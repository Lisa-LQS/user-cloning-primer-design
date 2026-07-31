"""Compare a template plasmid to a desired target plasmid, treating both as circles.

The output is a small number of `EditSite` objects: contiguous stretches where the
target differs from the template. Everything outside those sites is identical and can
therefore be amplified straight off the template, which is what makes USER cloning
possible in the first place.
"""

from __future__ import annotations

import difflib
from dataclasses import dataclass, field
from typing import List, Optional, Tuple

from .sequences import kmer_index, revcomp, rotate


@dataclass
class EditSite:
    """One contiguous difference: template[tpl_start:tpl_end] becomes target[tgt_start:tgt_end].

    Coordinates are on the *rotated* (common-frame) linear sequences produced by
    `compare_plasmids`, not on the original spreadsheet sequences.
    """

    tpl_start: int
    tpl_end: int
    tgt_start: int
    tgt_end: int
    tpl_seq: str
    tgt_seq: str

    @property
    def kind(self) -> str:
        if not self.tpl_seq:
            return "insertion"
        if not self.tgt_seq:
            return "deletion"
        if len(self.tpl_seq) == len(self.tgt_seq):
            return "substitution" if len(self.tpl_seq) > 1 else "point mutation"
        return "replacement"

    @property
    def new_length(self) -> int:
        return len(self.tgt_seq)

    @property
    def removed_length(self) -> int:
        return len(self.tpl_seq)

    def describe(self) -> str:
        def show(s: str, limit: int = 30) -> str:
            if not s:
                return "-"
            return s if len(s) <= limit else f"{s[:12]}...{s[-12:]} ({len(s)} nt)"

        return (
            f"{self.kind}: template {self.tpl_start + 1}..{self.tpl_end} "
            f"[{show(self.tpl_seq)}] -> target {self.tgt_start + 1}..{self.tgt_end} [{show(self.tgt_seq)}]"
        )


@dataclass
class PlasmidComparison:
    """Result of aligning target to template in a shared circular frame."""

    template_name: str
    target_name: str
    template: str          # rotated template, common frame
    target: str            # rotated target, common frame
    edits: List[EditSite] = field(default_factory=list)
    template_rotation: int = 0   # left-rotation applied to the original template
    target_rotation: int = 0     # left-rotation applied to the original target
    target_reverse_complemented: bool = False
    identity: float = 0.0
    notes: List[str] = field(default_factory=list)

    @property
    def identical(self) -> bool:
        return not self.edits

    def to_template_input(self, frame_pos: int) -> int:
        """Convert an aligned-frame template coordinate back to the input file's numbering."""
        return (frame_pos + self.template_rotation) % len(self.template)

    def to_target_input(self, frame_pos: int) -> int:
        """Convert an aligned-frame target coordinate back to the input file's numbering.

        If the target was reverse-complemented to align it, the returned coordinate refers
        to the reverse complement of the input sequence, not the input sequence itself.
        """
        return (frame_pos + self.target_rotation) % len(self.target)

    def conserved_blocks(self) -> List[Tuple[int, int, int]]:
        """Maximal (tpl_start, tgt_start, size) blocks that are identical in both."""
        blocks: List[Tuple[int, int, int]] = []
        tpl_pos = tgt_pos = 0
        for edit in sorted(self.edits, key=lambda e: e.tgt_start):
            size = edit.tgt_start - tgt_pos
            if size > 0:
                blocks.append((tpl_pos, tgt_pos, size))
            tpl_pos, tgt_pos = edit.tpl_end, edit.tgt_end
        tail = len(self.target) - tgt_pos
        if tail > 0:
            blocks.append((tpl_pos, tgt_pos, tail))
        return blocks


def compare_plasmids(
    template_name: str,
    template: str,
    target_name: str,
    target: str,
    kmer: int = 25,
    merge_gap: int = 15,
    allow_reverse_complement: bool = True,
) -> PlasmidComparison:
    """Align two circular plasmids and return their differences.

    Steps: (1) pick the target orientation, (2) rotate the target onto the template,
    (3) move the linear origin into a long conserved block so no edit straddles it,
    (4) diff and merge nearby differences into single edit sites.
    """
    notes: List[str] = []
    orientation_flipped = False
    forward_offset, forward_votes = _best_offset(template, target, kmer)

    if allow_reverse_complement:
        rc = revcomp(target)
        rc_offset, rc_votes = _best_offset(template, rc, kmer)
        if rc_votes > forward_votes:
            notes.append(
                "Target matches the template better on the opposite strand; "
                "it was reverse-complemented before designing."
            )
            target, forward_offset, forward_votes = rc, rc_offset, rc_votes
            orientation_flipped = True

    if forward_votes == 0:
        notes.append(
            f"No shared {kmer}-mer between the two plasmids: they may be unrelated, in which "
            "case single-template USER cloning is not applicable."
        )
        rotation = 0
    else:
        rotation = (-forward_offset) % len(target)
    target_rot = rotate(target, rotation)

    # Put the linear origin inside a long conserved block so that no edit wraps around it.
    # The two shifts differ whenever the plasmids are of unequal length: each sequence is
    # rotated to the *same base of the shared anchor*, not by the same number of bases.
    shift_tpl, shift_tgt = _origin_shift(template, target_rot)
    template_frame = rotate(template, shift_tpl)
    target_frame = rotate(target_rot, shift_tgt)
    template_rotation = shift_tpl % len(template)
    target_rotation = (rotation + shift_tgt) % len(target)

    edits = _diff_edits(template_frame, target_frame, merge_gap)
    matched = len(target_frame) - sum(max(e.removed_length, e.new_length) for e in edits)
    identity = matched / max(len(template_frame), len(target_frame)) if target_frame else 0.0

    for edit in edits:
        if edit.tgt_start == 0 or edit.tgt_end == len(target_frame):
            notes.append(
                "An edit sits at the linear origin of the common frame; junction placement "
                "there is handled circularly but double-check the report diagram."
            )
            break

    return PlasmidComparison(
        template_name=template_name,
        target_name=target_name,
        template=template_frame,
        target=target_frame,
        edits=edits,
        template_rotation=template_rotation,
        target_rotation=target_rotation,
        target_reverse_complemented=orientation_flipped,
        identity=identity,
        notes=notes,
    )


def _best_offset(template: str, target: str, kmer: int) -> Tuple[int, int]:
    """Most common (template_pos - target_pos) offset over shared k-mers, and its vote count."""
    if len(template) < kmer or len(target) < kmer:
        return 0, 0
    index = kmer_index(template, kmer, circular=True)
    votes: dict = {}
    step = max(1, len(target) // 2000)  # sample; plasmids are near-identical so this is ample
    for i in range(0, len(target) - kmer + 1, step):
        positions = index.get(target[i:i + kmer])
        if not positions or len(positions) > 4:  # skip repeats
            continue
        for p in positions:
            offset = (p - i) % len(template)
            votes[offset] = votes.get(offset, 0) + 1
    if not votes:
        return 0, 0
    best = max(votes.items(), key=lambda kv: kv[1])
    return best[0], best[1]


def _origin_shift(template: str, target: str, min_anchor: int = 60) -> Tuple[int, int]:
    """Left-rotations that land both linear origins on the same base of a long shared match.

    Returns (template_shift, target_shift). The two differ by the net indel offset between
    the start of each sequence and the anchor, which is exactly what keeps the rotated
    sequences in register.
    """
    matcher = difflib.SequenceMatcher(None, template, target, autojunk=False)
    best = matcher.find_longest_match(0, len(template), 0, len(target))
    if best.size < min_anchor:
        return 0, 0
    midpoint = best.size // 2
    return best.a + midpoint, best.b + midpoint


def _diff_edits(template: str, target: str, merge_gap: int) -> List[EditSite]:
    """Diff two aligned linear sequences, merging differences separated by short matches."""
    opcodes = difflib.SequenceMatcher(None, template, target, autojunk=False).get_opcodes()
    raw = [(a1, a2, b1, b2) for tag, a1, a2, b1, b2 in opcodes if tag != "equal"]
    if not raw:
        return []

    merged: List[List[int]] = [list(raw[0])]
    for a1, a2, b1, b2 in raw[1:]:
        prev = merged[-1]
        if a1 - prev[1] <= merge_gap and b1 - prev[3] <= merge_gap:
            prev[1], prev[3] = a2, b2
        else:
            merged.append([a1, a2, b1, b2])

    return [
        EditSite(
            tpl_start=a1, tpl_end=a2, tgt_start=b1, tgt_end=b2,
            tpl_seq=template[a1:a2], tgt_seq=target[b1:b2],
        )
        for a1, a2, b1, b2 in merged
    ]
