"""Work out which source plasmid contributes which part of the target.

Given one or more source plasmids and a desired target, decide how to carve the target into
PCR fragments. The result is a list of `JunctionSpec`s, which is all the primer designer
needs to know.

The backbone is whichever source explains the most of the target. Comparing the target to
it gives a small number of differences; each difference is then looked up in the *other*
sources. A difference whose new sequence is a contiguous block of another plasmid becomes
its own amplified fragment (so a 600 bp cassette is PCR'd, not written into primer tails),
while anything not found elsewhere stays as non-templated sequence carried on the primers.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

from .design import DesignError, DesignParams, JunctionSpec
from .plasmid_diff import EditSite, PlasmidComparison, compare_plasmids
from .sequences import SeqRecord, find_all, revcomp

MIN_AMPLIFIED_INSERT = 60
"""Below this, new sequence is cheaper to put on primer tails than to amplify separately."""


@dataclass
class Segment:
    """A contiguous stretch of the target contributed by one source plasmid."""

    index: int
    source_name: str
    source_seq: str
    tgt_start: int           # target coordinate, may exceed len(target) when wrapping
    tgt_end: int
    source_offset: int       # source_pos = (target_pos + source_offset) % len(source_seq)
    reverse_complemented: bool = False

    @property
    def length(self) -> int:
        return self.tgt_end - self.tgt_start


@dataclass
class AssemblyPlan:
    """How the target is built: segments amplified from sources, joined at junctions."""

    target_name: str
    target: str
    backbone_name: str
    comparison: PlasmidComparison
    segments: List[Segment] = field(default_factory=list)
    specs: List[JunctionSpec] = field(default_factory=list)
    sources_used: List[str] = field(default_factory=list)
    notes: List[str] = field(default_factory=list)

    @property
    def fragment_count(self) -> int:
        return len(self.segments)


def choose_backbone(
    sources: Sequence[SeqRecord],
    target: SeqRecord,
) -> Tuple[SeqRecord, Dict[str, PlasmidComparison]]:
    """The source that explains most of the target becomes the backbone."""
    comparisons: Dict[str, PlasmidComparison] = {}
    best: Optional[SeqRecord] = None
    best_identity = -1.0
    for source in sources:
        comparison = compare_plasmids(source.name, source.seq, target.name, target.seq)
        comparisons[source.name] = comparison
        if comparison.identity > best_identity:
            best_identity, best = comparison.identity, source
    if best is None:
        raise DesignError("No source plasmids were supplied")
    return best, comparisons


def plan_assembly(
    sources: Sequence[SeqRecord],
    target: SeqRecord,
    params: DesignParams = DesignParams(),
    min_amplified_insert: int = MIN_AMPLIFIED_INSERT,
) -> AssemblyPlan:
    """Decide the fragment layout for building `target` from `sources`."""
    backbone, comparisons = choose_backbone(sources, target)
    comparison = comparisons[backbone.name]
    plan = AssemblyPlan(
        target_name=target.name,
        target=comparison.target,
        backbone_name=backbone.name,
        comparison=comparison,
        sources_used=[backbone.name],
    )
    if len(sources) > 1:
        plan.notes.append(
            f"{backbone.name} explains {100 * comparison.identity:.1f}% of {target.name} and is "
            "used as the backbone."
        )
    if comparison.identical:
        return plan

    others = [s for s in sources if s.name != backbone.name]
    resolved = []
    for edit in comparison.edits:
        resolved.extend(_resolve_edit(edit, others, comparison, min_amplified_insert))
    for res in resolved:
        if res.source is not None and res.source.name not in plan.sources_used:
            plan.sources_used.append(res.source.name)
        plan.notes.extend(res.notes)

    plan.segments = _build_segments(comparison, backbone, resolved)
    plan.specs = _build_specs(comparison, plan.segments, resolved, params)
    return plan


@dataclass
class _ResolvedEdit:
    """An edit, plus where its new sequence comes from (another plasmid, or nowhere)."""

    edit: EditSite
    source: Optional[SeqRecord] = None
    source_start: int = 0            # 0-based position in that source
    reverse_complemented: bool = False
    notes: List[str] = field(default_factory=list)

    @property
    def amplified(self) -> bool:
        return self.source is not None


def _resolve_edit(
    edit: EditSite,
    others: Sequence[SeqRecord],
    comparison: PlasmidComparison,
    min_amplified_insert: int,
) -> List[_ResolvedEdit]:
    """Work out where this edit's new sequence comes from.

    Usually the whole block sits in one other source. When it does not, the block may still
    be a join of two or more donor stretches (a cassette assembled from several plasmids),
    so fall back to tiling it; whatever no donor covers stays on the primer tails.
    """
    new = edit.tgt_seq
    if not new or not others:
        return [_ResolvedEdit(edit=edit)]

    for source in others:
        hit = _locate_block(source.seq, new)
        if hit is None:
            continue
        start, flipped = hit
        if len(new) < min_amplified_insert:
            return [_ResolvedEdit(
                edit=edit,
                notes=[
                    f"The {len(new)} nt of new sequence is present in {source.name} but is short "
                    "enough to carry on the primer tails, so it is not amplified separately."
                ],
            )]
        strand = " (reverse complement)" if flipped else ""
        return [_ResolvedEdit(
            edit=edit, source=source, source_start=start, reverse_complemented=flipped,
            notes=[
                f"The {len(new)} nt insert matches {source.name} {start + 1}.."
                f"{start + len(new)}{strand} exactly, so it is amplified from {source.name} as "
                "its own fragment."
            ],
        )]

    if len(new) >= min_amplified_insert:
        tiled = _tile_edit(edit, others, min_amplified_insert)
        if tiled is not None:
            return tiled
        return [_ResolvedEdit(
            edit=edit,
            notes=[
                f"The {len(new)} nt of new sequence was not found as a contiguous block in "
                + ", ".join(s.name for s in others)
                + ", and could not be tiled across them either; it will be carried on the "
                "primers instead."
            ],
        )]
    return [_ResolvedEdit(edit=edit)]


def _longest_donor_block(
    motif: str,
    others: Sequence[SeqRecord],
) -> Optional[Tuple[int, SeqRecord, int, bool]]:
    """Longest prefix of `motif` that sits uniquely in one donor. (length, donor, start, rc)"""
    best: Optional[Tuple[int, SeqRecord, int, bool]] = None
    for source in others:
        # Presence shrinks with length and uniqueness grows with it, so the longest prefix
        # that occurs at all is also the one most likely to occur exactly once.
        lo, hi = 0, min(len(motif), len(source.seq))
        while lo < hi:
            mid = (lo + hi + 1) // 2
            if _locate_block(source.seq, motif[:mid]) is not None:
                lo = mid
            else:
                hi = mid - 1
        if lo == 0:
            continue
        hit = _locate_block(source.seq, motif[:lo])
        if hit is None:
            continue
        if best is None or lo > best[0]:
            best = (lo, source, hit[0], hit[1])
    return best


def _tile_edit(
    edit: EditSite,
    others: Sequence[SeqRecord],
    min_amplified_insert: int,
) -> Optional[List[_ResolvedEdit]]:
    """Split one insertion across several donors, longest donor stretch first.

    Returns None when tiling buys nothing -- fewer than two donor stretches, which the
    single-block search above has already ruled out.
    """
    new = edit.tgt_seq
    pieces: List[Tuple[int, int, Optional[SeqRecord], int, bool]] = []  # start,end,src,pos,rc
    pos = 0
    while pos < len(new):
        found = _longest_donor_block(new[pos:], others)
        if found is None:
            # No donor covers the next base; walk forward to where one picks up again.
            nxt = pos + 1
            while nxt < len(new) and _longest_donor_block(new[nxt:], others) is None:
                nxt += 1
            pieces.append((pos, nxt, None, 0, False))
            pos = nxt
            continue
        length, source, start, flipped = found
        pieces.append((pos, pos + length, source, start, flipped))
        pos += length

    # A donor stretch too short to be worth its own PCR goes onto the tails instead.
    pieces = [
        p if p[2] is None or (p[1] - p[0]) >= min_amplified_insert
        else (p[0], p[1], None, 0, False)
        for p in pieces
    ]
    merged: List[Tuple[int, int, Optional[SeqRecord], int, bool]] = []
    for piece in pieces:
        if merged and piece[2] is None and merged[-1][2] is None:
            prev = merged[-1]
            merged[-1] = (prev[0], piece[1], None, 0, False)
        else:
            merged.append(piece)

    if sum(1 for p in merged if p[2] is not None) < 2:
        return None

    resolved: List[_ResolvedEdit] = []
    for start, end, source, src_pos, flipped in merged:
        sub = EditSite(
            # An insertion does not consume template, so every piece shares the one point.
            tpl_start=edit.tpl_start,
            tpl_end=edit.tpl_end,
            tgt_start=edit.tgt_start + start,
            tgt_end=edit.tgt_start + end,
            tpl_seq="",
            tgt_seq=new[start:end],
        )
        if source is None:
            resolved.append(_ResolvedEdit(edit=sub, notes=[
                f"{end - start} nt of the insert ({start + 1}..{end} of {len(new)}) is in no "
                "donor and is carried on the primer tails."
            ]))
        else:
            strand = " (reverse complement)" if flipped else ""
            resolved.append(_ResolvedEdit(
                edit=sub, source=source, source_start=src_pos, reverse_complemented=flipped,
                notes=[
                    f"Insert {start + 1}..{end} of {len(new)} nt matches {source.name} "
                    f"{src_pos + 1}..{src_pos + (end - start)}{strand}, so it is amplified "
                    f"from {source.name} as its own fragment."
                ],
            ))
    names = ", ".join(p[2].name for p in merged if p[2] is not None)
    resolved[0].notes.insert(
        0,
        f"The {len(new)} nt insert is not a single block of any one donor; it is tiled "
        f"across {names}.",
    )
    return resolved


def _locate_block(source: str, motif: str) -> Optional[Tuple[int, bool]]:
    """Find `motif` exactly once in a circular source, on either strand."""
    if len(motif) > len(source):
        return None
    extended = source + source[:len(motif) - 1]
    forward = sorted({h % len(source) for h in find_all(extended, motif)})
    reverse = sorted({h % len(source) for h in find_all(extended, revcomp(motif))})
    if len(forward) == 1 and not reverse:
        return forward[0], False
    if len(reverse) == 1 and not forward:
        return reverse[0], True
    return None


def _build_segments(
    comparison: PlasmidComparison,
    backbone: SeqRecord,
    resolved: Sequence[_ResolvedEdit],
) -> List[Segment]:
    """Carve the target into alternating backbone and amplified-insert segments."""
    target = comparison.target
    n = len(target)
    amplified = [r for r in resolved if r.amplified]
    if not amplified:
        # Single fragment: the whole target comes off the backbone, primers carry the edits.
        first = comparison.edits[0]
        start = first.tgt_end
        return [
            Segment(
                index=0, source_name=backbone.name, source_seq=comparison.template,
                tgt_start=start, tgt_end=start + n,
                source_offset=(first.tpl_end - first.tgt_end) % len(comparison.template),
            )
        ]

    segments: List[Segment] = []
    # Walk the circle: a backbone stretch, then an insert, then backbone, and so on.
    order = sorted(amplified, key=lambda r: r.edit.tgt_start)
    # New sequence no donor carries rides on the tails, so the backbone picks up only after
    # it. Measure how much of that trails each amplified piece.
    by_start = sorted(resolved, key=lambda r: r.edit.tgt_start)
    tail_after: Dict[int, int] = {}
    head_before: Dict[int, int] = {}
    for i, res in enumerate(by_start):
        if not res.amplified:
            continue
        tail = 0
        j = i + 1
        while j < len(by_start) and not by_start[j].amplified \
                and by_start[j].edit.tgt_start == res.edit.tgt_end + tail:
            tail += by_start[j].edit.new_length
            j += 1
        tail_after[res.edit.tgt_start] = tail
        head = 0
        j = i - 1
        while j >= 0 and not by_start[j].amplified \
                and by_start[j].edit.tgt_end == res.edit.tgt_start - head:
            head += by_start[j].edit.new_length
            j -= 1
        head_before[res.edit.tgt_start] = head
    for i, res in enumerate(order):
        edit = res.edit
        nxt = order[(i + 1) % len(order)]
        insert_seq = res.source.seq if res.source else ""
        if res.reverse_complemented:
            # Work on the strand that reads the same way as the target.
            insert_seq = revcomp(insert_seq)
            start = len(res.source.seq) - res.source_start - len(edit.tgt_seq)
        else:
            start = res.source_start
        segments.append(
            Segment(
                index=len(segments),
                source_name=res.source.name if res.source else "",
                source_seq=insert_seq,
                tgt_start=edit.tgt_start,
                tgt_end=edit.tgt_end,
                source_offset=(start - edit.tgt_start) % len(insert_seq),
                reverse_complemented=res.reverse_complemented,
            )
        )
        # Backbone stretch from the end of this insert to the start of the next one,
        # picking up after any non-templated sequence that trails this insert.
        back_start = edit.tgt_end + tail_after.get(edit.tgt_start, 0)
        # ...and hands over before any non-templated sequence leading the next insert.
        back_end = nxt.edit.tgt_start - head_before.get(nxt.edit.tgt_start, 0)
        if back_end == back_start:
            # Two donor pieces butt together: no backbone between them.
            continue
        if back_end < back_start:
            back_end += n
        segments.append(
            Segment(
                index=len(segments),
                source_name=backbone.name,
                source_seq=comparison.template,
                tgt_start=back_start,
                tgt_end=back_end,
                # The template resumes at tpl_end, which now sits at back_start.
                source_offset=(edit.tpl_end - back_start) % len(comparison.template),
            )
        )
    return segments


def _build_specs(
    comparison: PlasmidComparison,
    segments: Sequence[Segment],
    resolved: Sequence[_ResolvedEdit],
    params: DesignParams,
) -> List[JunctionSpec]:
    """One junction per segment boundary, carrying any primer-borne sequence with it."""
    target = comparison.target
    primer_borne = {r.edit.tgt_start: r.edit for r in resolved if not r.amplified}
    # A primer-borne piece that trails an amplified one ends where the next segment starts,
    # so it is reached by its end coordinate rather than its start.
    primer_borne_by_end = {r.edit.tgt_end: r.edit for r in resolved if not r.amplified}

    if len(segments) == 1:
        # Whole-plasmid amplification; each edit becomes its own junction on the backbone.
        specs: List[JunctionSpec] = []
        edits = comparison.edits
        for i, edit in enumerate(edits):
            prev_edit = edits[(i - 1) % len(edits)]
            next_edit = edits[(i + 1) % len(edits)]
            specs.append(
                JunctionSpec(
                    index=i,
                    target=target,
                    boundary=edit.tgt_start,
                    new_length=edit.new_length,
                    up_source_name=comparison.template_name,
                    up_source=comparison.template,
                    up_offset=(edit.tpl_start - edit.tgt_start) % len(comparison.template),
                    up_available=(edit.tgt_start - prev_edit.tgt_end) % len(target) or len(target),
                    down_source_name=comparison.template_name,
                    down_source=comparison.template,
                    down_offset=(edit.tpl_end - edit.tgt_end) % len(comparison.template),
                    down_available=(next_edit.tgt_start - edit.tgt_end) % len(target) or len(target),
                    label=_edit_label(edit),
                    edit=edit,
                )
            )
        return specs

    specs = []
    count = len(segments)
    for i, segment in enumerate(segments):
        upstream = segments[(i - 1) % count]
        boundary = segment.tgt_start
        edit = primer_borne_by_end.get(boundary % len(target)) \
            or primer_borne.get(boundary % len(target))
        new_length = edit.new_length if edit is not None else 0
        # The upstream segment's own frame ends at `upstream.tgt_end`, which is the same
        # target position as `boundary` but possibly a whole turn of the circle away.
        up_shift = upstream.tgt_end - (boundary - new_length)
        specs.append(
            JunctionSpec(
                index=i,
                target=target,
                boundary=boundary - new_length,
                new_length=new_length,
                up_source_name=upstream.source_name,
                up_source=upstream.source_seq,
                up_offset=upstream.source_offset,
                up_shift=up_shift,
                up_available=upstream.length,
                down_source_name=segment.source_name,
                down_source=segment.source_seq,
                down_offset=segment.source_offset,
                down_available=segment.length,
                label=f"{upstream.source_name} → {segment.source_name}",
                edit=edit,
            )
        )
    return specs


def _edit_label(edit: EditSite) -> str:
    return f"{edit.kind} of {max(edit.removed_length, edit.new_length)} nt"


def verify_plan(plan: AssemblyPlan) -> None:
    """Check the segment layout really does spell out the target, before designing primers."""
    if not plan.segments or len(plan.segments) == 1:
        return
    target = plan.target
    rebuilt = []
    for segment in sorted(plan.segments, key=lambda s: s.tgt_start):
        expected = "".join(
            target[(segment.tgt_start + k) % len(target)] for k in range(segment.length)
        )
        actual = "".join(
            segment.source_seq[(segment.tgt_start + k + segment.source_offset)
                               % len(segment.source_seq)]
            for k in range(segment.length)
        )
        if expected != actual:
            raise DesignError(
                f"Segment {segment.index + 1} from {segment.source_name} does not match the "
                f"target over {segment.length} nt starting at target "
                f"{segment.tgt_start % len(target) + 1}"
            )
        rebuilt.append(expected)
    # Sequence no source carries is not in any segment; the primer tails spell it out.
    on_tails = sum(spec.new_length for spec in plan.specs)
    total = sum(s.length for s in plan.segments)
    if total + on_tails != len(target):
        raise DesignError(
            f"Segments cover {total} nt and primer tails add {on_tails} nt, but the target "
            f"is {len(target)} nt"
        )
