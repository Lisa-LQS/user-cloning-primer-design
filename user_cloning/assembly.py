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
    resolved = [
        _resolve_edit(edit, others, comparison, min_amplified_insert)
        for edit in comparison.edits
    ]
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
) -> _ResolvedEdit:
    """Look for this edit's new sequence as a contiguous block in one of the other sources."""
    new = edit.tgt_seq
    if not new or not others:
        return _ResolvedEdit(edit=edit)

    for source in others:
        hit = _locate_block(source.seq, new)
        if hit is None:
            continue
        start, flipped = hit
        if len(new) < min_amplified_insert:
            return _ResolvedEdit(
                edit=edit,
                notes=[
                    f"The {len(new)} nt of new sequence is present in {source.name} but is short "
                    "enough to carry on the primer tails, so it is not amplified separately."
                ],
            )
        strand = " (reverse complement)" if flipped else ""
        return _ResolvedEdit(
            edit=edit, source=source, source_start=start, reverse_complemented=flipped,
            notes=[
                f"The {len(new)} nt insert matches {source.name} {start + 1}.."
                f"{start + len(new)}{strand} exactly, so it is amplified from {source.name} as "
                "its own fragment."
            ],
        )

    if len(new) >= min_amplified_insert and others:
        return _ResolvedEdit(
            edit=edit,
            notes=[
                f"The {len(new)} nt of new sequence was not found as a contiguous block in "
                + ", ".join(s.name for s in others)
                + "; it will be carried on the primers instead."
            ],
        )
    return _ResolvedEdit(edit=edit)


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
        # Backbone stretch from the end of this insert to the start of the next one.
        back_start = edit.tgt_end
        back_end = nxt.edit.tgt_start
        if back_end <= back_start:
            back_end += n
        segments.append(
            Segment(
                index=len(segments),
                source_name=backbone.name,
                source_seq=comparison.template,
                tgt_start=back_start,
                tgt_end=back_end,
                source_offset=(edit.tpl_end - edit.tgt_end) % len(comparison.template),
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
        edit = primer_borne.get(boundary % len(target))
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
    total = sum(s.length for s in plan.segments)
    if total != len(target):
        raise DesignError(
            f"Segments cover {total} nt but the target is {len(target)} nt"
        )
