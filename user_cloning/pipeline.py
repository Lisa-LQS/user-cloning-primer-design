"""Top-level orchestration: spreadsheet in, verified USER cloning design out."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

from .design import (
    DesignError,
    DesignParams,
    Fragment,
    Junction,
    Primer,
    build_fragments,
    design_junctions,
)
from .assembly import plan_assembly, verify_plan
from .plasmid_diff import PlasmidComparison
from .sequences import RepeatBlock, SeqRecord, repeat_blocks, revcomp, rotate
from .simulate import Verification, verify_design

STATUS_OK = "ok"
STATUS_NO_DIFFERENCE = "no_difference"
STATUS_FAILED = "failed"



@dataclass
class PcrProtocol:
    """Suggested cycling conditions for one fragment."""

    fragment: str
    annealing_temp_C: float
    extension_seconds: int
    product_length: int


@dataclass
class DesignResult:
    """Everything produced for one template -> target request."""

    design_id: str
    design_date: str
    template_name: str
    target_name: str
    template_length: int
    target_length: int
    status: str
    batch: str = ""
    comparison: Optional[PlasmidComparison] = None
    junctions: List[Junction] = field(default_factory=list)
    fragments: List[Fragment] = field(default_factory=list)
    verification: Optional[Verification] = None
    protocols: List[PcrProtocol] = field(default_factory=list)
    params: DesignParams = DesignParams()
    template_repeats: List[RepeatBlock] = field(default_factory=list)
    source_names: List[str] = field(default_factory=list)
    backbone_name: str = ""
    warnings: List[str] = field(default_factory=list)
    messages: List[str] = field(default_factory=list)

    @property
    def primers(self) -> List[Primer]:
        out: List[Primer] = []
        for junction in self.junctions:
            out.extend([junction.forward, junction.reverse])
        return out

    @property
    def verified(self) -> bool:
        return self.verification is not None and self.verification.ok

    @property
    def stamp(self) -> str:
        """Date plus batch label, as it appears in every file and oligo name."""
        base = self.design_date.replace("-", "")
        return f"{base}_{self.batch}" if self.batch else base

    def summary_row(self) -> Dict[str, object]:
        return {
            "design_id": self.design_id,
            "design_date": self.design_date,
            "batch": self.batch,
            "target": self.target_name,
            "template": self.template_name,
            "status": self.status,
            "verified": "yes" if self.verified else "no",
            "edits": len(self.comparison.edits) if self.comparison else 0,
            "junctions": len(self.junctions),
            "fragments": len(self.fragments),
            "primers": len(self.primers),
            "warnings": len(self.warnings),
            "notes": "; ".join(self.messages) or "-",
        }


def resolve_records(
    records: Sequence[SeqRecord],
    template_name: str,
    target_name: str,
) -> Tuple[SeqRecord, SeqRecord]:
    by_name = {rec.name: rec for rec in records}
    lowered = {rec.name.lower(): rec for rec in records}
    missing = [n for n in (template_name, target_name) if n not in by_name and n.lower() not in lowered]
    if missing:
        raise DesignError(
            f"Sequence(s) {', '.join(missing)} not found in the input file. "
            f"Available: {', '.join(rec.name for rec in records)}"
        )
    pick = lambda n: by_name.get(n) or lowered[n.lower()]  # noqa: E731
    return pick(template_name), pick(target_name)


def design_one(
    template: SeqRecord,
    target: SeqRecord,
    design_date: str,
    params: DesignParams = DesignParams(),
    design_index: int = 1,
    batch: str = "",
) -> DesignResult:
    """Design (and verify) the USER cloning primers for one template -> target pair."""
    return design_assembly([template], target, design_date, params, design_index, batch)


def design_assembly(
    sources: Sequence[SeqRecord],
    target: SeqRecord,
    design_date: str,
    params: DesignParams = DesignParams(),
    design_index: int = 1,
    batch: str = "",
) -> DesignResult:
    """Design (and verify) a USER assembly of `target` from one or more source plasmids.

    `batch` labels this run so that several batches requested on the same day stay separate;
    it appears in the design ID, the output paths, the file names and the oligo names.
    """
    if not sources:
        raise DesignError("At least one source plasmid is required")
    stamp = design_date.replace("-", "")
    batch_part = f"-{batch}" if batch else ""
    design_id = f"USER-{stamp}{batch_part}-{target.name}-{design_index:02d}"
    result = DesignResult(
        design_id=design_id,
        design_date=design_date,
        batch=batch,
        template_name=", ".join(s.name for s in sources),
        target_name=target.name,
        template_length=sources[0].length,
        target_length=target.length,
        status=STATUS_FAILED,
        params=params,
        source_names=[s.name for s in sources],
    )

    try:
        plan = plan_assembly(sources, target, params)
        verify_plan(plan)
    except DesignError as exc:
        result.messages.append(str(exc))
        return result

    comparison = plan.comparison
    backbone = next(s for s in sources if s.name == plan.backbone_name)
    result.comparison = comparison
    result.backbone_name = plan.backbone_name
    result.template_length = backbone.length
    result.warnings.extend(comparison.notes)
    result.messages.extend(plan.notes)

    if comparison.identical:
        result.status = STATUS_NO_DIFFERENCE
        result.messages.append(
            f"{target.name} and {plan.backbone_name} are identical over all "
            f"{backbone.length} bp (as circular sequences), so there is nothing to clone. "
            "Check that the input file holds the intended target sequence."
        )
        return result

    # Repeats are found on the backbone as supplied, so the coordinates reported to the user
    # match their own file rather than the internal aligned frame.
    result.template_repeats = repeat_blocks(backbone.seq)
    for i, edit in enumerate(comparison.edits, start=1):
        start = comparison.to_template_input(edit.tpl_start)
        end = start + max(edit.removed_length, 1)
        for block in result.template_repeats:
            wrapped_end = min(end, backbone.length)
            overlaps = block.overlaps(start, wrapped_end) or (
                end > backbone.length and block.overlaps(0, end - backbone.length)
            )
            if overlaps:
                result.warnings.append(
                    f"Edit site {i} lies in a {block.length} nt sequence that occurs "
                    f"{block.copies}x in {plan.backbone_name} "
                    f"({block.start + 1}..{block.end + 1}); primers there have to be long to "
                    "become specific, and may not be possible at all."
                )

    if comparison.identity < 0.5 and len(plan.segments) <= 1:
        result.messages.append(
            f"Only {100 * comparison.identity:.1f}% of the target is present in "
            f"{plan.backbone_name}; a single-template design will need very long primers. "
            "Supplying the plasmid that carries the rest would let it be amplified as a "
            "second fragment."
        )

    try:
        junctions, junction_warnings = design_junctions(plan.specs, params)
        result.junctions = junctions
        result.warnings.extend(junction_warnings)
        result.fragments = build_fragments(plan.specs, junctions)
    except DesignError as exc:
        result.messages.append(str(exc))
        return result

    for junction in junctions:
        result.warnings.extend(f"J{junction.index + 1}: {w}" for w in junction.warnings)
        carried = junction.block_len + len(junction.forward.extra) + len(junction.reverse.extra)
        if carried > 80:
            result.messages.append(
                f"Junction {junction.index + 1} carries {carried} nt of new sequence on its "
                "primer pair. That works, but ordering the new region as a synthesised dsDNA "
                "fragment (with matching USER overhangs) and assembling it as a second "
                "fragment is usually cheaper and more accurate than two long Ultramers."
            )

    result.verification = verify_design(comparison.target, junctions, result.fragments)
    result.protocols = [_protocol(frag) for frag in result.fragments]

    if result.verification.ok:
        result.status = STATUS_OK
        result.messages.append(
            "In-silico PCR + USER digestion + annealing reproduces the requested target "
            "plasmid exactly."
        )
    else:
        result.messages.append(
            "Verification failed: " + "; ".join(result.verification.errors)
        )
    return result


def _protocol(fragment: Fragment, extension_s_per_kb: int = 30) -> PcrProtocol:
    """Phusion-style cycling suggestion. Ta uses the annealing regions only, because the
    non-templated 5' tails contribute nothing in the first cycles."""
    limiting_tm = min(fragment.forward.anneal_tm, fragment.reverse.anneal_tm)
    ta = min(72.0, max(55.0, limiting_tm + 3.0))
    extension = max(30, int(round(extension_s_per_kb * fragment.expected_length / 1000.0)))
    return PcrProtocol(
        fragment=fragment.name,
        annealing_temp_C=round(ta, 1),
        extension_seconds=extension,
        product_length=fragment.expected_length,
    )


def predicted_plasmid(result: DesignResult) -> str:
    """The assembled plasmid, rotated back to the frame of the input target sequence."""
    if result.verification is None or not result.verification.assembled:
        return ""
    assembled = result.verification.assembled
    rotation = result.verification.rotation_vs_target
    if rotation is None:
        return assembled
    # `assembled` equals the common-frame target rotated left by `rotation`; undo that,
    # then undo the rotation (and any flip) applied when the target entered that frame,
    # so the output is directly comparable to the sequence in the input file.
    in_frame = rotate(assembled, -rotation)
    if result.comparison is None:
        return in_frame
    oriented = rotate(in_frame, -result.comparison.target_rotation)
    if result.comparison.target_reverse_complemented:
        return revcomp(oriented)
    return oriented
