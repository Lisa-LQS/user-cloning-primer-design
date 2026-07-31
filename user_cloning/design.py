"""USER cloning primer design.

Design model
------------
USER (Uracil-Specific Excision Reagent) cloning: PCR primers carry a single
deoxyuridine a few bases in from the 5' end and are amplified with a uracil-tolerant
proofreading polymerase (Phusion U, Q5U, PfuTurbo Cx). USER enzyme (UDG + Endo VIII)
excises the uracil and nicks the backbone, releasing the short 5' flap of each strand
and leaving a 3' single-stranded overhang on the opposite strand. Fragments whose
overhangs are complementary anneal and are transformed directly; E. coli seals the nicks.

For a fragment amplified with

    forward primer  =  F + (extra) + anneal_f          (U replaces the last base of F)
    reverse primer  =  revcomp(F') + (extra) + anneal_r

the post-USER 3' overhangs are `revcomp(tail_R)` on the top strand and `revcomp(tail_F)`
on the bottom strand. They can only anneal if

    tail_R == revcomp(tail_F)

so each junction is defined by a single overhang block F, and because both tails must
end in the U-substituted T, **F must begin with A and end with T**. That single
constraint drives the whole junction search below.

A pleasant consequence: every primer is just a contiguous (circular) slice of the
*desired target* sequence, with one T swapped for U --

    forward primer  =  target[a : a + L + extra_f + anneal_f]
    reverse primer  =  revcomp(target[a - extra_r - anneal_r : a + L])

where the overhang block is target[a : a+L].
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

from . import thermo
from .plasmid_diff import EditSite, PlasmidComparison
from .sequences import (
    circular_slice,
    count_circular_occurrences,
    gc_fraction,
    revcomp,
)
from .thermo import DEFAULT_CONDITIONS, PcrConditions

U_MARK = "/ideoxyU/"  # IDT ordering notation for deoxyuridine


@dataclass(frozen=True)
class DesignParams:
    """Tunable design constraints. Defaults suit Phusion U on an 8 kb plasmid."""

    overhang_min: int = 8
    overhang_max: int = 12
    overhang_preferred: int = 9
    overhang_gc_min: float = 20.0
    overhang_gc_max: float = 80.0

    anneal_min: int = 18
    anneal_max: int = 32
    anneal_hard_max: int = 60
    """Only used as a fallback: plasmids often carry repeated elements (duplicated
    promoters, LTRs), and inside those a 32 nt primer can have several binding sites. The
    annealing region is then extended until it reaches sequence unique in the template."""
    tm_target: float = 62.0
    tm_min: float = 57.0
    tm_max: float = 69.0
    max_pair_tm_diff: float = 4.0

    gc_min: float = 30.0
    gc_max: float = 70.0
    max_homopolymer: int = 5

    soft_max_primer_len: int = 60   # standard desalted oligo limit at most vendors
    hard_max_primer_len: int = 120  # beyond this an Ultramer/gene fragment is saner

    junction_window_slack: int = 30  # how far from the edit the overhang block may sit
    candidates_per_junction: int = 40
    beam_width: int = 60
    conditions: PcrConditions = DEFAULT_CONDITIONS


@dataclass
class Primer:
    """One orderable oligo."""

    name: str
    junction: int
    direction: str            # "forward" | "reverse"
    sequence: str             # plain DNA, U written as T
    u_index: int              # 0-based position of the deoxyuridine
    tail: str                 # 5' flap removed by USER, including the U position
    extra: str                # non-templated sequence retained in the product
    anneal: str               # 3' region that base-pairs with the template
    anneal_tm: float
    full_tm: float
    gc_percent: float
    template_hits: int
    warnings: List[str] = field(default_factory=list)

    @property
    def length(self) -> int:
        return len(self.sequence)

    @property
    def order_sequence(self) -> str:
        """Sequence in IDT notation, e.g. ACGT/ideoxyU/GGCA."""
        return f"{self.sequence[:self.u_index]}{U_MARK}{self.sequence[self.u_index + 1:]}"

    @property
    def sequence_with_u(self) -> str:
        """Sequence with a plain U at the deoxyuridine, e.g. ACGTUGGCA.

        Easier to read than the vendor notation, and what most sequence editors expect.
        """
        return f"{self.sequence[:self.u_index]}U{self.sequence[self.u_index + 1:]}"

    @property
    def overhang(self) -> str:
        """The 3' single-stranded overhang exposed at this primer's end of the fragment.

        USER removes this primer's 5' flap, so the overhang sits on the opposite strand
        and reads as the reverse complement of the flap.
        """
        return revcomp(self.tail)

    def stats_row(self) -> Dict[str, object]:
        return {
            "primer_name": self.name,
            "junction": self.junction + 1,
            "direction": self.direction,
            "order_sequence": self.order_sequence,
            "sequence_with_U": self.sequence_with_u,
            "plain_sequence": self.sequence,
            "length_nt": self.length,
            "u_position_from_5prime": self.u_index + 1,
            "tail_removed_by_USER": self.tail,
            "non_templated_extra": self.extra or "-",
            "annealing_region": self.anneal,
            "annealing_length_nt": len(self.anneal),
            "annealing_tm_C": round(self.anneal_tm, 1),
            "full_length_tm_C": round(self.full_tm, 1),
            "gc_percent": round(self.gc_percent, 1),
            "template_binding_sites": self.template_hits,
            "warnings": "; ".join(self.warnings) or "-",
        }


@dataclass
class JunctionCandidate:
    """One possible overhang-block placement for a junction, with its two primers."""

    edit_index: int
    block_start: int          # target coordinate of the overhang block
    block_len: int
    overhang: str             # the block itself, i.e. the top-strand 3' overhang
    forward: Primer
    reverse: Primer
    overhang_tm: float
    penalty: float
    penalty_detail: Dict[str, float]
    warnings: List[str] = field(default_factory=list)


@dataclass
class Junction:
    """A chosen junction: where two fragment ends are joined by complementary overhangs."""

    index: int
    edit: Optional[EditSite]
    block_start: int
    block_len: int
    overhang: str
    forward: Primer
    reverse: Primer
    overhang_tm: float
    penalty: float
    penalty_detail: Dict[str, float]
    warnings: List[str] = field(default_factory=list)
    label: str = ""
    up_source_name: str = ""
    down_source_name: str = ""

    @property
    def crosses_sources(self) -> bool:
        return bool(self.up_source_name) and self.up_source_name != self.down_source_name

    @property
    def top_overhang(self) -> str:
        return self.overhang

    @property
    def bottom_overhang(self) -> str:
        return revcomp(self.overhang)


@dataclass
class Fragment:
    """One PCR product, amplified from one source plasmid and joined at two junctions."""

    index: int
    name: str
    source_name: str
    source_seq: str
    forward: Primer
    reverse: Primer
    left_junction: int
    right_junction: int
    template_start: int      # source coordinate where the annealing region starts
    template_end: int        # source coordinate just past the reverse annealing region
    expected_length: int     # length of the PCR product before USER treatment


@dataclass
class JunctionSpec:
    """Where and how one junction must be built, independent of how it arose.

    A junction is the point in the target where the sequence stops being contributed by
    one source plasmid and starts being contributed by the next. Between the two there may
    be a stretch that no source provides (a mutation or a small insertion); that stretch is
    `target[boundary : boundary + new_length]` and has to be carried on the primers.

    The upstream source supplies target sequence up to `boundary` and is amplified by the
    reverse primer; the downstream source supplies target sequence from
    `boundary + new_length` onwards and is amplified by the forward primer. Coordinates map
    as `source_pos = (target_pos + offset) % len(source)`.

    This one shape covers everything: a point mutation on a single template (both sources
    the same plasmid, `new_length` = 1), a deletion (`new_length` = 0, offsets differing by
    the deleted length), and a cassette swap between two plasmids (different sources).

    The two sides may need different representatives of the same target position: target and
    source have different lengths, so reducing a target coordinate modulo the target length
    shifts the mapped source position by a whole turn of the circle. `up_shift` (a multiple
    of the target length, usually 0) is added on the upstream side to undo that.
    """

    index: int
    target: str
    boundary: int
    new_length: int
    up_source_name: str
    up_source: str
    up_offset: int
    up_available: int        # how far back the reverse primer may reach, from `boundary`
    down_source_name: str
    down_source: str
    down_offset: int
    down_available: int      # how far forward the forward primer may reach
    label: str
    up_shift: int = 0
    edit: Optional[EditSite] = None

    @property
    def new_sequence(self) -> str:
        return circular_slice(self.target, self.boundary, self.boundary + self.new_length)

    @property
    def crosses_sources(self) -> bool:
        return self.up_source_name != self.down_source_name


class DesignError(Exception):
    """Raised when no valid USER design can be produced."""


def design_junctions(
    specs: Sequence[JunctionSpec],
    params: DesignParams = DesignParams(),
) -> Tuple[List[Junction], List[str]]:
    """Pick one overhang block (and its primer pair) per junction."""
    if not specs:
        return [], []

    warnings: List[str] = []
    candidate_sets: List[List[JunctionCandidate]] = []
    for spec in specs:
        candidates = _junction_candidates(spec, params)
        if not candidates:
            raise DesignError(
                f"No usable USER junction found at junction {spec.index + 1} ({spec.label}). "
                + _diagnose_no_junction(spec, params)
            )
        candidate_sets.append(candidates[: params.candidates_per_junction])

    chosen = _select_compatible(candidate_sets, params)
    junctions = [
        Junction(
            index=i,
            edit=specs[i].edit,
            block_start=c.block_start,
            block_len=c.block_len,
            overhang=c.overhang,
            forward=c.forward,
            reverse=c.reverse,
            overhang_tm=c.overhang_tm,
            penalty=c.penalty,
            penalty_detail=c.penalty_detail,
            warnings=list(c.warnings),
            label=specs[i].label,
            up_source_name=specs[i].up_source_name,
            down_source_name=specs[i].down_source_name,
        )
        for i, c in enumerate(chosen)
    ]
    if len(junctions) > 1:
        warnings.extend(_check_overhang_orthogonality(junctions, params))
    return junctions, warnings


def build_fragments(
    specs: Sequence[JunctionSpec],
    junctions: Sequence[Junction],
) -> List[Fragment]:
    """Pair up primers into PCR fragments.

    Fragment i runs from junction i to junction i+1, so it is amplified from the source
    that lies downstream of junction i -- which must be the same source that lies upstream
    of junction i+1.
    """
    fragments: List[Fragment] = []
    n = len(junctions)
    for i, left in enumerate(junctions):
        right = junctions[(i + 1) % n]
        left_spec, right_spec = specs[i], specs[(i + 1) % n]
        if left_spec.down_source_name != right_spec.up_source_name:
            raise DesignError(
                f"Fragment {i + 1} is inconsistent: junction {left.index + 1} continues into "
                f"{left_spec.down_source_name} but junction {right.index + 1} expects "
                f"{right_spec.up_source_name}"
            )
        source_name, source = left_spec.down_source_name, left_spec.down_source
        fwd, rev = left.forward, right.reverse
        t_start = _locate_unique(source, fwd.anneal)
        rev_site = _locate_unique(source, revcomp(rev.anneal))
        t_end = rev_site + len(rev.anneal)
        span = (t_end - t_start) % len(source)
        if span == 0:
            span = len(source)
        product_len = span + len(fwd.tail) + len(fwd.extra) + len(rev.tail) + len(rev.extra)
        fragments.append(
            Fragment(
                index=i,
                name=f"F{i + 1}",
                source_name=source_name,
                source_seq=source,
                forward=fwd,
                reverse=rev,
                left_junction=left.index,
                right_junction=right.index,
                template_start=t_start,
                template_end=t_end,
                expected_length=product_len,
            )
        )
    return fragments


# --------------------------------------------------------------------------------------
# Junction candidate enumeration
# --------------------------------------------------------------------------------------


def _junction_candidates(
    spec: JunctionSpec,
    params: DesignParams,
) -> List[JunctionCandidate]:
    target = spec.target
    gs, ge = spec.boundary, spec.boundary + spec.new_length

    candidates: List[JunctionCandidate] = []
    slack = params.junction_window_slack
    # `a` is deliberately left un-normalised (it may be negative or exceed the sequence
    # length); circular_slice does the wrapping, and plain min/max then works below.
    for block_len in range(params.overhang_min, params.overhang_max + 1):
        for a in range(gs - block_len - slack, ge + slack + 1):
            block = circular_slice(target, a, a + block_len)
            if block[0] != "A" or block[-1] != "T":
                continue  # both tails must end in a T that becomes the deoxyuridine
            gc = 100.0 * gc_fraction(block)
            if not params.overhang_gc_min <= gc <= params.overhang_gc_max:
                continue
            cand = _build_candidate(spec, a, block_len, block, params)
            if cand is not None:
                candidates.append(cand)

    candidates.sort(key=lambda c: c.penalty)
    return candidates


def _diagnose_no_junction(spec: JunctionSpec, params: DesignParams) -> str:
    """Explain *why* a junction could not be primed, so the message is actionable."""
    up_end = spec.boundary + spec.up_shift + spec.up_offset
    down_start = spec.boundary + spec.new_length + spec.down_offset
    reasons: List[str] = []

    for label, source_name, source, motif_of in (
        ("upstream", spec.up_source_name, spec.up_source,
         lambda n: circular_slice(spec.up_source, up_end - n, up_end)),
        ("downstream", spec.down_source_name, spec.down_source,
         lambda n: circular_slice(spec.down_source, down_start, down_start + n)),
    ):
        needed = None
        for n in range(params.anneal_min, params.anneal_hard_max + 1):
            if count_circular_occurrences(source, motif_of(n), both_strands=True) == 1:
                needed = n
                break
        if needed is None:
            copies = count_circular_occurrences(
                source, motif_of(params.anneal_hard_max), both_strands=True
            )
            reasons.append(
                f"the {label} flank is still repeated {copies}x in {source_name} at "
                f"{params.anneal_hard_max} nt, so no primer there can be specific"
            )
        elif needed > params.anneal_max:
            reasons.append(
                f"the {label} flank needs {needed} nt to become unique in {source_name}"
            )

    if min(spec.up_available, spec.down_available) < params.anneal_min:
        reasons.append(
            f"only {spec.up_available} nt upstream and {spec.down_available} nt downstream are "
            f"available before the neighbouring junction, against a {params.anneal_min} nt "
            "minimum annealing length"
        )

    if reasons:
        return (
            "Cause: " + "; ".join(reasons) + ". Move the junction, raise --anneal-max, or "
            "build this region as a separate synthesised fragment."
        )
    return (
        "No overhang block of the form A...T with acceptable GC could be placed near this "
        "site. Try widening --overhang-min/--overhang-max or --junction-slack."
    )


def _build_candidate(
    spec: JunctionSpec,
    a: int,
    block_len: int,
    block: str,
    params: DesignParams,
) -> Optional[JunctionCandidate]:
    target = spec.target
    gs, ge = spec.boundary, spec.boundary + spec.new_length

    # Where the templated (annealing) regions must start / end, in target coordinates.
    # Anything between those bounds and the overhang block is sequence no source provides,
    # so it has to be carried on the primer as non-templated "extra".
    ar_end = min(a, gs)
    af_start = max(a + block_len, ge)
    extra_r = circular_slice(target, ar_end, a)
    extra_f = circular_slice(target, a + block_len, af_start)

    # Moving the block away from the boundary eats into the reach of the annealing region,
    # because each primer may only anneal within its own source's contribution.
    avail_left = spec.up_available - (gs - ar_end)
    avail_right = spec.down_available - (af_start - ge)
    if avail_left < params.anneal_min or avail_right < params.anneal_min:
        return None

    fwd = _pick_primer(
        target, spec.down_source, "forward", block, extra_f, af_start,
        spec.down_offset, avail_right, params,
    )
    # The reverse primer reads along the bottom strand, so both its tail and its extra
    # segment are the reverse complements of the target-strand sequence they encode.
    rev = _pick_primer(
        target, spec.up_source, "reverse", revcomp(block), revcomp(extra_r),
        ar_end + spec.up_shift, spec.up_offset, avail_left, params,
    )
    if fwd is None or rev is None:
        return None

    # Both primers must come out as exact circular slices of the desired target. Anything
    # else is a bug in the coordinate arithmetic, not an infeasible candidate.
    fwd_expected = circular_slice(target, a, a + block_len + len(extra_f) + len(fwd.anneal))
    rev_expected = revcomp(circular_slice(target, ar_end - len(rev.anneal), a + block_len))
    if fwd.sequence != fwd_expected:
        raise DesignError(
            f"internal error: forward primer {fwd.sequence} is not the target slice "
            f"{fwd_expected} (block at {a}, length {block_len})"
        )
    if rev.sequence != rev_expected:
        raise DesignError(
            f"internal error: reverse primer {rev.sequence} is not the target slice "
            f"{rev_expected} (block at {a}, length {block_len})"
        )

    overhang_tm = thermo.duplex_tm(block, revcomp(block), params.conditions)
    penalty, detail, warns = _score(fwd, rev, block, block_len, overhang_tm, params)

    fwd.name = f"J{spec.index + 1}_F"
    rev.name = f"J{spec.index + 1}_R"
    fwd.junction = rev.junction = spec.index
    return JunctionCandidate(
        edit_index=spec.index,
        block_start=a,
        block_len=block_len,
        overhang=block,
        forward=fwd,
        reverse=rev,
        overhang_tm=overhang_tm,
        penalty=penalty,
        penalty_detail=detail,
        warnings=warns,
    )


def _forward_distance(start: int, end: int, n: int) -> int:
    return (end - start) % n


def _pick_primer(
    target: str,
    source: str,
    direction: str,
    tail: str,
    extra: str,
    anchor: int,
    source_offset: int,
    available: int,
    params: DesignParams,
) -> Optional[Primer]:
    """Choose the annealing-region length that best hits the Tm target.

    `anchor` is the target coordinate where the annealing region starts (forward) or
    ends (reverse); `source_offset` converts target coordinates to coordinates on the source
    plasmid this primer amplifies; `available` is how much of that source's contribution is
    reachable before the next junction.

    Only primers with a single binding site in the source are ever returned. If nothing in
    the normal length range is unique -- which happens inside duplicated elements -- the
    region is extended up to `anneal_hard_max` and the Tm ceiling is dropped, because a long
    primer is the only way to reach out of a repeat.
    """
    normal = _scan_anneal_lengths(
        target, source, direction, tail, extra, anchor, source_offset,
        params.anneal_min, min(params.anneal_max, available), params, enforce_tm_ceiling=True,
    )
    if normal is not None:
        return normal
    extended = _scan_anneal_lengths(
        target, source, direction, tail, extra, anchor, source_offset,
        params.anneal_max + 1, min(params.anneal_hard_max, available), params,
        enforce_tm_ceiling=False,
    )
    if extended is not None:
        extended.warnings.append(
            f"annealing region extended to {len(extended.anneal)} nt to reach a sequence that "
            "is unique in the template (a repeated element lies closer in)"
        )
    return extended


def _scan_anneal_lengths(
    target: str,
    template: str,
    direction: str,
    tail: str,
    extra: str,
    anchor: int,
    template_offset: int,
    length_min: int,
    length_max: int,
    params: DesignParams,
    enforce_tm_ceiling: bool,
) -> Optional[Primer]:
    best: Optional[Primer] = None
    best_key: Optional[Tuple[float, ...]] = None

    for n in range(length_min, length_max + 1):
        if direction == "forward":
            tgt_region = circular_slice(target, anchor, anchor + n)
            tpl_region = circular_slice(template, anchor + template_offset, anchor + template_offset + n)
            anneal = tgt_region
        else:
            tgt_region = circular_slice(target, anchor - n, anchor)
            tpl_region = circular_slice(
                template, anchor - n + template_offset, anchor + template_offset
            )
            anneal = revcomp(tgt_region)
        if tgt_region != tpl_region:
            continue  # annealing region must match the template exactly

        anneal_tm = thermo.tm(anneal, params.conditions)
        if anneal_tm < params.tm_min - 4:
            continue
        if enforce_tm_ceiling and anneal_tm > params.tm_max + 4:
            continue

        sequence = tail + extra + anneal
        if len(sequence) > params.hard_max_primer_len:
            continue
        u_index = len(tail) - 1
        if sequence[u_index] != "T":
            return None  # should not happen: tails are constructed to end in T

        # A primer that binds more than once would amplify several products, so it is
        # rejected outright rather than merely penalised.
        hits = count_circular_occurrences(template, anneal, both_strands=True)
        if hits != 1:
            continue

        stats = thermo.OligoStats.of(sequence, params.conditions)
        warns: List[str] = []
        if stats.max_homopolymer > params.max_homopolymer:
            warns.append(f"{stats.max_homopolymer}-nt homopolymer run")
        if len(sequence) > params.soft_max_primer_len:
            warns.append(f"{len(sequence)} nt: order as an Ultramer / PAGE-purified oligo")
        anneal_gc = 100.0 * gc_fraction(anneal)
        if not params.gc_min <= anneal_gc <= params.gc_max:
            warns.append(f"annealing-region GC {anneal_gc:.0f}%")

        primer = Primer(
            name="", junction=-1, direction=direction, sequence=sequence, u_index=u_index,
            tail=tail, extra=extra, anneal=anneal, anneal_tm=anneal_tm, full_tm=stats.tm,
            gc_percent=stats.gc_percent, template_hits=hits, warnings=warns,
        )
        key = (
            abs(anneal_tm - params.tm_target),
            _three_prime_penalty(anneal),
            float(len(sequence)),
        )
        if best_key is None or key < best_key:
            best, best_key = primer, key

    return best


def _three_prime_penalty(anneal: str) -> float:
    """Prefer a 1-3 G/C clamp and avoid GC-rich or A/T-only 3' ends."""
    clamp = thermo.three_prime_gc_clamp(anneal, window=5)
    penalty = 0.0
    if clamp == 0:
        penalty += 2.0
    elif clamp >= 4:
        penalty += 1.0
    if anneal[-1] in "AT":
        penalty += 0.5
    if anneal[-2:] in ("GG", "CC", "GC", "CG") and clamp >= 4:
        penalty += 0.5
    return penalty


def _score(
    fwd: Primer,
    rev: Primer,
    block: str,
    block_len: int,
    overhang_tm: float,
    params: DesignParams,
) -> Tuple[float, Dict[str, float], List[str]]:
    detail: Dict[str, float] = {}
    warns: List[str] = []

    detail["tm_offset"] = abs(fwd.anneal_tm - params.tm_target) + abs(rev.anneal_tm - params.tm_target)
    detail["tm_mismatch"] = 2.0 * abs(fwd.anneal_tm - rev.anneal_tm)
    if abs(fwd.anneal_tm - rev.anneal_tm) > params.max_pair_tm_diff:
        warns.append(
            f"annealing Tm differs by {abs(fwd.anneal_tm - rev.anneal_tm):.1f} C between the pair"
        )

    # 8-10 nt overhangs are all routine, so only penalise drifting outside that plateau.
    detail["overhang_length"] = 1.5 * max(0, abs(block_len - params.overhang_preferred) - 1)
    detail["overhang_gc"] = _range_penalty(100.0 * gc_fraction(block), 30.0, 70.0, weight=0.15)

    # The overhang has to anneal during the post-USER incubation, so it needs some stability;
    # very long or very GC-rich overhangs buy nothing, hence a window rather than a floor.
    # (The two primer tails are complementary by construction -- that duplex *is* the
    # overhang -- and being a 5'-5' overlap it cannot be extended, so it is not a dimer risk.)
    detail["overhang_stability"] = _range_penalty(overhang_tm, 20.0, 45.0, weight=0.4)
    if overhang_tm < 15.0:
        warns.append(
            f"overhang duplex Tm is only {overhang_tm:.0f} C, so anneal at 10 C rather than "
            "room temperature after the USER digest"
        )
    detail["three_prime"] = _three_prime_penalty(fwd.anneal) + _three_prime_penalty(rev.anneal)
    # Long oligos cost more and are synthesised less accurately, and crossing the standard
    # 60 nt limit forces an Ultramer order, so the penalty steepens there. Charging both
    # primers separately also pushes new sequence to be split evenly between the pair.
    detail["length"] = sum(
        0.15 * max(0, primer.length - 45)
        + 0.6 * max(0, primer.length - params.soft_max_primer_len)
        for primer in (fwd, rev)
    )
    # Total primer length is fixed once the junction has to carry a given amount of new
    # sequence, so the penalty above is blind to how that sequence is divided. Charging the
    # imbalance separately makes the pair converge on two medium oligos rather than one very
    # long one -- cheaper to order and synthesised more accurately.
    detail["length_imbalance"] = 0.15 * abs(fwd.length - rev.length)
    detail["specificity"] = 50.0 * ((fwd.template_hits != 1) + (rev.template_hits != 1))
    detail["homopolymer"] = 2.0 * (
        max(0, thermo.longest_homopolymer(fwd.sequence) - params.max_homopolymer)
        + max(0, thermo.longest_homopolymer(rev.sequence) - params.max_homopolymer)
    )
    warns.extend(fwd.warnings)
    warns.extend(rev.warnings)
    return sum(detail.values()), detail, sorted(set(warns))


def _range_penalty(value: float, low: float, high: float, weight: float = 1.0) -> float:
    if value < low:
        return weight * (low - value)
    if value > high:
        return weight * (value - high)
    return 0.0


# --------------------------------------------------------------------------------------
# Multi-junction selection
# --------------------------------------------------------------------------------------


def _select_compatible(
    candidate_sets: List[List[JunctionCandidate]],
    params: DesignParams,
) -> List[JunctionCandidate]:
    """Beam search for a set of junctions with mutually orthogonal overhangs."""
    if len(candidate_sets) == 1:
        return [candidate_sets[0][0]]

    beam: List[Tuple[float, List[JunctionCandidate]]] = [(0.0, [])]
    for candidates in candidate_sets:
        nxt: List[Tuple[float, List[JunctionCandidate]]] = []
        for score, chosen in beam:
            for cand in candidates:
                if not _orthogonal(cand, chosen, params):
                    continue
                nxt.append((score + cand.penalty, chosen + [cand]))
        if not nxt:  # nothing orthogonal: fall back to best-scoring and warn later
            for score, chosen in beam:
                nxt.append((score + candidates[0].penalty, chosen + [candidates[0]]))
        nxt.sort(key=lambda item: item[0])
        beam = nxt[: params.beam_width]
    return beam[0][1]


def _orthogonal(
    cand: JunctionCandidate,
    chosen: Sequence[JunctionCandidate],
    params: DesignParams,
) -> bool:
    for other in chosen:
        if cand.overhang == other.overhang:
            return False
        cross = max(
            thermo.duplex_tm(cand.overhang, revcomp(other.overhang), params.conditions),
            thermo.duplex_tm(other.overhang, revcomp(cand.overhang), params.conditions),
        )
        cognate = min(cand.overhang_tm, other.overhang_tm)
        if cross > cognate - 8.0:
            return False
    return True


def _check_overhang_orthogonality(junctions: Sequence[Junction], params: DesignParams) -> List[str]:
    warnings: List[str] = []
    for i, a in enumerate(junctions):
        for b in junctions[i + 1:]:
            cross = max(
                thermo.duplex_tm(a.overhang, revcomp(b.overhang), params.conditions),
                thermo.duplex_tm(b.overhang, revcomp(a.overhang), params.conditions),
            )
            if cross > min(a.overhang_tm, b.overhang_tm) - 8.0:
                warnings.append(
                    f"Junction {a.index + 1} and {b.index + 1} overhangs cross-anneal "
                    f"({cross:.0f} C vs cognate {min(a.overhang_tm, b.overhang_tm):.0f} C): "
                    "mis-assembly is possible."
                )
            if a.overhang == b.overhang:
                warnings.append(
                    f"Junctions {a.index + 1} and {b.index + 1} share an identical overhang; "
                    "assembly will not be directional."
                )
    return warnings


def _locate_unique(template: str, motif: str) -> int:
    from .sequences import find_all

    extended = template + template[:max(0, len(motif) - 1)]
    hits = find_all(extended, motif)
    if not hits:
        raise DesignError(f"Annealing region {motif} not found in the template")
    return hits[0] % len(template)
