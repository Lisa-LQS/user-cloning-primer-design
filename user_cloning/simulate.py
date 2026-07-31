"""In-silico PCR + USER digestion + assembly, used to verify every design.

Nothing is reported to the user unless the simulated assembly reproduces the requested
target plasmid base for base (as a circle).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Sequence, Tuple

from .design import Fragment, Junction
from .sequences import circular_equal, circular_slice, find_all, revcomp


@dataclass
class SimulatedProduct:
    """One PCR product, before and after USER treatment."""

    fragment: str
    top_strand: str            # full PCR product, top strand
    length: int
    top_after_user: str        # top strand after the 5' flap is excised
    bottom_after_user: str     # bottom strand after its 5' flap is excised
    left_overhang: str         # 3' overhang on the bottom strand, read 5'->3'
    right_overhang: str        # 3' overhang on the top strand, read 5'->3'


@dataclass
class Verification:
    """Outcome of simulating the whole design."""

    ok: bool
    products: List[SimulatedProduct] = field(default_factory=list)
    assembled: str = ""
    assembled_length: int = 0
    rotation_vs_target: Optional[int] = None
    errors: List[str] = field(default_factory=list)
    checks: List[Tuple[str, bool, str]] = field(default_factory=list)

    def record(self, name: str, passed: bool, detail: str = "") -> bool:
        self.checks.append((name, passed, detail))
        if not passed:
            self.errors.append(f"{name}: {detail}" if detail else name)
        return passed


def simulate_pcr(template: str, forward: str, reverse: str, anneal_f: str, anneal_r: str) -> str:
    """Return the top strand of the product amplified from a circular template."""
    f_site = _unique_site(template, anneal_f, "forward annealing region")
    r_site = _unique_site(template, revcomp(anneal_r), "reverse annealing region")
    f_end = (f_site + len(anneal_f)) % len(template)
    span_len = (r_site - f_end) % len(template)
    interior = circular_slice(template, f_end, f_end + span_len)
    return forward + interior + revcomp(reverse)


def user_digest(top_strand: str, forward_u_index: int, reverse_u_index: int) -> Tuple[str, str]:
    """Excise both 5' flaps. Returns (top strand, bottom strand), each read 5'->3'."""
    bottom = revcomp(top_strand)
    return top_strand[forward_u_index + 1:], bottom[reverse_u_index + 1:]


def verify_design(
    target: str,
    junctions: Sequence[Junction],
    fragments: Sequence[Fragment],
) -> Verification:
    """Simulate the full workflow and check it reconstitutes the target plasmid.

    Each fragment is amplified off its own source plasmid, so multi-template assemblies are
    handled the same way as single-template ones.
    """
    v = Verification(ok=True)
    if not junctions:
        v.record("has junctions", False, "no junctions were designed")
        v.ok = False
        return v

    for frag in fragments:
        fwd, rev = frag.forward, frag.reverse
        try:
            top = simulate_pcr(
                frag.source_seq, fwd.sequence, rev.sequence, fwd.anneal, rev.anneal
            )
        except ValueError as exc:
            v.record(f"PCR of {frag.name} from {frag.source_name}", False, str(exc))
            v.ok = False
            return v
        top_cut, bottom_cut = user_digest(top, fwd.u_index, rev.u_index)
        # The flap removed from the bottom strand was paired with the top strand's 3' end,
        # and vice versa, so each excision exposes a 3' overhang on the opposite strand.
        left_overhang = bottom_cut[-(fwd.u_index + 1):]
        right_overhang = top_cut[-(rev.u_index + 1):]
        v.products.append(
            SimulatedProduct(
                fragment=frag.name, top_strand=top, length=len(top),
                top_after_user=top_cut, bottom_after_user=bottom_cut,
                left_overhang=left_overhang, right_overhang=right_overhang,
            )
        )
        v.record(
            f"{frag.name} amplifies from {frag.source_name} at the predicted length",
            len(top) == frag.expected_length,
            f"simulated {len(top)} nt vs predicted {frag.expected_length} nt",
        )

    # Each junction joins the right end of one fragment to the left end of the next.
    n = len(fragments)
    for i, frag in enumerate(fragments):
        nxt = fragments[(i + 1) % n]
        right = v.products[i].right_overhang
        left = v.products[(i + 1) % n].left_overhang
        where = (
            f"{frag.name} self-circularises" if nxt is frag
            else f"{frag.name} joins {nxt.name}"
        )
        v.record(
            f"junction {frag.right_junction + 1} overhangs are complementary",
            right == revcomp(left),
            f"{where}: right-end overhang {right} vs left-end overhang {left}",
        )
        junction = junctions[frag.right_junction]
        v.record(
            f"junction {frag.right_junction + 1} overhang matches design",
            right == junction.top_overhang,
            f"simulated {right} vs designed {junction.top_overhang}",
        )

    # Annealing the fragments in order regenerates each strand of the circle intact,
    # nicked once per junction, so concatenating the post-USER top strands gives the plasmid.
    assembled = "".join(p.top_after_user for p in v.products)
    v.assembled = assembled
    v.assembled_length = len(assembled)
    v.record(
        "assembled length equals target length",
        len(assembled) == len(target),
        f"assembled {len(assembled)} nt vs target {len(target)} nt",
    )

    equal, rotation = circular_equal(assembled, target)
    v.rotation_vs_target = rotation
    v.record(
        "assembled plasmid is identical to the requested target (as a circle)",
        equal,
        "sequences differ" if not equal else f"target rotation +{rotation}",
    )

    bottom_join = "".join(p.bottom_after_user for p in reversed(v.products))
    bottom_equal, _ = circular_equal(bottom_join, revcomp(target))
    v.record(
        "bottom strand also reconstitutes the target",
        bottom_equal,
        "bottom strands do not close the circle" if not bottom_equal else "",
    )

    v.ok = all(passed for _, passed, _ in v.checks)
    return v


def _unique_site(template: str, motif: str, label: str) -> int:
    extended = template + template[:max(0, len(motif) - 1)]
    hits = [h % len(template) for h in find_all(extended, motif)]
    unique = sorted(set(hits))
    if not unique:
        raise ValueError(f"{label} {motif} does not occur in the template")
    if len(unique) > 1:
        raise ValueError(f"{label} {motif} occurs {len(unique)}x in the template at {unique}")
    return unique[0]
