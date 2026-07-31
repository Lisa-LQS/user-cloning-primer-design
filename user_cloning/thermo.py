"""Oligo thermodynamics via primer3-py (SantaLucia nearest-neighbour + salt correction).

Deoxyuridine is treated as thymidine for every calculation, which is the standard
approximation: dU:dA and dT:dA nearest-neighbour parameters are close enough that the
difference is well inside the error of the model.

Only duplex melting temperatures are computed here. Hairpin and self-dimer predictions are
deliberately not used anywhere in the pipeline: for USER cloning they flag far more primers
than actually fail, so they were dropped rather than left to distort primer choice.
"""

from __future__ import annotations

from dataclasses import dataclass

import primer3

from .sequences import gc_fraction


@dataclass(frozen=True)
class PcrConditions:
    """Reaction conditions used for every Tm calculation.

    Defaults match a standard 50 uL Phusion U / Q5U reaction: 1x HF buffer is ~50 mM
    monovalent salt, 2 mM Mg2+ (1.5 mM free after dNTP chelation), 200 uM each dNTP,
    500 nM each primer.
    """

    monovalent_mM: float = 50.0
    divalent_mM: float = 1.5
    dntp_mM: float = 0.8
    primer_nM: float = 500.0

    def as_primer3_kwargs(self) -> dict:
        return dict(
            mv_conc=self.monovalent_mM,
            dv_conc=self.divalent_mM,
            dntp_conc=self.dntp_mM,
            dna_conc=self.primer_nM,
        )


DEFAULT_CONDITIONS = PcrConditions()

MAX_DUPLEX_LENGTH = 60
"""primer3 refuses two-sequence duplex calculations above this length."""


def _dna(seq: str) -> str:
    return seq.upper().replace("U", "T")


def _duplex_window(seq: str) -> str:
    """Clamp a sequence to primer3's duplex limit, keeping the 3' end."""
    return _dna(seq)[-MAX_DUPLEX_LENGTH:]


def tm(seq: str, conditions: PcrConditions = DEFAULT_CONDITIONS) -> float:
    """Melting temperature (deg C) of an oligo:template duplex."""
    seq = _dna(seq)
    if len(seq) < 2:
        return 0.0
    return primer3.calc_tm(seq, **conditions.as_primer3_kwargs())


def duplex_tm(seq_a: str, seq_b: str, conditions: PcrConditions = DEFAULT_CONDITIONS) -> float:
    """Tm of the best duplex between two oligos (used for USER overhang annealing)."""
    result = primer3.calc_heterodimer(
        _duplex_window(seq_a), _duplex_window(seq_b), **conditions.as_primer3_kwargs()
    )
    return result.tm if result.structure_found else -99.0


def longest_homopolymer(seq: str) -> int:
    best = run = 1
    for i in range(1, len(seq)):
        run = run + 1 if seq[i] == seq[i - 1] else 1
        best = max(best, run)
    return best if seq else 0


def three_prime_gc_clamp(seq: str, window: int = 5) -> int:
    """Number of G/C in the last `window` bases."""
    tail = seq[-window:].upper()
    return sum(1 for c in tail if c in "GC")


@dataclass
class OligoStats:
    """Summary statistics for one oligo (or one sub-region of one)."""

    seq: str
    length: int
    gc_percent: float
    tm: float
    max_homopolymer: int

    @classmethod
    def of(cls, seq: str, conditions: PcrConditions = DEFAULT_CONDITIONS) -> "OligoStats":
        dna = _dna(seq)
        return cls(
            seq=seq,
            length=len(seq),
            gc_percent=100.0 * gc_fraction(dna),
            tm=tm(dna, conditions),
            max_homopolymer=longest_homopolymer(dna),
        )
