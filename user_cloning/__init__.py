"""USER cloning primer design pipeline.

Given a spreadsheet of named circular plasmid sequences, work out how to build one
plasmid from another by uracil-excision (USER) cloning, and emit verified, orderable
primers.
"""

from .design import DesignError, DesignParams, Junction, Primer
from .pipeline import DesignResult, design_one, resolve_records
from .plasmid_diff import EditSite, PlasmidComparison, compare_plasmids
from .sequences import SeqRecord, read_sequence_table, revcomp

__all__ = [
    "DesignError",
    "DesignParams",
    "DesignResult",
    "EditSite",
    "Junction",
    "PlasmidComparison",
    "Primer",
    "SeqRecord",
    "compare_plasmids",
    "design_one",
    "read_sequence_table",
    "resolve_records",
    "revcomp",
]

__version__ = "1.0.0"
