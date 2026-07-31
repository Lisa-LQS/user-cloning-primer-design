"""Sequence helpers and spreadsheet/FASTA I/O for the USER cloning pipeline."""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

DNA_LETTERS = set("ACGT")
IUPAC_LETTERS = set("ACGTNRYSWKMBDHV")
_COMPLEMENT = str.maketrans("ACGTNRYSWKMBDHVacgtnryswkmbdhvUu", "TGCANYRSWMKVHDBtgcanyrswmkvhdbAa")


def clean_seq(raw: str) -> str:
    """Uppercase a sequence and strip whitespace, digits and FASTA-style decorations."""
    return re.sub(r"[^A-Za-z]", "", str(raw)).upper()


def is_dna(seq: str, min_purity: float = 0.9) -> bool:
    if not seq:
        return False
    good = sum(1 for c in seq if c in IUPAC_LETTERS)
    acgt = sum(1 for c in seq if c in DNA_LETTERS)
    return good / len(seq) >= min_purity and acgt / len(seq) >= min_purity


def revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def gc_fraction(seq: str) -> float:
    if not seq:
        return 0.0
    return sum(1 for c in seq.upper() if c in "GCS") / len(seq)


def rotate(seq: str, n: int) -> str:
    """Left-rotate a circular sequence by n (n may be negative or > len)."""
    if not seq:
        return seq
    n %= len(seq)
    return seq[n:] + seq[:n]


def circular_slice(seq: str, start: int, end: int) -> str:
    """Slice a circular sequence; `end` may wrap past the origin. Half-open [start, end)."""
    n = len(seq)
    if n == 0:
        return ""
    length = end - start
    if length < 0:
        raise ValueError(f"circular_slice needs end >= start, got {start}..{end}")
    if length > n:
        raise ValueError(f"circular_slice of {length} nt exceeds sequence length {n}")
    start %= n
    if start + length <= n:
        return seq[start:start + length]
    return seq[start:] + seq[:start + length - n]


def find_all(haystack: str, needle: str) -> List[int]:
    """All (possibly overlapping) start positions of `needle` in `haystack`."""
    if not needle:
        return []
    hits, i = [], haystack.find(needle)
    while i != -1:
        hits.append(i)
        i = haystack.find(needle, i + 1)
    return hits


def count_circular_occurrences(seq: str, motif: str, both_strands: bool = True) -> int:
    """Occurrences of `motif` in a circular sequence, optionally on both strands."""
    if not motif or len(motif) > len(seq):
        return 0
    extended = seq + seq[:len(motif) - 1]
    n = len(find_all(extended, motif))
    if both_strands:
        n += len(find_all(extended, revcomp(motif)))
    return n


def circular_equal(a: str, b: str) -> Tuple[bool, Optional[int]]:
    """Are two sequences identical as circles? Returns (equal, rotation of b matching a)."""
    if len(a) != len(b) or not a:
        return (len(a) == len(b), 0 if len(a) == len(b) else None)
    idx = (b + b).find(a)
    if idx == -1 or idx >= len(b):
        return (False, None)
    return (True, idx)


def kmer_index(seq: str, k: int, circular: bool = True) -> Dict[str, List[int]]:
    space = seq + seq[:k - 1] if circular and len(seq) >= k else seq
    idx: Dict[str, List[int]] = {}
    for i in range(len(space) - k + 1):
        idx.setdefault(space[i:i + k], []).append(i % len(seq))
    return idx


@dataclass(frozen=True)
class RepeatBlock:
    """A stretch of sequence that is not unique in the plasmid.

    Primers cannot be made specific inside one of these, so USER junctions have to be
    placed outside them.
    """

    start: int          # 0-based, inclusive
    end: int            # 0-based, inclusive
    copies: int         # highest copy number of any k-mer in the block

    @property
    def length(self) -> int:
        return self.end - self.start + 1

    def overlaps(self, start: int, end: int) -> bool:
        """Does this block overlap the half-open interval [start, end)?"""
        return start <= self.end and end > self.start


def repeat_blocks(seq: str, k: int = 40, circular: bool = True) -> List[RepeatBlock]:
    """Find maximal stretches covered by k-mers that occur more than once.

    Only exact, same-strand repeats are reported, which is what matters for priming: a
    duplicated promoter or LTR will bind the same primer twice.
    """
    if len(seq) < k:
        return []
    index = kmer_index(seq, k, circular=circular)
    flagged = sorted({p for positions in index.values() if len(positions) > 1 for p in positions})
    if not flagged:
        return []
    counts = {p: len(positions) for positions in index.values()
              if len(positions) > 1 for p in positions}

    blocks: List[RepeatBlock] = []
    start, end, copies = flagged[0], flagged[0] + k - 1, counts[flagged[0]]
    for p in flagged[1:]:
        if p <= end + 1:
            end = max(end, p + k - 1)
            copies = max(copies, counts[p])
        else:
            blocks.append(RepeatBlock(start, min(end, len(seq) - 1), copies))
            start, end, copies = p, p + k - 1, counts[p]
    blocks.append(RepeatBlock(start, min(end, len(seq) - 1), copies))
    return blocks


def unique_stretches(seq: str, blocks: Sequence[RepeatBlock]) -> List[Tuple[int, int]]:
    """The complement of `blocks`: 0-based half-open intervals that are unique."""
    stretches: List[Tuple[int, int]] = []
    cursor = 0
    for block in sorted(blocks, key=lambda b: b.start):
        if block.start > cursor:
            stretches.append((cursor, block.start))
        cursor = max(cursor, block.end + 1)
    if cursor < len(seq):
        stretches.append((cursor, len(seq)))
    return stretches


@dataclass
class SeqRecord:
    """A named plasmid sequence as read from the input spreadsheet."""

    name: str
    seq: str
    source_row: int

    @property
    def length(self) -> int:
        return len(self.seq)


def read_sequence_table(path: str, sheet: Optional[str] = None) -> List[SeqRecord]:
    """Read (name, sequence) pairs from an .xlsx/.csv/.tsv/.fasta file.

    Column order is not assumed: the sequence column is the one that looks like DNA.
    A header row, if present, is skipped automatically.
    """
    lower = str(path).lower()
    if lower.endswith((".fa", ".fasta", ".fna", ".seq", ".txt")):
        return _read_fasta_records(path)
    rows = _read_rows_xlsx(path, sheet) if lower.endswith((".xlsx", ".xlsm")) else _read_rows_delimited(path)
    return _rows_to_records(rows)


def _read_rows_xlsx(path: str, sheet: Optional[str]) -> List[List[Optional[str]]]:
    import openpyxl

    wb = openpyxl.load_workbook(path, data_only=True, read_only=True)
    ws = wb[sheet] if sheet else wb.worksheets[0]
    rows = [[None if c is None else str(c) for c in row] for row in ws.iter_rows(values_only=True)]
    wb.close()
    return rows


def _read_rows_delimited(path: str) -> List[List[Optional[str]]]:
    import csv

    delim = "\t" if str(path).lower().endswith((".tsv", ".tab")) else ","
    with open(path, newline="") as fh:
        return [list(row) for row in csv.reader(fh, delimiter=delim)]


def _rows_to_records(rows: List[List[Optional[str]]]) -> List[SeqRecord]:
    records: List[SeqRecord] = []
    for row_no, row in enumerate(rows, start=1):
        cells = [("" if c is None else str(c).strip()) for c in row]
        seq_col = None
        for i, cell in enumerate(cells):
            candidate = clean_seq(cell)
            if len(candidate) >= 50 and is_dna(candidate):
                seq_col = i
                break
        if seq_col is None:
            continue  # header row, blank row, or notes
        name_cells = [c for i, c in enumerate(cells) if i != seq_col and c]
        name = name_cells[0] if name_cells else f"row{row_no}"
        records.append(SeqRecord(name=name, seq=clean_seq(cells[seq_col]), source_row=row_no))
    if not records:
        raise ValueError("No rows with a DNA-like sequence column were found")
    _check_unique_names(records)
    return records


def _read_fasta_records(path: str) -> List[SeqRecord]:
    records: List[SeqRecord] = []
    name, chunks, first_line = None, [], 0
    with open(path) as fh:
        for line_no, line in enumerate(fh, start=1):
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    records.append(SeqRecord(name, clean_seq("".join(chunks)), first_line))
                name, chunks, first_line = line[1:].split()[0] or f"seq{line_no}", [], line_no
            elif line:
                chunks.append(line)
    if name is not None:
        records.append(SeqRecord(name, clean_seq("".join(chunks)), first_line))
    if not records:
        raise ValueError(f"No FASTA records found in {path}")
    _check_unique_names(records)
    return records


def _check_unique_names(records: List[SeqRecord]) -> None:
    seen: Dict[str, int] = {}
    for rec in records:
        if rec.name in seen:
            raise ValueError(
                f"Duplicate sequence name {rec.name!r} (rows {seen[rec.name]} and {rec.source_row}); "
                "names must be unique so template/target can be resolved unambiguously"
            )
        seen[rec.name] = rec.source_row


def write_fasta(path: str, name: str, seq: str, description: str = "", width: int = 70) -> None:
    header = f">{name} {description}".rstrip()
    with open(path, "w") as fh:
        fh.write(header + "\n")
        for i in range(0, len(seq), width):
            fh.write(seq[i:i + width] + "\n")
