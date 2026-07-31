#!/usr/bin/env python3
"""Design USER cloning primers from a spreadsheet of plasmid sequences.

Examples
--------
List what is in a spreadsheet:
    python3 design_user_primers.py USER_design_1.xlsx --list

Design one construct (the case this pipeline was built for):
    python3 design_user_primers.py USER_design_1.xlsx --template pLL057P --target pLL082P

A whole batch, all cloned from the same backbone:
    python3 design_user_primers.py 2026-08_designs.xlsx --all-from pLL057P

Explicit pairs:
    python3 design_user_primers.py in.xlsx --pairs pLL057P:pLL082P pLL057P:pLL083P

Every run writes into designs/<date>/, so repeat batches never overwrite each other.
"""

from __future__ import annotations

import argparse
import os
import sys
from dataclasses import replace
from datetime import date
from typing import List, Optional, Sequence, Tuple

from user_cloning.design import DesignError, DesignParams
from user_cloning.pipeline import (
    STATUS_NO_DIFFERENCE,
    STATUS_OK,
    DesignResult,
    design_assembly,
    resolve_records,
)
from user_cloning.report import allocate_batch, write_all, write_batch_tables
from user_cloning.sequences import SeqRecord, read_sequence_table, repeat_blocks


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="design_user_primers.py",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("input", help="input .xlsx / .csv / .tsv / .fasta of named plasmid sequences")
    p.add_argument("--sheet", help="worksheet name (default: first sheet)")
    p.add_argument("--list", action="store_true", dest="list_only",
                   help="list the sequences in the input file and exit")

    sel = p.add_argument_group("what to design")
    sel.add_argument("--template", help="name of the plasmid to amplify from")
    sel.add_argument("--sources", nargs="+", metavar="NAME",
                     help="two or more source plasmids to assemble the target from")
    sel.add_argument("--target", help="name of the plasmid to build")
    sel.add_argument("--pairs", nargs="+", metavar="TEMPLATE:TARGET",
                     help="one or more explicit template:target pairs")
    sel.add_argument("--all-from", metavar="TEMPLATE",
                     help="design every other sequence in the file from this template")

    out = p.add_argument_group("output")
    out.add_argument("--outdir", default="designs", help="output root (default: designs)")
    out.add_argument("--date", default=date.today().isoformat(),
                     help="date stamp for this batch, YYYY-MM-DD (default: today)")
    out.add_argument("--batch", metavar="LABEL",
                     help="label for this batch, e.g. b3 or a project name. Several batches can be "
                          "requested on the same day, so results are keyed on date AND batch. "
                          "Defaults to the next free b1, b2, ... for the date. Existing "
                          "results are never overwritten or deleted.")

    tune = p.add_argument_group("design constraints")
    tune.add_argument("--overhang-min", type=int, default=DesignParams.overhang_min)
    tune.add_argument("--overhang-max", type=int, default=DesignParams.overhang_max)
    tune.add_argument("--anneal-min", type=int, default=DesignParams.anneal_min)
    tune.add_argument("--anneal-max", type=int, default=DesignParams.anneal_max)
    tune.add_argument("--tm-target", type=float, default=DesignParams.tm_target,
                      help="target Tm (C) for the template-binding region")
    tune.add_argument("--junction-slack", type=int, default=DesignParams.junction_window_slack,
                      help="how far from an edit the overhang block may be placed (nt)")
    tune.add_argument("--max-primer-len", type=int, default=DesignParams.soft_max_primer_len,
                      help="length above which a primer is flagged as an Ultramer order")
    return p


def params_from_args(args: argparse.Namespace) -> DesignParams:
    if args.overhang_min > args.overhang_max:
        raise SystemExit("--overhang-min must not exceed --overhang-max")
    if args.anneal_min > args.anneal_max:
        raise SystemExit("--anneal-min must not exceed --anneal-max")
    return replace(
        DesignParams(),
        overhang_min=args.overhang_min,
        overhang_max=args.overhang_max,
        overhang_preferred=min(max(9, args.overhang_min), args.overhang_max),
        anneal_min=args.anneal_min,
        anneal_max=args.anneal_max,
        tm_target=args.tm_target,
        junction_window_slack=args.junction_slack,
        soft_max_primer_len=args.max_primer_len,
    )


def resolve_pairs(args: argparse.Namespace, records: Sequence[SeqRecord]) -> List[Tuple[str, str]]:
    pairs: List[Tuple[str, str]] = []
    if args.pairs:
        for item in args.pairs:
            sep = ":" if ":" in item else (">" if ">" in item else None)
            if sep is None:
                raise SystemExit(f"--pairs entry {item!r} must look like TEMPLATE:TARGET")
            tpl, tgt = item.split(sep, 1)
            pairs.append((tpl.strip(), tgt.strip()))
    if args.all_from:
        for rec in records:
            if rec.name.lower() != args.all_from.lower():
                pairs.append((args.all_from, rec.name))
    if args.sources:
        if not args.target:
            raise SystemExit("--sources requires --target")
        if len(args.sources) < 2:
            raise SystemExit("--sources needs at least two plasmids; use --template for one")
    if args.template and args.target:
        pairs.append((args.template, args.target))
    elif args.template and not args.sources:
        raise SystemExit("--template and --target must be given together")
    elif args.target and not args.sources and not args.template:
        raise SystemExit("--target needs --template or --sources")

    if not pairs and not args.sources:
        raise SystemExit(
            "Nothing to design. Give --template/--target, --pairs, or --all-from "
            "(use --list to see what is in the file)."
        )
    seen, unique = set(), []
    for pair in pairs:
        if pair not in seen:
            seen.add(pair)
            unique.append(pair)
    return unique


def resolve_assembly(
    records: Sequence[SeqRecord],
    source_names: Sequence[str],
    target_name: str,
) -> Tuple[List[SeqRecord], SeqRecord]:
    """Look up several sources plus a target, tolerating case and punctuation differences."""
    index = {_key(rec.name): rec for rec in records}
    missing = [n for n in list(source_names) + [target_name] if _key(n) not in index]
    if missing:
        raise DesignError(
            f"Sequence(s) {', '.join(missing)} not found in the input file. "
            f"Available: {', '.join(rec.name for rec in records)}"
        )
    return [index[_key(n)] for n in source_names], index[_key(target_name)]


def _key(name: str) -> str:
    """Names get typed with and without hyphens, so compare on alphanumerics only."""
    return "".join(c for c in name.lower() if c.isalnum())


def print_listing(records: Sequence[SeqRecord], path: str) -> None:
    print(f"{len(records)} sequence(s) in {os.path.basename(path)}:")
    for rec in records:
        blocks = repeat_blocks(rec.seq)
        repetitive = sum(b.length for b in blocks)
        note = ""
        if blocks:
            note = (f"  ({len(blocks)} repeated block(s), {repetitive:,} nt = "
                    f"{100.0 * repetitive / rec.length:.0f}% not uniquely primable)")
        print(f"  row {rec.source_row:>3}  {rec.name:<24} {rec.length:>8,} bp{note}")
        for b in blocks:
            print(f"           repeat {b.start + 1:,}..{b.end + 1:,} "
                  f"({b.length:,} nt, {b.copies} copies)")
    for group in _identical_groups(records):
        print(f"  ! identical sequences: {', '.join(group)}")


def _identical_groups(records: Sequence[SeqRecord]) -> List[List[str]]:
    by_seq: dict = {}
    for rec in records:
        by_seq.setdefault(rec.seq, []).append(rec.name)
    return [names for names in by_seq.values() if len(names) > 1]


def _show(seq: str, limit: int = 40) -> str:
    if not seq:
        return "(nothing)"
    return seq if len(seq) <= limit else f"{seq[:18]}...{seq[-18:]} ({len(seq)} nt)"


def report_result(result: DesignResult, paths: dict) -> None:
    header = f"{result.target_name} from {result.template_name}  [{result.design_id}]"
    print()
    print(header)
    print("-" * len(header))
    if result.status == STATUS_NO_DIFFERENCE:
        for msg in result.messages:
            print(f"  SKIPPED: {msg}")
    elif result.status == STATUS_OK:
        comparison = result.comparison
        n_edits = len(comparison.edits) if comparison else 0
        print(f"  {n_edits} difference site(s) -> {len(result.junctions)} junction(s), "
              f"{len(result.fragments)} PCR fragment(s), {len(result.primers)} primer(s)")
        for e in comparison.edits if comparison else []:
            # Report positions on the sequences as supplied, not on the internal frame.
            start = comparison.to_template_input(e.tpl_start)
            where = (f"at {result.backbone_name} {start + 1}..{start + e.removed_length}"
                     if e.removed_length
                     else f"between {result.backbone_name} {start} and {start + 1}")
            print(f"    - {e.kind} {where}: {_show(e.tpl_seq)} -> {_show(e.tgt_seq)}")
        for primer in result.primers:
            print(f"  {primer.name:<8} {primer.order_sequence}")
            print(f"           {primer.length} nt, dU at {primer.u_index + 1}, "
                  f"annealing Tm {primer.anneal_tm:.1f} C")
        source_of = {f.name: f.source_name for f in result.fragments}
        for p in result.protocols:
            print(f"  {p.fragment} from {source_of.get(p.fragment, '?')}: "
                  f"{p.product_length:,} bp product, Ta {p.annealing_temp_C} C, "
                  f"extension {p.extension_seconds} s")
        print("  VERIFIED: in-silico assembly reproduces the target exactly")
    else:
        for msg in result.messages:
            print(f"  FAILED: {msg}")
    for warning in dict.fromkeys(result.warnings):
        print(f"  warning: {warning}")
    for kind, path in paths.items():
        print(f"  wrote {kind}: {path}")


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_arg_parser().parse_args(argv)
    try:
        records = read_sequence_table(args.input, args.sheet)
    except (OSError, ValueError) as exc:
        print(f"error reading {args.input}: {exc}", file=sys.stderr)
        return 2

    if args.list_only:
        print_listing(records, args.input)
        return 0

    params = params_from_args(args)
    pairs = resolve_pairs(args, records)

    try:
        batch = allocate_batch(args.outdir, args.date, args.batch)
    except (FileExistsError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    print(f"batch {batch} ({args.date}) -> {os.path.join(args.outdir, args.date, batch)}")

    results: List[DesignResult] = []
    exit_code = 0

    if args.sources:
        try:
            sources, target = resolve_assembly(records, args.sources, args.target)
            result = design_assembly(sources, target, args.date, params,
                                     design_index=1, batch=batch)
        except DesignError as exc:
            print(f"\n{args.target} from {', '.join(args.sources)}\n  FAILED: {exc}",
                  file=sys.stderr)
            return 1
        paths = write_all(args.outdir, result, args.input)
        report_result(result, paths)
        results.append(result)
        if not result.verified:
            exit_code = 3 if result.status == STATUS_NO_DIFFERENCE else 1
        pairs = []

    for index, (tpl_name, tgt_name) in enumerate(pairs, start=1):
        try:
            template, target = resolve_records(records, tpl_name, tgt_name)
            result = design_assembly([template], target, args.date, params,
                                     design_index=index, batch=batch)
        except DesignError as exc:
            print(f"\n{tgt_name} from {tpl_name}\n  FAILED: {exc}", file=sys.stderr)
            exit_code = 1
            continue
        paths = write_all(args.outdir, result, args.input)
        report_result(result, paths)
        results.append(result)
        if result.status == STATUS_NO_DIFFERENCE:
            exit_code = max(exit_code, 3)
        elif not result.verified:
            exit_code = max(exit_code, 1)

    if results:
        for kind, path in write_batch_tables(args.outdir, args.date, batch, results).items():
            print(f"\n{kind}: {path}" if kind == "batch_summary" else f"{kind}: {path}")
    return exit_code


if __name__ == "__main__":
    sys.exit(main())
