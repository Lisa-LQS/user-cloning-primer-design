"""Write date-stamped design outputs: primer table, order form, report, FASTA, JSON."""

from __future__ import annotations

import csv
import json
import os
from dataclasses import asdict
from typing import Dict, List, Optional, Sequence

from .design import U_MARK, Junction
from .pipeline import (
    STATUS_NO_DIFFERENCE,
    STATUS_OK,
    DesignResult,
    predicted_plasmid,
)
from .sequences import revcomp, unique_stretches, write_fasta

PRIMER_COLUMNS = [
    "design_id", "design_date", "batch", "target", "template", "primer_name", "junction",
    "direction", "order_sequence", "sequence_with_U", "plain_sequence", "length_nt",
    "u_position_from_5prime",
    "user_junction_sequence", "user_junction_length_nt", "user_junction_duplex_tm_C",
    "three_prime_overhang_this_end", "amplified_from",
    "tail_removed_by_USER", "non_templated_extra", "annealing_region", "annealing_length_nt",
    "annealing_tm_C", "full_length_tm_C", "gc_percent",
    "template_binding_sites", "warnings",
]

USER_JUNCTION_COLUMN_NOTE = (
    "user_junction_sequence is the single-stranded 3' overhang that USER exposes at this "
    "junction, written as the top strand reads it 5'->3'. Both primers of a junction share "
    "it: the forward primer's fragment presents its reverse complement, and the two anneal."
)

SUMMARY_COLUMNS = [
    "design_id", "design_date", "batch", "target", "template", "status", "verified",
    "edits", "junctions", "fragments", "primers", "warnings", "notes",
]


def batch_dir(outdir: str, design_date: str, batch: str) -> str:
    """`<outdir>/<date>/<batch>`, or `<outdir>/<date>` when no batch label is used."""
    parts = [outdir, design_date]
    if batch:
        parts.append(batch)
    return os.path.join(*parts)


def allocate_batch(outdir: str, design_date: str, requested: Optional[str] = None) -> str:
    """Pick a batch label that does not already hold results.

    Several batches can be requested on the same day, so results are never keyed on the date
    alone. Without an explicit label this returns the next free `b1`, `b2`, ... for that date;
    with one, it checks that the label is unused and raises rather than overwrite anything.
    """
    day = os.path.join(outdir, design_date)
    if requested:
        label = sanitise_label(requested)
        if os.path.isdir(os.path.join(day, label)):
            raise FileExistsError(
                f"{os.path.join(day, label)} already holds results. Choose a different "
                f"--batch label; existing designs are never overwritten."
            )
        return label
    existing = set(os.listdir(day)) if os.path.isdir(day) else set()
    n = 1
    while f"b{n}" in existing:
        n += 1
    return f"b{n}"


def sanitise_label(label: str) -> str:
    safe = "".join(c if c.isalnum() or c in "-_." else "_" for c in label.strip())
    if not safe:
        raise ValueError(f"batch label {label!r} contains no usable characters")
    return safe


def design_dir(outdir: str, design_date: str, result: DesignResult) -> str:
    path = os.path.join(batch_dir(outdir, design_date, result.batch), design_folder_name(result))
    os.makedirs(path, exist_ok=True)
    return path


def design_folder_name(result: DesignResult) -> str:
    """`<target>_from_<sources>`, joined with + and stripped of anything path-unsafe."""
    sources = result.source_names or [result.template_name]
    joined = "+".join(sources)
    safe = "".join(c if c.isalnum() or c in "-_.+" else "_" for c in joined)
    return f"{result.target_name}_from_{safe}"


def primer_rows(result: DesignResult) -> List[Dict[str, object]]:
    """One row per orderable oligo, tagged with the USER junction it builds."""
    amplified_by = {}
    for frag in result.fragments:
        amplified_by[id(frag.forward)] = frag.source_name
        amplified_by[id(frag.reverse)] = frag.source_name

    rows = []
    for junction in result.junctions:
        for primer in (junction.forward, junction.reverse):
            row: Dict[str, object] = {
                "design_id": result.design_id,
                "design_date": result.design_date,
                "batch": result.batch,
                "target": result.target_name,
                "template": result.template_name,
            }
            row.update(primer.stats_row())
            row["primer_name"] = f"{result.target_name}_{primer.name}_{result.stamp}"
            row["user_junction_sequence"] = junction.top_overhang
            row["user_junction_length_nt"] = junction.block_len
            row["user_junction_duplex_tm_C"] = round(junction.overhang_tm, 1)
            row["three_prime_overhang_this_end"] = primer.overhang
            row["amplified_from"] = amplified_by.get(id(primer), result.backbone_name)
            rows.append(row)
    return rows


def write_primer_csv(path: str, rows: Sequence[Dict[str, object]]) -> None:
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=PRIMER_COLUMNS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_order_form(path: str, rows: Sequence[Dict[str, object]]) -> None:
    """Tab-separated name/sequence/scale/purification, ready to paste into an IDT bulk order."""
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["Name", "Sequence", "Scale", "Purification"])
        for row in rows:
            long_oligo = int(row["length_nt"]) > 60
            scale = "4nm Ultramer" if long_oligo else "25nm"
            writer.writerow([row["primer_name"], row["order_sequence"], scale, "STD"])


def write_batch_tables(outdir: str, design_date: str, batch: str,
                       results: Sequence[DesignResult]) -> Dict[str, str]:
    """Combined primer table, order form and summary covering every design in the batch."""
    folder = batch_dir(outdir, design_date, batch)
    os.makedirs(folder, exist_ok=True)
    stamp = results[0].stamp if results else design_date.replace("-", "")
    rows = [row for result in results for row in primer_rows(result)]
    paths: Dict[str, str] = {}
    paths["batch_summary"] = os.path.join(folder, f"summary_{stamp}.csv")
    write_summary_csv(paths["batch_summary"], results)
    if rows:
        paths["batch_primers"] = os.path.join(folder, f"primers_{stamp}.csv")
        write_primer_csv(paths["batch_primers"], rows)
        paths["batch_order_form"] = os.path.join(folder, f"order_{stamp}.tsv")
        write_order_form(paths["batch_order_form"], rows)
    return paths


def write_summary_csv(path: str, results: Sequence[DesignResult]) -> None:
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=SUMMARY_COLUMNS, extrasaction="ignore")
        writer.writeheader()
        for result in results:
            writer.writerow(result.summary_row())


def write_json(path: str, result: DesignResult) -> None:
    comparison = result.comparison
    payload: Dict[str, object] = {
        "design_id": result.design_id,
        "design_date": result.design_date,
        "sources": result.source_names,
        "backbone": {"name": result.backbone_name, "length": result.template_length},
        "target": {"name": result.target_name, "length": result.target_length},
        "status": result.status,
        "verified": result.verified,
        "messages": result.messages,
        "warnings": result.warnings,
        "parameters": _params_dict(result),
    }
    if comparison is not None:
        payload["comparison"] = {
            "identity": round(comparison.identity, 6),
            "template_rotation": comparison.template_rotation,
            "target_rotation": comparison.target_rotation,
            "target_reverse_complemented": comparison.target_reverse_complemented,
            "edits": [
                {
                    "kind": e.kind,
                    "template_start_1based": e.tpl_start + 1,
                    "template_end": e.tpl_end,
                    "target_start_1based": e.tgt_start + 1,
                    "target_end": e.tgt_end,
                    "template_sequence": e.tpl_seq,
                    "target_sequence": e.tgt_seq,
                }
                for e in comparison.edits
            ],
        }
    payload["template_repeats"] = [
        {
            "start_1based": b.start + 1,
            "end_1based": b.end + 1,
            "length": b.length,
            "copies": b.copies,
        }
        for b in result.template_repeats
    ]
    payload["junctions"] = [
        {
            "junction": j.index + 1,
            "overhang_top_strand": j.top_overhang,
            "overhang_bottom_strand": j.bottom_overhang,
            "overhang_length": j.block_len,
            "overhang_duplex_tm_C": round(j.overhang_tm, 1),
            "target_position_1based": (j.block_start % max(1, result.target_length)) + 1,
            "penalty": round(j.penalty, 3),
            "penalty_detail": {k: round(v, 3) for k, v in j.penalty_detail.items()},
            "warnings": j.warnings,
        }
        for j in result.junctions
    ]
    payload["fragments"] = [
        {
            "fragment": f.name,
            "amplified_from": f.source_name,
            "forward_primer": f.forward.name,
            "reverse_primer": f.reverse.name,
            "left_junction": f.left_junction + 1,
            "right_junction": f.right_junction + 1,
            "template_start_1based": f.template_start + 1,
            "template_end": f.template_end,
            "expected_product_length": f.expected_length,
        }
        for f in result.fragments
    ]
    payload["primers"] = primer_rows(result)
    payload["pcr"] = [asdict(p) for p in result.protocols]
    if result.verification is not None:
        payload["verification"] = {
            "ok": result.verification.ok,
            "assembled_length": result.verification.assembled_length,
            "checks": [
                {"check": name, "passed": passed, "detail": detail}
                for name, passed, detail in result.verification.checks
            ],
            "errors": result.verification.errors,
        }
    with open(path, "w") as fh:
        json.dump(payload, fh, indent=2)
        fh.write("\n")


def _params_dict(result: DesignResult) -> Dict[str, object]:
    params = asdict(result.params)
    params["conditions"] = asdict(result.params.conditions)
    return params


def write_report(path: str, result: DesignResult, input_file: str) -> None:
    with open(path, "w") as fh:
        fh.write("\n".join(_report_lines(result, input_file)) + "\n")


def _report_lines(result: DesignResult, input_file: str) -> List[str]:
    L: List[str] = []
    L.append(f"# USER cloning design: {result.target_name} from {result.template_name}")
    L.append("")
    L.append(f"- **Design ID:** `{result.design_id}`")
    L.append(f"- **Design date:** {result.design_date}"
             + (f"  •  **batch {result.batch}**" if result.batch else ""))
    L.append(f"- **Input file:** `{os.path.basename(input_file)}`")
    if len(result.source_names) > 1:
        L.append(f"- **Source plasmids:** {', '.join(result.source_names)}")
        L.append(f"- **Backbone:** {result.backbone_name} "
                 f"({result.template_length:,} bp, circular)")
    else:
        L.append(f"- **Template:** {result.template_name} "
                 f"({result.template_length:,} bp, circular)")
    L.append(f"- **Target:** {result.target_name} ({result.target_length:,} bp, circular)")
    L.append(f"- **Status:** {result.status}"
             + ("  •  **in-silico assembly verified**" if result.verified else ""))
    L.append("")

    if result.status == STATUS_NO_DIFFERENCE:
        L.append("## No design produced")
        L.append("")
        for msg in result.messages:
            L.append(f"> {msg}")
        L.append("")
        return L

    if result.status != STATUS_OK:
        L.append("## Design incomplete")
        L.append("")
        for msg in result.messages:
            L.append(f"- {msg}")
        L.append("")

    L.extend(_differences_section(result))
    L.extend(_repeat_section(result))
    L.extend(_primer_section(result))
    L.extend(_junction_section(result))
    L.extend(_protocol_section(result))
    L.extend(_verification_section(result))

    if result.warnings:
        L.append("## Warnings")
        L.append("")
        for w in dict.fromkeys(result.warnings):
            L.append(f"- {w}")
        L.append("")

    L.append("## How to read the primer sequences")
    L.append("")
    L.append(f"`{U_MARK}` marks the deoxyuridine (order as dU / deoxyuridine; IDT writes it "
             f"exactly this way). Each primer is:")
    L.append("")
    L.append("```")
    L.append("5'-[ tail ][U][ extra ][ annealing region ]-3'")
    L.append("    ^^^^^^^^^^  removed by USER enzyme      ^ binds the template")
    L.append("```")
    L.append("")
    L.append("The tail plus the U position are excised, leaving a 3' overhang on the opposite "
             "strand. Paired fragment ends carry complementary overhangs, so they anneal "
             "directionally without ligase.")
    L.append("")
    return L


def _differences_section(result: DesignResult) -> List[str]:
    comparison = result.comparison
    if comparison is None:
        return []
    n = len(comparison.edits)
    L = ["## What changes between the two plasmids", ""]
    found = f"{n} difference site" if n == 1 else f"{n} difference sites"
    handled = "handled by one USER junction" if n == 1 else "each handled by its own USER junction"
    L.append(f"Sequence identity in the aligned frame: **{100 * comparison.identity:.2f}%**. "
             f"{found} found, {handled}.")
    L.append("")
    L.append("| Site | Type | Position in template (input file) | Position in target (input file) "
             "| Template sequence | Target sequence |")
    L.append("| --- | --- | --- | --- | --- | --- |")
    for i, e in enumerate(comparison.edits, start=1):
        L.append(
            f"| {i} | {e.kind} "
            f"| {_input_span(comparison.to_template_input(e.tpl_start), e.removed_length)} "
            f"| {_input_span(comparison.to_target_input(e.tgt_start), e.new_length)} "
            f"| `{_short(e.tpl_seq)}` | `{_short(e.tgt_seq)}` |"
        )
    L.append("")
    L.append("Positions are 1-based on the sequences exactly as they appear in the input file.")
    if comparison.target_reverse_complemented:
        L.append("")
        L.append("**Note:** the target sequence in the input file is on the opposite strand to "
                 "the template, so it was reverse-complemented before designing. Target "
                 "positions above refer to that reverse complement.")
    L.append("")
    return L


def _input_span(start: int, length: int) -> str:
    """A 1-based coordinate span; insertions have zero length and sit between two bases."""
    if length == 0:
        return f"between {start} and {start + 1} (insertion point)"
    return f"{start + 1}–{start + length}"


def _repeat_section(result: DesignResult) -> List[str]:
    if not result.template_repeats:
        return []
    length = result.template_length
    L = ["## Repeated regions of the template", ""]
    total = sum(b.length for b in result.template_repeats)
    L.append(
        f"{total:,} nt ({100.0 * total / max(1, length):.1f}% of {result.backbone_name}) lies in "
        "exactly repeated blocks of 40 nt or more. A primer inside one of these binds more than "
        "once, so USER junctions cannot normally be placed there. Positions are 1-based on the "
        "template as supplied."
    )
    L.append("")
    L.append("| Repeated block | Length | Copies |")
    L.append("| --- | --- | --- |")
    for b in result.template_repeats:
        L.append(f"| {b.start + 1:,}–{b.end + 1:,} | {b.length:,} nt | {b.copies} |")
    L.append("")
    stretches = [(s, e) for s, e in unique_stretches("N" * length, result.template_repeats)
                 if e - s >= 100]
    if stretches:
        readable = ", ".join(f"{s + 1:,}–{e:,}" for s, e in stretches)
        L.append(f"Junctions can be placed freely in the uniquely primable stretches: {readable}.")
        L.append("")
    return L


def _primer_section(result: DesignResult) -> List[str]:
    rows = primer_rows(result)
    L = ["## Primers to order", ""]
    L.append("| Name | Sequence (5'→3') | nt | U pos | Anneal Tm | GC% | USER junction |")
    L.append("| --- | --- | --- | --- | --- | --- | --- |")
    for row in rows:
        L.append(
            f"| {row['primer_name']} | `{row['order_sequence']}` | {row['length_nt']} "
            f"| {row['u_position_from_5prime']} | {row['annealing_tm_C']} °C "
            f"| {row['gc_percent']} | `{row['user_junction_sequence']}` |"
        )
    L.append("")
    L.append(f"_{USER_JUNCTION_COLUMN_NOTE}_")
    L.append("")
    L.append("### Primer anatomy")
    L.append("")
    for primer, row in zip(result.primers, rows):
        L.append(f"**{row['primer_name']}** ({primer.direction}, junction {primer.junction + 1})")
        L.append("")
        L.append("```")
        L.append(f"5'-{primer.tail[:-1]}[{primer.tail[-1]}=dU]{primer.extra}{primer.anneal}-3'")
        L.append(f"   tail            = {primer.tail}   ({len(primer.tail)} nt, dU at position "
                 f"{primer.u_index + 1})")
        L.append(f"   non-templated   = {primer.extra or '(none)'}")
        L.append(f"   annealing       = {primer.anneal}   ({len(primer.anneal)} nt, "
                 f"Tm {primer.anneal_tm:.1f} C, {primer.template_hits} site in template)")
        L.append(f"   3' overhang     = {primer.overhang}   (exposed after USER)")
        L.append("```")
        L.append("")
    return L


def _junction_section(result: DesignResult) -> List[str]:
    L = ["## Junctions and fragments", ""]
    L.append("| Junction | Overhang (top strand 5'→3') | nt | Duplex Tm | Joins |")
    L.append("| --- | --- | --- | --- | --- |")
    for j in result.junctions:
        joins = j.label or (j.edit.kind if j.edit else "-")
        if j.edit is not None and j.crosses_sources:
            joins = f"{j.label} ({j.edit.kind} carried on the primers)"
        L.append(f"| J{j.index + 1} | `{j.top_overhang}` | {j.block_len} "
                 f"| {j.overhang_tm:.1f} °C | {joins} |")
    L.append("")
    L.append("How each junction anneals after USER treatment:")
    L.append("")
    L.append("```")
    L.append("\n\n".join(junction_diagram(j) for j in result.junctions))
    L.append("```")
    L.append("")
    L.append("| Fragment | Amplified from | Forward | Reverse | Template span | Product length |")
    L.append("| --- | --- | --- | --- | --- | --- |")
    comparison = result.comparison
    for f in result.fragments:
        on_backbone = comparison is not None and f.source_name == result.backbone_name
        start = comparison.to_template_input(f.template_start) if on_backbone else f.template_start
        end = (comparison.to_template_input(f.template_end - 1) if on_backbone
               else f.template_end - 1)
        span = f"{start + 1:,}..{end + 1:,}"
        if len(result.fragments) == 1:
            span += " (whole plasmid, wraps the origin)"
        L.append(
            f"| {f.name} | {f.source_name} | {result.target_name}_{f.forward.name} "
            f"| {result.target_name}_{f.reverse.name} | {span} | {f.expected_length:,} bp |"
        )
    L.append("")
    if len(result.fragments) == 1:
        L.append("A single product spanning the whole plasmid is amplified; USER treatment "
                 "leaves complementary overhangs at its two ends so it self-circularises.")
        L.append("")
    else:
        sources = ", ".join(dict.fromkeys(f.source_name for f in result.fragments))
        L.append(f"{len(result.fragments)} products are amplified (from {sources}) and joined at "
                 f"{len(result.junctions)} junctions with distinct overhangs, so assembly is "
                 "directional: each pair of ends can only anneal to its intended partner.")
        L.append("")
    return L


def _protocol_section(result: DesignResult) -> List[str]:
    if not result.protocols:
        return []
    cond = result.params.conditions
    L = ["## Suggested wet-lab protocol", ""]
    L.append("### 1. PCR (uracil-tolerant polymerase)")
    L.append("")
    L.append("Use **Phusion U Hot Start**, **Q5U**, or **PfuTurbo Cx**. A standard proofreading "
             "polymerase will stall at the dU and must not be used.")
    L.append("")
    L.append("| Fragment | Template | Annealing temp | Extension | Product |")
    L.append("| --- | --- | --- | --- | --- |")
    source_of = {f.name: f.source_name for f in result.fragments}
    for p in result.protocols:
        L.append(f"| {p.fragment} | {source_of.get(p.fragment, '-')} | {p.annealing_temp_C} °C "
                 f"| {p.extension_seconds} s | {p.product_length:,} bp |")
    L.append("")
    L.append(f"Tm values assume {cond.monovalent_mM:.0f} mM monovalent salt, "
             f"{cond.divalent_mM:.1f} mM Mg²⁺, {cond.dntp_mM:.1f} mM dNTP, "
             f"{cond.primer_nM:.0f} nM primer. Annealing temperature is derived from the "
             "template-binding region only, because the 5' tails do not pair in the first cycles.")
    L.append("")
    L.append("Cycling: 98 °C 30 s; 30 cycles of [98 °C 10 s, Ta 20 s, 72 °C extension]; "
             "72 °C 5 min.")
    L.append("")
    L.append("### 2. USER treatment and annealing")
    L.append("")
    L.append("- Column- or gel-purify the PCR product, then combine ~0.1–0.2 pmol with 1 U USER "
             "enzyme in 10 µL 1× buffer (T4 ligase buffer or CutSmart both work).")
    L.append("- 37 °C for 25 min to excise the uracils.")
    weakest = min((j.overhang_tm for j in result.junctions), default=99.0)
    if weakest < 25.0:
        L.append(f"- Then anneal by ramping 37 → 10 °C at 0.1 °C/s and holding at 10 °C for "
                 f"20 min. The least stable overhang here has a duplex Tm of "
                 f"{weakest:.0f} °C, so a room-temperature hold would leave it largely "
                 "single-stranded.")
    else:
        L.append("- Then anneal at 25 °C for 25 min (a 37 → 10 °C ramp works too).")
    L.append("- Transform 1–2 µL directly into chemically competent cells. **No ligase is "
             "needed** — the nicks are sealed in vivo.")
    L.append("")
    L.append("### 3. Screening")
    L.append("")
    L.append("Sanger sequencing across each junction confirms the design; the junction sequences "
             "to look for are in the table above.")
    L.append("")
    L.append("**No restriction enzyme is used anywhere in this workflow.** The only enzymes are "
             "the uracil-tolerant polymerase and USER; the overhangs come from uracil excision, "
             "not from any restriction site, so no site needs to be present, added or preserved.")
    L.append("")
    return L


def _verification_section(result: DesignResult) -> List[str]:
    v = result.verification
    if v is None:
        return []
    L = ["## In-silico verification", ""]
    for name, passed, detail in v.checks:
        mark = "✅" if passed else "❌"
        L.append(f"- {mark} {name}" + (f" — {detail}" if detail else ""))
    L.append("")
    if v.ok:
        L.append(f"The simulated assembly is **{v.assembled_length:,} bp** and matches "
                 f"{result.target_name} exactly as a circular sequence.")
    else:
        L.append("**The design did not verify — do not order these primers.**")
    L.append("")
    return L


def write_all(
    outdir: str,
    result: DesignResult,
    input_file: str,
) -> Dict[str, str]:
    """Write every output file for one design and return {kind: path}."""
    folder = design_dir(outdir, result.design_date, result)
    base = f"{design_folder_name(result)}_{result.stamp}"
    paths: Dict[str, str] = {}

    paths["report"] = os.path.join(folder, f"report_{base}.md")
    write_report(paths["report"], result, input_file)

    paths["json"] = os.path.join(folder, f"design_{base}.json")
    write_json(paths["json"], result)

    rows = primer_rows(result)
    if rows:
        paths["primers"] = os.path.join(folder, f"primers_{base}.csv")
        write_primer_csv(paths["primers"], rows)
        paths["order_form"] = os.path.join(folder, f"order_{base}.tsv")
        write_order_form(paths["order_form"], rows)

    predicted = predicted_plasmid(result)
    if predicted:
        paths["predicted_fasta"] = os.path.join(folder, f"predicted_{base}.fasta")
        write_fasta(
            paths["predicted_fasta"],
            f"{result.target_name}_predicted",
            predicted,
            f"assembled in silico from {result.template_name} | {result.design_id} | "
            f"{result.design_date} | {len(predicted)} bp circular",
        )
    return paths


def junction_diagram(junction: Junction) -> str:
    """Text diagram of one annealed junction: the two 3' overhangs base-paired."""
    top = junction.top_overhang
    label = f"J{junction.index + 1}"
    pad = " " * len(label)
    return "\n".join([
        f"{label}  upstream fragment   3'-overhang  5'-{top}-3'",
        f"{pad}                                     {'|' * len(top)}",
        f"{pad}  downstream fragment 3'-overhang  3'-{revcomp(top)[::-1]}-5'",
        f"{pad}  duplex Tm {junction.overhang_tm:.1f} C",
    ])


def _short(seq: str, limit: int = 36) -> str:
    if not seq:
        return "—"
    if len(seq) <= limit:
        return seq
    return f"{seq[:16]}…{seq[-16:]} ({len(seq)} nt)"
