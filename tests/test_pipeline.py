"""Tests for the USER cloning pipeline.

Every design case is checked two ways: the pipeline's own in-silico verification must
pass, and the predicted plasmid must independently match the requested target as a
circular sequence. The real 8,116 bp pLL057P backbone is used as the template so the
tests exercise realistic sequence composition, repeats included.
"""

from __future__ import annotations

import os
import subprocess
import sys
import tempfile
import unittest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)

from user_cloning import DesignParams, SeqRecord, compare_plasmids, read_sequence_table  # noqa: E402
from user_cloning.design import U_MARK  # noqa: E402
from user_cloning.report import allocate_batch, primer_rows  # noqa: E402
from user_cloning.pipeline import (  # noqa: E402
    STATUS_NO_DIFFERENCE,
    STATUS_OK,
    design_assembly,
    design_one,
    predicted_plasmid,
)
from user_cloning.sequences import (  # noqa: E402
    circular_equal,
    count_circular_occurrences,
    repeat_blocks,
    revcomp,
    rotate,
    unique_stretches,
)

def _first_existing(*names: str) -> str:
    for name in names:
        path = os.path.join(REPO, name)
        if os.path.exists(path):
            return path
    return os.path.join(REPO, names[0])


INPUT_XLSX = _first_existing("USER_design_1.xlsx", "USER_design.xlsx")
DATE = "2026-07-30"


def load_backbone() -> str:
    records = read_sequence_table(INPUT_XLSX)
    return records[0].seq


BACKBONE = load_backbone()


def rec(name: str, seq: str) -> SeqRecord:
    return SeqRecord(name=name, seq=seq, source_row=1)


class DesignCaseMixin:
    """Shared assertions for a template -> target design."""

    def run_design(self, target_seq: str, name: str = "pTGT", params: DesignParams = None):
        template = rec("pLL057P", BACKBONE)
        target = rec(name, target_seq)
        result = design_one(template, target, DATE, params or DesignParams())
        return result

    def assert_good_design(self, target_seq: str, expected_fragments: int = 1, name: str = "pTGT"):
        result = self.run_design(target_seq, name=name)
        self.assertEqual(result.status, STATUS_OK, msg="; ".join(result.messages))
        self.assertTrue(result.verified, msg="; ".join(result.verification.errors))
        self.assertEqual(len(result.fragments), expected_fragments)
        self.assertEqual(len(result.primers), 2 * expected_fragments)

        predicted = predicted_plasmid(result)
        equal, _ = circular_equal(predicted, target_seq)
        self.assertTrue(equal, "predicted plasmid does not match the requested target")
        self.assertEqual(len(predicted), len(target_seq))

        for primer in result.primers:
            self.assertEqual(primer.sequence[primer.u_index], "T",
                             "the deoxyuridine must replace a T")
            self.assertEqual(primer.tail[-1], "T")
            self.assertEqual(primer.tail[0] if primer.direction == "forward" else
                             revcomp(primer.tail)[0], "A")
            self.assertGreaterEqual(len(primer.anneal), DesignParams().anneal_min)
            # The annealing region must occur exactly once in the circular template. Forward
            # primers match the top strand, reverse primers the bottom one.
            self.assertEqual(
                count_circular_occurrences(BACKBONE, primer.anneal, both_strands=True), 1,
                f"{primer.name} annealing region is not a unique site in the template",
            )
            self.assertEqual(primer.template_hits, 1)
            self.assertEqual(primer.order_sequence.replace(U_MARK, "T"), primer.sequence)

        for junction in result.junctions:
            self.assertEqual(junction.top_overhang, revcomp(junction.bottom_overhang))
            self.assertTrue(junction.top_overhang.startswith("A"))
            self.assertTrue(junction.top_overhang.endswith("T"))
        return result


class TestPointMutation(DesignCaseMixin, unittest.TestCase):
    def test_single_base_change(self):
        pos = 4000
        original = BACKBONE[pos]
        new = "A" if original != "A" else "C"
        target = BACKBONE[:pos] + new + BACKBONE[pos + 1:]
        result = self.assert_good_design(target)
        self.assertEqual(len(result.comparison.edits), 1)
        self.assertEqual(len(target), len(BACKBONE))

    def test_three_codon_substitution(self):
        pos = 2500
        block = "GCTGCAGCTGCAG"
        target = BACKBONE[:pos] + block + BACKBONE[pos + len(block):]
        self.assert_good_design(target)


class TestInsertions(DesignCaseMixin, unittest.TestCase):
    def test_short_tag_insertion(self):
        pos = 3300
        insert = "GGTGGTAGCGGTGGTAGC"  # GS linker
        target = BACKBONE[:pos] + insert + BACKBONE[pos:]
        result = self.assert_good_design(target)
        self.assertEqual(len(predicted_plasmid(result)), len(BACKBONE) + len(insert))
        self.assertEqual(result.comparison.edits[0].kind, "insertion")

    def test_long_insertion_is_split_across_both_primers(self):
        pos = 5000
        insert = "GATTACA" * 9  # 63 nt: too long for one overhang block
        target = BACKBONE[:pos] + insert + BACKBONE[pos:]
        result = self.assert_good_design(target)
        junction = result.junctions[0]
        carried = junction.block_len + len(junction.forward.extra) + len(junction.reverse.extra)
        self.assertEqual(carried, len(insert),
                         "the whole insertion must be encoded on the primer pair exactly once")
        self.assertTrue(junction.forward.extra and junction.reverse.extra,
                        "the insertion should be shared between the two primers")
        lengths = sorted(p.length for p in result.primers)
        self.assertLess(lengths[-1] - lengths[0], 30,
                        f"primer lengths should be roughly balanced, got {lengths}")

    def test_very_long_insertion_advises_a_synthesised_fragment(self):
        pos = 5000
        insert = "GGTGGTAGCGGTGGTAGCGGT" * 6 + "GGTGGT"  # 132 nt
        target = BACKBONE[:pos] + insert + BACKBONE[pos:]
        result = self.assert_good_design(target)
        self.assertTrue(any("synthesised dsDNA fragment" in m for m in result.messages))
        self.assertTrue(any(p.length > 60 for p in result.primers))


class TestDeletions(DesignCaseMixin, unittest.TestCase):
    def test_small_deletion(self):
        target = BACKBONE[:1500] + BACKBONE[1533:]
        result = self.assert_good_design(target)
        self.assertEqual(result.comparison.edits[0].kind, "deletion")

    def test_large_deletion(self):
        target = BACKBONE[:5000] + BACKBONE[5900:]
        result = self.assert_good_design(target)
        self.assertEqual(len(predicted_plasmid(result)), len(BACKBONE) - 900)

    def test_replacement_of_unequal_length(self):
        target = BACKBONE[:6000] + "ACGTACGTACGTACGT" + BACKBONE[6200:]
        self.assert_good_design(target)


class TestMultipleEdits(DesignCaseMixin, unittest.TestCase):
    def test_two_distant_edits_give_two_fragments(self):
        target = (
            BACKBONE[:1000] + "CCCGGGCCCGGG" + BACKBONE[1012:5000]
            + "TTAATTAATTAA" + BACKBONE[5012:]
        )
        result = self.assert_good_design(target, expected_fragments=2)
        self.assertEqual(len(result.junctions), 2)
        overhangs = [j.top_overhang for j in result.junctions]
        self.assertEqual(len(set(overhangs)), 2, "junction overhangs must differ")

    def test_three_edits_give_three_fragments(self):
        target = (
            BACKBONE[:800] + "AGGCCTAGGCCT" + BACKBONE[812:3000]
            + "CATCATCATCAT" + BACKBONE[3012:6000]
            + "TGGTGGTGGTGG" + BACKBONE[6012:]
        )
        result = self.assert_good_design(target, expected_fragments=3)
        self.assertEqual(len(set(j.top_overhang for j in result.junctions)), 3)


class TestCircularAndOrientationHandling(DesignCaseMixin, unittest.TestCase):
    def test_edit_at_sequence_origin(self):
        target = "GCGCGCGC" + BACKBONE[8:]
        self.assert_good_design(target)

    def test_edit_wrapping_the_origin(self):
        target = "TTTTT" + BACKBONE[5:-5] + "AAAAA"
        result = self.assert_good_design(target)
        self.assertTrue(result.verified)

    def test_target_supplied_as_rotation(self):
        pos = 4000
        mutant = BACKBONE[:pos] + "GGGGCCCC" + BACKBONE[pos + 8:]
        rotated = rotate(mutant, 3712)
        result = self.assert_good_design(rotated)
        self.assertEqual(len(result.comparison.edits), 1)

    def test_target_supplied_reverse_complemented(self):
        pos = 3000
        mutant = BACKBONE[:pos] + "TTAACCGGTTAA" + BACKBONE[pos + 12:]
        flipped = revcomp(rotate(mutant, 1234))
        result = self.assert_good_design(flipped)
        self.assertTrue(result.comparison.target_reverse_complemented)
        self.assertEqual(len(result.comparison.edits), 1)


class TestMultiTemplateAssembly(unittest.TestCase):
    """Building a target from two source plasmids: backbone from one, insert from the other."""

    def assemble(self, sources, target_seq, name="pASM"):
        target = rec(name, target_seq)
        result = design_assembly(sources, target, DATE)
        self.assertEqual(result.status, STATUS_OK, "; ".join(result.messages))
        self.assertTrue(result.verified, "; ".join(result.verification.errors))
        predicted = predicted_plasmid(result)
        equal, _ = circular_equal(predicted, target_seq)
        self.assertTrue(equal, "predicted plasmid does not match the requested target")
        for frag in result.fragments:
            # Each fragment must genuinely amplify off the plasmid it is assigned to.
            self.assertEqual(
                count_circular_occurrences(frag.source_seq, frag.forward.anneal), 1)
            self.assertEqual(
                count_circular_occurrences(frag.source_seq, frag.reverse.anneal), 1)
        return result

    def test_cassette_moved_between_plasmids(self):
        donor_payload = "".join(BACKBONE[i] for i in range(600, 1200))  # 600 nt from elsewhere
        donor = rec("pDONOR", BACKBONE[3000:4000] + donor_payload + BACKBONE[4000:5000])
        backbone = rec("pBB", BACKBONE)
        target_seq = BACKBONE[:5300] + donor_payload + BACKBONE[5300:]
        result = self.assemble([backbone, donor], target_seq)
        self.assertEqual(len(result.fragments), 2)
        self.assertEqual(len(result.junctions), 2)
        self.assertEqual({f.source_name for f in result.fragments}, {"pBB", "pDONOR"})
        self.assertEqual(
            sum(1 for f in result.fragments if f.source_name == "pDONOR"), 1)

    def test_backbone_is_the_better_matching_source(self):
        payload = BACKBONE[200:900]
        donor = rec("pSMALL", BACKBONE[3000:3600] + payload)
        backbone = rec("pBIG", BACKBONE)
        target_seq = BACKBONE[:5300] + payload + BACKBONE[5300:]
        result = self.assemble([donor, backbone], target_seq)
        self.assertEqual(result.backbone_name, "pBIG",
                         "the source explaining most of the target must be the backbone")

    def test_source_order_does_not_change_the_design(self):
        payload = BACKBONE[300:1000]
        donor = rec("pD", BACKBONE[2600:3200] + payload)
        backbone = rec("pB", BACKBONE)
        target_seq = BACKBONE[:5300] + payload + BACKBONE[5300:]
        first = self.assemble([backbone, donor], target_seq)
        second = self.assemble([donor, backbone], target_seq)
        self.assertEqual([p.sequence for p in first.primers],
                         [p.sequence for p in second.primers])

    def test_short_insert_stays_on_the_primers(self):
        """A 20 nt insert present in the donor is cheaper to write into the tails."""
        payload = "GGTGGTAGCGGTGGTAGCGG"
        donor = rec("pD2", BACKBONE[1000:1600] + payload + BACKBONE[1600:2000])
        backbone = rec("pB2", BACKBONE)
        target_seq = BACKBONE[:5300] + payload + BACKBONE[5300:]
        result = self.assemble([backbone, donor], target_seq)
        self.assertEqual(len(result.fragments), 1, "should stay a single-fragment design")
        self.assertEqual({f.source_name for f in result.fragments}, {"pB2"})

    def test_single_source_still_works_through_design_assembly(self):
        target_seq = BACKBONE[:3400] + "GGTTAACCGG" + BACKBONE[3410:]
        result = self.assemble([rec("pLL057P", BACKBONE)], target_seq)
        self.assertEqual(len(result.fragments), 1)



class TestBatching(unittest.TestCase):
    """Several batches can land on the same day, and results are never destroyed."""

    def test_auto_labels_increment_per_date(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.assertEqual(allocate_batch(tmp, DATE), "b1")
            os.makedirs(os.path.join(tmp, DATE, "b1"))
            self.assertEqual(allocate_batch(tmp, DATE), "b2")
            os.makedirs(os.path.join(tmp, DATE, "b2"))
            self.assertEqual(allocate_batch(tmp, DATE), "b3")
            # A different date starts over.
            self.assertEqual(allocate_batch(tmp, "2026-08-01"), "b1")

    def test_explicit_label_is_refused_when_already_used(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.assertEqual(allocate_batch(tmp, DATE, "cassette"), "cassette")
            os.makedirs(os.path.join(tmp, DATE, "cassette"))
            with self.assertRaises(FileExistsError):
                allocate_batch(tmp, DATE, "cassette")

    def test_labels_are_made_path_safe(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.assertEqual(allocate_batch(tmp, DATE, "batch 2/final"), "batch_2_final")
            with self.assertRaises(ValueError):
                allocate_batch(tmp, DATE, "  ")

    def test_two_runs_on_one_day_both_survive(self):
        with tempfile.TemporaryDirectory() as tmp:
            src = os.path.join(tmp, "in.csv")
            write_pair_csv(src, [
                ("pLL057P", BACKBONE),
                ("pA", BACKBONE[:3400] + "GGTTAACCGG" + BACKBONE[3410:]),
                ("pB", BACKBONE[:5200] + BACKBONE[5260:]),
            ])
            out = os.path.join(tmp, "designs")
            for target in ("pA", "pB"):
                proc = subprocess.run(
                    [sys.executable, "design_user_primers.py", src,
                     "--template", "pLL057P", "--target", target,
                     "--outdir", out, "--date", DATE],
                    cwd=REPO, capture_output=True, text=True,
                )
                self.assertEqual(proc.returncode, 0, proc.stdout + proc.stderr)
            self.assertTrue(os.path.isdir(os.path.join(out, DATE, "b1", "pA_from_pLL057P")))
            self.assertTrue(os.path.isdir(os.path.join(out, DATE, "b2", "pB_from_pLL057P")))
            # The first batch's files must still be intact.
            stamp1 = f"{DATE.replace('-', '')}_b1"
            self.assertTrue(os.path.exists(os.path.join(
                out, DATE, "b1", "pA_from_pLL057P", f"primers_pA_from_pLL057P_{stamp1}.csv")))


class TestUserJunctionColumn(unittest.TestCase):
    def test_junction_sequence_is_shared_by_the_pair_and_matches_the_overhang(self):
        target_seq = BACKBONE[:1000] + "CCCGGGCCCGGG" + BACKBONE[1012:5000] + \
            "TTAATTAATTAA" + BACKBONE[5012:]
        result = design_one(rec("pLL057P", BACKBONE), rec("pJ", target_seq), DATE, batch="b7")
        self.assertTrue(result.verified)
        rows = primer_rows(result)
        self.assertEqual(len(rows), 4)
        by_junction = {}
        for row in rows:
            by_junction.setdefault(row["junction"], set()).add(row["user_junction_sequence"])
        for junction, seqs in by_junction.items():
            self.assertEqual(len(seqs), 1,
                             f"both primers of junction {junction} must report one sequence")
        designed = {j.top_overhang for j in result.junctions}
        self.assertEqual({s for seqs in by_junction.values() for s in seqs}, designed)
        for row in rows:
            # The overhang this end presents is the junction sequence or its complement.
            self.assertIn(row["three_prime_overhang_this_end"],
                          {row["user_junction_sequence"],
                           revcomp(str(row["user_junction_sequence"]))})
            self.assertEqual(row["batch"], "b7")
            self.assertTrue(str(row["primer_name"]).endswith("_b7"))


class TestRepeatedRegions(unittest.TestCase):
    """pLL057P carries two CMV enhancer/promoter copies and two LTR blocks. Primers inside
    those cannot be specific, and the pipeline must say so rather than emit bad oligos."""

    def test_backbone_repeat_map(self):
        blocks = repeat_blocks(BACKBONE)
        self.assertTrue(blocks, "pLL057P is known to contain repeated elements")
        self.assertTrue(any(b.length > 400 for b in blocks),
                        "expected the duplicated ~500 nt promoter block")
        for block in blocks:
            self.assertGreaterEqual(block.copies, 2)
            motif = BACKBONE[block.start:block.start + 40]
            self.assertGreater(count_circular_occurrences(BACKBONE, motif, both_strands=False), 1)

    def test_unique_stretches_partition_the_sequence(self):
        blocks = repeat_blocks(BACKBONE)
        stretches = unique_stretches(BACKBONE, blocks)
        for start, end in stretches:
            for block in blocks:
                self.assertFalse(block.overlaps(start, end),
                                 "unique stretches must not overlap repeat blocks")

    def test_edit_inside_a_repeat_is_refused_with_a_reason(self):
        blocks = [b for b in repeat_blocks(BACKBONE) if b.length > 400]
        middle = (blocks[0].start + blocks[0].end) // 2
        target = BACKBONE[:middle] + "AGCTAGCTAGCT" + BACKBONE[middle + 12:]
        result = design_one(rec("pLL057P", BACKBONE), rec("pBAD", target), DATE)
        self.assertNotEqual(result.status, STATUS_OK)
        self.assertFalse(result.verified)
        combined = " ".join(result.messages + result.warnings)
        self.assertIn("repeat", combined.lower())
        self.assertTrue(any("No usable USER junction" in m for m in result.messages))

    def test_repeat_overlap_is_warned_even_when_design_succeeds(self):
        # Just outside the repeat, so a design is possible, but the flank still touches it.
        blocks = [b for b in repeat_blocks(BACKBONE) if b.length > 400]
        pos = blocks[0].end + 6
        target = BACKBONE[:pos] + "GGGTTTAAACCC" + BACKBONE[pos + 12:]
        result = design_one(rec("pLL057P", BACKBONE), rec("pEDGE", target), DATE)
        self.assertTrue(result.template_repeats)
        if result.status == STATUS_OK:
            self.assertTrue(result.verified)


class TestNoDifference(unittest.TestCase):
    def test_identical_sequences_are_reported_not_designed(self):
        template = rec("pLL057P", BACKBONE)
        target = rec("pLL082P", BACKBONE)
        result = design_one(template, target, DATE)
        self.assertEqual(result.status, STATUS_NO_DIFFERENCE)
        self.assertEqual(result.primers, [])
        self.assertFalse(result.verified)
        self.assertTrue(any("identical" in m for m in result.messages))

    def test_rotated_identical_sequences_also_report_no_difference(self):
        template = rec("a", BACKBONE)
        target = rec("b", rotate(BACKBONE, 2000))
        result = design_one(template, target, DATE)
        self.assertEqual(result.status, STATUS_NO_DIFFERENCE)


class TestComparison(unittest.TestCase):
    def test_edits_are_merged_not_fragmented(self):
        pos = 4000
        target = BACKBONE[:pos] + "GATTACAGATTACA" + BACKBONE[pos + 20:]
        comparison = compare_plasmids("t", BACKBONE, "g", target)
        self.assertEqual(len(comparison.edits), 1,
                         f"expected one merged edit, got {[e.describe() for e in comparison.edits]}")

    def test_conserved_blocks_cover_everything_outside_edits(self):
        target = BACKBONE[:3000] + "AAGGTTCC" + BACKBONE[3008:]
        comparison = compare_plasmids("t", BACKBONE, "g", target)
        covered = sum(size for _, _, size in comparison.conserved_blocks())
        edited = sum(e.new_length for e in comparison.edits)
        self.assertEqual(covered + edited, len(comparison.target))
        for tpl_start, tgt_start, size in comparison.conserved_blocks():
            self.assertEqual(
                comparison.template[tpl_start:tpl_start + size],
                comparison.target[tgt_start:tgt_start + size],
            )


class TestSequenceHelpers(unittest.TestCase):
    def test_revcomp_roundtrip(self):
        self.assertEqual(revcomp(revcomp(BACKBONE)), BACKBONE)

    def test_circular_equal_detects_rotation(self):
        equal, rot = circular_equal(rotate(BACKBONE, 777), BACKBONE)
        self.assertTrue(equal)
        self.assertEqual(rot, 777)

    def test_circular_equal_rejects_different_sequences(self):
        equal, _ = circular_equal(BACKBONE, BACKBONE[:-1] + "A" if BACKBONE[-1] != "A" else BACKBONE[:-1] + "C")
        self.assertFalse(equal)

    def test_reader_rejects_duplicate_names(self):
        import csv

        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "dup.csv")
            with open(path, "w", newline="") as fh:
                w = csv.writer(fh)
                w.writerow(["pA", BACKBONE[:200]])
                w.writerow(["pA", BACKBONE[:200]])
            with self.assertRaises(ValueError):
                read_sequence_table(path)

    def test_reader_handles_header_and_swapped_columns(self):
        import csv

        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "swapped.csv")
            with open(path, "w", newline="") as fh:
                w = csv.writer(fh)
                w.writerow(["sequence", "plasmid"])
                w.writerow([BACKBONE[:300], "pX"])
            records = read_sequence_table(path)
            self.assertEqual([r.name for r in records], ["pX"])
            self.assertEqual(records[0].seq, BACKBONE[:300])


def write_pair_csv(path: str, rows: "list") -> None:
    import csv

    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh)
        for name, seq in rows:
            writer.writerow([name, seq])


class TestCli(unittest.TestCase):
    def test_list_shows_names_lengths_and_repeat_map(self):
        proc = subprocess.run(
            [sys.executable, "design_user_primers.py", os.path.basename(INPUT_XLSX), "--list"],
            cwd=REPO, capture_output=True, text=True,
        )
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("pLL057P", proc.stdout)
        self.assertIn("8,116 bp", proc.stdout)
        self.assertIn("repeat", proc.stdout)

    def test_list_flags_duplicate_sequences(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "dupes.csv")
            write_pair_csv(path, [("pA", BACKBONE), ("pB", BACKBONE)])
            proc = subprocess.run(
                [sys.executable, "design_user_primers.py", path, "--list"],
                cwd=REPO, capture_output=True, text=True,
            )
            self.assertEqual(proc.returncode, 0, proc.stderr)
            self.assertIn("identical sequences", proc.stdout)

    def test_real_input_file_designs_and_verifies(self):
        """Guards the actual deliverable, without pinning the file's current contents."""
        records = read_sequence_table(INPUT_XLSX)
        self.assertGreaterEqual(len(records), 2)
        template, target = records[0], records[1]
        with tempfile.TemporaryDirectory() as tmp:
            proc = subprocess.run(
                [sys.executable, "design_user_primers.py", INPUT_XLSX,  # noqa: E501
                 "--template", template.name, "--target", target.name,
                 "--outdir", tmp, "--date", DATE],
                cwd=REPO, capture_output=True, text=True,
            )
            output = proc.stdout + proc.stderr
            if template.seq == target.seq:
                self.assertEqual(proc.returncode, 3, output)
                self.assertIn("SKIPPED", proc.stdout)
            else:
                self.assertEqual(proc.returncode, 0, output)
                self.assertIn("VERIFIED", proc.stdout)
                stamp = f"{DATE.replace('-', '')}_b1"
                predicted = os.path.join(
                    tmp, DATE, "b1", f"{target.name}_from_{template.name}",
                    f"predicted_{target.name}_from_{template.name}_{stamp}.fasta",
                )
                with open(predicted) as fh:
                    seq = "".join(l.strip() for l in fh if not l.startswith(">"))
                equal, _ = circular_equal(seq, target.seq)
                self.assertTrue(equal, "predicted plasmid must match the requested target")

    def test_end_to_end_writes_outputs(self):
        import csv

        with tempfile.TemporaryDirectory() as tmp:
            src = os.path.join(tmp, "batch.csv")
            mutant = BACKBONE[:5500] + "AGCTAGCTAGCT" + BACKBONE[5512:]
            with open(src, "w", newline="") as fh:
                w = csv.writer(fh)
                w.writerow(["pLL057P", BACKBONE])
                w.writerow(["pTEST", mutant])
            out = os.path.join(tmp, "designs")
            proc = subprocess.run(
                [sys.executable, "design_user_primers.py", src,
                 "--template", "pLL057P", "--target", "pTEST",
                 "--outdir", out, "--date", DATE],
                cwd=REPO, capture_output=True, text=True,
            )
            self.assertEqual(proc.returncode, 0, proc.stdout + proc.stderr)
            self.assertIn("VERIFIED", proc.stdout)
            stamp = f"{DATE.replace('-', '')}_b1"
            folder = os.path.join(out, DATE, "b1", "pTEST_from_pLL057P")
            for expected in [
                f"report_pTEST_from_pLL057P_{stamp}.md",
                f"primers_pTEST_from_pLL057P_{stamp}.csv",
                f"order_pTEST_from_pLL057P_{stamp}.tsv",
                f"design_pTEST_from_pLL057P_{stamp}.json",
                f"predicted_pTEST_from_pLL057P_{stamp}.fasta",
            ]:
                self.assertTrue(os.path.exists(os.path.join(folder, expected)), expected)
            for expected in [f"summary_{stamp}.csv", f"primers_{stamp}.csv", f"order_{stamp}.tsv"]:
                self.assertTrue(
                    os.path.exists(os.path.join(out, DATE, "b1", expected)), expected)
            with open(os.path.join(folder, f"primers_pTEST_from_pLL057P_{stamp}.csv")) as fh:
                rows = list(csv.DictReader(fh))
            self.assertEqual(len(rows), 2)
            for row in rows:
                self.assertEqual(row["design_date"], DATE)
                self.assertEqual(row["batch"], "b1")
                self.assertIn(stamp, row["primer_name"])
                self.assertIn(U_MARK, row["order_sequence"])
                self.assertTrue(row["user_junction_sequence"])
                # sequence_with_U is the same oligo with a plain U at the dU position.
                self.assertEqual(row["order_sequence"].replace(U_MARK, "U"),
                                 row["sequence_with_U"])
                self.assertEqual(row["sequence_with_U"].count("U"), 1)
                self.assertEqual(row["sequence_with_U"].replace("U", "T"),
                                 row["plain_sequence"])

    def test_identical_input_exits_with_code_three(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "same.csv")
            write_pair_csv(path, [("pSAME_A", BACKBONE), ("pSAME_B", BACKBONE)])
            proc = subprocess.run(
                [sys.executable, "design_user_primers.py", path,
                 "--template", "pSAME_A", "--target", "pSAME_B",
                 "--outdir", os.path.join(tmp, "out"), "--date", DATE],
                cwd=REPO, capture_output=True, text=True,
            )
            self.assertEqual(proc.returncode, 3, proc.stdout + proc.stderr)
            self.assertIn("SKIPPED", proc.stdout)


if __name__ == "__main__":
    unittest.main(verbosity=2)
