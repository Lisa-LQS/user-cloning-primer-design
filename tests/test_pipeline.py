"""Tests for the USER cloning pipeline.

Every design case is checked two ways: the pipeline's own in-silico verification must
pass, and the predicted plasmid must independently match the requested target as a
circular sequence. The real 8,116 bp pLL057P backbone is used as the template so the
tests exercise realistic sequence composition, repeats included.
"""

from __future__ import annotations

import os
import random
import subprocess
import sys
import tempfile
import unittest
from dataclasses import replace

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)

from user_cloning import DesignParams, SeqRecord, compare_plasmids, read_sequence_table  # noqa: E402
from user_cloning.design import U_MARK, pair_tm_penalty  # noqa: E402
from user_cloning.report import allocate_batch, primer_rows  # noqa: E402
from user_cloning.pipeline import (  # noqa: E402
    ANNEAL_TEMP_C,
    EXTENSION_TEMP_C,
    STATUS_NO_DIFFERENCE,
    STATUS_OK,
    _protocol,
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


def point_mutant(pos: int = 4000) -> str:
    """The backbone with one base changed. The replacement is chosen against the base
    actually at `pos`, because picking a fixed letter can silently reproduce the backbone
    and leave a test asserting over an empty design."""
    original = BACKBONE[pos]
    new = "A" if original != "A" else "C"
    return BACKBONE[:pos] + new + BACKBONE[pos + 1:]


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



class TestPrimingRegion(DesignCaseMixin, unittest.TestCase):
    """Tm is scored on everything that pairs in cycle 1, not on the 3' region alone.

    A USER junction is placed in sequence the two parents share, so the 5' block carrying
    the dU is normally plain template: it base-pairs from the first cycle and has to be
    counted, or the reported Tm understates the real duplex and Ta comes out too low.
    """

    def designed(self, target_seq: str):
        """Run a design and refuse to hand back an empty one.

        Most assertions here loop over `result.primers`, so a design that quietly produces
        nothing would make them all pass without testing anything.
        """
        result = self.run_design(target_seq)
        self.assertEqual(result.status, STATUS_OK, msg="; ".join(result.messages))
        self.assertTrue(result.primers, "design produced no primers")
        self.assertTrue(result.fragments, "design produced no fragments")
        return result

    def test_priming_region_is_a_templated_suffix_of_the_primer(self):
        result = self.designed(point_mutant())
        for primer in result.primers:
            self.assertTrue(primer.prime_region.endswith(primer.anneal))
            self.assertEqual(primer.sequence[len(primer.sequence) - len(primer.prime_region):],
                             primer.prime_region)
            self.assertGreaterEqual(len(primer.prime_region), len(primer.anneal))
            # Whatever we call templated must really occur in the template.
            self.assertGreaterEqual(
                count_circular_occurrences(BACKBONE, primer.prime_region, both_strands=True), 1,
                f"{primer.name} priming region does not occur in the template",
            )

    def test_priming_tm_is_at_least_the_annealing_tm(self):
        for primer in self.designed(point_mutant()).primers:
            self.assertGreaterEqual(primer.prime_tm, primer.anneal_tm - 1e-9)
            if len(primer.prime_region) > len(primer.anneal):
                self.assertGreater(primer.prime_tm, primer.anneal_tm)

    def test_tail_templated_means_the_du_base_itself_pairs(self):
        """`tail_templated` is about the dU base, not the whole 5' block.

        The flag says the templated run reaches back past every non-templated base and
        covers the dU, which sits at the end of the tail. The bases 5' of the dU may still
        hang off unpaired -- in a deletion, for instance, the overhang block is shared
        sequence in the *target* but comes from only one side of the *template*, so one
        primer of the pair carries a block the template only partly matches.
        """
        result = self.designed(BACKBONE[:3000] + BACKBONE[3120:])
        for primer in result.primers:
            self.assertEqual(primer.extra, "", "a deletion needs no non-templated insert")
            overhang_bases = len(primer.prime_region) - len(primer.anneal)
            if primer.tail_templated:
                self.assertGreaterEqual(
                    overhang_bases, len(primer.extra) + 1,
                    f"{primer.name} is flagged templated but its dU is outside the duplex",
                )
            else:
                self.assertLess(overhang_bases, len(primer.extra) + 1)
        self.assertTrue(any(p.tail_templated for p in result.primers),
                        "at least one primer of a deletion pair should pair through its dU")

    def test_priming_tm_lands_in_the_target_window(self):
        params = DesignParams()
        for primer in self.designed(point_mutant()).primers:
            self.assertGreaterEqual(primer.prime_tm, params.tm_min - 4)
            self.assertLessEqual(primer.prime_tm, params.tm_max + 4)

    def test_protocol_anneals_at_64_and_extends_at_68(self):
        """Primers are designed to 64 C and the bench program anneals there, while the
        mastermix still extends at its own 68 C -- so a routine design is 3-step."""
        params = DesignParams()
        self.assertEqual(params.tm_target, 64.0)
        self.assertEqual((params.tm_min, params.tm_max), (61.0, 67.0))

        result = self.designed(point_mutant())
        by_name = {f.name: f for f in result.fragments}
        for protocol in result.protocols:
            fragment = by_name[protocol.fragment]
            limiting = min(fragment.forward.prime_tm, fragment.reverse.prime_tm)
            self.assertAlmostEqual(protocol.limiting_primer_tm_C, round(limiting, 1), places=6)
            self.assertFalse(protocol.two_step)
            self.assertEqual(protocol.annealing_temp_C, ANNEAL_TEMP_C)
            self.assertLess(protocol.annealing_temp_C, EXTENSION_TEMP_C)

    def test_protocol_collapses_to_two_step_if_annealing_at_the_extension_temp(self):
        """The 2-step branch is not dead code: it is what annealing at 68 C would give."""
        fragment = self.designed(point_mutant()).fragments[0]
        two = _protocol(fragment, anneal_temp_C=EXTENSION_TEMP_C)
        self.assertTrue(two.two_step)
        self.assertEqual(two.annealing_temp_C, EXTENSION_TEMP_C)
        three = _protocol(fragment, anneal_temp_C=ANNEAL_TEMP_C)
        self.assertFalse(three.two_step)
        # The limiting primer is a property of the pair, not of the cycling choice.
        self.assertEqual(two.limiting_primer_tm_C, three.limiting_primer_tm_C)


class TestPrimerPairingAcrossJunctions(unittest.TestCase):
    """Tm is matched between the primers that share a tube, which span two junctions.

    Fragment i is amplified by junction i's forward primer and junction i+1's reverse
    (`build_fragments`), the last wrapping back to the first. So with two junctions the
    tubes hold J1_F + J2_R and J2_F + J1_R -- never J1_F + J1_R. Matching Tm within a
    junction optimises a pair that never meets and leaves the real pairs unconstrained.
    """

    def two_junction_design(self):
        donor_payload = "".join(BACKBONE[i] for i in range(600, 1200))
        donor = rec("pDONOR", BACKBONE[3000:4000] + donor_payload + BACKBONE[4000:5000])
        backbone = rec("pBB", BACKBONE)
        target_seq = BACKBONE[:5300] + donor_payload + BACKBONE[5300:]
        result = design_assembly([backbone, donor], rec("pASM", target_seq), DATE)
        self.assertEqual(result.status, STATUS_OK, "; ".join(result.messages))
        self.assertEqual(len(result.fragments), 2)
        return result

    def test_a_fragment_takes_its_primers_from_two_different_junctions(self):
        for fragment in self.two_junction_design().fragments:
            self.assertNotEqual(
                fragment.forward.junction, fragment.reverse.junction,
                f"{fragment.name} would be a junction-internal pair, which cannot happen "
                "in a multi-junction assembly",
            )

    def test_tm_is_matched_within_each_pcr_not_within_each_junction(self):
        result = self.two_junction_design()
        params = DesignParams()
        for fragment in result.fragments:
            diff = abs(fragment.forward.prime_tm - fragment.reverse.prime_tm)
            self.assertLessEqual(
                diff, params.max_pair_tm_diff,
                f"{fragment.name} ({fragment.forward.name} + {fragment.reverse.name}) "
                f"differs by {diff:.1f} C; these two share an annealing step",
            )

    def test_mismatch_warnings_name_a_fragment_not_a_junction(self):
        result = self.two_junction_design()
        for warning in result.warnings:
            if "priming Tm differs" in warning:
                self.assertRegex(warning, r"^F\d+:")
                self.assertIn("amplifying", warning)

    def test_pair_penalty_is_symmetric_and_zero_for_equal_tm(self):
        result = self.two_junction_design()
        a = result.fragments[0].forward
        b = result.fragments[0].reverse
        self.assertEqual(pair_tm_penalty(a, b), pair_tm_penalty(b, a))
        self.assertEqual(pair_tm_penalty(a, a), 0.0)


class TestReportedRegionsMatchReality(unittest.TestCase):
    """The reported priming region must be the duplex that really forms.

    `template_block_3prime` is only where the design put the junction boundary. The match
    with the source normally runs on past it, so a column claiming those 5' bases are
    non-templated is wrong -- and a Tm computed on the block alone understates the duplex.
    """

    def two_junction_design(self):
        donor_payload = "".join(BACKBONE[i] for i in range(600, 1200))
        donor = rec("pDONOR", BACKBONE[3000:4000] + donor_payload + BACKBONE[4000:5000])
        backbone = rec("pBB", BACKBONE)
        target_seq = BACKBONE[:5300] + donor_payload + BACKBONE[5300:]
        result = design_assembly([backbone, donor], rec("pASM", target_seq), DATE)
        self.assertEqual(result.status, STATUS_OK, "; ".join(result.messages))
        return result

    @staticmethod
    def longest_3prime_match(primer_seq: str, source: str, direction: str) -> int:
        """How far back from the 3' end the primer really matches the source, circularly."""
        doubled = source + source
        best = 0
        for length in range(1, len(primer_seq) + 1):
            probe = primer_seq[-length:]
            hay = doubled if direction == "forward" else revcomp(doubled)
            if probe in hay:
                best = length
        return best

    def test_priming_region_is_the_maximal_true_match(self):
        for fragment in self.two_junction_design().fragments:
            for primer in (fragment.forward, fragment.reverse):
                true_len = self.longest_3prime_match(
                    primer.sequence, fragment.source_seq, primer.direction)
                self.assertEqual(
                    len(primer.prime_region), true_len,
                    f"{primer.name}: reported priming region is "
                    f"{len(primer.prime_region)} nt but the primer really matches "
                    f"{fragment.source_name} over {true_len} nt",
                )
                self.assertEqual(primer.prime_region, primer.sequence[-true_len:])

    def test_the_layout_block_never_overstates_the_duplex(self):
        for fragment in self.two_junction_design().fragments:
            for primer in (fragment.forward, fragment.reverse):
                self.assertGreaterEqual(len(primer.prime_region), len(primer.anneal))
                self.assertGreaterEqual(primer.prime_tm, primer.anneal_tm - 1e-9)

    def test_csv_reports_how_far_the_match_runs_past_the_block(self):
        result = self.two_junction_design()
        for row, primer in zip(primer_rows(result), result.primers):
            self.assertEqual(row["templated_beyond_block_nt"],
                             len(primer.prime_region) - len(primer.anneal))
            # The old column name asserted these bases were non-templated. Whatever the
            # column is called now, it must not make that claim when they do pair.
            self.assertNotIn("non_templated", " ".join(row.keys()))


class TestLengthPenalty(unittest.TestCase):
    """Primer length is charged superlinearly, so very long oligos are not nearly free."""

    @staticmethod
    def cost(length: int, params: DesignParams = None) -> float:
        p = params or DesignParams()
        return (p.length_weight * max(0, length - p.length_free_upto) ** 2 / 10.0
                + 0.6 * max(0, length - p.soft_max_primer_len))

    def test_short_primers_are_free(self):
        for length in (20, 30, 40):
            self.assertEqual(self.cost(length), 0.0)

    def test_cost_accelerates_with_length(self):
        """Each extra 5 nt must cost more than the previous 5 did -- that is the point of
        a squared term over a linear one."""
        steps = [self.cost(n + 5) - self.cost(n) for n in range(40, 60, 5)]
        for earlier, later in zip(steps, steps[1:]):
            self.assertGreater(later, earlier)

    def test_a_long_oligo_outweighs_a_small_tm_mismatch(self):
        """The regression this guards: at the old linear 0.15/nt a 60 nt primer cost 2.25,
        less than a 2 C mismatch within a PCR pair, so length never influenced a choice."""
        two_degrees = 2.0 * 2.0  # pair_tm_penalty for a 2 C difference
        self.assertGreater(self.cost(60), two_degrees)
        self.assertGreater(self.cost(54), two_degrees)

    def test_crossing_the_ultramer_limit_adds_a_step(self):
        below, above = self.cost(60), self.cost(61)
        smooth = (DesignParams().length_weight
                  * (61 - DesignParams().length_free_upto) ** 2 / 10.0)
        self.assertGreater(above - below, 0.0)
        self.assertAlmostEqual(above - smooth, 0.6, places=6)

    def test_weight_zero_disables_the_term(self):
        params = replace(DesignParams(), length_weight=0.0)
        self.assertEqual(self.cost(55, params), 0.0)
        self.assertAlmostEqual(self.cost(70, params), 0.6 * 10)


class TestPartnerAwareResizing(unittest.TestCase):
    """After junctions are fixed, each primer is re-sized against its tube-mate.

    Sizing happens before anything knows which primers share a PCR, so each is first sized
    toward the global Tm target. `refine_pairs` revisits that once the fragments are known.
    """

    def design(self, params: DesignParams = None):
        donor_payload = "".join(BACKBONE[i] for i in range(600, 1200))
        donor = rec("pDONOR", BACKBONE[3000:4000] + donor_payload + BACKBONE[4000:5000])
        target_seq = BACKBONE[:5300] + donor_payload + BACKBONE[5300:]
        result = design_assembly([rec("pBB", BACKBONE), donor], rec("pASM", target_seq),
                                 DATE, params or DesignParams())
        self.assertEqual(result.status, STATUS_OK, "; ".join(result.messages))
        return result

    def test_resizing_cannot_change_the_assembled_plasmid(self):
        """Only the 3' annealing region may move. The overhang block and the carried extra
        define the junction, and the primer's 5' end is where the product begins."""
        result = self.design()
        self.assertTrue(result.verified, "; ".join(result.verification.errors))
        for junction in result.junctions:
            for primer in (junction.forward, junction.reverse):
                self.assertTrue(primer.sequence.startswith(primer.tail + primer.extra),
                                f"{primer.name} 5' end no longer matches its junction")
        # And the independent check: the predicted plasmid is still the requested one.
        donor_payload = "".join(BACKBONE[i] for i in range(600, 1200))
        target_seq = BACKBONE[:5300] + donor_payload + BACKBONE[5300:]
        equal, _ = circular_equal(predicted_plasmid(result), target_seq)
        self.assertTrue(equal)

    def test_pairs_end_up_matched(self):
        params = DesignParams()
        for fragment in self.design().fragments:
            self.assertLessEqual(
                abs(fragment.forward.prime_tm - fragment.reverse.prime_tm),
                params.max_pair_tm_diff,
                f"{fragment.name} pair is still split after re-sizing",
            )

    def test_matching_never_drags_a_primer_out_of_the_usable_window(self):
        """The regression this guards: scored on Tm gap and length alone, the pass once
        picked a pair agreeing at 59.5/60.3 C -- well matched, but below the 61 C floor
        and below the bench annealing step, so neither primer would hold."""
        params = DesignParams()
        for fragment in self.design().fragments:
            for primer in (fragment.forward, fragment.reverse):
                self.assertGreaterEqual(
                    primer.prime_tm, params.tm_min - 0.1,
                    f"{primer.name} at {primer.prime_tm:.1f} C is below the "
                    f"{params.tm_min} C floor",
                )
                self.assertGreaterEqual(
                    primer.prime_tm, ANNEAL_TEMP_C - params.max_pair_tm_diff,
                    f"{primer.name} at {primer.prime_tm:.1f} C cannot hold at the "
                    f"{ANNEAL_TEMP_C} C annealing step",
                )

    def test_pass_is_reported_when_it_acts(self):
        notes = [w for w in self.design().warnings if "re-sized" in w]
        for note in notes:
            self.assertRegex(note, r"^F\d+:")
            self.assertIn("Tm gap", note)


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


class TestInsertTiledAcrossDonors(unittest.TestCase):
    """An insert that is not one block of any single donor, but a join of two.

    pLL057Z is built this way: 607 nt of the cassette comes from pLL217, the next 372 nt
    from pNB0153c_3, and the last 6 nt from neither. Payloads here are synthetic rather
    than backbone slices, so that each one primes uniquely.
    """

    @staticmethod
    def payload(seed: int, length: int) -> str:
        rng = random.Random(seed)
        # Balanced GC and no long homopolymers, so junction placement is never the
        # thing under test here.
        return "".join(rng.choice("ACGT") for _ in range(length))

    def assemble(self, sources, target_seq, name="pTILED"):
        target = rec(name, target_seq)
        result = design_assembly(sources, target, DATE)
        self.assertEqual(result.status, STATUS_OK, "; ".join(result.messages))
        self.assertTrue(result.verified, "; ".join(result.verification.errors))
        equal, _ = circular_equal(predicted_plasmid(result), target_seq)
        self.assertTrue(equal, "predicted plasmid does not match the requested target")
        return result

    def test_two_donors_make_one_insert(self):
        left = self.payload(1, 700)
        right = self.payload(2, 600)
        donor_a = rec("pDA", self.payload(3, 800) + left)
        donor_b = rec("pDB", self.payload(4, 600) + right)
        backbone = rec("pBB", BACKBONE)
        target_seq = BACKBONE[:5300] + left + right + BACKBONE[5300:]
        result = self.assemble([backbone, donor_a, donor_b], target_seq)
        self.assertEqual(
            {f.source_name for f in result.fragments}, {"pBB", "pDA", "pDB"},
            "each donor should contribute its own PCR fragment")
        self.assertEqual(len(result.fragments), 3)

    def test_tiled_insert_with_untemplated_tail(self):
        """The real pLL057Z shape: two donor blocks plus a short non-templated tail."""
        left = self.payload(5, 700)
        right = self.payload(6, 600)
        tail = "GATCTAGA"
        donor_a = rec("pDA", self.payload(7, 800) + left)
        donor_b = rec("pDB", self.payload(8, 600) + right)
        backbone = rec("pBB", BACKBONE)
        target_seq = BACKBONE[:5300] + left + right + tail + BACKBONE[5300:]
        result = self.assemble([backbone, donor_a, donor_b], target_seq)
        self.assertEqual(len(result.fragments), 3)
        # The tail is on nobody's template, so it has to be written into a primer.
        self.assertTrue(
            any(tail in p.sequence.replace(U_MARK, "") for p in result.primers),
            "the non-templated tail must be carried on a primer")

    def test_single_block_insert_is_unchanged(self):
        """Tiling must not disturb the ordinary one-donor case."""
        load = self.payload(9, 700)
        donor = rec("pD", self.payload(10, 800) + load)
        backbone = rec("pBB", BACKBONE)
        target_seq = BACKBONE[:5300] + load + BACKBONE[5300:]
        result = self.assemble([backbone, donor], target_seq)
        self.assertEqual(len(result.fragments), 2)
