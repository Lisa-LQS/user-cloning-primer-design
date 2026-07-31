# USER cloning primer design pipeline

Takes a spreadsheet of named circular plasmid sequences, works out how to build one
plasmid from another by uracil-excision (USER) cloning, and writes date-stamped,
orderable primers — but only after simulating the whole workflow and confirming it
reproduces the requested plasmid base for base.

## Install

```bash
git clone git@github.com:Lisa-LQS/user-cloning-primer-design.git
cd user-cloning-primer-design
pip install -r requirements.txt      # primer3-py, openpyxl
```

Nothing else is needed (Python 3.8+).

## Example inputs

Two worked examples ship with the repo, one per design mode:

`USER_design_1.xlsx` ships with the repo: two plasmids, pLL057P and pLL082P, differing by a
27 bp deletion that removes an ATG-NLS-stop cassette.

```bash
python3 design_user_primers.py USER_design_1.xlsx --list
python3 design_user_primers.py USER_design_1.xlsx --template pLL057P --target pLL082P
```

That should end with `VERIFIED: in-silico assembly reproduces the target exactly`.

To assemble a target from two or more plasmids, put them all in one sheet and name them:

```bash
python3 design_user_primers.py your_designs.xlsx \
    --sources pBACKBONE pDONOR --target pNEW
```

## Use

```bash
# What is in the file?
python3 design_user_primers.py USER_design_1.xlsx --list

# One construct
python3 design_user_primers.py USER_design_1.xlsx --template pLL057P --target pLL082P

# A later batch: everything in the file built from the same backbone
python3 design_user_primers.py 2026-08_designs.xlsx --all-from pLL057P

# Explicit pairs
python3 design_user_primers.py in.xlsx --pairs pLL057P:pLL082P pLL057P:pLL083P

# Assemble one target from several source plasmids
python3 design_user_primers.py your_designs.xlsx \
    --sources pBACKBONE pDONOR --target pNEW
```

Plasmid names are matched loosely on alphanumerics, so `pMY-v5p` finds `pMYV5P`.

The input may be `.xlsx`, `.csv`, `.tsv` or `.fasta`. For spreadsheets, one row per
plasmid; column order does not matter (the sequence column is detected as the DNA-looking
one) and a header row is skipped automatically.

### Output: date **and** batch

Results are keyed on `designs/<date>/<batch>/`. The date alone is not enough — several
batches often arrive on the same day — so every run gets a batch label as well.

- Without `--batch`, the next free `b1`, `b2`, `b3`, … for that date is used.
- With `--batch mycassette` (any label), that label is used.
- **Existing results are never overwritten or deleted.** Naming an existing batch is an
  error (exit 2) rather than a silent overwrite; without a label the run moves to the next
  free one.

```
designs/2026-07-30/
├── b1/
│   ├── summary_20260730_b1.csv                one row per design in the batch
│   ├── primers_20260730_b1.csv                every primer in the batch, one table
│   ├── order_20260730_b1.tsv                  paste-into-IDT form for the whole batch
│   └── pLL082P_from_pLL057P/
│       ├── report_pLL082P_from_pLL057P_20260730_b1.md     design + protocol
│       ├── primers_pLL082P_from_pLL057P_20260730_b1.csv   full table with all metrics
│       ├── order_pLL082P_from_pLL057P_20260730_b1.tsv     order form
│       ├── design_pLL082P_from_pLL057P_20260730_b1.json   machine-readable everything
│       └── predicted_..._20260730_b1.fasta                the plasmid this should give
└── b2/
    └── pNEW_from_pBACKBONE+pDONOR/ ...
```

Date and batch appear in the directory, every file name, the design ID
(`USER-20260730-b1-pLL082P-01`) and every oligo name (`pLL082P_J1_F_20260730_b1`), so tubes
from different batches cannot be confused.

Each primer appears in three notations: `order_sequence` (`ACGT/ideoxyU/GGCA`, for pasting
into a vendor order), `sequence_with_U` (`ACGTUGGCA`, for reading and for sequence editors)
and `plain_sequence` (all T, for BLAST and alignment).

`primers_<date>_<batch>.csv` includes a **`user_junction_sequence`** column: the
single-stranded 3' overhang USER exposes at that junction, written as the top strand reads
it 5'→3'. Both primers of a junction share the value, so it identifies which oligos pair up;
`three_prime_overhang_this_end` gives the overhang that particular primer's fragment end
presents, and `amplified_from` names the source plasmid for that primer's PCR.

Exit codes: `0` success, `1` a design failed or did not verify, `2` bad input file or a
batch label that is already in use, `3` template and target were identical so there was
nothing to design.

## How the design works

USER cloning primers carry a single deoxyuridine a few bases in from the 5' end. After PCR
with a uracil-tolerant polymerase, USER enzyme (UDG + Endo VIII) excises the uracil and
nicks the backbone, releasing the short 5' flap of each strand and leaving a 3'
single-stranded overhang on the opposite strand. Fragments with complementary overhangs
anneal and are transformed directly — no ligase, no restriction sites.

For a fragment amplified with

```
forward primer = tail_F + extra + annealing region      (U replaces the last base of tail_F)
reverse primer = tail_R + extra + annealing region
```

the post-USER 3' overhangs are `revcomp(tail_R)` on the top strand and `revcomp(tail_F)` on
the bottom. Those can only anneal if `tail_R == revcomp(tail_F)`, so each junction is
defined by a single overhang block, and — since both tails must end in the T that becomes
the uracil — **the overhang block must begin with A and end with T**. That one constraint
drives the junction search.

A useful consequence: every primer is just a contiguous circular slice of the *desired
target* sequence with one T swapped for U. The pipeline exploits this both to build primers
and to check them.

### Steps

1. **Plan the assembly.** With several sources, the one explaining most of the target
   becomes the backbone; each difference against it is then looked up in the other
   plasmids. A difference whose new sequence is a contiguous block of another plasmid is
   amplified from it as its own fragment, rather than being written into primer tails.
2. **Align** backbone and target as circles: pick the target's orientation, rotate it onto
   the template, then move the linear origin into a long conserved block so no difference
   straddles it. Unequal-length plasmids get different rotations, chosen so both land on
   the same base of the shared anchor.
3. **Diff** them and merge nearby differences into single edit sites.
4. **Place one junction per edit site.** Every `(position, length)` overhang block of the
   form `A…T` within a window around the edit is enumerated, scored, and the best kept.
   New sequence too long for the block is carried as non-templated "extra" on the primers,
   split between the pair.
5. **Size each annealing region** to hit the Tm target. Only regions with exactly one
   binding site in the circular template (both strands) are accepted; if nothing in the
   normal 18–32 nt range is unique — which happens inside duplicated promoters or LTRs —
   the region is extended up to 60 nt to reach out of the repeat, and if that fails the
   site is refused with an explanation rather than given unusable primers.
6. **Pair primers into fragments.** One edit site gives a single whole-plasmid product that
   self-circularises; several give one fragment per junction, assembled directionally. Each
   fragment carries its own source plasmid, so a fragment amplified off a donor is primed
   and verified against that donor. Overhangs of different junctions are checked for
   cross-annealing so the assembly cannot go together the wrong way.
7. **Verify.** Simulate PCR off the template, excise the uracils, check that adjoining
   overhangs are complementary, anneal the fragments, and confirm the resulting circle is
   identical to the requested target on both strands. Nothing is reported as usable unless
   this passes.

Scoring balances annealing Tm against the target and against its partner, overhang length
and duplex stability, 3'-end composition, homopolymer runs, primer length, and the length
difference between the pair.

Hairpin and self-dimer predictions are deliberately **not** used. They flag far more
primers than actually fail in a USER reaction, and letting them into the score pushed the
design towards worse choices on everything that does matter. The one complementarity check
that remains is between the overhangs of *different* junctions, because that determines
whether a multi-fragment assembly can go together the wrong way.

### Repeat awareness

Plasmids with duplicated elements have regions where no specific primer exists. The
pipeline maps them (`--list` prints them) and reports which stretches are uniquely primable,
so an impossible request is recognised as impossible rather than answered with primers that
would amplify two products.

## Tuning

`--overhang-min/--overhang-max` (default 8–12), `--anneal-min/--anneal-max` (18–32),
`--tm-target` (62 °C), `--junction-slack` (how far from an edit the overhang block may sit,
30 nt), `--max-primer-len` (60 nt, the length above which an Ultramer is flagged).
Programmatic users can pass any `DesignParams` field.

## Layout

| File | Role |
| --- | --- |
| `design_user_primers.py` | CLI |
| `user_cloning/sequences.py` | sequence helpers, circular arithmetic, file reading, repeat finding |
| `user_cloning/plasmid_diff.py` | circular alignment and edit detection |
| `user_cloning/thermo.py` | melting temperatures via primer3 |
| `user_cloning/assembly.py` | which source plasmid contributes which part of the target |
| `user_cloning/design.py` | junction placement and primer construction |
| `user_cloning/simulate.py` | in-silico PCR, USER digestion, assembly, verification |
| `user_cloning/pipeline.py` | orchestration |
| `user_cloning/report.py` | CSV, TSV, JSON, FASTA and Markdown output |
| `tests/test_pipeline.py` | test suite |

## Tests

```bash
python3 -m unittest discover -s tests
```

The suite designs point mutations, substitutions, insertions (up to 132 nt), small and
large deletions, unequal-length replacements, and two- and three-junction assemblies
against the real 8,116 bp pLL057P backbone, plus synthetic two-source assemblies (cassette
moved between plasmids, backbone selection, source ordering), — checking each time both that the pipeline's
own verification passes and that the predicted plasmid independently matches the requested
target. It also covers rotated and reverse-complemented input, edits sitting on the
sequence origin, identical input, and edits inside repeated elements.

## Wet-lab notes

- **No restriction enzyme is involved at any point.** The overhangs come from uracil
  excision, so the design never needs a restriction site to be present, introduced or
  preserved. The only enzymes are the polymerase and USER.
- Use **Phusion U Hot Start**, **Q5U** or **PfuTurbo Cx**. Ordinary proofreading enzymes
  stall at the deoxyuridine.
- USER: 1 U per ~0.1–0.2 pmol purified product, 37 °C 25 min, then anneal. The report picks
  25 °C or a ramp to 10 °C based on the actual overhang duplex Tm.
- Transform directly; the nicks are sealed in vivo. No ligase.
- Sequence across each junction to confirm.
