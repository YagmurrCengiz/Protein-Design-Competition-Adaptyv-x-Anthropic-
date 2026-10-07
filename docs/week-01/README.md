# Week 1 recap and handoff

## Goal

Build a first evidence-based view of human and mouse EGFR and prepare an initial competition submission hypothesis.

## Work completed

1. Collected and aligned human and mouse EGFR sequence references; used the high conservation of Domain III to motivate a shared-species target hypothesis.
2. Examined EGFR structure and the EGF-facing Domain III surface using the available 1IVO and AlphaFold-derived structural files.
3. Analyzed candidate residues, mutation environments, interface contacts, and local mechanistic context.
4. Generated and validated computational side-chain mutant models for selected EGFR positions. These models are geometric hypotheses, not experimental observations.
5. Prepared the original competition panel and a v2 sequence panel with a diversity comparison report.

## Main files

- Human and mouse sequences: `data/sequences/human_egfr.fasta`, `data/sequences/mouse_egfr.fasta`
- Structure and derived analysis files: `data/structures/`
- Script workflow: [scripts/README.md](../../scripts/README.md)
- Design rationale: `data/structures/EGFR_competition_methodology.md`
- Submission sequences: `data/structures/EGFR_competition_submission_v2.fasta` and `data/structures/EGFR_competition_submission_v2.csv`
- Sequence diversity report: `data/structures/EGFR_competition_diversity_report.md`

## Interpretation limits

The sequence panel is computational and unvalidated. Domain III conservation supports the cross-species design rationale but does not demonstrate binding to either species. Histidine placement gives a protonation-based pH hypothesis, not a guarantee of binding at pH 6.5 or loss of detectable binding at pH 7.4. Local sequence identity comparisons cannot determine the competition platform's global novelty score.

## Handoff to Week 2

Week 2 has not started. Before adding new analyses, record its objective, required inputs, expected outputs, and success criteria here or in a separate Week 2 note. Preserve Week 1 source files and outputs so new work can be traced back to its inputs.
