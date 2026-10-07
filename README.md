# EGFR Protein Design Competition

Computational project for the Proteinbase competition, focused on EGFR Domain III and conditional binder design.

## Project status

- **Week 1 complete:** sequence comparison, structure/interface analysis, mutation prioritization, candidate evidence, and an initial sequence submission panel.
- **Week 2:** not started. Start only after reviewing the Week 1 handoff and confirming the next module's requirements.
- **Evidence level:** sequence and structure analyses are computational. Candidate binding, folding, cross-species recognition, and pH response remain unvalidated unless explicitly stated otherwise.

## Start here

- [Week 1 recap and handoff](docs/week-01/README.md)
- [Script map and execution notes](scripts/README.md)
- [Original competition methodology](data/structures/EGFR_competition_methodology.md)
- [Current submission panel (v2)](data/structures/EGFR_competition_submission_v2.fasta)
- [Submission diversity report](data/structures/EGFR_competition_diversity_report.md)

## Repository layout

| Path | Purpose |
|---|---|
| `data/sequences/` | Human and mouse EGFR FASTA inputs |
| `data/structures/` | Structural inputs, analysis tables, visualizations, and submission artifacts |
| `scripts/` | Analysis and modeling scripts; see its README for workflow order |
| `designs/` | Reserved for generated, filtered, and final design files |
| `predictions/` | Reserved for prediction outputs |
| `results/` | Reserved for consolidated results |
| `logs/` | Reserved for run logs |

Keep running scripts from the repository root unless a script explicitly documents another working directory; several existing scripts use root-relative input/output paths.

## Environment

Python dependencies are listed in `requirements.txt`. PyMOL is only needed for the optional mutant-modeling workflow. Do not interpret a generated model as experimental structural evidence.
