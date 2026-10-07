# Script map

Run commands from the repository root, for example `python scripts/align_egfr.py`. A number of scripts refer to `data/...` with paths relative to the current working directory.

## Week 1 workflows

### 1. EGFR sequence inputs and conservation

- `download_egfr.py` — fetches the EGFR reference sequences.
- `align_egfr.py` — aligns human and mouse EGFR.
- `analyze_conservation.py` — summarizes sequence conservation.

Inputs and sequence outputs belong in `data/sequences/`.

### 2. Structure and ligand-interface analysis

- `analyze_structure.py` — compares EGFR structures and produces the initial residue-level structure comparison.
- `analyze_1ivo_chains.py` — inspects chains in the 1IVO structure.
- `analyze_egf_contacts.py` — measures EGFR–EGF contacts.
- `analyze_egfr_egf_interface.py` — analyzes the EGFR–EGF interface table.
- `analyze_candidate_domains.py` — assigns candidate residues to domains.
- `analyze_candidate_interactions.py` — summarizes candidate/interface interactions.
- `analyze_high_priority_candidates.py` — prioritizes candidate positions.

### 3. Mutation environment and prioritization

Typical dependency order:

1. `analyze_mutations.py`
2. `analyze_mutation_properties.py` and `prioritize_mutations.py`
3. `analyze_mutation_neighbors.py`
4. `analyze_mutation_environment.py`
5. `integrate_mutation_analysis.py`

`analyze_interface_differences.py` and `analyze_interface_environment.py` add interface-specific context. Several outputs are inputs to later scripts, so keep their filenames and root-relative paths stable.

### 4. Evidence review and final candidate summaries

- `analyze_candidate_mechanistic_evidence.py`
- `build_candidate_mechanistic_comparison.py`
- `review_candidate_mechanistic_evidence.py`
- `build_final_candidate_evidence.py`
- `build_final_candidate_evidence_synthesis.py`

These scripts combine and review analysis outputs. Their generated evidence tables are not binding predictions.

### 5. Optional mutant modeling

- `model_egfr_mutants_pymol.py` — generates side-chain mutant models using PyMOL.
- `validate_egfr_mutant_models.py` — checks the generated structures and records validation results.

These are computational models, not experimentally observed structures.

## Practical notes

- Keep source inputs, derived tables, models, and final submission files distinguishable by their names; generated structure artifacts currently live together under `data/structures/` because scripts expect those paths.
- Do not move or rename analysis outputs without updating every downstream path listed in the scripts.
- Before rerunning an analysis, check the input/output constants near the top and bottom of its script. Some older scripts use fixed filenames and may overwrite a previous output.
- Week 2 scripts and analyses have not been added to this map yet.
