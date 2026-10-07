# EGFR conditional-binder sequence panel: methodology and status

## Scope and sequence sources

This is a sequence-level candidate panel for the EGFR conditional-binder challenge, not a set of validated binders. The repository contains full-length human and mouse EGFR FASTAs (`data/sequences/human_egfr.fasta`, `data/sequences/mouse_egfr.fasta`). The existing `scripts/align_egfr.py` aligns their extracellular segments (residues 25?645). The working Domain III window used here is residues 311?480 in the aligned full-length sequences:

- Human EGFR 311?480 (170 aa): `CGADSYEMEEDGVRKCKKCEGPCRKVCNGIGIGEFKDSLSINATNIKHFKNCTSISGDLHILPVAFRGDSFTHTPPLDPQELDILKTVKEITGFLLIQAWPENRTDLHAFENLEIIRGRTKQHGQFSLAVVSLNITSLGLRSLKEISDGDVIISGNKNLCYANTINWKKL`
- Mouse EGFR homologous 311?480 (170 aa): `CGPDYYEVEEDGIRKCKKCDGPCRKVCNGIGIGEFKDTLSINATNIKHFKYCTAISGDLHILPVAFKGDSFTRTPPLDPRELEILKTVKEITGFLLIQAWPDNWTDLHAFENLEIIRGRTKQHGQFSLAVVGLNITSLGLRSLKEISDGDVIISGNRNLCYANTINWKKL`
- Exact identity in this gap-free window: 154/170 (90.6%); the 16 sequence differences include human S348?mouse T and human S442?mouse G.

## Target region and cross-species rationale

The target hypothesis is the Domain III face contacted by EGF in the existing human EGFR?EGF 1IVO structure. The repository's contact table and validated mappings identify conserved human/mouse interface residues including D379, S380, F381, T382, L406, Q408, Q432, H433, Q435, F436, A439 and V441. The human 1IVO structural mapping places human EGFR 348 at chain A SER 324, near the interface (closest listed EGF context: chain C TYR 44), and human 442 at chain A SER 418, with a direct geometric contact record to chain C ARG 45. Human 348 and 442 themselves differ in mouse, so the intended shared-species recognition surface is the conserved surrounding patch; the design should not rely solely on the divergent S442 side chain. These mapping and contact statements come from the existing repository outputs; no mouse structure is implied.

## Binder format and sequence strategy

Eight unique 59-aa sequences are supplied in one coherent format: a manually authored single-chain, three-helix mini-protein hypothesis. The proposed scaffold uses short gly/pro/ser turns and a leucine/isoleucine-rich amphipathic core pattern; helix 3 carries a varied exposed interaction face. These sequences were authored de novo for this panel and are not copied therapeutic antibody sequences. The fold itself has not been predicted or experimentally checked.

The exposed helix-3 face varies histidine number/placement and aromatic/acidic residues. The proposed pH mechanism is conditional: at pH 6.5, partial histidine protonation could strengthen electrostatic/H-bond interactions with the conserved acidic EGFR D379 region; at pH 7.4, a smaller protonated fraction could weaken that contribution. Conserved Q408/H433/Q435/F436/V441 are candidate context for shared human/mouse recognition, while Tyr/Trp/Phe are included as possible affinity anchors. No docking establishes that these residues face or contact the proposed paratope. Histidine incorporation alone does not guarantee pH selectivity; a sharp ?detectable at pH 6.5 / undetectable at pH 7.4? window requires empirical screening and likely optimization.

## Validation and unresolved requirements

The CSV and FASTA contain 8 distinct candidates, all 59 aa, using only the 20 canonical amino acids, with one molecule class. There are no exact full-sequence matches to the two local EGFR FASTAs. The competition novelty score is platform-specific and cannot be calculated from repository data; its required threshold (?3/4) is therefore **unverified**. Local duplication checks do not establish global novelty.

No binder-design, folding, docking, pH-dependent scoring, or affinity-prediction model is available in this repository environment. Consequently, these sequences have **no computational binding prediction** and no experimental validation. Do not describe them as binders or claim human/mouse cross-reactivity, pH selectivity, affinity, or therapeutic effect. They are design hypotheses only. Given the one-submission-per-24-hours limit, verify novelty and run the competition's supported prediction/scoring workflow before submitting; the current panel is not platform-validated or submission-certified.
