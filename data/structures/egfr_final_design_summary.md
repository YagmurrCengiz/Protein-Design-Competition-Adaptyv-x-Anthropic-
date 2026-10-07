# EGFR final design evidence summary

## 1. Design rationale

Positions 442 and 348 were selected for structural modeling because existing 1IVO mapping places human EGFR 442 at chain A SER 418, with a direct geometric contact record to EGF chain C ARG 45 (4.439 A in WT; source interaction typing is LOW confidence), while human EGFR 348 maps to chain A SER 324 and is classified as near a 3D interface. The 348 site is 5.640 A from EGF chain C TYR 44 in WT and has reported nearby EGFR interface-environment partners, but no direct EGF contact record. These are geometric/contextual observations, not proof of a functional effect.
Evidence strength and interface class in the table are inherited from the existing site-level synthesis; they were not rescored per substitution.

## 2. 442 mutation panel

The S>G, S>A, and S>T variants form a controlled side-chain chemistry and geometry perturbation. G removes the hydroxyl and side-chain bulk; A removes the hydroxyl but retains a small methyl group; T retains the hydroxyl and adds methyl bulk. Comparing the three can test whether observed effects track hydroxyl presence, side-chain volume, or both. Existing model distances are 4.439 A WT to 5.626 A for G, 5.626 A for A, and 4.286 A for T. These are minimum heavy-atom distances to EGF chain C and do not measure binding.

## 3. 348 S>T

Human EGFR 348 S>T is a near-interface candidate. In the selected THR rotamer model, the minimum heavy-atom distance to EGF chain C remained 5.640 A, equal to WT at the reported precision. This does not rule out local packing or functional changes and does not establish a direct EGF contact.

## 4. 255 R>Q comparator

Human EGFR 255 R>Q is retained as a local-environment comparator: the source reports 14 1IVO neighbors within 5 A, including 7 charged and 5 polar neighbors, and a loss of positive charge. Four potential salt-bridge contexts are heuristic counts, not confirmed salt bridges. The site is classified NONE for interface relationship and no tested interface evidence is reported. No mutant geometry was modeled, so there is no WT-vs-mutant distance comparison.

## 5. Structural-model interpretation

The models are computational side-chain replacement models, not experimentally observed structures. No molecular dynamics was performed, no energetic stability was calculated, no binding affinity was calculated, and no biological causality can be inferred. Distance changes are geometric observations only. The 348 S>T model uses the first of three THR library rotamers; PyMOL's optional bump-scoring step was disabled after a batch selection-parser error, so the chosen rotamer was not ranked by that step. The 442 S>A and S>G PDBs are side-chain truncation/identity models without rotamer repacking. All four passed structure parsing, residue identity, duplicate-atom, and repository clash checks.

## 6. Final testable hypotheses

- **442 S>G:** Test whether removal of both the serine hydroxyl and side-chain bulk changes EGF binding or ligand-stimulated EGFR phosphorylation relative to WT.
- **442 S>A:** Test whether loss of the hydroxyl while retaining a methyl group changes EGF binding or ligand-stimulated phosphorylation relative to WT and the S>G/S>T variants.
- **442 S>T:** Test whether retaining a hydroxyl while adding methyl bulk changes EGF binding or ligand-stimulated phosphorylation relative to WT and the S>G/S>A variants.
- **348 S>T:** Test whether added methyl bulk at the near-interface site changes receptor behavior or ligand-stimulated phosphorylation relative to WT, despite unchanged minimum EGF-C distance in the selected model.
- **255 R>Q:** Test whether neutralizing Arg255 changes receptor surface abundance or ligand-stimulated phosphorylation independently of direct interface effects.

These are proposed experiments, not claims that any hypothesis is true.
