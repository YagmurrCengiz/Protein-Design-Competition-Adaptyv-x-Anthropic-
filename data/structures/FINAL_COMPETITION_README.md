# EGFR competition handoff

## 1. Problem and design objective

This work compares human and mouse EGFR sequence differences, then uses validated human EGFR structural context to define a compact, experimentally testable five-variant shortlist. The objective is to distinguish direct interface geometry from near-interface and local-environment hypotheses without claiming functional direction.

## 2. Sequence comparison

Human EGFR was compared with mouse EGFR. The selected human positions differ in mouse as shown below. ?Human?mouse? describes the sequence comparison; panel mutations such as 442 S?A and S?T are separate substitutions modeled on the corresponding **human** structure. 1IVO is a human EGFR?EGF structure, not a mouse structure.

| Human EGFR position | Human residue | Mouse residue | Human?mouse substitution | Design mutation | PDB chain | PDB residue | Structural model available? |
|---:|:---:|:---:|:---:|:---:|:---:|---:|:---:|
| 442 | S | G | S?G | S?G, S?A, S?T | A | 418 | Yes, all three human-structure variants |
| 348 | S | T | S?T | S?T | A | 324 | Yes, human-structure variant |
| 255 | R | Q | R?Q | R?Q | A | 231 | No |

PDB mappings and human?mouse substitutions are carried forward from the validated repository evidence; they were not recalculated here.

## 3. Structural mapping

- Human EGFR 348 maps to **human 1IVO chain A SER 324**; the closest listed EGF contact context is **chain C TYR 44**.
- Human EGFR 442 maps to **human 1IVO chain A SER 418**; its direct geometric contact record is **chain C ARG 45**.

The contact distances and classifications describe the human 1IVO structure. The 348 site is classified as near-interface, not a direct EGF contact. The 442 contact has LOW source interaction confidence.

## 4. Final five candidates

| Candidate | Why selected | Interface evidence | Model availability | Key limitation |
|---|---|---|---|---|
| **442 S?G** | Human?mouse sequence substitution; starts the controlled 442 chemistry panel. | Direct geometric WT contact record to EGF C:45 ARG (4.439 ?; LOW interaction confidence). | Yes; human A:418 model. | Side-chain truncation/identity model, no rotamer repacking; geometry does not establish function. |
| **442 S?A** | Panel control removes hydroxyl while retaining a methyl group. | Same site-level direct contact evidence as 442; modeled distance is variant-specific. | Yes; human A:418 model. | Side-chain truncation/identity model, no rotamer repacking; WT contact typing is LOW confidence. |
| **442 S?T** | Panel control retains hydroxyl while adding methyl bulk. | Same site-level direct contact evidence as 442. | Yes; human A:418 model. | Static side-chain model; WT contact typing is LOW confidence. |
| **348 S?T** | Human?mouse substitution at a near-interface site. | NEAR_3D_INTERFACE; EGF C:44 TYR is the closest listed contact context (WT 5.640 ?); no direct candidate?EGF contact record. | Yes; human A:324 model. | Uses the first of three standard THR rotamers; PyMOL's optional bump-scoring/ranking step had a parser limitation. |
| **255 R?Q** | Human?mouse substitution selected as the non-interface/local-environment comparator; nearby environment is charged and polar. | Interface class NONE; WT minimum to EGF chain C is 18.732 ?. Interface-environment output is NOT_AVAILABLE, which is not negative evidence. | No mutant structure. | No mutant geometry or WT-vs-mutant distance comparison; salt-bridge context is heuristic only. |

## 5. 442 mechanistic panel

Comparing S?G, S?A, and S?T is a controlled side-chain chemistry/geometry perturbation. S?G removes the hydroxyl and almost all side-chain bulk; S?A removes the hydroxyl while retaining a methyl; S?T retains the hydroxyl and adds methyl bulk. The panel can test which chemical and geometric features track with experimental readouts; it does not predetermine an outcome.

## 6. Structural results

Minimum heavy-atom distances to EGF chain C from the validated human-structure models:

| Human-structure variant | WT distance | Mutant distance | Clash count |
|---|---:|---:|---:|
| 442 S?G | 4.439 ? | 5.626 ? | 0 |
| 442 S?A | 4.439 ? | 5.626 ? | 0 |
| 442 S?T | 4.439 ? | 4.286 ? | 0 |
| 348 S?T | 5.640 ? | 5.640 ? | 0 |

255 R?Q has a WT distance of 18.732 ? but no mutant model, so there is no WT-vs-mutant distance comparison. These distances are geometric observations, not binding-energy predictions.

## 7. Limitations

- The mutant structures are computational side-chain models on the human EGFR 1IVO structure; they are not experimentally observed structures.
- No molecular dynamics, energetic stability calculation, or binding-affinity calculation was performed.
- No biological causality or beneficial effect is claimed.
- The 348 S?T model uses a standard library rotamer. PyMOL's optional bump-scoring/ranking step had a parser limitation, so the first of three rotamers was used without that ranking step.
- The 442 S?A and S?G models are side-chain truncation/identity models without rotamer repacking.
- Heuristic salt-bridge counts are not confirmed interactions.
- NOT_AVAILABLE data is not negative evidence.
- Sequence-near status is not treated as 3D interface evidence.

## 8. Final conclusion

Human EGFR 442 is the strongest mechanistically testable structural site in this shortlist because it has a direct geometric interface contact record and a controlled three-variant side-chain panel. Human EGFR 348 provides a near-interface comparison, while human EGFR 255 R?Q provides a non-interface/local-environment comparator. The models help define testable experiments; they do not show whether any candidate changes EGFR function or whether any mutation is beneficial.
