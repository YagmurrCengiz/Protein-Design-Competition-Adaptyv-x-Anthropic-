# EGFR competition submission v2: diversity and sequence report

## Scope and limits

Twelve de novo sequence hypotheses across several length and architecture families were internally drafted. Candidates with obvious local redundancy were removed; the exported panel contains seven new candidates and the unchanged EGFR_pHmini_08 anchor. Local comparisons cannot establish global sequence novelty. Only the competition platform determines novelty; a score of 3/4 or higher is not guaranteed.

## Biological and pH rationale

The objective is a single-chain de novo conditional-binder hypothesis for human and mouse EGFR. The existing Domain III alignment is 154/170 identical (90.6%). The designs are intended to recognize the conserved EGF-facing patch (D379, S380, F381, T382, L406, Q408, Q432, H433, Q435, F436, A439 and V441), avoiding sole dependence on divergent S348 or S442. Histidines occur at varied positions as candidate titratable interface residues. Partial protonation at pH 6.5 could alter electrostatic or hydrogen-bond contributions near acidic D379; lower protonation at pH 7.4 could weaken them. Histidine alone does not guarantee selectivity. Folding, binding, cross-species recognition and pH response are unpredicted and unvalidated. No therapeutic antibody sequence was intentionally copied.

## Exported candidates

| Candidate | Length | Architecture/family | Closest old-panel identity |
|---|---:|---|---:|
| EGFR_pHmini_v2_01 | 64 | compact three-helix bundle | 36.2% (EGFR_pHmini_07) |
| EGFR_pHmini_v2_02 | 45 | two-helix hairpin | 36.7% (EGFR_pHmini_07) |
| EGFR_pHmini_v2_03 | 63 | beta-alpha mini-domain | 21.7% (EGFR_pHmini_08) |
| EGFR_pHmini_v2_04 | 61 | beta-sheet with helical cap | 23.0% (EGFR_pHmini_08) |
| EGFR_pHmini_v2_05 | 68 | compact four-helix bundle | 31.4% (EGFR_pHmini_03) |
| EGFR_pHmini_v2_06 | 93 | alpha/beta domain with terminal helix | 26.9% (EGFR_pHmini_08) |
| EGFR_pHmini_v2_07 | 65 | beta-hairpin-rich mini-domain with acidic/polar turns | 16.9% (EGFR_pHmini_07) |
| EGFR_pHmini_08 | 59 | unchanged original three-helix anchor | 100.0% (EGFR_pHmini_08) |

## Pairwise identity among exported candidates

Global Needleman–Wunsch identity: matches divided by aligned columns, including gaps. Closest pair: **EGFR_pHmini_v2_04 / EGFR_pHmini_v2_06 = 65.6%**. This measures sequence identity only.

| Candidate A | Candidate B | Identity |
|---|---|---:|
| EGFR_pHmini_v2_04 | EGFR_pHmini_v2_06 | 65.6% |
| EGFR_pHmini_v2_03 | EGFR_pHmini_v2_06 | 63.4% |
| EGFR_pHmini_v2_03 | EGFR_pHmini_v2_04 | 61.8% |
| EGFR_pHmini_v2_01 | EGFR_pHmini_v2_05 | 45.5% |
| EGFR_pHmini_v2_02 | EGFR_pHmini_v2_05 | 44.1% |
| EGFR_pHmini_v2_01 | EGFR_pHmini_v2_02 | 43.1% |
| EGFR_pHmini_v2_02 | EGFR_pHmini_08 | 35.0% |
| EGFR_pHmini_v2_05 | EGFR_pHmini_v2_06 | 33.3% |
| EGFR_pHmini_v2_01 | EGFR_pHmini_08 | 33.3% |
| EGFR_pHmini_v2_01 | EGFR_pHmini_v2_07 | 31.9% |
| EGFR_pHmini_v2_03 | EGFR_pHmini_v2_05 | 31.1% |
| EGFR_pHmini_v2_04 | EGFR_pHmini_v2_05 | 30.9% |
| EGFR_pHmini_v2_01 | EGFR_pHmini_v2_04 | 30.9% |
| EGFR_pHmini_v2_01 | EGFR_pHmini_v2_06 | 30.5% |
| EGFR_pHmini_v2_06 | EGFR_pHmini_v2_07 | 30.1% |
| EGFR_pHmini_v2_05 | EGFR_pHmini_08 | 30.0% |
| EGFR_pHmini_v2_03 | EGFR_pHmini_v2_07 | 29.9% |
| EGFR_pHmini_v2_02 | EGFR_pHmini_v2_04 | 28.1% |
| EGFR_pHmini_v2_01 | EGFR_pHmini_v2_03 | 27.9% |
| EGFR_pHmini_v2_06 | EGFR_pHmini_08 | 26.9% |
| EGFR_pHmini_v2_02 | EGFR_pHmini_v2_06 | 26.6% |
| EGFR_pHmini_v2_04 | EGFR_pHmini_v2_07 | 26.1% |
| EGFR_pHmini_v2_02 | EGFR_pHmini_v2_03 | 25.0% |
| EGFR_pHmini_v2_05 | EGFR_pHmini_v2_07 | 24.7% |
| EGFR_pHmini_v2_02 | EGFR_pHmini_v2_07 | 23.1% |
| EGFR_pHmini_v2_04 | EGFR_pHmini_08 | 23.0% |
| EGFR_pHmini_v2_03 | EGFR_pHmini_08 | 21.7% |
| EGFR_pHmini_v2_07 | EGFR_pHmini_08 | 15.4% |

## Comparison against old panel

| New candidate | Closest old candidate | Identity |
|---|---|---:|
| EGFR_pHmini_v2_01 | EGFR_pHmini_07 | 36.2% |
| EGFR_pHmini_v2_02 | EGFR_pHmini_07 | 36.7% |
| EGFR_pHmini_v2_03 | EGFR_pHmini_08 | 21.7% |
| EGFR_pHmini_v2_04 | EGFR_pHmini_08 | 23.0% |
| EGFR_pHmini_v2_05 | EGFR_pHmini_03 | 31.4% |
| EGFR_pHmini_v2_06 | EGFR_pHmini_08 | 26.9% |
| EGFR_pHmini_v2_07 | EGFR_pHmini_07 | 16.9% |
| EGFR_pHmini_08 | EGFR_pHmini_08 | 100.0% |

Alignment scoring: match +2, mismatch -1, gap -2. Substitution-matrix similarity is omitted because the competition does not provide a local novelty metric. Local identity does not establish global novelty.
