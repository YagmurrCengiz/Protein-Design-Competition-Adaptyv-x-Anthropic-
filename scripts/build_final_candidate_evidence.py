import pandas as pd
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

HIGH_PRIORITY_PATH = "data/structures/egfr_high_priority_analysis.csv"
DOMAIN_PATH = "data/structures/egfr_candidate_domains.csv"
EGF_PATH = "data/structures/egfr_egf_contacts.csv"
INTERACTION_PATH = "data/structures/egfr_candidate_interactions.csv"

OUTPUT_PATH = "data/structures/egfr_final_candidate_evidence.csv"

high_priority = pd.read_csv(ROOT / HIGH_PRIORITY_PATH)
domains = pd.read_csv(ROOT / DOMAIN_PATH)
egf = pd.read_csv(ROOT / EGF_PATH)
interactions = pd.read_csv(ROOT / INTERACTION_PATH)

print("=" * 90)
print("EGFR FINAL CANDIDATE EVIDENCE")
print("=" * 90)

print()
print("High-priority candidates:", len(high_priority))
print("Domain records:", len(domains))
print("EGF records:", len(egf))
print("Interaction records:", len(interactions))

if "exposure_changed" not in high_priority.columns:

    high_priority["exposure_changed"] = (
        high_priority["1IVO_exposure"].astype(str)
        != high_priority["AlphaFold_exposure"].astype(str)
    )

high_priority_columns = [
    "position",
    "human",
    "mouse",
    "mutation_type",
    "1IVO_PDB_number",
    "1IVO_PDB_residue",
    "1IVO_SASA",
    "1IVO_exposure",
    "AlphaFold_PDB_number",
    "AlphaFold_PDB_residue",
    "AlphaFold_SASA",
    "AlphaFold_exposure",
    "AlphaFold_pLDDT",
    "SASA_difference",
    "exposure_changed",
    "pLDDT_category",
    "priority_score",
    "impact_score",
    "evidence_category"
]

high_priority_subset = high_priority[
    [
        column
        for column in high_priority_columns
        if column in high_priority.columns
    ]
].copy()

domain_columns = [
    column
    for column in domains.columns
    if column != "position"
    and column not in high_priority_subset.columns
]

domain_subset = domains[
    ["position"] + domain_columns
].copy()

egf_columns = [
    column
    for column in egf.columns
    if column != "position"
    and column not in high_priority_subset.columns
    and column not in domain_subset.columns
]

egf_subset = egf[
    ["position"] + egf_columns
].copy()

interaction_columns = [
    column
    for column in interactions.columns
    if column != "position"
    and column not in high_priority_subset.columns
    and column not in domain_subset.columns
    and column not in egf_subset.columns
]

interaction_subset = interactions[
    ["position"] + interaction_columns
].copy()

final = high_priority_subset.merge(
    domain_subset,
    on="position",
    how="left"
)

final = final.merge(
    egf_subset,
    on="position",
    how="left"
)

final = final.merge(
    interaction_subset,
    on="position",
    how="left"
)

if "exposure_changed" not in final.columns:

    final["exposure_changed"] = (
        final["1IVO_exposure"].astype(str)
        != final["AlphaFold_exposure"].astype(str)
    )

print()
print("=" * 90)
print("FINAL CANDIDATES")
print("=" * 90)

for _, row in final.iterrows():

    print()
    print("-" * 90)

    print(
        f"{int(row['position'])}: "
        f"{row['human']}->{row['mouse']}"
    )

    print(
        "Mutation type:",
        row["mutation_type"]
    )

    domain_column = "domain" if "domain" in final.columns else "EGFR_domain"
    if domain_column in final.columns:
        print(
            "Domain:",
            row[domain_column]
        )

    for label, column in [
        ("pLDDT", "AlphaFold_pLDDT"),
        ("SASA difference", "SASA_difference"),
        ("Exposure changed", "exposure_changed"),
        ("Neighbor count", "neighbor_count"),
        ("Positive neighbors", "positive_neighbors"),
        ("Negative neighbors", "negative_neighbors"),
        ("Polar neighbors", "polar_neighbors"),
        ("Hydrophobic neighbors", "hydrophobic_neighbors"),
        ("Potential salt bridges", "potential_salt_bridges"),
        ("EGF contact", "EGF_contact"),
        ("EGF minimum distance", "EGF_min_distance"),
        ("Impact score", "impact_score"),
        ("Evidence category", "evidence_category"),
    ]:
        if column in final.columns:
            print(f"{label}:", row[column])

os.makedirs(
    (ROOT / OUTPUT_PATH).parent,
    exist_ok=True
)

final.to_csv(
    ROOT / OUTPUT_PATH,
    index=False
)

print()
print("=" * 90)
print("FINAL EVIDENCE TABLE SAVED")
print("=" * 90)

print(
    "Output:",
    ROOT / OUTPUT_PATH
)

print(
    "Candidates:",
    len(final)
)

print()
print(
    final.to_string(index = False)
)
