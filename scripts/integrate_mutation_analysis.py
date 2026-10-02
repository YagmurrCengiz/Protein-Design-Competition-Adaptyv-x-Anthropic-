import pandas as pd
import os


# ============================================================
# FILE PATHS
# ============================================================

priority_path = (
    "data/structures/"
    "egfr_prioritized_mutations.csv"
)

environment_path = (
    "data/structures/"
    "egfr_mutation_environment.csv"
)

output_path = (
    "data/structures/"
    "egfr_integrated_analysis.csv"
)


# ============================================================
# LOAD DATA
# ============================================================

priority_df = pd.read_csv(priority_path)

environment_df = pd.read_csv(environment_path)

print("=" * 90)
print("EGFR MUTATION INTEGRATED ANALYSIS")
print("=" * 90)

print()

print(
    "Priority mutations:",
    len(priority_df)
)

print(
    "Chemical environments:",
    len(environment_df)
)

required_priority_columns = [
    "position",
    "human",
    "mouse",
    "mutation_type",
    "1IVO_exposure",
    "AlphaFold_exposure",
    "AlphaFold_pLDDT"
]

required_environment_columns = [
    "position",
    "neighbor_count",
    "closest_neighbor",
    "charged_neighbors",
    "polar_neighbors",
    "hydrophobic_neighbors",
    "special_neighbors",
    "chemical_contacts"
]

for column in required_priority_columns:

    if column not in priority_df.columns:

        raise ValueError(
            f"Missing column in priority file: {column}"
        )

for column in required_environment_columns:

    if column not in environment_df.columns:

        raise ValueError(
            f"Missing column in environment file: {column}"
        )

environment_features = environment_df[
    required_environment_columns
].copy()

priority_df["position"] = pd.to_numeric(
    priority_df["position"],
    errors="coerce"
)

environment_features["position"] = pd.to_numeric(
    environment_features["position"],
    errors="coerce"
)

integrated = priority_df.merge(
    environment_features,
    on="position",
    how="left"
)

print()

print(
    "Integrated rows:",
    len(integrated)
)

environment_numeric_columns = [
    "neighbor_count",
    "closest_neighbor",
    "charged_neighbors",
    "polar_neighbors",
    "hydrophobic_neighbors",
    "special_neighbors",
    "chemical_contacts"
]

for column in environment_numeric_columns:

    integrated[column] = pd.to_numeric(
        integrated[column],
        errors="coerce"
    )

    integrated[column] = integrated[column].fillna(0)

integrated["experimental_structure"] = (
    integrated["1IVO_exposure"]
    != "NOT MAPPED"
)

integrated["exposure_changed"] = False

mapped_mask = (
    integrated["1IVO_exposure"]
    != "NOT MAPPED"
)

integrated.loc[
    mapped_mask,
    "exposure_changed"
] = (
    integrated.loc[
        mapped_mask,
        "1IVO_exposure"
    ]
    !=
    integrated.loc[
        mapped_mask,
        "AlphaFold_exposure"
    ]
)

integrated["chemical_environment_score"] = (

    integrated["charged_neighbors"]

    +

    integrated["polar_neighbors"]

    +

    integrated["chemical_contacts"]

)

integrated["structural_change_score"] = 0

integrated.loc[
    integrated["exposure_changed"],
    "structural_change_score"
] += 2

integrated.loc[
    integrated["experimental_structure"],
    "structural_change_score"
] += 1

integrated.loc[
    integrated["chemical_contacts"] >= 7,
    "structural_change_score"
] += 2

integrated.loc[
    integrated["charged_neighbors"] >= 4,
    "structural_change_score"
] += 2

integrated.loc[
    integrated["polar_neighbors"] >= 5,
    "structural_change_score"
] += 1

mutation_type_scores = {

    "LOSS_OF_CHARGE": 3,

    "GAIN_OF_CHARGE": 3,

    "POLAR_TO_HYDROPHOBIC": 3,

    "HYDROPHOBIC_TO_POLAR": 2,

    "OTHER": 0
}

integrated["mutation_type_score"] = (
    integrated["mutation_type"]
    .map(mutation_type_scores)
    .fillna(0)
)

integrated["confidence_score"] = 0

integrated.loc[
    integrated["AlphaFold_pLDDT"] >= 90,
    "confidence_score"
] = 2

integrated.loc[
    (
        integrated["AlphaFold_pLDDT"] >= 70
    )
    &
    (
        integrated["AlphaFold_pLDDT"] < 90
    ),
    "confidence_score"
] = 1

integrated["impact_score"] = (

    integrated["structural_change_score"]

    +

    integrated["mutation_type_score"]

    +

    integrated["confidence_score"]

)

def classify_evidence(row):

    score = row["impact_score"]

    if score >= 9:

        return "HIGH_STRUCTURAL_INTEREST"

    elif score >= 6:

        return "MODERATE_STRUCTURAL_INTEREST"

    else:

        return "LOW_STRUCTURAL_INTEREST"

integrated["evidence_category"] = (
    integrated.apply(
        classify_evidence,
        axis=1
    )
)

integrated = integrated.sort_values(

    by=[
        "impact_score",
        "chemical_environment_score",
        "AlphaFold_pLDDT"
    ],

    ascending=False
)

display_columns = [

    "position",

    "human",

    "mouse",

    "mutation_type",

    "1IVO_exposure",

    "AlphaFold_exposure",

    "AlphaFold_pLDDT",

    "neighbor_count",

    "charged_neighbors",

    "polar_neighbors",

    "hydrophobic_neighbors",

    "chemical_contacts",

    "exposure_changed",

    "experimental_structure",

    "mutation_type_score",

    "structural_change_score",

    "confidence_score",

    "chemical_environment_score",

    "impact_score",

    "evidence_category"
]

print()

print("=" * 120)
print("TOP STRUCTURAL CANDIDATES")
print("=" * 120)

print()

print(
    integrated[
        display_columns
    ]
    .head(20)
    .to_string(index=False)
)

print()

print("=" * 90)
print("EVIDENCE CATEGORY SUMMARY")
print("=" * 90)

print()

print(
    integrated[
        "evidence_category"
    ].value_counts()
)

integrated.to_csv(
    output_path,
    index=False
)

print()

print("=" * 90)
print("INTEGRATED ANALYSIS SAVED")
print("=" * 90)

print()

print(
    "Output:",
    output_path
)

print(
    "Rows:",
    len(integrated)
)