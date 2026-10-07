import pandas as pd
import os

INPUT_PATH = "data/structures/egfr_egf_interface.csv"
OUTPUT_PATH = "data/structures/egfr_egf_interface_interactions.csv"

CONTACT_DISTANCE = 5.0
SALT_BRIDGE_DISTANCE = 4.0

charged_positive = {
    "LYS",
    "ARG",
    "HIS"
}

charged_negative = {
    "ASP",
    "GLU"
}

polar_residues = {
    "SER",
    "THR",
    "ASN",
    "GLN",
    "TYR",
    "HIS",
    "CYS"
}

hydrophobic_residues = {
    "ALA",
    "VAL",
    "ILE",
    "LEU",
    "MET",
    "PHE",
    "TRP",
    "PRO"
}

aromatic_residues = {
    "PHE",
    "TYR",
    "TRP",
    "HIS"
}


def classify_residue(residue):

    if residue in charged_positive:
        return "POSITIVE"

    if residue in charged_negative:
        return "NEGATIVE"

    if residue in polar_residues:
        return "POLAR"

    if residue in hydrophobic_residues:
        return "HYDROPHOBIC"

    return "OTHER"


def classify_interaction(
    egfr_residue,
    egf_residue,
    distance
):

    egfr_class = classify_residue(
        egfr_residue
    )

    egf_class = classify_residue(
        egf_residue
    )

    if (
        distance <= SALT_BRIDGE_DISTANCE
        and (
            (
                egfr_class == "POSITIVE"
                and egf_class == "NEGATIVE"
            )
            or
            (
                egfr_class == "NEGATIVE"
                and egf_class == "POSITIVE"
            )
        )
    ):

        return (
            "POSSIBLE_SALT_BRIDGE",
            "HIGH"
        )

    if (
        egfr_residue in aromatic_residues
        and egf_residue in aromatic_residues
        and distance <= 5.0
    ):

        return (
            "AROMATIC_CONTACT",
            "MEDIUM"
        )

    if (
        egfr_class == "HYDROPHOBIC"
        and egf_class == "HYDROPHOBIC"
        and distance <= 4.0
    ):

        return (
            "HYDROPHOBIC_CONTACT",
            "MEDIUM"
        )

    if (
        (
            egfr_class in {
                "POSITIVE",
                "NEGATIVE"
            }
            and egf_class in {
                "POSITIVE",
                "NEGATIVE"
            }
        )
        and distance <= 5.0
    ):

        return (
            "CHARGED_CONTACT",
            "MEDIUM"
        )

    if (
        (
            egfr_class == "POLAR"
            or egf_class == "POLAR"
        )
        and distance <= 4.0
    ):

        return (
            "POLAR_CONTACT",
            "MEDIUM"
        )

    if (
        (
            egfr_class in {
                "POLAR",
                "POSITIVE",
                "NEGATIVE"
            }
            or
            egf_class in {
                "POLAR",
                "POSITIVE",
                "NEGATIVE"
            }
        )
        and distance <= 3.5
    ):

        return (
            "POLAR_CONTACT",
            "LOW"
        )

    return (
        "OTHER",
        "LOW"
    )


df = pd.read_csv(
    INPUT_PATH
)

print("=" * 90)
print("EGFR-EGF INTERACTION CLASSIFICATION")
print("=" * 90)

print()
print(
    "Interface contacts:",
    len(df)
)

results = []

for _, row in df.iterrows():

    egfr_residue = row["EGFR_residue"]
    egf_residue = row["EGF_residue"]
    distance = float(
        row["minimum_distance"]
    )

    egfr_class = classify_residue(
        egfr_residue
    )

    egf_class = classify_residue(
        egf_residue
    )

    interaction_type, confidence = (
        classify_interaction(
            egfr_residue,
            egf_residue,
            distance
        )
    )

    result = row.to_dict()

    result["EGFR_class"] = egfr_class
    result["EGF_class"] = egf_class
    result["interaction_type"] = interaction_type
    result["interaction_confidence"] = confidence

    results.append(result)

result_df = pd.DataFrame(
    results
)

result_df = result_df.sort_values(
    by="minimum_distance"
)

print()
print("=" * 90)
print("INTERACTION SUMMARY")
print("=" * 90)

print()

print(
    result_df[
        [
            "EGFR_position",
            "EGFR_residue",
            "EGF_chain",
            "EGF_position",
            "EGF_residue",
            "minimum_distance",
            "EGFR_class",
            "EGF_class",
            "interaction_type",
            "interaction_confidence"
        ]
    ].to_string(
        index=False
    )
)

print()
print("=" * 90)
print("INTERACTION COUNTS")
print("=" * 90)

print()

print(
    result_df[
        "interaction_type"
    ].value_counts().to_string()
)

print()
print("=" * 90)
print("HIGH-CONFIDENCE INTERACTIONS")
print("=" * 90)

high_confidence = result_df[
    result_df[
        "interaction_confidence"
    ] == "HIGH"
]

if len(high_confidence) > 0:

    print()

    print(
        high_confidence[
            [
                "EGFR_position",
                "EGFR_residue",
                "EGF_chain",
                "EGF_position",
                "EGF_residue",
                "minimum_distance",
                "interaction_type"
            ]
        ].to_string(
            index=False
        )
    )

else:

    print()
    print(
        "No high-confidence interactions detected."
    )

os.makedirs(
    os.path.dirname(OUTPUT_PATH),
    exist_ok=True
)

result_df.to_csv(
    OUTPUT_PATH,
    index=False
)

print()
print("=" * 90)
print("INTERACTION ANALYSIS SAVED")
print("=" * 90)

print(
    "Output:",
    OUTPUT_PATH
)

print(
    "Interactions:",
    len(result_df)
)