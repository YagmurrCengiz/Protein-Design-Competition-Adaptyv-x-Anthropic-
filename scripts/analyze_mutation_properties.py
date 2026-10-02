import pandas as pd

input_path = "data/structures/egfr_structure_comparison.csv"

output_path = (
    "data/structures/"
    "egfr_mutation_properties.csv"
)

amino_acid_properties = {

    "A": {
        "charge": "neutral",
        "polarity": "nonpolar",
        "hydrophobicity": 1.8,
        "size": "small"
    },

    "R": {
        "charge": "positive",
        "polarity": "polar",
        "hydrophobicity": -4.5,
        "size": "large"
    },

    "N": {
        "charge": "neutral",
        "polarity": "polar",
        "hydrophobicity": -3.5,
        "size": "medium"
    },

    "D": {
        "charge": "negative",
        "polarity": "polar",
        "hydrophobicity": -3.5,
        "size": "medium"
    },

    "C": {
        "charge": "neutral",
        "polarity": "polar",
        "hydrophobicity": 2.5,
        "size": "medium"
    },

    "Q": {
        "charge": "neutral",
        "polarity": "polar",
        "hydrophobicity": -3.5,
        "size": "medium"
    },

    "E": {
        "charge": "negative",
        "polarity": "polar",
        "hydrophobicity": -3.5,
        "size": "medium"
    },

    "G": {
        "charge": "neutral",
        "polarity": "nonpolar",
        "hydrophobicity": -0.4,
        "size": "small"
    },

    "H": {
        "charge": "positive/neutral",
        "polarity": "polar",
        "hydrophobicity": -3.2,
        "size": "medium"
    },

    "I": {
        "charge": "neutral",
        "polarity": "nonpolar",
        "hydrophobicity": 4.5,
        "size": "medium"
    },

    "L": {
        "charge": "neutral",
        "polarity": "nonpolar",
        "hydrophobicity": 3.8,
        "size": "medium"
    },

    "K": {
        "charge": "positive",
        "polarity": "polar",
        "hydrophobicity": -3.9,
        "size": "large"
    },

    "M": {
        "charge": "neutral",
        "polarity": "nonpolar",
        "hydrophobicity": 1.9,
        "size": "large"
    },

    "F": {
        "charge": "neutral",
        "polarity": "nonpolar",
        "hydrophobicity": 2.8,
        "size": "large"
    },

    "P": {
        "charge": "neutral",
        "polarity": "nonpolar",
        "hydrophobicity": -1.6,
        "size": "medium"
    },

    "S": {
        "charge": "neutral",
        "polarity": "polar",
        "hydrophobicity": -0.8,
        "size": "small"
    },

    "T": {
        "charge": "neutral",
        "polarity": "polar",
        "hydrophobicity": -0.7,
        "size": "medium"
    },

    "W": {
        "charge": "neutral",
        "polarity": "nonpolar",
        "hydrophobicity": -0.9,
        "size": "large"
    },

    "Y": {
        "charge": "neutral",
        "polarity": "polar",
        "hydrophobicity": -1.3,
        "size": "large"
    },

    "V": {
        "charge": "neutral",
        "polarity": "nonpolar",
        "hydrophobicity": 4.2,
        "size": "medium"
    }
}

df = pd.read_csv(input_path)


for index, row in df.iterrows():

    human_aa = row["human"]
    mouse_aa = row["mouse"]

    human = amino_acid_properties[human_aa]
    mouse = amino_acid_properties[mouse_aa]

    df.loc[index, "human_charge"] = human["charge"]
    df.loc[index, "mouse_charge"] = mouse["charge"]

    df.loc[index, "charge_change"] = (
        f"{human['charge']} → {mouse['charge']}"
    )

    df.loc[index, "human_polarity"] = human["polarity"]
    df.loc[index, "mouse_polarity"] = mouse["polarity"]

    df.loc[index, "polarity_change"] = (
        f"{human['polarity']} → {mouse['polarity']}"
    )

    df.loc[index, "human_hydrophobicity"] = (
        human["hydrophobicity"]
    )

    df.loc[index, "mouse_hydrophobicity"] = (
        mouse["hydrophobicity"]
    )

    df.loc[index, "hydrophobicity_change"] = round(
        mouse["hydrophobicity"]
        - human["hydrophobicity"],
        2
    )

    df.loc[index, "size_change"] = (
        f"{human['size']} → {mouse['size']}"
    )

    if pd.notna(row["1IVO_SASA"]) and pd.notna(
        row["AlphaFold_SASA"]
    ):

        df.loc[index, "SASA_difference"] = round(
            row["AlphaFold_SASA"]
            - row["1IVO_SASA"],
            2
        )

    else:

        df.loc[index, "SASA_difference"] = None

    ivo_exposure = row["1IVO_exposure"]
    af_exposure = row["AlphaFold_exposure"]

    if (
        ivo_exposure != "NOT MAPPED"
        and pd.notna(af_exposure)
    ):

        if ivo_exposure == af_exposure:

            df.loc[index, "exposure_change"] = "NO CHANGE"

        else:

            df.loc[index, "exposure_change"] = (
                f"{ivo_exposure} → {af_exposure}"
            )

    else:

        df.loc[index, "exposure_change"] = "NOT COMPARABLE"


df.to_csv(
    output_path,
    index=False
)

print("=" * 80)
print("MUTATION PROPERTY ANALYSIS")
print("=" * 80)

print(
    "Input:",
    input_path
)

print(
    "Output:",
    output_path
)

print(
    "Rows:",
    len(df)
)

print("\nAnalysis completed.")