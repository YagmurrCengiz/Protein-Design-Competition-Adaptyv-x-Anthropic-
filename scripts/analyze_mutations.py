import pandas as pd

input_path = "data/structures/egfr_structure_comparison.csv"
output_path = "data/structures/egfr_mutation_analysis.csv"

df = pd.read_csv(input_path)

def classify_mutation(human, mouse):
    charged = {
        "D", "E", "K", "R"
    }
    polar = {
        "S", "T", "N", "Q"
    }
    hydrophobic = {
        "A", "V", "I", "L", "M", "F", "W", "Y"
    }

    if human in charged and mouse not in charged:
        return "LOSS_OF_CHARGE"

    if human not in charged and mouse in charged:
        return "GAIN_OF_CHARGE"

    if human in polar and mouse in hydrophobic:
        return "POLAR_TO_HYDROPHOBIC"

    if human in hydrophobic and mouse in polar:
        return "HYDROPHOBIC_TO_POLAR"

    if human in polar and mouse in charged:
        return "POLAR_TO_CHARGED"

    if human in charged and mouse in polar:
        return "CHARGED_TO_POLAR"

    if human in hydrophobic and mouse in charged:
        return "HYDROPHOBIC_TO_CHARGED"

    if human in charged and mouse in hydrophobic:
        return "CHARGED_TO_HYDROPHOBIC"

    return "OTHER"

df["mutation_type"] = df.apply(
    lambda row: classify_mutation(
        row["human"],
        row["mouse"]
    ),
    axis=1
)


df["SASA_difference"] = (
    df["AlphaFold_SASA"] -
    df["1IVO_SASA"]
)


df["exposure_changed"] = (
    df["1IVO_exposure"] !=
    df["AlphaFold_exposure"]
)


df["pLDDT_category"] = pd.cut(
    df["AlphaFold_pLDDT"],
    bins=[0, 50, 70, 90, 100],
    labels=[
        "Very Low",
        "Low",
        "Confident",
        "Very High"
    ],
    include_lowest=True
)


df.to_csv(
    output_path,
    index=False
)


print("=" * 80)
print("EGFR MUTATION ANALYSIS")
print("=" * 80)

print()

print(
    "Total mutations:",
    len(df)
)

print()

print("Mutation types:")
print(
    df["mutation_type"].value_counts()
)

print()

print("Exposure changes:")
print(
    df["exposure_changed"].value_counts()
)

print()

print("pLDDT categories:")
print(
    df["pLDDT_category"].value_counts()
)

print()

print(
    "Saved:",
    output_path
)