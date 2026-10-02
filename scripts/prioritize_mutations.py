import pandas as pd

input_path = "data/structures/egfr_mutation_analysis.csv"
output_path = "data/structures/egfr_prioritized_mutations.csv"

df = pd.read_csv(input_path)

def calculate_priority(row):

    score = 0

    if row["mutation_type"] != "OTHER":
        score += 2

    if row["mutation_type"] in [
        "LOSS_OF_CHARGE",
        "GAIN_OF_CHARGE"
    ]:
        score += 1

    if row["mutation_type"] in [
        "POLAR_TO_HYDROPHOBIC",
        "HYDROPHOBIC_TO_POLAR"
    ]:
        score += 1

    if row["exposure_changed"]:
        score += 2

    if row["AlphaFold_exposure"] == "EXPOSED":
        score += 1

    if row["AlphaFold_exposure"] == "PARTIALLY EXPOSED":
        score += 1

    if row["AlphaFold_pLDDT"] >= 90:
        score += 2

    elif row["AlphaFold_pLDDT"] >= 70:
        score += 1

    return score


df["priority_score"] = df.apply(
    calculate_priority,
    axis=1
)


df = df.sort_values(
    by=[
        "priority_score",
        "AlphaFold_pLDDT"
    ],
    ascending=[
        False,
        False
    ]
)

df.to_csv(
    output_path,
    index = False
)

print("=" * 80)
print("EGFR MUTATION PRIORITIZATION")
print("=" * 80)

print()

print(
    df[
        [
            "position",
            "human",
            "mouse",
            "mutation_type",
            "1IVO_exposure",
            "AlphaFold_exposure",
            "AlphaFold_pLDDT",
            "priority_score"
        ]
    ].head(20).to_string(index = False)
)

print()

print(
    "Saved:",
    output_path
)