import pandas as pd

input_path = "data/structures/egfr_mutation_neighbors.csv"
output_path = "data/structures/egfr_mutation_environment.csv"

df = pd.read_csv(input_path)

charged = {
    "ASP",
    "GLU",
    "LYS",
    "ARG",
    "HIS"
}

polar = {
    "ASN",
    "GLN",
    "SER",
    "THR",
    "TYR"
}

hydrophobic = {
    "ALA",
    "VAL",
    "ILE",
    "LEU",
    "MET",
    "PHE",
    "TRP"
}

special = {
    "GLY",
    "PRO",
    "CYS"
}

def classify_residue(residue):
    if residue in charged:
        return "CHARGED"

    if residue in polar:
        return "POLAR"

    if residue in hydrophobic:
        return "HYDROPHOBIC"

    if residue in special:
        return "SPECIAL"

    return "OTHER"

df["neighbor_class"] = df[
    "neighbor_residue"
].apply(classify_residue)


environment = df.groupby(
    "position"
).agg(

    neighbor_count = (
        "neighbor_position",
        "count"
    ),

    closest_neighbor=(
        "distance",
        "min"
    ),

    charged_neighbors=(
        "neighbor_class",
        lambda x: (x == "CHARGED").sum()
    ),

    polar_neighbors=(
        "neighbor_class",
        lambda x: (x == "POLAR").sum()
    ),

    hydrophobic_neighbors=(
        "neighbor_class",
        lambda x: (x == "HYDROPHOBIC").sum()
    ),

    special_neighbors=(
        "neighbor_class",
        lambda x: (x == "SPECIAL").sum()
    )

).reset_index()

mutation_info = df[
    [
        "position",
        "human",
        "mouse",
        "mutation_type",
        "AlphaFold_pLDDT"
    ]
].drop_duplicates(
    subset="position"
)

environment = environment.merge(
    mutation_info,
    on = "position",
    how = "left"
)

environment["chemical_contacts"] = (
    environment["charged_neighbors"]
    +
    environment["polar_neighbors"]
)

environment = environment.sort_values(
    by=[
        "chemical_contacts",
        "neighbor_count"
    ],
    ascending=[
        False,
        False
    ]
)

environment.to_csv(
    output_path,
    index=False
)

print("=" * 80)
print("EGFR MUTATION CHEMICAL ENVIRONMENT")
print("=" * 80)

print()

print(
    environment.head(25).to_string(
        index=False
    )
)

print()

print(
    "Total mutation environments:",
    len(environment)
)

print()

print(
    "Saved:",
    output_path
)