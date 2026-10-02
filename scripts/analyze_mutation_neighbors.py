from Bio.PDB import PDBParser
import pandas as pd
import math


pdb_path = "data/structures/alphafold_egfr.pdb"

input_path = "data/structures/egfr_prioritized_mutations.csv"

output_path = "data/structures/egfr_mutation_neighbors.csv"


parser = PDBParser(QUIET=True)

structure = parser.get_structure(
    "EGFR_AlphaFold",
    pdb_path
)

model = structure[0]

chain = next(iter(model))


df = pd.read_csv(input_path)


candidate_df = df[
    (df["AlphaFold_pLDDT"] >= 70)
    &
    (df["AlphaFold_exposure"] != "NOT MAPPED")
].copy()


residue_lookup = {}

for residue in chain:

    resname = residue.get_resname()

    residue_number = residue.id[1]

    residue_lookup[residue_number] = residue


def distance(atom1, atom2):

    return math.sqrt(
        sum(
            (a - b) ** 2
            for a, b in zip(
                atom1.coord,
                atom2.coord
            )
        )
    )


results = []


for _, row in candidate_df.iterrows():

    position = int(row["position"])

    target_residue = residue_lookup.get(position)

    if target_residue is None:
        continue


    target_atoms = list(
        target_residue.get_atoms()
    )


    neighbors = []


    for residue_number, residue in residue_lookup.items():

        if residue_number == position:
            continue


        residue_atoms = list(
            residue.get_atoms()
        )


        minimum_distance = min(
            distance(
                atom1,
                atom2
            )
            for atom1 in target_atoms
            for atom2 in residue_atoms
        )


        if minimum_distance <= 5.0:

            neighbors.append({
                "position": position,
                "neighbor_position": residue_number,
                "neighbor_residue": residue.get_resname(),
                "distance": round(
                    minimum_distance,
                    2
                )
            })


    if neighbors:

        for neighbor in neighbors:

            neighbor["human"] = row["human"]
            neighbor["mouse"] = row["mouse"]
            neighbor["mutation_type"] = row["mutation_type"]
            neighbor["AlphaFold_pLDDT"] = row["AlphaFold_pLDDT"]

            results.append(neighbor)


result_df = pd.DataFrame(results)


result_df = result_df.sort_values(
    by=[
        "position",
        "distance"
    ]
)


result_df.to_csv(
    output_path,
    index=False
)


print("=" * 80)
print("EGFR MUTATION NEIGHBOR ANALYSIS")
print("=" * 80)

print()

print(
    "Candidates analyzed:",
    candidate_df.shape[0]
)

print()

print(
    "Neighbor interactions:",
    result_df.shape[0]
)

print()

print(
    result_df.head(30).to_string(
        index=False
    )
)

print()

print(
    "Saved:",
    output_path
)