import pandas as pd

integrated_path = (
    "data/structures/"
    "egfr_integrated_analysis.csv"
)

neighbors_path = (
    "data/structures/"
    "egfr_mutation_neighbors.csv"
)

output_path = (
    "data/structures/"
    "egfr_high_priority_analysis.csv"
)

integrated = pd.read_csv(
    integrated_path
)

neighbors = pd.read_csv(
    neighbors_path
)

high = integrated[
    integrated["evidence_category"]
    ==
    "HIGH_STRUCTURAL_INTEREST"
].copy()

positions = high["position"].tolist()

print("=" * 90)
print("EGFR HIGH-PRIORITY CANDIDATE ANALYSIS")
print("=" * 90)

print()

print(
    "Candidates:",
    len(positions)
)

print(
    "Positions:",
    positions
)

results = []

for _, mutation in high.iterrows():

    position = int(
        mutation["position"]
    )

    human = mutation["human"]

    mouse = mutation["mouse"]

    mutation_type = mutation[
        "mutation_type"
    ]

    print()
    print("=" * 90)
    print(
        f"{position}: "
        f"{human}->{mouse}"
    )
    print("=" * 90)

    candidate_neighbors = neighbors[
        neighbors["position"]
        ==
        position
    ].copy()

    candidate_neighbors = (
        candidate_neighbors
        .sort_values(
            "distance"
        )
    )

    print()

    print("Mutation type:")
    print(
        mutation_type
    )

    print()

    print("1IVO:")
    print(
        mutation["1IVO_PDB_residue"],
        mutation["1IVO_PDB_number"]
    )

    print(
        "SASA:",
        mutation["1IVO_SASA"]
    )

    print(
        "Exposure:",
        mutation["1IVO_exposure"]
    )

    print()

    print("AlphaFold:")

    print(
        mutation["AlphaFold_PDB_residue"],
        mutation["AlphaFold_PDB_number"]
    )

    print(
        "SASA:",
        mutation["AlphaFold_SASA"]
    )

    print(
        "Exposure:",
        mutation["AlphaFold_exposure"]
    )

    print(
        "pLDDT:",
        mutation["AlphaFold_pLDDT"]
    )

    print()

    print("Chemical environment:")

    print(
        "Neighbors:",
        mutation["neighbor_count"]
    )

    print(
        "Charged:",
        mutation["charged_neighbors"]
    )

    print(
        "Polar:",
        mutation["polar_neighbors"]
    )

    print(
        "Hydrophobic:",
        mutation["hydrophobic_neighbors"]
    )

    print(
        "Special:",
        mutation["special_neighbors"]
    )

    print(
        "Chemical contacts:",
        mutation["chemical_contacts"]
    )

    print()

    print("Closest residues:")

    print(
        candidate_neighbors[
            [
                "neighbor_position",
                "neighbor_residue",
                "distance"
            ]
        ]
        .head(10)
        .to_string(index=False)
    )

    residue_counts = (
        candidate_neighbors[
            "neighbor_residue"
        ]
        .value_counts()
    )

    residue_summary = ", ".join(
        [
            f"{residue}:{count}"
            for residue, count
            in residue_counts.items()
        ]
    )

    if len(candidate_neighbors) > 0:

        mean_distance = round(
            candidate_neighbors[
                "distance"
            ].mean(),
            2
        )

        minimum_distance = round(
            candidate_neighbors[
                "distance"
            ].min(),
            2
        )

    else:

        mean_distance = None

        minimum_distance = None

    results.append({

        "position": position,

        "human": human,

        "mouse": mouse,

        "mutation_type":
            mutation_type,

        "1IVO_exposure":
            mutation["1IVO_exposure"],

        "AlphaFold_exposure":
            mutation["AlphaFold_exposure"],

        "1IVO_SASA":
            mutation["1IVO_SASA"],

        "AlphaFold_SASA":
            mutation["AlphaFold_SASA"],

        "SASA_difference":
            mutation["SASA_difference"],

        "AlphaFold_pLDDT":
            mutation["AlphaFold_pLDDT"],

        "neighbor_count":
            mutation["neighbor_count"],

        "closest_neighbor":
            mutation["closest_neighbor"],

        "charged_neighbors":
            mutation["charged_neighbors"],

        "polar_neighbors":
            mutation["polar_neighbors"],

        "hydrophobic_neighbors":
            mutation["hydrophobic_neighbors"],

        "special_neighbors":
            mutation["special_neighbors"],

        "chemical_contacts":
            mutation["chemical_contacts"],

        "mean_neighbor_distance":
            mean_distance,

        "minimum_neighbor_distance":
            minimum_distance,

        "neighbor_residue_summary":
            residue_summary,

        "impact_score":
            mutation["impact_score"]

    })

result_df = pd.DataFrame(
    results
)

result_df.to_csv(
    output_path,
    index=False
)

print()
print("=" * 90)
print("HIGH-PRIORITY ANALYSIS SAVED")
print("=" * 90)

print()

print(
    "Output:",
    output_path
)

print(
    "Candidates:",
    len(result_df)
)

print()

print(
    result_df.to_string(
        index=False
    )
)