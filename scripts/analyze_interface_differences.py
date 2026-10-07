import pandas as pd
import os

MUTATION_PATH = "data/structures/egfr_mutation_analysis.csv"
INTERFACE_PATH = "data/structures/egfr_egf_interface.csv"
OUTPUT_PATH = "data/structures/egfr_interface_differences.csv"

INTERFACE_DISTANCE = 5.0

mutations = pd.read_csv(MUTATION_PATH)
interface = pd.read_csv(INTERFACE_PATH)

print("=" * 90)
print("EGFR INTERFACE DIFFERENCE ANALYSIS")
print("=" * 90)

print()
print("Mutation records:", len(mutations))
print("Interface records:", len(interface))

mutation_positions = set(
    mutations["position"].astype(int)
)

interface_positions = set(
    interface["EGFR_position"].astype(int)
)

direct_matches = mutation_positions.intersection(
    interface_positions
)

print()
print("=" * 90)
print("DIRECT INTERFACE MATCHES")
print("=" * 90)

print()
print(
    "Human-mouse differences on EGFR-EGF interface:",
    len(direct_matches)
)

print(
    "Positions:",
    sorted(direct_matches)
)

interface_residues = (
    interface[
        [
            "EGFR_position",
            "EGFR_residue",
            "EGF_chain",
            "EGF_position",
            "EGF_residue",
            "minimum_distance",
            "EGFR_atom",
            "EGF_atom"
        ]
    ]
    .drop_duplicates()
)

results = []

for _, mutation in mutations.iterrows():

    position = int(mutation["position"])

    mutation_interface = interface_residues[
        interface_residues["EGFR_position"] == position
    ]

    if len(mutation_interface) > 0:

        for _, contact in mutation_interface.iterrows():

            results.append({

                "position": position,

                "human": mutation["human"],

                "mouse": mutation["mouse"],

                "mutation_type": mutation.get(
                    "mutation_type",
                    None
                ),

                "EGFR_residue": contact[
                    "EGFR_residue"
                ],

                "EGF_chain": contact[
                    "EGF_chain"
                ],

                "EGF_position": contact[
                    "EGF_position"
                ],

                "EGF_residue": contact[
                    "EGF_residue"
                ],

                "interface_distance": contact[
                    "minimum_distance"
                ],

                "EGFR_atom": contact[
                    "EGFR_atom"
                ],

                "EGF_atom": contact[
                    "EGF_atom"
                ],

                "direct_interface_match": True

            })

interface_positions_list = sorted(
    interface_positions
)

for _, mutation in mutations.iterrows():

    position = int(mutation["position"])

    if position in interface_positions:
        continue

    closest_interface_position = None
    closest_distance = float("inf")

    for interface_position in interface_positions_list:

        distance = abs(
            position - interface_position
        )

        if distance < closest_distance:

            closest_distance = distance
            closest_interface_position = interface_position

    if closest_distance <= INTERFACE_DISTANCE:

        matching_interface = interface_residues[
            interface_residues[
                "EGFR_position"
            ] == closest_interface_position
        ]

        for _, contact in matching_interface.iterrows():

            results.append({

                "position": position,

                "human": mutation["human"],

                "mouse": mutation["mouse"],

                "mutation_type": mutation.get(
                    "mutation_type",
                    None
                ),

                "EGFR_residue": None,

                "EGF_chain": contact[
                    "EGF_chain"
                ],

                "EGF_position": contact[
                    "EGF_position"
                ],

                "EGF_residue": contact[
                    "EGF_residue"
                ],

                "interface_distance": contact[
                    "minimum_distance"
                ],

                "EGFR_atom": None,

                "EGF_atom": contact[
                    "EGF_atom"
                ],

                "direct_interface_match": False

            })

results_df = pd.DataFrame(results)

if len(results_df) > 0:

    results_df = results_df.sort_values(
        [
            "direct_interface_match",
            "interface_distance"
        ],
        ascending=[
            False,
            True
        ]
    )

os.makedirs(
    os.path.dirname(OUTPUT_PATH),
    exist_ok=True
)

results_df.to_csv(
    OUTPUT_PATH,
    index=False
)

print()
print("=" * 90)
print("INTERFACE DIFFERENCE SUMMARY")
print("=" * 90)

print()

print(
    "Direct interface mutations:",
    int(
        results_df["direct_interface_match"].sum()
    )
    if len(results_df) > 0
    else 0
)

print(
    "Interface-near mutations:",
    int(
        (~results_df["direct_interface_match"]).sum()
    )
    if len(results_df) > 0
    else 0
)

print()

if len(results_df) > 0:

    print(
        results_df.to_string(
            index=False
        )
    )

else:

    print(
        "No human-mouse differences were found "
        "on or near the EGFR-EGF interface."
    )

print()
print("=" * 90)
print("INTERFACE DIFFERENCE ANALYSIS SAVED")
print("=" * 90)

print(
    "Output:",
    OUTPUT_PATH
)