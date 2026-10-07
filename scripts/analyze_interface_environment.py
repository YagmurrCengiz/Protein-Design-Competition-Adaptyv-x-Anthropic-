from Bio.PDB import PDBParser
import pandas as pd
import os

PDB_PATH = "data/structures/1IVO.pdb"
INTEGRATED_PATH = "data/structures/egfr_integrated_analysis.csv"
INTERFACE_PATH = "data/structures/egfr_egf_interface.csv"

OUTPUT_PATH = "data/structures/egfr_interface_environment.csv"

CONTACT_DISTANCE = 5.0

TARGET_POSITIONS = [
    412,
    348,
    467
]

parser = PDBParser(QUIET=True)

structure = parser.get_structure(
    "EGFR",
    PDB_PATH
)

model = structure[0]

egfr_chain = model["A"]

egf_chains = []

for chain_id in ["C", "D"]:

    if chain_id in model:

        egf_chains.append(
            model[chain_id]
        )

integrated = pd.read_csv(
    INTEGRATED_PATH
)

interface = pd.read_csv(
    INTERFACE_PATH
)

integrated_map = {}

for _, row in integrated.iterrows():

    position = int(row["position"])

    integrated_map[position] = row

interface_map = {}

for _, row in interface.iterrows():

    position = int(row["EGFR_position"])

    if position not in interface_map:

        interface_map[position] = []

    interface_map[position].append(row)

egfr_residue_map = {}

for residue in egfr_chain:

    hetflag, resseq, icode = residue.id

    if hetflag == " ":

        egfr_residue_map[resseq] = residue


def calculate_distance(atom1, atom2):

    return atom1 - atom2


def residue_distance(residue1, residue2):

    minimum_distance = float("inf")

    closest_atom1 = None
    closest_atom2 = None

    for atom1 in residue1.get_atoms():

        for atom2 in residue2.get_atoms():

            distance = calculate_distance(
                atom1,
                atom2
            )

            if distance < minimum_distance:

                minimum_distance = distance

                closest_atom1 = atom1
                closest_atom2 = atom2

    return (
        minimum_distance,
        closest_atom1,
        closest_atom2
    )


results = []

print("=" * 90)
print("EGFR-EGF INTERFACE ENVIRONMENT ANALYSIS")
print("=" * 90)

print()

print(
    "Target positions:",
    TARGET_POSITIONS
)

for position in TARGET_POSITIONS:

    print()
    print("-" * 90)

    integrated_row = integrated_map.get(position)

    if integrated_row is None:

        print(
            "WARNING: position not found:",
            position
        )

        continue

    human = integrated_row["human"]
    mouse = integrated_row["mouse"]
    mutation_type = integrated_row["mutation_type"]

    pdb_position = int(
        integrated_row["1IVO_PDB_number"]
    )

    residue = egfr_residue_map.get(
        pdb_position
    )

    print(
        f"{position}: {human}->{mouse}"
    )

    print(
        "Mutation type:",
        mutation_type
    )

    print(
        "PDB position:",
        pdb_position
    )

    if residue is None:

        print(
            "WARNING: PDB residue not found"
        )

        continue

    print(
        "PDB residue:",
        residue.get_resname()
    )

    interface_rows = interface_map.get(
        position,
        []
    )

    if interface_rows:

        print()
        print(
            "DIRECT INTERFACE CONTACTS"
        )

        for interface_row in interface_rows:

            print(
                f"EGF {interface_row['EGF_chain']}:"
                f"{int(interface_row['EGF_position'])} "
                f"{interface_row['EGF_residue']} "
                f"{interface_row['minimum_distance']:.2f} Å"
            )

    else:

        print(
            "No direct interface contact"
        )

    egfr_neighbors = []

    for other_residue in egfr_chain:

        hetflag, resseq, icode = other_residue.id

        if hetflag != " ":

            continue

        if resseq == pdb_position:

            continue

        distance, atom1, atom2 = residue_distance(
            residue,
            other_residue
        )

        if distance <= CONTACT_DISTANCE:

            egfr_neighbors.append(
                (
                    other_residue,
                    distance,
                    atom1,
                    atom2
                )
            )

    egfr_neighbors.sort(
        key=lambda x: x[1]
    )

    print()
    print(
        "EGFR LOCAL ENVIRONMENT"
    )

    for (
        neighbor,
        distance,
        atom1,
        atom2
    ) in egfr_neighbors[:15]:

        print(
            f"{neighbor.id[1]:>4} "
            f"{neighbor.get_resname():>3} "
            f"{distance:>6.2f} Å "
            f"{atom1.get_name():>4} "
            f"{atom2.get_name():>4}"
        )

        results.append({

            "position": position,
            "human": human,
            "mouse": mouse,
            "mutation_type": mutation_type,
            "PDB_position": pdb_position,
            "PDB_residue": residue.get_resname(),
            "interaction_group": "EGFR",
            "partner_chain": "A",
            "partner_position": neighbor.id[1],
            "partner_residue": neighbor.get_resname(),
            "distance": round(distance, 2),
            "candidate_atom": atom1.get_name(),
            "partner_atom": atom2.get_name()

        })

    egf_contacts = []

    for chain in egf_chains:

        for egf_residue in chain:

            hetflag, resseq, icode = egf_residue.id

            if hetflag != " ":

                continue

            distance, atom1, atom2 = residue_distance(
                residue,
                egf_residue
            )

            if distance <= CONTACT_DISTANCE:

                egf_contacts.append(
                    (
                        chain.id,
                        egf_residue,
                        distance,
                        atom1,
                        atom2
                    )
                )

    egf_contacts.sort(
        key=lambda x: x[2]
    )

    print()
    print(
        "EGF LOCAL ENVIRONMENT"
    )

    for (
        chain_id,
        egf_residue,
        distance,
        atom1,
        atom2
    ) in egf_contacts:

        print(
            f"Chain {chain_id} "
            f"{egf_residue.id[1]:>4} "
            f"{egf_residue.get_resname():>3} "
            f"{distance:>6.2f} Å "
            f"{atom1.get_name():>4} "
            f"{atom2.get_name():>4}"
        )

        results.append({

            "position": position,
            "human": human,
            "mouse": mouse,
            "mutation_type": mutation_type,
            "PDB_position": pdb_position,
            "PDB_residue": residue.get_resname(),
            "interaction_group": "EGF",
            "partner_chain": chain_id,
            "partner_position": egf_residue.id[1],
            "partner_residue": egf_residue.get_resname(),
            "distance": round(distance, 2),
            "candidate_atom": atom1.get_name(),
            "partner_atom": atom2.get_name()

        })

df = pd.DataFrame(
    results
)

os.makedirs(
    os.path.dirname(OUTPUT_PATH),
    exist_ok=True
)

df.to_csv(
    OUTPUT_PATH,
    index=False
)

print()
print("=" * 90)
print("INTERFACE ENVIRONMENT ANALYSIS SAVED")
print("=" * 90)

print(
    "Output:",
    OUTPUT_PATH
)

print(
    "Records:",
    len(df)
)