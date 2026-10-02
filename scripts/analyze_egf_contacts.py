from Bio.PDB import PDBParser
import pandas as pd
import os

PDB_PATH = "data/structures/1IVO.pdb"
INPUT_PATH = "data/structures/egfr_high_priority_analysis.csv"
OUTPUT_PATH = "data/structures/egfr_egf_contacts.csv"
CONTACT_DISTANCE = 5.0

parser = PDBParser(QUIET=True)

structure = parser.get_structure(
    "EGFR_EGF",
    PDB_PATH
)

model = structure[0]

candidates = pd.read_csv(INPUT_PATH)

position_to_pdb = {
    364: 340,
    255: 231,
    414: 390,
    126: 102
}

print("=" * 90)
print("EGFR-EGF CONTACT ANALYSIS")
print("=" * 90)

print()

print(
    "Candidates:",
    len(candidates)
)

print(
    "Positions:",
    candidates["position"].tolist()
)

print()
print("=" * 90)
print("STRUCTURE CHAINS")
print("=" * 90)

for chain in model:

    standard_residues = []

    for residue in chain:

        hetflag, resseq, icode = residue.id

        if hetflag == " ":

            standard_residues.append(residue)

    print(
        f"Chain {chain.id}: "
        f"{len(standard_residues)} standard residues"
    )

egfr_chain = model["A"]

egf_chains = []

for chain_id in ["C", "D"]:

    if chain_id in model:

        egf_chains.append(
            model[chain_id]
        )

print()
print(
    "EGF chains:",
    [chain.id for chain in egf_chains]
)

print()
print("=" * 90)
print("BUILDING EGFR RESIDUE MAP")
print("=" * 90)

egfr_residues = []

for residue in egfr_chain:

    hetflag, resseq, icode = residue.id

    if hetflag == " ":

        egfr_residues.append(residue)

print(
    "EGFR standard residues:",
    len(egfr_residues)
)

candidate_residue_map = {}

for residue in egfr_residues:

    pdb_number = residue.id[1]

    candidate_residue_map[pdb_number] = residue

egf_atoms = []

for chain in egf_chains:

    for residue in chain:

        hetflag, resseq, icode = residue.id

        if hetflag == " ":

            for atom in residue:

                egf_atoms.append(
                    (
                        chain.id,
                        residue,
                        atom
                    )
                )

print(
    "EGF atoms:",
    len(egf_atoms)
)

def calculate_distance(atom1, atom2):

    return atom1 - atom2

results = []

print()
print("=" * 90)
print("CANDIDATE CONTACT ANALYSIS")
print("=" * 90)

for _, candidate in candidates.iterrows():

    position = int(candidate["position"])
    pdb_position = position_to_pdb.get(position)

    human_aa = candidate["human"]
    mouse_aa = candidate["mouse"]
    mutation_type = candidate["mutation_type"]

    print()
    print("-" * 90)

    print(
        f"{position}: "
        f"{human_aa}->{mouse_aa}"
    )

    print(
        "Mutation type:",
        mutation_type
    )

    if pdb_position is None:

        print(
            "WARNING: PDB position mapping not found"
        )

        results.append({

            "position": position,

            "human": human_aa,

            "mouse": mouse_aa,

            "mutation_type": mutation_type,

            "PDB_position": None,

            "PDB_residue": None,

            "EGF_contact": False,

            "EGF_min_distance": None,

            "EGF_chain": None,

            "EGF_residue_number": None,

            "EGF_residue": None,

            "EGF_atom": None,

            "EGFR_atom": None

        })

        continue

    residue = candidate_residue_map.get(pdb_position)

    if residue is None:

        print(
            "WARNING: EGFR residue not found in PDB"
        )

        results.append({

            "position": position,

            "human": human_aa,

            "mouse": mouse_aa,

            "mutation_type": mutation_type,

            "PDB_position": pdb_position,

            "PDB_residue": None,

            "EGF_contact": False,

            "EGF_min_distance": None,

            "EGF_chain": None,

            "EGF_residue_number": None,

            "EGF_residue": None,

            "EGF_atom": None,

            "EGFR_atom": None

        })

        continue

    print(
        "EGFR PDB residue:",
        residue.get_resname(),
        residue.id[1]
    )

    print(
        "Human EGFR position:",
        position
    )

    print(
        "PDB position:",
        pdb_position
    )

    minimum_distance = float("inf")

    closest_egf_chain = None
    closest_egf_residue = None
    closest_egf_atom = None
    closest_egfr_atom = None

    for egfr_atom in residue.get_atoms():

        for (
            egf_chain_id,
            egf_residue,
            egf_atom
        ) in egf_atoms:

            distance = calculate_distance(
                egfr_atom,
                egf_atom
            )

            if distance < minimum_distance:

                minimum_distance = distance

                closest_egf_chain = egf_chain_id

                closest_egf_residue = egf_residue

                closest_egf_atom = egf_atom

                closest_egfr_atom = egfr_atom

    has_contact = (
        minimum_distance <= CONTACT_DISTANCE
    )

    if has_contact:

        print(
            f"EGF CONTACT FOUND: "
            f"{minimum_distance:.2f} Å"
        )

        print(
            "EGF chain:",
            closest_egf_chain
        )

        print(
            "EGF residue:",
            closest_egf_residue.get_resname(),
            closest_egf_residue.id[1]
        )

        print(
            "EGF atom:",
            closest_egf_atom.get_name()
        )

        print(
            "EGFR atom:",
            closest_egfr_atom.get_name()
        )

    else:

        print(
            f"No EGF contact within "
            f"{CONTACT_DISTANCE} Å"
        )

        print(
            f"Closest distance: "
            f"{minimum_distance:.2f} Å"
        )

    results.append({

        "position": position,

        "human": human_aa,

        "mouse": mouse_aa,

        "mutation_type": mutation_type,

        "PDB_position": pdb_position,

        "PDB_residue": residue.get_resname(),

        "EGF_contact": has_contact,

        "EGF_min_distance": round(
            minimum_distance,
            2
        ),

        "EGF_chain": closest_egf_chain,

        "EGF_residue_number": (
            closest_egf_residue.id[1]
            if closest_egf_residue
            else None
        ),

        "EGF_residue": (
            closest_egf_residue.get_resname()
            if closest_egf_residue
            else None
        ),

        "EGF_atom": (
            closest_egf_atom.get_name()
            if closest_egf_atom
            else None
        ),

        "EGFR_atom": (
            closest_egfr_atom.get_name()
            if closest_egfr_atom
            else None
        )

    })

df = pd.DataFrame(results)

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
print("EGF CONTACT SUMMARY")
print("=" * 90)

print()

print(
    "Total candidates:",
    len(df)
)

print(
    "Candidates contacting EGF:",
    int(df["EGF_contact"].sum())
)

print(
    "Candidates without EGF contact:",
    int((~df["EGF_contact"]).sum())
)

print()

print(df[[
    "position",
    "human",
    "mouse",
    "PDB_position",
    "PDB_residue",
    "EGF_contact",
    "EGF_min_distance",
    "EGF_chain",
    "EGF_residue_number",
    "EGF_residue"
]].to_string(index=False))

print()
print("=" * 90)
print("EGF CONTACT ANALYSIS SAVED")
print("=" * 90)

print(
    "Output:",
    OUTPUT_PATH
)