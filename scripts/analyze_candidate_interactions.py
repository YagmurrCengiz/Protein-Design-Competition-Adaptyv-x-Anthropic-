from Bio.PDB import PDBParser
import pandas as pd
import os
import math

PDB_PATH = "data/structures/1IVO.pdb"
INPUT_PATH = "data/structures/egfr_integrated_analysis.csv"
OUTPUT_PATH = "data/structures/egfr_candidate_interactions.csv"

CONTACT_DISTANCE = 5.0
SALT_BRIDGE_DISTANCE = 4.0

parser = PDBParser(QUIET=True)

structure = parser.get_structure(
    "EGFR",
    PDB_PATH
)

model = structure[0]
chain = model["A"]

candidates = pd.read_csv(INPUT_PATH)

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

candidate_map = {}

for residue in chain:

    hetflag, resseq, icode = residue.id

    if hetflag == " ":

        candidate_map[resseq] = residue


def atom_distance(atom1, atom2):

    return atom1 - atom2


def classify_residue(resname):

    if resname in charged_positive:
        return "POSITIVE"

    if resname in charged_negative:
        return "NEGATIVE"

    if resname in polar_residues:
        return "POLAR"

    if resname in hydrophobic_residues:
        return "HYDROPHOBIC"

    return "OTHER"


def is_salt_bridge(residue1, residue2, distance):

    name1 = residue1.get_resname()
    name2 = residue2.get_resname()

    if distance > SALT_BRIDGE_DISTANCE:
        return False

    if name1 in charged_positive and name2 in charged_negative:
        return True

    if name1 in charged_negative and name2 in charged_positive:
        return True

    return False


results = []

print("=" * 90)
print("EGFR CANDIDATE INTERACTION ANALYSIS")
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

for _, candidate in candidates.iterrows():

    position = int(candidate["position"])

    human = candidate["human"]

    mouse = candidate["mouse"]

    mutation_type = candidate["mutation_type"]

    if "PDB_position" in candidates.columns:

        pdb_position = int(candidate["PDB_position"])

    else:

        pdb_position = position - 24

    residue = candidate_map.get(pdb_position)

    print()
    print("-" * 90)

    print(
        f"{position}: "
        f"{human}->{mouse}"
    )

    print(
        "Mutation type:",
        mutation_type
    )

    if residue is None:

        print(
            "WARNING: residue not found"
        )

        continue

    print(
        "PDB residue:",
        residue.get_resname(),
        residue.id[1]
    )

    interactions = []

    for other_residue in chain:

        hetflag, resseq, icode = other_residue.id

        if hetflag != " ":
            continue

        if resseq == pdb_position:
            continue

        minimum_distance = float("inf")

        for atom1 in residue.get_atoms():

            for atom2 in other_residue.get_atoms():

                distance = atom_distance(
                    atom1,
                    atom2
                )

                if distance < minimum_distance:

                    minimum_distance = distance

        if minimum_distance <= CONTACT_DISTANCE:

            interactions.append(
                (
                    other_residue,
                    minimum_distance
                )
            )

    interactions.sort(
        key=lambda x: x[1]
    )

    positive_contacts = 0
    negative_contacts = 0
    polar_contacts = 0
    hydrophobic_contacts = 0
    salt_bridges = 0

    closest_positive_distance = None
    closest_negative_distance = None

    for other_residue, distance in interactions:

        resname = other_residue.get_resname()

        category = classify_residue(
            resname
        )

        if category == "POSITIVE":

            positive_contacts += 1

            if closest_positive_distance is None:

                closest_positive_distance = distance

        elif category == "NEGATIVE":

            negative_contacts += 1

            if closest_negative_distance is None:

                closest_negative_distance = distance

        elif category == "POLAR":

            polar_contacts += 1

        elif category == "HYDROPHOBIC":

            hydrophobic_contacts += 1

        if is_salt_bridge(
            residue,
            other_residue,
            distance
        ):

            salt_bridges += 1

    candidate_type = classify_residue(
        residue.get_resname()
    )

    print(
        "Candidate chemical class:",
        candidate_type
    )

    print(
        "Positive neighbors:",
        positive_contacts
    )

    print(
        "Negative neighbors:",
        negative_contacts
    )

    print(
        "Polar neighbors:",
        polar_contacts
    )

    print(
        "Hydrophobic neighbors:",
        hydrophobic_contacts
    )

    print(
        "Potential salt bridges:",
        salt_bridges
    )

    print()
    print(
        "Closest interactions:"
    )

    for other_residue, distance in interactions[:10]:

        print(
            f"{other_residue.id[1]:>4} "
            f"{other_residue.get_resname():>3} "
            f"{distance:>6.2f} Å "
            f"{classify_residue(other_residue.get_resname())}"
        )

    results.append({

        "position": position,

        "human": human,

        "mouse": mouse,

        "mutation_type": mutation_type,

        "PDB_position": pdb_position,

        "PDB_residue": residue.get_resname(),

        "candidate_class": candidate_type,

        "positive_neighbors": positive_contacts,

        "negative_neighbors": negative_contacts,

        "polar_neighbors": polar_contacts,

        "hydrophobic_neighbors": hydrophobic_contacts,

        "potential_salt_bridges": salt_bridges,

        "closest_positive_distance": (
            round(
                closest_positive_distance,
                2
            )
            if closest_positive_distance is not None
            else None
        ),

        "closest_negative_distance": (
            round(
                closest_negative_distance,
                2
            )
            if closest_negative_distance is not None
            else None
        ),

        "neighbor_count": len(interactions)

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
print("INTERACTION ANALYSIS SAVED")
print("=" * 90)

print(
    "Output:",
    OUTPUT_PATH
)

print(
    "Candidates analyzed:",
    len(df)
)