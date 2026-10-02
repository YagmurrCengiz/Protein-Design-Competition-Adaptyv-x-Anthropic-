from Bio.PDB import PDBParser
import pandas as pd
import os

PDB_PATH = "data/structures/1IVO.pdb"
OUTPUT_PATH = "data/structures/egfr_egf_interface.csv"

CONTACT_DISTANCE = 5.0

parser = PDBParser(QUIET=True)

structure = parser.get_structure(
    "EGFR_EGF",
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

print("=" * 90)
print("EGFR-EGF INTERFACE ANALYSIS")
print("=" * 90)

print()

print(
    "EGFR chain:",
    egfr_chain.id
)

print(
    "EGF chains:",
    [chain.id for chain in egf_chains]
)

egfr_residues = []

for residue in egfr_chain:

    hetflag, resseq, icode = residue.id

    if hetflag == " ":

        egfr_residues.append(
            residue
        )

egf_residues = []

for chain in egf_chains:

    for residue in chain:

        hetflag, resseq, icode = residue.id

        if hetflag == " ":

            egf_residues.append(
                (
                    chain.id,
                    residue
                )
            )

print()
print(
    "EGFR standard residues:",
    len(egfr_residues)
)

print(
    "EGF standard residues:",
    len(egf_residues)
)

print()
print("=" * 90)
print("SEARCHING FOR INTERFACE CONTACTS")
print("=" * 90)

results = []

for egfr_residue in egfr_residues:

    best_distance = float("inf")

    best_egf_chain = None
    best_egf_residue = None
    best_egfr_atom = None
    best_egf_atom = None

    for egf_chain_id, egf_residue in egf_residues:

        for egfr_atom in egfr_residue.get_atoms():

            for egf_atom in egf_residue.get_atoms():

                distance = egfr_atom - egf_atom

                if distance < best_distance:

                    best_distance = distance

                    best_egf_chain = egf_chain_id

                    best_egf_residue = egf_residue

                    best_egfr_atom = egfr_atom

                    best_egf_atom = egf_atom

    if best_distance <= CONTACT_DISTANCE:

        results.append({

            "EGFR_position":
                egfr_residue.id[1],

            "EGFR_residue":
                egfr_residue.get_resname(),

            "EGF_chain":
                best_egf_chain,

            "EGF_position":
                best_egf_residue.id[1],

            "EGF_residue":
                best_egf_residue.get_resname(),

            "minimum_distance":
                round(
                    best_distance,
                    2
                ),

            "EGFR_atom":
                best_egfr_atom.get_name(),

            "EGF_atom":
                best_egf_atom.get_name()

        })

print()
print(
    "Interface residue pairs:",
    len(results)
)

df = pd.DataFrame(results)

if not df.empty:

    df = df.sort_values(
        by="minimum_distance"
    )

print()
print("=" * 90)
print("EGFR-EGF INTERFACE SUMMARY")
print("=" * 90)

if not df.empty:

    print()

    print(
        "Unique EGFR interface residues:",
        df["EGFR_position"].nunique()
    )

    print(
        "Unique EGF interface residues:",
        df[
            [
                "EGF_chain",
                "EGF_position"
            ]
        ].drop_duplicates().shape[0]
    )

    print()

    print(
        df.to_string(
            index=False
        )
    )

else:

    print()
    print(
        "No EGFR-EGF interface contacts detected."
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
print("EGFR-EGF INTERFACE ANALYSIS SAVED")
print("=" * 90)

print(
    "Output:",
    OUTPUT_PATH
)