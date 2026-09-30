from Bio import SeqIO
from Bio.PDB import PDBParser, ShrakeRupley
from Bio.Align import PairwiseAligner
import pandas as pd
import requests
import os

os.makedirs("data/structures", exist_ok=True)

pdb_path = "data/structures/1IVO.pdb"
alphafold_path = "data/structures/alphafold_egfr.pdb"

if not os.path.exists(pdb_path):
    url = "https://files.rcsb.org/download/1IVO.pdb"
    r = requests.get(url)
    r.raise_for_status()

    with open(pdb_path, "wb") as f:
        f.write(r.content)

print("PDB loaded:", pdb_path)

human = SeqIO.read(
    "data/sequences/human_egfr.fasta",
    "fasta"
)

# Python starts with 0, that's why 24:645 since we need 25 to 645
human_extra = human.seq[24:645]

print("Human extracellular length:", len(human_extra))

mouse = SeqIO.read(
    "data/sequences/mouse_egfr.fasta",
    "fasta"
)

mouse_extra = mouse.seq[24:645]

differences = []

for position, (human_aa, mouse_aa) in enumerate(
    zip(human_extra, mouse_extra),
    start=25
):

    if human_aa != mouse_aa:

        differences.append({
            "position": position,
            "human": human_aa,
            "mouse": mouse_aa
        })

print("Total differences:", len(differences))

parser = PDBParser(QUIET=True)

three_to_one = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V"
}

# 1IVO experimental structure

structure = parser.get_structure(
    "EGFR",
    pdb_path
)

model = structure[0]

print("\nChains in 1IVO:")

for chain in model:

    print(
        chain.id,
        len(list(chain.get_residues()))
    )

# Extract EGFR chain A sequence
chain = model["A"]

pdb_sequence = ""
pdb_residues = []

for residue in chain:

    resname = residue.get_resname()

    if resname in three_to_one:

        pdb_sequence += three_to_one[resname]

        pdb_residues.append(residue)

print(
    "PDB Chain A sequence length:",
    len(pdb_sequence)
)

aligner = PairwiseAligner()

pdb_alignments = aligner.align(
    human_extra,
    pdb_sequence
)

pdb_alignment = pdb_alignments[0]

print("\nHuman EGFR ↔ PDB alignment:")
print(pdb_alignment)

# Create residue mapping
mapping = {}

human_position = 25
pdb_index = 0

aligned_human = pdb_alignment[0]
aligned_pdb = pdb_alignment[1]

for human_char, pdb_char in zip(
    aligned_human,
    aligned_pdb
):

    current_residue = None

    if pdb_char != "-":

        current_residue = pdb_residues[pdb_index]

        pdb_index += 1

    if human_char != "-":

        if current_residue is not None:

            mapping[human_position] = current_residue

        human_position += 1


# Calculate solvent accessibility
sr = ShrakeRupley()

sr.compute(
    structure,
    level="R"
)

print("\n")
print("=" * 90)
print("1IVO EXPERIMENTAL STRUCTURAL ANALYSIS")
print("=" * 90)

ivo_results = []

for difference in differences:

    position = difference["position"]

    human_aa = difference["human"]

    mouse_aa = difference["mouse"]

    residue = mapping.get(position)

    if residue is None:

        print(
            f"{position}: "
            f"{human_aa}->{mouse_aa} | "
            f"NOT MAPPED"
        )

        ivo_results.append({
            "position": position,
            "human": human_aa,
            "mouse": mouse_aa,
            "1IVO_PDB_number": None,
            "1IVO_PDB_residue": None,
            "1IVO_SASA": None,
            "1IVO_exposure": "NOT MAPPED"
        })

        continue

    pdb_number = residue.id[1]

    pdb_name = residue.get_resname()

    sasa = residue.sasa

    if sasa > 40:

        exposure = "EXPOSED"

    elif sasa > 15:

        exposure = "PARTIALLY EXPOSED"

    else:

        exposure = "BURIED"

    print(
        f"{position}: "
        f"{human_aa}->{mouse_aa} | "
        f"PDB {pdb_name} {pdb_number} | "
        f"SASA={sasa:.1f} | "
        f"{exposure}"
    )

    ivo_results.append({
        "position": position,
        "human": human_aa,
        "mouse": mouse_aa,
        "1IVO_PDB_number": pdb_number,
        "1IVO_PDB_residue": pdb_name,
        "1IVO_SASA": round(sasa, 2),
        "1IVO_exposure": exposure
    })

# ALPHAFOLD STRUCTURE

print("\n")
print("=" * 90)
print("ALPHAFOLD STRUCTURE")
print("=" * 90)

af_structure = parser.get_structure(
    "EGFR_AlphaFold",
    alphafold_path
)

af_model = af_structure[0]

print("\nChains in AlphaFold:")

for af_chain in af_model:

    print(
        af_chain.id,
        len(list(af_chain.get_residues()))
    )

# Use the first chain
af_chain = next(iter(af_model))

af_sequence = ""
af_residues = []

for residue in af_chain:

    resname = residue.get_resname()

    if resname in three_to_one:

        af_sequence += three_to_one[resname]

        af_residues.append(residue)

print(
    "AlphaFold sequence length:",
    len(af_sequence)
)

af_alignments = aligner.align(
    human_extra,
    af_sequence
)

af_alignment = af_alignments[0]

print("\nHuman EGFR ↔ AlphaFold alignment:")
print(af_alignment)

# Create AlphaFold residue mapping
af_mapping = {}

human_position = 25
af_index = 0

aligned_human = af_alignment[0]
aligned_af = af_alignment[1]

for human_char, af_char in zip(
    aligned_human,
    aligned_af
):

    current_residue = None

    if af_char != "-":

        current_residue = af_residues[af_index]

        af_index += 1

    if human_char != "-":

        if current_residue is not None:

            af_mapping[human_position] = current_residue

        human_position += 1


# Calculate AlphaFold solvent accessibility
af_sr = ShrakeRupley()

af_sr.compute(
    af_structure,
    level="R"
)

# compare 1IVO and ALPHAFOLD

print("\n")
print("=" * 100)
print("1IVO vs ALPHAFOLD COMPARISON")
print("=" * 100)

comparison_results = []

for difference in differences:

    position = difference["position"]

    human_aa = difference["human"]

    mouse_aa = difference["mouse"]

    # 1IVO

    ivo_residue = mapping.get(position)

    if ivo_residue is None:

        ivo_number = None
        ivo_name = None
        ivo_sasa = None
        ivo_exposure = "NOT MAPPED"

    else:

        ivo_number = ivo_residue.id[1]

        ivo_name = ivo_residue.get_resname()

        ivo_sasa = round(ivo_residue.sasa, 2)

        if ivo_sasa > 40:

            ivo_exposure = "EXPOSED"

        elif ivo_sasa > 15:

            ivo_exposure = "PARTIALLY EXPOSED"

        else:

            ivo_exposure = "BURIED"

    # AlphaFold

    af_residue = af_mapping.get(position)

    if af_residue is None:

        af_number = None
        af_name = None
        af_sasa = None
        af_exposure = "NOT MAPPED"
        af_plddt = None

    else:

        af_number = af_residue.id[1]

        af_name = af_residue.get_resname()

        af_sasa = round(af_residue.sasa, 2)

        af_plddt_values = [
            atom.bfactor
            for atom in af_residue.get_atoms()
        ]

        if af_plddt_values:

            af_plddt = round(
                sum(af_plddt_values) / len(af_plddt_values),
                2
            )

        else:

            af_plddt = None

        if af_sasa > 40:

            af_exposure = "EXPOSED"

        elif af_sasa > 15:

            af_exposure = "PARTIALLY EXPOSED"

        else:

            af_exposure = "BURIED"

    comparison_results.append({

        "position": position,

        "human": human_aa,

        "mouse": mouse_aa,

        "1IVO_PDB_number": ivo_number,

        "1IVO_PDB_residue": ivo_name,

        "1IVO_SASA": ivo_sasa,

        "1IVO_exposure": ivo_exposure,

        "AlphaFold_PDB_number": af_number,

        "AlphaFold_PDB_residue": af_name,

        "AlphaFold_SASA": af_sasa,

        "AlphaFold_exposure": af_exposure,

        "AlphaFold_pLDDT": af_plddt

    })

    print(
        f"{position}: "
        f"{human_aa}->{mouse_aa} | "
        f"1IVO={ivo_exposure} | "
        f"AlphaFold={af_exposure} | "
        f"pLDDT={af_plddt}"
    )

# save comparison results

df = pd.DataFrame(comparison_results)

output_path = (
    "data/structures/"
    "egfr_structure_comparison.csv"
)

df.to_csv(
    output_path,
    index=False
)

print("\n")
print("=" * 90)
print("COMPARISON SAVED")
print("=" * 90)

print(
    "Output:",
    output_path
)

print(
    "Rows:",
    len(df)
)