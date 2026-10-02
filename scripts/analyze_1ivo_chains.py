from Bio.PDB import PDBParser

pdb_path = "data/structures/1IVO.pdb"

parser = PDBParser(QUIET=True)
structure = parser.get_structure("EGFR", pdb_path)

model = structure[0]

print("=" * 80)
print("1IVO CHAIN CONTENT ANALYSIS")
print("=" * 80)

for chain in model:

    print(f"\nCHAIN {chain.id}")
    print("-" * 50)

    residues = []

    for residue in chain:

        hetflag, resseq, icode = residue.id

        # Standard amino-acid residues only
        if hetflag == " ":
            residues.append(residue)

    print("Number of standard residues:", len(residues))

    if residues:

        print(
            "First residues:",
            [
                f"{r.get_resname()} {r.id[1]}"
                for r in residues[:10]
            ]
        )

        print(
            "Last residues:",
            [
                f"{r.get_resname()} {r.id[1]}"
                for r in residues[-10:]
            ]
        )

    # Show all non-standard residues
    nonstandard = []

    for residue in chain:

        hetflag, resseq, icode = residue.id

        if hetflag != " ":

            nonstandard.append(
                (
                    residue.get_resname(),
                    resseq,
                    hetflag
                )
            )

    if nonstandard:

        print("Non-standard / hetero residues:")
        for item in nonstandard[:30]:
            print(item)

print("\n")
print("=" * 80)
print("DONE")
print("=" * 80)