from Bio import SeqIO
from Bio.Align import PairwiseAligner

human = SeqIO.read(
    "data/sequences/human_egfr.fasta",
    "fasta"
)

mouse = SeqIO.read(
    "data/sequences/mouse_egfr.fasta",
    "fasta"
)

EXTRACELLULAR_START = 25
EXTRACELLULAR_END = 645

human_extra = human.seq[
    EXTRACELLULAR_START - 1 : EXTRACELLULAR_END
]

mouse_extra = mouse.seq[
    EXTRACELLULAR_START - 1 : EXTRACELLULAR_END
]

aligner = PairwiseAligner()

alignments = aligner.align(
    human_extra,
    mouse_extra
)

alignment = alignments[0]

differences = [] 

for i, (human_aa, mouse_aa) in enumerate(
    zip(human_extra, mouse_extra),
    start = EXTRACELLULAR_START 
):

    if human_aa != mouse_aa:
        differences.append({
            "position" : i,
            "human" : human_aa,
            "mouse" : mouse_aa
        })

print("\nHuman/Mouse differences")
print("=" * 40)

for diff in differences:
    print( 
        f"EGFR residue {diff['position']}: "
        f"Human = {diff['human']} "
        f"Mouse = {diff['mouse']}"
    )

print("\nTotal differences:", len(differences))

identity = (
    1 - len(differences) / len(human_extra)
) * 100

print(f"Sequence identity: {identity:.2f}%")