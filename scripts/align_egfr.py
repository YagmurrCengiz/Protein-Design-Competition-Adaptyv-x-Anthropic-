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


print("Human EGFR length:", len(human.seq))
print("Mouse EGFR length:", len(mouse.seq))

print(
    "Human extracellular length:",
    len(human_extra)
)

print(
    "Mouse extracellular length:",
    len(mouse_extra)
)


aligner = PairwiseAligner()

alignments = aligner.align(
    human_extra,
    mouse_extra
)

alignment = alignments[0]


print("\nAlignment:")
print(alignment)