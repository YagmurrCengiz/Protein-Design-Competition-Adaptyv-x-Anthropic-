import requests
from pathlib import Path

OUT_DIR = Path("data/tnf_alpha")
OUT_DIR.mkdir(parents=True, exist_ok=True)

proteins = {
    "human_tnf_alpha": "P01375",
    "mouse_tnf_alpha": "P06804",
}

for name, accession in proteins.items():
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"

    response = requests.get(url, timeout=30)
    response.raise_for_status()

    fasta = response.text
    sequence = "".join(
        line.strip()
        for line in fasta.splitlines()
        if line.strip() and not line.startswith(">")
    )
    if not sequence:
        raise ValueError(f"UniProt returned no sequence for {accession}")

    output_file = OUT_DIR / f"{name}.fasta"
    output_file.write_text(fasta)

    print(f"Accession: {accession}")
    print(f"Saved to: {output_file}")
    print(f"Sequence length: {len(sequence)}")
