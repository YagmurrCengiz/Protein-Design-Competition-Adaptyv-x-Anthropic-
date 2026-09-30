import requests
from pathlib import Path

OUT_DIR = Path("data/sequences")
OUT_DIR.mkdir(parents = True, exist_ok = True)

proteins = {
    "human_egfr": "P00533",
    "mouse_egfr": "Q01279",
}

for name, accession in proteins.items():

    url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"

    response = requests.get(url)
    response.raise_for_status()

    fasta = response.text

    output_file = OUT_DIR / f"{name}.fasta"
    output_file.write_text(fasta)

    print(f"Downloaded {name}: {accession}")
    print(f"Saved to: {output_file}")