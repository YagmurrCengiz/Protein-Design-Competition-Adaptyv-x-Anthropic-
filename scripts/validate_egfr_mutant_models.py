"""Validate saved EGFR mutant models against 1IVO and write a compact CSV.

Clash criterion: count nonbonded heavy-atom pairs closer than 2.0 A where
one atom is in the mutated residue. Atoms in the mutated residue and its two
immediate sequence neighbors (same chain, PDB numbering +/-1) are excluded
from the opposing set, avoiding covalent/near-covalent local backbone pairs.
This is a geometric flag only; it is not an energy or stability calculation.
"""
from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser

ROOT = Path(__file__).resolve().parents[1]
STRUCTURES = ROOT / "data" / "structures"
WT_PATH = STRUCTURES / "1IVO.pdb"
CLASH_CUTOFF_A = 2.0
TARGET_BACKBONE_TOLERANCE_A = 0.15
CASES = [
    (442, 418, "S", "G", "EGFR_442_S_G.pdb"),
    (442, 418, "S", "A", "EGFR_442_S_A.pdb"),
    (442, 418, "S", "T", "EGFR_442_S_T.pdb"),
    (348, 324, "S", "T", "EGFR_348_S_T.pdb"),
]
THREE_TO_ONE = {"SER": "S", "GLY": "G", "ALA": "A", "THR": "T"}
parser = PDBParser(QUIET=True)
wt = parser.get_structure("WT", str(WT_PATH))[0]


def heavy_atoms(residue):
    return [atom for atom in residue if atom.element.strip().upper() not in {"H", "D"}]


def min_egf_distance(residue, chain_c):
    pairs = [(a, b, float(np.linalg.norm(a.coord - b.coord)))
             for a in heavy_atoms(residue)
             for r in chain_c
             for b in heavy_atoms(r)]
    atom_a, atom_b, distance = min(pairs, key=lambda x: x[2])
    r = atom_b.get_parent()
    return distance, f"{r.resname} {r.id[1]} {atom_b.name}", f"{atom_a.name}-{atom_b.name}"


def atom_key(atom):
    residue = atom.get_parent()
    chain = residue.get_parent()
    return chain.id, residue.id[1], residue.id[2], atom.name, atom.get_altloc()


def validate_coordinates(model, pdb_name, pdb_residue):
    target = model["A"][(' ', pdb_residue, ' ')]
    assert target.resname in THREE_TO_ONE, f"Unexpected target residue {target.resname}"
    backbone = {a.name for a in target}
    missing_bb = {"N", "CA", "C", "O"} - backbone
    if missing_bb:
        raise ValueError(f"{pdb_name}: missing backbone atoms {sorted(missing_bb)}")

    keys = [atom_key(a) for a in model.get_atoms()]
    duplicates = [k for k, count in Counter(keys).items() if count > 1]
    if duplicates:
        raise ValueError(f"{pdb_name}: duplicate atom keys {duplicates[:5]}")

    # Every atom outside the targeted residue must match WT by identity and
    # coordinate (1e-3 A tolerance). The target backbone must also be retained.
    actual = {atom_key(a): a for a in model.get_atoms()}
    changed = []
    for a in wt.get_atoms():
        k = atom_key(a)
        if k[0] == "A" and k[1] == pdb_residue and k[2] == " ":
            continue
        b = actual.get(k)
        if b is None or np.max(np.abs(a.coord - b.coord)) > 0.001:
            changed.append(k)
    if changed:
        raise ValueError(f"{pdb_name}: {len(changed)} non-target WT atoms missing/changed; first {changed[:5]}")
    wt_target = wt["A"][(' ', pdb_residue, ' ')]
    backbone_shift = max(float(np.linalg.norm(wt_target[name].coord - target[name].coord))
                         for name in ("N", "CA", "C", "O"))
    if backbone_shift > TARGET_BACKBONE_TOLERANCE_A:
        raise ValueError(f"{pdb_name}: target backbone moved {backbone_shift:.3f} A")
    return target, backbone_shift


def clash_count(target, model, pdb_residue):
    excluded = {pdb_residue - 1, pdb_residue, pdb_residue + 1}
    other = [a for chain in model for residue in chain
             if not (chain.id == "A" and residue.id[0] == " " and residue.id[1] in excluded)
             for a in heavy_atoms(residue)]
    n = 0
    for a in heavy_atoms(target):
        for b in other:
            if np.linalg.norm(a.coord - b.coord) < CLASH_CUTOFF_A:
                n += 1
    return n


rows = []
for human_position, pdb_residue, wt_one, mutant_one, filename in CASES:
    path = STRUCTURES / filename
    if not path.is_file():
        raise FileNotFoundError(path)
    mutant_model = parser.get_structure(filename, str(path))[0]
    target, backbone_shift = validate_coordinates(mutant_model, filename, pdb_residue)
    wt_target = wt["A"][(' ', pdb_residue, ' ')]
    if THREE_TO_ONE[wt_target.resname] != wt_one:
        raise ValueError(f"WT A:{pdb_residue} expected {wt_one}, observed {wt_target.resname}")
    if THREE_TO_ONE[target.resname] != mutant_one:
        raise ValueError(f"{filename}: expected {mutant_one}, observed {target.resname}")
    wt_distance, wt_contact, wt_atoms = min_egf_distance(wt_target, wt["C"])
    mutant_distance, mutant_contact, mutant_atoms = min_egf_distance(target, mutant_model["C"])
    clashes = clash_count(target, mutant_model, pdb_residue)
    rows.append({
        "human_position": human_position,
        "pdb_residue": f"A:{pdb_residue}",
        "mutation": f"{wt_one}>{mutant_one}",
        "mutant_pdb": f"data/structures/{filename}",
        "wt_residue": wt_target.resname,
        "mutant_residue": target.resname,
        "wt_min_egfC_distance_A": f"{wt_distance:.3f}",
        "mutant_min_egfC_distance_A": f"{mutant_distance:.3f}",
        "wt_closest_egf_residue": wt_contact,
        "mutant_closest_egf_residue": mutant_contact,
        "wt_closest_atom_pair": wt_atoms,
        "mutant_closest_atom_pair": mutant_atoms,
        "clash_count": clashes,
        "validation_status": "PASS",
        "notes": f"< {CLASH_CUTOFF_A:.1f} A heavy-atom pairs; excludes A:{pdb_residue-1}, A:{pdb_residue}, A:{pdb_residue+1}; non-target coordinates match WT; max target backbone shift {backbone_shift:.3f} A",
    })

out = STRUCTURES / "egfr_mutant_model_validation.csv"
with out.open("w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
for row in rows:
    print(f"{row['mutant_pdb']}: {row['wt_min_egfC_distance_A']} -> {row['mutant_min_egfC_distance_A']} A; clashes={row['clash_count']}; PASS")
print(f"Wrote {out}")
