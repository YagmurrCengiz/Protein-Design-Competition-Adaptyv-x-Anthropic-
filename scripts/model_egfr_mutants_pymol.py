"""Build the EGFR A324 S>T model and a PyMOL comparison session.

Run with the PyMOL 3.1.8 Python environment, with PYMOL_DATA pointing to the
installation's share/pymol/data directory. Every mutant is made from a fresh
load of data/structures/1IVO.pdb using PyMOL's Mutagenesis Wizard.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
STRUCTURES = ROOT / "data" / "structures"
WT_PDB = STRUCTURES / "1IVO.pdb"
MUTANTS = {
    "EGFR_442_S_G": ("418", "GLY"),
    "EGFR_442_S_A": ("418", "ALA"),
    "EGFR_442_S_T": ("418", "THR"),
    "EGFR_348_S_T": ("324", "THR"),
}

if not os.environ.get("PYMOL_SITE_PACKAGES"):
    raise RuntimeError("Set PYMOL_SITE_PACKAGES to the installed PyMOL site-packages directory")
sys.path.insert(0, os.environ["PYMOL_SITE_PACKAGES"])

import pymol  # noqa: E402
pymol.finish_launching(["pymol", "-cq"])
from pymol import cmd  # noqa: E402
from pymol.wizard import mutagenesis  # noqa: E402


def mutate_fresh(name: str, resi: str, target: str) -> None:
    cmd.reinitialize()
    cmd.load(str(WT_PDB), "wt")
    selection = f"wt and chain A and resi {resi}"
    if cmd.count_atoms(selection) == 0:
        raise RuntimeError(f"No atoms found for {selection}")
    # The wizard deletes its selection argument after constructing the
    # rotamers, so pass a named selection rather than a selection expression.
    cmd.select("target_residue", selection)
    cmd.set_wizard(mutagenesis.Mutagenesis())
    wizard = cmd.get_wizard()
    # Rotamer selection does not need PyMOL's optional sculpt bump scoring.
    # In this headless session that path raises a selection-parser exception;
    # the mutation itself still comes from the installed sidechain library.
    wizard.bump_check = 0
    wizard.set_mode(target)
    wizard.do_select("target_residue")
    if cmd.count_atoms("mutation") == 0:
        raise RuntimeError(f"Wizard did not create a rotamer object for {name}")
    # With sculpt scoring disabled, explicitly choose the first standard
    # library rotamer (wizard state 1), rather than inventing coordinates.
    state = 1
    cmd.frame(state)
    if cmd.count_atoms("mutation and chain A and resi " + resi) == 0:
        raise RuntimeError(f"Wizard rotamer is missing at {resi} for {name}")
    wizard.apply()
    cmd.set_wizard()
    cmd.save(str(STRUCTURES / f"{name}.pdb"), "wt", state=1)
    cmd.delete("all")


# Preserve all existing model files; only create A324 S>T when absent.
if not (STRUCTURES / "EGFR_348_S_T.pdb").is_file():
    mutate_fresh("EGFR_348_S_T", *MUTANTS["EGFR_348_S_T"])

# Assemble the visualization from independent objects. Existing mutant PDBs
# are loaded as-is and never overwritten here.
cmd.reinitialize()
cmd.load(str(WT_PDB), "wt")
for name in MUTANTS:
    cmd.load(str(STRUCTURES / f"{name}.pdb"), name.lower())

cmd.hide("everything", "all")
for name in ["wt", *(n.lower() for n in MUTANTS)]:
    cmd.show("cartoon", f"{name} and chain A+B")
    cmd.show("cartoon", f"{name} and chain C+D")
    cmd.color("slate", f"{name} and chain A")
    cmd.color("gray70", f"{name} and chain B")
    cmd.color("teal", f"{name} and chain C")
    cmd.color("cyan", f"{name} and chain D")
    cmd.show("sticks", f"{name} and chain A and resi 324+418 and not name N+C+O+CA")
    cmd.show("sticks", f"{name} and chain C and resi 44+45 and not name N+C+O+CA")
    cmd.color("orange", f"{name} and chain A and resi 324+418")
    cmd.color("yellow", f"{name} and chain C and resi 44+45")
    cmd.show("surface", f"{name} and chain A+C")
    cmd.set("transparency", 0.65, name)

# Keep all variants in the session, with WT visible initially.
for name in MUTANTS:
    cmd.disable(name.lower())
cmd.enable("wt")
cmd.bg_color("white")
cmd.set("ray_opaque_background", 1)
cmd.set("antialias", 2)
cmd.set("label_size", 18)
cmd.set("label_color", "black")
cmd.label("wt and chain A and resi 324 and name CA", '"EGFR 348 / PDB A:324 SER"')
cmd.label("wt and chain A and resi 418 and name CA", '"EGFR 442 / PDB A:418 SER"')
cmd.label("wt and chain C and resi 44 and name CA", '"EGF C:44 TYR"')
cmd.label("wt and chain C and resi 45 and name CA", '"EGF C:45 ARG"')

cmd.orient("wt and chain A+C")
cmd.zoom("wt and chain A+C", 8)
cmd.scene("full_complex", "store")
cmd.png(str(STRUCTURES / "EGFR_EGF_full.png"), width=2400, height=1800, dpi=300, ray=1)

cmd.disable("wt")
cmd.enable("egfr_348_s_t")
cmd.hide("surface", "all")
cmd.set("label_size", 12)
focus_348 = "(egfr_348_s_t and chain A and resi 324) or (egfr_348_s_t and chain C and resi 44)"
context_348 = f"(egfr_348_s_t and chain A+C) within 5 of ({focus_348})"
cmd.label("egfr_348_s_t and chain A and resi 324 and name CA", '"EGFR 348 A:324 THR"')
cmd.label("egfr_348_s_t and chain C and resi 44 and name CA", '"EGF C44 TYR"')
cmd.orient(focus_348)
cmd.zoom(f"({focus_348}) or ({context_348})", 3)
cmd.scene("A324_EGF_C44", "store")
cmd.png(str(STRUCTURES / "EGFR_EGF_348_interface.png"), width=2000, height=1600, dpi=300, ray=1)

cmd.disable("egfr_348_s_t")
cmd.enable("egfr_442_s_t")
focus_442 = "(egfr_442_s_t and chain A and resi 418) or (egfr_442_s_t and chain C and resi 45)"
context_442 = f"(egfr_442_s_t and chain A+C) within 5 of ({focus_442})"
cmd.label("egfr_442_s_t and chain A and resi 418 and name CA", '"EGFR 442 A:418 THR"')
cmd.label("egfr_442_s_t and chain C and resi 45 and name CA", '"EGF C45 ARG"')
cmd.orient(focus_442)
cmd.zoom(f"({focus_442}) or ({context_442})", 3)
cmd.scene("A418_EGF_C45", "store")
cmd.png(str(STRUCTURES / "EGFR_EGF_442_interface.png"), width=2000, height=1600, dpi=300, ray=1)

cmd.save(str(STRUCTURES / "EGFR_EGF_mutant_models.pse"))
print("Generated A324 S>T and saved mutant comparison session/images.")
