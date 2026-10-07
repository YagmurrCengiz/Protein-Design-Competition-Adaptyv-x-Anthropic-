from pathlib import Path

import pandas as pd
from Bio.PDB import PDBParser


ROOT = Path(__file__).resolve().parents[1]
STRUCTURES = ROOT / "data/structures"

HIGH_PRIORITY_PATH = STRUCTURES / "egfr_high_priority_analysis.csv"
INTEGRATED_PATH = STRUCTURES / "egfr_integrated_analysis.csv"
INTERFACE_PATH = STRUCTURES / "egfr_egf_interface.csv"
INTERFACE_DIFFERENCES_PATH = STRUCTURES / "egfr_interface_differences.csv"
ALPHAFOLD_NEIGHBORS_PATH = STRUCTURES / "egfr_mutation_neighbors.csv"
OUTPUT_PATH = STRUCTURES / "egfr_final_candidate_evidence.csv"
PDB_PATH = STRUCTURES / "1IVO.pdb"

CONTACT_DISTANCE = 5.0
SALT_BRIDGE_DISTANCE = 4.0

POSITIVE = {"LYS", "ARG", "HIS"}
NEGATIVE = {"ASP", "GLU"}
POLAR = {"SER", "THR", "ASN", "GLN", "TYR", "HIS", "CYS"}
HYDROPHOBIC = {"ALA", "VAL", "ILE", "LEU", "MET", "PHE", "TRP", "PRO"}
SPECIAL = {"GLY", "PRO", "CYS"}


def residue_class(resname):
    if resname in POSITIVE:
        return "positive"
    if resname in NEGATIVE:
        return "negative"
    if resname in POLAR:
        return "polar"
    if resname in HYDROPHOBIC:
        return "hydrophobic"
    if resname in SPECIAL:
        return "special"
    return "other"


def alphafold_residue_class(resname):
    # Match scripts/analyze_mutation_environment.py's residue groups.
    if resname in POSITIVE | NEGATIVE:
        return "charged"
    if resname in {"ASN", "GLN", "SER", "THR", "TYR"}:
        return "polar"
    if resname in {"ALA", "VAL", "ILE", "LEU", "MET", "PHE", "TRP"}:
        return "hydrophobic"
    if resname in {"GLY", "PRO", "CYS"}:
        return "special"
    return "other"


def residue_distance(residue_a, residue_b):
    minimum = float("inf")
    closest_atoms = (None, None)
    for atom_a in residue_a.get_atoms():
        for atom_b in residue_b.get_atoms():
            distance = atom_a - atom_b
            if distance < minimum:
                minimum = distance
                closest_atoms = (atom_a, atom_b)
    return minimum, closest_atoms


high_priority = pd.read_csv(HIGH_PRIORITY_PATH)
integrated = pd.read_csv(INTEGRATED_PATH).set_index("position", drop=False)
interface = pd.read_csv(INTERFACE_PATH)
legacy_differences = pd.read_csv(INTERFACE_DIFFERENCES_PATH)
alphafold_neighbors = pd.read_csv(ALPHAFOLD_NEIGHBORS_PATH)

# Keep the established universe: high-priority candidates plus every position
# in the legacy interface-difference output.
candidate_positions = sorted(
    set(pd.to_numeric(high_priority["position"], errors="coerce").dropna().astype(int))
    | set(pd.to_numeric(legacy_differences["position"], errors="coerce").dropna().astype(int))
)

parser = PDBParser(QUIET=True)
structure = parser.get_structure("EGFR_1IVO", PDB_PATH)
model = structure[0]
egfr_chain = model["A"]
egf_chains = [model[chain_id] for chain_id in ("C", "D") if chain_id in model]

egfr_residues = {
    residue.id[1]: residue
    for residue in egfr_chain
    if residue.id[0] == " "
}
egf_residues = [
    (chain.id, residue)
    for chain in egf_chains
    for residue in chain
    if residue.id[0] == " "
]

# The existing interface file records PDB chain A residue numbers. Validate
# those identities against chain A before using them to define the interface.
interface_pdb_positions = set()
for _, contact in interface.iterrows():
    pdb_position = int(contact["EGFR_position"])
    residue = egfr_residues.get(pdb_position)
    if residue is not None and residue.get_resname() == contact["EGFR_residue"]:
        interface_pdb_positions.add(pdb_position)

high_priority_positions = set(high_priority["position"].astype(int))
legacy_by_position = legacy_differences.groupby("position")
af_by_position = alphafold_neighbors.groupby("position")

rows = []
for position in candidate_positions:
    if position not in integrated.index:
        raise ValueError(f"Candidate position {position} missing from integrated analysis")
    source = integrated.loc[position]
    if isinstance(source, pd.DataFrame):
        source = source.iloc[0]

    pdb_position = source.get("1IVO_PDB_number")
    pdb_residue = None
    pdb_target = None
    pdb_mapping_status = "NOT_MAPPED"
    if pd.notna(pdb_position):
        pdb_position = int(pdb_position)
        pdb_target = egfr_residues.get(pdb_position)
        if pdb_target is not None:
            pdb_residue = pdb_target.get_resname()
            pdb_mapping_status = (
                "MAPPED_MATCH" if pdb_residue == source.get("1IVO_PDB_residue") else "MAPPED_RESIDUE_MISMATCH"
            )

    direct_distance = float("nan")
    closest_egf_chain = None
    closest_egf_position = None
    closest_egf_residue = None
    closest_egf_atoms = (None, None)
    egf_distance_status = "NOT_ANALYZED"
    egfr_local = []
    egf_local = []
    nearest_interface_distance = float("nan")
    nearest_interface_position = None

    if pdb_target is not None:
        egf_distance_status = "ANALYZED"
        for chain_id, egf_residue in egf_residues:
            distance, atoms = residue_distance(pdb_target, egf_residue)
            if distance < direct_distance or pd.isna(direct_distance):
                direct_distance = distance
                closest_egf_chain = chain_id
                closest_egf_position = egf_residue.id[1]
                closest_egf_residue = egf_residue.get_resname()
                closest_egf_atoms = atoms
            if distance <= CONTACT_DISTANCE:
                egf_local.append((chain_id, egf_residue, distance))

        for neighbor_position, neighbor in egfr_residues.items():
            if neighbor_position == pdb_position:
                continue
            distance, _ = residue_distance(pdb_target, neighbor)
            if distance <= CONTACT_DISTANCE:
                egfr_local.append((neighbor, distance))

        for interface_position in interface_pdb_positions:
            if interface_position == pdb_position:
                continue
            interface_residue = egfr_residues[interface_position]
            distance, _ = residue_distance(pdb_target, interface_residue)
            if distance < nearest_interface_distance or pd.isna(nearest_interface_distance):
                nearest_interface_distance = distance
                nearest_interface_position = interface_position

    direct_3d = bool(pd.notna(direct_distance) and direct_distance <= CONTACT_DISTANCE)
    near_3d = bool(
        not direct_3d
        and pd.notna(nearest_interface_distance)
        and nearest_interface_distance <= CONTACT_DISTANCE
    )

    legacy_rows = legacy_by_position.get_group(position) if position in legacy_by_position.groups else None
    old_position_match = bool(
        legacy_rows is not None
        and legacy_rows["direct_interface_match"].fillna(False).astype(bool).any()
    )
    old_sequence_near = bool(
        legacy_rows is not None
        and (~legacy_rows["direct_interface_match"].fillna(False).astype(bool)).any()
    )

    if direct_3d:
        interface_status = "DIRECT_3D_INTERFACE"
    elif near_3d:
        interface_status = "NEAR_3D_INTERFACE"
    elif old_sequence_near:
        interface_status = "SEQUENCE_NEAR_INTERFACE"
    else:
        interface_status = "NONE"

    env_classes = [residue_class(residue.get_resname()) for residue, _ in egfr_local]
    target_class = residue_class(pdb_residue) if pdb_residue else None
    closest_egfr_neighbor = min(egfr_local, key=lambda item: item[1]) if egfr_local else None
    nonadjacent_egfr_local = [
        item for item in egfr_local
        if abs(item[0].id[1] - pdb_position) > 1
    ] if pdb_target is not None else []
    closest_nonadjacent_egfr_neighbor = (
        min(nonadjacent_egfr_local, key=lambda item: item[1])
        if nonadjacent_egfr_local else None
    )
    closest_positive_egfr_neighbor = next(
        (item for item in sorted(egfr_local, key=lambda item: item[1]) if item[0].get_resname() in POSITIVE),
        None,
    )
    closest_negative_egfr_neighbor = next(
        (item for item in sorted(egfr_local, key=lambda item: item[1]) if item[0].get_resname() in NEGATIVE),
        None,
    )
    possible_salt_bridges = 0
    if pdb_residue:
        for neighbor, distance in egfr_local:
            neighbor_class = residue_class(neighbor.get_resname())
            if distance <= SALT_BRIDGE_DISTANCE and (
                (target_class == "positive" and neighbor_class == "negative")
                or (target_class == "negative" and neighbor_class == "positive")
            ):
                possible_salt_bridges += 1

    af_rows = af_by_position.get_group(position) if position in af_by_position.groups else None
    af_counts = {name: 0 for name in ("charged", "polar", "hydrophobic", "special", "other")}
    af_min_distance = float("nan")
    if af_rows is not None:
        for residue_name, count in af_rows["neighbor_residue"].map(alphafold_residue_class).value_counts().items():
            af_counts[residue_name] = int(count)
        af_min_distance = pd.to_numeric(af_rows["distance"], errors="coerce").min()

    af_position = source.get("AlphaFold_PDB_number")
    af_residue = source.get("AlphaFold_PDB_residue")
    af_metrics_complete = all(pd.notna(source.get(column)) for column in (
        "AlphaFold_PDB_number", "AlphaFold_PDB_residue", "AlphaFold_pLDDT",
        "AlphaFold_SASA", "AlphaFold_exposure", "SASA_difference"
    ))
    local_af_status = "ANALYZED" if af_rows is not None else "NOT_ANALYZED"
    local_1ivo_status = "ANALYZED" if pdb_target is not None else "NOT_ANALYZED"

    if pdb_mapping_status == "MAPPED_MATCH" and af_metrics_complete and local_af_status == "ANALYZED" and local_1ivo_status == "ANALYZED" and egf_distance_status == "ANALYZED":
        overall_coverage = "COMPLETE"
    else:
        overall_coverage = "INCOMPLETE"

    if not direct_3d and not near_3d and not old_sequence_near:
        interface_status_detail = "ANALYZED_NO_INTERFACE_RELATIONSHIP"
    else:
        interface_status_detail = "RELATIONSHIP_DETECTED"

    rows.append({
        "position": position,
        "human": source["human"],
        "mouse": source["mouse"],
        "mutation_type": source["mutation_type"],
        "domain": (
            "Domain I" if 25 <= position <= 190 else
            "Domain II" if 191 <= position <= 292 else
            "Domain III" if 293 <= position <= 444 else
            "Domain IV" if 445 <= position <= 645 else
            "Outside extracellular region"
        ),
        "PDB_position": pdb_position,
        "PDB_residue": pdb_residue,
        "AlphaFold_position": af_position,
        "AlphaFold_residue": af_residue,
        "AlphaFold_pLDDT": source.get("AlphaFold_pLDDT"),
        "SASA_difference": source.get("SASA_difference"),
        "exposure_changed": source.get("exposure_changed"),
        "high_priority_impact_score": source.get("impact_score") if position in high_priority_positions else pd.NA,
        "AlphaFold_local_neighbor_count": len(af_rows) if af_rows is not None else pd.NA,
        "AlphaFold_local_min_neighbor_distance": af_min_distance,
        "AlphaFold_charged_neighbors": af_counts["charged"] if af_rows is not None else pd.NA,
        "AlphaFold_polar_neighbors": af_counts["polar"] if af_rows is not None else pd.NA,
        "AlphaFold_hydrophobic_neighbors": af_counts["hydrophobic"] if af_rows is not None else pd.NA,
        "AlphaFold_special_neighbors": af_counts["special"] if af_rows is not None else pd.NA,
        "1IVO_local_EGFR_neighbor_count_5A": len(egfr_local) if pdb_target is not None else pd.NA,
        "1IVO_local_EGF_neighbor_count_5A": len(egf_local) if pdb_target is not None else pd.NA,
        "1IVO_local_min_neighbor_distance_A": closest_egfr_neighbor[1] if closest_egfr_neighbor else pd.NA,
        "1IVO_closest_neighbor_PDB_position": closest_egfr_neighbor[0].id[1] if closest_egfr_neighbor else pd.NA,
        "1IVO_closest_neighbor_residue": closest_egfr_neighbor[0].get_resname() if closest_egfr_neighbor else pd.NA,
        "1IVO_local_min_nonadjacent_neighbor_distance_A": closest_nonadjacent_egfr_neighbor[1] if closest_nonadjacent_egfr_neighbor else pd.NA,
        "1IVO_closest_nonadjacent_neighbor_PDB_position": closest_nonadjacent_egfr_neighbor[0].id[1] if closest_nonadjacent_egfr_neighbor else pd.NA,
        "1IVO_closest_nonadjacent_neighbor_residue": closest_nonadjacent_egfr_neighbor[0].get_resname() if closest_nonadjacent_egfr_neighbor else pd.NA,
        "1IVO_local_neighbor_scope_note": "Existing <=5 A EGFR-chain analysis includes sequence-adjacent residues; nonadjacent minimum is also reported separately.",
        "1IVO_closest_positive_neighbor_PDB_position": closest_positive_egfr_neighbor[0].id[1] if closest_positive_egfr_neighbor else pd.NA,
        "1IVO_closest_positive_neighbor_residue": closest_positive_egfr_neighbor[0].get_resname() if closest_positive_egfr_neighbor else pd.NA,
        "1IVO_closest_positive_neighbor_distance_A": closest_positive_egfr_neighbor[1] if closest_positive_egfr_neighbor else pd.NA,
        "1IVO_closest_negative_neighbor_PDB_position": closest_negative_egfr_neighbor[0].id[1] if closest_negative_egfr_neighbor else pd.NA,
        "1IVO_closest_negative_neighbor_residue": closest_negative_egfr_neighbor[0].get_resname() if closest_negative_egfr_neighbor else pd.NA,
        "1IVO_closest_negative_neighbor_distance_A": closest_negative_egfr_neighbor[1] if closest_negative_egfr_neighbor else pd.NA,
        "1IVO_local_neighbor_residue_identities_within_5A": "; ".join(
            f"{neighbor.id[1]}:{neighbor.get_resname()} ({distance:.2f} A)"
            for neighbor, distance in sorted(egfr_local, key=lambda item: item[1])
        ) if pdb_target is not None else pd.NA,
        "1IVO_local_positive_neighbors": env_classes.count("positive") if pdb_target is not None else pd.NA,
        "1IVO_local_negative_neighbors": env_classes.count("negative") if pdb_target is not None else pd.NA,
        "1IVO_local_polar_neighbors": env_classes.count("polar") if pdb_target is not None else pd.NA,
        "1IVO_local_hydrophobic_neighbors": env_classes.count("hydrophobic") if pdb_target is not None else pd.NA,
        "1IVO_local_potential_salt_bridges": possible_salt_bridges if pdb_target is not None else pd.NA,
        "actual_min_3D_distance_to_EGF_A": direct_distance,
        "closest_EGF_chain": closest_egf_chain,
        "closest_EGF_position": closest_egf_position,
        "closest_EGF_residue": closest_egf_residue,
        "closest_candidate_atom": closest_egf_atoms[0].get_name() if closest_egf_atoms[0] else None,
        "closest_EGF_atom": closest_egf_atoms[1].get_name() if closest_egf_atoms[1] else None,
        "min_3D_distance_to_EGFR_interface_residue_A": nearest_interface_distance,
        "closest_EGFR_interface_PDB_position": nearest_interface_position,
        "direct_3D_interface": direct_3d,
        "near_3D_interface": near_3d,
        "old_sequence_near_interface": old_sequence_near,
        "old_position_match_from_legacy_analysis": old_position_match,
        "interface_evidence_status": interface_status,
        "interface_analysis_status": interface_status_detail,
        "PDB_mapping_status": pdb_mapping_status,
        "AlphaFold_metric_status": "ANALYZED" if af_metrics_complete else "NOT_ANALYZED",
        "AlphaFold_local_environment_status": local_af_status,
        "1IVO_local_environment_status": local_1ivo_status,
        "EGF_distance_analysis_status": egf_distance_status,
        "sequence_proximity_analysis_status": "ANALYZED",
        "overall_evidence_coverage": overall_coverage,
    })

final = pd.DataFrame(rows).sort_values("position").reset_index(drop=True)
final.to_csv(OUTPUT_PATH, index=False)

print("EGFR candidate evidence / positioning audit")
print(f"Output: {OUTPUT_PATH}")
print(f"Candidates: {len(final)}; unique positions: {final['position'].nunique()}")
print(f"Direct 3D interface cutoff: <= {CONTACT_DISTANCE:.1f} A")
print(f"Near 3D interface cutoff: <= {CONTACT_DISTANCE:.1f} A to a separate EGFR interface residue")
print()
print(final[[
    "position", "human", "mouse", "actual_min_3D_distance_to_EGF_A",
    "closest_EGF_chain", "closest_EGF_position", "closest_EGF_residue",
    "direct_3D_interface", "near_3D_interface", "old_sequence_near_interface",
    "interface_evidence_status", "overall_evidence_coverage",
]].to_string(index=False))
print()
print("Coverage:")
print(final["overall_evidence_coverage"].value_counts(dropna=False).to_string())
