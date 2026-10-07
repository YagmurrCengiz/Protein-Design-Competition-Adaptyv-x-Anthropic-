from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
STRUCTURES = ROOT / "data/structures"

MECHANISTIC_PATH = STRUCTURES / "egfr_candidate_mechanistic_evidence.csv"
FINAL_EVIDENCE_PATH = STRUCTURES / "egfr_final_candidate_evidence.csv"
INTERFACE_DIFFERENCES_PATH = STRUCTURES / "egfr_interface_differences.csv"
INTERFACE_ENVIRONMENT_PATH = STRUCTURES / "egfr_interface_environment.csv"
MUTATION_ANALYSIS_PATH = STRUCTURES / "egfr_mutation_analysis.csv"
INTEGRATED_ANALYSIS_PATH = STRUCTURES / "egfr_integrated_analysis.csv"
CANDIDATE_INTERACTIONS_PATH = STRUCTURES / "egfr_candidate_interactions.csv"
MUTATION_ENVIRONMENT_PATH = STRUCTURES / "egfr_mutation_environment.csv"
EGF_INTERFACE_PATH = STRUCTURES / "egfr_egf_interface.csv"
EGF_INTERACTIONS_PATH = STRUCTURES / "egfr_egf_interface_interactions.csv"
OUTPUT_PATH = STRUCTURES / "egfr_candidate_mechanistic_comparison.csv"

EXPECTED_POSITIONS = {
    126, 255, 323, 330, 348, 361, 364, 377, 383, 412, 414, 442, 467
}


def partner_summary(group):
    return "; ".join(
        f"{row.partner_chain}:{int(row.partner_position)} {row.partner_residue} "
        f"({float(row.distance):.2f} A; {row.interaction_group})"
        for row in group.sort_values("distance").itertuples()
    )


mechanistic = pd.read_csv(MECHANISTIC_PATH)
final_evidence = pd.read_csv(FINAL_EVIDENCE_PATH)
mutation_analysis = pd.read_csv(MUTATION_ANALYSIS_PATH).set_index("position", drop=False)
integrated = pd.read_csv(INTEGRATED_ANALYSIS_PATH).set_index("position", drop=False)
candidate_interactions = pd.read_csv(CANDIDATE_INTERACTIONS_PATH).set_index("position", drop=False)
mutation_environment = pd.read_csv(MUTATION_ENVIRONMENT_PATH).set_index("position", drop=False)
alphafold_neighbor_records = pd.read_csv(STRUCTURES / "egfr_mutation_neighbors.csv")
interface_differences = pd.read_csv(INTERFACE_DIFFERENCES_PATH)
interface_environment = pd.read_csv(INTERFACE_ENVIRONMENT_PATH)
egf_interface = pd.read_csv(EGF_INTERFACE_PATH)
egf_interactions = pd.read_csv(EGF_INTERACTIONS_PATH)

positions = set(pd.to_numeric(mechanistic["position"], errors="coerce").dropna().astype(int))
if positions != EXPECTED_POSITIONS or len(mechanistic) != 13 or mechanistic["position"].nunique() != 13:
    raise ValueError(
        "Candidate universe mismatch: expected exactly the 13 validated positions "
        f"{sorted(EXPECTED_POSITIONS)}, got {sorted(positions)} in {len(mechanistic)} rows"
    )
if set(pd.to_numeric(final_evidence["position"], errors="coerce").dropna().astype(int)) != EXPECTED_POSITIONS:
    raise ValueError("Final evidence candidate positions do not match the validated universe")

interface_environment_groups = {
    int(position): group.copy()
    for position, group in interface_environment.groupby("position")
}
legacy_groups = {
    int(position): group.copy()
    for position, group in interface_differences.groupby("position")
}

rows = []
for _, candidate in mechanistic.sort_values("position").iterrows():
    position = int(candidate["position"])
    mutation = mutation_analysis.loc[position]
    integrated_row = integrated.loc[position]
    pdb_local = candidate_interactions.loc[position]
    final_row = final_evidence.loc[final_evidence["position"] == position].iloc[0]
    af_local = mutation_environment.loc[position] if position in mutation_environment.index else None

    # All contact records are matched using the candidate's mapped 1IVO chain A
    # residue number, never by comparing a human sequence position to a PDB id.
    pdb_position = int(mutation["1IVO_PDB_number"])
    mapped_residue = str(mutation["1IVO_PDB_residue"])
    mapped_contacts = egf_interface[
        (egf_interface["EGFR_position"] == pdb_position)
        & (egf_interface["EGFR_residue"] == mapped_residue)
        & (pd.to_numeric(egf_interface["minimum_distance"], errors="coerce") <= 5.0)
    ].copy()
    direct_contact_details = []
    for contact in mapped_contacts.itertuples():
        typed = egf_interactions[
            (egf_interactions["EGFR_position"] == pdb_position)
            & (egf_interactions["EGFR_residue"] == mapped_residue)
            & (egf_interactions["EGF_chain"] == contact.EGF_chain)
            & (egf_interactions["EGF_position"] == contact.EGF_position)
            & (egf_interactions["EGFR_atom"] == contact.EGFR_atom)
            & (egf_interactions["EGF_atom"] == contact.EGF_atom)
        ]
        if len(typed):
            interaction_type = str(typed.iloc[0]["interaction_type"])
            confidence = str(typed.iloc[0]["interaction_confidence"])
        else:
            interaction_type = "NOT_AVAILABLE"
            confidence = "NOT_AVAILABLE"
        direct_contact_details.append(
            f"EGFR {mapped_residue} {pdb_position}:{contact.EGFR_atom} - "
            f"EGF {contact.EGF_chain}:{int(contact.EGF_position)} "
            f"{contact.EGF_residue}:{contact.EGF_atom} "
            f"({float(contact.minimum_distance):.2f} A; {interaction_type}; {confidence})"
        )

    nearest_interface_position = candidate["nearest_EGFR_interface_residue_PDB_position"]
    nearest_interface_position = int(nearest_interface_position) if pd.notna(nearest_interface_position) else None
    nearest_interface_rows = (
        egf_interface[egf_interface["EGFR_position"] == nearest_interface_position]
        if nearest_interface_position is not None
        else egf_interface.iloc[0:0]
    )
    nearest_interface_partners = "; ".join(sorted({
        f"EGF {row.EGF_chain}:{int(row.EGF_position)} {row.EGF_residue} "
        f"({float(row.minimum_distance):.2f} A from interface residue)"
        for row in nearest_interface_rows.itertuples()
    }))

    environment_rows = interface_environment_groups.get(position)
    if environment_rows is None or environment_rows.empty:
        env_status = "NOT_AVAILABLE"
        env_record_count = "NOT_AVAILABLE"
        env_egfr_count = "NOT_AVAILABLE"
        env_egf_count = "NOT_AVAILABLE"
        env_partner_details = "NOT_AVAILABLE"
    else:
        env_status = "OBSERVED"
        env_record_count = len(environment_rows)
        env_egfr_count = int((environment_rows["interaction_group"] == "EGFR").sum())
        env_egf_count = int((environment_rows["interaction_group"] == "EGF").sum())
        env_partner_details = partner_summary(environment_rows)

    legacy_rows = legacy_groups.get(position)
    legacy_record_count = len(legacy_rows) if legacy_rows is not None else 0

    af_positive = "NOT_AVAILABLE"
    af_negative = "NOT_AVAILABLE"
    if af_local is not None:
        # Derive sign counts from the existing AlphaFold neighbor identities;
        # the source environment table stores charged neighbors as one total.
        candidate_neighbors = alphafold_neighbor_records[
            alphafold_neighbor_records["position"] == position
        ]
        af_positive = int(candidate_neighbors["neighbor_residue"].isin(["ARG", "LYS", "HIS"]).sum())
        af_negative = int(candidate_neighbors["neighbor_residue"].isin(["ASP", "GLU"]).sum())

    row = {
        "position": position,
        "human": candidate["human"],
        "mouse": candidate["mouse"],
        "mutation_type": mutation["mutation_type"],
        "charge_change": candidate["charge_change"],
        "polarity_change": candidate["polarity_change"],
        "hydrophobicity_delta_mouse_minus_human": candidate["hydrophobicity_delta_mouse_minus_human"],
        "side_chain_size_change": candidate["side_chain_size_change"],
        "chemical_change_summary": candidate["chemical_change_summary"],
        "domain": candidate["domain"],
        "PDB_chain": "A",
        "PDB_position": mutation["1IVO_PDB_number"],
        "PDB_residue": mutation["1IVO_PDB_residue"],
        "1IVO_SASA": mutation["1IVO_SASA"],
        "1IVO_exposure": mutation["1IVO_exposure"],
        "AlphaFold_position": mutation["AlphaFold_PDB_number"],
        "AlphaFold_residue": mutation["AlphaFold_PDB_residue"],
        "AlphaFold_SASA": mutation["AlphaFold_SASA"],
        "AlphaFold_exposure": mutation["AlphaFold_exposure"],
        "SASA_difference_AlphaFold_minus_1IVO": mutation["SASA_difference"],
        "exposure_changed": mutation["exposure_changed"],
        "AlphaFold_pLDDT": mutation["AlphaFold_pLDDT"],
        "charge_change": candidate["charge_change"],
        "polarity_change": candidate["polarity_change"],
        "hydrophobicity_delta_mouse_minus_human": candidate["hydrophobicity_delta_mouse_minus_human"],
        "side_chain_size_change": candidate["side_chain_size_change"],
        "chemical_change_summary": candidate["chemical_change_summary"],
        "1IVO_neighbor_count_within_5A": pdb_local["neighbor_count"],
        "1IVO_charged_neighbors_within_5A": int(pdb_local["positive_neighbors"] + pdb_local["negative_neighbors"]),
        "1IVO_positive_neighbors_within_5A": pdb_local["positive_neighbors"],
        "1IVO_negative_neighbors_within_5A": pdb_local["negative_neighbors"],
        "1IVO_polar_neighbors_within_5A": pdb_local["polar_neighbors"],
        "1IVO_hydrophobic_neighbors_within_5A": pdb_local["hydrophobic_neighbors"],
        "1IVO_chemical_contacts_charged_plus_polar": candidate["1IVO_chemical_contacts_charged_plus_polar"],
        "1IVO_potential_salt_bridges_heuristic": pdb_local["potential_salt_bridges"],
        "1IVO_closest_positive_neighbor": candidate["1IVO_closest_positive_neighbor_identity"],
        "1IVO_closest_positive_neighbor_distance_A": candidate["1IVO_closest_positive_neighbor_distance_A"],
        "1IVO_closest_negative_neighbor": candidate["1IVO_closest_negative_neighbor_identity"],
        "1IVO_closest_negative_neighbor_distance_A": candidate["1IVO_closest_negative_neighbor_distance_A"],
        "1IVO_closest_neighbor": candidate["1IVO_closest_nonadjacent_neighbor_identity"],
        "1IVO_min_nonadjacent_local_distance_A": candidate["1IVO_closest_nonadjacent_neighbor_distance_A"],
        "1IVO_local_partner_residues_within_5A": candidate["1IVO_neighbor_residue_identities_within_5A"],
        "AlphaFold_neighbor_count_within_5A": af_local["neighbor_count"] if af_local is not None else "NOT_AVAILABLE",
        "AlphaFold_charged_neighbors_within_5A": af_local["charged_neighbors"] if af_local is not None else "NOT_AVAILABLE",
        "AlphaFold_positive_neighbors_within_5A": af_positive,
        "AlphaFold_negative_neighbors_within_5A": af_negative,
        "AlphaFold_polar_neighbors_within_5A": af_local["polar_neighbors"] if af_local is not None else "NOT_AVAILABLE",
        "AlphaFold_hydrophobic_neighbors_within_5A": af_local["hydrophobic_neighbors"] if af_local is not None else "NOT_AVAILABLE",
        "AlphaFold_chemical_contacts": af_local["chemical_contacts"] if af_local is not None else "NOT_AVAILABLE",
        "AlphaFold_potential_salt_bridges": "NOT_AVAILABLE",
        "AlphaFold_closest_neighbor": candidate["AlphaFold_closest_nonadjacent_neighbor_identity"],
        "AlphaFold_min_nonadjacent_local_distance_A": candidate["AlphaFold_closest_nonadjacent_neighbor_distance_A"],
        "AlphaFold_local_partner_residues_within_5A": candidate["AlphaFold_neighbor_residue_identities_within_5A"],
        "actual_min_3D_distance_candidate_to_EGF_A": candidate["actual_min_3D_distance_to_EGF_A"],
        "closest_EGF_chain": candidate["closest_EGF_chain"],
        "closest_EGF_position": candidate["closest_EGF_position"],
        "closest_EGF_residue": candidate["closest_EGF_residue"],
        "closest_candidate_atom_to_EGF": final_row["closest_candidate_atom"],
        "closest_EGF_atom": final_row["closest_EGF_atom"],
        "direct_3D_interface": candidate["direct_3D_interface"],
        "near_3D_interface": candidate["near_3D_interface"],
        "legacy_sequence_near_flag_sequence_only": candidate["old_sequence_near_interface_sequence_only"],
        "interface_classification": candidate["interface_evidence_status"],
        "nearest_EGFR_interface_residue_PDB_position": nearest_interface_position if nearest_interface_position is not None else "NOT_AVAILABLE",
        "nearest_EGFR_interface_residue": candidate["nearest_EGFR_interface_residue"] if pd.notna(candidate["nearest_EGFR_interface_residue"]) else "NOT_AVAILABLE",
        "distance_candidate_to_nearest_EGFR_interface_residue_A": candidate["distance_to_nearest_EGFR_interface_residue_A"],
        "nearest_interface_residue_EGF_contacts": nearest_interface_partners or "NOT_AVAILABLE",
        "direct_EGF_contact_record_count": len(mapped_contacts),
        "direct_EGF_contact_details_mapped_PDB_residue": "; ".join(direct_contact_details) if direct_contact_details else "NO_DIRECT_CONTACT_WITHIN_5A",
        "interface_environment_status": env_status,
        "interface_environment_record_count": env_record_count,
        "interface_environment_EGFR_partner_count": env_egfr_count,
        "interface_environment_EGF_partner_count": env_egf_count,
        "interface_environment_partner_details": env_partner_details,
        "legacy_interface_difference_record_count": legacy_record_count,
        "existing_impact_score": integrated_row["impact_score"],
        "impact_score_interpretation": "EXISTING_PIPELINE_HEURISTIC_ONLY",
    }
    rows.append(row)

comparison = pd.DataFrame(rows).sort_values("position").reset_index(drop=True)
if len(comparison) != 13 or comparison["position"].nunique() != 13:
    raise ValueError("Output must contain exactly 13 unique candidate rows")
if comparison.columns.duplicated().any():
    raise ValueError("Duplicate output column names are not allowed")

comparison.to_csv(OUTPUT_PATH, index=False)

print(f"Output: {OUTPUT_PATH}")
print(f"Rows: {len(comparison)}; unique candidates: {comparison['position'].nunique()}")
print("Interface classification counts:")
print(comparison["interface_classification"].value_counts().to_string())
print()

display = comparison.copy()
display["mutation"] = display["human"].astype(str) + "->" + display["mouse"].astype(str)
display["local_environment_evidence"] = display.apply(
    lambda row: (
        f"1IVO n={row['1IVO_neighbor_count_within_5A']} "
        f"(+{row['1IVO_positive_neighbors_within_5A']}/-{row['1IVO_negative_neighbors_within_5A']}, "
        f"polar={row['1IVO_polar_neighbors_within_5A']}, hydrophobic={row['1IVO_hydrophobic_neighbors_within_5A']}); "
        f"AF n={row['AlphaFold_neighbor_count_within_5A']} "
        f"(charged={row['AlphaFold_charged_neighbors_within_5A']}, polar={row['AlphaFold_polar_neighbors_within_5A']}, "
        f"hydrophobic={row['AlphaFold_hydrophobic_neighbors_within_5A']})"
    ),
    axis=1,
)
display["salt_bridge_evidence"] = display["1IVO_potential_salt_bridges_heuristic"].map(
    lambda value: f"1IVO potential={value}; AF=NOT_AVAILABLE"
)
print(display[[
    "position", "mutation", "domain", "AlphaFold_pLDDT",
    "SASA_difference_AlphaFold_minus_1IVO", "mutation_type",
    "actual_min_3D_distance_candidate_to_EGF_A", "interface_classification",
    "local_environment_evidence", "salt_bridge_evidence", "existing_impact_score",
]].to_string(index=False))

unavailable_fields = [
    "AlphaFold_potential_salt_bridges: NOT_AVAILABLE in the existing AlphaFold environment outputs",
    "interface_environment summary: NOT_AVAILABLE for candidates without rows in egfr_interface_environment.csv",
]
print("\nUnavailable fields:")
for field in unavailable_fields:
    print(f"- {field}")
