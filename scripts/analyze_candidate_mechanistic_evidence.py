from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
STRUCTURES = ROOT / "data/structures"

AUDIT_PATH = STRUCTURES / "egfr_final_candidate_evidence.csv"
PROPERTIES_PATH = STRUCTURES / "egfr_mutation_properties.csv"
INTERACTIONS_PATH = STRUCTURES / "egfr_candidate_interactions.csv"
AF_NEIGHBORS_PATH = STRUCTURES / "egfr_mutation_neighbors.csv"
INTERFACE_PATH = STRUCTURES / "egfr_egf_interface.csv"
INTERFACE_INTERACTIONS_PATH = STRUCTURES / "egfr_egf_interface_interactions.csv"
OUTPUT_PATH = STRUCTURES / "egfr_candidate_mechanistic_evidence.csv"

EXPECTED_CANDIDATES = 13

# Qualitative side-chain chemistry based on residue identity. This describes
# capability, not a validated hydrogen bond in the structure.
SIDECHAIN_CHEMISTRY = {
    "A": ("methyl", "no side-chain donor", "no side-chain acceptor"),
    "R": ("guanidinium", "donor-rich", "no side-chain acceptor"),
    "N": ("carboxamide", "donor", "acceptor"),
    "D": ("carboxylate", "no donor in carboxylate form", "acceptor"),
    "C": ("thiol", "weak donor capability", "weak acceptor capability"),
    "Q": ("carboxamide", "donor", "acceptor"),
    "E": ("carboxylate", "no donor in carboxylate form", "acceptor"),
    "G": ("no side-chain functional group", "no side-chain donor", "no side-chain acceptor"),
    "H": ("imidazole", "donor capability", "acceptor capability; protonation-dependent"),
    "I": ("branched alkyl", "no side-chain donor", "no side-chain acceptor"),
    "L": ("branched alkyl", "no side-chain donor", "no side-chain acceptor"),
    "K": ("primary amino group", "donor", "no side-chain acceptor"),
    "M": ("thioether", "no side-chain donor", "weak acceptor capability"),
    "F": ("phenyl", "no side-chain donor", "no side-chain acceptor"),
    "P": ("cyclic pyrrolidine", "no side-chain donor", "no side-chain acceptor"),
    "S": ("hydroxymethyl", "donor", "weak acceptor capability"),
    "T": ("hydroxyethyl", "donor", "weak acceptor capability"),
    "W": ("indole", "donor", "no side-chain acceptor"),
    "Y": ("phenol", "donor", "weak/context-dependent acceptor capability"),
    "V": ("branched alkyl", "no side-chain donor", "no side-chain acceptor"),
}

CHEMICAL_SUMMARIES = {
    ("D", "G"): "Loses the acidic carboxylate side chain; glycine has no side-chain donor/acceptor group and is much smaller.",
    ("R", "Q"): "Loses the positively charged guanidinium group; gains a neutral carboxamide with donor and acceptor capability.",
    ("V", "I"): "Retains a neutral nonpolar branched side chain; isoleucine is larger and changes branching geometry.",
    ("E", "D"): "Retains the acidic carboxylate and acceptor capability; shortens the side chain by one methylene group.",
    ("S", "T"): "Retains a hydroxyl donor/weak acceptor; adds a methyl group and increases side-chain size.",
    ("N", "Y"): "Changes a carboxamide to an aromatic phenol; donor capability remains, while acceptor character and side-chain geometry change.",
    ("S", "A"): "Loses the hydroxyl group and its donor/weak acceptor capability; becomes a small methyl side chain.",
    ("R", "K"): "Retains positive charge and donor capability; replaces guanidinium with a primary amino group and changes side-chain geometry.",
    ("H", "R"): "Replaces imidazole with guanidinium; changes a protonation-dependent donor/acceptor group to a donor-rich positive group.",
    ("R", "W"): "Loses the positively charged guanidinium group; gains a bulky indole side chain with a donor but no side-chain acceptor.",
    ("S", "G"): "Loses the hydroxymethyl side chain and its hydroxyl donor/weak acceptor capability; glycine has no side-chain functional group.",
    ("K", "R"): "Retains positive charge and donor capability; replaces a primary amino group with guanidinium.",
}


def clean(value):
    return None if pd.isna(value) else value


def fmt_partner(row):
    return f"{int(row.neighbor_position)}:{row.neighbor_residue} ({float(row.distance):.2f} A)"


audit = pd.read_csv(AUDIT_PATH)
properties = pd.read_csv(PROPERTIES_PATH).set_index("position", drop=False)
interactions = pd.read_csv(INTERACTIONS_PATH).set_index("position", drop=False)
af_neighbors = pd.read_csv(AF_NEIGHBORS_PATH)
interface = pd.read_csv(INTERFACE_PATH)
interface_interactions = pd.read_csv(INTERFACE_INTERACTIONS_PATH)

if len(audit) != EXPECTED_CANDIDATES or audit["position"].nunique() != EXPECTED_CANDIDATES:
    raise ValueError(
        f"Expected {EXPECTED_CANDIDATES} unique audit candidates; got "
        f"{len(audit)} rows and {audit['position'].nunique()} positions"
    )

af_groups = {int(position): group.sort_values("distance") for position, group in af_neighbors.groupby("position")}
interface_classes = [
    "DIRECT_3D_INTERFACE",
    "NEAR_3D_INTERFACE",
    "SEQUENCE_NEAR_INTERFACE",
    "NONE",
]

results = []
for _, candidate in audit.iterrows():
    position = int(candidate.position)
    prop = properties.loc[position]
    if isinstance(prop, pd.DataFrame):
        prop = prop.iloc[0]
    interaction = interactions.loc[position]
    if isinstance(interaction, pd.DataFrame):
        interaction = interaction.iloc[0]

    human = str(candidate.human)
    mouse = str(candidate.mouse)
    human_group, human_hbd, human_hba = SIDECHAIN_CHEMISTRY[human]
    mouse_group, mouse_hbd, mouse_hba = SIDECHAIN_CHEMISTRY[mouse]
    chemical_summary = CHEMICAL_SUMMARIES.get(
        (human, mouse),
        f"Side-chain group changes from {human_group} to {mouse_group}.",
    )

    candidate_neighbors = af_groups.get(position)
    if candidate_neighbors is not None:
        neighbor_records = [fmt_partner(row) for row in candidate_neighbors.itertuples()]
        nonadjacent_neighbors = candidate_neighbors[
            (candidate_neighbors["neighbor_position"] - position).abs() > 1
        ]
        positive = candidate_neighbors[candidate_neighbors["neighbor_residue"].isin(["ARG", "LYS", "HIS"])]
        negative = candidate_neighbors[candidate_neighbors["neighbor_residue"].isin(["ASP", "GLU"])]
        closest_positive = fmt_partner(next(positive.head(1).itertuples())) if len(positive) else "NO POSITIVE NEIGHBOR WITHIN 5.0 A"
        closest_negative = fmt_partner(next(negative.head(1).itertuples())) if len(negative) else "NO NEGATIVE NEIGHBOR WITHIN 5.0 A"
        closest_af_neighbor = neighbor_records[0] if neighbor_records else None
        closest_nonadjacent_neighbor = (
            fmt_partner(next(nonadjacent_neighbors.head(1).itertuples()))
            if len(nonadjacent_neighbors)
            else None
        )
        closest_nonadjacent_distance = (
            float(nonadjacent_neighbors.iloc[0]["distance"])
            if len(nonadjacent_neighbors)
            else pd.NA
        )
        af_neighbor_status = "OBSERVED"
    else:
        neighbor_records = []
        closest_positive = closest_negative = closest_af_neighbor = "NOT AVAILABLE: AlphaFold local environment was not analyzed"
        closest_nonadjacent_neighbor = None
        closest_nonadjacent_distance = pd.NA
        af_neighbor_status = "NOT AVAILABLE"

    pdb_position = int(candidate.PDB_position) if pd.notna(candidate.PDB_position) else None
    direct_rows = (
        interface[interface["EGFR_position"] == pdb_position]
        if pdb_position is not None
        else interface.iloc[0:0]
    )
    direct_contact_descriptions = []
    interaction_descriptions = []
    for row in direct_rows.itertuples():
        direct_contact_descriptions.append(
            f"EGFR {row.EGFR_residue} {row.EGFR_position}:{row.EGFR_atom} — "
            f"EGF {row.EGF_chain}:{row.EGF_position} {row.EGF_residue}:{row.EGF_atom} "
            f"({float(row.minimum_distance):.2f} A)"
        )
        typed = interface_interactions[
            (interface_interactions["EGFR_position"] == row.EGFR_position)
            & (interface_interactions["EGF_chain"] == row.EGF_chain)
            & (interface_interactions["EGF_position"] == row.EGF_position)
            & (interface_interactions["EGFR_atom"] == row.EGFR_atom)
            & (interface_interactions["EGF_atom"] == row.EGF_atom)
        ]
        if len(typed):
            hit = typed.iloc[0]
            interaction_descriptions.append(
                f"{hit['EGFR_class']}–{hit['EGF_class']}; "
                f"existing classifier: {hit['interaction_type']} ({hit['interaction_confidence']})"
            )

    near_pdb_position = clean(candidate.closest_EGFR_interface_PDB_position)
    near_pdb_position = int(near_pdb_position) if near_pdb_position is not None else None
    near_interface_rows = (
        interface[interface["EGFR_position"] == near_pdb_position]
        if near_pdb_position is not None
        else interface.iloc[0:0]
    )
    near_partners = sorted({
        f"EGF {row.EGF_chain}:{int(row.EGF_position)} {row.EGF_residue}"
        for row in near_interface_rows.itertuples()
    })
    near_residue_name = (
        str(near_interface_rows.iloc[0]["EGFR_residue"])
        if len(near_interface_rows)
        else None
    )

    status = str(candidate.interface_evidence_status)
    if status == "DIRECT_3D_INTERFACE":
        contact_observation = "; ".join(direct_contact_descriptions) or "Direct contact observed from mapped 3D coordinates; atom-pair detail unavailable."
        if (human, mouse) == ("S", "G") and any("OG" in text for text in direct_contact_descriptions):
            direct_consequence = (
                "Observed: the human Ser side-chain OG is the listed atom contacting EGF Arg NH1. "
                "Derived: Ser-to-Gly removes that side-chain OG, so this exact side-chain atom contact cannot be retained; "
                "a hydrogen bond is not established by the available data."
            )
        else:
            direct_consequence = "The mapped contact is observed; the mutation-specific interaction change is not resolved by the available data."
        interface_interpretation = "Direct atom-level proximity is observed in 1IVO; chemical consequences are limited to the stated side-chain capability change."
    elif status == "NEAR_3D_INTERFACE":
        contact_observation = "No direct candidate-to-EGF contact within 5.0 A; candidate is within 5.0 A of a separate EGFR interface residue."
        near_target = f"EGFR chain A {near_residue_name} {near_pdb_position}" if near_residue_name else "EGFR interface residue"
        partner_text = ", ".join(near_partners) if near_partners else "EGF partner identity not available"
        direct_consequence = "No direct EGF contact is attributed to this candidate."
        interface_interpretation = (
            f"Observed: {float(candidate.min_3D_distance_to_EGFR_interface_residue_A):.2f} A to {near_target}, "
            f"which contacts {partner_text}. The mutation may affect nearby packing or local chemistry; a binding effect is not established."
        )
    elif status == "SEQUENCE_NEAR_INTERFACE":
        contact_observation = "No direct or tested 3D-near EGF relationship; only the legacy sequence-position proximity flag is present."
        direct_consequence = "No 3D interface consequence is inferred from the sequence-near flag."
        interface_interpretation = "Sequence-position proximity only; this is not 3D interface evidence."
    else:
        contact_observation = "No tested direct EGF contact or 3D-near interface relationship was detected."
        direct_consequence = "No EGF/interface interaction consequence is inferred; interpret the local structural and chemical evidence separately."
        interface_interpretation = "No tested direct or 3D-near EGF relationship was detected under the stated criteria."

    hbond_status = "NOT AVAILABLE: no hydrogen-bond validation or protonation analysis is present in the repository data."
    if status == "DIRECT_3D_INTERFACE" and direct_rows.empty:
        hbond_status += " The direct geometric contact has no matching atom-pair detail in the interface contact table."

    chemistry_basis = "DERIVED_FROM_OBSERVED_RESIDUE_IDENTITIES_AND_EXISTING_RESIDUE_PROPERTY_TABLE"
    results.append({
        "position": position,
        "human": human,
        "mouse": mouse,
        "mutation_type": candidate.mutation_type,
        "domain": candidate.domain,
        "PDB_chain": "A",
        "PDB_position": candidate.PDB_position,
        "PDB_residue": candidate.PDB_residue,
        "AlphaFold_position": candidate.AlphaFold_position,
        "AlphaFold_residue": candidate.AlphaFold_residue,
        "AlphaFold_pLDDT": candidate.AlphaFold_pLDDT,
        "1IVO_SASA": prop["1IVO_SASA"],
        "AlphaFold_SASA": prop["AlphaFold_SASA"],
        "SASA_difference": candidate.SASA_difference,
        "1IVO_exposure": prop["1IVO_exposure"],
        "AlphaFold_exposure": prop["AlphaFold_exposure"],
        "exposure_changed": candidate.exposure_changed,
        "human_charge": prop["human_charge"],
        "mouse_charge": prop["mouse_charge"],
        "charge_change": prop["charge_change"],
        "human_polarity": prop["human_polarity"],
        "mouse_polarity": prop["mouse_polarity"],
        "polarity_change": prop["polarity_change"],
        "hydrophobicity_delta_mouse_minus_human": prop["hydrophobicity_change"],
        "side_chain_size_change": prop["size_change"],
        "human_side_chain_group": human_group,
        "mouse_side_chain_group": mouse_group,
        "human_side_chain_H_bond_donor_capability": human_hbd,
        "mouse_side_chain_H_bond_donor_capability": mouse_hbd,
        "human_side_chain_H_bond_acceptor_capability": human_hba,
        "mouse_side_chain_H_bond_acceptor_capability": mouse_hba,
        "functional_group_change": f"{human_group} -> {mouse_group}",
        "chemical_change_summary": chemical_summary,
        "hydrogen_bond_interaction_status": hbond_status,
        "mutant_3D_geometry_status": "NOT AVAILABLE: no human-to-mouse mutant structure was modeled.",
        "AlphaFold_neighbor_count_within_5A": len(candidate_neighbors) if candidate_neighbors is not None else pd.NA,
        "AlphaFold_charged_neighbors": int(candidate.AlphaFold_charged_neighbors) if pd.notna(candidate.AlphaFold_charged_neighbors) else pd.NA,
        "AlphaFold_polar_neighbors": int(candidate.AlphaFold_polar_neighbors) if pd.notna(candidate.AlphaFold_polar_neighbors) else pd.NA,
        "AlphaFold_hydrophobic_neighbors": int(candidate.AlphaFold_hydrophobic_neighbors) if pd.notna(candidate.AlphaFold_hydrophobic_neighbors) else pd.NA,
        "AlphaFold_special_neighbors": int(candidate.AlphaFold_special_neighbors) if pd.notna(candidate.AlphaFold_special_neighbors) else pd.NA,
        "AlphaFold_chemical_contacts_charged_plus_polar": (
            int(candidate.AlphaFold_charged_neighbors + candidate.AlphaFold_polar_neighbors)
            if pd.notna(candidate.AlphaFold_charged_neighbors) and pd.notna(candidate.AlphaFold_polar_neighbors)
            else pd.NA
        ),
        "AlphaFold_closest_neighbor_distance_A": candidate.AlphaFold_local_min_neighbor_distance,
        "AlphaFold_closest_neighbor_identity": closest_af_neighbor,
        "AlphaFold_closest_nonadjacent_neighbor_distance_A": closest_nonadjacent_distance,
        "AlphaFold_closest_nonadjacent_neighbor_identity": closest_nonadjacent_neighbor,
        "AlphaFold_neighbor_scope_note": "Existing <=5 A analysis includes sequence-adjacent residues; nonadjacent minimum is also reported separately.",
        "AlphaFold_closest_positive_neighbor": closest_positive,
        "AlphaFold_closest_negative_neighbor": closest_negative,
        "AlphaFold_neighbor_residue_identities_within_5A": "; ".join(neighbor_records),
        "1IVO_local_neighbor_count_within_5A": candidate["1IVO_local_EGFR_neighbor_count_5A"],
        "1IVO_charged_neighbors_within_5A": int(interaction.positive_neighbors + interaction.negative_neighbors),
        "1IVO_positive_neighbors_within_5A": interaction.positive_neighbors,
        "1IVO_negative_neighbors_within_5A": interaction.negative_neighbors,
        "1IVO_polar_neighbors_within_5A": interaction.polar_neighbors,
        "1IVO_hydrophobic_neighbors_within_5A": interaction.hydrophobic_neighbors,
        "1IVO_chemical_contacts_charged_plus_polar": int(interaction.positive_neighbors + interaction.negative_neighbors + interaction.polar_neighbors),
        "1IVO_potential_salt_bridges": interaction.potential_salt_bridges,
        "1IVO_closest_neighbor_distance_A": candidate["1IVO_local_min_neighbor_distance_A"],
        "1IVO_closest_neighbor_identity": (
            f"{int(candidate['1IVO_closest_neighbor_PDB_position'])}:{candidate['1IVO_closest_neighbor_residue']}"
            if pd.notna(candidate["1IVO_closest_neighbor_PDB_position"])
            else None
        ),
        "1IVO_closest_nonadjacent_neighbor_identity": (
            f"{int(candidate['1IVO_closest_nonadjacent_neighbor_PDB_position'])}:"
            f"{candidate['1IVO_closest_nonadjacent_neighbor_residue']}"
            if pd.notna(candidate["1IVO_closest_nonadjacent_neighbor_PDB_position"])
            else None
        ),
        "1IVO_closest_nonadjacent_neighbor_distance_A": candidate["1IVO_local_min_nonadjacent_neighbor_distance_A"],
        "1IVO_neighbor_scope_note": candidate["1IVO_local_neighbor_scope_note"],
        "1IVO_closest_positive_neighbor_identity": (
            f"{int(candidate['1IVO_closest_positive_neighbor_PDB_position'])}:"
            f"{candidate['1IVO_closest_positive_neighbor_residue']}"
            if pd.notna(candidate["1IVO_closest_positive_neighbor_PDB_position"])
            else "NO POSITIVE NEIGHBOR WITHIN 5.0 A"
        ),
        "1IVO_closest_positive_neighbor_distance_A": candidate["1IVO_closest_positive_neighbor_distance_A"],
        "1IVO_closest_negative_neighbor_identity": (
            f"{int(candidate['1IVO_closest_negative_neighbor_PDB_position'])}:"
            f"{candidate['1IVO_closest_negative_neighbor_residue']}"
            if pd.notna(candidate["1IVO_closest_negative_neighbor_PDB_position"])
            else "NO NEGATIVE NEIGHBOR WITHIN 5.0 A"
        ),
        "1IVO_closest_negative_neighbor_distance_A": candidate["1IVO_closest_negative_neighbor_distance_A"],
        "1IVO_neighbor_residue_identities_within_5A": candidate["1IVO_local_neighbor_residue_identities_within_5A"],
        "1IVO_EGF_neighbor_count_within_5A": candidate["1IVO_local_EGF_neighbor_count_5A"],
        "actual_min_3D_distance_to_EGF_A": candidate.actual_min_3D_distance_to_EGF_A,
        "closest_EGF_chain": candidate.closest_EGF_chain,
        "closest_EGF_position": candidate.closest_EGF_position,
        "closest_EGF_residue": candidate.closest_EGF_residue,
        "direct_3D_interface": candidate.direct_3D_interface,
        "near_3D_interface": candidate.near_3D_interface,
        "old_sequence_near_interface_sequence_only": candidate.old_sequence_near_interface,
        "interface_evidence_status": status,
        "nearest_EGFR_interface_residue_PDB_position": near_pdb_position,
        "nearest_EGFR_interface_residue": near_residue_name,
        "distance_to_nearest_EGFR_interface_residue_A": candidate.min_3D_distance_to_EGFR_interface_residue_A,
        "nearest_interface_residue_EGF_partners": "; ".join(near_partners),
        "direct_contact_atom_pairs_observed": (
            "; ".join(direct_contact_descriptions)
            if direct_contact_descriptions
            else "NOT APPLICABLE: no direct 3D contact observed"
        ),
        "direct_contact_type_summary": (
            "; ".join(sorted(set(interaction_descriptions)))
            if interaction_descriptions
            else "NOT APPLICABLE: no direct 3D contact observed"
        ),
        "direct_contact_mutation_consequence": direct_consequence,
        "interface_mechanistic_interpretation": interface_interpretation,
        "interface_contact_observation": contact_observation,
        "basic_mutation_evidence_status": "OBSERVED_IN_EXISTING_MUTATION_ANALYSIS",
        "structural_mapping_evidence_status": candidate.PDB_mapping_status,
        "solvent_accessibility_evidence_status": "OBSERVED_IN_EXISTING_STRUCTURE_ANALYSIS",
        "chemical_change_evidence_status": chemistry_basis,
        "local_environment_1IVO_evidence_status": "OBSERVED_IN_EXISTING_1IVO_CANDIDATE_INTERACTION_ANALYSIS",
        "local_environment_AlphaFold_evidence_status": af_neighbor_status,
        "interface_3D_evidence_status": candidate.EGF_distance_analysis_status,
        "sequence_proximity_evidence_status": "OBSERVED_LEGACY_SEQUENCE_ANALYSIS_ONLY",
        "overall_evidence_coverage": candidate.overall_evidence_coverage,
    })

result = pd.DataFrame(results).sort_values("position").reset_index(drop=True)
if len(result) != EXPECTED_CANDIDATES or result["position"].nunique() != EXPECTED_CANDIDATES:
    raise ValueError("Mechanistic evidence output must contain exactly 13 unique candidate rows")

result.to_csv(OUTPUT_PATH, index=False)

print(f"Mechanistic evidence output: {OUTPUT_PATH}")
print(f"Rows: {len(result)}; unique positions: {result['position'].nunique()}")
print("Interface-class distribution:")
print(result["interface_evidence_status"].value_counts().reindex(interface_classes, fill_value=0).to_string())

for interface_class in interface_classes:
    print(f"\n{interface_class}")
    subset = result[result["interface_evidence_status"] == interface_class]
    if subset.empty:
        print("  (none)")
        continue
    for _, row in subset.iterrows():
        print(
            f"  {int(row['position'])} | {row['human']}->{row['mouse']} | {row['domain']} | "
            f"pLDDT {float(row['AlphaFold_pLDDT']):.2f} | SASA delta {float(row['SASA_difference']):.2f} | "
            f"{row['chemical_change_summary']} | EGF {float(row['actual_min_3D_distance_to_EGF_A']):.2f} A | "
            f"{row['interface_evidence_status']} | AF nonadjacent min "
            f"{float(row['AlphaFold_closest_nonadjacent_neighbor_distance_A']):.2f} A, "
            f"AF charged/polar neighbors {int(row['AlphaFold_charged_neighbors'])}/"
            f"{int(row['AlphaFold_polar_neighbors'])}; 1IVO nonadjacent closest "
            f"{row['1IVO_closest_nonadjacent_neighbor_identity']} at "
            f"{float(row['1IVO_closest_nonadjacent_neighbor_distance_A']):.2f} A"
        )
