"""Create a descriptive review of the existing 13-candidate evidence table.

This script reads only the completed mechanistic comparison CSV and writes a
separate review CSV. It does not modify upstream data or compute a score.
"""

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
INPUT_PATH = ROOT / "data/structures/egfr_candidate_mechanistic_comparison.csv"
OUTPUT_PATH = ROOT / "data/structures/egfr_candidate_mechanistic_review.csv"
EXPECTED_POSITIONS = {
    126, 255, 323, 330, 348, 361, 364, 377, 383, 412, 414, 442, 467
}


def available(value):
    return pd.notna(value) and str(value).strip() not in {"", "NOT_AVAILABLE"}


def display(value):
    return str(value) if available(value) else "NOT_AVAILABLE"


def local_evidence(row):
    keys = [
        "1IVO_neighbor_count_within_5A",
        "1IVO_charged_neighbors_within_5A",
        "1IVO_polar_neighbors_within_5A",
        "1IVO_hydrophobic_neighbors_within_5A",
        "1IVO_chemical_contacts_charged_plus_polar",
        "1IVO_potential_salt_bridges_heuristic",
    ]
    if not all(available(row[key]) for key in keys):
        return "NOT_AVAILABLE"
    return (
        f"1IVO within 5 A: neighbors={int(row[keys[0]])}, "
        f"charged={int(row[keys[1]])}, polar={int(row[keys[2]])}, "
        f"hydrophobic={int(row[keys[3]])}, chemical contacts={int(row[keys[4]])}; "
        f"potential salt-bridge heuristic={int(row[keys[5]])}"
    )


def charge_evidence(row):
    change = display(row["charge_change"])
    count = row["1IVO_potential_salt_bridges_heuristic"]
    if not available(count):
        bridge_context = "potential salt-bridge context NOT_AVAILABLE"
    elif int(count) > 0:
        bridge_context = f"{int(count)} potential salt-bridge context(s), heuristic only"
    else:
        bridge_context = "no potential salt-bridge context reported by the heuristic"
    return f"{change}; {bridge_context}"


def interface_class(row):
    value = display(row["interface_classification"])
    mapping = {
        "DIRECT_3D_INTERFACE": "DIRECT_3D_INTERFACE",
        "NEAR_3D_INTERFACE": "NEAR_3D_INTERFACE",
        "SEQUENCE_NEAR_INTERFACE": "SEQUENCE_NEAR_INTERFACE",
        "NONE": "NONE",
    }
    if value in mapping:
        return mapping[value]
    raise ValueError(f"Unrecognized or unavailable interface classification: {value}")


def rationale(row, local):
    cls = interface_class(row)
    chemistry = display(row["chemical_change_summary"])
    exposure = display(row["exposure_changed"])
    delta = row["SASA_difference_AlphaFold_minus_1IVO"]
    parts = [chemistry]
    if available(delta):
        parts.append(f"cross-structure SASA difference={float(delta):+.2f} A2")
    if exposure == "True":
        parts.append("exposure category differs between 1IVO and AlphaFold")
    elif exposure == "False":
        parts.append("exposure category is unchanged between 1IVO and AlphaFold")
    if local != "NOT_AVAILABLE":
        parts.append(local)
    if cls == "DIRECT_3D_INTERFACE":
        parts.append("candidate has a direct EGF contact under the existing 5 A criterion")
        details = display(row["direct_EGF_contact_details_mapped_PDB_residue"])
        if details != "NO_DIRECT_CONTACT_WITHIN_5A":
            parts.append(f"mapped contact record: {details}")
    elif cls == "NEAR_3D_INTERFACE":
        parts.append("candidate is 3D-near an EGFR interface residue, without direct EGF contact")
    elif cls == "SEQUENCE_NEAR_INTERFACE":
        parts.append("sequence-near flag only; this is not 3D interface evidence")
    else:
        parts.append("no tested interface relationship detected")
    if display(row["interface_environment_status"]) == "NOT_AVAILABLE":
        parts.append("interface-environment records NOT_AVAILABLE")
    return "; ".join(parts) + "."


def strength(row, local):
    """Descriptive evidence label; deliberately independent of impact_score."""
    cls = interface_class(row)
    if cls == "DIRECT_3D_INTERFACE" and int(row["direct_EGF_contact_record_count"]) > 0:
        details = display(row["direct_EGF_contact_details_mapped_PDB_residue"])
        # The source's sole direct contact is explicitly typed OTHER/LOW;
        # retain it as direct geometric evidence without overstating certainty.
        if "; HIGH)" in details or "; MEDIUM)" in details:
            return "STRONG"
        return "MODERATE"
    if cls == "NEAR_3D_INTERFACE" and local != "NOT_AVAILABLE":
        return "MODERATE"
    if cls in {"SEQUENCE_NEAR_INTERFACE", "NONE"} and local != "NOT_AVAILABLE":
        return "LIMITED"
    return "INSUFFICIENT"


def main():
    comparison = pd.read_csv(INPUT_PATH)
    positions = pd.to_numeric(comparison["position"], errors="coerce")
    observed = set(positions.dropna().astype(int))
    if len(comparison) != 13 or comparison["position"].nunique() != 13 or observed != EXPECTED_POSITIONS:
        raise ValueError(
            "Expected exactly 13 unique validated candidate positions; "
            f"found {len(comparison)} rows and positions {sorted(observed)}"
        )

    rows = []
    for _, source in comparison.sort_values("position").iterrows():
        local = local_evidence(source)
        row = {
            "position": int(source["position"]),
            "human": source["human"],
            "mouse": source["mouse"],
            "mutation_type": display(source["mutation_type"]),
            "EGFR_domain": display(source["domain"]),
            "AlphaFold_pLDDT": display(source["AlphaFold_pLDDT"]),
            "SASA_difference": display(source["SASA_difference_AlphaFold_minus_1IVO"]),
            "exposure_changed": display(source["exposure_changed"]),
            "mutation_chemistry": display(source["chemical_change_summary"]),
            "local_structural_evidence": local,
            "charge_evidence": charge_evidence(source),
            "interface_class": interface_class(source),
            "minimum_EGF_distance": display(source["actual_min_3D_distance_candidate_to_EGF_A"]),
            "closest_EGF_residue": (
                f"{source['closest_EGF_chain']}:{source['closest_EGF_position']} "
                f"{source['closest_EGF_residue']}"
                if all(available(source[key]) for key in ["closest_EGF_chain", "closest_EGF_position", "closest_EGF_residue"])
                else "NOT_AVAILABLE"
            ),
            "interface_environment_status": display(source["interface_environment_status"]),
            "evidence_strength": strength(source, local),
            "mechanistic_rationale": rationale(source, local),
        }
        rows.append(row)

    review = pd.DataFrame(rows)
    if len(review) != 13 or review["position"].nunique() != 13:
        raise ValueError("Review output must contain exactly 13 unique candidates")
    review.to_csv(OUTPUT_PATH, index=False)

    print(f"Output: {OUTPUT_PATH}")
    print(f"Validated rows: {len(review)}; unique positions: {review['position'].nunique()}")
    print("Evidence is descriptive and does not use existing_impact_score.")
    for cls in ["DIRECT_3D_INTERFACE", "NEAR_3D_INTERFACE", "SEQUENCE_NEAR_INTERFACE", "NONE"]:
        group = review[review["interface_class"] == cls]
        print(f"\n{cls}")
        if group.empty:
            print("(no candidates)")
            continue
        print(group[["position", "human", "mouse", "evidence_strength", "local_structural_evidence"]].to_string(index=False))

    print("\nUnavailable evidence coverage:")
    print("- interface environment: NOT_AVAILABLE for candidates lacking interface-environment records")
    print("- AlphaFold potential salt-bridge analysis: NOT_AVAILABLE in source comparison")


if __name__ == "__main__":
    main()
