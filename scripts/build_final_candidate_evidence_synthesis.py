"""Synthesize existing candidate review evidence into descriptive follow-up groups.

The review CSV is the sole source. This script creates no score and does not
change upstream data, candidate selection, or analysis outputs.
"""

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
INPUT_PATH = ROOT / "data/structures/egfr_candidate_mechanistic_review.csv"
OUTPUT_PATH = ROOT / "data/structures/egfr_final_candidate_evidence_synthesis.csv"
EXPECTED_POSITIONS = {
    126, 255, 323, 330, 348, 361, 364, 377, 383, 412, 414, 442, 467
}

# Descriptive synthesis categories, based on independent interface, chemical,
# and local-structure evidence in the existing review; these are not scores.
GROUPS = {
    126: "PRIORITY FOR FOLLOW-UP",  # charge loss + 3D-near + dense local contacts
    255: "PRIORITY FOR FOLLOW-UP",  # charge loss + dense charged local environment
    348: "PRIORITY FOR FOLLOW-UP",  # 3D-near + exposure change + observed environment
    383: "PRIORITY FOR FOLLOW-UP",  # charge gain + 3D-near + local chemical context
    442: "PRIORITY FOR FOLLOW-UP",  # direct geometry + side-chain loss + exposure change
    364: "SECONDARY FOLLOW-UP",
    377: "SECONDARY FOLLOW-UP",
    412: "SECONDARY FOLLOW-UP",
    414: "SECONDARY FOLLOW-UP",
    467: "SECONDARY FOLLOW-UP",
    323: "LOW CURRENT EVIDENCE",
    330: "LOW CURRENT EVIDENCE",
    361: "LOW CURRENT EVIDENCE",
}


def interface_evidence(row):
    cls = row["interface_class"]
    if cls == "DIRECT_3D_INTERFACE":
        return (
            "Direct geometric EGF contact is recorded; interaction typing/confidence is "
            "reported in the source contact detail and should be considered."
        )
    if cls == "NEAR_3D_INTERFACE":
        return "3D-near an EGFR interface residue; no direct EGF contact is reported."
    if cls == "SEQUENCE_NEAR_INTERFACE":
        return "Legacy sequence-near flag only; this is not 3D interface evidence."
    if cls == "NONE":
        return "No tested interface relationship detected; local evidence is considered separately."
    return "NOT_AVAILABLE"


def limitations(row):
    items = []
    cls = row["interface_class"]
    if cls == "SEQUENCE_NEAR_INTERFACE":
        items.append("interface association is sequence-near only, not structural proximity")
    elif cls == "NONE":
        items.append("no tested interface relationship")
    elif cls == "NEAR_3D_INTERFACE":
        items.append("3D proximity does not establish direct EGF contact")
    elif cls == "DIRECT_3D_INTERFACE":
        detail = str(row["mechanistic_rationale"])
        if "LOW" in detail:
            items.append("direct contact is geometrically observed but source interaction confidence is LOW")
    if row["interface_environment_status"] == "NOT_AVAILABLE":
        items.append("interface-environment data NOT_AVAILABLE, not evidence of absence")
    charge = str(row["charge_evidence"])
    if "heuristic only" in charge:
        items.append("potential salt-bridge count is heuristic/contextual, not a confirmed interaction")
    items.append("SASA difference compares 1IVO with AlphaFold and is not a measured mutation effect")
    items.append("no experimental causality or functional effect is established by these fields")
    return "; ".join(items) + "."


def local_evidence(row):
    return str(row["local_structural_evidence"])


def summary(row):
    chem = str(row["mutation_chemistry"])
    cls = row["interface_class"]
    struct = str(row["structural_evidence"])
    # Keep the terminal/CSV summary readable while retaining the observed fields.
    struct = struct.replace("A2", "A²")
    import re
    struct = re.sub(r"cross-structure SASA difference=([+-]?\d+\.\d+) A²", lambda m: f"cross-structure SASA difference={float(m.group(1)):+.2f} A²", struct)
    struct = re.sub(r"AlphaFold pLDDT=([0-9]+(?:\.\d+)?)\.", lambda m: f"AlphaFold pLDDT={float(m.group(1)):.2f}", struct)
    if cls == "DIRECT_3D_INTERFACE":
        interface = "a direct geometric EGF contact is recorded, with LOW interaction confidence"
    elif cls == "NEAR_3D_INTERFACE":
        interface = "the residue is 3D-near an interface residue, without direct EGF contact"
    elif cls == "SEQUENCE_NEAR_INTERFACE":
        interface = "only legacy sequence-near evidence is present"
    else:
        interface = "no tested interface relationship is reported"
    return f"{chem} Local evidence: {struct}; {interface}."


def main():
    review = pd.read_csv(INPUT_PATH)
    positions = set(pd.to_numeric(review["position"], errors="coerce").dropna().astype(int))
    if len(review) != 13 or review["position"].nunique() != 13 or positions != EXPECTED_POSITIONS:
        raise ValueError(
            "Expected exactly 13 unique review candidates; "
            f"found {len(review)} rows with positions {sorted(positions)}"
        )
    if review.columns.duplicated().any():
        raise ValueError("Input review contains duplicate columns")

    records = []
    for _, source in review.sort_values("position").iterrows():
        record = {
            "position": int(source["position"]),
            "human": source["human"],
            "mouse": source["mouse"],
            "mutation_type": source["mutation_type"],
            "EGFR_domain": source["EGFR_domain"],
            "interface_class": source["interface_class"],
            "mutation_chemistry": source["mutation_chemistry"],
            "structural_evidence": str(source["local_structural_evidence"])
            + f"; cross-structure SASA difference={source['SASA_difference']} A2; "
            + f"exposure_changed={source['exposure_changed']}; AlphaFold pLDDT={source['AlphaFold_pLDDT']}.",
            "interface_evidence": interface_evidence(source),
            "local_chemical_evidence": (
                f"{source['local_structural_evidence']}; {source['charge_evidence']}"
            ),
            "evidence_limitations": limitations(source),
            "evidence_strength": source["evidence_strength"],
            "follow_up_group": GROUPS[int(source["position"])],
            "mechanistic_summary": "",
        }
        record["mechanistic_summary"] = summary(record)
        records.append(record)

    synthesis = pd.DataFrame(records)
    expected_columns = [
        "position", "human", "mouse", "mutation_type", "EGFR_domain",
        "interface_class", "mutation_chemistry", "structural_evidence",
        "interface_evidence", "local_chemical_evidence", "evidence_limitations",
        "evidence_strength", "follow_up_group", "mechanistic_summary",
    ]
    synthesis = synthesis[expected_columns]
    if len(synthesis) != 13 or synthesis["position"].nunique() != 13:
        raise ValueError("Synthesis output must contain exactly 13 unique candidates")
    if synthesis.columns.duplicated().any():
        raise ValueError("Synthesis output contains duplicate columns")
    if not synthesis.astype(str).apply(
        lambda column: column.str.contains("NOT_AVAILABLE", regex=False).any()
    ).any():
        raise ValueError("Unavailable evidence must remain explicitly marked")

    synthesis.to_csv(OUTPUT_PATH, index=False)
    print(f"Output: {OUTPUT_PATH}")
    print(f"Validated rows: {len(synthesis)}; unique positions: {synthesis['position'].nunique()}")
    print(f"Duplicate columns: {synthesis.columns.duplicated().any()}")
    print("No numerical ranking score was created; existing impact_score was not read.")
    for group in ["PRIORITY FOR FOLLOW-UP", "SECONDARY FOLLOW-UP", "LOW CURRENT EVIDENCE"]:
        print(f"\n{group}")
        subset = synthesis[synthesis["follow_up_group"] == group]
        for row in subset.itertuples(index=False):
            mutation = f"{row.human}->{row.mouse}"
            rationale = row.mechanistic_summary
            print(
                f"{row.position} {mutation} | {row.interface_class} | "
                f"{row.evidence_strength} | {rationale}"
            )
    print("\nUnavailable evidence remains marked in the synthesis CSV, including interface-environment gaps.")


if __name__ == "__main__":
    main()
