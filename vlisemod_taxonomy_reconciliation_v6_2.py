#!/usr/bin/env python3
"""
V-LiSEMOD taxonomy reconciliation — Pass 6.2
Manual adjudication of the 10 Pass-6 source-conflict structures (READ ONLY)

Purpose
-------
Convert the corrected Pass-6.1 forensic result into explicit occurrence-level
adjudications for the 10 structures in the source-conflict cohort.

This pass DOES NOT modify:
- production SQLite
- Stage-09 / Stage-12 / Stage-14
- structure_classifications
- PROTACability scores
- source CIF/mmCIF files
- API/UI

The adjudications are exact PDB-level rules derived from manual structure review.
No fuzzy matching and no folder labels are used as deciding evidence.

Inputs
------
--pass6-csv
    taxonomy_source_conflict_forensics_pass6.csv from Pass 6.1

Outputs
-------
taxonomy_source_conflict_adjudication_pass6_2.csv
taxonomy_source_conflict_adjudication_summary.json
taxonomy_source_conflict_retained_viral_pass6_2.csv
taxonomy_source_conflict_excluded_nonviral_pass6_2.csv
taxonomy_source_conflict_unexpected_review_pass6_2.csv
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

import pandas as pd


def txt(v) -> str:
    if v is None or pd.isna(v):
        return ""
    return str(v).strip()


def resolve_col(df, candidates, required=True):
    by = {c.casefold(): c for c in df.columns}
    for c in candidates:
        if c.casefold() in by:
            return by[c.casefold()]
    if required:
        raise KeyError(
            f"Missing required column; expected one of {candidates}. "
            f"Available columns: {list(df.columns)}"
        )
    return None


# Exact reviewed adjudications.
#
# IMPORTANT:
# These are occurrence-context decisions, not assertions that the entire PDB
# belongs to only one organism. Mixed/chimeric structures are handled according
# to the ligand-contacting target context identified upstream.
ADJUDICATIONS = {
    # Retain as viral targets.
    "5W1X": {
        "decision": "RETAIN_VIRAL_TARGET",
        "target_browser_eligible": "YES",
        "canonical_target_id": "hpv_18_l1",
        "canonical_target_name": "Major capsid protein L1",
        "target_family": "capsid_protein",
        "entity_role": "VIRAL_TARGET",
        "reason": (
            "HPV18 major capsid protein L1. Earlier conflict was caused by "
            "virus/source alias handling rather than a nonviral target."
        ),
        "evidence_url": "https://www.rcsb.org/structure/5W1X",
    },
    "3ISN": {
        "decision": "RETAIN_VIRAL_TARGET",
        "target_browser_eligible": "YES",
        "canonical_target_id": "hiv_1_reverse_transcriptase",
        "canonical_target_name": "Reverse transcriptase",
        "target_family": "reverse_transcriptase",
        "entity_role": "VIRAL_TARGET",
        "reason": "HIV-1 reverse transcriptase inhibitor complex.",
        "evidence_url": "https://www.rcsb.org/structure/3ISN",
    },
    "3ITH": {
        "decision": "RETAIN_VIRAL_TARGET",
        "target_browser_eligible": "YES",
        "canonical_target_id": "hiv_1_reverse_transcriptase",
        "canonical_target_name": "Reverse transcriptase",
        "target_family": "reverse_transcriptase",
        "entity_role": "VIRAL_TARGET",
        "reason": "HIV-1 reverse transcriptase inhibitor complex.",
        "evidence_url": "https://www.rcsb.org/structure/3ITH",
    },
    "1GZL": {
        "decision": "RETAIN_VIRAL_TARGET",
        "target_browser_eligible": "YES",
        "canonical_target_id": "hiv_1_gp41",
        "canonical_target_name": "gp41",
        "target_family": "envelope_glycoprotein",
        "entity_role": "VIRAL_TARGET_WITH_ENGINEERED_COMPONENT",
        "reason": (
            "Structure targets the HIV-1 gp41 hydrophobic pocket; the entry also "
            "contains engineered/chimeric sequence context."
        ),
        "evidence_url": "https://www.rcsb.org/structure/1GZL",
    },

    # Exclude from viral target grouping for these ligand-contact occurrences.
    "6VQY": {
        "decision": "EXCLUDE_NONVIRAL_CONTACT_CONTEXT",
        "target_browser_eligible": "NO",
        "canonical_target_id": "",
        "canonical_target_name": "",
        "target_family": "",
        "entity_role": "HOST_PARTNER",
        "reason": (
            "HLA-B*27:05 presenting an HIV-1 peptide. The contacted protein "
            "context is human MHC, not HIV reverse transcriptase."
        ),
        "evidence_url": "https://www.rcsb.org/structure/6VQY",
    },
    "6VQZ": {
        "decision": "EXCLUDE_NONVIRAL_CONTACT_CONTEXT",
        "target_browser_eligible": "NO",
        "canonical_target_id": "",
        "canonical_target_name": "",
        "target_family": "",
        "entity_role": "HOST_PARTNER",
        "reason": (
            "HLA-B*27:05 presenting an HIV-1 peptide. The contacted protein "
            "context is human MHC, not HIV reverse transcriptase."
        ),
        "evidence_url": "https://www.rcsb.org/structure/6VQZ",
    },
    "7KGO": {
        "decision": "EXCLUDE_NONVIRAL_CONTACT_CONTEXT",
        "target_browser_eligible": "NO",
        "canonical_target_id": "",
        "canonical_target_name": "",
        "target_family": "",
        "entity_role": "HOST_PARTNER",
        "reason": (
            "Human HLA-A*02:01 presenting a SARS-CoV-2 nucleoprotein peptide. "
            "The contacted protein entity is human MHC."
        ),
        "evidence_url": "https://www.rcsb.org/structure/7KGO",
    },
    "8CMF": {
        "decision": "EXCLUDE_NONVIRAL_CONTACT_CONTEXT",
        "target_browser_eligible": "NO",
        "canonical_target_id": "",
        "canonical_target_name": "",
        "target_family": "",
        "entity_role": "HOST_PARTNER",
        "reason": (
            "Human HLA-DR1 presenting a SARS-CoV-2 nsp3 epitope. The contacted "
            "protein entity is human HLA rather than nsp3."
        ),
        "evidence_url": "https://www.rcsb.org/structure/8CMF",
    },
    "2X2D": {
        "decision": "EXCLUDE_NONVIRAL_CONTACT_CONTEXT",
        "target_browser_eligible": "NO",
        "canonical_target_id": "",
        "canonical_target_name": "",
        "target_family": "",
        "entity_role": "HOST_PARTNER",
        "reason": (
            "Human cyclophilin A/HIV-1 capsid complex; the conflicting contacted "
            "entity is human peptidyl-prolyl cis-trans isomerase A."
        ),
        "evidence_url": "https://www.rcsb.org/structure/2X2D",
    },
    "1K5M": {
        "decision": "EXCLUDE_NONVIRAL_CONTACT_CONTEXT",
        "target_browser_eligible": "NO",
        "canonical_target_id": "",
        "canonical_target_name": "",
        "target_family": "",
        "entity_role": "OTHER_VIRUS_CHIMERIC_SCAFFOLD",
        "reason": (
            "Human rhinovirus 14 capsid displaying an HIV-1 V3-loop insert. "
            "The conflicting contacted coat-protein entity belongs to rhinovirus, "
            "not an HIV protein target."
        ),
        "evidence_url": "https://www.rcsb.org/structure/1K5M",
    },
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--pass6-csv",
        type=Path,
        required=True,
        help="Pass-6.1 taxonomy_source_conflict_forensics_pass6.csv",
    )
    ap.add_argument("--outdir", type=Path, required=True)
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.pass6_csv, dtype=str, low_memory=False).fillna("")
    c_pdb = resolve_col(df, ["pdb_id", "pdb"])

    rows = []
    for _, r in df.iterrows():
        rec = dict(r)
        pdb = txt(r[c_pdb]).upper()

        rule = ADJUDICATIONS.get(pdb)
        if rule is None:
            rec.update({
                "pass6_2_decision": "UNEXPECTED_PDB_REVIEW",
                "pass6_2_target_browser_eligible": "REVIEW",
                "pass6_2_canonical_target_id": "",
                "pass6_2_canonical_target_name": "",
                "pass6_2_target_family": "",
                "pass6_2_entity_role": "UNRESOLVED",
                "pass6_2_reason": "PDB was not in the reviewed 10-structure adjudication set.",
                "pass6_2_evidence_url": "",
            })
        else:
            rec.update({
                "pass6_2_decision": rule["decision"],
                "pass6_2_target_browser_eligible": rule["target_browser_eligible"],
                "pass6_2_canonical_target_id": rule["canonical_target_id"],
                "pass6_2_canonical_target_name": rule["canonical_target_name"],
                "pass6_2_target_family": rule["target_family"],
                "pass6_2_entity_role": rule["entity_role"],
                "pass6_2_reason": rule["reason"],
                "pass6_2_evidence_url": rule["evidence_url"],
            })
        rows.append(rec)

    out = pd.DataFrame(rows)
    out = out.sort_values(
        ["pass6_2_decision", c_pdb],
        kind="mergesort",
    )

    out.to_csv(
        args.outdir / "taxonomy_source_conflict_adjudication_pass6_2.csv",
        index=False,
    )

    retained = out[out["pass6_2_decision"] == "RETAIN_VIRAL_TARGET"].copy()
    excluded = out[
        out["pass6_2_decision"] == "EXCLUDE_NONVIRAL_CONTACT_CONTEXT"
    ].copy()
    unexpected = out[
        out["pass6_2_decision"] == "UNEXPECTED_PDB_REVIEW"
    ].copy()

    retained.to_csv(
        args.outdir / "taxonomy_source_conflict_retained_viral_pass6_2.csv",
        index=False,
    )
    excluded.to_csv(
        args.outdir / "taxonomy_source_conflict_excluded_nonviral_pass6_2.csv",
        index=False,
    )
    unexpected.to_csv(
        args.outdir / "taxonomy_source_conflict_unexpected_review_pass6_2.csv",
        index=False,
    )

    counts = Counter(out["pass6_2_decision"].tolist())
    pdb_counts = {
        decision: int(grp[c_pdb].nunique())
        for decision, grp in out.groupby("pass6_2_decision")
    }

    summary = {
        "input_conflict_occurrences": int(len(out)),
        "distinct_input_pdbs": int(out[c_pdb].nunique()),
        "decision_counts": dict(sorted(counts.items())),
        "distinct_pdb_counts_by_decision": dict(sorted(pdb_counts.items())),
        "retained_viral_occurrences": int(len(retained)),
        "excluded_nonviral_occurrences": int(len(excluded)),
        "unexpected_review_occurrences": int(len(unexpected)),
        "reviewed_pdb_rule_count": len(ADJUDICATIONS),
        "production_data_modified": False,
        "policy": {
            "exact_pdb_adjudication_only": True,
            "folder_label_used_as_deciding_evidence": False,
            "fuzzy_matching": False,
            "scientific_scores_modified": False,
            "production_target_browser_modified": False,
        },
    }

    (args.outdir / "taxonomy_source_conflict_adjudication_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    print(json.dumps(summary, indent=2, sort_keys=True))
    print(f"\nWrote audit-only Pass-6.2 outputs to: {args.outdir}")


if __name__ == "__main__":
    main()
