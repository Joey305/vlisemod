#!/usr/bin/env python3
"""
V-LiSEMOD taxonomy reconciliation — Pass 8
Exact adjudication of the final 40 polyprotein/domain occurrences (READ ONLY)

Purpose
-------
Close the Pass-7 40-occurrence / 21-PDB polyprotein review branch using exact
structure-level adjudications after review of the deposited PDB context.

This pass does NOT modify:
- production SQLite
- Stage-09 / Stage-12 / Stage-14
- structure_classifications
- PROTACability scores
- source CIF/mmCIF files
- API/UI

No fuzzy matching.
No folder-label deciding evidence.
Rules apply only to the exact reviewed PDB IDs below.

Input
-----
--pass7-csv
    taxonomy_polyprotein_forensics_pass7.csv

Outputs
-------
taxonomy_polyprotein_final_adjudication_pass8.csv
taxonomy_polyprotein_final_adjudication_summary.json
taxonomy_polyprotein_resolved_viral_pass8.csv
taxonomy_polyprotein_resolved_interface_pass8.csv
taxonomy_polyprotein_excluded_context_pass8.csv
taxonomy_polyprotein_unexpected_review_pass8.csv
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


def viral(
    target_id,
    target_name,
    family,
    role="VIRAL_TARGET",
    reason="",
    evidence_url="",
):
    return {
        "decision": "RESOLVE_VIRAL_TARGET",
        "target_browser_eligible": "YES",
        "canonical_target_id": target_id,
        "canonical_target_name": target_name,
        "target_family": family,
        "entity_role": role,
        "reason": reason,
        "evidence_url": evidence_url,
    }


def interface(
    target_id,
    target_name,
    family,
    reason="",
    evidence_url="",
):
    return {
        "decision": "RESOLVE_VIRAL_INTERFACE",
        "target_browser_eligible": "YES",
        "canonical_target_id": target_id,
        "canonical_target_name": target_name,
        "target_family": family,
        "entity_role": "VIRAL_MULTIPROTEIN_INTERFACE",
        "reason": reason,
        "evidence_url": evidence_url,
    }


def exclude(role, reason, evidence_url=""):
    return {
        "decision": "EXCLUDE_CONTEXT_FROM_MATURE_VIRAL_TARGET_BROWSER",
        "target_browser_eligible": "NO",
        "canonical_target_id": "",
        "canonical_target_name": "",
        "target_family": "",
        "entity_role": role,
        "reason": reason,
        "evidence_url": evidence_url,
    }


# ---------------------------------------------------------------------------
# Exact reviewed PDB adjudications
# ---------------------------------------------------------------------------

RULES = {
    # HIV-1 CACTD-SP1 maturation-site structures.
    "7R7P": viral(
        "hiv_1_capsid_sp1_maturation_site",
        "Capsid–SP1 maturation site",
        "capsid_protein",
        role="VIRAL_MATURATION_SITE",
        reason=(
            "Immature HIV-1 CA C-terminal-domain/SP1 lattice. Bevirimat binds "
            "the CACTD-SP1 maturation region; do not force this junction into "
            "one mature-product coordinate interval."
        ),
        evidence_url="https://www.rcsb.org/structure/7R7P",
    ),
    "7R7Q": viral(
        "hiv_1_capsid_sp1_maturation_site",
        "Capsid–SP1 maturation site",
        "capsid_protein",
        role="VIRAL_MATURATION_SITE",
        reason=(
            "Immature HIV-1 CA C-terminal-domain/SP1 lattice. The deposited "
            "construct explicitly spans the capsid-SP1 maturation junction."
        ),
        evidence_url="https://www.rcsb.org/structure/7R7Q",
    ),

    # HIV capsid.
    "2XDE": viral(
        "hiv_1_capsid",
        "Capsid (p24)",
        "capsid_protein",
        reason="Engineered HIV capsid N-terminal domain bound to PF-3450074.",
        evidence_url="https://www.rcsb.org/structure/2XDE",
    ),

    # HIV protease.
    "2WHH": viral(
        "hiv_1_protease",
        "Protease",
        "protease",
        reason="HIV-1 protease tethered-dimer product complex.",
        evidence_url="https://www.rcsb.org/structure/2WHH",
    ),
    "2FLE": viral(
        "hiv_1_protease",
        "Protease",
        "protease",
        reason="HIV-1 protease V82A inhibitor complex.",
        evidence_url="https://www.rcsb.org/structure/2FLE",
    ),
    "2PWC": viral(
        "hiv_1_protease",
        "Protease",
        "protease",
        reason="HIV-1 protease inhibitor complex; deposited entity is a Gag-Pol-derived protease construct.",
        evidence_url="https://www.rcsb.org/structure/2PWC",
    ),
    "2PWR": viral(
        "hiv_1_protease",
        "Protease",
        "protease",
        reason="HIV-1 protease inhibitor complex; deposited 99-aa Gag-Pol-derived entity is protease.",
        evidence_url="https://www.rcsb.org/structure/2PWR",
    ),
    "2QNP": viral(
        "hiv_1_protease",
        "Protease",
        "protease",
        reason="HIV-1 protease inhibitor complex.",
        evidence_url="https://www.rcsb.org/structure/2QNP",
    ),
    "2QNQ": viral(
        "hiv_1_protease",
        "Protease",
        "protease",
        reason="HIV-1 protease inhibitor complex.",
        evidence_url="https://www.rcsb.org/structure/2QNQ",
    ),

    # HIV reverse transcriptase.
    "5XN1": viral(
        "hiv_1_reverse_transcriptase",
        "Reverse transcriptase",
        "reverse_transcriptase",
        reason="HIV-1 reverse transcriptase Q151M:DNA:entecavir-triphosphate ternary complex.",
        evidence_url="https://www.rcsb.org/structure/5XN1",
    ),
    "5XN2": viral(
        "hiv_1_reverse_transcriptase",
        "Reverse transcriptase",
        "reverse_transcriptase",
        reason="HIV-1 reverse transcriptase Q151M:DNA:dGTP ternary complex.",
        evidence_url="https://www.rcsb.org/structure/5XN2",
    ),

    # HIV integrase.
    "5OI8": viral(
        "hiv_1_integrase",
        "Integrase",
        "integrase",
        reason="HIV-1 integrase allosteric-inhibitor structure.",
        evidence_url="https://www.rcsb.org/structure/5OI8",
    ),
    "5OIA": viral(
        "hiv_1_integrase",
        "Integrase",
        "integrase",
        reason="HIV-1 integrase allosteric-inhibitor structure.",
        evidence_url="https://www.rcsb.org/structure/5OIA",
    ),

    # Host TSG101 structures using HIV p6-derived peptide ligands.
    "3P9G": exclude(
        "HOST_TARGET_WITH_VIRAL_PEPTIDE",
        (
            "Human TSG101 UEV domain bound to a modified HIV-1 Gag-p6-derived "
            "peptide. Preserve as a viral-host interaction context, but do not "
            "present TSG101 as an HIV viral protein target."
        ),
        "https://www.rcsb.org/structure/3P9G",
    ),
    "3P9H": exclude(
        "HOST_TARGET_WITH_VIRAL_PEPTIDE",
        (
            "Human TSG101 UEV domain bound to a modified HIV-1 Gag-p6-derived "
            "peptide. Preserve as a viral-host interaction context, but do not "
            "present TSG101 as an HIV viral protein target."
        ),
        "https://www.rcsb.org/structure/3P9H",
    ),

    # SARS-CoV-2 main protease miniprecursor.
    "8E4R": viral(
        "sars_cov_2_nsp5",
        "Main protease (nsp5)",
        "protease",
        reason=(
            "SARS-CoV-2 main-protease H41A miniprecursor in complex with GC373; "
            "the few precursor-flanking residues should not obscure the nsp5 identity."
        ),
        evidence_url="https://www.rcsb.org/structure/8E4R",
    ),

    # SARS-CoV-2 frameshifting/ribosome contexts.
    "7O7Y": exclude(
        "VIRAL_NASCENT_PEPTIDE_TRANSLATION_CONTEXT",
        (
            "Rabbit 80S ribosome stalled on the SARS-CoV-2 frameshift element. "
            "The viral ORF1ab segment is a nascent translation product rather "
            "than a mature viral protein target."
        ),
        "https://www.rcsb.org/structure/7O7Y",
    ),
    "7O7Z": exclude(
        "VIRAL_NASCENT_PEPTIDE_TRANSLATION_CONTEXT",
        (
            "Rabbit 80S ribosome stalled on the SARS-CoV-2 frameshift element. "
            "The viral ORF1ab segment is a nascent translation product rather "
            "than a mature viral protein target."
        ),
        "https://www.rcsb.org/structure/7O7Z",
    ),
    "7O80": exclude(
        "VIRAL_NASCENT_PEPTIDE_TRANSLATION_CONTEXT",
        (
            "Rabbit 80S termination complex at the mutated SARS-CoV-2 slippery "
            "site. Retain for provenance but exclude from mature viral-protein grouping."
        ),
        "https://www.rcsb.org/structure/7O80",
    ),
    "7O81": exclude(
        "VIRAL_NASCENT_PEPTIDE_TRANSLATION_CONTEXT",
        (
            "Collided rabbit ribosome complex stalled by the SARS-CoV-2 pseudoknot. "
            "The viral segment is a nascent translation context, not a mature target."
        ),
        "https://www.rcsb.org/structure/7O81",
    ),

    # SARS-CoV-2 nsp12/nsp9 capping complex.
    "8GW1": interface(
        "sars_cov_2_nsp12_nsp9_capping_interface",
        "nsp12–nsp9 capping interface",
        "polymerase,nsp_proteins",
        reason=(
            "SARS-CoV-2 RNA-capping complex. Pass-7 direct-contact evidence "
            "contains both ORF1ab/nsp12 and nsp9 entities; preserve as a genuine "
            "viral multi-protein interface rather than collapsing it to one protein."
        ),
        evidence_url="https://www.rcsb.org/structure/8GW1",
    ),
}


EXPECTED_PDBS = set(RULES)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--pass7-csv",
        type=Path,
        required=True,
        help="taxonomy_polyprotein_forensics_pass7.csv",
    )
    ap.add_argument("--outdir", type=Path, required=True)
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.pass7_csv, dtype=str, low_memory=False).fillna("")
    c_pdb = resolve_col(df, ["pdb_id", "pdb"])
    c_virus = resolve_col(df, ["virus_name", "virus"], required=False)

    rows = []
    for _, r in df.iterrows():
        rec = dict(r)
        pdb = txt(r[c_pdb]).upper()
        rule = RULES.get(pdb)

        if rule is None:
            rec.update({
                "pass8_decision": "UNEXPECTED_PDB_REVIEW",
                "pass8_target_browser_eligible": "REVIEW",
                "pass8_canonical_target_id": "",
                "pass8_canonical_target_name": "",
                "pass8_target_family": "",
                "pass8_entity_role": "UNRESOLVED",
                "pass8_reason": "PDB is outside the exact reviewed Pass-8 rule set.",
                "pass8_evidence_url": "",
            })
        else:
            rec.update({
                "pass8_decision": rule["decision"],
                "pass8_target_browser_eligible": rule["target_browser_eligible"],
                "pass8_canonical_target_id": rule["canonical_target_id"],
                "pass8_canonical_target_name": rule["canonical_target_name"],
                "pass8_target_family": rule["target_family"],
                "pass8_entity_role": rule["entity_role"],
                "pass8_reason": rule["reason"],
                "pass8_evidence_url": rule["evidence_url"],
            })
        rows.append(rec)

    out = pd.DataFrame(rows)
    out = out.sort_values(["pass8_decision", c_pdb], kind="mergesort")

    out.to_csv(
        args.outdir / "taxonomy_polyprotein_final_adjudication_pass8.csv",
        index=False,
    )

    resolved = out[
        out["pass8_decision"] == "RESOLVE_VIRAL_TARGET"
    ].copy()
    interfaces = out[
        out["pass8_decision"] == "RESOLVE_VIRAL_INTERFACE"
    ].copy()
    excluded = out[
        out["pass8_decision"] == "EXCLUDE_CONTEXT_FROM_MATURE_VIRAL_TARGET_BROWSER"
    ].copy()
    unexpected = out[
        out["pass8_decision"] == "UNEXPECTED_PDB_REVIEW"
    ].copy()

    resolved.to_csv(
        args.outdir / "taxonomy_polyprotein_resolved_viral_pass8.csv",
        index=False,
    )
    interfaces.to_csv(
        args.outdir / "taxonomy_polyprotein_resolved_interface_pass8.csv",
        index=False,
    )
    excluded.to_csv(
        args.outdir / "taxonomy_polyprotein_excluded_context_pass8.csv",
        index=False,
    )
    unexpected.to_csv(
        args.outdir / "taxonomy_polyprotein_unexpected_review_pass8.csv",
        index=False,
    )

    input_pdbs = set(out[c_pdb].astype(str).str.upper())
    missing_reviewed_pdbs = sorted(EXPECTED_PDBS - input_pdbs)
    unexpected_input_pdbs = sorted(input_pdbs - EXPECTED_PDBS)

    decision_counts = Counter(out["pass8_decision"].tolist())
    target_counts = Counter(
        x for x in out["pass8_canonical_target_name"].tolist() if txt(x)
    )
    role_counts = Counter(out["pass8_entity_role"].tolist())

    summary = {
        "input_occurrences": int(len(out)),
        "distinct_input_pdbs": int(out[c_pdb].nunique()),
        "reviewed_exact_pdb_rules": len(RULES),
        "resolved_viral_target_occurrences": int(len(resolved)),
        "resolved_viral_interface_occurrences": int(len(interfaces)),
        "excluded_context_occurrences": int(len(excluded)),
        "unexpected_review_occurrences": int(len(unexpected)),
        "decision_counts": dict(sorted(decision_counts.items())),
        "resolved_target_counts": dict(sorted(target_counts.items())),
        "entity_role_counts": dict(sorted(role_counts.items())),
        "missing_reviewed_pdbs_from_input": missing_reviewed_pdbs,
        "unexpected_input_pdbs": unexpected_input_pdbs,
        "production_data_modified": False,
        "policy": {
            "exact_pdb_adjudication_only": True,
            "folder_label_used_as_deciding_evidence": False,
            "fuzzy_matching": False,
            "mature_target_forced_for_translation_context": False,
            "genuine_viral_interface_collapsed": False,
            "scientific_scores_modified": False,
            "production_target_browser_modified": False,
        },
    }

    (args.outdir / "taxonomy_polyprotein_final_adjudication_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    print(json.dumps(summary, indent=2, sort_keys=True))
    print(f"\nWrote audit-only Pass-8 outputs to: {args.outdir}")


if __name__ == "__main__":
    main()
