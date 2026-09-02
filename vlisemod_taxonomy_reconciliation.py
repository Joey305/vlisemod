#!/usr/bin/env python3
"""
V-LiSEMOD taxonomy reconciliation pass.

Purpose
-------
Build a conservative, audit-only mapping from ligand-contacting mmCIF entity
descriptions to canonical V-LiSEMOD target labels.

This script DOES NOT modify the database, Stage-09, Stage-12, Stage-14, API,
website, source classifications, or scientific scores.

Strategy
--------
1. Read classification_reconciliation_audit.csv from the prior audit.
2. Learn exact entity-description -> canonical-target mappings only from
   already trusted occurrence categories where the target is unambiguous.
3. Apply optional reviewed manual overrides.
4. Reconcile disagreement rows only when every contacted entity description
   maps deterministically.
5. Preserve genuine multi-protein interfaces and unresolved cases for review.
6. Write a taxonomy candidate table, review queue, reconciled occurrence audit,
   grouping simulation, and summary JSON.

No fuzzy matching is used.
"""

from __future__ import annotations

import argparse
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable, Optional

import pandas as pd


TRUSTED_CATEGORIES = {
    "CLEAN_AGREEMENT",
    "CLASSIFICATION_PROVENANCE_CONTAMINATION",
    "MULTIPROTEIN_STRUCTURE_SINGLE_CONTACT_TARGET",
}

GENUINE_MULTI_CATEGORIES = {
    "GENUINE_MULTIPROTEIN_LIGAND_INTERFACE",
}

NO_CONTACT_CATEGORIES = {
    "UNRESOLVED_NO_CONTACT_EVIDENCE",
}

MULTI_SPLIT_RE = re.compile(r"\s*[;,|]\s*")


def norm_text(value) -> str:
    if value is None or pd.isna(value):
        return ""
    s = str(value).strip()
    s = re.sub(r"\s+", " ", s)
    return s


def norm_key(value) -> str:
    """Stable comparison key; intentionally not biologically fuzzy."""
    s = norm_text(value).casefold()
    s = s.replace("_", " ")
    s = re.sub(r"[\[\]\(\)\{\}]", " ", s)
    s = re.sub(r"[^a-z0-9.+/\- ]+", " ", s)
    s = re.sub(r"\s+", " ", s).strip()
    return s


def split_multi(value) -> list[str]:
    s = norm_text(value)
    if not s:
        return []
    return [x.strip() for x in MULTI_SPLIT_RE.split(s) if x.strip()]


def unique_preserve(values: Iterable[str]) -> list[str]:
    out, seen = [], set()
    for v in values:
        if not v:
            continue
        k = norm_key(v)
        if k and k not in seen:
            seen.add(k)
            out.append(v)
    return out


def resolve_column(df: pd.DataFrame, candidates: list[str], required: bool = True) -> Optional[str]:
    lower = {c.casefold(): c for c in df.columns}
    for c in candidates:
        if c.casefold() in lower:
            return lower[c.casefold()]
    if required:
        raise KeyError(
            f"Could not find any of {candidates}. Available columns: {list(df.columns)}"
        )
    return None


def load_overrides(path: Optional[Path]) -> dict[tuple[str, str], str]:
    """
    CSV schema:
      virus_name,entity_description,canonical_target,note
    virus_name may be '*' for a global exact override.
    """
    if not path:
        return {}
    odf = pd.read_csv(path, dtype=str).fillna("")
    vcol = resolve_column(odf, ["virus_name", "virus"])
    ecol = resolve_column(odf, ["entity_description", "contacting_entity_description"])
    tcol = resolve_column(odf, ["canonical_target", "target", "protein_label"])
    out = {}
    for _, r in odf.iterrows():
        virus = norm_text(r[vcol])
        entity = norm_text(r[ecol])
        target = norm_text(r[tcol])
        if not entity or not target:
            continue
        out[(virus.casefold(), norm_key(entity))] = target
    return out


def single_target(value) -> str:
    vals = unique_preserve(split_multi(value))
    return vals[0] if len(vals) == 1 else ""


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--audit-csv",
        required=True,
        type=Path,
        help="Prior classification_reconciliation_audit.csv",
    )
    ap.add_argument(
        "--outdir",
        required=True,
        type=Path,
        help="Directory for audit-only taxonomy outputs",
    )
    ap.add_argument(
        "--overrides",
        type=Path,
        default=None,
        help="Optional reviewed exact-mapping CSV",
    )
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.audit_csv, dtype=str, low_memory=False).fillna("")

    # Resolve likely columns from the audit schema.
    c_pdb = resolve_column(df, ["pdb_id", "pdb"])
    c_virus = resolve_column(df, ["virus_name", "virus"])
    c_ligand = resolve_column(df, ["ligand_id", "ligand", "ligand_code"], required=False)
    c_instance = resolve_column(
        df, ["ligand_instance", "ligand_instance_id", "occurrence_id"], required=False
    )
    c_folder = resolve_column(
        df, ["folder_protein_labels", "structure_protein_labels", "classification_labels"]
    )
    c_entity_desc = resolve_column(
        df, ["contacting_entity_descriptions", "contact_entity_descriptions", "entity_descriptions"]
    )
    c_entity_ids = resolve_column(
        df, ["contacting_entity_ids", "contact_entity_ids", "entity_ids"], required=False
    )
    c_norm_contact = resolve_column(
        df,
        ["normalized_contact_protein_labels", "normalized_contact_labels"],
        required=False,
    )
    c_current = resolve_column(
        df, ["current_stage14_protein_type", "current_target_label", "protein_type"]
    )
    c_proposed = resolve_column(
        df, ["proposed_canonical_target", "proposed_target"], required=False
    )
    c_category = resolve_column(df, ["audit_category", "category"])

    overrides = load_overrides(args.overrides)

    # ------------------------------------------------------------------
    # 1) Build an entity-description evidence table.
    # ------------------------------------------------------------------
    evidence = []
    for idx, row in df.iterrows():
        virus = norm_text(row[c_virus])
        category = norm_text(row[c_category])
        descriptions = unique_preserve(split_multi(row[c_entity_desc]))

        proposed = single_target(row[c_proposed]) if c_proposed else ""
        norm_contact = single_target(row[c_norm_contact]) if c_norm_contact else ""
        trusted_target = proposed or norm_contact

        for desc in descriptions:
            evidence.append(
                {
                    "row_index": idx,
                    "pdb_id": norm_text(row[c_pdb]),
                    "virus_name": virus,
                    "entity_description": desc,
                    "entity_description_key": norm_key(desc),
                    "audit_category": category,
                    "trusted_target": trusted_target if category in TRUSTED_CATEGORIES else "",
                    "current_target": norm_text(row[c_current]),
                    "folder_labels": norm_text(row[c_folder]),
                }
            )

    ev = pd.DataFrame(evidence)
    if ev.empty:
        raise RuntimeError("No contacting entity descriptions found in the audit CSV.")

    # Virus-specific exact mappings learned from already trusted rows.
    learned_specific: dict[tuple[str, str], str] = {}
    specific_conflicts: set[tuple[str, str]] = set()

    grouped = defaultdict(set)
    for _, r in ev.iterrows():
        if r["trusted_target"]:
            grouped[(r["virus_name"].casefold(), r["entity_description_key"])].add(
                norm_text(r["trusted_target"])
            )
    for key, targets in grouped.items():
        nt = unique_preserve(sorted(targets, key=str.casefold))
        if len(nt) == 1:
            learned_specific[key] = nt[0]
        else:
            specific_conflicts.add(key)

    # Global exact mappings are allowed only if the description resolves to the
    # same target across all viruses in trusted evidence.
    global_grouped = defaultdict(set)
    for _, r in ev.iterrows():
        if r["trusted_target"]:
            global_grouped[r["entity_description_key"]].add(norm_text(r["trusted_target"]))

    learned_global: dict[str, str] = {}
    global_conflicts: set[str] = set()
    for key, targets in global_grouped.items():
        nt = unique_preserve(sorted(targets, key=str.casefold))
        if len(nt) == 1:
            learned_global[key] = nt[0]
        else:
            global_conflicts.add(key)

    # ------------------------------------------------------------------
    # 2) Candidate taxonomy summary.
    # ------------------------------------------------------------------
    candidate_rows = []
    for (virus, desc_key), grp in ev.groupby(["virus_name", "entity_description_key"], sort=True):
        display_descs = unique_preserve(grp["entity_description"].tolist())
        trusted = unique_preserve([x for x in grp["trusted_target"].tolist() if x])
        current = unique_preserve(grp["current_target"].tolist())
        folders = unique_preserve(grp["folder_labels"].tolist())
        override = overrides.get((virus.casefold(), desc_key)) or overrides.get(("*", desc_key), "")
        learned = learned_specific.get((virus.casefold(), desc_key), "")
        global_learned = learned_global.get(desc_key, "")

        if override:
            suggested = override
            source = "MANUAL_OVERRIDE"
        elif learned:
            suggested = learned
            source = "TRUSTED_VIRUS_SPECIFIC_EXACT"
        elif global_learned and desc_key not in global_conflicts:
            suggested = global_learned
            source = "TRUSTED_GLOBAL_EXACT"
        else:
            suggested = ""
            source = "UNRESOLVED"

        candidate_rows.append(
            {
                "virus_name": virus,
                "entity_description": display_descs[0] if display_descs else "",
                "occurrence_evidence_rows": int(len(grp)),
                "distinct_pdbs": int(grp["pdb_id"].nunique()),
                "trusted_targets_observed": ";".join(trusted),
                "current_targets_observed": ";".join(current),
                "folder_labels_observed": ";".join(folders),
                "suggested_canonical_target": suggested,
                "suggestion_source": source,
                "specific_mapping_conflict": (virus.casefold(), desc_key) in specific_conflicts,
                "global_mapping_conflict": desc_key in global_conflicts,
            }
        )

    candidates = pd.DataFrame(candidate_rows).sort_values(
        ["suggestion_source", "occurrence_evidence_rows", "virus_name", "entity_description"],
        ascending=[True, False, True, True],
        kind="mergesort",
    )
    candidates.to_csv(args.outdir / "canonical_target_taxonomy_candidates.csv", index=False)

    # ------------------------------------------------------------------
    # 3) Reconcile each occurrence.
    # ------------------------------------------------------------------
    def map_entity(virus: str, desc: str) -> tuple[str, str]:
        k = norm_key(desc)
        specific_override = overrides.get((virus.casefold(), k))
        global_override = overrides.get(("*", k))
        if specific_override:
            return specific_override, "MANUAL_OVERRIDE_SPECIFIC"
        if global_override:
            return global_override, "MANUAL_OVERRIDE_GLOBAL"

        specific = learned_specific.get((virus.casefold(), k))
        if specific:
            return specific, "TRUSTED_VIRUS_SPECIFIC_EXACT"

        global_target = learned_global.get(k)
        if global_target and k not in global_conflicts:
            return global_target, "TRUSTED_GLOBAL_EXACT"

        return "", "UNRESOLVED"

    out_rows = []
    for idx, row in df.iterrows():
        virus = norm_text(row[c_virus])
        category = norm_text(row[c_category])
        descriptions = unique_preserve(split_multi(row[c_entity_desc]))
        mapped = []
        mapping_sources = []
        unmapped = []

        for desc in descriptions:
            target, source = map_entity(virus, desc)
            if target:
                mapped.append(target)
                mapping_sources.append(f"{desc}=>{target} [{source}]")
            else:
                unmapped.append(desc)

        mapped_unique = unique_preserve(mapped)

        if category in NO_CONTACT_CATEGORIES or not descriptions:
            reconciled = ""
            status = "UNRESOLVED_NO_CONTACT_EVIDENCE"
        elif category in GENUINE_MULTI_CATEGORIES:
            # Preserve as reviewable multi-protein context. If all entities map,
            # record them explicitly rather than collapsing them.
            reconciled = ";".join(mapped_unique) if mapped_unique and not unmapped else ""
            status = (
                "PRESERVED_GENUINE_MULTIPROTEIN_INTERFACE"
                if reconciled and len(mapped_unique) > 1
                else "MANUAL_REVIEW_MULTIPROTEIN_INTERFACE"
            )
        elif unmapped:
            reconciled = ""
            status = "MANUAL_REVIEW_UNMAPPED_ENTITY"
        elif len(mapped_unique) == 1:
            reconciled = mapped_unique[0]
            status = "AUTO_RESOLVED_EXACT"
        elif len(mapped_unique) > 1:
            reconciled = ";".join(mapped_unique)
            status = "MANUAL_REVIEW_MULTIPLE_CONTACT_TARGETS"
        else:
            reconciled = ""
            status = "MANUAL_REVIEW_UNRESOLVED"

        rec = dict(row)
        rec.update(
            {
                "taxonomy_reconciled_target": reconciled,
                "taxonomy_reconciliation_status": status,
                "taxonomy_mapping_evidence": " | ".join(mapping_sources),
                "taxonomy_unmapped_entity_descriptions": ";".join(unmapped),
            }
        )
        out_rows.append(rec)

    rdf = pd.DataFrame(out_rows)

    # Stable ordering for deterministic output.
    sort_cols = [x for x in [c_virus, c_pdb, c_instance, c_ligand] if x]
    if sort_cols:
        rdf = rdf.sort_values(sort_cols, kind="mergesort")

    rdf.to_csv(args.outdir / "classification_reconciliation_audit_taxonomy_pass.csv", index=False)

    # ------------------------------------------------------------------
    # 4) Review queue, ranked by impact.
    # ------------------------------------------------------------------
    review = rdf[
        rdf["taxonomy_reconciliation_status"].str.startswith("MANUAL_REVIEW", na=False)
    ].copy()

    queue_records = []
    if not review.empty:
        review["_descs"] = review[c_entity_desc].map(lambda x: ";".join(unique_preserve(split_multi(x))))
        for (virus, descs, status), grp in review.groupby(
            [c_virus, "_descs", "taxonomy_reconciliation_status"], dropna=False, sort=True
        ):
            queue_records.append(
                {
                    "virus_name": virus,
                    "entity_descriptions": descs,
                    "review_status": status,
                    "occurrence_count": int(len(grp)),
                    "distinct_pdbs": int(grp[c_pdb].nunique()),
                    "pdb_examples": ";".join(sorted(grp[c_pdb].dropna().unique())[:20]),
                    "current_targets": ";".join(
                        unique_preserve(sorted(grp[c_current].dropna().tolist(), key=str.casefold))
                    ),
                    "folder_labels": ";".join(
                        unique_preserve(sorted(grp[c_folder].dropna().tolist(), key=str.casefold))
                    ),
                }
            )

    qdf = pd.DataFrame(queue_records)
    if not qdf.empty:
        qdf = qdf.sort_values(
            ["occurrence_count", "distinct_pdbs", "virus_name", "entity_descriptions"],
            ascending=[False, False, True, True],
            kind="mergesort",
        )
    qdf.to_csv(args.outdir / "taxonomy_manual_review_queue.csv", index=False)

    # A ready-to-edit override template from the unresolved queue.
    template_rows = []
    seen = set()
    for _, r in qdf.iterrows() if not qdf.empty else []:
        for desc in split_multi(r["entity_descriptions"]):
            key = (norm_text(r["virus_name"]), norm_key(desc))
            if key in seen:
                continue
            seen.add(key)
            template_rows.append(
                {
                    "virus_name": r["virus_name"],
                    "entity_description": desc,
                    "canonical_target": "",
                    "note": f"{r['occurrence_count']} unresolved occurrences in this queue bucket",
                }
            )
    pd.DataFrame(
        template_rows,
        columns=["virus_name", "entity_description", "canonical_target", "note"],
    ).to_csv(args.outdir / "taxonomy_manual_overrides_template.csv", index=False)

    # ------------------------------------------------------------------
    # 5) Audit-only Target Browser grouping simulation.
    # ------------------------------------------------------------------
    sim = rdf.copy()
    sim["simulated_target"] = sim["taxonomy_reconciled_target"]
    resolved = sim[sim["simulated_target"].astype(str).str.len() > 0].copy()

    diff_rows = []
    if not resolved.empty:
        group_cols = [c_virus, c_current, "simulated_target"]
        for keys, grp in resolved.groupby(group_cols, dropna=False, sort=True):
            virus, current_target, simulated_target = keys
            if norm_key(current_target) == norm_key(simulated_target):
                change = "UNCHANGED"
            else:
                change = "MERGE_OR_RELABEL"
            diff_rows.append(
                {
                    "virus_name": virus,
                    "current_target_label": current_target,
                    "simulated_target_label": simulated_target,
                    "change_type": change,
                    "ligand_occurrence_count": int(len(grp)),
                    "structure_count": int(grp[c_pdb].nunique()),
                    "affected_pdb_ids": ";".join(sorted(grp[c_pdb].dropna().unique())),
                }
            )

    ddf = pd.DataFrame(diff_rows)
    if not ddf.empty:
        ddf = ddf.sort_values(
            ["change_type", "ligand_occurrence_count", "virus_name", "current_target_label"],
            ascending=[True, False, True, True],
            kind="mergesort",
        )
    ddf.to_csv(args.outdir / "target_browser_grouping_taxonomy_pass.csv", index=False)

    # ------------------------------------------------------------------
    # 6) Summary.
    # ------------------------------------------------------------------
    status_counts = Counter(rdf["taxonomy_reconciliation_status"].tolist())
    summary = {
        "input_occurrences": int(len(rdf)),
        "distinct_structures": int(rdf[c_pdb].nunique()),
        "exact_taxonomy_candidates": int(
            (candidates["suggested_canonical_target"].astype(str).str.len() > 0).sum()
        ),
        "unresolved_taxonomy_candidates": int(
            (candidates["suggested_canonical_target"].astype(str).str.len() == 0).sum()
        ),
        "status_counts": dict(sorted(status_counts.items())),
        "auto_resolved_occurrences": int(
            (rdf["taxonomy_reconciliation_status"] == "AUTO_RESOLVED_EXACT").sum()
        ),
        "manual_review_occurrences": int(
            rdf["taxonomy_reconciliation_status"].str.startswith("MANUAL_REVIEW", na=False).sum()
        ),
        "preserved_genuine_multiprotein_occurrences": int(
            (rdf["taxonomy_reconciliation_status"] == "PRESERVED_GENUINE_MULTIPROTEIN_INTERFACE").sum()
        ),
        "overrides_loaded": int(len(overrides)),
        "mapping_policy": {
            "fuzzy_matching": False,
            "virus_specific_exact_trusted_mapping": True,
            "global_exact_mapping_only_if_unambiguous_across_trusted_evidence": True,
            "manual_override_support": True,
            "genuine_multiprotein_interfaces_auto_collapsed": False,
            "production_data_modified": False,
        },
    }
    with open(args.outdir / "taxonomy_reconciliation_summary.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2, sort_keys=True)
        fh.write("\n")

    print(json.dumps(summary, indent=2, sort_keys=True))
    print(f"\nWrote audit-only outputs to: {args.outdir}")


if __name__ == "__main__":
    main()
