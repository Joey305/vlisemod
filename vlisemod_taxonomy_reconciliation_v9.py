#!/usr/bin/env python3
"""
V-LiSEMOD taxonomy reconciliation — Pass 9
Consolidate Passes 3–8 into one authoritative occurrence-level target table
and simulate the corrected Target Browser (READ ONLY).

Purpose
-------
Create the single audit artifact that production can eventually consume.

Precedence
----------
Base: Pass-5 full 7,355-occurrence table (contains Pass-3 and Pass-4 columns).

1. Pass-4 source-organism decision establishes YES / NO / REVIEW.
2. Pass-5 unambiguous polyprotein mature-product resolutions override Pass-4.
3. Pass-6.2 exact source-conflict adjudications override the relevant 87 rows.
4. Pass-8 exact final polyprotein adjudications override the final 40 rows.

Nothing is written to the production database or website.

Inputs
------
--pass5-csv
    taxonomy_polyprotein_domain_qc_pass5.csv   [all 7,355 occurrences]

--pass6-2-csv
    taxonomy_source_conflict_adjudication_pass6_2.csv

--pass8-csv
    taxonomy_polyprotein_final_adjudication_pass8.csv

--outdir
    output directory

Optional
--------
--key-cols
    Comma-separated occurrence key columns. Normally unnecessary: the script
    auto-detects a key and verifies that it matches the 87 Pass-6.2 rows and
    40 Pass-8 rows exactly, with no accidental extra base-row matches.

Outputs
-------
canonical_target_occurrences_pass9.csv
canonical_target_browser_groups_pass9.csv
target_browser_grouping_diff_pass9.csv
canonical_target_review_queue_pass9.csv
canonical_target_excluded_contexts_pass9.csv
canonical_target_authority_summary.json

Design rules
------------
- Only final_target_browser_eligible == YES may create a viral Target Browser card.
- NO rows remain in provenance but do not create mature viral-target cards.
- REVIEW rows remain visible to curation/audit but are not production-authoritative.
- Genuine viral interfaces are preserved as explicit canonical interface targets.
- Folder-derived combined labels are retained only as old/provenance fields.
- Scientific scores are never changed.
"""

from __future__ import annotations

import argparse
import itertools
import json
import re
from collections import Counter
from pathlib import Path
from typing import Iterable

import pandas as pd


def txt(v) -> str:
    if v is None or pd.isna(v):
        return ""
    return re.sub(r"\s+", " ", str(v).strip())


def norm(v) -> str:
    s = txt(v).casefold().replace("_", " ")
    s = re.sub(r"[^a-z0-9]+", " ", s)
    return re.sub(r"\s+", " ", s).strip()


def resolve_col(df: pd.DataFrame, candidates: list[str], required: bool = True):
    lower = {c.casefold(): c for c in df.columns}
    for c in candidates:
        if c.casefold() in lower:
            return lower[c.casefold()]
    if required:
        raise KeyError(
            f"Missing required column; expected one of {candidates}. "
            f"Available columns: {list(df.columns)}"
        )
    return None


def first_present(df: pd.DataFrame, candidates: list[str]):
    return resolve_col(df, candidates, required=False)


def series_or_blank(df: pd.DataFrame, col: str | None) -> pd.Series:
    if col:
        return df[col].fillna("").astype(str)
    return pd.Series([""] * len(df), index=df.index, dtype=str)


def make_keys(df: pd.DataFrame, cols: list[str]) -> pd.Series:
    parts = []
    for c in cols:
        parts.append(df[c].fillna("").astype(str).map(txt))
    if not parts:
        raise ValueError("No key columns.")
    out = parts[0]
    for p in parts[1:]:
        out = out + "\x1f" + p
    return out


def exact_subset_match(base: pd.DataFrame, sub: pd.DataFrame, cols: list[str]) -> bool:
    if any(c not in base.columns or c not in sub.columns for c in cols):
        return False
    b = Counter(make_keys(base, cols))
    s = Counter(make_keys(sub, cols))
    # Every override-key multiplicity must match the base exactly. This prevents
    # a PDB-level rule from accidentally touching non-override occurrences.
    return all(b.get(k, 0) == n for k, n in s.items())


def choose_key_cols(
    base: pd.DataFrame,
    p62: pd.DataFrame,
    p8: pd.DataFrame,
    pdb_col: str,
    explicit: str | None,
) -> list[str]:
    if explicit:
        cols = [x.strip() for x in explicit.split(",") if x.strip()]
        missing = [
            c for c in cols
            if c not in base.columns or c not in p62.columns or c not in p8.columns
        ]
        if missing:
            raise KeyError(f"--key-cols missing from one or more inputs: {missing}")
        if not exact_subset_match(base, p62, cols):
            raise RuntimeError("--key-cols do not exactly identify the Pass-6.2 subset.")
        if not exact_subset_match(base, p8, cols):
            raise RuntimeError("--key-cols do not exactly identify the Pass-8 subset.")
        return cols

    common = set(base.columns) & set(p62.columns) & set(p8.columns)

    priority_names = [
        "occurrence_id",
        "ligand_occurrence_id",
        "ligand_instance",
        "ligand_instance_id",
        "row_index",
        "ligand_id",
        "ligand",
        "ligand_code",
        "ligand_chain",
        "ligand_auth_asym_id",
        "ligand_resseq",
        "ligand_residue_number",
        "ligand_auth_seq_id",
        "ligand_label_seq_id",
        "model_id",
        "altloc",
        "contacting_protein_chains",
        "target_chains",
        "stage09_target_chains",
        "contacting_entity_ids",
        "contact_entity_ids",
        "entity_ids",
        "contacting_entity_descriptions",
        "entity_descriptions",
        "current_stage14_protein_type",
        "current_target_label",
        "protein_type",
        "pass5_status",
    ]

    by_lower = {c.casefold(): c for c in common}
    candidates = []
    for name in priority_names:
        c = by_lower.get(name.casefold())
        if c and c != pdb_col and c not in candidates:
            candidates.append(c)

    # Fast-path likely keys first.
    likely = []
    for name in [
        "occurrence_id", "ligand_occurrence_id",
        "ligand_instance", "ligand_instance_id", "row_index"
    ]:
        c = by_lower.get(name.casefold())
        if c:
            likely.append([pdb_col, c])

    for cols in likely:
        if exact_subset_match(base, p62, cols) and exact_subset_match(base, p8, cols):
            return cols

    # Greedily test informative combinations. We cap at 5 added columns to
    # avoid excessive search while still handling legacy audit schemas.
    max_r = min(5, len(candidates))
    for r in range(1, max_r + 1):
        for combo in itertools.combinations(candidates, r):
            cols = [pdb_col, *combo]
            if exact_subset_match(base, p62, cols) and exact_subset_match(base, p8, cols):
                return cols

    available = [pdb_col, *candidates]
    raise RuntimeError(
        "Could not auto-detect an occurrence key that matches both override "
        "subsets exactly. Re-run with --key-cols. Candidate columns were:\n"
        + ", ".join(available)
    )


def ensure_consistent_override(
    df: pd.DataFrame,
    key_col: str,
    value_cols: list[str],
    label: str,
) -> pd.DataFrame:
    """
    Multiple rows may share a key only if all override fields are identical.
    """
    keep = [key_col] + value_cols
    tmp = df[keep].copy().fillna("")
    for k, grp in tmp.groupby(key_col, sort=False):
        unique = grp[value_cols].drop_duplicates()
        if len(unique) > 1:
            raise RuntimeError(
                f"{label} has inconsistent override values for occurrence key {k!r}."
            )
    return tmp.drop_duplicates(subset=[key_col], keep="first")


def canonical_id_fallback(virus: str, target_id: str, target_name: str) -> str:
    if txt(target_id):
        return txt(target_id)
    if not txt(target_name):
        return ""
    slug = norm(target_name).replace(" ", "_")
    vslug = norm(virus).replace(" ", "_")
    return f"{vslug}_{slug}" if vslug else slug


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pass5-csv", type=Path, required=True)
    ap.add_argument("--pass6-2-csv", type=Path, required=True)
    ap.add_argument("--pass8-csv", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument(
        "--key-cols",
        default=None,
        help="Optional comma-separated occurrence key columns.",
    )
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    base = pd.read_csv(args.pass5_csv, dtype=str, low_memory=False).fillna("")
    p62 = pd.read_csv(args.pass6_2_csv, dtype=str, low_memory=False).fillna("")
    p8 = pd.read_csv(args.pass8_csv, dtype=str, low_memory=False).fillna("")

    # Core base columns.
    c_pdb = resolve_col(base, ["pdb_id", "pdb"])
    c_virus = resolve_col(base, ["virus_name", "virus"])
    c_ligand = first_present(base, ["ligand_id", "ligand", "ligand_code"])
    c_instance = first_present(
        base, ["occurrence_id", "ligand_occurrence_id", "ligand_instance", "ligand_instance_id"]
    )
    c_current = first_present(
        base, ["current_stage14_protein_type", "current_target_label", "protein_type"]
    )

    c_p3_id = first_present(base, ["pass3_canonical_target_id"])
    c_p3_name = first_present(base, ["pass3_canonical_target_name"])
    c_p3_family = first_present(base, ["pass3_target_family"])
    c_p3_role = first_present(base, ["pass3_entity_role"])
    c_p3_status = first_present(base, ["pass3_status"])

    c_p4_status = resolve_col(base, ["pass4_status"])
    c_p4_eligible = resolve_col(base, ["pass4_target_browser_eligible"])
    c_p4_source = first_present(base, ["pass4_source_aggregate"])

    c_p5_status = resolve_col(base, ["pass5_status"])
    c_p5_id = first_present(base, ["pass5_canonical_target_id"])
    c_p5_name = first_present(base, ["pass5_canonical_target_name"])
    c_p5_family = first_present(base, ["pass5_target_family"])

    # Verify input cardinalities before doing anything.
    if len(base) == 0:
        raise RuntimeError("Pass-5 base table is empty.")
    if len(p62) == 0:
        raise RuntimeError("Pass-6.2 adjudication table is empty.")
    if len(p8) == 0:
        raise RuntimeError("Pass-8 adjudication table is empty.")

    key_cols = choose_key_cols(base, p62, p8, c_pdb, args.key_cols)

    base = base.copy()
    p62 = p62.copy()
    p8 = p8.copy()

    base["pass9_occurrence_key"] = make_keys(base, key_cols)
    p62["pass9_occurrence_key"] = make_keys(p62, key_cols)
    p8["pass9_occurrence_key"] = make_keys(p8, key_cols)

    # ------------------------------------------------------------------
    # Start from Pass-4 authority.
    # ------------------------------------------------------------------
    out = base.copy()

    out["candidate_canonical_target_id"] = series_or_blank(out, c_p3_id)
    out["candidate_canonical_target_name"] = series_or_blank(out, c_p3_name)
    out["candidate_target_family"] = series_or_blank(out, c_p3_family)
    out["candidate_entity_role"] = series_or_blank(out, c_p3_role)

    out["final_canonical_target_id"] = ""
    out["final_canonical_target_name"] = ""
    out["final_target_family"] = ""
    out["final_entity_role"] = ""
    out["final_target_browser_eligible"] = "REVIEW"
    out["final_decision"] = "REVIEW_PENDING"
    out["final_authority_basis"] = "PASS4_REVIEW_PENDING"
    out["final_authority_note"] = ""

    # Pass 4 YES.
    yes4 = out[c_p4_eligible].eq("YES")
    out.loc[yes4, "final_canonical_target_id"] = series_or_blank(out, c_p3_id)[yes4]
    out.loc[yes4, "final_canonical_target_name"] = series_or_blank(out, c_p3_name)[yes4]
    out.loc[yes4, "final_target_family"] = series_or_blank(out, c_p3_family)[yes4]
    out.loc[yes4, "final_entity_role"] = series_or_blank(out, c_p3_role)[yes4].replace(
        {"VIRAL_TARGET_PROVISIONAL": "VIRAL_TARGET"}
    )
    out.loc[yes4, "final_target_browser_eligible"] = "YES"
    out.loc[yes4, "final_decision"] = "AUTHORITATIVE_VIRAL_TARGET"
    out.loc[yes4, "final_authority_basis"] = "PASS4_SOURCE_CONFIRMED"

    # Pass 4 NO.
    no4 = out[c_p4_eligible].eq("NO")
    out.loc[no4, "final_target_browser_eligible"] = "NO"
    out.loc[no4, "final_decision"] = "EXCLUDE_CONTEXT_FROM_MATURE_VIRAL_TARGET_BROWSER"
    out.loc[no4, "final_authority_basis"] = "PASS4_NONVIRAL_CONTEXT"
    out.loc[no4, "final_entity_role"] = series_or_blank(out, c_p3_role)[no4]

    # ------------------------------------------------------------------
    # Pass 5: 133 unambiguous mature-product resolutions.
    # ------------------------------------------------------------------
    p5resolved = out[c_p5_status].eq("RESOLVED_POLYPROTEIN_MATURE_PRODUCT")
    out.loc[p5resolved, "final_canonical_target_id"] = series_or_blank(out, c_p5_id)[p5resolved]
    out.loc[p5resolved, "final_canonical_target_name"] = series_or_blank(out, c_p5_name)[p5resolved]
    out.loc[p5resolved, "final_target_family"] = series_or_blank(out, c_p5_family)[p5resolved]
    out.loc[p5resolved, "final_entity_role"] = "VIRAL_TARGET"
    out.loc[p5resolved, "final_target_browser_eligible"] = "YES"
    out.loc[p5resolved, "final_decision"] = "AUTHORITATIVE_VIRAL_TARGET"
    out.loc[p5resolved, "final_authority_basis"] = "PASS5_UNAMBIGUOUS_POLYPROTEIN_ALIGNMENT"

    # ------------------------------------------------------------------
    # Pass 6.2 exact adjudication overrides.
    # ------------------------------------------------------------------
    p62_dec = resolve_col(p62, ["pass6_2_decision"])
    p62_el = resolve_col(p62, ["pass6_2_target_browser_eligible"])
    p62_id = resolve_col(p62, ["pass6_2_canonical_target_id"])
    p62_name = resolve_col(p62, ["pass6_2_canonical_target_name"])
    p62_family = resolve_col(p62, ["pass6_2_target_family"])
    p62_role = resolve_col(p62, ["pass6_2_entity_role"])
    p62_reason = resolve_col(p62, ["pass6_2_reason"])

    ov62_cols = [p62_dec, p62_el, p62_id, p62_name, p62_family, p62_role, p62_reason]
    ov62 = ensure_consistent_override(
        p62, "pass9_occurrence_key", ov62_cols, "Pass 6.2"
    ).set_index("pass9_occurrence_key")

    for idx in out.index[out["pass9_occurrence_key"].isin(ov62.index)]:
        k = out.at[idx, "pass9_occurrence_key"]
        r = ov62.loc[k]
        decision = txt(r[p62_dec])

        if decision == "RETAIN_VIRAL_TARGET":
            out.at[idx, "final_canonical_target_id"] = txt(r[p62_id])
            out.at[idx, "final_canonical_target_name"] = txt(r[p62_name])
            out.at[idx, "final_target_family"] = txt(r[p62_family])
            out.at[idx, "final_entity_role"] = txt(r[p62_role])
            out.at[idx, "final_target_browser_eligible"] = "YES"
            out.at[idx, "final_decision"] = "AUTHORITATIVE_VIRAL_TARGET"
            out.at[idx, "final_authority_basis"] = "PASS6_2_EXACT_SOURCE_CONFLICT_ADJUDICATION"
            out.at[idx, "final_authority_note"] = txt(r[p62_reason])

        elif decision == "EXCLUDE_NONVIRAL_CONTACT_CONTEXT":
            out.at[idx, "final_canonical_target_id"] = ""
            out.at[idx, "final_canonical_target_name"] = ""
            out.at[idx, "final_target_family"] = ""
            out.at[idx, "final_entity_role"] = txt(r[p62_role])
            out.at[idx, "final_target_browser_eligible"] = "NO"
            out.at[idx, "final_decision"] = "EXCLUDE_CONTEXT_FROM_MATURE_VIRAL_TARGET_BROWSER"
            out.at[idx, "final_authority_basis"] = "PASS6_2_EXACT_NONVIRAL_ADJUDICATION"
            out.at[idx, "final_authority_note"] = txt(r[p62_reason])
        else:
            raise RuntimeError(f"Unexpected Pass-6.2 decision: {decision!r}")

    # ------------------------------------------------------------------
    # Pass 8 exact final polyprotein adjudication overrides.
    # ------------------------------------------------------------------
    p8_dec = resolve_col(p8, ["pass8_decision"])
    p8_el = resolve_col(p8, ["pass8_target_browser_eligible"])
    p8_id = resolve_col(p8, ["pass8_canonical_target_id"])
    p8_name = resolve_col(p8, ["pass8_canonical_target_name"])
    p8_family = resolve_col(p8, ["pass8_target_family"])
    p8_role = resolve_col(p8, ["pass8_entity_role"])
    p8_reason = resolve_col(p8, ["pass8_reason"])

    ov8_cols = [p8_dec, p8_el, p8_id, p8_name, p8_family, p8_role, p8_reason]
    ov8 = ensure_consistent_override(
        p8, "pass9_occurrence_key", ov8_cols, "Pass 8"
    ).set_index("pass9_occurrence_key")

    for idx in out.index[out["pass9_occurrence_key"].isin(ov8.index)]:
        k = out.at[idx, "pass9_occurrence_key"]
        r = ov8.loc[k]
        decision = txt(r[p8_dec])

        if decision in {"RESOLVE_VIRAL_TARGET", "RESOLVE_VIRAL_INTERFACE"}:
            out.at[idx, "final_canonical_target_id"] = txt(r[p8_id])
            out.at[idx, "final_canonical_target_name"] = txt(r[p8_name])
            out.at[idx, "final_target_family"] = txt(r[p8_family])
            out.at[idx, "final_entity_role"] = txt(r[p8_role])
            out.at[idx, "final_target_browser_eligible"] = "YES"
            out.at[idx, "final_decision"] = (
                "AUTHORITATIVE_VIRAL_INTERFACE"
                if decision == "RESOLVE_VIRAL_INTERFACE"
                else "AUTHORITATIVE_VIRAL_TARGET"
            )
            out.at[idx, "final_authority_basis"] = "PASS8_EXACT_POLYPROTEIN_ADJUDICATION"
            out.at[idx, "final_authority_note"] = txt(r[p8_reason])

        elif decision == "EXCLUDE_CONTEXT_FROM_MATURE_VIRAL_TARGET_BROWSER":
            out.at[idx, "final_canonical_target_id"] = ""
            out.at[idx, "final_canonical_target_name"] = ""
            out.at[idx, "final_target_family"] = ""
            out.at[idx, "final_entity_role"] = txt(r[p8_role])
            out.at[idx, "final_target_browser_eligible"] = "NO"
            out.at[idx, "final_decision"] = "EXCLUDE_CONTEXT_FROM_MATURE_VIRAL_TARGET_BROWSER"
            out.at[idx, "final_authority_basis"] = "PASS8_EXACT_CONTEXT_ADJUDICATION"
            out.at[idx, "final_authority_note"] = txt(r[p8_reason])
        else:
            raise RuntimeError(f"Unexpected Pass-8 decision: {decision!r}")

    # Canonical ID fallback only for authoritative YES rows.
    for idx in out.index[out["final_target_browser_eligible"].eq("YES")]:
        out.at[idx, "final_canonical_target_id"] = canonical_id_fallback(
            txt(out.at[idx, c_virus]),
            txt(out.at[idx, "final_canonical_target_id"]),
            txt(out.at[idx, "final_canonical_target_name"]),
        )

    # ------------------------------------------------------------------
    # Validation.
    # ------------------------------------------------------------------
    eligibility_counts = Counter(out["final_target_browser_eligible"])
    if sum(eligibility_counts.values()) != len(out):
        raise RuntimeError("Final eligibility partition does not equal input rows.")

    yes = out["final_target_browser_eligible"].eq("YES")
    missing_target = out.loc[
        yes & (
            out["final_canonical_target_id"].map(txt).eq("")
            | out["final_canonical_target_name"].map(txt).eq("")
        )
    ]
    if len(missing_target):
        raise RuntimeError(
            f"{len(missing_target)} authoritative YES rows lack a canonical target ID/name."
        )

    # Override subset coverage must be exact.
    matched62 = int(out["pass9_occurrence_key"].isin(ov62.index).sum())
    matched8 = int(out["pass9_occurrence_key"].isin(ov8.index).sum())
    if matched62 != len(p62):
        raise RuntimeError(
            f"Pass-6.2 coverage mismatch: expected {len(p62)}, matched {matched62}."
        )
    if matched8 != len(p8):
        raise RuntimeError(
            f"Pass-8 coverage mismatch: expected {len(p8)}, matched {matched8}."
        )

    # A canonical target ID itself should never be the old delimiter-joined
    # provenance label. Explicit interface IDs use underscores, not commas/semicolons.
    bad_canonical_ids = out.loc[
        yes & out["final_canonical_target_id"].str.contains(r"[;,]", regex=True, na=False)
    ]

    # ------------------------------------------------------------------
    # Authoritative occurrence table.
    # ------------------------------------------------------------------
    occurrence_path = args.outdir / "canonical_target_occurrences_pass9.csv"
    out.to_csv(occurrence_path, index=False)

    # Review / excluded subsets.
    out[out["final_target_browser_eligible"].eq("REVIEW")].to_csv(
        args.outdir / "canonical_target_review_queue_pass9.csv",
        index=False,
    )
    out[out["final_target_browser_eligible"].eq("NO")].to_csv(
        args.outdir / "canonical_target_excluded_contexts_pass9.csv",
        index=False,
    )

    # ------------------------------------------------------------------
    # Corrected Target Browser card simulation: YES rows only.
    # ------------------------------------------------------------------
    yesdf = out[yes].copy()

    def examples(values: Iterable[str], n=20):
        vals = sorted({txt(x) for x in values if txt(x)})
        return ";".join(vals[:n])

    group_cols = [
        c_virus,
        "final_canonical_target_id",
        "final_canonical_target_name",
        "final_target_family",
        "final_entity_role",
    ]

    group_rows = []
    for keys, grp in yesdf.groupby(group_cols, dropna=False, sort=True):
        d = dict(zip(group_cols, keys if isinstance(keys, tuple) else (keys,)))
        d["occurrence_count"] = int(len(grp))
        d["distinct_pdbs"] = int(grp[c_pdb].nunique())
        d["pdb_examples"] = examples(grp[c_pdb])
        if c_ligand:
            d["distinct_ligands"] = int(grp[c_ligand].nunique())
            d["ligand_examples"] = examples(grp[c_ligand])
        if c_instance:
            d["distinct_ligand_instances"] = int(grp[c_instance].nunique())
        group_rows.append(d)

    groups = pd.DataFrame(group_rows)
    if not groups.empty:
        groups = groups.sort_values(
            ["occurrence_count", c_virus, "final_canonical_target_name"],
            ascending=[False, True, True],
            kind="mergesort",
        )
    groups.to_csv(
        args.outdir / "canonical_target_browser_groups_pass9.csv",
        index=False,
    )

    # ------------------------------------------------------------------
    # Old Stage-14 grouping -> final grouping mapping matrix.
    # ------------------------------------------------------------------
    old_target_series = series_or_blank(out, c_current)
    out["_pass9_old_target"] = old_target_series

    diff_cols = [
        c_virus,
        "_pass9_old_target",
        "final_target_browser_eligible",
        "final_canonical_target_id",
        "final_canonical_target_name",
        "final_target_family",
        "final_entity_role",
        "final_authority_basis",
    ]

    diff_rows = []
    for keys, grp in out.groupby(diff_cols, dropna=False, sort=True):
        d = dict(zip(diff_cols, keys if isinstance(keys, tuple) else (keys,)))
        d["occurrence_count"] = int(len(grp))
        d["distinct_pdbs"] = int(grp[c_pdb].nunique())
        d["pdb_examples"] = examples(grp[c_pdb])
        if c_ligand:
            d["distinct_ligands"] = int(grp[c_ligand].nunique())
        diff_rows.append(d)

    diff = pd.DataFrame(diff_rows)
    if not diff.empty:
        diff = diff.rename(columns={"_pass9_old_target": "current_stage14_target"})
        diff = diff.sort_values(
            ["occurrence_count", c_virus, "current_stage14_target"],
            ascending=[False, True, True],
            kind="mergesort",
        )
    diff.to_csv(
        args.outdir / "target_browser_grouping_diff_pass9.csv",
        index=False,
    )

    # ------------------------------------------------------------------
    # Summary / spot checks.
    # ------------------------------------------------------------------
    current_yes_groups = 0
    if c_current:
        current_yes_groups = int(
            yesdf[[c_virus, c_current]].drop_duplicates().shape[0]
        )

    final_groups = int(len(groups))

    current_compound_occ = 0
    if c_current:
        current_compound_occ = int(
            out[c_current].str.contains(r"[;,]", regex=True, na=False).sum()
        )

    authority_counts = Counter(out["final_authority_basis"])
    decision_counts = Counter(out["final_decision"])

    spot_checks = {}
    for pdb in ["2O4K"]:
        sub = out[out[c_pdb].astype(str).str.upper().eq(pdb)]
        if not sub.empty:
            spot_checks[pdb] = {
                "occurrences": int(len(sub)),
                "current_stage14_targets": sorted(
                    {txt(x) for x in series_or_blank(sub, c_current) if txt(x)}
                ),
                "final_eligibilities": sorted(
                    {txt(x) for x in sub["final_target_browser_eligible"] if txt(x)}
                ),
                "final_targets": sorted(
                    {txt(x) for x in sub["final_canonical_target_name"] if txt(x)}
                ),
                "final_target_ids": sorted(
                    {txt(x) for x in sub["final_canonical_target_id"] if txt(x)}
                ),
            }

    summary = {
        "input_occurrences": int(len(out)),
        "distinct_structures": int(out[c_pdb].nunique()),
        "occurrence_key_columns": key_cols,
        "override_coverage": {
            "pass6_2_rows": int(len(p62)),
            "pass6_2_rows_matched": matched62,
            "pass8_rows": int(len(p8)),
            "pass8_rows_matched": matched8,
        },
        "final_target_browser_eligibility_counts": dict(
            sorted(eligibility_counts.items())
        ),
        "final_decision_counts": dict(sorted(decision_counts.items())),
        "final_authority_basis_counts": dict(sorted(authority_counts.items())),
        "authoritative_yes_occurrences": int(yes.sum()),
        "excluded_no_occurrences": int(
            out["final_target_browser_eligible"].eq("NO").sum()
        ),
        "review_occurrences": int(
            out["final_target_browser_eligible"].eq("REVIEW").sum()
        ),
        "current_stage14_group_count_on_final_yes_rows": current_yes_groups,
        "final_canonical_target_browser_group_count": final_groups,
        "current_compound_label_occurrences": current_compound_occ,
        "authoritative_canonical_ids_with_semicolon_or_comma": int(
            len(bad_canonical_ids)
        ),
        "spot_checks": spot_checks,
        "validation": {
            "all_rows_partitioned": True,
            "all_authoritative_yes_rows_have_target_id_and_name": True,
            "pass6_2_exact_subset_match": matched62 == len(p62),
            "pass8_exact_subset_match": matched8 == len(p8),
        },
        "production_data_modified": False,
        "scientific_scores_modified": False,
    }

    (args.outdir / "canonical_target_authority_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    print(json.dumps(summary, indent=2, sort_keys=True))
    print(f"\nWrote audit-only Pass-9 outputs to: {args.outdir}")


if __name__ == "__main__":
    main()
