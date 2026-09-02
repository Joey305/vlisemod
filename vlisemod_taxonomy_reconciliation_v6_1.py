#!/usr/bin/env python3
"""
V-LiSEMOD taxonomy reconciliation — Pass 6
Source-conflict forensics (READ ONLY)

Purpose
-------
Investigate the Pass-4 rows with:

    CONFLICT_TARGET_ASSIGNED_BUT_SOURCE_NONVIRAL

These are potentially important because a canonical viral target was assigned
earlier, while the contacted entity's mmCIF gene/natural source is nonviral.

Pass 6 does NOT modify production data. It determines why the conflict exists
and places each conflict into a conservative evidence bucket.

Inputs
------
--pass5-csv
    taxonomy_polyprotein_domain_qc_pass5.csv
--seed-csv
    reviewed_taxonomy_seed_v3.csv

Outputs
-------
taxonomy_source_conflict_forensics_pass6.csv
taxonomy_source_conflict_signatures_pass6.csv
taxonomy_source_conflict_summary.json
taxonomy_source_conflict_review_queue_pass6.csv

Important:
- No fuzzy matching.
- No folder label is used as deciding evidence.
- Exact reviewed entity-role mappings from the Pass-3 seed are allowed.
- mmCIF _struct_ref metadata are inspected as independent reference evidence.
- Nothing is automatically written back to V-LiSEMOD.
"""

from __future__ import annotations

import argparse
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import pandas as pd

try:
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
except Exception as exc:
    raise SystemExit(
        "Biopython is required. Use the same V-LiSEMOD environment as the "
        "rebuild pipeline. Import error: " + repr(exc)
    )


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

NULLISH = {"", ".", "?", "none", "null", "na", "n/a", "unknown"}


def txt(v: Any) -> str:
    if v is None or pd.isna(v):
        return ""
    return re.sub(r"\s+", " ", str(v).strip())


def clean(v: Any) -> str:
    s = txt(v).strip("'\"")
    if s.casefold() in NULLISH:
        return ""
    return s


def norm(v: Any) -> str:
    s = clean(v).casefold().replace("_", " ")
    s = re.sub(r"[^a-z0-9]+", " ", s)
    return re.sub(r"\s+", " ", s).strip()


def as_list(v: Any) -> list[str]:
    if v is None:
        return []
    if isinstance(v, (list, tuple)):
        return [clean(x) for x in v]
    return [clean(v)]


def split_ids(v: Any) -> list[str]:
    s = txt(v)
    if not s:
        return []
    vals = [x.strip() for x in re.split(r"[;,|]", s) if x.strip()]
    out, seen = [], set()
    for x in vals:
        k = x.casefold()
        if k not in seen:
            seen.add(k)
            out.append(x)
    return out


def resolve_col(df: pd.DataFrame, candidates: list[str], required: bool = True) -> str | None:
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


def value_at(vals: list[str], i: int) -> str:
    return vals[i] if i < len(vals) else ""


# ---------------------------------------------------------------------------
# Reviewed seed
# ---------------------------------------------------------------------------

def load_seed(path: Path) -> dict[tuple[str, str], dict[str, str]]:
    df = pd.read_csv(path, dtype=str).fillna("")
    out = {}
    for _, r in df.iterrows():
        virus = clean(r.get("virus_name"))
        desc = clean(r.get("entity_description"))
        if not virus or not desc:
            continue
        out[(virus.casefold(), norm(desc))] = {
            "canonical_target_id": clean(r.get("canonical_target_id")),
            "canonical_target_name": clean(r.get("canonical_target_name")),
            "target_family": clean(r.get("target_family")),
            "entity_role": clean(r.get("entity_role")),
            "note": clean(r.get("note")),
        }
    return out


NONVIRAL_SEED_ROLES = {
    "HOST_PARTNER",
    "ENGINEERED_COMPONENT",
    "EXTERNAL_BINDER",
    "NONPOLYMER_ARTIFACT",
}

VIRAL_SEED_ROLE_PREFIX = "VIRAL_TARGET"


# ---------------------------------------------------------------------------
# Virus-reference evidence
# ---------------------------------------------------------------------------

VIRUS_TEXT_PATTERNS = {
    "sars_cov_2": [
        r"\bsevere acute respiratory syndrome coronavirus 2\b",
        r"\bsars cov 2\b",
        r"\bsars2\b",
        r"\bsars 2\b",
        r"\b2019 ncov\b",
    ],
    "hiv_1": [
        r"\bhuman immunodeficiency virus 1\b",
        r"\bhuman immunodeficiency virus type 1\b",
        r"\bhiv 1\b",
        # UniProt/PDB reference codes frequently look like HIV1... or HV1H2.
        r"\bhiv1[a-z0-9]*\b",
        r"\bhv1[a-z0-9]*\b",
    ],
    "hpv_18": [
        r"\bhuman papillomavirus type 18\b",
        r"\bhuman papillomavirus 18\b",
        r"\bhpv 18\b",
        r"\bhpv18[a-z0-9]*\b",
    ],
}


def virus_key(virus: str) -> str:
    s = norm(virus)
    if "sars" in s and ("cov" in s or "coronavirus" in s) and "2" in s:
        return "sars_cov_2"
    if ("hiv" in s and "1" in s) or "human immunodeficiency virus 1" in s:
        return "hiv_1"
    if ("hpv" in s and "18" in s) or ("human papillomavirus" in s and "18" in s):
        return "hpv_18"
    return re.sub(r"\s+", "_", s)


def text_supports_expected_virus(virus: str, value: str) -> bool:
    s = norm(value)
    if not s:
        return False
    for pat in VIRUS_TEXT_PATTERNS.get(virus_key(virus), []):
        if re.search(pat, s):
            return True
    return False


# ---------------------------------------------------------------------------
# mmCIF reference metadata
# ---------------------------------------------------------------------------

def parse_entity_desc(cif: dict[str, Any]) -> dict[str, str]:
    eids = as_list(cif.get("_entity.id"))
    descs = as_list(cif.get("_entity.pdbx_description"))
    return {eid: value_at(descs, i) for i, eid in enumerate(eids) if eid}


def parse_struct_refs(cif: dict[str, Any]) -> dict[str, dict[str, str]]:
    ids = as_list(cif.get("_struct_ref.id"))
    entity_ids = as_list(cif.get("_struct_ref.entity_id"))
    accessions = as_list(cif.get("_struct_ref.pdbx_db_accession"))
    codes = as_list(cif.get("_struct_ref.db_code"))
    names = as_list(cif.get("_struct_ref.db_name"))

    out = {}
    for i, rid in enumerate(ids):
        if not rid:
            continue
        out[rid] = {
            "entity_id": value_at(entity_ids, i),
            "accession": value_at(accessions, i),
            "db_code": value_at(codes, i),
            "db_name": value_at(names, i),
        }
    return out


def refs_for_entity(cif: dict[str, Any], entity_id: str) -> list[dict[str, str]]:
    refs = parse_struct_refs(cif)
    rows = []
    for rid, meta in refs.items():
        if clean(meta.get("entity_id")) == clean(entity_id):
            x = dict(meta)
            x["ref_id"] = rid
            rows.append(x)
    return rows


# ---------------------------------------------------------------------------
# Conflict evidence classification
# ---------------------------------------------------------------------------

def parse_pass4_source_details(value: str) -> list[dict[str, Any]]:
    s = clean(value)
    if not s:
        return []
    try:
        data = json.loads(s)
        return data if isinstance(data, list) else []
    except Exception:
        return []


def classify_entity(
    virus: str,
    entity_id: str,
    description: str,
    refs: list[dict[str, str]],
    seed: dict[tuple[str, str], dict[str, str]],
    source_detail: dict[str, Any] | None,
) -> dict[str, Any]:
    seed_hit = seed.get((virus.casefold(), norm(description)))

    ref_texts = []
    viral_ref_support = False
    for r in refs:
        t = " ".join(
            clean(r.get(x))
            for x in ["db_name", "db_code", "accession"]
            if clean(r.get(x))
        )
        if t:
            ref_texts.append(t)
            if text_supports_expected_virus(virus, t):
                viral_ref_support = True

    desc_viral_support = text_supports_expected_virus(virus, description)

    source_detail = source_detail or {}
    gene_sources = source_detail.get("gene_source", []) or []
    natural_sources = source_detail.get("natural_source", []) or []
    synthetic_sources = source_detail.get("synthetic_source", []) or []
    source_verdict = clean(source_detail.get("source_verdict"))

    seed_role = clean(seed_hit.get("entity_role")) if seed_hit else ""
    seed_target = clean(seed_hit.get("canonical_target_name")) if seed_hit else ""

    if seed_role in NONVIRAL_SEED_ROLES:
        evidence_class = "REVIEWED_NONVIRAL_IDENTITY"
        recommendation = "EXCLUDE_FROM_VIRAL_TARGET_BROWSER"
        confidence = "HIGH"

    elif seed_role.startswith(VIRAL_SEED_ROLE_PREFIX):
        # Exact reviewed viral identity conflicts with nonviral source metadata.
        evidence_class = "CURATED_VIRAL_IDENTITY_SOURCE_METADATA_CONFLICT"
        recommendation = "RETAIN_TARGET_CANDIDATE_REVIEW_SOURCE_METADATA"
        confidence = "HIGH"

    elif viral_ref_support and desc_viral_support:
        evidence_class = "VIRAL_DESCRIPTION_AND_REFERENCE_SOURCE_CONFLICT"
        recommendation = "RETAIN_TARGET_CANDIDATE_REVIEW_SOURCE_METADATA"
        confidence = "HIGH"

    elif viral_ref_support:
        evidence_class = "VIRAL_DATABASE_REFERENCE_SOURCE_CONFLICT"
        recommendation = "RETAIN_TARGET_CANDIDATE_REVIEW_SOURCE_METADATA"
        confidence = "MEDIUM"

    elif desc_viral_support:
        evidence_class = "VIRAL_DESCRIPTION_ONLY_SOURCE_CONFLICT"
        recommendation = "MANUAL_REVIEW"
        confidence = "MEDIUM"

    else:
        evidence_class = "NONVIRAL_SOURCE_WITHOUT_INDEPENDENT_VIRAL_SUPPORT"
        recommendation = "EXCLUDE_CANDIDATE_REQUIRES_REVIEW"
        confidence = "MEDIUM"

    return {
        "entity_id": entity_id,
        "entity_description": description,
        "seed_role": seed_role,
        "seed_target": seed_target,
        "gene_sources": gene_sources,
        "natural_sources": natural_sources,
        "synthetic_sources": synthetic_sources,
        "source_verdict": source_verdict,
        "struct_refs": refs,
        "viral_reference_support": viral_ref_support,
        "viral_description_support": desc_viral_support,
        "evidence_class": evidence_class,
        "recommendation": recommendation,
        "confidence": confidence,
    }


def aggregate_entities(entity_results: list[dict[str, Any]]) -> tuple[str, str, str]:
    if not entity_results:
        return (
            "NO_ENTITY_EVIDENCE",
            "MANUAL_REVIEW",
            "LOW",
        )

    classes = {x["evidence_class"] for x in entity_results}
    recs = {x["recommendation"] for x in entity_results}

    if classes == {"REVIEWED_NONVIRAL_IDENTITY"}:
        return (
            "CONFIRMED_NONVIRAL_CONTACT_CONTEXT",
            "EXCLUDE_FROM_VIRAL_TARGET_BROWSER",
            "HIGH",
        )

    if classes <= {
        "CURATED_VIRAL_IDENTITY_SOURCE_METADATA_CONFLICT",
        "VIRAL_DESCRIPTION_AND_REFERENCE_SOURCE_CONFLICT",
        "VIRAL_DATABASE_REFERENCE_SOURCE_CONFLICT",
    }:
        return (
            "LIKELY_VIRAL_TARGET_WITH_SOURCE_METADATA_CONFLICT",
            "RETAIN_TARGET_CANDIDATE_REVIEW_SOURCE_METADATA",
            "HIGH" if all(x["confidence"] == "HIGH" for x in entity_results) else "MEDIUM",
        )

    if "REVIEWED_NONVIRAL_IDENTITY" in classes and any(
        c.startswith("VIRAL_") or c.startswith("CURATED_VIRAL_")
        for c in classes
    ):
        return (
            "MIXED_VIRAL_NONVIRAL_IDENTITY_CONFLICT",
            "MANUAL_INTERFACE_REVIEW",
            "HIGH",
        )

    if recs == {"EXCLUDE_CANDIDATE_REQUIRES_REVIEW"}:
        return (
            "PROBABLE_NONVIRAL_CONTACT_CONTEXT",
            "EXCLUDE_CANDIDATE_REQUIRES_REVIEW",
            "MEDIUM",
        )

    return (
        "AMBIGUOUS_SOURCE_CONFLICT",
        "MANUAL_REVIEW",
        "LOW",
    )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--pass5-csv",
        type=Path,
        required=True,
        help="taxonomy_polyprotein_domain_qc_pass5.csv",
    )
    ap.add_argument(
        "--seed-csv",
        type=Path,
        required=True,
        help="reviewed_taxonomy_seed_v3.csv",
    )
    ap.add_argument("--outdir", type=Path, required=True)
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.pass5_csv, dtype=str, low_memory=False).fillna("")
    seed = load_seed(args.seed_csv)

    c_pdb = resolve_col(df, ["pdb_id", "pdb"])
    c_virus = resolve_col(df, ["virus_name", "virus"])
    c_pass4 = resolve_col(df, ["pass4_status"])
    c_cif = resolve_col(df, ["pass4_selected_cif_path"])
    c_entity_ids = resolve_col(
        df, ["contacting_entity_ids", "contact_entity_ids", "entity_ids"]
    )
    c_entity_descs = resolve_col(
        df, ["contacting_entity_descriptions", "entity_descriptions"], required=False
    )
    c_source_details = resolve_col(
        df, ["pass4_contact_entity_source_details"], required=False
    )
    c_target_id = resolve_col(
        df, ["pass3_canonical_target_id"], required=False
    )
    c_target_name = resolve_col(
        df, ["pass3_canonical_target_name"], required=False
    )
    c_target_family = resolve_col(
        df, ["pass3_target_family"], required=False
    )
    c_pass3_status = resolve_col(df, ["pass3_status"], required=False)

    conflicts = df[
        df[c_pass4] == "CONFLICT_TARGET_ASSIGNED_BUT_SOURCE_NONVIRAL"
    ].copy()

    cif_cache: dict[str, dict[str, Any]] = {}
    rows = []

    for _, r in conflicts.iterrows():
        pdb = clean(r[c_pdb]).upper()
        virus = clean(r[c_virus])
        cif_path = Path(clean(r[c_cif]))
        entity_ids = split_ids(r[c_entity_ids])

        rec = dict(r)
        rec.update({
            "pass6_entity_forensics": "",
            "pass6_conflict_class": "",
            "pass6_recommendation": "",
            "pass6_confidence": "",
            "pass6_error": "",
        })

        if not cif_path.exists():
            rec.update({
                "pass6_conflict_class": "CIF_NOT_FOUND",
                "pass6_recommendation": "MANUAL_REVIEW",
                "pass6_confidence": "LOW",
                "pass6_error": str(cif_path),
            })
            rows.append(rec)
            continue

        key = str(cif_path.resolve())
        if key not in cif_cache:
            try:
                cif_cache[key] = MMCIF2Dict(str(cif_path))
            except Exception as exc:
                cif_cache[key] = {"__ERROR__": f"{type(exc).__name__}: {exc}"}

        cif = cif_cache[key]
        if "__ERROR__" in cif:
            rec.update({
                "pass6_conflict_class": "CIF_PARSE_ERROR",
                "pass6_recommendation": "MANUAL_REVIEW",
                "pass6_confidence": "LOW",
                "pass6_error": cif["__ERROR__"],
            })
            rows.append(rec)
            continue

        entity_desc_map = parse_entity_desc(cif)
        source_rows = parse_pass4_source_details(
            r[c_source_details] if c_source_details else ""
        )
        source_by_entity = {
            clean(x.get("entity_id")): x
            for x in source_rows
            if clean(x.get("entity_id"))
        }

        results = []
        for eid in entity_ids:
            desc = clean(entity_desc_map.get(eid))
            refs = refs_for_entity(cif, eid)
            result = classify_entity(
                virus=virus,
                entity_id=eid,
                description=desc,
                refs=refs,
                seed=seed,
                source_detail=source_by_entity.get(eid),
            )
            results.append(result)

        conflict_class, recommendation, confidence = aggregate_entities(results)

        def uniq(vals):
            seen = set()
            out_vals = []
            for x in vals:
                x = clean(x)
                k = x.casefold()
                if x and k not in seen:
                    seen.add(k)
                    out_vals.append(x)
            return out_vals

        entity_signature = ";".join(
            uniq([x.get("entity_description", "") for x in results])
        )
        source_signature_values = []
        for x in results:
            source_signature_values.extend(x.get("gene_sources", []) or [])
            source_signature_values.extend(x.get("natural_sources", []) or [])
        source_signature = ";".join(uniq(source_signature_values))

        rec.update({
            "pass6_entity_forensics": json.dumps(
                results, ensure_ascii=False, sort_keys=True
            ),
            "pass6_entity_description_signature": entity_signature,
            "pass6_biological_source_signature": source_signature,
            "pass6_conflict_class": conflict_class,
            "pass6_recommendation": recommendation,
            "pass6_confidence": confidence,
        })
        rows.append(rec)

    out = pd.DataFrame(rows)
    if not out.empty:
        out = out.sort_values(
            ["pass6_conflict_class", c_virus, c_pdb],
            kind="mergesort",
        )
    out.to_csv(
        args.outdir / "taxonomy_source_conflict_forensics_pass6.csv",
        index=False,
    )

    # Signature-level report to make manual review manageable.
    signature_rows = []
    if not out.empty:
        group_cols = [
            c_virus,
            "pass6_conflict_class",
            "pass6_recommendation",
            "pass6_confidence",
            "pass6_entity_description_signature",
            "pass6_biological_source_signature",
        ]
        if c_target_name:
            group_cols.append(c_target_name)
        elif c_target_id:
            group_cols.append(c_target_id)

        for keys, grp in out.groupby(group_cols, dropna=False, sort=True):
            if not isinstance(keys, tuple):
                keys = (keys,)
            d = dict(zip(group_cols, keys))

            descriptions = []
            source_names = []
            for raw in grp["pass6_entity_forensics"]:
                try:
                    ents = json.loads(raw) if raw else []
                except Exception:
                    ents = []
                for e in ents:
                    if clean(e.get("entity_description")):
                        descriptions.append(clean(e["entity_description"]))
                    source_names.extend(e.get("gene_sources", []) or [])
                    source_names.extend(e.get("natural_sources", []) or [])

            def uniq(vals):
                seen = set()
                result = []
                for x in vals:
                    x = clean(x)
                    k = x.casefold()
                    if x and k not in seen:
                        seen.add(k)
                        result.append(x)
                return result

            d.update({
                "occurrence_count": int(len(grp)),
                "distinct_pdbs": int(grp[c_pdb].nunique()),
                "entity_descriptions": clean(grp["pass6_entity_description_signature"].iloc[0]),
                "biological_sources": clean(grp["pass6_biological_source_signature"].iloc[0]),
                "pdb_examples": ";".join(
                    sorted(set(grp[c_pdb].astype(str)))[:30]
                ),
            })
            signature_rows.append(d)

    sig = pd.DataFrame(signature_rows)
    if not sig.empty:
        sig = sig.sort_values(
            ["occurrence_count", "distinct_pdbs"],
            ascending=[False, False],
            kind="mergesort",
        )
    sig.to_csv(
        args.outdir / "taxonomy_source_conflict_signatures_pass6.csv",
        index=False,
    )

    # Everything remains auditable; "recommended" is not a production write.
    review = out.copy()
    review.to_csv(
        args.outdir / "taxonomy_source_conflict_review_queue_pass6.csv",
        index=False,
    )

    class_counts = Counter(out["pass6_conflict_class"]) if not out.empty else Counter()
    rec_counts = Counter(out["pass6_recommendation"]) if not out.empty else Counter()
    confidence_counts = Counter(out["pass6_confidence"]) if not out.empty else Counter()

    summary = {
        "input_conflict_occurrences": int(len(out)),
        "distinct_conflict_structures": int(out[c_pdb].nunique()) if not out.empty else 0,
        "conflict_class_counts": dict(sorted(class_counts.items())),
        "recommendation_counts": dict(sorted(rec_counts.items())),
        "confidence_counts": dict(sorted(confidence_counts.items())),
        "production_data_modified": False,
        "policy": {
            "folder_label_used_as_deciding_evidence": False,
            "fuzzy_matching": False,
            "reviewed_seed_exact_roles_used": True,
            "mmcif_struct_ref_used_as_independent_evidence": True,
            "recommendations_written_to_production": False,
        },
    }
    (args.outdir / "taxonomy_source_conflict_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    print(json.dumps(summary, indent=2, sort_keys=True))
    print(f"\nWrote audit-only Pass-6 outputs to: {args.outdir}")


if __name__ == "__main__":
    main()
