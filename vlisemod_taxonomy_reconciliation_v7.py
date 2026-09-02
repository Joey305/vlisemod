#!/usr/bin/env python3
"""
V-LiSEMOD taxonomy reconciliation — Pass 7
Forensic audit of the 40 remaining polyprotein/domain cases (READ ONLY)

Purpose
-------
Take the Pass-5 unresolved polyprotein review queue and explain exactly why each
occurrence could not be resolved automatically.

This pass DOES NOT modify:
- production SQLite
- Stage-09 / Stage-12 / Stage-14
- structure_classifications
- PROTACability scores
- source CIF/mmCIF files
- API/UI

It is an audit/triage pass only.

Input
-----
--pass5-review-csv
    taxonomy_polyprotein_domain_review_queue.csv

Outputs
-------
taxonomy_polyprotein_forensics_pass7.csv
taxonomy_polyprotein_forensic_signatures_pass7.csv
taxonomy_polyprotein_forensics_summary.json
taxonomy_polyprotein_manual_review_pass7.csv

Goals
-----
For each unresolved occurrence:
1. inspect contacted entity IDs/chains/descriptions;
2. inspect _struct_ref and _struct_ref_seq metadata;
3. classify the failure reason;
4. summarize unique PDB/evidence signatures;
5. make NO automatic mature-product assignment.

No fuzzy matching.
No folder-label resolution.
"""

from __future__ import annotations

import argparse
import json
import re
from collections import Counter
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


def split_multi(v: Any) -> list[str]:
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


def to_int(v: Any) -> int | None:
    s = clean(v)
    if not s:
        return None
    m = re.search(r"-?\d+", s)
    return int(m.group()) if m else None


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


def uniq(vals):
    out, seen = [], set()
    for x in vals:
        x = clean(x)
        k = x.casefold()
        if x and k not in seen:
            seen.add(k)
            out.append(x)
    return out


# ---------------------------------------------------------------------------
# mmCIF parsing
# ---------------------------------------------------------------------------

def entity_desc_map(cif: dict[str, Any]) -> dict[str, str]:
    ids = as_list(cif.get("_entity.id"))
    desc = as_list(cif.get("_entity.pdbx_description"))
    return {eid: value_at(desc, i) for i, eid in enumerate(ids) if eid}


def struct_asym_map(cif: dict[str, Any]) -> dict[str, str]:
    ids = as_list(cif.get("_struct_asym.id"))
    entity_ids = as_list(cif.get("_struct_asym.entity_id"))
    return {a: value_at(entity_ids, i) for i, a in enumerate(ids) if a}


def parse_refs(cif: dict[str, Any]) -> dict[str, dict[str, str]]:
    ids = as_list(cif.get("_struct_ref.id"))
    entity_ids = as_list(cif.get("_struct_ref.entity_id"))
    db_names = as_list(cif.get("_struct_ref.db_name"))
    db_codes = as_list(cif.get("_struct_ref.db_code"))
    accessions = as_list(cif.get("_struct_ref.pdbx_db_accession"))

    out = {}
    for i, rid in enumerate(ids):
        if not rid:
            continue
        out[rid] = {
            "ref_id": rid,
            "entity_id": value_at(entity_ids, i),
            "db_name": value_at(db_names, i),
            "db_code": value_at(db_codes, i),
            "accession": value_at(accessions, i),
        }
    return out


def parse_ref_seq(cif: dict[str, Any]) -> list[dict[str, Any]]:
    fields = {
        "align_id": as_list(cif.get("_struct_ref_seq.align_id")),
        "ref_id": as_list(cif.get("_struct_ref_seq.ref_id")),
        "strand": as_list(cif.get("_struct_ref_seq.pdbx_strand_id")),
        "seq_beg": as_list(cif.get("_struct_ref_seq.seq_align_beg")),
        "seq_end": as_list(cif.get("_struct_ref_seq.seq_align_end")),
        "db_beg": as_list(cif.get("_struct_ref_seq.db_align_beg")),
        "db_end": as_list(cif.get("_struct_ref_seq.db_align_end")),
        "accession": as_list(cif.get("_struct_ref_seq.pdbx_db_accession")),
    }
    n = max((len(v) for v in fields.values()), default=0)
    rows = []
    for i in range(n):
        rows.append({
            "align_id": value_at(fields["align_id"], i),
            "ref_id": value_at(fields["ref_id"], i),
            "strand_ids": split_multi(value_at(fields["strand"], i)),
            "seq_align_beg": to_int(value_at(fields["seq_beg"], i)),
            "seq_align_end": to_int(value_at(fields["seq_end"], i)),
            "db_align_beg": to_int(value_at(fields["db_beg"], i)),
            "db_align_end": to_int(value_at(fields["db_end"], i)),
            "accession": value_at(fields["accession"], i),
        })
    return rows


def entity_alignments(
    cif: dict[str, Any],
    entity_id: str,
    contact_chains: list[str],
) -> list[dict[str, Any]]:
    refs = parse_refs(cif)
    ref_seq = parse_ref_seq(cif)
    asym = struct_asym_map(cif)

    target_ref_ids = {
        rid for rid, r in refs.items()
        if clean(r.get("entity_id")) == clean(entity_id)
    }

    out = []
    for a in ref_seq:
        via_ref = bool(a["ref_id"] and a["ref_id"] in target_ref_ids)
        via_chain = False

        for strand in a.get("strand_ids", []):
            if strand in contact_chains:
                via_chain = True
            if asym.get(strand) == entity_id:
                via_chain = True

        if not via_ref and not via_chain:
            continue

        meta = refs.get(a["ref_id"], {})
        rec = dict(a)
        rec.update({
            "ref_entity_id": meta.get("entity_id", ""),
            "ref_db_name": meta.get("db_name", ""),
            "ref_db_code": meta.get("db_code", ""),
            "ref_accession": a.get("accession") or meta.get("accession", ""),
        })
        out.append(rec)

    seen = set()
    dedup = []
    for r in out:
        key = (
            r.get("ref_id"),
            r.get("db_align_beg"),
            r.get("db_align_end"),
            tuple(r.get("strand_ids", [])),
            r.get("ref_accession"),
        )
        if key not in seen:
            seen.add(key)
            dedup.append(r)
    return dedup


# ---------------------------------------------------------------------------
# Failure classification
# ---------------------------------------------------------------------------

def is_true_gag_pol(desc: str, refs: list[dict[str, Any]]) -> bool:
    d = norm(desc)
    ref_text = " ".join(
        norm(" ".join([
            clean(r.get("ref_db_code")),
            clean(r.get("ref_accession")),
            clean(r.get("ref_db_name")),
        ]))
        for r in refs
    )
    return (
        bool(re.search(r"\bgag pol\b", d))
        or "pr160" in d
        or bool(re.search(r"\bgag pol\b", ref_text))
        or "pr160" in ref_text
    )


def looks_like_sars_replicase(desc: str, refs: list[dict[str, Any]]) -> bool:
    s = norm(desc)
    ref = " ".join(norm(r.get("ref_db_code")) for r in refs)
    return (
        "replicase polyprotein" in s
        or "orf1a" in s
        or "orf1ab" in s
        or "r1ab sars2" in ref
        or "r1a sars2" in ref
    )


def looks_like_hiv_gag(desc: str, refs: list[dict[str, Any]]) -> bool:
    if is_true_gag_pol(desc, refs):
        return False
    s = norm(desc)
    ref = " ".join(norm(r.get("ref_db_code")) for r in refs)
    return "gag polyprotein" in s or s == "gag protein" or re.search(r"\bgag\b", ref) is not None


def looks_like_hiv_pol(desc: str, refs: list[dict[str, Any]]) -> bool:
    if is_true_gag_pol(desc, refs):
        return False
    s = norm(desc)
    ref = " ".join(norm(r.get("ref_db_code")) for r in refs)
    return (
        s in {"pol", "pol protein"}
        or "pol protein" in s
        or re.search(r"\bpol\b", ref) is not None
    )


def alignment_range_class(records: list[dict[str, Any]]) -> str:
    if not records:
        return "NO_ALIGNMENT_RECORDS"

    usable = [
        r for r in records
        if r.get("db_align_beg") is not None and r.get("db_align_end") is not None
    ]
    if not usable:
        return "ALIGNMENT_WITHOUT_DB_COORDINATES"

    ranges = {
        (min(r["db_align_beg"], r["db_align_end"]),
         max(r["db_align_beg"], r["db_align_end"]))
        for r in usable
    }
    if len(ranges) > 1:
        return "MULTIPLE_DISTINCT_DB_RANGES"

    return "SINGLE_DB_RANGE"


def classify_failure(
    virus: str,
    desc: str,
    records: list[dict[str, Any]],
    entity_count: int,
) -> tuple[str, str]:
    """
    Return (failure_class, recommended_next_action).
    No target is assigned here.
    """
    if entity_count > 1:
        return (
            "MULTIPLE_CONTACTED_ENTITIES",
            "REVIEW_ENTITY_BY_ENTITY",
        )

    if not records:
        return (
            "NO_STRUCT_REF_ALIGNMENT",
            "USE_SEQUENCE_ALIGNMENT_OR_MANUAL_PDB_REVIEW",
        )

    range_class = alignment_range_class(records)
    if range_class == "ALIGNMENT_WITHOUT_DB_COORDINATES":
        return (
            "STRUCT_REF_PRESENT_NO_DB_COORDINATES",
            "USE_SEQUENCE_ALIGNMENT_OR_MANUAL_PDB_REVIEW",
        )
    if range_class == "MULTIPLE_DISTINCT_DB_RANGES":
        return (
            "MULTIPLE_DATABASE_ALIGNMENT_RANGES",
            "REVIEW_ALIGNMENT_CHOICE",
        )

    if is_true_gag_pol(desc, records):
        return (
            "TRUE_HIV_GAG_POL_PRECURSOR",
            "MAP_GAG_POL_COORDINATE_SYSTEM_OR_SEQUENCE_ALIGN",
        )

    if looks_like_sars_replicase(desc, records):
        return (
            "SARS_REPLICASE_RANGE_NOT_SINGLE_MATURE_PRODUCT",
            "INSPECT_RANGE_AND_CONTACT_RESIDUES",
        )

    if looks_like_hiv_gag(desc, records):
        return (
            "HIV_GAG_RANGE_NOT_SINGLE_MATURE_PRODUCT",
            "INSPECT_RANGE_AND_CONTACT_RESIDUES",
        )

    if looks_like_hiv_pol(desc, records):
        return (
            "HIV_POL_RANGE_NOT_SINGLE_MATURE_PRODUCT",
            "INSPECT_RANGE_AND_CONTACT_RESIDUES",
        )

    return (
        "UNRECOGNIZED_POLYPROTEIN_ALIGNMENT_CONTEXT",
        "MANUAL_PDB_AND_SEQUENCE_REVIEW",
    )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--pass5-review-csv",
        type=Path,
        required=True,
        help="taxonomy_polyprotein_domain_review_queue.csv",
    )
    ap.add_argument("--outdir", type=Path, required=True)
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.pass5_review_csv, dtype=str, low_memory=False).fillna("")

    c_pdb = resolve_col(df, ["pdb_id", "pdb"])
    c_virus = resolve_col(df, ["virus_name", "virus"])
    c_cif = resolve_col(df, ["pass4_selected_cif_path"])
    c_entity_ids = resolve_col(
        df, ["contacting_entity_ids", "contact_entity_ids", "entity_ids"], required=False
    )
    c_desc = resolve_col(
        df, ["contacting_entity_descriptions", "entity_descriptions"], required=False
    )
    c_chains = resolve_col(
        df, ["contacting_protein_chains", "target_chains", "protein_contact_chains"],
        required=False,
    )
    c_pass5 = resolve_col(df, ["pass5_status"])

    rows = []
    cache: dict[str, dict[str, Any]] = {}

    for _, r in df.iterrows():
        rec = dict(r)

        pdb = clean(r[c_pdb]).upper()
        virus = clean(r[c_virus])
        cif_path = Path(clean(r[c_cif]))
        entity_ids = split_multi(r[c_entity_ids]) if c_entity_ids else []
        chains = split_multi(r[c_chains]) if c_chains else []
        desc_raw = clean(r[c_desc]) if c_desc else ""

        rec.update({
            "pass7_failure_class": "",
            "pass7_recommended_next_action": "",
            "pass7_entity_forensics": "",
            "pass7_error": "",
        })

        if not cif_path.exists():
            rec.update({
                "pass7_failure_class": "CIF_NOT_FOUND",
                "pass7_recommended_next_action": "MANUAL_PDB_REVIEW",
                "pass7_error": str(cif_path),
            })
            rows.append(rec)
            continue

        key = str(cif_path.resolve())
        if key not in cache:
            try:
                cache[key] = MMCIF2Dict(str(cif_path))
            except Exception as exc:
                cache[key] = {"__ERROR__": f"{type(exc).__name__}: {exc}"}

        cif = cache[key]
        if "__ERROR__" in cif:
            rec.update({
                "pass7_failure_class": "CIF_PARSE_ERROR",
                "pass7_recommended_next_action": "MANUAL_PDB_REVIEW",
                "pass7_error": cif["__ERROR__"],
            })
            rows.append(rec)
            continue

        desc_map = entity_desc_map(cif)

        entity_results = []
        failure_classes = []
        next_actions = []

        if not entity_ids:
            failure_classes.append("NO_CONTACTING_ENTITY_ID")
            next_actions.append("MANUAL_PDB_REVIEW")
        else:
            for eid in entity_ids:
                desc = clean(desc_map.get(eid)) or desc_raw
                aligns = entity_alignments(cif, eid, chains)

                fclass, action = classify_failure(
                    virus=virus,
                    desc=desc,
                    records=aligns,
                    entity_count=len(entity_ids),
                )
                failure_classes.append(fclass)
                next_actions.append(action)

                entity_results.append({
                    "entity_id": eid,
                    "entity_description": desc,
                    "contacting_chains": chains,
                    "alignment_records": aligns,
                    "failure_class": fclass,
                    "recommended_next_action": action,
                })

        rec.update({
            "pass7_failure_class": ";".join(uniq(failure_classes)),
            "pass7_recommended_next_action": ";".join(uniq(next_actions)),
            "pass7_entity_forensics": json.dumps(
                entity_results, ensure_ascii=False, sort_keys=True
            ),
        })
        rows.append(rec)

    out = pd.DataFrame(rows)
    out = out.sort_values(
        ["pass7_failure_class", c_virus, c_pdb],
        kind="mergesort",
    )
    out.to_csv(
        args.outdir / "taxonomy_polyprotein_forensics_pass7.csv",
        index=False,
    )

    # Compact signatures: keep PDB, description, failure type and alignment
    # evidence distinct so unrelated situations never get collapsed together.
    signatures = []
    if not out.empty:
        group_cols = [
            c_virus,
            c_pdb,
            c_pass5,
            "pass7_failure_class",
            "pass7_recommended_next_action",
        ]
        for keys, grp in out.groupby(group_cols, dropna=False, sort=True):
            if not isinstance(keys, tuple):
                keys = (keys,)
            d = dict(zip(group_cols, keys))

            descriptions = []
            accessions = []
            db_codes = []
            ranges = []
            entity_ids_seen = []

            for raw in grp["pass7_entity_forensics"]:
                try:
                    ents = json.loads(raw) if raw else []
                except Exception:
                    ents = []
                for e in ents:
                    entity_ids_seen.append(clean(e.get("entity_id")))
                    descriptions.append(clean(e.get("entity_description")))
                    for a in e.get("alignment_records", []) or []:
                        accessions.append(clean(a.get("ref_accession")))
                        db_codes.append(clean(a.get("ref_db_code")))
                        b = a.get("db_align_beg")
                        en = a.get("db_align_end")
                        if b is not None and en is not None:
                            ranges.append(f"{min(b,en)}-{max(b,en)}")

            d.update({
                "occurrence_count": int(len(grp)),
                "entity_ids": ";".join(uniq(entity_ids_seen)),
                "entity_descriptions": ";".join(uniq(descriptions)),
                "reference_accessions": ";".join(uniq(accessions)),
                "reference_db_codes": ";".join(uniq(db_codes)),
                "db_alignment_ranges": ";".join(uniq(ranges)),
            })
            signatures.append(d)

    sig = pd.DataFrame(signatures)
    if not sig.empty:
        sig = sig.sort_values(
            ["occurrence_count", "pass7_failure_class", c_pdb],
            ascending=[False, True, True],
            kind="mergesort",
        )
    sig.to_csv(
        args.outdir / "taxonomy_polyprotein_forensic_signatures_pass7.csv",
        index=False,
    )

    out.to_csv(
        args.outdir / "taxonomy_polyprotein_manual_review_pass7.csv",
        index=False,
    )

    failure_counts = Counter(out["pass7_failure_class"].tolist())
    action_counts = Counter(out["pass7_recommended_next_action"].tolist())

    summary = {
        "input_review_occurrences": int(len(out)),
        "distinct_review_structures": int(out[c_pdb].nunique()) if not out.empty else 0,
        "failure_class_counts": dict(sorted(failure_counts.items())),
        "recommended_action_counts": dict(sorted(action_counts.items())),
        "production_data_modified": False,
        "policy": {
            "automatic_target_assignment": False,
            "folder_label_used_as_deciding_evidence": False,
            "fuzzy_matching": False,
            "scientific_scores_modified": False,
        },
    }

    (args.outdir / "taxonomy_polyprotein_forensics_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    print(json.dumps(summary, indent=2, sort_keys=True))
    print(f"\nWrote audit-only Pass-7 outputs to: {args.outdir}")


if __name__ == "__main__":
    main()
