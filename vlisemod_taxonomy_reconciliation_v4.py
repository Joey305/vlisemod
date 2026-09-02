#!/usr/bin/env python3
"""
V-LiSEMOD taxonomy reconciliation — Pass 4
Source-organism / entity-role QC (READ ONLY)

Purpose
-------
Validate Pass-2/Pass-3 canonical target assignments against the source organism
of the ligand-contacting mmCIF entities.

This script does NOT modify:
- the production SQLite database
- structure_classifications
- Stage-09
- Stage-12
- Stage-14
- PROTACability scores
- API/UI code
- source CIF files

It only writes new audit artifacts.

Inputs
------
--pass3-csv
    taxonomy_role_reconciliation_pass3.csv

--repo-root
    VLISEMOD project root. The script indexes *.cif / *.mmcif files below it.

Optional:
--cif-root
    Repeatable explicit CIF root. If supplied, these are indexed instead of
    recursively scanning --repo-root.

Outputs
-------
taxonomy_source_organism_qc_pass4.csv
taxonomy_source_organism_qc_summary.json
taxonomy_source_organism_conflicts.csv
taxonomy_source_organism_review_queue.csv
taxonomy_target_browser_eligibility_pass4.csv

Design principles
-----------------
1. Use ligand-contacting entity IDs already established by the prior audit.
2. Use mmCIF gene/natural source organism as the biological source.
3. DO NOT treat the recombinant expression host as the protein source.
4. Viral source match confirms target identity; a nonviral gene/natural source
   contradicts a provisional viral-target assignment.
5. Synthetic/unknown source is reviewable, not automatically wrong.
6. Polyprotein/domain ambiguity remains unresolved even when viral source is
   confirmed.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import pandas as pd

try:
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
except Exception as exc:
    raise SystemExit(
        "Biopython is required for Pass 4. Install/activate the same environment "
        "used by the V-LiSEMOD rebuild pipeline. Import error: " + repr(exc)
    )


# ---------------------------------------------------------------------------
# Text helpers
# ---------------------------------------------------------------------------

NULLISH = {"", ".", "?", "none", "null", "n/a", "na", "unknown", "not applicable"}


def txt(value: Any) -> str:
    if value is None or pd.isna(value):
        return ""
    return re.sub(r"\s+", " ", str(value).strip())


def clean_cif_value(value: Any) -> str:
    s = txt(value).strip("'\"")
    if s.casefold() in NULLISH:
        return ""
    return s


def as_list(value: Any) -> list[str]:
    if value is None:
        return []
    if isinstance(value, (list, tuple)):
        return [clean_cif_value(x) for x in value]
    return [clean_cif_value(value)]


def split_ids(value: Any) -> list[str]:
    s = txt(value)
    if not s:
        return []
    vals = [x.strip() for x in re.split(r"[;,|]", s) if x.strip()]
    out = []
    seen = set()
    for v in vals:
        k = v.casefold()
        if k not in seen:
            seen.add(k)
            out.append(v)
    return out


def resolve_col(df: pd.DataFrame, candidates: list[str], required: bool = True) -> str | None:
    by_lower = {c.casefold(): c for c in df.columns}
    for name in candidates:
        if name.casefold() in by_lower:
            return by_lower[name.casefold()]
    if required:
        raise KeyError(
            f"Missing required column; expected one of {candidates}. "
            f"Available columns: {list(df.columns)}"
        )
    return None


def norm_for_match(value: Any) -> str:
    s = clean_cif_value(value).casefold()
    s = s.replace("_", " ")
    s = re.sub(r"[^a-z0-9]+", " ", s)
    return re.sub(r"\s+", " ", s).strip()


# ---------------------------------------------------------------------------
# Virus/source matching
# ---------------------------------------------------------------------------

VIRUS_SOURCE_PATTERNS = {
    "sars_cov_2": [
        r"\bsevere acute respiratory syndrome coronavirus 2\b",
        r"\bsars cov 2\b",
        r"\b2019 ncov\b",
    ],
    "hiv_1": [
        r"\bhuman immunodeficiency virus 1\b",
        r"\bhuman immunodeficiency virus type 1\b",
        r"\bhiv 1\b",
    ],
}

NONVIRAL_SOURCE_HINTS = [
    # Common host / expression / binder sources. A non-match is not called
    # nonviral solely because it is absent from this list; this list strengthens
    # the verdict for common cases.
    "homo sapiens",
    "oryctolagus cuniculus",
    "mus musculus",
    "escherichia coli",
    "spodoptera frugiperda",
    "trichoplusia ni",
    "bos taurus",
    "sus scrofa",
    "gallus gallus",
    "cricetulus griseus",
    "synthetic construct",
    "enterobacteria phage t4",
]


def expected_virus_key(virus_name: str) -> str:
    s = norm_for_match(virus_name)
    if "sars" in s and ("cov" in s or "coronavirus" in s) and "2" in s:
        return "sars_cov_2"
    if ("hiv" in s and "1" in s) or "human immunodeficiency virus 1" in s:
        return "hiv_1"
    return re.sub(r"\s+", "_", s)


def source_matches_virus(virus_name: str, source_name: str) -> bool | None:
    """
    True  = source clearly matches the expected virus
    False = source is present and clearly does not match
    None  = missing/ambiguous/synthetic source
    """
    src = norm_for_match(source_name)
    if not src:
        return None

    # Synthetic source is not sufficient for organism confirmation.
    if "synthetic" in src:
        return None

    key = expected_virus_key(virus_name)
    patterns = VIRUS_SOURCE_PATTERNS.get(key, [])
    for pat in patterns:
        if re.search(pat, src):
            return True

    # Generic fallback for viruses not explicitly listed:
    # require at least two informative tokens from the V-LiSEMOD virus label.
    if not patterns:
        expected_tokens = [
            t for t in norm_for_match(virus_name).split()
            if len(t) >= 3 and t not in {"virus", "protein"}
        ]
        if len(expected_tokens) >= 2 and sum(t in src.split() for t in expected_tokens) >= 2:
            return True

    # A populated gene/natural source that does not match the expected virus
    # is treated as a nonviral source for this virus context.
    return False


# ---------------------------------------------------------------------------
# mmCIF extraction
# ---------------------------------------------------------------------------

def category_rows(
    cif: dict[str, Any],
    entity_tag: str,
    value_tags: list[str],
) -> list[dict[str, str]]:
    entity_ids = as_list(cif.get(entity_tag))
    if not entity_ids:
        return []

    columns = {tag: as_list(cif.get(tag)) for tag in value_tags}
    n = len(entity_ids)
    rows = []
    for i in range(n):
        row = {"entity_id": clean_cif_value(entity_ids[i])}
        for tag, vals in columns.items():
            row[tag] = clean_cif_value(vals[i]) if i < len(vals) else ""
        rows.append(row)
    return rows


def parse_cif_entity_sources(path: Path) -> dict[str, Any]:
    cif = MMCIF2Dict(str(path))

    entity_desc = {}
    entity_ids = as_list(cif.get("_entity.id"))
    entity_descs = as_list(cif.get("_entity.pdbx_description"))
    for i, eid in enumerate(entity_ids):
        if not eid:
            continue
        entity_desc[eid] = entity_descs[i] if i < len(entity_descs) else ""

    result: dict[str, Any] = {
        "entities": defaultdict(lambda: {
            "description": "",
            "gene_source": [],
            "natural_source": [],
            "synthetic_source": [],
            "expression_host": [],
            "tax_ids": [],
        })
    }

    for eid, desc in entity_desc.items():
        result["entities"][eid]["description"] = desc

    # Genetically engineered/recombinant source.
    gen_rows = category_rows(
        cif,
        "_entity_src_gen.entity_id",
        [
            "_entity_src_gen.pdbx_gene_src_scientific_name",
            "_entity_src_gen.pdbx_gene_src_ncbi_taxonomy_id",
            "_entity_src_gen.pdbx_host_org_scientific_name",
            "_entity_src_gen.pdbx_host_org_ncbi_taxonomy_id",
        ],
    )
    for r in gen_rows:
        eid = r["entity_id"]
        if not eid:
            continue
        g = clean_cif_value(r.get("_entity_src_gen.pdbx_gene_src_scientific_name"))
        gt = clean_cif_value(r.get("_entity_src_gen.pdbx_gene_src_ncbi_taxonomy_id"))
        h = clean_cif_value(r.get("_entity_src_gen.pdbx_host_org_scientific_name"))
        ht = clean_cif_value(r.get("_entity_src_gen.pdbx_host_org_ncbi_taxonomy_id"))
        if g:
            result["entities"][eid]["gene_source"].append(g)
        if gt:
            result["entities"][eid]["tax_ids"].append(gt)
        if h:
            result["entities"][eid]["expression_host"].append(h)
        if ht:
            result["entities"][eid]["tax_ids"].append("host:" + ht)

    # Natural source.
    nat_rows = category_rows(
        cif,
        "_entity_src_nat.entity_id",
        [
            "_entity_src_nat.pdbx_organism_scientific",
            "_entity_src_nat.pdbx_ncbi_taxonomy_id",
        ],
    )
    for r in nat_rows:
        eid = r["entity_id"]
        if not eid:
            continue
        org = clean_cif_value(r.get("_entity_src_nat.pdbx_organism_scientific"))
        tax = clean_cif_value(r.get("_entity_src_nat.pdbx_ncbi_taxonomy_id"))
        if org:
            result["entities"][eid]["natural_source"].append(org)
        if tax:
            result["entities"][eid]["tax_ids"].append(tax)

    # Synthetic source.
    syn_rows = category_rows(
        cif,
        "_pdbx_entity_src_syn.entity_id",
        [
            "_pdbx_entity_src_syn.organism_scientific",
            "_pdbx_entity_src_syn.ncbi_taxonomy_id",
        ],
    )
    for r in syn_rows:
        eid = r["entity_id"]
        if not eid:
            continue
        org = clean_cif_value(r.get("_pdbx_entity_src_syn.organism_scientific"))
        tax = clean_cif_value(r.get("_pdbx_entity_src_syn.ncbi_taxonomy_id"))
        if org:
            result["entities"][eid]["synthetic_source"].append(org)
        if tax:
            result["entities"][eid]["tax_ids"].append(tax)

    # De-duplicate while preserving order.
    for eid, info in result["entities"].items():
        for field in ["gene_source", "natural_source", "synthetic_source", "expression_host", "tax_ids"]:
            seen = set()
            vals = []
            for x in info[field]:
                k = x.casefold()
                if x and k not in seen:
                    seen.add(k)
                    vals.append(x)
            info[field] = vals

    result["entities"] = dict(result["entities"])
    return result


def primary_biological_sources(info: dict[str, Any]) -> list[str]:
    """
    Gene source and natural source are biological-origin evidence.
    Expression host is intentionally excluded.
    Synthetic source is included only as weak/unknown evidence.
    """
    values = []
    for field in ["gene_source", "natural_source"]:
        values.extend(info.get(field, []))
    if not values:
        values.extend(info.get("synthetic_source", []))

    out = []
    seen = set()
    for x in values:
        if x and x.casefold() not in seen:
            seen.add(x.casefold())
            out.append(x)
    return out


def entity_source_verdict(virus_name: str, info: dict[str, Any]) -> str:
    sources = primary_biological_sources(info)
    if not sources:
        return "SOURCE_UNKNOWN"

    verdicts = [source_matches_virus(virus_name, s) for s in sources]
    if any(v is True for v in verdicts):
        if any(v is False for v in verdicts):
            return "SOURCE_MIXED"
        return "VIRAL_SOURCE_MATCH"
    if all(v is False for v in verdicts):
        return "NONVIRAL_SOURCE"
    return "SOURCE_UNKNOWN"


# ---------------------------------------------------------------------------
# CIF discovery
# ---------------------------------------------------------------------------

def build_cif_index(roots: list[Path]) -> dict[str, list[Path]]:
    index: dict[str, list[Path]] = defaultdict(list)
    seen_paths = set()

    for root in roots:
        root = root.resolve()
        if not root.exists():
            print(f"WARNING: CIF root does not exist: {root}", file=sys.stderr)
            continue

        patterns = ("*.cif", "*.mmcif")
        for pattern in patterns:
            for path in root.rglob(pattern):
                rp = path.resolve()
                if rp in seen_paths:
                    continue
                seen_paths.add(rp)
                index[path.stem.upper()].append(rp)

    for pdb in index:
        index[pdb] = sorted(index[pdb], key=lambda p: str(p).casefold())
    return dict(index)


def choose_cif(paths: list[Path]) -> Path | None:
    if not paths:
        return None
    # Deterministic only. Classification provenance is not used to select a
    # biologically "preferred" folder because that would reintroduce the bias
    # we are auditing.
    return sorted(paths, key=lambda p: str(p).casefold())[0]


# ---------------------------------------------------------------------------
# Occurrence-level QC
# ---------------------------------------------------------------------------

def aggregate_source_verdict(entity_verdicts: list[str]) -> str:
    vals = set(entity_verdicts)
    if not vals:
        return "NO_ENTITY_SOURCE_EVIDENCE"
    if vals == {"VIRAL_SOURCE_MATCH"}:
        return "ALL_CONTACT_ENTITIES_VIRAL"
    if vals == {"NONVIRAL_SOURCE"}:
        return "ALL_CONTACT_ENTITIES_NONVIRAL"
    if "VIRAL_SOURCE_MATCH" in vals and "NONVIRAL_SOURCE" in vals:
        return "MIXED_VIRAL_NONVIRAL_CONTACT"
    if "SOURCE_MIXED" in vals:
        return "MIXED_SOURCE_WITHIN_ENTITY"
    if "SOURCE_UNKNOWN" in vals and len(vals) == 1:
        return "ALL_CONTACT_ENTITY_SOURCES_UNKNOWN"
    if "VIRAL_SOURCE_MATCH" in vals and "SOURCE_UNKNOWN" in vals:
        return "VIRAL_PLUS_UNKNOWN_CONTACT"
    if "NONVIRAL_SOURCE" in vals and "SOURCE_UNKNOWN" in vals:
        return "NONVIRAL_PLUS_UNKNOWN_CONTACT"
    return "SOURCE_REVIEW"


def decide_pass4_status(pass3_status: str, aggregate: str) -> tuple[str, str]:
    """
    Returns (pass4_status, target_browser_eligibility)
    eligibility: YES / NO / REVIEW
    """

    # Nothing to validate from contact evidence.
    if pass3_status == "UNRESOLVED_NO_CONTACT_EVIDENCE":
        return "UNRESOLVED_NO_CONTACT_EVIDENCE", "REVIEW"

    # Missing entity metadata cannot be upgraded from PDB-level hint alone.
    if pass3_status == "PDB_TARGET_HINT_REQUIRES_OCCURRENCE_REVIEW":
        return "PDB_TARGET_HINT_STILL_REQUIRES_OCCURRENCE_REVIEW", "REVIEW"

    if pass3_status == "CONTACT_WITHOUT_ENTITY_METADATA_REVIEW":
        return "CONTACT_WITHOUT_ENTITY_METADATA_REVIEW", "REVIEW"

    # Explicitly excluded nonviral contexts from Pass 3.
    if pass3_status == "EXCLUDE_NONVIRAL_CONTACT_CONTEXT":
        if aggregate in {"ALL_CONTACT_ENTITIES_NONVIRAL", "NONVIRAL_PLUS_UNKNOWN_CONTACT"}:
            return "CONFIRMED_NONVIRAL_EXCLUSION", "NO"
        if aggregate == "ALL_CONTACT_ENTITIES_VIRAL":
            return "CONFLICT_PASS3_NONVIRAL_BUT_SOURCE_IS_VIRAL", "REVIEW"
        return "NONVIRAL_EXCLUSION_SOURCE_REVIEW", "NO"

    # Polyproteins: source can confirm viral origin, not the mature product.
    if pass3_status == "POLYPROTEIN_DOMAIN_REVIEW":
        if aggregate in {"ALL_CONTACT_ENTITIES_VIRAL", "VIRAL_PLUS_UNKNOWN_CONTACT"}:
            return "VIRAL_POLYPROTEIN_CONFIRMED_DOMAIN_UNRESOLVED", "REVIEW"
        if aggregate == "ALL_CONTACT_ENTITIES_NONVIRAL":
            return "CONFLICT_POLYPROTEIN_ENTITY_SOURCE_NONVIRAL", "REVIEW"
        return "POLYPROTEIN_SOURCE_OR_DOMAIN_REVIEW", "REVIEW"

    # Genuine/mixed interfaces remain reviewable by design.
    if pass3_status in {
        "MULTIVIRAL_INTERFACE_REVIEW",
        "VIRAL_NONVIRAL_INTERFACE_REVIEW",
        "VIRAL_TARGET_COMPONENT_REVIEW",
        "RESOLVED_TARGET_BUT_AUDIT_ENTITY_MAPPING_REVIEW",
    }:
        return f"{pass3_status}_SOURCE_QC", "REVIEW"

    # Core target confirmations.
    if pass3_status in {
        "KEEP_PASS2_EXACT_PENDING_ENTITY_SOURCE_QC",
        "RESOLVED_CURATED_EXACT",
        "RESOLVED_CURATED_WITH_ENGINEERED_COMPONENT",
    }:
        if aggregate == "ALL_CONTACT_ENTITIES_VIRAL":
            return "CONFIRMED_VIRAL_TARGET", "YES"
        if aggregate == "VIRAL_PLUS_UNKNOWN_CONTACT":
            return "VIRAL_TARGET_WITH_UNKNOWN_SOURCE_COMPONENT", "REVIEW"
        if aggregate == "ALL_CONTACT_ENTITIES_NONVIRAL":
            return "CONFLICT_TARGET_ASSIGNED_BUT_SOURCE_NONVIRAL", "REVIEW"
        if aggregate == "MIXED_VIRAL_NONVIRAL_CONTACT":
            return "VIRAL_NONVIRAL_INTERFACE_REVIEW", "REVIEW"
        if aggregate == "NONVIRAL_PLUS_UNKNOWN_CONTACT":
            return "NONVIRAL_PLUS_UNKNOWN_TARGET_REVIEW", "REVIEW"
        return "TARGET_SOURCE_UNKNOWN_REVIEW", "REVIEW"

    if pass3_status == "MANUAL_REVIEW_AFTER_ROLE_PASS":
        if aggregate == "ALL_CONTACT_ENTITIES_NONVIRAL":
            return "MANUAL_CASE_CONFIRMED_NONVIRAL_CONTEXT", "NO"
        if aggregate == "ALL_CONTACT_ENTITIES_VIRAL":
            return "MANUAL_CASE_VIRAL_SOURCE_CONFIRMED_TARGET_UNRESOLVED", "REVIEW"
        return "MANUAL_CASE_SOURCE_REVIEW", "REVIEW"

    # Conservative fallback.
    return "PASS4_REVIEW_UNHANDLED_PASS3_STATUS", "REVIEW"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--pass3-csv",
        type=Path,
        required=True,
        help="taxonomy_role_reconciliation_pass3.csv",
    )
    ap.add_argument(
        "--repo-root",
        type=Path,
        default=Path("."),
        help="VLISEMOD root; default current directory",
    )
    ap.add_argument(
        "--cif-root",
        type=Path,
        action="append",
        default=[],
        help="Explicit CIF root; repeatable. If omitted, indexes --repo-root.",
    )
    ap.add_argument(
        "--outdir",
        type=Path,
        required=True,
    )
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.pass3_csv, dtype=str, low_memory=False).fillna("")

    c_pdb = resolve_col(df, ["pdb_id", "pdb"])
    c_virus = resolve_col(df, ["virus_name", "virus"])
    c_entity_ids = resolve_col(
        df,
        ["contacting_entity_ids", "contact_entity_ids", "entity_ids"],
        required=False,
    )
    c_entity_desc = resolve_col(
        df,
        ["contacting_entity_descriptions", "entity_descriptions"],
        required=False,
    )
    c_pass3_status = resolve_col(df, ["pass3_status"])
    c_pass3_role = resolve_col(df, ["pass3_entity_role"], required=False)
    c_target_id = resolve_col(df, ["pass3_canonical_target_id"], required=False)
    c_target_name = resolve_col(df, ["pass3_canonical_target_name"], required=False)
    c_family = resolve_col(df, ["pass3_target_family"], required=False)

    roots = args.cif_root if args.cif_root else [args.repo_root]
    print("Indexing CIF files under:")
    for r in roots:
        print("  ", r.resolve())

    cif_index = build_cif_index(roots)
    print(f"Indexed {sum(len(v) for v in cif_index.values())} CIF/mmCIF files "
          f"for {len(cif_index)} unique PDB IDs.")

    # Parse each PDB once.
    cif_cache: dict[str, dict[str, Any]] = {}
    cif_errors: dict[str, str] = {}

    def get_cif(pdb: str) -> tuple[Path | None, dict[str, Any] | None, str]:
        pdb = pdb.upper()
        candidates = cif_index.get(pdb, [])
        selected = choose_cif(candidates)
        if not selected:
            return None, None, "CIF_NOT_FOUND"
        if pdb in cif_cache:
            return selected, cif_cache[pdb], ""
        if pdb in cif_errors:
            return selected, None, cif_errors[pdb]
        try:
            parsed = parse_cif_entity_sources(selected)
            cif_cache[pdb] = parsed
            return selected, parsed, ""
        except Exception as exc:
            msg = f"{type(exc).__name__}: {exc}"
            cif_errors[pdb] = msg
            return selected, None, msg

    output_rows = []

    for _, row in df.iterrows():
        rec = dict(row)
        pdb = txt(row[c_pdb]).upper()
        virus = txt(row[c_virus])
        pass3_status = txt(row[c_pass3_status])
        entity_ids = split_ids(row[c_entity_ids]) if c_entity_ids else []

        selected, parsed, cif_error = get_cif(pdb)
        candidates = cif_index.get(pdb, [])

        details = []
        entity_verdicts = []

        if parsed is not None and entity_ids:
            entities = parsed.get("entities", {})
            for eid in entity_ids:
                info = entities.get(eid)
                if info is None:
                    detail = {
                        "entity_id": eid,
                        "description": "",
                        "gene_source": [],
                        "natural_source": [],
                        "synthetic_source": [],
                        "expression_host": [],
                        "tax_ids": [],
                        "source_verdict": "ENTITY_ID_NOT_FOUND_IN_CIF",
                    }
                    details.append(detail)
                    entity_verdicts.append("SOURCE_UNKNOWN")
                    continue

                verdict = entity_source_verdict(virus, info)
                detail = {
                    "entity_id": eid,
                    "description": info.get("description", ""),
                    "gene_source": info.get("gene_source", []),
                    "natural_source": info.get("natural_source", []),
                    "synthetic_source": info.get("synthetic_source", []),
                    "expression_host": info.get("expression_host", []),
                    "tax_ids": info.get("tax_ids", []),
                    "source_verdict": verdict,
                }
                details.append(detail)
                entity_verdicts.append(verdict)

        elif entity_ids and parsed is None:
            entity_verdicts = ["SOURCE_UNKNOWN"] * len(entity_ids)

        aggregate = aggregate_source_verdict(entity_verdicts)
        pass4_status, eligibility = decide_pass4_status(pass3_status, aggregate)

        if cif_error == "CIF_NOT_FOUND":
            pass4_status = "CIF_NOT_FOUND_REVIEW"
            eligibility = "REVIEW"
        elif cif_error:
            pass4_status = "CIF_PARSE_ERROR_REVIEW"
            eligibility = "REVIEW"
        elif not entity_ids and pass3_status not in {
            "UNRESOLVED_NO_CONTACT_EVIDENCE",
            "PDB_TARGET_HINT_REQUIRES_OCCURRENCE_REVIEW",
            "CONTACT_WITHOUT_ENTITY_METADATA_REVIEW",
        }:
            pass4_status = "NO_CONTACTING_ENTITY_IDS_REVIEW"
            eligibility = "REVIEW"

        rec.update({
            "pass4_selected_cif_path": str(selected) if selected else "",
            "pass4_cif_candidate_count": len(candidates),
            "pass4_cif_error": cif_error,
            "pass4_contact_entity_source_details": json.dumps(
                details, sort_keys=True, ensure_ascii=False
            ),
            "pass4_source_aggregate": aggregate,
            "pass4_status": pass4_status,
            "pass4_target_browser_eligible": eligibility,
        })
        output_rows.append(rec)

    out = pd.DataFrame(output_rows)
    out = out.sort_values([c_virus, c_pdb], kind="mergesort")

    main_path = args.outdir / "taxonomy_source_organism_qc_pass4.csv"
    out.to_csv(main_path, index=False)

    status_counts = Counter(out["pass4_status"].tolist())
    source_counts = Counter(out["pass4_source_aggregate"].tolist())
    eligibility_counts = Counter(out["pass4_target_browser_eligible"].tolist())

    confirmed_viral = int((out["pass4_status"] == "CONFIRMED_VIRAL_TARGET").sum())
    confirmed_nonviral = int(
        out["pass4_status"].isin(
            ["CONFIRMED_NONVIRAL_EXCLUSION", "MANUAL_CASE_CONFIRMED_NONVIRAL_CONTEXT"]
        ).sum()
    )

    summary = {
        "input_occurrences": int(len(out)),
        "distinct_structures": int(out[c_pdb].nunique()),
        "cif_files_indexed": int(sum(len(v) for v in cif_index.values())),
        "unique_pdbs_with_cif": int(len(cif_index)),
        "parsed_unique_pdbs": int(len(cif_cache)),
        "cif_parse_or_lookup_errors": int(
            out["pass4_status"].isin(["CIF_NOT_FOUND_REVIEW", "CIF_PARSE_ERROR_REVIEW"]).sum()
        ),
        "confirmed_viral_target_occurrences": confirmed_viral,
        "confirmed_nonviral_context_occurrences": confirmed_nonviral,
        "target_browser_eligibility_counts": dict(sorted(eligibility_counts.items())),
        "source_aggregate_counts": dict(sorted(source_counts.items())),
        "status_counts": dict(sorted(status_counts.items())),
        "production_data_modified": False,
        "policy": {
            "uses_gene_or_natural_source_as_biological_origin": True,
            "expression_host_used_as_target_organism": False,
            "synthetic_unknown_auto_confirmed": False,
            "nonmatching_populated_gene_or_natural_source_flags_conflict": True,
            "polyprotein_domain_auto_resolved": False,
        },
    }

    (args.outdir / "taxonomy_source_organism_qc_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    conflict_statuses = {
        "CONFLICT_PASS3_NONVIRAL_BUT_SOURCE_IS_VIRAL",
        "CONFLICT_POLYPROTEIN_ENTITY_SOURCE_NONVIRAL",
        "CONFLICT_TARGET_ASSIGNED_BUT_SOURCE_NONVIRAL",
        "CIF_NOT_FOUND_REVIEW",
        "CIF_PARSE_ERROR_REVIEW",
        "NO_CONTACTING_ENTITY_IDS_REVIEW",
    }
    conflicts = out[out["pass4_status"].isin(conflict_statuses)].copy()
    conflicts.to_csv(args.outdir / "taxonomy_source_organism_conflicts.csv", index=False)

    review = out[out["pass4_target_browser_eligible"] == "REVIEW"].copy()
    review.to_csv(args.outdir / "taxonomy_source_organism_review_queue.csv", index=False)

    eligible_cols = [
        c for c in [
            c_pdb,
            c_virus,
            c_target_id,
            c_target_name,
            c_family,
            c_pass3_role,
            c_pass3_status,
            "pass4_source_aggregate",
            "pass4_status",
            "pass4_target_browser_eligible",
        ] if c
    ]
    out[eligible_cols].to_csv(
        args.outdir / "taxonomy_target_browser_eligibility_pass4.csv",
        index=False,
    )

    print(json.dumps(summary, indent=2, sort_keys=True))
    print(f"\nWrote audit-only Pass-4 outputs to: {args.outdir}")


if __name__ == "__main__":
    main()
