#!/usr/bin/env python3
"""
V-LiSEMOD taxonomy reconciliation — Pass 5
Polyprotein mature-product/domain reconciliation (READ ONLY)

Purpose
-------
Resolve Pass-4 viral polyprotein cases by mapping the ligand-contacting mmCIF
entity to its database-reference alignment range and intersecting that range
with conservative, reviewed mature-product boundaries.

This pass does NOT modify:
- production SQLite
- Stage-09 / Stage-12 / Stage-14
- structure_classifications
- PROTACability scores
- source CIF/mmCIF files
- API/UI

It writes audit artifacts only.

Input
-----
--pass4-csv
    taxonomy_source_organism_qc_pass4.csv

Outputs
-------
taxonomy_polyprotein_domain_qc_pass5.csv
taxonomy_polyprotein_domain_qc_summary.json
taxonomy_polyprotein_domain_resolved.csv
taxonomy_polyprotein_domain_review_queue.csv
taxonomy_target_browser_eligibility_pass5.csv

Safety policy
-------------
- Only rows already confirmed or strongly supported as viral polyprotein contexts
  are considered.
- No folder-derived target label is used as the resolving evidence.
- No fuzzy sequence/name matching is used.
- A target is auto-resolved only when the mmCIF database-alignment interval(s)
  for the contacted entity map unambiguously to exactly one reviewed mature
  product/domain.
- Gag-Pol precursor coordinates are not auto-resolved unless a safe mapping can
  be established; they remain review.
- Missing/contradictory alignment data remain review.
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

NULLISH = {"", ".", "?", "none", "null", "n/a", "na", "unknown"}


def txt(value: Any) -> str:
    if value is None or pd.isna(value):
        return ""
    return re.sub(r"\s+", " ", str(value).strip())


def clean(value: Any) -> str:
    s = txt(value).strip("'\"")
    if s.casefold() in NULLISH:
        return ""
    return s


def as_list(value: Any) -> list[str]:
    if value is None:
        return []
    if isinstance(value, (list, tuple)):
        return [clean(x) for x in value]
    return [clean(value)]


def to_int(value: Any) -> int | None:
    s = clean(value)
    if not s:
        return None
    m = re.search(r"-?\d+", s)
    return int(m.group()) if m else None


def split_ids(value: Any) -> list[str]:
    s = txt(value)
    if not s:
        return []
    vals = [x.strip() for x in re.split(r"[;,|]", s) if x.strip()]
    out, seen = [], set()
    for v in vals:
        k = v.casefold()
        if k not in seen:
            seen.add(k)
            out.append(v)
    return out


def norm(value: Any) -> str:
    s = clean(value).casefold().replace("_", " ")
    s = re.sub(r"[^a-z0-9]+", " ", s)
    return re.sub(r"\s+", " ", s).strip()


def resolve_col(df: pd.DataFrame, candidates: list[str], required: bool = True) -> str | None:
    by = {c.casefold(): c for c in df.columns}
    for c in candidates:
        if c.casefold() in by:
            return by[c.casefold()]
    if required:
        raise KeyError(
            f"Missing required column; expected one of {candidates}. "
            f"Available: {list(df.columns)}"
        )
    return None


def col_values(cif: dict[str, Any], tag: str) -> list[str]:
    return as_list(cif.get(tag))


def value_at(vals: list[str], i: int) -> str:
    return vals[i] if i < len(vals) else ""


# ---------------------------------------------------------------------------
# Reviewed mature-product boundaries
#
# Coordinates are precursor-protein amino-acid coordinates and are used only
# when the database alignment is compatible with that precursor coordinate
# system.
# ---------------------------------------------------------------------------

SARS_COV2_ORF1AB = [
    # id, display name, family, start, end
    ("sars_cov_2_nsp1",  "Non-structural protein 1 (nsp1)",  "nsp_proteins", 1, 180),
    ("sars_cov_2_nsp2",  "Non-structural protein 2 (nsp2)",  "nsp_proteins", 181, 818),
    ("sars_cov_2_nsp3",  "Non-structural protein 3 (nsp3)",  "nsp_proteins", 819, 2763),
    ("sars_cov_2_nsp4",  "Non-structural protein 4 (nsp4)",  "nsp_proteins", 2764, 3263),
    ("sars_cov_2_nsp5",  "Main protease (nsp5)",             "protease",     3264, 3569),
    ("sars_cov_2_nsp6",  "Non-structural protein 6 (nsp6)",  "nsp_proteins", 3570, 3859),
    ("sars_cov_2_nsp7",  "Non-structural protein 7 (nsp7)",  "nsp_proteins", 3860, 3942),
    ("sars_cov_2_nsp8",  "Non-structural protein 8 (nsp8)",  "nsp_proteins", 3943, 4140),
    ("sars_cov_2_nsp9",  "Non-structural protein 9 (nsp9)",  "nsp_proteins", 4141, 4253),
    ("sars_cov_2_nsp10", "Non-structural protein 10 (nsp10)","nsp_proteins", 4254, 4392),
    ("sars_cov_2_nsp12", "RNA-dependent RNA polymerase (nsp12)", "polymerase", 4393, 5324),
    ("sars_cov_2_nsp13", "Helicase (nsp13)",                "helicase",      5325, 5925),
    ("sars_cov_2_nsp14", "Exoribonuclease/N7-methyltransferase (nsp14)", "nsp_proteins", 5926, 6452),
    ("sars_cov_2_nsp15", "Uridylate-specific endoribonuclease (nsp15)", "nsp_proteins", 6453, 6798),
    ("sars_cov_2_nsp16", "2'-O-methyltransferase (nsp16)",   "nsp_proteins", 6799, 7096),
]

SARS_COV2_ORF1A = [
    ("sars_cov_2_nsp1",  "Non-structural protein 1 (nsp1)",  "nsp_proteins", 1, 180),
    ("sars_cov_2_nsp2",  "Non-structural protein 2 (nsp2)",  "nsp_proteins", 181, 818),
    ("sars_cov_2_nsp3",  "Non-structural protein 3 (nsp3)",  "nsp_proteins", 819, 2763),
    ("sars_cov_2_nsp4",  "Non-structural protein 4 (nsp4)",  "nsp_proteins", 2764, 3263),
    ("sars_cov_2_nsp5",  "Main protease (nsp5)",             "protease",     3264, 3569),
    ("sars_cov_2_nsp6",  "Non-structural protein 6 (nsp6)",  "nsp_proteins", 3570, 3859),
    ("sars_cov_2_nsp7",  "Non-structural protein 7 (nsp7)",  "nsp_proteins", 3860, 3942),
    ("sars_cov_2_nsp8",  "Non-structural protein 8 (nsp8)",  "nsp_proteins", 3943, 4140),
    ("sars_cov_2_nsp9",  "Non-structural protein 9 (nsp9)",  "nsp_proteins", 4141, 4253),
    ("sars_cov_2_nsp10", "Non-structural protein 10 (nsp10)","nsp_proteins", 4254, 4392),
    ("sars_cov_2_nsp11", "Non-structural protein 11 (nsp11)","nsp_proteins", 4393, 4405),
]

# HIV-1 Gag precursor mature products (HXB2-like coordinate convention).
HIV1_GAG = [
    ("hiv_1_matrix",        "Matrix (p17)",        "matrix_protein", 1, 132),
    ("hiv_1_capsid",        "Capsid (p24)",        "capsid_protein", 133, 363),
    ("hiv_1_sp1",           "Spacer peptide 1",    "other_peptides_rna_complexes", 364, 377),
    ("hiv_1_nucleocapsid",  "Nucleocapsid (p7)",  "nucleoprotein", 378, 432),
    ("hiv_1_sp2",           "Spacer peptide 2",    "other_peptides_rna_complexes", 433, 448),
    ("hiv_1_p6",            "p6",                  "other_peptides_rna_complexes", 449, 500),
]

# HIV-1 Pol precursor coordinate convention used only when the aligned reference
# is clearly a Pol (not Gag-Pol) precursor. Splitting RT polymerase and RNase H
# is intentional because V-LiSEMOD already distinguishes rnase_h.
HIV1_POL = [
    ("hiv_1_protease",              "Protease",              "protease",               1, 99),
    ("hiv_1_reverse_transcriptase", "Reverse transcriptase", "reverse_transcriptase", 100, 440),
    ("hiv_1_rnase_h",               "RNase H",               "rnase_h",               441, 560),
    ("hiv_1_integrase",             "Integrase",              "integrase",             561, 848),
]


def interval_within(start: int, end: int, product: tuple) -> bool:
    _, _, _, pstart, pend = product
    return start >= pstart and end <= pend


def products_overlapping(start: int, end: int, products: list[tuple]) -> list[tuple]:
    out = []
    for p in products:
        _, _, _, a, b = p
        if max(start, a) <= min(end, b):
            out.append(p)
    return out


# ---------------------------------------------------------------------------
# mmCIF reference/alignment parsing
# ---------------------------------------------------------------------------

def parse_entity_descriptions(cif: dict[str, Any]) -> dict[str, str]:
    eids = col_values(cif, "_entity.id")
    descs = col_values(cif, "_entity.pdbx_description")
    return {eid: value_at(descs, i) for i, eid in enumerate(eids) if eid}


def parse_struct_refs(cif: dict[str, Any]) -> dict[str, dict[str, str]]:
    """
    _struct_ref.id -> metadata including entity_id/accession/db_code/db_name.
    """
    ids = col_values(cif, "_struct_ref.id")
    entity_ids = col_values(cif, "_struct_ref.entity_id")
    accessions = col_values(cif, "_struct_ref.pdbx_db_accession")
    db_codes = col_values(cif, "_struct_ref.db_code")
    db_names = col_values(cif, "_struct_ref.db_name")

    out = {}
    for i, rid in enumerate(ids):
        if not rid:
            continue
        out[rid] = {
            "entity_id": value_at(entity_ids, i),
            "accession": value_at(accessions, i),
            "db_code": value_at(db_codes, i),
            "db_name": value_at(db_names, i),
        }
    return out


def parse_ref_alignments(cif: dict[str, Any]) -> list[dict[str, Any]]:
    """
    Parse _struct_ref_seq rows. These are database-reference coordinate ranges.
    """
    ref_ids = col_values(cif, "_struct_ref_seq.ref_id")
    align_ids = col_values(cif, "_struct_ref_seq.align_id")
    db_beg = col_values(cif, "_struct_ref_seq.db_align_beg")
    db_end = col_values(cif, "_struct_ref_seq.db_align_end")
    seq_beg = col_values(cif, "_struct_ref_seq.seq_align_beg")
    seq_end = col_values(cif, "_struct_ref_seq.seq_align_end")
    strands = col_values(cif, "_struct_ref_seq.pdbx_strand_id")
    accessions = col_values(cif, "_struct_ref_seq.pdbx_db_accession")

    n = max(
        len(ref_ids), len(align_ids), len(db_beg), len(db_end),
        len(seq_beg), len(seq_end), len(strands), len(accessions), 0
    )
    rows = []
    for i in range(n):
        rows.append({
            "ref_id": value_at(ref_ids, i),
            "align_id": value_at(align_ids, i),
            "db_align_beg": to_int(value_at(db_beg, i)),
            "db_align_end": to_int(value_at(db_end, i)),
            "seq_align_beg": to_int(value_at(seq_beg, i)),
            "seq_align_end": to_int(value_at(seq_end, i)),
            "strand_ids": split_ids(value_at(strands, i)),
            "accession": value_at(accessions, i),
        })
    return rows


def parse_struct_asym(cif: dict[str, Any]) -> dict[str, str]:
    asym = col_values(cif, "_struct_asym.id")
    entity = col_values(cif, "_struct_asym.entity_id")
    return {a: value_at(entity, i) for i, a in enumerate(asym) if a}


def entity_alignment_records(
    cif: dict[str, Any],
    entity_id: str,
    contacting_chains: list[str],
) -> list[dict[str, Any]]:
    refs = parse_struct_refs(cif)
    aligns = parse_ref_alignments(cif)
    asym_to_entity = parse_struct_asym(cif)

    target_ref_ids = {
        rid for rid, meta in refs.items()
        if clean(meta.get("entity_id")) == clean(entity_id)
    }

    rows = []
    for a in aligns:
        via_ref = a["ref_id"] in target_ref_ids if a["ref_id"] else False

        via_chain = False
        if a["strand_ids"] and contacting_chains:
            for strand in a["strand_ids"]:
                # strand IDs may be author chains, so first direct compare;
                # also accept struct_asym -> entity if available.
                if strand in contacting_chains:
                    via_chain = True
                if asym_to_entity.get(strand) == entity_id:
                    via_chain = True

        if not via_ref and not via_chain:
            continue

        meta = refs.get(a["ref_id"], {})
        rec = dict(a)
        rec.update({
            "entity_id": entity_id,
            "ref_entity_id": meta.get("entity_id", ""),
            "ref_accession": a.get("accession") or meta.get("accession", ""),
            "ref_db_code": meta.get("db_code", ""),
            "ref_db_name": meta.get("db_name", ""),
        })
        rows.append(rec)

    # Deduplicate deterministic records.
    seen = set()
    out = []
    for r in rows:
        k = (
            r.get("ref_id"), r.get("db_align_beg"), r.get("db_align_end"),
            tuple(r.get("strand_ids", [])), r.get("ref_accession")
        )
        if k not in seen:
            seen.add(k)
            out.append(r)
    return out


# ---------------------------------------------------------------------------
# Domain resolution
# ---------------------------------------------------------------------------

def choose_coordinate_system(
    virus_name: str,
    entity_description: str,
    ref_meta_text: str,
    start: int,
    end: int,
) -> tuple[str, list[tuple]] | tuple[None, None]:
    v = norm(virus_name)
    desc = norm(entity_description)
    ref = norm(ref_meta_text)

    if "sars" in v and ("cov" in v or "coronavirus" in v):
        if "replicase polyprotein 1ab" in desc or "orf1ab" in desc or "1ab" in ref or "orf1ab" in ref:
            return "SARS_COV2_ORF1AB", SARS_COV2_ORF1AB
        if "replicase polyprotein 1a" in desc or "orf1a" in desc or "1a" in ref or "orf1a" in ref:
            return "SARS_COV2_ORF1A", SARS_COV2_ORF1A

        # Conservative fallback for clearly SARS-CoV-2 replicase precursor
        # coordinates: ranges >4405 necessarily imply ORF1ab.
        if end > 4405:
            return "SARS_COV2_ORF1AB", SARS_COV2_ORF1AB

    if ("hiv" in v and "1" in v) or "human immunodeficiency virus 1" in v:
        # "gag polyprotein" contains the character substring "gag pol", so
        # Gag-Pol detection must require a word boundary after the standalone
        # token "pol".
        gag_pol_desc = bool(re.search(r"\\bgag pol\\b", desc)) or "pr160" in desc
        gag_pol_ref = bool(re.search(r"\\bgag pol\\b", ref)) or "pr160" in ref

        # Do not auto-map true Gag-Pol coordinates from the simple Gag/Pol tables.
        if gag_pol_desc or gag_pol_ref:
            return None, None

        if "gag polyprotein" in desc or desc == "gag protein" or re.search(r"\\bgag\\b", ref):
            return "HIV1_GAG", HIV1_GAG

        if desc in {"pol protein", "pol"} or "pol protein" in desc:
            if not gag_pol_ref and end <= 1000:
                return "HIV1_POL", HIV1_POL

    return None, None


def resolve_records(
    virus_name: str,
    entity_description: str,
    records: list[dict[str, Any]],
) -> dict[str, Any]:
    """
    Auto-resolve only if all usable alignment intervals independently resolve
    to the same mature product.
    """
    usable = []
    mapped = []

    for r in records:
        start = r.get("db_align_beg")
        end = r.get("db_align_end")
        if start is None or end is None:
            continue
        if start > end:
            start, end = end, start

        ref_text = " ".join([
            clean(r.get("ref_accession")),
            clean(r.get("ref_db_code")),
            clean(r.get("ref_db_name")),
        ])
        coord_name, products = choose_coordinate_system(
            virus_name, entity_description, ref_text, start, end
        )
        if not products:
            usable.append({
                "range": [start, end],
                "coordinate_system": "",
                "status": "COORDINATE_SYSTEM_UNRESOLVED",
                "ref": ref_text,
            })
            continue

        overlaps = products_overlapping(start, end, products)
        contained = [p for p in products if interval_within(start, end, p)]

        if len(contained) == 1:
            p = contained[0]
            usable.append({
                "range": [start, end],
                "coordinate_system": coord_name,
                "status": "RANGE_WITHIN_ONE_PRODUCT",
                "product_id": p[0],
                "product_name": p[1],
                "target_family": p[2],
                "product_range": [p[3], p[4]],
                "ref": ref_text,
            })
            mapped.append(p)
        else:
            usable.append({
                "range": [start, end],
                "coordinate_system": coord_name,
                "status": "RANGE_SPANS_OR_MISSES_PRODUCTS",
                "overlap_products": [p[0] for p in overlaps],
                "ref": ref_text,
            })

    if not usable:
        return {
            "status": "NO_USABLE_DB_ALIGNMENT",
            "evidence": [],
        }

    if not mapped:
        return {
            "status": "ALIGNMENT_NOT_UNAMBIGUOUS",
            "evidence": usable,
        }

    unique = {}
    for p in mapped:
        unique[p[0]] = p

    # Every usable alignment record must map cleanly; otherwise keep review.
    all_clean = all(x.get("status") == "RANGE_WITHIN_ONE_PRODUCT" for x in usable)

    if all_clean and len(unique) == 1:
        p = next(iter(unique.values()))
        return {
            "status": "RESOLVED_UNAMBIGUOUS",
            "canonical_target_id": p[0],
            "canonical_target_name": p[1],
            "target_family": p[2],
            "evidence": usable,
        }

    if len(unique) > 1:
        return {
            "status": "MULTIPLE_MATURE_PRODUCTS_PRESENT",
            "evidence": usable,
        }

    return {
        "status": "ALIGNMENT_PARTIALLY_AMBIGUOUS",
        "evidence": usable,
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

ELIGIBLE_PASS4_STATUSES = {
    "VIRAL_POLYPROTEIN_CONFIRMED_DOMAIN_UNRESOLVED",
    "POLYPROTEIN_SOURCE_OR_DOMAIN_REVIEW",
}

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--pass4-csv",
        required=True,
        type=Path,
        help="taxonomy_source_organism_qc_pass4.csv",
    )
    ap.add_argument(
        "--outdir",
        required=True,
        type=Path,
    )
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.pass4_csv, dtype=str, low_memory=False).fillna("")

    c_pdb = resolve_col(df, ["pdb_id", "pdb"])
    c_virus = resolve_col(df, ["virus_name", "virus"])
    c_pass4 = resolve_col(df, ["pass4_status"])
    c_cif = resolve_col(df, ["pass4_selected_cif_path"])
    c_entity_ids = resolve_col(
        df, ["contacting_entity_ids", "contact_entity_ids", "entity_ids"], required=False
    )
    c_entity_desc = resolve_col(
        df, ["contacting_entity_descriptions", "entity_descriptions"], required=False
    )
    c_chains = resolve_col(
        df, ["contacting_protein_chains", "target_chains", "protein_contact_chains"],
        required=False,
    )

    rows = []
    cif_cache: dict[str, dict[str, Any]] = {}

    for _, r in df.iterrows():
        rec = dict(r)
        p4 = clean(r[c_pass4])

        rec.update({
            "pass5_status": "NOT_APPLICABLE",
            "pass5_canonical_target_id": "",
            "pass5_canonical_target_name": "",
            "pass5_target_family": "",
            "pass5_coordinate_system": "",
            "pass5_alignment_evidence": "",
            "pass5_note": "",
            "pass5_target_browser_eligible": clean(r.get("pass4_target_browser_eligible", "")),
        })

        if p4 not in ELIGIBLE_PASS4_STATUSES:
            rows.append(rec)
            continue

        pdb = clean(r[c_pdb]).upper()
        virus = clean(r[c_virus])
        cif_path = Path(clean(r[c_cif]))
        entity_ids = split_ids(r[c_entity_ids]) if c_entity_ids else []
        chains = split_ids(r[c_chains]) if c_chains else []
        desc_raw = clean(r[c_entity_desc]) if c_entity_desc else ""

        if not cif_path.exists():
            rec.update({
                "pass5_status": "CIF_NOT_FOUND_REVIEW",
                "pass5_note": str(cif_path),
                "pass5_target_browser_eligible": "REVIEW",
            })
            rows.append(rec)
            continue

        cache_key = str(cif_path.resolve())
        if cache_key not in cif_cache:
            try:
                cif_cache[cache_key] = MMCIF2Dict(str(cif_path))
            except Exception as exc:
                cif_cache[cache_key] = {"__ERROR__": f"{type(exc).__name__}: {exc}"}

        cif = cif_cache[cache_key]
        if "__ERROR__" in cif:
            rec.update({
                "pass5_status": "CIF_PARSE_ERROR_REVIEW",
                "pass5_note": cif["__ERROR__"],
                "pass5_target_browser_eligible": "REVIEW",
            })
            rows.append(rec)
            continue

        descriptions = parse_entity_descriptions(cif)

        if not entity_ids:
            rec.update({
                "pass5_status": "NO_CONTACTING_ENTITY_ID_REVIEW",
                "pass5_target_browser_eligible": "REVIEW",
            })
            rows.append(rec)
            continue

        entity_results = []
        for eid in entity_ids:
            desc = clean(descriptions.get(eid)) or desc_raw
            alignment_records = entity_alignment_records(cif, eid, chains)
            result = resolve_records(virus, desc, alignment_records)
            result["entity_id"] = eid
            result["entity_description"] = desc
            result["alignment_records"] = alignment_records
            entity_results.append(result)

        resolved = [
            x for x in entity_results
            if x.get("status") == "RESOLVED_UNAMBIGUOUS"
        ]

        # Require every contacted entity in the polyprotein case to resolve,
        # and require all of them to the same mature product.
        if len(resolved) == len(entity_results) and resolved:
            target_ids = {x["canonical_target_id"] for x in resolved}
            if len(target_ids) == 1:
                x = resolved[0]
                coord_systems = set()
                for er in resolved:
                    for ev in er.get("evidence", []):
                        if ev.get("coordinate_system"):
                            coord_systems.add(ev["coordinate_system"])
                rec.update({
                    "pass5_status": "RESOLVED_POLYPROTEIN_MATURE_PRODUCT",
                    "pass5_canonical_target_id": x["canonical_target_id"],
                    "pass5_canonical_target_name": x["canonical_target_name"],
                    "pass5_target_family": x["target_family"],
                    "pass5_coordinate_system": ";".join(sorted(coord_systems)),
                    "pass5_alignment_evidence": json.dumps(
                        entity_results, ensure_ascii=False, sort_keys=True
                    ),
                    "pass5_target_browser_eligible": "YES",
                })
            else:
                rec.update({
                    "pass5_status": "MULTIPLE_MATURE_PRODUCTS_REVIEW",
                    "pass5_alignment_evidence": json.dumps(
                        entity_results, ensure_ascii=False, sort_keys=True
                    ),
                    "pass5_target_browser_eligible": "REVIEW",
                })
        else:
            statuses = sorted({x.get("status", "") for x in entity_results})
            # Explicitly distinguish Gag-Pol/non-mappable coordinate cases.
            if statuses == ["ALIGNMENT_NOT_UNAMBIGUOUS"] or \
               "NO_USABLE_DB_ALIGNMENT" in statuses:
                p5 = "POLYPROTEIN_ALIGNMENT_REVIEW"
            else:
                p5 = "POLYPROTEIN_DOMAIN_REVIEW_REMAINS"

            rec.update({
                "pass5_status": p5,
                "pass5_alignment_evidence": json.dumps(
                    entity_results, ensure_ascii=False, sort_keys=True
                ),
                "pass5_target_browser_eligible": "REVIEW",
            })

        rows.append(rec)

    out = pd.DataFrame(rows)
    out = out.sort_values([c_virus, c_pdb], kind="mergesort")

    out.to_csv(args.outdir / "taxonomy_polyprotein_domain_qc_pass5.csv", index=False)

    applicable = out[out["pass5_status"] != "NOT_APPLICABLE"].copy()
    resolved = out[
        out["pass5_status"] == "RESOLVED_POLYPROTEIN_MATURE_PRODUCT"
    ].copy()
    review = applicable[
        applicable["pass5_status"] != "RESOLVED_POLYPROTEIN_MATURE_PRODUCT"
    ].copy()

    resolved.to_csv(args.outdir / "taxonomy_polyprotein_domain_resolved.csv", index=False)
    review.to_csv(args.outdir / "taxonomy_polyprotein_domain_review_queue.csv", index=False)

    # Final eligibility view: Pass 5 overrides Pass 4 only for applicable rows.
    eligibility = out.copy()
    eligibility.to_csv(
        args.outdir / "taxonomy_target_browser_eligibility_pass5.csv",
        index=False,
    )

    status_counts = Counter(applicable["pass5_status"].tolist())
    target_counts = Counter(
        x for x in resolved["pass5_canonical_target_name"].tolist() if clean(x)
    )

    summary = {
        "input_occurrences": int(len(out)),
        "pass5_applicable_polyprotein_occurrences": int(len(applicable)),
        "resolved_polyprotein_occurrences": int(len(resolved)),
        "remaining_polyprotein_review_occurrences": int(len(review)),
        "status_counts": dict(sorted(status_counts.items())),
        "resolved_target_counts": dict(sorted(target_counts.items())),
        "production_data_modified": False,
        "policy": {
            "folder_label_used_as_resolving_evidence": False,
            "fuzzy_matching": False,
            "requires_unambiguous_db_alignment_to_single_mature_product": True,
            "gag_pol_precursor_auto_resolved": False,
            "scientific_scores_modified": False,
        },
    }

    (args.outdir / "taxonomy_polyprotein_domain_qc_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    print(json.dumps(summary, indent=2, sort_keys=True))
    print(f"\nWrote audit-only Pass-5 outputs to: {args.outdir}")


if __name__ == "__main__":
    main()
