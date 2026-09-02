#!/usr/bin/env python3
"""
V-LiSEMOD taxonomy reconciliation pass 3 (role-aware, audit-only).

Input:
  Pass-2 classification_reconciliation_audit_taxonomy_pass.csv
  reviewed_taxonomy_seed_v3.csv

Output:
  taxonomy_role_reconciliation_pass3.csv
  taxonomy_role_reconciliation_summary.json
  taxonomy_role_manual_review_queue.csv

This does NOT modify the production database or scientific scores.

Important design change:
  An entity can be a VIRAL_TARGET, HOST_PARTNER, ENGINEERED_COMPONENT,
  EXTERNAL_BINDER, NONPOLYMER_ARTIFACT, or POLYPROTEIN_AMBIGUOUS.
  Host/binder/engineered entities are never silently converted to a viral target.
"""

from __future__ import annotations
import argparse, csv, json, re
from collections import Counter
from pathlib import Path
import pandas as pd

def txt(x):
    if x is None or pd.isna(x):
        return ""
    return re.sub(r"\s+", " ", str(x).strip())

def key(x):
    s = txt(x).casefold().replace("_", " ")
    # Remove only wrapping quote noise; preserve meaningful internal punctuation.
    s = s.strip(" '\"")
    s = re.sub(r"\s+", " ", s)
    return s

def split_atomic_descriptions(s):
    """
    Split semicolon-delimited descriptions only at top level.
    Do not split semicolons inside quotes or (), [], {}.
    """
    s = txt(s)
    if not s:
        return []
    out, buf = [], []
    quote = None
    depth = 0
    pairs_open = "([{"
    pairs_close = ")]}"
    for ch in s:
        if quote:
            buf.append(ch)
            if ch == quote:
                quote = None
            continue
        if ch in ("'", '"'):
            quote = ch
            buf.append(ch)
            continue
        if ch in pairs_open:
            depth += 1
            buf.append(ch)
            continue
        if ch in pairs_close:
            depth = max(0, depth - 1)
            buf.append(ch)
            continue
        if ch == ";" and depth == 0:
            part = "".join(buf).strip()
            if part:
                out.append(part.strip(" '\""))
            buf = []
        else:
            buf.append(ch)
    part = "".join(buf).strip()
    if part:
        out.append(part.strip(" '\""))
    return [x for x in out if x]

def resolve_col(df, names, required=True):
    by = {c.casefold(): c for c in df.columns}
    for n in names:
        if n.casefold() in by:
            return by[n.casefold()]
    if required:
        raise KeyError(f"Missing any of {names}; columns={list(df.columns)}")
    return None

def load_seed(path):
    sdf = pd.read_csv(path, dtype=str).fillna("")
    out = {}
    for _, r in sdf.iterrows():
        out[(txt(r["virus_name"]).casefold(), key(r["entity_description"]))] = {
            "canonical_target_id": txt(r["canonical_target_id"]),
            "canonical_target_name": txt(r["canonical_target_name"]),
            "target_family": txt(r["target_family"]),
            "entity_role": txt(r["entity_role"]),
            "note": txt(r["note"]),
        }
    return out

# The 248 "contact without entity metadata" rows need a separate structural
# metadata pass. These hints are verified target-context hints, NOT automatic
# occurrence-level resolution.
PDB_TARGET_HINTS = {
    # gp41 loop analogue structures from one RCSB-linked study
    "1J8N": ("hiv_1_gp41", "gp41", "envelope_glycoprotein"),
    "1J8Z": ("hiv_1_gp41", "gp41", "envelope_glycoprotein"),
    "1J9V": ("hiv_1_gp41", "gp41", "envelope_glycoprotein"),
    "1JAA": ("hiv_1_gp41", "gp41", "envelope_glycoprotein"),
    "1JAR": ("hiv_1_gp41", "gp41", "envelope_glycoprotein"),
    "1JC8": ("hiv_1_gp41", "gp41", "envelope_glycoprotein"),
    "1JCP": ("hiv_1_gp41", "gp41", "envelope_glycoprotein"),
    "1RPV": ("hiv_1_rev", "Rev", "rev_protein"),
}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pass2-csv", type=Path, required=True)
    ap.add_argument("--seed-csv", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.pass2_csv, dtype=str, low_memory=False).fillna("")
    seed = load_seed(args.seed_csv)

    cpdb = resolve_col(df, ["pdb_id", "pdb"])
    cvirus = resolve_col(df, ["virus_name", "virus"])
    cdesc = resolve_col(df, ["contacting_entity_descriptions", "entity_descriptions"])
    cstatus = resolve_col(df, ["taxonomy_reconciliation_status"])
    ctarget = resolve_col(df, ["taxonomy_reconciled_target"], required=False)
    ccurrent = resolve_col(df, ["current_stage14_protein_type", "current_target_label", "protein_type"], required=False)

    rows = []
    for _, r in df.iterrows():
        pdb = txt(r[cpdb]).upper()
        virus = txt(r[cvirus])
        status2 = txt(r[cstatus])
        desc_raw = txt(r[cdesc])
        existing = txt(r[ctarget]) if ctarget else ""

        rec = dict(r)
        rec.update({
            "pass3_canonical_target_id": "",
            "pass3_canonical_target_name": "",
            "pass3_target_family": "",
            "pass3_entity_role": "",
            "pass3_status": "",
            "pass3_mapping_evidence": "",
            "pass3_unmapped_descriptions": "",
        })

        # Keep pass-2 exact resolution, but label it provisional until the
        # separate organism/source-entity QC confirms that it is viral.
        if status2 == "AUTO_RESOLVED_EXACT" and existing:
            rec.update({
                "pass3_canonical_target_id": existing,
                "pass3_canonical_target_name": existing,
                "pass3_entity_role": "VIRAL_TARGET_PROVISIONAL",
                "pass3_status": "KEEP_PASS2_EXACT_PENDING_ENTITY_SOURCE_QC",
                "pass3_mapping_evidence": "Pass-2 exact reconciliation",
            })
            rows.append(rec)
            continue

        if status2 == "UNRESOLVED_NO_CONTACT_EVIDENCE":
            rec.update({
                "pass3_status": "UNRESOLVED_NO_CONTACT_EVIDENCE",
                "pass3_entity_role": "UNRESOLVED",
            })
            rows.append(rec)
            continue

        if status2 == "MANUAL_REVIEW_CONTACT_WITHOUT_ENTITY_METADATA":
            if pdb in PDB_TARGET_HINTS:
                tid, tname, fam = PDB_TARGET_HINTS[pdb]
                rec.update({
                    "pass3_canonical_target_id": tid,
                    "pass3_canonical_target_name": tname,
                    "pass3_target_family": fam,
                    "pass3_entity_role": "VIRAL_TARGET_HINT",
                    "pass3_status": "PDB_TARGET_HINT_REQUIRES_OCCURRENCE_REVIEW",
                    "pass3_mapping_evidence": "PDB-level verified target context; entity-level metadata absent",
                })
            else:
                rec.update({
                    "pass3_entity_role": "UNRESOLVED",
                    "pass3_status": "CONTACT_WITHOUT_ENTITY_METADATA_REVIEW",
                })
            rows.append(rec)
            continue

        # First try the full atomic description as a curated seed key.
        full = seed.get((virus.casefold(), key(desc_raw)))
        if full:
            mappings = [(desc_raw, full)]
            unmapped = []
        else:
            mappings, unmapped = [], []
            for d in split_atomic_descriptions(desc_raw):
                m = seed.get((virus.casefold(), key(d)))
                if m:
                    mappings.append((d, m))
                else:
                    unmapped.append(d)

        viral = []
        nonviral_roles = []
        ambiguous = []
        evidence = []
        for d, m in mappings:
            role = m["entity_role"]
            evidence.append(f"{d}=>{role}:{m['canonical_target_id'] or m['canonical_target_name']}")
            if role.startswith("VIRAL_TARGET"):
                viral.append(m)
            elif role == "POLYPROTEIN_AMBIGUOUS":
                ambiguous.append(m)
            else:
                nonviral_roles.append(role)

        # De-duplicate viral target IDs.
        uniq_viral = {}
        for m in viral:
            tid = m["canonical_target_id"]
            if tid:
                uniq_viral[tid] = m

        if ambiguous and not uniq_viral:
            rec.update({
                "pass3_entity_role": "POLYPROTEIN_AMBIGUOUS",
                "pass3_status": "POLYPROTEIN_DOMAIN_REVIEW",
                "pass3_mapping_evidence": " | ".join(evidence),
                "pass3_unmapped_descriptions": ";".join(unmapped),
            })
        elif len(uniq_viral) > 1:
            rec.update({
                "pass3_entity_role": "MULTIPLE_VIRAL_TARGETS",
                "pass3_status": "MULTIVIRAL_INTERFACE_REVIEW",
                "pass3_mapping_evidence": " | ".join(evidence),
                "pass3_unmapped_descriptions": ";".join(unmapped),
            })
        elif len(uniq_viral) == 1:
            m = next(iter(uniq_viral.values()))
            rec.update({
                "pass3_canonical_target_id": m["canonical_target_id"],
                "pass3_canonical_target_name": m["canonical_target_name"],
                "pass3_target_family": m["target_family"],
                "pass3_mapping_evidence": " | ".join(evidence),
                "pass3_unmapped_descriptions": ";".join(unmapped),
            })
            if unmapped:
                rec["pass3_entity_role"] = "VIRAL_TARGET_WITH_UNRESOLVED_COMPONENT"
                rec["pass3_status"] = "VIRAL_TARGET_COMPONENT_REVIEW"
            elif any(x in ("HOST_PARTNER", "EXTERNAL_BINDER") for x in nonviral_roles):
                rec["pass3_entity_role"] = "VIRAL_TARGET_WITH_NONVIRAL_PARTNER"
                rec["pass3_status"] = "VIRAL_NONVIRAL_INTERFACE_REVIEW"
            elif any(x == "ENGINEERED_COMPONENT" for x in nonviral_roles):
                rec["pass3_entity_role"] = "VIRAL_TARGET_WITH_ENGINEERED_COMPONENT"
                rec["pass3_status"] = "RESOLVED_CURATED_WITH_ENGINEERED_COMPONENT"
            elif any(x == "NONPOLYMER_ARTIFACT" for x in nonviral_roles):
                rec["pass3_entity_role"] = "VIRAL_TARGET_WITH_NONPOLYMER_ARTIFACT"
                rec["pass3_status"] = "RESOLVED_TARGET_BUT_AUDIT_ENTITY_MAPPING_REVIEW"
            else:
                rec["pass3_entity_role"] = m["entity_role"]
                rec["pass3_status"] = "RESOLVED_CURATED_EXACT"
        elif mappings and not unmapped and nonviral_roles:
            # Contact context contains no viral target at all.
            roles = sorted(set(nonviral_roles))
            rec.update({
                "pass3_entity_role": ";".join(roles),
                "pass3_status": "EXCLUDE_NONVIRAL_CONTACT_CONTEXT",
                "pass3_mapping_evidence": " | ".join(evidence),
            })
        else:
            rec.update({
                "pass3_entity_role": "UNRESOLVED",
                "pass3_status": "MANUAL_REVIEW_AFTER_ROLE_PASS",
                "pass3_mapping_evidence": " | ".join(evidence),
                "pass3_unmapped_descriptions": ";".join(unmapped),
            })

        rows.append(rec)

    out = pd.DataFrame(rows)
    sort_cols = [cvirus, cpdb]
    out = out.sort_values(sort_cols, kind="mergesort")
    out_path = args.outdir / "taxonomy_role_reconciliation_pass3.csv"
    out.to_csv(out_path, index=False)

    counts = Counter(out["pass3_status"])
    role_counts = Counter(out["pass3_entity_role"])
    summary = {
        "input_occurrences": int(len(out)),
        "status_counts": dict(sorted(counts.items())),
        "role_counts": dict(sorted(role_counts.items())),
        "resolved_curated_exact": int((out["pass3_status"] == "RESOLVED_CURATED_EXACT").sum()),
        "exclude_nonviral_contact_context": int((out["pass3_status"] == "EXCLUDE_NONVIRAL_CONTACT_CONTEXT").sum()),
        "polyprotein_domain_review": int((out["pass3_status"] == "POLYPROTEIN_DOMAIN_REVIEW").sum()),
        "pdb_target_hint_requires_occurrence_review": int((out["pass3_status"] == "PDB_TARGET_HINT_REQUIRES_OCCURRENCE_REVIEW").sum()),
        "production_data_modified": False,
    }
    (args.outdir / "taxonomy_role_reconciliation_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    review_statuses = {
        "MANUAL_REVIEW_AFTER_ROLE_PASS",
        "VIRAL_TARGET_COMPONENT_REVIEW",
        "VIRAL_NONVIRAL_INTERFACE_REVIEW",
        "MULTIVIRAL_INTERFACE_REVIEW",
        "POLYPROTEIN_DOMAIN_REVIEW",
        "PDB_TARGET_HINT_REQUIRES_OCCURRENCE_REVIEW",
        "CONTACT_WITHOUT_ENTITY_METADATA_REVIEW",
        "RESOLVED_TARGET_BUT_AUDIT_ENTITY_MAPPING_REVIEW",
    }
    review = out[out["pass3_status"].isin(review_statuses)].copy()
    review.to_csv(args.outdir / "taxonomy_role_manual_review_queue.csv", index=False)

    print(json.dumps(summary, indent=2, sort_keys=True))
    print(f"\nWrote: {out_path}")

if __name__ == "__main__":
    main()
