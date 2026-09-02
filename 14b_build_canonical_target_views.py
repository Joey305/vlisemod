#!/usr/bin/env python3
"""Build canonical occurrence-level Target Browser views for V-LiSEMOD.

This is intentionally separate from the frozen Stage-14 compatibility views.
It does not rewrite scientific outputs or structure_classifications.

Required upstream table:
    canonical_ligand_targets

Created views:
    v2_target_browser_ligand_context
    v2_target_browser_groups

Only target_browser_eligible='YES' occurrences are exposed to the Target Browser.
Folder-derived Stage-14 labels are retained as source_* provenance fields.
"""
from __future__ import annotations

import argparse
import json
import sqlite3
from pathlib import Path

VERSION = "canonical-target-browser-views-pass9-v1"
TARGET_AUTHORITY_VERSION = "canonical-target-authority-pass9-v1"

EXPECTED_TOTAL = 7355
EXPECTED_ELIGIBILITY = {"NO": 796, "REVIEW": 1145, "YES": 5414}
EXPECTED_GROUPS = 32


def connect(database: Path) -> sqlite3.Connection:
    db = sqlite3.connect(str(database), timeout=60.0)
    db.row_factory = sqlite3.Row
    db.execute("PRAGMA foreign_keys=ON")
    db.execute("PRAGMA busy_timeout=60000")
    return db


def object_type(db, name: str):
    r = db.execute(
        "SELECT type FROM sqlite_master WHERE lower(name)=lower(?)", (name,)
    ).fetchone()
    return r[0] if r else None


def replace_view(db, name: str, sql: str):
    typ = object_type(db, name)
    if typ == "table":
        raise RuntimeError(f"{name} exists as a table; refusing to replace it.")
    if typ == "view":
        db.execute(f'DROP VIEW "{name}"')
    db.execute(f'CREATE VIEW "{name}" AS {sql}')


def authority_validation(db):
    if object_type(db, "canonical_ligand_targets") != "table":
        raise RuntimeError(
            "canonical_ligand_targets is missing. Run the validated Pass-9 loader first."
        )
    if object_type(db, "v2_ligand_context") != "view":
        raise RuntimeError(
            "v2_ligand_context is missing. Run the existing Stage 14 first."
        )

    total = db.execute(
        "SELECT count(*) FROM canonical_ligand_targets"
    ).fetchone()[0]
    if total != EXPECTED_TOTAL:
        raise RuntimeError(f"authority rows={total}; expected {EXPECTED_TOTAL}")

    eligibility = {
        str(r[0]): int(r[1])
        for r in db.execute(
            """SELECT target_browser_eligible,count(*)
               FROM canonical_ligand_targets
               GROUP BY target_browser_eligible"""
        )
    }
    if eligibility != EXPECTED_ELIGIBILITY:
        raise RuntimeError(
            f"authority eligibility={eligibility}; expected={EXPECTED_ELIGIBILITY}"
        )

    versions = {
        str(r[0]): int(r[1])
        for r in db.execute(
            """SELECT authority_version,count(*)
               FROM canonical_ligand_targets
               GROUP BY authority_version"""
        )
    }
    if versions != {TARGET_AUTHORITY_VERSION: EXPECTED_TOTAL}:
        raise RuntimeError(
            f"authority version distribution={versions}; "
            f"expected only {TARGET_AUTHORITY_VERSION}"
        )

    missing = db.execute(
        """SELECT count(*)
           FROM ligand_instances i
           LEFT JOIN canonical_ligand_targets c
             ON c.ligand_instance_id=i.ligand_instance_id
           WHERE i.curation_status='included'
             AND c.ligand_instance_id IS NULL"""
    ).fetchone()[0]
    if missing:
        raise RuntimeError(f"{missing} included ligand occurrences lack authority rows")

    bad_yes = db.execute(
        """SELECT count(*)
           FROM canonical_ligand_targets
           WHERE target_browser_eligible='YES'
             AND (
                canonical_target_id IS NULL
                OR trim(canonical_target_id)=''
                OR canonical_target_name IS NULL
                OR trim(canonical_target_name)=''
                OR canonical_target_id LIKE '%,%'
                OR canonical_target_id LIKE '%;%'
             )"""
    ).fetchone()[0]
    if bad_yes:
        raise RuntimeError(f"{bad_yes} authoritative YES rows have invalid target IDs")

    metadata_conflicts = db.execute(
        """SELECT count(*) FROM (
             SELECT virus_name, canonical_target_id
             FROM canonical_ligand_targets
             WHERE target_browser_eligible='YES'
             GROUP BY virus_name, canonical_target_id
             HAVING count(DISTINCT canonical_target_name) > 1
                OR count(DISTINCT coalesce(target_family, '')) > 1
                OR count(DISTINCT entity_role) > 1
           )"""
    ).fetchone()[0]
    if metadata_conflicts:
        raise RuntimeError(
            f"{metadata_conflicts} canonical target identities have inconsistent display metadata"
        )

    group_count = db.execute(
        """SELECT count(*) FROM (
             SELECT virus_name, canonical_target_id
             FROM canonical_ligand_targets
             WHERE target_browser_eligible='YES'
             GROUP BY virus_name, canonical_target_id
           )"""
    ).fetchone()[0]
    if group_count != EXPECTED_GROUPS:
        raise RuntimeError(
            f"canonical group count={group_count}; expected={EXPECTED_GROUPS}"
        )

    reg = db.execute(
        """SELECT c.ligand_instance_id,c.canonical_target_id,
                  c.canonical_target_name,c.target_browser_eligible
           FROM canonical_ligand_targets c
           JOIN ligand_instances i
             ON i.ligand_instance_id=c.ligand_instance_id
           JOIN structures s ON s.structure_id=i.structure_id
           WHERE upper(s.entry_id)='2O4K'
             AND upper(i.label_comp_id)='DR7'
             AND i.curation_status='included'"""
    ).fetchall()

    if len(reg) != 1:
        raise RuntimeError(f"2O4K/DR7 expected one authority row; observed {len(reg)}")

    r = reg[0]
    observed_reg = (
        int(r["ligand_instance_id"]),
        str(r["canonical_target_id"]),
        str(r["canonical_target_name"]),
        str(r["target_browser_eligible"]),
    )
    expected_reg = (36170, "protease", "protease", "YES")
    if observed_reg != expected_reg:
        raise RuntimeError(
            f"2O4K/DR7 regression failed: observed={observed_reg}, "
            f"expected={expected_reg}"
        )

    return {
        "authority_rows": total,
        "eligibility_counts": eligibility,
        "canonical_groups": group_count,
        "authority_version": TARGET_AUTHORITY_VERSION,
    }


def protected_counts(db):
    names = [
        "structures",
        "structure_classifications",
        "ligand_instances",
        "protacability_target_context",
        "protacability_assessment",
        "protacability_degrader_readiness",
        "protacability_attachment_sites",
        "protacability_attachment_site_summary",
    ]
    out = {}
    for name in names:
        if object_type(db, name) == "table":
            out[name] = db.execute(f'SELECT count(*) FROM "{name}"').fetchone()[0]
    return out


def run(database: Path):
    db = connect(database)
    try:
        validation = authority_validation(db)
        before = protected_counts(db)

        db.execute("BEGIN IMMEDIATE")

        replace_view(
            db,
            "v2_target_browser_ligand_context",
            """
            SELECT
                c.virus_name AS virus_name,
                c.canonical_target_name AS protein_type,
                c.canonical_target_id,
                c.canonical_target_name,
                c.target_family,
                c.entity_role,
                c.decision AS target_authority_decision,
                c.authority_basis AS target_authority_basis,
                c.authority_note AS target_authority_note,
                c.authority_version AS target_authority_version,

                lc.virus_name AS source_virus_name,
                lc.protein_type AS source_protein_type,

                lc.pdb_id,
                lc.ligand_instance_id,
                lc.model_id,
                lc.ligand,
                lc.chain,
                lc.ligand_residue_id,
                lc.ligand_insertion_code,
                lc.smiles,
                lc.canonical_smiles,
                lc.smiles_source,
                lc.chemical_status,
                lc.curation_status,
                lc.curation_reason

            FROM v2_ligand_context lc
            JOIN canonical_ligand_targets c
              ON c.ligand_instance_id=lc.ligand_instance_id

            WHERE lc.curation_status='included'
              AND c.target_browser_eligible='YES'
            """
        )

        replace_view(
            db,
            "v2_target_browser_groups",
            """
            SELECT
                virus_name,
                canonical_target_id,
                min(canonical_target_name) AS canonical_target_name,
                min(protein_type) AS protein_type,
                min(target_family) AS target_family,
                min(entity_role) AS entity_role,
                count(*) AS occurrence_count,
                count(DISTINCT pdb_id) AS structure_count,
                count(DISTINCT ligand) AS ligand_count

            FROM v2_target_browser_ligand_context

            GROUP BY
                virus_name,
                canonical_target_id
            """
        )

        occurrence_rows = db.execute(
            "SELECT count(*) FROM v2_target_browser_ligand_context"
        ).fetchone()[0]
        group_rows = db.execute(
            "SELECT count(*) FROM v2_target_browser_groups"
        ).fetchone()[0]

        if occurrence_rows != EXPECTED_ELIGIBILITY["YES"]:
            raise RuntimeError(
                f"Target Browser occurrence view has {occurrence_rows} rows; "
                f"expected {EXPECTED_ELIGIBILITY['YES']}"
            )
        if group_rows != EXPECTED_GROUPS:
            raise RuntimeError(
                f"Target Browser group view has {group_rows} rows; expected {EXPECTED_GROUPS}"
            )

        reg = db.execute(
            """SELECT ligand_instance_id,virus_name,protein_type,
                      canonical_target_id,source_protein_type
               FROM v2_target_browser_ligand_context
               WHERE upper(pdb_id)='2O4K'
                 AND upper(ligand)='DR7'"""
        ).fetchall()
        if len(reg) != 1:
            raise RuntimeError(
                f"Target Browser 2O4K/DR7 expected one row; observed {len(reg)}"
            )
        rr = reg[0]
        if (
            int(rr["ligand_instance_id"]) != 36170
            or str(rr["canonical_target_id"]).casefold() != "protease"
            or str(rr["protein_type"]).casefold() != "protease"
        ):
            raise RuntimeError(f"Target Browser 2O4K/DR7 regression failed: {dict(rr)}")

        after = protected_counts(db)
        if before != after:
            raise RuntimeError(
                f"Protected scientific/provenance table counts changed: "
                f"before={before}, after={after}"
            )

        fk = db.execute("PRAGMA foreign_key_check").fetchall()
        if fk:
            raise RuntimeError(f"foreign key errors: {fk[:10]}")

        db.commit()

        return {
            "version": VERSION,
            "database": str(database.resolve()),
            "authority": validation,
            "target_browser_occurrences": occurrence_rows,
            "target_browser_groups": group_rows,
            "protected_table_counts_unchanged": True,
            "foreign_key_errors": 0,
            "2O4K_DR7": {
                "ligand_instance_id": 36170,
                "canonical_target_id": "protease",
                "protein_type": "protease",
                "source_protein_type": str(rr["source_protein_type"]),
            },
            "production_scientific_tables_modified": False,
        }
    except Exception:
        db.rollback()
        raise
    finally:
        db.close()


def main():
    global VERSION, TARGET_AUTHORITY_VERSION, EXPECTED_GROUPS
    ap = argparse.ArgumentParser(
        description="Build occurrence-resolved canonical Target Browser views."
    )
    ap.add_argument("--database", type=Path, required=True)
    ap.add_argument("--authority-version", default=TARGET_AUTHORITY_VERSION)
    ap.add_argument("--expected-groups", type=int, default=EXPECTED_GROUPS)
    ap.add_argument("--view-version", default=VERSION)
    args = ap.parse_args()

    TARGET_AUTHORITY_VERSION = args.authority_version
    EXPECTED_GROUPS = args.expected_groups
    VERSION = args.view_version

    if not args.database.exists():
        raise SystemExit(f"database not found: {args.database}")

    result = run(args.database)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
