#!/usr/bin/env python3
"""
WEB_02_build_ligand_arp_diagram.py

Create the legacy-facing `Ligand_Arp_Diagram` compatibility VIEW used by the
V-LiSEMOD Flask application, using the already-validated Stage 14
`v2_ligand_context` compatibility view.

This is a WEBSITE INTEGRATION step, not a scientific-analysis stage.

Legacy columns exposed:
    virus_name
    pdb_id
    ligand
    chain
    ligand_id

Why v2_ligand_context is the source:
- structure_context can legitimately contain more than one classification for a PDB.
- Stage 14 already resolves/projects those classification contexts into
  v2_ligand_context.
- Re-deriving a single virus label per PDB would silently discard valid context.

Only curation_status='included' ligand instances are exposed.

The authoritative occurrence-level identity remains in ligand_instances and the
v2_* views. The five-column legacy view is only a website-compatibility surface.
"""

from __future__ import annotations

import argparse
import sqlite3
import sys
from pathlib import Path

VERSION = "vlisemod-web-ligand-arp-diagram-v1.2"


def object_type(conn: sqlite3.Connection, name: str):
    row = conn.execute(
        "SELECT type FROM sqlite_master WHERE lower(name)=lower(?) LIMIT 1",
        (name,),
    ).fetchone()
    return row[0] if row else None


def object_columns(conn: sqlite3.Connection, name: str):
    return {row[1] for row in conn.execute(f"PRAGMA table_info({name})").fetchall()}


def require_columns(conn: sqlite3.Connection, obj: str, required):
    typ = object_type(conn, obj)
    if typ not in {"table", "view"}:
        raise RuntimeError(f"Required database object is missing: {obj}")
    cols = object_columns(conn, obj)
    missing = [c for c in required if c not in cols]
    if missing:
        raise RuntimeError(
            f"{obj} is missing required columns: {', '.join(missing)}"
        )


def parse_args():
    p = argparse.ArgumentParser(
        description="Build the Ligand_Arp_Diagram website compatibility view."
    )
    p.add_argument(
        "--database",
        default="./viral_data.db",
        help="Website/deployment SQLite database (default: ./viral_data.db)",
    )
    p.add_argument(
        "--full-integrity",
        action="store_true",
        help="Run PRAGMA integrity_check; otherwise run faster PRAGMA quick_check.",
    )
    return p.parse_args()


def main():
    args = parse_args()
    db_path = Path(args.database).expanduser().resolve()

    if not db_path.exists():
        print(f"ERROR: database not found: {db_path}", file=sys.stderr)
        return 2

    print(f"V-LiSEMOD Ligand_Arp_Diagram builder: {VERSION}")
    print(f"database: {db_path}")

    with sqlite3.connect(str(db_path)) as conn:
        conn.row_factory = sqlite3.Row

        require_columns(
            conn,
            "ligand_instances",
            ["ligand_instance_id", "label_comp_id", "curation_status"],
        )
        require_columns(
            conn,
            "v2_ligand_context",
            [
                "virus_name",
                "pdb_id",
                "ligand_instance_id",
                "ligand",
                "chain",
                "ligand_residue_id",
                "curation_status",
            ],
        )

        included_instances = conn.execute(
            """
            SELECT COUNT(*)
            FROM ligand_instances
            WHERE curation_status='included'
            """
        ).fetchone()[0]

        included_components = conn.execute(
            """
            SELECT COUNT(DISTINCT label_comp_id)
            FROM ligand_instances
            WHERE curation_status='included'
            """
        ).fetchone()[0]

        represented_instances = conn.execute(
            """
            SELECT COUNT(DISTINCT ligand_instance_id)
            FROM v2_ligand_context
            WHERE curation_status='included'
            """
        ).fetchone()[0]

        missing_instance_ids = conn.execute(
            """
            SELECT COUNT(*)
            FROM (
                SELECT ligand_instance_id
                FROM ligand_instances
                WHERE curation_status='included'
                EXCEPT
                SELECT ligand_instance_id
                FROM v2_ligand_context
                WHERE curation_status='included'
            )
            """
        ).fetchone()[0]

        extra_instance_ids = conn.execute(
            """
            SELECT COUNT(*)
            FROM (
                SELECT ligand_instance_id
                FROM v2_ligand_context
                WHERE curation_status='included'
                EXCEPT
                SELECT ligand_instance_id
                FROM ligand_instances
                WHERE curation_status='included'
            )
            """
        ).fetchone()[0]

        multi_virus_instances = conn.execute(
            """
            SELECT COUNT(*)
            FROM (
                SELECT ligand_instance_id
                FROM v2_ligand_context
                WHERE curation_status='included'
                GROUP BY ligand_instance_id
                HAVING COUNT(DISTINCT COALESCE(virus_name,'')) > 1
            )
            """
        ).fetchone()[0]

        multi_context_instances = conn.execute(
            """
            SELECT COUNT(*)
            FROM (
                SELECT ligand_instance_id
                FROM v2_ligand_context
                WHERE curation_status='included'
                GROUP BY ligand_instance_id
                HAVING COUNT(*) > 1
            )
            """
        ).fetchone()[0]

        blank_required_source_rows = conn.execute(
            """
            SELECT COUNT(*)
            FROM v2_ligand_context
            WHERE curation_status='included'
              AND (
                    virus_name IS NULL OR TRIM(virus_name)=''
                 OR pdb_id IS NULL OR TRIM(pdb_id)=''
                 OR ligand IS NULL OR TRIM(ligand)=''
                 OR chain IS NULL OR TRIM(chain)=''
                 OR ligand_residue_id IS NULL OR TRIM(CAST(ligand_residue_id AS TEXT))=''
              )
            """
        ).fetchone()[0]

        # Preview the legacy projection before creating anything.
        legacy_rows = conn.execute(
            """
            SELECT COUNT(*)
            FROM (
                SELECT DISTINCT
                    virus_name,
                    pdb_id,
                    ligand,
                    chain,
                    CAST(ligand_residue_id AS TEXT) AS ligand_id
                FROM v2_ligand_context
                WHERE curation_status='included'
            )
            """
        ).fetchone()[0]

        # Count legacy keys that merge more than one authoritative occurrence.
        legacy_key_collision_groups = conn.execute(
            """
            SELECT COUNT(*)
            FROM (
                SELECT
                    virus_name,
                    pdb_id,
                    ligand,
                    chain,
                    CAST(ligand_residue_id AS TEXT) AS ligand_id
                FROM v2_ligand_context
                WHERE curation_status='included'
                GROUP BY
                    virus_name,
                    pdb_id,
                    ligand,
                    chain,
                    CAST(ligand_residue_id AS TEXT)
                HAVING COUNT(DISTINCT ligand_instance_id) > 1
            )
            """
        ).fetchone()[0]

        print()
        print("PRE-BUILD AUDIT")
        print(f"authoritative included ligand instances: {included_instances}")
        print(f"authoritative unique included components: {included_components}")
        print(f"included instances represented in v2_ligand_context: {represented_instances}")
        print(f"missing included instance IDs from v2_ligand_context: {missing_instance_ids}")
        print(f"extra included instance IDs in v2_ligand_context: {extra_instance_ids}")
        print(f"included instances with >1 v2 context row: {multi_context_instances}")
        print(f"included instances with >1 virus_name: {multi_virus_instances}")
        print(f"blank required source rows: {blank_required_source_rows}")
        print(f"legacy-facing distinct rows to create: {legacy_rows}")
        print(f"legacy five-column key collision groups: {legacy_key_collision_groups}")

        if represented_instances != included_instances:
            raise RuntimeError(
                "Refusing to build view: v2_ligand_context does not represent every "
                "authoritative included ligand instance."
            )
        if missing_instance_ids or extra_instance_ids:
            raise RuntimeError(
                "Refusing to build view: occurrence-ID mismatch between ligand_instances "
                "and v2_ligand_context."
            )
        if blank_required_source_rows:
            raise RuntimeError(
                f"Refusing to build view: {blank_required_source_rows} included context "
                "rows have blank required website fields."
            )

        existing = object_type(conn, "Ligand_Arp_Diagram")
        if existing == "table":
            raise RuntimeError(
                "Ligand_Arp_Diagram already exists as a TABLE. Refusing to drop it "
                "automatically."
            )
        if existing == "view":
            conn.execute("DROP VIEW Ligand_Arp_Diagram")

        conn.execute(
            """
            CREATE VIEW Ligand_Arp_Diagram AS
            SELECT DISTINCT
                virus_name,
                pdb_id,
                ligand,
                chain,
                CAST(ligand_residue_id AS TEXT) AS ligand_id
            FROM v2_ligand_context
            WHERE curation_status='included'
            """
        )
        conn.commit()

        view_rows = conn.execute(
            "SELECT COUNT(*) FROM Ligand_Arp_Diagram"
        ).fetchone()[0]

        view_components = conn.execute(
            "SELECT COUNT(DISTINCT ligand) FROM Ligand_Arp_Diagram"
        ).fetchone()[0]

        blank_view_rows = conn.execute(
            """
            SELECT COUNT(*)
            FROM Ligand_Arp_Diagram
            WHERE virus_name IS NULL OR TRIM(virus_name)=''
               OR pdb_id IS NULL OR TRIM(pdb_id)=''
               OR ligand IS NULL OR TRIM(ligand)=''
               OR chain IS NULL OR TRIM(chain)=''
               OR ligand_id IS NULL OR TRIM(ligand_id)=''
            """
        ).fetchone()[0]

        dr7 = conn.execute(
            """
            SELECT virus_name, pdb_id, ligand, chain, ligand_id
            FROM Ligand_Arp_Diagram
            WHERE UPPER(pdb_id)='3EKY'
              AND UPPER(ligand)='DR7'
            ORDER BY virus_name, chain, ligand_id
            """
        ).fetchall()

        check_name = "integrity_check" if args.full_integrity else "quick_check"
        print(f"running PRAGMA {check_name}...")
        integrity = conn.execute(f"PRAGMA {check_name}").fetchone()[0]
        fk_rows = conn.execute("PRAGMA foreign_key_check").fetchall()

    print()
    print("BUILD COMPLETE")
    print("object type: view")
    print(f"Ligand_Arp_Diagram rows: {view_rows}")
    print(f"unique ligands in compatibility view: {view_components}")
    print(f"blank required-field rows: {blank_view_rows}")
    print(f"3EKY/DR7 rows: {len(dr7)}")
    for row in dr7:
        print(
            "  3EKY/DR7 -> "
            f"virus={row['virus_name']} chain={row['chain']} ligand_id={row['ligand_id']}"
        )
    print(f"{check_name}: {integrity}")
    print(f"foreign_key_check rows: {len(fk_rows)}")

    fail = False
    if view_rows != legacy_rows:
        print(
            f"ERROR: created view row count {view_rows} != audited projection {legacy_rows}",
            file=sys.stderr,
        )
        fail = True
    if view_components != included_components:
        print(
            f"ERROR: view unique ligand count {view_components} != authoritative "
            f"included component count {included_components}",
            file=sys.stderr,
        )
        fail = True
    if blank_view_rows:
        print("ERROR: blank required fields exist in compatibility view", file=sys.stderr)
        fail = True
    if integrity != "ok":
        print(f"ERROR: SQLite {check_name} returned {integrity}", file=sys.stderr)
        fail = True
    if fk_rows:
        print("ERROR: foreign_key_check returned rows", file=sys.stderr)
        fail = True

    if fail:
        print("FINAL STATUS: FAIL", file=sys.stderr)
        return 1

    print("FINAL STATUS: PASS")
    if multi_virus_instances:
        print(
            "NOTE: multiple virus classifications are preserved from the validated "
            "Stage 14 context view rather than forcing one virus label per PDB."
        )
    if legacy_key_collision_groups:
        print(
            "NOTE: the old five-column website key merges some richer occurrence "
            "identities. The authoritative occurrence-level identity remains in "
            "ligand_instances / v2_ligand_context."
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
