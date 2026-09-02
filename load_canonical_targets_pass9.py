#!/usr/bin/env python3
"""
Load the reviewed Pass-9 canonical target authority into V-LiSEMOD.

Default mode is VALIDATE ONLY. No database writes occur unless --apply is given.

Authoritative occurrence key:
    ligand_instance_id

The CSV pdb_id is retained only as an audit cross-check and is NOT part of the
database primary key because ligand_instance_id is already globally unique.

Expected frozen Pass-9 partition:
    YES     5414
    NO       796
    REVIEW  1145
    TOTAL   7355
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sqlite3
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

VERSION = "canonical-target-authority-pass9-v1"

EXPECTED_TOTAL = 7355
EXPECTED_ELIGIBILITY = {"YES": 5414, "NO": 796, "REVIEW": 1145}
EXPECTED_GROUPS = 32

REQUIRED_CSV_COLUMNS = {
    "pdb_id",
    "ligand_instance_id",
    "virus_name",
    "final_canonical_target_id",
    "final_canonical_target_name",
    "final_target_family",
    "final_entity_role",
    "final_target_browser_eligible",
    "final_decision",
    "final_authority_basis",
    "final_authority_note",
}


def clean(v) -> str:
    return "" if v is None else str(v).strip()


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def connect(database: Path, readonly: bool):
    if readonly:
        uri = database.resolve().as_uri() + "?mode=ro"
        db = sqlite3.connect(uri, uri=True, timeout=60)
        db.execute("PRAGMA query_only=ON")
    else:
        db = sqlite3.connect(str(database), timeout=60)
    db.row_factory = sqlite3.Row
    db.execute("PRAGMA foreign_keys=ON")
    db.execute("PRAGMA busy_timeout=60000")
    return db


def read_csv(path: Path):
    with path.open(newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        fields = set(reader.fieldnames or [])
        missing = sorted(REQUIRED_CSV_COLUMNS - fields)
        if missing:
            raise RuntimeError(f"CSV missing required columns: {missing}")
        rows = [dict(r) for r in reader]
    return rows


def db_included_inventory(db):
    rows = db.execute("""
        SELECT
            i.ligand_instance_id,
            s.entry_id AS pdb_id,
            i.label_comp_id AS ligand
        FROM ligand_instances i
        JOIN structures s ON s.structure_id=i.structure_id
        WHERE i.curation_status='included'
    """).fetchall()
    return {
        int(r["ligand_instance_id"]): {
            "pdb_id": clean(r["pdb_id"]).upper(),
            "ligand": clean(r["ligand"]).upper(),
        }
        for r in rows
    }


def validate_csv_against_db(rows, db):
    problems = []

    if len(rows) != EXPECTED_TOTAL:
        problems.append(f"CSV row count {len(rows)} != expected {EXPECTED_TOTAL}")

    ids = []
    parsed = {}

    for n, r in enumerate(rows, start=2):
        try:
            lid = int(clean(r["ligand_instance_id"]))
        except Exception:
            problems.append(
                f"CSV line {n}: invalid ligand_instance_id={r.get('ligand_instance_id')!r}"
            )
            continue

        if lid in parsed:
            problems.append(f"duplicate ligand_instance_id in CSV: {lid}")

        elig = clean(r["final_target_browser_eligible"]).upper()
        if elig not in {"YES", "NO", "REVIEW"}:
            problems.append(f"{lid}: invalid eligibility {elig!r}")

        cid = clean(r["final_canonical_target_id"])
        cname = clean(r["final_canonical_target_name"])

        if elig == "YES":
            if not cid or not cname:
                problems.append(f"{lid}: YES row lacks canonical ID/name")
            if "," in cid or ";" in cid:
                problems.append(f"{lid}: canonical ID contains comma/semicolon: {cid!r}")

        parsed[lid] = r
        ids.append(lid)

    eligibility_counts = Counter(
        clean(r["final_target_browser_eligible"]).upper() for r in rows
    )
    eligibility_counts = dict(sorted(eligibility_counts.items()))

    if eligibility_counts != EXPECTED_ELIGIBILITY:
        problems.append(
            f"eligibility counts {eligibility_counts} != expected {EXPECTED_ELIGIBILITY}"
        )

    inventory = db_included_inventory(db)
    db_ids = set(inventory)
    csv_ids = set(parsed)

    missing_from_csv = sorted(db_ids - csv_ids)
    extra_in_csv = sorted(csv_ids - db_ids)

    if missing_from_csv:
        problems.append(
            f"{len(missing_from_csv)} included DB ligand_instance_id values missing from CSV; "
            f"examples={missing_from_csv[:20]}"
        )
    if extra_in_csv:
        problems.append(
            f"{len(extra_in_csv)} CSV ligand_instance_id values are not included DB occurrences; "
            f"examples={extra_in_csv[:20]}"
        )

    pdb_mismatches = []
    for lid in sorted(db_ids & csv_ids):
        csv_pdb = clean(parsed[lid]["pdb_id"]).upper()
        db_pdb = inventory[lid]["pdb_id"]
        if csv_pdb != db_pdb:
            pdb_mismatches.append((lid, csv_pdb, db_pdb))

    if pdb_mismatches:
        problems.append(
            f"{len(pdb_mismatches)} CSV/DB PDB mismatches; examples={pdb_mismatches[:10]}"
        )

    # Exact 2O4K/DR7 regression from the real database.
    r204k = db.execute("""
        SELECT i.ligand_instance_id
        FROM ligand_instances i
        JOIN structures s ON s.structure_id=i.structure_id
        WHERE upper(s.entry_id)='2O4K'
          AND upper(i.label_comp_id)='DR7'
          AND i.curation_status='included'
    """).fetchall()

    if len(r204k) != 1:
        problems.append(f"expected one included 2O4K/DR7 occurrence; observed {len(r204k)}")
        r204k_id = None
    else:
        r204k_id = int(r204k[0]["ligand_instance_id"])
        rr = parsed.get(r204k_id)
        if rr is None:
            problems.append(f"2O4K/DR7 ligand_instance_id {r204k_id} missing from CSV")
        else:
            if clean(rr["final_target_browser_eligible"]).upper() != "YES":
                problems.append("2O4K/DR7 is not YES in Pass-9 CSV")
            if clean(rr["final_canonical_target_id"]).casefold() != "protease":
                problems.append(
                    f"2O4K/DR7 canonical ID is "
                    f"{rr['final_canonical_target_id']!r}, expected 'protease'"
                )
            if clean(rr["final_canonical_target_name"]).casefold() != "protease":
                problems.append(
                    f"2O4K/DR7 canonical name is "
                    f"{rr['final_canonical_target_name']!r}, expected 'protease'"
                )

    yes_groups = {
        (
            clean(r["virus_name"]),
            clean(r["final_canonical_target_id"]),
            clean(r["final_canonical_target_name"]),
            clean(r["final_target_family"]),
            clean(r["final_entity_role"]),
        )
        for r in rows
        if clean(r["final_target_browser_eligible"]).upper() == "YES"
    }

    if len(yes_groups) != EXPECTED_GROUPS:
        problems.append(
            f"authoritative YES canonical group count {len(yes_groups)} "
            f"!= expected {EXPECTED_GROUPS}"
        )

    return {
        "csv_rows": len(rows),
        "database_included_occurrences": len(inventory),
        "eligibility_counts": eligibility_counts,
        "canonical_yes_group_count": len(yes_groups),
        "2O4K_DR7_ligand_instance_id": r204k_id,
        "problems": problems,
    }


SCHEMA_SQL = """
CREATE TABLE canonical_ligand_targets__new (
    ligand_instance_id INTEGER PRIMARY KEY
        REFERENCES ligand_instances(ligand_instance_id),

    source_pdb_id TEXT NOT NULL,
    virus_name TEXT NOT NULL,

    canonical_target_id TEXT,
    canonical_target_name TEXT,
    target_family TEXT,
    entity_role TEXT NOT NULL,

    target_browser_eligible TEXT NOT NULL
        CHECK(target_browser_eligible IN ('YES','NO','REVIEW')),

    decision TEXT NOT NULL,
    authority_basis TEXT NOT NULL,
    authority_note TEXT NOT NULL DEFAULT '',

    authority_version TEXT NOT NULL,
    source_csv_sha256 TEXT NOT NULL,
    loaded_at TEXT NOT NULL,

    CHECK(
        target_browser_eligible <> 'YES'
        OR (
            canonical_target_id IS NOT NULL
            AND trim(canonical_target_id) <> ''
            AND canonical_target_name IS NOT NULL
            AND trim(canonical_target_name) <> ''
        )
    )
)
"""


def apply_rows(rows, database: Path, csv_sha: str):
    db = connect(database, readonly=False)
    try:
        before = {
            "structures": db.execute("SELECT count(*) FROM structures").fetchone()[0],
            "included": db.execute(
                "SELECT count(*) FROM ligand_instances WHERE curation_status='included'"
            ).fetchone()[0],
            "assessments": db.execute(
                "SELECT count(*) FROM protacability_assessment"
            ).fetchone()[0],
            "structure_classifications": db.execute(
                "SELECT count(*) FROM structure_classifications"
            ).fetchone()[0],
        }

        validation = validate_csv_against_db(rows, db)
        if validation["problems"]:
            raise RuntimeError(
                "Refusing to write because validation failed:\n- "
                + "\n- ".join(validation["problems"])
            )

        loaded_at = datetime.now(timezone.utc).isoformat()

        db.execute("BEGIN IMMEDIATE")
        db.execute("DROP TABLE IF EXISTS canonical_ligand_targets__new")
        db.execute(SCHEMA_SQL)

        payload = []
        for r in rows:
            elig = clean(r["final_target_browser_eligible"]).upper()
            cid = clean(r["final_canonical_target_id"]) or None
            cname = clean(r["final_canonical_target_name"]) or None
            family = clean(r["final_target_family"]) or None

            payload.append((
                int(clean(r["ligand_instance_id"])),
                clean(r["pdb_id"]).upper(),
                clean(r["virus_name"]),
                cid,
                cname,
                family,
                clean(r["final_entity_role"]) or "UNRESOLVED",
                elig,
                clean(r["final_decision"]),
                clean(r["final_authority_basis"]),
                clean(r["final_authority_note"]),
                VERSION,
                csv_sha,
                loaded_at,
            ))

        db.executemany("""
            INSERT INTO canonical_ligand_targets__new(
                ligand_instance_id,
                source_pdb_id,
                virus_name,
                canonical_target_id,
                canonical_target_name,
                target_family,
                entity_role,
                target_browser_eligible,
                decision,
                authority_basis,
                authority_note,
                authority_version,
                source_csv_sha256,
                loaded_at
            ) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)
        """, payload)

        observed = dict(db.execute("""
            SELECT target_browser_eligible,count(*)
            FROM canonical_ligand_targets__new
            GROUP BY target_browser_eligible
        """).fetchall())

        if observed != EXPECTED_ELIGIBILITY:
            raise RuntimeError(
                f"staging eligibility counts {observed} != {EXPECTED_ELIGIBILITY}"
            )

        missing = db.execute("""
            SELECT count(*)
            FROM ligand_instances i
            LEFT JOIN canonical_ligand_targets__new c
              ON c.ligand_instance_id=i.ligand_instance_id
            WHERE i.curation_status='included'
              AND c.ligand_instance_id IS NULL
        """).fetchone()[0]
        if missing != 0:
            raise RuntimeError(f"{missing} included ligand occurrences missing authority rows")

        extra = db.execute("""
            SELECT count(*)
            FROM canonical_ligand_targets__new c
            JOIN ligand_instances i
              ON i.ligand_instance_id=c.ligand_instance_id
            WHERE i.curation_status<>'included'
        """).fetchone()[0]
        if extra != 0:
            raise RuntimeError(f"{extra} authority rows refer to non-included occurrences")

        fk = db.execute("PRAGMA foreign_key_check").fetchall()
        if fk:
            raise RuntimeError(f"foreign key check failed before swap: {fk[:10]}")

        # These are derived authority views and are rebuilt by
        # 14b_build_canonical_target_views.py after the atomic table swap.  They
        # must be dropped first: SQLite otherwise rejects the table rename
        # because their stored SQL still resolves the old table name.
        for view_name in ("v2_target_browser_groups", "v2_target_browser_ligand_context"):
            view_type = db.execute(
                "SELECT type FROM sqlite_master WHERE name=?", (view_name,)
            ).fetchone()
            if view_type and view_type[0] == "view":
                db.execute(f'DROP VIEW "{view_name}"')

        db.execute("DROP TABLE IF EXISTS canonical_ligand_targets")
        db.execute(
            "ALTER TABLE canonical_ligand_targets__new "
            "RENAME TO canonical_ligand_targets"
        )

        db.execute("""
            CREATE INDEX IF NOT EXISTS idx_canonical_ligand_targets_browser
            ON canonical_ligand_targets(
                target_browser_eligible,
                virus_name,
                canonical_target_id
            )
        """)
        db.execute("""
            CREATE INDEX IF NOT EXISTS idx_canonical_ligand_targets_pdb
            ON canonical_ligand_targets(source_pdb_id)
        """)

        after = {
            "structures": db.execute("SELECT count(*) FROM structures").fetchone()[0],
            "included": db.execute(
                "SELECT count(*) FROM ligand_instances WHERE curation_status='included'"
            ).fetchone()[0],
            "assessments": db.execute(
                "SELECT count(*) FROM protacability_assessment"
            ).fetchone()[0],
            "structure_classifications": db.execute(
                "SELECT count(*) FROM structure_classifications"
            ).fetchone()[0],
        }

        if before != after:
            raise RuntimeError(
                f"protected table counts changed unexpectedly: before={before} after={after}"
            )

        db.commit()

        final = {
            "table_rows": db.execute(
                "SELECT count(*) FROM canonical_ligand_targets"
            ).fetchone()[0],
            "eligibility_counts": dict(db.execute("""
                SELECT target_browser_eligible,count(*)
                FROM canonical_ligand_targets
                GROUP BY target_browser_eligible
            """).fetchall()),
            "canonical_yes_groups": db.execute("""
                SELECT count(*) FROM (
                    SELECT virus_name,canonical_target_id,canonical_target_name,
                           target_family,entity_role
                    FROM canonical_ligand_targets
                    WHERE target_browser_eligible='YES'
                    GROUP BY virus_name,canonical_target_id,canonical_target_name,
                             target_family,entity_role
                )
            """).fetchone()[0],
            "protected_counts_unchanged": before == after,
            "foreign_key_errors": len(db.execute("PRAGMA foreign_key_check").fetchall()),
        }
        return final
    except Exception:
        db.rollback()
        raise
    finally:
        db.close()


def main():
    global VERSION, EXPECTED_GROUPS
    ap = argparse.ArgumentParser()
    ap.add_argument("--database", type=Path, required=True)
    ap.add_argument("--authority-csv", type=Path, required=True)
    ap.add_argument(
        "--authority-version",
        default=VERSION,
        help="Authority version to stamp after validating the supplied artifact.",
    )
    ap.add_argument(
        "--expected-groups",
        type=int,
        default=EXPECTED_GROUPS,
        help="Expected YES (virus, canonical target) group count for this artifact.",
    )
    ap.add_argument(
        "--apply",
        action="store_true",
        help="Actually create/replace canonical_ligand_targets. "
             "Without this flag, validation is read-only.",
    )
    args = ap.parse_args()

    VERSION = args.authority_version
    EXPECTED_GROUPS = args.expected_groups

    if not args.database.exists():
        raise SystemExit(f"database not found: {args.database}")
    if not args.authority_csv.exists():
        raise SystemExit(f"authority CSV not found: {args.authority_csv}")

    rows = read_csv(args.authority_csv)
    csv_sha = sha256_file(args.authority_csv)

    db = connect(args.database, readonly=True)
    try:
        validation = validate_csv_against_db(rows, db)
    finally:
        db.close()

    report = {
        "mode": "APPLY" if args.apply else "VALIDATE_ONLY",
        "database": str(args.database.resolve()),
        "authority_csv": str(args.authority_csv.resolve()),
        "authority_csv_sha256": csv_sha,
        "authority_version": VERSION,
        "validation": validation,
        "production_data_modified": False,
    }

    if validation["problems"]:
        print(json.dumps(report, indent=2, sort_keys=True))
        raise SystemExit(2)

    if not args.apply:
        print(json.dumps(report, indent=2, sort_keys=True))
        print("\nVALIDATION PASSED. No database writes were performed.")
        print("Re-run with --apply only after reviewing this output.")
        return

    final = apply_rows(rows, args.database, csv_sha)
    report["apply_result"] = final
    report["production_data_modified"] = True
    print(json.dumps(report, indent=2, sort_keys=True))
    print("\nCanonical target authority table loaded successfully.")


if __name__ == "__main__":
    main()
