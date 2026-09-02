#!/usr/bin/env python3
"""Read-only release gate for V-LiSEMOD canonical Target Browser authority."""
from __future__ import annotations

import argparse
import json
import sqlite3
from pathlib import Path

VERSION = "canonical-target-browser-validation-pass9-v1"
TARGET_AUTHORITY_VERSION = "canonical-target-authority-pass9-v1"

EXPECTED_TOTAL = 7355
EXPECTED_ELIGIBILITY = {"NO": 796, "REVIEW": 1145, "YES": 5414}
EXPECTED_GROUPS = 32
EXPECTED_PROTACABILITY_ASSESSMENT_ROWS = 9462


def connect_ro(database: Path):
    uri = database.resolve().as_uri() + "?mode=ro"
    db = sqlite3.connect(uri, uri=True, timeout=60)
    db.row_factory = sqlite3.Row
    db.execute("PRAGMA query_only=ON")
    db.execute("PRAGMA foreign_keys=ON")
    db.execute("PRAGMA busy_timeout=60000")
    return db


def exists(db, typ, name):
    return bool(db.execute(
        "SELECT 1 FROM sqlite_master WHERE type=? AND name=?", (typ, name)
    ).fetchone())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--database", type=Path, required=True)
    args = ap.parse_args()

    if not args.database.exists():
        raise SystemExit(f"database not found: {args.database}")

    db = connect_ro(args.database)
    checks = []

    def check(name, observed, expected):
        ok = observed == expected
        checks.append({
            "name": name,
            "status": "PASS" if ok else "FAIL",
            "observed": observed,
            "expected": expected,
        })
        return ok

    try:
        check(
            "canonical_ligand_targets exists",
            exists(db, "table", "canonical_ligand_targets"),
            True,
        )
        check(
            "v2_target_browser_ligand_context exists",
            exists(db, "view", "v2_target_browser_ligand_context"),
            True,
        )
        check(
            "v2_target_browser_groups exists",
            exists(db, "view", "v2_target_browser_groups"),
            True,
        )

        if exists(db, "table", "canonical_ligand_targets"):
            check(
                "authority row count",
                db.execute(
                    "SELECT count(*) FROM canonical_ligand_targets"
                ).fetchone()[0],
                EXPECTED_TOTAL,
            )
            elig = {
                str(r[0]): int(r[1])
                for r in db.execute(
                    """SELECT target_browser_eligible,count(*)
                       FROM canonical_ligand_targets
                       GROUP BY target_browser_eligible"""
                )
            }
            check("authority eligibility partition", elig, EXPECTED_ELIGIBILITY)

            versions = {
                str(r[0]): int(r[1])
                for r in db.execute(
                    """SELECT authority_version,count(*)
                       FROM canonical_ligand_targets
                       GROUP BY authority_version"""
                )
            }
            check(
                "authority version",
                versions,
                {TARGET_AUTHORITY_VERSION: EXPECTED_TOTAL},
            )

            bad_ids = db.execute(
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
            check("invalid authoritative canonical IDs", bad_ids, 0)

            missing = db.execute(
                """SELECT count(*)
                   FROM ligand_instances i
                   LEFT JOIN canonical_ligand_targets c
                     ON c.ligand_instance_id=i.ligand_instance_id
                   WHERE i.curation_status='included'
                     AND c.ligand_instance_id IS NULL"""
            ).fetchone()[0]
            check("included occurrences missing authority", missing, 0)

        if exists(db, "view", "v2_target_browser_ligand_context"):
            check(
                "Target Browser occurrence rows",
                db.execute(
                    "SELECT count(*) FROM v2_target_browser_ligand_context"
                ).fetchone()[0],
                EXPECTED_ELIGIBILITY["YES"],
            )

            bad_source = db.execute(
                """SELECT count(*)
                   FROM v2_target_browser_ligand_context
                   WHERE canonical_target_id IS NULL
                      OR trim(canonical_target_id)=''
                      OR canonical_target_name IS NULL
                      OR trim(canonical_target_name)=''"""
            ).fetchone()[0]
            check("Target Browser rows without canonical identity", bad_source, 0)

        if exists(db, "view", "v2_target_browser_groups"):
            check(
                "Target Browser canonical group rows",
                db.execute(
                    "SELECT count(*) FROM v2_target_browser_groups"
                ).fetchone()[0],
                EXPECTED_GROUPS,
            )

        # 2O4K regression: source contamination remains visible only as provenance,
        # while the web-facing target is protease.
        if exists(db, "view", "v2_target_browser_ligand_context"):
            reg = db.execute(
                """SELECT ligand_instance_id,virus_name,protein_type,
                          canonical_target_id,canonical_target_name,
                          source_protein_type
                   FROM v2_target_browser_ligand_context
                   WHERE upper(pdb_id)='2O4K'
                     AND upper(ligand)='DR7'"""
            ).fetchall()

            check("2O4K/DR7 Target Browser row count", len(reg), 1)
            if len(reg) == 1:
                rr = reg[0]
                check("2O4K/DR7 ligand_instance_id", int(rr["ligand_instance_id"]), 36170)
                check("2O4K/DR7 canonical ID", str(rr["canonical_target_id"]), "protease")
                check("2O4K/DR7 canonical name", str(rr["canonical_target_name"]), "protease")
                check("2O4K/DR7 web protein_type", str(rr["protein_type"]), "protease")
                check(
                    "2O4K/DR7 source provenance retained",
                    str(rr["source_protein_type"]),
                    "capsid_protein,protease",
                )

        # Scientific release-denominator guardrail.
        if exists(db, "table", "protacability_assessment"):
            check(
                "Stage-12 assessment rows unchanged",
                db.execute(
                    "SELECT count(*) FROM protacability_assessment"
                ).fetchone()[0],
                EXPECTED_PROTACABILITY_ASSESSMENT_ROWS,
            )

        check(
            "PRAGMA integrity_check",
            db.execute("PRAGMA integrity_check").fetchone()[0],
            "ok",
        )
        check(
            "PRAGMA foreign_key_check",
            len(db.execute("PRAGMA foreign_key_check").fetchall()),
            0,
        )

    finally:
        db.close()

    failed = [x for x in checks if x["status"] == "FAIL"]
    result = {
        "validator_version": VERSION,
        "database": str(args.database.resolve()),
        "passed": not failed,
        "failure_count": len(failed),
        "checks": checks,
    }
    print(json.dumps(result, indent=2, sort_keys=True))
    raise SystemExit(0 if not failed else 1)


if __name__ == "__main__":
    main()
