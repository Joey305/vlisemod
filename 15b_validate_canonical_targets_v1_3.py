#!/usr/bin/env python3
"""
Validate canonical-target authority/browser integration without modifying data.

Pass 9 / validator v1.3

This validator checks:
  * canonical authority table and browser views exist
  * authority row counts, eligibility partition, and authority version
  * canonical IDs are valid and all included occurrences have authority rows
  * Target Browser row/group counts
  * the 2O4K/DR7 regression fixture
  * current Stage-12 PROTACability v2.8 counts and scientific digests are unchanged
  * SQLite integrity and foreign-key integrity

The Stage-12 baseline in this version reflects the current LIVE database after the
intentional 2026-09-01 single-ligand reruns for ligand_instance_id 18005 and 18093.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sqlite3
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, List, Tuple


VALIDATOR_VERSION = "canonical-target-browser-validation-pass9-v1.3"
AUTHORITY_VERSION = "canonical-target-authority-pass9-v1"
PROTACABILITY_METHOD_VERSION = "protacability-cif-v2.8"

EXPECTED_AUTHORITY_ROWS = 7355
EXPECTED_ELIGIBILITY = {"NO": 796, "REVIEW": 1145, "YES": 5414}
EXPECTED_TARGET_BROWSER_OCCURRENCES = 5414
EXPECTED_TARGET_BROWSER_GROUPS = 32

EXPECTED_STAGE12_ASSESSMENT_ROWS = 9466
EXPECTED_STAGE12_ASSESSED_INSTANCES = 6786
EXPECTED_STAGE12_READINESS_ROWS = 9466
EXPECTED_STAGE12_READINESS_INSTANCES = 6786
EXPECTED_STAGE12_DIRECT_CONTACT_BAD = 0

EXPECTED_STAGE12_ASSESSMENT_SHA256 = (
    "599aa641af86b5a83dfaf29360464ef85ae2e2713f9344c282f38ae2f15bbfbd"
)
EXPECTED_STAGE12_READINESS_SHA256 = (
    "7ffe2afbf79ef2a90ae014804cd38c2be39bc1707683b9ec28e517f45240b067"
)


class ValidationError(RuntimeError):
    pass


def connect_ro(path: Path) -> sqlite3.Connection:
    if not path.exists():
        raise ValidationError(f"Database does not exist: {path}")
    uri = path.resolve().as_uri() + "?mode=ro"
    db = sqlite3.connect(uri, uri=True, timeout=60)
    db.row_factory = sqlite3.Row
    db.execute("PRAGMA query_only=ON")
    db.execute("PRAGMA foreign_keys=ON")
    return db


def object_exists(db: sqlite3.Connection, name: str) -> bool:
    row = db.execute(
        """
        SELECT 1
        FROM sqlite_master
        WHERE name=? AND type IN ('table','view')
        LIMIT 1
        """,
        (name,),
    ).fetchone()
    return bool(row)


def table_exists(db: sqlite3.Connection, name: str) -> bool:
    row = db.execute(
        "SELECT 1 FROM sqlite_master WHERE type='table' AND name=? LIMIT 1",
        (name,),
    ).fetchone()
    return bool(row)


def grouped_counts(
    db: sqlite3.Connection, sql: str, params: Iterable[Any] = ()
) -> Dict[str, int]:
    return {str(r[0]): int(r[1]) for r in db.execute(sql, tuple(params)).fetchall()}


def add_check(
    checks: List[Dict[str, Any]], name: str, expected: Any, observed: Any
) -> None:
    checks.append(
        {
            "expected": expected,
            "name": name,
            "observed": observed,
            "status": "PASS" if observed == expected else "FAIL",
        }
    )


def table_columns(db: sqlite3.Connection, table: str) -> List[str]:
    return [str(r["name"]) for r in db.execute(f'PRAGMA table_info("{table}")')]


def scientific_digest(
    db: sqlite3.Connection, table: str, method_version: str
) -> str:
    cols = table_columns(db, table)
    if "method_version" not in cols:
        raise ValidationError(f"{table} has no method_version column")

    order: List[str] = []
    for candidate in (
        "ligand_instance_id",
        "chain_id",
        "method_version",
        "assessment_id",
        "readiness_id",
    ):
        if candidate in cols:
            order.append(candidate)

    sql = f'SELECT * FROM "{table}" WHERE method_version=?'
    if order:
        sql += " ORDER BY " + ",".join(order)

    h = hashlib.sha256()
    for row in db.execute(sql, (method_version,)):
        payload = {c: row[c] for c in cols}
        h.update(
            json.dumps(
                payload,
                sort_keys=True,
                default=str,
                separators=(",", ":"),
            ).encode("utf-8")
        )
        h.update(b"\n")
    return h.hexdigest()


def scalar(db: sqlite3.Connection, sql: str, params: Iterable[Any] = ()) -> Any:
    row = db.execute(sql, tuple(params)).fetchone()
    return None if row is None else row[0]


def first_existing_column(
    db: sqlite3.Connection, object_name: str, candidates: Iterable[str]
) -> str | None:
    cols = set(table_columns(db, object_name))
    for candidate in candidates:
        if candidate in cols:
            return candidate
    return None


def quote_ident(name: str) -> str:
    return '"' + name.replace('"', '""') + '"'


def detect_eligibility_column(
    db: sqlite3.Connection, table: str
) -> tuple[str | None, Dict[str, int] | None]:
    """Find the authority eligibility column robustly.

    Prefer known names, then fall back to detecting the column whose normalized
    values partition the table exactly into YES/REVIEW/NO with the expected
    counts. This avoids coupling the validator to one CSV/schema field name.
    """
    cols = table_columns(db, table)
    preferred = (
        "eligibility",
        "eligibility_status",
        "canonical_eligibility",
        "include_status",
        "eligibility_decision",
        "target_eligibility",
        "canonical_target_eligibility",
        "canonical_target_status",
    )

    ordered = [c for c in preferred if c in cols] + [c for c in cols if c not in preferred]
    for col in ordered:
        q = quote_ident(col)
        try:
            counts = {
                str(r[0]).strip().upper(): int(r[1])
                for r in db.execute(
                    f"""
                    SELECT upper(trim(CAST({q} AS TEXT))), count(*)
                    FROM {quote_ident(table)}
                    WHERE {q} IS NOT NULL
                    GROUP BY upper(trim(CAST({q} AS TEXT)))
                    ORDER BY upper(trim(CAST({q} AS TEXT)))
                    """
                ).fetchall()
            }
        except sqlite3.Error:
            continue
        if counts == EXPECTED_ELIGIBILITY:
            return col, counts
    return None, None


def validate(db: sqlite3.Connection) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    checks: List[Dict[str, Any]] = []
    metrics: Dict[str, Any] = {}

    # ------------------------------------------------------------------
    # Required canonical objects
    # ------------------------------------------------------------------
    for name in (
        "canonical_ligand_targets",
        "v2_target_browser_ligand_context",
        "v2_target_browser_groups",
    ):
        add_check(checks, f"{name} exists", True, object_exists(db, name))

    if not table_exists(db, "canonical_ligand_targets"):
        return checks, metrics

    authority_cols = set(table_columns(db, "canonical_ligand_targets"))

    # ------------------------------------------------------------------
    # Authority table counts/version
    # ------------------------------------------------------------------
    authority_rows = int(scalar(db, "SELECT count(*) FROM canonical_ligand_targets") or 0)
    add_check(checks, "authority row count", EXPECTED_AUTHORITY_ROWS, authority_rows)

    eligibility_col, eligibility = detect_eligibility_column(
        db, "canonical_ligand_targets"
    )
    if eligibility_col is None or eligibility is None:
        add_check(
            checks,
            "authority eligibility partition",
            EXPECTED_ELIGIBILITY,
            "MISSING_OR_UNRECOGNIZED_COLUMN",
        )
    else:
        add_check(
            checks,
            "authority eligibility partition",
            EXPECTED_ELIGIBILITY,
            eligibility,
        )
        metrics["authority_eligibility_column"] = eligibility_col

    version_col = next(
        (c for c in ("authority_version", "version") if c in authority_cols),
        None,
    )
    if version_col is None:
        add_check(
            checks,
            "authority version",
            {AUTHORITY_VERSION: EXPECTED_AUTHORITY_ROWS},
            "MISSING_COLUMN",
        )
    else:
        versions = grouped_counts(
            db,
            f"""
            SELECT {version_col}, count(*)
            FROM canonical_ligand_targets
            GROUP BY {version_col}
            ORDER BY {version_col}
            """,
        )
        add_check(
            checks,
            "authority version",
            {AUTHORITY_VERSION: EXPECTED_AUTHORITY_ROWS},
            versions,
        )

    # ------------------------------------------------------------------
    # Canonical ID validity
    # ------------------------------------------------------------------
    canonical_id_col = first_existing_column(
        db,
        "canonical_ligand_targets",
        ("canonical_target_id", "canonical_id"),
    )
    if canonical_id_col is None:
        add_check(checks, "invalid authoritative canonical IDs", 0, "MISSING_COLUMN")
    elif eligibility_col is None:
        add_check(
            checks,
            "invalid authoritative canonical IDs",
            0,
            "ELIGIBILITY_COLUMN_UNAVAILABLE",
        )
    else:
        q_id = quote_ident(canonical_id_col)
        q_elig = quote_ident(eligibility_col)
        invalid_ids = int(
            scalar(
                db,
                f"""
                SELECT count(*)
                FROM canonical_ligand_targets
                WHERE upper(trim(CAST({q_elig} AS TEXT)))='YES'
                  AND (
                       {q_id} IS NULL
                    OR trim(CAST({q_id} AS TEXT))=''
                    OR CAST({q_id} AS TEXT) GLOB '*[^a-z0-9_]*'
                  )
                """,
            )
            or 0
        )
        add_check(checks, "invalid authoritative canonical IDs", 0, invalid_ids)

    # ------------------------------------------------------------------
    # Included occurrence coverage.
    # We prefer ligand_instances.is_included when present; otherwise the
    # authority-table total itself is the validated included universe.
    # ------------------------------------------------------------------
    missing_authority = 0
    if table_exists(db, "ligand_instances"):
        li_cols = set(table_columns(db, "ligand_instances"))
        if "is_included" in li_cols and "ligand_instance_id" in li_cols:
            missing_authority = int(
                scalar(
                    db,
                    """
                    SELECT count(*)
                    FROM ligand_instances li
                    LEFT JOIN canonical_ligand_targets c
                      ON c.ligand_instance_id=li.ligand_instance_id
                    WHERE li.is_included=1
                      AND c.ligand_instance_id IS NULL
                    """,
                )
                or 0
            )
        else:
            missing_authority = max(0, EXPECTED_AUTHORITY_ROWS - authority_rows)
    add_check(checks, "included occurrences missing authority", 0, missing_authority)

    # ------------------------------------------------------------------
    # Target Browser views
    # ------------------------------------------------------------------
    if object_exists(db, "v2_target_browser_ligand_context"):
        browser_occurrences = int(
            scalar(db, "SELECT count(*) FROM v2_target_browser_ligand_context") or 0
        )
        add_check(
            checks,
            "Target Browser occurrence rows",
            EXPECTED_TARGET_BROWSER_OCCURRENCES,
            browser_occurrences,
        )

        browser_cols = set(table_columns(db, "v2_target_browser_ligand_context"))
        browser_canonical_col = next(
            (c for c in ("canonical_target_id", "canonical_id") if c in browser_cols),
            None,
        )
        if browser_canonical_col is None:
            add_check(
                checks,
                "Target Browser rows without canonical identity",
                0,
                "MISSING_COLUMN",
            )
        else:
            missing_identity = int(
                scalar(
                    db,
                    f"""
                    SELECT count(*)
                    FROM v2_target_browser_ligand_context
                    WHERE {browser_canonical_col} IS NULL
                       OR trim({browser_canonical_col})=''
                    """,
                )
                or 0
            )
            add_check(
                checks,
                "Target Browser rows without canonical identity",
                0,
                missing_identity,
            )

    if object_exists(db, "v2_target_browser_groups"):
        browser_groups = int(scalar(db, "SELECT count(*) FROM v2_target_browser_groups") or 0)
        add_check(
            checks,
            "Target Browser canonical group rows",
            EXPECTED_TARGET_BROWSER_GROUPS,
            browser_groups,
        )

    # ------------------------------------------------------------------
    # 2O4K / DR7 regression fixture
    # ------------------------------------------------------------------
    if object_exists(db, "v2_target_browser_ligand_context"):
        cols = set(table_columns(db, "v2_target_browser_ligand_context"))
        pdb_col = next((c for c in ("pdb_id", "entry_id") if c in cols), None)
        ligand_col = next(
            (c for c in ("ligand", "label_comp_id", "ligand_code") if c in cols),
            None,
        )
        id_col = "ligand_instance_id" if "ligand_instance_id" in cols else None
        canonical_col = next(
            (c for c in ("canonical_target_id", "canonical_id") if c in cols),
            None,
        )
        name_col = next(
            (
                c
                for c in (
                    "canonical_target_name",
                    "canonical_name",
                    "protein_type",
                )
                if c in cols
            ),
            None,
        )
        protein_type_col = "protein_type" if "protein_type" in cols else None
        source_col = next(
            (
                c
                for c in (
                    "source_protein_type",
                    "source_protein_types",
                    "protein_type_source",
                )
                if c in cols
            ),
            None,
        )

        if pdb_col and ligand_col:
            rows = db.execute(
                f"""
                SELECT *
                FROM v2_target_browser_ligand_context
                WHERE upper({pdb_col})='2O4K'
                  AND upper({ligand_col})='DR7'
                """
            ).fetchall()
        else:
            rows = []

        add_check(checks, "2O4K/DR7 Target Browser row count", 1, len(rows))
        if len(rows) == 1:
            r = rows[0]
            if id_col:
                add_check(checks, "2O4K/DR7 ligand_instance_id", 36170, int(r[id_col]))
            if canonical_col:
                add_check(checks, "2O4K/DR7 canonical ID", "protease", r[canonical_col])
            if name_col:
                add_check(checks, "2O4K/DR7 canonical name", "protease", r[name_col])
            if protein_type_col:
                add_check(checks, "2O4K/DR7 web protein_type", "protease", r[protein_type_col])
            if source_col:
                add_check(
                    checks,
                    "2O4K/DR7 source provenance retained",
                    "capsid_protein,protease",
                    r[source_col],
                )

    # ------------------------------------------------------------------
    # Stage-12 preservation gate
    # ------------------------------------------------------------------
    required_stage12_tables = (
        "protacability_assessment",
        "protacability_degrader_readiness",
    )
    stage12_available = all(table_exists(db, t) for t in required_stage12_tables)

    if not stage12_available:
        add_check(checks, "Stage-12 tables available", True, False)
    else:
        assessment_rows = int(
            scalar(
                db,
                """
                SELECT count(*)
                FROM protacability_assessment
                WHERE method_version=?
                """,
                (PROTACABILITY_METHOD_VERSION,),
            )
            or 0
        )
        assessed_instances = int(
            scalar(
                db,
                """
                SELECT count(DISTINCT ligand_instance_id)
                FROM protacability_assessment
                WHERE method_version=?
                """,
                (PROTACABILITY_METHOD_VERSION,),
            )
            or 0
        )
        readiness_rows = int(
            scalar(
                db,
                """
                SELECT count(*)
                FROM protacability_degrader_readiness
                WHERE method_version=?
                """,
                (PROTACABILITY_METHOD_VERSION,),
            )
            or 0
        )
        readiness_instances = int(
            scalar(
                db,
                """
                SELECT count(DISTINCT ligand_instance_id)
                FROM protacability_degrader_readiness
                WHERE method_version=?
                """,
                (PROTACABILITY_METHOD_VERSION,),
            )
            or 0
        )

        assessment_cols = set(table_columns(db, "protacability_assessment"))
        direct_contact_bad = 0
        if {
            "target_chain_selection_basis",
            "ligand_target_contact_pair_count",
        }.issubset(assessment_cols):
            direct_contact_bad = int(
                scalar(
                    db,
                    """
                    SELECT count(*)
                    FROM protacability_assessment
                    WHERE method_version=?
                      AND (
                        target_chain_selection_basis
                          <> 'stage09_arpeggio_direct_protein_contact'
                        OR ligand_target_contact_pair_count <= 0
                      )
                    """,
                    (PROTACABILITY_METHOD_VERSION,),
                )
                or 0
            )

        assessment_sha = scientific_digest(
            db, "protacability_assessment", PROTACABILITY_METHOD_VERSION
        )
        readiness_sha = scientific_digest(
            db, "protacability_degrader_readiness", PROTACABILITY_METHOD_VERSION
        )

        metrics.update(
            {
                "protacability_method_version": PROTACABILITY_METHOD_VERSION,
                "stage12_assessment_rows": assessment_rows,
                "stage12_assessed_ligand_instances": assessed_instances,
                "stage12_readiness_rows": readiness_rows,
                "stage12_readiness_ligand_instances": readiness_instances,
                "stage12_assessment_sha256": assessment_sha,
                "stage12_readiness_sha256": readiness_sha,
            }
        )

        add_check(
            checks,
            f"Stage-12 assessment rows ({PROTACABILITY_METHOD_VERSION}) unchanged",
            EXPECTED_STAGE12_ASSESSMENT_ROWS,
            assessment_rows,
        )
        add_check(
            checks,
            f"Stage-12 assessed ligand instances ({PROTACABILITY_METHOD_VERSION}) unchanged",
            EXPECTED_STAGE12_ASSESSED_INSTANCES,
            assessed_instances,
        )
        add_check(
            checks,
            f"Stage-12 readiness rows ({PROTACABILITY_METHOD_VERSION}) unchanged",
            EXPECTED_STAGE12_READINESS_ROWS,
            readiness_rows,
        )
        add_check(
            checks,
            f"Stage-12 readiness ligand instances ({PROTACABILITY_METHOD_VERSION}) unchanged",
            EXPECTED_STAGE12_READINESS_INSTANCES,
            readiness_instances,
        )
        add_check(
            checks,
            f"Stage-12 direct-contact basis ({PROTACABILITY_METHOD_VERSION}) unchanged",
            EXPECTED_STAGE12_DIRECT_CONTACT_BAD,
            direct_contact_bad,
        )
        add_check(
            checks,
            f"Stage-12 assessment SHA256 ({PROTACABILITY_METHOD_VERSION}) unchanged",
            EXPECTED_STAGE12_ASSESSMENT_SHA256,
            assessment_sha,
        )
        add_check(
            checks,
            f"Stage-12 readiness SHA256 ({PROTACABILITY_METHOD_VERSION}) unchanged",
            EXPECTED_STAGE12_READINESS_SHA256,
            readiness_sha,
        )

    # ------------------------------------------------------------------
    # Database integrity
    # ------------------------------------------------------------------
    integrity = str(scalar(db, "PRAGMA integrity_check") or "")
    add_check(checks, "PRAGMA integrity_check", "ok", integrity)

    fk_errors = len(db.execute("PRAGMA foreign_key_check").fetchall())
    add_check(checks, "PRAGMA foreign_key_check", 0, fk_errors)

    return checks, metrics


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Validate Pass-9 canonical target authority/browser integration."
    )
    parser.add_argument("--database", required=True, help="SQLite database to validate")
    args = parser.parse_args()

    db_path = Path(args.database).expanduser().resolve()

    try:
        db = connect_ro(db_path)
        try:
            checks, metrics = validate(db)
        finally:
            db.close()
    except Exception as exc:
        result = {
            "database": str(db_path),
            "error": f"{type(exc).__name__}: {exc}",
            "passed": False,
            "validator_version": VALIDATOR_VERSION,
        }
        print(json.dumps(result, indent=2, sort_keys=True))
        return 2

    failures = [c for c in checks if c["status"] != "PASS"]
    result = {
        "checks": checks,
        "database": str(db_path),
        "failure_count": len(failures),
        "metrics": metrics,
        "passed": len(failures) == 0,
        "validator_version": VALIDATOR_VERSION,
    }
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if not failures else 1


if __name__ == "__main__":
    sys.exit(main())
