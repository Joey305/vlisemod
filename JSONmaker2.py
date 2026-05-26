#!/usr/bin/env python3
"""
JSONmaker2.py

Builds static/data/recruiter_pdb_map.json directly from the SQLite table:

    Recruiter_Code_Crosswalk

This replaces the older manually generated recruiter_pdb_map.json that contained
Windows-specific absolute paths. The new output uses the current Linux filesystem
layout on STAN/KYLE and preserves the frontend-friendly schema:

    {
        "LR00001": {
            "ligase": "CHIP",
            "pdb_id": "9QEU",
            "ligand": "A90",
            "variant": 1,
            "pdb_file": "9QEU_A90.pdb",
            "absolute_path": "/home/.../Ligases/CHIP/PDB/9QEU_A90.pdb"
        }
    }

Default expected layouts supported:

1) STAN / web-tool layout:

    /home/jxs794/WebTools/E3Recruiter_Ligandalyzer/
        JSONmaker2.py
        Ligases/
            Ligase_Recruiter.db
            CRBN/PDB/*.pdb
            VHL/PDB/*.pdb
        static/data/

2) KYLE / VLISEMOD nested module layout:

    /home/jxs794/VLISEMOD/
        JSONmaker2.py
        Ligases/MODULE/e3-recruiter-mod/Ligases/
            Ligase_Recruiter.db
            CRBN/PDB/*.pdb
            VHL/PDB/*.pdb
        static/data/
"""

from __future__ import annotations

import argparse
import json
import sqlite3
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple


# -----------------------------------------------------------------------------
# Path discovery
# -----------------------------------------------------------------------------

SCRIPT_DIR = Path(__file__).resolve().parent


def find_project_root(start: Path) -> Path:
    """
    Walk upward from the script location and return the first directory that
    looks like the Flask project root, meaning it contains a static/ folder.

    If no static/ folder is found, use the script directory.
    """

    for candidate in [start, *start.parents]:
        if (candidate / "static").is_dir():
            return candidate

    return start


PROJECT_ROOT = find_project_root(SCRIPT_DIR)


DEFAULT_DB_CANDIDATES = [
    # KYLE / VLISEMOD nested module layout
    PROJECT_ROOT / "Ligases" / "MODULE" / "e3-recruiter-mod" / "Ligases" / "Ligase_Recruiter.db",

    # STAN / E3Recruiter_Ligandalyzer root layout
    PROJECT_ROOT / "Ligases" / "Ligase_Recruiter.db",

    # If this script is placed directly inside the Ligases directory
    SCRIPT_DIR / "Ligase_Recruiter.db",
]

DEFAULT_OUTPUT_JSON = PROJECT_ROOT / "static" / "data" / "recruiter_pdb_map.json"
DEFAULT_TABLE = "Recruiter_Code_Crosswalk"


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------


def clean_text(value: Any) -> Optional[str]:
    """
    Normalize text values from SQLite.
    """

    if value is None:
        return None

    text = str(value).strip()
    if text == "" or text.lower() in {"none", "null", "nan"}:
        return None

    return text



def to_int_or_none(value: Any) -> Optional[int]:
    """
    Convert SQLite value to int when possible.
    """

    if value is None:
        return None

    try:
        return int(value)
    except Exception:
        return None



def pick_existing_path(candidates: Iterable[Path]) -> Optional[Path]:
    """
    Return the first existing path from a candidate list.
    """

    for path in candidates:
        if path.exists():
            return path.resolve()

    return None



def discover_default_db() -> Path:
    """
    Locate Ligase_Recruiter.db from known project layouts.
    """

    db = pick_existing_path(DEFAULT_DB_CANDIDATES)
    if db is not None:
        return db

    # Conservative fallback: search below the project root for a matching DB.
    matches = list(PROJECT_ROOT.rglob("Ligase_Recruiter.db"))
    if len(matches) == 1:
        return matches[0].resolve()

    if len(matches) > 1:
        msg = "Multiple Ligase_Recruiter.db files found. Please specify --db:\n"
        msg += "\n".join(f"  {m}" for m in matches)
        raise FileExistsError(msg)

    checked = "\n".join(f"  {p}" for p in DEFAULT_DB_CANDIDATES)
    raise FileNotFoundError(
        "Could not find Ligase_Recruiter.db. Checked:\n"
        f"{checked}\n\n"
        "Pass it explicitly with:\n"
        "  python JSONmaker2.py --db /path/to/Ligase_Recruiter.db"
    )



def resolve_case_insensitive(path: Path) -> Path:
    """
    Resolve a path even when one or more components differ in case.

    Useful for cases like DB ligase name PARKIN while the folder is Parkin.
    Returns the corrected path if it can be walked; otherwise returns input path.
    """

    if path.exists():
        return path

    parts = path.parts

    if not parts:
        return path

    if path.is_absolute():
        current = Path(parts[0])
        remaining = parts[1:]
    else:
        current = Path(parts[0])
        remaining = parts[1:]

    for part in remaining:
        candidate = current / part
        if candidate.exists():
            current = candidate
            continue

        if not current.is_dir():
            return path

        try:
            children = list(current.iterdir())
        except Exception:
            return path

        match = next((child for child in children if child.name.lower() == part.lower()), None)
        if match is None:
            return path

        current = match

    return current



def get_repo_root_for_db(db_path: Path) -> Path:
    """
    The DB usually lives inside a Ligases/ directory.

    PDB_File values in Recruiter_Code_Crosswalk look like:

        Ligases/CHIP/PDB/9QEU_A90.pdb

    Therefore, those paths are relative to the parent of the Ligases/ directory:

        STAN: /home/jxs794/WebTools/E3Recruiter_Ligandalyzer
        KYLE: /home/jxs794/VLISEMOD/Ligases/MODULE/e3-recruiter-mod
    """

    db_parent = db_path.resolve().parent

    if db_parent.name == "Ligases":
        return db_parent.parent

    return PROJECT_ROOT



def get_ligases_dir_for_db(db_path: Path) -> Path:
    """
    Return the directory containing ligase folders such as CRBN/PDB, VHL/PDB, etc.
    """

    db_parent = db_path.resolve().parent

    if db_parent.name == "Ligases":
        return db_parent

    # Fallbacks for unusual layouts
    candidates = [
        PROJECT_ROOT / "Ligases" / "MODULE" / "e3-recruiter-mod" / "Ligases",
        PROJECT_ROOT / "Ligases",
        SCRIPT_DIR,
    ]

    existing = pick_existing_path(candidates)
    if existing is not None:
        return existing

    return db_parent



def normalize_relative_path(path_value: Optional[str]) -> Optional[str]:
    """
    Normalize DB PDB_File paths to slash-separated relative paths for JSON.
    """

    if not path_value:
        return None

    return path_value.replace("\\", "/").lstrip("/")



def resolve_pdb_paths(
    *,
    db_path: Path,
    ligase: Optional[str],
    pdb_id: Optional[str],
    ligand: Optional[str],
    variant: Optional[int],
    pdb_file_from_db: Optional[str],
) -> Tuple[Optional[str], Optional[str], Optional[str], bool]:
    """
    Resolve PDB filename, relative path, expected absolute path, and existence.

    Returns:
        pdb_file
        relative_path
        expected_absolute_path
        exists
    """

    repo_root = get_repo_root_for_db(db_path)
    ligases_dir = get_ligases_dir_for_db(db_path)

    rel_path = normalize_relative_path(pdb_file_from_db)

    candidate_paths: List[Path] = []

    if rel_path:
        candidate_paths.append(repo_root / rel_path)

    # Fallback filename guesses if PDB_File is missing or stale.
    guessed_files: List[str] = []
    if pdb_id and ligand:
        guessed_files.append(f"{pdb_id}_{ligand}.pdb")
        if variant is not None:
            guessed_files.append(f"{pdb_id}_{ligand}_{variant}.pdb")

    # Also try the file basename from the DB path.
    if rel_path:
        guessed_files.insert(0, Path(rel_path).name)

    # Remove duplicate guesses while preserving order.
    seen = set()
    guessed_files = [x for x in guessed_files if not (x in seen or seen.add(x))]

    if ligase:
        for fname in guessed_files:
            candidate_paths.append(ligases_dir / ligase / "PDB" / fname)

        # Last-resort glob for PDB files with variants, e.g. 8OG5_LQH_1.pdb.
        pdb_folder = resolve_case_insensitive(ligases_dir / ligase / "PDB")
        if pdb_folder.is_dir() and pdb_id and ligand:
            matches = sorted(pdb_folder.glob(f"{pdb_id}_{ligand}*.pdb"))
            candidate_paths.extend(matches)

    corrected_candidates = [resolve_case_insensitive(p) for p in candidate_paths]
    existing = pick_existing_path(corrected_candidates)

    if existing is not None:
        pdb_file = existing.name

        # Prefer DB-provided relative path if it resolves; otherwise make one
        # relative to repo_root when possible.
        try:
            resolved_rel = existing.relative_to(repo_root).as_posix()
        except ValueError:
            try:
                resolved_rel = existing.relative_to(PROJECT_ROOT).as_posix()
            except ValueError:
                resolved_rel = existing.name

        return pdb_file, resolved_rel, str(existing), True

    # If no physical file exists, still preserve the expected metadata.
    if rel_path:
        expected = resolve_case_insensitive(repo_root / rel_path)
        return Path(rel_path).name, rel_path, str(expected), False

    if guessed_files and ligase:
        expected = resolve_case_insensitive(ligases_dir / ligase / "PDB" / guessed_files[0])
        try:
            expected_rel = expected.relative_to(repo_root).as_posix()
        except ValueError:
            expected_rel = str(expected)
        return guessed_files[0], expected_rel, str(expected), False

    return None, None, None, False


# -----------------------------------------------------------------------------
# SQLite reading
# -----------------------------------------------------------------------------


def connect_db(db_path: Path) -> sqlite3.Connection:
    """
    Connect to SQLite DB with row access by column name.
    """

    if not db_path.exists():
        raise FileNotFoundError(f"SQLite DB does not exist: {db_path}")

    conn = sqlite3.connect(str(db_path))
    conn.row_factory = sqlite3.Row
    return conn



def table_exists(conn: sqlite3.Connection, table: str) -> bool:
    cur = conn.execute(
        "SELECT name FROM sqlite_master WHERE type='table' AND name = ?",
        (table,),
    )
    return cur.fetchone() is not None



def get_table_columns(conn: sqlite3.Connection, table: str) -> List[str]:
    cur = conn.execute(f'PRAGMA table_info("{table}")')
    return [row[1] for row in cur.fetchall()]



def require_columns(columns: List[str], required: List[str], table: str) -> None:
    missing = [col for col in required if col not in columns]
    if missing:
        raise ValueError(
            f"Table {table!r} is missing required columns: {missing}\n"
            f"Available columns: {columns}"
        )



def fetch_crosswalk_rows(conn: sqlite3.Connection, table: str) -> List[sqlite3.Row]:
    """
    Fetch rows from Recruiter_Code_Crosswalk in a stable order.
    """

    columns = get_table_columns(conn, table)
    order_col = "RECRUITER_CODE" if "RECRUITER_CODE" in columns else columns[0]

    cur = conn.execute(f'SELECT * FROM "{table}" ORDER BY "{order_col}"')
    return cur.fetchall()


# -----------------------------------------------------------------------------
# JSON building
# -----------------------------------------------------------------------------


def build_recruiter_pdb_map(
    *,
    db_path: Path,
    table: str = DEFAULT_TABLE,
    only_existing: bool = False,
) -> Tuple[Dict[str, Dict[str, Any]], Dict[str, int], List[str]]:
    """
    Build recruiter_pdb_map.json directly from Recruiter_Code_Crosswalk.
    """

    conn = connect_db(db_path)

    try:
        if not table_exists(conn, table):
            available = [
                row[0]
                for row in conn.execute(
                    "SELECT name FROM sqlite_master WHERE type='table' ORDER BY name"
                ).fetchall()
            ]
            raise ValueError(
                f"Table {table!r} not found in DB: {db_path}\n"
                f"Available tables: {available}"
            )

        columns = get_table_columns(conn, table)
        require_columns(
            columns,
            required=["RECRUITER_CODE", "Ligase", "pdb_id", "Ligand", "Variant", "PDB_File"],
            table=table,
        )

        rows = fetch_crosswalk_rows(conn, table)

    finally:
        conn.close()

    output: Dict[str, Dict[str, Any]] = {}
    warnings: List[str] = []

    stats = {
        "rows_read": len(rows),
        "rows_written": 0,
        "missing_recruiter_code": 0,
        "existing_pdb_files": 0,
        "missing_pdb_files": 0,
        "duplicates": 0,
    }

    for row in rows:
        recruiter_code = clean_text(row["RECRUITER_CODE"])
        if not recruiter_code:
            stats["missing_recruiter_code"] += 1
            warnings.append("Skipped row with missing RECRUITER_CODE")
            continue

        ligase = clean_text(row["Ligase"])
        pdb_id = clean_text(row["pdb_id"])
        ligand = clean_text(row["Ligand"])
        variant = to_int_or_none(row["Variant"])
        pdb_file_from_db = clean_text(row["PDB_File"])

        pdb_file, relative_path, expected_abs_path, exists = resolve_pdb_paths(
            db_path=db_path,
            ligase=ligase,
            pdb_id=pdb_id,
            ligand=ligand,
            variant=variant,
            pdb_file_from_db=pdb_file_from_db,
        )

        if exists:
            stats["existing_pdb_files"] += 1
            absolute_path = expected_abs_path
        else:
            stats["missing_pdb_files"] += 1
            absolute_path = None
            if pdb_file or expected_abs_path:
                warnings.append(
                    f"Missing PDB for {recruiter_code}: "
                    f"ligase={ligase}, pdb_id={pdb_id}, ligand={ligand}, expected={expected_abs_path}"
                )

        if only_existing and not exists:
            continue

        if recruiter_code in output:
            stats["duplicates"] += 1
            warnings.append(f"Duplicate RECRUITER_CODE encountered and overwritten: {recruiter_code}")

        # Preserve old frontend-compatible fields while adding useful DB metadata.
        output[recruiter_code] = {
            "ligase": ligase,
            "pdb_id": pdb_id,
            "ligand": ligand,
            "variant": variant,
            "pdb_file": pdb_file,
            "absolute_path": absolute_path,

            # New helpful fields from the DB / resolver.
            "relative_path": relative_path,
            "expected_absolute_path": expected_abs_path,
            "exists": exists,
            "recruiter_code": recruiter_code,
            "original_ligase_ligand_code": clean_text(row["Original_Ligase_Ligand_Code"]) if "Original_Ligase_Ligand_Code" in row.keys() else None,
            "smiles": clean_text(row["SMILES"]) if "SMILES" in row.keys() else None,
            "canonical_smiles": clean_text(row["Canonical_SMILES"]) if "Canonical_SMILES" in row.keys() else None,
            "smiles_source": clean_text(row["SMILES_Source"]) if "SMILES_Source" in row.keys() else None,
            "instance_key": clean_text(row["Instance_Key"]) if "Instance_Key" in row.keys() else None,
        }

    stats["rows_written"] = len(output)

    # Stable key sort.
    output = dict(sorted(output.items(), key=lambda kv: kv[0]))

    return output, stats, warnings



def write_json(data: Dict[str, Any], output_json: Path) -> None:
    output_json.parent.mkdir(parents=True, exist_ok=True)

    with output_json.open("w", encoding="utf-8") as f:
        json.dump(data, f, indent=4, ensure_ascii=False)
        f.write("\n")



def print_summary(
    *,
    db_path: Path,
    output_json: Path,
    table: str,
    stats: Dict[str, int],
    warnings: List[str],
    dry_run: bool,
    show_missing: bool,
) -> None:
    print("\n============================================================")
    print("JSONmaker2.py — Recruiter → PDB Map Generator")
    print("============================================================")
    print(f"Project root: {PROJECT_ROOT}")
    print(f"Database:     {db_path}")
    print(f"Table:        {table}")
    print(f"Output JSON:  {output_json}")
    print("------------------------------------------------------------")
    print(f"Rows read:              {stats['rows_read']}")
    print(f"Rows written:           {stats['rows_written']}")
    print(f"Existing PDB files:     {stats['existing_pdb_files']}")
    print(f"Missing PDB files:      {stats['missing_pdb_files']}")
    print(f"Missing recruiter code: {stats['missing_recruiter_code']}")
    print(f"Duplicate codes:        {stats['duplicates']}")
    print("============================================================\n")

    if dry_run:
        print("Dry run complete. No file was written.\n")

    if show_missing and warnings:
        print("Warnings / missing files:")
        for warning in warnings:
            print(f"  - {warning}")
        print()


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate static/data/recruiter_pdb_map.json from "
            "Ligase_Recruiter.db::Recruiter_Code_Crosswalk."
        )
    )

    parser.add_argument(
        "--db",
        default=None,
        help=(
            "Path to Ligase_Recruiter.db. If omitted, the script searches known "
            "STAN/KYLE project layouts."
        ),
    )

    parser.add_argument(
        "--table",
        default=DEFAULT_TABLE,
        help=f"SQLite table to read. Default: {DEFAULT_TABLE}",
    )

    parser.add_argument(
        "--output",
        default=str(DEFAULT_OUTPUT_JSON),
        help="Output JSON path. Default: static/data/recruiter_pdb_map.json",
    )

    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Read DB and print summary, but do not write JSON.",
    )

    parser.add_argument(
        "--show-missing",
        action="store_true",
        help="Print warnings for missing/unresolved PDB files.",
    )

    parser.add_argument(
        "--only-existing",
        action="store_true",
        help="Only include rows whose PDB files exist on disk.",
    )

    return parser.parse_args()



def main() -> int:
    args = parse_args()

    try:
        db_path = Path(args.db).expanduser().resolve() if args.db else discover_default_db()
        output_json = Path(args.output).expanduser().resolve()

        data, stats, warnings = build_recruiter_pdb_map(
            db_path=db_path,
            table=args.table,
            only_existing=args.only_existing,
        )

        print_summary(
            db_path=db_path,
            output_json=output_json,
            table=args.table,
            stats=stats,
            warnings=warnings,
            dry_run=args.dry_run,
            show_missing=args.show_missing,
        )

        if not args.dry_run:
            write_json(data, output_json)
            print("✅ Successfully generated recruiter_pdb_map.json")
            print(f"   {output_json}\n")

        return 0

    except Exception as exc:
        print("\n❌ JSONmaker2.py failed.", file=sys.stderr)
        print(f"Error: {exc}\n", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
