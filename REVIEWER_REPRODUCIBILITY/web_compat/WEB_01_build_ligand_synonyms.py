#!/usr/bin/env python3
"""
WEB_01_build_ligand_synonyms.py

Build the V-LiSEMOD website compatibility table `Ligand_Synonyms`
from the current included ligand population in the rebuilt CIF-native database.

Authoritative external source:
    RCSB PDB Data API core chemical-component endpoint
    https://data.rcsb.org/rest/v1/core/chemcomp/{CCD_ID}

No third-party data file is required.

The script:
  1. Reads distinct CCD/component IDs from ligand_instances where curation_status='included'.
  2. Fetches current chemical-component metadata from the RCSB PDB Data API.
  3. Extracts the primary chemical name plus all available synonym fields.
  4. Caches raw JSON responses so interrupted runs can resume without re-fetching.
  5. Builds `Ligand_Synonyms` using a staging table and atomically replaces it only
     after all rows have been prepared.
  6. Builds `Ligand_Synonym_Status` for provenance/coverage auditing.
  7. Optionally merges synonym pairs from a legacy V-LiSEMOD database, but only for
     component IDs that are present in the current included release population.

This is a WEBSITE-COMPATIBILITY augmentation. It should be run against the deployment
copy (for example ./viral_data.db), not against the frozen scientific release master.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import datetime as dt
import json
import os
import random
import sqlite3
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

VERSION = "vlisemod-web-ligand-synonyms-v1.0"
RCSB_CHEMCOMP_URL = "https://data.rcsb.org/rest/v1/core/chemcomp/{comp_id}"
DEFAULT_CACHE_DIR = "outputs/ligand_synonyms_rcsb_cache"


def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


def normalize_component_id(value: object) -> str:
    return str(value or "").strip().upper()


def clean_name(value: object) -> Optional[str]:
    if value is None:
        return None
    text = " ".join(str(value).strip().split())
    if not text or text in {"?", "."}:
        return None
    return text


def split_semicolon_synonyms(value: object) -> List[str]:
    text = clean_name(value)
    if not text:
        return []
    return [part.strip() for part in text.split(";") if part.strip()]


def dedupe_names(items: Iterable[Tuple[str, str, int]]) -> List[Tuple[str, str, int]]:
    """
    Return unique (name, source, is_primary_name) rows, deduplicated case-insensitively.
    Primary-name provenance wins if the same text appears in multiple source fields.
    """
    best: Dict[str, Tuple[str, str, int]] = {}
    for name, source, is_primary in items:
        cleaned = clean_name(name)
        if not cleaned:
            continue
        key = cleaned.casefold()
        current = best.get(key)
        candidate = (cleaned, source, int(bool(is_primary)))
        if current is None or candidate[2] > current[2]:
            best[key] = candidate
    return sorted(best.values(), key=lambda x: (-x[2], x[0].casefold()))


def table_exists(conn: sqlite3.Connection, table_name: str) -> bool:
    row = conn.execute(
        """
        SELECT 1
        FROM sqlite_master
        WHERE type='table' AND lower(name)=lower(?)
        LIMIT 1
        """,
        (table_name,),
    ).fetchone()
    return row is not None


def current_component_ids(conn: sqlite3.Connection) -> List[str]:
    columns = {
        row[1]
        for row in conn.execute("PRAGMA table_info(ligand_instances)").fetchall()
    }
    if "label_comp_id" not in columns:
        raise RuntimeError(
            "ligand_instances.label_comp_id was not found; this script is pinned to "
            "the rebuilt CIF-native schema."
        )

    rows = conn.execute(
        """
        SELECT DISTINCT UPPER(TRIM(label_comp_id))
        FROM ligand_instances
        WHERE curation_status='included'
          AND label_comp_id IS NOT NULL
          AND TRIM(label_comp_id) <> ''
        ORDER BY 1
        """
    ).fetchall()
    return [normalize_component_id(row[0]) for row in rows if normalize_component_id(row[0])]


def legacy_pairs(
    legacy_db: Path,
    allowed_components: Sequence[str],
) -> Dict[str, List[str]]:
    allowed = set(allowed_components)
    out: Dict[str, List[str]] = {}
    with sqlite3.connect(str(legacy_db)) as conn:
        if not table_exists(conn, "Ligand_Synonyms"):
            raise RuntimeError(f"{legacy_db} does not contain Ligand_Synonyms")
        for ligand, synonym in conn.execute(
            """
            SELECT ligand, synonym
            FROM Ligand_Synonyms
            WHERE ligand IS NOT NULL
              AND synonym IS NOT NULL
              AND TRIM(synonym) <> ''
            """
        ):
            comp_id = normalize_component_id(ligand)
            syn = clean_name(synonym)
            if comp_id in allowed and syn:
                out.setdefault(comp_id, []).append(syn)
    return out


def cache_path(cache_dir: Path, comp_id: str) -> Path:
    safe = "".join(ch for ch in comp_id if ch.isalnum() or ch in {"_", "-"})
    return cache_dir / f"{safe}.json"


def load_cached_payload(cache_dir: Path, comp_id: str) -> Optional[dict]:
    path = cache_path(cache_dir, comp_id)
    if not path.exists():
        return None
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return None
    if isinstance(payload, dict):
        return payload
    return None


def save_cached_payload(cache_dir: Path, comp_id: str, payload: dict) -> None:
    cache_dir.mkdir(parents=True, exist_ok=True)
    path = cache_path(cache_dir, comp_id)
    tmp = path.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    os.replace(tmp, path)


def fetch_rcsb_payload(
    comp_id: str,
    *,
    timeout: int,
    retries: int,
    user_agent: str,
) -> Tuple[str, Optional[dict], Optional[str], int]:
    url = RCSB_CHEMCOMP_URL.format(comp_id=comp_id)
    last_error: Optional[str] = None
    for attempt in range(retries + 1):
        request = urllib.request.Request(
            url,
            headers={
                "Accept": "application/json",
                "User-Agent": user_agent,
            },
        )
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                raw = response.read().decode("utf-8", errors="replace")
                payload = json.loads(raw)
                if not isinstance(payload, dict):
                    return "invalid_json_shape", None, "RCSB returned non-object JSON", response.status
                return "ok", payload, None, response.status
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                return "not_found", None, "RCSB chemical component not found", exc.code
            last_error = f"HTTP {exc.code}: {exc.reason}"
            if exc.code not in {408, 429, 500, 502, 503, 504}:
                return "http_error", None, last_error, exc.code
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as exc:
            last_error = f"{type(exc).__name__}: {exc}"

        if attempt < retries:
            delay = min(12.0, (2 ** attempt) + random.random())
            time.sleep(delay)

    return "fetch_error", None, last_error or "unknown fetch error", 0


def extract_names(payload: dict) -> Tuple[Optional[str], List[Tuple[str, str, int]]]:
    rows: List[Tuple[str, str, int]] = []

    chem_comp = payload.get("chem_comp")
    if not isinstance(chem_comp, dict):
        chem_comp = {}

    primary_name = clean_name(chem_comp.get("name"))
    if primary_name:
        rows.append((primary_name, "RCSB:chem_comp.name", 1))

    # The wwPDB CCD commonly exposes a semicolon-delimited synonym field here.
    for synonym in split_semicolon_synonyms(chem_comp.get("pdbx_synonyms")):
        rows.append((synonym, "RCSB:chem_comp.pdbx_synonyms", 0))

    # RCSB Data API's dedicated synonym category.
    synonym_categories = [
        "rcsb_chem_comp_synonyms",
        "pdbx_chem_comp_synonyms",
    ]
    for category_name in synonym_categories:
        category = payload.get(category_name)
        if isinstance(category, dict):
            category = [category]
        if not isinstance(category, list):
            continue
        for item in category:
            if not isinstance(item, dict):
                continue
            synonym = clean_name(item.get("name"))
            if synonym:
                rows.append((synonym, f"RCSB:{category_name}.name", 0))

    return primary_name, dedupe_names(rows)


def fetch_or_cache_one(
    comp_id: str,
    *,
    cache_dir: Path,
    refresh: bool,
    offline: bool,
    timeout: int,
    retries: int,
) -> dict:
    if not refresh:
        cached = load_cached_payload(cache_dir, comp_id)
        if cached is not None:
            return {
                "component_id": comp_id,
                "status": "cached",
                "payload": cached,
                "error": None,
                "http_status": 200,
                "source_url": RCSB_CHEMCOMP_URL.format(comp_id=comp_id),
            }

    if offline:
        return {
            "component_id": comp_id,
            "status": "offline_cache_miss",
            "payload": None,
            "error": "No cached RCSB response available",
            "http_status": 0,
            "source_url": RCSB_CHEMCOMP_URL.format(comp_id=comp_id),
        }

    status, payload, error, http_status = fetch_rcsb_payload(
        comp_id,
        timeout=timeout,
        retries=retries,
        user_agent=f"V-LiSEMOD/{VERSION}",
    )
    if payload is not None:
        save_cached_payload(cache_dir, comp_id, payload)

    return {
        "component_id": comp_id,
        "status": status,
        "payload": payload,
        "error": error,
        "http_status": http_status,
        "source_url": RCSB_CHEMCOMP_URL.format(comp_id=comp_id),
    }


def build_rows(
    component_ids: Sequence[str],
    results: Sequence[dict],
    legacy: Optional[Dict[str, List[str]]] = None,
    merge_legacy: bool = False,
) -> Tuple[List[Tuple], List[Tuple]]:
    result_by_id = {row["component_id"]: row for row in results}
    synonym_rows: List[Tuple] = []
    status_rows: List[Tuple] = []
    timestamp = utc_now()

    for comp_id in component_ids:
        result = result_by_id.get(comp_id) or {
            "status": "missing_result",
            "payload": None,
            "error": "No result returned",
            "http_status": 0,
            "source_url": RCSB_CHEMCOMP_URL.format(comp_id=comp_id),
        }

        primary_name: Optional[str] = None
        names: List[Tuple[str, str, int]] = []
        payload = result.get("payload")
        if isinstance(payload, dict):
            primary_name, names = extract_names(payload)

        legacy_names = (legacy or {}).get(comp_id, [])
        if merge_legacy or (not names and legacy_names):
            for synonym in legacy_names:
                names.append((synonym, "legacy_vlisemod_synonym_table", 0))
            names = dedupe_names(names)

        for name, source, is_primary in names:
            synonym_rows.append(
                (comp_id, name, source, is_primary, timestamp)
            )

        status_rows.append(
            (
                comp_id,
                primary_name,
                result.get("status"),
                int(result.get("http_status") or 0),
                len(names),
                result.get("source_url"),
                result.get("error"),
                timestamp,
            )
        )

    return synonym_rows, status_rows


def install_tables(
    conn: sqlite3.Connection,
    synonym_rows: Sequence[Tuple],
    status_rows: Sequence[Tuple],
) -> None:
    conn.execute("DROP TABLE IF EXISTS Ligand_Synonyms__staging")
    conn.execute("DROP TABLE IF EXISTS Ligand_Synonym_Status__staging")

    conn.execute(
        """
        CREATE TABLE Ligand_Synonyms__staging (
            ligand TEXT NOT NULL,
            synonym TEXT NOT NULL,
            source TEXT NOT NULL,
            is_primary_name INTEGER NOT NULL DEFAULT 0 CHECK(is_primary_name IN (0,1)),
            retrieved_at_utc TEXT NOT NULL,
            UNIQUE(ligand, synonym COLLATE NOCASE)
        )
        """
    )
    conn.execute(
        """
        CREATE TABLE Ligand_Synonym_Status__staging (
            ligand TEXT PRIMARY KEY,
            primary_name TEXT,
            fetch_status TEXT NOT NULL,
            http_status INTEGER NOT NULL DEFAULT 0,
            synonym_count INTEGER NOT NULL DEFAULT 0,
            source_url TEXT,
            error_message TEXT,
            updated_at_utc TEXT NOT NULL
        )
        """
    )

    conn.executemany(
        """
        INSERT OR IGNORE INTO Ligand_Synonyms__staging
            (ligand, synonym, source, is_primary_name, retrieved_at_utc)
        VALUES (?, ?, ?, ?, ?)
        """,
        synonym_rows,
    )
    conn.executemany(
        """
        INSERT INTO Ligand_Synonym_Status__staging
            (ligand, primary_name, fetch_status, http_status, synonym_count,
             source_url, error_message, updated_at_utc)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?)
        """,
        status_rows,
    )

    # Every installed synonym ligand must exist in the current included release population.
    invalid = conn.execute(
        """
        SELECT COUNT(*)
        FROM Ligand_Synonyms__staging s
        WHERE NOT EXISTS (
            SELECT 1
            FROM ligand_instances i
            WHERE i.curation_status='included'
              AND UPPER(TRIM(i.label_comp_id)) = s.ligand
        )
        """
    ).fetchone()[0]
    if invalid:
        raise RuntimeError(
            f"Refusing to install synonym table: {invalid} synonym rows refer to "
            "components outside the included release population."
        )

    # Replace only after staging tables have been fully populated and validated.
    conn.execute("DROP TABLE IF EXISTS Ligand_Synonyms")
    conn.execute("ALTER TABLE Ligand_Synonyms__staging RENAME TO Ligand_Synonyms")
    conn.execute("DROP TABLE IF EXISTS Ligand_Synonym_Status")
    conn.execute("ALTER TABLE Ligand_Synonym_Status__staging RENAME TO Ligand_Synonym_Status")

    conn.execute(
        "CREATE INDEX IF NOT EXISTS idx_ligand_synonyms_ligand ON Ligand_Synonyms(ligand)"
    )
    conn.execute(
        "CREATE INDEX IF NOT EXISTS idx_ligand_synonyms_synonym_nocase "
        "ON Ligand_Synonyms(synonym COLLATE NOCASE)"
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build Ligand_Synonyms for the V-LiSEMOD website from RCSB CCD metadata."
    )
    parser.add_argument(
        "--database",
        default="./viral_data.db",
        help="Deployment/website SQLite database (default: ./viral_data.db)",
    )
    parser.add_argument(
        "--cache-dir",
        default=DEFAULT_CACHE_DIR,
        help=f"Raw RCSB JSON cache directory (default: {DEFAULT_CACHE_DIR})",
    )
    parser.add_argument(
        "--legacy-database",
        default=None,
        help="Optional old V-LiSEMOD DB used only as fallback/seed provenance.",
    )
    parser.add_argument(
        "--merge-legacy",
        action="store_true",
        help="Merge legacy synonym pairs even when current RCSB names are available.",
    )
    parser.add_argument(
        "--refresh",
        action="store_true",
        help="Ignore existing cached RCSB JSON and fetch again.",
    )
    parser.add_argument(
        "--offline",
        action="store_true",
        help="Do not access the network; use cache and optional legacy DB only.",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=4,
        help="Concurrent RCSB requests (default: 4; max: 12)",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=30,
        help="Per-request timeout seconds (default: 30)",
    )
    parser.add_argument(
        "--retries",
        type=int,
        default=3,
        help="Retries for transient RCSB failures (default: 3)",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    db_path = Path(args.database).expanduser().resolve()
    cache_dir = Path(args.cache_dir).expanduser().resolve()
    legacy_db_path = (
        Path(args.legacy_database).expanduser().resolve()
        if args.legacy_database
        else None
    )

    if not db_path.exists():
        print(f"ERROR: database not found: {db_path}", file=sys.stderr)
        return 2

    workers = min(max(int(args.workers), 1), 12)

    print(f"V-LiSEMOD ligand synonym builder: {VERSION}")
    print(f"database: {db_path}")
    print(f"cache: {cache_dir}")
    print(f"offline: {bool(args.offline)}")
    print(f"workers: {workers}")

    with sqlite3.connect(str(db_path)) as conn:
        component_ids = current_component_ids(conn)

    print(f"included unique component IDs: {len(component_ids)}")
    if not component_ids:
        print("ERROR: no included ligand component IDs found.", file=sys.stderr)
        return 2

    legacy: Dict[str, List[str]] = {}
    if legacy_db_path:
        if not legacy_db_path.exists():
            print(f"ERROR: legacy database not found: {legacy_db_path}", file=sys.stderr)
            return 2
        legacy = legacy_pairs(legacy_db_path, component_ids)
        legacy_pair_count = sum(len(v) for v in legacy.values())
        print(
            f"legacy fallback loaded: {len(legacy)} current component IDs / "
            f"{legacy_pair_count} synonym pairs"
        )

    results: List[dict] = []
    completed = 0

    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        future_map = {
            pool.submit(
                fetch_or_cache_one,
                comp_id,
                cache_dir=cache_dir,
                refresh=bool(args.refresh),
                offline=bool(args.offline),
                timeout=max(1, int(args.timeout)),
                retries=max(0, int(args.retries)),
            ): comp_id
            for comp_id in component_ids
        }

        for future in concurrent.futures.as_completed(future_map):
            comp_id = future_map[future]
            try:
                result = future.result()
            except Exception as exc:
                result = {
                    "component_id": comp_id,
                    "status": "worker_exception",
                    "payload": None,
                    "error": f"{type(exc).__name__}: {exc}",
                    "http_status": 0,
                    "source_url": RCSB_CHEMCOMP_URL.format(comp_id=comp_id),
                }
            results.append(result)
            completed += 1
            if completed == 1 or completed % 100 == 0 or completed == len(component_ids):
                print(f"RCSB metadata: {completed}/{len(component_ids)}")

    synonym_rows, status_rows = build_rows(
        component_ids,
        results,
        legacy=legacy,
        merge_legacy=bool(args.merge_legacy),
    )

    with sqlite3.connect(str(db_path)) as conn:
        conn.execute("PRAGMA foreign_keys=ON")
        conn.execute("BEGIN IMMEDIATE")
        try:
            install_tables(conn, synonym_rows, status_rows)
            conn.commit()
        except Exception:
            conn.rollback()
            raise

        total_pairs = conn.execute(
            "SELECT COUNT(*) FROM Ligand_Synonyms"
        ).fetchone()[0]
        components_with_names = conn.execute(
            "SELECT COUNT(DISTINCT ligand) FROM Ligand_Synonyms"
        ).fetchone()[0]
        primary_names = conn.execute(
            "SELECT COUNT(*) FROM Ligand_Synonyms WHERE is_primary_name=1"
        ).fetchone()[0]
        status_counts = conn.execute(
            """
            SELECT fetch_status, COUNT(*)
            FROM Ligand_Synonym_Status
            GROUP BY fetch_status
            ORDER BY fetch_status
            """
        ).fetchall()
        empty_components = conn.execute(
            """
            SELECT COUNT(*)
            FROM Ligand_Synonym_Status
            WHERE synonym_count=0
            """
        ).fetchone()[0]
        invalid = conn.execute(
            """
            SELECT COUNT(*)
            FROM Ligand_Synonyms s
            WHERE NOT EXISTS (
                SELECT 1
                FROM ligand_instances i
                WHERE i.curation_status='included'
                  AND UPPER(TRIM(i.label_comp_id))=s.ligand
            )
            """
        ).fetchone()[0]
        integrity = conn.execute("PRAGMA integrity_check").fetchone()[0]
        fk = conn.execute("PRAGMA foreign_key_check").fetchall()

    print()
    print("BUILD COMPLETE")
    print(f"unique current component IDs: {len(component_ids)}")
    print(f"components with >=1 name/synonym: {components_with_names}")
    print(f"components with zero names/synonyms: {empty_components}")
    print(f"primary chemical names: {primary_names}")
    print(f"total Ligand_Synonyms rows: {total_pairs}")
    print(f"out-of-release synonym rows: {invalid}")
    print(f"integrity_check: {integrity}")
    print(f"foreign_key_check rows: {len(fk)}")
    print("fetch status:")
    for status, count in status_counts:
        print(f"  {status}: {count}")

    if invalid != 0 or integrity != "ok" or fk:
        print("FINAL STATUS: FAIL", file=sys.stderr)
        return 1

    print("FINAL STATUS: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
