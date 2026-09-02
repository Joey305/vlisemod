#!/usr/bin/env python3
"""Manifest-driven public RCSB mmCIF acquisition for reviewer reproduction.

This adapts the project's proven cache/retry/atomic-download design while using
the frozen release manifest—not a virus/PDB CSV—as the sole source of truth.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import threading
import time
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Callable, Iterable
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

import gemmi


RCSB_URL = "https://files.rcsb.org/download/{entry_id}.cif"
USER_AGENT = "V-LiSEMOD-reviewer-reproducibility/1.0 (+https://github.com/Joey305/vlisemod)"
REQUIRED_COLUMNS = {"entry_id", "filename", "relative_path", "sha256"}


@dataclass(frozen=True)
class ManifestRecord:
    entry_id: str
    filename: str
    relative_path: str
    sha256: str
    file_size: int | None = None


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_manifest(path: Path) -> list[ManifestRecord]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or not REQUIRED_COLUMNS.issubset(reader.fieldnames):
            raise ValueError(f"Manifest must contain {sorted(REQUIRED_COLUMNS)}; found {reader.fieldnames}")
        rows = []
        paths = set()
        for line, row in enumerate(reader, start=2):
            entry_id = (row.get("entry_id") or "").strip().upper()
            relative = (row.get("relative_path") or "").strip().replace("\\", "/")
            filename = (row.get("filename") or "").strip()
            expected = (row.get("sha256") or "").strip().lower()
            relpath = Path(relative)
            if not entry_id or not filename or not relative or len(expected) != 64 or any(c not in "0123456789abcdef" for c in expected):
                raise ValueError(f"Invalid required manifest value at line {line}")
            if relpath.is_absolute() or ".." in relpath.parts or relpath.name != filename:
                raise ValueError(f"Unsafe or inconsistent relative_path at line {line}: {relative}")
            if relative in paths:
                raise ValueError(f"Duplicate manifest relative_path: {relative}")
            paths.add(relative)
            size = row.get("file_size") or ""
            rows.append(ManifestRecord(entry_id, filename, relative, expected, int(size) if size.isdigit() else None))
    if not rows:
        raise ValueError("Manifest contains no records")
    return rows


def group_by_entry(rows: Iterable[ManifestRecord]) -> dict[str, list[ManifestRecord]]:
    grouped: dict[str, list[ManifestRecord]] = defaultdict(list)
    for row in rows:
        grouped[row.entry_id].append(row)
    for entry_id, records in grouped.items():
        if len({r.sha256 for r in records}) != 1:
            raise ValueError(f"Entry {entry_id} has inconsistent frozen checksums")
    return dict(grouped)


def inspect_mmcif(path: Path, entry_id: str) -> tuple[str, str | None, str | None]:
    """Return status, observed hash, and parsed entry id without network access."""
    if not path.exists() or not path.is_file() or path.stat().st_size == 0:
        return "MISSING", None, None
    observed = sha256(path)
    try:
        document = gemmi.cif.read_file(str(path))
        parsed = document.sole_block().find_value("_entry.id").strip().upper()
    except Exception:
        return "PARSE_FAILED", observed, None
    if parsed != entry_id.upper():
        return "ENTRY_ID_MISMATCH", observed, parsed
    return "PARSED", observed, parsed


def classify_local(path: Path, record: ManifestRecord) -> tuple[str, str | None, str | None]:
    status, observed, parsed = inspect_mmcif(path, record.entry_id)
    if status != "PARSED":
        return status, observed, parsed
    return ("VERIFIED_FROZEN_INPUT" if observed == record.sha256 else "LOCAL_CHECKSUM_MISMATCH"), observed, parsed


def fetch_to_cache(entry_id: str, cache_dir: Path, timeout: float, retries: int,
                   opener: Callable = urlopen, sleeper: Callable[[float], None] = time.sleep) -> tuple[str, str]:
    """Fetch one entry atomically. A leftover .part is never a cache hit."""
    url = RCSB_URL.format(entry_id=entry_id)
    cache_dir.mkdir(parents=True, exist_ok=True)
    part = cache_dir / f"{entry_id}.{threading.get_ident()}.part"
    last_error = ""
    for attempt in range(1, retries + 1):
        try:
            request = Request(url, headers={"User-Agent": USER_AGENT})
            with opener(request, timeout=timeout) as response, part.open("wb") as output:
                shutil.copyfileobj(response, output)
            part.replace(cache_dir / f"{entry_id}.cif")
            return "DOWNLOADED", url
        except HTTPError as exc:
            last_error = f"HTTP {exc.code}: {exc.reason}"
            if exc.code == 404:
                return "NOT_FOUND", last_error
        except (URLError, TimeoutError, OSError) as exc:
            last_error = str(exc)
        if attempt < retries:
            sleeper(min(8.0, float(attempt)))
    return "DOWNLOAD_FAILED", last_error or "unknown download error"


def acquire_entry(entry_id: str, records: list[ManifestRecord], cache_dir: Path, timeout: float, retries: int,
                  opener: Callable = urlopen, sleeper: Callable[[float], None] = time.sleep) -> dict:
    expected = records[0].sha256
    cache = cache_dir / f"{entry_id}.cif"
    cached_status, observed, parsed = inspect_mmcif(cache, entry_id)
    cache_hit = cached_status == "PARSED"
    downloaded = False
    message = ""
    if cached_status != "PARSED":
        fetched, message = fetch_to_cache(entry_id, cache_dir, timeout, retries, opener, sleeper)
        if fetched != "DOWNLOADED":
            return {"entry_id": entry_id, "status": fetched, "expected_sha256": expected,
                    "observed_sha256": observed, "parsed_entry_id": parsed, "cache_hit": False,
                    "download_attempted": True, "downloaded": False, "message": message, "url": RCSB_URL.format(entry_id=entry_id),
                    "relative_paths": [r.relative_path for r in records]}
        downloaded = True
        cached_status, observed, parsed = inspect_mmcif(cache, entry_id)
    if cached_status == "PARSE_FAILED":
        status = "PARSE_FAILED"
    elif cached_status == "ENTRY_ID_MISMATCH":
        status = "ENTRY_ID_MISMATCH"
    elif cached_status != "PARSED":
        status = "DOWNLOAD_FAILED"
    elif observed == expected:
        status = "VERIFIED_FROZEN_INPUT"
    else:
        status = "UPSTREAM_REVISION_CHANGED"
    return {"entry_id": entry_id, "status": status, "expected_sha256": expected,
            "observed_sha256": observed, "parsed_entry_id": parsed, "cache_hit": cache_hit,
            "download_attempted": not cache_hit, "downloaded": downloaded, "message": message, "url": RCSB_URL.format(entry_id=entry_id),
            "relative_paths": [r.relative_path for r in records]}


def materialize(rows: list[ManifestRecord], source: Path, entry_results: dict[str, dict], allow_current_upstream: bool) -> list[dict]:
    file_results = []
    for record in rows:
        target = source / record.relative_path
        target_status, observed, parsed = classify_local(target, record)
        if target_status == "VERIFIED_FROZEN_INPUT":
            file_results.append({"entry_id": record.entry_id, "relative_path": record.relative_path, "status": "EXISTING_VALID_HIERARCHY_FILE"})
            continue
        entry = entry_results[record.entry_id]
        accepted = entry["status"] == "VERIFIED_FROZEN_INPUT" or (allow_current_upstream and entry["status"] == "UPSTREAM_REVISION_CHANGED")
        if not accepted:
            file_results.append({"entry_id": record.entry_id, "relative_path": record.relative_path,
                                 "status": target_status if target_status != "MISSING" else entry["status"]})
            continue
        target.parent.mkdir(parents=True, exist_ok=True)
        cache = source / ".download_cache" / f"{record.entry_id}.cif"
        try:
            if target.exists():
                target.unlink()
            try:
                target.hardlink_to(cache)
                action = "MATERIALIZED_HARDLINK"
            except OSError:
                shutil.copy2(cache, target)
                action = "MATERIALIZED_COPY"
            final_status, _, _ = classify_local(target, record)
            if final_status != "VERIFIED_FROZEN_INPUT" and not allow_current_upstream:
                raise RuntimeError(final_status)
            file_results.append({"entry_id": record.entry_id, "relative_path": record.relative_path, "status": action})
        except Exception as exc:
            file_results.append({"entry_id": record.entry_id, "relative_path": record.relative_path, "status": "MATERIALIZE_FAILED", "message": str(exc)})
    return file_results


def summarize(entry_results: Iterable[dict], file_results: Iterable[dict], rows: list[ManifestRecord]) -> dict:
    entries, files = list(entry_results), list(file_results)
    return {
        "manifest_rows": len(rows), "unique_pdb_entries": len(group_by_entry(rows)),
        "downloads_attempted": sum(1 for r in entries if r.get("download_attempted")),
        "downloaded": sum(1 for r in entries if r.get("downloaded")),
        "cache_hits": sum(1 for r in entries if r.get("cache_hit")),
        "exact_frozen_checksum_matches": sum(r["status"] == "VERIFIED_FROZEN_INPUT" for r in entries),
        "upstream_revision_mismatches": sum(r["status"] == "UPSTREAM_REVISION_CHANGED" for r in entries),
        "http_failures": sum(r["status"] == "DOWNLOAD_FAILED" for r in entries),
        "not_found": sum(r["status"] == "NOT_FOUND" for r in entries),
        "parse_failures": sum(r["status"] == "PARSE_FAILED" for r in entries),
        "entry_id_mismatches": sum(r["status"] == "ENTRY_ID_MISMATCH" for r in entries),
        "materialized_hierarchy_files": sum(r["status"].startswith("MATERIALIZED_") for r in files),
        "existing_valid_hierarchy_files": sum(r["status"] == "EXISTING_VALID_HIERARCHY_FILE" for r in files),
        "invalid_local_files": sum(r["status"] in {"LOCAL_CHECKSUM_MISMATCH", "PARSE_FAILED", "ENTRY_ID_MISMATCH"} for r in files),
        "entry_statuses": dict(Counter(r["status"] for r in entries)),
        "file_statuses": dict(Counter(r["status"] for r in files)),
    }


def write_reports(outputs: Path, prefix: str, summary: dict, entries: list[dict], files: list[dict]) -> None:
    outputs.mkdir(parents=True, exist_ok=True)
    payload = {"summary": summary, "entries": entries, "files": files}
    (outputs / f"{prefix}.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    columns = ["entry_id", "status", "expected_sha256", "observed_sha256", "parsed_entry_id", "cache_hit", "downloaded", "url", "relative_paths", "message"]
    with (outputs / f"{prefix}.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for entry in entries:
            row = {key: entry.get(key, "") for key in columns}
            row["relative_paths"] = ";".join(entry.get("relative_paths", []))
            writer.writerow(row)
    lines = [f"# {prefix.replace('_', ' ').title()}", "", "## Summary", ""]
    lines.extend(f"- {key.replace('_', ' ')}: **{value}**" for key, value in summary.items() if not isinstance(value, dict))
    changed = [r for r in entries if r["status"] != "VERIFIED_FROZEN_INPUT"]
    if changed:
        lines += ["", "## Non-frozen or failed entries", "", "| PDB | Status | Expected SHA-256 | Observed SHA-256 | Paths |", "| --- | --- | --- | --- | --- |"]
        lines.extend(f"| {r['entry_id']} | {r['status']} | {r.get('expected_sha256','')} | {r.get('observed_sha256','') or ''} | {'; '.join(r['relative_paths'])} |" for r in changed)
    (outputs / f"{prefix}.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def download_inputs(manifest: Path, source: Path, outputs: Path, workers: int = 6, timeout: float = 30.0,
                    retries: int = 3, allow_current_upstream: bool = False, entry_ids: set[str] | None = None,
                    opener: Callable = urlopen, sleeper: Callable[[float], None] = time.sleep,
                    dry_run: bool = False) -> dict:
    all_rows = read_manifest(manifest)
    rows = [r for r in all_rows if entry_ids is None or r.entry_id in entry_ids]
    if not rows:
        raise ValueError("No manifest rows selected for download")
    grouped = group_by_entry(rows)
    source = source.resolve()
    print(f"Input source: {source}")
    print(f"Manifest rows: {len(rows)}")
    print(f"Unique PDB entries: {len(grouped)}")
    if dry_run:
        print("DOWNLOAD PLAN: dry run; no network, cache, hierarchy, or report writes")
        return {"summary": {"manifest_rows": len(rows), "unique_pdb_entries": len(grouped), "planned_downloads": len(grouped)}}
    cache = source / ".download_cache"
    entries = []
    with ThreadPoolExecutor(max_workers=max(1, workers)) as executor:
        futures = {executor.submit(acquire_entry, entry_id, records, cache, timeout, retries, opener, sleeper): entry_id for entry_id, records in grouped.items()}
        for completed, future in enumerate(as_completed(futures), start=1):
            result = future.result()
            entries.append(result)
            print(f"[{completed}/{len(futures)}] {result['entry_id']} {result['status']}")
    entries.sort(key=lambda r: r["entry_id"])
    entry_results = {r["entry_id"]: r for r in entries}
    files = materialize(rows, source, entry_results, allow_current_upstream)
    summary = summarize(entries, files, rows)
    write_reports(outputs, "INPUT_DOWNLOAD_REPORT", summary, entries, files)
    return {"summary": summary, "entries": entries, "files": files}


def verify_inputs(manifest: Path, source: Path, outputs: Path, allow_current_upstream: bool = False,
                  entry_ids: set[str] | None = None) -> dict:
    all_rows = read_manifest(manifest)
    rows = [r for r in all_rows if entry_ids is None or r.entry_id in entry_ids]
    source = source.resolve()
    expected = {r.relative_path for r in rows}
    files = []
    for record in rows:
        status, observed, parsed = classify_local(source / record.relative_path, record)
        files.append({"entry_id": record.entry_id, "relative_path": record.relative_path, "status": status,
                      "expected_sha256": record.sha256, "observed_sha256": observed, "parsed_entry_id": parsed})
    observed_paths = {p.relative_to(source).as_posix() for p in source.rglob("*.cif") if ".download_cache" not in p.parts}
    unexpected = sorted(observed_paths - expected)
    status_counts = Counter(r["status"] for r in files)
    summary = {
        "manifest_rows": len(rows), "unique_pdb_entries": len(group_by_entry(rows)),
        "verified_frozen_files": status_counts["VERIFIED_FROZEN_INPUT"], "missing": status_counts["MISSING"],
        "checksum_mismatch": status_counts["LOCAL_CHECKSUM_MISMATCH"], "parse_failure": status_counts["PARSE_FAILED"],
        "entry_mismatch": status_counts["ENTRY_ID_MISMATCH"], "unexpected_files": len(unexpected),
        "verification_mode": "CURRENT_UPSTREAM_RECONSTRUCTION" if allow_current_upstream else "FROZEN_RELEASE_REPRODUCTION",
    }
    strict_bad = summary["missing"] or summary["parse_failure"] or summary["entry_mismatch"] or summary["unexpected_files"]
    if not allow_current_upstream:
        strict_bad = strict_bad or summary["checksum_mismatch"]
    summary["passed"] = not strict_bad
    entries = []
    for entry_id, records in group_by_entry(rows).items():
        relevant = [item for item in files if item["entry_id"] == entry_id]
        entries.append({
            "entry_id": entry_id,
            "status": "VERIFIED_FROZEN_INPUT" if all(item["status"] == "VERIFIED_FROZEN_INPUT" for item in relevant) else "LOCAL_VERIFICATION_ISSUES",
            "relative_paths": [record.relative_path for record in records],
        })
    write_reports(outputs, "INPUT_VERIFICATION_REPORT", summary, entries, files)
    print("INPUT VERIFICATION: " + ("PASS" if summary["passed"] else "FAIL"))
    for key in ("manifest_rows", "unique_pdb_entries", "verified_frozen_files", "missing", "checksum_mismatch", "parse_failure", "entry_mismatch", "unexpected_files"):
        print(f"{key.replace('_', ' ').title()}: {summary[key]}")
    if not summary["passed"]:
        raise RuntimeError("Frozen input verification failed" if not allow_current_upstream else "Input verification failed beyond allowed upstream checksum revisions")
    return {"summary": summary, "entries": entries, "files": files, "unexpected": unexpected}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", default=Path(__file__).resolve().parents[1] / "manifests/FROZEN_CIF_CORPUS_MANIFEST.csv")
    parser.add_argument("--source", default="PDB_FILES")
    parser.add_argument("--outputs", default=Path(__file__).resolve().parents[1] / "outputs")
    parser.add_argument("--workers", type=int, default=6)
    parser.add_argument("--timeout", type=float, default=30)
    parser.add_argument("--retries", type=int, default=3)
    parser.add_argument("--allow-current-upstream", action="store_true")
    parser.add_argument("--pdb-id", action="append", help="Bounded acquisition/verification selector; repeatable")
    parser.add_argument("--verify", action="store_true", help="Verify locally only; never downloads")
    parser.add_argument("--dry-run", action="store_true", help="Print selected manifest denominators without downloading or writing")
    args = parser.parse_args()
    selected = {x.upper() for x in args.pdb_id} if args.pdb_id else None
    if args.verify:
        verify_inputs(Path(args.manifest), Path(args.source), Path(args.outputs), args.allow_current_upstream, selected)
    else:
        download_inputs(Path(args.manifest), Path(args.source), Path(args.outputs), args.workers, args.timeout, args.retries, args.allow_current_upstream, selected, dry_run=args.dry_run)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
