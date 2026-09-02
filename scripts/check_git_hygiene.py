#!/usr/bin/env python3
"""Fail on forbidden or oversized files currently staged for commit.

Run ``python scripts/check_git_hygiene.py --staged`` before committing.  It is
repository-local and has no external dependencies or global Git side effects.
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import PurePosixPath


DEFAULT_MAX_BYTES = 50 * 1024 * 1024


def forbidden_reason(path: str) -> str | None:
    parts = PurePosixPath(path).parts
    if ".download_cache" in parts:
        return "download cache"
    if path.endswith(".part"):
        return "partial download"
    if path.startswith("PDB_FILES/") or path.startswith("REVIEWER_REPRODUCIBILITY/PDB_FILES/"):
        return "runtime-downloaded PDB corpus"
    if "PDB_FILES" in parts and not path.startswith("REVIEWER_REPRODUCIBILITY/fixture/PDB_FILES/"):
        return "non-fixture PDB corpus"
    if path.startswith(("CIF_DATABASE_REBUILD/outputs/", "CIF_DATABASE_REBUILD/logs/", "CIF_DATABASE_REBUILD/temp/")):
        return "CIF rebuild runtime output"
    if path.startswith("REVIEWER_REPRODUCIBILITY/outputs/") and not path.endswith(".gitkeep"):
        return "reviewer runtime output"
    if path.startswith("VLiSEMOD_Reviewer_Reproducibility_") and path.endswith(".zip"):
        return "generated reviewer release ZIP"
    if path.endswith((".db-wal", ".db-shm", ".db-journal")):
        return "SQLite transient file"
    if any(marker in path.lower() for marker in ("_checkpoint", "_backup", "_pre_", "_post_", "_web_")) and path.endswith(".db"):
        return "generated database checkpoint or backup"
    return None


def staged_paths() -> list[str]:
    completed = subprocess.run(
        ["git", "diff", "--cached", "--name-only", "--diff-filter=ACMR", "-z"],
        check=True, stdout=subprocess.PIPE,
    )
    return [os.fsdecode(item) for item in completed.stdout.split(b"\0") if item]


def staged_size(path: str) -> int:
    completed = subprocess.run(["git", "cat-file", "-s", f":{path}"], check=True, text=True, stdout=subprocess.PIPE)
    return int(completed.stdout.strip())


def check(paths: list[str], max_bytes: int, allow: set[str]) -> list[str]:
    failures = []
    for path in paths:
        reason = forbidden_reason(path)
        if reason:
            failures.append(f"forbidden: {path} ({reason})")
            continue
        size = staged_size(path)
        if size > max_bytes and path not in allow:
            failures.append(f"oversize: {path} ({size} bytes; limit {max_bytes}; add an explicit --allow only after review)")
    return failures


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--staged", action="store_true", help="Check staged added/copied/modified files (default)")
    parser.add_argument("--max-mb", type=float, default=50, help="Maximum staged file size in MiB")
    parser.add_argument("--allow", action="append", default=[], help="Exact reviewed staged path allowed above the size limit")
    args = parser.parse_args()
    max_bytes = int(args.max_mb * 1024 * 1024)
    paths = staged_paths()
    failures = check(paths, max_bytes, set(args.allow))
    if failures:
        print("Git hygiene: FAIL")
        print("\n".join(f"- {failure}" for failure in failures))
        return 1
    print(f"Git hygiene: PASS ({len(paths)} staged added/copied/modified files; {args.max_mb:g} MiB limit)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
