#!/usr/bin/env python3
"""
SQLite backup utility (verbose + safe).
- Uses sqlite3 backup API to make a consistent .db snapshot
- Gzips the snapshot to backups/sqlite/<db>.<YYYYmmdd-HHMMSS>.sqlite.gz
- Keeps last N snapshots per DB and prunes older ones
"""

import os
import sqlite3
import time
import gzip
import shutil
from pathlib import Path
import glob
import sys

# === Config ===
PROJECT_ROOT = Path(__file__).resolve().parents[1]  # repo root (.. from TOOLS/)
DBS = [
    PROJECT_ROOT / "viral_data.db",
]
BACKUP_DIR = PROJECT_ROOT / "backups" / "sqlite"
KEEP_PER_DB = int(os.environ.get("SQLITE_BACKUP_KEEP", "14"))  # keep last 14 by default

def log(msg: str):
    print(msg, flush=True)

def ensure_dir(p: Path):
    p.mkdir(parents=True, exist_ok=True)

def snapshot_sqlite(src_db: Path, tmp_out: Path):
    """Create a consistent snapshot of src_db into tmp_out using sqlite backup API."""
    # open source database (readonly)
    src = sqlite3.connect(f"file:{src_db}?mode=ro", uri=True)
    # create destination/new DB file
    dst = sqlite3.connect(str(tmp_out))
    try:
        log(f"   - backing up with sqlite3 API …")
        src.backup(dst)  # consistent copy
    finally:
        dst.close()
        src.close()

def gzip_file(src: Path, dest_gz: Path):
    with open(src, "rb") as fin, gzip.open(dest_gz, "wb", compresslevel=6) as fout:
        shutil.copyfileobj(fin, fout)

def prune_old(backups_dir: Path, stem: str):
    pattern = str(backups_dir / f"{stem}.*.sqlite.gz")
    files = sorted(glob.glob(pattern))
    excess = max(0, len(files) - KEEP_PER_DB)
    for old in files[:excess]:
        os.remove(old)
        log(f"   🗑 pruned: {old}")

def main():
    ensure_dir(BACKUP_DIR)
    ts = time.strftime("%Y%m%d-%H%M%S")

    any_found = False
    for db_path in DBS:
        if not db_path.exists():
            log(f"↷ skip (missing): {db_path}")
            continue

        any_found = True
        stem = db_path.stem
        tmp_snap = BACKUP_DIR / f"{stem}.{ts}.snapshot.db"
        out_gz = BACKUP_DIR / f"{stem}.{ts}.sqlite.gz"

        log(f"✅ backing up: {db_path.name}")
        try:
            snapshot_sqlite(db_path, tmp_snap)
            gzip_file(tmp_snap, out_gz)
            log(f"   📦 wrote: {out_gz}")
        except Exception as e:
            log(f"   ❌ error while backing up {db_path.name}: {e}")
        finally:
            if tmp_snap.exists():
                tmp_snap.unlink(missing_ok=True)

        prune_old(BACKUP_DIR, stem)

    if not any_found:
        log("⚠️  No databases found. Checked:")
        for p in DBS:
            log(f"   - {p}")
        # Exit non-zero so cron can alert you if desired
        sys.exit(1)

if __name__ == "__main__":
    log(f"Working dir: {os.getcwd()}")
    log(f"Project root: {PROJECT_ROOT}")
    log(f"Backup dir  : {BACKUP_DIR}")
    main()
