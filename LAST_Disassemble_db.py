#!/usr/bin/env python3
"""
Safe disassembler/exporter (repo root auto-detected).

- viral_data.db  -> output_csvs/<YYYYmmdd-HHMMSS>/*.csv   (chunked)  [DB is deleted after export]
- users.db       -> users_info/<YYYYmmdd-HHMMSS>/*.csv    (chunked)  [DB is KEPT by default]

Env (optional):
  SANITIZE_USERS_EXPORT=1  -> drop users.password_hash and skip 'sessions'
  DELETE_USERS_DB=1        -> also delete users.db after export
  ROW_LIMIT=500000         -> chunk size for viral_data.db
  USERS_ROW_LIMIT=100000   -> chunk size for users.db
"""

import os
import re
import time
import sqlite3
import pandas as pd
from pathlib import Path
import shutil

# ---------- locate repo root robustly ----------
def find_repo_root(start: Path) -> Path:
    for p in [start, *start.parents]:
        if (p / "app.py").exists() or (p / "templates").is_dir() or (p / "requirements.txt").exists():
            return p
    return start  # fallback

HERE = Path(__file__).resolve().parent
ROOT = find_repo_root(HERE)

VIRAL_DB = ROOT / "viral_data.db"
USERS_DB = ROOT / "users.db"
OUT_VIRAL_BASE = ROOT / "output_csvs"
OUT_USERS_BASE = ROOT / "users_info"
LINKER_IMAGES_DIR = ROOT / "static" / "linker_images"

SANITIZE = os.getenv("SANITIZE_USERS_EXPORT", "0") == "1"
DELETE_USERS = os.getenv("DELETE_USERS_DB", "0") == "1"
ROW_LIMIT = int(os.getenv("ROW_LIMIT", "500000"))
USERS_ROW_LIMIT = int(os.getenv("USERS_ROW_LIMIT", "100000"))

def log(m): print(m, flush=True)

# Windows-safe filename sanitizer (remove illegal/reserved chars and trim)
RESERVED_BASENAMES = {"CON","PRN","AUX","NUL","COM1","COM2","COM3","COM4","COM5","COM6","COM7","COM8","COM9",
                      "LPT1","LPT2","LPT3","LPT4","LPT5","LPT6","LPT7","LPT8","LPT9"}
def sanitize_basename(name: str) -> str:
    # keep alnum, dash, underscore; replace others with _
    safe = re.sub(r"[^A-Za-z0-9_\-]+", "_", name).strip("._ ")
    if not safe:
        safe = "table"
    if safe.upper() in RESERVED_BASENAMES:
        safe = f"_{safe}"
    # cap length to something reasonable
    return safe[:100]

def new_stamp_dir(base: Path) -> Path:
    base.mkdir(parents=True, exist_ok=True)
    stamp = time.strftime("%Y%m%d-%H%M%S")
    out = base / stamp
    out.mkdir(parents=True, exist_ok=True)
    # write a small "LATEST" marker
    try:
        (base / "LATEST.txt").write_text(stamp + "\n")
    except Exception:
        pass
    return out

def list_tables(con, exclude=("sqlite_sequence",)):
    rows = con.execute(
        "SELECT name FROM sqlite_master WHERE type='table' AND name NOT LIKE 'sqlite_%';"
    ).fetchall()
    return [r[0] for r in rows if r[0] not in set(exclude)]

def table_columns(con, table):
    return [r[1] for r in con.execute(f'PRAGMA table_info("{table}")').fetchall()]

def export_table(con, table, out_dir: Path, row_limit: int, drop_cols=None):
    drop_cols = set(drop_cols or [])
    cols_all = table_columns(con, table)
    keep = [c for c in cols_all if c not in drop_cols]
    if not keep:
        log(f"↷ skip {table}: no columns after exclusion")
        return

    # SQL select with quoted column names
    col_list = ", ".join([f'"{c}"' for c in keep])
    sql = f'SELECT {col_list} FROM "{table}"'

    safe_base = sanitize_basename(table)
    log(f"→ Export {table} -> {safe_base}_part*.csv (cols {len(keep)}/{len(cols_all)})")

    part = 0
    for chunk in pd.read_sql_query(sql, con, chunksize=row_limit):
        part += 1
        out = out_dir / f"{safe_base}_part{part}.csv"
        # If Windows/WSL ever complains, write to a temp then atomically move
        tmp = out.with_suffix(".csv.tmp")
        chunk.to_csv(tmp, index=False)
        tmp.replace(out)
        log(f"   ✅ wrote {out} ({len(chunk)} rows)")

def export_db(db_path: Path, out_base: Path, row_limit: int,
              table_exclude=(), column_exclude_map=None, delete_db=False):
    if not db_path.exists():
        log(f"⚠️  DB missing, skip: {db_path}")
        return
    out_dir = new_stamp_dir(out_base)

    con = sqlite3.connect(db_path)
    try:
        for t in list_tables(con, exclude=("sqlite_sequence", *table_exclude)):
            drops = (column_exclude_map or {}).get(t, [])
            export_table(con, t, out_dir, row_limit, drop_cols=drops)
    finally:
        con.close()

    if delete_db:
        try:
            db_path.unlink(missing_ok=True)
            log(f"🗑 deleted DB: {db_path}")
        except Exception as e:
            log(f"⚠️  could not delete {db_path}: {e}")

def main():
    log(f"Script: {Path(__file__).name}")
    log(f"Detected repo root: {ROOT}")
    log(f"Looking for: {VIRAL_DB} and {USERS_DB}")

    # 1) viral_data.db -> output_csvs/<stamp> (delete after)
    if VIRAL_DB.exists():
        export_db(
            db_path=VIRAL_DB,
            out_base=OUT_VIRAL_BASE,
            row_limit=ROW_LIMIT,
            delete_db=True
        )
    else:
        log(f"↷ viral_data.db not found at {VIRAL_DB}")

    # prune generated heavy dir if present
    if LINKER_IMAGES_DIR.exists():
        try:
            shutil.rmtree(LINKER_IMAGES_DIR, ignore_errors=True)
            log(f"🗑 removed {LINKER_IMAGES_DIR}")
        except Exception as e:
            log(f"⚠️  could not remove {LINKER_IMAGES_DIR}: {e}")
    else:
        log(f"↷ no {LINKER_IMAGES_DIR} to delete")

    # 2) users.db -> users_info/<stamp> (keep by default; optional sanitize/delete)
    if USERS_DB.exists():
        table_exclude = ()
        col_exclude = {}
        if SANITIZE:
            table_exclude = (*table_exclude, "sessions")
            col_exclude["users"] = ["password_hash"]
        export_db(
            db_path=USERS_DB,
            out_base=OUT_USERS_BASE,
            row_limit=USERS_ROW_LIMIT,
            table_exclude=table_exclude,
            column_exclude_map=col_exclude,
            delete_db=DELETE_USERS
        )
    else:
        log(f"↷ users.db not found at {USERS_DB}")

    log("✅ Done.")

if __name__ == "__main__":
    main()



# import sqlite3
# import pandas as pd
# import os
# import shutil  # ✅ Added to remove directories

# def export_to_csv_and_delete_db(db_path, output_dir, linker_images_dir="static/linker_images", row_limit=500000):
#     # Ensure the output directory exists
#     if not os.path.exists(output_dir):
#         os.makedirs(output_dir)

#     # Connect to the SQLite database
#     conn = sqlite3.connect(db_path)
#     cursor = conn.cursor()

#     # Get the list of all tables
#     cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
#     tables = cursor.fetchall()

#     for table_name in tables:
#         table_name = table_name[0]
#         # Query all data from the table
#         query = f"SELECT * FROM {table_name};"
#         df = pd.read_sql_query(query, conn)

#         # Split the table into multiple CSVs if it's larger than the row limit
#         num_chunks = (len(df) // row_limit) + 1
        
#         for i, chunk in enumerate(range(0, len(df), row_limit)):
#             chunk_df = df.iloc[chunk:chunk + row_limit]
#             csv_path = f"{output_dir}/{table_name}_part{i+1}.csv"
#             chunk_df.to_csv(csv_path, index=False)
#             print(f"✅ Exported {table_name} chunk {i+1} to {csv_path}")

#     # Close the database connection
#     conn.close()

#     # ✅ Delete the database file
#     if os.path.exists(db_path):
#         os.remove(db_path)
#         print(f"🗑 Deleted the database file: {db_path}")
#     else:
#         print(f"⚠️ Database file not found: {db_path}")

#     # ✅ Remove linker_images directory if it exists
#     if os.path.exists(linker_images_dir):
#         shutil.rmtree(linker_images_dir)  # Recursively delete the folder
#         print(f"🗑 Deleted linker_images directory: {linker_images_dir}")
#     else:
#         print(f"⚠️ linker_images directory not found, skipping.")

# # Example usage
# export_to_csv_and_delete_db('viral_data.db', 'output_csvs')
