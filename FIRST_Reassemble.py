#!/usr/bin/env python3
"""
Reassemble SQLite DBs from the newest export folder (repo-root auto-detected).

- Rebuild viral_data.db from output_csvs/<latest>/
- Rebuild users.db from users_info/<latest>/
Falls back to root-level CSVs if no subfolders are found.
"""

import re
import sqlite3
import pandas as pd
from pathlib import Path
from collections import defaultdict

def find_repo_root(start: Path) -> Path:
    for p in [start, *start.parents]:
        if (p / "app.py").exists() or (p / "templates").is_dir() or (p / "requirements.txt").exists():
            return p
    return start

HERE = Path(__file__).resolve().parent
ROOT = find_repo_root(HERE)

VIRAL_DB = ROOT / "viral_data.db"
USERS_DB = ROOT / "users.db"
OUT_VIRAL_BASE = ROOT / "output_csvs"
OUT_USERS_BASE = ROOT / "users_info"

def log(m): print(m, flush=True)

def newest_subdir(base: Path) -> Path | None:
    if not base.exists():
        return None
    latest_txt = base / "LATEST.txt"
    if latest_txt.exists():
        try:
            stamp = latest_txt.read_text().strip()
            candidate = base / stamp
            if candidate.is_dir():
                return candidate
        except Exception:
            pass
    subs = [p for p in base.iterdir() if p.is_dir()]
    return max(subs, key=lambda p: p.stat().st_mtime) if subs else None

def group_csvs(input_dir: Path) -> dict[str, list[Path]]:
    groups = defaultdict(list)
    if not input_dir or not input_dir.exists():
        return groups
    for p in input_dir.glob("*.csv"):
        m = re.match(r"^(?P<base>.+?)(?:_part(?P<idx>\d+))?\.csv$", p.name)
        base = m.group("base") if m else p.stem
        idx = int(m.group("idx") or 1) if m else 1
        groups[base].append((idx, p))
    return {k: [p for _, p in sorted(v, key=lambda t: t[0])] for k, v in groups.items()}

def import_dir_to_db(db_path: Path, input_dir: Path):
    groups = group_csvs(input_dir)
    if not groups:
        log(f"↷ no CSVs in {input_dir}, skipping {db_path.name}")
        return
    con = sqlite3.connect(db_path)
    try:
        for table, files in groups.items():
            dfs = []
            for f in files:
                df = pd.read_csv(f)
                dfs.append(df)
                log(f"   loaded {f} ({len(df)} rows)")
            merged = pd.concat(dfs, ignore_index=True) if len(dfs) > 1 else dfs[0]
            merged.to_sql(table, con, if_exists="replace", index=False)
            log(f"→ imported {table} into {db_path.name} (rows: {len(merged)})")
    finally:
        con.close()

def rebuild_one(db_path: Path, base: Path, label: str):
    sub = newest_subdir(base)
    chosen = sub if sub else base
    log(f"{label}: using CSVs from {chosen}")
    import_dir_to_db(db_path, chosen)

def main():
    log(f"Script: {Path(__file__).name}")
    log(f"Detected repo root: {ROOT}")

    rebuild_one(VIRAL_DB, OUT_VIRAL_BASE, "viral_data.db")
    rebuild_one(USERS_DB, OUT_USERS_BASE, "users.db")

    log("✅ Done.")

if __name__ == "__main__":
    main()


# import sqlite3
# import pandas as pd
# import os

# def import_from_csv(db_path, input_dir):
#     # Connect to the SQLite database
#     conn = sqlite3.connect(db_path)

#     # Dictionary to store dataframes for merging
#     tables_data = {}

#     # Get all CSV files from the input directory
#     for file_name in os.listdir(input_dir):
#         if file_name.endswith(".csv"):
#             table_name = "_".join(file_name.split("_part")[:-1])  # Remove the part index
#             csv_path = os.path.join(input_dir, file_name)
            
#             # Read the CSV file into a pandas DataFrame
#             df = pd.read_csv(csv_path)

#             # Append to the respective table's DataFrame
#             if table_name not in tables_data:
#                 tables_data[table_name] = []
#             tables_data[table_name].append(df)
#             print(f"Loaded {csv_path} for {table_name}")

#     # Merge and write each table back to the SQLite database
#     for table_name, df_list in tables_data.items():
#         merged_df = pd.concat(df_list, ignore_index=True)
#         merged_df.to_sql(table_name, conn, if_exists='replace', index=False)
#         print(f"Merged and imported {table_name} into the database")

#     # Close the database connection
#     conn.close()

# # Example usage
# import_from_csv('viral_data.db', 'output_csvs')
