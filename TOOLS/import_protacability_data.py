#!/usr/bin/env python3
import sqlite3
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DB_PATH = ROOT / "viral_data.db"
PDB_DIR = ROOT / "PDB_FILES"

CSV_TABLES = {
    "protacability_assessment": {"path": PDB_DIR / "PROTACability_Assessment.csv", "required": True},
    "protacability_lysine_proximity": {"path": PDB_DIR / "PROTACability_Lysine_Ligand_Proximity.csv", "required": True},
    "protacability_ligand_inventory": {"path": PDB_DIR / "PROTACability_Ligand_Inventory.csv", "required": True},
    "protacability_warhead_linkability": {"path": PDB_DIR / "PROTACability_Warhead_Linkability.csv", "required": False},
    "protacability_degrader_readiness": {"path": PDB_DIR / "PROTACability_Degrader_Readiness.csv", "required": False},
}

INDEX_STATEMENTS = [
    "CREATE INDEX IF NOT EXISTS idx_protac_assess_pdb ON protacability_assessment(pdb_code);",
    "CREATE INDEX IF NOT EXISTS idx_protac_assess_virus ON protacability_assessment(virus_name);",
    "CREATE INDEX IF NOT EXISTS idx_protac_assess_protein ON protacability_assessment(protein_type);",
    "CREATE INDEX IF NOT EXISTS idx_protac_assess_score ON protacability_assessment(protacability_proxy_score);",
    "CREATE INDEX IF NOT EXISTS idx_protac_assess_tier ON protacability_assessment(protacability_tier);",
    "CREATE INDEX IF NOT EXISTS idx_protac_lys_pdb_chain ON protacability_lysine_proximity(pdb_code, chain_id);",
    "CREATE INDEX IF NOT EXISTS idx_protac_lig_pdb ON protacability_ligand_inventory(pdb_code);",
    "CREATE INDEX IF NOT EXISTS idx_protac_warhead_pdb ON protacability_warhead_linkability(pdb_code);",
    "CREATE INDEX IF NOT EXISTS idx_protac_warhead_ligand ON protacability_warhead_linkability(pdb_code, ligand_resname, ligand_chain, ligand_residue_id);",
    "CREATE INDEX IF NOT EXISTS idx_protac_warhead_score ON protacability_warhead_linkability(warhead_linkability_score);",
    "CREATE INDEX IF NOT EXISTS idx_protac_readiness_chain ON protacability_degrader_readiness(virus_name, protein_type, pdb_code, chain_id);",
    "CREATE INDEX IF NOT EXISTS idx_protac_readiness_score ON protacability_degrader_readiness(degrader_design_readiness_score);",
]


def sqlite_type_for_series(series):
    if pd.api.types.is_bool_dtype(series):
        return "INTEGER"
    if pd.api.types.is_integer_dtype(series):
        return "INTEGER"
    if pd.api.types.is_float_dtype(series):
        return "REAL"
    return "TEXT"


def table_exists(conn, table_name):
    row = conn.execute(
        "SELECT 1 FROM sqlite_master WHERE type = 'table' AND name = ?",
        (table_name,),
    ).fetchone()
    return row is not None


def existing_columns(conn, table_name):
    return {row[1] for row in conn.execute(f"PRAGMA table_info({table_name})").fetchall()}


def ensure_table_schema(conn, table_name, df):
    if not table_exists(conn, table_name):
        empty_df = df.iloc[0:0]
        empty_df.to_sql(table_name, conn, if_exists="append", index=False)
        return

    current_columns = existing_columns(conn, table_name)
    for column in df.columns:
        if column in current_columns:
            continue
        column_type = sqlite_type_for_series(df[column])
        conn.execute(f'ALTER TABLE "{table_name}" ADD COLUMN "{column}" {column_type}')


def import_table(conn, table_name, csv_path):
    df = pd.read_csv(csv_path, low_memory=False)
    df = df.where(pd.notna(df), None)
    ensure_table_schema(conn, table_name, df)
    conn.execute(f'DELETE FROM "{table_name}"')
    df.to_sql(table_name, conn, if_exists="append", index=False)
    print(f"Importing {csv_path.name} -> {table_name}")
    print(f"  rows: {len(df)}")


def main():
    missing_required = [
        str(meta["path"])
        for meta in CSV_TABLES.values()
        if meta["required"] and not meta["path"].exists()
    ]
    if missing_required:
        raise SystemExit(f"Missing required CSV files: {missing_required}")

    with sqlite3.connect(DB_PATH) as conn:
        for table_name, meta in CSV_TABLES.items():
            csv_path = meta["path"]
            if not csv_path.exists():
                print(f"Skipping missing optional CSV: {csv_path.name}")
                continue
            import_table(conn, table_name, csv_path)

        for statement in INDEX_STATEMENTS:
            conn.execute(statement)
        conn.commit()

    print(f"Imported PROTACability tables into {DB_PATH}")


if __name__ == "__main__":
    main()
