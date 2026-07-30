#!/usr/bin/env python3
"""Apply or roll back ligand SMILES graph enrichment tables.

Usage:
    python TOOLS/migrate_ligand_smiles_graphs.py --database /tmp/viral_data.db --upgrade
    python TOOLS/migrate_ligand_smiles_graphs.py --database /tmp/viral_data.db --downgrade
"""

from __future__ import annotations

import argparse
import sqlite3
from pathlib import Path

from ligand_smiles_graph_schema import apply_schema, rollback_schema


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--database", required=True, help="Explicit SQLite database path. Use a copied DB for testing.")
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--upgrade", action="store_true", help="Create graph tables and indexes.")
    action.add_argument("--downgrade", action="store_true", help="Drop graph tables.")
    parser.add_argument("--dry-run", action="store_true", help="Print the action without changing the database.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    db_path = Path(args.database).expanduser().resolve()
    if not db_path.exists():
        raise SystemExit(f"Database does not exist: {db_path}")

    action = "upgrade" if args.upgrade else "downgrade"
    if args.dry_run:
        print(f"Dry run: would {action} ligand SMILES graph schema in {db_path}")
        return 0

    with sqlite3.connect(db_path) as conn:
        if args.upgrade:
            apply_schema(conn)
        else:
            rollback_schema(conn)
        conn.commit()

    print(f"Completed {action} for ligand SMILES graph schema in {db_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

