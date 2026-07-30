#!/usr/bin/env python3
"""Apply or roll back PROTACability attachment-site enrichment tables.

Usage:
    python TOOLS/migrate_protacability_attachment_sites.py --database /tmp/viral_data.db --upgrade
    python TOOLS/migrate_protacability_attachment_sites.py --database /tmp/viral_data.db --downgrade
"""

from __future__ import annotations

import argparse
import sqlite3
from pathlib import Path

from protacability_attachment_schema import apply_schema, rollback_schema


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--database", required=True, help="Explicit SQLite database path. Use a copied DB for testing.")
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--upgrade", action="store_true", help="Create attachment-site tables and indexes.")
    action.add_argument("--downgrade", action="store_true", help="Drop attachment-site tables.")
    parser.add_argument("--dry-run", action="store_true", help="Print the action without changing the database.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    db_path = Path(args.database).expanduser().resolve()
    if not db_path.exists():
        raise SystemExit(f"Database does not exist: {db_path}")

    action = "upgrade" if args.upgrade else "downgrade"
    if args.dry_run:
        print(f"Dry run: would {action} PROTACability attachment-site schema in {db_path}")
        return 0

    with sqlite3.connect(db_path) as conn:
        if args.upgrade:
            apply_schema(conn)
        else:
            rollback_schema(conn)
        conn.commit()

    print(f"Completed {action} for PROTACability attachment-site schema in {db_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
