#!/usr/bin/env python3
"""Remove machine-specific paths from a frozen CIF manifest."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    args = parser.parse_args()
    path = Path(args.manifest)
    with path.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
        fields = [field for field in (handle.seek(0) or csv.DictReader(handle).fieldnames or []) if field != "absolute_path"]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows({field: row.get(field, "") for field in fields} for row in rows)
    print(f"sanitized {len(rows)} rows")


if __name__ == "__main__":
    main()
