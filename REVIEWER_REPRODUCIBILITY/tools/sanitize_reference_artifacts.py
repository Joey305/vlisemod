#!/usr/bin/env python3
"""Remove machine-specific paths from copied release-reference artifacts."""
from __future__ import annotations

import argparse
import re
from pathlib import Path


AUTHOR_WORKSPACE = re.compile("/" + "Users/[^/]+/Documents/VLISEMOD")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    root = Path(args.root).resolve() / "reference"
    for path in root.iterdir():
        if not path.is_file() or path.suffix.lower() not in {".csv", ".json", ".md", ".txt"}:
            continue
        text = path.read_text(encoding="utf-8", errors="replace")
        text = AUTHOR_WORKSPACE.sub("<author-workspace>", text)
        path.write_text(text, encoding="utf-8")


if __name__ == "__main__":
    main()
