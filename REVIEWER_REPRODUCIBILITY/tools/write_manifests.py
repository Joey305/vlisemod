#!/usr/bin/env python3
"""Write package-relative SHA-256 manifests for frozen inputs and scripts."""
from __future__ import annotations

import argparse
import hashlib
from pathlib import Path


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write(root: Path, output: Path, paths: list[Path]) -> None:
    output.write_text("".join(f"{digest(path)}  {path.relative_to(root).as_posix()}\n" for path in sorted(paths)), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    root = Path(args.root).resolve()
    write(root, root / "INPUT_SHA256.txt", [
        root / "manifests/FROZEN_CIF_CORPUS_MANIFEST.csv",
        root / "inputs/chemistry/frozen_component_chemistry.csv",
        root / "inputs/mapping/frozen_mapping_remediation_registry.csv",
        *sorted((root / "fixture").rglob("*.cif")),
    ])
    write(root, root / "SCRIPT_SHA256.txt", [root / "reproduce.py", root / "tools/download_inputs.py", *sorted((root / "pipeline").glob("*.py"))])


if __name__ == "__main__":
    main()
