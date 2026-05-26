#!/usr/bin/env python3
"""
JSONmaker.py

Regenerates the frontend ligase JSON file for the E3 Recruiter / Ligandalyzer dataset.

Expected project layout:

    ~/VLISEMOD/
        JSONmaker.py
        Ligases/
            MODULE/
                e3-recruiter-mod/
                    Ligases/
                        CRBN/
                            PDB/
                                4CI1_EF2.pdb
                                ...
                        VHL/
                            PDB/
                                ...
        static/
            data/
                ligases.json

Output JSON format:

    {
        "CRBN": [
            "4CI1_EF2.pdb",
            "4CI2_LVY.pdb"
        ],
        "VHL": [
            "3ZRC_L8B_1.pdb"
        ]
    }
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Tuple


# -------------------------------------------------------------------------
# Default paths
# -------------------------------------------------------------------------

SCRIPT_DIR = Path(__file__).resolve().parent

DEFAULT_BASE_DIR = (
    SCRIPT_DIR
    / "Ligases"
    / "MODULE"
    / "e3-recruiter-mod"
    / "Ligases"
)

DEFAULT_OUTPUT_JSON = (
    SCRIPT_DIR
    / "static"
    / "data"
    / "ligases.json"
)


# -------------------------------------------------------------------------
# Core logic
# -------------------------------------------------------------------------

def find_pdb_folder(ligase_path: Path) -> Path | None:
    """
    Find the PDB folder inside a ligase directory.

    Most folders should use:
        LigaseName/PDB/

    This also tolerates lowercase:
        LigaseName/pdb/
    """

    possible_names = ["PDB", "pdb"]

    for name in possible_names:
        candidate = ligase_path / name
        if candidate.is_dir():
            return candidate

    return None


def collect_pdb_files(pdb_folder: Path) -> List[str]:
    """
    Collect .pdb files from a PDB folder.

    Only direct files inside the PDB folder are included.
    Subdirectories are ignored.
    """

    pdb_files = []

    for item in pdb_folder.iterdir():
        if item.is_file() and item.suffix.lower() == ".pdb":
            pdb_files.append(item.name)

    # Remove duplicates just in case, then sort consistently
    pdb_files = sorted(set(pdb_files), key=lambda x: x.lower())

    return pdb_files


def build_ligase_json(base_dir: Path) -> Tuple[Dict[str, List[str]], List[str]]:
    """
    Build the ligase-to-PDB-file mapping.

    Returns:
        ligase_map:
            Dictionary of ligase name -> list of PDB filenames

        skipped:
            Human-readable list of skipped folders/files
    """

    ligase_map: Dict[str, List[str]] = {}
    skipped: List[str] = []

    if not base_dir.exists():
        raise FileNotFoundError(f"Base ligase directory does not exist:\n{base_dir}")

    if not base_dir.is_dir():
        raise NotADirectoryError(f"Base path is not a directory:\n{base_dir}")

    for item in sorted(base_dir.iterdir(), key=lambda p: p.name.lower()):
        # Skip files like .db, .csv, .md, .py, etc.
        if not item.is_dir():
            skipped.append(f"SKIP file: {item.name}")
            continue

        # Skip hidden/system folders
        if item.name.startswith("."):
            skipped.append(f"SKIP hidden/system folder: {item.name}")
            continue

        # Skip Python cache
        if item.name == "__pycache__":
            skipped.append(f"SKIP cache folder: {item.name}")
            continue

        pdb_folder = find_pdb_folder(item)

        if pdb_folder is None:
            skipped.append(f"SKIP no PDB folder: {item.name}")
            continue

        pdb_files = collect_pdb_files(pdb_folder)

        if not pdb_files:
            skipped.append(f"SKIP empty PDB folder: {item.name}")
            continue

        ligase_map[item.name] = pdb_files

    # Sort ligases by name for stable output
    ligase_map = dict(sorted(ligase_map.items(), key=lambda kv: kv[0].lower()))

    return ligase_map, skipped


def write_json(ligase_map: Dict[str, List[str]], output_json: Path) -> None:
    """
    Write ligase map to JSON.
    """

    output_json.parent.mkdir(parents=True, exist_ok=True)

    with output_json.open("w", encoding="utf-8") as f:
        json.dump(ligase_map, f, indent=4, ensure_ascii=False)

    # Add trailing newline for cleaner git diffs
    with output_json.open("a", encoding="utf-8") as f:
        f.write("\n")


def print_summary(
    ligase_map: Dict[str, List[str]],
    base_dir: Path,
    output_json: Path,
    skipped: List[str],
    show_skipped: bool = False,
) -> None:
    """
    Print a clean summary of what was found and written.
    """

    total_ligases = len(ligase_map)
    total_pdbs = sum(len(files) for files in ligase_map.values())

    print("\n============================================================")
    print("JSONmaker.py — Ligase JSON Regenerator")
    print("============================================================")
    print(f"Input ligase directory:")
    print(f"  {base_dir}")
    print()
    print(f"Output JSON:")
    print(f"  {output_json}")
    print()
    print(f"Ligases found: {total_ligases}")
    print(f"Total PDB files: {total_pdbs}")
    print("============================================================\n")

    if ligase_map:
        print("Ligase counts:")
        for ligase, files in ligase_map.items():
            print(f"  {ligase}: {len(files)} PDB files")
        print()

    if show_skipped and skipped:
        print("Skipped items:")
        for msg in skipped:
            print(f"  {msg}")
        print()


# -------------------------------------------------------------------------
# CLI
# -------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Regenerate static/data/ligases.json from the nested Ligases dataset."
    )

    parser.add_argument(
        "--base-dir",
        default=str(DEFAULT_BASE_DIR),
        help=(
            "Path to the ligase dataset folder. "
            "Default: Ligases/MODULE/e3-recruiter-mod/Ligases relative to this script."
        ),
    )

    parser.add_argument(
        "--output",
        default=str(DEFAULT_OUTPUT_JSON),
        help=(
            "Output JSON path. "
            "Default: static/data/ligases.json relative to this script."
        ),
    )

    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Scan and print summary, but do not write the JSON file.",
    )

    parser.add_argument(
        "--show-skipped",
        action="store_true",
        help="Print files/folders that were skipped.",
    )

    return parser.parse_args()


def main() -> int:
    args = parse_args()

    base_dir = Path(args.base_dir).expanduser().resolve()
    output_json = Path(args.output).expanduser().resolve()

    try:
        ligase_map, skipped = build_ligase_json(base_dir)

        print_summary(
            ligase_map=ligase_map,
            base_dir=base_dir,
            output_json=output_json,
            skipped=skipped,
            show_skipped=args.show_skipped,
        )

        if args.dry_run:
            print("Dry run complete. No JSON file was written.\n")
            return 0

        write_json(ligase_map, output_json)

        print("✅ Successfully generated ligases.json")
        print(f"   {output_json}\n")

        return 0

    except Exception as exc:
        print("\n❌ JSONmaker.py failed.")
        print(f"Error: {exc}\n", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())



