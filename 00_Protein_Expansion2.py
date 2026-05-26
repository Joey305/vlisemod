#!/usr/bin/env python3

import csv
import os
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Tuple

from Bio.Data.IUPACData import protein_letters_3to1
from Bio.PDB.MMCIFParser import MMCIFParser
try:
    from Bio.PDB.MMCIFParser import FastMMCIFParser
except Exception:
    FastMMCIFParser = None
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB.SASA import ShrakeRupley
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# =========================================================
# CONFIG
# =========================================================

PDB_ROOT = Path("PDB_FILES")
MANIFEST_CSV = PDB_ROOT / "download_manifest.csv"
CACHE_DIR = PDB_ROOT / "_CACHE_MMCIF"

OUTPUT_CSV = PDB_ROOT / "Lysine_ISO.csv"
FAILURE_CSV = PDB_ROOT / "Lysine_ISO_failures.csv"

# Use all but 1 CPU, minimum 1
# Use all but 3 CPUs, minimum 1
# MAX_WORKERS = max(1, (os.cpu_count() or 4) - 3)

MAX_WORKERS = 8
# Parsing / calculation options
USE_FAST_MMCIF_PARSER = True
USE_FIRST_MODEL_ONLY = True
SASA_PROBE_RADIUS = 1.40
SASA_N_POINTS = 100

# pI calculation needs a sequence; below this length we leave it blank
MIN_AA_FOR_PI = 2

# =========================================================
# CONSTANTS
# =========================================================

AA3_TO_1 = {k.upper(): v for k, v in protein_letters_3to1.items()}

OUTPUT_FIELDS = [
    "virus_name",
    "protein_type",
    "pdb_code",
    "target_path",
    "cache_path",
    "model_id",
    "chain_id",
    "chain_length_aa",
    "lys_count",
    "total_sasa_a2",
    "lysine_sasa_a2",
    "lysine_surface_fraction",
    "isoelectric_point",
    "chain_sequence",
]

FAILURE_FIELDS = [
    "pdb_code",
    "cache_path",
    "reason",
]

# =========================================================
# HELPERS
# =========================================================

def get_parser():
    """
    Use FastMMCIFParser when available for speed, otherwise fall back.
    """
    if USE_FAST_MMCIF_PARSER and FastMMCIFParser is not None:
        return FastMMCIFParser(QUIET=True)
    return MMCIFParser(QUIET=True)


def residue_to_one_letter(residue) -> str:
    """
    Convert a standard amino-acid residue to a one-letter code.
    Non-standard residues are skipped upstream.
    """
    resname = residue.get_resname().strip().upper()
    return AA3_TO_1.get(resname, "")


def read_manifest(path: Path) -> List[Dict]:
    with path.open("r", newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        return [row for row in reader]


def pick_model(structure):
    """
    For multi-model entries (e.g. NMR), analyze only the first model.
    """
    try:
        return next(structure.get_models())
    except StopIteration:
        return None


def extract_chain_metrics(structure, pdb_code: str) -> List[Dict]:
    """
    Compute per-chain:
      - total protein SASA
      - total LYS SASA
      - LYS SASA fraction of total protein SASA
      - isoelectric point from observed amino-acid sequence
    """
    model = pick_model(structure)
    if model is None:
        return []

    sr = ShrakeRupley(
        probe_radius=SASA_PROBE_RADIUS,
        n_points=SASA_N_POINTS
    )
    sr.compute(model, level="R")

    chain_rows: List[Dict] = []

    for chain in model:
        sequence_chars: List[str] = []
        aa_count = 0
        lys_count = 0
        total_sasa = 0.0
        lysine_sasa = 0.0

        for residue in chain:
            # Only standard amino-acid residues count toward protein surface/pI
            if not is_aa(residue, standard=True):
                continue

            one_letter = residue_to_one_letter(residue)
            if not one_letter:
                continue

            aa_count += 1
            sequence_chars.append(one_letter)

            residue_sasa = float(getattr(residue, "sasa", 0.0) or 0.0)
            total_sasa += residue_sasa

            if residue.get_resname().strip().upper() == "LYS":
                lys_count += 1
                lysine_sasa += residue_sasa

        sequence = "".join(sequence_chars)

        # Skip chains with no standard amino-acid sequence
        if not sequence:
            continue

        if len(sequence) >= MIN_AA_FOR_PI:
            try:
                pI_value = round(float(ProteinAnalysis(sequence).isoelectric_point()), 4)
            except Exception:
                pI_value = ""
        else:
            pI_value = ""

        lys_fraction = round(lysine_sasa / total_sasa, 6) if total_sasa > 0 else ""

        chain_rows.append(
            {
                "pdb_code": pdb_code,
                "model_id": getattr(model, "id", 0),
                "chain_id": chain.id if str(chain.id).strip() else "_blank_",
                "chain_length_aa": aa_count,
                "lys_count": lys_count,
                "total_sasa_a2": round(total_sasa, 4),
                "lysine_sasa_a2": round(lysine_sasa, 4),
                "lysine_surface_fraction": lys_fraction,
                "isoelectric_point": pI_value,
                "chain_sequence": sequence,
            }
        )

    return chain_rows


def analyze_cached_cif(task: Tuple[str, str]) -> Tuple[str, List[Dict], str]:
    """
    Worker function:
      input  -> (pdb_code, cif_path)
      output -> (pdb_code, chain_rows, error_message)
    """
    pdb_code, cif_path_str = task
    cif_path = Path(cif_path_str)

    try:
        parser = get_parser()
        structure = parser.get_structure(pdb_code, str(cif_path))
        chain_rows = extract_chain_metrics(structure, pdb_code)
        return pdb_code, chain_rows, ""
    except Exception as e:
        return pdb_code, [], f"{type(e).__name__}: {e}"


def read_completed_pdbs(output_csv: Path) -> set:
    """
    Read already-completed PDB codes from an existing output CSV.
    This lets the script resume after a crash or manual stop.
    """
    completed = set()

    if not output_csv.exists() or output_csv.stat().st_size == 0:
        return completed

    try:
        with output_csv.open("r", newline="", encoding="utf-8") as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                pdb_code = row.get("pdb_code", "").strip().upper()
                if pdb_code:
                    completed.add(pdb_code)
    except Exception as e:
        print(f"WARNING: Could not read existing output CSV for resume: {e}")

    return completed


def append_rows_csv(path: Path, fieldnames: List[str], rows: List[Dict]) -> None:
    """
    Append rows to a CSV.
    Writes the header only if the file does not exist yet or is empty.
    """
    if not rows:
        return

    write_header = not path.exists() or path.stat().st_size == 0

    with path.open("a", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        if write_header:
            writer.writeheader()
        writer.writerows(rows)
        fh.flush()
        os.fsync(fh.fileno())


def main():
    if not MANIFEST_CSV.exists():
        raise FileNotFoundError(f"Missing manifest CSV: {MANIFEST_CSV}")

    if not CACHE_DIR.exists():
        raise FileNotFoundError(f"Missing cache directory: {CACHE_DIR}")

    manifest_rows = [
        row for row in read_manifest(MANIFEST_CSV)
        if row.get("status", "").strip() in {"copied", "already_exists"}
    ]

    if not manifest_rows:
        raise ValueError("No usable rows found in download_manifest.csv")

    manifest_by_pdb: Dict[str, List[Dict]] = defaultdict(list)
    for row in manifest_rows:
        pdb_code = row["pdb_code"].strip().upper()
        row["pdb_code"] = pdb_code
        manifest_by_pdb[pdb_code].append(row)

    completed_pdbs = read_completed_pdbs(OUTPUT_CSV)

    unique_tasks: List[Tuple[str, str]] = []
    skipped_count = 0

    for pdb_code in sorted(manifest_by_pdb):
        if pdb_code in completed_pdbs:
            skipped_count += 1
            continue

        cache_file = CACHE_DIR / f"{pdb_code}.cif"
        if cache_file.exists():
            unique_tasks.append((pdb_code, str(cache_file)))
        else:
            fallback_path = manifest_by_pdb[pdb_code][0].get("target_path", "")
            unique_tasks.append((pdb_code, fallback_path))

    print(f"Manifest rows: {len(manifest_rows)}")
    print(f"Unique PDBs in manifest: {len(manifest_by_pdb)}")
    print(f"Already completed PDBs skipped: {skipped_count}")
    print(f"Remaining PDBs to analyze: {len(unique_tasks)}")
    print(f"Using {MAX_WORKERS} worker processes")

    if not unique_tasks:
        print("\nNothing left to process.")
        print(f"Existing output: {OUTPUT_CSV}")
        print(f"Existing failures: {FAILURE_CSV}")
        return

    total = len(unique_tasks)
    completed = 0
    successful_pdbs = 0
    failed_pdbs = 0
    empty_pdbs = 0
    written_rows = 0

    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as executor:
        future_to_pdb = {
            executor.submit(analyze_cached_cif, task): task[0]
            for task in unique_tasks
        }

        for future in as_completed(future_to_pdb):
            completed += 1
            pdb_code = future_to_pdb[future]

            try:
                result_pdb, chain_rows, error_message = future.result()
            except Exception as e:
                result_pdb = pdb_code
                chain_rows = []
                error_message = f"{type(e).__name__}: {e}"

            result_pdb = result_pdb.strip().upper()
            cache_path = str(CACHE_DIR / f"{result_pdb}.cif")

            if error_message:
                failed_pdbs += 1
                print(f"[{completed}/{total}] FAILED: {result_pdb} -> {error_message}")

                append_rows_csv(
                    FAILURE_CSV,
                    FAILURE_FIELDS,
                    [
                        {
                            "pdb_code": result_pdb,
                            "cache_path": cache_path,
                            "reason": error_message,
                        }
                    ],
                )
                continue

            if not chain_rows:
                empty_pdbs += 1
                reason = "No standard protein chains found in analyzed first model"
                print(f"[{completed}/{total}] EMPTY: {result_pdb} -> {reason}")

                append_rows_csv(
                    FAILURE_CSV,
                    FAILURE_FIELDS,
                    [
                        {
                            "pdb_code": result_pdb,
                            "cache_path": cache_path,
                            "reason": reason,
                        }
                    ],
                )
                continue

            output_rows_for_pdb: List[Dict] = []

            for manifest_row in manifest_by_pdb.get(result_pdb, []):
                for chain_row in chain_rows:
                    output_rows_for_pdb.append(
                        {
                            "virus_name": manifest_row.get("virus_name", ""),
                            "protein_type": manifest_row.get("protein_type", ""),
                            "pdb_code": result_pdb,
                            "target_path": manifest_row.get("target_path", ""),
                            "cache_path": cache_path,
                            "model_id": chain_row["model_id"],
                            "chain_id": chain_row["chain_id"],
                            "chain_length_aa": chain_row["chain_length_aa"],
                            "lys_count": chain_row["lys_count"],
                            "total_sasa_a2": chain_row["total_sasa_a2"],
                            "lysine_sasa_a2": chain_row["lysine_sasa_a2"],
                            "lysine_surface_fraction": chain_row["lysine_surface_fraction"],
                            "isoelectric_point": chain_row["isoelectric_point"],
                            "chain_sequence": chain_row["chain_sequence"],
                        }
                    )

            append_rows_csv(OUTPUT_CSV, OUTPUT_FIELDS, output_rows_for_pdb)

            successful_pdbs += 1
            written_rows += len(output_rows_for_pdb)

            print(
                f"[{completed}/{total}] OK: {result_pdb} "
                f"({len(chain_rows)} chains, {len(output_rows_for_pdb)} rows written)"
            )

    print("\nDone.")
    print(f"Chain metrics written/appended to: {OUTPUT_CSV}")
    print(f"Failures written/appended to: {FAILURE_CSV}")
    print(f"Skipped already-completed PDBs: {skipped_count}")
    print(f"Successful PDBs this run: {successful_pdbs}")
    print(f"Empty PDBs this run: {empty_pdbs}")
    print(f"Failed PDBs this run: {failed_pdbs}")
    print(f"Rows written this run: {written_rows}")


if __name__ == "__main__":
    main()