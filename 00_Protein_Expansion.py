#!/usr/bin/env python3

import csv
import os
import re
import shutil
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Tuple
from urllib.request import urlopen, Request
from urllib.error import HTTPError, URLError

# =========================================================
# CONFIG
# =========================================================

INPUT_CSV = "output_csvs/Virus_Proteins_part1.csv"
OUTPUT_ROOT = "PDB_FILES"
CACHE_DIRNAME = "_CACHE_MMCIF"
FAILED_LOG = "failed_downloads.csv"
MANIFEST_LOG = "download_manifest.csv"

OVERWRITE_EXISTING = False
REQUEST_TIMEOUT = 30
MAX_RETRIES = 3

# Use all available CPU cores minus 1, but never below 1
MAX_WORKERS = max(1, (os.cpu_count() or 2) - 1)

# Optional aliases for cleaner folder names
VIRUS_ALIASES = {
    "Human immunodeficiency virus 1": "HIV_1",
    "Severe acute respiratory syndrome coronavirus 2": "SARS_CoV_2",
    "Human papillomavirus 16": "HPV_16",
    "Human papillomavirus 18": "HPV_18",
}

PROTEIN_ALIASES = {
    "Protease": "protease",
    "Reverse Transcriptase": "reverse_transcriptase",
    "Integrase": "integrase",
    "Polymerase": "polymerase",
    "Capsid Protein": "capsid_protein",
    "Envelope Protein/Glycoprotein": "envelope_glycoprotein",
    "Accessory Proteins": "accessory_proteins",
    "Other (Peptides & RNA Complexes)": "other_peptides_rna_complexes",
}

# =========================================================
# HELPERS
# =========================================================

def safe_slug(text: str) -> str:
    """Convert arbitrary text into a safe folder name."""
    text = text.strip()
    text = text.replace("&", "and")
    text = text.replace("/", "_")
    text = re.sub(r"[^A-Za-z0-9._ -]+", "", text)
    text = re.sub(r"\s+", "_", text)
    text = re.sub(r"_+", "_", text)
    return text.strip("._-") or "unknown"


def normalize_virus_name(name: str) -> str:
    return VIRUS_ALIASES.get(name.strip(), safe_slug(name))


def normalize_protein_name(name: str) -> str:
    return PROTEIN_ALIASES.get(name.strip(), safe_slug(name).lower())


def normalize_pdb_code(pdb: str) -> str:
    return pdb.strip().upper()


def guess_has_header(first_row: List[str]) -> bool:
    joined = ",".join(first_row).lower()
    return any(
        token in joined
        for token in ["virus", "pdb", "protein", "virus_name", "protein_type", "pdb_code"]
    )


def read_virus_protein_csv(csv_path: Path) -> List[Tuple[str, str, str]]:
    """
    Reads either:
      - headerless rows: virus,pdb,protein
      - headered rows with flexible column names
    """
    rows_out: List[Tuple[str, str, str]] = []

    with csv_path.open("r", newline="", encoding="utf-8-sig") as fh:
        reader = csv.reader(fh)
        raw_rows = [row for row in reader if row and any(cell.strip() for cell in row)]

    if not raw_rows:
        raise ValueError(f"No rows found in {csv_path}")

    first_row = [c.strip() for c in raw_rows[0]]

    if guess_has_header(first_row):
        header = [h.strip().lower() for h in raw_rows[0]]
        data_rows = raw_rows[1:]

        def find_col(possible_names: List[str]) -> int:
            for i, col in enumerate(header):
                if col in possible_names:
                    return i
            raise ValueError(
                f"Could not find any of these columns in header: {possible_names}\nHeader found: {header}"
            )

        virus_idx = find_col(["virus", "virus_name", "virusname"])
        pdb_idx = find_col(["pdb", "pdb_code", "pdbid", "pdb_id"])
        protein_idx = find_col(["protein", "protein_type", "proteintype", "protein_name"])

        for row in data_rows:
            if len(row) <= max(virus_idx, pdb_idx, protein_idx):
                continue
            virus = row[virus_idx].strip()
            pdb = row[pdb_idx].strip()
            protein = row[protein_idx].strip()
            if virus and pdb and protein:
                rows_out.append((virus, normalize_pdb_code(pdb), protein))
    else:
        for row in raw_rows:
            if len(row) < 3:
                continue
            virus = row[0].strip()
            pdb = row[1].strip()
            protein = row[2].strip()
            if virus and pdb and protein:
                rows_out.append((virus, normalize_pdb_code(pdb), protein))

    return rows_out


def download_file(url: str, destination: Path, timeout: int = 30, retries: int = 3) -> None:
    """
    Download a file with retry logic.
    """
    headers = {
        "User-Agent": "Mozilla/5.0 (compatible; VLISEMOD-downloader/1.0)"
    }

    last_err = None
    for attempt in range(1, retries + 1):
        try:
            req = Request(url, headers=headers)
            with urlopen(req, timeout=timeout) as response, destination.open("wb") as out_f:
                shutil.copyfileobj(response, out_f)
            return
        except (HTTPError, URLError, TimeoutError) as err:
            last_err = err
            print(f"[retry {attempt}/{retries}] Failed: {url} -> {err}")
            time.sleep(1.0 * attempt)

    raise RuntimeError(f"Download failed after {retries} attempts: {url}\nLast error: {last_err}")


def download_one_pdb(pdb: str, cache_dir: Path) -> Tuple[str, str, str]:
    """
    Download one PDB mmCIF into the cache.
    Returns: (pdb, status, message)
    status in {"downloaded", "exists", "failed"}
    """
    cache_file = cache_dir / f"{pdb}.cif"
    url = f"https://files.rcsb.org/download/{pdb}.cif"

    if cache_file.exists() and not OVERWRITE_EXISTING:
        return (pdb, "exists", str(cache_file))

    part_file = cache_dir / f"{pdb}.{threading.get_ident()}.part"

    try:
        download_file(
            url=url,
            destination=part_file,
            timeout=REQUEST_TIMEOUT,
            retries=MAX_RETRIES
        )
        part_file.replace(cache_file)
        return (pdb, "downloaded", str(cache_file))
    except Exception as e:
        if part_file.exists():
            try:
                part_file.unlink()
            except Exception:
                pass
        return (pdb, "failed", str(e))


# =========================================================
# MAIN
# =========================================================

def main():
    input_csv = Path(INPUT_CSV)
    output_root = Path(OUTPUT_ROOT)
    cache_dir = output_root / CACHE_DIRNAME

    output_root.mkdir(parents=True, exist_ok=True)
    cache_dir.mkdir(parents=True, exist_ok=True)

    if not input_csv.exists():
        raise FileNotFoundError(f"Could not find input CSV: {input_csv}")

    records = read_virus_protein_csv(input_csv)
    if not records:
        raise ValueError("No valid virus/PDB/protein records found.")

    print(f"Loaded {len(records)} virus/PDB/protein associations from {input_csv}")

    # Unique PDB cache so the same CIF is only downloaded once
    unique_pdbs = sorted({pdb for _, pdb, _ in records})
    print(f"Unique PDB IDs to download: {len(unique_pdbs)}")
    print(f"Using {MAX_WORKERS} worker threads")

    failed_rows = []
    manifest_rows = []

    downloaded_cache: Dict[str, Path] = {}

    # -----------------------------------------------------
    # Step A: download each unique mmCIF once into cache
    #         in parallel
    # -----------------------------------------------------
    total = len(unique_pdbs)
    completed = 0
    downloaded_count = 0
    exists_count = 0
    failed_count = 0

    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        future_to_pdb = {
            executor.submit(download_one_pdb, pdb, cache_dir): pdb
            for pdb in unique_pdbs
        }

        for future in as_completed(future_to_pdb):
            completed += 1
            pdb = future_to_pdb[future]

            try:
                pdb_code, status, message = future.result()
            except Exception as e:
                pdb_code = pdb
                status = "failed"
                message = str(e)

            if status in {"downloaded", "exists"}:
                downloaded_cache[pdb_code] = cache_dir / f"{pdb_code}.cif"

            if status == "exists":
                exists_count += 1
                print(f"[{completed}/{total}] Cache exists: {pdb_code}")
            elif status == "downloaded":
                downloaded_count += 1
                print(f"[{completed}/{total}] Downloaded: {pdb_code}")
            else:
                failed_count += 1
                print(f"[{completed}/{total}] ERROR downloading {pdb_code}: {message}")
                failed_rows.append({
                    "virus_name": "",
                    "protein_type": "",
                    "pdb_code": pdb_code,
                    "reason": message
                })

    print("\nDownload stage complete.")
    print(f"  Downloaded new: {downloaded_count}")
    print(f"  Already cached: {exists_count}")
    print(f"  Failed:         {failed_count}")

    # -----------------------------------------------------
    # Step B: copy cached CIFs into virus/protein folders
    # -----------------------------------------------------
    copy_count = 0
    already_exists_count = 0
    missing_cache_count = 0

    for virus_name, pdb, protein_type in records:
        virus_folder = normalize_virus_name(virus_name)
        protein_folder = normalize_protein_name(protein_type)

        out_dir = output_root / virus_folder / protein_folder
        out_dir.mkdir(parents=True, exist_ok=True)

        target_file = out_dir / f"{pdb}.cif"

        if target_file.exists() and not OVERWRITE_EXISTING:
            already_exists_count += 1
            manifest_rows.append({
                "virus_name": virus_name,
                "protein_type": protein_type,
                "pdb_code": pdb,
                "target_path": str(target_file),
                "status": "already_exists"
            })
            continue

        cache_file = downloaded_cache.get(pdb)
        if cache_file is None or not cache_file.exists():
            missing_cache_count += 1
            failed_rows.append({
                "virus_name": virus_name,
                "protein_type": protein_type,
                "pdb_code": pdb,
                "reason": "download_missing_in_cache"
            })
            continue

        shutil.copy2(cache_file, target_file)
        copy_count += 1

        manifest_rows.append({
            "virus_name": virus_name,
            "protein_type": protein_type,
            "pdb_code": pdb,
            "target_path": str(target_file),
            "status": "copied"
        })

    print("\nCopy stage complete.")
    print(f"  Copied into virus/protein folders: {copy_count}")
    print(f"  Already existed in target folders: {already_exists_count}")
    print(f"  Missing cache entries:             {missing_cache_count}")

    # -----------------------------------------------------
    # Step C: write logs
    # -----------------------------------------------------
    manifest_path = output_root / MANIFEST_LOG
    with manifest_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["virus_name", "protein_type", "pdb_code", "target_path", "status"]
        )
        writer.writeheader()
        writer.writerows(manifest_rows)

    failed_path = output_root / FAILED_LOG
    with failed_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["virus_name", "protein_type", "pdb_code", "reason"]
        )
        writer.writeheader()
        writer.writerows(failed_rows)

    print("\nDone.")
    print(f"Manifest written to: {manifest_path}")
    print(f"Failures written to: {failed_path}")
    print(f"Output root: {output_root}")


if __name__ == "__main__":
    main()