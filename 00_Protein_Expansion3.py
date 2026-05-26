#!/usr/bin/env python3

"""
00_Protein_Expansion3.py

Builds PROTACability-style structural assessment data from the existing
V-LiSEMOD PDB_FILES directory.

Inputs expected:
    PDB_FILES/download_manifest.csv
    PDB_FILES/_CACHE_MMCIF/*.cif

Outputs:
    PDB_FILES/PROTACability_Assessment.csv
    PDB_FILES/PROTACability_Lysine_Ligand_Proximity.csv
    PDB_FILES/PROTACability_Ligand_Inventory.csv
    PDB_FILES/PROTACability_failures.csv

Important scientific framing:
    This script produces structural-priority / hypothesis-generation metrics.
    It does NOT produce a true experimental PROTACability score.
"""

import argparse
import csv
import math
import os
import shutil
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path
from statistics import median
from typing import Any, Dict, List, Optional, Tuple

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
# DEFAULT CONFIG
# =========================================================

PDB_ROOT = Path("PDB_FILES")
MANIFEST_CSV = PDB_ROOT / "download_manifest.csv"
CACHE_DIR = PDB_ROOT / "_CACHE_MMCIF"

ASSESSMENT_CSV = PDB_ROOT / "PROTACability_Assessment.csv"
LYSINE_PROXIMITY_CSV = PDB_ROOT / "PROTACability_Lysine_Ligand_Proximity.csv"
LIGAND_INVENTORY_CSV = PDB_ROOT / "PROTACability_Ligand_Inventory.csv"
FAILURE_CSV = PDB_ROOT / "PROTACability_failures.csv"

DEFAULT_MAX_WORKERS = 6
DEFAULT_SASA_N_POINTS = 100
DEFAULT_SASA_PROBE_RADIUS = 1.40

LYSINE_SURFACE_SASA_THRESHOLD = 30.0
LIGAND_PROXIMITY_THRESHOLD_A = 15.0
STRONG_LIGAND_PROXIMITY_THRESHOLD_A = 8.0

USE_FAST_MMCIF_PARSER = True


# =========================================================
# CONSTANTS
# =========================================================

AA3_TO_1 = {k.upper(): v for k, v in protein_letters_3to1.items()}

COMMON_SOLVENTS_IONS = {
    "HOH", "WAT", "DOD",
    "NA", "CL", "K", "MG", "CA", "ZN", "MN", "FE", "CU", "CO", "NI",
    "SO4", "PO4", "NO3", "ACT", "ACE",
    "GOL", "EDO", "PEG", "PG4", "PGE", "DMS", "DMSO",
    "TRS", "MES", "HEP", "BME",
}

BASIC_AA = {"K", "R", "H"}
ACIDIC_AA = {"D", "E"}
POLAR_AA = {"S", "T", "N", "Q", "Y", "C", "K", "R", "H", "D", "E"}
HYDROPHOBIC_AA = {"A", "V", "I", "L", "M", "F", "W", "P", "G"}

ASSESSMENT_FIELDS = [
    "virus_name",
    "protein_type",
    "pdb_code",
    "chain_id",
    "model_id",
    "chain_length_aa",
    "candidate_ligand_count",
    "candidate_ligand_resnames",
    "lys_count",
    "exposed_lys_count",
    "exposed_lys_fraction",
    "near_ligand_lys_count",
    "near_ligand_exposed_lys_count",
    "min_lys_ligand_distance_a",
    "median_lys_ligand_distance_a",
    "total_sasa_a2",
    "lysine_sasa_a2",
    "lysine_surface_fraction",
    "isoelectric_point",
    "basic_fraction",
    "acidic_fraction",
    "polar_fraction",
    "hydrophobic_fraction",
    "has_candidate_ligand",
    "has_exposed_lysine",
    "has_ligand_proximal_exposed_lysine",
    "linker_docking_site_annotation",
    "protein_ligand_druggability_proxy_score",
    "protacability_proxy_score",
    "protacability_tier",
    "notes",
]

LYSINE_PROXIMITY_FIELDS = [
    "virus_name",
    "protein_type",
    "pdb_code",
    "chain_id",
    "model_id",
    "lys_residue_id",
    "lys_insertion_code",
    "lys_observed_index",
    "lysine_sasa_a2",
    "is_surface_exposed",
    "distance_atom",
    "distance_atom_x",
    "distance_atom_y",
    "distance_atom_z",
    "nearest_ligand_resname",
    "nearest_ligand_chain",
    "nearest_ligand_residue_id",
    "nearest_ligand_insertion_code",
    "nearest_ligand_atom",
    "nearest_ligand_distance_a",
    "is_ligand_proximal",
    "linker_site_class",
]

LIGAND_INVENTORY_FIELDS = [
    "virus_name",
    "protein_type",
    "pdb_code",
    "model_id",
    "ligand_resname",
    "ligand_chain",
    "ligand_residue_id",
    "ligand_insertion_code",
    "ligand_atom_count",
    "ligand_heavy_atom_count",
    "centroid_x",
    "centroid_y",
    "centroid_z",
]

FAILURE_FIELDS = [
    "pdb_code",
    "stage",
    "reason",
]


# =========================================================
# CSV / FILE HELPERS
# =========================================================

def read_csv_rows(path: Path) -> List[Dict[str, str]]:
    if not path.exists():
        return []

    with path.open("r", newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def append_rows_csv(path: Path, fieldnames: List[str], rows: List[Dict[str, Any]]) -> None:
    if not rows:
        return

    write_header = not path.exists() or path.stat().st_size == 0

    with path.open("a", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        if write_header:
            writer.writeheader()
        writer.writerows(rows)
        fh.flush()
        os.fsync(fh.fileno())


def read_completed_pdbs(path: Path) -> set:
    completed = set()

    if not path.exists() or path.stat().st_size == 0:
        return completed

    try:
        with path.open("r", newline="", encoding="utf-8") as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                pdb_code = row.get("pdb_code", "").strip().upper()
                if pdb_code:
                    completed.add(pdb_code)
    except Exception as e:
        print(f"WARNING: Could not read completed PDBs from {path}: {e}")

    return completed


def archive_existing_outputs() -> None:
    output_files = [
        ASSESSMENT_CSV,
        LYSINE_PROXIMITY_CSV,
        LIGAND_INVENTORY_CSV,
        FAILURE_CSV,
    ]

    existing = [p for p in output_files if p.exists()]
    if not existing:
        print("No previous PROTACability output files found to archive.")
        return

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    archive_dir = PDB_ROOT / "_old_protacability_runs" / timestamp
    archive_dir.mkdir(parents=True, exist_ok=True)

    for path in existing:
        target = archive_dir / path.name
        shutil.move(str(path), str(target))
        print(f"Archived existing output: {path} -> {target}")


def load_manifest_by_pdb() -> Dict[str, List[Dict[str, str]]]:
    if not MANIFEST_CSV.exists():
        raise FileNotFoundError(f"Missing manifest CSV: {MANIFEST_CSV}")

    rows = [
        row for row in read_csv_rows(MANIFEST_CSV)
        if row.get("status", "").strip() in {"copied", "already_exists"}
    ]

    manifest_by_pdb: Dict[str, List[Dict[str, str]]] = defaultdict(list)

    for row in rows:
        pdb_code = row.get("pdb_code", "").strip().upper()
        if not pdb_code:
            continue
        row["pdb_code"] = pdb_code
        manifest_by_pdb[pdb_code].append(row)

    return manifest_by_pdb


# =========================================================
# STRUCTURE HELPERS
# =========================================================

def get_parser():
    if USE_FAST_MMCIF_PARSER and FastMMCIFParser is not None:
        return FastMMCIFParser(QUIET=True)
    return MMCIFParser(QUIET=True)


def pick_model(structure):
    try:
        return next(structure.get_models())
    except StopIteration:
        return None


def residue_to_one_letter(residue) -> str:
    resname = residue.get_resname().strip().upper()
    return AA3_TO_1.get(resname, "")


def residue_number_and_icode(residue) -> Tuple[str, str]:
    try:
        _, seq_id, insertion_code = residue.id
        insertion_code = "" if insertion_code == " " else str(insertion_code).strip()
        return str(seq_id), insertion_code
    except Exception:
        return "", ""


def atom_coord(atom) -> Tuple[float, float, float]:
    coord = atom.get_coord()
    return float(coord[0]), float(coord[1]), float(coord[2])


def distance_a(coord1: Tuple[float, float, float], coord2: Tuple[float, float, float]) -> float:
    return math.sqrt(
        (coord1[0] - coord2[0]) ** 2
        + (coord1[1] - coord2[1]) ** 2
        + (coord1[2] - coord2[2]) ** 2
    )


def get_preferred_residue_atom(residue) -> Tuple[str, Optional[Tuple[float, float, float]]]:
    for atom_name in ["NZ", "CA"]:
        if atom_name in residue:
            return atom_name, atom_coord(residue[atom_name])

    for atom in residue:
        return atom.get_name(), atom_coord(atom)

    return "", None


def is_candidate_ligand_residue(residue) -> bool:
    resname = residue.get_resname().strip().upper()

    if not resname or resname in COMMON_SOLVENTS_IONS:
        return False

    if is_aa(residue, standard=True):
        return False

    hetero_flag = residue.id[0]
    if hetero_flag == " ":
        return False

    atom_count = sum(1 for _ in residue.get_atoms())
    if atom_count < 3:
        return False

    return True


def extract_candidate_ligands(model) -> List[Dict[str, Any]]:
    ligands = []

    for chain in model:
        for residue in chain:
            if not is_candidate_ligand_residue(residue):
                continue

            residue_id, insertion_code = residue_number_and_icode(residue)
            atoms = []
            heavy_atoms = []

            for atom in residue.get_atoms():
                atom_name = atom.get_name()
                element = (getattr(atom, "element", "") or "").strip().upper()
                coord = atom_coord(atom)

                atom_data = {
                    "atom_name": atom_name,
                    "element": element,
                    "coord": coord,
                }

                atoms.append(atom_data)

                if element != "H":
                    heavy_atoms.append(atom_data)

            coords = [a["coord"] for a in atoms]
            if coords:
                centroid = (
                    round(sum(c[0] for c in coords) / len(coords), 4),
                    round(sum(c[1] for c in coords) / len(coords), 4),
                    round(sum(c[2] for c in coords) / len(coords), 4),
                )
            else:
                centroid = ("", "", "")

            ligands.append(
                {
                    "ligand_resname": residue.get_resname().strip().upper(),
                    "ligand_chain": chain.id if str(chain.id).strip() else "_blank_",
                    "ligand_residue_id": residue_id,
                    "ligand_insertion_code": insertion_code,
                    "ligand_atom_count": len(atoms),
                    "ligand_heavy_atom_count": len(heavy_atoms),
                    "centroid_x": centroid[0],
                    "centroid_y": centroid[1],
                    "centroid_z": centroid[2],
                    "atoms": heavy_atoms if heavy_atoms else atoms,
                }
            )

    return ligands


def find_nearest_ligand(
    query_coord: Tuple[float, float, float],
    ligands: List[Dict[str, Any]],
) -> Dict[str, Any]:
    best = {
        "nearest_ligand_resname": "",
        "nearest_ligand_chain": "",
        "nearest_ligand_residue_id": "",
        "nearest_ligand_insertion_code": "",
        "nearest_ligand_atom": "",
        "nearest_ligand_distance_a": "",
    }

    best_distance = None

    for ligand in ligands:
        for atom in ligand.get("atoms", []):
            d = distance_a(query_coord, atom["coord"])
            if best_distance is None or d < best_distance:
                best_distance = d
                best = {
                    "nearest_ligand_resname": ligand["ligand_resname"],
                    "nearest_ligand_chain": ligand["ligand_chain"],
                    "nearest_ligand_residue_id": ligand["ligand_residue_id"],
                    "nearest_ligand_insertion_code": ligand["ligand_insertion_code"],
                    "nearest_ligand_atom": atom["atom_name"],
                    "nearest_ligand_distance_a": round(d, 4),
                }

    return best


# =========================================================
# SCORING
# =========================================================

def safe_float(value: Any, default: float = 0.0) -> float:
    try:
        if value in {"", None}:
            return default
        return float(value)
    except Exception:
        return default


def classify_linker_site(has_ligand: bool, is_exposed: bool, nearest_distance: Any) -> str:
    """
    This is a geometry annotation only.
    It should not be interpreted as proof of PROTACability.
    """
    if not has_ligand:
        return "No candidate ligand detected"

    distance = safe_float(nearest_distance, default=9999.0)

    if is_exposed and distance <= STRONG_LIGAND_PROXIMITY_THRESHOLD_A:
        return "Strong ligand-proximal exposed lysine geometry cue"
    if is_exposed and distance <= LIGAND_PROXIMITY_THRESHOLD_A:
        return "Ligand-proximal exposed lysine geometry cue"
    if is_exposed:
        return "Surface-exposed lysine distal from ligand"
    if distance <= LIGAND_PROXIMITY_THRESHOLD_A:
        return "Ligand-proximal low-SASA lysine"

    return "No favorable lysine-linker geometry cue detected"


def score_chain_proxy(
    chain_length: int,
    ligand_count: int,
    lys_count: int,
    exposed_lys_count: int,
    near_ligand_exposed_lys_count: int,
    min_lys_ligand_distance: Any,
    pI_available: bool,
) -> Tuple[int, str, str, int, str]:
    """
    Transparent structural-priority heuristic.

    This is NOT a true experimental PROTACability score.

    Main scoring emphasis:
      - Candidate ligand context
      - Surface-exposed lysine availability
      - Surface lysine fraction
      - Absolute exposed lysine count
      - Reasonable modeled chain length
      - pI availability

    Ligand-proximal lysines are kept as a small secondary geometry cue only.
    """
    score = 0
    notes = []

    if ligand_count > 0:
        score += 30
        notes.append("candidate ligand context present")
    else:
        notes.append("no candidate ligand detected")

    if lys_count > 0:
        exposed_fraction = exposed_lys_count / lys_count

        if exposed_lys_count > 0:
            score += 20
            notes.append("surface-exposed lysines detected")
        else:
            notes.append("lysines detected but low surface exposure by current threshold")

        score += min(20, int(round(exposed_fraction * 20)))
        score += min(10, exposed_lys_count * 2)

        notes.append(f"{exposed_lys_count}/{lys_count} lysines exposed")
    else:
        notes.append("no lysines detected")

    if near_ligand_exposed_lys_count > 0:
        score += 5
        notes.append("ligand-proximal exposed lysine recorded as secondary geometry cue")

    if 50 <= chain_length <= 1500:
        score += 10
        notes.append("chain length in reasonable modeled protein range")
    elif chain_length > 0:
        score += 3
        notes.append("chain length available but outside preferred range")

    if pI_available:
        score += 5
        notes.append("pI available")

    score = max(0, min(100, score))

    if ligand_count == 0 and lys_count == 0:
        tier = "Insufficient structural context"
    elif score >= 70:
        tier = "High structural priority"
    elif score >= 45:
        tier = "Moderate structural priority"
    else:
        tier = "Low structural priority"

    if ligand_count > 0 and exposed_lys_count > 0:
        annotation = "Candidate ligand context with surface-exposed lysine accessibility"
    elif ligand_count > 0:
        annotation = "Candidate ligand context detected, but limited exposed lysine signal"
    elif exposed_lys_count > 0:
        annotation = "Surface-exposed lysines detected without candidate ligand context"
    else:
        annotation = "Limited ligand or lysine accessibility signal"

    druggability_score = score
    note_string = "; ".join(notes)

    return score, tier, annotation, druggability_score, note_string


# =========================================================
# WORKER
# =========================================================

def analyze_one_pdb(
    task: Tuple[str, str, float, int],
) -> Tuple[str, List[Dict], List[Dict], List[Dict], str]:
    pdb_code, cif_path_str, probe_radius, sasa_n_points = task
    cif_path = Path(cif_path_str)

    try:
        parser = get_parser()
        structure = parser.get_structure(pdb_code, str(cif_path))
        model = pick_model(structure)

        if model is None:
            return pdb_code, [], [], [], "No model found"

        sr = ShrakeRupley(
            probe_radius=probe_radius,
            n_points=sasa_n_points,
        )
        sr.compute(model, level="R")

        ligands = extract_candidate_ligands(model)

        ligand_rows = []
        for ligand in ligands:
            ligand_rows.append(
                {
                    "pdb_code": pdb_code,
                    "model_id": getattr(model, "id", 0),
                    "ligand_resname": ligand["ligand_resname"],
                    "ligand_chain": ligand["ligand_chain"],
                    "ligand_residue_id": ligand["ligand_residue_id"],
                    "ligand_insertion_code": ligand["ligand_insertion_code"],
                    "ligand_atom_count": ligand["ligand_atom_count"],
                    "ligand_heavy_atom_count": ligand["ligand_heavy_atom_count"],
                    "centroid_x": ligand["centroid_x"],
                    "centroid_y": ligand["centroid_y"],
                    "centroid_z": ligand["centroid_z"],
                }
            )

        assessment_rows = []
        lysine_rows = []

        for chain in model:
            chain_id = chain.id if str(chain.id).strip() else "_blank_"

            sequence_chars = []
            aa_count = 0
            total_sasa = 0.0

            lys_count = 0
            exposed_lys_count = 0
            near_ligand_lys_count = 0
            near_ligand_exposed_lys_count = 0
            lys_ligand_distances = []

            basic_count = 0
            acidic_count = 0
            polar_count = 0
            hydrophobic_count = 0

            chain_lysine_sasa = 0.0
            observed_index = 0

            for residue in chain:
                if not is_aa(residue, standard=True):
                    continue

                one_letter = residue_to_one_letter(residue)
                if not one_letter:
                    continue

                observed_index += 1
                aa_count += 1
                sequence_chars.append(one_letter)

                if one_letter in BASIC_AA:
                    basic_count += 1
                if one_letter in ACIDIC_AA:
                    acidic_count += 1
                if one_letter in POLAR_AA:
                    polar_count += 1
                if one_letter in HYDROPHOBIC_AA:
                    hydrophobic_count += 1

                residue_sasa = float(getattr(residue, "sasa", 0.0) or 0.0)
                total_sasa += residue_sasa

                if residue.get_resname().strip().upper() != "LYS":
                    continue

                lys_count += 1
                chain_lysine_sasa += residue_sasa

                residue_id, insertion_code = residue_number_and_icode(residue)
                distance_atom, query_coord = get_preferred_residue_atom(residue)
                is_surface_exposed = residue_sasa >= LYSINE_SURFACE_SASA_THRESHOLD

                if query_coord is None:
                    nearest = {
                        "nearest_ligand_resname": "",
                        "nearest_ligand_chain": "",
                        "nearest_ligand_residue_id": "",
                        "nearest_ligand_insertion_code": "",
                        "nearest_ligand_atom": "",
                        "nearest_ligand_distance_a": "",
                    }
                    x, y, z = "", "", ""
                else:
                    nearest = find_nearest_ligand(query_coord, ligands) if ligands else {
                        "nearest_ligand_resname": "",
                        "nearest_ligand_chain": "",
                        "nearest_ligand_residue_id": "",
                        "nearest_ligand_insertion_code": "",
                        "nearest_ligand_atom": "",
                        "nearest_ligand_distance_a": "",
                    }
                    x, y, z = (
                        round(query_coord[0], 4),
                        round(query_coord[1], 4),
                        round(query_coord[2], 4),
                    )

                nearest_distance = nearest.get("nearest_ligand_distance_a", "")
                is_ligand_proximal = (
                    nearest_distance != ""
                    and safe_float(nearest_distance, 9999.0) <= LIGAND_PROXIMITY_THRESHOLD_A
                )

                if nearest_distance != "":
                    lys_ligand_distances.append(safe_float(nearest_distance))

                if is_surface_exposed:
                    exposed_lys_count += 1

                if is_ligand_proximal:
                    near_ligand_lys_count += 1

                if is_surface_exposed and is_ligand_proximal:
                    near_ligand_exposed_lys_count += 1

                linker_site_class = classify_linker_site(
                    has_ligand=len(ligands) > 0,
                    is_exposed=is_surface_exposed,
                    nearest_distance=nearest_distance,
                )

                lysine_rows.append(
                    {
                        "pdb_code": pdb_code,
                        "chain_id": chain_id,
                        "model_id": getattr(model, "id", 0),
                        "lys_residue_id": residue_id,
                        "lys_insertion_code": insertion_code,
                        "lys_observed_index": observed_index,
                        "lysine_sasa_a2": round(residue_sasa, 4),
                        "is_surface_exposed": int(is_surface_exposed),
                        "distance_atom": distance_atom,
                        "distance_atom_x": x,
                        "distance_atom_y": y,
                        "distance_atom_z": z,
                        **nearest,
                        "is_ligand_proximal": int(is_ligand_proximal),
                        "linker_site_class": linker_site_class,
                    }
                )

            if aa_count == 0:
                continue

            sequence = "".join(sequence_chars)

            try:
                pI_value = round(float(ProteinAnalysis(sequence).isoelectric_point()), 4)
            except Exception:
                pI_value = ""

            exposed_fraction = round(exposed_lys_count / lys_count, 6) if lys_count else 0.0
            min_distance = round(min(lys_ligand_distances), 4) if lys_ligand_distances else ""
            median_distance = round(median(lys_ligand_distances), 4) if lys_ligand_distances else ""

            ligand_resnames = sorted({lig["ligand_resname"] for lig in ligands})
            ligand_count = len(ligands)

            lysine_surface_fraction = (
                round(chain_lysine_sasa / total_sasa, 6) if total_sasa > 0 else ""
            )

            score, tier, annotation, druggability_score, score_notes = score_chain_proxy(
                chain_length=aa_count,
                ligand_count=ligand_count,
                lys_count=lys_count,
                exposed_lys_count=exposed_lys_count,
                near_ligand_exposed_lys_count=near_ligand_exposed_lys_count,
                min_lys_ligand_distance=min_distance,
                pI_available=(pI_value != ""),
            )

            assessment_rows.append(
                {
                    "pdb_code": pdb_code,
                    "chain_id": chain_id,
                    "model_id": getattr(model, "id", 0),
                    "chain_length_aa": aa_count,
                    "candidate_ligand_count": ligand_count,
                    "candidate_ligand_resnames": ";".join(ligand_resnames),
                    "lys_count": lys_count,
                    "exposed_lys_count": exposed_lys_count,
                    "exposed_lys_fraction": exposed_fraction,
                    "near_ligand_lys_count": near_ligand_lys_count,
                    "near_ligand_exposed_lys_count": near_ligand_exposed_lys_count,
                    "min_lys_ligand_distance_a": min_distance,
                    "median_lys_ligand_distance_a": median_distance,
                    "total_sasa_a2": round(total_sasa, 4),
                    "lysine_sasa_a2": round(chain_lysine_sasa, 4),
                    "lysine_surface_fraction": lysine_surface_fraction,
                    "isoelectric_point": pI_value,
                    "basic_fraction": round(basic_count / aa_count, 6),
                    "acidic_fraction": round(acidic_count / aa_count, 6),
                    "polar_fraction": round(polar_count / aa_count, 6),
                    "hydrophobic_fraction": round(hydrophobic_count / aa_count, 6),
                    "has_candidate_ligand": int(ligand_count > 0),
                    "has_exposed_lysine": int(exposed_lys_count > 0),
                    "has_ligand_proximal_exposed_lysine": int(
                        near_ligand_exposed_lys_count > 0
                    ),
                    "linker_docking_site_annotation": annotation,
                    "protein_ligand_druggability_proxy_score": druggability_score,
                    "protacability_proxy_score": score,
                    "protacability_tier": tier,
                    "notes": (
                        score_notes
                        + "; structural proxy only, not an experimental PROTACability score"
                    ),
                }
            )

        if not assessment_rows:
            return pdb_code, [], [], ligand_rows, "No standard protein chains found"

        return pdb_code, assessment_rows, lysine_rows, ligand_rows, ""

    except Exception as e:
        return pdb_code, [], [], [], f"{type(e).__name__}: {e}"


# =========================================================
# ENRICHMENT / MAIN
# =========================================================

def enrich_rows_with_manifest(
    pdb_code: str,
    base_rows: List[Dict[str, Any]],
    manifest_by_pdb: Dict[str, List[Dict[str, str]]],
) -> List[Dict[str, Any]]:
    enriched = []

    manifest_entries = manifest_by_pdb.get(pdb_code, [])
    if not manifest_entries:
        manifest_entries = [
            {
                "virus_name": "",
                "protein_type": "",
                "pdb_code": pdb_code,
            }
        ]

    for manifest_row in manifest_entries:
        for row in base_rows:
            new_row = dict(row)
            new_row["virus_name"] = manifest_row.get("virus_name", "")
            new_row["protein_type"] = manifest_row.get("protein_type", "")
            enriched.append(new_row)

    return enriched


def build_tasks(
    manifest_by_pdb: Dict[str, List[Dict[str, str]]],
    resume: bool,
    skip_failures: bool,
    limit: Optional[int],
    probe_radius: float,
    sasa_n_points: int,
) -> List[Tuple[str, str, float, int]]:
    completed_pdbs = read_completed_pdbs(ASSESSMENT_CSV) if resume else set()
    failed_pdbs = read_completed_pdbs(FAILURE_CSV) if skip_failures else set()

    tasks = []

    for pdb_code in sorted(manifest_by_pdb):
        if pdb_code in completed_pdbs:
            continue
        if pdb_code in failed_pdbs:
            continue

        cache_file = CACHE_DIR / f"{pdb_code}.cif"
        if cache_file.exists():
            cif_path = str(cache_file)
        else:
            cif_path = manifest_by_pdb[pdb_code][0].get("target_path", "")

        tasks.append((pdb_code, cif_path, probe_radius, sasa_n_points))

    if limit is not None:
        tasks = tasks[:limit]

    print(f"Already completed PDBs skipped: {len(completed_pdbs)}")
    print(f"Previous failures skipped: {len(failed_pdbs)}")

    return tasks


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate PROTACability-style structural assessment tables."
    )

    parser.add_argument(
        "--fresh",
        action="store_true",
        help="Archive existing PROTACability output CSVs and rerun all PDBs from scratch.",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Resume from existing output CSVs. Ignored if --fresh is used.",
    )
    parser.add_argument(
        "--skip-failures",
        action="store_true",
        help="Skip PDBs already present in PROTACability_failures.csv.",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=DEFAULT_MAX_WORKERS,
        help=f"Number of worker processes. Default: {DEFAULT_MAX_WORKERS}",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Limit number of PDBs for testing.",
    )
    parser.add_argument(
        "--sasa-n-points",
        type=int,
        default=DEFAULT_SASA_N_POINTS,
        help=f"Shrake-Rupley points. Default: {DEFAULT_SASA_N_POINTS}",
    )
    parser.add_argument(
        "--probe-radius",
        type=float,
        default=DEFAULT_SASA_PROBE_RADIUS,
        help=f"SASA probe radius. Default: {DEFAULT_SASA_PROBE_RADIUS}",
    )

    return parser.parse_args()


def main():
    args = parse_args()

    if not CACHE_DIR.exists():
        raise FileNotFoundError(f"Missing cache directory: {CACHE_DIR}")

    if args.fresh:
        archive_existing_outputs()
        resume = False
        skip_failures = False
    else:
        resume = args.resume
        skip_failures = args.skip_failures

    manifest_by_pdb = load_manifest_by_pdb()

    tasks = build_tasks(
        manifest_by_pdb=manifest_by_pdb,
        resume=resume,
        skip_failures=skip_failures,
        limit=args.limit,
        probe_radius=args.probe_radius,
        sasa_n_points=args.sasa_n_points,
    )

    print("\nPROTACability Expansion 3")
    print(f"Manifest PDBs: {len(manifest_by_pdb)}")
    print(f"Remaining PDBs to analyze: {len(tasks)}")
    print(f"Using {args.workers} worker processes")
    print(f"SASA_N_POINTS: {args.sasa_n_points}")
    print(f"SASA probe radius: {args.probe_radius}")
    print(f"Lysine exposure threshold: {LYSINE_SURFACE_SASA_THRESHOLD} A^2")
    print(f"Ligand proximity threshold: {LIGAND_PROXIMITY_THRESHOLD_A} A")
    print(
        "Scoring note: ligand-proximal lysines are secondary geometry cues, "
        "not primary PROTACability determinants."
    )

    if not tasks:
        print("\nNothing left to process.")
        return

    total = len(tasks)
    completed = 0
    success_count = 0
    failure_count = 0
    assessment_row_count = 0
    lysine_row_count = 0
    ligand_row_count = 0

    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        future_to_pdb = {
            executor.submit(analyze_one_pdb, task): task[0]
            for task in tasks
        }

        for future in as_completed(future_to_pdb):
            completed += 1
            pdb_code = future_to_pdb[future]

            try:
                result_pdb, assessment_rows, lysine_rows, ligand_rows, error_message = future.result()
            except Exception as e:
                result_pdb = pdb_code
                assessment_rows = []
                lysine_rows = []
                ligand_rows = []
                error_message = f"{type(e).__name__}: {e}"

            result_pdb = result_pdb.strip().upper()

            if error_message:
                failure_count += 1
                print(f"[{completed}/{total}] FAILED: {result_pdb} -> {error_message}")

                append_rows_csv(
                    FAILURE_CSV,
                    FAILURE_FIELDS,
                    [
                        {
                            "pdb_code": result_pdb,
                            "stage": "protein_expansion3",
                            "reason": error_message,
                        }
                    ],
                )
                continue

            enriched_assessment_rows = enrich_rows_with_manifest(
                result_pdb,
                assessment_rows,
                manifest_by_pdb,
            )

            enriched_lysine_rows = enrich_rows_with_manifest(
                result_pdb,
                lysine_rows,
                manifest_by_pdb,
            )

            enriched_ligand_rows = enrich_rows_with_manifest(
                result_pdb,
                ligand_rows,
                manifest_by_pdb,
            )

            append_rows_csv(ASSESSMENT_CSV, ASSESSMENT_FIELDS, enriched_assessment_rows)
            append_rows_csv(LYSINE_PROXIMITY_CSV, LYSINE_PROXIMITY_FIELDS, enriched_lysine_rows)
            append_rows_csv(LIGAND_INVENTORY_CSV, LIGAND_INVENTORY_FIELDS, enriched_ligand_rows)

            success_count += 1
            assessment_row_count += len(enriched_assessment_rows)
            lysine_row_count += len(enriched_lysine_rows)
            ligand_row_count += len(enriched_ligand_rows)

            print(
                f"[{completed}/{total}] OK: {result_pdb} "
                f"assessment={len(enriched_assessment_rows)}, "
                f"lysines={len(enriched_lysine_rows)}, "
                f"ligands={len(enriched_ligand_rows)}"
            )

    print("\nDone.")
    print(f"Assessment CSV: {ASSESSMENT_CSV}")
    print(f"Lysine proximity CSV: {LYSINE_PROXIMITY_CSV}")
    print(f"Ligand inventory CSV: {LIGAND_INVENTORY_CSV}")
    print(f"Failures CSV: {FAILURE_CSV}")
    print(f"Successful PDBs this run: {success_count}")
    print(f"Failed PDBs this run: {failure_count}")
    print(f"Assessment rows written this run: {assessment_row_count}")
    print(f"Lysine rows written this run: {lysine_row_count}")
    print(f"Ligand rows written this run: {ligand_row_count}")


if __name__ == "__main__":
    main()