#!/usr/bin/env python3
"""
01_PROTACability_Warhead_Linkability_Enrichment.py

FAST ligand-centered enrichment layer for the V-LiSEMOD PROTACability workflow.

Run this AFTER 00_Protein_Expansion3.py.

Scientific framing
-------------------
This script deliberately separates:

1) Warhead / ligand linkability
   - Ligand-centered.
   - Evaluates whether the bound ligand has evidence for chemically modifiable,
     solvent-exposed atoms that may tolerate linker attachment while preserving
     binding.

2) Target lysine accessibility
   - Protein-centered.
   - Evaluates whether the target has exposed lysines that may be available for
     ubiquitination after ternary complex formation.
   - This does NOT mean the linker attaches to target lysines.

3) Ternary geometry cue
   - Weak, hypothesis-generating structural cue.
   - Based on exposed target lysines near the ligand-binding region.
   - Not proof of productive degradation.

Inputs expected
---------------
    PDB_FILES/PROTACability_Assessment.csv
    PDB_FILES/PROTACability_Ligand_Inventory.csv
    viral_data.db                                  optional, strongly preferred

Outputs
-------
    PDB_FILES/PROTACability_Warhead_Linkability.csv
    PDB_FILES/PROTACability_Degrader_Readiness.csv
    PDB_FILES/PROTACability_Warhead_Linkability_failures.csv

Why this version is faster
--------------------------
    - Deduplicates ligand inventory before scoring.
    - Bulk-loads SQLite evidence once, instead of repeated per-ligand queries.
    - Uses multiprocessing for RDKit/scoring work.
    - Adds progress logging.
    - Caches RDKit descriptor calculations per worker.

Example runs
------------
    python 01_PROTACability_Warhead_Linkability_Enrichment.py --limit 500 --workers 4
    python 01_PROTACability_Warhead_Linkability_Enrichment.py --workers 8
    python 01_PROTACability_Warhead_Linkability_Enrichment.py --serial --limit 100
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import sqlite3
import sys
import time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Set, Tuple


# =========================================================
# OPTIONAL RDKIT SUPPORT
# =========================================================

try:
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors
except Exception:
    Chem = None
    RDLogger = None
    Descriptors = None
    Lipinski = None
    rdMolDescriptors = None


# =========================================================
# DEFAULT CONFIG
# =========================================================

PDB_ROOT = Path("PDB_FILES")
DEFAULT_DB_PATH = Path("viral_data.db")

WARHEAD_LINKABILITY_CSV = PDB_ROOT / "PROTACability_Warhead_Linkability.csv"
DEGRADER_READINESS_CSV = PDB_ROOT / "PROTACability_Degrader_Readiness.csv"
FAILURE_CSV = PDB_ROOT / "PROTACability_Warhead_Linkability_failures.csv"

MEANINGFUL_CONTACTS = {
    "hbond", "weak_hbond", "ionic", "metal", "xbond", "covalent",
    "hydrophobic", "aromatic", "carbonyl", "polar", "weak_polar",
    "CARBONPI", "CATIONPI", "DONORPI", "HALOGENPI", "METSULPHURPI",
    "vdw", "vdw_clash",
}

STRONG_CONTACTS = {
    "hbond", "ionic", "metal", "xbond", "covalent", "aromatic",
    "CATIONPI", "DONORPI", "HALOGENPI",
}

NON_INFORMATIVE_CONTACTS = {"proximal"}

COMMON_CONTEXT_LIGANDS = {
    "HOH", "WAT", "DOD", "NA", "CL", "K", "MG", "CA", "ZN", "MN", "FE",
    "CU", "CO", "NI", "SO4", "PO4", "NO3", "ACT", "ACE", "GOL", "EDO",
    "PEG", "PG4", "PGE", "DMS", "DMSO", "TRS", "MES", "HEP", "BME",
    "MPD", "IPA", "SCN", "FMT", "MSE", "NH4", "SOR",
}

GLYCAN_CONTEXT_LIGANDS = {
    "NAG", "BMA", "MAN", "FUC", "GAL", "GLC", "SIA", "NDG", "BGC", "GLA",
    "GLCN", "A2G", "GCU", "XYL", "FUL", "FRU", "GME", "G7L", "G7O",
}

WARHEAD_FIELDS = [
    "virus_name",
    "protein_type",
    "pdb_code",
    "model_id",
    "ligand_resname",
    "ligand_chain",
    "ligand_residue_id",
    "ligand_insertion_code",
    "ligand_context_class",
    "source_inventory_row_count",
    "smiles_available",
    "representative_smiles",
    "smiles_source",
    "rdkit_available",
    "rdkit_valid_smiles",
    "mol_weight",
    "tpsa",
    "hbd",
    "hba",
    "rotatable_bonds",
    "heavy_atom_count_from_smiles",
    "pdb_ligand_heavy_atom_count",
    "pdb_to_smiles_mapped_atom_count",
    "functional_group_count",
    "functional_group_types",
    "solvent_exposed_ligand_atom_count",
    "solvent_exposed_mapped_atom_count",
    "meaningful_contact_count",
    "strong_contact_count",
    "contact_atom_count",
    "strong_contact_atom_count",
    "candidate_linker_atom_count",
    "candidate_linker_atom_ids",
    "interaction_preservation_score",
    "warhead_linkability_score",
    "warhead_linkability_tier",
    "warhead_linkability_label",
    "warhead_flags",
    "warhead_notes",
]

READINESS_FIELDS = [
    "virus_name",
    "protein_type",
    "pdb_code",
    "chain_id",
    "model_id",
    "best_ligand_resname",
    "best_ligand_chain",
    "best_ligand_residue_id",
    "protein_structural_priority_score",
    "warhead_linkability_score",
    "target_lysine_accessibility_score",
    "ternary_geometry_cue_score",
    "degrader_design_readiness_score",
    "degrader_design_readiness_tier",
    "evidence_level",
    "best_linker_geometry_class",
    "short_linker_geometry_feasible",
    "medium_linker_geometry_feasible",
    "long_linker_geometry_feasible",
    "exposed_lys_count",
    "lys_count",
    "exposed_lys_fraction",
    "lysine_surface_fraction",
    "min_lys_ligand_distance_a",
    "near_ligand_exposed_lys_count",
    "candidate_ligand_resnames",
    "readiness_flags",
    "readiness_notes",
]

FAILURE_FIELDS = ["stage", "pdb_code", "ligand_resname", "reason"]


# =========================================================
# BASIC HELPERS
# =========================================================

def norm_str(value: Any) -> str:
    return str(value or "").strip()


def norm_upper(value: Any) -> str:
    return norm_str(value).upper()


def safe_float(value: Any, default: float = 0.0) -> float:
    try:
        if value in (None, ""):
            return default
        return float(value)
    except Exception:
        return default


def safe_int(value: Any, default: int = 0) -> int:
    try:
        if value in (None, ""):
            return default
        return int(float(value))
    except Exception:
        return default


def clamp(value: float, lo: float = 0.0, hi: float = 100.0) -> float:
    return max(lo, min(hi, value))


def sorted_join(values: Iterable[Any], sep: str = ";") -> str:
    return sep.join(sorted({norm_str(v) for v in values if norm_str(v)}))


def split_ligands(value: Any) -> List[str]:
    if value is None:
        return []
    raw = str(value).replace("|", ";").replace(",", ";")
    out: List[str] = []
    seen: Set[str] = set()
    for part in raw.split(";"):
        item = part.strip().upper()
        if item and item not in seen:
            out.append(item)
            seen.add(item)
    return out


def classify_ligand_context(ligand_resname: str) -> str:
    lig = norm_upper(ligand_resname)
    if not lig:
        return "no_ligand_context"
    if lig in GLYCAN_CONTEXT_LIGANDS:
        return "glycan_only"
    if lig in COMMON_CONTEXT_LIGANDS:
        return "common_buffer_only"
    return "candidate_small_molecule"


def read_csv_rows(path: Path) -> List[Dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def write_csv(path: Path, fieldnames: Sequence[str], rows: Sequence[Dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = path.with_suffix(path.suffix + ".tmp")
    with tmp_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    tmp_path.replace(path)


def now_s(start: float) -> str:
    elapsed = time.time() - start
    if elapsed < 60:
        return f"{elapsed:.1f}s"
    return f"{elapsed/60:.1f}m"


def print_flush(message: str) -> None:
    print(message, flush=True)


def ligand_key_from_row(row: Dict[str, Any]) -> Tuple[str, str, str, str, str, str]:
    """Physical-ish ligand instance key used for scoring/deduplication."""
    return (
        norm_upper(row.get("pdb_code")),
        norm_upper(row.get("ligand_resname")),
        norm_str(row.get("ligand_chain")),
        norm_str(row.get("ligand_residue_id")),
        norm_str(row.get("ligand_insertion_code")),
        norm_str(row.get("model_id")),
    )


def db_key(pdb_code: Any, ligand: Any, chain: Any = "", residue_id: Any = "") -> Tuple[str, str, str, str]:
    return (norm_upper(pdb_code), norm_upper(ligand), norm_str(chain), norm_str(residue_id))


def db_key_fallbacks(pdb_code: Any, ligand: Any, chain: Any = "", residue_id: Any = "") -> List[Tuple[str, str, str, str]]:
    p = norm_upper(pdb_code)
    l = norm_upper(ligand)
    c = norm_str(chain)
    r = norm_str(residue_id)
    return [
        (p, l, c, r),
        (p, l, c, ""),
        (p, l, "", r),
        (p, l, "", ""),
        ("", l, "", ""),  # PDB Chemical Component Dictionary fallback
    ]


# =========================================================
# SQLITE BULK LOAD HELPERS
# =========================================================

def connect_db(path: Path) -> Optional[sqlite3.Connection]:
    if not path.exists():
        return None
    conn = sqlite3.connect(str(path))
    conn.row_factory = sqlite3.Row
    return conn


def table_exists(conn: sqlite3.Connection, table: str) -> bool:
    row = conn.execute(
        "SELECT 1 FROM sqlite_master WHERE type='table' AND name=?",
        (table,),
    ).fetchone()
    return row is not None


def table_columns(conn: sqlite3.Connection, table: str) -> Set[str]:
    if not table_exists(conn, table):
        return set()
    return {row[1] for row in conn.execute(f"PRAGMA table_info({table})")}


def first_existing(columns: Set[str], candidates: Sequence[str]) -> Optional[str]:
    for col in candidates:
        if col in columns:
            return col
    return None


def safe_select_rows(
    conn: Optional[sqlite3.Connection],
    table: str,
    wanted_columns: Sequence[str],
    pdb_filter: Optional[Set[str]] = None,
    progress_name: Optional[str] = None,
) -> Tuple[List[sqlite3.Row], Set[str]]:
    """Bulk-load a table with flexible column selection and optional PDB filter."""
    if conn is None or not table_exists(conn, table):
        return [], set()

    columns = table_columns(conn, table)
    selected = [c for c in wanted_columns if c in columns]
    if not selected:
        return [], columns

    pdb_col = first_existing(columns, ["pdb_id", "PDB_ID", "pdb_code", "PDB_Code"])
    where_sql = ""
    params: List[Any] = []

    # For speed and memory, restrict large tables to relevant PDB IDs.
    if pdb_filter and pdb_col:
        pdb_values = sorted({norm_upper(x) for x in pdb_filter if norm_upper(x)})
        if pdb_values:
            placeholders = ",".join(["?"] * len(pdb_values))
            where_sql = f" WHERE UPPER(CAST({pdb_col} AS TEXT)) IN ({placeholders})"
            params.extend(pdb_values)

    sql = f"SELECT {', '.join(selected)} FROM {table}{where_sql}"
    try:
        rows = list(conn.execute(sql, params).fetchall())
        if progress_name:
            print_flush(f"Loaded {len(rows):,} rows from {table}")
        return rows, columns
    except sqlite3.Error as exc:
        print_flush(f"WARNING: could not bulk-load {table}: {exc}")
        return [], columns


def row_get(row: Any, names: Sequence[str], default: Any = "") -> Any:
    if row is None:
        return default
    if isinstance(row, dict):
        for name in names:
            if name in row and row[name] not in (None, ""):
                return row[name]
        return default
    keys = set(row.keys()) if hasattr(row, "keys") else set()
    for name in names:
        if name in keys and row[name] not in (None, ""):
            return row[name]
    return default


def add_to_map_all_keys(mapping: Dict[Tuple[str, str, str, str], Any], key: Tuple[str, str, str, str], value: Any, reducer: str) -> None:
    """Add evidence to exact and fallback keys to tolerate incomplete DB rows."""
    p, l, c, r = key
    keys = [(p, l, c, r), (p, l, c, ""), (p, l, "", r), (p, l, "", "")]
    for k in keys:
        if reducer == "set":
            mapping[k].update(value)
        elif reducer == "list":
            mapping[k].extend(value)



def load_component_smiles(path: Optional[Path]) -> Dict[str, str]:
    """Load RCSB/PDB Chemical Component SMILES file.

    Expected tab-delimited format:
        SMILES<TAB>CCD_ID<TAB>component name

    Example:
        COC(=O)O    000    methyl hydrogen carbonate

    Returns CCD_ID -> SMILES. If duplicate CCD IDs exist, the first non-empty
    SMILES is kept.
    """
    if path is None or not path.exists():
        return {}

    out: Dict[str, str] = {}
    bad = 0
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                # Some files may be whitespace separated; keep this as fallback.
                parts = line.split(None, 2)
            if len(parts) < 2:
                bad += 1
                continue
            smiles = parts[0].strip()
            ligand_id = parts[1].strip().upper()
            if smiles and ligand_id and ligand_id not in out:
                out[ligand_id] = smiles
    print_flush(f"Loaded {len(out):,} component SMILES from {path}" + (f" ({bad:,} skipped lines)" if bad else ""))
    return out

def build_evidence_maps(conn: Optional[sqlite3.Connection], ligand_tasks: List[Dict[str, Any]], component_smiles: Optional[Dict[str, str]] = None) -> Dict[str, Any]:
    """Preload SQLite evidence into compact dictionaries keyed by ligand instance."""
    pdb_filter = {task["pdb_code"] for task in ligand_tasks}
    maps: Dict[str, Any] = {
        "smiles": defaultdict(set),
        "smiles_source": defaultdict(set),
        "functional_groups": defaultdict(set),
        "sasa_atom_ids": defaultdict(set),
        "mapped_atom_ids": defaultdict(set),
        "contact_atom_ids": defaultdict(set),
        "strong_contact_atom_ids": defaultdict(set),
        "meaningful_contact_count": defaultdict(int),
        "strong_contact_count": defaultdict(int),
        "tables_loaded": {},
    }

    component_smiles = component_smiles or {}

    if component_smiles:
        # Add ligand-code fallback SMILES from the PDB Chemical Component Dictionary.
        # This is intentionally keyed as (any PDB, ligand code, any chain/residue)
        # so it fills missing local DB SMILES without overwriting exact DB evidence.
        ligands_in_tasks = {task["ligand_resname"] for task in ligand_tasks if task.get("ligand_resname")}
        added = 0
        for ligand in ligands_in_tasks:
            smiles = component_smiles.get(ligand)
            if smiles:
                key = ("", ligand, "", "")
                maps["smiles"][key].add(smiles)
                maps["smiles_source"][key].add("component_smiles")
                added += 1
        print_flush(f"Component SMILES fallback matched {added:,} ligand codes in ligand inventory")

    if conn is None:
        return maps

    print_flush("\nBulk-loading SQLite evidence...")

    # Ligand_Atoms_Smiles
    rows, cols = safe_select_rows(
        conn,
        "Ligand_Atoms_Smiles",
        [
            "pdb_id", "PDB_ID", "pdb_code", "PDB_Code",
            "ligand", "Ligand", "ligand_resname",
            "chain", "Chain", "ligand_chain",
            "ligand_id", "Ligand_ID", "residue_id", "ligand_residue_id",
            "smiles", "SMILES",
        ],
        pdb_filter=pdb_filter,
        progress_name="Ligand_Atoms_Smiles",
    )
    maps["tables_loaded"]["Ligand_Atoms_Smiles"] = len(rows)
    for r in rows:
        key = db_key(
            row_get(r, ["pdb_id", "PDB_ID", "pdb_code", "PDB_Code"]),
            row_get(r, ["ligand", "Ligand", "ligand_resname"]),
            row_get(r, ["chain", "Chain", "ligand_chain"]),
            row_get(r, ["ligand_id", "Ligand_ID", "residue_id", "ligand_residue_id"]),
        )
        smiles = norm_str(row_get(r, ["smiles", "SMILES"]))
        if key[0] and key[1] and smiles:
            add_to_map_all_keys(maps["smiles"], key, {smiles}, "set")
            add_to_map_all_keys(maps["smiles_source"], key, {"viral_data.db"}, "set")

    # Functional_Group_Atoms
    rows, cols = safe_select_rows(
        conn,
        "Functional_Group_Atoms",
        [
            "pdb_id", "PDB_ID", "pdb_code", "PDB_Code",
            "ligand", "Ligand", "ligand_resname",
            "chain", "Chain", "ligand_chain",
            "ligand_id", "Ligand_ID", "residue_id", "ligand_residue_id",
            "functional_group", "Functional_Group", "functional_group_name",
        ],
        pdb_filter=pdb_filter,
        progress_name="Functional_Group_Atoms",
    )
    maps["tables_loaded"]["Functional_Group_Atoms"] = len(rows)
    for r in rows:
        key = db_key(
            row_get(r, ["pdb_id", "PDB_ID", "pdb_code", "PDB_Code"]),
            row_get(r, ["ligand", "Ligand", "ligand_resname"]),
            row_get(r, ["chain", "Chain", "ligand_chain"]),
            row_get(r, ["ligand_id", "Ligand_ID", "residue_id", "ligand_residue_id"]),
        )
        fg = norm_str(row_get(r, ["functional_group", "Functional_Group", "functional_group_name"]))
        if key[0] and key[1] and fg:
            add_to_map_all_keys(maps["functional_groups"], key, {fg}, "set")

    # RUPLEY_SASA_DATA
    rows, cols = safe_select_rows(
        conn,
        "RUPLEY_SASA_DATA",
        [
            "pdb_id", "PDB_ID", "pdb_code", "PDB_Code",
            "ligand", "Ligand", "ligand_resname",
            "chain", "Chain", "ligand_chain",
            "ligand_id", "Ligand_ID", "residue_id", "ligand_residue_id",
            "atom_id", "Atom_ID", "ligand_atom_id",
        ],
        pdb_filter=pdb_filter,
        progress_name="RUPLEY_SASA_DATA",
    )
    maps["tables_loaded"]["RUPLEY_SASA_DATA"] = len(rows)
    for r in rows:
        key = db_key(
            row_get(r, ["pdb_id", "PDB_ID", "pdb_code", "PDB_Code"]),
            row_get(r, ["ligand", "Ligand", "ligand_resname"]),
            row_get(r, ["chain", "Chain", "ligand_chain"]),
            row_get(r, ["ligand_id", "Ligand_ID", "residue_id", "ligand_residue_id"]),
        )
        atom = norm_str(row_get(r, ["atom_id", "Atom_ID", "ligand_atom_id"]))
        if key[0] and key[1] and atom:
            add_to_map_all_keys(maps["sasa_atom_ids"], key, {atom}, "set")

    # SMILES_MAP_PDB
    rows, cols = safe_select_rows(
        conn,
        "SMILES_MAP_PDB",
        [
            "pdb_id", "PDB_ID", "pdb_code", "PDB_Code",
            "ligand", "Ligand", "ligand_resname",
            "chain", "Chain", "ligand_chain",
            "ligand_id", "Ligand_ID", "residue_id", "ligand_residue_id",
            "atom_id", "Atom_ID", "ligand_atom_id",
            "smiles_atom_index",
        ],
        pdb_filter=pdb_filter,
        progress_name="SMILES_MAP_PDB",
    )
    maps["tables_loaded"]["SMILES_MAP_PDB"] = len(rows)
    for r in rows:
        key = db_key(
            row_get(r, ["pdb_id", "PDB_ID", "pdb_code", "PDB_Code"]),
            row_get(r, ["ligand", "Ligand", "ligand_resname"]),
            row_get(r, ["chain", "Chain", "ligand_chain"]),
            row_get(r, ["ligand_id", "Ligand_ID", "residue_id", "ligand_residue_id"]),
        )
        atom = norm_str(row_get(r, ["atom_id", "Atom_ID", "ligand_atom_id"]))
        if key[0] and key[1] and atom:
            add_to_map_all_keys(maps["mapped_atom_ids"], key, {atom}, "set")

    # Arpeggio_Contacts_Data
    rows, cols = safe_select_rows(
        conn,
        "Arpeggio_Contacts_Data",
        [
            "pdb_id", "PDB_ID", "pdb_code", "PDB_Code",
            "ligand", "Ligand", "ligand_resname",
            "chain", "Chain", "ligand_chain",
            "ligand_id", "Ligand_ID", "residue_id", "ligand_residue_id",
            "Contact", "contact", "Interaction", "interaction",
            "atom_id", "Atom_ID", "ligand_atom_id",
        ],
        pdb_filter=pdb_filter,
        progress_name="Arpeggio_Contacts_Data",
    )
    maps["tables_loaded"]["Arpeggio_Contacts_Data"] = len(rows)
    meaningful_lower = {x.lower() for x in MEANINGFUL_CONTACTS}
    strong_lower = {x.lower() for x in STRONG_CONTACTS}
    non_info_lower = {x.lower() for x in NON_INFORMATIVE_CONTACTS}

    for r in rows:
        key = db_key(
            row_get(r, ["pdb_id", "PDB_ID", "pdb_code", "PDB_Code"]),
            row_get(r, ["ligand", "Ligand", "ligand_resname"]),
            row_get(r, ["chain", "Chain", "ligand_chain"]),
            row_get(r, ["ligand_id", "Ligand_ID", "residue_id", "ligand_residue_id"]),
        )
        contact = norm_str(row_get(r, ["Contact", "contact", "Interaction", "interaction"]))
        contact_l = contact.lower()
        atom = norm_str(row_get(r, ["atom_id", "Atom_ID", "ligand_atom_id"]))
        if not key[0] or not key[1] or not contact_l or contact_l in non_info_lower:
            continue

        for k in db_key_fallbacks(*key):
            if contact_l in meaningful_lower:
                maps["meaningful_contact_count"][k] += 1
                if atom:
                    maps["contact_atom_ids"][k].add(atom)
            if contact_l in strong_lower:
                maps["strong_contact_count"][k] += 1
                if atom:
                    maps["strong_contact_atom_ids"][k].add(atom)

    print_flush("SQLite evidence load complete.")
    return maps


def evidence_for_task(task: Dict[str, Any], maps: Dict[str, Any]) -> Dict[str, Any]:
    p = task["pdb_code"]
    l = task["ligand_resname"]
    c = task["ligand_chain"]
    r = task["ligand_residue_id"]

    def lookup_set(name: str) -> Set[str]:
        table = maps.get(name, {})
        for k in db_key_fallbacks(p, l, c, r):
            val = table.get(k)
            if val:
                return set(val)
        return set()

    def lookup_int(name: str) -> int:
        table = maps.get(name, {})
        for k in db_key_fallbacks(p, l, c, r):
            val = table.get(k)
            if val:
                return int(val)
        return 0

    smiles_values = lookup_set("smiles")
    smiles_sources = lookup_set("smiles_source")
    fg_types = lookup_set("functional_groups")
    exposed_atom_ids = lookup_set("sasa_atom_ids")
    mapped_atom_ids = lookup_set("mapped_atom_ids")
    contact_atom_ids = lookup_set("contact_atom_ids")
    strong_contact_atom_ids = lookup_set("strong_contact_atom_ids")
    meaningful_contact_count = lookup_int("meaningful_contact_count")
    strong_contact_count = lookup_int("strong_contact_count")

    mapped_exposed_atom_ids = exposed_atom_ids & mapped_atom_ids if mapped_atom_ids else set(exposed_atom_ids)
    candidate_linker_atoms = mapped_exposed_atom_ids - strong_contact_atom_ids if mapped_exposed_atom_ids else exposed_atom_ids - strong_contact_atom_ids

    return {
        "representative_smiles": sorted(smiles_values)[0] if smiles_values else "",
        "smiles_source": sorted_join(smiles_sources) if smiles_sources else "",
        "smiles_available": int(bool(smiles_values)),
        "functional_group_count": len(fg_types),
        "functional_group_types": sorted_join(fg_types),
        "solvent_exposed_ligand_atom_ids": exposed_atom_ids,
        "solvent_exposed_ligand_atom_count": len(exposed_atom_ids),
        "pdb_to_smiles_mapped_atom_ids": mapped_atom_ids,
        "pdb_to_smiles_mapped_atom_count": len(mapped_atom_ids),
        "solvent_exposed_mapped_atom_count": len(mapped_exposed_atom_ids),
        "meaningful_contact_count": meaningful_contact_count,
        "strong_contact_count": strong_contact_count,
        "contact_atom_count": len(contact_atom_ids),
        "strong_contact_atom_count": len(strong_contact_atom_ids),
        "candidate_linker_atom_ids": candidate_linker_atoms,
        "candidate_linker_atom_count": len(candidate_linker_atoms),
    }


# =========================================================
# DEDUPLICATION
# =========================================================

def deduplicate_ligand_inventory(ligand_rows: List[Dict[str, str]], limit: Optional[int] = None) -> Tuple[List[Dict[str, Any]], Dict[Tuple[str, str, str, str, str, str], List[Dict[str, str]]]]:
    grouped: Dict[Tuple[str, str, str, str, str, str], List[Dict[str, str]]] = defaultdict(list)
    for row in ligand_rows:
        key = ligand_key_from_row(row)
        if not key[0] or not key[1]:
            continue
        grouped[key].append(row)

    tasks: List[Dict[str, Any]] = []
    for key, rows in grouped.items():
        pdb_code, ligand, chain, residue_id, insertion_code, model_id = key
        virus_names = sorted({norm_str(r.get("virus_name")) for r in rows if norm_str(r.get("virus_name"))})
        protein_types = sorted({norm_str(r.get("protein_type")) for r in rows if norm_str(r.get("protein_type"))})
        max_heavy = max((safe_int(r.get("ligand_heavy_atom_count")) for r in rows), default=0)

        tasks.append({
            "pdb_code": pdb_code,
            "ligand_resname": ligand,
            "ligand_chain": chain,
            "ligand_residue_id": residue_id,
            "ligand_insertion_code": insertion_code,
            "model_id": model_id,
            "virus_name": ";".join(virus_names),
            "protein_type": ";".join(protein_types),
            "source_inventory_row_count": len(rows),
            "pdb_ligand_heavy_atom_count": max_heavy,
        })

    tasks.sort(key=lambda x: (x["pdb_code"], x["ligand_resname"], x["ligand_chain"], x["ligand_residue_id"]))

    if limit is not None:
        tasks = tasks[:limit]

    return tasks, grouped


# =========================================================
# RDKit and scoring
# =========================================================

_RDKIT_CACHE: Dict[str, Dict[str, Any]] = {}


def setup_worker(rdkit_warnings: bool = False) -> None:
    if RDLogger is not None and not rdkit_warnings:
        RDLogger.DisableLog("rdApp.*")


def rdkit_descriptors(smiles: str) -> Dict[str, Any]:
    global _RDKIT_CACHE

    if smiles in _RDKIT_CACHE:
        return dict(_RDKIT_CACHE[smiles])

    result = {
        "rdkit_available": int(Chem is not None),
        "rdkit_valid_smiles": 0,
        "mol_weight": "",
        "tpsa": "",
        "hbd": "",
        "hba": "",
        "rotatable_bonds": "",
        "heavy_atom_count_from_smiles": "",
    }

    if not smiles or Chem is None:
        _RDKIT_CACHE[smiles] = result
        return dict(result)

    try:
        mol = Chem.MolFromSmiles(smiles)
    except Exception:
        mol = None

    if mol is None:
        _RDKIT_CACHE[smiles] = result
        return dict(result)

    result["rdkit_valid_smiles"] = 1
    try:
        result["mol_weight"] = round(float(Descriptors.MolWt(mol)), 4) if Descriptors else ""
        result["tpsa"] = round(float(rdMolDescriptors.CalcTPSA(mol)), 4) if rdMolDescriptors else ""
        result["hbd"] = int(Lipinski.NumHDonors(mol)) if Lipinski else ""
        result["hba"] = int(Lipinski.NumHAcceptors(mol)) if Lipinski else ""
        result["rotatable_bonds"] = int(Lipinski.NumRotatableBonds(mol)) if Lipinski else ""
        result["heavy_atom_count_from_smiles"] = int(mol.GetNumHeavyAtoms())
    except Exception:
        # Keep valid flag but leave descriptors blank if a descriptor fails.
        pass

    _RDKIT_CACHE[smiles] = result
    return dict(result)


def interaction_preservation_score(evidence: Dict[str, Any]) -> float:
    exposed = safe_int(evidence.get("solvent_exposed_ligand_atom_count"))
    candidates = safe_int(evidence.get("candidate_linker_atom_count"))
    strong_contact_atoms = safe_int(evidence.get("strong_contact_atom_count"))

    if exposed <= 0:
        return 0.0

    candidate_fraction = candidates / exposed
    score = 100.0 * candidate_fraction

    if strong_contact_atoms >= exposed and candidates == 0:
        score -= 30
    elif strong_contact_atoms > 0:
        score -= min(20, strong_contact_atoms * 4)

    return round(clamp(score), 2)


def score_warhead_linkability(task: Dict[str, Any], evidence: Dict[str, Any], rdkit_info: Dict[str, Any]) -> Tuple[float, str, str, str, str]:
    ligand = norm_upper(task.get("ligand_resname"))
    context_class = classify_ligand_context(ligand)
    flags: List[str] = []
    notes: List[str] = []
    score = 0.0

    if context_class == "candidate_small_molecule":
        score += 10
        notes.append("candidate small-molecule context")
    elif context_class == "glycan_only":
        flags.append("glycan_context_not_warhead_like")
        notes.append("glycan context; weak warhead evidence")
    elif context_class == "common_buffer_only":
        flags.append("common_buffer_or_crystal_ligand")
        notes.append("common buffer/crystallization ligand; not treated as warhead evidence")

    if evidence.get("smiles_available"):
        score += 12
        source = norm_str(evidence.get("smiles_source"))
        if source:
            notes.append(f"SMILES available from {source}")
        else:
            notes.append("SMILES available")
    else:
        flags.append("missing_smiles")

    if rdkit_info.get("rdkit_valid_smiles"):
        score += 8
        notes.append("valid RDKit molecule")
    elif evidence.get("smiles_available"):
        flags.append("invalid_smiles_for_rdkit")

    mapped = safe_int(evidence.get("pdb_to_smiles_mapped_atom_count"))
    if mapped > 0:
        score += min(15, 5 + mapped * 0.5)
        notes.append("PDB-to-SMILES atom mapping available")
    else:
        flags.append("missing_pdb_to_smiles_mapping")

    exposed = safe_int(evidence.get("solvent_exposed_ligand_atom_count"))
    exposed_mapped = safe_int(evidence.get("solvent_exposed_mapped_atom_count"))
    if exposed_mapped > 0:
        score += min(20, 8 + exposed_mapped * 3)
        notes.append("mapped solvent-exposed ligand atoms detected")
    elif exposed > 0:
        score += min(12, 4 + exposed * 2)
        notes.append("solvent-exposed ligand atoms detected, mapping incomplete")
    else:
        flags.append("no_solvent_exposed_ligand_atoms_detected")

    fg_count = safe_int(evidence.get("functional_group_count"))
    if fg_count > 0:
        score += min(15, 5 + fg_count * 2)
        notes.append("functional group annotations available")
    else:
        flags.append("no_functional_group_annotations")

    candidate_atoms = safe_int(evidence.get("candidate_linker_atom_count"))
    if candidate_atoms >= 5:
        score += 20
        notes.append("multiple candidate exposed non-critical linker atoms")
    elif candidate_atoms >= 2:
        score += 16
        notes.append("several candidate exposed non-critical linker atoms")
    elif candidate_atoms == 1:
        score += 10
        notes.append("one candidate exposed non-critical linker atom")
    else:
        flags.append("no_candidate_linker_atoms_after_contact_filtering")

    preservation = interaction_preservation_score(evidence)
    score += preservation * 0.10
    if preservation < 25:
        flags.append("high_risk_of_disrupting_binding_contacts")

    # Mild chemistry sanity penalties so scores distribute better.
    if evidence.get("smiles_available") and not rdkit_info.get("rdkit_valid_smiles"):
        score -= 10
    if not evidence.get("smiles_available"):
        score -= 8
    if safe_int(evidence.get("candidate_linker_atom_count")) == 0:
        score -= 10

    # Strong cap for contexts that are not warhead-like.
    if context_class in {"glycan_only", "common_buffer_only"}:
        score = min(score, 35)

    score = round(clamp(score), 2)

    if score >= 75:
        tier = "High warhead-linkability evidence"
        label = "Linkerable warhead candidate"
    elif score >= 55:
        tier = "Moderate warhead-linkability evidence"
        label = "Plausible warhead candidate"
    elif score >= 35:
        tier = "Exploratory warhead-linkability evidence"
        label = "Needs manual linker-site review"
    else:
        tier = "Weak warhead-linkability evidence"
        label = "Poor or incomplete warhead evidence"

    notes.append("ligand-centered linkability; target lysines are evaluated separately")
    return score, tier, label, ";".join(dict.fromkeys(flags)), "; ".join(notes)


def score_one_ligand(payload: Tuple[Dict[str, Any], Dict[str, Any]]) -> Dict[str, Any]:
    task, evidence = payload
    smiles = evidence.get("representative_smiles", "")
    rdkit_info = rdkit_descriptors(smiles)
    score, tier, label, flags, notes = score_warhead_linkability(task, evidence, rdkit_info)

    return {
        "virus_name": task.get("virus_name", ""),
        "protein_type": task.get("protein_type", ""),
        "pdb_code": task.get("pdb_code", ""),
        "model_id": task.get("model_id", ""),
        "ligand_resname": task.get("ligand_resname", ""),
        "ligand_chain": task.get("ligand_chain", ""),
        "ligand_residue_id": task.get("ligand_residue_id", ""),
        "ligand_insertion_code": task.get("ligand_insertion_code", ""),
        "ligand_context_class": classify_ligand_context(task.get("ligand_resname", "")),
        "source_inventory_row_count": task.get("source_inventory_row_count", 1),
        "smiles_available": evidence.get("smiles_available", 0),
        "representative_smiles": evidence.get("representative_smiles", ""),
        "smiles_source": evidence.get("smiles_source", ""),
        "rdkit_available": rdkit_info.get("rdkit_available", 0),
        "rdkit_valid_smiles": rdkit_info.get("rdkit_valid_smiles", 0),
        "mol_weight": rdkit_info.get("mol_weight", ""),
        "tpsa": rdkit_info.get("tpsa", ""),
        "hbd": rdkit_info.get("hbd", ""),
        "hba": rdkit_info.get("hba", ""),
        "rotatable_bonds": rdkit_info.get("rotatable_bonds", ""),
        "heavy_atom_count_from_smiles": rdkit_info.get("heavy_atom_count_from_smiles", ""),
        "pdb_ligand_heavy_atom_count": task.get("pdb_ligand_heavy_atom_count", ""),
        "pdb_to_smiles_mapped_atom_count": evidence.get("pdb_to_smiles_mapped_atom_count", 0),
        "functional_group_count": evidence.get("functional_group_count", 0),
        "functional_group_types": evidence.get("functional_group_types", ""),
        "solvent_exposed_ligand_atom_count": evidence.get("solvent_exposed_ligand_atom_count", 0),
        "solvent_exposed_mapped_atom_count": evidence.get("solvent_exposed_mapped_atom_count", 0),
        "meaningful_contact_count": evidence.get("meaningful_contact_count", 0),
        "strong_contact_count": evidence.get("strong_contact_count", 0),
        "contact_atom_count": evidence.get("contact_atom_count", 0),
        "strong_contact_atom_count": evidence.get("strong_contact_atom_count", 0),
        "candidate_linker_atom_count": evidence.get("candidate_linker_atom_count", 0),
        "candidate_linker_atom_ids": sorted_join(evidence.get("candidate_linker_atom_ids", [])),
        "interaction_preservation_score": interaction_preservation_score(evidence),
        "warhead_linkability_score": score,
        "warhead_linkability_tier": tier,
        "warhead_linkability_label": label,
        "warhead_flags": flags,
        "warhead_notes": notes,
    }


def build_warhead_rows(
    ligand_tasks: List[Dict[str, Any]],
    evidence_maps: Dict[str, Any],
    workers: int,
    chunksize: int,
    serial: bool,
    progress_every: int,
    rdkit_warnings: bool,
) -> List[Dict[str, Any]]:
    print_flush("\nPreparing ligand scoring payloads...")
    payloads = [(task, evidence_for_task(task, evidence_maps)) for task in ligand_tasks]
    total = len(payloads)

    smiles_count = sum(1 for _, ev in payloads if ev.get("smiles_available"))
    mapping_count = sum(1 for _, ev in payloads if safe_int(ev.get("pdb_to_smiles_mapped_atom_count")) > 0)
    sasa_count = sum(1 for _, ev in payloads if safe_int(ev.get("solvent_exposed_ligand_atom_count")) > 0)
    candidate_count = sum(1 for _, ev in payloads if safe_int(ev.get("candidate_linker_atom_count")) > 0)

    print_flush(f"Unique ligand instances to score: {total:,}")
    print_flush(f"  With SMILES evidence:              {smiles_count:,}")
    print_flush(f"  With PDB-to-SMILES mapping:        {mapping_count:,}")
    print_flush(f"  With solvent-exposed ligand atoms: {sasa_count:,}")
    print_flush(f"  With candidate linker atoms:       {candidate_count:,}")

    start = time.time()
    rows: List[Dict[str, Any]] = []

    if serial or workers <= 1:
        setup_worker(rdkit_warnings)
        for i, payload in enumerate(payloads, start=1):
            rows.append(score_one_ligand(payload))
            if progress_every and (i % progress_every == 0 or i == total):
                print_flush(f"[{i:,}/{total:,}] scored ligand instances ({now_s(start)})")
        return rows

    print_flush(f"Scoring with {workers} worker processes, chunksize={chunksize}...")
    completed = 0
    with ProcessPoolExecutor(
        max_workers=workers,
        initializer=setup_worker,
        initargs=(rdkit_warnings,),
    ) as executor:
        futures = [executor.submit(score_one_ligand, payload) for payload in payloads]
        for fut in as_completed(futures):
            rows.append(fut.result())
            completed += 1
            if progress_every and (completed % progress_every == 0 or completed == total):
                print_flush(f"[{completed:,}/{total:,}] scored ligand instances ({now_s(start)})")

    rows.sort(key=lambda x: (
        norm_upper(x.get("pdb_code")),
        norm_upper(x.get("ligand_resname")),
        norm_str(x.get("ligand_chain")),
        norm_str(x.get("ligand_residue_id")),
    ))
    return rows


# =========================================================
# READINESS SCORING
# =========================================================

def score_target_lysine_accessibility(assessment: Dict[str, Any]) -> float:
    exposed_count = safe_int(assessment.get("exposed_lys_count"))
    lys_count = safe_int(assessment.get("lys_count"))
    exposed_fraction = safe_float(assessment.get("exposed_lys_fraction"))
    lys_surface_fraction = safe_float(assessment.get("lysine_surface_fraction"))

    if lys_count <= 0:
        return 0.0

    score = 0.0
    score += min(40.0, exposed_fraction * 70.0)
    score += min(35.0, exposed_count * 5.0)
    score += min(25.0, lys_surface_fraction * 250.0)
    return round(clamp(score), 2)


def score_ternary_geometry_cue(assessment: Dict[str, Any]) -> Tuple[float, str, int, int, int]:
    near_exposed = safe_int(assessment.get("near_ligand_exposed_lys_count"))
    d = safe_float(assessment.get("min_lys_ligand_distance_a"), default=9999.0)

    short_ok = int(near_exposed > 0 and d <= 8.0)
    medium_ok = int(near_exposed > 0 and d <= 15.0)
    long_ok = int(near_exposed > 0 and d <= 25.0)

    if near_exposed <= 0:
        return 0.0, "No ligand-proximal exposed target lysine cue", 0, 0, 0
    if d <= 8.0:
        return 90.0, "Short-linker geometry cue", short_ok, medium_ok, long_ok
    if d <= 15.0:
        return 70.0, "Medium-linker geometry cue", short_ok, medium_ok, long_ok
    if d <= 25.0:
        return 40.0, "Long-linker exploratory geometry cue", short_ok, medium_ok, long_ok
    return 15.0, "Distant exposed lysine geometry cue", short_ok, medium_ok, long_ok


def readiness_tier(score: float) -> str:
    if score >= 75:
        return "High degrader-design readiness"
    if score >= 55:
        return "Moderate degrader-design readiness"
    if score >= 35:
        return "Exploratory degrader-design readiness"
    return "Weak degrader-design readiness"


def evidence_level(row: Dict[str, Any]) -> str:
    warhead = safe_float(row.get("warhead_linkability_score"))
    lys = safe_float(row.get("target_lysine_accessibility_score"))
    protein = safe_float(row.get("protein_structural_priority_score"))
    ternary = safe_float(row.get("ternary_geometry_cue_score"))

    if warhead >= 75 and lys >= 55 and protein >= 60:
        return "Strong: candidate warhead plus target lysine accessibility"
    if warhead >= 55 and lys >= 45:
        return "Moderate: plausible warhead with accessible lysine signal"
    if warhead >= 55 and lys < 45:
        return "Warhead-supported but lysine accessibility is weak"
    if warhead < 55 and lys >= 55:
        return "Surface-lysine-supported but warhead evidence is weak"
    if ternary >= 70:
        return "Geometry-supported exploratory case"
    return "Limited degrader-readiness evidence"


def index_warheads(warhead_rows: List[Dict[str, Any]]) -> Dict[Tuple[str, str], List[Dict[str, Any]]]:
    index: Dict[Tuple[str, str], List[Dict[str, Any]]] = defaultdict(list)
    for row in warhead_rows:
        pdb = norm_upper(row.get("pdb_code"))
        lig = norm_upper(row.get("ligand_resname"))
        if pdb and lig:
            index[(pdb, lig)].append(row)
    return index


def pick_best_warhead(assessment: Dict[str, Any], warhead_index: Dict[Tuple[str, str], List[Dict[str, Any]]]) -> Optional[Dict[str, Any]]:
    pdb = norm_upper(assessment.get("pdb_code"))
    ligands = split_ligands(assessment.get("candidate_ligand_resnames"))
    candidates: List[Dict[str, Any]] = []

    for lig in ligands:
        candidates.extend(warhead_index.get((pdb, lig), []))

    if not candidates:
        # Fallback to any ligand in this PDB.
        for (wpdb, _lig), rows in warhead_index.items():
            if wpdb == pdb:
                candidates.extend(rows)

    if not candidates:
        return None

    return max(candidates, key=lambda row: (
        safe_float(row.get("warhead_linkability_score")),
        safe_int(row.get("candidate_linker_atom_count")),
        safe_int(row.get("solvent_exposed_mapped_atom_count")),
        safe_int(row.get("functional_group_count")),
    ))


def readiness_flags(assessment: Dict[str, Any], best_warhead: Optional[Dict[str, Any]]) -> str:
    flags: List[str] = []
    if safe_int(assessment.get("lys_count")) <= 0:
        flags.append("no_target_lysines_detected")
    if safe_int(assessment.get("exposed_lys_count")) <= 0:
        flags.append("no_surface_exposed_target_lysines")
    if safe_int(assessment.get("candidate_ligand_count")) <= 0:
        flags.append("no_candidate_ligand_context")
    if safe_int(assessment.get("near_ligand_exposed_lys_count")) <= 0:
        flags.append("no_ligand_proximal_exposed_lysine_geometry_cue")
    if best_warhead is None:
        flags.append("no_warhead_linkability_row")
    else:
        wf = norm_str(best_warhead.get("warhead_flags"))
        if wf:
            flags.extend([f"warhead_{x}" for x in wf.split(";") if x])
    return ";".join(dict.fromkeys(flags))


def build_readiness_rows(assessment_rows: List[Dict[str, str]], warhead_rows: List[Dict[str, Any]], progress_every: int = 5000) -> List[Dict[str, Any]]:
    print_flush("\nBuilding degrader-readiness rows...")
    start = time.time()
    warhead_index = index_warheads(warhead_rows)
    readiness_rows: List[Dict[str, Any]] = []
    total = len(assessment_rows)

    for i, assessment in enumerate(assessment_rows, start=1):
        best = pick_best_warhead(assessment, warhead_index)
        protein_score = safe_float(assessment.get("protacability_proxy_score"))
        warhead_score = safe_float(best.get("warhead_linkability_score")) if best else 0.0
        lys_score = score_target_lysine_accessibility(assessment)
        ternary_score, geometry_class, short_ok, medium_ok, long_ok = score_ternary_geometry_cue(assessment)

        readiness = round(clamp(
            0.40 * warhead_score
            + 0.30 * lys_score
            + 0.20 * protein_score
            + 0.10 * ternary_score
        ), 2)

        row = {
            "virus_name": assessment.get("virus_name", ""),
            "protein_type": assessment.get("protein_type", ""),
            "pdb_code": assessment.get("pdb_code", ""),
            "chain_id": assessment.get("chain_id", ""),
            "model_id": assessment.get("model_id", ""),
            "best_ligand_resname": best.get("ligand_resname", "") if best else "",
            "best_ligand_chain": best.get("ligand_chain", "") if best else "",
            "best_ligand_residue_id": best.get("ligand_residue_id", "") if best else "",
            "protein_structural_priority_score": protein_score,
            "warhead_linkability_score": warhead_score,
            "target_lysine_accessibility_score": lys_score,
            "ternary_geometry_cue_score": ternary_score,
            "degrader_design_readiness_score": readiness,
            "degrader_design_readiness_tier": readiness_tier(readiness),
            "best_linker_geometry_class": geometry_class,
            "short_linker_geometry_feasible": short_ok,
            "medium_linker_geometry_feasible": medium_ok,
            "long_linker_geometry_feasible": long_ok,
            "exposed_lys_count": assessment.get("exposed_lys_count", ""),
            "lys_count": assessment.get("lys_count", ""),
            "exposed_lys_fraction": assessment.get("exposed_lys_fraction", ""),
            "lysine_surface_fraction": assessment.get("lysine_surface_fraction", ""),
            "min_lys_ligand_distance_a": assessment.get("min_lys_ligand_distance_a", ""),
            "near_ligand_exposed_lys_count": assessment.get("near_ligand_exposed_lys_count", ""),
            "candidate_ligand_resnames": assessment.get("candidate_ligand_resnames", ""),
            "readiness_notes": (
                "warhead linkability is ligand-centered; target lysine accessibility is a separate "
                "ubiquitination-readiness cue; ternary geometry is secondary/hypothesis-only"
            ),
        }
        row["evidence_level"] = evidence_level(row)
        row["readiness_flags"] = readiness_flags(assessment, best)
        readiness_rows.append(row)

        if progress_every and (i % progress_every == 0 or i == total):
            print_flush(f"[{i:,}/{total:,}] built readiness rows ({now_s(start)})")

    return readiness_rows


# =========================================================
# SUMMARY / VALIDATION
# =========================================================

def summarize_scores(rows: List[Dict[str, Any]], field: str, label: str) -> None:
    values = [safe_float(r.get(field), default=math.nan) for r in rows]
    values = [v for v in values if not math.isnan(v)]
    if not values:
        print_flush(f"{label}: no values")
        return
    values_sorted = sorted(values)
    n = len(values_sorted)
    mean = sum(values_sorted) / n
    p50 = values_sorted[n // 2]
    p90 = values_sorted[int(n * 0.9)] if n > 1 else values_sorted[0]
    print_flush(f"{label}: n={n:,}, min={values_sorted[0]:.2f}, median={p50:.2f}, mean={mean:.2f}, p90={p90:.2f}, max={values_sorted[-1]:.2f}")


# =========================================================
# CLI / MAIN
# =========================================================

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Fast ligand-centered warhead linkability and degrader-readiness enrichment."
    )
    parser.add_argument("--pdb-root", default=str(PDB_ROOT), help="PDB_FILES root directory. Default: PDB_FILES")
    parser.add_argument("--db", default=str(DEFAULT_DB_PATH), help="SQLite database path. Default: viral_data.db")
    parser.add_argument("--assessment-csv", default=None, help="Override PROTACability_Assessment.csv path")
    parser.add_argument("--ligand-inventory-csv", default=None, help="Override PROTACability_Ligand_Inventory.csv path")
    parser.add_argument("--warhead-output", default=None, help="Override warhead linkability output CSV")
    parser.add_argument("--readiness-output", default=None, help="Override degrader readiness output CSV")
    parser.add_argument("--failures-output", default=None, help="Override failures output CSV")
    parser.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 4) - 2), help="Worker processes. Default: CPU count - 2")
    parser.add_argument("--chunksize", type=int, default=100, help="Reserved for future batched map behavior. Default: 100")
    parser.add_argument("--limit", type=int, default=None, help="Limit unique ligand instances for testing.")
    parser.add_argument("--serial", action="store_true", help="Run ligand scoring serially for debugging.")
    parser.add_argument("--progress-every", type=int, default=500, help="Progress log interval. Default: 500")
    parser.add_argument("--rdkit-warnings", action="store_true", help="Show RDKit parser warnings instead of suppressing them.")
    parser.add_argument("--component-smiles", default="Components-smiles-stereo-oe.smi", help="Optional PDB Chemical Component SMILES .smi file. Default: Components-smiles-stereo-oe.smi")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    run_start = time.time()

    if RDLogger is not None and not args.rdkit_warnings:
        RDLogger.DisableLog("rdApp.*")

    pdb_root = Path(args.pdb_root)
    assessment_csv = Path(args.assessment_csv) if args.assessment_csv else pdb_root / "PROTACability_Assessment.csv"
    ligand_inventory_csv = Path(args.ligand_inventory_csv) if args.ligand_inventory_csv else pdb_root / "PROTACability_Ligand_Inventory.csv"

    warhead_output = Path(args.warhead_output) if args.warhead_output else pdb_root / "PROTACability_Warhead_Linkability.csv"
    readiness_output = Path(args.readiness_output) if args.readiness_output else pdb_root / "PROTACability_Degrader_Readiness.csv"
    failures_output = Path(args.failures_output) if args.failures_output else pdb_root / "PROTACability_Warhead_Linkability_failures.csv"

    if not assessment_csv.exists():
        raise FileNotFoundError(f"Missing assessment CSV: {assessment_csv}")
    if not ligand_inventory_csv.exists():
        raise FileNotFoundError(f"Missing ligand inventory CSV: {ligand_inventory_csv}")

    print_flush("\nPROTACability Warhead Linkability Enrichment FAST")
    print_flush(f"Assessment CSV:       {assessment_csv}")
    print_flush(f"Ligand inventory CSV: {ligand_inventory_csv}")
    print_flush(f"Component SMILES:     {args.component_smiles}")
    print_flush(f"Workers:              {args.workers}")
    print_flush(f"Serial mode:          {args.serial}")
    print_flush(f"Limit:                {args.limit if args.limit is not None else 'none'}")
    print_flush("Scoring note: linkability is ligand/warhead-centered; target lysines are separate accessibility cues.")

    print_flush("\nLoading CSV inputs...")
    assessment_rows = read_csv_rows(assessment_csv)
    ligand_inventory_rows = read_csv_rows(ligand_inventory_csv)
    print_flush(f"Assessment rows:       {len(assessment_rows):,}")
    print_flush(f"Ligand inventory rows: {len(ligand_inventory_rows):,}")

    ligand_tasks, grouped_inventory = deduplicate_ligand_inventory(ligand_inventory_rows, limit=args.limit)
    duplicate_count = len(ligand_inventory_rows) - len(grouped_inventory)
    print_flush("\nDeduplicated ligand inventory:")
    print_flush(f"  Unique ligand instances: {len(ligand_tasks):,}")
    print_flush(f"  Duplicate rows skipped for scoring: {duplicate_count:,}")

    component_path = Path(args.component_smiles) if args.component_smiles else None
    component_smiles = load_component_smiles(component_path)

    db_path = Path(args.db)
    conn = connect_db(db_path)
    if conn is None:
        print_flush(f"\nWARNING: SQLite database not found at {db_path}. Running CSV-only fallback.")
    else:
        print_flush(f"\nUsing SQLite database: {db_path}")

    evidence_maps = build_evidence_maps(conn, ligand_tasks, component_smiles=component_smiles)
    if conn is not None:
        conn.close()

    warhead_rows = build_warhead_rows(
        ligand_tasks=ligand_tasks,
        evidence_maps=evidence_maps,
        workers=args.workers,
        chunksize=args.chunksize,
        serial=args.serial,
        progress_every=args.progress_every,
        rdkit_warnings=args.rdkit_warnings,
    )

    readiness_rows = build_readiness_rows(assessment_rows, warhead_rows, progress_every=max(args.progress_every, 5000))

    failure_rows: List[Dict[str, Any]] = []
    # Add lightweight failures/flags table for missing essentials. Do not duplicate all flags.
    for row in warhead_rows:
        if not safe_int(row.get("smiles_available")):
            failure_rows.append({
                "stage": "warhead_evidence",
                "pdb_code": row.get("pdb_code", ""),
                "ligand_resname": row.get("ligand_resname", ""),
                "reason": "missing_smiles_after_db_and_component_fallback",
            })
        elif not safe_int(row.get("rdkit_valid_smiles")):
            failure_rows.append({
                "stage": "warhead_evidence",
                "pdb_code": row.get("pdb_code", ""),
                "ligand_resname": row.get("ligand_resname", ""),
                "reason": "invalid_smiles_for_rdkit",
            })

    print_flush("\nWriting outputs...")
    write_csv(warhead_output, WARHEAD_FIELDS, warhead_rows)
    write_csv(readiness_output, READINESS_FIELDS, readiness_rows)
    write_csv(failures_output, FAILURE_FIELDS, failure_rows)

    print_flush("\nScore summaries:")
    summarize_scores(warhead_rows, "warhead_linkability_score", "Warhead linkability")
    summarize_scores(readiness_rows, "degrader_design_readiness_score", "Degrader readiness")

    print_flush("\nDone.")
    print_flush(f"Warhead linkability CSV: {warhead_output}")
    print_flush(f"Degrader readiness CSV:  {readiness_output}")
    print_flush(f"Failures/flags CSV:      {failures_output}")
    print_flush(f"Warhead rows written:    {len(warhead_rows):,}")
    print_flush(f"Readiness rows written:  {len(readiness_rows):,}")
    print_flush(f"Failures/flags written:  {len(failure_rows):,}")
    print_flush(f"Total runtime:           {now_s(run_start)}")


if __name__ == "__main__":
    main()
