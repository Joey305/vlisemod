#!/usr/bin/env python3
"""Occurrence-resolved ligand, pocket, target-chain, and lysine geometry.

This stage replaces the legacy PDB_FILES/PROTACability geometry pass with a
CIF-native implementation that reuses the exact ligand occurrences and coherent
conformer policy established in stages 03-09.

Scientific scope
----------------
* Ligand centroids and atom-to-protein distances are computed from the selected
  deposited ligand occurrence.
* Binding-pocket atoms are protein atoms within 5 A (configurable) of any
  selected ligand heavy atom.
* A simple outward-facing cue is recorded for each ligand atom by comparing the
  ligand-centroid->atom vector with the ligand-centroid->pocket-centroid vector.
  Positive ``outward_score`` means the atom points away from the pocket centroid.
  This is a geometric cue, not path finding or proof of linker tolerance.
* Target-side lysine analysis measures whole-residue SASA, heavy-side-chain SASA,
  and terminal NZ-atom SASA. NZ exposure is summarized at >0.1, >1, and >5 A^2.
  Ligand-to-lysine distance is deliberately not used as a ubiquitination metric.
* Atom-level exit vectors are recorded from the ligand centroid through each
  ligand atom, together with a local forward-clearance cue. A clear corridor at
  the configured ray length means clearance reaches that sampling cap; it is not
  an exact distance beyond the cap. These are geometric descriptors for
  modification-site review, not linker-path predictions.

The authoritative source mmCIF is never edited.  A deterministic analysis-only
single-model mmCIF is produced with the same residue-wide conformer policy used
by stage 09.
"""
from __future__ import annotations

import argparse
import importlib
import json
import math
import sqlite3
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from statistics import median

import numpy as np
from Bio.PDB import MMCIFParser
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB.Polypeptide import is_aa
from Bio.SeqUtils.ProtParam import ProteinAnalysis

try:
    from scipy.spatial import cKDTree
except Exception:  # pragma: no cover - fallback used only when scipy is absent
    cKDTree = None

c = importlib.import_module("00_common")
arpeggio = importlib.import_module("09_run_arpeggio")

VERSION = "cif-ligand-geometry-v2.4"
DEFAULT_POCKET_RADIUS = 5.0
DEFAULT_LYSINE_SASA_THRESHOLD = 30.0
DEFAULT_SIDECHAIN_ACCESSIBILITY_THRESHOLD = 0.1
NZ_EXPOSURE_THRESHOLDS = (0.1, 1.0, 5.0)
LYS_BACKBONE_ATOMS = {"N", "CA", "C", "O", "OXT"}
DEFAULT_EXIT_RAY_LENGTH = 8.0
DEFAULT_EXIT_RAY_RADIUS = 2.0
DEFAULT_PROBE_RADIUS = 1.40
DEFAULT_SASA_POINTS = 100
SQLITE_BUSY_TIMEOUT_MS = 60000


def _prepare_sqlite_concurrency(database: str) -> None:
    """Enable one-writer/many-reader operation for the geometry worker pool.

    Stage 10 workers are read-only, while the parent process materializes completed
    results. WAL prevents those normal concurrent reads from being rejected while
    the parent commits. busy_timeout handles brief checkpoint/commit contention.
    This changes only SQLite journaling behavior, not scientific data semantics.
    """
    db = sqlite3.connect(str(Path(database).resolve()), timeout=60.0)
    try:
        db.execute(f"PRAGMA busy_timeout={SQLITE_BUSY_TIMEOUT_MS}")
        db.execute("PRAGMA journal_mode=WAL")
        db.execute("PRAGMA synchronous=NORMAL")
        db.commit()
    finally:
        db.close()


def _is_protein_component(component_id: str) -> bool:
    """Recognize standard and Biopython-known modified amino-acid residues."""
    comp = c.norm(component_id).upper()
    return bool(comp and is_aa(comp, standard=False))

AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}
BASIC_AA = {"K", "R", "H"}
ACIDIC_AA = {"D", "E"}
POLAR_AA = {"S", "T", "N", "Q", "Y", "C", "K", "R", "H", "D", "E"}
HYDROPHOBIC_AA = {"A", "V", "I", "L", "M", "F", "W", "P", "G"}
HYDROGENS = {"H", "D", "T"}


def ensure_schema(database: str) -> None:
    c.create_schema(database)
    with c.dbconn(database) as db:
        db.executescript(
            """
            CREATE TABLE IF NOT EXISTS protacability_ligand_inventory (
                inventory_id INTEGER PRIMARY KEY,
                run_id INTEGER NOT NULL REFERENCES analysis_runs(run_id),
                ligand_instance_id INTEGER NOT NULL REFERENCES ligand_instances(ligand_instance_id),
                virus_name TEXT,
                protein_type TEXT,
                pdb_code TEXT NOT NULL,
                model_id TEXT NOT NULL,
                ligand_resname TEXT NOT NULL,
                ligand_chain TEXT NOT NULL,
                ligand_residue_id TEXT NOT NULL,
                ligand_insertion_code TEXT,
                ligand_atom_count INTEGER NOT NULL,
                ligand_heavy_atom_count INTEGER NOT NULL,
                centroid_x REAL, centroid_y REAL, centroid_z REAL,
                pocket_radius_a REAL NOT NULL,
                pocket_protein_atom_count INTEGER NOT NULL,
                pocket_centroid_x REAL, pocket_centroid_y REAL, pocket_centroid_z REAL,
                min_ligand_protein_distance_a REAL,
                method_version TEXT NOT NULL,
                status TEXT NOT NULL,
                UNIQUE(ligand_instance_id, method_version)
            );
            CREATE TABLE IF NOT EXISTS ligand_atom_geometry (
                geometry_atom_id INTEGER PRIMARY KEY,
                run_id INTEGER NOT NULL REFERENCES analysis_runs(run_id),
                ligand_instance_id INTEGER NOT NULL REFERENCES ligand_instances(ligand_instance_id),
                ligand_instance_atom_id INTEGER NOT NULL REFERENCES ligand_instance_atoms(ligand_instance_atom_id),
                atom_site_id TEXT,
                exact_atom TEXT,
                element TEXT,
                x REAL, y REAL, z REAL,
                radial_distance_from_ligand_centroid_a REAL,
                nearest_protein_distance_a REAL,
                outward_cosine REAL,
                outward_score REAL,
                points_away_from_pocket INTEGER,
                pocket_centroid_available INTEGER NOT NULL,
                method_version TEXT NOT NULL,
                UNIQUE(ligand_instance_atom_id, method_version)
            );
            CREATE TABLE IF NOT EXISTS target_chain_geometry (
                chain_geometry_id INTEGER PRIMARY KEY,
                run_id INTEGER NOT NULL REFERENCES analysis_runs(run_id),
                ligand_instance_id INTEGER NOT NULL REFERENCES ligand_instances(ligand_instance_id),
                virus_name TEXT,
                protein_type TEXT,
                pdb_code TEXT NOT NULL,
                chain_id TEXT NOT NULL,
                model_id TEXT NOT NULL,
                chain_length_aa INTEGER NOT NULL,
                lys_count INTEGER NOT NULL,
                exposed_lys_count INTEGER NOT NULL,
                exposed_lys_fraction REAL NOT NULL,
                near_ligand_lys_count INTEGER NOT NULL,
                near_ligand_exposed_lys_count INTEGER NOT NULL,
                min_lys_ligand_distance_a REAL,
                median_lys_ligand_distance_a REAL,
                total_sasa_a2 REAL NOT NULL,
                lysine_sasa_a2 REAL NOT NULL,
                lysine_surface_fraction REAL,
                isoelectric_point REAL,
                basic_fraction REAL NOT NULL,
                acidic_fraction REAL NOT NULL,
                polar_fraction REAL NOT NULL,
                hydrophobic_fraction REAL NOT NULL,
                has_exposed_lysine INTEGER NOT NULL,
                has_ligand_proximal_exposed_lysine INTEGER NOT NULL,
                lysine_sidechain_sasa_a2 REAL NOT NULL DEFAULT 0,
                lysine_nz_sasa_a2 REAL NOT NULL DEFAULT 0,
                nz_observed_lys_count INTEGER NOT NULL DEFAULT 0,
                nz_observed_lys_fraction REAL NOT NULL DEFAULT 0,
                sidechain_accessible_lys_count INTEGER NOT NULL DEFAULT 0,
                sidechain_accessible_lys_fraction REAL NOT NULL DEFAULT 0,
                nz_exposed_lys_count_gt_0_1 INTEGER NOT NULL DEFAULT 0,
                nz_exposed_lys_fraction_gt_0_1 REAL NOT NULL DEFAULT 0,
                nz_exposed_lys_count_gt_1 INTEGER NOT NULL DEFAULT 0,
                nz_exposed_lys_fraction_gt_1 REAL NOT NULL DEFAULT 0,
                nz_exposed_lys_count_gt_5 INTEGER NOT NULL DEFAULT 0,
                nz_exposed_lys_fraction_gt_5 REAL NOT NULL DEFAULT 0,
                method_version TEXT NOT NULL,
                UNIQUE(ligand_instance_id, chain_id, method_version)
            );
            CREATE TABLE IF NOT EXISTS target_surface_lysines (
                surface_lysine_id INTEGER PRIMARY KEY,
                run_id INTEGER NOT NULL REFERENCES analysis_runs(run_id),
                ligand_instance_id INTEGER NOT NULL REFERENCES ligand_instances(ligand_instance_id),
                virus_name TEXT,
                protein_type TEXT,
                pdb_code TEXT NOT NULL,
                chain_id TEXT NOT NULL,
                model_id TEXT NOT NULL,
                lys_residue_id TEXT NOT NULL,
                lys_insertion_code TEXT,
                lys_observed_index INTEGER NOT NULL,
                lysine_sasa_a2 REAL NOT NULL,
                is_surface_exposed INTEGER NOT NULL,
                lysine_sidechain_sasa_a2 REAL NOT NULL DEFAULT 0,
                nz_sasa_a2 REAL,
                nz_observed INTEGER NOT NULL DEFAULT 0,
                sidechain_accessible INTEGER NOT NULL DEFAULT 0,
                nz_exposed_gt_0_1 INTEGER NOT NULL DEFAULT 0,
                nz_exposed_gt_1 INTEGER NOT NULL DEFAULT 0,
                nz_exposed_gt_5 INTEGER NOT NULL DEFAULT 0,
                reference_atom TEXT,
                reference_atom_x REAL, reference_atom_y REAL, reference_atom_z REAL,
                method_version TEXT NOT NULL,
                UNIQUE(ligand_instance_id, chain_id, lys_residue_id, lys_insertion_code, method_version)
            );
            CREATE TABLE IF NOT EXISTS protacability_lysine_proximity (
                lysine_geometry_id INTEGER PRIMARY KEY,
                run_id INTEGER NOT NULL REFERENCES analysis_runs(run_id),
                ligand_instance_id INTEGER NOT NULL REFERENCES ligand_instances(ligand_instance_id),
                virus_name TEXT,
                protein_type TEXT,
                pdb_code TEXT NOT NULL,
                chain_id TEXT NOT NULL,
                model_id TEXT NOT NULL,
                lys_residue_id TEXT NOT NULL,
                lys_insertion_code TEXT,
                lys_observed_index INTEGER NOT NULL,
                lysine_sasa_a2 REAL NOT NULL,
                is_surface_exposed INTEGER NOT NULL,
                distance_atom TEXT,
                distance_atom_x REAL, distance_atom_y REAL, distance_atom_z REAL,
                nearest_ligand_resname TEXT,
                nearest_ligand_chain TEXT,
                nearest_ligand_residue_id TEXT,
                nearest_ligand_insertion_code TEXT,
                nearest_ligand_atom TEXT,
                nearest_ligand_instance_atom_id INTEGER REFERENCES ligand_instance_atoms(ligand_instance_atom_id),
                nearest_ligand_distance_a REAL,
                is_ligand_proximal INTEGER NOT NULL,
                linker_site_class TEXT NOT NULL,
                method_version TEXT NOT NULL,
                UNIQUE(ligand_instance_id, chain_id, lys_residue_id, lys_insertion_code, method_version)
            );
            CREATE TABLE IF NOT EXISTS ligand_binding_pocket_atoms (
                ligand_pocket_atom_id INTEGER PRIMARY KEY,
                run_id INTEGER NOT NULL REFERENCES analysis_runs(run_id),
                ligand_instance_id INTEGER NOT NULL REFERENCES ligand_instances(ligand_instance_id),
                source_atom_site_id TEXT,
                label_asym_id TEXT,
                auth_asym_id TEXT,
                label_seq_id TEXT,
                auth_seq_id TEXT,
                insertion_code TEXT,
                residue_name TEXT,
                label_atom_id TEXT,
                auth_atom_id TEXT,
                element TEXT,
                x REAL, y REAL, z REAL,
                distance_a REAL NOT NULL,
                pocket_radius_a REAL NOT NULL,
                method_version TEXT NOT NULL,
                UNIQUE(ligand_instance_id, source_atom_site_id, method_version)
            );
            """
        )
        # Additive migrations.  v2.3 keeps the old exit_vector_clear name for
        # compatibility but introduces the clearer local_corridor_clear name and
        # explicit capped-clearance provenance.
        migrations = {
            "ligand_atom_geometry": (
                ("exit_vector_x", "REAL"), ("exit_vector_y", "REAL"), ("exit_vector_z", "REAL"),
                ("forward_clearance_a", "REAL"), ("forward_obstruction_count", "INTEGER"),
                ("exit_vector_clear", "INTEGER"), ("local_corridor_clear", "INTEGER"),
                ("forward_clearance_reaches_cap", "INTEGER"),
                ("exit_ray_length_a", "REAL"), ("exit_ray_radius_a", "REAL"),
            ),
            "target_chain_geometry": (
                ("lysine_sidechain_sasa_a2", "REAL NOT NULL DEFAULT 0"),
                ("lysine_nz_sasa_a2", "REAL NOT NULL DEFAULT 0"),
                ("nz_observed_lys_count", "INTEGER NOT NULL DEFAULT 0"),
                ("nz_observed_lys_fraction", "REAL NOT NULL DEFAULT 0"),
                ("sidechain_accessible_lys_count", "INTEGER NOT NULL DEFAULT 0"),
                ("sidechain_accessible_lys_fraction", "REAL NOT NULL DEFAULT 0"),
                ("nz_exposed_lys_count_gt_0_1", "INTEGER NOT NULL DEFAULT 0"),
                ("nz_exposed_lys_fraction_gt_0_1", "REAL NOT NULL DEFAULT 0"),
                ("nz_exposed_lys_count_gt_1", "INTEGER NOT NULL DEFAULT 0"),
                ("nz_exposed_lys_fraction_gt_1", "REAL NOT NULL DEFAULT 0"),
                ("nz_exposed_lys_count_gt_5", "INTEGER NOT NULL DEFAULT 0"),
                ("nz_exposed_lys_fraction_gt_5", "REAL NOT NULL DEFAULT 0"),
            ),
            "target_surface_lysines": (
                ("lysine_sidechain_sasa_a2", "REAL NOT NULL DEFAULT 0"),
                ("nz_sasa_a2", "REAL"),
                ("nz_observed", "INTEGER NOT NULL DEFAULT 0"),
                ("sidechain_accessible", "INTEGER NOT NULL DEFAULT 0"),
                ("nz_exposed_gt_0_1", "INTEGER NOT NULL DEFAULT 0"),
                ("nz_exposed_gt_1", "INTEGER NOT NULL DEFAULT 0"),
                ("nz_exposed_gt_5", "INTEGER NOT NULL DEFAULT 0"),
            ),
        }
        for table, columns in migrations.items():
            existing = {r[1] for r in db.execute(f"PRAGMA table_info({table})")}
            for column, kind in columns:
                if column not in existing:
                    db.execute(f"ALTER TABLE {table} ADD COLUMN {column} {kind}")


def _num_sort(value):
    text = str(value or "")
    try:
        return (0, float(text))
    except Exception:
        return (1, text)


def _safe_slug(value: str) -> str:
    return "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in str(value))


def _distance_matrix_min(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    if not len(a) or not len(b):
        return np.full(len(a), np.nan)
    if cKDTree is not None:
        tree = cKDTree(b)
        return np.asarray(tree.query(a, k=1)[0], dtype=float)
    # Bounded fallback; ligand atom sets are small.
    return np.sqrt(((a[:, None, :] - b[None, :, :]) ** 2).sum(axis=2)).min(axis=1)


def _protein_neighbors(protein_coords: np.ndarray, ligand_coords: np.ndarray, radius: float):
    if not len(protein_coords) or not len(ligand_coords):
        return np.array([], dtype=int)
    if cKDTree is not None:
        tree = cKDTree(protein_coords)
        hits = set()
        for found in tree.query_ball_point(ligand_coords, r=radius):
            hits.update(found)
        return np.asarray(sorted(hits), dtype=int)
    d2 = ((protein_coords[:, None, :] - ligand_coords[None, :, :]) ** 2).sum(axis=2)
    return np.where((d2 <= radius * radius).any(axis=1))[0]


def _exit_ray_metrics(atom_coord: np.ndarray, unit_vector: np.ndarray | None,
                      protein_coords: np.ndarray, ray_length: float, ray_radius: float):
    """Return local forward-clearance descriptors for one candidate exit vector.

    Protein atoms are counted only when they lie in front of the ligand atom,
    within ``ray_length`` along the vector and within ``ray_radius`` of the ray.
    This is a local steric descriptor, not a linker-routing or accessibility proof.
    """
    if unit_vector is None or not len(protein_coords):
        return None, None, None, None
    rel = protein_coords - atom_coord
    forward = rel @ unit_vector
    rel_sq = np.einsum("ij,ij->i", rel, rel)
    lateral_sq = np.maximum(0.0, rel_sq - forward * forward)
    mask = (forward > 1e-6) & (forward <= ray_length) & (lateral_sq <= ray_radius * ray_radius)
    count = int(np.count_nonzero(mask))
    reaches_cap = int(count == 0)
    # When reaches_cap=1, clearance means "at least ray_length", not an
    # exact steric distance beyond the sampled corridor.
    clearance = float(np.min(forward[mask])) if count else float(ray_length)
    return clearance, count, reaches_cap, reaches_cap


def _chain_labels(db, structure_id: int):
    rows = db.execute(
        "SELECT virus_label,protein_label FROM structure_classifications WHERE structure_id=?",
        (structure_id,),
    ).fetchall()
    viruses = sorted({c.norm(r["virus_label"]) for r in rows if c.norm(r["virus_label"])})
    proteins = sorted({c.norm(r["protein_label"]) for r in rows if c.norm(r["protein_label"])})
    return ";".join(viruses), ";".join(proteins)


def _parse_context(derived_path: Path, atom_map, probe_radius: float, sasa_points: int):
    parser = MMCIFParser(QUIET=True, auth_chains=True, auth_residues=True)
    structure = parser.get_structure("geometry_context", str(derived_path))
    sr = ShrakeRupley(probe_radius=probe_radius, n_points=sasa_points)
    sr.compute(structure, level="A")

    by_serial = {str(m["derived_atom_id"]): m for m in atom_map}
    records = []
    for atom in structure.get_atoms():
        serial = str(atom.get_serial_number())
        m = by_serial.get(serial)
        if not m:
            continue
        coord = atom.get_coord()
        rec = dict(m)
        rec["coord"] = np.asarray([float(coord[0]), float(coord[1]), float(coord[2])], dtype=float)
        rec["sasa"] = float(getattr(atom, "sasa", 0.0) or 0.0)
        records.append(rec)
    return records


def _build_model_context(db, structure_id: int, model_id: str, representative_iid: int,
                         output_path: Path, probe_radius: float, sasa_points: int):
    row = db.execute(
        """SELECT i.*,s.entry_id,s.source_cif_path,s.source_cif_sha256
           FROM ligand_instances i JOIN structures s ON s.structure_id=i.structure_id
           WHERE i.ligand_instance_id=?""",
        (representative_iid,),
    ).fetchone()
    if row is None:
        raise ValueError(f"missing representative ligand instance {representative_iid}")
    derived_path, atom_map, operations, manifest = arpeggio.build_derived_input(
        db, row, output_path, strategy="sanitized_full", context_radius=12.0
    )
    records = _parse_context(derived_path, atom_map, probe_radius, sasa_points)

    # Use Biopython's amino-acid recognition for pocket geometry so common
    # modified amino acids (for example MSE) remain part of the protein surface.
    # A model containing no recognized amino-acid polymer is not an execution
    # failure: target-side protein/lysine geometry is simply not applicable.
    protein_atoms = [r for r in records if _is_protein_component(r["source_component_id"])]
    target_context_status = "complete" if protein_atoms else "not_applicable_no_protein_atoms"

    residues = defaultdict(list)
    for rec in protein_atoms:
        key = (
            c.norm(rec["source_auth_asym_id"]), c.norm(rec["source_auth_seq_id"]),
            c.norm(rec["source_insertion_code"]), c.norm(rec["source_component_id"]).upper(),
            c.norm(rec["source_label_seq_id"]),
        )
        residues[key].append(rec)

    chains = defaultdict(list)
    for key, atoms in residues.items():
        chains[key[0]].append((key, atoms))

    chain_base = {}
    lysines = []
    for chain_id, items in chains.items():
        items = sorted(items, key=lambda x: (_num_sort(x[0][4]), _num_sort(x[0][1]), x[0][2]))
        # Keep the geometric chain length inclusive of recognized modified amino
        # acids. Sequence-derived descriptors use only standard residues, so a
        # modified residue is never silently rewritten to a different chemistry.
        sequence = "".join(AA3_TO_1[item[0][3]] for item in items if item[0][3] in AA3_TO_1)
        aa_count = len(items)
        total_sasa = sum(a["sasa"] for _, atoms in items for a in atoms)
        lys_items = [(k, atoms) for k, atoms in items if k[3] == "LYS"]
        chain_lysines = []
        for observed_index, (key, atoms) in enumerate(lys_items, 1):
            by_name = {c.norm(a["source_auth_atom_id"] or a["source_label_atom_id"]): a for a in atoms}
            nz = by_name.get("NZ")
            ref = nz or by_name.get("CA")
            if ref is None:
                ref = next((a for a in atoms if c.norm(a["source_element"]).upper() not in HYDROGENS), atoms[0])
            residue_sasa = sum(a["sasa"] for a in atoms)
            sidechain_atoms = [
                a for a in atoms
                if c.norm(a["source_element"]).upper() not in HYDROGENS
                and c.norm(a["source_auth_atom_id"] or a["source_label_atom_id"]).upper() not in LYS_BACKBONE_ATOMS
            ]
            sidechain_sasa = sum(a["sasa"] for a in sidechain_atoms)
            nz_sasa = float(nz["sasa"]) if nz is not None else None
            metric = {
                "chain_id": chain_id,
                "auth_seq_id": key[1],
                "insertion_code": key[2],
                "observed_index": observed_index,
                "sasa": residue_sasa,
                "sidechain_sasa": sidechain_sasa,
                "nz_sasa": nz_sasa,
                "nz_observed": int(nz is not None),
                "sidechain_accessible": int(sidechain_sasa > DEFAULT_SIDECHAIN_ACCESSIBILITY_THRESHOLD),
                "nz_exposed_gt_0_1": int(nz_sasa is not None and nz_sasa > 0.1),
                "nz_exposed_gt_1": int(nz_sasa is not None and nz_sasa > 1.0),
                "nz_exposed_gt_5": int(nz_sasa is not None and nz_sasa > 5.0),
                "distance_atom": c.norm(ref["source_auth_atom_id"] or ref["source_label_atom_id"]),
                "coord": ref["coord"],
            }
            chain_lysines.append(metric)
            lysines.append(metric)

        lys_count = len(chain_lysines)
        lysine_sasa = sum(x["sasa"] for x in chain_lysines)
        sidechain_sasa_total = sum(x["sidechain_sasa"] for x in chain_lysines)
        nz_sasa_total = sum((x["nz_sasa"] or 0.0) for x in chain_lysines)
        nz_observed = sum(x["nz_observed"] for x in chain_lysines)
        sidechain_accessible = sum(x["sidechain_accessible"] for x in chain_lysines)
        nz_0_1 = sum(x["nz_exposed_gt_0_1"] for x in chain_lysines)
        nz_1 = sum(x["nz_exposed_gt_1"] for x in chain_lysines)
        nz_5 = sum(x["nz_exposed_gt_5"] for x in chain_lysines)
        try:
            pi = float(ProteinAnalysis(sequence).isoelectric_point()) if len(sequence) >= 2 else None
        except Exception:
            pi = None
        counts = {name: sum(aa in pool for aa in sequence) for name, pool in (
            ("basic", BASIC_AA), ("acidic", ACIDIC_AA), ("polar", POLAR_AA), ("hydrophobic", HYDROPHOBIC_AA)
        )}
        chain_base[chain_id] = {
            "chain_length_aa": aa_count,
            "lys_count": lys_count,
            "total_sasa_a2": total_sasa,
            "lysine_sasa_a2": lysine_sasa,
            "lysine_surface_fraction": lysine_sasa / total_sasa if total_sasa > 0 else None,
            "lysine_sidechain_sasa_a2": sidechain_sasa_total,
            "lysine_nz_sasa_a2": nz_sasa_total,
            "nz_observed_lys_count": nz_observed,
            "nz_observed_lys_fraction": nz_observed / lys_count if lys_count else 0.0,
            "sidechain_accessible_lys_count": sidechain_accessible,
            "sidechain_accessible_lys_fraction": sidechain_accessible / lys_count if lys_count else 0.0,
            "nz_exposed_lys_count_gt_0_1": nz_0_1,
            "nz_exposed_lys_fraction_gt_0_1": nz_0_1 / lys_count if lys_count else 0.0,
            "nz_exposed_lys_count_gt_1": nz_1,
            "nz_exposed_lys_fraction_gt_1": nz_1 / lys_count if lys_count else 0.0,
            "nz_exposed_lys_count_gt_5": nz_5,
            "nz_exposed_lys_fraction_gt_5": nz_5 / lys_count if lys_count else 0.0,
            "isoelectric_point": pi,
            "basic_fraction": counts["basic"] / aa_count if aa_count else 0.0,
            "acidic_fraction": counts["acidic"] / aa_count if aa_count else 0.0,
            "polar_fraction": counts["polar"] / aa_count if aa_count else 0.0,
            "hydrophobic_fraction": counts["hydrophobic"] / aa_count if aa_count else 0.0,
        }

    pcoords = np.vstack([r["coord"] for r in protein_atoms]) if protein_atoms else np.empty((0, 3), dtype=float)
    return {
        "protein_atoms": protein_atoms,
        "protein_coords": pcoords,
        "chain_base": chain_base,
        "lysines": lysines,
        "target_context_status": target_context_status,
        "derived_path": str(derived_path),
        "manifest": manifest,
    }


def _instance_atoms(db, iid: int):
    return db.execute(
        """SELECT ligand_instance_atom_id,atom_site_id,label_atom_id,auth_atom_id,element,x,y,z
           FROM ligand_instance_atoms
           WHERE ligand_instance_id=? AND selected_conformer=1
           ORDER BY ligand_instance_atom_id""",
        (iid,),
    ).fetchall()


def _analyze_group(payload):
    (database, run_id, structure_id, model_id, iids, representative_iid, context_path,
     pocket_radius, lysine_sasa_threshold, exit_ray_length, exit_ray_radius,
     probe_radius, sasa_points) = payload
    db = sqlite3.connect(
        f"file:{Path(database).resolve()}?mode=ro", uri=True, timeout=60.0
    )
    db.row_factory = sqlite3.Row
    db.execute(f"PRAGMA busy_timeout={SQLITE_BUSY_TIMEOUT_MS}")
    db.execute("PRAGMA query_only=ON")
    try:
        model = _build_model_context(
            db, structure_id, model_id, representative_iid, Path(context_path), probe_radius, sasa_points
        )
        virus_name, protein_type = _chain_labels(db, structure_id)
        structure = db.execute("SELECT entry_id FROM structures WHERE structure_id=?", (structure_id,)).fetchone()
        pdb_code = structure["entry_id"] if structure else ""
        protein_atoms = model["protein_atoms"]
        protein_coords = model["protein_coords"]
        chain_base = model["chain_base"]
        lysines = model["lysines"]

        inventory_rows, atom_rows, chain_rows, lysine_rows, pocket_rows = [], [], [], [], []
        for iid in iids:
            inst = db.execute(
                """SELECT i.*,l.component_id FROM ligand_instances i JOIN ligands l ON l.ligand_id=i.ligand_id
                   WHERE i.ligand_instance_id=?""", (iid,)
            ).fetchone()
            atoms = _instance_atoms(db, iid)
            if not atoms:
                raise ValueError(f"ligand instance {iid} has no selected conformer atoms")
            heavy = [a for a in atoms if c.norm(a["element"]).upper() not in HYDROGENS]
            basis = heavy or list(atoms)
            ligand_coords = np.asarray([[a["x"], a["y"], a["z"]] for a in basis], dtype=float)
            centroid = ligand_coords.mean(axis=0)
            pocket_idx = _protein_neighbors(protein_coords, ligand_coords, pocket_radius)
            pocket_coords = protein_coords[pocket_idx] if len(pocket_idx) else np.empty((0, 3))
            pocket_centroid = pocket_coords.mean(axis=0) if len(pocket_coords) else None
            min_lp = float(np.nanmin(_distance_matrix_min(ligand_coords, protein_coords))) if len(protein_coords) else None

            inventory_rows.append((
                run_id, iid, virus_name, protein_type, pdb_code, str(model_id), inst["label_comp_id"],
                inst["auth_asym_id"], inst["auth_seq_id"], inst["insertion_code_normalized"],
                len(atoms), len(heavy), float(centroid[0]), float(centroid[1]), float(centroid[2]),
                pocket_radius, int(len(pocket_idx)),
                float(pocket_centroid[0]) if pocket_centroid is not None else None,
                float(pocket_centroid[1]) if pocket_centroid is not None else None,
                float(pocket_centroid[2]) if pocket_centroid is not None else None,
                min_lp, VERSION, model["target_context_status"]
            ))

            all_coords = np.asarray([[a["x"], a["y"], a["z"]] for a in atoms], dtype=float)
            nearest = _distance_matrix_min(all_coords, protein_coords)
            inward = (pocket_centroid - centroid) if pocket_centroid is not None else None
            inward_norm = float(np.linalg.norm(inward)) if inward is not None else 0.0
            for idx, atom in enumerate(atoms):
                vec = all_coords[idx] - centroid
                vec_norm = float(np.linalg.norm(vec))
                unit = (vec / vec_norm) if vec_norm > 1e-12 else None
                cosine = None
                outward_score = None
                points_away = None
                if pocket_centroid is not None and inward_norm > 1e-12 and unit is not None:
                    cosine = float(np.dot(unit, inward / inward_norm))
                    cosine = max(-1.0, min(1.0, cosine))
                    outward_score = -cosine
                    points_away = int(outward_score > 0.0)
                clearance, obstruction_count, local_corridor_clear, reaches_cap = _exit_ray_metrics(
                    all_coords[idx], unit, protein_coords, exit_ray_length, exit_ray_radius
                )
                atom_rows.append((
                    run_id, iid, atom["ligand_instance_atom_id"], atom["atom_site_id"],
                    atom["auth_atom_id"] or atom["label_atom_id"], atom["element"],
                    atom["x"], atom["y"], atom["z"], vec_norm,
                    float(nearest[idx]) if math.isfinite(float(nearest[idx])) else None,
                    cosine, outward_score, points_away, int(pocket_centroid is not None),
                    float(unit[0]) if unit is not None else None,
                    float(unit[1]) if unit is not None else None,
                    float(unit[2]) if unit is not None else None,
                    clearance, obstruction_count, local_corridor_clear, local_corridor_clear, reaches_cap,
                    exit_ray_length, exit_ray_radius, VERSION
                ))

            for pi in pocket_idx:
                patom = protein_atoms[int(pi)]
                d = float(_distance_matrix_min(np.asarray([patom["coord"]]), ligand_coords)[0])
                pocket_rows.append((
                    run_id, iid, patom["source_atom_site_id"], patom["source_label_asym_id"],
                    patom["source_auth_asym_id"], patom["source_label_seq_id"], patom["source_auth_seq_id"],
                    patom["source_insertion_code"], patom["source_component_id"], patom["source_label_atom_id"],
                    patom["source_auth_atom_id"], patom["source_element"],
                    float(patom["coord"][0]), float(patom["coord"][1]), float(patom["coord"][2]),
                    d, pocket_radius, VERSION
                ))

            per_chain_exposed = defaultdict(int)
            for lys in lysines:
                exposed = int(lys["sasa"] >= lysine_sasa_threshold)
                per_chain_exposed[lys["chain_id"]] += exposed
                lysine_rows.append((
                    run_id, iid, virus_name, protein_type, pdb_code, lys["chain_id"], str(model_id),
                    lys["auth_seq_id"], lys["insertion_code"], lys["observed_index"], lys["sasa"], exposed,
                    lys["sidechain_sasa"], lys["nz_sasa"], lys["nz_observed"], lys["sidechain_accessible"],
                    lys["nz_exposed_gt_0_1"], lys["nz_exposed_gt_1"], lys["nz_exposed_gt_5"],
                    lys["distance_atom"], float(lys["coord"][0]), float(lys["coord"][1]), float(lys["coord"][2]), VERSION
                ))

            for chain_id, base in chain_base.items():
                lys_count = base["lys_count"]
                exposed_count = per_chain_exposed.get(chain_id, 0)
                chain_rows.append((
                    run_id, iid, virus_name, protein_type, pdb_code, chain_id, str(model_id),
                    base["chain_length_aa"], lys_count, exposed_count,
                    exposed_count / lys_count if lys_count else 0.0,
                    0, 0, None, None,
                    base["total_sasa_a2"], base["lysine_sasa_a2"], base["lysine_surface_fraction"],
                    base["isoelectric_point"], base["basic_fraction"], base["acidic_fraction"],
                    base["polar_fraction"], base["hydrophobic_fraction"], int(exposed_count > 0),
                    0, base["lysine_sidechain_sasa_a2"], base["lysine_nz_sasa_a2"],
                    base["nz_observed_lys_count"], base["nz_observed_lys_fraction"],
                    base["sidechain_accessible_lys_count"], base["sidechain_accessible_lys_fraction"],
                    base["nz_exposed_lys_count_gt_0_1"], base["nz_exposed_lys_fraction_gt_0_1"],
                    base["nz_exposed_lys_count_gt_1"], base["nz_exposed_lys_fraction_gt_1"],
                    base["nz_exposed_lys_count_gt_5"], base["nz_exposed_lys_fraction_gt_5"], VERSION
                ))
        return {
            "status": "completed", "structure_id": structure_id, "model_id": model_id, "iids": iids,
            "target_context_status": model["target_context_status"],
            "inventory": inventory_rows, "atoms": atom_rows, "chains": chain_rows,
            "surface_lysines": lysine_rows, "pocket": pocket_rows,
        }
    except Exception as exc:
        return {
            "status": "failed", "structure_id": structure_id, "model_id": model_id, "iids": iids,
            "error": f"{type(exc).__name__}: {exc}",
        }
    finally:
        db.close()


def _write_outputs(database: str):
    c.dirs()
    with c.dbconn(database) as db:
        counts = {
            "ligand_inventory": db.execute("SELECT count(*) FROM protacability_ligand_inventory WHERE method_version=?", (VERSION,)).fetchone()[0],
            "protein_applicable_instances": db.execute("SELECT count(*) FROM protacability_ligand_inventory WHERE method_version=? AND status='complete'", (VERSION,)).fetchone()[0],
            "no_protein_target_instances": db.execute("SELECT count(*) FROM protacability_ligand_inventory WHERE method_version=? AND status='not_applicable_no_protein_atoms'", (VERSION,)).fetchone()[0],
            "ligand_atom_geometry": db.execute("SELECT count(*) FROM ligand_atom_geometry WHERE method_version=?", (VERSION,)).fetchone()[0],
            "target_chain_geometry": db.execute("SELECT count(*) FROM target_chain_geometry WHERE method_version=?", (VERSION,)).fetchone()[0],
            "surface_lysines": db.execute("SELECT count(*) FROM target_surface_lysines WHERE method_version=?", (VERSION,)).fetchone()[0],
            "binding_pocket_atoms": db.execute("SELECT count(*) FROM ligand_binding_pocket_atoms WHERE method_version=?", (VERSION,)).fetchone()[0],
        }
        lys_summary = db.execute(
            """SELECT COALESCE(sum(lys_count),0),COALESCE(sum(nz_observed_lys_count),0),
                      COALESCE(sum(nz_exposed_lys_count_gt_0_1),0),COALESCE(sum(nz_exposed_lys_count_gt_1),0),
                      COALESCE(sum(nz_exposed_lys_count_gt_5),0)
               FROM target_chain_geometry WHERE method_version=?""", (VERSION,)
        ).fetchone()
    report = ["# Stage 10 geometry report", "", f"* Method: {VERSION}"] + [f"* {k.replace('_',' ')}: {v}" for k, v in counts.items()]
    report += ["", "## Target lysine terminal-amine accessibility",
               f"* Lysine residues represented across ligand/chain contexts: {lys_summary[0]}",
               f"* NZ atoms observed: {lys_summary[1]}",
               f"* NZ SASA >0.1 A^2: {lys_summary[2]}",
               f"* NZ SASA >1 A^2: {lys_summary[3]}",
               f"* NZ SASA >5 A^2: {lys_summary[4]}",
               "", "Forward clearance is sampled only to the configured ray length. A clear row at 8 A means >=8 A within the tested local corridor, not an exact clearance beyond 8 A."]
    (c.ROOT / "outputs" / "GEOMETRY_STAGE_REPORT.md").write_text("\n".join(report) + "\n")


def run(database: str, limit=None, pdb_id=None, instance_id=None, workers=4, resume=False,
        pocket_radius=DEFAULT_POCKET_RADIUS, lysine_sasa_threshold=DEFAULT_LYSINE_SASA_THRESHOLD,
        exit_ray_length=DEFAULT_EXIT_RAY_LENGTH, exit_ray_radius=DEFAULT_EXIT_RAY_RADIUS,
        probe_radius=DEFAULT_PROBE_RADIUS, sasa_points=DEFAULT_SASA_POINTS, progress_every=25):
    ensure_schema(database); c.dirs(); _prepare_sqlite_concurrency(database)
    with c.dbconn(database) as db:
        db.execute(f"PRAGMA busy_timeout={SQLITE_BUSY_TIMEOUT_MS}")
        db.execute("UPDATE analysis_runs SET status='interrupted',end_time=? WHERE stage='geometry' AND status='running'", (c.now(),))
        q = """SELECT i.ligand_instance_id,i.structure_id,i.deposited_model_num
               FROM ligand_instances i JOIN structures s ON s.structure_id=i.structure_id
               WHERE i.curation_status='included'"""
        args = []
        if pdb_id:
            q += " AND UPPER(s.entry_id)=UPPER(?)"; args.append(pdb_id)
        if instance_id:
            q += " AND i.ligand_instance_id=?"; args.append(instance_id)
        if resume:
            q += " AND NOT EXISTS (SELECT 1 FROM protacability_ligand_inventory g WHERE g.ligand_instance_id=i.ligand_instance_id AND g.method_version=? AND g.status IN ('complete','not_applicable_no_protein_atoms'))"
            args.append(VERSION)
        rows = db.execute(q + " ORDER BY i.structure_id,i.deposited_model_num,i.ligand_instance_id", args).fetchall()
        if limit:
            rows = rows[:limit]
        ids = [r["ligand_instance_id"] for r in rows]
        rid = c.run_start(db, "geometry", {
            "method": VERSION, "limit": limit, "pdb_id": pdb_id, "ligand_instance_id": instance_id,
            "workers": workers, "resume": resume, "pocket_radius_a": pocket_radius,
            "lysine_sasa_threshold_a2": lysine_sasa_threshold, "exit_ray_length_a": exit_ray_length,
            "exit_ray_radius_a": exit_ray_radius, "ligand_to_lysine_distance_used": False,
            "probe_radius_a": probe_radius, "sasa_points": sasa_points,
        })
        if ids and not resume:
            marks = ",".join("?" for _ in ids)
            for table in ("protacability_ligand_inventory", "ligand_atom_geometry", "target_chain_geometry",
                          "target_surface_lysines", "ligand_binding_pocket_atoms"):
                db.execute(f"DELETE FROM {table} WHERE method_version=? AND ligand_instance_id IN ({marks})", [VERSION, *ids])
            db.execute(f"DELETE FROM receptor_binding_pocket_atoms WHERE ligand_instance_id IN ({marks})", ids)
        db.commit()

    groups = defaultdict(list)
    for r in rows:
        groups[(r["structure_id"], str(r["deposited_model_num"]))].append(r["ligand_instance_id"])
    payloads = []
    for (sid, model_id), iids in groups.items():
        context_dir = c.ROOT / "outputs" / "geometry" / "contexts" / str(rid) / f"{sid}_{_safe_slug(model_id)}"
        context_dir.mkdir(parents=True, exist_ok=True)
        payloads.append((
            str(database), rid, sid, model_id, iids, iids[0], str(context_dir / "derived_input.cif"),
            pocket_radius, lysine_sasa_threshold, exit_ray_length, exit_ray_radius, probe_radius, sasa_points
        ))

    print(f"geometry selection: pdb_id={pdb_id or 'ALL'} instances={len(ids)} model_groups={len(payloads)} method={VERSION}", flush=True)
    processed = success = not_applicable = failures = 0
    insert_sql = {
        "inventory": """INSERT OR REPLACE INTO protacability_ligand_inventory(run_id,ligand_instance_id,virus_name,protein_type,pdb_code,model_id,ligand_resname,ligand_chain,ligand_residue_id,ligand_insertion_code,ligand_atom_count,ligand_heavy_atom_count,centroid_x,centroid_y,centroid_z,pocket_radius_a,pocket_protein_atom_count,pocket_centroid_x,pocket_centroid_y,pocket_centroid_z,min_ligand_protein_distance_a,method_version,status) VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)""",
        "atoms": """INSERT OR REPLACE INTO ligand_atom_geometry(run_id,ligand_instance_id,ligand_instance_atom_id,atom_site_id,exact_atom,element,x,y,z,radial_distance_from_ligand_centroid_a,nearest_protein_distance_a,outward_cosine,outward_score,points_away_from_pocket,pocket_centroid_available,exit_vector_x,exit_vector_y,exit_vector_z,forward_clearance_a,forward_obstruction_count,exit_vector_clear,local_corridor_clear,forward_clearance_reaches_cap,exit_ray_length_a,exit_ray_radius_a,method_version) VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)""",
        "chains": """INSERT OR REPLACE INTO target_chain_geometry(run_id,ligand_instance_id,virus_name,protein_type,pdb_code,chain_id,model_id,chain_length_aa,lys_count,exposed_lys_count,exposed_lys_fraction,near_ligand_lys_count,near_ligand_exposed_lys_count,min_lys_ligand_distance_a,median_lys_ligand_distance_a,total_sasa_a2,lysine_sasa_a2,lysine_surface_fraction,isoelectric_point,basic_fraction,acidic_fraction,polar_fraction,hydrophobic_fraction,has_exposed_lysine,has_ligand_proximal_exposed_lysine,lysine_sidechain_sasa_a2,lysine_nz_sasa_a2,nz_observed_lys_count,nz_observed_lys_fraction,sidechain_accessible_lys_count,sidechain_accessible_lys_fraction,nz_exposed_lys_count_gt_0_1,nz_exposed_lys_fraction_gt_0_1,nz_exposed_lys_count_gt_1,nz_exposed_lys_fraction_gt_1,nz_exposed_lys_count_gt_5,nz_exposed_lys_fraction_gt_5,method_version) VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)""",
        "surface_lysines": """INSERT OR REPLACE INTO target_surface_lysines(run_id,ligand_instance_id,virus_name,protein_type,pdb_code,chain_id,model_id,lys_residue_id,lys_insertion_code,lys_observed_index,lysine_sasa_a2,is_surface_exposed,lysine_sidechain_sasa_a2,nz_sasa_a2,nz_observed,sidechain_accessible,nz_exposed_gt_0_1,nz_exposed_gt_1,nz_exposed_gt_5,reference_atom,reference_atom_x,reference_atom_y,reference_atom_z,method_version) VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)""",
        "pocket": """INSERT OR REPLACE INTO ligand_binding_pocket_atoms(run_id,ligand_instance_id,source_atom_site_id,label_asym_id,auth_asym_id,label_seq_id,auth_seq_id,insertion_code,residue_name,label_atom_id,auth_atom_id,element,x,y,z,distance_a,pocket_radius_a,method_version) VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)""",
    }

    with c.dbconn(database) as db:
        db.execute(f"PRAGMA busy_timeout={SQLITE_BUSY_TIMEOUT_MS}")
        executor = ProcessPoolExecutor(max_workers=max(1, workers)) if workers and workers > 1 else None
        futures = [executor.submit(_analyze_group, p) for p in payloads] if executor else None
        iterable = (
            (future.result() for future in as_completed(futures))
            if futures is not None
            else map(_analyze_group, payloads)
        )
        total_groups = len(payloads)
        for n, result in enumerate(iterable, 1):
            processed += len(result["iids"])
            if result["status"] == "completed":
                db.executemany(insert_sql["inventory"], result["inventory"])
                db.executemany(insert_sql["atoms"], result["atoms"])
                db.executemany(insert_sql["chains"], result["chains"])
                db.executemany(insert_sql["surface_lysines"], result["surface_lysines"])
                db.executemany(insert_sql["pocket"], result["pocket"])
                # Populate the compact foundation pocket table as well.
                db.executemany(
                    "INSERT INTO receptor_binding_pocket_atoms(run_id,ligand_instance_id,partner_atom_site_id,partner_label_asym_id,partner_auth_asym_id,partner_auth_seq_id,distance) VALUES(?,?,?,?,?,?,?)",
                    [(rid, r[1], r[2], r[3], r[4], r[6], r[15]) for r in result["pocket"]],
                )
                if result.get("target_context_status") == "not_applicable_no_protein_atoms":
                    not_applicable += len(result["iids"])
                    print(
                        f"geometry target not applicable: structure_id={result['structure_id']} "
                        f"model={result['model_id']} instances={len(result['iids'])} "
                        "reason=no recognized protein atoms",
                        flush=True,
                    )
                else:
                    success += len(result["iids"])
            else:
                failures += len(result["iids"])
                print(f"geometry group failed: structure_id={result['structure_id']} model={result['model_id']} instances={len(result['iids'])} error={result['error']}", flush=True)
                for iid in result["iids"]:
                    c.fail(db, rid, "geometry", result["error"], instance_id=iid, code="geometry_group_failure")
            if n % max(1, progress_every) == 0 or n == total_groups:
                db.commit()
                print(f"geometry progress: {n}/{total_groups} model-groups; instances={processed}/{len(ids)} success={success} not_applicable={not_applicable} failures={failures}", flush=True)
        if executor:
            executor.shutdown(wait=True)
        c.run_end(db, rid, "completed" if failures == 0 else "partial", processed, success, not_applicable, failures)

    _write_outputs(database)
    return {"run_id": rid, "processed": processed, "success": success, "not_applicable": not_applicable, "failures": failures, "model_groups": len(payloads)}


def main():
    p = argparse.ArgumentParser(description="Calculate occurrence-resolved ligand/target geometry and lysine accessibility.")
    p.add_argument("--database", default=str(c.ROOT / "viral_data_cif_v2.db"))
    p.add_argument("--limit", type=int)
    p.add_argument("--pdb-id")
    p.add_argument("--ligand-instance-id", type=int)
    p.add_argument("--workers", type=int, default=4)
    p.add_argument("--resume", action="store_true")
    p.add_argument("--pocket-radius", type=float, default=DEFAULT_POCKET_RADIUS)
    p.add_argument("--lysine-sasa-threshold", type=float, default=DEFAULT_LYSINE_SASA_THRESHOLD)
    p.add_argument("--exit-ray-length", type=float, default=DEFAULT_EXIT_RAY_LENGTH, help="Forward distance (A) used only for the local exit-vector clearance descriptor.")
    p.add_argument("--exit-ray-radius", type=float, default=DEFAULT_EXIT_RAY_RADIUS, help="Cylinder radius (A) used only for the local exit-vector clearance descriptor.")
    p.add_argument("--probe-radius", type=float, default=DEFAULT_PROBE_RADIUS)
    p.add_argument("--sasa-points", type=int, default=DEFAULT_SASA_POINTS)
    p.add_argument("--progress-every", type=int, default=25)
    a = p.parse_args()
    print(json.dumps(run(a.database, a.limit, a.pdb_id, a.ligand_instance_id, a.workers, a.resume,
                         a.pocket_radius, a.lysine_sasa_threshold, a.exit_ray_length,
                         a.exit_ray_radius, a.probe_radius, a.sasa_points, a.progress_every), indent=2))


if __name__ == "__main__":
    main()
