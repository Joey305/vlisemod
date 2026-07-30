#!/usr/bin/env python3
"""Build PROTACability attachment-site evidence for bound ligand instances.

The default mode is a dry run. Writes require an explicit --write flag and an
explicit database path, so the script can be run safely against copied SQLite
databases during development.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import shlex
import sqlite3
import sys
import time
from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from protacability_attachment_schema import apply_schema


METHOD_VERSION = "attachment_v1_1"
GRAPH_VERSION = "v1"

STRONG_CONTACT_TYPES = {
    "polar",
    "hydrophobic",
    "vdw",
    "vdw_clash",
    "hbond",
    "ionic",
    "aromatic",
    "carbonyl",
    "halogen",
}

VDW_RADII = {
    "H": 1.20,
    "C": 1.70,
    "N": 1.55,
    "O": 1.52,
    "S": 1.80,
    "P": 1.80,
    "F": 1.47,
    "CL": 1.75,
    "BR": 1.85,
    "I": 1.98,
    "B": 1.92,
}


@dataclass(frozen=True)
class LigandInstance:
    source_row_id: int
    pdb_code: str
    model_id: int
    ligand_chain: str
    ligand_residue_id: int
    ligand_insertion_code: str
    ligand_resname: str
    graph_id: int
    mapping_status: str
    source_smiles: str
    atom_count: int
    heavy_atom_count: int
    bond_count: int
    collision_group_size: int = 1

    @property
    def key(self) -> tuple[str, int, str, int, str, str]:
        return (
            self.pdb_code,
            self.model_id,
            self.ligand_chain,
            self.ligand_residue_id,
            self.ligand_insertion_code,
            self.ligand_resname,
        )


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def json_text(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def norm_text(value: Any, default: str = "") -> str:
    if value is None:
        return default
    text = str(value).strip()
    return text if text else default


def safe_int(value: Any, default: int = 0) -> int:
    if value in (None, ""):
        return default
    try:
        return int(float(value))
    except Exception:
        return default


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--database", required=True, help="Explicit SQLite database path. Use a copied DB for writes.")
    parser.add_argument("--write", action="store_true", help="Actually write attachment-site tables. Omit for dry-run.")
    parser.add_argument("--dry-run", action="store_true", help="Explicit dry-run flag; this is also the default.")
    parser.add_argument("--pdb-code", help="Restrict to one PDB code.")
    parser.add_argument("--ligand-resname", help="Restrict to one ligand residue name.")
    parser.add_argument("--limit", type=int, help="Maximum unique ligand instances to process.")
    parser.add_argument("--method-version", default=METHOD_VERSION)
    parser.add_argument("--graph-version", default=GRAPH_VERSION)
    parser.add_argument(
        "--replace-version",
        nargs="?",
        const="__CURRENT__",
        help="In write mode, delete rows for this method version before rebuilding. If no value is supplied, uses --method-version.",
    )
    parser.add_argument("--report-path", help="Write JSON summary report to this path.")
    parser.add_argument("--failure-report", help="Write CSV failure report to this path.")
    parser.add_argument("--inventory-report", help="Write JSON inventory/classification report to this path.")
    parser.add_argument("--progress-path", help="Write machine-readable progress JSON during processing.")
    parser.add_argument("--batch-size", type=int, default=250, help="Commit interval in write mode.")
    parser.add_argument("--cif-root", default="PDB_FILES", help="Local mmCIF/PDB structure root for collision resolution.")
    parser.add_argument("--routine-density", type=int, default=960)
    parser.add_argument("--validation-density", type=int, default=1920)
    parser.add_argument("--use-validation-density", action="store_true")
    parser.add_argument("--region-graph-distance", type=int, default=2)
    parser.add_argument("--region-spatial-distance", type=float, default=5.0)
    parser.add_argument("--region-vector-dot-min", type=float, default=-0.15)
    parser.add_argument("--minimum-surface-sasa", type=float, default=2.0)
    return parser.parse_args()


def confidence(score: float | None, completeness: dict[str, Any]) -> str:
    if score is None:
        return "none"
    missing = int(completeness.get("missing_complex_sasa_atoms", 0))
    if score >= 0.70 and missing == 0:
        return "high"
    if score >= 0.45:
        return "medium"
    if score > 0:
        return "low"
    return "none"


def element_radius(element: str | None) -> float:
    return VDW_RADII.get(norm_text(element).upper(), 1.70)


def fibonacci_sphere(count: int) -> list[tuple[float, float, float]]:
    if count <= 1:
        return [(0.0, 0.0, 1.0)]
    points = []
    golden = math.pi * (3.0 - math.sqrt(5.0))
    for idx in range(count):
        y = 1.0 - (idx / float(count - 1)) * 2.0
        radius = math.sqrt(max(0.0, 1.0 - y * y))
        theta = golden * idx
        points.append((math.cos(theta) * radius, y, math.sin(theta) * radius))
    return points


def unit(vector: tuple[float, float, float]) -> tuple[float, float, float] | None:
    norm = math.sqrt(sum(v * v for v in vector))
    if norm == 0:
        return None
    return tuple(v / norm for v in vector)  # type: ignore[return-value]


def distance(a: dict[str, Any], b: dict[str, Any]) -> float:
    return math.sqrt((a["x"] - b["x"]) ** 2 + (a["y"] - b["y"]) ** 2 + (a["z"] - b["z"]) ** 2)


def discover_instances(conn: sqlite3.Connection, args: argparse.Namespace) -> list[LigandInstance]:
    clauses = ["a.smiles_source_table='Ligand_Atoms_Smiles'", "a.graph_version=?", "g.rdkit_valid=1"]
    params: list[Any] = [args.graph_version]
    if args.pdb_code:
        clauses.append("las.pdb_id=?")
        params.append(args.pdb_code)
    if args.ligand_resname:
        clauses.append("las.ligand=?")
        params.append(args.ligand_resname)
    rows = conn.execute(
        f"""
        SELECT
            MIN(las.rowid) AS source_row_id,
            las.pdb_id AS pdb_code,
            0 AS model_id,
            las.chain AS ligand_chain,
            las.ligand_id AS ligand_residue_id,
            '' AS ligand_insertion_code,
            las.ligand AS ligand_resname,
            a.graph_id,
            a.mapping_status,
            las.smiles AS source_smiles,
            g.atom_count,
            g.heavy_atom_count,
            g.bond_count
        FROM Ligand_Atoms_Smiles las
        JOIN Ligand_SMILES_Graph_Assignments a
          ON a.smiles_source_table='Ligand_Atoms_Smiles'
         AND a.smiles_source_row_id=las.rowid
        JOIN Ligand_SMILES_Graphs g ON g.graph_id=a.graph_id
        WHERE {" AND ".join(clauses)}
        GROUP BY las.pdb_id, las.chain, las.ligand_id, las.ligand, a.graph_id
        ORDER BY
            CASE WHEN las.pdb_id='3EKY' AND las.ligand='DR7' THEN 0 ELSE 1 END,
            CASE WHEN a.mapping_status='complete' THEN 0 ELSE 1 END,
            las.pdb_id, las.ligand, las.chain, las.ligand_id
        """,
        params,
    ).fetchall()
    instances = [
        LigandInstance(
            source_row_id=int(row["source_row_id"]),
            pdb_code=norm_text(row["pdb_code"]),
            model_id=safe_int(row["model_id"]),
            ligand_chain=norm_text(row["ligand_chain"]),
            ligand_residue_id=safe_int(row["ligand_residue_id"]),
            ligand_insertion_code=norm_text(row["ligand_insertion_code"]),
            ligand_resname=norm_text(row["ligand_resname"]),
            graph_id=int(row["graph_id"]),
            mapping_status=norm_text(row["mapping_status"], "unknown"),
            source_smiles=norm_text(row["source_smiles"]),
            atom_count=int(row["atom_count"]),
            heavy_atom_count=int(row["heavy_atom_count"]),
            bond_count=int(row["bond_count"]),
        )
        for row in rows
    ]
    collision_counts = Counter((inst.pdb_code, inst.model_id, inst.ligand_chain, inst.ligand_resname) for inst in instances)
    instances = [
        LigandInstance(
            source_row_id=inst.source_row_id,
            pdb_code=inst.pdb_code,
            model_id=inst.model_id,
            ligand_chain=inst.ligand_chain,
            ligand_residue_id=inst.ligand_residue_id,
            ligand_insertion_code=inst.ligand_insertion_code,
            ligand_resname=inst.ligand_resname,
            graph_id=inst.graph_id,
            mapping_status=inst.mapping_status,
            source_smiles=inst.source_smiles,
            atom_count=inst.atom_count,
            heavy_atom_count=inst.heavy_atom_count,
            bond_count=inst.bond_count,
            collision_group_size=collision_counts[(inst.pdb_code, inst.model_id, inst.ligand_chain, inst.ligand_resname)],
        )
        for inst in instances
    ]
    if not args.pdb_code and not args.ligand_resname and args.limit:
        return select_representative_pilot(conn, instances, args.limit)
    if args.limit:
        return instances[: args.limit]
    return instances


def select_representative_pilot(
    conn: sqlite3.Connection,
    instances: list[LigandInstance],
    limit: int,
) -> list[LigandInstance]:
    """Pick a deterministic diverse pilot set from primary ligand instances."""
    by_key = {inst.key: inst for inst in instances}
    selected: list[LigandInstance] = []
    selected_keys: set[tuple[str, int, str, int, str, str]] = set()

    def add(inst: LigandInstance | None) -> None:
        if inst and inst.key not in selected_keys and len(selected) < limit:
            selected.append(inst)
            selected_keys.add(inst.key)

    add(next((inst for inst in instances if inst.pdb_code == "3EKY" and inst.ligand_resname == "DR7"), None))

    def best(predicate, key):
        pool = [inst for inst in instances if inst.key not in selected_keys and predicate(inst)]
        return sorted(pool, key=key)[0] if pool else None

    add(best(lambda i: i.mapping_status == "complete" and i.atom_count >= 35, lambda i: (-i.atom_count, i.pdb_code)))
    add(best(lambda i: i.mapping_status == "complete" and i.atom_count <= 8, lambda i: (i.atom_count, i.pdb_code)))
    add(best(lambda i: i.mapping_status != "complete", lambda i: (i.mapping_status, i.pdb_code)))

    aromatic_counts = {
        int(row["graph_id"]): int(row["aromatic_count"])
        for row in conn.execute(
            """
            SELECT graph_id, SUM(is_aromatic) AS aromatic_count
            FROM Ligand_SMILES_Atoms
            GROUP BY graph_id
            """
        )
    }
    hetero_counts = {
        int(row["graph_id"]): int(row["hetero_count"])
        for row in conn.execute(
            """
            SELECT graph_id, SUM(CASE WHEN element NOT IN ('C','H') THEN 1 ELSE 0 END) AS hetero_count
            FROM Ligand_SMILES_Atoms
            GROUP BY graph_id
            """
        )
    }
    add(best(lambda i: aromatic_counts.get(i.graph_id, 0) >= 6, lambda i: (-aromatic_counts.get(i.graph_id, 0), i.pdb_code)))
    add(best(lambda i: hetero_counts.get(i.graph_id, 0) >= 5, lambda i: (-hetero_counts.get(i.graph_id, 0), i.pdb_code)))

    exposed_counts = defaultdict(int)
    for row in conn.execute(
        """
        SELECT pdb_id, ligand, chain, COUNT(*) AS exposed_count
        FROM solvent_exposed_atoms
        GROUP BY pdb_id, ligand, chain
        """
    ):
        exposed_counts[(row["pdb_id"], row["ligand"], row["chain"])] = int(row["exposed_count"])
    add(best(lambda i: exposed_counts.get((i.pdb_code, i.ligand_resname, i.ligand_chain), 0) >= 8,
             lambda i: (-exposed_counts.get((i.pdb_code, i.ligand_resname, i.ligand_chain), 0), i.pdb_code)))
    add(best(lambda i: exposed_counts.get((i.pdb_code, i.ligand_resname, i.ligand_chain), 0) == 0,
             lambda i: (i.atom_count, i.pdb_code)))

    chain_counts = Counter((i.pdb_code, i.ligand_resname) for i in instances)
    add(best(lambda i: chain_counts[(i.pdb_code, i.ligand_resname)] > 1, lambda i: (i.pdb_code, i.ligand_resname, i.ligand_chain)))

    for inst in instances:
        add(inst)
        if len(selected) >= limit:
            break
    return selected


def load_graph_adjacency(conn: sqlite3.Connection, graph_id: int) -> dict[int, set[int]]:
    adj = {
        int(row["smiles_atom_index"]): set()
        for row in conn.execute(
            "SELECT smiles_atom_index FROM Ligand_SMILES_Atoms WHERE graph_id=?",
            (graph_id,),
        )
    }
    for row in conn.execute(
        "SELECT begin_atom_index, end_atom_index FROM Ligand_SMILES_Bonds WHERE graph_id=?",
        (graph_id,),
    ):
        a = int(row["begin_atom_index"])
        b = int(row["end_atom_index"])
        adj.setdefault(a, set()).add(b)
        adj.setdefault(b, set()).add(a)
    return adj


def graph_distance(adj: dict[int, set[int]], source: int, target: int) -> int | None:
    if source == target:
        return 0
    queue = deque([(source, 0)])
    seen = {source}
    while queue:
        atom, dist = queue.popleft()
        for neighbor in sorted(adj.get(atom, set())):
            if neighbor == target:
                return dist + 1
            if neighbor not in seen:
                seen.add(neighbor)
                queue.append((neighbor, dist + 1))
    return None


@dataclass
class InstanceResolution:
    status: str
    method: str
    ambiguity_flag: int
    coordinate_source: str
    mapping_source: str
    notes: str
    structure_path: str | None = None
    atom_records: list[dict[str, Any]] | None = None


def find_structure_file(cif_root: Path, pdb_code: str) -> Path | None:
    pdb = pdb_code.upper()
    matches = sorted(cif_root.glob(f"**/{pdb}.cif"))
    if matches:
        return matches[0]
    matches = sorted(cif_root.glob(f"**/{pdb}.pdb"))
    return matches[0] if matches else None


def parse_mmcif_atom_site(cif_path: Path) -> list[dict[str, Any]]:
    lines = cif_path.read_text(errors="ignore").splitlines()
    atoms: list[dict[str, Any]] = []
    i = 0
    while i < len(lines):
        if lines[i].strip() != "loop_":
            i += 1
            continue
        i += 1
        headers = []
        while i < len(lines) and lines[i].strip().startswith("_"):
            headers.append(lines[i].strip())
            i += 1
        if not headers or not any(header.startswith("_atom_site.") for header in headers):
            continue
        header_index = {name: idx for idx, name in enumerate(headers)}
        while i < len(lines):
            line = lines[i].strip()
            if not line:
                i += 1
                continue
            if line == "#" or line == "loop_" or line.startswith("_"):
                break
            parts = shlex.split(line)
            if len(parts) >= len(headers):
                def get(name: str, default: str = "") -> str:
                    idx = header_index.get(name)
                    if idx is None or idx >= len(parts):
                        return default
                    value = parts[idx]
                    return "" if value in {".", "?"} else value
                atoms.append({
                    "group_PDB": get("_atom_site.group_PDB"),
                    "atom_id": safe_int(get("_atom_site.id"), -1),
                    "element": get("_atom_site.type_symbol"),
                    "atom_name": get("_atom_site.auth_atom_id") or get("_atom_site.label_atom_id"),
                    "resname": get("_atom_site.auth_comp_id") or get("_atom_site.label_comp_id"),
                    "chain": get("_atom_site.auth_asym_id") or get("_atom_site.label_asym_id"),
                    "residue_id": safe_int(get("_atom_site.auth_seq_id") or get("_atom_site.label_seq_id"), -999999),
                    "insertion_code": get("_atom_site.pdbx_PDB_ins_code"),
                    "model_id": safe_int(get("_atom_site.pdbx_PDB_model_num"), 1) - 1,
                    "x": float(get("_atom_site.Cartn_x", "nan")),
                    "y": float(get("_atom_site.Cartn_y", "nan")),
                    "z": float(get("_atom_site.Cartn_z", "nan")),
                })
            i += 1
    return atoms


def resolve_instance(
    conn: sqlite3.Connection,
    inst: LigandInstance,
    cif_root: Path,
    structure_cache: dict[str, tuple[Path | None, list[dict[str, Any]]]],
) -> InstanceResolution:
    if inst.collision_group_size <= 1:
        return InstanceResolution(
            status="resolved",
            method="existing_database_identity",
            ambiguity_flag=0,
            coordinate_source="ligand_atoms_by_unique_pdb_ligand_chain",
            mapping_source="SMILES_MAP_PDB_by_pdb_ligand_chain_atom_id",
            notes="No same-PDB/model/chain/ligand-resname collision in Ligand_Atoms_Smiles inventory.",
        )
    if inst.mapping_status == "no_pdb_to_smiles_mapping":
        return InstanceResolution(
            status="unresolved",
            method="safe_exclusion_no_mapping",
            ambiguity_flag=1,
            coordinate_source="none",
            mapping_source="none",
            notes="Repeated ligand instance has no PDB-to-SMILES mapping; excluded before atom-level analysis.",
        )

    pdb = inst.pdb_code.upper()
    if pdb not in structure_cache:
        path = find_structure_file(cif_root, pdb)
        structure_cache[pdb] = (path, parse_mmcif_atom_site(path) if path and path.suffix.lower() == ".cif" else [])
    path, atom_site = structure_cache[pdb]
    if path is None or not atom_site:
        return InstanceResolution(
            status="unresolved",
            method="safe_exclusion_missing_structure_file",
            ambiguity_flag=1,
            coordinate_source="none",
            mapping_source="none",
            notes="Same-ligand/same-chain collision cannot be separated without a local mmCIF atom_site record.",
        )

    residue_atoms = [
        atom for atom in atom_site
        if atom["chain"] == inst.ligand_chain
        and atom["resname"] == inst.ligand_resname
        and atom["residue_id"] == inst.ligand_residue_id
        and atom["insertion_code"] == inst.ligand_insertion_code
        and atom["element"].upper() != "H"
        and atom["atom_id"] >= 0
    ]
    if not residue_atoms:
        return InstanceResolution(
            status="unresolved",
            method="safe_exclusion_structure_residue_not_found",
            ambiguity_flag=1,
            coordinate_source=str(path),
            mapping_source="none",
            notes="No matching ligand residue found in mmCIF atom_site for the full instance key.",
            structure_path=str(path),
        )

    name_to_smiles: dict[str, int] = {}
    for row in conn.execute(
        """
        SELECT exact_atom, smiles_atom_index
        FROM SMILES_MAP_PDB
        WHERE pdb_id=? AND ligand=? AND chain=?
          AND smiles_atom_index IS NOT NULL
        """,
        (inst.pdb_code, inst.ligand_resname, inst.ligand_chain),
    ):
        name = norm_text(row["exact_atom"])
        if name and name not in name_to_smiles:
            name_to_smiles[name] = safe_int(row["smiles_atom_index"], -1)
    if not name_to_smiles:
        return InstanceResolution(
            status="unresolved",
            method="safe_exclusion_no_atom_name_mapping",
            ambiguity_flag=1,
            coordinate_source=str(path),
            mapping_source="none",
            notes="Structure residue exists but SMILES_MAP_PDB has no atom-name mapping for this ligand.",
            structure_path=str(path),
        )
    resolved_atoms = []
    for atom in residue_atoms:
        smiles_idx = name_to_smiles.get(atom["atom_name"])
        if smiles_idx is None or smiles_idx < 0:
            continue
        resolved_atoms.append({**atom, "smiles_atom_index": smiles_idx})
    if not resolved_atoms:
        return InstanceResolution(
            status="unresolved",
            method="safe_exclusion_structure_mapping_mismatch",
            ambiguity_flag=1,
            coordinate_source=str(path),
            mapping_source="SMILES_MAP_PDB_atom_name_to_smiles_index",
            notes="No structure atom names matched the PDB-to-SMILES atom-name mapping.",
            structure_path=str(path),
        )
    return InstanceResolution(
        status="resolved",
        method="structure_file_atom_site",
        ambiguity_flag=0,
        coordinate_source=str(path),
        mapping_source="SMILES_MAP_PDB_atom_name_to_smiles_index",
        notes=f"Resolved repeated ligand instance using mmCIF atom_site residue {inst.ligand_chain}/{inst.ligand_resname}/{inst.ligand_residue_id}.",
        structure_path=str(path),
        atom_records=resolved_atoms,
    )


def load_atom_context(
    conn: sqlite3.Connection,
    inst: LigandInstance,
    resolution: InstanceResolution,
) -> list[dict[str, Any]]:
    if resolution.atom_records is not None:
        atom_rows = []
        for atom in resolution.atom_records:
            serial = int(atom["atom_id"])
            graph_atom = conn.execute(
                """
                SELECT is_aromatic, is_in_ring, degree, element
                FROM Ligand_SMILES_Atoms
                WHERE graph_id=? AND smiles_atom_index=?
                """,
                (inst.graph_id, atom["smiles_atom_index"]),
            ).fetchone()
            if graph_atom is None:
                continue
            sasa = conn.execute(
                """
                SELECT SASA_Area
                FROM RUPLEY_SASA_DATA
                WHERE pdb_id=? AND ligand=? AND chain=? AND atom_id=?
                LIMIT 1
                """,
                (inst.pdb_code, inst.ligand_resname, inst.ligand_chain, serial),
            ).fetchone()
            exposed = conn.execute(
                """
                SELECT 1
                FROM solvent_exposed_atoms
                WHERE pdb_id=? AND ligand=? AND chain=? AND atom_id=?
                LIMIT 1
                """,
                (inst.pdb_code, inst.ligand_resname, inst.ligand_chain, serial),
            ).fetchone()
            db_atom = conn.execute(
                """
                SELECT rowid
                FROM ligand_atoms
                WHERE pdb_id=? AND ligand=? AND chain=? AND atom_id=?
                LIMIT 1
                """,
                (inst.pdb_code, inst.ligand_resname, inst.ligand_chain, serial),
            ).fetchone()
            atom_rows.append({
                "database_atom_id": None if db_atom is None else db_atom["rowid"],
                "pdb_atom_serial": serial,
                "pdb_atom_name": atom["atom_name"],
                "element": atom["element"] or graph_atom["element"],
                "x": atom["x"],
                "y": atom["y"],
                "z": atom["z"],
                "smiles_atom_index": atom["smiles_atom_index"],
                "is_aromatic": graph_atom["is_aromatic"],
                "is_in_ring": graph_atom["is_in_ring"],
                "degree": graph_atom["degree"],
                "complex_sasa": None if sasa is None else sasa["SASA_Area"],
                "existing_exposed_atom_status": 0 if exposed is None else 1,
            })
        return atom_rows

    rows = conn.execute(
        """
        SELECT
            la.rowid AS database_atom_id,
            la.atom_id AS pdb_atom_serial,
            la.exact_atom AS pdb_atom_name,
            la.atom_type AS element,
            la.x, la.y, la.z,
            m.smiles_atom_index,
            ga.is_aromatic,
            ga.is_in_ring,
            ga.degree,
            COALESCE(s.SASA_Area, NULL) AS complex_sasa,
            CASE WHEN se.atom_id IS NULL THEN 0 ELSE 1 END AS existing_exposed_atom_status
        FROM ligand_atoms la
        JOIN SMILES_MAP_PDB m
          ON m.pdb_id=la.pdb_id
         AND m.ligand=la.ligand
         AND m.chain=la.chain
         AND m.atom_id=la.atom_id
        JOIN Ligand_SMILES_Atoms ga
          ON ga.graph_id=?
         AND ga.smiles_atom_index=m.smiles_atom_index
        LEFT JOIN RUPLEY_SASA_DATA s
          ON s.pdb_id=la.pdb_id
         AND s.ligand=la.ligand
         AND s.chain=la.chain
         AND s.atom_id=la.atom_id
        LEFT JOIN solvent_exposed_atoms se
          ON se.pdb_id=la.pdb_id
         AND se.ligand=la.ligand
         AND se.chain=la.chain
         AND se.atom_id=la.atom_id
        WHERE la.pdb_id=? AND la.ligand=? AND la.chain=?
          AND UPPER(COALESCE(la.atom_type, '')) != 'H'
        ORDER BY m.smiles_atom_index
        """,
        (inst.graph_id, inst.pdb_code, inst.ligand_resname, inst.ligand_chain),
    ).fetchall()
    atoms = []
    seen_serials = set()
    for row in rows:
        serial = int(row["pdb_atom_serial"])
        if serial in seen_serials:
            continue
        seen_serials.add(serial)
        atoms.append(dict(row))
    return atoms


def contact_cache_key(inst: LigandInstance) -> tuple[str, str, str, int]:
    return (inst.pdb_code, inst.ligand_resname, inst.ligand_chain, inst.ligand_residue_id)


def group_cache_key(inst: LigandInstance) -> tuple[str, str, str]:
    return (inst.pdb_code, inst.ligand_resname, inst.ligand_chain)


def empty_contact() -> dict[str, Any]:
    return {"strong": 0, "weak": 0, "types": set(), "nearest": None}


def add_contact(contact: dict[str, Any], contact_type: str, distance_value: Any) -> None:
    contact_types = {part.strip().lower() for part in contact_type.split(",") if part.strip()}
    is_strong = bool(contact_types & STRONG_CONTACT_TYPES)
    if is_strong:
        contact["strong"] += 1
    else:
        contact["weak"] += 1
    contact["types"].update(sorted(contact_types))
    if distance_value is not None:
        current = contact["nearest"]
        contact["nearest"] = float(distance_value) if current is None else min(float(distance_value), float(current))


def build_contact_cache(
    conn: sqlite3.Connection,
    instances: list[LigandInstance],
) -> dict[tuple[str, str, str, int], dict[int, dict[str, Any]]]:
    wanted = {contact_cache_key(inst) for inst in instances}
    cache: dict[tuple[str, str, str, int], dict[int, dict[str, Any]]] = {
        key: defaultdict(empty_contact) for key in wanted
    }
    if not wanted:
        return {}
    pdb_codes = sorted({key[0] for key in wanted})
    placeholders = ",".join("?" for _ in pdb_codes)
    for row in conn.execute(
        f"""
        SELECT pdb_id, ligand, chain, ligand_id, atom_id, Contact, Distance
        FROM Arpeggio_Contacts_Data
        WHERE pdb_id IN ({placeholders})
        """,
        pdb_codes,
    ):
        key = (row["pdb_id"], row["ligand"], row["chain"], safe_int(row["ligand_id"]))
        if key not in wanted:
            continue
        add_contact(cache[key][int(row["atom_id"])], norm_text(row["Contact"], "unknown"), row["Distance"])
    return {key: dict(value) for key, value in cache.items()}


def build_functional_group_cache(
    conn: sqlite3.Connection,
    instances: list[LigandInstance],
) -> dict[tuple[str, str, str], dict[int, list[str]]]:
    wanted = {group_cache_key(inst) for inst in instances}
    groups: dict[tuple[str, str, str], dict[int, set[str]]] = {
        key: defaultdict(set) for key in wanted
    }
    if not wanted:
        return {}
    pdb_codes = sorted({key[0] for key in wanted})
    placeholders = ",".join("?" for _ in pdb_codes)
    for row in conn.execute(
        f"""
        SELECT pdb_id, ligand, chain, atom_id, functional_group
        FROM Functional_Group_Atoms
        WHERE pdb_id IN ({placeholders})
        """,
        pdb_codes,
    ):
        key = (row["pdb_id"], row["ligand"], row["chain"])
        if key in wanted and row["functional_group"]:
            groups[key][int(row["atom_id"])].add(str(row["functional_group"]))
    return {
        key: {atom_id: sorted(values) for atom_id, values in atom_groups.items()}
        for key, atom_groups in groups.items()
    }


def load_contacts(conn: sqlite3.Connection, inst: LigandInstance) -> dict[int, dict[str, Any]]:
    contacts: dict[int, dict[str, Any]] = defaultdict(lambda: {"strong": 0, "weak": 0, "types": set(), "nearest": None})
    for row in conn.execute(
        """
        SELECT atom_id, Contact, Distance
        FROM Arpeggio_Contacts_Data
        WHERE pdb_id=? AND ligand=? AND chain=? AND ligand_id=?
        """,
        (inst.pdb_code, inst.ligand_resname, inst.ligand_chain, inst.ligand_residue_id),
    ):
        atom_id = int(row["atom_id"])
        contact_type = norm_text(row["Contact"], "unknown")
        contact_types = {part.strip().lower() for part in contact_type.split(",") if part.strip()}
        is_strong = bool(contact_types & STRONG_CONTACT_TYPES)
        if is_strong:
            contacts[atom_id]["strong"] += 1
        else:
            contacts[atom_id]["weak"] += 1
        contacts[atom_id]["types"].update(sorted(contact_types))
        dist = row["Distance"]
        if dist is not None:
            current = contacts[atom_id]["nearest"]
            contacts[atom_id]["nearest"] = float(dist) if current is None else min(float(dist), float(current))
    return contacts


def load_functional_groups(conn: sqlite3.Connection, inst: LigandInstance) -> dict[int, list[str]]:
    groups: dict[int, set[str]] = defaultdict(set)
    for row in conn.execute(
        """
        SELECT atom_id, functional_group
        FROM Functional_Group_Atoms
        WHERE pdb_id=? AND ligand=? AND chain=?
        """,
        (inst.pdb_code, inst.ligand_resname, inst.ligand_chain),
    ):
        if row["functional_group"]:
            groups[int(row["atom_id"])].add(str(row["functional_group"]))
    return {atom_id: sorted(values) for atom_id, values in groups.items()}


def calculate_isolated_surface(
    atoms: list[dict[str, Any]],
    density: int,
    probe_radius: float = 1.4,
) -> None:
    points = fibonacci_sphere(density)
    radii = [element_radius(atom.get("element")) + probe_radius for atom in atoms]
    occluders: list[list[tuple[dict[str, Any], float]]] = []
    for idx, atom in enumerate(atoms):
        possible = []
        radius = radii[idx]
        for other_idx, other in enumerate(atoms):
            if other_idx == idx:
                continue
            other_radius = radii[other_idx]
            center_cutoff = radius + other_radius
            if distance(atom, other) <= center_cutoff:
                possible.append((other, other_radius * other_radius))
        occluders.append(possible)

    for idx, atom in enumerate(atoms):
        radius = radii[idx]
        exposed = 0
        vector_sum = [0.0, 0.0, 0.0]
        for px, py, pz in points:
            sx = atom["x"] + px * radius
            sy = atom["y"] + py * radius
            sz = atom["z"] + pz * radius
            occluded = False
            for other, other_radius_sq in occluders[idx]:
                if (sx - other["x"]) ** 2 + (sy - other["y"]) ** 2 + (sz - other["z"]) ** 2 < other_radius_sq:
                    occluded = True
                    break
            if not occluded:
                exposed += 1
                vector_sum[0] += px
                vector_sum[1] += py
                vector_sum[2] += pz
        fraction = exposed / float(density)
        atom["surface_point_count"] = exposed
        atom["isolated_ligand_sasa"] = fraction * 4.0 * math.pi * radius * radius
        atom["regenerated_exposed_atom_status"] = int(fraction >= 0.05)
        direction = unit((vector_sum[0], vector_sum[1], vector_sum[2]))
        if direction is None:
            atom["outward_vector"] = None
        else:
            atom["outward_vector"] = direction


def score_atoms(
    atoms: list[dict[str, Any]],
    contacts: dict[int, dict[str, Any]],
    groups: dict[int, list[str]],
    params: dict[str, Any],
) -> None:
    for atom in atoms:
        serial = int(atom["pdb_atom_serial"])
        contact = contacts.get(serial, {"strong": 0, "weak": 0, "types": set(), "nearest": None})
        atom["strong_contact_count"] = int(contact["strong"])
        atom["weak_contact_count"] = int(contact["weak"])
        atom["interaction_types"] = sorted(contact["types"])
        atom["nearest_protein_distance"] = contact["nearest"]
        atom["functional_group_annotations"] = groups.get(serial, [])
        atom["terminal_atom_flag"] = int(safe_int(atom.get("degree")) <= 1)
        complex_sasa = atom.get("complex_sasa")
        isolated_sasa = atom.get("isolated_ligand_sasa")
        if complex_sasa is not None and isolated_sasa and isolated_sasa > 0:
            relative = max(0.0, min(float(complex_sasa) / float(isolated_sasa), 1.0))
        else:
            relative = None
        atom["relative_exposure"] = relative

        surface_support = atom["surface_point_count"] / float(params["surface_density"])
        complex_norm = 0.0 if complex_sasa is None else max(0.0, min(float(complex_sasa) / 30.0, 1.0))
        relative_norm = 0.0 if relative is None else relative
        strong_penalty = min(atom["strong_contact_count"] / 8.0, 1.0)
        weak_penalty = min(atom["weak_contact_count"] / 20.0, 0.4)

        score = (
            0.34 * relative_norm
            + 0.24 * complex_norm
            + 0.14 * surface_support
            + 0.10 * atom["terminal_atom_flag"]
            + 0.10 * (1.0 - strong_penalty)
            + 0.08 * int(atom.get("existing_exposed_atom_status", 0))
            - 0.04 * weak_penalty
        )
        cautions = []
        reasons = []
        if complex_sasa is None:
            cautions.append("missing_complex_sasa")
        elif complex_sasa <= 0:
            cautions.append("zero_complex_sasa")
        else:
            reasons.append(f"complex_sasa={round(float(complex_sasa), 3)}")
        if isolated_sasa is None or isolated_sasa <= 0:
            cautions.append("zero_isolated_ligand_sasa")
        elif relative is not None:
            reasons.append(f"relative_exposure={round(relative, 3)}")
        if atom["surface_point_count"] > 0:
            reasons.append(f"surface_points={atom['surface_point_count']}")
        if atom["strong_contact_count"] > 0:
            cautions.append(f"strong_contacts={atom['strong_contact_count']}")
        if atom.get("is_aromatic"):
            cautions.append("aromatic_atom_candidate_caution")
        if atom.get("is_in_ring"):
            cautions.append("ring_atom_candidate_caution")
        if atom["terminal_atom_flag"]:
            reasons.append("terminal_graph_atom")
        if atom.get("functional_group_annotations"):
            reasons.append("functional_group_annotated")

        surface_defining = bool(
            (complex_sasa is not None and float(complex_sasa) >= float(params["minimum_surface_sasa"]))
            or (
                complex_sasa is not None
                and float(complex_sasa) > 0
                and int(atom.get("existing_exposed_atom_status", 0)) == 1
            )
        )
        candidate = bool(
            surface_defining
            and score >= 0.18
            and atom["strong_contact_count"] <= 10
            and isolated_sasa
            and isolated_sasa > 0
        )
        atom["attachment_score"] = round(max(0.0, min(score, 1.0)), 4)
        atom["surface_defining_flag"] = int(surface_defining)
        atom["candidate_attachment_flag"] = int(candidate)
        atom["reasons"] = reasons
        atom["cautions"] = cautions


def connected_surface_regions(
    atoms: list[dict[str, Any]],
    adj: dict[int, set[int]],
    params: dict[str, Any],
) -> list[list[dict[str, Any]]]:
    surface_atoms = [atom for atom in atoms if atom["surface_defining_flag"]]
    by_idx = {int(atom["smiles_atom_index"]): atom for atom in surface_atoms}
    region_adj: dict[int, set[int]] = {idx: set() for idx in by_idx}
    for left_idx, left in by_idx.items():
        for right_idx, right in by_idx.items():
            if left_idx >= right_idx:
                continue
            gd = graph_distance(adj, left_idx, right_idx)
            if gd is None or gd > int(params["region_graph_distance"]):
                continue
            if distance(left, right) > float(params["region_spatial_distance"]):
                continue
            left_vec = left.get("outward_vector")
            right_vec = right.get("outward_vector")
            if left_vec is not None and right_vec is not None:
                dot = sum(left_vec[i] * right_vec[i] for i in range(3))
                if dot < float(params["region_vector_dot_min"]):
                    continue
            region_adj[left_idx].add(right_idx)
            region_adj[right_idx].add(left_idx)

    seen: set[int] = set()
    regions = []
    for idx in sorted(region_adj):
        if idx in seen:
            continue
        queue = deque([idx])
        seen.add(idx)
        member_indices = []
        while queue:
            current = queue.popleft()
            member_indices.append(current)
            for neighbor in sorted(region_adj[current]):
                if neighbor not in seen:
                    seen.add(neighbor)
                    queue.append(neighbor)
        regions.append([by_idx[item] for item in sorted(member_indices)])

    atom_by_idx = {int(atom["smiles_atom_index"]): atom for atom in atoms}
    expanded_regions = []
    assigned = set()
    for region in regions:
        member = {int(atom["smiles_atom_index"]): atom for atom in region}
        for atom in atoms:
            idx = int(atom["smiles_atom_index"])
            if idx in member or atom["surface_defining_flag"]:
                continue
            if atom["attachment_score"] <= 0:
                continue
            distances = [graph_distance(adj, idx, int(seed["smiles_atom_index"])) for seed in region]
            finite = [d for d in distances if d is not None]
            if finite and min(finite) <= 1 and distance(atom, min(region, key=lambda seed: distance(atom, seed))) <= 3.2:
                member[idx] = atom
        expanded = [member[idx] for idx in sorted(member)]
        expanded_regions.append(expanded)
        assigned.update(int(atom["smiles_atom_index"]) for atom in expanded)
    return expanded_regions


def summarize_region(region: list[dict[str, Any]], region_id: str, completeness: dict[str, Any]) -> dict[str, Any]:
    candidates = [atom for atom in region if atom["candidate_attachment_flag"]]
    best = max(candidates or region, key=lambda atom: atom["attachment_score"])
    total_complex = sum(float(atom["complex_sasa"] or 0.0) for atom in region)
    total_isolated = sum(float(atom["isolated_ligand_sasa"] or 0.0) for atom in region)
    weighted_relative = total_complex / total_isolated if total_isolated > 0 else None
    centroid = tuple(sum(float(atom[axis]) for atom in region) / len(region) for axis in ("x", "y", "z"))
    vectors = [atom["outward_vector"] for atom in region if atom.get("outward_vector") is not None]
    if vectors:
        avg = tuple(sum(vec[i] for vec in vectors) / len(vectors) for i in range(3))
        avg_unit = unit(avg)
        coherence = math.sqrt(sum(v * v for v in avg))
    else:
        avg_unit = None
        coherence = None
    extent = max(math.sqrt(sum((float(atom[axis]) - centroid[i]) ** 2 for i, axis in enumerate(("x", "y", "z")))) for atom in region)
    nearest_values = [atom["nearest_protein_distance"] for atom in region if atom.get("nearest_protein_distance") is not None]
    strong_total = sum(atom["strong_contact_count"] for atom in region)
    weak_total = sum(atom["weak_contact_count"] for atom in region)
    interaction_types = sorted({item for atom in region for item in atom["interaction_types"]})
    score = (
        0.60 * max(atom["attachment_score"] for atom in region)
        + 0.20 * min(total_complex / 40.0, 1.0)
        + 0.10 * min(sum(atom["surface_point_count"] for atom in region) / 1000.0, 1.0)
        + 0.10 * (0.0 if coherence is None else coherence)
    )
    reasons = [
        f"member_atoms={len(region)}",
        f"candidate_atoms={len(candidates)}",
        f"total_complex_sasa={round(total_complex, 3)}",
    ]
    cautions = []
    if any(atom["is_aromatic"] for atom in region):
        cautions.append("region_contains_aromatic_atoms")
    if strong_total:
        cautions.append(f"region_strong_contacts={strong_total}")
    if completeness.get("missing_complex_sasa_atoms"):
        cautions.append("analysis_has_missing_complex_sasa_atoms")
    return {
        "region_id": region_id,
        "member_atom_ids": [int(atom["pdb_atom_serial"]) for atom in region],
        "member_smiles_indices": [int(atom["smiles_atom_index"]) for atom in region],
        "candidate_atom_ids": [int(atom["pdb_atom_serial"]) for atom in candidates],
        "best_candidate_atom_id": int(best["pdb_atom_serial"]) if best else None,
        "total_complex_sasa": round(total_complex, 4),
        "total_isolated_ligand_sasa": round(total_isolated, 4),
        "weighted_relative_exposure": None if weighted_relative is None else round(weighted_relative, 4),
        "surface_point_count": sum(atom["surface_point_count"] for atom in region),
        "centroid": centroid,
        "outward_vector": avg_unit,
        "vector_coherence": None if coherence is None else round(coherence, 4),
        "spatial_extent": round(extent, 4),
        "nearest_protein_distance": None if not nearest_values else round(min(nearest_values), 4),
        "interaction_summary": {"strong": strong_total, "weak": weak_total, "types": interaction_types},
        "region_score": round(max(0.0, min(score, 1.0)), 4),
        "confidence": confidence(score, completeness),
        "reasons": reasons,
        "cautions": cautions,
        "atoms": region,
    }


def assign_region_ids(regions: list[list[dict[str, Any]]], completeness: dict[str, Any]) -> list[dict[str, Any]]:
    summaries = [summarize_region(region, "pending", completeness) for region in regions]
    summaries.sort(
        key=lambda item: (
            -item["region_score"],
            -item["total_complex_sasa"],
            item["member_smiles_indices"],
        )
    )
    for idx, summary in enumerate(summaries, start=1):
        summary["region_id"] = f"ASR-{idx:02d}"
    return summaries


def prototype_region_mapping(region_summaries: list[dict[str, Any]]) -> dict[str, str]:
    out = {}
    for region in region_summaries:
        names = {atom["pdb_atom_name"] for atom in region["atoms"]}
        if {"CAA", "OBI"} & names:
            out[region["region_id"]] = "prototype_R01_like_terminal_ester_branch"
        elif {"CAB", "OBJ"} & names:
            out[region["region_id"]] = "prototype_R02_like_terminal_ester_branch"
        elif {"CAN", "CAQ", "CAU"} & names:
            out[region["region_id"]] = "prototype_R03_like_phenyl_edge"
        elif {"CAO", "CAS", "CAR", "NBD"} & names:
            out[region["region_id"]] = "prototype_R04_like_heteroaromatic_facing_region"
    return out


def evaluate_instance(
    conn: sqlite3.Connection,
    inst: LigandInstance,
    method_version: str,
    params: dict[str, Any],
    resolution: InstanceResolution,
    contact_cache: dict[tuple[str, str, str, int], dict[int, dict[str, Any]]] | None = None,
    functional_group_cache: dict[tuple[str, str, str], dict[int, list[str]]] | None = None,
) -> dict[str, Any]:
    if resolution.status != "resolved":
        return skipped_result(
            inst,
            method_version,
            params,
            resolution,
            "ambiguous_ligand_instance" if resolution.ambiguity_flag else "no_mapping",
            resolution.notes,
        )
    atoms = load_atom_context(conn, inst, resolution)
    completeness = {
        "graph_assignment_mapping_status": inst.mapping_status,
        "mapped_heavy_atoms": len(atoms),
        "graph_heavy_atoms": inst.heavy_atom_count,
        "complete_pdb_to_smiles_mapping": inst.mapping_status == "complete" and len(atoms) == inst.heavy_atom_count,
        "complex_sasa_atoms": sum(1 for atom in atoms if atom.get("complex_sasa") is not None),
        "missing_complex_sasa_atoms": sum(1 for atom in atoms if atom.get("complex_sasa") is None),
        "contact_source_rows": 0,
        "functional_group_atoms": 0,
        "coordinate_source": resolution.coordinate_source,
        "mapping_source": resolution.mapping_source,
        "instance_resolution_status": resolution.status,
        "instance_resolution_method": resolution.method,
        "instance_ambiguity_flag": resolution.ambiguity_flag,
        "resolution_notes": resolution.notes,
    }
    if not atoms:
        return skipped_result(inst, method_version, params, resolution, "no_mapping", "No mapped ligand heavy atoms were found for the instance.")
    contacts = (
        contact_cache.get(contact_cache_key(inst), {})
        if contact_cache is not None
        else load_contacts(conn, inst)
    )
    groups = (
        functional_group_cache.get(group_cache_key(inst), {})
        if functional_group_cache is not None
        else load_functional_groups(conn, inst)
    )
    completeness["contact_source_rows"] = sum(data["strong"] + data["weak"] for data in contacts.values())
    completeness["functional_group_atoms"] = sum(1 for atom in atoms if int(atom["pdb_atom_serial"]) in groups)

    calculate_isolated_surface(atoms, int(params["surface_density"]))
    score_atoms(atoms, contacts, groups, params)
    adj = load_graph_adjacency(conn, inst.graph_id)
    regions = connected_surface_regions(atoms, adj, params)
    region_summaries = assign_region_ids(regions, completeness)
    atom_region = {}
    for region in region_summaries:
        for atom in region["atoms"]:
            atom_region[int(atom["pdb_atom_serial"])] = region["region_id"]
            atom["vector_coherence"] = region["vector_coherence"]

    candidates = [atom for atom in atoms if atom["candidate_attachment_flag"]]
    best_atom = max(candidates, key=lambda atom: atom["attachment_score"]) if candidates else None
    best_score = best_atom["attachment_score"] if best_atom else None
    return {
        "instance": inst,
        "analysis_status": "completed",
        "eligibility_status": "full_analysis_eligible" if completeness["complete_pdb_to_smiles_mapping"] else "limited_analysis_eligible",
        "skip_reason": None,
        "mapping_status": inst.mapping_status,
        "source_data_completeness": completeness,
        "has_attachment_site_evidence": int(bool(region_summaries and candidates)),
        "attachment_region_count": len(region_summaries),
        "candidate_atom_count": len(candidates),
        "best_attachment_score": best_score,
        "best_attachment_confidence": confidence(best_score, completeness),
        "bond_source": "Ligand_SMILES_Bonds",
        "sasa_source": "RUPLEY_SASA_DATA plus deterministic ligand-only Shrake-Rupley estimate",
        "contact_source": "Arpeggio_Contacts_Data",
        "method_version": method_version,
        "calculation_parameters": params,
        "resolution": resolution,
        "atoms": atoms,
        "regions": region_summaries,
        "prototype_region_mapping": prototype_region_mapping(region_summaries)
        if inst.pdb_code == "3EKY" and inst.ligand_resname == "DR7"
        else {},
    }


def skipped_result(
    inst: LigandInstance,
    method_version: str,
    params: dict[str, Any],
    resolution: InstanceResolution,
    eligibility_status: str,
    skip_reason: str,
) -> dict[str, Any]:
    return {
        "instance": inst,
        "analysis_status": "skipped",
        "eligibility_status": eligibility_status,
        "skip_reason": skip_reason,
        "mapping_status": inst.mapping_status,
        "source_data_completeness": {
            "graph_assignment_mapping_status": inst.mapping_status,
            "mapped_heavy_atoms": 0,
            "graph_heavy_atoms": inst.heavy_atom_count,
            "complete_pdb_to_smiles_mapping": False,
            "complex_sasa_atoms": 0,
            "missing_complex_sasa_atoms": inst.heavy_atom_count,
            "contact_source_rows": 0,
            "functional_group_atoms": 0,
            "coordinate_source": resolution.coordinate_source,
            "mapping_source": resolution.mapping_source,
            "instance_resolution_status": resolution.status,
            "instance_resolution_method": resolution.method,
            "instance_ambiguity_flag": resolution.ambiguity_flag,
            "resolution_notes": resolution.notes,
        },
        "has_attachment_site_evidence": 0,
        "attachment_region_count": 0,
        "candidate_atom_count": 0,
        "best_attachment_score": None,
        "best_attachment_confidence": "none",
        "bond_source": "Ligand_SMILES_Bonds",
        "sasa_source": "not_evaluated",
        "contact_source": "not_evaluated",
        "method_version": method_version,
        "calculation_parameters": params,
        "resolution": resolution,
        "atoms": [],
        "regions": [],
        "prototype_region_mapping": {},
    }


def delete_version(conn: sqlite3.Connection, method_version: str) -> dict[str, int]:
    counts = {}
    for table in (
        "protacability_attachment_atoms",
        "protacability_attachment_regions",
        "protacability_attachment_analysis",
    ):
        counts[table] = conn.execute(
            f"SELECT COUNT(*) FROM {table} WHERE method_version=?",
            (method_version,),
        ).fetchone()[0]
    conn.execute("DELETE FROM protacability_attachment_analysis WHERE method_version=?", (method_version,))
    return counts


def write_result(conn: sqlite3.Connection, result: dict[str, Any], now: str) -> int:
    inst: LigandInstance = result["instance"]
    resolution: InstanceResolution = result["resolution"]
    cur = conn.execute(
        """
        INSERT OR REPLACE INTO protacability_attachment_analysis (
            pdb_code, model_id, ligand_chain, ligand_residue_id, ligand_insertion_code,
            ligand_resname, graph_id, analysis_status, mapping_status,
            source_data_completeness_json, has_attachment_site_evidence,
            attachment_region_count, candidate_atom_count, best_attachment_score,
            best_attachment_confidence, bond_source, sasa_source, contact_source,
            instance_resolution_status, instance_resolution_method,
            instance_ambiguity_flag, coordinate_source, mapping_source,
            resolution_notes, eligibility_status, skip_reason,
            software_versions_json, method_version, calculation_parameters_json,
            generated_at
        )
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            inst.pdb_code,
            inst.model_id,
            inst.ligand_chain,
            inst.ligand_residue_id,
            inst.ligand_insertion_code,
            inst.ligand_resname,
            inst.graph_id,
            result["analysis_status"],
            result["mapping_status"],
            json_text(result["source_data_completeness"]),
            result["has_attachment_site_evidence"],
            result["attachment_region_count"],
            result["candidate_atom_count"],
            result["best_attachment_score"],
            result["best_attachment_confidence"],
            result["bond_source"],
            result["sasa_source"],
            result["contact_source"],
            resolution.status,
            resolution.method,
            resolution.ambiguity_flag,
            resolution.coordinate_source,
            resolution.mapping_source,
            resolution.notes,
            result["eligibility_status"],
            result["skip_reason"],
            json_text({
                "python": sys.version.split()[0],
                "sqlite": sqlite3.sqlite_version,
            }),
            result["method_version"],
            json_text(result["calculation_parameters"]),
            now,
        ),
    )
    analysis_id = int(cur.lastrowid)
    for region in result["regions"]:
        centroid = region["centroid"]
        vector = region["outward_vector"]
        conn.execute(
            """
            INSERT INTO protacability_attachment_regions (
                analysis_id, region_id, member_atom_ids_json,
                member_smiles_indices_json, candidate_atom_ids_json,
                best_candidate_atom_id, total_complex_sasa,
                total_isolated_ligand_sasa, weighted_relative_exposure,
                surface_point_count, centroid_x, centroid_y, centroid_z,
                outward_vector_x, outward_vector_y, outward_vector_z,
                vector_coherence, spatial_extent, nearest_protein_distance,
                interaction_summary_json, region_score, confidence,
                reasons_json, cautions_json, method_version
            )
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                analysis_id,
                region["region_id"],
                json_text(region["member_atom_ids"]),
                json_text(region["member_smiles_indices"]),
                json_text(region["candidate_atom_ids"]),
                region["best_candidate_atom_id"],
                region["total_complex_sasa"],
                region["total_isolated_ligand_sasa"],
                region["weighted_relative_exposure"],
                region["surface_point_count"],
                centroid[0],
                centroid[1],
                centroid[2],
                None if vector is None else vector[0],
                None if vector is None else vector[1],
                None if vector is None else vector[2],
                region["vector_coherence"],
                region["spatial_extent"],
                region["nearest_protein_distance"],
                json_text(region["interaction_summary"]),
                region["region_score"],
                region["confidence"],
                json_text(region["reasons"]),
                json_text(region["cautions"]),
                result["method_version"],
            ),
        )
    region_by_atom = {
        int(atom["pdb_atom_serial"]): region["region_id"]
        for region in result["regions"]
        for atom in region["atoms"]
    }
    for atom in result["atoms"]:
        vector = atom.get("outward_vector")
        conn.execute(
            """
            INSERT INTO protacability_attachment_atoms (
                analysis_id, region_id, database_atom_id, pdb_atom_serial,
                pdb_atom_name, element, smiles_atom_index, complex_sasa,
                isolated_ligand_sasa, relative_exposure,
                existing_exposed_atom_status, regenerated_exposed_atom_status,
                surface_point_count, outward_vector_x, outward_vector_y,
                outward_vector_z, vector_coherence, strong_contact_count,
                weak_contact_count, interaction_types_json,
                functional_group_annotations_json, aromatic_flag, ring_flag,
                terminal_atom_flag, graph_degree, candidate_attachment_flag,
                surface_defining_flag, attachment_score, confidence,
                reasons_json, cautions_json, method_version
            )
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                analysis_id,
                region_by_atom.get(int(atom["pdb_atom_serial"])),
                atom["database_atom_id"],
                atom["pdb_atom_serial"],
                atom["pdb_atom_name"],
                atom["element"],
                atom["smiles_atom_index"],
                atom.get("complex_sasa"),
                atom.get("isolated_ligand_sasa"),
                atom.get("relative_exposure"),
                atom["existing_exposed_atom_status"],
                atom["regenerated_exposed_atom_status"],
                atom["surface_point_count"],
                None if vector is None else vector[0],
                None if vector is None else vector[1],
                None if vector is None else vector[2],
                atom.get("vector_coherence"),
                atom["strong_contact_count"],
                atom["weak_contact_count"],
                json_text(atom["interaction_types"]),
                json_text(atom["functional_group_annotations"]),
                atom["is_aromatic"],
                atom["is_in_ring"],
                atom["terminal_atom_flag"],
                atom["degree"],
                atom["candidate_attachment_flag"],
                atom["surface_defining_flag"],
                atom["attachment_score"],
                confidence(atom["attachment_score"], result["source_data_completeness"]),
                json_text(atom["reasons"]),
                json_text(atom["cautions"]),
                result["method_version"],
            ),
        )
    return analysis_id


def write_failure_report(path: str | None, failures: list[dict[str, Any]]) -> None:
    if not path:
        return
    fields = [
        "pdb_code",
        "model_id",
        "ligand_chain",
        "ligand_residue_id",
        "ligand_insertion_code",
        "ligand_resname",
        "graph_id",
        "mapping_status",
        "error",
    ]
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for failure in failures:
            writer.writerow({field: failure.get(field) for field in fields})


def audit_counts(conn: sqlite3.Connection, graph_version: str) -> dict[str, Any]:
    eligible = conn.execute(
        """
        SELECT COUNT(*) FROM (
            SELECT las.pdb_id, 0, las.chain, las.ligand_id, '', las.ligand
            FROM Ligand_Atoms_Smiles las
            JOIN Ligand_SMILES_Graph_Assignments a
              ON a.smiles_source_table='Ligand_Atoms_Smiles'
             AND a.smiles_source_row_id=las.rowid
             AND a.graph_version=?
            JOIN Ligand_SMILES_Graphs g ON g.graph_id=a.graph_id AND g.rdkit_valid=1
            GROUP BY las.pdb_id, las.chain, las.ligand_id, las.ligand
        )
        """,
        (graph_version,),
    ).fetchone()[0]
    complete = conn.execute(
        """
        SELECT COUNT(*) FROM (
            SELECT las.pdb_id, 0, las.chain, las.ligand_id, '', las.ligand
            FROM Ligand_Atoms_Smiles las
            JOIN Ligand_SMILES_Graph_Assignments a
              ON a.smiles_source_table='Ligand_Atoms_Smiles'
             AND a.smiles_source_row_id=las.rowid
             AND a.graph_version=?
             AND a.mapping_status='complete'
            JOIN Ligand_SMILES_Graphs g ON g.graph_id=a.graph_id AND g.rdkit_valid=1
            GROUP BY las.pdb_id, las.chain, las.ligand_id, las.ligand
        )
        """,
        (graph_version,),
    ).fetchone()[0]
    return {
        "instance_key": [
            "pdb_code",
            "model_id",
            "ligand_chain",
            "ligand_residue_id",
            "ligand_insertion_code",
            "ligand_resname",
        ],
        "primary_source": "Ligand_Atoms_Smiles joined to Ligand_SMILES_Graph_Assignments",
        "eligible_unique_ligand_instances": int(eligible),
        "eligible_complete_mapping_instances": int(complete),
        "duplicate_assignment_policy": "Functional_GROUPED and protacability summary assignments are graph references only; they are not attachment-analysis units.",
    }


def inventory_record(result: dict[str, Any]) -> dict[str, Any]:
    inst: LigandInstance = result["instance"]
    completeness = result["source_data_completeness"]
    return {
        "pdb_code": inst.pdb_code,
        "model_id": inst.model_id,
        "ligand_chain": inst.ligand_chain,
        "ligand_residue_id": inst.ligand_residue_id,
        "ligand_insertion_code": inst.ligand_insertion_code,
        "ligand_resname": inst.ligand_resname,
        "graph_id": inst.graph_id,
        "mapping_status": inst.mapping_status,
        "coordinate_available": bool(completeness.get("mapped_heavy_atoms", 0)),
        "stored_sasa_available": bool(completeness.get("complex_sasa_atoms", 0)),
        "contact_data_available": bool(completeness.get("contact_source_rows", 0)),
        "functional_group_available": bool(completeness.get("functional_group_atoms", 0)),
        "instance_resolution_status": result["resolution"].status,
        "instance_resolution_method": result["resolution"].method,
        "instance_ambiguity_flag": result["resolution"].ambiguity_flag,
        "eligibility_status": result["eligibility_status"],
        "analysis_status": result["analysis_status"],
        "skip_or_limitation_reason": result["skip_reason"],
        "attachment_region_count": result["attachment_region_count"],
        "candidate_atom_count": result["candidate_atom_count"],
        "best_attachment_score": result["best_attachment_score"],
        "best_attachment_confidence": result["best_attachment_confidence"],
    }


def main() -> int:
    args = parse_args()
    db_path = Path(args.database).expanduser().resolve()
    if not db_path.exists():
        raise SystemExit(f"Database does not exist: {db_path}")
    write = bool(args.write)
    method_version = args.method_version
    replace_version = None
    if args.replace_version:
        replace_version = method_version if args.replace_version == "__CURRENT__" else args.replace_version
    density = args.validation_density if args.use_validation_density else args.routine_density
    cif_root = Path(args.cif_root).expanduser().resolve()
    params = {
        "method_version": method_version,
        "graph_version": args.graph_version,
        "surface_density": density,
        "routine_density": args.routine_density,
        "validation_density": args.validation_density,
        "probe_radius": 1.4,
        "region_graph_distance": args.region_graph_distance,
        "region_spatial_distance": args.region_spatial_distance,
        "region_vector_dot_min": args.region_vector_dot_min,
        "minimum_surface_sasa": args.minimum_surface_sasa,
        "surface_region_method": "region_seeded_hybrid_graph_v1",
        "instance_resolution_method_version": "instance_resolution_v1_1",
    }
    now = utc_now()
    start_time = time.monotonic()
    failures: list[dict[str, Any]] = []
    results: list[dict[str, Any]] = []
    inventory: list[dict[str, Any]] = []

    with sqlite3.connect(db_path) as conn:
        conn.row_factory = sqlite3.Row
        conn.execute("PRAGMA foreign_keys = ON")
        if write:
            apply_schema(conn)
        if write and replace_version:
            replace_counts = delete_version(conn, replace_version)
        else:
            replace_counts = {}
        instances = discover_instances(conn, args)
        audit = audit_counts(conn, args.graph_version)
        contact_cache = build_contact_cache(conn, instances)
        functional_group_cache = build_functional_group_cache(conn, instances)
        structure_cache: dict[str, tuple[Path | None, list[dict[str, Any]]]] = {}
        total = len(instances)
        for pos, inst in enumerate(instances, start=1):
            savepoint_name = f"attachment_instance_{len(results) + len(failures)}"
            try:
                resolution = resolve_instance(conn, inst, cif_root, structure_cache)
                if write:
                    conn.execute(f"SAVEPOINT {savepoint_name}")
                result = evaluate_instance(
                    conn,
                    inst,
                    method_version,
                    params,
                    resolution,
                    contact_cache,
                    functional_group_cache,
                )
                if write:
                    write_result(conn, result, now)
                    conn.execute(f"RELEASE {savepoint_name}")
                results.append(result)
                inventory.append(inventory_record(result))
            except Exception as exc:
                if write:
                    conn.execute(f"ROLLBACK TO {savepoint_name}")
                    conn.execute(f"RELEASE {savepoint_name}")
                failure = {
                    "pdb_code": inst.pdb_code,
                    "model_id": inst.model_id,
                    "ligand_chain": inst.ligand_chain,
                    "ligand_residue_id": inst.ligand_residue_id,
                    "ligand_insertion_code": inst.ligand_insertion_code,
                    "ligand_resname": inst.ligand_resname,
                    "graph_id": inst.graph_id,
                    "mapping_status": inst.mapping_status,
                    "error": str(exc),
                }
                failures.append(failure)
                inventory.append({
                    **{k: failure[k] for k in ("pdb_code", "model_id", "ligand_chain", "ligand_residue_id", "ligand_insertion_code", "ligand_resname", "graph_id", "mapping_status")},
                    "eligibility_status": "failed_validation",
                    "analysis_status": "failed",
                    "skip_or_limitation_reason": str(exc),
                    "instance_resolution_status": "unknown",
                    "instance_resolution_method": "unknown",
                    "instance_ambiguity_flag": None,
                    "coordinate_available": False,
                    "stored_sasa_available": False,
                    "contact_data_available": False,
                    "functional_group_available": False,
                })
            if args.progress_path and (pos == total or pos % max(1, min(args.batch_size, 100)) == 0):
                elapsed = time.monotonic() - start_time
                rate = pos / elapsed if elapsed else 0
                progress = {
                    "total_eligible": total,
                    "completed": pos,
                    "successful": sum(1 for result in results if result["analysis_status"] == "completed"),
                    "limited": sum(1 for result in results if result["eligibility_status"] == "limited_analysis_eligible"),
                    "failed": len(failures),
                    "skipped": sum(1 for result in results if result["analysis_status"] == "skipped"),
                    "current_instance": inst.key,
                    "elapsed_seconds": elapsed,
                    "estimated_time_remaining_seconds": None if not rate else (total - pos) / rate,
                    "checkpoint_position": pos,
                    "current_database_size": db_path.stat().st_size,
                }
                Path(args.progress_path).write_text(json.dumps(progress, indent=2, sort_keys=True))
            if write and pos % args.batch_size == 0:
                conn.commit()
        if write:
            conn.commit()

    status_counts = Counter(result["analysis_status"] for result in results)
    eligibility_counts = Counter(result["eligibility_status"] for result in results)
    resolution_counts = Counter(result["resolution"].method for result in results)
    region_counts = Counter(result["attachment_region_count"] for result in results)
    candidate_counts = Counter(result["candidate_atom_count"] for result in results)
    pilot_selection = [
        {
            "pdb_code": result["instance"].pdb_code,
            "model_id": result["instance"].model_id,
            "ligand_chain": result["instance"].ligand_chain,
            "ligand_residue_id": result["instance"].ligand_residue_id,
            "ligand_insertion_code": result["instance"].ligand_insertion_code,
            "ligand_resname": result["instance"].ligand_resname,
            "graph_id": result["instance"].graph_id,
            "mapping_status": result["mapping_status"],
            "eligibility_status": result["eligibility_status"],
            "resolution_method": result["resolution"].method,
            "ambiguity_flag": result["resolution"].ambiguity_flag,
            "atom_count": result["instance"].atom_count,
            "region_count": result["attachment_region_count"],
            "candidate_atom_count": result["candidate_atom_count"],
            "best_attachment_score": result["best_attachment_score"],
        }
        for result in results
    ]
    three = next(
        (
            result for result in results
            if result["instance"].pdb_code == "3EKY" and result["instance"].ligand_resname == "DR7"
        ),
        None,
    )
    report = {
        "mode": "write" if write else "dry-run",
        "database": str(db_path),
        "method_version": method_version,
        "generated_at": now,
        "audit": audit,
        "replace_version": replace_version,
        "replace_counts": replace_counts,
        "instances_selected": len(results) + len(failures),
        "instances_completed": len(results),
        "instances_failed": len(failures),
        "analysis_status_counts": dict(status_counts),
        "eligibility_status_counts": dict(eligibility_counts),
        "resolution_method_counts": dict(resolution_counts),
        "attachment_region_count_distribution": dict(region_counts),
        "candidate_atom_count_distribution": dict(candidate_counts),
        "pilot_selection": pilot_selection,
        "failure_count": len(failures),
        "runtime_seconds": time.monotonic() - start_time,
        "calculation_parameters": params,
        "validation_3EKY_DR7": None if three is None else {
            "atom_count": three["instance"].atom_count,
            "heavy_atom_count": three["instance"].heavy_atom_count,
            "bond_count": three["instance"].bond_count,
            "mapped_heavy_atoms": three["source_data_completeness"]["mapped_heavy_atoms"],
            "attachment_region_count": three["attachment_region_count"],
            "candidate_atom_count": three["candidate_atom_count"],
            "best_attachment_score": three["best_attachment_score"],
            "prototype_region_mapping": three["prototype_region_mapping"],
            "atom_scores": {
                atom["pdb_atom_name"]: {
                    "smiles_atom_index": atom["smiles_atom_index"],
                    "complex_sasa": atom["complex_sasa"],
                    "relative_exposure": atom["relative_exposure"],
                    "surface_points": atom["surface_point_count"],
                    "score": atom["attachment_score"],
                    "region_id": next((region["region_id"] for region in three["regions"] if atom in region["atoms"]), None),
                    "surface_defining": atom["surface_defining_flag"],
                    "candidate": atom["candidate_attachment_flag"],
                    "cautions": atom["cautions"],
                }
                for atom in three["atoms"]
                if atom["pdb_atom_name"] in {"CAO", "CAS", "CAE", "CAG", "CAR", "CAA", "CAB", "OBI", "OBJ"}
            },
            "regions": [
                {
                    key: region[key]
                    for key in (
                        "region_id",
                        "member_atom_ids",
                        "member_smiles_indices",
                        "candidate_atom_ids",
                        "best_candidate_atom_id",
                        "total_complex_sasa",
                        "weighted_relative_exposure",
                        "region_score",
                        "confidence",
                        "cautions",
                    )
                }
                for region in three["regions"]
            ],
        },
        "example_atom_json": None if not results else {
            key: value
            for key, value in max(results[0]["atoms"], key=lambda atom: atom["attachment_score"]).items()
            if key not in {"outward_vector"}
        },
        "example_region_json": None if not results or not results[0]["regions"] else {
            key: value
            for key, value in results[0]["regions"][0].items()
            if key != "atoms"
        },
    }
    if args.report_path:
        Path(args.report_path).write_text(json.dumps(report, indent=2, sort_keys=True))
    if args.inventory_report:
        Path(args.inventory_report).write_text(json.dumps({
            "audit": audit,
            "method_version": method_version,
            "generated_at": now,
            "inventory": inventory,
            "summary": {
                "total": len(inventory),
                "eligibility_status_counts": dict(Counter(row["eligibility_status"] for row in inventory)),
                "analysis_status_counts": dict(Counter(row["analysis_status"] for row in inventory)),
            },
        }, indent=2, sort_keys=True))
    write_failure_report(args.failure_report, failures)
    print(json.dumps(report, indent=2, sort_keys=True))
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
