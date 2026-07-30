#!/usr/bin/env python3
"""Build normalized RDKit molecular graphs for stored V-LiSEMOD SMILES.

The script is safe by default:

* an explicit --database path is required
* dry-run is the default
* --write is required to create or replace rows
* no production or RANDY connection is inferred
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sqlite3
import sys
from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

from ligand_smiles_graph_schema import apply_schema

try:
    from rdkit import Chem, rdBase
except Exception as exc:  # pragma: no cover - exercised by users without RDKit.
    Chem = None
    rdBase = None
    RDKIT_IMPORT_ERROR = exc
else:
    RDKIT_IMPORT_ERROR = None


GRAPH_METHOD = "rdkit_source_smiles_graph"
DEFAULT_GRAPH_VERSION = "v1"


@dataclass(frozen=True)
class SourceRecord:
    source_table: str
    source_row_id: int
    source_column: str
    source_smiles: str
    pdb_code: str | None = None
    model_id: int | None = None
    ligand_chain: str | None = None
    ligand_residue_id: int | None = None
    ligand_insertion_code: str | None = None
    ligand_resname: str | None = None

    @property
    def smiles_hash(self) -> str:
        return sha256_text(self.source_smiles)


def sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def norm_text(value: Any) -> str | None:
    if value is None:
        return None
    text = str(value).strip()
    return text or None


def safe_int(value: Any) -> int | None:
    if value in (None, ""):
        return None
    try:
        return int(float(value))
    except Exception:
        return None


def table_exists(conn: sqlite3.Connection, table_name: str) -> bool:
    row = conn.execute(
        "SELECT 1 FROM sqlite_master WHERE type='table' AND name=?",
        (table_name,),
    ).fetchone()
    return row is not None


def table_columns(conn: sqlite3.Connection, table_name: str) -> set[str]:
    return {row[1] for row in conn.execute(f'PRAGMA table_info("{table_name}")')}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--database", required=True, help="Explicit SQLite database path. Use a copied RANDY/local DB for writes.")
    parser.add_argument("--write", action="store_true", help="Actually write graph tables. Omit for dry-run.")
    parser.add_argument("--dry-run", action="store_true", help="Explicit dry-run flag; this is also the default when --write is omitted.")
    parser.add_argument("--limit", type=int, help="Maximum source SMILES records to process.")
    parser.add_argument("--pdb-code", help="Restrict to a PDB code.")
    parser.add_argument("--ligand-resname", help="Restrict to a ligand residue code.")
    parser.add_argument("--source-table", help="Restrict to one source table.")
    parser.add_argument("--graph-version", default=DEFAULT_GRAPH_VERSION)
    parser.add_argument(
        "--replace-version",
        nargs="?",
        const="__CURRENT__",
        help="In write mode, delete rows for the supplied graph version before rebuilding. "
        "If no value is supplied, uses --graph-version.",
    )
    parser.add_argument("--report-path", help="Write JSON summary report to this path.")
    parser.add_argument("--failure-report", help="Write CSV failure report to this path.")
    parser.add_argument("--batch-size", type=int, default=500)
    parser.add_argument(
        "--validate-pdb-components",
        action="store_true",
        help="Optionally compare mapped SMILES graph connectivity to local PDB CIF _chem_comp_bond records.",
    )
    parser.add_argument("--component-cif-root", default="PDB_FILES", help="Repository-relative or absolute CIF root for optional validation.")
    return parser.parse_args()


def require_rdkit() -> None:
    if Chem is None:
        raise SystemExit(
            "RDKit is required for ligand SMILES graph generation. "
            f"Import error: {RDKIT_IMPORT_ERROR}"
        )


def discover_source_records(conn: sqlite3.Connection, args: argparse.Namespace) -> list[SourceRecord]:
    """Discover current SMILES records without assuming every table exists."""
    sources: list[SourceRecord] = []
    allowed = {args.source_table} if args.source_table else None
    pdb_filter = norm_text(args.pdb_code)
    ligand_filter = norm_text(args.ligand_resname)

    def source_allowed(table: str) -> bool:
        return allowed is None or table == args.source_table

    if source_allowed("Ligand_Atoms_Smiles") and table_exists(conn, "Ligand_Atoms_Smiles"):
        rows = conn.execute(
            """
            SELECT rowid AS source_row_id, pdb_id, ligand, chain, ligand_id, smiles
            FROM Ligand_Atoms_Smiles
            WHERE smiles IS NOT NULL AND trim(smiles) <> ''
            ORDER BY pdb_id, ligand, chain, ligand_id, rowid
            """
        ).fetchall()
        for row in rows:
            if pdb_filter and row["pdb_id"] != pdb_filter:
                continue
            if ligand_filter and row["ligand"] != ligand_filter:
                continue
            sources.append(SourceRecord(
                source_table="Ligand_Atoms_Smiles",
                source_row_id=int(row["source_row_id"]),
                source_column="smiles",
                source_smiles=str(row["smiles"]).strip(),
                pdb_code=norm_text(row["pdb_id"]),
                model_id=None,
                ligand_chain=norm_text(row["chain"]),
                ligand_residue_id=safe_int(row["ligand_id"]),
                ligand_resname=norm_text(row["ligand"]),
            ))

    if source_allowed("Functional_GROUPED") and table_exists(conn, "Functional_GROUPED"):
        rows = conn.execute(
            """
            SELECT rowid AS source_row_id, pdb_id, ligand, smiles
            FROM Functional_GROUPED
            WHERE smiles IS NOT NULL AND trim(smiles) <> ''
            ORDER BY pdb_id, ligand, rowid
            """
        ).fetchall()
        for row in rows:
            if pdb_filter and row["pdb_id"] != pdb_filter:
                continue
            if ligand_filter and row["ligand"] != ligand_filter:
                continue
            sources.append(SourceRecord(
                source_table="Functional_GROUPED",
                source_row_id=int(row["source_row_id"]),
                source_column="smiles",
                source_smiles=str(row["smiles"]).strip(),
                pdb_code=norm_text(row["pdb_id"]),
                ligand_resname=norm_text(row["ligand"]),
            ))

    if source_allowed("protacability_warhead_linkability") and table_exists(conn, "protacability_warhead_linkability"):
        rows = conn.execute(
            """
            SELECT rowid AS source_row_id, pdb_code, model_id, ligand_resname,
                   ligand_chain, ligand_residue_id, ligand_insertion_code,
                   representative_smiles
            FROM protacability_warhead_linkability
            WHERE representative_smiles IS NOT NULL AND trim(representative_smiles) <> ''
            ORDER BY pdb_code, ligand_resname, ligand_chain, ligand_residue_id, rowid
            """
        ).fetchall()
        for row in rows:
            if pdb_filter and row["pdb_code"] != pdb_filter:
                continue
            if ligand_filter and row["ligand_resname"] != ligand_filter:
                continue
            sources.append(SourceRecord(
                source_table="protacability_warhead_linkability",
                source_row_id=int(row["source_row_id"]),
                source_column="representative_smiles",
                source_smiles=str(row["representative_smiles"]).strip(),
                pdb_code=norm_text(row["pdb_code"]),
                model_id=safe_int(row["model_id"]),
                ligand_chain=norm_text(row["ligand_chain"]),
                ligand_residue_id=safe_int(row["ligand_residue_id"]),
                ligand_insertion_code=norm_text(row["ligand_insertion_code"]),
                ligand_resname=norm_text(row["ligand_resname"]),
            ))

    if args.limit:
        sources = sources[: args.limit]
    return sources


def parse_graph(source_smiles: str) -> dict[str, Any]:
    """Parse source SMILES while preserving RDKit source atom indices."""
    require_rdkit()
    parse_message = ""
    mol = None
    try:
        mol = Chem.MolFromSmiles(source_smiles, sanitize=True)
    except Exception as exc:
        parse_message = str(exc)

    if mol is None:
        return {
            "mol": None,
            "rdkit_valid": 0,
            "parse_status": "failed",
            "parse_message": parse_message or "RDKit returned no molecule",
            "canonical_smiles": None,
            "isomeric_smiles": None,
            "atom_count": 0,
            "heavy_atom_count": 0,
            "bond_count": 0,
            "formal_charge": 0,
            "atoms": [],
            "bonds": [],
        }

    atoms = []
    formal_charge = 0
    for atom in mol.GetAtoms():
        formal_charge += int(atom.GetFormalCharge())
        atoms.append({
            "smiles_atom_index": int(atom.GetIdx()),
            "element": atom.GetSymbol(),
            "atomic_number": int(atom.GetAtomicNum()),
            "formal_charge": int(atom.GetFormalCharge()),
            "isotope": int(atom.GetIsotope()),
            "is_aromatic": int(atom.GetIsAromatic()),
            "is_in_ring": int(atom.IsInRing()),
            "hybridization": str(atom.GetHybridization()),
            "chiral_tag": str(atom.GetChiralTag()),
            "degree": int(atom.GetDegree()),
            "total_valence": int(atom.GetTotalValence()),
            "explicit_h_count": int(atom.GetNumExplicitHs()),
            "implicit_h_count": int(atom.GetNumImplicitHs()),
            "atom_map_number": int(atom.GetAtomMapNum()),
        })

    bonds = []
    for bond in mol.GetBonds():
        begin = int(bond.GetBeginAtomIdx())
        end = int(bond.GetEndAtomIdx())
        lo, hi = sorted((begin, end))
        bonds.append({
            "smiles_bond_index": int(bond.GetIdx()),
            "begin_atom_index": lo,
            "end_atom_index": hi,
            "bond_type": str(bond.GetBondType()),
            "bond_order": float(bond.GetBondTypeAsDouble()),
            "is_aromatic": int(bond.GetIsAromatic()),
            "is_conjugated": int(bond.GetIsConjugated()),
            "is_in_ring": int(bond.IsInRing()),
            "stereo": str(bond.GetStereo()),
            "bond_direction": str(bond.GetBondDir()),
        })

    return {
        "mol": mol,
        "rdkit_valid": 1,
        "parse_status": "parsed",
        "parse_message": "",
        "canonical_smiles": Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False),
        "isomeric_smiles": Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True),
        "atom_count": int(mol.GetNumAtoms()),
        "heavy_atom_count": int(mol.GetNumHeavyAtoms()),
        "bond_count": int(mol.GetNumBonds()),
        "formal_charge": formal_charge,
        "atoms": atoms,
        "bonds": bonds,
    }


def get_or_create_graph(
    conn: sqlite3.Connection,
    source_smiles: str,
    graph_data: dict[str, Any],
    graph_version: str,
    write: bool,
    now: str,
) -> int | None:
    smiles_hash = sha256_text(source_smiles)
    existing = conn.execute(
        """
        SELECT graph_id
        FROM Ligand_SMILES_Graphs
        WHERE smiles_hash=? AND graph_version=?
        """,
        (smiles_hash, graph_version),
    ).fetchone()
    if existing:
        return int(existing["graph_id"])
    if not write:
        return None

    cur = conn.execute(
        """
        INSERT INTO Ligand_SMILES_Graphs (
            source_smiles, canonical_smiles, isomeric_smiles, smiles_hash,
            atom_count, heavy_atom_count, bond_count, formal_charge,
            rdkit_valid, parse_status, parse_message, graph_method,
            graph_version, rdkit_version, generated_at
        )
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            source_smiles,
            graph_data["canonical_smiles"],
            graph_data["isomeric_smiles"],
            smiles_hash,
            graph_data["atom_count"],
            graph_data["heavy_atom_count"],
            graph_data["bond_count"],
            graph_data["formal_charge"],
            graph_data["rdkit_valid"],
            graph_data["parse_status"],
            graph_data["parse_message"],
            GRAPH_METHOD,
            graph_version,
            rdBase.rdkitVersion if rdBase is not None else None,
            now,
        ),
    )
    graph_id = int(cur.lastrowid)

    for atom in graph_data["atoms"]:
        conn.execute(
            """
            INSERT INTO Ligand_SMILES_Atoms (
                graph_id, smiles_atom_index, element, atomic_number, formal_charge,
                isotope, is_aromatic, is_in_ring, hybridization, chiral_tag,
                degree, total_valence, explicit_h_count, implicit_h_count,
                atom_map_number
            )
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                graph_id,
                atom["smiles_atom_index"],
                atom["element"],
                atom["atomic_number"],
                atom["formal_charge"],
                atom["isotope"],
                atom["is_aromatic"],
                atom["is_in_ring"],
                atom["hybridization"],
                atom["chiral_tag"],
                atom["degree"],
                atom["total_valence"],
                atom["explicit_h_count"],
                atom["implicit_h_count"],
                atom["atom_map_number"],
            ),
        )

    for bond in graph_data["bonds"]:
        conn.execute(
            """
            INSERT INTO Ligand_SMILES_Bonds (
                graph_id, smiles_bond_index, begin_atom_index, end_atom_index,
                bond_type, bond_order, is_aromatic, is_conjugated,
                is_in_ring, stereo, bond_direction
            )
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                graph_id,
                bond["smiles_bond_index"],
                bond["begin_atom_index"],
                bond["end_atom_index"],
                bond["bond_type"],
                bond["bond_order"],
                bond["is_aromatic"],
                bond["is_conjugated"],
                bond["is_in_ring"],
                bond["stereo"],
                bond["bond_direction"],
            ),
        )

    return graph_id


def build_mapping_caches(conn: sqlite3.Connection) -> tuple[dict[tuple[str, str, str], set[int]], dict[tuple[str, str, str], int]]:
    mapping_cache: dict[tuple[str, str, str], set[int]] = defaultdict(set)
    if table_exists(conn, "SMILES_MAP_PDB"):
        for row in conn.execute(
            """
            SELECT pdb_id, ligand, chain, smiles_atom_index
            FROM SMILES_MAP_PDB
            WHERE smiles_atom_index IS NOT NULL
            """
        ):
            idx = safe_int(row["smiles_atom_index"])
            if idx is not None:
                mapping_cache[(row["pdb_id"], row["ligand"], row["chain"])].add(idx)

    ligand_count_cache: dict[tuple[str, str, str], int] = {}
    if table_exists(conn, "ligand_atoms"):
        for row in conn.execute(
            """
            SELECT pdb_id, ligand, chain, COUNT(*) AS atom_count
            FROM ligand_atoms
            GROUP BY pdb_id, ligand, chain
            """
        ):
            ligand_count_cache[(row["pdb_id"], row["ligand"], row["chain"])] = int(row["atom_count"])
    return dict(mapping_cache), ligand_count_cache


def mapping_status(
    conn: sqlite3.Connection,
    source: SourceRecord,
    graph_data: dict[str, Any],
    mapping_cache: dict[tuple[str, str, str], set[int]] | None = None,
    ligand_count_cache: dict[tuple[str, str, str], int] | None = None,
) -> tuple[str, int | None, int | None, int | None]:
    if not source.pdb_code or not source.ligand_resname or not source.ligand_chain:
        return "not_instance_specific", None, None, None
    key = (source.pdb_code, source.ligand_resname, source.ligand_chain)
    if mapping_cache is None:
        rows = conn.execute(
            """
            SELECT smiles_atom_index
            FROM SMILES_MAP_PDB
            WHERE pdb_id=? AND ligand=? AND chain=?
            ORDER BY smiles_atom_index
            """,
            key,
        ).fetchall()
        mapped = {safe_int(row["smiles_atom_index"]) for row in rows}
        mapped.discard(None)
    else:
        mapped = mapping_cache.get(key, set())
    if ligand_count_cache is None:
        ligand_count_row = conn.execute(
            """
            SELECT COUNT(*) AS atom_count
            FROM ligand_atoms
            WHERE pdb_id=? AND ligand=? AND chain=?
            """,
            key,
        ).fetchone()
        pdb_atom_count = int(ligand_count_row["atom_count"]) if ligand_count_row else None
    else:
        pdb_atom_count = ligand_count_cache.get(key)
    graph_atom_count = int(graph_data["atom_count"])
    if not mapped:
        status = "no_pdb_to_smiles_mapping"
    elif graph_atom_count and all(0 <= idx < graph_atom_count for idx in mapped):
        status = "complete" if len(mapped) == graph_atom_count else "partial"
    else:
        status = "index_out_of_range"
    return status, len(mapped), pdb_atom_count, pdb_atom_count


def parse_cif_component_bonds(cif_path: Path, ligand_resname: str) -> set[tuple[str, str]]:
    bonds: set[tuple[str, str]] = set()
    if not cif_path.exists():
        return bonds
    in_loop = False
    with cif_path.open(errors="ignore") as handle:
        for raw in handle:
            line = raw.strip()
            if line == "loop_":
                in_loop = False
                continue
            if line.startswith("_chem_comp_bond."):
                in_loop = True
                continue
            if in_loop:
                if line.startswith("_"):
                    break
                parts = line.split()
                if len(parts) >= 4 and parts[0] == ligand_resname:
                    a1, a2 = parts[1], parts[2]
                    if not a1.startswith("H") and not a2.startswith("H"):
                        bonds.add(tuple(sorted((a1, a2))))
    return bonds


def find_cif(root: Path, pdb_code: str) -> Path | None:
    pdb_upper = pdb_code.upper()
    matches = sorted(root.glob(f"**/{pdb_upper}.cif"))
    return matches[0] if matches else None


def validate_pdb_component(
    conn: sqlite3.Connection,
    source: SourceRecord,
    graph_data: dict[str, Any],
    cif_root: Path,
) -> tuple[str, str]:
    if not source.pdb_code or not source.ligand_resname or not source.ligand_chain:
        return "insufficient_mapping", "Source record is not ligand-instance-specific."
    cif_path = find_cif(cif_root, source.pdb_code)
    if cif_path is None:
        return "not_checked", "No local CIF found."
    component_bonds = parse_cif_component_bonds(cif_path, source.ligand_resname)
    if not component_bonds:
        return "insufficient_mapping", "No _chem_comp_bond records found for ligand."

    map_rows = conn.execute(
        """
        SELECT exact_atom, smiles_atom_index
        FROM SMILES_MAP_PDB
        WHERE pdb_id=? AND ligand=? AND chain=?
        """,
        (source.pdb_code, source.ligand_resname, source.ligand_chain),
    ).fetchall()
    idx_to_name = {
        safe_int(row["smiles_atom_index"]): str(row["exact_atom"]).strip()
        for row in map_rows
        if safe_int(row["smiles_atom_index"]) is not None
    }
    if not idx_to_name:
        return "insufficient_mapping", "No PDB-to-SMILES atom mapping."

    graph_bonds = set()
    for bond in graph_data["bonds"]:
        a = idx_to_name.get(bond["begin_atom_index"])
        b = idx_to_name.get(bond["end_atom_index"])
        if a and b:
            graph_bonds.add(tuple(sorted((a, b))))

    missing_from_graph = sorted(component_bonds - graph_bonds)
    missing_from_component = sorted(graph_bonds - component_bonds)
    if not missing_from_graph and not missing_from_component:
        return "consistent", f"Compared {len(graph_bonds)} mapped heavy-atom bonds against {cif_path}."
    if len(missing_from_graph) <= 2 and len(missing_from_component) <= 2:
        return "minor_representation_difference", (
            f"Small mapped-bond difference; missing_from_graph={missing_from_graph[:3]}, "
            f"missing_from_component={missing_from_component[:3]}"
        )
    return "connectivity_conflict", (
        f"Mapped-bond conflict; missing_from_graph={missing_from_graph[:6]}, "
        f"missing_from_component={missing_from_component[:6]}"
    )


def insert_assignment(
    conn: sqlite3.Connection,
    graph_id: int,
    source: SourceRecord,
    graph_data: dict[str, Any],
    graph_version: str,
    now: str,
    write: bool,
    validation: tuple[str, str],
    status_tuple: tuple[str, int | None, int | None, int | None] | None = None,
) -> None:
    status, mapped_count, pdb_atom_count, pdb_heavy_count = (
        status_tuple if status_tuple is not None else mapping_status(conn, source, graph_data)
    )
    if not write:
        return
    conn.execute(
        """
        INSERT OR REPLACE INTO Ligand_SMILES_Graph_Assignments (
            graph_id, pdb_code, model_id, ligand_chain, ligand_residue_id,
            ligand_resname, ligand_insertion_code, smiles_source_table,
            smiles_source_row_id, smiles_source_column, source_smiles_hash,
            assignment_status, mapping_status, pdb_to_smiles_mapped_atom_count,
            pdb_ligand_atom_count, pdb_ligand_heavy_atom_count,
            pdb_component_validation_status, pdb_component_validation_message,
            graph_method, graph_version, generated_at
        )
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            graph_id,
            source.pdb_code,
            source.model_id,
            source.ligand_chain,
            source.ligand_residue_id,
            source.ligand_resname,
            source.ligand_insertion_code,
            source.source_table,
            source.source_row_id,
            source.source_column,
            source.smiles_hash,
            "assigned" if graph_data["rdkit_valid"] else "parse_failed",
            status,
            mapped_count,
            pdb_atom_count,
            pdb_heavy_count,
            validation[0],
            validation[1],
            GRAPH_METHOD,
            graph_version,
            now,
        ),
    )


def replace_version(conn: sqlite3.Connection, graph_version: str) -> dict[str, int]:
    counts = {}
    for table_name in (
        "Ligand_SMILES_Graph_Assignments",
        "Ligand_SMILES_Bonds",
        "Ligand_SMILES_Atoms",
        "Ligand_SMILES_Graphs",
    ):
        row = conn.execute(
            f'SELECT COUNT(*) AS n FROM "{table_name}" WHERE graph_version=?'
            if table_name in {"Ligand_SMILES_Graph_Assignments", "Ligand_SMILES_Graphs"}
            else f"""
                SELECT COUNT(*) AS n
                FROM "{table_name}"
                WHERE graph_id IN (
                    SELECT graph_id FROM Ligand_SMILES_Graphs WHERE graph_version=?
                )
            """,
            (graph_version,),
        ).fetchone()
        counts[table_name] = int(row["n"])
    conn.execute("DELETE FROM Ligand_SMILES_Graph_Assignments WHERE graph_version=?", (graph_version,))
    conn.execute(
        """
        DELETE FROM Ligand_SMILES_Bonds
        WHERE graph_id IN (SELECT graph_id FROM Ligand_SMILES_Graphs WHERE graph_version=?)
        """,
        (graph_version,),
    )
    conn.execute(
        """
        DELETE FROM Ligand_SMILES_Atoms
        WHERE graph_id IN (SELECT graph_id FROM Ligand_SMILES_Graphs WHERE graph_version=?)
        """,
        (graph_version,),
    )
    conn.execute("DELETE FROM Ligand_SMILES_Graphs WHERE graph_version=?", (graph_version,))
    return counts


def connected_components_for_graph(graph_data: dict[str, Any]) -> list[list[int]]:
    nodes = {atom["smiles_atom_index"] for atom in graph_data["atoms"]}
    adj: dict[int, set[int]] = {node: set() for node in nodes}
    for bond in graph_data["bonds"]:
        a = bond["begin_atom_index"]
        b = bond["end_atom_index"]
        adj.setdefault(a, set()).add(b)
        adj.setdefault(b, set()).add(a)
    seen = set()
    comps = []
    for node in sorted(nodes):
        if node in seen:
            continue
        comp = []
        queue = deque([node])
        seen.add(node)
        while queue:
            cur = queue.popleft()
            comp.append(cur)
            for nb in sorted(adj.get(cur, set())):
                if nb not in seen:
                    seen.add(nb)
                    queue.append(nb)
        comps.append(comp)
    return comps


def write_failure_report(path: str | None, failures: list[dict[str, Any]]) -> None:
    if not path:
        return
    fieldnames = [
        "source_table", "source_row_id", "source_column", "pdb_code",
        "ligand_resname", "ligand_chain", "ligand_residue_id",
        "source_smiles", "parse_status", "parse_message",
    ]
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for failure in failures:
            writer.writerow({name: failure.get(name) for name in fieldnames})


def main() -> int:
    args = parse_args()
    require_rdkit()
    if args.write and args.dry_run:
        raise SystemExit("Choose either --write or --dry-run, not both.")
    db_path = Path(args.database).expanduser().resolve()
    if not db_path.exists():
        raise SystemExit(f"Database does not exist: {db_path}")
    if args.replace_version and not args.write:
        raise SystemExit("--replace-version requires --write")

    graph_version = args.graph_version
    replace_target = graph_version if args.replace_version == "__CURRENT__" else args.replace_version
    now = utc_now()
    summary: dict[str, Any] = {
        "database": str(db_path),
        "mode": "write" if args.write else "dry-run",
        "graph_method": GRAPH_METHOD,
        "graph_version": graph_version,
        "rdkit_version": rdBase.rdkitVersion if rdBase is not None else None,
        "atom_index_convention": "zero-based RDKit atom indices, matching existing SMILES_MAP_PDB.smiles_atom_index convention",
        "replace_version": replace_target,
    }

    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA foreign_keys = ON")
    try:
        if args.write:
            apply_schema(conn)
            if replace_target:
                summary["replace_counts"] = replace_version(conn, replace_target)
        else:
            summary["dry_run_note"] = "No schema or data writes were performed."

        sources = discover_source_records(conn, args)
        mapping_cache, ligand_count_cache = build_mapping_caches(conn)
        unique_by_hash: dict[tuple[str, str], str] = {}
        graph_cache: dict[tuple[str, str], dict[str, Any]] = {}
        failures = []
        assignments_seen = 0
        validation_counts: Counter[str] = Counter()
        mapping_counts: Counter[str] = Counter()
        source_counts: Counter[str] = Counter()
        duplicate_exact_count = 0

        for source in sources:
            source_counts[source.source_table] += 1
            key = (source.smiles_hash, graph_version)
            if key in unique_by_hash:
                duplicate_exact_count += 1
            else:
                unique_by_hash[key] = source.source_smiles
                graph_cache[key] = parse_graph(source.source_smiles)

        atom_total = sum(data["atom_count"] for data in graph_cache.values() if data["rdkit_valid"])
        bond_total = sum(data["bond_count"] for data in graph_cache.values() if data["rdkit_valid"])

        if args.write:
            for index, source in enumerate(sources, 1):
                key = (source.smiles_hash, graph_version)
                graph_data = graph_cache[key]
                graph_id = get_or_create_graph(conn, source.source_smiles, graph_data, graph_version, True, now)
                validation = ("not_checked", "")
                if args.validate_pdb_components and graph_data["rdkit_valid"]:
                    cif_root = Path(args.component_cif_root)
                    if not cif_root.is_absolute():
                        cif_root = Path.cwd() / cif_root
                    validation = validate_pdb_component(conn, source, graph_data, cif_root)
                validation_counts[validation[0]] += 1
                status_tuple = mapping_status(conn, source, graph_data, mapping_cache, ligand_count_cache)
                status, _, _, _ = status_tuple
                mapping_counts[status] += 1
                insert_assignment(conn, graph_id, source, graph_data, graph_version, now, True, validation, status_tuple)
                assignments_seen += 1
                if not graph_data["rdkit_valid"]:
                    failures.append({
                        "source_table": source.source_table,
                        "source_row_id": source.source_row_id,
                        "source_column": source.source_column,
                        "pdb_code": source.pdb_code,
                        "ligand_resname": source.ligand_resname,
                        "ligand_chain": source.ligand_chain,
                        "ligand_residue_id": source.ligand_residue_id,
                        "source_smiles": source.source_smiles,
                        "parse_status": graph_data["parse_status"],
                        "parse_message": graph_data["parse_message"],
                    })
                if index % max(args.batch_size, 1) == 0:
                    conn.commit()
            conn.commit()
        else:
            for source in sources:
                graph_data = graph_cache[(source.smiles_hash, graph_version)]
                status, _, _, _ = mapping_status(conn, source, graph_data, mapping_cache, ligand_count_cache)
                mapping_counts[status] += 1
                validation_counts["not_checked"] += 1
                if not graph_data["rdkit_valid"]:
                    failures.append({
                        "source_table": source.source_table,
                        "source_row_id": source.source_row_id,
                        "source_column": source.source_column,
                        "pdb_code": source.pdb_code,
                        "ligand_resname": source.ligand_resname,
                        "ligand_chain": source.ligand_chain,
                        "ligand_residue_id": source.ligand_residue_id,
                        "source_smiles": source.source_smiles,
                        "parse_status": graph_data["parse_status"],
                        "parse_message": graph_data["parse_message"],
                    })

        valid_graphs = sum(1 for data in graph_cache.values() if data["rdkit_valid"])
        invalid_graphs = len(graph_cache) - valid_graphs
        summary.update({
            "source_records_discovered": len(sources),
            "source_records_by_table": dict(source_counts),
            "unique_exact_smiles_graphs": len(graph_cache),
            "valid_graphs": valid_graphs,
            "invalid_or_failed_graphs": invalid_graphs,
            "duplicate_exact_smiles_source_records": duplicate_exact_count,
            "atom_total_valid_graphs": atom_total,
            "bond_total_valid_graphs": bond_total,
            "assignment_rows_processed": assignments_seen if args.write else len(sources),
            "mapping_status_counts": dict(mapping_counts),
            "pdb_component_validation_counts": dict(validation_counts),
            "failure_count": len(failures),
        })

        # 3EKY/DR7 validation is reported when present in the processed set.
        dr7_source = next(
            (
                s for s in sources
                if s.source_table == "Ligand_Atoms_Smiles"
                and s.pdb_code == "3EKY"
                and s.ligand_resname == "DR7"
                and s.ligand_chain == "A"
                and s.ligand_residue_id == 100
            ),
            None,
        )
        if dr7_source:
            dr7_graph = graph_cache[(dr7_source.smiles_hash, graph_version)]
            comps = connected_components_for_graph(dr7_graph) if dr7_graph["rdkit_valid"] else []
            dr7_status, mapped_count, pdb_atom_count, _ = mapping_status(conn, dr7_source, dr7_graph, mapping_cache, ligand_count_cache)
            aromatic_atoms = sum(1 for atom in dr7_graph["atoms"] if atom["is_aromatic"])
            ring_atoms = sum(1 for atom in dr7_graph["atoms"] if atom["is_in_ring"])
            summary["validation_3EKY_DR7"] = {
                "graph_valid": bool(dr7_graph["rdkit_valid"]),
                "atom_count": dr7_graph["atom_count"],
                "heavy_atom_count": dr7_graph["heavy_atom_count"],
                "bond_count": dr7_graph["bond_count"],
                "mapped_atom_count": mapped_count,
                "pdb_atom_count": pdb_atom_count,
                "mapping_status": dr7_status,
                "connected_component_count": len(comps),
                "connected_components": comps,
                "aromatic_atom_count": aromatic_atoms,
                "ring_atom_count": ring_atoms,
                "supports_all_51_mapped_heavy_atoms": bool(dr7_graph["heavy_atom_count"] == 51 and mapped_count == 51),
            }

        write_failure_report(args.failure_report, failures)
        if args.report_path:
            Path(args.report_path).write_text(json.dumps(summary, indent=2, sort_keys=True))
        print(json.dumps(summary, indent=2, sort_keys=True))
        return 0 if not failures else 2
    finally:
        conn.close()


if __name__ == "__main__":
    raise SystemExit(main())
