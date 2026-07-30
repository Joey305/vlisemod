#!/usr/bin/env python3
"""Query helpers for normalized ligand SMILES graph tables."""

from __future__ import annotations

import sqlite3
from collections import deque
from typing import Any


def _rows(conn: sqlite3.Connection, query: str, params: tuple[Any, ...]) -> list[sqlite3.Row]:
    old_factory = conn.row_factory
    conn.row_factory = sqlite3.Row
    try:
        return conn.execute(query, params).fetchall()
    finally:
        conn.row_factory = old_factory


def get_bonds(conn: sqlite3.Connection, graph_id: int) -> list[dict[str, Any]]:
    return [
        dict(row)
        for row in _rows(
            conn,
            """
            SELECT *
            FROM Ligand_SMILES_Bonds
            WHERE graph_id=?
            ORDER BY smiles_bond_index
            """,
            (graph_id,),
        )
    ]


def get_neighbors(conn: sqlite3.Connection, graph_id: int, smiles_atom_index: int) -> list[dict[str, Any]]:
    rows = _rows(
        conn,
        """
        SELECT b.*, a.smiles_atom_index AS neighbor_atom_index, a.element AS neighbor_element
        FROM Ligand_SMILES_Bonds b
        JOIN Ligand_SMILES_Atoms a
          ON a.graph_id = b.graph_id
         AND a.smiles_atom_index = CASE
             WHEN b.begin_atom_index = ? THEN b.end_atom_index
             ELSE b.begin_atom_index
         END
        WHERE b.graph_id = ?
          AND (? IN (b.begin_atom_index, b.end_atom_index))
        ORDER BY neighbor_atom_index
        """,
        (smiles_atom_index, graph_id, smiles_atom_index),
    )
    return [dict(row) for row in rows]


def _adjacency(conn: sqlite3.Connection, graph_id: int) -> dict[int, set[int]]:
    atom_rows = _rows(
        conn,
        "SELECT smiles_atom_index FROM Ligand_SMILES_Atoms WHERE graph_id=?",
        (graph_id,),
    )
    adj = {int(row["smiles_atom_index"]): set() for row in atom_rows}
    for bond in get_bonds(conn, graph_id):
        a = int(bond["begin_atom_index"])
        b = int(bond["end_atom_index"])
        adj.setdefault(a, set()).add(b)
        adj.setdefault(b, set()).add(a)
    return adj


def shortest_graph_distance(
    conn: sqlite3.Connection,
    graph_id: int,
    begin_atom_index: int,
    end_atom_index: int,
) -> int | None:
    if begin_atom_index == end_atom_index:
        return 0
    adj = _adjacency(conn, graph_id)
    if begin_atom_index not in adj or end_atom_index not in adj:
        return None
    queue = deque([(begin_atom_index, 0)])
    seen = {begin_atom_index}
    while queue:
        atom, distance = queue.popleft()
        for neighbor in adj.get(atom, set()):
            if neighbor == end_atom_index:
                return distance + 1
            if neighbor not in seen:
                seen.add(neighbor)
                queue.append((neighbor, distance + 1))
    return None


def atoms_within_n_bonds(
    conn: sqlite3.Connection,
    graph_id: int,
    smiles_atom_index: int,
    max_distance: int,
) -> list[int]:
    adj = _adjacency(conn, graph_id)
    if smiles_atom_index not in adj:
        return []
    queue = deque([(smiles_atom_index, 0)])
    seen = {smiles_atom_index}
    out = []
    while queue:
        atom, distance = queue.popleft()
        if distance <= max_distance:
            out.append(atom)
        if distance == max_distance:
            continue
        for neighbor in sorted(adj.get(atom, set())):
            if neighbor not in seen:
                seen.add(neighbor)
                queue.append((neighbor, distance + 1))
    return sorted(out)


def connected_components(conn: sqlite3.Connection, graph_id: int) -> list[list[int]]:
    adj = _adjacency(conn, graph_id)
    seen = set()
    components = []
    for atom in sorted(adj):
        if atom in seen:
            continue
        queue = deque([atom])
        seen.add(atom)
        component = []
        while queue:
            current = queue.popleft()
            component.append(current)
            for neighbor in sorted(adj.get(current, set())):
                if neighbor not in seen:
                    seen.add(neighbor)
                    queue.append(neighbor)
        components.append(component)
    return components


def ring_atoms(conn: sqlite3.Connection, graph_id: int) -> list[int]:
    rows = _rows(
        conn,
        """
        SELECT smiles_atom_index
        FROM Ligand_SMILES_Atoms
        WHERE graph_id=? AND is_in_ring=1
        ORDER BY smiles_atom_index
        """,
        (graph_id,),
    )
    return [int(row["smiles_atom_index"]) for row in rows]


def terminal_atoms(conn: sqlite3.Connection, graph_id: int) -> list[int]:
    rows = _rows(
        conn,
        """
        SELECT smiles_atom_index
        FROM Ligand_SMILES_Atoms
        WHERE graph_id=? AND degree <= 1
        ORDER BY smiles_atom_index
        """,
        (graph_id,),
    )
    return [int(row["smiles_atom_index"]) for row in rows]


def map_pdb_atom_to_graph_atom(
    conn: sqlite3.Connection,
    *,
    pdb_code: str,
    ligand_resname: str,
    ligand_chain: str,
    ligand_residue_id: int | None,
    pdb_atom_id: int,
    graph_version: str = "v1",
) -> dict[str, Any] | None:
    assignment = _rows(
        conn,
        """
        SELECT graph_id
        FROM Ligand_SMILES_Graph_Assignments
        WHERE pdb_code=?
          AND ligand_resname=?
          AND ligand_chain=?
          AND (? IS NULL OR ligand_residue_id=?)
          AND graph_version=?
        ORDER BY assignment_id
        LIMIT 1
        """,
        (pdb_code, ligand_resname, ligand_chain, ligand_residue_id, ligand_residue_id, graph_version),
    )
    if not assignment:
        return None
    graph_id = int(assignment[0]["graph_id"])
    mapped = _rows(
        conn,
        """
        SELECT exact_atom, atom_id, smiles_atom_index
        FROM SMILES_MAP_PDB
        WHERE pdb_id=?
          AND ligand=?
          AND chain=?
          AND atom_id=?
        LIMIT 1
        """,
        (pdb_code, ligand_resname, ligand_chain, pdb_atom_id),
    )
    if not mapped:
        return None
    return {"graph_id": graph_id, **dict(mapped[0])}

