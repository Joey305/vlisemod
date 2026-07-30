#!/usr/bin/env python3
"""Query helpers for PROTACability attachment-site enrichment tables."""

from __future__ import annotations

import json
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


def get_analysis(
    conn: sqlite3.Connection,
    *,
    pdb_code: str,
    ligand_resname: str,
    ligand_chain: str,
    ligand_residue_id: int,
    model_id: int = 0,
    ligand_insertion_code: str = "",
    method_version: str = "attachment_v1",
) -> dict[str, Any] | None:
    rows = _rows(
        conn,
        """
        SELECT *
        FROM protacability_attachment_analysis
        WHERE pdb_code=?
          AND model_id=?
          AND ligand_chain=?
          AND ligand_residue_id=?
          AND ligand_insertion_code=?
          AND ligand_resname=?
          AND method_version=?
        """,
        (pdb_code, model_id, ligand_chain, ligand_residue_id, ligand_insertion_code, ligand_resname, method_version),
    )
    return dict(rows[0]) if rows else None


def get_attachment_atoms(conn: sqlite3.Connection, analysis_id: int) -> list[dict[str, Any]]:
    return [
        dict(row)
        for row in _rows(
            conn,
            """
            SELECT *
            FROM protacability_attachment_atoms
            WHERE analysis_id=?
            ORDER BY attachment_score DESC, pdb_atom_name
            """,
            (analysis_id,),
        )
    ]


def get_attachment_regions(conn: sqlite3.Connection, analysis_id: int) -> list[dict[str, Any]]:
    return [
        dict(row)
        for row in _rows(
            conn,
            """
            SELECT *
            FROM protacability_attachment_regions
            WHERE analysis_id=?
            ORDER BY region_score DESC, region_id
            """,
            (analysis_id,),
        )
    ]


def get_best_attachment_atoms(conn: sqlite3.Connection, analysis_id: int, limit: int = 10) -> list[dict[str, Any]]:
    return [
        dict(row)
        for row in _rows(
            conn,
            """
            SELECT *
            FROM protacability_attachment_atoms
            WHERE analysis_id=? AND candidate_attachment_flag=1
            ORDER BY attachment_score DESC, pdb_atom_name
            LIMIT ?
            """,
            (analysis_id, limit),
        )
    ]


def atoms_in_region(conn: sqlite3.Connection, analysis_id: int, region_id: str) -> list[dict[str, Any]]:
    return [
        dict(row)
        for row in _rows(
            conn,
            """
            SELECT *
            FROM protacability_attachment_atoms
            WHERE analysis_id=? AND region_id=?
            ORDER BY attachment_score DESC, pdb_atom_name
            """,
            (analysis_id, region_id),
        )
    ]


def parse_json_fields(row: dict[str, Any], *field_names: str) -> dict[str, Any]:
    parsed = dict(row)
    for field_name in field_names:
        value = parsed.get(field_name)
        if isinstance(value, str) and value:
            parsed[field_name] = json.loads(value)
    return parsed


def region_adjacency(conn: sqlite3.Connection, analysis_id: int) -> dict[str, set[str]]:
    rows = _rows(
        conn,
        """
        SELECT region_id, member_smiles_indices_json
        FROM protacability_attachment_regions
        WHERE analysis_id=?
        """,
        (analysis_id,),
    )
    members = {row["region_id"]: set(json.loads(row["member_smiles_indices_json"])) for row in rows}
    analysis = _rows(
        conn,
        "SELECT graph_id FROM protacability_attachment_analysis WHERE analysis_id=?",
        (analysis_id,),
    )
    if not analysis:
        return {}
    graph_id = analysis[0]["graph_id"]
    bonds = _rows(
        conn,
        """
        SELECT begin_atom_index, end_atom_index
        FROM Ligand_SMILES_Bonds
        WHERE graph_id=?
        """,
        (graph_id,),
    )
    out = {region_id: set() for region_id in members}
    for left_id, left_atoms in members.items():
        for right_id, right_atoms in members.items():
            if left_id >= right_id:
                continue
            linked = any(
                (bond["begin_atom_index"] in left_atoms and bond["end_atom_index"] in right_atoms)
                or (bond["end_atom_index"] in left_atoms and bond["begin_atom_index"] in right_atoms)
                for bond in bonds
            )
            if linked:
                out[left_id].add(right_id)
                out[right_id].add(left_id)
    return out


def shortest_region_distance(
    conn: sqlite3.Connection,
    analysis_id: int,
    source_region_id: str,
    target_region_id: str,
) -> int | None:
    if source_region_id == target_region_id:
        return 0
    adj = region_adjacency(conn, analysis_id)
    if source_region_id not in adj or target_region_id not in adj:
        return None
    queue = deque([(source_region_id, 0)])
    seen = {source_region_id}
    while queue:
        region_id, distance = queue.popleft()
        for neighbor in adj.get(region_id, set()):
            if neighbor == target_region_id:
                return distance + 1
            if neighbor not in seen:
                seen.add(neighbor)
                queue.append((neighbor, distance + 1))
    return None
