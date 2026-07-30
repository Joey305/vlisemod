#!/usr/bin/env python3
"""Schema helpers for V-LiSEMOD ligand SMILES graph enrichment.

The application database is provisioned outside the source tree. These helpers
apply or roll back the graph-enrichment tables on an explicitly supplied SQLite
database copy; they do not know about production or RANDY connection details.
"""

from __future__ import annotations

import sqlite3


GRAPH_TABLES = (
    "Ligand_SMILES_Graph_Assignments",
    "Ligand_SMILES_Bonds",
    "Ligand_SMILES_Atoms",
    "Ligand_SMILES_Graphs",
)


CREATE_STATEMENTS = (
    """
    CREATE TABLE IF NOT EXISTS Ligand_SMILES_Graphs (
        graph_id INTEGER PRIMARY KEY AUTOINCREMENT,
        source_smiles TEXT NOT NULL,
        canonical_smiles TEXT,
        isomeric_smiles TEXT,
        smiles_hash TEXT NOT NULL,
        atom_count INTEGER NOT NULL DEFAULT 0,
        heavy_atom_count INTEGER NOT NULL DEFAULT 0,
        bond_count INTEGER NOT NULL DEFAULT 0,
        formal_charge INTEGER NOT NULL DEFAULT 0,
        rdkit_valid INTEGER NOT NULL DEFAULT 0,
        parse_status TEXT NOT NULL,
        parse_message TEXT,
        graph_method TEXT NOT NULL,
        graph_version TEXT NOT NULL,
        rdkit_version TEXT,
        generated_at TEXT NOT NULL,
        UNIQUE (smiles_hash, graph_version)
    )
    """,
    """
    CREATE TABLE IF NOT EXISTS Ligand_SMILES_Atoms (
        graph_atom_id INTEGER PRIMARY KEY AUTOINCREMENT,
        graph_id INTEGER NOT NULL,
        smiles_atom_index INTEGER NOT NULL,
        element TEXT,
        atomic_number INTEGER,
        formal_charge INTEGER,
        isotope INTEGER,
        is_aromatic INTEGER NOT NULL DEFAULT 0,
        is_in_ring INTEGER NOT NULL DEFAULT 0,
        hybridization TEXT,
        chiral_tag TEXT,
        degree INTEGER,
        total_valence INTEGER,
        explicit_h_count INTEGER,
        implicit_h_count INTEGER,
        atom_map_number INTEGER,
        FOREIGN KEY (graph_id) REFERENCES Ligand_SMILES_Graphs(graph_id) ON DELETE CASCADE,
        UNIQUE (graph_id, smiles_atom_index)
    )
    """,
    """
    CREATE TABLE IF NOT EXISTS Ligand_SMILES_Bonds (
        graph_bond_id INTEGER PRIMARY KEY AUTOINCREMENT,
        graph_id INTEGER NOT NULL,
        smiles_bond_index INTEGER NOT NULL,
        begin_atom_index INTEGER NOT NULL,
        end_atom_index INTEGER NOT NULL,
        bond_type TEXT,
        bond_order REAL,
        is_aromatic INTEGER NOT NULL DEFAULT 0,
        is_conjugated INTEGER NOT NULL DEFAULT 0,
        is_in_ring INTEGER NOT NULL DEFAULT 0,
        stereo TEXT,
        bond_direction TEXT,
        FOREIGN KEY (graph_id) REFERENCES Ligand_SMILES_Graphs(graph_id) ON DELETE CASCADE,
        UNIQUE (graph_id, smiles_bond_index),
        UNIQUE (graph_id, begin_atom_index, end_atom_index),
        CHECK (begin_atom_index < end_atom_index)
    )
    """,
    """
    CREATE TABLE IF NOT EXISTS Ligand_SMILES_Graph_Assignments (
        assignment_id INTEGER PRIMARY KEY AUTOINCREMENT,
        graph_id INTEGER NOT NULL,
        pdb_code TEXT,
        model_id INTEGER,
        ligand_chain TEXT,
        ligand_residue_id INTEGER,
        ligand_resname TEXT,
        ligand_insertion_code TEXT,
        smiles_source_table TEXT NOT NULL,
        smiles_source_row_id INTEGER,
        smiles_source_column TEXT NOT NULL,
        source_smiles_hash TEXT NOT NULL,
        assignment_status TEXT NOT NULL,
        mapping_status TEXT NOT NULL,
        pdb_to_smiles_mapped_atom_count INTEGER,
        pdb_ligand_atom_count INTEGER,
        pdb_ligand_heavy_atom_count INTEGER,
        pdb_component_validation_status TEXT NOT NULL DEFAULT 'not_checked',
        pdb_component_validation_message TEXT,
        graph_method TEXT NOT NULL,
        graph_version TEXT NOT NULL,
        generated_at TEXT NOT NULL,
        FOREIGN KEY (graph_id) REFERENCES Ligand_SMILES_Graphs(graph_id) ON DELETE CASCADE,
        UNIQUE (
            smiles_source_table,
            smiles_source_row_id,
            smiles_source_column,
            graph_version
        )
    )
    """,
)


INDEX_STATEMENTS = (
    "CREATE INDEX IF NOT EXISTS idx_smiles_graph_hash ON Ligand_SMILES_Graphs(smiles_hash);",
    "CREATE INDEX IF NOT EXISTS idx_smiles_atoms_graph ON Ligand_SMILES_Atoms(graph_id);",
    "CREATE INDEX IF NOT EXISTS idx_smiles_atoms_element ON Ligand_SMILES_Atoms(graph_id, element);",
    "CREATE INDEX IF NOT EXISTS idx_smiles_bonds_graph ON Ligand_SMILES_Bonds(graph_id);",
    "CREATE INDEX IF NOT EXISTS idx_smiles_bonds_begin ON Ligand_SMILES_Bonds(graph_id, begin_atom_index);",
    "CREATE INDEX IF NOT EXISTS idx_smiles_bonds_end ON Ligand_SMILES_Bonds(graph_id, end_atom_index);",
    "CREATE INDEX IF NOT EXISTS idx_smiles_assign_graph ON Ligand_SMILES_Graph_Assignments(graph_id);",
    """
    CREATE INDEX IF NOT EXISTS idx_smiles_assign_ligand_instance
    ON Ligand_SMILES_Graph_Assignments(
        pdb_code,
        model_id,
        ligand_resname,
        ligand_chain,
        ligand_residue_id
    );
    """,
    """
    CREATE INDEX IF NOT EXISTS idx_smiles_assign_source
    ON Ligand_SMILES_Graph_Assignments(
        smiles_source_table,
        smiles_source_row_id,
        smiles_source_column
    );
    """,
)


def apply_schema(conn: sqlite3.Connection) -> None:
    """Create graph-enrichment tables and indexes."""
    conn.execute("PRAGMA foreign_keys = ON")
    for statement in CREATE_STATEMENTS:
        conn.execute(statement)
    for statement in INDEX_STATEMENTS:
        conn.execute(statement)


def rollback_schema(conn: sqlite3.Connection) -> None:
    """Drop graph-enrichment tables in dependency order."""
    conn.execute("PRAGMA foreign_keys = OFF")
    for table_name in GRAPH_TABLES:
        conn.execute(f'DROP TABLE IF EXISTS "{table_name}"')

