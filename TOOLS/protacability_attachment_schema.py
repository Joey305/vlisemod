#!/usr/bin/env python3
"""Schema helpers for PROTACability attachment-site enrichment.

The tables created here are intended for copied SQLite databases during
development. They store one analysis per bound ligand instance and method
version, evaluated atoms, and deterministic attachment regions.
"""

from __future__ import annotations

import sqlite3


ATTACHMENT_TABLES = (
    "protacability_attachment_atoms",
    "protacability_attachment_regions",
    "protacability_attachment_analysis",
)


CREATE_STATEMENTS = (
    """
    CREATE TABLE IF NOT EXISTS protacability_attachment_analysis (
        analysis_id INTEGER PRIMARY KEY AUTOINCREMENT,
        pdb_code TEXT NOT NULL,
        model_id INTEGER NOT NULL DEFAULT 0,
        ligand_chain TEXT NOT NULL,
        ligand_residue_id INTEGER NOT NULL,
        ligand_insertion_code TEXT NOT NULL DEFAULT '',
        ligand_resname TEXT NOT NULL,
        graph_id INTEGER,
        analysis_status TEXT NOT NULL,
        mapping_status TEXT NOT NULL,
        source_data_completeness_json TEXT NOT NULL,
        has_attachment_site_evidence INTEGER NOT NULL DEFAULT 0,
        attachment_region_count INTEGER NOT NULL DEFAULT 0,
        candidate_atom_count INTEGER NOT NULL DEFAULT 0,
        best_attachment_score REAL,
        best_attachment_confidence TEXT,
        bond_source TEXT NOT NULL,
        sasa_source TEXT NOT NULL,
        contact_source TEXT NOT NULL,
        method_version TEXT NOT NULL,
        calculation_parameters_json TEXT NOT NULL,
        generated_at TEXT NOT NULL,
        FOREIGN KEY (graph_id) REFERENCES Ligand_SMILES_Graphs(graph_id) ON DELETE SET NULL,
        UNIQUE (
            pdb_code,
            model_id,
            ligand_chain,
            ligand_residue_id,
            ligand_insertion_code,
            ligand_resname,
            method_version
        )
    )
    """,
    """
    CREATE TABLE IF NOT EXISTS protacability_attachment_regions (
        attachment_region_id INTEGER PRIMARY KEY AUTOINCREMENT,
        analysis_id INTEGER NOT NULL,
        region_id TEXT NOT NULL,
        member_atom_ids_json TEXT NOT NULL,
        member_smiles_indices_json TEXT NOT NULL,
        candidate_atom_ids_json TEXT NOT NULL,
        best_candidate_atom_id INTEGER,
        total_complex_sasa REAL,
        total_isolated_ligand_sasa REAL,
        weighted_relative_exposure REAL,
        surface_point_count INTEGER NOT NULL DEFAULT 0,
        centroid_x REAL,
        centroid_y REAL,
        centroid_z REAL,
        outward_vector_x REAL,
        outward_vector_y REAL,
        outward_vector_z REAL,
        vector_coherence REAL,
        spatial_extent REAL,
        nearest_protein_distance REAL,
        interaction_summary_json TEXT NOT NULL,
        region_score REAL,
        confidence TEXT,
        reasons_json TEXT NOT NULL,
        cautions_json TEXT NOT NULL,
        method_version TEXT NOT NULL,
        FOREIGN KEY (analysis_id) REFERENCES protacability_attachment_analysis(analysis_id) ON DELETE CASCADE,
        UNIQUE (analysis_id, region_id)
    )
    """,
    """
    CREATE TABLE IF NOT EXISTS protacability_attachment_atoms (
        attachment_atom_id INTEGER PRIMARY KEY AUTOINCREMENT,
        analysis_id INTEGER NOT NULL,
        region_id TEXT,
        database_atom_id INTEGER,
        pdb_atom_serial INTEGER,
        pdb_atom_name TEXT,
        element TEXT,
        smiles_atom_index INTEGER,
        complex_sasa REAL,
        isolated_ligand_sasa REAL,
        relative_exposure REAL,
        existing_exposed_atom_status INTEGER NOT NULL DEFAULT 0,
        regenerated_exposed_atom_status INTEGER NOT NULL DEFAULT 0,
        surface_point_count INTEGER NOT NULL DEFAULT 0,
        outward_vector_x REAL,
        outward_vector_y REAL,
        outward_vector_z REAL,
        vector_coherence REAL,
        strong_contact_count INTEGER NOT NULL DEFAULT 0,
        weak_contact_count INTEGER NOT NULL DEFAULT 0,
        interaction_types_json TEXT NOT NULL,
        functional_group_annotations_json TEXT NOT NULL,
        aromatic_flag INTEGER NOT NULL DEFAULT 0,
        ring_flag INTEGER NOT NULL DEFAULT 0,
        terminal_atom_flag INTEGER NOT NULL DEFAULT 0,
        graph_degree INTEGER,
        candidate_attachment_flag INTEGER NOT NULL DEFAULT 0,
        surface_defining_flag INTEGER NOT NULL DEFAULT 0,
        attachment_score REAL,
        confidence TEXT,
        reasons_json TEXT NOT NULL,
        cautions_json TEXT NOT NULL,
        method_version TEXT NOT NULL,
        FOREIGN KEY (analysis_id) REFERENCES protacability_attachment_analysis(analysis_id) ON DELETE CASCADE,
        UNIQUE (analysis_id, pdb_atom_serial)
    )
    """,
)


INDEX_STATEMENTS = (
    """
    CREATE INDEX IF NOT EXISTS idx_attachment_analysis_instance
    ON protacability_attachment_analysis (
        pdb_code,
        model_id,
        ligand_chain,
        ligand_residue_id,
        ligand_insertion_code,
        ligand_resname
    )
    """,
    "CREATE INDEX IF NOT EXISTS idx_attachment_analysis_graph ON protacability_attachment_analysis(graph_id)",
    "CREATE INDEX IF NOT EXISTS idx_attachment_atoms_analysis ON protacability_attachment_atoms(analysis_id)",
    "CREATE INDEX IF NOT EXISTS idx_attachment_atoms_smiles ON protacability_attachment_atoms(analysis_id, smiles_atom_index)",
    "CREATE INDEX IF NOT EXISTS idx_attachment_atoms_candidate ON protacability_attachment_atoms(analysis_id, candidate_attachment_flag)",
    "CREATE INDEX IF NOT EXISTS idx_attachment_regions_analysis ON protacability_attachment_regions(analysis_id)",
)


def apply_schema(conn: sqlite3.Connection) -> None:
    """Create attachment-site tables and indexes."""
    conn.execute("PRAGMA foreign_keys = ON")
    for statement in CREATE_STATEMENTS:
        conn.execute(statement)
    for statement in INDEX_STATEMENTS:
        conn.execute(statement)


def rollback_schema(conn: sqlite3.Connection) -> None:
    """Drop attachment-site tables in dependency order."""
    conn.execute("PRAGMA foreign_keys = OFF")
    for table_name in ATTACHMENT_TABLES:
        conn.execute(f'DROP TABLE IF EXISTS "{table_name}"')
