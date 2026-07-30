import json
import sqlite3
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT / "TOOLS"
BUILDER = TOOLS / "build_protacability_attachment_sites.py"
MIGRATION = TOOLS / "migrate_protacability_attachment_sites.py"

sys.path.insert(0, str(TOOLS))

from build_protacability_attachment_sites import score_atoms  # noqa: E402
from ligand_smiles_graph_schema import apply_schema as apply_graph_schema  # noqa: E402
from protacability_attachment_queries import (  # noqa: E402
    get_analysis,
    get_attachment_summary,
    get_attachment_atoms,
    get_attachment_regions,
)


class ProtacabilityAttachmentSiteTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.db_path = Path(self.tmp.name) / "sample.db"
        self._create_database()

    def tearDown(self):
        self.tmp.cleanup()

    def _create_database(self):
        with sqlite3.connect(self.db_path) as conn:
            conn.row_factory = sqlite3.Row
            apply_graph_schema(conn)
            conn.executescript(
                """
                CREATE TABLE Ligand_Atoms_Smiles (
                    virus_name TEXT,
                    pdb_id TEXT,
                    ligand TEXT,
                    chain TEXT,
                    ligand_id INTEGER,
                    smiles TEXT,
                    functional_groups TEXT,
                    molecular_weight REAL
                );
                CREATE TABLE ligand_atoms (
                    virus_name TEXT,
                    pdb_id TEXT,
                    ligand TEXT,
                    chain TEXT,
                    atom_id INTEGER,
                    exact_atom TEXT,
                    atom_type TEXT,
                    x REAL,
                    y REAL,
                    z REAL
                );
                CREATE TABLE SMILES_MAP_PDB (
                    virus_name TEXT,
                    pdb_id TEXT,
                    ligand TEXT,
                    chain TEXT,
                    exact_atom TEXT,
                    atom_id INTEGER,
                    atom_index INTEGER,
                    smiles_atom_index INTEGER
                );
                CREATE TABLE RUPLEY_SASA_DATA (
                    virus_name TEXT,
                    pdb_id TEXT,
                    ligand TEXT,
                    chain TEXT,
                    exact_atom TEXT,
                    atom_id INTEGER,
                    SASA_Area REAL
                );
                CREATE TABLE solvent_exposed_atoms (
                    virus_name TEXT,
                    pdb_id TEXT,
                    ligand TEXT,
                    chain TEXT,
                    atom_id INTEGER,
                    exact_atom TEXT,
                    atom_type TEXT,
                    x REAL,
                    y REAL,
                    z REAL
                );
                CREATE TABLE Arpeggio_Contacts_Data (
                    virus_name TEXT,
                    pdb_id TEXT,
                    ligand TEXT,
                    ligand_id INTEGER,
                    chain TEXT,
                    Contact TEXT,
                    Distance REAL,
                    exact_atom TEXT,
                    atom_id INTEGER,
                    residue TEXT,
                    residue_number INTEGER,
                    residue_atom TEXT,
                    residue_chain TEXT
                );
                CREATE TABLE Functional_Group_Atoms (
                    virus_name TEXT,
                    pdb_id TEXT,
                    ligand TEXT,
                    chain TEXT,
                    functional_group TEXT,
                    atom_id INTEGER,
                    exact_atom TEXT,
                    atom_type TEXT
                );
                """
            )
            conn.execute(
                """
                INSERT INTO Ligand_SMILES_Graphs (
                    graph_id, source_smiles, canonical_smiles, isomeric_smiles,
                    smiles_hash, atom_count, heavy_atom_count, bond_count,
                    formal_charge, rdkit_valid, parse_status, graph_method,
                    graph_version, rdkit_version, generated_at
                )
                VALUES (1, 'CC(C)N', 'CC(C)N', 'CC(C)N', 'hash1', 4, 4, 3, 0, 1, 'parsed', 'test', 'v1', 'test', '2026-01-01T00:00:00+00:00')
                """
            )
            atoms = [
                (1, 0, "C", 6, 0, 0, 0, 0, "SP3", "CHI_UNSPECIFIED", 3, 4, 0, 0, 0),
                (1, 1, "C", 6, 0, 0, 0, 0, "SP3", "CHI_UNSPECIFIED", 1, 4, 0, 0, 0),
                (1, 2, "C", 6, 0, 0, 0, 0, "SP3", "CHI_UNSPECIFIED", 1, 4, 0, 0, 0),
                (1, 3, "N", 7, 0, 0, 0, 0, "SP3", "CHI_UNSPECIFIED", 1, 3, 0, 1, 0),
            ]
            conn.executemany(
                """
                INSERT INTO Ligand_SMILES_Atoms (
                    graph_id, smiles_atom_index, element, atomic_number,
                    formal_charge, isotope, is_aromatic, is_in_ring,
                    hybridization, chiral_tag, degree, total_valence,
                    explicit_h_count, implicit_h_count, atom_map_number
                )
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                atoms,
            )
            conn.executemany(
                """
                INSERT INTO Ligand_SMILES_Bonds (
                    graph_id, smiles_bond_index, begin_atom_index, end_atom_index,
                    bond_type, bond_order
                )
                VALUES (1, ?, ?, ?, 'SINGLE', 1.0)
                """,
                [(0, 0, 1), (1, 0, 2), (2, 0, 3)],
            )
            conn.executemany(
                """
                INSERT INTO Ligand_Atoms_Smiles (
                    virus_name, pdb_id, ligand, chain, ligand_id, smiles
                )
                VALUES ('Virus', ?, ?, ?, ?, 'CC(C)N')
                """,
                [("T001", "LIG", "A", 1), ("T002", "BAD", "A", 1)],
            )
            conn.executemany(
                """
                INSERT INTO Ligand_SMILES_Graph_Assignments (
                    graph_id, pdb_code, model_id, ligand_chain,
                    ligand_residue_id, ligand_resname, ligand_insertion_code,
                    smiles_source_table, smiles_source_row_id,
                    smiles_source_column, source_smiles_hash, assignment_status,
                    mapping_status, pdb_to_smiles_mapped_atom_count,
                    pdb_ligand_atom_count, pdb_ligand_heavy_atom_count,
                    graph_method, graph_version, generated_at
                )
                VALUES (1, ?, 0, 'A', 1, ?, '', 'Ligand_Atoms_Smiles', ?, 'smiles', 'hash1', 'assigned', ?, ?, ?, ?, 'test', 'v1', '2026-01-01T00:00:00+00:00')
                """,
                [
                    ("T001", "LIG", 1, "complete", 4, 4, 4),
                    ("T002", "BAD", 2, "no_pdb_to_smiles_mapping", 0, 0, 0),
                ],
            )
            atom_rows = [
                (100, "CEN", "C", 0.0, 0.0, 0.0, 0, 0.0, False),
                (101, "LEFT", "C", -2.2, 0.0, 0.0, 1, 24.0, True),
                (102, "RGHT", "C", 2.2, 0.0, 0.0, 2, 22.0, True),
                (103, "NSTR", "N", 0.0, 2.2, 0.0, 3, 12.0, True),
            ]
            for atom_id, name, element, x, y, z, smiles_idx, sasa, exposed in atom_rows:
                conn.execute(
                    """
                    INSERT INTO ligand_atoms (
                        virus_name, pdb_id, ligand, chain, atom_id, exact_atom,
                        atom_type, x, y, z
                    )
                    VALUES ('Virus', 'T001', 'LIG', 'A', ?, ?, ?, ?, ?, ?)
                    """,
                    (atom_id, name, element, x, y, z),
                )
                conn.execute(
                    """
                    INSERT INTO SMILES_MAP_PDB (
                        virus_name, pdb_id, ligand, chain, exact_atom,
                        atom_id, atom_index, smiles_atom_index
                    )
                    VALUES ('Virus', 'T001', 'LIG', 'A', ?, ?, ?, ?)
                    """,
                    (name, atom_id, smiles_idx, smiles_idx),
                )
                conn.execute(
                    """
                    INSERT INTO RUPLEY_SASA_DATA (
                        virus_name, pdb_id, ligand, chain, exact_atom, atom_id,
                        SASA_Area
                    )
                    VALUES ('Virus', 'T001', 'LIG', 'A', ?, ?, ?)
                    """,
                    (name, atom_id, sasa),
                )
                if exposed:
                    conn.execute(
                        """
                        INSERT INTO solvent_exposed_atoms (
                            virus_name, pdb_id, ligand, chain, atom_id,
                            exact_atom, atom_type, x, y, z
                        )
                        VALUES ('Virus', 'T001', 'LIG', 'A', ?, ?, ?, ?, ?, ?)
                        """,
                        (atom_id, name, element, x, y, z),
                    )
            conn.execute(
                """
                INSERT INTO Arpeggio_Contacts_Data (
                    virus_name, pdb_id, ligand, ligand_id, chain, Contact,
                    Distance, exact_atom, atom_id, residue, residue_number,
                    residue_atom, residue_chain
                )
                VALUES ('Virus', 'T001', 'LIG', 1, 'A', 'polar', 2.8, 'NSTR', 103, 'ASP', 10, 'OD1', 'B')
                """
            )
            conn.execute(
                """
                INSERT INTO Functional_Group_Atoms (
                    virus_name, pdb_id, ligand, chain, functional_group,
                    atom_id, exact_atom, atom_type
                )
                VALUES ('Virus', 'T001', 'LIG', 'A', 'amine', 103, 'NSTR', 'N')
                """
            )
            conn.commit()

    def _run_builder(self, *args):
        return subprocess.run(
            [sys.executable, str(BUILDER), "--database", str(self.db_path), *args],
            cwd=ROOT,
            text=True,
            capture_output=True,
            check=False,
        )

    def test_migration_upgrade_and_rollback(self):
        up = subprocess.run(
            [sys.executable, str(MIGRATION), "--database", str(self.db_path), "--upgrade"],
            cwd=ROOT,
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertEqual(up.returncode, 0, up.stdout + up.stderr)
        with sqlite3.connect(self.db_path) as conn:
            self.assertIsNotNone(
                conn.execute(
                    "SELECT 1 FROM sqlite_master WHERE type='table' AND name='protacability_attachment_analysis'"
                ).fetchone()
            )
        down = subprocess.run(
            [sys.executable, str(MIGRATION), "--database", str(self.db_path), "--downgrade"],
            cwd=ROOT,
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertEqual(down.returncode, 0, down.stdout + down.stderr)
        with sqlite3.connect(self.db_path) as conn:
            self.assertIsNone(
                conn.execute(
                    "SELECT 1 FROM sqlite_master WHERE type='table' AND name='protacability_attachment_analysis'"
                ).fetchone()
            )

    def test_builder_dry_run_does_not_create_tables(self):
        result = self._run_builder("--pdb-code", "T001", "--ligand-resname", "LIG", "--limit", "1")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        with sqlite3.connect(self.db_path) as conn:
            self.assertIsNone(
                conn.execute(
                    "SELECT 1 FROM sqlite_master WHERE type='table' AND name='protacability_attachment_analysis'"
                ).fetchone()
            )

    def test_write_regions_are_deterministic_and_rerunnable(self):
        result = self._run_builder(
            "--write",
            "--replace-version",
            "attachment_v1",
            "--limit",
            "1",
            "--routine-density",
            "160",
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        with sqlite3.connect(self.db_path) as conn:
            conn.row_factory = sqlite3.Row
            analysis = get_analysis(
                conn,
                pdb_code="T001",
                ligand_resname="LIG",
                ligand_chain="A",
                ligand_residue_id=1,
            )
            self.assertIsNotNone(analysis)
            self.assertEqual(analysis["analysis_status"], "completed")
            summary = get_attachment_summary(
                conn,
                pdb_code="T001",
                ligand_resname="LIG",
                ligand_chain="A",
                ligand_residue_id=1,
            )
            self.assertEqual(summary["attachment_method_version"], "attachment_v1_1")
            self.assertEqual(summary["instance_resolution_status"], "resolved")
            regions = get_attachment_regions(conn, analysis["analysis_id"])
            self.assertGreaterEqual(len(regions), 2)
            atom_rows = get_attachment_atoms(conn, analysis["analysis_id"])
            left = next(row for row in atom_rows if row["pdb_atom_name"] == "LEFT")
            right = next(row for row in atom_rows if row["pdb_atom_name"] == "RGHT")
            self.assertNotEqual(left["region_id"], right["region_id"])
            nstr = next(row for row in atom_rows if row["pdb_atom_name"] == "NSTR")
            self.assertLess(nstr["attachment_score"], left["attachment_score"])
            self.assertIn("strong_contacts", nstr["cautions_json"])
            json.loads(left["reasons_json"])
            counts = {
                table: conn.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
                for table in (
                    "protacability_attachment_analysis",
                    "protacability_attachment_atoms",
                    "protacability_attachment_regions",
                )
            }
        rerun = self._run_builder("--write", "--limit", "1", "--routine-density", "160")
        self.assertEqual(rerun.returncode, 0, rerun.stdout + rerun.stderr)
        with sqlite3.connect(self.db_path) as conn:
            self.assertEqual(
                conn.execute("SELECT COUNT(*) FROM protacability_attachment_analysis").fetchone()[0],
                counts["protacability_attachment_analysis"],
            )
            self.assertEqual(
                conn.execute("SELECT COUNT(*) FROM protacability_attachment_atoms").fetchone()[0],
                counts["protacability_attachment_atoms"],
            )
            self.assertEqual(
                conn.execute("SELECT COUNT(*) FROM protacability_attachment_regions").fetchone()[0],
                counts["protacability_attachment_regions"],
            )

    def test_batch_failure_isolation(self):
        result = self._run_builder("--limit", "2", "--routine-density", "80")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        report = json.loads(result.stdout)
        self.assertEqual(report["instances_completed"], 2)
        self.assertEqual(report["instances_failed"], 0)
        self.assertEqual(report["analysis_status_counts"]["completed"], 1)
        self.assertEqual(report["analysis_status_counts"]["skipped"], 1)
        self.assertEqual(report["eligibility_status_counts"]["no_mapping"], 1)

    def test_zero_isolated_sasa_and_missing_sasa_are_not_positive_evidence(self):
        atoms = [
            {
                "pdb_atom_serial": 1,
                "complex_sasa": 5.0,
                "isolated_ligand_sasa": 0.0,
                "surface_point_count": 0,
                "existing_exposed_atom_status": 1,
                "degree": 1,
                "is_aromatic": 0,
                "is_in_ring": 0,
            },
            {
                "pdb_atom_serial": 2,
                "complex_sasa": None,
                "isolated_ligand_sasa": 30.0,
                "surface_point_count": 20,
                "existing_exposed_atom_status": 0,
                "degree": 1,
                "is_aromatic": 1,
                "is_in_ring": 1,
            },
        ]
        score_atoms(
            atoms,
            {},
            {},
            {"surface_density": 80, "minimum_surface_sasa": 2.0},
        )
        self.assertIn("zero_isolated_ligand_sasa", atoms[0]["cautions"])
        self.assertEqual(atoms[0]["candidate_attachment_flag"], 0)
        self.assertIn("missing_complex_sasa", atoms[1]["cautions"])
        self.assertEqual(atoms[1]["surface_defining_flag"], 0)


if __name__ == "__main__":
    unittest.main()
