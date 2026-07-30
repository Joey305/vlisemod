import sqlite3
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT / "TOOLS"
BUILDER = TOOLS / "build_ligand_smiles_graphs.py"
MIGRATION = TOOLS / "migrate_ligand_smiles_graphs.py"

sys.path.insert(0, str(TOOLS))

from ligand_smiles_graph_queries import (  # noqa: E402
    atoms_within_n_bonds,
    connected_components,
    get_neighbors,
    map_pdb_atom_to_graph_atom,
    ring_atoms,
    shortest_graph_distance,
    terminal_atoms,
)


def rdkit_available():
    try:
        import rdkit  # noqa: F401
        return True
    except Exception:
        return False


class LigandSmilesGraphTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.db_path = Path(self.tmp.name) / "sample.db"
        self._create_source_database()

    def tearDown(self):
        self.tmp.cleanup()

    def _create_source_database(self):
        with sqlite3.connect(self.db_path) as conn:
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
                CREATE TABLE Functional_GROUPED (
                    virus_name TEXT,
                    pdb_id TEXT,
                    ligand TEXT,
                    smiles TEXT,
                    functional_groups TEXT
                );
                CREATE TABLE protacability_warhead_linkability (
                    pdb_code TEXT,
                    model_id INTEGER,
                    ligand_resname TEXT,
                    ligand_chain TEXT,
                    ligand_residue_id INTEGER,
                    ligand_insertion_code TEXT,
                    representative_smiles TEXT
                );
                CREATE TABLE SMILES_MAP_PDB (
                    pdb_id TEXT,
                    ligand TEXT,
                    chain TEXT,
                    exact_atom TEXT,
                    atom_id INTEGER,
                    atom_index INTEGER,
                    smiles_atom_index INTEGER
                );
                CREATE TABLE ligand_atoms (
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
                """
            )
            rows = [
                ("Virus", "T001", "ALP", "A", 1, "CCO"),
                ("Virus", "T002", "ARO", "A", 1, "c1ccccc1"),
                ("Virus", "T003", "CHG", "A", 1, "C[NH3+]"),
                ("Virus", "T004", "STC", "A", 1, "F[C@H](Cl)Br"),
                ("Virus", "T005", "FRG", "A", 1, "CC.O"),
                ("Virus", "T006", "SLT", "A", 1, "C[NH3+].[Cl-]"),
                ("Virus", "T007", "BAD", "A", 1, "C1CC"),
                ("Virus", "T008", "DUP", "A", 1, "CCO"),
                ("Virus", "T009", "ORD", "A", 1, "OCC"),
            ]
            conn.executemany(
                """
                INSERT INTO Ligand_Atoms_Smiles (
                    virus_name, pdb_id, ligand, chain, ligand_id, smiles
                )
                VALUES (?, ?, ?, ?, ?, ?)
                """,
                rows,
            )
            conn.execute(
                """
                INSERT INTO Functional_GROUPED (virus_name, pdb_id, ligand, smiles)
                VALUES ('Virus', 'T010', 'FG1', 'CC(=O)O')
                """
            )
            conn.execute(
                """
                INSERT INTO protacability_warhead_linkability (
                    pdb_code, model_id, ligand_resname, ligand_chain,
                    ligand_residue_id, representative_smiles
                )
                VALUES ('T011', 0, 'PRT', 'A', 1, 'N#CC')
                """
            )
            for idx, name in enumerate(["C1", "C2", "O1"]):
                conn.execute(
                    """
                    INSERT INTO SMILES_MAP_PDB (
                        pdb_id, ligand, chain, exact_atom, atom_id, atom_index,
                        smiles_atom_index
                    )
                    VALUES ('T001', 'ALP', 'A', ?, ?, ?, ?)
                    """,
                    (name, 100 + idx, idx, idx),
                )
                conn.execute(
                    """
                    INSERT INTO ligand_atoms (
                        pdb_id, ligand, chain, atom_id, exact_atom, atom_type, x, y, z
                    )
                    VALUES ('T001', 'ALP', 'A', ?, ?, ?, ?, 0, 0)
                    """,
                    (100 + idx, name, name[0], float(idx)),
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

    @unittest.skipUnless(rdkit_available(), "RDKit is required for graph generation tests")
    def test_builder_writes_graphs_and_is_rerunnable(self):
        result = self._run_builder("--write")
        self.assertEqual(result.returncode, 2, result.stdout + result.stderr)
        with sqlite3.connect(self.db_path) as conn:
            graph_count = conn.execute("SELECT COUNT(*) FROM Ligand_SMILES_Graphs").fetchone()[0]
            assignment_count = conn.execute("SELECT COUNT(*) FROM Ligand_SMILES_Graph_Assignments").fetchone()[0]
            atom_count = conn.execute("SELECT COUNT(*) FROM Ligand_SMILES_Atoms").fetchone()[0]
            bond_count = conn.execute("SELECT COUNT(*) FROM Ligand_SMILES_Bonds").fetchone()[0]
        self.assertEqual(graph_count, 10)
        self.assertEqual(assignment_count, 11)
        self.assertGreater(atom_count, 0)
        self.assertGreater(bond_count, 0)

        rerun = self._run_builder("--write")
        self.assertEqual(rerun.returncode, 2, rerun.stdout + rerun.stderr)
        with sqlite3.connect(self.db_path) as conn:
            self.assertEqual(conn.execute("SELECT COUNT(*) FROM Ligand_SMILES_Graphs").fetchone()[0], graph_count)
            self.assertEqual(conn.execute("SELECT COUNT(*) FROM Ligand_SMILES_Graph_Assignments").fetchone()[0], assignment_count)

    @unittest.skipUnless(rdkit_available(), "RDKit is required for graph generation tests")
    def test_atom_index_stability_bonds_rings_and_helpers(self):
        result = self._run_builder("--write")
        self.assertEqual(result.returncode, 2, result.stdout + result.stderr)
        with sqlite3.connect(self.db_path) as conn:
            conn.row_factory = sqlite3.Row
            graph_id = conn.execute(
                """
                SELECT graph_id
                FROM Ligand_SMILES_Graph_Assignments
                WHERE pdb_code='T001' AND ligand_resname='ALP'
                """
            ).fetchone()["graph_id"]
            neighbors = get_neighbors(conn, graph_id, 1)
            self.assertEqual({row["neighbor_atom_index"] for row in neighbors}, {0, 2})
            self.assertEqual(shortest_graph_distance(conn, graph_id, 0, 2), 2)
            self.assertEqual(atoms_within_n_bonds(conn, graph_id, 0, 1), [0, 1])
            self.assertEqual(connected_components(conn, graph_id), [[0, 1, 2]])
            self.assertIn(0, terminal_atoms(conn, graph_id))
            mapped = map_pdb_atom_to_graph_atom(
                conn,
                pdb_code="T001",
                ligand_resname="ALP",
                ligand_chain="A",
                ligand_residue_id=1,
                pdb_atom_id=100,
            )
            self.assertEqual(mapped["smiles_atom_index"], 0)

            aromatic_graph = conn.execute(
                """
                SELECT graph_id
                FROM Ligand_SMILES_Graph_Assignments
                WHERE pdb_code='T002' AND ligand_resname='ARO'
                """
            ).fetchone()["graph_id"]
            self.assertEqual(len(ring_atoms(conn, aromatic_graph)), 6)
            aromatic_bonds = conn.execute(
                "SELECT COUNT(*) FROM Ligand_SMILES_Bonds WHERE graph_id=? AND is_aromatic=1",
                (aromatic_graph,),
            ).fetchone()[0]
            self.assertEqual(aromatic_bonds, 6)

    @unittest.skipUnless(rdkit_available(), "RDKit is required for graph generation tests")
    def test_duplicate_exact_smiles_share_graph_but_atom_order_variants_do_not(self):
        result = self._run_builder("--write")
        self.assertEqual(result.returncode, 2, result.stdout + result.stderr)
        with sqlite3.connect(self.db_path) as conn:
            conn.row_factory = sqlite3.Row
            graph_alp = conn.execute(
                "SELECT graph_id FROM Ligand_SMILES_Graph_Assignments WHERE ligand_resname='ALP'"
            ).fetchone()["graph_id"]
            graph_dup = conn.execute(
                "SELECT graph_id FROM Ligand_SMILES_Graph_Assignments WHERE ligand_resname='DUP'"
            ).fetchone()["graph_id"]
            graph_ord = conn.execute(
                "SELECT graph_id FROM Ligand_SMILES_Graph_Assignments WHERE ligand_resname='ORD'"
            ).fetchone()["graph_id"]
            self.assertEqual(graph_alp, graph_dup)
            self.assertNotEqual(graph_alp, graph_ord)

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
                    "SELECT 1 FROM sqlite_master WHERE type='table' AND name='Ligand_SMILES_Graphs'"
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
                    "SELECT 1 FROM sqlite_master WHERE type='table' AND name='Ligand_SMILES_Graphs'"
                ).fetchone()
            )


if __name__ == "__main__":
    unittest.main()

