import importlib.util
import json
import sqlite3
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
APP_SPEC = importlib.util.spec_from_file_location("vlismod_ligand_indexer_app", ROOT / "app.py")
vlismod_app = importlib.util.module_from_spec(APP_SPEC)
sys.modules[APP_SPEC.name] = vlismod_app
APP_SPEC.loader.exec_module(vlismod_app)


class LigandIndexerInteractionTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.db_path = Path(self.tmp.name) / "ligand-indexer.db"
        self.chart_dir = Path(self.tmp.name) / "charts"
        self._create_database()
        self.original_db_path = vlismod_app.LOCAL_DB_PATH
        self.original_charts_dir = vlismod_app.CHARTS_DIR
        vlismod_app.LOCAL_DB_PATH = self.db_path
        vlismod_app.CHARTS_DIR = str(self.chart_dir)
        vlismod_app.app.config.update(TESTING=True)
        self.client = vlismod_app.app.test_client()

    def tearDown(self):
        vlismod_app.LOCAL_DB_PATH = self.original_db_path
        vlismod_app.CHARTS_DIR = self.original_charts_dir
        self.tmp.cleanup()

    def _create_database(self):
        partner = json.dumps(
            {"label_comp_id": "PRO", "auth_seq_id": "42", "auth_atom_id": "CA", "auth_asym_id": "P"}
        )
        with sqlite3.connect(self.db_path) as conn:
            conn.executescript(
                """
                CREATE TABLE structures (structure_id INTEGER PRIMARY KEY, entry_id TEXT UNIQUE);
                CREATE TABLE ligand_instances (
                    ligand_instance_id INTEGER PRIMARY KEY, structure_id INTEGER,
                    label_comp_id TEXT, auth_asym_id TEXT, auth_seq_id TEXT,
                    insertion_code_normalized TEXT, deposited_model_num TEXT,
                    curation_status TEXT
                );
                CREATE TABLE ligand_instance_atoms (
                    ligand_instance_atom_id INTEGER PRIMARY KEY,
                    ligand_instance_id INTEGER, atom_site_id TEXT,
                    label_atom_id TEXT, auth_atom_id TEXT
                );
                CREATE TABLE ligand_arpeggio_runs (
                    run_id INTEGER PRIMARY KEY, ligand_instance_id INTEGER, status TEXT
                );
                CREATE TABLE arpeggio_raw_contact_labels (
                    raw_contact_id INTEGER PRIMARY KEY, run_id INTEGER,
                    ligand_instance_id INTEGER, raw_contact_index INTEGER,
                    interaction_label TEXT, distance REAL, bgn_json TEXT, end_json TEXT,
                    ligand_instance_atom_id INTEGER, partner_identity_json TEXT,
                    filter_class TEXT, partner_source_atom_site_id TEXT,
                    partner_mapping_status TEXT
                );
                CREATE TABLE Ligand_Synonyms (ligand TEXT, synonym TEXT);
                """
            )
            conn.execute("INSERT INTO structures VALUES (1, '3EL1')")
            conn.executemany(
                "INSERT INTO ligand_instances VALUES (?, 1, 'DR7', 'A', '100', '', '1', 'included')",
                [(11,), (12,), (13,)],
            )
            conn.executemany(
                "INSERT INTO ligand_instance_atoms VALUES (?, ?, ?, ?, ?)",
                [(101, 11, '1', 'C1', 'C1'), (102, 12, '2', 'N1', 'N1')],
            )
            conn.executemany(
                "INSERT INTO ligand_arpeggio_runs VALUES (?, ?, 'completed')",
                [(201, 11), (202, 12)],
            )
            conn.executemany(
                """
                INSERT INTO arpeggio_raw_contact_labels VALUES (?, ?, ?, ?, ?, ?, '{}', '{}', ?, ?, 'raw_environment', NULL, NULL)
                """,
                [
                    (301, 201, 11, 1, 'hbond', 2.8, 101, partner),
                    (302, 201, 11, 2, 'vdw', 3.6, 101, partner),
                    (303, 202, 12, 1, 'proximal', 4.0, 102, partner),
                ],
            )
            conn.execute("INSERT INTO Ligand_Synonyms VALUES ('DR7', 'ATAZANAVIR')")

    def _chart_payload(self, ligand_instance_id=11):
        return {
            'pdb_id': '3EL1', 'ligand': 'DR7', 'ligand_id': '100',
            'chain': 'A', 'ligand_instance_id': ligand_instance_id,
        }

    def test_ligand_indexer_page_loads(self):
        self.assertEqual(self.client.get('/ligand_indexer').status_code, 200)

    def test_ligand_indexer_ligand_options_load(self):
        response = self.client.get('/get_ligands_with_synonyms')
        self.assertEqual(response.status_code, 200)
        self.assertIn({'ligand_code': 'DR7', 'synonym': 'ATAZANAVIR'}, response.get_json())

    def test_ligand_indexer_occurrences_load_for_dr7(self):
        response = self.client.get('/get_pdb_residue_by_ligand/DR7')
        self.assertEqual(response.status_code, 200)
        pairs = response.get_json()['pairs']
        self.assertEqual({pair['ligand_instance_id'] for pair in pairs}, {11, 12})
        self.assertTrue(all(pair['pdb_id'] == '3EL1' for pair in pairs))

    def test_generate_charts_dr7_3el1_returns_data(self):
        response = self.client.post('/generate_charts', json=self._chart_payload())
        self.assertEqual(response.status_code, 200)
        payload = response.get_json()
        self.assertEqual(payload['interaction_label_count'], 2)
        self.assertEqual(payload['ligand_instance_id'], 11)
        self.assertEqual(len(payload['chart_urls']), 5)
        self.assertEqual(len(list(self.chart_dir.glob('*.png'))), 5)

    def test_generate_charts_no_data_returns_clear_response(self):
        response = self.client.post('/generate_charts', json=self._chart_payload(13))
        self.assertEqual(response.status_code, 404)
        self.assertEqual(response.get_json()['error'], 'no_interaction_data')

    def test_ligand_indexer_synonym_submits_component_id(self):
        page = self.client.get('/ligand_indexer').get_data(as_text=True)
        self.assertIn('value="${item.ligand_code}"', page)
        self.assertIn('ligand: selectedLigand', page)


if __name__ == '__main__':
    unittest.main()
