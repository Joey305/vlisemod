import importlib.util
import json
import os
import sqlite3
import sys
import tempfile
import time
import unittest
from pathlib import Path
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
APP_SPEC = importlib.util.spec_from_file_location("vlismod_ligand_comparison_app", ROOT / "app.py")
vlismod_app = importlib.util.module_from_spec(APP_SPEC)
sys.modules[APP_SPEC.name] = vlismod_app
APP_SPEC.loader.exec_module(vlismod_app)


class LigandComparisonTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.db_path = Path(self.tmp.name) / "ligand-comparison.db"
        self._create_database()
        self.original_db_path = vlismod_app.LOCAL_DB_PATH
        vlismod_app.LOCAL_DB_PATH = self.db_path
        vlismod_app.app.config.update(TESTING=True)
        self.client = vlismod_app.app.test_client()

    def tearDown(self):
        vlismod_app.LOCAL_DB_PATH = self.original_db_path
        self.tmp.cleanup()

    def _create_database(self):
        partner = json.dumps(
            {"label_comp_id": "PRO", "auth_seq_id": "42", "auth_atom_id": "CA", "auth_asym_id": "P"}
        )
        with sqlite3.connect(self.db_path) as conn:
            conn.executescript(
                """
                CREATE TABLE structures (structure_id INTEGER PRIMARY KEY, entry_id TEXT UNIQUE);
                CREATE TABLE structure_classifications (structure_id INTEGER, virus_label TEXT);
                CREATE TABLE ligand_instances (
                    ligand_instance_id INTEGER PRIMARY KEY, structure_id INTEGER, ligand_id INTEGER,
                    deposited_model_num TEXT, label_asym_id TEXT, label_comp_id TEXT,
                    auth_asym_id TEXT, auth_seq_id TEXT, insertion_code_normalized TEXT,
                    curation_status TEXT
                );
                CREATE TABLE ligand_instance_atoms (
                    ligand_instance_atom_id INTEGER PRIMARY KEY, ligand_instance_id INTEGER,
                    atom_site_id TEXT, label_atom_id TEXT, auth_atom_id TEXT
                );
                CREATE TABLE ligand_smiles_atom_mapping (
                    mapping_id INTEGER PRIMARY KEY, run_id INTEGER, ligand_instance_id INTEGER,
                    ligand_instance_atom_id INTEGER, smiles_atom_index INTEGER,
                    mapping_status TEXT, method_version TEXT
                );
                CREATE TABLE ligand_sasa_atoms (
                    sasa_id INTEGER PRIMARY KEY, run_id INTEGER, ligand_instance_id INTEGER,
                    ligand_instance_atom_id INTEGER, sasa_area REAL, legacy_exposed INTEGER,
                    status TEXT, method_version TEXT
                );
                CREATE TABLE ligand_arpeggio_runs (
                    arpeggio_run_id INTEGER PRIMARY KEY, run_id INTEGER,
                    ligand_instance_id INTEGER, status TEXT
                );
                CREATE TABLE arpeggio_raw_contact_labels (
                    raw_contact_id INTEGER PRIMARY KEY, run_id INTEGER,
                    ligand_instance_id INTEGER, raw_contact_index INTEGER,
                    interaction_label TEXT, distance REAL, bgn_json TEXT, end_json TEXT,
                    ligand_instance_atom_id INTEGER, partner_identity_json TEXT,
                    filter_class TEXT
                );
                CREATE TABLE Ligand_Synonyms (ligand TEXT, synonym TEXT);
                CREATE TABLE Ligand_Atoms_Smiles (
                    ligand TEXT, pdb_id TEXT, smiles TEXT, molecular_weight REAL
                );
                CREATE TABLE ligands (
                    ligand_id INTEGER PRIMARY KEY, component_id TEXT, smiles TEXT
                );
                """
            )
            conn.executemany("INSERT INTO structures VALUES (?, ?)", [(1, "3EKY"), (2, "4PHV")])
            conn.executemany(
                "INSERT INTO structure_classifications VALUES (?, ?)", [(1, "HIV_1"), (2, "HIV_1")]
            )
            conn.executemany(
                """INSERT INTO ligand_instances VALUES (?, ?, 1, '1', 'A', 'DR7', 'A', '100', '', 'included')""",
                [(11, 1), (12, 2), (13, 2)],
            )
            conn.executemany(
                "INSERT INTO ligand_instance_atoms VALUES (?, ?, ?, ?, ?)",
                [(101, 11, "1", "C1", "C1"), (102, 11, "2", "N1", "N1"), (201, 12, "1", "C1", "C1"), (301, 13, "1", "O1", "O1")],
            )
            conn.executemany(
                """INSERT INTO ligand_smiles_atom_mapping VALUES (?, ?, ?, ?, ?, 'mapped', 'legacy_mcs_etkdg_uff_cif_v2.5')""",
                [(1, 501, 11, 101, 0), (2, 501, 11, 102, 1), (3, 502, 12, 201, 0), (4, 503, 13, 301, 0)],
            )
            conn.executemany(
                """INSERT INTO ligand_sasa_atoms VALUES (?, ?, ?, ?, 10.0, ?, 'complete', 'biopython-shrake_rupley-1.40-cif-v2.1')""",
                [(1, 601, 11, 101, 1), (2, 601, 11, 102, 0), (3, 602, 12, 201, 1), (4, 603, 13, 301, 0)],
            )
            conn.executemany(
                "INSERT INTO ligand_arpeggio_runs VALUES (?, ?, ?, 'completed')", [(1, 701, 11), (2, 702, 12)])
            conn.executemany(
                """INSERT INTO arpeggio_raw_contact_labels VALUES (?, ?, ?, ?, ?, ?, '{}', '{}', ?, ?, 'raw_environment')""",
                [(1, 701, 11, 1, "hbond", 2.8, 101, partner), (2, 701, 11, 2, "proximal", 4.0, 102, partner), (3, 702, 12, 1, "vdw", 3.5, 201, partner)],
            )
            conn.execute("INSERT INTO Ligand_Synonyms VALUES ('DR7', 'ATAZANAVIR')")
            conn.execute("INSERT INTO Ligand_Atoms_Smiles VALUES ('DR7', '3EKY', 'COC', 46.0)")
            conn.execute("INSERT INTO ligands VALUES (1, 'DR7', 'CCO')")

    def _payload(self, occurrence_ids=(11, 12)):
        return {"ligand": "DR7", "occurrence_ids": list(occurrence_ids), "pdb_ids": ["3EKY-100-A", "4PHV-100-A"]}

    def test_compare_ligands_page_loads(self):
        page = self.client.get("/compare_ligands")
        self.assertEqual(page.status_code, 200)
        html = page.get_data(as_text=True)
        self.assertIn('id="comparison-rail"', html)
        self.assertIn('id="comparison-rail-content"', html)
        self.assertNotIn('id="comparison-protacability-occurrences"', html)

    def test_compare_ligands_list_loads(self):
        response = self.client.get("/get_ligands_with_synonyms")
        self.assertEqual(response.status_code, 200)
        self.assertIn({"ligand_code": "DR7", "synonym": "ATAZANAVIR"}, response.get_json())

    def test_compare_ligands_smiles_svg_dr7(self):
        response = self.client.get("/get_smiles_svg/DR7")
        self.assertEqual(response.status_code, 200)
        self.assertIn("<svg", response.get_json()["svg"])

    def test_comparison_svg_uses_the_mapping_smiles_atom_order(self):
        payload = vlismod_app._local_get_smiles_payload("DR7")
        self.assertEqual(payload["smiles"], "CCO")

    def test_compare_ligands_occurrence_mapping_dr7(self):
        response = self.client.get("/get_pdb_mapping/DR7")
        self.assertEqual(response.status_code, 200)
        mappings = response.get_json()["pdb_mapping"]
        self.assertEqual(set(mappings), {"11", "12", "13"})
        self.assertEqual(mappings["11"]["legacy_key"], "3EKY-100-A")
        self.assertEqual(mappings["11"]["ligand_instance_id"], 11)

    def test_compare_ligands_contact_query_is_occurrence_scoped(self):
        response = self.client.post("/compare_ligand_interactions", json=self._payload((11,)))
        self.assertEqual(response.status_code, 200)
        rows = response.get_json()["interactions_data"]
        self.assertEqual([row["ligand_instance_id"] for row in rows], [11])
        self.assertEqual(rows[0]["mapped_atom_count"], 2)
        self.assertEqual(rows[0]["solvent_exposed_atom_count"], 1)
        self.assertEqual(len(rows[0]["interactions"]), 1)

    def test_compare_ligand_interactions_route_accepts_post(self):
        rule = next(rule for rule in vlismod_app.app.url_map.iter_rules() if rule.rule == "/compare_ligand_interactions")
        self.assertIn("POST", rule.methods)

    def test_compare_ligands_dr7_returns_within_reasonable_time(self):
        started = time.perf_counter()
        response = self.client.post("/compare_ligand_interactions", json=self._payload())
        elapsed = time.perf_counter() - started
        self.assertEqual(response.status_code, 200)
        self.assertLess(elapsed, 2.0)

    def test_compare_ligands_no_data_response(self):
        response = self.client.post("/compare_ligand_interactions", json=self._payload((13,)))
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.get_json()["status"], "no_data")
        self.assertEqual(response.get_json()["interactions_data"], [])

    def test_compare_ligands_randy_payload_is_proxied(self):
        expected = {"interactions_data": [{"pdb_id": "3EKY", "interactions": []}], "smiles_interactions_data": []}
        with patch.dict(os.environ, {"VLISMOD_DATA_BACKEND": "randy"}, clear=False), patch.object(
            vlismod_app, "randy_post", return_value=expected
        ) as randy_post:
            response = self.client.post("/compare_ligand_interactions", json=self._payload((11,)))
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.get_json(), expected)
        randy_post.assert_called_once_with(
            "ligand-interactions/compare",
            json={"ligand": "DR7", "ligand_instance_ids": [11]},
        )

    def test_compare_ligands_randy_failure_is_not_a_zero_result(self):
        with patch.dict(os.environ, {"VLISMOD_DATA_BACKEND": "randy"}, clear=False), patch.object(
            vlismod_app, "randy_post", side_effect=vlismod_app.RandyBackendError("upstream route unavailable", status_code=404)
        ):
            response = self.client.post("/compare_ligand_interactions", json=self._payload((11,)))
        self.assertEqual(response.status_code, 404)
        self.assertEqual(response.get_json()["error"], "upstream route unavailable")

    def test_comparison_frontend_keeps_retrieval_failures_distinct_from_zero_results(self):
        template = (ROOT / "templates" / "compare_ligands.html").read_text()
        self.assertIn("Unable to retrieve interaction data for the selected ligand occurrences.", template)
        self.assertIn("Interaction retrieval failed", template)
        self.assertIn("data.status === 'no_data'", template)

    def test_compare_ligands_synonym_uses_component_id(self):
        page = self.client.get("/compare_ligands").get_data(as_text=True)
        self.assertIn('value="${item.ligand_code}"', page)
        self.assertIn("ligand_instance_ids: occurrenceIds", page)


if __name__ == "__main__":
    unittest.main()
