"""Regression coverage for Randy-backed Protein Query filter options."""

import os
import sqlite3
import unittest
from pathlib import Path

from flask import Flask

from RANDY import vlismod_data_routes as routes


ROOT = Path(__file__).resolve().parents[1]
DATABASE = Path(os.environ.get("VLISMOD_DB_PATH", ROOT / "viral_data.db"))


@unittest.skipUnless(DATABASE.exists(), "requires the checked-in V-LiSEMOD database")
class RandyProteinQueryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.previous_db = os.environ.get("VLISMOD_DB_PATH")
        cls.previous_token = os.environ.get("VLISMOD_API_TOKEN")
        os.environ["VLISMOD_DB_PATH"] = str(DATABASE)
        os.environ["VLISMOD_API_TOKEN"] = "test-token"
        cls.app = Flask(__name__)
        cls.app.register_blueprint(routes.vlismod_data_bp)
        cls.client = cls.app.test_client()

    @classmethod
    def tearDownClass(cls):
        if cls.previous_db is None:
            os.environ.pop("VLISMOD_DB_PATH", None)
        else:
            os.environ["VLISMOD_DB_PATH"] = cls.previous_db
        if cls.previous_token is None:
            os.environ.pop("VLISMOD_API_TOKEN", None)
        else:
            os.environ["VLISMOD_API_TOKEN"] = cls.previous_token

    def test_filter_options_supports_initial_and_hiv_cascade(self):
        headers = {"Authorization": "Bearer test-token"}
        response = self.client.get("/api/vlismod/virus-proteins/filter-options", headers=headers)
        self.assertEqual(response.status_code, 200)
        self.assertIn("HIV_1", response.get_json()["virus_names"])

        response = self.client.get(
            "/api/vlismod/virus-proteins/filter-options?virus_name=HIV_1",
            headers=headers,
        )
        self.assertEqual(response.status_code, 200)
        payload = response.get_json()
        self.assertIn("protease", payload["protein_types"])
        self.assertTrue(payload["ligands"])

    def test_ligand_comparison_mapping_preserves_occurrence_identity(self):
        response = self.client.get(
            "/api/vlismod/pdb-mapping?ligand_code=DR7",
            headers={"Authorization": "Bearer test-token"},
        )
        self.assertEqual(response.status_code, 200)
        records = response.get_json()["pdb_mapping"].values()
        self.assertTrue(records)
        self.assertTrue(all(record.get("ligand_instance_id") for record in records))
        self.assertTrue(all(record.get("legacy_key") for record in records))

    def test_ligand_smiles_uses_the_v2_mapping_atom_order(self):
        response = self.client.get(
            "/api/vlismod/ligand-smiles?ligand_id=DR7",
            headers={"Authorization": "Bearer test-token"},
        )
        self.assertEqual(response.status_code, 200)
        payload = response.get_json()
        self.assertTrue(payload["smiles"])
        with sqlite3.connect(DATABASE) as conn:
            expected = conn.execute(
                "SELECT smiles FROM ligands WHERE component_id = ?", ("DR7",)
            ).fetchone()[0]
        self.assertEqual(payload["smiles"], expected)

    def test_ligand_indexer_pairs_are_occurrence_resolved(self):
        headers = {"Authorization": "Bearer test-token"}
        response = self.client.get("/api/vlismod/pdb-residues/by-ligand?ligand_code=DR7", headers=headers)
        self.assertEqual(response.status_code, 200)
        pairs = response.get_json()["pairs"]
        self.assertTrue(pairs)
        self.assertEqual(len({pair["ligand_instance_id"] for pair in pairs}), len(pairs))
        self.assertTrue(all(pair["ligand_instance_id"] for pair in pairs))
        self.assertTrue(all(pair.get("ligand_id") for pair in pairs))
        exact = next(pair for pair in pairs if pair["pdb_id"] == "3EM4" and pair["chain"] == "V" and str(pair["ligand_id"]) == "100")
        self.assertEqual(exact["ligand_instance_id"], 60443)

    def test_ligand_indexer_interactions_honor_selected_occurrence(self):
        headers = {"Authorization": "Bearer test-token"}
        response = self.client.get(
            "/api/vlismod/interaction-records?pdb_id=3EM4&ligand=DR7&ligand_id=100&chain=V&ligand_instance_id=60443",
            headers=headers,
        )
        self.assertEqual(response.status_code, 200)
        records = response.get_json()["records"]
        self.assertTrue(records)
        self.assertEqual({record["ligand_instance_id"] for record in records}, {60443})
