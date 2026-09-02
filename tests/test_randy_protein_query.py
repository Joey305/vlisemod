"""Regression coverage for Randy-backed Protein Query filter options."""

import os
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

