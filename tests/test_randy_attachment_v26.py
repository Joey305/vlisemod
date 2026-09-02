"""Regression coverage for the Randy v2.6 attachment-site API contract."""

import os
import sqlite3
import unittest
from pathlib import Path

from flask import Flask

try:
    from RANDY import vlismod_data_routes as routes
except ModuleNotFoundError:  # deployed beside Randy's service module
    import vlismod_data_routes as routes


ROOT = Path(__file__).resolve().parents[1]
DATABASE = Path(os.environ.get("VLISMOD_DB_PATH", ROOT / "viral_data.db"))


@unittest.skipUnless(DATABASE.exists(), "requires the checked-in V-LiSEMOD database")
class RandyAttachmentV26Tests(unittest.TestCase):
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

    def test_current_database_exposes_v26_views_and_3tkg_roc_records(self):
        with sqlite3.connect(DATABASE) as conn:
            names = {
                row[0]
                for row in conn.execute(
                    "SELECT name FROM sqlite_master WHERE name IN (?, ?)",
                    ("v2_attachment_site_summary", "v2_attachment_site_candidates"),
                )
            }
            self.assertEqual(names, {"v2_attachment_site_summary", "v2_attachment_site_candidates"})
            rows = conn.execute(
                """
                SELECT ligand_instance_id, mapped_atom_count, exposed_mapped_atom_count,
                       candidate_attachment_atom_count,
                       chemically_supported_candidate_count, top_attachment_exact_atom,
                       top_attachment_site_score
                FROM protacability_attachment_site_summary
                WHERE pdb_code='3TKG' AND ligand_resname='ROC'
                ORDER BY ligand_instance_id
                """
            ).fetchall()
        self.assertEqual(rows, [
            (74791, 49, 18, 18, 4, "N1", 79.0),
            (74793, 49, 17, 17, 4, "N1", 79.0),
        ])

    def test_structure_detail_returns_current_atom_level_attachment_evidence(self):
        response = self.client.get(
            "/api/vlismod/protacability/structure-detail/3TKG?collapse_labels=1",
            headers={"Authorization": "Bearer test-token"},
        )
        self.assertEqual(response.status_code, 200)
        attachment = response.get_json()["attachment_sites"]
        self.assertTrue(attachment["data_available"])
        self.assertEqual(attachment["ligand_instance_id"], 74793)
        self.assertEqual(attachment["summary"]["mapped_atom_count"], 49)
        self.assertEqual(attachment["summary"]["attachment_exposed_mapped_atom_count"], 17)
        self.assertEqual(attachment["summary"]["attachment_candidate_atom_count"], 17)
        self.assertEqual(attachment["summary"]["best_attachment_confidence"], "Moderate")
        self.assertEqual(len(attachment["atoms"]), 17)

    def test_structure_detail_honors_exact_ligand_occurrence_identity(self):
        response = self.client.get(
            "/api/vlismod/protacability/structure-detail/3TKG?collapse_labels=1&ligand_instance_id=74791",
            headers={"Authorization": "Bearer test-token"},
        )
        self.assertEqual(response.status_code, 200)
        body = response.get_json()
        self.assertEqual(body["representative_ligand"]["ligand_instance_id"], 74791)
        self.assertEqual(body["attachment_sites"]["ligand_instance_id"], 74791)
        self.assertEqual(body["attachment_sites"]["summary"]["attachment_candidate_atom_count"], 18)
        self.assertTrue(all(context.get("ligand_instance_id") for context in body["ligand_contexts"]))

    def test_search_summary_preserves_attachment_evidence(self):
        response = self.client.get(
            "/api/vlismod/protacability/search?pdb_code=3TKG&collapse_labels=1",
            headers={"Authorization": "Bearer test-token"},
        )
        self.assertEqual(response.status_code, 200)
        row = response.get_json()["rows"][0]
        self.assertEqual(row["has_attachment_site_evidence"], 1)
        self.assertEqual(row["attachment_candidate_atom_count"], 18)
        self.assertEqual(row["attachment_display_site_count"], 1)
        self.assertEqual(row["best_attachment_score"], 79.0)
        self.assertEqual(row["best_attachment_confidence"], "Moderate")

    def test_missing_compatibility_views_are_explicitly_unavailable(self):
        with sqlite3.connect(":memory:") as conn:
            payload = routes._attachment_detail_payload(conn, {})
        self.assertFalse(payload["data_available"])
        self.assertIn("compatibility views", payload["message"])
