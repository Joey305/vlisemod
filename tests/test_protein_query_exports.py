import io
import importlib.util
import json
import sqlite3
import sys
import tempfile
import unittest
import zipfile
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
APP_SPEC = importlib.util.spec_from_file_location("vlismod_web_app", ROOT / "app.py")
vlismod_app = importlib.util.module_from_spec(APP_SPEC)
sys.modules[APP_SPEC.name] = vlismod_app
APP_SPEC.loader.exec_module(vlismod_app)


class ProteinQueryExportTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.db_path = Path(self.tmp.name) / "protein-query.db"
        self._create_database()
        self.original_db_path = vlismod_app.LOCAL_DB_PATH
        vlismod_app.LOCAL_DB_PATH = self.db_path
        vlismod_app.app.config.update(TESTING=True)
        self.client = vlismod_app.app.test_client()

    def tearDown(self):
        vlismod_app.LOCAL_DB_PATH = self.original_db_path
        self.tmp.cleanup()

    def _create_database(self):
        with sqlite3.connect(self.db_path) as conn:
            conn.executescript(
                """
                CREATE TABLE Virus_Proteins (virus_name TEXT, pdb_id TEXT, protein TEXT);
                CREATE TABLE structures (structure_id INTEGER PRIMARY KEY, entry_id TEXT UNIQUE);
                CREATE TABLE ligand_instances (
                    ligand_instance_id INTEGER PRIMARY KEY,
                    structure_id INTEGER,
                    label_comp_id TEXT,
                    curation_status TEXT
                );
                CREATE TABLE ligand_instance_atoms (
                    ligand_instance_id INTEGER,
                    atom_site_id INTEGER,
                    selected_conformer INTEGER
                );
                CREATE TABLE ligand_atoms (pdb_id TEXT, ligand TEXT, atom_id INTEGER);
                CREATE TABLE Ligand_Synonyms (ligand TEXT, synonym TEXT);
                CREATE TABLE Ligand_Arp_Diagram (pdb_id TEXT, ligand TEXT);
                """
            )
            conn.executemany(
                "INSERT INTO Virus_Proteins VALUES (?, ?, ?)",
                [
                    ("HIV_1", "GOOD1", "integrase"),
                    ("HIV_1", "GOOD2", "integrase"),
                    ("HIV_1", "EMPTY1", "integrase"),
                    ("HPV_16", "HPVGOOD", "capsid_protein"),
                ],
            )
            conn.executemany(
                "INSERT INTO structures VALUES (?, ?)",
                [(1, "GOOD1"), (2, "GOOD2"), (3, "EMPTY1"), (4, "HPVGOOD")],
            )
            conn.executemany(
                "INSERT INTO ligand_instances VALUES (?, ?, ?, ?)",
                [
                    (101, 1, "LIG", "included"),
                    (102, 2, "LIG", "included"),
                    (103, 3, "LIG", "excluded"),
                    (104, 4, "HPV", "included"),
                ],
            )
            conn.executemany(
                "INSERT INTO ligand_instance_atoms VALUES (?, ?, ?)",
                [(101, 1, 1), (102, 2, 1), (103, 3, 1), (104, 4, 1)],
            )
            conn.executemany(
                "INSERT INTO ligand_atoms VALUES (?, ?, ?)",
                [("GOOD1", "LIG", 1), ("GOOD2", "LIG", 2)],
            )

    def _export(self, pdb_codes):
        return self.client.post(
            "/export_data_to_excel",
            json={"pdb_codes": pdb_codes, "data_sets": ["Ligand Atoms"]},
        )

    def test_protein_query_excludes_nonexportable_pdbs(self):
        response = self.client.post(
            "/get_pdbs_for_virus_protein",
            json={"virus_name": "HIV_1", "protein_types": ["integrase"]},
        )
        self.assertEqual(response.status_code, 200)
        self.assertNotIn("EMPTY1", response.get_json()["pdb_codes"])

    def test_protein_query_includes_exportable_pdbs(self):
        with sqlite3.connect(self.db_path) as conn:
            pdb_codes = vlismod_app.get_exportable_protein_query_pdbs(
                conn, "HIV_1", ["integrase"]
            )
        self.assertEqual(pdb_codes, ["GOOD1", "GOOD2"])

    def test_filter_options_cascade_from_virus_to_protein_to_ligand(self):
        viruses = self.client.get("/get_protein_query_filter_options").get_json()
        self.assertEqual(viruses["virus_names"], ["HIV_1", "HPV_16"])

        proteins = self.client.get(
            "/get_protein_query_filter_options?virus_name=HPV_16"
        ).get_json()
        self.assertEqual(proteins["protein_types"], ["capsid_protein"])

        ligands = self.client.get(
            "/get_protein_query_filter_options?virus_name=HIV_1&protein_type=integrase"
        ).get_json()
        self.assertEqual(ligands["ligands"], [{"ligand_code": "LIG", "synonyms": []}])

    def test_empty_export_returns_no_data_error(self):
        response = self._export(["EMPTY1"])
        self.assertEqual(response.status_code, 422)
        payload = response.get_json()
        self.assertEqual(payload["error"], "no_exportable_data")
        self.assertEqual(payload["non_exportable_pdb_codes"], ["EMPTY1"])

    def test_valid_export_contains_rows(self):
        response = self._export(["GOOD1"])
        self.assertEqual(response.status_code, 200)
        with zipfile.ZipFile(io.BytesIO(response.data)) as archive:
            self.assertIn("Ligand_Atoms.csv", archive.namelist())
            csv_text = archive.read("Ligand_Atoms.csv").decode()
            self.assertIn("GOOD1", csv_text)
            self.assertNotIn("GOOD2", csv_text)

    def test_valid_export_workbook_not_blank(self):
        response = self._export(["GOOD1"])
        with zipfile.ZipFile(io.BytesIO(response.data)) as archive:
            workbook = pd.ExcelFile(io.BytesIO(archive.read("combined_data.xlsx")))
            self.assertEqual(workbook.sheet_names, ["Ligand Atoms"])
            self.assertEqual(len(pd.read_excel(workbook, sheet_name="Ligand Atoms")), 1)

    def test_mixed_export_filters_nonexportable_pdb_and_records_it(self):
        response = self._export(["GOOD1", "EMPTY1"])
        self.assertEqual(response.status_code, 200)
        with zipfile.ZipFile(io.BytesIO(response.data)) as archive:
            summary = json.loads(archive.read("export_summary.json"))
            self.assertEqual(summary["exportable_pdb_codes"], ["GOOD1"])
            self.assertEqual(summary["non_exportable_pdb_codes"], ["EMPTY1"])


if __name__ == "__main__":
    unittest.main()
