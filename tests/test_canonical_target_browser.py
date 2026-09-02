import importlib.util
import sqlite3
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATABASE = ROOT / "viral_data.db"
PASS10_AUTHORITY = ROOT / "outputs" / "protein_taxonomy_reconciliation_pass10" / "canonical_target_occurrences_pass10.csv"
APP_SPEC = importlib.util.spec_from_file_location("vlismod_canonical_target_app", ROOT / "app.py")
vlismod_app = importlib.util.module_from_spec(APP_SPEC)
sys.modules[APP_SPEC.name] = vlismod_app
APP_SPEC.loader.exec_module(vlismod_app)


@unittest.skipUnless(DATABASE.exists(), "live V-LiSEMOD database is required")
class CanonicalTargetBrowserTests(unittest.TestCase):
    def test_normalized_composite_virus_target_identity_has_no_synonym_fragmentation(self):
        db = sqlite3.connect(DATABASE)
        try:
            occurrence_count, distinct_target_ids, composite_groups = db.execute(
                """
                SELECT count(*), count(DISTINCT canonical_target_id),
                       count(DISTINCT virus_name || char(31) || canonical_target_id)
                FROM v2_target_browser_ligand_context
                """
            ).fetchone()
            group_count = db.execute("SELECT count(*) FROM v2_target_browser_groups").fetchone()[0]
            non_eligible_leaks = db.execute(
                """
                SELECT count(*)
                FROM v2_target_browser_ligand_context AS v
                JOIN canonical_ligand_targets AS c
                  ON c.ligand_instance_id = v.ligand_instance_id
                WHERE c.target_browser_eligible <> 'YES'
                """
            ).fetchone()[0]
            duplicate_occurrences = db.execute(
                """
                SELECT count(*)
                FROM (
                    SELECT ligand_instance_id
                    FROM v2_target_browser_ligand_context
                    GROUP BY ligand_instance_id
                    HAVING count(*) > 1
                )
                """
            ).fetchone()[0]
            group_only = db.execute(
                """
                SELECT virus_name, canonical_target_id, canonical_target_name
                FROM v2_target_browser_groups
                EXCEPT
                SELECT virus_name, canonical_target_id, canonical_target_name
                FROM v2_target_browser_ligand_context
                """
            ).fetchall()
            occurrence_only = db.execute(
                """
                SELECT virus_name, canonical_target_id, canonical_target_name
                FROM v2_target_browser_ligand_context
                GROUP BY virus_name, canonical_target_id, canonical_target_name
                EXCEPT
                SELECT virus_name, canonical_target_id, canonical_target_name
                FROM v2_target_browser_groups
                """
            ).fetchall()
        finally:
            db.close()

        self.assertEqual((occurrence_count, distinct_target_ids, composite_groups), (5414, 31, 34))
        self.assertEqual(group_count, 34)
        self.assertEqual(non_eligible_leaks, 0)
        self.assertEqual(duplicate_occurrences, 0)
        self.assertEqual(group_only, [])
        self.assertEqual(occurrence_only, [])

    def test_target_rows_keep_stable_canonical_id_and_source_provenance(self):
        db = vlismod_app.connect_db_row()
        try:
            assessment_rows = vlismod_app._load_canonical_target_browser_assessment_rows(db)
            readiness_rows, warhead_rows, attachment_rows = vlismod_app._load_protacability_enrichment_tables(db)
        finally:
            db.close()

        self.assertEqual(len(assessment_rows), 5414)
        fixture = next(row for row in assessment_rows if row["ligand_instance_id"] == 36170)
        self.assertEqual(fixture["canonical_target_id"], "protease")
        self.assertEqual(fixture["source_protein_type"], "capsid_protein,protease")

        decorated = vlismod_app._decorate_protacability_rows(
            assessment_rows,
            readiness_rows=readiness_rows,
            warhead_rows=warhead_rows,
            attachment_rows=attachment_rows,
        )
        groups = vlismod_app._group_target_rows(decorated)
        self.assertEqual(len(groups), 34)
        self.assertEqual(
            len({(row["virus_name"], row["canonical_target_id"]) for row in groups}),
            34,
        )
        hiv_protease = next(
            row for row in groups
            if row["virus_name"] == "HIV_1" and row["canonical_target_id"] == "protease"
        )
        self.assertEqual(hiv_protease["target_key"], "HIV_1::protease")
        self.assertIn("capsid_protein,protease", hiv_protease["source_protein_types"])

    def test_pass10_merges_hiv_protease_alias_without_losing_provenance(self):
        self.assertTrue(PASS10_AUTHORITY.exists())
        db = sqlite3.connect(DATABASE)
        try:
            duplicate_aliases = db.execute(
                """
                SELECT count(*) FROM canonical_ligand_targets
                WHERE virus_name='HIV_1'
                  AND target_browser_eligible='YES'
                  AND canonical_target_id='hiv_1_protease'
                """
            ).fetchone()[0]
            protease_rows = db.execute(
                """
                SELECT count(*) FROM canonical_ligand_targets
                WHERE virus_name='HIV_1'
                  AND target_browser_eligible='YES'
                  AND canonical_target_id='protease'
                """
            ).fetchone()[0]
            generic_nsp_rows = db.execute(
                """
                SELECT count(*) FROM canonical_ligand_targets
                WHERE target_browser_eligible='YES'
                  AND canonical_target_id='nsp_proteins'
                """
            ).fetchone()[0]
        finally:
            db.close()

        self.assertEqual(duplicate_aliases, 0)
        self.assertEqual(protease_rows, 748)
        self.assertEqual(generic_nsp_rows, 0)

    def test_target_browser_links_send_canonical_target_id(self):
        template = (ROOT / "templates" / "protacability_assessment.html").read_text()
        self.assertIn('data-target-id="${escapeHtml(row.canonical_target_id)}"', template)
        self.assertIn("params.set('canonical_target_id', button.data('target-id'))", template)
        self.assertIn("params.set('canonical_target_id', targetId)", template)


if __name__ == "__main__":
    unittest.main()
