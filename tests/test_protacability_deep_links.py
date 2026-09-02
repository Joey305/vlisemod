import importlib.util
import sqlite3
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
APP_SPEC = importlib.util.spec_from_file_location("vlismod_protacability_deep_links_app", ROOT / "app.py")
vlismod_app = importlib.util.module_from_spec(APP_SPEC)
sys.modules[APP_SPEC.name] = vlismod_app
APP_SPEC.loader.exec_module(vlismod_app)


class QueryArgs(dict):
    def get(self, key, default=None, type=None):
        value = super().get(key, default)
        if type is not None and value not in (None, ""):
            return type(value)
        return value

    def getlist(self, key):
        value = self.get(key)
        return [] if value is None else (value if isinstance(value, list) else [value])


class ProtacabilityDeepLinkTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.db_path = Path(self.tmp.name) / "protacability-deep-links.db"
        with sqlite3.connect(self.db_path) as conn:
            conn.executescript(
                """
                CREATE TABLE protacability_ligand_inventory (
                    inventory_id INTEGER PRIMARY KEY,
                    ligand_instance_id INTEGER, pdb_code TEXT, model_id TEXT,
                    ligand_resname TEXT, ligand_chain TEXT, ligand_residue_id TEXT,
                    ligand_insertion_code TEXT
                );
                CREATE TABLE protacability_assessment (
                    assessment_id INTEGER PRIMARY KEY, ligand_instance_id INTEGER,
                    pdb_code TEXT, chain_id TEXT, virus_name TEXT, protein_type TEXT,
                    method_version TEXT
                );
                """
            )
            conn.executemany(
                "INSERT INTO protacability_ligand_inventory VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
                [
                    (1, 59418, "3EKY", "1", "DR7", "A", "100", ""),
                    (2, 59757, "3EL1", "1", "DR7", "A", "100", ""),
                    (3, 70000, "8CYB", "1", "AH2", "A", "1312", ""),
                ],
            )
            conn.executemany(
                "INSERT INTO protacability_assessment VALUES (?, ?, ?, ?, ?, ?, ?)",
                [
                    (1, 59418, "3EKY", "A", "HIV_1", "protease", vlismod_app.PROTACABILITY_METHOD_VERSION),
                    (2, 59418, "3EKY", "B", "HIV_1", "protease", vlismod_app.PROTACABILITY_METHOD_VERSION),
                    (3, 59757, "3EL1", "A", "HIV_1", "protease", vlismod_app.PROTACABILITY_METHOD_VERSION),
                ],
            )
        self.conn = sqlite3.connect(self.db_path)
        self.conn.row_factory = sqlite3.Row

    def tearDown(self):
        self.conn.close()
        self.tmp.cleanup()

    def test_protacability_exact_occurrence_deep_link(self):
        context, error = vlismod_app._resolve_protacability_occurrence_context(
            self.conn,
            QueryArgs({"ligand_instance_id": "59418", "ligand": "DR7", "pdb_code": "3EKY"}),
        )
        self.assertIsNone(error)
        self.assertEqual(context["ligand_resname"], "DR7")
        self.assertEqual(context["pdb_code"], "3EKY")

    def test_exact_occurrence_metadata_mismatch_is_rejected(self):
        _, error = vlismod_app._resolve_protacability_occurrence_context(
            self.conn, QueryArgs({"ligand_instance_id": "59418", "ligand": "BAD"})
        )
        self.assertIn("do not match", error)

    def test_exact_occurrence_rejects_pdb_mismatch(self):
        _, error = vlismod_app._resolve_protacability_occurrence_context(
            self.conn, QueryArgs({"ligand_instance_id": "59418", "pdb_code": "3EL1"})
        )
        self.assertIn("do not match", error)

    def test_protacability_invalid_occurrence_deep_link(self):
        _, error = vlismod_app._resolve_protacability_occurrence_context(
            self.conn, QueryArgs({"ligand_instance_id": "999999"})
        )
        self.assertEqual(error, "The requested ligand occurrence could not be found.")

    def test_invalid_occurrence_returns_clear_error(self):
        _, error = vlismod_app._resolve_protacability_occurrence_context(
            self.conn, QueryArgs({"ligand_instance_id": "not-an-id"})
        )
        self.assertEqual(error, "The requested ligand occurrence could not be found.")

    def test_protacability_ligand_deep_link(self):
        rows = vlismod_app._load_protacability_assessment_rows(self.conn, ligand="DR7")
        self.assertEqual({row["ligand_instance_id"] for row in rows}, {59418, 59757})

    def test_exact_assessment_load_is_sql_constrained(self):
        rows = vlismod_app._load_protacability_assessment_rows(
            self.conn, ligand="DR7", pdb_code="3EL1", ligand_instance_id=59757
        )
        self.assertEqual([(row["pdb_code"], row["ligand_instance_id"]) for row in rows], [("3EL1", 59757)])

    def test_ligand_indexer_protacability_exact_link(self):
        indexer = (ROOT / "templates" / "ligand_query.html").read_text()
        self.assertIn("ligand_instance_id: selectedOccurrence.val()", indexer)
        self.assertIn("target=\"_blank\"", indexer)

    def test_ligand_indexer_protacability_all_ligand_link(self):
        indexer = (ROOT / "templates" / "ligand_query.html").read_text()
        self.assertIn("new URLSearchParams({ ligand: selectedLigand })", indexer)

    def test_synonym_deep_link_uses_canonical_component_id(self):
        indexer = (ROOT / "templates" / "ligand_query.html").read_text()
        self.assertIn('value="${item.ligand_code}"', indexer)
        self.assertIn("ligand: selectedLigand", indexer)

    def test_ligand_comparison_protacability_exact_links(self):
        comparison = (ROOT / "templates" / "compare_ligands.html").read_text()
        self.assertIn("ligand_instance_id: occurrenceId", comparison)
        self.assertIn("noopener noreferrer", comparison)

    def test_ligand_indexer_chart_modal_keeps_exact_occurrence_link(self):
        template = (ROOT / "templates" / "ligand_query.html").read_text()
        self.assertIn('id="carouselProtacabilityLink"', template)
        self.assertIn("ligand_instance_id: selectedOccurrence.val()", template)
        self.assertIn("displayCarouselModal(payload.chart_urls, selectedLigand, selectedOccurrence)", template)

    def test_ligand_image_workspace_can_open_the_exact_occurrence(self):
        template = (ROOT / "templates" / "display_images.html").read_text()
        self.assertIn('id="assess-protacability-link"', template)
        self.assertIn("data-ligand-instance-id", template)
        self.assertIn("updateProtacabilityLink(checkbox)", template)

    def test_ligand_image_workspace_can_highlight_exact_attachment_sites_in_3d(self):
        template = (ROOT / "templates" / "display_images.html").read_text()
        self.assertIn('id="highlight-sasa-attachment-button"', template)
        self.assertIn('id="sasa-attachment-sites-table"', template)
        self.assertIn("/api/protacability/structure_detail/", template)
        self.assertIn("mapAttachmentAtomsToLigandIndices('viewport', atoms)", template)
        self.assertIn("highlightAttachmentRegionSets('viewport', regions", template)
        self.assertIn("displayLigandSpinRequested = true", template)
        self.assertIn("const hasLoaded3DViewer = viewport", template)
        self.assertIn("focusSasaAttachmentSite", template)
        self.assertIn('class="attachment-site-row"', template)
        self.assertLess(
            template.index("Display Ligand Interaction Diagram"),
            template.index('id="sasa-attachment-panel"'),
        )

    def test_ligand_comparison_protacability_all_ligand_link(self):
        comparison = (ROOT / "templates" / "compare_ligands.html").read_text()
        self.assertIn("Explore ${ligand} across PROTACability", comparison)

    def test_protacability_deep_link_filter_state_visible(self):
        template = (ROOT / "templates" / "protacability_assessment.html").read_text()
        self.assertIn("$('#protac-ligand-instance-id').val(params.get('ligand_instance_id') || '');", template)
        self.assertIn("return refreshDynamicFilterOptions({ skipLoader: true }).then(function()", template)
        self.assertIn("Viewing PROTACability for", template)

    def test_exact_deep_link_opens_the_selected_structure_summary(self):
        template = (ROOT / "templates" / "protacability_assessment.html").read_text()
        self.assertIn("openDeepLinkedOccurrenceDetail(data);", template)
        self.assertIn("ligand_instance_id: occurrenceId", template)
        self.assertIn("Selected ligand:", template)

    def test_attachment_region_badge_is_kept_on_one_line(self):
        template = (ROOT / "templates" / "protacability_assessment.html").read_text()
        stylesheet = (ROOT / "static" / "css" / "styles.css").read_text()
        self.assertIn("attachment-region-badge", template)
        self.assertIn("attachment-region-cell", template)
        self.assertIn("white-space: nowrap", stylesheet)

    def test_detail_popup_has_selected_ligand_handoff_actions(self):
        template = (ROOT / "templates" / "protacability_assessment.html").read_text()
        stylesheet = (ROOT / "static" / "css" / "styles.css").read_text()
        self.assertIn('id="protac-detail-footer"', template)
        self.assertIn("Use as Warhead", template)
        self.assertIn("Open PROTAC Builder", template)
        self.assertIn("function updateProtacDetailFooter", template)
        self.assertIn("#protac-detail-content", stylesheet)
        self.assertIn("overflow-y: auto", stylesheet)

    def test_no_protacability_record_stays_empty_for_the_exact_occurrence(self):
        context, error = vlismod_app._resolve_protacability_occurrence_context(
            self.conn, QueryArgs({"ligand_instance_id": "70000", "ligand": "AH2", "pdb_code": "8CYB"})
        )
        self.assertIsNone(error)
        rows = vlismod_app._load_protacability_assessment_rows(
            self.conn,
            ligand=context["ligand_resname"],
            pdb_code=context["pdb_code"],
            ligand_instance_id=context["ligand_instance_id"],
        )
        self.assertEqual(rows, [])


if __name__ == "__main__":
    unittest.main()
