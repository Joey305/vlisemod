import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
WORKSPACE_JS = (ROOT / "static" / "js" / "workspace-ui.js").read_text()
INDEXER = (ROOT / "templates" / "ligand_query.html").read_text()
COMPARISON = (ROOT / "templates" / "compare_ligands.html").read_text()
PROTAC = (ROOT / "templates" / "protacability_assessment.html").read_text()


class WorkspaceUiTests(unittest.TestCase):
    def test_occurrence_card_uses_canonical_identity_attributes(self):
        self.assertIn("data-ligand-instance-id", WORKSPACE_JS)
        self.assertIn("data-ligand-code", WORKSPACE_JS)
        self.assertIn("data-pdb-code", WORKSPACE_JS)

    def test_selection_updates_canonical_shareable_urls(self):
        self.assertIn("ligand_instance_id: occurrenceId", INDEXER)
        self.assertIn("instances: ids.join(',')", COMPARISON)
        self.assertIn("history.replaceState({}, '', canonicalUrl)", INDEXER)
        self.assertIn("history.replaceState({}, '', canonicalUrl)", COMPARISON)

    def test_synonym_display_does_not_define_scientific_url_identity(self):
        self.assertIn('value="${item.ligand_code}"', INDEXER)
        self.assertIn('value="${item.ligand_code}"', COMPARISON)

    def test_shared_no_data_and_copy_link_behaviors_exist(self):
        self.assertIn("function noData", WORKSPACE_JS)
        self.assertIn("Link copied", WORKSPACE_JS)
        self.assertIn("Interaction evidence unavailable", INDEXER)

    def test_current_method_generations_are_centralized(self):
        for version in (
            "legacy_mcs_etkdg_uff_cif_v2.5",
            "biopython-shrake_rupley-1.40-cif-v2.1",
            "arpeggio-cif-v2.2",
            "cif-ligand-geometry-v2.4",
            "rdkit-smarts-functional-groups-v2.3",
            "protacability-cif-v2.8",
            "attachment-sites-cif-v2.6",
        ):
            self.assertIn(version, WORKSPACE_JS)

    def test_protacability_uses_attachment_site_terminology_in_new_ui(self):
        self.assertIn("Attachment sites", PROTAC)
        self.assertIn("Attachment-site map", PROTAC)
        self.assertIn("Download PROTACability evidence", PROTAC)

    def test_attachment_site_selection_is_occurrence_resolved(self):
        self.assertIn("ligand_instance_id: nextInstanceId", PROTAC)
        self.assertIn("renderAttachmentSitesPanel(detail.attachment_sites || {}, detail.summary || {})", PROTAC)
        self.assertIn("ligand_instance_id", PROTAC)


if __name__ == "__main__":
    unittest.main()
