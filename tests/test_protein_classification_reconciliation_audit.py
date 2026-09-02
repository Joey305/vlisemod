import importlib.util
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "protein_classification_reconciliation_audit",
    ROOT / "scripts" / "audit_protein_classification_reconciliation.py",
)
audit = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(audit)


class ProteinClassificationReconciliationAuditTests(unittest.TestCase):
    def test_entity_normalizer_prefers_specific_protein_terms(self):
        self.assertEqual(audit.normalize_description("Exoribonuclease H, p66 RT"), "rnase_h")
        self.assertEqual(audit.normalize_description("RNA-directed RNA polymerase"), "polymerase")
        self.assertEqual(audit.normalize_description("HIV-1 protease"), "protease")

    def test_local_2o4k_entity_metadata_maps_both_contact_chains_to_protease(self):
        metadata = audit.cif_entity_metadata(str(ROOT / "PDB_FILES" / "HIV_1" / "protease" / "2O4K.cif"))
        self.assertFalse(metadata["error"])
        self.assertEqual(metadata["chains"]["A"]["entity_id"], "1")
        self.assertEqual(metadata["chains"]["A"]["normalized_label"], "protease")
        self.assertEqual(metadata["chains"]["B"]["normalized_label"], "protease")


if __name__ == "__main__":
    unittest.main()
