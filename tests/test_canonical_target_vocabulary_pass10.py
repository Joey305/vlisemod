import importlib.util
import unittest
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
PASS9 = ROOT / "outputs" / "protein_taxonomy_reconciliation_pass9" / "canonical_target_occurrences_pass9.csv"
PASS10 = ROOT / "outputs" / "protein_taxonomy_reconciliation_pass10" / "canonical_target_occurrences_pass10.csv"
SPEC = importlib.util.spec_from_file_location(
    "canonical_target_vocabulary_pass10", ROOT / "normalize_canonical_target_vocabulary_pass10.py"
)
normalizer = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(normalizer)


class CanonicalTargetVocabularyPass10Tests(unittest.TestCase):
    def setUp(self):
        self.pass9 = pd.read_csv(PASS9, dtype=str).fillna("")
        self.pass10 = pd.read_csv(PASS10, dtype=str).fillna("")

    def test_hiv_protease_aliases_merge_to_one_canonical_identity(self):
        rows = self.pass10[
            (self.pass10["virus_name"] == "HIV_1")
            & (self.pass10["final_target_browser_eligible"] == "YES")
            & self.pass10["final_canonical_target_id"].isin(["protease", "hiv_1_protease"])
        ]
        self.assertEqual(set(rows["final_canonical_target_id"]), {"protease"})
        self.assertEqual(len(rows), 748)

    def test_generic_nsp_bucket_is_resolved_by_entity_description(self):
        rows = self.pass10[
            (self.pass10["virus_name"] == "SARS_CoV_2")
            & (self.pass10["final_target_browser_eligible"] == "YES")
        ]
        self.assertNotIn("nsp_proteins", set(rows["final_canonical_target_id"]))
        self.assertIn("nsp3", set(rows["final_canonical_target_id"]))
        self.assertIn("nsp14", set(rows["final_canonical_target_id"]))
        self.assertIn("nsp16", set(rows["final_canonical_target_id"]))

    def test_source_aliases_do_not_split_a_same_virus_display_concept(self):
        yes = self.pass10[self.pass10["final_target_browser_eligible"] == "YES"]
        conflicts = (
            yes.groupby(["virus_name", "final_canonical_target_name"])["final_canonical_target_id"]
            .nunique()
        )
        self.assertTrue(conflicts[conflicts > 1].empty)

    def test_normalizer_keeps_dr7_protease_provenance(self):
        row = self.pass10[self.pass10["ligand_instance_id"].astype(str).eq("36170")].iloc[0]
        self.assertEqual(row["final_canonical_target_id"], "protease")
        self.assertEqual(row["current_stage14_protein_type"], "capsid_protein,protease")


if __name__ == "__main__":
    unittest.main()
