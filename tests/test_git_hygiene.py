from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("check_git_hygiene", ROOT / "scripts/check_git_hygiene.py")
module = importlib.util.module_from_spec(spec)
assert spec.loader is not None
sys.modules[spec.name] = module
spec.loader.exec_module(module)


class GitHygieneTests(unittest.TestCase):
    def test_runtime_corpus_and_outputs_are_rejected(self):
        self.assertIsNotNone(module.forbidden_reason("PDB_FILES/HIV_1/protease/3EKY.cif"))
        self.assertIsNotNone(module.forbidden_reason("REVIEWER_REPRODUCIBILITY/PDB_FILES/.download_cache/3EKY.cif"))
        self.assertIsNotNone(module.forbidden_reason("CIF_DATABASE_REBUILD/outputs/arpeggio/x.log"))
        self.assertIsNotNone(module.forbidden_reason("foo.db-wal"))

    def test_reviewer_source_inputs_and_fixture_are_allowed(self):
        self.assertIsNone(module.forbidden_reason("REVIEWER_REPRODUCIBILITY/fixture/PDB_FILES/HIV_1/protease/3EKY.cif"))
        self.assertIsNone(module.forbidden_reason("REVIEWER_REPRODUCIBILITY/manifests/FROZEN_CIF_CORPUS_MANIFEST.csv"))
        self.assertIsNone(module.forbidden_reason("REVIEWER_REPRODUCIBILITY/inputs/chemistry/frozen_component_chemistry.csv"))
        self.assertIsNone(module.forbidden_reason("app.py"))


if __name__ == "__main__":
    unittest.main()
