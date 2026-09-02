from __future__ import annotations

import importlib.util
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("reproduce", ROOT / "reproduce.py")
reproduce = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(reproduce)


class ReproductionTests(unittest.TestCase):
    def test_canonical_stage_order_and_exclusions(self):
        keys = [key for _, key, _ in reproduce.STAGES]
        scripts = [script for _, _, script in reproduce.STAGES]
        self.assertEqual(keys, ["inventory", "create", "ingest", "curate", "chemistry", "remediation-registry", "mapping", "sasa", "arpeggio", "geometry", "functional-groups", "protacability", "attachment-sites", "compatibility-views", "validate"])
        self.assertNotIn("04_identify_ligand_instances.py", scripts)
        self.assertFalse(any("OLD" in script or "attachment_site.py" in script for script in scripts))
        self.assertEqual(reproduce.selected_stages("mapping", "validate")[0][1], "mapping")
        self.assertEqual(reproduce.selected_stages("mapping", "validate")[-1][1], "validate")

    def test_fixture_dry_run_creates_no_database(self):
        with tempfile.TemporaryDirectory() as tmp:
            database = Path(tmp) / "must_not_exist.db"
            completed = subprocess.run([sys.executable, str(ROOT / "reproduce.py"), "--fixture", "--dry-run", "--database", str(database)], text=True, capture_output=True)
            self.assertEqual(completed.returncode, 0, completed.stderr + completed.stdout)
            self.assertFalse(database.exists())
            self.assertIn("Stage 13: attachment-sites", completed.stdout)
            self.assertNotIn("Stage 15", completed.stdout)

    def test_fixture_manifest_is_exact(self):
        report = reproduce.verify_inputs(reproduce.FIXTURE_SOURCE, fixture=True)
        self.assertEqual(report["summary"]["manifest_rows"], 2)
        self.assertEqual(report["summary"]["missing"], 0)
        self.assertEqual(report["summary"]["checksum_mismatch"], 0)

    def test_release_fingerprint_reference_is_valid_json(self):
        import json
        reference = json.loads((ROOT / "reference/SCIENTIFIC_CONTENT_FINGERPRINTS.json").read_text())
        self.assertEqual(reference["version"], "scientific-content-fingerprint-v1")
        self.assertIn("attachment_sites", reference["components"])


if __name__ == "__main__":
    unittest.main()
