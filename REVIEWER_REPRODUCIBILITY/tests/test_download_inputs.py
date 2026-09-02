from __future__ import annotations

import csv
import hashlib
import importlib.util
import io
import sys
import tempfile
import unittest
from pathlib import Path
from urllib.error import HTTPError, URLError


ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("download_inputs", ROOT / "tools/download_inputs.py")
downloader = importlib.util.module_from_spec(spec)
assert spec.loader is not None
sys.modules[spec.name] = downloader
spec.loader.exec_module(downloader)


def cif(entry_id: str) -> bytes:
    return f"data_{entry_id}\n_entry.id {entry_id}\n".encode()


class DownloadInputTests(unittest.TestCase):
    def make_manifest(self, directory: Path, content: bytes, duplicate: bool = True) -> Path:
        path = directory / "manifest.csv"
        checksum = hashlib.sha256(content).hexdigest()
        rows = [{"entry_id": "3EKY", "filename": "3EKY.cif", "relative_path": "HIV_1/capsid_protein/3EKY.cif", "sha256": checksum, "file_size": str(len(content))}]
        if duplicate:
            rows.append({"entry_id": "3EKY", "filename": "3EKY.cif", "relative_path": "HIV_1/protease/3EKY.cif", "sha256": checksum, "file_size": str(len(content))})
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
            writer.writeheader(); writer.writerows(rows)
        return path

    def test_manifest_parsing_and_duplicate_entry_grouping(self):
        with tempfile.TemporaryDirectory() as tmp:
            manifest = self.make_manifest(Path(tmp), cif("3EKY"))
            rows = downloader.read_manifest(manifest)
            grouped = downloader.group_by_entry(rows)
            self.assertEqual(len(rows), 2)
            self.assertEqual(list(grouped), ["3EKY"])
            self.assertEqual(len(grouped["3EKY"]), 2)

    def test_dry_run_never_calls_network_or_writes_source(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); manifest = self.make_manifest(root, cif("3EKY"))
            def opener(*_, **__): raise AssertionError("dry run must not download")
            result = downloader.download_inputs(manifest, root / "PDB_FILES", root / "reports", opener=opener, dry_run=True)
            self.assertEqual(result["summary"]["planned_downloads"], 1)
            self.assertFalse((root / "PDB_FILES").exists())
            self.assertFalse((root / "reports").exists())

    def test_download_once_and_materialize_exact_hierarchy(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); content = cif("3EKY"); manifest = self.make_manifest(root, content)
            calls = []
            def opener(request, timeout): calls.append(request.full_url); return io.BytesIO(content)
            result = downloader.download_inputs(manifest, root / "PDB_FILES", root / "reports", workers=1, opener=opener, sleeper=lambda _: None)
            self.assertEqual(calls, ["https://files.rcsb.org/download/3EKY.cif"])
            self.assertEqual(result["summary"]["materialized_hierarchy_files"], 2)
            for rel in ("HIV_1/capsid_protein/3EKY.cif", "HIV_1/protease/3EKY.cif"):
                self.assertEqual((root / "PDB_FILES" / rel).read_bytes(), content)
            verified = downloader.verify_inputs(manifest, root / "PDB_FILES", root / "reports")
            self.assertTrue(verified["summary"]["passed"])

    def test_valid_cache_skips_network_and_part_file_is_ignored(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); content = cif("3EKY"); manifest = self.make_manifest(root, content, duplicate=False)
            cache = root / "PDB_FILES/.download_cache"; cache.mkdir(parents=True)
            (cache / "3EKY.cif").write_bytes(content)
            (cache / "3EKY.interrupted.part").write_bytes(b"partial")
            def opener(*_): raise AssertionError("network must not be called for valid cache")
            result = downloader.download_inputs(manifest, root / "PDB_FILES", root / "reports", workers=1, opener=opener)
            self.assertEqual(result["summary"]["cache_hits"], 1)
            self.assertTrue((cache / "3EKY.interrupted.part").exists())

    def test_corrupt_cache_retries_then_replaces(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); content = cif("3EKY"); manifest = self.make_manifest(root, content, duplicate=False)
            cache = root / "PDB_FILES/.download_cache"; cache.mkdir(parents=True); (cache / "3EKY.cif").write_text("not cif")
            calls = []
            def opener(*_, **__): calls.append(1); return io.BytesIO(content)
            result = downloader.download_inputs(manifest, root / "PDB_FILES", root / "reports", workers=1, opener=opener, sleeper=lambda _: None)
            self.assertEqual(len(calls), 1)
            self.assertEqual(result["entries"][0]["status"], "VERIFIED_FROZEN_INPUT")

    def test_retry_not_found_parse_entry_and_checksum_classifications(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); content = cif("3EKY"); manifest = self.make_manifest(root, content, duplicate=False)
            attempts = []
            def retry_opener(*_, **__):
                attempts.append(1)
                if len(attempts) == 1: raise URLError("temporary")
                return io.BytesIO(content)
            result = downloader.download_inputs(manifest, root / "ok", root / "reports", workers=1, retries=2, opener=retry_opener, sleeper=lambda _: None)
            self.assertEqual(len(attempts), 2); self.assertEqual(result["entries"][0]["status"], "VERIFIED_FROZEN_INPUT")
            def missing_opener(request, **_): raise HTTPError(request.full_url, 404, "missing", {}, None)
            missing = downloader.download_inputs(manifest, root / "missing", root / "reports2", workers=1, opener=missing_opener)
            self.assertEqual(missing["entries"][0]["status"], "NOT_FOUND")
            def bad_parse(*_, **__): return io.BytesIO(b"not a CIF")
            parsed = downloader.download_inputs(manifest, root / "bad", root / "reports3", workers=1, opener=bad_parse)
            self.assertEqual(parsed["entries"][0]["status"], "PARSE_FAILED")
            def bad_entry(*_, **__): return io.BytesIO(cif("1ABC"))
            wrong = downloader.download_inputs(manifest, root / "wrong", root / "reports4", workers=1, opener=bad_entry)
            self.assertEqual(wrong["entries"][0]["status"], "ENTRY_ID_MISMATCH")
            def upstream(*_, **__): return io.BytesIO(cif("3EKY") + b"# upstream bytes\n")
            changed = downloader.download_inputs(manifest, root / "changed", root / "reports5", workers=1, opener=upstream, allow_current_upstream=True)
            self.assertEqual(changed["entries"][0]["status"], "UPSTREAM_REVISION_CHANGED")
            with self.assertRaises(RuntimeError): downloader.verify_inputs(manifest, root / "changed", root / "reports5")
            allowed = downloader.verify_inputs(manifest, root / "changed", root / "reports5", allow_current_upstream=True)
            self.assertTrue(allowed["summary"]["passed"])


if __name__ == "__main__":
    unittest.main()
