import importlib.util
import sys
import unittest
from pathlib import Path
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("vlismod_remote_proxy_app", ROOT / "app.py")
app_module = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = app_module
SPEC.loader.exec_module(app_module)


class RemoteProtacabilityProxyTests(unittest.TestCase):
    def setUp(self):
        app_module.app.config.update(TESTING=True)
        self.client = app_module.app.test_client()

    def _remote_payload(self, path, *, params=None, max_bytes=None):
        if path == "protacability/filter-options":
            return {"data_available": True, "virus_names": ["HIV_1"], "protein_types": ["protease"], "ligands": []}
        if path == "protacability/search":
            return {"data_available": True, "view": "targets", "rows": [], "summary": {}, "limit": 50, "offset": 0, "total_rows": 0, "has_more": False, "sort": "ligand_priority_desc"}
        raise AssertionError(path)

    def test_filter_options_never_fetches_unfiltered_source_in_randy_mode(self):
        with patch.object(app_module, "_normalized_backend_mode", return_value="randy"), \
             patch.object(app_module, "_remote_protacability_get", side_effect=self._remote_payload) as remote, \
             patch.object(app_module, "_load_protacability_source_payload", side_effect=AssertionError("raw source must not be fetched")):
            response = self.client.get("/api/protacability/filter_options?view=targets")
        self.assertEqual(response.status_code, 200)
        self.assertEqual(remote.call_args.args[0], "protacability/filter-options")

    def test_target_search_proxies_randy_not_local_source_or_database(self):
        with patch.object(app_module, "_normalized_backend_mode", return_value="randy"), \
             patch.object(app_module, "_remote_protacability_get", side_effect=self._remote_payload) as remote, \
             patch.object(app_module, "_load_protacability_source_payload", side_effect=AssertionError("raw source must not be fetched")), \
             patch.object(app_module, "connect_db_row", side_effect=AssertionError("target browser must proxy Randy")):
            response = self.client.get("/api/protacability/search?view=targets&page_size=25&canonical_target_id=protease")
        self.assertEqual(response.status_code, 200)
        self.assertEqual(remote.call_args.args[0], "protacability/search")

