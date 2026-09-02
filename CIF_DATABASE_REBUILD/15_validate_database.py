#!/usr/bin/env python3
"""Final read-only release gate for the CIF-native V-LiSEMOD rebuild.

This replaces the original foundation-only Stage 15 validator. It validates the
frozen corpus, occurrence-level traceability, and the final materialized outputs
from Stages 07-14 without recalculating scientific data.

The validator is intentionally strict for the frozen release represented by the
expected counts below. A failed check exits non-zero and should block release.
It never writes to the SQLite database. Reports are written under outputs/.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sqlite3
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable

try:
    from importlib import import_module
    c = import_module("00_common")
    ROOT = Path(c.ROOT)
except Exception:
    c = None
    ROOT = Path(__file__).resolve().parent

VERSION = "final-release-validation-cif-v2.4"

# Frozen evidence generations.
MAPPING_VERSION = "legacy_mcs_etkdg_uff_cif_v2.5"
SASA_VERSION = "biopython-shrake_rupley-1.40-cif-v2.1"
ARPEGGIO_VERSION = "arpeggio-cif-v2.2"
GEOMETRY_VERSION = "cif-ligand-geometry-v2.4"
FUNCTIONAL_GROUP_VERSION = "rdkit-smarts-functional-groups-v2.3"
PROTACABILITY_VERSION = "protacability-cif-v2.8"
ATTACHMENT_VERSION = "attachment-sites-cif-v2.6"
COMPATIBILITY_VERSION = "compatibility-views-cif-v2.7"

EXPECTED = {
    "manifest_cifs": 11533,
    "structures": 7610,
    "included_ligands": 7355,
    "resolved_chemistry": 7335,
    "unresolved_chemistry": 20,
    "mapping_complete": 3864,
    "mapping_altloc_complete": 570,
    "mapping_partial": 2182,
    "mapping_remediation": 719,
    "mapping_downstream_eligible": 6616,
    "arpeggio_complete": 7355,
    "geometry_complete": 7308,
    "geometry_no_protein": 47,
    "functional_group_complete": 7335,
    "functional_group_matches": 53224,
    "functional_group_mapped_atoms": 134419,
    "functional_group_unmapped_atoms": 24472,
    "protacability_warhead_rows": 7355,
    "target_context_applicable": 6784,
    "target_context_no_contacting_chain": 524,
    "target_context_no_protein": 47,
    "target_assessed_instances": 6784,
    "target_chain_rows": 9462,
    "attachment_summary_instances": 7355,
    "attachment_atom_rows": 227080,
    "attachment_candidate_atoms": 62392,
    "attachment_high_atoms": 1544,
    "attachment_instances_with_candidates": 5814,
    "attachment_instances_with_high": 1120,
    "compatibility_views": 23,
}

EXPECTED_MAPPING_STATUSES = {
    "complete": 3864,
    "complete_altloc_resolved": 570,
    "partial_ccd_difference": 2182,
    "skipped_pending_remediation": 719,
}
EXPECTED_TARGET_CONTEXT = {
    "applicable_contacting_protein_chain": 6784,
    "not_applicable_no_contacting_protein_chain": 524,
    "not_applicable_no_protein_atoms": 47,
}
EXPECTED_READINESS_TIERS = {
    "High degrader-design readiness": 6579,
    "Moderate degrader-design readiness": 2289,
    "Weak degrader-design readiness": 421,
    "Exploratory degrader-design readiness": 173,
}
EXPECTED_ATTACHMENT_TIERS = {
    "Low attachment-site priority": 164688,
    "Exploratory attachment-site priority": 53606,
    "Moderate attachment-site priority": 7242,
    "High attachment-site priority": 1544,
}
EXPECTED_CHEMICAL_ROLES = {
    "unclassified_atom_context": 116334,
    "functional_group_context_only": 83144,
    "conditional_substitution_site": 18530,
    "direct_attachment_atom": 9072,
}
EXPECTED_VIEWS = [
    "v2_structure_context",
    "v2_ligand_context",
    "v2_ligand_atom_evidence",
    "v2_ligand_comparison_atom_contacts",
    "v2_all_chain_lysine_geometry",
    "v2_target_lysine_accessibility",
    "v2_protacability_target_context",
    "v2_protacability_best",
    "v2_attachment_site_candidates",
    "v2_attachment_site_high_priority",
    "v2_attachment_site_summary",
    "Virus_Proteins",
    "Ligand_Atoms_Smiles",
    "Functional_GROUPED",
    "ligand_atoms",
    "solvent_exposed_atoms",
    "RUPLEY_SASA_DATA",
    "SMILES_MAP_PDB",
    "Functional_Group_Atoms",
    "Arpeggio_Contacts_Data",
    "receptor_binding_pocket",
    "Covalent_Noncovalent",
    "distal_atoms",
]
EXPECTED_SCRIPT_VERSIONS = {
    "07_map_cif_atoms_to_smiles.py": MAPPING_VERSION,
    "08_calculate_ligand_sasa.py": SASA_VERSION,
    "09_run_arpeggio.py": ARPEGGIO_VERSION,
    "10_calculate_ligand_geometry.py": GEOMETRY_VERSION,
    "11_assign_functional_groups.py": FUNCTIONAL_GROUP_VERSION,
    "12_build_protacability.py": PROTACABILITY_VERSION,
    "13_build_attachment_sites.py": ATTACHMENT_VERSION,
    "14_build_compatibility_views.py": COMPATIBILITY_VERSION,
}

REQUIRED_TABLES = [
    "structures", "ligands", "ligand_instances", "ligand_instance_atoms",
    "analysis_runs", "pipeline_failures", "mapping_remediation_queue",
    "ligand_mapping_runs", "ligand_smiles_atom_mapping", "ligand_sasa_atoms",
    "ligand_arpeggio_runs", "arpeggio_raw_contact_labels", "arpeggio_unique_atom_pairs",
    "protacability_ligand_inventory", "ligand_atom_geometry", "target_chain_geometry",
    "target_surface_lysines", "ligand_functional_group_matches",
    "ligand_functional_group_atoms", "ligand_functional_group_summary",
    "protacability_warhead_linkability", "protacability_target_context",
    "protacability_assessment", "protacability_degrader_readiness",
    "protacability_attachment_sites", "protacability_attachment_site_summary",
]


@dataclass
class Check:
    section: str
    name: str
    status: str
    observed: Any = None
    expected: Any = None
    detail: str = ""


class Gate:
    def __init__(self) -> None:
        self.checks: list[Check] = []

    def add(self, section: str, name: str, ok: bool, observed: Any = None,
            expected: Any = None, detail: str = "") -> None:
        self.checks.append(Check(section, name, "PASS" if ok else "FAIL", observed, expected, detail))

    def warn(self, section: str, name: str, observed: Any = None,
             expected: Any = None, detail: str = "") -> None:
        self.checks.append(Check(section, name, "WARN", observed, expected, detail))

    @property
    def failures(self) -> list[Check]:
        return [x for x in self.checks if x.status == "FAIL"]

    @property
    def warnings(self) -> list[Check]:
        return [x for x in self.checks if x.status == "WARN"]

    @property
    def passed(self) -> bool:
        return not self.failures


def progress(message: str) -> None:
    stamp = time.strftime("%H:%M:%S")
    print(f"[Stage15 {stamp}] {message}", flush=True)


def connect_readonly(database: str) -> sqlite3.Connection:
    # URI read-only mode prevents accidental mutation by the release gate.
    uri = Path(database).resolve().as_uri() + "?mode=ro"
    db = sqlite3.connect(uri, uri=True, timeout=60.0)
    db.row_factory = sqlite3.Row
    db.execute("PRAGMA foreign_keys=ON")
    db.execute("PRAGMA busy_timeout=60000")
    db.execute("PRAGMA query_only=ON")
    db.execute("PRAGMA automatic_index=ON")
    db.execute("PRAGMA temp_store=MEMORY")
    return db


def scalar(db: sqlite3.Connection, sql: str, params: Iterable[Any] = ()) -> Any:
    row = db.execute(sql, tuple(params)).fetchone()
    return row[0] if row is not None else None


def grouped_counts(db: sqlite3.Connection, sql: str, params: Iterable[Any] = ()) -> dict[str, int]:
    return {str(r[0]): int(r[1]) for r in db.execute(sql, tuple(params)).fetchall()}


def table_exists(db: sqlite3.Connection, name: str) -> bool:
    return bool(scalar(db, "SELECT count(*) FROM sqlite_master WHERE type='table' AND name=?", (name,)))


def view_exists(db: sqlite3.Connection, name: str) -> bool:
    return bool(scalar(db, "SELECT count(*) FROM sqlite_master WHERE type='view' AND name=?", (name,)))


def safe_check(gate: Gate, section: str, name: str, fn, expected: Any = None,
               predicate=None, detail: str = "") -> Any:
    try:
        observed = fn()
        ok = predicate(observed) if predicate else observed == expected
        gate.add(section, name, ok, observed, expected, detail)
        return observed
    except Exception as exc:
        gate.add(section, name, False, f"{type(exc).__name__}: {exc}", expected, detail)
        return None


def manifest_count(path: Path) -> int:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return max(0, sum(1 for _ in csv.reader(handle)) - 1)


def script_version(path: Path) -> str | None:
    if not path.exists():
        return None
    text = path.read_text(encoding="utf-8", errors="replace")
    match = re.search(r"\bVERSION\s*=\s*['\"]([^'\"]+)['\"]", text)
    return match.group(1) if match else None


def view_sql(db: sqlite3.Connection, name: str) -> str:
    row = db.execute("SELECT sql FROM sqlite_master WHERE type='view' AND name=?", (name,)).fetchone()
    return (row[0] or "") if row else ""


def validate(database: str, root: Path = ROOT) -> dict[str, Any]:
    gate = Gate()
    metrics: dict[str, Any] = {}
    db = connect_readonly(database)
    try:
        progress("Schema, compatibility-view metadata, and code-version locks")
        # ------------------------------------------------------------------
        # Schema and script-generation locks
        # ------------------------------------------------------------------
        for name in REQUIRED_TABLES:
            gate.add("Schema", f"table exists: {name}", table_exists(db, name), int(table_exists(db, name)), 1)
        for name in EXPECTED_VIEWS:
            gate.add("Compatibility", f"view exists: {name}", view_exists(db, name), int(view_exists(db, name)), 1)

        existing_views = sum(1 for name in EXPECTED_VIEWS if view_exists(db, name))
        metrics["compatibility_views_present"] = existing_views
        gate.add("Compatibility", "all expected compatibility views present",
                 existing_views == EXPECTED["compatibility_views"], existing_views, EXPECTED["compatibility_views"])

        for filename, expected_version in EXPECTED_SCRIPT_VERSIONS.items():
            observed_version = script_version(root / filename)
            gate.add("Code freeze", f"{filename} VERSION", observed_version == expected_version,
                     observed_version, expected_version)

        progress("Foundation corpus and denominator checks")
        # ------------------------------------------------------------------
        # Foundation
        # ------------------------------------------------------------------
        manifest = root / "manifests" / "FROZEN_CIF_CORPUS_MANIFEST.csv"
        if manifest.exists():
            n_manifest = manifest_count(manifest)
            metrics["manifest_cif_files"] = n_manifest
            gate.add("Foundation", "frozen mmCIF manifest rows", n_manifest == EXPECTED["manifest_cifs"],
                     n_manifest, EXPECTED["manifest_cifs"])
        else:
            gate.add("Foundation", "frozen mmCIF manifest exists", False, str(manifest), "existing file")

        n_structures = scalar(db, "SELECT count(*) FROM structures")
        n_pdb = scalar(db, "SELECT count(DISTINCT entry_id) FROM structures")
        included = scalar(db, "SELECT count(*) FROM ligand_instances WHERE curation_status='included'")
        resolved = scalar(db, """SELECT count(*) FROM ligand_instances i JOIN ligands l ON l.ligand_id=i.ligand_id
                                WHERE i.curation_status='included' AND l.chemical_status='resolved'""")
        unresolved = included - resolved
        metrics.update({
            "structures": n_structures,
            "unique_pdb_entries": n_pdb,
            "included_ligand_instances": included,
            "resolved_chemistry_instances": resolved,
            "unresolved_chemistry_instances": unresolved,
        })
        gate.add("Foundation", "structure rows", n_structures == EXPECTED["structures"], n_structures, EXPECTED["structures"])
        gate.add("Foundation", "unique PDB entries", n_pdb == EXPECTED["structures"], n_pdb, EXPECTED["structures"])
        gate.add("Foundation", "retained ligand instances", included == EXPECTED["included_ligands"], included, EXPECTED["included_ligands"])
        gate.add("Foundation", "resolved-chemistry denominator", resolved == EXPECTED["resolved_chemistry"], resolved, EXPECTED["resolved_chemistry"])
        gate.add("Foundation", "included instances without resolved chemistry", unresolved == EXPECTED["unresolved_chemistry"], unresolved, EXPECTED["unresolved_chemistry"])

        progress("Stage 07 mapping checks")
        # ------------------------------------------------------------------
        # Stage 07 mapping
        # ------------------------------------------------------------------
        mapping_total = scalar(db, "SELECT count(*) FROM ligand_mapping_runs WHERE method_version=?", (MAPPING_VERSION,))
        mapping_status = grouped_counts(db,
            "SELECT mapping_status,count(*) FROM ligand_mapping_runs WHERE method_version=? GROUP BY mapping_status",
            (MAPPING_VERSION,))
        mapping_eligible = scalar(db, "SELECT count(*) FROM ligand_mapping_runs WHERE method_version=? AND downstream_mapping_eligibility=1", (MAPPING_VERSION,))
        pending_remediation = scalar(db, "SELECT count(*) FROM mapping_remediation_queue WHERE remediation_status='pending'")
        metrics.update({
            "mapping_total": mapping_total,
            "mapping_statuses": mapping_status,
            "mapping_downstream_eligible": mapping_eligible,
            "mapping_pending_remediation": pending_remediation,
        })
        gate.add("Stage 07 mapping", "mapping population", mapping_total == EXPECTED["resolved_chemistry"], mapping_total, EXPECTED["resolved_chemistry"])
        gate.add("Stage 07 mapping", "mapping status distribution", mapping_status == EXPECTED_MAPPING_STATUSES, mapping_status, EXPECTED_MAPPING_STATUSES)
        gate.add("Stage 07 mapping", "downstream-eligible mappings", mapping_eligible == EXPECTED["mapping_downstream_eligible"], mapping_eligible, EXPECTED["mapping_downstream_eligible"])
        gate.add("Stage 07 mapping", "pending remediation population", pending_remediation == EXPECTED["mapping_remediation"], pending_remediation, EXPECTED["mapping_remediation"])
        map_bad = scalar(db, "SELECT count(*) FROM ligand_mapping_runs WHERE method_version=? AND mapping_status IN ('failed','mapping_timeout')", (MAPPING_VERSION,))
        gate.add("Stage 07 mapping", "failed/time-out mappings", map_bad == 0, map_bad, 0)
        map_dup_smiles = scalar(db, """SELECT count(*) FROM (
            SELECT ligand_instance_id,smiles_atom_index,count(*) n
            FROM ligand_smiles_atom_mapping
            WHERE method_version=? AND ligand_instance_atom_id IS NOT NULL AND smiles_atom_index IS NOT NULL
            GROUP BY ligand_instance_id,smiles_atom_index HAVING n>1)""", (MAPPING_VERSION,))
        gate.add("Stage 07 mapping", "unique SMILES index per occurrence", map_dup_smiles == 0, map_dup_smiles, 0)

        progress("Stage 08 SASA checks")
        # ------------------------------------------------------------------
        # Stage 08 SASA
        # ------------------------------------------------------------------
        sasa_instances = scalar(db, "SELECT count(DISTINCT ligand_instance_id) FROM ligand_sasa_atoms WHERE method_version=? AND status='complete'", (SASA_VERSION,))
        sasa_rows = scalar(db, "SELECT count(*) FROM ligand_sasa_atoms WHERE method_version=? AND status='complete'", (SASA_VERSION,))
        metrics.update({"sasa_instances": sasa_instances, "sasa_atom_rows": sasa_rows})
        gate.add("Stage 08 SASA", "instance coverage", sasa_instances == EXPECTED["included_ligands"], sasa_instances, EXPECTED["included_ligands"])
        sasa_param_bad = scalar(db, """SELECT count(*) FROM ligand_sasa_atoms WHERE method_version=? AND
            (status<>'complete' OR abs(probe_radius-1.40)>0.000001 OR point_density<>100 OR water_treatment<>'HOH_removed')""", (SASA_VERSION,))
        gate.add("Stage 08 SASA", "1.40 A / 100-point / HOH-removed parameter consistency", sasa_param_bad == 0, sasa_param_bad, 0)
        progress("Stage 08: selected-atom SASA coverage (set-difference scan)")
        sasa_missing = scalar(db, """SELECT count(*) FROM (
            SELECT a.ligand_instance_atom_id
            FROM ligand_instance_atoms a
            JOIN ligand_instances i ON i.ligand_instance_id=a.ligand_instance_id
            WHERE i.curation_status='included' AND a.selected_conformer=1
            EXCEPT
            SELECT s.ligand_instance_atom_id
            FROM ligand_sasa_atoms s
            WHERE s.method_version=? AND s.status='complete'
        )""", (SASA_VERSION,))
        gate.add("Stage 08 SASA", "selected ligand atoms missing SASA", sasa_missing == 0, sasa_missing, 0)

        progress("Stage 09 Arpeggio checks")
        # ------------------------------------------------------------------
        # Stage 09 Arpeggio
        # ------------------------------------------------------------------
        progress("Stage 09: latest Arpeggio outcome distribution")
        arpeggio_latest_status = grouped_counts(db, """WITH latest_id AS (
            SELECT r.ligand_instance_id,MAX(r.run_id) AS run_id
            FROM ligand_arpeggio_runs r
            JOIN ligand_instances i ON i.ligand_instance_id=r.ligand_instance_id
            WHERE i.curation_status='included'
            GROUP BY r.ligand_instance_id
        )
        SELECT r.status,count(DISTINCT r.ligand_instance_id)
        FROM ligand_arpeggio_runs r
        JOIN latest_id x ON x.ligand_instance_id=r.ligand_instance_id AND x.run_id=r.run_id
        GROUP BY r.status""")
        arpeggio_latest_completed = arpeggio_latest_status.get('completed', 0)
        arpeggio_latest_noncomplete = sum(v for k,v in arpeggio_latest_status.items() if k != 'completed')
        arpeggio_latest_duplicate_keys = scalar(db, """SELECT count(*) FROM (
            WITH latest_id AS (
                SELECT r.ligand_instance_id,MAX(r.run_id) AS run_id
                FROM ligand_arpeggio_runs r
                JOIN ligand_instances i ON i.ligand_instance_id=r.ligand_instance_id
                WHERE i.curation_status='included'
                GROUP BY r.ligand_instance_id
            )
            SELECT r.ligand_instance_id,r.run_id
            FROM ligand_arpeggio_runs r
            JOIN latest_id x ON x.ligand_instance_id=r.ligand_instance_id AND x.run_id=r.run_id
            GROUP BY r.ligand_instance_id,r.run_id
            HAVING count(*)>1
        )""")
        arpeggio_excluded_latest_completed = scalar(db, """SELECT count(DISTINCT r.ligand_instance_id)
            FROM ligand_arpeggio_runs r
            JOIN ligand_instances i ON i.ligand_instance_id=r.ligand_instance_id
            JOIN (
                SELECT r2.ligand_instance_id,MAX(r2.run_id) AS run_id
                FROM ligand_arpeggio_runs r2
                JOIN ligand_instances i2 ON i2.ligand_instance_id=r2.ligand_instance_id
                WHERE i2.curation_status<>'included'
                GROUP BY r2.ligand_instance_id
            ) x ON x.ligand_instance_id=r.ligand_instance_id AND x.run_id=r.run_id
            WHERE i.curation_status<>'included' AND r.status='completed'""")
        progress("Stage 09: included-occurrence coverage (set-difference scan)")
        arpeggio_missing = scalar(db, """SELECT count(*) FROM (
            SELECT ligand_instance_id FROM ligand_instances WHERE curation_status='included'
            EXCEPT
            SELECT ligand_instance_id FROM ligand_arpeggio_runs
        )""")
        metrics.update({
            "arpeggio_latest_completed": arpeggio_latest_completed,
            "arpeggio_latest_noncomplete": arpeggio_latest_noncomplete,
            "arpeggio_latest_duplicate_keys": arpeggio_latest_duplicate_keys,
            "arpeggio_excluded_latest_completed": arpeggio_excluded_latest_completed,
            "arpeggio_missing": arpeggio_missing,
        })
        gate.add("Stage 09 Arpeggio", "latest completed outcomes in included release population", arpeggio_latest_completed == EXPECTED["arpeggio_complete"], arpeggio_latest_completed, EXPECTED["arpeggio_complete"])
        gate.add("Stage 09 Arpeggio", "latest unresolved outcomes in included release population", arpeggio_latest_noncomplete == 0, arpeggio_latest_noncomplete, 0)
        gate.add("Stage 09 Arpeggio", "duplicate latest outcome keys in included release population", arpeggio_latest_duplicate_keys == 0, arpeggio_latest_duplicate_keys, 0)
        gate.add("Stage 09 Arpeggio", "included instances without Arpeggio outcome", arpeggio_missing == 0, arpeggio_missing, 0)

        progress("Stage 10 geometry checks")
        # ------------------------------------------------------------------
        # Stage 10 geometry / surface lysines
        # ------------------------------------------------------------------
        geometry_status = grouped_counts(db,
            "SELECT status,count(*) FROM protacability_ligand_inventory WHERE method_version=? GROUP BY status",
            (GEOMETRY_VERSION,))
        geometry_total = sum(geometry_status.values())
        metrics.update({"geometry_statuses": geometry_status, "geometry_total": geometry_total})
        gate.add("Stage 10 geometry", "inventory coverage", geometry_total == EXPECTED["included_ligands"], geometry_total, EXPECTED["included_ligands"])
        gate.add("Stage 10 geometry", "protein-applicable instances", geometry_status.get("complete", 0) == EXPECTED["geometry_complete"], geometry_status.get("complete", 0), EXPECTED["geometry_complete"])
        gate.add("Stage 10 geometry", "no-protein not-applicable instances", geometry_status.get("not_applicable_no_protein_atoms", 0) == EXPECTED["geometry_no_protein"], geometry_status.get("not_applicable_no_protein_atoms", 0), EXPECTED["geometry_no_protein"])
        geometry_other = geometry_total - geometry_status.get("complete", 0) - geometry_status.get("not_applicable_no_protein_atoms", 0)
        gate.add("Stage 10 geometry", "unexpected geometry statuses", geometry_other == 0, geometry_other, 0)

        progress("Stage 11 functional-group checks")
        # ------------------------------------------------------------------
        # Stage 11 functional groups
        # ------------------------------------------------------------------
        fg_summary = db.execute("""SELECT count(*) n,COALESCE(sum(functional_group_match_count),0) matches,
            COALESCE(sum(mapped_functional_group_atom_count),0) mapped_atoms,
            COALESCE(sum(unmapped_functional_group_atom_count),0) unmapped_atoms
            FROM ligand_functional_group_summary WHERE method_version=? AND status='complete'""", (FUNCTIONAL_GROUP_VERSION,)).fetchone()
        metrics.update({
            "functional_group_complete": int(fg_summary["n"]),
            "functional_group_matches": int(fg_summary["matches"]),
            "functional_group_mapped_atoms": int(fg_summary["mapped_atoms"]),
            "functional_group_unmapped_atoms": int(fg_summary["unmapped_atoms"]),
        })
        gate.add("Stage 11 functional groups", "complete summaries", fg_summary["n"] == EXPECTED["functional_group_complete"], int(fg_summary["n"]), EXPECTED["functional_group_complete"])
        gate.add("Stage 11 functional groups", "SMARTS match count", fg_summary["matches"] == EXPECTED["functional_group_matches"], int(fg_summary["matches"]), EXPECTED["functional_group_matches"])
        gate.add("Stage 11 functional groups", "mapped functional-group atom occurrences", fg_summary["mapped_atoms"] == EXPECTED["functional_group_mapped_atoms"], int(fg_summary["mapped_atoms"]), EXPECTED["functional_group_mapped_atoms"])
        gate.add("Stage 11 functional groups", "unmapped functional-group atom occurrences", fg_summary["unmapped_atoms"] == EXPECTED["functional_group_unmapped_atoms"], int(fg_summary["unmapped_atoms"]), EXPECTED["functional_group_unmapped_atoms"])
        fg_full_unmapped = scalar(db, """SELECT COALESCE(sum(s.unmapped_functional_group_atom_count),0)
            FROM ligand_functional_group_summary s JOIN ligand_mapping_runs r ON r.ligand_instance_id=s.ligand_instance_id
            WHERE s.method_version=? AND r.method_version=? AND r.mapping_status IN ('complete','complete_altloc_resolved')""",
            (FUNCTIONAL_GROUP_VERSION, MAPPING_VERSION))
        gate.add("Stage 11 functional groups", "unmapped FG atoms on complete/altloc-complete mappings", fg_full_unmapped == 0, fg_full_unmapped, 0)
        fg_bad_mapping_label = scalar(db, """SELECT count(*) FROM ligand_functional_group_atoms f
            WHERE f.method_version=? AND f.ligand_instance_atom_id IS NOT NULL AND f.mapping_status<>'mapped_element_validated'""", (FUNCTIONAL_GROUP_VERSION,))
        gate.add("Stage 11 functional groups", "mapped FG atoms lacking element validation", fg_bad_mapping_label == 0, fg_bad_mapping_label, 0)

        progress("Stage 12 PROTACability checks")
        # ------------------------------------------------------------------
        # Stage 12 PROTACability target context and surface-lysine scoring
        # ------------------------------------------------------------------
        warhead_rows = scalar(db, "SELECT count(*) FROM protacability_warhead_linkability WHERE method_version=?", (PROTACABILITY_VERSION,))
        context_counts = grouped_counts(db,
            "SELECT target_context_status,count(*) FROM protacability_target_context WHERE method_version=? GROUP BY target_context_status",
            (PROTACABILITY_VERSION,))
        assessed_instances = scalar(db, "SELECT count(DISTINCT ligand_instance_id) FROM protacability_assessment WHERE method_version=?", (PROTACABILITY_VERSION,))
        assessment_rows = scalar(db, "SELECT count(*) FROM protacability_assessment WHERE method_version=?", (PROTACABILITY_VERSION,))
        readiness_instances = scalar(db, "SELECT count(DISTINCT ligand_instance_id) FROM protacability_degrader_readiness WHERE method_version=?", (PROTACABILITY_VERSION,))
        readiness_rows = scalar(db, "SELECT count(*) FROM protacability_degrader_readiness WHERE method_version=?", (PROTACABILITY_VERSION,))
        readiness_tiers = grouped_counts(db,
            "SELECT degrader_design_readiness_tier,count(*) FROM protacability_degrader_readiness WHERE method_version=? GROUP BY degrader_design_readiness_tier",
            (PROTACABILITY_VERSION,))
        metrics.update({
            "warhead_rows": warhead_rows,
            "target_context": context_counts,
            "target_assessed_instances": assessed_instances,
            "target_chain_rows": assessment_rows,
            "readiness_instances": readiness_instances,
            "readiness_rows": readiness_rows,
            "readiness_tiers": readiness_tiers,
        })
        gate.add("Stage 12 PROTACability", "warhead/linkability coverage", warhead_rows == EXPECTED["protacability_warhead_rows"], warhead_rows, EXPECTED["protacability_warhead_rows"])
        gate.add("Stage 12 PROTACability", "target-context distribution", context_counts == EXPECTED_TARGET_CONTEXT, context_counts, EXPECTED_TARGET_CONTEXT)
        gate.add("Stage 12 PROTACability", "assessed ligand instances", assessed_instances == EXPECTED["target_assessed_instances"], assessed_instances, EXPECTED["target_assessed_instances"])
        gate.add("Stage 12 PROTACability", "target-chain assessment rows", assessment_rows == EXPECTED["target_chain_rows"], assessment_rows, EXPECTED["target_chain_rows"])
        gate.add("Stage 12 PROTACability", "readiness instance coverage", readiness_instances == EXPECTED["target_assessed_instances"], readiness_instances, EXPECTED["target_assessed_instances"])
        gate.add("Stage 12 PROTACability", "readiness row count", readiness_rows == EXPECTED["target_chain_rows"], readiness_rows, EXPECTED["target_chain_rows"])
        gate.add("Stage 12 PROTACability", "readiness tier distribution", readiness_tiers == EXPECTED_READINESS_TIERS, readiness_tiers, EXPECTED_READINESS_TIERS)

        target_basis_bad = scalar(db, """SELECT count(*) FROM protacability_assessment WHERE method_version=?
            AND (target_chain_selection_basis<>'stage09_arpeggio_direct_protein_contact' OR ligand_target_contact_pair_count<=0)""", (PROTACABILITY_VERSION,))
        gate.add("Stage 12 PROTACability", "all target rows supported by direct Stage-09 protein contact", target_basis_bad == 0, target_basis_bad, 0)

        # Recompute the target-side score only from ALL-chain lysine NZ coverage and
        # surface exposure. If this passes, ligand-to-lysine distance cannot be
        # contributing to the stored target score.
        target_formula_bad = scalar(db, """SELECT count(*) FROM protacability_assessment WHERE method_version=? AND
            abs(protacability_proxy_score - round(CASE WHEN lys_count<=0 THEN 0.0 ELSE
                10.0*nz_observed_lys_count/lys_count +
                60.0*nz_exposed_lys_count_gt_1/lys_count +
                30.0*nz_exposed_lys_count_gt_5/lys_count END,2)) > 0.011""", (PROTACABILITY_VERSION,))
        gate.add("Stage 12 PROTACability", "surface-lysine score recomputes from NZ coverage/SASA only", target_formula_bad == 0, target_formula_bad, 0)

        readiness_formula_bad = scalar(db, """SELECT count(*) FROM protacability_degrader_readiness WHERE method_version=? AND
            abs(degrader_design_readiness_score -
                CASE WHEN nz_exposed_lys_count_gt_1<=0 THEN
                    min(round(0.60*warhead_linkability_score + 0.40*target_lysine_accessibility_score,2),34.99)
                ELSE round(0.60*warhead_linkability_score + 0.40*target_lysine_accessibility_score,2) END) > 0.011""", (PROTACABILITY_VERSION,))
        gate.add("Stage 12 PROTACability", "readiness score recomputes from warhead + surface lysines only", readiness_formula_bad == 0, readiness_formula_bad, 0)

        proximity_scoring_bad = scalar(db, """SELECT count(*) FROM protacability_degrader_readiness WHERE method_version=? AND
            (min_lys_ligand_distance_a IS NOT NULL OR near_ligand_exposed_lys_count<>0 OR ternary_geometry_cue_score<>0)""", (PROTACABILITY_VERSION,))
        gate.add("Stage 12 PROTACability", "ligand-to-lysine proximity disabled in readiness rows", proximity_scoring_bad == 0, proximity_scoring_bad, 0)

        zero_nz_high = scalar(db, """SELECT count(*) FROM protacability_degrader_readiness WHERE method_version=?
            AND nz_exposed_lys_count_gt_1=0 AND degrader_design_readiness_tier IN
            ('High degrader-design readiness','Moderate degrader-design readiness')""", (PROTACABILITY_VERSION,))
        gate.add("Stage 12 PROTACability", "zero accessible NZ cannot be Moderate/High", zero_nz_high == 0, zero_nz_high, 0)

        best_chain_bad = scalar(db, """SELECT count(*) FROM (
            SELECT ligand_instance_id,sum(best_chain_for_instance) n FROM protacability_degrader_readiness
            WHERE method_version=? GROUP BY ligand_instance_id HAVING n<>1)""", (PROTACABILITY_VERSION,))
        gate.add("Stage 12 PROTACability", "exactly one best target chain per assessed ligand", best_chain_bad == 0, best_chain_bad, 0)

        progress("Stage 13 attachment-site checks")
        # ------------------------------------------------------------------
        # Stage 13 attachment sites
        # ------------------------------------------------------------------
        att_summary_rows = scalar(db, "SELECT count(*) FROM protacability_attachment_site_summary WHERE method_version=? AND status='complete'", (ATTACHMENT_VERSION,))
        att_atom_rows = scalar(db, "SELECT count(*) FROM protacability_attachment_sites WHERE method_version=?", (ATTACHMENT_VERSION,))
        att_candidate = scalar(db, "SELECT COALESCE(sum(candidate_attachment_atom),0) FROM protacability_attachment_sites WHERE method_version=?", (ATTACHMENT_VERSION,))
        att_high = scalar(db, "SELECT COALESCE(sum(high_priority_attachment_atom),0) FROM protacability_attachment_sites WHERE method_version=?", (ATTACHMENT_VERSION,))
        att_candidate_instances = scalar(db, "SELECT count(DISTINCT CASE WHEN candidate_attachment_atom=1 THEN ligand_instance_id END) FROM protacability_attachment_sites WHERE method_version=?", (ATTACHMENT_VERSION,))
        att_high_instances = scalar(db, "SELECT count(DISTINCT CASE WHEN high_priority_attachment_atom=1 THEN ligand_instance_id END) FROM protacability_attachment_sites WHERE method_version=?", (ATTACHMENT_VERSION,))
        att_tiers = grouped_counts(db, "SELECT attachment_priority_tier,count(*) FROM protacability_attachment_sites WHERE method_version=? GROUP BY attachment_priority_tier", (ATTACHMENT_VERSION,))
        att_roles = grouped_counts(db, "SELECT atom_chemical_role,count(*) FROM protacability_attachment_sites WHERE method_version=? GROUP BY atom_chemical_role", (ATTACHMENT_VERSION,))
        metrics.update({
            "attachment_summary_instances": att_summary_rows,
            "attachment_atom_rows": att_atom_rows,
            "attachment_candidate_atoms": att_candidate,
            "attachment_high_atoms": att_high,
            "attachment_instances_with_candidates": att_candidate_instances,
            "attachment_instances_with_high": att_high_instances,
            "attachment_tiers": att_tiers,
            "attachment_chemical_roles": att_roles,
        })
        gate.add("Stage 13 attachment sites", "complete instance summaries", att_summary_rows == EXPECTED["attachment_summary_instances"], att_summary_rows, EXPECTED["attachment_summary_instances"])
        gate.add("Stage 13 attachment sites", "atom rows", att_atom_rows == EXPECTED["attachment_atom_rows"], att_atom_rows, EXPECTED["attachment_atom_rows"])
        gate.add("Stage 13 attachment sites", "candidate atoms", att_candidate == EXPECTED["attachment_candidate_atoms"], att_candidate, EXPECTED["attachment_candidate_atoms"])
        gate.add("Stage 13 attachment sites", "High-priority atoms", att_high == EXPECTED["attachment_high_atoms"], att_high, EXPECTED["attachment_high_atoms"])
        gate.add("Stage 13 attachment sites", "instances with candidates", att_candidate_instances == EXPECTED["attachment_instances_with_candidates"], att_candidate_instances, EXPECTED["attachment_instances_with_candidates"])
        gate.add("Stage 13 attachment sites", "instances with High sites", att_high_instances == EXPECTED["attachment_instances_with_high"], att_high_instances, EXPECTED["attachment_instances_with_high"])
        gate.add("Stage 13 attachment sites", "priority tier distribution", att_tiers == EXPECTED_ATTACHMENT_TIERS, att_tiers, EXPECTED_ATTACHMENT_TIERS)
        gate.add("Stage 13 attachment sites", "atom chemical-role distribution", att_roles == EXPECTED_CHEMICAL_ROLES, att_roles, EXPECTED_CHEMICAL_ROLES)

        candidate_rule_bad = scalar(db, """SELECT count(*) FROM protacability_attachment_sites WHERE method_version=? AND
            candidate_attachment_atom <> CASE WHEN upper(COALESCE(element,'')) NOT IN ('H','D','T')
            AND mapped=1 AND solvent_exposed=1 AND COALESCE(strong_contact_count,0)=0 THEN 1 ELSE 0 END""", (ATTACHMENT_VERSION,))
        gate.add("Stage 13 attachment sites", "candidate-core rule equivalence", candidate_rule_bad == 0, candidate_rule_bad, 0)

        high_rule_bad = scalar(db, """SELECT count(*) FROM protacability_attachment_sites WHERE method_version=?
            AND high_priority_attachment_atom=1 AND (mapped<>1 OR solvent_exposed<>1 OR COALESCE(strong_contact_count,0)<>0
            OR direct_attachment_support<>1 OR atom_chemical_role<>'direct_attachment_atom'
            OR points_away_from_pocket<>1 OR local_corridor_clear<>1)""", (ATTACHMENT_VERSION,))
        gate.add("Stage 13 attachment sites", "High-rule violations", high_rule_bad == 0, high_rule_bad, 0)

        high_tier_flag_bad = scalar(db, """SELECT count(*) FROM protacability_attachment_sites WHERE method_version=? AND
            ((attachment_priority_tier='High attachment-site priority' AND high_priority_attachment_atom<>1)
             OR (high_priority_attachment_atom=1 AND attachment_priority_tier<>'High attachment-site priority'))""", (ATTACHMENT_VERSION,))
        gate.add("Stage 13 attachment sites", "High tier/flag consistency", high_tier_flag_bad == 0, high_tier_flag_bad, 0)

        progress("Occurrence traceability checks")
        # ------------------------------------------------------------------
        # Occurrence scoping / cross-instance joins
        # ------------------------------------------------------------------
        cross_queries = {
            "mapping atom": ("ligand_smiles_atom_mapping", MAPPING_VERSION, "m.method_version=?", "m"),
            "SASA atom": ("ligand_sasa_atoms", SASA_VERSION, "m.method_version=?", "m"),
            "geometry atom": ("ligand_atom_geometry", GEOMETRY_VERSION, "m.method_version=?", "m"),
            "functional-group atom": ("ligand_functional_group_atoms", FUNCTIONAL_GROUP_VERSION, "m.method_version=?", "m"),
            "attachment-site atom": ("protacability_attachment_sites", ATTACHMENT_VERSION, "m.method_version=?", "m"),
        }
        for label, (table, version, condition, alias) in cross_queries.items():
            progress(f"Occurrence traceability: {label}")
            bad = scalar(db, f"""SELECT count(*) FROM {table} {alias}
                JOIN ligand_instance_atoms a ON a.ligand_instance_atom_id={alias}.ligand_instance_atom_id
                WHERE {condition} AND {alias}.ligand_instance_id<>a.ligand_instance_id""", (version,))
            gate.add("Occurrence traceability", f"cross-instance {label} joins", bad == 0, bad, 0)
        progress("Occurrence traceability: raw Arpeggio ligand-atom joins")
        arpeggio_cross = scalar(db, """SELECT count(*) FROM arpeggio_raw_contact_labels r
            JOIN ligand_instance_atoms a ON a.ligand_instance_atom_id=r.ligand_instance_atom_id
            WHERE r.ligand_instance_atom_id IS NOT NULL AND r.ligand_instance_id<>a.ligand_instance_id""")
        gate.add("Occurrence traceability", "cross-instance Arpeggio ligand-atom joins", arpeggio_cross == 0, arpeggio_cross, 0)

        progress("Stage 14 compatibility-view metadata checks")
        # ------------------------------------------------------------------
        # Stage 14 compatibility views: metadata-only checks (never COUNT(*)
        # every complex view; that caused the multi-hour report bottleneck).
        # ------------------------------------------------------------------
        latest_compat = db.execute("""SELECT run_id,status,processed_count,success_count,failure_count,parameters_json
            FROM analysis_runs WHERE stage='compatibility_views' ORDER BY run_id DESC LIMIT 1""").fetchone()
        if latest_compat:
            compat_tuple = (latest_compat["status"], latest_compat["processed_count"], latest_compat["success_count"], latest_compat["failure_count"])
            gate.add("Stage 14 compatibility", "latest build completed 23/23/0",
                     compat_tuple == ("completed", 23, 23, 0), compat_tuple, ("completed", 23, 23, 0))
        else:
            gate.add("Stage 14 compatibility", "latest compatibility build exists", False, None, "completed run")

        view_requirements = {
            "v2_target_lysine_accessibility": [PROTACABILITY_VERSION, GEOMETRY_VERSION, "protacability_assessment"],
            "v2_protacability_target_context": [PROTACABILITY_VERSION],
            "v2_protacability_best": [PROTACABILITY_VERSION, "best_chain_for_instance=1"],
            "v2_attachment_site_candidates": [ATTACHMENT_VERSION, "candidate_attachment_atom=1"],
            "v2_attachment_site_high_priority": [ATTACHMENT_VERSION, "high_priority_attachment_atom=1", "direct_attachment_support=1", "direct_attachment_atom"],
            "v2_attachment_site_summary": [ATTACHMENT_VERSION],
        }
        for name, tokens in view_requirements.items():
            sql = re.sub(r"\s+", "", view_sql(db, name)).lower()
            missing_tokens = [token for token in tokens if re.sub(r"\s+", "", token).lower() not in sql]
            gate.add("Stage 14 compatibility", f"{name} pinned to current evidence", not missing_tokens,
                     "all tokens present" if not missing_tokens else f"missing: {missing_tokens}", "all tokens present")

        progress("3EKY/DR7 regression fixture")
        # ------------------------------------------------------------------
        # 3EKY / DR7 regression fixture
        # ------------------------------------------------------------------
        dr7_rows = db.execute("""SELECT i.ligand_instance_id FROM ligand_instances i
            JOIN structures s ON s.structure_id=i.structure_id
            WHERE s.entry_id='3EKY' AND i.label_comp_id='DR7' AND i.curation_status='included'""").fetchall()
        gate.add("3EKY/DR7 regression", "unique included DR7 occurrence", len(dr7_rows) == 1, len(dr7_rows), 1)
        if len(dr7_rows) == 1:
            dr7 = int(dr7_rows[0][0])
            metrics["3EKY_DR7_ligand_instance_id"] = dr7
            mr = db.execute("""SELECT mapping_status,mapped_count,structural_atom_count,smiles_atom_count,
                heavy_atom_mapping_fraction,downstream_mapping_eligibility FROM ligand_mapping_runs
                WHERE ligand_instance_id=? AND method_version=?""", (dr7, MAPPING_VERSION)).fetchone()
            expected_map = ("complete", 51, 51, 51, 1.0, 1)
            observed_map = tuple(mr) if mr else None
            map_ok = bool(mr and mr["mapping_status"] == "complete" and mr["mapped_count"] == 51 and
                          mr["structural_atom_count"] == 51 and mr["smiles_atom_count"] == 51 and
                          abs(float(mr["heavy_atom_mapping_fraction"] or 0)-1.0) < 1e-9 and
                          mr["downstream_mapping_eligibility"] == 1)
            gate.add("3EKY/DR7 regression", "51/51 authoritative mapping", map_ok, observed_map, expected_map)

            dr7_fg_bad = scalar(db, """SELECT count(*) FROM ligand_functional_group_atoms
                WHERE ligand_instance_id=? AND method_version=? AND ligand_instance_atom_id IS NOT NULL
                  AND mapping_status<>'mapped_element_validated'""", (dr7, FUNCTIONAL_GROUP_VERSION))
            gate.add("3EKY/DR7 regression", "functional-group mapped atoms element-validated", dr7_fg_bad == 0, dr7_fg_bad, 0)

            tc = db.execute("""SELECT target_context_status,contacting_protein_chain_count,contacting_protein_chain_ids,
                target_chain_selection_basis FROM protacability_target_context
                WHERE ligand_instance_id=? AND method_version=?""", (dr7, PROTACABILITY_VERSION)).fetchone()
            target_ok = bool(tc and tc["target_context_status"] == "applicable_contacting_protein_chain" and
                             tc["contacting_protein_chain_count"] == 2 and tc["contacting_protein_chain_ids"] == "A;B" and
                             tc["target_chain_selection_basis"] == "stage09_arpeggio_direct_protein_contact")
            gate.add("3EKY/DR7 regression", "direct target chains A;B", target_ok, dict(tc) if tc else None,
                     {"status": "applicable_contacting_protein_chain", "count": 2, "chains": "A;B"})

            scores = {r["chain_id"]: float(r["protacability_proxy_score"]) for r in db.execute(
                "SELECT chain_id,protacability_proxy_score FROM protacability_assessment WHERE ligand_instance_id=? AND method_version=?",
                (dr7, PROTACABILITY_VERSION)).fetchall()}
            score_ok = set(scores) == {"A", "B"} and abs(scores["A"]-52.86) <= 0.01 and abs(scores["B"]-81.43) <= 0.01
            gate.add("3EKY/DR7 regression", "surface-lysine scores A=52.86, B=81.43", score_ok, scores, {"A": 52.86, "B": 81.43})

            dr7_high = scalar(db, "SELECT count(*) FROM protacability_attachment_sites WHERE ligand_instance_id=? AND method_version=? AND high_priority_attachment_atom=1", (dr7, ATTACHMENT_VERSION))
            gate.add("3EKY/DR7 regression", "no forced High direct-attachment site", dr7_high == 0, dr7_high, 0)
            dr7_moderate = [r[0] for r in db.execute("""SELECT exact_atom FROM protacability_attachment_sites
                WHERE ligand_instance_id=? AND method_version=? AND attachment_priority_tier='Moderate attachment-site priority'
                ORDER BY exact_atom""", (dr7, ATTACHMENT_VERSION)).fetchall()]
            gate.add("3EKY/DR7 regression", "conditional Moderate sites", dr7_moderate == ["CAO", "CAR", "CAS", "NBD"],
                     dr7_moderate, ["CAO", "CAR", "CAS", "NBD"])

        progress("SQLite integrity and foreign-key checks")
        # ------------------------------------------------------------------
        # SQLite health
        # ------------------------------------------------------------------
        integrity = scalar(db, "PRAGMA integrity_check")
        fk_rows = db.execute("PRAGMA foreign_key_check").fetchall()
        metrics["integrity_check"] = integrity
        metrics["foreign_key_errors"] = len(fk_rows)
        gate.add("Database health", "PRAGMA integrity_check", integrity == "ok", integrity, "ok")
        gate.add("Database health", "PRAGMA foreign_key_check", len(fk_rows) == 0, len(fk_rows), 0)

    finally:
        db.close()

    progress("Validation queries complete; writing reports")
    return {
        "validator_version": VERSION,
        "database": str(Path(database).resolve()),
        "passed": gate.passed,
        "failure_count": len(gate.failures),
        "warning_count": len(gate.warnings),
        "metrics": metrics,
        "checks": [asdict(x) for x in gate.checks],
    }


def markdown_report(result: dict[str, Any]) -> str:
    checks = [Check(**x) for x in result["checks"]]
    sections: list[str] = []
    overall = "PASS" if result["passed"] else "FAIL"
    sections.extend([
        "# V-LiSEMOD final release validation",
        "",
        f"* Validator: `{result['validator_version']}`",
        f"* Database: `{result['database']}`",
        f"* Final release validation: **{overall}**",
        f"* Failed checks: {result['failure_count']}",
        f"* Warnings: {result['warning_count']}",
        "",
        "This is a read-only release gate. It validates the frozen evidence generations and does not recalculate structural analyses.",
    ])
    order: list[str] = []
    for item in checks:
        if item.section not in order:
            order.append(item.section)
    for section in order:
        sections.extend(["", f"## {section}", ""])
        for item in (x for x in checks if x.section == section):
            symbol = "PASS" if item.status == "PASS" else item.status
            exp = f"; expected={item.expected}" if item.expected is not None else ""
            det = f"; {item.detail}" if item.detail else ""
            sections.append(f"* **{symbol}** - {item.name}: observed={item.observed}{exp}{det}")

    m = result["metrics"]
    sections.extend([
        "",
        "## Frozen release snapshot for manuscript / reviewer response",
        "",
        f"* Frozen mmCIF files: {m.get('manifest_cif_files', 'NA')}",
        f"* Unique PDB entries: {m.get('unique_pdb_entries', 'NA')}",
        f"* Retained ligand instances: {m.get('included_ligand_instances', 'NA')}",
        f"* Resolved-chemistry instances: {m.get('resolved_chemistry_instances', 'NA')}",
        f"* Stage-07 downstream-usable mappings: {m.get('mapping_downstream_eligible', 'NA')}",
        f"* Stage-07 pending remediation instances: {m.get('mapping_pending_remediation', 'NA')}",
        f"* Stage-09 completed Arpeggio outcomes in included release population: {m.get('arpeggio_latest_completed', 'NA')}",
        f"* Excluded/non-release ligand instances with completed Arpeggio history: {m.get('arpeggio_excluded_latest_completed', 'NA')} (reported for provenance; outside the 7,355 release denominator)",
        f"* Stage-10 protein-applicable geometry: {m.get('geometry_statuses', {}).get('complete', 'NA')}",
        f"* Stage-10 no-protein target N/A: {m.get('geometry_statuses', {}).get('not_applicable_no_protein_atoms', 'NA')}",
        f"* Stage-12 directly contacting target contexts: {m.get('target_context', {}).get('applicable_contacting_protein_chain', 'NA')}",
        f"* Stage-12 target-chain assessments: {m.get('target_chain_rows', 'NA')}",
        f"* Stage-13 ligand instances with candidate attachment sites: {m.get('attachment_instances_with_candidates', 'NA')}",
        f"* Stage-13 ligand instances with High direct attachment sites: {m.get('attachment_instances_with_high', 'NA')}",
        f"* Stage-13 candidate atoms: {m.get('attachment_candidate_atoms', 'NA')}",
        f"* Stage-13 High-priority direct attachment atoms: {m.get('attachment_high_atoms', 'NA')}",
        "",
        "### Interpretation guardrails",
        "",
        "* The 7,355 retained-ligand denominator covers the included structure-derived release population; Arpeggio history for excluded/non-release occurrences is retained for provenance but is not counted in that denominator. The 7,335 denominator is the resolved-chemistry population eligible for SMILES-based mapping and SMARTS annotation.",
        "* Target PROTACability uses Stage-09-confirmed ligand-contacting protein chains, then evaluates lysine NZ solvent accessibility across the entire selected chain surface.",
        "* Ligand-to-lysine distance is retained only in upstream descriptive geometry where present and is not used in the v2.8 target-accessibility or degrader-readiness score.",
        "* High attachment-site priority requires a mapped, solvent-exposed, low-strong-contact atom with direct atom-level attachment chemistry, an outward-facing vector, and a locally clear corridor.",
        "* These are structural-priority heuristics for follow-up design, not experimentally calibrated predictors of linker tolerance, ubiquitination, or degradation.",
        "",
    ])
    return "\n".join(sections)


def write_reports(result: dict[str, Any], root: Path, report_path: str | None = None,
                  json_path: str | None = None) -> tuple[Path, Path]:
    outdir = root / "outputs"
    outdir.mkdir(parents=True, exist_ok=True)
    md = Path(report_path) if report_path else outdir / "FINAL_RELEASE_VALIDATION_REPORT.md"
    js = Path(json_path) if json_path else outdir / "FINAL_RELEASE_VALIDATION_REPORT.json"
    if not md.is_absolute():
        md = root / md
    if not js.is_absolute():
        js = root / js
    md.parent.mkdir(parents=True, exist_ok=True)
    js.parent.mkdir(parents=True, exist_ok=True)
    md.write_text(markdown_report(result) + "\n", encoding="utf-8")
    js.write_text(json.dumps(result, indent=2, sort_keys=True, default=str) + "\n", encoding="utf-8")
    return md, js


def print_summary(result: dict[str, Any], md: Path, js: Path) -> None:
    status = "PASS" if result["passed"] else "FAIL"
    print(f"FINAL RELEASE VALIDATION: {status}")
    print(f"validator: {result['validator_version']}")
    print(f"failed_checks: {result['failure_count']}")
    print(f"warnings: {result['warning_count']}")
    print(f"report: {md}")
    print(f"json_report: {js}")
    if not result["passed"]:
        print("\nFailed checks:")
        for item in result["checks"]:
            if item["status"] == "FAIL":
                print(f"- [{item['section']}] {item['name']}: observed={item['observed']} expected={item['expected']}")


def main() -> int:
    parser = argparse.ArgumentParser(description="Final read-only V-LiSEMOD release validation gate")
    parser.add_argument("--database", default=str(ROOT / "viral_data_cif_v2.db"))
    parser.add_argument("--root", default=str(ROOT), help="Pipeline root containing scripts/manifests/outputs")
    parser.add_argument("--report")
    parser.add_argument("--json-report")
    parser.add_argument("--no-fail-exit", action="store_true", help="Return exit code 0 even when release checks fail")
    args = parser.parse_args()

    root = Path(args.root).resolve()
    result = validate(args.database, root)
    md, js = write_reports(result, root, args.report, args.json_report)
    print_summary(result, md, js)
    if result["passed"] or args.no_fail_exit:
        return 0
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
