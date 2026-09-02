#!/usr/bin/env python3
"""Reviewer-facing orchestrator for the frozen CIF-native V-LiSEMOD release."""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import shutil
import sqlite3
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent
PIPELINE = ROOT / "pipeline"
OUTPUTS = ROOT / "outputs"
MANIFEST = ROOT / "manifests" / "FROZEN_CIF_CORPUS_MANIFEST.csv"
CHEMISTRY = ROOT / "inputs" / "chemistry" / "frozen_component_chemistry.csv"
REGISTRY = ROOT / "inputs" / "mapping" / "frozen_mapping_remediation_registry.csv"
FIXTURE_SOURCE = ROOT / "fixture" / "PDB_FILES"
REFERENCE_FINGERPRINTS = ROOT / "reference" / "SCIENTIFIC_CONTENT_FINGERPRINTS.json"

STAGES = [
    ("01", "inventory", "01_inventory_cifs.py"),
    ("02", "create", "02_create_database.py"),
    ("03", "ingest", "03_ingest_structures.py"),
    ("05", "curate", "05_curate_ligands.py"),
    ("06", "chemistry", "06_load_ligand_chemistry.py"),
    ("06a", "remediation-registry", "06a_import_mapping_remediation_registry.py"),
    ("07", "mapping", "07_map_cif_atoms_to_smiles.py"),
    ("08", "sasa", "08_calculate_ligand_sasa.py"),
    ("09", "arpeggio", "09_run_arpeggio.py"),
    ("10", "geometry", "10_calculate_ligand_geometry.py"),
    ("11", "functional-groups", "11_assign_functional_groups.py"),
    ("12", "protacability", "12_build_protacability.py"),
    ("13", "attachment-sites", "13_build_attachment_sites.py"),
    ("14", "compatibility-views", "14_build_compatibility_views.py"),
    ("15", "validate", "15_validate_database.py"),
]
STAGE_INDEX = {key: index for index, (_, key, _) in enumerate(STAGES)}
REQUIRED_MODULES = ("gemmi", "Bio", "numpy", "scipy", "rdkit", "pandas")
sys.path.insert(0, str(ROOT / "tools"))
from download_inputs import download_inputs as acquire_public_inputs, verify_inputs as verify_local_inputs  # noqa: E402


def now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def check_environment(require_arpeggio: bool = True) -> list[str]:
    missing = [name for name in REQUIRED_MODULES if importlib.util.find_spec(name) is None]
    executable = shutil.which("pdbe-arpeggio") or shutil.which("arpeggio")
    if require_arpeggio and not executable:
        missing.append("pdbe-arpeggio executable")
    print("Environment check")
    print(f"  Python: {sys.version.split()[0]}")
    print(f"  SQLite: {sqlite3.sqlite_version}")
    print(f"  pdbe-arpeggio: {executable or 'MISSING'}")
    print(f"  required modules: {'OK' if not [x for x in missing if x != 'pdbe-arpeggio executable'] else 'MISSING'}")
    if missing:
        print("  missing: " + ", ".join(missing))
    return missing


def verify_inputs(source: Path, fixture: bool = False, allow_current_upstream: bool = False,
                  entry_ids: set[str] | None = None) -> dict:
    """Local-only manifest validation; cache CIFs are intentionally ignored."""
    selected = {"3EKY"} if fixture else entry_ids
    return verify_local_inputs(MANIFEST, source, OUTPUTS, allow_current_upstream, selected)


def selected_stages(from_stage: str | None, to_stage: str | None):
    start = STAGE_INDEX[from_stage] if from_stage else 0
    end = STAGE_INDEX[to_stage] if to_stage else len(STAGES) - 1
    if start > end:
        raise ValueError("--from-stage must precede --to-stage")
    return STAGES[start : end + 1]


def run_command(command: list[str], dry_run: bool, event: dict) -> None:
    print("$ " + " ".join(command), flush=True)
    event["command"] = command
    if dry_run:
        event["status"] = "dry-run"
        return
    started = time.monotonic()
    completed = subprocess.run(command, cwd=ROOT, text=True)
    event["elapsed_seconds"] = round(time.monotonic() - started, 3)
    event["returncode"] = completed.returncode
    event["status"] = "completed" if completed.returncode == 0 else "failed"
    if completed.returncode:
        raise RuntimeError(f"Required stage failed ({completed.returncode}): {' '.join(command)}")


def current_mapping_timeouts(database: Path) -> list[int]:
    if not database.exists():
        return []
    with sqlite3.connect(database) as db:
        return [row[0] for row in db.execute("SELECT ligand_instance_id FROM ligand_mapping_runs WHERE method_version='legacy_mcs_etkdg_uff_cif_v2.5' AND mapping_status='mapping_timeout'")]


def current_arpeggio_unresolved(database: Path) -> int:
    if not database.exists():
        return 0
    with sqlite3.connect(database) as db:
        row = db.execute("""
            SELECT count(*) FROM ligand_instances i
            WHERE i.curation_status='included' AND NOT EXISTS (
              SELECT 1 FROM ligand_arpeggio_runs a
              WHERE a.ligand_instance_id=i.ligand_instance_id
                AND a.status='completed'
            )
        """).fetchone()
    return int(row[0])


def stage_command(key: str, script: str, database: Path, source: Path, args, fixture: bool) -> list[str]:
    command = [sys.executable, str(PIPELINE / script)]
    if key == "inventory":
        return command + ["--source", str(source), "--output-dir", str(OUTPUTS / "run_inventory")] + (["--pdb-id", "3EKY"] if fixture else [])
    if key == "create":
        return command + ["--database", str(database)]
    if key == "ingest":
        command += ["--database", str(database), "--source", str(source)]
    elif key == "curate":
        command += ["--database", str(database)]
    elif key == "chemistry":
        command += ["--database", str(database), "--chemistry-input", str(CHEMISTRY)]
    elif key == "remediation-registry":
        command += ["--database", str(database), "--registry", str(REGISTRY)]
        if fixture:
            command.append("--allow-subset")
    elif key == "validate":
        return command + ["--database", str(database), "--root", str(ROOT)]
    else:
        command += ["--database", str(database)]
    if key not in {"curate", "chemistry", "remediation-registry", "validate"} and fixture:
        command += ["--pdb-id", "3EKY"]
    if key in {"mapping", "sasa", "arpeggio", "geometry", "functional-groups", "protacability", "attachment-sites"} and args.limit:
        command += ["--limit", str(args.limit)]
    if key == "arpeggio":
        command += ["--workers", str(args.workers), "--per-instance-timeout", str(args.per_instance_timeout), "--retry-timeout", str(args.retry_timeout)]
        if args.resume:
            command.append("--resume")
    if key == "mapping":
        command += ["--workers", str(args.workers), "--per-instance-timeout", str(args.per_instance_timeout)]
        if args.resume:
            command.append("--resume")
    return command


def validate_fixture(database: Path) -> None:
    with sqlite3.connect(database) as db:
        dr7 = db.execute("""SELECT i.ligand_instance_id FROM ligand_instances i
            JOIN structures s ON s.structure_id=i.structure_id
            WHERE s.entry_id='3EKY' AND i.label_comp_id='DR7' AND i.auth_asym_id='A' AND i.auth_seq_id='100'""").fetchone()
        if not dr7:
            raise AssertionError("3EKY/DR7 chain A residue 100 is absent")
        instance_id = dr7[0]
        mapped = db.execute("SELECT count(*) FROM ligand_smiles_atom_mapping WHERE ligand_instance_id=? AND method_version='legacy_mcs_etkdg_uff_cif_v2.5' AND mapping_status IN ('mapped_element_validated','complete')", (instance_id,)).fetchone()[0]
        exposed = db.execute("SELECT count(*) FROM ligand_sasa_atoms WHERE ligand_instance_id=? AND method_version='biopython-shrake_rupley-1.40-cif-v2.1' AND legacy_exposed=1", (instance_id,)).fetchone()[0]
        scores = dict(db.execute("SELECT chain_id, protacability_proxy_score FROM protacability_assessment WHERE ligand_instance_id=? AND method_version='protacability-cif-v2.8'", (instance_id,)).fetchall())
        high = db.execute("SELECT count(*) FROM protacability_attachment_sites WHERE ligand_instance_id=? AND method_version='attachment-sites-cif-v2.6' AND high_priority_attachment_atom=1", (instance_id,)).fetchone()[0]
        moderate = {row[0] for row in db.execute("SELECT exact_atom FROM protacability_attachment_sites WHERE ligand_instance_id=? AND method_version='attachment-sites-cif-v2.6' AND attachment_priority_score=79 AND atom_chemical_role='conditional_substitution_site'", (instance_id,))}
    assert mapped == 51, f"mapping {mapped}/51"
    assert exposed == 12, f"exposed atoms {exposed}/12"
    assert scores == {"A": 52.86, "B": 81.43}, scores
    assert high == 0, f"high sites {high}"
    assert moderate == {"CAO", "CAR", "CAS", "NBD"}, moderate
    print("3EKY/DR7 fixture: PASS")


def write_release_artifacts(database: Path) -> None:
    """Fingerprint a Stage-15-passing full reconstruction and compare it to release."""
    fingerprints = OUTPUTS / "SCIENTIFIC_CONTENT_FINGERPRINTS.json"
    command = [sys.executable, str(ROOT / "tools" / "release_fingerprints.py"), "--database", str(database),
               "--output", str(fingerprints), "--reference", str(REFERENCE_FINGERPRINTS)]
    event = {"stage": "release-fingerprints", "key": "fingerprints", "start": now()}
    run_command(command, False, event)
    event["end"] = now()
    checksums = OUTPUTS / "REPRODUCED_RELEASE_SHA256.txt"
    targets = [database, OUTPUTS / "FINAL_RELEASE_VALIDATION_REPORT.md", OUTPUTS / "FINAL_RELEASE_VALIDATION_REPORT.json", fingerprints]
    checksums.write_text("".join(f"{sha256(path)}  {path.name}\n" for path in targets), encoding="utf-8")
    print(f"Release artifacts: {fingerprints.name}, {checksums.name}")


def run_pipeline(args, fixture: bool) -> int:
    source = FIXTURE_SOURCE if fixture else Path(args.source or "PDB_FILES").resolve()
    database = Path(args.database).resolve()
    if not fixture:
        verify_inputs(source, allow_current_upstream=args.allow_current_upstream)
    else:
        verify_inputs(source, fixture=True)
    if not args.dry_run:
        missing = check_environment(require_arpeggio=True)
        if missing:
            raise RuntimeError("Environment check failed before expensive work: " + ", ".join(missing))
        OUTPUTS.mkdir(exist_ok=True)
    events = []
    stage_plan = selected_stages(args.from_stage, args.to_stage)
    # The fixture validates real Stage 01--13 science against its compact
    # corpus.  Stage 15 is intentionally a full-release gate with full-corpus
    # denominators and therefore runs only under --full.
    if fixture and args.to_stage is None:
        stage_plan = stage_plan[:STAGE_INDEX["attachment-sites"] + 1]
    for number, key, script in stage_plan:
        event = {"stage": number, "key": key, "start": now()}
        print(f"\n=== Stage {number}: {key} ===", flush=True)
        command = stage_command(key, script, database, source, args, fixture)
        run_command(command, args.dry_run, event)
        if not args.dry_run and key == "mapping":
            timeouts = current_mapping_timeouts(database)
            for instance_id in timeouts:
                retry = [sys.executable, str(PIPELINE / script), "--database", str(database), "--ligand-instance-id", str(instance_id), "--per-instance-timeout", str(args.mapping_retry_timeout)]
                retry_event = {"stage": "07-retry", "key": key, "ligand_instance_id": instance_id, "start": now()}
                run_command(retry, False, retry_event); events.append(retry_event)
            if current_mapping_timeouts(database):
                raise RuntimeError("Stage 07 mapping_timeout rows remain after deterministic retry")
        if not args.dry_run and key == "arpeggio":
            previous = current_arpeggio_unresolved(database)
            for round_number in range(1, args.arpeggio_recovery_rounds + 1):
                if previous == 0:
                    break
                retry = command + ["--resume"] if "--resume" not in command else command
                retry_event = {"stage": f"09-resume-{round_number}", "key": key, "start": now(), "unresolved_before": previous}
                run_command(retry, False, retry_event); events.append(retry_event)
                current = current_arpeggio_unresolved(database)
                retry_event["unresolved_after"] = current
                if current >= previous:
                    raise RuntimeError(f"Stage 09 recovery made no progress ({previous} unresolved)")
                previous = current
            if previous:
                raise RuntimeError(f"Stage 09 unresolved included cases remain: {previous}")
        event["end"] = now(); events.append(event)
    if fixture and not args.dry_run:
        validate_fixture(database)
    if not fixture and not args.dry_run and stage_plan and stage_plan[-1][1] == "validate":
        write_release_artifacts(database)
    if not args.dry_run:
        run_manifest = {"started": events[0].get("start") if events else now(), "ended": now(), "fixture": fixture,
                        "reproduction_mode": "CURRENT_UPSTREAM_RECONSTRUCTION" if args.allow_current_upstream else "FROZEN_RELEASE_REPRODUCTION",
                        "database": str(database), "events": events}
        (OUTPUTS / "RUN_MANIFEST.json").write_text(json.dumps(run_manifest, indent=2) + "\n", encoding="utf-8")
    return 0


def build_web_compat(args) -> int:
    scientific = Path(args.database).resolve()
    if not scientific.exists():
        raise FileNotFoundError(f"Validated scientific database is required for --web-compat: {scientific}")
    target = OUTPUTS / "viral_data_cif_v2_WEB_COMPAT.db"
    commands = [
        [sys.executable, str(ROOT / "web_compat" / "WEB_01_build_ligand_synonyms.py"), "--database", str(target), "--offline"],
        [sys.executable, str(ROOT / "web_compat" / "WEB_02_build_ligand_arp_diagram.py"), "--database", str(target)],
    ]
    print(f"Scientific input (never modified): {scientific}")
    print(f"Website-compatibility copy: {target}")
    if args.dry_run:
        for command in commands:
            print("$ " + " ".join(command))
        return 0
    OUTPUTS.mkdir(exist_ok=True)
    shutil.copy2(scientific, target)
    for command in commands:
        completed = subprocess.run(command, cwd=ROOT, text=True)
        if completed.returncode:
            raise RuntimeError(f"Website compatibility step failed: {' '.join(command)}")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--check-environment", action="store_true")
    mode.add_argument("--download-inputs", action="store_true")
    mode.add_argument("--verify-inputs", action="store_true")
    mode.add_argument("--fixture", action="store_true")
    mode.add_argument("--full", action="store_true")
    parser.add_argument("--source", help="Input root; defaults to ./PDB_FILES for public download/verification")
    parser.add_argument("--database", default=str(OUTPUTS / "viral_data_cif_v2_REPRODUCED_RELEASE.db"))
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--from-stage", choices=STAGE_INDEX)
    parser.add_argument("--to-stage", choices=STAGE_INDEX)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--per-instance-timeout", type=float, default=60)
    parser.add_argument("--mapping-retry-timeout", type=float, default=300)
    parser.add_argument("--retry-timeout", type=float, default=600)
    parser.add_argument("--arpeggio-recovery-rounds", type=int, default=3)
    parser.add_argument("--web-compat", action="store_true")
    parser.add_argument("--allow-current-upstream", action="store_true", help="Permit parsed RCSB revisions whose bytes differ from the frozen release")
    parser.add_argument("--pdb-id", action="append", help="Bounded public-input selector for download/verification; repeatable")
    args = parser.parse_args()
    if args.check_environment:
        return 1 if check_environment() else 0
    source = Path(args.source or "PDB_FILES")
    if args.download_inputs:
        acquire_public_inputs(MANIFEST, source, OUTPUTS, args.workers, allow_current_upstream=args.allow_current_upstream,
                              entry_ids={entry.upper() for entry in args.pdb_id} if args.pdb_id else None,
                              dry_run=args.dry_run)
        return 0
    if args.verify_inputs:
        verify_inputs(source, allow_current_upstream=args.allow_current_upstream,
                      entry_ids={entry.upper() for entry in args.pdb_id} if args.pdb_id else None); return 0
    if args.web_compat:
        return build_web_compat(args)
    if not args.fixture and not args.full:
        parser.error("choose --check-environment, --verify-inputs, --fixture, or --full")
    return run_pipeline(args, fixture=args.fixture)


if __name__ == "__main__":
    raise SystemExit(main())
