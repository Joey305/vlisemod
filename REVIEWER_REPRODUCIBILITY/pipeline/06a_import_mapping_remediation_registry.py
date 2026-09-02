"""Import the frozen, occurrence-keyed Stage-07 remediation registry."""
from __future__ import annotations

import argparse
import csv
from collections import Counter
from importlib import import_module
from pathlib import Path

c = import_module("00_common")
DEFAULT_REGISTRY = c.ROOT / "inputs" / "mapping" / "frozen_mapping_remediation_registry.csv"
EXPECTED = {
    ("P11_OTHER", "QUEUE_BOND_TEMPLATE"): 561,
    ("P1_EXACT_HEAVY_TOPOLOGY_MCS_PARTIAL", "QUEUE_EXACT_GRAPH"): 99,
    ("P3_SMALL_ATOMSET_DIFFERENCE", "QUEUE_LOCAL_RECOVERY"): 30,
    ("P7_MCS_TIMEOUT", "QUEUE_COMPLEX_GRAPH"): 22,
    ("P8_NO_MCS_MAPPING", "QUEUE_FRAGMENT_REVIEW"): 5,
    ("P5_LARGE_GRAPH_DIFFERENCE", "QUEUE_CURATION_REVIEW"): 2,
}


def load_rows(path: Path):
    if not path.exists():
        raise FileNotFoundError(f"Frozen remediation registry is required: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    counts = Counter((row["mapping_reason_class"], row["algorithm_queue"]) for row in rows)
    if len(rows) != 719 or counts != Counter(EXPECTED):
        raise ValueError(f"Invalid frozen remediation registry: rows={len(rows)} distribution={dict(counts)}")
    return rows


def import_registry(database: str, registry: Path = DEFAULT_REGISTRY, allow_subset: bool = False):
    rows = load_rows(Path(registry))
    c.create_schema(database)
    imported, skipped = 0, 0
    with c.dbconn(database) as db:
        run_id = c.run_start(db, "mapping_remediation_registry", {"registry": Path(registry).name, "frozen_rows": len(rows)})
        for row in rows:
            match = db.execute(
                """
                SELECT i.ligand_instance_id, s.entry_id
                FROM ligand_instances i JOIN structures s ON s.structure_id=i.structure_id
                WHERE s.entry_id=? AND i.deposited_model_num=? AND i.label_asym_id=?
                  AND i.label_comp_id=? AND i.auth_asym_id=? AND i.auth_seq_id=?
                  AND i.insertion_code_normalized=?
                """,
                (row["entry_id"], row["deposited_model_num"], row["label_asym_id"], row["component_id"], row["auth_asym_id"], row["auth_seq_id"], row["insertion_code_normalized"]),
            ).fetchall()
            if not match:
                skipped += 1
                continue
            if len(match) != 1:
                raise RuntimeError(f"Non-unique frozen remediation occurrence: {row}")
            instance_id = match[0]["ligand_instance_id"]
            db.execute(
                """
                INSERT INTO mapping_remediation_queue(
                  ligand_instance_id,pdb_id,component_id,mapping_reason_class,algorithm_queue,
                  preflight_mapping_status,heavy_atoms_structural,heavy_atoms_reference,
                  heavy_atoms_mapped,heavy_atom_mapping_fraction,remediation_status,
                  first_identified_run_id,notes)
                VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?)
                ON CONFLICT(ligand_instance_id) DO UPDATE SET
                  mapping_reason_class=excluded.mapping_reason_class,
                  algorithm_queue=excluded.algorithm_queue,
                  preflight_mapping_status=excluded.preflight_mapping_status,
                  remediation_status=excluded.remediation_status,
                  notes=excluded.notes
                """,
                (instance_id, row["entry_id"], row["component_id"], row["mapping_reason_class"], row["algorithm_queue"], row["preflight_mapping_status"], row["heavy_atoms_structural"] or None, row["heavy_atoms_reference"] or None, row["heavy_atoms_mapped"] or None, row["heavy_atom_mapping_fraction"] or None, "pending", run_id, row["notes"]),
            )
            imported += 1
        if skipped and not allow_subset:
            raise RuntimeError(f"Frozen remediation registry did not resolve {skipped} occurrences; refusing a full Stage-07 run")
        c.run_end(db, run_id, "completed", imported, imported, skipped, 0)
    return {"frozen_rows": len(rows), "imported": imported, "skipped_not_in_subset": skipped}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--database", required=True)
    parser.add_argument("--registry", default=str(DEFAULT_REGISTRY))
    parser.add_argument("--allow-subset", action="store_true")
    args = parser.parse_args()
    print(import_registry(args.database, Path(args.registry), args.allow_subset))
