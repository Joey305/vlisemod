#!/usr/bin/env python3
"""Author-only utility to derive small frozen reviewer inputs from a release DB.

This utility is not used by reproduction.  It makes the provenance of the
checked-in chemistry and remediation inputs explicit instead of requiring the
website database at reviewer runtime.
"""
from __future__ import annotations

import argparse
import csv
import sqlite3
from pathlib import Path


def write_csv(path: Path, columns: list[str], rows) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--database", required=True)
    parser.add_argument("--package-root", required=True)
    args = parser.parse_args()
    root = Path(args.package_root).resolve()
    db = sqlite3.connect(Path(args.database).resolve().as_uri() + "?mode=ro", uri=True)
    db.row_factory = sqlite3.Row

    chemistry = db.execute(
        """
        SELECT component_id, smiles AS source_smiles,
               COALESCE(smiles_source, 'frozen_final_release') AS source_name,
               canonical_smiles, chemical_status
        FROM ligands
        WHERE smiles IS NOT NULL AND TRIM(smiles) <> ''
        ORDER BY component_id
        """
    ).fetchall()
    write_csv(
        root / "inputs/chemistry/frozen_component_chemistry.csv",
        ["component_id", "source_smiles", "source_name", "canonical_smiles", "chemical_status"],
        (dict(row) for row in chemistry),
    )

    registry = db.execute(
        """
        SELECT s.entry_id, i.deposited_model_num, i.label_asym_id,
               i.label_comp_id AS component_id, i.auth_asym_id, i.auth_seq_id,
               i.insertion_code_normalized, q.mapping_reason_class,
               q.algorithm_queue, q.preflight_mapping_status,
               q.heavy_atoms_structural, q.heavy_atoms_reference,
               q.heavy_atoms_mapped, q.heavy_atom_mapping_fraction,
               q.remediation_status, q.notes
        FROM mapping_remediation_queue q
        JOIN ligand_instances i ON i.ligand_instance_id = q.ligand_instance_id
        JOIN structures s ON s.structure_id = i.structure_id
        WHERE q.remediation_status = 'pending'
        ORDER BY s.entry_id, i.deposited_model_num, i.label_asym_id,
                 i.label_comp_id, i.auth_asym_id, i.auth_seq_id,
                 i.insertion_code_normalized
        """
    ).fetchall()
    columns = list(registry[0].keys()) if registry else []
    write_csv(root / "inputs/mapping/frozen_mapping_remediation_registry.csv", columns, (dict(row) for row in registry))
    print(f"chemistry rows={len(chemistry)} remediation rows={len(registry)}")


if __name__ == "__main__":
    main()
