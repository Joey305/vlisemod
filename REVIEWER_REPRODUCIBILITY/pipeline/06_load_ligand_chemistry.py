"""Stage 06: resolve ligand chemistry from a frozen package input only."""
from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from importlib import import_module
from pathlib import Path

from rdkit import Chem

c = import_module("00_common")
VERSION = "local-chemistry-v2.1"
DEFAULT_INPUT = c.ROOT / "inputs" / "chemistry" / "frozen_component_chemistry.csv"


def canonical(smiles: str):
    molecule = Chem.MolFromSmiles(smiles)
    return (
        Chem.MolToSmiles(molecule, canonical=True),
        Chem.MolToInchiKey(molecule),
    ) if molecule else (None, None)


def load_frozen_chemistry(path: Path):
    if not path.exists():
        raise FileNotFoundError(
            f"Frozen chemistry input is required and was not found: {path}. "
            "Supply --chemistry-input; Stage 06 never reads a website database."
        )
    rows = defaultdict(list)
    with path.open(newline="", encoding="utf-8") as handle:
        required = {"component_id", "source_smiles"}
        if not required.issubset(set(handle.readline().strip().split(","))):
            raise ValueError(f"Chemistry input is missing required columns {sorted(required)}: {path}")
        handle.seek(0)
        for row in csv.DictReader(handle):
            component = (row.get("component_id") or "").strip().upper()
            smiles = (row.get("source_smiles") or "").strip()
            if component and smiles:
                rows[component].append((
                    smiles,
                    (row.get("source_name") or "frozen_component_chemistry").strip(),
                    path.name,
                    1,
                ))
    if not rows:
        raise ValueError(f"Frozen chemistry input is empty: {path}")
    return rows


def resolve(database: str, chemistry_input: Path = DEFAULT_INPUT, limit=None, pdb_id=None):
    c.create_schema(database)
    c.dirs()
    sources = load_frozen_chemistry(Path(chemistry_input))
    with c.dbconn(database) as db:
        where = """WHERE EXISTS (
            SELECT 1 FROM ligand_instances i JOIN structures s ON s.structure_id=i.structure_id
            WHERE i.ligand_id=l.ligand_id AND i.curation_status='included'
        )"""
        args = []
        if pdb_id:
            where += " AND EXISTS (SELECT 1 FROM ligand_instances i JOIN structures s ON s.structure_id=i.structure_id WHERE i.ligand_id=l.ligand_id AND s.entry_id=?)"
            args.append(pdb_id)
        ligands = db.execute(f"SELECT ligand_id,component_id FROM ligands l {where} ORDER BY component_id", args).fetchall()
        ligands = ligands[:limit] if limit else ligands
        run_id = c.run_start(db, "chemistry", {"method": VERSION, "chemistry_input": Path(chemistry_input).name, "network": False})
        report, counts = [], defaultdict(int)
        for ligand in ligands:
            ligand_id, component = ligand["ligand_id"], ligand["component_id"]
            candidates = sources.get(component.upper(), [])
            parsed = []
            for smiles, name, locator, priority in candidates:
                canonical_smiles, inchikey = canonical(smiles)
                status = "valid" if canonical_smiles else "invalid"
                db.execute("INSERT OR IGNORE INTO ligand_chemistry_sources(ligand_id,source_name,source_smiles,source_locator,source_priority,parse_status) VALUES(?,?,?,?,?,?)", (ligand_id, name, smiles, locator, priority, status))
                if canonical_smiles:
                    parsed.append((canonical_smiles, inchikey, smiles, name, locator, priority))
            unique = {item[0] for item in parsed}
            if not candidates:
                status, chosen = "missing_smiles", (None,) * 6
            elif not parsed:
                status, chosen = "invalid_smiles", (None,) * 6
            else:
                status = "conflicting_sources" if len(unique) > 1 else "resolved"
                chosen = sorted(parsed, key=lambda item: (item[5], item[3], item[2]))[0]
            canonical_smiles, inchikey, smiles, name, locator, _ = chosen
            db.execute("UPDATE ligands SET smiles=?,smiles_source=?,source_version=?,canonical_smiles=?,inchikey=?,chemical_definition_version=?,chemical_status=? WHERE ligand_id=?", (smiles, name, VERSION, canonical_smiles, inchikey, VERSION, status, ligand_id))
            counts[status] += 1
            report.append({"ligand_id": ligand_id, "component_id": component, "status": status, "source_smiles": smiles or "", "canonical_smiles": canonical_smiles or "", "smiles_source": name or "", "source_locator": locator or "", "source_count": len(candidates)})
        c.run_end(db, run_id, "completed", len(ligands), len(ligands) - counts["missing_smiles"] - counts["invalid_smiles"], counts["conflicting_sources"], 0)
    output = c.ROOT / "outputs" / "CHEMISTRY_RESOLUTION.csv"
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(report[0]) if report else ["ligand_id"])
        writer.writeheader(); writer.writerows(report)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--database", default=str(c.ROOT / "outputs" / "viral_data_cif_v2_reproduced.db"))
    parser.add_argument("--chemistry-input", default=str(DEFAULT_INPUT))
    parser.add_argument("--limit", type=int)
    parser.add_argument("--pdb-id")
    args = parser.parse_args()
    print(f"resolved {len(resolve(args.database, Path(args.chemistry_input), args.limit, args.pdb_id))} identities")
