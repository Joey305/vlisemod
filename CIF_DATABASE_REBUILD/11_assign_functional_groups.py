#!/usr/bin/env python3
"""Assign SMARTS functional groups to occurrence-resolved ligand atoms.

The legacy V-LiSEMOD functional-group workflow identified SMARTS matches on a
ligand SMILES and then joined those matches back to deposited ligand atoms.  In
this CIF rebuild, stage 07 already supplies the authoritative occurrence-scoped
PDB/mmCIF-to-SMILES atom mapping.  Stage 11 therefore performs SMARTS matching
in the exact ``ligands.smiles`` atom-index namespace used by Stage 07; canonical
SMILES is deliberately not used for atom-index joins.  Every mapped index is
validated for uniqueness and element identity before functional-group rows are
materialized.  No coordinate-nearest or atom-number repair is performed here.

An optional ``functional_groups.txt`` file using ``Name: SMARTS`` lines can be
supplied.  If absent, a built-in library covering the historical group names is
used.  Functional-group labels are interpretive chemical context, not claims of
reactivity or synthetic feasibility.
"""
from __future__ import annotations

import argparse
import csv
import importlib
import json
from collections import defaultdict
from pathlib import Path

from rdkit import Chem, RDLogger

c = importlib.import_module("00_common")
VERSION = "rdkit-smarts-functional-groups-v2.3"
MAPPING_VERSION = "legacy_mcs_etkdg_uff_cif_v2.5"

# Historical names retained where practical.  Patterns are intentionally
# interpretable rather than mutually exclusive; one atom may belong to several
# functional-group matches.
DEFAULT_SMARTS = {
    "Hydroxyl/Alcohol": "[OX2H;!$(O-[c,n,o,s,p])]",
    "Ether": "[OD2]([#6])[#6]",
    "Aldehyde": "[CX3H1](=O)[#6,#1]",
    "Ketone": "[#6][CX3](=O)[#6]",
    "Carboxylic-Acid": "[CX3](=O)[OX2H1]",
    "Ester": "[CX3](=O)[OX2][#6]",
    "Terminal-Methyl-Ester": "[CX3](=O)O[CH3]",
    "Amine": "[NX3;H2,H1,H0;!$(N-C=O);!$(N-S(=O)=O)]",
    "Amide": "[NX3][CX3](=O)",
    "Benzene": "c1ccccc1",
    "Sulfonamide": "[SX4](=O)(=O)[NX3]",
    "Thiol": "[SX2H]",
    "Nitrile": "[CX2]#N",
    "Nitro": "[$([N+](=O)[O-]),$([NX3](=O)=O)]",
    "Azo": "[NX2]=[NX2]",
    "Hydrazine": "[NX3][NX3]",
    "Isonitrile": "[N+]#[C-]",
    "Imide": "[NX3]([CX3](=O))[CX3](=O)",
    "Sulfoxide": "[#16X3](=O)",
    "Sulfone": "[#16X4](=O)(=O)",
    "Phosphoric-Acid": "[PX4](=O)([OX2H])([OX2H])[OX2H]",
    "Sulfenic-Acid": "[SX2][OX2H]",
    "Terminal-Alkyne": "[CX2H]#[CX2]",
    "Aryl-Halide": "[c][F,Cl,Br,I]",
    "Terminal-Alkene": "[CX3H2]=[CX3]",
    "Thioester": "[CX3](=O)[SX2][#6]",
    "Isocyanate": "[NX2]=[CX2]=O",
    "Aldoxime": "[CX3H1]=[NX2][OX2H]",
    "Carbamate": "[NX3][CX3](=O)[OX2][#6]",
    "Isothiocyanate": "[NX2]=[CX2]=S",
    "Phosphonate-Terminal": "[PX4](=O)([OX1-,OX2H])[OX1-,OX2H]",
    "Haloalkane": "[CX4][F,Cl,Br,I]",
    "Diazo": "[$([#6][N+]=[N-]),$([#6]=[N+]=[N-])]",
    "Azide": "[$([NX1-]=[NX2+]=[NX1]),$([NX1]#[NX2+][NX1-])]",
    "Phenol": "[c][OX2H]",
    "Thioether": "[#6][SX2][#6]",
    "Carbamic-Acid": "[NX3][CX3](=O)[OX2H]",
    "Urea": "[NX3][CX3](=O)[NX3]",
    "Amidine": "[NX3][CX3](=[NX2])[#6,#7]",
    "Guanidine": "[NX3][CX3](=[NX2])[NX3]",
    "Sulfonic-Acid": "[SX4](=O)(=O)[OX2H]",
    "Thiol-Ester": "[CX3](=O)[SX2][#6]",
    "Phosphine": "[PX3]",
    "Imine": "[CX3]=[NX2]",
    "Hydrazone": "[NX3][NX2]=[CX3]",
    "Sulfonic-Acid-Ester": "[SX4](=O)(=O)[OX2][#6]",
    "Thiocyanate": "[SX1][CX2]#N",
    "Carbodiimide": "[NX2]=[CX2]=[NX2]",
    "Selenide": "[#6][SeX2][#6]",
    "Acyl-Halide": "[CX3](=O)[F,Cl,Br,I]",
    "Enamine": "[NX3][CX3]=[CX3]",
    "Boronic-Acid": "[BX3]([OX2H])[OX2H]",
    "Pyridine": "n1ccccc1",
    "Pyrrole": "[nH]1cccc1",
    "Oxime": "[CX3]=[NX2][OX2H]",
    "Cyclohexane": "[C;R]1[C;R][C;R][C;R][C;R][C;R]1",
}

DIRECT_HANDLE_GROUPS = {
    "Hydroxyl/Alcohol", "Phenol", "Amine", "Thiol", "Carboxylic-Acid",
    "Terminal-Alkyne", "Terminal-Alkene", "Azide", "Aryl-Halide",
    "Haloalkane", "Boronic-Acid",
}
CONDITIONAL_HANDLE_GROUPS = {
    "Ether", "Ester", "Terminal-Methyl-Ester", "Amide", "Sulfonamide",
    "Nitrile", "Carbamate", "Thioether", "Urea", "Imine", "Pyridine",
    "Pyrrole", "Oxime", "Thioester", "Thiol-Ester", "Phosphonate-Terminal",
}


def tractability_role(group_name: str) -> str:
    if group_name in DIRECT_HANDLE_GROUPS:
        return "direct_handle_context"
    if group_name in CONDITIONAL_HANDLE_GROUPS:
        return "conditional_handle_context"
    return "structural_context_only"


def ensure_schema(database: str) -> None:
    c.create_schema(database)
    with c.dbconn(database) as db:
        db.executescript(
            """
            CREATE TABLE IF NOT EXISTS ligand_functional_group_matches (
                functional_group_match_id INTEGER PRIMARY KEY,
                run_id INTEGER NOT NULL REFERENCES analysis_runs(run_id),
                ligand_instance_id INTEGER NOT NULL REFERENCES ligand_instances(ligand_instance_id),
                ligand_id INTEGER NOT NULL REFERENCES ligands(ligand_id),
                functional_group TEXT NOT NULL,
                occurrence_index INTEGER NOT NULL,
                smarts TEXT NOT NULL,
                smiles_atom_indices_json TEXT NOT NULL,
                tractability_role TEXT NOT NULL,
                method_version TEXT NOT NULL,
                UNIQUE(ligand_instance_id,functional_group,occurrence_index,method_version)
            );
            CREATE TABLE IF NOT EXISTS ligand_functional_group_atoms (
                functional_group_atom_id INTEGER PRIMARY KEY,
                run_id INTEGER NOT NULL REFERENCES analysis_runs(run_id),
                ligand_instance_id INTEGER NOT NULL REFERENCES ligand_instances(ligand_instance_id),
                functional_group_match_id INTEGER NOT NULL REFERENCES ligand_functional_group_matches(functional_group_match_id),
                functional_group TEXT NOT NULL,
                occurrence_index INTEGER NOT NULL,
                smiles_atom_index INTEGER NOT NULL,
                ligand_instance_atom_id INTEGER REFERENCES ligand_instance_atoms(ligand_instance_atom_id),
                atom_site_id TEXT,
                exact_atom TEXT,
                element TEXT,
                tractability_role TEXT NOT NULL,
                mapping_status TEXT NOT NULL,
                method_version TEXT NOT NULL,
                UNIQUE(functional_group_match_id,smiles_atom_index)
            );
            CREATE TABLE IF NOT EXISTS ligand_functional_group_summary (
                functional_group_summary_id INTEGER PRIMARY KEY,
                run_id INTEGER NOT NULL REFERENCES analysis_runs(run_id),
                ligand_instance_id INTEGER NOT NULL REFERENCES ligand_instances(ligand_instance_id),
                functional_group_match_count INTEGER NOT NULL,
                functional_group_type_count INTEGER NOT NULL,
                functional_group_types TEXT NOT NULL,
                direct_handle_group_count INTEGER NOT NULL,
                conditional_handle_group_count INTEGER NOT NULL,
                mapped_functional_group_atom_count INTEGER NOT NULL,
                unmapped_functional_group_atom_count INTEGER NOT NULL,
                library_source TEXT NOT NULL,
                method_version TEXT NOT NULL,
                status TEXT NOT NULL,
                UNIQUE(ligand_instance_id,method_version)
            );
            """
        )


def load_library(path: str | None):
    source = "builtin_historical_name_library"
    raw = dict(DEFAULT_SMARTS)
    if path:
        p = Path(path)
        if not p.exists():
            raise FileNotFoundError(p)
        raw = {}
        for line in p.read_text(encoding="utf8", errors="replace").splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if ":" not in line:
                raise ValueError(f"Invalid functional-group line: {line}")
            name, smarts = line.split(":", 1)
            raw[name.strip()] = smarts.strip()
        source = str(p.resolve())
    compiled = {}
    invalid = []
    for name, smarts in raw.items():
        patt = Chem.MolFromSmarts(smarts)
        if patt is None:
            invalid.append((name, smarts))
        else:
            compiled[name] = (smarts, patt)
    if invalid:
        raise ValueError("Invalid SMARTS definitions: " + "; ".join(f"{n}={s}" for n, s in invalid))
    return compiled, source


def latest_mapping(db, iid: int):
    """Return only the validated Stage-07 v2.5 mapping for one occurrence.

    Do not fall back to older mapping generations.  In particular, an
    occurrence that is pending in the remediation registry receives a v2.5
    ``skipped_pending_remediation`` row and must remain structurally unmapped
    here rather than inheriting an older adapter-limited mapping.
    """
    return db.execute(
        """SELECT run_id,smiles_atom_count,mapping_status,method_version
           FROM ligand_mapping_runs
           WHERE ligand_instance_id=? AND method_version=?
             AND downstream_mapping_eligibility=1
             AND mapping_status IN ('complete','complete_altloc_resolved','partial_ccd_difference')
           ORDER BY run_id DESC LIMIT 1""",
        (iid, MAPPING_VERSION),
    ).fetchone()


def _normalized_element(value: str | None) -> str:
    """Normalize deposited element tokens for comparison with RDKit symbols."""
    token = c.norm(value).upper() if value is not None else ""
    if token in {"D", "T"}:
        return "H"
    return token


def mapping_index(db, iid: int, mol):
    """Build and validate the Stage-07 source-SMILES-index -> deposited-atom map.

    Stage 07 numbers ``smiles_atom_index`` in the atom order created from
    ``ligands.smiles``.  This function therefore validates against the *same*
    source-SMILES RDKit molecule used by this stage.  Canonical SMILES atom
    ordering must never be substituted for this index namespace.

    Partial mappings are allowed: a missing index is handled downstream as an
    unmapped functional-group atom.  Ambiguous mappings, out-of-range indices,
    missing deposited atoms, or element disagreements are hard validation
    failures because silently accepting them would attach chemistry to the
    wrong deposited atom.
    """
    mapping_run = latest_mapping(db, iid)
    if mapping_run is None:
        # Some chemistry-resolved occurrences are intentionally ineligible for
        # downstream atom mapping.  We can still assign chemical SMARTS groups;
        # their deposited-atom links remain explicitly unmapped.
        return {}, None, []

    expected_count = mapping_run["smiles_atom_count"]
    if expected_count is not None and int(expected_count) != mol.GetNumAtoms():
        return {}, mapping_run["run_id"], [
            f"source_smiles_atom_count_mismatch:stage07={expected_count}:stage11={mol.GetNumAtoms()}"
        ]

    rows = db.execute(
        """SELECT m.smiles_atom_index,m.ligand_instance_atom_id,a.atom_site_id,
                  COALESCE(a.auth_atom_id,a.label_atom_id) exact_atom,a.element
           FROM ligand_smiles_atom_mapping m
           LEFT JOIN ligand_instance_atoms a ON a.ligand_instance_atom_id=m.ligand_instance_atom_id
           WHERE m.ligand_instance_id=? AND m.run_id=? AND m.ligand_instance_atom_id IS NOT NULL
           ORDER BY m.smiles_atom_index,m.ligand_instance_atom_id""",
        (iid, mapping_run["run_id"]),
    ).fetchall()

    grouped = defaultdict(list)
    errors = []
    for r in rows:
        if r["smiles_atom_index"] is None:
            continue
        idx = int(r["smiles_atom_index"])
        if idx < 0 or idx >= mol.GetNumAtoms():
            errors.append(f"smiles_atom_index_out_of_range:{idx}")
            continue
        grouped[idx].append(r)

    out = {}
    for idx, candidates in grouped.items():
        # Collapse exact duplicate rows, but never choose between distinct
        # deposited atoms for the same source-SMILES atom index.
        unique = {}
        for r in candidates:
            unique[r["ligand_instance_atom_id"]] = r
        candidates = list(unique.values())
        if len(candidates) != 1:
            atom_ids = ",".join(str(r["ligand_instance_atom_id"]) for r in candidates)
            errors.append(f"ambiguous_smiles_index:{idx}:ligand_instance_atom_ids={atom_ids}")
            continue
        atom = candidates[0]
        if atom["ligand_instance_atom_id"] is None or atom["element"] is None:
            errors.append(f"missing_deposited_atom_for_smiles_index:{idx}")
            continue
        expected = mol.GetAtomWithIdx(idx).GetSymbol().upper()
        observed = _normalized_element(atom["element"])
        if expected != observed:
            errors.append(
                f"element_mismatch:smiles_index={idx}:smiles_element={expected}:"
                f"deposited_element={observed}:ligand_instance_atom_id={atom['ligand_instance_atom_id']}"
            )
            continue
        out[idx] = atom
    return out, mapping_run["run_id"], errors


def annotate_instance(db, iid: int, library):
    row = db.execute(
        """SELECT i.ligand_instance_id,i.ligand_id,l.smiles,l.canonical_smiles,l.chemical_status
           FROM ligand_instances i JOIN ligands l ON l.ligand_id=i.ligand_id
           WHERE i.ligand_instance_id=?""", (iid,)
    ).fetchone()
    if row is None:
        return {"status": "failed", "reason": "unknown_ligand_instance"}

    # CRITICAL: Stage 07's smiles_atom_index namespace comes from l.smiles, not
    # canonical_smiles.  Using canonical_smiles here can preserve connectivity
    # while reordering atoms and therefore attach correct SMARTS matches to the
    # wrong deposited atoms.
    source_smiles = row["smiles"]
    if row["chemical_status"] != "resolved" or not source_smiles:
        return {"status": "skipped", "reason": "missing_or_unresolved_smiles", "ligand_id": row["ligand_id"]}
    mol = Chem.MolFromSmiles(source_smiles)
    if mol is None:
        return {"status": "failed", "reason": "rdkit_invalid_source_smiles", "ligand_id": row["ligand_id"]}

    mapping, mapping_run_id, validation_errors = mapping_index(db, iid, mol)
    if validation_errors:
        return {
            "status": "failed",
            "reason": "mapping_namespace_validation_failed:" + "|".join(validation_errors[:10]),
            "ligand_id": row["ligand_id"],
            "mapping_run_id": mapping_run_id,
        }

    matches = []
    for group_name, (smarts, patt) in library.items():
        found = sorted(set(tuple(int(x) for x in match) for match in mol.GetSubstructMatches(patt, uniquify=True)))
        for occurrence_index, atom_indices in enumerate(found, 1):
            matches.append({
                "group": group_name,
                "occurrence": occurrence_index,
                "smarts": smarts,
                "indices": atom_indices,
                "role": tractability_role(group_name),
            })
    matches.sort(key=lambda x: (x["group"], x["occurrence"], x["indices"]))
    return {
        "status": "complete", "ligand_id": row["ligand_id"], "matches": matches,
        "mapping": mapping, "mapping_run_id": mapping_run_id, "source_smiles": source_smiles,
    }


def _write_report(database: str, library_source: str):
    c.dirs()
    with c.dbconn(database) as db:
        vals = dict(db.execute(
            """SELECT
                 count(*) instances,
                 sum(functional_group_match_count) matches,
                 sum(mapped_functional_group_atom_count) mapped_atoms,
                 sum(unmapped_functional_group_atom_count) unmapped_atoms
               FROM ligand_functional_group_summary WHERE method_version=?""", (VERSION,)
        ).fetchone())
        direct = db.execute(
            "SELECT sum(direct_handle_group_count) FROM ligand_functional_group_summary WHERE method_version=?", (VERSION,)
        ).fetchone()[0] or 0
    lines = ["# Stage 11 functional-group report", "", f"* Method: {VERSION}", f"* Library: {library_source}", "* Atom index space: ligands.smiles (same source-SMILES ordering used by Stage 07)", "* Mapping validation: unique source-SMILES index plus RDKit/deposited element identity", f"* Required Stage-07 mapping: {MAPPING_VERSION}"]
    lines += [f"* {k.replace('_',' ')}: {v or 0}" for k, v in vals.items()]
    lines.append(f"* direct handle group matches: {direct}")
    (c.ROOT / "outputs" / "FUNCTIONAL_GROUP_STAGE_REPORT.md").write_text("\n".join(lines) + "\n")


def run(database: str, limit=None, pdb_id=None, instance_id=None, functional_groups_file=None,
        resume=False, progress_every=250):
    ensure_schema(database); c.dirs(); RDLogger.DisableLog("rdApp.*")
    library, library_source = load_library(functional_groups_file)
    with c.dbconn(database) as db:
        q = """SELECT i.ligand_instance_id FROM ligand_instances i
               JOIN ligands l ON l.ligand_id=i.ligand_id
               JOIN structures s ON s.structure_id=i.structure_id
               WHERE i.curation_status='included' AND l.chemical_status='resolved'"""
        args = []
        if pdb_id:
            q += " AND s.entry_id=?"; args.append(pdb_id)
        if instance_id:
            q += " AND i.ligand_instance_id=?"; args.append(instance_id)
        if resume:
            q += " AND NOT EXISTS (SELECT 1 FROM ligand_functional_group_summary f WHERE f.ligand_instance_id=i.ligand_instance_id AND f.method_version=? AND f.status='complete')"
            args.append(VERSION)
        ids = [r[0] for r in db.execute(q + " ORDER BY i.ligand_instance_id", args)]
        if limit:
            ids = ids[:limit]
        rid = c.run_start(db, "functional_groups", {
            "method": VERSION, "limit": limit, "pdb_id": pdb_id, "ligand_instance_id": instance_id,
            "functional_groups_file": functional_groups_file, "library_source": library_source, "resume": resume,
            "atom_index_space": "ligands.smiles / Stage-07 smiles_atom_index", "mapping_element_validation": True, "mapping_method_version": MAPPING_VERSION,
        })
        if ids and not resume:
            marks = ",".join("?" for _ in ids)
            # Children first because of FK to match rows.
            db.execute(f"DELETE FROM ligand_functional_group_atoms WHERE method_version=? AND ligand_instance_id IN ({marks})", [VERSION, *ids])
            db.execute(f"DELETE FROM ligand_functional_group_matches WHERE method_version=? AND ligand_instance_id IN ({marks})", [VERSION, *ids])
            db.execute(f"DELETE FROM ligand_functional_group_summary WHERE method_version=? AND ligand_instance_id IN ({marks})", [VERSION, *ids])

        success = skipped = failures = 0
        for n, iid in enumerate(ids, 1):
            try:
                result = annotate_instance(db, iid, library)
                if result["status"] == "failed":
                    failures += 1
                    c.fail(db, rid, "functional_groups", result["reason"], instance_id=iid, code=result["reason"])
                    continue
                if result["status"] == "skipped":
                    skipped += 1
                    db.execute(
                        """INSERT OR REPLACE INTO ligand_functional_group_summary(
                             run_id,ligand_instance_id,functional_group_match_count,functional_group_type_count,
                             functional_group_types,direct_handle_group_count,conditional_handle_group_count,
                             mapped_functional_group_atom_count,unmapped_functional_group_atom_count,
                             library_source,method_version,status) VALUES(?,?,?,?,?,?,?,?,?,?,?,?)""",
                        (rid, iid, 0, 0, "", 0, 0, 0, 0, library_source, VERSION, "skipped_unresolved_chemistry")
                    )
                    continue

                mapped_count = unmapped_count = 0
                types = []
                direct_count = conditional_count = 0
                for match in result["matches"]:
                    types.append(match["group"])
                    direct_count += int(match["role"] == "direct_handle_context")
                    conditional_count += int(match["role"] == "conditional_handle_context")
                    cur = db.execute(
                        """INSERT INTO ligand_functional_group_matches(
                             run_id,ligand_instance_id,ligand_id,functional_group,occurrence_index,smarts,
                             smiles_atom_indices_json,tractability_role,method_version)
                           VALUES(?,?,?,?,?,?,?,?,?)""",
                        (rid, iid, result["ligand_id"], match["group"], match["occurrence"], match["smarts"],
                         json.dumps(match["indices"]), match["role"], VERSION)
                    )
                    match_id = cur.lastrowid
                    for smiles_idx in match["indices"]:
                        mapped = result["mapping"].get(int(smiles_idx))
                        if mapped is None:
                            unmapped_count += 1
                            db.execute(
                                """INSERT INTO ligand_functional_group_atoms(
                                     run_id,ligand_instance_id,functional_group_match_id,functional_group,occurrence_index,
                                     smiles_atom_index,ligand_instance_atom_id,atom_site_id,exact_atom,element,
                                     tractability_role,mapping_status,method_version)
                                   VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?)""",
                                (rid, iid, match_id, match["group"], match["occurrence"], smiles_idx, None, None, None, None,
                                 match["role"], "unmapped_smiles_atom", VERSION)
                            )
                        else:
                            atom = mapped
                            mapped_count += 1
                            db.execute(
                                """INSERT INTO ligand_functional_group_atoms(
                                     run_id,ligand_instance_id,functional_group_match_id,functional_group,occurrence_index,
                                     smiles_atom_index,ligand_instance_atom_id,atom_site_id,exact_atom,element,
                                     tractability_role,mapping_status,method_version)
                                   VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?)""",
                                (rid, iid, match_id, match["group"], match["occurrence"], smiles_idx,
                                 atom["ligand_instance_atom_id"], atom["atom_site_id"], atom["exact_atom"], atom["element"],
                                 match["role"], "mapped_element_validated", VERSION)
                            )
                unique_types = sorted(set(types))
                db.execute(
                    """INSERT OR REPLACE INTO ligand_functional_group_summary(
                         run_id,ligand_instance_id,functional_group_match_count,functional_group_type_count,
                         functional_group_types,direct_handle_group_count,conditional_handle_group_count,
                         mapped_functional_group_atom_count,unmapped_functional_group_atom_count,
                         library_source,method_version,status) VALUES(?,?,?,?,?,?,?,?,?,?,?,?)""",
                    (rid, iid, len(result["matches"]), len(unique_types), ";".join(unique_types), direct_count,
                     conditional_count, mapped_count, unmapped_count, library_source, VERSION, "complete")
                )
                success += 1
            except Exception as exc:
                failures += 1
                c.fail(db, rid, "functional_groups", f"{type(exc).__name__}: {exc}", instance_id=iid, code="functional_group_exception")
            if n % max(1, progress_every) == 0 or n == len(ids):
                db.commit(); print(f"functional-group progress: {n}/{len(ids)} success={success} skipped={skipped} failures={failures}", flush=True)
        c.run_end(db, rid, "completed" if failures == 0 else "partial", len(ids), success, skipped, failures)
    _write_report(database, library_source)
    return {"run_id": rid, "processed": len(ids), "success": success, "skipped": skipped, "failures": failures,
            "functional_group_definitions": len(library), "library_source": library_source}


def main():
    p = argparse.ArgumentParser(description="Assign SMARTS functional groups and map them to occurrence-resolved deposited atoms.")
    p.add_argument("--database", default=str(c.ROOT / "viral_data_cif_v2.db"))
    p.add_argument("--limit", type=int)
    p.add_argument("--pdb-id")
    p.add_argument("--ligand-instance-id", type=int)
    p.add_argument("--functional-groups-file")
    p.add_argument("--resume", action="store_true")
    p.add_argument("--progress-every", type=int, default=250)
    a = p.parse_args()
    print(json.dumps(run(a.database, a.limit, a.pdb_id, a.ligand_instance_id,
                         a.functional_groups_file, a.resume, a.progress_every), indent=2))


if __name__ == "__main__":
    main()
