"""Occurrence-resolved, provenance-preserving pdbe-arpeggio runner.

Canonical mmCIF files are immutable inputs.  When pdbe-arpeggio cannot consume
one directly, this stage writes a deterministic analysis-only mmCIF plus an
atom-level source map under ``outputs/arpeggio/raw/<run_id>/``.  Completed runs
from older pipeline versions remain valid and are skipped by ``--resume``.
"""
from __future__ import annotations

import argparse
import importlib.metadata
import json
import math
import shutil
import subprocess
import time
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from functools import lru_cache
from importlib import import_module
from pathlib import Path

c = import_module("00_common")
VERSION = "arpeggio-cif-v2.2"
DEFAULT_RETRY_TIMEOUT = 600.0
DEFAULT_CONTEXT_RADIUS = 12.0


class DerivedInputError(RuntimeError):
    """A derived input could not be created without ambiguous provenance."""


class OutputValidationError(RuntimeError):
    """Arpeggio output exists but is not safe to import."""


def executable_provenance():
    exe = shutil.which("pdbe-arpeggio")
    try:
        version = importlib.metadata.version("pdbe-arpeggio")
    except Exception:
        version = "unknown"
    return exe, version


def log_text(value):
    return value.decode(errors="replace") if isinstance(value, bytes) else (value or "")


def occurrence(db, iid):
    return db.execute(
        """SELECT i.*,s.entry_id,s.source_cif_path,s.source_cif_sha256
           FROM ligand_instances i JOIN structures s ON s.structure_id=i.structure_id
           WHERE i.ligand_instance_id=?""",
        (iid,),
    ).fetchone()


def selector_guard(db, row):
    count = db.execute(
        """SELECT count(*) FROM ligand_instances i
           WHERE i.structure_id=? AND i.deposited_model_num=?
             AND i.auth_asym_id=? AND i.auth_seq_id=?
             AND i.insertion_code_normalized=?""",
        (
            row["structure_id"],
            row["deposited_model_num"],
            row["auth_asym_id"],
            row["auth_seq_id"],
            row["insertion_code_normalized"],
        ),
    ).fetchone()[0]
    if count != 1:
        return None, "selector_not_occurrence_unique"
    return selector_for(row), None


def selector_for(row):
    return f'/{row["auth_asym_id"]}/{row["auth_seq_id"]}{row["insertion_code_normalized"]}/'


def atom_name(value):
    """Remove one genuine CIF quoting pair while preserving a chemical prime."""
    value = (value or "").strip()
    return value[1:-1] if len(value) >= 2 and value[0] in "\"'" and value[-1] == value[0] else value


def _token(value, missing="?"):
    value = c.raw(value)
    try:
        value = c.gemmi.cif.as_string(value)
    except Exception:
        value = atom_name(value)
    return c.gemmi.cif.quote(value) if value else missing


def _residue_key(atom):
    return (
        c.norm(atom["_atom_site.label_asym_id"]),
        c.norm(atom["_atom_site.label_comp_id"]),
        c.norm(atom["_atom_site.label_seq_id"]),
        c.norm(atom["_atom_site.auth_asym_id"]),
        c.norm(atom["_atom_site.auth_seq_id"]),
        c.norm(atom["_atom_site.pdbx_PDB_ins_code"]),
    )


def _auth_residue_key(atom):
    return (
        c.norm(atom["_atom_site.auth_asym_id"]),
        c.norm(atom["_atom_site.auth_seq_id"]),
        c.norm(atom["_atom_site.pdbx_PDB_ins_code"]),
        c.norm(atom["_atom_site.auth_comp_id"]),
    )


def _xyz(atom):
    values = tuple(c.fnum(atom[tag]) for tag in ("_atom_site.Cartn_x", "_atom_site.Cartn_y", "_atom_site.Cartn_z"))
    return values if all(value is not None and math.isfinite(value) for value in values) else None


def _choose_conformer(rows, target_source_ids):
    """Keep blank/shared atoms plus one named altloc, chosen residue-wide."""
    if any(c.norm(row["_atom_site.id"]) in target_source_ids for row in rows):
        return [row for row in rows if c.norm(row["_atom_site.id"]) in target_source_ids]
    shared = [row for row in rows if not c.norm(row["_atom_site.label_alt_id"])]
    named = defaultdict(list)
    for row in rows:
        altloc = c.norm(row["_atom_site.label_alt_id"])
        if altloc:
            named[altloc].append(row)
    if not named:
        return shared
    scores = []
    for altloc, candidates in named.items():
        heavy = sum(c.norm(row["_atom_site.type_symbol"]).upper() not in {"H", "D", "T"} for row in candidates)
        occupancy = sum(c.fnum(row["_atom_site.occupancy"]) or 0.0 for row in candidates)
        scores.append((-heavy, -occupancy, altloc))
    chosen = sorted(scores)[0][2]
    return shared + named[chosen]


def _within_context(rows, target_rows, radius):
    target_xyz = [point for point in (_xyz(row) for row in target_rows) if point]
    if not target_xyz:
        raise DerivedInputError("target_missing_finite_coordinates")
    radius_squared = radius * radius
    for row in rows:
        point = _xyz(row)
        if point and any(sum((point[i] - target[i]) ** 2 for i in range(3)) <= radius_squared for target in target_xyz):
            return True
    return False


@lru_cache(maxsize=1)
def _source_payload(source_path, source_sha256):
    """Cache one source structure while consecutive occurrences are prepared."""
    path = Path(source_path)
    block = c.cif_doc(path).sole_block()
    atoms = c.loop_rows(block, c.ATOM_TAGS)
    comp_tags = ["_chem_comp.id", "_chem_comp.type", "_chem_comp.name"]
    return atoms, c.loop_rows(block, comp_tags)


ATOM_OUTPUT_TAGS = [
    "_atom_site.group_PDB",
    "_atom_site.id",
    "_atom_site.type_symbol",
    "_atom_site.label_atom_id",
    "_atom_site.label_comp_id",
    "_atom_site.label_asym_id",
    "_atom_site.label_seq_id",
    "_atom_site.Cartn_x",
    "_atom_site.Cartn_y",
    "_atom_site.Cartn_z",
    "_atom_site.occupancy",
    "_atom_site.B_iso_or_equiv",
    "_atom_site.pdbx_PDB_model_num",
    "_atom_site.auth_atom_id",
    "_atom_site.auth_comp_id",
    "_atom_site.auth_asym_id",
    "_atom_site.auth_seq_id",
    "_atom_site.pdbx_PDB_ins_code",
    "_atom_site.label_alt_id",
]


def build_derived_input(db, row, path, strategy="sanitized_full", context_radius=DEFAULT_CONTEXT_RADIUS):
    """Build a deterministic mmCIF and return its complete retained-atom map.

    The derived file contains one canonical deposited model, one coherent
    conformer per residue, sequential atom IDs, and blank derived altlocs.  A
    pocket strategy additionally retains whole residues having any atom within
    ``context_radius`` Angstrom of a canonical target atom.
    """
    source_path = Path(row["source_cif_path"])
    source_atoms, comps = _source_payload(str(source_path), row["source_cif_sha256"])
    model = str(row["deposited_model_num"])
    selected = {
        str(record["atom_site_id"]): record["ligand_instance_atom_id"]
        for record in db.execute(
            """SELECT ligand_instance_atom_id,atom_site_id FROM ligand_instance_atoms
               WHERE ligand_instance_id=? AND selected_conformer=1""",
            (row["ligand_instance_id"],),
        )
    }
    if not selected:
        raise DerivedInputError("target_has_no_selected_conformer_atoms")

    model_rows = [atom for atom in source_atoms if c.norm(atom["_atom_site.pdbx_PDB_model_num"]) == model]
    residue_groups = defaultdict(list)
    for atom in model_rows:
        residue_groups[_residue_key(atom)].append(atom)
    retained_groups = {key: _choose_conformer(atoms, set(selected)) for key, atoms in residue_groups.items()}
    target_rows = [atom for atoms in retained_groups.values() for atom in atoms if c.norm(atom["_atom_site.id"]) in selected]
    if {c.norm(atom["_atom_site.id"]) for atom in target_rows} != set(selected):
        raise DerivedInputError("target_conformer_atom_set_not_preserved")

    operations = {
        "single_model_selected": model,
        "coherent_altloc_selection": True,
        "derived_altloc_normalized_to_blank": True,
        "atom_ids_renumbered_sequentially": True,
        "source_atom_count_in_model": len(model_rows),
        "context_radius_angstrom": context_radius if strategy == "sanitized_pocket" else None,
    }
    if strategy == "sanitized_pocket":
        retained_groups = {
            key: atoms for key, atoms in retained_groups.items() if _within_context(atoms, target_rows, context_radius)
        }
        operations["whole_residue_context_crop"] = True
    else:
        operations["whole_residue_context_crop"] = False

    # BioPython keys residues by the author namespace.  When distinct label
    # residues collapse onto one author residue, retain the target identity and
    # assign deterministic derived-only chain aliases to the other residues.
    auth_collisions = defaultdict(list)
    for key, atoms in retained_groups.items():
        if atoms:
            auth_collisions[_auth_residue_key(atoms[0])].append(key)
    target_keys = {key for key, atoms in retained_groups.items() if any(c.norm(a["_atom_site.id"]) in selected for a in atoms)}
    used_chains = {c.norm(atom["_atom_site.auth_asym_id"]) for atoms in retained_groups.values() for atom in atoms}
    aliases = {}
    alias_number = 1
    for keys in auth_collisions.values():
        if len(keys) < 2:
            continue
        ordered = sorted(keys, key=lambda key: (key not in target_keys, key))
        for key in ordered[1:]:
            while f"AR{alias_number}" in used_chains:
                alias_number += 1
            aliases[key] = f"AR{alias_number}"
            used_chains.add(aliases[key])
            alias_number += 1
    operations["author_residue_chain_aliases"] = len(aliases)

    prepared = []
    missing_coordinates = 0
    exact_duplicates = 0
    seen_atom_keys = {}
    for residue_key in sorted(retained_groups):
        for atom in sorted(retained_groups[residue_key], key=lambda item: (c.norm(item["_atom_site.id"]), c.norm(item["_atom_site.label_atom_id"]))):
            if _xyz(atom) is None:
                missing_coordinates += 1
                continue
            derived = dict(atom)
            if residue_key in aliases:
                derived["_atom_site.auth_asym_id"] = aliases[residue_key]
            derived["_atom_site.label_alt_id"] = "."
            derived["_atom_site.pdbx_PDB_model_num"] = "1"
            identity_key = (
                c.norm(derived["_atom_site.auth_asym_id"]),
                c.norm(derived["_atom_site.auth_seq_id"]),
                c.norm(derived["_atom_site.pdbx_PDB_ins_code"]),
                c.norm(derived["_atom_site.auth_comp_id"]),
                atom_name(derived["_atom_site.auth_atom_id"]),
            )
            signature = (_xyz(derived), c.norm(derived["_atom_site.type_symbol"]).upper())
            if identity_key in seen_atom_keys:
                previous_signature, previous_source_id = seen_atom_keys[identity_key]
                source_id = c.norm(atom["_atom_site.id"])
                if signature == previous_signature and source_id not in selected and previous_source_id not in selected:
                    exact_duplicates += 1
                    continue
                raise DerivedInputError(f"ambiguous_derived_atom_identity:{identity_key}")
            seen_atom_keys[identity_key] = (signature, c.norm(atom["_atom_site.id"]))
            prepared.append((atom, derived, residue_key in aliases))

    operations["missing_coordinate_atoms_excluded"] = missing_coordinates
    operations["exact_duplicate_atoms_excluded"] = exact_duplicates
    operations["derived_atom_count"] = len(prepared)
    target_retained = {c.norm(source["_atom_site.id"]) for source, _, _ in prepared if c.norm(source["_atom_site.id"]) in selected}
    if target_retained != set(selected):
        raise DerivedInputError("derived_target_atom_set_incomplete")

    lines = [f"data_{c.norm(row['entry_id']) or 'arpeggio_derived'}", "#", f"_entry.id {_token(row['entry_id'])}", "#"]
    comp_tags = ["_chem_comp.id", "_chem_comp.type", "_chem_comp.name"]
    if comps:
        lines.extend(["loop_", *comp_tags])
        lines.extend(" ".join(_token(comp[tag]) for tag in comp_tags) for comp in comps)
        lines.append("#")
    lines.extend(["loop_", *ATOM_OUTPUT_TAGS])
    atom_map = []
    for serial, (source, derived, aliased) in enumerate(prepared, 1):
        derived["_atom_site.id"] = str(serial)
        values = []
        for tag in ATOM_OUTPUT_TAGS:
            missing = "." if tag in {"_atom_site.label_alt_id", "_atom_site.pdbx_PDB_ins_code"} else "?"
            values.append(_token(derived[tag], missing))
        lines.append(" ".join(values))
        source_id = c.norm(source["_atom_site.id"])
        atom_map.append(
            {
                "derived_atom_id": str(serial),
                "source_atom_site_id": source_id,
                "ligand_instance_atom_id": selected.get(source_id),
                "source_model_num": model,
                "source_label_asym_id": c.norm(source["_atom_site.label_asym_id"]),
                "source_auth_asym_id": c.norm(source["_atom_site.auth_asym_id"]),
                "source_label_seq_id": c.norm(source["_atom_site.label_seq_id"]),
                "source_auth_seq_id": c.norm(source["_atom_site.auth_seq_id"]),
                "source_component_id": c.norm(source["_atom_site.label_comp_id"]),
                "source_label_atom_id": atom_name(source["_atom_site.label_atom_id"]),
                "source_auth_atom_id": atom_name(source["_atom_site.auth_atom_id"]),
                "source_element": c.norm(source["_atom_site.type_symbol"]),
                "source_altloc": c.norm(source["_atom_site.label_alt_id"]),
                "source_insertion_code": c.norm(source["_atom_site.pdbx_PDB_ins_code"]),
                "derived_model_num": "1",
                "derived_label_asym_id": c.norm(derived["_atom_site.label_asym_id"]),
                "derived_auth_asym_id": c.norm(derived["_atom_site.auth_asym_id"]),
                "derived_label_seq_id": c.norm(derived["_atom_site.label_seq_id"]),
                "derived_auth_seq_id": c.norm(derived["_atom_site.auth_seq_id"]),
                "derived_component_id": c.norm(derived["_atom_site.label_comp_id"]),
                "derived_label_atom_id": atom_name(derived["_atom_site.label_atom_id"]),
                "derived_auth_atom_id": atom_name(derived["_atom_site.auth_atom_id"]),
                "derived_element": c.norm(derived["_atom_site.type_symbol"]),
                "derived_altloc": "",
                "derived_insertion_code": c.norm(derived["_atom_site.pdbx_PDB_ins_code"]),
                "mapping_status": "exact_with_derived_chain_alias" if aliased else "exact",
            }
        )
    lines.append("#")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")
    # Parse our own artifact before allowing it to reach external software.
    c.cif_doc(path)
    manifest = {
        "format_version": 1,
        "strategy": strategy,
        "authoritative_source_path": str(source_path),
        "authoritative_source_sha256": row["source_cif_sha256"],
        "ligand_instance_id": row["ligand_instance_id"],
        "pdb_id": row["entry_id"],
        "selector": selector_for(row),
        "operations": operations,
        "derived_input_path": str(path),
        "derived_input_sha256": c.sha256(path),
        "atom_map_path": str(path.with_name("atom_provenance.jsonl")),
        "provenance_validation_status": "target_exact_all_retained_atoms_source_mapped",
    }
    with path.with_name("atom_provenance.jsonl").open("w") as handle:
        for mapping in atom_map:
            handle.write(json.dumps(mapping, sort_keys=True) + "\n")
    path.with_name("derivation_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return path, atom_map, operations, manifest


def derived_cif(db, row, path):
    """Backward-compatible wrapper used by older fixtures and callers."""
    return build_derived_input(db, row, path, "canonical_derived")[0]


def persist_derived_map(db, run_id, strategy, atom_map):
    columns = [
        "derived_atom_id", "source_atom_site_id", "ligand_instance_atom_id", "source_model_num",
        "source_label_asym_id", "source_auth_asym_id", "source_label_seq_id", "source_auth_seq_id",
        "source_component_id", "source_label_atom_id", "source_auth_atom_id", "source_element",
        "source_altloc", "source_insertion_code", "derived_model_num", "derived_label_asym_id",
        "derived_auth_asym_id", "derived_label_seq_id", "derived_auth_seq_id", "derived_component_id",
        "derived_label_atom_id", "derived_auth_atom_id", "derived_element", "derived_altloc",
        "derived_insertion_code", "mapping_status",
    ]
    sql = f"INSERT INTO arpeggio_derived_atom_map(run_id,input_strategy,{','.join(columns)}) VALUES({','.join('?' for _ in range(len(columns)+2))})"
    db.executemany(sql, [(run_id, strategy, *(mapping[column] for column in columns)) for mapping in atom_map])


def classify_failure(stderr_text, status=None, fallback="arpeggio_subprocess_failed"):
    text = stderr_text or ""
    if status == "timed_out":
        return "timeout"
    patterns = (
        ("unterminated 'string'", "derived_cif_quoting_error"),
        ("defined twice in residue", "duplicate_atom_identity"),
        ("OBBioMatchError", "openbabel_biopython_mapping"),
        ("could not be matched to a BioPython counterpart", "openbabel_biopython_mapping"),
        ("PDBConstructionException", "biopython_parse_failure"),
        ("read_mmcif_to_biopython", "biopython_parse_failure"),
        ("OpenBabel", "openbabel_parse_failure"),
        ("Open Babel", "openbabel_parse_failure"),
        ("No such file", "missing_input_or_output"),
    )
    for needle, category in patterns:
        if needle in text:
            return category
    return fallback


RECOVERABLE_INPUT_FAILURES = {
    "derived_cif_quoting_error",
    "duplicate_atom_identity",
    "openbabel_biopython_mapping",
    "biopython_parse_failure",
    "openbabel_parse_failure",
    "arpeggio_subprocess_failed",
    "missing_or_invalid_output",
}


def _endpoint_is_ligand(endpoint, row):
    if not endpoint:
        return False
    insertion = c.norm(endpoint.get("pdbx_PDB_ins_code"))
    return (
        str(endpoint.get("auth_asym_id", "")) == str(row["auth_asym_id"])
        and str(endpoint.get("auth_seq_id", "")) == str(row["auth_seq_id"])
        and str(endpoint.get("label_comp_id", "")) == str(row["label_comp_id"])
        and insertion == str(row["insertion_code_normalized"])
    )


def validate_output_file(path, row):
    if not path.exists():
        raise OutputValidationError("expected_json_missing")
    try:
        raw = json.loads(path.read_text())
    except Exception as exc:
        raise OutputValidationError(f"invalid_json:{type(exc).__name__}:{exc}") from exc
    if not isinstance(raw, list):
        raise OutputValidationError("json_root_not_list")
    malformed = sum(not isinstance(item, dict) or "bgn" not in item or "end" not in item for item in raw)
    if malformed:
        raise OutputValidationError(f"malformed_contact_records:{malformed}")
    selected_endpoints = sum(
        _endpoint_is_ligand(item.get(side), row) for item in raw for side in ("bgn", "end")
    )
    if raw and not selected_endpoints:
        raise OutputValidationError("selected_ligand_absent_from_contacts")
    return raw, {"contact_records": len(raw), "selected_ligand_endpoints": selected_endpoints}


def _endpoint_key(endpoint):
    return (
        str(endpoint.get("auth_asym_id", "")),
        str(endpoint.get("auth_seq_id", "")),
        c.norm(endpoint.get("pdbx_PDB_ins_code")),
        str(endpoint.get("label_comp_id", "")),
        atom_name(endpoint.get("auth_atom_id") or endpoint.get("label_atom_id")),
    )


def reconcile(db, iid, endpoint, row=None):
    if not endpoint:
        return None, "unmatched"
    raw_name = (endpoint.get("auth_atom_id") or endpoint.get("label_atom_id") or "").strip()
    name = atom_name(raw_name)
    if not name:
        return None, "unmatched"
    atoms = db.execute(
        """SELECT ligand_instance_atom_id,auth_atom_id,label_atom_id
           FROM ligand_instance_atoms WHERE ligand_instance_id=? AND selected_conformer=1""",
        (iid,),
    ).fetchall()
    if "," in name:
        grouped_names = [part.strip() for part in name.split(",") if part.strip()]
        grouped_matches = []
        for grouped_name in grouped_names:
            ids = {
                atom[0] for atom in atoms
                if atom_name(atom["auth_atom_id"]) == grouped_name or atom_name(atom["label_atom_id"]) == grouped_name
            }
            if len(ids) != 1:
                return None, "ambiguous" if ids else "unmatched"
            grouped_matches.extend(ids)
        return (None, "grouped") if grouped_matches else (None, "unmatched")
    matches = list(
        {
            atom[0]
            for atom in atoms
            if (atom["auth_atom_id"] or "").strip() == raw_name or (atom["label_atom_id"] or "").strip() == raw_name
        }
    )
    if not matches:
        matches = list(
            {
                atom[0]
                for atom in atoms
                if atom_name(atom["auth_atom_id"]) == name or atom_name(atom["label_atom_id"]) == name
            }
        )
    if len(matches) == 1:
        return matches[0], None
    return None, "ambiguous" if matches else "unmatched"


def _derived_endpoint_map(db, run_id, strategy):
    mapped = defaultdict(list)
    if not strategy or strategy == "direct_original":
        return mapped
    for row in db.execute(
        """SELECT source_atom_site_id,derived_auth_asym_id,derived_auth_seq_id,
                  derived_insertion_code,derived_component_id,derived_auth_atom_id
           FROM arpeggio_derived_atom_map WHERE run_id=? AND input_strategy=?""",
        (run_id, strategy),
    ):
        key = (
            str(row["derived_auth_asym_id"]), str(row["derived_auth_seq_id"]),
            c.norm(row["derived_insertion_code"]), str(row["derived_component_id"]),
            atom_name(row["derived_auth_atom_id"]),
        )
        mapped[key].append(row["source_atom_site_id"])
    return mapped


def parse_contacts(db, run_id, iid, path, row, strategy=None):
    raw, _ = validate_output_file(path, row)
    pairs = Counter()
    endpoint_stats = Counter()
    parsed = []
    provenance = _derived_endpoint_map(db, run_id, strategy)
    for index, contact in enumerate(raw):
        beginning, end = contact.get("bgn", {}), contact.get("end", {})
        ligand_b = _endpoint_is_ligand(beginning, row)
        ligand_e = _endpoint_is_ligand(end, row)
        atom_b, reason_b = reconcile(db, iid, beginning) if ligand_b else (None, None)
        atom_e, reason_e = reconcile(db, iid, end) if ligand_e else (None, None)
        for selected_endpoint, atom_id, reason in ((ligand_b, atom_b, reason_b), (ligand_e, atom_e, reason_e)):
            if selected_endpoint:
                endpoint_stats["observed"] += 1
                endpoint_stats["reconciled" if atom_id else reason] += 1
        atom_id = atom_b or atom_e
        partner = end if ligand_b else beginning
        component = partner.get("label_comp_id", "")
        filter_class = (
            "raw_environment" if atom_id and component not in c.WATER
            else "raw_water" if atom_id and component in c.WATER
            else "raw_environment_grouped_ligand" if (ligand_b or ligand_e) and (reason_b == "grouped" or reason_e == "grouped")
            else "not_selected_ligand"
        )
        partner_key = json.dumps(partner, sort_keys=True)
        source_matches = provenance.get(_endpoint_key(partner), [])
        source_atom = source_matches[0] if len(source_matches) == 1 else None
        partner_mapping = (
            "exact" if len(source_matches) == 1 else "ambiguous" if source_matches
            else "not_available_direct_input" if strategy == "direct_original" else "unmatched"
        )
        parsed.append((index, contact, atom_id, partner, partner_key, filter_class, source_atom, partner_mapping))

    if endpoint_stats["observed"] and endpoint_stats["reconciled"] + endpoint_stats["grouped"] != endpoint_stats["observed"]:
        raise OutputValidationError(
            f"ligand_endpoint_mapping_invalid:observed={endpoint_stats['observed']};"
            f"reconciled={endpoint_stats['reconciled']};unmatched={endpoint_stats['unmatched']};"
            f"ambiguous={endpoint_stats['ambiguous']}"
        )

    inserted = 0
    for index, contact, atom_id, partner, partner_key, filter_class, source_atom, partner_mapping in parsed:
        labels = contact.get("contact", []) or ["UNLABELLED"]
        for label in labels:
            db.execute(
                """INSERT INTO arpeggio_raw_contact_labels(
                     run_id,ligand_instance_id,raw_contact_index,interaction_label,distance,
                     bgn_json,end_json,ligand_instance_atom_id,partner_identity_json,filter_class,
                     partner_source_atom_site_id,partner_mapping_status)
                   VALUES(?,?,?,?,?,?,?,?,?,?,?,?)""",
                (
                    run_id, iid, index, label, contact.get("distance"),
                    json.dumps(contact.get("bgn", {}), sort_keys=True),
                    json.dumps(contact.get("end", {}), sort_keys=True), atom_id,
                    partner_key, filter_class, source_atom, partner_mapping,
                ),
            )
            inserted += 1
            if atom_id and filter_class == "raw_environment":
                pairs[(atom_id, partner_key)] += 1
    for (atom_id, partner_key), count in pairs.items():
        db.execute(
            """INSERT INTO arpeggio_unique_atom_pairs(
                 run_id,ligand_instance_id,ligand_instance_atom_id,partner_identity_json,raw_label_count)
               VALUES(?,?,?,?,?)""",
            (run_id, iid, atom_id, partner_key, count),
        )
    return inserted, len(pairs), dict(endpoint_stats)


def execute_attempt(plan, row):
    started_wall = c.now()
    started = time.monotonic()
    stdout_path, stderr_path = Path(plan["stdout"]), Path(plan["stderr"])
    stdout_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        process = subprocess.run(plan["command"], text=True, capture_output=True, timeout=plan["timeout"])
        stdout_path.write_text(log_text(process.stdout))
        stderr_path.write_text(log_text(process.stderr))
        elapsed = time.monotonic() - started
        output = Path(plan["output"])
        if process.returncode != 0:
            failure_class = classify_failure(log_text(process.stderr))
            return {**plan, "status": "failed", "returncode": process.returncode, "failure_class": failure_class,
                    "error_class": "ArpeggioSubprocessError", "error_message": log_text(process.stderr)[-2000:],
                    "elapsed": elapsed, "start_time": started_wall, "end_time": c.now(), "raw": None,
                    "output_validation_status": "not_validated", "provenance_validation_status": plan["provenance_status"]}
        try:
            _, validation = validate_output_file(output, row)
        except OutputValidationError as exc:
            return {**plan, "status": "failed", "returncode": process.returncode,
                    "failure_class": "missing_or_invalid_output", "error_class": type(exc).__name__,
                    "error_message": str(exc), "elapsed": elapsed, "start_time": started_wall,
                    "end_time": c.now(), "raw": None, "output_validation_status": str(exc),
                    "provenance_validation_status": plan["provenance_status"]}
        return {**plan, "status": "completed", "returncode": process.returncode, "failure_class": None,
                "error_class": None, "error_message": None, "elapsed": elapsed, "start_time": started_wall,
                "end_time": c.now(), "raw": output, "output_validation_status": "valid_parseable_ligand_selected",
                "provenance_validation_status": plan["provenance_status"], "validation": validation}
    except subprocess.TimeoutExpired as exc:
        stdout_path.write_text(log_text(exc.stdout))
        stderr_path.write_text(log_text(exc.stderr))
        return {**plan, "status": "timed_out", "returncode": None, "failure_class": "timeout",
                "error_class": type(exc).__name__, "error_message": f"timeout_after_{plan['timeout']}_seconds",
                "elapsed": time.monotonic() - started, "start_time": started_wall, "end_time": c.now(),
                "raw": None, "output_validation_status": "not_validated",
                "provenance_validation_status": plan["provenance_status"]}
    except Exception as exc:
        stderr_path.write_text(repr(exc))
        return {**plan, "status": "failed", "returncode": None, "failure_class": classify_failure(str(exc), fallback="runner_exception"),
                "error_class": type(exc).__name__, "error_message": str(exc), "elapsed": time.monotonic() - started,
                "start_time": started_wall, "end_time": c.now(), "raw": None,
                "output_validation_status": "not_validated", "provenance_validation_status": plan["provenance_status"]}


def execute_task(task):
    attempts = []
    for plan in task["plans"]:
        if attempts:
            previous = attempts[-1]
            trigger = plan.get("trigger")
            if trigger == "recoverable_failure" and not (
                previous["status"] == "failed" and previous["failure_class"] in RECOVERABLE_INPUT_FAILURES
            ):
                continue
            if trigger == "timeout" and previous["status"] != "timed_out":
                continue
        result = execute_attempt(plan, task["row"])
        attempts.append(result)
        if result["status"] == "completed":
            break
    return {**task, "attempts": attempts, "final": attempts[-1]}


def _latest_outcome(db, iid):
    return db.execute(
        """SELECT * FROM ligand_arpeggio_runs WHERE ligand_instance_id=?
           ORDER BY run_id DESC LIMIT 1""",
        (iid,),
    ).fetchone()


def _legacy_failure_class(previous):
    if not previous:
        return None
    if previous["status"] == "timed_out":
        return "timeout"
    path = Path(previous["stderr_path"]) if previous["stderr_path"] else None
    text = path.read_text(errors="replace") if path and path.exists() else ""
    return classify_failure(text, previous["status"])


def _needs_canonical_derived(db, row):
    return (
        str(row["deposited_model_num"]) != "1"
        or bool(row["insertion_code_normalized"])
        or bool(
            db.execute(
                """SELECT 1 FROM ligand_instance_atoms
                   WHERE ligand_instance_id=? AND selected_conformer=1 AND altloc<>'' LIMIT 1""",
                (row["ligand_instance_id"],),
            ).fetchone()
        )
    )


def _make_plan(exe, row, run_id, base, strategy, timeout, is_fallback, fallback_reason, db,
               context_radius=DEFAULT_CONTEXT_RADIUS, trigger=None):
    attempt_dir = base / strategy
    attempt_dir.mkdir(parents=True, exist_ok=True)
    operations = {}
    if strategy == "direct_original":
        input_path = Path(row["source_cif_path"])
        provenance_status = "authoritative_source_input"
    else:
        input_path, atom_map, operations, manifest = build_derived_input(
            db, row, attempt_dir / "derived_input.cif", strategy, context_radius
        )
        persist_derived_map(db, run_id, strategy, atom_map)
        provenance_status = manifest["provenance_validation_status"]
    selector = selector_for(row)
    output = attempt_dir / f"{input_path.stem}.json"
    command = [exe, "-s", selector, "-o", str(attempt_dir), str(input_path)]
    return {
        "strategy": strategy,
        "input": str(input_path),
        "input_sha256": c.sha256(input_path),
        "selector": selector,
        "command": command,
        "command_json": json.dumps(command),
        "timeout": float(timeout),
        "is_fallback": bool(is_fallback),
        "fallback_reason": fallback_reason,
        "operations": operations,
        "stdout": str(attempt_dir / "stdout.log"),
        "stderr": str(attempt_dir / "stderr.log"),
        "output": str(output),
        "provenance_status": provenance_status,
        "trigger": trigger,
    }


def _persist_attempts(db, run_id, iid, attempts):
    for number, attempt in enumerate(attempts, 1):
        output_sha = c.sha256(attempt["raw"]) if attempt.get("raw") else None
        db.execute(
            """INSERT INTO arpeggio_attempts(
                 run_id,ligand_instance_id,attempt_number,input_strategy,input_path,input_sha256,
                 is_fallback,fallback_reason,sanitization_operations_json,timeout_seconds,start_time,
                 end_time,runtime_seconds,command,return_code,status,failure_class,error_class,error_message,
                 stdout_path,stderr_path,output_path,output_sha256,output_validation_status,
                 provenance_validation_status)
               VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)""",
            (
                run_id, iid, number, attempt["strategy"], attempt["input"], attempt["input_sha256"],
                int(attempt["is_fallback"]), attempt["fallback_reason"], json.dumps(attempt["operations"], sort_keys=True),
                attempt["timeout"], attempt["start_time"], attempt["end_time"], attempt["elapsed"],
                attempt["command_json"], attempt["returncode"], attempt["status"], attempt["failure_class"],
                attempt["error_class"], attempt["error_message"], attempt["stdout"], attempt["stderr"],
                str(attempt["raw"]) if attempt.get("raw") else attempt["output"], output_sha,
                attempt["output_validation_status"], attempt["provenance_validation_status"],
            ),
        )


def run(database, limit=None, pdb_id=None, instance_id=None, resume=False, per_instance_timeout=300,
        progress_every=25, workers=1, retry_timeout=DEFAULT_RETRY_TIMEOUT,
        fallback_radius=DEFAULT_CONTEXT_RADIUS):
    """Run bounded direct/fallback attempts, committing one final record per occurrence."""
    c.create_schema(database)
    exe, arpeggio_version = executable_provenance()
    if not exe:
        raise RuntimeError("pdbe-arpeggio executable not found")
    results = []
    with c.dbconn(database) as db:
        query = """SELECT i.ligand_instance_id FROM ligand_instances i
                   JOIN structures s ON s.structure_id=i.structure_id WHERE 1=1"""
        args = []
        if instance_id is None:
            query += " AND i.curation_status='included'"
        if pdb_id:
            query += " AND s.entry_id=?"
            args.append(pdb_id)
        if instance_id:
            query += " AND i.ligand_instance_id=?"
            args.append(instance_id)
        if resume:
            query += """ AND NOT EXISTS (
                         SELECT 1 FROM ligand_arpeggio_runs r
                         WHERE r.ligand_instance_id=i.ligand_instance_id AND r.status='completed')"""
        ids = [record[0] for record in db.execute(query + " ORDER BY i.ligand_instance_id", args)]
        if limit:
            ids = ids[:limit]
        db.execute("UPDATE ligand_arpeggio_runs SET status='interrupted' WHERE status='running'")

    tasks = []
    overall_started = time.monotonic()
    for preparation_sequence, iid in enumerate(ids, 1):
        with c.dbconn(database) as db:
            row = occurrence(db, iid)
            previous = _latest_outcome(db, iid)
            previous_status = previous["status"] if previous else None
            prior_failure = _legacy_failure_class(previous)
            run_id = c.run_start(
                db,
                "arpeggio",
                {
                    "ligand_instance_id": iid,
                    "method": VERSION,
                    "per_instance_timeout": per_instance_timeout,
                    "retry_timeout": retry_timeout,
                    "fallback_radius": fallback_radius,
                    "resume": resume,
                },
            )
            base = c.ROOT / "outputs" / "arpeggio" / "raw" / str(run_id)
            base.mkdir(parents=True, exist_ok=True)
            plans = []
            try:
                if resume and previous_status == "timed_out":
                    plans.append(_make_plan(exe, row, run_id, base, "sanitized_pocket", retry_timeout, True,
                                            prior_failure, db, fallback_radius))
                elif resume and previous_status in {"failed", "blocked", "interrupted"}:
                    plans.append(_make_plan(exe, row, run_id, base, "sanitized_full", per_instance_timeout, True,
                                            prior_failure, db, fallback_radius))
                else:
                    selector, blocked = selector_guard(db, row)
                    if blocked:
                        plans.append(_make_plan(exe, row, run_id, base, "sanitized_full", per_instance_timeout, True,
                                                blocked, db, fallback_radius))
                    else:
                        first_strategy = "canonical_derived" if _needs_canonical_derived(db, row) else "direct_original"
                        plans.append(_make_plan(exe, row, run_id, base, first_strategy, per_instance_timeout, False,
                                                None, db, fallback_radius))
            except Exception as exc:
                error_class = type(exc).__name__
                message = str(exc)
                db.execute(
                    """INSERT INTO ligand_arpeggio_runs(
                         run_id,ligand_instance_id,source_cif_sha256,selector,selector_namespace,
                         arpeggio_version,status,canonical_deposited_model_num,original_attempt_status,
                         fallback_attempted,fallback_reason,attempt_count,error_class,error_message,
                         provenance_validation_status,output_validation_status)
                       VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)""",
                    (run_id, iid, row["source_cif_sha256"], selector_for(row), "auth_asym_id/auth_seq_id",
                     arpeggio_version, "failed", row["deposited_model_num"], previous_status, 1, prior_failure,
                     0, error_class, message, "failed_during_derivation", "not_run"),
                )
                c.fail(db, run_id, "arpeggio", message, instance_id=iid, code="derived_input_failure")
                c.run_end(db, run_id, "failed", 1, 0, 0, 1)
                results.append({"iid": iid, "status": "failed", "reason": "derived_input_failure", "message": message})
                continue

            first = plans[0]
            db.execute(
                """INSERT INTO ligand_arpeggio_runs(
                     run_id,ligand_instance_id,source_cif_sha256,selector,selector_namespace,command,
                     arpeggio_version,stdout_path,stderr_path,status,derived_input_path,
                     derived_input_sha256,canonical_deposited_model_num,input_strategy,
                     original_attempt_status,fallback_attempted,fallback_reason,sanitization_operations_json,
                     attempt_count,timeout_seconds,provenance_validation_status,output_validation_status)
                   VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)""",
                (
                    run_id, iid, row["source_cif_sha256"], first["selector"], "auth_asym_id/auth_seq_id",
                    first["command_json"], arpeggio_version, first["stdout"], first["stderr"], "running",
                    first["input"] if first["strategy"] != "direct_original" else None,
                    first["input_sha256"] if first["strategy"] != "direct_original" else None,
                    row["deposited_model_num"], first["strategy"], previous_status,
                    int(any(plan["is_fallback"] for plan in plans)), prior_failure,
                    json.dumps(first["operations"], sort_keys=True), 0, first["timeout"],
                    first["provenance_status"], "pending",
                ),
            )
            tasks.append({"iid": iid, "rid": run_id, "row": dict(row), "plans": plans})
            if preparation_sequence % progress_every == 0 or preparation_sequence == len(ids):
                print(
                    f"arpeggio preparation: {preparation_sequence}/{len(ids)} "
                    f"last_ligand_instance_id={iid} first_mode={plans[0]['strategy']}",
                    flush=True,
                )

    total = len(tasks)
    with ThreadPoolExecutor(max_workers=max(1, workers)) as pool:
        futures = [pool.submit(execute_task, task) for task in tasks]
        for sequence, future in enumerate(as_completed(futures), 1):
            task = future.result()
            iid, run_id, row, attempts, final = task["iid"], task["rid"], task["row"], task["attempts"], task["final"]
            status = final["status"]
            result = {"iid": iid, "status": status, "completion_mode": final["strategy"], "attempts": len(attempts)}
            with c.dbconn(database) as db:
                _persist_attempts(db, run_id, iid, attempts)
                if status == "completed":
                    try:
                        raw_count, pairs, endpoints = parse_contacts(
                            db, run_id, iid, final["raw"], row, final["strategy"]
                        )
                        result.update(
                            raw_labels=raw_count,
                            unique_pairs=pairs,
                            selected_ligand_endpoints_observed=endpoints.get("observed", 0),
                            endpoints_reconciled=endpoints.get("reconciled", 0),
                            endpoints_unmatched=endpoints.get("unmatched", 0),
                            endpoints_ambiguous=endpoints.get("ambiguous", 0),
                        )
                    except Exception as exc:
                        status = "failed"
                        final["failure_class"] = "output_provenance_validation_failed"
                        final["error_class"] = type(exc).__name__
                        final["error_message"] = str(exc)
                        final["output_validation_status"] = str(exc)
                        db.execute("DELETE FROM arpeggio_raw_contact_labels WHERE run_id=?", (run_id,))
                        db.execute("DELETE FROM arpeggio_unique_atom_pairs WHERE run_id=?", (run_id,))
                        db.execute(
                            """UPDATE arpeggio_attempts SET status='failed',failure_class=?,error_class=?,
                               error_message=?,output_validation_status=? WHERE run_id=? AND attempt_number=?""",
                            (final["failure_class"], final["error_class"], final["error_message"],
                             final["output_validation_status"], run_id, len(attempts)),
                        )
                        result.update(status=status, reason=final["failure_class"], message=str(exc))
                output_sha = c.sha256(final["raw"]) if status == "completed" and final.get("raw") else None
                fallback_attempts = [attempt for attempt in attempts if attempt["is_fallback"]]
                db.execute(
                    """UPDATE ligand_arpeggio_runs SET command=?,stdout_path=?,stderr_path=?,exit_status=?,
                       output_sha256=?,status=?,end_time=?,derived_input_path=?,derived_input_sha256=?,
                       input_strategy=?,fallback_attempted=?,fallback_reason=?,sanitization_operations_json=?,
                       attempt_count=?,timeout_seconds=?,runtime_seconds=?,error_class=?,error_message=?,
                       provenance_validation_status=?,output_validation_status=?,completion_mode=? WHERE run_id=?""",
                    (
                        final["command_json"], final["stdout"], final["stderr"], final["returncode"], output_sha,
                        status, c.now(), final["input"] if final["strategy"] != "direct_original" else None,
                        final["input_sha256"] if final["strategy"] != "direct_original" else None,
                        final["strategy"], int(bool(fallback_attempts)), final["fallback_reason"],
                        json.dumps(final["operations"], sort_keys=True), len(attempts), final["timeout"],
                        sum(attempt["elapsed"] for attempt in attempts), final["error_class"], final["error_message"],
                        final["provenance_validation_status"], final["output_validation_status"],
                        final["strategy"] if status == "completed" else None, run_id,
                    ),
                )
                if status == "completed":
                    c.run_end(db, run_id, "completed", 1, 1, 0, 0)
                else:
                    reason = final["failure_class"] or "arpeggio_failed"
                    c.fail(db, run_id, "arpeggio", final["error_message"] or reason, instance_id=iid, code=reason)
                    c.run_end(db, run_id, status, 1, 0, 0, 1)
                    result.setdefault("reason", reason)
                    result.setdefault("elapsed_seconds", round(sum(a["elapsed"] for a in attempts), 3))
            results.append(result)
            if sequence % progress_every == 0 or sequence == total:
                print(
                    f"arpeggio progress: {sequence}/{total} last_ligand_instance_id={iid} "
                    f"status={status} mode={final['strategy']} attempts={len(attempts)} "
                    f"elapsed_seconds={time.monotonic()-overall_started:.1f} "
                    f"last_instance_seconds={sum(a['elapsed'] for a in attempts):.1f}",
                    flush=True,
                )
    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--database", default=str(c.ROOT / "viral_data_cif_v2.db"))
    parser.add_argument("--limit", type=int)
    parser.add_argument("--pdb-id")
    parser.add_argument("--ligand-instance-id", type=int)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--per-instance-timeout", type=float, default=300)
    parser.add_argument("--retry-timeout", type=float, default=DEFAULT_RETRY_TIMEOUT)
    parser.add_argument("--fallback-radius", type=float, default=DEFAULT_CONTEXT_RADIUS)
    parser.add_argument("--progress-every", type=int, default=25)
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args()
    print(
        run(
            args.database, args.limit, args.pdb_id, args.ligand_instance_id, args.resume,
            args.per_instance_timeout, args.progress_every, args.workers,
            args.retry_timeout, args.fallback_radius,
        )
    )
