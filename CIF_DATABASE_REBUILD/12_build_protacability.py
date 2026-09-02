#!/usr/bin/env python3
"""Build occurrence-resolved PROTACability triage tables.

This stage ports V-LiSEMOD PROTACability triage onto the CIF-native evidence
produced by stages 07-11. Ligand linkability and target-side lysine terminal-amine
surface accessibility are kept separate. Ligand-to-lysine distance is deliberately not
used as an ubiquitination or degrader-readiness metric. Atom-level exit-vector
descriptors from stage 10 are retained as ligand-side structural cues. Scores
are transparent structural-priority heuristics for follow-up; they are not
predictions of retained affinity, ubiquitination, or degradation.

Target-side scoring is restricted to protein chains that directly contact the
selected ligand in the completed Stage-09 Arpeggio result. Once a target chain
is selected, lysine accessibility is evaluated across the entire chain surface;
ligand-to-lysine distance is not used.
"""
from __future__ import annotations

import argparse
import importlib
import json
import math
import re
from collections import defaultdict

from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors
from Bio.PDB.Polypeptide import is_aa

c = importlib.import_module("00_common")
VERSION = "protacability-cif-v2.8"
GEOMETRY_VERSION = "cif-ligand-geometry-v2.4"
FUNCTIONAL_GROUP_VERSION = "rdkit-smarts-functional-groups-v2.3"
MAPPING_VERSION = "legacy_mcs_etkdg_uff_cif_v2.5"
HYDROGENS = {"H", "D", "T"}

# Preserved from the prior warhead-linkability implementation, normalized below.
STRONG_CONTACTS = {
    "hbond", "hydrogenbond", "ionic", "metal", "metalcomplex", "xbond",
    "halogenbond", "covalent", "aromatic", "cationpi", "donorpi", "halogenpi",
}
NON_INFORMATIVE_CONTACTS = {"proximal"}


def ensure_schema(database: str) -> None:
    c.create_schema(database)
    with c.dbconn(database) as db:
        db.executescript(
            """
            CREATE TABLE IF NOT EXISTS protacability_warhead_linkability (
                warhead_id INTEGER PRIMARY KEY,
                run_id INTEGER NOT NULL REFERENCES analysis_runs(run_id),
                ligand_instance_id INTEGER NOT NULL REFERENCES ligand_instances(ligand_instance_id),
                virus_name TEXT, protein_type TEXT,
                pdb_code TEXT NOT NULL, model_id TEXT NOT NULL,
                ligand_resname TEXT NOT NULL, ligand_chain TEXT NOT NULL,
                ligand_residue_id TEXT NOT NULL, ligand_insertion_code TEXT,
                ligand_context_class TEXT NOT NULL,
                source_inventory_row_count INTEGER NOT NULL,
                smiles_available INTEGER NOT NULL,
                representative_smiles TEXT, smiles_source TEXT,
                rdkit_available INTEGER NOT NULL, rdkit_valid_smiles INTEGER NOT NULL,
                mol_weight REAL, tpsa REAL, hbd INTEGER, hba INTEGER, rotatable_bonds INTEGER,
                heavy_atom_count_from_smiles INTEGER,
                pdb_ligand_heavy_atom_count INTEGER,
                pdb_to_smiles_mapped_atom_count INTEGER NOT NULL,
                functional_group_count INTEGER NOT NULL,
                functional_group_types TEXT NOT NULL,
                solvent_exposed_ligand_atom_count INTEGER NOT NULL,
                solvent_exposed_mapped_atom_count INTEGER NOT NULL,
                meaningful_contact_count INTEGER NOT NULL,
                strong_contact_count INTEGER NOT NULL,
                contact_atom_count INTEGER NOT NULL,
                strong_contact_atom_count INTEGER NOT NULL,
                candidate_linker_atom_count INTEGER NOT NULL,
                candidate_linker_atom_ids TEXT NOT NULL,
                outward_supported_candidate_count INTEGER NOT NULL,
                clear_exit_candidate_count INTEGER NOT NULL DEFAULT 0,
                ligand_exit_geometry_score REAL NOT NULL DEFAULT 0,
                direct_handle_candidate_count INTEGER NOT NULL,
                interaction_preservation_score REAL NOT NULL,
                warhead_linkability_score REAL NOT NULL,
                warhead_linkability_tier TEXT NOT NULL,
                warhead_linkability_label TEXT NOT NULL,
                warhead_flags TEXT NOT NULL,
                warhead_notes TEXT NOT NULL,
                method_version TEXT NOT NULL,
                UNIQUE(ligand_instance_id,method_version)
            );
            CREATE TABLE IF NOT EXISTS protacability_assessment (
                assessment_id INTEGER PRIMARY KEY,
                run_id INTEGER NOT NULL REFERENCES analysis_runs(run_id),
                ligand_instance_id INTEGER NOT NULL REFERENCES ligand_instances(ligand_instance_id),
                virus_name TEXT, protein_type TEXT,
                pdb_code TEXT NOT NULL, chain_id TEXT NOT NULL, model_id TEXT NOT NULL,
                chain_length_aa INTEGER NOT NULL,
                candidate_ligand_count INTEGER NOT NULL,
                candidate_ligand_resnames TEXT NOT NULL,
                lys_count INTEGER NOT NULL,
                exposed_lys_count INTEGER NOT NULL,
                exposed_lys_fraction REAL NOT NULL,
                near_ligand_lys_count INTEGER NOT NULL,
                near_ligand_exposed_lys_count INTEGER NOT NULL,
                min_lys_ligand_distance_a REAL,
                median_lys_ligand_distance_a REAL,
                total_sasa_a2 REAL NOT NULL,
                lysine_sasa_a2 REAL NOT NULL,
                lysine_surface_fraction REAL,
                isoelectric_point REAL,
                basic_fraction REAL NOT NULL,
                acidic_fraction REAL NOT NULL,
                polar_fraction REAL NOT NULL,
                hydrophobic_fraction REAL NOT NULL,
                has_candidate_ligand INTEGER NOT NULL,
                has_exposed_lysine INTEGER NOT NULL,
                has_ligand_proximal_exposed_lysine INTEGER NOT NULL,
                lysine_sidechain_sasa_a2 REAL NOT NULL DEFAULT 0,
                lysine_nz_sasa_a2 REAL NOT NULL DEFAULT 0,
                nz_observed_lys_count INTEGER NOT NULL DEFAULT 0,
                nz_observed_lys_fraction REAL NOT NULL DEFAULT 0,
                nz_exposed_lys_count_gt_1 INTEGER NOT NULL DEFAULT 0,
                nz_exposed_lys_fraction_gt_1 REAL NOT NULL DEFAULT 0,
                nz_exposed_lys_count_gt_5 INTEGER NOT NULL DEFAULT 0,
                nz_exposed_lys_fraction_gt_5 REAL NOT NULL DEFAULT 0,
                linker_docking_site_annotation TEXT NOT NULL,
                protein_ligand_druggability_proxy_score REAL NOT NULL,
                protacability_proxy_score REAL NOT NULL,
                protacability_tier TEXT NOT NULL,
                notes TEXT NOT NULL,
                target_chain_selection_basis TEXT NOT NULL DEFAULT '',
                ligand_target_contact_pair_count INTEGER NOT NULL DEFAULT 0,
                method_version TEXT NOT NULL,
                UNIQUE(ligand_instance_id,chain_id,method_version)
            );
            CREATE TABLE IF NOT EXISTS protacability_target_context (
                context_id INTEGER PRIMARY KEY,
                run_id INTEGER NOT NULL REFERENCES analysis_runs(run_id),
                ligand_instance_id INTEGER NOT NULL REFERENCES ligand_instances(ligand_instance_id),
                arpeggio_run_id INTEGER REFERENCES analysis_runs(run_id),
                target_context_status TEXT NOT NULL,
                contacting_protein_chain_count INTEGER NOT NULL,
                contacting_protein_chain_ids TEXT NOT NULL,
                target_chain_selection_basis TEXT NOT NULL,
                notes TEXT NOT NULL,
                method_version TEXT NOT NULL,
                UNIQUE(ligand_instance_id,method_version)
            );
            CREATE TABLE IF NOT EXISTS protacability_degrader_readiness (
                readiness_id INTEGER PRIMARY KEY,
                run_id INTEGER NOT NULL REFERENCES analysis_runs(run_id),
                ligand_instance_id INTEGER NOT NULL REFERENCES ligand_instances(ligand_instance_id),
                virus_name TEXT, protein_type TEXT,
                pdb_code TEXT NOT NULL, chain_id TEXT NOT NULL, model_id TEXT NOT NULL,
                best_ligand_resname TEXT NOT NULL,
                best_ligand_chain TEXT NOT NULL,
                best_ligand_residue_id TEXT NOT NULL,
                protein_structural_priority_score REAL NOT NULL,
                warhead_linkability_score REAL NOT NULL,
                target_lysine_accessibility_score REAL NOT NULL,
                ternary_geometry_cue_score REAL NOT NULL,
                ligand_exit_geometry_score REAL NOT NULL DEFAULT 0,
                clear_exit_candidate_count INTEGER NOT NULL DEFAULT 0,
                degrader_design_readiness_score REAL NOT NULL,
                degrader_design_readiness_tier TEXT NOT NULL,
                evidence_level TEXT NOT NULL,
                best_linker_geometry_class TEXT NOT NULL,
                short_linker_geometry_feasible INTEGER NOT NULL,
                medium_linker_geometry_feasible INTEGER NOT NULL,
                long_linker_geometry_feasible INTEGER NOT NULL,
                exposed_lys_count INTEGER NOT NULL,
                lys_count INTEGER NOT NULL,
                exposed_lys_fraction REAL NOT NULL,
                lysine_surface_fraction REAL,
                min_lys_ligand_distance_a REAL,
                near_ligand_exposed_lys_count INTEGER NOT NULL,
                nz_observed_lys_count INTEGER NOT NULL DEFAULT 0,
                nz_observed_lys_fraction REAL NOT NULL DEFAULT 0,
                nz_exposed_lys_count_gt_1 INTEGER NOT NULL DEFAULT 0,
                nz_exposed_lys_fraction_gt_1 REAL NOT NULL DEFAULT 0,
                nz_exposed_lys_count_gt_5 INTEGER NOT NULL DEFAULT 0,
                nz_exposed_lys_fraction_gt_5 REAL NOT NULL DEFAULT 0,
                candidate_ligand_resnames TEXT NOT NULL,
                readiness_flags TEXT NOT NULL,
                readiness_notes TEXT NOT NULL,
                best_chain_for_instance INTEGER NOT NULL DEFAULT 0,
                target_chain_selection_basis TEXT NOT NULL DEFAULT '',
                ligand_target_contact_pair_count INTEGER NOT NULL DEFAULT 0,
                method_version TEXT NOT NULL,
                UNIQUE(ligand_instance_id,chain_id,method_version)
            );
            """
        )
        migrations = {
            "protacability_warhead_linkability": (("clear_exit_candidate_count", "INTEGER NOT NULL DEFAULT 0"), ("ligand_exit_geometry_score", "REAL NOT NULL DEFAULT 0")),
            "protacability_assessment": (
                ("lysine_sidechain_sasa_a2", "REAL NOT NULL DEFAULT 0"),
                ("lysine_nz_sasa_a2", "REAL NOT NULL DEFAULT 0"),
                ("nz_observed_lys_count", "INTEGER NOT NULL DEFAULT 0"),
                ("nz_observed_lys_fraction", "REAL NOT NULL DEFAULT 0"),
                ("nz_exposed_lys_count_gt_1", "INTEGER NOT NULL DEFAULT 0"),
                ("nz_exposed_lys_fraction_gt_1", "REAL NOT NULL DEFAULT 0"),
                ("nz_exposed_lys_count_gt_5", "INTEGER NOT NULL DEFAULT 0"),
                ("nz_exposed_lys_fraction_gt_5", "REAL NOT NULL DEFAULT 0"),
                ("target_chain_selection_basis", "TEXT NOT NULL DEFAULT ''"),
                ("ligand_target_contact_pair_count", "INTEGER NOT NULL DEFAULT 0"),
            ),
            "protacability_degrader_readiness": (
                ("ligand_exit_geometry_score", "REAL NOT NULL DEFAULT 0"), ("clear_exit_candidate_count", "INTEGER NOT NULL DEFAULT 0"),
                ("nz_observed_lys_count", "INTEGER NOT NULL DEFAULT 0"),
                ("nz_observed_lys_fraction", "REAL NOT NULL DEFAULT 0"),
                ("nz_exposed_lys_count_gt_1", "INTEGER NOT NULL DEFAULT 0"),
                ("nz_exposed_lys_fraction_gt_1", "REAL NOT NULL DEFAULT 0"),
                ("nz_exposed_lys_count_gt_5", "INTEGER NOT NULL DEFAULT 0"),
                ("nz_exposed_lys_fraction_gt_5", "REAL NOT NULL DEFAULT 0"),
                ("target_chain_selection_basis", "TEXT NOT NULL DEFAULT ''"),
                ("ligand_target_contact_pair_count", "INTEGER NOT NULL DEFAULT 0"),
            ),
        }
        for table, columns in migrations.items():
            existing = {r[1] for r in db.execute(f"PRAGMA table_info({table})")}
            for column, kind in columns:
                if column not in existing:
                    db.execute(f"ALTER TABLE {table} ADD COLUMN {column} {kind}")


def clamp(value, lo=0.0, hi=100.0):
    return max(lo, min(hi, float(value)))


def normalize_contact(label: str) -> str:
    token = re.sub(r"[^a-z0-9]+", "", str(label or "").lower())
    aliases = {
        "hbond": "hbond", "hydrogenbond": "hydrogenbond", "weakhbond": "weakhbond",
        "metalcomplex": "metalcomplex", "metal": "metal", "halogenbond": "halogenbond",
        "xbond": "xbond", "cationpi": "cationpi", "donorpi": "donorpi", "halogenpi": "halogenpi",
        "carbonpi": "carbonpi", "metsulphurpi": "metsulphurpi", "vdwclash": "vdwclash",
    }
    return aliases.get(token, token)


def is_strong_contact(label: str) -> bool:
    token = normalize_contact(label)
    if token in STRONG_CONTACTS:
        return True
    return any(x in token for x in ("covalent", "ionic", "aromatic", "cationpi", "donorpi", "halogenpi"))


def is_meaningful_contact(label: str) -> bool:
    return normalize_contact(label) not in NON_INFORMATIVE_CONTACTS


def _latest_mapping_run(db, iid):
    """Use only the repaired Stage-07 v2.5 mapping generation."""
    r = db.execute(
        """SELECT MAX(run_id) FROM ligand_mapping_runs WHERE ligand_instance_id=?
           AND method_version=? AND downstream_mapping_eligibility=1
           AND mapping_status IN ('complete','complete_altloc_resolved','partial_ccd_difference')""",
        (iid, MAPPING_VERSION),
    ).fetchone()
    return r[0] if r and r[0] is not None else None


def _latest_sasa_run(db, iid):
    r = db.execute("SELECT MAX(run_id) FROM ligand_sasa_atoms WHERE ligand_instance_id=? AND status='complete'", (iid,)).fetchone()
    return r[0] if r and r[0] is not None else None


def _latest_arpeggio_run(db, iid):
    r = db.execute("SELECT MAX(run_id) FROM ligand_arpeggio_runs WHERE ligand_instance_id=? AND status='completed'", (iid,)).fetchone()
    return r[0] if r and r[0] is not None else None


def _protein_component(component_id: str) -> bool:
    comp = c.norm(component_id).upper()
    return bool(comp and is_aa(comp, standard=False))


def contacting_target_chains(db, iid: int):
    """Return direct ligand-contacting protein chains from the latest Stage-09 run.

    Arpeggio is used only to establish which protein chain(s) constitute the
    ligand-binding target context.  Lysine scoring is then performed across the
    full surface of those chains.  Ligand-to-lysine distance is never used.

    For fallback/sanitized Arpeggio inputs, partner_source_atom_site_id is used
    to restore the deposited source-chain identity when provenance is available.
    Direct-input contacts fall back to the endpoint auth_asym_id/label_asym_id.
    Contact-pair counts are deduplicated across Arpeggio interaction labels.
    """
    arp_run = _latest_arpeggio_run(db, iid)
    if arp_run is None:
        raise ValueError("no completed Arpeggio run; run stage 09 first")

    available = {
        str(r["chain_id"]): r
        for r in db.execute(
            """SELECT * FROM target_chain_geometry
               WHERE ligand_instance_id=? AND method_version=?""",
            (iid, GEOMETRY_VERSION),
        )
    }
    if not available:
        return arp_run, {}, "not_applicable_no_protein_atoms"

    # Restore source identities for derived/sanitized Arpeggio inputs.
    provenance = {}
    for r in db.execute(
        """SELECT source_atom_site_id,source_auth_asym_id,source_label_asym_id,
                  source_component_id
           FROM arpeggio_derived_atom_map WHERE run_id=?""",
        (arp_run,),
    ):
        provenance[str(r["source_atom_site_id"])] = {
            "auth_asym_id": c.norm(r["source_auth_asym_id"]),
            "label_asym_id": c.norm(r["source_label_asym_id"]),
            "component_id": c.norm(r["source_component_id"]).upper(),
        }

    per_chain_pairs = defaultdict(set)
    rows = db.execute(
        """SELECT ligand_instance_atom_id,partner_identity_json,
                  partner_source_atom_site_id,partner_mapping_status
           FROM arpeggio_raw_contact_labels
           WHERE ligand_instance_id=? AND run_id=?
             AND filter_class='raw_environment'
             AND ligand_instance_atom_id IS NOT NULL""",
        (iid, arp_run),
    ).fetchall()

    for r in rows:
        try:
            partner = json.loads(r["partner_identity_json"] or "{}")
        except Exception:
            partner = {}

        source_atom_site_id = c.norm(r["partner_source_atom_site_id"])
        source = provenance.get(source_atom_site_id) if source_atom_site_id else None
        if source:
            chain_id = source["auth_asym_id"] or source["label_asym_id"]
            component_id = source["component_id"]
            partner_key = f"source_atom_site:{source_atom_site_id}"
        else:
            chain_id = c.norm(partner.get("auth_asym_id")) or c.norm(partner.get("label_asym_id"))
            component_id = c.norm(partner.get("label_comp_id") or partner.get("auth_comp_id")).upper()
            # partner_identity_json includes residue/atom identity, so this is
            # stable enough to deduplicate the same pair across interaction labels.
            partner_key = r["partner_identity_json"] or "{}"

        if not chain_id or chain_id not in available:
            continue
        if not _protein_component(component_id):
            continue

        per_chain_pairs[chain_id].add((int(r["ligand_instance_atom_id"]), partner_key))

    evidence = {
        chain_id: {"chain": available[chain_id], "contact_pair_count": len(pairs)}
        for chain_id, pairs in per_chain_pairs.items()
        if pairs
    }
    if not evidence:
        return arp_run, {}, "not_applicable_no_contacting_protein_chain"
    return arp_run, evidence, "applicable_contacting_protein_chain"


def _latest_geometry_run(db, iid):
    r = db.execute(
        "SELECT MAX(run_id) FROM ligand_atom_geometry WHERE ligand_instance_id=? AND method_version=?",
        (iid, GEOMETRY_VERSION),
    ).fetchone()
    return r[0] if r and r[0] is not None else None


def _latest_functional_group_run(db, iid):
    """Use only the corrected Stage-11 source-SMILES atom namespace."""
    r = db.execute(
        """SELECT run_id FROM ligand_functional_group_summary
           WHERE ligand_instance_id=? AND method_version=? AND status='complete'
           ORDER BY run_id DESC LIMIT 1""",
        (iid, FUNCTIONAL_GROUP_VERSION),
    ).fetchone()
    return r[0] if r and r[0] is not None else None


def load_atom_evidence(db, iid: int):
    atoms = {
        r["ligand_instance_atom_id"]: {
            "ligand_instance_atom_id": r["ligand_instance_atom_id"],
            "atom_site_id": r["atom_site_id"],
            "exact_atom": r["auth_atom_id"] or r["label_atom_id"],
            "element": c.norm(r["element"]).upper(),
            "mapped": False, "smiles_atom_indices": set(),
            "sasa_area": None, "exposed": False,
            "labels": [], "meaningful_contact_count": 0, "strong_contact_count": 0,
            "unique_partner_count": 0,
            "functional_groups": set(), "tractability_roles": set(),
            "outward_score": None, "points_away": None, "nearest_protein_distance_a": None,
            "exit_vector_clear": None, "local_corridor_clear": None, "forward_clearance_a": None,
            "forward_obstruction_count": None, "forward_clearance_reaches_cap": None,
        }
        for r in db.execute(
            """SELECT ligand_instance_atom_id,atom_site_id,label_atom_id,auth_atom_id,element
               FROM ligand_instance_atoms WHERE ligand_instance_id=? AND selected_conformer=1""", (iid,)
        )
    }
    mapping_run = _latest_mapping_run(db, iid)
    if mapping_run is not None:
        for r in db.execute(
            """SELECT ligand_instance_atom_id,smiles_atom_index FROM ligand_smiles_atom_mapping
               WHERE ligand_instance_id=? AND run_id=? AND ligand_instance_atom_id IS NOT NULL""", (iid, mapping_run)
        ):
            if r["ligand_instance_atom_id"] in atoms:
                atoms[r["ligand_instance_atom_id"]]["mapped"] = True
                if r["smiles_atom_index"] is not None:
                    atoms[r["ligand_instance_atom_id"]]["smiles_atom_indices"].add(int(r["smiles_atom_index"]))
    sasa_run = _latest_sasa_run(db, iid)
    if sasa_run is not None:
        for r in db.execute(
            """SELECT ligand_instance_atom_id,sasa_area,legacy_exposed FROM ligand_sasa_atoms
               WHERE ligand_instance_id=? AND run_id=?""", (iid, sasa_run)
        ):
            if r["ligand_instance_atom_id"] in atoms:
                atoms[r["ligand_instance_atom_id"]]["sasa_area"] = float(r["sasa_area"])
                atoms[r["ligand_instance_atom_id"]]["exposed"] = bool(r["legacy_exposed"])
    arp_run = _latest_arpeggio_run(db, iid)
    if arp_run is not None:
        partner_counts = defaultdict(set)
        for r in db.execute(
            """SELECT ligand_instance_atom_id,interaction_label,partner_identity_json
               FROM arpeggio_raw_contact_labels
               WHERE ligand_instance_id=? AND run_id=? AND filter_class='raw_environment'
                 AND ligand_instance_atom_id IS NOT NULL""", (iid, arp_run)
        ):
            aid = r["ligand_instance_atom_id"]
            if aid not in atoms:
                continue
            label = r["interaction_label"]
            atoms[aid]["labels"].append(label)
            if is_meaningful_contact(label):
                atoms[aid]["meaningful_contact_count"] += 1
            if is_strong_contact(label):
                atoms[aid]["strong_contact_count"] += 1
            partner_counts[aid].add(r["partner_identity_json"])
        for aid, partners in partner_counts.items():
            atoms[aid]["unique_partner_count"] = len(partners)
    # Functional groups from the corrected Stage 11 only.  Older Stage-11
    # runs used canonical-SMILES atom ordering for these joins and must never
    # be mixed into atom-level evidence.
    if db.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND name='ligand_functional_group_atoms'").fetchone():
        fg_run = _latest_functional_group_run(db, iid)
        if fg_run is not None:
            for r in db.execute(
                """SELECT ligand_instance_atom_id,functional_group,tractability_role
                   FROM ligand_functional_group_atoms
                   WHERE ligand_instance_id=? AND run_id=? AND method_version=?
                     AND mapping_status='mapped_element_validated'
                     AND ligand_instance_atom_id IS NOT NULL""",
                (iid, fg_run, FUNCTIONAL_GROUP_VERSION),
            ):
                aid = r["ligand_instance_atom_id"]
                if aid in atoms:
                    atoms[aid]["functional_groups"].add(r["functional_group"])
                    atoms[aid]["tractability_roles"].add(r["tractability_role"])
    if db.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND name='ligand_atom_geometry'").fetchone():
        geometry_run = _latest_geometry_run(db, iid)
        if geometry_run is not None:
            for r in db.execute(
                """SELECT ligand_instance_atom_id,outward_score,points_away_from_pocket,nearest_protein_distance_a,
                          exit_vector_clear,local_corridor_clear,forward_clearance_a,forward_obstruction_count,
                          forward_clearance_reaches_cap
                   FROM ligand_atom_geometry WHERE ligand_instance_id=? AND run_id=? AND method_version=?""",
                (iid, geometry_run, GEOMETRY_VERSION)
            ):
                aid = r["ligand_instance_atom_id"]
                if aid in atoms:
                    atoms[aid]["outward_score"] = r["outward_score"]
                    atoms[aid]["points_away"] = r["points_away_from_pocket"]
                    atoms[aid]["nearest_protein_distance_a"] = r["nearest_protein_distance_a"]
                    atoms[aid]["exit_vector_clear"] = r["exit_vector_clear"]
                    atoms[aid]["local_corridor_clear"] = r["local_corridor_clear"]
                    atoms[aid]["forward_clearance_a"] = r["forward_clearance_a"]
                    atoms[aid]["forward_obstruction_count"] = r["forward_obstruction_count"]
                    atoms[aid]["forward_clearance_reaches_cap"] = r["forward_clearance_reaches_cap"]
    return atoms


def chemical_context(atom):
    roles = atom["tractability_roles"]
    if "direct_handle_context" in roles:
        return "direct_handle_context"
    if "conditional_handle_context" in roles:
        return "conditional_handle_context"
    if roles:
        return "structural_context_only"
    return "unclassified_atom_context"


def candidate_core(atom):
    return (
        atom["element"] not in HYDROGENS
        and atom["mapped"] and atom["exposed"] and atom["strong_contact_count"] == 0
    )


def attachment_atom_rows_for_instance(db, iid: int):
    rows = []
    for atom in load_atom_evidence(db, iid).values():
        core = candidate_core(atom)
        chem = chemical_context(atom)
        outward = atom["points_away"]
        chemical_support = chem in {"direct_handle_context", "conditional_handle_context"}
        clear_exit = atom["local_corridor_clear"]
        high_priority = bool(core and outward == 1 and clear_exit == 1 and chemical_support)
        score = 0.0
        if atom["mapped"]: score += 20
        if atom["exposed"]: score += 30
        if atom["strong_contact_count"] == 0: score += 20
        if atom["unique_partner_count"] == 0: score += 10
        elif atom["unique_partner_count"] <= 2: score += 7
        elif atom["unique_partner_count"] <= 5: score += 3
        if chem == "direct_handle_context": score += 10
        elif chem == "conditional_handle_context": score += 6
        elif chem == "structural_context_only": score += 3
        else: score += 1
        if outward == 1:
            score += 8 if (atom["outward_score"] or 0) >= 0.25 else 5
        elif outward is None:
            score += 1
        if outward == 1 and clear_exit == 1:
            score += 7
        if not core:
            score = min(score, 49.0)
        score = round(clamp(score), 2)
        tier = "High attachment-site priority" if score >= 80 else (
            "Moderate attachment-site priority" if score >= 60 else (
                "Exploratory attachment-site priority" if score >= 40 else "Low attachment-site priority"
            )
        )
        rows.append({**atom, "candidate_core": int(core), "chemical_context": chem,
                     "chemical_support": int(chemical_support), "high_priority": int(high_priority),
                     "attachment_priority_score": score, "attachment_priority_tier": tier})
    return rows


def rdkit_info(smiles):
    out = {"rdkit_available": 1, "rdkit_valid_smiles": 0, "mol_weight": None, "tpsa": None,
           "hbd": None, "hba": None, "rotatable_bonds": None, "heavy_atom_count_from_smiles": None}
    if not smiles:
        return out
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return out
    out.update(
        rdkit_valid_smiles=1,
        mol_weight=round(float(Descriptors.MolWt(mol)), 4),
        tpsa=round(float(rdMolDescriptors.CalcTPSA(mol)), 4),
        hbd=int(Lipinski.NumHDonors(mol)), hba=int(Lipinski.NumHAcceptors(mol)),
        rotatable_bonds=int(Lipinski.NumRotatableBonds(mol)),
        heavy_atom_count_from_smiles=int(mol.GetNumHeavyAtoms()),
    )
    return out


def interaction_preservation_score(atom_rows):
    exposed = [a for a in atom_rows if a["exposed"] and a["mapped"]]
    if not exposed:
        return 0.0
    candidates = sum(a["candidate_core"] for a in exposed)
    strong_atoms = sum(a["strong_contact_count"] > 0 for a in exposed)
    score = 100.0 * candidates / len(exposed)
    if strong_atoms >= len(exposed) and candidates == 0:
        score -= 30
    elif strong_atoms:
        score -= min(20, strong_atoms * 4)
    return round(clamp(score), 2)


def warhead_for_instance(db, iid: int):
    inst = db.execute(
        """SELECT i.*,l.smiles,l.canonical_smiles,l.smiles_source,l.chemical_status,s.entry_id
           FROM ligand_instances i JOIN ligands l ON l.ligand_id=i.ligand_id
           JOIN structures s ON s.structure_id=i.structure_id WHERE i.ligand_instance_id=?""", (iid,)
    ).fetchone()
    inv = db.execute("SELECT * FROM protacability_ligand_inventory WHERE ligand_instance_id=? AND method_version=? ORDER BY run_id DESC LIMIT 1", (iid, GEOMETRY_VERSION)).fetchone()
    if inst is None or inv is None:
        raise ValueError("missing stage-10 ligand inventory or ligand instance")
    atoms = attachment_atom_rows_for_instance(db, iid)
    smiles = inst["canonical_smiles"] or inst["smiles"] or ""
    rdk = rdkit_info(smiles)
    mapped = sum(a["mapped"] for a in atoms)
    exposed = sum(a["exposed"] for a in atoms)
    exposed_mapped = sum(a["exposed"] and a["mapped"] for a in atoms)
    meaningful = sum(a["meaningful_contact_count"] for a in atoms)
    strong = sum(a["strong_contact_count"] for a in atoms)
    contact_atoms = sum(a["meaningful_contact_count"] > 0 for a in atoms)
    strong_atoms = sum(a["strong_contact_count"] > 0 for a in atoms)
    candidates = [a for a in atoms if a["candidate_core"]]
    outward_candidates = sum(a["candidate_core"] and a["points_away"] == 1 for a in atoms)
    clear_exit_candidates = sum(a["candidate_core"] and a["points_away"] == 1 and a["local_corridor_clear"] == 1 for a in atoms)
    direct_candidates = sum(a["candidate_core"] and a["chemical_context"] == "direct_handle_context" for a in atoms)
    fg_summary = db.execute("SELECT * FROM ligand_functional_group_summary WHERE ligand_instance_id=? AND method_version=? AND status='complete' ORDER BY run_id DESC LIMIT 1", (iid, FUNCTIONAL_GROUP_VERSION)).fetchone()
    fg_count = int(fg_summary["functional_group_type_count"]) if fg_summary else 0
    fg_types = fg_summary["functional_group_types"] if fg_summary else ""

    score = 10.0  # included ligand candidate / small-molecule context, matching legacy logic
    flags, notes = [], ["included ligand-candidate context"]
    if inv["status"] == "not_applicable_no_protein_atoms":
        flags.append("no_protein_target_geometry")
        notes.append("no recognized protein atoms in the deposited model; target-side PROTACability is not applicable")
    if smiles:
        score += 12; notes.append("SMILES available")
    else:
        flags.append("missing_smiles"); score -= 8
    if rdk["rdkit_valid_smiles"]:
        score += 8; notes.append("valid RDKit molecule")
    elif smiles:
        flags.append("invalid_smiles_for_rdkit"); score -= 10
    if mapped:
        score += min(15, 5 + mapped * 0.5); notes.append("PDB-to-SMILES atom mapping available")
    else:
        flags.append("missing_pdb_to_smiles_mapping")
    if exposed_mapped:
        score += min(20, 8 + exposed_mapped * 3); notes.append("mapped solvent-exposed ligand atoms detected")
    elif exposed:
        score += min(12, 4 + exposed * 2); notes.append("solvent-exposed atoms detected but mapping incomplete")
    else:
        flags.append("no_solvent_exposed_ligand_atoms_detected")
    if fg_count:
        score += min(15, 5 + fg_count * 2); notes.append("functional-group annotations available")
    else:
        flags.append("no_functional_group_annotations")
    ncan = len(candidates)
    exit_geometry_score = round(clamp(
        (40.0 * outward_candidates / ncan + 60.0 * clear_exit_candidates / ncan) if ncan else 0.0
    ), 2)
    if ncan >= 5:
        score += 20; notes.append("multiple mapped exposed atoms without strong contacts")
    elif ncan >= 2:
        score += 16; notes.append("several mapped exposed atoms without strong contacts")
    elif ncan == 1:
        score += 10; notes.append("one mapped exposed atom without strong contacts")
    else:
        flags.append("no_candidate_linker_atoms_after_contact_filtering"); score -= 10
    preservation = interaction_preservation_score(atoms)
    score += preservation * 0.10
    if preservation < 25:
        flags.append("high_risk_of_disrupting_binding_contacts")
    score = round(clamp(score), 2)
    if score >= 75:
        tier, label = "High warhead-linkability evidence", "Linkerable warhead candidate"
    elif score >= 55:
        tier, label = "Moderate warhead-linkability evidence", "Plausible warhead candidate"
    elif score >= 35:
        tier, label = "Exploratory warhead-linkability evidence", "Needs manual linker-site review"
    else:
        tier, label = "Weak warhead-linkability evidence", "Poor or incomplete warhead evidence"
    if outward_candidates:
        notes.append(f"{outward_candidates} candidate atoms also point away from the pocket centroid")
    if clear_exit_candidates:
        notes.append(f"{clear_exit_candidates} candidate atoms have a locally clear forward exit-vector corridor")
    if direct_candidates:
        notes.append(f"{direct_candidates} candidate atoms occur in direct-handle functional-group context")
    notes.append("ligand-centered linkability; target lysines are evaluated separately")
    return {
        "iid": iid, "inv": inv, "inst": inst, "smiles": smiles, "rdk": rdk,
        "mapped": mapped, "fg_count": fg_count, "fg_types": fg_types,
        "exposed": exposed, "exposed_mapped": exposed_mapped, "meaningful": meaningful, "strong": strong,
        "contact_atoms": contact_atoms, "strong_atoms": strong_atoms, "candidates": candidates,
        "outward_candidates": outward_candidates, "clear_exit_candidates": clear_exit_candidates,
        "exit_geometry_score": exit_geometry_score, "direct_candidates": direct_candidates,
        "preservation": preservation, "score": score, "tier": tier, "label": label,
        "flags": ";".join(dict.fromkeys(flags)), "notes": "; ".join(notes),
    }


def target_lysine_accessibility(chain):
    """Score target-side *surface lysine availability* only.

    Ligand-to-lysine distance is deliberately excluded. The relevant binary-
    structure question is whether the target presents lysine terminal amines on
    the solvent-accessible protein surface that could in principle be reached by
    ubiquitination machinery after ternary-complex formation.

    The denominator for exposure is ALL lysines in the modeled chain, not only
    lysines whose NZ atom happens to be observed. This prevents sparse coordinate
    coverage from making one exposed observed lysine look like near-complete target
    accessibility. Missing NZ coordinates remain unknown, but they lower confidence
    and therefore the structural-priority score rather than being called buried.

    Transparent heuristic (0-100):
      10%  NZ coordinate coverage among all lysines
      60%  fraction of all lysines with NZ SASA > 1 A^2
      30%  fraction of all lysines with NZ SASA > 5 A^2

    The >1 A^2 term is the primary surface-accessibility criterion; >5 A^2 rewards
    clearly exposed terminal amines. No ligand-proximity term is present.
    """
    lys_count = int(chain["lys_count"] or 0)
    if lys_count <= 0:
        return 0.0

    nz_observed = int(chain["nz_observed_lys_count"] or 0)
    gt1_count = int(chain["nz_exposed_lys_count_gt_1"] or 0)
    gt5_count = int(chain["nz_exposed_lys_count_gt_5"] or 0)

    coverage = nz_observed / lys_count
    exposed_gt1_all_lys = gt1_count / lys_count
    exposed_gt5_all_lys = gt5_count / lys_count

    score = (
        10.0 * coverage
        + 60.0 * exposed_gt1_all_lys
        + 30.0 * exposed_gt5_all_lys
    )
    return round(clamp(score), 2)


def protein_priority(chain):
    score = target_lysine_accessibility(chain)
    lys_count = int(chain["lys_count"] or 0)
    nz_observed = int(chain["nz_observed_lys_count"] or 0)
    gt1 = int(chain["nz_exposed_lys_count_gt_1"] or 0)
    gt5 = int(chain["nz_exposed_lys_count_gt_5"] or 0)
    notes = [
        f"NZ observed for {nz_observed}/{lys_count} lysines" if lys_count else "no lysines detected",
        f"{gt1}/{lys_count} total lysines have observed NZ SASA >1 A^2" if lys_count else "",
        f"{gt5}/{lys_count} total lysines have observed NZ SASA >5 A^2" if lys_count else "",
        "target-side score measures surface-accessible lysine NZ availability across the chain; ligand-to-lysine distance is not scored",
        "missing NZ coordinates reduce structural confidence and are not classified as buried",
    ]
    tier = "High structural priority" if score >= 70 else ("Moderate structural priority" if score >= 45 else "Low structural priority")
    return score, tier, "; ".join(x for x in notes if x)


def ligand_exit_geometry(warhead):
    ncan = len(warhead["candidates"])
    outward = int(warhead["outward_candidates"])
    clear = int(warhead["clear_exit_candidates"])
    score = float(warhead["exit_geometry_score"])
    if clear > 0:
        label = "Open outward ligand exit-vector cue"
    elif outward > 0:
        label = "Outward ligand atom cue with local forward obstruction"
    elif ncan > 0:
        label = "Candidate attachment atoms without an outward exit-vector cue"
    else:
        label = "No candidate ligand attachment atom"
    return score, label


def readiness_tier(score):
    return "High degrader-design readiness" if score >= 75 else (
        "Moderate degrader-design readiness" if score >= 55 else (
            "Exploratory degrader-design readiness" if score >= 35 else "Weak degrader-design readiness"
        )
    )


def evidence_level(warhead, lys, exit_geometry):
    if warhead >= 75 and lys >= 55:
        return "Strong: candidate warhead plus target lysine surface accessibility"
    if warhead >= 55 and lys >= 45:
        return "Moderate: plausible warhead with accessible target lysine surface"
    if warhead >= 55 and lys < 45:
        return "Warhead-supported but target lysine surface accessibility is weak"
    if warhead < 55 and lys >= 55:
        return "Target-surface-supported but warhead evidence is weak"
    if exit_geometry >= 60:
        return "Exit-vector-supported exploratory case"
    return "Limited degrader-readiness evidence"


def _write_report(database):
    c.dirs()
    with c.dbconn(database) as db:
        wc = db.execute("SELECT count(*) FROM protacability_warhead_linkability WHERE method_version=?", (VERSION,)).fetchone()[0]
        ac = db.execute("SELECT count(*) FROM protacability_assessment WHERE method_version=?", (VERSION,)).fetchone()[0]
        rc = db.execute("SELECT count(*) FROM protacability_degrader_readiness WHERE method_version=?", (VERSION,)).fetchone()[0]
        assessed_instances = db.execute("SELECT count(DISTINCT ligand_instance_id) FROM protacability_assessment WHERE method_version=?", (VERSION,)).fetchone()[0]
        contexts = db.execute("SELECT target_context_status,count(*) n FROM protacability_target_context WHERE method_version=? GROUP BY 1 ORDER BY n DESC", (VERSION,)).fetchall()
        tiers = db.execute("SELECT degrader_design_readiness_tier,count(*) n FROM protacability_degrader_readiness WHERE method_version=? GROUP BY 1 ORDER BY n DESC", (VERSION,)).fetchall()
    lines = ["# Stage 12 PROTACability report", "", f"* Method: {VERSION}", f"* Warhead rows: {wc}", f"* Target-assessed ligand instances: {assessed_instances}", f"* Assessment rows: {ac}", f"* Readiness rows: {rc}", "", "## Target-context status"]
    lines += [f"* {r['target_context_status']}: {r['n']}" for r in contexts]
    lines += ["", "## Readiness tiers"]
    lines += [f"* {r['degrader_design_readiness_tier']}: {r['n']}" for r in tiers]
    lines += ["", "Scores are structural-priority heuristics for hypothesis generation and are not experimental degradation predictions.",
              "Ligand-to-lysine distance is not used in the target accessibility or degrader-readiness scores; target-side lysine evidence is the fraction of ALL chain lysines with observed solvent-exposed NZ terminal amines plus NZ coordinate coverage.",
              "Target-side PROTACability/readiness is calculated only for protein chains that directly contact the selected ligand in the completed Stage-09 Arpeggio result.",
              "Once a target chain is selected, lysine accessibility is evaluated across the entire selected chain surface; ligand-to-lysine distance is not used.",
              "Models with no recognized protein atoms, or with no ligand-contacting protein chain, retain ligand-side warhead evidence but are explicitly not applicable for target-side PROTACability/readiness."]
    (c.ROOT / "outputs" / "PROTACABILITY_STAGE_REPORT.md").write_text("\n".join(lines) + "\n")


def run(database: str, limit=None, pdb_id=None, instance_id=None, resume=False, progress_every=100):
    ensure_schema(database); RDLogger.DisableLog("rdApp.*"); c.dirs()
    required = {"protacability_ligand_inventory", "target_chain_geometry", "ligand_functional_group_summary",
                "ligand_arpeggio_runs", "arpeggio_raw_contact_labels"}
    with c.dbconn(database) as db:
        existing = {r[0] for r in db.execute("SELECT name FROM sqlite_master WHERE type='table'")}
        missing = sorted(required - existing)
        if missing:
            raise RuntimeError("Required upstream stage tables are missing: " + ", ".join(missing))
        q = """SELECT i.ligand_instance_id FROM ligand_instances i JOIN structures s ON s.structure_id=i.structure_id
               WHERE i.curation_status='included'"""
        args = []
        if pdb_id: q += " AND UPPER(s.entry_id)=UPPER(?)"; args.append(pdb_id)
        if instance_id: q += " AND i.ligand_instance_id=?"; args.append(instance_id)
        if resume:
            q += " AND NOT EXISTS (SELECT 1 FROM protacability_warhead_linkability w WHERE w.ligand_instance_id=i.ligand_instance_id AND w.method_version=?)"
            args.append(VERSION)
        ids = [r[0] for r in db.execute(q + " ORDER BY i.ligand_instance_id", args)]
        if limit: ids = ids[:limit]
        rid = c.run_start(db, "protacability", {
            "method": VERSION, "limit": limit, "pdb_id": pdb_id, "ligand_instance_id": instance_id, "resume": resume,
            "mapping_method_version": MAPPING_VERSION,
            "functional_group_method_version": FUNCTIONAL_GROUP_VERSION,
            "geometry_method_version": GEOMETRY_VERSION,
            "target_chain_selection": "stage09_arpeggio_direct_protein_contact",
            "target_lysine_scope": "entire_selected_chain_surface",
            "ligand_to_lysine_distance_scored": False,
        })
        if ids and not resume:
            marks = ",".join("?" for _ in ids)
            for table in ("protacability_degrader_readiness", "protacability_assessment", "protacability_target_context", "protacability_warhead_linkability"):
                db.execute(f"DELETE FROM {table} WHERE method_version=? AND ligand_instance_id IN ({marks})", [VERSION, *ids])

        success = not_applicable = failures = 0
        for n, iid in enumerate(ids, 1):
            try:
                w = warhead_for_instance(db, iid)
                inv, inst, rdk = w["inv"], w["inst"], w["rdk"]
                db.execute(
                    """INSERT OR REPLACE INTO protacability_warhead_linkability(
                         run_id,ligand_instance_id,virus_name,protein_type,pdb_code,model_id,ligand_resname,ligand_chain,
                         ligand_residue_id,ligand_insertion_code,ligand_context_class,source_inventory_row_count,
                         smiles_available,representative_smiles,smiles_source,rdkit_available,rdkit_valid_smiles,
                         mol_weight,tpsa,hbd,hba,rotatable_bonds,heavy_atom_count_from_smiles,pdb_ligand_heavy_atom_count,
                         pdb_to_smiles_mapped_atom_count,functional_group_count,functional_group_types,
                         solvent_exposed_ligand_atom_count,solvent_exposed_mapped_atom_count,meaningful_contact_count,
                         strong_contact_count,contact_atom_count,strong_contact_atom_count,candidate_linker_atom_count,
                         candidate_linker_atom_ids,outward_supported_candidate_count,clear_exit_candidate_count,ligand_exit_geometry_score,direct_handle_candidate_count,
                         interaction_preservation_score,warhead_linkability_score,warhead_linkability_tier,
                         warhead_linkability_label,warhead_flags,warhead_notes,method_version)
                       VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)""",
                    (rid, iid, inv["virus_name"], inv["protein_type"], inv["pdb_code"], inv["model_id"], inv["ligand_resname"],
                     inv["ligand_chain"], inv["ligand_residue_id"], inv["ligand_insertion_code"], "candidate_small_molecule", 1,
                     int(bool(w["smiles"])), w["smiles"], inst["smiles_source"] or "", 1, rdk["rdkit_valid_smiles"],
                     rdk["mol_weight"], rdk["tpsa"], rdk["hbd"], rdk["hba"], rdk["rotatable_bonds"], rdk["heavy_atom_count_from_smiles"],
                     inv["ligand_heavy_atom_count"], w["mapped"], w["fg_count"], w["fg_types"], w["exposed"], w["exposed_mapped"],
                     w["meaningful"], w["strong"], w["contact_atoms"], w["strong_atoms"], len(w["candidates"]),
                     ";".join(str(a["atom_site_id"]) for a in w["candidates"]), w["outward_candidates"], w["clear_exit_candidates"],
                     w["exit_geometry_score"], w["direct_candidates"], w["preservation"], w["score"], w["tier"], w["label"], w["flags"], w["notes"], VERSION)
                )
                arp_run, contacting, target_context_status = contacting_target_chains(db, iid)
                chain_ids = sorted(contacting)
                target_context_notes = (
                    "Arpeggio direct protein contacts select the ligand-binding target chain(s); "
                    "lysine accessibility is then evaluated across each selected chain's entire surface; "
                    "ligand-to-lysine distance is not used."
                )
                if target_context_status == "not_applicable_no_protein_atoms":
                    target_context_notes = (
                        "No recognized protein atoms are present in the deposited model; target-side "
                        "PROTACability/readiness is not applicable. Ligand-side warhead evidence is retained."
                    )
                elif target_context_status == "not_applicable_no_contacting_protein_chain":
                    target_context_notes = (
                        "Stage 10 contains protein geometry, but the completed Stage-09 Arpeggio result does "
                        "not identify a direct ligand-contacting protein chain. Target-side PROTACability/readiness "
                        "is therefore not assigned; ligand-side warhead evidence is retained."
                    )
                db.execute(
                    """INSERT OR REPLACE INTO protacability_target_context(
                         run_id,ligand_instance_id,arpeggio_run_id,target_context_status,
                         contacting_protein_chain_count,contacting_protein_chain_ids,
                         target_chain_selection_basis,notes,method_version)
                       VALUES(?,?,?,?,?,?,?,?,?)""",
                    (rid, iid, arp_run, target_context_status, len(chain_ids), ";".join(chain_ids),
                     "stage09_arpeggio_direct_protein_contact", target_context_notes, VERSION),
                )
                if not contacting:
                    not_applicable += 1
                    if n % max(1, progress_every) == 0 or n == len(ids):
                        db.commit(); print(f"PROTACability progress: {n}/{len(ids)} success={success} target_not_applicable={not_applicable} failures={failures}", flush=True)
                    continue

                readiness_ids = []
                for chain_id in chain_ids:
                    target_evidence = contacting[chain_id]
                    chain = target_evidence["chain"]
                    contact_pair_count = int(target_evidence["contact_pair_count"])
                    protein_score, ptier, pnotes = protein_priority(chain)
                    annotation = (
                        "Candidate ligand context with solvent-accessible lysine terminal amines" if int(chain["nz_exposed_lys_count_gt_1"]) > 0
                        else "Candidate ligand context with limited terminal-amine accessibility signal"
                    )
                    db.execute(
                        """INSERT OR REPLACE INTO protacability_assessment(
                             run_id,ligand_instance_id,virus_name,protein_type,pdb_code,chain_id,model_id,chain_length_aa,
                             candidate_ligand_count,candidate_ligand_resnames,lys_count,exposed_lys_count,exposed_lys_fraction,
                             near_ligand_lys_count,near_ligand_exposed_lys_count,min_lys_ligand_distance_a,median_lys_ligand_distance_a,
                             total_sasa_a2,lysine_sasa_a2,lysine_surface_fraction,isoelectric_point,basic_fraction,acidic_fraction,
                             polar_fraction,hydrophobic_fraction,has_candidate_ligand,has_exposed_lysine,
                             has_ligand_proximal_exposed_lysine,lysine_sidechain_sasa_a2,lysine_nz_sasa_a2,nz_observed_lys_count,
                             nz_observed_lys_fraction,nz_exposed_lys_count_gt_1,nz_exposed_lys_fraction_gt_1,
                             nz_exposed_lys_count_gt_5,nz_exposed_lys_fraction_gt_5,linker_docking_site_annotation,protein_ligand_druggability_proxy_score,
                             protacability_proxy_score,protacability_tier,notes,target_chain_selection_basis,
                             ligand_target_contact_pair_count,method_version)
                           VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)""",
                        (rid, iid, chain["virus_name"], chain["protein_type"], chain["pdb_code"], chain["chain_id"], chain["model_id"],
                         chain["chain_length_aa"], 1, inv["ligand_resname"], chain["lys_count"], chain["exposed_lys_count"],
                         chain["exposed_lys_fraction"], 0, 0, None, None, chain["total_sasa_a2"],
                         chain["lysine_sasa_a2"], chain["lysine_surface_fraction"], chain["isoelectric_point"], chain["basic_fraction"],
                         chain["acidic_fraction"], chain["polar_fraction"], chain["hydrophobic_fraction"], 1,
                         int(int(chain["nz_exposed_lys_count_gt_1"] or 0) > 0), 0, chain["lysine_sidechain_sasa_a2"], chain["lysine_nz_sasa_a2"],
                         chain["nz_observed_lys_count"], chain["nz_observed_lys_fraction"], chain["nz_exposed_lys_count_gt_1"],
                         chain["nz_exposed_lys_fraction_gt_1"], chain["nz_exposed_lys_count_gt_5"], chain["nz_exposed_lys_fraction_gt_5"], annotation,
                         protein_score, protein_score, ptier,
                         pnotes + "; target chain selected by direct Stage-09 Arpeggio protein contact; structural proxy only, not an experimental PROTACability score",
                         "stage09_arpeggio_direct_protein_contact", contact_pair_count, VERSION)
                    )
                    lys_score = target_lysine_accessibility(chain)
                    exit_geometry_score, geometry_class = ligand_exit_geometry(w)
                    ternary_score = 0.0  # retained only as a compatibility field; ligand-to-lysine distance is not scored
                    short_ok = medium_ok = long_ok = 0

                    # Both pieces are required for degrader-oriented structural
                    # prioritization: a plausible ligand exit/attachment context AND
                    # an accessible target-lysine surface. Give target lysine
                    # availability substantial weight and prevent a ligand-perfect
                    # structure with no exposed target lysines from receiving a high
                    # readiness classification.
                    readiness = round(clamp(0.60 * w["score"] + 0.40 * lys_score), 2)
                    if int(chain["nz_exposed_lys_count_gt_1"] or 0) <= 0:
                        readiness = min(readiness, 34.99)
                    flags = []
                    if int(chain["lys_count"]) <= 0: flags.append("no_target_lysines_detected")
                    if int(chain["nz_observed_lys_count"]) <= 0: flags.append("no_observed_lysine_nz_atoms")
                    elif int(chain["nz_exposed_lys_count_gt_1"]) <= 0: flags.append("no_lysine_nz_atoms_with_sasa_gt_1_a2")
                    if int(w["clear_exit_candidates"]) <= 0: flags.append("no_locally_clear_outward_ligand_exit_vector")
                    if w["flags"]: flags.extend("warhead_" + x for x in w["flags"].split(";") if x)
                    cur = db.execute(
                        """INSERT OR REPLACE INTO protacability_degrader_readiness(
                             run_id,ligand_instance_id,virus_name,protein_type,pdb_code,chain_id,model_id,best_ligand_resname,
                             best_ligand_chain,best_ligand_residue_id,protein_structural_priority_score,warhead_linkability_score,
                             target_lysine_accessibility_score,ternary_geometry_cue_score,ligand_exit_geometry_score,clear_exit_candidate_count,degrader_design_readiness_score,
                             degrader_design_readiness_tier,evidence_level,best_linker_geometry_class,short_linker_geometry_feasible,
                             medium_linker_geometry_feasible,long_linker_geometry_feasible,exposed_lys_count,lys_count,exposed_lys_fraction,
                             lysine_surface_fraction,min_lys_ligand_distance_a,near_ligand_exposed_lys_count,nz_observed_lys_count,
                             nz_observed_lys_fraction,nz_exposed_lys_count_gt_1,nz_exposed_lys_fraction_gt_1,
                             nz_exposed_lys_count_gt_5,nz_exposed_lys_fraction_gt_5,candidate_ligand_resnames,
                             readiness_flags,readiness_notes,best_chain_for_instance,target_chain_selection_basis,
                             ligand_target_contact_pair_count,method_version)
                           VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)""",
                        (rid, iid, chain["virus_name"], chain["protein_type"], chain["pdb_code"], chain["chain_id"], chain["model_id"],
                         inv["ligand_resname"], inv["ligand_chain"], inv["ligand_residue_id"], protein_score, w["score"], lys_score,
                         ternary_score, exit_geometry_score, w["clear_exit_candidates"], readiness, readiness_tier(readiness), evidence_level(w["score"], lys_score, exit_geometry_score),
                         geometry_class, short_ok, medium_ok, long_ok, chain["exposed_lys_count"], chain["lys_count"], chain["exposed_lys_fraction"],
                         chain["lysine_surface_fraction"], None, 0, chain["nz_observed_lys_count"], chain["nz_observed_lys_fraction"],
                         chain["nz_exposed_lys_count_gt_1"], chain["nz_exposed_lys_fraction_gt_1"], chain["nz_exposed_lys_count_gt_5"],
                         chain["nz_exposed_lys_fraction_gt_5"], inv["ligand_resname"], ";".join(dict.fromkeys(flags)),
                         "warhead linkability is ligand-centered; target chain is selected by direct Stage-09 Arpeggio protein contact; target lysine accessibility then measures the fraction of ALL lysines across that chain surface with observed solvent-exposed NZ terminal amines plus NZ observation coverage; ligand-to-lysine distance is not used; exit-vector geometry is ligand-side and descriptive",
                         0, "stage09_arpeggio_direct_protein_contact", contact_pair_count, VERSION)
                    )
                    readiness_ids.append((cur.lastrowid, readiness))
                # Mark the highest-scoring chain for convenient structure-level summaries.
                best = db.execute(
                    """SELECT readiness_id FROM protacability_degrader_readiness
                       WHERE ligand_instance_id=? AND method_version=? ORDER BY degrader_design_readiness_score DESC,chain_id LIMIT 1""",
                    (iid, VERSION)
                ).fetchone()
                if best:
                    db.execute("UPDATE protacability_degrader_readiness SET best_chain_for_instance=1 WHERE readiness_id=?", (best[0],))
                success += 1
            except Exception as exc:
                failures += 1
                c.fail(db, rid, "protacability", f"{type(exc).__name__}: {exc}", instance_id=iid, code="protacability_exception")
            if n % max(1, progress_every) == 0 or n == len(ids):
                db.commit(); print(f"PROTACability progress: {n}/{len(ids)} success={success} target_not_applicable={not_applicable} failures={failures}", flush=True)
        c.run_end(db, rid, "completed" if failures == 0 else "partial", len(ids), success, not_applicable, failures)
    _write_report(database)
    return {"run_id": rid, "processed": len(ids), "success": success, "target_not_applicable": not_applicable, "failures": failures}


def main():
    p = argparse.ArgumentParser(description="Build occurrence-resolved PROTACability structural-priority tables.")
    p.add_argument("--database", default=str(c.ROOT / "viral_data_cif_v2.db"))
    p.add_argument("--limit", type=int)
    p.add_argument("--pdb-id")
    p.add_argument("--ligand-instance-id", type=int)
    p.add_argument("--resume", action="store_true")
    p.add_argument("--progress-every", type=int, default=100)
    a = p.parse_args()
    print(json.dumps(run(a.database, a.limit, a.pdb_id, a.ligand_instance_id, a.resume, a.progress_every), indent=2))


if __name__ == "__main__":
    main()
