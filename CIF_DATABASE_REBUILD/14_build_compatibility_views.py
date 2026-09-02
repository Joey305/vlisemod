#!/usr/bin/env python3
"""Build explicit occurrence-resolved compatibility views for V-LiSEMOD.

The new CIF-native tables remain authoritative.  This stage exposes familiar
legacy table names only where their meaning can be projected without discarding
the ligand occurrence key.  It also creates explicit ``v2_*`` views for new
code.  No ambiguous legacy rows are copied back into the canonical tables.
"""
from __future__ import annotations

import argparse
import csv
import importlib
import json
import sqlite3
from pathlib import Path

c = importlib.import_module("00_common")
VERSION = "compatibility-views-cif-v2.7"
MAPPING_VERSION = "legacy_mcs_etkdg_uff_cif_v2.5"
FUNCTIONAL_GROUP_VERSION = "rdkit-smarts-functional-groups-v2.3"
GEOMETRY_VERSION = "cif-ligand-geometry-v2.4"
PROTACABILITY_VERSION = "protacability-cif-v2.8"
ATTACHMENT_VERSION = "attachment-sites-cif-v2.6"


def validate_dependency_scripts():
    """Refuse to build views against a mismatched local analysis stack."""
    stage12 = importlib.import_module("12_build_protacability")
    stage13 = importlib.import_module("13_build_attachment_sites")
    checks = {
        "Stage 12 VERSION": (getattr(stage12, "VERSION", None), PROTACABILITY_VERSION),
        "Stage 12 MAPPING_VERSION": (getattr(stage12, "MAPPING_VERSION", None), MAPPING_VERSION),
        "Stage 12 FUNCTIONAL_GROUP_VERSION": (getattr(stage12, "FUNCTIONAL_GROUP_VERSION", None), FUNCTIONAL_GROUP_VERSION),
        "Stage 12 GEOMETRY_VERSION": (getattr(stage12, "GEOMETRY_VERSION", None), GEOMETRY_VERSION),
        "Stage 13 VERSION": (getattr(stage13, "VERSION", None), ATTACHMENT_VERSION),
        "Stage 13 expected PROTACability": (getattr(stage13, "EXPECTED_PROTACABILITY_VERSION", None), PROTACABILITY_VERSION),
        "Stage 13 expected mapping": (getattr(stage13, "EXPECTED_MAPPING_VERSION", None), MAPPING_VERSION),
        "Stage 13 expected functional groups": (getattr(stage13, "EXPECTED_FUNCTIONAL_GROUP_VERSION", None), FUNCTIONAL_GROUP_VERSION),
        "Stage 13 expected geometry": (getattr(stage13, "EXPECTED_GEOMETRY_VERSION", None), GEOMETRY_VERSION),
    }
    bad = [f"{label}: found={found!r}, expected={expected!r}" for label, (found, expected) in checks.items() if found != expected]
    if bad:
        raise RuntimeError("Stage 14 dependency mismatch; refusing to project stale views: " + "; ".join(bad))


def object_type(db, name):
    r = db.execute("SELECT type FROM sqlite_master WHERE lower(name)=lower(?)", (name,)).fetchone()
    return r[0] if r else None


def create_or_replace_view(db, name, sql, log):
    typ = object_type(db, name)
    if typ == "table":
        log.append((name, "skipped_existing_table"))
        return
    if typ == "view":
        db.execute(f'DROP VIEW "{name}"')
    db.execute(f'CREATE VIEW "{name}" AS {sql}')
    log.append((name, "created_view"))


def ensure_synonyms(db, project_root: Path):
    db.execute(
        """CREATE TABLE IF NOT EXISTS ligand_synonyms(
             ligand TEXT NOT NULL, synonym TEXT NOT NULL, source TEXT NOT NULL,
             UNIQUE(ligand,synonym))"""
    )
    before = db.execute("SELECT count(*) FROM ligand_synonyms").fetchone()[0]
    sources = []
    legacy = project_root / "viral_data.db"
    if legacy.exists():
        try:
            old = sqlite3.connect(f"file:{legacy}?mode=ro", uri=True)
            if old.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND lower(name)='ligand_synonyms'").fetchone():
                for lig, syn in old.execute("SELECT ligand,synonym FROM ligand_synonyms"):
                    lig, syn = str(lig or "").strip(), str(syn or "").strip()
                    if lig and syn:
                        db.execute("INSERT OR IGNORE INTO ligand_synonyms(ligand,synonym,source) VALUES(?,?,?)", (lig, syn, "legacy viral_data.db"))
                sources.append(str(legacy))
            old.close()
        except Exception:
            pass
    csv_candidates = [
        project_root / "LigandSynonym" / "Ligand_Synonyms_Clean.csv",
        project_root / "Ligand_Synonyms_Clean.csv",
        c.ROOT / "Ligand_Synonyms_Clean.csv",
    ]
    for path in csv_candidates:
        if not path.exists():
            continue
        try:
            with path.open(newline="", encoding="utf8", errors="replace") as fh:
                for row in csv.DictReader(fh):
                    lig = (row.get("Ligand") or row.get("ligand") or "").strip()
                    syn = (row.get("Synonym") or row.get("synonym") or "").strip()
                    if lig and syn:
                        db.execute("INSERT OR IGNORE INTO ligand_synonyms(ligand,synonym,source) VALUES(?,?,?)", (lig, syn, str(path)))
            sources.append(str(path))
        except Exception:
            pass
    after = db.execute("SELECT count(*) FROM ligand_synonyms").fetchone()[0]
    return {"rows_before": before, "rows_after": after, "sources": sources}


def run(database: str):
    validate_dependency_scripts()
    c.create_schema(database); c.dirs()
    created = []
    with c.dbconn(database) as db:
        rid = c.run_start(db, "compatibility_views", {"method": VERSION, "mapping_method_version": MAPPING_VERSION, "functional_group_method_version": FUNCTIONAL_GROUP_VERSION, "geometry_method_version": GEOMETRY_VERSION, "protacability_method_version": PROTACABILITY_VERSION, "attachment_method_version": ATTACHMENT_VERSION})
        synonym_info = ensure_synonyms(db, c.ROOT.parent)

        context_sql = """
        SELECT s.structure_id,s.entry_id AS pdb_id,
               GROUP_CONCAT(DISTINCT sc.virus_label) AS virus_name,
               GROUP_CONCAT(DISTINCT sc.protein_label) AS protein_type,
               s.source_cif_path,s.source_cif_sha256
        FROM structures s LEFT JOIN structure_classifications sc ON sc.structure_id=s.structure_id
        GROUP BY s.structure_id,s.entry_id,s.source_cif_path,s.source_cif_sha256
        """
        create_or_replace_view(db, "v2_structure_context", context_sql, created)

        create_or_replace_view(db, "v2_ligand_context", """
        SELECT ctx.virus_name,ctx.protein_type,ctx.pdb_id,i.ligand_instance_id,i.deposited_model_num AS model_id,
               i.label_comp_id AS ligand,i.auth_asym_id AS chain,i.auth_seq_id AS ligand_residue_id,
               i.insertion_code_normalized AS ligand_insertion_code,l.smiles,l.canonical_smiles,l.smiles_source,l.chemical_status,
               i.curation_status,i.curation_reason
        FROM ligand_instances i JOIN ligands l ON l.ligand_id=i.ligand_id
        JOIN v2_structure_context ctx ON ctx.structure_id=i.structure_id
        """, created)

        create_or_replace_view(db, "v2_ligand_atom_evidence", f"""
        SELECT lc.*,a.ligand_instance_atom_id,a.atom_site_id,
               COALESCE(a.auth_atom_id,a.label_atom_id) AS exact_atom,a.element,a.x,a.y,a.z,
               m.smiles_atom_index,sa.sasa_area,sa.legacy_exposed AS solvent_exposed,
               g.outward_score,g.points_away_from_pocket,g.nearest_protein_distance_a,
               g.exit_vector_x,g.exit_vector_y,g.exit_vector_z,g.forward_clearance_a,
               g.forward_obstruction_count,g.local_corridor_clear,g.forward_clearance_reaches_cap
        FROM v2_ligand_context lc
        JOIN ligand_instance_atoms a ON a.ligand_instance_id=lc.ligand_instance_id AND a.selected_conformer=1
        LEFT JOIN ligand_smiles_atom_mapping m ON m.ligand_instance_atom_id=a.ligand_instance_atom_id
          AND m.method_version='{MAPPING_VERSION}'
          AND m.run_id=(SELECT MAX(m2.run_id) FROM ligand_smiles_atom_mapping m2 WHERE m2.ligand_instance_id=lc.ligand_instance_id AND m2.method_version='{MAPPING_VERSION}')
        LEFT JOIN ligand_sasa_atoms sa ON sa.ligand_instance_atom_id=a.ligand_instance_atom_id
          AND sa.run_id=(SELECT MAX(s2.run_id) FROM ligand_sasa_atoms s2 WHERE s2.ligand_instance_id=lc.ligand_instance_id AND s2.status='complete')
        LEFT JOIN ligand_atom_geometry g ON g.ligand_instance_atom_id=a.ligand_instance_atom_id
          AND g.method_version='{GEOMETRY_VERSION}'
          AND g.run_id=(SELECT MAX(g2.run_id) FROM ligand_atom_geometry g2 WHERE g2.ligand_instance_id=lc.ligand_instance_id AND g2.method_version='{GEOMETRY_VERSION}')
        """, created)

        create_or_replace_view(db, "v2_ligand_comparison_atom_contacts", f"""
        SELECT lc.ligand,lc.canonical_smiles,lc.pdb_id,lc.model_id,lc.chain,lc.ligand_residue_id,lc.ligand_insertion_code,
               lc.ligand_instance_id,m.smiles_atom_index,
               COUNT(r.raw_contact_id) AS contact_label_count,
               COUNT(DISTINCT r.partner_identity_json) AS unique_partner_count
        FROM v2_ligand_context lc
        JOIN ligand_smiles_atom_mapping m ON m.ligand_instance_id=lc.ligand_instance_id
          AND m.method_version='{MAPPING_VERSION}'
          AND m.run_id=(SELECT MAX(m2.run_id) FROM ligand_smiles_atom_mapping m2 WHERE m2.ligand_instance_id=lc.ligand_instance_id AND m2.method_version='{MAPPING_VERSION}')
        LEFT JOIN arpeggio_raw_contact_labels r ON r.ligand_instance_atom_id=m.ligand_instance_atom_id
          AND r.filter_class='raw_environment'
          AND r.run_id=(SELECT MAX(ar.run_id) FROM ligand_arpeggio_runs ar WHERE ar.ligand_instance_id=lc.ligand_instance_id AND ar.status='completed')
        WHERE m.smiles_atom_index IS NOT NULL
        GROUP BY lc.ligand,lc.canonical_smiles,lc.pdb_id,lc.model_id,lc.chain,lc.ligand_residue_id,
                 lc.ligand_insertion_code,lc.ligand_instance_id,m.smiles_atom_index
        """, created)

        # Preserve the complete Stage-10 chain inventory separately.  This is
        # descriptive geometry and intentionally includes non-contacting chains.
        create_or_replace_view(db, "v2_all_chain_lysine_geometry", f"""
        SELECT g.*
        FROM target_chain_geometry g
        WHERE g.method_version='{GEOMETRY_VERSION}' AND g.run_id=(
            SELECT MAX(g2.run_id) FROM target_chain_geometry g2
            WHERE g2.ligand_instance_id=g.ligand_instance_id AND g2.chain_id=g.chain_id AND g2.method_version='{GEOMETRY_VERSION}'
        )
        """, created)

        # PROTACability lysine accessibility is restricted to Stage-09-confirmed
        # ligand-contacting protein chains by joining the authoritative v2.8
        # assessment rather than exposing every chain in the deposited model.
        create_or_replace_view(db, "v2_target_lysine_accessibility", f"""
        SELECT g.*,
               a.protacability_proxy_score AS surface_lysine_accessibility_score,
               a.protacability_tier AS surface_lysine_accessibility_tier,
               a.target_chain_selection_basis,
               a.ligand_target_contact_pair_count
        FROM target_chain_geometry g
        JOIN protacability_assessment a
          ON a.ligand_instance_id=g.ligand_instance_id
         AND a.chain_id=g.chain_id
         AND a.method_version='{PROTACABILITY_VERSION}'
         AND a.run_id=(SELECT MAX(a2.run_id) FROM protacability_assessment a2
                       WHERE a2.ligand_instance_id=a.ligand_instance_id
                         AND a2.chain_id=a.chain_id
                         AND a2.method_version='{PROTACABILITY_VERSION}')
        WHERE g.method_version='{GEOMETRY_VERSION}' AND g.run_id=(
            SELECT MAX(g2.run_id) FROM target_chain_geometry g2
            WHERE g2.ligand_instance_id=g.ligand_instance_id AND g2.chain_id=g.chain_id AND g2.method_version='{GEOMETRY_VERSION}'
        )
        """, created)

        create_or_replace_view(db, "v2_protacability_target_context", f"""
        SELECT t.*
        FROM protacability_target_context t
        WHERE t.method_version='{PROTACABILITY_VERSION}'
          AND t.run_id=(SELECT MAX(t2.run_id) FROM protacability_target_context t2
                        WHERE t2.ligand_instance_id=t.ligand_instance_id
                          AND t2.method_version='{PROTACABILITY_VERSION}')
        """, created)

        create_or_replace_view(db, "v2_protacability_best", f"""
        SELECT r.* FROM protacability_degrader_readiness r
        WHERE r.best_chain_for_instance=1 AND r.method_version='{PROTACABILITY_VERSION}'
          AND r.run_id=(SELECT MAX(r2.run_id) FROM protacability_degrader_readiness r2 WHERE r2.ligand_instance_id=r.ligand_instance_id AND r2.method_version='{PROTACABILITY_VERSION}')
        """, created)

        # Familiar legacy names.  Full occurrence identifiers remain available as
        # additional columns where useful; old clients can select only the columns
        # they historically used.
        create_or_replace_view(db, "Virus_Proteins", """
        SELECT sc.virus_label AS virus_name,s.entry_id AS pdb_id,sc.protein_label AS protein
        FROM structure_classifications sc JOIN structures s ON s.structure_id=sc.structure_id
        """, created)

        create_or_replace_view(db, "Ligand_Atoms_Smiles", f"""
        SELECT lc.virus_name,lc.pdb_id,lc.ligand,lc.chain,lc.ligand_residue_id AS ligand_id,
               COALESCE(lc.canonical_smiles,lc.smiles) AS smiles,
               COALESCE(f.functional_group_types,'') AS functional_groups,
               w.mol_weight AS molecular_weight,lc.ligand_instance_id,lc.model_id
        FROM v2_ligand_context lc
        LEFT JOIN ligand_functional_group_summary f ON f.ligand_instance_id=lc.ligand_instance_id
          AND f.method_version='{FUNCTIONAL_GROUP_VERSION}'
          AND f.run_id=(SELECT MAX(f2.run_id) FROM ligand_functional_group_summary f2 WHERE f2.ligand_instance_id=lc.ligand_instance_id AND f2.method_version='{FUNCTIONAL_GROUP_VERSION}' AND f2.status='complete')
        LEFT JOIN protacability_warhead_linkability w ON w.ligand_instance_id=lc.ligand_instance_id AND w.method_version='{PROTACABILITY_VERSION}'
          AND w.run_id=(SELECT MAX(w2.run_id) FROM protacability_warhead_linkability w2 WHERE w2.ligand_instance_id=lc.ligand_instance_id AND w2.method_version='{PROTACABILITY_VERSION}')
        WHERE lc.curation_status='included'
        """, created)

        create_or_replace_view(db, "Functional_GROUPED", f"""
        SELECT lc.virus_name,lc.pdb_id,lc.ligand,COALESCE(lc.canonical_smiles,lc.smiles) AS smiles,
               COALESCE(f.functional_group_types,'') AS functional_groups,lc.ligand_instance_id
        FROM v2_ligand_context lc
        LEFT JOIN ligand_functional_group_summary f ON f.ligand_instance_id=lc.ligand_instance_id
          AND f.method_version='{FUNCTIONAL_GROUP_VERSION}'
          AND f.run_id=(SELECT MAX(f2.run_id) FROM ligand_functional_group_summary f2 WHERE f2.ligand_instance_id=lc.ligand_instance_id AND f2.method_version='{FUNCTIONAL_GROUP_VERSION}' AND f2.status='complete')
        WHERE lc.curation_status='included'
        """, created)

        create_or_replace_view(db, "ligand_atoms", """
        SELECT lc.virus_name,lc.pdb_id,lc.ligand,lc.chain,a.atom_site_id AS atom_id,
               COALESCE(a.auth_atom_id,a.label_atom_id) AS exact_atom,a.element AS atom_type,a.x,a.y,a.z,
               lc.ligand_instance_id,lc.ligand_residue_id AS ligand_id,lc.model_id
        FROM v2_ligand_context lc JOIN ligand_instance_atoms a ON a.ligand_instance_id=lc.ligand_instance_id
        WHERE a.selected_conformer=1 AND lc.curation_status='included'
        """, created)

        create_or_replace_view(db, "solvent_exposed_atoms", """
        SELECT lc.virus_name,lc.pdb_id,lc.ligand,lc.chain,a.atom_site_id AS atom_id,
               COALESCE(a.auth_atom_id,a.label_atom_id) AS exact_atom,a.element AS atom_type,a.x,a.y,a.z,
               lc.ligand_instance_id,lc.ligand_residue_id AS ligand_id,lc.model_id
        FROM v2_ligand_context lc JOIN ligand_instance_atoms a ON a.ligand_instance_id=lc.ligand_instance_id
        JOIN ligand_sasa_atoms s ON s.ligand_instance_atom_id=a.ligand_instance_atom_id
          AND s.run_id=(SELECT MAX(s2.run_id) FROM ligand_sasa_atoms s2 WHERE s2.ligand_instance_id=lc.ligand_instance_id AND s2.status='complete')
        WHERE a.selected_conformer=1 AND s.legacy_exposed=1 AND lc.curation_status='included'
        """, created)

        create_or_replace_view(db, "RUPLEY_SASA_DATA", """
        SELECT lc.virus_name,lc.pdb_id,lc.ligand,lc.chain,COALESCE(a.auth_atom_id,a.label_atom_id) AS exact_atom,
               a.atom_site_id AS atom_id,s.sasa_area AS SASA_Area,lc.ligand_instance_id,lc.ligand_residue_id AS ligand_id,lc.model_id
        FROM v2_ligand_context lc JOIN ligand_instance_atoms a ON a.ligand_instance_id=lc.ligand_instance_id
        JOIN ligand_sasa_atoms s ON s.ligand_instance_atom_id=a.ligand_instance_atom_id
          AND s.run_id=(SELECT MAX(s2.run_id) FROM ligand_sasa_atoms s2 WHERE s2.ligand_instance_id=lc.ligand_instance_id AND s2.status='complete')
        WHERE a.selected_conformer=1 AND lc.curation_status='included'
        """, created)

        create_or_replace_view(db, "SMILES_MAP_PDB", f"""
        SELECT lc.virus_name,lc.pdb_id,lc.ligand,lc.chain,COALESCE(a.auth_atom_id,a.label_atom_id) AS exact_atom,
               a.atom_site_id AS atom_id,a.ligand_instance_atom_id AS atom_index,m.smiles_atom_index,
               lc.ligand_instance_id,lc.ligand_residue_id AS ligand_id,lc.model_id
        FROM v2_ligand_context lc JOIN ligand_smiles_atom_mapping m ON m.ligand_instance_id=lc.ligand_instance_id
          AND m.method_version='{MAPPING_VERSION}'
          AND m.run_id=(SELECT MAX(m2.run_id) FROM ligand_smiles_atom_mapping m2 WHERE m2.ligand_instance_id=lc.ligand_instance_id AND m2.method_version='{MAPPING_VERSION}')
        LEFT JOIN ligand_instance_atoms a ON a.ligand_instance_atom_id=m.ligand_instance_atom_id
        WHERE m.ligand_instance_atom_id IS NOT NULL AND lc.curation_status='included'
        """, created)

        create_or_replace_view(db, "Functional_Group_Atoms", f"""
        SELECT lc.virus_name,lc.pdb_id,lc.ligand,lc.chain,f.functional_group,
               f.atom_site_id AS atom_id,f.exact_atom,f.element AS atom_type,
               lc.ligand_instance_id,lc.ligand_residue_id AS ligand_id,lc.model_id,f.smiles_atom_index
        FROM v2_ligand_context lc JOIN ligand_functional_group_atoms f ON f.ligand_instance_id=lc.ligand_instance_id
          AND f.method_version='{FUNCTIONAL_GROUP_VERSION}'
          AND f.run_id=(SELECT MAX(s.run_id) FROM ligand_functional_group_summary s
                        WHERE s.ligand_instance_id=lc.ligand_instance_id
                          AND s.method_version='{FUNCTIONAL_GROUP_VERSION}'
                          AND s.status='complete')
        WHERE f.mapping_status='mapped_element_validated' AND lc.curation_status='included'
        """, created)

        create_or_replace_view(db, "Arpeggio_Contacts_Data", """
        SELECT lc.virus_name,lc.pdb_id,lc.ligand,lc.ligand_residue_id AS ligand_id,lc.chain,
               r.interaction_label AS Contact,r.distance AS Distance,
               COALESCE(a.auth_atom_id,a.label_atom_id) AS exact_atom,a.atom_site_id AS atom_id,
               json_extract(r.partner_identity_json,'$.label_comp_id') AS residue,
               json_extract(r.partner_identity_json,'$.auth_seq_id') AS residue_number,
               COALESCE(json_extract(r.partner_identity_json,'$.auth_atom_id'),json_extract(r.partner_identity_json,'$.label_atom_id')) AS residue_atom,
               json_extract(r.partner_identity_json,'$.auth_asym_id') AS residue_chain,
               lc.ligand_instance_id,lc.model_id
        FROM v2_ligand_context lc JOIN arpeggio_raw_contact_labels r ON r.ligand_instance_id=lc.ligand_instance_id
          AND r.run_id=(SELECT MAX(ar.run_id) FROM ligand_arpeggio_runs ar WHERE ar.ligand_instance_id=lc.ligand_instance_id AND ar.status='completed')
        LEFT JOIN ligand_instance_atoms a ON a.ligand_instance_atom_id=r.ligand_instance_atom_id
        WHERE r.filter_class='raw_environment' AND lc.curation_status='included'
        """, created)

        create_or_replace_view(db, "receptor_binding_pocket", f"""
        SELECT lc.virus_name,lc.pdb_id,p.residue_name AS residue,p.auth_asym_id AS residue_chain,
               p.auth_seq_id AS residue_number,COALESCE(p.auth_atom_id,p.label_atom_id) AS residue_atom,p.element AS atom_type,
               p.x,p.y,p.z,lc.ligand_instance_id,lc.ligand AS ligand,lc.chain AS ligand_chain,
               lc.ligand_residue_id AS ligand_id,p.distance_a AS distance
        FROM ligand_binding_pocket_atoms p JOIN v2_ligand_context lc ON lc.ligand_instance_id=p.ligand_instance_id
        WHERE p.method_version='{GEOMETRY_VERSION}'
          AND p.run_id=(SELECT MAX(p2.run_id) FROM ligand_binding_pocket_atoms p2 WHERE p2.ligand_instance_id=p.ligand_instance_id AND p2.method_version='{GEOMETRY_VERSION}')
        """, created)

        create_or_replace_view(db, "Covalent_Noncovalent", """
        SELECT lc.virus_name,lc.pdb_id,lc.ligand,lc.ligand_residue_id AS ligand_id,lc.chain,
               CASE WHEN EXISTS(
                 SELECT 1 FROM arpeggio_raw_contact_labels r
                 WHERE r.ligand_instance_id=lc.ligand_instance_id
                   AND r.run_id=(SELECT MAX(ar.run_id) FROM ligand_arpeggio_runs ar WHERE ar.ligand_instance_id=lc.ligand_instance_id AND ar.status='completed')
                   AND lower(r.interaction_label) LIKE '%covalent%'
               ) THEN 'Covalent' ELSE 'Noncovalent' END AS Inhibitor_Type,
               lc.ligand_instance_id,lc.model_id
        FROM v2_ligand_context lc WHERE lc.curation_status='included'
        """, created)

        # Explicit geometry-derived view; unlike the historical name this has a
        # documented meaning: outward-facing according to the pocket-centroid cue.
        create_or_replace_view(db, "distal_atoms", f"""
        SELECT lc.virus_name,lc.pdb_id,lc.ligand,lc.chain,g.atom_site_id AS atom_id,g.exact_atom,
               g.element AS atom_type,g.x,g.y,g.z,lc.ligand_instance_id,lc.ligand_residue_id AS ligand_id,
               g.outward_score
        FROM ligand_atom_geometry g JOIN v2_ligand_context lc ON lc.ligand_instance_id=g.ligand_instance_id
        WHERE g.points_away_from_pocket=1 AND g.method_version='{GEOMETRY_VERSION}'
          AND g.run_id=(SELECT MAX(g2.run_id) FROM ligand_atom_geometry g2 WHERE g2.ligand_instance_id=g.ligand_instance_id AND g2.method_version='{GEOMETRY_VERSION}')
        """, created)

        # Atom-level attachment evidence from Stage 13 v2.6.  Candidate rows
        # remain broad enough for medicinal-chemistry review; High rows are the
        # strict direct-handle + exposure + contact + outward/clear subset.
        create_or_replace_view(db, "v2_attachment_site_candidates", f"""
        SELECT a.*,lc.virus_name,lc.protein_type,lc.model_id
        FROM protacability_attachment_sites a JOIN v2_ligand_context lc ON lc.ligand_instance_id=a.ligand_instance_id
        WHERE a.candidate_attachment_atom=1 AND a.method_version='{ATTACHMENT_VERSION}'
          AND a.run_id=(SELECT MAX(a2.run_id) FROM protacability_attachment_sites a2 WHERE a2.ligand_instance_id=a.ligand_instance_id AND a2.method_version='{ATTACHMENT_VERSION}')
        """, created)

        create_or_replace_view(db, "v2_attachment_site_high_priority", f"""
        SELECT a.*,lc.virus_name,lc.protein_type,lc.model_id
        FROM protacability_attachment_sites a JOIN v2_ligand_context lc ON lc.ligand_instance_id=a.ligand_instance_id
        WHERE a.high_priority_attachment_atom=1
          AND a.direct_attachment_support=1
          AND a.atom_chemical_role='direct_attachment_atom'
          AND a.method_version='{ATTACHMENT_VERSION}'
          AND a.run_id=(SELECT MAX(a2.run_id) FROM protacability_attachment_sites a2 WHERE a2.ligand_instance_id=a.ligand_instance_id AND a2.method_version='{ATTACHMENT_VERSION}')
        """, created)

        create_or_replace_view(db, "v2_attachment_site_summary", f"""
        SELECT s.*,lc.virus_name,lc.protein_type,lc.model_id
        FROM protacability_attachment_site_summary s
        JOIN v2_ligand_context lc ON lc.ligand_instance_id=s.ligand_instance_id
        WHERE s.method_version='{ATTACHMENT_VERSION}'
          AND s.run_id=(SELECT MAX(s2.run_id) FROM protacability_attachment_site_summary s2
                        WHERE s2.ligand_instance_id=s.ligand_instance_id
                          AND s2.method_version='{ATTACHMENT_VERSION}'
                          AND s2.status='complete')
        """, created)

        c.run_end(db, rid, "completed", len(created), len(created), 0, 0)
        db.commit()

        report_rows = []
        for name, status in created:
            try:
                n = db.execute(f'SELECT count(*) FROM "{name}"').fetchone()[0] if status == "created_view" else None
            except Exception:
                n = None
            report_rows.append((name, status, n))

    lines = ["# Stage 14 compatibility-view report", "", f"* Method: {VERSION}",
             f"* Synonym rows: {synonym_info['rows_after']}",
             f"* Synonym sources: {', '.join(synonym_info['sources']) if synonym_info['sources'] else 'none found; table left available for later load'}",
             "", "## Views"]
    lines += [f"* {name}: {status}" + (f" ({n} rows)" if n is not None else "") for name, status, n in report_rows]
    lines += ["", "The CIF-native normalized tables remain authoritative. Compatibility views retain ligand_instance_id/model identifiers wherever possible."]
    (c.ROOT / "outputs" / "COMPATIBILITY_VIEW_MANIFEST.md").write_text("\n".join(lines) + "\n")
    return {"run_id": rid, "views": len(created), "synonym_rows": synonym_info["rows_after"], "details": created}


def main():
    ap = argparse.ArgumentParser(description="Create occurrence-resolved compatibility and web-facing views.")
    ap.add_argument("--database", default=str(c.ROOT / "viral_data_cif_v2.db"))
    a = ap.parse_args()
    print(json.dumps(run(a.database), indent=2))


if __name__ == "__main__":
    main()
