#!/usr/bin/env python3
"""Create order-independent fingerprints for released scientific outputs.

The fingerprints intentionally omit run IDs, timestamps, local filenames, and
other execution provenance.  They are comparisons of scientific content only.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sqlite3
from pathlib import Path


VERSION = "scientific-content-fingerprint-v1"


QUERIES = {
    "ligand_instances": """
        SELECT s.entry_id, i.deposited_model_num, i.label_asym_id, i.label_comp_id,
               i.auth_asym_id, i.auth_seq_id, i.insertion_code_normalized,
               i.identity_status, i.curation_status, i.curation_reason,
               l.canonical_smiles, l.chemical_status
        FROM ligand_instances i JOIN structures s ON s.structure_id=i.structure_id
        JOIN ligands l ON l.ligand_id=i.ligand_id
        WHERE i.curation_status='included'
    """,
    "mapping": """
        SELECT s.entry_id, i.deposited_model_num, i.label_asym_id, i.label_comp_id,
               i.auth_asym_id, i.auth_seq_id, i.insertion_code_normalized,
               m.mapping_status, m.mapped_count, m.structural_atom_count,
               m.smiles_atom_count, m.mcs_atom_count, m.mapping_outcome,
               m.mapping_reason_class, m.algorithm_queue, m.heavy_atoms_structural,
               m.heavy_atoms_reference, m.heavy_atoms_mapped,
               m.heavy_atom_mapping_fraction, m.mapping_complete,
               m.downstream_mapping_eligibility
        FROM ligand_mapping_runs m JOIN ligand_instances i ON i.ligand_instance_id=m.ligand_instance_id
        JOIN structures s ON s.structure_id=i.structure_id
        WHERE m.method_version='legacy_mcs_etkdg_uff_cif_v2.5'
    """,
    "sasa": """
        SELECT s.entry_id, i.deposited_model_num, i.label_asym_id, i.label_comp_id,
               i.auth_asym_id, i.auth_seq_id, i.insertion_code_normalized,
               a.atom_site_id, a.label_atom_id, sa.sasa_area, sa.legacy_exposed,
               sa.probe_radius, sa.point_density, sa.water_treatment,
               sa.conformer_policy, sa.status
        FROM ligand_sasa_atoms sa JOIN ligand_instance_atoms a ON a.ligand_instance_atom_id=sa.ligand_instance_atom_id
        JOIN ligand_instances i ON i.ligand_instance_id=sa.ligand_instance_id
        JOIN structures s ON s.structure_id=i.structure_id
        WHERE sa.method_version='biopython-shrake_rupley-1.40-cif-v2.1'
    """,
    "functional_groups": """
        SELECT s.entry_id, i.deposited_model_num, i.label_asym_id, i.label_comp_id,
               i.auth_asym_id, i.auth_seq_id, i.insertion_code_normalized,
               f.functional_group_match_count, f.functional_group_type_count,
               f.functional_group_types, f.direct_handle_group_count,
               f.conditional_handle_group_count, f.mapped_functional_group_atom_count,
               f.unmapped_functional_group_atom_count, f.library_source, f.status
        FROM ligand_functional_group_summary f JOIN ligand_instances i ON i.ligand_instance_id=f.ligand_instance_id
        JOIN structures s ON s.structure_id=i.structure_id
        WHERE f.method_version='rdkit-smarts-functional-groups-v2.3'
    """,
    "target_context": """
        SELECT s.entry_id, i.deposited_model_num, i.label_asym_id, i.label_comp_id,
               i.auth_asym_id, i.auth_seq_id, i.insertion_code_normalized,
               t.target_context_status, t.contacting_protein_chain_count,
               t.contacting_protein_chain_ids, t.target_chain_selection_basis, t.notes
        FROM protacability_target_context t JOIN ligand_instances i ON i.ligand_instance_id=t.ligand_instance_id
        JOIN structures s ON s.structure_id=i.structure_id
        WHERE t.method_version='protacability-cif-v2.8'
    """,
    "protacability": """
        SELECT pdb_code, chain_id, model_id, candidate_ligand_resnames,
               exposed_lys_count, nz_exposed_lys_count_gt_1,
               protein_ligand_druggability_proxy_score, protacability_proxy_score,
               protacability_tier, target_chain_selection_basis,
               ligand_target_contact_pair_count
        FROM protacability_assessment WHERE method_version='protacability-cif-v2.8'
    """,
    "attachment_sites": """
        SELECT pdb_code, ligand_resname, ligand_chain, ligand_residue_id,
               ligand_insertion_code, atom_site_id, exact_atom, element,
               smiles_atom_indices, mapped, sasa_area_a2, solvent_exposed,
               meaningful_contact_count, strong_contact_count, functional_groups,
               chemical_context, chemical_support, nearest_protein_distance_a,
               outward_score, exit_vector_clear, local_corridor_clear,
               candidate_attachment_atom, high_priority_attachment_atom,
               attachment_priority_score, attachment_priority_tier,
               atom_chemical_role, direct_attachment_support,
               conditional_substitution_support, chemical_rule_labels
        FROM protacability_attachment_sites WHERE method_version='attachment-sites-cif-v2.6'
    """,
}


def canonical(value):
    return json.dumps(list(value), ensure_ascii=False, separators=(",", ":"), default=str)


def fingerprint_rows(db: sqlite3.Connection, query: str) -> dict[str, object]:
    row_hashes = []
    for row in db.execute(query):
        row_hashes.append(hashlib.sha256(canonical(row).encode("utf-8")).hexdigest())
    aggregate = hashlib.sha256("\n".join(sorted(row_hashes)).encode("ascii")).hexdigest()
    return {"rows": len(row_hashes), "sha256": aggregate}


def build(database: Path) -> dict[str, object]:
    uri = database.resolve().as_uri() + "?mode=ro"
    with sqlite3.connect(uri, uri=True) as db:
        return {"version": VERSION, "components": {name: fingerprint_rows(db, query) for name, query in QUERIES.items()}}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--database", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--reference", help="Fail if fingerprints differ from this JSON file")
    args = parser.parse_args()
    result = build(Path(args.database))
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if args.reference:
        expected = json.loads(Path(args.reference).read_text(encoding="utf-8"))
        if result != expected:
            print("Scientific content fingerprints: DIFFER")
            return 1
    print("Scientific content fingerprints: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
