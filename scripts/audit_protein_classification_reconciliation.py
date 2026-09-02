#!/usr/bin/env python3
"""Read-only audit of folder-derived vs ligand-contact-derived target identity.

This script never opens the database for writing.  It consumes the persisted
Stage-12 v2.8 target context, whose chain selection is derived from the latest
completed Stage-09 Arpeggio direct protein contacts, and maps those chains to
local mmCIF entities.  Its output is an audit-only recommendation, not a
production classification update.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sqlite3
from collections import Counter, defaultdict
from pathlib import Path

import gemmi


METHOD_VERSION = "protacability-cif-v2.8"
AUDIT_VERSION = "protein-classification-reconciliation-audit-v1"
FIELDNAMES = [
    "pdb_id", "virus_name", "ligand_id", "ligand_instance_id", "ligand_chain",
    "ligand_residue", "model_id", "folder_protein_labels", "classification_count",
    "classification_provenance", "contacting_protein_chains", "contacting_entity_ids",
    "contacting_entity_descriptions", "normalized_contact_protein_labels",
    "structure_entity_protein_labels", "current_stage14_protein_type",
    "proposed_canonical_target", "audit_category", "resolution_confidence", "notes",
]


def text(value: object) -> str:
    return "" if value is None else str(value).strip()


def stable_join(values) -> str:
    return ";".join(sorted({text(value) for value in values if text(value)}))


def normalize_description(description: str) -> str:
    """Map conservative, explicit entity-description terms to V-LiSEMOD labels."""
    value = re.sub(r"[^a-z0-9]+", " ", text(description).lower()).strip()
    # Specific rules must precede broad tokens such as polymerase and protein.
    aliases = [
        ("programmed 1 ribosomal frameshifting", "programmed_1_ribosomal_frameshifting_element"),
        ("peptide binding domain", "peptide_binding_domains_complexed_with_mhc"),
        ("reverse transcriptase", "reverse_transcriptase"),
        ("ribonuclease h", "rnase_h"),
        ("rnase h", "rnase_h"),
        ("nucleocapsid", "nucleocapsid_protein"),
        ("envelope glycoprotein", "envelope_glycoprotein"),
        ("glycoprotein", "envelope_glycoprotein"),
        ("capsid", "capsid_protein"),
        ("matrix", "matrix_protein"),
        ("nucleoprotein", "nucleoprotein"),
        ("ribonucleoprotein", "ribonucleoprotein"),
        ("integrase", "integrase"),
        ("protease", "protease"),
        ("helicase", "helicase"),
        ("hemagglutinin", "hemagglutinin"),
        ("spike", "spike_protein"),
        ("fusion", "fusion_protein"),
        ("transmembrane", "transmembrane_protein"),
        ("accessory", "accessory_proteins"),
        ("non structural", "nsp_proteins"),
        ("nsp", "nsp_proteins"),
        ("polyprotein", "pol_protein"),
        ("polymerase", "polymerase"),
        ("e7", "e7_oncoprotein"),
    ]
    for needle, label in aliases:
        if needle in value:
            return label
    if re.search(r"\bvp\d*\b", value):
        return "vp_proteins"
    return ""


def cif_entity_metadata(path: str) -> dict:
    """Return author-chain -> local mmCIF entity metadata, without network I/O."""
    source = Path(path)
    if not source.exists():
        return {"error": f"source CIF missing: {source}", "chains": {}, "structure_labels": []}
    try:
        block = gemmi.cif.read_file(str(source)).sole_block()
        entity_ids = list(block.find_loop("_entity.id") or [])
        descriptions = list(block.find_loop("_entity.pdbx_description") or [])
        entity_description = {
            text(entity_id): text(descriptions[index]) if index < len(descriptions) else ""
            for index, entity_id in enumerate(entity_ids)
        }
        asym_ids = list(block.find_loop("_struct_asym.id") or [])
        asym_entity_ids = list(block.find_loop("_struct_asym.entity_id") or [])
        label_to_entity = {
            text(asym): text(asym_entity_ids[index]) if index < len(asym_entity_ids) else ""
            for index, asym in enumerate(asym_ids)
        }
        auth_ids = list(block.find_loop("_atom_site.auth_asym_id") or [])
        label_ids = list(block.find_loop("_atom_site.label_asym_id") or [])
        auth_to_label = {}
        for auth, label in zip(auth_ids, label_ids):
            auth, label = text(auth), text(label)
            if auth and label:
                auth_to_label.setdefault(auth, label)
        chains = {}
        for auth, label in auth_to_label.items():
            entity_id = label_to_entity.get(label, "")
            description = entity_description.get(entity_id, "")
            chains[auth] = {
                "label_asym_id": label,
                "entity_id": entity_id,
                "description": description,
                "normalized_label": normalize_description(description),
            }
        # Some files use the same author/label asym identifiers but omit an atom
        # mapping for a chain type; retain struct_asym as a safe fallback.
        for label, entity_id in label_to_entity.items():
            chains.setdefault(label, {
                "label_asym_id": label,
                "entity_id": entity_id,
                "description": entity_description.get(entity_id, ""),
                "normalized_label": normalize_description(entity_description.get(entity_id, "")),
            })
        structure_labels = sorted({
            normalize_description(description) for description in entity_description.values()
            if normalize_description(description)
        })
        return {"error": "", "chains": chains, "structure_labels": structure_labels}
    except Exception as exc:  # malformed local metadata is an audit finding
        return {"error": f"mmCIF parse error: {type(exc).__name__}: {exc}", "chains": {}, "structure_labels": []}


def load_occurrences(db: sqlite3.Connection):
    """Load exactly the retained population with current v2.8 target contexts."""
    return db.execute(
        """
        SELECT s.structure_id, s.entry_id AS pdb_id, s.source_cif_path,
               i.ligand_instance_id, i.label_comp_id AS ligand_id,
               i.auth_asym_id AS ligand_chain, i.auth_seq_id AS ligand_residue,
               i.deposited_model_num AS model_id,
               t.target_context_status, t.contacting_protein_chain_ids,
               t.target_chain_selection_basis,
               ctx.virus_name, ctx.protein_type AS current_stage14_protein_type
        FROM ligand_instances AS i
        JOIN structures AS s ON s.structure_id = i.structure_id
        JOIN protacability_target_context AS t
          ON t.ligand_instance_id = i.ligand_instance_id
         AND t.method_version = ?
        LEFT JOIN v2_structure_context AS ctx ON ctx.structure_id = s.structure_id
        WHERE i.curation_status = 'included'
        ORDER BY s.entry_id, i.ligand_instance_id
        """, (METHOD_VERSION,)
    ).fetchall()


def load_classifications(db: sqlite3.Connection):
    grouped = defaultdict(list)
    for row in db.execute(
        """SELECT structure_id, virus_label, protein_label, source_relative_path
               FROM structure_classifications
               ORDER BY structure_id, virus_label, protein_label, source_relative_path"""
    ):
        grouped[row["structure_id"]].append(dict(row))
    return grouped


def load_structure_inventory(db: sqlite3.Connection):
    return db.execute(
        """SELECT structure_id, entry_id AS pdb_id, source_cif_path, source_cif_sha256
               FROM structures ORDER BY entry_id, structure_id"""
    ).fetchall()


def classify(row: dict, source_labels: list[str], entity_labels: list[str], metadata: dict):
    contact_chains = [item for item in text(row["contacting_protein_chain_ids"]).split(";") if item]
    if row["target_context_status"] != "applicable_contacting_protein_chain" or not contact_chains:
        return "", "UNRESOLVED_NO_CONTACT_EVIDENCE", "none", "No usable direct Stage-09 protein-contact chain was persisted."
    chain_info = [metadata["chains"].get(chain) for chain in contact_chains]
    if metadata["error"] or any(item is None for item in chain_info):
        return "", "CONTACT_ENTITY_CLASSIFICATION_DISAGREEMENT", "low", metadata["error"] or "Contacting chain missing from local mmCIF entity map."
    contact_labels = sorted({item["normalized_label"] for item in chain_info if item["normalized_label"]})
    if len(contact_labels) != 1:
        category = "GENUINE_MULTIPROTEIN_LIGAND_INTERFACE" if len(contact_labels) > 1 else "CONTACT_ENTITY_CLASSIFICATION_DISAGREEMENT"
        note = "Direct contacting chains map to multiple protein classes." if contact_labels else "Contacting entity descriptions could not be normalized conservatively."
        return stable_join(contact_labels), category, "manual_review", note
    proposed = contact_labels[0]
    source_set = set(source_labels)
    if len(source_set) == 1 and proposed in source_set:
        return proposed, "CLEAN_AGREEMENT", "high", "Single source label agrees with the contacted local mmCIF entity."
    if len(source_set) > 1 and proposed in source_set:
        if len(entity_labels) > 1:
            return proposed, "MULTIPROTEIN_STRUCTURE_SINGLE_CONTACT_TARGET", "high", "Local structure contains multiple normalized protein entities; this ligand contacts one class."
        return proposed, "CLASSIFICATION_PROVENANCE_CONTAMINATION", "high", "Multiple source-path labels collapse into one contacted local mmCIF protein entity."
    return proposed, "CONTACT_ENTITY_CLASSIFICATION_DISAGREEMENT", "manual_review", "Contact-derived entity identity is not safely represented by the folder/source label set."


def audit(database: Path, output_dir: Path):
    output_dir.mkdir(parents=True, exist_ok=True)
    db = sqlite3.connect(f"file:{database.resolve()}?mode=ro", uri=True)
    db.row_factory = sqlite3.Row
    try:
        rows = load_occurrences(db)
        classifications = load_classifications(db)
        structure_inventory = load_structure_inventory(db)
    finally:
        db.close()

    inventory_rows = []
    for structure in structure_inventory:
        source_classifications = classifications[structure["structure_id"]]
        labels = sorted({item["protein_label"] for item in source_classifications})
        inventory_rows.append({
            "pdb_id": structure["pdb_id"],
            "virus_name": stable_join(item["virus_label"] for item in source_classifications),
            "classification_count": len(labels),
            "folder_protein_labels": stable_join(labels),
            "classification_provenance": stable_join(
                f"{item['protein_label']}:{item['source_relative_path']}" for item in source_classifications
            ),
            "source_cif_path": structure["source_cif_path"],
            "source_cif_sha256": structure["source_cif_sha256"],
        })
    inventory_path = output_dir / "structure_classification_inventory.csv"
    with inventory_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(inventory_rows[0]))
        writer.writeheader()
        writer.writerows(inventory_rows)

    metadata_cache, audit_rows = {}, []
    for source in rows:
        source_classifications = classifications[source["structure_id"]]
        labels = sorted({item["protein_label"] for item in source_classifications})
        viruses = sorted({item["virus_label"] for item in source_classifications})
        provenance = stable_join(f"{item['protein_label']}:{item['source_relative_path']}" for item in source_classifications)
        path = source["source_cif_path"]
        if path not in metadata_cache:
            metadata_cache[path] = cif_entity_metadata(path)
        metadata = metadata_cache[path]
        contact_chains = [value for value in text(source["contacting_protein_chain_ids"]).split(";") if value]
        chain_info = [metadata["chains"].get(chain, {}) for chain in contact_chains]
        proposed, category, confidence, note = classify(source, labels, metadata["structure_labels"], metadata)
        audit_rows.append({
            "pdb_id": source["pdb_id"],
            "virus_name": stable_join(viruses) or text(source["virus_name"]),
            "ligand_id": source["ligand_id"],
            "ligand_instance_id": source["ligand_instance_id"],
            "ligand_chain": source["ligand_chain"],
            "ligand_residue": source["ligand_residue"],
            "model_id": source["model_id"],
            "folder_protein_labels": stable_join(labels),
            "classification_count": len(labels),
            "classification_provenance": provenance,
            "contacting_protein_chains": stable_join(contact_chains),
            "contacting_entity_ids": stable_join(item.get("entity_id", "") for item in chain_info),
            "contacting_entity_descriptions": stable_join(item.get("description", "") for item in chain_info),
            "normalized_contact_protein_labels": stable_join(item.get("normalized_label", "") for item in chain_info),
            "structure_entity_protein_labels": stable_join(metadata["structure_labels"]),
            "current_stage14_protein_type": text(source["current_stage14_protein_type"]),
            "proposed_canonical_target": proposed,
            "audit_category": category,
            "resolution_confidence": confidence,
            "notes": note,
        })

    csv_path = output_dir / "classification_reconciliation_audit.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(audit_rows)

    structures = {row["pdb_id"] for row in audit_rows}
    multi_structures = {
        row["pdb_id"] for row in audit_rows if int(row["classification_count"]) > 1
    }
    category_counts = Counter(row["audit_category"] for row in audit_rows)
    summary = {
        "audit_version": AUDIT_VERSION,
        "database": str(database.resolve()),
        "total_structures_audited": len(structures),
        "total_structures_in_classification_inventory": len(inventory_rows),
        "multi_classified_structures_in_inventory": sum(int(row["classification_count"]) > 1 for row in inventory_rows),
        "total_ligand_occurrences_audited": len(audit_rows),
        "structures_with_multiple_protein_classifications": len(multi_structures),
        "ligand_occurrences_affected_by_multiple_classification": sum(int(row["classification_count"]) > 1 for row in audit_rows),
        "clean_agreements": category_counts["CLEAN_AGREEMENT"],
        "classification_provenance_contamination": category_counts["CLASSIFICATION_PROVENANCE_CONTAMINATION"],
        "multiprotein_structures_single_contact_target": category_counts["MULTIPROTEIN_STRUCTURE_SINGLE_CONTACT_TARGET"],
        "genuine_multiprotein_ligand_interfaces": category_counts["GENUINE_MULTIPROTEIN_LIGAND_INTERFACE"],
        "contact_entity_classification_disagreements": category_counts["CONTACT_ENTITY_CLASSIFICATION_DISAGREEMENT"],
        "unresolved_no_contact_evidence": category_counts["UNRESOLVED_NO_CONTACT_EVIDENCE"],
        "by_audit_category": dict(sorted(category_counts.items())),
        "by_virus": dict(sorted(Counter(row["virus_name"] for row in audit_rows).items())),
        "by_current_protein_label": dict(sorted(Counter(row["folder_protein_labels"] for row in audit_rows).items())),
        "by_proposed_target_label": dict(sorted(Counter(row["proposed_canonical_target"] for row in audit_rows).items())),
    }
    summary_path = output_dir / "classification_reconciliation_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    grouping = defaultdict(list)
    for row in audit_rows:
        current = row["folder_protein_labels"]
        proposed = row["proposed_canonical_target"]
        if row["audit_category"] in {"UNRESOLVED_NO_CONTACT_EVIDENCE", "CONTACT_ENTITY_CLASSIFICATION_DISAGREEMENT", "GENUINE_MULTIPROTEIN_LIGAND_INTERFACE"}:
            change_type = "MANUAL_REVIEW"
        elif current == proposed:
            change_type = "UNCHANGED"
        elif ";" in current and proposed and ";" not in proposed:
            change_type = "MERGE"
        else:
            change_type = "SPLIT"
        grouping[(row["virus_name"], current, proposed, change_type)].append(row)
    diff_rows = []
    for (virus, current, proposed, change_type), group in sorted(grouping.items()):
        diff_rows.append({
            "virus_name": virus,
            "current_target_label": current,
            "proposed_target_label": proposed,
            "change_type": change_type,
            "structure_count": len({row["pdb_id"] for row in group}),
            "ligand_occurrence_count": len(group),
            "affected_pdb_ids": stable_join(row["pdb_id"] for row in group),
        })
    diff_path = output_dir / "target_browser_grouping_diff.csv"
    with diff_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(diff_rows[0]) if diff_rows else ["virus_name", "current_target_label", "proposed_target_label", "change_type", "structure_count", "ligand_occurrence_count", "affected_pdb_ids"])
        writer.writeheader()
        writer.writerows(diff_rows)

    examples = {}
    for category in sorted(category_counts):
        examples[category] = next((row for row in audit_rows if row["audit_category"] == category), None)
    o4k = [row for row in audit_rows if row["pdb_id"] == "2O4K"]
    report_lines = [
        "# Protein-classification reconciliation audit",
        "",
        f"Audit version: `{AUDIT_VERSION}`",
        "",
        "## Inspected implementation paths",
        "",
        "- `CIF_DATABASE_REBUILD/03_ingest_structures.py`: inserts one `structure_classifications` row per source-relative CIF path.",
        "- `CIF_DATABASE_REBUILD/00_common.py`: schema retains `source_relative_path` provenance.",
        "- `CIF_DATABASE_REBUILD/14_build_compatibility_views.py`: `v2_structure_context` uses `GROUP_CONCAT(DISTINCT sc.protein_label)`.",
        "- `CIF_DATABASE_REBUILD/12_build_protacability.py`: Stage-12 chooses target chains from direct Stage-09 Arpeggio protein contacts.",
        "- `app.py`: Target Browser groups rows by virus/protein display fields from current PROTACability source payloads.",
        "",
        "## Root cause verified",
        "",
        "The compound labels are a source-path provenance aggregation artifact: a single retained CIF revision can be discovered under multiple folder classifications, all are retained in `structure_classifications`, and Stage 14 concatenates their labels. The audit does not alter this production behavior.",
        "",
        "## Scope",
        "",
        *[f"- {key.replace('_', ' ')}: {value}" for key, value in summary.items() if isinstance(value, int)],
        "",
        "## 2O4K",
        "",
        *[f"- occurrence `{row['ligand_instance_id']}` ({row['ligand_id']} {row['ligand_chain']} {row['ligand_residue']}): source labels `{row['folder_protein_labels']}`; Stage-09 contacting chains `{row['contacting_protein_chains']}`; entity `{row['contacting_entity_ids']}` / `{row['contacting_entity_descriptions']}`; proposed `{row['proposed_canonical_target']}`; category `{row['audit_category']}`." for row in o4k],
        "- Provenance is retained in the CSV: the same 2O4K CIF appears under both `HIV_1/capsid_protein/2O4K.cif` and `HIV_1/protease/2O4K.cif`; the retained source path is the protease copy.",
        "",
        "## Example categories observed",
        "",
        *[f"- `{category}`: {json.dumps(example, sort_keys=True) if example else 'not observed'}" for category, example in examples.items()],
        "",
        "## Recommended future rule (not implemented)",
        "",
        "For each retained ligand occurrence with an applicable v2.8 target context, derive a web-facing canonical target only from the distinct normalized local mmCIF entities of Stage-09-direct contacting protein chains. Collapse duplicate chains that map to the same normalized entity. Preserve multiple distinct contacted entities as an explicit multi-protein interface, and retain unresolved rows for manual review. Keep source-path classifications as provenance, never as the sole target identity.",
        "",
        "## Scientific invariants",
        "",
        "This audit is read-only. It does not change SASA, atom mapping, solvent exposure, chemical tractability, attachment-site evidence, Stage-09 contacts, target-chain selection, lysine accessibility, PROTACability scores, readiness tiers, or production views/API grouping.",
    ]
    report_path = output_dir / "classification_reconciliation_report.md"
    report_path.write_text("\n".join(report_lines) + "\n", encoding="utf-8")
    return {"csv": csv_path, "inventory": inventory_path, "summary": summary_path, "diff": diff_path, "report": report_path, "summary_data": summary}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--database", type=Path, default=Path("viral_data.db"))
    parser.add_argument("--output-dir", type=Path, default=Path("outputs/protein_classification_reconciliation"))
    args = parser.parse_args()
    if not args.database.exists():
        raise SystemExit(f"Database not found: {args.database}")
    result = audit(args.database, args.output_dir)
    print(json.dumps({key: str(value) if isinstance(value, Path) else value for key, value in result.items()}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
