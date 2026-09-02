#!/usr/bin/env python3
"""Create the Pass-10 canonical-target authority from the frozen Pass-9 table.

Pass 9 correctly resolved occurrence eligibility, but it retained a mixture of
generic fallback identifiers and exact virus-prefixed curator identifiers.  This
stage is deliberately upstream of the database/web application: it reconciles
those identifiers into a semantic vocabulary while preserving every source
label, eligibility decision, and scientific result as provenance.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd


VERSION = "canonical-target-authority-pass10-v1"
EXPECTED_TOTAL = 7355
EXPECTED_ELIGIBILITY = {"YES": 5414, "REVIEW": 1145, "NO": 796}


def clean(value: object) -> str:
    return "" if value is None or pd.isna(value) else str(value).strip()


def normalized(value: object) -> str:
    return re.sub(r"\s+", " ", re.sub(r"[^a-z0-9]+", " ", clean(value).casefold())).strip()


TARGET_METADATA = {
    "capsid_protein": ("capsid_protein", "capsid_protein"),
    "envelope_glycoprotein_precursor": ("Envelope glycoprotein precursor (gp160)", "envelope_glycoprotein"),
    "gp120": ("gp120", "envelope_glycoprotein"),
    "gp41": ("gp41", "envelope_glycoprotein"),
    "integrase": ("integrase", "integrase"),
    "matrix_protein": ("matrix_protein", "matrix_protein"),
    "protease": ("protease", "protease"),
    "papain_like_protease": ("Papain-like protease (nsp3)", "protease"),
    "reverse_transcriptase": ("reverse_transcriptase", "reverse_transcriptase"),
    "tat": ("Tat", "accessory_proteins"),
    "spacer_peptide_1": ("Spacer peptide 1", "maturation_peptide"),
    "hiv_1_capsid_sp1_maturation_site": ("Capsid–SP1 maturation site", "capsid_protein"),
    "spike_glycoprotein": ("Spike glycoprotein", "spike_protein"),
    "spike_protein_s1": ("Spike protein S1", "spike_protein"),
    "polymerase": ("RNA-dependent RNA polymerase (nsp12)", "polymerase"),
    "nsp1": ("Non-structural protein 1 (nsp1)", "nsp_proteins"),
    "nsp3": ("Non-structural protein 3 (nsp3)", "nsp_proteins"),
    "nsp7": ("Non-structural protein 7 (nsp7)", "nsp_proteins"),
    "nsp9": ("Non-structural protein 9 (nsp9)", "nsp_proteins"),
    "nsp10": ("Non-structural protein 10 (nsp10)", "nsp_proteins"),
    "nsp11": ("Non-structural protein 11 (nsp11)", "nsp_proteins"),
    "nsp14": ("Non-structural protein 14 (nsp14)", "nsp_proteins"),
    "nsp15": ("Uridylate-specific endoribonuclease (nsp15)", "nsp_proteins"),
    "nsp16": ("Non-structural protein 16 (nsp16)", "nsp_proteins"),
    "orf3a": ("ORF3a protein", "accessory_proteins"),
}

TARGET_ROLES = {
    **{target_id: "VIRAL_TARGET" for target_id in TARGET_METADATA},
    "hiv_1_capsid_sp1_maturation_site": "VIRAL_MATURATION_SITE",
}


DIRECT_MERGES = {
    "HIV_1": {
        "hiv_1_capsid": "capsid_protein",
        "hiv_1_gp120": "gp120",
        "hiv_1_gp41": "gp41",
        "hiv_1_integrase": "integrase",
        "hiv_1_matrix": "matrix_protein",
        "hiv_1_protease": "protease",
        "hiv_1_reverse_transcriptase": "reverse_transcriptase",
        "hiv_1_tat": "tat",
    },
    "HPV_18": {"hpv_18_l1": "capsid_protein"},
    "SARS_CoV_2": {
        "sars_cov_2_nsp3": "nsp3",
        "sars_cov_2_nsp5": "protease",
        "sars_cov_2_nsp12": "polymerase",
        "sars_cov_2_nsp15": "nsp15",
        "sars_cov_2_orf3a": "orf3a",
        "sars_cov_2_spike": "spike_glycoprotein",
    },
}


def _sars_nsp_target(description: str) -> str:
    """Resolve the former broad ``nsp_proteins`` fallback using entity text."""
    if "3c like proteinase nsp5" in description:
        return "protease"
    if "non structural protein 3" in description or "nsp3 macrodomain" in description:
        return "nsp3"
    if "nsp14" in description or "proofreading exoribonuclease" in description or "guanine n7 methyltransferase" in description:
        return "nsp14"
    if "nsp16" in description or "2 o methyltransferase" in description:
        return "nsp16"
    if "nsp1" in description or "host translation inhibitor" in description:
        return "nsp1"
    if "non structural protein 9" in description:
        return "nsp9"
    if "non structural protein 7" in description:
        return "nsp7"
    if "non structural protein 11" in description:
        return "nsp11"
    if "non structural protein 10" in description:
        return "nsp10"
    if "non structural protein 12" in description or "rna directed rna polymerase nsp12" in description:
        return "polymerase"
    raise ValueError(f"No Pass-10 SARS-CoV-2 NSP mapping for description: {description!r}")


def canonicalize_row(row: pd.Series) -> tuple[str, str]:
    """Return the Pass-10 target ID and a precise explanation of the change."""
    prior = clean(row["final_canonical_target_id"])
    virus = clean(row["virus_name"])
    description = normalized(row.get("contacting_entity_descriptions", ""))
    target_id = DIRECT_MERGES.get(virus, {}).get(prior, prior)
    reason = "direct semantic alias merge" if target_id != prior else "retained canonical identity"

    if virus == "HIV_1" and target_id == "envelope_glycoprotein":
        if "gp120" in description or "glycoprotein 120" in description:
            target_id, reason = "gp120", "split generic envelope label using gp120 entity description"
        elif "gp41" in description or "transmembrane protein gp41" in description:
            target_id, reason = "gp41", "split generic envelope label using gp41 entity description"
        elif "gp160" in description:
            target_id, reason = "envelope_glycoprotein_precursor", "retain gp160 precursor as a related, distinct target"

    if virus == "HIV_1" and target_id == "capsid_protein" and "spacer peptide 1" in description:
        if "capsid" in description or "p24" in description:
            target_id, reason = "hiv_1_capsid_sp1_maturation_site", "retain capsid-SP1 maturation junction as a distinct target"
        else:
            target_id, reason = "spacer_peptide_1", "split mature spacer peptide from capsid target"

    if virus == "SARS_CoV_2":
        if target_id == "protease":
            if "papain like protease" in description:
                target_id, reason = "papain_like_protease", "split PLpro/nsp3 from main/HIV protease class"
            elif "main protease" in description or "3c like proteinase" in description:
                target_id, reason = "protease", "merge SARS-CoV-2 main protease/nsp5 alias"
        elif target_id == "nsp_proteins":
            target_id, reason = _sars_nsp_target(description), "replace broad NSP fallback using entity description"
        elif target_id == "spike_protein" and "s1" in description:
            target_id, reason = "spike_protein_s1", "retain spike S1 subunit as related but distinct from full spike"

    return target_id, reason


def metadata_for(target_id: str, row: pd.Series) -> tuple[str, str]:
    if target_id in TARGET_METADATA:
        return TARGET_METADATA[target_id]
    return clean(row["final_canonical_target_name"]), clean(row["final_target_family"])


def vocabulary_rows(frame: pd.DataFrame) -> list[dict[str, object]]:
    result: list[dict[str, object]] = []
    for target_id, group in frame.groupby("final_canonical_target_id", sort=True):
        if not clean(target_id):
            continue
        counts = Counter(group["final_target_browser_eligible"])
        yes = group[group["final_target_browser_eligible"].eq("YES")]
        result.append({
            "canonical_target_id": target_id,
            "canonical_target_name": "; ".join(sorted({clean(v) for v in group["final_canonical_target_name"] if clean(v)})),
            "viruses": "; ".join(sorted({clean(v) for v in group["virus_name"] if clean(v)})),
            "authority_rows": len(group),
            "yes_count": counts["YES"],
            "review_count": counts["REVIEW"],
            "no_count": counts["NO"],
            "source_protein_types": "; ".join(sorted({clean(v) for v in group["current_stage14_protein_type"] if clean(v)})),
            "representative_pdbs": "; ".join(sorted({clean(v) for v in yes["pdb_id"] if clean(v)})[:10]),
            "representative_ligand_instances": "; ".join(sorted({clean(v) for v in yes["ligand_instance_id"] if clean(v)})[:10]),
        })
    return result


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def run(pass9_csv: Path, outdir: Path) -> dict[str, object]:
    frame = pd.read_csv(pass9_csv, dtype=str, low_memory=False).fillna("")
    if len(frame) != EXPECTED_TOTAL:
        raise RuntimeError(f"Pass-9 input rows={len(frame)}; expected {EXPECTED_TOTAL}")
    if Counter(frame["final_target_browser_eligible"]) != EXPECTED_ELIGIBILITY:
        raise RuntimeError("Pass-9 eligibility partition is not the protected baseline")

    out = frame.copy()
    yes_mask = out["final_target_browser_eligible"].eq("YES")
    prior_ids = out.loc[yes_mask, "final_canonical_target_id"].copy()
    changes: list[dict[str, str]] = []
    for idx in out.index[yes_mask]:
        prior = clean(out.at[idx, "final_canonical_target_id"])
        target_id, reason = canonicalize_row(out.loc[idx])
        target_name, target_family = metadata_for(target_id, out.loc[idx])
        out.at[idx, "final_canonical_target_id"] = target_id
        out.at[idx, "final_canonical_target_name"] = target_name
        out.at[idx, "final_target_family"] = target_family
        if target_id in TARGET_ROLES:
            out.at[idx, "final_entity_role"] = TARGET_ROLES[target_id]
        if target_id != prior:
            out.at[idx, "final_authority_basis"] = clean(out.at[idx, "final_authority_basis"]) + ";PASS10_SEMANTIC_VOCABULARY_NORMALIZATION"
            note = clean(out.at[idx, "final_authority_note"])
            out.at[idx, "final_authority_note"] = (note + " | " if note else "") + reason
            changes.append({
                "virus_name": clean(out.at[idx, "virus_name"]),
                "from_canonical_target_id": prior,
                "to_canonical_target_id": target_id,
                "reason": reason,
                "pdb_id": clean(out.at[idx, "pdb_id"]),
                "ligand_instance_id": clean(out.at[idx, "ligand_instance_id"]),
            })

    eligibility = Counter(out["final_target_browser_eligible"])
    if eligibility != EXPECTED_ELIGIBILITY:
        raise RuntimeError(f"normalization changed eligibility: {dict(eligibility)}")
    if out.loc[yes_mask, "final_canonical_target_id"].map(clean).eq("").any():
        raise RuntimeError("Pass-10 created a YES row without a canonical target ID")
    if out.loc[yes_mask, "final_canonical_target_id"].eq("nsp_proteins").any():
        raise RuntimeError("broad nsp_proteins fallback survived Pass-10 normalization")

    # A name is display-only, but duplicate IDs for the same virus and display
    # concept are synonym fragmentation and must be rejected at generation time.
    identity_conflicts = (
        out.loc[yes_mask]
        .groupby(["virus_name", "final_canonical_target_name"])["final_canonical_target_id"]
        .nunique()
    )
    conflicts = identity_conflicts[identity_conflicts.gt(1)]
    if not conflicts.empty:
        raise RuntimeError(f"same virus/display concept maps to multiple IDs: {conflicts.to_dict()}")

    outdir.mkdir(parents=True, exist_ok=True)
    authority_csv = outdir / "canonical_target_occurrences_pass10.csv"
    out.to_csv(authority_csv, index=False)
    out[out["final_target_browser_eligible"].eq("REVIEW")].to_csv(outdir / "canonical_target_review_queue_pass10.csv", index=False)
    out[out["final_target_browser_eligible"].eq("NO")].to_csv(outdir / "canonical_target_excluded_contexts_pass10.csv", index=False)

    yes = out.loc[yes_mask].copy()
    group_columns = ["virus_name", "final_canonical_target_id"]
    groups = yes.groupby(group_columns, dropna=False).agg(
        final_canonical_target_name=("final_canonical_target_name", "first"),
        final_target_family=("final_target_family", "first"),
        final_entity_role=("final_entity_role", "first"),
        occurrence_count=("ligand_instance_id", "size"),
        structure_count=("pdb_id", "nunique"),
        ligand_count=("ligand_id", "nunique"),
    ).reset_index().sort_values(group_columns)
    groups.to_csv(outdir / "canonical_target_browser_groups_pass10.csv", index=False)

    change_frame = pd.DataFrame(changes)
    if not change_frame.empty:
        change_frame = change_frame.sort_values(["virus_name", "from_canonical_target_id", "to_canonical_target_id", "ligand_instance_id"])
    change_frame.to_csv(outdir / "canonical_target_semantic_changes_pass10.csv", index=False)
    pd.DataFrame(vocabulary_rows(out)).to_csv(outdir / "canonical_target_vocabulary_audit_pass10.csv", index=False)

    old_group_count = int(frame.loc[frame["final_target_browser_eligible"].eq("YES"), ["virus_name", "final_canonical_target_id"]].drop_duplicates().shape[0])
    summary = {
        "authority_version": VERSION,
        "input_pass9_csv": str(pass9_csv.resolve()),
        "input_pass9_sha256": sha256_file(pass9_csv),
        "authority_rows": len(out),
        "eligibility_counts": dict(sorted(eligibility.items())),
        "eligible_occurrence_rows": int(yes_mask.sum()),
        "old_distinct_canonical_target_ids": int(prior_ids.nunique()),
        "new_distinct_canonical_target_ids": int(yes["final_canonical_target_id"].nunique()),
        "old_virus_target_groups": old_group_count,
        "new_virus_target_groups": int(len(groups)),
        "semantic_change_rows": len(changes),
        "semantic_changes_by_pair": {
            f"{key[0]}::{key[1]}->{key[2]}": int(value)
            for key, value in Counter((c["virus_name"], c["from_canonical_target_id"], c["to_canonical_target_id"]) for c in changes).items()
        },
        "protected_scientific_data_modified": False,
    }
    (outdir / "canonical_target_authority_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {**summary, "authority_csv": str(authority_csv.resolve()), "authority_csv_sha256": sha256_file(authority_csv)}


def main() -> None:
    parser = argparse.ArgumentParser(description="Normalize the Pass-9 canonical target vocabulary into Pass 10.")
    parser.add_argument("--pass9-csv", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(run(args.pass9_csv, args.outdir), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
