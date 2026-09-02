#!/usr/bin/env python3
"""Materialize atom-level attachment-site evidence with atom-specific chemistry.

Stage 13 is the ligand-side modification-site synthesis layer. It combines the
validated occurrence-resolved evidence from Stages 07-11/12, but it does not
use target-lysine proximity or degrader-readiness scores.

Why v2.6
--------
Stage 11 deliberately assigns broad functional-group *context*. A SMARTS match
such as pyridine, ester, carbamate, ether, or aryl halide can contain several
atoms, but not every atom in that match is itself a defensible linker-attachment
reaction center. Stage 13 v2.6 therefore keeps the Stage-11 context and adds an
explicit atom-specific medicinal-chemistry role derived in the same source-
SMILES atom-index namespace used by Stages 07 and 11.

Atom-specific chemical roles
----------------------------
``direct_attachment_atom``
    A common, identifiable atom-level handle for single-step derivatization
    (for example a free primary/secondary amine N, alcohol/phenol O, thiol S,
    carboxylic-acid acyl carbon, terminal-alkyne carbon, or terminal azide N).

``conditional_substitution_site``
    A chemically plausible atom-level modification site that normally requires
    substitution/coupling, precursor redesign, or a charge-/function-changing
    transformation (for example pyridine N, an exposed C-H position within a
    Stage-11 pyridine/pyrrole context, C-X substitution/coupling sites,
    boron-bearing carbon, amide/sulfonamide N-H, or aldehyde/ketone carbonyl C).

``functional_group_context_only``
    The atom belongs to an informative Stage-11 functional group, but the atom
    itself is not automatically treated as a reaction center. Saturated ether
    O atoms and generic atoms inside ester/carbamate/amide SMARTS matches fall
    here unless an explicit atom-specific rule applies.

``unclassified_atom_context``
    No mapped Stage-11 functional-group context and no atom-specific rule.

Scientific gate
---------------
A *High attachment-site priority* atom must satisfy all of the following:
  1. deposited heavy atom;
  2. mapped through Stage 07 v2.5;
  3. solvent exposed by Stage 08;
  4. no strong non-water Arpeggio interaction at that ligand atom;
  5. atom-specific role == ``direct_attachment_atom``;
  6. points away from the local binding-pocket centroid in Stage 10 v2.4; and
  7. has a locally clear forward corridor in Stage 10 v2.4.

Conditional substitution sites remain useful medicinal-chemistry hypotheses,
but are capped at Moderate priority so that "High" is reserved for a direct
atom-level handle plus favorable structural geometry. Functional-group context
alone can never create a Moderate/High chemical site.

The forward corridor is only a local geometric cue. It is not a linker path,
ternary-complex model, SAR result, or proof that modification will preserve
binding or cause degradation.
"""
from __future__ import annotations

import argparse
import importlib
import json
from collections import defaultdict

from rdkit import Chem, RDLogger

c = importlib.import_module("00_common")
p = importlib.import_module("12_build_protacability")

VERSION = "attachment-sites-cif-v2.6"
EXPECTED_PROTACABILITY_VERSION = "protacability-cif-v2.8"
EXPECTED_MAPPING_VERSION = "legacy_mcs_etkdg_uff_cif_v2.5"
EXPECTED_FUNCTIONAL_GROUP_VERSION = "rdkit-smarts-functional-groups-v2.3"
EXPECTED_GEOMETRY_VERSION = "cif-ligand-geometry-v2.4"

ROLE_DIRECT = "direct_attachment_atom"
ROLE_CONDITIONAL = "conditional_substitution_site"
ROLE_CONTEXT = "functional_group_context_only"
ROLE_UNCLASSIFIED = "unclassified_atom_context"

ROLE_RANK = {
    ROLE_UNCLASSIFIED: 0,
    ROLE_CONTEXT: 1,
    ROLE_CONDITIONAL: 2,
    ROLE_DIRECT: 3,
}

# Each SMARTS identifies a *reaction-center atom*, not an entire functional
# group. ``anchor`` is the atom position within the SMARTS match that receives
# the role. These are deliberately conservative and auditable; they are not a
# reaction-prediction engine or a guarantee of synthetic feasibility.
ATOM_TRACTABILITY_RULES = (
    # Direct atom-level handles.
    {
        "role": ROLE_DIRECT,
        "label": "free_primary_secondary_amine_n",
        "smarts": "[N;X3;H1,H2;+0;!$(N-C=O);!$(N-S(=O)=O);!$(N-C(=N)-N)]",
        "anchor": 0,
        "rationale": "free primary/secondary amine N; common N-derivatization handle",
    },
    {
        "role": ROLE_DIRECT,
        "label": "alcohol_or_phenol_o",
        "smarts": "[O;X2;H1;+0;!$([O;H1]-[C,S,P]=O)]",
        "anchor": 0,
        "rationale": "neutral O-H atom outside carboxylic/sulfonic/phosphoric acid; common O-derivatization handle",
    },
    {
        "role": ROLE_DIRECT,
        "label": "thiol_s",
        "smarts": "[S;X2;H1;+0]",
        "anchor": 0,
        "rationale": "thiol S-H atom; common S-derivatization handle",
    },
    {
        "role": ROLE_DIRECT,
        "label": "carboxylic_acid_acyl_c",
        "smarts": "[C;X3](=[O;X1])[O;X2;H1]",
        "anchor": 0,
        "rationale": "carboxylic-acid acyl carbon; linker installation can be tested by amide/ester coupling",
    },
    {
        "role": ROLE_DIRECT,
        "label": "terminal_alkyne_c",
        "smarts": "[C;X2;H1]#[C;X2]",
        "anchor": 0,
        "rationale": "terminal alkyne carbon; common click-compatible handle",
    },

    # Conditional atom-level sites.
    {
        "role": ROLE_CONDITIONAL,
        "label": "tertiary_amine_n_charge_changing",
        "smarts": "[N;X3;H0;+0;!$(N-C=O);!$(N-S(=O)=O)]",
        "anchor": 0,
        "rationale": "neutral tertiary amine N; N-alkylation is possible but changes charge/substitution state",
    },
    {
        "role": ROLE_CONDITIONAL,
        "label": "pyridine_like_aromatic_n_charge_changing",
        "smarts": "[n;H0;+0]",
        "anchor": 0,
        "rationale": "pyridine-like aromatic N; N-functionalization is possible but normally changes charge/electronics",
    },
    {
        "role": ROLE_CONDITIONAL,
        "label": "amide_or_carbamate_nh",
        "smarts": "[N;X3;H1,H2][C;X3](=[O;X1])",
        "anchor": 0,
        "rationale": "amide/carbamate N-H; N-derivatization is possible but can alter resonance and binding",
    },
    {
        "role": ROLE_CONDITIONAL,
        "label": "sulfonamide_nh",
        "smarts": "[N;X3;H1,H2][S;X4](=[O;X1])(=[O;X1])",
        "anchor": 0,
        "rationale": "sulfonamide N-H; N-derivatization is possible but may alter acidity/resonance",
    },
    {
        "role": ROLE_CONDITIONAL,
        "label": "aryl_halide_carbon",
        "smarts": "[c]-[F,Cl,Br,I]",
        "anchor": 0,
        "rationale": "aryl carbon bearing halogen; atom-specific substitution/cross-coupling site",
    },
    {
        "role": ROLE_CONDITIONAL,
        "label": "alkyl_halide_carbon",
        "smarts": "[C;X4]-[F,Cl,Br,I]",
        "anchor": 0,
        "rationale": "sp3 carbon bearing halogen; atom-specific substitution site",
    },
    {
        "role": ROLE_CONDITIONAL,
        "label": "boron_bearing_carbon",
        "smarts": "[#6]-[B;X3]([O])[O]",
        "anchor": 0,
        "rationale": "carbon bonded to boron; coupling can replace the C-B relationship with a linker-bearing C-C/C-X bond",
    },
    {
        "role": ROLE_CONDITIONAL,
        "label": "terminal_alkene_c",
        "smarts": "[C;X3;H2]=[C;X3]",
        "anchor": 0,
        "rationale": "terminal alkene carbon; chemically addressable but transformation-dependent",
    },
    {
        "role": ROLE_CONDITIONAL,
        "label": "aldehyde_carbonyl_c",
        "smarts": "[C;X3;H1](=[O;X1])",
        "anchor": 0,
        "rationale": "aldehyde carbonyl carbon; condensation/reductive-amination style linker installation is transformation-dependent",
    },
    {
        "role": ROLE_CONDITIONAL,
        "label": "ketone_carbonyl_c",
        "smarts": "[#6][C;X3](=[O;X1])[#6]",
        "anchor": 1,
        "rationale": "ketone carbonyl carbon; condensation/reductive-amination style linker installation is transformation-dependent",
    },
    {
        "role": ROLE_CONDITIONAL,
        "label": "oxime_oh_o",
        "smarts": "[C;X3]=[N;X2][O;X2;H1]",
        "anchor": 2,
        "rationale": "oxime O-H atom; O-functionalization is possible but transformation/context dependent",
    },
    {
        "role": ROLE_CONDITIONAL,
        "label": "acyl_halide_acyl_c",
        "smarts": "[C;X3](=[O;X1])[F,Cl,Br,I]",
        "anchor": 0,
        "rationale": "acyl-halide carbonyl carbon; reactive substitution site requiring chemical-context review",
    },
    {
        "role": ROLE_CONDITIONAL,
        "label": "isocyanate_c",
        "smarts": "[N;X2]=[C;X2]=[O;X1]",
        "anchor": 1,
        "rationale": "isocyanate carbon; electrophilic derivatization site requiring chemical-context review",
    },
    {
        "role": ROLE_CONDITIONAL,
        "label": "isothiocyanate_c",
        "smarts": "[N;X2]=[C;X2]=[S;X1]",
        "anchor": 1,
        "rationale": "isothiocyanate carbon; electrophilic derivatization site requiring chemical-context review",
    },
)


def _compile_rules():
    compiled = []
    bad = []
    for rule in ATOM_TRACTABILITY_RULES:
        patt = Chem.MolFromSmarts(rule["smarts"])
        if patt is None:
            bad.append(f"{rule['label']}={rule['smarts']}")
            continue
        if rule["anchor"] < 0 or rule["anchor"] >= patt.GetNumAtoms():
            bad.append(f"{rule['label']}:anchor={rule['anchor']}:pattern_atoms={patt.GetNumAtoms()}")
            continue
        compiled.append((rule, patt))
    if bad:
        raise RuntimeError("Invalid Stage-13 atom tractability rules: " + "; ".join(bad))
    return tuple(compiled)


COMPILED_ATOM_TRACTABILITY_RULES = _compile_rules()


def _validate_dependencies() -> None:
    observed = {
        "Stage12": getattr(p, "VERSION", None),
        "mapping": getattr(p, "MAPPING_VERSION", None),
        "functional_groups": getattr(p, "FUNCTIONAL_GROUP_VERSION", None),
        "geometry": getattr(p, "GEOMETRY_VERSION", None),
    }
    expected = {
        "Stage12": EXPECTED_PROTACABILITY_VERSION,
        "mapping": EXPECTED_MAPPING_VERSION,
        "functional_groups": EXPECTED_FUNCTIONAL_GROUP_VERSION,
        "geometry": EXPECTED_GEOMETRY_VERSION,
    }
    mismatches = [
        f"{name}: expected={expected[name]!r} observed={observed[name]!r}"
        for name in expected
        if observed[name] != expected[name]
    ]
    if mismatches:
        raise RuntimeError(
            "Stage 13 dependency mismatch; refusing to mix evidence generations: "
            + "; ".join(mismatches)
        )


def ensure_schema(database: str) -> None:
    c.create_schema(database)
    with c.dbconn(database) as db:
        db.executescript(
            """
            CREATE TABLE IF NOT EXISTS protacability_attachment_sites (
                attachment_site_id INTEGER PRIMARY KEY,
                run_id INTEGER NOT NULL REFERENCES analysis_runs(run_id),
                ligand_instance_id INTEGER NOT NULL REFERENCES ligand_instances(ligand_instance_id),
                ligand_instance_atom_id INTEGER NOT NULL REFERENCES ligand_instance_atoms(ligand_instance_atom_id),
                pdb_code TEXT NOT NULL,
                ligand_resname TEXT NOT NULL,
                ligand_chain TEXT NOT NULL,
                ligand_residue_id TEXT NOT NULL,
                ligand_insertion_code TEXT,
                atom_site_id TEXT,
                exact_atom TEXT,
                element TEXT,
                smiles_atom_indices TEXT NOT NULL,
                mapped INTEGER NOT NULL,
                sasa_area_a2 REAL,
                solvent_exposed INTEGER NOT NULL,
                meaningful_contact_count INTEGER NOT NULL,
                strong_contact_count INTEGER NOT NULL,
                unique_partner_count INTEGER NOT NULL,
                contact_labels TEXT NOT NULL,
                functional_groups TEXT NOT NULL,
                chemical_context TEXT NOT NULL,
                chemical_support INTEGER NOT NULL,
                atom_chemical_role TEXT NOT NULL DEFAULT 'unclassified_atom_context',
                direct_attachment_support INTEGER NOT NULL DEFAULT 0,
                conditional_substitution_support INTEGER NOT NULL DEFAULT 0,
                chemical_rule_labels TEXT NOT NULL DEFAULT '',
                chemical_rationale TEXT NOT NULL DEFAULT '',
                nearest_protein_distance_a REAL,
                outward_score REAL,
                points_away_from_pocket INTEGER,
                forward_clearance_a REAL,
                forward_obstruction_count INTEGER,
                exit_vector_clear INTEGER,
                local_corridor_clear INTEGER,
                forward_clearance_reaches_cap INTEGER,
                candidate_attachment_atom INTEGER NOT NULL,
                high_priority_attachment_atom INTEGER NOT NULL,
                attachment_priority_score REAL NOT NULL,
                attachment_priority_tier TEXT NOT NULL,
                method_version TEXT NOT NULL,
                UNIQUE(ligand_instance_atom_id,method_version)
            );
            CREATE TABLE IF NOT EXISTS protacability_attachment_site_summary (
                attachment_summary_id INTEGER PRIMARY KEY,
                run_id INTEGER NOT NULL REFERENCES analysis_runs(run_id),
                ligand_instance_id INTEGER NOT NULL REFERENCES ligand_instances(ligand_instance_id),
                pdb_code TEXT NOT NULL,
                ligand_resname TEXT NOT NULL,
                mapped_atom_count INTEGER NOT NULL,
                exposed_mapped_atom_count INTEGER NOT NULL,
                candidate_attachment_atom_count INTEGER NOT NULL,
                outward_supported_candidate_count INTEGER NOT NULL,
                clear_exit_supported_candidate_count INTEGER NOT NULL DEFAULT 0,
                chemically_supported_candidate_count INTEGER NOT NULL,
                direct_attachment_candidate_count INTEGER NOT NULL DEFAULT 0,
                conditional_substitution_candidate_count INTEGER NOT NULL DEFAULT 0,
                conditional_clear_exit_candidate_count INTEGER NOT NULL DEFAULT 0,
                high_priority_attachment_atom_count INTEGER NOT NULL,
                high_priority_direct_attachment_atom_count INTEGER NOT NULL DEFAULT 0,
                top_attachment_site_score REAL,
                top_attachment_atom_site_id TEXT,
                top_attachment_exact_atom TEXT,
                method_version TEXT NOT NULL,
                status TEXT NOT NULL,
                UNIQUE(ligand_instance_id,method_version)
            );
            """
        )
        migrations = {
            "protacability_attachment_sites": (
                ("forward_clearance_a", "REAL"),
                ("forward_obstruction_count", "INTEGER"),
                ("exit_vector_clear", "INTEGER"),
                ("local_corridor_clear", "INTEGER"),
                ("forward_clearance_reaches_cap", "INTEGER"),
                ("atom_chemical_role", "TEXT NOT NULL DEFAULT 'unclassified_atom_context'"),
                ("direct_attachment_support", "INTEGER NOT NULL DEFAULT 0"),
                ("conditional_substitution_support", "INTEGER NOT NULL DEFAULT 0"),
                ("chemical_rule_labels", "TEXT NOT NULL DEFAULT ''"),
                ("chemical_rationale", "TEXT NOT NULL DEFAULT ''"),
            ),
            "protacability_attachment_site_summary": (
                ("clear_exit_supported_candidate_count", "INTEGER NOT NULL DEFAULT 0"),
                ("direct_attachment_candidate_count", "INTEGER NOT NULL DEFAULT 0"),
                ("conditional_substitution_candidate_count", "INTEGER NOT NULL DEFAULT 0"),
                ("conditional_clear_exit_candidate_count", "INTEGER NOT NULL DEFAULT 0"),
                ("high_priority_direct_attachment_atom_count", "INTEGER NOT NULL DEFAULT 0"),
            ),
        }
        for table, columns in migrations.items():
            existing = {r[1] for r in db.execute(f"PRAGMA table_info({table})")}
            for column, kind in columns:
                if column not in existing:
                    db.execute(f"ALTER TABLE {table} ADD COLUMN {column} {kind}")


def _promote_role(record, role: str, label: str, rationale: str) -> None:
    if ROLE_RANK[role] > ROLE_RANK[record["role"]]:
        record["role"] = role
    record["labels"].add(label)
    record["rationales"].add(rationale)


def _atom_roles_for_molecule(mol, broad_groups_by_index: dict[int, set[str]]):
    """Return atom-specific chemical roles in the source-SMILES index space."""
    records = {
        idx: {
            "role": ROLE_CONTEXT if broad_groups_by_index.get(idx) else ROLE_UNCLASSIFIED,
            "labels": set(),
            "rationales": set(),
        }
        for idx in range(mol.GetNumAtoms())
    }

    # Explicit reaction-center SMARTS.
    for rule, patt in COMPILED_ATOM_TRACTABILITY_RULES:
        for match in mol.GetSubstructMatches(patt, uniquify=True):
            idx = int(match[rule["anchor"]])
            _promote_role(records[idx], rule["role"], rule["label"], rule["rationale"])

    # Organic azide: the terminal anionic N is the atom that forms the new bond
    # in common azide-alkyne click products.  Select it by graph/charge rather
    # than relying on one particular resonance SMARTS ordering.
    for atom in mol.GetAtoms():
        if (
            atom.GetSymbol() == "N"
            and atom.GetFormalCharge() == -1
            and atom.GetDegree() == 1
            and atom.GetNeighbors()
            and atom.GetNeighbors()[0].GetSymbol() == "N"
        ):
            _promote_role(
                records[atom.GetIdx()], ROLE_DIRECT, "terminal_azide_n",
                "terminal azide N; common click-compatible atom-level handle",
            )

    # Stage 11 classifies pyridine/pyrrole as conditional handle *contexts*.
    # At atom level, an aromatic C-H position in those heteroaromatic rings is a
    # plausible analogue-design/substitution site, not a direct existing handle.
    # This is intentionally not generalized to every benzene C-H atom.
    for idx, groups in broad_groups_by_index.items():
        if idx < 0 or idx >= mol.GetNumAtoms():
            continue
        atom = mol.GetAtomWithIdx(idx)
        if (
            atom.GetSymbol() == "C"
            and atom.GetIsAromatic()
            and atom.GetTotalNumHs(includeNeighbors=True) > 0
            and ({"Pyridine", "Pyrrole"} & groups)
        ):
            _promote_role(
                records[idx], ROLE_CONDITIONAL, "heteroaromatic_c_h_precursor_substitution",
                "heteroaromatic C-H position; plausible analogue-design/substitution site but not a pre-existing direct handle",
            )

        # Pyrrolic N-H can be derivatized, but its aromatic/electronic role makes
        # this conditional rather than a direct-handle call.
        if (
            atom.GetSymbol() == "N"
            and atom.GetIsAromatic()
            and atom.GetTotalNumHs(includeNeighbors=True) > 0
            and "Pyrrole" in groups
        ):
            _promote_role(
                records[idx], ROLE_CONDITIONAL, "pyrrolic_n_h_derivatization",
                "pyrrolic N-H; N-derivatization is possible but changes heteroaromatic electronics",
            )

    return records


def _load_atom_specific_chemistry(db, iid: int, atom_evidence: dict[int, dict]):
    """Map source-SMILES reaction-center roles back to deposited ligand atoms."""
    row = db.execute(
        """SELECT l.smiles,l.chemical_status
           FROM ligand_instances i JOIN ligands l ON l.ligand_id=i.ligand_id
           WHERE i.ligand_instance_id=?""",
        (iid,),
    ).fetchone()

    # Preserve broad Stage-11 context even for intentionally unresolved
    # chemistry. Such atoms cannot receive direct/conditional support here.
    default = {}
    for aid, atom in atom_evidence.items():
        default[aid] = {
            "atom_chemical_role": ROLE_CONTEXT if atom["functional_groups"] else ROLE_UNCLASSIFIED,
            "direct_attachment_support": 0,
            "conditional_substitution_support": 0,
            "chemical_rule_labels": "",
            "chemical_rationale": "",
        }

    if row is None or row["chemical_status"] != "resolved" or not row["smiles"]:
        return default

    mol = Chem.MolFromSmiles(row["smiles"])
    if mol is None:
        raise ValueError("resolved ligand has invalid source SMILES")

    broad_groups_by_index = defaultdict(set)
    for atom in atom_evidence.values():
        for idx in atom["smiles_atom_indices"]:
            if idx < 0 or idx >= mol.GetNumAtoms():
                raise ValueError(f"source smiles atom index out of range: {idx}/{mol.GetNumAtoms()}")
            broad_groups_by_index[int(idx)].update(atom["functional_groups"])

    by_index = _atom_roles_for_molecule(mol, broad_groups_by_index)
    result = {}
    for aid, atom in atom_evidence.items():
        merged_role = ROLE_CONTEXT if atom["functional_groups"] else ROLE_UNCLASSIFIED
        labels = set()
        rationales = set()
        for idx in sorted(atom["smiles_atom_indices"]):
            rec = by_index[int(idx)]
            if ROLE_RANK[rec["role"]] > ROLE_RANK[merged_role]:
                merged_role = rec["role"]
            labels.update(rec["labels"])
            rationales.update(rec["rationales"])
        result[aid] = {
            "atom_chemical_role": merged_role,
            "direct_attachment_support": int(merged_role == ROLE_DIRECT),
            "conditional_substitution_support": int(merged_role == ROLE_CONDITIONAL),
            "chemical_rule_labels": ";".join(sorted(labels)),
            "chemical_rationale": "; ".join(sorted(rationales)),
        }
    return result


def _attachment_rows_for_instance(db, iid: int):
    """Build Stage-13 atom rows using atom-specific chemistry plus structural evidence."""
    evidence = p.load_atom_evidence(db, iid)
    chemistry = _load_atom_specific_chemistry(db, iid, evidence)
    rows = []

    for atom in evidence.values():
        core = bool(p.candidate_core(atom))
        broad_chem = p.chemical_context(atom)
        chem = chemistry[atom["ligand_instance_atom_id"]]
        atom_role = chem["atom_chemical_role"]
        chemical_support = atom_role in {ROLE_DIRECT, ROLE_CONDITIONAL}
        outward = atom["points_away"]
        clear_exit = atom["local_corridor_clear"]

        # "High" is intentionally reserved for a direct atom-level handle.
        high_priority = bool(
            core
            and atom_role == ROLE_DIRECT
            and outward == 1
            and clear_exit == 1
        )
        conditional_clear = bool(
            core
            and atom_role == ROLE_CONDITIONAL
            and outward == 1
            and clear_exit == 1
        )

        # Transparent additive structural evidence score. Tier caps are the
        # semantic safeguard; they prevent context-only or conditional sites
        # from being mislabeled High.
        score = 0.0
        if atom["mapped"]:
            score += 20
        if atom["exposed"]:
            score += 30
        if atom["strong_contact_count"] == 0:
            score += 20

        if atom["unique_partner_count"] == 0:
            score += 10
        elif atom["unique_partner_count"] <= 2:
            score += 7
        elif atom["unique_partner_count"] <= 5:
            score += 3

        if atom_role == ROLE_DIRECT:
            score += 10
        elif atom_role == ROLE_CONDITIONAL:
            score += 5
        elif atom_role == ROLE_CONTEXT:
            score += 1

        if outward == 1:
            score += 8 if (atom["outward_score"] or 0.0) >= 0.25 else 5
        elif outward is None:
            score += 1

        if outward == 1 and clear_exit == 1:
            score += 7

        # Gate-aware caps.
        if not core:
            score = min(score, 39.0)
        elif atom_role in {ROLE_UNCLASSIFIED, ROLE_CONTEXT}:
            score = min(score, 59.0)
        elif atom_role == ROLE_CONDITIONAL:
            # Conditional sites are useful but never labeled High by this stage.
            score = min(score, 79.0)
        elif atom_role == ROLE_DIRECT and (outward != 1 or clear_exit != 1):
            score = min(score, 79.0)

        score = round(p.clamp(score), 2)
        tier = (
            "High attachment-site priority" if high_priority and score >= 80 else
            "Moderate attachment-site priority" if core and chemical_support and score >= 60 else
            "Exploratory attachment-site priority" if score >= 40 else
            "Low attachment-site priority"
        )

        rows.append({
            **atom,
            "candidate_core": int(core),
            "chemical_context": broad_chem,
            "chemical_support": int(chemical_support),
            **chem,
            "conditional_clear_exit": int(conditional_clear),
            "high_priority": int(high_priority),
            "attachment_priority_score": score,
            "attachment_priority_tier": tier,
        })
    return rows


def _write_report(database: str) -> None:
    c.dirs()
    with c.dbconn(database) as db:
        summary = dict(db.execute(
            """SELECT count(*) instances,
                      sum(mapped_atom_count) mapped_atoms,
                      sum(exposed_mapped_atom_count) exposed_mapped_atoms,
                      sum(candidate_attachment_atom_count) candidate_atoms,
                      sum(outward_supported_candidate_count) outward_supported_candidates,
                      sum(clear_exit_supported_candidate_count) clear_exit_supported_candidates,
                      sum(chemically_supported_candidate_count) chemically_supported_candidates,
                      sum(direct_attachment_candidate_count) direct_attachment_candidates,
                      sum(conditional_substitution_candidate_count) conditional_substitution_candidates,
                      sum(conditional_clear_exit_candidate_count) conditional_clear_exit_candidates,
                      sum(high_priority_attachment_atom_count) high_priority_atoms
               FROM protacability_attachment_site_summary WHERE method_version=?""",
            (VERSION,),
        ).fetchone())
        tiers = db.execute(
            "SELECT attachment_priority_tier,count(*) n FROM protacability_attachment_sites WHERE method_version=? GROUP BY 1 ORDER BY n DESC",
            (VERSION,),
        ).fetchall()
        roles = db.execute(
            "SELECT atom_chemical_role,count(*) n FROM protacability_attachment_sites WHERE method_version=? GROUP BY 1 ORDER BY n DESC",
            (VERSION,),
        ).fetchall()
    lines = [
        "# Stage 13 attachment-site report",
        "",
        f"* Method: {VERSION}",
        f"* Stage 12 evidence helper: {EXPECTED_PROTACABILITY_VERSION}",
        f"* Mapping evidence: {EXPECTED_MAPPING_VERSION}",
        f"* Functional-group evidence: {EXPECTED_FUNCTIONAL_GROUP_VERSION}",
        f"* Geometry evidence: {EXPECTED_GEOMETRY_VERSION}",
        "* Atom-specific chemistry: curated source-SMILES reaction-center rules; Stage-11 broad functional-group context retained separately",
    ]
    lines += [f"* {k.replace('_',' ')}: {v or 0}" for k, v in summary.items()]
    lines += ["", "## Atom chemical roles"] + [f"* {r['atom_chemical_role']}: {r['n']}" for r in roles]
    lines += ["", "## Atom priority tiers"] + [f"* {r['attachment_priority_tier']}: {r['n']}" for r in tiers]
    lines += [
        "",
        "High attachment-site priority requires mapping, solvent exposure, absence of a strong non-water environment contact, an atom-specific direct attachment handle, outward orientation, and a locally clear forward corridor.",
        "Conditional substitution sites can reach Moderate priority but are never labeled High solely from structural geometry.",
        "Functional-group context alone is not treated as atom-specific chemical support.",
        "Candidate and priority labels are medicinal-chemistry/structural hypothesis-generation outputs, not experimentally calibrated modification-tolerance predictions.",
    ]
    (c.ROOT / "outputs" / "ATTACHMENT_SITE_STAGE_REPORT.md").write_text("\n".join(lines) + "\n")


def run(database: str, limit=None, pdb_id=None, instance_id=None, resume=False, progress_every=100):
    _validate_dependencies()
    ensure_schema(database)
    c.dirs()
    RDLogger.DisableLog("rdApp.*")
    with c.dbconn(database) as db:
        q = """SELECT i.ligand_instance_id,s.entry_id,i.label_comp_id,i.auth_asym_id,i.auth_seq_id,i.insertion_code_normalized
               FROM ligand_instances i JOIN structures s ON s.structure_id=i.structure_id
               WHERE i.curation_status='included'"""
        args = []
        if pdb_id:
            q += " AND UPPER(s.entry_id)=UPPER(?)"
            args.append(pdb_id)
        if instance_id:
            q += " AND i.ligand_instance_id=?"
            args.append(instance_id)
        if resume:
            q += " AND NOT EXISTS (SELECT 1 FROM protacability_attachment_site_summary x WHERE x.ligand_instance_id=i.ligand_instance_id AND x.method_version=? AND x.status='complete')"
            args.append(VERSION)
        instances = db.execute(q + " ORDER BY i.ligand_instance_id", args).fetchall()
        if limit:
            instances = instances[:limit]
        ids = [r["ligand_instance_id"] for r in instances]
        rid = c.run_start(db, "attachment_sites", {
            "method": VERSION,
            "limit": limit,
            "pdb_id": pdb_id,
            "ligand_instance_id": instance_id,
            "resume": resume,
            "protacability_helper_version": EXPECTED_PROTACABILITY_VERSION,
            "mapping_method_version": EXPECTED_MAPPING_VERSION,
            "functional_group_method_version": EXPECTED_FUNCTIONAL_GROUP_VERSION,
            "geometry_method_version": EXPECTED_GEOMETRY_VERSION,
            "chemical_role_model": "atom-specific-source-smiles-reaction-centers-v1",
            "high_priority_requires_direct_handle": True,
        })
        if ids and not resume:
            marks = ",".join("?" for _ in ids)
            db.execute(
                f"DELETE FROM protacability_attachment_sites WHERE method_version=? AND ligand_instance_id IN ({marks})",
                [VERSION, *ids],
            )
            db.execute(
                f"DELETE FROM protacability_attachment_site_summary WHERE method_version=? AND ligand_instance_id IN ({marks})",
                [VERSION, *ids],
            )

        success = failures = 0
        for n, inst in enumerate(instances, 1):
            iid = inst["ligand_instance_id"]
            try:
                atom_rows = _attachment_rows_for_instance(db, iid)
                if not atom_rows:
                    raise ValueError("no selected ligand atoms")
                for a in atom_rows:
                    db.execute(
                        """INSERT OR REPLACE INTO protacability_attachment_sites(
                             run_id,ligand_instance_id,ligand_instance_atom_id,pdb_code,ligand_resname,ligand_chain,
                             ligand_residue_id,ligand_insertion_code,atom_site_id,exact_atom,element,smiles_atom_indices,
                             mapped,sasa_area_a2,solvent_exposed,meaningful_contact_count,strong_contact_count,
                             unique_partner_count,contact_labels,functional_groups,chemical_context,chemical_support,
                             atom_chemical_role,direct_attachment_support,conditional_substitution_support,chemical_rule_labels,chemical_rationale,
                             nearest_protein_distance_a,outward_score,points_away_from_pocket,forward_clearance_a,forward_obstruction_count,
                             exit_vector_clear,local_corridor_clear,forward_clearance_reaches_cap,candidate_attachment_atom,high_priority_attachment_atom,
                             attachment_priority_score,attachment_priority_tier,method_version)
                           VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)""",
                        (
                            rid, iid, a["ligand_instance_atom_id"], inst["entry_id"], inst["label_comp_id"], inst["auth_asym_id"],
                            inst["auth_seq_id"], inst["insertion_code_normalized"], a["atom_site_id"], a["exact_atom"], a["element"],
                            ";".join(str(x) for x in sorted(a["smiles_atom_indices"])), int(a["mapped"]), a["sasa_area"], int(a["exposed"]),
                            a["meaningful_contact_count"], a["strong_contact_count"], a["unique_partner_count"],
                            ";".join(sorted(set(a["labels"]))), ";".join(sorted(a["functional_groups"])), a["chemical_context"],
                            a["chemical_support"], a["atom_chemical_role"], a["direct_attachment_support"],
                            a["conditional_substitution_support"], a["chemical_rule_labels"], a["chemical_rationale"],
                            a["nearest_protein_distance_a"], a["outward_score"], a["points_away"], a["forward_clearance_a"],
                            a["forward_obstruction_count"], a["exit_vector_clear"], a["local_corridor_clear"],
                            a["forward_clearance_reaches_cap"], a["candidate_core"], a["high_priority"],
                            a["attachment_priority_score"], a["attachment_priority_tier"], VERSION,
                        ),
                    )

                mapped_count = sum(a["mapped"] for a in atom_rows)
                exposed_mapped = sum(a["mapped"] and a["exposed"] for a in atom_rows)
                candidates = [a for a in atom_rows if a["candidate_core"]]
                outward = sum(a["candidate_core"] and a["points_away"] == 1 for a in atom_rows)
                clear_exit = sum(
                    a["candidate_core"] and a["points_away"] == 1 and a["local_corridor_clear"] == 1
                    for a in atom_rows
                )
                chemical = sum(a["candidate_core"] and a["chemical_support"] for a in atom_rows)
                direct = sum(a["candidate_core"] and a["direct_attachment_support"] for a in atom_rows)
                conditional = sum(a["candidate_core"] and a["conditional_substitution_support"] for a in atom_rows)
                conditional_clear = sum(a["conditional_clear_exit"] for a in atom_rows)
                high = sum(a["high_priority"] for a in atom_rows)
                high_direct = sum(a["high_priority"] and a["direct_attachment_support"] for a in atom_rows)
                top = max(
                    atom_rows,
                    key=lambda x: (
                        x["attachment_priority_score"],
                        x["high_priority"],
                        x["direct_attachment_support"],
                        x["conditional_substitution_support"],
                        x["candidate_core"],
                        -int(x["ligand_instance_atom_id"]),
                    ),
                )
                db.execute(
                    """INSERT OR REPLACE INTO protacability_attachment_site_summary(
                         run_id,ligand_instance_id,pdb_code,ligand_resname,mapped_atom_count,exposed_mapped_atom_count,
                         candidate_attachment_atom_count,outward_supported_candidate_count,clear_exit_supported_candidate_count,
                         chemically_supported_candidate_count,direct_attachment_candidate_count,conditional_substitution_candidate_count,
                         conditional_clear_exit_candidate_count,high_priority_attachment_atom_count,high_priority_direct_attachment_atom_count,
                         top_attachment_site_score,top_attachment_atom_site_id,top_attachment_exact_atom,method_version,status)
                       VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)""",
                    (
                        rid, iid, inst["entry_id"], inst["label_comp_id"], mapped_count, exposed_mapped, len(candidates), outward,
                        clear_exit, chemical, direct, conditional, conditional_clear, high, high_direct,
                        top["attachment_priority_score"], top["atom_site_id"], top["exact_atom"], VERSION, "complete",
                    ),
                )
                success += 1
            except Exception as exc:
                failures += 1
                c.fail(
                    db, rid, "attachment_sites", f"{type(exc).__name__}: {exc}",
                    instance_id=iid, code="attachment_site_exception",
                )
            if n % max(1, progress_every) == 0 or n == len(instances):
                db.commit()
                print(
                    f"attachment-site progress: {n}/{len(instances)} success={success} failures={failures}",
                    flush=True,
                )
        c.run_end(
            db, rid, "completed" if failures == 0 else "partial",
            len(instances), success, 0, failures,
        )
    _write_report(database)
    return {"run_id": rid, "processed": len(instances), "success": success, "failures": failures}


def main():
    ap = argparse.ArgumentParser(description="Build atom-level attachment-site evidence with atom-specific chemical tractability.")
    ap.add_argument("--database", default=str(c.ROOT / "viral_data_cif_v2.db"))
    ap.add_argument("--limit", type=int)
    ap.add_argument("--pdb-id")
    ap.add_argument("--ligand-instance-id", type=int)
    ap.add_argument("--resume", action="store_true")
    ap.add_argument("--progress-every", type=int, default=100)
    a = ap.parse_args()
    print(json.dumps(run(a.database, a.limit, a.pdb_id, a.ligand_instance_id, a.resume, a.progress_every), indent=2))


if __name__ == "__main__":
    main()
