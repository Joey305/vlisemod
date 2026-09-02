#!/usr/bin/env python3
"""Materialize atom-level attachment-site evidence and priorities.

The candidate flag preserves the historical ligand-centered logic: a deposited
heavy atom must be mapped to the ligand SMILES, solvent exposed, and free of a
strong Arpeggio contact.  Functional-group context and the stage-10 outward cue
are retained as additional prioritization evidence rather than being treated as
experimental proof of linker tolerance.
"""
from __future__ import annotations

import argparse
import importlib
import json

c = importlib.import_module("00_common")
p = importlib.import_module("12_build_protacability")
VERSION = "attachment-sites-cif-v2.1"


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
                nearest_protein_distance_a REAL,
                outward_score REAL,
                points_away_from_pocket INTEGER,
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
                chemically_supported_candidate_count INTEGER NOT NULL,
                high_priority_attachment_atom_count INTEGER NOT NULL,
                top_attachment_site_score REAL,
                top_attachment_atom_site_id TEXT,
                top_attachment_exact_atom TEXT,
                method_version TEXT NOT NULL,
                status TEXT NOT NULL,
                UNIQUE(ligand_instance_id,method_version)
            );
            """
        )


def _write_report(database):
    c.dirs()
    with c.dbconn(database) as db:
        summary = dict(db.execute(
            """SELECT count(*) instances,
                      sum(mapped_atom_count) mapped_atoms,
                      sum(exposed_mapped_atom_count) exposed_mapped_atoms,
                      sum(candidate_attachment_atom_count) candidate_atoms,
                      sum(outward_supported_candidate_count) outward_supported_candidates,
                      sum(chemically_supported_candidate_count) chemically_supported_candidates,
                      sum(high_priority_attachment_atom_count) high_priority_atoms
               FROM protacability_attachment_site_summary WHERE method_version=?""", (VERSION,)
        ).fetchone())
        tiers = db.execute(
            "SELECT attachment_priority_tier,count(*) n FROM protacability_attachment_sites WHERE method_version=? GROUP BY 1 ORDER BY n DESC",
            (VERSION,),
        ).fetchall()
    lines = ["# Stage 13 attachment-site report", "", f"* Method: {VERSION}"]
    lines += [f"* {k.replace('_',' ')}: {v or 0}" for k, v in summary.items()]
    lines += ["", "## Atom priority tiers"] + [f"* {r['attachment_priority_tier']}: {r['n']}" for r in tiers]
    lines += ["", "Candidate and priority labels are structural hypothesis-generation outputs, not experimental modification-tolerance predictions."]
    (c.ROOT / "outputs" / "ATTACHMENT_SITE_STAGE_REPORT.md").write_text("\n".join(lines) + "\n")


def run(database: str, limit=None, pdb_id=None, instance_id=None, resume=False, progress_every=100):
    ensure_schema(database); c.dirs()
    with c.dbconn(database) as db:
        q = """SELECT i.ligand_instance_id,s.entry_id,i.label_comp_id,i.auth_asym_id,i.auth_seq_id,i.insertion_code_normalized
               FROM ligand_instances i JOIN structures s ON s.structure_id=i.structure_id
               WHERE i.curation_status='included'"""
        args = []
        if pdb_id: q += " AND s.entry_id=?"; args.append(pdb_id)
        if instance_id: q += " AND i.ligand_instance_id=?"; args.append(instance_id)
        if resume:
            q += " AND NOT EXISTS (SELECT 1 FROM protacability_attachment_site_summary x WHERE x.ligand_instance_id=i.ligand_instance_id AND x.method_version=? AND x.status='complete')"
            args.append(VERSION)
        instances = db.execute(q + " ORDER BY i.ligand_instance_id", args).fetchall()
        if limit: instances = instances[:limit]
        ids = [r["ligand_instance_id"] for r in instances]
        rid = c.run_start(db, "attachment_sites", {"method": VERSION, "limit": limit, "pdb_id": pdb_id, "ligand_instance_id": instance_id, "resume": resume})
        if ids and not resume:
            marks = ",".join("?" for _ in ids)
            db.execute(f"DELETE FROM protacability_attachment_sites WHERE method_version=? AND ligand_instance_id IN ({marks})", [VERSION, *ids])
            db.execute(f"DELETE FROM protacability_attachment_site_summary WHERE method_version=? AND ligand_instance_id IN ({marks})", [VERSION, *ids])

        success = failures = 0
        for n, inst in enumerate(instances, 1):
            iid = inst["ligand_instance_id"]
            try:
                atom_rows = p.attachment_atom_rows_for_instance(db, iid)
                if not atom_rows:
                    raise ValueError("no selected ligand atoms")
                for a in atom_rows:
                    db.execute(
                        """INSERT OR REPLACE INTO protacability_attachment_sites(
                             run_id,ligand_instance_id,ligand_instance_atom_id,pdb_code,ligand_resname,ligand_chain,
                             ligand_residue_id,ligand_insertion_code,atom_site_id,exact_atom,element,smiles_atom_indices,
                             mapped,sasa_area_a2,solvent_exposed,meaningful_contact_count,strong_contact_count,
                             unique_partner_count,contact_labels,functional_groups,chemical_context,chemical_support,
                             nearest_protein_distance_a,outward_score,points_away_from_pocket,candidate_attachment_atom,
                             high_priority_attachment_atom,attachment_priority_score,attachment_priority_tier,method_version)
                           VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)""",
                        (rid, iid, a["ligand_instance_atom_id"], inst["entry_id"], inst["label_comp_id"], inst["auth_asym_id"],
                         inst["auth_seq_id"], inst["insertion_code_normalized"], a["atom_site_id"], a["exact_atom"], a["element"],
                         ";".join(str(x) for x in sorted(a["smiles_atom_indices"])), int(a["mapped"]), a["sasa_area"], int(a["exposed"]),
                         a["meaningful_contact_count"], a["strong_contact_count"], a["unique_partner_count"],
                         ";".join(sorted(set(a["labels"]))), ";".join(sorted(a["functional_groups"])), a["chemical_context"],
                         a["chemical_support"], a["nearest_protein_distance_a"], a["outward_score"], a["points_away"],
                         a["candidate_core"], a["high_priority"], a["attachment_priority_score"], a["attachment_priority_tier"], VERSION)
                    )
                mapped_count = sum(a["mapped"] for a in atom_rows)
                exposed_mapped = sum(a["mapped"] and a["exposed"] for a in atom_rows)
                candidates = [a for a in atom_rows if a["candidate_core"]]
                outward = sum(a["candidate_core"] and a["points_away"] == 1 for a in atom_rows)
                chemical = sum(a["candidate_core"] and a["chemical_support"] for a in atom_rows)
                high = sum(a["high_priority"] for a in atom_rows)
                top = max(atom_rows, key=lambda x: (x["attachment_priority_score"], x["candidate_core"], -int(x["ligand_instance_atom_id"])))
                db.execute(
                    """INSERT OR REPLACE INTO protacability_attachment_site_summary(
                         run_id,ligand_instance_id,pdb_code,ligand_resname,mapped_atom_count,exposed_mapped_atom_count,
                         candidate_attachment_atom_count,outward_supported_candidate_count,chemically_supported_candidate_count,
                         high_priority_attachment_atom_count,top_attachment_site_score,top_attachment_atom_site_id,
                         top_attachment_exact_atom,method_version,status)
                       VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)""",
                    (rid, iid, inst["entry_id"], inst["label_comp_id"], mapped_count, exposed_mapped, len(candidates), outward,
                     chemical, high, top["attachment_priority_score"], top["atom_site_id"], top["exact_atom"], VERSION, "complete")
                )
                success += 1
            except Exception as exc:
                failures += 1
                c.fail(db, rid, "attachment_sites", f"{type(exc).__name__}: {exc}", instance_id=iid, code="attachment_site_exception")
            if n % max(1, progress_every) == 0 or n == len(instances):
                db.commit(); print(f"attachment-site progress: {n}/{len(instances)} success={success} failures={failures}", flush=True)
        c.run_end(db, rid, "completed" if failures == 0 else "partial", len(instances), success, 0, failures)
    _write_report(database)
    return {"run_id": rid, "processed": len(instances), "success": success, "failures": failures}


def main():
    ap = argparse.ArgumentParser(description="Build atom-level attachment-site evidence and priority records.")
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
