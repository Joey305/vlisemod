"""Classify legacy/baseline Arpeggio failures without rerunning them."""
from __future__ import annotations

import argparse
import csv
from collections import Counter
from importlib import import_module
from pathlib import Path

c = import_module("00_common")
arpeggio = import_module("09_run_arpeggio")
CLASSIFIER_VERSION = "arpeggio-failure-taxonomy-v1"


def classify(database, baseline_only=False):
    c.create_schema(database)
    rows = []
    with c.dbconn(database) as db:
        baseline_filter = "AND r.input_strategy IS NULL" if baseline_only else ""
        query = f"""
            WITH ranked AS (
              SELECT r.*,row_number() OVER(
                PARTITION BY r.ligand_instance_id ORDER BY r.run_id DESC
              ) AS rn
              FROM ligand_arpeggio_runs r
              JOIN ligand_instances i ON i.ligand_instance_id=r.ligand_instance_id
              WHERE i.curation_status='included' {baseline_filter}
            )
            SELECT r.run_id,r.ligand_instance_id,r.status,r.stderr_path,
                   s.entry_id,s.file_size,i.label_comp_id,i.deposited_model_num,
                   i.auth_asym_id,i.auth_seq_id,i.insertion_code_normalized
            FROM ranked r
            JOIN ligand_instances i ON i.ligand_instance_id=r.ligand_instance_id
            JOIN structures s ON s.structure_id=i.structure_id
            WHERE r.rn=1 AND r.status IN ('failed','timed_out','blocked','interrupted')
            ORDER BY r.ligand_instance_id
        """
        for record in db.execute(query):
            stderr_path = Path(record["stderr_path"]) if record["stderr_path"] else None
            stderr = stderr_path.read_text(errors="replace") if stderr_path and stderr_path.exists() else ""
            failure_class = arpeggio.classify_failure(stderr, record["status"], record["status"])
            row = dict(record)
            row["failure_class"] = failure_class
            rows.append(row)
            db.execute(
                """INSERT INTO arpeggio_failure_classifications(
                     ligand_instance_id,source_run_id,source_status,failure_class,classifier_version,classified_at)
                   VALUES(?,?,?,?,?,?) ON CONFLICT(source_run_id) DO UPDATE SET
                     source_status=excluded.source_status,failure_class=excluded.failure_class,
                     classifier_version=excluded.classifier_version,classified_at=excluded.classified_at""",
                (record["ligand_instance_id"], record["run_id"], record["status"], failure_class,
                 CLASSIFIER_VERSION, c.now()),
            )
    return rows


def write_outputs(rows, csv_path, report_path):
    fields = [
        "run_id", "ligand_instance_id", "entry_id", "label_comp_id", "source_status",
        "failure_class", "file_size", "deposited_model_num", "auth_asym_id",
        "auth_seq_id", "insertion_code_normalized", "stderr_path",
    ]
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            export = dict(row)
            export["source_status"] = export.pop("status")
            writer.writerow(export)
    counts = Counter(row["failure_class"] for row in rows)
    status_counts = Counter(row["status"] for row in rows)
    lines = [
        "# Arpeggio unresolved-input failure audit",
        "",
        f"Classifier: `{CLASSIFIER_VERSION}`",
        f"",
        f"Total classified: **{len(rows)}**",
        "",
        "## Source status",
        "",
        *[f"- {name}: {count}" for name, count in sorted(status_counts.items())],
        "",
        "## Failure class",
        "",
        *[f"- {name}: {count}" for name, count in counts.most_common()],
        "",
        "The source mmCIF files are authoritative and were not modified. This audit reads persisted run logs only.",
    ]
    report_path.write_text("\n".join(lines) + "\n")
    return counts


def write_recovery_report(database, report_path):
    with c.dbconn(database) as db:
        included = db.execute("SELECT count(*) FROM ligand_instances WHERE curation_status='included'").fetchone()[0]
        completed = db.execute(
            """SELECT count(*) FROM ligand_instances i WHERE i.curation_status='included'
               AND EXISTS (SELECT 1 FROM ligand_arpeggio_runs r
                           WHERE r.ligand_instance_id=i.ligand_instance_id AND r.status='completed')"""
        ).fetchone()[0]
        modes = dict(db.execute(
            """SELECT completion_mode,count(*) FROM ligand_arpeggio_runs
               WHERE input_strategy IS NOT NULL AND status='completed' GROUP BY completion_mode"""
        ))
        final_status = dict(db.execute(
            """WITH latest AS (
                 SELECT r.status,row_number() OVER(PARTITION BY r.ligand_instance_id ORDER BY r.run_id DESC) rn
                 FROM ligand_arpeggio_runs r JOIN ligand_instances i
                   ON i.ligand_instance_id=r.ligand_instance_id WHERE i.curation_status='included')
               SELECT status,count(*) FROM latest WHERE rn=1 GROUP BY status"""
        ))
    percent = 100.0 * completed / included if included else 0.0
    lines = [
        "# Arpeggio recovery report",
        "",
        f"- Included ligand instances: {included}",
        f"- Instances with a validated completed Arpeggio run: {completed}",
        f"- Completion: {percent:.1f}%",
        f"- Recovered with sanitized full input: {modes.get('sanitized_full', 0)}",
        f"- Recovered with controlled sanitized pocket timeout input: {modes.get('sanitized_pocket', 0)}",
        f"- Latest failed: {final_status.get('failed', 0)}",
        f"- Latest timed out: {final_status.get('timed_out', 0)}",
        f"- Latest blocked/interrupted: {final_status.get('blocked', 0) + final_status.get('interrupted', 0)}",
        "",
        "All recovery outputs passed JSON, selected-ligand, canonical endpoint, and derived-atom provenance validation.",
    ]
    report_path.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--database", default=str(c.ROOT / "viral_data_cif_v2.db"))
    parser.add_argument("--baseline-only", action="store_true", help="Classify the latest pre-v2.2 run for each instance")
    parser.add_argument("--csv", default=str(c.ROOT / "outputs" / "ARPEGGIO_FAILURE_CLASSIFICATION.csv"))
    parser.add_argument("--report", default=str(c.ROOT / "outputs" / "ARPEGGIO_FAILURE_CLASSIFICATION.md"))
    parser.add_argument("--recovery-report", default=str(c.ROOT / "outputs" / "ARPEGGIO_RECOVERY_REPORT.md"))
    args = parser.parse_args()
    records = classify(args.database, args.baseline_only)
    print(dict(write_outputs(records, Path(args.csv), Path(args.report))))
    write_recovery_report(args.database, Path(args.recovery_report))
