"""Privacy-preserving, V-LiSEMOD-only analytics receiver routes.

This module deliberately stores no request IP, molecular input, search terms, or
free-form user content.  The V-LiSEMOD web app enriches safe events and sends
them through the existing bearer-token channel.
"""
from __future__ import annotations

import json
import os
import sqlite3
import re
from collections import Counter
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any

from flask import Blueprint, jsonify, request


bp = Blueprint("vlismod_analytics", __name__)
_ALLOWED_EVENTS = {"page_view", "workflow_started", "upload_started", "analysis_submitted", "analysis_completed", "analysis_failed", "results_viewed", "export_generated", "companion_handoff"}
_ALLOWED_FEATURES = {"analysis_builder", "protein_query", "ligand_query", "ligand_comparison", "protacability", "ligand_interactions", "protacability_search", "pymol_session", "ligand_images", "protein_query_results", "ligand_interaction_results", "ligand_comparison_results", "protacability_target_detail", "protacability_structure_detail", "protacability_ligand_detail", "protein_query_export", "protacability_evidence_export", "builder_from_landing", "builder_from_ligand_query", "builder_from_ligand_comparison", "builder_from_protacability", "builder_from_viral_protac_design"}
_ALLOWED_STAGES = {"selection", "input_validation", "interaction_analysis", "comparison", "protacability_search", "image_generation", "pymol_generation", "export", "handoff", "unknown"}
_HANDOFF_ID_RE = re.compile(r"^[0-9a-f]{8}-[0-9a-f]{4}-4[0-9a-f]{3}-[89ab][0-9a-f]{3}-[0-9a-f]{12}$", re.I)


def _token() -> str:
    return os.environ.get("VLISMOD_API_TOKEN", "").strip() or os.environ.get("RANDY_BACKUP_TOKEN", "").strip()


def _authorized() -> bool:
    token = _token()
    return bool(token) and request.headers.get("Authorization", "") == f"Bearer {token}"


def _db_path() -> Path:
    root = Path(os.environ.get("PROTAC_BACKUP_DIR", str(Path(__file__).resolve().parents[1] / "data"))).expanduser()
    root.mkdir(parents=True, exist_ok=True)
    return root / "vlismod_analytics.sqlite3"


def _connect() -> sqlite3.Connection:
    conn = sqlite3.connect(_db_path())
    conn.row_factory = sqlite3.Row
    conn.execute("""CREATE TABLE IF NOT EXISTS vlismod_events (
        event_id TEXT PRIMARY KEY, occurred_at TEXT NOT NULL, event_type TEXT NOT NULL,
        visitor_id TEXT NOT NULL, session_id TEXT NOT NULL, path TEXT NOT NULL,
        referrer TEXT NOT NULL, device_class TEXT NOT NULL, country_code TEXT,
        country_name TEXT, latitude REAL, longitude REAL, feature TEXT, failure_stage TEXT, handoff_id TEXT
    )""")
    columns = {row[1] for row in conn.execute("PRAGMA table_info(vlismod_events)")}
    if "handoff_id" not in columns:
        conn.execute("ALTER TABLE vlismod_events ADD COLUMN handoff_id TEXT")
    conn.execute("CREATE INDEX IF NOT EXISTS idx_vlismod_events_time ON vlismod_events(occurred_at)")
    conn.execute("CREATE INDEX IF NOT EXISTS idx_vlismod_events_type ON vlismod_events(event_type)")
    return conn


def _safe(value: Any, limit: int = 120) -> str:
    return str(value or "").strip()[:limit]


@bp.post("/analytics/events")
def ingest_event():
    if not _authorized():
        return jsonify({"ok": False, "error": "Unauthorized."}), 401
    data = request.get_json(silent=True) or {}
    event_type = _safe(data.get("event_type"), 40)
    if event_type not in _ALLOWED_EVENTS:
        return jsonify({"ok": False, "error": "Unsupported analytics event."}), 400
    feature = _safe(data.get("feature"), 80)
    if event_type != "page_view" and feature not in _ALLOWED_FEATURES:
        return jsonify({"ok": False, "error": "Unsupported analytics feature."}), 400
    handoff_id = _safe(data.get("handoff_id"), 64)
    if event_type == "companion_handoff" and not _HANDOFF_ID_RE.fullmatch(handoff_id):
        return jsonify({"ok": False, "error": "Invalid handoff identifier."}), 400
    if handoff_id and (event_type != "companion_handoff" or not _HANDOFF_ID_RE.fullmatch(handoff_id)):
        return jsonify({"ok": False, "error": "Invalid handoff identifier."}), 400
    event_id, visitor_id, session_id = (_safe(data.get(k), 80) for k in ("event_id", "visitor_id", "session_id"))
    if not event_id or not visitor_id or not session_id:
        return jsonify({"ok": False, "error": "Missing anonymous event identifier."}), 400
    path = _safe(data.get("path"), 240)
    if not path.startswith("/") or "?" in path or "#" in path:
        return jsonify({"ok": False, "error": "Invalid normalized path."}), 400
    occurred = _safe(data.get("occurred_at"), 40) or datetime.now(timezone.utc).isoformat()
    stage = _safe(data.get("failure_stage"), 40)
    if stage not in _ALLOWED_STAGES:
        stage = ""
    row = (event_id, occurred, event_type, visitor_id, session_id, path, _safe(data.get("referrer"), 100) or "direct", _safe(data.get("device_class"), 20) or "desktop", _safe(data.get("country_code"), 3), _safe(data.get("country_name"), 80), data.get("latitude"), data.get("longitude"), feature, stage, handoff_id)
    with _connect() as conn:
        before = conn.total_changes
        conn.execute("""INSERT OR IGNORE INTO vlismod_events
            (event_id, occurred_at, event_type, visitor_id, session_id, path, referrer, device_class,
             country_code, country_name, latitude, longitude, feature, failure_stage, handoff_id)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""", row)
        inserted = conn.total_changes > before
    return jsonify({"ok": True, "inserted": inserted})


def _start_date(period: str) -> str | None:
    days = {"7d": 7, "30d": 30, "90d": 90, "1y": 365}.get(period)
    return (datetime.now(timezone.utc) - timedelta(days=days)).isoformat() if days else None


@bp.get("/analytics/rollup")
def rollup():
    if not _authorized():
        return jsonify({"ok": False, "error": "Unauthorized."}), 401
    start = _start_date(request.args.get("period", "30d"))
    clause, params = ("WHERE occurred_at >= ?", (start,)) if start else ("", ())
    with _connect() as conn:
        rows = conn.execute(f"SELECT * FROM vlismod_events {clause} ORDER BY occurred_at DESC", params).fetchall()
    events = [dict(r) for r in rows]
    page_views = [e for e in events if e["event_type"] == "page_view"]
    visitors, sessions = {e["visitor_id"] for e in events}, {e["session_id"] for e in events}
    daily = Counter(e["occurred_at"][:10] for e in page_views)
    by_referrer, by_device = Counter(e["referrer"] for e in page_views), Counter(e["device_class"] for e in page_views)
    pages: dict[str, dict[str, Any]] = {}
    for e in page_views:
        item = pages.setdefault(e["path"], {"path": e["path"], "views": 0, "visitors": set()})
        item["views"] += 1; item["visitors"].add(e["visitor_id"])
    countries: dict[str, dict[str, Any]] = {}
    for e in page_views:
        if not e["country_code"]: continue
        item = countries.setdefault(e["country_code"], {"code": e["country_code"], "name": e["country_name"] or e["country_code"], "count": 0, "latitude": e["latitude"], "longitude": e["longitude"]})
        item["count"] += 1
    funnel_types = [("Landing page", "page_view"), ("Workflow started", "workflow_started"), ("Analysis submitted", "analysis_submitted"), ("Analysis completed", "analysis_completed"), ("Results viewed", "results_viewed"), ("Export or Builder handoff", "terminal")]
    funnel = []
    previous = None
    for label, typ in funnel_types:
        count = len({e["session_id"] for e in events if e["event_type"] in {"export_generated", "companion_handoff"}}) if typ == "terminal" else len({e["session_id"] for e in events if e["event_type"] == typ})
        funnel.append({"label": label, "count": count, "conversion": round((count / previous * 100), 1) if previous else None})
        previous = count
    submitted = sum(e["event_type"] == "analysis_submitted" for e in events); completed = sum(e["event_type"] == "analysis_completed" for e in events); failed = sum(e["event_type"] == "analysis_failed" for e in events)
    workflows = []
    workflow_specs = [("Protein Query", "protein_query", "protein_query_export", ()), ("Ligand Interactions", "ligand_interactions", "", ("builder_from_ligand_query",)), ("Ligand Comparison", "ligand_comparison", "", ("builder_from_ligand_comparison",)), ("PROTACability", "protacability_search", "protacability_evidence_export", ("builder_from_protacability",)), ("PyMOL Session", "pymol_session", "", ()), ("Ligand Images", "ligand_images", "", ())]
    for label, feature, export_feature, handoff_origins in workflow_specs:
        relevant = [e for e in events if e["feature"] == feature]
        starts = sum(e["event_type"] in {"workflow_started", "analysis_submitted"} for e in relevant)
        done = sum(e["event_type"] == "analysis_completed" for e in relevant)
        failures = sum(e["event_type"] == "analysis_failed" for e in relevant)
        workflows.append({"label": label, "starts": starts, "completed": done, "failed": failures, "completion_rate": round(done / starts * 100, 1) if starts else None, "exports": sum(e["event_type"] == "export_generated" and e["feature"] == export_feature for e in events), "handoffs": sum(e["event_type"] == "companion_handoff" and e["feature"] in handoff_origins for e in events)})
    handoffs = Counter(e["feature"] for e in events if e["event_type"] == "companion_handoff")
    return jsonify({"ok": True, "period": request.args.get("period", "30d"), "legacy": {"available": False, "note": "No inferred historical visitor, session, referrer, or location data is displayed."}, "summary": {"unique_visitors": len(visitors), "sessions": len(sessions), "page_views": len(page_views), "workflow_starts": sum(e["event_type"] == "workflow_started" for e in events), "completed_analyses": completed, "exports": sum(e["event_type"] == "export_generated" for e in events), "builder_handoffs": sum(handoffs.values())}, "daily": [{"date": k, "views": v} for k, v in sorted(daily.items())], "referrers": by_referrer.most_common(10), "devices": by_device.most_common(), "pages": sorted(({**x, "visitors": len(x["visitors"])} for x in pages.values()), key=lambda x: x["views"], reverse=True)[:25], "countries": sorted(countries.values(), key=lambda x: x["count"], reverse=True), "funnel": funnel, "workflows": workflows, "handoffs": {"total": sum(handoffs.values()), "session_rate": round(len({e["session_id"] for e in events if e["event_type"] == "companion_handoff"}) / len(sessions) * 100, 1) if sessions else None, "by_origin": handoffs.most_common()}, "operations": {"submitted": submitted, "completed": completed, "failed": failed, "success_rate": round(completed / submitted * 100, 1) if submitted else None, "failures": [{"when": e["occurred_at"], "stage": e["failure_stage"] or "unknown", "feature": e["feature"] or "analysis"} for e in events if e["event_type"] == "analysis_failed"][:10], "top_features": Counter(e["feature"] or e["path"] for e in events if e["event_type"] != "page_view").most_common(10)}})
