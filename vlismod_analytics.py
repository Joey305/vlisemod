"""V-LiSEMOD client-side analytics plumbing with privacy-safe server enrichment."""
from __future__ import annotations

import os
import re
import secrets
import threading
from datetime import datetime, timezone
from functools import wraps
from urllib.parse import urlparse

import requests
from flask import Blueprint, current_app, jsonify, redirect, render_template, request, session, url_for


bp = Blueprint("vlismod_analytics_ui", __name__)
COOKIE_VISITOR, COOKIE_SESSION = "vlismod_vid", "vlismod_sid"
EVENTS = {"workflow_started", "upload_started", "analysis_submitted", "analysis_completed", "analysis_failed", "results_viewed", "export_generated", "companion_handoff"}
FEATURES = {
    "analysis_builder", "protein_query", "ligand_query", "ligand_comparison", "protacability",
    "ligand_interactions", "protacability_search", "pymol_session", "ligand_images",
    "protein_query_results", "ligand_interaction_results", "ligand_comparison_results",
    "protacability_target_detail", "protacability_structure_detail", "protacability_ligand_detail",
    "protein_query_export", "protacability_evidence_export", "builder_from_landing",
    "builder_from_ligand_query", "builder_from_ligand_comparison", "builder_from_protacability",
    "builder_from_viral_protac_design",
}
FAILURE_STAGES = {"selection", "input_validation", "interaction_analysis", "comparison", "protacability_search", "image_generation", "pymol_generation", "export", "handoff", "unknown"}
HANDOFF_ID_RE = re.compile(r"^[0-9a-f]{8}-[0-9a-f]{4}-4[0-9a-f]{3}-[89ab][0-9a-f]{3}-[0-9a-f]{12}$", re.I)
SERVER_OUTCOMES = {
    "/generate_pymol_session": ("pymol_session", "pymol_generation", True),
    "/generate_ligand_images": ("ligand_images", "image_generation", True),
    "/generate_charts": ("ligand_interactions", "interaction_analysis", False),
    "/compare_ligand_interactions": ("ligand_comparison", "comparison", False),
    "/get_pdbs_for_virus_protein": ("protein_query", "selection", False),
    "/api/protacability/search": ("protacability_search", "protacability_search", False),
    "/api/protacability/export": ("protacability_evidence_export", "export", False),
    "/export_data_to_excel": ("protein_query_export", "export", False),
}


def _ids():
    return request.cookies.get(COOKIE_VISITOR) or secrets.token_urlsafe(24), request.cookies.get(COOKIE_SESSION) or secrets.token_urlsafe(24)


def _device_class() -> str:
    ua = request.user_agent.string.lower()
    if "ipad" in ua or "tablet" in ua: return "tablet"
    return "mobile" if any(x in ua for x in ("mobile", "iphone", "android")) else "desktop"


def _referrer() -> str:
    host = urlparse(request.referrer or "").hostname
    return host[:100] if host else "direct"


def _geoip_for_ip(client_ip: str) -> dict:
    if os.environ.get("VLISMOD_USAGE_GEOIP", "0") != "1": return {}
    # This transient address is used for a country lookup only, never stored or logged.
    try:
        reply = requests.get(f"https://ipwho.is/{client_ip}", timeout=1.5).json()
        if reply.get("success") is False: return {}
        return {"country_code": str(reply.get("country_code") or "")[:3], "country_name": str(reply.get("country") or "")[:80], "latitude": reply.get("latitude"), "longitude": reply.get("longitude")}
    except (requests.RequestException, ValueError):
        return {}


def _geoip() -> dict:
    return _geoip_for_ip(request.headers.get("X-Forwarded-For", "").split(",", 1)[0].strip() or request.remote_addr)


def _deliver(payload, client_ip, sender, logger):
    try:
        payload.update(_geoip_for_ip(client_ip))
        sender("analytics/events", json=payload, max_bytes=1024 * 1024)
    except Exception as exc:  # Analytics cannot affect the user workflow.
        logger.info("VLiSEMOD analytics receiver unavailable: %s", type(exc).__name__)


def _emit(event_type: str, *, path: str | None = None, feature: str = "", failure_stage: str = "", event_id: str | None = None, visitor_id: str | None = None, session_id: str | None = None, handoff_id: str = "") -> bool:
    if event_type not in EVENTS | {"page_view"} or (feature and feature not in FEATURES):
        return False
    if failure_stage and failure_stage not in FAILURE_STAGES:
        failure_stage = "unknown"
    if event_type == "companion_handoff" and not HANDOFF_ID_RE.fullmatch(handoff_id):
        return False
    if handoff_id and (event_type != "companion_handoff" or not HANDOFF_ID_RE.fullmatch(handoff_id)):
        return False
    if not current_app.config.get("RANDY_ANALYTICS_ENABLED", True): return False
    visitor_id, session_id = visitor_id or _ids()[0], session_id or _ids()[1]
    event_id = event_id or secrets.token_urlsafe(24)
    payload = {"event_id": event_id, "event_type": event_type, "occurred_at": datetime.now(timezone.utc).isoformat(), "visitor_id": visitor_id, "session_id": session_id, "path": path or request.path, "referrer": _referrer(), "device_class": _device_class(), "feature": feature[:80], "failure_stage": failure_stage[:40], "handoff_id": handoff_id}
    client_ip = request.headers.get("X-Forwarded-For", "").split(",", 1)[0].strip() or request.remote_addr
    sender, logger = current_app.config["RANDY_POST"], current_app.logger
    if current_app.config.get("RANDY_ANALYTICS_SYNC", False):
        _deliver(payload, client_ip, sender, logger)
    else:
        threading.Thread(target=_deliver, args=(payload, client_ip, sender, logger), daemon=True, name="vlismod-analytics").start()
    return True


def track_page_view(response):
    if request.method != "GET" or request.path.startswith(("/static/", "/admin", "/analytics/", "/api/", "/get_", "/health")): return response
    if response.mimetype != "text/html" or response.status_code >= 400: return response
    visitor_id, session_id = _ids()
    _emit("page_view", visitor_id=visitor_id, session_id=session_id)
    if not request.cookies.get(COOKIE_VISITOR): response.set_cookie(COOKIE_VISITOR, visitor_id, max_age=60 * 60 * 24 * 400, secure=request.is_secure, httponly=True, samesite="Lax")
    if not request.cookies.get(COOKIE_SESSION): response.set_cookie(COOKIE_SESSION, session_id, secure=request.is_secure, httponly=True, samesite="Lax")
    return response


def track_server_outcome(response):
    """Record only safe route outcomes; request payloads never enter analytics."""
    outcome = SERVER_OUTCOMES.get(request.path)
    if not outcome or request.method not in {"POST", "GET"}:
        return response
    feature, failure_stage, is_output = outcome
    if 200 <= response.status_code < 400:
        if is_output:
            _emit("analysis_completed", feature=feature)
        elif feature.endswith("_export"):
            _emit("export_generated", feature=feature)
        else:
            _emit("analysis_completed", feature=feature)
            result_feature = {
                "protein_query": "protein_query_results", "ligand_interactions": "ligand_interaction_results",
                "ligand_comparison": "ligand_comparison_results", "protacability_search": "protacability_target_detail",
            }.get(feature)
            if result_feature:
                _emit("results_viewed", feature=result_feature)
    elif response.status_code >= 400:
        _emit("analysis_failed", feature=feature, failure_stage=failure_stage)
    return response


def _admin_required(view):
    @wraps(view)
    def wrapped(*args, **kwargs):
        if not session.get("vlismod_admin"):
            return redirect(url_for("vlismod_analytics_ui.admin_login", next=request.path))
        return view(*args, **kwargs)
    return wrapped


@bp.route("/admin/login", methods=["GET", "POST"])
def admin_login():
    error = None
    if request.method == "POST":
        import hmac
        expected_email, expected_password = os.environ.get("ADMIN_EMAIL", ""), os.environ.get("ADMIN_PASSWORD", "")
        if expected_email and expected_password and hmac.compare_digest(request.form.get("email", ""), expected_email) and hmac.compare_digest(request.form.get("password", ""), expected_password):
            session.clear(); session["vlismod_admin"] = True
            return redirect(request.args.get("next") or url_for("vlismod_analytics_ui.admin_analytics"))
        error = "Invalid administrator credentials."
    return render_template("admin_login.html", error=error)


@bp.post("/admin/logout")
def admin_logout():
    session.clear()
    return redirect(url_for("vlismod_analytics_ui.admin_login"))


@bp.get("/admin/analytics")
@_admin_required
def admin_analytics():
    period = request.args.get("period", "30d")
    if period not in {"7d", "30d", "90d", "1y", "all"}: period = "30d"
    try:
        data = current_app.config["RANDY_GET"]("analytics/rollup", params={"period": period})
        error = None
    except Exception as exc:
        current_app.logger.warning("Analytics rollup unavailable: %s", exc)
        data, error = None, "The analytics receiver is temporarily unavailable. Existing data was not changed."
    return render_template("admin_analytics.html", data=data, period=period, error=error)


@bp.post("/analytics/event")
def browser_event():
    payload = request.get_json(silent=True) or {}
    event_type = str(payload.get("event_type") or "")
    feature = str(payload.get("feature") or "")
    if event_type not in EVENTS or feature not in FEATURES: return jsonify({"ok": False}), 400
    event_id = str(payload.get("event_id") or "")[:80]
    handoff_id = str(payload.get("handoff_id") or "")[:64]
    if event_type == "companion_handoff" and not HANDOFF_ID_RE.fullmatch(handoff_id):
        return jsonify({"ok": False}), 400
    # Only a fixed, non-sensitive feature label and safe failure stage are accepted.
    ok = _emit(event_type, feature=feature, failure_stage=str(payload.get("failure_stage") or "")[:40], event_id=event_id or None, handoff_id=handoff_id)
    return jsonify({"ok": ok})
