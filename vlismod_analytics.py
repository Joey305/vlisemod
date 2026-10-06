"""V-LiSEMOD client-side analytics plumbing with privacy-safe server enrichment."""
from __future__ import annotations

import hashlib
import os
import secrets
from datetime import datetime, timezone
from functools import wraps
from urllib.parse import urlparse

import requests
from flask import Blueprint, current_app, jsonify, redirect, render_template, request, session, url_for


bp = Blueprint("vlismod_analytics_ui", __name__)
COOKIE_VISITOR, COOKIE_SESSION = "vlismod_vid", "vlismod_sid"
EVENTS = {"workflow_started", "upload_started", "analysis_submitted", "analysis_completed", "analysis_failed", "results_viewed", "export_generated", "companion_handoff"}


def _ids():
    return request.cookies.get(COOKIE_VISITOR) or secrets.token_urlsafe(24), request.cookies.get(COOKIE_SESSION) or secrets.token_urlsafe(24)


def _device_class() -> str:
    ua = request.user_agent.string.lower()
    if "ipad" in ua or "tablet" in ua: return "tablet"
    return "mobile" if any(x in ua for x in ("mobile", "iphone", "android")) else "desktop"


def _referrer() -> str:
    host = urlparse(request.referrer or "").hostname
    return host[:100] if host else "direct"


def _geoip() -> dict:
    if os.environ.get("VLISMOD_USAGE_GEOIP", "0") != "1": return {}
    # The address exists only in this outbound lookup and is never persisted or logged.
    try:
        # Heroku terminates TLS before Flask, so remote_addr is its private router
        # address.  The left-most forwarded value is the original browser address.
        client_ip = (request.headers.get("X-Forwarded-For", "").split(",", 1)[0].strip() or request.remote_addr)
        reply = requests.get(f"https://ipwho.is/{client_ip}", timeout=1.5).json()
        if reply.get("success") is False: return {}
        return {"country_code": str(reply.get("country_code") or "")[:3], "country_name": str(reply.get("country") or "")[:80], "latitude": reply.get("latitude"), "longitude": reply.get("longitude")}
    except (requests.RequestException, ValueError):
        return {}


def _emit(event_type: str, *, path: str | None = None, feature: str = "", failure_stage: str = "", event_id: str | None = None, visitor_id: str | None = None, session_id: str | None = None) -> bool:
    if not current_app.config.get("RANDY_ANALYTICS_ENABLED", True): return False
    visitor_id, session_id = visitor_id or _ids()[0], session_id or _ids()[1]
    event_id = event_id or secrets.token_urlsafe(24)
    payload = {"event_id": event_id, "event_type": event_type, "occurred_at": datetime.now(timezone.utc).isoformat(), "visitor_id": visitor_id, "session_id": session_id, "path": path or request.path, "referrer": _referrer(), "device_class": _device_class(), "feature": feature[:80], "failure_stage": failure_stage[:40], **_geoip()}
    try:
        current_app.config["RANDY_POST"]("analytics/events", json=payload, max_bytes=1024 * 1024)
        return True
    except Exception as exc:  # Usage collection never interrupts scientific work.
        current_app.logger.info("VLiSEMOD analytics receiver unavailable: %s", exc)
        return False


def track_page_view(response):
    if request.method != "GET" or request.path.startswith(("/static/", "/admin", "/analytics/", "/api/", "/get_", "/health")): return response
    if response.mimetype != "text/html" or response.status_code >= 400: return response
    visitor_id, session_id = _ids()
    _emit("page_view", visitor_id=visitor_id, session_id=session_id)
    if not request.cookies.get(COOKIE_VISITOR): response.set_cookie(COOKIE_VISITOR, visitor_id, max_age=60 * 60 * 24 * 400, secure=request.is_secure, httponly=True, samesite="Lax")
    if not request.cookies.get(COOKIE_SESSION): response.set_cookie(COOKIE_SESSION, session_id, secure=request.is_secure, httponly=True, samesite="Lax")
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
    if event_type not in EVENTS: return jsonify({"ok": False}), 400
    event_id = str(payload.get("event_id") or "")[:80]
    # Only a fixed, non-sensitive feature label and safe failure stage are accepted.
    ok = _emit(event_type, feature=str(payload.get("feature") or "")[:80], failure_stage=str(payload.get("failure_stage") or "")[:40], event_id=event_id or None)
    return jsonify({"ok": ok})
