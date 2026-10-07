import os

from flask import Flask

os.environ.setdefault("ADMIN_EMAIL", "admin@example.test")
os.environ.setdefault("ADMIN_PASSWORD", "correct-password")

from vlismod_analytics import bp, track_page_view, track_server_outcome


def make_app(events):
    app = Flask(__name__)
    app.secret_key = "test-secret"
    app.config.update(RANDY_ANALYTICS_ENABLED=True, RANDY_POST=lambda *args, **kwargs: events.append(kwargs["json"]) or {"ok": True}, RANDY_GET=lambda *args, **kwargs: {"summary": {"unique_visitors": 1, "sessions": 1, "page_views": 1, "workflow_starts": 0, "completed_analyses": 0, "exports": 0, "builder_handoffs": 0}, "referrers": [], "devices": [], "pages": [], "countries": [], "funnel": [], "workflows": [], "handoffs": {"session_rate": None, "by_origin": []}, "operations": {"submitted": 0, "completed": 0, "failed": 0, "success_rate": None, "failures": []}, "legacy": {"note": "Unavailable"}})
    app.register_blueprint(bp)
    app.after_request(track_page_view)
    app.after_request(track_server_outcome)
    @app.get("/")
    def home(): return "<html>public</html>"
    return app


def test_anonymous_visitor_and_session_are_deduplicated():
    events = []; client = make_app(events).test_client()
    client.get("/"); client.get("/")
    assert len(events) == 2
    assert events[0]["visitor_id"] == events[1]["visitor_id"]
    assert events[0]["session_id"] == events[1]["session_id"]
    assert all(e["event_type"] == "page_view" for e in events)


def test_admin_login_protection_and_dashboard_rendering():
    events = []; client = make_app(events).test_client()
    assert client.get("/admin/analytics").status_code == 302
    assert client.post("/admin/login", data={"email": "admin@example.test", "password": "wrong"}).status_code == 200
    response = client.post("/admin/login", data={"email": "admin@example.test", "password": "correct-password"}, follow_redirects=True)
    assert response.status_code == 200
    assert b"Product analytics" in response.data


def test_geoip_uses_forwarded_browser_address(monkeypatch):
    from vlismod_analytics import _geoip
    captured = {}
    class Reply:
        def json(self): return {"country_code": "US", "country": "United States", "latitude": 38, "longitude": -97}
    monkeypatch.setenv("VLISMOD_USAGE_GEOIP", "1")
    monkeypatch.setattr("vlismod_analytics.requests.get", lambda url, timeout: captured.setdefault("url", url) and Reply())
    app = Flask(__name__)
    with app.test_request_context("/", headers={"X-Forwarded-For": "203.0.113.8, 10.0.0.1"}):
        assert _geoip()["country_code"] == "US"
    assert captured["url"].endswith("/203.0.113.8")


def test_browser_events_accept_only_safe_feature_labels():
    events = []; client = make_app(events).test_client()
    assert client.post("/analytics/event", json={"event_type": "workflow_started", "feature": "protein_query", "event_id": "safe"}).status_code == 200
    assert client.post("/analytics/event", json={"event_type": "workflow_started", "feature": "PDB-1ABC", "event_id": "unsafe"}).status_code == 400
