import os
import sqlite3

from flask import Flask


def test_event_ingestion_is_idempotent_and_rolls_up(tmp_path, monkeypatch):
    monkeypatch.setenv("RANDY_BACKUP_TOKEN", "test-token")
    monkeypatch.setenv("PROTAC_BACKUP_DIR", str(tmp_path))
    from RANDY.vlismod_analytics_routes import bp
    app = Flask(__name__); app.register_blueprint(bp, url_prefix="/backup/vlismod")
    client, headers = app.test_client(), {"Authorization": "Bearer test-token"}
    event = {"event_id": "id-1", "event_type": "page_view", "occurred_at": "2026-10-06T12:00:00+00:00", "visitor_id": "visitor-1", "session_id": "session-1", "path": "/", "referrer": "direct", "device_class": "desktop"}
    assert client.post("/backup/vlismod/analytics/events", json=event, headers=headers).json["inserted"]
    assert not client.post("/backup/vlismod/analytics/events", json=event, headers=headers).json["inserted"]
    payload = client.get("/backup/vlismod/analytics/rollup?period=all", headers=headers).json
    assert payload["summary"]["unique_visitors"] == 1
    assert payload["summary"]["page_views"] == 1


def test_receiver_rejects_scientific_or_arbitrary_feature_values(tmp_path, monkeypatch):
    monkeypatch.setenv("RANDY_BACKUP_TOKEN", "test-token")
    monkeypatch.setenv("PROTAC_BACKUP_DIR", str(tmp_path))
    from RANDY.vlismod_analytics_routes import bp
    app = Flask(__name__); app.register_blueprint(bp, url_prefix="/backup/vlismod")
    payload = {"event_id": "id-unsafe", "event_type": "workflow_started", "feature": "PDB: 6M0J", "occurred_at": "2026-10-06T12:00:00+00:00", "visitor_id": "visitor-1", "session_id": "session-1", "path": "/", "referrer": "direct", "device_class": "desktop"}
    response = app.test_client().post("/backup/vlismod/analytics/events", json=payload, headers={"Authorization": "Bearer test-token"})
    assert response.status_code == 400


def test_receiver_persists_only_valid_opaque_handoff_id(tmp_path, monkeypatch):
    monkeypatch.setenv("RANDY_BACKUP_TOKEN", "test-token")
    monkeypatch.setenv("PROTAC_BACKUP_DIR", str(tmp_path))
    from RANDY.vlismod_analytics_routes import bp
    app = Flask(__name__); app.register_blueprint(bp, url_prefix="/backup/vlismod")
    handoff_id = "1b99df97-dde1-4b03-8ba7-c92af1b43da1"
    payload = {"event_id": "id-handoff", "event_type": "companion_handoff", "feature": "builder_from_landing", "handoff_id": handoff_id, "occurred_at": "2026-10-06T12:00:00+00:00", "visitor_id": "visitor-1", "session_id": "session-1", "path": "/", "referrer": "direct", "device_class": "desktop"}
    client = app.test_client(); headers = {"Authorization": "Bearer test-token"}
    assert client.post("/backup/vlismod/analytics/events", json=payload, headers=headers).status_code == 200
    assert sqlite3.connect(tmp_path / "vlismod_analytics.sqlite3").execute("select handoff_id from vlismod_events").fetchone()[0] == handoff_id
    assert client.post("/backup/vlismod/analytics/events", json={**payload, "event_id": "id-bad", "handoff_id": "PDB-1ABC"}, headers=headers).status_code == 400
