import os

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
