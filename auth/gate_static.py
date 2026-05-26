from flask import Blueprint, send_from_directory, abort
from flask_login import login_required
from pathlib import Path

protected_static = Blueprint("protected_static", __name__)
STATIC_DIR = Path(__file__).resolve().parent / "static"

@protected_static.route("/<path:filename>")
@login_required
def static_file(filename):
    target = STATIC_DIR / filename
    if target.exists():
        return send_from_directory(STATIC_DIR, filename)
    return abort(404)
