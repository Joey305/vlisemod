from __future__ import annotations

import os
import sqlite3
from functools import wraps
from pathlib import Path
from typing import Any, Callable

from flask import Blueprint, current_app, jsonify, request
from werkzeug.exceptions import HTTPException

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DB_PATH = PROJECT_ROOT / "viral_data.db"

vlismod_data_bp = Blueprint("vlismod_data", __name__, url_prefix="/api/vlismod")

REQUIRED_TABLES = (
    "ligand_atoms",
    "Ligand_Atoms_Smiles",
    "Functional_GROUPED",
    "Ligand_Arp_Diagram",
    "Functional_Group_Atoms",
    "Arpeggio_Contacts_Data",
    "RUPLEY_SASA_DATA",
    "SMILES_MAP_PDB",
)


def _configured_token() -> str:
    return os.environ.get("VLISMOD_API_TOKEN", "").strip()


def _db_path() -> Path:
    configured = os.environ.get("VLISMOD_DB_PATH", "").strip()
    if configured:
        return Path(configured).expanduser()
    return DEFAULT_DB_PATH


def _row_dict(row: sqlite3.Row) -> dict[str, Any]:
    return {key: row[key] for key in row.keys()}


def _connect() -> sqlite3.Connection:
    db_path = _db_path()
    if not db_path.exists():
        raise FileNotFoundError(f"Database not found: {db_path}")
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    return conn


def _json_error(message: str, status_code: int):
    return jsonify({"ok": False, "error": message}), status_code


def require_token(view: Callable):
    @wraps(view)
    def wrapped(*args, **kwargs):
        token = _configured_token()
        if not token:
            return _json_error("VLISMOD API token is not configured.", 500)

        auth_header = request.headers.get("Authorization", "")
        if not auth_header.startswith("Bearer "):
            return _json_error("Unauthorized.", 401)

        provided = auth_header.split(" ", 1)[1].strip()
        if provided != token:
            return _json_error("Unauthorized.", 401)

        return view(*args, **kwargs)

    return wrapped


def _required_arg(name: str) -> str:
    value = str(request.args.get(name, "")).strip()
    if not value:
        raise ValueError(f"Missing required parameter: {name}")
    return value


def _fetch_scalar_list(query: str, params: tuple[Any, ...]) -> list[Any]:
    with _connect() as conn:
        rows = conn.execute(query, params).fetchall()
    return [row[0] for row in rows]


def _fetch_rows(query: str, params: tuple[Any, ...]) -> list[sqlite3.Row]:
    with _connect() as conn:
        return conn.execute(query, params).fetchall()


@vlismod_data_bp.errorhandler(HTTPException)
def _json_http_error(error: HTTPException):
    return _json_error(error.description or error.name, error.code or 500)


@vlismod_data_bp.errorhandler(FileNotFoundError)
def _json_missing_file(error: FileNotFoundError):
    return _json_error(str(error), 404)


@vlismod_data_bp.errorhandler(sqlite3.OperationalError)
def _json_sqlite_error(error: sqlite3.OperationalError):
    message = str(error)
    status = 404 if "no such table" in message.lower() else 500
    return _json_error(message, status)


@vlismod_data_bp.errorhandler(ValueError)
def _json_value_error(error: ValueError):
    return _json_error(str(error), 400)


@vlismod_data_bp.errorhandler(Exception)
def _json_unhandled_error(error: Exception):
    current_app.logger.exception("Unhandled V-LiSEMOD RANDY route error")
    return _json_error("Internal server error.", 500)


@vlismod_data_bp.get("/health")
def vlismod_health():
    db_path = _db_path()
    return jsonify(
        {
            "ok": True,
            "service": "vlismod-data",
            "db_path": str(db_path),
            "db_exists": db_path.exists(),
            "auth_configured": bool(_configured_token()),
        }
    )


@vlismod_data_bp.get("/db-health")
@require_token
def vlismod_db_health():
    db_path = _db_path()
    table_status: dict[str, bool] = {}

    with _connect() as conn:
        existing = {
            row["name"]
            for row in conn.execute(
                "SELECT name FROM sqlite_master WHERE type = 'table' AND name IN ({})".format(
                    ",".join("?" for _ in REQUIRED_TABLES)
                ),
                REQUIRED_TABLES,
            ).fetchall()
        }
    for table in REQUIRED_TABLES:
        table_status[table] = table in existing

    return jsonify(
        {
            "ok": True,
            "db_path": str(db_path),
            "db_exists": db_path.exists(),
            "required_tables": table_status,
        }
    )


@vlismod_data_bp.get("/viruses")
@require_token
def get_viruses():
    viruses = _fetch_scalar_list(
        "SELECT DISTINCT virus_name FROM ligand_atoms ORDER BY virus_name",
        (),
    )
    return jsonify(viruses=viruses)


@vlismod_data_bp.get("/pdb-codes")
@require_token
def get_pdb_codes():
    virus_name = _required_arg("virus_name")
    pdb_codes = _fetch_scalar_list(
        """
        SELECT DISTINCT pdb_id
        FROM ligand_Atoms_Smiles
        WHERE virus_name = ?
        ORDER BY pdb_id
        """,
        (virus_name,),
    )
    return jsonify(pdb_codes=pdb_codes)


@vlismod_data_bp.get("/ligands")
@require_token
def get_ligands():
    pdb_code = _required_arg("pdb_code")
    rows = _fetch_rows(
        """
        SELECT ligand, MIN(chain) AS chain,
               MAX(
                   CASE
                       WHEN EXISTS (
                           SELECT 1
                           FROM Functional_GROUPED fg
                           WHERE fg.pdb_id = ligand_atoms.pdb_id
                             AND fg.ligand = ligand_atoms.ligand
                             AND fg.smiles IS NOT NULL
                             AND fg.smiles != ''
                       ) THEN 1 ELSE 0
                   END
               ) AS has_smiles
        FROM ligand_atoms
        WHERE pdb_id = ?
        GROUP BY ligand
        ORDER BY ligand
        """,
        (pdb_code,),
    )
    ligands = [
        {"ligand": row["ligand"], "chain": row["chain"], "has_smiles": row["has_smiles"]}
        for row in rows
    ]
    return jsonify(ligands=ligands)


@vlismod_data_bp.get("/functional-groups/check")
@require_token
def check_functional_groups():
    pdb_code = _required_arg("pdb_code")
    rows = _fetch_rows(
        """
        SELECT COUNT(*) AS functional_group_count
        FROM Functional_Group_Atoms
        WHERE pdb_id = ?
        """,
        (pdb_code,),
    )
    count = int(rows[0]["functional_group_count"]) if rows else 0
    return jsonify({"has_functional_groups": count > 0})


@vlismod_data_bp.get("/ligands/list")
@require_token
def get_ligands_list():
    ligands = _fetch_scalar_list(
        "SELECT DISTINCT ligand FROM Ligand_Atoms_Smiles ORDER BY ligand",
        (),
    )
    return jsonify(ligands=ligands)


@vlismod_data_bp.get("/viruses/by-ligand")
@require_token
def get_viruses_by_ligand():
    ligand_code = _required_arg("ligand_code")
    viruses = _fetch_scalar_list(
        """
        SELECT DISTINCT Virus_Name
        FROM Arpeggio_Contacts_Data
        WHERE Ligand = ?
        ORDER BY Virus_Name
        """,
        (ligand_code,),
    )
    return jsonify(viruses=viruses)


@vlismod_data_bp.get("/pdb-residues/by-ligand")
@require_token
def get_pdb_residue_by_ligand():
    ligand_code = _required_arg("ligand_code")
    rows = _fetch_rows(
        """
        SELECT DISTINCT pdb_id, chain, ligand_id
        FROM Ligand_Arp_Diagram
        WHERE ligand = ?
        ORDER BY pdb_id, chain, ligand_id
        """,
        (ligand_code,),
    )
    pairs = [
        {"pdb_id": row["pdb_id"], "chain": row["chain"], "ligand_id": row["ligand_id"]}
        for row in rows
    ]
    return jsonify(pairs=pairs)


@vlismod_data_bp.get("/pdb-mapping")
@require_token
def get_pdb_mapping():
    ligand_code = _required_arg("ligand_code")
    rows = _fetch_rows(
        """
        SELECT DISTINCT pdb_id, chain, ligand_id, virus_name, ligand
        FROM Ligand_Atoms_Smiles
        WHERE ligand = ?
        ORDER BY pdb_id, chain, ligand_id
        """,
        (ligand_code,),
    )

    pdb_mapping: dict[str, dict[str, Any]] = {}
    for row in rows:
        unique_key = f"{row['pdb_id']}-{row['ligand_id']}-{row['chain']}"
        if unique_key not in pdb_mapping:
            pdb_mapping[unique_key] = {
                "pdb_id": row["pdb_id"],
                "ligand_id": row["ligand_id"],
                "chain": row["chain"],
                "virus_name": row["virus_name"],
                "ligand": row["ligand"],
            }

    return jsonify(pdb_mapping=pdb_mapping)


@vlismod_data_bp.get("/sasa-chains")
@require_token
def get_sasa_chains():
    pdb_code = _required_arg("pdb_code")
    ligand_name = _required_arg("ligand_name")
    rows = _fetch_rows(
        """
        SELECT atom_id, chain
        FROM RUPLEY_SASA_DATA
        WHERE pdb_id = ? AND ligand = ?
        """,
        (pdb_code, ligand_name),
    )
    sasa_chains = [[row["atom_id"], row["chain"]] for row in rows]
    return jsonify(sasa_chains)
