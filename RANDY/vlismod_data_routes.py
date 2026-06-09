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

REQUIRED_TABLES = (
    "ligand_atoms",
    "Ligand_Atoms_Smiles",
    "Ligand_Synonyms",
    "Functional_GROUPED",
    "Ligand_Arp_Diagram",
    "Functional_Group_Atoms",
    "Arpeggio_Contacts_Data",
    "RUPLEY_SASA_DATA",
    "SMILES_MAP_PDB",
    "Virus_Proteins",
)


def _configured_token() -> str:
    return (
        os.environ.get("VLISMOD_API_TOKEN", "").strip()
        or os.environ.get("RANDY_BACKUP_TOKEN", "").strip()
        or os.environ.get("PROTAC_BACKUP_TOKEN", "").strip()
    )


def _db_path() -> Path:
    configured = os.environ.get("VLISMOD_DB_PATH", "").strip()
    if configured:
        return Path(configured).expanduser()
    return DEFAULT_DB_PATH


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


def _required_json_arg(payload: dict[str, Any], name: str) -> str:
    value = str(payload.get(name, "")).strip()
    if not value:
        raise ValueError(f"Missing required JSON field: {name}")
    return value


def _json_flag(payload: dict[str, Any], name: str) -> bool:
    return bool(payload.get(name))


def _fetch_scalar_list(query: str, params: tuple[Any, ...]) -> list[Any]:
    with _connect() as conn:
        rows = conn.execute(query, params).fetchall()
    return [row[0] for row in rows]


def _fetch_rows(query: str, params: tuple[Any, ...]) -> list[sqlite3.Row]:
    with _connect() as conn:
        return conn.execute(query, params).fetchall()


def _table_exists(conn: sqlite3.Connection, table_name: str) -> bool:
    row = conn.execute(
        "SELECT 1 FROM sqlite_master WHERE type='table' AND name = ?",
        (table_name,),
    ).fetchone()
    return row is not None


def _protacability_tables_available(conn: sqlite3.Connection) -> bool:
    required = (
        "protacability_assessment",
        "protacability_lysine_proximity",
        "protacability_ligand_inventory",
    )
    return all(_table_exists(conn, table_name) for table_name in required)


def _available_export_data_sets(conn: sqlite3.Connection) -> list[str]:
    data_sets = [
        "Solvent Exposed Atoms",
        "Ligand Atoms",
        "Binding Pocket",
        "Smiles and Functional Groups",
        "Interatomic Interactions",
        "Functional Group Atoms",
        "Smiles & PDB Mapping",
    ]
    if _protacability_tables_available(conn):
        data_sets.extend([
            "PROTACability Assessment",
            "PROTACability Lysine Proximity",
            "PROTACability Ligand Inventory",
        ])
        if _table_exists(conn, "protacability_warhead_linkability"):
            data_sets.append("PROTACability Warhead Linkability")
        if _table_exists(conn, "protacability_degrader_readiness"):
            data_sets.append("PROTACability Degrader Readiness")
    return data_sets


def create_vlismod_blueprint(blueprint_name: str, url_prefix: str) -> Blueprint:
    bp = Blueprint(blueprint_name, __name__, url_prefix=url_prefix)

    @bp.errorhandler(HTTPException)
    def _json_http_error(error: HTTPException):
        return _json_error(error.description or error.name, error.code or 500)

    @bp.errorhandler(FileNotFoundError)
    def _json_missing_file(error: FileNotFoundError):
        return _json_error(str(error), 404)

    @bp.errorhandler(sqlite3.OperationalError)
    def _json_sqlite_error(error: sqlite3.OperationalError):
        message = str(error)
        status = 404 if "no such table" in message.lower() else 500
        return _json_error(message, status)

    @bp.errorhandler(ValueError)
    def _json_value_error(error: ValueError):
        return _json_error(str(error), 400)

    @bp.errorhandler(Exception)
    def _json_unhandled_error(error: Exception):
        current_app.logger.exception("Unhandled V-LiSEMOD RANDY route error")
        return _json_error("Internal server error.", 500)

    @bp.get("/health")
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

    @bp.get("/db-health")
    @require_token
    def vlismod_db_health():
        db_path = _db_path()
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

        table_status = {table: table in existing for table in REQUIRED_TABLES}
        return jsonify(
            {
                "ok": True,
                "db_path": str(db_path),
                "db_exists": db_path.exists(),
                "required_tables": table_status,
            }
        )

    @bp.get("/viruses")
    @require_token
    def get_viruses():
        viruses = _fetch_scalar_list(
            "SELECT DISTINCT virus_name FROM ligand_atoms ORDER BY virus_name",
            (),
        )
        return jsonify(viruses=viruses)

    @bp.get("/pdb-codes")
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

    @bp.get("/ligands")
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

    @bp.get("/functional-groups/check")
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

    @bp.get("/ligands/list")
    @require_token
    def get_ligands_list():
        ligands = _fetch_scalar_list(
            "SELECT DISTINCT ligand FROM Ligand_Atoms_Smiles ORDER BY ligand",
            (),
        )
        return jsonify(ligands=ligands)

    @bp.get("/viruses/by-ligand")
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

    @bp.get("/pdb-residues/by-ligand")
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

    @bp.get("/pdb-mapping")
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

    @bp.get("/sasa-chains")
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

    @bp.get("/page-metadata")
    @require_token
    def get_page_metadata():
        with _connect() as conn:
            return jsonify(
                {
                    "available_export_data_sets": _available_export_data_sets(conn),
                    "protacability_data_available": _protacability_tables_available(conn),
                }
            )

    @bp.get("/ligands/with-synonyms")
    @require_token
    def get_ligands_with_synonyms():
        rows = _fetch_rows(
            """
            SELECT ligand, synonym
            FROM Ligand_Synonyms
            ORDER BY ligand, synonym
            """,
            (),
        )
        return jsonify(
            [
                {"ligand_code": row["ligand"], "synonym": row["synonym"]}
                for row in rows
            ]
        )

    @bp.get("/ligand-info")
    @require_token
    def get_ligand_info():
        ligand_code = _required_arg("ligand_code")
        rows = _fetch_rows(
            """
            SELECT ligand, pdb_id, smiles, molecular_weight
            FROM Ligand_Atoms_Smiles
            WHERE ligand = ?
            LIMIT 1
            """,
            (ligand_code,),
        )
        if not rows:
            return _json_error("Ligand not found", 404)
        row = rows[0]
        return jsonify(
            {
                "ligand": row["ligand"],
                "pdb_id": row["pdb_id"],
                "smiles": row["smiles"],
                "molecular_weight": row["molecular_weight"],
            }
        )

    @bp.get("/ligand-options")
    @require_token
    def get_ligand_options():
        rows = _fetch_rows(
            """
            SELECT virus_name, pdb_id, ligand, chain, ligand_id
            FROM Ligand_Arp_Diagram
            ORDER BY virus_name, pdb_id, ligand, chain, ligand_id
            """,
            (),
        )
        return jsonify(
            {
                "options": [
                    {
                        "value": f"{row['pdb_id']}-{row['ligand']}",
                        "text": f"{row['virus_name']}, {row['pdb_id']}, {row['ligand']}, {row['chain']}, {row['ligand_id']}",
                    }
                    for row in rows
                ]
            }
        )

    @bp.get("/ligand-smiles")
    @require_token
    def get_ligand_smiles():
        ligand_code = _required_arg("ligand_id")
        rows = _fetch_rows(
            """
            SELECT smiles
            FROM Ligand_Atoms_Smiles
            WHERE ligand = ?
            LIMIT 1
            """,
            (ligand_code,),
        )
        if not rows or not rows[0]["smiles"]:
            return _json_error("SMILES not found", 404)
        return jsonify({"ligand_id": ligand_code, "smiles": rows[0]["smiles"]})

    @bp.get("/interaction-records")
    @require_token
    def get_interaction_records():
        pdb_id = _required_arg("pdb_id")
        ligand = _required_arg("ligand")
        ligand_id = _required_arg("ligand_id")
        chain = _required_arg("chain")
        rows = _fetch_rows(
            """
            SELECT A.pdb_id, A.ligand, A.chain, A.Contact, A.Distance, A.exact_atom, A.atom_id,
                   A.residue, A.residue_number, A.residue_atom, A.residue_chain,
                   S.smiles_atom_index, A.virus_name, A.ligand_id
            FROM Arpeggio_Contacts_Data A
            LEFT JOIN SMILES_MAP_PDB S
              ON A.atom_id = S.atom_id
             AND A.pdb_id = S.pdb_id
             AND A.chain = S.chain
             AND A.exact_atom = S.exact_atom
            WHERE A.pdb_id = ?
              AND A.ligand_id = ?
              AND A.chain = ?
              AND A.ligand = ?
            """,
            (pdb_id, ligand_id, chain, ligand),
        )
        return jsonify({"records": [dict(row) for row in rows]})

    @bp.post("/ligand-interactions/compare")
    @require_token
    def compare_ligand_interactions():
        payload = request.get_json(silent=True) or {}
        selected_pdbs = payload.get("pdb_ids") or []
        ligand = _required_json_arg(payload, "ligand")

        interactions_data = []
        smiles_interactions_data = []
        with _connect() as conn:
            for unique_key in selected_pdbs:
                pdb_id, ligand_id, chain = str(unique_key).split("-")
                rows = conn.execute(
                    """
                    SELECT A.pdb_id, A.ligand, A.chain, A.Contact, A.Distance, A.exact_atom, A.atom_id,
                           A.residue, A.residue_number, A.residue_atom, A.residue_chain,
                           S.smiles_atom_index, A.virus_name, A.ligand_id
                    FROM Arpeggio_Contacts_Data A
                    LEFT JOIN SMILES_MAP_PDB S
                      ON A.atom_id = S.atom_id
                     AND A.pdb_id = S.pdb_id
                     AND A.chain = S.chain
                     AND A.exact_atom = S.exact_atom
                    WHERE A.pdb_id = ?
                      AND A.ligand_id = ?
                      AND A.chain = ?
                      AND A.ligand = ?
                    """,
                    (pdb_id, ligand_id, chain, ligand),
                ).fetchall()

                cleaned_rows = []
                smiles_rows = []
                for row in rows:
                    row_dict = dict(row)
                    contact = row_dict.get("Contact")
                    if contact == "proximal":
                        pass
                    else:
                        if contact == "weak_polar":
                            row_dict["Contact"] = "polar"
                        elif contact == "vdw_clash":
                            row_dict["Contact"] = "vdw"
                        atom_id = row_dict.get("atom_id")
                        if atom_id not in (None, "", "N/A"):
                            try:
                                row_dict["atom_id"] = int(atom_id)
                                cleaned_rows.append(row_dict)
                            except (TypeError, ValueError):
                                pass

                    smiles_atom_index = row_dict.get("smiles_atom_index")
                    if smiles_atom_index not in (None, "", "N/A"):
                        try:
                            smiles_rows.append(
                                {
                                    "pdb_id": row_dict.get("pdb_id"),
                                    "Contact": row_dict.get("Contact"),
                                    "smiles_atom_index": int(float(smiles_atom_index)),
                                }
                            )
                        except (TypeError, ValueError):
                            pass

                if cleaned_rows:
                    interactions_data.append(
                        {
                            "pdb_id": pdb_id,
                            "virus_name": cleaned_rows[0].get("virus_name"),
                            "ligand_id": int(cleaned_rows[0].get("ligand_id")),
                            "interactions": cleaned_rows,
                        }
                    )
                if smiles_rows:
                    smiles_interactions_data.append(
                        {"pdb_id": pdb_id, "interactions": smiles_rows}
                    )

        return jsonify(
            {
                "interactions_data": interactions_data,
                "smiles_interactions_data": smiles_interactions_data,
            }
        )

    @bp.get("/virus-proteins/virus-names")
    @require_token
    def get_virus_names():
        return jsonify(
            _fetch_scalar_list(
                "SELECT DISTINCT virus_name FROM Virus_Proteins ORDER BY virus_name",
                (),
            )
        )

    @bp.get("/virus-proteins/protein-types")
    @require_token
    def get_protein_types():
        return jsonify(
            _fetch_scalar_list(
                "SELECT DISTINCT protein FROM Virus_Proteins ORDER BY protein",
                (),
            )
        )

    @bp.post("/virus-proteins/pdbs")
    @require_token
    def get_pdbs_for_virus_protein():
        payload = request.get_json(silent=True) or {}
        virus_name = _required_json_arg(payload, "virus_name")
        protein_types = payload.get("protein_types") or []
        ligand_filter = str(payload.get("ligand", "")).strip()
        if not protein_types:
            raise ValueError("Missing required JSON field: protein_types")

        placeholders = ", ".join(["?"] * len(protein_types))
        query = f"SELECT DISTINCT pdb_id FROM Virus_Proteins WHERE virus_name = ? AND protein IN ({placeholders})"
        params: list[Any] = [virus_name, *protein_types]
        if ligand_filter:
            query += """
                AND pdb_id IN (
                    SELECT pdb_id FROM Ligand_Arp_Diagram WHERE ligand = ? OR ligand IN (
                        SELECT synonym FROM Ligand_Synonyms WHERE ligand = ?
                    )
                )
            """
            params.extend([ligand_filter, ligand_filter])

        rows = _fetch_rows(query, tuple(params))
        return jsonify({"pdb_codes": [row["pdb_id"] for row in rows]})

    @bp.post("/export-data")
    @require_token
    def export_data():
        payload = request.get_json(silent=True) or {}
        pdb_codes = payload.get("pdb_codes") or []
        data_sets = payload.get("data_sets") or []
        if not pdb_codes or not data_sets:
            raise ValueError("Missing required JSON fields: pdb_codes and data_sets")

        query_map = {
            "Solvent Exposed Atoms": ("SELECT * FROM RUPLEY_SASA_data WHERE pdb_id IN ({placeholders})", "pdb_id"),
            "Ligand Atoms": ("SELECT * FROM ligand_atoms WHERE pdb_id IN ({placeholders})", "pdb_id"),
            "Binding Pocket": ("SELECT * FROM receptor_binding_pocket WHERE pdb_id IN ({placeholders})", "pdb_id"),
            "Smiles and Functional Groups": ("SELECT * FROM Ligand_Atoms_Smiles WHERE pdb_id IN ({placeholders})", "pdb_id"),
            "Interatomic Interactions": ("SELECT * FROM Arpeggio_Contacts_Data WHERE pdb_id IN ({placeholders})", "pdb_id"),
            "Functional Group Atoms": ("SELECT * FROM Functional_Group_Atoms WHERE pdb_id IN ({placeholders})", "pdb_id"),
            "Smiles & PDB Mapping": ("SELECT * FROM SMILES_MAP_PDB WHERE pdb_id IN ({placeholders})", "pdb_id"),
            "PROTACability Assessment": ("SELECT * FROM protacability_assessment WHERE pdb_code IN ({placeholders})", "pdb_code"),
            "PROTACability Lysine Proximity": ("SELECT * FROM protacability_lysine_proximity WHERE pdb_code IN ({placeholders})", "pdb_code"),
            "PROTACability Ligand Inventory": ("SELECT * FROM protacability_ligand_inventory WHERE pdb_code IN ({placeholders})", "pdb_code"),
            "PROTACability Warhead Linkability": ("SELECT * FROM protacability_warhead_linkability WHERE pdb_code IN ({placeholders})", "pdb_code"),
            "PROTACability Degrader Readiness": ("SELECT * FROM protacability_degrader_readiness WHERE pdb_code IN ({placeholders})", "pdb_code"),
        }

        with _connect() as conn:
            allowed_sets = set(_available_export_data_sets(conn))
            invalid = [name for name in data_sets if name not in allowed_sets or name not in query_map]
            if invalid:
                raise ValueError(f"Unsupported data sets requested: {invalid}")

            result: dict[str, Any] = {"data_sets": {}}
            for data_set in data_sets:
                query_template, _column = query_map[data_set]
                placeholders = ", ".join(["?"] * len(pdb_codes))
                query = query_template.format(placeholders=placeholders)
                rows = conn.execute(query, tuple(pdb_codes)).fetchall()
                result["data_sets"][data_set] = [dict(row) for row in rows]

        return jsonify(result)

    @bp.get("/protacability/source")
    @require_token
    def get_protacability_source():
        with _connect() as conn:
            data_available = _protacability_tables_available(conn)
            if not data_available:
                return jsonify(
                    {
                        "data_available": False,
                        "assessment_rows": [],
                        "readiness_rows": [],
                        "warhead_rows": [],
                        "lysine_rows": [],
                        "ligand_inventory": [],
                    }
                )

            pdb_code = str(request.args.get("pdb_code", "")).strip()
            virus_name = str(request.args.get("virus_name", "")).strip()
            protein_type = str(request.args.get("protein_type", "")).strip()
            include_lysine = str(request.args.get("include_lysine", "")).strip().lower() in {"1", "true", "yes"}
            include_inventory = str(request.args.get("include_inventory", "")).strip().lower() in {"1", "true", "yes"}

            clauses = []
            params: list[Any] = []
            if pdb_code:
                clauses.append("pdb_code = ?")
                params.append(pdb_code)
            if virus_name:
                clauses.append("virus_name = ?")
                params.append(virus_name)
            if protein_type:
                clauses.append("protein_type = ?")
                params.append(protein_type)
            where_sql = f" WHERE {' AND '.join(clauses)}" if clauses else ""

            assessment_rows = [
                dict(row)
                for row in conn.execute(
                    f"SELECT * FROM protacability_assessment{where_sql}",
                    tuple(params),
                ).fetchall()
            ]
            pdb_codes = sorted({row.get("pdb_code") for row in assessment_rows if row.get("pdb_code")})

            def _load_optional_rows(table_name: str) -> list[dict[str, Any]]:
                if not _table_exists(conn, table_name):
                    return []
                if not pdb_codes:
                    return [dict(row) for row in conn.execute(f"SELECT * FROM {table_name}").fetchall()]
                placeholders = ", ".join(["?"] * len(pdb_codes))
                return [
                    dict(row)
                    for row in conn.execute(
                        f"SELECT * FROM {table_name} WHERE pdb_code IN ({placeholders})",
                        tuple(pdb_codes),
                    ).fetchall()
                ]

            readiness_rows = _load_optional_rows("protacability_degrader_readiness")
            warhead_rows = _load_optional_rows("protacability_warhead_linkability")
            lysine_rows = _load_optional_rows("protacability_lysine_proximity") if include_lysine else []
            ligand_inventory = _load_optional_rows("protacability_ligand_inventory") if include_inventory else []

            return jsonify(
                {
                    "data_available": True,
                    "assessment_rows": assessment_rows,
                    "readiness_rows": readiness_rows,
                    "warhead_rows": warhead_rows,
                    "lysine_rows": lysine_rows,
                    "ligand_inventory": ligand_inventory,
                }
            )

    @bp.get("/protacability/raw-table")
    @require_token
    def get_protacability_raw_table():
        raw_export = _required_arg("raw_export")
        table_map = {
            "PROTACability Assessment": "protacability_assessment",
            "PROTACability Lysine Proximity": "protacability_lysine_proximity",
            "PROTACability Ligand Inventory": "protacability_ligand_inventory",
            "PROTACability Warhead Linkability": "protacability_warhead_linkability",
            "PROTACability Degrader Readiness": "protacability_degrader_readiness",
        }
        table_name = table_map.get(raw_export)
        if not table_name:
            raise ValueError("Unknown PROTACability export selection.")
        with _connect() as conn:
            if not _table_exists(conn, table_name):
                return _json_error(f"{raw_export} has not been imported yet.", 404)
            rows = [dict(row) for row in conn.execute(f"SELECT * FROM {table_name}").fetchall()]
        return jsonify({"table_name": table_name, "rows": rows})

    @bp.post("/ligand-images-data")
    @require_token
    def get_ligand_images_data():
        payload = request.get_json(silent=True) or {}
        virus_name = _required_json_arg(payload, "virus_name")
        pdb_code = _required_json_arg(payload, "pdb_code")
        ligand_name = _required_json_arg(payload, "ligand_name")
        requested_chain = str(payload.get("chain", "")).strip() or None

        with _connect() as conn:
            smiles_rows = conn.execute(
                """
                SELECT DISTINCT virus_name, pdb_id, ligand, smiles
                FROM Functional_GROUPED
                WHERE virus_name = ? AND pdb_id = ? AND ligand = ?
                """,
                (virus_name, pdb_code, ligand_name),
            ).fetchall()

            chain_rows = conn.execute(
                """
                SELECT DISTINCT chain, ligand_id
                FROM Ligand_Arp_Diagram
                WHERE virus_name = ? AND pdb_id = ? AND ligand = ?
                ORDER BY chain, ligand_id
                """,
                (virus_name, pdb_code, ligand_name),
            ).fetchall()

            if not smiles_rows:
                return _json_error("No ligand image data found for the selected virus, PDB code, and ligand.", 404)

            chains_to_map = [requested_chain] if requested_chain else [
                row["chain"] for row in chain_rows if row["chain"]
            ]
            solvent_exposed_atom_map: dict[str, list[int]] = {}

            for chain in chains_to_map:
                sasa_rows = conn.execute(
                    """
                    SELECT atom_id, exact_atom
                    FROM RUPLEY_SASA_DATA
                    WHERE virus_name = ? AND pdb_id = ? AND ligand = ? AND chain = ?
                    """,
                    (virus_name, pdb_code, ligand_name, chain),
                ).fetchall()

                smiles_indices: list[int] = []
                for sasa_row in sasa_rows:
                    mapped_row = conn.execute(
                        """
                        SELECT smiles_atom_index
                        FROM SMILES_MAP_PDB
                        WHERE virus_name = ?
                          AND pdb_id = ?
                          AND ligand = ?
                          AND chain = ?
                          AND atom_id = ?
                          AND exact_atom = ?
                        """,
                        (
                            virus_name,
                            pdb_code,
                            ligand_name,
                            chain,
                            sasa_row["atom_id"],
                            sasa_row["exact_atom"],
                        ),
                    ).fetchone()

                    if mapped_row and mapped_row["smiles_atom_index"] is not None:
                        smiles_indices.append(int(mapped_row["smiles_atom_index"]))

                solvent_exposed_atom_map[f"{pdb_code}|{ligand_name}|{chain}"] = sorted(set(smiles_indices))

        return jsonify(
            {
                "ok": True,
                "smiles_data": [
                    {
                        "virus_name": row["virus_name"],
                        "pdb_id": row["pdb_id"],
                        "ligand": row["ligand"],
                        "smiles": row["smiles"],
                    }
                    for row in smiles_rows
                ],
                "chain_residue_data": [
                    {"chain": row["chain"], "ligand_id": row["ligand_id"]}
                    for row in chain_rows
                ],
                "solvent_exposed_atom_map": solvent_exposed_atom_map,
            }
        )

    @bp.post("/pymol-session-data")
    @require_token
    def get_pymol_session_data():
        payload = request.get_json(silent=True) or {}
        pdb_code = _required_json_arg(payload, "pdb_code")
        ligand_name = _required_json_arg(payload, "ligand_name")
        requested_chain = str(payload.get("chain", "")).strip() or None
        options = payload.get("options") or {}

        with _connect() as conn:
            chain_rows = conn.execute(
                """
                SELECT DISTINCT chain
                FROM ligand_atoms
                WHERE pdb_id = ? AND ligand = ?
                ORDER BY chain
                """,
                (pdb_code, ligand_name),
            ).fetchall()
            chains = [row["chain"] for row in chain_rows if row["chain"]]

            if requested_chain:
                if requested_chain not in chains:
                    return _json_error("Specified ligand chain was not found.", 404)
                ligand_chain = requested_chain
            else:
                if not chains:
                    return _json_error("No ligand chain found for the selected PDB and ligand.", 404)
                if len(chains) > 1:
                    return _json_error("Multiple chains found; please specify the chain.", 400)
                ligand_chain = chains[0]

            virus_row = conn.execute(
                """
                SELECT virus_name
                FROM ligand_atoms
                WHERE pdb_id = ? AND ligand = ? AND chain = ?
                LIMIT 1
                """,
                (pdb_code, ligand_name, ligand_chain),
            ).fetchone()
            virus_name = virus_row["virus_name"] if virus_row else None

            response_payload: dict[str, Any] = {
                "ok": True,
                "ligand_chain": ligand_chain,
                "functional_groups": {},
                "binding_pocket": [],
                "distal_atoms": [],
                "solvent_exposed_atoms": [],
                "hydrated_atoms": [],
                "rupley_sasa": [],
            }

            if _json_flag(options, "functional_groups") and virus_name:
                fg_rows = conn.execute(
                    """
                    SELECT functional_group, atom_id, exact_atom, atom_type
                    FROM Functional_Group_Atoms
                    WHERE virus_name = ?
                      AND pdb_id = ?
                      AND ligand = ?
                      AND chain = ?
                    ORDER BY functional_group, atom_id
                    """,
                    (virus_name, pdb_code, ligand_name, ligand_chain),
                ).fetchall()
                functional_groups: dict[str, list[dict[str, Any]]] = {}
                for row in fg_rows:
                    functional_groups.setdefault(row["functional_group"], []).append(
                        {
                            "atom_id": row["atom_id"],
                            "exact_atom": row["exact_atom"],
                            "atom_type": row["atom_type"],
                        }
                    )
                response_payload["functional_groups"] = functional_groups

            if _json_flag(options, "binding_pocket"):
                binding_rows = conn.execute(
                    """
                    SELECT residue_chain, residue_number
                    FROM receptor_binding_pocket
                    WHERE pdb_id = ?
                    ORDER BY residue_chain, residue_number
                    """,
                    (pdb_code,),
                ).fetchall()
                response_payload["binding_pocket"] = [
                    {
                        "residue_chain": row["residue_chain"],
                        "residue_number": row["residue_number"],
                    }
                    for row in binding_rows
                ]

            if _json_flag(options, "distal_atoms"):
                distal_rows = conn.execute(
                    """
                    SELECT chain, atom_id
                    FROM distal_atoms
                    WHERE pdb_id = ? AND ligand = ?
                    ORDER BY chain, atom_id
                    """,
                    (pdb_code, ligand_name),
                ).fetchall()
                response_payload["distal_atoms"] = [
                    {"chain": row["chain"], "atom_id": row["atom_id"]}
                    for row in distal_rows
                ]

            if _json_flag(options, "solvent_exposed_atoms"):
                solvent_rows = conn.execute(
                    """
                    SELECT atom_id, chain
                    FROM RUPLEY_SASA_DATA
                    WHERE pdb_id = ? AND ligand = ? AND chain = ?
                    ORDER BY chain, atom_id
                    """,
                    (pdb_code, ligand_name, ligand_chain),
                ).fetchall()
                response_payload["solvent_exposed_atoms"] = [
                    {"atom_id": row["atom_id"], "chain": row["chain"]}
                    for row in solvent_rows
                ]

            if _json_flag(options, "hydrated_atoms"):
                hydrated_rows = conn.execute(
                    """
                    SELECT chain, atom_id
                    FROM ligand_atoms
                    WHERE pdb_id = ? AND ligand = ? AND chain = ?
                    ORDER BY chain, atom_id
                    """,
                    (pdb_code, ligand_name, ligand_chain),
                ).fetchall()
                response_payload["hydrated_atoms"] = [
                    {"chain": row["chain"], "atom_id": row["atom_id"]}
                    for row in hydrated_rows
                ]

            if _json_flag(options, "rupley_sasa"):
                rupley_rows = conn.execute(
                    """
                    SELECT atom_id, chain
                    FROM RUPLEY_SASA_DATA
                    WHERE pdb_id = ? AND ligand = ? AND chain = ?
                    ORDER BY chain, atom_id
                    """,
                    (pdb_code, ligand_name, ligand_chain),
                ).fetchall()
                response_payload["rupley_sasa"] = [
                    {"atom_id": row["atom_id"], "chain": row["chain"]}
                    for row in rupley_rows
                ]

        return jsonify(response_payload)

    return bp


vlismod_data_bp = create_vlismod_blueprint("vlismod_data", "/api/vlismod")
vlismod_backup_bp = create_vlismod_blueprint("vlismod_backup", "/backup/vlismod")
