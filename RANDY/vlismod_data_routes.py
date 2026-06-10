from __future__ import annotations

import os
import re
import sqlite3
from collections import Counter, defaultdict
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


def _protacability_source_payload_local(
    *,
    pdb_code: str | None = None,
    virus_name: str | None = None,
    protein_type: str | None = None,
    include_lysine: bool = False,
    include_inventory: bool = False,
) -> dict[str, Any]:
    with _connect() as conn:
        data_available = _protacability_tables_available(conn)
        if not data_available:
            return {
                "data_available": False,
                "assessment_rows": [],
                "readiness_rows": [],
                "warhead_rows": [],
                "lysine_rows": [],
                "ligand_inventory": [],
            }

        pdb_code = str(pdb_code or "").strip()
        virus_name = str(virus_name or "").strip()
        protein_type = str(protein_type or "").strip()

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

        return {
            "data_available": True,
            "assessment_rows": assessment_rows,
            "readiness_rows": readiness_rows,
            "warhead_rows": warhead_rows,
            "lysine_rows": lysine_rows,
            "ligand_inventory": ligand_inventory,
        }


#SCRIPTS FOR PROTEIN QUERYING

data_set_queries = {
    "Solvent Exposed Atoms": "SELECT * FROM RUPLEY_SASA_data WHERE pdb_id IN ({placeholders})",
    "Ligand Atoms": "SELECT * FROM ligand_atoms WHERE pdb_id IN ({placeholders})",
    "Binding Pocket": "SELECT * FROM receptor_binding_pocket WHERE pdb_id IN ({placeholders})",
    "Smiles and Functional Groups": "SELECT * FROM Ligand_Atoms_Smiles WHERE pdb_id IN ({placeholders})",
    "Interatomic Interactions": "SELECT * FROM Arpeggio_Contacts_Data WHERE pdb_id IN ({placeholders})",
    "Functional Group Atoms": "SELECT * FROM Functional_Group_Atoms WHERE pdb_id IN ({placeholders})",
    "Smiles & PDB Mapping": "SELECT * FROM SMILES_MAP_PDB WHERE pdb_id IN ({placeholders})",
    "PROTACability Assessment": "SELECT * FROM protacability_assessment WHERE pdb_code IN ({placeholders})",
    "PROTACability Lysine Proximity": "SELECT * FROM protacability_lysine_proximity WHERE pdb_code IN ({placeholders})",
    "PROTACability Ligand Inventory": "SELECT * FROM protacability_ligand_inventory WHERE pdb_code IN ({placeholders})",
    "PROTACability Warhead Linkability": "SELECT * FROM protacability_warhead_linkability WHERE pdb_code IN ({placeholders})",
    "PROTACability Degrader Readiness": "SELECT * FROM protacability_degrader_readiness WHERE pdb_code IN ({placeholders})",
}

def connect_db():
    return _connect()


PROTACABILITY_REQUIRED_TABLES = (
    "protacability_assessment",
    "protacability_lysine_proximity",
    "protacability_ligand_inventory",
)

PROTACABILITY_OPTIONAL_TABLES = (
    "protacability_warhead_linkability",
    "protacability_degrader_readiness",
)

PROTACABILITY_SORT_COLUMNS = {
    "protacability_proxy_score": "protacability_proxy_score",
    "virus_name": "virus_name",
    "protein_type": "protein_type",
    "pdb_code": "pdb_code",
    "chain_id": "chain_id",
    "isoelectric_point": "isoelectric_point",
    "exposed_lys_fraction": "exposed_lys_fraction",
    "candidate_ligand_count": "candidate_ligand_count",
    "chain_length_aa": "chain_length_aa",
    "protacability_tier": "protacability_tier",
}

PROTACABILITY_GLYCAN_LIGANDS = {
    "NAG", "BMA", "MAN", "FUC", "GAL", "GLC", "SIA", "NDG", "BGC", "GLA",
    "GLCN", "A2G", "GCU", "XYL", "FUL", "FRU", "GME", "G7L", "G7O",
}

PROTACABILITY_COMMON_CONTEXT_LIGANDS = {
    "GOL", "PEG", "PG4", "PGE", "SO4", "PO4", "EDO", "DMS", "DMSO", "ACT",
    "CL", "NA", "MG", "CA", "ZN", "MN", "FE", "NH4", "MPD", "IPA", "SCN",
    "FMT", "MSE", "ACE", "TRS", "MES", "HEP", "SOR",
}

PROTACABILITY_VIEW_MODES = {"targets", "summary", "chains", "protein"}

PROTACABILITY_CHAIN_COLUMNS = [
    "virus_name",
    "protein_type",
    "pdb_code",
    "chain_id",
    "chain_length_aa",
    "candidate_ligand_count",
    "candidate_ligand_resnames",
    "lys_count",
    "exposed_lys_count",
    "exposed_lys_fraction",
    "lysine_surface_fraction",
    "isoelectric_point",
    "protein_ligand_druggability_proxy_score",
    "protacability_proxy_score",
    "protacability_tier",
    "linker_docking_site_annotation",
    "has_candidate_ligand",
    "has_exposed_lysine",
    "has_ligand_proximal_exposed_lysine",
    "notes",
]

PROTACABILITY_TIER_ORDER = {
    "High structural priority": 0,
    "Moderate structural priority": 1,
    "Low structural priority": 2,
    "Insufficient structural context": 3,
}

PROTACABILITY_SORT_OPTIONS = {
    "targets": {
        "readiness_priority_desc": ("readiness_rank_score", True),
        "ligand_priority_desc": ("ligand_priority", True),
        "best_score_desc": ("best_score", True),
        "warhead_score_desc": ("warhead_linkability_score", True),
        "avg_score_desc": ("avg_score", True),
        "pdb_count_desc": ("pdb_count", True),
        "protein_type_asc": ("protein_type", False),
    },
    "summary": {
        "readiness_priority_desc": ("readiness_rank_score", True),
        "ligand_priority_desc": ("ligand_priority", True),
        "best_score_desc": ("best_score", True),
        "warhead_score_desc": ("warhead_linkability_score", True),
        "avg_score_desc": ("avg_score", True),
        "exposed_fraction_desc": ("best_exposed_lys_fraction", True),
        "candidate_ligand_count_desc": ("informative_ligand_count", True),
        "pdb_code_asc": ("pdb_code", False),
        "protein_type_asc": ("protein_type", False),
    },
    "protein": {
        "readiness_priority_desc": ("readiness_rank_score", True),
        "ligand_priority_desc": ("ligand_priority", True),
        "best_score_desc": ("best_score", True),
        "warhead_score_desc": ("warhead_linkability_score", True),
        "avg_score_desc": ("avg_score", True),
        "ligand_context_count_desc": ("ligand_context_count", True),
        "pdb_count_desc": ("pdb_count", True),
        "protein_type_asc": ("protein_type", False),
    },
    "chains": {
        "readiness_priority_desc": ("readiness_rank_score", True),
        "protacability_proxy_score_desc": ("protacability_proxy_score", True),
        "ligand_priority_desc": ("ligand_priority", True),
        "warhead_score_desc": ("warhead_linkability_score", True),
        "exposed_fraction_desc": ("exposed_lys_fraction", True),
        "candidate_ligand_count_desc": ("informative_ligand_count", True),
        "pdb_code_asc": ("pdb_code", False),
        "protein_type_asc": ("display_protein_type", False),
    },
}


def connect_db_row():
    return _connect()


def protacability_tables_available(conn=None):
    owns_conn = False
    if conn is None:
        conn = connect_db()
        owns_conn = True
    try:
        cursor = conn.cursor()
        cursor.execute(
            "SELECT name FROM sqlite_master WHERE type='table' AND name IN (?, ?, ?)",
            PROTACABILITY_REQUIRED_TABLES
        )
        found = {row[0] for row in cursor.fetchall()}
        return all(table in found for table in PROTACABILITY_REQUIRED_TABLES)
    finally:
        if owns_conn:
            conn.close()


def _table_exists(conn, table_name):
    row = conn.execute(
        "SELECT 1 FROM sqlite_master WHERE type='table' AND name = ?",
        (table_name,),
    ).fetchone()
    return row is not None


def _table_columns(conn, table_name):
    if not _table_exists(conn, table_name):
        return set()
    return {row[1] for row in conn.execute(f"PRAGMA table_info({table_name})").fetchall()}


def _protacability_optional_table_names(conn):
    return {name for name in PROTACABILITY_OPTIONAL_TABLES if _table_exists(conn, name)}


def get_available_export_data_sets():
    data_sets = [
        "Solvent Exposed Atoms",
        "Ligand Atoms",
        "Binding Pocket",
        "Smiles and Functional Groups",
        "Interatomic Interactions",
        "Functional Group Atoms",
        "Smiles & PDB Mapping",
    ]
    if protacability_tables_available():
        data_sets.extend([
            "PROTACability Assessment",
            "PROTACability Lysine Proximity",
            "PROTACability Ligand Inventory",
        ])
        with connect_db() as conn:
            optional_tables = _protacability_optional_table_names(conn)
        if "protacability_warhead_linkability" in optional_tables:
            data_sets.append("PROTACability Warhead Linkability")
        if "protacability_degrader_readiness" in optional_tables:
            data_sets.append("PROTACability Degrader Readiness")
    return data_sets


def _normalize_multi_values(values):
    normalized = []
    for raw in values:
        if raw is None:
            continue
        if isinstance(raw, str):
            parts = [part.strip() for part in raw.split(",")]
        else:
            parts = [str(raw).strip()]
        for part in parts:
            if part:
                normalized.append(part)
    return normalized


def _parse_bool_filter(value):
    if value is None or value == "":
        return None
    return 1 if str(value).strip().lower() in {"1", "true", "yes", "y"} else 0


def _build_protacability_filters(args):
    filters = {}
    filters["virus_names"] = _normalize_multi_values(args.getlist("virus_name"))
    filters["protein_types"] = _normalize_multi_values(args.getlist("protein_type"))
    filters["tiers"] = _normalize_multi_values(args.getlist("tier"))
    filters["warhead_tiers"] = _normalize_multi_values(args.getlist("warhead_tier"))
    filters["readiness_tiers"] = _normalize_multi_values(args.getlist("readiness_tier"))
    filters["evidence_levels"] = _normalize_multi_values(args.getlist("evidence_level"))
    filters["smiles_sources"] = _normalize_multi_values(args.getlist("smiles_source"))
    filters["ligand"] = (args.get("ligand") or "").strip()
    filters["ligand_presence"] = (args.get("ligand_presence") or "").strip()
    filters["pdb_code"] = (args.get("pdb_code") or "").strip()
    filters["ligand_context_class"] = (args.get("ligand_context_class") or "").strip()
    filters["min_score"] = args.get("min_score", type=float)
    filters["min_warhead_score"] = args.get("min_warhead_score", type=float)
    filters["min_readiness_score"] = args.get("min_readiness_score", type=float)
    filters["has_candidate_ligand"] = _parse_bool_filter(args.get("has_candidate_ligand"))
    filters["has_exposed_lysine"] = _parse_bool_filter(args.get("has_exposed_lysine"))
    filters["has_ligand_proximal_exposed_lysine"] = _parse_bool_filter(args.get("has_ligand_proximal_exposed_lysine"))
    filters["has_candidate_linker_atoms"] = _parse_bool_filter(args.get("has_candidate_linker_atoms"))
    filters["has_solvent_exposed_ligand_atoms"] = _parse_bool_filter(args.get("has_solvent_exposed_ligand_atoms"))
    filters["has_mapped_atoms"] = _parse_bool_filter(args.get("has_mapped_atoms"))
    filters["has_valid_rdkit_smiles"] = _parse_bool_filter(args.get("has_valid_rdkit_smiles"))
    filters["has_exposed_target_lysines"] = _parse_bool_filter(args.get("has_exposed_target_lysines"))
    return filters


def _copy_protacability_filters(filters):
    return {
        "virus_names": list(filters.get("virus_names") or []),
        "protein_types": list(filters.get("protein_types") or []),
        "tiers": list(filters.get("tiers") or []),
        "warhead_tiers": list(filters.get("warhead_tiers") or []),
        "readiness_tiers": list(filters.get("readiness_tiers") or []),
        "evidence_levels": list(filters.get("evidence_levels") or []),
        "smiles_sources": list(filters.get("smiles_sources") or []),
        "ligand": filters.get("ligand") or "",
        "ligand_presence": filters.get("ligand_presence") or "",
        "pdb_code": filters.get("pdb_code") or "",
        "ligand_context_class": filters.get("ligand_context_class") or "",
        "min_score": filters.get("min_score"),
        "min_warhead_score": filters.get("min_warhead_score"),
        "min_readiness_score": filters.get("min_readiness_score"),
        "has_candidate_ligand": filters.get("has_candidate_ligand"),
        "has_exposed_lysine": filters.get("has_exposed_lysine"),
        "has_ligand_proximal_exposed_lysine": filters.get("has_ligand_proximal_exposed_lysine"),
        "has_candidate_linker_atoms": filters.get("has_candidate_linker_atoms"),
        "has_solvent_exposed_ligand_atoms": filters.get("has_solvent_exposed_ligand_atoms"),
        "has_mapped_atoms": filters.get("has_mapped_atoms"),
        "has_valid_rdkit_smiles": filters.get("has_valid_rdkit_smiles"),
        "has_exposed_target_lysines": filters.get("has_exposed_target_lysines"),
    }


def _clear_protacability_filter_dimension(filters, key):
    if key == "virus_names":
        filters["virus_names"] = []
    elif key == "protein_types":
        filters["protein_types"] = []
    elif key == "tiers":
        filters["tiers"] = []
    elif key == "warhead_tiers":
        filters["warhead_tiers"] = []
    elif key == "readiness_tiers":
        filters["readiness_tiers"] = []
    elif key == "evidence_levels":
        filters["evidence_levels"] = []
    elif key == "smiles_sources":
        filters["smiles_sources"] = []
    elif key in {"ligand", "ligand_presence", "pdb_code", "ligand_context_class"}:
        filters[key] = ""
    elif key in {
        "min_score",
        "min_warhead_score",
        "min_readiness_score",
        "has_candidate_ligand",
        "has_exposed_lysine",
        "has_ligand_proximal_exposed_lysine",
        "has_candidate_linker_atoms",
        "has_solvent_exposed_ligand_atoms",
        "has_mapped_atoms",
        "has_valid_rdkit_smiles",
        "has_exposed_target_lysines",
    }:
        filters[key] = None


def _filter_options_for_context(rows, filters, collapse_labels=True, ignore_key=None):
    scoped_filters = _copy_protacability_filters(filters)
    if ignore_key:
        _clear_protacability_filter_dimension(scoped_filters, ignore_key)

    filtered_rows = _filter_protacability_rows(rows, scoped_filters, collapse_labels=collapse_labels)
    if ignore_key != "ligand_context_class":
        filtered_rows = _apply_ligand_context_filter(filtered_rows, scoped_filters.get("ligand_context_class"))
    return filtered_rows


def _build_protacability_filter_options_payload(conn, args):
    collapse_labels = _protacability_collapse_labels(args.get("collapse_labels"))
    filters = _build_protacability_filters(args)
    readiness_rows, warhead_rows = _load_protacability_enrichment_tables(conn)
    rows = _decorate_protacability_rows(
        _load_protacability_assessment_rows(conn),
        collapse_labels=collapse_labels,
        readiness_rows=readiness_rows,
        warhead_rows=warhead_rows,
    )

    virus_rows = _filter_options_for_context(rows, filters, collapse_labels=collapse_labels, ignore_key="virus_names")
    protein_rows = _filter_options_for_context(rows, filters, collapse_labels=collapse_labels, ignore_key="protein_types")
    tier_rows = _filter_options_for_context(rows, filters, collapse_labels=collapse_labels, ignore_key="tiers")
    ligand_rows = _filter_options_for_context(rows, filters, collapse_labels=collapse_labels, ignore_key="ligand")
    context_rows = _filter_options_for_context(rows, filters, collapse_labels=collapse_labels, ignore_key="ligand_context_class")

    protein_field = "display_protein_type" if collapse_labels else "protein_type"
    ligand_context_classes = [
        {"value": "candidate_small_molecule", "label": "Candidate small-molecule context"},
        {"value": "candidate_plus_glycan", "label": "Candidate + glycan context"},
        {"value": "glycan_only", "label": "Glycan-only context"},
        {"value": "glycan_common_mixed", "label": "Glycan/common structural context"},
        {"value": "common_buffer_only", "label": "Common buffer/crystal context"},
        {"value": "no_ligand_context", "label": "No ligand context"},
    ]
    available_contexts = {row.get("ligand_context_class") for row in context_rows if row.get("ligand_context_class")}

    return {
        "data_available": True,
        "virus_names": sorted({row.get("virus_name") for row in virus_rows if row.get("virus_name")}),
        "protein_types": sorted({row.get(protein_field) for row in protein_rows if row.get(protein_field)}),
        "tiers": sorted(
            {row.get("protacability_tier") for row in tier_rows if row.get("protacability_tier")},
            key=lambda tier: (_tier_rank(tier), tier),
        ),
        "warhead_tiers": sorted({row.get("warhead_linkability_tier") for row in rows if row.get("warhead_linkability_tier")}),
        "readiness_tiers": sorted({row.get("degrader_design_readiness_tier") for row in rows if row.get("degrader_design_readiness_tier")}),
        "evidence_levels": sorted({row.get("evidence_level") for row in rows if row.get("evidence_level")}),
        "smiles_sources": sorted({row.get("smiles_source") for row in rows if row.get("smiles_source")}),
        "ligands": sorted({ligand for row in ligand_rows for ligand in (row.get("ligand_names") or []) if ligand}),
        "ligand_context_classes": [item for item in ligand_context_classes if item["value"] in available_contexts],
    }


def _protacability_where_clause(filters):
    clauses = []
    params = []

    if filters["virus_names"]:
        placeholders = ", ".join(["?"] * len(filters["virus_names"]))
        clauses.append(f"virus_name IN ({placeholders})")
        params.extend(filters["virus_names"])

    if filters["protein_types"]:
        placeholders = ", ".join(["?"] * len(filters["protein_types"]))
        clauses.append(f"protein_type IN ({placeholders})")
        params.extend(filters["protein_types"])

    if filters["tiers"]:
        placeholders = ", ".join(["?"] * len(filters["tiers"]))
        clauses.append(f"protacability_tier IN ({placeholders})")
        params.extend(filters["tiers"])

    if filters["ligand"]:
        clauses.append("candidate_ligand_resnames LIKE ?")
        params.append(f"%{filters['ligand']}%")

    if filters["pdb_code"]:
        clauses.append("pdb_code LIKE ?")
        params.append(f"%{filters['pdb_code']}%")

    if filters["min_score"] is not None:
        clauses.append("protacability_proxy_score >= ?")
        params.append(filters["min_score"])

    for key in (
        "has_candidate_ligand",
        "has_exposed_lysine",
        "has_ligand_proximal_exposed_lysine",
    ):
        if filters[key] is not None:
            clauses.append(f"{key} = ?")
            params.append(filters[key])

    where_sql = ""
    if clauses:
        where_sql = "WHERE " + " AND ".join(clauses)
    return where_sql, params


def _protacability_sort_clause(sort_value):
    sort_raw = (sort_value or "protacability_proxy_score_desc").strip().lower()
    direction = "DESC"
    column_key = sort_raw
    if sort_raw.endswith("_asc"):
        column_key = sort_raw[:-4]
        direction = "ASC"
    elif sort_raw.endswith("_desc"):
        column_key = sort_raw[:-5]
        direction = "DESC"

    column = PROTACABILITY_SORT_COLUMNS.get(column_key, "protacability_proxy_score")
    return f"{column} {direction}"


def _protacability_view_mode(value):
    view = (value or "targets").strip().lower()
    return view if view in PROTACABILITY_VIEW_MODES else "targets"


def _protacability_collapse_labels(value):
    if value is None:
        return True
    return str(value).strip().lower() not in {"0", "false", "no", "off"}


def _protacability_sort_token(view, sort_value):
    options = PROTACABILITY_SORT_OPTIONS[view]
    if sort_value in options:
        return sort_value
    return "readiness_priority_desc"


def _split_candidate_ligands(value):
    if value is None:
        return []
    if isinstance(value, list):
        raw_parts = value
    else:
        raw_parts = re.split(r"[;,|]", str(value))
    ligands = []
    seen = set()
    for raw in raw_parts:
        name = str(raw).strip().upper()
        if name and name not in seen:
            ligands.append(name)
            seen.add(name)
    return ligands


def classify_ligand_context(ligand_resnames):
    ligands = _split_candidate_ligands(ligand_resnames)
    if not ligands:
        return {
            "ligand_context_class": "no_ligand_context",
            "ligand_context_label": "No ligand context",
            "informative_ligand_count": 0,
            "glycan_ligand_count": 0,
            "common_ligand_count": 0,
            "candidate_ligand_count_display": 0,
            "glycan_ligands": [],
            "common_ligands": [],
            "candidate_ligands": [],
            "ligand_priority": 0,
        }

    glycan_ligands = sorted([ligand for ligand in ligands if ligand in PROTACABILITY_GLYCAN_LIGANDS])
    common_ligands = sorted([ligand for ligand in ligands if ligand in PROTACABILITY_COMMON_CONTEXT_LIGANDS])
    candidate_ligands = sorted([
        ligand for ligand in ligands
        if ligand not in PROTACABILITY_GLYCAN_LIGANDS and ligand not in PROTACABILITY_COMMON_CONTEXT_LIGANDS
    ])

    glycan_count = len(glycan_ligands)
    common_count = len(common_ligands)
    candidate_count = len(candidate_ligands)

    if candidate_count > 0 and glycan_count > 0:
        context_class = "candidate_plus_glycan"
        context_label = "Candidate + glycan context"
        priority = 4
    elif candidate_count > 0:
        context_class = "candidate_small_molecule"
        context_label = "Candidate small-molecule context"
        priority = 5
    elif glycan_count > 0 and common_count > 0:
        context_class = "glycan_common_mixed"
        context_label = "Glycan/common structural context"
        priority = 2
    elif glycan_count > 0:
        context_class = "glycan_only"
        context_label = "Glycan-only context"
        priority = 1
    elif common_count > 0:
        context_class = "common_buffer_only"
        context_label = "Common buffer/crystal context"
        priority = 0
    else:
        context_class = "no_ligand_context"
        context_label = "No ligand context"
        priority = 0

    return {
        "ligand_context_class": context_class,
        "ligand_context_label": context_label,
        "informative_ligand_count": candidate_count if candidate_count > 0 else 0,
        "glycan_ligand_count": glycan_count,
        "common_ligand_count": common_count,
        "candidate_ligand_count_display": candidate_count,
        "glycan_ligands": glycan_ligands,
        "common_ligands": common_ligands,
        "candidate_ligands": candidate_ligands,
        "ligand_priority": priority,
    }


def _display_protein_priority(label):
    label_normalized = (label or "").strip().lower()
    priority_map = {
        "protease": 0,
        "reverse transcriptase": 1,
        "integrase": 2,
        "rnase h": 3,
        "spike protein": 4,
        "envelope protein/glycoprotein": 5,
        "polymerase": 10,
        "helicase": 11,
        "pol protein": 80,
        "nsp proteins": 90,
        "fusion protein": 91,
        "transmembrane protein": 92,
    }
    return priority_map.get(label_normalized, 40)


def _choose_display_protein_type(labels):
    cleaned = [label for label in labels if str(label).strip()]
    if not cleaned:
        return ""
    unique_labels = sorted(set(cleaned), key=lambda label: (_display_protein_priority(label), len(label), label))
    return unique_labels[0]


def _numeric_value(value, default=0.0):
    try:
        if value is None or value == "":
            return default
        return float(value)
    except (TypeError, ValueError):
        return default


def _tier_rank(tier):
    return PROTACABILITY_TIER_ORDER.get(tier or "", 9)


def _normalize_text_key(value):
    return str(value or "").strip().upper()


def _normalize_residue_key(value):
    text = str(value or "").strip()
    if not text:
        return ""
    return text.upper()


def _boolish(value):
    if value in (None, "", "nan"):
        return 0
    if isinstance(value, bool):
        return int(value)
    text = str(value).strip().lower()
    return 1 if text in {"1", "true", "yes", "y"} else 0


def _has_positive_value(value):
    return _numeric_value(value, 0) > 0


def _safe_int(value):
    try:
        if value in (None, "", "nan"):
            return None
        return int(float(value))
    except (TypeError, ValueError):
        return None


def _candidate_linker_ids_list(value):
    if value in (None, ""):
        return []
    return [part.strip() for part in str(value).split(";") if part.strip()]


def _warhead_evidence_rank(row):
    return (
        _has_positive_value(row.get("candidate_linker_atom_count")),
        _has_positive_value(row.get("solvent_exposed_mapped_atom_count")),
        _has_positive_value(row.get("pdb_to_smiles_mapped_atom_count")),
        _has_positive_value(row.get("solvent_exposed_ligand_atom_count")),
        _boolish(row.get("rdkit_valid_smiles")),
        _boolish(row.get("smiles_available")),
        _numeric_value(row.get("interaction_preservation_score")),
        _numeric_value(row.get("warhead_linkability_score")),
    )


def _readiness_evidence_rank(row):
    return (
        _numeric_value(row.get("degrader_design_readiness_score")),
        _numeric_value(row.get("warhead_linkability_score")),
        _numeric_value(row.get("target_lysine_accessibility_score")),
        _numeric_value(row.get("protein_structural_priority_score")),
        _numeric_value(row.get("ternary_geometry_cue_score")),
    )


def _choose_best_record(rows, rank_fn):
    rows = list(rows or [])
    if not rows:
        return None
    return max(rows, key=rank_fn)


def _load_optional_table_rows(conn, table_name):
    if not _table_exists(conn, table_name):
        return []
    return [dict(row) for row in conn.execute(f"SELECT * FROM {table_name}").fetchall()]


def _build_warhead_indexes(rows):
    exact = defaultdict(list)
    by_resname = defaultdict(list)
    for row in rows:
        exact_key = (
            _normalize_text_key(row.get("pdb_code")),
            _normalize_text_key(row.get("ligand_resname")),
            _normalize_text_key(row.get("ligand_chain")),
            _normalize_residue_key(row.get("ligand_residue_id")),
        )
        if exact_key[0] and exact_key[1]:
            exact[exact_key].append(row)
            by_resname[(exact_key[0], exact_key[1])].append(row)
    return {
        "exact": {key: _choose_best_record(value, _warhead_evidence_rank) for key, value in exact.items()},
        "by_resname": {key: _choose_best_record(value, _warhead_evidence_rank) for key, value in by_resname.items()},
    }


def _build_readiness_indexes(rows):
    by_full = defaultdict(list)
    by_chain = defaultdict(list)
    by_target = defaultdict(list)
    for row in rows:
        full_key = (
            _normalize_text_key(row.get("virus_name")),
            _normalize_text_key(row.get("protein_type")),
            _normalize_text_key(row.get("pdb_code")),
            _normalize_text_key(row.get("chain_id")),
        )
        chain_key = (full_key[2], full_key[3])
        target_key = (full_key[0], full_key[1], full_key[2])
        if full_key[2] and full_key[3]:
            by_full[full_key].append(row)
            by_chain[chain_key].append(row)
        if full_key[2]:
            by_target[target_key].append(row)
    return {
        "by_full": {key: _choose_best_record(value, _readiness_evidence_rank) for key, value in by_full.items()},
        "by_chain": {key: _choose_best_record(value, _readiness_evidence_rank) for key, value in by_chain.items()},
        "by_target": {key: _choose_best_record(value, _readiness_evidence_rank) for key, value in by_target.items()},
    }


def _row_representative_ligands(row):
    ligands = list(row.get("ligand_names") or [])
    best_resname = _normalize_text_key(row.get("best_ligand_resname"))
    if best_resname and best_resname not in ligands:
        ligands.insert(0, best_resname)
    return [ligand for ligand in ligands if ligand]


def _find_matching_warhead_row(row, warhead_indexes):
    if not warhead_indexes:
        return None
    pdb_code = _normalize_text_key(row.get("pdb_code"))
    ligand_candidates = _row_representative_ligands(row)
    best_resname = _normalize_text_key(row.get("best_ligand_resname"))
    best_chain = _normalize_text_key(row.get("best_ligand_chain"))
    best_residue_id = _normalize_residue_key(row.get("best_ligand_residue_id"))
    if pdb_code and best_resname:
        exact = warhead_indexes["exact"].get((pdb_code, best_resname, best_chain, best_residue_id))
        if exact:
            return exact
        fallback = warhead_indexes["by_resname"].get((pdb_code, best_resname))
        if fallback:
            return fallback
    for ligand in ligand_candidates:
        fallback = warhead_indexes["by_resname"].get((pdb_code, _normalize_text_key(ligand)))
        if fallback:
            return fallback
    return None


def _find_matching_readiness_row(row, readiness_indexes):
    if not readiness_indexes:
        return None
    full_key = (
        _normalize_text_key(row.get("virus_name")),
        _normalize_text_key(row.get("protein_type")),
        _normalize_text_key(row.get("pdb_code")),
        _normalize_text_key(row.get("chain_id")),
    )
    match = readiness_indexes["by_full"].get(full_key)
    if match:
        return match
    chain_match = readiness_indexes["by_chain"].get((full_key[2], full_key[3]))
    if chain_match:
        return chain_match
    return readiness_indexes["by_target"].get((full_key[0], full_key[1], full_key[2]))


def _merge_optional_protacability_data(rows, readiness_rows=None, warhead_rows=None):
    readiness_indexes = _build_readiness_indexes(readiness_rows or []) if readiness_rows else None
    warhead_indexes = _build_warhead_indexes(warhead_rows or []) if warhead_rows else None

    merged = []
    for row in rows:
        current = dict(row)
        readiness_match = _find_matching_readiness_row(current, readiness_indexes)
        if readiness_match:
            for key, value in readiness_match.items():
                if key not in current or current.get(key) in (None, "", []):
                    current[key] = value
                elif key in {
                    "protein_structural_priority_score",
                    "warhead_linkability_score",
                    "target_lysine_accessibility_score",
                    "ternary_geometry_cue_score",
                    "degrader_design_readiness_score",
                    "degrader_design_readiness_tier",
                    "evidence_level",
                    "best_linker_geometry_class",
                    "short_linker_geometry_feasible",
                    "medium_linker_geometry_feasible",
                    "long_linker_geometry_feasible",
                    "min_lys_ligand_distance_a",
                    "near_ligand_exposed_lys_count",
                    "readiness_flags",
                    "readiness_notes",
                    "best_ligand_resname",
                    "best_ligand_chain",
                    "best_ligand_residue_id",
                }:
                    current[key] = value

        warhead_match = _find_matching_warhead_row(current, warhead_indexes)
        if warhead_match:
            for key, value in warhead_match.items():
                if key not in current or current.get(key) in (None, "", []):
                    current[key] = value
                elif key in {
                    "candidate_linker_atom_count",
                    "candidate_linker_atom_ids",
                    "solvent_exposed_ligand_atom_count",
                    "solvent_exposed_mapped_atom_count",
                    "pdb_to_smiles_mapped_atom_count",
                    "functional_group_count",
                    "functional_group_types",
                    "interaction_preservation_score",
                    "smiles_available",
                    "smiles_source",
                    "representative_smiles",
                    "rdkit_available",
                    "rdkit_valid_smiles",
                    "mol_weight",
                    "tpsa",
                    "hbd",
                    "hba",
                    "rotatable_bonds",
                    "meaningful_contact_count",
                    "strong_contact_count",
                    "contact_atom_count",
                    "strong_contact_atom_count",
                    "warhead_flags",
                    "warhead_notes",
                    "warhead_linkability_tier",
                    "warhead_linkability_label",
                }:
                    current[key] = value

        current["candidate_linker_atom_ids_list"] = _candidate_linker_ids_list(current.get("candidate_linker_atom_ids"))
        current["has_candidate_linker_atoms"] = int(_has_positive_value(current.get("candidate_linker_atom_count")))
        current["has_solvent_exposed_ligand_atoms"] = int(_has_positive_value(current.get("solvent_exposed_ligand_atom_count")))
        current["has_mapped_atoms"] = int(_has_positive_value(current.get("pdb_to_smiles_mapped_atom_count")))
        current["has_valid_rdkit_smiles"] = int(_boolish(current.get("rdkit_valid_smiles")))
        current["has_exposed_target_lysines"] = int(_has_positive_value(current.get("exposed_lys_count")))
        current["warhead_evidence_score"] = (
            50 * current["has_candidate_linker_atoms"]
            + 25 * int(_has_positive_value(current.get("solvent_exposed_mapped_atom_count")))
            + 15 * current["has_mapped_atoms"]
            + 10 * current["has_valid_rdkit_smiles"]
        )
        current["readiness_rank_score"] = (
            _numeric_value(current.get("degrader_design_readiness_score"))
            + current["warhead_evidence_score"] / 10.0
        )
        merged.append(current)
    return merged


def _ligand_context_rank(context_class):
    rank_map = {
        "candidate_small_molecule": 5,
        "candidate_plus_glycan": 4,
        "glycan_common_mixed": 2,
        "glycan_only": 1,
        "common_buffer_only": 0,
        "no_ligand_context": 0,
    }
    return rank_map.get(context_class or "no_ligand_context", 0)


def _row_has_any_ligand_context(row):
    ligand_context_class = row.get("ligand_context_class") or "no_ligand_context"
    candidate_names = (
        row.get("candidate_ligand_resnames")
        or row.get("candidate_ligand_resnames_display")
        or row.get("candidate_ligand_resnames_full")
        or row.get("distinct_ligands")
        or row.get("distinct_ligands_full")
        or row.get("top_candidate_ligands_full")
        or ""
    )
    informative_count = (
        row.get("informative_ligand_count")
        or row.get("distinct_ligand_count")
        or row.get("ligand_context_count")
        or row.get("candidate_ligand_count")
        or 0
    )
    return (
        ligand_context_class not in ("", "no_ligand_context", None)
        or bool(str(candidate_names).strip())
        or int(_numeric_value(informative_count, 0)) > 0
    )


def _row_has_candidate_small_molecule(row):
    ligand_context_class = row.get("ligand_context_class") or row.get("best_ligand_context_class") or ""
    return ligand_context_class in ("candidate_small_molecule", "candidate_plus_glycan")


def _apply_ligand_presence_filter(rows, ligand_presence):
    target = (ligand_presence or "").strip()
    if not target:
        return rows

    filtered = []
    for row in rows:
        has_any = _row_has_any_ligand_context(row)
        has_candidate = _row_has_candidate_small_molecule(row)
        if target == "has_any_ligand" and not has_any:
            continue
        if target == "no_ligand" and has_any:
            continue
        if target == "has_candidate_ligand" and not has_candidate:
            continue
        if target == "no_candidate_ligand" and has_candidate:
            continue
        filtered.append(row)
    return filtered


def _row_priority_key(row):
    return (
        _numeric_value(row.get("degrader_design_readiness_score")),
        _numeric_value(row.get("warhead_evidence_score")),
        _numeric_value(row.get("warhead_linkability_score")),
        _numeric_value(row.get("target_lysine_accessibility_score")),
        _numeric_value(row.get("protein_structural_priority_score") or row.get("protacability_proxy_score")),
        _numeric_value(row.get("ligand_priority")),
        _numeric_value(row.get("informative_ligand_count")),
        _numeric_value(row.get("exposed_lys_fraction")),
        _numeric_value(row.get("candidate_ligand_count_display") or row.get("candidate_ligand_count")),
        _numeric_value(row.get("lysine_surface_fraction")),
        _numeric_value(row.get("chain_length_aa")),
        -_tier_rank(row.get("protacability_tier")),
    )


def _decorate_protacability_rows(rows, collapse_labels=True, readiness_rows=None, warhead_rows=None):
    decorated = _merge_optional_protacability_data(rows, readiness_rows=readiness_rows, warhead_rows=warhead_rows)
    display_labels = {}
    if collapse_labels:
        labels_by_chain = defaultdict(list)
        for row in decorated:
            key = (row.get("virus_name"), row.get("pdb_code"), row.get("chain_id"))
            labels_by_chain[key].append(row.get("protein_type") or "")
        display_labels = {
            key: _choose_display_protein_type(values)
            for key, values in labels_by_chain.items()
        }

    for row in decorated:
        key = (row.get("virus_name"), row.get("pdb_code"), row.get("chain_id"))
        ligand_names = _split_candidate_ligands(row.get("candidate_ligand_resnames"))
        context = classify_ligand_context(ligand_names)
        display_protein_type = display_labels.get(key, row.get("protein_type")) if collapse_labels else row.get("protein_type")
        row["display_protein_type"] = display_protein_type
        row["ligand_names"] = ligand_names
        row["candidate_ligand_resnames_full"] = "; ".join(ligand_names) if ligand_names else ""
        row["candidate_ligand_resnames_display"] = row["candidate_ligand_resnames_full"] or "None"
        row.update(context)
        row["has_small_molecule_context"] = 1 if row["ligand_context_class"] in {"candidate_small_molecule", "candidate_plus_glycan"} else 0
        row["has_glycan_context"] = 1 if row["glycan_ligand_count"] > 0 else 0
        row["has_common_context"] = 1 if row["common_ligand_count"] > 0 else 0
    return decorated


def _filter_protacability_rows(rows, filters, collapse_labels=True):
    protein_field = "display_protein_type" if collapse_labels else "protein_type"
    virus_filter = set(filters["virus_names"])
    protein_filter = set(filters["protein_types"])
    tier_filter = set(filters["tiers"])
    warhead_tier_filter = set(filters["warhead_tiers"])
    readiness_tier_filter = set(filters["readiness_tiers"])
    evidence_level_filter = set(filters["evidence_levels"])
    smiles_source_filter = set(filters["smiles_sources"])
    ligand_query = (filters["ligand"] or "").strip().upper()
    pdb_query = (filters["pdb_code"] or "").strip().upper()

    filtered = []
    for row in rows:
        if virus_filter and row.get("virus_name") not in virus_filter:
            continue
        if protein_filter and row.get(protein_field) not in protein_filter:
            continue
        if tier_filter and row.get("protacability_tier") not in tier_filter:
            continue
        if warhead_tier_filter and row.get("warhead_linkability_tier") not in warhead_tier_filter:
            continue
        if readiness_tier_filter and row.get("degrader_design_readiness_tier") not in readiness_tier_filter:
            continue
        if evidence_level_filter and row.get("evidence_level") not in evidence_level_filter:
            continue
        if smiles_source_filter and (row.get("smiles_source") or "") not in smiles_source_filter:
            continue
        if ligand_query:
            ligand_names = row.get("ligand_names") or []
            ligand_blob = ";".join(ligand_names)
            if ligand_query not in ligand_blob and ligand_query not in (row.get("candidate_ligand_resnames_display") or "").upper():
                continue
        ligand_presence = filters.get("ligand_presence") or ""
        if ligand_presence == "has_any_ligand" and not _row_has_any_ligand_context(row):
            continue
        if ligand_presence == "no_ligand" and _row_has_any_ligand_context(row):
            continue
        if ligand_presence == "has_candidate_ligand" and not _row_has_candidate_small_molecule(row):
            continue
        if ligand_presence == "no_candidate_ligand" and _row_has_candidate_small_molecule(row):
            continue
        if pdb_query and pdb_query not in str(row.get("pdb_code") or "").upper():
            continue
        if filters["min_score"] is not None and _numeric_value(row.get("protacability_proxy_score"), -1) < filters["min_score"]:
            continue
        if filters["min_warhead_score"] is not None and _numeric_value(row.get("warhead_linkability_score"), -1) < filters["min_warhead_score"]:
            continue
        if filters["min_readiness_score"] is not None and _numeric_value(row.get("degrader_design_readiness_score"), -1) < filters["min_readiness_score"]:
            continue
        if filters["has_candidate_ligand"] is not None and int(bool(row.get("has_candidate_ligand"))) != filters["has_candidate_ligand"]:
            continue
        if filters["has_exposed_lysine"] is not None and int(bool(row.get("has_exposed_lysine"))) != filters["has_exposed_lysine"]:
            continue
        if filters["has_ligand_proximal_exposed_lysine"] is not None and int(bool(row.get("has_ligand_proximal_exposed_lysine"))) != filters["has_ligand_proximal_exposed_lysine"]:
            continue
        if filters["has_candidate_linker_atoms"] is not None and int(bool(row.get("has_candidate_linker_atoms"))) != filters["has_candidate_linker_atoms"]:
            continue
        if filters["has_solvent_exposed_ligand_atoms"] is not None and int(bool(row.get("has_solvent_exposed_ligand_atoms"))) != filters["has_solvent_exposed_ligand_atoms"]:
            continue
        if filters["has_mapped_atoms"] is not None and int(bool(row.get("has_mapped_atoms"))) != filters["has_mapped_atoms"]:
            continue
        if filters["has_valid_rdkit_smiles"] is not None and int(bool(row.get("has_valid_rdkit_smiles"))) != filters["has_valid_rdkit_smiles"]:
            continue
        if filters["has_exposed_target_lysines"] is not None and int(bool(row.get("has_exposed_target_lysines"))) != filters["has_exposed_target_lysines"]:
            continue
        filtered.append(row)
    return filtered


def _apply_ligand_context_filter(rows, ligand_context_class):
    target = (ligand_context_class or "").strip()
    alias_map = {
        "mixed_candidate_and_glycan": "candidate_plus_glycan",
        "glycan_mixed": "glycan_common_mixed",
        "none": "no_ligand_context",
    }
    target = alias_map.get(target, target)
    if not target:
        return rows
    return [row for row in rows if row.get("ligand_context_class") == target]


def _dedupe_display_chain_rows(rows):
    grouped = defaultdict(list)
    for row in rows:
        key = (
            row.get("virus_name"),
            row.get("pdb_code"),
            row.get("chain_id"),
            row.get("display_protein_type") or row.get("protein_type"),
        )
        grouped[key].append(row)
    deduped = []
    for group_rows in grouped.values():
        deduped.append(max(group_rows, key=_row_priority_key))
    return deduped


def _combine_ligands(rows):
    ligands = []
    seen = set()
    for row in rows:
        for ligand in row.get("ligand_names") or []:
            if ligand not in seen:
                ligands.append(ligand)
                seen.add(ligand)
    return ligands


def _ligand_preview(ligands, max_items=8):
    ligands = list(ligands or [])
    if not ligands:
        return "None"
    preview_items = ligands[:max_items]
    preview = "; ".join(preview_items)
    remaining = len(ligands) - len(preview_items)
    if remaining > 0:
        preview = f"{preview} + {remaining} more"
    return preview


def _round_or_none(value, digits=2):
    if value is None:
        return None
    return round(float(value), digits)


def _protacability_enrichment_snapshot(row):
    fields = [
        "protein_structural_priority_score",
        "warhead_linkability_score",
        "warhead_linkability_tier",
        "warhead_linkability_label",
        "target_lysine_accessibility_score",
        "ternary_geometry_cue_score",
        "degrader_design_readiness_score",
        "degrader_design_readiness_tier",
        "evidence_level",
        "readiness_flags",
        "readiness_notes",
        "best_ligand_resname",
        "best_ligand_chain",
        "best_ligand_residue_id",
        "candidate_linker_atom_count",
        "candidate_linker_atom_ids",
        "candidate_linker_atom_ids_list",
        "solvent_exposed_ligand_atom_count",
        "solvent_exposed_mapped_atom_count",
        "pdb_to_smiles_mapped_atom_count",
        "functional_group_count",
        "interaction_preservation_score",
        "smiles_available",
        "smiles_source",
        "rdkit_valid_smiles",
        "warhead_flags",
        "warhead_notes",
        "best_linker_geometry_class",
        "short_linker_geometry_feasible",
        "medium_linker_geometry_feasible",
        "long_linker_geometry_feasible",
        "min_lys_ligand_distance_a",
        "near_ligand_exposed_lys_count",
        "has_candidate_linker_atoms",
        "has_solvent_exposed_ligand_atoms",
        "has_mapped_atoms",
        "has_valid_rdkit_smiles",
        "has_exposed_target_lysines",
        "warhead_evidence_score",
        "readiness_rank_score",
    ]
    return {field: row.get(field) for field in fields}


def _group_structure_rows(rows):
    chain_rows = _dedupe_display_chain_rows(rows)
    grouped = defaultdict(list)
    for row in chain_rows:
        key = (row.get("virus_name"), row.get("display_protein_type") or row.get("protein_type"), row.get("pdb_code"))
        grouped[key].append(row)

    grouped_rows = []
    for (virus_name, protein_type, pdb_code), group_rows in grouped.items():
        representative = max(group_rows, key=_row_priority_key)
        ligands = _combine_ligands(group_rows)
        context = classify_ligand_context(ligands)
        scores = [_numeric_value(row.get("protacability_proxy_score")) for row in group_rows]
        p_is = [_numeric_value(row.get("isoelectric_point")) for row in group_rows if row.get("isoelectric_point") not in (None, "")]
        exposed_fracs = [_numeric_value(row.get("exposed_lys_fraction")) for row in group_rows]
        lys_surface_fracs = [_numeric_value(row.get("lysine_surface_fraction")) for row in group_rows]
        grouped_rows.append({
            "view_type": "summary",
            "virus_name": virus_name,
            "protein_type": protein_type,
            "pdb_code": pdb_code,
            "representative_chain_id": representative.get("chain_id"),
            "chain_count": len({row.get("chain_id") for row in group_rows}),
            "best_score": _round_or_none(max(scores), 2),
            "avg_score": _round_or_none(sum(scores) / len(scores), 2),
            "best_tier": representative.get("protacability_tier"),
            "best_pI": _round_or_none(_numeric_value(representative.get("isoelectric_point"), None), 2) if representative.get("isoelectric_point") not in (None, "") else None,
            "min_pI": _round_or_none(min(p_is), 2) if p_is else None,
            "max_pI": _round_or_none(max(p_is), 2) if p_is else None,
            "total_lys_count": int(sum(_numeric_value(row.get("lys_count")) for row in group_rows)),
            "max_chain_lys_count": int(max(_numeric_value(row.get("lys_count")) for row in group_rows)) if group_rows else 0,
            "total_exposed_lys_count": int(sum(_numeric_value(row.get("exposed_lys_count")) for row in group_rows)),
            "best_exposed_lys_fraction": _round_or_none(max(exposed_fracs), 3) if exposed_fracs else None,
            "avg_exposed_lys_fraction": _round_or_none(sum(exposed_fracs) / len(exposed_fracs), 3) if exposed_fracs else None,
            "max_lysine_surface_fraction": _round_or_none(max(lys_surface_fracs), 3) if lys_surface_fracs else None,
            "has_candidate_ligand": int(any(int(bool(row.get("has_candidate_ligand"))) for row in group_rows)),
            "candidate_ligand_resnames_full": "; ".join(ligands) if ligands else "",
            "candidate_ligand_resnames": _ligand_preview(ligands, max_items=8),
            "candidate_ligand_count": len(ligands),
            "has_exposed_lysine": int(any(int(bool(row.get("has_exposed_lysine"))) for row in group_rows)),
            "has_ligand_proximal_exposed_lysine": int(any(int(bool(row.get("has_ligand_proximal_exposed_lysine"))) for row in group_rows)),
            "best_annotation": representative.get("linker_docking_site_annotation") or "—",
            "grouped_notes": f"Grouped structure summary across {len(group_rows)} chain-level rows.",
            "raw_chain_rows_count": len(group_rows),
            "distinct_ligand_count": len(ligands),
            **context,
            **_protacability_enrichment_snapshot(representative),
        })
    return grouped_rows


def _group_protein_rows(rows):
    structure_rows = _group_structure_rows(rows)
    grouped = defaultdict(list)
    for row in structure_rows:
        key = (row.get("virus_name"), row.get("protein_type"))
        grouped[key].append(row)

    grouped_rows = []
    for (virus_name, protein_type), group_rows in grouped.items():
        representative = max(
            group_rows,
            key=lambda row: (
                _numeric_value(row.get("ligand_priority")),
                _numeric_value(row.get("best_score")),
                _numeric_value(row.get("informative_ligand_count")),
                _numeric_value(row.get("best_exposed_lys_fraction")),
                _numeric_value(row.get("candidate_ligand_count")),
            ),
        )
        ligands = _combine_ligands([{ "ligand_names": _split_candidate_ligands(row.get("candidate_ligand_resnames")) } for row in group_rows])
        context = classify_ligand_context(ligands)
        grouped_rows.append({
            "view_type": "protein",
            "virus_name": virus_name,
            "protein_type": protein_type,
            "pdb_count": len({row.get("pdb_code") for row in group_rows}),
            "chain_count": int(sum(_numeric_value(row.get("chain_count")) for row in group_rows)),
            "high_priority_count": sum(1 for row in group_rows if row.get("best_tier") == "High structural priority"),
            "moderate_priority_count": sum(1 for row in group_rows if row.get("best_tier") == "Moderate structural priority"),
            "low_priority_count": sum(1 for row in group_rows if row.get("best_tier") == "Low structural priority"),
            "best_score": representative.get("best_score"),
            "avg_score": _round_or_none(sum(_numeric_value(row.get("avg_score")) for row in group_rows) / len(group_rows), 2),
            "top_pdb_code": representative.get("pdb_code"),
            "top_chain_id": representative.get("representative_chain_id"),
            "ligand_context_count": sum(1 for row in group_rows if row.get("ligand_context_class") in {"candidate_small_molecule", "candidate_plus_glycan"}),
            "glycan_only_context_count": sum(1 for row in group_rows if row.get("ligand_context_class") == "glycan_only"),
            "common_only_context_count": sum(1 for row in group_rows if row.get("ligand_context_class") == "common_buffer_only"),
            "glycan_context_count": sum(1 for row in group_rows if row.get("ligand_context_class") in {"glycan_only", "glycan_common_mixed", "candidate_plus_glycan"}),
            "exposed_lysine_context_count": sum(1 for row in group_rows if row.get("has_exposed_lysine")),
            "distinct_ligands_full": "; ".join(ligands) if ligands else "",
            "distinct_ligands": _ligand_preview(ligands, max_items=8),
            "distinct_ligand_count": len(ligands),
            "candidate_ligand_count": len(ligands),
            "best_annotation": representative.get("best_annotation") or "—",
            "best_tier": representative.get("best_tier"),
            **context,
            "structure_group_count": len(group_rows),
            **_protacability_enrichment_snapshot(representative),
        })
    return grouped_rows


def _target_interpretation(row):
    context_class = row.get("ligand_context_class")
    candidate_structures = int(_numeric_value(row.get("candidate_small_molecule_structure_count")))
    mixed_structures = int(_numeric_value(row.get("mixed_candidate_glycan_structure_count")))
    glycan_only = int(_numeric_value(row.get("glycan_only_structure_count")))
    exposed_fraction = _numeric_value(row.get("best_exposed_lys_fraction"))

    if candidate_structures > 0:
        return "Candidate ligand-supported target", "Contains candidate small-molecule contexts suitable for warhead triage."
    if mixed_structures > 0 or context_class == "candidate_plus_glycan":
        return "Candidate ligand-supported target", "Contains mixed candidate + glycan contexts that may still support ligand-first triage."
    if glycan_only > 0 and exposed_fraction >= 0.2:
        return "Surface-lysine-rich target", "Limited ligand evidence but strong exposed lysine signal for linker geometry hypotheses."
    if glycan_only > 0 or context_class in {"glycan_only", "glycan_common_mixed"}:
        return "Glycan/context-heavy target", "Dominated by glycan/common structural context; useful for surface biology but weaker drug-like evidence."
    return "Limited ligand evidence", "Current structures are mostly common-buffer or no-ligand contexts."


def _group_target_rows(rows):
    protein_rows = _group_protein_rows(rows)
    structure_rows = _group_structure_rows(rows)
    grouped = defaultdict(list)
    structures_by_target = defaultdict(list)
    for row in protein_rows:
        grouped[(row.get("virus_name"), row.get("protein_type"))].append(row)
    for srow in structure_rows:
        structures_by_target[(srow.get("virus_name"), srow.get("protein_type"))].append(srow)

    target_rows = []
    for (virus_name, protein_type), group_rows in grouped.items():
        representative = max(group_rows, key=lambda row: (_numeric_value(row.get("ligand_priority")), _numeric_value(row.get("best_score"))))
        target_structures = structures_by_target.get((virus_name, protein_type), [])
        distinct_ligands = _combine_ligands([{"ligand_names": _split_candidate_ligands(row.get("distinct_ligands_full"))} for row in group_rows])
        p_is = []
        exposed_fracs = []
        for srow in target_structures:
            if srow.get("min_pI") is not None:
                p_is.append(_numeric_value(srow.get("min_pI")))
            if srow.get("max_pI") is not None:
                p_is.append(_numeric_value(srow.get("max_pI")))
            if srow.get("best_exposed_lys_fraction") is not None:
                exposed_fracs.append(_numeric_value(srow.get("best_exposed_lys_fraction")))

        row = {
            "view_type": "targets",
            "virus_name": virus_name,
            "protein_type": protein_type,
            "target_key": f"{virus_name}::{protein_type}",
            "pdb_count": int(sum(_numeric_value(r.get("pdb_count")) for r in group_rows)),
            "chain_count": int(sum(_numeric_value(r.get("chain_count")) for r in group_rows)),
            "best_score": representative.get("best_score"),
            "avg_score": _round_or_none(sum(_numeric_value(r.get("avg_score")) for r in group_rows) / len(group_rows), 2),
            "best_tier": representative.get("best_tier"),
            "high_priority_count": int(sum(_numeric_value(r.get("high_priority_count")) for r in group_rows)),
            "candidate_small_molecule_structure_count": sum(1 for s in target_structures if s.get("ligand_context_class") == "candidate_small_molecule"),
            "mixed_candidate_glycan_structure_count": sum(1 for s in target_structures if s.get("ligand_context_class") == "candidate_plus_glycan"),
            "glycan_only_structure_count": sum(1 for s in target_structures if s.get("ligand_context_class") == "glycan_only"),
            "common_buffer_only_structure_count": sum(1 for s in target_structures if s.get("ligand_context_class") == "common_buffer_only"),
            "no_ligand_structure_count": max(0, int(sum(_numeric_value(r.get("pdb_count")) for r in group_rows)) - int(sum(_numeric_value(r.get("ligand_context_count")) for r in group_rows))),
            "best_pdb_code": representative.get("top_pdb_code"),
            "best_chain_id": representative.get("top_chain_id"),
            "best_ligand_context_class": representative.get("ligand_context_class"),
            "best_ligand_context_label": representative.get("ligand_context_label"),
            "top_candidate_ligands_preview": _ligand_preview(distinct_ligands, max_items=6),
            "top_candidate_ligands_full": "; ".join(distinct_ligands) if distinct_ligands else "",
            "distinct_candidate_ligand_count": len(distinct_ligands),
            "pI_min": _round_or_none(min(p_is), 2) if p_is else None,
            "pI_max": _round_or_none(max(p_is), 2) if p_is else None,
            "exposed_lys_total": int(sum(_numeric_value(r.get("exposed_lysine_context_count")) for r in group_rows)),
            "best_exposed_lys_fraction": _round_or_none(max(exposed_fracs), 3) if exposed_fracs else None,
            "ligand_context_class": representative.get("ligand_context_class"),
            "ligand_context_label": representative.get("ligand_context_label"),
            "ligand_priority": _numeric_value(representative.get("ligand_priority")),
            **_protacability_enrichment_snapshot(representative),
        }
        interpretation_label, interpretation_note = _target_interpretation(row)
        row["interpretation_label"] = interpretation_label
        row["interpretation_note"] = interpretation_note
        target_rows.append(row)
    return target_rows


def _sort_protacability_rows(rows, view, sort_value):
    sort_token = _protacability_sort_token(view, sort_value)
    field_name, reverse = PROTACABILITY_SORT_OPTIONS[view][sort_token]

    if sort_token in {"ligand_priority_desc", "readiness_priority_desc"}:
        pre_sorted = sorted(
            rows,
            key=lambda row: (
                str(row.get("virus_name") or "").lower(),
                str(row.get("protein_type") or row.get("display_protein_type") or "").lower(),
                str(row.get("pdb_code") or "").lower(),
                str(row.get("chain_id") or row.get("representative_chain_id") or "").lower(),
            ),
        )
        return sorted(
            pre_sorted,
            key=lambda row: (
                _numeric_value(row.get("has_candidate_linker_atoms")),
                _numeric_value(row.get("has_mapped_atoms")),
                _numeric_value(row.get("has_solvent_exposed_ligand_atoms")),
                _numeric_value(row.get("ligand_priority")),
                _numeric_value(row.get("degrader_design_readiness_score")),
                _numeric_value(row.get("warhead_linkability_score")),
                _numeric_value(row.get("target_lysine_accessibility_score")),
                _numeric_value(row.get("protein_structural_priority_score") or row.get("best_score") or row.get("protacability_proxy_score")),
                _numeric_value(row.get("best_score") or row.get("protacability_proxy_score")),
                _numeric_value(row.get("informative_ligand_count")),
                _numeric_value(row.get("best_exposed_lys_fraction") or row.get("exposed_lys_fraction")),
            ),
            reverse=True,
        ), sort_token

    def sort_key(row):
        value = row.get(field_name)
        if field_name in {"pdb_code", "protein_type", "display_protein_type"}:
            primary = str(value or "").lower()
        else:
            primary = _numeric_value(value)
        tie_breaker = (
            _numeric_value(row.get("has_candidate_linker_atoms")),
            _numeric_value(row.get("has_mapped_atoms")),
            _numeric_value(row.get("has_solvent_exposed_ligand_atoms")),
            _numeric_value(row.get("degrader_design_readiness_score")),
            _numeric_value(row.get("warhead_linkability_score")),
            _numeric_value(row.get("target_lysine_accessibility_score")),
            _numeric_value(row.get("protein_structural_priority_score") or row.get("best_score") or row.get("protacability_proxy_score")),
            _numeric_value(row.get("ligand_priority")),
            _numeric_value(row.get("best_score") or row.get("protacability_proxy_score")),
            _numeric_value(row.get("informative_ligand_count")),
            _numeric_value(row.get("best_exposed_lys_fraction") or row.get("exposed_lys_fraction")),
            -_numeric_value(row.get("candidate_ligand_count") or row.get("candidate_ligand_count_display")),
            str(row.get("virus_name") or "").lower(),
            str(row.get("protein_type") or row.get("display_protein_type") or "").lower(),
            str(row.get("pdb_code") or "").lower(),
            str(row.get("chain_id") or row.get("representative_chain_id") or "").lower(),
        )
        return (primary, tie_breaker)

    return sorted(rows, key=sort_key, reverse=reverse), sort_token


def _paginate_rows(rows, limit, offset):
    total_rows = len(rows)
    paged_rows = rows[offset: offset + limit]
    return paged_rows, total_rows, (offset + limit) < total_rows


def _build_summary_cards(view, rows):
    if view == "targets":
        score_values = [_numeric_value(row.get("degrader_design_readiness_score"), None) for row in rows if row.get("degrader_design_readiness_score") not in (None, "")]
        return {
            "targets_assessed": len(rows),
            "candidate_warheads_with_exposed_mapped_atoms": sum(
                1 for row in rows
                if _has_positive_value(row.get("solvent_exposed_mapped_atom_count")) and _has_positive_value(row.get("pdb_to_smiles_mapped_atom_count"))
            ),
            "high_warhead_linkability_ligands": sum(1 for row in rows if _numeric_value(row.get("warhead_linkability_score")) >= 70),
            "high_degrader_readiness_targets": sum(1 for row in rows if _numeric_value(row.get("degrader_design_readiness_score")) >= 70),
            "average_degrader_readiness": _round_or_none(sum(score_values) / len(score_values), 2) if score_values else None,
            "total_rows": len(rows),
        }
    if view == "protein":
        score_values = [_numeric_value(row.get("degrader_design_readiness_score"), None) for row in rows if row.get("degrader_design_readiness_score") not in (None, "")]
        return {
            "targets_assessed": len(rows),
            "candidate_warheads_with_exposed_mapped_atoms": sum(
                1 for row in rows
                if _has_positive_value(row.get("solvent_exposed_mapped_atom_count")) and _has_positive_value(row.get("pdb_to_smiles_mapped_atom_count"))
            ),
            "high_warhead_linkability_ligands": sum(1 for row in rows if _numeric_value(row.get("warhead_linkability_score")) >= 70),
            "high_degrader_readiness_targets": sum(1 for row in rows if _numeric_value(row.get("degrader_design_readiness_score")) >= 70),
            "average_degrader_readiness": _round_or_none(sum(score_values) / len(score_values), 2) if score_values else None,
            "total_rows": len(rows),
        }
    if view == "summary":
        score_values = [_numeric_value(row.get("degrader_design_readiness_score"), None) for row in rows if row.get("degrader_design_readiness_score") not in (None, "")]
        return {
            "targets_assessed": len(rows),
            "candidate_warheads_with_exposed_mapped_atoms": sum(
                1 for row in rows
                if _has_positive_value(row.get("solvent_exposed_mapped_atom_count")) and _has_positive_value(row.get("pdb_to_smiles_mapped_atom_count"))
            ),
            "high_warhead_linkability_ligands": sum(1 for row in rows if _numeric_value(row.get("warhead_linkability_score")) >= 70),
            "high_degrader_readiness_targets": sum(1 for row in rows if _numeric_value(row.get("degrader_design_readiness_score")) >= 70),
            "average_degrader_readiness": _round_or_none(sum(score_values) / len(score_values), 2) if score_values else None,
            "total_rows": len(rows),
        }
    score_values = [_numeric_value(row.get("degrader_design_readiness_score"), None) for row in rows if row.get("degrader_design_readiness_score") not in (None, "")]
    return {
        "targets_assessed": len(rows),
        "candidate_warheads_with_exposed_mapped_atoms": sum(
            1 for row in rows
            if _has_positive_value(row.get("solvent_exposed_mapped_atom_count")) and _has_positive_value(row.get("pdb_to_smiles_mapped_atom_count"))
        ),
        "high_warhead_linkability_ligands": sum(1 for row in rows if _numeric_value(row.get("warhead_linkability_score")) >= 70),
        "high_degrader_readiness_targets": sum(1 for row in rows if _numeric_value(row.get("degrader_design_readiness_score")) >= 70),
        "average_degrader_readiness": _round_or_none(sum(score_values) / len(score_values), 2) if score_values else None,
        "total_rows": len(rows),
    }


def _load_protacability_assessment_rows(conn, pdb_code=None):
    query = f"SELECT {', '.join(PROTACABILITY_CHAIN_COLUMNS)} FROM protacability_assessment"
    params = []
    if pdb_code:
        query += " WHERE pdb_code = ?"
        params.append(pdb_code)
    return conn.execute(query, params).fetchall()


def _load_protacability_enrichment_tables(conn):
    optional_tables = _protacability_optional_table_names(conn)
    readiness_rows = _load_optional_table_rows(conn, "protacability_degrader_readiness") if "protacability_degrader_readiness" in optional_tables else []
    warhead_rows = _load_optional_table_rows(conn, "protacability_warhead_linkability") if "protacability_warhead_linkability" in optional_tables else []
    return readiness_rows, warhead_rows


def _prepare_protacability_result_set(conn, args, export_all=False):
    view = _protacability_view_mode(args.get("view"))
    collapse_labels = _protacability_collapse_labels(args.get("collapse_labels"))
    requested_limit = args.get("page_size", type=int) or args.get("limit", type=int) or 50
    limit = min(max(requested_limit, 1), 100)
    page = args.get("page", type=int)
    if page and page > 0 and args.get("offset", type=int) is None:
        offset = (page - 1) * limit
    else:
        offset = max(args.get("offset", type=int) or 0, 0)
    filters = _build_protacability_filters(args)

    readiness_rows, warhead_rows = _load_protacability_enrichment_tables(conn)
    rows = _decorate_protacability_rows(
        _load_protacability_assessment_rows(conn),
        collapse_labels=collapse_labels,
        readiness_rows=readiness_rows,
        warhead_rows=warhead_rows,
    )
    filtered_rows = _filter_protacability_rows(rows, filters, collapse_labels=collapse_labels)

    if view == "targets":
        result_rows = _group_target_rows(filtered_rows)
    elif view == "protein":
        result_rows = _group_protein_rows(filtered_rows)
    elif view == "chains":
        result_rows = filtered_rows
    else:
        result_rows = _group_structure_rows(filtered_rows)

    result_rows = _apply_ligand_context_filter(result_rows, filters.get("ligand_context_class"))
    result_rows = _apply_ligand_presence_filter(result_rows, filters.get("ligand_presence"))

    sorted_rows, sort_token = _sort_protacability_rows(result_rows, view, args.get("sort"))
    summary = _build_summary_cards(view, sorted_rows)

    if export_all:
        page_rows = sorted_rows
        total_rows = len(sorted_rows)
        has_more = False
        page_offset = 0
    else:
        page_rows, total_rows, has_more = _paginate_rows(sorted_rows, limit, offset)
        page_offset = offset

    return {
        "view": view,
        "collapse_labels": collapse_labels,
        "rows": page_rows,
        "summary": summary,
        "total_rows": total_rows,
        "has_more": has_more,
        "limit": limit,
        "offset": page_offset,
        "sort": sort_token,
        "all_rows": sorted_rows,
        "filtered_chain_rows": filtered_rows,
    }


def _build_protacability_filter_options_payload_from_rows(assessment_rows, readiness_rows, warhead_rows, args):
    collapse_labels = _protacability_collapse_labels(args.get("collapse_labels"))
    filters = _build_protacability_filters(args)
    rows = _decorate_protacability_rows(
        assessment_rows,
        collapse_labels=collapse_labels,
        readiness_rows=readiness_rows,
        warhead_rows=warhead_rows,
    )

    virus_rows = _filter_options_for_context(rows, filters, collapse_labels=collapse_labels, ignore_key="virus_names")
    protein_rows = _filter_options_for_context(rows, filters, collapse_labels=collapse_labels, ignore_key="protein_types")
    tier_rows = _filter_options_for_context(rows, filters, collapse_labels=collapse_labels, ignore_key="tiers")
    ligand_rows = _filter_options_for_context(rows, filters, collapse_labels=collapse_labels, ignore_key="ligand")
    context_rows = _filter_options_for_context(rows, filters, collapse_labels=collapse_labels, ignore_key="ligand_context_class")

    protein_field = "display_protein_type" if collapse_labels else "protein_type"
    ligand_context_classes = [
        {"value": "candidate_small_molecule", "label": "Candidate small-molecule context"},
        {"value": "candidate_plus_glycan", "label": "Candidate + glycan context"},
        {"value": "glycan_only", "label": "Glycan-only context"},
        {"value": "glycan_common_mixed", "label": "Glycan/common structural context"},
        {"value": "common_buffer_only", "label": "Common buffer/crystal context"},
        {"value": "no_ligand_context", "label": "No ligand context"},
    ]
    available_contexts = {row.get("ligand_context_class") for row in context_rows if row.get("ligand_context_class")}

    return {
        "data_available": True,
        "virus_names": sorted({row.get("virus_name") for row in virus_rows if row.get("virus_name")}),
        "protein_types": sorted({row.get(protein_field) for row in protein_rows if row.get(protein_field)}),
        "tiers": sorted(
            {row.get("protacability_tier") for row in tier_rows if row.get("protacability_tier")},
            key=lambda tier: (_tier_rank(tier), tier),
        ),
        "warhead_tiers": sorted({row.get("warhead_linkability_tier") for row in rows if row.get("warhead_linkability_tier")}),
        "readiness_tiers": sorted({row.get("degrader_design_readiness_tier") for row in rows if row.get("degrader_design_readiness_tier")}),
        "evidence_levels": sorted({row.get("evidence_level") for row in rows if row.get("evidence_level")}),
        "smiles_sources": sorted({row.get("smiles_source") for row in rows if row.get("smiles_source")}),
        "ligands": sorted({ligand for row in ligand_rows for ligand in (row.get("ligand_names") or []) if ligand}),
        "ligand_context_classes": [item for item in ligand_context_classes if item["value"] in available_contexts],
    }


def _prepare_protacability_result_set_from_rows(assessment_rows, readiness_rows, warhead_rows, args, export_all=False):
    view = _protacability_view_mode(args.get("view"))
    collapse_labels = _protacability_collapse_labels(args.get("collapse_labels"))
    requested_limit = args.get("page_size", type=int) or args.get("limit", type=int) or 50
    limit = min(max(requested_limit, 1), 100)
    page = args.get("page", type=int)
    if page and page > 0 and args.get("offset", type=int) is None:
        offset = (page - 1) * limit
    else:
        offset = max(args.get("offset", type=int) or 0, 0)
    filters = _build_protacability_filters(args)

    rows = _decorate_protacability_rows(
        assessment_rows,
        collapse_labels=collapse_labels,
        readiness_rows=readiness_rows,
        warhead_rows=warhead_rows,
    )
    filtered_rows = _filter_protacability_rows(rows, filters, collapse_labels=collapse_labels)

    if view == "targets":
        result_rows = _group_target_rows(filtered_rows)
    elif view == "protein":
        result_rows = _group_protein_rows(filtered_rows)
    elif view == "chains":
        result_rows = filtered_rows
    else:
        result_rows = _group_structure_rows(filtered_rows)

    result_rows = _apply_ligand_context_filter(result_rows, filters.get("ligand_context_class"))
    result_rows = _apply_ligand_presence_filter(result_rows, filters.get("ligand_presence"))

    sorted_rows, sort_token = _sort_protacability_rows(result_rows, view, args.get("sort"))
    summary = _build_summary_cards(view, sorted_rows)

    if export_all:
        page_rows = sorted_rows
        total_rows = len(sorted_rows)
        has_more = False
        page_offset = 0
    else:
        page_rows, total_rows, has_more = _paginate_rows(sorted_rows, limit, offset)
        page_offset = offset

    return {
        "view": view,
        "collapse_labels": collapse_labels,
        "rows": page_rows,
        "summary": summary,
        "total_rows": total_rows,
        "has_more": has_more,
        "limit": limit,
        "offset": page_offset,
        "sort": sort_token,
        "all_rows": sorted_rows,
        "filtered_chain_rows": filtered_rows,
    }


def _pick_representative_ligand_record(ligand_inventory, preferred_ligands=None, allow_glycan=False, preferred_chain=None):
    preferred = {str(name).strip().upper() for name in (preferred_ligands or []) if str(name).strip()}
    excluded = set(PROTACABILITY_COMMON_CONTEXT_LIGANDS) | (set() if allow_glycan else set(PROTACABILITY_GLYCAN_LIGANDS))
    candidates = []
    for record in ligand_inventory or []:
        resname = str(record.get("ligand_resname") or "").strip().upper()
        if not resname:
            continue
        if resname in excluded:
            continue
        if preferred and resname in preferred:
            candidates.append(record)
    if not candidates:
        candidates = [
            record for record in (ligand_inventory or [])
            if str(record.get("ligand_resname") or "").strip()
            and str(record.get("ligand_resname") or "").strip().upper() not in excluded
        ]
    if not candidates and allow_glycan:
        candidates = [record for record in (ligand_inventory or []) if str(record.get("ligand_resname") or "").strip()]
    if not candidates:
        return None
    return sorted(
        candidates,
        key=lambda r: (
            0 if str(r.get("ligand_chain") or "").strip().upper() == str(preferred_chain or "").strip().upper() and preferred_chain else 1,
            -_numeric_value(r.get("ligand_atom_count")),
            -_numeric_value(r.get("ligand_heavy_atom_count")),
            str(r.get("ligand_resname") or ""),
            str(r.get("ligand_chain") or ""),
            str(r.get("ligand_residue_id") or ""),
        ),
    )[0]


def _serialize_ligand_contexts(ligand_inventory):
    contexts = []
    for row in ligand_inventory or []:
        contexts.append({
            "label": " ".join(part for part in [
                str(row.get("ligand_resname") or "").strip(),
                str(row.get("ligand_chain") or "").strip(),
                str(row.get("ligand_residue_id") or "").strip(),
            ] if part),
            "ligand_resname": row.get("ligand_resname"),
            "ligand_chain": row.get("ligand_chain"),
            "ligand_residue_id": row.get("ligand_residue_id"),
            "ligand_atom_count": row.get("ligand_atom_count"),
            "ligand_heavy_atom_count": row.get("ligand_heavy_atom_count"),
        })
    return contexts

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

    @bp.get("/protacability/filter-options")
    @require_token
    def get_protacability_filter_options():
        payload = _protacability_source_payload_local()
        if not payload.get("data_available"):
            return jsonify({
                "data_available": False,
                "message": "PROTACability data has not been imported yet. Run tools/import_protacability_data.py after generating the expansion outputs.",
                "virus_names": [],
                "protein_types": [],
                "tiers": [],
                "warhead_tiers": [],
                "readiness_tiers": [],
                "evidence_levels": [],
                "smiles_sources": [],
                "ligands": [],
                "ligand_context_classes": [],
            })
        return jsonify(
            _build_protacability_filter_options_payload_from_rows(
                payload.get("assessment_rows", []),
                payload.get("readiness_rows", []),
                payload.get("warhead_rows", []),
                request.args,
            )
        )

    @bp.get("/protacability/search")
    @require_token
    def get_protacability_search():
        payload = _protacability_source_payload_local()
        if not payload.get("data_available"):
            return jsonify({
                "data_available": False,
                "message": "PROTACability data has not been imported yet. Run tools/import_protacability_data.py after generating the expansion outputs.",
                "rows": [],
                "summary": {},
            })
        result = _prepare_protacability_result_set_from_rows(
            payload.get("assessment_rows", []),
            payload.get("readiness_rows", []),
            payload.get("warhead_rows", []),
            request.args,
        )
        return jsonify({
            "data_available": True,
            "view": result["view"],
            "collapse_labels": result["collapse_labels"],
            "rows": result["rows"],
            "summary": result["summary"],
            "limit": result["limit"],
            "offset": result["offset"],
            "total_rows": result["total_rows"],
            "has_more": result["has_more"],
            "sort": result["sort"],
        })

    @bp.get("/protacability/detail/<pdb_code>/<chain_id>")
    @require_token
    def get_protacability_detail(pdb_code: str, chain_id: str):
        payload = _protacability_source_payload_local(
            pdb_code=pdb_code,
            include_lysine=True,
            include_inventory=True,
        )
        if not payload.get("data_available"):
            return _json_error("PROTACability data has not been imported yet. Run tools/import_protacability_data.py after generating the expansion outputs.", 404)

        readiness_rows = payload.get("readiness_rows", [])
        warhead_rows = payload.get("warhead_rows", [])
        raw_assessment_rows = [
            row for row in payload.get("assessment_rows", [])
            if row.get("pdb_code") == pdb_code and row.get("chain_id") == chain_id
        ]
        assessment_rows = _decorate_protacability_rows(
            raw_assessment_rows,
            collapse_labels=True,
            readiness_rows=readiness_rows,
            warhead_rows=warhead_rows,
        )
        assessment = max(assessment_rows, key=_row_priority_key) if assessment_rows else None
        if assessment is None:
            return _json_error("Assessment row not found", 404)

        lysine_rows = [
            {
                "lys_residue_id": row.get("lys_residue_id"),
                "lys_observed_index": row.get("lys_observed_index"),
                "lysine_sasa_a2": row.get("lysine_sasa_a2"),
                "is_surface_exposed": row.get("is_surface_exposed"),
                "nearest_ligand_resname": row.get("nearest_ligand_resname"),
                "nearest_ligand_distance_a": row.get("nearest_ligand_distance_a"),
                "linker_site_class": row.get("linker_site_class"),
            }
            for row in payload.get("lysine_rows", [])
            if row.get("pdb_code") == pdb_code and row.get("chain_id") == chain_id
        ]
        ligand_inventory = [
            {
                "ligand_resname": row.get("ligand_resname"),
                "ligand_chain": row.get("ligand_chain"),
                "ligand_residue_id": row.get("ligand_residue_id"),
                "ligand_atom_count": row.get("ligand_atom_count"),
                "ligand_heavy_atom_count": row.get("ligand_heavy_atom_count"),
                "centroid_x": row.get("centroid_x"),
                "centroid_y": row.get("centroid_y"),
                "centroid_z": row.get("centroid_z"),
            }
            for row in payload.get("ligand_inventory", [])
            if row.get("pdb_code") == pdb_code
        ]
        related_chains = [
            {
                "chain_id": row.get("chain_id"),
                "protacability_proxy_score": row.get("protacability_proxy_score"),
                "protacability_tier": row.get("protacability_tier"),
                "candidate_ligand_count": row.get("candidate_ligand_count"),
                "exposed_lys_count": row.get("exposed_lys_count"),
            }
            for row in payload.get("assessment_rows", [])
            if row.get("pdb_code") == pdb_code
        ]
        return jsonify({
            "data_available": True,
            "assessment": dict(assessment),
            "lysine_rows": lysine_rows,
            "ligand_inventory": ligand_inventory,
            "ligand_contexts": _serialize_ligand_contexts(ligand_inventory),
            "related_chains": related_chains,
        })

    @bp.get("/protacability/structure-detail/<pdb_code>")
    @require_token
    def get_protacability_structure_detail(pdb_code: str):
        payload = _protacability_source_payload_local(
            pdb_code=pdb_code,
            include_lysine=True,
            include_inventory=True,
        )
        if not payload.get("data_available"):
            return _json_error("PROTACability data has not been imported yet. Run tools/import_protacability_data.py after generating the expansion outputs.", 404)

        collapse_labels = _protacability_collapse_labels(request.args.get("collapse_labels"))
        virus_name = str(request.args.get("virus_name", "") or "").strip()
        protein_type = str(request.args.get("protein_type", "") or "").strip()
        readiness_rows = payload.get("readiness_rows", [])
        warhead_rows = payload.get("warhead_rows", [])
        decorated_rows = _decorate_protacability_rows(
            payload.get("assessment_rows", []),
            collapse_labels=collapse_labels,
            readiness_rows=readiness_rows,
            warhead_rows=warhead_rows,
        )
        if virus_name:
            decorated_rows = [row for row in decorated_rows if row.get("virus_name") == virus_name]
        if protein_type:
            decorated_rows = [row for row in decorated_rows if (row.get("display_protein_type") if collapse_labels else row.get("protein_type")) == protein_type]
        if not decorated_rows:
            return _json_error("Structure summary not found", 404)

        summary_rows = _group_structure_rows(decorated_rows)
        summary_row = max(summary_rows, key=lambda row: (_numeric_value(row.get("best_score")), _numeric_value(row.get("best_exposed_lys_fraction"))))
        chain_rows = _dedupe_display_chain_rows(decorated_rows)
        chain_rows, _ = _sort_protacability_rows(chain_rows, "chains", "protacability_proxy_score_desc")
        representative_chain = summary_row.get("representative_chain_id")
        lysine_rows = [
            {
                "lys_residue_id": row.get("lys_residue_id"),
                "lys_observed_index": row.get("lys_observed_index"),
                "lysine_sasa_a2": row.get("lysine_sasa_a2"),
                "is_surface_exposed": row.get("is_surface_exposed"),
                "nearest_ligand_resname": row.get("nearest_ligand_resname"),
                "nearest_ligand_distance_a": row.get("nearest_ligand_distance_a"),
                "linker_site_class": row.get("linker_site_class"),
            }
            for row in payload.get("lysine_rows", [])
            if row.get("pdb_code") == pdb_code and row.get("chain_id") == representative_chain
        ]
        ligand_inventory = [
            {
                "ligand_resname": row.get("ligand_resname"),
                "ligand_chain": row.get("ligand_chain"),
                "ligand_residue_id": row.get("ligand_residue_id"),
                "ligand_atom_count": row.get("ligand_atom_count"),
                "ligand_heavy_atom_count": row.get("ligand_heavy_atom_count"),
                "centroid_x": row.get("centroid_x"),
                "centroid_y": row.get("centroid_y"),
                "centroid_z": row.get("centroid_z"),
            }
            for row in payload.get("ligand_inventory", [])
            if row.get("pdb_code") == pdb_code
        ]
        preferred_ligands = _split_candidate_ligands(summary_row.get("candidate_ligand_resnames_full"))
        representative_ligand = _pick_representative_ligand_record(
            ligand_inventory,
            preferred_ligands=preferred_ligands,
            allow_glycan=summary_row.get("ligand_context_class") == "glycan_only",
            preferred_chain=representative_chain,
        )
        return jsonify({
            "data_available": True,
            "summary": summary_row,
            "chain_rows": chain_rows,
            "representative_chain_id": representative_chain,
            "representative_ligand": representative_ligand,
            "representative_ligand_resname": (representative_ligand or {}).get("ligand_resname"),
            "representative_ligand_chain": (representative_ligand or {}).get("ligand_chain"),
            "representative_ligand_residue_id": (representative_ligand or {}).get("ligand_residue_id"),
            "lysine_rows": lysine_rows,
            "ligand_inventory": ligand_inventory,
            "ligand_contexts": _serialize_ligand_contexts(ligand_inventory),
        })

    @bp.get("/protacability/protein-detail")
    @require_token
    def get_protacability_protein_detail():
        virus_name = str(request.args.get("virus_name", "") or "").strip()
        protein_type = str(request.args.get("protein_type", "") or "").strip()
        payload = _protacability_source_payload_local(virus_name=virus_name, protein_type=protein_type)
        if not payload.get("data_available"):
            return _json_error("PROTACability data has not been imported yet. Run tools/import_protacability_data.py after generating the expansion outputs.", 404)

        collapse_labels = _protacability_collapse_labels(request.args.get("collapse_labels"))
        readiness_rows = payload.get("readiness_rows", [])
        warhead_rows = payload.get("warhead_rows", [])
        decorated_rows = _decorate_protacability_rows(
            payload.get("assessment_rows", []),
            collapse_labels=collapse_labels,
            readiness_rows=readiness_rows,
            warhead_rows=warhead_rows,
        )
        decorated_rows = [
            row for row in decorated_rows
            if row.get("virus_name") == virus_name
            and ((row.get("display_protein_type") if collapse_labels else row.get("protein_type")) == protein_type)
        ]
        if not decorated_rows:
            return _json_error("Protein summary not found", 404)

        protein_rows = _group_protein_rows(decorated_rows)
        protein_row = protein_rows[0]
        structure_rows = _group_structure_rows(decorated_rows)
        structure_rows, _ = _sort_protacability_rows(structure_rows, "summary", "best_score_desc")
        tier_distribution = dict(Counter(row.get("best_tier") or "Unknown" for row in structure_rows))
        return jsonify({
            "data_available": True,
            "summary": protein_row,
            "top_structures": structure_rows[:10],
            "tier_distribution": tier_distribution,
            "explanation": "This view groups multiple structures and chains into a single protein-level summary so repeated biological contexts do not dominate the table.",
        })

    @bp.get("/protacability/target-detail")
    @require_token
    def get_protacability_target_detail():
        collapse_labels = _protacability_collapse_labels(request.args.get("collapse_labels"))
        virus_name = str(request.args.get("virus_name", "") or "").strip()
        protein_type = str(request.args.get("protein_type", "") or "").strip()
        ligand_context_class = str(request.args.get("ligand_context_class", "") or "").strip()
        min_score = request.args.get("min_score", type=float)
        if not virus_name or not protein_type:
            return _json_error("virus_name and protein_type are required", 400)

        payload = _protacability_source_payload_local(
            virus_name=virus_name,
            protein_type=protein_type,
            include_inventory=True,
        )
        if not payload.get("data_available"):
            return _json_error("PROTACability data is not available.", 404)

        readiness_rows = payload.get("readiness_rows", [])
        warhead_rows = payload.get("warhead_rows", [])
        rows = _decorate_protacability_rows(
            payload.get("assessment_rows", []),
            collapse_labels=collapse_labels,
            readiness_rows=readiness_rows,
            warhead_rows=warhead_rows,
        )
        rows = [row for row in rows if row.get("virus_name") == virus_name and (row.get("display_protein_type") or row.get("protein_type")) == protein_type]
        if min_score is not None:
            rows = [row for row in rows if _numeric_value(row.get("protacability_proxy_score"), -1) >= min_score]
        if ligand_context_class:
            rows = _apply_ligand_context_filter(rows, ligand_context_class)
        if not rows:
            return _json_error("Target detail not found", 404)

        target_summary = _group_target_rows(rows)[0]
        structure_rows = _group_structure_rows(rows)
        structure_rows, _ = _sort_protacability_rows(structure_rows, "summary", "best_score_desc")
        ligand_groups_map = defaultdict(list)
        for srow in structure_rows:
            ligands = _split_candidate_ligands(srow.get("candidate_ligand_resnames_full"))
            for ligand in ligands:
                ligand_groups_map[ligand].append(srow)
        ligand_groups = []
        for ligand, grows in ligand_groups_map.items():
            top = max(grows, key=lambda r: _numeric_value(r.get("best_score")))
            ligand_groups.append({
                "ligand_resname": ligand,
                "context_class": top.get("ligand_context_class"),
                "pdb_count": len({r.get("pdb_code") for r in grows}),
                "top_pdb_code": top.get("pdb_code"),
                "associated_protein_chains": sorted({r.get("representative_chain_id") for r in grows if r.get("representative_chain_id")}),
                "best_score": top.get("best_score"),
            })
        ligand_groups = sorted(ligand_groups, key=lambda r: (_ligand_context_rank(r.get("context_class")), _numeric_value(r.get("best_score"))), reverse=True)

        def _pick_context(context_class: str):
            matches = [row for row in structure_rows if row.get("ligand_context_class") == context_class]
            if not matches:
                return None
            best = max(matches, key=lambda r: _numeric_value(r.get("best_score")))
            ligand_rows = [
                {
                    "ligand_resname": row.get("ligand_resname"),
                    "ligand_chain": row.get("ligand_chain"),
                    "ligand_residue_id": row.get("ligand_residue_id"),
                    "ligand_atom_count": row.get("ligand_atom_count"),
                    "ligand_heavy_atom_count": row.get("ligand_heavy_atom_count"),
                }
                for row in payload.get("ligand_inventory", [])
                if row.get("pdb_code") == best.get("pdb_code")
            ]
            preferred = _split_candidate_ligands(best.get("candidate_ligand_resnames_full"))
            ligand_record = _pick_representative_ligand_record(
                ligand_rows,
                preferred_ligands=preferred,
                allow_glycan=context_class == "glycan_only",
                preferred_chain=best.get("representative_chain_id"),
            )
            return {
                "pdb_code": best.get("pdb_code"),
                "chain_id": best.get("representative_chain_id"),
                "ligand_resname": (ligand_record or {}).get("ligand_resname"),
                "ligand_chain": (ligand_record or {}).get("ligand_chain"),
                "ligand_residue_id": (ligand_record or {}).get("ligand_residue_id"),
                "ligand_context_class": best.get("ligand_context_class"),
            }

        representative_contexts = {
            "best_candidate_small_molecule_context": _pick_context("candidate_small_molecule"),
            "best_candidate_plus_glycan_context": _pick_context("candidate_plus_glycan"),
            "best_glycan_only_context": _pick_context("glycan_only"),
        }
        representative_ligand = (
            representative_contexts.get("best_candidate_small_molecule_context")
            or representative_contexts.get("best_candidate_plus_glycan_context")
            or representative_contexts.get("best_glycan_only_context")
        )
        active_pdb_code = (representative_ligand or {}).get("pdb_code") or target_summary.get("best_pdb_code")
        ligand_inventory = [
            {
                "ligand_resname": row.get("ligand_resname"),
                "ligand_chain": row.get("ligand_chain"),
                "ligand_residue_id": row.get("ligand_residue_id"),
                "ligand_atom_count": row.get("ligand_atom_count"),
                "ligand_heavy_atom_count": row.get("ligand_heavy_atom_count"),
            }
            for row in payload.get("ligand_inventory", [])
            if row.get("pdb_code") == active_pdb_code
        ]
        return jsonify({
            "data_available": True,
            "target_summary": target_summary,
            "structure_rows": structure_rows,
            "ligand_groups": ligand_groups,
            "representative_contexts": representative_contexts,
            "representative_ligand": representative_ligand,
            "ligand_contexts": _serialize_ligand_contexts(ligand_inventory),
        })

    @bp.get("/protacability/export-filtered")
    @require_token
    def get_protacability_export_filtered():
        payload = _protacability_source_payload_local()
        if not payload.get("data_available"):
            return _json_error("PROTACability data has not been imported yet. Run tools/import_protacability_data.py after generating the expansion outputs.", 404)
        result = _prepare_protacability_result_set_from_rows(
            payload.get("assessment_rows", []),
            payload.get("readiness_rows", []),
            payload.get("warhead_rows", []),
            request.args,
            export_all=True,
        )
        return jsonify({
            "data_available": True,
            "view": result["view"],
            "rows": result["rows"],
            "sort": result["sort"],
            "total_rows": result["total_rows"],
        })

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
