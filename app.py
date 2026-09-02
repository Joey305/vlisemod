import os 
import io
import json
import re
import glob
import urllib.request
import urllib.error
import requests
from pathlib import Path
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from PIL import Image, ImageDraw, ImageFont
from flask import Flask, render_template, request, jsonify, redirect, url_for, send_from_directory, send_file, render_template_string
import sqlite3
import logging
import time
from collections import Counter, defaultdict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import zipfile
import shutil
import threading
import uuid
from io import BytesIO
import zipfile
from datetime import timedelta
from flask import request, jsonify, render_template
try:
    from Bio.PDB import MMCIFParser, PDBIO
except Exception:  # optional at runtime
    MMCIFParser = None
    PDBIO = None

LIGAND_INSTANCE_MAPPING_CACHE = {}

from rdkit import Chem
''' 
How to Start the File:

waitress-serve --listen=127.0.0.1:5002 app:app

'''

# Set up logging
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

# Initialize the Flask app
app = Flask(__name__)


class RandyBackendError(Exception):
    def __init__(self, message, status_code=502):
        super().__init__(message)
        self.status_code = status_code


LOCAL_DB_PATH = Path(os.environ.get("VLISMOD_LOCAL_DB_PATH", "viral_data.db")).expanduser()
LIGAND_IMAGE_REQUIRED_TABLES = (
    "Functional_GROUPED",
    "Ligand_Arp_Diagram",
    "solvent_exposed_atoms",
    "SMILES_MAP_PDB",
)
COMPARE_LIGAND_REQUIRED_TABLES = (
    "structures",
    "structure_classifications",
    "ligand_instances",
    "ligand_instance_atoms",
    "ligand_smiles_atom_mapping",
    "ligand_sasa_atoms",
    "ligand_arpeggio_runs",
    "arpeggio_raw_contact_labels",
)
LIGAND_INTERACTION_REQUIRED_TABLES = (
    "structures",
    "ligand_instances",
    "ligand_instance_atoms",
    "ligand_arpeggio_runs",
    "arpeggio_raw_contact_labels",
)
LIGAND_SYNONYM_REQUIRED_TABLES = ("Ligand_Synonyms",)
LIGAND_INFO_REQUIRED_TABLES = ("Ligand_Atoms_Smiles", "ligands")
LIGAND_OPTIONS_REQUIRED_TABLES = ("Ligand_Arp_Diagram",)
QUERY_PROTEIN_REQUIRED_TABLES = (
    "Virus_Proteins",
    "Ligand_Synonyms",
    "Ligand_Arp_Diagram",
    "structures",
    "ligand_instances",
    "ligand_instance_atoms",
)
PYMOL_REQUIRED_TABLES = (
    "ligand_atoms",
    "Functional_Group_Atoms",
    "receptor_binding_pocket",
    "distal_atoms",
    "solvent_exposed_atoms",
    "RUPLEY_SASA_DATA",
)
PROTACABILITY_ALL_TABLES = (
    "protacability_assessment",
    "protacability_lysine_proximity",
    "protacability_ligand_inventory",
    "protacability_warhead_linkability",
    "protacability_degrader_readiness",
)
DATA_SET_REQUIRED_TABLES = {
    "Solvent Exposed Atoms": ("solvent_exposed_atoms",),
    "Ligand Atoms": ("ligand_atoms",),
    "Binding Pocket": ("receptor_binding_pocket",),
    "Smiles and Functional Groups": ("Ligand_Atoms_Smiles",),
    "Interatomic Interactions": ("Arpeggio_Contacts_Data",),
    "Functional Group Atoms": ("Functional_Group_Atoms",),
    "Smiles & PDB Mapping": ("SMILES_MAP_PDB",),
    "PROTACability Assessment": ("protacability_assessment",),
    "PROTACability Lysine Proximity": ("protacability_lysine_proximity",),
    "PROTACability Ligand Inventory": ("protacability_ligand_inventory",),
    "PROTACability Warhead Linkability": ("protacability_warhead_linkability",),
    "PROTACability Degrader Readiness": ("protacability_degrader_readiness",),
}


def _normalized_backend_mode():
    mode = os.environ.get("VLISMOD_DATA_BACKEND", "local").strip().lower()
    if mode not in {"local", "randy", "auto"}:
        return "local"
    return mode


def _vlismod_backup_url():
    return os.environ.get("VLISMOD_BACKUP_URL", "").strip().rstrip("/")


def _randy_api_base_url():
    return os.environ.get("RANDY_API_BASE_URL", "").strip().rstrip("/")


def _randy_base_url():
    backup_url = _vlismod_backup_url()
    if backup_url:
        return backup_url

    legacy_base = _randy_api_base_url()
    if not legacy_base:
        return ""

    if legacy_base.endswith("/api/vlismod") or legacy_base.endswith("/backup/vlismod"):
        return legacy_base

    return f"{legacy_base}/api/vlismod"


def _randy_api_token():
    return (
        os.environ.get("RANDY_API_TOKEN", "").strip()
        or os.environ.get("VLISMOD_API_TOKEN", "").strip()
    )


def _randy_api_timeout_seconds(default=30):
    raw_value = os.environ.get("RANDY_API_TIMEOUT_SECONDS", "").strip()
    if not raw_value:
        return default
    try:
        return max(1, int(raw_value))
    except ValueError:
        return default


def use_randy_backend():
    return _normalized_backend_mode() == "randy"


def randy_available():
    return bool(_randy_base_url() and _randy_api_token())


def _enforce_randy_json_size(response, *, max_bytes):
    if max_bytes is None:
        return
    raw_length = response.headers.get("Content-Length", "").strip()
    if not raw_length:
        return
    try:
        content_length = int(raw_length)
    except ValueError:
        return
    if content_length > max_bytes:
        raise RandyBackendError(
            f"RANDY API payload too large for UI route ({content_length} bytes). Use a paginated or targeted endpoint instead.",
            status_code=502,
        )


def randy_get(path, params=None, *, max_bytes=10 * 1024 * 1024):
    if not randy_available():
        raise RandyBackendError("RANDY API is not configured.", status_code=500)

    url = f"{_randy_base_url()}/{str(path or '').lstrip('/')}"
    headers = {"Authorization": f"Bearer {_randy_api_token()}"}

    try:
        response = requests.get(
            url,
            params=params,
            headers=headers,
            timeout=_randy_api_timeout_seconds(30),
        )
    except requests.RequestException as exc:
        raise RandyBackendError(f"RANDY API request failed: {exc}", status_code=502) from exc

    _enforce_randy_json_size(response, max_bytes=max_bytes)

    try:
        payload = response.json()
    except ValueError:
        payload = None

    if response.status_code >= 400:
        message = payload.get("error") if isinstance(payload, dict) else None
        raise RandyBackendError(
            message or f"RANDY API request failed with status {response.status_code}.",
            status_code=response.status_code,
        )

    if payload is None:
        raise RandyBackendError("RANDY API returned non-JSON response.", status_code=502)

    return payload


def randy_post(path, json=None, *, max_bytes=10 * 1024 * 1024):
    if not randy_available():
        raise RandyBackendError("RANDY API is not configured.", status_code=500)

    url = f"{_randy_base_url()}/{str(path or '').lstrip('/')}"
    headers = {"Authorization": f"Bearer {_randy_api_token()}"}

    try:
        response = requests.post(
            url,
            json=json,
            headers=headers,
            timeout=_randy_api_timeout_seconds(60),
        )
    except requests.RequestException as exc:
        raise RandyBackendError(f"RANDY API request failed: {exc}", status_code=502) from exc

    _enforce_randy_json_size(response, max_bytes=max_bytes)

    try:
        payload = response.json()
    except ValueError:
        payload = None

    if response.status_code >= 400:
        message = payload.get("error") if isinstance(payload, dict) else None
        raise RandyBackendError(
            message or f"RANDY API request failed with status {response.status_code}.",
            status_code=response.status_code,
        )

    if payload is None:
        raise RandyBackendError("RANDY API returned non-JSON response.", status_code=502)

    return payload


def local_tables_available(required_tables):
    db_path = LOCAL_DB_PATH
    if not db_path.exists():
        return False, f"Local V-LiSEMOD database not found at {db_path}."

    with sqlite3.connect(str(db_path)) as conn:
        rows = conn.execute(
            "SELECT name FROM sqlite_master WHERE type IN ('table','view') AND name IN ({})".format(
                ",".join("?" for _ in required_tables)
            ),
            tuple(required_tables),
        ).fetchall()

    existing = {row[0] for row in rows}
    missing = [table for table in required_tables if table not in existing]
    if missing:
        return False, f"Local V-LiSEMOD database is missing required tables: {', '.join(missing)}."

    return True, None


def _connect_local_db(required_tables=None):
    if required_tables:
        tables_ok, error_message = local_tables_available(required_tables)
        if not tables_ok:
            raise RandyBackendError(error_message)
    return sqlite3.connect(str(LOCAL_DB_PATH))


def _form_flag(name):
    return request.form.get(name) is not None


def _remote_db_enabled():
    mode = _normalized_backend_mode()
    return mode == "randy" or (mode == "auto" and randy_available())


def _remote_page_metadata():
    return randy_get("page-metadata")


def _local_page_metadata():
    return {
        "available_export_data_sets": get_available_export_data_sets(),
        "protacability_data_available": protacability_tables_available(),
    }


def _required_tables_for_datasets(data_sets):
    required = []
    for data_set in data_sets:
        for table_name in DATA_SET_REQUIRED_TABLES.get(data_set, ()):
            if table_name not in required:
                required.append(table_name)
    return tuple(required)


def _local_get_ligands_with_synonyms_payload():
    with _connect_local_db(LIGAND_SYNONYM_REQUIRED_TABLES) as conn:
        cursor = conn.cursor()
        cursor.execute("SELECT ligand, synonym FROM Ligand_Synonyms")
        return [{'ligand_code': row[0], 'synonym': row[1]} for row in cursor.fetchall()]


def _local_get_ligand_info_payload(ligand_code):
    with _connect_local_db(LIGAND_INFO_REQUIRED_TABLES) as conn:
        cursor = conn.cursor()
        cursor.execute(
            """
            SELECT ligand, pdb_id, smiles, molecular_weight
            FROM Ligand_Atoms_Smiles
            WHERE ligand = ?
            """,
            (ligand_code,),
        )
        ligand_data = cursor.fetchone()

    if not ligand_data:
        return None

    return {
        "ligand": ligand_data[0],
        "pdb_id": ligand_data[1],
        "smiles": ligand_data[2],
        "molecular_weight": ligand_data[3],
    }


def _local_get_ligand_options_payload():
    with _connect_local_db(LIGAND_OPTIONS_REQUIRED_TABLES) as conn:
        cursor = conn.cursor()
        cursor.execute("SELECT virus_name, pdb_id, ligand, chain, ligand_id FROM Ligand_Arp_Diagram")
        rows = cursor.fetchall()
    return {
        "options": [
            {'value': f"{row[1]}-{row[2]}", 'text': f"{row[0]}, {row[1]}, {row[2]}, {row[3]}, {row[4]}"}
            for row in rows
        ]
    }


def _local_get_smiles_payload(ligand_code):
    with _connect_local_db(LIGAND_INFO_REQUIRED_TABLES) as conn:
        cursor = conn.cursor()
        # Stage 07 stores smiles_atom_index values against ligands.smiles.  Do
        # not substitute the canonicalized display SMILES from
        # Ligand_Atoms_Smiles here: it can describe the same molecule with a
        # different atom order, causing hover highlights to land on the wrong
        # rendered atom.
        cursor.execute(
            "SELECT smiles FROM ligands WHERE component_id = ? AND smiles IS NOT NULL LIMIT 1",
            (ligand_code,),
        )
        result = cursor.fetchone()
    return {"ligand_id": ligand_code, "smiles": result[0]} if result and result[0] else None


def _remote_smiles_payload(ligand_code):
    return randy_get("ligand-smiles", params={"ligand_id": ligand_code})


def _load_smiles_payload(ligand_code):
    mode = _normalized_backend_mode()
    if mode == "randy":
        return _remote_smiles_payload(ligand_code)
    if mode == "auto" and randy_available():
        try:
            return _remote_smiles_payload(ligand_code)
        except RandyBackendError:
            logging.warning("Falling back to local SMILES lookup for ligand %s", ligand_code)
    return _local_get_smiles_payload(ligand_code)


def _load_smiles_for_ligand(ligand_code):
    payload = _load_smiles_payload(ligand_code)
    return payload.get("smiles") if payload else None


def _local_interaction_records_payload(pdb_id, ligand, ligand_id, chain, ligand_instance_id):
    """Load one retained occurrence's latest Stage-09 interaction labels."""
    with _connect_local_db(LIGAND_INTERACTION_REQUIRED_TABLES) as conn:
        query = '''
            SELECT
                s.entry_id AS pdb_id,
                i.label_comp_id AS ligand,
                i.auth_asym_id AS chain,
                r.interaction_label AS Contact,
                r.distance AS Distance,
                COALESCE(a.auth_atom_id, a.label_atom_id) AS exact_atom,
                a.atom_site_id AS atom_id,
                json_extract(r.partner_identity_json, '$.label_comp_id') AS residue,
                json_extract(r.partner_identity_json, '$.auth_seq_id') AS residue_number,
                COALESCE(
                    json_extract(r.partner_identity_json, '$.auth_atom_id'),
                    json_extract(r.partner_identity_json, '$.label_atom_id')
                ) AS residue_atom,
                json_extract(r.partner_identity_json, '$.auth_asym_id') AS residue_chain,
                i.auth_seq_id AS ligand_id,
                i.ligand_instance_id,
                i.deposited_model_num AS model_id
            FROM ligand_instances AS i
            JOIN structures AS s ON s.structure_id = i.structure_id
            JOIN arpeggio_raw_contact_labels AS r
              ON r.ligand_instance_id = i.ligand_instance_id
             AND r.filter_class = 'raw_environment'
             AND r.run_id = (
                 SELECT MAX(ar.run_id)
                 FROM ligand_arpeggio_runs AS ar
                 WHERE ar.ligand_instance_id = i.ligand_instance_id
                   AND ar.status = 'completed'
             )
            LEFT JOIN ligand_instance_atoms AS a
              ON a.ligand_instance_atom_id = r.ligand_instance_atom_id
            WHERE i.ligand_instance_id = ?
              AND s.entry_id = ?
              AND i.label_comp_id = ?
              AND i.auth_seq_id = ?
              AND i.auth_asym_id = ?
              AND i.curation_status = 'included'
        '''
        df = pd.read_sql(
            query,
            conn,
            params=(ligand_instance_id, pdb_id, ligand, ligand_id, chain),
        )
    return {"records": df.replace({np.nan: None}).to_dict(orient="records")}


def _records_to_interaction_dataframe(records):
    if not records:
        return pd.DataFrame()
    return pd.DataFrame(records)


def _dispatch_supported_lookup(remote_path, *, params=None, local_loader):
    mode = _normalized_backend_mode()

    if mode == "local":
        return local_loader()

    if mode == "randy":
        try:
            payload = randy_get(remote_path, params=params)
        except RandyBackendError as exc:
            return jsonify({"error": str(exc)}), exc.status_code
        return jsonify(payload)

    if randy_available():
        try:
            payload = randy_get(remote_path, params=params)
            return jsonify(payload)
        except RandyBackendError:
            logging.warning("Falling back to local V-LiSEMOD database for %s", remote_path)

    return local_loader()


def _env_flag(name, default=False):
    raw = os.environ.get(name)
    if raw is None:
        return default
    return raw.strip().lower() in {"1", "true", "yes", "on"}


SHOW_DRUG_GPT_NAV = _env_flag("SHOW_DRUG_GPT_NAV", default=False)
ENABLE_DRUG_GPT = _env_flag("ENABLE_DRUG_GPT", default=False)
ENABLE_LOCAL_LLM = _env_flag("ENABLE_LOCAL_LLM", default=False)

configure_local_llm = None

if ENABLE_DRUG_GPT:
    from DRUGapp import dp, configure_local_llm
    app.register_blueprint(dp, url_prefix='/drugapp')




app.config.update(
    SECRET_KEY=os.environ.get("FLASK_SECRET_KEY", "dev-secret-change-me"),
)
app.config["SHOW_DRUG_GPT_NAV"] = SHOW_DRUG_GPT_NAV
app.config["ENABLE_DRUG_GPT"] = ENABLE_DRUG_GPT
app.config["ENABLE_LOCAL_LLM"] = ENABLE_LOCAL_LLM

if configure_local_llm is not None:
    configure_local_llm(app)
else:
    app.config["LOCAL_LLM_ENABLED"] = False
    app.config["LOCAL_LLM_MODEL_ID"] = os.environ.get("MODEL_ID") or os.environ.get("LLM_MODEL_ID")
    app.config["LLAMA_MODEL"] = None
    app.config["LLAMA_TOKENIZER"] = None
    app.config["LOCAL_LLM_ERROR"] = None


@app.context_processor
def inject_feature_flags():
    return {
        "show_drug_gpt_nav": app.config.get("SHOW_DRUG_GPT_NAV", False),
        "enable_drug_gpt": app.config.get("ENABLE_DRUG_GPT", False),
    }


@app.route("/healthz")
def healthz():
    return {"ok": True}, 200


@app.route("/backend-health")
def backend_health():
    mode = _normalized_backend_mode()
    status = {
        "ok": True,
        "backend_mode": mode,
        "backup_url_configured": bool(_vlismod_backup_url() or _randy_api_base_url()),
        "randy_token_configured": bool(_randy_api_token()),
        "strict_randy_disables_local_fallback": mode == "randy",
        "randy_health": None,
        "randy_db_health": None,
    }
    if randy_available():
        try:
            status["randy_health"] = randy_get("health")
        except RandyBackendError as exc:
            status["ok"] = False
            status["randy_health"] = {"error": str(exc)}
        try:
            status["randy_db_health"] = randy_get("db-health")
        except RandyBackendError as exc:
            status["ok"] = False
            status["randy_db_health"] = {"error": str(exc)}
    elif mode == "randy":
        status["ok"] = False
        status["randy_health"] = {"error": "RANDY API is not configured."}
    return jsonify(status)

# Define a consistent color palette for interaction types
INTERACTION_COLORS = {
    'clash': '#FF6347',           # Tomato (Red)
    'covalent': '#4682B4',        # Steel Blue
    'vdw_clash': '#FF4500',       # Orange Red
    'vdw': '#2196F3',             # Blue
    'proximal': '#8A2BE2',        # Blue Violet 
    'hbond': '#32CD32',           # Lime Green
    'weak_hbond': '#ADFF2F',      # Green Yellow
    'xbond': '#FF8C00',           # Dark Orange
    'ionic': '#00CED1',           # Dark Turquoise
    'metal': '#FFC107',        # Yellow
    'aromatic': '#DA70D6',        # Orchid
    'hydrophobic': '#2E8B57',     # Sea Green
    'carbonyl': '#D2691E',        # Chocolate
    'polar':    '#FFD700',        # Gold
    'weak_polar': '#9C27B0',      # Purple
    'CARBONPI': '#4CAF50',           # Green
    'CATIONPI': '#FFB6C1',        # Light Pink
    'DONORPI': '#BA55D3',         # Medium Orchid
    'HALOGENPI': '#8B0000',       # Dark Red
    'METSULPHURPI': '#FF9800'     # Orange
}


# Function to load the virus list from the text file
def load_virus_list(file_path):
    virus_list = []
    if os.path.exists(file_path):
        with open(file_path, 'r') as file:
            virus_list = [line.strip() for line in file.readlines()]
    return virus_list

# Function to dynamically generate green shades based on the number of viruses
def generate_green_shades(num_colors):
    base_hue = 120  # Hue for green
    saturation = 80  # Saturation for vivid colors
    start_lightness = 20  # Starting lightness percentage
    lightness_step = 2  # Step size for increasing lightness

    colors = []
    for i in range(num_colors):
        lightness = start_lightness + (i * lightness_step)
        if lightness > 90:  # Ensure lightness doesn't exceed the max value
            lightness = 90
        color = f'hsl({base_hue}, {saturation}%, {lightness}%)'
        colors.append(color)
    
    return colors

# Load the virus list from the text file
virus_file_path = os.path.join('static', 'virus_list', 'virus_list.txt')
virus_list = load_virus_list(virus_file_path)

# Generate green shades dynamically based on the number of viruses
colors = generate_green_shades(len(virus_list))

# Create a color map for the viruses
virusColorMap = {}
for i, virus in enumerate(virus_list):
    virusColorMap[virus] = colors[i]




@app.errorhandler(404)
def page_not_found(e):
    return render_template("404.html"), 404



@app.route('/get_virus_color_map', methods=['GET'])
def get_virus_color_map():
    return jsonify(virusColorMap)


CHARTS_DIR = 'static/charts/'
CHART_GENERATION_LOCK = threading.Lock()


@app.route("/", endpoint="index")
def index_page():
    return render_template("index.html")

@app.route('/about')
def about_page():
    return render_template('about.html')

@app.route('/contact')
def contact_page():
    return render_template('contact.html')

@app.route('/use-cases')
def use_cases_page():
    return render_template('use_cases.html')

@app.route('/viral-protac-design')
def viral_protac_design_page():
    return render_template('viral_protac_design.html')

@app.route('/viral-drug-targets')
def viral_drug_targets_page():
    return render_template('viral_drug_targets.html')

@app.route('/in-silico-virology-tools')
def in_silico_virology_tools_page():
    return render_template('in_silico_virology_tools.html')

@app.route('/methods')
def methods_page():
    return render_template('methods.html')

@app.route('/faq')
def faq_page():
    return render_template('faq.html')

@app.route('/citation')
def citation_page():
    return render_template('citation.html')

@app.route('/generate_ligand_images', methods=['POST'])
def generate_ligand_images():
    virus_name = str(request.form.get('virus') or '').strip()
    pdb_code = str(request.form.get('pdb_code') or '').strip()
    ligand_name = str(request.form.get('ligand') or '').strip()
    selected_chain = str(request.form.get('chain') or '').strip() or None

    if not virus_name or not pdb_code or not ligand_name:
        return jsonify({"error": "Missing required ligand image parameters."}), 400

    mode = _normalized_backend_mode()
    remote_payload = None

    if mode == "randy":
        try:
            remote_payload = randy_post(
                "ligand-images-data",
                json={
                    "virus_name": virus_name,
                    "pdb_code": pdb_code,
                    "ligand_name": ligand_name,
                    "chain": selected_chain,
                },
            )
        except RandyBackendError as exc:
            return jsonify({"error": str(exc)}), exc.status_code
    elif mode == "auto" and randy_available():
        try:
            remote_payload = randy_post(
                "ligand-images-data",
                json={
                    "virus_name": virus_name,
                    "pdb_code": pdb_code,
                    "ligand_name": ligand_name,
                    "chain": selected_chain,
                },
            )
        except RandyBackendError:
            logging.warning(
                "Falling back to local V-LiSEMOD database for ligand image generation: %s / %s / %s",
                virus_name,
                pdb_code,
                ligand_name,
            )

    if remote_payload is not None:
        smiles_data = [
            (
                row.get("virus_name"),
                row.get("pdb_id"),
                row.get("ligand"),
                row.get("smiles"),
            )
            for row in remote_payload.get("smiles_data", [])
        ]
        chain_residue_data = [
            (row.get("chain"), row.get("ligand_id"))
            for row in remote_payload.get("chain_residue_data", [])
        ]
        solvent_exposed_atom_map = remote_payload.get("solvent_exposed_atom_map") or {}
    else:
        try:
            with _connect_local_db(LIGAND_IMAGE_REQUIRED_TABLES) as conn:
                cursor = conn.cursor()
                cursor.execute(
                    '''
                    SELECT DISTINCT virus_name, pdb_id, ligand, smiles
                    FROM Functional_GROUPED
                    WHERE virus_name = ? AND pdb_id = ? AND ligand = ?
                    ''',
                    (virus_name, pdb_code, ligand_name),
                )
                smiles_data = cursor.fetchall()
                cursor.execute(
                    '''
                    SELECT DISTINCT chain, ligand_id
                    FROM Ligand_Arp_Diagram
                    WHERE virus_name = ? AND pdb_id = ? AND ligand = ?
                    ''',
                    (virus_name, pdb_code, ligand_name),
                )
                chain_residue_data = cursor.fetchall()
        except RandyBackendError as exc:
            status_code = 500 if mode == "local" else exc.status_code
            return jsonify({"error": str(exc)}), status_code
        solvent_exposed_atom_map = None

    if not smiles_data:
        return "No SMILES data found for the selected virus, PDB code, and ligand.", 404

    effective_chain = selected_chain or next(
        (str(chain).strip() for chain, _ligand_id in chain_residue_data if chain),
        None,
    )

    images = generate_images_from_smiles(
        smiles_data,
        effective_chain,
        "static/ligand_images",
        solvent_exposed_atom_map=solvent_exposed_atom_map,
    )
    return render_template(
        'display_images.html',
        images=images,
        chain_residues=chain_residue_data,
        structure_url=None
    )




from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
import os


def generate_images_from_smiles(smiles_data, selected_chain, output_folder, solvent_exposed_atom_map=None):
    images = []
    os.makedirs(output_folder, exist_ok=True)

    for virus_name, pdb_id, ligand_code, smiles in smiles_data:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            print(f"❌ Invalid SMILES: {ligand_code}")
            continue

        AllChem.Compute2DCoords(mol)

        # -------------------------
        # STANDARD SVG
        # -------------------------
        drawer = rdMolDraw2D.MolDraw2DSVG(600, 600)
        opts = drawer.drawOptions()
        opts.bondLineWidth = 2.0
        opts.addStereoAnnotation = True

        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()

        svg = drawer.GetDrawingText()

        # Inject ligand label (SVG-safe)
        svg = svg.replace(
            "</svg>",
            f"""
            <text x="300" y="580"
                  text-anchor="middle"
                  font-size="48"
                  font-family="Arial"
                  fill="black">{ligand_code}</text>
            </svg>
            """
        )

        svg_path = f"{output_folder}/{pdb_id}_{ligand_code}.svg"
        with open(svg_path, "w") as f:
            f.write(svg)

        # -------------------------
        # SOLVENT-EXPOSED SVG
        # -------------------------
        solvent_svg_path = f"{output_folder}/{pdb_id}_{ligand_code}_solvent_exposed.svg"
        if solvent_exposed_atom_map is not None:
            solvent_key = f"{pdb_id}|{ligand_code}|{selected_chain}"
            highlight_solvent_exposed_atoms_from_indices(
                mol,
                solvent_exposed_atom_map.get(solvent_key, []),
                solvent_svg_path,
            )
        else:
            highlight_solvent_exposed_atoms(
                mol,
                virus_name,
                pdb_id,
                ligand_code,
                selected_chain,
                solvent_svg_path
            )

        images.append({
            "path": svg_path,
            "solvent_exposed_path": solvent_svg_path,
            "virus_name": virus_name,
            "pdb_id": pdb_id,
            "ligand_code": ligand_code,
            "chain": selected_chain,
            "filename": os.path.basename(svg_path)
        })

        print(f"✅ SVG saved: {svg_path}")
        print(f"✅ Solvent SVG saved: {solvent_svg_path}")

    return images


# Function to write PyMOL script
def _pymol_ligand_atom_selection(ligand_name, ligand_chain, ligand_residue_id, atom):
    """Build a stable PyMOL selection for a ligand atom.

    V2 data records use mmCIF atom_site IDs, while ``fetch ..., type=pdb``
    assigns PDB serial IDs.  Atom serials therefore must never be used to
    select a rebuilt ligand atom in PyMOL.  The author atom name together
    with the resolved ligand residue identity is stable across the formats.
    """
    atom_name = str(atom.get("exact_atom") or "").strip()
    if not atom_name:
        return "none"
    parts = [f"resn {ligand_name}", f"chain {ligand_chain}"]
    if ligand_residue_id not in (None, ""):
        parts.append(f"resi {ligand_residue_id}")
    parts.append(f"name {atom_name}")
    return " and ".join(parts)


def _pymol_protein_atom_selection(atom):
    """Select a pocket atom by residue identity, not CIF/PDB serial number."""
    chain = str(atom.get("residue_chain") or "").strip()
    residue_id = str(atom.get("residue_number") or "").strip()
    atom_name = str(atom.get("residue_atom") or "").strip()
    if not (chain and residue_id and atom_name):
        return "none"
    return f"polymer and chain {chain} and resi {residue_id} and name {atom_name}"


def write_pymol_script(pdb_code, ligand_name, ligand_chain, ligand_residue_id, options):
    output_dir = './pml_sessions'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    output_script = f"{output_dir}/{pdb_code}.pml"
    
    with open(output_script, 'w') as script:
        script.write(f"fetch {pdb_code}, async=0, type=pdb\n")
        script.write(f"select {ligand_name}_Ligand, resn {ligand_name} and chain {ligand_chain}\n")
        script.write(f"create Ligand_Object, {ligand_name}_Ligand\n")
        script.write(f"show sticks, Ligand_Object\n")
        script.write(f"color gray, Ligand_Object and elem C\n")
        script.write(f"color red, Ligand_Object and elem O\n")
        script.write(f"color blue, Ligand_Object and elem N\n")
        script.write(f"color yellow, Ligand_Object and elem S\n")
        
        # Functional groups
        if options.get('functional_groups'):
            functional_groups = options['functional_groups']
            for fg_name, atoms in functional_groups.items():
                atom_selection = " or ".join(
                    _pymol_ligand_atom_selection(ligand_name, ligand_chain, ligand_residue_id, atom)
                    for atom in atoms
                )
                script.write(f"select {fg_name}, ({atom_selection})\n")
                script.write(f"create {fg_name}_Object, {fg_name}\n")
                script.write(f"show sticks, {fg_name}_Object\n")
                script.write(f"color magenta, {fg_name}_Object\n")
        
        # Binding pocket
        if options.get('binding_pocket'):
            binding_pocket_atoms = options['binding_pocket']
            script.write("create Binding_Pocket, none\n")
            for atom in binding_pocket_atoms:
                script.write(f"select temp, {_pymol_protein_atom_selection(atom)}\n")
                script.write("create Binding_Pocket, Binding_Pocket or temp\n")
                script.write("delete temp\n")
            script.write("show surface, Binding_Pocket\n")
            script.write("color yellow, Binding_Pocket\n")
            script.write("set transparency, 0.1, Binding_Pocket\n")
        
        # Distal atoms
        if options.get('distal_atoms'):
            distal_atoms = options['distal_atoms']
            script.write("create Distal_Atoms, none\n")
            for atom in distal_atoms:
                script.write(f"select temp, {_pymol_ligand_atom_selection(ligand_name, ligand_chain, ligand_residue_id, atom)}\n")
                script.write("create Distal_Atoms, Distal_Atoms or temp\n")
                script.write("delete temp\n")
            script.write("show spheres, Distal_Atoms\n")
            script.write("color blue, Distal_Atoms\n")
        
        # Solvent exposed atoms
        if options.get('solvent_exposed_atoms'):
            solvent_exposed_atoms = options['solvent_exposed_atoms']
            logging.debug(f"Solvent-exposed atoms (SASA): {solvent_exposed_atoms}")
            script.write("create Solvent_Exposed_Atoms, none\n")
            for atom in solvent_exposed_atoms:
                script.write(f"select temp, {_pymol_ligand_atom_selection(ligand_name, ligand_chain, ligand_residue_id, atom)}\n")
                script.write("create Solvent_Exposed_Atoms, Solvent_Exposed_Atoms or temp\n")
                script.write("delete temp\n")
            script.write("show spheres, Solvent_Exposed_Atoms\n")
            script.write("color firebrick, Solvent_Exposed_Atoms\n")

        # Hydrated atoms
        if options.get('hydrated_atoms'):
            hydrated_atoms = options['hydrated_atoms']
            script.write("create Hydrated_Atoms, none\n")
            for atom in hydrated_atoms:
                script.write(f"select temp, {_pymol_ligand_atom_selection(ligand_name, ligand_chain, ligand_residue_id, atom)}\n")
                script.write("create Hydrated_Atoms, Hydrated_Atoms or temp\n")
                script.write("delete temp\n")
            script.write("show sticks, Hydrated_Atoms\n")
            script.write("color cyan, Hydrated_Atoms\n")

        # Rupley SASA atoms
        if options.get('rupley_sasa'):
            rupley_sasa_atoms = options['rupley_sasa']
            script.write("create RUPLEY_SASA, none\n")
            for atom in rupley_sasa_atoms:
                script.write(f"select temp, {_pymol_ligand_atom_selection(ligand_name, ligand_chain, ligand_residue_id, atom)}\n")
                script.write("create RUPLEY_SASA, RUPLEY_SASA or temp\n")
                script.write("delete temp\n")
            script.write("show spheres, RUPLEY_SASA\n")
            script.write("color purple, RUPLEY_SASA\n")
        
        script.write("show cartoon, all\n")
        script.write("hide lines, all\n")
        session_file = os.path.join(output_dir, f"{pdb_code}.pse")
        script.write(f"save {session_file}\n")
    
    return output_script

@app.route('/download_image/<filename>')
def download_image(filename):
    return send_from_directory('static/ligand_images', filename)

@app.route('/download_SASA_image/<filename>')
def download_SASA_image(filename):
    # Ensure the filename refers to the SASA image file
    return send_from_directory('static/ligand_images', filename)


@app.route('/generate_pymol_session', methods=['POST'])
def generate_pymol_session():
    pdb_code = str(request.form.get('pdb_code') or '').strip()
    ligand_name = str(request.form.get('ligand') or '').strip()
    requested_chain = str(request.form.get('chain') or '').strip() or None
    requested_ligand_instance_id = str(request.form.get('ligand_instance_id') or '').strip() or None

    if not pdb_code or not ligand_name:
        return jsonify({"error": "Missing required PyMOL session parameters."}), 400

    option_flags = {
        'functional_groups': _form_flag('functional_groups'),
        'binding_pocket': _form_flag('binding_pocket'),
        'distal_atoms': _form_flag('distal_atoms'),
        'solvent_exposed_atoms': _form_flag('solvent_exposed_atoms'),
        'hydrated_atoms': _form_flag('hydrated_atoms'),
        'rupley_sasa': _form_flag('rupley_sasa'),
    }

    mode = _normalized_backend_mode()
    remote_payload = None

    if mode == "randy":
        try:
            remote_payload = randy_post(
                "pymol-session-data",
                json={
                    "pdb_code": pdb_code,
                    "ligand_name": ligand_name,
                    "chain": requested_chain,
                    "ligand_instance_id": requested_ligand_instance_id,
                    "options": option_flags,
                },
            )
        except RandyBackendError as exc:
            return jsonify({"error": str(exc)}), exc.status_code
    elif mode == "auto" and randy_available():
        try:
            remote_payload = randy_post(
                "pymol-session-data",
                json={
                    "pdb_code": pdb_code,
                    "ligand_name": ligand_name,
                    "chain": requested_chain,
                    "ligand_instance_id": requested_ligand_instance_id,
                    "options": option_flags,
                },
            )
        except RandyBackendError:
            logging.warning(
                "Falling back to local V-LiSEMOD database for PyMOL session generation: %s / %s",
                pdb_code,
                ligand_name,
            )

    if remote_payload is not None:
        ligand_chain = remote_payload.get("ligand_chain")
        if not ligand_chain:
            return jsonify({"error": "RANDY response did not include a ligand chain."}), 502

        ligand_residue_id = remote_payload.get("ligand_residue_id")
        options = {}
        if option_flags['functional_groups']:
            options['functional_groups'] = {
                group_name: list(atoms)
                for group_name, atoms in (remote_payload.get("functional_groups") or {}).items()
            }
        if option_flags['binding_pocket']:
            options['binding_pocket'] = list(remote_payload.get("binding_pocket", []))
        if option_flags['distal_atoms']:
            options['distal_atoms'] = list(remote_payload.get("distal_atoms", []))
        if option_flags['solvent_exposed_atoms']:
            options['solvent_exposed_atoms'] = list(remote_payload.get("solvent_exposed_atoms", []))
        if option_flags['hydrated_atoms']:
            options['hydrated_atoms'] = list(remote_payload.get("hydrated_atoms", []))
        if option_flags['rupley_sasa']:
            options['rupley_sasa'] = list(remote_payload.get("rupley_sasa", []))
    else:
        try:
            with _connect_local_db(PYMOL_REQUIRED_TABLES) as conn:
                cursor = conn.cursor()
                cursor.execute(
                    '''
                    SELECT DISTINCT chain
                    FROM ligand_atoms
                    WHERE pdb_id = ? AND ligand = ?
                    ORDER BY chain
                    ''',
                    (pdb_code, ligand_name),
                )
                chains = [row[0] for row in cursor.fetchall() if row[0]]

                if requested_chain:
                    if requested_chain not in chains:
                        return jsonify({"error": "Specified ligand chain was not found."}), 404
                    ligand_chain = requested_chain
                else:
                    if not chains:
                        return jsonify({"error": "No ligand chain found for the selected PDB and ligand."}), 404
                    if len(chains) > 1:
                        return jsonify({"error": "Multiple chains found; please specify the chain."}), 400
                    ligand_chain = chains[0]

                options = {}
                if option_flags['functional_groups']:
                    cursor.execute(
                        '''
                        SELECT functional_group, atom_id, exact_atom, atom_type, chain
                        FROM Functional_Group_Atoms
                        WHERE virus_name = (
                              SELECT virus_name
                              FROM ligand_atoms
                              WHERE pdb_id = ? AND ligand = ? AND chain = ?
                              LIMIT 1
                          )
                          AND pdb_id = ?
                          AND ligand = ?
                          AND chain = ?
                        ''',
                        (pdb_code, ligand_name, ligand_chain, pdb_code, ligand_name, ligand_chain),
                    )
                    rows = cursor.fetchall()
                    functional_groups = {}
                    for fg_name, atom_id, exact_atom, atom_type, _chain in rows:
                        functional_groups.setdefault(fg_name, []).append((atom_id, exact_atom, atom_type))
                    options['functional_groups'] = functional_groups

                if option_flags['binding_pocket']:
                    cursor.execute(
                        '''
                        SELECT residue_chain, residue_number
                        FROM receptor_binding_pocket
                        WHERE pdb_id = ?
                        ''',
                        (pdb_code,),
                    )
                    options['binding_pocket'] = cursor.fetchall()

                if option_flags['distal_atoms']:
                    cursor.execute(
                        '''
                        SELECT chain, atom_id
                        FROM distal_atoms
                        WHERE pdb_id = ? AND ligand = ?
                        ''',
                        (pdb_code, ligand_name),
                    )
                    options['distal_atoms'] = cursor.fetchall()

                if option_flags['solvent_exposed_atoms']:
                    cursor.execute(
                        '''
                        SELECT DISTINCT atom_id, chain
                        FROM solvent_exposed_atoms
                        WHERE pdb_id = ? AND ligand = ? AND (chain = ? OR ? IS NULL)
                        ''',
                        (pdb_code, ligand_name, ligand_chain, ligand_chain),
                    )
                    options['solvent_exposed_atoms'] = cursor.fetchall()

                if option_flags['hydrated_atoms']:
                    cursor.execute(
                        '''
                        SELECT chain, atom_id
                        FROM ligand_atoms
                        WHERE pdb_id = ? AND ligand = ? AND (chain = ? OR ? IS NULL)
                        ''',
                        (pdb_code, ligand_name, ligand_chain, ligand_chain),
                    )
                    options['hydrated_atoms'] = cursor.fetchall()

                if option_flags['rupley_sasa']:
                    cursor.execute(
                        '''
                        SELECT atom_id, chain
                        FROM RUPLEY_SASA_DATA
                        WHERE pdb_id = ? AND ligand = ? AND (chain = ? OR ? IS NULL)
                        ''',
                        (pdb_code, ligand_name, ligand_chain, ligand_chain),
                    )
                    options['rupley_sasa'] = cursor.fetchall()
        except RandyBackendError as exc:
            status_code = 500 if mode == "local" else exc.status_code
            return jsonify({"error": str(exc)}), status_code

    pymol_script_path = write_pymol_script(
        pdb_code, ligand_name, ligand_chain, ligand_residue_id if remote_payload is not None else None, options
    )
    return send_file(pymol_script_path, as_attachment=True)


def _local_get_viruses_payload():
    conn = sqlite3.connect(str(LOCAL_DB_PATH))
    cursor = conn.cursor()
    cursor.execute('SELECT DISTINCT virus_name FROM ligand_atoms')
    viruses = [row[0] for row in cursor.fetchall()]
    conn.close()
    return {'viruses': viruses}


def _local_get_pdb_codes_payload(virus_name):
    conn = sqlite3.connect(str(LOCAL_DB_PATH))
    cursor = conn.cursor()
    cursor.execute(
        '''
        SELECT DISTINCT pdb_id
        FROM ligand_Atoms_Smiles
        WHERE virus_name = ?
        ''',
        (virus_name,),
    )
    pdb_codes = [row[0] for row in cursor.fetchall()]
    conn.close()
    return {'pdb_codes': pdb_codes}


def _local_get_ligands_payload(pdb_code):
    conn = sqlite3.connect(str(LOCAL_DB_PATH))
    cursor = conn.cursor()
    cursor.execute(
        '''
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
        ''',
        (pdb_code,),
    )
    ligands = [{'ligand': row[0], 'chain': row[1], 'has_smiles': row[2]} for row in cursor.fetchall()]
    conn.close()
    return {'ligands': ligands}


def _local_check_functional_groups_payload(pdb_code):
    conn = sqlite3.connect(str(LOCAL_DB_PATH))
    cursor = conn.cursor()
    cursor.execute(
        '''
        SELECT COUNT(*)
        FROM Functional_Group_Atoms
        WHERE pdb_id = ?
        ''',
        (pdb_code,),
    )
    count = cursor.fetchone()[0]
    conn.close()
    return {'has_functional_groups': count > 0}


def _local_get_ligands_list_payload():
    conn = sqlite3.connect(str(LOCAL_DB_PATH))
    cursor = conn.cursor()
    cursor.execute('SELECT DISTINCT ligand FROM Ligand_Atoms_Smiles')
    ligands = [row[0] for row in cursor.fetchall()]
    conn.close()
    return {'ligands': ligands}


def _local_get_viruses_by_ligand_payload(ligand_code):
    conn = sqlite3.connect(str(LOCAL_DB_PATH))
    cursor = conn.cursor()
    cursor.execute('SELECT DISTINCT Virus_Name FROM Arpeggio_Contacts_Data WHERE Ligand = ?', (ligand_code,))
    viruses = [row[0] for row in cursor.fetchall()]
    conn.close()
    return {'viruses': viruses}


def _local_get_pdb_residue_by_ligand_payload(ligand_code):
    with _connect_local_db(LIGAND_INTERACTION_REQUIRED_TABLES) as conn:
        rows = conn.execute(
            '''
            SELECT DISTINCT
                s.entry_id AS pdb_id,
                i.auth_asym_id AS chain,
                i.auth_seq_id AS ligand_id,
                i.insertion_code_normalized AS insertion_code,
                i.deposited_model_num AS model_id,
                i.ligand_instance_id
            FROM ligand_instances AS i
            JOIN structures AS s ON s.structure_id = i.structure_id
            WHERE i.label_comp_id = ?
              AND i.curation_status = 'included'
              AND EXISTS (
                  SELECT 1
                  FROM arpeggio_raw_contact_labels AS r
                  WHERE r.ligand_instance_id = i.ligand_instance_id
                    AND r.filter_class = 'raw_environment'
                    AND r.run_id = (
                        SELECT MAX(ar.run_id)
                        FROM ligand_arpeggio_runs AS ar
                        WHERE ar.ligand_instance_id = i.ligand_instance_id
                          AND ar.status = 'completed'
                    )
              )
            ORDER BY s.entry_id, i.deposited_model_num, i.auth_asym_id,
                     i.auth_seq_id, i.insertion_code_normalized, i.ligand_instance_id
            ''',
            (ligand_code,),
        ).fetchall()
    pairs = [
        {
            'pdb_id': row[0],
            'chain': row[1],
            'ligand_id': row[2],
            'insertion_code': row[3],
            'model_id': row[4],
            'ligand_instance_id': row[5],
        }
        for row in rows
    ]
    return {'pairs': pairs}


def _local_get_sasa_chains_payload(pdb_code, ligand_name):
    conn = sqlite3.connect(str(LOCAL_DB_PATH))
    cursor = conn.cursor()
    cursor.execute(
        '''
        SELECT DISTINCT atom_id, chain
        FROM solvent_exposed_atoms
        WHERE pdb_id = ? AND ligand = ?
        ''',
        (pdb_code, ligand_name),
    )
    sasa_chains = cursor.fetchall()
    conn.close()
    return sasa_chains


def _local_get_pdb_mapping_payload(ligand_code):
    """Return mapped, included ligand occurrences without collapsing identity."""
    with _connect_local_db(COMPARE_LIGAND_REQUIRED_TABLES) as conn:
        rows = conn.execute(
            '''
            WITH selected_occurrences AS (
                SELECT
                    i.ligand_instance_id,
                    s.entry_id AS pdb_id,
                    i.deposited_model_num AS model_id,
                    i.auth_asym_id AS chain,
                    i.auth_seq_id AS ligand_id,
                    i.insertion_code_normalized AS insertion_code,
                    COALESCE(sc.virus_label, 'Unknown') AS virus_name,
                    i.label_comp_id AS ligand
                FROM ligand_instances AS i
                JOIN structures AS s ON s.structure_id = i.structure_id
                LEFT JOIN structure_classifications AS sc ON sc.structure_id = i.structure_id
                WHERE i.label_comp_id = ?
                  AND i.curation_status = 'included'
            ), latest_mapping_runs AS (
                SELECT m.ligand_instance_id, MAX(m.run_id) AS run_id
                FROM ligand_smiles_atom_mapping AS m
                JOIN selected_occurrences AS so
                  ON so.ligand_instance_id = m.ligand_instance_id
                WHERE m.method_version = 'legacy_mcs_etkdg_uff_cif_v2.5'
                GROUP BY m.ligand_instance_id
            )
            SELECT so.*
            FROM selected_occurrences AS so
            JOIN latest_mapping_runs AS lmr
              ON lmr.ligand_instance_id = so.ligand_instance_id
            ORDER BY so.pdb_id, so.model_id, so.chain, so.ligand_id,
                     so.insertion_code, so.ligand_instance_id
            ''',
            (ligand_code,),
        ).fetchall()

    pdb_mapping = {}
    for row in rows:
        ligand_instance_id, pdb_id, model_id, chain, ligand_id, insertion_code, virus_name, ligand = row
        unique_key = str(ligand_instance_id)
        pdb_mapping[unique_key] = {
            'ligand_instance_id': ligand_instance_id,
            'pdb_id': pdb_id,
            'model_id': model_id,
            'chain': chain,
            'ligand_id': ligand_id,
            'insertion_code': insertion_code,
            'virus_name': virus_name,
            'ligand': ligand,
            'legacy_key': f"{pdb_id}-{ligand_id}-{chain}",
        }
    return {'pdb_mapping': pdb_mapping}

@app.route('/get_viruses')
def get_viruses():
    return _dispatch_supported_lookup(
        'viruses',
        local_loader=lambda: jsonify(_local_get_viruses_payload()),
    )

@app.route('/get_pdb_codes/<virus_name>')
def get_pdb_codes(virus_name):
    return _dispatch_supported_lookup(
        'pdb-codes',
        params={'virus_name': virus_name},
        local_loader=lambda: jsonify(_local_get_pdb_codes_payload(virus_name)),
    )



@app.route('/get_ligands/<pdb_code>')
def get_ligands(pdb_code):
    return _dispatch_supported_lookup(
        'ligands',
        params={'pdb_code': pdb_code},
        local_loader=lambda: jsonify(_local_get_ligands_payload(pdb_code)),
    )




@app.route('/check_functional_groups/<pdb_code>')
def check_functional_groups(pdb_code):
    return _dispatch_supported_lookup(
        'functional-groups/check',
        params={'pdb_code': pdb_code},
        local_loader=lambda: jsonify(_local_check_functional_groups_payload(pdb_code)),
    )


@app.route('/coming-soon')
def coming_soon():
    return render_template('coming-soon.html')







@app.route('/view_ligand_3d/<ligand_code>/<pdb_id>')
def view_ligand_3d(ligand_code, pdb_id):
    return render_template('ligand_3d_viewer.html', ligand_code=ligand_code, pdb_id=pdb_id)









@app.route('/get_ligands_list', methods=['GET'])
def get_ligands_list():
    return _dispatch_supported_lookup(
        'ligands/list',
        local_loader=lambda: jsonify(_local_get_ligands_list_payload()),
    )



@app.route('/get_viruses_by_ligand/<ligand_code>')
def get_viruses_by_ligand(ligand_code):
    return _dispatch_supported_lookup(
        'viruses/by-ligand',
        params={'ligand_code': ligand_code},
        local_loader=lambda: jsonify(_local_get_viruses_by_ligand_payload(ligand_code)),
    )





# Function to clean and preprocess interactions
def preprocess_interactions(df):
    # Filter out 'proximal' interactions
    df_filtered = df[df['Interaction'] != 'proximal']

    # Merge 'weak polar' and 'polar' into 'polar'
    df_filtered['Interaction'] = df_filtered['Interaction'].replace({'weak_polar': 'polar'})

    # Merge 'vdw' and 'vdw clash' into 'vdw'
    df_filtered['Interaction'] = df_filtered['Interaction'].replace({'vdw_clash': 'vdw'})
    
    # Merge 'vdw' and 'vdw clash' into 'vdw'
    df_filtered['Interaction'] = df_filtered['Interaction'].replace({'weak_hbond': 'hbond'})
    

    return df_filtered


def filter_valid_atom_ids(df):
    # Check if either 'Atom_ID' or 'atom_id' exists in the DataFrame
    if 'atom_id' in df.columns:
        column_name = 'atom_id'
    elif 'atom_id' in df.columns:
        column_name = 'atom_id'
    else:
        print("Warning: Neither 'atom_id' nor 'atom_id' column is present in the DataFrame")
        return df  # Return the unfiltered DataFrame if neither column is present
    
    # Filter based on the found column name
    df_filtered = df[df[column_name].notna() & (df[column_name] != "N/A")]
    df_filtered[column_name] = df_filtered[column_name].astype(int)
    return df_filtered





@app.route('/ligand_indexer')
def ligand_indexer():
    return render_template('ligand_query.html')


@app.route('/get_pdb_residue_by_ligand/<ligand_code>')
def get_pdb_residue_by_ligand(ligand_code):
    return _dispatch_supported_lookup(
        'pdb-residues/by-ligand',
        params={'ligand_code': ligand_code},
        local_loader=lambda: jsonify(_local_get_pdb_residue_by_ligand_payload(ligand_code)),
    )


# Function to generate the DataFrame from the SQLite database
def get_interaction_data(pdb_id, ligand, ligand_id, chain, ligand_instance_id=None):
    mode = _normalized_backend_mode()

    if mode == "randy":
        payload = randy_get(
            "interaction-records",
            params={"pdb_id": pdb_id, "ligand": ligand, "ligand_id": ligand_id, "chain": chain, "ligand_instance_id": ligand_instance_id},
        )
        return _records_to_interaction_dataframe(payload.get("records", []))

    if mode == "auto" and randy_available():
        try:
            payload = randy_get(
                "interaction-records",
                params={"pdb_id": pdb_id, "ligand": ligand, "ligand_id": ligand_id, "chain": chain, "ligand_instance_id": ligand_instance_id},
            )
            return _records_to_interaction_dataframe(payload.get("records", []))
        except RandyBackendError:
            logging.warning(
                "Falling back to local interaction data for %s / %s / %s / %s",
                pdb_id,
                ligand,
                ligand_id,
                chain,
            )

    if ligand_instance_id in (None, ""):
        raise ValueError("ligand_instance_id is required for occurrence-resolved interaction lookup.")
    payload = _local_interaction_records_payload(
        pdb_id, ligand, ligand_id, chain, int(ligand_instance_id)
    )
    return _records_to_interaction_dataframe(payload.get("records", []))

# Function to clean and preprocess interactions (now applies to 'Contact')
def preprocess_interactions(df):
    # Filter out 'proximal' interactions
    df_filtered = df[df['Contact'] != 'proximal']

    # Merge 'weak polar' and 'polar' into 'polar'
    df_filtered.loc[:, 'Contact'] = df_filtered['Contact'].replace({'weak_polar': 'polar'})

    # Merge 'vdw' and 'vdw clash' into 'vdw'
    df_filtered.loc[:, 'Contact'] = df_filtered['Contact'].replace({'vdw_clash': 'vdw'})

    return df_filtered

# Function to remove rows without valid atom_id and ensure atom_ids are integers
def filter_valid_atom_ids(df):
    # Check if either 'atom_id' or 'atom_id' exists in the DataFrame
    if 'atom_id' in df.columns:
        column_name = 'atom_id'
    elif 'atom_id' in df.columns:
        column_name = 'atom_id'
    else:
        print("Warning: Neither 'atom_id' nor 'atom_id' column is present in the DataFrame")
        return df  # Return the unfiltered DataFrame if neither column is present
    
    # Filter rows where the column is not 'N/A' or missing
    df_filtered = df[df[column_name].notna() & (df[column_name] != "N/A")]
    
    # Ensure the column is treated as an integer
    df_filtered[column_name] = df_filtered[column_name].astype(int)
    
    return df_filtered


# Function to generate a pie chart for interaction types, with filtering
def plot_pie_chart(df, output_file, pdb_id, ligand_code, small_threshold=5):
    # Filter and preprocess interactions using the existing preprocessing logic
    df_filtered = preprocess_interactions(df)  # Cleaned and processed interaction data
    
    interaction_counts = df_filtered['Contact'].value_counts()

    # Separate out small interaction types below the threshold
    small_interactions = interaction_counts[interaction_counts < small_threshold]
    main_interactions = interaction_counts[interaction_counts >= small_threshold]

    # Only create an "Others" category if there are two or more small interactions
    if len(small_interactions) > 1:
        others_sum = small_interactions.sum()
        main_interactions['Others'] = others_sum

    # Otherwise, keep the small category as its own label
    elif len(small_interactions) == 1:
        main_interactions = pd.concat([main_interactions, small_interactions])

    # Get the color mapping for interactions (use INTERACTION_COLORS defined earlier)
    colors = [INTERACTION_COLORS.get(interaction, '#808080') for interaction in main_interactions.index]

    # Use explode to separate the smallest slices slightly (for better readability)
    explode = [0.1 if interaction == 'Others' else 0 for interaction in main_interactions.index]

    # Plot the pie chart as a donut
    plt.figure(figsize=(8, 8))
    plt.pie(
        main_interactions, 
        labels=main_interactions.index, 
        autopct='%1.1f%%', 
        startangle=90,
        colors=colors, 
        wedgeprops={'edgecolor': 'white', 'linewidth': 1.5},
        explode=explode
    )

    # Draw a circle at the center to turn the pie into a donut
    center_circle = plt.Circle((0, 0), 0.70, fc='white')
    plt.gca().add_artist(center_circle)

    plt.title(f'Distribution of Interaction Types for PDB: {pdb_id} & Ligand: {ligand_code}')
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()





# Function to generate a bar chart for interaction types
def plot_bar_chart(df, output_file, pdb_id, ligand_code):
    # Filter and preprocess interactions using the existing preprocessing logic
    df_filtered = preprocess_interactions(df)  # Cleaned and processed interaction data
    interaction_counts = df_filtered['Contact'].value_counts()

    # Get the color mapping for interactions
    colors = [INTERACTION_COLORS[interaction] for interaction in interaction_counts.index]

    # Plot the bar chart
    plt.figure(figsize=(10, 6))
    plt.bar(interaction_counts.index, interaction_counts.values, color=colors)
    
    plt.title(f'Count of Interaction Types for PDB: {pdb_id} & Ligand: {ligand_code} ')
    plt.xlabel('Interaction Type')
    plt.ylabel('Count')
    plt.xticks(rotation=45)
# Add counts at the top of each bar
    for index, value in enumerate(interaction_counts.values):
        plt.text(index, value, str(value), ha='center', va='bottom')

    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()


# Function to generate a scatter plot to show distances of interactions per type
def plot_scatter_chart(df, output_file, pdb_id, ligand_code):
    # Preprocess the dataframe to filter out unwanted interactions
    df_filtered = preprocess_interactions(df)

    # Plot the scatter chart with the filtered data
    plt.figure(figsize=(10, 6))
    sns.scatterplot(x='Distance', y='Contact', data=df_filtered, hue='Contact', palette=INTERACTION_COLORS)
    
    plt.title(f'Interaction Distance by Type for PDB: {pdb_id} & Ligand: {ligand_code} ')
    plt.xlabel('Distance')
    plt.ylabel('Interaction Type')
    
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()


# Function to generate a bar chart for interactions per atom, excluding 'proximal'
def plot_interactions_per_atom(df, output_file, pdb_id, ligand_code, exclude_proximal=True):
    # If exclude_proximal is True, filter out 'proximal' interactions
    if exclude_proximal:
        df = df[df['Contact'] != 'proximal'].copy()
    else:
        df = df.copy()

    # Merge 'weak_polar' into 'polar' and 'vdw_clash' into 'vdw'
    df.loc[:, 'Contact'] = df['Contact'].replace({'weak_polar': 'polar', 'vdw_clash': 'vdw', 'weak_hbond': 'hbond'})

    # Create a new column that combines atom_id and Ligand_Atom for labeling
    df.loc[:, 'Atom_Label'] = df['atom_id'].astype(str) + ' (' + df['exact_atom'].astype(str) + ')'

    # Filter out rows where Atom_Label contains 'nan'
    df = df.loc[df['Atom_Label'].notna() & (df['Atom_Label'] != 'nan')].copy()

    # Group by Atom_Label and Contact and count occurrences
    df_atom_interactions = df.groupby(['Atom_Label', 'Contact']).size().unstack(fill_value=0)

    # Define colors in the same order as interaction columns
    colors = [INTERACTION_COLORS.get(col, '#808080') for col in df_atom_interactions.columns]

    # Plotting
    plt.figure(figsize=(12, 6))
    ax = df_atom_interactions.plot(kind='bar', stacked=True, color=colors, width=0.8, figsize=(12, 6))

    # Set title and labels
    if exclude_proximal:
        plt.title(f'Interactions per Atom (Excluding Proximal) for PDB: {pdb_id} & Ligand: {ligand_code}')
    else:
        plt.title(f'Interactions per Atom (Including Proximal) for PDB: {pdb_id} & Ligand: {ligand_code}')

    plt.xlabel('Atom Number (Ligand Atom)')
    plt.ylabel('Number of Interactions')

    # Manually set the X-axis tick labels to be the actual Atom_Label values
    ax.set_xticks(range(len(df_atom_interactions.index)))
    ax.set_xticklabels(df_atom_interactions.index, rotation=45, ha='right')

    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()


# Update generate_charts to include this new plot
@app.route('/generate_charts', methods=['POST'])
def generate_charts():
    data = request.get_json(silent=True) or {}
    required_fields = ('pdb_id', 'ligand', 'ligand_id', 'chain', 'ligand_instance_id')
    missing = [field for field in required_fields if data.get(field) in (None, '')]
    if missing:
        return jsonify({
            'error': 'invalid_interaction_request',
            'message': f"Missing required interaction field(s): {', '.join(missing)}.",
        }), 400

    pdb_id = str(data['pdb_id']).strip()
    ligand = str(data['ligand']).strip()
    ligand_id = str(data['ligand_id']).strip()
    chain = str(data['chain']).strip()
    ligand_instance_id = data['ligand_instance_id']

    try:
        df = get_interaction_data(
            pdb_id, ligand, ligand_id, chain, ligand_instance_id
        )
    except RandyBackendError as exc:
        return jsonify({'error': str(exc)}), exc.status_code
    except (TypeError, ValueError) as exc:
        return jsonify({'error': 'invalid_interaction_request', 'message': str(exc)}), 400

    if df.empty:
        return jsonify({
            'error': 'no_interaction_data',
            'message': 'No recorded protein–ligand interactions are available for this selected ligand occurrence.',
        }), 404

    df_clean = filter_valid_atom_ids(preprocess_interactions(df).copy())
    if df_clean.empty:
        return jsonify({
            'error': 'no_interaction_data',
            'message': 'No non-proximal protein–ligand interactions are available for this selected ligand occurrence.',
        }), 422
    
    # Filter the 'proximal' interactions and ensure they have valid atom_ids
    df_proximal = df[df['Contact'] == 'proximal']
    df_proximal_clean = filter_valid_atom_ids(df_proximal)

    # Add 'proximal' interactions back to the dataset for specific plots
    df_with_proximal = pd.concat([df_clean, df_proximal_clean])

    os.makedirs(CHARTS_DIR, exist_ok=True)
    chart_request_id = uuid.uuid4().hex
    chart_specs = [
        ('pie_chart', plot_pie_chart, df_clean),
        ('bar_chart', plot_bar_chart, df_clean),
        ('scatter_chart', plot_scatter_chart, df_clean),
        ('interactions_per_atom_with_proximal', plot_interactions_per_atom, df_with_proximal),
        ('interactions_per_atom_without_proximal', plot_interactions_per_atom, df_clean),
    ]
    chart_filenames = [f"{chart_request_id}_{name}.png" for name, _, _ in chart_specs]

    # Matplotlib is process-global. Serializing only rendering prevents a
    # concurrent request from corrupting figures, while request-specific names
    # keep one user's charts from replacing another's.
    with CHART_GENERATION_LOCK:
        for (name, plotter, chart_df), filename in zip(chart_specs, chart_filenames):
            output_file = os.path.join(CHARTS_DIR, filename)
            if name.startswith('interactions_per_atom'):
                plotter(
                    chart_df,
                    output_file,
                    pdb_id,
                    ligand,
                    exclude_proximal=name.endswith('without_proximal'),
                )
            else:
                plotter(chart_df, output_file, pdb_id, ligand)

    return jsonify({
        'chart_urls': [
            f"{url_for('static', filename=f'charts/{filename}')}?v={chart_request_id}"
            for filename in chart_filenames
        ],
        'interaction_label_count': int(len(df.index)),
        'unique_partner_count': int(
            df[['residue_chain', 'residue_number', 'residue_atom']].drop_duplicates().shape[0]
        ),
        'ligand_instance_id': int(ligand_instance_id),
    })



@app.route('/get_sasa_chains/<pdb_code>/<ligand_name>', methods=['GET'])
def get_sasa_chains(pdb_code, ligand_name):
    return _dispatch_supported_lookup(
        'sasa-chains',
        params={'pdb_code': pdb_code, 'ligand_name': ligand_name},
        local_loader=lambda: jsonify(_local_get_sasa_chains_payload(pdb_code, ligand_name)),
    )

def highlight_solvent_exposed_atoms(
    molecule,
    virus_name,
    pdb_id,
    ligand_id,
    chain,
    output_svg
):
    """
    Highlight only atoms classified as solvent exposed by the validated
    `solvent_exposed_atoms` compatibility view.

    `RUPLEY_SASA_DATA` contains the full quantitative per-atom SASA population,
    including buried / zero-SASA atoms, so it must not be used as the exposed
    subset.
    """
    with _connect_local_db(LIGAND_IMAGE_REQUIRED_TABLES) as conn:
        cur = conn.cursor()
        cur.execute(
            """
            SELECT DISTINCT m.smiles_atom_index
            FROM solvent_exposed_atoms e
            JOIN SMILES_MAP_PDB m
              ON m.virus_name = e.virus_name
             AND m.pdb_id = e.pdb_id
             AND m.ligand = e.ligand
             AND m.chain = e.chain
             AND m.atom_id = e.atom_id
            WHERE e.virus_name = ?
              AND e.pdb_id = ?
              AND e.ligand = ?
              AND e.chain = ?
              AND m.smiles_atom_index IS NOT NULL
            ORDER BY m.smiles_atom_index
            """,
            (virus_name, pdb_id, ligand_id, chain),
        )
        sasa_smiles_indices = [int(row[0]) for row in cur.fetchall()]

    logging.debug(
        "Solvent-exposed 2D highlight %s/%s/%s/%s: %d mapped exposed atoms",
        virus_name,
        pdb_id,
        ligand_id,
        chain,
        len(sasa_smiles_indices),
    )

    highlight_solvent_exposed_atoms_from_indices(
        molecule,
        sasa_smiles_indices,
        output_svg,
    )

def highlight_solvent_exposed_atoms_from_indices(molecule, smiles_indices, output_svg):
    AllChem.Compute2DCoords(molecule)

    drawer = rdMolDraw2D.MolDraw2DSVG(600, 600)
    opts = drawer.drawOptions()
    opts.bondLineWidth = 2.0

    drawer.DrawMolecule(
        molecule,
        highlightAtoms=smiles_indices,
        highlightAtomColors={i: (1.0, 0.0, 0.0) for i in smiles_indices}
    )

    drawer.FinishDrawing()

    with open(output_svg, "w") as f:
        f.write(drawer.GetDrawingText())



@app.route('/compare_ligands')
def compare_ligands():
    return render_template('compare_ligands.html')



def build_query(filters):
    base_query = "SELECT * FROM Arpeggio_Contacts_Data WHERE 1=1"
    
    query_params = []
    
    # Add dynamic filters to the query
    if 'pdb_id' in filters and filters['pdb_id']:
        base_query += " AND PDB_ID = ?"
        query_params.append(filters['pdb_id'])

    if 'ligand_id' in filters and filters['ligand_id']:
        base_query += " AND Ligand = ?"
        query_params.append(filters['ligand_id'])

    if 'chain' in filters and filters['chain']:
        base_query += " AND Chain = ?"
        query_params.append(filters['chain'])
    
    if 'ligand_id' in filters and filters['ligand_id']:
        base_query += " AND Residue_Number = ?"
        query_params.append(filters['ligand_id'])
    
    return base_query, query_params

@app.route('/compare_ligand_interactions', methods=['POST'])
def compare_ligand_interactions():
    data = request.get_json(silent=True) or {}
    ligand = str(data.get('ligand') or '').strip()
    # `ligand_instance_ids` is the public occurrence-resolved contract.
    # Keep the previous spelling only for clients deployed during the migration.
    occurrence_ids = data.get('ligand_instance_ids') or data.get('occurrence_ids') or []
    legacy_pdb_ids = data.get('pdb_ids') or []
    if not ligand or not isinstance(occurrence_ids, list) or not occurrence_ids:
        return jsonify({
            'error': 'invalid_comparison_request',
            'message': 'Select a ligand and one or more mapped ligand occurrences before comparing.',
        }), 400

    mode = _normalized_backend_mode()
    if mode == "randy":
        try:
            payload = randy_post(
                'ligand-interactions/compare',
                json={
                    'ligand': ligand,
                    'ligand_instance_ids': occurrence_ids,
                },
            )
        except RandyBackendError as exc:
            return jsonify({'error': str(exc)}), exc.status_code
        return jsonify(payload)

    if mode == "auto" and randy_available():
        try:
            payload = randy_post(
                'ligand-interactions/compare',
                json={
                    'ligand': ligand,
                    'ligand_instance_ids': occurrence_ids,
                },
            )
            return jsonify(payload)
        except RandyBackendError:
            logging.warning("Falling back to local ligand interaction comparison for %s", ligand)

    try:
        payload = _local_compare_ligand_interactions_payload(ligand, occurrence_ids)
    except (TypeError, ValueError) as exc:
        return jsonify({'error': 'invalid_comparison_request', 'message': str(exc)}), 400

    if not payload['interactions_data']:
        return jsonify({
            **payload,
            'status': 'no_data',
            'message': 'No recorded non-proximal protein–ligand interactions are available for the selected ligand occurrence(s).',
        })
    return jsonify(payload)


def _local_compare_ligand_interactions_payload(ligand, occurrence_ids):
    """Return Stage-07 mapped, latest Stage-09 contacts for explicit occurrences.

    The query begins with the selected occurrence IDs.  This intentionally
    avoids materializing broad compatibility views before applying the user's
    ligand/structure selection.
    """
    try:
        normalized_ids = list(dict.fromkeys(int(value) for value in occurrence_ids))
    except (TypeError, ValueError) as exc:
        raise ValueError('Selected ligand occurrence identifiers must be integers.') from exc
    if not normalized_ids:
        raise ValueError('Select one or more mapped ligand occurrences before comparing.')

    placeholders = ', '.join('?' for _ in normalized_ids)
    query = f'''
        WITH selected_occurrences AS (
            SELECT
                i.ligand_instance_id,
                s.entry_id AS pdb_id,
                i.deposited_model_num AS model_id,
                i.auth_asym_id AS chain,
                i.auth_seq_id AS ligand_id,
                i.insertion_code_normalized AS insertion_code,
                COALESCE(sc.virus_label, 'Unknown') AS virus_name
            FROM ligand_instances AS i
            JOIN structures AS s ON s.structure_id = i.structure_id
            LEFT JOIN structure_classifications AS sc ON sc.structure_id = i.structure_id
            WHERE i.ligand_instance_id IN ({placeholders})
              AND i.label_comp_id = ?
              AND i.curation_status = 'included'
        ), latest_mapping_runs AS (
            SELECT m.ligand_instance_id, MAX(m.run_id) AS run_id
            FROM ligand_smiles_atom_mapping AS m
            JOIN selected_occurrences AS so
              ON so.ligand_instance_id = m.ligand_instance_id
            WHERE m.method_version = 'legacy_mcs_etkdg_uff_cif_v2.5'
            GROUP BY m.ligand_instance_id
        ), latest_contact_runs AS (
            SELECT ar.ligand_instance_id, MAX(ar.run_id) AS run_id
            FROM ligand_arpeggio_runs AS ar
            JOIN selected_occurrences AS so
              ON so.ligand_instance_id = ar.ligand_instance_id
            WHERE ar.status = 'completed'
            GROUP BY ar.ligand_instance_id
        ), latest_sasa_runs AS (
            SELECT sa.ligand_instance_id, MAX(sa.run_id) AS run_id
            FROM ligand_sasa_atoms AS sa
            JOIN selected_occurrences AS so
              ON so.ligand_instance_id = sa.ligand_instance_id
            WHERE sa.method_version = 'biopython-shrake_rupley-1.40-cif-v2.1'
              AND sa.status = 'complete'
            GROUP BY sa.ligand_instance_id
        ), occurrence_metrics AS (
            SELECT
                m.ligand_instance_id,
                COUNT(DISTINCT m.ligand_instance_atom_id) AS mapped_atom_count,
                COUNT(DISTINCT CASE WHEN sa.legacy_exposed = 1 THEN m.ligand_instance_atom_id END)
                    AS solvent_exposed_atom_count
            FROM ligand_smiles_atom_mapping AS m
            JOIN latest_mapping_runs AS lmr
              ON lmr.ligand_instance_id = m.ligand_instance_id
             AND lmr.run_id = m.run_id
            LEFT JOIN latest_sasa_runs AS lsr
              ON lsr.ligand_instance_id = m.ligand_instance_id
            LEFT JOIN ligand_sasa_atoms AS sa
              ON sa.ligand_instance_id = m.ligand_instance_id
             AND sa.run_id = lsr.run_id
             AND sa.ligand_instance_atom_id = m.ligand_instance_atom_id
            WHERE m.smiles_atom_index IS NOT NULL
            GROUP BY m.ligand_instance_id
        )
        SELECT
            so.ligand_instance_id,
            so.pdb_id,
            so.model_id,
            so.chain,
            so.ligand_id,
            so.insertion_code,
            so.virus_name,
            om.mapped_atom_count,
            om.solvent_exposed_atom_count,
            m.smiles_atom_index,
            r.interaction_label AS Contact,
            r.distance AS Distance,
            COALESCE(a.auth_atom_id, a.label_atom_id) AS exact_atom,
            a.atom_site_id AS atom_id,
            json_extract(r.partner_identity_json, '$.label_comp_id') AS residue,
            json_extract(r.partner_identity_json, '$.auth_seq_id') AS residue_number,
            COALESCE(
                json_extract(r.partner_identity_json, '$.auth_atom_id'),
                json_extract(r.partner_identity_json, '$.label_atom_id')
            ) AS residue_atom,
            json_extract(r.partner_identity_json, '$.auth_asym_id') AS residue_chain
        FROM selected_occurrences AS so
        JOIN latest_mapping_runs AS lmr
          ON lmr.ligand_instance_id = so.ligand_instance_id
        JOIN ligand_smiles_atom_mapping AS m
          ON m.ligand_instance_id = so.ligand_instance_id
         AND m.run_id = lmr.run_id
        JOIN latest_contact_runs AS lcr
          ON lcr.ligand_instance_id = so.ligand_instance_id
        JOIN arpeggio_raw_contact_labels AS r
          ON r.ligand_instance_id = so.ligand_instance_id
         AND r.run_id = lcr.run_id
         AND r.filter_class = 'raw_environment'
         AND r.ligand_instance_atom_id = m.ligand_instance_atom_id
        JOIN occurrence_metrics AS om
          ON om.ligand_instance_id = so.ligand_instance_id
        LEFT JOIN ligand_instance_atoms AS a
          ON a.ligand_instance_atom_id = m.ligand_instance_atom_id
        WHERE m.smiles_atom_index IS NOT NULL
        ORDER BY so.pdb_id, so.model_id, so.chain, so.ligand_id,
                 so.insertion_code, so.ligand_instance_id, m.smiles_atom_index
    '''
    with _connect_local_db(COMPARE_LIGAND_REQUIRED_TABLES) as conn:
        frame = pd.read_sql(query, conn, params=(*normalized_ids, ligand))

    if frame.empty:
        return {'interactions_data': [], 'smiles_interactions_data': []}

    frame = frame.replace({np.nan: None})
    interactions_data = []
    smiles_interactions_data = []
    for occurrence_id, occurrence_frame in frame.groupby('ligand_instance_id', sort=False):
        clean_frame = filter_valid_atom_ids(preprocess_interactions(occurrence_frame.copy()))
        if clean_frame.empty:
            continue
        first = clean_frame.iloc[0]
        insertion = first['insertion_code'] or ''
        occurrence_label = (
            f"{first['pdb_id']} · model {first['model_id']} · {first['chain']}:{first['ligand_id']}"
            f"{insertion}"
        )
        interaction_records = clean_frame.to_dict(orient='records')
        interactions_data.append({
            'pdb_id': first['pdb_id'],
            'occurrence_label': occurrence_label,
            'virus_name': first['virus_name'],
            'ligand_id': str(first['ligand_id']),
            'ligand_instance_id': int(occurrence_id),
            'model_id': first['model_id'],
            'chain': first['chain'],
            'insertion_code': first['insertion_code'],
            'mapped_atom_count': int(first['mapped_atom_count']),
            'solvent_exposed_atom_count': int(first['solvent_exposed_atom_count']),
            'interactions': interaction_records,
        })
        smiles_interactions_data.append({
            'pdb_id': first['pdb_id'],
            'occurrence_label': occurrence_label,
            'ligand_instance_id': int(occurrence_id),
            'interactions': [
                {
                    'pdb_id': row['pdb_id'],
                    'Contact': row['Contact'],
                    'smiles_atom_index': int(row['smiles_atom_index']),
                }
                for row in interaction_records
                if row['smiles_atom_index'] is not None
            ],
        })
    return {
        'interactions_data': interactions_data,
        'smiles_interactions_data': smiles_interactions_data,
    }

# Function to generate the atom bar chart that compares multiple PDB structures
def plot_atom_interactions_comparison(df_list, output_file, pdb_ids, ligand_code):
    atom_labels = []
    interaction_map = {}

    # Iterate through the DataFrame for each PDB structure
    for df, pdb_id in zip(df_list, pdb_ids):
        # Remove rows with NaN in 'smiles_atom_index' or 'Contact'
        df = df.dropna(subset=['smiles_atom_index', 'Contact'])

        df['Atom_Label'] = df['smiles_atom_index'].astype(str) + f' ({pdb_id})'
        for _, row in df.iterrows():
            atom = row['Atom_Label']
            contact = row['Contact']

            if atom not in interaction_map:
                interaction_map[atom] = {}
            if contact not in interaction_map[atom]:
                interaction_map[atom][contact] = 0

            interaction_map[atom][contact] += 1

    # Prepare the data for plotting
    traces = []
    for atom in interaction_map:
        for contact, count in interaction_map[atom].items():
            traces.append({
                'x': [atom],
                'y': [count],
                'name': contact,
                'type': 'bar'
            })

    # Plot using Plotly or Matplotlib
    plt.figure(figsize=(14, 7))
    sns.barplot(x=list(interaction_map.keys()), y=[sum(atom.values()) for atom in interaction_map.values()], palette="Set1")

    plt.title(f'Interaction Comparison Across PDBs for Ligand {ligand_code}')
    plt.xlabel('Atoms (SMILES Atom Index)')
    plt.ylabel('Interaction Count')
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()


@app.route('/get_pdb_mapping/<ligand_code>')
def get_pdb_mapping(ligand_code):
    return _dispatch_supported_lookup(
        'pdb-mapping',
        params={'ligand_code': ligand_code},
        local_loader=lambda: jsonify(_local_get_pdb_mapping_payload(ligand_code)),
    )

@app.route('/rdkittest')
def rdkit_test():
    # Serve the RDKit test page under the '/rdkittest' route
    return render_template('rdkittest.html')

@app.route('/TESTPAGE')
def TESTPAGE():
    # Serve the RDKit test page under the '/rdkittest' route
    return render_template('TESTPAGE.html')
    
@app.route('/get_ligand_options')
def get_ligand_options():
    mode = _normalized_backend_mode()
    if mode == "randy":
        try:
            payload = randy_get('ligand-options')
        except RandyBackendError as exc:
            return jsonify({'error': str(exc)}), exc.status_code
        return jsonify(payload)
    if mode == "auto" and randy_available():
        try:
            return jsonify(randy_get('ligand-options'))
        except RandyBackendError:
            logging.warning("Falling back to local ligand options")
    try:
        payload = _local_get_ligand_options_payload()
    except RandyBackendError as exc:
        return jsonify({'error': str(exc)}), 500
    return jsonify(payload)



# Function to get the SMILES string based on ligand_id
def get_smiles_from_identifier(ligand):
    return _load_smiles_for_ligand(ligand)

# Function to generate SVG from SMILES
def generate_svg_from_smiles(smiles):
    molecule = Chem.MolFromSmiles(smiles)
    if molecule:
        AllChem.Compute2DCoords(molecule)
        drawer = rdMolDraw2D.MolDraw2DSVG(400, 400)
        drawer.DrawMolecule(molecule)
        drawer.FinishDrawing()
        return drawer.GetDrawingText().replace('\n', '')
    return None


def generate_svg_with_highlight(smiles, atom_indices):
    """
    Generates SVG with specific atoms highlighted.
    :param smiles: The SMILES string to generate the molecule.
    :param atom_indices: The list of atom indices to highlight.
    :return: The SVG string with highlighted atoms.
    """
    molecule = Chem.MolFromSmiles(smiles)
    if molecule:
        AllChem.Compute2DCoords(molecule)
        drawer = rdMolDraw2D.MolDraw2DSVG(400, 400)
        try:
            drawer.DrawMolecule(molecule, highlightAtoms=atom_indices)
        except Exception as e:
            logging.error("Error highlighting atoms: %s", str(e))
            return None
        drawer.FinishDrawing()
        return drawer.GetDrawingText().replace('\n', '')
    return None


# Get SMILES and render SVG using only ligand_id
@app.route('/get_smiles_svg/<ligand_id>', methods=['GET'])
def get_smiles_svg(ligand_id):
    smiles = get_smiles_from_identifier(ligand_id)
    if smiles:
        svg = generate_svg_from_smiles(smiles)
        if svg:
            return jsonify({'svg': svg})
    return jsonify({'error': 'Unable to generate SVG for the provided SMILES'}), 404





# Get atom count for the given ligand
@app.route('/get_atom_count/<ligand_id>', methods=['GET'])
def get_atom_count(ligand_id):
    smiles = get_smiles_from_identifier(ligand_id)
    if smiles:
        molecule = Chem.MolFromSmiles(smiles)
        if molecule:
            return jsonify({'atom_count': molecule.GetNumAtoms()})
    return jsonify({'error': 'SMILES not found'}), 404



# Generate highlighted SVG for selected atoms
@app.route('/highlight_atoms', methods=['POST'])
def highlight_atoms():
    try:
        data = request.json
        print(data)  # Debug: print data to see what's being received

        atom_indices = data.get('atom_indices', [])
        ligand_id = data.get('ligand_id')  # Using only ligand_id

        # Fetch the SMILES string from the database using ligand_id
        smiles = get_smiles_from_identifier(ligand_id)
        if not smiles:
            logging.error("SMILES string not found for ligand_id: %s", ligand_id)
            return jsonify({'error': 'SMILES string not found'}), 404

        # Generate SVG with highlighted atoms
        svg = generate_svg_with_highlight(smiles, atom_indices)
        if not svg:
            return jsonify({'error': 'Failed to generate highlighted image'}), 500

        # Return the image as SVG
        return jsonify({'image': svg})

    except Exception as e:
        logging.exception("Unexpected error in /highlight_atoms")
        return jsonify({'error': 'Unexpected error occurred'}), 500





#SCRIPTS FOR PROTEIN QUERYING

data_set_queries = {
    "Solvent Exposed Atoms": "SELECT * FROM solvent_exposed_atoms WHERE pdb_id IN ({placeholders})",
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


def _protein_query_export_eligible_pdbs(
    conn,
    *,
    pdb_codes=None,
    virus_name=None,
    protein_types=None,
    ligand_filter=None,
):
    """Return PDB codes with retained, exportable ligand-centered data.

    ``Virus_Proteins`` is deliberately used only for virus/protein
    classification.  Eligibility is based on the current occurrence-level
    population: an included ``ligand_instances`` row with a selected ligand
    atom.  This is the source population for the occurrence-preserving
    ``ligand_atoms`` compatibility view used by Protein Query exports.
    """
    requested_codes = []
    for pdb_code in pdb_codes or []:
        normalized = str(pdb_code or "").strip().upper()
        if normalized and normalized not in requested_codes:
            requested_codes.append(normalized)

    normalized_proteins = [
        str(protein_type).strip()
        for protein_type in (protein_types or [])
        if str(protein_type).strip()
    ]
    normalized_virus = str(virus_name or "").strip()
    normalized_ligand = str(ligand_filter or "").strip()

    params = []
    if normalized_virus or normalized_proteins:
        if not (normalized_virus and normalized_proteins):
            return []
        placeholders = ", ".join(["?"] * len(normalized_proteins))
        query = f"""
            SELECT DISTINCT s.entry_id AS pdb_id
            FROM structures AS s
            JOIN ligand_instances AS li ON li.structure_id = s.structure_id
            JOIN ligand_instance_atoms AS lia
                ON lia.ligand_instance_id = li.ligand_instance_id
               AND lia.selected_conformer = 1
            WHERE li.curation_status = 'included'
              AND s.entry_id IN (
                  SELECT DISTINCT pdb_id
                  FROM Virus_Proteins
                  WHERE virus_name = ? AND protein IN ({placeholders})
              )
        """
        params.extend([normalized_virus, *normalized_proteins])
    else:
        query = """
            SELECT DISTINCT s.entry_id AS pdb_id
            FROM structures AS s
            JOIN ligand_instances AS li ON li.structure_id = s.structure_id
            JOIN ligand_instance_atoms AS lia
                ON lia.ligand_instance_id = li.ligand_instance_id
               AND lia.selected_conformer = 1
            WHERE li.curation_status = 'included'
        """

    if requested_codes:
        placeholders = ", ".join(["?"] * len(requested_codes))
        query += f" AND s.entry_id IN ({placeholders})"
        params.extend(requested_codes)

    if normalized_ligand:
        query += """
            AND (
                li.label_comp_id = ?
                OR li.label_comp_id IN (
                    SELECT synonym FROM Ligand_Synonyms WHERE ligand = ?
                )
                OR li.label_comp_id IN (
                    SELECT ligand FROM Ligand_Synonyms WHERE synonym = ?
                )
            )
        """
        params.extend([normalized_ligand, normalized_ligand, normalized_ligand])

    query += " ORDER BY pdb_id"
    return [row[0] for row in conn.execute(query, tuple(params)).fetchall()]


def get_exportable_protein_query_pdbs(conn, virus_name, protein_types, ligand_filter=None):
    """Central Protein Query eligibility API for selector population."""
    return _protein_query_export_eligible_pdbs(
        conn,
        virus_name=virus_name,
        protein_types=protein_types,
        ligand_filter=ligand_filter,
    )


def get_protein_query_filter_options(conn, virus_names=None, protein_types=None):
    """Return cascading Protein Query options from the exportable population.

    Every option is backed by an included ligand instance with a selected
    conformer atom, exactly as the PDB picker is.  This keeps a user from
    reaching a filter combination that cannot produce an export.
    """
    normalized_viruses = sorted({
        str(virus_name).strip()
        for virus_name in (virus_names or [])
        if str(virus_name).strip()
    })
    normalized_proteins = sorted({
        str(protein_type).strip()
        for protein_type in (protein_types or [])
        if str(protein_type).strip()
    })

    clauses = [
        "li.curation_status = 'included'",
        "lia.selected_conformer = 1",
    ]
    params = []
    if normalized_viruses:
        placeholders = ", ".join(["?"] * len(normalized_viruses))
        clauses.append(f"vp.virus_name IN ({placeholders})")
        params.extend(normalized_viruses)
    if normalized_proteins:
        placeholders = ", ".join(["?"] * len(normalized_proteins))
        clauses.append(f"vp.protein IN ({placeholders})")
        params.extend(normalized_proteins)

    if normalized_viruses:
        # The PDB eligibility helper has an efficient query plan for a
        # virus/protein slice.  Use it first, then derive the three option
        # lists from that small, known-exportable PDB set.  A direct four-way
        # join is much slower for a selected virus on the production database.
        protein_scope = {}
        placeholders = ", ".join(["?"] * len(normalized_viruses))
        for virus_name, protein_type in conn.execute(
            f"""
            SELECT DISTINCT virus_name, protein
            FROM Virus_Proteins
            WHERE virus_name IN ({placeholders})
            """,
            tuple(normalized_viruses),
        ):
            protein_scope.setdefault(virus_name, []).append(protein_type)

        eligible_pdbs = []
        for virus_name, available_proteins in protein_scope.items():
            scoped_proteins = normalized_proteins or available_proteins
            for pdb_code in get_exportable_protein_query_pdbs(
                conn, virus_name, scoped_proteins
            ):
                if pdb_code not in eligible_pdbs:
                    eligible_pdbs.append(pdb_code)

        if not eligible_pdbs:
            return {"virus_names": [], "protein_types": [], "ligands": []}

        pdb_placeholders = ", ".join(["?"] * len(eligible_pdbs))
        classification_clauses = [
            f"pdb_id IN ({pdb_placeholders})",
            f"virus_name IN ({placeholders})",
        ]
        classification_params = [*eligible_pdbs, *normalized_viruses]
        if normalized_proteins:
            protein_placeholders = ", ".join(["?"] * len(normalized_proteins))
            classification_clauses.append(f"protein IN ({protein_placeholders})")
            classification_params.extend(normalized_proteins)
        classification_rows = conn.execute(
            f"""
            SELECT DISTINCT virus_name, protein
            FROM Virus_Proteins
            WHERE {' AND '.join(classification_clauses)}
            ORDER BY virus_name, protein
            """,
            tuple(classification_params),
        ).fetchall()
        ligand_codes = {
            row[0]
            for row in conn.execute(
                f"""
                SELECT DISTINCT li.label_comp_id
                FROM structures AS s
                JOIN ligand_instances AS li ON li.structure_id = s.structure_id
                JOIN ligand_instance_atoms AS lia
                    ON lia.ligand_instance_id = li.ligand_instance_id
                WHERE s.entry_id IN ({pdb_placeholders})
                  AND li.curation_status = 'included'
                  AND lia.selected_conformer = 1
                """,
                tuple(eligible_pdbs),
            )
        }
        virus_options = sorted({row[0] for row in classification_rows})
        protein_options = sorted({row[1] for row in classification_rows})
    else:
        rows = conn.execute(
            f"""
            SELECT DISTINCT vp.virus_name, vp.protein, li.label_comp_id
            FROM Virus_Proteins AS vp
            JOIN structures AS s ON s.entry_id = vp.pdb_id
            JOIN ligand_instances AS li ON li.structure_id = s.structure_id
            JOIN ligand_instance_atoms AS lia
                ON lia.ligand_instance_id = li.ligand_instance_id
            WHERE {' AND '.join(clauses)}
            ORDER BY vp.virus_name, vp.protein, li.label_comp_id
            """,
            tuple(params),
        ).fetchall()
        ligand_codes = {row[2] for row in rows}
        virus_options = sorted({row[0] for row in rows})
        protein_options = sorted({row[1] for row in rows})
    synonym_rows = []
    if ligand_codes:
        placeholders = ", ".join(["?"] * len(ligand_codes))
        synonym_rows = conn.execute(
            f"""
            SELECT ligand, synonym
            FROM Ligand_Synonyms
            WHERE ligand IN ({placeholders}) OR synonym IN ({placeholders})
            """,
            tuple(ligand_codes) * 2,
        ).fetchall()

    canonical_for_code = {ligand_code: ligand_code for ligand_code in ligand_codes}
    synonyms_by_ligand = {}
    for ligand, synonym in synonym_rows:
        if synonym in ligand_codes:
            canonical_for_code[synonym] = ligand
        if ligand in ligand_codes:
            canonical_for_code[ligand] = ligand
        synonyms_by_ligand.setdefault(ligand, set()).add(synonym)

    available_ligands = {}
    for ligand_code in ligand_codes:
        canonical = canonical_for_code[ligand_code]
        available_ligands.setdefault(canonical, set()).update(
            synonyms_by_ligand.get(canonical, set())
        )

    return {
        "virus_names": virus_options,
        "protein_types": protein_options,
        "ligands": [
            {
                "ligand_code": ligand_code,
                "synonyms": sorted(synonym for synonym in synonyms if synonym != ligand_code),
            }
            for ligand_code, synonyms in sorted(available_ligands.items())
        ],
    }

def connect_db():
    return sqlite3.connect(str(LOCAL_DB_PATH))


PROTACABILITY_REQUIRED_TABLES = (
    "protacability_assessment",
    "protacability_lysine_proximity",
    "protacability_ligand_inventory",
)

PROTACABILITY_OPTIONAL_TABLES = (
    "protacability_warhead_linkability",
    "protacability_degrader_readiness",
    "v2_attachment_site_summary",
    "v2_attachment_site_candidates",
    "v2_attachment_site_high_priority",
)

ATTACHMENT_METHOD_VERSION = "attachment-sites-cif-v2.6"
# Website-only presentation cutoff retained from the prior attachment map.
# It groups nearby atom candidates for human review; it neither changes Stage-13
# atom scores nor persists a new scientific region model.
ATTACHMENT_DISPLAY_SITE_DISTANCE_A = 5.0
PROTACABILITY_METHOD_VERSION = "protacability-cif-v2.8"
PROTACABILITY_ATTACHMENT_TABLES = (
    "v2_attachment_site_summary",
    "v2_attachment_site_candidates",
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
    conn = sqlite3.connect(str(LOCAL_DB_PATH))
    conn.row_factory = sqlite3.Row
    return conn


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
        "SELECT 1 FROM sqlite_master WHERE type IN ('table','view') AND name = ?",
        (table_name,),
    ).fetchone()
    return row is not None


def _table_columns(conn, table_name):
    if not _table_exists(conn, table_name):
        return set()
    return {row[1] for row in conn.execute(f"PRAGMA table_info({table_name})").fetchall()}


def _attachment_tables_available(conn):
    return all(_table_exists(conn, table_name) for table_name in PROTACABILITY_ATTACHMENT_TABLES)


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
        if PROTACABILITY_ATTACHMENT_TABLES[0] in optional_tables and _attachment_tables_available(conn):
            data_sets.extend([
                "PROTACability Attachment-Site Summary",
                "PROTACability Candidate Attachment Atoms",
                "PROTACability High-Priority Attachment Sites",
            ])
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
    filters["canonical_target_ids"] = _normalize_multi_values(args.getlist("canonical_target_id"))
    filters["tiers"] = _normalize_multi_values(args.getlist("tier"))
    filters["warhead_tiers"] = _normalize_multi_values(args.getlist("warhead_tier"))
    filters["readiness_tiers"] = _normalize_multi_values(args.getlist("readiness_tier"))
    filters["evidence_levels"] = _normalize_multi_values(args.getlist("evidence_level"))
    filters["smiles_sources"] = _normalize_multi_values(args.getlist("smiles_source"))
    filters["ligand"] = (args.get("ligand") or "").strip()
    filters["ligand_instance_id"] = (args.get("ligand_instance_id") or "").strip()
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
    filters["has_attachment_sites"] = _parse_bool_filter(args.get("has_attachment_sites"))
    return filters


def _copy_protacability_filters(filters):
    return {
        "virus_names": list(filters.get("virus_names") or []),
        "protein_types": list(filters.get("protein_types") or []),
        "canonical_target_ids": list(filters.get("canonical_target_ids") or []),
        "tiers": list(filters.get("tiers") or []),
        "warhead_tiers": list(filters.get("warhead_tiers") or []),
        "readiness_tiers": list(filters.get("readiness_tiers") or []),
        "evidence_levels": list(filters.get("evidence_levels") or []),
        "smiles_sources": list(filters.get("smiles_sources") or []),
        "ligand": filters.get("ligand") or "",
        "ligand_instance_id": filters.get("ligand_instance_id") or "",
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
        "has_attachment_sites": filters.get("has_attachment_sites"),
    }


def _clear_protacability_filter_dimension(filters, key):
    if key == "virus_names":
        filters["virus_names"] = []
    elif key == "protein_types":
        filters["protein_types"] = []
    elif key == "canonical_target_ids":
        filters["canonical_target_ids"] = []
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
        "has_attachment_sites",
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
    readiness_rows, warhead_rows, attachment_rows = _load_protacability_enrichment_tables(conn)
    rows = _decorate_protacability_rows(
        _load_protacability_assessment_rows(conn),
        collapse_labels=collapse_labels,
        readiness_rows=readiness_rows,
        warhead_rows=warhead_rows,
        attachment_rows=attachment_rows,
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


def _load_optional_table_rows(conn, table_name, method_version=None):
    if not _table_exists(conn, table_name):
        return []
    query = f"SELECT * FROM {table_name}"
    params = []
    if method_version is not None and "method_version" in _table_columns(conn, table_name):
        query += " WHERE method_version = ?"
        params.append(method_version)
    return [dict(row) for row in conn.execute(query, tuple(params)).fetchall()]


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
    for ligand in _split_candidate_ligands(row.get("candidate_ligand_resnames")):
        if ligand not in ligands:
            ligands.append(ligand)
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


def _normalize_attachment_key(row):
    if not row:
        return None

    ligand_instance_id = row.get("ligand_instance_id")
    try:
        if ligand_instance_id not in (None, ""):
            return ("ligand_instance_id", int(ligand_instance_id))
    except (TypeError, ValueError):
        pass

    pdb_code = _normalize_text_key(row.get("pdb_code"))
    ligand_resname = _normalize_text_key(row.get("ligand_resname") or row.get("best_ligand_resname"))
    ligand_chain = _normalize_text_key(row.get("ligand_chain") or row.get("best_ligand_chain"))
    ligand_residue_id = row.get("ligand_residue_id")
    if ligand_residue_id in (None, ""):
        ligand_residue_id = row.get("best_ligand_residue_id")
    ligand_insertion_code = _normalize_text_key(row.get("ligand_insertion_code"))
    model_id = row.get("model_id")
    try:
        model_id = int(model_id or 1)
    except (TypeError, ValueError):
        model_id = 1
    try:
        ligand_residue_id = int(ligand_residue_id)
    except (TypeError, ValueError):
        return None
    if not (pdb_code and ligand_resname and ligand_chain):
        return None
    return (
        "legacy_occurrence",
        pdb_code,
        model_id,
        ligand_chain,
        ligand_residue_id,
        ligand_insertion_code,
        ligand_resname,
    )

def _build_attachment_index(attachment_rows):
    index = {}
    for row in attachment_rows or []:
        key = _normalize_attachment_key(row)
        if key:
            index[key] = dict(row)
    return index


def _attachment_defaults():
    return {
        "attachment_analysis_id": None,
        "attachment_method_version": None,
        "attachment_analysis_status": None,
        "attachment_eligibility_status": None,
        "attachment_mapping_status": None,
        "attachment_region_count": 0,
        "attachment_candidate_atom_count": 0,
        "best_attachment_score": None,
        "best_attachment_confidence": None,
        "has_attachment_site_evidence": 0,
        "has_candidate_attachment_regions": 0,
        "attachment_instance_resolution_status": None,
        "attachment_instance_ambiguity_flag": 0,
        "attachment_site_model": "atom-level-v2.6",
        "attachment_region_semantics_available": False,
        "attachment_display_site_count": 0,
        "mapped_atom_count": 0,
        "attachment_mapped_atom_count": 0,
        "attachment_exposed_mapped_atom_count": 0,
        "attachment_outward_supported_candidate_count": 0,
        "attachment_clear_exit_supported_candidate_count": 0,
        "attachment_chemically_supported_candidate_count": 0,
        "attachment_direct_candidate_count": 0,
        "attachment_conditional_candidate_count": 0,
        "attachment_conditional_clear_exit_candidate_count": 0,
        "attachment_high_priority_atom_count": 0,
        "attachment_high_priority_direct_atom_count": 0,
        "best_attachment_priority_tier": None,
        "best_attachment_atom": None,
    }


def _short_attachment_tier(value):
    text_value = str(value or "").strip()
    if not text_value:
        return None
    for prefix in ("High", "Moderate", "Exploratory", "Low"):
        if text_value.lower().startswith(prefix.lower()):
            return prefix
    return text_value


def _attachment_summary_from_match(attachment_match):
    summary = _attachment_defaults()
    if not attachment_match:
        return summary

    candidate_count = int(_numeric_value(
        attachment_match.get("candidate_attachment_atom_count")
        or attachment_match.get("attachment_candidate_atom_count")
    ))
    chemically_supported = int(_numeric_value(
        attachment_match.get("chemically_supported_candidate_count")
    ))
    high_count = int(_numeric_value(
        attachment_match.get("high_priority_attachment_atom_count")
    ))
    best_score = attachment_match.get("top_attachment_site_score")
    if best_score in (None, ""):
        best_score = attachment_match.get("best_attachment_score")

    best_tier = (
        attachment_match.get("top_attachment_priority_tier")
        or attachment_match.get("best_attachment_priority_tier")
    )
    short_tier = _short_attachment_tier(best_tier)
    has_evidence = int(candidate_count > 0)

    summary.update({
        "attachment_analysis_id": (
            attachment_match.get("attachment_summary_id")
            or attachment_match.get("analysis_id")
        ),
        "attachment_method_version": (
            attachment_match.get("method_version")
            or attachment_match.get("attachment_method_version")
        ),
        "attachment_analysis_status": (
            attachment_match.get("status")
            or attachment_match.get("analysis_status")
        ),
        "attachment_eligibility_status": (
            "candidate_atoms_present" if has_evidence else "no_candidate_atoms"
        ),
        "attachment_mapping_status": (
            "mapped_atoms_present"
            if _has_positive_value(attachment_match.get("mapped_atom_count"))
            else "no_mapped_atoms"
        ),
        "attachment_region_count": 0,
        "attachment_candidate_atom_count": candidate_count,
        "best_attachment_score": best_score,
        "best_attachment_confidence": short_tier,
        "has_attachment_site_evidence": has_evidence,
        "has_candidate_attachment_regions": 0,
        "attachment_instance_resolution_status": "ligand_instance_id",
        "attachment_instance_ambiguity_flag": 0,
        "attachment_site_model": "atom-level-v2.6",
        "attachment_region_semantics_available": False,
        "attachment_display_site_count": int(_numeric_value(
            attachment_match.get("attachment_display_site_count")
        )),
        "attachment_mapped_atom_count": int(_numeric_value(
            attachment_match.get("mapped_atom_count")
        )),
        # Kept as a concise API-facing alias for consumers that render the
        # attachment summary independently of the parent assessment row.
        "mapped_atom_count": int(_numeric_value(
            attachment_match.get("mapped_atom_count")
        )),
        "attachment_exposed_mapped_atom_count": int(_numeric_value(
            attachment_match.get("exposed_mapped_atom_count")
        )),
        "attachment_outward_supported_candidate_count": int(_numeric_value(
            attachment_match.get("outward_supported_candidate_count")
        )),
        "attachment_clear_exit_supported_candidate_count": int(_numeric_value(
            attachment_match.get("clear_exit_supported_candidate_count")
        )),
        "attachment_chemically_supported_candidate_count": chemically_supported,
        "attachment_direct_candidate_count": int(_numeric_value(
            attachment_match.get("direct_attachment_candidate_count")
        )),
        "attachment_conditional_candidate_count": int(_numeric_value(
            attachment_match.get("conditional_substitution_candidate_count")
        )),
        "attachment_conditional_clear_exit_candidate_count": int(_numeric_value(
            attachment_match.get("conditional_clear_exit_candidate_count")
        )),
        "attachment_high_priority_atom_count": high_count,
        "attachment_high_priority_direct_atom_count": int(_numeric_value(
            attachment_match.get("high_priority_direct_attachment_atom_count")
        )),
        "best_attachment_priority_tier": best_tier,
        "best_attachment_atom": (
            attachment_match.get("top_attachment_exact_atom")
            or attachment_match.get("best_attachment_atom")
        ),
    })
    return summary

def _json_list(value):
    if value in (None, ""):
        return []
    if isinstance(value, list):
        return value
    try:
        parsed = json.loads(value)
    except (TypeError, ValueError):
        return []
    return parsed if isinstance(parsed, list) else []


def _json_dict(value):
    if value in (None, ""):
        return {}
    if isinstance(value, dict):
        return value
    try:
        parsed = json.loads(value)
    except (TypeError, ValueError):
        return {}
    return parsed if isinstance(parsed, dict) else {}


def _load_attachment_analysis_rows(conn):
    # Do not query the v2 compatibility views here.  Each view expands a
    # correlated latest-run subquery and, when nested in the former per-summary
    # candidate lookup below, made the initial PROTACability request effectively
    # quadratic.  The release tables are immutable; this is just a set-based
    # read of their same latest complete records.
    if not (
        _table_exists(conn, "protacability_attachment_site_summary")
        and _table_exists(conn, "protacability_attachment_sites")
    ):
        return []

    query = """
        WITH current_summaries AS (
            SELECT
                s.*,
                ROW_NUMBER() OVER (
                    PARTITION BY s.ligand_instance_id
                    ORDER BY s.run_id DESC
                ) AS summary_rank
            FROM protacability_attachment_site_summary s
            WHERE s.method_version = ?
              AND s.status = 'complete'
        ),
        current_candidate_tiers AS (
            SELECT
                a.ligand_instance_id,
                a.attachment_priority_tier,
                ROW_NUMBER() OVER (
                    PARTITION BY a.ligand_instance_id
                    ORDER BY
                        a.run_id DESC,
                        a.attachment_priority_score DESC,
                        a.attachment_site_id
                ) AS candidate_rank
            FROM protacability_attachment_sites a
            WHERE a.method_version = ?
              AND a.candidate_attachment_atom = 1
        )
        SELECT
            s.*,
            c.attachment_priority_tier AS top_attachment_priority_tier
        FROM current_summaries s
        LEFT JOIN current_candidate_tiers c
          ON c.ligand_instance_id = s.ligand_instance_id
         AND c.candidate_rank = 1
        WHERE s.summary_rank = 1
    """
    attachment_rows = [
        dict(row)
        for row in conn.execute(
            query,
            (ATTACHMENT_METHOD_VERSION, ATTACHMENT_METHOD_VERSION),
        ).fetchall()
    ]
    # Compute the retained website display grouping in one batch.  This keeps
    # summary-table badges responsive without changing the atom-level release
    # tables or treating the clusters as scientific Stage-13 regions.
    candidate_rows = [
        dict(row)
        for row in conn.execute(
            """
            WITH ranked_candidates AS (
                SELECT
                    s.ligand_instance_id,
                    s.ligand_instance_atom_id,
                    s.run_id,
                    s.atom_site_id,
                    s.exact_atom,
                    s.smiles_atom_indices,
                    s.attachment_priority_score,
                    s.candidate_attachment_atom,
                    MAX(s.run_id) OVER (
                        PARTITION BY s.ligand_instance_id
                    ) AS current_run_id
                FROM protacability_attachment_sites s
                WHERE s.method_version = ?
                  AND s.candidate_attachment_atom = 1
            )
            SELECT
                s.ligand_instance_id,
                s.atom_site_id,
                s.exact_atom,
                s.smiles_atom_indices,
                s.attachment_priority_score,
                s.candidate_attachment_atom,
                a.x, a.y, a.z,
                l.canonical_smiles
            FROM ranked_candidates s
            LEFT JOIN ligand_instance_atoms a
              ON a.ligand_instance_atom_id = s.ligand_instance_atom_id
            LEFT JOIN ligand_instances li
              ON li.ligand_instance_id = s.ligand_instance_id
            LEFT JOIN ligands l
              ON l.ligand_id = li.ligand_id
            WHERE s.run_id = s.current_run_id
            """,
            (ATTACHMENT_METHOD_VERSION,),
        ).fetchall()
    ]
    candidates_by_instance = defaultdict(list)
    for candidate in candidate_rows:
        candidate["pdb_atom_serial"] = candidate.get("atom_site_id")
        candidate["pdb_atom_name"] = candidate.get("exact_atom")
        candidate["attachment_score"] = candidate.get("attachment_priority_score")
        candidate["candidate_attachment_flag"] = candidate.get("candidate_attachment_atom")
        candidate["smiles_atom_indices_list"] = _parse_smiles_atom_indices(
            candidate.get("smiles_atom_indices")
        )
        candidates_by_instance[candidate["ligand_instance_id"]].append(candidate)
    for attachment in attachment_rows:
        attachment["attachment_display_site_count"] = len(
            _attachment_display_site_clusters(
                candidates_by_instance.get(attachment["ligand_instance_id"], [])
            )
        )
    return attachment_rows


def _local_protacability_raw_export(conn, raw_export):
    """Return a human-labelled current-release PROTACability raw export.

    Attachment exports deliberately avoid the correlated compatibility views.
    They expose the same latest Stage-13 records used by the application,
    without carrying the retired region vocabulary forward.
    """
    table_exports = {
        "PROTACability Assessment": "protacability_assessment",
        "PROTACability Lysine Proximity": "protacability_lysine_proximity",
        "PROTACability Ligand Inventory": "protacability_ligand_inventory",
        "PROTACability Warhead Linkability": "protacability_warhead_linkability",
        "PROTACability Degrader Readiness": "protacability_degrader_readiness",
    }
    if raw_export in table_exports:
        table_name = table_exports[raw_export]
        if not _table_exists(conn, table_name):
            raise KeyError(raw_export)
        return table_name, pd.read_sql_query(f"SELECT * FROM {table_name}", conn)

    attachment_queries = {
        "PROTACability Attachment-Site Summary": """
            WITH ranked AS (
                SELECT s.*, ROW_NUMBER() OVER (
                    PARTITION BY s.ligand_instance_id ORDER BY s.run_id DESC
                ) AS current_rank
                FROM protacability_attachment_site_summary s
                WHERE s.method_version = ? AND s.status = 'complete'
            )
            SELECT * FROM ranked WHERE current_rank = 1
        """,
        "PROTACability Candidate Attachment Atoms": """
            WITH ranked AS (
                SELECT a.*, MAX(a.run_id) OVER (
                    PARTITION BY a.ligand_instance_id
                ) AS current_run_id
                FROM protacability_attachment_sites a
                WHERE a.method_version = ? AND a.candidate_attachment_atom = 1
            )
            SELECT * FROM ranked WHERE run_id = current_run_id
        """,
        "PROTACability High-Priority Attachment Sites": """
            WITH ranked AS (
                SELECT a.*, MAX(a.run_id) OVER (
                    PARTITION BY a.ligand_instance_id
                ) AS current_run_id
                FROM protacability_attachment_sites a
                WHERE a.method_version = ? AND a.high_priority_attachment_atom = 1
            )
            SELECT * FROM ranked WHERE run_id = current_run_id
        """,
    }
    query = attachment_queries.get(raw_export)
    if query is None:
        raise KeyError(raw_export)
    if not (
        _table_exists(conn, "protacability_attachment_site_summary")
        and _table_exists(conn, "protacability_attachment_sites")
    ):
        raise KeyError(raw_export)
    return raw_export.replace("PROTACability ", "").lower().replace(" ", "_"), pd.read_sql_query(
        query,
        conn,
        params=(ATTACHMENT_METHOD_VERSION,),
    )


def _empty_attachment_graph_payload():
    return {"graph_id": None, "nodes": [], "bonds": []}


def _attachment_graph_payload(conn, analysis_dict, atoms):
    graph_id = analysis_dict.get("graph_id")
    payload = {"graph_id": graph_id, "nodes": [], "bonds": []}
    try:
        graph_id_int = int(graph_id)
    except (TypeError, ValueError):
        return _empty_attachment_graph_payload()

    atom_by_index = {}
    for atom in atoms:
        try:
            smiles_index = int(atom.get("smiles_atom_index"))
        except (TypeError, ValueError):
            continue
        current = atom_by_index.get(smiles_index)
        if (
            current is None
            or int(bool(atom.get("candidate_attachment_flag"))) > int(bool(current.get("candidate_attachment_flag")))
            or (atom.get("attachment_score") or 0) > (current.get("attachment_score") or 0)
        ):
            atom_by_index[smiles_index] = atom

    if _table_exists(conn, "Ligand_SMILES_Atoms"):
        graph_atoms = conn.execute(
            """
            SELECT smiles_atom_index, element, atomic_number, is_aromatic, is_in_ring
            FROM Ligand_SMILES_Atoms
            WHERE graph_id=?
            ORDER BY smiles_atom_index
            """,
            (graph_id_int,),
        ).fetchall()
        for graph_atom in graph_atoms:
            node = dict(graph_atom)
            detail = atom_by_index.get(node.get("smiles_atom_index"))
            if detail:
                node.update(
                    {
                        "pdb_atom_serial": detail.get("pdb_atom_serial"),
                        "pdb_atom_name": detail.get("pdb_atom_name"),
                        "region_id": detail.get("region_id"),
                        "candidate_attachment_flag": detail.get("candidate_attachment_flag"),
                        "surface_defining_flag": detail.get("surface_defining_flag"),
                        "attachment_score": detail.get("attachment_score"),
                        "confidence": detail.get("confidence"),
                    }
                )
            payload["nodes"].append(node)
    else:
        for smiles_index, detail in sorted(atom_by_index.items()):
            payload["nodes"].append(
                {
                    "smiles_atom_index": smiles_index,
                    "element": detail.get("element"),
                    "pdb_atom_serial": detail.get("pdb_atom_serial"),
                    "pdb_atom_name": detail.get("pdb_atom_name"),
                    "region_id": detail.get("region_id"),
                    "candidate_attachment_flag": detail.get("candidate_attachment_flag"),
                    "surface_defining_flag": detail.get("surface_defining_flag"),
                    "attachment_score": detail.get("attachment_score"),
                    "confidence": detail.get("confidence"),
                }
            )

    if _table_exists(conn, "Ligand_SMILES_Bonds"):
        payload["bonds"] = [
            dict(row)
            for row in conn.execute(
                """
                SELECT smiles_bond_index, begin_atom_index, end_atom_index, bond_type, bond_order, is_aromatic, is_in_ring
                FROM Ligand_SMILES_Bonds
                WHERE graph_id=?
                ORDER BY smiles_bond_index
                """,
                (graph_id_int,),
            ).fetchall()
        ]
    return payload


def _parse_smiles_atom_indices(value):
    if value in (None, ""):
        return []
    indices = []
    for token in re.split(r"[;,|\s]+", str(value).strip()):
        if not token:
            continue
        try:
            idx = int(float(token))
        except (TypeError, ValueError):
            continue
        if idx not in indices:
            indices.append(idx)
    return indices


def _attachment_display_site_clusters(atoms):
    """Group nearby v2.6 candidate atoms for website display only.

    Stage 13 remains atom-specific.  This recreates the previous map's useful
    human-review rule (within 5 Å and within two ligand bonds) without writing
    ASR rows or substituting a cluster score for any atom score.
    """
    candidates = [
        atom for atom in (atoms or [])
        if int(_numeric_value(atom.get("candidate_attachment_flag") or atom.get("candidate_attachment_atom")))
    ]
    for candidate in candidates:
        candidate["display_site_id"] = None
    if not candidates:
        return []

    parent = list(range(len(candidates)))

    def find(index):
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    def union(left, right):
        left_root, right_root = find(left), find(right)
        if left_root != right_root:
            parent[right_root] = left_root

    def coordinates(atom):
        values = []
        for axis in ("x", "y", "z"):
            value = atom.get(axis, atom.get(f"atom_{axis}"))
            try:
                values.append(float(value))
            except (TypeError, ValueError):
                return None
        return values

    # The historic display cutoff combined coordinate proximity with a local
    # ligand-bond neighborhood.  Coordinate proximity alone would join folded
    # but chemically distant portions of DR7, so retain the bond condition.
    canonical_smiles = next(
        (atom.get("canonical_smiles") for atom in candidates if atom.get("canonical_smiles")),
        None,
    )
    molecule = Chem.MolFromSmiles(canonical_smiles) if canonical_smiles else None

    candidate_smiles_indices = [
        _parse_smiles_atom_indices(candidate.get("smiles_atom_indices"))
        for candidate in candidates
    ]
    bond_neighborhoods = {}

    def nearby_bond_indices(atom_index):
        """Return indices up to two bonds away without repeated RDKit paths."""
        if molecule is None or atom_index < 0 or atom_index >= molecule.GetNumAtoms():
            return set()
        if atom_index not in bond_neighborhoods:
            direct = [neighbor.GetIdx() for neighbor in molecule.GetAtomWithIdx(atom_index).GetNeighbors()]
            bond_neighborhoods[atom_index] = {atom_index, *direct}
            for neighbor_index in direct:
                bond_neighborhoods[atom_index].update(
                    neighbor.GetIdx()
                    for neighbor in molecule.GetAtomWithIdx(neighbor_index).GetNeighbors()
                )
        return bond_neighborhoods[atom_index]

    def close_in_bond_graph(left_index, right_index):
        if molecule is None:
            return False
        right_indices = set(candidate_smiles_indices[right_index])
        return any(
            right_indices.intersection(nearby_bond_indices(atom_index))
            for atom_index in candidate_smiles_indices[left_index]
        )

    candidate_coordinates = [coordinates(atom) for atom in candidates]
    for left, left_coordinates in enumerate(candidate_coordinates):
        if left_coordinates is None:
            continue
        for right in range(left + 1, len(candidates)):
            right_coordinates = candidate_coordinates[right]
            if right_coordinates is None:
                continue
            distance_squared = sum(
                (left_coordinates[axis] - right_coordinates[axis]) ** 2
                for axis in range(3)
            )
            if (
                distance_squared <= ATTACHMENT_DISPLAY_SITE_DISTANCE_A ** 2
                and close_in_bond_graph(left, right)
            ):
                union(left, right)

    grouped = defaultdict(list)
    for index, atom in enumerate(candidates):
        grouped[find(index)].append(atom)

    # A display *site* has to be a bonded local group.  Retain isolated v2.6
    # candidate atoms in the raw atom table/map, but do not relabel them as a
    # legacy-style site cluster.
    clusters = [members for members in grouped.values() if len(members) > 1]
    clusters.sort(key=lambda members: (
        -max(_numeric_value(member.get("attachment_score") or member.get("attachment_priority_score")) for member in members),
        min(_numeric_value(member.get("pdb_atom_serial") or member.get("atom_site_id")) for member in members),
        sorted(str(member.get("pdb_atom_name") or member.get("exact_atom") or "") for member in members),
    ))

    display_clusters = []
    for number, members in enumerate(clusters, start=1):
        ordered_members = sorted(
            members,
            key=lambda member: (
                -_numeric_value(member.get("attachment_score") or member.get("attachment_priority_score")),
                str(member.get("pdb_atom_name") or member.get("exact_atom") or ""),
            ),
        )
        site_id = f"Site {number}"
        for member in ordered_members:
            member["display_site_id"] = site_id
        best = ordered_members[0]
        display_clusters.append({
            "display_site_id": site_id,
            "candidate_atom_serials": [
                int(member["pdb_atom_serial"])
                for member in ordered_members
                if member.get("pdb_atom_serial") not in (None, "")
            ],
            "candidate_atom_names": [
                member.get("pdb_atom_name") or member.get("exact_atom")
                for member in ordered_members
                if member.get("pdb_atom_name") or member.get("exact_atom")
            ],
            "surface_atom_serials": [
                int(member["pdb_atom_serial"])
                for member in ordered_members
                if member.get("pdb_atom_serial") not in (None, "")
                and int(_numeric_value(member.get("surface_defining_flag") or member.get("solvent_exposed")))
            ],
            "chemically_supported_candidate_count": sum(
                int(bool(member.get("chemically_supported"))) for member in ordered_members
            ),
            "high_priority_atom_count": sum(
                int(bool(member.get("high_priority_attachment_atom"))) for member in ordered_members
            ),
            "best_attachment_score": best.get("attachment_score") or best.get("attachment_priority_score"),
            "best_attachment_priority_tier": best.get("priority_tier_short") or _short_attachment_tier(best.get("attachment_priority_tier")),
        })
    return display_clusters


def _attachment_graph_from_canonical_smiles(ligand_instance_id, atoms):
    """Build the complete ligand graph used by the attachment-site display.

    The release tables map only attachment candidates.  The map must still
    include every atom and bond in the canonical ligand structure so reviewers
    can see candidate clusters in their actual chemical context.
    """
    canonical_smiles = next(
        (atom.get("canonical_smiles") for atom in (atoms or []) if atom.get("canonical_smiles")),
        None,
    )
    molecule = Chem.MolFromSmiles(canonical_smiles) if canonical_smiles else None
    payload = {
        "graph_id": f"ligand_instance:{ligand_instance_id}",
        "nodes": [],
        "bonds": [],
        "site_model": "atom-level-v2.6",
        "region_semantics_available": False,
    }
    if molecule is None:
        return payload

    details_by_index = {}
    for atom in atoms or []:
        for smiles_index in atom.get("smiles_atom_indices_list") or _parse_smiles_atom_indices(atom.get("smiles_atom_indices")):
            existing = details_by_index.get(smiles_index)
            if existing is None or _numeric_value(atom.get("attachment_score")) > _numeric_value(existing.get("attachment_score")):
                details_by_index[smiles_index] = atom

    for rdkit_atom in molecule.GetAtoms():
        smiles_index = rdkit_atom.GetIdx()
        detail = details_by_index.get(smiles_index, {})
        payload["nodes"].append({
            "smiles_atom_index": smiles_index,
            "element": rdkit_atom.GetSymbol(),
            "atomic_number": rdkit_atom.GetAtomicNum(),
            "is_aromatic": int(rdkit_atom.GetIsAromatic()),
            "is_in_ring": int(rdkit_atom.IsInRing()),
            "pdb_atom_serial": detail.get("pdb_atom_serial"),
            "pdb_atom_name": detail.get("pdb_atom_name"),
            "display_site_id": detail.get("display_site_id"),
            "candidate_attachment_flag": int(bool(detail.get("candidate_attachment_flag"))),
            "surface_defining_flag": int(bool(detail.get("surface_defining_flag"))),
            "attachment_score": detail.get("attachment_score"),
            "confidence": detail.get("confidence"),
            "priority_tier_short": detail.get("priority_tier_short"),
            "atom_chemical_role": detail.get("atom_chemical_role"),
            "direct_attachment_support": detail.get("direct_attachment_support"),
            "conditional_substitution_support": detail.get("conditional_substitution_support"),
            "chemically_supported": detail.get("chemically_supported"),
        })

    for bond_index, bond in enumerate(molecule.GetBonds()):
        payload["bonds"].append({
            "smiles_bond_index": bond_index,
            "begin_atom_index": bond.GetBeginAtomIdx(),
            "end_atom_index": bond.GetEndAtomIdx(),
            "bond_type": str(bond.GetBondType()),
            "bond_order": bond.GetBondTypeAsDouble(),
            "is_aromatic": int(bond.GetIsAromatic()),
            "is_in_ring": int(bond.IsInRing()),
        })
    return payload


def _resolve_attachment_ligand_instance_id(conn, row):
    if not row:
        return None

    try:
        ligand_instance_id = row.get("ligand_instance_id")
        if ligand_instance_id not in (None, ""):
            return int(ligand_instance_id)
    except (TypeError, ValueError):
        pass

    pdb_code = _normalize_text_key(row.get("pdb_code"))
    ligand_resname = _normalize_text_key(
        row.get("ligand_resname") or row.get("best_ligand_resname")
    )
    ligand_chain = _normalize_text_key(
        row.get("ligand_chain") or row.get("best_ligand_chain")
    )
    ligand_residue_id = row.get("ligand_residue_id")
    if ligand_residue_id in (None, ""):
        ligand_residue_id = row.get("best_ligand_residue_id")
    ligand_insertion_code = _normalize_text_key(row.get("ligand_insertion_code"))
    model_id = row.get("model_id")

    if not (pdb_code and ligand_resname and ligand_chain and ligand_residue_id not in (None, "")):
        return None

    params = [pdb_code, ligand_resname, ligand_chain, str(ligand_residue_id)]
    query = """
        SELECT ligand_instance_id
        FROM protacability_ligand_inventory
        WHERE UPPER(TRIM(pdb_code)) = ?
          AND UPPER(TRIM(ligand_resname)) = ?
          AND UPPER(TRIM(ligand_chain)) = ?
          AND CAST(ligand_residue_id AS TEXT) = ?
    """

    if ligand_insertion_code:
        query += " AND UPPER(COALESCE(TRIM(ligand_insertion_code), '')) = ?"
        params.append(ligand_insertion_code)

    if model_id not in (None, "", 0, "0"):
        query += " AND CAST(model_id AS TEXT) = ?"
        params.append(str(model_id))

    query += " ORDER BY ligand_instance_id LIMIT 2"
    matches = conn.execute(query, tuple(params)).fetchall()
    if len(matches) == 1:
        return int(matches[0][0])

    if not matches and (model_id not in (None, "", 0, "0") or ligand_insertion_code):
        fallback = conn.execute(
            """
            SELECT ligand_instance_id
            FROM protacability_ligand_inventory
            WHERE UPPER(TRIM(pdb_code)) = ?
              AND UPPER(TRIM(ligand_resname)) = ?
              AND UPPER(TRIM(ligand_chain)) = ?
              AND CAST(ligand_residue_id AS TEXT) = ?
            ORDER BY ligand_instance_id
            LIMIT 2
            """,
            (pdb_code, ligand_resname, ligand_chain, str(ligand_residue_id)),
        ).fetchall()
        if len(fallback) == 1:
            return int(fallback[0][0])

    return None


def _attachment_detail_payload(conn, row):
    if not _attachment_tables_available(conn):
        return {
            "data_available": False,
            "summary": _attachment_defaults(),
            "regions": [],
            "display_site_clusters": [],
            "atoms": [],
            "candidate_atom_serials": [],
            "chemically_supported_atom_serials": [],
            "priority_atom_serials": [],
            "high_priority_atom_serials": [],
            "surface_atom_serials": [],
            "graph": _empty_attachment_graph_payload(),
            "site_model": "atom-level-v2.6",
            "region_semantics_available": False,
        }

    ligand_instance_id = _resolve_attachment_ligand_instance_id(conn, row)
    if ligand_instance_id is None:
        return {
            "data_available": True,
            "summary": _attachment_defaults(),
            "regions": [],
            "display_site_clusters": [],
            "atoms": [],
            "candidate_atom_serials": [],
            "chemically_supported_atom_serials": [],
            "priority_atom_serials": [],
            "high_priority_atom_serials": [],
            "surface_atom_serials": [],
            "graph": _empty_attachment_graph_payload(),
            "site_model": "atom-level-v2.6",
            "region_semantics_available": False,
            "message": "No unique ligand occurrence could be resolved for attachment-site lookup.",
        }

    summary_row = conn.execute(
        """
        SELECT *
        FROM protacability_attachment_site_summary
        WHERE ligand_instance_id = ?
          AND method_version = ?
          AND status = 'complete'
        ORDER BY run_id DESC
        LIMIT 1
        """,
        (ligand_instance_id, ATTACHMENT_METHOD_VERSION),
    ).fetchone()

    if not summary_row:
        return {
            "data_available": True,
            "summary": _attachment_defaults(),
            "regions": [],
            "display_site_clusters": [],
            "atoms": [],
            "candidate_atom_serials": [],
            "chemically_supported_atom_serials": [],
            "priority_atom_serials": [],
            "high_priority_atom_serials": [],
            "surface_atom_serials": [],
            "graph": _empty_attachment_graph_payload(),
            "site_model": "atom-level-v2.6",
            "region_semantics_available": False,
            "ligand_instance_id": ligand_instance_id,
        }

    summary_dict = dict(summary_row)

    site_rows = [
        dict(site)
        for site in conn.execute(
            """
            SELECT s.*, a.x, a.y, a.z, l.canonical_smiles
            FROM protacability_attachment_sites s
            LEFT JOIN ligand_instance_atoms a
              ON a.ligand_instance_atom_id = s.ligand_instance_atom_id
            LEFT JOIN ligand_instances li
              ON li.ligand_instance_id = s.ligand_instance_id
            LEFT JOIN ligands l
              ON l.ligand_id = li.ligand_id
            WHERE s.ligand_instance_id = ?
              AND s.method_version = ?
              AND s.run_id = ?
              AND s.candidate_attachment_atom = 1
            ORDER BY
                s.attachment_priority_score DESC,
                s.chemical_support DESC,
                s.exact_atom,
                s.attachment_site_id
            """,
            (
                ligand_instance_id,
                ATTACHMENT_METHOD_VERSION,
                summary_dict["run_id"],
            ),
        ).fetchall()
    ]

    # The best tier is a presentation label for the current atom-specific
    # candidates.  It is derived from the same run as the authoritative summary
    # instead of re-entering the expensive compatibility view.
    summary_dict["top_attachment_priority_tier"] = (
        site_rows[0].get("attachment_priority_tier") if site_rows else None
    )
    summary = _attachment_summary_from_match(summary_dict)

    atoms = []
    candidate_atom_serials = []
    chemically_supported_atom_serials = []
    priority_atom_serials = []
    high_priority_atom_serials = []
    surface_atom_serials = []
    graph_nodes_by_index = {}

    for site in site_rows:
        serial = site.get("atom_site_id")
        try:
            serial_int = int(serial) if serial not in (None, "") else None
        except (TypeError, ValueError):
            serial_int = None

        smiles_indices = _parse_smiles_atom_indices(site.get("smiles_atom_indices"))
        tier_short = _short_attachment_tier(site.get("attachment_priority_tier"))
        chemically_supported = bool(
            site.get("direct_attachment_support")
            or site.get("conditional_substitution_support")
            or site.get("chemical_support")
        )

        atom = {
            **site,
            "pdb_atom_serial": serial_int,
            "pdb_atom_name": site.get("exact_atom"),
            "smiles_atom_index": smiles_indices[0] if smiles_indices else None,
            "smiles_atom_indices_list": smiles_indices,
            "candidate_attachment_flag": int(bool(site.get("candidate_attachment_atom"))),
            "surface_defining_flag": int(bool(site.get("solvent_exposed"))),
            "attachment_score": site.get("attachment_priority_score"),
            "confidence": tier_short,
            "region_id": None,
            "priority_tier_short": tier_short,
            "chemically_supported": int(chemically_supported),
            "interaction_types": [
                part.strip()
                for part in str(site.get("contact_labels") or "").split(";")
                if part.strip()
            ],
            "functional_group_annotations": [
                part.strip()
                for part in str(site.get("functional_groups") or "").split(";")
                if part.strip()
            ],
            "reasons": [
                value for value in [
                    site.get("chemical_rationale"),
                    "Mapped" if site.get("mapped") else None,
                    "Solvent exposed" if site.get("solvent_exposed") else None,
                    "Points away from pocket" if site.get("points_away_from_pocket") else None,
                    "Local corridor clear" if site.get("local_corridor_clear") else None,
                ] if value
            ],
            "cautions": [
                value for value in [
                    "Strong protein contact present"
                    if _numeric_value(site.get("strong_contact_count")) > 0 else None,
                    "Chemical context only"
                    if site.get("atom_chemical_role") == "functional_group_context_only" else None,
                    "No direct or conditional attachment chemistry"
                    if not chemically_supported else None,
                ] if value
            ],
        }
        atoms.append(atom)

        if serial_int is not None and site.get("candidate_attachment_atom"):
            candidate_atom_serials.append(serial_int)
        if serial_int is not None and chemically_supported:
            chemically_supported_atom_serials.append(serial_int)
        if serial_int is not None and tier_short in {"High", "Moderate"}:
            priority_atom_serials.append(serial_int)
        if serial_int is not None and site.get("high_priority_attachment_atom"):
            high_priority_atom_serials.append(serial_int)
        if serial_int is not None and site.get("solvent_exposed"):
            surface_atom_serials.append(serial_int)

        for smiles_index in smiles_indices:
            existing = graph_nodes_by_index.get(smiles_index)
            candidate = {
                "smiles_atom_index": smiles_index,
                "element": site.get("element"),
                "pdb_atom_serial": serial_int,
                "pdb_atom_name": site.get("exact_atom"),
                "region_id": None,
                "candidate_attachment_flag": int(bool(site.get("candidate_attachment_atom"))),
                "surface_defining_flag": int(bool(site.get("solvent_exposed"))),
                "attachment_score": site.get("attachment_priority_score"),
                "confidence": tier_short,
                "atom_chemical_role": site.get("atom_chemical_role"),
                "direct_attachment_support": site.get("direct_attachment_support"),
                "conditional_substitution_support": site.get("conditional_substitution_support"),
            }
            if existing is None or _numeric_value(candidate.get("attachment_score")) > _numeric_value(existing.get("attachment_score")):
                graph_nodes_by_index[smiles_index] = candidate

    # Display-site grouping is intentionally performed after the atom payload
    # is built.  It assigns a presentation-only label to each candidate atom;
    # the underlying Stage 13 values and the empty legacy-region payload stay
    # unchanged.
    display_site_clusters = _attachment_display_site_clusters(atoms)
    summary["attachment_display_site_count"] = len(display_site_clusters)
    graph = _attachment_graph_from_canonical_smiles(ligand_instance_id, atoms)

    logging.info(
        "[attachment-v2.6] ligand_instance_id=%s pdb=%s ligand=%s candidates=%s chemically_supported=%s priority=%s high=%s",
        ligand_instance_id,
        summary_dict.get("pdb_code"),
        summary_dict.get("ligand_resname"),
        len(candidate_atom_serials),
        len(chemically_supported_atom_serials),
        len(priority_atom_serials),
        len(high_priority_atom_serials),
    )

    return {
        "data_available": True,
        "summary": summary,
        "regions": [],
        "display_site_clusters": display_site_clusters,
        "display_site_grouping": {
            "method": "candidate-coordinate-neighborhood",
            "distance_cutoff_a": ATTACHMENT_DISPLAY_SITE_DISTANCE_A,
            "presentation_only": True,
        },
        "atoms": atoms,
        "candidate_atom_serials": sorted(set(candidate_atom_serials)),
        "chemically_supported_atom_serials": sorted(set(chemically_supported_atom_serials)),
        "priority_atom_serials": sorted(set(priority_atom_serials)),
        "high_priority_atom_serials": sorted(set(high_priority_atom_serials)),
        "surface_atom_serials": sorted(set(surface_atom_serials)),
        "graph": graph,
        "site_model": "atom-level-v2.6",
        "region_semantics_available": False,
        "ligand_instance_id": ligand_instance_id,
        "method_version": ATTACHMENT_METHOD_VERSION,
    }


def _normalize_remote_attachment_site_display(payload):
    """Apply the presentation-only bonded/SASA grouping to Randy atom payloads.

    Randy is intentionally RDKit-free; it sends the deposited coordinates and
    canonical SMILES.  The public app then uses the same ≤2-bond clustering
    rule as local detail routes, so both paths show identical display regions.
    """
    attachment = (payload or {}).get("attachment_sites")
    if not isinstance(attachment, dict) or not attachment.get("data_available"):
        return payload
    atoms = attachment.get("atoms") or []
    clusters = _attachment_display_site_clusters(atoms)
    attachment["display_site_clusters"] = clusters
    summary = attachment.setdefault("summary", {})
    summary["attachment_display_site_count"] = len(clusters)
    attachment["display_site_grouping"] = {
        "method": "candidate-coordinate-and-two-bond-neighborhood",
        "distance_cutoff_a": ATTACHMENT_DISPLAY_SITE_DISTANCE_A,
        "presentation_only": True,
    }
    return payload

def _merge_optional_protacability_data(rows, readiness_rows=None, warhead_rows=None, attachment_rows=None):
    readiness_indexes = _build_readiness_indexes(readiness_rows or []) if readiness_rows else None
    warhead_indexes = _build_warhead_indexes(warhead_rows or []) if warhead_rows else None
    attachment_index = _build_attachment_index(attachment_rows or []) if attachment_rows else None

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

        attachment_match = attachment_index.get(_normalize_attachment_key(current)) if attachment_index else None
        current.update(_attachment_summary_from_match(attachment_match))
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


def _decorate_protacability_rows(rows, collapse_labels=True, readiness_rows=None, warhead_rows=None, attachment_rows=None):
    decorated = _merge_optional_protacability_data(
        rows,
        readiness_rows=readiness_rows,
        warhead_rows=warhead_rows,
        attachment_rows=attachment_rows,
    )
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
    canonical_target_filter = set(filters.get("canonical_target_ids") or [])
    tier_filter = set(filters["tiers"])
    warhead_tier_filter = set(filters["warhead_tiers"])
    readiness_tier_filter = set(filters["readiness_tiers"])
    evidence_level_filter = set(filters["evidence_levels"])
    smiles_source_filter = set(filters["smiles_sources"])
    ligand_query = (filters["ligand"] or "").strip().upper()
    ligand_instance_query = (filters.get("ligand_instance_id") or "").strip()
    pdb_query = (filters["pdb_code"] or "").strip().upper()

    filtered = []
    for row in rows:
        if ligand_instance_query and str(row.get("ligand_instance_id") or "") != ligand_instance_query:
            continue
        if virus_filter and row.get("virus_name") not in virus_filter:
            continue
        if protein_filter and row.get(protein_field) not in protein_filter:
            continue
        if canonical_target_filter and row.get("canonical_target_id") not in canonical_target_filter:
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
        if filters["has_attachment_sites"] is not None and int(bool(row.get("has_attachment_site_evidence"))) != filters["has_attachment_sites"]:
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
        "has_attachment_site_evidence",
        "has_candidate_attachment_regions",
        "attachment_analysis_id",
        "attachment_method_version",
        "attachment_analysis_status",
        "attachment_eligibility_status",
        "attachment_mapping_status",
        "attachment_region_count",
        "attachment_candidate_atom_count",
        "attachment_display_site_count",
        "best_attachment_score",
        "best_attachment_confidence",
        "attachment_instance_resolution_status",
        "attachment_instance_ambiguity_flag",
        "attachment_site_model",
        "attachment_region_semantics_available",
        "attachment_exposed_mapped_atom_count",
        "attachment_outward_supported_candidate_count",
        "attachment_clear_exit_supported_candidate_count",
        "attachment_chemically_supported_candidate_count",
        "attachment_direct_candidate_count",
        "attachment_conditional_candidate_count",
        "attachment_conditional_clear_exit_candidate_count",
        "attachment_high_priority_atom_count",
        "attachment_high_priority_direct_atom_count",
        "best_attachment_priority_tier",
        "best_attachment_atom",
        "warhead_evidence_score",
        "readiness_rank_score",
    ]
    return {field: row.get(field) for field in fields}


def _row_has_mapped_exposed_warhead_evidence(row):
    explicit_flag = row.get("has_mapped_exposed_warhead_evidence")
    if explicit_flag not in (None, ""):
        return _has_positive_value(explicit_flag)
    return (
        _has_positive_value(row.get("solvent_exposed_mapped_atom_count"))
        and _has_positive_value(row.get("pdb_to_smiles_mapped_atom_count"))
    )


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
            "mapped_exposed_chain_count": sum(1 for row in group_rows if _row_has_mapped_exposed_warhead_evidence(row)),
            "has_mapped_exposed_warhead_evidence": int(
                any(_row_has_mapped_exposed_warhead_evidence(row) for row in group_rows)
            ),
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
            "mapped_exposed_structure_count": sum(
                1 for row in group_rows if _row_has_mapped_exposed_warhead_evidence(row)
            ),
            "has_mapped_exposed_warhead_evidence": int(
                any(_row_has_mapped_exposed_warhead_evidence(row) for row in group_rows)
            ),
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
    # Target Browser identity is the canonical target ID scoped to the virus.
    # The display name is deliberately kept separate: legacy source labels and
    # canonical names are not safe grouping keys.
    canonical_mode = any(row.get("canonical_target_id") for row in rows)
    canonical_metadata = {}
    grouping_rows = rows
    if canonical_mode:
        grouping_rows = []
        for source_row in rows:
            row = dict(source_row)
            canonical_target_id = row.get("canonical_target_id")
            if not canonical_target_id:
                # A canonical Target Browser response must never quietly fall
                # back to a legacy protein label.
                continue
            key = (row.get("virus_name"), canonical_target_id)
            metadata = canonical_metadata.setdefault(key, {
                "canonical_target_id": canonical_target_id,
                "canonical_target_name": row.get("canonical_target_name") or row.get("protein_type"),
                "target_family": row.get("target_family"),
                "entity_role": row.get("entity_role"),
                "source_protein_types": set(),
            })
            if row.get("source_protein_type"):
                metadata["source_protein_types"].add(row["source_protein_type"])
            # Reuse the established structure/protein aggregation with the
            # stable identifier, then restore the human-facing name below.
            row["protein_type"] = canonical_target_id
            row["display_protein_type"] = canonical_target_id
            grouping_rows.append(row)

    protein_rows = _group_protein_rows(grouping_rows)
    structure_rows = _group_structure_rows(grouping_rows)
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

        metadata = canonical_metadata.get((virus_name, protein_type), {}) if canonical_mode else {}
        display_protein_type = metadata.get("canonical_target_name") or protein_type
        row = {
            "view_type": "targets",
            "virus_name": virus_name,
            "protein_type": display_protein_type,
            "canonical_target_id": metadata.get("canonical_target_id"),
            "canonical_target_name": metadata.get("canonical_target_name"),
            "target_family": metadata.get("target_family"),
            "entity_role": metadata.get("entity_role"),
            "source_protein_types": sorted(metadata.get("source_protein_types") or []),
            "source_protein_type": "; ".join(sorted(metadata.get("source_protein_types") or [])),
            "target_key": f"{virus_name}::{metadata.get('canonical_target_id') or protein_type}",
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
            "mapped_exposed_structure_count": sum(
                1 for row in target_structures if _row_has_mapped_exposed_warhead_evidence(row)
            ),
            "has_mapped_exposed_warhead_evidence": int(
                any(_row_has_mapped_exposed_warhead_evidence(row) for row in target_structures)
            ),
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
                1 for row in rows if _row_has_mapped_exposed_warhead_evidence(row)
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
                1 for row in rows if _row_has_mapped_exposed_warhead_evidence(row)
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
                1 for row in rows if _row_has_mapped_exposed_warhead_evidence(row)
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
            1 for row in rows if _row_has_mapped_exposed_warhead_evidence(row)
        ),
        "high_warhead_linkability_ligands": sum(1 for row in rows if _numeric_value(row.get("warhead_linkability_score")) >= 70),
        "high_degrader_readiness_targets": sum(1 for row in rows if _numeric_value(row.get("degrader_design_readiness_score")) >= 70),
        "average_degrader_readiness": _round_or_none(sum(score_values) / len(score_values), 2) if score_values else None,
        "total_rows": len(rows),
    }


def _load_protacability_assessment_rows(conn, pdb_code=None, virus_name=None, protein_type=None, ligand=None, ligand_instance_id=None):
    # The table retains older releases for provenance.  Browser ranking must
    # read the final v2.8 release only, not decorate every historical row.
    query = "SELECT * FROM protacability_assessment WHERE method_version = ?"
    params = [PROTACABILITY_METHOD_VERSION]
    if pdb_code:
        query += " AND pdb_code = ?"
        params.append(pdb_code)
    if virus_name:
        query += " AND virus_name = ?"
        params.append(virus_name)
    if protein_type:
        query += " AND protein_type = ?"
        params.append(protein_type)
    if ligand:
        query += """
            AND EXISTS (
                SELECT 1
                FROM protacability_ligand_inventory AS inventory
                WHERE inventory.ligand_instance_id = protacability_assessment.ligand_instance_id
                  AND UPPER(inventory.ligand_resname) = UPPER(?)
            )
        """
        params.append(ligand)
    if ligand_instance_id is not None:
        query += " AND ligand_instance_id = ?"
        params.append(ligand_instance_id)
    return conn.execute(query, params).fetchall()


def _resolve_protacability_occurrence_context(conn, args):
    """Resolve and validate an occurrence deep link before searching assessment rows."""
    raw_occurrence_id = (args.get("ligand_instance_id") or "").strip()
    if not raw_occurrence_id:
        return None, None
    try:
        occurrence_id = int(raw_occurrence_id)
    except ValueError:
        return None, "The requested ligand occurrence could not be found."

    occurrence = conn.execute(
        """
        SELECT ligand_instance_id, pdb_code, model_id, ligand_resname,
               ligand_chain, ligand_residue_id, ligand_insertion_code
        FROM protacability_ligand_inventory
        WHERE ligand_instance_id = ?
        ORDER BY inventory_id DESC
        LIMIT 1
        """,
        (occurrence_id,),
    ).fetchone()
    if occurrence is None:
        return None, "The requested ligand occurrence could not be found."

    context = dict(occurrence)
    supplied_values = {
        "ligand": (args.get("ligand") or "").strip().upper(),
        "pdb_code": (args.get("pdb_code") or "").strip().upper(),
        "ligand_chain": (args.get("ligand_chain") or "").strip(),
        "ligand_residue_id": (args.get("ligand_residue_id") or "").strip(),
        "model_id": (args.get("model_id") or "").strip(),
    }
    canonical_values = {
        "ligand": str(context["ligand_resname"] or "").upper(),
        "pdb_code": str(context["pdb_code"] or "").upper(),
        "ligand_chain": str(context["ligand_chain"] or ""),
        "ligand_residue_id": str(context["ligand_residue_id"] or ""),
        "model_id": str(context["model_id"] or ""),
    }
    for key, supplied_value in supplied_values.items():
        if supplied_value and supplied_value != canonical_values[key]:
            return None, "The supplied ligand occurrence details do not match the current database record."
    return context, None


def _load_protacability_enrichment_tables(conn):
    optional_tables = _protacability_optional_table_names(conn)
    readiness_rows = _load_optional_table_rows(
        conn,
        "protacability_degrader_readiness",
        PROTACABILITY_METHOD_VERSION,
    ) if "protacability_degrader_readiness" in optional_tables else []
    warhead_rows = _load_optional_table_rows(
        conn,
        "protacability_warhead_linkability",
        PROTACABILITY_METHOD_VERSION,
    ) if "protacability_warhead_linkability" in optional_tables else []
    attachment_rows = _load_attachment_analysis_rows(conn) if "v2_attachment_site_summary" in optional_tables else []
    return readiness_rows, warhead_rows, attachment_rows


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

    readiness_rows, warhead_rows, attachment_rows = _load_protacability_enrichment_tables(conn)
    assessment_rows = _load_protacability_assessment_rows(conn)
    if view == "targets":
        assessment_rows = _load_canonical_target_browser_assessment_rows(conn)
    rows = _decorate_protacability_rows(
        assessment_rows,
        collapse_labels=collapse_labels,
        readiness_rows=readiness_rows,
        warhead_rows=warhead_rows,
        attachment_rows=attachment_rows,
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


def _load_canonical_target_browser_assessment_rows(conn):
    """One current assessment per canonical, public ligand occurrence.

    Canonical target views own public target eligibility and identity.  The
    selected assessment remains scientific provenance and is not altered.
    """
    rows = conn.execute(
        """
        SELECT a.*, v.canonical_target_id, v.canonical_target_name,
               v.source_protein_type, v.target_family, v.entity_role
        FROM protacability_assessment AS a
        JOIN v2_target_browser_ligand_context AS v
          ON v.ligand_instance_id = a.ligand_instance_id
        WHERE a.method_version = ?
        ORDER BY a.ligand_instance_id,
                 a.protacability_proxy_score DESC,
                 a.assessment_id ASC
        """,
        (PROTACABILITY_METHOD_VERSION,),
    ).fetchall()
    selected = []
    seen_occurrences = set()
    for raw_row in rows:
        row = dict(raw_row)
        occurrence_id = row["ligand_instance_id"]
        if occurrence_id in seen_occurrences:
            continue
        seen_occurrences.add(occurrence_id)
        row["protein_type"] = row["canonical_target_name"]
        row["canonical_target_id"] = row["canonical_target_id"]
        row["source_protein_type"] = row["source_protein_type"]
        selected.append(row)
    return selected


def _canonical_target_browser_rows_from_assessments(assessment_rows):
    conn = connect_db_row()
    try:
        authority_rows = conn.execute(
            """
            SELECT ligand_instance_id, canonical_target_id, canonical_target_name,
                   source_protein_type, target_family, entity_role
            FROM v2_target_browser_ligand_context
            """
        ).fetchall()
    finally:
        conn.close()
    authority = {row["ligand_instance_id"]: dict(row) for row in authority_rows}
    selected = {}
    for raw_row in assessment_rows or []:
        row = dict(raw_row)
        occurrence_id = row.get("ligand_instance_id")
        canonical = authority.get(occurrence_id)
        if not canonical:
            continue
        prior = selected.get(occurrence_id)
        if prior and _numeric_value(prior.get("protacability_proxy_score")) >= _numeric_value(row.get("protacability_proxy_score")):
            continue
        row.update(canonical)
        row["protein_type"] = canonical["canonical_target_name"]
        selected[occurrence_id] = row
    return list(selected.values())


def _build_protacability_filter_options_payload_from_rows(assessment_rows, readiness_rows, warhead_rows, args, attachment_rows=None):
    collapse_labels = _protacability_collapse_labels(args.get("collapse_labels"))
    filters = _build_protacability_filters(args)
    if _protacability_view_mode(args.get("view")) == "targets":
        assessment_rows = _canonical_target_browser_rows_from_assessments(assessment_rows)
    rows = _decorate_protacability_rows(
        assessment_rows,
        collapse_labels=collapse_labels,
        readiness_rows=readiness_rows,
        warhead_rows=warhead_rows,
        attachment_rows=attachment_rows,
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


def _prepare_protacability_result_set_from_rows(assessment_rows, readiness_rows, warhead_rows, args, export_all=False, attachment_rows=None):
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

    if view == "targets":
        assessment_rows = _canonical_target_browser_rows_from_assessments(assessment_rows)
    rows = _decorate_protacability_rows(
        assessment_rows,
        collapse_labels=collapse_labels,
        readiness_rows=readiness_rows,
        warhead_rows=warhead_rows,
        attachment_rows=attachment_rows,
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


def _local_protacability_source_payload(*, pdb_code=None, virus_name=None, protein_type=None, ligand=None, ligand_instance_id=None, include_lysine=False, include_inventory=False, request_args=None):
    tables_ok, error_message = local_tables_available(PROTACABILITY_REQUIRED_TABLES)
    if not tables_ok:
        raise RandyBackendError(error_message)

    with connect_db_row() as conn:
        if not protacability_tables_available(conn):
            return {
                "data_available": False,
                "assessment_rows": [],
                "readiness_rows": [],
                "warhead_rows": [],
                "attachment_rows": [],
                "lysine_rows": [],
                "ligand_inventory": [],
            }

        occurrence_context = None
        if ligand_instance_id not in (None, ""):
            occurrence_context, deep_link_error = _resolve_protacability_occurrence_context(conn, request_args or {})
            if deep_link_error:
                return {
                    "data_available": True,
                    "assessment_rows": [],
                    "readiness_rows": [],
                    "warhead_rows": [],
                    "attachment_rows": [],
                    "lysine_rows": [],
                    "ligand_inventory": [],
                    "deep_link_error": deep_link_error,
                }
            ligand_instance_id = occurrence_context["ligand_instance_id"]
            ligand = occurrence_context["ligand_resname"]
            pdb_code = occurrence_context["pdb_code"]

        assessment_rows = [
            dict(row)
            for row in _load_protacability_assessment_rows(
                conn,
                pdb_code=pdb_code,
                virus_name=virus_name,
                protein_type=protein_type,
                ligand=ligand,
                ligand_instance_id=ligand_instance_id,
            )
        ]
        readiness_rows, warhead_rows, attachment_rows = _load_protacability_enrichment_tables(conn)
        lysine_rows = []
        ligand_inventory = []
        if include_lysine:
            if pdb_code:
                lysine_rows = [dict(row) for row in conn.execute("SELECT * FROM protacability_lysine_proximity WHERE pdb_code = ?", (pdb_code,)).fetchall()]
            else:
                lysine_rows = [dict(row) for row in conn.execute("SELECT * FROM protacability_lysine_proximity").fetchall()]
        if include_inventory:
            if pdb_code:
                ligand_inventory = [dict(row) for row in conn.execute("SELECT * FROM protacability_ligand_inventory WHERE pdb_code = ?", (pdb_code,)).fetchall()]
            else:
                ligand_inventory = [dict(row) for row in conn.execute("SELECT * FROM protacability_ligand_inventory").fetchall()]

    return {
        "data_available": True,
        "assessment_rows": assessment_rows,
        "readiness_rows": readiness_rows,
        "warhead_rows": warhead_rows,
        "attachment_rows": attachment_rows,
        "lysine_rows": lysine_rows,
        "ligand_inventory": ligand_inventory,
        "deep_link_context": occurrence_context,
    }


def _load_protacability_source_payload(*, pdb_code=None, virus_name=None, protein_type=None, ligand=None, ligand_instance_id=None, include_lysine=False, include_inventory=False, request_args=None):
    mode = _normalized_backend_mode()
    params = {
        "pdb_code": pdb_code or "",
        "virus_name": virus_name or "",
        "protein_type": protein_type or "",
        "ligand": ligand or "",
        "ligand_instance_id": ligand_instance_id or "",
        "include_lysine": str(bool(include_lysine)).lower(),
        "include_inventory": str(bool(include_inventory)).lower(),
    }

    if mode == "randy":
        return randy_get("protacability/source", params=params)

    if mode == "auto" and randy_available():
        try:
            return randy_get("protacability/source", params=params)
        except RandyBackendError:
            logging.warning("Falling back to local PROTACability source payload")

    return _local_protacability_source_payload(
        pdb_code=pdb_code,
        virus_name=virus_name,
        protein_type=protein_type,
        ligand=ligand,
        ligand_instance_id=ligand_instance_id,
        include_lysine=include_lysine,
        include_inventory=include_inventory,
        request_args=request_args,
    )


def _remote_protacability_get(path, *, params=None, max_bytes=10 * 1024 * 1024):
    mode = _normalized_backend_mode()
    if mode == "randy":
        return randy_get(path, params=params, max_bytes=max_bytes)
    if mode == "auto" and randy_available():
        return randy_get(path, params=params, max_bytes=max_bytes)
    raise RandyBackendError("RANDY API is not configured.", status_code=500)


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


def _attach_ligand_label_asym_ids(ligand_inventory):
    """Attach stored occurrence metadata for browser-only RCSB ligand retrieval."""
    records = [dict(record) for record in (ligand_inventory or [])]
    instance_ids = [record.get("ligand_instance_id") for record in records if record.get("ligand_instance_id") not in (None, "")]
    if not instance_ids:
        return records
    try:
        with connect_db_row() as conn:
            placeholders = ",".join("?" for _ in instance_ids)
            rows = conn.execute(
                f"SELECT ligand_instance_id, label_asym_id FROM ligand_instances WHERE ligand_instance_id IN ({placeholders})",
                instance_ids,
            ).fetchall()
        labels = {str(row["ligand_instance_id"]): row["label_asym_id"] for row in rows}
        for record in records:
            record["label_asym_id"] = labels.get(str(record.get("ligand_instance_id")), "")
    except Exception:
        # Viewer-only metadata must never block stored analysis results.
        pass
    return records


def _serialize_ligand_contexts(ligand_inventory):
    contexts = []
    for row in ligand_inventory or []:
        contexts.append({
            "ligand_instance_id": row.get("ligand_instance_id"),
            "model_id": row.get("model_id"),
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


def make_ligand_code_aliases(code):
    raw = str(code or "").strip().upper()
    if not raw:
        return []
    aliases = [raw]
    o_to_zero = raw.replace("O", "0")
    zero_to_o = raw.replace("0", "O")
    for alias in (o_to_zero, zero_to_o):
        if alias and alias not in aliases:
            aliases.append(alias)
    return aliases


def _table_exists(conn, table_name):
    row = conn.execute(
        "SELECT name FROM sqlite_master WHERE type IN ('table','view') AND name = ?",
        (table_name,)
    ).fetchone()
    return row is not None


def _table_columns(conn, table_name):
    return [row[1] for row in conn.execute(f"PRAGMA table_info({table_name})").fetchall()]


def _first_present_column(columns, candidates):
    column_set = set(columns or [])
    for candidate in candidates:
        if candidate in column_set:
            return candidate
    return None


def _normalize_debug_row(row, table_name, columns):
    row_dict = dict(row)
    norm = {
        "source_table": table_name,
        "pdb_code": row_dict.get(_first_present_column(columns, ["pdb_code", "pdb_id"])),
        "pdb_id": row_dict.get(_first_present_column(columns, ["pdb_id", "pdb_code"])),
        "ligand": row_dict.get(_first_present_column(columns, ["ligand", "ligand_resname"])),
        "ligand_resname": row_dict.get(_first_present_column(columns, ["ligand_resname", "ligand"])),
        "chain": row_dict.get(_first_present_column(columns, ["chain", "chain_id", "ligand_chain"])),
        "ligand_chain": row_dict.get(_first_present_column(columns, ["ligand_chain", "chain", "chain_id"])),
        "ligand_id": row_dict.get(_first_present_column(columns, ["ligand_id", "ligand_residue_id", "residue_id", "ligand_residue"])),
        "residue_id": row_dict.get(_first_present_column(columns, ["residue_id", "ligand_id", "ligand_residue_id"])),
        "ligand_residue_id": row_dict.get(_first_present_column(columns, ["ligand_residue_id", "ligand_id", "residue_id"])),
        "atom_id": row_dict.get(_first_present_column(columns, ["atom_id"])),
        "exact_atom": row_dict.get(_first_present_column(columns, ["exact_atom"])),
        "smiles_atom_index": row_dict.get(_first_present_column(columns, ["smiles_atom_index"])),
        "residue_number": row_dict.get(_first_present_column(columns, ["residue_number", "ligand_residue_id", "residue_id", "ligand_id"])),
        "residue_chain": row_dict.get(_first_present_column(columns, ["residue_chain", "ligand_chain", "chain", "chain_id"])),
    }
    row_dict["normalized"] = norm
    return row_dict


def _build_debug_table_payload(conn, table_name, query_ctx):
    payload = {
        "table": table_name,
        "row_count": 0,
        "columns_present": [],
        "rows": [],
        "error": None,
    }
    try:
        if not _table_exists(conn, table_name):
            payload["error"] = "table does not exist"
            return payload
        columns = _table_columns(conn, table_name)
        payload["columns_present"] = columns
        where_clauses = []
        params = []

        pdb_col = _first_present_column(columns, ["pdb_code", "pdb_id"])
        if pdb_col and query_ctx["pdb_code_upper"]:
            where_clauses.append(f"UPPER(CAST({pdb_col} AS TEXT)) = ?")
            params.append(query_ctx["pdb_code_upper"])

        virus_col = _first_present_column(columns, ["virus_name"])
        if virus_col and query_ctx["virus_name"]:
            where_clauses.append(f"CAST({virus_col} AS TEXT) = ?")
            params.append(query_ctx["virus_name"])

        protein_col = _first_present_column(columns, ["protein_type", "protein"])
        if protein_col and query_ctx["protein_type"]:
            where_clauses.append(f"CAST({protein_col} AS TEXT) = ?")
            params.append(query_ctx["protein_type"])

        ligand_col = _first_present_column(columns, ["ligand", "ligand_resname"])
        if ligand_col and query_ctx["ligand_alias_candidates"]:
            placeholders = ",".join(["?"] * len(query_ctx["ligand_alias_candidates"]))
            where_clauses.append(f"UPPER(CAST({ligand_col} AS TEXT)) IN ({placeholders})")
            params.extend(query_ctx["ligand_alias_candidates"])

        chain_col = _first_present_column(columns, ["ligand_chain", "chain", "chain_id"])
        if chain_col and query_ctx["chain_upper"]:
            where_clauses.append(f"UPPER(CAST({chain_col} AS TEXT)) = ?")
            params.append(query_ctx["chain_upper"])

        resid_col = _first_present_column(columns, ["ligand_residue_id", "ligand_id", "residue_id"])
        resid_clauses = []
        resid_values = []
        if resid_col and query_ctx["residue_id_string"]:
            resid_clauses.append(f"CAST({resid_col} AS TEXT) = ?")
            resid_values.append(query_ctx["residue_id_string"])
        if resid_col and query_ctx["residue_id_int_if_possible"] is not None:
            resid_clauses.append(f"{resid_col} = ?")
            resid_values.append(query_ctx["residue_id_int_if_possible"])
        if resid_clauses:
            where_clauses.append("(" + " OR ".join(resid_clauses) + ")")
            params.extend(resid_values)

        where_sql = f" WHERE {' AND '.join(where_clauses)}" if where_clauses else ""
        limit = 100000 if query_ctx["debug_full"] else 50
        query = f"SELECT * FROM {table_name}{where_sql} LIMIT {limit}"
        rows = conn.execute(query, params).fetchall()
        payload["row_count"] = len(rows)
        payload["rows"] = [_normalize_debug_row(row, table_name, columns) for row in rows]
        return payload
    except Exception as exc:
        payload["error"] = str(exc)
        return payload


def _parse_pdb_coordinate_text(text):
    summary = defaultdict(int)
    matching_lines_preview = []
    if not text:
        return summary, matching_lines_preview
    for line in text.splitlines():
        rec = line[:6].strip().upper()
        if rec not in {"HETATM", "ATOM"}:
            continue
        resname = line[17:20].strip().upper() if len(line) >= 20 else ""
        chain = line[21:22].strip().upper() if len(line) >= 22 else ""
        resno = line[22:27].strip().upper() if len(line) >= 27 else ""
        if rec == "HETATM":
            summary[(resname, chain, resno)] += 1
        if len(matching_lines_preview) < 20:
            matching_lines_preview.append(line[:160])
    return summary, matching_lines_preview


def _parse_cif_coordinate_text(text):
    summary = defaultdict(int)
    matching_lines_preview = []
    if not text:
        return summary, matching_lines_preview
    headers = []
    in_atom_loop = False
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.lower() == "loop_":
            in_atom_loop = True
            headers = []
            continue
        if in_atom_loop and line.startswith("_atom_site."):
            headers.append(line)
            continue
        if in_atom_loop and headers and not line.startswith("_"):
            parts = line.split()
            if len(parts) < len(headers):
                continue
            header_map = {headers[i]: parts[i] for i in range(len(headers))}
            group_pdb = str(header_map.get("_atom_site.group_PDB", "")).upper()
            if group_pdb != "HETATM":
                continue
            resname = str(
                header_map.get("_atom_site.auth_comp_id")
                or header_map.get("_atom_site.label_comp_id")
                or ""
            ).upper()
            chain = str(
                header_map.get("_atom_site.auth_asym_id")
                or header_map.get("_atom_site.label_asym_id")
                or ""
            ).upper()
            resno = str(
                header_map.get("_atom_site.auth_seq_id")
                or header_map.get("_atom_site.label_seq_id")
                or ""
            ).upper()
            summary[(resname, chain, resno)] += 1
            if len(matching_lines_preview) < 20:
                matching_lines_preview.append(raw_line[:160])
    return summary, matching_lines_preview


def _scan_coordinate_source(source_name, text, ligand_aliases, chain=None, residue_id=None):
    chain_upper = str(chain or "").strip().upper()
    residue_text = str(residue_id or "").strip().upper()
    residue_digits = "".join(ch for ch in residue_text if ch.isdigit())
    contains_requested_ligand_text = any(alias in text.upper() for alias in ligand_aliases)
    contains_requested_ligand_hetatm = False
    residue_chain_match = False
    if source_name.lower().endswith(".cif"):
        het_summary, preview = _parse_cif_coordinate_text(text)
    else:
        het_summary, preview = _parse_pdb_coordinate_text(text)
    hetatm_residue_summary = []
    for (resname, res_chain, resno), atom_count in sorted(het_summary.items()):
        if resname in ligand_aliases:
            contains_requested_ligand_hetatm = True
        if chain_upper and residue_digits:
            this_digits = "".join(ch for ch in str(resno) if ch.isdigit())
            if res_chain == chain_upper and this_digits == residue_digits:
                residue_chain_match = True
        hetatm_residue_summary.append({
            "resname": resname,
            "chain": res_chain,
            "resno": resno,
            "atom_count": atom_count
        })
    matching_lines_preview = [line for line in preview if any(alias in line.upper() for alias in ligand_aliases)][:20]
    return {
        "source": source_name,
        "exists": True,
        "contains_requested_ligand_text": contains_requested_ligand_text,
        "contains_requested_ligand_hetatm": contains_requested_ligand_hetatm,
        "contains_chain_residue_match": residue_chain_match,
        "hetatm_residue_summary": hetatm_residue_summary,
        "matching_lines_preview": matching_lines_preview,
    }


def is_pdb_text(text):
    if not text:
        return False
    for line in text.splitlines():
        if line.startswith("ATOM  ") or line.startswith("HETATM"):
            return True
    return False


def pdb_has_protein_atoms(text):
    if not text:
        return False
    return any(line.startswith("ATOM  ") for line in text.splitlines())


def summarize_pdb_hetatm(text):
    summary = defaultdict(int)
    if not text:
        return []
    for line in text.splitlines():
        if not line.startswith("HETATM"):
            continue
        resname = line[17:20].strip().upper() if len(line) >= 20 else ""
        chain = line[21:22].strip().upper() if len(line) >= 22 else ""
        resno = line[22:26].strip().upper() if len(line) >= 26 else ""
        summary[(resname, chain, resno)] += 1
    return [
        {"resname": r, "chain": c, "resno": n, "atom_count": ct}
        for (r, c, n), ct in sorted(summary.items())
    ]


def pdb_contains_ligand(text, ligand_code, chain=None, residue_id=None):
    ligand = str(ligand_code or "").strip().upper()
    if not ligand:
        return False
    chain_upper = str(chain or "").strip().upper()
    residue_upper = str(residue_id or "").strip().upper()
    residue_digits = "".join(ch for ch in residue_upper if ch.isdigit())
    for line in text.splitlines():
        if not line.startswith("HETATM"):
            continue
        resname = line[17:20].strip().upper() if len(line) >= 20 else ""
        if resname != ligand:
            continue
        if chain_upper:
            line_chain = line[21:22].strip().upper() if len(line) >= 22 else ""
            if line_chain != chain_upper:
                continue
        if residue_digits:
            line_res = line[22:26].strip().upper() if len(line) >= 26 else ""
            if "".join(ch for ch in line_res if ch.isdigit()) != residue_digits:
                continue
        return True
    return False


def pdb_ligand_hetatm_count(text, ligand_code, chain=None, residue_id=None):
    ligand = str(ligand_code or "").strip().upper()
    if not ligand:
        return 0
    chain_upper = str(chain or "").strip().upper()
    residue_upper = str(residue_id or "").strip().upper()
    residue_digits = "".join(ch for ch in residue_upper if ch.isdigit())
    count = 0
    for line in text.splitlines():
        if not line.startswith("HETATM"):
            continue
        resname = line[17:20].strip().upper() if len(line) >= 20 else ""
        if resname != ligand:
            continue
        if chain_upper:
            line_chain = line[21:22].strip().upper() if len(line) >= 22 else ""
            if line_chain != chain_upper:
                continue
        if residue_digits:
            line_res = line[22:26].strip().upper() if len(line) >= 26 else ""
            if "".join(ch for ch in line_res if ch.isdigit()) != residue_digits:
                continue
        count += 1
    return count


def pdb_first_hetatm_lines(text, limit=20):
    lines = []
    for line in text.splitlines():
        if line.startswith("HETATM"):
            lines.append(line)
            if len(lines) >= limit:
                break
    return lines


def pdb_matching_ligand_lines(text, ligand_code, chain=None, residue_id=None, limit=50):
    ligand = str(ligand_code or "").strip().upper()
    chain_upper = str(chain or "").strip().upper()
    residue_upper = str(residue_id or "").strip().upper()
    residue_digits = "".join(ch for ch in residue_upper if ch.isdigit())
    matches = []
    for line in text.splitlines():
        if not line.startswith("HETATM"):
            continue
        resname = line[17:20].strip().upper() if len(line) >= 20 else ""
        if ligand and resname != ligand:
            continue
        if chain_upper:
            line_chain = line[21:22].strip().upper() if len(line) >= 22 else ""
            if line_chain != chain_upper:
                continue
        if residue_digits:
            line_res = line[22:26].strip().upper() if len(line) >= 26 else ""
            if "".join(ch for ch in line_res if ch.isdigit()) != residue_digits:
                continue
        matches.append(line)
        if len(matches) >= limit:
            break
    return matches


def _candidate_coordinate_files(pdb_code, extensions=None):
    pdb_upper = str(pdb_code or "").strip().upper()
    candidates = []
    extensions = extensions or [".pdb", ".cif", ".mmcif"]
    search_dirs = [
        os.getcwd(),
        os.path.join(os.getcwd(), "PDB_FILES"),
        os.path.join(os.getcwd(), "static"),
        os.path.join(os.getcwd(), "output_files"),
        os.path.join(os.getcwd(), "pml_sessions"),
        os.path.join(os.getcwd(), "temp"),
    ]
    patterns = []
    for ext in extensions:
        ext = ext.lower()
        patterns.extend([
            f"{pdb_upper}{ext}",
            f"{pdb_upper.lower()}{ext}",
            f"*{pdb_upper}*{ext}",
            f"*{pdb_upper.lower()}*{ext}",
        ])
    for base in search_dirs:
        if not os.path.isdir(base):
            continue
        for pattern in patterns:
            for match in glob.glob(os.path.join(base, "**", pattern), recursive=True):
                if os.path.isfile(match):
                    candidates.append(os.path.abspath(match))
    seen = set()
    deduped = []
    for item in candidates:
        if item in seen:
            continue
        seen.add(item)
        deduped.append(item)
    return deduped


def _coordinate_cache_dir():
    path = os.path.join(os.getcwd(), "static", "coordinate_cache")
    os.makedirs(path, exist_ok=True)
    return path


def _ligand_sdf_cache_dir():
    path = os.path.join(os.getcwd(), "static", "ligand_sdf_cache")
    os.makedirs(path, exist_ok=True)
    return path


def _cache_pdb_path(pdb_code, ligand_code=None, chain=None, residue_id=None):
    base = [str(pdb_code or "").strip().upper()]
    if ligand_code:
        base.append(str(ligand_code or "").strip().upper())
    if chain:
        base.append(str(chain or "").strip().upper())
    if residue_id:
        base.append(str(residue_id or "").strip().upper())
    safe = "_".join(re.sub(r"[^A-Z0-9]+", "", token) for token in base if token)
    safe = safe or "COORD"
    return os.path.join(_coordinate_cache_dir(), f"{safe}.pdb")


def _cache_sdf_path(pdb_code, ligand_code=None, chain=None, residue_id=None, label_asym_id=None):
    base = [str(pdb_code or "").strip().upper()]
    if ligand_code:
        base.append(str(ligand_code or "").strip().upper())
    if chain:
        base.append(str(chain or "").strip().upper())
    if residue_id:
        base.append(str(residue_id or "").strip().upper())
    if label_asym_id:
        base.append(str(label_asym_id or "").strip().upper())
    safe = "_".join(re.sub(r"[^A-Z0-9]+", "", token) for token in base if token)
    safe = safe or "LIGAND"
    return os.path.join(_ligand_sdf_cache_dir(), f"{safe}.sdf")


def _convert_cif_to_pdb(cif_path, output_pdb_path):
    if MMCIFParser is None or PDBIO is None:
        raise RuntimeError("Biopython MMCIFParser/PDBIO is unavailable")
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("structure", cif_path)
    io_writer = PDBIO()
    io_writer.set_structure(structure)
    io_writer.save(output_pdb_path)
    return output_pdb_path


def _fetch_url_text(url, timeout=15):
    with urllib.request.urlopen(url, timeout=timeout) as resp:
        raw = resp.read()
    return raw.decode("utf-8", errors="ignore")


def _fetch_url_bytes(url, timeout=15):
    with urllib.request.urlopen(url, timeout=timeout) as resp:
        return resp.read()


def _extract_cif_ligand_instances(text, ligand_code=None, auth_chain=None, auth_seq_id=None):
    ligand_upper = str(ligand_code or "").strip().upper()
    auth_chain_upper = str(auth_chain or "").strip().upper()
    auth_seq_upper = str(auth_seq_id or "").strip().upper()
    residue_digits = "".join(ch for ch in auth_seq_upper if ch.isdigit())
    headers = []
    in_atom_loop = False
    instance_counts = defaultdict(int)

    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.lower() == "loop_":
            in_atom_loop = True
            headers = []
            continue
        if in_atom_loop and line.startswith("_atom_site."):
            headers.append(line)
            continue
        if in_atom_loop and headers and not line.startswith("_"):
            parts = line.split()
            if len(parts) < len(headers):
                continue
            header_map = {headers[i]: parts[i] for i in range(len(headers))}
            group_pdb = str(header_map.get("_atom_site.group_PDB", "")).upper()
            if group_pdb != "HETATM":
                continue

            auth_comp_id = str(header_map.get("_atom_site.auth_comp_id") or header_map.get("_atom_site.label_comp_id") or "").upper()
            auth_asym_id = str(header_map.get("_atom_site.auth_asym_id") or "").upper()
            auth_seq = str(header_map.get("_atom_site.auth_seq_id") or "").upper()
            label_asym_id = str(header_map.get("_atom_site.label_asym_id") or "").upper()

            if ligand_upper and auth_comp_id != ligand_upper:
                continue
            if auth_chain_upper and auth_asym_id != auth_chain_upper:
                continue
            if residue_digits:
                this_digits = "".join(ch for ch in auth_seq if ch.isdigit())
                if this_digits != residue_digits:
                    continue

            key = (
                auth_comp_id,
                auth_asym_id,
                auth_seq,
                label_asym_id
            )
            instance_counts[key] += 1

    instances = []
    for (comp_id, auth_asym_id, auth_seq, label_asym_id), atom_count in sorted(instance_counts.items()):
        instances.append({
            "ligand_code": comp_id,
            "auth_asym_id": auth_asym_id,
            "auth_seq_id": auth_seq,
            "label_asym_id": label_asym_id,
            "atom_count": atom_count,
        })
    return instances


def _collect_ligand_mapping_db_rows(conn, pdb_code, ligand_code, auth_chain=None, auth_seq_id=None):
    query_ctx = {
        "pdb_code_upper": str(pdb_code or "").strip().upper(),
        "ligand_alias_candidates": make_ligand_code_aliases(ligand_code),
        "chain_upper": str(auth_chain or "").strip().upper(),
        "residue_id_string": str(auth_seq_id or "").strip(),
        "residue_id_int_if_possible": None,
        "virus_name": "",
        "protein_type": "",
        "debug_full": False,
    }
    if query_ctx["residue_id_string"]:
        try:
            query_ctx["residue_id_int_if_possible"] = int(query_ctx["residue_id_string"])
        except (TypeError, ValueError):
            query_ctx["residue_id_int_if_possible"] = None

    tables = [
        "protacability_ligand_inventory",
        "Ligand_Arp_Diagram",
        "ligand_atoms",
        "Ligand_Atoms_Smiles",
        "SMILES_MAP_PDB",
        "Arpeggio_Contacts_Data",
    ]
    rows = []
    for table_name in tables:
        payload = _build_debug_table_payload(conn, table_name, query_ctx)
        if payload.get("row_count"):
            rows.extend(payload.get("rows", []))
    return rows


def _build_ligand_instance_sdf_url(pdb_code, ligand_code, auth_seq_id, label_asym_id):
    pdb_lower = str(pdb_code or "").strip().lower()
    ligand_upper = str(ligand_code or "").strip().upper()
    label_upper = str(label_asym_id or "").strip().upper()
    auth_seq_text = str(auth_seq_id or "").strip()
    return (
        f"https://models.rcsb.org/v1/{pdb_lower}/ligand"
        f"?auth_seq_id={auth_seq_text}&label_asym_id={label_upper}"
        f"&encoding=sdf&filename={pdb_lower}_{label_upper}_{ligand_upper}.sdf"
    )


def _resolve_ligand_instance_mapping(pdb_code, ligand_code, auth_chain=None, auth_seq_id=None):
    pdb_upper = str(pdb_code or "").strip().upper()
    ligand_upper = str(ligand_code or "").strip().upper()
    auth_chain_upper = str(auth_chain or "").strip().upper()
    auth_seq_text = str(auth_seq_id or "").strip()
    cache_key = (pdb_upper, ligand_upper, auth_chain_upper, auth_seq_text)
    if cache_key in LIGAND_INSTANCE_MAPPING_CACHE:
        return dict(LIGAND_INSTANCE_MAPPING_CACHE[cache_key])

    diagnostics = {
        "success": False,
        "pdb_code": pdb_upper,
        "ligand_code": ligand_upper,
        "auth_chain": auth_chain_upper,
        "auth_seq_id": auth_seq_text,
        "db_rows_found": [],
        "local_cif_instances_found": [],
        "rcsb_cif_instances_found": [],
        "chosen_label_asym_id": None,
        "sdf_url": None,
        "source": None,
    }

    mode = _normalized_backend_mode()
    if mode != "randy":
        conn = connect_db_row()
        try:
            diagnostics["db_rows_found"] = _collect_ligand_mapping_db_rows(
                conn,
                pdb_upper,
                ligand_upper,
                auth_chain=auth_chain_upper,
                auth_seq_id=auth_seq_text
            )
        finally:
            conn.close()

    local_cif_paths = _candidate_coordinate_files(pdb_upper, extensions=[".cif", ".mmcif"])
    for path in local_cif_paths:
        try:
            text = open(path, "r", encoding="utf-8", errors="ignore").read()
            matches = _extract_cif_ligand_instances(
                text,
                ligand_code=ligand_upper,
                auth_chain=auth_chain_upper,
                auth_seq_id=auth_seq_text
            )
            if matches:
                diagnostics["local_cif_instances_found"].append({
                    "source": path,
                    "instances": matches
                })
        except Exception as exc:
            diagnostics["local_cif_instances_found"].append({
                "source": path,
                "error": str(exc),
                "instances": []
            })

    chosen = None
    if diagnostics["local_cif_instances_found"]:
        for entry in diagnostics["local_cif_instances_found"]:
            if entry.get("instances"):
                chosen = dict(entry["instances"][0])
                diagnostics["source"] = "local_mmcif_atom_site"
                diagnostics["source_file"] = entry.get("source")
                break

    if not chosen:
        remote_cif_url = f"https://files.rcsb.org/download/{pdb_upper}.cif"
        try:
            cif_text = _fetch_url_text(remote_cif_url)
            remote_instances = _extract_cif_ligand_instances(
                cif_text,
                ligand_code=ligand_upper,
                auth_chain=auth_chain_upper,
                auth_seq_id=auth_seq_text
            )
            diagnostics["rcsb_cif_instances_found"] = remote_instances
            if remote_instances:
                chosen = dict(remote_instances[0])
                diagnostics["source"] = "rcsb_mmcif_atom_site"
                diagnostics["source_file"] = remote_cif_url
        except Exception as exc:
            diagnostics["rcsb_cif_error"] = str(exc)

    if chosen and chosen.get("label_asym_id"):
        diagnostics["success"] = True
        diagnostics["chosen_label_asym_id"] = chosen["label_asym_id"]
        diagnostics["atom_count"] = chosen.get("atom_count")
        diagnostics["sdf_url"] = _build_ligand_instance_sdf_url(
            pdb_upper,
            ligand_upper,
            auth_seq_text,
            chosen["label_asym_id"]
        )

    LIGAND_INSTANCE_MAPPING_CACHE[cache_key] = dict(diagnostics)
    return diagnostics


def _count_sdf_atoms(sdf_text):
    lines = str(sdf_text or "").splitlines()
    if len(lines) < 4:
        return 0
    counts_line = lines[3]
    try:
        return int(counts_line[0:3].strip())
    except Exception:
        parts = counts_line.split()
        try:
            return int(parts[0])
        except Exception:
            return 0


def _protein_only_pdb_text(text):
    allowed = {"HEADER", "TITLE", "COMPND", "SOURCE", "KEYWDS", "EXPDTA", "AUTHOR", "REMARK", "DBREF", "SEQRES", "HELIX", "SHEET", "TURN", "SSBOND", "CRYST1", "MODEL", "ATOM", "ANISOU", "TER", "ENDMDL", "END"}
    lines = []
    for raw_line in str(text or "").splitlines():
        record = raw_line[:6].strip().upper()
        if record in allowed:
            lines.append(raw_line)
    if not lines or (lines and lines[-1].strip().upper() != "END"):
        lines.append("END")
    return "\n".join(lines) + "\n"


def _resolve_coordinate_pdb(pdb_code, ligand_code=None, chain=None, residue_id=None):
    pdb_upper = str(pdb_code or "").strip().upper()
    ligand_upper = str(ligand_code or "").strip().upper()
    chain_upper = str(chain or "").strip().upper()
    residue_upper = str(residue_id or "").strip().upper()
    diagnostics = {
        "requested": {
            "pdb": pdb_upper,
            "ligand": ligand_upper,
            "chain": chain_upper,
            "residue_id": residue_upper,
        },
        "candidate_evaluations": [],
        "local_pdb_candidates": [],
        "local_cif_candidates": [],
        "remote_attempts": [],
        "selected_source": None,
        "selected_source_format": None,
        "converted_to_pdb": None,
        "served_pdb_path": None,
        "served_pdb_url": None,
        "contains_ligand": False,
        "requested_ligand_hetatm_count": 0,
        "has_protein_atoms": False,
        "hetatm_summary": [],
    }

    def _validate_pdb_text(text):
        has_atom_records = pdb_has_protein_atoms(text)
        ligand_atom_count = pdb_ligand_hetatm_count(text, ligand_upper, chain=chain_upper, residue_id=residue_upper) if ligand_upper else 0
        contains_ligand = ligand_atom_count > 0 if ligand_upper else True
        return {
            "has_atom_records": has_atom_records,
            "contains_ligand": contains_ligand,
            "ligand_atom_count": ligand_atom_count,
            "hetatm_summary": summarize_pdb_hetatm(text),
            "usable": has_atom_records and contains_ligand
        }

    def _record_candidate(priority, source, source_format, validation, source_file_before_conversion=None, converted_to_pdb=None):
        diagnostics["candidate_evaluations"].append({
            "priority": priority,
            "source": source,
            "source_format": source_format,
            "source_file_before_conversion": source_file_before_conversion or source,
            "converted_to_pdb": converted_to_pdb,
            "has_atom_records": validation.get("has_atom_records"),
            "contains_requested_ligand_hetatm": validation.get("contains_ligand"),
            "requested_ligand_hetatm_count": validation.get("ligand_atom_count"),
            "usable": validation.get("usable"),
            "hetatm_summary": validation.get("hetatm_summary"),
        })

    for path in _candidate_coordinate_files(pdb_upper, extensions=[".pdb"]):
        try:
            text = open(path, "r", encoding="utf-8", errors="ignore").read()
            validation = _validate_pdb_text(text)
            _record_candidate("A", path, "pdb", validation)
            diagnostics["local_pdb_candidates"].append({"path": path, "usable": validation["usable"]})
            if validation["usable"]:
                diagnostics["selected_source"] = path
                diagnostics["selected_source_format"] = "pdb"
                diagnostics["served_pdb_path"] = path
                diagnostics["contains_ligand"] = validation["contains_ligand"]
                diagnostics["requested_ligand_hetatm_count"] = validation["ligand_atom_count"]
                diagnostics["has_protein_atoms"] = validation["has_atom_records"]
                diagnostics["hetatm_summary"] = validation["hetatm_summary"]
                return diagnostics
        except Exception as exc:
            diagnostics["local_pdb_candidates"].append({"path": path, "usable": False, "error": str(exc)})

    for path in _candidate_coordinate_files(pdb_upper, extensions=[".cif", ".mmcif"]):
        try:
            cif_text = open(path, "r", encoding="utf-8", errors="ignore").read()
            if ligand_upper and ligand_upper not in cif_text.upper():
                diagnostics["local_cif_candidates"].append({"path": path, "has_ligand_text": False, "usable": False})
                continue
            cache_pdb = _cache_pdb_path(pdb_upper, ligand_upper, chain_upper, residue_upper)
            converted_path = _convert_cif_to_pdb(path, cache_pdb)
            pdb_text = open(converted_path, "r", encoding="utf-8", errors="ignore").read()
            validation = _validate_pdb_text(pdb_text)
            _record_candidate("B", path, "cif", validation, source_file_before_conversion=path, converted_to_pdb=converted_path)
            diagnostics["local_cif_candidates"].append({
                "path": path,
                "has_ligand_text": True,
                "converted_pdb": converted_path,
                "usable": validation["usable"]
            })
            if validation["usable"]:
                diagnostics["selected_source"] = path
                diagnostics["selected_source_format"] = "cif"
                diagnostics["converted_to_pdb"] = converted_path
                diagnostics["served_pdb_path"] = converted_path
                diagnostics["contains_ligand"] = validation["contains_ligand"]
                diagnostics["requested_ligand_hetatm_count"] = validation["ligand_atom_count"]
                diagnostics["has_protein_atoms"] = validation["has_atom_records"]
                diagnostics["hetatm_summary"] = validation["hetatm_summary"]
                return diagnostics
        except Exception as exc:
            diagnostics["local_cif_candidates"].append({"path": path, "usable": False, "error": str(exc)})

    remote_pdb_urls = [
        f"https://files.rcsb.org/view/{pdb_upper}.pdb",
        f"https://files.rcsb.org/download/{pdb_upper}.pdb",
    ]
    for url in remote_pdb_urls:
        try:
            text = _fetch_url_text(url)
            cache_pdb = _cache_pdb_path(pdb_upper, ligand_upper, chain_upper, residue_upper or "RCSB")
            with open(cache_pdb, "w", encoding="utf-8") as handle:
                handle.write(text)
            validation = _validate_pdb_text(text)
            _record_candidate("C", url, "pdb", validation, source_file_before_conversion=url, converted_to_pdb=cache_pdb)
            diagnostics["remote_attempts"].append({"url": url, "format": "pdb", "usable": validation["usable"], "cached_path": cache_pdb})
            if validation["usable"]:
                diagnostics["selected_source"] = url
                diagnostics["selected_source_format"] = "pdb"
                diagnostics["served_pdb_path"] = cache_pdb
                diagnostics["contains_ligand"] = validation["contains_ligand"]
                diagnostics["requested_ligand_hetatm_count"] = validation["ligand_atom_count"]
                diagnostics["has_protein_atoms"] = validation["has_atom_records"]
                diagnostics["hetatm_summary"] = validation["hetatm_summary"]
                return diagnostics
        except Exception as exc:
            diagnostics["remote_attempts"].append({"url": url, "format": "pdb", "usable": False, "error": str(exc)})

    remote_cif_url = f"https://files.rcsb.org/download/{pdb_upper}.cif"
    try:
        cif_text = _fetch_url_text(remote_cif_url)
        if ligand_upper and ligand_upper not in cif_text.upper():
            diagnostics["remote_attempts"].append({"url": remote_cif_url, "format": "cif", "usable": False, "reason": "ligand text missing"})
            return diagnostics
        cache_cif = os.path.join(_coordinate_cache_dir(), f"{pdb_upper}_remote.cif")
        with open(cache_cif, "w", encoding="utf-8") as handle:
            handle.write(cif_text)
        cache_pdb = _cache_pdb_path(pdb_upper, ligand_upper, chain_upper, residue_upper or "CIF")
        converted_path = _convert_cif_to_pdb(cache_cif, cache_pdb)
        pdb_text = open(converted_path, "r", encoding="utf-8", errors="ignore").read()
        validation = _validate_pdb_text(pdb_text)
        _record_candidate("D", remote_cif_url, "cif", validation, source_file_before_conversion=cache_cif, converted_to_pdb=converted_path)
        diagnostics["remote_attempts"].append({
            "url": remote_cif_url,
            "format": "cif",
            "converted_pdb": converted_path,
            "usable": validation["usable"]
        })
        if validation["usable"]:
            diagnostics["selected_source"] = remote_cif_url
            diagnostics["selected_source_format"] = "cif"
            diagnostics["converted_to_pdb"] = converted_path
            diagnostics["served_pdb_path"] = converted_path
            diagnostics["contains_ligand"] = validation["contains_ligand"]
            diagnostics["requested_ligand_hetatm_count"] = validation["ligand_atom_count"]
            diagnostics["has_protein_atoms"] = validation["has_atom_records"]
            diagnostics["hetatm_summary"] = validation["hetatm_summary"]
            return diagnostics
    except Exception as exc:
        diagnostics["remote_attempts"].append({"url": remote_cif_url, "format": "cif", "usable": False, "error": str(exc)})

    return diagnostics

@app.route('/get_virus_names_list_distinct', methods=['GET'])
def get_virus_names_list_distinct():
    mode = _normalized_backend_mode()
    if mode == "randy":
        try:
            return jsonify(randy_get('virus-proteins/virus-names'))
        except RandyBackendError as exc:
            return jsonify({"error": str(exc)}), exc.status_code
    if mode == "auto" and randy_available():
        try:
            return jsonify(randy_get('virus-proteins/virus-names'))
        except RandyBackendError:
            logging.warning("Falling back to local virus names list")
    try:
        with _connect_local_db(QUERY_PROTEIN_REQUIRED_TABLES) as conn:
            cursor = conn.cursor()
            cursor.execute("SELECT DISTINCT virus_name FROM Virus_Proteins")
            virus_names = [row[0] for row in cursor.fetchall()]
        return jsonify(virus_names)
    except RandyBackendError as exc:
        return jsonify({"error": str(exc)}), 500

@app.route('/get_protein_types_list_distinct', methods=['GET'])
def get_protein_types_list_distinct():
    mode = _normalized_backend_mode()
    if mode == "randy":
        try:
            return jsonify(randy_get('virus-proteins/protein-types'))
        except RandyBackendError as exc:
            return jsonify({"error": str(exc)}), exc.status_code
    if mode == "auto" and randy_available():
        try:
            return jsonify(randy_get('virus-proteins/protein-types'))
        except RandyBackendError:
            logging.warning("Falling back to local protein types list")
    try:
        with _connect_local_db(QUERY_PROTEIN_REQUIRED_TABLES) as conn:
            cursor = conn.cursor()
            cursor.execute("SELECT DISTINCT protein FROM Virus_Proteins")
            protein_types = [row[0] for row in cursor.fetchall()]
        return jsonify(protein_types)
    except RandyBackendError as exc:
        return jsonify({"error": str(exc)}), 500



@app.route('/get_protein_query_filter_options', methods=['GET'])
def get_protein_query_filter_options_route():
    """Provide the dependent options for the Protein Query filter cascade."""
    virus_names = request.args.getlist('virus_name')
    protein_types = request.args.getlist('protein_type')
    mode = _normalized_backend_mode()
    remote_params = [
        *(('virus_name', virus_name) for virus_name in virus_names),
        *(('protein_type', protein_type) for protein_type in protein_types),
    ]
    if mode == "randy":
        try:
            return jsonify(randy_get('virus-proteins/filter-options', params=remote_params))
        except RandyBackendError as exc:
            return jsonify({"error": str(exc)}), exc.status_code
    if mode == "auto" and randy_available():
        try:
            return jsonify(randy_get('virus-proteins/filter-options', params=remote_params))
        except RandyBackendError:
            logging.warning("Falling back to local Protein Query filter options")
    try:
        with _connect_local_db(QUERY_PROTEIN_REQUIRED_TABLES) as conn:
            return jsonify(get_protein_query_filter_options(
                conn,
                virus_names=virus_names,
                protein_types=protein_types,
            ))
    except RandyBackendError as exc:
        return jsonify({"error": str(exc)}), 500




@app.route('/get_pdbs_for_virus_protein', methods=['POST'])
def get_pdbs_for_virus_protein():
    data = request.get_json(silent=True) or {}
    mode = _normalized_backend_mode()
    if mode == "randy":
        try:
            return jsonify(randy_post('virus-proteins/pdbs', json=data))
        except RandyBackendError as exc:
            return jsonify({"error": str(exc)}), exc.status_code
    if mode == "auto" and randy_available():
        try:
            return jsonify(randy_post('virus-proteins/pdbs', json=data))
        except RandyBackendError:
            logging.warning("Falling back to local pdb lookup for virus/protein query")

    virus_name = str(data.get('virus_name') or '').strip()
    protein_types = data.get('protein_types') or []
    ligand_filter = data.get('ligand', None)
    if not virus_name or not protein_types:
        return jsonify({'error': 'virus_name and protein_types are required.'}), 400
    try:
        with _connect_local_db(QUERY_PROTEIN_REQUIRED_TABLES) as conn:
            pdb_codes = get_exportable_protein_query_pdbs(
                conn, virus_name, protein_types, ligand_filter
            )
        return jsonify({'pdb_codes': pdb_codes})
    except RandyBackendError as exc:
        return jsonify({"error": str(exc)}), 500





@app.route('/export_data_to_excel', methods=['POST'])
def export_data_to_excel():
    data = request.get_json(silent=True) or {}
    requested_pdb_codes = []
    for pdb_code in data.get('pdb_codes', []):
        normalized = str(pdb_code or '').strip().upper()
        if normalized and normalized not in requested_pdb_codes:
            requested_pdb_codes.append(normalized)
    pdb_codes = requested_pdb_codes
    data_sets = data.get('data_sets', [])

    if not pdb_codes or not data_sets:
        return jsonify({'success': False, 'error': 'No PDB codes or datasets provided.'}), 400

    mode = _normalized_backend_mode()
    if mode == "randy":
        try:
            page_metadata = _remote_page_metadata()
        except RandyBackendError as exc:
            return jsonify({'success': False, 'error': str(exc)}), exc.status_code
    elif mode == "auto" and randy_available():
        try:
            page_metadata = _remote_page_metadata()
        except RandyBackendError:
            logging.warning("Falling back to local page metadata for export validation")
            page_metadata = _local_page_metadata()
    else:
        page_metadata = _local_page_metadata()

    allowed_data_sets = set(page_metadata.get("available_export_data_sets", []))
    invalid_data_sets = [data_set for data_set in data_sets if data_set not in allowed_data_sets]
    if invalid_data_sets:
        return jsonify({'success': False, 'error': f'Unsupported data sets requested: {invalid_data_sets}'}), 400
    try:
        remote_payload = None
        if mode == "randy":
            remote_payload = randy_post('export-data', json=data)
        elif mode == "auto" and randy_available():
            try:
                remote_payload = randy_post('export-data', json=data)
            except RandyBackendError:
                logging.warning("Falling back to local export data generation")

        non_exportable_pdb_codes = []
        if remote_payload is None:
            eligibility_tables = (
                "structures",
                "ligand_instances",
                "ligand_instance_atoms",
            )
            with _connect_local_db(eligibility_tables) as conn:
                eligible_pdb_codes = _protein_query_export_eligible_pdbs(
                    conn, pdb_codes=pdb_codes
                )
        else:
            # The remote API returns the datasets themselves, so the same
            # no-data guard below remains authoritative for RANDY responses.
            eligible_pdb_codes = list(pdb_codes)

        eligible_code_set = set(eligible_pdb_codes)
        non_exportable_pdb_codes = [
            code for code in pdb_codes if code not in eligible_code_set
        ]
        if not eligible_pdb_codes:
            return jsonify({
                'success': False,
                'error': 'no_exportable_data',
                'message': 'No retained V-LiSEMOD ligand data are available for the selected structure(s).',
                'requested_pdb_codes': pdb_codes,
                'non_exportable_pdb_codes': non_exportable_pdb_codes,
            }), 422

        populated_data_sets = {}
        dataset_row_counts = {}
        for data_set in data_sets:
            try:
                if remote_payload is not None:
                    dataset_rows = (remote_payload.get("data_sets") or {}).get(data_set, [])
                    df = pd.DataFrame(dataset_rows)
                else:
                    with _connect_local_db(_required_tables_for_datasets([data_set])) as conn:
                        placeholders = ', '.join(['?'] * len(eligible_pdb_codes))
                        query = data_set_queries[data_set].format(placeholders=placeholders)
                        df = pd.read_sql(query, conn, params=eligible_pdb_codes)
                dataset_row_counts[data_set] = int(len(df.index))
                if not df.empty:
                    populated_data_sets[data_set] = df
            except RandyBackendError as exc:
                return jsonify({'success': False, 'error': str(exc)}), exc.status_code

        if not populated_data_sets:
            return jsonify({
                'success': False,
                'error': 'no_exportable_data',
                'message': 'No rows are available for the selected data set(s) and retained structure(s).',
                'requested_pdb_codes': pdb_codes,
                'exportable_pdb_codes': eligible_pdb_codes,
                'non_exportable_pdb_codes': non_exportable_pdb_codes,
            }), 422

        excel_buffer = BytesIO()
        with pd.ExcelWriter(excel_buffer, engine='xlsxwriter') as writer:
            for data_set, df in populated_data_sets.items():
                df.to_excel(writer, sheet_name=data_set[:30], index=False)

        export_summary = {
            'requested_pdb_codes': pdb_codes,
            'exportable_pdb_codes': eligible_pdb_codes,
            'non_exportable_pdb_codes': non_exportable_pdb_codes,
            'dataset_row_counts': dataset_row_counts,
        }
        zip_buffer = BytesIO()
        with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zipf:
            for data_set, df in populated_data_sets.items():
                csv_buffer = io.StringIO()
                df.to_csv(csv_buffer, index=False)
                zipf.writestr(f"{data_set.replace(' ', '_')}.csv", csv_buffer.getvalue())
            zipf.writestr('export_summary.json', json.dumps(export_summary, indent=2))
            zipf.writestr('combined_data.xlsx', excel_buffer.getvalue())

        zip_buffer.seek(0)  # Move back to the start of the BytesIO buffer

        # Send the ZIP archive as a downloadable file
        return send_file(zip_buffer, as_attachment=True, download_name='data_sets.zip', mimetype='application/zip')

    except Exception as e:
        logging.error("Failed during file creation: %s", e)
        return jsonify({'success': False, 'error': str(e)}), 500



# @app.route('/export_data_to_excel', methods=['POST'])
# def export_data_to_excel():
#     data = request.json
#     pdb_codes = data.get('pdb_codes', [])
#     data_sets = data.get('data_sets', [])
#     output_dir = 'output_files'  # Directory where CSV files will be saved

#     # Ensure the output directory exists
#     os.makedirs(output_dir, exist_ok=True)

#     if not pdb_codes or not data_sets:
#         return jsonify({'success': False, 'error': 'No PDB codes or datasets provided.'}), 400

#     conn = connect_db()
#     response = {}
#     try:
#         for data_set in data_sets:
#             placeholders = ', '.join(['?'] * len(pdb_codes))
#             query = data_set_queries[data_set].format(placeholders=placeholders)
#             try:
#                 df = pd.read_sql(query, conn, params=pdb_codes)
#                 if not df.empty:
#                     csv_file_path = os.path.join(output_dir, f"{data_set.replace(' ', '_')}.csv")
#                     df.to_csv(csv_file_path, index=False)
#                     response[data_set] = f"File saved: {csv_file_path}"
#                 else:
#                     response[data_set] = "No data to save."
#             except Exception as e:
#                 response[data_set] = f"Error querying database: {str(e)}"
#                 logging.error("Error querying database: %s", e)

#         # Instead of zipping, send the folder as a downloadable
#         archive_path = shutil.make_archive(output_dir, 'zip', output_dir)
#         return send_file(archive_path, as_attachment=True, download_name='data_sets.zip')
#     except Exception as e:
#         logging.error("Failed during file creation: %s", e)
#         return jsonify({'success': False, 'error': str(e)}), 500
#     finally:
#         conn.close()
    

@app.route('/query_protein_virus_page')
def query_protein_virus_page():
    mode = _normalized_backend_mode()
    if mode == "randy":
        try:
            metadata = _remote_page_metadata()
        except RandyBackendError as exc:
            return render_template('query_protein_virus.html', available_export_data_sets=[], protacability_data_available=False, backend_error=str(exc)), exc.status_code
    elif mode == "auto" and randy_available():
        try:
            metadata = _remote_page_metadata()
        except RandyBackendError:
            logging.warning("Falling back to local page metadata for query_protein_virus_page")
            metadata = _local_page_metadata()
    else:
        metadata = _local_page_metadata()
    return render_template(
        'query_protein_virus.html',
        available_export_data_sets=metadata.get("available_export_data_sets", []),
        protacability_data_available=metadata.get("protacability_data_available", False)
    )


############################
###Scripts for Ligand Querying



@app.route('/get_ligands_with_synonyms', methods=['GET'])
def get_ligands_with_synonyms():
    mode = _normalized_backend_mode()
    if mode == "randy":
        try:
            return jsonify(randy_get('ligands/with-synonyms'))
        except RandyBackendError as exc:
            return jsonify({"error": str(exc)}), exc.status_code
    if mode == "auto" and randy_available():
        try:
            return jsonify(randy_get('ligands/with-synonyms'))
        except RandyBackendError:
            logging.warning("Falling back to local ligands with synonyms")
    try:
        return jsonify(_local_get_ligands_with_synonyms_payload())
    except RandyBackendError as exc:
        return jsonify({"error": str(exc)}), 500





@app.route("/get_ligand_info/<ligand_code>", methods=["GET"])
def get_ligand_info(ligand_code):
    mode = _normalized_backend_mode()
    if mode == "randy":
        try:
            payload = randy_get('ligand-info', params={"ligand_code": ligand_code})
        except RandyBackendError as exc:
            return jsonify({"error": str(exc)}), exc.status_code
    elif mode == "auto" and randy_available():
        try:
            payload = randy_get('ligand-info', params={"ligand_code": ligand_code})
        except RandyBackendError:
            logging.warning("Falling back to local ligand info for %s", ligand_code)
            payload = _local_get_ligand_info_payload(ligand_code)
    else:
        payload = _local_get_ligand_info_payload(ligand_code)

    if payload:
        return jsonify(payload)
    return jsonify({"error": "Ligand not found"}), 404









####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################






PROTAC_BUILDER_EXTERNAL_URL = os.environ.get("PROTAC_BUILDER_EXTERNAL_URL", "https://protacbuilder.com")


@app.route("/copy", defaults={"legacy_path": ""}, methods=["GET", "HEAD", "OPTIONS"])
@app.route("/copy/<path:legacy_path>", methods=["GET", "HEAD", "OPTIONS"])
def redirect_legacy_copy_routes(legacy_path):
    target = PROTAC_BUILDER_EXTERNAL_URL.rstrip("/")
    if legacy_path:
        target = f"{target}/{legacy_path.lstrip('/')}"
    query_string = request.query_string.decode("utf-8")
    if query_string:
        target = f"{target}?{query_string}"
    return redirect(target, code=302)

@app.route("/run-drug-analysis", methods=["POST", "GET", "OPTIONS"])
def run_drug_analysis():
    return jsonify({
        "success": False,
        "error": "Deprecated on VLISEMOD. Use standalone PROTAC Builder.",
        "protac_builder_url": PROTAC_BUILDER_EXTERNAL_URL,
        "recommended_endpoint": "/api/deeppk/run"
    }), 410

@app.route("/download/<filename>")
def download_file(filename):
    return jsonify({
        "success": False,
        "error": "Deprecated on VLISEMOD. Use standalone PROTAC Builder.",
        "protac_builder_url": PROTAC_BUILDER_EXTERNAL_URL
    }), 410


@app.route('/static/python/PrepFiles.py')
def serve_prep_script():
    return jsonify({
        "success": False,
        "error": "Deprecated on VLISEMOD. Use standalone PROTAC Builder.",
        "protac_builder_url": PROTAC_BUILDER_EXTERNAL_URL
    }), 410


@app.route('/process', methods=['POST', 'GET', 'OPTIONS'])
def process_files():
    return jsonify({
        "success": False,
        "error": "Deprecated on VLISEMOD. Use standalone PROTAC Builder.",
        "protac_builder_url": PROTAC_BUILDER_EXTERNAL_URL
    }), 410


@app.route('/ligase_details', methods=['GET', 'POST', 'OPTIONS'])
def ligase_details():
    return jsonify({
        "success": False,
        "error": "Deprecated on VLISEMOD. Use standalone PROTAC Builder.",
        "protac_builder_url": PROTAC_BUILDER_EXTERNAL_URL
    }), 410


if not ENABLE_DRUG_GPT:
    @app.route("/drugapp/", methods=["GET", "HEAD", "OPTIONS"])
    def drug_gpt_disabled_home():
        if request.method == "HEAD":
            return ("", 503)
        return render_template_string(
            """
            {% extends "base.html" %}
            {% block title %}V-LiSEMOD | Drug GPT Disabled{% endblock %}
            {% block content %}
            <section class="page-intro">
                <span class="page-kicker">Drug GPT</span>
                <h1>Drug GPT is currently disabled for this deployment.</h1>
                <p>V-LiSEMOD is running in its lightweight default mode, so no local LLM weights are loaded automatically.</p>
            </section>
            <section class="page-card">
                <p>Set <code>ENABLE_DRUG_GPT=1</code> to re-enable the Drug GPT route and <code>ENABLE_LOCAL_LLM=1</code> to allow local model loading.</p>
            </section>
            {% endblock %}
            """
        ), 503


    @app.route("/drugapp/query", methods=["POST", "GET", "OPTIONS"])
    def drug_gpt_disabled_query():
        return jsonify({
            "error": "Drug GPT is currently disabled for this deployment.",
            "enable_drug_gpt": app.config.get("ENABLE_DRUG_GPT", False),
            "enable_local_llm": app.config.get("ENABLE_LOCAL_LLM", False),
        }), 503


@app.route('/protacability_page', methods=['GET', 'POST', 'OPTIONS'])
def protacability_page():
    mode = _normalized_backend_mode()
    if mode == "randy":
        try:
            metadata = _remote_page_metadata()
        except RandyBackendError as exc:
            return render_template('protacability_assessment.html', protacability_data_available=False, backend_error=str(exc)), exc.status_code
    elif mode == "auto" and randy_available():
        try:
            metadata = _remote_page_metadata()
        except RandyBackendError:
            logging.warning("Falling back to local page metadata for protacability_page")
            metadata = _local_page_metadata()
    else:
        metadata = _local_page_metadata()
    return render_template(
        'protacability_assessment.html',
        protacability_data_available=metadata.get("protacability_data_available", False)
    )


@app.route('/api/protacability/filters')
def protacability_filters():
    mode = _normalized_backend_mode()
    if mode == "randy":
        try:
            return jsonify(_remote_protacability_get("protacability/filter-options", params=request.args, max_bytes=2 * 1024 * 1024))
        except RandyBackendError as exc:
            return jsonify({"data_available": False, "message": str(exc)}), exc.status_code
    if mode == "auto" and randy_available():
        try:
            return jsonify(_remote_protacability_get("protacability/filter-options", params=request.args, max_bytes=2 * 1024 * 1024))
        except RandyBackendError:
            logging.warning("Falling back to local PROTACability filters payload")
    try:
        payload = _load_protacability_source_payload()
    except RandyBackendError as exc:
        return jsonify({"data_available": False, "message": str(exc)}), exc.status_code

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
            attachment_rows=payload.get("attachment_rows", []),
        )
    )


@app.route('/api/protacability/filter_options')
def protacability_filter_options():
    mode = _normalized_backend_mode()
    started_at = time.perf_counter()
    # The normal browser must be served by Randy's aggregated endpoint.  Do
    # not special-case target mode into a local/raw source-data path.
    if mode == "randy":
        try:
            payload = _remote_protacability_get("protacability/filter-options", params=request.args, max_bytes=2 * 1024 * 1024)
            logging.info(
                "PROTAC filter_options backend=randy elapsed_ms=%.1f payload_bytes=%s ligand_count=%s protein_count=%s virus_count=%s",
                (time.perf_counter() - started_at) * 1000,
                len(json.dumps(payload)),
                len(payload.get("ligands") or []),
                len(payload.get("protein_types") or []),
                len(payload.get("virus_names") or []),
            )
            return jsonify(payload)
        except RandyBackendError as exc:
            return jsonify({"data_available": False, "message": str(exc)}), exc.status_code
    if mode == "auto" and randy_available():
        try:
            payload = _remote_protacability_get("protacability/filter-options", params=request.args, max_bytes=2 * 1024 * 1024)
            logging.info(
                "PROTAC filter_options backend=auto-randy elapsed_ms=%.1f payload_bytes=%s ligand_count=%s protein_count=%s virus_count=%s",
                (time.perf_counter() - started_at) * 1000,
                len(json.dumps(payload)),
                len(payload.get("ligands") or []),
                len(payload.get("protein_types") or []),
                len(payload.get("virus_names") or []),
            )
            return jsonify(payload)
        except RandyBackendError:
            logging.warning("Falling back to local PROTACability filter_options payload")
    try:
        payload = _load_protacability_source_payload()
    except RandyBackendError as exc:
        return jsonify({"data_available": False, "message": str(exc)}), exc.status_code

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
    response_payload = _build_protacability_filter_options_payload_from_rows(
        payload.get("assessment_rows", []),
        payload.get("readiness_rows", []),
        payload.get("warhead_rows", []),
        request.args,
        attachment_rows=payload.get("attachment_rows", []),
    )
    logging.info(
        "PROTAC filter_options backend=local elapsed_ms=%.1f payload_bytes=%s ligand_count=%s protein_count=%s virus_count=%s",
        (time.perf_counter() - started_at) * 1000,
        len(json.dumps(response_payload)),
        len(response_payload.get("ligands") or []),
        len(response_payload.get("protein_types") or []),
        len(response_payload.get("virus_names") or []),
    )
    return jsonify(response_payload)


@app.route('/api/protacability/search')
def protacability_search():
    mode = _normalized_backend_mode()
    started_at = time.perf_counter()
    if mode == "randy":
        try:
            payload = _remote_protacability_get("protacability/search", params=request.args)
            logging.info(
                "PROTAC search backend=randy view=%s offset=%s limit=%s total_rows=%s mapped_exposed=%s elapsed_ms=%.1f payload_bytes=%s",
                payload.get("view"),
                payload.get("offset"),
                payload.get("limit"),
                payload.get("total_rows"),
                (payload.get("summary") or {}).get("candidate_warheads_with_exposed_mapped_atoms"),
                (time.perf_counter() - started_at) * 1000,
                len(json.dumps(payload)),
            )
            return jsonify(payload)
        except RandyBackendError as exc:
            return jsonify({"data_available": False, "message": str(exc), "rows": [], "summary": {}}), exc.status_code
    if mode == "auto" and randy_available():
        try:
            payload = _remote_protacability_get("protacability/search", params=request.args)
            logging.info(
                "PROTAC search backend=auto-randy view=%s offset=%s limit=%s total_rows=%s mapped_exposed=%s elapsed_ms=%.1f payload_bytes=%s",
                payload.get("view"),
                payload.get("offset"),
                payload.get("limit"),
                payload.get("total_rows"),
                (payload.get("summary") or {}).get("candidate_warheads_with_exposed_mapped_atoms"),
                (time.perf_counter() - started_at) * 1000,
                len(json.dumps(payload)),
            )
            return jsonify(payload)
        except RandyBackendError:
            logging.warning("Falling back to local PROTACability search payload")
    filters = _build_protacability_filters(request.args)
    try:
        source_payload = _load_protacability_source_payload(
            ligand=filters.get("ligand") or None,
            ligand_instance_id=filters.get("ligand_instance_id") or None,
            request_args=request.args,
        )
    except RandyBackendError as exc:
        return jsonify({"data_available": False, "message": str(exc), "rows": [], "summary": {}}), exc.status_code

    if not source_payload.get("data_available"):
        return jsonify({
            "data_available": False,
            "message": "PROTACability data has not been imported yet. Run tools/import_protacability_data.py after generating the expansion outputs.",
            "rows": [],
            "summary": {}
        })

    if source_payload.get("deep_link_error"):
        return jsonify({
            "data_available": False,
            "message": source_payload["deep_link_error"],
            "rows": [],
            "summary": {},
        }), 400

    payload = _prepare_protacability_result_set_from_rows(
        source_payload.get("assessment_rows", []),
        source_payload.get("readiness_rows", []),
        source_payload.get("warhead_rows", []),
        request.args,
        attachment_rows=source_payload.get("attachment_rows", []),
    )
    response_payload = {
        "data_available": True,
        "view": payload["view"],
        "collapse_labels": payload["collapse_labels"],
        "rows": payload["rows"],
        "summary": payload["summary"],
        "limit": payload["limit"],
        "offset": payload["offset"],
        "total_rows": payload["total_rows"],
        "has_more": payload["has_more"],
        "sort": payload["sort"],
        "deep_link_context": source_payload.get("deep_link_context"),
    }
    logging.info(
        "PROTAC search backend=local view=%s offset=%s limit=%s total_rows=%s mapped_exposed=%s elapsed_ms=%.1f payload_bytes=%s",
        response_payload.get("view"),
        response_payload.get("offset"),
        response_payload.get("limit"),
        response_payload.get("total_rows"),
        (response_payload.get("summary") or {}).get("candidate_warheads_with_exposed_mapped_atoms"),
        (time.perf_counter() - started_at) * 1000,
        len(json.dumps(response_payload)),
    )
    return jsonify(response_payload)


@app.route('/api/protacability/detail/<pdb_code>/<chain_id>')
def protacability_detail(pdb_code, chain_id):
    mode = _normalized_backend_mode()
    if mode == "randy":
        try:
            return jsonify(_normalize_remote_attachment_site_display(_remote_protacability_get(f"protacability/detail/{pdb_code}/{chain_id}", max_bytes=2 * 1024 * 1024)))
        except RandyBackendError as exc:
            return jsonify({"data_available": False, "message": str(exc)}), exc.status_code
    if mode == "auto" and randy_available():
        try:
            return jsonify(_normalize_remote_attachment_site_display(_remote_protacability_get(f"protacability/detail/{pdb_code}/{chain_id}", max_bytes=2 * 1024 * 1024)))
        except RandyBackendError:
            logging.warning("Falling back to local PROTACability detail payload for %s/%s", pdb_code, chain_id)
    try:
        source_payload = _load_protacability_source_payload(pdb_code=pdb_code, include_lysine=True, include_inventory=True)
    except RandyBackendError as exc:
        return jsonify({"data_available": False, "message": str(exc)}), exc.status_code

    if not source_payload.get("data_available"):
        return jsonify({
            "data_available": False,
            "message": "PROTACability data has not been imported yet. Run tools/import_protacability_data.py after generating the expansion outputs."
        }), 404

    readiness_rows = source_payload.get("readiness_rows", [])
    warhead_rows = source_payload.get("warhead_rows", [])
    raw_assessment_rows = [
        row for row in source_payload.get("assessment_rows", [])
        if row.get("pdb_code") == pdb_code and row.get("chain_id") == chain_id
    ]

    assessment_rows = _decorate_protacability_rows(
        raw_assessment_rows,
        collapse_labels=True,
        readiness_rows=readiness_rows,
        warhead_rows=warhead_rows,
        attachment_rows=source_payload.get("attachment_rows", []),
    )
    assessment = max(assessment_rows, key=_row_priority_key) if assessment_rows else None

    if assessment is None:
        return jsonify({"error": "Assessment row not found"}), 404

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
        for row in source_payload.get("lysine_rows", [])
        if row.get("pdb_code") == pdb_code and row.get("chain_id") == chain_id
    ]

    ligand_inventory = [
        {
            "ligand_instance_id": row.get("ligand_instance_id"),
            "model_id": row.get("model_id"),
            "ligand_resname": row.get("ligand_resname"),
            "ligand_chain": row.get("ligand_chain"),
            "ligand_residue_id": row.get("ligand_residue_id"),
            "ligand_insertion_code": row.get("ligand_insertion_code"),
            "ligand_atom_count": row.get("ligand_atom_count"),
            "ligand_heavy_atom_count": row.get("ligand_heavy_atom_count"),
            "centroid_x": row.get("centroid_x"),
            "centroid_y": row.get("centroid_y"),
            "centroid_z": row.get("centroid_z"),
        }
        for row in source_payload.get("ligand_inventory", [])
        if row.get("pdb_code") == pdb_code
    ]
    ligand_inventory = _attach_ligand_label_asym_ids(ligand_inventory)

    related_chains = [
        {
            "chain_id": row.get("chain_id"),
            "protacability_proxy_score": row.get("protacability_proxy_score"),
            "protacability_tier": row.get("protacability_tier"),
            "candidate_ligand_count": row.get("candidate_ligand_count"),
            "exposed_lys_count": row.get("exposed_lys_count"),
        }
        for row in source_payload.get("assessment_rows", [])
        if row.get("pdb_code") == pdb_code
    ]
    with connect_db_row() as conn:
        attachment_sites = _attachment_detail_payload(conn, assessment)

    return jsonify({
        "data_available": True,
        "assessment": dict(assessment),
        "lysine_rows": lysine_rows,
        "ligand_inventory": ligand_inventory,
        "ligand_contexts": _serialize_ligand_contexts(ligand_inventory),
        "related_chains": related_chains,
        "attachment_sites": attachment_sites,
    })


@app.route('/api/protacability/structure_detail/<pdb_code>')
def protacability_structure_detail(pdb_code):
    mode = _normalized_backend_mode()
    if mode == "randy":
        try:
            return jsonify(_normalize_remote_attachment_site_display(_remote_protacability_get(f"protacability/structure-detail/{pdb_code}", params=request.args, max_bytes=2 * 1024 * 1024)))
        except RandyBackendError as exc:
            return jsonify({"data_available": False, "message": str(exc)}), exc.status_code
    if mode == "auto" and randy_available():
        try:
            return jsonify(_normalize_remote_attachment_site_display(_remote_protacability_get(f"protacability/structure-detail/{pdb_code}", params=request.args, max_bytes=2 * 1024 * 1024)))
        except RandyBackendError:
            logging.warning("Falling back to local PROTACability structure detail payload for %s", pdb_code)
    try:
        source_payload = _load_protacability_source_payload(pdb_code=pdb_code, include_lysine=True, include_inventory=True)
    except RandyBackendError as exc:
        return jsonify({"data_available": False, "message": str(exc)}), exc.status_code

    if not source_payload.get("data_available"):
        return jsonify({
            "data_available": False,
            "message": "PROTACability data has not been imported yet. Run tools/import_protacability_data.py after generating the expansion outputs."
        }), 404

    collapse_labels = _protacability_collapse_labels(request.args.get("collapse_labels"))
    virus_name = (request.args.get("virus_name") or "").strip()
    protein_type = (request.args.get("protein_type") or "").strip()
    requested_ligand_instance_id = str(request.args.get("ligand_instance_id") or "").strip()

    readiness_rows = source_payload.get("readiness_rows", [])
    warhead_rows = source_payload.get("warhead_rows", [])
    decorated_rows = _decorate_protacability_rows(
        source_payload.get("assessment_rows", []),
        collapse_labels=collapse_labels,
        readiness_rows=readiness_rows,
        warhead_rows=warhead_rows,
        attachment_rows=source_payload.get("attachment_rows", []),
    )
    if virus_name:
        decorated_rows = [row for row in decorated_rows if row.get("virus_name") == virus_name]
    if protein_type:
        decorated_rows = [row for row in decorated_rows if (row.get("display_protein_type") if collapse_labels else row.get("protein_type")) == protein_type]

    if not decorated_rows:
        return jsonify({"error": "Structure summary not found"}), 404

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
        for row in source_payload.get("lysine_rows", [])
        if row.get("pdb_code") == pdb_code and row.get("chain_id") == representative_chain
    ]

    ligand_inventory = [
        {
            "ligand_instance_id": row.get("ligand_instance_id"),
            "model_id": row.get("model_id"),
            "ligand_resname": row.get("ligand_resname"),
            "ligand_chain": row.get("ligand_chain"),
            "ligand_residue_id": row.get("ligand_residue_id"),
            "ligand_insertion_code": row.get("ligand_insertion_code"),
            "ligand_atom_count": row.get("ligand_atom_count"),
            "ligand_heavy_atom_count": row.get("ligand_heavy_atom_count"),
            "centroid_x": row.get("centroid_x"),
            "centroid_y": row.get("centroid_y"),
            "centroid_z": row.get("centroid_z"),
        }
        for row in source_payload.get("ligand_inventory", [])
        if row.get("pdb_code") == pdb_code
    ]
    ligand_inventory = _attach_ligand_label_asym_ids(ligand_inventory)
    selected_ligand_instance = None
    if requested_ligand_instance_id:
        selected_ligand_instance = next(
            (
                record for record in ligand_inventory
                if str(record.get("ligand_instance_id") or "") == requested_ligand_instance_id
            ),
            None,
        )
        if selected_ligand_instance is None:
            return jsonify({"error": "Ligand occurrence was not found for this structure"}), 404
    preferred_ligands = _split_candidate_ligands(summary_row.get("candidate_ligand_resnames_full"))
    representative_ligand = selected_ligand_instance or _pick_representative_ligand_record(
        ligand_inventory,
        preferred_ligands=preferred_ligands,
        allow_glycan=summary_row.get("ligand_context_class") == "glycan_only",
        preferred_chain=representative_chain
    )
    attachment_lookup_row = {
        **summary_row,
        **(representative_ligand or {}),
        "pdb_code": pdb_code,
        "model_id": (representative_ligand or {}).get("model_id") or summary_row.get("model_id") or 0,
    }
    with connect_db_row() as conn:
        attachment_sites = _attachment_detail_payload(conn, attachment_lookup_row)
    return jsonify({
        "data_available": True,
        "summary": summary_row,
        "chain_rows": chain_rows,
        "representative_chain_id": representative_chain,
        "representative_ligand": representative_ligand,
        "selected_ligand_instance": selected_ligand_instance,
        "representative_ligand_resname": (representative_ligand or {}).get("ligand_resname"),
        "representative_ligand_chain": (representative_ligand or {}).get("ligand_chain"),
        "representative_ligand_residue_id": (representative_ligand or {}).get("ligand_residue_id"),
        "lysine_rows": lysine_rows,
        "ligand_inventory": ligand_inventory,
        "ligand_contexts": _serialize_ligand_contexts(ligand_inventory),
        "attachment_sites": attachment_sites,
    })


@app.route('/api/protacability/protein_detail')
def protacability_protein_detail():
    virus_name = (request.args.get("virus_name") or "").strip()
    protein_type = (request.args.get("protein_type") or "").strip()
    mode = _normalized_backend_mode()
    if mode == "randy":
        try:
            return jsonify(_normalize_remote_attachment_site_display(_remote_protacability_get("protacability/protein-detail", params=request.args, max_bytes=2 * 1024 * 1024)))
        except RandyBackendError as exc:
            return jsonify({"data_available": False, "message": str(exc)}), exc.status_code
    if mode == "auto" and randy_available():
        try:
            return jsonify(_normalize_remote_attachment_site_display(_remote_protacability_get("protacability/protein-detail", params=request.args, max_bytes=2 * 1024 * 1024)))
        except RandyBackendError:
            logging.warning("Falling back to local PROTACability protein detail payload for %s / %s", virus_name, protein_type)
    try:
        source_payload = _load_protacability_source_payload(virus_name=virus_name, protein_type=protein_type)
    except RandyBackendError as exc:
        return jsonify({"data_available": False, "message": str(exc)}), exc.status_code

    if not source_payload.get("data_available"):
        return jsonify({
            "data_available": False,
            "message": "PROTACability data has not been imported yet. Run tools/import_protacability_data.py after generating the expansion outputs."
        }), 404

    collapse_labels = _protacability_collapse_labels(request.args.get("collapse_labels"))
    readiness_rows = source_payload.get("readiness_rows", [])
    warhead_rows = source_payload.get("warhead_rows", [])
    decorated_rows = _decorate_protacability_rows(
        source_payload.get("assessment_rows", []),
        collapse_labels=collapse_labels,
        readiness_rows=readiness_rows,
        warhead_rows=warhead_rows,
        attachment_rows=source_payload.get("attachment_rows", []),
    )
    decorated_rows = [
        row for row in decorated_rows
        if row.get("virus_name") == virus_name
        and ((row.get("display_protein_type") if collapse_labels else row.get("protein_type")) == protein_type)
    ]

    if not decorated_rows:
        return jsonify({"error": "Protein summary not found"}), 404

    protein_rows = _group_protein_rows(decorated_rows)
    protein_row = protein_rows[0]
    structure_rows = _group_structure_rows(decorated_rows)
    structure_rows, _ = _sort_protacability_rows(structure_rows, "summary", "best_score_desc")
    tier_distribution = dict(Counter(row.get("best_tier") or "Unknown" for row in structure_rows))
    attachment_lookup_row = {
        **protein_row,
        "pdb_code": protein_row.get("top_pdb_code"),
        "ligand_resname": protein_row.get("best_ligand_resname"),
        "ligand_chain": protein_row.get("best_ligand_chain"),
        "ligand_residue_id": protein_row.get("best_ligand_residue_id"),
    }
    with connect_db_row() as conn:
        attachment_sites = _attachment_detail_payload(conn, attachment_lookup_row)
    return jsonify({
        "data_available": True,
        "summary": protein_row,
        "top_structures": structure_rows[:10],
        "tier_distribution": tier_distribution,
        "explanation": "This view groups multiple structures and chains into a single protein-level summary so repeated biological contexts do not dominate the table.",
        "attachment_sites": attachment_sites,
    })


@app.route('/api/protacability/target_detail')
def protacability_target_detail():
    collapse_labels = _protacability_collapse_labels(request.args.get("collapse_labels"))
    virus_name = (request.args.get("virus_name") or "").strip()
    protein_type = (request.args.get("protein_type") or "").strip()
    canonical_target_id = (request.args.get("canonical_target_id") or "").strip()
    ligand_context_class = (request.args.get("ligand_context_class") or "").strip()
    min_score = request.args.get("min_score", type=float)

    if not virus_name or not (canonical_target_id or protein_type):
        return jsonify({"error": "virus_name and canonical_target_id or protein_type are required"}), 400

    mode = _normalized_backend_mode()
    # Randy owns the canonical occurrence authority and accepts both the
    # stable canonical_target_id and the legacy display parameter.  Heroku
    # must not attempt a local database read for canonical detail views.
    if mode == "randy":
        try:
            return jsonify(_normalize_remote_attachment_site_display(_remote_protacability_get("protacability/target-detail", params=request.args, max_bytes=10 * 1024 * 1024)))
        except RandyBackendError as exc:
            return jsonify({"data_available": False, "message": str(exc)}), exc.status_code
    if mode == "auto" and randy_available():
        try:
            return jsonify(_normalize_remote_attachment_site_display(_remote_protacability_get("protacability/target-detail", params=request.args, max_bytes=10 * 1024 * 1024)))
        except RandyBackendError:
            logging.warning("Falling back to local PROTACability target detail payload for %s / %s", virus_name, protein_type)

    if canonical_target_id:
        try:
            conn = connect_db_row()
            try:
                assessment_rows = _load_canonical_target_browser_assessment_rows(conn)
                readiness_rows, warhead_rows, attachment_rows = _load_protacability_enrichment_tables(conn)
            finally:
                conn.close()
        except sqlite3.Error as exc:
            return jsonify({"data_available": False, "message": str(exc)}), 500
        rows = _decorate_protacability_rows(
            assessment_rows,
            collapse_labels=collapse_labels,
            readiness_rows=readiness_rows,
            warhead_rows=warhead_rows,
            attachment_rows=attachment_rows,
        )
        rows = [
            row for row in rows
            if row.get("virus_name") == virus_name
            and row.get("canonical_target_id") == canonical_target_id
        ]
    else:
        try:
            source_payload = _load_protacability_source_payload(virus_name=virus_name, protein_type=protein_type)
        except RandyBackendError as exc:
            return jsonify({"data_available": False, "message": str(exc)}), exc.status_code

        if not source_payload.get("data_available"):
            return jsonify({"data_available": False, "message": "PROTACability data is not available."}), 404

        readiness_rows = source_payload.get("readiness_rows", [])
        warhead_rows = source_payload.get("warhead_rows", [])
        rows = _decorate_protacability_rows(
            source_payload.get("assessment_rows", []),
            collapse_labels=collapse_labels,
            readiness_rows=readiness_rows,
            warhead_rows=warhead_rows,
            attachment_rows=source_payload.get("attachment_rows", []),
        )
        rows = [row for row in rows if row.get("virus_name") == virus_name and (row.get("display_protein_type") or row.get("protein_type")) == protein_type]
    if min_score is not None:
        rows = [row for row in rows if _numeric_value(row.get("protacability_proxy_score"), -1) >= min_score]
    if ligand_context_class:
        rows = _apply_ligand_context_filter(rows, ligand_context_class)

    if not rows:
        return jsonify({"error": "Target detail not found"}), 404

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

    def _pick_context(context_class):
        matches = [row for row in structure_rows if row.get("ligand_context_class") == context_class]
        if not matches:
            return None
        best = max(matches, key=lambda r: _numeric_value(r.get("best_score")))
        ligand_payload = _load_protacability_source_payload(pdb_code=best.get("pdb_code"), include_inventory=True)
        ligand_rows = ligand_payload.get("ligand_inventory", []) if ligand_payload.get("data_available") else []
        preferred = _split_candidate_ligands(best.get("candidate_ligand_resnames_full"))
        ligand_record = _pick_representative_ligand_record(
            ligand_rows,
            preferred_ligands=preferred,
            allow_glycan=context_class == "glycan_only",
            preferred_chain=best.get("representative_chain_id")
        )
        ligand_resname = (ligand_record or {}).get("ligand_resname")
        return {
            "ligand_instance_id": (ligand_record or {}).get("ligand_instance_id"),
            "model_id": (ligand_record or {}).get("model_id"),
            "pdb_code": best.get("pdb_code"),
            "chain_id": best.get("representative_chain_id"),
            "ligand_resname": ligand_resname,
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
    ligand_inventory = []
    if active_pdb_code:
        try:
            ligand_payload = _load_protacability_source_payload(pdb_code=active_pdb_code, include_inventory=True)
            ligand_inventory = [
                {
                    "ligand_instance_id": row.get("ligand_instance_id"),
                    "model_id": row.get("model_id"),
                    "ligand_resname": row.get("ligand_resname"),
                    "ligand_chain": row.get("ligand_chain"),
                    "ligand_residue_id": row.get("ligand_residue_id"),
                    "ligand_atom_count": row.get("ligand_atom_count"),
                    "ligand_heavy_atom_count": row.get("ligand_heavy_atom_count"),
                }
                for row in ligand_payload.get("ligand_inventory", [])
                if row.get("pdb_code") == active_pdb_code
            ]
        except RandyBackendError:
            ligand_inventory = []
    attachment_lookup_row = {
        **target_summary,
        **(representative_ligand or {}),
        "pdb_code": active_pdb_code,
        "model_id": (representative_ligand or {}).get("model_id") or target_summary.get("model_id") or 0,
    }
    with connect_db_row() as conn:
        attachment_sites = _attachment_detail_payload(conn, attachment_lookup_row)
    return jsonify({
        "data_available": True,
        "target_summary": target_summary,
        "structure_rows": structure_rows,
        "ligand_groups": ligand_groups,
        "representative_contexts": representative_contexts,
        "representative_ligand": representative_ligand,
        "ligand_contexts": _serialize_ligand_contexts(ligand_inventory),
        "attachment_sites": attachment_sites,
    })


@app.route('/api/protacability/export')
def protacability_export():
    raw_export = (request.args.get("raw_export") or "").strip()
    if raw_export:
        mode = _normalized_backend_mode()
        if mode == "randy":
            try:
                raw_payload = randy_get("protacability/raw-table", params={"raw_export": raw_export}, max_bytes=None)
            except RandyBackendError as exc:
                return jsonify({"success": False, "message": str(exc)}), exc.status_code
            table_name = raw_payload.get("table_name")
            df = pd.DataFrame(raw_payload.get("rows", []))
        elif mode == "auto" and randy_available():
            try:
                raw_payload = randy_get("protacability/raw-table", params={"raw_export": raw_export}, max_bytes=None)
                table_name = raw_payload.get("table_name")
                df = pd.DataFrame(raw_payload.get("rows", []))
            except RandyBackendError:
                with connect_db_row() as conn:
                    try:
                        table_name, df = _local_protacability_raw_export(conn, raw_export)
                    except KeyError:
                        return jsonify({"success": False, "message": f"{raw_export} has not been imported yet."}), 404
        else:
            with connect_db_row() as conn:
                try:
                    table_name, df = _local_protacability_raw_export(conn, raw_export)
                except KeyError:
                    return jsonify({"success": False, "message": f"{raw_export} has not been imported yet."}), 404
        csv_buffer = io.StringIO()
        df.to_csv(csv_buffer, index=False)
        byte_buffer = BytesIO(csv_buffer.getvalue().encode("utf-8"))
        byte_buffer.seek(0)
        download_name = f"{table_name}.csv"
        return send_file(byte_buffer, mimetype="text/csv", as_attachment=True, download_name=download_name)

    mode = _normalized_backend_mode()
    if mode == "randy":
        try:
            payload = _remote_protacability_get("protacability/export-filtered", params=request.args, max_bytes=None)
        except RandyBackendError as exc:
            return jsonify({"success": False, "message": str(exc)}), exc.status_code
    elif mode == "auto" and randy_available():
        try:
            payload = _remote_protacability_get("protacability/export-filtered", params=request.args, max_bytes=None)
        except RandyBackendError:
            logging.warning("Falling back to local PROTACability filtered export payload")
            payload = None
    else:
        payload = None

    if payload is None:
        try:
            source_payload = _load_protacability_source_payload()
        except RandyBackendError as exc:
            return jsonify({"success": False, "message": str(exc)}), exc.status_code

        if not source_payload.get("data_available"):
            return jsonify({
                "success": False,
                "message": "PROTACability data has not been imported yet. Run tools/import_protacability_data.py after generating the expansion outputs."
            }), 404

        payload = _prepare_protacability_result_set_from_rows(
            source_payload.get("assessment_rows", []),
            source_payload.get("readiness_rows", []),
            source_payload.get("warhead_rows", []),
            request.args,
            export_all=True,
            attachment_rows=source_payload.get("attachment_rows", []),
        )

    df = pd.DataFrame(payload["rows"])

    csv_buffer = io.StringIO()
    df.to_csv(csv_buffer, index=False)
    byte_buffer = BytesIO(csv_buffer.getvalue().encode("utf-8"))
    byte_buffer.seek(0)
    download_name = f"protacability_{payload['view']}_filtered.csv"
    return send_file(
        byte_buffer,
        mimetype="text/csv",
        as_attachment=True,
        download_name=download_name
    )


@app.route('/api/debug/coordinate_ligand_presence/<pdb_code>/<ligand_code>')
def debug_coordinate_ligand_presence(pdb_code, ligand_code):
    pdb_upper = str(pdb_code or "").strip().upper()
    ligand_upper = str(ligand_code or "").strip().upper()
    chain = str(request.args.get("chain", "") or "").strip().upper()
    residue_id = str(request.args.get("residue_id", "") or "").strip()
    return jsonify({
        "pdb_code": pdb_upper,
        "requested_ligand": ligand_upper,
        "chain": chain,
        "residue_id": residue_id,
        "viewer_source": "direct_rcsb_mmcif",
        "viewer_coordinate_url": f"https://files.rcsb.org/download/{pdb_upper}.cif",
        "note": "Viewer coordinates are retrieved client-side and are not converted or cached by V-LiSEMOD.",
    })


@app.route('/api/debug/served_coordinate/<pdb_code>/<ligand_code>')
def debug_served_coordinate(pdb_code, ligand_code):
    pdb_upper = str(pdb_code or "").strip().upper()
    ligand_upper = str(ligand_code or "").strip().upper()
    chain = str(request.args.get("chain", "") or "").strip().upper()
    residue_id = str(request.args.get("residue_id", "") or "").strip().upper()

    return jsonify({
        "pdb_code": pdb_upper,
        "ligand_code": ligand_upper,
        "chain": chain,
        "residue_id": residue_id,
        "viewer_source": "direct_rcsb_mmcif",
        "viewer_coordinate_url": f"https://files.rcsb.org/download/{pdb_upper}.cif",
        "note": "No coordinate file is served, converted, or retained by V-LiSEMOD.",
    })


@app.route('/api/ligand_instance_metadata/<pdb_code>/<ligand_code>')
def ligand_instance_metadata(pdb_code, ligand_code):
    """Return stored occurrence metadata only; never retrieve or cache coordinates."""
    pdb_upper = str(pdb_code or "").strip().upper()
    ligand_upper = str(ligand_code or "").strip().upper()
    auth_chain = str(request.args.get("auth_chain", "") or "").strip().upper()
    auth_seq_id = str(request.args.get("auth_seq_id", "") or "").strip()
    if not re.fullmatch(r"[A-Z0-9]{4}", pdb_upper):
        return jsonify({"error": "Invalid PDB code"}), 400
    try:
        with connect_db_row() as conn:
            row = conn.execute(
                """
                SELECT i.label_asym_id, i.auth_asym_id, i.auth_seq_id, i.deposited_model_num
                FROM ligand_instances AS i
                JOIN structures AS s ON s.structure_id = i.structure_id
                WHERE UPPER(s.entry_id) = ?
                  AND UPPER(COALESCE(i.auth_comp_id, i.label_comp_id)) = ?
                  AND (? = '' OR UPPER(i.auth_asym_id) = ?)
                  AND (? = '' OR CAST(i.auth_seq_id AS TEXT) = ?)
                ORDER BY i.deposited_model_num, i.label_asym_id
                LIMIT 1
                """,
                (pdb_upper, ligand_upper, auth_chain, auth_chain, auth_seq_id, auth_seq_id),
            ).fetchone()
    except Exception:
        row = None
    if not row:
        return jsonify({"found": False}), 404
    return jsonify({
        "found": True,
        "label_asym_id": row["label_asym_id"],
        "auth_chain": row["auth_asym_id"],
        "auth_seq_id": row["auth_seq_id"],
        "model_id": row["deposited_model_num"],
    })


@app.route('/api/debug/ligand_context/<pdb_code>/<ligand_code>')
def debug_ligand_context(pdb_code, ligand_code):
    conn = connect_db_row()
    query_ctx = {
        "pdb_code_upper": str(pdb_code or "").strip().upper(),
        "ligand_code_upper": str(ligand_code or "").strip().upper(),
        "ligand_alias_candidates": make_ligand_code_aliases(ligand_code),
        "chain_upper": str(request.args.get("chain", "") or "").strip().upper(),
        "residue_id_string": str(request.args.get("residue_id", "") or "").strip(),
        "residue_id_int_if_possible": None,
        "virus_name": str(request.args.get("virus_name", "") or "").strip(),
        "protein_type": str(request.args.get("protein_type", "") or "").strip(),
        "debug_full": request.args.get("debug_full", "0") == "1",
    }
    if query_ctx["residue_id_string"]:
        try:
            query_ctx["residue_id_int_if_possible"] = int(query_ctx["residue_id_string"])
        except (TypeError, ValueError):
            query_ctx["residue_id_int_if_possible"] = None

    app.logger.info(
        "[ligand-debug] pdb=%s ligand=%s chain=%s residue=%s aliases=%s",
        query_ctx["pdb_code_upper"],
        query_ctx["ligand_code_upper"],
        query_ctx["chain_upper"] or "(none)",
        query_ctx["residue_id_string"] or "(none)",
        query_ctx["ligand_alias_candidates"],
    )

    table_names = [
        "Ligand_Arp_Diagram",
        "Ligand_Atoms_Smiles",
        "ligand_atoms",
        "SMILES_MAP_PDB",
        "solvent_exposed_atoms",
        "RUPLEY_SASA_DATA",
        "Arpeggio_Contacts_Data",
        "protacability_assessment",
        "protacability_ligand_inventory",
        "protacability_lysine_proximity",
        "v2_attachment_site_summary",
        "v2_attachment_site_candidates",
        "v2_attachment_site_high_priority",
    ]
    tables = []
    for table_name in table_names:
        payload = _build_debug_table_payload(conn, table_name, query_ctx)
        app.logger.info("[ligand-debug] %s rows=%s error=%s", table_name, payload.get("row_count"), payload.get("error"))
        tables.append(payload)

    coordinate_presence_summary = []
    coord_files = _candidate_coordinate_files(query_ctx["pdb_code_upper"])
    for path in coord_files[:40]:
        try:
            with open(path, "r", encoding="utf-8", errors="ignore") as handle:
                text = handle.read()
            coordinate_presence_summary.append(
                _scan_coordinate_source(
                    path,
                    text,
                    query_ctx["ligand_alias_candidates"],
                    chain=query_ctx["chain_upper"],
                    residue_id=query_ctx["residue_id_string"],
                )
            )
        except Exception as exc:
            coordinate_presence_summary.append({
                "source": path,
                "exists": os.path.isfile(path),
                "error": str(exc),
            })
    conn.close()
    return jsonify({
        "normalized_query": {
            "input": {
                "pdb_code": pdb_code,
                "ligand_code": ligand_code,
                "chain": request.args.get("chain"),
                "residue_id": request.args.get("residue_id"),
            },
            "normalized": {
                "pdb_code_upper": query_ctx["pdb_code_upper"],
                "ligand_code_upper": query_ctx["ligand_code_upper"],
                "ligand_alias_candidates": query_ctx["ligand_alias_candidates"],
                "chain_upper": query_ctx["chain_upper"],
                "residue_id_string": query_ctx["residue_id_string"],
                "residue_id_int_if_possible": query_ctx["residue_id_int_if_possible"],
            }
        },
        "filters": {
            "virus_name": query_ctx["virus_name"],
            "protein_type": query_ctx["protein_type"],
            "debug_full": query_ctx["debug_full"],
        },
        "tables": tables,
        "coordinate_presence_summary": coordinate_presence_summary,
        "recommended_structure_url": url_for(
            "serve_coordinate_for_viewer",
            pdb_code=query_ctx["pdb_code_upper"],
            ligand_code=query_ctx["ligand_code_upper"],
        ),
    })


@app.route('/api/debug/protacability_detail_payload/<pdb_code>')
def debug_protacability_detail_payload(pdb_code):
    conn = connect_db_row()
    chain_id = str(request.args.get("chain_id", "") or "").strip()
    virus_name = str(request.args.get("virus_name", "") or "").strip()
    protein_type = str(request.args.get("protein_type", "") or "").strip()
    view_type = str(request.args.get("view_type", "chain") or "chain").strip().lower()

    assessment_rows = _decorate_protacability_rows(_load_protacability_assessment_rows(conn, pdb_code=pdb_code), collapse_labels=True)
    if chain_id:
        assessment_rows = [row for row in assessment_rows if str(row.get("chain_id") or "").strip().upper() == chain_id.upper()]
    if virus_name:
        assessment_rows = [row for row in assessment_rows if str(row.get("virus_name") or "").strip() == virus_name]
    if protein_type:
        assessment_rows = [row for row in assessment_rows if str(row.get("display_protein_type") or row.get("protein_type") or "").strip() == protein_type]

    assessment = max(assessment_rows, key=_row_priority_key) if assessment_rows else None
    representative_chain_id = (assessment or {}).get("chain_id") or chain_id

    ligand_inventory = [
        dict(row) for row in conn.execute(
            """
            SELECT
                ligand_instance_id,
                model_id,
                ligand_resname,
                ligand_chain,
                ligand_residue_id,
                ligand_insertion_code,
                ligand_atom_count,
                ligand_heavy_atom_count,
                centroid_x,
                centroid_y,
                centroid_z
            FROM protacability_ligand_inventory
            WHERE pdb_code = ?
            ORDER BY ligand_resname, ligand_chain, ligand_residue_id
            """,
            (pdb_code,)
        ).fetchall()
    ] if _table_exists(conn, "protacability_ligand_inventory") else []

    preferred_ligands = _split_candidate_ligands((assessment or {}).get("candidate_ligand_resnames_full"))
    representative_ligand = _pick_representative_ligand_record(
        ligand_inventory,
        preferred_ligands=preferred_ligands,
        allow_glycan=(assessment or {}).get("ligand_context_class") == "glycan_only",
        preferred_chain=representative_chain_id
    )
    ligand_contexts = _serialize_ligand_contexts(ligand_inventory)
    selected_default_ligand_context = representative_ligand or (ligand_contexts[0] if ligand_contexts else None)
    viewer_payload_expected = {
        "pdbCode": str(pdb_code or "").strip().upper(),
        "chainId": representative_chain_id or "",
        "ligandResname": (selected_default_ligand_context or {}).get("ligand_resname", ""),
        "ligandChain": (selected_default_ligand_context or {}).get("ligand_chain", ""),
        "ligandResidueId": str((selected_default_ligand_context or {}).get("ligand_residue_id") or ""),
    }
    conn.close()
    return jsonify({
        "pdb_code": str(pdb_code or "").strip().upper(),
        "view_type": view_type,
        "chain_id": chain_id,
        "virus_name": virus_name,
        "protein_type": protein_type,
        "assessment_row": dict(assessment) if assessment else None,
        "representative_chain_id": representative_chain_id,
        "representative_ligand": representative_ligand,
        "ligand_contexts": ligand_contexts,
        "ligand_inventory": ligand_inventory,
        "selected_default_ligand_context": selected_default_ligand_context,
        "viewer_payload_expected": viewer_payload_expected,
    })

    
# if __name__ == '__main__':
#     app.run(debug=True)


# if __name__ == '__main__':
#     app.run(host='0.0.0.0', port=5002, debug=True)

if __name__ == "__main__":
    app.run(host="127.0.0.1", port=5003, debug=True)

    pass
