import os 
import io
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
    "RUPLEY_SASA_DATA",
    "SMILES_MAP_PDB",
)
PYMOL_REQUIRED_TABLES = (
    "ligand_atoms",
    "Functional_Group_Atoms",
    "receptor_binding_pocket",
    "distal_atoms",
    "RUPLEY_SASA_DATA",
)


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


def use_randy_backend():
    return _normalized_backend_mode() == "randy"


def randy_available():
    return bool(_randy_base_url() and _randy_api_token())


def randy_get(path, params=None):
    if not randy_available():
        raise RandyBackendError("RANDY API is not configured.", status_code=500)

    url = f"{_randy_base_url()}/{str(path or '').lstrip('/')}"
    headers = {"Authorization": f"Bearer {_randy_api_token()}"}

    try:
        response = requests.get(url, params=params, headers=headers, timeout=10)
    except requests.RequestException as exc:
        raise RandyBackendError(f"RANDY API request failed: {exc}", status_code=502) from exc

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


def randy_post(path, json=None):
    if not randy_available():
        raise RandyBackendError("RANDY API is not configured.", status_code=500)

    url = f"{_randy_base_url()}/{str(path or '').lstrip('/')}"
    headers = {"Authorization": f"Bearer {_randy_api_token()}"}

    try:
        response = requests.post(url, json=json, headers=headers, timeout=20)
    except requests.RequestException as exc:
        raise RandyBackendError(f"RANDY API request failed: {exc}", status_code=502) from exc

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
            "SELECT name FROM sqlite_master WHERE type = 'table' AND name IN ({})".format(
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
    structure_url = None
    if images:
        primary = images[0]
        structure_url = url_for(
            "serve_coordinate_for_viewer",
            pdb_code=str(primary.get("pdb_id") or "").strip().upper(),
            ligand_code=str(primary.get("ligand_code") or "").strip().upper()
        )
    return render_template(
        'display_images.html',
        images=images,
        chain_residues=chain_residue_data,
        structure_url=structure_url
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
def write_pymol_script(pdb_code, ligand_name, ligand_chain, options):
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
                atom_selection = " or ".join([f"id {atom[0]}" for atom in atoms])
                script.write(f"select {fg_name}, chain {ligand_chain} and ({atom_selection})\n")
                script.write(f"create {fg_name}_Object, {fg_name}\n")
                script.write(f"show sticks, {fg_name}_Object\n")
                script.write(f"color magenta, {fg_name}_Object\n")
        
        # Binding pocket
        if options.get('binding_pocket'):
            binding_pocket_atoms = options['binding_pocket']
            script.write("create Binding_Pocket, none\n")
            for atom in binding_pocket_atoms:
                script.write(f"select temp, chain {atom[0]} and id {atom[1]}\n")
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
                script.write(f"select temp, chain {atom[0]} and id {atom[1]}\n")
                script.write("create Distal_Atoms, Distal_Atoms or temp\n")
                script.write("delete temp\n")
            script.write("show spheres, Distal_Atoms\n")
            script.write("color blue, Distal_Atoms\n")
        
        # Solvent exposed atoms
        if options.get('solvent_exposed_atoms'):
            solvent_exposed_atoms = options['solvent_exposed_atoms']
            logging.debug(f"Solvent-exposed atoms (SASA): {solvent_exposed_atoms}")
            script.write("create Solvent_Exposed_Atoms, none\n")
            for atom_id, atom_chain in solvent_exposed_atoms:
                script.write(f"select temp, chain {atom_chain} and id {atom_id}\n")
                script.write("create Solvent_Exposed_Atoms, Solvent_Exposed_Atoms or temp\n")
                script.write("delete temp\n")
            script.write("show spheres, Solvent_Exposed_Atoms\n")
            script.write("color firebrick, Solvent_Exposed_Atoms\n")

        # Hydrated atoms
        if options.get('hydrated_atoms'):
            hydrated_atoms = options['hydrated_atoms']
            script.write("create Hydrated_Atoms, none\n")
            for atom in hydrated_atoms:
                script.write(f"select temp, chain {atom[0]} and id {atom[1]}\n")
                script.write("create Hydrated_Atoms, Hydrated_Atoms or temp\n")
                script.write("delete temp\n")
            script.write("show sticks, Hydrated_Atoms\n")
            script.write("color cyan, Hydrated_Atoms\n")

        # Rupley SASA atoms
        if options.get('rupley_sasa'):
            rupley_sasa_atoms = options['rupley_sasa']
            script.write("create RUPLEY_SASA, none\n")
            for atom in rupley_sasa_atoms:
                script.write(f"select temp, chain {atom[1]} and id {atom[0]}\n")
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

        options = {}
        if option_flags['functional_groups']:
            options['functional_groups'] = {
                group_name: [
                    (atom.get("atom_id"), atom.get("exact_atom"), atom.get("atom_type"))
                    for atom in atoms
                ]
                for group_name, atoms in (remote_payload.get("functional_groups") or {}).items()
            }
        if option_flags['binding_pocket']:
            options['binding_pocket'] = [
                (row.get("residue_chain"), row.get("residue_number"))
                for row in remote_payload.get("binding_pocket", [])
            ]
        if option_flags['distal_atoms']:
            options['distal_atoms'] = [
                (row.get("chain"), row.get("atom_id"))
                for row in remote_payload.get("distal_atoms", [])
            ]
        if option_flags['solvent_exposed_atoms']:
            options['solvent_exposed_atoms'] = [
                (row.get("atom_id"), row.get("chain"))
                for row in remote_payload.get("solvent_exposed_atoms", [])
            ]
        if option_flags['hydrated_atoms']:
            options['hydrated_atoms'] = [
                (row.get("chain"), row.get("atom_id"))
                for row in remote_payload.get("hydrated_atoms", [])
            ]
        if option_flags['rupley_sasa']:
            options['rupley_sasa'] = [
                (row.get("atom_id"), row.get("chain"))
                for row in remote_payload.get("rupley_sasa", [])
            ]
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
                        SELECT atom_id, chain
                        FROM RUPLEY_SASA_DATA
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

    pymol_script_path = write_pymol_script(pdb_code, ligand_name, ligand_chain, options)
    return send_file(pymol_script_path, as_attachment=True)


def _local_get_viruses_payload():
    conn = sqlite3.connect('viral_data.db')
    cursor = conn.cursor()
    cursor.execute('SELECT DISTINCT virus_name FROM ligand_atoms')
    viruses = [row[0] for row in cursor.fetchall()]
    conn.close()
    return {'viruses': viruses}


def _local_get_pdb_codes_payload(virus_name):
    conn = sqlite3.connect('viral_data.db')
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
    conn = sqlite3.connect('viral_data.db')
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
    conn = sqlite3.connect('viral_data.db')
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
    conn = sqlite3.connect('viral_data.db')
    cursor = conn.cursor()
    cursor.execute('SELECT DISTINCT ligand FROM Ligand_Atoms_Smiles')
    ligands = [row[0] for row in cursor.fetchall()]
    conn.close()
    return {'ligands': ligands}


def _local_get_viruses_by_ligand_payload(ligand_code):
    conn = sqlite3.connect('viral_data.db')
    cursor = conn.cursor()
    cursor.execute('SELECT DISTINCT Virus_Name FROM Arpeggio_Contacts_Data WHERE Ligand = ?', (ligand_code,))
    viruses = [row[0] for row in cursor.fetchall()]
    conn.close()
    return {'viruses': viruses}


def _local_get_pdb_residue_by_ligand_payload(ligand_code):
    conn = sqlite3.connect('viral_data.db')
    cursor = conn.cursor()
    cursor.execute(
        '''
        SELECT DISTINCT pdb_id, chain, ligand_id
        FROM Ligand_Arp_Diagram
        WHERE ligand = ?
        ''',
        (ligand_code,),
    )
    pairs = [{'pdb_id': row[0], 'chain': row[1], 'ligand_id': row[2]} for row in cursor.fetchall()]
    conn.close()
    return {'pairs': pairs}


def _local_get_sasa_chains_payload(pdb_code, ligand_name):
    conn = sqlite3.connect('viral_data.db')
    cursor = conn.cursor()
    cursor.execute(
        '''
        SELECT atom_id, chain
        FROM RUPLEY_SASA_DATA
        WHERE pdb_id = ? AND ligand = ?
        ''',
        (pdb_code, ligand_name),
    )
    sasa_chains = cursor.fetchall()
    conn.close()
    return sasa_chains


def _local_get_pdb_mapping_payload(ligand_code):
    conn = sqlite3.connect('viral_data.db')
    cursor = conn.cursor()
    cursor.execute(
        '''
        SELECT DISTINCT pdb_id, chain, ligand_id, virus_name, ligand
        FROM Ligand_Atoms_Smiles
        WHERE ligand = ?
        ''',
        (ligand_code,),
    )

    pdb_mapping = {}
    for row in cursor.fetchall():
        pdb_id = row[0]
        chain = row[1]
        ligand_id = row[2]
        virus_name = row[3]
        ligand = row[4]
        unique_key = f"{pdb_id}-{ligand_id}-{chain}"
        if unique_key not in pdb_mapping:
            pdb_mapping[unique_key] = {
                'pdb_id': pdb_id,
                'ligand_id': ligand_id,
                'chain': chain,
                'virus_name': virus_name,
                'ligand': ligand
            }

    conn.close()
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
def get_interaction_data(pdb_id, ligand, ligand_id, chain):
    conn = sqlite3.connect('viral_data.db')
    query = '''
        SELECT * FROM Arpeggio_Contacts_Data
        WHERE pdb_id = ? AND ligand = ? AND ligand_id = ? AND chain = ?
    '''
    df = pd.read_sql(query, conn, params=(pdb_id, ligand, ligand_id, chain))
    conn.close()
    return df

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
    sns.barplot(x=interaction_counts.index, y=interaction_counts.values, palette=colors)
    
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
        df = df[df['Contact'] != 'proximal']

    # Merge 'weak_polar' into 'polar' and 'vdw_clash' into 'vdw'
    df['Contact'] = df['Contact'].replace({'weak_polar': 'polar', 'vdw_clash': 'vdw', 'weak_hbond': 'hbond'})

    # Create a new column that combines atom_id and Ligand_Atom for labeling
    df['Atom_Label'] = df['atom_id'].astype(str) + ' (' + df['exact_atom'] + ')'

    # Filter out rows where Atom_Label contains 'nan'
    df = df[df['Atom_Label'].notna() & df['Atom_Label'] != 'nan']

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
    # Get the input from the form submission (AJAX POST request)
    data = request.json
    pdb_id = data['pdb_id']
    ligand = data['ligand']
    ligand_id = data['ligand_id']  # Corresponds to ligand_id
    chain = data['chain']
    
    # Query the database to get the relevant data
    df = get_interaction_data(pdb_id, ligand, ligand_id, chain)

    if df.empty:
        return jsonify({'error': 'No data found for the selected inputs'}), 404

    # Preprocess the dataframe for filtering
    df_clean = preprocess_interactions(df)  # Filtered interactions excluding 'proximal', merging 'vdw_clash' and 'weak_polar'
    df_clean = filter_valid_atom_ids(df_clean)  # Ensure valid atom_ids
    
    # Filter the 'proximal' interactions and ensure they have valid atom_ids
    df_proximal = df[df['Contact'] == 'proximal']
    df_proximal_clean = filter_valid_atom_ids(df_proximal)

    # Add 'proximal' interactions back to the dataset for specific plots
    df_with_proximal = pd.concat([df_clean, df_proximal_clean])

    # Ensure CHARTS_DIR exists
    if not os.path.exists(CHARTS_DIR):
        os.makedirs(CHARTS_DIR)

    # File paths for the generated charts
    pie_chart_file = os.path.join(CHARTS_DIR, 'pie_chart.png')
    bar_chart_file = os.path.join(CHARTS_DIR, 'bar_chart.png')
    scatter_chart_file = os.path.join(CHARTS_DIR, 'scatter_chart.png')
    interactions_per_atom_with_proximal_file = os.path.join(CHARTS_DIR, 'interactions_per_atom_with_proximal.png')
    interactions_per_atom_without_proximal_file = os.path.join(CHARTS_DIR, 'interactions_per_atom_without_proximal.png')

     # Generate charts with PDB and Ligand info in the titles
    plot_pie_chart(df_clean, pie_chart_file, pdb_id, ligand)
    plot_bar_chart(df_clean, bar_chart_file, pdb_id, ligand)
    plot_scatter_chart(df_clean, scatter_chart_file, pdb_id, ligand)
    plot_interactions_per_atom(df_with_proximal, interactions_per_atom_with_proximal_file,  pdb_id, ligand, exclude_proximal=False)
    plot_interactions_per_atom(df_clean, interactions_per_atom_without_proximal_file,  pdb_id, ligand, exclude_proximal=True)

    # Return the file paths of the charts to be displayed in the carousel
    return jsonify({
        'chart_paths': [
            pie_chart_file,
            bar_chart_file,
            scatter_chart_file,
            interactions_per_atom_with_proximal_file,
            interactions_per_atom_without_proximal_file
        ]
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
    with _connect_local_db(LIGAND_IMAGE_REQUIRED_TABLES) as conn:
        cur = conn.cursor()
        cur.execute(
            """
            SELECT atom_id, exact_atom
            FROM RUPLEY_SASA_DATA
            WHERE virus_name = ? AND pdb_id = ? AND ligand = ? AND chain = ?
            """,
            (virus_name, pdb_id, ligand_id, chain),
        )
        sasa_atoms = cur.fetchall()

        sasa_smiles_indices = []
        for atom_id, exact_atom in sasa_atoms:
            cur.execute(
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
                (virus_name, pdb_id, ligand_id, chain, atom_id, exact_atom),
            )

            row = cur.fetchone()
            if row and row[0] is not None:
                sasa_smiles_indices.append(int(row[0]))

    highlight_solvent_exposed_atoms_from_indices(molecule, sasa_smiles_indices, output_svg)


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
    data = request.json
    selected_pdbs = data['pdb_ids']  # Format like ['PDB_ID-Residue-Chain', ...]
    ligand = data['ligand']

    interactions_data = []
    smiles_interactions_data = []

    conn = sqlite3.connect('viral_data.db')
    cursor = conn.cursor()

    for unique_key in selected_pdbs:
        pdb_id, ligand_id, chain = unique_key.split('-')

        # Query Arpeggio_Contacts_Data and join with SMILES_MAP_PDB to get smiles_atom_index
        query = '''
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
    '''

        df = pd.read_sql(query, conn, params=(pdb_id, ligand_id, chain, ligand))

        # Replace NaN with None for proper JSON serialization
        df_clean = df.replace({np.nan: None})

        # Apply filtering and preprocessing to clean the data
        df_clean = preprocess_interactions(df_clean)
        df_clean = filter_valid_atom_ids(df_clean)

        if not df_clean.empty:
            interactions_data.append({
                'pdb_id': pdb_id,
                'virus_name': df_clean['virus_name'].iloc[0],
                'ligand_id': int(df_clean['ligand_id'].iloc[0]),
                'interactions': df_clean.to_dict(orient='records')
            })
        else:
            logging.warning(f"No data found for PDB ID: {pdb_id}, Residue ID: {ligand_id}, Chain: {chain}, Ligand: {ligand}")

        # Process SMILES Atom Index data
        if 'smiles_atom_index' in df.columns:
            smiles_df = df[['pdb_id', 'Contact', 'smiles_atom_index']].dropna(subset=['smiles_atom_index'])
            
            # Convert the SMILES atom index to integers and filter out NaN values
            smiles_df['smiles_atom_index'] = smiles_df['smiles_atom_index'].astype(float).astype('Int64').replace({pd.NA: None})
            
            smiles_interactions_data.append({
                'pdb_id': pdb_id,
                'interactions': smiles_df.to_dict(orient='records')
            })
        else:
            print(f"'smiles_atom_index' not found for PDB ID: {pdb_id}")

    conn.close()

    return jsonify({
        'interactions_data': interactions_data,
        'smiles_interactions_data': smiles_interactions_data
    })

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
    conn = sqlite3.connect('viral_data.db')
    cursor = conn.cursor()
    cursor.execute("SELECT virus_name, pdb_id, ligand, chain, ligand_id FROM Ligand_Arp_Diagram")
    rows = cursor.fetchall()
    options = [{'value': f"{row[1]}-{row[2]}", 'text': f"{row[0]}, {row[1]}, {row[2]}, {row[3]}, {row[4]}"} for row in rows]
    conn.close()
    return jsonify(options=options)



# Function to get the SMILES string based on ligand_id
def get_smiles_from_identifier(ligand):
    conn = sqlite3.connect('viral_data.db')
    cursor = conn.cursor()
    cursor.execute("SELECT smiles FROM Ligand_Atoms_Smiles WHERE ligand = ?", (ligand,))
    result = cursor.fetchone()
    conn.close()
    if result:
        return result[0]
    return None

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
    return sqlite3.connect('viral_data.db')


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
    conn = sqlite3.connect('viral_data.db')
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
    limit = min(max(args.get("limit", type=int) or 50, 1), 100)
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
        "SELECT name FROM sqlite_master WHERE type='table' AND name = ?",
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
    conn = connect_db()
    cursor = conn.cursor()
    cursor.execute("SELECT DISTINCT virus_name FROM Virus_Proteins")
    virus_names = [row[0] for row in cursor.fetchall()]
    conn.close()
    return jsonify(virus_names)

@app.route('/get_protein_types_list_distinct', methods=['GET'])
def get_protein_types_list_distinct():
    conn = connect_db()
    cursor = conn.cursor()
    cursor.execute("SELECT DISTINCT protein FROM Virus_Proteins")
    protein_types = [row[0] for row in cursor.fetchall()]
    conn.close()
    return jsonify(protein_types)




@app.route('/get_pdbs_for_virus_protein', methods=['POST'])
def get_pdbs_for_virus_protein():
    data = request.json
    virus_name = data['virus_name']
    protein_types = data['protein_types']
    ligand_filter = data.get('ligand', None)  # Optional ligand filter

    conn = sqlite3.connect('viral_data.db')
    cursor = conn.cursor()

    # Basic query for PDB codes based on virus and protein types
    placeholders = ', '.join(['?'] * len(protein_types))
    query = f"SELECT DISTINCT pdb_id FROM Virus_Proteins WHERE virus_name = ? AND protein IN ({placeholders})"
    params = [virus_name] + protein_types

    # If a ligand filter is provided, include it in the query
    if ligand_filter:
        # Subquery to filter PDBs based on ligand or synonym presence
        query += '''
            AND pdb_id IN (
                SELECT pdb_id FROM Ligand_ARP_Diagram WHERE ligand = ? OR ligand IN (
                    SELECT synonym FROM Ligand_Synonyms WHERE ligand = ?
                )
            )
        '''
        params += [ligand_filter, ligand_filter]

    cursor.execute(query, params)
    pdb_codes = [row[0] for row in cursor.fetchall()]
    conn.close()

    return jsonify({'pdb_codes': pdb_codes})





@app.route('/export_data_to_excel', methods=['POST'])
def export_data_to_excel():
    data = request.json
    pdb_codes = data.get('pdb_codes', [])
    data_sets = data.get('data_sets', [])
    output_dir = 'output_files'  # Directory where CSV files will be saved

    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    if not pdb_codes or not data_sets:
        return jsonify({'success': False, 'error': 'No PDB codes or datasets provided.'}), 400

    allowed_data_sets = set(get_available_export_data_sets())
    invalid_data_sets = [data_set for data_set in data_sets if data_set not in allowed_data_sets]
    if invalid_data_sets:
        return jsonify({'success': False, 'error': f'Unsupported data sets requested: {invalid_data_sets}'}), 400

    conn = connect_db()
    response = {}
    excel_file_path = os.path.join(output_dir, "combined_data.xlsx")

    try:
        # Create an ExcelWriter object to write multiple sheets
        with pd.ExcelWriter(excel_file_path, engine='xlsxwriter') as writer:
            for data_set in data_sets:
                placeholders = ', '.join(['?'] * len(pdb_codes))
                query = data_set_queries[data_set].format(placeholders=placeholders)
                try:
                    df = pd.read_sql(query, conn, params=pdb_codes)
                    if not df.empty:
                        # Write each dataset to a CSV file
                        csv_file_path = os.path.join(output_dir, f"{data_set.replace(' ', '_')}.csv")
                        df.to_csv(csv_file_path, index=False)
                        response[data_set] = f"File saved: {csv_file_path}"
                        
                        # Write each dataset as a sheet in the Excel file
                        df.to_excel(writer, sheet_name=data_set[:30], index=False)
                    else:
                        response[data_set] = "No data to save."
                except Exception as e:
                    response[data_set] = f"Error querying database: {str(e)}"
                    logging.error("Error querying database: %s", e)
        
        # Create a ZIP archive that includes the CSV files and the combined Excel file
        zip_buffer = BytesIO()  # Use an in-memory buffer to store the ZIP
        with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zipf:
            # Add each CSV file to the ZIP
            for data_set in data_sets:
                csv_file_path = os.path.join(output_dir, f"{data_set.replace(' ', '_')}.csv")
                if os.path.exists(csv_file_path):
                    zipf.write(csv_file_path, os.path.basename(csv_file_path))

            # Add the combined Excel file to the ZIP
            zipf.write(excel_file_path, os.path.basename(excel_file_path))

        zip_buffer.seek(0)  # Move back to the start of the BytesIO buffer

        # Clean up temporary files
        shutil.rmtree(output_dir)

        # Send the ZIP archive as a downloadable file
        return send_file(zip_buffer, as_attachment=True, download_name='data_sets.zip', mimetype='application/zip')

    except Exception as e:
        logging.error("Failed during file creation: %s", e)
        return jsonify({'success': False, 'error': str(e)}), 500
    finally:
        conn.close()



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
    return render_template(
        'query_protein_virus.html',
        available_export_data_sets=get_available_export_data_sets(),
        protacability_data_available=protacability_tables_available()
    )


############################
###Scripts for Ligand Querying



@app.route('/get_ligands_with_synonyms', methods=['GET'])
def get_ligands_with_synonyms():
    conn = sqlite3.connect('viral_data.db')
    cursor = conn.cursor()
    cursor.execute("SELECT ligand, synonym FROM Ligand_Synonyms")
    ligands = [{'ligand_code': row[0], 'synonym': row[1]} for row in cursor.fetchall()]
    conn.close()
    return jsonify(ligands)





@app.route("/get_ligand_info/<ligand_code>", methods=["GET"])
def get_ligand_info(ligand_code):
    """
    Fetch ligand details (SMILES, PDB ID, molecular weight) from the Ligand_Atoms_Smiles.csv database.
    """
    conn = sqlite3.connect("viral_data.db")  # Connect to your database
    cursor = conn.cursor()

    cursor.execute("""
        SELECT ligand, pdb_id, smiles, molecular_weight 
        FROM Ligand_Atoms_Smiles 
        WHERE ligand = ?
    """, (ligand_code,))
    
    ligand_data = cursor.fetchone()
    conn.close()

    if ligand_data:
        return jsonify({
            "ligand": ligand_data[0],
            "pdb_id": ligand_data[1],
            "smiles": ligand_data[2],
            "molecular_weight": ligand_data[3]
        })
    else:
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
    return render_template(
        'protacability_assessment.html',
        protacability_data_available=protacability_tables_available()
    )


@app.route('/api/protacability/filters')
def protacability_filters():
    conn = connect_db_row()
    if not protacability_tables_available(conn):
        conn.close()
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
    payload = _build_protacability_filter_options_payload(conn, request.args)
    conn.close()
    return jsonify(payload)


@app.route('/api/protacability/filter_options')
def protacability_filter_options():
    conn = connect_db_row()
    if not protacability_tables_available(conn):
        conn.close()
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

    payload = _build_protacability_filter_options_payload(conn, request.args)
    conn.close()
    return jsonify(payload)


@app.route('/api/protacability/search')
def protacability_search():
    conn = connect_db_row()
    if not protacability_tables_available(conn):
        conn.close()
        return jsonify({
            "data_available": False,
            "message": "PROTACability data has not been imported yet. Run tools/import_protacability_data.py after generating the expansion outputs.",
            "rows": [],
            "summary": {}
        })

    payload = _prepare_protacability_result_set(conn, request.args)
    conn.close()
    return jsonify({
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
    })


@app.route('/api/protacability/detail/<pdb_code>/<chain_id>')
def protacability_detail(pdb_code, chain_id):
    conn = connect_db_row()
    if not protacability_tables_available(conn):
        conn.close()
        return jsonify({
            "data_available": False,
            "message": "PROTACability data has not been imported yet. Run tools/import_protacability_data.py after generating the expansion outputs."
        }), 404

    readiness_rows, warhead_rows = _load_protacability_enrichment_tables(conn)
    raw_assessment_rows = conn.execute(
        "SELECT * FROM protacability_assessment WHERE pdb_code = ? AND chain_id = ?",
        (pdb_code, chain_id)
    ).fetchall()

    assessment_rows = _decorate_protacability_rows(
        raw_assessment_rows,
        collapse_labels=True,
        readiness_rows=readiness_rows,
        warhead_rows=warhead_rows,
    )
    assessment = max(assessment_rows, key=_row_priority_key) if assessment_rows else None

    if assessment is None:
        conn.close()
        return jsonify({"error": "Assessment row not found"}), 404

    lysine_rows = [
        dict(row) for row in conn.execute(
            """
            SELECT
                lys_residue_id,
                lys_observed_index,
                lysine_sasa_a2,
                is_surface_exposed,
                nearest_ligand_resname,
                nearest_ligand_distance_a,
                linker_site_class
            FROM protacability_lysine_proximity
            WHERE pdb_code = ? AND chain_id = ?
            ORDER BY lys_residue_id
            """,
            (pdb_code, chain_id)
        ).fetchall()
    ]

    ligand_inventory = [
        dict(row) for row in conn.execute(
            """
            SELECT
                ligand_resname,
                ligand_chain,
                ligand_residue_id,
                ligand_atom_count,
                NULL AS ligand_heavy_atom_count,
                centroid_x,
                centroid_y,
                centroid_z
            FROM protacability_ligand_inventory
            WHERE pdb_code = ?
            ORDER BY ligand_resname, ligand_chain, ligand_residue_id
            """,
            (pdb_code,)
        ).fetchall()
    ]

    related_chains = [
        dict(row) for row in conn.execute(
            """
            SELECT
                chain_id,
                protacability_proxy_score,
                protacability_tier,
                candidate_ligand_count,
                exposed_lys_count
            FROM protacability_assessment
            WHERE pdb_code = ?
            ORDER BY chain_id
            """,
            (pdb_code,)
        ).fetchall()
    ]
    conn.close()

    return jsonify({
        "data_available": True,
        "assessment": dict(assessment),
        "lysine_rows": lysine_rows,
        "ligand_inventory": ligand_inventory,
        "ligand_contexts": _serialize_ligand_contexts(ligand_inventory),
        "related_chains": related_chains
    })


@app.route('/api/protacability/structure_detail/<pdb_code>')
def protacability_structure_detail(pdb_code):
    conn = connect_db_row()
    if not protacability_tables_available(conn):
        conn.close()
        return jsonify({
            "data_available": False,
            "message": "PROTACability data has not been imported yet. Run tools/import_protacability_data.py after generating the expansion outputs."
        }), 404

    collapse_labels = _protacability_collapse_labels(request.args.get("collapse_labels"))
    virus_name = (request.args.get("virus_name") or "").strip()
    protein_type = (request.args.get("protein_type") or "").strip()

    readiness_rows, warhead_rows = _load_protacability_enrichment_tables(conn)
    decorated_rows = _decorate_protacability_rows(
        _load_protacability_assessment_rows(conn, pdb_code=pdb_code),
        collapse_labels=collapse_labels,
        readiness_rows=readiness_rows,
        warhead_rows=warhead_rows,
    )
    if virus_name:
        decorated_rows = [row for row in decorated_rows if row.get("virus_name") == virus_name]
    if protein_type:
        decorated_rows = [row for row in decorated_rows if (row.get("display_protein_type") if collapse_labels else row.get("protein_type")) == protein_type]

    if not decorated_rows:
        conn.close()
        return jsonify({"error": "Structure summary not found"}), 404

    summary_rows = _group_structure_rows(decorated_rows)
    summary_row = max(summary_rows, key=lambda row: (_numeric_value(row.get("best_score")), _numeric_value(row.get("best_exposed_lys_fraction"))))
    chain_rows = _dedupe_display_chain_rows(decorated_rows)
    chain_rows, _ = _sort_protacability_rows(chain_rows, "chains", "protacability_proxy_score_desc")
    representative_chain = summary_row.get("representative_chain_id")

    lysine_rows = [
        dict(row) for row in conn.execute(
            """
            SELECT
                lys_residue_id,
                lys_observed_index,
                lysine_sasa_a2,
                is_surface_exposed,
                nearest_ligand_resname,
                nearest_ligand_distance_a,
                linker_site_class
            FROM protacability_lysine_proximity
            WHERE pdb_code = ? AND chain_id = ?
            ORDER BY lys_residue_id
            """,
            (pdb_code, representative_chain)
        ).fetchall()
    ]

    ligand_inventory = [
        dict(row) for row in conn.execute(
            """
            SELECT
                ligand_resname,
                ligand_chain,
                ligand_residue_id,
                ligand_atom_count,
                NULL AS ligand_heavy_atom_count,
                centroid_x,
                centroid_y,
                centroid_z
            FROM protacability_ligand_inventory
            WHERE pdb_code = ?
            ORDER BY ligand_resname, ligand_chain, ligand_residue_id
            """,
            (pdb_code,)
        ).fetchall()
    ]
    preferred_ligands = _split_candidate_ligands(summary_row.get("candidate_ligand_resnames_full"))
    representative_ligand = _pick_representative_ligand_record(
        ligand_inventory,
        preferred_ligands=preferred_ligands,
        allow_glycan=summary_row.get("ligand_context_class") == "glycan_only",
        preferred_chain=representative_chain
    )
    conn.close()
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


@app.route('/api/protacability/protein_detail')
def protacability_protein_detail():
    conn = connect_db_row()
    if not protacability_tables_available(conn):
        conn.close()
        return jsonify({
            "data_available": False,
            "message": "PROTACability data has not been imported yet. Run tools/import_protacability_data.py after generating the expansion outputs."
        }), 404

    collapse_labels = _protacability_collapse_labels(request.args.get("collapse_labels"))
    virus_name = (request.args.get("virus_name") or "").strip()
    protein_type = (request.args.get("protein_type") or "").strip()

    readiness_rows, warhead_rows = _load_protacability_enrichment_tables(conn)
    decorated_rows = _decorate_protacability_rows(
        _load_protacability_assessment_rows(conn),
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
        conn.close()
        return jsonify({"error": "Protein summary not found"}), 404

    protein_rows = _group_protein_rows(decorated_rows)
    protein_row = protein_rows[0]
    structure_rows = _group_structure_rows(decorated_rows)
    structure_rows, _ = _sort_protacability_rows(structure_rows, "summary", "best_score_desc")
    tier_distribution = dict(Counter(row.get("best_tier") or "Unknown" for row in structure_rows))
    conn.close()
    return jsonify({
        "data_available": True,
        "summary": protein_row,
        "top_structures": structure_rows[:10],
        "tier_distribution": tier_distribution,
        "explanation": "This view groups multiple structures and chains into a single protein-level summary so repeated biological contexts do not dominate the table.",
    })


@app.route('/api/protacability/target_detail')
def protacability_target_detail():
    conn = connect_db_row()
    if not protacability_tables_available(conn):
        conn.close()
        return jsonify({"data_available": False, "message": "PROTACability data is not available."}), 404

    collapse_labels = _protacability_collapse_labels(request.args.get("collapse_labels"))
    virus_name = (request.args.get("virus_name") or "").strip()
    protein_type = (request.args.get("protein_type") or "").strip()
    ligand_context_class = (request.args.get("ligand_context_class") or "").strip()
    min_score = request.args.get("min_score", type=float)

    if not virus_name or not protein_type:
        conn.close()
        return jsonify({"error": "virus_name and protein_type are required"}), 400

    readiness_rows, warhead_rows = _load_protacability_enrichment_tables(conn)
    rows = _decorate_protacability_rows(
        _load_protacability_assessment_rows(conn),
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
        conn.close()
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
        ligand_rows = [
            dict(row) for row in conn.execute(
                """
                SELECT ligand_resname, ligand_chain, ligand_residue_id, ligand_atom_count
                FROM protacability_ligand_inventory
                WHERE pdb_code = ?
                ORDER BY ligand_atom_count DESC, ligand_chain, ligand_residue_id
                """,
                (best.get("pdb_code"),)
            ).fetchall()
        ]
        preferred = _split_candidate_ligands(best.get("candidate_ligand_resnames_full"))
        ligand_record = _pick_representative_ligand_record(
            ligand_rows,
            preferred_ligands=preferred,
            allow_glycan=context_class == "glycan_only",
            preferred_chain=best.get("representative_chain_id")
        )
        ligand_resname = (ligand_record or {}).get("ligand_resname")
        return {
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
        ligand_inventory = [
            dict(row) for row in conn.execute(
                """
                SELECT
                    ligand_resname,
                    ligand_chain,
                    ligand_residue_id,
                    ligand_atom_count,
                    NULL AS ligand_heavy_atom_count
                FROM protacability_ligand_inventory
                WHERE pdb_code = ?
                ORDER BY ligand_resname, ligand_chain, ligand_residue_id
                """,
                (active_pdb_code,)
            ).fetchall()
        ]
    conn.close()
    return jsonify({
        "data_available": True,
        "target_summary": target_summary,
        "structure_rows": structure_rows,
        "ligand_groups": ligand_groups,
        "representative_contexts": representative_contexts,
        "representative_ligand": representative_ligand,
        "ligand_contexts": _serialize_ligand_contexts(ligand_inventory),
    })


@app.route('/api/protacability/export')
def protacability_export():
    conn = connect_db_row()
    if not protacability_tables_available(conn):
        conn.close()
        return jsonify({
            "success": False,
            "message": "PROTACability data has not been imported yet. Run tools/import_protacability_data.py after generating the expansion outputs."
        }), 404

    raw_export = (request.args.get("raw_export") or "").strip()
    if raw_export:
        query = data_set_queries.get(raw_export)
        if not query:
            conn.close()
            return jsonify({"success": False, "message": "Unknown PROTACability export selection."}), 400
        table_name = {
            "PROTACability Assessment": "protacability_assessment",
            "PROTACability Lysine Proximity": "protacability_lysine_proximity",
            "PROTACability Ligand Inventory": "protacability_ligand_inventory",
            "PROTACability Warhead Linkability": "protacability_warhead_linkability",
            "PROTACability Degrader Readiness": "protacability_degrader_readiness",
        }.get(raw_export)
        if table_name and not _table_exists(conn, table_name):
            conn.close()
            return jsonify({"success": False, "message": f"{raw_export} has not been imported yet."}), 404
        df = pd.read_sql_query(f"SELECT * FROM {table_name}", conn)
        conn.close()
        csv_buffer = io.StringIO()
        df.to_csv(csv_buffer, index=False)
        byte_buffer = BytesIO(csv_buffer.getvalue().encode("utf-8"))
        byte_buffer.seek(0)
        download_name = f"{table_name}.csv"
        return send_file(byte_buffer, mimetype="text/csv", as_attachment=True, download_name=download_name)

    payload = _prepare_protacability_result_set(conn, request.args, export_all=True)
    conn.close()
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


@app.route('/api/coordinates/<pdb_code>.pdb')
def serve_coordinate_for_viewer(pdb_code):
    ligand_code = str(request.args.get("ligand_code", "") or "").strip().upper()
    chain = str(request.args.get("chain", "") or "").strip().upper()
    residue_id = str(request.args.get("residue_id", "") or "").strip().upper()
    protein_only = str(request.args.get("protein_only", "0") or "0").strip() == "1"
    pdb_upper = str(pdb_code or "").strip().upper()
    if not re.match(r"^[A-Z0-9]{4}$", pdb_upper):
        return jsonify({"error": "Invalid PDB code"}), 400
    diagnostics = _resolve_coordinate_pdb(
        pdb_upper,
        ligand_code="" if protein_only else ligand_code,
        chain=chain,
        residue_id=residue_id
    )
    diagnostics["served_pdb_url"] = request.url
    app.logger.info("[coordinates] requested pdb=%s ligand=%s chain=%s residue=%s", pdb_upper, ligand_code or "(none)", chain or "(none)", residue_id or "(none)")
    app.logger.info("[coordinates] selected source=%s", diagnostics.get("selected_source"))
    app.logger.info("[coordinates] selected source format=%s", diagnostics.get("selected_source_format"))
    app.logger.info("[coordinates] converted_to_pdb=%s", diagnostics.get("converted_to_pdb") or "none")
    app.logger.info("[coordinates] FINAL_SERVED_FILE=%s", diagnostics.get("served_pdb_path"))
    app.logger.info("[coordinates] has_ATOM_records=%s", diagnostics.get("has_protein_atoms"))
    app.logger.info("[coordinates] requested_ligand=%s", ligand_code or "(none)")
    app.logger.info("[coordinates] requested_chain=%s", chain or "(none)")
    app.logger.info("[coordinates] requested_residue=%s", residue_id or "(none)")
    app.logger.info("[coordinates] contains_requested_ligand_hetatm=%s", diagnostics.get("contains_ligand"))
    app.logger.info("[coordinates] requested_ligand_hetatm_count=%s", diagnostics.get("requested_ligand_hetatm_count", 0))
    app.logger.info("[coordinates] hetatm_summary=%s", diagnostics.get("hetatm_summary"))

    served_pdb_path = diagnostics.get("served_pdb_path")
    if not served_pdb_path or not os.path.isfile(served_pdb_path):
        return jsonify({
            "error": "Unable to resolve a valid PDB coordinate source",
            "viewer_served_format": "pdb",
            "diagnostics": diagnostics
        }), 404
    if (not protein_only) and ligand_code and not diagnostics.get("contains_ligand"):
        return jsonify({
            "error": "No coordinate source contains the requested ligand context in served PDB",
            "viewer_served_format": "pdb",
            "diagnostics": diagnostics
        }), 404
    final_path = served_pdb_path
    if protein_only:
        source_text = open(served_pdb_path, "r", encoding="utf-8", errors="ignore").read()
        protein_only_path = _cache_pdb_path(pdb_upper, "PROTEINONLY")
        with open(protein_only_path, "w", encoding="utf-8") as handle:
            handle.write(_protein_only_pdb_text(source_text))
        final_path = protein_only_path
    response = send_file(final_path, mimetype="chemical/x-pdb", as_attachment=False)
    response.headers["Cache-Control"] = "no-store"
    return response


@app.route('/api/debug/coordinate_ligand_presence/<pdb_code>/<ligand_code>')
def debug_coordinate_ligand_presence(pdb_code, ligand_code):
    pdb_upper = str(pdb_code or "").strip().upper()
    ligand_upper = str(ligand_code or "").strip().upper()
    chain = str(request.args.get("chain", "") or "").strip().upper()
    residue_id = str(request.args.get("residue_id", "") or "").strip()
    aliases = make_ligand_code_aliases(ligand_upper)

    local_pdb_candidates = _candidate_coordinate_files(pdb_upper, extensions=[".pdb"])
    local_cif_candidates = _candidate_coordinate_files(pdb_upper, extensions=[".cif", ".mmcif"])
    coordinate_sources = []
    for path in (local_pdb_candidates + local_cif_candidates):
        entry = {"source": path, "exists": os.path.isfile(path)}
        if not entry["exists"]:
            coordinate_sources.append(entry)
            continue
        try:
            with open(path, "r", encoding="utf-8", errors="ignore") as handle:
                text = handle.read()
            parsed = _scan_coordinate_source(path, text, aliases, chain=chain, residue_id=residue_id)
            parsed["source_format"] = "cif" if path.lower().endswith((".cif", ".mmcif")) else "pdb"
            coordinate_sources.append(parsed)
        except Exception as exc:
            entry["error"] = str(exc)
            coordinate_sources.append(entry)

    resolved = _resolve_coordinate_pdb(
        pdb_upper,
        ligand_code=ligand_upper,
        chain=chain,
        residue_id=residue_id
    )
    served_url = url_for("serve_coordinate_for_viewer", pdb_code=pdb_upper, ligand_code=ligand_upper, chain=chain, residue_id=residue_id)
    resolved["served_pdb_url"] = served_url

    return jsonify({
        "pdb_code": pdb_upper,
        "requested_ligand": ligand_upper,
        "requested_aliases": aliases,
        "chain": chain,
        "residue_id": residue_id,
        "viewer_served_format": "pdb",
        "searched_files": local_pdb_candidates + local_cif_candidates,
        "local_pdb_candidates": local_pdb_candidates,
        "local_cif_candidates": local_cif_candidates,
        "coordinate_sources": coordinate_sources,
        "selected_source": resolved.get("selected_source"),
        "selected_source_format": resolved.get("selected_source_format"),
        "converted_to_pdb": resolved.get("converted_to_pdb"),
        "final_served_pdb_path": resolved.get("served_pdb_path"),
        "final_served_pdb_url": served_url,
        "final_pdb_has_protein_atoms": resolved.get("has_protein_atoms"),
        "final_pdb_contains_ligand": resolved.get("contains_ligand"),
        "final_pdb_hetatm_summary": resolved.get("hetatm_summary"),
        "resolution_diagnostics": resolved,
    })


@app.route('/api/debug/served_coordinate/<pdb_code>/<ligand_code>')
def debug_served_coordinate(pdb_code, ligand_code):
    pdb_upper = str(pdb_code or "").strip().upper()
    ligand_upper = str(ligand_code or "").strip().upper()
    chain = str(request.args.get("chain", "") or "").strip().upper()
    residue_id = str(request.args.get("residue_id", "") or "").strip().upper()

    resolved = _resolve_coordinate_pdb(
        pdb_upper,
        ligand_code=ligand_upper,
        chain=chain,
        residue_id=residue_id
    )
    served_path = resolved.get("served_pdb_path")
    served_url = url_for(
        "serve_coordinate_for_viewer",
        pdb_code=pdb_upper,
        ligand_code=ligand_upper,
        chain=chain,
        residue_id=residue_id
    )
    if not served_path or not os.path.isfile(served_path):
        return jsonify({
            "error": "No served coordinate file could be resolved",
            "resolution_diagnostics": resolved,
            "final_served_url": served_url
        }), 404

    text = open(served_path, "r", encoding="utf-8", errors="ignore").read()
    requested_ligand_present = pdb_contains_ligand(text, ligand_upper, chain=chain, residue_id=residue_id)
    requested_ligand_atom_count = pdb_ligand_hetatm_count(text, ligand_upper, chain=chain, residue_id=residue_id)
    return jsonify({
        "final_served_file": served_path,
        "final_served_url": served_url,
        "source_file_before_conversion": resolved.get("selected_source"),
        "source_format": resolved.get("selected_source_format"),
        "converted": bool(resolved.get("converted_to_pdb")),
        "converted_to_pdb": resolved.get("converted_to_pdb"),
        "has_atom_records": pdb_has_protein_atoms(text),
        "requested_ligand_present": requested_ligand_present,
        "requested_ligand_atom_count": requested_ligand_atom_count,
        "hetatm_summary": summarize_pdb_hetatm(text),
        "first_20_hetatm_lines": pdb_first_hetatm_lines(text, limit=20),
        "matching_ligand_lines": pdb_matching_ligand_lines(text, ligand_upper, chain=chain, residue_id=residue_id, limit=50),
        "resolution_diagnostics": resolved
    })


@app.route('/api/ligand_instance_sdf_url/<pdb_code>/<ligand_code>')
def ligand_instance_sdf_url(pdb_code, ligand_code):
    auth_chain = str(request.args.get("auth_chain", "") or "").strip().upper()
    auth_seq_id = str(request.args.get("auth_seq_id", "") or "").strip()
    mapping = _resolve_ligand_instance_mapping(
        pdb_code,
        ligand_code,
        auth_chain=auth_chain,
        auth_seq_id=auth_seq_id
    )
    if not mapping.get("success"):
        return jsonify(mapping), 404
    return jsonify({
        "success": True,
        "pdb_code": mapping.get("pdb_code"),
        "ligand_code": mapping.get("ligand_code"),
        "auth_chain": mapping.get("auth_chain"),
        "auth_seq_id": mapping.get("auth_seq_id"),
        "label_asym_id": mapping.get("chosen_label_asym_id"),
        "sdf_url": mapping.get("sdf_url"),
        "source": "rcsb_model_server",
        "mapping_source": mapping.get("source"),
    })


@app.route('/api/debug/ligand_instance_mapping/<pdb_code>/<ligand_code>')
def debug_ligand_instance_mapping(pdb_code, ligand_code):
    auth_chain = str(request.args.get("auth_chain", "") or "").strip().upper()
    auth_seq_id = str(request.args.get("auth_seq_id", "") or "").strip()
    mapping = _resolve_ligand_instance_mapping(
        pdb_code,
        ligand_code,
        auth_chain=auth_chain,
        auth_seq_id=auth_seq_id
    )
    return jsonify(mapping), (200 if mapping.get("success") else 404)


@app.route('/api/ligand_instance_sdf/<pdb_code>/<ligand_code>.sdf')
def ligand_instance_sdf_proxy(pdb_code, ligand_code):
    auth_chain = str(request.args.get("auth_chain", "") or "").strip().upper()
    auth_seq_id = str(request.args.get("auth_seq_id", "") or "").strip()
    mapping = _resolve_ligand_instance_mapping(
        pdb_code,
        ligand_code,
        auth_chain=auth_chain,
        auth_seq_id=auth_seq_id
    )
    app.logger.info(
        "[ligand-sdf] requested pdb=%s ligand=%s auth_chain=%s auth_seq_id=%s",
        str(pdb_code or "").strip().upper(),
        str(ligand_code or "").strip().upper(),
        auth_chain or "(none)",
        auth_seq_id or "(none)",
    )
    if not mapping.get("success") or not mapping.get("sdf_url"):
        return jsonify({
            "success": False,
            "error": "Unable to resolve ligand instance SDF URL",
            "mapping": mapping
        }), 404

    label_asym_id = mapping.get("chosen_label_asym_id")
    sdf_url = mapping.get("sdf_url")
    app.logger.info("[ligand-sdf] resolved label_asym_id=%s", label_asym_id)
    app.logger.info("[ligand-sdf] sdf_url=%s", sdf_url)

    try:
        sdf_bytes = _fetch_url_bytes(sdf_url)
        sdf_text = sdf_bytes.decode("utf-8", errors="ignore")
        atom_count = _count_sdf_atoms(sdf_text)
        app.logger.info("[ligand-sdf] fetched bytes=%s", len(sdf_bytes))
        app.logger.info("[ligand-sdf] atom_count=%s", atom_count)

        cache_path = _cache_sdf_path(pdb_code, ligand_code, auth_chain, auth_seq_id, label_asym_id)
        with open(cache_path, "wb") as handle:
            handle.write(sdf_bytes)
    except Exception as exc:
        return jsonify({
            "success": False,
            "error": f"Failed to fetch ligand SDF: {exc}",
            "mapping": mapping
        }), 502

    response = send_file(cache_path, mimetype="chemical/x-mdl-sdfile", as_attachment=False)
    response.headers["Cache-Control"] = "no-store"
    return response


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
        "RUPLEY_SASA_DATA",
        "Arpeggio_Contacts_Data",
        "protacability_assessment",
        "protacability_ligand_inventory",
        "protacability_lysine_proximity",
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
                ligand_resname,
                ligand_chain,
                ligand_residue_id,
                ligand_atom_count,
                NULL AS ligand_heavy_atom_count,
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
