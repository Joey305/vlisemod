"""Shared CIF-native foundation functions.  Never parses legacy PDB columns."""
from __future__ import annotations
import csv, hashlib, json, os, platform, sqlite3, sys, traceback
from collections import Counter, defaultdict
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable
import gemmi

ROOT=Path(__file__).resolve().parent
DEFAULT_SOURCE=(ROOT.parent/'PDB_FILES').resolve()
PIPELINE_VERSION='0.1.0'; PARSER_VERSION=getattr(gemmi,'__version__','unknown')
MISSING={'.','?',''}
# Preserved comparison policy from FIRSTROUND_Python/2_PDBSORT_2.py.
EXCLUDED=set('03N 1GP 1PE 12P 2CO 2DA 2FX 2GL 2DX 3DR 3PA 41H 4BA 43X 4DX 73O 90A ABA ACA ACE ACY ACT ADN ADP ADE AIB AMP ANP ARF ARS ASN ATM AZT ATP B3E B3P BAM BDG BDF BEN BGC BIF BMA BME BOG BR BTE BU3 BYZ CAC CAF CAS CA CD CIT CL CME CMO CMM CS CSD CSO CO CO3 CU DAL DAS DBU DDG DHI DIV DIQ DLE DLY DMS DCY DGL DGN DPN DOC DOD DIL DPV DTD DTR DTP DPP DTT DTY DVA DPR DTV EDO EPE ESD FAD FLC FE FG7 FMT FRU FUC FUM G46 G47 GGL GAL G3P GLO GLC GLY GOA GVE GOL GLU HAI HEM HEP HEZ HG HOH HPH IDG IIL IMD IOD IPA IVA K KCX KF2 MAN LAC LEU LYS MES MEA MG MLA MN M7G MK8 MNK MPD MPT MRG MSE MYR NAD NA NAG NAO NEN NH2 NH4 NI NLE NHE NO3 NTB OAS OIL OMC OXY PB PC PCA PEG PGE PG2 PG3 PG4 P6G PH2 PUT PO4 P03 P6S PPI PTR PYZ QNC RIB SIA SLZ SNU SMC SO3 SO2 SO4 STA SRT TAR TAM TEO TPO TRS TYS U2X UB4 UZ1 UZ4 UZ7 UMP URE UNX VLM VME XCP XPC YCM ZN Z9N'.split())
WATER={'HOH','DOD','WAT'}

def now(): return datetime.now(timezone.utc).isoformat()
def norm(v): return '' if v is None or str(v).strip() in MISSING else str(v).strip()
def raw(v): return '' if v is None else str(v).strip()
def sha256(path):
 h=hashlib.sha256()
 with open(path,'rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
 return h.hexdigest()
def dirs():
 for n in ('manifests','outputs','logs','temp','sql','tests'):(ROOT/n).mkdir(parents=True,exist_ok=True)
def source_paths(source, limit=None, pdb_id=None):
 ps=sorted(Path(source).rglob('*.cif'))
 if pdb_id: ps=[p for p in ps if p.stem.upper()==pdb_id.upper()]
 return ps[:limit] if limit else ps
def classification(path,source):
 rel=path.relative_to(source); bits=rel.parts
 return (bits[0] if bits else '', bits[1] if len(bits)>2 else '', str(rel))
def cif_doc(path): return gemmi.cif.read_file(str(path))
def getv(block, tag, default=''):
 try: return raw(block.find_value(tag)) or default
 except Exception:return default
def loop_rows(block, tags):
 column=block.find_loop(tags[0])
 if not column: return []
 # Loop columns can be physically ordered differently from requested tags.
 # Read each tagged column explicitly rather than zipping raw loop storage.
 columns=[list(block.find_loop(tag)) for tag in tags]
 return [dict(zip(tags,values)) for values in zip(*columns)]
ATOM_TAGS=['_atom_site.group_PDB','_atom_site.id','_atom_site.pdbx_PDB_model_num','_atom_site.label_asym_id','_atom_site.auth_asym_id','_atom_site.label_comp_id','_atom_site.auth_comp_id','_atom_site.label_seq_id','_atom_site.auth_seq_id','_atom_site.pdbx_PDB_ins_code','_atom_site.label_atom_id','_atom_site.auth_atom_id','_atom_site.type_symbol','_atom_site.label_alt_id','_atom_site.occupancy','_atom_site.B_iso_or_equiv','_atom_site.Cartn_x','_atom_site.Cartn_y','_atom_site.Cartn_z']
def atoms(path):
 b=cif_doc(path).sole_block(); return getv(b,'_entry.id',path.stem), loop_rows(b,ATOM_TAGS)
def fnum(v):
 try:return float(v)
 except:return None
def curation(comp):
 c=norm(comp).upper()
 if c in WATER:return 'excluded','water'
 if c in EXCLUDED:return 'excluded','legacy_non_ligand_exclusion'
 return 'included','included_ligand_candidate'
@contextmanager
def dbconn(db):
 c=sqlite3.connect(db); c.row_factory=sqlite3.Row; c.execute('PRAGMA foreign_keys=ON')
 try: yield c; c.commit()
 except: c.rollback(); raise
 finally:c.close()
SCHEMA='''
CREATE TABLE IF NOT EXISTS structures (structure_id INTEGER PRIMARY KEY, entry_id TEXT NOT NULL, source_cif_path TEXT NOT NULL, source_cif_sha256 TEXT NOT NULL, file_size INTEGER NOT NULL, source_status TEXT, parser TEXT NOT NULL, parser_version TEXT NOT NULL, ingest_time TEXT NOT NULL, UNIQUE(entry_id,source_cif_sha256));
CREATE TABLE IF NOT EXISTS structure_classifications (classification_id INTEGER PRIMARY KEY, structure_id INTEGER NOT NULL REFERENCES structures(structure_id), virus_label TEXT NOT NULL, protein_label TEXT NOT NULL, source_relative_path TEXT NOT NULL, UNIQUE(structure_id,source_relative_path));
CREATE TABLE IF NOT EXISTS structure_context (context_id INTEGER PRIMARY KEY, pdb_id TEXT NOT NULL, virus_raw TEXT NOT NULL, virus_normalized TEXT, target_context_raw TEXT NOT NULL, target_context_normalized TEXT, source_cif_path TEXT NOT NULL UNIQUE, classification_source TEXT NOT NULL, classification_qc TEXT NOT NULL);
CREATE TABLE IF NOT EXISTS mapping_remediation_queue (remediation_id INTEGER PRIMARY KEY, ligand_instance_id INTEGER NOT NULL UNIQUE REFERENCES ligand_instances(ligand_instance_id), pdb_id TEXT NOT NULL, component_id TEXT NOT NULL, virus_normalized TEXT, target_context_normalized TEXT, mapping_reason_class TEXT NOT NULL, algorithm_queue TEXT NOT NULL, preflight_mapping_status TEXT NOT NULL, heavy_atoms_structural INTEGER, heavy_atoms_reference INTEGER, heavy_atoms_mapped INTEGER, heavy_atom_mapping_fraction REAL, remediation_status TEXT NOT NULL, first_identified_run_id INTEGER REFERENCES analysis_runs(run_id), latest_attempt_run_id INTEGER REFERENCES analysis_runs(run_id), notes TEXT NOT NULL DEFAULT '');
CREATE TABLE IF NOT EXISTS ligands (ligand_id INTEGER PRIMARY KEY, component_id TEXT NOT NULL, smiles TEXT, smiles_source TEXT, source_version TEXT, canonical_smiles TEXT, inchikey TEXT, chemical_definition_version TEXT NOT NULL DEFAULT 'local-unresolved-v1', chemical_status TEXT NOT NULL DEFAULT 'missing', UNIQUE(component_id,chemical_definition_version));
CREATE TABLE IF NOT EXISTS ligand_instances (ligand_instance_id INTEGER PRIMARY KEY, structure_id INTEGER NOT NULL REFERENCES structures(structure_id), ligand_id INTEGER NOT NULL REFERENCES ligands(ligand_id), deposited_model_num TEXT NOT NULL, parser_model_index INTEGER, label_asym_id TEXT NOT NULL, label_comp_id TEXT NOT NULL, auth_asym_id TEXT NOT NULL, auth_comp_id TEXT, label_seq_id TEXT, auth_seq_id TEXT NOT NULL, insertion_code_raw TEXT, insertion_code_normalized TEXT NOT NULL, identity_status TEXT NOT NULL, curation_status TEXT NOT NULL, curation_reason TEXT NOT NULL, curation_rule_version TEXT NOT NULL, UNIQUE(structure_id,deposited_model_num,label_asym_id,label_comp_id,auth_asym_id,auth_seq_id,insertion_code_normalized));
CREATE TABLE IF NOT EXISTS ligand_instance_atoms (ligand_instance_atom_id INTEGER PRIMARY KEY, ligand_instance_id INTEGER NOT NULL REFERENCES ligand_instances(ligand_instance_id), atom_site_id TEXT NOT NULL, label_atom_id TEXT, auth_atom_id TEXT, element TEXT, x REAL,y REAL,z REAL,altloc TEXT,occupancy REAL,b_factor REAL,deposited_model_num TEXT NOT NULL,selected_conformer INTEGER NOT NULL DEFAULT 0,conformer_selection_reason TEXT, UNIQUE(ligand_instance_id,atom_site_id));
CREATE TABLE IF NOT EXISTS analysis_runs (run_id INTEGER PRIMARY KEY,stage TEXT NOT NULL,method TEXT NOT NULL,method_version TEXT NOT NULL,pipeline_version TEXT NOT NULL,parameters_json TEXT NOT NULL,python_version TEXT NOT NULL,package_versions TEXT NOT NULL,start_time TEXT NOT NULL,end_time TEXT,source_manifest_sha256 TEXT,status TEXT NOT NULL,processed_count INTEGER DEFAULT 0,success_count INTEGER DEFAULT 0,partial_count INTEGER DEFAULT 0,failure_count INTEGER DEFAULT 0);
CREATE TABLE IF NOT EXISTS pipeline_failures (failure_id INTEGER PRIMARY KEY,run_id INTEGER REFERENCES analysis_runs(run_id),structure_id INTEGER REFERENCES structures(structure_id),ligand_instance_id INTEGER REFERENCES ligand_instances(ligand_instance_id),stage TEXT NOT NULL,reason_code TEXT NOT NULL,message TEXT NOT NULL,traceback_log_path TEXT,created_at TEXT NOT NULL);
CREATE TABLE IF NOT EXISTS ligand_sasa_atoms (sasa_id INTEGER PRIMARY KEY,run_id INTEGER NOT NULL REFERENCES analysis_runs(run_id),ligand_instance_id INTEGER NOT NULL REFERENCES ligand_instances(ligand_instance_id),ligand_instance_atom_id INTEGER NOT NULL REFERENCES ligand_instance_atoms(ligand_instance_atom_id),sasa_area REAL NOT NULL,legacy_exposed INTEGER NOT NULL,probe_radius REAL NOT NULL,point_density INTEGER NOT NULL,water_treatment TEXT NOT NULL,conformer_policy TEXT NOT NULL,status TEXT NOT NULL,UNIQUE(run_id,ligand_instance_atom_id));
CREATE TABLE IF NOT EXISTS ligand_smiles_atom_mapping (mapping_id INTEGER PRIMARY KEY,run_id INTEGER NOT NULL REFERENCES analysis_runs(run_id),ligand_instance_id INTEGER NOT NULL REFERENCES ligand_instances(ligand_instance_id),ligand_instance_atom_id INTEGER REFERENCES ligand_instance_atoms(ligand_instance_atom_id),smiles_atom_index INTEGER,mapping_status TEXT NOT NULL,method_version TEXT NOT NULL,error_reason TEXT,UNIQUE(run_id,ligand_instance_atom_id,smiles_atom_index));
CREATE TABLE IF NOT EXISTS ligand_chemistry_sources (chemistry_source_id INTEGER PRIMARY KEY,ligand_id INTEGER NOT NULL REFERENCES ligands(ligand_id),source_name TEXT NOT NULL,source_smiles TEXT,source_locator TEXT,source_priority INTEGER NOT NULL,parse_status TEXT NOT NULL,UNIQUE(ligand_id,source_name,source_smiles));
CREATE TABLE IF NOT EXISTS ligand_mapping_runs (mapping_run_id INTEGER PRIMARY KEY,run_id INTEGER NOT NULL UNIQUE REFERENCES analysis_runs(run_id),ligand_instance_id INTEGER NOT NULL REFERENCES ligand_instances(ligand_instance_id),mapping_status TEXT NOT NULL,mapped_count INTEGER NOT NULL,structural_atom_count INTEGER NOT NULL,smiles_atom_count INTEGER NOT NULL,mcs_atom_count INTEGER NOT NULL,unmatched_structural_atom_ids_json TEXT NOT NULL,unmatched_smiles_atom_indices_json TEXT NOT NULL,ambiguity_reason TEXT,method_version TEXT NOT NULL);
CREATE TABLE IF NOT EXISTS ligand_arpeggio_runs (arpeggio_run_id INTEGER PRIMARY KEY,run_id INTEGER NOT NULL REFERENCES analysis_runs(run_id),ligand_instance_id INTEGER NOT NULL REFERENCES ligand_instances(ligand_instance_id),source_cif_sha256 TEXT NOT NULL,selector TEXT,selector_namespace TEXT,command TEXT,arpeggio_version TEXT,start_time TEXT,end_time TEXT,exit_status INTEGER,stdout_path TEXT,stderr_path TEXT,output_sha256 TEXT,status TEXT NOT NULL);
CREATE TABLE IF NOT EXISTS arpeggio_attempts (attempt_id INTEGER PRIMARY KEY,run_id INTEGER NOT NULL REFERENCES analysis_runs(run_id),ligand_instance_id INTEGER NOT NULL REFERENCES ligand_instances(ligand_instance_id),attempt_number INTEGER NOT NULL,input_strategy TEXT NOT NULL,input_path TEXT NOT NULL,input_sha256 TEXT NOT NULL,is_fallback INTEGER NOT NULL,fallback_reason TEXT,sanitization_operations_json TEXT NOT NULL,timeout_seconds REAL NOT NULL,start_time TEXT NOT NULL,end_time TEXT,runtime_seconds REAL,command TEXT NOT NULL,return_code INTEGER,status TEXT NOT NULL,failure_class TEXT,error_class TEXT,error_message TEXT,stdout_path TEXT NOT NULL,stderr_path TEXT NOT NULL,output_path TEXT,output_sha256 TEXT,output_validation_status TEXT,provenance_validation_status TEXT,UNIQUE(run_id,attempt_number));
CREATE TABLE IF NOT EXISTS arpeggio_derived_atom_map (derived_atom_map_id INTEGER PRIMARY KEY,run_id INTEGER NOT NULL REFERENCES analysis_runs(run_id),input_strategy TEXT NOT NULL,derived_atom_id TEXT NOT NULL,source_atom_site_id TEXT NOT NULL,ligand_instance_atom_id INTEGER REFERENCES ligand_instance_atoms(ligand_instance_atom_id),source_model_num TEXT NOT NULL,source_label_asym_id TEXT,source_auth_asym_id TEXT,source_label_seq_id TEXT,source_auth_seq_id TEXT,source_component_id TEXT,source_label_atom_id TEXT,source_auth_atom_id TEXT,source_element TEXT,source_altloc TEXT,source_insertion_code TEXT,derived_model_num TEXT NOT NULL,derived_label_asym_id TEXT,derived_auth_asym_id TEXT,derived_label_seq_id TEXT,derived_auth_seq_id TEXT,derived_component_id TEXT,derived_label_atom_id TEXT,derived_auth_atom_id TEXT,derived_element TEXT,derived_altloc TEXT,derived_insertion_code TEXT,mapping_status TEXT NOT NULL,UNIQUE(run_id,input_strategy,derived_atom_id));
CREATE TABLE IF NOT EXISTS arpeggio_failure_classifications (classification_id INTEGER PRIMARY KEY,ligand_instance_id INTEGER NOT NULL REFERENCES ligand_instances(ligand_instance_id),source_run_id INTEGER NOT NULL REFERENCES analysis_runs(run_id),source_status TEXT NOT NULL,failure_class TEXT NOT NULL,classifier_version TEXT NOT NULL,classified_at TEXT NOT NULL,UNIQUE(source_run_id));
CREATE TABLE IF NOT EXISTS arpeggio_raw_contact_labels (raw_contact_id INTEGER PRIMARY KEY,run_id INTEGER NOT NULL REFERENCES analysis_runs(run_id),ligand_instance_id INTEGER NOT NULL REFERENCES ligand_instances(ligand_instance_id),raw_contact_index INTEGER NOT NULL,interaction_label TEXT NOT NULL,distance REAL,bgn_json TEXT NOT NULL,end_json TEXT NOT NULL,ligand_instance_atom_id INTEGER REFERENCES ligand_instance_atoms(ligand_instance_atom_id),partner_identity_json TEXT NOT NULL,filter_class TEXT NOT NULL,UNIQUE(run_id,raw_contact_index,interaction_label));
CREATE TABLE IF NOT EXISTS arpeggio_unique_atom_pairs (pair_id INTEGER PRIMARY KEY,run_id INTEGER NOT NULL REFERENCES analysis_runs(run_id),ligand_instance_id INTEGER NOT NULL REFERENCES ligand_instances(ligand_instance_id),ligand_instance_atom_id INTEGER REFERENCES ligand_instance_atoms(ligand_instance_atom_id),partner_identity_json TEXT NOT NULL,raw_label_count INTEGER NOT NULL,UNIQUE(run_id,ligand_instance_atom_id,partner_identity_json));
CREATE TABLE IF NOT EXISTS receptor_binding_pocket_atoms (pocket_atom_id INTEGER PRIMARY KEY,run_id INTEGER NOT NULL REFERENCES analysis_runs(run_id),ligand_instance_id INTEGER NOT NULL REFERENCES ligand_instances(ligand_instance_id),partner_atom_site_id TEXT,partner_label_asym_id TEXT,partner_auth_asym_id TEXT,partner_auth_seq_id TEXT,distance REAL NOT NULL);
'''
def create_schema(db):
 with dbconn(db) as c:
  c.executescript(SCHEMA)
  # Formal additive migration for databases created by the foundation milestone.
  for table,column,kind in [('ligand_smiles_atom_mapping','mapped_count','INTEGER'),('ligand_smiles_atom_mapping','structural_atom_count','INTEGER'),('ligand_smiles_atom_mapping','smiles_atom_count','INTEGER'),('ligand_smiles_atom_mapping','mcs_atom_count','INTEGER'),('ligand_smiles_atom_mapping','unmatched_structural_atom_ids_json','TEXT'),('ligand_smiles_atom_mapping','unmatched_smiles_atom_indices_json','TEXT'),('ligand_smiles_atom_mapping','mapping_method','TEXT'),('ligand_mapping_runs','mapping_outcome','TEXT'),('ligand_mapping_runs','mapping_reason_class','TEXT'),('ligand_mapping_runs','algorithm_queue','TEXT'),('ligand_mapping_runs','heavy_atoms_structural','INTEGER'),('ligand_mapping_runs','heavy_atoms_reference','INTEGER'),('ligand_mapping_runs','heavy_atoms_mapped','INTEGER'),('ligand_mapping_runs','heavy_atom_mapping_fraction','REAL'),('ligand_mapping_runs','mapping_complete','INTEGER'),('ligand_mapping_runs','downstream_mapping_eligibility','INTEGER'),('ligand_arpeggio_runs','derived_input_path','TEXT'),('ligand_arpeggio_runs','derived_input_sha256','TEXT'),('ligand_arpeggio_runs','canonical_deposited_model_num','TEXT'),('ligand_arpeggio_runs','input_strategy','TEXT'),('ligand_arpeggio_runs','original_attempt_status','TEXT'),('ligand_arpeggio_runs','fallback_attempted','INTEGER NOT NULL DEFAULT 0'),('ligand_arpeggio_runs','fallback_reason','TEXT'),('ligand_arpeggio_runs','sanitization_operations_json','TEXT'),('ligand_arpeggio_runs','attempt_count','INTEGER NOT NULL DEFAULT 0'),('ligand_arpeggio_runs','timeout_seconds','REAL'),('ligand_arpeggio_runs','runtime_seconds','REAL'),('ligand_arpeggio_runs','error_class','TEXT'),('ligand_arpeggio_runs','error_message','TEXT'),('ligand_arpeggio_runs','provenance_validation_status','TEXT'),('ligand_arpeggio_runs','output_validation_status','TEXT'),('ligand_arpeggio_runs','completion_mode','TEXT'),('arpeggio_raw_contact_labels','partner_source_atom_site_id','TEXT'),('arpeggio_raw_contact_labels','partner_mapping_status','TEXT'),('ligand_sasa_atoms','deposited_model_num','TEXT'),('ligand_sasa_atoms','method_version','TEXT')]:
   if column not in {r[1] for r in c.execute(f'PRAGMA table_info({table})')}:
    c.execute(f'ALTER TABLE {table} ADD COLUMN {column} {kind}')
def run_start(c,stage,params,manifest=''):
 return c.execute("INSERT INTO analysis_runs(stage,method,method_version,pipeline_version,parameters_json,python_version,package_versions,start_time,source_manifest_sha256,status) VALUES(?,?,?,?,?,?,?,?,?,?)",(stage,stage,'v2-foundation-1',PIPELINE_VERSION,json.dumps(params,sort_keys=True),sys.version.split()[0],json.dumps({'gemmi':PARSER_VERSION}),now(),manifest,'running')).lastrowid
def run_end(c,rid,status,processed=0,success=0,partial=0,failures=0): c.execute('UPDATE analysis_runs SET end_time=?,status=?,processed_count=?,success_count=?,partial_count=?,failure_count=? WHERE run_id=?',(now(),status,processed,success,partial,failures,rid))
def fail(c,rid,stage,msg,structure_id=None,instance_id=None,code='exception'):
 c.execute('INSERT INTO pipeline_failures(run_id,structure_id,ligand_instance_id,stage,reason_code,message,created_at) VALUES(?,?,?,?,?,?,?)',(rid,structure_id,instance_id,stage,code,msg,now()))
def set_conformer_flags(c, instance_id):
 rows=c.execute('SELECT ligand_instance_atom_id,altloc,occupancy,element FROM ligand_instance_atoms WHERE ligand_instance_id=?',(instance_id,)).fetchall()
 named=defaultdict(list); shared=[]
 for r in rows:
  a=norm(r['altloc']); (shared if not a else named[a]).append(r)
 if not named:
  c.execute('UPDATE ligand_instance_atoms SET selected_conformer=1,conformer_selection_reason=? WHERE ligand_instance_id=?',('shared_atoms_only',instance_id)); return
 # Heavy-atom completeness then summed occupancy then lexical tie break.
 scores=[]
 for alt, rs in named.items(): scores.append((sum(1 for r in rs if norm(r['element']).upper()!='H'),sum(r['occupancy'] or 0 for r in rs),alt))
 best=sorted(scores,key=lambda x:(-x[0],-x[1],x[2]))[0]; reason=f'coherent_altloc={best[2]};heavy_atoms={best[0]};occupancy_sum={best[1]:.3f}'
 ids=[r['ligand_instance_atom_id'] for r in shared+named[best[2]]]
 c.execute('UPDATE ligand_instance_atoms SET selected_conformer=0,conformer_selection_reason=? WHERE ligand_instance_id=?',('rejected_'+reason,instance_id))
 c.executemany('UPDATE ligand_instance_atoms SET selected_conformer=1,conformer_selection_reason=? WHERE ligand_instance_atom_id=?',[(reason,i) for i in ids])
