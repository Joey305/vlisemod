from __future__ import annotations
import argparse,importlib,sys
from pathlib import Path
c=importlib.import_module('00_common')
STAGES={'inventory':'01_inventory_cifs','create':'02_create_database','ingest':'03_ingest_structures','curate':'05_curate_ligands','chemistry':'06_load_ligand_chemistry','mapping':'07_map_cif_atoms_to_smiles','sasa':'08_calculate_ligand_sasa','arpeggio':'09_run_arpeggio','geometry':'10_calculate_ligand_geometry','functional-groups':'11_assign_functional_groups','protacability':'12_build_protacability','attachment':'13_build_attachment_sites','validate':'15_validate_database','compare':'16_compare_to_legacy'}
def main():
 p=argparse.ArgumentParser();p.add_argument('--all',action='store_true');p.add_argument('--stage',choices=STAGES);p.add_argument('--source',default=str(c.DEFAULT_SOURCE));p.add_argument('--database',default=str(c.ROOT/'viral_data_cif_v2.db'));p.add_argument('--limit',type=int);p.add_argument('--pdb-id');p.add_argument('--ligand-instance-id',type=int);p.add_argument('--workers',type=int,default=1);p.add_argument('--resume',action='store_true');p.add_argument('--per-instance-timeout',type=float,default=60);p.add_argument('--retry-timeout',type=float,default=600);p.add_argument('--fallback-radius',type=float,default=12);p.add_argument('--progress-every',type=int,default=25);p.add_argument('--dry-run',action='store_true');a=p.parse_args()
 print(f'Arguments: source={a.source} database={a.database} limit={a.limit} pdb_id={a.pdb_id}', flush=True)
 todo=['inventory','create','ingest','curate','validate'] if a.all else [a.stage]
 if a.dry_run: print('Would run:',', '.join(todo));return
 for s in todo:
  print(f'Running stage: {s}', flush=True)
  if s=='inventory': importlib.import_module(STAGES[s]).inventory(Path(a.source),a.limit,a.pdb_id)
  elif s=='create': c.create_schema(a.database)
  elif s=='ingest': importlib.import_module(STAGES[s]).ingest(Path(a.source),a.database,a.limit,a.pdb_id,a.resume)
  elif s=='curate':
   with c.dbconn(a.database) as db:
    for r in db.execute('SELECT ligand_instance_id,label_comp_id FROM ligand_instances'):
     q,reason=c.curation(r['label_comp_id']);db.execute('UPDATE ligand_instances SET curation_status=?,curation_reason=? WHERE ligand_instance_id=?',(q,reason,r['ligand_instance_id']))
  elif s=='validate': print(importlib.import_module(STAGES[s]).validate(a.database))
  elif s=='chemistry': print(f"chemistry identities: {len(importlib.import_module(STAGES[s]).resolve(a.database,a.limit,a.pdb_id))}")
  elif s=='mapping': print([(x['iid'],x['status']) for x in importlib.import_module(STAGES[s]).run(a.database,a.limit,a.pdb_id,a.ligand_instance_id,a.per_instance_timeout)])
  elif s=='arpeggio': print(importlib.import_module(STAGES[s]).run(a.database,a.limit,a.pdb_id,a.ligand_instance_id,a.resume,a.per_instance_timeout,a.progress_every,a.workers,a.retry_timeout,a.fallback_radius))
  elif s=='sasa': print(importlib.import_module(STAGES[s]).run(a.database,a.limit,a.pdb_id,a.ligand_instance_id))
  else: __import__('runpy').run_path(str(c.ROOT/(STAGES[s]+'.py')),run_name='__main__')
if __name__=='__main__':
 main()
