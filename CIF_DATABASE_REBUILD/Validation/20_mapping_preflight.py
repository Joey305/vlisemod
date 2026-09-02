"""Resumable, read-only, process-isolated mapping preflight."""
from __future__ import annotations
import argparse,csv,sqlite3,time
from pathlib import Path
from importlib import import_module
c=import_module('00_common'); mapping=import_module('07_map_cif_atoms_to_smiles')
FIELDS=['sequence_number','ligand_instance_id','pdb_id','component_id','mapping_status','reason_code','elapsed_seconds','exception_type','exception_message']
def metadata(db,iid):
 r=db.execute('''SELECT l.component_id,s.entry_id FROM ligand_instances i JOIN ligands l ON l.ligand_id=i.ligand_id JOIN structures s ON s.structure_id=i.structure_id WHERE i.ligand_instance_id=?''',(iid,)).fetchone()
 return (r['entry_id'],r['component_id']) if r else ('','')
def display_status(status):
 # Legacy stored name retained in DB; the preflight reports its scientific
 # classification required for bulk-execution review.
 return 'partial_valid_deposited_difference' if status=='partial_ccd_difference' else status
def resume_ids(path):
 if not path.exists(): return set()
 with path.open(newline='') as f:return {r['ligand_instance_id'] for r in csv.DictReader(f) if r.get('ligand_instance_id')}
def count_status(counts,status):
 if status in {'complete','complete_altloc_resolved'}:counts['complete']+=1
 elif status=='partial_valid_deposited_difference':counts['valid_partial']+=1
 elif status=='partial_adapter_limitation':counts['adapter']+=1
 elif status=='mapping_timeout':counts['timeout']+=1
 elif status=='failed':counts['failed']+=1
 else:counts['other']+=1
def append_result(path,row):
 new=not path.exists() or path.stat().st_size==0
 with path.open('a',newline='') as f:
  w=csv.DictWriter(f,fieldnames=FIELDS)
  if new:w.writeheader()
  w.writerow(row);f.flush()
def main(database,limit,resume,progress_every,per_instance_timeout,results_path):
 db=sqlite3.connect(f'file:{Path(database).resolve()}?mode=ro',uri=True);db.row_factory=sqlite3.Row
 ids=[r[0] for r in db.execute("SELECT i.ligand_instance_id FROM ligand_instances i JOIN ligands l USING(ligand_id) WHERE i.curation_status='included' AND l.chemical_status='resolved' ORDER BY i.ligand_instance_id LIMIT ?",(limit,))]
 counts={'complete':0,'valid_partial':0,'timeout':0,'failed':0,'adapter':0,'other':0};started=time.monotonic();scanned=0
 done=resume_ids(results_path) if resume else set()
 if resume and results_path.exists():
  with results_path.open(newline='') as f:
   for old in csv.DictReader(f): count_status(counts,old['mapping_status'])
 for sequence,iid in enumerate(ids,1):
  if str(iid) in done: continue
  pdb_id,component_id=metadata(db,iid);result=mapping.map_one_isolated(database,iid,per_instance_timeout);status='failed';reason='';etype='';emessage=''
  if result['outcome']=='result': status=display_status(result['result']['status']);reason=result['result'].get('reason','')
  elif result['outcome']=='timeout': status='mapping_timeout';reason='mapping_timeout';emessage=f"{result['elapsed_seconds']:.3f}s"
  else: reason='instance_exception';etype=result['exception_type'];emessage=result['exception_message']
  elapsed=result['elapsed_seconds'];row={'sequence_number':sequence,'ligand_instance_id':iid,'pdb_id':pdb_id,'component_id':component_id,'mapping_status':status,'reason_code':reason,'elapsed_seconds':f'{elapsed:.6f}','exception_type':etype,'exception_message':emessage};append_result(results_path,row);scanned+=1
  count_status(counts,status)
  if scanned%progress_every==0 or sequence==len(ids):
   print(f'preflight progress: {sequence}/{len(ids)} ligand_instance_id={iid} pdb_id={pdb_id} component_id={component_id} elapsed_seconds={time.monotonic()-started:.1f} last_instance_seconds={elapsed:.3f} complete={counts["complete"]} valid_partial={counts["valid_partial"]} timeout={counts["timeout"]} failure={counts["failed"]}',flush=True)
 print(f'preflight scanned={len(ids)} newly_scanned={scanned} complete={counts["complete"]} valid_partial={counts["valid_partial"]} partial_adapter_limitation={counts["adapter"]} mapping_timeout={counts["timeout"]} failed={counts["failed"]} other={counts["other"]}',flush=True)
 db.close()
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--database',default=str(c.ROOT/'viral_data_cif_v2.db'));p.add_argument('--limit',type=int,default=7335);p.add_argument('--resume',action='store_true');p.add_argument('--progress-every',type=int,default=25);p.add_argument('--per-instance-timeout',type=float,default=60);p.add_argument('--results',default=str(c.ROOT/'outputs/full_mapping_preflight_results.csv'));a=p.parse_args();main(a.database,a.limit,a.resume,a.progress_every,a.per_instance_timeout,Path(a.results))
