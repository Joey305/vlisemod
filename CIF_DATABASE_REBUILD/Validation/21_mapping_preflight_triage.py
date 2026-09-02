"""Read-only Stage 2E evidence audit for a completed mapping preflight."""
from __future__ import annotations
import argparse,csv,json,sqlite3
from collections import Counter,defaultdict
from pathlib import Path
from importlib import import_module
c=import_module('00_common'); mapping=import_module('07_map_cif_atoms_to_smiles')
ROOT=Path(__file__).resolve().parent; OUT=ROOT/'outputs'
STATUSES=['complete','complete_altloc_resolved','partial_valid_deposited_difference','partial_missing_deposited_atoms','partial_adapter_limitation','ambiguous','mapping_timeout','failed','other']
def read_csv(path):
 with path.open(newline='') as f:return list(csv.DictReader(f))
def write_csv(path,rows):
 rows=list(rows);fields=list(rows[0]) if rows else []
 with path.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
def metals(atoms): return sorted({(a['element'] or '').upper() for a in atoms if (a['element'] or '').upper() not in {'H','D','T','C','N','O','P','S','F','CL','BR','I','B','SI','SE'}})
def classify(r,atoms):
 if metals(atoms): return 'METAL_OR_UNSUPPORTED_ELEMENT'
 if r.get('mapped_heavy_count',0)==r.get('structural_heavy_count',-1)==r.get('smiles_heavy_count',-2): return 'STATUS_CLASSIFICATION_BUG'
 if r.get('structural_heavy_count',0)!=r.get('smiles_heavy_count',0): return 'CCD_SMILES_ATOM_COUNT_MISMATCH'
 return 'MCS_RING_OR_BOND_ORDER_LIMITATION'
def main(database,checkpoint,timeout):
 rows=read_csv(checkpoint);counts=Counter(r['mapping_status'] for r in rows);db=sqlite3.connect(f'file:{Path(database).resolve()}?mode=ro',uri=True);db.row_factory=sqlite3.Row
 audit=[{'mapping_status':s,'count':counts.get(s,0)} for s in STATUSES];write_csv(OUT/'FULL_MAPPING_PREFLIGHT_STATUS.csv',audit)
 successes=defaultdict(set)
 for x in rows:
  if x['mapping_status'] in {'complete','complete_altloc_resolved','partial_valid_deposited_difference'}:successes[x['component_id']].add(x['mapping_status'])
 adapter=[x for x in rows if x['mapping_status']=='partial_adapter_limitation']; detail=[];comparisons=[]
 for x in adapter:
  iid=int(x['ligand_instance_id']);atoms=db.execute('SELECT * FROM ligand_instance_atoms WHERE ligand_instance_id=? AND selected_conformer=1 ORDER BY ligand_instance_atom_id',(iid,)).fetchall();result=mapping.map_one_isolated(database,iid,timeout)
  if result['outcome']!='result':
   r={'structural_count':len(atoms),'structural_heavy_count':sum(mapping.is_heavy_element(a['element']) for a in atoms),'smiles_count':0,'smiles_heavy_count':0,'mapped_heavy_count':0,'mcs':0,'pairs':[],'unmatched':[],'unmatched_s':[],'reason':result['outcome']}
  else:r=result['result']
  unmatched=set(r.get('unmatched',[]));unstruct=[a['auth_atom_id'] or a['label_atom_id'] or a['atom_site_id'] for a in atoms if a['atom_site_id'] in unmatched and mapping.is_heavy_element(a['element'])]
  status=classify(r,atoms);methods=Counter(m for _,_,m in r.get('pairs',[]) if m)
  detail.append({'ligand_instance_id':iid,'pdb_id':x['pdb_id'],'component_id':x['component_id'],'model':db.execute('SELECT deposited_model_num FROM ligand_instances WHERE ligand_instance_id=?',(iid,)).fetchone()[0],'auth_chain':db.execute('SELECT auth_asym_id FROM ligand_instances WHERE ligand_instance_id=?',(iid,)).fetchone()[0],'auth_residue':db.execute('SELECT auth_seq_id FROM ligand_instances WHERE ligand_instance_id=?',(iid,)).fetchone()[0],'insertion_code':db.execute('SELECT insertion_code_normalized FROM ligand_instances WHERE ligand_instance_id=?',(iid,)).fetchone()[0],'selected_atom_count':r.get('structural_count',0),'selected_heavy_atom_count':r.get('structural_heavy_count',0),'smiles_heavy_atom_count':r.get('smiles_heavy_count',0),'mapped_heavy_atom_count':r.get('mapped_heavy_count',0),'heavy_atom_mapping_fraction':round(r.get('mapped_heavy_count',0)/max(1,r.get('smiles_heavy_count',0)),6),'unmatched_structural_heavy_atoms':json.dumps(unstruct),'unmatched_smiles_heavy_atoms':json.dumps([i for i in r.get('unmatched_s',[]) ]),'initial_mcs_count':r.get('mcs',0),'post_oxygen_recovery_count':len(r.get('pairs',[])),'mapping_provenance':json.dumps(methods,sort_keys=True),'reason_code':r.get('reason',''),'root_cause':status})
  comparisons.append({'component_id':x['component_id'],'limitation_instance_id':iid,'successful_statuses':';'.join(sorted(successes[x['component_id']])),'same_component_succeeds_elsewhere':bool(successes[x['component_id']]),'selected_atom_count':r.get('structural_count',0),'selected_heavy_atom_count':r.get('structural_heavy_count',0),'altloc_count':sum(bool(a['altloc']) for a in atoms),'occupancy_min':min((a['occupancy'] or 0) for a in atoms) if atoms else 0,'model':detail[-1]['model'],'unmatched_heavy_count':len(unstruct),'root_cause':status})
 write_csv(OUT/'ADAPTER_LIMITATION_COMPONENT_COMPARISON.csv',comparisons);write_csv(OUT/'FULL_MAPPING_PREFLIGHT_ADAPTER_DETAILS.csv',detail)
 timeouts=[]
 for x in (r for r in rows if r['mapping_status']=='mapping_timeout'):
  iid=int(x['ligand_instance_id']);atoms=db.execute('SELECT element FROM ligand_instance_atoms WHERE ligand_instance_id=? AND selected_conformer=1',(iid,)).fetchall();timeouts.append({'ligand_instance_id':iid,'pdb_id':x['pdb_id'],'component_id':x['component_id'],'atom_count':len(atoms),'heavy_atom_count':sum(mapping.is_heavy_element(a['element']) for a in atoms),'smiles_atom_count':db.execute('SELECT length(l.smiles) FROM ligand_instances i JOIN ligands l ON l.ligand_id=i.ligand_id WHERE i.ligand_instance_id=?',(iid,)).fetchone()[0],'same_component_succeeds_elsewhere':bool(successes[x['component_id']]),'timed_step':'watchdog_exceeded_unlocalized','elapsed_seconds':x['elapsed_seconds']})
 failures=[]
 for x in (r for r in rows if r['mapping_status']=='failed'):
  iid=int(x['ligand_instance_id']);z=mapping.map_one_isolated(database,iid,timeout);reason=z.get('result',{}).get('reason','') if z['outcome']=='result' else z.get('exception_message','');failures.append({'ligand_instance_id':iid,'pdb_id':x['pdb_id'],'component_id':x['component_id'],'exception_type':z.get('exception_type',''),'exception_message':reason,'stage_of_failure':'MCS_no_mapping' if reason=='failed' else 'adapter'})
 write_csv(OUT/'MAPPING_TIMEOUT_TRIAGE.csv',timeouts);write_csv(OUT/'MAPPING_FAILURE_TRIAGE.csv',failures)
 root=Counter(d['root_cause'] for d in detail);components=Counter(d['component_id'] for d in detail)
 md=[f'# Full Mapping Preflight Triage\n\nTotal checkpoint rows: {len(rows)}.', '## Status audit', *[f'- {s}: {counts.get(s,0)}' for s in STATUSES],f'\nAdapter limitations: {len(detail)} across {len(components)} component IDs.','## Root causes',*[f'- {k}: {v}' for k,v in root.most_common()], '\n## Timeouts',f'- {len(timeouts)} timeouts; component concentration: {Counter(x["component_id"] for x in timeouts).most_common()}. Each is recorded as `watchdog_exceeded_unlocalized`; no timeout limit was increased.', '\n## Failures',f'- {len(failures)} failures; components: {Counter(x["component_id"] for x in failures)}.','\n## Decision','- No production mapping was run or modified.','- No safe global adapter fix is approved from this classification-only pass.','- Bulk mapping is **not approved**: partial_adapter_limitation is not zero.']
 (OUT/'FULL_MAPPING_PREFLIGHT_TRIAGE.md').write_text('\n'.join(md)+'\n');(OUT/'STAGE2E_FULL_MAPPING_TRIAGE_DECISION.md').write_text('\n'.join(md)+'\n')
 print(f'triage complete: adapter_limitations={len(detail)} unique_components={len(components)} timeouts={len(timeouts)} failures={len(failures)}',flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--database',default=str(ROOT/'viral_data_cif_v2.db'));p.add_argument('--checkpoint',default=str(OUT/'full_mapping_preflight_results.csv'));p.add_argument('--per-instance-timeout',type=float,default=60);a=p.parse_args();main(a.database,Path(a.checkpoint),a.per_instance_timeout)
