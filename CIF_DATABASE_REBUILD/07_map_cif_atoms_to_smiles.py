"""CIF-instance-to-SMILES mapping; structural identity is always ligand_instance_id."""
from __future__ import annotations
import argparse,json,multiprocessing,queue,sqlite3,time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from importlib import import_module
from rdkit import Chem
from rdkit.Chem import AllChem,rdFMCS
c=import_module('00_common'); VERSION='legacy_mcs_etkdg_uff_cif_v2.5'
def selected_atoms(db,iid): return db.execute('SELECT * FROM ligand_instance_atoms WHERE ligand_instance_id=? AND selected_conformer=1 ORDER BY ligand_instance_atom_id',(iid,)).fetchall()
def rdkit_element(value):
 """Return the RDKit element and isotope without rewriting source mmCIF data."""
 token=(value or 'C').strip()
 # D and T are valid mmCIF type symbols, but RDKit represents them as
 # hydrogen isotopes rather than periodic-table elements.
 if token.upper()=='D': return ('H',2)
 if token.upper()=='T': return ('H',3)
 return (token[:1].upper()+token[1:].lower() if token else 'C',None)
def is_heavy_element(value): return (value or '').strip().upper() not in {'H','D','T'}
def atom_identity(row):
 """Stable deposited atom identity; never compare this to a legacy PDB serial."""
 return {'model':row['deposited_model_num'],'auth_chain':row['auth_asym_id'],'auth_seq_id':row['auth_seq_id'],'insertion_code':row['insertion_code_normalized'],'comp_id':row['label_comp_id']}
def cif_occurrence_mol(rows, source_cif, comp_id, occurrence):
 """Build one RDKit graph from exact selected atom rows and CIF CCD bonds."""
 block=c.cif_doc(Path(source_cif)).sole_block(); tags=['_chem_comp_bond.comp_id','_chem_comp_bond.atom_id_1','_chem_comp_bond.atom_id_2','_chem_comp_bond.value_order']
 bonds=c.loop_rows(block,tags); rw=Chem.RWMol(); by_name={}; order={'SING':Chem.BondType.SINGLE,'DOUB':Chem.BondType.DOUBLE,'TRIP':Chem.BondType.TRIPLE,'AROM':Chem.BondType.AROMATIC}
 for row in rows:
  symbol,isotope=rdkit_element(row['element']); a=Chem.Atom(symbol)
  if isotope is not None: a.SetIsotope(isotope)
  idx=rw.AddAtom(a); name=row['auth_atom_id'] or row['label_atom_id']; by_name[name]=idx
  for k,v in {'cif_atom_site_id':row['atom_site_id'],'cif_auth_atom_id':row['auth_atom_id'] or '','cif_label_atom_id':row['label_atom_id'] or '','cif_alt_id':row['altloc'] or '','cif_type_symbol_original':row['element'] or '','cif_occurrence_key':json.dumps(occurrence,sort_keys=True)}.items(): rw.GetAtomWithIdx(idx).SetProp(k,str(v))
 for b in bonds:
  if c.norm(b['_chem_comp_bond.comp_id'])!=comp_id:continue
  n1,n2=c.norm(b['_chem_comp_bond.atom_id_1']),c.norm(b['_chem_comp_bond.atom_id_2'])
  if n1 in by_name and n2 in by_name: rw.AddBond(by_name[n1],by_name[n2],order.get(c.norm(b['_chem_comp_bond.value_order'])[:4],Chem.BondType.SINGLE))
 mol=rw.GetMol()
 try: Chem.SanitizeMol(mol)
 except Exception: pass # Some deposited graphs have non-standard valence.
 # Substructure matching requires ring information even if full sanitization is
 # unavailable for a deposited/non-standard component.
 Chem.FastFindRings(mol)
 return mol
def pdb_block(rows):
 lines=[]
 for n,r in enumerate(rows,1):
  serial=int(r['atom_site_id']) if str(r['atom_site_id']).isdigit() and int(r['atom_site_id'])<100000 else n
  name=(r['auth_atom_id'] or r['label_atom_id'] or r['element'] or 'X')[:4].rjust(4); element=(r['element'] or 'C')[:2].rjust(2)
  lines.append(f'HETATM{serial:5d} {name} LIG A   1    {r["x"]:8.3f}{r["y"]:8.3f}{r["z"]:8.3f}{(r["occupancy"] if r["occupancy"] is not None else 1):6.2f}{(r["b_factor"] if r["b_factor"] is not None else 0):6.2f}          {element}')
 return '\n'.join(lines)+'\nEND\n'
def oxygen_recovery(struct,smiles,pairs):
 # pairs are (structural-RDKit-index, SMILES-RDKit-index). Keep those
 # namespaces separate; the previous reversed bookkeeping disabled the
 # historical carbonyl-oxygen recovery for 3EKY/DR7.
 mapped_struct={s for s,_ in pairs}; mapped_smiles={t for _,t in pairs}; extra=[]
 for si in set(range(struct.GetNumAtoms()))-mapped_struct:
  a=struct.GetAtomWithIdx(si)
  if a.GetSymbol()!='O':continue
  for nb in a.GetNeighbors():
   hit=next((t for s,t in pairs if s==nb.GetIdx()),None)
   if hit is None:continue
   for b in smiles.GetAtomWithIdx(hit).GetBonds():
    o=b.GetOtherAtom(smiles.GetAtomWithIdx(hit))
    if o.GetSymbol()=='O' and b.GetBondType()==Chem.BondType.DOUBLE and o.GetIdx() not in mapped_smiles: extra.append((si,o.GetIdx()));mapped_struct.add(si);mapped_smiles.add(o.GetIdx());break
 return extra
def map_one(db,iid):
 row=db.execute('''SELECT i.*,l.smiles,l.chemical_status,s.source_cif_path FROM ligand_instances i JOIN ligands l ON l.ligand_id=i.ligand_id JOIN structures s ON s.structure_id=i.structure_id WHERE i.ligand_instance_id=?''',(iid,)).fetchone()
 if not row: raise ValueError(f'unknown ligand_instance_id {iid}')
 atoms=selected_atoms(db,iid); structural_count=len(atoms); smi=row['smiles']
 if row['chemical_status']!='resolved' or not smi:return {'iid':iid,'status':'failed','reason':'missing_or_unresolved_smiles','atoms':atoms,'pairs':[],'structural_count':structural_count,'smiles_count':0,'mcs':0}
 smiles=Chem.MolFromSmiles(smi)
 if not smiles:return {'iid':iid,'status':'failed','reason':'invalid_smiles','atoms':atoms,'pairs':[],'structural_count':structural_count,'smiles_count':0,'mcs':0}
 # The historical pipeline generated a 3-D conformer before matching.  Its
 # MCS operation, however, is entirely graph-based.  Some valid lipid-like
 # CCD SMILES cannot be embedded by ETKDG; retain their validated 2-D graph
 # rather than turning a coordinate-generation limitation into a false
 # atom-mapping failure.  The fallback is retained in run provenance.
 preparation='etkdg_uff'
 conformer=Chem.AddHs(smiles)
 if AllChem.EmbedMolecule(conformer,AllChem.ETKDG())==0:
  try: AllChem.UFFOptimizeMolecule(conformer)
  except Exception: preparation='etkdg_uff_failed_2d_fallback'
  else: smiles=Chem.RemoveHs(conformer)
 else: preparation='etkdg_failed_2d_fallback'
 structural=cif_occurrence_mol(atoms,row['source_cif_path'],row['label_comp_id'],atom_identity(row))
 if not structural:return {'iid':iid,'status':'failed','reason':'occurrence_pdb_adapter_failed','atoms':atoms,'pairs':[],'structural_count':structural_count,'smiles_count':smiles.GetNumAtoms(),'mcs':0}
 try:
  m=rdFMCS.FindMCS([structural,smiles],completeRingsOnly=True,bondCompare=rdFMCS.BondCompare.CompareOrder)
  patt=Chem.MolFromSmarts(m.smartsString); a=structural.GetSubstructMatch(patt) if patt else (); b=smiles.GetSubstructMatch(patt) if patt else ()
 except Exception as exc:
  return {'iid':iid,'status':'failed','reason':f'rdkit_graph_exception:{type(exc).__name__}','atoms':atoms,'pairs':[],'structural_count':structural_count,'smiles_count':smiles.GetNumAtoms(),'mcs':0,'unmatched':[x['atom_site_id'] for x in atoms],'unmatched_s':list(range(smiles.GetNumAtoms()))}
 mcs_pairs=list(zip(a,b)); recovery_pairs=oxygen_recovery(structural,smiles,mcs_pairs)
 pairs=list(dict.fromkeys(mcs_pairs+recovery_pairs)); recovery_set=set(recovery_pairs)
 # RDKit indices are transient: provenance returns every mapping to its exact
 # occurrence-scoped database atom, without any PDB serial assumption.
 resolved=[(atoms[si],mi) for si,mi in pairs if si < len(atoms)]
 matched={x['ligand_instance_atom_id'] for x,_ in resolved}; unmatched=[x['atom_site_id'] for x in atoms if x['ligand_instance_atom_id'] not in matched]; unmatched_s=sorted(set(range(smiles.GetNumAtoms()))-{mi for _,mi in resolved})
 # Hydrogen isotopes are never heavy atoms.  Keep their original source token
 # and all-deposited-atom provenance, but assess chemical completeness on the
 # same heavy-atom basis used by the functional-group mapping.
 structural_heavy=sum(is_heavy_element(x['element']) for x in atoms); smiles_heavy=sum(a.GetAtomicNum()>1 for a in smiles.GetAtoms())
 mapped_heavy=sum(is_heavy_element(atom['element']) and smiles.GetAtomWithIdx(si).GetAtomicNum()>1 for atom,si in resolved)
 unmatched_heavy=[x['atom_site_id'] for x in atoms if x['ligand_instance_atom_id'] not in matched and is_heavy_element(x['element'])]
 unmatched_smiles_heavy=[si for si in unmatched_s if smiles.GetAtomWithIdx(si).GetAtomicNum()>1]
 heavy_status='complete' if not unmatched_heavy and not unmatched_smiles_heavy else ('partial_ccd_difference' if structural_heavy!=smiles_heavy else 'partial_adapter_limitation')
 if heavy_status=='complete': status='complete_altloc_resolved' if any(a['altloc'] for a in atoms) else 'complete'
 elif resolved: status='partial_ccd_difference' if structural_count!=smiles.GetNumAtoms() else 'partial_adapter_limitation'
 else: status='failed'
 reason=status if status!='complete' else (preparation if preparation!='etkdg_uff' else '')
 return {'iid':iid,'status':status,'reason':reason,'atoms':atoms,'pairs':resolved,'mapping_methods':{(atom['ligand_instance_atom_id'],si):('legacy_oxygen_recovery' if (atoms.index(atom),si) in recovery_set else 'MCS') for atom,si in resolved},'structural_count':structural_count,'smiles_count':smiles.GetNumAtoms(),'mcs':m.numAtoms,'unmatched':unmatched,'unmatched_s':unmatched_s,'structural_heavy_count':structural_heavy,'smiles_heavy_count':smiles_heavy,'mapped_heavy_count':mapped_heavy,'heavy_atom_mapping_status':heavy_status}
def exception_result(db,iid,exc):
 """Convert one unexpected adapter error into an auditable instance failure."""
 row=db.execute('''SELECT i.*,l.smiles,s.entry_id,l.component_id FROM ligand_instances i JOIN ligands l ON l.ligand_id=i.ligand_id JOIN structures s ON s.structure_id=i.structure_id WHERE i.ligand_instance_id=?''',(iid,)).fetchone()
 atoms=selected_atoms(db,iid) if row else []
 message=f'{type(exc).__name__}: {exc}'
 return {'iid':iid,'status':'failed','reason':f'instance_exception:{message}','atoms':atoms,'pairs':[],'structural_count':len(atoms),'smiles_count':0,'mcs':0,'unmatched':[x['atom_site_id'] for x in atoms],'unmatched_s':[],'component_id':row['component_id'] if row else None,'pdb_id':row['entry_id'] if row else None}
BUSY_TIMEOUT_MS=60000
DOWNSTREAM_ELIGIBLE_STATUSES={'complete','complete_altloc_resolved','partial_ccd_difference'}

def _configure_connection(db):
 """Apply bounded SQLite waiting without changing scientific data."""
 db.execute(f'PRAGMA busy_timeout={BUSY_TIMEOUT_MS}')
 return db

def _map_worker(database,iid,queue):
 """Child entry point: RDKit runs outside any parent write transaction."""
 db=None
 try:
  db=sqlite3.connect(f'file:{Path(database).resolve()}?mode=ro',uri=True,timeout=BUSY_TIMEOUT_MS/1000.0)
  db.row_factory=sqlite3.Row;_configure_connection(db);db.execute('PRAGMA query_only=ON')
  r=map_one(db,iid)
  payload={k:v for k,v in r.items() if k not in {'atoms','pairs','mapping_methods'}}
  payload['pairs']=[(atom['ligand_instance_atom_id'],si,r.get('mapping_methods',{}).get((atom['ligand_instance_atom_id'],si))) for atom,si in r['pairs']]
  queue.put(('result',payload))
 except Exception as exc:
  queue.put(('exception',type(exc).__name__,str(exc)))
 finally:
  if db is not None: db.close()

def run_killable(worker,args,timeout,context='spawn'):
 """Run one operation behind a terminable process boundary."""
 ctx=multiprocessing.get_context(context); result_queue=ctx.Queue(); proc=ctx.Process(target=worker,args=(*args,result_queue)); started=time.monotonic();proc.start();proc.join(timeout);elapsed=time.monotonic()-started
 if proc.is_alive():
  proc.terminate();proc.join()
  try: result_queue.close()
  except Exception: pass
  return {'outcome':'timeout','elapsed_seconds':elapsed}
 try: item=result_queue.get(timeout=1)
 except queue.Empty:
  try: result_queue.close()
  except Exception: pass
  return {'outcome':'exception','elapsed_seconds':elapsed,'exception_type':'WorkerExit','exception_message':f'worker exited with code {proc.exitcode}'}
 try: result_queue.close()
 except Exception: pass
 if item[0]=='exception': return {'outcome':'exception','elapsed_seconds':elapsed,'exception_type':item[1],'exception_message':item[2]}
 return {'outcome':'result','elapsed_seconds':elapsed,'result':item[1]}

def map_one_isolated(database,iid,timeout=60):
 """Map an occurrence without allowing a hung RDKit call to block the run."""
 return run_killable(_map_worker,(str(database),iid),timeout)

def _load_work(database,limit=None,pdb_id=None,instance_id=None,use_remediation_registry=True,resume=False):
 """Load the worklist and close SQLite before any worker process starts."""
 with c.dbconn(database) as db:
  _configure_connection(db)
  q="SELECT i.ligand_instance_id FROM ligand_instances i JOIN ligands l ON l.ligand_id=i.ligand_id JOIN structures s ON s.structure_id=i.structure_id WHERE i.curation_status='included' AND l.chemical_status='resolved'";args=[]
  if pdb_id:q+=' AND UPPER(s.entry_id)=UPPER(?)';args.append(pdb_id)
  if instance_id:q+=' AND i.ligand_instance_id=?';args.append(instance_id)
  if resume:
   q+=' AND NOT EXISTS (SELECT 1 FROM ligand_mapping_runs mr WHERE mr.ligand_instance_id=i.ligand_instance_id AND mr.method_version=?)';args.append(VERSION)
  ids=[r[0] for r in db.execute(q+' ORDER BY i.ligand_instance_id',args)]
  ids=ids[:limit] if limit else ids
  skips={r['ligand_instance_id']:dict(r) for r in db.execute('SELECT * FROM mapping_remediation_queue WHERE remediation_status="pending"')} if use_remediation_registry else {}
 return ids,skips

def _start_instance_run(database,iid,per_instance_timeout,use_remediation_registry):
 """Do parent bookkeeping in a short transaction, then release the lock."""
 with c.dbconn(database) as db:
  _configure_connection(db)
  db.execute('DELETE FROM ligand_smiles_atom_mapping WHERE ligand_instance_id=? AND method_version=?',(iid,VERSION))
  db.execute('DELETE FROM ligand_mapping_runs WHERE ligand_instance_id=? AND method_version=?',(iid,VERSION))
  rid=c.run_start(db,'mapping',{'ligand_instance_id':iid,'method':VERSION,'per_instance_timeout':per_instance_timeout,'remediation_registry':use_remediation_registry,'sqlite_lock_policy':'no_parent_write_transaction_during_worker;busy_timeout_ms=60000'})
 return rid

def _materialize_result(database,rid,iid,r,outcome,reason_class,algorithm_queue,failure_code=None):
 """Write one finished worker result only after the read-only child has exited."""
 eligible=1 if r['status'] in DOWNSTREAM_ELIGIBLE_STATUSES else 0
 with c.dbconn(database) as db:
  _configure_connection(db)
  if failure_code:
   c.fail(db,rid,'mapping',r['reason'],instance_id=iid,code=failure_code)
  db.execute('INSERT OR REPLACE INTO ligand_mapping_runs(run_id,ligand_instance_id,mapping_status,mapped_count,structural_atom_count,smiles_atom_count,mcs_atom_count,unmatched_structural_atom_ids_json,unmatched_smiles_atom_indices_json,ambiguity_reason,method_version,mapping_outcome,mapping_reason_class,algorithm_queue,heavy_atoms_structural,heavy_atoms_reference,heavy_atoms_mapped,heavy_atom_mapping_fraction,mapping_complete,downstream_mapping_eligibility) VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)',(rid,iid,r['status'],len(r['pairs']),r['structural_count'],r['smiles_count'],r['mcs'],json.dumps(r.get('unmatched',[])),json.dumps(r.get('unmatched_s',[])),r['reason'],VERSION,outcome,reason_class,algorithm_queue,r.get('structural_heavy_count'),r.get('smiles_heavy_count'),r.get('mapped_heavy_count'),(r.get('mapped_heavy_count',0)/max(1,r.get('smiles_heavy_count',1))) if r.get('smiles_heavy_count') is not None else None,1 if r['status'] in {'complete','complete_altloc_resolved'} else 0,eligible))
  for atom,si in r['pairs']:
   db.execute('INSERT OR IGNORE INTO ligand_smiles_atom_mapping(run_id,ligand_instance_id,ligand_instance_atom_id,smiles_atom_index,mapping_status,method_version,error_reason,mapped_count,structural_atom_count,smiles_atom_count,mcs_atom_count,unmatched_structural_atom_ids_json,unmatched_smiles_atom_indices_json,mapping_method) VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?)',(rid,iid,atom['ligand_instance_atom_id'],si,r['status'],VERSION,r['reason'],len(r['pairs']),r['structural_count'],r['smiles_count'],r['mcs'],json.dumps(r.get('unmatched',[])),json.dumps(r.get('unmatched_s',[])),r.get('mapping_methods',{}).get((atom['ligand_instance_atom_id'],si))))
  is_failure=r['status'] in {'failed','mapping_timeout'}
  c.run_end(db,rid,'failed' if is_failure else 'completed',1,1 if r['status'] in {'complete','complete_altloc_resolved'} else 0,1 if r['status'].startswith('partial') else 0,1 if is_failure else 0)

def _decode_isolated_result(database,iid,isolated):
 """Convert one child-process outcome into the standard mapping result shape."""
 failure_code=None
 if isolated['outcome']=='result':
  r=isolated['result'];pair_payload=r['pairs'];r['pairs']=[({'ligand_instance_atom_id':atom_id},si) for atom_id,si,_ in pair_payload];r['mapping_methods']={(atom_id,si):method for atom_id,si,method in pair_payload}
 elif isolated['outcome']=='timeout':
  r={'iid':iid,'status':'mapping_timeout','reason':f'mapping_timeout:{isolated["elapsed_seconds"]:.3f}s','pairs':[],'structural_count':0,'smiles_count':0,'mcs':0,'unmatched':[],'unmatched_s':[]};failure_code='mapping_timeout'
 else:
  with c.dbconn(database) as db:
   _configure_connection(db);exc=RuntimeError(f"{isolated['exception_type']}: {isolated['exception_message']}");r=exception_result(db,iid,exc)
  failure_code='instance_exception'
 outcome='timeout' if r['status']=='mapping_timeout' else ('failed' if r['status']=='failed' else ('partial' if r['status'].startswith('partial') else r['status']))
 return r,outcome,None,None,failure_code

def run(database,limit=None,pdb_id=None,instance_id=None,per_instance_timeout=60,use_remediation_registry=True,resume=False,progress_every=100,workers=4):
 """Run occurrence mapping with read-only workers separated from SQLite writes.

 Workers are launched in bounded waves.  All read-only child processes in a
 wave finish before the parent materializes any results, so SQLite never has to
 arbitrate a long-lived parent write transaction against mapping readers.
 """
 c.create_schema(database)
 ids,skips=_load_work(database,limit,pdb_id,instance_id,use_remediation_registry,resume)
 workers=max(1,int(workers or 1))
 print(f'mapping selection: pdb_id={pdb_id or "ALL"} instances={len(ids)} method={VERSION} remediation_registry={use_remediation_registry} workers={workers}',flush=True)
 results=[];complete=partial=skipped=failures=timeouts=0;processed=0
 for batch_start in range(0,len(ids),workers):
  batch_ids=ids[batch_start:batch_start+workers]
  run_ids={iid:_start_instance_run(database,iid,per_instance_timeout,use_remediation_registry) for iid in batch_ids}
  decoded={}
  active=[iid for iid in batch_ids if iid not in skips]
  # All active workers in this wave are read-only.  Do not write to SQLite
  # until the entire wave has returned or timed out.
  if active:
   with ThreadPoolExecutor(max_workers=min(workers,len(active))) as pool:
    futures={iid:pool.submit(map_one_isolated,database,iid,per_instance_timeout) for iid in active}
    isolated_results={iid:futures[iid].result() for iid in active}
   for iid in active:
    decoded[iid]=_decode_isolated_result(database,iid,isolated_results[iid])
  for iid in batch_ids:
   remediation=skips.get(iid)
   if remediation:
    r={'iid':iid,'status':'skipped_pending_remediation','reason':f"{remediation['mapping_reason_class']}:{remediation['algorithm_queue']}",'pairs':[],'structural_count':remediation['heavy_atoms_structural'] or 0,'smiles_count':remediation['heavy_atoms_reference'] or 0,'mcs':0,'unmatched':[],'unmatched_s':[],'structural_heavy_count':remediation['heavy_atoms_structural'],'smiles_heavy_count':remediation['heavy_atoms_reference'],'mapped_heavy_count':remediation['heavy_atoms_mapped']}
    outcome='skipped_pending_remediation';reason_class=remediation['mapping_reason_class'];algorithm_queue=remediation['algorithm_queue'];failure_code=None;skipped+=1
   else:
    r,outcome,reason_class,algorithm_queue,failure_code=decoded[iid]
   if r['status'] in {'complete','complete_altloc_resolved'}:complete+=1
   elif r['status'].startswith('partial'):partial+=1
   elif r['status']=='mapping_timeout':timeouts+=1
   elif r['status']=='failed':failures+=1
   _materialize_result(database,run_ids[iid],iid,r,outcome,reason_class,algorithm_queue,failure_code)
   results.append(r);processed+=1
   if processed%max(1,progress_every)==0 or processed==len(ids):
    print(f'mapping progress: {processed}/{len(ids)} complete={complete} partial={partial} skipped={skipped} timeouts={timeouts} failures={failures}',flush=True)
 return results

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--database',default=str(c.ROOT/'viral_data_cif_v2.db'));p.add_argument('--limit',type=int);p.add_argument('--pdb-id');p.add_argument('--ligand-instance-id',type=int);p.add_argument('--per-instance-timeout',type=float,default=60);p.add_argument('--ignore-remediation-registry',action='store_true');p.add_argument('--resume',action='store_true');p.add_argument('--progress-every',type=int,default=100);p.add_argument('--workers',type=int,default=4);a=p.parse_args()
 rows=run(a.database,a.limit,a.pdb_id,a.ligand_instance_id,a.per_instance_timeout,not a.ignore_remediation_registry,a.resume,a.progress_every,a.workers)
 counts={}
 for r in rows: counts[r['status']]=counts.get(r['status'],0)+1
 print(json.dumps({'method':VERSION,'processed':len(rows),'status_counts':counts},indent=2,sort_keys=True))
