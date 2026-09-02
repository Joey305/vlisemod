"""Controlled internal timing profile for one timeout representative/component."""
from __future__ import annotations
import argparse,csv,multiprocessing,sqlite3,time
from pathlib import Path
from importlib import import_module
from rdkit import Chem
from rdkit.Chem import AllChem,rdFMCS
c=import_module('00_common');m=import_module('07_map_cif_atoms_to_smiles');ROOT=Path(__file__).resolve().parent;OUT=ROOT/'outputs'
def worker(database,iid,q):
 try:
  db=sqlite3.connect(f'file:{Path(database).resolve()}?mode=ro',uri=True);db.row_factory=sqlite3.Row;r= m.meta if False else None
  row=db.execute('''SELECT i.*,l.smiles,s.source_cif_path FROM ligand_instances i JOIN ligands l ON l.ligand_id=i.ligand_id JOIN structures s ON s.structure_id=i.structure_id WHERE i.ligand_instance_id=?''',(iid,)).fetchone();atoms=m.selected_atoms(db,iid);times={};t=time.monotonic();smiles=Chem.MolFromSmiles(row['smiles']);times['smiles_parsing']=time.monotonic()-t;q.put(('step','smiles_parsing',times['smiles_parsing']));t=time.monotonic();conf=Chem.AddHs(smiles);times['AddHs']=time.monotonic()-t;q.put(('step','AddHs',times['AddHs']));t=time.monotonic();embed=AllChem.EmbedMolecule(conf,AllChem.ETKDG());times['ETKDG']=time.monotonic()-t;q.put(('step','ETKDG',times['ETKDG']));t=time.monotonic();
  if embed==0:
   try:AllChem.UFFOptimizeMolecule(conf)
   except Exception:pass
  times['UFF']=time.monotonic()-t;q.put(('step','UFF',times['UFF']));t=time.monotonic();smiles=Chem.RemoveHs(conf) if embed==0 else smiles;times['RemoveHs']=time.monotonic()-t;q.put(('step','RemoveHs',times['RemoveHs']));t=time.monotonic();struct=m.cif_occurrence_mol(atoms,row['source_cif_path'],row['label_comp_id'],m.atom_identity(row));times['structural_graph']=time.monotonic()-t;q.put(('step','structural_graph',times['structural_graph']));t=time.monotonic();x=rdFMCS.FindMCS([struct,smiles],completeRingsOnly=True,bondCompare=rdFMCS.BondCompare.CompareOrder);times['FindMCS']=time.monotonic()-t;q.put(('ok',times,x.numAtoms))
 except Exception as e:q.put(('error',{},f'{type(e).__name__}: {e}'))
def main(database,checkpoint,outer):
 rows=list(csv.DictReader(open(checkpoint)));seen={};
 for r in rows:
  if r['mapping_status']=='mapping_timeout':seen.setdefault(r['component_id'],r)
 out=[]
 for comp,r in seen.items():
  q=multiprocessing.get_context('spawn').Queue();p=multiprocessing.get_context('spawn').Process(target=worker,args=(database,int(r['ligand_instance_id']),q));p.start();started=time.monotonic();steps={};final=None
  while p.is_alive() and time.monotonic()-started<outer:
   try:z=q.get(timeout=1)
   except Exception:continue
   if z[0]=='step':steps[z[1]]=z[2]
   else:final=z;break
  if p.is_alive():p.terminate();p.join();out.append({'component_id':comp,'ligand_instance_id':r['ligand_instance_id'],'outcome':'outer_timeout','completed_steps':str(steps),'mcs_atoms':''})
  else:
   while not final and not q.empty():
    z=q.get();final=z if z[0]!='step' else final
   z=final or ('error',{},'worker exit');out.append({'component_id':comp,'ligand_instance_id':r['ligand_instance_id'],'outcome':z[0],'completed_steps':str(z[1]),'mcs_atoms':z[2] if len(z)>2 and isinstance(z[2],int) else ''})
 with open(OUT/'STAGE2F_TIMEOUT_PROFILE.csv','w',newline='') as f:w=csv.DictWriter(f,fieldnames=out[0]);w.writeheader();w.writerows(out)
 (OUT/'STAGE2F_TIMEOUT_PROFILE.md').write_text('# Timeout profile\n\n'+str(out)+'\n');print(f'timeout profiles={len(out)}',flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--database',default=str(ROOT/'viral_data_cif_v2.db'));p.add_argument('--checkpoint',default=str(OUT/'full_mapping_preflight_results.csv'));p.add_argument('--outer-timeout',type=float,default=180);a=p.parse_args();main(a.database,a.checkpoint,a.outer_timeout)
