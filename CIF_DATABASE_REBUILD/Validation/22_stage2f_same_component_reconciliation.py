"""Read-only paired graph reconciliation for Stage 2F."""
from __future__ import annotations
import argparse,csv,json,sqlite3
from collections import defaultdict
from pathlib import Path
from importlib import import_module
from rdkit import Chem
c=import_module('00_common'); m=import_module('07_map_cif_atoms_to_smiles'); ROOT=Path(__file__).resolve().parent;OUT=ROOT/'outputs'
SUCCESS={'complete','complete_altloc_resolved','partial_valid_deposited_difference'}
def rows(p):
 with open(p,newline='') as f:return list(csv.DictReader(f))
def save(name,data):
 data=list(data);fields=list(data[0]) if data else []
 with open(OUT/name,'w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(data)
def atoms(db,iid): return db.execute('SELECT * FROM ligand_instance_atoms WHERE ligand_instance_id=? AND selected_conformer=1 ORDER BY ligand_instance_atom_id',(iid,)).fetchall()
def meta(db,iid): return db.execute('''SELECT i.*,s.entry_id,s.source_cif_path,l.component_id,l.smiles FROM ligand_instances i JOIN structures s ON s.structure_id=i.structure_id JOIN ligands l ON l.ligand_id=i.ligand_id WHERE i.ligand_instance_id=?''',(iid,)).fetchone()
def graph(mol):
 a=[{'idx':x.GetIdx(),'name':x.GetProp('cif_auth_atom_id') if x.HasProp('cif_auth_atom_id') else str(x.GetIdx()),'element':x.GetSymbol(),'charge':x.GetFormalCharge(),'degree':x.GetDegree(),'aromatic':x.GetIsAromatic()} for x in mol.GetAtoms()]
 names={x['idx']:x['name'] for x in a}
 b=[(*sorted((names[x.GetBeginAtomIdx()],names[x.GetEndAtomIdx()])),str(x.GetBondType()),x.GetIsAromatic()) for x in mol.GetBonds()]
 return a,b
def main(database,checkpoint):
 db=sqlite3.connect(f'file:{Path(database).resolve()}?mode=ro',uri=True);db.row_factory=sqlite3.Row;allrows=rows(checkpoint);bycomp=defaultdict(list)
 for r in allrows:bycomp[r['component_id']].append(r)
 limits=[r for r in allrows if r['mapping_status']=='partial_adapter_limitation' and any(x['mapping_status'] in SUCCESS for x in bycomp[r['component_id']])]
 paired=[];atomcmp=[];graphcmp=[];notes=[]
 for lim in limits:
  iid=int(lim['ligand_instance_id']);la=atoms(db,iid);lm=meta(db,iid); candidates=[x for x in bycomp[lim['component_id']] if x['mapping_status'] in SUCCESS]
  def distance(x):
   sa=atoms(db,int(x['ligand_instance_id']));return abs(len(sa)-len(la))+abs(sum(bool(a['altloc']) for a in sa)-sum(bool(a['altloc']) for a in la))
  suc=min(candidates,key=distance);sid=int(suc['ligand_instance_id']);sa=atoms(db,sid);sm=meta(db,sid)
  paired.append({'component_id':lim['component_id'],'limitation_instance_id':iid,'limitation_pdb':lm['entry_id'],'success_instance_id':sid,'success_pdb':sm['entry_id'],'limitation_status':lim['mapping_status'],'success_status':suc['mapping_status'],'selected_atom_count':len(la),'selected_heavy_atom_count':sum(m.is_heavy_element(a['element']) for a in la),'success_selected_heavy_atom_count':sum(m.is_heavy_element(a['element']) for a in sa),'altloc_count':sum(bool(a['altloc']) for a in la),'success_altloc_count':sum(bool(a['altloc']) for a in sa),'occupancy_min':min((a['occupancy'] or 0) for a in la),'success_occupancy_min':min((a['occupancy'] or 0) for a in sa)})
  key=lambda a:(a['auth_atom_id'] or a['label_atom_id'],a['element'])
  lk={key(a):a for a in la};sk={key(a):a for a in sa};shared=set(lk)&set(sk)
  atomcmp.append({'component_id':lim['component_id'],'limitation_instance_id':iid,'success_instance_id':sid,'atoms_present_in_both':json.dumps(sorted(map(str,shared))),'atoms_only_in_success':json.dumps(sorted(map(str,set(sk)-set(lk)))),'atoms_only_in_limitation':json.dumps(sorted(map(str,set(lk)-set(sk)))),'element_mismatches':json.dumps([]),'atom_name_mismatches':json.dumps([])})
  lo=m.cif_occurrence_mol(la,lm['source_cif_path'],lm['label_comp_id'],m.atom_identity(lm));so=m.cif_occurrence_mol(sa,sm['source_cif_path'],sm['label_comp_id'],m.atom_identity(sm));lga,lgb=graph(lo);sga,sgb=graph(so)
  graphcmp.append({'component_id':lim['component_id'],'limitation_instance_id':iid,'success_instance_id':sid,'limitation_atoms_json':json.dumps(lga),'success_atoms_json':json.dumps(sga),'limitation_bonds_json':json.dumps(lgb),'success_bonds_json':json.dumps(sgb),'same_graph':lgb==sgb and [(x['name'],x['element']) for x in lga]==[(x['name'],x['element']) for x in sga],'missing_or_extra_bonds':json.dumps(sorted(map(str,set(lgb)^set(sgb))) )})
  notes.append(f"- {lim['component_id']} {iid} vs {sid}: atoms shared={len(shared)}/{len(lk)}; bonds differ={len(set(lgb)^set(sgb))}.")
 save('STAGE2F_PAIRED_OCCURRENCES.csv',paired);save('STAGE2F_ATOM_SET_COMPARISON.csv',atomcmp);save('STAGE2F_STRUCTURAL_GRAPH_COMPARISON.csv',graphcmp)
 (OUT/'STAGE2F_STRUCTURAL_GRAPH_DIFFERENCES.md').write_text('# Stage 2F structural graph differences\n\n'+'\n'.join(notes)+'\n')
 w3=[r for r in allrows if r['component_id']=='W3J'];(OUT/'STAGE2F_W3J_FAILURE_ANALYSIS.md').write_text('# W3J failure analysis\n\n'+json.dumps(w3,indent=2)+'\n')
 (OUT/'STAGE2F_MAPPING_RECONCILIATION_DECISION.md').write_text(f'# Stage 2F decision\n\nPaired limitation occurrences analysed: {len(paired)}. No production change was made. Inspect graph-difference output before any adapter fix. Bulk mapping is not approved.\n')
 print(f'stage2f pairs={len(paired)} components={len(set(x["component_id"] for x in paired))}',flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--database',default=str(ROOT/'viral_data_cif_v2.db'));p.add_argument('--checkpoint',default=str(OUT/'full_mapping_preflight_results.csv'));a=p.parse_args();main(a.database,a.checkpoint)
