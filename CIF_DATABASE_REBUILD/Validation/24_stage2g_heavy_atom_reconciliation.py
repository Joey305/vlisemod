"""Read-only Stage 2G heavy-atom and local-template reconciliation."""
from __future__ import annotations
import argparse,csv,json,sqlite3,time
from collections import Counter
from pathlib import Path
from importlib import import_module
from rdkit import Chem
from rdkit.Chem import rdFMCS
c=import_module('00_common');m=import_module('07_map_cif_atoms_to_smiles');ROOT=Path(__file__).resolve().parent;OUT=ROOT/'outputs'
def csvrows(p):
 with open(p,newline='') as f:return list(csv.DictReader(f))
def write(name,data):
 data=list(data);fields=list(data[0]) if data else []
 with open(OUT/name,'w',newline='') as f:w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(data)
def heavy_graph(mol):
 keep=[a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum()>1];rw=Chem.RWMol();oldnew={}
 for old in keep:
  a=Chem.Atom(mol.GetAtomWithIdx(old));oldnew[old]=rw.AddAtom(a)
 for b in mol.GetBonds():
  if b.GetBeginAtomIdx() in oldnew and b.GetEndAtomIdx() in oldnew:rw.AddBond(oldnew[b.GetBeginAtomIdx()],oldnew[b.GetEndAtomIdx()],b.GetBondType())
 result=rw.GetMol();Chem.FastFindRings(result);return result,oldnew
def named_heavy(mol):
 names={a.GetIdx():(a.GetProp('cif_auth_atom_id') if a.HasProp('cif_auth_atom_id') else str(a.GetIdx())) for a in mol.GetAtoms() if a.GetAtomicNum()>1}
 bonds={tuple(sorted((names[b.GetBeginAtomIdx()],names[b.GetEndAtomIdx()])))+(str(b.GetBondType()),) for b in mol.GetBonds() if b.GetBeginAtomIdx() in names and b.GetEndAtomIdx() in names}
 return set(names.values()),bonds
def main(database,pairsfile):
 db=sqlite3.connect(f'file:{Path(database).resolve()}?mode=ro',uri=True);db.row_factory=sqlite3.Row;pairs=csvrows(pairsfile);audit=[];before=[];templates=[]
 for p in pairs:
  li,si=int(p['limitation_instance_id']),int(p['success_instance_id']);lr=db.execute('''SELECT i.*,l.smiles,s.source_cif_path FROM ligand_instances i JOIN ligands l ON l.ligand_id=i.ligand_id JOIN structures s ON s.structure_id=i.structure_id WHERE i.ligand_instance_id=?''',(li,)).fetchone();sr=db.execute('''SELECT i.*,s.source_cif_path FROM ligand_instances i JOIN structures s ON s.structure_id=i.structure_id WHERE i.ligand_instance_id=?''',(si,)).fetchone();la,sa=m.selected_atoms(db,li),m.selected_atoms(db,si);lm=m.cif_occurrence_mol(la,lr['source_cif_path'],lr['label_comp_id'],m.atom_identity(lr));sm=m.cif_occurrence_mol(sa,sr['source_cif_path'],sr['label_comp_id'],m.atom_identity(sr));lset,lb=named_heavy(lm);sset,sb=named_heavy(sm);cur=m.map_one(db,li);hm,idx=heavy_graph(lm);ref=Chem.MolFromSmiles(lr['smiles']);t=time.monotonic();x=rdFMCS.FindMCS([hm,ref],completeRingsOnly=True,bondCompare=rdFMCS.BondCompare.CompareOrder,timeout=5);heavy_seconds=time.monotonic()-t;pat=Chem.MolFromSmarts(x.smartsString) if not x.canceled else None;a=hm.GetSubstructMatch(pat) if pat else ();b=ref.GetSubstructMatch(pat) if pat else ();hpairs=list(zip(a,b));hpairs.extend(m.oxygen_recovery(hm,ref,hpairs));inverse={new:old for old,new in idx.items()};mapped_ids=[la[inverse[hi]]['ligand_instance_atom_id'] for hi,_ in hpairs]
  audit.append({'component_id':p['component_id'],'limitation_instance_id':li,'success_instance_id':si,'heavy_atom_set_equal':lset==sset,'heavy_heavy_connectivity_equal':lb==sb,'heavy_atom_count_limitation':len(lset),'heavy_atom_count_comparator':len(sset),'heavy_bond_count_limitation':len(lb),'heavy_bond_count_comparator':len(sb),'mapped_heavy_atom_count':cur.get('mapped_heavy_count',0),'reference_smiles_heavy_atom_count':cur.get('smiles_heavy_count',0),'heavy_atom_fraction_mapped':cur.get('mapped_heavy_count',0)/max(1,cur.get('smiles_heavy_count',0)),'unmatched_structural_heavy_atoms':json.dumps([v for v in cur.get('unmatched',[])]),'unmatched_smiles_heavy_atoms':json.dumps(cur.get('unmatched_s',[])),'mcs_heavy_atom_count':x.numAtoms,'mapping_status':cur['status'],'reason_code':cur['reason']})
  before.append({'component_id':p['component_id'],'limitation_instance_id':li,'current_mcs_atoms':cur['mcs'],'current_mapped_heavy_atoms':cur.get('mapped_heavy_count',0),'heavy_mcs_atoms':x.numAtoms,'heavy_mapped_atoms':len(hpairs),'heavy_mapped_ligand_instance_atom_ids':json.dumps(mapped_ids),'current_status':cur['status'],'heavy_runtime_seconds':round(heavy_seconds,6),'method':'controlled_heavy_atom_mcs_plus_legacy_oxygen_recovery'})
 targets={x['component_id'] for x in pairs}|{'QNG','9VI','146','B0F','FSC','GTG','XIO','W3J'}
 for comp in sorted(targets):
  row=db.execute('''SELECT i.ligand_instance_id,s.source_cif_path,l.smiles FROM ligand_instances i JOIN ligands l ON l.ligand_id=i.ligand_id JOIN structures s ON s.structure_id=i.structure_id WHERE l.component_id=? LIMIT 1''',(comp,)).fetchone();bond_rows=0
  if row:
   block=c.cif_doc(Path(row['source_cif_path'])).sole_block();bond_rows=sum(1 for x in c.loop_rows(block,['_chem_comp_bond.comp_id','_chem_comp_bond.atom_id_1','_chem_comp_bond.atom_id_2','_chem_comp_bond.value_order']) if c.norm(x['_chem_comp_bond.comp_id'])==comp)
  legacy=sqlite3.connect(ROOT.parent/'viral_data.db').execute('SELECT count(*) FROM Ligand_Atoms_Smiles WHERE ligand=?',(comp,)).fetchone()[0]
  templates.append({'component_id':comp,'mmcif_component_bond_rows':bond_rows,'local_smiles_available':bool(row and row['smiles']),'legacy_occurrence_smiles_rows':legacy,'atom_name_aware_reference_smiles_template_available':False,'template_path_decision':'not_safe_to_implement_without_authoritative_CIF_atom_name_to_SMILES_index_mapping'})
 write('STAGE2G_HEAVY_ATOM_PAIR_AUDIT.csv',audit);write('STAGE2G_HEAVY_MCS_BEFORE_AFTER.csv',before);write('STAGE2G_COMPONENT_TEMPLATE_AUDIT.csv',templates)
 eq=sum(str(x['heavy_atom_set_equal'])=='True' and str(x['heavy_heavy_connectivity_equal'])=='True' for x in audit);status='Current mapper determines partial_adapter_limitation from heavy-atom MCS incompleteness when total atom counts match; explicit H is no longer counted as heavy, so the 99 hydrogen-only graph observation does not by itself justify reclassification.'
 (OUT/'STAGE2G_HEAVY_ATOM_PAIR_AUDIT.md').write_text(f'# Heavy atom audit\n\nIdentical named heavy topology: {eq}/{len(audit)}.\n')
 (OUT/'STAGE2G_STATUS_CLASSIFICATION_AUDIT.md').write_text('# Status classification audit\n\n'+status+'\n')
 (OUT/'STAGE2G_COMPONENT_TEMPLATE_AUDIT.md').write_text('# Component template audit\n\nmmCIF `_chem_comp_bond` definitions are present where reported, but no authoritative local CIF atom-name → SMILES-index template was found. An exact-template correspondence is therefore not safe to add yet.\n')
 (OUT/'STAGE2G_MAPPING_DECISION.md').write_text(f'# Stage 2G decision\n\nHeavy-topology-identical paired cases: {eq}/{len(audit)}. No production code changed. Exact-template path is not approved because atom-name-aware reference correspondence is unavailable locally. FindMCS remains the confirmed timeout step. Bulk mapping is not approved; a new full preflight is only warranted after a separately proven mapping/status fix.\n')
 print(f'heavy audit pairs={len(audit)} identical_heavy_topology={eq}',flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--database',default=str(ROOT/'viral_data_cif_v2.db'));p.add_argument('--pairs',default=str(OUT/'STAGE2F_PAIRED_OCCURRENCES.csv'));a=p.parse_args();main(a.database,a.pairs)
