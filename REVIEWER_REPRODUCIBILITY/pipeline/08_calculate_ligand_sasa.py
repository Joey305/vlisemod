"""Occurrence-specific atom-level Shrake--Rupley SASA using retained CIF coordinates."""
from __future__ import annotations
import argparse
from collections import defaultdict
from importlib import import_module
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
c=import_module('00_common'); VERSION='biopython-shrake_rupley-1.40-cif-v2.1'
def record_line(a):
 serial=int(a['_atom_site.id']) if str(a['_atom_site.id']).isdigit() and int(a['_atom_site.id'])<100000 else 0
 name=(c.norm(a['_atom_site.auth_atom_id']) or c.norm(a['_atom_site.label_atom_id']) or c.norm(a['_atom_site.type_symbol']) or 'X')[:4].rjust(4);res=(c.norm(a['_atom_site.auth_comp_id']) or c.norm(a['_atom_site.label_comp_id']) or 'LIG')[:3].rjust(3);chain=(c.norm(a['_atom_site.auth_asym_id']) or 'A')[:1];seq=c.norm(a['_atom_site.auth_seq_id']) or '1'; ins=c.norm(a['_atom_site.pdbx_PDB_ins_code'])[:1] or ' '
 try:seqi=int(seq)
 except:seqi=1
 return f'{c.norm(a["_atom_site.group_PDB"]).upper():<6}{serial:5d} {name} {res} {chain}{seqi:4d}{ins}   {float(a["_atom_site.Cartn_x"]):8.3f}{float(a["_atom_site.Cartn_y"]):8.3f}{float(a["_atom_site.Cartn_z"]):8.3f}{(c.fnum(a["_atom_site.occupancy"]) or 1):6.2f}{(c.fnum(a["_atom_site.B_iso_or_equiv"]) or 0):6.2f}          {(c.norm(a["_atom_site.type_symbol"]) or "C")[:2].rjust(2)}'
def context_pdb(db,iid):
 target=db.execute('''SELECT i.*,s.source_cif_path FROM ligand_instances i JOIN structures s ON s.structure_id=i.structure_id WHERE i.ligand_instance_id=?''',(iid,)).fetchone(); chosen={str(r[0]) for r in db.execute('SELECT atom_site_id FROM ligand_instance_atoms WHERE ligand_instance_id=? AND selected_conformer=1',(iid,))}
 entry,atoms=c.atoms(__import__('pathlib').Path(target['source_cif_path'])); model=target['deposited_model_num']; lines=['MODEL        1']; seen_alt=set()
 for a in atoms:
  if c.norm(a['_atom_site.pdbx_PDB_model_num'])!=model:continue
  if c.norm(a['_atom_site.label_comp_id']).upper()=='HOH':continue
  sid=c.norm(a['_atom_site.id']); is_target=sid in chosen
  alt=c.norm(a['_atom_site.label_alt_id'])
  # Target uses central coherent-conformer flags. Context has a deterministic
  # per-atom fallback only because it is not an analyzed ligand occurrence.
  if not is_target and alt:
   key=(c.norm(a['_atom_site.label_asym_id']),c.norm(a['_atom_site.auth_seq_id']),c.norm(a['_atom_site.label_atom_id']))
   if key in seen_alt:continue
   seen_alt.add(key)
  if is_target or not any(sid==x for x in chosen): lines.append(record_line(a))
 lines+=['ENDMDL','END'];return '\n'.join(lines)+'\n',chosen,target
def sasa_one(db,iid):
 text,chosen,target=context_pdb(db,iid); parser=PDBParser(QUIET=True); structure=parser.get_structure('cif_context',__import__('io').StringIO(text));sr=ShrakeRupley(probe_radius=1.40,n_points=100);sr.compute(structure,level='A')
 vals={str(a.get_serial_number()):float(getattr(a,'sasa',0.0)) for a in structure.get_atoms() if str(a.get_serial_number()) in chosen}
 rows=db.execute('SELECT ligand_instance_atom_id,atom_site_id FROM ligand_instance_atoms WHERE ligand_instance_id=? AND selected_conformer=1',(iid,)).fetchall()
 return [(r['ligand_instance_atom_id'],r['atom_site_id'],vals.get(str(r['atom_site_id']),0.0)) for r in rows],target
def run(database,limit=None,pdb_id=None,instance_id=None):
 c.create_schema(database)
 with c.dbconn(database) as db:
  q="SELECT i.ligand_instance_id FROM ligand_instances i JOIN structures s ON s.structure_id=i.structure_id WHERE i.curation_status='included'";args=[]
  if pdb_id:q+=' AND s.entry_id=?';args.append(pdb_id)
  if instance_id:q+=' AND i.ligand_instance_id=?';args.append(instance_id)
  ids=[x[0] for x in db.execute(q+' ORDER BY i.ligand_instance_id',args)];ids=ids[:limit] if limit else ids;result=[]
  for iid in ids:
   # Materialized V2 facts are idempotent at occurrence/method version.
   db.execute('DELETE FROM ligand_sasa_atoms WHERE ligand_instance_id=? AND method_version=?',(iid,VERSION))
   rid=c.run_start(db,'sasa',{'ligand_instance_id':iid,'probe_radius':1.40,'n_points':100,'water_treatment':'HOH_removed'}); rows,target=sasa_one(db,iid)
   for aid,sid,value in rows:db.execute('INSERT OR REPLACE INTO ligand_sasa_atoms(run_id,ligand_instance_id,ligand_instance_atom_id,sasa_area,legacy_exposed,probe_radius,point_density,water_treatment,conformer_policy,status,deposited_model_num,method_version) VALUES(?,?,?,?,?,?,?,?,?,?,?,?)',(rid,iid,aid,value,int(value>0.1),1.4,100,'HOH_removed','central_coherent_altloc_v1','complete',target['deposited_model_num'],VERSION))
   c.run_end(db,rid,'completed',1,1,0,0);result.append((iid,len(rows),sum(v>0.1 for _,_,v in rows)))
 return result
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--database',default=str(c.ROOT/'viral_data_cif_v2.db'));p.add_argument('--limit',type=int);p.add_argument('--pdb-id');p.add_argument('--ligand-instance-id',type=int);a=p.parse_args();print(run(a.database,a.limit,a.pdb_id,a.ligand_instance_id))
