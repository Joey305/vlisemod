"""Load exact CIF revisions, classifications, HET occurrences, and raw atom rows."""
from __future__ import annotations
import argparse,csv
from pathlib import Path
from importlib import import_module
c=import_module('00_common')
def ingest(source,database,limit=None,pdb_id=None,resume=False):
 c.create_schema(database); processed=ok=failed=instances=0
 manifest=c.ROOT/'manifests/FROZEN_CIF_CORPUS_MANIFEST.csv'; msha=c.sha256(manifest) if manifest.exists() else ''
 with c.dbconn(database) as db:
  rid=c.run_start(db,'ingest',{'source':str(source),'limit':limit,'pdb_id':pdb_id,'resume':resume},msha)
  for path in c.source_paths(source,limit,pdb_id):
   processed+=1
   try:
    entry,ats=c.atoms(path); checksum=c.sha256(path);virus,protein,rel=c.classification(path,source)
    db.execute('INSERT INTO structures(entry_id,source_cif_path,source_cif_sha256,file_size,source_status,parser,parser_version,ingest_time) VALUES(?,?,?,?,?,?,?,?) ON CONFLICT(entry_id,source_cif_sha256) DO UPDATE SET source_cif_path=excluded.source_cif_path',(entry,str(path.resolve()),checksum,path.stat().st_size,'retained','gemmi',c.PARSER_VERSION,c.now()))
    sid=db.execute('SELECT structure_id FROM structures WHERE entry_id=? AND source_cif_sha256=?',(entry,checksum)).fetchone()[0]
    db.execute('INSERT OR IGNORE INTO structure_classifications(structure_id,virus_label,protein_label,source_relative_path) VALUES(?,?,?,?)',(sid,virus,protein,rel))
    groups={}
    for a in ats:
     if c.norm(a['_atom_site.group_PDB']).upper()!='HETATM':continue
     key=tuple(c.norm(a[k]) for k in ['_atom_site.pdbx_PDB_model_num','_atom_site.label_asym_id','_atom_site.label_comp_id','_atom_site.auth_asym_id','_atom_site.auth_seq_id','_atom_site.pdbx_PDB_ins_code'])
     groups.setdefault(key,[]).append(a)
    for (model,label_asym,label_comp,auth_asym,auth_seq,ins),rs in groups.items():
     status,reason=c.curation(label_comp);identity='complete' if all([model,label_asym,label_comp,auth_asym,auth_seq]) else 'incomplete'
     db.execute("INSERT INTO ligands(component_id,chemical_definition_version,chemical_status) VALUES(?,?,?) ON CONFLICT(component_id,chemical_definition_version) DO NOTHING",(label_comp,'local-unresolved-v1','missing'))
     lid=db.execute('SELECT ligand_id FROM ligands WHERE component_id=? AND chemical_definition_version=?',(label_comp,'local-unresolved-v1')).fetchone()[0]
     one=rs[0]; db.execute('''INSERT INTO ligand_instances(structure_id,ligand_id,deposited_model_num,label_asym_id,label_comp_id,auth_asym_id,auth_comp_id,label_seq_id,auth_seq_id,insertion_code_raw,insertion_code_normalized,identity_status,curation_status,curation_reason,curation_rule_version) VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?) ON CONFLICT(structure_id,deposited_model_num,label_asym_id,label_comp_id,auth_asym_id,auth_seq_id,insertion_code_normalized) DO UPDATE SET curation_status=excluded.curation_status''',(sid,lid,model,label_asym,label_comp,auth_asym,c.norm(one['_atom_site.auth_comp_id']),c.norm(one['_atom_site.label_seq_id']),auth_seq,c.raw(one['_atom_site.pdbx_PDB_ins_code']),ins,identity,status,reason,'legacy_2_PDBSORT_2_cif_v2'))
     iid=db.execute('SELECT ligand_instance_id FROM ligand_instances WHERE structure_id=? AND deposited_model_num=? AND label_asym_id=? AND label_comp_id=? AND auth_asym_id=? AND auth_seq_id=? AND insertion_code_normalized=?',(sid,model,label_asym,label_comp,auth_asym,auth_seq,ins)).fetchone()[0]
     for a in rs:
      db.execute('''INSERT INTO ligand_instance_atoms(ligand_instance_id,atom_site_id,label_atom_id,auth_atom_id,element,x,y,z,altloc,occupancy,b_factor,deposited_model_num) VALUES(?,?,?,?,?,?,?,?,?,?,?,?) ON CONFLICT(ligand_instance_id,atom_site_id) DO NOTHING''',(iid,c.norm(a['_atom_site.id']),c.norm(a['_atom_site.label_atom_id']),c.norm(a['_atom_site.auth_atom_id']),c.norm(a['_atom_site.type_symbol']),c.fnum(a['_atom_site.Cartn_x']),c.fnum(a['_atom_site.Cartn_y']),c.fnum(a['_atom_site.Cartn_z']),c.norm(a['_atom_site.label_alt_id']),c.fnum(a['_atom_site.occupancy']),c.fnum(a['_atom_site.B_iso_or_equiv']),model))
     c.set_conformer_flags(db,iid);instances+=1
    ok+=1
   except Exception as e: failed+=1;c.fail(db,rid,'ingest',f'{type(e).__name__}: {e}',code='parse_or_ingest_failure')
  c.run_end(db,rid,'completed' if not failed else 'partial',processed,ok,0,failed)
 return processed,ok,failed,instances
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--source',default=str(c.DEFAULT_SOURCE));p.add_argument('--database',default=str(c.ROOT/'viral_data_cif_v2.db'));p.add_argument('--limit',type=int);p.add_argument('--pdb-id');p.add_argument('--resume',action='store_true');a=p.parse_args();print(ingest(Path(a.source),a.database,a.limit,a.pdb_id,a.resume))
