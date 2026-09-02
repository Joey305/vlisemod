"""Freeze and profile the retained CIF corpus without changing it."""
from __future__ import annotations
import argparse,csv,json,os
from collections import Counter,defaultdict
from pathlib import Path
from importlib import import_module
c=import_module('00_common')

FIELDS=['entry_id','filename','absolute_path','relative_path','virus_label','protein_label','sha256','file_size','modification_time_utc','parse_status','parse_error','deposited_model_count','atom_site_count','het_occurrence_count']
def inventory(source,limit=None,pdb_id=None):
 c.dirs(); rows=[]; lengths=defaultdict(Counter); examples=defaultdict(list)
 for path in c.source_paths(source,limit,pdb_id):
  virus,protein,rel=c.classification(path,source); base={'filename':path.name,'absolute_path':str(path.resolve()),'relative_path':rel,'virus_label':virus,'protein_label':protein,'sha256':c.sha256(path),'file_size':path.stat().st_size,'modification_time_utc':__import__('datetime').datetime.fromtimestamp(path.stat().st_mtime,__import__('datetime').timezone.utc).isoformat(),'parse_status':'success','parse_error':'','entry_id':'','deposited_model_count':0,'atom_site_count':0,'het_occurrence_count':0}
  try:
   # Inventory is deliberately columnar: do not materialize 19-key dictionaries
   # for every atom in the 11k-file corpus just to count/profiling fields.
   block=c.cif_doc(path).sole_block();entry=c.getv(block,'_entry.id',path.stem);base['entry_id']=entry
   col=lambda tag:list(block.find_loop(tag))
   group,model,lasym,lcomp,aasym,acomp,aseq,ins=(col(t) for t in ['_atom_site.group_PDB','_atom_site.pdbx_PDB_model_num','_atom_site.label_asym_id','_atom_site.label_comp_id','_atom_site.auth_asym_id','_atom_site.auth_comp_id','_atom_site.auth_seq_id','_atom_site.pdbx_PDB_ins_code'])
   base['atom_site_count']=len(group)
   base['deposited_model_count']=len({c.norm(x) for x in model})
   # Keep group/identity columns aligned without reliance on CIF loop order.
   het={tuple(c.norm(x) for x in values[1:]) for values in zip(group,model,lasym,lcomp,aasym,aseq,ins) if c.norm(values[0]).upper()=='HETATM'}
   base['het_occurrence_count']=len(het)
   for field,vals in [('_entry.id',[entry]),('_atom_site.label_comp_id',lcomp),('_atom_site.auth_comp_id',acomp),('_atom_site.auth_asym_id',aasym),('_atom_site.label_asym_id',lasym)]:
    vals=[c.norm(x) for x in vals]
    for v in vals:
     if v: lengths[field][len(v)]+=1; examples[field].append(v)
  except Exception as e: base.update(parse_status='failed',parse_error=f'{type(e).__name__}: {e}')
  rows.append(base)
 mdir=c.ROOT/'manifests'; out=mdir/'FROZEN_CIF_CORPUS_MANIFEST.csv'
 with out.open('w',newline='',encoding='utf8') as f:w=csv.DictWriter(f,fieldnames=FIELDS);w.writeheader();w.writerows(rows)
 with (mdir/'IDENTIFIER_LENGTH_REPORT.csv').open('w',newline='',encoding='utf8') as f:
  w=csv.DictWriter(f,fieldnames=['field','length','count','examples']);w.writeheader()
  for field,counts in lengths.items():
   for length,count in sorted(counts.items()):w.writerow({'field':field,'length':length,'count':count,'examples':' | '.join(sorted(set(x for x in examples[field] if len(x)==length))[:10])})
 byid=defaultdict(list);bysha=defaultdict(list)
 for r in rows:
  if r['parse_status']=='success':byid[r['entry_id']].append(r);bysha[r['sha256']].append(r)
 duplicate_ids={k:v for k,v in byid.items() if len(v)>1}; same_id_diff={k:v for k,v in duplicate_ids.items() if len({x['sha256'] for x in v})>1}
 longcomp=sorted({x for x in examples['_atom_site.label_comp_id'] if len(x)>3})
 report=['# Frozen CIF corpus report','',f'* Files: {len(rows)}',f'* Parse successes: {sum(r["parse_status"]=="success" for r in rows)}',f'* Parse failures: {sum(r["parse_status"]=="failed" for r in rows)}',f'* Unique entry IDs: {len(byid)}',f'* Duplicate entry IDs: {len(duplicate_ids)}',f'* Identical-checksum duplicate groups: {sum(len(v)>1 for v in bysha.values())}',f'* Same entry ID, different checksums: {len(same_id_diff)}',f'* Component IDs >3 characters: {len(longcomp)}',f'* Examples: {", ".join(longcomp[:20]) or "none"}','', 'Conflicting same-ID/different-checksum files are retained as distinct revisions; nothing was silently selected.']
 (mdir/'FROZEN_CIF_CORPUS_REPORT.md').write_text('\n'.join(report)+'\n')
 idreport=['# Identifier length report','', 'All values are textual mmCIF tokens. Schema columns are SQLite TEXT and no length restriction is applied.','',f'* Component IDs longer than 3: {", ".join(longcomp[:100]) or "none"}']
 (mdir/'IDENTIFIER_LENGTH_REPORT.md').write_text('\n'.join(idreport)+'\n')
 return rows
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--source',default=str(c.DEFAULT_SOURCE));p.add_argument('--limit',type=int);p.add_argument('--pdb-id');a=p.parse_args();r=inventory(Path(a.source),a.limit,a.pdb_id);print(f'inventoried {len(r)} CIF files')
