"""Resolve included component chemistry from retained local sources only."""
from __future__ import annotations
import argparse,csv,json,sqlite3
from collections import defaultdict
from pathlib import Path
from importlib import import_module
from rdkit import Chem
c=import_module('00_common')
PROJECT=c.ROOT.parent
VERSION='local-chemistry-v2.1'
def component_file():
 out=defaultdict(list); p=PROJECT/'Components-smiles-stereo-oe.smi'
 if p.exists():
  with p.open(encoding='utf8',errors='replace') as f:
   for n,line in enumerate(f,1):
    bits=line.rstrip('\n').split('\t')
    if len(bits)>=2 and bits[1].strip(): out[bits[1].strip().upper()].append((bits[0].strip(),'Components-smiles-stereo-oe.smi',str(p),1))
 return out
def legacy_sources():
 out=defaultdict(list); path=PROJECT/'viral_data.db'
 if not path.exists(): return out
 db=sqlite3.connect(f'file:{path}?mode=ro',uri=True)
 for table in ('Functional_GROUPED','Ligand_Atoms_Smiles'):
  try:
   for comp,smiles in db.execute(f'SELECT ligand,smiles FROM {table} WHERE smiles IS NOT NULL'):
    if comp and smiles: out[str(comp).upper()].append((str(smiles).strip(),f'viral_data.db:{table}',str(path),10))
  except sqlite3.Error: pass
 db.close();return out
def canonical(smiles):
 mol=Chem.MolFromSmiles(smiles)
 return (Chem.MolToSmiles(mol,canonical=True),Chem.MolToInchiKey(mol)) if mol else (None,None)
def resolve(database,limit=None,pdb_id=None):
 c.create_schema(database); c.dirs(); sources=component_file(); legacy=legacy_sources()
 for k,v in legacy.items():sources[k].extend(v)
 with c.dbconn(database) as db:
  where="WHERE EXISTS (SELECT 1 FROM ligand_instances i JOIN structures s ON s.structure_id=i.structure_id WHERE i.ligand_id=l.ligand_id AND i.curation_status='included')"; args=[]
  if pdb_id: where+=" AND EXISTS (SELECT 1 FROM ligand_instances i JOIN structures s ON s.structure_id=i.structure_id WHERE i.ligand_id=l.ligand_id AND s.entry_id=?)";args=[pdb_id]
  ligands=db.execute(f'SELECT ligand_id,component_id FROM ligands l {where} ORDER BY component_id',args).fetchall(); ligands=ligands[:limit] if limit else ligands
  rid=c.run_start(db,'chemistry',{'limit':limit,'pdb_id':pdb_id,'network':False})
  report=[]; counts=defaultdict(int)
  for row in ligands:
   lid,comp=row['ligand_id'],row['component_id']; candidates=sources.get(comp.upper(),[])
   parsed=[]
   for smi,name,loc,prio in candidates:
    can,key=canonical(smi); status='valid' if can else 'invalid'
    db.execute('INSERT OR IGNORE INTO ligand_chemistry_sources(ligand_id,source_name,source_smiles,source_locator,source_priority,parse_status) VALUES(?,?,?,?,?,?)',(lid,name,smi,loc,prio,status))
    if can: parsed.append((can,key,smi,name,loc,prio))
   distinct={x[0] for x in parsed}
   if not candidates: status='missing_smiles';chosen=(None,None,None,None,None,None)
   elif not parsed: status='invalid_smiles';chosen=(None,None,None,None,None,None)
   elif len(distinct)>1: status='conflicting_sources';chosen=sorted(parsed,key=lambda x:(x[5],x[3],x[2]))[0]
   else: status='resolved';chosen=sorted(parsed,key=lambda x:(x[5],x[3],x[2]))[0]
   can,key,smi,name,loc,prio=chosen
   db.execute('UPDATE ligands SET smiles=?,smiles_source=?,source_version=?,canonical_smiles=?,inchikey=?,chemical_definition_version=?,chemical_status=? WHERE ligand_id=?',(smi,name,VERSION,can,key,VERSION,status,lid))
   counts[status]+=1;report.append({'ligand_id':lid,'component_id':comp,'status':status,'source_smiles':smi or '','canonical_smiles':can or '','smiles_source':name or '','source_locator':loc or '','source_count':len(candidates)})
  c.run_end(db,rid,'completed',len(ligands),len(ligands)-counts['missing_smiles']-counts['invalid_smiles'],counts['conflicting_sources'],0)
 out=c.ROOT/'outputs';
 with (out/'CHEMISTRY_RESOLUTION.csv').open('w',newline='',encoding='utf8') as f:w=csv.DictWriter(f,fieldnames=report[0].keys() if report else ['ligand_id']);w.writeheader();w.writerows(report)
 five=[r for r in report if len(r['component_id'])==5]
 with (out/'FIVE_CHARACTER_COMPONENT_CHEMISTRY.csv').open('w',newline='',encoding='utf8') as f:w=csv.DictWriter(f,fieldnames=report[0].keys() if report else ['ligand_id']);w.writeheader();w.writerows(five)
 audit=['# Chemistry source audit','',f'* Unique included chemical identities processed: {len(report)}']+[f'* {k}: {v}' for k,v in sorted(counts.items())]+[f'* Five-character identities: {len(five)}',f"* Five-character identities resolved: {sum(x['status']=='resolved' for x in five)}",'', 'Sources: retained Components-smiles-stereo-oe.smi (precedence 1), then read-only legacy Functional_GROUPED and Ligand_Atoms_Smiles (precedence 10). Distinct canonical SMILES are recorded as conflicts, never silently erased.']
 (out/'CHEMISTRY_SOURCE_AUDIT.md').write_text('\n'.join(audit)+'\n')
 return report
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--database',default=str(c.ROOT/'viral_data_cif_v2.db'));p.add_argument('--limit',type=int);p.add_argument('--pdb-id');a=p.parse_args();print(f'resolved {len(resolve(a.database,a.limit,a.pdb_id))} identities')
