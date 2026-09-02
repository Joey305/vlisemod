"""Lightweight, source-path-derived virus/target context audit."""
from __future__ import annotations
import csv,re,sqlite3
from collections import Counter,defaultdict
from pathlib import Path
from importlib import import_module
ROOT=Path(__file__).resolve().parent;SOURCE=ROOT.parent/'PDB_FILES';OUT=ROOT/'outputs'
c=import_module('00_common')
def virus(raw):
 key=raw.lower().replace('-','_').replace(' ','_')
 fixed={'hiv_1':'HIV-1','sars_cov_2':'SARS-CoV-2','hpv_16':'Human papillomavirus type 16','human_papillomavirus_type_16':'Human papillomavirus type 16','human_papillomavirus_16':'Human papillomavirus type 16'}
 if key in fixed:return fixed[key]
 match=re.fullmatch(r'(?:hpv_|human_papillomavirus_(?:type_)?)?(\d+[a-z]?)',key)
 return f'Human papillomavirus type {match.group(1)}' if match else raw.replace('_',' ')
def target(raw):
 fixed={'protease':'Protease','capsid_protein':'Capsid protein','reverse_transcriptase':'Reverse transcriptase','integrase':'Integrase','polymerase':'Polymerase','spike_protein':'Spike protein'}
 return fixed.get(raw,raw.replace('_',' ').capitalize())
def write(name,rows):
 rows=list(rows);fields=list(rows[0]) if rows else []
 with open(OUT/name,'w',newline='') as f:w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
def main(database):
 c.create_schema(database)
 paths=sorted(SOURCE.rglob('*.cif'));records=[]
 for p in paths:
  rel=p.relative_to(SOURCE);parts=rel.parts
  if len(parts)<3:continue
  raw,ctx=parts[0],parts[1];records.append({'pdb_id':p.stem.upper(),'virus_raw':raw,'virus_normalized':virus(raw),'target_context_raw':ctx,'target_context_normalized':target(ctx),'source_cif_path':str(p),'classification_source':'PDB_FILES_directory_hierarchy','classification_qc':''})
 by=defaultdict(list)
 for r in records:by[r['pdb_id']].append(r)
 coverage=[];conflicts=[]
 for pid,xs in by.items():
  viruses={x['virus_normalized'] for x in xs};contexts={x['target_context_normalized'] for x in xs};qc='conflicting_virus' if len(viruses)>1 else ('multiple_target_contexts' if len(contexts)>1 else 'consistent')
  for x in xs:x['classification_qc']=qc
  coverage.append({'pdb_id':pid,'placement_count':len(xs),'virus_assignment_count':len(viruses),'target_context_assignment_count':len(contexts),'virus_normalized':';'.join(sorted(viruses)),'target_context_normalized':';'.join(sorted(contexts)),'classification_qc':qc})
  if qc!='consistent':conflicts.append(coverage[-1])
 with sqlite3.connect(database) as db:
  db.row_factory=sqlite3.Row
  db.execute('PRAGMA foreign_keys=ON');db.execute('DELETE FROM structure_context')
  db.executemany('INSERT INTO structure_context(pdb_id,virus_raw,virus_normalized,target_context_raw,target_context_normalized,source_cif_path,classification_source,classification_qc) VALUES(:pdb_id,:virus_raw,:virus_normalized,:target_context_raw,:target_context_normalized,:source_cif_path,:classification_source,:classification_qc)',records)
  retained={r[0] for r in db.execute('SELECT DISTINCT entry_id FROM structures')};lig_total=db.execute('SELECT count(*) FROM ligand_instances').fetchone()[0];lig_covered=db.execute('''SELECT count(*) FROM ligand_instances i JOIN structures s ON s.structure_id=i.structure_id WHERE EXISTS(SELECT 1 FROM structure_context c WHERE c.pdb_id=s.entry_id)''').fetchone()[0]
  taxonomy=OUT/'PARTIAL_MAPPING_TAXONOMY.csv';queues=[]
  if taxonomy.exists():
   tax=list(csv.DictReader(open(taxonomy)));db.execute('CREATE TEMP TABLE tax(iid TEXT,queue TEXT)');db.executemany('INSERT INTO tax VALUES(?,?)',[(x['ligand_instance_id'],x['algorithm_queue']) for x in tax]);queues=[dict(r) for r in db.execute('''SELECT c.virus_normalized,c.target_context_normalized,t.queue AS algorithm_queue,count(*) AS instance_count FROM tax t JOIN ligand_instances i ON i.ligand_instance_id=t.iid JOIN structures s ON s.structure_id=i.structure_id JOIN structure_context c ON c.pdb_id=s.entry_id GROUP BY 1,2,3 ORDER BY 1,2,3''').fetchall()]
 write('STRUCTURE_CONTEXT_COVERAGE.csv',coverage);write('STRUCTURE_CONTEXT_CONFLICTS.csv',conflicts);write('MAPPING_QC_BY_VIRUS_CONTEXT.csv',queues)
 multiple_contexts=sum(x['classification_qc']=='multiple_target_contexts' for x in coverage);virus_conflicts=sum(x['classification_qc']=='conflicting_virus' for x in coverage);three_paths='; '.join(x['source_cif_path'] for x in by['3EKY']);three_contexts=[(x['virus_normalized'],x['target_context_normalized']) for x in by['3EKY']]
 md=['# Structure Context QC','',f'- Unique retained PDB IDs: {len(retained)}',f'- PDBs with virus assignment: {len(retained & set(by))}',f'- PDBs with target-context assignment: {len(retained & set(by))}',f'- PDBs with multiple target contexts: {multiple_contexts}',f'- PDBs with conflicting virus assignments: {virus_conflicts}',f'- Ligand instances with context coverage: {lig_covered}/{lig_total}',f'- 3EKY source paths: {three_paths}',f'- 3EKY normalized virus/context: {three_contexts}']
 (OUT/'STRUCTURE_CONTEXT_QC.md').write_text('\n'.join(md)+'\n');(OUT/'STRUCTURE_CONTEXT_RELEASE_DECISION.md').write_text('\n'.join(md)+f'\n\nEvery retained PDB with a source placement can inherit virus/context through `structure_context`; multiple contexts are retained. Mapping queues are groupable by virus/context. No ligand, SASA, or mapping facts were changed.\n')
 print(f'context placements={len(records)} retained_pdbs={len(retained)} conflicts={sum(x["classification_qc"]=="conflicting_virus" for x in coverage)}',flush=True)
if __name__=='__main__':main(str(ROOT/'viral_data_cif_v2.db'))
