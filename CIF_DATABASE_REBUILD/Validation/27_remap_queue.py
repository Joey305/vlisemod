"""Create/query the Stage 2I remediation registry; never rebuilds foundation data."""
from __future__ import annotations
import argparse,csv,sqlite3
from pathlib import Path
from importlib import import_module
c=import_module('00_common');ROOT=Path(__file__).resolve().parent;OUT=ROOT/'outputs';SKIP={'QUEUE_EXACT_GRAPH','QUEUE_LOCAL_RECOVERY','QUEUE_BOND_TEMPLATE','QUEUE_COMPLEX_GRAPH','QUEUE_FRAGMENT_REVIEW','QUEUE_CURATION_REVIEW'}
def load_taxonomy():
 with open(OUT/'PARTIAL_MAPPING_TAXONOMY.csv',newline='') as f:return [r for r in csv.DictReader(f) if r['algorithm_queue'] in SKIP]
def provision(database):
 c.create_schema(database);tax=load_taxonomy()
 with c.dbconn(database) as db:
  run_id=c.run_start(db,'mapping_remediation_registry',{'source':'PARTIAL_MAPPING_TAXONOMY.csv','skip_queues':sorted(SKIP)})
  rows=[]
  for x in tax:
   context=db.execute('''SELECT group_concat(DISTINCT virus_normalized),group_concat(DISTINCT target_context_normalized) FROM ligand_instances i JOIN structures s ON s.structure_id=i.structure_id LEFT JOIN structure_context c ON c.pdb_id=s.entry_id WHERE i.ligand_instance_id=?''',(x['ligand_instance_id'],)).fetchone();rows.append((x['ligand_instance_id'],x['pdb_id'],x['component_id'],context[0],context[1],x['mapping_reason_class'],x['algorithm_queue'],x['preflight_status'],x['heavy_atoms_structural'] or None,x['heavy_atoms_reference'] or None,x['heavy_atoms_mapped'] or None,x['heavy_atom_mapping_fraction'] or None,'pending',run_id,'Stage 2I QC-controlled mapping deferral'))
  db.executemany('''INSERT INTO mapping_remediation_queue(ligand_instance_id,pdb_id,component_id,virus_normalized,target_context_normalized,mapping_reason_class,algorithm_queue,preflight_mapping_status,heavy_atoms_structural,heavy_atoms_reference,heavy_atoms_mapped,heavy_atom_mapping_fraction,remediation_status,first_identified_run_id,notes) VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?) ON CONFLICT(ligand_instance_id) DO UPDATE SET mapping_reason_class=excluded.mapping_reason_class,algorithm_queue=excluded.algorithm_queue,preflight_mapping_status=excluded.preflight_mapping_status,notes=excluded.notes''',rows)
  c.run_end(db,run_id,'completed',len(rows),len(rows));manifest=[dict(r) for r in db.execute('SELECT ligand_instance_id,pdb_id,component_id,virus_normalized,target_context_normalized,mapping_reason_class,algorithm_queue,heavy_atom_mapping_fraction,remediation_status,notes FROM mapping_remediation_queue WHERE remediation_status="pending" ORDER BY ligand_instance_id')]
 fields=list(manifest[0])
 with open(OUT/'PRODUCTION_MAPPING_SKIPPED_INSTANCES.csv','w',newline='') as f:w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(manifest)
 (OUT/'PRODUCTION_MAPPING_SKIPPED_INSTANCES.md').write_text(f'# Provisional production mapping skip manifest\n\nPending exact-instance deferrals: {len(manifest)}. These rows remain in ligand instances, SASA, context, and all mapping-independent layers.\n')
 print(f'remediation registry pending={len(manifest)}',flush=True)
def main():
 p=argparse.ArgumentParser();p.add_argument('--database',default=str(ROOT/'viral_data_cif_v2.db'));p.add_argument('--queue',choices=sorted(SKIP));p.add_argument('--provision-only',action='store_true');a=p.parse_args();provision(a.database)
if __name__=='__main__':main()
