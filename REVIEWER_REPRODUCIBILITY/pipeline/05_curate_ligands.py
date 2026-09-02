"""Curation is deterministic during ingestion; this entrypoint documents/reapplies it."""
from importlib import import_module
import argparse
c=import_module('00_common')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--database',default=str(c.ROOT/'viral_data_cif_v2.db'));a=p.parse_args()
 with c.dbconn(a.database) as db:
  for r in db.execute('SELECT ligand_instance_id,label_comp_id FROM ligand_instances'):
   s,reason=c.curation(r['label_comp_id']);db.execute('UPDATE ligand_instances SET curation_status=?,curation_reason=?,curation_rule_version=? WHERE ligand_instance_id=?',(s,reason,'legacy_2_PDBSORT_2_cif_v2',r['ligand_instance_id']))
 print('curation reapplied')
