"""Post-run QC summary for the Stage 2I provisional mapping build."""
from __future__ import annotations
import csv,sqlite3
from pathlib import Path
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'outputs';VERSION='legacy_mcs_etkdg_uff_cif_v2.4'
def main(database):
 db=sqlite3.connect(database);rows=db.execute('SELECT mapping_status,mapping_outcome,count(*) FROM ligand_mapping_runs WHERE method_version=? GROUP BY 1,2',(VERSION,)).fetchall();skip=db.execute('SELECT count(*) FROM mapping_remediation_queue WHERE remediation_status="pending"').fetchone()[0]
 with open(OUT/'PROVISIONAL_V2_MAPPING_QC.csv','w',newline='') as f:w=csv.writer(f);w.writerow(['mapping_status','mapping_outcome','count']);w.writerows(rows)
 (OUT/'PROVISIONAL_V2_MAPPING_QC.md').write_text(f'# V-LiSEMOD CIF V2 PROVISIONAL / QC-AWARE BUILD\n\nPending remediation skips: {skip}. Mapping coverage is not 100%; only analyses requiring trusted atom mapping exclude these explicit queue records.\n')
 print(f'provisional mapping QC pending_skips={skip}',flush=True)
if __name__=='__main__':main(str(ROOT/'viral_data_cif_v2.db'))
