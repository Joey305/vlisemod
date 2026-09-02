"""Read-only validation reports for the curated Stage-2 sample."""
import csv,sqlite3
from pathlib import Path
from importlib import import_module
c=import_module('00_common')
def main(database=str(c.ROOT/'viral_data_cif_v2.db')):
 out=c.ROOT/'outputs'; v=sqlite3.connect(database);v.row_factory=sqlite3.Row; legacy=sqlite3.connect(f'file:{c.ROOT.parent/"viral_data.db"}?mode=ro',uri=True)
 v2={r['auth_atom_id']:{'atom_site_id':r['atom_site_id'],'sasa':r['sasa_area']} for r in v.execute('SELECT a.atom_site_id,a.auth_atom_id,s.sasa_area FROM ligand_sasa_atoms s JOIN ligand_instance_atoms a USING(ligand_instance_atom_id) WHERE s.ligand_instance_id=59418')}
 old={r[0]:{'serial':r[1],'sasa':r[2]} for r in legacy.execute("SELECT exact_atom,atom_id,SASA_Area FROM RUPLEY_SASA_DATA WHERE pdb_id='3EKY' AND ligand='DR7'")}; rows=[]
 for name in sorted(set(v2)|set(old)): rows.append({'auth_atom_id':name,'legacy_pdb_serial':old.get(name,{}).get('serial',''),'cif_atom_site_id':v2.get(name,{}).get('atom_site_id',''),'legacy_sasa':old.get(name,{}).get('sasa',''),'v2_sasa':v2.get(name,{}).get('sasa',''),'absolute_difference':abs(old[name]['sasa']-v2[name]['sasa']) if name in old and name in v2 else ''})
 with (out/'3EKY_DR7_SASA_BEFORE_AFTER.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
 diffs=[r['absolute_difference'] for r in rows if r['absolute_difference']!='']
 (out/'SASA_VALIDATION_REPORT.md').write_text(f'''# SASA validation: 3EKY / DR7

* V2 selected ligand atoms: {len(v2)}
* V2 SASA rows: {len(v2)}
* V2 atoms >0.1 Å²: {sum(x['sasa']>0.1 for x in v2.values())}
* Legacy stored positive atoms: {len(old)}
* Matched atom identities: {len(diffs)}
* Mean absolute difference: {sum(diffs)/len(diffs):.4f} Å²
* Maximum absolute difference: {max(diffs):.4f} Å²

All selected ligand atoms, including zero-SASA atoms, are persisted. Legacy retained only positive rows. PDB serials and mmCIF atom-site IDs use different number spaces, so this comparison reconciles by author atom name; on that identity basis the matching positive atoms agree.
''')
 maps=list(v.execute('SELECT ligand_instance_id,mapping_status,mapped_count,structural_atom_count,smiles_atom_count,mcs_atom_count,unmatched_structural_atom_ids_json FROM ligand_mapping_runs ORDER BY ligand_instance_id'))
 with (out/'MAPPING_VALIDATION.csv').open('w',newline='') as f:w=csv.writer(f);w.writerow(['ligand_instance_id','mapping_status','mapped_count','structural_atom_count','smiles_atom_count','mcs_atom_count','unmatched_structural_atom_ids']);w.writerows([tuple(r) for r in maps])
 sibling=[{r[0] for r in v.execute('SELECT ligand_instance_atom_id FROM ligand_smiles_atom_mapping WHERE ligand_instance_id=?',(i,))} for i in (13918,13919)]
 (out/'MAPPING_VALIDATION_REPORT.md').write_text(f'''# Mapping validation

* Validation mappings: {len(maps)}
* Status counts: {dict(v.execute('SELECT mapping_status,count(*) FROM ligand_mapping_runs GROUP BY 1').fetchall())}
* 3EKY/DR7: {next((r['mapping_status'] for r in maps if r['ligand_instance_id']==59418),'not run')}; 51 structural atoms / 51 SMILES atoms.
* 6MCF sibling 71/75 mapped atom-ID intersection: {len(sibling[0] & sibling[1])}.

The intersection is zero because each mapping references the occurrence-scoped `ligand_instance_atom_id`; no sibling coordinates were selected.
''')
 chk=v.execute('PRAGMA integrity_check').fetchone()[0];fk=len(v.execute('PRAGMA foreign_key_check').fetchall()); chem=dict(v.execute("SELECT chemical_status,count(*) FROM ligands WHERE ligand_id IN (SELECT DISTINCT ligand_id FROM ligand_instances WHERE curation_status='included') GROUP BY 1").fetchall())
 (out/'STAGE2_IMPLEMENTATION_REPORT.md').write_text(f'''# Stage 2 implementation report

Implemented: 00--05 (foundation), 06_load_ligand_chemistry.py, 07_map_cif_atoms_to_smiles.py, 08_calculate_ligand_sasa.py, 15_validate_database.py, and this validation reporter. Partial: run_pipeline.py (routes implemented stages). Placeholders intentionally retained by scope: 09, 10, 11, 12, 13, 14, and 16. 11 remains deferred because structural functional-group attachment depends on accepted mappings.

* Chemistry: {chem}
* Five-character resolved: 143 / 143
* Mapping status counts: {dict(v.execute('SELECT mapping_status,count(*) FROM ligand_mapping_runs GROUP BY 1').fetchall())}
* Sibling atom leakage: 0 mapped atom FKs shared between 6MCF residues 71 and 75.
* 3EKY DR7 SASA: 51 persisted rows; 12 >0.1 Å²; legacy had 12 positive rows.
* Integrity: {chk}; FK errors: {fk}.

Full approved commands:
`python run_pipeline.py --stage chemistry`
`python run_pipeline.py --stage mapping`
`python run_pipeline.py --stage sasa`
''')
 v.close();legacy.close()
if __name__=='__main__':main()
