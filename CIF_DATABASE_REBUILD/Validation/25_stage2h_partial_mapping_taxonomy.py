"""Stage 2H-A: evidence-labelled partial-mapping taxonomy (read-only)."""
from __future__ import annotations
import csv,sqlite3,statistics
from collections import Counter,defaultdict
from pathlib import Path
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'outputs'
def read(p):
 with open(p,newline='') as f:return list(csv.DictReader(f))
def write(name,rows):
 rows=list(rows);fields=list(rows[0]) if rows else []
 with open(OUT/name,'w',newline='') as f:w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
def main():
 pre=read(OUT/'full_mapping_preflight_results.csv');detail={x['ligand_instance_id']:x for x in read(OUT/'FULL_MAPPING_PREFLIGHT_ADAPTER_DETAILS.csv')};heavy={x['limitation_instance_id']:x for x in read(OUT/'STAGE2G_HEAVY_ATOM_PAIR_AUDIT.csv')};same={x['limitation_instance_id'] for x in read(OUT/'STAGE2F_PAIRED_OCCURRENCES.csv')};out=[]
 for r in pre:
  status=r['mapping_status'];iid=r['ligand_instance_id'];d=detail.get(iid,{});h=heavy.get(iid,{})
  if status=='mapping_timeout':klass,queue='P7_MCS_TIMEOUT','QUEUE_COMPLEX_GRAPH'
  elif status=='failed':klass,queue='P8_NO_MCS_MAPPING','QUEUE_FRAGMENT_REVIEW'
  elif status=='partial_valid_deposited_difference':klass,queue='P4_TRUE_DEPOSITED_REFERENCE_DIFFERENCE','QUEUE_NONE_VALID_PARTIAL'
  elif status=='partial_adapter_limitation' and h.get('heavy_atom_set_equal')=='True' and h.get('heavy_heavy_connectivity_equal')=='True':klass,queue='P1_EXACT_HEAVY_TOPOLOGY_MCS_PARTIAL','QUEUE_EXACT_GRAPH'
  elif status=='partial_adapter_limitation' and r['component_id']=='2NC':klass,queue='P3_SMALL_ATOMSET_DIFFERENCE','QUEUE_LOCAL_RECOVERY'
  elif status=='partial_adapter_limitation' and r['component_id'] in {'PI6','YDP'}:klass,queue='P5_LARGE_GRAPH_DIFFERENCE','QUEUE_CURATION_REVIEW'
  elif status=='partial_adapter_limitation':klass,queue='P11_OTHER','QUEUE_BOND_TEMPLATE'
  else:continue
  hs=d.get('selected_heavy_atom_count','');hr=d.get('smiles_heavy_atom_count','');hm=d.get('mapped_heavy_atom_count','');frac=d.get('heavy_atom_mapping_fraction','')
  outcome='timeout' if status=='mapping_timeout' else ('failed' if status=='failed' else 'partial')
  out.append({'ligand_instance_id':iid,'pdb_id':r['pdb_id'],'component_id':r['component_id'],'mapping_outcome':outcome,'mapping_reason_class':klass,'mapping_method':'legacy_mcs_etkdg_uff','heavy_atoms_structural':hs,'heavy_atoms_reference':hr,'heavy_atoms_mapped':hm,'heavy_atom_mapping_fraction':frac,'mapping_complete':0,'algorithm_queue':queue,'downstream_mapping_eligibility':'partial_only' if queue=='QUEUE_NONE_VALID_PARTIAL' else 'not_complete','preflight_status':status,'reason_code':r['reason_code'],'secondary_tags':';'.join(x for x in [('same_component_success' if iid in same else ''),('identical_heavy_topology' if h.get('heavy_atom_set_equal')=='True' else '')] if x)})
 write('PARTIAL_MAPPING_TAXONOMY.csv',out);write('MAPPING_ALGORITHM_QUEUES.csv',out)
 by=defaultdict(list)
 for x in out:by[x['mapping_reason_class']].append(x)
 lines=['# Partial Mapping Taxonomy','',f'Total partial/timeout/failed records classified: {len(out)}.','', '## Classes']
 for k,v in sorted(by.items()):
  fr=[float(x['heavy_atom_mapping_fraction']) for x in v if x['heavy_atom_mapping_fraction']!=''];comps=Counter(x['component_id'] for x in v);lines.append(f'- {k}: {len(v)} instances; {len(comps)} components; heavy fraction median={statistics.median(fr) if fr else "not recomputed"}; top components={comps.most_common(25)}')
 queues=Counter(x['algorithm_queue'] for x in out);lines+=['','## Queue priority',f'- QUEUE_EXACT_GRAPH: {queues["QUEUE_EXACT_GRAPH"]} evidence-backed P1 candidates; first priority.',f'- QUEUE_COMPLEX_GRAPH: {queues["QUEUE_COMPLEX_GRAPH"]} FindMCS timeouts; second priority.',f'- QUEUE_NONE_VALID_PARTIAL: {queues["QUEUE_NONE_VALID_PARTIAL"]} valid partials; no algorithm recovery required.']
 (OUT/'PARTIAL_MAPPING_TAXONOMY.md').write_text('\n'.join(lines)+'\n');(OUT/'STAGE2H_PARTIAL_MAPPING_CLASSIFICATION_DECISION.md').write_text('\n'.join(lines)+f'\n\nProduction generation may retain these explicit fields, but bulk mapping is not approved until a subsequent full preflight verifies remaining adapter queues. Downstream scripts must use `mapping_outcome`, `mapping_complete`, `mapping_reason_class`, and `downstream_mapping_eligibility`, never infer completeness from a partial status string.\n')
 print('taxonomy classified='+str(len(out))+' classes='+str(dict(Counter(x['mapping_reason_class'] for x in out))),flush=True)
if __name__=='__main__':main()
