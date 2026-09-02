"""Compatibility entrypoint: occurrence identification is atom-site ingestion in 03."""
from importlib import import_module
from pathlib import Path
import argparse
i=import_module('03_ingest_structures'); c=import_module('00_common')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--source',default=str(c.DEFAULT_SOURCE));p.add_argument('--database',default=str(c.ROOT/'viral_data_cif_v2.db'));p.add_argument('--limit',type=int);p.add_argument('--pdb-id');a=p.parse_args();print(i.ingest(Path(a.source),a.database,a.limit,a.pdb_id,True))
