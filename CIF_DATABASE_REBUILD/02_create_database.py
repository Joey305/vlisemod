from pathlib import Path
import argparse
from importlib import import_module
c=import_module('00_common')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--database',default=str(c.ROOT/'viral_data_cif_v2.db'));a=p.parse_args();c.create_schema(a.database);print(Path(a.database).resolve())
