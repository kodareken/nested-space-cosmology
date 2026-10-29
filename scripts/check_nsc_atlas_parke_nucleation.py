"""Retain the radius/mass result or the source/reheating admission check."""
import argparse,hashlib,json
from pathlib import Path
from recursive_horizons.nsc_atlas_parke_nucleation import geometry_record,source_record
ROOT=Path(__file__).resolve().parents[1]

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--stage',choices=('geometry','source'),required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args();out=(ROOT/a.output).resolve()
    if not out.is_relative_to(ROOT) or out.exists():raise ValueError('new repository output required')
    r=geometry_record() if a.stage=='geometry' else source_record()
    r['source_sha256']={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in (
        'src/recursive_horizons/nsc_atlas_parke_nucleation.py','scripts/check_nsc_atlas_parke_nucleation.py')}
    out.write_text(json.dumps(r,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(r))
if __name__=='__main__':main()
