import argparse,hashlib,json
from pathlib import Path
from recursive_horizons.nsc_atlas_curvature_mechanisms import curvature_audit
ROOT=Path(__file__).resolve().parents[1]
def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True)
    out=(ROOT/p.parse_args().output).resolve()
    if not out.is_relative_to(ROOT) or out.exists():raise ValueError('new repository output required')
    r=curvature_audit();r['source_sha256']={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in (
        'src/recursive_horizons/nsc_atlas_curvature_mechanisms.py','scripts/check_nsc_atlas_curvature_mechanisms.py')}
    out.write_text(json.dumps(r,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n');print(json.dumps(r))
if __name__=='__main__':main()
