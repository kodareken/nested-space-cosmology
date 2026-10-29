import argparse,json,hashlib
from pathlib import Path
from recursive_horizons.nsc_atlas_wall_exit import calculation
ROOT=Path(__file__).resolve().parents[1]
def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True)
    out=(ROOT/p.parse_args().output).resolve()
    if not out.is_relative_to(ROOT) or out.exists():raise ValueError('new repository output required')
    r=calculation();r['source_sha256']={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in (
        'src/recursive_horizons/nsc_atlas_wall_exit.py','scripts/check_nsc_atlas_wall_exit.py')}
    out.write_text(json.dumps(r,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
    print(r['reference_projection']);print(r['HM'])
if __name__=='__main__':main()
