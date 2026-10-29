"""Compute five finite flows, relative horizon entropies and the HM saddle."""
import argparse,hashlib,json
from pathlib import Path
from recursive_horizons.nsc_atlas_volume_routes import calculation

ROOT=Path(__file__).resolve().parents[1]

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    out=(ROOT/p.parse_args().output).resolve()
    if not out.is_relative_to(ROOT) or out.exists():raise ValueError('new repository output required')
    r=calculation()
    r['source_sha256']={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in (
        'src/recursive_horizons/nsc_atlas_volume_routes.py','scripts/check_nsc_atlas_volume_routes.py',
        'src/recursive_horizons/nsc_atlas_parent_scale_selector.py')}
    out.write_text(json.dumps(r,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
    for row in r['frontier']:print(row['d'],row['j'],row['Lambda_Lref0_squared'],row['endpoint']['s'],row['entropy_over_reference'])
    print(r['same_action_Hawking_Moss'])

if __name__=='__main__':main()
