"""Record all five HM coefficients and the bounded CDL classification."""
import argparse,json,hashlib
from pathlib import Path
from recursive_horizons.nsc_atlas_volume_lifetime import calculation
ROOT=Path(__file__).resolve().parents[1]
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    out=(ROOT/p.parse_args().output).resolve()
    if not out.is_relative_to(ROOT) or out.exists():raise ValueError('new repository output required')
    r=calculation();r['source_sha256']={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in (
        'src/recursive_horizons/nsc_atlas_volume_lifetime.py','scripts/check_nsc_atlas_volume_lifetime.py')}
    out.write_text(json.dumps(r,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
    for row in r['rows']:print(row['d'],row['j'],row['HM_exponent_times_g'],row['CDL_status'])
    print(r['reference_search']['status'])
if __name__=='__main__':main()
