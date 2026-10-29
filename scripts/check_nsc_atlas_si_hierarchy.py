"""Retain the SI hierarchy audit without selecting g from an observation."""
import argparse,hashlib,json,math
from pathlib import Path
from recursive_horizons.nsc_atlas_si_hierarchy import calculation,A0,VSTAR,B_HM
ROOT=Path(__file__).resolve().parents[1]

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    out=(ROOT/p.parse_args().output).resolve()
    if not out.is_relative_to(ROOT) or out.exists():raise ValueError('new repository output required')
    pins=('results/summaries/nsc-atlas-volume-lifetime-v1.json',
          'results/summaries/nsc-atlas-parent-volume-selector-v1.json',
          'results/summaries/nsc-atlas-SI-scale-ladder-v1.json')
    lifetime,volume,scale=[json.loads((ROOT/n).read_text(encoding='utf-8')) for n in pins]
    minimum=next(v for v in volume['vacuum_extrema'] if v['kind']=='minimum')
    assert math.isclose(B_HM,lifetime['rows'][0]['HM_exponent_times_g'],rel_tol=1e-14)
    assert math.isclose(VSTAR,minimum['vacuum'],rel_tol=1e-14)
    assert math.isclose(A0,scale['fixed_parent']['Lref_over_model_unit'],rel_tol=1e-14)
    r=calculation();r['input_sha256']={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in pins}
    r['source_sha256']={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in (
        'src/recursive_horizons/nsc_atlas_si_hierarchy.py','scripts/check_nsc_atlas_si_hierarchy.py')}
    out.write_text(json.dumps(r,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
    print(r['targets']);print(r['forward_identities'])
    for row in r['exponential_inverse_diagnostics']:
        print(row['p'],row['required_g_using_displayed_log'],row['same_geometry']['log_hierarchy'],row['same_geometry_log_miss_percent'])
    print(r['target_inverse_consistency']);print(r['stopping_outcome'])

if __name__=='__main__':main()
