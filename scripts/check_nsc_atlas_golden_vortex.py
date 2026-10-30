"""Retain the vortex, source-stability and black-hole/BKL comparisons."""
import argparse,hashlib,json
from pathlib import Path
from recursive_horizons.nsc_atlas_golden_vortex import calculation
ROOT=Path(__file__).resolve().parents[1]


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    out=(ROOT/p.parse_args().output).resolve()
    if not out.is_relative_to(ROOT) or out.exists():raise ValueError('new repository output required')
    names=('results/summaries/nsc-atlas-pen-geometry-v1.json','results/summaries/nsc-atlas-parke-bubble-v1.json')
    r=calculation(*(json.loads((ROOT/n).read_text(encoding='utf-8')) for n in names))
    r['input_sha256']={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names}
    sources=('src/recursive_horizons/nsc_atlas_golden_vortex.py','src/recursive_horizons/nsc_atlas_pen_geometry.py',
             'scripts/check_nsc_atlas_golden_vortex.py')
    r['source_sha256']={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in sources}
    out.write_text(json.dumps(r,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
    print('Golden ratios:',r['golden']['angle'],r['golden']['per_full_turn_ratio'])
    print('Golden area:',r['golden_area']['measure']['area_fraction'],r['golden_area']['measure']['area_enclosure'])
    print('Sources:',[(x['omega0'],x['limit']['area_fraction']) for x in r['source_stability']])
    print('BKL fixed-point derivative:',r['bkl']['era_derivative_at_fixed_point'])
    print('No physical transition ratio or observer event is selected.')


if __name__=='__main__':main()
