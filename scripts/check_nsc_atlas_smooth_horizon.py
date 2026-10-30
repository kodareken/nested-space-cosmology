"""Save a compact smooth-horizon calculation; never overwrite historical output."""
import argparse
import hashlib
import json
from pathlib import Path
from recursive_horizons.nsc_atlas_smooth_horizon import calculation

ROOT=Path(__file__).resolve().parents[1]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    out=(ROOT/p.parse_args().output).resolve()
    if not out.is_relative_to(ROOT) or out.exists():
        raise ValueError('new repository output required')
    result=calculation()
    names=('src/recursive_horizons/nsc_atlas_smooth_horizon.py',
           'src/recursive_horizons/nsc_atlas_parke_nucleation.py',
           'src/recursive_horizons/nsc_atlas_parent_scale_selector.py',
           'src/recursive_horizons/nsc_atlas_pen_geometry.py',
           'scripts/check_nsc_atlas_smooth_horizon.py')
    result['source_sha256']={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names}
    out.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
    print('Parke limit:',result['shell_limit']['limit_radius'])
    print('Same-action vacuum:',result['admitted_vacuum'])
    print('Radius-inversion diagnostic:',result['radius_inversion_control'])
    print('Symmetric eta control: 0; physical baryon/photon eta: unselected')


if __name__=='__main__':
    main()
