"""Save the bounded charge-map and Parke-ruler calculation."""
import argparse
import hashlib
import json
from pathlib import Path
from recursive_horizons.nsc_atlas_charge_inheritance import calculation

ROOT=Path(__file__).resolve().parents[1]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    out=(ROOT/parser.parse_args().output).resolve()
    if not out.is_relative_to(ROOT) or out.exists():
        raise ValueError('new repository output required')
    result=calculation()
    names=('src/recursive_horizons/nsc_atlas_charge_inheritance.py',
           'src/recursive_horizons/nsc_atlas_flux_transition.py',
           'src/recursive_horizons/nsc_atlas_parent_scale_selector.py',
           'src/recursive_horizons/nsc_atlas_parke_nucleation.py',
           'src/recursive_horizons/nsc_atlas_si_hierarchy.py',
           'src/recursive_horizons/nsc_atlas_scale_selection.py',
           'scripts/check_nsc_atlas_charge_inheritance.py')
    result['source_sha256']={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in names}
    out.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
    print('Fixed flux-step r:',result['flux_step']['r_exact'])
    print('Parke ruler:',result['parke_ruler'])
    print('Targets, not predictions:',result['target_diagnostics'])
    print('Prediction:',result['prediction'])


if __name__=='__main__':main()
