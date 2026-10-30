"""Retain the contrast/Planck/horizon calculation without selecting an SI input."""
import argparse
import hashlib
import json
from pathlib import Path
from recursive_horizons.nsc_atlas_contrast_hierarchy import calculation
ROOT=Path(__file__).resolve().parents[1]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    out=(ROOT/p.parse_args().output).resolve()
    if not out.is_relative_to(ROOT) or out.exists():raise ValueError('new repository output required')
    r=calculation()
    names=('src/recursive_horizons/nsc_atlas_contrast_hierarchy.py',
           'src/recursive_horizons/nsc_atlas_parent_scale_selector.py',
           'src/recursive_horizons/nsc_atlas_parke_nucleation.py',
           'src/recursive_horizons/nsc_atlas_si_hierarchy.py',
           'src/recursive_horizons/nsc_atlas_scale_selection.py',
           'scripts/check_nsc_atlas_contrast_hierarchy.py')
    r['source_sha256']={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names}
    out.write_text(json.dumps(r,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
    print('Observed same-definition benchmark:',r['observed_benchmark'])
    print('Parke controls:',r['auxiliary_Parke_snapshot_chain'])
    for row in r['fixed_coupling_j2_sectors']:
        print('Fixed-coupling sector:',row)
    print('Prediction:',r['prediction'])


if __name__=='__main__':main()
