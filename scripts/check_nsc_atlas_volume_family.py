"""Record the analytic dimension/ancestor selection audit and finite witnesses."""
import argparse
import hashlib
import json
from pathlib import Path

from recursive_horizons.nsc_atlas_volume_family import calculation

ROOT=Path(__file__).resolve().parents[1]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    output=(ROOT/p.parse_args().output).resolve()
    if not output.is_relative_to(ROOT) or output.exists():
        raise ValueError('new repository output required')
    result=calculation()
    result['source_sha256']={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in (
        'src/recursive_horizons/nsc_atlas_volume_family.py',
        'src/recursive_horizons/nsc_atlas_parent_scale_selector.py',
        'scripts/check_nsc_atlas_volume_family.py')}
    output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
    for row in result['witnesses']:
        print(row['d'],row['ancestor_depth'],row['initial_stiffness'],row['Lambda_Lref0_squared'])
    print(result['conclusion'])


if __name__=='__main__':main()
