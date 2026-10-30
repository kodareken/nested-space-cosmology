"""Produce the compact constrained metric/cubic-matter response record."""
import argparse
import hashlib
import json
from pathlib import Path

from recursive_horizons.nsc_atlas_parent_backreaction import calculation

ROOT=Path(__file__).resolve().parents[1]


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    output=(ROOT/p.parse_args().output).resolve()
    if not output.is_relative_to(ROOT) or output.exists():raise ValueError('fresh repository output required')
    r=calculation()
    names=['src/recursive_horizons/nsc_atlas_parent_backreaction.py',
        'src/recursive_horizons/nsc_atlas_cubic_operators.py',
        'src/recursive_horizons/nsc_atlas_volume_feedback.py',
        'src/recursive_horizons/nsc_atlas_parent_circulation.py',
        'src/recursive_horizons/nsc_atlas_horizon_response.py',
        'src/recursive_horizons/nsc_atlas_parent_scale_selector.py',
        'src/recursive_horizons/nsc_atlas_s3_polynomials.py',
        'scripts/check_nsc_atlas_parent_backreaction.py']
    r['source_sha256']={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names}
    output.write_text(json.dumps(r,sort_keys=True,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(controls=r['numerical_controls'],energy=r['energy_balances'],
        observer=r['observer'],cross_observable=r['cross_observable'],
        harmonic_generation=r['harmonic_generation'],elapsed=r['elapsed_seconds']),indent=2))


if __name__=='__main__':main()
