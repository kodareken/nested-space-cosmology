"""Evaluate the finite backward bootstrap and two unchanged scalar clocks."""
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
import hashlib
import json
from pathlib import Path

from recursive_horizons.nsc_atlas_parent_past import past_bound,solve_past

ROOT=Path(__file__).resolve().parents[1]


def worker(clock):return solve_past(clock=clock)


def main():
    out=ROOT/'results/summaries/nsc-atlas-parent-past-v1.json'
    if out.exists():raise ValueError('refuse to replace retained proof record')
    with ProcessPoolExecutor(max_workers=2) as pool:
        controls=list(pool.map(worker,('proper','efolds')))
    b=past_bound()
    result=dict(schema='nsc-atlas-parent-past-v1',
                rational_bounds={k:str(v) for k,v in b.items()},
                approximate_values={k:float(v) for k,v in b.items()},
                unaveraged_controls=controls,
                clock_endpoint_difference=abs(controls[0]['endpoint_density_ratio']-controls[1]['endpoint_density_ratio']),
                scope='Finite homogeneous Cauchy continuation to A=1/3, not formation, global wall completion or initial-state selection.',
                source_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in (
                    'src/recursive_horizons/nsc_atlas_parent_past.py',
                    'src/recursive_horizons/nsc_atlas_inherited_matter.py',
                    'scripts/derive_nsc_atlas_parent_past.py')})
    out.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(result))


if __name__=='__main__':main()
