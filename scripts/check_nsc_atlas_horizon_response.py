"""Produce scalar-sourced horizon area, entropy, tail bounds and flux directions."""
import argparse
import hashlib
import json
from pathlib import Path
import time

from recursive_horizons.nsc_atlas_horizon_response import calculation

ROOT=Path(__file__).resolve().parents[1]


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path)
    output=(ROOT/p.parse_args().output).resolve()
    if not output.is_relative_to(ROOT) or output.exists():raise ValueError('fresh repository output required')
    start=time.perf_counter();r=calculation();r['elapsed_seconds']=time.perf_counter()-start
    if r['elapsed_seconds']>180:raise RuntimeError('finite horizon-response budget exhausted')
    names=['src/recursive_horizons/nsc_atlas_horizon_response.py',
           'src/recursive_horizons/nsc_atlas_parent_circulation.py',
           'src/recursive_horizons/nsc_atlas_parent_scale_selector.py',
           'scripts/check_nsc_atlas_horizon_response.py',
           'results/summaries/nsc-atlas-parent-circulation-v1.json']
    r['source_sha256']={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names}
    output.write_text(json.dumps(r,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({k:r[k] for k in ('exact','area_history','angular_transfer','controls','infinite_future_tail','elapsed_seconds')},indent=2))


if __name__=='__main__':main()
