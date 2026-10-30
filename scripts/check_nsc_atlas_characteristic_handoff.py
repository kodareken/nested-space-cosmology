"""Reproduce the inherited two-trace scalar continuation outside the parent horizon."""
import argparse
import hashlib
import json
from pathlib import Path
import time

from recursive_horizons.nsc_atlas_characteristic_handoff import calculation

ROOT=Path(__file__).resolve().parents[1]


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path)
    output=(ROOT/p.parse_args().output).resolve()
    if not output.is_relative_to(ROOT) or output.exists():raise ValueError('fresh repository output required')
    start=time.perf_counter();r=calculation();r['elapsed_seconds']=time.perf_counter()-start
    if r['elapsed_seconds']>180:raise RuntimeError('finite handoff budget exhausted')
    names=['src/recursive_horizons/nsc_atlas_characteristic_handoff.py',
           'src/recursive_horizons/nsc_atlas_horizon_response.py',
           'src/recursive_horizons/nsc_atlas_parent_circulation.py',
           'src/recursive_horizons/nsc_atlas_parent_scale_selector.py',
           'scripts/check_nsc_atlas_characteristic_handoff.py',
           'results/summaries/nsc-atlas-parent-circulation-v1.json',
           'results/summaries/nsc-atlas-horizon-response-v1.json']
    r['source_sha256']={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names}
    output.write_text(json.dumps(r,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(domain=r['domain'],cases=[{k:x[k] for k in ('cells_per_leg','relative_RMS_field_error','horizon_normal_derivative_relative_RMS','receiver_wrapped_phase_error','noether_flux')} for x in r['resolutions']],
        receiver=r['resolutions'][-1]['receiver'],orders=r['observed_orders'],control=r['one_trace_control'],elapsed_seconds=r['elapsed_seconds']),indent=2))


if __name__=='__main__':main()
