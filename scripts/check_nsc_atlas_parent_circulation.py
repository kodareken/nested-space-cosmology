"""Reproduce smooth-parent circulation and its finite guided-light readout."""
import argparse
import hashlib
import json
from pathlib import Path
import time

from recursive_horizons.nsc_atlas_parent_circulation import calculation

ROOT=Path(__file__).resolve().parents[1]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',required=True,type=Path)
    path=(ROOT/parser.parse_args().output).resolve()
    if not path.is_relative_to(ROOT) or path.exists():
        raise ValueError('fresh repository output required')
    start=time.perf_counter()
    record=calculation()
    record['elapsed_seconds']=time.perf_counter()-start
    if record['elapsed_seconds']>180:
        raise RuntimeError('finite calculation budget exceeded')
    names=['src/recursive_horizons/nsc_atlas_parent_circulation.py',
           'src/recursive_horizons/nsc_atlas_parent_scale_selector.py',
           'src/recursive_horizons/nsc_atlas_offset_generator.py',
           'scripts/check_nsc_atlas_parent_circulation.py']
    record['source_sha256']={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in names}
    path.write_text(json.dumps(record,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({key:record[key] for key in ('controls','light_comparisons','direct_null_controls','horizon_exchange','elapsed_seconds')},indent=2))


if __name__=='__main__':main()
