"""Solve three finite-amplitude S3 initial-data resolutions; retain summaries."""
import argparse
import hashlib
import json
from pathlib import Path

from recursive_horizons.nsc_atlas_finite_parent_data import calculation

ROOT=Path(__file__).resolve().parents[1]


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path)
    out=(ROOT/p.parse_args().output).resolve()
    if not out.is_relative_to(ROOT) or out.exists():raise ValueError('fresh repository output required')
    r=calculation()
    names=['src/recursive_horizons/nsc_atlas_finite_parent_data.py',
        'src/recursive_horizons/nsc_atlas_s3_spectral.py',
        'src/recursive_horizons/nsc_atlas_parent_circulation.py',
        'src/recursive_horizons/nsc_atlas_parent_scale_selector.py',
        'scripts/check_nsc_atlas_finite_parent_data.py']
    r['source_sha256']={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names}
    out.write_text(json.dumps(r,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(levels=r['levels'],convergence=r['convergence'],compute=r['compute']),indent=2))


if __name__=='__main__':main()
