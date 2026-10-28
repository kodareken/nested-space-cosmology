"""Evaluate the fixed-parent scale ladder without using the target as a selector."""
import argparse
import hashlib
import json
from pathlib import Path

from recursive_horizons.nsc_atlas_scale_selection import calculation

ROOT=Path(__file__).resolve().parents[1]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();output=(ROOT/args.output).resolve()
    if not output.is_relative_to(ROOT) or output.exists():raise ValueError('new repository output required')
    result=calculation()
    result['source_sha256']={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in (
        'src/recursive_horizons/nsc_atlas_scale_selection.py','scripts/check_nsc_atlas_scale_selection.py',
        'tests/test_nsc_atlas_scale_selection.py')}
    result['prior_art']=['https://arxiv.org/html/hep-th/0004134#S2.SS2',
                         'https://arxiv.org/html/1804.09985','https://physics.nist.gov/cuu/Constants/Table/allascii.txt']
    output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(result))


if __name__=='__main__':main()
