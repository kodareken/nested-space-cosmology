"""Reproduce the finite physical-volume candidate without a target Lambda."""
import argparse
import hashlib
import json
from pathlib import Path

from recursive_horizons.nsc_atlas_parent_scale_selector import calculation

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    output = (ROOT/args.output).resolve()
    if not output.is_relative_to(ROOT) or output.exists():
        raise ValueError('new repository output required')
    result = calculation()
    result['source_sha256'] = {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in (
        'scripts/check_nsc_atlas_parent_scale_selector.py',
        'src/recursive_horizons/nsc_atlas_parent_scale_selector.py')}
    output.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
