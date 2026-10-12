#!/usr/bin/env python3
"""Read saved nf128/256 T=.3/1 states for a whole-realization translation control."""
from __future__ import annotations
import argparse
import json
import os
from pathlib import Path
import sys

for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[name] = '1'

from recursive_horizons import nsc_discovery_translation as translation

OUTPUT_DIR = translation.LAB / 'results/development/nsc-discovery-translation-v1'


def destination(path):
    path = Path(path).resolve()
    if OUTPUT_DIR.resolve() not in path.parents or path.suffix != '.json':
        raise PermissionError('translation records belong under ' + str(OUTPUT_DIR))
    if path.exists():
        raise FileExistsError('refusing to replace ' + str(path))
    return path


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true', help='read-only assessment, print JSON')
    parser.add_argument('--write', type=Path, help='exclusive new JSON record')
    parser.add_argument('--shift', type=float, default=translation.DEFAULT_SHIFT)
    parser.add_argument('--step', type=float, help='optional single RK4 step, at most 0.001')
    args = parser.parse_args(argv)
    if not args.check and args.write is None:
        parser.error('choose --check or --write PATH')
    target = None if args.write is None else destination(args.write)
    record = translation.assessment(shift=args.shift, step=args.step)
    blob = json.dumps(record, sort_keys=True, indent=2, allow_nan=False) + '\n'
    if target is not None:
        if not record['ok']:
            raise RuntimeError('translation covariance checks failed; refusing production record')
        if len(blob.encode()) > 64 * 1024 * 1024:
            raise RuntimeError('translation record exceeds 64 MiB')
        target.parent.mkdir(parents=True, exist_ok=True)
        with target.open('x') as handle:
            handle.write(blob)
    sys.stdout.write(blob)
    return 0 if record['ok'] else 1


if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except (RuntimeError, ValueError, PermissionError, FileExistsError) as error:
        print(str(error), file=sys.stderr)
        raise SystemExit(1)
