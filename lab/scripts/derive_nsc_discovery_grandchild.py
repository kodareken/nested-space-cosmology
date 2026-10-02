#!/usr/bin/env python3
"""Preview, lock a grandchild forecast, then independently measure its future arm."""
from __future__ import annotations
import argparse
import json
import os
from pathlib import Path

for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[name]='1'

from recursive_horizons import nsc_discovery_grandchild as grandchild


def output_path(path):
    destination=Path(path).resolve()
    if grandchild.OUTPUT.resolve() not in destination.parents:
        raise PermissionError('grandchild records belong under a new child prefix of '+str(grandchild.OUTPUT))
    return destination


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    modes=parser.add_mutually_exclusive_group()
    modes.add_argument('--prepare',action='store_true',help='propagate baseline+tangent and lock forecast before future held-out measurement')
    modes.add_argument('--run',action='store_true',help='measure the future held-out arm using an existing locked forecast')
    modes.add_argument('--check',action='store_true',help='authenticate/replay new records without evolution')
    parser.add_argument('--nf',type=int,choices=(128,256),default=256)
    parser.add_argument('--output',type=Path,help='new record prefix directory; default grandchild-v1/nf256')
    parser.add_argument('--no-matched-step',action='store_true',help='prepare only cap .0005')
    args=parser.parse_args(argv)
    if args.prepare or args.run or args.check:
        destination=output_path(args.output or grandchild.OUTPUT/f'nf{args.nf}')
        if args.prepare:report=grandchild.prepare(destination,args.nf,matched=not args.no_matched_step)
        elif args.run:report=grandchild.run(destination)
        else:report=grandchild.check(destination)
    else:report=grandchild.preview(args.nf)
    print(json.dumps(report,sort_keys=True,indent=2,allow_nan=False))
    return 0


if __name__=='__main__':
    main()
