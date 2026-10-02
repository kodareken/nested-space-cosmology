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
    modes.add_argument('--confirm',action='store_true',help='one bounded half-step confirmation of an immutable reference forecast')
    parser.add_argument('--nf',type=int,choices=(128,256),default=256)
    parser.add_argument('--output',type=Path,help='new record prefix directory; default grandchild-v1/nf256')
    parser.add_argument('--no-matched-step',action='store_true',help='prepare only cap .0005')
    parser.add_argument('--proper-increment',type=float,default=None,help='future normal-clock increment in (0,1]; default .05; use .75 in a new successor prefix')
    parser.add_argument('--reference',type=Path,help='with --confirm, completed original forecast/measurement directory')
    parser.add_argument('--spatial',action='store_true',help='with --confirm, optionally add nf128 at cap .0005 if measured admission fits')
    args=parser.parse_args(argv)
    if args.proper_increment is not None and (args.run or args.check or args.confirm):
        parser.error('--run/--check use the already locked proper increment; set --proper-increment during --prepare or preview')
    increment=grandchild.validate_proper_increment(grandchild.DELTA_TAU if args.proper_increment is None else args.proper_increment)
    if args.reference is not None and not args.confirm:parser.error('--reference belongs to --confirm')
    if args.spatial and not args.confirm:parser.error('--spatial belongs to --confirm')
    if args.confirm and (args.reference is None or args.output is None):parser.error('--confirm requires explicit --reference and new --output')
    if args.confirm and (args.nf!=256 or args.no_matched_step):parser.error('--confirm has one fixed nf256 half-step pair; use --spatial for optional nf128')
    if args.prepare or args.run or args.check or args.confirm:
        destination=output_path(args.output or grandchild.OUTPUT/f'nf{args.nf}')
        if args.prepare:report=grandchild.prepare(destination,args.nf,proper_increment=increment,matched=not args.no_matched_step)
        elif args.run:report=grandchild.run(destination)
        elif args.confirm:report=grandchild.confirm(args.reference,destination,spatial=args.spatial)
        else:report=grandchild.check(destination)
    else:report=grandchild.preview(args.nf,proper_increment=increment)
    print(json.dumps(report,sort_keys=True,indent=2,allow_nan=False))
    return 0


if __name__=='__main__':
    main()
