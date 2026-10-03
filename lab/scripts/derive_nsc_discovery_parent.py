#!/usr/bin/env python3
"""Preview or create/check one frozen global-parent Cauchy preparation."""
from __future__ import annotations
import argparse
import json
import os
from pathlib import Path

for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[name]='1'
from recursive_horizons import nsc_discovery_parent as parent


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    modes=parser.add_mutually_exclusive_group()
    modes.add_argument('--prepare',action='store_true');modes.add_argument('--pilot',action='store_true');modes.add_argument('--check',action='store_true')
    parser.add_argument('--output',type=Path,default=parent.OUTPUT)
    parser.add_argument('--nf',type=int,default=64)
    parser.add_argument('--population',type=int,choices=parent.POPULATIONS,default=0)
    parser.add_argument('--sign',type=int,choices=parent.SIGNS,default=1)
    parser.add_argument('--k',type=float,default=None)
    parser.add_argument('--cpu-limit',type=float,default=parent.CPU_LIMIT)
    parser.add_argument('--producer-commit',default=None)
    args=parser.parse_args(argv);destination=args.output.resolve()
    if args.prepare or args.pilot:
        if destination!=parent.OUTPUT.resolve() and parent.OUTPUT.resolve() not in destination.parents:
            raise PermissionError('new parent evidence belongs under '+str(parent.OUTPUT))
        if args.pilot:
            if args.k is not None:raise ValueError('pilot computes conservative common k from both resolutions')
            result=parent.pilot(destination,population=args.population,sign=args.sign,cpu_limit=args.cpu_limit,producer_commit=args.producer_commit)
        else:result=parent.prepare(destination,nf=args.nf,population=args.population,sign=args.sign,
                k_override=args.k,cpu_limit=args.cpu_limit,producer_commit=args.producer_commit)
    elif args.check:result=parent.check(destination)
    else:result=parent.preview()
    print(json.dumps(parent._plain(result),sort_keys=True,indent=2,allow_nan=False))
    return 0


if __name__=='__main__':main()
