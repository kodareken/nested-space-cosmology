#!/usr/bin/env python3
"""Preview or prepare, lock a source-selected turn prediction, measure and check."""
from __future__ import annotations
import argparse
import json
import os
from pathlib import Path

for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[name]='1'

from recursive_horizons import nsc_discovery_homogeneous_turn as turn


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    modes=parser.add_mutually_exclusive_group()
    for mode in ('prepare','predict','measure','check'):
        modes.add_argument('--'+mode,action='store_true')
    parser.add_argument('--output',type=Path,default=turn.OUTPUT)
    args=parser.parse_args(argv)
    selected=next((mode for mode in ('prepare','predict','measure','check') if getattr(args,mode)),None)
    destination=args.output.resolve()
    if selected!='check' and selected is not None and destination!=turn.OUTPUT.resolve() and turn.OUTPUT.resolve() not in destination.parents:
        raise PermissionError('new turn records belong under '+str(turn.OUTPUT))
    result=turn.preview() if selected is None else getattr(turn,selected)(destination)
    print(json.dumps(turn._json(result),sort_keys=True,indent=2,allow_nan=False))
    return 0


if __name__=='__main__':
    main()
