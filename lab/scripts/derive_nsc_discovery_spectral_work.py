#!/usr/bin/env python3
"""Preview or exclusively prepare, predict, run and check the prescribed pulse."""
from __future__ import annotations
import argparse
import json
import os
from pathlib import Path

for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[name]='1'

from recursive_horizons import nsc_discovery_spectral_work as work


def output_path(path):
    path=Path(path).resolve()
    if path!=work.OUTPUT.resolve() and work.OUTPUT.resolve() not in path.parents:
        raise PermissionError('new pulse records belong under '+str(work.OUTPUT))
    return path


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    modes=parser.add_mutually_exclusive_group()
    for name in ('prepare','predict','run','check'):
        modes.add_argument('--'+name,action='store_true')
    parser.add_argument('--output',type=Path,default=work.OUTPUT)
    args=parser.parse_args(argv)
    selected=next((name for name in ('prepare','predict','run','check') if getattr(args,name)),None)
    result=work.preview() if selected is None else getattr(work,selected)(output_path(args.output))
    print(json.dumps(result,sort_keys=True,indent=2,allow_nan=False))
    return 0


if __name__=='__main__':
    main()
