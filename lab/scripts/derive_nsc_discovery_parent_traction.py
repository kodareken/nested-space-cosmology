#!/usr/bin/env python3
"""Preview, freeze, predict and measure a finite static source-responsive collar."""
from __future__ import annotations
import argparse
import json
import os
from pathlib import Path

for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[name]='1'

from recursive_horizons import nsc_discovery_parent_traction as collar


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    modes=parser.add_mutually_exclusive_group()
    for mode in ('prepare','predict','measure','check'):
        modes.add_argument('--'+mode,action='store_true')
    parser.add_argument('--output',type=Path,default=collar.OUTPUT)
    args=parser.parse_args(argv)
    selected=next((mode for mode in ('prepare','predict','measure','check') if getattr(args,mode)),None)
    destination=args.output.resolve()
    if selected not in (None,'check') and destination!=collar.OUTPUT.resolve() and collar.OUTPUT.resolve() not in destination.parents:
        raise PermissionError('new collar records belong under '+str(collar.OUTPUT))
    report=collar.preview() if selected is None else getattr(collar,selected)(destination)
    print(json.dumps(report,sort_keys=True,indent=2,allow_nan=False))
    return 0


if __name__=='__main__':
    main()
