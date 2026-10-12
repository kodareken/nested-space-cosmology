#!/usr/bin/env python3
"""Preview, prepare, predict or run the bounded parent incoming-cut response."""
import os
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[name]='1'
import argparse
import json
from pathlib import Path
from recursive_horizons import nsc_discovery_parent_cut_response as cut


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__);modes=parser.add_mutually_exclusive_group()
    for name in ('prepare','predict','run'):modes.add_argument('--'+name,action='store_true')
    parser.add_argument('--output',type=Path,default=cut.OUTPUT)
    parser.add_argument('--station',type=float,choices=cut.STATIONS,default=2.25)
    parser.add_argument('--amplitude',type=float,default=cut.ALPHA)
    args=parser.parse_args(argv);mode=next((n for n in ('prepare','predict','run') if getattr(args,n)),None)
    output=args.output.resolve()
    if mode and output!=cut.OUTPUT.resolve() and cut.OUTPUT.resolve() not in output.parents:raise PermissionError('new response stages require the designated fresh prefix')
    if mode=='run':result=cut.run(output,station=args.station,amplitude=args.amplitude)
    else:result=cut.preview() if mode is None else getattr(cut,mode)(output)
    print(json.dumps(result,indent=2,sort_keys=True,allow_nan=False))


if __name__=='__main__':main()
