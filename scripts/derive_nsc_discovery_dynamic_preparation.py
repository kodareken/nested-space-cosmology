#!/usr/bin/env python3
"""Preflight/prepare initial data or explicitly run the owned diagnostic episode."""
import argparse
import json
import os
from pathlib import Path

for name in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS","VECLIB_MAXIMUM_THREADS"):
    os.environ[name]="1"

from recursive_horizons import nsc_discovery_dynamic_preparation as preparation


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepare",action="store_true",help="explicit bounded initial-data construction")
    parser.add_argument("--coherent",type=Path,help="successor of frozen v1 including its declared rotated source")
    parser.add_argument("--episode",type=Path,help="prepare six native-basis cases from a three-source checkpoint")
    parser.add_argument("--run-episode",type=Path,help="explicit existing-runner execution of prepared cases")
    parser.add_argument("--execute",action="store_true",help="execute coherent/episode preparation rather than preflight")
    parser.add_argument("--observed-binding",type=Path,help="external post-run producer binding for preserved v1")
    parser.add_argument("--producer-commit",help="exact frozen producer revision for new records")
    parser.add_argument("--workers",type=int,default=4)
    parser.add_argument("--episode-cpu",type=float,default=300.)
    parser.add_argument("--accepted-six-hour",action="store_true")
    parser.add_argument("--max-steps",type=int,help="bounded diagnostic preview")
    parser.add_argument("--write",type=Path,help="create a new JSON+NPZ checkpoint after preparation")
    parser.add_argument("--check",type=Path,help="read-only authentication/replay of a saved checkpoint")
    parser.add_argument("--cpu-limit",type=float,default=30.)
    args=parser.parse_args(argv)
    if sum(bool(v) for v in (args.prepare,args.coherent,args.episode,args.run_episode,args.check))>1:
        parser.error("choose one mode")
    if args.check is not None:
        if args.prepare or args.write is not None:
            parser.error("check cannot prepare or write")
        report=preparation.check_record(args.check,observed_binding=args.observed_binding)
    elif args.run_episode:
        if args.write is not None:
            parser.error("run uses its already prepared directory")
        report=preparation.run_episode(args.run_episode,workers=args.workers,cpu_budget=args.episode_cpu,
                accepted_six_hour=args.accepted_six_hour,max_steps=args.max_steps)
    elif args.episode:
        if args.execute and args.write is None:
            parser.error("episode execution requires a new output directory")
        report=preparation.prepare_episode(args.episode,args.write,execute=args.execute,
            producer_commit=args.producer_commit,observed_binding=args.observed_binding)
    elif args.coherent:
        report,arrays=preparation.build_coherent_successor(args.coherent,execute=args.execute,cpu_limit=args.cpu_limit,
            producer_commit=args.producer_commit,observed_binding=args.observed_binding)
        if args.write is not None:
            if not args.execute:
                parser.error("successor write requires execute")
            report=preparation.write_record(report,arrays,args.write)
    else:
        if args.write is not None and not args.prepare:
            parser.error("write requires explicit prepare")
        if args.write is not None:
            if any(p.exists() for p in preparation.output_paths(args.write)):
                raise FileExistsError("checkpoint already exists")
        closure=None;commit=None
        if args.prepare:
            if not args.producer_commit:
                parser.error("new production preparation requires a frozen producer commit")
            closure=preparation.dependency_hashes()
            commit=preparation._git_hashes(args.producer_commit,closure)
        report,arrays=preparation.build_record(execute=args.prepare,cpu_limit=args.cpu_limit)
        if args.prepare:
            if closure!=preparation.dependency_hashes():
                raise ValueError("dependency closure changed during preparation")
            report["source_hashes"]=closure;report["producing_commit"]=commit
        if args.write is not None:
            report=preparation.write_record(report,arrays,args.write)
    print(json.dumps(report,indent=2,allow_nan=False))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
