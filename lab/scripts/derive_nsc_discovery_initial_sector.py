#!/usr/bin/env python3
"""Explicit preparation/episode commands for the bounded initial-sector controls."""
import argparse
import json
import os
from pathlib import Path

for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[name] = "1"

from recursive_horizons import nsc_discovery_initial_sector as sector


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--prepare", action="store_true", help="construct three new initial sectors; no evolution")
    modes.add_argument("--episode", type=Path, help="create six exact cap forks from the prepared sector record")
    modes.add_argument("--run-episode", type=Path, help="explicit existing-pipeline scientific execution")
    modes.add_argument("--check", type=Path, help="read-only prepared-record replay")
    modes.add_argument("--assess-episode", type=Path, help="read saved episode rows; no evolution")
    parser.add_argument("--initial", type=Path, default=sector.PREPARED)
    parser.add_argument("--balanced", action="store_true", help="conditional cubic chi-star and fixed +/-1 percent balanced-v2 controls")
    parser.add_argument("--comparison", type=Path, default=sector.COMPARISON)
    parser.add_argument("--producer-commit", help="frozen revision containing the complete new source closure")
    parser.add_argument("--execute", action="store_true", help="write episode preparation rather than preflight")
    parser.add_argument("--write", type=Path, help="creation-only preparation stem or episode directory")
    parser.add_argument("--cpu-limit", type=float, default=30.)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--episode-cpu", type=float, default=300.)
    parser.add_argument("--max-steps", type=int, help="explicit bounded preview; no implied trajectory reset")
    args = parser.parse_args(argv)
    if args.run_episode:
        if args.write:
            parser.error("run uses its already prepared successor directory")
        report = sector.run_episode(args.run_episode, workers=args.workers,
                                    cpu_budget=args.episode_cpu, max_steps=args.max_steps)
    elif args.check:
        if args.write:
            parser.error("check is read-only")
        report = sector.check_record(args.check)
    elif args.assess_episode:
        if args.write:
            parser.error("assessment is read-only")
        report = sector.assess_episode(args.assess_episode)
    elif args.episode:
        if args.execute and args.write is None:
            parser.error("episode creation requires a new output directory")
        report = sector.prepare_episode(args.episode, args.write,
                                        execute=args.execute, producer_commit=args.producer_commit)
    else:
        if args.write and not args.prepare:
            parser.error("preparation write requires --prepare")
        if args.prepare and not args.producer_commit:
            parser.error("production preparation requires its frozen producing commit")
        report, arrays = sector.build_record(args.initial, execute=args.prepare,
            cpu_limit=args.cpu_limit, producer_commit=args.producer_commit, comparison=args.comparison,
            balanced=args.balanced)
        if args.write:
            report = sector.write_record(report, arrays, args.write)
    print(json.dumps(episode_jsonable(report), indent=2, allow_nan=False))
    return 0


def episode_jsonable(value):
    from recursive_horizons.nsc_discovery_episode import jsonable
    return jsonable(value)


if __name__ == "__main__":
    raise SystemExit(main())
