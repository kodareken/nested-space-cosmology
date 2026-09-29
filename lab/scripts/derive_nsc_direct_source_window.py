#!/usr/bin/env python3
"""Resume or replay direct-vacuum bounds for group14 angular+1, 2<=E<16.

Pass an explicit row budget and CPU budget. This command does not launch the
87-row window unless those budgets ask for it. A partial run stays OPEN, and
finishing this finite window still does not close the upstream budget or the
local gate.
"""
import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src"))
from recursive_horizons.nsc_direct_source_window import cover


DEFAULT_OUTPUT = "results/development/nsc-direct-source-window-v1"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--resume", action="store_true",
                      help="validate saved rows, then capture at most --row-budget new rows")
    mode.add_argument("--check", action="store_true",
                      help="replay saved traces and do not solve")
    parser.add_argument("--row-budget", type=int,
                        help="maximum number of new positive rows to solve")
    parser.add_argument("--cpu-budget", type=float,
                        help="process-time budget in seconds for new solves")
    parser.add_argument("--output", default=DEFAULT_OUTPUT,
                        help="checkpoint directory, created only by --resume")
    args = parser.parse_args(argv)
    output = Path(args.output)
    if not output.is_absolute():
        output = ROOT/output
    if args.check:
        if args.row_budget is not None or args.cpu_budget is not None:
            parser.error("check replays saved witnesses and does not accept a solve budget")
        result = cover(ROOT, output, mode="check")
    else:
        if args.row_budget is None or args.cpu_budget is None:
            parser.error("resume requires explicit --row-budget and --cpu-budget")
        result = cover(ROOT, output, mode="resume", row_budget=args.row_budget,
                       cpu_budget=args.cpu_budget)
    print(json.dumps(result, sort_keys=True, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
