#!/usr/bin/env python3
"""Resume or replay actual18 direct-vacuum rows, group14/low32_1, pi/2<=E<2.

The nine positive rows are 32 through 40, each with its original negative
partner. Pass an explicit row budget and CPU budget. A partial run stays
OPEN. Missing rows are identities with null bounds. Finishing this finite
window does not close the upstream budget or the local gate, and it does
not claim the frozen 2<=E<16 direct window.
"""
import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from nsc_threshold18_source import BASE_COMMIT, cover


DEFAULT_OUTPUT = "results/development/nsc-threshold18-actual18-v1"


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
    parser.add_argument("--report",
                        help="optional JSON report path; omitted by tests that only parse status")
    args = parser.parse_args(argv)
    output = Path(args.output)
    if not output.is_absolute():
        output = ROOT / output
    if args.check:
        if args.row_budget is not None or args.cpu_budget is not None:
            parser.error("check replays saved witnesses and does not accept a solve budget")
        result = cover(output, mode="check")
    else:
        if args.row_budget is None or args.cpu_budget is None:
            parser.error("resume requires explicit --row-budget and --cpu-budget")
        result = cover(output, mode="resume", row_budget=args.row_budget,
                       cpu_budget=args.cpu_budget)
    result["retrieval"] = {
        "base_commit": BASE_COMMIT,
        "repository": "https://github.com/kodareken/nested-space-cosmology",
        "myrsa": "unavailable",
        "myrsa_used": False,
        "stale_nuc_checkout_used": False,
        "authoritative_transport": "github-pin",
    }
    text = json.dumps(result, sort_keys=True, indent=2, allow_nan=False)
    print(text)
    if args.report:
        report = Path(args.report)
        if not report.is_absolute():
            report = ROOT / report
        report.parent.mkdir(parents=True, exist_ok=True)
        report.write_text(text + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
