#!/usr/bin/env python3
"""Record or replay the compact local-observer correspondence v1 artifact."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from recursive_horizons.evidence_io import publish_exclusive_file
from recursive_horizons.nsc_local_observer_correspondence import (
    ARTIFACT_ID,
    compare_records,
    correspondence_record,
)

OUTPUT = ROOT / "results/development/nsc-local-observer-correspondence.json"


def compute(*, source_context=None, source_contexts=None):
    return correspondence_record(ROOT, source_context=source_context,
                                 source_contexts=source_contexts)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--record", action="store_true")
    modes.add_argument("--check", action="store_true")
    parser.add_argument("--source-context", help="full original source commit for an explicit saved-record check")
    parser.add_argument("--source-pin", action="append", default=[], metavar="PATH=COMMIT",
                        help="explicit original source commit for one declared source path")
    args = parser.parse_args()
    if (args.source_context is not None or args.source_pin) and not args.check:
        parser.error("original source contexts are only valid with --check")
    overrides = {}
    for declaration in args.source_pin:
        path, separator, commit = declaration.partition("=")
        if not separator or not path or not commit or path in overrides:
            parser.error("--source-pin requires a unique PATH=COMMIT")
        overrides[path] = commit
    if args.check and not OUTPUT.is_file():
        raise FileNotFoundError("local-observer correspondence record is missing; --check cannot create evidence")
    result = compute(source_context=args.source_context, source_contexts=overrides)
    if args.record:
        publish_exclusive_file(
            ROOT,
            str(OUTPUT.relative_to(ROOT)),
            (json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n").encode(),
        )
    else:
        compare_records(json.loads(OUTPUT.read_text()), result)
    print(
        json.dumps(
            {
                "artifact_id": ARTIFACT_ID,
                "status": result["status"],
                "claim_status": result["claim_status"],
                "physical_claim_verified": result["physical_claim_verified"],
                "missing_common_action_evidence": result["missing_common_action_evidence"],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
