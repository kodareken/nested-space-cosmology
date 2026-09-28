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
    correspondence_record,
)

OUTPUT = ROOT / "results/development/nsc-local-observer-correspondence.json"


def compute():
    return correspondence_record(ROOT)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--record", action="store_true")
    modes.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.check and not OUTPUT.is_file():
        raise FileNotFoundError("local-observer correspondence record is missing; --check cannot create evidence")
    result = compute()
    if args.record:
        publish_exclusive_file(
            ROOT,
            str(OUTPUT.relative_to(ROOT)),
            (json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n").encode(),
        )
    elif result != json.loads(OUTPUT.read_text()):
        raise ValueError("local-observer correspondence replay differs")
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
