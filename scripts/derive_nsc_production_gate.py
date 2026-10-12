"""The freeze gate's verdict on the t=40 launch attempt.

Round-70's preflight passed; round-71 executed the ceremony for real:
calibration written, prepare attempted against the frozen root commit.
The runner refused with 'producer bytes are not the frozen root commit'.
This record documents the gate outcome and its diagnosis: the working
tree's producer nsc_rn_pg.py differs from the frozen commit's blob - a
pre-existing dirty-tree condition, not an edit made by this programme.
The freeze gate did its job; the launch now awaits the root owner's
tree reconciliation.

    python scripts/lab.py scripts/derive_nsc_production_gate.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

LAB = Path(__file__).resolve().parents[1]
if str(LAB / "src") not in sys.path:
    sys.path.insert(0, str(LAB / "src"))

from recursive_horizons import nsc_influence_chain as chain  # noqa: E402

OUTPUT = LAB / "results" / "development" / "nsc-production-gate-v1.json"
SCHEMA = "NSC-PRODUCTION-GATE-v1"


def main():
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "attempted": {
            "calibration_written": True,
            "calibration_path": "nsc-rn-calibration-t40-v1.json",
            "prepare_refused": True,
            "refusal_text": "producer bytes are not the frozen root commit",
            "no_directory_created": True,
        },
        "diagnosis": {
            "modified_producer": "lab/src/recursive_horizons/nsc_rn_pg.py",
            "modification_pre_existing": True,
            "tree_was_dirty_before_this_programme": True,
            "no_producer_edit_made_here": True,
        },
        "meaning": (
            "the freeze gate is working as designed: no production run on "
            "a tree whose producer bytes differ from the frozen commit. "
            "Launching the priced t=40 campaign is now a root-owner act "
            "(tree reconciliation), not a discovery-side task."
        ),
        "notes": [
            "the calibration file remains as a legitimate dry artifact",
            "the preflight record (round 70) stays valid: the price and "
            "ceremony are unchanged; the gate is the missing step",
        ],
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_production_gate.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_production_gate.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    print("  gate:", record["attempted"]["refusal_text"])
    print("  diagnosis:", record["diagnosis"]["modified_producer"], "(pre-existing)")
    print("  meaning:", record["meaning"][:110])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())