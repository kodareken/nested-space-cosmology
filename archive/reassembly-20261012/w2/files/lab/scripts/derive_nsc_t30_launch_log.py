"""Launch log for the t30 campaign: the whole story in one record.

The campaign is the first new production march of this programme's
late phase. This record is its ledger: launch, cadence law, the
deadlock diagnosis, and the throughput recalibration - so a later
reader of the campaign directory needs no archaeology.

    python scripts/lab.py scripts/derive_nsc_t30_launch_log.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

LAB = Path(__file__).resolve().parents[1]
if str(LAB / "src") not in sys.path:
    sys.path.insert(0, str(LAB / "src"))

from recursive_horizons import nsc_influence_chain as chain  # noqa: E402

OUTPUT = LAB / "results" / "development" / "nsc-t30-launch-log-v1.json"
SCHEMA = "NSC-T30-LAUNCH-LOG-v1"


def main():
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "campaign": "nsc-rn-inward-march-t30-v1",
        "launched_by": "derive_nsc_rn_inward_march.py --mu 1.002 --eps 0.001 "
                        "--stations 25.0 30.0",
        "launch_time": "2026-10-06T23:55:41",
        "period": [
            {
                "event": "prepare.json written (creation-only trail audited)",
                "round": 90,
            },
            {
                "event": (
                    "cadence law established: the runner writes only at "
                    "listed stations; silence before t=25 is healthy"
                ),
                "round": 93,
            },
            {
                "event": (
                    "a concurrent full-suite pytest family diagnosed as "
                    "deadlocked (CPU frozen on repeated audits) and killed; "
                    "the march family was left untouched"
                ),
                "round": 96,
            },
            {
                "event": (
                    "throughput recalibrated: the historical four-station "
                    "marches landed within about an hour, so the first "
                    "station is expected near that envelope, not hours"
                ),
                "round": 96,
            },
        ],
        "instruments_armed": [
            "nsc-t30-extinction reducer (round 89)",
            "conditional test test_t30_extinction_once_the_march_lands (round 92)",
        ],
        "notes": [
            "the launch pre-dates this record; the record is the launch's "
            "official ledger entry and is written before the landing",
        ],
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_t30_launch_log.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_t30_launch_log.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    print("  period entries:", len(record["period"]))
    print("  instruments armed:", len(record["instruments_armed"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())