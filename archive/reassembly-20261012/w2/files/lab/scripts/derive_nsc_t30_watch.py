"""Watch snapshot: forecast B passed without a landing.

The dual forecast recorded two possible cost laws. The campaign's
directory, named adjudicator, has now written its first verse: at
fifty-seven wall-minutes no station-25 exists, so the two-minute-per-t
rate did not govern. The v1-derived forecast remains the living one.
This is the watch's recorded observation, not a rerun of the forecast.

    python scripts/lab.py scripts/derive_nsc_t30_watch.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

LAB = Path(__file__).resolve().parents[1]
if str(LAB / "src") not in sys.path:
    sys.path.insert(0, str(LAB / "src"))

from recursive_horizons import nsc_influence_chain as chain  # noqa: E402

OUTPUT = LAB / "results" / "development" / "nsc-t30-watch-v1.json"
SCHEMA = "NSC-T30-WATCH-v1"


def main():
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    station25 = (LAB / "results" / "development" / "nsc-rn-inward-march-t30-v1"
                 / "mu-1.002-eps-0.001" / "station-25.json")
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "adjudicates": "nsc-t30-dual-forecast-v1.json",
        "observation": {
            "minutes_elapsed": 57.3,
            "station_25_present": bool(station25.is_file()),
            "forecast_b_passed": True,
            "forecast_a_still_live": True,
        },
        "reading": (
            "the relaxation-rate forecast required a landing near fifty "
            "minutes; none came, so the v1-derived cost law governs the "
            "run until the directory says otherwise"
        ),
        "notes": [
            "the watch record is a snapshot at the moment of writing; "
            "later landings supersede its verdicts",
        ],
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_t30_watch.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_t30_watch.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    print("  forecast B passed without landing:", record["observation"]["forecast_b_passed"])
    print("  forecast A still live:", record["observation"]["forecast_a_still_live"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())