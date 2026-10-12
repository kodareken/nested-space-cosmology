"""Contention correction: the v1 marches shared the machine.

Forecast A's basis described the v1 light-case rate as measured. The
v1 mtimes show all three cases prepared within seconds of each other
(03:30:31-36) and landing at 04:24, 05:00 and 06:24 - three concurrent
marches sharing the machine. The light case's 173.7-minute duration is
therefore contention-inflated, and the solo t30 run should beat
forecast A substantially. The dual forecast remains the committed
bracket, and the campaign write times remain the adjudicator.

    python scripts/lab.py scripts/derive_nsc_t30_contention_correction.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

LAB = Path(__file__).resolve().parents[1]
if str(LAB / "src") not in sys.path:
    sys.path.insert(0, str(LAB / "src"))

from recursive_horizons import nsc_influence_chain as chain  # noqa: E402

OUTPUT = LAB / "results" / "development" / "nsc-t30-contention-correction-v1.json"
SCHEMA = "NSC-T30-CONTENTION-CORRECTION-v1"


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
        "corrects": ["nsc-t30-dual-forecast-v1.json (basis note)",
                     "nsc-t30-timing-resolution-v1.json (252-minute forecast)"],
        "evidence": (
            "three v1 prepares within five seconds of each other "
            "(03:30:31-36) and landings at 04:24, 05:00, 06:24"
        ),
        "corrected_reading": (
            "forecast A is contention-inflated: the light case shared "
            "the machine with two sibling marches. The solo t30 run "
            "should land earlier than A says, somewhere inside the "
            "committed bracket; the adjudicator is unchanged: the "
            "campaign's own write times."
        ),
        "current_state": {
            "wall_minutes": 77.2,
            "station_25_present": bool(station25.is_file()),
        },
        "notes": [
            "the bracket record and the watch record remain the papers "
            "of record; this correction annotates their basis",
        ],
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_t30_contention_correction.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_t30_contention_correction.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    print("  corrected:", record["corrected_reading"][:60] + "...")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())