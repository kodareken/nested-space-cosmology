"""Dual forecast: two recorded cost laws disagree on the t30 arrival.

The v1 march mtimes imply ~8.7 minutes per t for mu=1.002 (173.7 min
to t=20). The relaxation campaign's mtimes show a steady 2.0 minutes
per t across all six stations to t=80. Both are recorded bytes. This
record lays the two forecasts side by side and names the campaign
directory as the adjudicator - the arrival time itself will say which
law governed the run.

    python scripts/lab.py scripts/derive_nsc_t30_dual_forecast.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

LAB = Path(__file__).resolve().parents[1]
if str(LAB / "src") not in sys.path:
    sys.path.insert(0, str(LAB / "src"))

from recursive_horizons import nsc_influence_chain as chain  # noqa: E402

OUTPUT = LAB / "results" / "development" / "nsc-t30-dual-forecast-v1.json"
SCHEMA = "NSC-T30-DUAL-FORECAST-v1"


def main():
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "forecast_a_v1_mtimes": {
            "basis": "nsc-march-duration-table-v1.json (mu=1.002)",
            "minutes_per_t": 173.7 / 20.0,
            "station_25_minutes": 212.8,
            "station_30_minutes": 251.9,
        },
        "forecast_b_relaxation_mtimes": {
            "basis": "nsc-relaxation-duration-v1.json (steady 2.0 min/t)",
            "minutes_per_t": 2.0,
            "station_25_minutes": 50.0,
            "station_30_minutes": 60.0,
        },
        "adjudicator": "nsc-rn-inward-march-t30-v1 station write times",
        "current_watch": {
            "minutes_so_far": 55.3,
            "status": (
                "past forecast-B's station-25, well inside forecast-A - "
                "the arrival will pick the governing law"
            ),
        },
        "notes": [
            "the v1 light-case machine was historically contended; the "
            "relaxation campaign's steady rate is the cleanest recorded "
            "march rate",
        ],
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_t30_dual_forecast.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_t30_dual_forecast.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    print("  A (v1 mtimes): station-25 {:.1f} min".format(
        record["forecast_a_v1_mtimes"]["station_25_minutes"]))
    print("  B (relaxation): station-25 {:.0f} min".format(
        record["forecast_b_relaxation_mtimes"]["station_25_minutes"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())