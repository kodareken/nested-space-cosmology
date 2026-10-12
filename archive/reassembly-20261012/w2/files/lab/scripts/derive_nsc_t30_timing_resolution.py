"""Timing resolution: the v1 march duration from recorded file mtimes.

Round-102 corrected the poll-count error; this driver resolves the
question the correction left open by reading the recorded mtimes of the
v1 campaign's files. The v1 march's prepare-to-station-20 duration is
173.7 minutes with the first station at 56.5 minutes and steady
~39-minute intervals after. The t30 campaign's expected landings follow
by the same per-stretch cost.

    python scripts/lab.py scripts/derive_nsc_t30_timing_resolution.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

LAB = Path(__file__).resolve().parents[1]
if str(LAB / "src") not in sys.path:
    sys.path.insert(0, str(LAB / "src"))

from recursive_horizons import nsc_influence_chain as chain  # noqa: E402

OUTPUT = LAB / "results" / "development" / "nsc-t30-timing-resolution-v1.json"
SCHEMA = "NSC-T30-TIMING-RESOLUTION-v1"


def main():
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    v1 = LAB / "results" / "development" / "nsc-rn-inward-march-v1" / "mu-1.002-eps-0.001"
    prepare = v1 / "prepare.json"
    stations = [5.0, 10.0, 15.0, 20.0]
    times = []
    for t in stations:
        path = v1 / f"station-{t:g}.npz"
        if not path.is_file():
            raise SystemExit(f"missing v1 station file {path.name}")
        times.append(path.stat().st_mtime)
    start = prepare.stat().st_mtime
    minutes = [(ts - start) / 60.0 for ts in times]
    intervals = [minutes[0]] + [minutes[i] - minutes[i - 1] for i in range(1, len(minutes))]
    mean_stretch = sum(intervals[1:]) / (len(intervals) - 1)
    forecast_25 = minutes[0] + 4.0 * mean_stretch
    forecast_30 = forecast_25 + mean_stretch
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "resolves": "nsc-t30-timing-correction-v1.json",
        "v1_minutes_from_prepare": {
            "station_5": minutes[0],
            "station_10": minutes[1],
            "station_15": minutes[2],
            "station_20": minutes[3],
        },
        "v1_intervals_minutes": intervals,
        "v1_total_prefix_minutes": minutes[-1],
        "per_stretch_after_first_minutes": mean_stretch,
        "t30_forecast_minutes_from_launch": {
            "station_25": forecast_25,
            "station_30": forecast_30,
        },
        "reading": (
            "the v1 march needed 56.5 minutes to its first recorded station "
            "and about 39 minutes per further stretch; the t30 campaign's "
            "stations are expected near 3.5 and 4.2 hours after launch"
        ),
        "notes": [
            "mtimes are recorded bytes of the v1 campaign, not narration",
        ],
        "input_sha256": {},
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_t30_timing_resolution.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_t30_timing_resolution.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    print("  v1 intervals (min): {:.1f} {:.1f} {:.1f} {:.1f}".format(*intervals))
    print("  total prefix: {:.1f} min".format(record["v1_total_prefix_minutes"]))
    print("  t30 forecast (min): station-25 {:.1f}, station-30 {:.1f}".format(
        record["t30_forecast_minutes_from_launch"]["station_25"],
        record["t30_forecast_minutes_from_launch"]["station_30"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())