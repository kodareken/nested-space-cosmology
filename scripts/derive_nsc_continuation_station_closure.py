"""Station-level interior closure across the recorded continuation t=10..20.

Round-14 verified the run-level identity dM = d(excision) on the whole
continuation; the family stations only reach t=18.5. This driver mines
the continuation's own half-station snapshots (t10.0000 .. t20.0000) and
evaluates the SAME identity interval by interval, plus the outer-boundary
sign per interval, across the domain the family never covered (beyond
t=18.5). Pure reading of recorded npz bytes.

    python scripts/lab.py scripts/derive_nsc_continuation_station_closure.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

LAB = Path(__file__).resolve().parents[1]
if str(LAB / "src") not in sys.path:
    sys.path.insert(0, str(LAB / "src"))

import numpy as np

from recursive_horizons import nsc_influence_chain as chain  # noqa: E402

CONTINUATION = LAB / "results" / "development" / "nsc-rn-neutral-continuation-v1"
OUTPUT = LAB / "results" / "development" / "nsc-continuation-station-closure-v1.json"
SCHEMA = "NSC-CONTINUATION-STATION-CLOSURE-v1"
STATIONS = [10.0 + 0.5 * k for k in range(21)]
CLOSURE_TOLERANCE = 1.0e-9


def scalar(value):
    value = np.asarray(value)
    return float(value.reshape(-1)[0])


def main():
    if not CONTINUATION.is_dir():
        raise SystemExit("missing continuation directory")
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    rows = []
    for station in STATIONS:
        path = CONTINUATION / "t{:0.4f}.npz".format(station)
        with np.load(path) as payload:
            rows.append({
                "t": station,
                "inner_mass": scalar(payload["coupled_inner_mass"]),
                "excision": scalar(payload["coupled_excision_flux"]),
                "outer": scalar(payload["coupled_outer_flux"]),
            })
    intervals = []
    for previous, current in zip(rows, rows[1:]):
        dm = current["inner_mass"] - previous["inner_mass"]
        de = current["excision"] - previous["excision"]
        do = current["outer"] - previous["outer"]
        intervals.append({
            "interval": [previous["t"], current["t"]],
            "delta_mass": dm,
            "delta_excision": de,
            "delta_outer": do,
            "residual": dm - de,
        })
    residuals = [abs(item["residual"]) for item in intervals]
    outer_signs = [item["delta_outer"] for item in intervals]
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "stations": len(rows),
        "intervals": intervals,
        "maximum_residual": max(residuals),
        "all_intervals_close": bool(max(residuals) <= CLOSURE_TOLERANCE),
        "outer_leaks_every_interval": bool(all(delta < 0.0 for delta in outer_signs)),
        "verdict": (
            "LATE_CLOSURE_HOLDS"
            if max(residuals) <= CLOSURE_TOLERANCE and all(delta < 0.0 for delta in outer_signs)
            else "LATE_CLOSURE_OPEN"
        ),
        "beyond_family_coverage": {
            "family_last_station": 18.5,
            "new_intervals": [item for item in intervals if item["interval"][1] > 18.5],
            "max_residual_beyond": max(
                abs(item["residual"]) for item in intervals if item["interval"][1] > 18.5
            ),
        },
        "notes": [
            "the continuation snapshots are 0.5-station resolved; the "
            "identity is the round-14 one applied interval by interval",
            "outer boundary sign per interval retests the leak-only reading",
        ],
        "input_sha256": {
            "lab/results/development/nsc-rn-neutral-continuation-v1/summary.json":
                chain.path_sha256(CONTINUATION / "summary.json"),
        },
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_continuation_station_closure.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_continuation_station_closure.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    print("  stations {}  intervals {}  max residual {:.2e}  verdict: {}".format(
        record["stations"], len(intervals), record["maximum_residual"], record["verdict"]))
    print("  outer sign per interval: {} negative of {}".format(
        sum(delta < 0.0 for delta in outer_signs), len(outer_signs)))
    beyond = record["beyond_family_coverage"]
    print("  beyond t=18.5: {} intervals, max residual {:.2e}".format(
        len(beyond["new_intervals"]), beyond["max_residual_beyond"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())