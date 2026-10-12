"""Timing of the trap-margin response to the mass feed.

Round-64 recorded the trap's outward drift correlated 0.835 with the
same-interval mass step over 40 intervals. This driver completes the
mechanism with the lag profile - drift(i) against mass(i+k) for
k = -1, 0, +1 - and the simultaneous regression slope (trap displacement
per unit mass, dimensionless). The wall's balance responds one station
late; the trap margin's own timing is the question.

    python scripts/lab.py scripts/derive_nsc_trap_timing.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

LAB = Path(__file__).resolve().parents[1]
if str(LAB / "src") not in sys.path:
    sys.path.insert(0, str(LAB / "src"))

from recursive_horizons import nsc_influence_chain as chain  # noqa: E402

SCALARS = LAB / "results" / "development" / "nsc-rn-neutral-continuation-v1" / "scalars.jsonl"
OUTPUT = LAB / "results" / "development" / "nsc-trap-timing-v1.json"
SCHEMA = "NSC-TRAP-TIMING-v1"


def pearson(x, y):
    n = len(x)
    xm = sum(x) / n
    ym = sum(y) / n
    sxy = sum((a - xm) * (b - ym) for a, b in zip(x, y))
    sxx = sum((a - xm) ** 2 for a in x)
    syy = sum((b - ym) ** 2 for b in y)
    return None if (sxx * syy) == 0.0 else sxy / (sxx * syy) ** 0.5


def main():
    if not SCALARS.is_file():
        raise SystemExit("missing scalars jsonl")
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    lines = SCALARS.read_text(encoding="utf-8").splitlines()
    rows = []
    for line in lines:
        entry = json.loads(line)
        rows.append({
            "t": float(entry["coordinate_t_over_rm"]),
            "trap": float(entry["trapping_horizon"]),
            "mass": float(entry["coupled_inner_mass"]),
        })
    rows.sort(key=lambda row: row["t"])
    drifts = []
    masses = []
    for previous, current in zip(rows, rows[1:]):
        dt = current["t"] - previous["t"]
        if dt <= 0.0:
            continue
        drifts.append((current["trap"] - previous["trap"]) / dt)
        masses.append(current["mass"] - previous["mass"])
    lags = {}
    for lag in (-1, 0, 1):
        x = []
        y = []
        for index in range(len(drifts)):
            bound = index + lag
            if 0 <= bound < len(masses):
                x.append(drifts[index])
                y.append(masses[bound])
        lags[str(lag)] = {
            "correlation": pearson(x, y),
            "pairs": len(x),
        }
    xm = sum(masses) / len(masses)
    ym = sum(drifts) / len(drifts)
    sxy = sum((a - xm) * (b - ym) for a, b in zip(masses, drifts))
    sxx = sum((a - xm) ** 2 for a in masses)
    slope = sxy / sxx if sxx else None
    peak = max((-1, 0, 1), key=lambda lag: abs(lags[str(lag)]["correlation"]))
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "lag_profile": lags,
        "peak_lag": peak,
        "trap_displacement_per_unit_mass": slope,
        "verdict": (
            "TRAP_RESPONDS_SIMULTANEOUSLY" if peak == 0 else "TRAP_LAGGED"
        ),
        "reading": (
            "the trap margin moves with the same-interval feed if the "
            "profile peaks at lag 0; the wall's balance lags one station"
        ),
        "notes": [
            "interval cadence is the ledger's own half-station split",
            "slope units: trap drift in r_m per t over mass in r_m per t",
        ],
        "input_sha256": {
            "lab/results/development/nsc-rn-neutral-continuation-v1/scalars.jsonl":
                chain.path_sha256(SCALARS),
        },
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_trap_timing.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_trap_timing.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    for lag in (-1, 0, 1):
        print("  lag {:+d}: corr {!s} (n={})".format(
            lag, None if lags[str(lag)]["correlation"] is None
            else "{:.3f}".format(lags[str(lag)]["correlation"]),
            lags[str(lag)]["pairs"]))
    print("  slope {:.4f}  peak lag {:+d}  verdict: {}".format(
        slope, peak, record["verdict"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())