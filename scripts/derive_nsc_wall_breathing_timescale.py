"""Breathing timescale: is the wall rate slaved to the infall clock?

Round-24 left three monotone trends across mu. This driver tests one
parameter-free prediction: the wall's radial rise/fall rates in horizon
units should scale as 1/T_total(mu) if the breathing is set by the
infall clock of each chart (instrument B of the round-12/13 atlas).
It computes T_total on the vacuum chart of each mu and forms the
products rate*T_total; constancy across the family is the test.

    python scripts/lab.py scripts/derive_nsc_wall_breathing_timescale.py
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
from recursive_horizons import nsc_rn_source as source  # noqa: E402

BREATHING = LAB / "results" / "development" / "nsc-wall-breathing-family-v1.json"
OUTPUT = LAB / "results" / "development" / "nsc-wall-breathing-timescale-v1.json"
SCHEMA = "NSC-WALL-BREATHING-TIMESCALE-v1"
MASSES = {"mu-1.002-eps-0.001": 1.002, "mu-1.010-eps-0.001": 1.010,
          "mu-1.040-eps-0.001": 1.040}


def infall_time(mass_over_rm):
    model = source.action_reference(mass_over_rm)
    r_plus = float(model.horizons["r_plus"])
    radius = np.linspace(r_plus, 32.0, 4000)
    f = 1.0 - 2.0 * float(model.mass) / radius + float(model.magnetic_r2) / radius ** 2
    weight = 1.0 / (1.0 + np.sqrt(np.maximum(0.0, f)))
    return float(np.trapezoid(weight, radius)), r_plus


def spread(values):
    mean = sum(values) / len(values)
    return (sum((v - mean) ** 2 for v in values) / len(values)) ** 0.5 / abs(mean)


def main():
    if not BREATHING.is_file():
        raise SystemExit("missing breathing record: " + str(BREATHING))
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    breathing = json.loads(BREATHING.read_text(encoding="utf-8"))
    cases = {}
    for name, mass in MASSES.items():
        t_total, r_plus = infall_time(mass)
        entry = breathing["cases"][name]
        rise = entry["mean_rise_rate"]
        fall = entry["mean_fall_rate"]
        cases[name] = {
            "mass_over_rm": mass,
            "r_plus": r_plus,
            "T_total": t_total,
            "rise_rate": rise,
            "fall_rate": fall,
            "rise_times_T": None if rise is None else rise * t_total,
            "fall_times_T": None if fall is None else fall * t_total,
            "rise_per_horizon_free_fall": None if rise is None else rise * (t_total / r_plus),
        }
    rise_products = [entry["rise_times_T"] for entry in cases.values()
                     if entry["rise_times_T"] is not None]
    fall_products = [entry["fall_times_T"] for entry in cases.values()
                     if entry["fall_times_T"] is not None]
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "cases": cases,
        "rise_times_T_relative_spread": spread(rise_products),
        "fall_times_T_relative_spread": spread(fall_products),
        "reading": (
            "if the breathing is infall-clocked, rate*T_total is a family "
            "constant; the measured spread quantifies the departure"
        ),
        "notes": [
            "T_total is the vacuum-chart infall time of instrument B, "
            "omega = 0 quadrature from r_+ to 32 r_m",
            "parameter-free prediction; nothing is fitted",
        ],
        "input_sha256": {
            "lab/results/development/nsc-wall-breathing-family-v1.json":
                chain.path_sha256(BREATHING),
        },
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_wall_breathing_timescale.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_wall_breathing_timescale.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    for name, entry in cases.items():
        print(
            "  {:<20} T_total {: .4f} rise_tT {: .4f} fall_tT {: .4f}".format(
                name, entry["T_total"], entry["rise_times_T"], entry["fall_times_T"],
            )
        )
    print("  rise*T spread: {:.2%}   fall*T spread: {:.2%}".format(
        record["rise_times_T_relative_spread"], record["fall_times_T_relative_spread"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())