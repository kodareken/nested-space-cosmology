"""Fall-phase clock: quarter-power gap law versus half-power.

Round-26 excluded the half-power gap clock for the fall phase
(fall*sqrt(gap) spread unchanged at 15.6%). The pre-registered follow-up
is the quarter-power law: rise follows gap^-1/2, the natural partner is
fall following gap^-1/4. This driver evaluates fall*gap^1/4 across the
three cases and holds its spread against the fall*gap^1/2 spread; the
quarter-power exponent is declared BEFORE the numbers are formed. The
fitted exponent is also reported and flagged as post-hoc.

    python scripts/lab.py scripts/derive_nsc_wall_clock_fall.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

LAB = Path(__file__).resolve().parents[1]
if str(LAB / "src") not in sys.path:
    sys.path.insert(0, str(LAB / "src"))

from recursive_horizons import nsc_influence_chain as chain  # noqa: E402

CLOCK = LAB / "results" / "development" / "nsc-wall-clock-horizon-v1.json"
OUTPUT = LAB / "results" / "development" / "nsc-wall-clock-fall-v1.json"
SCHEMA = "NSC-WALL-CLOCK-FALL-v1"
PRE_REGISTERED_POWER = 0.25


def spread(values):
    value = [abs(v) for v in values]
    mean = sum(value) / len(value)
    return (sum((v - mean) ** 2 for v in value) / len(value)) ** 0.5 / mean


def main():
    if not CLOCK.is_file():
        raise SystemExit("missing horizon clock record: " + str(CLOCK))
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    clock = json.loads(CLOCK.read_text(encoding="utf-8"))
    rows = []
    for name, entry in clock["cases"].items():
        gap = entry["horizon_gap"]
        rise = entry["rise_rate"]
        fall = entry["fall_rate"]
        rows.append({
            "case": name,
            "mass_over_rm": entry["mass_over_rm"],
            "gap": gap,
            "rise": rise,
            "fall": fall,
            "fall_times_gap_quarter": fall * (gap ** PRE_REGISTERED_POWER),
            "fall_times_gap_half": fall * (gap ** 0.5),
            "rise_times_gap_half": rise * (gap ** 0.5),
        })
    quarter = [row["fall_times_gap_quarter"] for row in rows]
    half = [row["fall_times_gap_half"] for row in rows]
    half_rise = [row["rise_times_gap_half"] for row in rows]
    fitted_power = None
    largest = max(rows, key=lambda row: row["gap"])
    smallest = min(rows, key=lambda row: row["gap"])
    if largest["gap"] > 0.0 and smallest["gap"] > 0.0:
        import math
        fitted_power = math.log(abs(largest["fall"]) / abs(smallest["fall"])) \
            / math.log(largest["gap"] / smallest["gap"])
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "pre_registered_power": PRE_REGISTERED_POWER,
        "rows": rows,
        "fall_gap_half_spread": spread(half),
        "fall_gap_quarter_spread": spread(quarter),
        "rise_gap_half_spread": spread(half_rise),
        "fitted_fall_power_post_hoc": fitted_power,
        "reading": (
            "the quarter-power law was declared before the products were "
            "formed; the fitted power is reported separately as post-hoc"
        ),
        "notes": [
            "rise uses the round-26 half-power law for reference in the "
            "same file",
        ],
        "input_sha256": {
            "lab/results/development/nsc-wall-clock-horizon-v1.json":
                chain.path_sha256(CLOCK),
        },
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_wall_clock_fall.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_wall_clock_fall.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    for row in rows:
        print("  {:<16} fall*gap^1/4 {: .4f}   fall*gap^1/2 {: .4f}   rise*gap^1/2 {: .4f}".format(
            row["case"], row["fall_times_gap_quarter"], row["fall_times_gap_half"],
            row["rise_times_gap_half"],
        ))
    print("  spreads: fall*gap^1/2 {:.1%} -> fall*gap^1/4 {:.1%}   (rise*gap^1/2 {:.1%})".format(
        record["fall_gap_half_spread"], record["fall_gap_quarter_spread"],
        record["rise_gap_half_spread"],
    ))
    print("  post-hoc fitted fall power: {!s}".format(
        "none" if fitted_power is None else "{:.4f}".format(fitted_power)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())