"""Is the wall's own clock the near-horizon scale?

Round-25 refuted the infall clock (mu-blind) for the mu-tuned breathing.
The remaining candidate is built from the horizon gap Delta = r_+ - r_-.
For the audited P^2 = 1 the gap swings 4.5x across the family, the right
shape to absorb the 2.3x rise-rate swing. This driver tests the
parameter-free form rate ~ Delta^-1/2 by measuring the spread of
rate * sqrt(Delta) across the family, against the raw rate spread.

    python scripts/lab.py scripts/derive_nsc_wall_clock_horizon.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

LAB = Path(__file__).resolve().parents[1]
if str(LAB / "src") not in sys.path:
    sys.path.insert(0, str(LAB / "src"))

from recursive_horizons import nsc_influence_chain as chain  # noqa: E402
from recursive_horizons import nsc_rn_source as source  # noqa: E402

BREATHING = LAB / "results" / "development" / "nsc-wall-breathing-family-v1.json"
OUTPUT = LAB / "results" / "development" / "nsc-wall-clock-horizon-v1.json"
SCHEMA = "NSC-WALL-CLOCK-HORIZON-v1"
MASSES = {"mu-1.002-eps-0.001": 1.002, "mu-1.010-eps-0.001": 1.010,
          "mu-1.040-eps-0.001": 1.040}


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
        model = source.action_reference(mass)
        r_plus = float(model.horizons["r_plus"])
        r_minus = float(model.horizons["r_minus"])
        gap = r_plus - r_minus
        entry = breathing["cases"][name]
        rise = entry["mean_rise_rate"]
        fall = entry["mean_fall_rate"]
        cases[name] = {
            "mass_over_rm": mass,
            "r_plus": r_plus,
            "r_minus": r_minus,
            "horizon_gap": gap,
            "rise_rate": rise,
            "fall_rate": fall,
            "rise_times_sqrt_gap": None if rise is None else rise * (gap ** 0.5),
            "fall_times_sqrt_gap": None if fall is None else fall * (gap ** 0.5),
        }
    raw_rise = [entry["rise_rate"] for entry in cases.values()]
    raw_fall = [entry["fall_rate"] for entry in cases.values()]
    gapped_rise = [entry["rise_times_sqrt_gap"] for entry in cases.values()]
    gapped_fall = [entry["fall_times_sqrt_gap"] for entry in cases.values()]
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "cases": cases,
        "raw_rise_spread": spread(raw_rise),
        "raw_fall_spread": spread(raw_fall),
        "rise_sqrt_gap_spread": spread(gapped_rise),
        "fall_sqrt_gap_spread": spread(gapped_fall),
        "reading": (
            "if rate ~ Delta^-1/2 the gapped spread collapses below the raw "
            "spread; if not, the near-horizon clock is refuted too"
        ),
        "notes": [
            "Delta = r_+ - r_- of the vacuum chart at each mu; parameter-free",
        ],
        "input_sha256": {
            "lab/results/development/nsc-wall-breathing-family-v1.json":
                chain.path_sha256(BREATHING),
        },
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_wall_clock_horizon.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_wall_clock_horizon.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    for name, entry in cases.items():
        print(
            "  {:<20} gap {:.4f} rise {:.4f} rise*sqrt {:.4f} fall*sqrt {:.4f}".format(
                name, entry["horizon_gap"], entry["rise_rate"],
                entry["rise_times_sqrt_gap"], entry["fall_times_sqrt_gap"],
            )
        )
    print("  spreads: raw rise {:.1%} -> gapped {:.1%}   raw fall {:.1%} -> gapped {:.1%}".format(
        record["raw_rise_spread"], record["rise_sqrt_gap_spread"],
        record["raw_fall_spread"], record["fall_sqrt_gap_spread"],
    ))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())