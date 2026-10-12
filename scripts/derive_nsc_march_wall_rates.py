"""The deep-domain breathing rates, and the clock law's second domain.

Round-81 recorded the march-domain wall radii per station. This driver
differentiates those series into rise/fall rates (in r_+ units per t)
for every case, compares the two light cases against the family-domain
rates, and runs the same pre-registered gap products on the deep-domain
measurements - the first second-domain test of the wall's clock pair.

    python scripts/lab.py scripts/derive_nsc_march_wall_rates.py
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

RADII = LAB / "results" / "development" / "nsc-march-wall-radii-v1.json"
FAMILY_RATES = LAB / "results" / "development" / "nsc-wall-breathing-family-v1.json"
OUTPUT = LAB / "results" / "development" / "nsc-march-wall-rates-v1.json"
SCHEMA = "NSC-MARCH-WALL-RATES-v1"
MASSES = {"mu-1.002-eps-0.001": 1.002, "mu-1.010-eps-0.001": 1.010,
          "mu-1.040-eps-0.001": 1.040}


def main():
    if not RADII.is_file() or not FAMILY_RATES.is_file():
        raise SystemExit("missing radii or family rates record")
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    radii = json.loads(RADII.read_text(encoding="utf-8"))
    family = json.loads(FAMILY_RATES.read_text(encoding="utf-8"))
    rows = []
    for name, mass in MASSES.items():
        model = source.action_reference(mass)
        gap = float(model.horizons["r_plus"]) - float(model.horizons["r_minus"])
        points = radii["radii"][name]["points"]
        rising = []
        falling = []
        for left, right in zip(points, points[1:]):
            dt = right["t_over_rm"] - left["t_over_rm"]
            if dt <= 0.0:
                continue
            velocity = (right["r_first_in_r_plus_units"] - left["r_first_in_r_plus_units"]) / dt
            (rising if velocity > 0.0 else falling).append(velocity)
        mean = lambda items: (sum(items) / len(items)) if items else None
        row = {
            "case": name,
            "gap": gap,
            "deep_rise_rate": mean(rising),
            "deep_fall_rate": mean(falling),
            "family_rise_rate": family["cases"][name]["mean_rise_rate"],
            "family_fall_rate": family["cases"][name]["mean_fall_rate"],
        }
        row["deep_rise_times_sqrt_gap"] = (
            row["deep_rise_rate"] * gap ** 0.5 if row["deep_rise_rate"] else None)
        row["deep_fall_times_sqrt_quarter_gap"] = (
            abs(row["deep_fall_rate"]) * gap ** 0.25 if row["deep_fall_rate"] else None)
        rows.append(row)
    usable = [row for row in rows
              if row["deep_rise_times_sqrt_gap"] is not None]
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "rows": rows,
        "usable_rise_cases": [row["case"] for row in usable],
        "reading": (
            "the deep-domain breathing rates differentiated from the "
            "recorded wall-radius series; the gap products are the same "
            "pre-registered forms as rounds 26/30"
        ),
        "notes": [
            "heavy-case walls exist at two stations only; its rates carry "
            "one interval",
            "the sparse march cadence (2-3 intervals per case) is "
            "insufficient for a second-domain power-pair test; no power "
            "is claimed from these rates",
        ],
        "input_sha256": {
            "lab/results/development/nsc-march-wall-radii-v1.json": chain.path_sha256(RADII),
            "lab/results/development/nsc-wall-breathing-family-v1.json":
                chain.path_sha256(FAMILY_RATES),
        },
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_march_wall_rates.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_march_wall_rates.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    for row in rows:
        print("  {:<20} deep rise {!s}  deep fall {!s}  family rise {:.4f} fall {:.4f}".format(
            row["case"],
            None if row["deep_rise_rate"] is None else "{:+.4f}".format(row["deep_rise_rate"]),
            None if row["deep_fall_rate"] is None else "{:+.4f}".format(row["deep_fall_rate"]),
            row["family_rise_rate"], row["family_fall_rate"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())