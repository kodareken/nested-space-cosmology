"""The deep-grid wall radii: the breathing envelope's second domain.

Rounds 22-24 measured the wall's radial breathing on the family chart.
The march balance record stored the walls' first radii on the
horizon-reaching grid but the trajectory was never read. This driver
extracts the per-station r_first_wall for the wall-bearing march
stations, expresses them in r_+ units, and holds the ranges against
the recorded family envelopes - the first domain-parity test of the
wall radius law.

    python scripts/lab.py scripts/derive_nsc_march_wall_radii.py
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

MARCH_BALANCE = LAB / "results" / "development" / "nsc-march-balance-trajectory-v1.json"
BREATHING = LAB / "results" / "development" / "nsc-wall-breathing-family-v1.json"
OUTPUT = LAB / "results" / "development" / "nsc-march-wall-radii-v1.json"
SCHEMA = "NSC-MARCH-WALL-RADII-v1"
MASSES = {"mu-1.002-eps-0.001": 1.002, "mu-1.010-eps-0.001": 1.010,
          "mu-1.040-eps-0.001": 1.040}


def main():
    if not MARCH_BALANCE.is_file() or not BREATHING.is_file():
        raise SystemExit("missing march balance or breathing record")
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    march = json.loads(MARCH_BALANCE.read_text(encoding="utf-8"))
    breathing = json.loads(BREATHING.read_text(encoding="utf-8"))
    radii = {}
    for name, rows in march["series"].items():
        model = source.action_reference(MASSES[name])
        r_plus = float(model.horizons["r_plus"])
        points = []
        for row in rows:
            if row["r_first_wall"] is not None:
                points.append({
                    "t_over_rm": row["t_over_rm"],
                    "r_first_wall": row["r_first_wall"],
                    "r_first_in_r_plus_units": row["r_first_wall"] / r_plus,
                })
        values = [point["r_first_in_r_plus_units"] for point in points]
        radii[name] = {
            "r_plus": r_plus,
            "points": points,
            "range_r_plus_units": [min(values), max(values)] if values else None,
            "family_envelope": breathing["cases"][name]["radial_range"],
            "inside_family_envelope": (
                bool(values) and
                breathing["cases"][name]["radial_range"][0] <= min(values)
                and max(values) <= breathing["cases"][name]["radial_range"][1]),
        }
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "radii": radii,
        "reading": (
            "the deep-grid wall radii of the light cases, held against "
            "the family-chart breathing envelope - the wall radius law's "
            "second domain"
        ),
        "notes": [
            "the heavy case has no walls after t=10 on this domain; its "
            "radii appear only where walls exist",
        ],
        "input_sha256": {
            "lab/results/development/nsc-march-balance-trajectory-v1.json":
                chain.path_sha256(MARCH_BALANCE),
            "lab/results/development/nsc-wall-breathing-family-v1.json":
                chain.path_sha256(BREATHING),
        },
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_march_wall_radii.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_march_wall_radii.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    for name, entry in radii.items():
        print("  {:<20} r_+ {:.3f}  range {!s}  family envelope {} inside {}".format(
            name, entry["r_plus"],
            None if entry["range_r_plus_units"] is None
            else "[{:.3f},{:.3f}]".format(*entry["range_r_plus_units"]),
            entry["family_envelope"], entry["inside_family_envelope"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())