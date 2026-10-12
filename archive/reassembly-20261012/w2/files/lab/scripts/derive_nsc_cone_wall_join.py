"""Geometric join of the cone trap into the wall and window envelopes.

The cone-map v2 trap sits a fraction of a cell above the reference
horizon; the wall envelope breathes across [0.90, 1.92] r_+ and the
charged window shell spans [1.0, 2.0] r_+. This driver lays the three
recorded locations on one r_+,unit table per case and records the
nesting (cone trap inside both envelopes, the window shell bounding the
wall), plus the redshift-at-excision uniformity across cases.

    python scripts/lab.py scripts/derive_nsc_cone_wall_join.py
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

CONE = LAB / "results" / "development" / "nsc-rn-cone-map-v2.json"
BREATHING = LAB / "results" / "development" / "nsc-wall-breathing-family-v1.json"
OUTPUT = LAB / "results" / "development" / "nsc-cone-wall-join-v1.json"
SCHEMA = "NSC-CONE-WALL-JOIN-v1"
MASSES = {"mu-1.002-eps-0.001": 1.002, "mu-1.010-eps-0.001": 1.010,
          "mu-1.040-eps-0.001": 1.040}


def main():
    if not CONE.is_file() or not BREATHING.is_file():
        raise SystemExit("missing cone or breathing record")
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    cone = json.loads(CONE.read_text(encoding="utf-8"))
    breathing = json.loads(BREATHING.read_text(encoding="utf-8"))
    rows = []
    for case in cone["cases"]:
        name = case["case"]
        mass = MASSES[name]
        model = source.action_reference(mass)
        r_plus = float(model.horizons["r_plus"])
        trap = float(case["trap_radius_last"])
        wall = breathing["cases"][name]["radial_range"]
        rows.append({
            "case": name,
            "r_plus": r_plus,
            "trap_last_r_m": trap,
            "trap_in_r_plus_units": trap / r_plus,
            "wall_envelope_r_plus_units": wall,
            "wall_contains_trap": bool(wall[0] <= trap / r_plus <= wall[1]),
            "window_shell": [1.0, 2.0],
            "window_contains_trap": bool(1.0 <= trap / r_plus <= 2.0),
            "window_bounds_wall": bool(wall[0] >= 0.9 and wall[1] <= 2.0),
            "redshift_at_excision_last": float(
                cone["summary"]["redshift_at_excision_last_station"][name]),
        })
    redshifts = [row["redshift_at_excision_last"] for row in rows]
    mean = sum(redshifts) / len(redshifts)
    spread = (sum((v - mean) ** 2 for v in redshifts) / len(redshifts)) ** 0.5 / abs(mean)
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "rows": rows,
        "redshift_at_excision_spread": spread,
        "reading": (
            "three recorded locations nested on one r_+ unit table: the "
            "cone trap inside the breathing wall, the window shell "
            "bounding the wall, all per case"
        ),
        "notes": [
            "the property that the trap sits above r_+ is the recorded "
            "offset trap_minus_reference of the v2 family",
            "redshift uniformity is the sixth mu-blind scaffold quantity",
        ],
        "input_sha256": {
            "lab/results/development/nsc-rn-cone-map-v2.json": chain.path_sha256(CONE),
            "lab/results/development/nsc-wall-breathing-family-v1.json":
                chain.path_sha256(BREATHING),
        },
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_cone_wall_join.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_cone_wall_join.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    for row in rows:
        print("  {:<20} trap {:.4f} r_+  wall [{:.3f},{:.3f}] contains {}  window contains {}  bounds {}".format(
            row["case"], row["trap_in_r_plus_units"], row["wall_envelope_r_plus_units"][0],
            row["wall_envelope_r_plus_units"][1], row["wall_contains_trap"],
            row["window_contains_trap"], row["window_bounds_wall"]))
    print("  redshift-at-excision spread: {:.2e} (mu-uniform)".format(record["redshift_at_excision_spread"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())