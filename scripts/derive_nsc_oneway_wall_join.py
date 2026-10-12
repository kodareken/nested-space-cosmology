"""One-way v3 join: the +0.283 transient band versus the wall clock.

The one-way programme records, at mu=1.040, a single outward cell in
the signed-current book: the low band [0, 0.5] with +0.283; every other
band, and all bands at mu=1.010, point interior-ward. The wall
programme's kinetic rates (rise 0.131, fall 0.242, in 1/t units) live
inside exactly that band, and its mass book recorded zero positive
outer ticks in 20 continuation intervals. This driver lays the three
book statements on one table and pre-registers the containment reading:
the current book's only outward cell is the band that contains the
wall's own clock rates.

    python scripts/lab.py scripts/derive_nsc_oneway_wall_join.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

LAB = Path(__file__).resolve().parents[1]
if str(LAB / "src") not in sys.path:
    sys.path.insert(0, str(LAB / "src"))

from recursive_horizons import nsc_influence_chain as chain  # noqa: E402

ONEWAY = LAB / "results" / "development" / "nsc-rn-oneway-v3.json"
BREATHING = LAB / "results" / "development" / "nsc-wall-breathing-family-v1.json"
CLOSURE = LAB / "results" / "development" / "nsc-continuation-station-closure-v1.json"
OUTPUT = LAB / "results" / "development" / "nsc-oneway-wall-join-v1.json"
SCHEMA = "NSC-ONEWAY-WALL-JOIN-v1"


def main():
    if any(not p.is_file() for p in (ONEWAY, BREATHING, CLOSURE)):
        raise SystemExit("missing join input record")
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    oneway = json.loads(ONEWAY.read_text(encoding="utf-8"))
    breathing = json.loads(BREATHING.read_text(encoding="utf-8"))
    closure = json.loads(CLOSURE.read_text(encoding="utf-8"))
    case = next(c for c in oneway["cases"] if c["mass_over_rm"] == 1.04)
    bands = case["final_station"]["bands"]
    wall = breathing["cases"]["mu-1.040-eps-0.001"]
    rise = abs(wall["mean_rise_rate"])
    fall = abs(wall["mean_fall_rate"])
    positive_cells = [
        {"band": [band["band_low"], band["band_high"]],
         "interior_max_direction": band["interior_max_direction"],
         "exterior_min_direction": band["exterior_min_direction"]}
        for band in bands
        if band["interior_max_direction"] is not None
        and float(band["interior_max_direction"]) > 0.0
    ]
    low_band_positive = (
        0.0 < rise < 0.5 and 0.0 < fall < 0.5
        and any(cell["band"] == [0.0, 0.5] for cell in positive_cells)
    )
    outer_intervals = closure["intervals"]
    positive_mass_ticks = sum(delta > 0.0 for delta in
                              [item["delta_outer"] for item in outer_intervals])
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "oneway_positive_cells_at_mu_1040": positive_cells,
        "wall_rates": {"rise": rise, "fall": fall},
        "wall_rates_inside_low_band": bool(low_band_positive),
        "mass_book": {
            "outer_intervals": len(outer_intervals),
            "positive_outer_ticks": positive_mass_ticks,
            "verdict_from_closure": closure["verdict"],
        },
        "verdict": (
            "TRANSIENT_BAND_CONTAINS_WALL_CLOCK"
            if low_band_positive and positive_mass_ticks == 0 and positive_cells
            else "NO_JOIN"
        ),
        "notes": [
            "frequency bands and kinetic rates share 1/t units but are "
            "different instruments; the containment reading is hypothesis-"
            "labeled",
            "the mass book shows no positive tick in any of the recorded "
            "continuation intervals",
        ],
        "input_sha256": {
            "lab/results/development/nsc-rn-oneway-v3.json": chain.path_sha256(ONEWAY),
            "lab/results/development/nsc-wall-breathing-family-v1.json":
                chain.path_sha256(BREATHING),
            "lab/results/development/nsc-continuation-station-closure-v1.json":
                chain.path_sha256(CLOSURE),
        },
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_oneway_wall_join.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_oneway_wall_join.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    print("  one-way positive cells at mu=1.040:", positive_cells)
    print("  wall rates: rise {:.4f} fall {:.4f}  inside low band: {}".format(
        rise, fall, low_band_positive))
    print("  mass book: {} intervals, {} positive outer ticks".format(
        len(outer_intervals), positive_mass_ticks))
    print("  verdict:", record["verdict"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())