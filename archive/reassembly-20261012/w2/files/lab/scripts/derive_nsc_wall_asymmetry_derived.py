"""Derivation test: the wall asymmetry law composed from the two clocks.

Round-30 left the power pair rise ~ gap^-1/2, fall ~ gap^-1/4. The
implied asymmetry law is A(mu) = (C_f / C_r) * gap^1/4 with constants
C_r = rise*gap^1/2, C_f = |fall|*gap^1/4. This driver performs the
leave-one-out test: for each case, C_r and C_f are built from the OTHER
two cases only, then the asymmetry is predicted for the held-out case
and compared with the recorded round-24 empirical asymmetry. Nothing is
fitted to the held-out case.

    python scripts/lab.py scripts/derive_nsc_wall_asymmetry_derived.py
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
BREATHING = LAB / "results" / "development" / "nsc-wall-breathing-family-v1.json"
FALL = LAB / "results" / "development" / "nsc-wall-clock-fall-v1.json"
OUTPUT = LAB / "results" / "development" / "nsc-wall-asymmetry-derived-v1.json"
SCHEMA = "NSC-WALL-ASYMMETRY-DERIVED-v1"


def main():
    if not CLOCK.is_file() or not BREATHING.is_file() or not FALL.is_file():
        raise SystemExit("missing clock, breathing, or fall record")
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    clock = json.loads(CLOCK.read_text(encoding="utf-8"))
    breathing = json.loads(BREATHING.read_text(encoding="utf-8"))
    fall_record = json.loads(FALL.read_text(encoding="utf-8"))
    fall_quarter = {
        row["case"]: abs(row["fall_times_gap_quarter"]) for row in fall_record["rows"]
    }
    sites = {}
    for name, entry in clock["cases"].items():
        gap = entry["horizon_gap"]
        sites[name] = {
            "gap": gap,
            "c_r": entry["rise_times_sqrt_gap"],
            "c_f": fall_quarter[name],
            "measured_asymmetry": breathing["cases"][name]["asymmetry_fall_over_rise"],
        }
    predictions = {}
    names = list(sites)
    for index, held_out in enumerate(names):
        others = [n for j, n in enumerate(names) if j != index]
        c_r = 0.5 * (sites[others[0]]["c_r"] + sites[others[1]]["c_r"])
        c_f = 0.5 * (sites[others[0]]["c_f"] + sites[others[1]]["c_f"])
        predicted = (c_f / c_r) * (sites[held_out]["gap"] ** 0.25)
        measured = sites[held_out]["measured_asymmetry"]
        predictions[held_out] = {
            "held_out_case": held_out,
            "constants_from": others,
            "predicted_asymmetry": predicted,
            "measured_asymmetry": measured,
            "relative_error": (predicted - measured) / measured,
        }
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "law": "A = (C_f/C_r) * gap^1/4, constants from the other two cases",
        "predictions": predictions,
        "notes": [
            "leave-one-out: the held-out case contributes nothing to its "
            "own prediction",
            "C_f/C_r values carry the round-26/30 product spreads into the "
            "composition",
        ],
        "input_sha256": {
            "lab/results/development/nsc-wall-clock-horizon-v1.json":
                chain.path_sha256(CLOCK),
            "lab/results/development/nsc-wall-breathing-family-v1.json":
                chain.path_sha256(BREATHING),
            "lab/results/development/nsc-wall-clock-fall-v1.json":
                chain.path_sha256(FALL),
        },
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_wall_asymmetry_derived.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_wall_asymmetry_derived.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    for name, entry in predictions.items():
        print("  predict {} from {}: {:.4f} vs measured {:.4f}  ({:+.1%})".format(
            name, entry["constants_from"], entry["predicted_asymmetry"],
            entry["measured_asymmetry"], entry["relative_error"],
        ))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())