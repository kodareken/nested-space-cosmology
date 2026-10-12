"""Placement of every wall kinematic constant inside the outward cell.

Round-54 joined the one-way low band [0, 0.5] to the wall's rise/fall
rates. This driver completes the set with the round-39 oscillator
constant omega = 0.386 and records, for each of the three wall
kinematic numbers, its placement inside the single outward cell of the
mu=1.040 current book and its margin to the band edge. No new motion;
pure table over recorded bytes.

    python scripts/lab.py scripts/derive_nsc_lowband_placement.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

LAB = Path(__file__).resolve().parents[1]
if str(LAB / "src") not in sys.path:
    sys.path.insert(0, str(LAB / "src"))

from recursive_horizons import nsc_influence_chain as chain  # noqa: E402

ONEWAY_WALL = LAB / "results" / "development" / "nsc-oneway-wall-join-v1.json"
DERIVATION = LAB / "results" / "development" / "nsc-wall-clock-derivation-v1.json"
OUTPUT = LAB / "results" / "development" / "nsc-lowband-placement-v1.json"
SCHEMA = "NSC-LOWBAND-PLACEMENT-v1"
BAND = (0.0, 0.5)


def main():
    if not ONEWAY_WALL.is_file() or not DERIVATION.is_file():
        raise SystemExit("missing join or derivation record")
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    join = json.loads(ONEWAY_WALL.read_text(encoding="utf-8"))
    derivation = json.loads(DERIVATION.read_text(encoding="utf-8"))
    omega = derivation["omega_mean"]
    rise = join["wall_rates"]["rise"]
    fall = join["wall_rates"]["fall"]
    transient = join["oneway_positive_cells_at_mu_1040"][0]["interior_max_direction"]
    entries = [
        {"name": "rise_rate", "value": rise},
        {"name": "fall_rate", "value": fall},
        {"name": "wall_oscillator_omega", "value": omega},
    ]
    rows = []
    for entry in entries:
        value = entry["value"]
        rows.append({
            "name": entry["name"],
            "value": value,
            "band": list(BAND),
            "placement_fraction_of_band": value / BAND[1],
            "margin_to_lower_edge": value - BAND[0],
            "margin_to_upper_edge": BAND[1] - value,
            "ratio_to_transient_magnitude": value / transient
            if transient else None,
        })
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "band": list(BAND),
        "transient_magnitude": transient,
        "rows": rows,
        "all_inside": bool(all(row["value"] < BAND[1] for row in rows)),
        "verdict": (
            "ALL_WALL_KINEMATICS_INSIDE_OUTWARD_CELL"
            if all(row["value"] < BAND[1] for row in rows) else "PARTIAL"
        ),
        "notes": [
            "the three constants come from three different records "
            "(breathing, clock derivation) and share 1/t units with the "
            "one-way frequency bands",
            "placement table only; no dynamical mechanism is claimed",
        ],
        "input_sha256": {
            "lab/results/development/nsc-oneway-wall-join-v1.json":
                chain.path_sha256(ONEWAY_WALL),
            "lab/results/development/nsc-wall-clock-derivation-v1.json":
                chain.path_sha256(DERIVATION),
        },
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_lowband_placement.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_lowband_placement.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    for row in rows:
        print("  {:<24} value {:.4f}  band placement {:.1%}  upper margin {:.4f}".format(
            row["name"], row["value"], row["placement_fraction_of_band"],
            row["margin_to_upper_edge"]))
    print("  transient magnitude +{:.4f}  verdict: {}".format(
        transient, record["verdict"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())