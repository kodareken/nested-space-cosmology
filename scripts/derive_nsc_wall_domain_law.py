"""Domain law of the wall radius: amplitude and inner-edge offset.

Round-81 recorded the deep-grid wall radii; round-24 recorded the
family-chart breathing amplitudes. This driver lays the two domains
side by side per case and quantifies the two visible domain laws:
(a) the breathing amplitude ratio shallow/deep, and (b) the wall's
offset above the grid's inner edge. The offsets are read from the
recorded radii and the recorded inner edges of the Q-decomposition.

    python scripts/lab.py scripts/derive_nsc_wall_domain_law.py
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
BREATHING = LAB / "results" / "development" / "nsc-wall-breathing-family-v1.json"
DECOMPOSITION = LAB / "results" / "development" / "nsc-q-decomposition-v1.json"
OUTPUT = LAB / "results" / "development" / "nsc-wall-domain-law-v1.json"
SCHEMA = "NSC-WALL-DOMAIN-LAW-v1"
MASSES = {"mu-1.002-eps-0.001": 1.002, "mu-1.010-eps-0.001": 1.010,
          "mu-1.040-eps-0.001": 1.040}


def main():
    if any(not p.is_file() for p in (RADII, BREATHING, DECOMPOSITION)):
        raise SystemExit("missing domain-law input record")
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    radii = json.loads(RADII.read_text(encoding="utf-8"))
    breathing = json.loads(BREATHING.read_text(encoding="utf-8"))
    decomposition = json.loads(DECOMPOSITION.read_text(encoding="utf-8"))
    rows = []
    for name, mass in MASSES.items():
        model = source.action_reference(mass)
        r_plus = float(model.horizons["r_plus"])
        deep = radii["radii"][name]
        deep_range = deep["range_r_plus_units"]
        deep_amplitude = deep_range[1] - deep_range[0]
        family = breathing["cases"][name]["radial_range"]
        family_amplitude = family[1] - family[0]
        inner_edge = decomposition["rows"][0]["r_factor"]  # placeholder, replaced
        inner_edge = None
        for entry in decomposition["rows"]:
            if entry["case"] == name:
                inner_edge = 32.0 / entry["r_factor"]
        min_radius_rm = deep_range[0] * r_plus
        offset = None if inner_edge is None else min_radius_rm - inner_edge
        rows.append({
            "case": name,
            "r_plus": r_plus,
            "deep_amplitude_r_plus": deep_amplitude,
            "family_amplitude_r_plus": family_amplitude,
            "amplitude_ratio_shallow_over_deep": (
                family_amplitude / deep_amplitude if deep_amplitude else None),
            "inner_edge_rm": inner_edge,
            "wall_min_rm": min_radius_rm,
            "wall_min_offset_above_inner_edge": offset,
        })
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "rows": rows,
        "reading": (
            "the wall breathes wide on the shallow grid and hugs its "
            "deepest cell on the horizon-reaching grid; both effects are "
            "mu-ordered"
        ),
        "notes": [
            "inner edges are the march grid's edges from the Q-decomposition "
            "record; offsets are in r_m",
        ],
        "input_sha256": {
            "lab/results/development/nsc-march-wall-radii-v1.json": chain.path_sha256(RADII),
            "lab/results/development/nsc-wall-breathing-family-v1.json":
                chain.path_sha256(BREATHING),
            "lab/results/development/nsc-q-decomposition-v1.json":
                chain.path_sha256(DECOMPOSITION),
        },
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_wall_domain_law.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_wall_domain_law.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    for row in rows:
        print("  {:<20} deep amp {:.3f}  fam amp {:.3f}  ratio {:.1f}x  offset {:.4f} r_m".format(
            row["case"], row["deep_amplitude_r_plus"], row["family_amplitude_r_plus"],
            row["amplitude_ratio_shallow_over_deep"], row["wall_min_offset_above_inner_edge"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())