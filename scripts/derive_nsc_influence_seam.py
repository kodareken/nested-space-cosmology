"""Seam profile of the two-level influence funnel.

Extends the single-seam tangent probe to a seam sweep: for a direction
living in the outer level [24, 32] r_m, the lapse-response L1 fraction
inside seams 4, 8, 16 and 24 r_m is the inward influence profile. The
mirrored direction (inner-localized) reports the outward profile. This
profile is the fractional-level generalisation the chain closure needs,
and it can be compared with the resolvent ladder scale 1/Omega.

    python scripts/lab.py scripts/derive_nsc_influence_seam.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

LAB = Path(__file__).resolve().parents[1]
if str(LAB / "src") not in sys.path:
    sys.path.insert(0, str(LAB / "src"))
if str(LAB / "scripts") not in sys.path:
    sys.path.insert(0, str(LAB / "scripts"))

import numpy as np

import derive_nsc_rn_parent as driver
from recursive_horizons import nsc_influence_chain as chain
from recursive_horizons import nsc_rn_pg as pg
from recursive_horizons import nsc_rn_response as response
from recursive_horizons import nsc_rn_source as source

OUTER_START_OVER_RM = 24.0
SEAMS = (4.0, 8.0, 16.0, 24.0, 30.0)
OUTPUT = LAB / "results" / "development" / "nsc-influence-seam-v1.json"
SCHEMA = "NSC-INFLUENCE-SEAM-v1"


def _window(radius, left, right, sharpness=8.0):
    edge_left = 0.5 * (1.0 + np.tanh(sharpness * (radius - left)))
    edge_right = 0.5 * (1.0 + np.tanh(sharpness * (right - radius)))
    return edge_left * edge_right


def main():
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    geometry = driver.case_geometry(1.04)
    grid = pg.build_grid(
        geometry["mass"], r_m=geometry["r_m"], points=geometry["route"]["points"],
        r_out=geometry["route"]["outer"],
    )
    packet = source.prepare_radial_packet(
        mass_over_rm=1.04, energy_target=1.0e-3, rank=driver.PILOT_RANK,
        model=grid["model"], radii=grid["radius"], weights=grid["weights"],
        derivative=grid["derivative"],
    )
    state = pg.make_state(
        grid, np.asarray(packet["phi"]),
        inner_mass=geometry["mass"],
        killing_frequency=float(packet["radial_moments"]["requested_omega"]),
        occupations=np.asarray(packet["occupations"], dtype=float),
    )
    radius = grid["radius"]
    weights = grid["weights"]

    def tangent_for(left, right):
        mask = _window(radius, left, right)
        direction = {
            "dphi": np.asarray(packet["phi"]) * mask[None, :, None],
            "d_occupations": np.zeros_like(state["occupations"]),
            "d_inner_mass": 0.0,
        }
        return response.stage_tangent(state, grid, direction)

    forward = tangent_for(OUTER_START_OVER_RM, float(geometry["route"]["outer"]))
    backward = tangent_for(0.0, 8.0)
    rows = []
    for seam in SEAMS:
        f_fwd = chain.spatial_fee(radius, weights, forward["d_lapse"], seam)[
            "l1_mass_response_fraction_inside"]
        f_bwd = chain.spatial_fee(radius, weights, backward["d_lapse"], seam)[
            "l1_mass_response_fraction_inside"]
        rows.append({
            "seam_over_rm": float(seam),
            "inner_fraction_of_cap": float(seam / float(geometry["route"]["outer"])),
            "forward_inward_fraction": f_fwd,
            "backward_inward_fraction": f_bwd,
        })
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "outer_level": [OUTER_START_OVER_RM, float(geometry["route"]["outer"])],
        "seams": rows,
        "resolvent_ladder_scale": 1.0 / 3.973074368754331,
        "resolvent_ladder_scale_source": "round-1 scale binding Omega",
        "notes": [
            "the forward profile is the inward influence of one outer level; "
            "the backward profile is the mirrored locality margin",
            "fractions are L1 lapse-response fractions, same definition as "
            "the single-seam probe record",
        ],
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_influence_seam.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_influence_seam.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    for row in rows:
        print(
            "  seam={: .1f} (cap fraction {: .3f})  forward={: .6f}  backward={: .6f}".format(
                row["seam_over_rm"], row["inner_fraction_of_cap"],
                row["forward_inward_fraction"], row["backward_inward_fraction"],
            )
        )
    print("  resolvent ladder scale for comparison:", record["resolvent_ladder_scale"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())