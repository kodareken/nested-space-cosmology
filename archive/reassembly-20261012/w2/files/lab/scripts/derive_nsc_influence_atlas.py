"""Cross-instrument interval atlas of the radial influence.

Two independent instruments measure the same partition weights on the
mu = 1.04 chart, intervals [grid_0, 8], [8, 16], [16, 24], [24, 32]:

    A. the stage tangent's lapse response to an outer-localized direction
       (L1 share per interval), the instrument of the influence funnel;
    B. the infall-winding time integral on the vacuum chart of the same
       mass (share of T_total per interval), the wavefront clock.

Shared edges, independent code paths. The plateau constancy of A below
the source edge is also measured (the linear-inward law from round 11).

    python scripts/lab.py scripts/derive_nsc_influence_atlas.py
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

EDGES = (8.0, 16.0, 24.0, 32.0)
OUTER_START_OVER_RM = 24.0
OUTPUT = LAB / "results" / "development" / "nsc-influence-atlas-v1.json"
SCHEMA = "NSC-INFLUENCE-ATLAS-v1"


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
    direction = {
        "dphi": np.asarray(packet["phi"]) * _window(
            radius, OUTER_START_OVER_RM, float(geometry["route"]["outer"]))[None, :, None],
        "d_occupations": np.zeros_like(state["occupations"]),
        "d_inner_mass": 0.0,
    }
    tangent = response.stage_tangent(state, grid, direction)
    d_lapse = tangent["d_lapse"]
    full_edges = [float(np.min(radius))] + list(EDGES)
    shares_a = chain.interval_shares(radius, weights, d_lapse, full_edges)
    below_source = radius < OUTER_START_OVER_RM
    plateau = float(np.max(np.abs(d_lapse[below_source]
                                  - np.mean(d_lapse[below_source])))) \
        / float(np.mean(np.abs(d_lapse[below_source])))
    full_mass_shares = chain.interval_shares(radius, weights, tangent["d_mass"], full_edges)

    # Instrument B: infall time on the vacuum chart of the same mass.
    model = source.action_reference(geometry["mass"])
    r_plus_h = float(model.horizons["r_plus"])
    fine = np.linspace(r_plus_h, 32.0, 2400)
    fval = 1.0 - 2.0 * float(model.mass) / fine + float(model.magnetic_r2) / fine ** 2
    weight_b = 1.0 / (1.0 + np.sqrt(np.maximum(0.0, fval)))
    edges_b = [r_plus_h] + list(EDGES)
    shares_b = []
    for left, right in zip(edges_b, edges_b[1:]):
        inside = (fine >= left) & (fine < right)
        shares_b.append({
            "interval": [left, right],
            "l1_share": float(np.trapezoid(weight_b[inside], fine[inside]))
            / float(np.trapezoid(weight_b, fine)),
        })

    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "mass_over_rm": 1.04,
        "intervals": [[float(np.min(radius))] + list(EDGES)],
        "instrument_a_tangent_lapse": shares_a,
        "instrument_a_tangent_mass": full_mass_shares,
        "instrument_b_infall_time": shares_b,
        "plateau_variation_below_source": plateau,
        "plateau_reading": (
            "round-11 linear profile requires d_lapse nearly flat below the "
            "source edge; this is the measured relative variation"
        ),
        "notes": [
            "instrument A and B share edges but not code paths: A is the "
            "stage tangent of the sourced chart, B the time quadrature of "
            "the vacuum chart at the same mass",
            "comparison is observational; no equality is asserted in advance",
        ],
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_influence_atlas.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_influence_atlas.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    print("instrument A (tangent lapse) shares:")
    for item in shares_a:
        print("  {!s:24} {:.6f}".format(str(item["interval"]), item["l1_share"]))
    print("instrument B (infall time) shares:")
    for item in shares_b:
        print("  {!s:24} {:.6f}".format(str(item["interval"]), item["l1_share"]))
    print("plateau variation below source: {:.3e}".format(plateau))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())