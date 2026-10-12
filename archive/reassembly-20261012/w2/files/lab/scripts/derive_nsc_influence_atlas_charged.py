"""Cross-instrument atlas on the charged admitted slice.

The round-12 atlas is rebuilt on the admitted dyonic slice (audited A,
C_F; q_mag=1, Q_E=4, mass ratio 1.01). Instrument A is the stage
tangent's lapse response to an outer-localized direction (neutral kappa=1
state on the dyonic chart, so the tangent's own ownership holds);
instrument B is the infall-time quadrature of the same chart. Both are
partitioned at cap fractions 25/50/75% instead of absolute radii, which
is the scale-free version of the round-12 comparison.

    python scripts/lab.py scripts/derive_nsc_influence_atlas_charged.py
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

from recursive_horizons import nsc_influence_chain as chain
from recursive_horizons import nsc_rn_charged_lll as lll
from recursive_horizons import nsc_rn_pg as pg
from recursive_horizons import nsc_rn_response as response
from recursive_horizons import nsc_rn_source as source

A = 0.04501936182826115
C_F = 0.01125484045706527
Q_MAG = 1.0
Q_E = 4.0
MASS_RATIO = 1.01
OUTPUT = LAB / "results" / "development" / "nsc-influence-atlas-charged-v1.json"
SCHEMA = "NSC-INFLUENCE-ATLAS-CHARGED-v1"


def _window(radius, left, right, sharpness=8.0):
    edge_left = 0.5 * (1.0 + np.tanh(sharpness * (radius - left)))
    edge_right = 0.5 * (1.0 + np.tanh(sharpness * (right - radius)))
    return edge_left * edge_right


def main():
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    ref = lll.dyonic_reference(
        A=A, C_F=C_F, q_mag=Q_MAG, Q_E=Q_E, mass_over_total_radius=MASS_RATIO,
    )
    grid = lll.charged_grid(model=ref["model"], r_m_tot=ref["r_m_tot"], points=876)
    # The charged grid forces its last node onto the cap, compressing the
    # final cell. The tangent's jet contract wants strict uniformity, so
    # rebuild a strictly uniform SBP grid over the same span (the last node
    # then sits 0.05 h beyond the declared cap).
    cap = float(grid["outer"])
    spacing = float(grid["spacing"])
    radius = np.asarray(grid["inner_edge"]) + spacing * np.arange(876, dtype=float)
    weights, derivative, q_matrix = pg.sbp_operator(876, spacing, layout="csr")
    grid = dict(grid)
    grid.update({
        "radius": radius,
        "weights": weights,
        "derivative": derivative,
        "q_matrix": q_matrix,
        "outer": float(radius[-1]),
    })
    grid["charge"] = 0.0  # the neutral kappa=1 state carries no electric current
    packet = source.prepare_radial_packet(
        mass_over_rm=MASS_RATIO, energy_target=1.0e-3, rank=9,
        model=ref["model"], radii=grid["radius"], weights=grid["weights"],
        derivative=grid["derivative"],
    )
    state = pg.make_state(
        grid, np.asarray(packet["phi"]),
        inner_mass=float(ref["model"].mass),
        killing_frequency=float(packet["radial_moments"]["requested_omega"]),
        occupations=np.asarray(packet["occupations"], dtype=float),
    )
    radius = grid["radius"]
    weights = grid["weights"]
    cap = 32.0 * ref["r_m_tot"]
    source_edge = 0.75 * cap
    direction = {
        "dphi": np.asarray(packet["phi"]) * _window(
            radius, source_edge, cap)[None, :, None],
        "d_occupations": np.zeros_like(state["occupations"]),
        "d_inner_mass": 0.0,
    }
    tangent = response.stage_tangent(state, grid, direction)
    edges = [float(np.min(radius)), 0.25 * cap, 0.5 * cap, 0.75 * cap, cap]
    shares_a = chain.interval_shares(radius, weights, tangent["d_lapse"], edges)
    mass_shares = chain.interval_shares(radius, weights, tangent["d_mass"], edges)
    below_source = radius < source_edge
    plateau = float(np.max(np.abs(tangent["d_lapse"][below_source]
                                  - np.mean(tangent["d_lapse"][below_source])))) \
        / float(np.mean(np.abs(tangent["d_lapse"][below_source])))

    r_plus_h = float(ref["horizons"]["r_plus"])
    fine = np.linspace(r_plus_h, cap, 2400)
    fval = 1.0 - 2.0 * float(ref["model"].mass) / fine + float(ref["model"].magnetic_r2) / fine ** 2
    weight_b = 1.0 / (1.0 + np.sqrt(np.maximum(0.0, fval)))
    edges_b = [r_plus_h, 0.25 * cap, 0.5 * cap, 0.75 * cap, cap]
    total_b = float(np.trapezoid(weight_b, fine))
    shares_b = []
    for left, right in zip(edges_b, edges_b[1:]):
        inside = (fine >= left) & (fine < right)
        shares_b.append({
            "interval": [left, right],
            "l1_share": float(np.trapezoid(weight_b[inside], fine[inside])) / total_b,
        })

    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "slice": {
            "A": A, "C_F": C_F, "q_mag": Q_MAG, "Q_E": Q_E,
            "mass_over_total_radius": MASS_RATIO,
            "r_m_tot": ref["r_m_tot"], "r_plus": ref["horizons"]["r_plus"],
            "cap": cap,
        },
        "instrument_a_tangent_lapse": shares_a,
        "instrument_a_tangent_mass": mass_shares,
        "instrument_b_infall_time": shares_b,
        "plateau_variation_below_source": plateau,
        "neutral_round12_reference": {
            "lapse_shares": [0.2610, 0.2991, 0.2990, 0.1409],
            "infall_shares": [0.2415, 0.2564, 0.2514, 0.2491],
        },
        "notes": [
            "the tangent runs the neutral kappa=1 state on the dyonic chart; "
            "its ownership of the chart response is unchanged",
            "cap-fraction edges make the comparison scale-free across charts",
        ],
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_influence_atlas_charged.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_influence_atlas_charged.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    print("instrument A (tangent lapse) cap-fraction shares:")
    for item in shares_a:
        print("  {:<40} {:.6f}".format(
            "[" + ",".join("{:g}".format(round(v, 2)) for v in item["interval"]) + "]",
            item["l1_share"]))
    print("instrument B (infall time) cap-fraction shares:")
    for item in shares_b:
        print("  {:<40} {:.6f}".format(
            "[" + ",".join("{:g}".format(round(v, 2)) for v in item["interval"]) + "]",
            item["l1_share"]))
    print("mass shares:", [round(item["l1_share"], 6) for item in mass_shares])
    print("plateau variation: {:.3e}".format(plateau))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())