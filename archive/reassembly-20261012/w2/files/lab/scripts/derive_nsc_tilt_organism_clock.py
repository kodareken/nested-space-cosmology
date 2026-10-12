"""Does the fringe tilt step with the organism's sqrt(gap) clock?

Round-44 established the charged inner-seam eta-tilt as fringe physics.
With only three eta values its scaling was not decidable. The fixed-
probe charged tangent is family-free, so this driver extends the series
to eta = 1.01..1.05 (five points, seam 12.5% of cap, 876 nodes) and
pre-registers the organism-clock hypothesis: the tilt increment across
adjacent eta steps proportional to the step of sqrt(gap). The ratio
series delta_i / d sqrt(gap)_i is measured; constancy is the test.

    python scripts/lab.py scripts/derive_nsc_tilt_organism_clock.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

LAB = Path(__file__).resolve().parents[1]
if str(LAB / "src") not in sys.path:
    sys.path.insert(0, str(LAB / "src"))

import numpy as np

from recursive_horizons import nsc_influence_chain as chain
from recursive_horizons import nsc_rn_charged_lll as lll
from recursive_horizons import nsc_rn_pg as pg
from recursive_horizons import nsc_rn_response as response
from recursive_horizons import nsc_rn_source as source

OUTPUT = LAB / "results" / "development" / "nsc-tilt-organism-clock-v1.json"
SCHEMA = "NSC-TILT-ORGANISM-CLOCK-v1"
A = 0.04501936182826115
C_F = 0.01125484045706527
ETAS = (1.01, 1.02, 1.03, 1.04, 1.05)
SEAM_FRACTION = 0.125


def _window(radius, left, right, sharpness=8.0):
    return 0.5 * (1.0 + np.tanh(sharpness * (radius - left))) \
        * 0.5 * (1.0 + np.tanh(sharpness * (right - radius)))


def main():
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    template_ref = lll.dyonic_reference(
        A=A, C_F=C_F, q_mag=1.0, Q_E=4.0, mass_over_total_radius=1.01,
    )
    base_grid = lll.charged_grid(
        model=template_ref["model"], r_m_tot=template_ref["r_m_tot"], points=876)
    base_spacing = float(base_grid["spacing"])
    base_radius = np.asarray(base_grid["inner_edge"]) \
        + base_spacing * np.arange(876, dtype=float)
    weights, derivative, _ = pg.sbp_operator(876, base_spacing, layout="csr")
    packet = source.prepare_radial_packet(
        mass_over_rm=1.01, energy_target=1.0e-3, rank=9,
        model=template_ref["model"], radii=base_radius, weights=weights,
        derivative=derivative,
    )
    rows = []
    for eta in ETAS:
        ref = lll.dyonic_reference(
            A=A, C_F=C_F, q_mag=1.0, Q_E=4.0, mass_over_total_radius=eta,
        )
        grid = lll.charged_grid(model=ref["model"], r_m_tot=ref["r_m_tot"], points=876)
        spacing = float(grid["spacing"])
        radius = np.asarray(grid["inner_edge"]) + spacing * np.arange(876, dtype=float)
        weights_g, derivative_g, q_matrix = pg.sbp_operator(876, spacing, layout="csr")
        grid = dict(grid)
        grid.update({
            "radius": radius, "weights": weights_g, "derivative": derivative_g,
            "q_matrix": q_matrix, "outer": float(radius[-1]), "charge": 0.0,
        })
        state = pg.make_state(
            grid, np.asarray(packet["phi"]), inner_mass=float(ref["model"].mass),
            killing_frequency=float(packet["radial_moments"]["requested_omega"]),
            occupations=np.asarray(packet["occupations"], dtype=float),
        )
        cap = 32.0 * ref["r_m_tot"]
        direction = {
            "dphi": np.asarray(packet["phi"]) * _window(
                grid["radius"], 0.75 * cap, cap)[None, :, None],
            "d_occupations": np.zeros_like(state["occupations"]),
            "d_inner_mass": 0.0,
        }
        tangent = response.stage_tangent(state, grid, direction)
        fee = chain.spatial_fee(
            grid["radius"], grid["weights"], tangent["d_lapse"],
            SEAM_FRACTION * cap,
        )["l1_mass_response_fraction_inside"]
        gap = float(ref["horizons"]["r_plus"]) - float(ref["horizons"]["r_minus"])
        rows.append({
            "eta": eta,
            "gap": gap,
            "sqrt_gap": gap ** 0.5,
            "fee": fee,
        })
    steps = []
    for previous, current in zip(rows, rows[1:]):
        delta_fee = current["fee"] - previous["fee"]
        delta_sqrt = current["sqrt_gap"] - previous["sqrt_gap"]
        steps.append({
            "eta_from": previous["eta"],
            "eta_to": current["eta"],
            "delta_fee": delta_fee,
            "delta_sqrt_gap": delta_sqrt,
            "ratio": delta_fee / delta_sqrt if delta_sqrt else None,
        })
    ratios = [step["ratio"] for step in steps]
    spread = (sum((v - sum(ratios) / len(ratios)) ** 2 for v in ratios) / len(ratios)) ** 0.5 \
        / abs(sum(ratios) / len(ratios))
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "hypothesis_pre_registered": (
            "tilt increments across eta proportional to the sqrt(gap) step "
            "if the fringe marches on the organism clock"
        ),
        "rows": rows,
        "steps": steps,
        "ratio_spread": spread,
        "verdict": (
            "ORGANISM_CLOCK_HOLDS" if spread < 0.20 else "ORGANISM_CLOCK_REJECTED"
        ),
        "notes": [
            "five eta points give four steps; the ratio constancy is the test",
            "the fixed 1.01 probe is unchanged from round 43/44",
        ],
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_tilt_organism_clock.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_tilt_organism_clock.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    for row in rows:
        print("  eta {:<5} gap {:.4f} sqrt {:.4f} fee {:.5f}".format(
            row["eta"], row["gap"], row["sqrt_gap"], row["fee"]))
    for step in steps:
        print("  step {:<5}->{:<5} dF {:.5f} dSqrt {:.4f} ratio {:.4f}".format(
            step["eta_from"], step["eta_to"], step["delta_fee"],
            step["delta_sqrt_gap"], step["ratio"]))
    print("  ratio spread: {:.1%}  verdict: {}".format(record["ratio_spread"], record["verdict"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())