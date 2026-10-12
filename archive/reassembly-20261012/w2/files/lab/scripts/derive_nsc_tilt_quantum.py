"""The tilt staircase's quantum: one node at the seam?

Round-46 found the fringe fee advancing in exact quanta (0.00131) with
flat plateaus on the 0.005-eta lattice. Pre-registered reading: the
quantum is the L1 share of the single node nearest the seam (the only
candidate unit of this chart at fixed seam). This driver computes, per
eta, the single-node share of |d_lapse| at the 12.5%-cap seam and holds
it against the recorded quantum.

    python scripts/lab.py scripts/derive_nsc_tilt_quantum.py
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

CROSSINGS = LAB / "results" / "development" / "nsc-tilt-cell-crossings-v1.json"
OUTPUT = LAB / "results" / "development" / "nsc-tilt-quantum-v1.json"
SCHEMA = "NSC-TILT-QUANTUM-v1"
A = 0.04501936182826115
C_F = 0.01125484045706527
ETAS = (1.015, 1.020, 1.025, 1.030, 1.035, 1.040, 1.045)
SEAM_FRACTION = 0.125
RECORDED_QUANTUM = 0.00131


def _window(radius, left, right, sharpness=8.0):
    return 0.5 * (1.0 + np.tanh(sharpness * (radius - left))) \
        * 0.5 * (1.0 + np.tanh(sharpness * (right - radius)))


def main():
    if not CROSSINGS.is_file():
        raise SystemExit("missing crossings record: " + str(CROSSINGS))
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
        seam = SEAM_FRACTION * cap
        nearest = int(np.argmin(np.abs(radius - seam)))
        magnitude = np.abs(tangent["d_lapse"]) * grid["weights"]
        total = float(np.sum(magnitude))
        single_node_share = {
            "node_index": nearest,
            "seam": seam,
            "node_radius": float(radius[nearest]),
            "share": float(magnitude[nearest]) / total,
        }
        rows.append({"eta": eta, "single_node_share": single_node_share})
    shares = [row["single_node_share"]["share"] for row in rows]
    mean_share = sum(shares) / len(shares)
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "hypothesis_pre_registered": (
            "the staircase quantum equals the single-node L1 share at the seam"
        ),
        "recorded_quantum": RECORDED_QUANTUM,
        "rows": rows,
        "mean_single_node_share": mean_share,
        "ratio_share_to_quantum": mean_share / RECORDED_QUANTUM,
        "verdict": (
            "QUANTUM_IS_SINGLE_NODE"
            if 0.5 <= mean_share / RECORDED_QUANTUM <= 2.0 else "QUANTUM_UNEXPLAINED"
        ),
        "notes": [
            "the quantum 0.00131 is the round-46 recorded step size; the "
            "single-node share isolates the node nearest the seam",
        ],
        "input_sha256": {
            "lab/results/development/nsc-tilt-cell-crossings-v1.json":
                chain.path_sha256(CROSSINGS),
        },
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_tilt_quantum.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_tilt_quantum.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    for row in rows:
        print("  eta {:<5} share {:.5f} (node {} at {:.3f})".format(
            row["eta"], row["single_node_share"]["share"],
            row["single_node_share"]["node_index"], row["single_node_share"]["node_radius"]))
    print("  mean share {:.5f} vs quantum 0.00131  ratio {:.2f}  verdict: {}".format(
        record["mean_single_node_share"], record["ratio_share_to_quantum"],
        record["verdict"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())