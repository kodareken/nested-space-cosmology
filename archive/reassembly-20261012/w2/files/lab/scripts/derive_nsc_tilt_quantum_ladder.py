"""The quantum ladder: does the fringe deliver integer quanta per seam?

Round-46 resolved the 12.5%-seam staircase into a single-node quantum
(0.00131). The round-43 record hints at more: the per-seam tilts over
the eta span are ~3x, ~2x, ~1x, ~0x of that quantum (3 at 12.5%, 2 at
25%, ~1 at 50%, ~0 at 75%). Pre-registered: one tangent per eta on the
fine lattice, four seam fees each; each seam's active steps must equal
the same quantum (within 15%), and its net tilt must be an integer
multiple of it. The ladder {3,2,1,0} would be the fringe's delivery
profile across the scaffold.

    python scripts/lab.py scripts/derive_nsc_tilt_quantum_ladder.py
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

OUTPUT = LAB / "results" / "development" / "nsc-tilt-quantum-ladder-v1.json"
SCHEMA = "NSC-TILT-QUANTUM-LADDER-v1"
A = 0.04501936182826115
C_F = 0.01125484045706527
ETAS = (1.015, 1.020, 1.025, 1.030, 1.035, 1.040, 1.045)
SEAM_FRACTIONS = (0.125, 0.25, 0.5, 0.75)
QUANTUM = 0.00131


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
    series = {fraction: [] for fraction in SEAM_FRACTIONS}
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
        for fraction in SEAM_FRACTIONS:
            series[fraction].append(chain.spatial_fee(
                grid["radius"], grid["weights"], tangent["d_lapse"], fraction * cap,
            )["l1_mass_response_fraction_inside"])
    seams = {}
    for fraction in SEAM_FRACTIONS:
        values = series[fraction]
        steps = [values[i + 1] - values[i] for i in range(len(values) - 1)]
        active = [step for step in steps if abs(step) > 0.0005]
        quantum_ok = (
            all(abs(abs(step) - QUANTUM) <= 0.15 * QUANTUM for step in active)
            if active else True
        )
        net = values[-1] - values[0]
        ratio = net / QUANTUM
        seams[str(fraction)] = {
            "series": values,
            "steps": steps,
            "active_steps": active,
            "all_active_steps_equal_quantum": quantum_ok,
            "net_tilt": net,
            "net_quantum_ratio": ratio,
            "nearest_integer": round(ratio),
        }
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "hypothesis_pre_registered": (
            "every seam's active steps carry the same quantum, and each "
            "seam's net tilt is an integer multiple of it"
        ),
        "quantum": QUANTUM,
        "seams": seams,
        "verdict": (
            "QUANTUM_LADDER_HOLDS"
            if all(entry["all_active_steps_equal_quantum"] for entry in seams.values())
            and all(abs(entry["net_quantum_ratio"] - entry["nearest_integer"]) < 0.25
                    for entry in seams.values())
            else "QUANTUM_LADDER_BROKEN"
        ),
        "notes": [
            "one tangent per eta; all four seam fees from the same d_lapse",
            "active steps are those above half the quantum",
        ],
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_tilt_quantum_ladder.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_tilt_quantum_ladder.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    for fraction, entry in seams.items():
        print("  seam {:<5} steps {!s} active {!s} net {:.5f} ratio {:.2f} -> {}".format(
            fraction, [round(v, 5) for v in entry["steps"]],
            [round(v, 5) for v in entry["active_steps"]], entry["net_tilt"],
            entry["net_quantum_ratio"], entry["nearest_integer"]))
    print("  verdict:", record["verdict"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())