"""Does the tilt staircase step on grid-cell crossings of the gap?

Round-45 recorded a fee staircase with a plateau at eta 1.03..1.04.
Pre-registered trigger: the fringe steps when the horizon zone width
gap crosses an integer multiple of the grid cell h (the only discrete
unit of this chart). This driver samples eta finely across the
suspected crossing band and labels each pair as crossing/non-crossing
by whether an integer multiple of h(eta) lies inside [gap_prev, gap_next].
The comparison of mean |delta fee| between the classes is the test.

    python scripts/lab.py scripts/derive_nsc_tilt_cell_crossings.py
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

OUTPUT = LAB / "results" / "development" / "nsc-tilt-cell-crossings-v1.json"
SCHEMA = "NSC-TILT-CELL-CROSSINGS-v1"
A = 0.04501936182826115
C_F = 0.01125484045706527
ETAS = (1.015, 1.020, 1.025, 1.030, 1.035, 1.040, 1.045)
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
        h = spacing
        rows.append({"eta": eta, "gap": gap, "h": h, "fee": fee})
    pairs = []
    for previous, current in zip(rows, rows[1:]):
        low = min(previous["gap"], current["gap"])
        high = max(previous["gap"], current["gap"])
        h_mean = 0.5 * (previous["h"] + current["h"])
        crossing = any(
            low <= k * h_mean <= high for k in range(1, 40)
        )
        pairs.append({
            "eta_pair": [previous["eta"], current["eta"]],
            "delta_fee": abs(current["fee"] - previous["fee"]),
            "gap_crosses_integer_cells": crossing,
        })
    crossing_pairs = [p for p in pairs if p["gap_crosses_integer_cells"]]
    quiet_pairs = [p for p in pairs if not p["gap_crosses_integer_cells"]]
    mean = lambda items: (sum(item["delta_fee"] for item in items) / len(items)) if items else None
    crossing_mean = mean(crossing_pairs)
    quiet_mean = mean(quiet_pairs)
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "hypothesis_pre_registered": (
            "the fringe staircase steps when the gap crosses an integer "
            "multiple of the grid cell"
        ),
        "rows": rows,
        "pairs": pairs,
        "crossing_mean_delta_fee": crossing_mean,
        "quiet_mean_delta_fee": quiet_mean,
        "verdict": (
            "QUANTIZED" if crossing_mean is not None and quiet_mean is not None
            and crossing_mean >= 2.0 * quiet_mean else "NOT_QUANTIZED"
        ),
        "notes": [
            "integer range 1..39 covers the gap/h ratio band",
            "crossing label built from the mean cell size of the pair",
        ],
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_tilt_cell_crossings.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_tilt_cell_crossings.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    for pair in pairs:
        print("  {:<6}->{:<6} dF {:.5f} crossing {}".format(
            pair["eta_pair"][0], pair["eta_pair"][1], pair["delta_fee"],
            pair["gap_crosses_integer_cells"]))
    print("  crossing mean {!s}  quiet mean {!s}  verdict: {}".format(
        crossing_mean, quiet_mean, record["verdict"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())