"""Resolution discriminator for the charged eta-tilt.

Round-43 recorded a weak systematic eta-tilt at the inner seams of the
dyonic scaffold. Two honest readings remain: a grid-inner-cell artifact
or a horizon-fringe effect entering the outside. They separate by
resolution: the cell-driven tilt halves when the grid density doubles;
a gravitational fringe persists. This driver re-runs the 12.5% and 25%
cap-seam fees at 876 and 1752 points for eta = 1.01, 1.02, 1.04 with
the same interpolated probe.

    python scripts/lab.py scripts/derive_nsc_charged_tilt_resolution.py
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

OUTPUT = LAB / "results" / "development" / "nsc-charged-tilt-resolution-v1.json"
SCHEMA = "NSC-CHARGED-TILT-RESOLUTION-v1"
A = 0.04501936182826115
C_F = 0.01125484045706527
ETAS = (1.01, 1.02, 1.04)
SEAM_FRACTIONS = (0.125, 0.25)
POINT_SETS = (876, 1752)


def _window(radius, left, right, sharpness=8.0):
    return 0.5 * (1.0 + np.tanh(sharpness * (radius - left))) \
        * 0.5 * (1.0 + np.tanh(sharpness * (right - radius)))


def template_probe(points):
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
    occupations = packet["occupations"]
    omega = packet["radial_moments"]["requested_omega"]
    if points == 876:
        return np.asarray(packet["phi"]), occupations, omega
    fraction = np.minimum(
        1.0, (base_radius - base_radius[0]) / (base_radius[-1] - base_radius[0]))
    base_phi = np.asarray(packet["phi"])
    was_3d = base_phi.ndim == 3
    if was_3d:
        base_phi = base_phi[..., 0]
    if base_phi.ndim == 1:
        rows = base_phi[None, :]
    else:
        rows = base_phi
    phi = np.stack([
        np.interp(np.linspace(0.0, 1.0, points), fraction, rows[column])
        for column in range(rows.shape[0])
    ])
    if phi.shape[0] == 1 and not was_3d:
        phi = phi[0]
    if was_3d:
        phi = phi[..., None]
    return phi, occupations, omega


def fee(eta, points, phi, occupations, omega):
    ref = lll.dyonic_reference(
        A=A, C_F=C_F, q_mag=1.0, Q_E=4.0, mass_over_total_radius=eta,
    )
    grid = lll.charged_grid(model=ref["model"], r_m_tot=ref["r_m_tot"], points=points)
    spacing = float(grid["spacing"])
    radius = np.asarray(grid["inner_edge"]) + spacing * np.arange(points, dtype=float)
    weights, derivative, q_matrix = pg.sbp_operator(points, spacing, layout="csr")
    grid = dict(grid)
    grid.update({
        "radius": radius, "weights": weights, "derivative": derivative,
        "q_matrix": q_matrix, "outer": float(radius[-1]), "charge": 0.0,
    })
    state = pg.make_state(
        grid, np.asarray(phi), inner_mass=float(ref["model"].mass),
        killing_frequency=float(omega),
        occupations=np.asarray(occupations, dtype=float),
    )
    cap = 32.0 * ref["r_m_tot"]
    direction = {
        "dphi": np.asarray(phi) * _window(grid["radius"], 0.75 * cap, cap)[None, :, None],
        "d_occupations": np.zeros_like(state["occupations"]),
        "d_inner_mass": 0.0,
    }
    tangent = response.stage_tangent(state, grid, direction)
    return {
        fraction: chain.spatial_fee(
            grid["radius"], grid["weights"], tangent["d_lapse"], fraction * cap,
        )["l1_mass_response_fraction_inside"]
        for fraction in SEAM_FRACTIONS
    }


def main():
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    results = {}
    for points in POINT_SETS:
        phi, occupations, omega = template_probe(points)
        for eta in ETAS:
            results[("eta-{:.3f}".format(eta), points)] = fee(
                eta, points, phi, occupations, omega)
    tilts = {}
    for fraction in SEAM_FRACTIONS:
        tilt_876 = results[("eta-1.040", 876)][fraction] - results[("eta-1.010", 876)][fraction]
        tilt_1752 = results[("eta-1.040", 1752)][fraction] - results[("eta-1.010", 1752)][fraction]
        tilts[str(fraction)] = {
            "tilt_at_876": tilt_876,
            "tilt_at_1752": tilt_1752,
            "ratio_1752_over_876": tilt_1752 / tilt_876 if tilt_876 else None,
        }
    fraction_key = str(SEAM_FRACTIONS[0])
    ratio = tilts[fraction_key]["ratio_1752_over_876"]
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "results": {str(key): value for key, value in results.items()},
        "tilts": tilts,
        "verdict": (
            "GRID_ARTIFACT" if ratio is not None and ratio <= 0.6
            else ("FRINGE_PHYSICS" if ratio is not None and ratio >= 0.8 else "INCONCLUSIVE")
        ),
        "notes": [
            "the probe is the same interpolated 1.01 packet at both "
            "resolutions; only the chart grid density changes",
            "tilt = f(eta=1.04) - f(eta=1.01) at the inner seam",
        ],
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_charged_tilt_resolution.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_charged_tilt_resolution.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    for fraction, tilt in tilts.items():
        print("  seam {}: tilt876 {:.4f} tilt1752 {:.4f} ratio {!s}".format(
            fraction, tilt["tilt_at_876"], tilt["tilt_at_1752"],
            "none" if tilt["ratio_1752_over_876"] is None
            else "{:.2f}".format(tilt["ratio_1752_over_876"])))
    print("  verdict:", record["verdict"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())