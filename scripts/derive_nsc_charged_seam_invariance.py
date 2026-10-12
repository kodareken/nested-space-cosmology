"""Charged seam-profile invariance across eta (scaffold parity test).

Round-33 recorded the neutral seam profile as chart-invariant at 2.1e-3.
This driver runs the same tangent probe on the dyonic charts eta =
1.01, 1.02, 1.04 at cap-fraction seams 12.5/25/50/75% and measures the
pairwise maximum deviation. If the outside scaffold is one object, the
charged profile is eta-invariant like its neutral twin.

    python scripts/lab.py scripts/derive_nsc_charged_seam_invariance.py
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

OUTPUT = LAB / "results" / "development" / "nsc-charged-seam-invariance-v1.json"
SCHEMA = "NSC-CHARGED-SEAM-INVARIANCE-v1"
A = 0.04501936182826115
C_F = 0.01125484045706527
ETAS = (1.01, 1.02, 1.04)
SEAM_FRACTIONS = (0.125, 0.25, 0.5, 0.75)


def _window(radius, left, right, sharpness=8.0):
    return 0.5 * (1.0 + np.tanh(sharpness * (radius - left))) \
        * 0.5 * (1.0 + np.tanh(sharpness * (right - radius)))


def probe(eta, packet_template):
    ref = lll.dyonic_reference(
        A=A, C_F=C_F, q_mag=1.0, Q_E=4.0, mass_over_total_radius=eta,
    )
    grid = lll.charged_grid(model=ref["model"], r_m_tot=ref["r_m_tot"], points=876)
    spacing = float(grid["spacing"])
    radius = np.asarray(grid["inner_edge"]) + spacing * np.arange(876, dtype=float)
    weights, derivative, q_matrix = pg.sbp_operator(876, spacing, layout="csr")
    grid = dict(grid)
    grid.update({
        "radius": radius, "weights": weights, "derivative": derivative,
        "q_matrix": q_matrix, "outer": float(radius[-1]), "charge": 0.0,
    })
    state = pg.make_state(
        grid, np.asarray(packet_template["phi"]), inner_mass=float(ref["model"].mass),
        killing_frequency=float(packet_template["radial_moments"]["requested_omega"]),
        occupations=np.asarray(packet_template["occupations"], dtype=float),
    )
    cap = 32.0 * ref["r_m_tot"]
    direction = {
        "dphi": np.asarray(packet_template["phi"]) * _window(
            grid["radius"], 0.75 * cap, cap)[None, :, None],
        "d_occupations": np.zeros_like(state["occupations"]),
        "d_inner_mass": 0.0,
    }
    tangent = response.stage_tangent(state, grid, direction)
    profile = {}
    for fraction in SEAM_FRACTIONS:
        seam = fraction * cap
        profile[fraction] = chain.spatial_fee(
            grid["radius"], grid["weights"], tangent["d_lapse"], seam,
        )["l1_mass_response_fraction_inside"]
    return profile


def main():
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    template_ref = lll.dyonic_reference(
        A=A, C_F=C_F, q_mag=1.0, Q_E=4.0, mass_over_total_radius=1.01,
    )
    template_grid = lll.charged_grid(
        model=template_ref["model"], r_m_tot=template_ref["r_m_tot"], points=876)
    template_spacing = float(template_grid["spacing"])
    template_radius = np.asarray(template_grid["inner_edge"]) \
        + template_spacing * np.arange(876, dtype=float)
    weights, derivative, _ = pg.sbp_operator(876, template_spacing, layout="csr")
    template_packet = source.prepare_radial_packet(
        mass_over_rm=1.01, energy_target=1.0e-3, rank=9,
        model=template_ref["model"], radii=template_radius, weights=weights,
        derivative=derivative,
    )
    profiles = {}
    for eta in ETAS:
        profiles["eta-{:.3f}".format(eta)] = probe(eta, template_packet)
    deviations = []
    keys = list(profiles)
    for index in range(len(keys)):
        for other in range(index + 1, len(keys)):
            a = profiles[keys[index]]
            b = profiles[keys[other]]
            deviations.append({
                "pair": [keys[index], keys[other]],
                "max_abs_deviation": max(
                    abs(a[fraction] - b[fraction]) for fraction in SEAM_FRACTIONS),
            })
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "seam_fractions": list(SEAM_FRACTIONS),
        "profiles": profiles,
        "pairwise": deviations,
        "max_pairwise_deviation": max(item["max_abs_deviation"] for item in deviations),
        "verdict": (
            "PROFILES_COINCIDE"
            if max(item["max_abs_deviation"] for item in deviations) < 0.01
            else "PROFILES_DIVERGE"
        ),
        "notes": [
            "cap-fraction seams make the comparison scale-free across eta",
            "the source band is the outer quarter of each cap",
            "single fixed 1.01 probe field in fractional node coordinates "
            "on every chart (the packet preparer is family-pinned); the "
            "charts themselves differ",
        ],
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_charged_seam_invariance.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_charged_seam_invariance.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    for key, profile in profiles.items():
        print("  {:<10} ".format(key) + "  ".join(
            "f({:.0%})={:.4f}".format(fraction, profile[fraction])
            for fraction in SEAM_FRACTIONS))
    print("  max pairwise deviation: {:.4f}   verdict: {}".format(
        record["max_pairwise_deviation"], record["verdict"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())