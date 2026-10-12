"""How much of the chart's budget does the cosmic-shell band hold?

Round-35 located the pre-registered local band at seams [25.9, 29.0]
r_m. This driver gives the location a weight: the same mu=1.04 tangent
probe is partitioned at edges 26 and 29 with instrument A (L1 share of
the lapse response) and instrument B (infall-time share of the vacuum
quadrature), and the band's share is compared with its length fraction
of the spanned route.

    python scripts/lab.py scripts/derive_nsc_cosmic_shell_weight.py
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

OUTPUT = LAB / "results" / "development" / "nsc-cosmic-shell-weight-v1.json"
SCHEMA = "NSC-COSMIC-SHELL-WEIGHT-v1"
BAND = (26.0, 29.0)
CAP = 32.0
MASS = 1.040
OUTER_START = 24.0


def _window(radius, left, right, sharpness=8.0):
    return 0.5 * (1.0 + np.tanh(sharpness * (radius - left))) \
        * 0.5 * (1.0 + np.tanh(sharpness * (right - radius)))


def main():
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    geometry = driver.case_geometry(MASS)
    grid = pg.build_grid(
        geometry["mass"], r_m=geometry["r_m"], points=geometry["route"]["points"],
        r_out=geometry["route"]["outer"],
    )
    packet = source.prepare_radial_packet(
        mass_over_rm=MASS, energy_target=1.0e-3, rank=driver.PILOT_RANK,
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
            radius, OUTER_START, float(geometry["route"]["outer"]))[None, :, None],
        "d_occupations": np.zeros_like(state["occupations"]),
        "d_inner_mass": 0.0,
    }
    tangent = response.stage_tangent(state, grid, direction)
    edges_a = [float(np.min(radius)), BAND[0], BAND[1], CAP]
    shares_a = chain.interval_shares(radius, weights, tangent["d_lapse"], edges_a)
    mass_shares = chain.interval_shares(radius, weights, tangent["d_mass"], edges_a)

    model = source.action_reference(MASS)
    r_plus_h = float(model.horizons["r_plus"])
    fine = np.linspace(r_plus_h, CAP, 2400)
    fval = 1.0 - 2.0 * float(model.mass) / fine + float(model.magnetic_r2) / fine ** 2
    weight_b = 1.0 / (1.0 + np.sqrt(np.maximum(0.0, fval)))
    total_b = float(np.trapezoid(weight_b, fine))
    band_inside = (fine >= BAND[0]) & (fine < BAND[1])
    share_b_band = float(np.trapezoid(weight_b[band_inside], fine[band_inside])) / total_b
    inside_below = (fine >= r_plus_h) & (fine < BAND[0])
    share_b_below = float(np.trapezoid(weight_b[inside_below], fine[inside_below])) / total_b
    length_fraction = (BAND[1] - BAND[0]) / (CAP - r_plus_h)
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "band": list(BAND),
        "instrument_a_lapse_shares": shares_a,
        "instrument_a_mass_shares": mass_shares,
        "instrument_b": {
            "band_share": share_b_band,
            "below_band_share": share_b_below,
            "above_band_share": 1.0 - share_b_band - share_b_below,
        },
        "band_length_fraction_of_span": length_fraction,
        "concentration": {
            "a_band_over_length": shares_a[1]["l1_share"] / length_fraction
            if length_fraction else None,
            "b_band_over_length": share_b_band / length_fraction if length_fraction else None,
        },
        "reading": (
            "the cosmic-shell band's causal weight in both instruments, "
            "normalized by its length fraction of the route"
        ),
        "notes": [
            "instrument A localizes the source above 24 as in the seam "
            "records; B is the vacuum infall quadrature of the same chart",
        ],
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_cosmic_shell_weight.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_cosmic_shell_weight.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    print("  A shares: below26 {:.4f}  band[26,29] {:.4f}  above29 {:.4f}".format(
        shares_a[0]["l1_share"], shares_a[1]["l1_share"], shares_a[2]["l1_share"]))
    print("  A mass:   below26 {:.4f}  band {:.4f}  above29 {:.4f}".format(
        mass_shares[0]["l1_share"], mass_shares[1]["l1_share"], mass_shares[2]["l1_share"]))
    print("  B shares: below26 {:.4f}  band {:.4f}  above29 {:.4f}".format(
        share_b_below, share_b_band, 1.0 - share_b_band - share_b_below))
    print("  band length fraction: {:.2%}  concentration A {:.2f} B {:.2f}".format(
        length_fraction, record["concentration"]["a_band_over_length"],
        record["concentration"]["b_band_over_length"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())