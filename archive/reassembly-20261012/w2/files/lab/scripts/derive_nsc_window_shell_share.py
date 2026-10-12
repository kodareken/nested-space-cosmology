"""The window shell's causal weight in both instruments, both charts.

The containment chain wall < window shell [r_+, 2 r_+] < atlas cell 1
was established structurally (rounds 3, 22). Its quantitative link was
never measured: how much of each instrument's budget does the window
shell itself hold? This driver partitions both instruments at the window
edge (2 r_+) on the neutral mu=1.04 chart and on the charged admitted
slice, reporting the shell's share below the next atlas edge (8 r_m and
25% of the cap respectively).

    python scripts/lab.py scripts/derive_nsc_window_shell_share.py
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
from recursive_horizons import nsc_rn_charged_lll as lll
from recursive_horizons import nsc_rn_pg as pg
from recursive_horizons import nsc_rn_response as response
from recursive_horizons import nsc_rn_source as source

OUTPUT = LAB / "results" / "development" / "nsc-window-shell-share-v1.json"
SCHEMA = "NSC-WINDOW-SHELL-SHARE-v1"
A_CHARGED = 0.04501936182826115
C_F = 0.01125484045706527


def _window(radius, left, right, sharpness=8.0):
    return 0.5 * (1.0 + np.tanh(sharpness * (radius - left))) \
        * 0.5 * (1.0 + np.tanh(sharpness * (right - radius)))


def neutral_instrument_a():
    geometry = driver.case_geometry(1.040)
    grid = pg.build_grid(
        geometry["mass"], r_m=geometry["r_m"], points=geometry["route"]["points"],
        r_out=geometry["route"]["outer"],
    )
    packet = source.prepare_radial_packet(
        mass_over_rm=1.040, energy_target=1.0e-3, rank=driver.PILOT_RANK,
        model=grid["model"], radii=grid["radius"], weights=grid["weights"],
        derivative=grid["derivative"],
    )
    state = pg.make_state(
        grid, np.asarray(packet["phi"]), inner_mass=geometry["mass"],
        killing_frequency=float(packet["radial_moments"]["requested_omega"]),
        occupations=np.asarray(packet["occupations"], dtype=float),
    )
    radius = grid["radius"]
    direction = {
        "dphi": np.asarray(packet["phi"]) * _window(radius, 24.0, 32.0)[None, :, None],
        "d_occupations": np.zeros_like(state["occupations"]),
        "d_inner_mass": 0.0,
    }
    tangent = response.stage_tangent(state, grid, direction)
    model = source.action_reference(1.040)
    r_plus = float(model.horizons["r_plus"])
    edges = [float(np.min(radius)), 2.0 * r_plus, 8.0]
    shares = chain.interval_shares(radius, grid["weights"], tangent["d_lapse"], edges)
    return r_plus, shares


def charged_instrument_a():
    ref = lll.dyonic_reference(
        A=A_CHARGED, C_F=C_F, q_mag=1.0, Q_E=4.0, mass_over_total_radius=1.01,
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
    packet = source.prepare_radial_packet(
        mass_over_rm=1.01, energy_target=1.0e-3, rank=9,
        model=ref["model"], radii=grid["radius"], weights=grid["weights"],
        derivative=grid["derivative"],
    )
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
    r_plus = float(ref["horizons"]["r_plus"])
    edges = [float(np.min(grid["radius"])), 2.0 * r_plus, 0.25 * cap]
    shares = chain.interval_shares(grid["radius"], grid["weights"],
                                   tangent["d_lapse"], edges)
    return r_plus, cap, shares


def infall_shares(model, r_plus, inner_edge, cap):
    fine = np.linspace(inner_edge, cap, 2400)
    f = 1.0 - 2.0 * float(model.mass) / fine + float(model.magnetic_r2) / fine ** 2
    weight = 1.0 / (1.0 + np.sqrt(np.maximum(0.0, f)))
    total = float(np.trapezoid(weight, fine))
    shell = (fine >= r_plus) & (fine < 2.0 * r_plus)
    below_next = (fine >= 2.0 * r_plus) & (fine < (8.0 if cap == 32.0 else 0.25 * cap))
    return {
        "window_shell_share": float(np.trapezoid(weight[shell], fine[shell])) / total,
        "next_band_share": float(np.trapezoid(weight[below_next], fine[below_next])) / total,
    }


def main():
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    neutral_r_plus, neutral_a = neutral_instrument_a()
    charged_r_plus, charged_cap, charged_a = charged_instrument_a()
    neutral_model = source.action_reference(1.040)
    neutral_b = infall_shares(neutral_model, neutral_r_plus, neutral_r_plus, 32.0)
    ref = lll.dyonic_reference(
        A=A_CHARGED, C_F=C_F, q_mag=1.0, Q_E=4.0, mass_over_total_radius=1.01,
    )
    charged_b = infall_shares(ref["model"], charged_r_plus, charged_r_plus, charged_cap)
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "neutral_chart": {
            "r_plus": neutral_r_plus,
            "instrument_a": neutral_a,
            "instrument_b": neutral_b,
        },
        "charged_chart": {
            "r_plus": charged_r_plus,
            "cap": charged_cap,
            "instrument_a": charged_a,
            "instrument_b": charged_b,
        },
        "reading": (
            "window shell [r_+, 2 r_+] share of each instrument below the "
            "next atlas edge; the containment chain becomes quantitative"
        ),
        "notes": [
            "instrument A localizes the source in the outer band (>=24 r_m "
            "neutral, >=75% cap charged); B is the vacuum quadrature",
        ],
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_window_shell_share.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_window_shell_share.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    for chart in ("neutral_chart", "charged_chart"):
        entry = record[chart]
        a = entry["instrument_a"]
        b = entry["instrument_b"]
        print("  {:<14} A: shell {:.4f} next {:.4f}   B: shell {:.4f} next {:.4f}".format(
            chart, a[0]["l1_share"], a[1]["l1_share"],
            b["window_shell_share"], b["next_band_share"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())