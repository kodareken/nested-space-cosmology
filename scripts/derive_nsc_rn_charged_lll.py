"""Charged kappa=0 lowest-Landau preflight on the extended dyonic grid.

Builds the note's manufactured dyonic slice (A, C_F, q_mag=1, Q_E=4,
mass ratio 2), integrates the future-horizon regular kappa=0 column from
one cell above the inner horizon out to the declared parent cap
32 r_m,tot, normalizes it in the weighted Gram, and records:

    - the window gate reading (window open, cap admission false, as the
      note already declared for this sample);
    - the plus-identically-zero property and its one admitted rate (the
      horizon jet derivative);
    - the unit canonical charge carried by the normalized column;
    - the Killing-form energy, probability current, and the radial shell
      where omega - s Phi(r) < 0 lives (the window edge radius).

No time step is taken and no probe is evolved against the metric. Output
is creation-only.

    python scripts/lab.py scripts/derive_nsc_rn_charged_lll.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

LAB = Path(__file__).resolve().parents[1]
if str(LAB / "src") not in sys.path:
    sys.path.insert(0, str(LAB / "src"))

import numpy as np

from recursive_horizons import nsc_rn_charged_lll as lll  # noqa: E402
from recursive_horizons import nsc_rn_charged_window as window  # noqa: E402

A = 0.125
C_F = 0.5
Q_MAG = 1.0
Q_E = 4.0
MASS_RATIO = 2.0
SIGN = 1
PHI_H = 0.04063643378570953
OMEGA = PHI_H / 2.0
OUTPUT = LAB / "results" / "development" / "nsc-rn-charged-lll-preflight-v1.json"
SCHEMA = "NSC-RN-CHARGED-LLL-PREFLIGHT-v1"


def _plain(value):
    if isinstance(value, complex):
        return {"real": float(value.real), "imag": float(value.imag)}
    if isinstance(value, np.floating):
        return float(value)
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.ndarray):
        return [_plain(item) for item in value.tolist()]
    if isinstance(value, (list, tuple)):
        return [_plain(item) for item in value]
    if isinstance(value, dict):
        return {key: _plain(item) for key, item in value.items()}
    return value


def _path_sha(path):
    import hashlib
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def field_energy_normals(column):
    """Tags for the recorded scalars without inventing a full stress."""
    weight = np.abs(column[1]) ** 2
    return {
        "minus_carried_total": float(np.sum(weight)),
        "plus_carried_total": float(np.sum(np.abs(column[0]) ** 2)),
        "charge_row_nonzero": bool(np.any(weight > 0.0)),
    }


def window_edge_radius(phi_h, omega, q_e):
    """The radius where the effective frequency passes through zero."""
    import math as _math
    if omega <= 0.0 or phi_h <= 0.0:
        return None
    edge = q_e / (16.0 * _math.pi * C_F * omega)
    return None if not _math.isfinite(edge) else float(edge)


def main():
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))

    gate = window.window_gate(
        A=A, C_F=C_F, q_mag=Q_MAG, Q_E=Q_E,
        mass_over_total_radius=MASS_RATIO, sign=SIGN,
    )
    ref = lll.dyonic_reference(
        A=A, C_F=C_F, q_mag=Q_MAG, Q_E=Q_E, mass_over_total_radius=MASS_RATIO,
    )
    grid = lll.charged_grid(model=ref["model"], r_m_tot=ref["r_m_tot"], points=876)
    chart = lll.charged_chart(ref["model"], grid["radius"])
    result = lll.charged_lll_mode(
        model=ref["model"], radius=grid["radius"], omega=OMEGA, sign=SIGN,
        Q_E=Q_E, C_F=C_F,
    )
    normalized = lll.normalized_column(result["mode"], grid["weights"])
    column = normalized["column"]
    form = lll.charged_mode_form(
        column, chart=chart, weights=grid["weights"], derivative=grid["derivative"],
        radius=grid["radius"], Q_E=Q_E, C_F=C_F, sign=SIGN, omega=OMEGA,
    )
    support = field_energy_normals(column)
    charge = SIGN * form["weighted_norm"] ** 2
    edge_radius = window_edge_radius(PHI_H, OMEGA, Q_E)

    settings = {
        "window_gate": _plain(gate),
        "mode": {
            "kappa": result["kappa"],
            "plus_identically_zero": bool(result["plus_identically_zero"]),
            "nodes": int(result["node_count"]),
            "nodes_dropped": int(result["nodes_dropped"]),
            "ode_evaluations": int(result["evaluations"]),
            "resampling": bool(result["resampling"]),
            "horizon_jet": _plain(result["horizon_jet"]),
        },
        "one_admitted_rate": {
            "definition": (
                "the only admitted rate is the horizon-jet minus derivative; "
                "the kappa=0 system has no metric feedback stage"
            ),
            "effective_frequency_at_horizon": _plain(result["horizon_jet"]["effective_frequency_at_horizon"]),
            "minus_derivative_at_horizon": _plain(result["horizon_jet"]["minus_derivative"]),
            "phi_h": float(PHI_H),
            "omega": float(OMEGA),
        },
        "form": _plain(form),
        "gauss": {
            "canonical_charge": float(charge),
            "sign_convention": "Q_E increases outward; positive canonical charge forces A_t < 0",
            "charge_row_nonzero": bool(support["charge_row_nonzero"]),
            "minus_carried_total": float(support["minus_carried_total"]),
            "plus_carried_total": float(support["plus_carried_total"]),
        },
        "window_edge": {
            "edge_radius": edge_radius,
            "edge_in_r_m_tot_units": None if edge_radius is None else edge_radius / ref["r_m_tot"],
            "reading": (
                "the B_< window occupies the finite radial shell from r_+ to "
                "the radius where omega - s Phi(r) changes sign"
            ),
        },
        "grid": {
            "points": int(grid["points"]),
            "inner_edge": _plain(grid["inner_edge"]),
            "inner_minus_r_hminus": _plain(grid["inner_edge"] - ref["horizons"]["r_minus"]),
            "outer": _plain(grid["outer"]),
            "cap_multiple": _plain(grid["cap_multiple"]),
            "spacing": _plain(grid["spacing"]),
        },
    }
    record = {
        "schema": SCHEMA,
        "api_version": lll.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "admitted_gate": bool(gate["admitted"]),
        "window_open": bool(gate["window_open"]),
        "settings": settings,
        "input_sha256": {
            "lab/src/recursive_horizons/nsc_rn_charged_lll.py": _path_sha(
                LAB / "src" / "recursive_horizons" / "nsc_rn_charged_lll.py"
            ),
            "lab/src/recursive_horizons/nsc_rn_charged_window.py": _path_sha(
                LAB / "src" / "recursive_horizons" / "nsc_rn_charged_window.py"
            ),
            "lab/docs/nsc-rn-charged-source.md": _path_sha(
                LAB / "docs" / "nsc-rn-charged-source.md"
            ),
        },
        "notes": [
            "The manufactured sample itself is not cap-admitted (window open, "
            "budget width too narrow) exactly as the note declares; admitting "
            "production physics requires a changed declaration, not this file.",
            "The 'one admitted rate' is the horizon jet derivative. The charged "
            "branch has no metric-backreaction stage yet, so there is nothing "
            "like the neutral stage_rates to call and no field is evolved.",
        ],
    }
    written = json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n"
    OUTPUT.write_text(written, encoding="utf-8")
    print("wrote", OUTPUT.relative_to(LAB))
    print(
        "  slices: r_+={:.6f}  edge={}".format(
            ref["horizons"]["r_plus"],
            None if edge_radius is None else round(edge_radius, 6),
        )
    )
    print(
        "  mode: kappa={} plus_zero={} nodes={} ode_eval={} norm={:.6f}".format(
            result["kappa"], result["plus_identically_zero"], result["node_count"],
            result["evaluations"], form["weighted_norm"],
        )
    )
    print(
        "  charge={:+.3f}  response_energy={:.6f}  edge_in_r_m={}".format(
            charge, form["energy"],
            None if edge_radius is None else round(edge_radius / ref["r_m_tot"], 4),
        )
    )
    print(
        "  gate: window_open={} admitted={} reason={}".format(
            gate["window_open"], gate["admitted"], gate["reason"],
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())