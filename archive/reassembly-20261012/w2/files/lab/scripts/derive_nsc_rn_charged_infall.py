"""First cap-admitted charged slice: kappa=0 mode and infall phase integral.

Consumes the best admitted row of nsc-rn-charged-window-assessment-v1.json
(audited A, C_F; q_mag=1, Q_E=4, mass ratio 1.01, sign +1), integrates the
kappa=0 lowest-Landau column at the admitted band centre omega_central on
the extended grid, and computes the exact infall winding quadrature:

    W(omega) = omega T_total - s (Q_E/(16 pi C_F)) I_1,
    omega_star = s Phi_mean,    Phi_mean = time-averaged horizon potential.

Below omega_star the entire infall winds opposite to the exterior Killing
phase. No time step is taken. Output is creation-only.

    python scripts/lab.py scripts/derive_nsc_rn_charged_infall.py
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

ASSESSMENT = LAB / "results" / "development" / "nsc-rn-charged-window-assessment-v1.json"
OUTPUT = LAB / "results" / "development" / "nsc-rn-charged-infall-v1.json"
SCHEMA = "NSC-RN-CHARGED-INFALL-v1"


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


def main():
    if not ASSESSMENT.is_file():
        raise SystemExit("missing window assessment record: " + str(ASSESSMENT))
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))

    assessment = json.loads(ASSESSMENT.read_text(encoding="utf-8"))
    row = assessment["verdict"]["best_admitted"]
    coefficients = assessment["coefficients"]
    A = float(coefficients["A"])
    C_F = float(coefficients["C_F"])
    q_mag = float(coefficients["scan_flux"])
    Q_E = float(row["Q_E"])
    eta = float(row["mass_over_total_radius"])
    sign = int(row["fermion_sign"])
    omega = float(row["omega_central"])
    neutral = float(row["neutral_frequency"])

    gate = window.window_gate(
        A=A, C_F=C_F, q_mag=q_mag, Q_E=Q_E, mass_over_total_radius=eta, sign=sign,
    )
    if not gate["admitted"]:
        raise SystemExit("the best admitted row does not open the gate on this checkout")
    ref = lll.dyonic_reference(A=A, C_F=C_F, q_mag=q_mag, Q_E=Q_E,
                               mass_over_total_radius=eta)
    grid = lll.charged_grid(model=ref["model"], r_m_tot=ref["r_m_tot"], points=876)
    chart = lll.charged_chart(ref["model"], grid["radius"])
    result = lll.charged_lll_mode(
        model=ref["model"], radius=grid["radius"], omega=omega, sign=sign,
        Q_E=Q_E, C_F=C_F,
    )
    normalized = lll.normalized_column(result["mode"], grid["weights"])
    form = lll.charged_mode_form(
        normalized["column"], chart=chart, weights=grid["weights"],
        derivative=grid["derivative"], radius=grid["radius"],
        Q_E=Q_E, C_F=C_F, sign=sign, omega=omega,
    )
    winding = lll.infall_winding(
        model=ref["model"], omega=omega, sign=sign, Q_E=Q_E, C_F=C_F,
        r_outer=float(grid["outer"]), points=2400,
    )
    winding_neutral = lll.infall_winding(
        model=ref["model"], omega=neutral, sign=sign, Q_E=Q_E, C_F=C_F,
        r_outer=float(grid["outer"]), points=2400,
    )
    edge_radius = Q_E / (16.0 * 3.141592653589793 * C_F * omega)
    charge = float(sign) * float(form["weighted_norm"]) ** 2

    record = {
        "schema": SCHEMA,
        "api_version": lll.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "slice": {
            "A": A,
            "C_F": C_F,
            "q_mag": q_mag,
            "Q_E": Q_E,
            "mass_over_total_radius": eta,
            "fermion_sign": sign,
            "omega": omega,
            "omega_source": "admitted band centre from the window assessment",
            "neutral_frequency": neutral,
            "neutral_frequency_in_window": neutral < gate["window_upper"],
            "Phi_H": gate["Phi_H"],
            "r_plus": gate["r_plus"],
            "r_minus": gate["r_minus"],
            "r_m_tot": gate["charge_squares"]["total_charge_radius"],
            "cap_room": gate["cap_room"],
            "window_upper": gate["window_upper"],
            "omega_min": gate["omega_min"],
            "gate_admitted": bool(gate["admitted"]),
        },
        "mode": {
            "kappa": result["kappa"],
            "plus_identically_zero": bool(result["plus_identically_zero"]),
            "nodes": int(result["node_count"]),
            "ode_evaluations": int(result["evaluations"]),
            "norm": form["weighted_norm"],
            "canonical_charge": charge,
            "response_energy": form["energy"],
            "window_entered": bool(form["window_entered"]),
            "window_edge_radius": edge_radius,
            "window_edge_over_r_plus": edge_radius / gate["r_plus"],
        },
        "infall_winding_at_band_centre": _plain(winding),
        "infall_winding_at_neutral_frequency": _plain(winding_neutral),
        "band_centre_below_omega_star": bool(omega < winding["omega_star"]),
        "neutral_frequency_below_omega_star": bool(neutral < winding_neutral["omega_star"]),
        "assessment_sha256": _path_sha(ASSESSMENT),
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_rn_charged_lll.py": _path_sha(
                LAB / "src" / "recursive_horizons" / "nsc_rn_charged_lll.py"
            ),
            "lab/scripts/derive_nsc_rn_charged_infall.py": _path_sha(
                LAB / "scripts" / "derive_nsc_rn_charged_infall.py"
            ),
        },
        "notes": [
            "phi_mean is a phase integral along the throughgoing characteristic; "
            "it is not the horizon potential and not a fit.",
            "below omega_star the infall phase winds opposite to the exterior "
            "Killing phase everywhere along the fall.",
        ],
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    print("slice: r_+={:.6f} Phi_H={:.6f} omega_central={:.6f} neutral={:.6f}".format(
        gate["r_plus"], gate["Phi_H"], omega, neutral))
    print("mode: charge={:+.6f} norm={:.6f} edge/r_+={:.6f}".format(
        charge, form["weighted_norm"], edge_radius / gate["r_plus"]))
    print("winding: T_total={:.6f} Phi_mean={:.6f} omega_star={:.6f}".format(
        winding["T_total"], winding["Phi_mean"], winding["omega_star"]))
    print("  band centre  -> phase {:.6f} ({}) below_star={}".format(
        winding["total_phase_at_omega"], winding["net_winding"],
        bool(omega < winding["omega_star"])))
    print("  neutral freq -> phase {:.6f} ({}) below_star={}".format(
        winding_neutral["total_phase_at_omega"], winding_neutral["net_winding"],
        bool(neutral < winding_neutral["omega_star"])))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())