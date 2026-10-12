"""Two-sign reflection control for the twisted kappa=sqrt(2) column.

Same admitted dyonic slice, same omega scan, both fermion signs. For
s = +1 the electric window opens below Phi_H and the round-5 scan already
showed the reflection rising toward the floor. For s = -1 the effective
frequency omega + Phi(r) is positive everywhere, so no window exists and
the reflection curve is the control that isolates the window signature.

The s = +1 arm is extended below the round-5 floor toward omega_min to
test whether any channel crosses into radiation (reflection above unity).

    python scripts/lab.py scripts/derive_nsc_rn_charged_twist_control.py
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

A = 0.04501936182826115
C_F = 0.01125484045706527
Q_E = 4.0
MASS_RATIO = 1.01
KAPPA = 2.0 ** 0.5
OMEGA_COMMON = (1.5, 1.2, 0.8843561919019439, 0.5, 0.25, 0.1)
OMEGA_FLOOR = (0.05, 0.0366)
OUTPUT = LAB / "results" / "development" / "nsc-rn-charged-twist-control-v1.json"
SCHEMA = "NSC-RN-CHARGED-TWIST-CONTROL-v1"


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


def main():
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    ref = lll.dyonic_reference(
        A=A, C_F=C_F, q_mag=1.0, Q_E=Q_E, mass_over_total_radius=MASS_RATIO,
    )
    grid = lll.charged_grid(model=ref["model"], r_m_tot=ref["r_m_tot"], points=876)
    chart = lll.charged_chart(ref["model"], grid["radius"])
    arms = {}
    for sign in (1, -1):
        omegas = OMEGA_COMMON + (OMEGA_FLOOR if sign == 1 else ())
        rows = []
        for omega in omegas:
            result = lll.charged_twisted_column(
                model=ref["model"], radius=grid["radius"], omega=omega, sign=sign,
                Q_E=Q_E, C_F=C_F, kappa=KAPPA,
            )
            ratio = lll.reflection_ratio(result["mode"], chart=chart)
            pivot = result["horizon_pivot"]
            rows.append({
                "omega": float(omega),
                "omega_over_Phi_H": float(omega / ref["Phi_H"]),
                "effective_frequency_at_horizon": float(pivot["effective_frequency_at_horizon"]),
                "reflection_ratio": ratio["reflection_ratio"],
                "outgoing_flux": float(ratio["outgoing_flux"]),
                "ingoing_flux": float(ratio["ingoing_flux"]),
                "radiative": bool(ratio["reflection_ratio"] > 1.0 + 1e-12) if ratio["reflection_ratio"] is not None else None,
            })
        arms["s_plus_1" if sign == 1 else "s_minus_1"] = rows
    record = {
        "schema": SCHEMA,
        "api_version": lll.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "slice": {
            "A": A, "C_F": C_F, "q_mag": 1.0, "Q_E": Q_E,
            "mass_over_total_radius": MASS_RATIO, "kappa": KAPPA,
            "Phi_H": ref["Phi_H"], "r_plus": ref["horizons"]["r_plus"],
            "omega_min": 0.0365869161951395,
            "omega_min_source": "window assessment record",
        },
        "arms": arms,
        "floor_points_included_only_for_s_plus_1": OMEGA_FLOOR,
        "notes": [
            "s = -1 keeps E_eff = omega + Phi(r) positive everywhere, so its "
            "arm is the no-window control",
            "radiative = reflection above unity at the outer node; the floor "
            "extension asks whether the round-5 rising trend ever crosses it",
        ],
    }
    import hashlib
    sources = {}
    for relative, target in (
        ("lab/src/recursive_horizons/nsc_rn_charged_lll.py",
         LAB / "src" / "recursive_horizons" / "nsc_rn_charged_lll.py"),
        ("lab/scripts/derive_nsc_rn_charged_twist_control.py",
         LAB / "scripts" / "derive_nsc_rn_charged_twist_control.py"),
    ):
        digest = hashlib.sha256()
        with target.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1 << 20), b""):
                digest.update(chunk)
        sources[relative] = digest.hexdigest()
    record["source_sha256"] = sources
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    for arm, rows in arms.items():
        for row in rows:
            print(
                "  {:<10} omega={: .6f} (omega/Phi_H={: .4f}) reflection={: .6e}{}".format(
                    arm, row["omega"], row["omega_over_Phi_H"], row["reflection_ratio"],
                    "  RADIATIVE" if row["radiative"] else "",
                )
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())