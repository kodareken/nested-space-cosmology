"""Twisted (kappa = sqrt(2)) charged column and its reflection scan.

First nonzero-angular charged multiplet on the admitted dyonic slice.
Solves the coupled characteristic system whose regular horizon root gives
both rows nonzero; the plus row is the reflected/outgoing channel. Scans
the band centre and four other declared frequencies and records the
reflection ratio at the outermost node, plus the horizon pivot ratio.

No time step; a mode-level scan on the extended grid. Output is
creation-only.

    python scripts/lab.py scripts/derive_nsc_rn_charged_twist.py
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
SIGN = 1
KAPPA = 2.0 ** 0.5
OMEGA_CENTRE = 0.8843561919019439
OMEGA_SCAN = (1.5, 1.2, OMEGA_CENTRE, 0.5, 0.25, 0.1)
OUTPUT = LAB / "results" / "development" / "nsc-rn-charged-twist-v1.json"
SCHEMA = "NSC-RN-CHARGED-TWIST-v1"


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
    rows = []
    for omega in OMEGA_SCAN:
        result = lll.charged_twisted_column(
            model=ref["model"], radius=grid["radius"], omega=omega, sign=SIGN,
            Q_E=Q_E, C_F=C_F, kappa=KAPPA,
        )
        ratio = lll.reflection_ratio(result["mode"], chart=chart)
        pivot = result["horizon_pivot"]
        weight = np.asarray(grid["weights"], dtype=float)
        norm = float(np.sqrt(np.real(np.sum(weight * np.sum(np.abs(result["mode"]) ** 2, axis=0)))))
        rows.append({
            "omega": float(omega),
            "omega_over_Phi_H": float(omega / ref["Phi_H"]),
            "reflection_ratio": ratio["reflection_ratio"],
            "outgoing_flux": ratio["outgoing_flux"],
            "ingoing_flux": ratio["ingoing_flux"],
            "horizon_plus_over_minus_modulus": float(abs(pivot["plus"])),
            "horizon_plus_over_minus": _plain(pivot["plus"]),
            "weighted_norm": float(norm),
            "ode_evaluations": int(result["evaluations"]),
        })
    record = {
        "schema": SCHEMA,
        "api_version": lll.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "slice": {
            "A": A, "C_F": C_F, "q_mag": 1.0, "Q_E": Q_E,
            "mass_over_total_radius": MASS_RATIO, "fermion_sign": SIGN,
            "kappa": KAPPA, "kappa_note": "first nonzero-angular charged multiplet, n=1",
            "Phi_H": ref["Phi_H"], "r_plus": ref["horizons"]["r_plus"],
            "r_minus": ref["horizons"]["r_minus"], "r_m_tot": ref["r_m_tot"],
            "grid_points": int(grid["points"]), "grid_outer": float(grid["outer"]),
        },
        "scan": rows,
        "notes": [
            "the reflection ratio is the outgoing/ingoing characteristic flux "
            "ratio at the outermost node; sub-unity means absorption into the "
            "hole, above unity would be a radiative window",
            "the horizon pivot ratio |chi_+/chi_-| at r_+ is the regular-root "
            "coupling alpha / |i E_eff - a_+'/2|",
        ],
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_rn_charged_lll.py": None,
            "lab/scripts/derive_nsc_rn_charged_twist.py": None,
        },
    }
    import hashlib
    for relative, target in (
        ("lab/src/recursive_horizons/nsc_rn_charged_lll.py",
         LAB / "src" / "recursive_horizons" / "nsc_rn_charged_lll.py"),
        ("lab/scripts/derive_nsc_rn_charged_twist.py",
         LAB / "scripts" / "derive_nsc_rn_charged_twist.py"),
    ):
        digest = hashlib.sha256()
        with target.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1 << 20), b""):
                digest.update(chunk)
        record["source_sha256"][relative] = digest.hexdigest()
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    for row in rows:
        print(
            "  omega={: .6f} (omega/Phi_H={: .4f}) reflection={: .6e}".format(
                row["omega"], row["omega_over_Phi_H"], row["reflection_ratio"],
            )
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())