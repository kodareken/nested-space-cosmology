#!/usr/bin/env python3
"""Locate the nearest complex zeros of the exact flat graviton form factor."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

import mpmath as mp
import numpy as np
from scipy.special import dawsn


ROOT = Path(__file__).resolve().parents[1]
for value in (str(ROOT), str(ROOT / "src")):
    if value not in sys.path:
        sys.path.insert(0, value)

from recursive_horizons.evidence_io import canonical_json_bytes  # noqa: E402


OUTPUT = ROOT / "results/nsc-1-s-one-flat-spectral-poles.json"
FORM_FACTOR = ROOT / "results/nsc-1-s-one-finite-momentum-graviton.json"


def _h_array(value):
    argument = np.asarray(value, dtype=complex)
    root = np.sqrt(argument)
    with np.errstate(divide="ignore", invalid="ignore"):
        answer = 2.0 * dawsn(root / 2.0) / root
    small = np.abs(argument) < 1.0e-12
    return np.where(
        small,
        1.0 - argument / 6.0 + argument**2 / 60.0,
        answer,
    )


def _f_array(value):
    argument = np.asarray(value, dtype=complex)
    return (1.0 + argument / 4.0) * _h_array(argument) - 2.0


def _h_mp(value):
    if abs(value) < mp.mpf("1e-30"):
        return 1 - value / 6 + value**2 / 60
    root = mp.sqrt(value)
    return mp.sqrt(mp.pi) * mp.exp(-value / 4) * mp.erfi(root / 2) / root


def _f_mp(value):
    return (1 + value / 4) * _h_mp(value) - 2


def _winding(radius: float, count: int = 32768) -> dict[str, object]:
    angles = np.linspace(0.0, 2.0 * np.pi, count + 1)
    values = _f_array(radius * np.exp(1j * angles))
    phases = np.unwrap(np.angle(values))
    winding_float = float((phases[-1] - phases[0]) / (2.0 * np.pi))
    return {
        "radius": radius,
        "sample_count": count,
        "winding_float": winding_float,
        "winding_integer": int(round(winding_float)),
        "minimum_boundary_modulus": float(np.min(np.abs(values))),
    }


def run() -> dict[str, object]:
    if OUTPUT.exists():
        raise RuntimeError("flat spectral-pole output already exists")

    source_raw = FORM_FACTOR.read_bytes()
    source = json.loads(source_raw)
    if source.get("terminal") is not True:
        raise RuntimeError("finite-momentum form factor is not terminal")

    mp.mp.dps = 60
    upper = mp.findroot(
        _f_mp,
        (mp.mpc(-1, 6), mp.mpc(-1.5, 7)),
        tol=mp.mpf("1e-50"),
    )
    lower = mp.conj(upper)
    roots = (lower, upper) if lower.imag < upper.imag else (upper, lower)
    windings = [_winding(radius) for radius in (1.0, 6.8, 6.9, 10.0)]

    real_axis = np.linspace(-100.0, 100.0, 40001)
    real_values = _f_array(real_axis)
    maximum_real_axis_imaginary = float(np.max(np.abs(real_values.imag)))
    minimum_real_axis_modulus = float(np.min(np.abs(real_values)))

    record = {
        "artifact_id": "NSC-1-S-ONE-FLAT-SPECTRAL-POLES",
        "schema": "NSC-1-S-ONE-FLAT-SPECTRAL-POLES-v1",
        "classification": "the_fixed_full_graviton_form_factor_has_no_additional_pole_inside_the_unit_spectral_cutoff_disk",
        "authenticated_input": {
            "path": str(FORM_FACTOR.relative_to(ROOT)),
            "sha256": hashlib.sha256(source_raw).hexdigest(),
            "artifact_id": source.get("artifact_id"),
        },
        "analytic_continuation": {
            "form_factor": "F2(z)=(1+z/4)*integral_0^1 exp[-a*(1-a)*z] da-2",
            "entire": True,
            "Euclidean_spacelike_axis": "z>=0",
            "Lorentzian_timelike_axis": "z<=0_after_flat_Wick_continuation",
        },
        "nearest_complex_zeros": [
            {
                "real": mp.nstr(root.real, 30),
                "imaginary": mp.nstr(root.imag, 30),
                "modulus": mp.nstr(abs(root), 30),
                "residual": mp.nstr(abs(_f_mp(root)), 8),
            }
            for root in roots
        ],
        "argument_principle": windings,
        "real_axis_scan": {
            "minimum": -100.0,
            "maximum": 100.0,
            "sample_count": int(real_axis.size),
            "maximum_imaginary_roundoff": maximum_real_axis_imaginary,
            "minimum_modulus": minimum_real_axis_modulus,
            "zero_found": bool(np.any(np.abs(real_values) < 1.0e-10)),
        },
        "physical_result": {
            "zeros_inside_abs_z_one": 0,
            "additional_pole_inside_local_cutoff": False,
            "nearest_zero_modulus": float(abs(upper)),
            "nearest_zeros_beyond_cutoff": abs(upper) > 1,
            "real_timelike_extra_mass_pole_found": False,
            "real_spacelike_instability_pole_found": False,
        },
        "scope": {
            "flat_background": True,
            "full_curved_background_Hessian": False,
            "retarded_sheet_and_branch_cut_prescription_completed": False,
            "physical_massless_room_graviton_included": False,
        },
        "next_result": "evaluate_the_same_argument_principle_for_the_background_adjusted_two_sheet_retarded_Hessian_below_the_local_cutoff",
        "gate": {
            "input_authenticated": source.get("terminal") is True,
            "root_pair_residual_below_1e_40": all(
                abs(_f_mp(root)) < mp.mpf("1e-40") for root in roots
            ),
            "unit_disk_winding_zero": windings[0]["winding_integer"] == 0,
            "no_zeros_inside_radius_6_8": windings[1]["winding_integer"] == 0,
            "two_zeros_inside_radius_6_9": windings[2]["winding_integer"] == 2,
            "real_axis_has_no_zero": not bool(
                np.any(np.abs(real_values) < 1.0e-10)
            ),
            "curved_retarded_test_completed": False,
        },
        "terminal": True,
    }
    for key, value in record["gate"].items():
        if key != "curved_retarded_test_completed" and value is not True:
            raise RuntimeError(f"flat spectral-pole gate failed: {key}")
    OUTPUT.write_bytes(canonical_json_bytes(record))
    return record


def main() -> int:
    record = run()
    print(json.dumps(record, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
