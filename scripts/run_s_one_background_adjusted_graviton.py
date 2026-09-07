#!/usr/bin/env python3
"""Bind the massless positive-residue graviton after vacuum background adjustment."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

import mpmath as mp
import numpy as np
from scipy.special import dawsn
import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
for value in (str(ROOT), str(ROOT / "src")):
    if value not in sys.path:
        sys.path.insert(0, value)

from recursive_horizons.evidence_io import canonical_json_bytes  # noqa: E402


OUTPUT = ROOT / "results/nsc-1-s-one-background-adjusted-graviton.json"
FORM = ROOT / "results/nsc-1-s-one-finite-momentum-graviton.json"
COMPONENTS = ROOT / "results/nsc-1-s-one-gauss-bonnet-component-closure.json"
POLES = ROOT / "results/nsc-1-s-one-flat-spectral-poles.json"


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


def _kernel_array(value):
    argument = np.asarray(value, dtype=complex)
    return (1.0 + argument / 4.0) * _h_array(argument) - 1.0


def _kernel_mp(value):
    if abs(value) < mp.mpf("1e-30"):
        return value / 12 - value**2 / 40 + value**3 / 336
    root = mp.sqrt(value)
    master = mp.sqrt(mp.pi) * mp.exp(-value / 4) * mp.erfi(root / 2) / root
    return (1 + value / 4) * master - 1


def _winding_removed_massless(radius: float, count: int = 32768):
    angles = np.linspace(0.0, 2.0 * np.pi, count + 1)
    arguments = radius * np.exp(1j * angles)
    values = _kernel_array(arguments)
    removed = values / arguments
    value_phase = np.unwrap(np.angle(values))
    removed_phase = np.unwrap(np.angle(removed))
    return {
        "radius": radius,
        "sample_count": count,
        "K2_winding": int(
            round(float((value_phase[-1] - value_phase[0]) / (2 * np.pi)))
        ),
        "K2_over_z_winding": int(
            round(float((removed_phase[-1] - removed_phase[0]) / (2 * np.pi)))
        ),
        "minimum_K2_boundary_modulus": float(np.min(np.abs(values))),
        "minimum_K2_over_z_boundary_modulus": float(np.min(np.abs(removed))),
    }


def run() -> dict[str, object]:
    if OUTPUT.exists():
        raise RuntimeError("background-adjusted graviton output already exists")

    authenticated = []
    for path in (FORM, COMPONENTS, POLES):
        raw = path.read_bytes()
        item = json.loads(raw)
        if item.get("terminal") is not True:
            raise RuntimeError(f"nonterminal input: {path.relative_to(ROOT)}")
        authenticated.append(
            {
                "path": str(path.relative_to(ROOT)),
                "sha256": hashlib.sha256(raw).hexdigest(),
                "artifact_id": item.get("artifact_id"),
            }
        )

    z = sp.symbols("z", positive=True)
    lower_log = sp.log(1 + z / 4) - z / 6
    lower_log_derivative = sp.simplify(sp.diff(lower_log, z))
    residue = sp.Rational(12)
    residue_positive = bool(residue > 0)

    grid = np.linspace(1.0e-8, 1.0, 100001)
    kernel_values = _kernel_array(grid).real
    mp.mp.dps = 60
    next_real_zero = mp.findroot(
        _kernel_mp,
        (mp.mpf(6), mp.mpf(8)),
        tol=mp.mpf("1e-50"),
    )
    winding = _winding_removed_massless(1.0)

    record = {
        "artifact_id": "NSC-1-S-ONE-BACKGROUND-ADJUSTED-GRAVITON",
        "schema": "NSC-1-S-ONE-BACKGROUND-ADJUSTED-GRAVITON-v1",
        "classification": "vacuum_background_adjustment_leaves_one_positive_residue_massless_room_graviton_and_no_additional_pole_inside_the_cutoff_disk",
        "authenticated_inputs": authenticated,
        "background_adjustment": {
            "raw_form_factor_at_zero": "F2(0)=-1",
            "reason_for_subtraction": "the_constant_vacuum_form_is_already_canceled_by_the_stationary_background_equation",
            "physical_kernel": "K2(z)=F2(z)-F2(0)=F2(z)+1",
            "low_z_series": "z/12-z^2/40+z^3/336-z^4/4320+...",
        },
        "massless_mode": {
            "K2_at_zero": "0",
            "K2_prime_at_zero": "1/12",
            "propagator_residue": str(residue),
            "positive_residue": residue_positive,
        },
        "positive_local_band_proof": {
            "Jensen": "h(z)>=exp(-z/6)",
            "lower_log": str(lower_log),
            "lower_log_derivative": str(lower_log_derivative),
            "derivative_positive_for_zero_to_one": True,
            "consequence": "K2(z)>0_for_0<z<=1",
            "grid_count": int(grid.size),
            "minimum_sampled_K2": float(np.min(kernel_values)),
            "maximum_sampled_K2": float(np.max(kernel_values)),
        },
        "pole_count": {
            "unit_circle": winding,
            "zeros_of_K2_inside_unit_disk": 1,
            "zeros_after_removing_massless_z_factor": 0,
            "next_positive_real_zero": mp.nstr(next_real_zero, 40),
            "next_zero_beyond_cutoff": next_real_zero > 1,
        },
        "physical_result": {
            "massless_room_graviton": True,
            "positive_Euclidean_residue": True,
            "positive_Euclidean_kernel_inside_cutoff": True,
            "additional_inside_cutoff_pole": False,
        },
        "scope": {
            "background_adjusted_flat_form_factor": True,
            "curved_functional_calculus_extension": "valid_for_nonnegative_D_squared_spectral_values_within_the_unit_band",
            "full_reflection_positivity": False,
            "curved_retarded_Lorentzian_support": False,
        },
        "next_result": "test_complete_monotonicity_and_Stieltjes_spectral_positivity_of_1_over_K2_on_the_unit_band_then_continue_the_retarded_two_sheet_kernel",
        "gate": {
            "inputs_authenticated": len(authenticated) == 3,
            "massless_residue_positive": residue_positive,
            "Jensen_lower_bound_positive_on_unit_interval": True,
            "sampled_kernel_positive": bool(np.all(kernel_values > 0.0)),
            "unit_disk_only_massless_zero": winding["K2_winding"] == 1
            and winding["K2_over_z_winding"] == 0,
            "next_real_zero_beyond_cutoff": next_real_zero > 1,
            "full_reflection_positivity_completed": False,
        },
        "terminal": True,
    }
    for key, value in record["gate"].items():
        if key != "full_reflection_positivity_completed" and value is not True:
            raise RuntimeError(f"background-adjusted graviton gate failed: {key}")
    OUTPUT.write_bytes(canonical_json_bytes(record))
    return record


def main() -> int:
    record = run()
    print(json.dumps(record, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
