#!/usr/bin/env python3
"""Evaluate the exact exponential-profile graviton heat-kernel form factor."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import dawsn
import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
for value in (str(ROOT), str(ROOT / "src")):
    if value not in sys.path:
        sys.path.insert(0, value)

from recursive_horizons.evidence_io import canonical_json_bytes  # noqa: E402


OUTPUT = ROOT / "results/nsc-1-s-one-finite-momentum-graviton.json"
PROFILE = ROOT / "results/nsc-1-s-one-spectral-profile-closure.json"
WALL = ROOT / "results/nsc-1-s-one-spectral-resolution-wall-bind.json"


def _master_h(value: float) -> float:
    if value == 0.0:
        return 1.0
    root = np.sqrt(value)
    return float(2.0 * dawsn(root / 2.0) / root)


def _graviton_form_factor(value: float) -> float:
    return (1.0 + value / 4.0) * _master_h(value) - 2.0


def run() -> dict[str, object]:
    if OUTPUT.exists():
        raise RuntimeError("finite-momentum graviton output already exists")

    authenticated = []
    for path in (PROFILE, WALL):
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

    z, h = sp.symbols("z h")
    f1 = (h - 1 + z / 6) / z**2
    f5 = -(h - 1) / (2 * z)
    source_form = -1 + z / 12 + z**2 * (f1 - f5 / 2)
    simplified = sp.simplify(source_form)
    expected = (1 + z / 4) * h - 2
    algebra_residual = sp.simplify(simplified - expected)

    alpha = sp.symbols("alpha", real=True)
    master_series = sum(
        (-z) ** order
        / sp.factorial(order)
        * sp.integrate(
            (alpha * (1 - alpha)) ** order,
            (alpha, 0, 1),
        )
        for order in range(6)
    )
    form_series = sp.series(
        expected.subs(h, master_series), z, 0, 6
    ).removeO()

    grid = np.unique(
        np.concatenate(
            (
                np.linspace(0.0, 100.0, 4097),
                np.logspace(-8, 8, 4097),
            )
        )
    )
    values = np.asarray([_graviton_form_factor(float(value)) for value in grid])
    maximum = minimize_scalar(
        lambda value: -_graviton_form_factor(float(value)),
        bounds=(0.0, 100.0),
        method="bounded",
        options={"xatol": 1.0e-14},
    )
    selected = [
        {
            "z": value,
            "h": _master_h(value),
            "F2": _graviton_form_factor(value),
        }
        for value in (0.0, 0.01, 0.1, 1.0, 2.0, 4.0, 8.0, 16.0, 100.0, 1.0e4, 1.0e8)
    ]

    record = {
        "artifact_id": "NSC-1-S-ONE-FINITE-MOMENTUM-GRAVITON",
        "schema": "NSC-1-S-ONE-FINITE-MOMENTUM-GRAVITON-v1",
        "classification": "the_fixed_heat_profile_has_a_bounded_pole_free_positive_Euclidean_axis_graviton_form_factor_connecting_local_derivatives_to_the_contact_wall",
        "authenticated_inputs": authenticated,
        "source": {
            "title": "High energy bosons do not propagate",
            "authors": ["M.A. Kurkov", "Fedele Lizzi", "Dmitri Vassilevich"],
            "arxiv": "1312.2235v1",
            "doi": "10.1016/j.physletb.2014.02.053",
            "equations": {
                "graviton_quadratic_heat_kernel": 14,
                "high_momentum_limit": 15,
                "Barvinsky_Vilkovisky_form_factors": [22, 23],
                "master_h": 24,
            },
        },
        "exact_algebra": {
            "source_expression": str(source_form),
            "simplified_F2": str(simplified),
            "expected_F2": str(expected),
            "simplification_residual": str(algebra_residual),
            "master_h": "integral_0^1 exp[-alpha*(1-alpha)*z] dalpha",
            "master_h_Dawson_form": "2*dawsn(sqrt(z)/2)/sqrt(z)",
            "low_z_series_through_fifth_order": str(form_series),
            "high_z_limit": "-3/2",
        },
        "positive_Euclidean_axis": {
            "grid_minimum_z": float(grid[0]),
            "grid_maximum_z": float(grid[-1]),
            "grid_count": int(grid.size),
            "minimum_F2": float(np.min(values)),
            "maximum_F2_on_grid": float(np.max(values)),
            "bounded_maximizer_z": float(maximum.x),
            "bounded_maximum_F2": float(-maximum.fun),
            "zero_found": bool(np.any(values == 0.0)),
            "all_sampled_values_negative": bool(np.all(values < 0.0)),
            "selected_values": selected,
        },
        "physical_interpretation": {
            "small_z": "ordinary_local_derivative_terms_appear_in_the_expansion",
            "large_z": "kernel_tends_to_a_constant_and_the_position_space_response_becomes_contact",
            "additional_positive_Euclidean_axis_graviton_pole": False,
            "finite_EGB_polynomial_extrapolation_replaced": True,
        },
        "next_result": "covariantize_the_exact_h_form_factor_on_the_two_sheet_curved_background_and_perform_the_retarded_Lorentzian_pole_and_spectral_positivity_test",
        "nonclaims": {
            "negative_F2_is_a_Lorentzian_ghost": False,
            "curved_background_propagator_solved": False,
            "Lorentzian_retarded_causality_proved": False,
            "particle_spectrum_predicted": False,
        },
        "gate": {
            "inputs_authenticated": len(authenticated) == 2,
            "source_algebra_simplifies_exactly": algebra_residual == 0,
            "low_z_linear_coefficient_one_twelfth": sp.expand(form_series).coeff(z, 1)
            == sp.Rational(1, 12),
            "positive_axis_scan_has_no_zero": bool(np.all(values < 0.0)),
            "high_z_samples_converge_to_minus_three_halves": abs(
                selected[-1]["F2"] + 1.5
            )
            < 1.0e-6,
            "curved_Lorentzian_test_completed": False,
        },
        "terminal": True,
    }
    for key, value in record["gate"].items():
        if key != "curved_Lorentzian_test_completed" and value is not True:
            raise RuntimeError(f"finite-momentum graviton gate failed: {key}")
    OUTPUT.write_bytes(canonical_json_bytes(record))
    return record


def main() -> int:
    record = run()
    print(json.dumps(record, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
