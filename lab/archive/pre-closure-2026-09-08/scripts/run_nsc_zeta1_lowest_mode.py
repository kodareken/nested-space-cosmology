#!/usr/bin/env python3
"""Compute the first joined/disconnected ZETA1 Dirac spectral minimum."""

from __future__ import annotations

import hashlib
import json
from math import log, pi, sqrt
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.optimize import brentq


ROOT = Path(__file__).resolve().parents[1]
for value in (str(ROOT), str(ROOT / "src")):
    if value not in sys.path:
        sys.path.insert(0, value)

from recursive_horizons.evidence_io import canonical_json_bytes  # noqa: E402


OUTPUT = ROOT / "results/nsc-2-zeta1-lowest-mode.json"
INPUTS = (
    ROOT / "results/nsc-2-zeta1-foliation.json",
    ROOT / "results/nsc-2-zeta1-self-adjoint-domain.json",
    ROOT / "results/nsc-1-s-one-child-scale-correction.json",
)
CHILD_THRESHOLD = 3.0 * pi / 2.0
ROOT_UPPER = 10.0


def _authenticate() -> list[dict[str, object]]:
    authenticated = []
    for path in INPUTS:
        raw = path.read_bytes()
        item = json.loads(raw)
        if item.get("terminal") is not True:
            raise RuntimeError(f"nonterminal input: {path.relative_to(ROOT)}")
        authenticated.append(
            {
                "path": str(path.relative_to(ROOT)),
                "sha256": hashlib.sha256(raw).hexdigest(),
                "artifact_id": item["artifact_id"],
            }
        )
    return authenticated


def _partner_spectrum(
    box_radius: float,
    half_intervals: int,
    partner_sign: int,
) -> tuple[np.ndarray, np.ndarray]:
    spacing = box_radius / half_intervals

    def spectrum(left: float, right: float, count: int) -> np.ndarray:
        rho = np.linspace(left + spacing, right - spacing, count)
        superpotential = 1.0 / np.sqrt(1.0 + rho**2)
        derivative = -rho / (1.0 + rho**2) ** 1.5
        potential = superpotential**2 + partner_sign * derivative
        diagonal = 2.0 / spacing**2 + potential
        off_diagonal = np.full(count - 1, -1.0 / spacing**2)
        return eigh_tridiagonal(diagonal, off_diagonal, eigvals_only=True)

    joined = spectrum(
        -box_radius,
        box_radius,
        2 * half_intervals - 1,
    )
    parent = spectrum(-box_radius, 0.0, half_intervals - 1)
    child = spectrum(0.0, box_radius, half_intervals - 1)
    return joined, np.concatenate((parent, child))


def _complete_spectrum(
    box_radius: float, half_intervals: int
) -> tuple[np.ndarray, np.ndarray]:
    joined_parts = []
    disconnected_parts = []
    for partner_sign in (1, -1):
        joined, disconnected = _partner_spectrum(
            box_radius, half_intervals, partner_sign
        )
        joined_parts.append(joined)
        disconnected_parts.append(disconnected)
    return np.concatenate(joined_parts), np.concatenate(disconnected_parts)


def _spectral_quantities(
    zeta: float, eigenvalues: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    mu_squared = 1.0 - CHILD_THRESHOLD / zeta
    heat = np.exp(-(eigenvalues + mu_squared) / zeta)
    logarithmic_gradient = (
        (eigenvalues + 1.0) / zeta
        - 2.0 * CHILD_THRESHOLD / zeta**2
    )
    logarithmic_hessian = (
        -(eigenvalues + 1.0) / zeta
        + 4.0 * CHILD_THRESHOLD / zeta**2
    )
    first = heat * logarithmic_gradient
    second = heat * (
        logarithmic_gradient**2 + logarithmic_hessian
    )
    return heat, first, second


def _relative_values(
    zeta: float, joined: np.ndarray, disconnected: np.ndarray
) -> tuple[float, float, float]:
    joined_values = _spectral_quantities(zeta, joined)
    disconnected_values = _spectral_quantities(zeta, disconnected)
    return tuple(
        float(np.sum(left) - np.sum(right))
        for left, right in zip(joined_values, disconnected_values, strict=True)
    )


def _root_record(box_radius: float, half_intervals: int) -> dict[str, object]:
    joined, disconnected = _complete_spectrum(box_radius, half_intervals)

    def derivative(zeta: float) -> float:
        return _relative_values(zeta, joined, disconnected)[1]

    lower = CHILD_THRESHOLD + 1.0e-6
    lower_derivative = derivative(lower)
    upper_derivative = derivative(ROOT_UPPER)
    if not lower_derivative < 0.0 < upper_derivative:
        raise RuntimeError("lowest-mode scale derivative does not bracket a root")
    root = brentq(derivative, lower, ROOT_UPPER, xtol=1.0e-13, rtol=1.0e-13)
    heat, first, second = _relative_values(root, joined, disconnected)
    return {
        "box_radius": box_radius,
        "half_intervals": half_intervals,
        "spacing": box_radius / half_intervals,
        "joined_eigenvalue_count": int(joined.size),
        "disconnected_eigenvalue_count": int(disconnected.size),
        "minimum_joined_eigenvalue": float(np.min(joined)),
        "minimum_disconnected_eigenvalue": float(np.min(disconnected)),
        "lower_bracket_derivative": lower_derivative,
        "upper_bracket_derivative": upper_derivative,
        "zeta_root": root,
        "mu_squared_at_root": 1.0 - CHILD_THRESHOLD / root,
        "mu_at_root": sqrt(1.0 - CHILD_THRESHOLD / root),
        "relative_heat_trace": heat,
        "log_zeta_derivative": first,
        "log_zeta_hessian": second,
    }


def run() -> dict[str, object]:
    if OUTPUT.exists():
        raise RuntimeError("ZETA1 lowest-mode output already exists")
    authenticated = _authenticate()

    box_records = [
        _root_record(radius, int(round(radius / 0.04)))
        for radius in (20.0, 40.0, 60.0)
    ]
    resolution_records = [
        _root_record(40.0, intervals)
        for intervals in (500, 1000, 2000)
    ]
    coarse, medium, fine = resolution_records
    observed_order = log(
        abs(
            (coarse["zeta_root"] - medium["zeta_root"])
            / (medium["zeta_root"] - fine["zeta_root"])
        ),
        2.0,
    )
    extrapolated_root = (
        4.0 * fine["zeta_root"] - medium["zeta_root"]
    ) / 3.0
    extrapolation_error = abs(extrapolated_root - fine["zeta_root"])
    extrapolated_mu_squared = 1.0 - CHILD_THRESHOLD / extrapolated_root
    extrapolated_mu = sqrt(extrapolated_mu_squared)
    box_spread = max(item["zeta_root"] for item in box_records) - min(
        item["zeta_root"] for item in box_records
    )

    zeta_samples = []
    joined, disconnected = _complete_spectrum(40.0, 1000)
    for zeta in (5.0, 10.0, 20.0, 50.0, 100.0, 200.0):
        heat, first, second = _relative_values(zeta, joined, disconnected)
        zeta_samples.append(
            {
                "zeta": zeta,
                "mu_squared": 1.0 - CHILD_THRESHOLD / zeta,
                "relative_heat_trace": heat,
                "log_zeta_derivative": first,
                "log_zeta_hessian": second,
            }
        )

    record = {
        "artifact_id": "NSC-2-ZETA1-LOWEST-MODE",
        "schema": "NSC-2-ZETA1-LOWEST-MODE-v1",
        "classification": "the_lowest_joined_Dirac_partner_heat_spectrum_has_an_isolated_child_allowed_scale_minimum",
        "authenticated_inputs": authenticated,
        "method": {
            "angular_sector": "abs(kappa)=1",
            "y_sector": "unwarped_zero_mode_proxy",
            "squared_Dirac_partners": [
                "-d2/drho2+1/(1+rho2)-rho/(1+rho2)^(3/2)",
                "-d2/drho2+1/(1+rho2)+rho/(1+rho2)^(3/2)",
            ],
            "joined_domain": "rho_in[-R,R]_with_outer_Dirichlet",
            "disconnected_domain": "rho_in[-R,0]_direct_sum_[0,R]_with_throat_and_outer_Dirichlet",
            "profile": "exp[-(lambda_j+mu_squared(zeta))/zeta]",
            "child_constraint": "mu_squared(zeta)=1-3*pi/(2*zeta)",
            "stationarity": "d_DeltaK/d_log_zeta=0",
        },
        "box_convergence": box_records,
        "resolution_convergence": {
            "records": resolution_records,
            "observed_order": observed_order,
            "Richardson_zeta": extrapolated_root,
            "Richardson_error_estimate": extrapolation_error,
            "box_spread_at_spacing_0_04": box_spread,
        },
        "first_scale_estimate": {
            "zeta": extrapolated_root,
            "zeta_minus_child_threshold": extrapolated_root - CHILD_THRESHOLD,
            "mu_squared": extrapolated_mu_squared,
            "mu": extrapolated_mu,
            "minimum_is_above_3pi_over_2": extrapolated_root > CHILD_THRESHOLD,
            "minimum_Hessian": fine["log_zeta_hessian"],
        },
        "zeta_samples": zeta_samples,
        "next_result": "add_the_exact_y_warp_lapse_shift_and_higher_angular_sectors_then_compute_the_anomaly_consistent_relative_determinant",
        "gate": {
            "inputs_authenticated": len(authenticated) == len(INPUTS),
            "all_joined_and_disconnected_spectra_positive": all(
                item["minimum_joined_eigenvalue"] > 0.0
                and item["minimum_disconnected_eigenvalue"] > 0.0
                for item in box_records + resolution_records
            ),
            "all_roots_bracketed": all(
                item["lower_bracket_derivative"] < 0.0
                < item["upper_bracket_derivative"]
                for item in box_records + resolution_records
            ),
            "all_stationary_residuals_below_1e_10": all(
                abs(item["log_zeta_derivative"]) < 1.0e-10
                for item in box_records + resolution_records
            ),
            "all_stationary_Hessians_positive": all(
                item["log_zeta_hessian"] > 0.0
                for item in box_records + resolution_records
            ),
            "box_stability_below_1e_8": box_spread < 1.0e-8,
            "second_order_grid_convergence": observed_order > 1.9,
            "physical_side_minimum": extrapolated_root > CHILD_THRESHOLD,
            "full_transition_zeta_derived": False,
        },
        "nonclaims": {
            "lowest_sector_is_the_complete_spectrum": False,
            "y_warp_and_shift_included": False,
            "fermionic_determinant_and_anomaly_included": False,
            "extrapolated_zeta_is_the_physical_final_value": False,
        },
        "terminal": True,
    }
    for key, value in record["gate"].items():
        if key == "full_transition_zeta_derived":
            continue
        if value is not True:
            raise RuntimeError(f"ZETA1 lowest-mode gate failed: {key}")
    OUTPUT.write_bytes(canonical_json_bytes(record))
    return record


def main() -> int:
    record = run()
    print(json.dumps(record, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
