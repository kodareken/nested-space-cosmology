#!/usr/bin/env python3
"""Compute the angularly converged regulated-determinant ZETA1 minimum."""

from __future__ import annotations

import hashlib
import json
from math import log, pi, sqrt
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.optimize import brentq
from scipy.special import exp1


ROOT = Path(__file__).resolve().parents[1]
for value in (str(ROOT), str(ROOT / "src")):
    if value not in sys.path:
        sys.path.insert(0, value)

from recursive_horizons.evidence_io import canonical_json_bytes  # noqa: E402


OUTPUT = ROOT / "results/nsc-2-zeta1-regulated-determinant.json"
INPUTS = (
    ROOT / "results/nsc-2-zeta1-angular-tower.json",
    ROOT / "results/nsc-2-zeta1-self-adjoint-domain.json",
    ROOT / "results/nsc-1-s-one-invariant-partition-bind.json",
    ROOT / "results/nsc-1-s-one-spectral-profile-closure.json",
)
CHILD_THRESHOLD = 3.0 * pi / 2.0
ANGULAR_MAX = 16


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


def _sector_spectrum(
    box_radius: float, half_intervals: int, angular: int
) -> tuple[np.ndarray, np.ndarray]:
    spacing = box_radius / half_intervals
    joined_parts = []
    disconnected_parts = []

    def spectrum(
        left: float, right: float, count: int, partner_sign: int
    ) -> np.ndarray:
        rho = np.linspace(left + spacing, right - spacing, count)
        superpotential = angular / np.sqrt(1.0 + rho**2)
        derivative = -angular * rho / (1.0 + rho**2) ** 1.5
        potential = superpotential**2 + partner_sign * derivative
        diagonal = 2.0 / spacing**2 + potential
        off_diagonal = np.full(count - 1, -1.0 / spacing**2)
        return eigh_tridiagonal(diagonal, off_diagonal, eigvals_only=True)

    for partner_sign in (1, -1):
        joined_parts.append(
            spectrum(
                -box_radius,
                box_radius,
                2 * half_intervals - 1,
                partner_sign,
            )
        )
        disconnected_parts.extend(
            (
                spectrum(
                    -box_radius,
                    0.0,
                    half_intervals - 1,
                    partner_sign,
                ),
                spectrum(
                    0.0,
                    box_radius,
                    half_intervals - 1,
                    partner_sign,
                ),
            )
        )
    return np.concatenate(joined_parts), np.concatenate(disconnected_parts)


def _mode_values(zeta: float, values: np.ndarray) -> tuple[float, float, float]:
    quotient = (values + 1.0) / zeta - CHILD_THRESHOLD / zeta**2
    gradient_numerator = (
        (values + 1.0) / zeta - 2.0 * CHILD_THRESHOLD / zeta**2
    )
    gradient_numerator_y = (
        -(values + 1.0) / zeta + 4.0 * CHILD_THRESHOLD / zeta**2
    )
    action = 0.5 * exp1(quotient)
    exponential = np.exp(-quotient)
    first = 0.5 * exponential * gradient_numerator / quotient
    second = (
        0.5
        * exponential
        * (
            gradient_numerator**2 * (quotient + 1.0)
            + gradient_numerator_y * quotient
        )
        / quotient**2
    )
    return float(np.sum(action)), float(np.sum(first)), float(np.sum(second))


def _solve(
    box_radius: float, half_intervals: int, angular_max: int
) -> dict[str, object]:
    sectors = [
        _sector_spectrum(box_radius, half_intervals, angular)
        for angular in range(1, angular_max + 1)
    ]

    def relative(zeta: float) -> tuple[float, float, float]:
        total = np.zeros(3)
        for angular, (joined, disconnected) in enumerate(sectors, start=1):
            joined_values = _mode_values(zeta, joined)
            disconnected_values = _mode_values(zeta, disconnected)
            total += 4.0 * angular * (
                np.asarray(joined_values) - np.asarray(disconnected_values)
            )
        return tuple(float(value) for value in total)

    lower = CHILD_THRESHOLD + 1.0e-5
    upper = 8.0
    lower_derivative = relative(lower)[1]
    upper_derivative = relative(upper)[1]
    if not lower_derivative < 0.0 < upper_derivative:
        raise RuntimeError("regulated determinant does not bracket a scale root")
    root = brentq(
        lambda zeta: relative(zeta)[1],
        lower,
        upper,
        xtol=1.0e-13,
        rtol=1.0e-13,
    )
    action, derivative, hessian = relative(root)
    return {
        "box_radius": box_radius,
        "half_intervals": half_intervals,
        "spacing": box_radius / half_intervals,
        "angular_max": angular_max,
        "included_spinor_degeneracy": sum(
            4 * angular for angular in range(1, angular_max + 1)
        ),
        "lower_derivative": lower_derivative,
        "upper_derivative": upper_derivative,
        "zeta": root,
        "mu_squared": 1.0 - CHILD_THRESHOLD / root,
        "mu": sqrt(1.0 - CHILD_THRESHOLD / root),
        "relative_regulated_determinant": action,
        "log_zeta_derivative": derivative,
        "log_zeta_hessian": hessian,
    }


def run() -> dict[str, object]:
    if OUTPUT.exists():
        raise RuntimeError("ZETA1 regulated-determinant output already exists")
    authenticated = _authenticate()

    angular_records = [_solve(30.0, 750, value) for value in range(1, ANGULAR_MAX + 1)]
    box_records = [
        _solve(radius, int(round(radius / 0.04)), 12)
        for radius in (20.0, 30.0, 40.0, 60.0)
    ]
    resolution_records = [
        _solve(40.0, intervals, 12) for intervals in (500, 1000, 2000)
    ]

    coarse, medium, fine = resolution_records
    observed_order = log(
        abs(
            (coarse["zeta"] - medium["zeta"])
            / (medium["zeta"] - fine["zeta"])
        ),
        2.0,
    )
    grid_extrapolated = (4.0 * fine["zeta"] - medium["zeta"]) / 3.0
    finest_box = box_records[-1]["zeta"]
    medium_box = box_records[-2]["zeta"]
    box_correction = finest_box - medium_box
    best_zeta = grid_extrapolated + box_correction
    grid_error = abs(grid_extrapolated - fine["zeta"])
    box_error = abs(box_correction)
    combined_error = grid_error + box_error
    best_mu_squared = 1.0 - CHILD_THRESHOLD / best_zeta
    best_mu = sqrt(best_mu_squared)
    angular_tail = abs(
        angular_records[-1]["zeta"] - angular_records[-2]["zeta"]
    )

    record = {
        "artifact_id": "NSC-2-ZETA1-REGULATED-DETERMINANT",
        "schema": "NSC-2-ZETA1-REGULATED-DETERMINANT-v1",
        "classification": "the_exponentially_regulated_joined_Dirac_determinant_restores_an_angularly_converged_child_allowed_scale_minimum",
        "authenticated_inputs": authenticated,
        "functional": {
            "relative_action": "DeltaGamma_det=one_half*sum_joined E1[(lambda+mu_squared)/zeta]-one_half*sum_disconnected E1[(lambda+mu_squared)/zeta]",
            "proper_time_definition": "minus_log_det_Lambda=one_half*integral_(1/Lambda_squared)^infinity dt_over_t*Tr_exp(-t*H_squared)",
            "profile": "the_same_fixed_exponential_heat_regulator",
            "child_curve": "mu_squared=1-3*pi/(2*zeta)",
            "independent_bosonic_weight_added": False,
        },
        "angular_convergence": angular_records,
        "box_convergence": box_records,
        "grid_convergence": {
            "records": resolution_records,
            "observed_order": observed_order,
            "grid_Richardson_zeta": grid_extrapolated,
            "grid_error": grid_error,
            "box_correction": box_correction,
            "box_error": box_error,
            "angular_tail": angular_tail,
        },
        "current_scale_candidate": {
            "zeta": best_zeta,
            "error_enclosure": combined_error + angular_tail,
            "zeta_minus_3pi_over_2": best_zeta - CHILD_THRESHOLD,
            "mu_squared": best_mu_squared,
            "mu": best_mu,
            "physical_side": best_zeta > CHILD_THRESHOLD,
            "Hessian_positive": angular_records[-1]["log_zeta_hessian"] > 0.0,
        },
        "causal_result": {
            "heat_only_angular_nonpass_preserved": True,
            "regulated_fermionic_determinant_has_stable_root": True,
            "next_owner": "the_exact_y_warp_lapse_shift_and_local_compensating_anomaly",
            "next_result": "add_those_terms_and_require_the_root_to_remain_above_3pi_over_2_before_promoting_zeta",
        },
        "gate": {
            "inputs_authenticated": len(authenticated) == len(INPUTS),
            "all_roots_bracketed": all(
                item["lower_derivative"] < 0.0 < item["upper_derivative"]
                for item in angular_records + box_records + resolution_records
            ),
            "all_Hessians_positive": all(
                item["log_zeta_hessian"] > 0.0
                for item in angular_records + box_records + resolution_records
            ),
            "angular_tail_below_1e_8": angular_tail < 1.0e-8,
            "box_40_to_60_change_below_1e_6": abs(box_correction) < 1.0e-6,
            "second_order_grid_convergence": observed_order > 1.9,
            "candidate_above_child_threshold": best_zeta > CHILD_THRESHOLD,
            "full_warped_anomaly_consistent_zeta_derived": False,
        },
        "nonclaims": {
            "candidate_is_final_zeta": False,
            "y_warp_lapse_shift_included": False,
            "local_compensating_anomaly_included": False,
            "physical_gap_promoted": False,
        },
        "terminal": True,
    }
    for key, value in record["gate"].items():
        if key == "full_warped_anomaly_consistent_zeta_derived":
            continue
        if value is not True:
            raise RuntimeError(f"ZETA1 regulated determinant failed: {key}")
    OUTPUT.write_bytes(canonical_json_bytes(record))
    return record


def main() -> int:
    record = run()
    print(json.dumps(record, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
