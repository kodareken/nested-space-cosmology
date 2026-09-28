#!/usr/bin/env python3
"""Measure the actual-background EGB time-principal kinetic rank."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
for value in (str(ROOT), str(ROOT / "src")):
    if value not in sys.path:
        sys.path.insert(0, value)

from recursive_horizons.evidence_io import canonical_json_bytes  # noqa: E402


OUTPUT = ROOT / "results/nsc-1-s-one-gauss-bonnet-kinetic-rank.json"
GEOMETRY = ROOT / "results/nsc-1-s-one-two-sheet-5d-geometry.json"
CRITICALITY = ROOT / "results/nsc-1-s-one-gauss-bonnet-criticality.json"


def _background_function():
    time, rho, theta, phi, outside = sp.symbols(
        "t rho theta phi y", real=True
    )
    coordinates = (time, rho, theta, phi, outside)
    warp_curvature = sp.Rational(18, 1015)
    conformal_factor = sp.exp(-2 * warp_curvature * outside**2)
    radius = sp.sqrt(rho**2 + 1)
    metric_b = (
        -3 * sp.pi / 2
        + 1 / (1 + rho**2)
        + 3 * (rho / (1 + rho**2) + sp.atan(rho))
    )
    metric_a = sp.simplify(metric_b * radius**2)
    metric = sp.diag(
        conformal_factor * metric_a,
        -conformal_factor / metric_a,
        -conformal_factor * radius**2,
        -conformal_factor * radius**2 * sp.sin(theta) ** 2,
        -conformal_factor,
    )
    inverse = metric.inv()
    dimension = 5
    christoffel = [
        [
            [
                sp.simplify(
                    sp.Rational(1, 2)
                    * sum(
                        inverse[upper, lower]
                        * (
                            sp.diff(metric[lower, right], coordinates[left])
                            + sp.diff(metric[lower, left], coordinates[right])
                            - sp.diff(metric[left, right], coordinates[lower])
                        )
                        for lower in range(dimension)
                    )
                )
                for right in range(dimension)
            ]
            for left in range(dimension)
        ]
        for upper in range(dimension)
    ]
    riemann_mixed = [
        [
            [
                [
                    sp.simplify(
                        sp.diff(
                            christoffel[upper][lower][fourth],
                            coordinates[third],
                        )
                        - sp.diff(
                            christoffel[upper][lower][third],
                            coordinates[fourth],
                        )
                        + sum(
                            christoffel[upper][middle][third]
                            * christoffel[middle][lower][fourth]
                            - christoffel[upper][middle][fourth]
                            * christoffel[middle][lower][third]
                            for middle in range(dimension)
                        )
                    )
                    for fourth in range(dimension)
                ]
                for third in range(dimension)
            ]
            for lower in range(dimension)
        ]
        for upper in range(dimension)
    ]
    riemann_down = [
        [
            [
                [
                    sp.simplify(
                        metric[first, first]
                        * riemann_mixed[first][second][third][fourth]
                    )
                    for fourth in range(dimension)
                ]
                for third in range(dimension)
            ]
            for second in range(dimension)
        ]
        for first in range(dimension)
    ]
    return sp.lambdify(
        (rho, outside, theta),
        ([metric[index, index] for index in range(dimension)], riemann_down),
        "numpy",
    )


def _egb_equation(metric: np.ndarray, riemann: np.ndarray, alpha: float) -> np.ndarray:
    inverse = np.linalg.inv(metric)
    ricci = np.einsum("ac,abcd->bd", inverse, riemann)
    ricci_scalar = np.einsum("ab,ab", inverse, ricci)
    ricci_up = inverse @ ricci @ inverse
    ricci_square_tensor = ricci @ inverse @ ricci
    ricci_riemann = np.einsum("cd,acbd->ab", ricci_up, riemann)
    riemann_up_three = np.einsum(
        "cC,dD,eE,aCDE->acde", inverse, inverse, inverse, riemann
    )
    riemann_square_tensor = np.einsum(
        "acde,bcde->ab", riemann_up_three, riemann
    )
    ricci_squared = np.einsum("ab,aA,bB,AB", ricci, inverse, inverse, ricci)
    riemann_up_four = np.einsum(
        "aA,bB,cC,dD,ABCD->abcd",
        inverse,
        inverse,
        inverse,
        inverse,
        riemann,
    )
    riemann_squared = np.einsum("abcd,abcd", riemann, riemann_up_four)
    gauss_bonnet = (
        ricci_scalar**2 - 4 * ricci_squared + riemann_squared
    )
    lanczos = 2 * (
        ricci_scalar * ricci
        - 2 * ricci_square_tensor
        - 2 * ricci_riemann
        + riemann_square_tensor
    ) - 0.5 * metric * gauss_bonnet
    einstein = ricci - 0.5 * metric * ricci_scalar
    return einstein + alpha * lanczos


def _principal_riemann(covector: np.ndarray, perturbation: np.ndarray) -> np.ndarray:
    answer = np.zeros((5, 5, 5, 5), dtype=float)
    for first in range(5):
        for second in range(5):
            for third in range(5):
                for fourth in range(5):
                    answer[first, second, third, fourth] = 0.5 * (
                        covector[third]
                        * covector[second]
                        * perturbation[first, fourth]
                        + covector[fourth]
                        * covector[first]
                        * perturbation[second, third]
                        - covector[fourth]
                        * covector[second]
                        * perturbation[first, third]
                        - covector[third]
                        * covector[first]
                        * perturbation[second, fourth]
                    )
    return answer


def _principal_action(
    metric: np.ndarray,
    background_riemann: np.ndarray,
    covector: np.ndarray,
    perturbation: np.ndarray,
    alpha: float,
) -> np.ndarray:
    variation = _principal_riemann(covector, perturbation)
    # The EGB equation is quadratic in curvature, so this centered difference
    # extracts its linear principal variation exactly up to floating arithmetic.
    return 0.5 * (
        _egb_equation(metric, background_riemann + variation, alpha)
        - _egb_equation(metric, background_riemann - variation, alpha)
    )


def _kinetic_record(
    background,
    rho: float,
    outside: float,
    alpha: float,
) -> dict[str, object]:
    diagonal, riemann = background(rho, outside, np.pi / 2)
    metric = np.diag(np.asarray(diagonal, dtype=float).reshape(5))
    riemann = np.asarray(riemann, dtype=float)
    positive_directions = [
        index for index in range(5) if metric[index, index] > 0
    ]
    if len(positive_directions) != 1:
        raise RuntimeError("background does not have one local time direction")
    time_index = positive_directions[0]
    covector = np.zeros(5)
    covector[time_index] = 1.0
    spatial = [index for index in range(5) if index != time_index]
    pairs = [
        (left, right)
        for offset, left in enumerate(spatial)
        for right in spatial[offset:]
    ]
    matrix = np.zeros((10, 10))
    for column, (left, right) in enumerate(pairs):
        perturbation = np.zeros((5, 5))
        perturbation[left, right] = 1.0
        perturbation[right, left] = 1.0
        output = _principal_action(
            metric, riemann, covector, perturbation, alpha
        )
        for row, (first, second) in enumerate(pairs):
            matrix[row, column] = output[first, second]
    singular_values = np.linalg.svd(matrix, compute_uv=False)
    tolerance = float(np.max(singular_values) * 1.0e-10)
    rank = int(np.count_nonzero(singular_values > tolerance))

    gauge_residual = 0.0
    for direction in range(5):
        gauge_vector = np.zeros(5)
        gauge_vector[direction] = 1.0
        gauge_perturbation = np.outer(covector, gauge_vector) + np.outer(
            gauge_vector, covector
        )
        gauge_output = _principal_action(
            metric, riemann, covector, gauge_perturbation, alpha
        )
        gauge_residual = max(gauge_residual, float(np.max(np.abs(gauge_output))))

    return {
        "rho": rho,
        "outside_y": outside,
        "time_coordinate_index": time_index,
        "spatial_component_order": [f"{left}{right}" for left, right in pairs],
        "rank": rank,
        "dimension": 10,
        "singular_values": [float(value) for value in singular_values],
        "rank_tolerance": tolerance,
        "gauge_principal_residual": gauge_residual,
    }


def run() -> dict[str, object]:
    if OUTPUT.exists():
        raise RuntimeError("Gauss-Bonnet kinetic-rank output already exists")

    authenticated = []
    for path in (GEOMETRY, CRITICALITY):
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

    background = _background_function()
    rho_samples = (-2.0, -1.0, 0.0, 0.5, 1.0, 2.0)
    alpha = 1015.0 / 144.0
    gr_control = [
        _kinetic_record(background, rho, 0.0, 0.0) for rho in rho_samples
    ]
    local_critical = [
        _kinetic_record(background, rho, 0.0, alpha) for rho in rho_samples
    ]
    outside_diagnostic = [
        _kinetic_record(background, rho, 1.0, alpha)
        for rho in (-1.0, 0.5, 2.0)
    ]

    record = {
        "artifact_id": "NSC-1-S-ONE-GAUSS-BONNET-KINETIC-RANK",
        "schema": "NSC-1-S-ONE-GAUSS-BONNET-KINETIC-RANK-v1",
        "classification": "the_critical_EGB_bulk_loses_half_of_the_time_principal_rank_exactly_on_the_local_sheet_and_requires_induced_boundary_dynamics",
        "authenticated_inputs": authenticated,
        "method": {
            "background": "exact_symbolic_5D_Riemann_tensor",
            "principal_variation": "centered_exact_linear_part_of_the_quadratic_EGB_curvature_equation",
            "time_covector": "coordinate_with_positive_metric_sign_at_each_point",
            "physical_acceleration_space": "ten_symmetric_components_tangent_to_the_non_null_time_slice",
            "gauge_control": "five_xi_(a_X_b)_directions",
            "rank_tolerance": "largest_singular_value_times_1e-10",
        },
        "GR_control": gr_control,
        "critical_EGB_local_sheet": local_critical,
        "critical_EGB_outside_diagnostic": outside_diagnostic,
        "result": {
            "GR_rank_all_samples": sorted({item["rank"] for item in gr_control}),
            "critical_local_rank_all_samples": sorted(
                {item["rank"] for item in local_critical}
            ),
            "critical_outside_rank_all_samples": sorted(
                {item["rank"] for item in outside_diagnostic}
            ),
            "maximum_gauge_residual": max(
                item["gauge_principal_residual"]
                for item in gr_control + local_critical + outside_diagnostic
            ),
            "local_sheet_is_kinetically_critical": True,
        },
        "causal_interpretation": {
            "pure_bulk_EGB_completion_healthy": False,
            "failure_owner": "critical_bulk_principal_rank_on_y_equals_zero",
            "not_gauge": "all_tested_diffeomorphism_principal_directions_vanish",
            "not_global_rank_loss": "the_same_EGB_background_has_full_tested_rank_at_y_equals_one",
            "required_repair": "include_the_induced_local_spectral_Einstein_term_and_off_diagonal_sheet_modes_in_the_interface_principal_system",
            "theory_scope": "this_rejects_pure_critical_bulk_EGB_as_the_complete_propagation_law_not_the_one_operator_architecture",
        },
        "next_result": "add_the_fixed_local_spectral_Einstein_kinetic_block_to_the_critical_bulk_interface_and_test_whether_the_combined_principal_matrix_restores_full_local_rank",
        "nonclaims": {
            "complete_characteristic_polynomial_evaluated": False,
            "combined_bulk_boundary_system_tested": False,
            "full_hyperbolicity_proved": False,
        },
        "gate": {
            "inputs_authenticated": len(authenticated) == 2,
            "GR_control_full_rank": all(item["rank"] == 10 for item in gr_control),
            "critical_local_rank_is_five": all(
                item["rank"] == 5 for item in local_critical
            ),
            "gauge_residual_zero": all(
                item["gauge_principal_residual"] == 0.0
                for item in gr_control + local_critical + outside_diagnostic
            ),
            "outside_diagnostic_full_rank": all(
                item["rank"] == 10 for item in outside_diagnostic
            ),
            "combined_interface_rank_restored": False,
        },
        "terminal": True,
    }
    for key, value in record["gate"].items():
        if key != "combined_interface_rank_restored" and value is not True:
            raise RuntimeError(f"Gauss-Bonnet kinetic-rank gate failed: {key}")
    OUTPUT.write_bytes(canonical_json_bytes(record))
    return record


def main() -> int:
    record = run()
    print(json.dumps(record, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
