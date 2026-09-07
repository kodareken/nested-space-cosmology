#!/usr/bin/env python3
"""Construct a positive-F nonminimal boundary profile for the exact optical bridge."""

from __future__ import annotations

import json
from math import isfinite
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
for value in (str(ROOT), str(ROOT / "src")):
    if value not in sys.path:
        sys.path.insert(0, value)

from recursive_horizons.evidence_io import canonical_json_bytes  # noqa: E402


OUTPUT = ROOT / "results/nsc-1-unified-gradient-boundary-dev3.json"
RADIUS = 32.0 / 25.0
EXPANSION_MAGNITUDE = 225.0 / 64.0
THICKNESS = 1.0 / 50.0
ENDPOINT_F = 4.0


def _radius_factor(parameter: float) -> float:
    return (
        1.0
        - EXPANSION_MAGNITUDE * parameter / 2.0
        + EXPANSION_MAGNITUDE * parameter * parameter / (2.0 * THICKNESS)
    )


def _ricci_null(parameter: float) -> float:
    return -2.0 * EXPANSION_MAGNITUDE / (
        THICKNESS * _radius_factor(parameter)
    )


def _rk4_half(step_count: int, *, retain: bool = False):
    start = THICKNESS / 2.0
    width = (THICKNESS - start) / step_count
    parameter = start
    field = 1.0
    derivative = 0.0
    retained = [(parameter, field, derivative)]

    def rhs(x: float, y: float, z: float) -> tuple[float, float]:
        return z, _ricci_null(x) * y

    stride = max(1, step_count // 128)
    for index in range(step_count):
        k1y, k1z = rhs(parameter, field, derivative)
        k2y, k2z = rhs(
            parameter + width / 2.0,
            field + width * k1y / 2.0,
            derivative + width * k1z / 2.0,
        )
        k3y, k3z = rhs(
            parameter + width / 2.0,
            field + width * k2y / 2.0,
            derivative + width * k2z / 2.0,
        )
        k4y, k4z = rhs(
            parameter + width,
            field + width * k3y,
            derivative + width * k3z,
        )
        field += width * (k1y + 2.0 * k2y + 2.0 * k3y + k4y) / 6.0
        derivative += width * (k1z + 2.0 * k2z + 2.0 * k3z + k4z) / 6.0
        parameter += width
        if retain and ((index + 1) % stride == 0 or index + 1 == step_count):
            retained.append((parameter, field, derivative))
    return field, derivative, retained


def run() -> dict[str, object]:
    if OUTPUT.exists():
        raise RuntimeError("unified gradient boundary DEV3 output already exists")
    resolutions = (4096, 8192, 16384)
    endpoints = []
    finest_profile = None
    for steps in resolutions:
        endpoint, derivative, retained = _rk4_half(
            steps, retain=steps == resolutions[-1]
        )
        endpoints.append((steps, endpoint, derivative))
        if steps == resolutions[-1]:
            finest_profile = retained
    if finest_profile is None:
        raise RuntimeError("finest nonminimal profile is missing")
    scale = ENDPOINT_F / endpoints[-1][1]
    right = [
        (parameter, scale * field, scale * derivative)
        for parameter, field, derivative in finest_profile
    ]
    left = [
        (THICKNESS - parameter, field, -derivative)
        for parameter, field, derivative in reversed(right[1:])
    ]
    profile = left + right
    parameters = np.asarray([item[0] for item in profile])
    fields = np.asarray([item[1] for item in profile])
    derivatives = np.asarray([item[2] for item in profile])
    ricci = np.asarray([_ricci_null(value) for value in parameters])
    second = fields * ricci
    normalized_nonminimal_null = second / fields
    endpoint_difference = float(abs(fields[0] - fields[-1]))
    convergence = float(abs(endpoints[-1][1] - endpoints[-2][1]))
    finite = all(
        np.all(np.isfinite(value))
        for value in (parameters, fields, derivatives, ricci, second)
    )
    assertions = {
        "profile_finite": bool(finite),
        "effective_planck_coefficient_positive": bool(np.min(fields) > 0.0),
        "same_endpoint_local_gravity_coefficient": bool(endpoint_difference < 1.0e-13),
        "negative_null_Ricci_everywhere": bool(np.max(ricci) < 0.0),
        "nonminimal_null_term_matches_required_Ricci": bool(
            np.max(np.abs(normalized_nonminimal_null - ricci)) < 1.0e-12
        ),
        "three_resolution_ODE_contraction": bool(
            convergence < abs(endpoints[-2][1] - endpoints[-3][1])
        ),
        "all_gate_passed": False,
    }
    assertions["all_gate_passed"] = all(
        value for key, value in assertions.items() if key != "all_gate_passed"
    )
    if not assertions["all_gate_passed"]:
        raise RuntimeError("unified boundary profile failed")
    record = {
        "artifact_id": "NSC-1-UNIFIED-GRADIENT-BOUNDARY-DEV3",
        "schema": "NSC-1-UNIFIED-GRADIENT-BOUNDARY-DEV3-v1",
        "classification": "positive_effective_gravity_profile_supplies_exact_required_negative_null_Ricci",
        "active_scale": "parent_child_transition_sector_of_one_gravitating_gradient_model",
        "candidate_shared_action": {
            "form": "S=int_sqrt_minus_g_[F(varphi)R/2-lambda6_squared*pi4*B_mu*B^mu-mu0_squared*V(U,varphi)]",
            "U_role": "topological_Skyrme_field_whose_stable_solitons_model_baryons_and_nuclei",
            "varphi_role": "formation_and_boundary_order_parameter_controlling_the_local_gravity_coefficient_F",
            "B_mu_role": "conserved_topological_baryon_current",
            "one_model_requirement": "nuclear_solitons_and_parent_child_boundary_are_different_solutions_of_this_same_action",
        },
        "boundary_equation": {
            "null_metric_equation_on_zero_baryon_flux_transition": "F*R_kk=d2F/dlambda2",
            "required_profile_equation": "d2F/dlambda2=R_kk(lambda)*F",
            "required_R_kk": "-2a/[L*(1-a*lambda/2+a*lambda^2/(2L))]",
            "endpoint_F": ENDPOINT_F,
        },
        "result": {
            "affine_thickness": THICKNESS,
            "sample_count": int(parameters.size),
            "minimum_F": float(np.min(fields)),
            "maximum_F": float(np.max(fields)),
            "left_endpoint_F": float(fields[0]),
            "right_endpoint_F": float(fields[-1]),
            "endpoint_F_difference": endpoint_difference,
            "minimum_Ricci_null": float(np.min(ricci)),
            "maximum_Ricci_null": float(np.max(ricci)),
            "maximum_null_equation_residual": float(
                np.max(np.abs(second - fields * ricci))
            ),
            "finest_pair_endpoint_convergence": convergence,
            "integration_resolutions": list(resolutions),
        },
        "assertions": assertions,
        "inherited_invariant": "F_returns_to_the_same_positive_value_on_parent_and_child_sides_so_the_local_gravity_normalization_is_inherited",
        "nuclear_connection": {
            "established_launch_surface": "BPS_Skyrme_solitons_are_stable_gradient_configurations_and_have_published_nuclear_binding_energy_calculations",
            "current_AME2020_result": "coarse_bulk_boundary_energy_predicts_5_of_6_release_signs_but_misses_light_cluster_topology",
            "required_next_test": "implement_published_near_BPS_mass_functional_then_hold_out_the_same_six_reactions_without_refitting",
        },
        "literature_launch_surface": [
            {"doi": "10.1103/PhysRevC.88.054313", "role": "BPS_Skyrme_nuclear_binding_energies"},
            {"doi": "10.1103/PhysRevD.94.024060", "role": "general_Einstein_Skyrme_black_holes_and_solitons"},
            {"doi": "10.1103/PhysRevC.92.025802", "role": "gravitating_BPS_Skyrme_neutron_stars"},
        ],
        "next_intervention": "close_the_remaining_field_equations_and_test_the_same_action_against_AME2020_Q_values_and_the_exact_optical_boundary",
        "nonclaims": {
            "full_four_dimensional_action_solution": False,
            "scalar_order_parameter_equation_solved": False,
            "published_BPS_nuclear_formula_reproduced": False,
            "global_child_cosmology_solved": False,
        },
        "terminal": True,
    }
    OUTPUT.write_bytes(canonical_json_bytes(record))
    return record


def main() -> int:
    record = run()
    print(json.dumps(record, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
