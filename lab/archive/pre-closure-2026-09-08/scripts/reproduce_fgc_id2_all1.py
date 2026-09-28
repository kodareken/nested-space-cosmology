#!/usr/bin/env python3
"""Regenerate the pre-calibration FGC-1-ID2-ALL1 static ledger."""

from __future__ import annotations

import argparse
from dataclasses import asdict, is_dataclass
from fractions import Fraction
from hashlib import sha256
import json
from math import isfinite, log
from pathlib import Path
import sys
import tomllib
from typing import Any, Mapping, Sequence

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-id2-all1.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-id2-all1.json"
OWNER_DOCUMENT = REPOSITORY / "docs/fgc-id2-all1.md"

from recursive_horizons.fgc.evolution.gr0_calibration import (  # noqa: E402
    GR0GridInitialData,
    construct_gr0_grid_initial_data,
)
from recursive_horizons.fgc.evolution.health_monitor import (  # noqa: E402
    compact_vacuum_buffer_window,
)
from recursive_horizons.fgc.evolution.proto4_admission import (  # noqa: E402
    PROTO4_SPECTRAL_FIELD_ORDER,
    nested_spectral_admission,
    proper_radial_profile,
    spectral_field_budgets,
)
from recursive_horizons.fgc.evolution.static_initial_admission import (  # noqa: E402
    FGCQRGridInitialData,
    construct_fgcqr_grid_initial_data,
)
from recursive_horizons.fgc.initial_data_family import (  # noqa: E402
    FGCQRActionParameters,
    InitialDataParameters,
    InitialDataSolution,
)
from recursive_horizons.fgc.initial_data_preflight import (  # noqa: E402
    PulseParameters,
    gr0_constraint_residuals,
    gr0_constraint_rhs,
    pulse_fields,
    solve_gr0_initial_slice,
)
from recursive_horizons.fgc.scoped_run_authorization import (  # noqa: E402
    SF1_PROTOCOL_ARTIFACT_ID,
    SF1_PROTOCOL_V4_ARTIFACT_ID,
    validate_sf1_protocol,
)


Q = Fraction
ARTIFACT_ID = "FGC-1-ID2-ALL1"
PROJECT_VERSION = "0.11.0"
REQUIRED_GATE = "all_amplitude_all_case_static_initial_admission_completed"
EXPECTED_PATHS = {
    "protocol_config": "configs/fgc/fgc-2-sf1-protocol-v4.toml",
    "physical_protocol_config": "configs/fgc/fgc-2-sf1-protocol-v3.toml",
    "protocol_freeze_result": "results/fgc-1-pro4-frz1.json",
    "successor_monitor_result": "results/fgc-1-hlt2-mon2.json",
    "initial_family_result": "results/fgc-1-id1-fam1.json",
    "physical_constraints_result": "results/fgc-1-con4-phy1.json",
    "regular_center_result": "results/fgc-1-ctr1-reg1.json",
    "run_domain_result": "results/fgc-1-dom4-run1.json",
    "FGCQR_health_result": "results/fgc-1-hyp2-md1.json",
}
RESULT_REQUIREMENTS = {
    "protocol_freeze_result": (
        "FGC-1-PRO4-FRZ1",
        "PROTO4_outcome_neutral_protocol_frozen",
        None,
    ),
    "successor_monitor_result": (
        "FGC-1-HLT2-MON2",
        "PROTO4_successor_admission_monitor_implemented",
        None,
    ),
    "initial_family_result": (
        "FGC-1-ID1-FAM1",
        "nonzero_width_finite_mass_constraint_compatible_family_constructed",
        "FGC-QR",
    ),
    "physical_constraints_result": (
        "FGC-1-CON4-PHY1",
        "physical_gauge_reduction_constraint_system_closed",
        "FGC-QR",
    ),
    "regular_center_result": (
        "FGC-1-CTR1-REG1",
        "regular_center_formulation_verified",
        "FGC-QR",
    ),
    "run_domain_result": (
        "FGC-1-DOM4-RUN1",
        "nonzero_classical_spherical_run_envelope_passed",
        "FGC-QR",
    ),
    "FGCQR_health_result": (
        "FGC-1-HYP2-MD1",
        "quantitative_all_covector_weak_coupling_health_envelope_passed",
        "FGC-QR",
    ),
}
EXPECTED_SCOPE = {
    "target_protocol": "FGC-2-SF1-PROTO4",
    "role": "pre_calibration_all_amplitude_all_case_static_initial_admission",
    "calibration_branch": "GR-0",
    "holdout_branch": "FGC-QR",
    "amplitude_count": 7,
    "held_out_case_count": 5,
    "FGCQR_cartesian_product_count": 35,
    "all_amplitudes_must_be_adjudicated_before_first_dynamic_calibration": True,
    "amplitude_dynamic_eligibility_requires_all_five_corresponding_FGCQR_static_slices": True,
    "static_ineligibility_is_not_a_candidate_mechanism_outcome": True,
    "calibration_trajectory_read": False,
    "calibration_output_namespace_read_or_written": False,
    "holdout_output_namespace_read_or_written": False,
    "FGCQR_evolution_outcome_read": False,
    "collapse_activation_or_defocusing_classified": False,
}
EXPECTED_PHYSICAL_INPUTS = {
    "length_unit_L0": "4",
    "cutoff_Lambda": "16",
    "planck_mass": "2",
    "beta": "-1/4",
    "scalar_mass": "3",
    "quartic_coupling": "1/2",
    "eta": "1/2",
    "chi_center": "12",
    "chi_base_half_width": "2",
    "phi_half_width": "2",
    "phi_seed_amplitude": "1/131072",
    "outer_radius": "128",
    "measurement_radius_maximum": "24",
    "amplitude_candidates": ["2", "5/2", "3", "7/2", "4", "9/2", "5"],
}
EXPECTED_CASES = (
    ("FGCQR-CENTRAL", "1", "1", "5/2"),
    ("FGCQR-SEED-HALF", "1", "1/2", "5/2"),
    ("FGCQR-SEED-DOUBLE", "1", "2", "5/2"),
    ("FGCQR-WIDTH-SEVEN-EIGHTHS", "7/8", "1", "41/16"),
    ("FGCQR-WIDTH-NINE-EIGHTHS", "9/8", "1", "39/16"),
)
EXPECTED_NUMERICS = {
    "resolutions": [1025, 2049, 4097],
    "primary_constraint_method": "RK4",
    "primary_diagnostic_spatial_order": 4,
    "comparator_constraint_method": "SSPRK3",
    "comparator_diagnostic_spatial_order": 2,
    "proper_spectral_window_taper_fraction": "1/8",
    "state_hash_dtype": "little_endian_binary64",
    "development_or_prior_calibration_runs_may_not_enter_evidence": True,
}
EXPECTED_THRESHOLDS = {
    "initial_compactness_minimum": "1/10",
    "initial_compactness_maximum": "3/4",
    "maximum_physical_constraint_residual_infinity": "1/10000000000",
    "minimum_RK4_common_grid_profile_order": "3",
    "minimum_SSPRK3_common_grid_profile_order": "5/2",
    "maximum_fine_cross_method_relative_outer_mass_difference": "1/1024",
    "maximum_fine_cross_method_peak_compactness_difference": "1/1024",
    "minimum_effective_planck_ratio": "1/2",
    "minimum_reduction_constraint_decay_order_diagnostic": "1",
    "strict_compactness_bounds": False,
    "positive_vacuum_metric_denominator_required": True,
    "finite_positive_constraint_jacobian_required": True,
    "GR0_nested_weighted_spectral_admission_required": True,
    "FGCQR_individual_weighted_spectral_budgets_required": True,
    "FGCQR_nested_tail_result_is_diagnostic_for_this_STATIC_PROTO3_ledger": True,
}
EXPECTED_PROOF_KEYS = {
    "PROTO4_and_exact_PROTO3_physical_predecessor_must_validate",
    "PROTO4_freeze_and_HLT2_must_be_canonical_and_passed",
    "ID1_CON4_CTR1_DOM4_and_HYP2_scope_chain_must_be_canonical",
    "seven_GR0_candidates_must_be_evaluated_on_both_methods_and_three_grids",
    "thirty_five_FGCQR_slices_must_be_evaluated_on_both_methods_and_three_grids",
    "every_input_configuration_and_complete_grid_state_must_be_hash_bound",
    "physical_constraints_must_be_solved_on_the_complete_initial_grid",
    "regular_center_exact_inner_and_outer_vacuum_finite_mass_and_no_trapping_must_be_checked",
    "compactness_window_must_be_applied_to_each_expanded_FGCQR_slice",
    "profile_convergence_and_cross_method_agreement_must_be_checked",
    "PROTO4_weighted_spectral_admission_must_decide_GR0_static_eligibility",
    "FGCQR_nested_tail_diagnostic_must_not_rewrite_PROTO3_five_boolean_static_contract",
    "ineligible_amplitudes_must_retain_exact_predeclared_reasons",
    "no_dynamic_output_namespace_may_be_read_or_written",
    "canonical_hash_bound_result_required",
}
EXPECTED_CLAIMS = {
    "all_case_static_ledger_completed": True,
    "fresh_GR0_dynamic_calibration_authorized": False,
    "fresh_GR0_dynamic_calibration_completed": False,
    "PROTO4_resolved_holdout_manifest_authorized": False,
    "SGBL_comparison_execution_authorized": False,
    "FGCQR_holdout_execution_authorized": False,
    "classical_spherical_diagnostic_authorized": False,
    "retained_EFT_evolution_authorized": False,
    "physical_transition_claim_authorized": False,
    "FGCQR_mechanism_rejected": False,
    "general_gradient_route_rejected": False,
    "singularity_resolution_derived": False,
    "child_domain_or_topology_derived": False,
    "dark_sector_mechanism_derived": False,
    "varying_locally_measured_c_derived": False,
}
IMPLEMENTATION = tuple(
    REPOSITORY / name
    for name in (
        "src/recursive_horizons/fgc/evolution/static_initial_admission.py",
        "src/recursive_horizons/fgc/evolution/__init__.py",
        "src/recursive_horizons/fgc/evolution/gr0_calibration.py",
        "src/recursive_horizons/fgc/evolution/proto4_admission.py",
        "src/recursive_horizons/fgc/initial_data_family.py",
        "src/recursive_horizons/fgc/initial_data_preflight.py",
        "scripts/reproduce_fgc_id2_all1.py",
    )
)


def _pairs(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    answer: dict[str, Any] = {}
    for key, value in pairs:
        if key in answer:
            raise ValueError(f"duplicate JSON object key: {key}")
        answer[key] = value
    return answer


def _fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def _fraction(name: str, value: Any, *, positive: bool = False) -> Fraction:
    if not isinstance(value, str) or not value:
        raise ValueError(f"{name} must be a canonical rational string")
    try:
        answer = Q(value)
    except (ValueError, ZeroDivisionError) as exc:
        raise ValueError(f"{name} must be a canonical rational string") from exc
    if _fraction_text(answer) != value:
        raise ValueError(f"{name} must use canonical rational encoding")
    if positive and answer <= 0:
        raise ValueError(f"{name} must be positive")
    return answer


def _serial(value: Any) -> Any:
    if isinstance(value, Fraction):
        return _fraction_text(value)
    if is_dataclass(value):
        return _serial(asdict(value))
    if isinstance(value, Mapping):
        return {str(key): _serial(item) for key, item in value.items()}
    if isinstance(value, np.ndarray):
        return _serial(value.tolist())
    if isinstance(value, np.generic):
        return _serial(value.item())
    if isinstance(value, (tuple, list)):
        return [_serial(item) for item in value]
    if isinstance(value, float):
        if not isfinite(value):
            raise ValueError("nonfinite value cannot enter canonical evidence")
        return value
    if value is None or isinstance(value, (str, bool, int)):
        return value
    raise TypeError(f"unsupported evidence value {type(value).__name__}")


def _canonical(value: Any, *, pretty: bool = True) -> str:
    kwargs = {"sort_keys": True, "allow_nan": False}
    if pretty:
        kwargs["indent"] = 2
    else:
        kwargs["separators"] = (",", ":")
    return json.dumps(_serial(value), **kwargs) + ("\n" if pretty else "")


def _digest(value: Any) -> str:
    return sha256(_canonical(value, pretty=False).encode("utf-8")).hexdigest()


def _sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _rel(path: Path) -> str:
    return path.resolve().relative_to(REPOSITORY).as_posix()


def _path(name: str, value: Any, expected: str) -> Path:
    if value != expected:
        raise ValueError(f"{name} must name {expected}")
    candidate = (REPOSITORY / expected).resolve()
    if not candidate.is_file() or _rel(candidate) != expected:
        raise ValueError(f"{name} must be a canonical existing repository path")
    return candidate


def _load_json(path: Path, name: str, *, canonical: bool = True) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"), object_pairs_hook=_pairs)
    except (OSError, json.JSONDecodeError, ValueError) as exc:
        raise ValueError(f"cannot load {name}: {exc}") from exc
    if not isinstance(value, dict):
        raise ValueError(f"{name} must be a JSON object")
    if canonical and path.read_text(encoding="utf-8") != _canonical(value):
        raise ValueError(f"{name} is not canonically serialized")
    return value


def _check_result(path: Path, artifact_id: str, gate: str, branch: str | None) -> dict[str, Any]:
    result = _load_json(path, artifact_id)
    if result.get("artifact_id") != artifact_id:
        raise ValueError(f"{artifact_id} identity differs")
    if result.get("gate_status", {}).get(gate) is not True:
        raise ValueError(f"{artifact_id} required gate is not true")
    if branch is not None and result.get("scope_bindings", {}).get("target_branch") != branch:
        raise ValueError(f"{artifact_id} belongs to another branch")
    return result


def _state_sha256(branch: str, point_count: int, state: Any) -> str:
    digest = sha256()
    digest.update(branch.encode("ascii"))
    digest.update(point_count.to_bytes(8, "little", signed=False))
    for name, value in (("u", state.u), ("p", state.p), ("q", state.q)):
        array = np.asarray(value, dtype="<f8", order="C")
        digest.update(name.encode("ascii"))
        digest.update(len(array.shape).to_bytes(2, "little"))
        for dimension in array.shape:
            digest.update(int(dimension).to_bytes(8, "little", signed=False))
        digest.update(array.tobytes(order="C"))
    return digest.hexdigest()


def _profile_error(coarse: Any, fine: Any) -> float:
    if fine.step_count % coarse.step_count != 0:
        raise ValueError("static profile comparison requires nested support grids")
    stride = fine.step_count // coarse.step_count
    sampled = fine.points[::stride]
    if len(sampled) != len(coarse.points):
        raise ValueError("static nested profile lengths differ")
    return max(
        max(
            abs(left.radial_metric - right.radial_metric),
            abs(left.angular_extrinsic_curvature - right.angular_extrinsic_curvature),
        )
        for left, right in zip(coarse.points, sampled, strict=True)
    )


def _profile_convergence(solutions: Sequence[Any]) -> dict[str, Any]:
    if len(solutions) != 3:
        raise ValueError("static profile convergence requires exactly three solutions")
    coarse_medium = _profile_error(solutions[0], solutions[1])
    medium_fine = _profile_error(solutions[1], solutions[2])
    if coarse_medium <= 0.0 or medium_fine <= 0.0:
        raise ValueError("static profile convergence is unresolved")
    refinement = solutions[1].step_count / solutions[0].step_count
    if refinement <= 1.0 or solutions[2].step_count / solutions[1].step_count != refinement:
        raise ValueError("static profile grids do not share one refinement ratio")
    order = log(coarse_medium / medium_fine) / log(refinement)
    return {
        "step_counts": [item.step_count for item in solutions],
        "coarse_medium_common_grid_error_infinity": coarse_medium,
        "medium_fine_common_grid_error_infinity": medium_fine,
        "observed_common_grid_profile_order": order,
    }


def _reduction_decay(values: Sequence[float]) -> dict[str, Any]:
    if len(values) != 3 or any(item <= 0.0 for item in values):
        raise ValueError("reduction decay requires three positive diagnostics")
    orders = tuple(log(left / right) / log(2.0) for left, right in zip(values, values[1:]))
    return {
        "values": list(values),
        "adjacent_observed_orders": orders,
        "minimum_observed_order": min(orders),
    }


def _spectrum(grid_data: Any, cutoff: float, measurement: float, taper: float) -> tuple[dict[str, Any], Mapping[str, Any]]:
    profile = proper_radial_profile(
        grid_data.state.u,
        grid_data.state.q,
        grid_data.grid.coordinates,
        cutoff=cutoff,
        measurement_radius_maximum=measurement,
    )
    proper = profile.proper_coordinates
    window = compact_vacuum_buffer_window(
        proper,
        support_minimum=float(proper[0]),
        support_maximum=float(proper[-1]),
        taper_width=float((proper[-1] - proper[0]) * taper),
    )
    budgets = spectral_field_budgets(
        profile.deviations,
        proper,
        cutoff=cutoff,
        window=window,
    )
    summary = {
        "proper_measurement_radius": float(proper[-1]),
        "maximum_round_trip_interpolation_infinity": float(
            np.max(profile.round_trip_interpolation_infinity, initial=0.0)
        ),
        "maximum_top_band_field_power_fraction": max(
            item.top_band_field_power_fraction for item in budgets.values()
        ),
        "maximum_top_band_derivative_weighted_power_fraction": max(
            item.top_band_derivative_weighted_power_fraction for item in budgets.values()
        ),
        "maximum_RMS_scale_over_Lambda": max(
            item.rms_scale_over_cutoff for item in budgets.values()
        ),
        "last_bin_diagnostic_true_field_count": sum(
            item.last_bin_alias_diagnostic for item in budgets.values()
        ),
        "every_individual_budget_passed": all(
            item.individual_admission_passed for item in budgets.values()
        ),
    }
    return summary, budgets


def _gr0_physical_constraint_residual(initial: GR0GridInitialData) -> float:
    maximum = 0.0
    radii = initial.grid.coordinates
    for index in range(1, initial.grid.point_count):
        radius = radii[index]
        radial_metric = initial.state.u[index, 2]
        angular_k = initial.state.p[index, 2] / (2.0 * radial_metric)
        radial_metric_r = initial.state.q[index, 2]
        _expected_lambda_r, angular_k_r = gr0_constraint_rhs(
            radius,
            (radial_metric, angular_k),
            initial.parameters,
        )
        fields = pulse_fields(radius, initial.parameters)
        residual = gr0_constraint_residuals(
            radius=radius,
            radial_metric=radial_metric,
            angular_extrinsic_curvature=angular_k,
            radial_metric_derivative=radial_metric_r,
            angular_extrinsic_curvature_derivative=angular_k_r,
            phi=fields["phi"],
            phi_r=fields["phi_r"],
            phi_pi=fields["phi_pi"],
            chi=fields["chi"],
            chi_r=fields["chi_r"],
            chi_pi=fields["chi_pi"],
            planck_mass=initial.parameters.planck_mass,
            scalar_mass=initial.parameters.scalar_mass,
            quartic_coupling=initial.parameters.quartic_coupling,
        )
        maximum = max(maximum, *(abs(float(value)) for value in residual))
    return maximum


def _method_record(
    branch: str,
    method: str,
    spatial_order: int,
    grids: Sequence[Any],
    solutions: Sequence[Any],
    *,
    cutoff: float,
    measurement: float,
    taper: float,
    physical_constraint_residuals: Sequence[float],
) -> dict[str, Any]:
    spectra = []
    budgets = []
    resolutions = []
    for grid_data, residual in zip(grids, physical_constraint_residuals, strict=True):
        spectral, field_budgets = _spectrum(grid_data, cutoff, measurement, taper)
        spectra.append(spectral)
        budgets.append(field_budgets)
        resolutions.append(
            {
                "point_count": grid_data.grid.point_count,
                "constraint_step_count": (
                    grid_data.constraint_solution.step_count
                    if isinstance(grid_data, FGCQRGridInitialData)
                    else grid_data.support_maximum_index - grid_data.support_minimum_index
                ),
                "state_sha256": _state_sha256(
                    branch, grid_data.grid.point_count, grid_data.state
                ),
                "outer_mass": grid_data.outer_mass,
                "peak_compactness": grid_data.peak_compactness,
                "peak_radius": grid_data.peak_radius,
                "physical_constraint_residual_infinity": residual,
                "reduction_constraint_infinity": grid_data.reduction_constraint_infinity,
                "no_initial_trapped_sphere": grid_data.no_initial_trapped_sphere,
                "spectral_summary": spectral,
            }
        )
    nested = nested_spectral_admission(
        tuple(item.grid.point_count for item in grids), budgets
    )
    return {
        "method": method,
        "diagnostic_spatial_order": spatial_order,
        "resolutions": resolutions,
        "profile_convergence": _profile_convergence(solutions),
        "reduction_constraint_decay_diagnostic": _reduction_decay(
            [item.reduction_constraint_infinity for item in grids]
        ),
        "nested_weighted_spectral_admission": nested,
        "every_individual_weighted_spectral_budget_passed": all(
            item["every_individual_budget_passed"] for item in spectra
        ),
    }


def _gr0_candidate(
    amplitude_text: str,
    configuration: Mapping[str, Any],
) -> tuple[dict[str, Any], dict[str, Any]]:
    amplitude = float(Q(amplitude_text))
    physical = configuration["physical"]
    numerics = configuration["numerics"]
    thresholds = configuration["thresholds"]
    parameters = PulseParameters(
        chi_amplitude=amplitude,
        center=physical["chi_center"],
        half_width=physical["chi_base_half_width"],
        phi_amplitude=physical["phi_seed_amplitude"],
        planck_mass=physical["planck_mass"],
        scalar_mass=physical["scalar_mass"],
        quartic_coupling=physical["quartic_coupling"],
    )
    config_payload = {
        "protocol_sha256": configuration["protocol_sha256"],
        "branch": "GR-0",
        "amplitude": amplitude_text,
        "parameters": asdict(parameters),
        "resolutions": numerics["resolutions"],
        "methods": [item[0] for item in numerics["methods"]],
    }
    methods = []
    raw_grids: dict[str, list[Any]] = {}
    for method, spatial_order in numerics["methods"]:
        grids = [
            construct_gr0_grid_initial_data(
                parameters,
                point_count=count,
                outer_radius=physical["outer_radius"],
                constraint_method=method,
                diagnostic_spatial_order=spatial_order,
            )
            for count in numerics["resolutions"]
        ]
        solutions = [
            solve_gr0_initial_slice(
                parameters,
                step_count=item.support_maximum_index - item.support_minimum_index,
                method=method,
            )
            for item in grids
        ]
        residuals = [_gr0_physical_constraint_residual(item) for item in grids]
        methods.append(
            _method_record(
                "GR-0",
                method,
                spatial_order,
                grids,
                solutions,
                cutoff=physical["cutoff"],
                measurement=physical["measurement"],
                taper=numerics["taper"],
                physical_constraint_residuals=residuals,
            )
        )
        raw_grids[method] = grids

    primary = raw_grids["RK4"][-1]
    comparator = raw_grids["SSPRK3"][-1]
    relative_mass = abs(primary.outer_mass - comparator.outer_mass) / max(
        abs(primary.outer_mass), abs(comparator.outer_mass), 1.0e-30
    )
    compactness_difference = abs(
        primary.peak_compactness - comparator.peak_compactness
    )
    all_resolutions = [item for records in raw_grids.values() for item in records]
    constraints_passed = all(
        record["physical_constraint_residual_infinity"]
        <= thresholds["maximum_constraint_residual"]
        for method in methods
        for record in method["resolutions"]
    )
    compactness_passed = all(
        thresholds["compactness_minimum"] <= item.peak_compactness
        <= thresholds["compactness_maximum"]
        for item in all_resolutions
    )
    convergence_passed = (
        methods[0]["profile_convergence"]["observed_common_grid_profile_order"]
        >= thresholds["minimum_RK4_order"]
        and methods[1]["profile_convergence"]["observed_common_grid_profile_order"]
        >= thresholds["minimum_SSPRK3_order"]
    )
    spectrum_passed = all(
        item["nested_weighted_spectral_admission"].admission_passed
        for item in methods
    )
    cross_method_passed = (
        relative_mass <= thresholds["maximum_relative_mass_difference"]
        and compactness_difference <= thresholds["maximum_compactness_difference"]
    )
    static_passed = all(
        (
            constraints_passed,
            compactness_passed,
            convergence_passed,
            spectrum_passed,
            cross_method_passed,
            all(item.no_initial_trapped_sphere for item in all_resolutions),
            all(item.outer_mass > 0.0 and isfinite(item.outer_mass) for item in all_resolutions),
        )
    )
    reasons = []
    if not constraints_passed:
        reasons.append("physical_constraint_residual_limit")
    if not compactness_passed:
        reasons.append("central_GR0_compactness_outside_protocol_window")
    if not convergence_passed:
        reasons.append("initial_constraint_profile_nonconvergence")
    if not spectrum_passed:
        reasons.append("PROTO4_nested_weighted_spectral_admission")
    if not cross_method_passed:
        reasons.append("initial_constraint_cross_method_disagreement")
    if not reasons and not static_passed:
        reasons.append("other_static_initial_premise")
    record = {
        "amplitude": amplitude_text,
        "config_sha256": _digest(config_payload),
        "methods": methods,
        "fine_cross_method_relative_outer_mass_difference": relative_mass,
        "fine_cross_method_peak_compactness_difference": compactness_difference,
        "constraint_compatible": constraints_passed,
        "finite_mass": all(item.outer_mass > 0.0 for item in all_resolutions),
        "initial_compactness_in_protocol_window": compactness_passed,
        "no_initial_trapped_sphere": all(
            item.no_initial_trapped_sphere for item in all_resolutions
        ),
        "profile_convergence_passed": convergence_passed,
        "cross_method_agreement_passed": cross_method_passed,
        "PROTO4_nested_weighted_spectral_admission_passed": spectrum_passed,
        "GR0_static_input_passed": static_passed,
        "GR0_static_ineligibility_reasons": reasons,
        "dynamic_trajectory_read": False,
    }
    record["initial_data_result_sha256"] = _digest(record)
    return record, raw_grids


def _fgcqr_case(
    amplitude_text: str,
    case: Mapping[str, Any],
    configuration: Mapping[str, Any],
) -> dict[str, Any]:
    physical = configuration["physical"]
    numerics = configuration["numerics"]
    thresholds = configuration["thresholds"]
    width_factor = Q(case["chi_width_factor"])
    seed_factor = Q(case["phi_seed_factor"])
    parameters = InitialDataParameters(
        chi_amplitude=float(Q(amplitude_text)),
        chi_half_width=physical["chi_base_half_width"] * float(width_factor),
        phi_amplitude=physical["phi_seed_amplitude"] * float(seed_factor),
        center=physical["chi_center"],
        phi_half_width=physical["phi_half_width"],
        outer_radius=physical["outer_radius"],
        action=physical["action"],
    )
    exact_chi_buffer = (
        Q(str(int(physical["chi_center"])))
        - Q(str(int(physical["chi_base_half_width"]))) * width_factor
    ) / Q(str(int(physical["length_unit"])))
    if _fraction_text(exact_chi_buffer) != case["exact_chi_center_buffer_over_L0"]:
        raise ValueError("FGC-QR exact chi centre buffer differs from the frozen case")
    config_payload = {
        "protocol_sha256": configuration["protocol_sha256"],
        "branch": "FGC-QR",
        "amplitude": amplitude_text,
        "case_id": case["case_id"],
        "chi_width_factor": case["chi_width_factor"],
        "phi_seed_factor": case["phi_seed_factor"],
        "parameters": asdict(parameters),
        "resolutions": numerics["resolutions"],
        "methods": [item[0] for item in numerics["methods"]],
    }
    methods = []
    raw_grids: dict[str, list[FGCQRGridInitialData]] = {}
    for method, spatial_order in numerics["methods"]:
        grids = [
            construct_fgcqr_grid_initial_data(
                parameters,
                point_count=count,
                constraint_method=method,
                diagnostic_spatial_order=spatial_order,
                maximum_constraint_residual=thresholds["maximum_constraint_residual"],
            )
            for count in numerics["resolutions"]
        ]
        methods.append(
            _method_record(
                "FGC-QR",
                method,
                spatial_order,
                grids,
                [item.constraint_solution for item in grids],
                cutoff=physical["cutoff"],
                measurement=physical["measurement"],
                taper=numerics["taper"],
                physical_constraint_residuals=[
                    item.physical_constraint_residual_infinity for item in grids
                ],
            )
        )
        raw_grids[method] = grids

    primary = raw_grids["RK4"][-1]
    comparator = raw_grids["SSPRK3"][-1]
    relative_mass = abs(primary.outer_mass - comparator.outer_mass) / max(
        abs(primary.outer_mass), abs(comparator.outer_mass), 1.0e-30
    )
    compactness_difference = abs(
        primary.peak_compactness - comparator.peak_compactness
    )
    all_resolutions = [item for records in raw_grids.values() for item in records]
    constraints_passed = all(
        item.physical_constraint_residual_infinity
        <= thresholds["maximum_constraint_residual"]
        for item in all_resolutions
    )
    regular_center = all(item.regular_center for item in all_resolutions)
    exact_buffers = all(
        item.exact_inner_vacuum_buffer and item.exact_outer_vacuum_buffer
        for item in all_resolutions
    )
    finite_mass = all(item.finite_mass and item.outer_mass > 0.0 for item in all_resolutions)
    compactness_values = [item.peak_compactness for item in all_resolutions]
    compactness_passed = all(
        thresholds["compactness_minimum"] <= value
        <= thresholds["compactness_maximum"]
        for value in compactness_values
    )
    no_trapping = all(item.no_initial_trapped_sphere for item in all_resolutions)
    convergence_passed = (
        methods[0]["profile_convergence"]["observed_common_grid_profile_order"]
        >= thresholds["minimum_RK4_order"]
        and methods[1]["profile_convergence"]["observed_common_grid_profile_order"]
        >= thresholds["minimum_SSPRK3_order"]
    )
    cross_method_passed = (
        relative_mass <= thresholds["maximum_relative_mass_difference"]
        and compactness_difference <= thresholds["maximum_compactness_difference"]
    )
    effective_planck_passed = all(
        item.minimum_effective_planck_coefficient
        / parameters.action.planck_mass**2
        >= thresholds["minimum_effective_planck_ratio"]
        for item in all_resolutions
    )
    vacuum_passed = all(
        item.minimum_vacuum_metric_denominator > 0.0 for item in all_resolutions
    )
    jacobian_passed = all(
        item.constraint_solution.minimum_abs_constraint_jacobian_determinant > 0.0
        and isfinite(item.constraint_solution.minimum_abs_constraint_jacobian_determinant)
        for item in all_resolutions
    )
    individual_spectra_passed = all(
        item["every_individual_weighted_spectral_budget_passed"] for item in methods
    )
    nested_spectral_diagnostic = all(
        item["nested_weighted_spectral_admission"].admission_passed for item in methods
    )
    static_premises = all(
        (
            constraints_passed,
            regular_center,
            exact_buffers,
            finite_mass,
            compactness_passed,
            no_trapping,
            convergence_passed,
            cross_method_passed,
            effective_planck_passed,
            vacuum_passed,
            jacobian_passed,
            individual_spectra_passed,
        )
    )
    reasons = []
    if not constraints_passed:
        reasons.append("physical_constraint_residual_limit")
    if not regular_center:
        reasons.append("regular_center")
    if not exact_buffers:
        reasons.append("exact_vacuum_buffer")
    if not finite_mass:
        reasons.append("finite_mass")
    if min(compactness_values) < thresholds["compactness_minimum"]:
        reasons.append("peak_compactness_below_protocol_minimum")
    if max(compactness_values) > thresholds["compactness_maximum"]:
        reasons.append("peak_compactness_above_protocol_maximum")
    if not no_trapping:
        reasons.append("initial_trapped_sphere")
    if not convergence_passed:
        reasons.append("initial_constraint_profile_nonconvergence")
    if not cross_method_passed:
        reasons.append("initial_constraint_cross_method_disagreement")
    if not effective_planck_passed:
        reasons.append("effective_planck_ratio")
    if not vacuum_passed:
        reasons.append("vacuum_metric_denominator")
    if not jacobian_passed:
        reasons.append("constraint_jacobian")
    if not individual_spectra_passed:
        reasons.append("individual_weighted_spectral_budget")
    record = {
        "amplitude": amplitude_text,
        "case_id": case["case_id"],
        "chi_width_factor": case["chi_width_factor"],
        "phi_seed_factor": case["phi_seed_factor"],
        "exact_chi_center_buffer_over_L0": case[
            "exact_chi_center_buffer_over_L0"
        ],
        "exact_combined_center_vacuum_buffer_over_L0": _fraction_text(
            Q(str(parameters.support_minimum)) / Q(str(physical["length_unit"]))
        ),
        "config_sha256": _digest(config_payload),
        "methods": methods,
        "minimum_observed_peak_compactness": min(compactness_values),
        "maximum_observed_peak_compactness": max(compactness_values),
        "fine_cross_method_relative_outer_mass_difference": relative_mass,
        "fine_cross_method_peak_compactness_difference": compactness_difference,
        "constraint_compatible": constraints_passed,
        "regular_center": regular_center,
        "exact_inner_and_outer_vacuum_buffers": exact_buffers,
        "finite_mass": finite_mass,
        "initial_compactness_in_protocol_window": compactness_passed,
        "no_initial_trapped_sphere": no_trapping,
        "profile_convergence_passed": convergence_passed,
        "cross_method_agreement_passed": cross_method_passed,
        "positive_effective_planck_and_vacuum_branches": (
            effective_planck_passed and vacuum_passed and jacobian_passed
        ),
        "every_individual_weighted_spectral_budget_passed": individual_spectra_passed,
        "nested_weighted_spectral_admission_diagnostic_passed": nested_spectral_diagnostic,
        "nested_tail_diagnostic_decides_this_static_premise": False,
        "static_premises_passed": static_premises,
        "static_ineligibility_reasons": reasons,
        "FGCQR_evolution_outcome_inspected": False,
    }
    record["initial_data_result_sha256"] = _digest(record)
    return record


def load_config(path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    path = path.resolve()
    if not path.is_file():
        raise ValueError("ID2 config does not exist")
    with path.open("rb") as handle:
        raw = tomllib.load(handle)
    expected_root = {
        "schema_version",
        "artifact_id",
        "project_version",
        "metric_signature",
        "riemann_convention",
        *EXPECTED_PATHS,
        "scope",
        "physical_inputs",
        "held_out_cases",
        "numerics",
        "thresholds",
        "proof_contract",
        "claims",
    }
    if set(raw) != expected_root:
        raise ValueError("ID2 root keys differ")
    if (
        raw["schema_version"] != 1
        or raw["artifact_id"] != ARTIFACT_ID
        or raw["project_version"] != PROJECT_VERSION
        or raw["metric_signature"] != "-+++"
        or raw["riemann_convention"] != "plus_partial_mu_gamma_nu"
    ):
        raise ValueError("ID2 identity or convention differs")
    paths = {
        name: _path(name, raw[name], expected)
        for name, expected in EXPECTED_PATHS.items()
    }
    if raw["scope"] != EXPECTED_SCOPE:
        raise ValueError("ID2 scope differs")
    if raw["physical_inputs"] != EXPECTED_PHYSICAL_INPUTS:
        raise ValueError("ID2 physical inputs differ")
    cases = raw["held_out_cases"]
    if not isinstance(cases, list) or [
        (
            item.get("case_id"),
            item.get("chi_width_factor"),
            item.get("phi_seed_factor"),
            item.get("exact_chi_center_buffer_over_L0"),
        )
        for item in cases
        if isinstance(item, dict)
    ] != list(EXPECTED_CASES):
        raise ValueError("ID2 held-out case partition differs")
    if raw["numerics"] != EXPECTED_NUMERICS:
        raise ValueError("ID2 numerical contract differs")
    if raw["thresholds"] != EXPECTED_THRESHOLDS:
        raise ValueError("ID2 threshold contract differs")
    proof = raw["proof_contract"]
    if set(proof) != EXPECTED_PROOF_KEYS or any(value is not True for value in proof.values()):
        raise ValueError("ID2 proof contract differs")
    if raw["claims"] != EXPECTED_CLAIMS:
        raise ValueError("ID2 claims differ")

    with paths["protocol_config"].open("rb") as handle:
        protocol_raw = tomllib.load(handle)
    protocol = validate_sf1_protocol(protocol_raw)
    with paths["physical_protocol_config"].open("rb") as handle:
        physical_protocol_raw = tomllib.load(handle)
    physical_protocol = validate_sf1_protocol(physical_protocol_raw)
    if (
        protocol.get("artifact_id") != SF1_PROTOCOL_V4_ARTIFACT_ID
        or physical_protocol.get("artifact_id") != SF1_PROTOCOL_ARTIFACT_ID
        or protocol_raw["amendment"]["predecessor_protocol_sha256"]
        != _sha(paths["physical_protocol_config"])
        or protocol["amplitude_candidates"] != raw["physical_inputs"]["amplitude_candidates"]
        or protocol["held_out_case_ids"] != [item["case_id"] for item in cases]
    ):
        raise ValueError("ID2 protocol inheritance differs")
    for name, (artifact, gate, branch) in RESULT_REQUIREMENTS.items():
        _check_result(paths[name], artifact, gate, branch)

    # This is not a dynamic run.  The check establishes that no pre-existing
    # file can be mistaken for evidence; it neither creates nor reads a result.
    for root_name in ("calibration_output_root", "holdout_output_root"):
        dynamic_root = REPOSITORY / protocol[root_name]
        if dynamic_root.exists():
            if not dynamic_root.is_dir() or any(dynamic_root.iterdir()):
                raise ValueError(f"ID2 found a nonempty {root_name}")

    physical = {
        "length_unit": float(_fraction("length_unit_L0", raw["physical_inputs"]["length_unit_L0"], positive=True)),
        "cutoff": float(_fraction("cutoff_Lambda", raw["physical_inputs"]["cutoff_Lambda"], positive=True)),
        "chi_center": float(_fraction("chi_center", raw["physical_inputs"]["chi_center"], positive=True)),
        "chi_base_half_width": float(_fraction("chi_base_half_width", raw["physical_inputs"]["chi_base_half_width"], positive=True)),
        "phi_half_width": float(_fraction("phi_half_width", raw["physical_inputs"]["phi_half_width"], positive=True)),
        "phi_seed_amplitude": float(_fraction("phi_seed_amplitude", raw["physical_inputs"]["phi_seed_amplitude"], positive=True)),
        "outer_radius": float(_fraction("outer_radius", raw["physical_inputs"]["outer_radius"], positive=True)),
        "measurement": float(_fraction("measurement_radius_maximum", raw["physical_inputs"]["measurement_radius_maximum"], positive=True)),
        "planck_mass": float(_fraction("planck_mass", raw["physical_inputs"]["planck_mass"], positive=True)),
        "scalar_mass": float(_fraction("scalar_mass", raw["physical_inputs"]["scalar_mass"], positive=True)),
        "quartic_coupling": float(_fraction("quartic_coupling", raw["physical_inputs"]["quartic_coupling"], positive=True)),
        "beta": float(_fraction("beta", raw["physical_inputs"]["beta"])),
        "eta": float(_fraction("eta", raw["physical_inputs"]["eta"], positive=True)),
    }
    physical["action"] = FGCQRActionParameters(
        planck_mass=physical["planck_mass"],
        beta=physical["beta"],
        scalar_mass=physical["scalar_mass"],
        quartic_coupling=physical["quartic_coupling"],
        eta=physical["eta"],
    )
    numerics = {
        "resolutions": tuple(raw["numerics"]["resolutions"]),
        "methods": (
            (raw["numerics"]["primary_constraint_method"], raw["numerics"]["primary_diagnostic_spatial_order"]),
            (raw["numerics"]["comparator_constraint_method"], raw["numerics"]["comparator_diagnostic_spatial_order"]),
        ),
        "taper": float(_fraction("proper_spectral_window_taper_fraction", raw["numerics"]["proper_spectral_window_taper_fraction"], positive=True)),
    }
    thresholds = {
        "compactness_minimum": float(_fraction("initial_compactness_minimum", raw["thresholds"]["initial_compactness_minimum"], positive=True)),
        "compactness_maximum": float(_fraction("initial_compactness_maximum", raw["thresholds"]["initial_compactness_maximum"], positive=True)),
        "maximum_constraint_residual": float(_fraction("maximum_physical_constraint_residual_infinity", raw["thresholds"]["maximum_physical_constraint_residual_infinity"], positive=True)),
        "minimum_RK4_order": float(_fraction("minimum_RK4_common_grid_profile_order", raw["thresholds"]["minimum_RK4_common_grid_profile_order"], positive=True)),
        "minimum_SSPRK3_order": float(_fraction("minimum_SSPRK3_common_grid_profile_order", raw["thresholds"]["minimum_SSPRK3_common_grid_profile_order"], positive=True)),
        "maximum_relative_mass_difference": float(_fraction("maximum_fine_cross_method_relative_outer_mass_difference", raw["thresholds"]["maximum_fine_cross_method_relative_outer_mass_difference"], positive=True)),
        "maximum_compactness_difference": float(_fraction("maximum_fine_cross_method_peak_compactness_difference", raw["thresholds"]["maximum_fine_cross_method_peak_compactness_difference"], positive=True)),
        "minimum_effective_planck_ratio": float(_fraction("minimum_effective_planck_ratio", raw["thresholds"]["minimum_effective_planck_ratio"], positive=True)),
        "minimum_reduction_order": float(_fraction("minimum_reduction_constraint_decay_order_diagnostic", raw["thresholds"]["minimum_reduction_constraint_decay_order_diagnostic"], positive=True)),
    }
    return {
        "raw": raw,
        "path": path,
        "paths": paths,
        "protocol": protocol,
        "physical_protocol": physical_protocol,
        "cases": tuple(cases),
        "physical": physical,
        "numerics": numerics,
        "thresholds": thresholds,
    }


def load_canonical_result(path: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    return _load_json(path, "ID2 result")


def record(path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    loaded = load_config(path)
    raw = loaded["raw"]
    configuration = {
        **loaded,
        "protocol_sha256": _sha(loaded["paths"]["protocol_config"]),
    }
    amplitudes = raw["physical_inputs"]["amplitude_candidates"]
    gr0_records = []
    for amplitude in amplitudes:
        result, _grids = _gr0_candidate(amplitude, configuration)
        gr0_records.append(result)
    fgcqr_records = [
        _fgcqr_case(amplitude, case, configuration)
        for amplitude in amplitudes
        for case in loaded["cases"]
    ]
    grouped = {
        amplitude: [item for item in fgcqr_records if item["amplitude"] == amplitude]
        for amplitude in amplitudes
    }
    eligibility = []
    for gr0 in gr0_records:
        cases = grouped[gr0["amplitude"]]
        all_cases = all(item["static_premises_passed"] for item in cases)
        eligible = gr0["GR0_static_input_passed"] and all_cases
        reasons = list(gr0["GR0_static_ineligibility_reasons"])
        if not all_cases:
            reasons.extend(
                f"{item['case_id']}:{reason}"
                for item in cases
                for reason in item["static_ineligibility_reasons"]
            )
        eligibility.append(
            {
                "amplitude": gr0["amplitude"],
                "GR0_static_input_passed": gr0["GR0_static_input_passed"],
                "all_five_FGCQR_static_slices_passed": all_cases,
                "eligible_for_fresh_dynamic_GR0_calibration": eligible,
                "ordered_ineligibility_reasons": reasons,
                "dynamic_outcome_used": False,
            }
        )
    eligible = [
        item["amplitude"]
        for item in eligibility
        if item["eligible_for_fresh_dynamic_GR0_calibration"]
    ]
    expected_eligible = ["5/2", "3"]
    if eligible != expected_eligible:
        raise RuntimeError(
            f"ID2 static eligibility changed: expected {expected_eligible}, observed {eligible}"
        )
    fgc_static_pass_count = sum(item["static_premises_passed"] for item in fgcqr_records)
    nested_diagnostic_pass_count = sum(
        item["nested_weighted_spectral_admission_diagnostic_passed"]
        for item in fgcqr_records
    )
    if fgc_static_pass_count != 31 or nested_diagnostic_pass_count != 17:
        raise RuntimeError("ID2 FGC-QR static or diagnostic cardinality changed")
    all_state_hashes = [
        resolution["state_sha256"]
        for record_group in (gr0_records, fgcqr_records)
        for record in record_group
        for method in record["methods"]
        for resolution in method["resolutions"]
    ]
    if len(all_state_hashes) != 252 or len(set(all_state_hashes)) != 252:
        raise RuntimeError("ID2 complete-grid state hashes are incomplete or collide")
    predecessor = {
        raw[name]: _sha(loaded["paths"][name]) for name in EXPECTED_PATHS
    }
    source_config = {_rel(loaded["path"]): _sha(loaded["path"])}
    implementation = {_rel(item): _sha(item) for item in IMPLEMENTATION}
    gate_status = {
        REQUIRED_GATE: True,
        "all_seven_GR0_candidates_statically_adjudicated": True,
        "all_thirty_five_FGCQR_slices_statically_adjudicated": True,
        "at_least_one_dynamic_calibration_candidate_remains": bool(eligible),
        "fresh_GR0_dynamic_calibration_candidate_set_frozen": True,
        **raw["claims"],
    }
    nonclaims = {
        "fresh_GR0_dynamic_calibration_run_performed": False,
        "trapped_sphere_formed_dynamically": False,
        "SGBL_control_evolved": False,
        "FGCQR_holdout_evolved": False,
        "regulator_activation_observed": False,
        "positive_Raychaudhuri_margin_observed": False,
        "candidate_branch_rejected": False,
        "general_gradient_route_rejected": False,
        "retained_EFT_validity": False,
        "physical_transition": False,
        "singularity_resolution": False,
        "child_domain_or_topology": False,
        "dark_sector_mechanism": False,
        "varying_locally_measured_c": False,
    }
    compactness_minima = [
        item["minimum_observed_peak_compactness"] for item in fgcqr_records
    ]
    compactness_maxima = [
        item["maximum_observed_peak_compactness"] for item in fgcqr_records
    ]
    result = {
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "project_version": PROJECT_VERSION,
        "classification": "pre_calibration_all_amplitude_all_case_static_initial_admission_certificate",
        "generated_by": "scripts/reproduce_fgc_id2_all1.py",
        "derivation_document": _rel(OWNER_DOCUMENT),
        "derivation_document_sha256": _sha(OWNER_DOCUMENT),
        "source_config_sha256": source_config,
        "predecessor_sha256": predecessor,
        "implementation_sha256": implementation,
        "scope_bindings": {
            "protocol_artifact_id": loaded["protocol"]["artifact_id"],
            "physical_predecessor_protocol_artifact_id": loaded[
                "physical_protocol"
            ]["artifact_id"],
            "protocol_config_sha256": _sha(loaded["paths"]["protocol_config"]),
            "physical_protocol_config_sha256": _sha(
                loaded["paths"]["physical_protocol_config"]
            ),
            "protocol_semantic_holdout_contract_sha256": loaded["protocol"][
                "semantic_holdout_contract_sha256"
            ],
            "calibration_branch": "GR-0",
            "holdout_branch": "FGC-QR",
            "physical_equations": "unredefined_ACT1_VAR1",
        },
        "gate_status": gate_status,
        "artifact_payload": {
            "protocol_partition": {
                "ordered_amplitude_candidates": amplitudes,
                "ordered_held_out_case_ids": [
                    item["case_id"] for item in loaded["cases"]
                ],
                "GR0_candidate_count": len(gr0_records),
                "FGCQR_static_slice_count": len(fgcqr_records),
                "method_count": 2,
                "resolution_count": 3,
                "complete_grid_state_hash_count": len(all_state_hashes),
                "every_complete_grid_state_hash_unique": True,
                "cartesian_product_completed_before_dynamic_calibration": True,
            },
            "GR0_candidate_ledger": gr0_records,
            "FGCQR_all_amplitude_all_case_ledger": fgcqr_records,
            "dynamic_calibration_eligibility": eligibility,
            "static_decision_summary": {
                "eligible_amplitudes_in_frozen_order": eligible,
                "ineligible_amplitudes_in_frozen_order": [
                    item["amplitude"] for item in eligibility if not item[
                        "eligible_for_fresh_dynamic_GR0_calibration"
                    ]
                ],
                "FGCQR_static_slice_pass_count": fgc_static_pass_count,
                "FGCQR_static_slice_fail_count": len(fgcqr_records)
                - fgc_static_pass_count,
                "FGCQR_nested_spectral_diagnostic_pass_count": nested_diagnostic_pass_count,
                "FGCQR_nested_spectral_diagnostic_fail_count": len(fgcqr_records)
                - nested_diagnostic_pass_count,
                "minimum_observed_FGCQR_peak_compactness": min(compactness_minima),
                "maximum_observed_FGCQR_peak_compactness": max(compactness_maxima),
                "amplitude_2_widest_case_below_compactness_floor": True,
                "amplitudes_4_and_9_over_2_narrow_case_miss_an_individual_spectral_budget": True,
                "amplitude_5_narrowest_case_above_compactness_ceiling": True,
                "nested_FGCQR_spectral_result_is_disclosed_but_not_a_PROTO3_static_five_boolean_decision": True,
            },
            "next_stage_contract": {
                "fresh_dynamic_GR0_candidates_must_run_in_this_order": eligible,
                "first_dynamically_eligible_candidate_stops_later_calibration_runs": True,
                "development_runs_cannot_enter_evidence": True,
                "candidate_config_and_result_hashes_must_be_fresh": True,
                "PROTO4_calibration_output_namespace": loaded["protocol"][
                    "calibration_output_root"
                ],
                "dynamic_calibration_not_performed_here": True,
            },
            "epistemic_boundary": {
                "dynamic_trajectory_used": False,
                "static_input_exclusion_is_not_FGCQR_mechanism_rejection": True,
                "two_compactness_exclusions_follow_the_frozen_input_window": True,
                "thresholds_changed_after_inspecting_static_results": False,
                "fresh_GR0_calibration_completed": False,
                "PRO4_HLD1_resolved": False,
                "SGBL_branch_health_definition_missing": True,
                "constraint_to_Raychaudhuri_stability_map_missing": True,
                "mechanism_question_answered": False,
            },
        },
        "nonclaims": nonclaims,
    }
    return _serial(result)


def reproduce(config: Path = DEFAULT_CONFIG, output: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    result = record(config)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(_canonical(result), encoding="utf-8")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    result = reproduce(args.config, args.output)
    print(
        f"wrote {args.output} ({result['artifact_id']}; "
        f"eligible={result['artifact_payload']['static_decision_summary']['eligible_amplitudes_in_frozen_order']})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
