#!/usr/bin/env python3
"""Reproduce the prospective FGC-1-RSP1-FRZ1 authorization record.

RSP1 asks one narrow question left open by CAL8/PREF10: whether amplitude
three's evolved ``phi/Lambda`` derivative tail enters the unchanged direct
``< 1/4`` contraction regime on the genuinely new
``4097 -> 8193 -> 16385`` ladder.  This pre-trajectory record rebuilds and
hashes all six inputs, exercises the SRC3 source instrument, publishes the
time-zero ratios without treating them as the endpoint outcome, and binds the
runner.  It neither advances nor reads an RSP1 trajectory.
"""

from __future__ import annotations

import argparse
from copy import deepcopy
from hashlib import sha256
from pathlib import Path
import subprocess
import sys
import tomllib
from typing import Any, Mapping

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts import reproduce_fgc_cal8_pref10 as cal8  # noqa: E402
from scripts import reproduce_fgc_hlt4_mon4 as inherited  # noqa: E402
from scripts import reproduce_fgc_src3 as src3  # noqa: E402
from recursive_horizons.fgc.evolution.gr0_calibration import (  # noqa: E402
    construct_gr0_grid_initial_data,
)
from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    EvolutionState,
    SBPFirstDerivative,
    UniformRadialGrid,
    array_content_sha256,
)
from recursive_horizons.fgc.evolution.proto11_runtime import (  # noqa: E402
    exact_spherical_minkowski_reference,
    project_gr0_reference_balanced_state,
    reference_balanced_spatial_derivatives,
)
from recursive_horizons.fgc.evolution.rsp1_resolution_runtime import (  # noqa: E402
    RSP1GR0EvolutionOperator,
    RSP1_POINT_COUNTS,
    RSP1_SOURCE_POINT_BATCH_SIZE,
    diagnose_gr0_reference_balanced_accelerations_chunked,
    rsp1_gr0_common_event,
)
from recursive_horizons.fgc.evolution.src3_reference_balanced_source import (  # noqa: E402
    diagnose_gr0_reference_balanced_accelerations,
)
from recursive_horizons.fgc.initial_data_preflight import (  # noqa: E402
    PulseParameters,
)


Q = inherited.Q
ARTIFACT_ID = "FGC-1-RSP1-FRZ1"
PROJECT_VERSION = "0.11.0"
CHECKPOINT_COMMIT = "4feda557dba06ccee751ea8db5a86e1925e767e4"
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-rsp1-frz1.toml"
DEFAULT_PLAN = REPOSITORY / "configs/fgc/fgc-1-rsp1-run1.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-rsp1-frz1.json"
OWNER_DOCUMENT = REPOSITORY / "docs/fgc-rsp1-frz1.md"

EXPECTED_SCOPE = {
    "role": "pre_trajectory_amplitude_three_resolution_spectrum_freeze_and_runtime_authorization",
    "branch": "GR-0",
    "amplitude": "3",
    "trajectory_advanced": False,
    "PROTO11_raw_outcomes_reinterpreted": False,
    "SGBL_trajectory_read": False,
    "FGCQR_trajectory_read": False,
    "collapse_or_trapped_outcome_classified": False,
    "mechanism_question_answered": False,
}

EXPECTED_CLAIMS = {
    "RSP1_runtime_implemented": True,
    "RSP1_execution_authorized": True,
    "amplitude_three_spectral_veto_cleared": False,
    "PROTO12_frozen": False,
    "fresh_GR0_calibration_completed": False,
    "classical_spherical_diagnostic_authorized": False,
    "SGBL_execution_authorized": False,
    "FGCQR_holdout_execution_authorized": False,
    "FGCQR_mechanism_rejected": False,
    "retained_EFT_evolution_authorized": False,
    "physical_transition_claim_authorized": False,
    "general_gradient_route_rejected": False,
    "singularity_resolution_derived": False,
    "child_domain_or_topology_derived": False,
    "dark_sector_mechanism_derived": False,
    "varying_locally_measured_c_derived": False,
}

EXPECTED_PLAN_CLAIMS = dict(EXPECTED_CLAIMS)
EXPECTED_PLAN_CLAIMS.pop("RSP1_runtime_implemented")
EXPECTED_PLAN_CLAIMS["RSP1_execution_authorized"] = False

EXPECTED_LINEAGE = {
    "checkpoint_commit": CHECKPOINT_COMMIT,
    "src3_result": "results/fgc-1-src3.json",
    "src3_result_sha256": "dea22654cb70b555053941009037c8b5c58e98dfc5a692dc27620e83363765e9",
    "src3_config": "configs/fgc/fgc-1-src3.toml",
    "src3_config_sha256": "9a0034ef739f75674d622223f42586c1abc56abb8afdbbd580d9143e22af61ba",
    "cal8_result": "results/fgc-1-cal8-pref10.json",
    "cal8_result_sha256": "10c6898e1ec6e7c07793180c4a4f8fa26a1b2f660115a0b2bc54497ff8bbb89d",
    "cal8_config": "configs/fgc/fgc-1-cal8-pref10.toml",
    "cal8_config_sha256": "a7646fb7ef2b229c36981a2f4a87e4428de1bdc5cbe324baedc965097c674cb3",
    "checkpoint_must_be_ancestor_of_HEAD": True,
    "tracked_lineage_must_be_canonical_and_hash_matched": True,
}

EXPECTED_STUDY_FREEZE = {
    "target_coordinate_time": "1/16",
    "amplitude": "3",
    "methods_in_order": ["RK4", "SSPRK3"],
    "point_counts_in_order": list(RSP1_POINT_COUNTS),
    "new_point_count": 16385,
    "source_point_batch_size": RSP1_SOURCE_POINT_BATCH_SIZE,
    "target_field": "phi_over_Lambda",
    "target_quantity": "top_band_derivative_weighted_power_fraction",
    "unchanged_direct_ratio_ceiling": "1/4",
    "both_target_pairs_and_both_methods_must_pass": True,
    "target_saturation_or_epsilon_floor_forbidden": True,
    "all_source_constraint_absolute_budget_and_profile_guards_remain_vetoes": True,
    "t0_target_ratios_are_public_but_not_the_evolved_endpoint_outcome": True,
    "only_one_synchronized_evolved_endpoint_may_be_read": True,
    "no_parameter_may_change_after_the_first_trajectory": True,
}

EXPECTED_RUNTIME_BINDING_KEYS = {
    "SRC3_reference_covariant_source_required",
    "PROTO11_reference_state_derivative_map_required",
    "bounded_point_batching_changes_no_per_point_equation_or_root",
    "exact_Minkowski_must_remain_bitwise_stationary",
    "chunked_and_unchunked_source_must_agree_on_frozen_controls",
    "complete_unredefined_raw_source_residual_remains_serialized",
    "raw_source_threshold_remains_strictly_below_1e_minus_12",
    "condition_limit_remains_strictly_below_1e10",
    "both_methods_use_their_original_integrator_and_SBP_order",
    "transaction_retry_boundary_and_health_rules_unchanged",
}

EXPECTED_PROOF_KEYS = {
    "all_six_initial_states_must_be_rebuilt_and_hash_bound",
    "all_six_t0_SRC3_source_prechecks_must_pass",
    "both_t0_constraint_admissions_must_pass",
    "all_t0_medium_and_fine_absolute_spectral_budgets_must_pass",
    "t0_target_profile_must_contract_for_both_methods",
    "t0_direct_target_ratio_failures_must_remain_public",
    "one_bit_input_mutation_must_change_the_chunked_source",
    "invalid_batch_grid_threshold_and_claim_mutations_must_fail_closed",
    "runner_must_refuse_dirty_tracked_state_drift_or_namespace_overwrite",
    "canonical_hash_bound_authorization_required",
}

EXPECTED_NAMESPACE = {
    "output_root": "runs/fgc-2-sf1/rsp1/amplitude3-spectrum",
    "runner_id": "FGC-1-RSP1-RUN1-RUNNER",
    "fresh_output_root_required_at_first_launch": True,
    "authorization_creates_no_run_output": True,
    "raw_outputs_remain_ignored_and_separate_from_tracked_certificate": True,
}

IMPLEMENTATION = (
    REPOSITORY / "src/recursive_horizons/fgc/evolution/rsp1_resolution_runtime.py",
    REPOSITORY / "scripts/run_fgc_rsp1_resolution_study.py",
    Path(__file__).resolve(),
)


def _rel(path: Path) -> str:
    return path.resolve().relative_to(REPOSITORY).as_posix()


def _sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _strict_keys(name: str, value: Mapping[str, Any], expected: set[str]) -> None:
    if set(value) != expected:
        raise ValueError(f"{name} keys differ")


def _all_true(name: str, value: Mapping[str, Any], expected: set[str]) -> None:
    if set(value) != expected or any(item is not True for item in value.values()):
        raise ValueError(f"{name} differs")


def validate_run_plan(path: Path = DEFAULT_PLAN) -> dict[str, Any]:
    """Strictly validate the outcome-neutral RSP1 trajectory plan."""

    with path.open("rb") as handle:
        raw = tomllib.load(handle)
    _strict_keys(
        "RSP1 run plan",
        raw,
        {
            "schema_version",
            "artifact_id",
            "project_version",
            "metric_signature",
            "riemann_convention",
            "freeze_result",
            "scope",
            "question",
            "physical_inputs",
            "primary_method",
            "comparator_method",
            "numerics",
            "universal_thresholds",
            "endpoint_classification",
            "provenance",
            "claims",
        },
    )
    if (
        raw["schema_version"] != 1
        or raw["artifact_id"] != "FGC-1-RSP1-RUN1-PLAN"
        or raw["project_version"] != PROJECT_VERSION
        or raw["metric_signature"] != "-+++"
        or raw["riemann_convention"] != "plus_partial_mu_gamma_nu"
        or raw["freeze_result"] != _rel(DEFAULT_OUTPUT)
    ):
        raise ValueError("RSP1 run-plan identity differs")
    if raw["scope"] != {
        "branch": "GR-0",
        "amplitude": "3",
        "role": "prospective_single_endpoint_resolution_spectrum_adjudication",
        "continuum_equations": "unredefined_GR0_specialization_of_ACT1_VAR1_REF1",
        "source_instrument": "FGC-1-SRC3",
        "reference_state_map": "PROTO11",
        "target_coordinate_time": "1/16",
        "calibration_or_candidate_selection": False,
        "collapse_or_trapped_outcome_classified": False,
        "SGBL_or_FGCQR_state_read": False,
        "retained_EFT_interpretation": False,
    }:
        raise ValueError("RSP1 run-plan scope differs")
    if raw["question"] != {
        "target_field": "phi_over_Lambda",
        "target_quantity": "top_band_derivative_weighted_power_fraction",
        "point_counts": list(RSP1_POINT_COUNTS),
        "new_point_count": 16385,
        "adjacent_pairs": ["4097_to_8193", "8193_to_16385"],
        "unchanged_strict_nested_ratio_ceiling": "1/4",
        "both_adjacent_target_ratios_must_pass_directly": True,
        "diagnostic_saturation_for_target_ratios_forbidden": True,
        "medium_and_fine_target_absolute_budgets_required": True,
        "all_medium_and_fine_field_absolute_budgets_required": True,
        "target_whole_profile_contraction_required": True,
        "both_methods_required": True,
        "outcome_may_not_change_the_frozen_question": True,
    }:
        raise ValueError("RSP1 frozen question differs")
    if raw["physical_inputs"] != {
        "length_unit_L0": "4",
        "cutoff_Lambda": "16",
        "planck_mass": "2",
        "scalar_mass": "3",
        "quartic_coupling": "1/2",
        "chi_amplitude": "3",
        "chi_center": "12",
        "chi_half_width": "2",
        "phi_seed_amplitude": "1/131072",
        "outer_radius": "128",
        "measurement_radius_maximum": "24",
    }:
        raise ValueError("RSP1 physical inputs differ")
    if raw["primary_method"] != {
        "method_label": "RK4",
        "integrator_id": "fourth_order_diagonal_norm_SBP_with_second_order_boundary_closure_plus_RK4",
        "constraint_solve_method": "RK4",
        "spatial_order": 4,
        "coarsest_constraint_guard": "1/50",
        "finest_constraint_guard": "1/1000",
    } or raw["comparator_method"] != {
        "method_label": "SSPRK3",
        "integrator_id": "second_order_diagonal_norm_SBP_plus_SSPRK3",
        "constraint_solve_method": "SSPRK3",
        "spatial_order": 2,
        "coarsest_constraint_guard": "1/10",
        "finest_constraint_guard": "1/200",
    }:
        raise ValueError("RSP1 numerical methods differ")
    numerics = raw["numerics"]
    if (
        numerics.get("resolutions") != list(RSP1_POINT_COUNTS)
        or numerics.get("source_point_batch_size") != RSP1_SOURCE_POINT_BATCH_SIZE
        or numerics.get("final_coordinate_time") != "1/16"
        or numerics.get("common_output_interval") != "1/16"
        or numerics.get("checkpoint_interval") != "1/16"
        or numerics.get("minimum_constraint_finest_pair_order") != "3/2"
        or numerics.get("initial_q_reference_state_rule")
        != "D_h(u-u_ref)+q_ref"
        or numerics.get("source_q_r_reference_state_rule") != "D_h(q-q_ref)"
        or numerics.get("reduction_constraint_reference_state_rule")
        != "q-(D_h(u-u_ref)+q_ref)"
        or numerics.get("q_remains_independently_evolved") is not True
        or numerics.get(
            "interior_q_reprojection_after_initialization_or_Runge_Kutta_stage"
        )
        is not False
        or numerics.get("monitor_every_proposed_Runge_Kutta_stage") is not True
        or numerics.get("accepted_state_source_precheck_before_every_proposal")
        is not True
        or numerics.get("serialize_every_rejected_source_trial") is not True
        or numerics.get("common_event_endpoint_is_bitwise_target_time") is not True
        or numerics.get("no_excision_in_measurement_region") is not True
    ):
        raise ValueError("RSP1 numerical contract differs")
    if raw["universal_thresholds"] != {
        "source_residual_infinity_max": "1/1000000000000",
        "source_iteration_max": 16,
        "source_residual_must_decrease_monotonically": True,
        "kinetic_condition_number_max": "10000000000",
        "hat_normal_factor": "9",
        "minimum_boundary_causal_buffer": "16",
    }:
        raise ValueError("RSP1 source or health thresholds differ")
    if raw["endpoint_classification"] != {
        "positive": "completed_direct_amplitude_three_spectrum_pass",
        "negative": "completed_direct_amplitude_three_spectrum_obstruction",
        "runtime_stop": "stopped_before_resolution_endpoint",
        "invalid": "invalid_implementation_or_nonconverged_run",
        "positive_clears_only_the_CAL8_amplitude_three_veto_for_successor_design": True,
        "negative_preserves_the_veto_and_localizes_the_failed_pair": True,
        "neither_outcome_is_a_GR0_calibration_or_mechanism_result": True,
    }:
        raise ValueError("RSP1 endpoint classification differs")
    if raw["provenance"] != {
        "output_root": EXPECTED_NAMESPACE["output_root"],
        "campaign_manifest": EXPECTED_NAMESPACE["output_root"] + "/manifest.json",
        "append_only_event_log_name": "events.jsonl",
        "checkpoint_extension": ".npz",
        "fresh_output_root_required_at_study_start": True,
        "runner_refuses_overwrite": True,
        "resume_requires_matching_plan_authorization_and_manifest_hashes": True,
        "PROTO11_outputs_are_immutable_predecessor_evidence_not_RSP1_outcomes": True,
    } or raw["claims"] != EXPECTED_PLAN_CLAIMS:
        raise ValueError("RSP1 provenance or pre-authorization claims differ")
    return raw


def _validate_config_mapping(raw: Mapping[str, Any]) -> None:
    _strict_keys(
        "RSP1 freeze config",
        raw,
        {
            "schema_version",
            "artifact_id",
            "project_version",
            "metric_signature",
            "riemann_convention",
            "run_plan",
            "owner_document",
            "scope",
            "lineage",
            "study_freeze",
            "runtime_binding",
            "proof_contract",
            "namespace",
            "claims",
        },
    )
    if (
        raw["schema_version"] != 1
        or raw["artifact_id"] != ARTIFACT_ID
        or raw["project_version"] != PROJECT_VERSION
        or raw["metric_signature"] != "-+++"
        or raw["riemann_convention"] != "plus_partial_mu_gamma_nu"
        or raw["run_plan"] != _rel(DEFAULT_PLAN)
        or raw["owner_document"] != _rel(OWNER_DOCUMENT)
        or raw["scope"] != EXPECTED_SCOPE
        or raw["lineage"] != EXPECTED_LINEAGE
        or raw["study_freeze"] != EXPECTED_STUDY_FREEZE
        or raw["namespace"] != EXPECTED_NAMESPACE
        or raw["claims"] != EXPECTED_CLAIMS
    ):
        raise ValueError("RSP1 freeze identity, scope, lineage, or claims differ")
    _all_true("RSP1 runtime binding", raw["runtime_binding"], EXPECTED_RUNTIME_BINDING_KEYS)
    _all_true("RSP1 proof contract", raw["proof_contract"], EXPECTED_PROOF_KEYS)


def load_config(path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    with path.open("rb") as handle:
        raw = tomllib.load(handle)
    _validate_config_mapping(raw)
    return raw


def _load_canonical(path: Path, artifact_id: str) -> dict[str, Any]:
    value = inherited._load_json_bytes(path.read_bytes(), artifact_id)
    if value.get("artifact_id") != artifact_id:
        raise ValueError(f"{_rel(path)} artifact identity differs")
    return value


def _validate_lineage(config: Mapping[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    lineage = config["lineage"]
    process = subprocess.run(
        ["git", "merge-base", "--is-ancestor", CHECKPOINT_COMMIT, "HEAD"],
        cwd=REPOSITORY,
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if process.returncode != 0:
        raise ValueError("RSP1 checkpoint is not an ancestor of HEAD")
    for path_key, hash_key in (
        ("src3_result", "src3_result_sha256"),
        ("src3_config", "src3_config_sha256"),
        ("cal8_result", "cal8_result_sha256"),
        ("cal8_config", "cal8_config_sha256"),
    ):
        path = REPOSITORY / lineage[path_key]
        if not path.is_file() or _sha(path) != lineage[hash_key]:
            raise ValueError(f"RSP1 lineage hash differs: {lineage[path_key]}")
    src3_result = _load_canonical(REPOSITORY / lineage["src3_result"], "FGC-1-SRC3")
    cal8_result = _load_canonical(
        REPOSITORY / lineage["cal8_result"], "FGC-1-CAL8-PREF10"
    )
    src3_gates = src3_result["gate_status"]
    cal8_gates = cal8_result["gate_status"]
    diagnosis = cal8_result["artifact_payload"]["campaign_diagnosis"]
    amplitude_three = diagnosis["amplitude_three_direct_spectral_obstruction"]
    historical = amplitude_three["methods"]
    if (
        src3_gates.get(
            "SRC3_runtime_source_instrument_authorized_for_a_separately_frozen_GR0_protocol"
        )
        is not True
        or src3_gates.get("amplitude_three_spectral_veto_cleared") is not False
        or cal8_gates.get("PROTO11_campaign_terminated_normally") is not True
        or cal8_gates.get("PROTO11_GR0_case_eligible") is not False
        or cal8_gates.get("PROTO11_direct_coarse_phi_derivative_veto_localized")
        is not True
        or amplitude_three.get("amplitude") != "3"
        or amplitude_three.get("source_only_rejection_count") != 0
        or amplitude_three.get("sole_direct_coarse_to_medium_veto")
        != "phi_over_Lambda derivative tail"
        or historical["RK4"].get("phi_derivative_tail_ratios")
        != [0.326717281113685, 0.22812049822366257]
        or historical["SSPRK3"].get("phi_derivative_tail_ratios")
        != [0.26675704013240803, 0.2142551901915683]
        or any(value is not False for value in src3_result.get("nonclaims", {}).values())
        or any(value is not False for value in cal8_result.get("nonclaims", {}).values())
    ):
        raise ValueError("RSP1 predecessor decisions differ")
    return src3_result, cal8_result


def _parameters(plan: Mapping[str, Any]) -> PulseParameters:
    physical = plan["physical_inputs"]
    return PulseParameters(
        chi_amplitude=float(Q(physical["chi_amplitude"])),
        center=float(Q(physical["chi_center"])),
        half_width=float(Q(physical["chi_half_width"])),
        phi_amplitude=float(Q(physical["phi_seed_amplitude"])),
        planck_mass=float(Q(physical["planck_mass"])),
        scalar_mass=float(Q(physical["scalar_mass"])),
        quartic_coupling=float(Q(physical["quartic_coupling"])),
    )


def _operator(
    plan: Mapping[str, Any], grid: UniformRadialGrid, order: int
) -> RSP1GR0EvolutionOperator:
    numerics = plan["numerics"]
    thresholds = plan["universal_thresholds"]
    return RSP1GR0EvolutionOperator(
        grid,
        spatial_order=order,
        ko_dissipation=float(Q(numerics["ko_dissipation"])),
        raw_tolerance=float(Q(thresholds["source_residual_infinity_max"])),
        kinetic_condition_maximum=float(
            Q(thresholds["kinetic_condition_number_max"])
        ),
        point_batch_size=numerics["source_point_batch_size"],
    )


def _freeze_inputs(plan: Mapping[str, Any], plan_hash: str):
    records: list[dict[str, Any]] = []
    prechecks: list[dict[str, Any]] = []
    states: dict[str, list[tuple[EvolutionState, UniformRadialGrid]]] = {}
    for method_key in ("primary_method", "comparator_method"):
        method = plan[method_key]
        label = method["method_label"]
        order = method["spatial_order"]
        for point_count in RSP1_POINT_COUNTS:
            initial = construct_gr0_grid_initial_data(
                _parameters(plan),
                point_count=point_count,
                outer_radius=float(Q(plan["physical_inputs"]["outer_radius"])),
                constraint_method=method["constraint_solve_method"],
                diagnostic_spatial_order=order,
            )
            state = project_gr0_reference_balanced_state(
                initial, spatial_order=order
            )
            before = array_content_sha256(state.u, state.p, state.q)
            rhs = _operator(plan, initial.grid, order)(0.0, state)
            after = array_content_sha256(state.u, state.p, state.q)
            diagnostics = rhs.diagnostics
            if (
                before != after
                or diagnostics["source_raw_gate_passed"] is not True
                or diagnostics["SRC3_reference_covariant_source"] is not True
                or diagnostics["PROTO11_reference_state_map_retained"] is not True
                or diagnostics["PROTO11_exact_reference_equilibrium_applied"]
                is not False
                or diagnostics["PROTO11_interior_q_reprojected"] is not False
            ):
                raise ValueError("RSP1 initial source precheck failed")
            record = {
                "amplitude": "3",
                "method": label,
                "integrator_id": method["integrator_id"],
                "spatial_order": order,
                "point_count": point_count,
                "grid_spacing": initial.grid.spacing,
                "plan_sha256": plan_hash,
                "projected_state_sha256": before,
                "source_instrument": "FGC-1-SRC3",
                "reference_state_map": "PROTO11",
                "q_remains_independently_evolved": True,
                "interior_q_reprojected": False,
                "trajectory_advanced": False,
            }
            record["expanded_run_config_sha256"] = inherited._digest(record)
            records.append(record)
            prechecks.append(
                {
                    "amplitude": "3",
                    "method": label,
                    "point_count": point_count,
                    "source_residual_infinity": diagnostics[
                        "source_residual_infinity"
                    ],
                    "source_residual_maximum": float(
                        Q(
                            plan["universal_thresholds"][
                                "source_residual_infinity_max"
                            ]
                        )
                    ),
                    "kinetic_condition_infinity": diagnostics[
                        "kinetic_condition_infinity"
                    ],
                    "source_point_batch_size": diagnostics[
                        "RSP1_source_point_batch_size"
                    ],
                    "source_point_batch_count": diagnostics[
                        "RSP1_source_point_batch_count"
                    ],
                    "accepted_state_source_gate_passed": True,
                    "operator_input_state_unchanged": True,
                    "trajectory_advanced": False,
                }
            )
            states.setdefault(label, []).append((state, initial.grid))
    return records, prechecks, states


def _initial_events(plan: Mapping[str, Any], states):
    events = []
    for method in ("RK4", "SSPRK3"):
        records = states[method]
        assessment = rsp1_gr0_common_event(
            [item[0] for item in records],
            [item[1] for item in records],
            accepted_stage_counts=(0, 0, 0),
            method=method,
            coordinate_time=0.0,
            cutoff=float(Q(plan["physical_inputs"]["cutoff_Lambda"])),
            measurement_radius_maximum=float(
                Q(plan["physical_inputs"]["measurement_radius_maximum"])
            ),
            taper_fraction=float(
                Q(plan["numerics"]["proper_spectral_taper_fraction"])
            ),
            fixed_outer_rows=plan["numerics"]["fixed_outer_rows"],
            point_batch_size=plan["numerics"]["source_point_batch_size"],
        )
        ratios = assessment.target_derivative_tail_ratios
        if (
            assessment.constraint_admission.admission_passed is not True
            or assessment.target_medium_and_fine_absolute_budgets_passed is not True
            or assessment.all_medium_and_fine_absolute_budgets_passed is not True
            or assessment.target_complete_profile_contraction_passed is not True
            or assessment.target_direct_pair_passed != (False, False)
            or not all(value is not None and value > 0.25 for value in ratios)
            or assessment.study_admission_passed is not False
        ):
            raise ValueError("RSP1 time-zero premise classification differs")
        events.append(
            {
                "amplitude": "3",
                "method": method,
                "coordinate_time": 0.0,
                "point_counts": list(RSP1_POINT_COUNTS),
                "constraint_admission_passed": True,
                "all_medium_and_fine_absolute_budgets_passed": True,
                "target_complete_profile_contraction_passed": True,
                "target_derivative_tail_ratios": list(ratios),
                "target_direct_pair_passed": [False, False],
                "study_admission_passed": False,
                "t0_target_failure_is_the_evolved_endpoint_outcome": False,
                "complete_common_event": inherited._serial(assessment),
                "trajectory_advanced": False,
            }
        )
    return events


def _controls(plan: Mapping[str, Any], states) -> dict[str, Any]:
    state, grid = states["RK4"][0]
    derivative = SBPFirstDerivative(grid, 4)
    p_r, q_r, _ = reference_balanced_spatial_derivatives(state, derivative)
    sample = slice(1, 129)
    args = (
        state.u[sample],
        state.p[sample],
        state.q[sample],
        p_r[sample],
        q_r[sample],
        grid.coordinates[sample],
    )
    direct = diagnose_gr0_reference_balanced_accelerations(*args)
    chunked = diagnose_gr0_reference_balanced_accelerations_chunked(
        *args, point_batch_size=64
    )
    if (
        not np.array_equal(chunked.accelerations, direct.accelerations)
        or not np.array_equal(chunked.residuals, direct.residuals)
        or chunked.residual_infinity != direct.residual_infinity
    ):
        raise ValueError("RSP1 chunking changed a pointwise source solve")

    mutated_u = state.u[sample].copy()
    # The frozen first positive node has exactly zero shift.  Advancing that
    # binary64 value by one representable step supplies a literal one-bit-edge
    # perturbation without changing any configured tolerance.
    mutated_u[0, 1] = np.nextafter(mutated_u[0, 1], np.inf)
    mutated = diagnose_gr0_reference_balanced_accelerations_chunked(
        mutated_u, *args[1:], point_batch_size=64
    )
    one_bit_visible = not np.array_equal(
        chunked.accelerations, mutated.accelerations
    )
    if not one_bit_visible:
        raise ValueError("RSP1 one-bit source mutation was hidden")

    exact_records = []
    for method, order in (("RK4", 4), ("SSPRK3", 2)):
        exact_grid = UniformRadialGrid(0.0, 128.0, RSP1_POINT_COUNTS[0])
        exact = exact_spherical_minkowski_reference(exact_grid)
        before = array_content_sha256(exact.u, exact.p, exact.q)
        rhs = _operator(plan, exact_grid, order)(0.0, exact)
        after = array_content_sha256(exact.u, exact.p, exact.q)
        if (
            before != after
            or not np.array_equal(rhs.du, np.zeros_like(rhs.du))
            or not np.array_equal(rhs.dp, np.zeros_like(rhs.dp))
            or not np.array_equal(rhs.dq, np.zeros_like(rhs.dq))
            or rhs.diagnostics["source_residual_infinity"] != 0.0
            or rhs.diagnostics["PROTO11_exact_reference_equilibrium_applied"]
            is not True
        ):
            raise ValueError("RSP1 exact-reference control failed")
        exact_records.append(
            {
                "method": method,
                "point_count": exact_grid.point_count,
                "source_residual_infinity": 0.0,
                "du_dp_dq_bitwise_zero": True,
                "operator_input_unchanged": True,
            }
        )

    typed_failures: dict[str, bool] = {}
    try:
        diagnose_gr0_reference_balanced_accelerations_chunked(
            *args, point_batch_size=0
        )
    except ValueError:
        typed_failures["invalid_batch_size"] = True
    try:
        RSP1GR0EvolutionOperator(
            UniformRadialGrid(1.0, 2.0, 17),
            spatial_order=4,
            ko_dissipation=0.0,
            raw_tolerance=1.0e-12,
            kinetic_condition_maximum=1.0e10,
        )
    except ValueError:
        typed_failures["noncentre_grid"] = True
    try:
        RSP1GR0EvolutionOperator(
            grid,
            spatial_order=4,
            ko_dissipation=0.0,
            raw_tolerance=0.0,
            kinetic_condition_maximum=1.0e10,
        )
    except ValueError:
        typed_failures["invalid_source_threshold"] = True
    attacked = deepcopy(load_config(DEFAULT_CONFIG))
    attacked["claims"]["FGCQR_holdout_execution_authorized"] = True
    try:
        _validate_config_mapping(attacked)
    except ValueError:
        typed_failures["claim_promotion"] = True
    if typed_failures != {
        "invalid_batch_size": True,
        "noncentre_grid": True,
        "invalid_source_threshold": True,
        "claim_promotion": True,
    }:
        raise ValueError("RSP1 typed-failure controls are incomplete")
    return {
        "chunked_control_point_count": 128,
        "chunked_batch_size": 64,
        "chunked_and_direct_accelerations_bitwise_equal": True,
        "chunked_and_direct_residuals_bitwise_equal": True,
        "chunked_complete_residual_infinity": chunked.residual_infinity,
        "direct_complete_residual_infinity": direct.residual_infinity,
        "one_bit_shift_mutation_visible": True,
        "baseline_acceleration_sha256": array_content_sha256(
            chunked.accelerations
        ),
        "mutated_acceleration_sha256": array_content_sha256(
            mutated.accelerations
        ),
        "exact_reference_records": exact_records,
        "typed_failure_controls": typed_failures,
        "all_controls_passed": True,
    }


def _namespace_precondition(
    config: Mapping[str, Any], historical: Mapping[str, Any] | None
) -> dict[str, Any]:
    relative = config["namespace"]["output_root"]
    if historical is not None:
        if historical != {
            "path": relative,
            "empty_before_authorization": True,
            "authorization_created_no_namespace": True,
        }:
            raise ValueError("historical RSP1 namespace evidence differs")
        return dict(historical)
    path = REPOSITORY / relative
    if path.exists() and (not path.is_dir() or any(path.iterdir())):
        raise ValueError("RSP1 output namespace is not fresh")
    return {
        "path": relative,
        "empty_before_authorization": True,
        "authorization_created_no_namespace": True,
    }


def reproduce(
    config_path: Path = DEFAULT_CONFIG,
    *,
    historical_namespace_evidence: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    config = load_config(config_path)
    plan_path = REPOSITORY / config["run_plan"]
    plan = validate_run_plan(plan_path)
    src3_result, cal8_result = _validate_lineage(config)
    plan_hash = _sha(plan_path)
    inputs, prechecks, states = _freeze_inputs(plan, plan_hash)
    events = _initial_events(plan, states)
    controls = _controls(plan, states)
    namespace = _namespace_precondition(config, historical_namespace_evidence)
    if len(inputs) != 6 or len(prechecks) != 6 or len(events) != 2:
        raise ValueError("RSP1 authorization evidence is incomplete")
    input_manifest_hash = inherited._digest(inputs)
    study_id = inherited._digest(
        {
            "artifact_id": plan["artifact_id"],
            "run_plan_sha256": plan_hash,
            "input_manifest_sha256": input_manifest_hash,
            "amplitude": "3",
            "methods": ["RK4", "SSPRK3"],
            "point_counts": list(RSP1_POINT_COUNTS),
            "target_coordinate_time": "1/16",
            "output_root": config["namespace"]["output_root"],
        }
    )
    implementation = dict(cal8_result.get("implementation_sha256", {}))
    implementation.update(src3_result.get("implementation_sha256", {}))
    implementation.update({_rel(path): _sha(path) for path in IMPLEMENTATION})
    for relative, expected in implementation.items():
        path = REPOSITORY / relative
        if not path.is_file() or _sha(path) != expected:
            raise ValueError(f"RSP1 inherited implementation drifted: {relative}")
    return {
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "project_version": PROJECT_VERSION,
        "classification": "pre_trajectory_amplitude_three_resolution_spectrum_runtime_authorization",
        "generated_by": _rel(Path(__file__)),
        "derivation_document": _rel(OWNER_DOCUMENT),
        "derivation_document_sha256": _sha(OWNER_DOCUMENT),
        "source_config_sha256": {
            _rel(config_path): _sha(config_path),
            _rel(plan_path): plan_hash,
        },
        "predecessor_sha256": {
            config["lineage"]["src3_result"]: config["lineage"][
                "src3_result_sha256"
            ],
            config["lineage"]["cal8_result"]: config["lineage"][
                "cal8_result_sha256"
            ],
        },
        "implementation_sha256": implementation,
        "scope_bindings": dict(config["scope"]),
        "gate_status": dict(EXPECTED_CLAIMS),
        "artifact_payload": {
            "immutable_lineage": {
                "checkpoint_commit": CHECKPOINT_COMMIT,
                "checkpoint_is_ancestor_of_HEAD": True,
                "SRC3_canonical_hash_matched": True,
                "CAL8_canonical_hash_matched": True,
                "CAL8_amplitude_three_veto_preserved": True,
                "historical_RK4_phi_derivative_tail_ratios": [
                    0.326717281113685,
                    0.22812049822366257,
                ],
                "historical_SSPRK3_phi_derivative_tail_ratios": [
                    0.26675704013240803,
                    0.2142551901915683,
                ],
            },
            "frozen_run_plan": {
                "artifact_id": plan["artifact_id"],
                "run_plan_sha256": plan_hash,
                "study_id": study_id,
                "input_manifest_sha256": input_manifest_hash,
                "amplitude": "3",
                "methods": ["RK4", "SSPRK3"],
                "point_counts": list(RSP1_POINT_COUNTS),
                "target_coordinate_time": "1/16",
                "target_field": "phi_over_Lambda",
                "strict_direct_ratio_ceiling": 0.25,
                "expanded_input_count": 6,
                "source_precheck_passed_count": 6,
                "t0_premise_event_count": 2,
            },
            "frozen_run_inputs": inputs,
            "accepted_state_source_prechecks": prechecks,
            "initial_premise_events": events,
            "runtime_controls": controls,
            "namespace_precondition": namespace,
            "decision": (
                "RSP1_may_run_one_frozen_amplitude_three_three_resolution_study_"
                "while_every_calibration_candidate_and_physical_claim_remains_closed"
            ),
            "epistemic_boundary": {
                "trajectory_read": False,
                "trajectory_advanced": False,
                "t0_target_ratios_are_public": True,
                "t0_target_ratios_are_the_evolved_endpoint_outcome": False,
                "RSP1_execution_authorized": True,
                "RSP1_execution_completed": False,
                "amplitude_three_spectral_veto_cleared": False,
                "fresh_GR0_calibration_completed": False,
                "SGBL_or_FGCQR_outcome_read": False,
                "mechanism_question_answered": False,
            },
        },
        "nonclaims": {
            key: value for key, value in EXPECTED_CLAIMS.items() if value is False
        },
    }


def load_canonical_result(path: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    return _load_canonical(path, ARTIFACT_ID)


def verify_canonical(
    path: Path = DEFAULT_OUTPUT, config_path: Path = DEFAULT_CONFIG
) -> None:
    observed = load_canonical_result(path)
    historical = observed["artifact_payload"]["namespace_precondition"]
    expected = reproduce(
        config_path,
        historical_namespace_evidence=historical,
    )
    if observed != inherited._serial(expected):
        raise ValueError("RSP1 authorization differs from a fresh reproduction")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    config = args.config.resolve()
    output = args.output.resolve()
    if args.check:
        verify_canonical(output, config)
        print(f"verified {output}")
        return
    record = reproduce(config)
    output.write_text(inherited._canonical(record), encoding="utf-8")
    print(
        f"wrote {output} "
        "(RSP1 runtime=true; execution authorized=true; trajectory_read=false)"
    )


if __name__ == "__main__":
    main()
