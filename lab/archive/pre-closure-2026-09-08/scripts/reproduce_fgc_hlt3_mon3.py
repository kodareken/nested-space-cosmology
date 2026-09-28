#!/usr/bin/env python3
"""Regenerate the pre-trajectory FGC-1-HLT3-MON3 authorization record."""

from __future__ import annotations

import argparse
from dataclasses import asdict, is_dataclass
from fractions import Fraction
from hashlib import sha256
import json
from math import isfinite
from pathlib import Path
import subprocess
import sys
import tomllib
from typing import Any, Mapping, Sequence

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-hlt3-mon3.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-hlt3-mon3.json"
OWNER_DOCUMENT = REPOSITORY / "docs/fgc-hlt3-mon3.md"

from recursive_horizons.fgc.evolution.boundary_domain import (  # noqa: E402
    BoundaryGeometry,
    CausalBudgetState,
)
from recursive_horizons.fgc.evolution.calibration_runtime import (  # noqa: E402
    GR0_UNIVERSAL_RUNTIME_STOP_IDS,
    project_gr0_semidiscrete_state,
)
from recursive_horizons.fgc.evolution.gr0_calibration import (  # noqa: E402
    construct_gr0_grid_initial_data,
)
from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    COMPARATOR_METHOD,
    PRIMARY_METHOD,
    EvolutionRHS,
    EvolutionState,
    array_content_sha256,
    propose_step,
)
from recursive_horizons.fgc.evolution.proto4_admission import (  # noqa: E402
    PROTO4_SPECTRAL_FIELD_ORDER,
)
from recursive_horizons.fgc.evolution.proto5_runtime import (  # noqa: E402
    CFLRetryRequired,
    GR0RuntimeStageTransaction,
    GR0RuntimeStop,
    GR0UniversalThresholds,
    GR0_RUNTIME_STOP_PRIORITY,
    PROTO5_METHOD_CONTRACT,
    proto5_gr0_common_event,
    proto5_event_log_recovery_plan,
    proto5_temporal_spectral_admission,
)
from recursive_horizons.fgc.evolution.protocol_v5 import (  # noqa: E402
    SF1_PROTOCOL_V5_ARTIFACT_ID,
    validate_sf1_protocol_v5,
)
from recursive_horizons.fgc.initial_data_preflight import (  # noqa: E402
    PulseParameters,
)


Q = Fraction
ARTIFACT_ID = "FGC-1-HLT3-MON3"
PROJECT_VERSION = "0.11.0"
EXPECTED_PATHS = {
    "protocol_config": "configs/fgc/fgc-2-sf1-protocol-v5.toml",
    "protocol_freeze_result": "results/fgc-1-pro5-frz1.json",
    "diagnosis_result": "results/fgc-1-cal1-pref3.json",
    "static_input_result": "results/fgc-1-id2-all1.json",
    "successor_monitor_result": "results/fgc-1-hlt2-mon2.json",
    "run_plan_config": "configs/fgc/fgc-1-cal2-run1.toml",
}
EXPECTED_SCOPE = {
    "target_protocol": "FGC-2-SF1-PROTO5",
    "runtime_owner": ARTIFACT_ID,
    "calibration_branch": "GR-0",
    "role": "pre_trajectory_runtime_compositor_and_fresh_input_freeze",
    "fresh_GR0_trajectory_read": False,
    "SGBL_trajectory_read": False,
    "FGCQR_trajectory_read": False,
    "collapse_or_trapped_outcome_classified": False,
    "candidate_action_health_values_invented_for_GR0": False,
}
EXPECTED_LINEAGE = {
    "checkpoint_commit": "cf48cf729ff1c3a6effce6ecd3fdbeaf8d9645bf",
    "protocol_config_sha256": "516061e19ac14857eda26076fa18d9163db2fcc7c9d31aa6aebefd698377c9ab",
    "protocol_freeze_result_sha256": "7810557ed2eb82d9d484b0aeb65b75bbe3fd9565a4a52c5305c07144614add16",
    "diagnosis_result_sha256": "3ea43baf0fb9830b348fc49a79af5c535d8f7089aa03a957a22f6927528a4395",
    "static_input_result_sha256": "10e2661be521087467b970f1cd877311bd63a0cb6874a1ee745ca00ec74ab981",
    "successor_monitor_result_sha256": "61a9843a7fe7374c74373fea5580bd83e558d42faff3dd724f01e12b2411061a",
    "calibration_runtime_source_sha256": "e05c91116a64587dc763145117bcf6c35289f7c1a0f9e4eb7ad5c697aecef5f9",
    "checkpoint_must_be_ancestor_of_HEAD": True,
    "all_predecessor_blobs_must_be_read_from_checkpoint": True,
}
EXPECTED_RUNTIME = {
    "universal_runtime_stop_ids": list(GR0_UNIVERSAL_RUNTIME_STOP_IDS),
    "candidate_action_stop_count_for_GR0": 0,
    "replaced_numerical_decisions_are_common_event_constraints_and_weighted_spectra": True,
    "source_residual_and_iterative_refinement_map_to_historical_solver_stop_names": True,
    "hat_Lorentzian_margin_is_computed_from_actual_ADM_stage_metric": True,
    "causal_ledger_consumes_every_internal_stage_and_candidate_endpoint": True,
    "CFL_failure_is_retryable_and_not_a_scientific_stop": True,
    "scientific_failure_latches_first_premise_and_rolls_back_stage_and_causal_acceptance": True,
    "common_event_constraints_use_raw_magnitude_guards_and_roundoff_only_for_zero_classification": True,
    "spatial_spectrum_uses_proper_radial_distance_and_two_edge_compact_window": True,
    "temporal_spectrum_uses_causal_past_normal_flow_tracers_and_uniform_proper_time_resampling": True,
    "temporal_resampling_error_is_public_and_nonnegative": True,
    "trapped_sign_uses_both_metric_null_expansions_with_fixed_orientation": True,
    "trapped_sign_error_components_are_additive_only": True,
    "constraint_and_spectral_failures_are_vetoes_not_arbitrary_expansion_error_conversions": True,
}
EXPECTED_INPUT_FREEZE = {
    "eligible_amplitudes_in_order": ["5/2", "3"],
    "methods_in_order": ["RK4", "SSPRK3"],
    "resolutions_in_order": [1025, 2049, 4097],
    "expected_run_input_count": 12,
    "continuum_state_hashes_must_match_ID2": True,
    "physical_u_and_p_must_be_bitwise_preserved": True,
    "only_q_may_change_under_native_SBP_projection": True,
    "projected_state_and_expanded_run_config_hashes_must_be_serialized": True,
    "both_amplitudes_and_methods_must_pass_the_t0_common_event": True,
    "first_eligible_candidate_stops_later_candidates": True,
}
EXPECTED_NAMESPACE = {
    "calibration_output_root": "runs/fgc-2-sf1/proto5/calibration",
    "holdout_output_root": "runs/fgc-2-sf1/proto5/holdout",
    "both_roots_must_be_absent_or_empty_at_authorization": True,
    "authorization_must_not_create_either_root": True,
    "runner_must_refuse_overwrite": True,
}
EXPECTED_PROOF_KEYS = {
    "PROTO5_and_PRO5_FRZ1_must_validate_and_pass_from_immutable_checkpoint",
    "CAL1_ID2_and_HLT2_must_be_canonical_and_hash_matched_from_immutable_checkpoint",
    "runtime_stop_partition_must_match_PROTO5_inheritance_exactly",
    "candidate_only_health_fields_must_be_absent_from_GR0_runtime_API",
    "synthetic_runtime_accept_failure_rollback_CFL_and_boundary_controls_must_pass",
    "spatial_and_temporal_weighted_spectral_controls_must_pass_and_fail_when_injected",
    "all_twelve_projected_inputs_must_reconstruct_deterministically",
    "all_four_amplitude_method_t0_common_events_must_pass",
    "run_plan_must_bind_every_protocol_value_and_output_rule",
    "authorized_runner_must_be_hash_bound_and_refuse_drift_overwrite_or_dirty_tracked_state",
    "atomic_checkpoint_must_bind_members_prior_amplitude_records_and_exact_event_log_bytes",
    "resume_must_archive_uncommitted_log_tail_and_reject_divergent_history",
    "fresh_namespaces_must_be_observed_empty_without_creation",
    "canonical_hash_bound_result_required",
    "no_trajectory_may_be_read_or_advanced",
}
EXPECTED_CLAIMS = {
    "PROTO5_successor_runtime_compositor_implemented": True,
    "PROTO5_fresh_GR0_dynamic_calibration_authorized": True,
    "PROTO5_resolved_holdout_manifest_authorized": False,
    "classical_spherical_diagnostic_authorized": False,
    "FGCQR_holdout_execution_authorized": False,
    "retained_EFT_evolution_authorized": False,
    "physical_transition_claim_authorized": False,
    "FGCQR_mechanism_rejected": False,
    "general_gradient_route_rejected": False,
    "singularity_resolution_derived": False,
    "child_domain_or_topology_derived": False,
    "dark_sector_mechanism_derived": False,
    "varying_locally_measured_c_derived": False,
}
IMMUTABLE_PATHS = {
    "protocol_config": EXPECTED_PATHS["protocol_config"],
    "protocol_freeze_result": EXPECTED_PATHS["protocol_freeze_result"],
    "diagnosis_result": EXPECTED_PATHS["diagnosis_result"],
    "static_input_result": EXPECTED_PATHS["static_input_result"],
    "successor_monitor_result": EXPECTED_PATHS["successor_monitor_result"],
    "calibration_runtime_source": "src/recursive_horizons/fgc/evolution/calibration_runtime.py",
    "physical_protocol": "configs/fgc/fgc-2-sf1-protocol-v3.toml",
}
IMPLEMENTATION = (
    REPOSITORY / "src/recursive_horizons/fgc/evolution/calibration_runtime.py",
    REPOSITORY / "src/recursive_horizons/fgc/evolution/proto5_runtime.py",
    REPOSITORY / "scripts/reproduce_fgc_hlt3_mon3.py",
    REPOSITORY / "scripts/run_fgc_gr0_calibration.py",
)


def _serial(value: Any) -> Any:
    if isinstance(value, Fraction):
        return (
            str(value.numerator)
            if value.denominator == 1
            else f"{value.numerator}/{value.denominator}"
        )
    if is_dataclass(value):
        return _serial(asdict(value))
    if isinstance(value, Mapping):
        return {str(key): _serial(item) for key, item in value.items()}
    if isinstance(value, np.ndarray):
        return _serial(value.tolist())
    if isinstance(value, (tuple, list)):
        return [_serial(item) for item in value]
    if isinstance(value, np.generic):
        return _serial(value.item())
    if isinstance(value, float):
        if not isfinite(value):
            raise ValueError("nonfinite value cannot enter canonical evidence")
        return value
    if value is None or isinstance(value, (str, bool, int)):
        return value
    raise TypeError(f"unsupported evidence type {type(value).__name__}")


def _canonical(value: Any) -> str:
    return json.dumps(
        _serial(value),
        indent=2,
        sort_keys=True,
        ensure_ascii=True,
        allow_nan=False,
    ) + "\n"


def _pairs(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    answer: dict[str, Any] = {}
    for key, value in pairs:
        if key in answer:
            raise ValueError(f"duplicate JSON key: {key}")
        answer[key] = value
    return answer


def _digest(value: Any) -> str:
    return sha256(
        json.dumps(
            _serial(value),
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
            allow_nan=False,
        ).encode("utf-8")
    ).hexdigest()


def _sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _rel(path: Path) -> str:
    return path.resolve().relative_to(REPOSITORY).as_posix()


def _path(name: str, value: Any, expected: str) -> Path:
    if value != expected:
        raise ValueError(f"{name} must name {expected}")
    candidate = (REPOSITORY / expected).resolve()
    if not candidate.is_file() or _rel(candidate) != expected:
        raise ValueError(f"{name} must be a canonical existing repository file")
    return candidate


def _load_json_bytes(source: bytes, name: str) -> dict[str, Any]:
    try:
        value = json.loads(source.decode("utf-8"), object_pairs_hook=_pairs)
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as exc:
        raise ValueError(f"{name} is not valid unique-key JSON") from exc
    if not isinstance(value, dict) or source.decode("utf-8") != _canonical(value):
        raise ValueError(f"{name} is not canonical sorted JSON")
    return value


def _load_toml_bytes(source: bytes, name: str) -> dict[str, Any]:
    try:
        value = tomllib.loads(source.decode("utf-8"))
    except (UnicodeDecodeError, tomllib.TOMLDecodeError) as exc:
        raise ValueError(f"{name} is not valid UTF-8 TOML") from exc
    if not isinstance(value, dict):
        raise ValueError(f"{name} root must be a table")
    return value


def _fraction(name: str, value: Any, *, positive: bool = False) -> Fraction:
    if not isinstance(value, str):
        raise TypeError(f"{name} must be a rational string")
    answer = Q(value)
    if positive and answer <= 0:
        raise ValueError(f"{name} must be positive")
    return answer


def _git_blob(commit: str, relative: str, name: str) -> bytes:
    completed = subprocess.run(
        ["git", "show", f"{commit}:{relative}"],
        cwd=REPOSITORY,
        check=False,
        capture_output=True,
    )
    if completed.returncode != 0:
        raise ValueError(f"cannot read immutable {name} from {commit}")
    return completed.stdout


def _state_sha256(branch: str, point_count: int, state: EvolutionState) -> str:
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


def _validate_run_plan(path: Path) -> dict[str, Any]:
    with path.open("rb") as handle:
        raw = tomllib.load(handle)
    expected_root = {
        "schema_version",
        "artifact_id",
        "project_version",
        "metric_signature",
        "riemann_convention",
        "protocol_config",
        "protocol_freeze_result",
        "runtime_authorization_result",
        "static_input_result",
        "scope",
        "candidate_selection",
        "physical_inputs",
        "primary_method",
        "comparator_method",
        "numerics",
        "universal_thresholds",
        "trapped_sign_error",
        "provenance",
        "claims",
    }
    if set(raw) != expected_root:
        raise ValueError("CAL2 run-plan root keys differ")
    if (
        raw["schema_version"] != 1
        or raw["artifact_id"] != "FGC-1-CAL2-RUN1-PLAN"
        or raw["project_version"] != PROJECT_VERSION
        or raw["metric_signature"] != "-+++"
        or raw["riemann_convention"] != "plus_partial_mu_gamma_nu"
    ):
        raise ValueError("CAL2 run-plan identity differs")
    for name in (
        "protocol_config",
        "protocol_freeze_result",
        "static_input_result",
    ):
        if raw[name] != EXPECTED_PATHS[name]:
            raise ValueError(f"CAL2 {name} differs")
    if raw["runtime_authorization_result"] != "results/fgc-1-hlt3-mon3.json":
        raise ValueError("CAL2 runtime authorization path differs")
    if raw["scope"] != {
        "target_protocol": SF1_PROTOCOL_V5_ARTIFACT_ID,
        "branch": "GR-0",
        "role": "fresh_outcome_neutral_dynamic_amplitude_calibration",
        "physical_equations": "unredefined_GR0_specialization_of_REF1",
        "candidate_action_fields_or_health_stops_forbidden": True,
        "retained_EFT_interpretation": False,
    }:
        raise ValueError("CAL2 scope differs")
    selection = raw["candidate_selection"]
    if (
        selection["ordered_amplitudes"] != ["5/2", "3"]
        or selection["minimum_consecutive_trapped_common_events"] != 8
        or selection["trapped_sign_margin_over_combined_error_factor"] != 4
        or selection["primary_and_comparator_confirmation_required"] is not True
        or selection["all_three_nested_resolutions_required"] is not True
        or selection["first_eligible_candidate_stops_later_candidates"] is not True
    ):
        raise ValueError("CAL2 candidate-selection contract differs")
    physical = raw["physical_inputs"]
    expected_physical = {
        "length_unit_L0": "4",
        "cutoff_Lambda": "16",
        "planck_mass": "2",
        "scalar_mass": "3",
        "quartic_coupling": "1/2",
        "chi_center": "12",
        "chi_half_width": "2",
        "phi_seed_amplitude": "1/131072",
        "outer_radius": "128",
        "measurement_radius_maximum": "24",
    }
    if physical != expected_physical:
        raise ValueError("CAL2 physical inputs differ")
    methods = (raw["primary_method"], raw["comparator_method"])
    if methods != (
        {
            "method_label": "RK4",
            "integrator_id": PRIMARY_METHOD,
            "constraint_solve_method": "RK4",
            "spatial_order": 4,
            "coarsest_constraint_guard": "1/50",
            "finest_constraint_guard": "1/1000",
        },
        {
            "method_label": "SSPRK3",
            "integrator_id": COMPARATOR_METHOD,
            "constraint_solve_method": "SSPRK3",
            "spatial_order": 2,
            "coarsest_constraint_guard": "1/10",
            "finest_constraint_guard": "1/200",
        },
    ):
        raise ValueError("CAL2 method contract differs")
    numerics = raw["numerics"]
    required_numerics = {
        "resolutions": [1025, 2049, 4097],
        "cfl_maximum": "1/8",
        "retry_factor": "1/2",
        "maximum_CFL_retries_per_step": 32,
        "minimum_step_size": "1/1073741824",
        "ko_dissipation": "1/64",
        "fixed_outer_rows": 4,
        "final_coordinate_time": "32",
        "common_output_interval": "1/16",
        "checkpoint_interval": "1/4",
        "minimum_constraint_finest_pair_order": "3/2",
        "proper_spectral_taper_fraction": "1/8",
        "temporal_history_minimum_samples": 64,
        "normal_flow_tracer_radius_minimum": "1/2",
        "normal_flow_tracer_radius_maximum": "24",
        "normal_flow_tracer_spacing": "1/2",
        "normal_flow_tracer_integrator": "explicit_trapezoid_with_Euler_predictor_on_dr_dt_minus_shift_and_dtau_dt_lapse",
        "monitor_every_proposed_Runge_Kutta_stage": True,
        "common_event_constraints_and_spatial_spectra": True,
        "causal_past_temporal_spectra_after_minimum_samples": True,
        "common_event_endpoint_is_bitwise_target_time": True,
        "no_excision_in_measurement_region": True,
    }
    if numerics != required_numerics:
        raise ValueError("CAL2 numerical contract differs")
    thresholds = raw["universal_thresholds"]
    if thresholds != {
        "source_residual_infinity_max": "1/1000000000000",
        "source_iteration_max": 16,
        "source_residual_must_decrease_monotonically": True,
        "kinetic_condition_number_max": "10000000000",
        "hat_normal_factor": "9",
        "minimum_boundary_causal_buffer": "16",
    }:
        raise ValueError("CAL2 universal thresholds differ")
    if raw["provenance"] != {
        "output_root": "runs/fgc-2-sf1/proto5/calibration",
        "campaign_manifest": "runs/fgc-2-sf1/proto5/calibration/manifest.json",
        "append_only_event_log_name": "events.jsonl",
        "checkpoint_extension": ".npz",
        "fresh_output_root_required_at_campaign_start": True,
        "runner_refuses_overwrite": True,
        "resume_requires_matching_plan_and_manifest_hashes": True,
        "earlier_protocol_outputs_forbidden_as_evidence": True,
        "canonical_reduced_result_is_tracked_only_after_campaign_classification": True,
    }:
        raise ValueError("CAL2 provenance differs")
    if raw["claims"] != {key: False for key in EXPECTED_CLAIMS if key not in {
        "PROTO5_successor_runtime_compositor_implemented",
        "PROTO5_fresh_GR0_dynamic_calibration_authorized",
        "PROTO5_resolved_holdout_manifest_authorized",
    }}:
        # The run-plan omits protocol-sequencing booleans and carries only
        # downstream physical nonclaims.
        raise ValueError("CAL2 physical nonclaims differ")
    return raw


def load_config(path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
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
        "immutable_lineage",
        "runtime_composition",
        "input_freeze",
        "namespace",
        "proof_contract",
        "claims",
    }
    if set(raw) != expected_root:
        raise ValueError("HLT3 root keys differ")
    if (
        raw["schema_version"] != 1
        or raw["artifact_id"] != ARTIFACT_ID
        or raw["project_version"] != PROJECT_VERSION
        or raw["metric_signature"] != "-+++"
        or raw["riemann_convention"] != "plus_partial_mu_gamma_nu"
    ):
        raise ValueError("HLT3 identity or convention differs")
    paths = {
        name: _path(name, raw[name], expected)
        for name, expected in EXPECTED_PATHS.items()
    }
    if raw["scope"] != EXPECTED_SCOPE:
        raise ValueError("HLT3 scope differs")
    if raw["immutable_lineage"] != EXPECTED_LINEAGE:
        raise ValueError("HLT3 immutable lineage differs")
    if raw["runtime_composition"] != EXPECTED_RUNTIME:
        raise ValueError("HLT3 runtime composition differs")
    if raw["input_freeze"] != EXPECTED_INPUT_FREEZE:
        raise ValueError("HLT3 input freeze differs")
    if raw["namespace"] != EXPECTED_NAMESPACE:
        raise ValueError("HLT3 namespace contract differs")
    if set(raw["proof_contract"]) != EXPECTED_PROOF_KEYS or any(
        value is not True for value in raw["proof_contract"].values()
    ):
        raise ValueError("HLT3 proof contract differs")
    if raw["claims"] != EXPECTED_CLAIMS:
        raise ValueError("HLT3 claims differ")
    run_plan = _validate_run_plan(paths["run_plan_config"])
    return {"raw": raw, "paths": paths, "run_plan": run_plan}


def _immutable_lineage(loaded: Mapping[str, Any]) -> dict[str, Any]:
    lineage = loaded["raw"]["immutable_lineage"]
    commit = lineage["checkpoint_commit"]
    ancestor = subprocess.run(
        ["git", "merge-base", "--is-ancestor", commit, "HEAD"],
        cwd=REPOSITORY,
        check=False,
        capture_output=True,
    )
    if ancestor.returncode != 0:
        raise ValueError("HLT3 checkpoint is not an ancestor of HEAD")
    blobs = {
        name: _git_blob(commit, relative, name)
        for name, relative in IMMUTABLE_PATHS.items()
    }
    hash_keys = {
        "protocol_config": "protocol_config_sha256",
        "protocol_freeze_result": "protocol_freeze_result_sha256",
        "diagnosis_result": "diagnosis_result_sha256",
        "static_input_result": "static_input_result_sha256",
        "successor_monitor_result": "successor_monitor_result_sha256",
        "calibration_runtime_source": "calibration_runtime_source_sha256",
    }
    for name, key in hash_keys.items():
        if sha256(blobs[name]).hexdigest() != lineage[key]:
            raise ValueError(f"immutable {name} hash differs")

    protocol_raw = _load_toml_bytes(blobs["protocol_config"], "immutable PROTO5")
    protocol = validate_sf1_protocol_v5(protocol_raw)
    if protocol["artifact_id"] != SF1_PROTOCOL_V5_ARTIFACT_ID:
        raise ValueError("immutable PROTO5 identity differs")
    freeze = _load_json_bytes(blobs["protocol_freeze_result"], "immutable PRO5-FRZ1")
    diagnosis = _load_json_bytes(blobs["diagnosis_result"], "immutable CAL1")
    static = _load_json_bytes(blobs["static_input_result"], "immutable ID2")
    hlt2 = _load_json_bytes(blobs["successor_monitor_result"], "immutable HLT2")
    if freeze.get("gate_status", {}).get("PROTO5_outcome_neutral_protocol_frozen") is not True:
        raise ValueError("immutable PRO5-FRZ1 gate is not true")
    if diagnosis.get("gate_status", {}).get("PROTO5_premise_revision_required") is not True:
        raise ValueError("immutable CAL1 diagnosis gate differs")
    if static.get("gate_status", {}).get("all_case_static_ledger_completed") is not True:
        raise ValueError("immutable ID2 gate differs")
    if hlt2.get("gate_status", {}).get("PROTO4_successor_admission_monitor_implemented") is not True:
        raise ValueError("immutable HLT2 gate differs")
    eligible = static.get("artifact_payload", {}).get(
        "static_decision_summary", {}
    ).get("eligible_amplitudes_in_frozen_order")
    if eligible != ["5/2", "3"]:
        raise ValueError("immutable ID2 eligible amplitude order differs")
    physical_protocol = _load_toml_bytes(
        blobs["physical_protocol"], "immutable physical PROTO3"
    )
    return {
        "checkpoint_commit": commit,
        "checkpoint_is_ancestor_of_HEAD": True,
        "protocol": protocol,
        "protocol_blob_sha256": sha256(blobs["protocol_config"]).hexdigest(),
        "freeze_artifact_id": freeze["artifact_id"],
        "freeze_gate_passed": True,
        "diagnosis_artifact_id": diagnosis["artifact_id"],
        "diagnosis_is_pre_trajectory": diagnosis.get("artifact_payload", {})
        .get("epistemic_boundary", {})
        .get("trajectory_read")
        is False,
        "static_input_artifact_id": static["artifact_id"],
        "eligible_amplitudes": eligible,
        "successor_monitor_artifact_id": hlt2["artifact_id"],
        "HLT2_monitor_gate_passed": True,
        "static_result": static,
        "physical_protocol": physical_protocol,
    }


def _id2_method_record(
    static_result: Mapping[str, Any],
    amplitude: str,
    method: str,
) -> Mapping[str, Any]:
    ledger = static_result["artifact_payload"]["GR0_candidate_ledger"]
    candidate = next(item for item in ledger if item["amplitude"] == amplitude)
    if candidate["GR0_static_input_passed"] is not True:
        raise ValueError(f"ID2 amplitude {amplitude} is not statically eligible")
    return next(item for item in candidate["methods"] if item["method"] == method)


def _freeze_inputs(
    loaded: Mapping[str, Any],
    lineage: Mapping[str, Any],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    run_plan = loaded["run_plan"]
    physical = run_plan["physical_inputs"]
    numerics = run_plan["numerics"]
    plan_hash = _sha(loaded["paths"]["run_plan_config"])
    records: list[dict[str, Any]] = []
    common_events: list[dict[str, Any]] = []
    for amplitude in run_plan["candidate_selection"]["ordered_amplitudes"]:
        method_states: dict[str, list[EvolutionState]] = {}
        method_grids: dict[str, list[Any]] = {}
        for method_table in (
            run_plan["primary_method"],
            run_plan["comparator_method"],
        ):
            method = method_table["method_label"]
            order = method_table["spatial_order"]
            id2_method = _id2_method_record(
                lineage["static_result"], amplitude, method
            )
            states: list[EvolutionState] = []
            grids: list[Any] = []
            for point_count, id2_resolution in zip(
                numerics["resolutions"],
                id2_method["resolutions"],
                strict=True,
            ):
                if id2_resolution["point_count"] != point_count:
                    raise ValueError("ID2 resolution order differs from CAL2")
                parameters = PulseParameters(
                    chi_amplitude=float(Q(amplitude)),
                    center=float(Q(physical["chi_center"])),
                    half_width=float(Q(physical["chi_half_width"])),
                    phi_amplitude=float(Q(physical["phi_seed_amplitude"])),
                    planck_mass=float(Q(physical["planck_mass"])),
                    scalar_mass=float(Q(physical["scalar_mass"])),
                    quartic_coupling=float(Q(physical["quartic_coupling"])),
                )
                initial = construct_gr0_grid_initial_data(
                    parameters,
                    point_count=point_count,
                    outer_radius=float(Q(physical["outer_radius"])),
                    constraint_method=method,
                    diagnostic_spatial_order=order,
                )
                continuum_hash = _state_sha256("GR-0", point_count, initial.state)
                if continuum_hash != id2_resolution["state_sha256"]:
                    raise ValueError("reconstructed continuum state differs from ID2")
                projected = project_gr0_semidiscrete_state(
                    initial,
                    spatial_order=order,
                )
                if not np.array_equal(projected.u, initial.state.u) or not np.array_equal(
                    projected.p, initial.state.p
                ):
                    raise ValueError("PROTO5 projection changed physical u or p")
                projected_hash = array_content_sha256(
                    projected.u,
                    projected.p,
                    projected.q,
                )
                payload = {
                    "plan_sha256": plan_hash,
                    "protocol_sha256": lineage["protocol_blob_sha256"],
                    "branch": "GR-0",
                    "amplitude": amplitude,
                    "method": method,
                    "integrator_id": method_table["integrator_id"],
                    "spatial_order": order,
                    "point_count": point_count,
                    "grid_spacing": initial.grid.spacing,
                    "continuum_ID2_state_sha256": continuum_hash,
                    "projected_state_sha256": projected_hash,
                    "cfl_maximum": numerics["cfl_maximum"],
                    "ko_dissipation": numerics["ko_dissipation"],
                    "common_output_interval": numerics["common_output_interval"],
                    "final_coordinate_time": numerics["final_coordinate_time"],
                }
                records.append(
                    {
                        **payload,
                        "expanded_run_config_sha256": _digest(payload),
                        "u_bitwise_preserved": True,
                        "p_bitwise_preserved": True,
                        "q_initialized_by_native_SBP_operator": True,
                        "trajectory_advanced": False,
                    }
                )
                states.append(projected)
                grids.append(initial.grid)
            assessed = proto5_gr0_common_event(
                states,
                grids,
                method=method,
                coordinate_time=0.0,
                cutoff=float(Q(physical["cutoff_Lambda"])),
                measurement_radius_maximum=float(
                    Q(physical["measurement_radius_maximum"])
                ),
                taper_fraction=float(
                    Q(numerics["proper_spectral_taper_fraction"])
                ),
            )
            if not assessed.admission_passed:
                raise ValueError(f"{amplitude} {method} failed its frozen t0 event")
            common_events.append(
                {
                    "amplitude": amplitude,
                    "method": method,
                    "coordinate_time": 0.0,
                    "constraint_admission": assessed.constraint_admission,
                    "spatial_spectral_admission": assessed.spatial_spectral_admission,
                    "maximum_spatial_round_trip_interpolation_infinity": (
                        assessed.maximum_spatial_round_trip_interpolation_infinity
                    ),
                    "admission_passed": True,
                    "trajectory_advanced": False,
                }
            )
            method_states[method] = states
            method_grids[method] = grids
        if tuple(method_states) != ("RK4", "SSPRK3"):
            raise RuntimeError("HLT3 method reconstruction order differs")
    if len(records) != 12 or len(common_events) != 4:
        raise RuntimeError("HLT3 did not freeze twelve inputs and four t0 events")
    return records, common_events


def _synthetic_state() -> EvolutionState:
    radii = np.linspace(0.0, 1.0, 9)
    u = np.zeros((9, 6))
    u[:, 0] = 1.0
    u[:, 2] = 1.0
    u[:, 3] = radii
    p = np.zeros_like(u)
    q = np.zeros_like(u)
    q[:, 3] = 1.0
    return EvolutionState(u, p, q)


def _synthetic_rhs(*, bad_time: float | None = None):
    def rhs(time: float, state: EvolutionState) -> EvolutionRHS:
        zero = np.zeros(state.shape)
        return EvolutionRHS(
            zero,
            zero,
            zero,
            {
                "source_residual_infinity": 0.0,
                "source_refinement_iterations": 0,
                "source_residual_decreased_monotonically": True,
                "kinetic_condition_infinity": 1.0,
                "coordinate_speed_upper": 1.0,
                "minimum_lapse": (
                    0.0 if bad_time is not None and time >= bad_time else 1.0
                ),
                "minimum_radial_metric": 1.0,
                "minimum_areal_radius_away_from_center": 0.125,
            },
        )

    return rhs


def _transaction(
    *, boundary_geometry: BoundaryGeometry | None = None
) -> GR0RuntimeStageTransaction:
    return GR0RuntimeStageTransaction(
        thresholds=GR0UniversalThresholds(),
        causal_state=CausalBudgetState(previous_speed_upper=1.0),
        boundary_geometry=(
            BoundaryGeometry(128.0, 24.0, 16.0, 0.375)
            if boundary_geometry is None
            else boundary_geometry
        ),
        grid_spacing=0.125,
        cfl_maximum=0.125,
    )


def _runtime_controls() -> dict[str, Any]:
    state = _synthetic_state()
    healthy_proposal = propose_step(
        method=PRIMARY_METHOD,
        time=0.0,
        step_size=0.01,
        state=state,
        rhs=_synthetic_rhs(),
    )
    healthy = _transaction()
    receipt = dict(healthy(healthy_proposal))

    failure_proposal = propose_step(
        method=PRIMARY_METHOD,
        time=0.0,
        step_size=0.01,
        state=state,
        rhs=_synthetic_rhs(bad_time=0.01),
    )
    failed = _transaction()
    causal_before = failed.causal_state
    try:
        failed(failure_proposal)
    except GR0RuntimeStop as stop:
        failure_reason = stop.reason
    else:
        raise RuntimeError("HLT3 synthetic scientific failure was accepted")
    rollback = (
        failed.state.accepted_stage_count == 0
        and failed.causal_state == causal_before
    )
    try:
        failed(failure_proposal)
    except GR0RuntimeStop as repeated:
        immutable = repeated.reason == failure_reason == "nonpositive_lapse"
    else:
        immutable = False

    oversized = propose_step(
        method=COMPARATOR_METHOD,
        time=0.0,
        step_size=0.02,
        state=state,
        rhs=_synthetic_rhs(),
    )
    cfl = _transaction()
    try:
        cfl(oversized)
    except CFLRetryRequired as retry:
        cfl_record = {
            "retry_required": True,
            "observed_ratio": retry.observed_ratio,
            "maximum_ratio": retry.maximum_ratio,
            "scientific_failure_latched": cfl.state.first_failed_premise is not None,
            "causal_or_stage_acceptance_advanced": (
                cfl.state.accepted_stage_count != 0
                or cfl.causal_state.accepted_time != 0.0
            ),
        }
    else:
        raise RuntimeError("HLT3 oversized CFL proposal was accepted")

    boundary = _transaction(
        boundary_geometry=BoundaryGeometry(0.5, 0.1, 0.35, 0.05)
    )
    try:
        boundary(healthy_proposal)
    except GR0RuntimeStop as stop:
        boundary_reason = stop.reason
    else:
        raise RuntimeError("HLT3 boundary-exhaustion proposal was accepted")
    boundary_rollback = (
        boundary_reason == "boundary_causal_buffer"
        and boundary.state.accepted_stage_count == 0
        and boundary.causal_state.accepted_time == 0.0
    )
    if not rollback or not immutable or not boundary_rollback:
        raise RuntimeError("HLT3 transaction rollback or first-failure control failed")
    checkpoint_log = b'{"event_index":0}\n'
    uncommitted_tail = b'{"event_index":1}\n'
    recovery_plans = {
        "absent": proto5_event_log_recovery_plan(None, checkpoint_log),
        "matching": proto5_event_log_recovery_plan(
            checkpoint_log,
            checkpoint_log,
        ),
        "behind": proto5_event_log_recovery_plan(b"", checkpoint_log),
        "ahead": proto5_event_log_recovery_plan(
            checkpoint_log + uncommitted_tail,
            checkpoint_log,
        ),
    }
    expected_recovery_actions = {
        "absent": "restored_missing_log_from_checkpoint",
        "matching": "event_log_already_matches_checkpoint",
        "behind": "completed_checkpoint_owned_log_write",
        "ahead": "archived_uncommitted_tail_and_restored_checkpoint",
    }
    if any(
        recovery_plans[name].action != action
        or recovery_plans[name].restored_payload != checkpoint_log
        for name, action in expected_recovery_actions.items()
    ) or recovery_plans["ahead"].orphaned_tail != uncommitted_tail:
        raise RuntimeError("HLT3 event-log recovery controls failed")
    try:
        proto5_event_log_recovery_plan(
            b'{"event_index":99}\n',
            checkpoint_log,
        )
    except ValueError as error:
        divergent_rejected = "diverges" in str(error)
    else:
        divergent_rejected = False
    if not divergent_rejected:
        raise RuntimeError("HLT3 divergent event-log history was accepted")
    return {
        "universal_stop_ids": list(GR0_RUNTIME_STOP_PRIORITY),
        "candidate_action_stop_count": 0,
        "healthy_receipt": receipt,
        "late_scientific_failure_reason": failure_reason,
        "late_failure_rolled_back": rollback,
        "first_scientific_failure_immutable": immutable,
        "CFL_retry_control": cfl_record,
        "boundary_failure_reason": boundary_reason,
        "boundary_failure_rolled_back": boundary_rollback,
        "event_log_recovery_controls": {
            "actions": {
                name: plan.action for name, plan in recovery_plans.items()
            },
            "ahead_orphaned_tail_sha256": sha256(uncommitted_tail).hexdigest(),
            "divergent_history_rejected": divergent_rejected,
            "all_recovery_controls_passed": True,
        },
        "all_controls_passed": True,
    }


def _temporal_controls() -> dict[str, Any]:
    sample_count = 64
    tracer_count = 2
    times = np.repeat(
        np.linspace(0.0, 1.0, sample_count)[:, None], tracer_count, axis=1
    )
    zero = np.zeros(
        (sample_count, tracer_count, len(PROTO4_SPECTRAL_FIELD_ORDER))
    )
    admitted = proto5_temporal_spectral_admission(
        point_counts=(1025, 2049, 4097),
        proper_times_by_resolution=(times, times, times),
        field_histories_by_resolution=(zero, zero, zero),
    )
    alias = zero.copy()
    alias[:, :, 0] = (-1.0) ** np.arange(sample_count)[:, None]
    rejected = proto5_temporal_spectral_admission(
        point_counts=(1025, 2049, 4097),
        proper_times_by_resolution=(times, times, times),
        field_histories_by_resolution=(alias, alias, alias),
    )
    if not admitted.admission_passed or rejected.admission_passed:
        raise RuntimeError("HLT3 temporal spectral controls failed")
    return {
        "minimum_sample_count": 64,
        "zero_history_control": admitted,
        "Nyquist_alias_injection": rejected,
        "low_frequency_control_passed": True,
        "alias_injection_rejected": True,
    }


def _validate_historical_namespace_evidence(
    raw: Mapping[str, Any],
    evidence: Mapping[str, Any],
) -> dict[str, Any]:
    expected_paths = {
        "calibration_output_root": raw["namespace"]["calibration_output_root"],
        "holdout_output_root": raw["namespace"]["holdout_output_root"],
    }
    records = evidence.get("records")
    if not isinstance(records, list) or len(records) != 2:
        raise ValueError("historical HLT3 namespace evidence has wrong cardinality")
    by_name = {item.get("name"): item for item in records}
    if set(by_name) != set(expected_paths):
        raise ValueError("historical HLT3 namespace names differ")
    for name, expected_path in expected_paths.items():
        record = by_name[name]
        if (
            record.get("path") != expected_path
            or record.get("entry_count") != 0
            or record.get("empty") is not True
            or record.get("created_by_authorization") is not False
            or not isinstance(record.get("existed_before_authorization"), bool)
        ):
            raise ValueError("historical HLT3 namespace evidence differs")
    if (
        evidence.get("both_fresh_namespaces_empty") is not True
        or evidence.get("authorization_created_no_namespace") is not True
    ):
        raise ValueError("historical HLT3 namespace conclusion differs")
    return _serial(evidence)


def _namespace_precondition(raw: Mapping[str, Any]) -> dict[str, Any]:
    records = []
    for key in ("calibration_output_root", "holdout_output_root"):
        relative = raw["namespace"][key]
        path = REPOSITORY / relative
        existed = path.exists()
        if existed:
            if not path.is_dir():
                raise ValueError(f"{relative} exists and is not a directory")
            entries = sorted(item.name for item in path.iterdir())
        else:
            entries = []
        if entries:
            raise ValueError(f"HLT3 requires empty fresh namespace {relative}")
        records.append(
            {
                "name": key,
                "path": relative,
                "existed_before_authorization": existed,
                "entry_count": len(entries),
                "empty": True,
                "created_by_authorization": False,
            }
        )
    return {
        "records": records,
        "both_fresh_namespaces_empty": True,
        "authorization_created_no_namespace": True,
    }


def record(
    path: Path = DEFAULT_CONFIG,
    *,
    historical_namespace_evidence: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    loaded = load_config(path)
    raw = loaded["raw"]
    lineage = _immutable_lineage(loaded)
    if tuple(GR0_RUNTIME_STOP_PRIORITY) != tuple(GR0_UNIVERSAL_RUNTIME_STOP_IDS):
        raise RuntimeError("HLT3 runtime stop partition differs")
    if dict(PROTO5_METHOD_CONTRACT) != {
        "RK4": {
            "spatial_order": 4,
            "coarsest_constraint_guard": 1.0 / 50.0,
            "finest_constraint_guard": 1.0 / 1000.0,
        },
        "SSPRK3": {
            "spatial_order": 2,
            "coarsest_constraint_guard": 1.0 / 10.0,
            "finest_constraint_guard": 1.0 / 200.0,
        },
    }:
        raise RuntimeError("HLT3 method-owned guards differ from PROTO5")
    namespace = (
        _namespace_precondition(raw)
        if historical_namespace_evidence is None
        else _validate_historical_namespace_evidence(
            raw,
            historical_namespace_evidence,
        )
    )
    inputs, common_events = _freeze_inputs(loaded, lineage)
    runtime = _runtime_controls()
    temporal = _temporal_controls()
    plan_hash = _sha(loaded["paths"]["run_plan_config"])
    input_manifest_hash = _digest(inputs)
    campaign_id = _digest(
        {
            "artifact_id": "FGC-1-CAL2-RUN1-PLAN",
            "run_plan_sha256": plan_hash,
            "input_manifest_sha256": input_manifest_hash,
            "ordered_amplitudes": ["5/2", "3"],
            "output_root": raw["namespace"]["calibration_output_root"],
        }
    )
    gate_status = dict(raw["claims"])
    result = {
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "project_version": PROJECT_VERSION,
        "classification": "pre_trajectory_PROTO5_GR0_runtime_compositor_and_fresh_calibration_authorization",
        "generated_by": "scripts/reproduce_fgc_hlt3_mon3.py",
        "derivation_document": _rel(OWNER_DOCUMENT),
        "derivation_document_sha256": _sha(OWNER_DOCUMENT),
        "source_config_sha256": {
            _rel(path): _sha(path),
            _rel(loaded["paths"]["run_plan_config"]): plan_hash,
        },
        "predecessor_sha256": {
            EXPECTED_PATHS["protocol_config"]: raw["immutable_lineage"][
                "protocol_config_sha256"
            ],
            EXPECTED_PATHS["protocol_freeze_result"]: raw["immutable_lineage"][
                "protocol_freeze_result_sha256"
            ],
            EXPECTED_PATHS["diagnosis_result"]: raw["immutable_lineage"][
                "diagnosis_result_sha256"
            ],
            EXPECTED_PATHS["static_input_result"]: raw["immutable_lineage"][
                "static_input_result_sha256"
            ],
            EXPECTED_PATHS["successor_monitor_result"]: raw[
                "immutable_lineage"
            ]["successor_monitor_result_sha256"],
            "src/recursive_horizons/fgc/evolution/calibration_runtime.py": raw[
                "immutable_lineage"
            ]["calibration_runtime_source_sha256"],
        },
        "implementation_sha256": {
            _rel(item): _sha(item) for item in IMPLEMENTATION
        },
        "scope_bindings": dict(raw["scope"]),
        "gate_status": gate_status,
        "artifact_payload": {
            "immutable_lineage": {
                key: value
                for key, value in lineage.items()
                if key not in {"static_result", "physical_protocol"}
            },
            "runtime_composition": runtime,
            "temporal_spectral_controls": temporal,
            "frozen_run_plan": {
                "artifact_id": loaded["run_plan"]["artifact_id"],
                "run_plan_sha256": plan_hash,
                "campaign_id": campaign_id,
                "input_manifest_sha256": input_manifest_hash,
                "ordered_amplitudes": ["5/2", "3"],
                "first_eligible_candidate_stops_later_candidates": True,
                "expanded_input_count": len(inputs),
            },
            "frozen_run_inputs": inputs,
            "initial_common_events": common_events,
            "namespace_precondition": namespace,
            "decision": "HLT3_implements_the_PROTO5_GR0_runtime_contract_and_authorizes_only_the_fresh_calibration_campaign",
            "epistemic_boundary": {
                "fresh_GR0_trajectory_read": False,
                "fresh_GR0_trajectory_advanced": False,
                "synthetic_zero_dynamics_used_only_for_transaction_controls": True,
                "t0_physical_inputs_reconstructed": True,
                "trapped_or_collapse_outcome_classified": False,
                "SGBL_or_FGCQR_outcome_read": False,
                "fresh_GR0_calibration_is_now_authorized_but_not_completed": True,
                "resolved_holdout_manifest_exists": False,
                "mechanism_question_answered": False,
            },
        },
        "nonclaims": {
            "fresh_GR0_calibration_completed": False,
            "trapped_sphere_formed_dynamically": False,
            "SGBL_health_definition_supplied": False,
            "SGBL_comparison_completed": False,
            "FGCQR_trajectory_opened": False,
            "positive_Raychaudhuri_margin_measured": False,
            "constraint_to_DEF1_stability_map_supplied": False,
            "retained_EFT_validity": False,
            "physical_transition": False,
            "singularity_resolution": False,
            "child_domain_or_topology": False,
            "dark_sector_mechanism": False,
            "varying_locally_measured_c": False,
        },
    }
    return result


def load_canonical_result(path: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    return _load_json_bytes(path.read_bytes(), "HLT3 result")


def _verify_authorized_postlaunch_namespace(
    loaded: Mapping[str, Any],
    observed: Mapping[str, Any],
    result_path: Path,
) -> None:
    raw = loaded["raw"]
    historical = observed.get("artifact_payload", {}).get(
        "namespace_precondition", {}
    )
    _validate_historical_namespace_evidence(raw, historical)
    calibration_root = REPOSITORY / raw["namespace"]["calibration_output_root"]
    holdout_root = REPOSITORY / raw["namespace"]["holdout_output_root"]
    if not calibration_root.is_dir():
        raise ValueError("postlaunch HLT3 check requires its calibration directory")
    if holdout_root.exists():
        if not holdout_root.is_dir() or any(holdout_root.iterdir()):
            raise ValueError("HLT3 does not authorize a nonempty holdout namespace")
    manifest_path = calibration_root / "manifest.json"
    if not manifest_path.is_file():
        raise ValueError("postlaunch calibration manifest is absent")
    manifest = _load_json_bytes(manifest_path.read_bytes(), "CAL2 launch manifest")
    frozen = observed["artifact_payload"]["frozen_run_plan"]
    expected_manifest_fields = {
        "schema_version": 1,
        "runner_id": "FGC-1-CAL2-RUN1-RUNNER",
        "campaign_id": frozen["campaign_id"],
        "plan_path": EXPECTED_PATHS["run_plan_config"],
        "plan_sha256": frozen["run_plan_sha256"],
        "authorization_path": _rel(result_path),
        "authorization_sha256": _sha(result_path),
        "input_manifest_sha256": frozen["input_manifest_sha256"],
        "ordered_amplitudes": frozen["ordered_amplitudes"],
        "implementation_sha256": observed["implementation_sha256"],
        "outcome_fields_present_at_creation": False,
    }
    for key, value in expected_manifest_fields.items():
        if manifest.get(key) != value:
            raise ValueError(f"postlaunch CAL2 manifest field differs: {key}")
    if not isinstance(manifest.get("git_commit"), str):
        raise ValueError("postlaunch CAL2 manifest lacks its Git checkpoint")
    ancestor = subprocess.run(
        ["git", "merge-base", "--is-ancestor", manifest["git_commit"], "HEAD"],
        cwd=REPOSITORY,
        check=False,
        capture_output=True,
    )
    if ancestor.returncode != 0:
        raise ValueError("postlaunch CAL2 Git checkpoint is not an ancestor of HEAD")
    created = manifest.get("created_unix_time")
    if not isinstance(created, (int, float)) or not isfinite(float(created)):
        raise ValueError("postlaunch CAL2 manifest creation time is invalid")


def verify_canonical(path: Path, config_path: Path = DEFAULT_CONFIG) -> None:
    observed = load_canonical_result(path)
    loaded = load_config(config_path)
    calibration_root = (
        REPOSITORY
        / loaded["raw"]["namespace"]["calibration_output_root"]
    )
    if calibration_root.exists() and any(calibration_root.iterdir()):
        _verify_authorized_postlaunch_namespace(loaded, observed, path)
        expected = record(
            config_path,
            historical_namespace_evidence=observed["artifact_payload"][
                "namespace_precondition"
            ],
        )
    else:
        expected = record(config_path)
    if observed != _serial(expected):
        raise ValueError("HLT3 output differs from a fresh reproduction")


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
    result = record(config)
    output.write_text(_canonical(result), encoding="utf-8")
    print(f"wrote {output}")


if __name__ == "__main__":
    main()
