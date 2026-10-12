#!/usr/bin/env python3
"""Regenerate the pre-trajectory FGC-1-HLT8-MON8 authorization record."""

from __future__ import annotations

import argparse
from copy import deepcopy
from dataclasses import asdict, replace
from hashlib import sha256
from pathlib import Path
import subprocess
import sys
import tomllib
from typing import Any, Mapping

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-hlt8-mon8.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-hlt8-mon8.json"
OWNER_DOCUMENT = REPOSITORY / "docs/fgc-hlt8-mon8.md"

from scripts import reproduce_fgc_hlt4_mon4 as hlt4  # noqa: E402
from recursive_horizons.fgc.evolution.calibration_runtime import (  # noqa: E402
    project_gr0_semidiscrete_state,
)
from recursive_horizons.fgc.evolution.gr0_calibration import (  # noqa: E402
    construct_gr0_grid_initial_data,
)
from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    EvolutionState,
    SBPFirstDerivative,
    UniformRadialGrid,
    array_content_sha256,
)
from recursive_horizons.fgc.evolution.proto10_runtime import (  # noqa: E402
    Proto10GR0CommonEventAssessment,
    proto10_gr0_common_event,
)
from recursive_horizons.fgc.evolution.proto4_admission import (  # noqa: E402
    PROTO4_SPECTRAL_FIELD_ORDER,
    SpectralThresholds,
)
from recursive_horizons.fgc.evolution.proto7_runtime import (  # noqa: E402
    Proto7GR0EvolutionOperator,
)
from recursive_horizons.fgc.evolution.proto9_runtime import (  # noqa: E402
    PROTO9_POINT_COUNTS,
)
from recursive_horizons.fgc.evolution.protocol_v10 import (  # noqa: E402
    validate_sf1_protocol_v10,
)
from recursive_horizons.fgc.evolution.spectral_sensitivity import (  # noqa: E402
    resolved_or_saturated_admission,
)
from recursive_horizons.fgc.evolution.vectorized_source import (  # noqa: E402
    ADM_CENTER_PARITIES,
)
from recursive_horizons.fgc.initial_data_preflight import (  # noqa: E402
    PulseParameters,
)


Q = hlt4.Q
ARTIFACT_ID = "FGC-1-HLT8-MON8"
PROJECT_VERSION = "0.11.0"
CHECKPOINT_COMMIT = "3f7ac298a027209e5bdb2d59c8c82e04c25211ab"
EXPECTED_PATHS = {
    "protocol_config": "configs/fgc/fgc-2-sf1-protocol-v10.toml",
    "protocol_freeze_config": "configs/fgc/fgc-1-pro10-frz1.toml",
    "protocol_freeze_result": "results/fgc-1-pro10-frz1.json",
    "predecessor_runtime_result": "results/fgc-1-hlt7-mon7.json",
    "diagnosis_result": "results/fgc-1-cal6-pref8.json",
    "static_input_result": "results/fgc-1-id2-all1.json",
    "predecessor_run_plan": "configs/fgc/fgc-1-cal6-run1.toml",
    "run_plan_config": "configs/fgc/fgc-1-cal7-run1.toml",
}
EXPECTED_SCOPE = {
    "target_protocol": "FGC-2-SF1-PROTO10",
    "runtime_owner": ARTIFACT_ID,
    "calibration_branch": "GR-0",
    "role": "pre_trajectory_guarded_finest_pair_runtime_binding_and_fresh_input_freeze",
    "PROTO9_preflight_and_diagnosis_history_disclosed": True,
    "PROTO10_trajectory_read": False,
    "SGBL_trajectory_read": False,
    "FGCQR_trajectory_read": False,
    "collapse_or_trapped_outcome_classified": False,
    "diagnostic_saturation_promoted_to_physical_resolution": False,
}
EXPECTED_LINEAGE = {
    "checkpoint_commit": CHECKPOINT_COMMIT,
    "protocol_config_sha256": "f06ea93e82a4cdda2a0474623bc0db9ac679e4044657afb32be0d5342fe437fa",
    "protocol_freeze_config_sha256": "ed75ba9be8f9dbcb11a1a9e40a1ce3bcaa7e36688ac80b2ed2c13a3ab285b23a",
    "protocol_freeze_result_sha256": "f807eabf24e1bfec3dd9a30cacf489ac812356ddd60a08881d42563cd3e17ae3",
    "predecessor_runtime_result_sha256": "e6eb5832302eea39f7adf7fee2996cf659979ebbc58317b8f6d3b7e092f789c5",
    "diagnosis_result_sha256": "3d4624cca60363c9daf6f65f9393d959d7a7fafbcf6eac7f8dc9c2f237c24f5e",
    "static_input_result_sha256": "10e2661be521087467b970f1cd877311bd63a0cb6874a1ee745ca00ec74ab981",
    "predecessor_run_plan_sha256": "ed0daf5be8353274131ad20ea13534fed491aa238a4234f07648e2a3bd6cafb3",
    "checkpoint_must_be_ancestor_of_HEAD": True,
    "all_predecessor_blobs_must_be_read_from_checkpoint": True,
}
EXPECTED_RUNTIME = {
    "inherited_PROTO9_raw_common_event_compositor_byte_unchanged": True,
    "inherited_PROTO7_evolution_transaction_and_source_retry_byte_unchanged": True,
    "exact_point_counts_required": [2049, 4097, 8193],
    "raw_PROTO9_assessment_must_be_recomputed_and_serialized": True,
    "raw_ratio_and_absolute_budget_decisions_must_match_PROTO9_exactly": True,
    "coarse_to_medium_ratios_remain_direct_vetoes": True,
    "failed_finest_pair_requires_medium_and_fine_tail_erasure_map_witnesses": True,
    "failed_finest_pair_requires_two_pair_round_trip_contraction": True,
    "failed_finest_pair_requires_whole_profile_infinity_and_RMS_contraction": True,
    "round_trip_is_not_a_continuum_error_bound": True,
    "diagnostic_saturation_is_not_physical_resolution": True,
    "physical_inputs_equations_methods_CFL_retries_and_stops_unchanged": True,
    "runner_binds_and_restores_immutable_predecessor_engine": True,
}
EXPECTED_INPUT = {
    "eligible_amplitudes_in_order": ["5/2", "3"],
    "methods_in_order": ["RK4", "SSPRK3"],
    "resolutions_in_order": [2049, 4097, 8193],
    "expected_run_input_count": 12,
    "all_twelve_projected_state_hashes_must_match_HLT7": True,
    "physical_u_p_and_initial_profiles_unchanged": True,
    "q_must_be_initialized_by_native_SBP_operator": True,
    "expanded_PROTO10_run_config_hashes_must_be_serialized": True,
    "all_twelve_accepted_state_source_prechecks_must_pass": True,
    "both_amplitudes_and_methods_must_pass_the_PROTO10_t0_common_event": True,
    "all_raw_PROTO9_t0_failures_must_remain_public": True,
    "first_eligible_candidate_stops_later_candidates": True,
}
EXPECTED_NAMESPACE = {
    "calibration_output_root": "runs/fgc-2-sf1/proto10/calibration",
    "holdout_output_root": "runs/fgc-2-sf1/proto10/holdout",
    "both_roots_must_be_absent_or_empty_at_authorization": True,
    "authorization_must_not_create_either_root": True,
    "runner_must_refuse_overwrite": True,
}
EXPECTED_PROOF_KEYS = {
    "PROTO10_and_PRO10_FRZ1_must_validate_from_immutable_checkpoint",
    "HLT7_CAL6_ID2_and_CAL6_plan_must_be_canonical_and_hash_matched_from_immutable_checkpoint",
    "CAL7_plan_must_differ_from_CAL6_only_by_the_frozen_spectral_interpretation_and_namespace",
    "old_or_reordered_ladder_must_fail_the_PROTO10_runtime",
    "direct_finest_pair_route_must_remain_available_without_saturation",
    "coarse_pair_failure_must_not_use_saturation",
    "missing_tail_erasure_map_witness_must_veto_saturation",
    "missing_whole_profile_or_round_trip_contraction_must_veto_saturation",
    "all_twelve_projected_inputs_must_match_HLT7",
    "all_twelve_t0_accepted_source_prechecks_must_pass",
    "all_four_raw_PROTO9_t0_failures_must_reproduce_exactly",
    "all_four_guarded_PROTO10_t0_admissions_must_pass",
    "inherited_transaction_source_boundary_health_and_affine_stops_must_remain_hash_bound",
    "authorized_runner_must_be_hash_bound_and_refuse_drift_overwrite_or_dirty_tracked_state",
    "fresh_namespaces_must_be_observed_empty_without_creation",
    "canonical_hash_bound_result_required",
    "no_PROTO10_trajectory_may_be_read_or_advanced",
}
EXPECTED_CLAIMS = {
    "PROTO10_successor_runtime_implemented": True,
    "PROTO10_fresh_GR0_dynamic_calibration_authorized": True,
    "PROTO10_resolved_holdout_manifest_authorized": False,
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
IMPLEMENTATION = (
    REPOSITORY / "src/recursive_horizons/fgc/evolution/calibration_runtime.py",
    REPOSITORY / "src/recursive_horizons/fgc/evolution/cal2_source_diagnosis.py",
    REPOSITORY / "src/recursive_horizons/fgc/evolution/proto5_runtime.py",
    REPOSITORY / "src/recursive_horizons/fgc/evolution/proto6_runtime.py",
    REPOSITORY / "src/recursive_horizons/fgc/evolution/proto7_runtime.py",
    REPOSITORY / "src/recursive_horizons/fgc/evolution/cal4_common_event_diagnosis.py",
    REPOSITORY / "src/recursive_horizons/fgc/evolution/proto8_runtime.py",
    REPOSITORY / "src/recursive_horizons/fgc/evolution/proto9_runtime.py",
    REPOSITORY / "src/recursive_horizons/fgc/evolution/spectral_sensitivity.py",
    REPOSITORY / "src/recursive_horizons/fgc/evolution/proto10_runtime.py",
    REPOSITORY / "src/recursive_horizons/fgc/evolution/protocol_v10.py",
    REPOSITORY / "scripts/run_fgc_gr0_calibration.py",
    REPOSITORY / "scripts/run_fgc_gr0_calibration_v6.py",
    REPOSITORY / "scripts/run_fgc_gr0_calibration_v7.py",
    REPOSITORY / "scripts/run_fgc_gr0_calibration_v8.py",
    REPOSITORY / "scripts/run_fgc_gr0_calibration_v9.py",
    REPOSITORY / "scripts/run_fgc_gr0_calibration_v10.py",
    Path(__file__).resolve(),
)


def _path(name: str, value: Any, expected: str) -> Path:
    if value != expected:
        raise ValueError(f"HLT8 {name} differs")
    path = REPOSITORY / expected
    if not path.is_file():
        raise ValueError(f"HLT8 {name} is absent")
    return path


def validate_run_plan(path: Path) -> dict[str, Any]:
    """Prove CAL7 differs from immutable CAL6 only by PROTO10's freeze."""

    with path.open("rb") as handle:
        raw = tomllib.load(handle)
    if (
        raw.get("schema_version") != 1
        or raw.get("artifact_id") != "FGC-1-CAL7-RUN1-PLAN"
        or raw.get("project_version") != PROJECT_VERSION
        or raw.get("protocol_config") != EXPECTED_PATHS["protocol_config"]
        or raw.get("protocol_freeze_result")
        != EXPECTED_PATHS["protocol_freeze_result"]
        or raw.get("runtime_authorization_result")
        != "results/fgc-1-hlt8-mon8.json"
        or raw.get("static_input_result") != EXPECTED_PATHS["static_input_result"]
    ):
        raise ValueError("CAL7 run-plan identity or ownership differs")
    if raw.get("scope") != {
        "target_protocol": "FGC-2-SF1-PROTO10",
        "branch": "GR-0",
        "role": "fresh_outcome_neutral_dynamic_amplitude_recalibration_after_guarded_finest_pair_spectral_interpretation",
        "physical_equations": "unredefined_GR0_specialization_of_REF1",
        "candidate_action_fields_or_health_stops_forbidden": True,
        "retained_EFT_interpretation": False,
    }:
        raise ValueError("CAL7 scope differs")
    numerics = raw.get("numerics", {})
    required_spectral_delta = {
        "common_event_spatial_coarse_to_medium_tail_ratios_required_directly": True,
        "common_event_spatial_finest_pair_direct_or_conditioning_guarded_saturation_required": True,
        "common_event_spatial_guarded_saturation_requires_medium_and_fine_map_witnesses": True,
        "common_event_spatial_guarded_saturation_requires_two_pair_round_trip_contraction": True,
        "common_event_spatial_guarded_saturation_requires_whole_profile_infinity_and_RMS_contraction": True,
        "common_event_spatial_round_trip_is_not_a_continuum_error_bound": True,
        "common_event_spatial_saturation_is_not_physical_resolution": True,
    }
    if (
        numerics.get("resolutions") != list(PROTO9_POINT_COUNTS)
        or any(numerics.get(key) is not value for key, value in required_spectral_delta.items())
        or "common_event_spatial_all_adjacent_tail_ratios_required" in numerics
    ):
        raise ValueError("CAL7 spectral interpretation or ladder differs")
    if raw.get("provenance") != {
        "output_root": "runs/fgc-2-sf1/proto10/calibration",
        "campaign_manifest": "runs/fgc-2-sf1/proto10/calibration/manifest.json",
        "append_only_event_log_name": "events.jsonl",
        "checkpoint_extension": ".npz",
        "fresh_output_root_required_at_campaign_start": True,
        "runner_refuses_overwrite": True,
        "resume_requires_matching_plan_and_manifest_hashes": True,
        "PROTO9_outputs_forbidden_as_PROTO10_outcomes": True,
        "canonical_reduced_result_is_tracked_only_after_campaign_classification": True,
    }:
        raise ValueError("CAL7 provenance differs")

    normalized = deepcopy(raw)
    normalized["artifact_id"] = "FGC-1-CAL6-RUN1-PLAN"
    normalized["protocol_config"] = "configs/fgc/fgc-2-sf1-protocol-v9.toml"
    normalized["protocol_freeze_result"] = "results/fgc-1-pro9-frz1.json"
    normalized["runtime_authorization_result"] = "results/fgc-1-hlt7-mon7.json"
    normalized["scope"]["target_protocol"] = "FGC-2-SF1-PROTO9"
    normalized["scope"]["role"] = (
        "fresh_outcome_neutral_dynamic_amplitude_recalibration_after_resolution_ladder_shift"
    )
    for key in required_spectral_delta:
        normalized["numerics"].pop(key)
    normalized["numerics"][
        "common_event_spatial_all_adjacent_tail_ratios_required"
    ] = True
    normalized["provenance"] = {
        "output_root": "runs/fgc-2-sf1/proto9/calibration",
        "campaign_manifest": "runs/fgc-2-sf1/proto9/calibration/manifest.json",
        "append_only_event_log_name": "events.jsonl",
        "checkpoint_extension": ".npz",
        "fresh_output_root_required_at_campaign_start": True,
        "runner_refuses_overwrite": True,
        "resume_requires_matching_plan_and_manifest_hashes": True,
        "PROTO8_outputs_forbidden_as_PROTO9_outcomes": True,
        "canonical_reduced_result_is_tracked_only_after_campaign_classification": True,
    }
    baseline_blob = hlt4._git_blob(
        CHECKPOINT_COMMIT,
        EXPECTED_PATHS["predecessor_run_plan"],
        "immutable CAL6 plan",
    )
    if sha256(baseline_blob).hexdigest() != EXPECTED_LINEAGE[
        "predecessor_run_plan_sha256"
    ]:
        raise ValueError("immutable CAL6 plan hash differs")
    if normalized != tomllib.loads(baseline_blob.decode("utf-8")):
        raise ValueError("CAL7 changes more than the declared PROTO10 delta")
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
        "runtime_binding",
        "input_freeze",
        "namespace",
        "proof_contract",
        "claims",
    }
    if set(raw) != expected_root:
        raise ValueError("HLT8 root keys differ")
    if (
        raw["schema_version"] != 1
        or raw["artifact_id"] != ARTIFACT_ID
        or raw["project_version"] != PROJECT_VERSION
        or raw["metric_signature"] != "-+++"
        or raw["riemann_convention"] != "plus_partial_mu_gamma_nu"
    ):
        raise ValueError("HLT8 identity or convention differs")
    paths = {
        name: _path(name, raw[name], expected)
        for name, expected in EXPECTED_PATHS.items()
    }
    if raw["scope"] != EXPECTED_SCOPE:
        raise ValueError("HLT8 scope differs")
    if raw["immutable_lineage"] != EXPECTED_LINEAGE:
        raise ValueError("HLT8 immutable lineage differs")
    if raw["runtime_binding"] != EXPECTED_RUNTIME:
        raise ValueError("HLT8 runtime binding differs")
    if raw["input_freeze"] != EXPECTED_INPUT:
        raise ValueError("HLT8 input freeze differs")
    if raw["namespace"] != EXPECTED_NAMESPACE:
        raise ValueError("HLT8 namespace differs")
    if set(raw["proof_contract"]) != EXPECTED_PROOF_KEYS or any(
        value is not True for value in raw["proof_contract"].values()
    ):
        raise ValueError("HLT8 proof contract differs")
    if raw["claims"] != EXPECTED_CLAIMS:
        raise ValueError("HLT8 claims differ")
    return {
        "raw": raw,
        "paths": paths,
        "run_plan": validate_run_plan(paths["run_plan_config"]),
    }


def _immutable_lineage(loaded: Mapping[str, Any]) -> dict[str, Any]:
    lineage = loaded["raw"]["immutable_lineage"]
    commit = lineage["checkpoint_commit"]
    if subprocess.run(
        ["git", "merge-base", "--is-ancestor", commit, "HEAD"],
        cwd=REPOSITORY,
        check=False,
        capture_output=True,
    ).returncode != 0:
        raise ValueError("HLT8 checkpoint is not an ancestor of HEAD")
    immutable_names = (
        "protocol_config",
        "protocol_freeze_config",
        "protocol_freeze_result",
        "predecessor_runtime_result",
        "diagnosis_result",
        "static_input_result",
        "predecessor_run_plan",
    )
    blobs = {
        name: hlt4._git_blob(commit, EXPECTED_PATHS[name], name)
        for name in immutable_names
    }
    for name, blob in blobs.items():
        expected = lineage[f"{name}_sha256"]
        if sha256(blob).hexdigest() != expected:
            raise ValueError(f"immutable {name} hash differs")
        current = REPOSITORY / EXPECTED_PATHS[name]
        if sha256(current.read_bytes()).hexdigest() != expected:
            raise ValueError(f"current {name} differs from immutable checkpoint")

    protocol = validate_sf1_protocol_v10(
        tomllib.loads(blobs["protocol_config"].decode("utf-8"))
    )
    freeze = hlt4._load_json_bytes(
        blobs["protocol_freeze_result"], "immutable PRO10-FRZ1"
    )
    predecessor = hlt4._load_json_bytes(
        blobs["predecessor_runtime_result"], "immutable HLT7"
    )
    diagnosis = hlt4._load_json_bytes(
        blobs["diagnosis_result"], "immutable CAL6/PREF8"
    )
    static = hlt4._load_json_bytes(blobs["static_input_result"], "immutable ID2")
    inputs = predecessor.get("artifact_payload", {}).get("frozen_run_inputs", [])
    events = predecessor.get("artifact_payload", {}).get("initial_common_events", [])
    aggregate = diagnosis.get("artifact_payload", {}).get("aggregate", {})
    if (
        protocol.get("artifact_id") != "FGC-2-SF1-PROTO10"
        or protocol.get("frozen") is not True
        or protocol.get("point_counts") != list(PROTO9_POINT_COUNTS)
        or protocol.get("coarse_to_medium_direct_ratio_veto") is not True
        or freeze.get("artifact_id") != "FGC-1-PRO10-FRZ1"
        or any(value is not False for value in freeze.get("gate_status", {}).values())
        or predecessor.get("artifact_id") != "FGC-1-HLT7-MON7"
        or predecessor.get("gate_status", {}).get(
            "PROTO9_successor_runtime_implemented"
        )
        is not True
        or predecessor.get("gate_status", {}).get(
            "PROTO9_fresh_GR0_dynamic_calibration_authorized"
        )
        is not False
        or diagnosis.get("artifact_id") != "FGC-1-CAL6-PREF8"
        or aggregate.get("case_count") != 4
        or aggregate.get("raw_failed_finest_pair_metric_count") != 16
        or aggregate.get("prospective_diagnostically_saturated_metric_count") != 16
        or aggregate.get("all_prospective_resolved_or_saturated_admissions_passed")
        is not True
        or static.get("gate_status", {}).get("all_case_static_ledger_completed")
        is not True
        or len(inputs) != 12
        or len(events) != 4
        or any(item.get("admission_passed") is not False for item in events)
        or any(value is not False for value in predecessor.get("nonclaims", {}).values())
        or any(value is not False for value in diagnosis.get("nonclaims", {}).values())
    ):
        raise ValueError("immutable HLT8 predecessor gate differs")

    hash_sets = (
        freeze.get("implementation_sha256", {}),
        predecessor.get("implementation_sha256", {}),
        diagnosis.get("implementation_sha256", {}),
    )
    checked: set[str] = set()
    for hashes in hash_sets:
        for relative, expected in hashes.items():
            current = REPOSITORY / relative
            if not current.is_file() or hlt4._sha(current) != expected:
                raise ValueError(f"inherited PROTO10 premise drifted: {relative}")
            checked.add(relative)
    return {
        "checkpoint_commit": commit,
        "checkpoint_is_ancestor_of_HEAD": True,
        "protocol": {
            "artifact_id": protocol["artifact_id"],
            "protocol_version": protocol["protocol_version"],
            "point_counts": protocol["point_counts"],
            "maximum_nested_tail_ratio": protocol["maximum_nested_tail_ratio"],
            "semantic_holdout_contract_sha256": protocol[
                "semantic_holdout_contract_sha256"
            ],
        },
        "freeze": {
            "artifact_id": freeze["artifact_id"],
            "result_sha256": lineage["protocol_freeze_result_sha256"],
            "claims_all_false": True,
        },
        "predecessor_runtime": {
            "artifact_id": predecessor["artifact_id"],
            "result_sha256": lineage["predecessor_runtime_result_sha256"],
            "raw_t0_event_count": len(events),
        },
        "diagnosis": {
            "artifact_id": diagnosis["artifact_id"],
            "result_sha256": lineage["diagnosis_result_sha256"],
            "raw_failed_finest_pair_metric_count": aggregate[
                "raw_failed_finest_pair_metric_count"
            ],
            "prospective_saturated_metric_count": aggregate[
                "prospective_diagnostically_saturated_metric_count"
            ],
        },
        "static_input": {
            "artifact_id": static["artifact_id"],
            "result_sha256": lineage["static_input_result_sha256"],
        },
        "inherited_implementation_files_hash_matched": len(checked),
        "predecessor_run_plan_sha256": lineage["predecessor_run_plan_sha256"],
        "predecessor_frozen_inputs": inputs,
        "predecessor_raw_events": events,
    }


def _flat_state(grid: UniformRadialGrid, order: int) -> EvolutionState:
    u = np.zeros((grid.point_count, 6), dtype=np.float64)
    u[:, 0] = 1.0
    u[:, 2] = 1.0
    u[:, 3] = grid.coordinates
    p = np.zeros_like(u)
    q = SBPFirstDerivative(grid, order).differentiate(
        u, center_parities=ADM_CENTER_PARITIES
    )
    return EvolutionState(u, p, q)


def _adapter_controls() -> dict[str, Any]:
    old_grids = tuple(
        UniformRadialGrid(0.0, 128.0, count) for count in (1025, 2049, 4097)
    )
    old_states = tuple(_flat_state(grid, 4) for grid in old_grids)
    old_rejected = False
    try:
        proto10_gr0_common_event(
            old_states,
            old_grids,
            accepted_stage_counts=(0, 0, 0),
            method="RK4",
            coordinate_time=0.0,
        )
    except ValueError:
        old_rejected = True
    reordered_rejected = False
    try:
        proto10_gr0_common_event(
            tuple(reversed(old_states)),
            tuple(reversed(old_grids)),
            accepted_stage_counts=(0, 0, 0),
            method="RK4",
            coordinate_time=0.0,
        )
    except ValueError:
        reordered_rejected = True

    direct_grids = tuple(
        UniformRadialGrid(0.0, 128.0, count) for count in PROTO9_POINT_COUNTS
    )
    direct_states = tuple(_flat_state(grid, 4) for grid in direct_grids)
    direct = proto10_gr0_common_event(
        direct_states,
        direct_grids,
        accepted_stage_counts=(0, 0, 0),
        method="RK4",
        coordinate_time=0.0,
    )
    direct_classes = (
        *direct.spatial_spectral_admission.finest_pair_field_classification.values(),
        *direct.spatial_spectral_admission.finest_pair_derivative_classification.values(),
    )
    if (
        not old_rejected
        or not reordered_rejected
        or direct.raw_spatial_spectral_admission.admission_passed is not True
        or direct.spatial_spectral_admission.admission_passed is not True
        or direct.spatial_spectral_admission.saturation_used is not False
        or any(value != "directly_resolved" for value in direct_classes)
    ):
        raise RuntimeError("PROTO10 adapter controls failed")
    return {
        "old_PROTO8_ladder_rejected": True,
        "reordered_ladder_rejected": True,
        "correct_ladder_direct_spectral_control_passed": True,
        "direct_finest_pair_route_used_saturation": False,
        "synthetic_direct_control_asserted_as_physical_constraint_data": False,
        "all_controls_passed": True,
    }


def _raw_event_projection(
    amplitude: str,
    assessment: Proto10GR0CommonEventAssessment,
) -> dict[str, Any]:
    raw = assessment.raw_PROTO9_common_event
    return {
        "amplitude": amplitude,
        "method": raw.method,
        "coordinate_time": raw.coordinate_time,
        "point_counts": list(raw.point_counts),
        "accepted_stage_counts": list(raw.accepted_stage_counts),
        "excluded_outer_rows": list(raw.constraint_admission.excluded_outer_rows),
        "owned_raw_global_norms": list(
            raw.constraint_admission.owned_raw_global_norms
        ),
        "constraint_admission_passed": raw.constraint_admission.admission_passed,
        "constraint_coarsest_guard_passed": (
            raw.constraint_admission.coarsest_guard_passed
        ),
        "constraint_finest_guard_passed": (
            raw.constraint_admission.finest_guard_passed
        ),
        "constraint_monotone_refinement_passed": (
            raw.constraint_admission.monotone_refinement_passed
        ),
        "constraint_finest_pair_order_passed": (
            raw.constraint_admission.finest_pair_order_passed
        ),
        "constraint_minimum_finite_finest_pair_order": (
            raw.constraint_admission.minimum_finite_finest_pair_order
        ),
        "coarsest_spatial_individual_budget_passed": (
            raw.spatial_spectral_admission.coarsest_individual_budget_passed
        ),
        "finest_pair_spatial_individual_budgets_passed": (
            raw.spatial_spectral_admission.finest_pair_individual_budgets_passed
        ),
        "all_adjacent_spatial_tail_ratios_passed": (
            raw.spatial_spectral_admission.every_nested_tail_ratio_passed
        ),
        "field_power_tail_ratios": {
            name: list(values)
            for name, values in raw.spatial_spectral_admission.field_power_tail_ratios.items()
        },
        "derivative_power_tail_ratios": {
            name: list(values)
            for name, values in raw.spatial_spectral_admission.derivative_power_tail_ratios.items()
        },
        "maximum_spatial_round_trip_interpolation_infinity": (
            raw.maximum_spatial_round_trip_interpolation_infinity
        ),
        "new_8193_state_present": True,
        "conditional_projection_used": False,
        "admission_passed": raw.admission_passed,
        "trajectory_advanced": False,
    }


def _freeze_inputs(
    loaded: Mapping[str, Any],
    lineage: Mapping[str, Any],
) -> tuple[
    list[dict[str, Any]],
    list[dict[str, Any]],
    list[dict[str, Any]],
    list[Proto10GR0CommonEventAssessment],
]:
    plan = loaded["run_plan"]
    physical = plan["physical_inputs"]
    numerics = plan["numerics"]
    thresholds = plan["universal_thresholds"]
    plan_hash = hlt4._sha(loaded["paths"]["run_plan_config"])
    protocol_hash = loaded["raw"]["immutable_lineage"]["protocol_config_sha256"]
    predecessor_inputs = lineage["predecessor_frozen_inputs"]
    predecessor_by_key = {
        (item["amplitude"], item["method"], item["point_count"]): item
        for item in predecessor_inputs
    }
    expected_order = [
        (amplitude, method, count)
        for amplitude in ("5/2", "3")
        for method in ("RK4", "SSPRK3")
        for count in PROTO9_POINT_COUNTS
    ]
    if [
        (item["amplitude"], item["method"], item["point_count"])
        for item in predecessor_inputs
    ] != expected_order:
        raise ValueError("immutable HLT7 input order differs")
    predecessor_events = {
        (item["amplitude"], item["method"]): item
        for item in lineage["predecessor_raw_events"]
    }

    frozen: list[dict[str, Any]] = []
    source_prechecks: list[dict[str, Any]] = []
    states: dict[tuple[str, str], list[tuple[int, EvolutionState, Any]]] = {}
    for amplitude in ("5/2", "3"):
        for method in ("RK4", "SSPRK3"):
            order = 4 if method == "RK4" else 2
            for count in PROTO9_POINT_COUNTS:
                predecessor = predecessor_by_key[(amplitude, method, count)]
                initial = construct_gr0_grid_initial_data(
                    PulseParameters(
                        chi_amplitude=float(Q(amplitude)),
                        center=float(Q(physical["chi_center"])),
                        half_width=float(Q(physical["chi_half_width"])),
                        phi_amplitude=float(Q(physical["phi_seed_amplitude"])),
                        planck_mass=float(Q(physical["planck_mass"])),
                        scalar_mass=float(Q(physical["scalar_mass"])),
                        quartic_coupling=float(Q(physical["quartic_coupling"])),
                    ),
                    point_count=count,
                    outer_radius=float(Q(physical["outer_radius"])),
                    constraint_method=method,
                    diagnostic_spatial_order=order,
                )
                state = project_gr0_semidiscrete_state(initial, spatial_order=order)
                state_hash = array_content_sha256(state.u, state.p, state.q)
                native_q = SBPFirstDerivative(initial.grid, order).differentiate(
                    state.u,
                    center_parities=ADM_CENTER_PARITIES,
                )
                if (
                    state_hash != predecessor["projected_state_sha256"]
                    or not np.array_equal(state.u, initial.state.u)
                    or not np.array_equal(state.p, initial.state.p)
                    or not np.array_equal(state.q, native_q)
                ):
                    raise ValueError("PROTO10 projected input differs from immutable HLT7")

                operator = Proto7GR0EvolutionOperator(
                    initial.grid,
                    spatial_order=order,
                    ko_dissipation=float(Q(numerics["ko_dissipation"])),
                    raw_tolerance=float(
                        Q(thresholds["source_residual_infinity_max"])
                    ),
                    kinetic_condition_maximum=float(
                        Q(thresholds["kinetic_condition_number_max"])
                    ),
                )
                rhs = operator(0.0, state)
                residual = float(rhs.diagnostics["source_residual_infinity"])
                threshold = float(Q(thresholds["source_residual_infinity_max"]))
                if residual >= threshold:
                    raise ValueError("PROTO10 t0 accepted-state source precheck failed")
                source_prechecks.append(
                    {
                        "amplitude": amplitude,
                        "method": method,
                        "point_count": count,
                        "source_residual_infinity": residual,
                        "source_residual_maximum": threshold,
                        "accepted_state_source_gate_passed": True,
                        "trajectory_advanced": False,
                    }
                )

                record = dict(predecessor)
                record["predecessor_PROTO9_expanded_run_config_sha256"] = record.pop(
                    "expanded_run_config_sha256"
                )
                record["predecessor_PROTO9_resolution_role"] = record.pop(
                    "PROTO9_resolution_role"
                )
                record.update(
                    {
                        "plan_sha256": plan_hash,
                        "protocol_sha256": protocol_hash,
                        "PROTO10_resolution_role": (
                            "direct_state_with_guarded_finest_pair_diagnostic"
                        ),
                        "projected_state_matches_HLT7": True,
                        "diagnostic_saturation_is_state_mutation": False,
                        "trajectory_advanced": False,
                    }
                )
                record["expanded_run_config_sha256"] = hlt4._digest(record)
                frozen.append(record)
                states.setdefault((amplitude, method), []).append(
                    (count, state, initial.grid)
                )

    event_payloads: list[dict[str, Any]] = []
    assessments: list[Proto10GR0CommonEventAssessment] = []
    for amplitude in ("5/2", "3"):
        for method in ("RK4", "SSPRK3"):
            records = sorted(states[(amplitude, method)], key=lambda item: item[0])
            assessed = proto10_gr0_common_event(
                [item[1] for item in records],
                [item[2] for item in records],
                accepted_stage_counts=(0, 0, 0),
                method=method,
                coordinate_time=0.0,
                cutoff=float(Q(physical["cutoff_Lambda"])),
                measurement_radius_maximum=float(
                    Q(physical["measurement_radius_maximum"])
                ),
                taper_fraction=float(Q(numerics["proper_spectral_taper_fraction"])),
                fixed_outer_rows=numerics["fixed_outer_rows"],
            )
            raw_projection = _raw_event_projection(amplitude, assessed)
            if raw_projection != predecessor_events[(amplitude, method)]:
                raise ValueError("PROTO10 raw t0 event differs from immutable HLT7")
            guarded = assessed.spatial_spectral_admission
            if (
                assessed.raw_PROTO9_admission_passed is not False
                or assessed.constraint_admission.admission_passed is not True
                or guarded.finest_pair_individual_budgets_passed is not True
                or guarded.every_profile_contracted is not True
                or guarded.saturation_used is not True
                or guarded.every_tail_resolved_or_saturated is not True
                or assessed.admission_passed is not True
            ):
                raise ValueError("PROTO10 t0 guarded admission shape differs")
            event_payloads.append(
                {
                    "amplitude": amplitude,
                    "method": method,
                    "coordinate_time": 0.0,
                    "raw_PROTO9_event_projection": raw_projection,
                    "guarded_PROTO10_common_event": hlt4._serial(assessed),
                    "raw_PROTO9_admission_passed": False,
                    "guarded_PROTO10_admission_passed": True,
                    "diagnostic_saturation_used": True,
                    "diagnostic_saturation_is_physical_resolution": False,
                    "round_trip_is_a_continuum_error_bound": False,
                    "trajectory_advanced": False,
                }
            )
            assessments.append(assessed)
    return frozen, source_prechecks, event_payloads, assessments


def _guard_controls(
    assessments: list[Proto10GR0CommonEventAssessment],
) -> dict[str, Any]:
    event = assessments[0]
    guarded = event.spatial_spectral_admission
    saturated_name = next(
        name
        for name in PROTO4_SPECTRAL_FIELD_ORDER
        if (
            guarded.finest_pair_field_classification[name]
            == "diagnostically_saturated"
            or guarded.finest_pair_derivative_classification[name]
            == "diagnostically_saturated"
        )
    )

    budgets = [dict(record) for record in event.spectral_budgets]
    coarse_name = next(
        name
        for name in PROTO4_SPECTRAL_FIELD_ORDER
        if budgets[0][name].top_band_field_power_fraction > 0.0
    )
    coarse_value = budgets[0][coarse_name].top_band_field_power_fraction
    budgets[1][coarse_name] = replace(
        budgets[1][coarse_name],
        top_band_field_power_fraction=coarse_value / 2.0,
    )
    coarse_attack = resolved_or_saturated_admission(
        PROTO9_POINT_COUNTS,
        budgets,
        event.spectral_tail_sensitivities,
        event.profile_convergence,
    )

    sensitivities = [dict(record) for record in event.spectral_tail_sensitivities]
    sensitivities[2][saturated_name] = replace(
        sensitivities[2][saturated_name], diagnostically_saturated=False
    )
    sensitivity_attack = resolved_or_saturated_admission(
        PROTO9_POINT_COUNTS,
        event.spectral_budgets,
        sensitivities,
        event.profile_convergence,
    )

    complete = dict(event.profile_convergence.complete_profile_contraction_by_field)
    complete[saturated_name] = False
    profile_attack_record = replace(
        event.profile_convergence,
        complete_profile_contraction_by_field=complete,
        every_field_contracted=False,
    )
    profile_attack = resolved_or_saturated_admission(
        PROTO9_POINT_COUNTS,
        event.spectral_budgets,
        event.spectral_tail_sensitivities,
        profile_attack_record,
    )

    round_trip = dict(event.profile_convergence.round_trip_contracted_by_field)
    round_trip[saturated_name] = False
    round_trip_complete = dict(
        event.profile_convergence.complete_profile_contraction_by_field
    )
    round_trip_complete[saturated_name] = False
    round_trip_record = replace(
        event.profile_convergence,
        round_trip_contracted_by_field=round_trip,
        complete_profile_contraction_by_field=round_trip_complete,
        every_field_contracted=False,
    )
    round_trip_attack = resolved_or_saturated_admission(
        PROTO9_POINT_COUNTS,
        event.spectral_budgets,
        event.spectral_tail_sensitivities,
        round_trip_record,
    )
    if (
        coarse_attack.admission_passed
        or sensitivity_attack.admission_passed
        or profile_attack.admission_passed
        or round_trip_attack.admission_passed
    ):
        raise RuntimeError("PROTO10 adversarial guard control failed open")
    defaults = asdict(SpectralThresholds())
    return {
        "attacked_event_amplitude": "5/2",
        "attacked_event_method": "RK4",
        "coarse_direct_failure_field": coarse_name,
        "coarse_direct_failure_admitted": False,
        "missing_tail_erasure_map_witness_field": saturated_name,
        "missing_tail_erasure_map_witness_admitted": False,
        "missing_complete_profile_contraction_admitted": False,
        "missing_round_trip_contraction_admitted": False,
        "spectral_thresholds": defaults,
        "maximum_nested_tail_ratio": defaults["maximum_nested_tail_ratio"],
        "new_epsilon_floor_added": False,
        "all_adversarial_controls_passed": True,
    }


def _namespace_precondition(
    raw: Mapping[str, Any],
    historical: Mapping[str, Any] | None,
) -> dict[str, Any]:
    if historical is not None:
        records = historical.get("records", [])
        if (
            historical.get("both_fresh_namespaces_empty") is not True
            or historical.get("authorization_created_no_namespace") is not True
            or [item.get("path") for item in records]
            != [
                raw["namespace"]["calibration_output_root"],
                raw["namespace"]["holdout_output_root"],
            ]
            or any(item.get("empty") is not True for item in records)
        ):
            raise ValueError("historical HLT8 namespace evidence differs")
        return dict(historical)
    records = []
    for key in ("calibration_output_root", "holdout_output_root"):
        relative = raw["namespace"][key]
        path = REPOSITORY / relative
        exists = path.exists()
        if exists and (not path.is_dir() or any(path.iterdir())):
            raise ValueError(f"HLT8 fresh namespace is nonempty: {relative}")
        records.append(
            {
                "path": relative,
                "existed_before_authorization": exists,
                "entry_count": 0,
                "empty": True,
                "created_by_authorization": False,
            }
        )
    return {
        "both_fresh_namespaces_empty": True,
        "authorization_created_no_namespace": True,
        "records": records,
    }


def record(
    config_path: Path = DEFAULT_CONFIG,
    *,
    historical_namespace_evidence: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    loaded = load_config(config_path)
    raw = loaded["raw"]
    lineage = _immutable_lineage(loaded)
    adapter_controls = _adapter_controls()
    inputs, source_prechecks, events, assessments = _freeze_inputs(loaded, lineage)
    guard_controls = _guard_controls(assessments)
    if (
        len(events) != 4
        or any(item["raw_PROTO9_admission_passed"] for item in events)
        or any(not item["guarded_PROTO10_admission_passed"] for item in events)
        or any(not item["diagnostic_saturation_used"] for item in events)
    ):
        raise ValueError("HLT8 t0 authorization evidence is incomplete")
    namespace = _namespace_precondition(raw, historical_namespace_evidence)
    plan_hash = hlt4._sha(loaded["paths"]["run_plan_config"])
    input_manifest_hash = hlt4._digest(inputs)
    campaign_id = hlt4._digest(
        {
            "artifact_id": loaded["run_plan"]["artifact_id"],
            "run_plan_sha256": plan_hash,
            "input_manifest_sha256": input_manifest_hash,
            "ordered_amplitudes": ["5/2", "3"],
            "point_counts": list(PROTO9_POINT_COUNTS),
            "output_root": raw["namespace"]["calibration_output_root"],
        }
    )
    public_lineage = dict(lineage)
    public_lineage.pop("predecessor_frozen_inputs")
    public_lineage.pop("predecessor_raw_events")
    return {
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "project_version": PROJECT_VERSION,
        "classification": (
            "pre_trajectory_PROTO10_guarded_spectral_runtime_and_fresh_GR0_calibration_authorization"
        ),
        "generated_by": "scripts/reproduce_fgc_hlt8_mon8.py",
        "derivation_document": hlt4._rel(OWNER_DOCUMENT),
        "derivation_document_sha256": hlt4._sha(OWNER_DOCUMENT),
        "source_config_sha256": {
            hlt4._rel(config_path): hlt4._sha(config_path),
            hlt4._rel(loaded["paths"]["run_plan_config"]): plan_hash,
        },
        "predecessor_sha256": {
            EXPECTED_PATHS[name]: raw["immutable_lineage"][f"{name}_sha256"]
            for name in (
                "protocol_config",
                "protocol_freeze_config",
                "protocol_freeze_result",
                "predecessor_runtime_result",
                "diagnosis_result",
                "static_input_result",
                "predecessor_run_plan",
            )
        },
        "implementation_sha256": {
            hlt4._rel(path): hlt4._sha(path) for path in IMPLEMENTATION
        },
        "scope_bindings": dict(raw["scope"]),
        "gate_status": dict(raw["claims"]),
        "artifact_payload": {
            "immutable_lineage": public_lineage,
            "runtime_adapter_controls": adapter_controls,
            "adversarial_guard_controls": guard_controls,
            "inherited_runtime": {
                "PROTO9_raw_common_event_and_PROTO7_evolution_unchanged": True,
                "inherited_implementation_files_hash_matched": public_lineage[
                    "inherited_implementation_files_hash_matched"
                ],
                "process_local_runner_bindings_are_exception_safe": True,
            },
            "frozen_run_plan": {
                "artifact_id": loaded["run_plan"]["artifact_id"],
                "run_plan_sha256": plan_hash,
                "campaign_id": campaign_id,
                "input_manifest_sha256": input_manifest_hash,
                "ordered_amplitudes": ["5/2", "3"],
                "point_counts": list(PROTO9_POINT_COUNTS),
                "first_eligible_candidate_stops_later_candidates": True,
                "expanded_input_count": len(inputs),
                "projected_state_hash_match_count": sum(
                    item["projected_state_matches_HLT7"] for item in inputs
                ),
                "raw_PROTO9_t0_passed_count": sum(
                    item["raw_PROTO9_admission_passed"] for item in events
                ),
                "guarded_PROTO10_t0_passed_count": sum(
                    item["guarded_PROTO10_admission_passed"] for item in events
                ),
                "t0_common_event_evaluated_count": len(events),
            },
            "frozen_run_inputs": inputs,
            "accepted_state_source_prechecks": source_prechecks,
            "initial_common_events": events,
            "namespace_precondition": namespace,
            "decision": (
                "HLT8_implements_the_PROTO10_compositor_and_authorizes_one_fresh_"
                "GR0_calibration_but_keeps_every_candidate_and_physical_claim_closed"
            ),
            "epistemic_boundary": {
                "PROTO9_preflight_and_diagnosis_history_disclosed": True,
                "PROTO10_trajectory_read": False,
                "PROTO10_trajectory_advanced": False,
                "t0_physical_inputs_reconstructed": True,
                "all_twelve_source_prechecks_passed": True,
                "all_four_raw_PROTO9_admissions_failed_and_remain_public": True,
                "all_four_guarded_PROTO10_admissions_passed": True,
                "diagnostic_saturation_is_physical_resolution": False,
                "round_trip_is_a_continuum_error_bound": False,
                "fresh_GR0_recalibration_authorized": True,
                "fresh_GR0_recalibration_completed": False,
                "trapped_or_collapse_outcome_classified": False,
                "resolved_holdout_manifest_exists": False,
                "SGBL_or_FGCQR_outcome_read": False,
                "mechanism_question_answered": False,
            },
        },
        "nonclaims": {
            "fresh_GR0_recalibration_completed": False,
            "trapped_sphere_formed_dynamically": False,
            "SGBL_health_definition_supplied": False,
            "SGBL_comparison_completed": False,
            "FGCQR_trajectory_opened": False,
            "regulator_activation_observed": False,
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


def load_canonical_result(path: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    return hlt4._load_json_bytes(path.read_bytes(), "HLT8 result")


def verify_canonical(
    path: Path = DEFAULT_OUTPUT,
    config_path: Path = DEFAULT_CONFIG,
) -> None:
    observed = load_canonical_result(path)
    historical = observed.get("artifact_payload", {}).get("namespace_precondition")
    expected = record(
        config_path,
        historical_namespace_evidence=historical,
    )
    if observed != hlt4._serial(expected):
        raise ValueError("HLT8 output differs from a fresh reproduction")


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
    output.write_text(hlt4._canonical(result), encoding="utf-8")
    print(
        f"wrote {output} "
        "(PROTO10 runtime=true; fresh GR0 calibration authorized=true; "
        "trajectory_read=false)"
    )


if __name__ == "__main__":
    main()
