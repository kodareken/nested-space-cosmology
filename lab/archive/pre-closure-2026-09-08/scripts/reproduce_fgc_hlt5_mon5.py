#!/usr/bin/env python3
"""Regenerate the pre-trajectory FGC-1-HLT5-MON5 authorization record."""

from __future__ import annotations

import argparse
from copy import deepcopy
from dataclasses import asdict
from pathlib import Path
from types import SimpleNamespace
import subprocess
import sys
import tomllib
from typing import Any, Mapping

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-hlt5-mon5.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-hlt5-mon5.json"
OWNER_DOCUMENT = REPOSITORY / "docs/fgc-hlt5-mon5.md"

from scripts import reproduce_fgc_hlt4_mon4 as hlt4  # noqa: E402
from recursive_horizons.fgc.evolution.boundary_domain import (  # noqa: E402
    BoundaryGeometry,
    CausalBudgetState,
)
from recursive_horizons.fgc.evolution.calibration_runtime import (  # noqa: E402
    project_gr0_semidiscrete_state,
)
from recursive_horizons.fgc.evolution.gr0_calibration import (  # noqa: E402
    GR0EvolutionOperator,
    construct_gr0_grid_initial_data,
)
from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    EvolutionRHS,
    EvolutionState,
    PRIMARY_METHOD,
    UniformRadialGrid,
    array_content_sha256,
)
from recursive_horizons.fgc.evolution.proto5_runtime import (  # noqa: E402
    GR0RuntimeStageTransaction,
    GR0RuntimeStop,
    GR0UniversalThresholds,
    proto5_gr0_common_event,
)
from recursive_horizons.fgc.evolution.proto6_runtime import (  # noqa: E402
    restore_common_event_snapshots,
)
from recursive_horizons.fgc.evolution.proto7_runtime import (  # noqa: E402
    Proto7GR0EvolutionOperator,
    Proto7TerminalStop,
    attempt_proto7_step,
)
from recursive_horizons.fgc.evolution.protocol_v7 import (  # noqa: E402
    validate_sf1_protocol_v7,
)
from recursive_horizons.fgc.initial_data_preflight import (  # noqa: E402
    PulseParameters,
)


Q = hlt4.Q
ARTIFACT_ID = "FGC-1-HLT5-MON5"
PROJECT_VERSION = "0.11.0"
CHECKPOINT_COMMIT = "ec8ac628528b822c600de35eb2ca6b2b49506051"
EXPECTED_PATHS = {
    "protocol_config": "configs/fgc/fgc-2-sf1-protocol-v7.toml",
    "protocol_freeze_config": "configs/fgc/fgc-1-pro7-frz1.toml",
    "protocol_freeze_result": "results/fgc-1-pro7-frz1.json",
    "predecessor_runtime_result": "results/fgc-1-hlt4-mon4.json",
    "static_input_result": "results/fgc-1-id2-all1.json",
    "run_plan_config": "configs/fgc/fgc-1-cal4-run1.toml",
}
EXPECTED_SCOPE = {
    "target_protocol": "FGC-2-SF1-PROTO7",
    "runtime_owner": ARTIFACT_ID,
    "calibration_branch": "GR-0",
    "role": "pre_trajectory_general_unaccepted_proposal_runtime_compositor_and_fresh_input_freeze",
    "PROTO7_trajectory_read": False,
    "SGBL_trajectory_read": False,
    "FGCQR_trajectory_read": False,
    "collapse_or_trapped_outcome_classified": False,
    "candidate_action_health_values_invented_for_GR0": False,
}
EXPECTED_LINEAGE = {
    "checkpoint_commit": CHECKPOINT_COMMIT,
    "protocol_config_sha256": "ea0563b205bf5c6c90d8f728afc56a13f624953fffc67c6d6b9a97c03c35f3fb",
    "protocol_freeze_config_sha256": "130008687df18e1dd4fee3ef83ab6337b96c23bd253f69d776ef383cf30f3e5b",
    "protocol_freeze_result_sha256": "6cb16d2e80e27f69af2f8dab0a74b583447fb7f84b1877d5aed5b1cc7cf40800",
    "predecessor_runtime_result_sha256": "9b5082a5775263345d23f7f7b43b9c05fde1653b0d76b7a34553d1ac9d93f452",
    "static_input_result_sha256": "10e2661be521087467b970f1cd877311bd63a0cb6874a1ee745ca00ec74ab981",
    "checkpoint_must_be_ancestor_of_HEAD": True,
    "all_predecessor_blobs_must_be_read_from_checkpoint": True,
}
EXPECTED_RUNTIME = {
    "accepted_state_source_gate_evaluated_before_every_proposal": True,
    "accepted_state_source_failure_is_terminal_and_latched": True,
    "complete_shadow_transaction_precedes_retry_classification": True,
    "every_source_only_unaccepted_proposal_failure_is_retryable": True,
    "candidate_endpoint_is_part_of_the_unaccepted_proposal": True,
    "any_non_source_failure_vetoes_retry": True,
    "inherited_chronological_and_within_stage_priority_is_unchanged": True,
    "retry_restarts_from_bitwise_last_accepted_state": True,
    "retry_advances_no_fields_time_stage_transaction_causal_debit_or_tracer": True,
    "retry_factor": "1/2",
    "maximum_source_retries_per_step": 32,
    "minimum_step_size": "1/1073741824",
    "raw_source_threshold_solver_and_refinement_unchanged": True,
    "all_non_source_PROTO6_stops_unchanged": True,
    "every_rejected_proposal_is_fsynced_to_the_event_ledger": True,
    "terminal_and_exhaustion_records_retain_complete_failure_evidence": True,
    "cross_member_failure_rolls_back_to_last_common_event": True,
    "PROTO6_common_event_spectral_trapped_and_recovery_controls_inherited": True,
}
EXPECTED_INPUT = {
    "eligible_amplitudes_in_order": ["5/2", "3"],
    "methods_in_order": ["RK4", "SSPRK3"],
    "resolutions_in_order": [1025, 2049, 4097],
    "expected_run_input_count": 12,
    "projected_state_hashes_must_match_HLT4": True,
    "physical_u_p_q_and_initial_profiles_unchanged": True,
    "expanded_PROTO7_run_config_hashes_must_be_serialized": True,
    "all_twelve_accepted_state_source_prechecks_must_pass": True,
    "both_amplitudes_and_methods_must_pass_the_t0_common_event": True,
    "first_eligible_candidate_stops_later_candidates": True,
}
EXPECTED_NAMESPACE = {
    "calibration_output_root": "runs/fgc-2-sf1/proto7/calibration",
    "holdout_output_root": "runs/fgc-2-sf1/proto7/holdout",
    "both_roots_must_be_absent_or_empty_at_authorization": True,
    "authorization_must_not_create_either_root": True,
    "runner_must_refuse_overwrite": True,
}
EXPECTED_PROOF_KEYS = {
    "PROTO7_and_PRO7_FRZ1_must_validate_from_immutable_checkpoint",
    "HLT4_and_ID2_must_be_canonical_and_hash_matched_from_immutable_checkpoint",
    "PROTO7_run_plan_must_differ_from_PROTO6_only_by_declared_ownership_observability_and_namespace",
    "accepted_source_failure_injected_control_must_stop_and_latch",
    "internal_and_candidate_endpoint_source_only_controls_must_retry_without_mutation",
    "non_source_veto_controls_must_stop_at_every_proposal_position",
    "complete_failure_set_and_selected_reason_controls_must_pass",
    "healthy_transaction_and_preaccept_order_control_must_pass",
    "accepted_source_operator_must_match_original_GR0_operator",
    "all_twelve_projected_inputs_must_reconstruct_deterministically",
    "all_twelve_t0_accepted_source_prechecks_must_pass",
    "authorized_runner_must_be_hash_bound_and_refuse_drift_overwrite_or_dirty_tracked_state",
    "retry_ledger_checkpoint_resume_exhaustion_and_cross_member_rollback_controls_must_pass",
    "fresh_namespaces_must_be_observed_empty_without_creation",
    "canonical_hash_bound_result_required",
    "no_trajectory_may_be_read_or_advanced",
}
EXPECTED_CLAIMS = {
    "PROTO7_successor_runtime_compositor_implemented": True,
    "PROTO7_fresh_GR0_dynamic_calibration_authorized": True,
    "PROTO7_resolved_holdout_manifest_authorized": False,
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
    REPOSITORY / "scripts/run_fgc_gr0_calibration.py",
    REPOSITORY / "scripts/run_fgc_gr0_calibration_v6.py",
    REPOSITORY / "scripts/run_fgc_gr0_calibration_v7.py",
    Path(__file__).resolve(),
)


def _path(name: str, value: Any, expected: str) -> Path:
    if value != expected:
        raise ValueError(f"HLT5 {name} differs")
    path = REPOSITORY / expected
    if not path.is_file():
        raise ValueError(f"HLT5 {name} is absent")
    return path


def validate_run_plan(path: Path) -> dict[str, Any]:
    with path.open("rb") as handle:
        raw = tomllib.load(handle)
    if (
        raw.get("schema_version") != 1
        or raw.get("artifact_id") != "FGC-1-CAL4-RUN1-PLAN"
        or raw.get("project_version") != PROJECT_VERSION
        or raw.get("protocol_config") != EXPECTED_PATHS["protocol_config"]
        or raw.get("protocol_freeze_result") != EXPECTED_PATHS["protocol_freeze_result"]
        or raw.get("runtime_authorization_result") != "results/fgc-1-hlt5-mon5.json"
    ):
        raise ValueError("CAL4 run-plan identity or ownership differs")
    if raw.get("scope") != {
        "target_protocol": "FGC-2-SF1-PROTO7",
        "branch": "GR-0",
        "role": "fresh_outcome_neutral_dynamic_amplitude_recalibration_after_general_transaction_repair",
        "physical_equations": "unredefined_GR0_specialization_of_REF1",
        "candidate_action_fields_or_health_stops_forbidden": True,
        "retained_EFT_interpretation": False,
    }:
        raise ValueError("CAL4 scope differs")
    numerics = raw.get("numerics", {})
    if (
        numerics.get("maximum_source_retries_per_step") != 32
        or numerics.get("accepted_state_source_precheck_before_every_proposal") is not True
        or numerics.get("serialize_every_rejected_source_trial") is not True
        or numerics.get("candidate_endpoint_source_only_failure_retryable") is not True
        or numerics.get("complete_terminal_and_exhaustion_failure_evidence") is not True
    ):
        raise ValueError("CAL4 PROTO7 numerical ownership differs")
    if raw.get("provenance") != {
        "output_root": "runs/fgc-2-sf1/proto7/calibration",
        "campaign_manifest": "runs/fgc-2-sf1/proto7/calibration/manifest.json",
        "append_only_event_log_name": "events.jsonl",
        "checkpoint_extension": ".npz",
        "fresh_output_root_required_at_campaign_start": True,
        "runner_refuses_overwrite": True,
        "resume_requires_matching_plan_and_manifest_hashes": True,
        "PROTO6_outputs_forbidden_as_PROTO7_outcomes": True,
        "canonical_reduced_result_is_tracked_only_after_campaign_classification": True,
    }:
        raise ValueError("CAL4 provenance differs")

    normalized = deepcopy(raw)
    normalized["artifact_id"] = "FGC-1-CAL3-RUN1-PLAN"
    normalized["protocol_config"] = "configs/fgc/fgc-2-sf1-protocol-v6.toml"
    normalized["protocol_freeze_result"] = "results/fgc-1-pro6-frz1.json"
    normalized["runtime_authorization_result"] = "results/fgc-1-hlt4-mon4.json"
    normalized["scope"]["target_protocol"] = "FGC-2-SF1-PROTO6"
    normalized["scope"]["role"] = "fresh_outcome_neutral_dynamic_amplitude_recalibration"
    normalized["numerics"].pop("candidate_endpoint_source_only_failure_retryable")
    normalized["numerics"].pop("complete_terminal_and_exhaustion_failure_evidence")
    normalized["provenance"] = {
        "output_root": "runs/fgc-2-sf1/proto6/calibration",
        "campaign_manifest": "runs/fgc-2-sf1/proto6/calibration/manifest.json",
        "append_only_event_log_name": "events.jsonl",
        "checkpoint_extension": ".npz",
        "fresh_output_root_required_at_campaign_start": True,
        "runner_refuses_overwrite": True,
        "resume_requires_matching_plan_and_manifest_hashes": True,
        "PROTO5_outputs_forbidden_as_PROTO6_outcomes": True,
        "canonical_reduced_result_is_tracked_only_after_campaign_classification": True,
    }
    baseline_blob = hlt4._git_blob(
        CHECKPOINT_COMMIT,
        "configs/fgc/fgc-1-cal3-run1.toml",
        "immutable CAL3 plan",
    )
    baseline = tomllib.loads(baseline_blob.decode("utf-8"))
    if normalized != baseline:
        raise ValueError("CAL4 changes more than the declared PROTO7 delta")
    return raw


def load_config(path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    with path.open("rb") as handle:
        raw = tomllib.load(handle)
    expected_root = {
        "schema_version", "artifact_id", "project_version", "metric_signature",
        "riemann_convention", *EXPECTED_PATHS, "scope", "immutable_lineage",
        "runtime_composition", "input_freeze", "namespace", "proof_contract",
        "claims",
    }
    if set(raw) != expected_root:
        raise ValueError("HLT5 root keys differ")
    if (
        raw["schema_version"] != 1
        or raw["artifact_id"] != ARTIFACT_ID
        or raw["project_version"] != PROJECT_VERSION
        or raw["metric_signature"] != "-+++"
        or raw["riemann_convention"] != "plus_partial_mu_gamma_nu"
    ):
        raise ValueError("HLT5 identity or convention differs")
    paths = {
        name: _path(name, raw[name], expected)
        for name, expected in EXPECTED_PATHS.items()
    }
    if raw["scope"] != EXPECTED_SCOPE:
        raise ValueError("HLT5 scope differs")
    if raw["immutable_lineage"] != EXPECTED_LINEAGE:
        raise ValueError("HLT5 immutable lineage differs")
    if raw["runtime_composition"] != EXPECTED_RUNTIME:
        raise ValueError("HLT5 runtime composition differs")
    if raw["input_freeze"] != EXPECTED_INPUT:
        raise ValueError("HLT5 input freeze differs")
    if raw["namespace"] != EXPECTED_NAMESPACE:
        raise ValueError("HLT5 namespace differs")
    if set(raw["proof_contract"]) != EXPECTED_PROOF_KEYS or any(
        value is not True for value in raw["proof_contract"].values()
    ):
        raise ValueError("HLT5 proof contract differs")
    if raw["claims"] != EXPECTED_CLAIMS:
        raise ValueError("HLT5 claims differ")
    run_plan = validate_run_plan(paths["run_plan_config"])
    return {"raw": raw, "paths": paths, "run_plan": run_plan}


def _immutable_lineage(loaded: Mapping[str, Any]) -> dict[str, Any]:
    lineage = loaded["raw"]["immutable_lineage"]
    commit = lineage["checkpoint_commit"]
    if subprocess.run(
        ["git", "merge-base", "--is-ancestor", commit, "HEAD"],
        cwd=REPOSITORY,
        check=False,
        capture_output=True,
    ).returncode != 0:
        raise ValueError("HLT5 checkpoint is not an ancestor of HEAD")
    names = (
        "protocol_config", "protocol_freeze_config", "protocol_freeze_result",
        "predecessor_runtime_result", "static_input_result",
    )
    blobs = {
        name: hlt4._git_blob(commit, EXPECTED_PATHS[name], name) for name in names
    }
    for name, blob in blobs.items():
        if hlt4.sha256(blob).hexdigest() != lineage[f"{name}_sha256"]:
            raise ValueError(f"immutable {name} hash differs")
    protocol = validate_sf1_protocol_v7(
        tomllib.loads(blobs["protocol_config"].decode("utf-8"))
    )
    freeze = hlt4._load_json_bytes(
        blobs["protocol_freeze_result"], "immutable PRO7-FRZ1"
    )
    predecessor = hlt4._load_json_bytes(
        blobs["predecessor_runtime_result"], "immutable HLT4"
    )
    static = hlt4._load_json_bytes(blobs["static_input_result"], "immutable ID2")
    if (
        protocol.get("artifact_id") != "FGC-2-SF1-PROTO7"
        or freeze.get("artifact_id") != "FGC-1-PRO7-FRZ1"
        or freeze.get("gate_status", {}).get("PROTO7_outcome_neutral_protocol_frozen")
        is not True
        or predecessor.get("artifact_id") != "FGC-1-HLT4-MON4"
        or predecessor.get("gate_status", {}).get(
            "PROTO6_successor_runtime_compositor_implemented"
        )
        is not True
        or static.get("gate_status", {}).get("all_case_static_ledger_completed")
        is not True
    ):
        raise ValueError("immutable HLT5 predecessor gate differs")
    frozen_inputs = predecessor.get("artifact_payload", {}).get(
        "frozen_run_inputs", []
    )
    if len(frozen_inputs) != 12:
        raise ValueError("immutable HLT4 does not carry twelve inputs")
    return {
        "checkpoint_commit": commit,
        "checkpoint_is_ancestor_of_HEAD": True,
        "protocol": protocol,
        "protocol_freeze_artifact_id": freeze["artifact_id"],
        "protocol_freeze_gate_passed": True,
        "predecessor_runtime_artifact_id": predecessor["artifact_id"],
        "predecessor_runtime_gate_passed": True,
        "static_input_artifact_id": static["artifact_id"],
        "static_input_gate_passed": True,
        "predecessor_frozen_inputs": frozen_inputs,
    }


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


def _synthetic_rhs(
    *, source_calls: set[int] | None = None, lapse_calls: set[int] | None = None
):
    calls = 0
    sources = set() if source_calls is None else set(source_calls)
    lapses = set() if lapse_calls is None else set(lapse_calls)

    def rhs(_time: float, state: EvolutionState) -> EvolutionRHS:
        nonlocal calls
        index = calls
        calls += 1
        zero = np.zeros(state.shape)
        return EvolutionRHS(
            zero, zero, zero,
            {
                "source_residual_infinity": 2.0e-12 if index in sources else 0.0,
                "source_refinement_iterations": 0,
                "source_residual_decreased_monotonically": True,
                "kinetic_condition_infinity": 1.0,
                "coordinate_speed_upper": 1.0,
                "minimum_lapse": 0.0 if index in lapses else 1.0,
                "minimum_radial_metric": 1.0,
                "minimum_areal_radius_away_from_center": 0.125,
            },
        )

    return rhs


def _transaction() -> GR0RuntimeStageTransaction:
    return GR0RuntimeStageTransaction(
        thresholds=GR0UniversalThresholds(),
        causal_state=CausalBudgetState(previous_speed_upper=1.0),
        boundary_geometry=BoundaryGeometry(128.0, 24.0, 16.0, 0.375),
        grid_spacing=0.125,
        cfl_maximum=0.125,
    )


def _attempt(rhs, transaction, *, preaccept=None):
    return attempt_proto7_step(
        method=PRIMARY_METHOD,
        time=0.0,
        step_size=0.01,
        state=_synthetic_state(),
        rhs=rhs,
        projector=None,
        transaction=transaction,
        previous_step_index=0,
        previous_transaction_serial=0,
        preaccept=preaccept,
    )


def _runtime_controls() -> dict[str, Any]:
    accepted_transaction = _transaction()
    try:
        _attempt(_synthetic_rhs(source_calls={0}), accepted_transaction)
    except Proto7TerminalStop as stop:
        accepted = {
            "reason": stop.reason,
            "evidence": asdict(stop.evidence),
        }
    else:
        raise RuntimeError("accepted-state source miss became retryable")

    retry_positions = []
    for call, expected_stage in (
        (1, "rk4_k1"), (2, "rk4_k2"), (3, "rk4_k3"),
        (4, "rk4_k4"), (5, "candidate_endpoint"),
    ):
        transaction = _transaction()
        state_before = transaction.state
        causal_before = transaction.causal_state
        result = _attempt(_synthetic_rhs(source_calls={call}), transaction)
        retry_positions.append(
            {
                "stage": expected_stage,
                "observed_stage": result.retry.failed_evaluations[0].stage_name,
                "retryable": result.retry.retryable_source_only,
                "state_unchanged": transaction.state == state_before,
                "causal_debit_unchanged": transaction.causal_state == causal_before,
            }
        )

    veto_positions = []
    for call in range(1, 6):
        transaction = _transaction()
        try:
            _attempt(
                _synthetic_rhs(source_calls={call}, lapse_calls={call}),
                transaction,
            )
        except Proto7TerminalStop as stop:
            veto_positions.append(
                {
                    "call": call,
                    "reason": stop.reason,
                    "non_source_failure_vetoed_retry": (
                        stop.evidence.non_source_failure_vetoed_retry
                    ),
                    "complete_failure_set": list(stop.evidence.complete_failure_set),
                }
            )
        else:
            raise RuntimeError("non-source proposal failure became retryable")

    mixed_transaction = _transaction()
    try:
        _attempt(
            _synthetic_rhs(source_calls={2}, lapse_calls={4}),
            mixed_transaction,
        )
    except Proto7TerminalStop as stop:
        mixed = {
            "selected_reason": stop.reason,
            "complete_failure_set": list(stop.evidence.complete_failure_set),
            "failed_stage_names": [
                item.stage_name for item in stop.evidence.failed_evaluations
            ],
        }
    else:
        raise RuntimeError("later non-source failure did not veto source retry")

    healthy_transaction = _transaction()
    ordering: list[int] = []
    healthy = _attempt(
        _synthetic_rhs(),
        healthy_transaction,
        preaccept=lambda _proposal: ordering.append(
            healthy_transaction.state.accepted_stage_count
        ),
    )
    all_passed = (
        accepted["reason"] == "newton_residual_limit"
        and all(
            item["stage"] == item["observed_stage"]
            and item["retryable"]
            and item["state_unchanged"]
            and item["causal_debit_unchanged"]
            for item in retry_positions
        )
        and all(
            item["reason"] == "nonpositive_lapse"
            and item["non_source_failure_vetoed_retry"]
            for item in veto_positions
        )
        and mixed["selected_reason"] == "newton_residual_limit"
        and mixed["failed_stage_names"] == ["rk4_k2", "rk4_k4"]
        and healthy.accepted is not None
        and ordering == [0]
    )
    if not all_passed:
        raise RuntimeError("PROTO7 transaction controls failed")
    return {
        "accepted_state_terminal_control": accepted,
        "source_only_retry_controls_at_every_proposal_position": retry_positions,
        "non_source_veto_controls_at_every_proposal_position": veto_positions,
        "mixed_chronological_priority_control": mixed,
        "healthy_preaccept_then_commit_control": {
            "accepted": True,
            "preaccept_observed_before_commit": True,
            "accepted_stage_count": healthy_transaction.state.accepted_stage_count,
        },
        "all_controls_passed": True,
    }


def _runner_controls() -> dict[str, Any]:
    from scripts.run_fgc_gr0_calibration_v7 import (
        Proto7RunMember,
        Proto7SourceRetryExhausted,
    )

    class Tracers:
        def __init__(self) -> None:
            self.positions = np.array([0.5])
            self.proper_times = np.array([0.0])
            self.event_proper_times: list[np.ndarray] = []
            self.event_fields: list[np.ndarray] = []
            self.preview_count = 0
            self.commit_count = 0

        def preview_advance(self, **_kwargs):
            self.preview_count += 1
            return np.array([1.0]), np.array([2.0])

        def commit_advance(self, positions, proper_times) -> None:
            self.commit_count += 1
            self.positions = positions
            self.proper_times = proper_times

        def append_common_event(self, *_args) -> None:
            self.event_proper_times.append(self.proper_times.copy())
            self.event_fields.append(np.zeros((1, 6)))

    def member(rhs) -> Proto7RunMember:
        grid = UniformRadialGrid(0.0, 1.0, 9)
        return Proto7RunMember(
            amplitude="5/2", method_label="RK4", integrator_id=PRIMARY_METHOD,
            spatial_order=4, point_count=9, input_hash="0" * 64,
            initial=SimpleNamespace(grid=grid), state=_synthetic_state(),
            operator=rhs, projector=None, transaction=_transaction(),
            tracers=Tracers(),
        )

    accepted_member = member(_synthetic_rhs(source_calls={5}))
    retry_records: list[Mapping[str, Any]] = []
    accepted_member.advance_to(
        0.01,
        retry_factor=0.5,
        maximum_CFL_retries=32,
        maximum_source_retries=32,
        minimum_step_size=2.0**-30,
        rejected_trial_sink=retry_records.append,
    )
    exhausted_member = member(_synthetic_rhs(source_calls={5, 11}))
    exhaustion_records: list[Mapping[str, Any]] = []
    try:
        exhausted_member.advance_to(
            0.01,
            retry_factor=0.5,
            maximum_CFL_retries=32,
            maximum_source_retries=1,
            minimum_step_size=2.0**-30,
            rejected_trial_sink=exhaustion_records.append,
        )
    except Proto7SourceRetryExhausted as stop:
        exhaustion = stop.evidence
    else:
        raise RuntimeError("PROTO7 retry exhaustion control did not stop")

    rollback_member = member(_synthetic_rhs())
    snapshot = rollback_member.snapshot()
    rollback_member.time = 1.0
    rollback_member.step_index = 9
    restore_common_event_snapshots(
        {"member": rollback_member}, {"member": snapshot}
    )
    all_passed = (
        len(retry_records) == 1
        and retry_records[0]["failed_evaluations"][0]["stage_name"]
        == "candidate_endpoint"
        and retry_records[0]["complete_failure_evidence_serialized_before_retry"]
        is True
        and retry_records[0]["external_transaction_state_preserved"] is True
        and accepted_member.time == 0.01
        and len(exhaustion_records) == 2
        and exhaustion["exhaustion_kind"] == "maximum_source_retries"
        and "last_rejected_proposal" in exhaustion
        and rollback_member.time == 0.0
        and rollback_member.step_index == 0
    )
    if not all_passed:
        raise RuntimeError("PROTO7 runner controls failed")
    return {
        "candidate_endpoint_retry_then_accept": {
            "serialized_retry_count": len(retry_records),
            "final_time": accepted_member.time,
            "accepted_step_count": accepted_member.step_index,
        },
        "retry_exhaustion": {
            "serialized_retry_count": len(exhaustion_records),
            "evidence": exhaustion,
        },
        "cross_member_restore_adapter": {
            "time_restored": rollback_member.time == 0.0,
            "step_index_restored": rollback_member.step_index == 0,
        },
        "all_controls_passed": True,
    }


def _operator_equivalence() -> dict[str, Any]:
    records = []
    for order, method in ((4, "RK4"), (2, "SSPRK3")):
        initial = construct_gr0_grid_initial_data(
            PulseParameters(chi_amplitude=2.5),
            point_count=129,
            outer_radius=128.0,
            constraint_method=method,
            diagnostic_spatial_order=order,
        )
        state = project_gr0_semidiscrete_state(initial, spatial_order=order)
        original = GR0EvolutionOperator(
            initial.grid, spatial_order=order, ko_dissipation=1.0 / 64.0,
            residual_tolerance=1.0e-12, kinetic_condition_maximum=1.0e10,
        )(0.0, state)
        successor = Proto7GR0EvolutionOperator(
            initial.grid, spatial_order=order, ko_dissipation=1.0 / 64.0,
            raw_tolerance=1.0e-12, kinetic_condition_maximum=1.0e10,
        )(0.0, state)
        if not all(
            np.array_equal(left, right)
            for left, right in (
                (original.du, successor.du),
                (original.dp, successor.dp),
                (original.dq, successor.dq),
            )
        ):
            raise RuntimeError("PROTO7 accepted source differs from original GR-0")
        records.append(
            {
                "method": method,
                "point_count": 129,
                "binary64_RHS_arrays_identical": True,
                "source_residual_infinity": successor.diagnostics[
                    "source_residual_infinity"
                ],
            }
        )
    return {"controls": records, "all_controls_passed": True}


def _freeze_inputs(
    loaded: Mapping[str, Any],
    lineage: Mapping[str, Any],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    plan = loaded["run_plan"]
    physical = plan["physical_inputs"]
    numerics = plan["numerics"]
    thresholds = plan["universal_thresholds"]
    plan_hash = hlt4._sha(loaded["paths"]["run_plan_config"])
    protocol_hash = loaded["raw"]["immutable_lineage"]["protocol_config_sha256"]
    predecessor_inputs = lineage["predecessor_frozen_inputs"]
    expected_order = [
        (amplitude, method, count)
        for amplitude in ("5/2", "3")
        for method in ("RK4", "SSPRK3")
        for count in (1025, 2049, 4097)
    ]
    if [
        (item["amplitude"], item["method"], item["point_count"])
        for item in predecessor_inputs
    ] != expected_order:
        raise ValueError("immutable HLT4 input order differs")

    frozen: list[dict[str, Any]] = []
    source_prechecks: list[dict[str, Any]] = []
    states: dict[tuple[str, str], list[tuple[int, EvolutionState, Any]]] = {}
    for predecessor in predecessor_inputs:
        amplitude = predecessor["amplitude"]
        method = predecessor["method"]
        count = predecessor["point_count"]
        order = 4 if method == "RK4" else 2
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
        if state_hash != predecessor["projected_state_sha256"]:
            raise ValueError("PROTO7 projected state differs from immutable HLT4")
        operator = Proto7GR0EvolutionOperator(
            initial.grid,
            spatial_order=order,
            ko_dissipation=float(Q(numerics["ko_dissipation"])),
            raw_tolerance=float(Q(thresholds["source_residual_infinity_max"])),
            kinetic_condition_maximum=float(
                Q(thresholds["kinetic_condition_number_max"])
            ),
        )
        rhs = operator(0.0, state)
        residual = float(rhs.diagnostics["source_residual_infinity"])
        threshold = float(Q(thresholds["source_residual_infinity_max"]))
        if residual >= threshold:
            raise ValueError("PROTO7 t0 accepted-state source precheck failed")
        source_prechecks.append(
            {
                "amplitude": amplitude, "method": method, "point_count": count,
                "source_residual_infinity": residual,
                "source_residual_maximum": threshold,
                "accepted_state_source_gate_passed": True,
                "trajectory_advanced": False,
            }
        )
        record = dict(predecessor)
        record["plan_sha256"] = plan_hash
        record["protocol_sha256"] = protocol_hash
        record["source_stop_ownership"] = (
            "accepted_terminal_all_unaccepted_proposal_source_only_retry"
        )
        record["complete_failure_observability"] = True
        record["trajectory_advanced"] = False
        record.pop("expanded_run_config_sha256", None)
        record["expanded_run_config_sha256"] = hlt4._digest(record)
        frozen.append(record)
        states.setdefault((amplitude, method), []).append((count, state, initial.grid))

    common_events: list[dict[str, Any]] = []
    for amplitude in ("5/2", "3"):
        for method in ("RK4", "SSPRK3"):
            records = sorted(states[(amplitude, method)], key=lambda item: item[0])
            assessed = proto5_gr0_common_event(
                [item[1] for item in records],
                [item[2] for item in records],
                method=method,
                coordinate_time=0.0,
                cutoff=float(Q(physical["cutoff_Lambda"])),
                measurement_radius_maximum=float(
                    Q(physical["measurement_radius_maximum"])
                ),
                taper_fraction=float(Q(numerics["proper_spectral_taper_fraction"])),
            )
            if not assessed.admission_passed:
                raise ValueError("PROTO7 initial common event failed")
            common_events.append(
                {
                    "amplitude": amplitude, "method": method,
                    "coordinate_time": 0.0,
                    "point_counts": list(assessed.point_counts),
                    "admission_passed": True, "trajectory_advanced": False,
                }
            )
    return frozen, source_prechecks, common_events


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
            != [raw["namespace"]["calibration_output_root"], raw["namespace"]["holdout_output_root"]]
            or any(item.get("empty") is not True for item in records)
        ):
            raise ValueError("historical HLT5 namespace evidence differs")
        return dict(historical)
    records = []
    for key in ("calibration_output_root", "holdout_output_root"):
        relative = raw["namespace"][key]
        path = REPOSITORY / relative
        exists = path.exists()
        if exists and (not path.is_dir() or any(path.iterdir())):
            raise ValueError(f"HLT5 fresh namespace is nonempty: {relative}")
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
    runtime = _runtime_controls()
    runner = _runner_controls()
    equivalence = _operator_equivalence()
    inputs, source_prechecks, common_events = _freeze_inputs(loaded, lineage)
    namespace = _namespace_precondition(raw, historical_namespace_evidence)
    plan_hash = hlt4._sha(loaded["paths"]["run_plan_config"])
    input_manifest_hash = hlt4._digest(inputs)
    campaign_id = hlt4._digest(
        {
            "artifact_id": loaded["run_plan"]["artifact_id"],
            "run_plan_sha256": plan_hash,
            "input_manifest_sha256": input_manifest_hash,
            "ordered_amplitudes": ["5/2", "3"],
            "output_root": raw["namespace"]["calibration_output_root"],
        }
    )
    public_lineage = dict(lineage)
    public_lineage.pop("predecessor_frozen_inputs")
    return {
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "project_version": PROJECT_VERSION,
        "classification": "pre_trajectory_PROTO7_GR0_general_transaction_runtime_and_fresh_calibration_authorization",
        "generated_by": "scripts/reproduce_fgc_hlt5_mon5.py",
        "derivation_document": hlt4._rel(OWNER_DOCUMENT),
        "derivation_document_sha256": hlt4._sha(OWNER_DOCUMENT),
        "source_config_sha256": {
            hlt4._rel(config_path): hlt4._sha(config_path),
            hlt4._rel(loaded["paths"]["run_plan_config"]): plan_hash,
        },
        "predecessor_sha256": {
            EXPECTED_PATHS[name]: raw["immutable_lineage"][f"{name}_sha256"]
            for name in (
                "protocol_config", "protocol_freeze_config",
                "protocol_freeze_result", "predecessor_runtime_result",
                "static_input_result",
            )
        },
        "implementation_sha256": {
            hlt4._rel(path): hlt4._sha(path) for path in IMPLEMENTATION
        },
        "scope_bindings": dict(raw["scope"]),
        "gate_status": dict(raw["claims"]),
        "artifact_payload": {
            "immutable_lineage": public_lineage,
            "runtime_controls": runtime,
            "runner_controls": runner,
            "accepted_source_operator_equivalence": equivalence,
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
            "accepted_state_source_prechecks": source_prechecks,
            "initial_common_events": common_events,
            "namespace_precondition": namespace,
            "decision": "HLT5_implements_only_PROTO7_transaction_ownership_and_authorizes_one_fresh_GR0_recalibration",
            "epistemic_boundary": {
                "PROTO6_calibration_history_disclosed": True,
                "PROTO7_trajectory_read": False,
                "PROTO7_trajectory_advanced": False,
                "synthetic_failures_used_only_for_transaction_controls": True,
                "t0_physical_inputs_reconstructed": True,
                "trapped_or_collapse_outcome_classified": False,
                "SGBL_or_FGCQR_outcome_read": False,
                "fresh_GR0_recalibration_is_now_authorized_but_not_completed": True,
                "resolved_holdout_manifest_exists": False,
                "mechanism_question_answered": False,
            },
        },
        "nonclaims": {
            "fresh_GR0_recalibration_completed": False,
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


def load_canonical_result(path: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    return hlt4._load_json_bytes(path.read_bytes(), "HLT5 result")


def _verify_postlaunch(
    loaded: Mapping[str, Any],
    observed: Mapping[str, Any],
    result_path: Path,
) -> None:
    raw = loaded["raw"]
    historical = observed.get("artifact_payload", {}).get(
        "namespace_precondition", {}
    )
    _namespace_precondition(raw, historical)
    calibration_root = REPOSITORY / raw["namespace"]["calibration_output_root"]
    holdout_root = REPOSITORY / raw["namespace"]["holdout_output_root"]
    if not calibration_root.is_dir():
        raise ValueError("postlaunch HLT5 calibration directory is absent")
    if holdout_root.exists() and (
        not holdout_root.is_dir() or any(holdout_root.iterdir())
    ):
        raise ValueError("HLT5 does not authorize a nonempty holdout namespace")
    manifest = hlt4._load_json_bytes(
        (calibration_root / "manifest.json").read_bytes(), "CAL4 launch manifest"
    )
    frozen = observed["artifact_payload"]["frozen_run_plan"]
    expected = {
        "schema_version": 1,
        "runner_id": "FGC-1-CAL4-RUN1-RUNNER",
        "campaign_id": frozen["campaign_id"],
        "plan_path": EXPECTED_PATHS["run_plan_config"],
        "plan_sha256": frozen["run_plan_sha256"],
        "authorization_path": hlt4._rel(result_path),
        "authorization_sha256": hlt4._sha(result_path),
        "input_manifest_sha256": frozen["input_manifest_sha256"],
        "ordered_amplitudes": frozen["ordered_amplitudes"],
        "implementation_sha256": observed["implementation_sha256"],
        "outcome_fields_present_at_creation": False,
    }
    for key, value in expected.items():
        if manifest.get(key) != value:
            raise ValueError(f"postlaunch CAL4 manifest field differs: {key}")
    commit = manifest.get("git_commit")
    if not isinstance(commit, str) or subprocess.run(
        ["git", "merge-base", "--is-ancestor", commit, "HEAD"],
        cwd=REPOSITORY,
        check=False,
        capture_output=True,
    ).returncode != 0:
        raise ValueError("postlaunch CAL4 Git checkpoint is invalid")


def verify_canonical(path: Path, config_path: Path = DEFAULT_CONFIG) -> None:
    observed = load_canonical_result(path)
    loaded = load_config(config_path)
    calibration_root = REPOSITORY / loaded["raw"]["namespace"][
        "calibration_output_root"
    ]
    if calibration_root.exists() and any(calibration_root.iterdir()):
        _verify_postlaunch(loaded, observed, path)
        expected = record(
            config_path,
            historical_namespace_evidence=observed["artifact_payload"]
            ["namespace_precondition"],
        )
    else:
        expected = record(config_path)
    if observed != hlt4._serial(expected):
        raise ValueError("HLT5 output differs from a fresh reproduction")


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
    print(f"wrote {output}")


if __name__ == "__main__":
    main()
