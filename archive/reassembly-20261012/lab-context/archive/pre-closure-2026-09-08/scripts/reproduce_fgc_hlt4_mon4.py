#!/usr/bin/env python3
"""Regenerate the pre-trajectory FGC-1-HLT4-MON4 authorization record."""

from __future__ import annotations

import argparse
from copy import deepcopy
from dataclasses import asdict, is_dataclass
from fractions import Fraction
from hashlib import sha256
import json
from math import isfinite
from pathlib import Path
import subprocess
import sys
import tomllib
from typing import Any, Mapping

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-hlt4-mon4.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-hlt4-mon4.json"
OWNER_DOCUMENT = REPOSITORY / "docs/fgc-hlt4-mon4.md"

from scripts import reproduce_fgc_hlt3_mon3 as hlt3  # noqa: E402
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
    array_content_sha256,
)
from recursive_horizons.fgc.evolution.proto5_runtime import (  # noqa: E402
    GR0RuntimeStageTransaction,
    GR0RuntimeStop,
    GR0UniversalThresholds,
    proto5_gr0_common_event,
)
from recursive_horizons.fgc.evolution.proto6_runtime import (  # noqa: E402
    Proto6GR0EvolutionOperator,
    attempt_proto6_step,
    restore_common_event_snapshots,
)
from recursive_horizons.fgc.evolution.protocol_v6 import (  # noqa: E402
    validate_sf1_protocol_v6,
)
from recursive_horizons.fgc.initial_data_preflight import (  # noqa: E402
    PulseParameters,
)


Q = Fraction
ARTIFACT_ID = "FGC-1-HLT4-MON4"
PROJECT_VERSION = "0.11.0"
EXPECTED_PATHS = {
    "protocol_config": "configs/fgc/fgc-2-sf1-protocol-v6.toml",
    "protocol_freeze_config": "configs/fgc/fgc-1-pro6-frz1.toml",
    "protocol_freeze_result": "results/fgc-1-pro6-frz1.json",
    "predecessor_runtime_result": "results/fgc-1-hlt3-mon3.json",
    "static_input_result": "results/fgc-1-id2-all1.json",
    "run_plan_config": "configs/fgc/fgc-1-cal3-run1.toml",
}
EXPECTED_SCOPE = {
    "target_protocol": "FGC-2-SF1-PROTO6",
    "runtime_owner": ARTIFACT_ID,
    "calibration_branch": "GR-0",
    "role": "pre_trajectory_source_retry_runtime_compositor_and_fresh_input_freeze",
    "PROTO6_trajectory_read": False,
    "SGBL_trajectory_read": False,
    "FGCQR_trajectory_read": False,
    "collapse_or_trapped_outcome_classified": False,
    "candidate_action_health_values_invented_for_GR0": False,
}
EXPECTED_LINEAGE = {
    "checkpoint_commit": "f3a32ff9faddfa4dba0e0f5e192d58ef94d5e094",
    "protocol_config_sha256": "e50093914f7b84a88bd593361f8b3eb43942cf9a8502f637fc38d8e61dd82448",
    "protocol_freeze_config_sha256": "7c78f645d3e5974a94b0787b000dbc9257a36458da576ec2055f6f93ef2b7a06",
    "protocol_freeze_result_sha256": "ed9a8ce8df75a3ade38a0416eb2ecb2e09c0840d89c5115860eb536d48799b56",
    "predecessor_runtime_result_sha256": "80ba72924b0fc7509b97105ba881e1df9951f5e935bf5429d9b10672bb236ec6",
    "static_input_result_sha256": "10e2661be521087467b970f1cd877311bd63a0cb6874a1ee745ca00ec74ab981",
    "checkpoint_must_be_ancestor_of_HEAD": True,
    "all_predecessor_blobs_must_be_read_from_checkpoint": True,
}
EXPECTED_RUNTIME = {
    "accepted_state_source_gate_evaluated_before_every_proposal": True,
    "accepted_state_source_failure_is_terminal_and_latched": True,
    "complete_shadow_transaction_precedes_retry_classification": True,
    "only_source_only_internal_RK_stage_failure_is_retryable": True,
    "candidate_endpoint_source_failure_remains_terminal": True,
    "simultaneous_non_source_failure_remains_terminal": True,
    "retry_restarts_from_bitwise_last_accepted_state": True,
    "retry_advances_no_time_stage_transaction_causal_debit_or_tracer": True,
    "retry_factor": "1/2",
    "maximum_source_retries_per_step": 32,
    "minimum_step_size": "1/1073741824",
    "raw_source_threshold_solver_and_refinement_unchanged": True,
    "all_non_source_PROTO5_stops_unchanged": True,
    "every_rejected_trial_is_fsynced_to_the_event_ledger": True,
    "cross_member_failure_rolls_back_to_last_common_event": True,
    "PROTO5_common_event_spectral_trapped_and_recovery_controls_inherited": True,
}
EXPECTED_INPUT = {
    "eligible_amplitudes_in_order": ["5/2", "3"],
    "methods_in_order": ["RK4", "SSPRK3"],
    "resolutions_in_order": [1025, 2049, 4097],
    "expected_run_input_count": 12,
    "projected_state_hashes_must_match_HLT3": True,
    "physical_u_p_q_and_initial_profiles_unchanged": True,
    "expanded_PROTO6_run_config_hashes_must_be_serialized": True,
    "all_twelve_accepted_state_source_prechecks_must_pass": True,
    "both_amplitudes_and_methods_must_pass_the_t0_common_event": True,
    "first_eligible_candidate_stops_later_candidates": True,
}
EXPECTED_NAMESPACE = {
    "calibration_output_root": "runs/fgc-2-sf1/proto6/calibration",
    "holdout_output_root": "runs/fgc-2-sf1/proto6/holdout",
    "both_roots_must_be_absent_or_empty_at_authorization": True,
    "authorization_must_not_create_either_root": True,
    "runner_must_refuse_overwrite": True,
}
EXPECTED_PROOF_KEYS = {
    "PROTO6_and_PRO6_FRZ1_must_validate_from_immutable_checkpoint",
    "HLT3_and_ID2_must_be_canonical_and_hash_matched_from_immutable_checkpoint",
    "PROTO6_run_plan_must_differ_from_PROTO5_only_by_declared_ownership_and_namespace",
    "accepted_source_failure_injected_control_must_stop_and_latch",
    "internal_source_only_failure_injected_control_must_retry_without_mutation",
    "candidate_endpoint_and_simultaneous_non_source_controls_must_stop",
    "healthy_transaction_and_preaccept_order_control_must_pass",
    "accepted_source_operator_must_match_original_GR0_operator",
    "all_twelve_projected_inputs_must_reconstruct_deterministically",
    "all_twelve_t0_accepted_source_prechecks_must_pass",
    "authorized_runner_must_be_hash_bound_and_refuse_drift_overwrite_or_dirty_tracked_state",
    "retry_ledger_checkpoint_resume_and_cross_member_rollback_controls_must_pass",
    "fresh_namespaces_must_be_observed_empty_without_creation",
    "canonical_hash_bound_result_required",
    "no_trajectory_may_be_read_or_advanced",
}
EXPECTED_CLAIMS = {
    "PROTO6_successor_runtime_compositor_implemented": True,
    "PROTO6_fresh_GR0_dynamic_calibration_authorized": True,
    "PROTO6_resolved_holdout_manifest_authorized": False,
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
    REPOSITORY / "scripts/run_fgc_gr0_calibration.py",
    REPOSITORY / "scripts/run_fgc_gr0_calibration_v6.py",
    REPOSITORY / "scripts/reproduce_fgc_hlt4_mon4.py",
)


def _serial(value: Any) -> Any:
    if isinstance(value, Fraction):
        return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"
    if isinstance(value, Mapping):
        return {str(key): _serial(item) for key, item in value.items()}
    if isinstance(value, np.ndarray):
        return _serial(value.tolist())
    if isinstance(value, np.generic):
        return _serial(value.item())
    if is_dataclass(value):
        return _serial(asdict(value))
    if isinstance(value, (tuple, list)):
        return [_serial(item) for item in value]
    if isinstance(value, float):
        if not isfinite(value):
            raise ValueError("nonfinite HLT4 evidence")
        return value
    if value is None or isinstance(value, (str, bool, int)):
        return value
    raise TypeError(f"unsupported HLT4 evidence type {type(value).__name__}")


def _canonical(value: Any) -> str:
    return json.dumps(
        _serial(value),
        sort_keys=True,
        indent=2,
        ensure_ascii=True,
        allow_nan=False,
    ) + "\n"


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
    return path.relative_to(REPOSITORY).as_posix()


def _reject_pairs(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    answer: dict[str, Any] = {}
    for key, value in pairs:
        if key in answer:
            raise ValueError(f"duplicate JSON key {key}")
        answer[key] = value
    return answer


def _load_json_bytes(source: bytes, name: str) -> dict[str, Any]:
    value = json.loads(source.decode("utf-8"), object_pairs_hook=_reject_pairs)
    if not isinstance(value, dict) or source.decode("utf-8") != _canonical(value):
        raise ValueError(f"{name} is not canonical sorted JSON")
    return value


def _path(name: str, value: Any, expected: str) -> Path:
    if value != expected:
        raise ValueError(f"HLT4 {name} differs")
    path = REPOSITORY / expected
    if not path.is_file():
        raise ValueError(f"HLT4 {name} is absent")
    return path


def _git_blob(commit: str, relative: str, name: str) -> bytes:
    result = subprocess.run(
        ["git", "show", f"{commit}:{relative}"],
        cwd=REPOSITORY,
        check=False,
        capture_output=True,
    )
    if result.returncode != 0:
        raise ValueError(f"cannot read immutable {name} blob")
    return result.stdout


def _validate_run_plan(path: Path) -> dict[str, Any]:
    with path.open("rb") as handle:
        raw = tomllib.load(handle)
    if (
        raw.get("schema_version") != 1
        or raw.get("artifact_id") != "FGC-1-CAL3-RUN1-PLAN"
        or raw.get("project_version") != PROJECT_VERSION
        or raw.get("protocol_config") != EXPECTED_PATHS["protocol_config"]
        or raw.get("protocol_freeze_result") != EXPECTED_PATHS["protocol_freeze_result"]
        or raw.get("runtime_authorization_result") != "results/fgc-1-hlt4-mon4.json"
    ):
        raise ValueError("CAL3 run-plan identity or ownership differs")
    if raw.get("scope") != {
        "target_protocol": "FGC-2-SF1-PROTO6",
        "branch": "GR-0",
        "role": "fresh_outcome_neutral_dynamic_amplitude_recalibration",
        "physical_equations": "unredefined_GR0_specialization_of_REF1",
        "candidate_action_fields_or_health_stops_forbidden": True,
        "retained_EFT_interpretation": False,
    }:
        raise ValueError("CAL3 scope differs")
    numerics = raw.get("numerics", {})
    if (
        numerics.get("maximum_source_retries_per_step") != 32
        or numerics.get("accepted_state_source_precheck_before_every_proposal") is not True
        or numerics.get("serialize_every_rejected_source_trial") is not True
    ):
        raise ValueError("CAL3 PROTO6 numerical ownership differs")
    if raw.get("provenance") != {
        "output_root": "runs/fgc-2-sf1/proto6/calibration",
        "campaign_manifest": "runs/fgc-2-sf1/proto6/calibration/manifest.json",
        "append_only_event_log_name": "events.jsonl",
        "checkpoint_extension": ".npz",
        "fresh_output_root_required_at_campaign_start": True,
        "runner_refuses_overwrite": True,
        "resume_requires_matching_plan_and_manifest_hashes": True,
        "PROTO5_outputs_forbidden_as_PROTO6_outcomes": True,
        "canonical_reduced_result_is_tracked_only_after_campaign_classification": True,
    }:
        raise ValueError("CAL3 provenance differs")

    # Normalize only the declared PROTO6 ownership/namespace delta and demand
    # exact equality with HLT3's already validated PROTO5 run plan.
    normalized = deepcopy(raw)
    normalized["artifact_id"] = "FGC-1-CAL2-RUN1-PLAN"
    normalized["protocol_config"] = "configs/fgc/fgc-2-sf1-protocol-v5.toml"
    normalized["protocol_freeze_result"] = "results/fgc-1-pro5-frz1.json"
    normalized["runtime_authorization_result"] = "results/fgc-1-hlt3-mon3.json"
    normalized["scope"]["target_protocol"] = "FGC-2-SF1-PROTO5"
    normalized["scope"]["role"] = "fresh_outcome_neutral_dynamic_amplitude_calibration"
    for key in (
        "maximum_source_retries_per_step",
        "accepted_state_source_precheck_before_every_proposal",
        "serialize_every_rejected_source_trial",
    ):
        normalized["numerics"].pop(key)
    normalized["provenance"] = {
        "output_root": "runs/fgc-2-sf1/proto5/calibration",
        "campaign_manifest": "runs/fgc-2-sf1/proto5/calibration/manifest.json",
        "append_only_event_log_name": "events.jsonl",
        "checkpoint_extension": ".npz",
        "fresh_output_root_required_at_campaign_start": True,
        "runner_refuses_overwrite": True,
        "resume_requires_matching_plan_and_manifest_hashes": True,
        "earlier_protocol_outputs_forbidden_as_evidence": True,
        "canonical_reduced_result_is_tracked_only_after_campaign_classification": True,
    }
    baseline = hlt3._validate_run_plan(
        REPOSITORY / "configs/fgc/fgc-1-cal2-run1.toml"
    )
    if normalized != baseline:
        raise ValueError("CAL3 changes more than the declared PROTO6 delta")
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
        raise ValueError("HLT4 root keys differ")
    if (
        raw["schema_version"] != 1
        or raw["artifact_id"] != ARTIFACT_ID
        or raw["project_version"] != PROJECT_VERSION
        or raw["metric_signature"] != "-+++"
        or raw["riemann_convention"] != "plus_partial_mu_gamma_nu"
    ):
        raise ValueError("HLT4 identity or convention differs")
    paths = {
        name: _path(name, raw[name], expected)
        for name, expected in EXPECTED_PATHS.items()
    }
    if raw["scope"] != EXPECTED_SCOPE:
        raise ValueError("HLT4 scope differs")
    if raw["immutable_lineage"] != EXPECTED_LINEAGE:
        raise ValueError("HLT4 immutable lineage differs")
    if raw["runtime_composition"] != EXPECTED_RUNTIME:
        raise ValueError("HLT4 runtime composition differs")
    if raw["input_freeze"] != EXPECTED_INPUT:
        raise ValueError("HLT4 input freeze differs")
    if raw["namespace"] != EXPECTED_NAMESPACE:
        raise ValueError("HLT4 namespace differs")
    if set(raw["proof_contract"]) != EXPECTED_PROOF_KEYS or any(
        value is not True for value in raw["proof_contract"].values()
    ):
        raise ValueError("HLT4 proof contract differs")
    if raw["claims"] != EXPECTED_CLAIMS:
        raise ValueError("HLT4 claims differ")
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
        raise ValueError("HLT4 checkpoint is not an ancestor of HEAD")
    names = (
        "protocol_config",
        "protocol_freeze_config",
        "protocol_freeze_result",
        "predecessor_runtime_result",
        "static_input_result",
    )
    blobs = {
        name: _git_blob(commit, EXPECTED_PATHS[name], name) for name in names
    }
    for name in names:
        if sha256(blobs[name]).hexdigest() != lineage[f"{name}_sha256"]:
            raise ValueError(f"immutable {name} hash differs")
    protocol = validate_sf1_protocol_v6(tomllib.loads(blobs["protocol_config"].decode()))
    freeze = _load_json_bytes(blobs["protocol_freeze_result"], "immutable PRO6-FRZ1")
    predecessor = _load_json_bytes(
        blobs["predecessor_runtime_result"], "immutable HLT3"
    )
    static = _load_json_bytes(blobs["static_input_result"], "immutable ID2")
    if (
        freeze.get("artifact_id") != "FGC-1-PRO6-FRZ1"
        or freeze.get("gate_status", {}).get("PROTO6_outcome_neutral_protocol_frozen")
        is not True
        or predecessor.get("artifact_id") != "FGC-1-HLT3-MON3"
        or predecessor.get("gate_status", {}).get(
            "PROTO5_successor_runtime_compositor_implemented"
        )
        is not True
        or static.get("gate_status", {}).get("all_case_static_ledger_completed")
        is not True
    ):
        raise ValueError("immutable HLT4 predecessor gate differs")
    frozen_inputs = predecessor.get("artifact_payload", {}).get(
        "frozen_run_inputs", []
    )
    if len(frozen_inputs) != 12:
        raise ValueError("immutable HLT3 does not carry twelve inputs")
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


def _synthetic_rhs(*, failed_call: int | None, non_source: bool = False):
    calls = 0

    def rhs(_time: float, state: EvolutionState) -> EvolutionRHS:
        nonlocal calls
        failed = calls == failed_call
        calls += 1
        zero = np.zeros(state.shape)
        return EvolutionRHS(
            zero,
            zero,
            zero,
            {
                "source_residual_infinity": 2.0e-12 if failed else 0.0,
                "source_refinement_iterations": 0,
                "source_residual_decreased_monotonically": True,
                "kinetic_condition_infinity": 1.0,
                "coordinate_speed_upper": 1.0,
                "minimum_lapse": 0.0 if failed and non_source else 1.0,
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
    return attempt_proto6_step(
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
        _attempt(_synthetic_rhs(failed_call=0), accepted_transaction)
    except GR0RuntimeStop as stop:
        accepted_stop = stop.reason
    else:
        raise RuntimeError("accepted-state source failure did not stop")

    retry_transaction = _transaction()
    state_before = retry_transaction.state
    causal_before = retry_transaction.causal_state
    preaccept_count = 0

    def preaccept(_proposal):
        nonlocal preaccept_count
        preaccept_count += 1

    retry = _attempt(
        _synthetic_rhs(failed_call=4), retry_transaction, preaccept=preaccept
    )
    if retry.retry is None:
        raise RuntimeError("internal source-only control was not retryable")

    endpoint_transaction = _transaction()
    try:
        _attempt(_synthetic_rhs(failed_call=5), endpoint_transaction)
    except GR0RuntimeStop as stop:
        endpoint_stop = stop.reason
    else:
        raise RuntimeError("candidate endpoint source failure did not stop")

    simultaneous_transaction = _transaction()
    try:
        _attempt(
            _synthetic_rhs(failed_call=4, non_source=True),
            simultaneous_transaction,
        )
    except GR0RuntimeStop as stop:
        simultaneous_stop = {
            "reason": stop.reason,
            "failures": list(stop.failures),
        }
    else:
        raise RuntimeError("simultaneous non-source failure became retryable")

    healthy_transaction = _transaction()
    healthy_preaccept_counts: list[int] = []
    healthy = _attempt(
        _synthetic_rhs(failed_call=None),
        healthy_transaction,
        preaccept=lambda _proposal: healthy_preaccept_counts.append(
            healthy_transaction.state.accepted_stage_count
        ),
    )

    class Restorable:
        def __init__(self) -> None:
            self.value = 9

        def restore(self, snapshot: Mapping[str, Any]) -> None:
            self.value = int(snapshot["value"])

    member = Restorable()
    restore_common_event_snapshots({"member": member}, {"member": {"value": 3}})
    serialized_retry = _serial(asdict(retry.retry))
    return {
        "accepted_state_control": {
            "terminal_reason": accepted_stop,
            "first_failure_latched": accepted_transaction.state.first_failed_premise
            == "newton_residual_limit",
            "accepted_stage_count": accepted_transaction.state.accepted_stage_count,
        },
        "internal_trial_control": {
            "failed_stage_name": retry.retry.failed_stage_name,
            "source_residual_infinity": retry.retry.source_residual_infinity,
            "raw_source_threshold": retry.retry.source_residual_maximum,
            "state_unchanged": retry_transaction.state == state_before,
            "causal_debit_unchanged": retry_transaction.causal_state == causal_before,
            "preaccept_not_called": preaccept_count == 0,
            "serialized_retry_fields": sorted(serialized_retry),
        },
        "candidate_endpoint_control": {
            "terminal_reason": endpoint_stop,
            "first_failure_latched": endpoint_transaction.state.first_failed_premise
            == "newton_residual_limit",
        },
        "simultaneous_failure_control": simultaneous_stop,
        "healthy_control": {
            "accepted": healthy.accepted is not None,
            "accepted_stage_count": healthy_transaction.state.accepted_stage_count,
            "preaccept_observed_before_commit": healthy_preaccept_counts == [0],
        },
        "cross_member_rollback_control": {
            "restored_value": member.value,
            "incomplete_key_set_rejected": _incomplete_rollback_rejected(),
        },
        "all_controls_passed": (
            accepted_stop == "newton_residual_limit"
            and accepted_transaction.state.accepted_stage_count == 0
            and retry.retry.failed_stage_name == "rk4_k4"
            and retry_transaction.state == state_before
            and retry_transaction.causal_state == causal_before
            and preaccept_count == 0
            and endpoint_stop == "newton_residual_limit"
            and simultaneous_stop["reason"] == "nonpositive_lapse"
            and "newton_residual_limit" in simultaneous_stop["failures"]
            and healthy.accepted is not None
            and healthy_transaction.state.accepted_stage_count == 5
            and healthy_preaccept_counts == [0]
            and member.value == 3
        ),
    }


def _incomplete_rollback_rejected() -> bool:
    try:
        restore_common_event_snapshots({"a": object()}, {})
    except ValueError:
        return True
    return False


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
            initial.grid,
            spatial_order=order,
            ko_dissipation=1.0 / 64.0,
            residual_tolerance=1.0e-12,
            kinetic_condition_maximum=1.0e10,
        )(0.0, state)
        successor = Proto6GR0EvolutionOperator(
            initial.grid,
            spatial_order=order,
            ko_dissipation=1.0 / 64.0,
            raw_tolerance=1.0e-12,
            kinetic_condition_maximum=1.0e10,
        )(0.0, state)
        exact = all(
            np.array_equal(left, right)
            for left, right in (
                (original.du, successor.du),
                (original.dp, successor.dp),
                (original.dq, successor.dq),
            )
        )
        if not exact:
            raise RuntimeError("PROTO6 accepted source differs from original GR-0")
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
    plan_hash = _sha(loaded["paths"]["run_plan_config"])
    protocol_hash = loaded["raw"]["immutable_lineage"]["protocol_config_sha256"]
    frozen: list[dict[str, Any]] = []
    source_prechecks: list[dict[str, Any]] = []
    states: dict[tuple[str, str], list[tuple[int, EvolutionState, Any]]] = {}
    predecessor_inputs = lineage["predecessor_frozen_inputs"]
    expected_order = [
        (amplitude, method, count)
        for amplitude in ("5/2", "3")
        for method in ("RK4", "SSPRK3")
        for count in (1025, 2049, 4097)
    ]
    observed_order = [
        (item["amplitude"], item["method"], item["point_count"])
        for item in predecessor_inputs
    ]
    if observed_order != expected_order:
        raise ValueError("immutable HLT3 input order differs")
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
            raise ValueError("PROTO6 projected state differs from immutable HLT3")
        operator = Proto6GR0EvolutionOperator(
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
        passed = residual < threshold
        if not passed:
            raise ValueError("PROTO6 t0 accepted-state source precheck failed")
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
        record["plan_sha256"] = plan_hash
        record["protocol_sha256"] = protocol_hash
        record["source_stop_ownership"] = (
            "accepted_terminal_internal_source_only_trial_retry"
        )
        record["trajectory_advanced"] = False
        record.pop("expanded_run_config_sha256", None)
        record["expanded_run_config_sha256"] = _digest(record)
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
                raise ValueError("PROTO6 initial common event failed")
            common_events.append(
                {
                    "amplitude": amplitude,
                    "method": method,
                    "coordinate_time": 0.0,
                    "point_counts": list(assessed.point_counts),
                    "admission_passed": True,
                    "trajectory_advanced": False,
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
            raise ValueError("historical HLT4 namespace evidence differs")
        return dict(historical)
    records = []
    for key in ("calibration_output_root", "holdout_output_root"):
        relative = raw["namespace"][key]
        path = REPOSITORY / relative
        exists = path.exists()
        if exists and (not path.is_dir() or any(path.iterdir())):
            raise ValueError(f"HLT4 fresh namespace is nonempty: {relative}")
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
    equivalence = _operator_equivalence()
    inputs, source_prechecks, common_events = _freeze_inputs(loaded, lineage)
    namespace = _namespace_precondition(raw, historical_namespace_evidence)
    plan_hash = _sha(loaded["paths"]["run_plan_config"])
    input_manifest_hash = _digest(inputs)
    campaign_id = _digest(
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
        "classification": "pre_trajectory_PROTO6_GR0_source_retry_runtime_and_fresh_calibration_authorization",
        "generated_by": "scripts/reproduce_fgc_hlt4_mon4.py",
        "derivation_document": _rel(OWNER_DOCUMENT),
        "derivation_document_sha256": _sha(OWNER_DOCUMENT),
        "source_config_sha256": {
            _rel(config_path): _sha(config_path),
            _rel(loaded["paths"]["run_plan_config"]): plan_hash,
        },
        "predecessor_sha256": {
            EXPECTED_PATHS[name]: raw["immutable_lineage"][f"{name}_sha256"]
            for name in (
                "protocol_config",
                "protocol_freeze_config",
                "protocol_freeze_result",
                "predecessor_runtime_result",
                "static_input_result",
            )
        },
        "implementation_sha256": {_rel(path): _sha(path) for path in IMPLEMENTATION},
        "scope_bindings": dict(raw["scope"]),
        "gate_status": dict(raw["claims"]),
        "artifact_payload": {
            "immutable_lineage": public_lineage,
            "runtime_controls": runtime,
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
            "decision": "HLT4_implements_only_PROTO6_source_retry_ownership_and_authorizes_one_fresh_GR0_recalibration",
            "epistemic_boundary": {
                "PROTO5_calibration_history_disclosed": True,
                "PROTO6_trajectory_read": False,
                "PROTO6_trajectory_advanced": False,
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
    return _load_json_bytes(path.read_bytes(), "HLT4 result")


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
        raise ValueError("postlaunch HLT4 calibration directory is absent")
    if holdout_root.exists() and (
        not holdout_root.is_dir() or any(holdout_root.iterdir())
    ):
        raise ValueError("HLT4 does not authorize a nonempty holdout namespace")
    manifest_path = calibration_root / "manifest.json"
    manifest = _load_json_bytes(manifest_path.read_bytes(), "CAL3 launch manifest")
    frozen = observed["artifact_payload"]["frozen_run_plan"]
    expected = {
        "schema_version": 1,
        "runner_id": "FGC-1-CAL3-RUN1-RUNNER",
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
    for key, value in expected.items():
        if manifest.get(key) != value:
            raise ValueError(f"postlaunch CAL3 manifest field differs: {key}")
    commit = manifest.get("git_commit")
    if not isinstance(commit, str) or subprocess.run(
        ["git", "merge-base", "--is-ancestor", commit, "HEAD"],
        cwd=REPOSITORY,
        check=False,
        capture_output=True,
    ).returncode != 0:
        raise ValueError("postlaunch CAL3 Git checkpoint is invalid")


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
            historical_namespace_evidence=observed["artifact_payload"][
                "namespace_precondition"
            ],
        )
    else:
        expected = record(config_path)
    if observed != _serial(expected):
        raise ValueError("HLT4 output differs from a fresh reproduction")


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
