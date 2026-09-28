#!/usr/bin/env python3
"""Reduce CAL3 and replay PROTO6's terminal proposal-ownership boundary."""

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
from typing import Any, Mapping

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-cal3-pref5.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-cal3-pref5.json"
OWNER_DOCUMENT = REPOSITORY / "docs/fgc-cal3-pref5.md"

from scripts import run_fgc_gr0_calibration_v6 as campaign  # noqa: E402
from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    array_content_sha256,
)
from recursive_horizons.fgc.evolution.proto5_runtime import (  # noqa: E402
    CFLRetryRequired,
    GR0RuntimeStop,
)
from recursive_horizons.fgc.evolution.proto6_runtime import (  # noqa: E402
    attempt_proto6_step,
)


Q = Fraction
ARTIFACT_ID = "FGC-1-CAL3-PREF5"
PROJECT_VERSION = "0.11.0"
EXPECTED_PATHS = {
    "runtime_authorization_result": "results/fgc-1-hlt4-mon4.json",
    "run_plan_config": "configs/fgc/fgc-1-cal3-run1.toml",
    "campaign_manifest": "runs/fgc-2-sf1/proto6/calibration/manifest.json",
    "campaign_event_log": "runs/fgc-2-sf1/proto6/calibration/events.jsonl",
    "campaign_checkpoint": "runs/fgc-2-sf1/proto6/calibration/latest-checkpoint.npz",
    "campaign_result": "runs/fgc-2-sf1/proto6/calibration/campaign-result.json",
}
EXPECTED_SCOPE = {
    "target_protocol": "FGC-2-SF1-PROTO6",
    "role": "post_calibration_unaccepted_candidate_endpoint_source_stop_ownership_diagnosis",
    "calibration_branch": "GR-0",
    "GR0_calibration_trajectory_read": True,
    "SGBL_trajectory_read": False,
    "FGCQR_trajectory_read": False,
    "first_nonzero_common_event_completed": False,
    "collapse_or_trapped_outcome_classified": False,
    "mechanism_question_answered": False,
}
EXPECTED_CAMPAIGN = {
    "authorization_commit": "f8412dd202886dc60650acb887721be3fc62a737",
    "campaign_id": "62a0f63e0c381da444890bdb5b2fe7bfeeb2c9b8ba84b07a3f9f77669af6f1b8",
    "manifest_sha256": "8a9f92748e677c6a768ca61f38fd1390d6704d25ceb9421e6744e1a8b8490104",
    "event_log_sha256": "a4885ebe857f0b395012d0a09552cad775c670686c1d146c80b93afa97ec2bef",
    "checkpoint_sha256": "b5c0220ba5373cc182a852b071b724dc425af1c9914741e9ecd67eb4df34c9d2",
    "campaign_result_sha256": "fb50beb09821763964445f952a080215f7a8ba02a7e6dc70b3a0d4b76cfd3397",
    "terminal_classification": "calibration_failed_no_eligible_GR0_case",
    "amplitudes_in_order": ["5/2", "3"],
    "common_stop_reason": "newton_residual_limit",
    "event_log_line_count": 4,
    "raw_residual_tolerance": "1/1000000000000",
}
EXPECTED_REPLAY = {
    "amplitude": "3",
    "member": "RK4-4097",
    "target_common_event_time": "1/16",
    "accepted_attempt_count": 36,
    "CFL_retry_count": 37,
    "source_retry_count": 1,
    "last_accepted_time": "0.058593708693194996",
    "last_accepted_transaction_serial": 180,
    "accepted_state_residual": "5.370340537498194e-13",
    "retry_failed_stages": ["rk4_k2", "rk4_k3"],
    "terminal_failed_stages": ["rk4_k4", "candidate_endpoint"],
    "terminal_first_failure_serial": 183,
    "terminal_first_failure_time": "0.05940750925443598",
    "terminal_k4_residual": "1.3639071551646946e-12",
    "terminal_candidate_endpoint_residual": "1.3639179325880144e-12",
    "additional_retry_factor": "1/2",
    "additional_retry_must_pass": True,
    "maximum_relative_acceleration_correction": "1/1000000000000",
    "maximum_kinetic_condition": "10000000000",
    "failed_residual_point_index": 0,
    "failed_residual_radius": "1/32",
}
EXPECTED_REPAIR = {
    "physical_inputs_unchanged": True,
    "amplitude_order_unchanged": True,
    "raw_source_residual_tolerance_unchanged": True,
    "accepted_state_source_failure_remains_terminal": True,
    "every_source_only_failure_in_an_unaccepted_proposal_is_retryable": True,
    "candidate_endpoint_is_unaccepted_until_transaction_commit": True,
    "any_non_source_failure_prevents_retry": True,
    "retry_factor_inherited": "1/2",
    "maximum_retries_inherited": 32,
    "minimum_step_size_inherited": "1/1073741824",
    "terminal_and_retry_exhaustion_evidence_must_name_member_stage_and_residual": True,
    "constraint_spectral_trapped_boundary_and_physical_rules_unchanged": True,
    "fresh_output_namespace_and_new_protocol_version_required": True,
}
EXPECTED_CLAIMS = {
    "PROTO6_unaccepted_candidate_endpoint_stop_ownership_contract_obstructed": True,
    "PROTO7_general_transaction_revision_required": True,
    "fresh_GR0_dynamic_calibration_completed": False,
    "PROTO6_resolved_holdout_manifest_authorized": False,
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
EXPECTED_PROOF = {
    "raw_campaign_files_must_match_frozen_hashes",
    "manifest_must_bind_HLT4_plan_implementation_and_commit",
    "event_log_must_contain_two_t0_events_and_two_exact_retry_records",
    "terminal_checkpoint_must_bind_atomic_cross_member_rollback",
    "both_amplitudes_must_share_the_universal_source_stop",
    "diagnostic_member_must_reconstruct_from_the_HLT4_t0_hash",
    "replay_retry_must_match_the_public_event_record",
    "last_accepted_state_must_pass_the_original_raw_gate",
    "terminal_proposal_must_fail_only_at_unaccepted_k4_and_candidate_endpoint",
    "one_further_halving_must_pass_the_unchanged_complete_transaction",
    "kinetic_and_relative_correction_margins_must_be_serialized",
    "diagnosis_must_not_read_SGBL_or_FGCQR",
    "diagnosis_must_not_answer_the_mechanism_question",
    "canonical_hash_bound_result_required",
}
IMPLEMENTATION = (Path(__file__).resolve(),)
RK4_CALL_NAMES = (
    "accepted_state_source_precheck",
    "rk4_k1",
    "rk4_k2",
    "rk4_k3",
    "rk4_k4",
    "candidate_endpoint",
)


def _pairs(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    answer: dict[str, Any] = {}
    for key, value in pairs:
        if key in answer:
            raise ValueError(f"duplicate JSON key {key}")
        answer[key] = value
    return answer


def _serial(value: Any) -> Any:
    if isinstance(value, Fraction):
        return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"
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
            raise ValueError("nonfinite CAL3 evidence")
        return value
    if value is None or isinstance(value, (str, bool, int)):
        return value
    raise TypeError(f"unsupported CAL3 evidence type {type(value).__name__}")


def _canonical(value: Any) -> str:
    return json.dumps(_serial(value), sort_keys=True, indent=2, allow_nan=False) + "\n"


def _canonical_line(value: Any) -> str:
    return json.dumps(
        _serial(value), sort_keys=True, separators=(",", ":"), allow_nan=False
    ) + "\n"


def _sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _rel(path: Path) -> str:
    return path.resolve().relative_to(REPOSITORY).as_posix()


def _path(name: str, value: Any, expected: str) -> Path:
    if value != expected:
        raise ValueError(f"{name} must name {expected}")
    path = (REPOSITORY / expected).resolve()
    if not path.is_file() or _rel(path) != expected:
        raise ValueError(f"{name} must be the canonical source path")
    return path


def _load_json(path: Path, name: str) -> dict[str, Any]:
    source = path.read_text(encoding="utf-8")
    value = json.loads(source, object_pairs_hook=_pairs)
    if not isinstance(value, dict) or source != _canonical(value):
        raise ValueError(f"{name} must be canonical unique-key JSON")
    return value


def _load_events(path: Path) -> list[dict[str, Any]]:
    answer = []
    for line_number, line in enumerate(
        path.read_text(encoding="utf-8").splitlines(keepends=True), start=1
    ):
        value = json.loads(line, object_pairs_hook=_pairs)
        if not isinstance(value, dict) or line != _canonical_line(value):
            raise ValueError(f"CAL3 event line {line_number} is noncanonical")
        answer.append(value)
    return answer


def _is_ancestor(commit: str) -> bool:
    return subprocess.run(
        ["git", "merge-base", "--is-ancestor", commit, "HEAD"],
        cwd=REPOSITORY,
        check=False,
    ).returncode == 0


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
        "immutable_campaign",
        "diagnostic_replay",
        "prospective_repair",
        "proof_contract",
        "claims",
    }
    if set(raw) != expected_root:
        raise ValueError("CAL3 root keys differ")
    if (
        raw["schema_version"] != 1
        or raw["artifact_id"] != ARTIFACT_ID
        or raw["project_version"] != PROJECT_VERSION
        or raw["metric_signature"] != "-+++"
        or raw["riemann_convention"] != "plus_partial_mu_gamma_nu"
        or raw["scope"] != EXPECTED_SCOPE
        or raw["immutable_campaign"] != EXPECTED_CAMPAIGN
        or raw["diagnostic_replay"] != EXPECTED_REPLAY
        or raw["prospective_repair"] != EXPECTED_REPAIR
        or raw["claims"] != EXPECTED_CLAIMS
        or set(raw["proof_contract"]) != EXPECTED_PROOF
        or any(value is not True for value in raw["proof_contract"].values())
    ):
        raise ValueError("CAL3 frozen contract differs")
    paths = {
        name: _path(name, raw[name], expected)
        for name, expected in EXPECTED_PATHS.items()
    }
    for path_key, hash_key in (
        ("campaign_manifest", "manifest_sha256"),
        ("campaign_event_log", "event_log_sha256"),
        ("campaign_checkpoint", "checkpoint_sha256"),
        ("campaign_result", "campaign_result_sha256"),
    ):
        if _sha(paths[path_key]) != raw["immutable_campaign"][hash_key]:
            raise ValueError(f"CAL3 raw campaign hash differs: {path_key}")
    if not _is_ancestor(raw["immutable_campaign"]["authorization_commit"]):
        raise ValueError("CAL3 authorization commit is not an ancestor of HEAD")
    plan, authorization = campaign._validate_authorization(
        paths["run_plan_config"],
        paths["runtime_authorization_result"],
        require_fresh_namespace=False,
    )
    return {"raw": raw, "paths": paths, "plan": plan, "authorization": authorization}


def _validate_campaign(loaded: Mapping[str, Any]) -> dict[str, Any]:
    raw = loaded["raw"]
    paths = loaded["paths"]
    manifest = _load_json(paths["campaign_manifest"], "CAL3 manifest")
    result = _load_json(paths["campaign_result"], "CAL3 campaign result")
    events = _load_events(paths["campaign_event_log"])
    frozen = loaded["authorization"]["artifact_payload"]["frozen_run_plan"]
    if (
        manifest.get("campaign_id") != raw["immutable_campaign"]["campaign_id"]
        or manifest.get("git_commit") != raw["immutable_campaign"]["authorization_commit"]
        or manifest.get("plan_sha256") != _sha(paths["run_plan_config"])
        or manifest.get("authorization_sha256")
        != _sha(paths["runtime_authorization_result"])
        or manifest.get("input_manifest_sha256") != frozen["input_manifest_sha256"]
        or manifest.get("implementation_sha256")
        != loaded["authorization"]["implementation_sha256"]
        or manifest.get("outcome_fields_present_at_creation") is not False
    ):
        raise ValueError("CAL3 manifest does not bind HLT4 and its launch commit")
    if (
        result.get("campaign_id") != manifest["campaign_id"]
        or result.get("classification")
        != raw["immutable_campaign"]["terminal_classification"]
        or result.get("selected_amplitude") is not None
        or result.get("event_log_sha256") != _sha(paths["campaign_event_log"])
        or result.get("holdout_execution_authorized") is not False
        or result.get("mechanism_question_answered") is not False
        or result.get("SGBL_outcome_read") is not False
        or result.get("FGCQR_outcome_read") is not False
    ):
        raise ValueError("CAL3 terminal result or nonclaim boundary differs")
    if len(events) != raw["immutable_campaign"]["event_log_line_count"]:
        raise ValueError("CAL3 event-log length differs")
    retry_events: dict[str, dict[str, Any]] = {}
    for index, amplitude in enumerate(raw["immutable_campaign"]["amplitudes_in_order"]):
        initial = events[2 * index]
        retry = events[2 * index + 1]
        assessment = initial.get("assessment", {})
        if (
            initial.get("event_type") != "common_event"
            or initial.get("amplitude") != amplitude
            or initial.get("event_index") != 0
            or assessment.get("coordinate_time") != 0.0
            or assessment.get("primary_common_event", {}).get("admission_passed") is not True
            or assessment.get("comparator_common_event", {}).get("admission_passed") is not True
            or retry.get("event_type") != "rejected_internal_source_trial"
            or retry.get("amplitude") != amplitude
            or retry.get("member") != raw["diagnostic_replay"]["member"]
            or [item.get("stage_name") for item in retry.get("failed_stages", [])]
            != raw["diagnostic_replay"]["retry_failed_stages"]
            or retry.get("accepted_state_bitwise_preserved") is not True
            or retry.get("time_advanced") is not False
            or retry.get("transaction_advanced") is not False
            or retry.get("causal_debit_advanced") is not False
        ):
            raise ValueError("CAL3 initial/retry event sequence differs")
        retry_events[amplitude] = retry
    records = result.get("amplitude_records")
    if not isinstance(records, list) or len(records) != 2:
        raise ValueError("CAL3 terminal result must contain two amplitude records")
    for amplitude, record in zip(
        raw["immutable_campaign"]["amplitudes_in_order"], records, strict=True
    ):
        stop = record.get("stop", {})
        if (
            record.get("amplitude") != amplitude
            or record.get("eligible") is not False
            or record.get("last_completed_common_event_index") != 0
            or record.get("last_completed_coordinate_time") != 0.0
            or stop.get("classification") != "scientific_universal_runtime_stop"
            or stop.get("reason") != raw["immutable_campaign"]["common_stop_reason"]
            or stop.get("failures") != [raw["immutable_campaign"]["common_stop_reason"]]
            or stop.get("target_common_event_index") != 1
            or stop.get("cross_member_state_rolled_back") is not True
        ):
            raise ValueError("CAL3 amplitude stop differs")
    checkpoint, _members, checkpoint_log = campaign._restore_checkpoint(
        paths["campaign_checkpoint"],
        plan=loaded["plan"],
        authorization=loaded["authorization"],
    )
    if (
        checkpoint.get("terminal") is not True
        or checkpoint.get("campaign_id") != manifest["campaign_id"]
        or checkpoint.get("completed_common_event_index") != 0
        or checkpoint.get("completed_amplitude_records") != records
        or any(item.get("time") != 0.0 for item in checkpoint.get("members", {}).values())
        or sha256(checkpoint_log).hexdigest() != _sha(paths["campaign_event_log"])
    ):
        raise ValueError("CAL3 terminal checkpoint does not prove common-event rollback")
    return {
        "manifest": manifest,
        "result": result,
        "events": events,
        "retry_events": retry_events,
        "checkpoint": checkpoint,
    }


def _call_payload(name: str, time: float, rhs: Any) -> dict[str, Any]:
    diagnostics = rhs.diagnostics
    return {
        "stage_name": name,
        "coordinate_time": time,
        "source_residual_infinity": diagnostics["source_residual_infinity"],
        "source_raw_gate_passed": diagnostics["source_raw_gate_passed"],
        "source_refinement_iterations": diagnostics["source_refinement_iterations"],
        "source_roundoff_floor_reached": diagnostics["source_roundoff_floor_reached"],
        "relative_acceleration_correction": diagnostics[
            "source_relative_acceleration_correction"
        ],
        "kinetic_condition_infinity": diagnostics["kinetic_condition_infinity"],
        "minimum_lapse": diagnostics["minimum_lapse"],
        "maximum_residual_point_index": diagnostics[
            "source_maximum_residual_point_index"
        ],
        "maximum_residual_row_index": diagnostics[
            "source_maximum_residual_row_index"
        ],
        "maximum_residual_radius": diagnostics["source_maximum_residual_radius"],
    }


def _named_calls(calls: list[tuple[float, Any]]) -> list[dict[str, Any]]:
    if len(calls) != len(RK4_CALL_NAMES):
        raise ValueError("CAL3 RK4 attempt did not expose six source calls")
    return [
        _call_payload(name, time, rhs)
        for name, (time, rhs) in zip(RK4_CALL_NAMES, calls, strict=True)
    ]


def _attempt(
    member: Any,
    step_size: float,
) -> tuple[Any, list[dict[str, Any]]]:
    calls: list[tuple[float, Any]] = []
    original = member.operator

    def tracing_rhs(time: float, state: Any) -> Any:
        rhs = original(time, state)
        calls.append((time, rhs))
        return rhs

    def preview(proposal: Any) -> Any:
        return member.tracers.preview_advance(
            old_state=member.state,
            new_state=proposal.candidate_state,
            coordinates=member.initial.grid.coordinates,
            step_size=proposal.final_time - member.time,
        )

    try:
        result = attempt_proto6_step(
            method=member.integrator_id,
            time=member.time,
            step_size=step_size,
            state=member.state,
            rhs=tracing_rhs,
            projector=member.projector,
            transaction=member.transaction,
            previous_step_index=member.step_index,
            previous_transaction_serial=member.transaction_serial,
            preaccept=preview,
        )
    except (CFLRetryRequired, GR0RuntimeStop) as error:
        setattr(error, "cal3_calls", _named_calls(calls))
        raise
    return result, _named_calls(calls)


def _replay(loaded: Mapping[str, Any], public_retry: Mapping[str, Any]) -> dict[str, Any]:
    raw = loaded["raw"]
    plan = loaded["plan"]
    replay = raw["diagnostic_replay"]
    member = campaign._build_members(
        plan, loaded["authorization"], replay["amplitude"]
    )[replay["member"]]
    initial_hash = array_content_sha256(member.state.u, member.state.p, member.state.q)
    frozen_input = next(
        item
        for item in loaded["authorization"]["artifact_payload"]["frozen_run_inputs"]
        if item["amplitude"] == replay["amplitude"]
        and item["method"] == member.method_label
        and item["point_count"] == member.point_count
    )
    if initial_hash != frozen_input["projected_state_sha256"]:
        raise ValueError("CAL3 reconstructed member differs from HLT4")
    target = float(Q(replay["target_common_event_time"]))
    retry_factor = float(Q(plan["numerics"]["retry_factor"]))
    minimum_step = float(Q(plan["numerics"]["minimum_step_size"]))
    accepted_attempts = 0
    cfl_retries = 0
    source_retries = 0
    replay_retry: dict[str, Any] | None = None
    terminal: dict[str, Any] | None = None
    smaller_control: dict[str, Any] | None = None
    tolerance = 4096.0 * np.finfo(np.float64).eps
    while member.time < target - tolerance:
        previous_speed = max(
            member.transaction.causal_state.previous_speed_upper,
            np.finfo(np.float64).tiny,
        )
        proposed_size = min(
            target - member.time,
            member.transaction.cfl_maximum
            * member.initial.grid.spacing
            / previous_speed,
        )
        while True:
            if proposed_size < minimum_step:
                raise ValueError("CAL3 replay fell below the frozen minimum step")
            snapshot = member.snapshot()
            try:
                attempt, calls = _attempt(member, proposed_size)
            except CFLRetryRequired:
                member.restore(snapshot)
                cfl_retries += 1
                proposed_size *= retry_factor
                continue
            except GR0RuntimeStop as stop:
                calls = getattr(stop, "cal3_calls")
                terminal = {
                    "reason": stop.reason,
                    "failures": list(stop.failures),
                    "last_accepted_time": snapshot["time"],
                    "last_accepted_state_sha256": array_content_sha256(
                        snapshot["state"].u,
                        snapshot["state"].p,
                        snapshot["state"].q,
                    ),
                    "last_accepted_transaction_serial": snapshot[
                        "transaction_serial"
                    ],
                    "terminal_trial_step_size": proposed_size,
                    "latched_monitor_state": asdict(stop.state),
                    "calls": calls,
                }
                member.restore(snapshot)
                smaller_size = proposed_size * float(Q(replay["additional_retry_factor"]))
                try:
                    smaller, smaller_calls = _attempt(member, smaller_size)
                except (CFLRetryRequired, GR0RuntimeStop) as smaller_stop:
                    smaller_control = {
                        "accepted": False,
                        "error_type": type(smaller_stop).__name__,
                        "calls": getattr(smaller_stop, "cal3_calls"),
                        "step_size": smaller_size,
                    }
                else:
                    smaller_control = {
                        "accepted": smaller.accepted is not None,
                        "retry_returned": smaller.retry is not None,
                        "step_size": smaller_size,
                        "calls": smaller_calls,
                    }
                member.restore(snapshot)
                break
            if attempt.retry is not None:
                member.restore(snapshot)
                source_retries += 1
                replay_retry = {
                    **asdict(attempt.retry),
                    "member": member.key,
                    "amplitude": member.amplitude,
                    "calls": calls,
                }
                proposed_size *= retry_factor
                continue
            if attempt.accepted is None:
                raise RuntimeError("CAL3 replay produced no accepted step")
            member.tracers.commit_advance(*attempt.preaccept_payload)
            member.state = attempt.accepted.state
            member.time = attempt.accepted.time
            member.step_index = attempt.accepted.step_index
            member.transaction_serial = attempt.accepted.transaction_serial
            accepted_attempts += 1
            break
        if terminal is not None:
            break
    if terminal is None or replay_retry is None or smaller_control is None:
        raise ValueError("CAL3 replay did not reproduce retry then terminal stop")
    public_projection = {
        key: replay_retry[key]
        for key in (
            "initial_time",
            "proposed_step_size",
            "failed_stage_name",
            "failed_stage_time",
            "source_residual_infinity",
            "source_residual_maximum",
            "failed_transaction_serial",
            "failed_stages",
            "accepted_state_bitwise_preserved",
            "time_advanced",
            "transaction_advanced",
            "causal_debit_advanced",
            "member",
            "amplitude",
        )
    }
    public_expected = {key: public_retry[key] for key in public_projection}
    if _serial(public_projection) != _serial(public_expected):
        raise ValueError("CAL3 replay retry differs from the immutable event log")
    failed_terminal = [
        item["stage_name"]
        for item in terminal["calls"]
        if item["source_raw_gate_passed"] is False
    ]
    if (
        accepted_attempts != replay["accepted_attempt_count"]
        or cfl_retries != replay["CFL_retry_count"]
        or source_retries != replay["source_retry_count"]
        or terminal["last_accepted_time"] != float(replay["last_accepted_time"])
        or terminal["last_accepted_transaction_serial"]
        != replay["last_accepted_transaction_serial"]
        or terminal["calls"][0]["source_residual_infinity"]
        != float(replay["accepted_state_residual"])
        or terminal["calls"][0]["source_raw_gate_passed"] is not True
        or failed_terminal != replay["terminal_failed_stages"]
        or terminal["latched_monitor_state"]["first_failed_transaction_serial"]
        != replay["terminal_first_failure_serial"]
        or terminal["latched_monitor_state"]["first_failed_time"]
        != float(replay["terminal_first_failure_time"])
        or terminal["calls"][4]["source_residual_infinity"]
        != float(replay["terminal_k4_residual"])
        or terminal["calls"][5]["source_residual_infinity"]
        != float(replay["terminal_candidate_endpoint_residual"])
        or smaller_control["accepted"] is not replay["additional_retry_must_pass"]
        or smaller_control.get("retry_returned") is not False
        or any(
            item["source_raw_gate_passed"] is not True
            for item in smaller_control["calls"]
        )
    ):
        raise ValueError("CAL3 deterministic replay differs from the frozen diagnosis")
    maximum_condition = max(
        item["kinetic_condition_infinity"] for item in terminal["calls"]
    )
    maximum_correction = max(
        item["relative_acceleration_correction"] for item in terminal["calls"]
    )
    failed_calls = [
        item for item in terminal["calls"] if not item["source_raw_gate_passed"]
    ]
    if (
        maximum_condition >= float(Q(replay["maximum_kinetic_condition"]))
        or maximum_correction >= float(Q(replay["maximum_relative_acceleration_correction"]))
        or any(
            item["maximum_residual_point_index"]
            != replay["failed_residual_point_index"]
            or item["maximum_residual_radius"]
            != float(Q(replay["failed_residual_radius"]))
            for item in failed_calls
        )
    ):
        raise ValueError("CAL3 replay health or localization margin differs")
    return {
        "initial_state_sha256": initial_hash,
        "accepted_attempt_count": accepted_attempts,
        "CFL_retry_count": cfl_retries,
        "source_retry_count": source_retries,
        "public_retry_record_reproduced_exactly": True,
        "retry_evidence": replay_retry,
        "terminal_evidence": terminal,
        "one_further_halving_control": smaller_control,
        "maximum_terminal_kinetic_condition": maximum_condition,
        "maximum_terminal_relative_acceleration_correction": maximum_correction,
        "all_terminal_failures_are_source_only": terminal["failures"]
        == ["newton_residual_limit"],
        "last_accepted_state_passes_raw_gate": terminal["calls"][0][
            "source_raw_gate_passed"
        ],
        "failed_terminal_stages": failed_terminal,
    }


def record(config_path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    loaded = load_config(config_path)
    raw = loaded["raw"]
    campaign_evidence = _validate_campaign(loaded)
    replay = _replay(
        loaded,
        campaign_evidence["retry_events"][raw["diagnostic_replay"]["amplitude"]],
    )
    gates = dict(raw["claims"])
    return {
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "project_version": PROJECT_VERSION,
        "classification": "post_PROTO6_GR0_unaccepted_candidate_endpoint_source_stop_ownership_diagnosis",
        "generated_by": "scripts/reproduce_fgc_cal3_pref5.py",
        "derivation_document": _rel(OWNER_DOCUMENT),
        "derivation_document_sha256": _sha(OWNER_DOCUMENT),
        "source_config_sha256": {_rel(config_path): _sha(config_path)},
        "metric_signature": "-+++",
        "riemann_convention": "plus_partial_mu_gamma_nu",
        "scope_bindings": raw["scope"],
        "gate_status": gates,
        "artifact_payload": {
            "immutable_campaign": {
                "campaign_id": campaign_evidence["manifest"]["campaign_id"],
                "authorization_commit": campaign_evidence["manifest"]["git_commit"],
                "manifest_sha256": _sha(loaded["paths"]["campaign_manifest"]),
                "event_log_sha256": _sha(loaded["paths"]["campaign_event_log"]),
                "checkpoint_sha256": _sha(loaded["paths"]["campaign_checkpoint"]),
                "campaign_result_sha256": _sha(loaded["paths"]["campaign_result"]),
                "terminal_classification": campaign_evidence["result"]["classification"],
                "event_log_line_count": len(campaign_evidence["events"]),
                "amplitude_stop_records": campaign_evidence["result"]["amplitude_records"],
                "terminal_checkpoint_common_event_index": campaign_evidence[
                    "checkpoint"
                ]["completed_common_event_index"],
                "cross_member_rollback_proved": True,
            },
            "deterministic_terminal_replay": replay,
            "decision": {
                "PROTO6_unaccepted_candidate_endpoint_stop_ownership_contract_obstructed": True,
                "reason": "the_last_accepted_state_passes_but_a_source_only_unaccepted_proposal_is_terminal_only_because_its_candidate_endpoint_also_misses_the_raw_gate",
                "prospective_repair": raw["prospective_repair"],
                "new_protocol_required": "FGC-2-SF1-PROTO7",
            },
            "epistemic_boundary": {
                "GR0_calibration_completed": False,
                "GR0_amplitude_selected": False,
                "first_nonzero_common_event_completed": False,
                "dynamic_trapped_sphere_classified": False,
                "SGBL_or_FGCQR_trajectory_read": False,
                "mechanism_question_answered": False,
                "candidate_or_general_gradient_route_rejected": False,
            },
        },
        "implementation_sha256": {_rel(path): _sha(path) for path in IMPLEMENTATION},
        "nonclaims": {
            "fresh_GR0_calibration_completed": False,
            "trapped_sphere_observed": False,
            "SGBL_control_completed": False,
            "FGCQR_holdout_evolved": False,
            "regulator_activation_observed": False,
            "positive_Raychaudhuri_margin_observed": False,
            "FGCQR_mechanism_rejected": False,
            "general_gradient_route_rejected": False,
            "retained_EFT_validity": False,
            "physical_transition": False,
        },
    }


def _validate_stored(result: Mapping[str, Any], loaded: Mapping[str, Any]) -> None:
    raw = loaded["raw"]
    payload = result.get("artifact_payload", {})
    immutable = payload.get("immutable_campaign", {})
    replay = payload.get("deterministic_terminal_replay", {})
    terminal = replay.get("terminal_evidence", {})
    smaller = replay.get("one_further_halving_control", {})
    if (
        result.get("artifact_id") != ARTIFACT_ID
        or result.get("classification")
        != "post_PROTO6_GR0_unaccepted_candidate_endpoint_source_stop_ownership_diagnosis"
        or result.get("scope_bindings") != raw["scope"]
        or result.get("gate_status") != raw["claims"]
        or immutable.get("campaign_id") != raw["immutable_campaign"]["campaign_id"]
        or immutable.get("manifest_sha256") != raw["immutable_campaign"]["manifest_sha256"]
        or immutable.get("event_log_sha256") != raw["immutable_campaign"]["event_log_sha256"]
        or immutable.get("checkpoint_sha256") != raw["immutable_campaign"]["checkpoint_sha256"]
        or immutable.get("campaign_result_sha256")
        != raw["immutable_campaign"]["campaign_result_sha256"]
        or immutable.get("cross_member_rollback_proved") is not True
        or replay.get("accepted_attempt_count")
        != raw["diagnostic_replay"]["accepted_attempt_count"]
        or replay.get("CFL_retry_count") != raw["diagnostic_replay"]["CFL_retry_count"]
        or replay.get("source_retry_count") != raw["diagnostic_replay"]["source_retry_count"]
        or replay.get("public_retry_record_reproduced_exactly") is not True
        or replay.get("last_accepted_state_passes_raw_gate") is not True
        or replay.get("failed_terminal_stages")
        != raw["diagnostic_replay"]["terminal_failed_stages"]
        or replay.get("all_terminal_failures_are_source_only") is not True
        or terminal.get("last_accepted_time")
        != float(raw["diagnostic_replay"]["last_accepted_time"])
        or terminal.get("calls", [])[0].get("source_residual_infinity")
        != float(raw["diagnostic_replay"]["accepted_state_residual"])
        or smaller.get("accepted") is not True
        or smaller.get("retry_returned") is not False
        or any(not item.get("source_raw_gate_passed") for item in smaller.get("calls", []))
        or payload.get("decision", {}).get("new_protocol_required")
        != "FGC-2-SF1-PROTO7"
        or payload.get("epistemic_boundary", {}).get("mechanism_question_answered")
        is not False
    ):
        raise ValueError("stored CAL3 result violates its frozen semantic boundary")
    expected_impl = {_rel(path): _sha(path) for path in IMPLEMENTATION}
    if (
        result.get("implementation_sha256") != expected_impl
        or result.get("derivation_document") != _rel(OWNER_DOCUMENT)
        or result.get("derivation_document_sha256") != _sha(OWNER_DOCUMENT)
        or result.get("source_config_sha256")
        != {_rel(DEFAULT_CONFIG): _sha(DEFAULT_CONFIG)}
    ):
        raise ValueError("stored CAL3 implementation or derivation hash differs")


def load_canonical_result(path: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    return _load_json(path, ARTIFACT_ID)


def verify_canonical(
    path: Path = DEFAULT_OUTPUT,
    config_path: Path = DEFAULT_CONFIG,
    *,
    replay: bool = False,
) -> None:
    loaded = load_config(config_path)
    _validate_campaign(loaded)
    actual = load_canonical_result(path)
    _validate_stored(actual, loaded)
    if replay and actual != record(config_path):
        raise ValueError(f"{_rel(path)} differs from full CAL3 replay")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--replay", action="store_true")
    args = parser.parse_args()
    config = args.config.resolve()
    output = args.output.resolve()
    if args.check:
        verify_canonical(output, config, replay=args.replay)
        print(f"verified {_rel(output)} (full_replay={str(args.replay).lower()})")
        return
    payload = record(config)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(_canonical(payload), encoding="utf-8")
    print(
        f"wrote {_rel(output)} "
        "(PROTO6 candidate-endpoint ownership obstructed=true; mechanism answered=false)"
    )


if __name__ == "__main__":
    main()
