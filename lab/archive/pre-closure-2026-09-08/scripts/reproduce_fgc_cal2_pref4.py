#!/usr/bin/env python3
"""Reduce and diagnose the immutable PROTO5 GR-0 CAL2 campaign."""

from __future__ import annotations

import argparse
from dataclasses import asdict, is_dataclass
from fractions import Fraction
from hashlib import sha256
import json
from math import isfinite
from pathlib import Path
import re
import subprocess
import sys
import tomllib
from typing import Any, Mapping

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-cal2-pref4.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-cal2-pref4.json"
OWNER_DOCUMENT = REPOSITORY / "docs/fgc-cal2-pref4.md"

from recursive_horizons.fgc.evolution.cal2_source_diagnosis import (  # noqa: E402
    DiagnosticGR0EvolutionOperator,
)
from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    propose_step,
)
from scripts import run_fgc_gr0_calibration as campaign  # noqa: E402


Q = Fraction
ARTIFACT_ID = "FGC-1-CAL2-PREF4"
PROJECT_VERSION = "0.11.0"
EXPECTED_PATHS = {
    "runtime_authorization_result": "results/fgc-1-hlt3-mon3.json",
    "run_plan_config": "configs/fgc/fgc-1-cal2-run1.toml",
    "campaign_manifest": "runs/fgc-2-sf1/proto5/calibration/manifest.json",
    "campaign_event_log": "runs/fgc-2-sf1/proto5/calibration/events.jsonl",
    "campaign_checkpoint": "runs/fgc-2-sf1/proto5/calibration/latest-checkpoint.npz",
    "campaign_result": "runs/fgc-2-sf1/proto5/calibration/campaign-result.json",
}
EXPECTED_SCOPE = {
    "target_protocol": "FGC-2-SF1-PROTO5",
    "role": "post_calibration_source_stop_ownership_diagnosis",
    "calibration_branch": "GR-0",
    "GR0_calibration_trajectory_read": True,
    "SGBL_trajectory_read": False,
    "FGCQR_trajectory_read": False,
    "first_nonzero_common_event_completed": False,
    "collapse_or_trapped_outcome_classified": False,
    "mechanism_question_answered": False,
}
EXPECTED_CAMPAIGN = {
    "authorization_commit": "fb7d1dc2499438ac639e6a1686cd3f5561fa93bc",
    "campaign_id": "9df7046b47a044df920bfe192e6b63fbda7953ab97d9174850927afec10177be",
    "manifest_sha256": "f47c9c54c9ddccb9c1b97807512c3ed1ae93f6692030e0d583ac8b71c6ba55d0",
    "event_log_sha256": "022ca3629134eb710df931aeef03ccba91c288856473696bf7a6ddeaf9e1b2d1",
    "checkpoint_sha256": "8607800ada752c7571b00e975399b8c038c145fa8635b3ffaedf4a010c1e7f12",
    "campaign_result_sha256": "9cbd2cdd44a85bb7c3d65b2c13c7b6f192e4d38375f8e4efc642208717be2296",
    "terminal_classification": "calibration_failed_no_eligible_GR0_case",
    "amplitudes_in_order": ["5/2", "3"],
    "common_stop_reason": "newton_residual_limit",
    "raw_residual_tolerance": "1/1000000000000",
    "first_residual_observed": "1.0511642392237476e-12",
    "second_residual_observed": "1.0511640449698101e-12",
}
EXPECTED_REPLAY = {
    "checkpoint_amplitude": "3",
    "member": "RK4-4097",
    "target_common_event_time": "1/16",
    "proposal_step_factors": ["1", "1/2", "1/4", "1/8"],
    "accepted_state_source_gate_must_pass": True,
    "unit_step_factor_must_miss_raw_gate": True,
    "one_half_step_factor_must_pass_every_raw_stage_gate": True,
    "one_eighth_step_factor_must_pass_every_raw_stage_gate": True,
    "condition_number_must_remain_below_frozen_limit": True,
    "failed_residual_must_localize_to_first_annular_node": True,
    "unresolved_relative_acceleration_correction_max": "1/1000000000000",
}
EXPECTED_REPAIR = {
    "physical_inputs_unchanged": True,
    "amplitude_order_unchanged": True,
    "raw_source_residual_tolerance_unchanged": True,
    "accepted_state_source_failure_remains_terminal": True,
    "only_internal_unaccepted_stage_source_failure_becomes_retryable": True,
    "retry_factor_inherited": "1/2",
    "maximum_retries_inherited": 32,
    "minimum_step_size_inherited": "1/1073741824",
    "candidate_action_health_stops_unchanged": True,
    "constraint_spectral_trapped_and_boundary_rules_unchanged": True,
    "fresh_output_namespace_and_new_protocol_version_required": True,
}
EXPECTED_PROOF = {
    "raw_campaign_files_must_match_the_frozen_hashes",
    "manifest_must_bind_HLT3_plan_inputs_implementations_and_commit",
    "event_log_must_contain_only_two_passing_t0_common_events",
    "checkpoint_must_bind_the_terminal_mixed_member_state_without_claiming_a_common_event",
    "both_amplitudes_must_stop_on_the_same_universal_source_gate",
    "accepted_fine_state_must_pass_the_original_raw_source_gate",
    "full_unaccepted_proposal_must_reproduce_a_raw_source_gate_miss",
    "one_half_unaccepted_proposal_must_pass_the_unchanged_raw_source_gate",
    "one_eighth_unaccepted_proposal_must_pass_the_unchanged_raw_source_gate",
    "kinetic_condition_and_acceleration_correction_margins_must_be_serialized",
    "diagnosis_must_not_read_SGBL_or_FGCQR",
    "diagnosis_must_not_answer_the_mechanism_question",
    "canonical_hash_bound_result_required",
}
EXPECTED_CLAIMS = {
    "PROTO5_GR0_source_stop_ownership_contract_obstructed": True,
    "PROTO6_premise_revision_required": True,
    "fresh_GR0_dynamic_calibration_completed": False,
    "PROTO5_resolved_holdout_manifest_authorized": False,
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
        "src/recursive_horizons/fgc/evolution/cal2_source_diagnosis.py",
        "scripts/reproduce_fgc_cal2_pref4.py",
    )
)
RESIDUAL_PATTERN = re.compile(
    r"observed=([0-9.eE+-]+), refinements=([0-9]+)"
)


def _pairs(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    answer: dict[str, Any] = {}
    for key, value in pairs:
        if key in answer:
            raise ValueError(f"duplicate JSON key: {key}")
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
            raise ValueError("nonfinite value cannot enter canonical evidence")
        return value
    if value is None or isinstance(value, (str, bool, int)):
        return value
    raise TypeError(f"unsupported evidence value {type(value).__name__}")


def _canonical(value: Any) -> str:
    return json.dumps(_serial(value), indent=2, sort_keys=True, allow_nan=False) + "\n"


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


def _load_event_log(path: Path) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for line_number, line in enumerate(
        path.read_text(encoding="utf-8").splitlines(keepends=True), start=1
    ):
        value = json.loads(line, object_pairs_hook=_pairs)
        if not isinstance(value, dict) or line != _canonical_line(value):
            raise ValueError(f"campaign event line {line_number} is not canonical")
        records.append(value)
    return records


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
        raise ValueError("CAL2 root keys differ")
    if (
        raw["schema_version"] != 1
        or raw["artifact_id"] != ARTIFACT_ID
        or raw["project_version"] != PROJECT_VERSION
        or raw["metric_signature"] != "-+++"
        or raw["riemann_convention"] != "plus_partial_mu_gamma_nu"
    ):
        raise ValueError("CAL2 identity or convention differs")
    if raw["scope"] != EXPECTED_SCOPE:
        raise ValueError("CAL2 scope differs")
    if raw["immutable_campaign"] != EXPECTED_CAMPAIGN:
        raise ValueError("CAL2 immutable campaign binding differs")
    if raw["diagnostic_replay"] != EXPECTED_REPLAY:
        raise ValueError("CAL2 diagnostic replay contract differs")
    if raw["prospective_repair"] != EXPECTED_REPAIR:
        raise ValueError("CAL2 prospective repair differs")
    if set(raw["proof_contract"]) != EXPECTED_PROOF or any(
        value is not True for value in raw["proof_contract"].values()
    ):
        raise ValueError("CAL2 proof contract differs")
    if raw["claims"] != EXPECTED_CLAIMS:
        raise ValueError("CAL2 claims differ")
    paths = {
        name: _path(name, raw[name], expected)
        for name, expected in EXPECTED_PATHS.items()
    }
    hash_keys = {
        "campaign_manifest": "manifest_sha256",
        "campaign_event_log": "event_log_sha256",
        "campaign_checkpoint": "checkpoint_sha256",
        "campaign_result": "campaign_result_sha256",
    }
    for path_key, hash_key in hash_keys.items():
        if _sha(paths[path_key]) != raw["immutable_campaign"][hash_key]:
            raise ValueError(f"CAL2 raw campaign hash differs: {path_key}")
    if not _is_ancestor(raw["immutable_campaign"]["authorization_commit"]):
        raise ValueError("CAL2 authorization commit is not an ancestor of HEAD")
    plan, authorization = campaign._validate_authorization(
        paths["run_plan_config"],
        paths["runtime_authorization_result"],
        require_fresh_namespace=False,
    )
    return {
        "raw": raw,
        "paths": paths,
        "plan": plan,
        "authorization": authorization,
    }


def _validate_campaign(
    loaded: Mapping[str, Any],
) -> tuple[
    dict[str, Any],
    list[dict[str, Any]],
    dict[str, Any],
    dict[str, Any],
    Mapping[str, Any],
]:
    raw = loaded["raw"]
    paths = loaded["paths"]
    plan = loaded["plan"]
    authorization = loaded["authorization"]
    manifest = _load_json(paths["campaign_manifest"], "CAL2 manifest")
    result = _load_json(paths["campaign_result"], "CAL2 campaign result")
    events = _load_event_log(paths["campaign_event_log"])
    frozen_plan = authorization["artifact_payload"]["frozen_run_plan"]
    if (
        manifest.get("campaign_id") != raw["immutable_campaign"]["campaign_id"]
        or manifest.get("git_commit")
        != raw["immutable_campaign"]["authorization_commit"]
        or manifest.get("plan_sha256") != _sha(paths["run_plan_config"])
        or manifest.get("authorization_sha256")
        != _sha(paths["runtime_authorization_result"])
        or manifest.get("input_manifest_sha256")
        != frozen_plan["input_manifest_sha256"]
        or manifest.get("implementation_sha256")
        != authorization["implementation_sha256"]
        or manifest.get("outcome_fields_present_at_creation") is not False
    ):
        raise ValueError("CAL2 manifest does not bind the authorized campaign")
    if len(events) != 2:
        raise ValueError("CAL2 event log must contain exactly two initial events")
    for expected_amplitude, event in zip(
        raw["immutable_campaign"]["amplitudes_in_order"], events, strict=True
    ):
        assessment = event.get("assessment", {})
        if (
            event.get("event_type") != "common_event"
            or event.get("amplitude") != expected_amplitude
            or event.get("event_index") != 0
            or assessment.get("coordinate_time") != 0.0
            or assessment.get("primary_common_event", {}).get("admission_passed")
            is not True
            or assessment.get("comparator_common_event", {}).get(
                "admission_passed"
            )
            is not True
            or assessment.get("qualified_trapped_common_event") is not False
        ):
            raise ValueError("CAL2 initial common-event record differs")
    if (
        result.get("campaign_id") != manifest["campaign_id"]
        or result.get("classification")
        != raw["immutable_campaign"]["terminal_classification"]
        or result.get("selected_amplitude") is not None
        or result.get("holdout_execution_authorized") is not False
        or result.get("mechanism_question_answered") is not False
        or result.get("SGBL_outcome_read") is not False
        or result.get("FGCQR_outcome_read") is not False
        or result.get("event_log_sha256") != _sha(paths["campaign_event_log"])
    ):
        raise ValueError("CAL2 terminal result or nonclaim boundary differs")
    amplitude_records = result.get("amplitude_records")
    if not isinstance(amplitude_records, list) or len(amplitude_records) != 2:
        raise ValueError("CAL2 must contain two amplitude records")
    for expected_amplitude, record in zip(
        raw["immutable_campaign"]["amplitudes_in_order"],
        amplitude_records,
        strict=True,
    ):
        stop = record.get("stop", {})
        if (
            record.get("amplitude") != expected_amplitude
            or record.get("eligible") is not False
            or record.get("last_completed_common_event_index") != 0
            or record.get("last_completed_coordinate_time") != 0.0
            or stop.get("reason")
            != raw["immutable_campaign"]["common_stop_reason"]
            or stop.get("target_common_event_index") != 1
            or stop.get("last_fully_completed_common_event_index") != 0
        ):
            raise ValueError("CAL2 amplitude stop differs")
    checkpoint, members, checkpoint_log = campaign._restore_checkpoint(
        paths["campaign_checkpoint"],
        plan=plan,
        authorization=authorization,
    )
    if (
        checkpoint.get("terminal") is not True
        or checkpoint.get("campaign_id") != manifest["campaign_id"]
        or checkpoint.get("amplitude")
        != raw["diagnostic_replay"]["checkpoint_amplitude"]
        or checkpoint.get("completed_common_event_index") != 0
        or checkpoint.get("completed_amplitude_records") != amplitude_records
        or sha256(checkpoint_log).hexdigest() != _sha(paths["campaign_event_log"])
    ):
        raise ValueError("CAL2 terminal checkpoint binding differs")
    return manifest, events, result, checkpoint, members


def _stop_evidence(
    raw: Mapping[str, Any], result: Mapping[str, Any]
) -> list[dict[str, Any]]:
    expected_values = (
        raw["immutable_campaign"]["first_residual_observed"],
        raw["immutable_campaign"]["second_residual_observed"],
    )
    evidence = []
    for record, expected_text in zip(
        result["amplitude_records"], expected_values, strict=True
    ):
        stop = record["stop"]
        match = RESIDUAL_PATTERN.search(stop["error_message"])
        if match is None:
            raise ValueError("CAL2 source stop lacks residual/refinement evidence")
        observed_text, refinements_text = match.groups()
        if observed_text != expected_text:
            raise ValueError("CAL2 observed residual differs from frozen value")
        observed = float(observed_text)
        tolerance = float(Q(raw["immutable_campaign"]["raw_residual_tolerance"]))
        evidence.append(
            {
                "amplitude": record["amplitude"],
                "reason": stop["reason"],
                "observed_residual_infinity": observed,
                "raw_residual_tolerance": tolerance,
                "relative_excess_over_raw_tolerance": observed / tolerance - 1.0,
                "refinement_iterations": int(refinements_text),
                "last_fully_completed_common_event_index": 0,
            }
        )
    relative_difference = abs(
        evidence[0]["observed_residual_infinity"]
        - evidence[1]["observed_residual_infinity"]
    ) / max(item["observed_residual_infinity"] for item in evidence)
    for item in evidence:
        item["cross_amplitude_relative_residual_difference"] = relative_difference
    return evidence


def _stage_payload(record: Any) -> dict[str, Any]:
    diagnostics = record.rhs.diagnostics
    return {
        "stage_name": record.stage_name,
        "coordinate_time": record.time,
        "raw_source_residual_infinity": diagnostics[
            "source_residual_infinity"
        ],
        "raw_source_gate_passed": diagnostics["source_raw_gate_passed"],
        "refinement_iterations": diagnostics["source_refinement_iterations"],
        "roundoff_floor_reached": diagnostics["source_roundoff_floor_reached"],
        "relative_acceleration_correction": diagnostics[
            "source_relative_acceleration_correction"
        ],
        "kinetic_condition_infinity": diagnostics[
            "kinetic_condition_infinity"
        ],
        "maximum_residual_point_index": diagnostics[
            "source_maximum_residual_point_index"
        ],
        "maximum_residual_row_index": diagnostics[
            "source_maximum_residual_row_index"
        ],
        "maximum_residual_radius": diagnostics[
            "source_maximum_residual_radius"
        ],
    }


def _diagnostic_replay(
    loaded: Mapping[str, Any], members: Mapping[str, Any]
) -> dict[str, Any]:
    raw = loaded["raw"]
    plan = loaded["plan"]
    member = members.get(raw["diagnostic_replay"]["member"])
    if member is None:
        raise ValueError("CAL2 diagnostic member is absent from checkpoint")
    if member.amplitude != raw["diagnostic_replay"]["checkpoint_amplitude"]:
        raise ValueError("CAL2 diagnostic member amplitude differs")
    thresholds = plan["universal_thresholds"]
    numerics = plan["numerics"]
    operator = DiagnosticGR0EvolutionOperator(
        member.initial.grid,
        spatial_order=member.spatial_order,
        ko_dissipation=float(Q(numerics["ko_dissipation"])),
        raw_tolerance=float(Q(thresholds["source_residual_infinity_max"])),
        kinetic_condition_maximum=float(
            Q(thresholds["kinetic_condition_number_max"])
        ),
    )
    target = float(Q(raw["diagnostic_replay"]["target_common_event_time"]))
    previous_speed = member.transaction.causal_state.previous_speed_upper
    base_step = min(
        target - member.time,
        member.transaction.cfl_maximum
        * member.initial.grid.spacing
        / previous_speed,
    )
    proposals = []
    for factor_text in raw["diagnostic_replay"]["proposal_step_factors"]:
        factor = float(Q(factor_text))
        proposal = propose_step(
            method=member.integrator_id,
            time=member.time,
            step_size=base_step * factor,
            state=member.state,
            rhs=operator,
            projector=member.projector,
        )
        stages = [_stage_payload(stage) for stage in proposal.stages]
        proposals.append(
            {
                "step_factor": factor_text,
                "proposed_step_size": base_step * factor,
                "all_stage_raw_source_gates_passed": all(
                    stage["raw_source_gate_passed"] for stage in stages
                ),
                "maximum_raw_source_residual_infinity": max(
                    stage["raw_source_residual_infinity"] for stage in stages
                ),
                "maximum_relative_acceleration_correction": max(
                    stage["relative_acceleration_correction"] for stage in stages
                ),
                "maximum_kinetic_condition_infinity": max(
                    stage["kinetic_condition_infinity"] for stage in stages
                ),
                "stages": stages,
            }
        )
    accepted_stage = proposals[0]["stages"][0]
    failed_stages = [
        stage
        for stage in proposals[0]["stages"]
        if not stage["raw_source_gate_passed"]
    ]
    if not accepted_stage["raw_source_gate_passed"]:
        raise ValueError("CAL2 last accepted state misses the original source gate")
    if proposals[0]["all_stage_raw_source_gates_passed"]:
        raise ValueError("CAL2 unit proposal did not reproduce a source miss")
    if not proposals[1]["all_stage_raw_source_gates_passed"]:
        raise ValueError("CAL2 one-half proposal did not restore the raw gate")
    if not proposals[-1]["all_stage_raw_source_gates_passed"]:
        raise ValueError("CAL2 one-eighth proposal did not restore the raw gate")
    spacing = member.initial.grid.spacing
    if not failed_stages or any(
        abs(stage["maximum_residual_radius"] - spacing)
        > 4096.0 * np.finfo(np.float64).eps * max(1.0, spacing)
        for stage in failed_stages
    ):
        raise ValueError("CAL2 failed residual is not at the first annular node")
    condition_limit = float(Q(thresholds["kinetic_condition_number_max"]))
    if any(
        item["maximum_kinetic_condition_infinity"] >= condition_limit
        for item in proposals
    ):
        raise ValueError("CAL2 replay reaches the kinetic condition stop")
    correction_limit = float(
        Q(
            raw["diagnostic_replay"][
                "unresolved_relative_acceleration_correction_max"
            ]
        )
    )
    if any(
        stage["relative_acceleration_correction"] >= correction_limit
        for item in proposals
        for stage in item["stages"]
    ):
        raise ValueError("CAL2 replay correction exceeds the declared bound")
    return {
        "checkpoint_member": member.key,
        "last_accepted_coordinate_time": member.time,
        "target_common_event_time": target,
        "grid_spacing": spacing,
        "previous_accepted_speed_upper": previous_speed,
        "base_proposed_step_size": base_step,
        "accepted_state_source_evidence": accepted_stage,
        "proposals": proposals,
        "first_retry_factor_that_passes_every_stage": next(
            item["step_factor"]
            for item in proposals
            if item["all_stage_raw_source_gates_passed"]
        ),
        "failed_unit_proposal_stage_names": [
            stage["stage_name"] for stage in failed_stages
        ],
        "failure_localizes_to_first_annular_node": True,
        "kinetic_condition_limit": condition_limit,
        "relative_acceleration_correction_limit": correction_limit,
    }


def record(config_path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    loaded = load_config(config_path)
    raw = loaded["raw"]
    manifest, events, result, checkpoint, members = _validate_campaign(loaded)
    stops = _stop_evidence(raw, result)
    replay = _diagnostic_replay(loaded, members)
    gates = dict(raw["claims"])
    gates.update(
        {
            "both_amplitudes_reached_only_t0_common_event": True,
            "both_amplitudes_share_threshold_adjacent_source_stop": True,
            "last_accepted_fine_state_passes_original_source_gate": True,
            "unaccepted_full_proposal_misses_original_source_gate": True,
            "one_half_proposal_passes_unchanged_source_gate": True,
            "one_eighth_proposal_passes_unchanged_source_gate": True,
            "kinetic_condition_stop_reached": False,
        }
    )
    return {
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "project_version": PROJECT_VERSION,
        "classification": "post_calibration_numerical_trial_stop_ownership_obstruction_not_mechanism_result",
        "generated_by": "scripts/reproduce_fgc_cal2_pref4.py",
        "derivation_document": _rel(OWNER_DOCUMENT),
        "derivation_document_sha256": _sha(OWNER_DOCUMENT),
        "source_config_sha256": {_rel(config_path): _sha(config_path)},
        "predecessor_sha256": {
            _rel(loaded["paths"]["runtime_authorization_result"]): _sha(
                loaded["paths"]["runtime_authorization_result"]
            ),
            _rel(loaded["paths"]["run_plan_config"]): _sha(
                loaded["paths"]["run_plan_config"]
            ),
        },
        "raw_campaign_sha256": {
            _rel(loaded["paths"][name]): _sha(loaded["paths"][name])
            for name in (
                "campaign_manifest",
                "campaign_event_log",
                "campaign_checkpoint",
                "campaign_result",
            )
        },
        "implementation_sha256": {
            _rel(path): _sha(path) for path in IMPLEMENTATION
        },
        "certificate_contract": sorted(raw["proof_contract"]),
        "scope_bindings": {
            "target_protocol": raw["scope"]["target_protocol"],
            "calibration_branch": raw["scope"]["calibration_branch"],
            "GR0_calibration_trajectory_read": True,
            "SGBL_trajectory_read": False,
            "FGCQR_trajectory_read": False,
            "first_nonzero_common_event_completed": False,
        },
        "gate_status": gates,
        "artifact_payload": {
            "immutable_campaign": {
                "manifest": manifest,
                "terminal_classification": result["classification"],
                "event_log_line_count": len(events),
                "terminal_checkpoint": {
                    "amplitude": checkpoint["amplitude"],
                    "completed_common_event_index": checkpoint[
                        "completed_common_event_index"
                    ],
                    "terminal": checkpoint["terminal"],
                    "member_times": {
                        key: value["time"]
                        for key, value in checkpoint["members"].items()
                    },
                },
            },
            "amplitude_stop_evidence": stops,
            "diagnostic_replay": replay,
            "decision": {
                "PROTO5_GR0_source_stop_ownership_contract_obstructed": True,
                "reason": "an_unaccepted_internal_RK_stage_missed_an_absolute_binary64_source_residual_gate_while_the_last_accepted_state_passed_and_a_smaller_trial_restored_the_same_gate",
                "prospective_repair": raw["prospective_repair"],
                "new_protocol_required": "FGC-2-SF1-PROTO6",
            },
            "epistemic_boundary": {
                "GR0_calibration_completed": False,
                "GR0_amplitude_selected": False,
                "dynamic_trapped_sphere_classified": False,
                "SGBL_or_FGCQR_trajectory_read": False,
                "mechanism_question_answered": False,
                "candidate_or_general_gradient_route_rejected": False,
            },
        },
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


def load_canonical_result(path: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    return _load_json(path, ARTIFACT_ID)


def verify_canonical(
    path: Path = DEFAULT_OUTPUT, config_path: Path = DEFAULT_CONFIG
) -> None:
    expected = record(config_path)
    actual = load_canonical_result(path)
    if actual != expected:
        raise ValueError(f"{_rel(path)} differs from canonical CAL2 reproduction")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    config = args.config.resolve()
    output = args.output.resolve()
    if args.check:
        verify_canonical(output, config)
        print(f"verified {_rel(output)}")
        return
    payload = record(config)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(_canonical(payload), encoding="utf-8")
    print(
        f"wrote {_rel(output)} "
        "(PROTO5 source-stop ownership obstructed=true; mechanism answered=false)"
    )


if __name__ == "__main__":
    main()
