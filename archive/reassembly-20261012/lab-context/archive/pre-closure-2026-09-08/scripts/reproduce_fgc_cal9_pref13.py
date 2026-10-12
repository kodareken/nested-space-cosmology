#!/usr/bin/env python3
"""Reproduce CAL9/PREF13's outcome-neutral PROTO12 campaign diagnosis."""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile
import tomllib
from typing import Any, Mapping

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
import sys

sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from recursive_horizons.fgc.evolution.cal9_proto12_campaign_diagnosis import (  # noqa: E402
    EXPECTED_AMPLITUDES,
    EXPECTED_METHODS,
    EXPECTED_POINT_COUNTS,
    diagnose_proto12_campaign,
)
from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    array_content_sha256,
)


ARTIFACT_ID = "FGC-1-CAL9-PREF13"
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-cal9-pref13.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-cal9-pref13.json"
OWNER_DOCUMENT = REPOSITORY / "docs/fgc-cal9-pref13.md"
DIAGNOSIS_MODULE = (
    REPOSITORY
    / "src/recursive_horizons/fgc/evolution/cal9_proto12_campaign_diagnosis.py"
)
IMPLEMENTATION = (Path(__file__).resolve(), DIAGNOSIS_MODULE)
Q = Fraction


EXPECTED_CLAIMS = {
    "PROTO12_campaign_terminated_normally": True,
    "PROTO12_both_amplitudes_reached_late_evolved_common_events": True,
    "PROTO12_GR0_case_eligible": False,
    "SRC4_arithmetic_wall_cleared_for_both_amplitudes": True,
    "PROTO12_dual_amplitude_SSPRK3_radial_momentum_order_veto_localized": True,
    "RSP2_high_ladder_constraint_preflight_required": True,
    "PROTO13_frozen": False,
    "fresh_GR0_dynamic_calibration_completed": False,
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


def _rel(path: Path) -> str:
    return path.resolve().relative_to(REPOSITORY).as_posix()


def _sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _bytes_sha(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _canonical(value: Mapping[str, Any]) -> str:
    return json.dumps(
        value,
        sort_keys=True,
        indent=2,
        ensure_ascii=True,
        allow_nan=False,
    ) + "\n"


def _canonical_line(value: Mapping[str, Any]) -> bytes:
    return (
        json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
            allow_nan=False,
        )
        + "\n"
    ).encode("utf-8")


def _atomic_write(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        Path(temporary).replace(path)
    except BaseException:
        Path(temporary).unlink(missing_ok=True)
        raise


def _load_json(path: Path, artifact_id: str | None = None) -> dict[str, Any]:
    raw = path.read_bytes()
    value = json.loads(raw)
    if not isinstance(value, dict) or raw != _canonical(value).encode("utf-8"):
        raise ValueError(f"{_rel(path)} must contain canonical JSON")
    if artifact_id is not None and value.get("artifact_id") != artifact_id:
        raise ValueError(f"{_rel(path)} artifact identifier differs")
    return value


def _load_events(path: Path) -> tuple[list[dict[str, Any]], bytes]:
    raw = path.read_bytes()
    if not raw.endswith(b"\n"):
        raise ValueError("campaign event log must end in a newline")
    records: list[dict[str, Any]] = []
    rebuilt = bytearray()
    for line in raw.splitlines():
        value = json.loads(line)
        if not isinstance(value, dict):
            raise ValueError("campaign event records must be JSON objects")
        records.append(value)
        rebuilt.extend(_canonical_line(value))
    if bytes(rebuilt) != raw:
        raise ValueError("campaign event log is not canonical JSONL")
    return records, raw


def _git(*args: str, check: bool = True) -> subprocess.CompletedProcess[bytes]:
    return subprocess.run(
        ["git", *args],
        cwd=REPOSITORY,
        check=check,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )


def _validate_commit_lineage(commit: str, implementation: Mapping[str, str]) -> None:
    _git("cat-file", "-e", f"{commit}^{{commit}}")
    if _git("merge-base", "--is-ancestor", commit, "HEAD", check=False).returncode:
        raise ValueError("campaign authorization commit is not an ancestor of HEAD")
    if not implementation:
        raise ValueError("campaign implementation lineage is absent")
    for relative, expected in implementation.items():
        blob = _git("show", f"{commit}:{relative}").stdout
        if _bytes_sha(blob) != expected:
            raise ValueError(f"immutable campaign implementation differs: {relative}")


def _strict_keys(name: str, value: Mapping[str, Any], expected: set[str]) -> None:
    if set(value) != expected:
        raise ValueError(f"{name} keys differ")


def _as_float(value: str) -> float:
    return float(Q(value))


def _numeric_list(values: list[str]) -> list[float]:
    return [_as_float(value) for value in values]


def _configured_CFL_counts(value: Mapping[str, int]) -> dict[str, int]:
    expected_keys = {
        f"{method}_{count}"
        for method in EXPECTED_METHODS
        for count in EXPECTED_POINT_COUNTS
    }
    _strict_keys("terminal CFL counts", value, expected_keys)
    if any(
        isinstance(count, bool) or not isinstance(count, int) or count <= 0
        for count in value.values()
    ):
        raise ValueError("terminal CFL counts must be positive integers")
    return {key.replace("_", "-", 1): count for key, count in value.items()}


def _validate_amplitude_config(
    name: str,
    value: Mapping[str, Any],
    *,
    expected: Mapping[str, Any],
) -> None:
    _strict_keys(
        name,
        value,
        {
            "amplitude",
            "common_event_count",
            "last_fully_admitted_event_index",
            "last_fully_admitted_coordinate_time",
            "terminal_event_index",
            "terminal_coordinate_time",
            "terminal_classification",
            "terminal_reason",
            "serialized_source_retry_count",
            "primary_constraint_minimum_order",
            "comparator_previous_radial_momentum_order",
            "comparator_terminal_radial_momentum_order",
            "comparator_order_shortfall",
            "terminal_finest_primary_trapped_score",
            "terminal_finest_comparator_trapped_score",
            "terminal_radial_momentum_raw_owned_norms",
            "terminal_accumulated_roundoff_enclosures",
            "trailing_comparator_radial_momentum_orders",
            "trailing_event_indices",
            "terminal_CFL_retry_counts",
        },
    )
    for key, expected_value in expected.items():
        if value.get(key) != expected_value:
            raise ValueError(f"{name} {key} differs")
    if (
        value.get("terminal_classification")
        != "common_event_constraint_or_spectral_stop"
        or value.get("terminal_reason")
        != "common_event_constraint_or_spatial_spectral_admission"
        or value.get("serialized_source_retry_count") != 0
        or len(value.get("terminal_radial_momentum_raw_owned_norms", [])) != 3
        or len(value.get("terminal_accumulated_roundoff_enclosures", [])) != 3
        or len(value.get("trailing_comparator_radial_momentum_orders", [])) != 6
    ):
        raise ValueError(f"{name} numerical contract differs")
    _configured_CFL_counts(value["terminal_CFL_retry_counts"])


def load_config(path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    with path.open("rb") as handle:
        config = tomllib.load(handle)
    _strict_keys(
        "CAL9 config",
        config,
        {
            "schema_version",
            "artifact_id",
            "project_version",
            "metric_signature",
            "riemann_convention",
            "runtime_authorization_result",
            "run_plan_config",
            "campaign_manifest",
            "campaign_event_log",
            "campaign_checkpoint",
            "campaign_result",
            "scope",
            "immutable_campaign",
            "shared_veto",
            "amplitude_five_halves",
            "amplitude_three",
            "prospective_resolution_study",
            "proof_contract",
            "claims",
        },
    )
    scope = config["scope"]
    immutable = config["immutable_campaign"]
    shared = config["shared_veto"]
    revision = config["prospective_resolution_study"]
    _strict_keys(
        "CAL9 scope",
        scope,
        {
            "target_protocol",
            "role",
            "calibration_branch",
            "GR0_calibration_trajectory_read",
            "SGBL_trajectory_read",
            "FGCQR_trajectory_read",
            "collapse_or_trapped_outcome_classified",
            "mechanism_question_answered",
        },
    )
    _strict_keys(
        "CAL9 immutable campaign",
        immutable,
        {
            "authorization_commit",
            "campaign_id",
            "manifest_sha256",
            "event_log_sha256",
            "checkpoint_sha256",
            "campaign_result_sha256",
            "terminal_classification",
            "amplitudes_in_order",
            "event_log_line_count",
            "event_log_byte_count",
            "common_event_count",
            "serialized_source_retry_count",
            "selected_amplitude_absent",
            "holdout_execution_authorized",
        },
    )
    _strict_keys(
        "CAL9 shared veto",
        shared,
        {
            "method",
            "constraint_component",
            "gate",
            "minimum_constraint_finest_pair_order",
            "same_configured_gate_failed_for_both_amplitudes",
            "all_other_terminal_comparator_constraint_subgates_passed",
            "both_terminal_PROTO12_spatial_admissions_passed",
            "common_physical_cause_derived",
            "common_numerical_cause_derived",
            "preasymptotic_or_persistent_decided",
        },
    )
    if (
        config["schema_version"] != 1
        or config["artifact_id"] != ARTIFACT_ID
        or scope
        != {
            "target_protocol": "FGC-2-SF1-PROTO12",
            "role": "post_calibration_dual_amplitude_comparator_constraint_order_veto_diagnosis",
            "calibration_branch": "GR-0",
            "GR0_calibration_trajectory_read": True,
            "SGBL_trajectory_read": False,
            "FGCQR_trajectory_read": False,
            "collapse_or_trapped_outcome_classified": False,
            "mechanism_question_answered": False,
        }
        or immutable.get("authorization_commit")
        != "2f9bddbd4967e4a8998845513f345df3475a4bea"
        or immutable.get("campaign_id")
        != "32ac798e8bd7ccb984dd66262439d55ec01ad7c5cabfea785fe30d1fc60eff36"
        or immutable.get("terminal_classification")
        != "calibration_failed_no_eligible_GR0_case"
        or immutable.get("amplitudes_in_order") != list(EXPECTED_AMPLITUDES)
        or immutable.get("event_log_line_count") != 49
        or immutable.get("event_log_byte_count") != 21_096_450
        or immutable.get("common_event_count") != 49
        or immutable.get("serialized_source_retry_count") != 0
        or immutable.get("selected_amplitude_absent") is not True
        or immutable.get("holdout_execution_authorized") is not False
        or shared
        != {
            "method": "SSPRK3",
            "constraint_component": "radial_momentum",
            "gate": "finest_pair_order_passed",
            "minimum_constraint_finest_pair_order": "3/2",
            "same_configured_gate_failed_for_both_amplitudes": True,
            "all_other_terminal_comparator_constraint_subgates_passed": True,
            "both_terminal_PROTO12_spatial_admissions_passed": True,
            "common_physical_cause_derived": False,
            "common_numerical_cause_derived": False,
            "preasymptotic_or_persistent_decided": False,
        }
        or revision.get("next_gate") != "FGC-1-RSP2-FRZ1"
        or revision.get("purpose")
        != "prospectively_freeze_a_higher_ladder_late_event_constraint_convergence_test_before_any_PROTO13_design"
        or any(
            value is not True
            for key, value in revision.items()
            if key not in {"next_gate", "purpose"}
        )
        or any(value is not True for value in config["proof_contract"].values())
        or config["claims"] != EXPECTED_CLAIMS
    ):
        raise ValueError("CAL9/PREF13 configuration contract differs")

    _validate_amplitude_config(
        "amplitude five halves",
        config["amplitude_five_halves"],
        expected={
            "amplitude": "5/2",
            "common_event_count": 25,
            "last_fully_admitted_event_index": 23,
            "last_fully_admitted_coordinate_time": "23/16",
            "terminal_event_index": 24,
            "terminal_coordinate_time": "3/2",
            "trailing_event_indices": [19, 20, 21, 22, 23, 24],
        },
    )
    _validate_amplitude_config(
        "amplitude three",
        config["amplitude_three"],
        expected={
            "amplitude": "3",
            "common_event_count": 24,
            "last_fully_admitted_event_index": 22,
            "last_fully_admitted_coordinate_time": "11/8",
            "terminal_event_index": 23,
            "terminal_coordinate_time": "23/16",
            "trailing_event_indices": [18, 19, 20, 21, 22, 23],
        },
    )
    return config


def _validate_campaign(
    config: Mapping[str, Any],
) -> tuple[dict[str, Any], list[dict[str, Any]], bytes, dict[str, Any]]:
    immutable = config["immutable_campaign"]
    manifest_path = REPOSITORY / config["campaign_manifest"]
    event_path = REPOSITORY / config["campaign_event_log"]
    checkpoint_path = REPOSITORY / config["campaign_checkpoint"]
    result_path = REPOSITORY / config["campaign_result"]
    for path, expected in (
        (manifest_path, immutable["manifest_sha256"]),
        (event_path, immutable["event_log_sha256"]),
        (checkpoint_path, immutable["checkpoint_sha256"]),
        (result_path, immutable["campaign_result_sha256"]),
    ):
        if _sha(path) != expected:
            raise ValueError(f"immutable campaign hash differs: {_rel(path)}")
    manifest = _load_json(manifest_path)
    result = _load_json(result_path)
    events, event_bytes = _load_events(event_path)
    if (
        manifest.get("campaign_id") != immutable["campaign_id"]
        or manifest.get("git_commit") != immutable["authorization_commit"]
        or manifest.get("runner_id") != "FGC-1-CAL9-RUN1-RUNNER"
        or manifest.get("ordered_amplitudes") != list(EXPECTED_AMPLITUDES)
        or manifest.get("outcome_fields_present_at_creation") is not False
        or result.get("campaign_id") != immutable["campaign_id"]
        or result.get("runner_id") != "FGC-1-CAL9-RUN1-RUNNER"
        or result.get("classification") != immutable["terminal_classification"]
        or result.get("event_log_sha256") != immutable["event_log_sha256"]
        or result.get("selected_amplitude") is not None
        or result.get("holdout_execution_authorized") is not False
        or result.get("mechanism_question_answered") is not False
        or result.get("SGBL_outcome_read") is not False
        or result.get("FGCQR_outcome_read") is not False
        or result.get("retained_EFT_evolution_authorized") is not False
        or len(events) != immutable["event_log_line_count"]
        or len(event_bytes) != immutable["event_log_byte_count"]
    ):
        raise ValueError("immutable campaign manifest/result contract differs")

    authorization_path = REPOSITORY / config["runtime_authorization_result"]
    plan_path = REPOSITORY / config["run_plan_config"]
    authorization = _load_json(authorization_path, "FGC-1-HLT10-MON10")
    frozen_plan = authorization.get("artifact_payload", {}).get("frozen_run_plan", {})
    if (
        _sha(authorization_path) != manifest.get("authorization_sha256")
        or _sha(plan_path) != manifest.get("plan_sha256")
        or _rel(authorization_path) != manifest.get("authorization_path")
        or _rel(plan_path) != manifest.get("plan_path")
        or authorization.get("gate_status", {}).get(
            "PROTO12_fresh_GR0_dynamic_calibration_authorized"
        )
        is not True
        or authorization.get("gate_status", {}).get(
            "FGCQR_holdout_execution_authorized"
        )
        is not False
        or frozen_plan.get("campaign_id") != manifest.get("campaign_id")
        or frozen_plan.get("run_plan_sha256") != manifest.get("plan_sha256")
        or frozen_plan.get("input_manifest_sha256")
        != manifest.get("input_manifest_sha256")
        or frozen_plan.get("ordered_amplitudes") != list(EXPECTED_AMPLITUDES)
        or frozen_plan.get("point_counts") != list(EXPECTED_POINT_COUNTS)
    ):
        raise ValueError("campaign authorization or plan binding differs")
    _validate_commit_lineage(
        immutable["authorization_commit"],
        manifest.get("implementation_sha256", {}),
    )

    amplitude_three_CFL = _configured_CFL_counts(
        config["amplitude_three"]["terminal_CFL_retry_counts"]
    )
    with np.load(checkpoint_path, allow_pickle=False) as checkpoint:
        checkpoint_event_bytes = bytes(checkpoint["event_log_utf8"])
        metadata = json.loads(bytes(checkpoint["metadata_utf8"]).decode("utf-8"))
        members = metadata.get("members")
        if (
            checkpoint_event_bytes != event_bytes
            or metadata.get("schema_version") != 1
            or metadata.get("terminal") is not True
            or metadata.get("campaign_id") != immutable["campaign_id"]
            or metadata.get("amplitude") != "3"
            or metadata.get("amplitude_index") != 1
            or metadata.get("completed_common_event_index") != 23
            or metadata.get("consecutive_qualified_events") != 0
            or metadata.get("event_log", {}).get("sha256")
            != immutable["event_log_sha256"]
            or metadata.get("event_log", {}).get("line_count") != len(events)
            or metadata.get("event_log", {}).get("byte_count") != len(event_bytes)
            or metadata.get("completed_amplitude_records")
            != result.get("amplitude_records")
            or not isinstance(members, Mapping)
            or set(members) != set(amplitude_three_CFL)
        ):
            raise ValueError("terminal checkpoint metadata differs")
        final_record = result["amplitude_records"][-1]
        for member, expected_hash in final_record["member_final_state_hashes"].items():
            member_metadata = members[member]
            monitor = member_metadata.get("runtime_monitor_state", {})
            if (
                member_metadata.get("source_retry_count") != 0
                or member_metadata.get("CFL_retry_count")
                != amplitude_three_CFL[member]
                or member_metadata.get("time") != 1.4375
                or monitor.get("first_failed_premise") is not None
                or monitor.get("first_failed_time") is not None
                or monitor.get("first_failed_transaction_serial") is not None
                or monitor.get("last_accepted_time") != 1.4375
            ):
                raise ValueError(f"terminal checkpoint runtime differs: {member}")
            method, count = member.split("-")
            observed = array_content_sha256(
                checkpoint[f"{method}_{count}_u"],
                checkpoint[f"{method}_{count}_p"],
                checkpoint[f"{method}_{count}_q"],
            )
            if observed != expected_hash:
                raise ValueError(f"terminal checkpoint state differs: {member}")
    return manifest, events, event_bytes, result


def _validate_configured_amplitude(
    configured: Mapping[str, Any], diagnosis: Mapping[str, Any]
) -> None:
    expected_scalars = {
        "common_event_count": configured["common_event_count"],
        "last_fully_admitted_event_index": configured[
            "last_fully_admitted_event_index"
        ],
        "last_fully_admitted_coordinate_time": _as_float(
            configured["last_fully_admitted_coordinate_time"]
        ),
        "terminal_event_index": configured["terminal_event_index"],
        "terminal_coordinate_time": _as_float(configured["terminal_coordinate_time"]),
        "terminal_classification": configured["terminal_classification"],
        "terminal_reason": configured["terminal_reason"],
        "serialized_source_retry_count": configured["serialized_source_retry_count"],
        "primary_constraint_minimum_order": _as_float(
            configured["primary_constraint_minimum_order"]
        ),
        "comparator_previous_radial_momentum_order": _as_float(
            configured["comparator_previous_radial_momentum_order"]
        ),
        "comparator_terminal_radial_momentum_order": _as_float(
            configured["comparator_terminal_radial_momentum_order"]
        ),
        "comparator_order_shortfall": _as_float(
            configured["comparator_order_shortfall"]
        ),
        "terminal_finest_primary_trapped_score": _as_float(
            configured["terminal_finest_primary_trapped_score"]
        ),
        "terminal_finest_comparator_trapped_score": _as_float(
            configured["terminal_finest_comparator_trapped_score"]
        ),
    }
    if (
        diagnosis.get("amplitude") != configured["amplitude"]
        or any(diagnosis.get(key) != value for key, value in expected_scalars.items())
        or diagnosis.get("terminal_CFL_retry_count_by_member")
        != _configured_CFL_counts(configured["terminal_CFL_retry_counts"])
        or diagnosis.get("terminal_comparator_radial_momentum_raw_owned_norms")
        != _numeric_list(configured["terminal_radial_momentum_raw_owned_norms"])
        or diagnosis.get("terminal_comparator_accumulated_roundoff_enclosures")
        != _numeric_list(configured["terminal_accumulated_roundoff_enclosures"])
        or diagnosis.get("trailing_comparator_radial_momentum_orders")
        != _numeric_list(configured["trailing_comparator_radial_momentum_orders"])
        or diagnosis.get("trailing_event_indices")
        != configured["trailing_event_indices"]
        or diagnosis.get("comparator_only_component_below_minimum")
        != "radial_momentum"
        or diagnosis.get("both_PROTO12_spatial_spectral_admissions_passed")
        is not True
        or diagnosis.get("qualified_trapped_event_observed") is not False
    ):
        raise ValueError(
            f"configured CAL9 amplitude diagnosis differs: {configured['amplitude']}"
        )


def _validate_configured_diagnosis(
    config: Mapping[str, Any], diagnosis: Mapping[str, Any]
) -> None:
    immutable = config["immutable_campaign"]
    shared = diagnosis["shared_declared_veto"]
    configured_shared = config["shared_veto"]
    if (
        diagnosis["event_log_line_count"] != immutable["event_log_line_count"]
        or diagnosis["common_event_count"] != immutable["common_event_count"]
        or diagnosis["serialized_source_retry_count"]
        != immutable["serialized_source_retry_count"]
        or diagnosis["both_amplitudes_reached_late_evolved_common_events"]
        is not True
        or diagnosis["some_adaptive_CFL_step_reductions_occurred"] is not True
        or diagnosis["campaign_stop_was_not_CFL_retry_exhaustion"] is not True
        or diagnosis["SRC4_arithmetic_wall_cleared_for_both_amplitudes"] is not True
        or diagnosis["no_GR0_amplitude_eligible"] is not True
        or diagnosis["candidate_or_mechanism_outcome_read"] is not False
        or shared.get("method") != configured_shared["method"]
        or shared.get("constraint_component")
        != configured_shared["constraint_component"]
        or shared.get("gate") != configured_shared["gate"]
        or shared.get("frozen_minimum_order")
        != _as_float(configured_shared["minimum_constraint_finest_pair_order"])
        or any(
            shared.get(key) != configured_shared[key]
            for key in (
                "same_configured_gate_failed_for_both_amplitudes",
                "all_other_terminal_comparator_constraint_subgates_passed",
                "both_terminal_PROTO12_spatial_admissions_passed",
                "common_physical_cause_derived",
                "common_numerical_cause_derived",
                "preasymptotic_or_persistent_decided",
            )
        )
    ):
        raise ValueError("configured CAL9 campaign diagnosis differs")
    _validate_configured_amplitude(
        config["amplitude_five_halves"], diagnosis["amplitudes"]["5/2"]
    )
    _validate_configured_amplitude(
        config["amplitude_three"], diagnosis["amplitudes"]["3"]
    )


def reproduce(config_path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    config = load_config(config_path)
    manifest, events, event_bytes, campaign_result = _validate_campaign(config)
    minimum_order = _as_float(
        config["shared_veto"]["minimum_constraint_finest_pair_order"]
    )
    diagnosis = diagnose_proto12_campaign(
        events,
        campaign_result,
        minimum_constraint_order=minimum_order,
        trailing_event_count=6,
    )
    _validate_configured_diagnosis(config, diagnosis)
    runtime_authorization = REPOSITORY / config["runtime_authorization_result"]
    run_plan = REPOSITORY / config["run_plan_config"]
    record = {
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "project_version": config["project_version"],
        "classification": (
            "post_calibration_PROTO12_dual_amplitude_SSPRK3_radial_momentum_order_veto_diagnosis"
        ),
        "generated_by": "scripts/reproduce_fgc_cal9_pref13.py",
        "derivation_document": _rel(OWNER_DOCUMENT),
        "derivation_document_sha256": _sha(OWNER_DOCUMENT),
        "source_config_sha256": {_rel(config_path): _sha(config_path)},
        "predecessor_sha256": {
            _rel(runtime_authorization): _sha(runtime_authorization),
            _rel(run_plan): _sha(run_plan),
        },
        "implementation_sha256": {_rel(path): _sha(path) for path in IMPLEMENTATION},
        "scope_bindings": dict(config["scope"]),
        "gate_status": dict(EXPECTED_CLAIMS),
        "artifact_payload": {
            "immutable_campaign": {
                "authorization_commit": config["immutable_campaign"][
                    "authorization_commit"
                ],
                "campaign_id": manifest["campaign_id"],
                "manifest_sha256": config["immutable_campaign"]["manifest_sha256"],
                "event_log_sha256": _bytes_sha(event_bytes),
                "checkpoint_sha256": config["immutable_campaign"][
                    "checkpoint_sha256"
                ],
                "campaign_result_sha256": config["immutable_campaign"][
                    "campaign_result_sha256"
                ],
                "elapsed_wall_seconds": campaign_result["elapsed_wall_seconds"],
                "classification": campaign_result["classification"],
            },
            "campaign_diagnosis": diagnosis,
            "prospective_resolution_study": dict(
                config["prospective_resolution_study"]
            ),
            "epistemic_boundary": {
                "observed_results": [
                    "SRC4 advanced both frozen amplitudes to late evolved common events with zero serialized source retry",
                    "both amplitudes stopped when only the terminal SSPRK3 radial-momentum finest-pair constraint order fell below the unchanged three-halves gate",
                    "both terminal RK4 admissions and both terminal PROTO12 pairwise spatial-spectrum admissions passed",
                    "nonzero permitted adaptive CFL step reductions occurred and were not source retries or the terminal stop",
                ],
                "localized_numerical_question": (
                    "whether a prospectively frozen higher nested ladder clears or preserves the late SSPRK3 radial-momentum order veto"
                ),
                "not_derived": [
                    "one common physical or numerical cause for the two matching vetoes",
                    "an eligible GR-0 calibration case",
                    "a qualified trapped common event",
                    "SGB-L or FGC-QR dynamics",
                    "regulator activation or metric-null defocusing",
                    "a continuum branch or hyperbolicity loss",
                    "a retained-EFT or physical transition",
                ],
                "successor_status": (
                    "RSP2 prospective higher-ladder preflight is required; PROTO13 is not frozen"
                ),
            },
        },
        "nonclaims": {
            key: value for key, value in EXPECTED_CLAIMS.items() if value is False
        },
    }
    return record


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    record = reproduce(args.config.resolve())
    payload = _canonical(record).encode("utf-8")
    output = args.output.resolve()
    if args.check:
        if not output.is_file() or output.read_bytes() != payload:
            raise SystemExit(f"stored result differs: {output}")
        print(f"verified {output}")
        return
    _atomic_write(output, payload)
    print(f"wrote {output}")


if __name__ == "__main__":
    main()
