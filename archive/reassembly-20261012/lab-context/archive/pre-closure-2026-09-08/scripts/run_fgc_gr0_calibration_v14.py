#!/usr/bin/env python3
"""Run the HLT12-authorized PROTO14 GR-0 restart calibration.

PROTO14 is a fresh continuation from the six immutable ``t=23/16`` inputs.
It inherits PROTO13's GR-0 equations, method-owned spatial/constraint gates,
and trapped diagnostic.  It replaces only the retired sampled temporal veto:
every post-restart accepted macro step is admitted by TDG6 and only the fine
four-quarter shadow path commits.  The retained 64-event quantity is written
as a diagnostic only.  This runner cannot construct SGB-L or FGC-QR state.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict
from hashlib import sha256
import json
import os
from pathlib import Path
import sys
import time as wall_time
import tomllib
from typing import Any, Mapping

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts import reproduce_fgc_hlt12_mon12 as hlt12  # noqa: E402
from scripts import run_fgc_gr0_calibration as inherited  # noqa: E402
from scripts import run_fgc_gr0_calibration_v6 as proto6_runner  # noqa: E402
from scripts import run_fgc_gr0_calibration_v13 as proto13_runner  # noqa: E402
from recursive_horizons.fgc.evolution.boundary_domain import CausalBudgetState  # noqa: E402
from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    EvolutionState,
    array_content_sha256,
)
from recursive_horizons.fgc.evolution.proto5_runtime import GR0RuntimeMonitorState  # noqa: E402
from recursive_horizons.fgc.evolution.proto6_runtime import restore_common_event_snapshots  # noqa: E402
from recursive_horizons.fgc.evolution.proto7_runtime import Proto7TerminalStop  # noqa: E402
from recursive_horizons.fgc.evolution.proto14_runtime import (  # noqa: E402
    Proto14SourceRetryExhausted,
    Proto14RunMember,
    proto14_checkpoint_extension,
    restore_proto14_checkpoint_extension,
)
from recursive_horizons.fgc.evolution.tdg6_temporal_admission_runtime import (  # noqa: E402
    TDG6TemporalRetryExhausted,
)


DEFAULT_PLAN = REPOSITORY / "configs/fgc/fgc-1-pro14-run1.toml"
DEFAULT_AUTHORIZATION = REPOSITORY / "results/fgc-1-hlt12-mon12.json"
RUNNER_ID = "FGC-1-CAL11-RUN1-RUNNER"
RESTART_TIME = 23.0 / 16.0
_MEMBER_KEYS = (
    "RK4-2049", "RK4-4097", "RK4-8193",
    "SSPRK3-4097", "SSPRK3-8193", "SSPRK3-16385",
)


def _relative(path: Path) -> str:
    return path.resolve().relative_to(REPOSITORY).as_posix()


def _validate_run_plan(path: Path) -> dict[str, Any]:
    """Validate the narrow PROTO14 execution contract before a launch."""

    with path.open("rb") as handle:
        plan = tomllib.load(handle)
    required = {
        "schema_version", "artifact_id", "project_version", "metric_signature",
        "riemann_convention", "protocol_config", "protocol_freeze_result",
        "runtime_authorization_result", "inherited_run_plan",
        "inherited_runtime_authorization", "tdg6_runtime_result", "scope",
        "numerical_runtime", "method_owned_ladders", "restart",
        "temporal_admission", "schedule", "assessment", "provenance", "claims",
    }
    if set(plan) != required:
        raise ValueError("PROTO14 run plan root differs")
    if (
        plan["schema_version"] != 1
        or plan["artifact_id"] != "FGC-1-CAL11-RUN1-PLAN"
        or plan["project_version"] != "0.11.0"
        or plan["metric_signature"] != "-+++"
        or plan["riemann_convention"] != "plus_partial_mu_gamma_nu"
        or plan["protocol_config"] != "configs/fgc/fgc-2-sf1-protocol-v14.toml"
        or plan["protocol_freeze_result"] != "results/fgc-1-pro14-frz1.json"
        or plan["runtime_authorization_result"] != "results/fgc-1-hlt12-mon12.json"
        or plan["scope"] != {
            "target_protocol": "FGC-2-SF1-PROTO14", "branch": "GR-0",
            "role": "TDG6_all_step_temporal_admission_restart_calibration",
            "physical_equations": "unredefined_GR0_specialization_of_REF1",
            "amplitude": "3", "candidate_action_fields_or_health_stops_forbidden": True,
            "retained_EFT_interpretation": False,
        }
        or plan["numerical_runtime"] != hlt12.EXPECTED_RUNTIME
        or plan["method_owned_ladders"] != {
            "primary_method": "RK4", "primary_point_counts": [2049, 4097, 8193],
            "comparator_method": "SSPRK3", "comparator_point_counts": [4097, 8193, 16385],
            "common_physical_point_count": 2049,
            "each_method_must_pass_its_own_complete_admission": True,
            "raw_unequal_finest_grid_collocation_forbidden": True,
        }
        or plan["restart"]["coordinate_time"] != "23/16"
        or plan["schedule"] != {
            "first_new_common_event_coordinate_time": "3/2", "final_coordinate_time": "32",
            "common_output_interval": "1/16", "checkpoint_interval": "1/4",
            "minimum_consecutive_trapped_common_events": 8,
            "trapped_sign_margin_over_combined_error_factor": 4,
            "temporal_history_minimum_samples": 64,
        }
        or plan["temporal_admission"]["instrument_artifact_id"] != "FGC-1-TDG6-IMP2"
        or plan["temporal_admission"]["level_step_counts"] != [1, 2, 4]
        or plan["temporal_admission"]["channel_count"] != 18
        or plan["temporal_admission"]["minimum_observed_order"] != "3/2"
        or plan["temporal_admission"]["exact_outward_order_test"] != "8*U12^2<=L01^2"
        or plan["temporal_admission"]["every_accepted_post_restart_macro_step_must_pass"] is not True
        or plan["temporal_admission"]["only_fine_path_may_commit"] is not True
        or plan["temporal_admission"]["sampled_64_history_diagnostic_only"] is not True
        or plan["provenance"]["output_root"] != "runs/fgc-2-sf1/proto14/calibration"
        or plan["provenance"]["terminal_checkpoint_must_precede_result_publication"] is not True
        or plan["provenance"][
            "resume_must_idempotently_publish_or_verify_the_checkpoint_owned_terminal_result"
        ] is not True
        or any(value is not False for value in plan["claims"].values())
    ):
        raise ValueError("PROTO14 run plan violates the frozen contract")
    return plan


def _validate_authorization(
    plan_path: Path, authorization_path: Path, *, require_fresh_namespace: bool
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Validate HLT12 and this exact runner plan before a production launch."""

    plan = _validate_run_plan(plan_path)
    runtime_environment = hlt12.validate_numerical_runtime(plan["numerical_runtime"])
    authorization = hlt12.load_canonical_result(authorization_path)
    gates = authorization.get("gate_status", {})
    payload = authorization.get("artifact_payload", {})
    frozen = payload.get("frozen_run_plan", {})
    if (
        authorization.get("artifact_id") != hlt12.ARTIFACT_ID
        or gates.get("PROTO14_successor_runtime_implemented") is not True
        or gates.get("PROTO14_fresh_GR0_dynamic_calibration_authorized") is not True
        or gates.get("classical_spherical_diagnostic_authorized") is not False
        or gates.get("SGBL_execution_authorized") is not False
        or gates.get("FGCQR_holdout_execution_authorized") is not False
        or frozen.get("run_plan_sha256") != hlt12._sha(plan_path)
        or frozen.get("primary_point_counts") != [2049, 4097, 8193]
        or frozen.get("comparator_point_counts") != [4097, 8193, 16385]
        or frozen.get("restart_coordinate_time") != RESTART_TIME
        or frozen.get("TDG6_every_accepted_macro_step_required") is not True
        or frozen.get("sampled_64_history_is_nonveto_diagnostic_only") is not True
        or payload.get("runtime_environment_observed") != runtime_environment
    ):
        raise ValueError("HLT12 authorization is absent, incomplete, or over-broad")
    for ledger_name in ("source_config_sha256", "predecessor_sha256", "implementation_sha256"):
        ledger = authorization.get(ledger_name)
        if not isinstance(ledger, dict) or not ledger:
            raise ValueError(f"HLT12 {ledger_name} is absent")
        for relative, expected in ledger.items():
            owned = REPOSITORY / relative
            if not owned.is_file() or hlt12._sha(owned) != expected:
                raise ValueError(f"HLT12-authorized file drifted: {relative}")
    owner = REPOSITORY / authorization.get("derivation_document", "")
    if not owner.is_file() or hlt12._sha(owner) != authorization.get("derivation_document_sha256"):
        raise ValueError("HLT12 owner document drifted")
    if require_fresh_namespace:
        reproduced = hlt12.record(hlt12.DEFAULT_CONFIG)
        if inherited._serial(reproduced) != authorization:
            raise ValueError("HLT12 authorization does not reproduce before launch")
    return plan, authorization

def _restore_frozen_members(
    authorization: Mapping[str, Any],
) -> dict[str, Proto14RunMember]:
    """Restore the immutable six payloads and add fresh zero TDG6 ledgers."""

    # The PROTO13 restoration helper reads only its immutable PROTO12/RSP2 and
    # PROTO13-FRZ1 lineage.  PROTO14 deliberately changes no input array.
    base = proto13_runner._restore_frozen_members(authorization)
    members = {
        key: Proto14RunMember.from_proto7(member)
        for key, member in base.items()
    }
    if tuple(sorted(members)) != tuple(sorted(_MEMBER_KEYS)):
        raise RuntimeError("PROTO14 member shells differ from the frozen ladders")
    for member in members.values():
        ledger = member.temporal_ledger
        if (
            member.time != RESTART_TIME
            or ledger.last_accepted_time != RESTART_TIME
            or ledger.accepted_macro_step_count != 0
            or ledger.cumulative_temporal_retry_count != 0
            or ledger.serialized_temporal_rejections
            or any(ledger.accumulated_debit_vector)
        ):
            raise ValueError(f"PROTO14 restart TDG6 ledger differs: {member.key}")
    return members


def _method_members(
    members: Mapping[str, Proto14RunMember], method: str
) -> tuple[Proto14RunMember, ...]:
    expected = (2049, 4097, 8193) if method == "RK4" else (4097, 8193, 16385)
    records = tuple(sorted(
        (item for item in members.values() if item.method_label == method),
        key=lambda item: item.point_count,
    ))
    if tuple(item.point_count for item in records) != expected:
        raise ValueError(f"PROTO14 {method} member ownership differs")
    return records


def _temporal_ledger_summary(member: Proto14RunMember) -> dict[str, Any]:
    ledger = member.temporal_ledger
    return {
        "last_accepted_time": ledger.last_accepted_time,
        "accepted_macro_step_count": ledger.accepted_macro_step_count,
        "cumulative_temporal_retry_count": ledger.cumulative_temporal_retry_count,
        "current_macro_step_temporal_retry_count": ledger.current_macro_step_temporal_retry_count,
        "last_accepted_macro_step_temporal_retry_count": ledger.last_accepted_macro_step_temporal_retry_count,
        "accumulated_debit_vector": list(ledger.accumulated_debit_vector),
        "accumulated_debit_sha256": array_content_sha256(
            np.asarray(ledger.accumulated_debit_vector, dtype=np.float64)
        ),
        "serialized_temporal_rejection_count": len(ledger.serialized_temporal_rejections),
        "physical_classification": False,
    }


def _tdg6_ledger_consistent(
    member: Proto14RunMember, coordinate_time: float
) -> bool:
    """Check the complete post-restart TDG6 ledger at one common boundary.

    This does not turn TDG6's debit into a physical error estimate.  It only
    proves that the currently accepted member boundary is backed by an
    admitted fine-path ledger, with no unresolved temporal retry.
    """

    ledger = member.temporal_ledger
    debit = np.asarray(ledger.accumulated_debit_vector, dtype=np.float64)
    if debit.shape != (18,) or not np.all(np.isfinite(debit)) or np.any(debit < 0.0):
        return False
    if (
        np.float64(member.time).tobytes() != np.float64(coordinate_time).tobytes()
        or np.float64(ledger.last_accepted_time).tobytes()
        != np.float64(member.time).tobytes()
        or ledger.current_macro_step_temporal_retry_count != 0
        or ledger.accepted_macro_step_count < 1
        or ledger.cumulative_temporal_retry_count
        != len(ledger.serialized_temporal_rejections)
    ):
        return False
    try:
        for serialized in ledger.serialized_temporal_rejections:
            event = json.loads(serialized)
            if (
                event.get("event_type") != "rejected_TDG6_temporal_admission"
                or inherited._canonical(event).decode("utf-8").strip() != serialized
            ):
                return False
    except (TypeError, ValueError, json.JSONDecodeError):
        return False
    return True


def _event_assessment(
    plan: Mapping[str, Any], members: Mapping[str, Proto14RunMember], coordinate_time: float
) -> dict[str, Any]:
    """Reuse PROTO13 spatial gates and demote only its sampled temporal veto."""

    base = proto13_runner._event_assessment(plan, members, coordinate_time)
    sampled = dict(base["temporal_spectral_admission"])
    sampled["used_as_veto"] = False
    sampled["classification_scope"] = "diagnostic_only_not_temporal_admission"
    base["sampled_temporal_history_diagnostic"] = sampled
    base["temporal_spectral_admission"] = {
        "replaced_by": "FGC-1-TDG6-IMP2",
        "sampled_history_rule_used_as_veto": False,
        "all_post_restart_accepted_macro_steps_TDG6_admitted": True,
    }
    tdg6_consistent = all(
        _tdg6_ledger_consistent(item, coordinate_time)
        for item in members.values()
    )
    base["TDG6_ledger_consistent_at_common_event"] = tdg6_consistent
    base["qualified_trapped_common_event"] = (
        base["trapped_assessment"].trapped_sign_passed
        and base["primary_common_event"].admission_passed
        and base["comparator_common_event"].admission_passed
        and tdg6_consistent
    )
    base["member_temporal_ledgers"] = {
        key: _temporal_ledger_summary(item) for key, item in members.items()
    }
    base["accumulated_temporal_debit_is_Raychaudhuri_error_bound"] = False
    return base


def _manifest(
    plan_path: Path, authorization_path: Path, authorization: Mapping[str, Any]
) -> dict[str, Any]:
    frozen = authorization["artifact_payload"]["frozen_run_plan"]
    return {
        "schema_version": 1,
        "runner_id": RUNNER_ID,
        "campaign_id": frozen["campaign_id"],
        "authorization_checkpoint_commit": inherited._git_head(),
        "plan_path": _relative(plan_path),
        "plan_sha256": hlt12._sha(plan_path),
        "authorization_path": _relative(authorization_path),
        "authorization_sha256": hlt12._sha(authorization_path),
        "implementation_sha256": authorization["implementation_sha256"],
        "numerical_runtime_contract": dict(frozen["numerical_runtime_contract"]),
        "runtime_environment": hlt12.validate_numerical_runtime(_validate_run_plan(plan_path)["numerical_runtime"]),
        "restart_coordinate_time": RESTART_TIME,
        "tdg6_all_step_temporal_admission": True,
        "sampled_temporal_history_diagnostic_only": True,
        "created_unix_time": wall_time.time(),
        "outcome_fields_present_at_creation": False,
    }


def _member_metadata(member: Proto14RunMember) -> dict[str, Any]:
    value = proto6_runner._member_metadata(member)
    value["state_sha256"] = array_content_sha256(member.state.u, member.state.p, member.state.q)
    extension, _ = proto14_checkpoint_extension(member)
    value["tdg6_temporal_ledger"] = extension
    return value


def _checkpoint_metadata(
    *, manifest: Mapping[str, Any], event_index: int, consecutive: int,
    event_log: Path, members: Mapping[str, Proto14RunMember], terminal: bool,
    terminal_result: Mapping[str, Any] | None,
) -> dict[str, Any]:
    serialized_terminal_result = (
        None if terminal_result is None else inherited._serial(terminal_result)
    )
    return {
        "schema_version": 2,
        "campaign_id": manifest["campaign_id"],
        "manifest_sha256": sha256(inherited._canonical(manifest)).hexdigest(),
        "amplitude": "3",
        "amplitude_index": 0,
        "completed_common_event_index": event_index,
        "consecutive_qualified_events": consecutive,
        "completed_amplitude_records": [],
        "event_log": inherited._event_log_identity(inherited._event_log_bytes(event_log)),
        "terminal": terminal,
        "terminal_result": serialized_terminal_result,
        "terminal_result_sha256": (
            None
            if serialized_terminal_result is None
            else sha256(inherited._canonical(serialized_terminal_result)).hexdigest()
        ),
        "members": {key: _member_metadata(item) for key, item in members.items()},
        "TDG6_checkpoint_extension_present": True,
        "sampled_temporal_history_veto_enabled": False,
    }


def _checkpoint_arrays(
    *, metadata: Mapping[str, Any], members: Mapping[str, Proto14RunMember], event_log: Path
) -> dict[str, np.ndarray]:
    arrays = inherited._checkpoint_arrays(
        metadata=metadata, members=members, event_log_payload=inherited._event_log_bytes(event_log)
    )
    for key, member in members.items():
        _, debit = proto14_checkpoint_extension(member)
        arrays[f"{key.replace('-', '_')}_tdg6_accumulated_debit"] = debit
    return arrays


def _write_checkpoint(
    path: Path, *, metadata: Mapping[str, Any], members: Mapping[str, Proto14RunMember], event_log: Path
) -> None:
    if metadata.get("event_log") != inherited._event_log_identity(inherited._event_log_bytes(event_log)):
        raise ValueError("PROTO14 checkpoint metadata does not match event log")
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}.npz")
    if temporary.exists():
        raise FileExistsError(temporary)
    arrays = _checkpoint_arrays(metadata=metadata, members=members, event_log=event_log)
    with temporary.open("xb") as handle:
        np.savez_compressed(handle, **arrays)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def _restore_campaign_checkpoint(
    path: Path, *, manifest: Mapping[str, Any]
) -> tuple[dict[str, Any], dict[str, Proto14RunMember], bytes]:
    base = proto13_runner._build_frozen_member_shells()
    members = {key: Proto14RunMember.from_proto7(member) for key, member in base.items()}
    with np.load(path, allow_pickle=False) as archive:
        metadata = json.loads(bytes(archive["metadata_utf8"]).decode("utf-8"))
        event_log = bytes(archive["event_log_utf8"])
        if (
            metadata.get("schema_version") != 2
            or metadata.get("campaign_id") != manifest["campaign_id"]
            or metadata.get("manifest_sha256") != sha256(inherited._canonical(manifest)).hexdigest()
            or metadata.get("event_log") != inherited._event_log_identity(event_log)
            or metadata.get("TDG6_checkpoint_extension_present") is not True
            or metadata.get("sampled_temporal_history_veto_enabled") is not False
            or set(metadata.get("members", {})) != set(members)
        ):
            raise ValueError("PROTO14 checkpoint identity differs")
        _checkpoint_terminal_result(
            metadata=metadata,
            manifest=manifest,
            event_log_payload=event_log,
        )
        for key, member in members.items():
            prefix = key.replace("-", "_")
            stored = metadata["members"][key]
            if stored["input_hash"] != member.input_hash:
                raise ValueError(f"PROTO14 checkpoint input hash differs: {key}")
            member.state = EvolutionState(archive[f"{prefix}_u"].copy(), archive[f"{prefix}_p"].copy(), archive[f"{prefix}_q"].copy())
            if array_content_sha256(member.state.u, member.state.p, member.state.q) != stored["state_sha256"]:
                raise ValueError(f"PROTO14 checkpoint state hash differs: {key}")
            member.time = float(stored["time"])
            member.step_index = int(stored["step_index"])
            member.transaction_serial = int(stored["transaction_serial"])
            member.CFL_retry_count = int(stored["CFL_retry_count"])
            member.source_retry_count = int(stored["source_retry_count"])
            member.transaction.state = GR0RuntimeMonitorState(**stored["runtime_monitor_state"])
            member.transaction.causal_state = CausalBudgetState(**stored["causal_state"])
            member.tracers.positions = archive[f"{prefix}_tracer_positions"].copy()
            member.tracers.proper_times = archive[f"{prefix}_tracer_proper_times"].copy()
            member.tracers.event_proper_times = [row.copy() for row in archive[f"{prefix}_event_proper_times"]]
            member.tracers.event_fields = [row.copy() for row in archive[f"{prefix}_event_fields"]]
            restore_proto14_checkpoint_extension(
                member, stored["tdg6_temporal_ledger"], archive[f"{prefix}_tdg6_accumulated_debit"]
            )
    return metadata, members, event_log


def _validate_resume_manifest(
    manifest: Mapping[str, Any], plan_path: Path, authorization_path: Path,
    authorization: Mapping[str, Any],
) -> None:
    frozen = authorization["artifact_payload"]["frozen_run_plan"]
    if (
        manifest.get("runner_id") != RUNNER_ID
        or manifest.get("campaign_id") != frozen["campaign_id"]
        or manifest.get("plan_sha256") != hlt12._sha(plan_path)
        or manifest.get("authorization_sha256") != hlt12._sha(authorization_path)
        or manifest.get("implementation_sha256") != authorization["implementation_sha256"]
        or manifest.get("authorization_checkpoint_commit") != inherited._git_head()
        or manifest.get("numerical_runtime_contract") != frozen["numerical_runtime_contract"]
        or manifest.get("tdg6_all_step_temporal_admission") is not True
        or manifest.get("sampled_temporal_history_diagnostic_only") is not True
    ):
        raise ValueError("PROTO14 resume manifest differs from authorization")


def _load_canonical_campaign_result(path: Path) -> dict[str, Any]:
    """Load a terminal result only when its complete bytes are canonical."""

    try:
        payload = path.read_bytes()
        value = json.loads(payload.decode("utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError("PROTO14 terminal result is unreadable") from error
    if not isinstance(value, dict) or payload != inherited._canonical(value):
        raise ValueError("PROTO14 terminal result is not canonical")
    return value


def _checkpoint_terminal_result(
    *,
    metadata: Mapping[str, Any],
    manifest: Mapping[str, Any],
    event_log_payload: bytes,
) -> dict[str, Any] | None:
    """Validate the checkpoint-owned half of the terminal transaction."""

    terminal = metadata.get("terminal")
    result = metadata.get("terminal_result")
    result_sha256 = metadata.get("terminal_result_sha256")
    if terminal is False:
        if result is not None or result_sha256 is not None:
            raise ValueError("PROTO14 nonterminal checkpoint carries terminal evidence")
        return None
    if terminal is not True or not isinstance(result, dict):
        raise ValueError("PROTO14 checkpoint terminal contract differs")
    canonical_result = inherited._serial(result)
    canonical_sha256 = sha256(inherited._canonical(canonical_result)).hexdigest()
    records = canonical_result.get("amplitude_records")
    members = metadata.get("members")
    event_identity = inherited._event_log_identity(event_log_payload)
    allowed_classifications = {
        "GR0_calibration_completed_selected_amplitude",
        "calibration_failed_no_eligible_GR0_case",
        "invalid_implementation_or_nonconverged_run",
    }
    if (
        result_sha256 != canonical_sha256
        or canonical_result.get("schema_version") != 2
        or canonical_result.get("runner_id") != RUNNER_ID
        or canonical_result.get("campaign_id") != manifest.get("campaign_id")
        or canonical_result.get("classification") not in allowed_classifications
        or canonical_result.get("event_log_sha256") != event_identity["sha256"]
        or not isinstance(records, list)
        or len(records) != 1
        or not isinstance(records[0], dict)
        or records[0].get("last_completed_common_event_index")
        != metadata.get("completed_common_event_index")
        or records[0].get("consecutive_qualified_trapped_common_events")
        != metadata.get("consecutive_qualified_events")
        or not isinstance(members, dict)
        or records[0].get("member_final_state_hashes")
        != {key: value.get("state_sha256") for key, value in members.items()}
        or canonical_result.get("FGCQR_outcome_read") is not False
        or canonical_result.get("SGBL_outcome_read") is not False
        or canonical_result.get("holdout_execution_authorized") is not False
        or canonical_result.get("retained_EFT_evolution_authorized") is not False
        or canonical_result.get("mechanism_question_answered") is not False
    ):
        raise ValueError("PROTO14 terminal checkpoint result differs")
    return canonical_result


def _publish_or_verify_terminal_result(
    path: Path, result: Mapping[str, Any]
) -> dict[str, Any]:
    """Atomically publish once, or prove an already-published exact match."""

    canonical_result = inherited._serial(result)
    payload = inherited._canonical(canonical_result)
    if path.exists():
        observed = _load_canonical_campaign_result(path)
        if observed != canonical_result:
            raise ValueError("PROTO14 published result differs from terminal checkpoint")
        return observed

    temporary = path.with_name(
        f".{path.name}.pending-{os.getpid()}-{wall_time.time_ns()}"
    )
    try:
        inherited._write_exclusive(temporary, payload)
        try:
            os.link(temporary, path)
        except FileExistsError:
            observed = _load_canonical_campaign_result(path)
            if observed != canonical_result:
                raise ValueError(
                    "PROTO14 concurrently published result differs from terminal checkpoint"
                )
            return observed
        directory_descriptor = os.open(path.parent, os.O_RDONLY)
        try:
            os.fsync(directory_descriptor)
        finally:
            os.close(directory_descriptor)
    finally:
        temporary.unlink(missing_ok=True)
    return canonical_result


def _reconcile_terminal_publication(
    *,
    metadata: Mapping[str, Any],
    manifest: Mapping[str, Any],
    event_log_payload: bytes,
    event_log_relative: str,
    result_path: Path,
) -> dict[str, Any] | None:
    """Finish or verify an interrupted terminal publication idempotently."""

    result = _checkpoint_terminal_result(
        metadata=metadata,
        manifest=manifest,
        event_log_payload=event_log_payload,
    )
    if result is None:
        if result_path.exists():
            raise ValueError("PROTO14 result exists beside a nonterminal checkpoint")
        return None
    if result.get("event_log") != event_log_relative:
        raise ValueError("PROTO14 terminal result names a different event log")
    return _publish_or_verify_terminal_result(result_path, result)


def _durable_attempt_sink(
    event_log: Path, *, attempt_id: str, member: Proto14RunMember, kind: str,
):
    def sink(evidence: Mapping[str, Any]) -> None:
        inherited._append_fsync(event_log, inherited._canonical_line({
            **dict(evidence),
            "common_event_attempt_id": attempt_id,
            "member": member.key,
            "method": member.method_label,
            "point_count": member.point_count,
            "retry_kind": kind,
            "classification_is_numerical_not_physical": True,
        }))
    return sink


def run_campaign(*, plan_path: Path, authorization_path: Path, resume: bool) -> dict[str, Any]:
    inherited._require_clean_tracked_worktree()
    plan, authorization = _validate_authorization(
        plan_path, authorization_path, require_fresh_namespace=not resume
    )
    provenance = plan["provenance"]
    output_root = REPOSITORY / provenance["output_root"]
    manifest_path = REPOSITORY / provenance["campaign_manifest"]
    event_log = output_root / provenance["append_only_event_log_name"]
    checkpoint_path = output_root / provenance["checkpoint_name"]
    result_path = output_root / provenance["result_name"]
    resume_recovery = None
    if resume:
        if not output_root.is_dir() or not manifest_path.is_file():
            raise ValueError("PROTO14 resume requires an existing campaign")
        manifest = inherited._load_manifest(manifest_path)
        _validate_resume_manifest(manifest, plan_path, authorization_path, authorization)
        metadata, members, saved_log = _restore_campaign_checkpoint(checkpoint_path, manifest=manifest)
        resume_recovery = inherited._recover_event_log_from_checkpoint(event_log, saved_log)
        terminal_result = _reconcile_terminal_publication(
            metadata=metadata,
            manifest=manifest,
            event_log_payload=inherited._event_log_bytes(event_log),
            event_log_relative=_relative(event_log),
            result_path=result_path,
        )
        if terminal_result is not None:
            return terminal_result
        event_index = int(metadata["completed_common_event_index"])
        consecutive = int(metadata["consecutive_qualified_events"])
    else:
        if output_root.exists():
            raise ValueError("fresh PROTO14 output root already exists")
        members = _restore_frozen_members(authorization)
        initial = _event_assessment(plan, members, RESTART_TIME)
        if not initial["primary_common_event"].admission_passed or not initial["comparator_common_event"].admission_passed:
            raise ValueError("PROTO14 restart spatial or constraint admission changed after HLT12")
        output_root.mkdir(parents=True, exist_ok=False)
        manifest = _manifest(plan_path, authorization_path, authorization)
        inherited._write_exclusive(manifest_path, inherited._canonical(manifest))
        event_index, consecutive = 23, 0
        inherited._append_fsync(event_log, inherited._canonical_line({
            "event_type": "restart_common_event", "event_index": event_index,
            "assessment": initial, "trajectory_advanced_by_PROTO14": False,
            "counts_toward_new_consecutive_trapped_events": False,
        }))
        _write_checkpoint(checkpoint_path, metadata=_checkpoint_metadata(
            manifest=manifest, event_index=event_index, consecutive=consecutive,
            event_log=event_log, members=members, terminal=False, terminal_result=None,
        ), members=members, event_log=event_log)

    output_interval, checkpoint_every, final_event_index = 1.0 / 16.0, 4, 512
    required_consecutive, retry_factor, minimum_step = 8, 0.5, 1.0 / 1073741824.0
    candidate_stop: dict[str, Any] | None = None
    selected = False
    started = wall_time.monotonic()
    while event_index < final_event_index:
        target_index = event_index + 1
        target_time = target_index * output_interval
        attempt_id = f"proto14-common-{target_index:04d}"
        snapshots = {key: item.snapshot() for key, item in members.items()}
        active: Proto14RunMember | None = None
        try:
            for active in members.values():
                active.advance_to(
                    target_time, retry_factor=retry_factor, maximum_CFL_retries=32,
                    maximum_source_retries=32, minimum_step_size=minimum_step,
                    rejected_trial_sink=_durable_attempt_sink(event_log, attempt_id=attempt_id, member=active, kind="PROTO7_source"),
                    temporal_rejection_sink=_durable_attempt_sink(event_log, attempt_id=attempt_id, member=active, kind="TDG6_temporal"),
                )
            for member in members.values():
                member.tracers.append_common_event(
                    member.state, member.initial.grid.coordinates
                )
            assessment = _event_assessment(plan, members, target_time)
        except Proto7TerminalStop as stop:
            candidate_stop = {"classification": "scientific_universal_runtime_stop", "reason": stop.reason,
                "failures": list(stop.failures), "complete_failure_evidence": asdict(stop.evidence),
                "member": None if active is None else active.key, "target_common_event_index": target_index,
                "last_fully_completed_common_event_index": event_index, "cross_member_state_rolled_back": True}
        except Proto14SourceRetryExhausted as stop:
            candidate_stop = {"classification": "scientific_source_retry_exhausted", "reason": stop.reason,
                "complete_failure_evidence": stop.evidence, "member": None if active is None else active.key,
                "target_common_event_index": target_index, "last_fully_completed_common_event_index": event_index,
                "cross_member_state_rolled_back": True}
        except TDG6TemporalRetryExhausted as stop:
            candidate_stop = {"classification": "invalid_implementation_or_nonconverged_run",
                "reason": "TDG6_temporal_retry_exhausted", "TDG6_reason": stop.reason,
                "complete_failure_evidence": asdict(stop.evidence),
                "member": None if active is None else active.key,
                "target_common_event_index": target_index,
                "last_fully_completed_common_event_index": event_index,
                "cross_member_state_rolled_back": True}
        except BaseException as error:
            if isinstance(error, KeyboardInterrupt):
                raise
            candidate_stop = {"classification": "invalid_implementation_or_nonconverged_run",
                "reason": inherited._source_failure_reason(error), "error_type": type(error).__name__,
                "error_message": str(error), "member": None if active is None else active.key,
                "target_common_event_index": target_index, "last_fully_completed_common_event_index": event_index,
                "cross_member_state_rolled_back": True}
        if candidate_stop is not None:
            inherited._append_fsync(event_log, inherited._canonical_line({
                "event_type": "common_event_attempt_aborted", "common_event_attempt_id": attempt_id,
                "target_common_event_index": target_index, "last_fully_completed_common_event_index": event_index,
                "active_member": None if active is None else active.key, "stop": candidate_stop,
                "atomic_checkpoint_owns_accepted_state": True,
                "log_tail_is_recovery_orphan_evidence_not_physical_classification": True,
            }))
            restore_common_event_snapshots(members, snapshots)
            break

        event_index = target_index
        consecutive = consecutive + 1 if assessment["qualified_trapped_common_event"] else 0
        inherited._append_fsync(event_log, inherited._canonical_line({
            "event_type": "common_event", "common_event_attempt_id": attempt_id, "amplitude": "3",
            "event_index": event_index, "assessment": assessment,
            "consecutive_qualified_trapped_common_events": consecutive,
        }))
        if event_index % checkpoint_every == 0 or consecutive > 0:
            _write_checkpoint(checkpoint_path, metadata=_checkpoint_metadata(
                manifest=manifest, event_index=event_index, consecutive=consecutive,
                event_log=event_log, members=members, terminal=False, terminal_result=None,
            ), members=members, event_log=event_log)
        if event_index % 8 == 0 or consecutive > 0:
            trapped = assessment["trapped_assessment"]
            print(f"PROTO14 event={event_index}/{final_event_index} t={target_time:.6g} "
                  f"trapped_score={trapped.primary_trapped_scores[-1]:.6g}/"
                  f"{trapped.comparator_trapped_scores[-1]:.6g} qualified_run={consecutive}/{required_consecutive}", flush=True)
        if consecutive >= required_consecutive:
            selected = True
            break
        common_pass = assessment["primary_common_event"].admission_passed and assessment["comparator_common_event"].admission_passed
        if not common_pass:
            candidate_stop = {"classification": "common_event_constraint_or_spatial_stop",
                "reason": "common_event_constraint_or_spatial_spectral_admission",
                "event_index": event_index, "coordinate_time": target_time}
            break

    classification = ("GR0_calibration_completed_selected_amplitude" if selected else
        "invalid_implementation_or_nonconverged_run" if candidate_stop and candidate_stop["classification"] == "invalid_implementation_or_nonconverged_run"
        else "calibration_failed_no_eligible_GR0_case")
    record = {"amplitude": "3", "last_completed_common_event_index": event_index,
        "last_completed_coordinate_time": event_index * output_interval,
        "consecutive_qualified_trapped_common_events": consecutive, "eligible": selected,
        "stop": candidate_stop,
        "member_source_retry_counts": {key: item.source_retry_count for key, item in members.items()},
        "member_CFL_retry_counts": {key: item.CFL_retry_count for key, item in members.items()},
        "member_temporal_ledgers": {key: _temporal_ledger_summary(item) for key, item in members.items()},
        "member_final_state_hashes": {key: array_content_sha256(item.state.u, item.state.p, item.state.q) for key, item in members.items()}}
    result = {"schema_version": 2, "runner_id": RUNNER_ID, "campaign_id": manifest["campaign_id"],
        "classification": classification, "selected_amplitude": "3" if selected else None,
        "amplitude_records": [record], "event_log": _relative(event_log), "event_log_sha256": hlt12._sha(event_log),
        "elapsed_wall_seconds": wall_time.monotonic() - started, "resume_recovery": resume_recovery,
        "restart_coordinate_time": RESTART_TIME, "numerical_runtime_contract": manifest["numerical_runtime_contract"],
        "runtime_environment": manifest["runtime_environment"], "method_owned_ladders_enabled": True,
        "common_physical_node_comparison_enabled": True, "TDG6_all_step_temporal_admission_enabled": True,
        "sampled_temporal_history_veto_enabled": False, "FGCQR_outcome_read": False, "SGBL_outcome_read": False,
        "holdout_execution_authorized": False, "retained_EFT_evolution_authorized": False,
        "mechanism_question_answered": False}
    terminal_metadata = _checkpoint_metadata(
        manifest=manifest, event_index=event_index, consecutive=consecutive, event_log=event_log,
        members=members, terminal=True, terminal_result=result,
    )
    _write_checkpoint(
        checkpoint_path,
        metadata=terminal_metadata,
        members=members,
        event_log=event_log,
    )
    published = _reconcile_terminal_publication(
        metadata=terminal_metadata,
        manifest=manifest,
        event_log_payload=inherited._event_log_bytes(event_log),
        event_log_relative=_relative(event_log),
        result_path=result_path,
    )
    if published is None:
        raise RuntimeError("PROTO14 terminal checkpoint did not own a result")
    return published


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, default=DEFAULT_PLAN)
    parser.add_argument("--authorization", type=Path, default=DEFAULT_AUTHORIZATION)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--check-authorization-only", action="store_true")
    args = parser.parse_args()
    if args.check_authorization_only:
        inherited._require_clean_tracked_worktree()
        plan, authorization = _validate_authorization(args.plan.resolve(), args.authorization.resolve(), require_fresh_namespace=True)
        print(json.dumps({"authorized": True, "campaign_id": authorization["artifact_payload"]["frozen_run_plan"]["campaign_id"],
            "amplitude": plan["scope"]["amplitude"], "TDG6_all_step_temporal_admission": True,
            "output_created": False}, sort_keys=True))
        return
    print(json.dumps(run_campaign(plan_path=args.plan.resolve(), authorization_path=args.authorization.resolve(), resume=args.resume), sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
