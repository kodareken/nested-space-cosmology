"""Read-only reconstruction of the PROTO14 immediate pre-shadow terminal.

This module owns no runner, output namespace, checkpoint writer, or resumption
operation.  It proves only the recorded terminal transaction: the checkpoint
owns the accepted event-23 boundary, embeds the complete two-line event log and
the external terminal result, and records the deterministic TDG6 subdivision
validator fault before any new macro step was accepted.
"""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .numerical_engine import array_content_sha256


MEMBERS = (
    "RK4-2049", "RK4-4097", "RK4-8193",
    "SSPRK3-4097", "SSPRK3-8193", "SSPRK3-16385",
)
ZERO_DEBIT_SHA256 = "a2f69ec892ac5ccef01ceb3ffe952383a3c96eae9c538e1263dee545c12e7af9"
MEMBER_SUFFIXES = (
    "u", "p", "q", "tracer_positions", "tracer_proper_times",
    "event_proper_times", "event_fields",
)


def canonical(value: Mapping[str, Any]) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True, ensure_ascii=True,
                       allow_nan=False) + "\n").encode("utf-8")


def sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def load_canonical_json(path: Path) -> dict[str, Any]:
    payload = path.read_bytes()
    value = json.loads(payload)
    if not isinstance(value, dict) or payload != canonical(value):
        raise ValueError(f"{path} is not canonical JSON")
    return value


def load_canonical_jsonl(path: Path) -> tuple[list[dict[str, Any]], bytes]:
    payload = path.read_bytes()
    if not payload.endswith(b"\n"):
        raise ValueError("event log does not end in a newline")
    records = [json.loads(line) for line in payload.splitlines()]
    if any(not isinstance(record, dict) for record in records):
        raise ValueError("event log record is not an object")
    rebuilt = b"".join(
        (json.dumps(record, sort_keys=True, separators=(",", ":"),
                    ensure_ascii=True, allow_nan=False) + "\n").encode("utf-8")
        for record in records
    )
    if rebuilt != payload:
        raise ValueError("event log is not canonical JSONL")
    return records, payload


def _checkpoint_payload(path: Path) -> tuple[dict[str, Any], bytes, dict[str, Any]]:
    with np.load(path, allow_pickle=False) as archive:
        required = {"metadata_utf8", "event_log_utf8"}
        for member in MEMBERS:
            prefix = member.replace("-", "_")
            required.update(f"{prefix}_{suffix}" for suffix in MEMBER_SUFFIXES)
            required.add(f"{prefix}_tdg6_accumulated_debit")
        if set(archive.files) != required:
            raise ValueError("checkpoint array schema differs")
        metadata = json.loads(bytes(archive["metadata_utf8"]).decode("utf-8"))
        event_bytes = bytes(archive["event_log_utf8"])
        member_arrays: dict[str, Any] = {}
        for member in MEMBERS:
            prefix = member.replace("-", "_")
            names = tuple(f"{prefix}_{suffix}" for suffix in MEMBER_SUFFIXES)
            arrays = tuple(archive[name].copy() for name in names)
            debit = archive[f"{prefix}_tdg6_accumulated_debit"].copy()
            member_arrays[member] = {
                "state_sha256": array_content_sha256(*arrays[:3]),
                "restart_payload_sha256": array_content_sha256(*arrays),
                "state_shape": list(arrays[0].shape),
                "tracer_count": int(arrays[3].shape[0]),
                "event_sample_count": int(arrays[5].shape[0]),
                "debit": debit,
            }
    if not isinstance(metadata, dict):
        raise ValueError("checkpoint metadata is not an object")
    return metadata, event_bytes, member_arrays


def reconstruct_terminal(
    *, manifest_path: Path, event_log_path: Path, checkpoint_path: Path,
    result_path: Path, expected: Mapping[str, Any],
    expected_members: Mapping[str, Mapping[str, Any]],
) -> dict[str, Any]:
    """Verify the raw four-file identity and zero-progress terminal boundary."""

    manifest = load_canonical_json(manifest_path)
    result = load_canonical_json(result_path)
    events, event_bytes = load_canonical_jsonl(event_log_path)
    metadata, embedded_events, arrays = _checkpoint_payload(checkpoint_path)
    if len(events) != 2 or len(event_bytes) != int(expected["event_log_byte_count"]):
        raise ValueError("PROTO14 event-log identity differs")
    if embedded_events != event_bytes or sha256(event_bytes).hexdigest() != expected["event_log_sha256"]:
        raise ValueError("checkpoint does not embed the exact external event log")
    lineage = expected["authorization_lineage"]
    if (
        manifest.get("campaign_id") != expected["campaign_id"]
        or result.get("campaign_id") != expected["campaign_id"]
        or manifest.get("authorization_checkpoint_commit") != expected["authorization_commit"]
        or manifest.get("authorization_sha256") != lineage["runtime_authorization_sha256"]
        or manifest.get("plan_sha256") != lineage["run_plan_sha256"]
        or manifest.get("runner_id") != lineage["runner_id"]
        or manifest.get("implementation_sha256", {}).get(lineage["runner"])
        != lineage["runner_sha256"]
        or manifest.get("implementation_sha256", {}).get(lineage["runtime_module"])
        != lineage["runtime_module_sha256"]
        or manifest.get("implementation_sha256", {}).get(lineage["tdg6_runtime_module"])
        != lineage["tdg6_runtime_module_sha256"]
    ):
        raise ValueError("campaign identity differs")
    if metadata.get("terminal") is not True or metadata.get("terminal_result") != result:
        raise ValueError("terminal checkpoint does not own the external result")
    if metadata.get("terminal_result_sha256") != sha256(canonical(result)).hexdigest():
        raise ValueError("terminal checkpoint result hash differs")
    if (
        metadata.get("completed_common_event_index") != 23
        or metadata.get("campaign_id") != expected["campaign_id"]
        or metadata.get("manifest_sha256") != expected["manifest_sha256"]
        or metadata.get("event_log") != {
            "sha256": expected["event_log_sha256"],
            "line_count": expected["event_log_line_count"],
            "byte_count": expected["event_log_byte_count"],
        }
        or metadata.get("TDG6_checkpoint_extension_present") is not True
    ):
        raise ValueError("checkpoint does not retain event-23 boundary")
    if events[0].get("event_type") != "restart_common_event" or events[0].get("event_index") != 23:
        raise ValueError("restart event identity differs")
    if events[0].get("trajectory_advanced_by_PROTO14") is not False:
        raise ValueError("restart was misclassified as new trajectory progress")
    abort = events[1]
    stop = abort.get("stop", {})
    identity = expected["terminal_identity"]
    if (
        abort.get("event_type") != "common_event_attempt_aborted"
        or abort.get("common_event_attempt_id") != "proto14-common-0024"
        or abort.get("active_member") != identity["active_member"]
        or abort.get("target_common_event_index") != identity["target_common_event_index"]
        or abort.get("last_fully_completed_common_event_index") != identity["last_fully_completed_common_event_index"]
        or abort.get("atomic_checkpoint_owns_accepted_state") is not True
        or abort.get("log_tail_is_recovery_orphan_evidence_not_physical_classification") is not True
        or any(stop.get(key) != identity[key] for key in ("classification", "reason", "error_type", "error_message"))
        or stop.get("cross_member_state_rolled_back") is not True
    ):
        raise ValueError("aborted attempt identity differs")
    record = result.get("amplitude_records", [None])[0]
    if (
        not isinstance(record, dict)
        or result.get("classification") != identity["classification"]
        or result.get("runner_id") != lineage["runner_id"]
        or result.get("selected_amplitude") is not None
        or result.get("FGCQR_outcome_read") is not False
        or result.get("SGBL_outcome_read") is not False
        or result.get("holdout_execution_authorized") is not False
        or result.get("retained_EFT_evolution_authorized") is not False
        or result.get("mechanism_question_answered") is not False
        or record.get("last_completed_common_event_index") != 23
        or record.get("last_completed_coordinate_time") != 23 / 16
        or record.get("consecutive_qualified_trapped_common_events") != 0
        or record.get("eligible") is not False
        or record.get("stop") != stop
    ):
        raise ValueError("terminal result classification differs")
    member_summaries: dict[str, Any] = {}
    for member in MEMBERS:
        stored = metadata.get("members", {}).get(member, {})
        ledger = stored.get("tdg6_temporal_ledger", {})
        result_ledger = record.get("member_temporal_ledgers", {}).get(member, {})
        frozen = expected_members.get(member, {})
        if (
            not frozen
            or arrays[member]["state_sha256"] != frozen.get("state_sha256")
            or arrays[member]["restart_payload_sha256"] != frozen.get("restart_payload_sha256")
            or arrays[member]["state_shape"] != frozen.get("state_shape")
            or arrays[member]["tracer_count"] != frozen.get("tracer_count")
            or arrays[member]["event_sample_count"] != frozen.get("event_sample_count")
            or stored.get("state_sha256") != arrays[member]["state_sha256"]
            or stored.get("state_sha256") != record.get("member_final_state_hashes", {}).get(member)
            or stored.get("input_hash") != frozen.get("input_hash")
            or stored.get("time") != 23 / 16
            or stored.get("step_index") != frozen.get("step_index")
            or stored.get("transaction_serial") != frozen.get("transaction_serial")
            or stored.get("source_retry_count") != frozen.get("source_retry_count")
            or stored.get("CFL_retry_count") != frozen.get("CFL_retry_count")
            or stored.get("runtime_monitor_state", {}).get("accepted_stage_count")
            != frozen.get("accepted_stage_count")
            or stored.get("runtime_monitor_state", {}).get("last_accepted_time") != 23 / 16
            or stored.get("causal_state", {}).get("accepted_time") != 23 / 16
            or ledger.get("accepted_macro_step_count") != 0
            or ledger.get("cumulative_temporal_retry_count") != 0
            or ledger.get("current_macro_step_temporal_retry_count") != 0
            or ledger.get("last_accepted_macro_step_temporal_retry_count") != 0
            or ledger.get("serialized_temporal_rejections") != []
            or ledger.get("accumulated_debit_sha256") != ZERO_DEBIT_SHA256
            or ledger.get("last_accepted_time") != 23 / 16
            or ledger.get("physical_classification") is not False
            or arrays[member]["debit"].shape != (18,)
            or not np.array_equal(arrays[member]["debit"], np.zeros(18))
            or result_ledger.get("accepted_macro_step_count") != 0
            or result_ledger.get("cumulative_temporal_retry_count") != 0
            or result_ledger.get("serialized_temporal_rejection_count") != 0
            or result_ledger.get("accumulated_debit_sha256") != ZERO_DEBIT_SHA256
            or result_ledger.get("accumulated_debit_vector") != [0.0] * 18
        ):
            raise ValueError(f"{member} does not retain the zero-progress TDG6 boundary")
        member_summaries[member] = {
            "state_sha256": arrays[member]["state_sha256"],
            "restart_payload_sha256": arrays[member]["restart_payload_sha256"],
            "accepted_macro_step_count": 0,
            "temporal_retry_count": 0,
            "accumulated_debit_sha256": ZERO_DEBIT_SHA256,
        }
    return {
        "campaign_id": manifest["campaign_id"],
        "event_log_sha256": sha256(event_bytes).hexdigest(),
        "terminal_result_sha256": sha256(canonical(result)).hexdigest(),
        "completed_common_event_index": 23,
        "aborted_target_common_event_index": 24,
        "all_members_zero_progress": True,
        "checkpoint_embeds_event_log_and_result": True,
        "member_summaries": member_summaries,
    }


def source_pinned_static_replay(
    *, source_payload: bytes, expected: Mapping[str, Any]
) -> dict[str, Any]:
    """Verify the exact pre-shadow guard without executing a campaign or shadows.

    This is intentionally static.  The raw abort lacks a persisted proposed
    macro width, so it cannot support an exact arithmetic invocation without
    reconstructing and advancing mutable runner state.  The binder therefore
    proves the pinned source contains the recorded guard and leaves causal
    root-cause diagnosis to a prospective branch.
    """

    source = source_payload.decode("utf-8")
    if sha256(source_payload).hexdigest() != expected["runtime_source_sha256"]:
        raise ValueError("source-pinned replay source hash differs")
    for token in (expected["required_function"], expected["required_guard"], expected["required_exception"]):
        if token not in source:
            raise ValueError("source-pinned pre-shadow guard differs")
    return {
        "kind": expected["execution_scope"],
        "actual_arithmetic_invocation": False,
        "campaign_execution": False,
        "campaign_resume": False,
        "output_namespace_creation": False,
        "guard_present": True,
        "fault_cause_fully_derived": False,
    }
