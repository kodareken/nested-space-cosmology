"""Read-only, byte-exact legacy history replay for PREF26.

This layer receives already-read external/embedded journal bytes and decoded
checkpoint metadata.  It never normalizes a journal, opens a campaign path,
or writes a namespace.  A semantic match is deliberately insufficient: the
embedded journal must be the same byte sequence as its external counterpart.
"""

from __future__ import annotations

from hashlib import sha256
import json
import math
from typing import Any, Mapping, Sequence


class Proto18HistoricalReplayError(ValueError):
    """A legacy history is malformed, mismatched, or not the frozen grammar."""


PROTO12_EVENT_BYTES = 21_096_450
PROTO12_EVENT_SHA256 = "c2d35d07db2f6c983a449cdcf2db469360565cf392d78348dcafbd1d81234ba7"
PROTO12_EVENT_LINES = 49
RSP2_EVENT_BYTES = 15_060
RSP2_EVENT_SHA256 = "ca9168c36ae0592fc69d20bfb32f1b501100d162b50d500a2f6b9529cdac2937"
RSP2_EVENT_LINES = 25


def _reject_nonfinite(value: object, *, context: str) -> None:
    if isinstance(value, float) and not math.isfinite(value):
        raise Proto18HistoricalReplayError(f"{context}: nonfinite number")
    if isinstance(value, Mapping):
        for key, child in value.items():
            if not isinstance(key, str):
                raise Proto18HistoricalReplayError(f"{context}: non-text object key")
            _reject_nonfinite(child, context=context)
    elif isinstance(value, list):
        for child in value:
            _reject_nonfinite(child, context=context)


def _plain(value: object) -> object:
    """Copy immutable admitted metadata into ordinary JSON containers."""
    if isinstance(value, Mapping):
        return {key: _plain(child) for key, child in value.items()}
    if isinstance(value, tuple):
        return [_plain(child) for child in value]
    if isinstance(value, list):
        return [_plain(child) for child in value]
    return value


def _pairs(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    value: dict[str, Any] = {}
    for key, child in pairs:
        if key in value:
            raise Proto18HistoricalReplayError("duplicate JSON object key")
        value[key] = child
    return value


def _constant(token: str) -> object:
    raise Proto18HistoricalReplayError(f"nonfinite JSON constant: {token}")


def parse_jsonl_exact(
    payload: bytes,
    *,
    label: str,
    expected_bytes: int,
    expected_sha256: str,
    expected_lines: int,
) -> tuple[dict[str, Any], ...]:
    """Decode one frozen JSONL stream without changing or canonicalizing bytes."""
    if not isinstance(payload, bytes):
        raise Proto18HistoricalReplayError(f"{label}: payload is not bytes")
    if len(payload) != expected_bytes or sha256(payload).hexdigest() != expected_sha256:
        raise Proto18HistoricalReplayError(f"{label}: raw byte identity differs")
    if not payload.endswith(b"\n") or b"\r" in payload or payload.startswith(b"\xef\xbb\xbf"):
        raise Proto18HistoricalReplayError(f"{label}: JSONL framing differs")
    lines = payload[:-1].split(b"\n")
    if len(lines) != expected_lines or any(not line for line in lines):
        raise Proto18HistoricalReplayError(f"{label}: JSONL line count/framing differs")
    records: list[dict[str, Any]] = []
    for index, line in enumerate(lines):
        try:
            record = json.loads(
                line.decode("utf-8"), object_pairs_hook=_pairs, parse_constant=_constant,
            )
        except (UnicodeDecodeError, json.JSONDecodeError, Proto18HistoricalReplayError) as error:
            raise Proto18HistoricalReplayError(f"{label}: malformed JSONL record {index}") from error
        if not isinstance(record, dict):
            raise Proto18HistoricalReplayError(f"{label}: record {index} is not an object")
        _reject_nonfinite(record, context=f"{label}: record {index}")
        records.append(record)
    return tuple(records)


def _require_embedded_equal(external: bytes, embedded: bytes, *, label: str) -> None:
    if not isinstance(embedded, bytes) or embedded != external:
        raise Proto18HistoricalReplayError(f"{label}: embedded/external journal bytes differ")


def _event_identity(metadata: Mapping[str, Any], payload: bytes, lines: int, *, label: str) -> None:
    expected = {"byte_count": len(payload), "line_count": lines, "sha256": sha256(payload).hexdigest()}
    if metadata.get("event_log") != expected:
        raise Proto18HistoricalReplayError(f"{label}: checkpoint event-log identity differs")


def _common_time(record: Mapping[str, Any], index: int, *, label: str) -> None:
    assessment = record.get("assessment")
    if not isinstance(assessment, Mapping):
        raise Proto18HistoricalReplayError(f"{label}: common assessment is absent")
    for key in ("primary_common_event", "comparator_common_event"):
        event = assessment.get(key)
        if not isinstance(event, Mapping) or event.get("coordinate_time") != index / 16.0:
            raise Proto18HistoricalReplayError(f"{label}: common-event time differs")


def replay_proto12_history(
    external: bytes, embedded: bytes, metadata: Mapping[str, Any],
    selected_physical_state_sha256: Mapping[str, str],
) -> dict[str, Any]:
    """Verify PROTO12's 49-record two-amplitude legacy journal.

    The selected ``amplitude=3`` record at event 23 is the exact restart
    boundary inherited by PRO18; this function verifies it is cross-linked to
    the checkpoint's declared completed-amplitude and member metadata.
    """
    _require_embedded_equal(external, embedded, label="PROTO12")
    records = parse_jsonl_exact(
        external, label="PROTO12", expected_bytes=PROTO12_EVENT_BYTES,
        expected_sha256=PROTO12_EVENT_SHA256, expected_lines=PROTO12_EVENT_LINES,
    )
    if not isinstance(metadata, Mapping):
        raise Proto18HistoricalReplayError("PROTO12: checkpoint metadata is not an object")
    metadata = _plain(metadata)
    if not isinstance(metadata, Mapping):  # narrowed for static readers
        raise Proto18HistoricalReplayError("PROTO12: checkpoint metadata copy failed")
    _event_identity(metadata, external, PROTO12_EVENT_LINES, label="PROTO12")
    groups = (("5/2", 25), ("3", 24))
    offset = 0
    for amplitude, count in groups:
        for local, record in enumerate(records[offset: offset + count]):
            expected_keys = {"amplitude", "assessment", "event_index", "event_type"}
            if local:
                expected_keys.add("consecutive_qualified_trapped_common_events")
            if (
                set(record) != expected_keys or record.get("amplitude") != amplitude
                or record.get("event_type") != "common_event" or record.get("event_index") != local
            ):
                raise Proto18HistoricalReplayError("PROTO12: event grammar/order differs")
            if local and record.get("consecutive_qualified_trapped_common_events") != 0:
                raise Proto18HistoricalReplayError("PROTO12: qualified-event counter differs")
            _common_time(record, local, label="PROTO12")
        offset += count
    selected = records[48]
    selected_assessment = selected["assessment"]
    comparator = selected_assessment.get("comparator_common_event")
    if (
        selected.get("amplitude") != "3" or selected.get("event_index") != 23
        or not isinstance(comparator, Mapping) or comparator.get("coordinate_time") != 23 / 16
        or comparator.get("admission_passed") is not False
    ):
        raise Proto18HistoricalReplayError("PROTO12: selected amplitude-3 event 23 differs")
    completed = metadata.get("completed_amplitude_records")
    if (not isinstance(completed, Sequence)
            or isinstance(completed, (str, bytes, bytearray))
            or len(completed) != 2):
        raise Proto18HistoricalReplayError("PROTO12: completed-amplitude metadata differs")
    selected_meta = completed[1]
    stop = selected_meta.get("stop") if isinstance(selected_meta, Mapping) else None
    if not isinstance(selected_meta, Mapping) or not isinstance(stop, Mapping) or (
        selected_meta.get("amplitude") != "3"
        or selected_meta.get("last_completed_common_event_index") != 23
        or selected_meta.get("last_completed_coordinate_time") != 23 / 16
        or stop.get("classification") != "common_event_constraint_or_spectral_stop"
    ):
        raise Proto18HistoricalReplayError("PROTO12: selected checkpoint completion differs")
    if metadata.get("amplitude") != "3" or metadata.get("completed_common_event_index") != 23:
        raise Proto18HistoricalReplayError("PROTO12: checkpoint active completion differs")
    members = metadata.get("members")
    final_hashes = selected_meta.get("member_final_state_hashes") if isinstance(selected_meta, Mapping) else None
    diagnostics = selected_assessment.get("member_diagnostics")
    expected_selected = {"RK4-2049", "RK4-4097", "RK4-8193", "SSPRK3-4097", "SSPRK3-8193"}
    if (not isinstance(members, Mapping) or not isinstance(final_hashes, Mapping)
            or not isinstance(diagnostics, Mapping)
            or set(selected_physical_state_sha256) != expected_selected):
        raise Proto18HistoricalReplayError("PROTO12: selected member cross-links are absent")
    for key in sorted(expected_selected):
        member = members.get(key)
        diagnostic = diagnostics.get(key)
        expected_hash = selected_physical_state_sha256[key]
        if (not isinstance(member, Mapping) or not isinstance(diagnostic, Mapping)
                or member.get("time") != 23 / 16
                or final_hashes.get(key) != expected_hash
                or diagnostic.get("state_sha256") != expected_hash
                or diagnostic.get("step_index") != member.get("step_index")
                or diagnostic.get("transaction_serial") != member.get("transaction_serial")
                or diagnostic.get("source_retry_count") != member.get("source_retry_count")
                or diagnostic.get("CFL_retry_count") != member.get("CFL_retry_count")
                or diagnostic.get("accepted_stage_count")
                != member.get("runtime_monitor_state", {}).get("accepted_stage_count")):
            raise Proto18HistoricalReplayError("PROTO12: selected member cross-link differs")
    return {
        "raw_sha256": PROTO12_EVENT_SHA256,
        "byte_count": PROTO12_EVENT_BYTES,
        "line_count": PROTO12_EVENT_LINES,
        "embedded_external_bytes_equal": True,
        "selected_amplitude": "3",
        "selected_event_index": 23,
        "selected_coordinate_time": "23/16",
    }


def replay_rsp2_history(
    external: bytes, embedded: bytes, metadata: Mapping[str, Any],
    selected_physical_state_sha256: str,
) -> dict[str, Any]:
    """Verify RSP2's initial, 23-boundary, and endpoint record grammar."""
    _require_embedded_equal(external, embedded, label="RSP2")
    records = parse_jsonl_exact(
        external, label="RSP2", expected_bytes=RSP2_EVENT_BYTES,
        expected_sha256=RSP2_EVENT_SHA256, expected_lines=RSP2_EVENT_LINES,
    )
    if not isinstance(metadata, Mapping):
        raise Proto18HistoricalReplayError("RSP2: checkpoint metadata is not an object")
    metadata = _plain(metadata)
    if not isinstance(metadata, Mapping):  # narrowed for static readers
        raise Proto18HistoricalReplayError("RSP2: checkpoint metadata copy failed")
    _event_identity(metadata, external, RSP2_EVENT_LINES, label="RSP2")
    initial = records[0]
    if set(initial) != {"boundary_index", "coordinate_time", "event_type", "member", "trajectory_advanced"} or initial.get("event_type") != "initial_premise_event" or initial.get("boundary_index") != 0 or initial.get("coordinate_time") != 0.0 or initial.get("trajectory_advanced") is not False:
        raise Proto18HistoricalReplayError("RSP2: initial record differs")
    for index, record in enumerate(records[1:24], start=1):
        member = record.get("member")
        if (
            set(record) != {"boundary_index", "coordinate_time", "event_type", "member"}
            or record.get("event_type") != "accepted_progress_boundary"
            or record.get("boundary_index") != index or record.get("coordinate_time") != index / 16.0
            or not isinstance(member, Mapping) or member.get("time") != index / 16.0
            or member.get("source_retry_count") != 0 or member.get("CFL_retry_count") != 0
        ):
            raise Proto18HistoricalReplayError("RSP2: accepted-boundary grammar/order differs")
    final_progress = records[23]
    final_progress_member = final_progress.get("member")
    terminal = records[24]
    assessment = terminal.get("assessment")
    if (
        set(terminal) != {"assessment", "boundary_index", "classification", "coordinate_time", "event_type"}
        or terminal.get("event_type") != "evolved_RSP2_endpoint" or terminal.get("boundary_index") != 23
        or terminal.get("coordinate_time") != 23 / 16 or terminal.get("classification") != "completed_target_and_complete_constraint_pass"
        or not isinstance(assessment, Mapping) or assessment.get("classification") != terminal.get("classification")
    ):
        raise Proto18HistoricalReplayError("RSP2: endpoint grammar differs")
    member = metadata.get("members", {}).get("SSPRK3-16385") if isinstance(metadata.get("members"), Mapping) else None
    result = metadata.get("terminal_result")
    if (
        metadata.get("terminal") is not True or metadata.get("completed_common_event_index") != 23
        or not isinstance(member, Mapping) or member.get("time") != 23 / 16
        or member.get("step_index") != 1771 or member.get("transaction_serial") != 7084
        or member.get("state_sha256") != selected_physical_state_sha256
        or not isinstance(final_progress_member, Mapping)
        or final_progress_member.get("state_sha256") != selected_physical_state_sha256
        or final_progress_member.get("step_index") != member.get("step_index")
        or final_progress_member.get("transaction_serial") != member.get("transaction_serial")
        or not isinstance(result, Mapping) or result.get("endpoint_assessment") != assessment
        or result.get("classification") != terminal.get("classification")
    ):
        raise Proto18HistoricalReplayError("RSP2: checkpoint endpoint cross-link differs")
    return {
        "raw_sha256": RSP2_EVENT_SHA256,
        "byte_count": RSP2_EVENT_BYTES,
        "line_count": RSP2_EVENT_LINES,
        "embedded_external_bytes_equal": True,
        "endpoint_boundary_index": 23,
        "endpoint_coordinate_time": "23/16",
        "classification": "completed_target_and_complete_constraint_pass",
    }


__all__ = [
    "PROTO12_EVENT_BYTES", "PROTO12_EVENT_LINES", "PROTO12_EVENT_SHA256",
    "RSP2_EVENT_BYTES", "RSP2_EVENT_LINES", "RSP2_EVENT_SHA256",
    "Proto18HistoricalReplayError", "parse_jsonl_exact", "replay_proto12_history",
    "replay_rsp2_history",
]
