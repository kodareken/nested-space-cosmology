"""Immutable CAL9 checkpoint adapter for the RSP2 late-event study.

RSP2 does not evolve its two predecessor resolutions.  This module loads only
the terminal SSPRK3 4,097- and 8,193-point arrays from the hash-bound CAL9
checkpoint, verifies their campaign metadata and content hashes, and exposes
them as read-only evidence for the pure constraint classifier.
"""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
import json
from pathlib import Path
from typing import Mapping

import numpy as np

from .numerical_engine import (
    EvolutionState,
    UniformRadialGrid,
    array_content_sha256,
)


RSP2_CAL9_CAMPAIGN_ID = (
    "32ac798e8bd7ccb984dd66262439d55ec01ad7c5cabfea785fe30d1fc60eff36"
)
RSP2_CAL9_EVENT_LOG_SHA256 = (
    "c2d35d07db2f6c983a449cdcf2db469360565cf392d78348dcafbd1d81234ba7"
)
RSP2_CAL9_POINT_COUNTS = (4097, 8193)
RSP2_CAL9_TARGET_TIME = 23.0 / 16.0
RSP2_CAL9_STAGE_COUNTS = (2524, 3948)
RSP2_CAL9_STATE_HASHES = (
    "ab96602c0f05633375710efee2be01be84e24c32a31d07d6ef23390788a72162",
    "c46abc3f0f72fa2a13c1fee7bb329b94698cc2c26003406ae3eb1592f11b648c",
)


def file_sha256(path: Path) -> str:
    """Return a streaming SHA-256 digest without mutating the source file."""

    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _readonly(value: object) -> np.ndarray:
    answer = np.ascontiguousarray(np.asarray(value, dtype=np.float64)).copy()
    if answer.ndim != 2 or answer.shape[1] != 6 or not np.all(np.isfinite(answer)):
        raise ValueError("CAL9 predecessor arrays must be finite (N, 6) binary64")
    answer.setflags(write=False)
    return answer


@dataclass(frozen=True, slots=True)
class RSP2CAL9Predecessors:
    """Verified read-only terminal evidence consumed by RSP2."""

    checkpoint_sha256: str
    campaign_id: str
    coordinate_time: float
    event_log_sha256: str
    event_log_line_count: int
    point_counts: tuple[int, int]
    accepted_stage_counts: tuple[int, int]
    input_hashes: tuple[str, str]
    state_hashes: tuple[str, str]
    states: tuple[EvolutionState, EvolutionState]
    grids: tuple[UniformRadialGrid, UniformRadialGrid]


def _terminal_amplitude_record(metadata: Mapping[str, object]) -> Mapping[str, object]:
    records = metadata.get("completed_amplitude_records")
    if not isinstance(records, list):
        raise ValueError("CAL9 checkpoint amplitude records are absent")
    selected = [item for item in records if isinstance(item, dict) and item.get("amplitude") == "3"]
    if len(selected) != 1:
        raise ValueError("CAL9 checkpoint does not contain one amplitude-three record")
    return selected[0]


def load_rsp2_cal9_predecessors(
    path: Path,
    *,
    expected_sha256: str,
    outer_radius: float = 128.0,
) -> RSP2CAL9Predecessors:
    """Load and fully bind the two CAL9 predecessor states required by RSP2."""

    if not path.is_file():
        raise FileNotFoundError(path)
    observed_sha = file_sha256(path)
    if observed_sha != expected_sha256:
        raise ValueError("CAL9 checkpoint hash differs from the RSP2 freeze")

    with np.load(path, allow_pickle=False) as archive:
        metadata = json.loads(bytes(archive["metadata_utf8"]).decode("utf-8"))
        if not isinstance(metadata, dict):
            raise ValueError("CAL9 checkpoint metadata must be an object")
        event_log_payload = bytes(archive["event_log_utf8"])
        observed_event_log_sha256 = sha256(event_log_payload).hexdigest()
        observed_event_log_line_count = event_log_payload.count(b"\n")
        if event_log_payload and not event_log_payload.endswith(b"\n"):
            raise ValueError("CAL9 checkpoint event log is not newline terminated")
        amplitude_record = _terminal_amplitude_record(metadata)
        member_metadata = metadata.get("members")
        if not isinstance(member_metadata, dict):
            raise ValueError("CAL9 checkpoint member metadata are absent")

        states: list[EvolutionState] = []
        grids: list[UniformRadialGrid] = []
        stage_counts: list[int] = []
        input_hashes: list[str] = []
        state_hashes: list[str] = []
        expected_record_hashes = amplitude_record.get("member_final_state_hashes")
        if not isinstance(expected_record_hashes, dict):
            raise ValueError("CAL9 checkpoint terminal state hashes are absent")

        for point_count, expected_stage_count, expected_state_hash in zip(
            RSP2_CAL9_POINT_COUNTS,
            RSP2_CAL9_STAGE_COUNTS,
            RSP2_CAL9_STATE_HASHES,
            strict=True,
        ):
            key = f"SSPRK3-{point_count}"
            stored = member_metadata.get(key)
            if not isinstance(stored, dict):
                raise ValueError(f"CAL9 checkpoint member is absent: {key}")
            runtime_state = stored.get("runtime_monitor_state")
            if not isinstance(runtime_state, dict):
                raise ValueError(f"CAL9 runtime monitor is absent: {key}")
            if (
                np.float64(stored.get("time")).tobytes()
                != np.float64(RSP2_CAL9_TARGET_TIME).tobytes()
                or runtime_state.get("accepted_stage_count") != expected_stage_count
                or stored.get("source_retry_count") != 0
            ):
                raise ValueError(f"CAL9 terminal metadata differ: {key}")

            prefix = key.replace("-", "_")
            state = EvolutionState(
                _readonly(archive[f"{prefix}_u"]),
                _readonly(archive[f"{prefix}_p"]),
                _readonly(archive[f"{prefix}_q"]),
            )
            observed_state_hash = array_content_sha256(state.u, state.p, state.q)
            if (
                observed_state_hash != expected_state_hash
                or expected_record_hashes.get(key) != expected_state_hash
            ):
                raise ValueError(f"CAL9 terminal state hash differs: {key}")
            input_hash = stored.get("input_hash")
            if not isinstance(input_hash, str) or len(input_hash) != 64:
                raise ValueError(f"CAL9 input hash differs: {key}")
            states.append(state)
            grids.append(UniformRadialGrid(0.0, outer_radius, point_count))
            stage_counts.append(expected_stage_count)
            input_hashes.append(input_hash)
            state_hashes.append(observed_state_hash)

    event_log = metadata.get("event_log")
    amplitude_stop = amplitude_record.get("stop")
    source_retries = amplitude_record.get("member_source_retry_counts")
    if (
        metadata.get("schema_version") != 1
        or metadata.get("terminal") is not True
        or metadata.get("amplitude") != "3"
        or metadata.get("campaign_id") != RSP2_CAL9_CAMPAIGN_ID
        or metadata.get("completed_common_event_index") != 23
        or not isinstance(event_log, dict)
        or event_log.get("sha256") != RSP2_CAL9_EVENT_LOG_SHA256
        or event_log.get("line_count") != 49
        or observed_event_log_sha256 != RSP2_CAL9_EVENT_LOG_SHA256
        or observed_event_log_line_count != 49
        or amplitude_record.get("eligible") is not False
        or amplitude_record.get("last_completed_common_event_index") != 23
        or amplitude_record.get("last_completed_coordinate_time")
        != RSP2_CAL9_TARGET_TIME
        or not isinstance(source_retries, dict)
        or source_retries.get("SSPRK3-4097") != 0
        or source_retries.get("SSPRK3-8193") != 0
        or not isinstance(amplitude_stop, dict)
        or amplitude_stop.get("classification")
        != "common_event_constraint_or_spectral_stop"
        or amplitude_stop.get("coordinate_time") != RSP2_CAL9_TARGET_TIME
        or amplitude_stop.get("event_index") != 23
        or amplitude_stop.get("reason")
        != "common_event_constraint_or_spatial_spectral_admission"
    ):
        raise ValueError("CAL9 checkpoint campaign identity differs")

    return RSP2CAL9Predecessors(
        checkpoint_sha256=observed_sha,
        campaign_id=RSP2_CAL9_CAMPAIGN_ID,
        coordinate_time=RSP2_CAL9_TARGET_TIME,
        event_log_sha256=RSP2_CAL9_EVENT_LOG_SHA256,
        event_log_line_count=49,
        point_counts=RSP2_CAL9_POINT_COUNTS,
        accepted_stage_counts=tuple(stage_counts),
        input_hashes=tuple(input_hashes),
        state_hashes=tuple(state_hashes),
        states=tuple(states),
        grids=tuple(grids),
    )


__all__ = [
    "RSP2CAL9Predecessors",
    "RSP2_CAL9_CAMPAIGN_ID",
    "RSP2_CAL9_EVENT_LOG_SHA256",
    "RSP2_CAL9_POINT_COUNTS",
    "RSP2_CAL9_STAGE_COUNTS",
    "RSP2_CAL9_STATE_HASHES",
    "RSP2_CAL9_TARGET_TIME",
    "file_sha256",
    "load_rsp2_cal9_predecessors",
]
