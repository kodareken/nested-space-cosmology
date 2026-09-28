"""No-store child-process isolation for prospective PRO20 RSRC1.

The parent owns authority, scheduling, and any later durable publication.  A
fresh child reconstructs one scheduled physical member, restores one accepted
bundle, and returns immutable bytes.  The parent derives the protocol class
again from predecessor/successor bytes and typed evidence; it never trusts a
child-supplied kind label.  Child crashes, signals, timeouts, or malformed
handoffs are forensic uncertainty and cannot be converted into a resource or
scientific terminal.
"""

from __future__ import annotations

from base64 import b64decode, b64encode
from dataclasses import dataclass
from hashlib import sha256
import json
import os
from pathlib import Path
import resource
import subprocess
import time
from types import MappingProxyType
from typing import Callable, Mapping, Sequence

from recursive_horizons.evidence_io import (
    CanonicalJSONError,
    canonical_json_bytes,
    load_canonical_json,
    read_regular_file,
)

from .hlt17_member_checkpoint import HLT17GenerationBundle, HLT17InMemoryCheckpoint
from .pro20_rsrc1_attempt import (
    AttemptClassification,
    RSRC1AttemptPremiseError,
    classify_attempt,
    execute_scheduled_attempt,
    restore_checkpoint_from_bundle,
)
from .pro20_origin import MEMBER_KEYS, capture_pro20_historical_origin
from .pro20_rsrc1_member import build_rsrc1_member
from .proto19_gr0_static_factory import STATIC_INPUT_PATHS


ARTIFACT_ID = "FGC-1-PRO20-EV1-RSRC1"
REQUEST_SCHEMA = "FGC-1-PRO20-EV1-RSRC1-child-request-v1"
RESPONSE_SCHEMA = "FGC-1-PRO20-EV1-RSRC1-child-response-v1"
CHILD_MAX_PEAK_RSS_BYTES = 4 * 1024**3
PARENT_MAX_CURRENT_RSS_BYTES = 1 * 1024**3
HOST_RESERVE_BYTES = 4 * 1024**3
MIN_AVAILABLE_MEMORY_BYTES = (
    CHILD_MAX_PEAK_RSS_BYTES + PARENT_MAX_CURRENT_RSS_BYTES + HOST_RESERVE_BYTES
)
MAX_CHILD_WALL_SECONDS = 86400.0
MAX_HANDOFF_BYTES = 256 * 1024**2
MAX_STDERR_BYTES = 1024**2
CHILD_MONITOR_INTERVAL_SECONDS = 0.25
RESOURCE_POLICY = MappingProxyType({
    "child_max_peak_rss_bytes": CHILD_MAX_PEAK_RSS_BYTES,
    "parent_max_current_rss_bytes": PARENT_MAX_CURRENT_RSS_BYTES,
    "host_reserve_bytes": HOST_RESERVE_BYTES,
    "min_available_memory_bytes": MIN_AVAILABLE_MEMORY_BYTES,
    "max_child_wall_seconds": MAX_CHILD_WALL_SECONDS,
    "max_handoff_bytes": MAX_HANDOFF_BYTES,
    "max_stderr_bytes": MAX_STDERR_BYTES,
    "child_monitor_interval_seconds": CHILD_MONITOR_INTERVAL_SECONDS,
    "scientific_threshold_fitted": False,
    "fresh_child_per_scheduled_attempt": True,
    "parent_is_sole_future_publisher": True,
})


class PRO20RSRC1IsolationError(ValueError):
    """The isolated request, response, or resource premise differs."""


class PRO20RSRC1ForensicUncertainty(RuntimeError):
    """No valid child terminal exists; a future writer must remain unclosed."""

    publish_typed_terminal = False
    close_writer = False


def _sha(raw: bytes) -> str:
    if type(raw) is not bytes:
        raise PRO20RSRC1IsolationError("handoff payload is not immutable bytes")
    return sha256(raw).hexdigest()


def _require_sha(value: object, label: str) -> str:
    if type(value) is not str or len(value) != 64 or value.lower() != value:
        raise PRO20RSRC1IsolationError(f"{label} is not lowercase SHA-256")
    try:
        int(value, 16)
    except ValueError as error:
        raise PRO20RSRC1IsolationError(f"{label} is not lowercase SHA-256") from error
    return value


def _require_member_key(value: object) -> str:
    if type(value) is not str or value not in MEMBER_KEYS:
        raise PRO20RSRC1IsolationError("handoff member is outside the live cohort")
    return value


def _exact_mapping(
    value: object,
    keys: frozenset[str],
    label: str,
) -> dict[str, object]:
    if type(value) is not dict or set(value) != keys:
        raise PRO20RSRC1IsolationError(f"{label} fields differ")
    return dict(value)


def _parse_canonical(raw: bytes, label: str) -> dict[str, object]:
    if type(raw) is not bytes or len(raw) > MAX_HANDOFF_BYTES:
        raise PRO20RSRC1IsolationError(f"{label} exceeds the handoff boundary")
    try:
        value = load_canonical_json(raw)
    except (
        CanonicalJSONError,
        UnicodeDecodeError,
        json.JSONDecodeError,
        TypeError,
        ValueError,
    ) as error:
        raise PRO20RSRC1IsolationError(f"{label} is malformed canonical JSON") from error
    if type(value) is not dict or canonical_json_bytes(value) != raw:
        raise PRO20RSRC1IsolationError(f"{label} is not canonical JSON")
    return value


def _bundle_mapping(bundle: HLT17GenerationBundle) -> dict[str, object]:
    if not isinstance(bundle, HLT17GenerationBundle):
        raise TypeError("bundle must be HLT17GenerationBundle")
    bundle.__post_init__()
    return {
        "cursor_b64": b64encode(bundle.cursor_bytes).decode("ascii"),
        "cursor_sha256": bundle.cursor_sha256,
        "descriptor_b64": b64encode(bundle.descriptor).decode("ascii"),
        "descriptor_sha256": bundle.descriptor_sha256,
        "payload_b64": b64encode(bundle.payload).decode("ascii"),
        "payload_sha256": bundle.payload_sha256,
        "physical_state_sha256": bundle.physical_state_sha256,
    }


_BUNDLE_KEYS = frozenset(
    {
        "cursor_b64",
        "cursor_sha256",
        "descriptor_b64",
        "descriptor_sha256",
        "payload_b64",
        "payload_sha256",
        "physical_state_sha256",
    }
)


def _decode_b64(value: object, label: str) -> bytes:
    if type(value) is not str:
        raise PRO20RSRC1IsolationError(f"{label} is not base64 text")
    try:
        raw = b64decode(value.encode("ascii"), validate=True)
    except (UnicodeEncodeError, ValueError) as error:
        raise PRO20RSRC1IsolationError(f"{label} is not canonical base64") from error
    if b64encode(raw).decode("ascii") != value:
        raise PRO20RSRC1IsolationError(f"{label} is not canonical base64")
    return raw


def _bundle_from_mapping(value: object) -> HLT17GenerationBundle:
    mapping = _exact_mapping(value, _BUNDLE_KEYS, "bundle")
    descriptor = _decode_b64(mapping["descriptor_b64"], "descriptor")
    payload = _decode_b64(mapping["payload_b64"], "payload")
    cursor = _decode_b64(mapping["cursor_b64"], "cursor")
    if len(descriptor) + len(payload) + len(cursor) > MAX_HANDOFF_BYTES:
        raise PRO20RSRC1IsolationError("bundle exceeds the handoff boundary")
    try:
        bundle = HLT17GenerationBundle(descriptor, payload, cursor)
    except (TypeError, ValueError) as error:
        raise PRO20RSRC1IsolationError("bundle failed independent reauthentication") from error
    expected = _bundle_mapping(bundle)
    if mapping != expected:
        raise PRO20RSRC1IsolationError("bundle hashes or bytes differ")
    return bundle


def encode_bundle(bundle: HLT17GenerationBundle) -> dict[str, object]:
    """Public immutable bundle wire used by RSRC1 seed and attempt children."""

    return _bundle_mapping(bundle)


def decode_bundle(value: object) -> HLT17GenerationBundle:
    """Reauthenticate one RSRC1 bundle wire without opening a store."""

    return _bundle_from_mapping(value)


def parse_canonical_handoff(raw: bytes, label: str) -> dict[str, object]:
    """Parse one bounded exact-canonical RSRC1 handoff object."""

    return _parse_canonical(raw, label)


@dataclass(frozen=True, slots=True)
class IsolatedAttemptRequest:
    member_key: str
    predecessor_bundle: HLT17GenerationBundle
    authority_receipt_sha256: str
    origin_capture_sha256: str
    generation1_bridge_content_id: str
    source_closure_sha256: str
    source_configuration_sha256: str

    def __post_init__(self) -> None:
        _require_member_key(self.member_key)
        self.predecessor_bundle.__post_init__()
        for value, label in (
            (self.authority_receipt_sha256, "authority receipt"),
            (self.origin_capture_sha256, "origin capture"),
            (self.generation1_bridge_content_id, "generation-one bridge"),
            (self.source_closure_sha256, "source closure"),
            (self.source_configuration_sha256, "source configuration"),
        ):
            _require_sha(value, label)

    @property
    def mapping_without_id(self) -> dict[str, object]:
        return {
            "authority_receipt_sha256": self.authority_receipt_sha256,
            "generation1_bridge_content_id": self.generation1_bridge_content_id,
            "member_key": self.member_key,
            "origin_capture_sha256": self.origin_capture_sha256,
            "predecessor_bundle": _bundle_mapping(self.predecessor_bundle),
            "schema": REQUEST_SCHEMA,
            "source_closure_sha256": self.source_closure_sha256,
            "source_configuration_sha256": self.source_configuration_sha256,
        }

    @property
    def request_sha256(self) -> str:
        return _sha(canonical_json_bytes(self.mapping_without_id))

    @property
    def mapping(self) -> dict[str, object]:
        return {**self.mapping_without_id, "request_sha256": self.request_sha256}


_REQUEST_KEYS = frozenset(
    {
        "authority_receipt_sha256",
        "generation1_bridge_content_id",
        "member_key",
        "origin_capture_sha256",
        "predecessor_bundle",
        "request_sha256",
        "schema",
        "source_closure_sha256",
        "source_configuration_sha256",
    }
)


def encode_attempt_request(request: IsolatedAttemptRequest) -> bytes:
    request.__post_init__()
    return canonical_json_bytes(request.mapping)


def decode_attempt_request(raw: bytes) -> IsolatedAttemptRequest:
    mapping = _exact_mapping(_parse_canonical(raw, "child request"), _REQUEST_KEYS, "child request")
    if mapping["schema"] != REQUEST_SCHEMA:
        raise PRO20RSRC1IsolationError("child request schema differs")
    request = IsolatedAttemptRequest(
        member_key=_require_member_key(mapping["member_key"]),
        predecessor_bundle=_bundle_from_mapping(mapping["predecessor_bundle"]),
        authority_receipt_sha256=_require_sha(
            mapping["authority_receipt_sha256"], "authority receipt"
        ),
        origin_capture_sha256=_require_sha(
            mapping["origin_capture_sha256"], "origin capture"
        ),
        generation1_bridge_content_id=_require_sha(
            mapping["generation1_bridge_content_id"], "generation-one bridge"
        ),
        source_closure_sha256=_require_sha(
            mapping["source_closure_sha256"], "source closure"
        ),
        source_configuration_sha256=_require_sha(
            mapping["source_configuration_sha256"], "source configuration"
        ),
    )
    if mapping["request_sha256"] != request.request_sha256:
        raise PRO20RSRC1IsolationError("child request digest differs")
    return request


@dataclass(frozen=True, slots=True)
class ChildMetrics:
    peak_rss_bytes: int
    wall_seconds: float

    def __post_init__(self) -> None:
        if type(self.peak_rss_bytes) is not int or self.peak_rss_bytes < 0:
            raise PRO20RSRC1IsolationError("child peak RSS is invalid")
        if type(self.wall_seconds) is not float or not (
            0.0 <= self.wall_seconds <= MAX_CHILD_WALL_SECONDS
        ):
            raise PRO20RSRC1IsolationError("child wall observation is invalid")


_RESPONSE_KEYS = frozenset(
    {
        "accepted_state_advanced",
        "child_peak_rss_bytes",
        "child_wall_seconds_hex",
        "member_key",
        "request_sha256",
        "schema",
        "stop_evidence",
        "stop_kind",
        "successor_bundle",
    }
)


def _stop_payload(classification: AttemptClassification) -> tuple[str, object]:
    evidence = classification.terminal_evidence
    if classification.kind == "resource_exhausted":
        return "resource", dict(evidence or {})
    if classification.kind == "invalid_premise":
        return "premise", dict(evidence or {})
    if classification.kind == "temporal_retry_exhausted":
        return "temporal", dict(evidence or {})
    return "none", None


def encode_child_response(
    request: IsolatedAttemptRequest,
    classification: AttemptClassification,
    *,
    peak_rss_bytes: int,
    wall_seconds: float,
) -> bytes:
    """Encode no labels the parent cannot independently derive."""

    request.__post_init__()
    classification.__post_init__()
    metrics = ChildMetrics(int(peak_rss_bytes), float(wall_seconds))
    if classification.member_key != request.member_key:
        raise PRO20RSRC1IsolationError("child classification member differs")
    if metrics.peak_rss_bytes > CHILD_MAX_PEAK_RSS_BYTES:
        classification = classify_attempt(
            member_key=request.member_key,
            predecessor_bundle=request.predecessor_bundle,
            successor_bundle=request.predecessor_bundle,
            accepted_state_advanced=False,
            resource_stop={
                "reason": "isolated child peak RSS exceeded",
                "observed_peak_rss_bytes": metrics.peak_rss_bytes,
                "limit_bytes": CHILD_MAX_PEAK_RSS_BYTES,
                "scope": "one_scheduled_attempt",
            },
        )
    stop_kind, stop_evidence = _stop_payload(classification)
    response = {
        "accepted_state_advanced": classification.accepted_state_advanced,
        "child_peak_rss_bytes": metrics.peak_rss_bytes,
        "child_wall_seconds_hex": metrics.wall_seconds.hex(),
        "member_key": request.member_key,
        "request_sha256": request.request_sha256,
        "schema": RESPONSE_SCHEMA,
        "stop_evidence": stop_evidence,
        "stop_kind": stop_kind,
        "successor_bundle": _bundle_mapping(classification.successor_bundle),
    }
    return canonical_json_bytes(response)


def classify_child_response(
    request: IsolatedAttemptRequest,
    raw: bytes,
) -> tuple[AttemptClassification, ChildMetrics]:
    """Reauthenticate bytes and independently reduce the child response."""

    request.__post_init__()
    mapping = _exact_mapping(
        _parse_canonical(raw, "child response"), _RESPONSE_KEYS, "child response"
    )
    if mapping["schema"] != RESPONSE_SCHEMA:
        raise PRO20RSRC1IsolationError("child response schema differs")
    if mapping["request_sha256"] != request.request_sha256:
        raise PRO20RSRC1IsolationError("child response request digest differs")
    if mapping["member_key"] != request.member_key:
        raise PRO20RSRC1IsolationError("child response member differs")
    if type(mapping["accepted_state_advanced"]) is not bool:
        raise PRO20RSRC1IsolationError("child accepted-state flag differs")
    if type(mapping["child_peak_rss_bytes"]) is not int:
        raise PRO20RSRC1IsolationError("child peak RSS type differs")
    if type(mapping["child_wall_seconds_hex"]) is not str:
        raise PRO20RSRC1IsolationError("child wall type differs")
    try:
        wall = float.fromhex(mapping["child_wall_seconds_hex"])
    except ValueError as error:
        raise PRO20RSRC1IsolationError("child wall value differs") from error
    metrics = ChildMetrics(mapping["child_peak_rss_bytes"], wall)
    successor = _bundle_from_mapping(mapping["successor_bundle"])
    stop_kind = mapping["stop_kind"]
    evidence = mapping["stop_evidence"]
    if metrics.peak_rss_bytes > CHILD_MAX_PEAK_RSS_BYTES and stop_kind != "resource":
        raise PRO20RSRC1IsolationError(
            "over-ceiling child did not return a resource stop"
        )
    if stop_kind != "none" and mapping["accepted_state_advanced"] is not False:
        raise PRO20RSRC1IsolationError(
            "typed child stop advertised accepted-state advance"
        )
    common = {
        "member_key": request.member_key,
        "predecessor_bundle": request.predecessor_bundle,
        "successor_bundle": successor,
        "accepted_state_advanced": mapping["accepted_state_advanced"],
    }
    if stop_kind == "none":
        if evidence is not None:
            raise PRO20RSRC1IsolationError("unstopped child supplied stop evidence")
        classification = classify_attempt(**common)
    elif stop_kind == "resource":
        if type(evidence) is not dict:
            raise PRO20RSRC1IsolationError("resource child omitted evidence")
        classification = classify_attempt(**common, resource_stop=evidence)
    elif stop_kind == "premise":
        if type(evidence) is not dict:
            raise PRO20RSRC1IsolationError("premise child omitted evidence")
        classification = classify_attempt(**common, premise_stop=evidence)
    elif stop_kind == "temporal":
        if type(evidence) is not dict:
            raise PRO20RSRC1IsolationError("temporal child omitted evidence")
        classification = classify_attempt(
            **common, temporal_exhaustion_evidence=evidence
        )
    else:
        raise PRO20RSRC1IsolationError("child stop kind differs")
    return classification, metrics


def current_child_peak_rss_bytes() -> int:
    rss = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    if rss < 0:
        raise PRO20RSRC1IsolationError("child peak RSS is negative")
    return rss if os.uname().sysname == "Darwin" else rss * 1024


def child_process_rss_bytes(pid: int) -> int:
    if type(pid) is not int or pid <= 0:
        raise PRO20RSRC1IsolationError("child pid is invalid")
    try:
        completed = subprocess.run(
            ("ps", "-o", "rss=", "-p", str(pid)),
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
            timeout=5.0,
        )
    except (OSError, subprocess.TimeoutExpired) as error:
        raise PRO20RSRC1IsolationError("child current RSS cannot be observed") from error
    value = completed.stdout.strip()
    if completed.returncode != 0 or not value.isdecimal():
        raise PRO20RSRC1IsolationError("child current RSS cannot be observed")
    return int(value) * 1024


def _terminate_child(process: subprocess.Popen[bytes]) -> None:
    try:
        process.terminate()
    except ProcessLookupError:
        return
    try:
        process.wait(timeout=5.0)
    except subprocess.TimeoutExpired:
        try:
            process.kill()
        except ProcessLookupError:
            return
        process.wait()


def check_parent_resource_premises(
    *,
    parent_current_rss_bytes: int,
    available_memory_bytes: int,
) -> None:
    if type(parent_current_rss_bytes) is not int or parent_current_rss_bytes < 0:
        raise PRO20RSRC1IsolationError("parent current RSS is invalid")
    if type(available_memory_bytes) is not int or available_memory_bytes < 0:
        raise PRO20RSRC1IsolationError("available memory is invalid")
    if parent_current_rss_bytes > PARENT_MAX_CURRENT_RSS_BYTES:
        raise PRO20RSRC1IsolationError("parent current RSS exceeds RSRC1 ceiling")
    if available_memory_bytes < MIN_AVAILABLE_MEMORY_BYTES:
        raise PRO20RSRC1IsolationError("available memory is below RSRC1 reserve")


def _load_static_inputs(root: Path) -> dict[str, bytes]:
    return {relative: read_regular_file(root, relative) for relative in STATIC_INPUT_PATHS}


def execute_child_request(
    repository_root: Path,
    request_raw: bytes,
    *,
    origin_loader: Callable[[Path], object] = capture_pro20_historical_origin,
    static_loader: Callable[[Path], Mapping[str, bytes]] = _load_static_inputs,
    member_builder: Callable[..., object] = build_rsrc1_member,
    step_fn: Callable[[HLT17InMemoryCheckpoint], AttemptClassification] = (
        execute_scheduled_attempt
    ),
    peak_rss: Callable[[], int] = current_child_peak_rss_bytes,
    now: Callable[[], float] = time.monotonic,
) -> bytes:
    """Physical child seam. It reads no campaign store and writes nothing."""

    root = Path(repository_root)
    if not root.is_absolute() or root.resolve(strict=True) != root:
        raise PRO20RSRC1IsolationError("child repository root is not canonical")
    request = decode_attempt_request(request_raw)
    started = float(now())
    origin = origin_loader(root)
    static_inputs = static_loader(root)
    record = member_builder(
        origin,
        member_key=request.member_key,
        static_input_bytes=static_inputs,
        source_closure_sha256=request.source_closure_sha256,
    )
    for observed, expected, label in (
        (getattr(record, "captured_origin_sha256", None), request.origin_capture_sha256, "origin"),
        (
            getattr(record, "generation1_bridge_content_id", None),
            request.generation1_bridge_content_id,
            "generation-one bridge",
        ),
        (getattr(record, "source_closure_sha256", None), request.source_closure_sha256, "source closure"),
        (
            getattr(record, "source_configuration_sha256", None),
            request.source_configuration_sha256,
            "source configuration",
        ),
    ):
        if observed != expected:
            raise PRO20RSRC1IsolationError(f"child {label} identity differs")
    checkpoint = getattr(record, "checkpoint", None)
    if not isinstance(checkpoint, HLT17InMemoryCheckpoint):
        raise PRO20RSRC1IsolationError("child member omitted its live checkpoint")
    restore_checkpoint_from_bundle(checkpoint, request.predecessor_bundle)
    try:
        classification = step_fn(checkpoint)
    except RSRC1AttemptPremiseError as error:
        classification = classify_attempt(
            member_key=request.member_key,
            predecessor_bundle=request.predecessor_bundle,
            successor_bundle=request.predecessor_bundle,
            accepted_state_advanced=False,
            premise_stop={"reason": str(error) or "invalid_premise"},
        )
    elapsed = float(now()) - started
    return encode_child_response(
        request,
        classification,
        peak_rss_bytes=peak_rss(),
        wall_seconds=elapsed,
    )


@dataclass(frozen=True, slots=True)
class CompletedChild:
    returncode: int
    stdout: bytes
    stderr: bytes
    wall_seconds: float

    def __post_init__(self) -> None:
        if type(self.returncode) is not int:
            raise TypeError("child returncode must be an integer")
        if type(self.stdout) is not bytes or type(self.stderr) is not bytes:
            raise TypeError("child streams must be immutable bytes")
        if type(self.wall_seconds) is not float or self.wall_seconds < 0.0:
            raise TypeError("child wall time must be nonnegative float")


def launch_child_process(
    command: Sequence[str],
    *,
    repository_root: Path,
    request_raw: bytes,
    timeout_seconds: float = MAX_CHILD_WALL_SECONDS,
    child_rss: Callable[[int], int] = child_process_rss_bytes,
    now: Callable[[], float] = time.monotonic,
) -> CompletedChild:
    """Launch one fresh child. Any abnormal exit is handled by the caller."""

    if (
        not isinstance(command, (tuple, list))
        or not command
        or any(type(item) is not str or not item for item in command)
    ):
        raise PRO20RSRC1IsolationError("child command is invalid")
    root = Path(repository_root)
    if not root.is_absolute() or root.resolve(strict=True) != root:
        raise PRO20RSRC1IsolationError("child working directory is not canonical")
    if type(request_raw) is not bytes or len(request_raw) > MAX_HANDOFF_BYTES:
        raise PRO20RSRC1IsolationError("child request exceeds the handoff boundary")
    if (
        type(timeout_seconds) is not float
        or not 0.0 < timeout_seconds <= MAX_CHILD_WALL_SECONDS
    ):
        raise PRO20RSRC1IsolationError("child timeout differs from the resource policy")
    started = float(now())
    try:
        process = subprocess.Popen(
            tuple(command),
            cwd=root,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            close_fds=True,
            start_new_session=True,
        )
    except OSError as error:
        raise PRO20RSRC1ForensicUncertainty("isolated child could not start") from error
    first = True
    deadline = started + timeout_seconds
    while True:
        remaining = deadline - float(now())
        if remaining <= 0.0:
            _terminate_child(process)
            raise PRO20RSRC1ForensicUncertainty(
                "isolated child exceeded its wall allowance"
            )
        try:
            stdout, stderr = process.communicate(
                input=request_raw if first else None,
                timeout=min(CHILD_MONITOR_INTERVAL_SECONDS, remaining),
            )
            break
        except subprocess.TimeoutExpired:
            first = False
            if process.poll() is not None:
                continue
            try:
                observed = child_rss(int(process.pid))
            except (OSError, TypeError, ValueError) as error:
                _terminate_child(process)
                raise PRO20RSRC1ForensicUncertainty(
                    "isolated child current RSS could not be observed"
                ) from error
            if observed > CHILD_MAX_PEAK_RSS_BYTES:
                _terminate_child(process)
                raise PRO20RSRC1ForensicUncertainty(
                    "isolated child crossed the hard current-RSS ceiling"
                )
    return CompletedChild(
        returncode=int(process.returncode),
        stdout=bytes(stdout),
        stderr=bytes(stderr),
        wall_seconds=float(now()) - started,
    )


def reduce_completed_child(
    request: IsolatedAttemptRequest,
    completed: CompletedChild,
) -> tuple[AttemptClassification, ChildMetrics]:
    """Accept only one silent, bounded, zero-exit canonical response."""

    request.__post_init__()
    completed.__post_init__()
    if completed.returncode != 0:
        raise PRO20RSRC1ForensicUncertainty(
            "isolated child exited without a valid terminal"
        )
    if completed.stderr:
        if len(completed.stderr) > MAX_STDERR_BYTES:
            raise PRO20RSRC1ForensicUncertainty("isolated child stderr was unbounded")
        raise PRO20RSRC1ForensicUncertainty("isolated child emitted stderr")
    if len(completed.stdout) > MAX_HANDOFF_BYTES:
        raise PRO20RSRC1ForensicUncertainty("isolated child response was unbounded")
    if completed.wall_seconds > MAX_CHILD_WALL_SECONDS:
        raise PRO20RSRC1ForensicUncertainty(
            "isolated child exceeded its wall allowance"
        )
    try:
        return classify_child_response(request, completed.stdout)
    except (PRO20RSRC1IsolationError, TypeError, ValueError) as error:
        raise PRO20RSRC1ForensicUncertainty(
            "isolated child response failed parent reauthentication"
        ) from error


def plan_status() -> dict[str, object]:
    return {
        "artifact_id": ARTIFACT_ID,
        "classification": "prospective_resource_isolation_core_no_authority",
        "resource_policy": dict(RESOURCE_POLICY),
        "child_reads_campaign_store": False,
        "child_writes_campaign_store": False,
        "parent_publication_implemented": False,
        "fresh_namespace_selected": False,
        "freeze_implemented": False,
        "live_target_present": False,
        "campaign_execution_authorized": False,
        "resume_authorized": False,
        "takeover_authorized": False,
        "calibration_eligible": False,
        "mechanism_claimed": False,
        "physics_claimed": False,
    }


__all__ = [
    "ARTIFACT_ID",
    "CHILD_MAX_PEAK_RSS_BYTES",
    "CHILD_MONITOR_INTERVAL_SECONDS",
    "ChildMetrics",
    "CompletedChild",
    "HOST_RESERVE_BYTES",
    "IsolatedAttemptRequest",
    "MAX_CHILD_WALL_SECONDS",
    "MAX_HANDOFF_BYTES",
    "MIN_AVAILABLE_MEMORY_BYTES",
    "PARENT_MAX_CURRENT_RSS_BYTES",
    "PRO20RSRC1ForensicUncertainty",
    "PRO20RSRC1IsolationError",
    "REQUEST_SCHEMA",
    "RESOURCE_POLICY",
    "RESPONSE_SCHEMA",
    "check_parent_resource_premises",
    "child_process_rss_bytes",
    "classify_child_response",
    "current_child_peak_rss_bytes",
    "decode_attempt_request",
    "decode_bundle",
    "encode_attempt_request",
    "encode_bundle",
    "encode_child_response",
    "execute_child_request",
    "launch_child_process",
    "plan_status",
    "parse_canonical_handoff",
    "reduce_completed_child",
]
