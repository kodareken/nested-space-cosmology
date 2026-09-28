"""Closed first-event HLT17/PROTO19 checkpoint and journal schema.

This module is the wire contract for a synthetic first-event store.  It does
not write a filesystem, authenticate a physical origin or source, authorize a
campaign, or claim calibration, eligibility, or physics.  A checkpoint names
the journal-record hash that produced it; that record does not name a
successor checkpoint hash.  Kernel-transition digests, full-cursor hashes,
journal hashes, and checkpoint hashes are distinct identities.  Store
generation and journal sequence are not kernel counters.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from hashlib import sha256
import re
from types import MappingProxyType
from typing import Final

from recursive_horizons.evidence_io import CanonicalJSONError, canonical_json_bytes, load_canonical_json

from .hlt17_imp1_cursor import (
    CFL_OWNER,
    FRESH_READY,
    FUTURE_AUTHENTICATED_ORIGIN_SEAM,
    FUTURE_SOURCE_MANIFEST_SEAM,
    HLT17_GENESIS_DIGEST,
    HLT17_MEMBER_KEYS,
    HLT17_PRODUCTION_MEMBERS,
    HLT17_SYNTHETIC_MEMBER_KEYS,
    HLT17IMP1Cursor,
    IDENTITY_SCOPE_LIVE,
    IDENTITY_SCOPE_SYNTHETIC,
    RETRY_PENDING,
    SOURCE_OWNER,
    TEMPORAL_OWNER,
    cursor_sha256,
    hlt17_member_spec,
    ledger_sha256,
    restore_hlt17_cursor,
)
from .hlt17_admission_runtime import (
    IMPLEMENTATION_IDENTITY_KEYS,
    hlt17_implementation_identity,
)
from .hlt17_member_checkpoint import HLT17CheckpointError, HLT17GenerationBundle
from .hlt17_member_codec import ORIGIN_PROTOCOL, RUNTIME_PROTOCOL, SNAPSHOT_ROLE_ACCEPTED
from .numerical_engine import COMPARATOR_METHOD, PRIMARY_METHOD
from .tdg6_temporal_admission_design import (
    TDG6_MAXIMUM_TEMPORAL_RETRIES_PER_STEP,
    TDG6_MINIMUM_MACRO_STEP,
)
from .tdg11_imp1_ledger import (
    IMP1_ACCEPTED_SUBSTEP_COUNT,
    IMP1_METHOD_STAGE_RECORD_COUNT,
    append_imp1_rejection,
    restore_imp1_checkpoint_extension,
)


PROTOCOL_ARTIFACT_ID: Final[str] = "FGC-2-SF1-PROTO19"
SCHEMA_VERSION: Final[int] = 1
ROOT_SCHEMA: Final[str] = "FGC-1-HLT17-campaign-root-v1"
JOURNAL_SCHEMA: Final[str] = "FGC-1-HLT17-campaign-journal-v1"
CHECKPOINT_SCHEMA: Final[str] = "FGC-1-HLT17-campaign-checkpoint-v1"
GENERATION_MANIFEST_SCHEMA: Final[str] = "FGC-1-HLT17-generation-manifest-v1"
SESSION_OPEN_SCHEMA: Final[str] = "FGC-1-HLT17-writer-session-open-v1"
SESSION_CLOSE_SCHEMA: Final[str] = "FGC-1-HLT17-writer-session-close-v1"
TERMINAL_LOCK_SCHEMA: Final[str] = "FGC-1-HLT17-terminal-lock-v1"
TEMPORAL_EXHAUSTION_REASONS: Final[frozenset[str]] = frozenset(
    {"coordinate_lattice", "maximum_temporal_retries", "minimum_macro_step"}
)
STORE_KIND_SYNTHETIC: Final[str] = "synthetic_temporary"
GENESIS_DIGEST: Final[str] = HLT17_GENESIS_DIGEST
EVENT_ORIGIN: Final[float] = 23.0 / 16.0
EVENT_TARGET: Final[float] = 3.0 / 2.0
EVENT_ORIGIN_HEX: Final[str] = EVENT_ORIGIN.hex()
EVENT_TARGET_HEX: Final[str] = EVENT_TARGET.hex()
COMMON_EVENT_INDEX: Final[int] = 23
ACCEPTED_SUBSTEPS: Final[int] = IMP1_ACCEPTED_SUBSTEP_COUNT
RK4_STAGE_RECORDS: Final[int] = IMP1_METHOD_STAGE_RECORD_COUNT[PRIMARY_METHOD]
SSPRK3_STAGE_RECORDS: Final[int] = IMP1_METHOD_STAGE_RECORD_COUNT[COMPARATOR_METHOD]
_SHA256 = re.compile(r"\A[0-9a-f]{64}\Z")

RECORD_KINDS: Final[tuple[str, ...]] = (
    "seed",
    "accepted_fine",
    "source_retry",
    "cfl_retry",
    "temporal_retry",
    "source_retry_exhausted",
    "cfl_retry_exhausted",
    "temporal_retry_exhausted",
    "invalid_premise",
    "resource_exhausted",
)
TERMINAL_KINDS: Final[frozenset[str]] = frozenset(
    {
        "source_retry_exhausted",
        "cfl_retry_exhausted",
        "temporal_retry_exhausted",
        "invalid_premise",
        "resource_exhausted",
    }
)
BLOB_UNCHANGED_KINDS: Final[frozenset[str]] = frozenset(
    {"temporal_retry_exhausted", "invalid_premise", "resource_exhausted"}
)
CURSOR_CHANGED_RETRY_KINDS: Final[frozenset[str]] = frozenset(
    {
        "source_retry",
        "cfl_retry",
        "temporal_retry",
        "source_retry_exhausted",
        "cfl_retry_exhausted",
    }
)
DISPOSITIONS: Final[frozenset[str]] = frozenset(
    {"nonterminal", "first_event_complete"} | TERMINAL_KINDS
)

ROOT_MANIFEST_KEYS: Final[tuple[str, ...]] = (
    "schema",
    "schema_version",
    "protocol_artifact_id",
    "campaign_id",
    "store_kind",
    "member_keys",
    "event_origin_time_hex",
    "event_target_time_hex",
    "common_event_index",
    "historical_origin_authenticated",
    "source_authenticated",
    "environment_authenticated",
    "campaign_execution_authorized",
    "production_write_authorized",
    "live_runner_authorized",
    "calibration_claimed",
    "eligibility_claimed",
    "physical_classification",
    "physical_origin_seam",
    "source_manifest_seam",
    "environment_gate_seam",
    "resource_gate_seam",
    "independent_binder_seam",
    "live_runner_seam",
)
JOURNAL_ENVELOPE_KEYS: Final[tuple[str, ...]] = (
    "schema",
    "schema_version",
    "protocol_artifact_id",
    "campaign_id",
    "kind",
    "journal_sequence",
    "journal_parent_sha256",
    "store_generation",
    "predecessor_checkpoint_sha256",
    "payload",
    "campaign_execution_authorized",
    "historical_origin_authenticated",
    "production_write_authorized",
    "physical_classification",
    "calibration_claimed",
    "eligibility_claimed",
)
FORBIDDEN_JOURNAL_KEYS: Final[frozenset[str]] = frozenset(
    {
        "checkpoint_sha256",
        "successor_checkpoint_sha256",
        "generation_manifest_sha256",
        "latest",
        "journal_tip_sha256",
        "checkpoint_tip_sha256",
    }
)
SEED_PAYLOAD_KEYS: Final[tuple[str, ...]] = ("member_keys", "members")
ATTEMPT_PAYLOAD_KEYS: Final[tuple[str, ...]] = (
    "member_key",
    "predecessor_descriptor_sha256",
    "predecessor_payload_sha256",
    "predecessor_cursor_sha256",
    "successor_descriptor_sha256",
    "successor_payload_sha256",
    "successor_cursor_sha256",
    "kernel_transition_parent_sha256",
    "kernel_transition_digest",
    "accepted_generation",
    "cursor_generation",
    "attempt_serial",
    "kernel_checkpoint_generation",
    "source_current",
    "source_total",
    "cfl_current",
    "cfl_total",
    "mode",
    "latest_owner",
    "accepted_time_hex",
    "physical_state_sha256",
    "accepted_state_preserved",
    "tdg6_snapshot_sha256",
    "parked",
    "public_fresh",
    "executable_retry_cursor",
    "step_index",
    "transaction_serial",
    "terminal_evidence",
)
MEMBER_STATE_KEYS: Final[tuple[str, ...]] = (
    "descriptor_sha256",
    "payload_sha256",
    "cursor_sha256",
    "kernel_transition_digest",
    "kernel_transition_parent_sha256",
    "accepted_generation",
    "cursor_generation",
    "attempt_serial",
    "kernel_checkpoint_generation",
    "accepted_time_hex",
    "physical_state_sha256",
    "mode",
    "latest_owner",
    "source_current",
    "source_total",
    "cfl_current",
    "cfl_total",
    "tdg6_snapshot_sha256",
    "parked",
    "public_fresh",
    "step_index",
    "transaction_serial",
    *IMPLEMENTATION_IDENTITY_KEYS,
)
CHECKPOINT_KEYS: Final[tuple[str, ...]] = (
    "schema",
    "schema_version",
    "protocol_artifact_id",
    "campaign_id",
    "store_generation",
    "journal_sequence",
    "journal_record_sha256",
    "journal_parent_sha256",
    "predecessor_checkpoint_sha256",
    "event_origin_time_hex",
    "event_target_time_hex",
    "member_keys",
    "members",
    "kind",
    "changed_member_key",
    "disposition",
    "first_event_complete",
    "terminal",
    "campaign_execution_authorized",
    "historical_origin_authenticated",
    "production_write_authorized",
    "physical_classification",
    "calibration_claimed",
    "eligibility_claimed",
)
GENERATION_MANIFEST_KEYS: Final[tuple[str, ...]] = (
    "schema",
    "schema_version",
    "protocol_artifact_id",
    "campaign_id",
    "store_generation",
    "journal_sequence",
    "journal_record_sha256",
    "checkpoint_sha256",
    "predecessor_checkpoint_sha256",
    "kind",
    "changed_member_key",
    "new_blob_sha256s",
    "leaves",
    "production_write_authorized",
    "campaign_execution_authorized",
)
SESSION_OPEN_KEYS: Final[tuple[str, ...]] = (
    "schema",
    "schema_version",
    "protocol_artifact_id",
    "campaign_id",
    "session_index",
    "owner_token",
    "host",
    "pid",
    "thread_id",
    "exclusion_lock_sha256",
    "production_write_authorized",
    "campaign_execution_authorized",
)
SESSION_CLOSE_KEYS: Final[tuple[str, ...]] = (
    "schema",
    "schema_version",
    "protocol_artifact_id",
    "campaign_id",
    "session_index",
    "session_open_sha256",
    "owner_token",
)
TERMINAL_LOCK_KEYS: Final[tuple[str, ...]] = (
    "schema",
    "schema_version",
    "protocol_artifact_id",
    "campaign_id",
    "checkpoint_sha256",
    "journal_record_sha256",
    "store_generation",
    "disposition",
    "changed_member_key",
)
PRODUCTION_MEMBER_IDENTITIES: Final[tuple[tuple[str, int, int, int, int], ...]] = tuple(
    (
        item.member_key,
        item.point_count,
        item.owned_rows,
        item.fallback_d01,
        item.fallback_d12,
    )
    for item in HLT17_PRODUCTION_MEMBERS
)

FUTURE_PHYSICAL_ORIGIN_SEAM: Final[str] = FUTURE_AUTHENTICATED_ORIGIN_SEAM
FUTURE_ENVIRONMENT_GATE_SEAM: Final[str] = (
    "environment and numerical-platform qualification is an external PRO20 "
    "gate; this schema never authenticates a live machine or compiler image"
)
FUTURE_RESOURCE_GATE_SEAM: Final[str] = (
    "largest-grid and physical-RHS resource qualification is an external "
    "gate; this schema does not execute or authorize a production preparation"
)
FUTURE_INDEPENDENT_BINDER_SEAM: Final[str] = (
    "an independent PRO20 binder must reconstruct evidence without trusting "
    "runner labels; this schema is not that binder"
)
FUTURE_LIVE_RUNNER_SEAM: Final[str] = (
    "this slice provides no live runner, Make target, or production write "
    "authority"
)
FUTURE_SOURCE_MANIFEST_SEAM_TEXT: Final[str] = FUTURE_SOURCE_MANIFEST_SEAM

_ROOT_KEY_SET = frozenset(ROOT_MANIFEST_KEYS)
_JOURNAL_KEY_SET = frozenset(JOURNAL_ENVELOPE_KEYS)
_SEED_PAYLOAD_KEY_SET = frozenset(SEED_PAYLOAD_KEYS)
_ATTEMPT_PAYLOAD_KEY_SET = frozenset(ATTEMPT_PAYLOAD_KEYS)
_MEMBER_STATE_KEY_SET = frozenset(MEMBER_STATE_KEYS)
_CHECKPOINT_KEY_SET = frozenset(CHECKPOINT_KEYS)
_MANIFEST_KEY_SET = frozenset(GENERATION_MANIFEST_KEYS)
_SESSION_OPEN_KEY_SET = frozenset(SESSION_OPEN_KEYS)
_SESSION_CLOSE_KEY_SET = frozenset(SESSION_CLOSE_KEYS)
_TERMINAL_LOCK_KEY_SET = frozenset(TERMINAL_LOCK_KEYS)


class HLT17ProtocolError(ValueError):
    """Malformed PROTO19 checkpoint, journal, manifest, or transition data."""


def _fail(message: str) -> None:
    raise HLT17ProtocolError(message)


def _digest(payload: bytes) -> str:
    return sha256(payload).hexdigest()


def _text(value: object, label: str) -> str:
    if type(value) is not str or not value:
        _fail(f"{label} must be nonempty text")
    return value


def _digest_text(value: object, label: str) -> str:
    text = _text(value, label)
    if _SHA256.match(text) is None:
        _fail(f"{label} is not a lowercase SHA-256 digest")
    return text


def _bool_false(value: object, label: str) -> bool:
    if type(value) is not bool or value is not False:
        _fail(f"{label} must be false")
    return False


def _bool(value: object, label: str) -> bool:
    if type(value) is not bool:
        _fail(f"{label} must be a built-in bool")
    return value


def _nonneg_int(value: object, label: str) -> int:
    if type(value) is not int or value < 0:
        _fail(f"{label} must be a nonnegative built-in integer")
    return value


def _hex_time(value: object, label: str) -> str:
    text = _text(value, label)
    try:
        number = float.fromhex(text)
    except ValueError as error:
        raise HLT17ProtocolError(f"{label} is not a binary64 hex time") from error
    if not number == number or number.hex() != text:
        _fail(f"{label} is not a canonical binary64 hex time")
    return text


def _exact_mapping(value: object, names: frozenset[str], label: str) -> dict[str, object]:
    if type(value) is not dict:
        _fail(f"{label} must be a JSON object")
    if set(value) != names:
        _fail(f"{label} keys differ")
    return dict(value)


def _canonical(value: object, label: str) -> bytes:
    try:
        raw = canonical_json_bytes(value)
        loaded = load_canonical_json(raw)
    except (CanonicalJSONError, TypeError, ValueError) as error:
        raise HLT17ProtocolError(f"{label} is not canonical JSON") from error
    if loaded != value:
        _fail(f"{label} is not a redundant canonical mapping")
    return raw


def _same_bits(left: float, right: float) -> bool:
    return float(left).hex() == float(right).hex()


def _optional_owner(value: object) -> str | None:
    if value is None:
        return None
    text = _text(value, "latest_owner")
    if text not in {SOURCE_OWNER, CFL_OWNER, TEMPORAL_OWNER}:
        _fail("latest_owner differs")
    return text


def _optional_member_key(value: object, allowed: tuple[str, ...]) -> str | None:
    if value is None:
        return None
    text = _text(value, "member_key")
    if text not in allowed:
        _fail("member_key is outside the cohort")
    return text


def production_member_identities() -> tuple[tuple[str, int, int, int, int], ...]:
    """Exact six production metadata identities. Not a size or run claim."""

    expected = (
        ("RK4-2049", 2049, 2044, 16352, 32704),
        ("RK4-4097", 4097, 4092, 32736, 65472),
        ("RK4-8193", 8193, 8188, 65504, 131008),
        ("SSPRK3-4097", 4097, 4092, 32736, 65472),
        ("SSPRK3-8193", 8193, 8188, 65504, 131008),
        ("SSPRK3-16385", 16385, 16380, 131040, 262080),
    )
    if PRODUCTION_MEMBER_IDENTITIES != expected:
        _fail("HLT17 production metadata identities differ")
    if RK4_STAGE_RECORDS != 20 or SSPRK3_STAGE_RECORDS != 16:
        _fail("accepted fine stage-record counts differ")
    if ACCEPTED_SUBSTEPS != 4:
        _fail("accepted fine substep count differs")
    return expected


def cohort_member_keys(value: object) -> tuple[str, ...]:
    if type(value) is not list or not value:
        _fail("member_keys must be a nonempty JSON array")
    keys = tuple(_text(item, "member_key") for item in value)
    if keys == HLT17_SYNTHETIC_MEMBER_KEYS:
        for key in keys:
            if hlt17_member_spec(key).identity_scope != IDENTITY_SCOPE_SYNTHETIC:
                _fail("synthetic cohort used a live identity")
        return keys
    if keys == HLT17_MEMBER_KEYS:
        production_member_identities()
        for key in keys:
            if hlt17_member_spec(key).identity_scope != IDENTITY_SCOPE_LIVE:
                _fail("production cohort used a synthetic identity")
        return keys
    _fail("cohort must be the synthetic RK4-9/SSPRK3-9 pair or the six production keys")
    raise AssertionError("unreachable")


def generation_directory_name(store_generation: int) -> str:
    index = _nonneg_int(store_generation, "store_generation")
    return f"generation-{index:016d}"


def parse_generation_directory_name(name: object) -> int:
    text = _text(name, "generation directory")
    prefix = "generation-"
    if not text.startswith(prefix) or len(text) != len(prefix) + 16:
        _fail("generation directory name differs")
    digits = text[len(prefix) :]
    if not digits.isdigit() or digits != f"{int(digits):016d}":
        _fail("generation directory name differs")
    return int(digits)


def session_open_name(session_index: int) -> str:
    return f"session-{_nonneg_int(session_index, 'session_index'):016d}-open.json"


def session_close_name(session_index: int) -> str:
    return f"session-{_nonneg_int(session_index, 'session_index'):016d}-close.json"


def parse_session_name(name: object) -> tuple[int, str]:
    text = _text(name, "session file")
    for suffix, kind in (("-open.json", "open"), ("-close.json", "close")):
        prefix = "session-"
        if text.startswith(prefix) and text.endswith(suffix):
            digits = text[len(prefix) : -len(suffix)]
            if len(digits) == 16 and digits.isdigit() and digits == f"{int(digits):016d}":
                return int(digits), kind
    _fail("session file name differs")
    raise AssertionError("unreachable")


def revalidate_bundle(bundle: object) -> HLT17GenerationBundle:
    """Rebuild a member bundle from its three immutable byte-string components."""

    if not isinstance(bundle, HLT17GenerationBundle):
        raise TypeError("bundle must be HLT17GenerationBundle")
    try:
        rebuilt = HLT17GenerationBundle(
            bytes(bundle.descriptor),
            bytes(bundle.payload),
            bytes(bundle.cursor_bytes),
        )
    except (HLT17CheckpointError, TypeError, ValueError) as error:
        raise HLT17ProtocolError(
            "member bundle failed immutable-byte revalidation"
        ) from error
    if (
        rebuilt.descriptor != bytes(bundle.descriptor)
        or rebuilt.payload != bytes(bundle.payload)
        or rebuilt.cursor_bytes != bytes(bundle.cursor_bytes)
    ):
        _fail("revalidated bundle bytes differ from the supplied source")
    return rebuilt


def bundle_from_bytes(
    descriptor: bytes, payload: bytes, cursor_bytes: bytes
) -> HLT17GenerationBundle:
    if type(descriptor) is not bytes or type(payload) is not bytes or type(cursor_bytes) is not bytes:
        raise TypeError("bundle parts must be bytes")
    try:
        return HLT17GenerationBundle(descriptor, payload, cursor_bytes)
    except (HLT17CheckpointError, TypeError, ValueError) as error:
        raise HLT17ProtocolError(
            "member bundle failed immutable-byte revalidation"
        ) from error


def _cursor_from_bundle(bundle: HLT17GenerationBundle) -> HLT17IMP1Cursor:
    try:
        mapping = load_canonical_json(bundle.cursor_bytes)
        cursor = restore_hlt17_cursor(mapping)
    except (CanonicalJSONError, TypeError, ValueError) as error:
        raise HLT17ProtocolError("full-cursor bytes are not a closed HLT17 cursor") from error
    if _digest(bundle.cursor_bytes) != cursor_sha256(cursor):
        _fail("full-cursor hash is not the hash of the immutable cursor bytes")
    if cursor.descriptor_sha256 != _digest(bundle.descriptor):
        _fail("cursor does not name the rehashed descriptor bytes")
    if cursor.descriptor_sha256 != bundle.descriptor_sha256:
        _fail("bundle descriptor hash disagrees with the descriptor bytes")
    return cursor


def _descriptor_mapping(bundle: HLT17GenerationBundle) -> dict[str, object]:
    try:
        mapping = load_canonical_json(bundle.descriptor)
    except (CanonicalJSONError, TypeError, ValueError) as error:
        raise HLT17ProtocolError("accepted descriptor is not canonical JSON") from error
    if type(mapping) is not dict:
        _fail("accepted descriptor is not a JSON object")
    if mapping.get("snapshot_role") != SNAPSHOT_ROLE_ACCEPTED:
        _fail("published bundle must export an accepted physical-state descriptor")
    return mapping


def member_state_from_bundle(bundle: object) -> dict[str, object]:
    """Derive member checkpoint fields from immutable bundle bytes only."""

    bundle = revalidate_bundle(bundle)
    cursor = _cursor_from_bundle(bundle)
    descriptor = _descriptor_mapping(bundle)
    identity = descriptor["runtime_identity"]
    if type(identity) is not dict:
        _fail("runtime identity differs")
    accepted_generation = _nonneg_int(identity["accepted_generation"], "accepted_generation")
    if accepted_generation != bundle.accepted_generation:
        _fail("accepted generation label disagrees with the descriptor bytes")
    if cursor.checkpoint_generation != accepted_generation:
        _fail("kernel checkpoint generation differs from accepted physical generation")
    if cursor.ledger.accepted_macro_step_count != accepted_generation:
        _fail("IMP1 accepted macro-step count differs from accepted generation")
    physical = _digest_text(descriptor["physical_state_sha256"], "physical_state_sha256")
    if physical != bundle.physical_state_sha256 or physical != cursor.physical_state_sha256:
        _fail("physical u/p/q hash disagrees with immutable bytes")
    if physical == bundle.descriptor_sha256:
        _fail("descriptor and physical identities are not interchangeable")
    parked = _same_bits(cursor.accepted_time, cursor.event_target)
    mapping = {
        "descriptor_sha256": _digest(bundle.descriptor),
        "payload_sha256": _digest(bundle.payload),
        "cursor_sha256": _digest(bundle.cursor_bytes),
        "kernel_transition_digest": cursor.kernel_transition_digest,
        "kernel_transition_parent_sha256": cursor.kernel_transition_parent_sha256,
        "accepted_generation": accepted_generation,
        "cursor_generation": cursor.cursor_generation,
        "attempt_serial": cursor.attempt_serial,
        "kernel_checkpoint_generation": cursor.checkpoint_generation,
        "accepted_time_hex": cursor.accepted_time.hex(),
        "physical_state_sha256": physical,
        "mode": cursor.mode,
        "latest_owner": cursor.latest_owner,
        "source_current": cursor.source_current,
        "source_total": cursor.source_total,
        "cfl_current": cursor.cfl_current,
        "cfl_total": cursor.cfl_total,
        "tdg6_snapshot_sha256": _digest_text(
            load_canonical_json(bundle.cursor_bytes)["inherited_snapshot_sha256"],
            "inherited_snapshot_sha256",
        ),
        "parked": parked,
        "public_fresh": cursor.public_fresh_permitted(),
        "step_index": cursor.accepted_step_index,
        "transaction_serial": cursor.accepted_transaction_serial,
        **hlt17_implementation_identity(),
    }
    _exact_mapping(mapping, _MEMBER_STATE_KEY_SET, "member state")
    _canonical(mapping, "member state")
    return mapping


def _require_distinct_identities(
    *,
    checkpoint_sha256: str,
    journal_sha256: str,
    states: Mapping[str, Mapping[str, object]],
) -> None:
    if checkpoint_sha256 == journal_sha256:
        _fail("checkpoint hash collided with the journal hash")
    for key, state in states.items():
        cursor_hash = _digest_text(state["cursor_sha256"], "cursor_sha256")
        kernel = _digest_text(state["kernel_transition_digest"], "kernel_transition_digest")
        identities = (checkpoint_sha256, journal_sha256, cursor_hash, kernel)
        if len(set(identities)) != 4:
            _fail(f"{key} checkpoint/journal/cursor/kernel identities are not distinct")
        if cursor_hash == _digest_text(state["descriptor_sha256"], "descriptor_sha256"):
            _fail(f"{key} full-cursor hash collided with the descriptor hash")
        if kernel == _digest_text(state["physical_state_sha256"], "physical_state_sha256"):
            _fail(f"{key} kernel-transition hash collided with the physical hash")


def _stage_records(member_key: str) -> int:
    spec = hlt17_member_spec(member_key)
    return IMP1_METHOD_STAGE_RECORD_COUNT[spec.method]


def _cohort_complete(states: Mapping[str, Mapping[str, object]]) -> bool:
    for state in states.values():
        if state["mode"] != FRESH_READY:
            return False
        if state["latest_owner"] is not None:
            return False
        if state["parked"] is not True or state["public_fresh"] is not True:
            return False
        if state["accepted_time_hex"] != EVENT_TARGET_HEX:
            return False
        if state["source_current"] != 0 or state["cfl_current"] != 0:
            return False
    return True


def validate_root_manifest(value: object) -> dict[str, object]:
    item = _exact_mapping(value, _ROOT_KEY_SET, "root manifest")
    if item["schema"] != ROOT_SCHEMA or item["schema_version"] != SCHEMA_VERSION:
        _fail("root schema differs")
    if item["protocol_artifact_id"] != PROTOCOL_ARTIFACT_ID:
        _fail("root protocol differs")
    campaign_id = _text(item["campaign_id"], "campaign_id")
    if "/" in campaign_id or "\\" in campaign_id or ".." in campaign_id:
        _fail("campaign_id is not a safe identifier")
    if item["store_kind"] != STORE_KIND_SYNTHETIC:
        _fail("store_kind must be synthetic_temporary")
    keys = cohort_member_keys(item["member_keys"])
    if _hex_time(item["event_origin_time_hex"], "event_origin_time_hex") != EVENT_ORIGIN_HEX:
        _fail("first-event origin must be t=23/16")
    if _hex_time(item["event_target_time_hex"], "event_target_time_hex") != EVENT_TARGET_HEX:
        _fail("first-event target must be t=3/2")
    if _nonneg_int(item["common_event_index"], "common_event_index") != COMMON_EVENT_INDEX:
        _fail("first-event index must be 23")
    for name in (
        "historical_origin_authenticated",
        "source_authenticated",
        "environment_authenticated",
        "campaign_execution_authorized",
        "production_write_authorized",
        "live_runner_authorized",
        "calibration_claimed",
        "eligibility_claimed",
        "physical_classification",
    ):
        _bool_false(item[name], name)
    seams = {
        "physical_origin_seam": FUTURE_PHYSICAL_ORIGIN_SEAM,
        "source_manifest_seam": FUTURE_SOURCE_MANIFEST_SEAM_TEXT,
        "environment_gate_seam": FUTURE_ENVIRONMENT_GATE_SEAM,
        "resource_gate_seam": FUTURE_RESOURCE_GATE_SEAM,
        "independent_binder_seam": FUTURE_INDEPENDENT_BINDER_SEAM,
        "live_runner_seam": FUTURE_LIVE_RUNNER_SEAM,
    }
    for name, expected in seams.items():
        if item[name] != expected:
            _fail(f"{name} differs")
    item["member_keys"] = list(keys)
    _canonical(item, "root manifest")
    return item


def build_root_manifest(*, campaign_id: str, member_keys: object) -> dict[str, object]:
    mapping = {
        "schema": ROOT_SCHEMA,
        "schema_version": SCHEMA_VERSION,
        "protocol_artifact_id": PROTOCOL_ARTIFACT_ID,
        "campaign_id": campaign_id,
        "store_kind": STORE_KIND_SYNTHETIC,
        "member_keys": list(member_keys),
        "event_origin_time_hex": EVENT_ORIGIN_HEX,
        "event_target_time_hex": EVENT_TARGET_HEX,
        "common_event_index": COMMON_EVENT_INDEX,
        "historical_origin_authenticated": False,
        "source_authenticated": False,
        "environment_authenticated": False,
        "campaign_execution_authorized": False,
        "production_write_authorized": False,
        "live_runner_authorized": False,
        "calibration_claimed": False,
        "eligibility_claimed": False,
        "physical_classification": False,
        "physical_origin_seam": FUTURE_PHYSICAL_ORIGIN_SEAM,
        "source_manifest_seam": FUTURE_SOURCE_MANIFEST_SEAM_TEXT,
        "environment_gate_seam": FUTURE_ENVIRONMENT_GATE_SEAM,
        "resource_gate_seam": FUTURE_RESOURCE_GATE_SEAM,
        "independent_binder_seam": FUTURE_INDEPENDENT_BINDER_SEAM,
        "live_runner_seam": FUTURE_LIVE_RUNNER_SEAM,
    }
    return validate_root_manifest(mapping)


def encode_session_open(
    *,
    campaign_id: str,
    session_index: int,
    owner_token: str,
    host: str,
    pid: int,
    thread_id: int,
    exclusion_lock_sha256: str,
) -> dict[str, object]:
    mapping = {
        "schema": SESSION_OPEN_SCHEMA,
        "schema_version": SCHEMA_VERSION,
        "protocol_artifact_id": PROTOCOL_ARTIFACT_ID,
        "campaign_id": _text(campaign_id, "campaign_id"),
        "session_index": _nonneg_int(session_index, "session_index"),
        "owner_token": _text(owner_token, "owner_token"),
        "host": _text(host, "host"),
        "pid": _nonneg_int(pid, "pid"),
        "thread_id": _nonneg_int(thread_id, "thread_id"),
        "exclusion_lock_sha256": _digest_text(
            exclusion_lock_sha256, "exclusion_lock_sha256"
        ),
        "production_write_authorized": False,
        "campaign_execution_authorized": False,
    }
    _exact_mapping(mapping, _SESSION_OPEN_KEY_SET, "session open")
    _canonical(mapping, "session open")
    return mapping


def encode_session_close(
    *,
    campaign_id: str,
    session_index: int,
    session_open_sha256: str,
    owner_token: str,
) -> dict[str, object]:
    mapping = {
        "schema": SESSION_CLOSE_SCHEMA,
        "schema_version": SCHEMA_VERSION,
        "protocol_artifact_id": PROTOCOL_ARTIFACT_ID,
        "campaign_id": _text(campaign_id, "campaign_id"),
        "session_index": _nonneg_int(session_index, "session_index"),
        "session_open_sha256": _digest_text(session_open_sha256, "session_open_sha256"),
        "owner_token": _text(owner_token, "owner_token"),
    }
    _exact_mapping(mapping, _SESSION_CLOSE_KEY_SET, "session close")
    _canonical(mapping, "session close")
    return mapping


def encode_terminal_lock(
    *,
    campaign_id: str,
    checkpoint_sha256: str,
    journal_record_sha256: str,
    store_generation: int,
    disposition: str,
    changed_member_key: str | None,
) -> dict[str, object]:
    if disposition not in DISPOSITIONS or disposition == "nonterminal":
        _fail("terminal lock disposition differs")
    mapping = {
        "schema": TERMINAL_LOCK_SCHEMA,
        "schema_version": SCHEMA_VERSION,
        "protocol_artifact_id": PROTOCOL_ARTIFACT_ID,
        "campaign_id": _text(campaign_id, "campaign_id"),
        "checkpoint_sha256": _digest_text(checkpoint_sha256, "checkpoint_sha256"),
        "journal_record_sha256": _digest_text(
            journal_record_sha256, "journal_record_sha256"
        ),
        "store_generation": _nonneg_int(store_generation, "store_generation"),
        "disposition": _text(disposition, "disposition"),
        "changed_member_key": changed_member_key,
    }
    _exact_mapping(mapping, _TERMINAL_LOCK_KEY_SET, "terminal lock")
    _canonical(mapping, "terminal lock")
    return mapping


def validate_session_open(value: object) -> dict[str, object]:
    item = _exact_mapping(value, _SESSION_OPEN_KEY_SET, "session open")
    if item["schema"] != SESSION_OPEN_SCHEMA or item["schema_version"] != SCHEMA_VERSION:
        _fail("session open schema differs")
    encode_session_open(
        campaign_id=item["campaign_id"],
        session_index=item["session_index"],
        owner_token=item["owner_token"],
        host=item["host"],
        pid=item["pid"],
        thread_id=item["thread_id"],
        exclusion_lock_sha256=item["exclusion_lock_sha256"],
    )
    return item


def validate_session_close(value: object) -> dict[str, object]:
    item = _exact_mapping(value, _SESSION_CLOSE_KEY_SET, "session close")
    if item["schema"] != SESSION_CLOSE_SCHEMA or item["schema_version"] != SCHEMA_VERSION:
        _fail("session close schema differs")
    encode_session_close(
        campaign_id=item["campaign_id"],
        session_index=item["session_index"],
        session_open_sha256=item["session_open_sha256"],
        owner_token=item["owner_token"],
    )
    return item


def validate_terminal_lock(value: object) -> dict[str, object]:
    item = _exact_mapping(value, _TERMINAL_LOCK_KEY_SET, "terminal lock")
    if item["schema"] != TERMINAL_LOCK_SCHEMA or item["schema_version"] != SCHEMA_VERSION:
        _fail("terminal lock schema differs")
    encode_terminal_lock(
        campaign_id=item["campaign_id"],
        checkpoint_sha256=item["checkpoint_sha256"],
        journal_record_sha256=item["journal_record_sha256"],
        store_generation=item["store_generation"],
        disposition=item["disposition"],
        changed_member_key=item["changed_member_key"],
    )
    return item


@dataclass(frozen=True, slots=True)
class HLT17ValidatedGeneration:
    """Closed generation: journal, checkpoint, manifest, and referenced blobs."""

    kind: str
    store_generation: int
    journal_sequence: int
    journal_raw: bytes
    journal_sha256: str
    checkpoint_raw: bytes
    checkpoint_sha256: str
    manifest_raw: bytes
    files: Mapping[str, bytes]
    bundles: Mapping[str, HLT17GenerationBundle]
    changed_member_key: str | None
    first_event_complete: bool
    terminal: bool
    disposition: str

    def __post_init__(self) -> None:
        object.__setattr__(self, "files", MappingProxyType(dict(self.files)))
        object.__setattr__(
            self,
            "bundles",
            MappingProxyType(
                {key: revalidate_bundle(bundle) for key, bundle in dict(self.bundles).items()}
            ),
        )
        if type(self.journal_raw) is not bytes or type(self.checkpoint_raw) is not bytes:
            raise TypeError("generation documents must be bytes")
        if type(self.manifest_raw) is not bytes:
            raise TypeError("generation manifest must be bytes")
        derived = _derived_generation_fields(
            journal_raw=self.journal_raw,
            checkpoint_raw=self.checkpoint_raw,
            manifest_raw=self.manifest_raw,
            files=self.files,
            bundles=self.bundles,
        )
        if _digest(self.journal_raw) != derived["journal_sha256"]:
            _fail("journal hash is not the hash of the journal bytes")
        if _digest(self.checkpoint_raw) != derived["checkpoint_sha256"]:
            _fail("checkpoint hash is not the hash of the checkpoint bytes")
        if derived["checkpoint_sha256"] == derived["journal_sha256"]:
            _fail("checkpoint hash collided with the journal hash")
        supplied = {
            "kind": self.kind,
            "store_generation": self.store_generation,
            "journal_sequence": self.journal_sequence,
            "journal_sha256": self.journal_sha256,
            "checkpoint_sha256": self.checkpoint_sha256,
            "changed_member_key": self.changed_member_key,
            "first_event_complete": self.first_event_complete,
            "terminal": self.terminal,
            "disposition": self.disposition,
        }
        for name, value in supplied.items():
            if value != derived[name]:
                _fail(f"{name} disagrees with the immutable generation documents")
        for name, value in derived.items():
            if name in supplied:
                object.__setattr__(self, name, value)

    @property
    def directory_name(self) -> str:
        return generation_directory_name(self.store_generation)

    @property
    def checkpoint_mapping(self) -> dict[str, object]:
        mapping = load_canonical_json(self.checkpoint_raw)
        if type(mapping) is not dict:
            _fail("checkpoint is not a JSON object")
        return mapping

    @property
    def journal_mapping(self) -> dict[str, object]:
        mapping = load_canonical_json(self.journal_raw)
        if type(mapping) is not dict:
            _fail("journal is not a JSON object")
        return mapping


def _authority_false(mapping: Mapping[str, object]) -> None:
    for name in (
        "campaign_execution_authorized",
        "historical_origin_authenticated",
        "production_write_authorized",
        "physical_classification",
        "calibration_claimed",
        "eligibility_claimed",
    ):
        if name in mapping:
            _bool_false(mapping[name], name)


def _validate_journal_envelope(
    value: object, *, campaign_id: str, member_keys: tuple[str, ...]
) -> dict[str, object]:
    if type(value) is not dict:
        _fail("journal must be a JSON object")
    if FORBIDDEN_JOURNAL_KEYS.intersection(value):
        _fail("journal contains a successor-checkpoint or latest pointer")
    item = _exact_mapping(value, _JOURNAL_KEY_SET, "journal envelope")
    if item["schema"] != JOURNAL_SCHEMA or item["schema_version"] != SCHEMA_VERSION:
        _fail("journal schema differs")
    if item["protocol_artifact_id"] != PROTOCOL_ARTIFACT_ID:
        _fail("journal protocol differs")
    if item["campaign_id"] != campaign_id:
        _fail("journal campaign_id differs")
    kind = _text(item["kind"], "kind")
    if kind not in RECORD_KINDS:
        _fail("journal kind differs")
    _nonneg_int(item["journal_sequence"], "journal_sequence")
    _digest_text(item["journal_parent_sha256"], "journal_parent_sha256")
    _nonneg_int(item["store_generation"], "store_generation")
    _digest_text(item["predecessor_checkpoint_sha256"], "predecessor_checkpoint_sha256")
    _authority_false(item)
    payload = item["payload"]
    if kind == "seed":
        body = _exact_mapping(payload, _SEED_PAYLOAD_KEY_SET, "seed payload")
        if tuple(cohort_member_keys(body["member_keys"])) != member_keys:
            _fail("seed payload member_keys differ")
        members = body["members"]
        if type(members) is not dict or set(members) != set(member_keys):
            _fail("seed payload members differ")
        for key in member_keys:
            state = _exact_mapping(members[key], _MEMBER_STATE_KEY_SET, "member state")
            members[key] = state
    else:
        body = _exact_mapping(payload, _ATTEMPT_PAYLOAD_KEY_SET, "attempt payload")
        if body["member_key"] not in member_keys:
            _fail("attempt member_key is outside the cohort")
        for name in (
            "predecessor_descriptor_sha256",
            "predecessor_payload_sha256",
            "predecessor_cursor_sha256",
            "successor_descriptor_sha256",
            "successor_payload_sha256",
            "successor_cursor_sha256",
            "kernel_transition_parent_sha256",
            "kernel_transition_digest",
            "physical_state_sha256",
            "tdg6_snapshot_sha256",
        ):
            _digest_text(body[name], name)
        if body["kernel_transition_digest"] == body["successor_cursor_sha256"]:
            _fail("kernel-transition hash is not the full-cursor hash")
        _optional_owner(body["latest_owner"])
        _hex_time(body["accepted_time_hex"], "accepted_time_hex")
        if body["mode"] not in {FRESH_READY, RETRY_PENDING}:
            _fail("attempt mode differs")
        _bool(body["accepted_state_preserved"], "accepted_state_preserved")
        _bool(body["parked"], "parked")
        _bool(body["public_fresh"], "public_fresh")
        _bool(body["executable_retry_cursor"], "executable_retry_cursor")
        for name in (
            "accepted_generation",
            "cursor_generation",
            "attempt_serial",
            "kernel_checkpoint_generation",
            "source_current",
            "source_total",
            "cfl_current",
            "cfl_total",
            "step_index",
            "transaction_serial",
        ):
            _nonneg_int(body[name], name)
        if body["terminal_evidence"] is not None and type(body["terminal_evidence"]) is not dict:
            _fail("terminal_evidence must be a JSON object or null")
    _canonical(item, "journal")
    return item


def _validate_checkpoint_document(
    value: object, *, campaign_id: str, member_keys: tuple[str, ...]
) -> dict[str, object]:
    if type(value) is not dict:
        _fail("checkpoint must be a JSON object")
    if "checkpoint_sha256" in value or "successor_checkpoint_sha256" in value:
        _fail("checkpoint contains a circular successor or self hash")
    if "latest" in value:
        _fail("checkpoint contains a mutable latest pointer")
    item = _exact_mapping(value, _CHECKPOINT_KEY_SET, "checkpoint")
    if item["schema"] != CHECKPOINT_SCHEMA or item["schema_version"] != SCHEMA_VERSION:
        _fail("checkpoint schema differs")
    if item["protocol_artifact_id"] != PROTOCOL_ARTIFACT_ID:
        _fail("checkpoint protocol differs")
    if item["campaign_id"] != campaign_id:
        _fail("checkpoint campaign_id differs")
    _nonneg_int(item["store_generation"], "store_generation")
    _nonneg_int(item["journal_sequence"], "journal_sequence")
    _digest_text(item["journal_record_sha256"], "journal_record_sha256")
    _digest_text(item["journal_parent_sha256"], "journal_parent_sha256")
    _digest_text(item["predecessor_checkpoint_sha256"], "predecessor_checkpoint_sha256")
    if _hex_time(item["event_origin_time_hex"], "event_origin_time_hex") != EVENT_ORIGIN_HEX:
        _fail("checkpoint origin differs")
    if _hex_time(item["event_target_time_hex"], "event_target_time_hex") != EVENT_TARGET_HEX:
        _fail("checkpoint target differs")
    if tuple(cohort_member_keys(item["member_keys"])) != member_keys:
        _fail("checkpoint member_keys differ")
    members = item["members"]
    if type(members) is not dict or set(members) != set(member_keys):
        _fail("checkpoint members differ")
    for key in member_keys:
        members[key] = _exact_mapping(members[key], _MEMBER_STATE_KEY_SET, "member state")
    kind = _text(item["kind"], "kind")
    if kind not in RECORD_KINDS:
        _fail("checkpoint kind differs")
    changed = _optional_member_key(item["changed_member_key"], member_keys)
    if kind == "seed":
        if changed is not None:
            _fail("seed cannot name a single changed member")
    elif changed is None:
        _fail("attempt checkpoint omitted changed_member_key")
    disposition = _text(item["disposition"], "disposition")
    if disposition not in DISPOSITIONS:
        _fail("checkpoint disposition differs")
    complete = _bool(item["first_event_complete"], "first_event_complete")
    terminal = _bool(item["terminal"], "terminal")
    if complete is True:
        if disposition != "first_event_complete" or terminal is not True:
            _fail("first-event completion must be a terminal first_event_complete disposition")
        if kind != "accepted_fine":
            _fail("first-event completion must be the last accepted fine")
    elif kind in TERMINAL_KINDS:
        if disposition != kind or terminal is not True or complete is not False:
            _fail("typed terminal disposition differs")
    else:
        if disposition != "nonterminal" or terminal is not False:
            _fail("nonterminal disposition differs")
    _authority_false(item)
    _canonical(item, "checkpoint")
    return item


def _derived_generation_fields(
    *,
    journal_raw: bytes,
    checkpoint_raw: bytes,
    manifest_raw: bytes,
    files: Mapping[str, bytes],
    bundles: Mapping[str, HLT17GenerationBundle],
) -> dict[str, object]:
    """Recompute generation identity from immutable documents, not constructor labels."""

    if files.get("journal.json") != journal_raw:
        _fail("journal.json file bytes differ from journal_raw")
    if files.get("checkpoint.json") != checkpoint_raw:
        _fail("checkpoint.json file bytes differ from checkpoint_raw")
    if files.get("generation_manifest.json") != manifest_raw:
        _fail("generation_manifest.json file bytes differ from manifest_raw")
    try:
        journal = load_canonical_json(journal_raw)
        checkpoint = load_canonical_json(checkpoint_raw)
        manifest = load_canonical_json(manifest_raw)
    except (CanonicalJSONError, TypeError, ValueError) as error:
        raise HLT17ProtocolError("generation documents are not canonical JSON") from error
    if type(journal) is not dict or type(checkpoint) is not dict or type(manifest) is not dict:
        _fail("generation documents are not JSON objects")
    campaign_id = _text(checkpoint.get("campaign_id"), "campaign_id")
    member_keys = tuple(cohort_member_keys(checkpoint.get("member_keys")))
    checkpoint = _validate_checkpoint_document(
        checkpoint, campaign_id=campaign_id, member_keys=member_keys
    )
    journal = _validate_journal_envelope(
        journal, campaign_id=campaign_id, member_keys=member_keys
    )
    manifest = _exact_mapping(manifest, _MANIFEST_KEY_SET, "generation manifest")
    if manifest["schema"] != GENERATION_MANIFEST_SCHEMA:
        _fail("generation manifest schema differs")
    journal_sha256 = _digest(journal_raw)
    checkpoint_sha256 = _digest(checkpoint_raw)
    if journal_sha256 != checkpoint["journal_record_sha256"]:
        _fail("checkpoint journal hash disagrees with journal bytes")
    if journal_sha256 != manifest["journal_record_sha256"]:
        _fail("manifest journal hash disagrees with journal bytes")
    if checkpoint_sha256 != manifest["checkpoint_sha256"]:
        _fail("manifest checkpoint hash disagrees with checkpoint bytes")
    if set(bundles) != set(member_keys):
        _fail("generation bundles do not match the cohort")
    derived_states = {key: member_state_from_bundle(bundles[key]) for key in member_keys}
    if derived_states != checkpoint["members"]:
        _fail("checkpoint member states disagree with immutable bundle bytes")
    complete = _cohort_complete(derived_states)
    if complete != bool(checkpoint["first_event_complete"]):
        _fail("first-event completion disagrees with the reconstructed cohort")
    return {
        "kind": str(checkpoint["kind"]),
        "store_generation": int(checkpoint["store_generation"]),
        "journal_sequence": int(checkpoint["journal_sequence"]),
        "journal_sha256": journal_sha256,
        "checkpoint_sha256": checkpoint_sha256,
        "changed_member_key": checkpoint["changed_member_key"],
        "first_event_complete": complete,
        "terminal": bool(checkpoint["terminal"]),
        "disposition": str(checkpoint["disposition"]),
    }


def _rebind_generation(generation: object) -> HLT17ValidatedGeneration:
    if not isinstance(generation, HLT17ValidatedGeneration):
        raise TypeError("generation must be HLT17ValidatedGeneration")
    generation.__post_init__()
    return generation


def _generation_forbids_append(generation: HLT17ValidatedGeneration) -> str | None:
    """Refuse successors from immutable predecessor checkpoint bytes, not labels."""

    try:
        checkpoint = load_canonical_json(generation.checkpoint_raw)
    except (CanonicalJSONError, TypeError, ValueError) as error:
        raise HLT17ProtocolError("predecessor checkpoint is not canonical JSON") from error
    if type(checkpoint) is not dict:
        _fail("predecessor checkpoint is not a JSON object")
    complete = bool(checkpoint.get("first_event_complete"))
    terminal = bool(checkpoint.get("terminal"))
    disposition = checkpoint.get("disposition")
    if complete or disposition == "first_event_complete":
        return "first-event completion is a final boundary"
    if terminal or disposition in TERMINAL_KINDS:
        return "typed terminal forbids further appends"
    return None


def _blob_path(digest: str) -> str:
    return f"blobs/{_digest_text(digest, 'blob digest')}"


def _files_from_parts(
    *,
    journal_raw: bytes,
    checkpoint_raw: bytes,
    manifest_mapping: Mapping[str, object],
    new_blobs: Mapping[str, bytes],
) -> dict[str, bytes]:
    files: dict[str, bytes] = {
        "journal.json": journal_raw,
        "checkpoint.json": checkpoint_raw,
    }
    for digest, raw in new_blobs.items():
        if type(raw) is not bytes:
            raise TypeError("blob bytes differ")
        if _digest(raw) != digest:
            _fail("content-addressed blob name disagrees with bytes")
        files[_blob_path(digest)] = raw
    leaves = {
        path: _digest(payload)
        for path, payload in files.items()
    }
    expected_leaves = manifest_mapping.get("leaves")
    if expected_leaves != leaves:
        _fail("generation manifest leaves disagree with the file inventory")
    manifest_raw = _canonical(dict(manifest_mapping), "generation manifest")
    files["generation_manifest.json"] = manifest_raw
    return files


def _new_blobs_for_cohort(
    *,
    bundles: Mapping[str, HLT17GenerationBundle],
    predecessor_bundles: Mapping[str, HLT17GenerationBundle] | None,
    changed_member_key: str | None,
) -> dict[str, bytes]:
    blobs: dict[str, bytes] = {}

    def add(raw: bytes) -> None:
        blobs[_digest(raw)] = raw

    if predecessor_bundles is None:
        for bundle in bundles.values():
            add(bundle.descriptor)
            add(bundle.payload)
            add(bundle.cursor_bytes)
        return blobs
    if changed_member_key is None:
        return blobs
    previous = predecessor_bundles[changed_member_key]
    current = bundles[changed_member_key]
    if current.descriptor != previous.descriptor:
        add(current.descriptor)
    if current.payload != previous.payload:
        add(current.payload)
    if current.cursor_bytes != previous.cursor_bytes:
        add(current.cursor_bytes)
    return blobs


def _require_seed_bundle(bundle: HLT17GenerationBundle, member_key: str) -> None:
    cursor = _cursor_from_bundle(bundle)
    descriptor = _descriptor_mapping(bundle)
    identity = descriptor["runtime_identity"]
    if identity["member_key"] != member_key or cursor.member_key != member_key:
        _fail("seed bundle member_key differs")
    if hlt17_member_spec(member_key).identity_scope != IDENTITY_SCOPE_SYNTHETIC:
        _fail("seed refuses production-size members")
    if cursor.event_target.hex() != EVENT_TARGET_HEX:
        _fail("seed event target must be t=3/2")
    if cursor.accepted_time.hex() != EVENT_ORIGIN_HEX:
        _fail("seed accepted time must be t=23/16")
    if bundle.accepted_generation != 0 or cursor.checkpoint_generation != 0:
        _fail("seed accepted generation must be zero")
    if cursor.cursor_generation != 0 or cursor.attempt_serial != 0:
        _fail("seed cursor/attempt revisions must be zero")
    if cursor.kernel_transition_parent_sha256 != GENESIS_DIGEST:
        _fail("seed kernel parent must be the genesis digest")
    predecessor = descriptor["predecessor_cursor_identity"]
    if type(predecessor) is not dict:
        _fail("seed predecessor identity differs")
    if predecessor.get("protocol_artifact_id") != ORIGIN_PROTOCOL:
        _fail("seed descriptor must name the sealed PROTO17 origin identity")
    if cursor.mode != FRESH_READY or cursor.latest_owner is not None:
        _fail("seed cursor must be fresh")
    if _same_bits(cursor.accepted_time, cursor.event_target):
        _fail("seed cannot already be parked at the event target")


def _require_fine_transition(
    previous: Mapping[str, object],
    current: Mapping[str, object],
    *,
    previous_bundle: HLT17GenerationBundle,
    current_bundle: HLT17GenerationBundle,
    member_key: str,
) -> None:
    if current["accepted_generation"] != previous["accepted_generation"] + 1:
        _fail("accepted fine must increment accepted generation")
    if current["kernel_checkpoint_generation"] != previous["kernel_checkpoint_generation"] + 1:
        _fail("accepted fine must increment the kernel checkpoint generation")
    if current["cursor_generation"] != previous["cursor_generation"] + 1:
        _fail("accepted fine must increment cursor generation")
    if current["attempt_serial"] != previous["attempt_serial"] + 1:
        _fail("accepted fine must increment attempt serial")
    if current["step_index"] != previous["step_index"] + ACCEPTED_SUBSTEPS:
        _fail("accepted fine must advance exactly four substeps")
    expected_serial = _stage_records(member_key)
    if current["transaction_serial"] != previous["transaction_serial"] + expected_serial:
        _fail("accepted fine must advance 20 RK4 or 16 SSPRK3 records")
    if current["descriptor_sha256"] == previous["descriptor_sha256"]:
        _fail("accepted fine must publish a new descriptor")
    if current["payload_sha256"] == previous["payload_sha256"]:
        _fail("accepted fine must publish a new payload")
    if current["physical_state_sha256"] == previous["physical_state_sha256"]:
        _fail("accepted fine must advance accepted physical state")
    if current["cursor_sha256"] == previous["cursor_sha256"]:
        _fail("accepted fine must publish a new full cursor")
    if current["mode"] != FRESH_READY or current["latest_owner"] is not None:
        _fail("accepted fine must return a fresh cursor")
    if current["source_current"] != 0 or current["cfl_current"] != 0:
        _fail("accepted fine must reset per-macro source/CFL counters")
    if current["source_total"] != previous["source_total"] or current["cfl_total"] != previous["cfl_total"]:
        _fail("accepted fine must retain inherited source/CFL totals")
    if current["tdg6_snapshot_sha256"] != previous["tdg6_snapshot_sha256"]:
        _fail("accepted fine must retain the immutable TDG6 sibling")
    descriptor = _descriptor_mapping(current_bundle)
    if descriptor.get("protocol_artifact_id") != RUNTIME_PROTOCOL:
        _fail("fine descriptor protocol differs")
    if descriptor.get("predecessor_cursor_sha256") != previous["cursor_sha256"]:
        _fail("new descriptor must name the old full-cursor hash")
    if descriptor.get("predecessor_descriptor_sha256") != previous["descriptor_sha256"]:
        _fail("new descriptor must name the old accepted descriptor")
    new_cursor = _cursor_from_bundle(current_bundle)
    old_cursor = _cursor_from_bundle(previous_bundle)
    if new_cursor.kernel_transition_parent_sha256 != old_cursor.kernel_transition_digest:
        _fail("fine kernel parent is not the predecessor kernel digest")
    if new_cursor.descriptor_sha256 == old_cursor.descriptor_sha256:
        _fail("fine cursor retained the old descriptor")


def _require_retry_transition(
    previous: Mapping[str, object],
    current: Mapping[str, object],
    *,
    previous_bundle: HLT17GenerationBundle,
    current_bundle: HLT17GenerationBundle,
    kind: str,
) -> None:
    if current["accepted_generation"] != previous["accepted_generation"]:
        _fail("retry cannot increment accepted generation")
    if current["kernel_checkpoint_generation"] != previous["kernel_checkpoint_generation"]:
        _fail("retry cannot increment the kernel checkpoint generation")
    if current["descriptor_sha256"] != previous["descriptor_sha256"]:
        _fail("retry must preserve the accepted descriptor blob")
    if current["payload_sha256"] != previous["payload_sha256"]:
        _fail("retry must preserve the accepted payload blob")
    if current["physical_state_sha256"] != previous["physical_state_sha256"]:
        _fail("retry must preserve accepted physical state")
    if current["accepted_time_hex"] != previous["accepted_time_hex"]:
        _fail("retry must preserve accepted time")
    if current["step_index"] != previous["step_index"]:
        _fail("retry must preserve the accepted step index")
    if current["transaction_serial"] != previous["transaction_serial"]:
        _fail("retry must preserve the accepted transaction serial")
    if current["tdg6_snapshot_sha256"] != previous["tdg6_snapshot_sha256"]:
        _fail("retry must preserve the immutable TDG6 sibling")
    if previous_bundle.descriptor != current_bundle.descriptor:
        _fail("retry descriptor bytes differ")
    if previous_bundle.payload != current_bundle.payload:
        _fail("retry payload bytes differ")
    cursor = _cursor_from_bundle(current_bundle)
    old = _cursor_from_bundle(previous_bundle)
    if kind in CURSOR_CHANGED_RETRY_KINDS:
        if current["cursor_sha256"] == previous["cursor_sha256"]:
            _fail("retry cursor bytes must change")
        if current["cursor_generation"] != previous["cursor_generation"] + 1:
            _fail("retry must increment cursor generation")
        if current["attempt_serial"] != previous["attempt_serial"] + 1:
            _fail("retry must increment attempt serial")
        if cursor.kernel_transition_parent_sha256 != old.kernel_transition_digest:
            _fail("retry kernel parent is not the predecessor kernel digest")
    if kind == "source_retry":
        if cursor.overlay is None or cursor.overlay.owner != SOURCE_OWNER or cursor.overlay.exhausted:
            _fail("source retry cursor overlay differs")
        if current["source_current"] != previous["source_current"] + 1:
            _fail("source retry must increment the source current counter")
        if current["source_total"] != previous["source_total"] + 1:
            _fail("source retry must increment the inherited source total")
        if current["cfl_current"] != previous["cfl_current"] or current["cfl_total"] != previous["cfl_total"]:
            _fail("source retry cannot rewrite CFL counters")
        if cursor.next_attempt_plan() is None:
            _fail("source retry must remain executable")
    elif kind == "cfl_retry":
        if cursor.overlay is None or cursor.overlay.owner != CFL_OWNER or cursor.overlay.exhausted:
            _fail("CFL retry cursor overlay differs")
        if current["cfl_current"] != previous["cfl_current"] + 1:
            _fail("CFL retry must increment the CFL current counter")
        if current["cfl_total"] != previous["cfl_total"] + 1:
            _fail("CFL retry must increment the inherited CFL total")
        if current["source_current"] != previous["source_current"] or current["source_total"] != previous["source_total"]:
            _fail("CFL retry cannot rewrite source counters")
        if cursor.next_attempt_plan() is None:
            _fail("CFL retry must remain executable")
    elif kind == "temporal_retry":
        if cursor.mode != RETRY_PENDING or cursor.temporal is None:
            _fail("temporal retry must retain IMP1 temporal identity")
        if current["source_total"] != previous["source_total"] or current["cfl_total"] != previous["cfl_total"]:
            _fail("temporal retry must keep inherited source/CFL totals")
        if cursor.next_attempt_plan() is None:
            _fail("temporal retry omitted an executable successor plan")
    elif kind == "source_retry_exhausted":
        if cursor.overlay is None or cursor.overlay.owner != SOURCE_OWNER or not cursor.overlay.exhausted:
            _fail("source exhaustion overlay differs")
        if cursor.next_attempt_plan() is not None:
            _fail("source exhaustion must not keep an executable retry cursor")
    elif kind == "cfl_retry_exhausted":
        if cursor.overlay is None or cursor.overlay.owner != CFL_OWNER or not cursor.overlay.exhausted:
            _fail("CFL exhaustion overlay differs")
        if cursor.next_attempt_plan() is not None:
            _fail("CFL exhaustion must not keep an executable retry cursor")


def _require_unchanged_attempt(
    previous: Mapping[str, object],
    current: Mapping[str, object],
    *,
    previous_bundle: HLT17GenerationBundle,
    current_bundle: HLT17GenerationBundle,
    kind: str,
    terminal_evidence: object,
) -> None:
    if current != previous:
        _fail(f"{kind} must not change the accepted member checkpoint state")
    if (
        previous_bundle.descriptor != current_bundle.descriptor
        or previous_bundle.payload != current_bundle.payload
        or previous_bundle.cursor_bytes != current_bundle.cursor_bytes
    ):
        _fail(f"{kind} must not publish a fabricated successor cursor")
    if type(terminal_evidence) is not dict:
        _fail(f"{kind} must carry terminal evidence")
    if terminal_evidence.get("executable_retry_cursor") is not False:
        _fail(f"{kind} evidence cannot advertise an executable retry cursor")
    if terminal_evidence.get("kind") != kind:
        _fail(f"{kind} evidence kind differs")
    if kind == "temporal_retry_exhausted":
        _require_temporal_exhaustion_successor(
            previous_bundle=previous_bundle,
            terminal_evidence=terminal_evidence,
        )
    elif kind == "invalid_premise":
        if terminal_evidence.get("spends_retry") is not False:
            _fail("invalid premise must not spend a retry")
    elif kind == "resource_exhausted":
        if not _text(terminal_evidence.get("reason"), "resource exhaustion reason"):
            _fail("resource exhaustion omitted a reason")


def _require_temporal_exhaustion_successor(
    *,
    previous_bundle: HLT17GenerationBundle,
    terminal_evidence: Mapping[str, object],
) -> None:
    """The evidence ledger must be the exact one-rejection exhausting successor."""

    reason = _text(terminal_evidence.get("reason"), "temporal exhaustion reason")
    if reason not in TEMPORAL_EXHAUSTION_REASONS:
        _fail("temporal exhaustion reason differs")
    ledger_mapping = terminal_evidence.get("updated_imp1_ledger")
    digest = terminal_evidence.get("updated_ledger_sha256")
    if type(ledger_mapping) is not dict:
        _fail("temporal exhaustion must retain the updated IMP1 ledger as evidence")
    try:
        updated = restore_imp1_checkpoint_extension(ledger_mapping)
    except (TypeError, ValueError) as error:
        raise HLT17ProtocolError(
            "temporal exhaustion ledger is not a closed IMP1 ledger"
        ) from error
    encoded = _canonical(ledger_mapping, "temporal exhaustion ledger")
    if _digest(encoded) != _digest_text(digest, "updated_ledger_sha256"):
        _fail("temporal exhaustion ledger hash differs")
    cursor = _cursor_from_bundle(previous_bundle)
    current = cursor.ledger
    if updated.method != current.method:
        _fail("temporal exhaustion ledger method differs")
    if updated.origin_receipt_sha256 != current.origin_receipt_sha256:
        _fail("temporal exhaustion origin receipt differs")
    if updated.origin_state_sha256 != current.origin_state_sha256:
        _fail("temporal exhaustion origin state differs")
    if updated.origin_step_index != current.origin_step_index:
        _fail("temporal exhaustion origin step differs")
    if updated.origin_transaction_serial != current.origin_transaction_serial:
        _fail("temporal exhaustion origin serial differs")
    if updated.inherited_snapshot != current.inherited_snapshot:
        _fail("temporal exhaustion cannot rewrite the inherited TDG6 sibling")
    if not _same_bits(updated.last_accepted_time, current.last_accepted_time):
        _fail("temporal exhaustion cannot advance accepted time")
    if updated.current_state_sha256 != current.current_state_sha256:
        _fail("temporal exhaustion cannot replace the accepted physical state")
    if updated.current_step_index != current.current_step_index:
        _fail("temporal exhaustion cannot advance the step index")
    if updated.current_transaction_serial != current.current_transaction_serial:
        _fail("temporal exhaustion cannot advance the transaction serial")
    if updated.accepted_macro_step_count != current.accepted_macro_step_count:
        _fail("temporal exhaustion cannot change accepted macro-step count")
    if updated.accumulated_debit_vector != current.accumulated_debit_vector:
        _fail("temporal exhaustion cannot rewrite accumulated debit")
    if (
        updated.last_accepted_macro_step_temporal_retry_count
        != current.last_accepted_macro_step_temporal_retry_count
    ):
        _fail("temporal exhaustion cannot rewrite last-accepted temporal retry history")
    current_rejections = current.serialized_temporal_rejections
    updated_rejections = updated.serialized_temporal_rejections
    if updated_rejections[: len(current_rejections)] != current_rejections:
        _fail("temporal exhaustion ledger is not a prefix-preserving successor")
    appended = updated_rejections[len(current_rejections) :]
    if len(appended) != 1:
        _fail("temporal exhaustion must be the exact one-rejection successor")
    if updated.current_macro_step_temporal_retry_count != (
        current.current_macro_step_temporal_retry_count + 1
    ):
        _fail("temporal exhaustion current retry count is not the successor")
    if updated.cumulative_temporal_retry_count != (
        current.cumulative_temporal_retry_count + 1
    ):
        _fail("temporal exhaustion cumulative retry count is not the successor")
    try:
        expected = append_imp1_rejection(current, load_canonical_json(appended[0].encode("ascii")))
    except (CanonicalJSONError, TypeError, ValueError) as error:
        raise HLT17ProtocolError(
            "temporal exhaustion rejection is not a closed successor of this ledger"
        ) from error
    if ledger_sha256(expected) != ledger_sha256(updated):
        _fail("updated IMP1 ledger is not the exact successor of this member cursor ledger")
    last = load_canonical_json(appended[0].encode("ascii"))
    if type(last) is not dict:
        _fail("temporal exhaustion rejection is not a JSON object")
    if reason == "maximum_temporal_retries":
        if current.current_macro_step_temporal_retry_count != TDG6_MAXIMUM_TEMPORAL_RETRIES_PER_STEP:
            _fail("premature temporal exhaustion")
        if updated.current_macro_step_temporal_retry_count <= TDG6_MAXIMUM_TEMPORAL_RETRIES_PER_STEP:
            _fail("premature temporal exhaustion")
    elif reason == "minimum_macro_step":
        attempted = float.fromhex(_hex_time(last.get("attempted_width_hex"), "attempted_width_hex"))
        if attempted / 2.0 >= float(TDG6_MINIMUM_MACRO_STEP):
            _fail("premature temporal exhaustion")
    elif reason != "coordinate_lattice":
        _fail("temporal exhaustion reason differs")


def _attempt_payload(
    *,
    member_key: str,
    previous: Mapping[str, object],
    current: Mapping[str, object],
    kind: str,
    terminal_evidence: object,
) -> dict[str, object]:
    preserved = kind != "accepted_fine"
    executable = kind in {"source_retry", "cfl_retry", "temporal_retry"}
    payload = {
        "member_key": member_key,
        "predecessor_descriptor_sha256": previous["descriptor_sha256"],
        "predecessor_payload_sha256": previous["payload_sha256"],
        "predecessor_cursor_sha256": previous["cursor_sha256"],
        "successor_descriptor_sha256": current["descriptor_sha256"],
        "successor_payload_sha256": current["payload_sha256"],
        "successor_cursor_sha256": current["cursor_sha256"],
        "kernel_transition_parent_sha256": current["kernel_transition_parent_sha256"],
        "kernel_transition_digest": current["kernel_transition_digest"],
        "accepted_generation": current["accepted_generation"],
        "cursor_generation": current["cursor_generation"],
        "attempt_serial": current["attempt_serial"],
        "kernel_checkpoint_generation": current["kernel_checkpoint_generation"],
        "source_current": current["source_current"],
        "source_total": current["source_total"],
        "cfl_current": current["cfl_current"],
        "cfl_total": current["cfl_total"],
        "mode": current["mode"],
        "latest_owner": current["latest_owner"],
        "accepted_time_hex": current["accepted_time_hex"],
        "physical_state_sha256": current["physical_state_sha256"],
        "accepted_state_preserved": preserved,
        "tdg6_snapshot_sha256": current["tdg6_snapshot_sha256"],
        "parked": current["parked"],
        "public_fresh": current["public_fresh"],
        "executable_retry_cursor": executable,
        "step_index": current["step_index"],
        "transaction_serial": current["transaction_serial"],
        "terminal_evidence": None if terminal_evidence is None else dict(terminal_evidence),
    }
    if kind in TERMINAL_KINDS:
        payload["executable_retry_cursor"] = False
        if kind in BLOB_UNCHANGED_KINDS:
            payload["accepted_state_preserved"] = True
    return payload


def _journal_envelope(
    *,
    campaign_id: str,
    kind: str,
    journal_sequence: int,
    journal_parent_sha256: str,
    store_generation: int,
    predecessor_checkpoint_sha256: str,
    payload: Mapping[str, object],
) -> dict[str, object]:
    mapping = {
        "schema": JOURNAL_SCHEMA,
        "schema_version": SCHEMA_VERSION,
        "protocol_artifact_id": PROTOCOL_ARTIFACT_ID,
        "campaign_id": campaign_id,
        "kind": kind,
        "journal_sequence": journal_sequence,
        "journal_parent_sha256": journal_parent_sha256,
        "store_generation": store_generation,
        "predecessor_checkpoint_sha256": predecessor_checkpoint_sha256,
        "payload": dict(payload),
        "campaign_execution_authorized": False,
        "historical_origin_authenticated": False,
        "production_write_authorized": False,
        "physical_classification": False,
        "calibration_claimed": False,
        "eligibility_claimed": False,
    }
    return mapping


def _checkpoint_document(
    *,
    campaign_id: str,
    member_keys: tuple[str, ...],
    store_generation: int,
    journal_sequence: int,
    journal_sha256: str,
    journal_parent_sha256: str,
    predecessor_checkpoint_sha256: str,
    members: Mapping[str, Mapping[str, object]],
    kind: str,
    changed_member_key: str | None,
) -> dict[str, object]:
    complete = _cohort_complete(members)
    if complete:
        if kind != "accepted_fine":
            _fail("first-event completion is not available for this kind")
        disposition = "first_event_complete"
        terminal = True
    elif kind in TERMINAL_KINDS:
        disposition = kind
        terminal = True
    else:
        disposition = "nonterminal"
        terminal = False
        complete = False
    mapping = {
        "schema": CHECKPOINT_SCHEMA,
        "schema_version": SCHEMA_VERSION,
        "protocol_artifact_id": PROTOCOL_ARTIFACT_ID,
        "campaign_id": campaign_id,
        "store_generation": store_generation,
        "journal_sequence": journal_sequence,
        "journal_record_sha256": journal_sha256,
        "journal_parent_sha256": journal_parent_sha256,
        "predecessor_checkpoint_sha256": predecessor_checkpoint_sha256,
        "event_origin_time_hex": EVENT_ORIGIN_HEX,
        "event_target_time_hex": EVENT_TARGET_HEX,
        "member_keys": list(member_keys),
        "members": {key: dict(members[key]) for key in member_keys},
        "kind": kind,
        "changed_member_key": changed_member_key,
        "disposition": disposition,
        "first_event_complete": complete,
        "terminal": terminal,
        "campaign_execution_authorized": False,
        "historical_origin_authenticated": False,
        "production_write_authorized": False,
        "physical_classification": False,
        "calibration_claimed": False,
        "eligibility_claimed": False,
    }
    return mapping


def _manifest_document(
    *,
    campaign_id: str,
    store_generation: int,
    journal_sequence: int,
    journal_sha256: str,
    checkpoint_sha256: str,
    predecessor_checkpoint_sha256: str,
    kind: str,
    changed_member_key: str | None,
    new_blobs: Mapping[str, bytes],
    data_files: Mapping[str, bytes],
) -> dict[str, object]:
    leaves = {path: _digest(payload) for path, payload in data_files.items()}
    mapping = {
        "schema": GENERATION_MANIFEST_SCHEMA,
        "schema_version": SCHEMA_VERSION,
        "protocol_artifact_id": PROTOCOL_ARTIFACT_ID,
        "campaign_id": campaign_id,
        "store_generation": store_generation,
        "journal_sequence": journal_sequence,
        "journal_record_sha256": journal_sha256,
        "checkpoint_sha256": checkpoint_sha256,
        "predecessor_checkpoint_sha256": predecessor_checkpoint_sha256,
        "kind": kind,
        "changed_member_key": changed_member_key,
        "new_blob_sha256s": sorted(new_blobs),
        "leaves": leaves,
        "production_write_authorized": False,
        "campaign_execution_authorized": False,
    }
    _exact_mapping(mapping, _MANIFEST_KEY_SET, "generation manifest")
    return mapping


def build_seed_generation(
    *,
    root: Mapping[str, object],
    bundles: Mapping[str, HLT17GenerationBundle],
) -> HLT17ValidatedGeneration:
    """Validate a complete cohort seed before any directory is visible."""

    root = validate_root_manifest(dict(root))
    member_keys = tuple(root["member_keys"])
    if member_keys != HLT17_SYNTHETIC_MEMBER_KEYS:
        _fail("actual seed bundles must use the synthetic RK4-9/SSPRK3-9 cohort")
    if set(bundles) != set(member_keys):
        _fail("seed bundles do not match the cohort")
    rebuilt = {key: revalidate_bundle(bundles[key]) for key in member_keys}
    states = {key: member_state_from_bundle(rebuilt[key]) for key in member_keys}
    for key in member_keys:
        _require_seed_bundle(rebuilt[key], key)
        if states[key]["source_current"] != 0 or states[key]["cfl_current"] != 0:
            _fail("seed per-macro retry counters must start at zero")
    payload = {"member_keys": list(member_keys), "members": {key: dict(states[key]) for key in member_keys}}
    new_blobs = _new_blobs_for_cohort(
        bundles=rebuilt, predecessor_bundles=None, changed_member_key=None
    )
    return _assemble_generation(
        root=root,
        kind="seed",
        store_generation=0,
        journal_sequence=0,
        journal_parent_sha256=GENESIS_DIGEST,
        predecessor_checkpoint_sha256=GENESIS_DIGEST,
        bundles=rebuilt,
        states=states,
        changed_member_key=None,
        payload=payload,
        new_blobs=new_blobs,
        predecessor=None,
        predecessor_bundles=None,
    )


def build_successor_generation(
    *,
    root: Mapping[str, object],
    predecessor: HLT17ValidatedGeneration,
    predecessor_bundles: Mapping[str, HLT17GenerationBundle],
    kind: str,
    member_key: str,
    successor_bundle: HLT17GenerationBundle,
    terminal_evidence: Mapping[str, object] | None = None,
) -> HLT17ValidatedGeneration:
    """Validate one member attempt against the unique predecessor before write."""

    root = validate_root_manifest(dict(root))
    member_keys = tuple(root["member_keys"])
    if kind not in RECORD_KINDS or kind == "seed":
        _fail("successor kind differs")
    predecessor = _rebind_generation(predecessor)
    closed = _generation_forbids_append(predecessor)
    if closed is not None:
        _fail(closed)
    if member_key not in member_keys:
        _fail("changed member is outside the cohort")
    if set(predecessor_bundles) != set(member_keys):
        _fail("predecessor bundles do not match the cohort")
    prior = {key: revalidate_bundle(predecessor_bundles[key]) for key in member_keys}
    for key in member_keys:
        expected = predecessor.bundles[key]
        if (
            prior[key].descriptor != expected.descriptor
            or prior[key].payload != expected.payload
            or prior[key].cursor_bytes != expected.cursor_bytes
        ):
            _fail("predecessor bundle bytes disagree with the authenticated generation")
    successor = revalidate_bundle(successor_bundle)
    rebuilt = dict(prior)
    rebuilt[member_key] = successor
    prior_states = {key: member_state_from_bundle(prior[key]) for key in member_keys}
    new_states = {key: member_state_from_bundle(rebuilt[key]) for key in member_keys}
    changed = [
        key
        for key in member_keys
        if (
            prior[key].descriptor != rebuilt[key].descriptor
            or prior[key].payload != rebuilt[key].payload
            or prior[key].cursor_bytes != rebuilt[key].cursor_bytes
        )
    ]
    if kind in BLOB_UNCHANGED_KINDS:
        if changed:
            _fail(f"{kind} must not fabricate a successor cursor or new accepted blobs")
    else:
        if changed != [member_key]:
            _fail("exactly one cohort member may change per attempt")
    previous = prior_states[member_key]
    current = new_states[member_key]
    evidence: object
    if kind == "accepted_fine":
        _require_fine_transition(
            previous,
            current,
            previous_bundle=prior[member_key],
            current_bundle=successor,
            member_key=member_key,
        )
        evidence = None
    elif kind in CURSOR_CHANGED_RETRY_KINDS:
        _require_retry_transition(
            previous,
            current,
            previous_bundle=prior[member_key],
            current_bundle=successor,
            kind=kind,
        )
        evidence = None if kind not in TERMINAL_KINDS else {
            "kind": kind,
            "exhausted": True,
            "executable_retry_cursor": False,
            "next_plan": None,
        }
    elif kind in BLOB_UNCHANGED_KINDS:
        evidence = None if terminal_evidence is None else dict(terminal_evidence)
        _require_unchanged_attempt(
            previous,
            current,
            previous_bundle=prior[member_key],
            current_bundle=successor,
            kind=kind,
            terminal_evidence=evidence,
        )
    else:
        _fail("successor kind differs")
        raise AssertionError("unreachable")
    payload = _attempt_payload(
        member_key=member_key,
        previous=previous,
        current=previous if kind in BLOB_UNCHANGED_KINDS else current,
        kind=kind,
        terminal_evidence=evidence,
    )
    new_blobs = _new_blobs_for_cohort(
        bundles=rebuilt, predecessor_bundles=prior, changed_member_key=member_key
    )
    return _assemble_generation(
        root=root,
        kind=kind,
        store_generation=predecessor.store_generation + 1,
        journal_sequence=predecessor.journal_sequence + 1,
        journal_parent_sha256=predecessor.journal_sha256,
        predecessor_checkpoint_sha256=predecessor.checkpoint_sha256,
        bundles=rebuilt,
        states=new_states,
        changed_member_key=member_key,
        payload=payload,
        new_blobs=new_blobs,
        predecessor=predecessor,
        predecessor_bundles=prior,
    )


def _assemble_generation(
    *,
    root: Mapping[str, object],
    kind: str,
    store_generation: int,
    journal_sequence: int,
    journal_parent_sha256: str,
    predecessor_checkpoint_sha256: str,
    bundles: Mapping[str, HLT17GenerationBundle],
    states: Mapping[str, Mapping[str, object]],
    changed_member_key: str | None,
    payload: Mapping[str, object],
    new_blobs: Mapping[str, bytes],
    predecessor: HLT17ValidatedGeneration | None,
    predecessor_bundles: Mapping[str, HLT17GenerationBundle] | None,
) -> HLT17ValidatedGeneration:
    campaign_id = str(root["campaign_id"])
    member_keys = tuple(root["member_keys"])
    if store_generation != journal_sequence:
        _fail("store generation and journal sequence must advance together from zero")
    if predecessor is None:
        if store_generation != 0 or journal_sequence != 0:
            _fail("seed store generation and journal sequence must be zero")
        if journal_parent_sha256 != GENESIS_DIGEST or predecessor_checkpoint_sha256 != GENESIS_DIGEST:
            _fail("seed parents must be the genesis digest")
    else:
        if store_generation != predecessor.store_generation + 1:
            _fail("store generation is not the unique successor")
        if journal_sequence != predecessor.journal_sequence + 1:
            _fail("journal sequence is not the unique successor")
        if journal_parent_sha256 != predecessor.journal_sha256:
            _fail("journal parent is not the predecessor journal hash")
        if predecessor_checkpoint_sha256 != predecessor.checkpoint_sha256:
            _fail("checkpoint parent is not the predecessor checkpoint hash")
        predecessor = _rebind_generation(predecessor)
        closed = _generation_forbids_append(predecessor)
        if closed is not None:
            _fail(closed)
    journal = _validate_journal_envelope(
        _journal_envelope(
            campaign_id=campaign_id,
            kind=kind,
            journal_sequence=journal_sequence,
            journal_parent_sha256=journal_parent_sha256,
            store_generation=store_generation,
            predecessor_checkpoint_sha256=predecessor_checkpoint_sha256,
            payload=payload,
        ),
        campaign_id=campaign_id,
        member_keys=member_keys,
    )
    journal_raw = _canonical(journal, "journal")
    journal_sha256 = _digest(journal_raw)
    checkpoint = _validate_checkpoint_document(
        _checkpoint_document(
            campaign_id=campaign_id,
            member_keys=member_keys,
            store_generation=store_generation,
            journal_sequence=journal_sequence,
            journal_sha256=journal_sha256,
            journal_parent_sha256=journal_parent_sha256,
            predecessor_checkpoint_sha256=predecessor_checkpoint_sha256,
            members=states,
            kind=kind,
            changed_member_key=changed_member_key,
        ),
        campaign_id=campaign_id,
        member_keys=member_keys,
    )
    checkpoint_raw = _canonical(checkpoint, "checkpoint")
    checkpoint_sha256 = _digest(checkpoint_raw)
    if checkpoint["journal_record_sha256"] != journal_sha256:
        _fail("checkpoint does not point to the actual journal record hash")
    checkpoint_text = checkpoint_raw.decode("ascii")
    if "successor_checkpoint_sha256" in journal or checkpoint_sha256 in journal_raw.decode("ascii"):
        _fail("journal record contains the successor checkpoint hash")
    if "checkpoint_sha256" in journal:
        _fail("journal record contains a checkpoint hash")
    del checkpoint_text
    _require_distinct_identities(
        checkpoint_sha256=checkpoint_sha256,
        journal_sha256=journal_sha256,
        states=states,
    )
    data_files = {
        "journal.json": journal_raw,
        "checkpoint.json": checkpoint_raw,
    }
    for digest, raw in new_blobs.items():
        if _digest(raw) != digest:
            _fail("content-addressed blob name disagrees with bytes")
        data_files[_blob_path(digest)] = raw
    manifest = _manifest_document(
        campaign_id=campaign_id,
        store_generation=store_generation,
        journal_sequence=journal_sequence,
        journal_sha256=journal_sha256,
        checkpoint_sha256=checkpoint_sha256,
        predecessor_checkpoint_sha256=predecessor_checkpoint_sha256,
        kind=kind,
        changed_member_key=changed_member_key,
        new_blobs=new_blobs,
        data_files=data_files,
    )
    files = _files_from_parts(
        journal_raw=journal_raw,
        checkpoint_raw=checkpoint_raw,
        manifest_mapping=manifest,
        new_blobs=new_blobs,
    )
    generation = HLT17ValidatedGeneration(
        kind=kind,
        store_generation=store_generation,
        journal_sequence=journal_sequence,
        journal_raw=journal_raw,
        journal_sha256=journal_sha256,
        checkpoint_raw=checkpoint_raw,
        checkpoint_sha256=checkpoint_sha256,
        manifest_raw=files["generation_manifest.json"],
        files=files,
        bundles={key: bundles[key] for key in member_keys},
        changed_member_key=changed_member_key,
        first_event_complete=bool(checkpoint["first_event_complete"]),
        terminal=bool(checkpoint["terminal"]),
        disposition=str(checkpoint["disposition"]),
    )
    return parse_validated_generation(
        generation.files,
        root=root,
        predecessor=predecessor,
        predecessor_bundles=predecessor_bundles,
    )


def parse_validated_generation(
    files: Mapping[str, bytes],
    *,
    root: Mapping[str, object],
    predecessor: HLT17ValidatedGeneration | None,
    predecessor_bundles: Mapping[str, HLT17GenerationBundle] | None,
    expected: HLT17ValidatedGeneration | None = None,
) -> HLT17ValidatedGeneration:
    """Reconstruct one generation from immutable file bytes."""

    root = validate_root_manifest(dict(root))
    member_keys = tuple(root["member_keys"])
    campaign_id = str(root["campaign_id"])
    if type(files) is not dict and not isinstance(files, Mapping):
        _fail("generation files differ")
    inventory = {path: bytes(payload) for path, payload in files.items()}
    if "generation_manifest.json" not in inventory or "checkpoint.json" not in inventory or "journal.json" not in inventory:
        _fail("generation is missing a required document")
    try:
        manifest = load_canonical_json(inventory["generation_manifest.json"])
        checkpoint = load_canonical_json(inventory["checkpoint.json"])
        journal = load_canonical_json(inventory["journal.json"])
    except (CanonicalJSONError, TypeError, ValueError) as error:
        raise HLT17ProtocolError("generation documents are not canonical JSON") from error
    manifest = _exact_mapping(manifest, _MANIFEST_KEY_SET, "generation manifest")
    if manifest["schema"] != GENERATION_MANIFEST_SCHEMA:
        _fail("generation manifest schema differs")
    checkpoint = _validate_checkpoint_document(
        checkpoint, campaign_id=campaign_id, member_keys=member_keys
    )
    journal = _validate_journal_envelope(
        journal, campaign_id=campaign_id, member_keys=member_keys
    )
    journal_raw = inventory["journal.json"]
    checkpoint_raw = inventory["checkpoint.json"]
    journal_sha256 = _digest(journal_raw)
    checkpoint_sha256 = _digest(checkpoint_raw)
    if journal_sha256 != checkpoint["journal_record_sha256"]:
        _fail("checkpoint journal hash disagrees with journal bytes")
    if journal_sha256 != manifest["journal_record_sha256"]:
        _fail("manifest journal hash disagrees with journal bytes")
    if checkpoint_sha256 != manifest["checkpoint_sha256"]:
        _fail("manifest checkpoint hash disagrees with checkpoint bytes")
    if checkpoint_sha256.encode("ascii") in journal_raw:
        _fail("journal record contains the successor checkpoint hash")
    if journal["kind"] != checkpoint["kind"] or journal["kind"] != manifest["kind"]:
        _fail("generation kind disagrees")
    if journal["store_generation"] != checkpoint["store_generation"]:
        _fail("store generation disagrees between journal and checkpoint")
    if journal["journal_sequence"] != checkpoint["journal_sequence"]:
        _fail("journal sequence disagrees between journal and checkpoint")
    if journal["journal_parent_sha256"] != checkpoint["journal_parent_sha256"]:
        _fail("journal parent disagrees between journal and checkpoint")
    if journal["predecessor_checkpoint_sha256"] != checkpoint["predecessor_checkpoint_sha256"]:
        _fail("checkpoint parent disagrees between journal and checkpoint")
    if manifest["predecessor_checkpoint_sha256"] != checkpoint["predecessor_checkpoint_sha256"]:
        _fail("manifest checkpoint parent disagrees")
    leaves = manifest["leaves"]
    if type(leaves) is not dict:
        _fail("manifest leaves differ")
    expected_paths = set(leaves) | {"generation_manifest.json"}
    if set(inventory) != expected_paths:
        _fail("generation file inventory differs from the manifest")
    for path, digest in leaves.items():
        if path not in inventory or _digest(inventory[path]) != _digest_text(digest, "leaf digest"):
            _fail("generation leaf digest differs")
        if path.startswith("blobs/"):
            name = path[len("blobs/") :]
            if name != digest:
                _fail("blob filename is not its content hash")
    kind = str(journal["kind"])
    changed = checkpoint["changed_member_key"]
    new_blob_ids = manifest["new_blob_sha256s"]
    if type(new_blob_ids) is not list or sorted(new_blob_ids) != new_blob_ids:
        _fail("new_blob_sha256s must be a sorted JSON array")
    new_blobs = {
        _digest_text(item, "new blob"): inventory[_blob_path(item)]
        for item in new_blob_ids
    }
    for digest, raw in new_blobs.items():
        if _digest(raw) != digest:
            _fail("new blob bytes disagree with the content hash")
    bundles: dict[str, HLT17GenerationBundle] = {}
    if kind == "seed":
        if predecessor is not None:
            _fail("seed cannot have a predecessor generation")
        if journal["journal_parent_sha256"] != GENESIS_DIGEST:
            _fail("seed journal parent differs")
        members_payload = journal["payload"]["members"]
        for key in member_keys:
            state = checkpoint["members"][key]
            descriptor = new_blobs[str(state["descriptor_sha256"])]
            payload = new_blobs[str(state["payload_sha256"])]
            cursor_bytes = new_blobs[str(state["cursor_sha256"])]
            bundle = bundle_from_bytes(descriptor, payload, cursor_bytes)
            derived = member_state_from_bundle(bundle)
            if derived != state or derived != members_payload[key]:
                _fail("seed member labels disagree with immutable bytes")
            _require_seed_bundle(bundle, key)
            bundles[key] = bundle
    else:
        if predecessor is None or predecessor_bundles is None:
            _fail("successor generation omitted its predecessor")
        predecessor = _rebind_generation(predecessor)
        closed = _generation_forbids_append(predecessor)
        if closed is not None:
            _fail(closed)
        if journal["journal_parent_sha256"] != predecessor.journal_sha256:
            _fail("journal parent is not the unique predecessor journal")
        if checkpoint["predecessor_checkpoint_sha256"] != predecessor.checkpoint_sha256:
            _fail("checkpoint parent is not the unique predecessor checkpoint")
        if journal["store_generation"] != predecessor.store_generation + 1:
            _fail("store generation hole or fork")
        if journal["journal_sequence"] != predecessor.journal_sequence + 1:
            _fail("journal sequence hole or fork")
        prior = {key: revalidate_bundle(predecessor_bundles[key]) for key in member_keys}
        for key in member_keys:
            previous_state = member_state_from_bundle(prior[key])
            listed = checkpoint["members"][key]
            if key != changed:
                if previous_state != listed:
                    _fail("unchanged member checkpoint state differs")
                bundles[key] = prior[key]
                continue
            descriptor_digest = str(listed["descriptor_sha256"])
            payload_digest = str(listed["payload_sha256"])
            cursor_digest = str(listed["cursor_sha256"])
            descriptor = (
                new_blobs[descriptor_digest]
                if descriptor_digest in new_blobs
                else prior[key].descriptor
            )
            payload_bytes = (
                new_blobs[payload_digest]
                if payload_digest in new_blobs
                else prior[key].payload
            )
            cursor_bytes = (
                new_blobs[cursor_digest]
                if cursor_digest in new_blobs
                else prior[key].cursor_bytes
            )
            if _digest(descriptor) != descriptor_digest or _digest(payload_bytes) != payload_digest:
                _fail("changed member blob digest disagrees")
            if _digest(cursor_bytes) != cursor_digest:
                _fail("changed member cursor digest disagrees")
            bundle = bundle_from_bytes(descriptor, payload_bytes, cursor_bytes)
            derived = member_state_from_bundle(bundle)
            if derived != listed:
                _fail("changed member labels disagree with immutable bytes")
            bundles[key] = bundle
        payload = journal["payload"]
        previous = member_state_from_bundle(prior[str(changed)])
        current = member_state_from_bundle(bundles[str(changed)])
        if kind == "accepted_fine":
            _require_fine_transition(
                previous,
                current,
                previous_bundle=prior[str(changed)],
                current_bundle=bundles[str(changed)],
                member_key=str(changed),
            )
        elif kind in CURSOR_CHANGED_RETRY_KINDS:
            _require_retry_transition(
                previous,
                current,
                previous_bundle=prior[str(changed)],
                current_bundle=bundles[str(changed)],
                kind=kind,
            )
        elif kind in BLOB_UNCHANGED_KINDS:
            _require_unchanged_attempt(
                previous,
                current,
                previous_bundle=prior[str(changed)],
                current_bundle=bundles[str(changed)],
                kind=kind,
                terminal_evidence=payload["terminal_evidence"],
            )
        expected_payload = _attempt_payload(
            member_key=str(changed),
            previous=previous,
            current=current if kind not in BLOB_UNCHANGED_KINDS else previous,
            kind=kind,
            terminal_evidence=payload["terminal_evidence"],
        )
        if kind == "accepted_fine":
            expected_payload["accepted_state_preserved"] = False
            expected_payload["executable_retry_cursor"] = False
        if expected_payload != payload:
            _fail("attempt payload disagrees with immutable predecessor/successor bytes")
    _require_distinct_identities(
        checkpoint_sha256=checkpoint_sha256,
        journal_sha256=journal_sha256,
        states=checkpoint["members"],
    )
    derived_states = {
        key: member_state_from_bundle(bundles[key]) for key in member_keys
    }
    if derived_states != checkpoint["members"]:
        _fail("checkpoint member states disagree with immutable bundle bytes")
    derived_complete = _cohort_complete(derived_states)
    if derived_complete != bool(checkpoint["first_event_complete"]):
        _fail("first-event completion disagrees with the reconstructed cohort")
    if derived_complete and (
        kind != "accepted_fine" or checkpoint["disposition"] != "first_event_complete"
    ):
        _fail("first-event completion must be the last accepted fine")
    if (not derived_complete) and checkpoint["disposition"] == "first_event_complete":
        _fail("first-event completion is not available")
    generation = HLT17ValidatedGeneration(
        kind=kind,
        store_generation=int(checkpoint["store_generation"]),
        journal_sequence=int(checkpoint["journal_sequence"]),
        journal_raw=journal_raw,
        journal_sha256=journal_sha256,
        checkpoint_raw=checkpoint_raw,
        checkpoint_sha256=checkpoint_sha256,
        manifest_raw=inventory["generation_manifest.json"],
        files=inventory,
        bundles=bundles,
        changed_member_key=None if changed is None else str(changed),
        first_event_complete=bool(checkpoint["first_event_complete"]),
        terminal=bool(checkpoint["terminal"]),
        disposition=str(checkpoint["disposition"]),
    )
    if expected is not None:
        if generation.files != expected.files:
            _fail("assembled generation files are not redundant")
    return generation


def authenticate_generation_chain(
    *,
    root: Mapping[str, object],
    generations: Sequence[Mapping[str, bytes]],
) -> tuple[HLT17ValidatedGeneration, ...]:
    """Authenticate a unique complete chain with no holes or forks."""

    root = validate_root_manifest(dict(root))
    if not generations:
        _fail("campaign store has no published generations")
    authenticated: list[HLT17ValidatedGeneration] = []
    predecessor: HLT17ValidatedGeneration | None = None
    predecessor_bundles: dict[str, HLT17GenerationBundle] | None = None
    seen_checkpoint: set[str] = set()
    seen_journal: set[str] = set()
    for index, files in enumerate(generations):
        generation = parse_validated_generation(
            files,
            root=root,
            predecessor=predecessor,
            predecessor_bundles=predecessor_bundles,
        )
        if generation.store_generation != index:
            _fail("generation directory order is not the unique store sequence")
        if generation.checkpoint_sha256 in seen_checkpoint or generation.journal_sha256 in seen_journal:
            _fail("generation chain forked")
        seen_checkpoint.add(generation.checkpoint_sha256)
        seen_journal.add(generation.journal_sha256)
        authenticated.append(generation)
        predecessor = generation
        predecessor_bundles = dict(generation.bundles)
        if generation.terminal and index != len(generations) - 1:
            _fail("typed terminal is not the last published generation")
    return tuple(authenticated)


__all__ = [
    "ACCEPTED_SUBSTEPS",
    "ATTEMPT_PAYLOAD_KEYS",
    "BLOB_UNCHANGED_KINDS",
    "CHECKPOINT_KEYS",
    "CHECKPOINT_SCHEMA",
    "COMMON_EVENT_INDEX",
    "CURSOR_CHANGED_RETRY_KINDS",
    "DISPOSITIONS",
    "EVENT_ORIGIN",
    "EVENT_ORIGIN_HEX",
    "EVENT_TARGET",
    "EVENT_TARGET_HEX",
    "FORBIDDEN_JOURNAL_KEYS",
    "FUTURE_ENVIRONMENT_GATE_SEAM",
    "FUTURE_INDEPENDENT_BINDER_SEAM",
    "FUTURE_LIVE_RUNNER_SEAM",
    "FUTURE_PHYSICAL_ORIGIN_SEAM",
    "FUTURE_RESOURCE_GATE_SEAM",
    "FUTURE_SOURCE_MANIFEST_SEAM_TEXT",
    "GENERATION_MANIFEST_KEYS",
    "GENERATION_MANIFEST_SCHEMA",
    "GENESIS_DIGEST",
    "HLT17ProtocolError",
    "HLT17ValidatedGeneration",
    "JOURNAL_ENVELOPE_KEYS",
    "JOURNAL_SCHEMA",
    "MEMBER_STATE_KEYS",
    "PRODUCTION_MEMBER_IDENTITIES",
    "PROTOCOL_ARTIFACT_ID",
    "RECORD_KINDS",
    "RK4_STAGE_RECORDS",
    "ROOT_MANIFEST_KEYS",
    "ROOT_SCHEMA",
    "SCHEMA_VERSION",
    "SEED_PAYLOAD_KEYS",
    "SESSION_CLOSE_KEYS",
    "SESSION_CLOSE_SCHEMA",
    "SESSION_OPEN_KEYS",
    "SESSION_OPEN_SCHEMA",
    "SSPRK3_STAGE_RECORDS",
    "STORE_KIND_SYNTHETIC",
    "TERMINAL_KINDS",
    "TERMINAL_LOCK_KEYS",
    "TERMINAL_LOCK_SCHEMA",
    "authenticate_generation_chain",
    "build_root_manifest",
    "build_seed_generation",
    "build_successor_generation",
    "bundle_from_bytes",
    "cohort_member_keys",
    "encode_session_close",
    "encode_session_open",
    "encode_terminal_lock",
    "generation_directory_name",
    "member_state_from_bundle",
    "parse_generation_directory_name",
    "parse_session_name",
    "parse_validated_generation",
    "production_member_identities",
    "revalidate_bundle",
    "session_close_name",
    "session_open_name",
    "validate_root_manifest",
    "validate_session_close",
    "validate_session_open",
    "validate_terminal_lock",
]
