"""Production-sibling PROTO19/PRO20-EV1 first-event journal and checkpoint schema.

This module is not a boolean flip of the synthetic HLT17/PROTO19 protocol.
It validates the six live C1R1 members, binds an externally supplied
immutable authority receipt, and keeps every generation-level authority or
claim flag false.  Construction of a receipt in tests is not live
authorization.  The module does not write a filesystem, run a campaign, or
treat runner labels as evidence.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from hashlib import sha256
import re
from types import MappingProxyType
from typing import Final

from recursive_horizons.evidence_io import CanonicalJSONError, canonical_json_bytes, load_canonical_json

from .hlt17_admission_runtime import (
    IMPLEMENTATION_IDENTITY_KEYS,
    hlt17_implementation_identity,
)
from .hlt17_imp1_cursor import (
    FRESH_READY,
    HLT17_MEMBER_KEYS,
    IDENTITY_SCOPE_LIVE,
    RETRY_PENDING,
    hlt17_member_spec,
)
from .hlt17_member_checkpoint import HLT17GenerationBundle
from .hlt17_member_codec import ORIGIN_PROTOCOL
from .protocol_v19 import (
    ATTEMPT_PAYLOAD_KEYS,
    BLOB_UNCHANGED_KINDS,
    COMMON_EVENT_INDEX,
    CURSOR_CHANGED_RETRY_KINDS,
    DISPOSITIONS,
    EVENT_ORIGIN,
    EVENT_ORIGIN_HEX,
    EVENT_TARGET,
    EVENT_TARGET_HEX,
    FORBIDDEN_JOURNAL_KEYS,
    GENESIS_DIGEST,
    HLT17ProtocolError,
    MEMBER_STATE_KEYS,
    RECORD_KINDS,
    SEED_PAYLOAD_KEYS,
    TERMINAL_KINDS,
    _attempt_payload,
    _blob_path,
    _cohort_complete,
    _cursor_from_bundle,
    _descriptor_mapping,
    _new_blobs_for_cohort,
    _require_fine_transition,
    _require_retry_transition,
    _require_unchanged_attempt,
    bundle_from_bytes,
    generation_directory_name,
    member_state_from_bundle,
    parse_generation_directory_name as _parse_generation_directory_name,
    parse_session_name as _parse_session_name,
    production_member_identities,
    revalidate_bundle,
    session_close_name,
    session_open_name,
)

if EVENT_ORIGIN.hex() != EVENT_ORIGIN_HEX or EVENT_TARGET.hex() != EVENT_TARGET_HEX:
    raise RuntimeError("PRO20-EV1 first-event times disagree with t=23/16 and t=3/2")


PROTOCOL_ARTIFACT_ID: Final[str] = "FGC-2-SF1-PROTO19"
FREEZE_ARTIFACT_ID: Final[str] = "FGC-1-PRO20-EV1-FRZ1"
ARTIFACT_ID: Final[str] = "FGC-1-PRO20-EV1-FRZ1"
SCHEMA_VERSION: Final[int] = 1
STORE_KIND_PRODUCTION: Final[str] = "production_event1"
PRODUCTION_NAMESPACE: Final[str] = "runs/fgc-2-sf1/pro20-event1/calibration"
RECEIPT_SCHEMA: Final[str] = "FGC-1-PRO20-EV1-authority-receipt-v1"
ROOT_SCHEMA: Final[str] = "FGC-1-PRO20-EV1-campaign-root-v1"
JOURNAL_SCHEMA: Final[str] = "FGC-1-PRO20-EV1-campaign-journal-v1"
CHECKPOINT_SCHEMA: Final[str] = "FGC-1-PRO20-EV1-campaign-checkpoint-v1"
GENERATION_MANIFEST_SCHEMA: Final[str] = "FGC-1-PRO20-EV1-generation-manifest-v1"
SESSION_OPEN_SCHEMA: Final[str] = "FGC-1-PRO20-EV1-writer-session-open-v1"
SESSION_CLOSE_SCHEMA: Final[str] = "FGC-1-PRO20-EV1-writer-session-close-v1"
TERMINAL_LOCK_SCHEMA: Final[str] = "FGC-1-PRO20-EV1-terminal-lock-v1"
_SHA256 = re.compile(r"\A[0-9a-f]{64}\Z")

PINNED_PROTOCOL_V19_PRIVATE_DEPENDENCIES: Final[tuple[str, ...]] = (
    "protocol_v19._attempt_payload",
    "protocol_v19._blob_path",
    "protocol_v19._cohort_complete",
    "protocol_v19._cursor_from_bundle",
    "protocol_v19._descriptor_mapping",
    "protocol_v19._new_blobs_for_cohort",
    "protocol_v19._require_fine_transition",
    "protocol_v19._require_retry_transition",
    "protocol_v19._require_unchanged_attempt",
)

RUNNER_SEAM: Final[str] = (
    "a committed live runner is an external owner; this slice does not execute "
    "or treat runner labels as evidence"
)
EXACT_DELTA_AUTHORITY_SEAM: Final[str] = (
    "exact-delta freeze authority is an external owner; this receipt type is "
    "not that freeze and construction is not live authorization"
)
RESOURCE_GATE_SEAM: Final[str] = (
    "largest-grid and physical-RHS resource qualification is an external gate; "
    "this slice does not execute or authorize a production preparation"
)
INDEPENDENT_BINDER_SEAM: Final[str] = (
    "an independent PREF1 binder must reconstruct evidence without trusting "
    "runner labels; this slice is not that binder"
)

NONCLAIM_FLAG_NAMES: Final[tuple[str, ...]] = (
    "historical_origin_authenticated",
    "source_authenticated",
    "environment_authenticated",
    "campaign_execution_authorized",
    "production_write_authorized",
    "live_runner_authorized",
    "calibration_claimed",
    "eligibility_claimed",
    "physical_classification",
    "calibration_eligible",
    "sgbl_claimed",
    "def1_claimed",
    "candidate_claimed",
    "mechanism_claimed",
    "physics_claimed",
)

RECEIPT_KEYS: Final[tuple[str, ...]] = (
    "schema",
    "schema_version",
    "protocol_artifact_id",
    "freeze_artifact_id",
    "campaign_id",
    "namespace",
    "member_keys",
    "event_origin_time_hex",
    "event_target_time_hex",
    "common_event_index",
    "authority_sha256",
    "implementation_sha256",
    "config_sha256",
    "source_sha256",
    "origin_sha256",
    "environment_sha256",
    "implementation_identity",
    *NONCLAIM_FLAG_NAMES,
    "runner_seam",
    "exact_delta_authority_seam",
    "resource_gate_seam",
    "independent_binder_seam",
)
ROOT_MANIFEST_KEYS: Final[tuple[str, ...]] = (
    "schema",
    "schema_version",
    "protocol_artifact_id",
    "freeze_artifact_id",
    "campaign_id",
    "store_kind",
    "namespace",
    "member_keys",
    "event_origin_time_hex",
    "event_target_time_hex",
    "common_event_index",
    "authority_receipt_sha256",
    "authority_sha256",
    "implementation_sha256",
    "config_sha256",
    "source_sha256",
    "origin_sha256",
    "environment_sha256",
    "implementation_identity",
    *NONCLAIM_FLAG_NAMES,
    "runner_seam",
    "exact_delta_authority_seam",
    "resource_gate_seam",
    "independent_binder_seam",
)
JOURNAL_ENVELOPE_KEYS: Final[tuple[str, ...]] = (
    "schema",
    "schema_version",
    "protocol_artifact_id",
    "freeze_artifact_id",
    "campaign_id",
    "kind",
    "journal_sequence",
    "journal_parent_sha256",
    "store_generation",
    "predecessor_checkpoint_sha256",
    "payload",
    *NONCLAIM_FLAG_NAMES,
)
CHECKPOINT_KEYS: Final[tuple[str, ...]] = (
    "schema",
    "schema_version",
    "protocol_artifact_id",
    "freeze_artifact_id",
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
    *NONCLAIM_FLAG_NAMES,
)
GENERATION_MANIFEST_KEYS: Final[tuple[str, ...]] = (
    "schema",
    "schema_version",
    "protocol_artifact_id",
    "freeze_artifact_id",
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
    "calibration_claimed",
    "eligibility_claimed",
    "physical_classification",
)
SESSION_OPEN_KEYS: Final[tuple[str, ...]] = (
    "schema",
    "schema_version",
    "protocol_artifact_id",
    "freeze_artifact_id",
    "campaign_id",
    "session_index",
    "owner_token",
    "host",
    "pid",
    "thread_id",
    "exclusion_lock_sha256",
    "authority_receipt_sha256",
    "namespace",
    "production_write_authorized",
    "campaign_execution_authorized",
)
SESSION_CLOSE_KEYS: Final[tuple[str, ...]] = (
    "schema",
    "schema_version",
    "protocol_artifact_id",
    "freeze_artifact_id",
    "campaign_id",
    "session_index",
    "session_open_sha256",
    "owner_token",
    "authority_receipt_sha256",
)
TERMINAL_LOCK_KEYS: Final[tuple[str, ...]] = (
    "schema",
    "schema_version",
    "protocol_artifact_id",
    "freeze_artifact_id",
    "campaign_id",
    "checkpoint_sha256",
    "journal_record_sha256",
    "store_generation",
    "disposition",
    "changed_member_key",
    "first_event_complete",
    "calibration_claimed",
    "eligibility_claimed",
    "physical_classification",
    "calibration_eligible",
    "sgbl_claimed",
    "def1_claimed",
    "candidate_claimed",
    "mechanism_claimed",
    "physics_claimed",
)

_RECEIPT_KEY_SET = frozenset(RECEIPT_KEYS)
_ROOT_KEY_SET = frozenset(ROOT_MANIFEST_KEYS)
_JOURNAL_KEY_SET = frozenset(JOURNAL_ENVELOPE_KEYS)
_CHECKPOINT_KEY_SET = frozenset(CHECKPOINT_KEYS)
_MANIFEST_KEY_SET = frozenset(GENERATION_MANIFEST_KEYS)
_SESSION_OPEN_KEY_SET = frozenset(SESSION_OPEN_KEYS)
_SESSION_CLOSE_KEY_SET = frozenset(SESSION_CLOSE_KEYS)
_TERMINAL_LOCK_KEY_SET = frozenset(TERMINAL_LOCK_KEYS)
_SEED_PAYLOAD_KEY_SET = frozenset(SEED_PAYLOAD_KEYS)
_ATTEMPT_PAYLOAD_KEY_SET = frozenset(ATTEMPT_PAYLOAD_KEYS)
_MEMBER_STATE_KEY_SET = frozenset(MEMBER_STATE_KEYS)


class PRO20EV1ProtocolError(ValueError):
    """Malformed PRO20-EV1 receipt, checkpoint, journal, manifest, or transition."""


def _fail(message: str) -> None:
    raise PRO20EV1ProtocolError(message)


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
        raise PRO20EV1ProtocolError(f"{label} is not a binary64 hex time") from error
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
        raise PRO20EV1ProtocolError(f"{label} is not canonical JSON") from error
    if loaded != value:
        _fail(f"{label} is not a redundant canonical mapping")
    return raw


def _wrap_hlt17(error: Exception) -> None:
    raise PRO20EV1ProtocolError(str(error)) from error


def _member_state(bundle: object) -> dict[str, object]:
    try:
        return member_state_from_bundle(bundle)
    except (HLT17ProtocolError, TypeError, ValueError) as error:
        _wrap_hlt17(error)
        raise AssertionError("unreachable")


def _revalidate(bundle: object) -> HLT17GenerationBundle:
    try:
        return revalidate_bundle(bundle)
    except (HLT17ProtocolError, TypeError, ValueError) as error:
        _wrap_hlt17(error)
        raise AssertionError("unreachable")


def _from_bytes(descriptor: bytes, payload: bytes, cursor_bytes: bytes) -> HLT17GenerationBundle:
    try:
        return bundle_from_bytes(descriptor, payload, cursor_bytes)
    except (HLT17ProtocolError, TypeError, ValueError) as error:
        _wrap_hlt17(error)
        raise AssertionError("unreachable")


def parse_generation_directory_name(name: object) -> int:
    try:
        return _parse_generation_directory_name(name)
    except HLT17ProtocolError as error:
        raise PRO20EV1ProtocolError(str(error)) from error


def parse_session_name(name: object) -> tuple[int, str]:
    try:
        return _parse_session_name(name)
    except HLT17ProtocolError as error:
        raise PRO20EV1ProtocolError(str(error)) from error


def live_member_keys(value: object) -> tuple[str, ...]:
    """Exact six live production keys in canonical cohort order."""

    if type(value) is not list or not value:
        _fail("member_keys must be a nonempty JSON array")
    keys = tuple(_text(item, "member_key") for item in value)
    try:
        production_member_identities()
    except HLT17ProtocolError as error:
        _wrap_hlt17(error)
    if keys != HLT17_MEMBER_KEYS:
        _fail("cohort must be the exact six live production members in canonical order")
    for key in keys:
        spec = hlt17_member_spec(key)
        if spec.identity_scope != IDENTITY_SCOPE_LIVE:
            _fail("production cohort used a synthetic identity")
    return keys


def _campaign_id(value: object) -> str:
    campaign_id = _text(value, "campaign_id")
    if "/" in campaign_id or "\\" in campaign_id or ".." in campaign_id:
        _fail("campaign_id is not a safe identifier")
    return campaign_id


def _implementation_identity(value: object) -> dict[str, str]:
    expected = hlt17_implementation_identity()
    if type(value) is not dict:
        _fail("implementation_identity must be a JSON object")
    if set(value) != set(IMPLEMENTATION_IDENTITY_KEYS):
        _fail("implementation_identity keys differ")
    mapping = {name: _text(value[name], name) for name in IMPLEMENTATION_IDENTITY_KEYS}
    if mapping != expected:
        _fail("implementation_identity is not the fixed C1R1 identity")
    return expected


def _nonclaim_flags(mapping: Mapping[str, object]) -> None:
    for name in NONCLAIM_FLAG_NAMES:
        if name in mapping:
            _bool_false(mapping[name], name)


def _seams(mapping: Mapping[str, object]) -> None:
    expected = {
        "runner_seam": RUNNER_SEAM,
        "exact_delta_authority_seam": EXACT_DELTA_AUTHORITY_SEAM,
        "resource_gate_seam": RESOURCE_GATE_SEAM,
        "independent_binder_seam": INDEPENDENT_BINDER_SEAM,
    }
    for name, text in expected.items():
        if mapping[name] != text:
            _fail(f"{name} differs")


def _event_bounds(mapping: Mapping[str, object]) -> None:
    if _hex_time(mapping["event_origin_time_hex"], "event_origin_time_hex") != EVENT_ORIGIN_HEX:
        _fail("first-event origin must be t=23/16")
    if _hex_time(mapping["event_target_time_hex"], "event_target_time_hex") != EVENT_TARGET_HEX:
        _fail("first-event target must be t=3/2")
    if _nonneg_int(mapping["common_event_index"], "common_event_index") != COMMON_EVENT_INDEX:
        _fail("first-event index must be 23")


def _identity_hashes(mapping: Mapping[str, object]) -> None:
    for name in (
        "authority_sha256",
        "implementation_sha256",
        "config_sha256",
        "source_sha256",
        "origin_sha256",
        "environment_sha256",
    ):
        _digest_text(mapping[name], name)
    identities = [
        mapping["authority_sha256"],
        mapping["implementation_sha256"],
        mapping["config_sha256"],
        mapping["source_sha256"],
        mapping["origin_sha256"],
        mapping["environment_sha256"],
    ]
    if len(set(identities)) != 6:
        _fail("authority/implementation/config/source/origin/environment identities must be distinct")


def _protocol_freeze(mapping: Mapping[str, object]) -> None:
    if mapping["protocol_artifact_id"] != PROTOCOL_ARTIFACT_ID:
        _fail("protocol differs")
    if mapping["freeze_artifact_id"] != FREEZE_ARTIFACT_ID:
        _fail("freeze artifact differs")
    if mapping["namespace"] != PRODUCTION_NAMESPACE:
        _fail("namespace differs")


def _optional_member_key(value: object, allowed: tuple[str, ...]) -> str | None:
    if value is None:
        return None
    text = _text(value, "member_key")
    if text not in allowed:
        _fail("member_key is outside the cohort")
    return text


def validate_authority_receipt(value: object) -> dict[str, object]:
    """Validate an external receipt mapping. This is not live authorization."""

    item = _exact_mapping(value, _RECEIPT_KEY_SET, "authority receipt")
    if item["schema"] != RECEIPT_SCHEMA or item["schema_version"] != SCHEMA_VERSION:
        _fail("authority receipt schema differs")
    _protocol_freeze(item)
    item["campaign_id"] = _campaign_id(item["campaign_id"])
    item["member_keys"] = list(live_member_keys(item["member_keys"]))
    _event_bounds(item)
    _identity_hashes(item)
    item["implementation_identity"] = _implementation_identity(item["implementation_identity"])
    _nonclaim_flags(item)
    _seams(item)
    _canonical(item, "authority receipt")
    return item


def build_authority_receipt(
    *,
    campaign_id: str,
    authority_sha256: str,
    implementation_sha256: str,
    config_sha256: str,
    source_sha256: str,
    origin_sha256: str,
    environment_sha256: str,
    implementation_identity: Mapping[str, str] | None = None,
) -> "PRO20EV1AuthorityReceipt":
    """Construct a validated receipt type. Construction is not live authorization."""

    mapping = {
        "schema": RECEIPT_SCHEMA,
        "schema_version": SCHEMA_VERSION,
        "protocol_artifact_id": PROTOCOL_ARTIFACT_ID,
        "freeze_artifact_id": FREEZE_ARTIFACT_ID,
        "campaign_id": campaign_id,
        "namespace": PRODUCTION_NAMESPACE,
        "member_keys": list(HLT17_MEMBER_KEYS),
        "event_origin_time_hex": EVENT_ORIGIN_HEX,
        "event_target_time_hex": EVENT_TARGET_HEX,
        "common_event_index": COMMON_EVENT_INDEX,
        "authority_sha256": authority_sha256,
        "implementation_sha256": implementation_sha256,
        "config_sha256": config_sha256,
        "source_sha256": source_sha256,
        "origin_sha256": origin_sha256,
        "environment_sha256": environment_sha256,
        "implementation_identity": dict(
            implementation_identity or hlt17_implementation_identity()
        ),
        **{name: False for name in NONCLAIM_FLAG_NAMES},
        "runner_seam": RUNNER_SEAM,
        "exact_delta_authority_seam": EXACT_DELTA_AUTHORITY_SEAM,
        "resource_gate_seam": RESOURCE_GATE_SEAM,
        "independent_binder_seam": INDEPENDENT_BINDER_SEAM,
    }
    return PRO20EV1AuthorityReceipt.from_mapping(mapping)


@dataclass(frozen=True, slots=True)
class PRO20EV1AuthorityReceipt:
    """Externally supplied immutable receipt. Construction is not live authorization."""

    mapping: Mapping[str, object]
    raw: bytes
    sha256: str

    def __post_init__(self) -> None:
        if type(self.raw) is not bytes:
            raise TypeError("authority receipt raw must be bytes")
        item = validate_authority_receipt(dict(self.mapping))
        raw = _canonical(item, "authority receipt")
        digest = _digest(raw)
        if raw != self.raw or digest != self.sha256:
            _fail("authority receipt bytes disagree with the validated mapping")
        object.__setattr__(self, "mapping", MappingProxyType(item))
        object.__setattr__(self, "raw", raw)
        object.__setattr__(self, "sha256", digest)

    @classmethod
    def from_mapping(cls, value: object) -> "PRO20EV1AuthorityReceipt":
        item = validate_authority_receipt(value)
        raw = _canonical(item, "authority receipt")
        return cls(mapping=item, raw=raw, sha256=_digest(raw))

    def as_mapping(self) -> dict[str, object]:
        return validate_authority_receipt(dict(self.mapping))


def _receipt_from_root_fields(item: Mapping[str, object]) -> dict[str, object]:
    mapping = {
        "schema": RECEIPT_SCHEMA,
        "schema_version": SCHEMA_VERSION,
        "protocol_artifact_id": item["protocol_artifact_id"],
        "freeze_artifact_id": item["freeze_artifact_id"],
        "campaign_id": item["campaign_id"],
        "namespace": item["namespace"],
        "member_keys": list(item["member_keys"]),
        "event_origin_time_hex": item["event_origin_time_hex"],
        "event_target_time_hex": item["event_target_time_hex"],
        "common_event_index": item["common_event_index"],
        "authority_sha256": item["authority_sha256"],
        "implementation_sha256": item["implementation_sha256"],
        "config_sha256": item["config_sha256"],
        "source_sha256": item["source_sha256"],
        "origin_sha256": item["origin_sha256"],
        "environment_sha256": item["environment_sha256"],
        "implementation_identity": dict(item["implementation_identity"]),
        **{name: item[name] for name in NONCLAIM_FLAG_NAMES},
        "runner_seam": item["runner_seam"],
        "exact_delta_authority_seam": item["exact_delta_authority_seam"],
        "resource_gate_seam": item["resource_gate_seam"],
        "independent_binder_seam": item["independent_binder_seam"],
    }
    return validate_authority_receipt(mapping)


def validate_root_manifest(value: object) -> dict[str, object]:
    item = _exact_mapping(value, _ROOT_KEY_SET, "root manifest")
    if item["schema"] != ROOT_SCHEMA or item["schema_version"] != SCHEMA_VERSION:
        _fail("root schema differs")
    if item["store_kind"] != STORE_KIND_PRODUCTION:
        _fail("store_kind must be production_event1")
    _protocol_freeze(item)
    item["campaign_id"] = _campaign_id(item["campaign_id"])
    item["member_keys"] = list(live_member_keys(item["member_keys"]))
    _event_bounds(item)
    _identity_hashes(item)
    item["implementation_identity"] = _implementation_identity(item["implementation_identity"])
    _nonclaim_flags(item)
    _seams(item)
    receipt = _receipt_from_root_fields(item)
    expected = _digest(_canonical(receipt, "authority receipt"))
    if _digest_text(item["authority_receipt_sha256"], "authority_receipt_sha256") != expected:
        _fail("root authority_receipt_sha256 disagrees with the bound receipt")
    _canonical(item, "root manifest")
    return item


def build_root_manifest(receipt: object) -> dict[str, object]:
    if not isinstance(receipt, PRO20EV1AuthorityReceipt):
        receipt = PRO20EV1AuthorityReceipt.from_mapping(receipt)
    receipt.__post_init__()
    body = receipt.as_mapping()
    mapping = {
        "schema": ROOT_SCHEMA,
        "schema_version": SCHEMA_VERSION,
        "protocol_artifact_id": body["protocol_artifact_id"],
        "freeze_artifact_id": body["freeze_artifact_id"],
        "campaign_id": body["campaign_id"],
        "store_kind": STORE_KIND_PRODUCTION,
        "namespace": body["namespace"],
        "member_keys": list(body["member_keys"]),
        "event_origin_time_hex": body["event_origin_time_hex"],
        "event_target_time_hex": body["event_target_time_hex"],
        "common_event_index": body["common_event_index"],
        "authority_receipt_sha256": receipt.sha256,
        "authority_sha256": body["authority_sha256"],
        "implementation_sha256": body["implementation_sha256"],
        "config_sha256": body["config_sha256"],
        "source_sha256": body["source_sha256"],
        "origin_sha256": body["origin_sha256"],
        "environment_sha256": body["environment_sha256"],
        "implementation_identity": dict(body["implementation_identity"]),
        **{name: False for name in NONCLAIM_FLAG_NAMES},
        "runner_seam": RUNNER_SEAM,
        "exact_delta_authority_seam": EXACT_DELTA_AUTHORITY_SEAM,
        "resource_gate_seam": RESOURCE_GATE_SEAM,
        "independent_binder_seam": INDEPENDENT_BINDER_SEAM,
    }
    return validate_root_manifest(mapping)


def authority_receipt_from_root(root: Mapping[str, object]) -> PRO20EV1AuthorityReceipt:
    """Rebuild the bound receipt from a validated production root manifest."""

    item = validate_root_manifest(dict(root))
    return PRO20EV1AuthorityReceipt.from_mapping(_receipt_from_root_fields(item))


def encode_session_open(
    *,
    campaign_id: str,
    session_index: int,
    owner_token: str,
    host: str,
    pid: int,
    thread_id: int,
    exclusion_lock_sha256: str,
    authority_receipt_sha256: str,
) -> dict[str, object]:
    mapping = {
        "schema": SESSION_OPEN_SCHEMA,
        "schema_version": SCHEMA_VERSION,
        "protocol_artifact_id": PROTOCOL_ARTIFACT_ID,
        "freeze_artifact_id": FREEZE_ARTIFACT_ID,
        "campaign_id": _campaign_id(campaign_id),
        "session_index": _nonneg_int(session_index, "session_index"),
        "owner_token": _text(owner_token, "owner_token"),
        "host": _text(host, "host"),
        "pid": _nonneg_int(pid, "pid"),
        "thread_id": _nonneg_int(thread_id, "thread_id"),
        "exclusion_lock_sha256": _digest_text(exclusion_lock_sha256, "exclusion_lock_sha256"),
        "authority_receipt_sha256": _digest_text(
            authority_receipt_sha256, "authority_receipt_sha256"
        ),
        "namespace": PRODUCTION_NAMESPACE,
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
    authority_receipt_sha256: str,
) -> dict[str, object]:
    mapping = {
        "schema": SESSION_CLOSE_SCHEMA,
        "schema_version": SCHEMA_VERSION,
        "protocol_artifact_id": PROTOCOL_ARTIFACT_ID,
        "freeze_artifact_id": FREEZE_ARTIFACT_ID,
        "campaign_id": _campaign_id(campaign_id),
        "session_index": _nonneg_int(session_index, "session_index"),
        "session_open_sha256": _digest_text(session_open_sha256, "session_open_sha256"),
        "owner_token": _text(owner_token, "owner_token"),
        "authority_receipt_sha256": _digest_text(
            authority_receipt_sha256, "authority_receipt_sha256"
        ),
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
    first_event_complete: bool,
) -> dict[str, object]:
    if disposition not in DISPOSITIONS or disposition == "nonterminal":
        _fail("terminal lock disposition differs")
    complete = _bool(first_event_complete, "first_event_complete")
    if complete is True and disposition != "first_event_complete":
        _fail("successful event completion lock disposition differs")
    if disposition == "first_event_complete" and complete is not True:
        _fail("first_event_complete lock omitted the completion flag")
    mapping = {
        "schema": TERMINAL_LOCK_SCHEMA,
        "schema_version": SCHEMA_VERSION,
        "protocol_artifact_id": PROTOCOL_ARTIFACT_ID,
        "freeze_artifact_id": FREEZE_ARTIFACT_ID,
        "campaign_id": _campaign_id(campaign_id),
        "checkpoint_sha256": _digest_text(checkpoint_sha256, "checkpoint_sha256"),
        "journal_record_sha256": _digest_text(
            journal_record_sha256, "journal_record_sha256"
        ),
        "store_generation": _nonneg_int(store_generation, "store_generation"),
        "disposition": _text(disposition, "disposition"),
        "changed_member_key": changed_member_key,
        "first_event_complete": complete,
        "calibration_claimed": False,
        "eligibility_claimed": False,
        "physical_classification": False,
        "calibration_eligible": False,
        "sgbl_claimed": False,
        "def1_claimed": False,
        "candidate_claimed": False,
        "mechanism_claimed": False,
        "physics_claimed": False,
    }
    _exact_mapping(mapping, _TERMINAL_LOCK_KEY_SET, "terminal lock")
    _canonical(mapping, "terminal lock")
    return mapping


def validate_session_open(value: object) -> dict[str, object]:
    item = _exact_mapping(value, _SESSION_OPEN_KEY_SET, "session open")
    if item["schema"] != SESSION_OPEN_SCHEMA or item["schema_version"] != SCHEMA_VERSION:
        _fail("session open schema differs")
    expected = encode_session_open(
        campaign_id=item["campaign_id"],
        session_index=item["session_index"],
        owner_token=item["owner_token"],
        host=item["host"],
        pid=item["pid"],
        thread_id=item["thread_id"],
        exclusion_lock_sha256=item["exclusion_lock_sha256"],
        authority_receipt_sha256=item["authority_receipt_sha256"],
    )
    if item != expected:
        _fail("session open differs from its canonical non-authorizing encoding")
    return item


def validate_session_close(value: object) -> dict[str, object]:
    item = _exact_mapping(value, _SESSION_CLOSE_KEY_SET, "session close")
    if item["schema"] != SESSION_CLOSE_SCHEMA or item["schema_version"] != SCHEMA_VERSION:
        _fail("session close schema differs")
    expected = encode_session_close(
        campaign_id=item["campaign_id"],
        session_index=item["session_index"],
        session_open_sha256=item["session_open_sha256"],
        owner_token=item["owner_token"],
        authority_receipt_sha256=item["authority_receipt_sha256"],
    )
    if item != expected:
        _fail("session close differs from its canonical receipt binding")
    return item


def validate_terminal_lock(value: object) -> dict[str, object]:
    item = _exact_mapping(value, _TERMINAL_LOCK_KEY_SET, "terminal lock")
    if item["schema"] != TERMINAL_LOCK_SCHEMA or item["schema_version"] != SCHEMA_VERSION:
        _fail("terminal lock schema differs")
    expected = encode_terminal_lock(
        campaign_id=item["campaign_id"],
        checkpoint_sha256=item["checkpoint_sha256"],
        journal_record_sha256=item["journal_record_sha256"],
        store_generation=item["store_generation"],
        disposition=item["disposition"],
        changed_member_key=item["changed_member_key"],
        first_event_complete=item["first_event_complete"],
    )
    if item != expected:
        _fail("terminal lock differs from its canonical nonclaim encoding")
    return item


def _require_live_seed_bundle(bundle: HLT17GenerationBundle, member_key: str) -> None:
    try:
        cursor = _cursor_from_bundle(bundle)
        descriptor = _descriptor_mapping(bundle)
    except HLT17ProtocolError as error:
        _wrap_hlt17(error)
        raise AssertionError("unreachable")
    identity = descriptor["runtime_identity"]
    if type(identity) is not dict:
        _fail("runtime identity differs")
    if identity["member_key"] != member_key or cursor.member_key != member_key:
        _fail("seed bundle member_key differs")
    spec = hlt17_member_spec(member_key)
    if spec.identity_scope != IDENTITY_SCOPE_LIVE or member_key not in HLT17_MEMBER_KEYS:
        _fail("seed refuses synthetic or diagnostic members")
    if identity.get("scope") != IDENTITY_SCOPE_LIVE:
        _fail("seed bundle is not a live C1R1 member")
    if identity.get("amplitude") != "3":
        _fail("live seed amplitude differs")
    if identity.get("method") != spec.method_label:
        _fail("live method identity differs")
    if identity.get("point_count") != spec.point_count:
        _fail("live point count differs")
    if identity.get("spatial_order") != (4 if spec.method_label == "RK4" else 2):
        _fail("live spatial operator identity differs")
    if cursor.member.identity_scope != IDENTITY_SCOPE_LIVE:
        _fail("seed cursor is not a live C1R1 member")
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
    if cursor.accepted_time.hex() == cursor.event_target.hex():
        _fail("seed cannot already be parked at the event target")


def _generation_forbids_append(generation: "PRO20EV1ValidatedGeneration") -> str | None:
    try:
        checkpoint = load_canonical_json(generation.checkpoint_raw)
    except (CanonicalJSONError, TypeError, ValueError) as error:
        raise PRO20EV1ProtocolError("predecessor checkpoint is not canonical JSON") from error
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


def next_required_member_key(generation: "PRO20EV1ValidatedGeneration") -> str:
    """Canonical cohort order; pending retry ownership may force the same member."""

    generation = _rebind_generation(generation)
    closed = _generation_forbids_append(generation)
    if closed is not None:
        _fail(closed)
    member_keys = tuple(live_member_keys(generation.checkpoint_mapping["member_keys"]))
    states = {key: _member_state(generation.bundles[key]) for key in member_keys}
    pending = [
        key
        for key in member_keys
        if states[key]["mode"] != FRESH_READY or states[key]["latest_owner"] is not None
    ]
    if pending:
        if len(pending) != 1:
            _fail("multiple members hold pending retry ownership")
        return pending[0]
    for key in member_keys:
        if states[key]["parked"] is not True:
            return key
    _fail("nonterminal generation has no required member")
    raise AssertionError("unreachable")


@dataclass(frozen=True, slots=True)
class PRO20EV1ValidatedGeneration:
    """Closed production generation: journal, checkpoint, manifest, and blobs."""

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
                {key: _revalidate(bundle) for key, bundle in dict(self.bundles).items()}
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


def _rebind_generation(generation: object) -> PRO20EV1ValidatedGeneration:
    if not isinstance(generation, PRO20EV1ValidatedGeneration):
        raise TypeError("generation must be PRO20EV1ValidatedGeneration")
    generation.__post_init__()
    return generation


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
    if item["freeze_artifact_id"] != FREEZE_ARTIFACT_ID:
        _fail("journal freeze artifact differs")
    if item["campaign_id"] != campaign_id:
        _fail("journal campaign_id differs")
    kind = _text(item["kind"], "kind")
    if kind not in RECORD_KINDS:
        _fail("journal kind differs")
    _nonneg_int(item["journal_sequence"], "journal_sequence")
    _digest_text(item["journal_parent_sha256"], "journal_parent_sha256")
    _nonneg_int(item["store_generation"], "store_generation")
    _digest_text(item["predecessor_checkpoint_sha256"], "predecessor_checkpoint_sha256")
    _nonclaim_flags(item)
    payload = item["payload"]
    if kind == "seed":
        body = _exact_mapping(payload, _SEED_PAYLOAD_KEY_SET, "seed payload")
        if tuple(live_member_keys(body["member_keys"])) != member_keys:
            _fail("seed payload member_keys differ")
        members = body["members"]
        if type(members) is not dict or set(members) != set(member_keys):
            _fail("seed payload members differ")
        for key in member_keys:
            members[key] = _exact_mapping(members[key], _MEMBER_STATE_KEY_SET, "member state")
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
    if item["freeze_artifact_id"] != FREEZE_ARTIFACT_ID:
        _fail("checkpoint freeze artifact differs")
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
    if tuple(live_member_keys(item["member_keys"])) != member_keys:
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
    _nonclaim_flags(item)
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
        raise PRO20EV1ProtocolError("generation documents are not canonical JSON") from error
    if type(journal) is not dict or type(checkpoint) is not dict or type(manifest) is not dict:
        _fail("generation documents are not JSON objects")
    campaign_id = _text(checkpoint.get("campaign_id"), "campaign_id")
    member_keys = tuple(live_member_keys(checkpoint.get("member_keys")))
    checkpoint = _validate_checkpoint_document(
        checkpoint, campaign_id=campaign_id, member_keys=member_keys
    )
    journal = _validate_journal_envelope(
        journal, campaign_id=campaign_id, member_keys=member_keys
    )
    manifest = _validate_generation_manifest(
        manifest,
        campaign_id=campaign_id,
    )
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
    derived_states = {key: _member_state(bundles[key]) for key in member_keys}
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
    leaves = {path: _digest(payload) for path, payload in files.items()}
    if manifest_mapping.get("leaves") != leaves:
        _fail("generation manifest leaves disagree with the file inventory")
    files["generation_manifest.json"] = _canonical(dict(manifest_mapping), "generation manifest")
    return files


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
    return {
        "schema": JOURNAL_SCHEMA,
        "schema_version": SCHEMA_VERSION,
        "protocol_artifact_id": PROTOCOL_ARTIFACT_ID,
        "freeze_artifact_id": FREEZE_ARTIFACT_ID,
        "campaign_id": campaign_id,
        "kind": kind,
        "journal_sequence": journal_sequence,
        "journal_parent_sha256": journal_parent_sha256,
        "store_generation": store_generation,
        "predecessor_checkpoint_sha256": predecessor_checkpoint_sha256,
        "payload": dict(payload),
        **{name: False for name in NONCLAIM_FLAG_NAMES},
    }


def production_cohort_complete(
    members: Mapping[str, Mapping[str, object]],
) -> bool:
    """Classify the exact six-member event boundary without evolving a member.

    This is the pure production wrapper around PROTO19's already-tested
    all-member predicate.  It validates the complete state-record shape and
    canonical member order; it does not manufacture a generation, terminal,
    authority, or accepted endpoint.
    """

    if not isinstance(members, Mapping) or tuple(members) != HLT17_MEMBER_KEYS:
        _fail("completion states must use the canonical six live members")
    checked = {
        key: _exact_mapping(
            dict(members[key]), _MEMBER_STATE_KEY_SET, f"{key} completion state"
        )
        for key in HLT17_MEMBER_KEYS
    }
    return bool(_cohort_complete(checked))


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
    complete = production_cohort_complete(members)
    if complete:
        if kind != "accepted_fine":
            _fail("first-event completion is not available for this kind")
        disposition = "first_event_complete"
        terminal = True
    elif kind in TERMINAL_KINDS:
        disposition = kind
        terminal = True
        complete = False
    else:
        disposition = "nonterminal"
        terminal = False
        complete = False
    return {
        "schema": CHECKPOINT_SCHEMA,
        "schema_version": SCHEMA_VERSION,
        "protocol_artifact_id": PROTOCOL_ARTIFACT_ID,
        "freeze_artifact_id": FREEZE_ARTIFACT_ID,
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
        **{name: False for name in NONCLAIM_FLAG_NAMES},
    }


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
        "freeze_artifact_id": FREEZE_ARTIFACT_ID,
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
        "calibration_claimed": False,
        "eligibility_claimed": False,
        "physical_classification": False,
    }
    return _validate_generation_manifest(mapping, campaign_id=campaign_id)


def _validate_generation_manifest(
    value: object,
    *,
    campaign_id: str,
) -> dict[str, object]:
    item = _exact_mapping(value, _MANIFEST_KEY_SET, "generation manifest")
    if (
        item["schema"] != GENERATION_MANIFEST_SCHEMA
        or item["schema_version"] != SCHEMA_VERSION
    ):
        _fail("generation manifest schema differs")
    if item["protocol_artifact_id"] != PROTOCOL_ARTIFACT_ID:
        _fail("generation manifest protocol differs")
    if item["freeze_artifact_id"] != FREEZE_ARTIFACT_ID:
        _fail("generation manifest freeze artifact differs")
    if item["campaign_id"] != _campaign_id(campaign_id):
        _fail("generation manifest campaign_id differs")
    item["store_generation"] = _nonneg_int(
        item["store_generation"], "manifest store_generation"
    )
    item["journal_sequence"] = _nonneg_int(
        item["journal_sequence"], "manifest journal_sequence"
    )
    for name in (
        "journal_record_sha256",
        "checkpoint_sha256",
        "predecessor_checkpoint_sha256",
    ):
        item[name] = _digest_text(item[name], f"manifest {name}")
    item["kind"] = _text(item["kind"], "generation manifest kind")
    if item["kind"] not in RECORD_KINDS:
        _fail("generation manifest kind differs")
    item["changed_member_key"] = _optional_member_key(
        item["changed_member_key"], HLT17_MEMBER_KEYS
    )
    for name in (
        "production_write_authorized",
        "campaign_execution_authorized",
        "calibration_claimed",
        "eligibility_claimed",
        "physical_classification",
    ):
        _bool_false(item[name], f"generation manifest {name}")
    _canonical(item, "generation manifest")
    return item


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


def build_seed_generation(
    *,
    root: Mapping[str, object],
    bundles: Mapping[str, HLT17GenerationBundle],
) -> PRO20EV1ValidatedGeneration:
    """Validate a complete live cohort seed before any directory is visible."""

    root = validate_root_manifest(dict(root))
    member_keys = tuple(root["member_keys"])
    if member_keys != HLT17_MEMBER_KEYS:
        _fail("actual seed bundles must use the six live production members")
    if set(bundles) != set(member_keys):
        _fail("seed bundles do not match the cohort")
    rebuilt = {key: _revalidate(bundles[key]) for key in member_keys}
    states = {key: _member_state(rebuilt[key]) for key in member_keys}
    for key in member_keys:
        _require_live_seed_bundle(rebuilt[key], key)
        if states[key]["source_current"] != 0 or states[key]["cfl_current"] != 0:
            _fail("seed per-macro retry counters must start at zero")
    payload = {
        "member_keys": list(member_keys),
        "members": {key: dict(states[key]) for key in member_keys},
    }
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
    predecessor: PRO20EV1ValidatedGeneration,
    predecessor_bundles: Mapping[str, HLT17GenerationBundle],
    kind: str,
    member_key: str,
    successor_bundle: HLT17GenerationBundle,
    terminal_evidence: Mapping[str, object] | None = None,
) -> PRO20EV1ValidatedGeneration:
    """Validate one live member attempt in canonical order before write."""

    root = validate_root_manifest(dict(root))
    member_keys = tuple(root["member_keys"])
    if kind not in RECORD_KINDS or kind == "seed":
        _fail("successor kind differs")
    predecessor = _rebind_generation(predecessor)
    closed = _generation_forbids_append(predecessor)
    if closed is not None:
        _fail(closed)
    required = next_required_member_key(predecessor)
    if member_key != required:
        _fail("changed member is not the canonical next live member")
    if member_key not in member_keys:
        _fail("changed member is outside the cohort")
    if set(predecessor_bundles) != set(member_keys):
        _fail("predecessor bundles do not match the cohort")
    prior = {key: _revalidate(predecessor_bundles[key]) for key in member_keys}
    for key in member_keys:
        expected = predecessor.bundles[key]
        if (
            prior[key].descriptor != expected.descriptor
            or prior[key].payload != expected.payload
            or prior[key].cursor_bytes != expected.cursor_bytes
        ):
            _fail("predecessor bundle bytes disagree with the authenticated generation")
    successor = _revalidate(successor_bundle)
    rebuilt = dict(prior)
    rebuilt[member_key] = successor
    prior_states = {key: _member_state(prior[key]) for key in member_keys}
    new_states = {key: _member_state(rebuilt[key]) for key in member_keys}
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
    try:
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
    except HLT17ProtocolError as error:
        _wrap_hlt17(error)
        raise AssertionError("unreachable")
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
    predecessor: PRO20EV1ValidatedGeneration | None,
    predecessor_bundles: Mapping[str, HLT17GenerationBundle] | None,
) -> PRO20EV1ValidatedGeneration:
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
    if checkpoint_sha256.encode("ascii") in journal_raw:
        _fail("journal record contains the successor checkpoint hash")
    if "checkpoint_sha256" in journal:
        _fail("journal record contains a checkpoint hash")
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
    generation = PRO20EV1ValidatedGeneration(
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
    predecessor: PRO20EV1ValidatedGeneration | None,
    predecessor_bundles: Mapping[str, HLT17GenerationBundle] | None,
    expected: PRO20EV1ValidatedGeneration | None = None,
) -> PRO20EV1ValidatedGeneration:
    """Reconstruct one production generation from immutable file bytes."""

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
        raise PRO20EV1ProtocolError("generation documents are not canonical JSON") from error
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
    if manifest["store_generation"] != checkpoint["store_generation"]:
        _fail("manifest store generation disagrees with checkpoint")
    if manifest["journal_sequence"] != checkpoint["journal_sequence"]:
        _fail("manifest journal sequence disagrees with checkpoint")
    if (
        manifest["predecessor_checkpoint_sha256"]
        != checkpoint["predecessor_checkpoint_sha256"]
    ):
        _fail("manifest checkpoint parent disagrees with checkpoint")
    if manifest["changed_member_key"] != checkpoint["changed_member_key"]:
        _fail("manifest changed member disagrees with checkpoint")
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
    bundles: dict[str, HLT17GenerationBundle] = {}
    if kind == "seed":
        if predecessor is not None:
            _fail("seed cannot have a predecessor generation")
        if journal["journal_parent_sha256"] != GENESIS_DIGEST:
            _fail("seed journal parent differs")
        members_payload = journal["payload"]["members"]
        for key in member_keys:
            state = checkpoint["members"][key]
            bundle = _from_bytes(
                new_blobs[str(state["descriptor_sha256"])],
                new_blobs[str(state["payload_sha256"])],
                new_blobs[str(state["cursor_sha256"])],
            )
            derived = _member_state(bundle)
            if derived != state or derived != members_payload[key]:
                _fail("seed member labels disagree with immutable bytes")
            _require_live_seed_bundle(bundle, key)
            bundles[key] = bundle
    else:
        if predecessor is None or predecessor_bundles is None:
            _fail("successor generation omitted its predecessor")
        predecessor = _rebind_generation(predecessor)
        closed = _generation_forbids_append(predecessor)
        if closed is not None:
            _fail(closed)
        required = next_required_member_key(predecessor)
        if str(changed) != required:
            _fail("changed member is not the canonical next live member")
        if journal["journal_parent_sha256"] != predecessor.journal_sha256:
            _fail("journal parent is not the unique predecessor journal")
        if checkpoint["predecessor_checkpoint_sha256"] != predecessor.checkpoint_sha256:
            _fail("checkpoint parent is not the unique predecessor checkpoint")
        if journal["store_generation"] != predecessor.store_generation + 1:
            _fail("store generation hole or fork")
        if journal["journal_sequence"] != predecessor.journal_sequence + 1:
            _fail("journal sequence hole or fork")
        prior = {key: _revalidate(predecessor_bundles[key]) for key in member_keys}
        for key in member_keys:
            previous_state = _member_state(prior[key])
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
            bundle = _from_bytes(descriptor, payload_bytes, cursor_bytes)
            derived = _member_state(bundle)
            if derived != listed:
                _fail("changed member labels disagree with immutable bytes")
            bundles[key] = bundle
        payload = journal["payload"]
        previous = _member_state(prior[str(changed)])
        current = _member_state(bundles[str(changed)])
        try:
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
        except HLT17ProtocolError as error:
            _wrap_hlt17(error)
            raise AssertionError("unreachable")
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
    derived_states = {key: _member_state(bundles[key]) for key in member_keys}
    if derived_states != checkpoint["members"]:
        _fail("checkpoint member states disagree with immutable bundle bytes")
    derived_complete = _cohort_complete(derived_states)
    if derived_complete != bool(checkpoint["first_event_complete"]):
        _fail("first-event completion disagrees with the reconstructed cohort")
    if derived_complete and (
        kind != "accepted_fine" or checkpoint["disposition"] != "first_event_complete"
    ):
        _fail("first-event completion must be the last accepted fine")
    generation = PRO20EV1ValidatedGeneration(
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
    if expected is not None and generation.files != expected.files:
        _fail("assembled generation files are not redundant")
    return generation


def authenticate_generation_chain(
    *,
    root: Mapping[str, object],
    generations: Sequence[Mapping[str, bytes]],
) -> tuple[PRO20EV1ValidatedGeneration, ...]:
    """Authenticate a unique complete live chain with no holes or forks."""

    root = validate_root_manifest(dict(root))
    if not generations:
        _fail("campaign store has no published generations")
    authenticated: list[PRO20EV1ValidatedGeneration] = []
    predecessor: PRO20EV1ValidatedGeneration | None = None
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
    "CHECKPOINT_KEYS",
    "CHECKPOINT_SCHEMA",
    "COMMON_EVENT_INDEX",
    "EVENT_ORIGIN",
    "EVENT_ORIGIN_HEX",
    "EVENT_TARGET",
    "EVENT_TARGET_HEX",
    "EXACT_DELTA_AUTHORITY_SEAM",
    "FREEZE_ARTIFACT_ID",
    "GENERATION_MANIFEST_SCHEMA",
    "GENESIS_DIGEST",
    "INDEPENDENT_BINDER_SEAM",
    "JOURNAL_SCHEMA",
    "NONCLAIM_FLAG_NAMES",
    "PRODUCTION_NAMESPACE",
    "PROTOCOL_ARTIFACT_ID",
    "PRO20EV1AuthorityReceipt",
    "PRO20EV1ProtocolError",
    "PRO20EV1ValidatedGeneration",
    "RECEIPT_KEYS",
    "RECEIPT_SCHEMA",
    "RESOURCE_GATE_SEAM",
    "ROOT_MANIFEST_KEYS",
    "ROOT_SCHEMA",
    "RUNNER_SEAM",
    "STORE_KIND_PRODUCTION",
    "TERMINAL_LOCK_SCHEMA",
    "authenticate_generation_chain",
    "authority_receipt_from_root",
    "build_authority_receipt",
    "build_root_manifest",
    "build_seed_generation",
    "build_successor_generation",
    "encode_session_close",
    "encode_session_open",
    "encode_terminal_lock",
    "generation_directory_name",
    "live_member_keys",
    "next_required_member_key",
    "parse_generation_directory_name",
    "parse_session_name",
    "parse_validated_generation",
    "production_cohort_complete",
    "session_close_name",
    "session_open_name",
    "validate_authority_receipt",
    "validate_root_manifest",
    "validate_session_close",
    "validate_session_open",
    "validate_terminal_lock",
]
