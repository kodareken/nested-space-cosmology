"""Independent, outcome-neutral PRO20 first-event evidence binder core.

This module deliberately does not import the PRO20 runner, authority,
protocol, or store modules.  It reads a completed production namespace
through no-follow descriptors, authenticates its immutable byte grammar, and
derives transition classes from accepted bundle and cursor bytes.  Runner
``kind``/``disposition`` fields are compared only after that derivation.

The tracked PREF1 config/result and the expensive physical replay are added
only after the one authorized event has a valid terminal.  Keeping this core
outcome-neutral lets its hostile-store and label-independence controls be
tested before any result is known.
"""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
import os
from pathlib import Path, PurePosixPath
import re
import stat
from types import MappingProxyType
from typing import Any, Mapping, Sequence

from recursive_horizons.evidence_io import (
    CanonicalJSONError,
    DEFAULT_MAX_BYTES,
    canonical_json_bytes,
    load_canonical_json,
)

from .hlt17_admission_runtime import (
    IMPLEMENTATION_IDENTITY_KEYS,
    hlt17_implementation_identity,
)
from .hlt17_imp1_cursor import (
    CFL_OWNER,
    FRESH_READY,
    HLT17_CFL_RETRY_CAP,
    HLT17_MEMBER_KEYS,
    HLT17_SOURCE_RETRY_CAP,
    IDENTITY_SCOPE_LIVE,
    RETRY_PENDING,
    SOURCE_OWNER,
    TEMPORAL_OWNER,
    cursor_sha256,
    hlt17_member_spec,
    ledger_sha256,
    restore_hlt17_cursor,
)
from .hlt17_member_checkpoint import HLT17GenerationBundle
from .hlt17_member_codec import ORIGIN_PROTOCOL, RUNTIME_PROTOCOL
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


ARTIFACT_ID = "FGC-1-PRO20-EV1-PREF1"
PROTOCOL_ARTIFACT_ID = "FGC-2-SF1-PROTO19"
FREEZE_ARTIFACT_ID = "FGC-1-PRO20-EV1-FRZ1"
PRODUCTION_NAMESPACE = "runs/fgc-2-sf1/pro20-event1/calibration"
PRODUCTION_CONTAINER = "runs/fgc-2-sf1/pro20-event1"
SCHEMA_VERSION = 1
STORE_KIND = "production_event1"
EVENT_ORIGIN_HEX = (23.0 / 16.0).hex()
EVENT_TARGET_HEX = (3.0 / 2.0).hex()
COMMON_EVENT_INDEX = 23
GENESIS_DIGEST = "0" * 64
LIVE_AMPLITUDE = "3"
ACCEPTED_SUBSTEPS = IMP1_ACCEPTED_SUBSTEP_COUNT
TEMPORAL_EXHAUSTION_REASONS = frozenset(
    {"coordinate_lattice", "maximum_temporal_retries", "minimum_macro_step"}
)
CONTAINER_STAGING_PREFIX = ".pro20-event1.evidence-io-stage-"
RUNNER_SEAM = (
    "a committed live runner is an external owner; this slice does not execute "
    "or treat runner labels as evidence"
)
EXACT_DELTA_AUTHORITY_SEAM = (
    "exact-delta freeze authority is an external owner; this receipt type is "
    "not that freeze and construction is not live authorization"
)
RESOURCE_GATE_SEAM = (
    "largest-grid and physical-RHS resource qualification is an external gate; "
    "this slice does not execute or authorize a production preparation"
)
INDEPENDENT_BINDER_SEAM = (
    "an independent PREF1 binder must reconstruct evidence without trusting "
    "runner labels; this slice is not that binder"
)
SEAM_VALUES = {
    "runner_seam": RUNNER_SEAM,
    "exact_delta_authority_seam": EXACT_DELTA_AUTHORITY_SEAM,
    "resource_gate_seam": RESOURCE_GATE_SEAM,
    "independent_binder_seam": INDEPENDENT_BINDER_SEAM,
}
ROOT_IDENTITY_HASH_NAMES = (
    "authority_sha256",
    "implementation_sha256",
    "config_sha256",
    "source_sha256",
    "origin_sha256",
    "environment_sha256",
)

ROOT_SCHEMA = "FGC-1-PRO20-EV1-campaign-root-v1"
RECEIPT_SCHEMA = "FGC-1-PRO20-EV1-authority-receipt-v1"
JOURNAL_SCHEMA = "FGC-1-PRO20-EV1-campaign-journal-v1"
CHECKPOINT_SCHEMA = "FGC-1-PRO20-EV1-campaign-checkpoint-v1"
MANIFEST_SCHEMA = "FGC-1-PRO20-EV1-generation-manifest-v1"
SESSION_OPEN_SCHEMA = "FGC-1-PRO20-EV1-writer-session-open-v1"
SESSION_CLOSE_SCHEMA = "FGC-1-PRO20-EV1-writer-session-close-v1"
TERMINAL_LOCK_SCHEMA = "FGC-1-PRO20-EV1-terminal-lock-v1"
WRITER_EXCLUSION_SCHEMA = "FGC-1-PRO20-EV1-writer-exclusion-v1"

ROOT_MANIFEST_NAME = "root_manifest.json"
WRITER_EXCLUSION_NAME = "writer_exclusion.lock"
TERMINAL_LOCK_NAME = "terminal.lock"

RECORD_KINDS = (
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
TERMINAL_KINDS = frozenset(
    {
        "source_retry_exhausted",
        "cfl_retry_exhausted",
        "temporal_retry_exhausted",
        "invalid_premise",
        "resource_exhausted",
    }
)
BLOB_UNCHANGED_KINDS = frozenset(
    {"temporal_retry_exhausted", "invalid_premise", "resource_exhausted"}
)
CURSOR_RETRY_KINDS = frozenset(
    {
        "source_retry",
        "cfl_retry",
        "temporal_retry",
        "source_retry_exhausted",
        "cfl_retry_exhausted",
    }
)

NONCLAIM_NAMES = (
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

MEMBER_STATE_KEYS = (
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

ROOT_KEYS = (
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
    *NONCLAIM_NAMES,
    "runner_seam",
    "exact_delta_authority_seam",
    "resource_gate_seam",
    "independent_binder_seam",
)

RECEIPT_KEYS = (
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
    *NONCLAIM_NAMES,
    "runner_seam",
    "exact_delta_authority_seam",
    "resource_gate_seam",
    "independent_binder_seam",
)

JOURNAL_KEYS = (
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
    *NONCLAIM_NAMES,
)

CHECKPOINT_KEYS = (
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
    *NONCLAIM_NAMES,
)

MANIFEST_KEYS = (
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

SESSION_OPEN_KEYS = (
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

SESSION_CLOSE_KEYS = (
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

TERMINAL_LOCK_KEYS = (
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

ATTEMPT_PAYLOAD_KEYS = (
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

_SHA = re.compile(r"\A[0-9a-f]{64}\Z")
_GENERATION = re.compile(r"\Ageneration-([0-9]{16})\Z")
_SESSION = re.compile(r"\Asession-([0-9]{16})-(open|close)\.json\Z")
_MAX_LEAF_BYTES = DEFAULT_MAX_BYTES
_MAX_TOTAL_BYTES = 16 * 1024 * 1024 * 1024
_MAX_LEAVES = 32768
_MAX_GENERATIONS = 4097
_DIRECTORY_FLAGS = os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) | getattr(os, "O_NOFOLLOW", 0)
_LEAF_FLAGS = os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0)


class PRO20EV1PREF1Error(ValueError):
    """The raw store or independent transition reduction differs."""

    def __init__(self, stop_id: str, detail: str) -> None:
        super().__init__(f"{stop_id}: {detail}")
        self.stop_id = stop_id
        self.detail = detail


def _stop(stop_id: str, detail: str) -> None:
    raise PRO20EV1PREF1Error(stop_id, detail)


def _digest(raw: bytes) -> str:
    if type(raw) is not bytes:
        _stop("PREF1_BYTES", "payload is not immutable bytes")
    return sha256(raw).hexdigest()


def _sha(value: object, label: str) -> str:
    if type(value) is not str or _SHA.fullmatch(value) is None:
        _stop("PREF1_SCHEMA", f"{label} is not a SHA-256 digest")
    return value


def _exact_keys(value: object, keys: Sequence[str], label: str) -> dict[str, Any]:
    if type(value) is not dict or set(value) != set(keys):
        _stop("PREF1_SCHEMA", f"{label} fields differ")
    return dict(value)


def _json(raw: bytes, label: str) -> dict[str, Any]:
    try:
        value = load_canonical_json(raw)
    except (CanonicalJSONError, TypeError, ValueError) as error:
        raise PRO20EV1PREF1Error(
            "PREF1_CANONICAL", f"{label} is not canonical JSON"
        ) from error
    if type(value) is not dict:
        _stop("PREF1_SCHEMA", f"{label} must be a JSON object")
    return value


def _nonnegative_int(value: object, label: str) -> int:
    if type(value) is not int or value < 0:
        _stop("PREF1_SCHEMA", f"{label} must be a nonnegative integer")
    return value


def _false(value: object, label: str) -> None:
    if value is not False:
        _stop("PREF1_PROMOTION", f"{label} must remain false")


def _nonclaims(mapping: Mapping[str, object], *, names: Sequence[str] = NONCLAIM_NAMES) -> None:
    for name in names:
        _false(mapping.get(name), name)


def _text(value: object, label: str) -> str:
    if type(value) is not str or not value:
        _stop("PREF1_SCHEMA", f"{label} must be nonempty text")
    return value


def _campaign_id(value: object) -> str:
    campaign_id = _text(value, "campaign_id")
    if "/" in campaign_id or "\\" in campaign_id or ".." in campaign_id:
        _stop("PREF1_ROOT", "campaign_id is not a safe identifier")
    return campaign_id


def _seams(mapping: Mapping[str, object], *, label: str) -> None:
    for name, expected in SEAM_VALUES.items():
        if mapping.get(name) != expected:
            _stop("PREF1_ROOT", f"{label}.{name} differs")


def _identity_hashes(mapping: Mapping[str, object]) -> None:
    identities: list[str] = []
    for name in ROOT_IDENTITY_HASH_NAMES:
        identities.append(_sha(mapping.get(name), f"root.{name}"))
    if len(set(identities)) != len(ROOT_IDENTITY_HASH_NAMES):
        _stop(
            "PREF1_ROOT",
            "authority/implementation/config/source/origin/environment identities must be distinct",
        )


def _directory_names(root: Path, relative: str) -> tuple[str, ...]:
    parts = _path_parts(relative) if relative else ()
    descriptor = _open_directory_chain(root, parts)
    try:
        with os.scandir(descriptor) as scan:
            return tuple(sorted(entry.name for entry in scan))
    finally:
        os.close(descriptor)


def _reject_container_staging(root: Path) -> None:
    parent = str(PurePosixPath(*_path_parts(PRODUCTION_CONTAINER)[:-1]))
    names = _directory_names(root, parent)
    if any(name.startswith(CONTAINER_STAGING_PREFIX) for name in names):
        _stop("PREF1_TREE", "publication staging sibling is present")
    container = _path_parts(PRODUCTION_CONTAINER)[-1]
    if container not in names:
        _stop("PREF1_TREE", "production container is absent")


def _require_calibration_container(root: Path, store_relative: str) -> None:
    parts = _path_parts(store_relative)
    parent = str(PurePosixPath(*parts[:-1]))
    names = _directory_names(root, parent)
    final = parts[-1]
    if any(name.startswith(f".{final}.") for name in names):
        _stop("PREF1_TREE", "publication staging sibling is present")
    if names != (final,):
        _stop("PREF1_TREE", "production container has a foreign sibling")


def _spatial_order(method_label: str) -> int:
    if method_label == "RK4":
        return 4
    if method_label == "SSPRK3":
        return 2
    _stop("PREF1_SEED", "live spatial operator identity differs")
    raise AssertionError("unreachable")


def _stage_records(member_key: str) -> int:
    spec = hlt17_member_spec(member_key)
    return IMP1_METHOD_STAGE_RECORD_COUNT[spec.method]


def _overlay_exhausted(cursor: object) -> bool:
    overlay = getattr(cursor, "overlay", None)
    if overlay is None:
        return False
    if overlay.exhausted is True or overlay.next_plan is None:
        return True
    cap = HLT17_SOURCE_RETRY_CAP if overlay.owner == SOURCE_OWNER else HLT17_CFL_RETRY_CAP
    return overlay.owner_current_count > cap


@dataclass(frozen=True, slots=True)
class RawLeaf:
    relative: str
    raw: bytes
    sha256: str
    byte_count: int
    identity: tuple[int, int, int, int, int, int]


@dataclass(frozen=True, slots=True)
class StoreSnapshot:
    repository: Path
    store_relative: str
    leaves: Mapping[str, RawLeaf]
    directories: tuple[str, ...]
    lexical_sha256: str
    byte_count: int

    def __post_init__(self) -> None:
        object.__setattr__(self, "leaves", MappingProxyType(dict(self.leaves)))


@dataclass(frozen=True, slots=True)
class IndependentGeneration:
    index: int
    kind: str
    derived_kind: str
    changed_member_key: str | None
    checkpoint_sha256: str
    journal_sha256: str
    disposition: str
    first_event_complete: bool
    terminal: bool
    member_states: Mapping[str, Mapping[str, object]]
    bundles: Mapping[str, HLT17GenerationBundle]

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "member_states",
            MappingProxyType(
                {key: MappingProxyType(dict(value)) for key, value in self.member_states.items()}
            ),
        )
        object.__setattr__(self, "bundles", MappingProxyType(dict(self.bundles)))


@dataclass(frozen=True, slots=True)
class IndependentTerminal:
    campaign_id: str
    generations: tuple[IndependentGeneration, ...]
    disposition: str
    first_event_complete: bool
    terminal_kind: str
    checkpoint_sha256: str
    journal_sha256: str
    terminal_lock_sha256: str
    session_count: int
    store_leaf_count: int
    store_byte_count: int
    store_lexical_sha256: str
    runner_label_sufficient: bool = False
    calibration_eligible: bool = False

    def __post_init__(self) -> None:
        if not self.generations:
            _stop("PREF1_LINEAGE", "terminal has no generations")
        if self.runner_label_sufficient or self.calibration_eligible:
            _stop("PREF1_PROMOTION", "terminal cannot trust labels or open calibration")


def _path_parts(relative: str) -> tuple[str, ...]:
    path = PurePosixPath(relative)
    if path.is_absolute() or not path.parts or any(part in {"", ".", ".."} for part in path.parts):
        _stop("PREF1_PATH", "repository-relative path is unsafe")
    return tuple(path.parts)


def _identity(info: os.stat_result) -> tuple[int, int, int, int, int, int]:
    return (
        info.st_dev,
        info.st_ino,
        info.st_size,
        info.st_mtime_ns,
        info.st_ctime_ns,
        info.st_nlink,
    )


def _open_directory_chain(root: Path, parts: Sequence[str]) -> int:
    try:
        root_info = root.lstat()
    except OSError as error:
        raise PRO20EV1PREF1Error("PREF1_PATH", "repository root is absent") from error
    if stat.S_ISLNK(root_info.st_mode) or not stat.S_ISDIR(root_info.st_mode):
        _stop("PREF1_PATH", "repository root is unsafe")
    current = os.open(root, _DIRECTORY_FLAGS)
    try:
        opened = os.fstat(current)
        if (root_info.st_dev, root_info.st_ino) != (opened.st_dev, opened.st_ino):
            _stop("PREF1_PATH", "repository root raced")
        for component in parts:
            before = os.stat(component, dir_fd=current, follow_symlinks=False)
            if stat.S_ISLNK(before.st_mode) or not stat.S_ISDIR(before.st_mode):
                _stop("PREF1_PATH", f"unsafe directory component: {component}")
            child = os.open(component, _DIRECTORY_FLAGS, dir_fd=current)
            after = os.fstat(child)
            if (before.st_dev, before.st_ino) != (after.st_dev, after.st_ino):
                os.close(child)
                _stop("PREF1_PATH", f"directory component raced: {component}")
            os.close(current)
            current = child
        return current
    except BaseException:
        os.close(current)
        raise


def _read_leaf(parent: int, name: str, relative: str) -> RawLeaf:
    before = os.stat(name, dir_fd=parent, follow_symlinks=False)
    if (
        stat.S_ISLNK(before.st_mode)
        or not stat.S_ISREG(before.st_mode)
        or before.st_nlink != 1
        or before.st_size < 0
        or before.st_size > _MAX_LEAF_BYTES
    ):
        _stop("PREF1_PATH", f"unsafe or oversized leaf: {relative}")
    descriptor = os.open(name, _LEAF_FLAGS, dir_fd=parent)
    try:
        opened = os.fstat(descriptor)
        expected = _identity(before)
        if expected != _identity(opened):
            _stop("PREF1_PATH", f"leaf raced before read: {relative}")
        remaining = before.st_size
        chunks: list[bytes] = []
        while remaining:
            block = os.read(descriptor, min(1024 * 1024, remaining))
            if not block:
                _stop("PREF1_PATH", f"leaf ended early: {relative}")
            chunks.append(block)
            remaining -= len(block)
        if os.read(descriptor, 1):
            _stop("PREF1_PATH", f"leaf grew during read: {relative}")
        if expected != _identity(os.fstat(descriptor)):
            _stop("PREF1_PATH", f"leaf changed during read: {relative}")
        if expected != _identity(os.stat(name, dir_fd=parent, follow_symlinks=False)):
            _stop("PREF1_PATH", f"leaf path changed during read: {relative}")
        raw = b"".join(chunks)
        return RawLeaf(relative, raw, _digest(raw), len(raw), expected)
    finally:
        os.close(descriptor)


def _scan_fd(
    directory_fd: int,
    *,
    prefix: str,
    leaves: dict[str, RawLeaf],
    directories: list[str],
    byte_total: list[int],
    depth: int,
) -> None:
    if depth > 2:
        _stop("PREF1_TREE", "store nesting exceeds the generation/blob grammar")
    with os.scandir(directory_fd) as entries:
        ordered = sorted(entries, key=lambda item: item.name)
    for entry in ordered:
        name = entry.name
        if name in {"", ".", ".."} or "/" in name or "\x00" in name:
            _stop("PREF1_PATH", "unsafe directory entry")
        relative = name if not prefix else f"{prefix}/{name}"
        info = os.stat(name, dir_fd=directory_fd, follow_symlinks=False)
        if stat.S_ISLNK(info.st_mode):
            _stop("PREF1_PATH", f"symlink in store: {relative}")
        if stat.S_ISREG(info.st_mode):
            if len(leaves) >= _MAX_LEAVES:
                _stop("PREF1_RESOURCE", "store leaf-count ceiling exceeded")
            leaf = _read_leaf(directory_fd, name, relative)
            byte_total[0] += leaf.byte_count
            if byte_total[0] > _MAX_TOTAL_BYTES:
                _stop("PREF1_RESOURCE", "store byte ceiling exceeded")
            leaves[relative] = leaf
            continue
        if not stat.S_ISDIR(info.st_mode):
            _stop("PREF1_PATH", f"special object in store: {relative}")
        child = os.open(name, _DIRECTORY_FLAGS, dir_fd=directory_fd)
        try:
            opened = os.fstat(child)
            if (info.st_dev, info.st_ino) != (opened.st_dev, opened.st_ino):
                _stop("PREF1_PATH", f"directory raced: {relative}")
            directories.append(relative)
            if len(directories) > 2 * _MAX_GENERATIONS:
                _stop("PREF1_RESOURCE", "store directory-count ceiling exceeded")
            _scan_fd(
                child,
                prefix=relative,
                leaves=leaves,
                directories=directories,
                byte_total=byte_total,
                depth=depth + 1,
            )
            after = os.stat(name, dir_fd=directory_fd, follow_symlinks=False)
            if (info.st_dev, info.st_ino) != (after.st_dev, after.st_ino):
                _stop("PREF1_PATH", f"directory changed: {relative}")
        finally:
            os.close(child)


def _scan_once(repository: Path, store_relative: str) -> tuple[dict[str, RawLeaf], tuple[str, ...]]:
    descriptor = _open_directory_chain(repository, _path_parts(store_relative))
    try:
        leaves: dict[str, RawLeaf] = {}
        directories: list[str] = []
        byte_total = [0]
        _scan_fd(
            descriptor,
            prefix="",
            leaves=leaves,
            directories=directories,
            byte_total=byte_total,
            depth=0,
        )
        return leaves, tuple(sorted(directories))
    finally:
        os.close(descriptor)


def _lexical_digest(leaves: Mapping[str, RawLeaf]) -> str:
    digestor = sha256()
    for relative in sorted(leaves):
        leaf = leaves[relative]
        digestor.update(relative.encode("utf-8"))
        digestor.update(b"\0")
        digestor.update(str(leaf.byte_count).encode("ascii"))
        digestor.update(b"\0")
        digestor.update(leaf.sha256.encode("ascii"))
        digestor.update(b"\n")
    return digestor.hexdigest()


def snapshot_store(
    repository: Path,
    *,
    store_relative: str = PRODUCTION_NAMESPACE,
) -> StoreSnapshot:
    """Read one stable complete store snapshot without importing store code."""

    root = Path(repository)
    if not root.is_absolute():
        _stop("PREF1_PATH", "repository path must be absolute")
    try:
        resolved = root.resolve(strict=True)
    except OSError as error:
        raise PRO20EV1PREF1Error("PREF1_PATH", "repository cannot be resolved") from error
    if resolved != root:
        _stop("PREF1_PATH", "repository path must already be canonical")
    _reject_container_staging(root)
    _require_calibration_container(root, store_relative)
    first = _scan_once(root, store_relative)
    second = _scan_once(root, store_relative)
    if first != second:
        _stop("PREF1_PATH", "store changed during independent snapshot")
    leaves, directories = first
    return StoreSnapshot(
        repository=root,
        store_relative=store_relative,
        leaves=leaves,
        directories=directories,
        lexical_sha256=_lexical_digest(leaves),
        byte_count=sum(item.byte_count for item in leaves.values()),
    )


def _schema_header(mapping: Mapping[str, object], *, schema: str, label: str) -> None:
    if (
        mapping.get("schema") != schema
        or mapping.get("schema_version") != SCHEMA_VERSION
        or mapping.get("protocol_artifact_id") != PROTOCOL_ARTIFACT_ID
        or mapping.get("freeze_artifact_id") != FREEZE_ARTIFACT_ID
    ):
        _stop("PREF1_SCHEMA", f"{label} identity differs")


def _root_manifest(raw: bytes) -> dict[str, Any]:
    item = _exact_keys(_json(raw, "root manifest"), ROOT_KEYS, "root manifest")
    _schema_header(item, schema=ROOT_SCHEMA, label="root manifest")
    if (
        item["store_kind"] != STORE_KIND
        or item["namespace"] != PRODUCTION_NAMESPACE
        or tuple(item["member_keys"]) != tuple(HLT17_MEMBER_KEYS)
        or item["event_origin_time_hex"] != EVENT_ORIGIN_HEX
        or item["event_target_time_hex"] != EVENT_TARGET_HEX
        or item["common_event_index"] != COMMON_EVENT_INDEX
    ):
        _stop("PREF1_ROOT", "production root contract differs")
    item["campaign_id"] = _campaign_id(item["campaign_id"])
    _sha(item["authority_receipt_sha256"], "root.authority_receipt_sha256")
    _identity_hashes(item)
    identity = item["implementation_identity"]
    if type(identity) is not dict or identity != hlt17_implementation_identity():
        _stop("PREF1_ROOT", "implementation identity differs")
    _nonclaims(item)
    _seams(item, label="root")
    receipt = {
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
        "implementation_identity": dict(identity),
        **{name: item[name] for name in NONCLAIM_NAMES},
        "runner_seam": item["runner_seam"],
        "exact_delta_authority_seam": item["exact_delta_authority_seam"],
        "resource_gate_seam": item["resource_gate_seam"],
        "independent_binder_seam": item["independent_binder_seam"],
    }
    _exact_keys(receipt, RECEIPT_KEYS, "authority receipt")
    if _digest(canonical_json_bytes(receipt)) != item["authority_receipt_sha256"]:
        _stop("PREF1_ROOT", "authority receipt digest differs")
    return item


def _writer_exclusion(raw: bytes, root: Mapping[str, object]) -> dict[str, Any]:
    keys = (
        "schema",
        "schema_version",
        "protocol_artifact_id",
        "freeze_artifact_id",
        "role",
        "namespace",
        "authority_receipt_sha256",
        "production_write_authorized",
        "campaign_execution_authorized",
        "automatic_stale_writer_recovery",
    )
    item = _exact_keys(_json(raw, "writer exclusion"), keys, "writer exclusion")
    _schema_header(item, schema=WRITER_EXCLUSION_SCHEMA, label="writer exclusion")
    if (
        item["role"] != "process_thread_exclusion"
        or item["namespace"] != PRODUCTION_NAMESPACE
        or item["authority_receipt_sha256"] != root["authority_receipt_sha256"]
    ):
        _stop("PREF1_SESSION", "writer exclusion identity differs")
    for name in (
        "production_write_authorized",
        "campaign_execution_authorized",
        "automatic_stale_writer_recovery",
    ):
        _false(item[name], f"writer_exclusion.{name}")
    return item


def _generation_index(name: str) -> int:
    match = _GENERATION.fullmatch(name)
    if match is None:
        _stop("PREF1_TREE", f"generation directory name differs: {name}")
    value = int(match.group(1))
    if name != f"generation-{value:016d}":
        _stop("PREF1_TREE", "generation directory is noncanonical")
    return value


def _generation_names(snapshot: StoreSnapshot) -> tuple[str, ...]:
    roots = sorted(path for path in snapshot.directories if "/" not in path)
    names = tuple(path for path in roots if path.startswith("generation-"))
    if len(names) > _MAX_GENERATIONS:
        _stop("PREF1_RESOURCE", "generation ceiling exceeded")
    indices = [_generation_index(name) for name in names]
    if indices != list(range(len(indices))):
        _stop("PREF1_LINEAGE", "generation sequence has a gap or fork")
    required_directories = set(names)
    optional_blob_directories = {f"{name}/blobs" for name in names}
    observed = set(snapshot.directories)
    if not required_directories.issubset(observed) or not observed.issubset(
        required_directories | optional_blob_directories
    ):
        _stop("PREF1_TREE", "store directory grammar differs")
    return names


def _root_leaf_grammar(snapshot: StoreSnapshot, generation_names: Sequence[str]) -> None:
    generation_prefixes = tuple(f"{name}/" for name in generation_names)
    root_names = {
        relative
        for relative in snapshot.leaves
        if not relative.startswith(generation_prefixes)
    }
    required = {ROOT_MANIFEST_NAME, WRITER_EXCLUSION_NAME, TERMINAL_LOCK_NAME}
    if not required.issubset(root_names):
        _stop("PREF1_TREE", "completed store lacks root/lock evidence")
    for name in root_names - required:
        if _SESSION.fullmatch(name) is None:
            _stop("PREF1_TREE", f"foreign root leaf: {name}")


def _generation_leaf_map(snapshot: StoreSnapshot, name: str) -> dict[str, RawLeaf]:
    prefix = f"{name}/"
    answer = {
        relative[len(prefix) :]: leaf
        for relative, leaf in snapshot.leaves.items()
        if relative.startswith(prefix)
    }
    required = {"journal.json", "checkpoint.json", "generation_manifest.json"}
    if not required.issubset(answer):
        _stop("PREF1_TREE", f"{name} lacks required leaves")
    for relative in answer:
        if relative in required:
            continue
        parts = PurePosixPath(relative).parts
        if len(parts) != 2 or parts[0] != "blobs" or _SHA.fullmatch(parts[1]) is None:
            _stop("PREF1_TREE", f"foreign generation leaf: {name}/{relative}")
    return answer


def _session_index(name: str) -> tuple[int, str]:
    match = _SESSION.fullmatch(name)
    if match is None:
        _stop("PREF1_SESSION", "session filename differs")
    value = int(match.group(1))
    if name != f"session-{value:016d}-{match.group(2)}.json":
        _stop("PREF1_SESSION", "session filename is noncanonical")
    return value, match.group(2)


def _closed_sessions(
    snapshot: StoreSnapshot,
    root: Mapping[str, object],
    exclusion_raw: bytes,
) -> int:
    opened: dict[int, tuple[RawLeaf, dict[str, Any]]] = {}
    closed: dict[int, tuple[RawLeaf, dict[str, Any]]] = {}
    for relative, leaf in snapshot.leaves.items():
        if _SESSION.fullmatch(relative) is None:
            continue
        index, kind = _session_index(relative)
        if kind == "open":
            item = _exact_keys(_json(leaf.raw, relative), SESSION_OPEN_KEYS, relative)
            _schema_header(item, schema=SESSION_OPEN_SCHEMA, label=relative)
            if (
                item["campaign_id"] != root["campaign_id"]
                or item["session_index"] != index
                or item["namespace"] != PRODUCTION_NAMESPACE
                or item["authority_receipt_sha256"] != root["authority_receipt_sha256"]
                or item["exclusion_lock_sha256"] != _digest(exclusion_raw)
                or type(item["owner_token"]) is not str
                or not item["owner_token"]
                or type(item["host"]) is not str
                or not item["host"]
                or type(item["pid"]) is not int
                or item["pid"] < 0
                or type(item["thread_id"]) is not int
                or item["thread_id"] < 0
            ):
                _stop("PREF1_SESSION", "session open differs")
            _false(item["production_write_authorized"], "session.production_write_authorized")
            _false(item["campaign_execution_authorized"], "session.campaign_execution_authorized")
            opened[index] = (leaf, item)
        else:
            item = _exact_keys(_json(leaf.raw, relative), SESSION_CLOSE_KEYS, relative)
            _schema_header(item, schema=SESSION_CLOSE_SCHEMA, label=relative)
            closed[index] = (leaf, item)
    if not opened or set(opened) != set(range(max(opened) + 1)):
        _stop("PREF1_SESSION", "writer-session sequence is absent or has a hole")
    if set(closed) != set(opened):
        _stop("PREF1_SESSION", "completed terminal must have no unclosed writer")
    for index in sorted(opened):
        open_leaf, opening = opened[index]
        _close_leaf, closing = closed[index]
        if (
            closing["campaign_id"] != root["campaign_id"]
            or closing["session_index"] != index
            or closing["session_open_sha256"] != open_leaf.sha256
            or closing["owner_token"] != opening["owner_token"]
            or closing["authority_receipt_sha256"] != root["authority_receipt_sha256"]
        ):
            _stop("PREF1_SESSION", "session close does not bind its open")
    return len(opened)


def _bundle_state(bundle: HLT17GenerationBundle) -> dict[str, object]:
    try:
        descriptor = _json(bundle.descriptor, "accepted descriptor")
        cursor_mapping = _json(bundle.cursor_bytes, "full cursor")
        cursor = restore_hlt17_cursor(cursor_mapping)
    except (TypeError, ValueError) as error:
        raise PRO20EV1PREF1Error("PREF1_BUNDLE", "bundle reauthentication failed") from error
    identity = descriptor.get("runtime_identity")
    if type(identity) is not dict:
        _stop("PREF1_BUNDLE", "descriptor runtime identity differs")
    accepted_generation = _nonnegative_int(
        identity.get("accepted_generation"), "accepted_generation"
    )
    physical = _sha(descriptor.get("physical_state_sha256"), "physical_state_sha256")
    if (
        accepted_generation != bundle.accepted_generation
        or cursor.checkpoint_generation != accepted_generation
        or cursor.ledger.accepted_macro_step_count != accepted_generation
        or physical != bundle.physical_state_sha256
        or physical != cursor.physical_state_sha256
        or cursor_sha256(cursor) != _digest(bundle.cursor_bytes)
    ):
        _stop("PREF1_BUNDLE", "accepted bundle identities differ")
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
        "tdg6_snapshot_sha256": _sha(
            cursor_mapping.get("inherited_snapshot_sha256"),
            "inherited_snapshot_sha256",
        ),
        "parked": cursor.accepted_time.hex() == cursor.event_target.hex(),
        "public_fresh": cursor.public_fresh_permitted(),
        "step_index": cursor.accepted_step_index,
        "transaction_serial": cursor.accepted_transaction_serial,
        **hlt17_implementation_identity(),
    }
    _exact_keys(mapping, MEMBER_STATE_KEYS, "derived member state")
    if mapping["parked"] is True and (
        mapping["source_current"] != 0 or mapping["cfl_current"] != 0
    ):
        _stop("PREF1_STATE", "parked member retained a nonzero retry current")
    if mapping["kernel_transition_digest"] == mapping["cursor_sha256"]:
        _stop("PREF1_STATE", "kernel-transition hash is not the full-cursor hash")
    if mapping["cursor_sha256"] == mapping["descriptor_sha256"]:
        _stop("PREF1_STATE", "full-cursor hash collided with the descriptor hash")
    if mapping["kernel_transition_digest"] == mapping["physical_state_sha256"]:
        _stop("PREF1_STATE", "kernel-transition hash collided with the physical hash")
    return mapping


def _generation_blobs(files: Mapping[str, RawLeaf]) -> dict[str, bytes]:
    blobs: dict[str, bytes] = {}
    for relative, leaf in files.items():
        if not relative.startswith("blobs/"):
            continue
        digest = PurePosixPath(relative).name
        if leaf.sha256 != digest:
            _stop("PREF1_BLOB", "content-addressed blob name differs")
        blobs[digest] = leaf.raw
    return blobs


def _resolve_bundle_part(
    digest: str,
    *,
    new_blobs: Mapping[str, bytes],
    predecessor: bytes | None,
    seed: bool,
) -> tuple[bytes, bool]:
    if digest in new_blobs:
        return new_blobs[digest], True
    if predecessor is not None and _digest(predecessor) == digest:
        return predecessor, False
    if seed:
        _stop("PREF1_BLOB", "seed must resolve state bytes from current generation blobs")
    _stop("PREF1_BLOB", "generation references a future or rollback blob")
    raise AssertionError("unreachable")


def _bundle_from_state(
    state: Mapping[str, object],
    new_blobs: Mapping[str, bytes],
    *,
    predecessor: HLT17GenerationBundle | None,
) -> tuple[HLT17GenerationBundle, frozenset[str]]:
    used: set[str] = set()
    descriptor_digest = _sha(state.get("descriptor_sha256"), "descriptor_sha256")
    payload_digest = _sha(state.get("payload_sha256"), "payload_sha256")
    cursor_digest = _sha(state.get("cursor_sha256"), "cursor_sha256")
    if predecessor is not None:
        for digest, raw in (
            (descriptor_digest, predecessor.descriptor),
            (payload_digest, predecessor.payload),
            (cursor_digest, predecessor.cursor_bytes),
        ):
            if digest in new_blobs and _digest(raw) == digest:
                _stop("PREF1_BLOB", "content-addressed blob was reintroduced")
    descriptor, used_descriptor = _resolve_bundle_part(
        descriptor_digest,
        new_blobs=new_blobs,
        predecessor=None if predecessor is None else predecessor.descriptor,
        seed=predecessor is None,
    )
    payload, used_payload = _resolve_bundle_part(
        payload_digest,
        new_blobs=new_blobs,
        predecessor=None if predecessor is None else predecessor.payload,
        seed=predecessor is None,
    )
    cursor, used_cursor = _resolve_bundle_part(
        cursor_digest,
        new_blobs=new_blobs,
        predecessor=None if predecessor is None else predecessor.cursor_bytes,
        seed=predecessor is None,
    )
    if used_descriptor:
        used.add(descriptor_digest)
    if used_payload:
        used.add(payload_digest)
    if used_cursor:
        used.add(cursor_digest)
    try:
        bundle = HLT17GenerationBundle(descriptor, payload, cursor)
    except (TypeError, ValueError) as error:
        raise PRO20EV1PREF1Error("PREF1_BUNDLE", "bundle bytes are invalid") from error
    derived = _bundle_state(bundle)
    if dict(state) != derived:
        _stop("PREF1_STATE", "checkpoint member state differs from bundle bytes")
    return bundle, frozenset(used)


def _required_member(states: Mapping[str, Mapping[str, object]]) -> str:
    pending = [
        key
        for key in HLT17_MEMBER_KEYS
        if states[key]["mode"] != FRESH_READY or states[key]["latest_owner"] is not None
    ]
    if pending:
        if len(pending) != 1:
            _stop("PREF1_SCHEDULE", "multiple members hold retry ownership")
        return pending[0]
    for key in HLT17_MEMBER_KEYS:
        if states[key]["parked"] is not True:
            return key
    _stop("PREF1_SCHEDULE", "nonterminal cohort has no required member")
    raise AssertionError("unreachable")


def _cohort_complete(states: Mapping[str, Mapping[str, object]]) -> bool:
    return all(
        state["mode"] == FRESH_READY
        and state["latest_owner"] is None
        and state["parked"] is True
        and state["public_fresh"] is True
        and state["accepted_time_hex"] == EVENT_TARGET_HEX
        and state["source_current"] == 0
        and state["cfl_current"] == 0
        for state in states.values()
    )


def _derive_changed_kind(
    previous_bundle: HLT17GenerationBundle,
    current_bundle: HLT17GenerationBundle,
) -> str:
    previous = _bundle_state(previous_bundle)
    current = _bundle_state(current_bundle)
    cursor = restore_hlt17_cursor(_json(current_bundle.cursor_bytes, "successor cursor"))
    encoding_changed = (
        previous_bundle.descriptor != current_bundle.descriptor
        or previous_bundle.payload != current_bundle.payload
    )
    physical_preserved = (
        previous["physical_state_sha256"] == current["physical_state_sha256"]
        and previous["accepted_time_hex"] == current["accepted_time_hex"]
    )
    if encoding_changed:
        if physical_preserved or cursor.mode != FRESH_READY or not cursor.public_fresh_permitted():
            _stop("PREF1_CLASS", "physical advance has an invalid fresh successor")
        return "accepted_fine"
    if not physical_preserved:
        _stop("PREF1_CLASS", "cursor-only transition changed accepted physical state")
    if cursor.overlay is not None:
        exhausted = _overlay_exhausted(cursor)
        if cursor.overlay.owner == SOURCE_OWNER:
            return "source_retry_exhausted" if exhausted else "source_retry"
        if cursor.overlay.owner == CFL_OWNER:
            return "cfl_retry_exhausted" if exhausted else "cfl_retry"
        _stop("PREF1_CLASS", "retry overlay has a foreign owner")
    if (
        cursor.mode == RETRY_PENDING
        and cursor.latest_owner == TEMPORAL_OWNER
        and cursor.temporal is not None
        and cursor.next_attempt_plan() is not None
    ):
        return "temporal_retry"
    _stop("PREF1_CLASS", "changed cursor has no independently executable owner")
    raise AssertionError("unreachable")


def _require_fine_transition(
    previous: Mapping[str, object],
    current: Mapping[str, object],
    *,
    previous_bundle: HLT17GenerationBundle,
    current_bundle: HLT17GenerationBundle,
    member_key: str,
) -> None:
    if current["accepted_generation"] != previous["accepted_generation"] + 1:
        _stop("PREF1_CLASS", "accepted fine must increment accepted generation")
    if current["kernel_checkpoint_generation"] != previous["kernel_checkpoint_generation"] + 1:
        _stop("PREF1_CLASS", "accepted fine must increment the kernel checkpoint generation")
    if current["cursor_generation"] != previous["cursor_generation"] + 1:
        _stop("PREF1_CLASS", "accepted fine must increment cursor generation")
    if current["attempt_serial"] != previous["attempt_serial"] + 1:
        _stop("PREF1_CLASS", "accepted fine must increment attempt serial")
    if current["step_index"] != previous["step_index"] + ACCEPTED_SUBSTEPS:
        _stop("PREF1_CLASS", "accepted fine must advance exactly four substeps")
    expected_serial = _stage_records(member_key)
    if current["transaction_serial"] != previous["transaction_serial"] + expected_serial:
        _stop("PREF1_CLASS", "accepted fine must advance 20 RK4 or 16 SSPRK3 records")
    if current["descriptor_sha256"] == previous["descriptor_sha256"]:
        _stop("PREF1_CLASS", "accepted fine must publish a new descriptor")
    if current["payload_sha256"] == previous["payload_sha256"]:
        _stop("PREF1_CLASS", "accepted fine must publish a new payload")
    if current["physical_state_sha256"] == previous["physical_state_sha256"]:
        _stop("PREF1_CLASS", "accepted fine must advance accepted physical state")
    if current["cursor_sha256"] == previous["cursor_sha256"]:
        _stop("PREF1_CLASS", "accepted fine must publish a new full cursor")
    if current["mode"] != FRESH_READY or current["latest_owner"] is not None:
        _stop("PREF1_CLASS", "accepted fine must return a fresh cursor")
    if current["source_current"] != 0 or current["cfl_current"] != 0:
        _stop("PREF1_CLASS", "accepted fine must reset per-macro source/CFL counters")
    if current["source_total"] != previous["source_total"] or current["cfl_total"] != previous["cfl_total"]:
        _stop("PREF1_CLASS", "accepted fine must retain inherited source/CFL totals")
    if current["tdg6_snapshot_sha256"] != previous["tdg6_snapshot_sha256"]:
        _stop("PREF1_CLASS", "accepted fine must retain the immutable TDG6 sibling")
    descriptor = _json(current_bundle.descriptor, "fine descriptor")
    if descriptor.get("protocol_artifact_id") != RUNTIME_PROTOCOL:
        _stop("PREF1_CLASS", "fine descriptor protocol differs")
    if descriptor.get("predecessor_cursor_sha256") != previous["cursor_sha256"]:
        _stop("PREF1_CLASS", "new descriptor must name the old full-cursor hash")
    if descriptor.get("predecessor_descriptor_sha256") != previous["descriptor_sha256"]:
        _stop("PREF1_CLASS", "new descriptor must name the old accepted descriptor")
    new_cursor = restore_hlt17_cursor(_json(current_bundle.cursor_bytes, "fine cursor"))
    old_cursor = restore_hlt17_cursor(_json(previous_bundle.cursor_bytes, "predecessor cursor"))
    if new_cursor.kernel_transition_parent_sha256 != old_cursor.kernel_transition_digest:
        _stop("PREF1_CLASS", "fine kernel parent is not the predecessor kernel digest")
    if new_cursor.descriptor_sha256 == old_cursor.descriptor_sha256:
        _stop("PREF1_CLASS", "fine cursor retained the old descriptor")


def _require_retry_transition(
    previous: Mapping[str, object],
    current: Mapping[str, object],
    *,
    previous_bundle: HLT17GenerationBundle,
    current_bundle: HLT17GenerationBundle,
    kind: str,
) -> None:
    if current["accepted_generation"] != previous["accepted_generation"]:
        _stop("PREF1_CLASS", "retry cannot increment accepted generation")
    if current["kernel_checkpoint_generation"] != previous["kernel_checkpoint_generation"]:
        _stop("PREF1_CLASS", "retry cannot increment the kernel checkpoint generation")
    if current["descriptor_sha256"] != previous["descriptor_sha256"]:
        _stop("PREF1_CLASS", "retry must preserve the accepted descriptor blob")
    if current["payload_sha256"] != previous["payload_sha256"]:
        _stop("PREF1_CLASS", "retry must preserve the accepted payload blob")
    if current["physical_state_sha256"] != previous["physical_state_sha256"]:
        _stop("PREF1_CLASS", "retry must preserve accepted physical state")
    if current["accepted_time_hex"] != previous["accepted_time_hex"]:
        _stop("PREF1_CLASS", "retry must preserve accepted time")
    if current["step_index"] != previous["step_index"]:
        _stop("PREF1_CLASS", "retry must preserve the accepted step index")
    if current["transaction_serial"] != previous["transaction_serial"]:
        _stop("PREF1_CLASS", "retry must preserve the accepted transaction serial")
    if current["tdg6_snapshot_sha256"] != previous["tdg6_snapshot_sha256"]:
        _stop("PREF1_CLASS", "retry must preserve the immutable TDG6 sibling")
    if previous_bundle.descriptor != current_bundle.descriptor:
        _stop("PREF1_CLASS", "retry descriptor bytes differ")
    if previous_bundle.payload != current_bundle.payload:
        _stop("PREF1_CLASS", "retry payload bytes differ")
    cursor = restore_hlt17_cursor(_json(current_bundle.cursor_bytes, "retry cursor"))
    old = restore_hlt17_cursor(_json(previous_bundle.cursor_bytes, "predecessor cursor"))
    if kind in CURSOR_RETRY_KINDS:
        if current["cursor_sha256"] == previous["cursor_sha256"]:
            _stop("PREF1_CLASS", "retry cursor bytes must change")
        if current["cursor_generation"] != previous["cursor_generation"] + 1:
            _stop("PREF1_CLASS", "retry must increment cursor generation")
        if current["attempt_serial"] != previous["attempt_serial"] + 1:
            _stop("PREF1_CLASS", "retry must increment attempt serial")
        if cursor.kernel_transition_parent_sha256 != old.kernel_transition_digest:
            _stop("PREF1_CLASS", "retry kernel parent is not the predecessor kernel digest")
    if kind == "source_retry":
        if (
            cursor.overlay is None
            or cursor.overlay.owner != SOURCE_OWNER
            or cursor.overlay.exhausted
            or _overlay_exhausted(cursor)
        ):
            _stop("PREF1_CLASS", "source retry cursor overlay differs")
        if current["source_current"] != previous["source_current"] + 1:
            _stop("PREF1_CLASS", "source retry must increment the source current counter")
        if current["source_total"] != previous["source_total"] + 1:
            _stop("PREF1_CLASS", "source retry must increment the inherited source total")
        if current["cfl_current"] != previous["cfl_current"] or current["cfl_total"] != previous["cfl_total"]:
            _stop("PREF1_CLASS", "source retry cannot rewrite CFL counters")
        if cursor.next_attempt_plan() is None:
            _stop("PREF1_CLASS", "source retry must remain executable")
    elif kind == "cfl_retry":
        if (
            cursor.overlay is None
            or cursor.overlay.owner != CFL_OWNER
            or cursor.overlay.exhausted
            or _overlay_exhausted(cursor)
        ):
            _stop("PREF1_CLASS", "CFL retry cursor overlay differs")
        if current["cfl_current"] != previous["cfl_current"] + 1:
            _stop("PREF1_CLASS", "CFL retry must increment the CFL current counter")
        if current["cfl_total"] != previous["cfl_total"] + 1:
            _stop("PREF1_CLASS", "CFL retry must increment the inherited CFL total")
        if (
            current["source_current"] != previous["source_current"]
            or current["source_total"] != previous["source_total"]
        ):
            _stop("PREF1_CLASS", "CFL retry cannot rewrite source counters")
        if cursor.next_attempt_plan() is None:
            _stop("PREF1_CLASS", "CFL retry must remain executable")
    elif kind == "temporal_retry":
        if cursor.mode != RETRY_PENDING or cursor.temporal is None:
            _stop("PREF1_CLASS", "temporal retry must retain IMP1 temporal identity")
        if current["source_total"] != previous["source_total"] or current["cfl_total"] != previous["cfl_total"]:
            _stop("PREF1_CLASS", "temporal retry must keep inherited source/CFL totals")
        if cursor.next_attempt_plan() is None:
            _stop("PREF1_CLASS", "temporal retry omitted an executable successor plan")
    elif kind == "source_retry_exhausted":
        if cursor.overlay is None or cursor.overlay.owner != SOURCE_OWNER or not cursor.overlay.exhausted:
            _stop("PREF1_CLASS", "source exhaustion overlay differs")
        if cursor.next_attempt_plan() is not None:
            _stop("PREF1_CLASS", "source exhaustion must not keep an executable retry cursor")
    elif kind == "cfl_retry_exhausted":
        if cursor.overlay is None or cursor.overlay.owner != CFL_OWNER or not cursor.overlay.exhausted:
            _stop("PREF1_CLASS", "CFL exhaustion overlay differs")
        if cursor.next_attempt_plan() is not None:
            _stop("PREF1_CLASS", "CFL exhaustion must not keep an executable retry cursor")
    else:
        _stop("PREF1_CLASS", "retry kind differs")


def _require_temporal_exhaustion_successor(
    *,
    previous_bundle: HLT17GenerationBundle,
    terminal_evidence: Mapping[str, object],
) -> None:
    reason = _text(terminal_evidence.get("reason"), "temporal exhaustion reason")
    if reason not in TEMPORAL_EXHAUSTION_REASONS:
        _stop("PREF1_TERMINAL", "temporal exhaustion reason differs")
    ledger_mapping = terminal_evidence.get("updated_imp1_ledger")
    digest = terminal_evidence.get("updated_ledger_sha256")
    if type(ledger_mapping) is not dict:
        _stop("PREF1_TERMINAL", "temporal exhaustion must retain the updated IMP1 ledger as evidence")
    try:
        updated = restore_imp1_checkpoint_extension(ledger_mapping)
    except (TypeError, ValueError) as error:
        raise PRO20EV1PREF1Error(
            "PREF1_TERMINAL", "temporal exhaustion ledger is not a closed IMP1 ledger"
        ) from error
    encoded = canonical_json_bytes(ledger_mapping)
    if _digest(encoded) != _sha(digest, "updated_ledger_sha256"):
        _stop("PREF1_TERMINAL", "temporal exhaustion ledger hash differs")
    cursor = restore_hlt17_cursor(_json(previous_bundle.cursor_bytes, "exhaustion predecessor"))
    current = cursor.ledger
    if updated.method != current.method:
        _stop("PREF1_TERMINAL", "temporal exhaustion ledger method differs")
    if updated.origin_receipt_sha256 != current.origin_receipt_sha256:
        _stop("PREF1_TERMINAL", "temporal exhaustion origin receipt differs")
    if updated.origin_state_sha256 != current.origin_state_sha256:
        _stop("PREF1_TERMINAL", "temporal exhaustion origin state differs")
    if updated.origin_step_index != current.origin_step_index:
        _stop("PREF1_TERMINAL", "temporal exhaustion origin step differs")
    if updated.origin_transaction_serial != current.origin_transaction_serial:
        _stop("PREF1_TERMINAL", "temporal exhaustion origin serial differs")
    if updated.inherited_snapshot != current.inherited_snapshot:
        _stop("PREF1_TERMINAL", "temporal exhaustion cannot rewrite the inherited TDG6 sibling")
    if updated.last_accepted_time.hex() != current.last_accepted_time.hex():
        _stop("PREF1_TERMINAL", "temporal exhaustion cannot advance accepted time")
    if updated.current_state_sha256 != current.current_state_sha256:
        _stop("PREF1_TERMINAL", "temporal exhaustion cannot replace the accepted physical state")
    if updated.current_step_index != current.current_step_index:
        _stop("PREF1_TERMINAL", "temporal exhaustion cannot advance the step index")
    if updated.current_transaction_serial != current.current_transaction_serial:
        _stop("PREF1_TERMINAL", "temporal exhaustion cannot advance the transaction serial")
    if updated.accepted_macro_step_count != current.accepted_macro_step_count:
        _stop("PREF1_TERMINAL", "temporal exhaustion cannot change accepted macro-step count")
    if updated.accumulated_debit_vector != current.accumulated_debit_vector:
        _stop("PREF1_TERMINAL", "temporal exhaustion cannot rewrite accumulated debit")
    if (
        updated.last_accepted_macro_step_temporal_retry_count
        != current.last_accepted_macro_step_temporal_retry_count
    ):
        _stop("PREF1_TERMINAL", "temporal exhaustion cannot rewrite last-accepted temporal retry history")
    current_rejections = current.serialized_temporal_rejections
    updated_rejections = updated.serialized_temporal_rejections
    if updated_rejections[: len(current_rejections)] != current_rejections:
        _stop("PREF1_TERMINAL", "temporal exhaustion ledger is not a prefix-preserving successor")
    appended = updated_rejections[len(current_rejections) :]
    if len(appended) != 1:
        _stop("PREF1_TERMINAL", "temporal exhaustion must be the exact one-rejection successor")
    if updated.current_macro_step_temporal_retry_count != (
        current.current_macro_step_temporal_retry_count + 1
    ):
        _stop("PREF1_TERMINAL", "temporal exhaustion current retry count is not the successor")
    if updated.cumulative_temporal_retry_count != (
        current.cumulative_temporal_retry_count + 1
    ):
        _stop("PREF1_TERMINAL", "temporal exhaustion cumulative retry count is not the successor")
    try:
        expected = append_imp1_rejection(
            current, load_canonical_json(appended[0].encode("ascii"))
        )
    except (CanonicalJSONError, TypeError, ValueError) as error:
        raise PRO20EV1PREF1Error(
            "PREF1_TERMINAL",
            "temporal exhaustion rejection is not a closed successor of this ledger",
        ) from error
    if ledger_sha256(expected) != ledger_sha256(updated):
        _stop(
            "PREF1_TERMINAL",
            "updated IMP1 ledger is not the exact successor of this member cursor ledger",
        )
    last = load_canonical_json(appended[0].encode("ascii"))
    if type(last) is not dict:
        _stop("PREF1_TERMINAL", "temporal exhaustion rejection is not a JSON object")
    if reason == "maximum_temporal_retries":
        if current.current_macro_step_temporal_retry_count != TDG6_MAXIMUM_TEMPORAL_RETRIES_PER_STEP:
            _stop("PREF1_TERMINAL", "premature temporal exhaustion")
        if updated.current_macro_step_temporal_retry_count <= TDG6_MAXIMUM_TEMPORAL_RETRIES_PER_STEP:
            _stop("PREF1_TERMINAL", "premature temporal exhaustion")
    elif reason == "minimum_macro_step":
        attempted = float.fromhex(
            _text(last.get("attempted_width_hex"), "attempted_width_hex")
        )
        if attempted / 2.0 >= float(TDG6_MINIMUM_MACRO_STEP):
            _stop("PREF1_TERMINAL", "premature temporal exhaustion")
    elif reason != "coordinate_lattice":
        _stop("PREF1_TERMINAL", "temporal exhaustion reason differs")


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
        _stop("PREF1_TERMINAL", f"{kind} must not change the accepted member checkpoint state")
    if (
        previous_bundle.descriptor != current_bundle.descriptor
        or previous_bundle.payload != current_bundle.payload
        or previous_bundle.cursor_bytes != current_bundle.cursor_bytes
    ):
        _stop("PREF1_TERMINAL", f"{kind} must not publish a fabricated successor cursor")
    if type(terminal_evidence) is not dict:
        _stop("PREF1_TERMINAL", f"{kind} must carry terminal evidence")
    if terminal_evidence.get("executable_retry_cursor") is not False:
        _stop("PREF1_TERMINAL", f"{kind} evidence cannot advertise an executable retry cursor")
    if terminal_evidence.get("kind") != kind:
        _stop("PREF1_TERMINAL", f"{kind} evidence kind differs")
    if kind == "temporal_retry_exhausted":
        _require_temporal_exhaustion_successor(
            previous_bundle=previous_bundle,
            terminal_evidence=terminal_evidence,
        )
    elif kind == "invalid_premise":
        if terminal_evidence.get("spends_retry") is not False:
            _stop("PREF1_TERMINAL", "invalid premise must not spend a retry")
        _text(terminal_evidence.get("reason"), "invalid premise reason")
    elif kind == "resource_exhausted":
        _text(terminal_evidence.get("reason"), "resource exhaustion reason")
    else:
        _stop("PREF1_TERMINAL", "unchanged terminal kind differs")


def _typed_unchanged_kind(payload: Mapping[str, object]) -> str:
    evidence = payload.get("terminal_evidence")
    if type(evidence) is not dict:
        _stop("PREF1_TERMINAL", "unchanged terminal lacks typed evidence")
    kind = evidence.get("kind")
    if kind not in BLOB_UNCHANGED_KINDS:
        _stop("PREF1_TERMINAL", "unchanged terminal kind differs")
    return str(kind)


def _exhausted_retry_evidence(kind: str) -> dict[str, object]:
    return {
        "kind": kind,
        "exhausted": True,
        "executable_retry_cursor": False,
        "next_plan": None,
    }


def _require_live_seed_bundle(bundle: HLT17GenerationBundle, member_key: str) -> None:
    cursor = restore_hlt17_cursor(_json(bundle.cursor_bytes, "seed cursor"))
    descriptor = _json(bundle.descriptor, "seed descriptor")
    identity = descriptor.get("runtime_identity")
    if type(identity) is not dict:
        _stop("PREF1_SEED", "runtime identity differs")
    if identity.get("member_key") != member_key or cursor.member_key != member_key:
        _stop("PREF1_SEED", "seed bundle member_key differs")
    spec = hlt17_member_spec(member_key)
    if spec.identity_scope != IDENTITY_SCOPE_LIVE or member_key not in HLT17_MEMBER_KEYS:
        _stop("PREF1_SEED", "seed refuses synthetic or diagnostic members")
    if identity.get("scope") != IDENTITY_SCOPE_LIVE:
        _stop("PREF1_SEED", "seed bundle is not a live C1R1 member")
    if identity.get("amplitude") != LIVE_AMPLITUDE:
        _stop("PREF1_SEED", "live seed amplitude differs")
    if identity.get("method") != spec.method_label:
        _stop("PREF1_SEED", "live method identity differs")
    if identity.get("point_count") != spec.point_count:
        _stop("PREF1_SEED", "live point count differs")
    if identity.get("spatial_order") != _spatial_order(spec.method_label):
        _stop("PREF1_SEED", "live spatial operator identity differs")
    if cursor.member.identity_scope != IDENTITY_SCOPE_LIVE:
        _stop("PREF1_SEED", "seed cursor is not a live C1R1 member")
    if cursor.event_target.hex() != EVENT_TARGET_HEX:
        _stop("PREF1_SEED", "seed event target must be t=3/2")
    if cursor.accepted_time.hex() != EVENT_ORIGIN_HEX:
        _stop("PREF1_SEED", "seed accepted time must be t=23/16")
    if bundle.accepted_generation != 0 or cursor.checkpoint_generation != 0:
        _stop("PREF1_SEED", "seed accepted generation must be zero")
    if cursor.cursor_generation != 0 or cursor.attempt_serial != 0:
        _stop("PREF1_SEED", "seed cursor/attempt revisions must be zero")
    if cursor.kernel_transition_parent_sha256 != GENESIS_DIGEST:
        _stop("PREF1_SEED", "seed kernel parent must be the genesis digest")
    predecessor = descriptor.get("predecessor_cursor_identity")
    if type(predecessor) is not dict:
        _stop("PREF1_SEED", "seed predecessor identity differs")
    if predecessor.get("protocol_artifact_id") != ORIGIN_PROTOCOL:
        _stop("PREF1_SEED", "seed descriptor must name the sealed PROTO17 origin identity")
    if cursor.mode != FRESH_READY or cursor.latest_owner is not None:
        _stop("PREF1_SEED", "seed cursor must be fresh")
    if cursor.accepted_time.hex() == cursor.event_target.hex():
        _stop("PREF1_SEED", "seed cannot already be parked at the event target")
    if cursor.source_current != 0 or cursor.cfl_current != 0:
        _stop("PREF1_SEED", "seed per-macro retry counters must start at zero")


def _require_distinct_identities(
    *,
    checkpoint_sha256: str,
    journal_sha256: str,
    states: Mapping[str, Mapping[str, object]],
) -> None:
    if checkpoint_sha256 == journal_sha256:
        _stop("PREF1_LINEAGE", "checkpoint hash collided with the journal hash")
    for key, state in states.items():
        cursor_hash = _sha(state["cursor_sha256"], "cursor_sha256")
        kernel = _sha(state["kernel_transition_digest"], "kernel_transition_digest")
        identities = (checkpoint_sha256, journal_sha256, cursor_hash, kernel)
        if len(set(identities)) != 4:
            _stop(
                "PREF1_LINEAGE",
                f"{key} checkpoint/journal/cursor/kernel identities are not distinct",
            )


def _attempt_payload(
    member_key: str,
    previous: Mapping[str, object],
    current: Mapping[str, object],
    *,
    kind: str,
    terminal_evidence: object,
) -> dict[str, object]:
    result = {
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
        "accepted_state_preserved": kind != "accepted_fine",
        "tdg6_snapshot_sha256": current["tdg6_snapshot_sha256"],
        "parked": current["parked"],
        "public_fresh": current["public_fresh"],
        "executable_retry_cursor": kind in {"source_retry", "cfl_retry", "temporal_retry"},
        "step_index": current["step_index"],
        "transaction_serial": current["transaction_serial"],
        "terminal_evidence": terminal_evidence,
    }
    _exact_keys(result, ATTEMPT_PAYLOAD_KEYS, "derived attempt payload")
    return result


def _parse_generation(
    *,
    name: str,
    files: Mapping[str, RawLeaf],
    root: Mapping[str, object],
    predecessor: IndependentGeneration | None,
    seen_blobs: set[str],
) -> IndependentGeneration:
    index = _generation_index(name)
    manifest = _exact_keys(
        _json(files["generation_manifest.json"].raw, f"{name} manifest"),
        MANIFEST_KEYS,
        f"{name} manifest",
    )
    journal = _exact_keys(
        _json(files["journal.json"].raw, f"{name} journal"),
        JOURNAL_KEYS,
        f"{name} journal",
    )
    checkpoint = _exact_keys(
        _json(files["checkpoint.json"].raw, f"{name} checkpoint"),
        CHECKPOINT_KEYS,
        f"{name} checkpoint",
    )
    _schema_header(manifest, schema=MANIFEST_SCHEMA, label=f"{name} manifest")
    _schema_header(journal, schema=JOURNAL_SCHEMA, label=f"{name} journal")
    _schema_header(checkpoint, schema=CHECKPOINT_SCHEMA, label=f"{name} checkpoint")
    if any(
        _campaign_id(item["campaign_id"]) != root["campaign_id"]
        for item in (manifest, journal, checkpoint)
    ):
        _stop("PREF1_LINEAGE", "generation campaign identity differs")
    journal_sha = files["journal.json"].sha256
    checkpoint_sha = files["checkpoint.json"].sha256
    actual_leaves = {
        relative: leaf.sha256
        for relative, leaf in files.items()
        if relative != "generation_manifest.json"
    }
    new_blobs = _generation_blobs(files)
    new_blob_names = sorted(new_blobs)
    if any(digest in seen_blobs for digest in new_blobs):
        _stop("PREF1_BLOB", "content-addressed blob was reintroduced")
    if (
        manifest["store_generation"] != index
        or manifest["journal_sequence"] != index
        or manifest["journal_record_sha256"] != journal_sha
        or manifest["checkpoint_sha256"] != checkpoint_sha
        or manifest["leaves"] != actual_leaves
        or manifest["new_blob_sha256s"] != new_blob_names
    ):
        _stop("PREF1_MANIFEST", "generation manifest does not address its bytes")
    for name_false in (
        "production_write_authorized",
        "campaign_execution_authorized",
        "calibration_claimed",
        "eligibility_claimed",
        "physical_classification",
    ):
        _false(manifest[name_false], f"manifest.{name_false}")
    _nonclaims(journal)
    _nonclaims(checkpoint)
    if (
        journal["store_generation"] != index
        or journal["journal_sequence"] != index
        or checkpoint["store_generation"] != index
        or checkpoint["journal_sequence"] != index
        or checkpoint["journal_record_sha256"] != journal_sha
        or checkpoint["journal_parent_sha256"] != journal["journal_parent_sha256"]
        or checkpoint["predecessor_checkpoint_sha256"]
        != journal["predecessor_checkpoint_sha256"]
        or tuple(checkpoint["member_keys"]) != tuple(HLT17_MEMBER_KEYS)
        or checkpoint["event_origin_time_hex"] != EVENT_ORIGIN_HEX
        or checkpoint["event_target_time_hex"] != EVENT_TARGET_HEX
    ):
        _stop("PREF1_LINEAGE", "journal/checkpoint addressing differs")
    if type(checkpoint["members"]) is not dict or set(checkpoint["members"]) != set(HLT17_MEMBER_KEYS):
        _stop("PREF1_STATE", "checkpoint member order differs")
    if type(journal["payload"]) is not dict:
        _stop("PREF1_JOURNAL", "journal payload must be an object")
    states: dict[str, Mapping[str, object]] = {}
    bundles: dict[str, HLT17GenerationBundle] = {}
    addressed: set[str] = set()
    for member_key in HLT17_MEMBER_KEYS:
        state = _exact_keys(
            checkpoint["members"][member_key],
            MEMBER_STATE_KEYS,
            f"{name}.{member_key}",
        )
        prior_bundle = None if predecessor is None else predecessor.bundles[member_key]
        bundle, used = _bundle_from_state(state, new_blobs, predecessor=prior_bundle)
        addressed.update(used)
        states[member_key] = state
        bundles[member_key] = bundle
    if addressed != set(new_blobs):
        _stop("PREF1_BLOB", "generation contains a non-addressed blob")
    complete = _cohort_complete(states)
    if predecessor is None:
        _exact_keys(
            journal["payload"],
            ("member_keys", "members"),
            "seed journal payload",
        )
        if index != 0:
            _stop("PREF1_LINEAGE", "first generation is not zero")
        if (
            journal["journal_parent_sha256"] != GENESIS_DIGEST
            or journal["predecessor_checkpoint_sha256"] != GENESIS_DIGEST
            or manifest["predecessor_checkpoint_sha256"] != GENESIS_DIGEST
            or journal["kind"] != "seed"
            or checkpoint["kind"] != "seed"
            or manifest["kind"] != "seed"
            or checkpoint["changed_member_key"] is not None
            or manifest["changed_member_key"] is not None
            or journal["payload"]
            != {"member_keys": list(HLT17_MEMBER_KEYS), "members": states}
        ):
            _stop("PREF1_LINEAGE", "seed generation differs")
        for member_key, state in states.items():
            _require_live_seed_bundle(bundles[member_key], member_key)
            if (
                state["accepted_time_hex"] != EVENT_ORIGIN_HEX
                or state["accepted_generation"] != 0
                or state["cursor_generation"] != 0
                or state["attempt_serial"] != 0
                or state["mode"] != FRESH_READY
                or state["latest_owner"] is not None
                or state["parked"] is not False
                or state["source_current"] != 0
                or state["cfl_current"] != 0
            ):
                _stop("PREF1_STATE", "seed member is not the sealed origin")
        derived_kind = "seed"
        changed_member = None
    else:
        _exact_keys(
            journal["payload"],
            ATTEMPT_PAYLOAD_KEYS,
            "attempt journal payload",
        )
        if predecessor.terminal or predecessor.first_event_complete:
            _stop("PREF1_LINEAGE", "generation appears after a terminal")
        if (
            index != predecessor.index + 1
            or journal["journal_parent_sha256"] != predecessor.journal_sha256
            or journal["predecessor_checkpoint_sha256"] != predecessor.checkpoint_sha256
            or manifest["predecessor_checkpoint_sha256"] != predecessor.checkpoint_sha256
        ):
            _stop("PREF1_LINEAGE", "generation parent relation differs")
        changed_member = checkpoint["changed_member_key"]
        if changed_member != manifest["changed_member_key"] or changed_member != journal["payload"].get("member_key"):
            _stop("PREF1_LINEAGE", "changed-member identity differs")
        required = _required_member(predecessor.member_states)
        if changed_member != required:
            _stop("PREF1_SCHEDULE", "generation advanced a noncanonical member")
        changed = [
            key
            for key in HLT17_MEMBER_KEYS
            if (
                predecessor.bundles[key].descriptor != bundles[key].descriptor
                or predecessor.bundles[key].payload != bundles[key].payload
                or predecessor.bundles[key].cursor_bytes != bundles[key].cursor_bytes
            )
        ]
        current_state = states[changed_member]
        previous_state = predecessor.member_states[changed_member]
        previous_bundle = predecessor.bundles[changed_member]
        current_bundle = bundles[changed_member]
        if not changed:
            derived_kind = _typed_unchanged_kind(journal["payload"])
            terminal_evidence = journal["payload"].get("terminal_evidence")
            _require_unchanged_attempt(
                previous_state,
                current_state,
                previous_bundle=previous_bundle,
                current_bundle=current_bundle,
                kind=derived_kind,
                terminal_evidence=terminal_evidence,
            )
        elif changed == [changed_member]:
            derived_kind = _derive_changed_kind(previous_bundle, current_bundle)
            if derived_kind == "accepted_fine":
                _require_fine_transition(
                    previous_state,
                    current_state,
                    previous_bundle=previous_bundle,
                    current_bundle=current_bundle,
                    member_key=changed_member,
                )
                terminal_evidence = None
            elif derived_kind in CURSOR_RETRY_KINDS:
                _require_retry_transition(
                    previous_state,
                    current_state,
                    previous_bundle=previous_bundle,
                    current_bundle=current_bundle,
                    kind=derived_kind,
                )
                terminal_evidence = (
                    None
                    if derived_kind not in TERMINAL_KINDS
                    else _exhausted_retry_evidence(derived_kind)
                )
            else:
                _stop("PREF1_CLASS", "changed cursor has no independently executable owner")
                raise AssertionError("unreachable")
            if journal["payload"].get("terminal_evidence") != terminal_evidence:
                _stop("PREF1_JOURNAL", "attempt terminal evidence differs")
        else:
            _stop("PREF1_STATE", "attempt changed zero/multiple foreign members")
            raise AssertionError("unreachable")
        expected_payload = _attempt_payload(
            changed_member,
            previous_state,
            current_state if derived_kind not in BLOB_UNCHANGED_KINDS else previous_state,
            kind=derived_kind,
            terminal_evidence=journal["payload"].get("terminal_evidence"),
        )
        if journal["payload"] != expected_payload:
            _stop("PREF1_JOURNAL", "attempt journal differs from independent state reduction")
    if journal["kind"] != derived_kind or checkpoint["kind"] != derived_kind or manifest["kind"] != derived_kind:
        _stop("PREF1_CLASS", "stored kind differs from independent reduction")
    if complete and derived_kind != "accepted_fine":
        _stop("PREF1_CLASS", "first-event completion must be the last accepted fine")
    expected_disposition = (
        "first_event_complete"
        if complete
        else derived_kind
        if derived_kind in TERMINAL_KINDS
        else "nonterminal"
    )
    expected_terminal = complete or derived_kind in TERMINAL_KINDS
    if (
        checkpoint["first_event_complete"] is not complete
        or checkpoint["terminal"] is not expected_terminal
        or checkpoint["disposition"] != expected_disposition
    ):
        _stop("PREF1_CLASS", "checkpoint terminal reduction differs")
    _require_distinct_identities(
        checkpoint_sha256=checkpoint_sha,
        journal_sha256=journal_sha,
        states=states,
    )
    seen_blobs.update(new_blobs)
    return IndependentGeneration(
        index=index,
        kind=str(journal["kind"]),
        derived_kind=derived_kind,
        changed_member_key=changed_member,
        checkpoint_sha256=checkpoint_sha,
        journal_sha256=journal_sha,
        disposition=expected_disposition,
        first_event_complete=complete,
        terminal=expected_terminal,
        member_states=states,
        bundles=bundles,
    )


def _terminal_lock(
    raw: bytes,
    *,
    root: Mapping[str, object],
    last: IndependentGeneration,
) -> dict[str, Any]:
    item = _exact_keys(_json(raw, "terminal lock"), TERMINAL_LOCK_KEYS, "terminal lock")
    _schema_header(item, schema=TERMINAL_LOCK_SCHEMA, label="terminal lock")
    if (
        item["campaign_id"] != root["campaign_id"]
        or item["checkpoint_sha256"] != last.checkpoint_sha256
        or item["journal_record_sha256"] != last.journal_sha256
        or item["store_generation"] != last.index
        or item["disposition"] != last.disposition
        or item["changed_member_key"] != last.changed_member_key
        or item["first_event_complete"] is not last.first_event_complete
    ):
        _stop("PREF1_TERMINAL", "terminal lock does not bind the final generation")
    for name in (
        "calibration_claimed",
        "eligibility_claimed",
        "physical_classification",
        "calibration_eligible",
        "sgbl_claimed",
        "def1_claimed",
        "candidate_claimed",
        "mechanism_claimed",
        "physics_claimed",
    ):
        _false(item[name], f"terminal_lock.{name}")
    return item


def bind_terminal_snapshot(snapshot: StoreSnapshot) -> IndependentTerminal:
    """Authenticate and independently reduce one already-stable terminal snapshot."""

    generation_names = _generation_names(snapshot)
    if not generation_names:
        _stop("PREF1_LINEAGE", "store has no generation")
    _root_leaf_grammar(snapshot, generation_names)
    root_leaf = snapshot.leaves[ROOT_MANIFEST_NAME]
    root = _root_manifest(root_leaf.raw)
    exclusion_leaf = snapshot.leaves[WRITER_EXCLUSION_NAME]
    _writer_exclusion(exclusion_leaf.raw, root)
    session_count = _closed_sessions(snapshot, root, exclusion_leaf.raw)
    generation_files: dict[str, dict[str, RawLeaf]] = {}
    for name in generation_names:
        generation_files[name] = _generation_leaf_map(snapshot, name)
    generations: list[IndependentGeneration] = []
    predecessor: IndependentGeneration | None = None
    seen_blobs: set[str] = set()
    for name in generation_names:
        generation = _parse_generation(
            name=name,
            files=generation_files[name],
            root=root,
            predecessor=predecessor,
            seen_blobs=seen_blobs,
        )
        generations.append(generation)
        predecessor = generation
    last = generations[-1]
    if not last.terminal:
        _stop("PREF1_TERMINAL", "last generation is nonterminal")
    lock_leaf = snapshot.leaves[TERMINAL_LOCK_NAME]
    _terminal_lock(lock_leaf.raw, root=root, last=last)
    terminal_kind = "first_event_complete" if last.first_event_complete else last.derived_kind
    return IndependentTerminal(
        campaign_id=str(root["campaign_id"]),
        generations=tuple(generations),
        disposition=last.disposition,
        first_event_complete=last.first_event_complete,
        terminal_kind=terminal_kind,
        checkpoint_sha256=last.checkpoint_sha256,
        journal_sha256=last.journal_sha256,
        terminal_lock_sha256=lock_leaf.sha256,
        session_count=session_count,
        store_leaf_count=len(snapshot.leaves),
        store_byte_count=snapshot.byte_count,
        store_lexical_sha256=snapshot.lexical_sha256,
    )


def bind_terminal_store(repository: Path) -> IndependentTerminal:
    """Read and bind one complete terminal store.  Never writes or replays."""

    return bind_terminal_snapshot(snapshot_store(Path(repository)))


def bind_independent_seed(
    terminal: IndependentTerminal,
    expected_bundles: Mapping[str, HLT17GenerationBundle],
) -> dict[str, object]:
    """Bind generation zero to separately reconstructed physical-origin bundles.

    The caller owns construction of the historical origin and physical source
    image.  This function accepts only the resulting immutable bundles; it
    neither imports the PRO20 factory nor treats root-manifest authentication
    labels as proof.
    """

    if not isinstance(terminal, IndependentTerminal):
        raise TypeError("terminal must be an IndependentTerminal")
    if not isinstance(expected_bundles, Mapping) or set(expected_bundles) != set(
        HLT17_MEMBER_KEYS
    ):
        _stop("PREF1_SEED", "independent seed must contain exactly six live members")
    seed = terminal.generations[0]
    if seed.derived_kind != "seed" or seed.index != 0:
        _stop("PREF1_SEED", "terminal generation zero is not a seed")
    identities: dict[str, dict[str, str]] = {}
    for key in HLT17_MEMBER_KEYS:
        supplied = expected_bundles[key]
        if not isinstance(supplied, HLT17GenerationBundle):
            raise TypeError(f"expected_bundles[{key!r}] must be HLT17GenerationBundle")
        try:
            rebuilt = HLT17GenerationBundle(
                bytes(supplied.descriptor),
                bytes(supplied.payload),
                bytes(supplied.cursor_bytes),
            )
        except (TypeError, ValueError) as error:
            raise PRO20EV1PREF1Error(
                "PREF1_SEED", f"independent seed bundle failed for {key}"
            ) from error
        raw = seed.bundles[key]
        if (
            raw.descriptor != rebuilt.descriptor
            or raw.payload != rebuilt.payload
            or raw.cursor_bytes != rebuilt.cursor_bytes
        ):
            _stop("PREF1_SEED", f"raw seed differs from independent origin for {key}")
        identities[key] = {
            "descriptor_sha256": _digest(rebuilt.descriptor),
            "payload_sha256": _digest(rebuilt.payload),
            "cursor_sha256": _digest(rebuilt.cursor_bytes),
            "physical_state_sha256": rebuilt.physical_state_sha256,
        }
    stream = canonical_json_bytes(
        {"member_keys": list(HLT17_MEMBER_KEYS), "identities": identities}
    )
    return {
        "member_keys": list(HLT17_MEMBER_KEYS),
        "identities": identities,
        "identity_stream_sha256": _digest(stream),
        "independently_reconstructed_origin_matches_raw_seed": True,
        "root_manifest_origin_label_sufficient": False,
        "state_advanced_by_seed_binding": False,
    }


def independent_outcome(terminal: IndependentTerminal) -> dict[str, object]:
    """Reduce the authenticated terminal without opening the next campaign gate."""

    if not isinstance(terminal, IndependentTerminal):
        raise TypeError("terminal must be an IndependentTerminal")
    if terminal.first_event_complete:
        if terminal.terminal_kind != "first_event_complete":
            _stop("PREF1_CLASS", "successful terminal kind differs")
        classification = "independently_bound_six_member_first_event_complete"
        common_event = True
    else:
        if terminal.terminal_kind not in TERMINAL_KINDS:
            _stop("PREF1_CLASS", "typed terminal kind differs")
        classification = (
            f"independently_bound_{terminal.terminal_kind}_no_common_event"
        )
        common_event = False
    return {
        "classification": classification,
        "terminal_kind": terminal.terminal_kind,
        "first_event_complete": common_event,
        "common_event_completed": common_event,
        "runner_label_sufficient": False,
        "calibration_eligible": False,
        "calibration_authorized": False,
        "sgbl_health_established": False,
        "def1_error_map_passed": False,
        "holdout_authorized": False,
        "candidate_execution_authorized": False,
        "mechanism_claimed": False,
        "physics_claimed": False,
    }


def independent_progression(terminal: IndependentTerminal) -> dict[str, object]:
    """Summarize authenticated progression without interpreting it as calibration."""

    if not isinstance(terminal, IndependentTerminal):
        raise TypeError("terminal must be an IndependentTerminal")
    counts = {kind: 0 for kind in RECORD_KINDS}
    for generation in terminal.generations:
        counts[generation.derived_kind] += 1
    final = terminal.generations[-1]
    parked = tuple(
        key for key in HLT17_MEMBER_KEYS if final.member_states[key]["parked"] is True
    )
    accepted_generations = {
        key: int(final.member_states[key]["accepted_generation"])
        for key in HLT17_MEMBER_KEYS
    }
    accepted_times = {
        key: str(final.member_states[key]["accepted_time_hex"])
        for key in HLT17_MEMBER_KEYS
    }
    accepted_fine_count = counts["accepted_fine"]
    return {
        "generation_count_including_seed": len(terminal.generations),
        "attempt_count": len(terminal.generations) - 1,
        "kind_counts": counts,
        "accepted_fine_count": accepted_fine_count,
        "accepted_state_advanced": accepted_fine_count > 0,
        "parked_member_keys": list(parked),
        "parked_member_count": len(parked),
        "accepted_generation_by_member": accepted_generations,
        "accepted_time_hex_by_member": accepted_times,
        "all_six_parked": len(parked) == len(HLT17_MEMBER_KEYS),
        "common_event_completed": terminal.first_event_complete,
        "calibration_eligible": False,
        "binder_mutated_store": False,
        "runner_label_sufficient": False,
    }


__all__ = [
    "ARTIFACT_ID",
    "IndependentGeneration",
    "IndependentTerminal",
    "PRO20EV1PREF1Error",
    "PRODUCTION_NAMESPACE",
    "RawLeaf",
    "StoreSnapshot",
    "bind_terminal_snapshot",
    "bind_terminal_store",
    "bind_independent_seed",
    "independent_outcome",
    "independent_progression",
    "snapshot_store",
]
