"""Read-only historical-origin capture for the prospective PRO20 first event.

This module authenticates the sealed PREF27 eight-leaf identities and the
expanded proto17/calibration HLT16 tree, then pins generation-1 bridge state
as immutable captured bytes.  It does not construct GR-0 shells, run a
source/RHS, repair a store, publish, or authorize a campaign.
"""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
import json
import os
from pathlib import Path
import stat
from types import MappingProxyType, SimpleNamespace
from typing import Any, Mapping

import numpy as np

from recursive_horizons.evidence_io import (
    CanonicalJSONError,
    UnsafePathError,
    canonical_json_bytes,
    load_canonical_json,
    read_regular_file,
)

from .hlt16_campaign_schema import CampaignCheckpoint, HLT16CampaignSchemaError
from .hlt16_campaign_store import HLT16CampaignStore, HLT16CampaignStoreError
from .hlt16_member_codec import (
    ARRAY_NAMES,
    GEN0_PROTOCOL,
    HLT16MemberSnapshot,
    validate_codec_metadata,
)
from .hlt16_state_store import (
    STATE_SCHEMA,
    HLT16StateStoreError,
    decode_payload,
    digest as hlt16_digest,
)
from .protocol_v17 import MEMBER_KEYS, STATE_ARRAY_NAMES


SELECTED_METHOD = "TDG11-IMP1"
SOURCE_TREE_RELATIVE = "runs/fgc-2-sf1/proto17"
SOURCE_STORE_RELATIVE = "runs/fgc-2-sf1/proto17/calibration"
PREF27_COMPACT_RELATIVE = "results/fgc-1-pro18-pref27.json"
CLEANUP_DOMAIN = "CLEANUP-BASELINE-TREE-v1"
SELECTED_GENERATION = 1
ORIGIN_EVENT = 23
ORIGIN_TIME_RATIONAL = "23/16"
ORIGIN_TIME_HEX = "0x1.7000000000000p+0"
ORIGIN_TIME = MappingProxyType(
    {"rational": ORIGIN_TIME_RATIONAL, "binary64_hex": ORIGIN_TIME_HEX}
)
CAPTURE_SCHEMA = "FGC-1-PRO20-historical-origin-capture-v1"
HLT16_DESCRIPTOR_SCHEMA = STATE_SCHEMA
GEN0_BRIDGES = frozenset(
    {
        "authenticated_GEN0_predecessor",
        "authenticated_GEN0_campaign_rebase",
    }
)
_HLT16_CODEC_KEYS = frozenset(
    {
        "runtime_identity",
        "template",
        "template_sha256",
        "accepted_plan",
        "counters",
        "tdg6",
        "monitor",
        "causal",
        "tracer",
        "provenance_cursor",
    }
)


class Pro20OriginError(ValueError):
    """The historical origin cannot be captured as a closed read-only slice."""


@dataclass(frozen=True, slots=True)
class _OwnerPins:
    selected_method: str
    selected_generation: int
    origin_event: int
    origin_time_rational: str
    origin_time_hex: str
    generation1_bridge_content_id: str
    expanded_leaf_count: int
    expanded_byte_count: int
    expanded_cleanup_digest: str
    pref27_compact_raw_sha256: str
    pref27_tree_sha256: str
    hlt15_checkpoint_content_id: str
    hlt15_checkpoint_raw_sha256: str
    hlt15_receipt_content_id: str
    hlt15_receipt_raw_sha256: str
    member_counters: Mapping[str, Mapping[str, int]]


_OWNER_PINS = _OwnerPins(
    selected_method=SELECTED_METHOD,
    selected_generation=SELECTED_GENERATION,
    origin_event=ORIGIN_EVENT,
    origin_time_rational=ORIGIN_TIME_RATIONAL,
    origin_time_hex=ORIGIN_TIME_HEX,
    generation1_bridge_content_id=(
        "23f4a1ac64af5e1463862bf67ca240980921611a4a77d31050b0f92c81d44895"
    ),
    expanded_leaf_count=56,
    expanded_byte_count=6054860,
    expanded_cleanup_digest=(
        "8e0b1f3ba70a4de7061d322ed041ff0328ae4e6ca7ed78da79421883f47ad777"
    ),
    pref27_compact_raw_sha256=(
        "102b385d058066d668e9e8541c974c736c583f6d5c4afaf9b99e95a16aabeb83"
    ),
    pref27_tree_sha256=(
        "95a1bb96570216b8007e32a0cc05f4c31f420052f1d90c1b7134fa27ce958585"
    ),
    hlt15_checkpoint_content_id=(
        "7c059dcc2197a9f57baf8012f4113c240917090171ee06cc639e9e8734603d2d"
    ),
    hlt15_checkpoint_raw_sha256=(
        "9bcbc68ca2f168d9f7ae87c37346d9ebf90aeaf2568f26546441634734511b3b"
    ),
    hlt15_receipt_content_id=(
        "78ae89a9f74193658003cf51f8567c4088db80fd6fb86867c2d581a0c16d5baf"
    ),
    hlt15_receipt_raw_sha256=(
        "7e99fb5fb70b6d2d590232d21222e9aa510b823dfdb3342aa86489b21d4d2f4a"
    ),
    member_counters=MappingProxyType(
        {
            "RK4-2049": MappingProxyType(
                {
                    "step_index": 342,
                    "transaction_serial": 1710,
                    "source_retry_count": 0,
                    "CFL_retry_count": 217,
                }
            ),
            "RK4-4097": MappingProxyType(
                {
                    "step_index": 673,
                    "transaction_serial": 3365,
                    "source_retry_count": 0,
                    "CFL_retry_count": 443,
                }
            ),
            "RK4-8193": MappingProxyType(
                {
                    "step_index": 976,
                    "transaction_serial": 4880,
                    "source_retry_count": 0,
                    "CFL_retry_count": 161,
                }
            ),
            "SSPRK3-4097": MappingProxyType(
                {
                    "step_index": 631,
                    "transaction_serial": 2524,
                    "source_retry_count": 0,
                    "CFL_retry_count": 358,
                }
            ),
            "SSPRK3-8193": MappingProxyType(
                {
                    "step_index": 987,
                    "transaction_serial": 3948,
                    "source_retry_count": 0,
                    "CFL_retry_count": 187,
                }
            ),
            "SSPRK3-16385": MappingProxyType(
                {
                    "step_index": 1771,
                    "transaction_serial": 7084,
                    "source_retry_count": 0,
                    "CFL_retry_count": 0,
                }
            ),
        }
    ),
)


def _sha256_hex(raw: bytes) -> str:
    if type(raw) is not bytes:
        raise Pro20OriginError("captured payload is not immutable bytes")
    return sha256(raw).hexdigest()


def _require_sha(value: object, label: str) -> str:
    if not isinstance(value, str) or len(value) != 64 or value.lower() != value:
        raise Pro20OriginError(f"{label} is not SHA-256")
    try:
        int(value, 16)
    except ValueError as error:
        raise Pro20OriginError(f"{label} is not SHA-256") from error
    return value


def _mapping(value: object, label: str) -> dict[str, Any]:
    if not isinstance(value, Mapping):
        raise Pro20OriginError(f"{label} is not a mapping")
    return dict(value)


def _load_json_object(raw: bytes, label: str, *, canonical: bool) -> dict[str, Any]:
    if type(raw) is not bytes:
        raise Pro20OriginError(f"{label} is not immutable bytes")

    def reject_duplicates(
        pairs: list[tuple[str, Any]],
    ) -> dict[str, Any]:
        answer: dict[str, Any] = {}
        for key, item in pairs:
            if key in answer:
                raise Pro20OriginError(f"{label} has duplicate keys")
            answer[key] = item
        return answer

    try:
        if canonical:
            value = load_canonical_json(raw)
        else:
            value = json.loads(
                raw.decode("utf-8"),
                object_pairs_hook=reject_duplicates,
            )
    except (CanonicalJSONError, UnicodeDecodeError, json.JSONDecodeError, ValueError) as error:
        raise Pro20OriginError(f"{label} is malformed JSON") from error
    if not isinstance(value, dict):
        raise Pro20OriginError(f"{label} is not a JSON object")
    return value


def physical_array_bytes(array: object) -> bytes:
    """Return little-endian C-order float64 bytes; identity is the bytes."""

    if type(array) is bytes:
        return array
    if not isinstance(array, np.ndarray):
        raise Pro20OriginError("physical array cannot be captured as bytes")
    normalized = np.asarray(array, dtype=np.dtype("<f8"), order="C")
    if (
        normalized.dtype != np.dtype("<f8")
        or not normalized.flags.c_contiguous
        or not np.isfinite(normalized).all()
    ):
        raise Pro20OriginError("physical array is not finite little-endian C float64")
    return bytes(normalized.tobytes(order="C"))


def physical_array_sha256(array: object) -> str:
    """Hash only little-endian C bytes; this is not a content ID."""

    return _sha256_hex(physical_array_bytes(array))


def full_cursor_bytes(cursor: Mapping[str, Any]) -> bytes:
    """Canonical bytes of the complete original cursor mapping."""

    return canonical_json_bytes(_mapping(cursor, "original cursor"))


def full_cursor_sha256(cursor: Mapping[str, Any]) -> str:
    """Hash of the full original cursor; distinct from cursor_chain_sha256."""

    return _sha256_hex(full_cursor_bytes(cursor))


def cleanup_baseline_tree_digest(
    records: Mapping[str, Mapping[str, Any]] | tuple[tuple[str, int, str], ...]
) -> str:
    """Recompute the owner CLEANUP-BASELINE-TREE-v1 digest.

    The stream is the domain line followed by lexically sorted
    ``path NUL byte-count NUL leaf-SHA-256 LF`` records.
    """

    rows: list[tuple[str, int, str]]
    if isinstance(records, Mapping):
        rows = []
        for path, item in records.items():
            leaf = _mapping(item, "inventory leaf")
            rows.append(
                (
                    str(path),
                    int(leaf["byte_count"]),
                    _require_sha(leaf["raw_sha256"], "inventory raw hash"),
                )
            )
    else:
        rows = list(records)
    stream = bytearray(f"{CLEANUP_DOMAIN}\n".encode("ascii"))
    for path, byte_count, digest in sorted(rows, key=lambda item: item[0]):
        if not isinstance(path, str) or not path or type(byte_count) is not int:
            raise Pro20OriginError("cleanup inventory record differs")
        _require_sha(digest, "cleanup leaf hash")
        stream.extend(path.encode("ascii"))
        stream.append(0)
        stream.extend(str(byte_count).encode("ascii"))
        stream.append(0)
        stream.extend(digest.encode("ascii"))
        stream.append(0x0A)
    return _sha256_hex(bytes(stream))


def _require_canonical_root(repository_root: Path | str) -> Path:
    """Refuse a relative or symlink-traversing root; do not rewrite it."""

    try:
        root = Path(repository_root)
    except TypeError as error:
        raise Pro20OriginError("repository root is absent") from error
    if not root.is_absolute() or any(part in {".", ".."} for part in root.parts):
        raise Pro20OriginError("repository root must be an absolute canonical path")
    try:
        info = root.lstat()
    except OSError as error:
        raise Pro20OriginError("repository root is absent") from error
    if stat.S_ISLNK(info.st_mode):
        raise Pro20OriginError("repository root traverses a symlink")
    if not stat.S_ISDIR(info.st_mode):
        raise Pro20OriginError("repository root is unsafe")
    current = root
    while True:
        try:
            current_info = current.lstat()
        except OSError as error:
            raise Pro20OriginError("repository root is unsafe") from error
        if stat.S_ISLNK(current_info.st_mode):
            raise Pro20OriginError("repository root traverses a symlink")
        if current == current.parent:
            break
        current = current.parent
    try:
        resolved = root.resolve(strict=True)
    except OSError as error:
        raise Pro20OriginError("repository root is absent") from error
    if resolved != root:
        raise Pro20OriginError("repository root traverses a symlink")
    return root


def _read_bytes(root: Path, relative: str, label: str) -> bytes:
    try:
        return read_regular_file(root, relative)
    except UnsafePathError as error:
        raise Pro20OriginError(f"{label} is absent or unsafe") from error


def _proxy(value: Mapping[str, Any]) -> Mapping[str, Any]:
    return MappingProxyType(dict(value))


def _is_checkpoint_path(path: str) -> bool:
    return "/checkpoints/" in path and path.endswith(".json")


def _is_receipt_path(path: str) -> bool:
    return path.endswith("/receipts/generation-zero.json")


def _is_state_path(path: str) -> bool:
    return "/states/" in path and path.endswith(".json")


def _content_id_for_leaf(path: str, payload: Mapping[str, Any], raw: bytes) -> str:
    """Own content IDs by path/schema, not the first digest-looking key."""

    if _is_receipt_path(path):
        return _require_sha(payload.get("receipt_sha256"), "receipt content ID")
    if _is_checkpoint_path(path):
        return _require_sha(payload.get("checkpoint_sha256"), "checkpoint content ID")
    if _is_state_path(path):
        schema = payload.get("schema")
        if schema == HLT16_DESCRIPTOR_SCHEMA or "descriptor_sha256" in payload:
            return _require_sha(payload.get("descriptor_sha256"), "descriptor content ID")
        return _sha256_hex(raw)
    raise Pro20OriginError(f"historical leaf schema is unknown: {path}")


@dataclass(frozen=True, slots=True)
class CapturedLeaf:
    path: str
    raw_sha256: str
    content_id: str
    raw_bytes: bytes

    def __post_init__(self) -> None:
        if type(self.raw_bytes) is not bytes:
            raise Pro20OriginError(f"{self.path} is not immutable bytes")
        if _sha256_hex(self.raw_bytes) != self.raw_sha256:
            raise Pro20OriginError(f"{self.path} raw identity drifted")
        _require_sha(self.raw_sha256, f"{self.path} raw hash")
        _require_sha(self.content_id, f"{self.path} content ID")


@dataclass(frozen=True, slots=True)
class SourceInventory:
    relative_root: str
    leaves: tuple[CapturedLeaf, ...]
    byte_count: int
    cleanup_digest: str

    def __post_init__(self) -> None:
        object.__setattr__(self, "leaves", tuple(self.leaves))
        records: list[tuple[str, int, str]] = []
        seen: set[str] = set()
        for leaf in self.leaves:
            if leaf.path in seen:
                raise Pro20OriginError("source inventory path is duplicated")
            seen.add(leaf.path)
            records.append((leaf.path, len(leaf.raw_bytes), leaf.raw_sha256))
        byte_count = sum(len(leaf.raw_bytes) for leaf in self.leaves)
        digest = cleanup_baseline_tree_digest(tuple(records))
        if byte_count != self.byte_count or digest != self.cleanup_digest:
            raise Pro20OriginError("source inventory identity drifted")

    @property
    def leaf_count(self) -> int:
        return len(self.leaves)

    def by_path(self) -> dict[str, CapturedLeaf]:
        return {leaf.path: leaf for leaf in self.leaves}


@dataclass(frozen=True, slots=True)
class CapturedMemberOrigin:
    member_key: str
    descriptor_content_id: str
    payload_raw_sha256: str
    hlt15_state_content_id: str
    descriptor_bytes: bytes
    payload_bytes: bytes
    physical_array_bytes: Mapping[str, bytes]
    physical_array_sha256: Mapping[str, str]
    physical_array_shapes: Mapping[str, tuple[int, ...]]
    original_cursor_bytes: bytes
    full_cursor_sha256: str
    cursor_chain_sha256: str
    input_hash: str
    monitor_bytes: bytes
    monitor_sha256: str
    causal_bytes: bytes
    causal_sha256: str
    step_index: int
    transaction_serial: int
    source_retry_count: int
    cfl_retry_count: int
    accepted_time: Mapping[str, str]

    def __post_init__(self) -> None:
        arrays = MappingProxyType(
            {name: bytes(raw) for name, raw in dict(self.physical_array_bytes).items()}
        )
        hashes = MappingProxyType(dict(self.physical_array_sha256))
        shapes = MappingProxyType(
            {
                name: tuple(int(dimension) for dimension in shape)
                for name, shape in dict(self.physical_array_shapes).items()
            }
        )
        object.__setattr__(self, "physical_array_bytes", arrays)
        object.__setattr__(self, "physical_array_sha256", hashes)
        object.__setattr__(self, "physical_array_shapes", shapes)
        object.__setattr__(self, "accepted_time", _proxy(self.accepted_time))
        if type(self.descriptor_bytes) is not bytes or type(self.payload_bytes) is not bytes:
            raise Pro20OriginError(f"{self.member_key} bridge bytes are not immutable")
        if tuple(arrays) != ARRAY_NAMES or tuple(shapes) != ARRAY_NAMES:
            raise Pro20OriginError(f"{self.member_key} bridge array order differs")
        for name in ARRAY_NAMES:
            expected_size = 8
            for dimension in shapes[name]:
                expected_size *= dimension
            if _sha256_hex(arrays[name]) != hashes[name] or len(arrays[name]) != expected_size:
                raise Pro20OriginError(
                    f"{self.member_key} {name} physical-array identity drifted"
                )
        if _sha256_hex(self.payload_bytes) != self.payload_raw_sha256:
            raise Pro20OriginError(f"{self.member_key} payload identity drifted")
        descriptor = _load_json_object(
            self.descriptor_bytes, f"{self.member_key} descriptor", canonical=True
        )
        descriptor_id = _content_id_for_leaf(
            f"{SOURCE_STORE_RELATIVE}/states/{self.descriptor_content_id}.json",
            descriptor,
            self.descriptor_bytes,
        )
        if descriptor_id != self.descriptor_content_id:
            raise Pro20OriginError(f"{self.member_key} descriptor content ID drifted")
        if type(self.original_cursor_bytes) is not bytes:
            raise Pro20OriginError(f"{self.member_key} cursor is not immutable bytes")
        if _sha256_hex(self.original_cursor_bytes) != self.full_cursor_sha256:
            raise Pro20OriginError(f"{self.member_key} full cursor identity drifted")
        if _sha256_hex(self.monitor_bytes) != self.monitor_sha256:
            raise Pro20OriginError(f"{self.member_key} monitor identity drifted")
        if _sha256_hex(self.causal_bytes) != self.causal_sha256:
            raise Pro20OriginError(f"{self.member_key} causal identity drifted")
        _require_sha(self.descriptor_content_id, f"{self.member_key} descriptor")
        _require_sha(self.payload_raw_sha256, f"{self.member_key} payload raw hash")
        _require_sha(self.hlt15_state_content_id, f"{self.member_key} HLT15 state")
        _require_sha(self.cursor_chain_sha256, f"{self.member_key} cursor chain")
        _require_sha(self.input_hash, f"{self.member_key} input hash")
        if self.full_cursor_sha256 == self.cursor_chain_sha256:
            raise Pro20OriginError(
                f"{self.member_key} full cursor hash collapsed into cursor chain ID"
            )
        if self.descriptor_content_id == self.payload_raw_sha256:
            raise Pro20OriginError(
                f"{self.member_key} descriptor content ID collapsed into payload raw hash"
            )


@dataclass(frozen=True, slots=True)
class Pro20OriginCapture:
    selected_method: str
    origin_event: int
    origin_time: Mapping[str, str]
    generation1_bridge_content_id: str
    generation1_checkpoint: CapturedLeaf
    pref27_compact: CapturedLeaf
    pref27_tree_sha256: str
    original_leaves: tuple[CapturedLeaf, ...]
    source_inventory: SourceInventory
    members: tuple[CapturedMemberOrigin, ...]
    captured_bytes: bytes

    def __post_init__(self) -> None:
        object.__setattr__(self, "origin_time", _proxy(self.origin_time))
        object.__setattr__(self, "original_leaves", tuple(self.original_leaves))
        object.__setattr__(self, "members", tuple(self.members))
        if type(self.captured_bytes) is not bytes:
            raise Pro20OriginError("capture record is not immutable bytes")
        if self.selected_method != SELECTED_METHOD:
            raise Pro20OriginError("selected method is not sealed TDG11-IMP1")
        if tuple(item.member_key for item in self.members) != MEMBER_KEYS:
            raise Pro20OriginError("captured members are not in canonical six-member order")
        if load_canonical_json(self.captured_bytes) != self.identity_record():
            raise Pro20OriginError("captured identity drifted behind cached metadata")
        reconstructed = _identity_from_retained_bytes(self)
        if canonical_json_bytes(reconstructed) != self.captured_bytes:
            raise Pro20OriginError("captured labels are not bound to retained origin bytes")

    def identity_record(self) -> dict[str, Any]:
        return _identity_record(
            selected_method=self.selected_method,
            origin_event=self.origin_event,
            origin_time=self.origin_time,
            generation1_bridge_content_id=self.generation1_bridge_content_id,
            generation1_checkpoint=self.generation1_checkpoint,
            pref27_compact=self.pref27_compact,
            pref27_tree_sha256=self.pref27_tree_sha256,
            original_leaves=self.original_leaves,
            source_inventory=self.source_inventory,
            members=self.members,
        )

    @property
    def capture_sha256(self) -> str:
        return _sha256_hex(self.captured_bytes)

    def fresh_physical_arrays(self, member_key: str) -> dict[str, np.ndarray]:
        """Decoded HLT15 physical arrays; callers receive copies."""

        return {
            name: array
            for name, array in self.fresh_bridge_arrays(member_key).items()
            if name in STATE_ARRAY_NAMES
        }

    def fresh_bridge_arrays(self, member_key: str) -> dict[str, np.ndarray]:
        """Decoded HLT16 snapshot arrays, including static grid and labels."""

        member = {item.member_key: item for item in self.members}.get(member_key)
        if member is None:
            raise Pro20OriginError("captured member is absent")
        return {
            name: np.frombuffer(raw, dtype="<f8")
            .reshape(member.physical_array_shapes[name])
            .copy()
            for name, raw in member.physical_array_bytes.items()
        }

    def retained_bridge_payload(self, member_key: str) -> tuple[bytes, bytes]:
        """Immutable descriptor JSON and payload NPZ for a later HLT16 restore."""

        member = {item.member_key: item for item in self.members}.get(member_key)
        if member is None:
            raise Pro20OriginError("captured member is absent")
        return bytes(member.descriptor_bytes), bytes(member.payload_bytes)

    def fresh_hlt16_snapshot(self, member_key: str) -> HLT16MemberSnapshot:
        """Reconstruct one fresh codec snapshot from retained authenticated bytes."""

        member = {item.member_key: item for item in self.members}.get(member_key)
        if member is None:
            raise Pro20OriginError("captured member is absent")
        descriptor, metadata, address = _decode_hlt16_descriptor(
            member.descriptor_bytes, member_key
        )
        if address != member.descriptor_content_id:
            raise Pro20OriginError("captured descriptor address differs")
        arrays = _decode_hlt16_payload(member.payload_bytes, descriptor, member_key)
        return HLT16MemberSnapshot(
            arrays={name: np.asarray(arrays[name]).copy() for name in ARRAY_NAMES},
            metadata=dict(metadata),
        )


def _identity_record(
    *,
    selected_method: str,
    origin_event: int,
    origin_time: Mapping[str, str],
    generation1_bridge_content_id: str,
    generation1_checkpoint: CapturedLeaf,
    pref27_compact: CapturedLeaf,
    pref27_tree_sha256: str,
    original_leaves: tuple[CapturedLeaf, ...],
    source_inventory: SourceInventory,
    members: tuple[CapturedMemberOrigin, ...],
) -> dict[str, Any]:
    return {
        "generation1_bridge_content_id": generation1_bridge_content_id,
        "generation1_checkpoint": {
            "content_id": generation1_checkpoint.content_id,
            "path": generation1_checkpoint.path,
            "raw_sha256": generation1_checkpoint.raw_sha256,
        },
        "members": [
            {
                "CFL_retry_count": item.cfl_retry_count,
                "accepted_time": dict(item.accepted_time),
                "causal_sha256": item.causal_sha256,
                "cursor_chain_sha256": item.cursor_chain_sha256,
                "descriptor_content_id": item.descriptor_content_id,
                "descriptor_raw_sha256": _sha256_hex(item.descriptor_bytes),
                "full_cursor_sha256": item.full_cursor_sha256,
                "hlt15_state_content_id": item.hlt15_state_content_id,
                "input_hash": item.input_hash,
                "member_key": item.member_key,
                "monitor_sha256": item.monitor_sha256,
                "payload_raw_sha256": item.payload_raw_sha256,
                "physical_array_sha256": {
                    name: item.physical_array_sha256[name] for name in STATE_ARRAY_NAMES
                },
                "physical_array_shapes": {
                    name: list(item.physical_array_shapes[name])
                    for name in STATE_ARRAY_NAMES
                },
                "bridge_array_sha256": dict(item.physical_array_sha256),
                "bridge_array_shapes": {
                    name: list(shape)
                    for name, shape in item.physical_array_shapes.items()
                },
                "source_retry_count": item.source_retry_count,
                "step_index": item.step_index,
                "transaction_serial": item.transaction_serial,
            }
            for item in members
        ],
        "origin_event": origin_event,
        "origin_time": dict(origin_time),
        "pref27": {
            "compact_content_id": pref27_compact.content_id,
            "compact_path": pref27_compact.path,
            "compact_raw_sha256": pref27_compact.raw_sha256,
            "leaves": [
                {
                    "content_id": leaf.content_id,
                    "path": leaf.path,
                    "raw_sha256": leaf.raw_sha256,
                }
                for leaf in original_leaves
            ],
            "tree_sha256": pref27_tree_sha256,
        },
        "schema": CAPTURE_SCHEMA,
        "selected_method": selected_method,
        "source_inventory": {
            "byte_count": source_inventory.byte_count,
            "cleanup_digest": source_inventory.cleanup_digest,
            "leaf_count": source_inventory.leaf_count,
            "leaves": [
                {
                    "byte_count": len(leaf.raw_bytes),
                    "path": leaf.path,
                    "raw_sha256": leaf.raw_sha256,
                }
                for leaf in source_inventory.leaves
            ],
            "relative_root": source_inventory.relative_root,
        },
    }


def _pref27_leaf_records(compact: Mapping[str, Any]) -> tuple[dict[str, str], ...]:
    if compact.get("artifact_id") != "FGC-1-PRO18-PREF27":
        raise Pro20OriginError("PREF27 compact identity differs")
    payload = compact.get("artifact_payload")
    if not isinstance(payload, Mapping):
        raise Pro20OriginError("PREF27 compact payload differs")
    rows = payload.get("store_leaves")
    tree = payload.get("store_tree_sha256")
    if not isinstance(rows, list) or len(rows) != 8:
        raise Pro20OriginError("PREF27 compact does not bind eight original leaves")
    leaves: list[dict[str, str]] = []
    for row in rows:
        item = _mapping(row, "PREF27 leaf")
        path = item.get("path")
        if not isinstance(path, str) or not path.startswith(SOURCE_STORE_RELATIVE + "/"):
            raise Pro20OriginError("PREF27 leaf path differs")
        leaves.append(
            {
                "path": path,
                "raw_sha256": _require_sha(item.get("raw_sha256"), "PREF27 raw hash"),
                "content_sha256": _require_sha(
                    item.get("content_sha256"), "PREF27 content ID"
                ),
            }
        )
    ordered = tuple(sorted(leaves, key=lambda item: item["path"]))
    observed = _sha256_hex(canonical_json_bytes(list(ordered)))
    if observed != _require_sha(tree, "PREF27 tree digest"):
        raise Pro20OriginError("PREF27 eight-leaf projection digest differs")
    return ordered


def _authenticate_pref27_original_leaves(
    compact_raw: bytes, inventory: SourceInventory
) -> tuple[CapturedLeaf, dict[str, Any], tuple[CapturedLeaf, ...]]:
    compact = _load_json_object(compact_raw, "PREF27 compact", canonical=False)
    records = _pref27_leaf_records(compact)
    inventory_map = inventory.by_path()
    captured: list[CapturedLeaf] = []
    parsed: dict[str, dict[str, Any]] = {}
    for record in records:
        try:
            inventory_leaf = inventory_map[record["path"]]
        except KeyError as error:
            raise Pro20OriginError(
                f"expanded source is missing original PREF27 leaf: {record['path']}"
            ) from error
        raw = inventory_leaf.raw_bytes
        if _sha256_hex(raw) != record["raw_sha256"]:
            raise Pro20OriginError(f"PREF27 original leaf raw hash differs: {record['path']}")
        payload = _load_json_object(raw, record["path"], canonical=True)
        content_id = _content_id_for_leaf(record["path"], payload, raw)
        if content_id != record["content_sha256"]:
            raise Pro20OriginError(
                f"PREF27 original leaf content ID differs: {record['path']}"
            )
        leaf = CapturedLeaf(
            path=record["path"],
            raw_sha256=record["raw_sha256"],
            content_id=content_id,
            raw_bytes=raw,
        )
        captured.append(leaf)
        parsed[record["path"]] = payload
    compact_leaf = CapturedLeaf(
        path=PREF27_COMPACT_RELATIVE,
        raw_sha256=_sha256_hex(compact_raw),
        content_id=_sha256_hex(canonical_json_bytes(compact)),
        raw_bytes=compact_raw,
    )
    return compact_leaf, parsed, tuple(captured)


def _original_hlt15_objects(
    leaves: tuple[CapturedLeaf, ...], parsed: Mapping[str, dict[str, Any]]
) -> tuple[CapturedLeaf, CapturedLeaf, dict[str, Any], dict[str, CapturedLeaf]]:
    checkpoint_leaf = [item for item in leaves if _is_checkpoint_path(item.path)]
    receipt_leaf = [item for item in leaves if _is_receipt_path(item.path)]
    state_leaves = [item for item in leaves if _is_state_path(item.path)]
    if len(checkpoint_leaf) != 1 or len(receipt_leaf) != 1 or len(state_leaves) != 6:
        raise Pro20OriginError("PREF27 original leaf grammar differs")
    checkpoint = parsed[checkpoint_leaf[0].path]
    if checkpoint.get("committed_common_event_index") not in {None, ORIGIN_EVENT}:
        raise Pro20OriginError("HLT15 checkpoint is not common event 23")
    states = checkpoint.get("states")
    cursors = checkpoint.get("cursors")
    if not isinstance(states, Mapping) or not isinstance(cursors, Mapping):
        raise Pro20OriginError("HLT15 checkpoint member maps differ")
    if set(states) != set(MEMBER_KEYS) or set(cursors) != set(MEMBER_KEYS):
        raise Pro20OriginError("HLT15 checkpoint member maps differ")
    by_content = {item.content_id: item for item in state_leaves}
    if len(by_content) != 6:
        raise Pro20OriginError("HLT15 original state content IDs collide")
    member_leaves: dict[str, CapturedLeaf] = {}
    for key in MEMBER_KEYS:
        state_id = _require_sha(states[key], f"{key} HLT15 state")
        if state_id not in by_content:
            raise Pro20OriginError(f"{key} original HLT15 state leaf is absent")
        member_leaves[key] = by_content[state_id]
    return checkpoint_leaf[0], receipt_leaf[0], checkpoint, member_leaves


def _nonfollowing_relatives(root: Path, relative: str) -> tuple[str, ...]:
    base = root / relative
    try:
        info = base.lstat()
    except OSError as error:
        raise Pro20OriginError("historical source tree is absent") from error
    if stat.S_ISLNK(info.st_mode) or not stat.S_ISDIR(info.st_mode):
        raise Pro20OriginError("historical source tree is unsafe")
    pending = [relative]
    files: list[str] = []
    while pending:
        current_relative = pending.pop()
        current = root / current_relative
        try:
            with os.scandir(current) as iterator:
                children = sorted(iterator, key=lambda item: item.name)
        except OSError as error:
            raise Pro20OriginError("historical source tree cannot be scanned") from error
        for child in children:
            child_relative = f"{current_relative}/{child.name}"
            try:
                child_info = child.stat(follow_symlinks=False)
            except OSError as error:
                raise Pro20OriginError(
                    f"historical source path is unsafe: {child_relative}"
                ) from error
            if stat.S_ISLNK(child_info.st_mode):
                raise Pro20OriginError(
                    f"historical source path is unsafe: {child_relative}"
                )
            if stat.S_ISDIR(child_info.st_mode):
                pending.append(child_relative)
                continue
            if not stat.S_ISREG(child_info.st_mode) or child_info.st_nlink != 1:
                raise Pro20OriginError(
                    f"historical source path is unsafe: {child_relative}"
                )
            files.append(child_relative)
    return tuple(sorted(files))


def inventory_source_tree(root: Path, relative: str = SOURCE_TREE_RELATIVE) -> SourceInventory:
    """Safe no-follow inventory of the expanded historical source tree."""

    root = _require_canonical_root(root)
    store_prefix = SOURCE_STORE_RELATIVE + "/"
    paths = _nonfollowing_relatives(root, relative)
    if any(not path.startswith(store_prefix) for path in paths):
        raise Pro20OriginError("historical source tree is forked or has extra paths")
    leaves: list[CapturedLeaf] = []
    for path in paths:
        raw = _read_bytes(root, path, "historical source leaf")
        leaves.append(
            CapturedLeaf(
                path=path,
                raw_sha256=_sha256_hex(raw),
                content_id=_sha256_hex(raw),
                raw_bytes=raw,
            )
        )
    records = tuple((item.path, len(item.raw_bytes), item.raw_sha256) for item in leaves)
    return SourceInventory(
        relative_root=relative,
        leaves=tuple(leaves),
        byte_count=sum(len(item.raw_bytes) for item in leaves),
        cleanup_digest=cleanup_baseline_tree_digest(records),
    )


def _inventories_equal(first: SourceInventory, second: SourceInventory) -> bool:
    if (
        first.cleanup_digest != second.cleanup_digest
        or first.byte_count != second.byte_count
        or first.leaf_count != second.leaf_count
    ):
        return False
    left = first.by_path()
    right = second.by_path()
    if set(left) != set(right):
        return False
    return all(left[path].raw_bytes == right[path].raw_bytes for path in left)


def _accepted_time(value: object, label: str) -> dict[str, str]:
    mapping = _mapping(value, label)
    rational = mapping.get("rational")
    binary = mapping.get("binary64_hex")
    if rational != ORIGIN_TIME_RATIONAL or binary != ORIGIN_TIME_HEX:
        raise Pro20OriginError(f"{label} is not t=23/16")
    return {"rational": ORIGIN_TIME_RATIONAL, "binary64_hex": ORIGIN_TIME_HEX}


def _hlt15_state_object(raw: bytes, key: str) -> dict[str, Any]:
    state = _load_json_object(raw, f"{key} HLT15 state", canonical=True)
    if state.get("member_key") != key:
        raise Pro20OriginError(f"{key} original HLT15 state member differs")
    arrays = state.get("physical_arrays")
    if not isinstance(arrays, list) or len(arrays) != len(STATE_ARRAY_NAMES):
        raise Pro20OriginError(f"{key} original physical array inventory differs")
    for expected, item in zip(STATE_ARRAY_NAMES, arrays, strict=True):
        descriptor = _mapping(item, f"{key} {expected} array descriptor")
        if descriptor.get("logical_name") != expected:
            raise Pro20OriginError(f"{key} original physical array order differs")
        _require_sha(
            descriptor.get("little_endian_c_bytes_sha256"),
            f"{key} {expected} physical-array hash",
        )
    return state


def _decode_hlt16_descriptor(raw: bytes, key: str) -> tuple[dict[str, Any], dict[str, Any], str]:
    descriptor = _load_json_object(raw, f"{key} HLT16 descriptor", canonical=True)
    if descriptor.get("schema") != HLT16_DESCRIPTOR_SCHEMA:
        raise Pro20OriginError(f"{key} HLT16 descriptor schema differs")
    bare = dict(descriptor)
    observed = bare.pop("descriptor_sha256", None)
    address = hlt16_digest(bare)
    if observed != address:
        raise Pro20OriginError(f"{key} HLT16 descriptor address differs")
    metadata = _mapping(descriptor.get("metadata"), f"{key} HLT16 metadata")
    if _HLT16_CODEC_KEYS <= set(metadata):
        try:
            metadata = validate_codec_metadata(metadata)
        except Exception as error:
            raise Pro20OriginError(f"{key} HLT16 codec metadata differs") from error
    return descriptor, metadata, address


def _decode_hlt16_payload(
    payload_bytes: bytes, descriptor: Mapping[str, Any], key: str
) -> dict[str, np.ndarray]:
    try:
        arrays, semantic, _manifest = decode_payload(
            payload_bytes,
            expected_semantic_sha256=descriptor.get("semantic_sha256"),
            expected_manifest=descriptor.get("arrays"),
        )
    except HLT16StateStoreError as error:
        raise Pro20OriginError(f"{key} HLT16 payload cannot be decoded") from error
    if semantic != descriptor.get("semantic_sha256"):
        raise Pro20OriginError(f"{key} HLT16 payload semantic identity differs")
    if _sha256_hex(payload_bytes) != descriptor.get("raw_archive_sha256"):
        raise Pro20OriginError(f"{key} HLT16 payload raw hash differs")
    return arrays


def _loaded_metadata(loaded: Any, key: str) -> dict[str, Any]:
    descriptor = getattr(loaded, "descriptor", None)
    if isinstance(descriptor, Mapping) and isinstance(descriptor.get("metadata"), Mapping):
        return dict(descriptor["metadata"])
    metadata = getattr(loaded, "metadata", None)
    if isinstance(metadata, Mapping):
        return dict(metadata)
    raise Pro20OriginError(f"{key} loaded bridge metadata differs")


def _loaded_arrays(loaded: Any, key: str) -> Mapping[str, Any]:
    arrays = getattr(loaded, "arrays", None)
    if not isinstance(arrays, Mapping) or any(name not in arrays for name in ARRAY_NAMES):
        raise Pro20OriginError(f"{key} loaded physical arrays differ")
    extra = set(arrays) - set(ARRAY_NAMES)
    if extra:
        raise Pro20OriginError(f"{key} loaded array inventory has extra names")
    return arrays


def compare_bridge_member_to_original(
    member_key: str,
    *,
    loaded: Any,
    original_state: Mapping[str, Any],
    original_cursor: Mapping[str, Any],
    original_state_content_id: str,
    descriptor_bytes: bytes,
    payload_bytes: bytes,
) -> CapturedMemberOrigin:
    """Compare one loaded generation-1 member against its original HLT15 object."""

    if member_key not in MEMBER_KEYS:
        raise Pro20OriginError("member key differs")
    metadata = _loaded_metadata(loaded, member_key)
    runtime = _mapping(metadata.get("runtime_identity"), f"{member_key} runtime identity")
    if runtime.get("member_key") != member_key:
        raise Pro20OriginError(f"{member_key} loaded member identity differs")
    if runtime.get("generation") != 0:
        raise Pro20OriginError(f"{member_key} loaded generation is not the import bridge")
    provenance = _mapping(
        metadata.get("provenance_cursor"), f"{member_key} provenance cursor"
    )
    if (
        provenance.get("source_protocol") != GEN0_PROTOCOL
        or provenance.get("bridge") not in GEN0_BRIDGES
    ):
        raise Pro20OriginError(f"{member_key} PROTO17 predecessor provenance differs")
    original_cursor = _mapping(original_cursor, f"{member_key} original cursor")
    loaded_cursor = _mapping(provenance.get("payload"), f"{member_key} loaded cursor")
    cursor_bytes = full_cursor_bytes(original_cursor)
    if canonical_json_bytes(loaded_cursor) != cursor_bytes:
        raise Pro20OriginError(f"{member_key} original cursor differs")
    chain = _require_sha(
        original_cursor.get("cursor_chain_sha256"), f"{member_key} cursor chain"
    )
    accepted = _accepted_time(
        original_state.get("accepted_boundary_time"), f"{member_key} HLT15 time"
    )
    runtime_time = runtime.get("accepted_time")
    if isinstance(runtime_time, Mapping):
        _accepted_time(runtime_time, f"{member_key} loaded time")
    elif not isinstance(runtime_time, str) or runtime_time != ORIGIN_TIME_HEX:
        raise Pro20OriginError(f"{member_key} loaded time is not t=23/16")
    counters = _mapping(metadata.get("counters"), f"{member_key} loaded counters")
    for field, original_field in (
        ("step_index", "step_index"),
        ("transaction_serial", "transaction_serial"),
        ("source_retry_count", "source_retry_count"),
        ("CFL_retry_count", "CFL_retry_count"),
    ):
        if counters.get(field) != original_state.get(original_field):
            raise Pro20OriginError(f"{member_key} inherited {field} differs")
    template = metadata.get("template")
    input_hash = (
        template.get("input_hash")
        if isinstance(template, Mapping)
        else metadata.get("input_hash")
    )
    original_input = original_state.get("input_hash")
    if input_hash != original_input:
        raise Pro20OriginError(f"{member_key} input hash differs")
    _require_sha(input_hash, f"{member_key} input hash")
    monitor = metadata.get("monitor")
    causal = metadata.get("causal")
    if monitor != original_state.get("runtime_monitor_state"):
        raise Pro20OriginError(f"{member_key} monitor state differs")
    if causal != original_state.get("causal_state"):
        raise Pro20OriginError(f"{member_key} causal state differs")
    monitor_bytes = canonical_json_bytes(_mapping(monitor, f"{member_key} monitor"))
    causal_bytes = canonical_json_bytes(_mapping(causal, f"{member_key} causal"))
    loaded_arrays = _loaded_arrays(loaded, member_key)
    original_arrays = {
        item["logical_name"]: item for item in original_state["physical_arrays"]
    }
    captured_arrays: dict[str, bytes] = {}
    captured_hashes: dict[str, str] = {}
    captured_shapes: dict[str, tuple[int, ...]] = {}
    for name in ARRAY_NAMES:
        array = loaded_arrays[name]
        raw = physical_array_bytes(array)
        digest = _sha256_hex(raw)
        if name in STATE_ARRAY_NAMES:
            expected = original_arrays[name]["little_endian_c_bytes_sha256"]
            if digest != expected:
                raise Pro20OriginError(f"{member_key} {name} physical-array hash differs")
        captured_arrays[name] = raw
        captured_hashes[name] = digest
        captured_shapes[name] = tuple(int(dimension) for dimension in np.asarray(array).shape)
    descriptor = _require_sha(
        getattr(loaded, "descriptor_sha256", None), f"{member_key} descriptor"
    )
    payload_raw = _require_sha(
        getattr(loaded, "raw_archive_sha256", None), f"{member_key} payload raw hash"
    )
    if _sha256_hex(payload_bytes) != payload_raw:
        raise Pro20OriginError(f"{member_key} retained payload bytes differ")
    return CapturedMemberOrigin(
        member_key=member_key,
        descriptor_content_id=descriptor,
        payload_raw_sha256=payload_raw,
        hlt15_state_content_id=_require_sha(
            original_state_content_id, f"{member_key} HLT15 state"
        ),
        descriptor_bytes=bytes(descriptor_bytes),
        payload_bytes=bytes(payload_bytes),
        physical_array_bytes=captured_arrays,
        physical_array_sha256=captured_hashes,
        physical_array_shapes=captured_shapes,
        original_cursor_bytes=cursor_bytes,
        full_cursor_sha256=_sha256_hex(cursor_bytes),
        cursor_chain_sha256=chain,
        input_hash=str(input_hash),
        monitor_bytes=monitor_bytes,
        monitor_sha256=_sha256_hex(monitor_bytes),
        causal_bytes=causal_bytes,
        causal_sha256=_sha256_hex(causal_bytes),
        step_index=int(counters["step_index"]),
        transaction_serial=int(counters["transaction_serial"]),
        source_retry_count=int(counters["source_retry_count"]),
        cfl_retry_count=int(counters["CFL_retry_count"]),
        accepted_time=accepted,
    )


def _open_historical_store(source_root: Path) -> HLT16CampaignStore:
    return HLT16CampaignStore(source_root)


def _select_generation_one_bridge(store: Any) -> Any:
    try:
        checkpoint = store.authenticated_checkpoint_at_generation(SELECTED_GENERATION)
    except HLT16CampaignStoreError as error:
        raise Pro20OriginError(
            "expanded historical store is not a closed authenticated tree"
        ) from error
    generation = getattr(checkpoint, "generation", None)
    if generation != SELECTED_GENERATION:
        raise Pro20OriginError("selected checkpoint is not generation 1")
    event = getattr(checkpoint, "event", None)
    if event != ORIGIN_EVENT:
        raise Pro20OriginError("selected checkpoint is not common event 23")
    return checkpoint


def _generation1_leaf(inventory: SourceInventory, content_id: str) -> CapturedLeaf:
    suffix = f"/checkpoints/{SELECTED_GENERATION:020d}-{content_id}.json"
    matches = [leaf for leaf in inventory.leaves if leaf.path.endswith(suffix)]
    if len(matches) != 1:
        raise Pro20OriginError("generation-1 checkpoint leaf is absent")
    return matches[0]


def _generation1_descriptors(raw: bytes, path: str) -> tuple[str, dict[str, str]]:
    payload = _load_json_object(raw, "generation-1 checkpoint", canonical=True)
    content_id = _content_id_for_leaf(path, payload, raw)
    if payload.get("generation") != SELECTED_GENERATION:
        raise Pro20OriginError("selected checkpoint is not generation 1")
    if payload.get("event") != ORIGIN_EVENT:
        raise Pro20OriginError("selected checkpoint is not common event 23")
    members = _mapping(payload.get("members"), "generation-1 members")
    if set(members) != set(MEMBER_KEYS):
        raise Pro20OriginError("generation-1 members are not the sealed six-member set")
    descriptors: dict[str, str] = {}
    complete = True
    for key in MEMBER_KEYS:
        item = _mapping(members[key], f"{key} generation-1 member")
        descriptors[key] = _require_sha(item.get("descriptor_sha256"), f"{key} descriptor")
        if not {"cursor", "ledger", "descriptor_sha256"}.issubset(item):
            complete = False
    if complete:
        try:
            checked = CampaignCheckpoint.from_mapping(payload)
        except (HLT16CampaignSchemaError, TypeError, ValueError) as error:
            raise Pro20OriginError("generation-1 checkpoint schema differs") from error
        if checked.sha256 != content_id:
            raise Pro20OriginError("generation-1 checkpoint hash differs")
    return content_id, descriptors


def _require_owner_pins(
    *,
    compact_leaf: CapturedLeaf,
    compact_tree: str,
    checkpoint_leaf: CapturedLeaf,
    receipt_leaf: CapturedLeaf,
    inventory: SourceInventory,
    bridge_content_id: str,
    members: tuple[CapturedMemberOrigin, ...],
) -> None:
    pins = _OWNER_PINS
    if pins.selected_method != SELECTED_METHOD:
        raise Pro20OriginError("owner method pin differs")
    if compact_leaf.raw_sha256 != pins.pref27_compact_raw_sha256:
        raise Pro20OriginError("PREF27 compact raw hash differs from owner")
    if compact_tree != pins.pref27_tree_sha256:
        raise Pro20OriginError("PREF27 eight-leaf projection digest differs from owner")
    if (
        checkpoint_leaf.content_id != pins.hlt15_checkpoint_content_id
        or checkpoint_leaf.raw_sha256 != pins.hlt15_checkpoint_raw_sha256
    ):
        raise Pro20OriginError("HLT15 checkpoint identities differ from owner")
    if (
        receipt_leaf.content_id != pins.hlt15_receipt_content_id
        or receipt_leaf.raw_sha256 != pins.hlt15_receipt_raw_sha256
    ):
        raise Pro20OriginError("HLT15 receipt identities differ from owner")
    if (
        inventory.leaf_count != pins.expanded_leaf_count
        or inventory.byte_count != pins.expanded_byte_count
        or inventory.cleanup_digest != pins.expanded_cleanup_digest
    ):
        raise Pro20OriginError("expanded historical source inventory differs from owner")
    if bridge_content_id != pins.generation1_bridge_content_id:
        raise Pro20OriginError("generation-1 bridge content ID differs from owner")
    for member in members:
        expected = pins.member_counters[member.member_key]
        if (
            member.step_index != expected["step_index"]
            or member.transaction_serial != expected["transaction_serial"]
            or member.source_retry_count != expected["source_retry_count"]
            or member.cfl_retry_count != expected["CFL_retry_count"]
        ):
            raise Pro20OriginError(
                f"{member.member_key} inherited counters differ from owner"
            )


def _reconstruct_member_from_bytes(
    key: str,
    *,
    member: CapturedMemberOrigin,
    original_state: Mapping[str, Any],
    original_cursor: Mapping[str, Any],
    original_state_content_id: str,
) -> CapturedMemberOrigin:
    descriptor, metadata, address = _decode_hlt16_descriptor(member.descriptor_bytes, key)
    arrays = _decode_hlt16_payload(member.payload_bytes, descriptor, key)
    loaded = SimpleNamespace(
        arrays=arrays,
        descriptor={**descriptor, "metadata": metadata},
        metadata=metadata,
        descriptor_sha256=address,
        raw_archive_sha256=_sha256_hex(member.payload_bytes),
    )
    return compare_bridge_member_to_original(
        key,
        loaded=loaded,
        original_state=original_state,
        original_cursor=original_cursor,
        original_state_content_id=original_state_content_id,
        descriptor_bytes=member.descriptor_bytes,
        payload_bytes=member.payload_bytes,
    )


def _identity_from_retained_bytes(capture: Pro20OriginCapture) -> dict[str, Any]:
    compact_leaf, parsed, original_leaves = _authenticate_pref27_original_leaves(
        capture.pref27_compact.raw_bytes, capture.source_inventory
    )
    if compact_leaf.raw_bytes != capture.pref27_compact.raw_bytes:
        raise Pro20OriginError("PREF27 compact identity drifted")
    checkpoint_leaf, receipt_leaf, checkpoint, member_leaves = _original_hlt15_objects(
        original_leaves, parsed
    )
    bridge_id, descriptors = _generation1_descriptors(
        capture.generation1_checkpoint.raw_bytes,
        capture.generation1_checkpoint.path,
    )
    source_checkpoint = _generation1_leaf(capture.source_inventory, bridge_id)
    if (
        source_checkpoint.path != capture.generation1_checkpoint.path
        or source_checkpoint.raw_sha256 != capture.generation1_checkpoint.raw_sha256
        or source_checkpoint.raw_bytes != capture.generation1_checkpoint.raw_bytes
    ):
        raise Pro20OriginError(
            "retained generation-1 checkpoint is not its source-inventory leaf"
        )
    reconstructed: list[CapturedMemberOrigin] = []
    if tuple(item.member_key for item in capture.members) != MEMBER_KEYS:
        raise Pro20OriginError("captured members are not in canonical six-member order")
    for expected, member in zip(MEMBER_KEYS, capture.members, strict=True):
        if member.member_key != expected:
            raise Pro20OriginError("captured member order differs")
        if descriptors[expected] != _decode_hlt16_descriptor(member.descriptor_bytes, expected)[2]:
            raise Pro20OriginError(f"{expected} retained descriptor is not the generation-1 member")
        original_state = _hlt15_state_object(member_leaves[expected].raw_bytes, expected)
        reconstructed.append(
            _reconstruct_member_from_bytes(
                expected,
                member=member,
                original_state=original_state,
                original_cursor=checkpoint["cursors"][expected],
                original_state_content_id=member_leaves[expected].content_id,
            )
        )
    compact_payload = _load_json_object(
        capture.pref27_compact.raw_bytes, "PREF27 compact", canonical=False
    )["artifact_payload"]
    compact_tree = _require_sha(
        compact_payload.get("store_tree_sha256"), "PREF27 tree digest"
    )
    _require_owner_pins(
        compact_leaf=compact_leaf,
        compact_tree=compact_tree,
        checkpoint_leaf=checkpoint_leaf,
        receipt_leaf=receipt_leaf,
        inventory=capture.source_inventory,
        bridge_content_id=bridge_id,
        members=tuple(reconstructed),
    )
    return _identity_record(
        selected_method=SELECTED_METHOD,
        origin_event=ORIGIN_EVENT,
        origin_time=dict(ORIGIN_TIME),
        generation1_bridge_content_id=bridge_id,
        generation1_checkpoint=capture.generation1_checkpoint,
        pref27_compact=compact_leaf,
        pref27_tree_sha256=compact_tree,
        original_leaves=original_leaves,
        source_inventory=capture.source_inventory,
        members=tuple(reconstructed),
    )


def revalidate_pro20_origin_capture(capture: Pro20OriginCapture) -> Pro20OriginCapture:
    """Rebuild every derived claim from retained bytes; ignore success labels."""

    if not isinstance(capture, Pro20OriginCapture):
        raise Pro20OriginError("origin capture type differs")
    reconstructed = _identity_from_retained_bytes(capture)
    if canonical_json_bytes(reconstructed) != capture.captured_bytes:
        raise Pro20OriginError("captured labels are not bound to retained origin bytes")
    if reconstructed != capture.identity_record():
        raise Pro20OriginError("captured identity record is not bound to retained bytes")
    return capture


def capture_pro20_historical_origin(repository_root: Path | str) -> Pro20OriginCapture:
    """Capture the verified HLT15 t=23/16 origin as an immutable input slice."""

    root = _require_canonical_root(repository_root)
    compact_before = _read_bytes(root, PREF27_COMPACT_RELATIVE, "PREF27 compact")
    inventory_before = inventory_source_tree(root, SOURCE_TREE_RELATIVE)
    compact_leaf, parsed_leaves, original_leaves = _authenticate_pref27_original_leaves(
        compact_before, inventory_before
    )
    checkpoint_leaf, receipt_leaf, checkpoint, member_leaves = _original_hlt15_objects(
        original_leaves, parsed_leaves
    )
    compact_payload = _load_json_object(
        compact_before, "PREF27 compact", canonical=False
    )["artifact_payload"]
    compact_tree = _require_sha(
        compact_payload.get("store_tree_sha256"), "PREF27 tree digest"
    )
    store = _open_historical_store(root / SOURCE_STORE_RELATIVE)
    checkpoint_bridge = _select_generation_one_bridge(store)
    store_bridge_id = _require_sha(
        getattr(checkpoint_bridge, "sha256", None), "generation-1 bridge"
    )
    inventory_gen1 = _generation1_leaf(inventory_before, store_bridge_id)
    bridge_content_id, gen1_descriptors = _generation1_descriptors(
        inventory_gen1.raw_bytes, inventory_gen1.path
    )
    generation1_leaf = CapturedLeaf(
        path=inventory_gen1.path,
        raw_sha256=inventory_gen1.raw_sha256,
        content_id=bridge_content_id,
        raw_bytes=inventory_gen1.raw_bytes,
    )
    if bridge_content_id != store_bridge_id:
        raise Pro20OriginError("generation-1 store selection differs from retained checkpoint")
    members_map = getattr(checkpoint_bridge, "members", None)
    if not isinstance(members_map, Mapping) or set(members_map) != set(MEMBER_KEYS):
        raise Pro20OriginError("generation-1 members are not the sealed six-member set")
    inventory_map = inventory_before.by_path()
    captured_members: list[CapturedMemberOrigin] = []
    for key in MEMBER_KEYS:
        descriptor = getattr(members_map[key], "descriptor_sha256", None)
        if descriptor != gen1_descriptors[key]:
            raise Pro20OriginError(f"{key} generation-1 descriptor address differs")
        try:
            loaded = store.load_state(descriptor)
        except HLT16CampaignStoreError as error:
            raise Pro20OriginError(f"{key} generation-1 state cannot be loaded") from error
        descriptor_path = f"{SOURCE_STORE_RELATIVE}/states/{descriptor}.json"
        try:
            descriptor_leaf = inventory_map[descriptor_path]
        except KeyError as error:
            raise Pro20OriginError(f"{key} generation-1 descriptor is not in the source inventory") from error
        descriptor_obj, metadata, address = _decode_hlt16_descriptor(
            descriptor_leaf.raw_bytes, key
        )
        if address != descriptor:
            raise Pro20OriginError(f"{key} retained descriptor address differs")
        semantic = descriptor_obj.get("semantic_sha256")
        payload_path = f"{SOURCE_STORE_RELATIVE}/payloads/{semantic}.npz"
        try:
            payload_leaf = inventory_map[payload_path]
        except KeyError as error:
            raise Pro20OriginError(f"{key} generation-1 payload is not in the source inventory") from error
        decoded_arrays = _decode_hlt16_payload(payload_leaf.raw_bytes, descriptor_obj, key)
        loaded_arrays = _loaded_arrays(loaded, key)
        for name in ARRAY_NAMES:
            if physical_array_bytes(loaded_arrays[name]) != physical_array_bytes(decoded_arrays[name]):
                raise Pro20OriginError(f"{key} loaded {name} differs from retained payload")
        loaded = SimpleNamespace(
            arrays=decoded_arrays,
            descriptor={**descriptor_obj, "metadata": metadata},
            metadata=metadata,
            descriptor_sha256=address,
            raw_archive_sha256=_sha256_hex(payload_leaf.raw_bytes),
        )
        original_state = _hlt15_state_object(member_leaves[key].raw_bytes, key)
        captured_members.append(
            compare_bridge_member_to_original(
                key,
                loaded=loaded,
                original_state=original_state,
                original_cursor=checkpoint["cursors"][key],
                original_state_content_id=member_leaves[key].content_id,
                descriptor_bytes=descriptor_leaf.raw_bytes,
                payload_bytes=payload_leaf.raw_bytes,
            )
        )
    compact_after = _read_bytes(root, PREF27_COMPACT_RELATIVE, "PREF27 compact")
    inventory_after = inventory_source_tree(root, SOURCE_TREE_RELATIVE)
    if compact_after != compact_before:
        raise Pro20OriginError("PREF27 compact drifted during capture")
    if not _inventories_equal(inventory_before, inventory_after):
        raise Pro20OriginError("historical source inventory drifted during capture")
    member_tuple = tuple(captured_members)
    _require_owner_pins(
        compact_leaf=compact_leaf,
        compact_tree=compact_tree,
        checkpoint_leaf=checkpoint_leaf,
        receipt_leaf=receipt_leaf,
        inventory=inventory_before,
        bridge_content_id=bridge_content_id,
        members=member_tuple,
    )
    identity = _identity_record(
        selected_method=SELECTED_METHOD,
        origin_event=ORIGIN_EVENT,
        origin_time=dict(ORIGIN_TIME),
        generation1_bridge_content_id=bridge_content_id,
        generation1_checkpoint=generation1_leaf,
        pref27_compact=compact_leaf,
        pref27_tree_sha256=compact_tree,
        original_leaves=original_leaves,
        source_inventory=inventory_before,
        members=member_tuple,
    )
    capture = Pro20OriginCapture(
        selected_method=SELECTED_METHOD,
        origin_event=ORIGIN_EVENT,
        origin_time=dict(ORIGIN_TIME),
        generation1_bridge_content_id=bridge_content_id,
        generation1_checkpoint=generation1_leaf,
        pref27_compact=compact_leaf,
        pref27_tree_sha256=compact_tree,
        original_leaves=original_leaves,
        source_inventory=inventory_before,
        members=member_tuple,
        captured_bytes=canonical_json_bytes(identity),
    )
    return revalidate_pro20_origin_capture(capture)


__all__ = [
    "CAPTURE_SCHEMA",
    "CLEANUP_DOMAIN",
    "CapturedLeaf",
    "CapturedMemberOrigin",
    "ORIGIN_EVENT",
    "ORIGIN_TIME",
    "ORIGIN_TIME_HEX",
    "ORIGIN_TIME_RATIONAL",
    "PREF27_COMPACT_RELATIVE",
    "Pro20OriginCapture",
    "Pro20OriginError",
    "SELECTED_GENERATION",
    "SELECTED_METHOD",
    "SOURCE_STORE_RELATIVE",
    "SOURCE_TREE_RELATIVE",
    "SourceInventory",
    "capture_pro20_historical_origin",
    "cleanup_baseline_tree_digest",
    "compare_bridge_member_to_original",
    "full_cursor_sha256",
    "inventory_source_tree",
    "physical_array_sha256",
    "revalidate_pro20_origin_capture",
]
