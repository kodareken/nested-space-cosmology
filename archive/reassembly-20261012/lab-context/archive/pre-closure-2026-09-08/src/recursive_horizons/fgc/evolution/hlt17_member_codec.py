"""Pure in-memory HLT17 member codec: canonical descriptor plus binary payload.

The wire object is not a store writer, journal, or scientific gate.  Descriptor
SHA-256 hashes the canonical JSON bytes.  Physical state SHA-256 hashes the
``u``, ``p``, and ``q`` arrays and is checked separately.  Restore always
reauthenticates those original encoding bytes; a mutable decoded snapshot is
never a restore authority.

Accepted physical-state generation is not a cursor/checkpoint/journal
revision.  Accepted descriptors are either the origin same-boundary or a
fine successor (+4 substeps and 20 RK4 / 16 SSPRK3 records) whose executed
TDG7 plan current/endpoint match predecessor time and accepted time.
Retry-ledger overlays are a distinct snapshot role on that same accepted
boundary; they do not advance state or time and do not relax successor
checks.  Attempted retry plans remain cursor-owned.

The sealed HLT15 origin protocol is ``FGC-2-SF1-PROTO17``.  Synthetic member
scope stays explicit.  A successor cursor is forbidden because it would make
the hash circular.  ``predecessor_cursor_identity_sha256`` hashes the small
identity projection; ``predecessor_cursor_sha256`` is a supplied full-cursor
digest reference.  Matching either digest is not authentication.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from hashlib import sha256
import math
import re
import struct
from types import MappingProxyType
from typing import Any, Callable, Mapping

import numpy as np

from recursive_horizons.evidence_io import (
    CanonicalJSONError,
    canonical_json_bytes,
    load_canonical_json,
)

from .boundary_domain import CausalBudgetState
from .hlt17_admission_runtime import (
    IMPLEMENTATION_IDENTITY_KEYS,
    hlt17_implementation_identity,
    require_hlt17_implementation_identity,
)
from .hlt17_runtime_member import (
    HLT17_LIVE_MEMBER_KEYS,
    HLT17MemberCapture,
    HLT17RuntimeMember,
    HLT17RuntimeMemberError,
    IDENTITY_SCOPE_LIVE,
    IDENTITY_SCOPE_SYNTHETIC,
    INTEGRATOR_BY_LABEL,
    RUNTIME_PROTOCOL,
    SPATIAL_ORDER_BY_LABEL,
)
from .numerical_engine import EvolutionState, array_content_sha256
from .proto5_runtime import GR0RuntimeMonitorState
from .tdg11_imp1_ledger import (
    IMP1_ACCEPTED_SUBSTEP_COUNT,
    IMP1_METHOD_STAGE_RECORD_COUNT,
    imp1_checkpoint_extension,
    restore_imp1_checkpoint_extension,
)
from .tdg7_binary64_subdivision_lattice import (
    TDG7Binary64SubdivisionPlan,
    validate_tdg7_binary64_subdivision_plan,
)


CODEC_SCHEMA = "FGC-1-HLT17-C1R1-member-state-v1"
CODEC_SCHEMA_VERSION = 1
CODEC_LEGACY_IMP1_SCHEMA = "FGC-1-HLT17-member-state-v1"
ORIGIN_PROTOCOL = "FGC-2-SF1-PROTO17"
SNAPSHOT_ROLE_ACCEPTED = "accepted_physical_state"
SNAPSHOT_ROLE_RETRY = "retry_ledger_overlay"
SNAPSHOT_ROLES = (SNAPSHOT_ROLE_ACCEPTED, SNAPSHOT_ROLE_RETRY)
ARRAY_NAMES = (
    "u",
    "p",
    "q",
    "grid_coordinates",
    "tracer_labels",
    "tracer_positions",
    "tracer_proper_times",
    "event_proper_times",
    "event_fields",
)
PAYLOAD_MAGIC = b"HLT17ARR"
PAYLOAD_VERSION = 1
NAME_WIDTH = 32
SHAPE_SLOTS = 4
HEADER_STRUCT = struct.Struct("<8sII")
ENTRY_STRUCT = struct.Struct("<32sI4x6Q32s")
HEADER_BYTES = HEADER_STRUCT.size
ENTRY_BYTES = ENTRY_STRUCT.size
TABLE_BYTES = ENTRY_BYTES * len(ARRAY_NAMES)
DATA_OFFSET = HEADER_BYTES + TABLE_BYTES
MAX_DESCRIPTOR_BYTES = 16 * 1024 * 1024
MAX_PAYLOAD_BYTES = 96 * 1024 * 1024
MAX_ARRAY_BYTES = 32 * 1024 * 1024
MAX_LABELS = 16385
MAX_EVENTS = 65536
_SHA256 = re.compile(r"\A[0-9a-f]{64}\Z")

DESCRIPTOR_KEYS = (
    "accepted_time_hex",
    "array_inventory",
    "campaign_execution_authorized",
    "causal",
    "counters",
    "durable_store_bound",
    "executed_plan",
    "historical_origin_authenticated",
    "imp1_checkpoint",
    "monitor",
    "origin_reference",
    "payload_byte_count",
    "payload_sha256",
    "physical_state_sha256",
    "predecessor_cursor_authenticated",
    "predecessor_cursor_identity",
    "predecessor_cursor_identity_sha256",
    "predecessor_cursor_sha256",
    "predecessor_descriptor_sha256",
    "protocol_artifact_id",
    "runtime_identity",
    "snapshot_role",
    "runtime_template",
    "runtime_template_sha256",
    "schema",
    "schema_version",
    "source_operator_closure",
    "source_operator_closure_authenticated",
    "tracer",
)
RUNTIME_IDENTITY_KEYS = (
    "accepted_generation",
    "admission_runtime_id",
    "amplitude",
    "campaign_id",
    "implementation_id",
    "integrator_id",
    "mathematical_object",
    "member_key",
    "method",
    "point_count",
    "reference_wire_artifact_id",
    "reference_wire_evaluator_id",
    "scope",
    "spatial_order",
)
COUNTER_KEYS = (
    "CFL_retry_count",
    "accepted_macro_steps",
    "source_retry_count",
    "step_index",
    "transaction_serial",
)
MONITOR_KEYS = (
    "accepted_stage_count",
    "first_failed_premise",
    "first_failed_time_hex",
    "first_failed_transaction_serial",
    "last_accepted_time_hex",
    "last_transaction_serial",
)
CAUSAL_KEYS = (
    "accepted_time_hex",
    "accumulated_characteristic_distance_hex",
    "previous_speed_upper_hex",
)
TRACER_KEYS = (
    "cutoff_hex",
    "event_count",
    "label_count",
    "outer_radius_hex",
)
INVENTORY_KEYS = (
    "byte_count",
    "bytes_sha256",
    "dtype",
    "name",
    "offset",
    "order",
    "shape",
)
PREDECESSOR_CURSOR_IDENTITY_KEYS = (
    "accepted_boundary_time_hex",
    "accepted_generation",
    "accepted_state_sha256",
    "campaign_id",
    "member_key",
    "method",
    "point_count",
    "previous_step_index",
    "previous_transaction_serial",
    "protocol_artifact_id",
)
PLAN_KEYS = (
    "boundaries_hex",
    "conservative_reduction_hex",
    "current_hex",
    "event_target_hex",
    "fine_width_hex",
    "macro_width_hex",
    "minimum_width_hex",
    "quantum_hex",
    "requested_cap_hex",
    "selection_budget_hex",
    "target_limited",
)
REFERENCE_KEYS = ("authenticated", "kind", "sha256")
FORBIDDEN_DESCRIPTOR_KEYS = (
    "successor_cursor",
    "successor_cursor_sha256",
    "cursor_chain_sha256",
)

HLT17_CURSOR_BRIDGE_SEAMS = (
    "Accepted physical-state generation is independent of cursor_generation, "
    "checkpoint_generation, attempt_serial, and kernel-transition digests.  "
    "Those revisions belong to the HLT17 cursor; FUTURE store journal/"
    "checkpoint revisions are a later exclusive writer.",
    "predecessor_cursor_identity_sha256 hashes the small closed identity "
    "projection.  predecessor_cursor_sha256 is a supplied full-cursor digest "
    "reference.  Matching either digest is not cursor or origin authentication.",
    "The sealed HLT15 origin protocol is FGC-2-SF1-PROTO17.  Synthetic scope "
    "stays explicit; live keys are never inferred from nine-node arrays.",
    "A successor cursor is never embedded.  Retry plans and mixed-owner "
    "overlays remain cursor-owned; this codec only round-trips the IMP1 "
    "retry ledger on the same accepted boundary.",
    "source_operator_closure and origin_reference remain unauthenticated "
    "supplied digests until the durable authority binds them.",
    "No journal, campaign store, atomic publication, or live production "
    "cursor is implemented here.",
)


class HLT17MemberCodecError(ValueError):
    """A member encoding is malformed, noncanonical, or resource-unbounded."""


class HLT17MemberCodecResourceError(HLT17MemberCodecError):
    """A codec bound was exceeded before any member mutation."""


def _fail(message: str) -> None:
    raise HLT17MemberCodecError(message)


def _text(value: object, label: str) -> str:
    if type(value) is not str or not value:
        _fail(f"{label} must be nonempty text")
    return value


def _digest_text(value: object, label: str) -> str:
    text = _text(value, label)
    if _SHA256.fullmatch(text) is None:
        _fail(f"{label} must be a SHA-256 digest")
    return text


def _bool_false(value: object, label: str) -> bool:
    if type(value) is not bool or value is not False:
        _fail(f"{label} must be false")
    return False


def _nonnegative_int(value: object, label: str) -> int:
    if type(value) is not int or value < 0:
        _fail(f"{label} must be a nonnegative built-in integer")
    return value


def _integer(value: object, label: str, *, minimum: int) -> int:
    if type(value) is not int or value < minimum:
        _fail(f"{label} must be a built-in integer >= {minimum}")
    return value


def _finite_hex(value: object, label: str, *, positive: bool = False) -> str:
    if type(value) is not str:
        _fail(f"{label} is not canonical binary64 text")
    try:
        number = float.fromhex(value)
    except (OverflowError, ValueError) as error:
        raise HLT17MemberCodecError(f"{label} is not canonical binary64 text") from error
    if not math.isfinite(number) or number.hex() != value:
        _fail(f"{label} is not canonical binary64 text")
    if positive and not (number > 0.0):
        _fail(f"{label} must be positive")
    return value


def _from_hex(value: object, label: str) -> float:
    return float.fromhex(_finite_hex(value, label))


def _exact_mapping(value: object, names: tuple[str, ...], label: str) -> dict[str, object]:
    if type(value) is not dict or set(value) != set(names):
        _fail(f"{label} is incomplete or broadened")
    return value


def _canonical(value: object, label: str) -> bytes:
    try:
        raw = canonical_json_bytes(value)
    except CanonicalJSONError as error:
        if "byte bound" in str(error):
            raise HLT17MemberCodecResourceError(
                f"{label} exceeds MAX_DESCRIPTOR_BYTES"
            ) from error
        raise HLT17MemberCodecError(f"{label} is not canonical JSON") from error
    if len(raw) > MAX_DESCRIPTOR_BYTES:
        raise HLT17MemberCodecResourceError(f"{label} exceeds MAX_DESCRIPTOR_BYTES")
    return raw


def _load_descriptor(raw: object) -> dict[str, object]:
    if type(raw) not in (bytes, bytearray):
        raise TypeError("descriptor must be bytes")
    payload = bytes(raw)
    if len(payload) > MAX_DESCRIPTOR_BYTES:
        raise HLT17MemberCodecResourceError("descriptor exceeds MAX_DESCRIPTOR_BYTES")
    try:
        value = load_canonical_json(payload)
    except CanonicalJSONError as error:
        raise HLT17MemberCodecError("descriptor is not canonical JSON") from error
    if type(value) is not dict:
        _fail("descriptor is not a JSON object")
    return value


def _digest_bytes(payload: bytes) -> str:
    return sha256(payload).hexdigest()


def _strict_array(value: object, name: str) -> np.ndarray:
    if not isinstance(value, np.ndarray):
        _fail(f"{name} is not a NumPy array")
    if (
        value.dtype != np.dtype("<f8")
        or not value.flags.c_contiguous
        or not value.shape
        or not np.isfinite(value).all()
    ):
        _fail(f"{name} is not finite little-endian C-order float64")
    if value.nbytes > MAX_ARRAY_BYTES:
        raise HLT17MemberCodecResourceError(f"{name} exceeds MAX_ARRAY_BYTES")
    return np.array(value, dtype=np.dtype("<f8"), order="C", copy=True)


def _shape_list(shape: tuple[int, ...]) -> list[int]:
    if any(type(item) is not int or item < 1 for item in shape):
        _fail("array shape is not a positive integer tuple")
    return [int(item) for item in shape]


def _reference_mapping(value: object, expected_kind: str) -> dict[str, object]:
    item = _exact_mapping(value, REFERENCE_KEYS, expected_kind)
    kind = _text(item["kind"], "reference kind")
    if kind != expected_kind:
        _fail(f"{expected_kind} kind differs")
    return {
        "kind": kind,
        "sha256": _digest_text(item["sha256"], f"{expected_kind} sha256"),
        "authenticated": _bool_false(item["authenticated"], f"{expected_kind} authenticated"),
    }


def encode_plan(plan: TDG7Binary64SubdivisionPlan) -> dict[str, object]:
    validate_tdg7_binary64_subdivision_plan(plan)
    minimum = plan.minimum_width
    mapping = {
        "current_hex": plan.current.hex(),
        "event_target_hex": plan.event_target.hex(),
        "requested_cap_hex": plan.requested_cap.hex(),
        "quantum_hex": plan.quantum.hex(),
        "selection_budget_hex": plan.selection_budget.hex(),
        "macro_width_hex": plan.macro_width.hex(),
        "fine_width_hex": plan.fine_width.hex(),
        "conservative_reduction_hex": plan.conservative_reduction.hex(),
        "minimum_width_hex": None if minimum is None else minimum.hex(),
        "boundaries_hex": [item.hex() for item in plan.boundaries],
        "target_limited": plan.target_limited,
    }
    if set(mapping) != set(PLAN_KEYS):
        _fail("executed TDG7 plan schema differs")
    _canonical(mapping, "executed TDG7 plan")
    return mapping


def decode_plan(value: object) -> TDG7Binary64SubdivisionPlan:
    item = _exact_mapping(value, PLAN_KEYS, "executed TDG7 plan")
    if type(item["target_limited"]) is not bool:
        _fail("TDG7 target_limited is not a boolean")
    boundaries = item["boundaries_hex"]
    if type(boundaries) is not list or len(boundaries) != 5:
        _fail("TDG7 boundaries differ")
    minimum = item["minimum_width_hex"]
    try:
        plan = TDG7Binary64SubdivisionPlan(
            current=_from_hex(item["current_hex"], "current_hex"),
            event_target=_from_hex(item["event_target_hex"], "event_target_hex"),
            requested_cap=_from_hex(item["requested_cap_hex"], "requested_cap_hex"),
            quantum=_from_hex(item["quantum_hex"], "quantum_hex"),
            selection_budget=_from_hex(item["selection_budget_hex"], "selection_budget_hex"),
            macro_width=_from_hex(item["macro_width_hex"], "macro_width_hex"),
            fine_width=_from_hex(item["fine_width_hex"], "fine_width_hex"),
            conservative_reduction=_from_hex(
                item["conservative_reduction_hex"], "conservative_reduction_hex"
            ),
            minimum_width=(
                None if minimum is None else _from_hex(minimum, "minimum_width_hex")
            ),
            boundaries=tuple(
                _from_hex(entry, f"boundaries_hex[{index}]")
                for index, entry in enumerate(boundaries)
            ),
            target_limited=item["target_limited"],
        )
        validate_tdg7_binary64_subdivision_plan(plan)
    except (TypeError, ValueError) as error:
        raise HLT17MemberCodecError("executed TDG7 plan differs") from error
    if encode_plan(plan) != item:
        _fail("executed TDG7 plan is not a closed encoding")
    return plan


def _validated_predecessor_cursor(value: object) -> dict[str, object]:
    item = _exact_mapping(value, PREDECESSOR_CURSOR_IDENTITY_KEYS, "predecessor cursor")
    point_count = _nonnegative_int(item["point_count"], "predecessor point_count")
    method = _text(item["method"], "predecessor method")
    member_key = _text(item["member_key"], "predecessor member_key")
    if member_key != f"{method}-{point_count}":
        _fail("predecessor member key disagrees with method and point count")
    checked = {
        "protocol_artifact_id": _text(
            item["protocol_artifact_id"], "predecessor protocol"
        ),
        "campaign_id": _text(item["campaign_id"], "predecessor campaign_id"),
        "member_key": member_key,
        "method": method,
        "point_count": point_count,
        "accepted_boundary_time_hex": _finite_hex(
            item["accepted_boundary_time_hex"], "predecessor accepted time"
        ),
        "accepted_state_sha256": _digest_text(
            item["accepted_state_sha256"], "predecessor accepted_state_sha256"
        ),
        "previous_step_index": _nonnegative_int(
            item["previous_step_index"], "predecessor previous_step_index"
        ),
        "previous_transaction_serial": _nonnegative_int(
            item["previous_transaction_serial"],
            "predecessor previous_transaction_serial",
        ),
        "accepted_generation": _nonnegative_int(
            item["accepted_generation"], "predecessor accepted_generation"
        ),
    }
    _canonical(checked, "predecessor cursor identity")
    return checked


@dataclass(frozen=True, slots=True)
class HLT17PredecessorIdentity:
    """Supplied predecessor references.  Not an HLT17 cursor object.

    ``identity_sha256`` hashes the small closed identity projection.
    ``cursor_sha256`` is a supplied digest of the complete predecessor cursor
    bytes.  Validating either hex string is not authentication.
    """

    descriptor_sha256: str
    cursor_identity: Mapping[str, object]
    cursor_sha256: str
    cursor_authenticated: bool = False

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "descriptor_sha256",
            _digest_text(self.descriptor_sha256, "predecessor descriptor_sha256"),
        )
        object.__setattr__(
            self,
            "cursor_identity",
            _validated_predecessor_cursor(dict(self.cursor_identity)),
        )
        object.__setattr__(
            self,
            "cursor_sha256",
            _digest_text(self.cursor_sha256, "predecessor cursor_sha256"),
        )
        object.__setattr__(
            self,
            "cursor_authenticated",
            _bool_false(self.cursor_authenticated, "predecessor cursor_authenticated"),
        )

    @property
    def identity_sha256(self) -> str:
        return _digest_bytes(
            _canonical(self.cursor_identity, "predecessor cursor identity")
        )


@dataclass(frozen=True, slots=True)
class HLT17MemberEncoding:
    """Authenticated descriptor+payload bytes. Identity hashes are derived."""

    descriptor: bytes
    payload: bytes
    _decoded: HLT17DecodedSnapshot = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        if type(self.descriptor) is not bytes or type(self.payload) is not bytes:
            _fail("encoding bytes differ")
        decoded = decode_member(self.descriptor, self.payload)
        object.__setattr__(self, "descriptor", bytes(decoded.descriptor))
        object.__setattr__(self, "payload", bytes(decoded.payload))
        object.__setattr__(self, "_decoded", decoded)

    @property
    def descriptor_sha256(self) -> str:
        return self._decoded.descriptor_sha256

    @property
    def payload_sha256(self) -> str:
        return self._decoded.payload_sha256

    @property
    def physical_state_sha256(self) -> str:
        return self._decoded.physical_state_sha256

    @property
    def decoded(self) -> HLT17DecodedSnapshot:
        return self._decoded

    @property
    def executed_plan(self) -> TDG7Binary64SubdivisionPlan | None:
        return self._decoded.executed_plan

    @property
    def metadata(self) -> Mapping[str, Any]:
        # A caller may edit its projection, but not the cached validation
        # source used by identity properties. The authoritative bytes remain
        # immutable and every projection is reconstructed from them.
        return _load_descriptor(self.descriptor)


@dataclass(frozen=True, slots=True)
class HLT17DecodedSnapshot:
    arrays: Mapping[str, np.ndarray]
    metadata: Mapping[str, Any]
    ledger: object
    executed_plan: TDG7Binary64SubdivisionPlan | None
    descriptor: bytes
    payload: bytes
    descriptor_sha256: str
    physical_state_sha256: str
    payload_sha256: str


def _inventory_and_payload(
    arrays: Mapping[str, np.ndarray],
) -> tuple[list[dict[str, object]], bytes]:
    if not isinstance(arrays, Mapping) or tuple(arrays) != ARRAY_NAMES:
        _fail("snapshot array inventory/order differs")
    checked = {name: _strict_array(arrays[name], name) for name in ARRAY_NAMES}
    blobs: list[bytes] = []
    inventory: list[dict[str, object]] = []
    offset = DATA_OFFSET
    total = 0
    for name in ARRAY_NAMES:
        array = checked[name]
        blob = array.tobytes(order="C")
        size = len(blob)
        if size != array.nbytes or size % 8 != 0:
            _fail(f"{name} byte image is not packed binary64")
        if offset % 8 != 0:
            _fail("payload offset is not binary64 aligned")
        total += size
        if total > MAX_PAYLOAD_BYTES:
            raise HLT17MemberCodecResourceError("payload exceeds MAX_PAYLOAD_BYTES")
        inventory.append(
            {
                "name": name,
                "dtype": "<f8",
                "shape": _shape_list(array.shape),
                "order": "C",
                "offset": offset,
                "byte_count": size,
                "bytes_sha256": _digest_bytes(blob),
            }
        )
        blobs.append(blob)
        offset += size
    table = bytearray()
    for item, array in zip(inventory, (checked[name] for name in ARRAY_NAMES), strict=True):
        shape = list(array.shape) + [0] * (SHAPE_SLOTS - array.ndim)
        if array.ndim > SHAPE_SLOTS:
            _fail("array rank exceeds the codec bound")
        name_bytes = str(item["name"]).encode("ascii")
        if len(name_bytes) >= NAME_WIDTH:
            _fail("array name exceeds the codec bound")
        table.extend(
            ENTRY_STRUCT.pack(
                name_bytes.ljust(NAME_WIDTH, b"\0"),
                array.ndim,
                shape[0],
                shape[1],
                shape[2],
                shape[3],
                int(item["offset"]),
                int(item["byte_count"]),
                bytes.fromhex(str(item["bytes_sha256"])),
            )
        )
    payload = (
        HEADER_STRUCT.pack(PAYLOAD_MAGIC, PAYLOAD_VERSION, len(ARRAY_NAMES))
        + bytes(table)
        + b"".join(blobs)
    )
    if len(payload) != offset:
        _fail("payload length disagrees with the inventory")
    if len(payload) > MAX_PAYLOAD_BYTES:
        raise HLT17MemberCodecResourceError("payload exceeds MAX_PAYLOAD_BYTES")
    return inventory, payload


def _parse_payload_table(payload: bytes) -> list[dict[str, object]]:
    if type(payload) is not bytes:
        raise TypeError("payload must be bytes")
    if len(payload) > MAX_PAYLOAD_BYTES:
        raise HLT17MemberCodecResourceError("payload exceeds MAX_PAYLOAD_BYTES")
    if len(payload) < DATA_OFFSET:
        _fail("payload is truncated before its array table")
    magic, version, count = HEADER_STRUCT.unpack_from(payload, 0)
    if magic != PAYLOAD_MAGIC or version != PAYLOAD_VERSION or count != len(ARRAY_NAMES):
        _fail("payload header differs")
    inventory: list[dict[str, object]] = []
    view = memoryview(payload)
    cursor = HEADER_BYTES
    for expected in ARRAY_NAMES:
        name_raw, ndim, d0, d1, d2, d3, offset, byte_count, digest = ENTRY_STRUCT.unpack_from(
            payload, cursor
        )
        name = name_raw.split(b"\0", 1)[0].decode("ascii")
        if name != expected or type(ndim) is not int or ndim < 1 or ndim > SHAPE_SLOTS:
            _fail("payload array inventory/order differs")
        dims = (d0, d1, d2, d3)[:ndim]
        unused = (d0, d1, d2, d3)[ndim:]
        if any(item != 0 for item in unused):
            _fail("payload shape padding is not zero")
        if any(type(item) is not int or item < 1 for item in dims):
            _fail("payload shape is not a positive integer tuple")
        if (
            type(offset) is not int
            or type(byte_count) is not int
            or offset < DATA_OFFSET
            or byte_count < 8
            or byte_count > MAX_ARRAY_BYTES
            or offset % 8 != 0
            or byte_count % 8 != 0
            or offset + byte_count > len(payload)
        ):
            _fail("payload offset/length is unsafe")
        expected_bytes = 8
        for item in dims:
            expected_bytes *= item
            if expected_bytes > MAX_ARRAY_BYTES:
                raise HLT17MemberCodecResourceError(
                    f"{name} exceeds MAX_ARRAY_BYTES"
                )
        if expected_bytes != byte_count:
            _fail(f"{name} byte count disagrees with its shape")
        digest_hex = digest.hex()
        if _digest_bytes(bytes(view[offset : offset + byte_count])) != digest_hex:
            _fail(f"{name} checksum differs")
        inventory.append(
            {
                "name": name,
                "dtype": "<f8",
                "shape": list(dims),
                "order": "C",
                "offset": offset,
                "byte_count": byte_count,
                "bytes_sha256": digest_hex,
            }
        )
        cursor += ENTRY_BYTES
    expected_offset = DATA_OFFSET
    for item in inventory:
        if item["offset"] != expected_offset:
            _fail("payload arrays are not tightly packed")
        expected_offset += int(item["byte_count"])
    if expected_offset != len(payload):
        _fail("payload has truncated or trailing bytes")
    return inventory


def _allocate_arrays(
    payload: bytes, inventory: list[dict[str, object]]
) -> dict[str, np.ndarray]:
    arrays: dict[str, np.ndarray] = {}
    for item in inventory:
        name = str(item["name"])
        offset = int(item["offset"])
        byte_count = int(item["byte_count"])
        shape = tuple(int(dim) for dim in item["shape"])
        raw = np.frombuffer(payload, dtype="<f8", count=byte_count // 8, offset=offset)
        array = np.array(raw.reshape(shape), dtype=np.dtype("<f8"), order="C", copy=True)
        if not np.isfinite(array).all():
            _fail(f"{name} contains a nonfinite value")
        if _digest_bytes(array.tobytes(order="C")) != item["bytes_sha256"]:
            _fail(f"{name} changed during allocation")
        arrays[name] = array
    if tuple(arrays) != ARRAY_NAMES:
        _fail("decoded array inventory/order differs")
    return arrays


def _expected_shapes(points: int, labels: int, events: int) -> dict[str, tuple[int, ...]]:
    if points < 9 or points > 16385:
        _fail("point count is outside the codec bound")
    if labels < 1 or labels > MAX_LABELS:
        _fail("tracer label count is outside the codec bound")
    if events < 1 or events > MAX_EVENTS:
        _fail("tracer event count is outside the codec bound")
    return {
        "u": (points, 6),
        "p": (points, 6),
        "q": (points, 6),
        "grid_coordinates": (points,),
        "tracer_labels": (labels,),
        "tracer_positions": (labels,),
        "tracer_proper_times": (labels,),
        "event_proper_times": (events, labels),
        "event_fields": (events, labels, 6),
    }


def _validate_inventory(
    inventory: object, metadata: Mapping[str, Any]
) -> list[dict[str, object]]:
    if type(inventory) is not list or len(inventory) != len(ARRAY_NAMES):
        _fail("array inventory/order differs")
    identity = metadata["runtime_identity"]
    tracer = metadata["tracer"]
    expected = _expected_shapes(
        int(identity["point_count"]),
        int(tracer["label_count"]),
        int(tracer["event_count"]),
    )
    checked: list[dict[str, object]] = []
    for name, item in zip(ARRAY_NAMES, inventory, strict=True):
        entry = _exact_mapping(item, INVENTORY_KEYS, f"{name} inventory")
        if (
            entry["name"] != name
            or entry["dtype"] != "<f8"
            or entry["order"] != "C"
        ):
            _fail(f"{name} inventory identity differs")
        shape = entry["shape"]
        if type(shape) is not list or tuple(shape) != expected[name]:
            _fail(f"{name} shape differs")
        checked.append(
            {
                "name": name,
                "dtype": "<f8",
                "shape": [int(dim) for dim in shape],
                "order": "C",
                "offset": _nonnegative_int(entry["offset"], f"{name} offset"),
                "byte_count": _nonnegative_int(entry["byte_count"], f"{name} byte_count"),
                "bytes_sha256": _digest_text(entry["bytes_sha256"], f"{name} bytes_sha256"),
            }
        )
    return checked


def _validate_monitor(value: object) -> dict[str, object]:
    item = _exact_mapping(value, MONITOR_KEYS, "monitor")
    failed = item["first_failed_premise"]
    failed_serial = item["first_failed_transaction_serial"]
    failed_time = item["first_failed_time_hex"]
    if failed is not None and type(failed) is not str:
        _fail("monitor failed premise differs")
    if failed_serial is not None:
        failed_serial = _nonnegative_int(
            failed_serial, "first failed transaction serial"
        )
    if failed_time is not None:
        failed_time = _finite_hex(failed_time, "first failed time")
    if (failed is None) != (failed_serial is None) or (failed is None) != (
        failed_time is None
    ):
        _fail("monitor failure tuple differs")
    return {
        "accepted_stage_count": _nonnegative_int(
            item["accepted_stage_count"], "accepted stage count"
        ),
        "last_transaction_serial": _integer(
            item["last_transaction_serial"], "monitor serial", minimum=-1
        ),
        "last_accepted_time_hex": _finite_hex(
            item["last_accepted_time_hex"], "monitor last accepted time"
        ),
        "first_failed_premise": failed,
        "first_failed_transaction_serial": failed_serial,
        "first_failed_time_hex": failed_time,
    }


def _validate_causal(value: object) -> dict[str, object]:
    item = _exact_mapping(value, CAUSAL_KEYS, "causal")
    return {
        "accepted_time_hex": _finite_hex(item["accepted_time_hex"], "causal accepted time"),
        "accumulated_characteristic_distance_hex": _finite_hex(
            item["accumulated_characteristic_distance_hex"],
            "accumulated characteristic distance",
        ),
        "previous_speed_upper_hex": _finite_hex(
            item["previous_speed_upper_hex"], "previous speed upper"
        ),
    }


def _validate_runtime_identity(value: object) -> dict[str, object]:
    item = _exact_mapping(value, RUNTIME_IDENTITY_KEYS, "runtime identity")
    scope = _text(item["scope"], "scope")
    method = _text(item["method"], "method")
    if method not in INTEGRATOR_BY_LABEL:
        _fail("runtime method label differs")
    point_count = _nonnegative_int(item["point_count"], "point_count")
    member_key = _text(item["member_key"], "member_key")
    amplitude = _text(item["amplitude"], "amplitude")
    integrator = _text(item["integrator_id"], "integrator_id")
    order = _nonnegative_int(item["spatial_order"], "spatial_order")
    if member_key != f"{method}-{point_count}":
        _fail("runtime member key disagrees with method and point count")
    if integrator != INTEGRATOR_BY_LABEL[method] or order != SPATIAL_ORDER_BY_LABEL[method]:
        _fail("runtime integrator or spatial order disagrees with the method")
    if scope == IDENTITY_SCOPE_LIVE:
        if member_key not in HLT17_LIVE_MEMBER_KEYS or amplitude != "3":
            _fail("live HLT17 identity is outside the future production set")
    elif scope == IDENTITY_SCOPE_SYNTHETIC:
        if member_key in HLT17_LIVE_MEMBER_KEYS or amplitude != "synthetic":
            _fail("synthetic identity cannot use a live production key")
    else:
        _fail("runtime identity scope differs")
    try:
        identity = require_hlt17_implementation_identity(
            {name: item[name] for name in IMPLEMENTATION_IDENTITY_KEYS},
            label="runtime identity",
        )
    except (TypeError, ValueError) as error:
        raise HLT17MemberCodecError(str(error)) from error
    return {
        "scope": scope,
        "campaign_id": _text(item["campaign_id"], "campaign_id"),
        "amplitude": amplitude,
        "member_key": member_key,
        "method": method,
        "point_count": point_count,
        "integrator_id": integrator,
        "spatial_order": order,
        "accepted_generation": _nonnegative_int(
            item["accepted_generation"], "accepted_generation"
        ),
        **identity,
    }


def validate_descriptor(metadata: Mapping[str, Any]) -> dict[str, Any]:
    """Validate one complete, non-circular HLT17 descriptor mapping."""

    if not isinstance(metadata, Mapping):
        _fail("codec metadata schema differs")
    if any(key in metadata for key in FORBIDDEN_DESCRIPTOR_KEYS):
        _fail("descriptor embeds a successor cursor")
    item = _exact_mapping(dict(metadata), DESCRIPTOR_KEYS, "codec metadata")
    if item["schema"] == CODEC_LEGACY_IMP1_SCHEMA:
        _fail("legacy IMP1 HLT17 member-state schema cannot be restored as C1R1")
    if item["schema"] != CODEC_SCHEMA or item["schema_version"] != CODEC_SCHEMA_VERSION:
        _fail("codec schema differs")
    if item["protocol_artifact_id"] != RUNTIME_PROTOCOL:
        _fail("runtime protocol differs")
    for name in (
        "campaign_execution_authorized",
        "historical_origin_authenticated",
        "source_operator_closure_authenticated",
        "predecessor_cursor_authenticated",
        "durable_store_bound",
    ):
        _bool_false(item[name], name)
    identity = _validate_runtime_identity(item["runtime_identity"])
    template = item["runtime_template"]
    if type(template) is not dict:
        _fail("runtime template differs")
    template_bytes = _canonical(template, "runtime template")
    if _digest_bytes(template_bytes) != _digest_text(
        item["runtime_template_sha256"], "runtime_template_sha256"
    ):
        _fail("runtime template digest differs")
    if (
        template.get("scope") != identity["scope"]
        or template.get("method") != identity["method"]
        or template.get("point_count") != identity["point_count"]
        or template.get("amplitude") != identity["amplitude"]
    ):
        _fail("template/runtime identity differs")
    closure = _reference_mapping(item["source_operator_closure"], "source_operator_closure")
    origin = _reference_mapping(item["origin_reference"], "origin_receipt")
    accepted_time = _finite_hex(item["accepted_time_hex"], "accepted_time_hex")
    physical = _digest_text(item["physical_state_sha256"], "physical_state_sha256")
    counters = _exact_mapping(item["counters"], COUNTER_KEYS, "counters")
    checked_counters = {
        name: _nonnegative_int(counters[name], name) for name in COUNTER_KEYS
    }
    try:
        ledger = restore_imp1_checkpoint_extension(item["imp1_checkpoint"])
        checkpoint = imp1_checkpoint_extension(ledger)
    except (TypeError, ValueError) as error:
        raise HLT17MemberCodecError("IMP1 checkpoint extension differs") from error
    if _canonical(checkpoint, "IMP1 checkpoint") != _canonical(
        item["imp1_checkpoint"], "IMP1 checkpoint"
    ):
        _fail("IMP1 checkpoint is not a closed encoding")
    if (
        ledger.current_state_sha256 != physical
        or ledger.last_accepted_time.hex() != accepted_time
        or ledger.current_step_index != checked_counters["step_index"]
        or ledger.current_transaction_serial != checked_counters["transaction_serial"]
        or ledger.accepted_macro_step_count != checked_counters["accepted_macro_steps"]
        or ledger.origin_receipt_sha256 != origin["sha256"]
        or ledger.historical_origin_authenticated is not False
    ):
        _fail("IMP1 checkpoint disagrees with the member identity")
    monitor = _validate_monitor(item["monitor"])
    causal = _validate_causal(item["causal"])
    if causal["accepted_time_hex"] != accepted_time:
        _fail("causal/runtime time relation differs")
    if monitor["last_transaction_serial"] + 1 != checked_counters["transaction_serial"]:
        _fail("monitor/runtime serial relation differs")
    if monitor["accepted_stage_count"] != checked_counters["transaction_serial"]:
        _fail("monitor stage count disagrees with the transaction serial")
    if monitor["accepted_stage_count"] == 0:
        if monitor["last_transaction_serial"] != -1:
            _fail("empty monitor serial is not the origin sentinel")
    else:
        if monitor["last_accepted_time_hex"] != accepted_time:
            _fail("nonempty monitor accepted time differs from member time")
    tracer = _exact_mapping(item["tracer"], TRACER_KEYS, "tracer")
    checked_tracer = {
        "label_count": _nonnegative_int(tracer["label_count"], "tracer label count"),
        "event_count": _nonnegative_int(tracer["event_count"], "tracer event count"),
        "cutoff_hex": _finite_hex(tracer["cutoff_hex"], "tracer cutoff", positive=True),
        "outer_radius_hex": _finite_hex(
            tracer["outer_radius_hex"], "tracer outer radius", positive=True
        ),
    }
    snapshot_role = _text(item["snapshot_role"], "snapshot_role")
    if snapshot_role not in SNAPSHOT_ROLES:
        _fail("snapshot role differs")
    accepted_generation = identity["accepted_generation"]
    plan_value = item["executed_plan"]
    predecessor_cursor = _validated_predecessor_cursor(item["predecessor_cursor_identity"])
    predecessor_identity_sha256 = _digest_text(
        item["predecessor_cursor_identity_sha256"],
        "predecessor_cursor_identity_sha256",
    )
    if predecessor_identity_sha256 != _digest_bytes(
        _canonical(predecessor_cursor, "predecessor cursor identity")
    ):
        _fail("predecessor cursor identity digest differs")
    predecessor_cursor_sha256 = _digest_text(
        item["predecessor_cursor_sha256"], "predecessor_cursor_sha256"
    )
    if predecessor_identity_sha256 == predecessor_cursor_sha256:
        _fail(
            "projected cursor-identity digest and supplied full-cursor digest "
            "must remain distinct"
        )
    predecessor_descriptor = _digest_text(
        item["predecessor_descriptor_sha256"], "predecessor_descriptor_sha256"
    )
    if predecessor_cursor["member_key"] != identity["member_key"]:
        _fail("predecessor/runtime member relation differs")
    if predecessor_cursor["method"] != identity["method"]:
        _fail("predecessor/runtime method relation differs")
    if predecessor_cursor["point_count"] != identity["point_count"]:
        _fail("predecessor/runtime point-count relation differs")
    if predecessor_cursor["campaign_id"] != identity["campaign_id"]:
        _fail("predecessor/runtime campaign relation differs")
    stage_records = IMP1_METHOD_STAGE_RECORD_COUNT[identity["integrator_id"]]
    if snapshot_role == SNAPSHOT_ROLE_ACCEPTED:
        if ledger.current_macro_step_temporal_retry_count != 0:
            _fail("accepted descriptor cannot carry an active retry ledger")
        if accepted_generation == 0:
            if plan_value is not None:
                _fail("origin accepted encoding cannot carry an executed plan")
            plan = None
            if predecessor_cursor["protocol_artifact_id"] != ORIGIN_PROTOCOL:
                _fail("origin predecessor protocol differs")
            if predecessor_cursor["accepted_generation"] != 0:
                _fail("origin predecessor accepted generation differs")
            if predecessor_cursor["accepted_boundary_time_hex"] != accepted_time:
                _fail("origin predecessor time differs")
            if predecessor_cursor["accepted_state_sha256"] != physical:
                _fail("origin predecessor physical hash differs")
            if predecessor_cursor["previous_step_index"] != checked_counters["step_index"]:
                _fail("origin predecessor step differs")
            if (
                predecessor_cursor["previous_transaction_serial"]
                != checked_counters["transaction_serial"]
            ):
                _fail("origin predecessor serial differs")
            if checked_counters["accepted_macro_steps"] != 0:
                _fail("origin accepted macro-step count differs")
        else:
            plan = decode_plan(plan_value)
            if predecessor_cursor["protocol_artifact_id"] != RUNTIME_PROTOCOL:
                _fail("accepted-successor predecessor protocol differs")
            if predecessor_cursor["accepted_generation"] != accepted_generation - 1:
                _fail("predecessor accepted generation is not the immediate ancestor")
            if plan.current.hex() != predecessor_cursor["accepted_boundary_time_hex"]:
                _fail("executed plan current is not the predecessor accepted time")
            if plan.boundaries[4].hex() != accepted_time:
                _fail("executed plan endpoint is not the accepted time")
            predecessor_time = _from_hex(
                predecessor_cursor["accepted_boundary_time_hex"],
                "predecessor accepted time",
            )
            current_time = _from_hex(accepted_time, "accepted_time_hex")
            if not (predecessor_time < current_time):
                _fail("accepted successor time is not strictly later")
            if (
                checked_counters["step_index"]
                != predecessor_cursor["previous_step_index"] + IMP1_ACCEPTED_SUBSTEP_COUNT
            ):
                _fail("accepted successor does not advance exactly four substeps")
            if (
                checked_counters["transaction_serial"]
                != predecessor_cursor["previous_transaction_serial"] + stage_records
            ):
                _fail("accepted successor serial is not the method-owned stage increment")
            if checked_counters["accepted_macro_steps"] != accepted_generation:
                _fail("accepted generation disagrees with accepted macro-step count")
    else:
        if ledger.current_macro_step_temporal_retry_count < 1:
            _fail("retry overlay omits an active IMP1 retry ledger")
        if plan_value is not None:
            _fail("retry overlay cannot carry an executed accepted plan")
        plan = None
        if predecessor_cursor["accepted_generation"] != accepted_generation:
            _fail("retry overlay accepted generation differs from its accepted boundary")
        if predecessor_cursor["accepted_boundary_time_hex"] != accepted_time:
            _fail("retry overlay cannot advance accepted time")
        if predecessor_cursor["accepted_state_sha256"] != physical:
            _fail("retry overlay cannot replace the accepted physical state")
        if predecessor_cursor["previous_step_index"] != checked_counters["step_index"]:
            _fail("retry overlay cannot advance the step index")
        if (
            predecessor_cursor["previous_transaction_serial"]
            != checked_counters["transaction_serial"]
        ):
            _fail("retry overlay cannot advance the transaction serial")
        if accepted_generation == 0:
            if predecessor_cursor["protocol_artifact_id"] not in {
                ORIGIN_PROTOCOL,
                RUNTIME_PROTOCOL,
            }:
                _fail("origin retry predecessor protocol differs")
            if checked_counters["accepted_macro_steps"] != 0:
                _fail("origin retry accepted macro-step count differs")
        else:
            if predecessor_cursor["protocol_artifact_id"] != RUNTIME_PROTOCOL:
                _fail("retry overlay predecessor protocol differs")
            if checked_counters["accepted_macro_steps"] != accepted_generation:
                _fail("retry overlay accepted generation disagrees with macro-step count")
    inventory = _validate_inventory(item["array_inventory"], {"runtime_identity": identity, "tracer": checked_tracer})
    payload_sha256 = _digest_text(item["payload_sha256"], "payload_sha256")
    payload_byte_count = _nonnegative_int(item["payload_byte_count"], "payload_byte_count")
    if payload_byte_count > MAX_PAYLOAD_BYTES:
        raise HLT17MemberCodecResourceError("payload exceeds MAX_PAYLOAD_BYTES")
    if physical == payload_sha256 or physical == predecessor_descriptor:
        _fail("physical state hash is not distinct from codec/provenance digests")
    result = {
        "schema": CODEC_SCHEMA,
        "schema_version": CODEC_SCHEMA_VERSION,
        "protocol_artifact_id": RUNTIME_PROTOCOL,
        "snapshot_role": snapshot_role,
        "runtime_identity": identity,
        "runtime_template": dict(template),
        "runtime_template_sha256": str(item["runtime_template_sha256"]),
        "source_operator_closure": closure,
        "origin_reference": origin,
        "accepted_time_hex": accepted_time,
        "physical_state_sha256": physical,
        "counters": checked_counters,
        "imp1_checkpoint": dict(item["imp1_checkpoint"]),
        "monitor": monitor,
        "causal": causal,
        "tracer": checked_tracer,
        "executed_plan": None if plan is None else encode_plan(plan),
        "predecessor_descriptor_sha256": predecessor_descriptor,
        "predecessor_cursor_identity": predecessor_cursor,
        "predecessor_cursor_identity_sha256": predecessor_identity_sha256,
        "predecessor_cursor_sha256": predecessor_cursor_sha256,
        "array_inventory": inventory,
        "payload_sha256": payload_sha256,
        "payload_byte_count": payload_byte_count,
        "campaign_execution_authorized": False,
        "historical_origin_authenticated": False,
        "source_operator_closure_authenticated": False,
        "predecessor_cursor_authenticated": False,
        "durable_store_bound": False,
    }
    _canonical(result, "codec metadata")
    return result


def _member_arrays(member: HLT17RuntimeMember) -> dict[str, np.ndarray]:
    tracers = member.tracers
    return {
        "u": _strict_array(member.state.u, "u"),
        "p": _strict_array(member.state.p, "p"),
        "q": _strict_array(member.state.q, "q"),
        "grid_coordinates": _strict_array(member.coordinates, "grid_coordinates"),
        "tracer_labels": _strict_array(tracers.labels, "tracer_labels"),
        "tracer_positions": _strict_array(tracers.positions, "tracer_positions"),
        "tracer_proper_times": _strict_array(tracers.proper_times, "tracer_proper_times"),
        "event_proper_times": _strict_array(
            np.stack(tuple(tracers.event_proper_times)), "event_proper_times"
        ),
        "event_fields": _strict_array(
            np.stack(tuple(tracers.event_fields)), "event_fields"
        ),
    }


def _encode_snapshot(
    member: HLT17RuntimeMember,
    predecessor: HLT17PredecessorIdentity,
    *,
    accepted_generation: int,
    snapshot_role: str,
) -> HLT17MemberEncoding:
    if not isinstance(member, HLT17RuntimeMember):
        raise TypeError("member must be HLT17RuntimeMember")
    if not isinstance(predecessor, HLT17PredecessorIdentity):
        raise TypeError("predecessor must be HLT17PredecessorIdentity")
    member.__post_init__()
    accepted_generation = _nonnegative_int(accepted_generation, "accepted_generation")
    if snapshot_role not in SNAPSHOT_ROLES:
        _fail("snapshot role differs")
    if member.identity.scope == IDENTITY_SCOPE_LIVE:
        if member.key not in HLT17_LIVE_MEMBER_KEYS:
            _fail("live HLT17 identity is outside the future production set")
    elif member.key in HLT17_LIVE_MEMBER_KEYS:
        _fail("synthetic identity cannot use a live production key")
    arrays = _member_arrays(member)
    physical = array_content_sha256(arrays["u"], arrays["p"], arrays["q"])
    if physical != member.physical_state_sha256():
        _fail("physical u/p/q hash differs from the live member")
    if np.any(np.diff(arrays["tracer_labels"]) <= 0.0):
        _fail("tracer labels are not strictly increasing")
    if array_content_sha256(arrays["grid_coordinates"]) != member.coordinates_sha256:
        _fail("encoded coordinates differ from the live grid identity")
    inventory, payload = _inventory_and_payload(arrays)
    template = member.runtime_template()
    checkpoint = imp1_checkpoint_extension(member.temporal_ledger)
    plan: object
    if snapshot_role == SNAPSHOT_ROLE_ACCEPTED:
        if member.temporal_ledger.current_macro_step_temporal_retry_count != 0:
            _fail("accepted encoding cannot carry an active retry ledger")
        if accepted_generation == 0:
            if member.executed_plan is not None:
                _fail("origin accepted member cannot carry an executed plan")
            plan = None
        else:
            if member.executed_plan is None:
                _fail("accepted successor requires the exact executed TDG7 plan")
            plan = encode_plan(member.executed_plan)
    else:
        if member.temporal_ledger.current_macro_step_temporal_retry_count < 1:
            _fail("retry overlay requires an active IMP1 retry ledger")
        plan = None
    monitor = member.transaction.state
    causal = member.transaction.causal_state
    metadata = {
        "schema": CODEC_SCHEMA,
        "schema_version": CODEC_SCHEMA_VERSION,
        "protocol_artifact_id": RUNTIME_PROTOCOL,
        "snapshot_role": snapshot_role,
        "runtime_identity": {
            "scope": member.identity.scope,
            "campaign_id": member.identity.campaign_id,
            "amplitude": member.identity.amplitude,
            "member_key": member.key,
            "method": member.identity.method_label,
            "point_count": member.identity.point_count,
            "integrator_id": member.identity.integrator_id,
            "spatial_order": member.identity.spatial_order,
            "accepted_generation": accepted_generation,
            **hlt17_implementation_identity(),
        },
        "runtime_template": template,
        "runtime_template_sha256": member.runtime_template_sha256(),
        "source_operator_closure": member.source_operator_closure.as_mapping(),
        "origin_reference": member.origin_reference.as_mapping(),
        "accepted_time_hex": member.time.hex(),
        "physical_state_sha256": physical,
        "counters": {
            "step_index": member.step_index,
            "transaction_serial": member.transaction_serial,
            "accepted_macro_steps": member.temporal_ledger.accepted_macro_step_count,
            "source_retry_count": member.source_retry_count,
            "CFL_retry_count": member.CFL_retry_count,
        },
        "imp1_checkpoint": checkpoint,
        "monitor": {
            "accepted_stage_count": monitor.accepted_stage_count,
            "last_transaction_serial": monitor.last_transaction_serial,
            "last_accepted_time_hex": float(monitor.last_accepted_time).hex(),
            "first_failed_premise": monitor.first_failed_premise,
            "first_failed_transaction_serial": monitor.first_failed_transaction_serial,
            "first_failed_time_hex": (
                None
                if monitor.first_failed_time is None
                else float(monitor.first_failed_time).hex()
            ),
        },
        "causal": {
            "accepted_time_hex": float(causal.accepted_time).hex(),
            "accumulated_characteristic_distance_hex": float(
                causal.accumulated_characteristic_distance
            ).hex(),
            "previous_speed_upper_hex": float(causal.previous_speed_upper).hex(),
        },
        "tracer": {
            "label_count": int(arrays["tracer_labels"].size),
            "event_count": int(arrays["event_fields"].shape[0]),
            "cutoff_hex": float(member.tracers.cutoff).hex(),
            "outer_radius_hex": float(member.tracers.outer_radius).hex(),
        },
        "executed_plan": plan,
        "predecessor_descriptor_sha256": predecessor.descriptor_sha256,
        "predecessor_cursor_identity": dict(predecessor.cursor_identity),
        "predecessor_cursor_identity_sha256": predecessor.identity_sha256,
        "predecessor_cursor_sha256": predecessor.cursor_sha256,
        "array_inventory": inventory,
        "payload_sha256": _digest_bytes(payload),
        "payload_byte_count": len(payload),
        "campaign_execution_authorized": False,
        "historical_origin_authenticated": False,
        "source_operator_closure_authenticated": False,
        "predecessor_cursor_authenticated": False,
        "durable_store_bound": False,
    }
    checked = validate_descriptor(metadata)
    descriptor = _canonical(checked, "codec metadata")
    descriptor_sha256 = _digest_bytes(descriptor)
    payload_sha256 = _digest_bytes(payload)
    if descriptor_sha256 == physical or payload_sha256 == physical:
        _fail("descriptor or payload hash collided with the physical u/p/q hash")
    return HLT17MemberEncoding(descriptor=descriptor, payload=payload)


def encode_member(
    member: HLT17RuntimeMember,
    predecessor: HLT17PredecessorIdentity,
    *,
    generation: int,
) -> HLT17MemberEncoding:
    """Encode one accepted physical-state boundary.

    ``generation`` is accepted physical-state generation, not a cursor or
    journal revision.  Active retry ledgers must use encode_retry_overlay.
    """

    return _encode_snapshot(
        member,
        predecessor,
        accepted_generation=generation,
        snapshot_role=SNAPSHOT_ROLE_ACCEPTED,
    )


def encode_retry_overlay(
    member: HLT17RuntimeMember,
    predecessor: HLT17PredecessorIdentity,
    *,
    accepted_generation: int,
) -> HLT17MemberEncoding:
    """Encode an IMP1 retry ledger on an unchanged accepted physical boundary."""

    return _encode_snapshot(
        member,
        predecessor,
        accepted_generation=accepted_generation,
        snapshot_role=SNAPSHOT_ROLE_RETRY,
    )


def decode_member(descriptor: bytes, payload: bytes) -> HLT17DecodedSnapshot:
    """Decode only after schema, checksum, and bound checks succeed."""

    metadata = validate_descriptor(_load_descriptor(descriptor))
    descriptor_sha256 = _digest_bytes(bytes(descriptor))
    if type(payload) is not bytes:
        raise TypeError("payload must be bytes")
    if len(payload) != metadata["payload_byte_count"]:
        _fail("payload length differs from the descriptor")
    if len(payload) > MAX_PAYLOAD_BYTES:
        raise HLT17MemberCodecResourceError("payload exceeds MAX_PAYLOAD_BYTES")
    payload_sha256 = _digest_bytes(payload)
    if payload_sha256 != metadata["payload_sha256"]:
        _fail("payload checksum differs")
    table = _parse_payload_table(payload)
    if table != metadata["array_inventory"]:
        _fail("payload inventory disagrees with the descriptor")
    allocated = _allocate_arrays(payload, table)
    arrays = MappingProxyType(
        {
            name: np.array(allocated[name], dtype=np.dtype("<f8"), order="C", copy=True)
            for name in ARRAY_NAMES
        }
    )
    for array in arrays.values():
        array.setflags(write=False)
    identity = metadata["runtime_identity"]
    tracer = metadata["tracer"]
    expected = _expected_shapes(
        int(identity["point_count"]),
        int(tracer["label_count"]),
        int(tracer["event_count"]),
    )
    if any(arrays[name].shape != shape for name, shape in expected.items()):
        _fail("decoded array shapes differ")
    template_grid = metadata["runtime_template"].get("grid")
    if (
        type(template_grid) is not dict
        or array_content_sha256(arrays["grid_coordinates"])
        != template_grid.get("coordinates_sha256")
        or np.any(np.diff(arrays["tracer_labels"]) <= 0.0)
    ):
        _fail("grid or tracer labels differ from the template")
    physical = array_content_sha256(arrays["u"], arrays["p"], arrays["q"])
    if physical != metadata["physical_state_sha256"]:
        _fail("physical u/p/q hash differs from the descriptor")
    if physical == descriptor_sha256 or physical == payload_sha256:
        _fail("physical state hash is not distinct")
    ledger = restore_imp1_checkpoint_extension(metadata["imp1_checkpoint"])
    if ledger.current_state_sha256 != physical:
        _fail("IMP1 ledger physical hash differs from decoded u/p/q")
    plan = (
        None
        if metadata["executed_plan"] is None
        else decode_plan(metadata["executed_plan"])
    )
    return HLT17DecodedSnapshot(
        arrays=arrays,
        metadata=MappingProxyType(metadata),
        ledger=ledger,
        executed_plan=plan,
        descriptor=bytes(descriptor),
        payload=bytes(payload),
        descriptor_sha256=descriptor_sha256,
        physical_state_sha256=physical,
        payload_sha256=payload_sha256,
    )


def _capture_from_decoded(
    member: HLT17RuntimeMember, decoded: HLT17DecodedSnapshot
) -> HLT17MemberCapture:
    metadata = decoded.metadata
    identity = metadata["runtime_identity"]
    if (
        identity["member_key"] != member.key
        or identity["scope"] != member.identity.scope
        or identity["campaign_id"] != member.identity.campaign_id
        or identity["method"] != member.identity.method_label
        or identity["point_count"] != member.identity.point_count
        or identity["integrator_id"] != member.identity.integrator_id
    ):
        _fail("restore branch/member differs")
    if metadata["runtime_template"] != member.runtime_template():
        _fail("restore shell template differs")
    if metadata["runtime_template_sha256"] != member.runtime_template_sha256():
        _fail("restore template digest differs")
    if metadata["source_operator_closure"] != member.source_operator_closure.as_mapping():
        _fail("restore source/operator closure reference differs")
    if metadata["origin_reference"] != member.origin_reference.as_mapping():
        _fail("restore origin reference differs")
    if array_content_sha256(decoded.arrays["grid_coordinates"]) != member.coordinates_sha256:
        _fail("restore coordinates differ from the live grid identity")
    monitor = metadata["monitor"]
    causal = metadata["causal"]
    tracer = metadata["tracer"]
    arrays = decoded.arrays
    failed_time = monitor["first_failed_time_hex"]
    return HLT17MemberCapture(
        state=EvolutionState(arrays["u"], arrays["p"], arrays["q"]),
        time=_from_hex(metadata["accepted_time_hex"], "accepted_time_hex"),
        step_index=int(metadata["counters"]["step_index"]),
        transaction_serial=int(metadata["counters"]["transaction_serial"]),
        source_retry_count=int(metadata["counters"]["source_retry_count"]),
        CFL_retry_count=int(metadata["counters"]["CFL_retry_count"]),
        monitor_state=GR0RuntimeMonitorState(
            accepted_stage_count=int(monitor["accepted_stage_count"]),
            last_transaction_serial=int(monitor["last_transaction_serial"]),
            last_accepted_time=_from_hex(
                monitor["last_accepted_time_hex"], "monitor last accepted time"
            ),
            first_failed_premise=monitor["first_failed_premise"],
            first_failed_transaction_serial=monitor["first_failed_transaction_serial"],
            first_failed_time=(
                None if failed_time is None else _from_hex(failed_time, "first failed time")
            ),
        ),
        causal_state=CausalBudgetState(
            accepted_time=_from_hex(causal["accepted_time_hex"], "causal accepted time"),
            accumulated_characteristic_distance=_from_hex(
                causal["accumulated_characteristic_distance_hex"],
                "accumulated characteristic distance",
            ),
            previous_speed_upper=_from_hex(
                causal["previous_speed_upper_hex"], "previous speed upper"
            ),
        ),
        tracer_labels=arrays["tracer_labels"],
        tracer_positions=arrays["tracer_positions"],
        tracer_proper_times=arrays["tracer_proper_times"],
        event_proper_times=tuple(row.copy() for row in arrays["event_proper_times"]),
        event_fields=tuple(row.copy() for row in arrays["event_fields"]),
        tracer_cutoff=_from_hex(tracer["cutoff_hex"], "tracer cutoff"),
        tracer_outer_radius=_from_hex(tracer["outer_radius_hex"], "tracer outer radius"),
        temporal_ledger=decoded.ledger,
        executed_plan=decoded.executed_plan,
        coordinates_sha256=member.coordinates_sha256,
        inherited_snapshot=member.temporal_ledger.inherited_snapshot,
        operator=member.operator,
        projector=member.projector,
        initial=member.initial,
    )


def _check_encoding_hashes(encoding: HLT17MemberEncoding) -> None:
    if not isinstance(encoding, HLT17MemberEncoding):
        raise TypeError("encoding must be HLT17MemberEncoding")
    encoding.__post_init__()


def _encoding_from(
    encoded: object,
) -> tuple[HLT17MemberEncoding, HLT17DecodedSnapshot | None]:
    if isinstance(encoded, HLT17MemberEncoding):
        return encoded, None
    if isinstance(encoded, HLT17DecodedSnapshot):
        if type(encoded.descriptor) is not bytes or type(encoded.payload) is not bytes:
            _fail("restore requires original encoding bytes")
        encoding = HLT17MemberEncoding(
            descriptor=encoded.descriptor,
            payload=encoded.payload,
        )
        return encoding, encoded
    _fail("restore requires original encoding bytes")
    raise AssertionError("unreachable")


def _require_snapshot_matches(
    claimed: HLT17DecodedSnapshot, authentic: HLT17DecodedSnapshot
) -> None:
    if claimed.descriptor != authentic.descriptor or claimed.payload != authentic.payload:
        _fail("decoded snapshot original bytes differ")
    if (
        claimed.descriptor_sha256 != authentic.descriptor_sha256
        or claimed.payload_sha256 != authentic.payload_sha256
        or claimed.physical_state_sha256 != authentic.physical_state_sha256
    ):
        _fail("decoded snapshot hashes differ from original bytes")
    if dict(claimed.metadata) != dict(authentic.metadata):
        _fail("decoded snapshot metadata was altered")
    if not isinstance(claimed.arrays, Mapping) or tuple(claimed.arrays) != ARRAY_NAMES:
        _fail("decoded snapshot array inventory/order differs")
    for name in ARRAY_NAMES:
        claimed_array = claimed.arrays[name]
        if (
            not isinstance(claimed_array, np.ndarray)
            or claimed_array.tobytes() != authentic.arrays[name].tobytes()
        ):
            _fail(f"decoded snapshot {name} was altered")
    try:
        claimed_ledger = imp1_checkpoint_extension(claimed.ledger)
        authentic_ledger = imp1_checkpoint_extension(authentic.ledger)
    except (TypeError, ValueError) as error:
        raise HLT17MemberCodecError("decoded snapshot ledger was altered") from error
    if _canonical(claimed_ledger, "claimed ledger") != _canonical(
        authentic_ledger, "authentic ledger"
    ):
        _fail("decoded snapshot ledger was altered")
    try:
        claimed_plan = (
            None if claimed.executed_plan is None else encode_plan(claimed.executed_plan)
        )
        authentic_plan = (
            None
            if authentic.executed_plan is None
            else encode_plan(authentic.executed_plan)
        )
    except (TypeError, ValueError) as error:
        raise HLT17MemberCodecError("decoded snapshot plan was altered") from error
    if claimed_plan != authentic_plan:
        _fail("decoded snapshot plan was altered")


def restore_member(
    member: HLT17RuntimeMember,
    encoded: HLT17MemberEncoding | HLT17DecodedSnapshot,
    *,
    fault_hook: Callable[[str], None] | None = None,
) -> HLT17DecodedSnapshot:
    """Restore only after original bytes, inventories, and hashes reauthenticate."""

    if not isinstance(member, HLT17RuntimeMember):
        raise TypeError("member must be HLT17RuntimeMember")
    encoding, claimed = _encoding_from(encoded)
    _check_encoding_hashes(encoding)
    decoded = decode_member(encoding.descriptor, encoding.payload)
    if (
        decoded.descriptor_sha256 != encoding.descriptor_sha256
        or decoded.physical_state_sha256 != encoding.physical_state_sha256
        or decoded.payload_sha256 != encoding.payload_sha256
    ):
        _fail("encoding hashes disagree with the decoded object")
    if claimed is not None:
        _require_snapshot_matches(claimed, decoded)
    capture = _capture_from_decoded(member, decoded)
    try:
        member.restore(capture, fault_hook=fault_hook)
    except HLT17RuntimeMemberError as error:
        raise HLT17MemberCodecError(str(error)) from error
    return decoded


__all__ = [
    "ARRAY_NAMES",
    "CODEC_LEGACY_IMP1_SCHEMA",
    "CODEC_SCHEMA",
    "CODEC_SCHEMA_VERSION",
    "HLT17_CURSOR_BRIDGE_SEAMS",
    "HLT17DecodedSnapshot",
    "HLT17MemberCodecError",
    "HLT17MemberCodecResourceError",
    "HLT17MemberEncoding",
    "HLT17PredecessorIdentity",
    "MAX_ARRAY_BYTES",
    "MAX_DESCRIPTOR_BYTES",
    "MAX_PAYLOAD_BYTES",
    "ORIGIN_PROTOCOL",
    "PREDECESSOR_CURSOR_IDENTITY_KEYS",
    "RUNTIME_PROTOCOL",
    "SNAPSHOT_ROLE_ACCEPTED",
    "SNAPSHOT_ROLE_RETRY",
    "decode_member",
    "decode_plan",
    "encode_member",
    "encode_plan",
    "encode_retry_overlay",
    "restore_member",
    "validate_descriptor",
]
