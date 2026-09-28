"""Immutable HLT17/PROTO19 cursor, codec, and pure retry transitions.

This module is the first synthetic HLT17 slice: a new schema for durable
fresh/retry identity, complete TDG7 plans, independent source/CFL/temporal
counters, and a composite overlay that does not overwrite temporal identity.
The in-memory path binds the fixed C1R1 admission implementation. IMP1 remains
the sealed ledger/rejection wire. This module does not read a store,
authenticate an origin or physical source, prepare a production grid, or
authorize a campaign.
"""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
import math
from numbers import Real
import re
from typing import Final, Mapping

from recursive_horizons.evidence_io import CanonicalJSONError, canonical_json_bytes

from .hlt17_admission_runtime import (
    HLT17_ADMISSION_RUNTIME_ID,
    HLT17_IMPLEMENTATION_ID,
    HLT17_MATHEMATICAL_OBJECT,
    HLT17_REFERENCE_WIRE_ARTIFACT_ID,
    HLT17_REFERENCE_WIRE_EVALUATOR_ID,
    IMPLEMENTATION_IDENTITY_KEYS,
    hlt17_implementation_identity,
    require_hlt17_implementation_identity,
)
from .numerical_engine import COMPARATOR_METHOD, PRIMARY_METHOD
from .tdg11_imp1_enclosure import IMP1EnclosureLimits
from .tdg11_imp1_ledger import (
    IMP1_METHOD_LABEL,
    TDG11IMP1Error,
    TDG11IMP1Ledger,
    TDG11IMP1TemporalRejection,
    append_imp1_rejection,
    imp1_checkpoint_extension,
    restore_imp1_checkpoint_extension,
)
from .tdg6_temporal_admission_design import TDG6_MINIMUM_MACRO_STEP
from .tdg7_binary64_subdivision_lattice import (
    CoordinateLatticeLimitReached,
    TDG7Binary64SubdivisionPlan,
    plan_forward_proto14_subdivision,
    validate_tdg7_binary64_subdivision_plan,
)


HLT17_CURSOR_SCHEMA: Final[str] = "FGC-1-HLT17-C1R1-cursor-v1"
HLT17_CURSOR_SCHEMA_VERSION: Final[int] = 1
HLT17_LEGACY_IMP1_CURSOR_SCHEMA: Final[str] = "FGC-1-HLT17-IMP1-cursor-v1"
HLT17_SOURCE_RETRY_CAP: Final[int] = 32
HLT17_CFL_RETRY_CAP: Final[int] = 32
HLT17_GENESIS_DIGEST: Final[str] = "0" * 64
MAX_CURSOR_BYTES: Final[int] = 16 * 1024 * 1024
FRESH_READY: Final[str] = "FRESH_READY"
RETRY_PENDING: Final[str] = "RETRY_PENDING"
SOURCE_OWNER: Final[str] = "source"
CFL_OWNER: Final[str] = "cfl"
TEMPORAL_OWNER: Final[str] = "temporal"
IDENTITY_SCOPE_SYNTHETIC: Final[str] = "synthetic"
IDENTITY_SCOPE_LIVE: Final[str] = "live"
IDENTITY_SCOPES: Final[tuple[str, ...]] = (
    IDENTITY_SCOPE_SYNTHETIC,
    IDENTITY_SCOPE_LIVE,
)

FUTURE_MEMBER_CODEC_SEAM: Final[str] = (
    "protocol_v19 member codec must bind complete physical u/p/q arrays, the "
    "serialized descriptor, operator/source closure bytes, monitor/causal/"
    "tracer/guard identities, this cursor, and both IMP1 ledgers"
)
FUTURE_SOURCE_MANIFEST_SEAM: Final[str] = (
    "physical source/operator full-code and configuration authentication is an "
    "external manifest; the synthetic callable-binding digest is not that proof"
)
FUTURE_CONTENT_ADDRESSED_STORE_SEAM: Final[str] = (
    "a new content-addressed store must persist cursor+ledger+physical state "
    "without treating this in-memory codec as a writer or origin proof"
)
FUTURE_CHECKPOINT_JOURNAL_SEAM: Final[str] = (
    "checkpoint/journal schema must validate a proposed transition before "
    "publishing its journal; kernel_transition_parent_sha256 and "
    "kernel_transition_digest are hashes of this in-memory cursor "
    "transition, not a store-journal parent/tip or durable publication. "
    "Actual FUTURE store journal and store checkpoint revisions stay unbound "
    "until the external journal owner exists"
)
FUTURE_AUTHENTICATED_ORIGIN_SEAM: Final[str] = (
    "HLT17 production seeding must independently verify the synchronized "
    "HLT15 t=23/16 origin; this cursor never sets historical origin true"
)
FUTURE_CAMPAIGN_TRANSACTION_SEAM: Final[str] = (
    "campaign recovery must reject malformed tails, forks, partial terminals "
    "and informal takeover; mixed source/CFL/temporal state is represented "
    "here and must not be collapsed to an exclusive pending_owner"
)

PLAN_KEYS: Final[tuple[str, ...]] = (
    "current_hex",
    "event_target_hex",
    "requested_cap_hex",
    "quantum_hex",
    "selection_budget_hex",
    "macro_width_hex",
    "fine_width_hex",
    "conservative_reduction_hex",
    "minimum_width_hex_or_none",
    "boundaries_hex",
    "target_limited",
)
TEMPORAL_KEYS: Final[tuple[str, ...]] = (
    "rejection",
    "predecessor_plan",
    "successor_plan",
    "predecessor_imp1_ledger",
    "predecessor_imp1_ledger_sha256",
    "preparation_sha256",
    "assessment_sha256",
    "rejection_sha256",
    "imp1_half_cap_hex",
)
OVERLAY_KEYS: Final[tuple[str, ...]] = (
    "owner",
    "attempted_plan",
    "next_plan",
    "evidence",
    "owner_current_count",
    "exhausted",
)
OVERLAY_EVIDENCE_COMMON: Final[tuple[str, ...]] = (
    "kind",
    "path",
    "attempted_requested_cap_hex",
    "attempted_macro_width_hex",
    "accepted_boundary_preserved",
    "physical_classification",
)
SOURCE_EVIDENCE_KEYS: Final[tuple[str, ...]] = OVERLAY_EVIDENCE_COMMON + (
    "source_retry",
)
CFL_EVIDENCE_KEYS: Final[tuple[str, ...]] = OVERLAY_EVIDENCE_COMMON + (
    "observed_ratio_hex",
    "maximum_ratio_hex",
)
CURSOR_KEYS: Final[tuple[str, ...]] = (
    "schema",
    "schema_version",
    "implementation_id",
    "mathematical_object",
    "reference_wire_artifact_id",
    "reference_wire_evaluator_id",
    "admission_runtime_id",
    "member_key",
    "identity_scope",
    "method",
    "method_label",
    "point_count",
    "owned_rows",
    "fallback_d01",
    "fallback_d12",
    "mode",
    "event_target_hex",
    "accepted_time_hex",
    "accepted_step_index",
    "accepted_transaction_serial",
    "descriptor_sha256",
    "physical_state_sha256",
    "coordinates_sha256",
    "synthetic_callable_binding_sha256",
    "monitor_sha256",
    "causal_state_sha256",
    "tracer_sha256",
    "guard_configuration_sha256",
    "imp1_ledger",
    "imp1_ledger_sha256",
    "inherited_snapshot_sha256",
    "temporal",
    "overlay",
    "source_current",
    "source_total",
    "cfl_current",
    "cfl_total",
    "source_cap",
    "cfl_cap",
    "cursor_generation",
    "attempt_serial",
    "checkpoint_generation",
    "kernel_transition_parent_sha256",
    "kernel_transition_digest",
    "historical_origin_authenticated",
    "physical_classification",
    "campaign_execution_authorized",
    "origin_authenticated_by_this_cursor",
    "source_authenticated_by_this_cursor",
    "store_authenticated_by_this_cursor",
)

_CURSOR_KEY_SET = frozenset(CURSOR_KEYS)
_PLAN_KEY_SET = frozenset(PLAN_KEYS)
_TEMPORAL_KEY_SET = frozenset(TEMPORAL_KEYS)
_OVERLAY_KEY_SET = frozenset(OVERLAY_KEYS)
_SHA256 = re.compile(r"\A[0-9a-f]{64}\Z")
_METHOD_BY_LABEL: Final[dict[str, str]] = {
    "RK4": PRIMARY_METHOD,
    "SSPRK3": COMPARATOR_METHOD,
}


class HLT17IMP1Error(ValueError):
    """Malformed HLT17 cursor, plan, overlay, or transition data."""


class HLT17InvalidPremise(HLT17IMP1Error):
    """Replay or authentication failed; not a new source/CFL/temporal attempt."""

    spends_retry = False
    spends_debit = False
    disposition = "invalid_premise"


class HLT17OwnerRetryExhausted(RuntimeError):
    """Typed source/CFL/temporal lattice or budget exhaustion."""

    physical_classification = False
    accepted_state_advanced = False

    def __init__(
        self,
        owner: str,
        reason: str,
        *,
        attempted_plan: TDG7Binary64SubdivisionPlan | None = None,
        ledger: TDG11IMP1Ledger | None = None,
    ) -> None:
        if owner not in {SOURCE_OWNER, CFL_OWNER, TEMPORAL_OWNER}:
            raise HLT17IMP1Error("unknown exhausted retry owner")
        self.owner = owner
        self.reason = reason
        self.attempted_plan = attempted_plan
        self.ledger = ledger
        super().__init__(f"HLT17 {owner} retry exhausted: {reason}")


def _fail(message: str) -> None:
    raise HLT17IMP1Error(message)


def _implementation_fields(
    cursor: "HLT17IMP1Cursor | None" = None,
) -> dict[str, str]:
    if cursor is None:
        return dict(hlt17_implementation_identity())
    return {
        "implementation_id": cursor.implementation_id,
        "mathematical_object": cursor.mathematical_object,
        "reference_wire_artifact_id": cursor.reference_wire_artifact_id,
        "reference_wire_evaluator_id": cursor.reference_wire_evaluator_id,
        "admission_runtime_id": cursor.admission_runtime_id,
    }


def _digest(payload: bytes) -> str:
    return sha256(payload).hexdigest()


def _text(value: object, label: str) -> str:
    if type(value) is not str:
        raise TypeError(f"{label} must be text")
    return value


def _bool(value: object, expected: bool | None, label: str) -> bool:
    if type(value) is not bool:
        _fail(f"{label} must be a built-in bool")
    if expected is not None and value is not expected:
        _fail(f"{label} crossed its numerical scope")
    return value


def _digest_text(value: object, label: str) -> str:
    text = _text(value, label)
    if _SHA256.fullmatch(text) is None:
        _fail(f"{label} must be a SHA-256 digest")
    return text


def _nonnegative_int(value: object, label: str) -> int:
    if type(value) is not int or value < 0:
        raise TypeError(f"{label} must be a nonnegative built-in integer")
    return value


def _canonical_bytes(value: object, *, label: str) -> bytes:
    try:
        raw = canonical_json_bytes(value)
    except CanonicalJSONError as error:
        if "byte bound" in str(error):
            raise HLT17IMP1Error(f"{label} exceeds MAX_CURSOR_BYTES") from error
        raise HLT17IMP1Error(f"{label} is not canonical JSON") from error
    if len(raw) > MAX_CURSOR_BYTES:
        raise HLT17IMP1Error(f"{label} exceeds MAX_CURSOR_BYTES")
    return raw


def _parse_hex(value: object, label: str, *, positive: bool = False) -> float:
    text = _text(value, label)
    try:
        number = float.fromhex(text)
    except (OverflowError, ValueError) as error:
        raise HLT17IMP1Error(f"{label} is not canonical binary64 hex") from error
    if not math.isfinite(number) or number.hex() != text:
        _fail(f"{label} is not canonical binary64 hex")
    if positive and not (number > 0.0):
        _fail(f"{label} must be positive")
    return number


def _canonical_hex(value: object, label: str) -> str:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{label} must be a finite real scalar")
    try:
        number = float(value)
    except (OverflowError, ValueError, TypeError) as error:
        raise HLT17IMP1Error(f"{label} must be finite") from error
    if not math.isfinite(number):
        _fail(f"{label} must be finite")
    text = number.hex()
    if float.fromhex(text).hex() != text:
        _fail(f"{label} is not canonical binary64 hex")
    return text


def _same_binary64(left: float, right: float) -> bool:
    return left.hex() == right.hex()


def _exact_mapping(value: object, names: frozenset[str], label: str) -> dict[str, object]:
    if type(value) is not dict or set(value) != names:
        _fail(f"{label} is incomplete or broadened")
    return value


def _binary64_half(width: float) -> float:
    half = width / 2.0
    if not math.isfinite(half) or not (half > 0.0):
        _fail("IMP1 half-cap is not a positive finite binary64 value")
    return half


def ledger_sha256(ledger: TDG11IMP1Ledger) -> str:
    """SHA-256 of the closed IMP1 checkpoint mapping; same formula as IMP1."""

    if not isinstance(ledger, TDG11IMP1Ledger):
        raise TypeError("ledger must be TDG11IMP1Ledger")
    return _digest(_canonical_bytes(imp1_checkpoint_extension(ledger), label="IMP1 ledger"))


@dataclass(frozen=True, slots=True)
class HLT17MemberSpec:
    """Declared member identity. Live keys are never inferred from array size."""

    member_key: str
    method_label: str
    point_count: int
    identity_scope: str

    def __post_init__(self) -> None:
        if type(self.member_key) is not str or not self.member_key:
            _fail("member_key must be nonempty text")
        if self.method_label not in _METHOD_BY_LABEL:
            _fail("unknown HLT17 method label")
        if type(self.point_count) is not int or self.point_count < 9:
            raise TypeError("point_count must be a built-in integer at least 9")
        if self.identity_scope not in IDENTITY_SCOPES:
            _fail("HLT17 identity scope differs")
        expected_key = f"{self.method_label}-{self.point_count}"
        if self.member_key != expected_key:
            _fail("member_key differs from method-label-point_count")
        live_key = self.member_key in {
            "RK4-2049",
            "RK4-4097",
            "RK4-8193",
            "SSPRK3-4097",
            "SSPRK3-8193",
            "SSPRK3-16385",
        }
        if self.identity_scope == IDENTITY_SCOPE_LIVE and not live_key:
            _fail("live HLT17 identity is outside the production set")
        if self.identity_scope == IDENTITY_SCOPE_SYNTHETIC and live_key:
            _fail("synthetic HLT17 identity cannot use a live production key")

    @property
    def method(self) -> str:
        return _METHOD_BY_LABEL[self.method_label]

    @property
    def owned_rows(self) -> int:
        return self.point_count - 5

    @property
    def fallback_d01(self) -> int:
        return 8 * self.owned_rows

    @property
    def fallback_d12(self) -> int:
        return 16 * self.owned_rows

    def enclosure_limits(self) -> IMP1EnclosureLimits:
        return IMP1EnclosureLimits(
            maximum_owned_rows=self.owned_rows,
            fallback_maximum_candidates_D01=self.fallback_d01,
            fallback_maximum_candidates_D12=self.fallback_d12,
        )

    def as_mapping(self) -> dict[str, object]:
        return {
            "member_key": self.member_key,
            "identity_scope": self.identity_scope,
            "method": self.method,
            "method_label": self.method_label,
            "point_count": self.point_count,
            "owned_rows": self.owned_rows,
            "fallback_d01": self.fallback_d01,
            "fallback_d12": self.fallback_d12,
        }


HLT17ProductionMember = HLT17MemberSpec

HLT17_PRODUCTION_MEMBERS: Final[tuple[HLT17MemberSpec, ...]] = (
    HLT17MemberSpec("RK4-2049", "RK4", 2049, IDENTITY_SCOPE_LIVE),
    HLT17MemberSpec("RK4-4097", "RK4", 4097, IDENTITY_SCOPE_LIVE),
    HLT17MemberSpec("RK4-8193", "RK4", 8193, IDENTITY_SCOPE_LIVE),
    HLT17MemberSpec("SSPRK3-4097", "SSPRK3", 4097, IDENTITY_SCOPE_LIVE),
    HLT17MemberSpec("SSPRK3-8193", "SSPRK3", 8193, IDENTITY_SCOPE_LIVE),
    HLT17MemberSpec("SSPRK3-16385", "SSPRK3", 16385, IDENTITY_SCOPE_LIVE),
)
HLT17_SYNTHETIC_MEMBERS: Final[tuple[HLT17MemberSpec, ...]] = (
    HLT17MemberSpec("RK4-9", "RK4", 9, IDENTITY_SCOPE_SYNTHETIC),
    HLT17MemberSpec("SSPRK3-9", "SSPRK3", 9, IDENTITY_SCOPE_SYNTHETIC),
)
HLT17_MEMBER_KEYS: Final[tuple[str, ...]] = tuple(
    item.member_key for item in HLT17_PRODUCTION_MEMBERS
)
HLT17_SYNTHETIC_MEMBER_KEYS: Final[tuple[str, ...]] = tuple(
    item.member_key for item in HLT17_SYNTHETIC_MEMBERS
)
HLT17_MEMBER_BY_KEY: Final[dict[str, HLT17MemberSpec]] = {
    item.member_key: item
    for item in (*HLT17_PRODUCTION_MEMBERS, *HLT17_SYNTHETIC_MEMBERS)
}


def _validate_production_table() -> None:
    if len(HLT17_PRODUCTION_MEMBERS) != 6:
        _fail("HLT17 production member table must contain the six actual members")
    expected = (
        ("RK4-2049", 2049, 2044, 16352, 32704),
        ("RK4-4097", 4097, 4092, 32736, 65472),
        ("RK4-8193", 8193, 8188, 65504, 131008),
        ("SSPRK3-4097", 4097, 4092, 32736, 65472),
        ("SSPRK3-8193", 8193, 8188, 65504, 131008),
        ("SSPRK3-16385", 16385, 16380, 131040, 262080),
    )
    actual = tuple(
        (
            item.member_key,
            item.point_count,
            item.owned_rows,
            item.fallback_d01,
            item.fallback_d12,
        )
        for item in HLT17_PRODUCTION_MEMBERS
    )
    if actual != expected:
        _fail("HLT17 production sizes or 8N/16N ceilings differ")
    if actual[-1][2:] != (16380, 131040, 262080):
        _fail("largest HLT17 member ceilings differ")
    if HLT17_SYNTHETIC_MEMBER_KEYS != ("RK4-9", "SSPRK3-9"):
        _fail("HLT17 synthetic member keys differ")
    for item in HLT17_SYNTHETIC_MEMBERS:
        if item.identity_scope != IDENTITY_SCOPE_SYNTHETIC:
            _fail("synthetic member must use the synthetic identity scope")
        if item.member_key in HLT17_MEMBER_KEYS:
            _fail("synthetic member collided with a production key")


_validate_production_table()


def hlt17_member_spec(member_key: object) -> HLT17MemberSpec:
    key = _text(member_key, "member_key")
    try:
        return HLT17_MEMBER_BY_KEY[key]
    except KeyError as error:
        raise HLT17IMP1Error("unknown HLT17 member identity") from error


def production_member(member_key: object) -> HLT17MemberSpec:
    spec = hlt17_member_spec(member_key)
    if spec.identity_scope != IDENTITY_SCOPE_LIVE:
        _fail("member_key is not a live HLT17 production member")
    return spec


def imp1_enclosure_limits_for_member(member_key: object) -> IMP1EnclosureLimits:
    """Exact per-member IMP1 injection; production ceilings stay unshrunk."""

    return hlt17_member_spec(member_key).enclosure_limits()


def encode_tdg7_plan(plan: TDG7Binary64SubdivisionPlan) -> dict[str, object]:
    """Persist the complete TDG7 plan, including requested_cap ≠ macro_width."""

    if not isinstance(plan, TDG7Binary64SubdivisionPlan):
        raise TypeError("plan must be TDG7Binary64SubdivisionPlan")
    validate_tdg7_binary64_subdivision_plan(plan)
    mapping = {
        "current_hex": plan.current.hex(),
        "event_target_hex": plan.event_target.hex(),
        "requested_cap_hex": plan.requested_cap.hex(),
        "quantum_hex": plan.quantum.hex(),
        "selection_budget_hex": plan.selection_budget.hex(),
        "macro_width_hex": plan.macro_width.hex(),
        "fine_width_hex": plan.fine_width.hex(),
        "conservative_reduction_hex": plan.conservative_reduction.hex(),
        "minimum_width_hex_or_none": (
            None if plan.minimum_width is None else plan.minimum_width.hex()
        ),
        "boundaries_hex": [item.hex() for item in plan.boundaries],
        "target_limited": plan.target_limited,
    }
    if set(mapping) != _PLAN_KEY_SET:
        _fail("TDG7 plan mapping is incomplete or broadened")
    _canonical_bytes(mapping, label="TDG7 plan")
    return mapping


def decode_tdg7_plan(value: object) -> TDG7Binary64SubdivisionPlan:
    item = _exact_mapping(value, _PLAN_KEY_SET, "TDG7 plan")
    boundaries = item["boundaries_hex"]
    if type(boundaries) is not list or len(boundaries) != 5:
        _fail("TDG7 plan must retain all five boundaries")
    minimum = item["minimum_width_hex_or_none"]
    plan = TDG7Binary64SubdivisionPlan(
        current=_parse_hex(item["current_hex"], "current_hex"),
        event_target=_parse_hex(item["event_target_hex"], "event_target_hex"),
        requested_cap=_parse_hex(
            item["requested_cap_hex"], "requested_cap_hex", positive=True
        ),
        quantum=_parse_hex(item["quantum_hex"], "quantum_hex", positive=True),
        selection_budget=_parse_hex(
            item["selection_budget_hex"], "selection_budget_hex", positive=True
        ),
        macro_width=_parse_hex(
            item["macro_width_hex"], "macro_width_hex", positive=True
        ),
        fine_width=_parse_hex(item["fine_width_hex"], "fine_width_hex", positive=True),
        conservative_reduction=_parse_hex(
            item["conservative_reduction_hex"], "conservative_reduction_hex"
        ),
        minimum_width=(
            None
            if minimum is None
            else _parse_hex(minimum, "minimum_width_hex_or_none", positive=True)
        ),
        boundaries=tuple(
            _parse_hex(entry, f"boundaries_hex[{index}]")
            for index, entry in enumerate(boundaries)
        ),
        target_limited=_bool(item["target_limited"], None, "target_limited"),
    )
    if encode_tdg7_plan(plan) != item:
        _fail("TDG7 plan identities are not redundant")
    return plan


def plan_imp1_lattice(
    current: object, event_target: object, requested_cap: object
) -> TDG7Binary64SubdivisionPlan:
    """Certified IMP1/TDG7 lattice; requested_cap is not silently replaced by width."""

    try:
        return plan_forward_proto14_subdivision(
            current,
            event_target,
            requested_cap,
            minimum_width=float(TDG6_MINIMUM_MACRO_STEP),
        )
    except CoordinateLatticeLimitReached as error:
        raise HLT17OwnerRetryExhausted(
            TEMPORAL_OWNER, error.reason
        ) from error


def temporal_successor_plan(
    predecessor: TDG7Binary64SubdivisionPlan,
) -> TDG7Binary64SubdivisionPlan:
    """Exact IMP1 half-cap successor of one rejected proposal."""

    validate_tdg7_binary64_subdivision_plan(predecessor)
    try:
        successor = plan_imp1_lattice(
            predecessor.current,
            predecessor.event_target,
            _binary64_half(predecessor.macro_width),
        )
    except HLT17OwnerRetryExhausted as error:
        raise HLT17OwnerRetryExhausted(
            TEMPORAL_OWNER,
            error.reason,
            attempted_plan=predecessor,
        ) from error
    if not (successor.macro_width < predecessor.macro_width):
        raise HLT17OwnerRetryExhausted(
            TEMPORAL_OWNER,
            "not_strictly_smaller",
            attempted_plan=predecessor,
        )
    if not _same_binary64(
        successor.requested_cap, _binary64_half(predecessor.macro_width)
    ):
        _fail("IMP1 temporal successor requested_cap is not the exact half-cap")
    return successor


def strictly_smaller_tdg7_plan(
    plan: TDG7Binary64SubdivisionPlan, *, owner: str
) -> TDG7Binary64SubdivisionPlan:
    """Next source/CFL lattice plan; exhaustion belongs to that owner."""

    if owner not in {SOURCE_OWNER, CFL_OWNER}:
        _fail("strictly smaller lattice plans are owned by source or CFL")
    validate_tdg7_binary64_subdivision_plan(plan)
    try:
        successor = plan_imp1_lattice(
            plan.current,
            plan.event_target,
            _binary64_half(plan.macro_width),
        )
    except HLT17OwnerRetryExhausted as error:
        raise HLT17OwnerRetryExhausted(
            owner, error.reason, attempted_plan=plan
        ) from error
    if not (successor.macro_width < plan.macro_width):
        raise HLT17OwnerRetryExhausted(
            owner, "not_strictly_smaller", attempted_plan=plan
        )
    return successor


def _plan_or_none(value: object) -> TDG7Binary64SubdivisionPlan | None:
    if value is None:
        return None
    return decode_tdg7_plan(value)


def _rejection_from_mapping(value: object) -> TDG11IMP1TemporalRejection:
    if not isinstance(value, Mapping):
        raise TypeError("IMP1 rejection must be a mapping")
    try:
        return TDG11IMP1TemporalRejection.from_mapping(dict(value))
    except (TypeError, TDG11IMP1Error, ValueError) as error:
        raise HLT17IMP1Error("IMP1 rejection evidence is malformed") from error


def _ledger_from_checkpoint(value: object, label: str) -> TDG11IMP1Ledger:
    if not isinstance(value, Mapping):
        raise TypeError(f"{label} must be an IMP1 checkpoint mapping")
    try:
        return restore_imp1_checkpoint_extension(dict(value))
    except (TypeError, TDG11IMP1Error, ValueError) as error:
        raise HLT17IMP1Error(f"{label} is not a closed IMP1 checkpoint") from error


def _overlay_evidence(value: object, owner: str) -> dict[str, object]:
    names = frozenset(SOURCE_EVIDENCE_KEYS if owner == SOURCE_OWNER else CFL_EVIDENCE_KEYS)
    item = _exact_mapping(value, names, "overlay evidence")
    kind = _text(item["kind"], "overlay evidence kind")
    expected = "source_retry" if owner == SOURCE_OWNER else "cfl_retry"
    if kind != expected:
        _fail("overlay evidence kind differs from owner")
    _text(item["path"], "overlay evidence path")
    attempted_cap = _parse_hex(
        item["attempted_requested_cap_hex"],
        "attempted_requested_cap_hex",
        positive=True,
    )
    attempted_width = _parse_hex(
        item["attempted_macro_width_hex"],
        "attempted_macro_width_hex",
        positive=True,
    )
    if attempted_cap < attempted_width:
        _fail("overlay attempted requested_cap is smaller than executed width")
    _bool(item["accepted_boundary_preserved"], True, "accepted_boundary_preserved")
    _bool(item["physical_classification"], False, "physical_classification")
    if owner == SOURCE_OWNER:
        if type(item["source_retry"]) is not dict:
            _fail("source overlay must retain wired source-retry evidence")
        _canonical_bytes(item["source_retry"], label="source retry evidence")
    else:
        _parse_hex(item["observed_ratio_hex"], "observed_ratio_hex", positive=True)
        _parse_hex(item["maximum_ratio_hex"], "maximum_ratio_hex", positive=True)
    return item


@dataclass(frozen=True, slots=True)
class HLT17TemporalIdentity:
    """Exact preceding IMP1 rejection and its half-cap successor identity."""

    rejection: TDG11IMP1TemporalRejection
    predecessor_plan: TDG7Binary64SubdivisionPlan
    successor_plan: TDG7Binary64SubdivisionPlan
    predecessor_ledger: TDG11IMP1Ledger
    preparation_sha256: str
    assessment_sha256: str

    def __post_init__(self) -> None:
        if not isinstance(self.rejection, TDG11IMP1TemporalRejection):
            raise TypeError("rejection must be TDG11IMP1TemporalRejection")
        validate_tdg7_binary64_subdivision_plan(self.predecessor_plan)
        validate_tdg7_binary64_subdivision_plan(self.successor_plan)
        if not isinstance(self.predecessor_ledger, TDG11IMP1Ledger):
            raise TypeError("predecessor_ledger must be TDG11IMP1Ledger")
        object.__setattr__(
            self,
            "preparation_sha256",
            _digest_text(self.preparation_sha256, "preparation_sha256"),
        )
        object.__setattr__(
            self,
            "assessment_sha256",
            _digest_text(self.assessment_sha256, "assessment_sha256"),
        )
        if self.rejection.preparation_sha256 != self.preparation_sha256:
            _fail("temporal rejection preparation hash differs")
        if self.rejection.assessment_sha256 != self.assessment_sha256:
            _fail("temporal rejection assessment hash differs")
        if not _same_binary64(
            self.predecessor_plan.current,
            _parse_hex(self.rejection.time_hex, "time_hex"),
        ):
            _fail("predecessor plan current time differs from the IMP1 rejection")
        if not _same_binary64(
            self.predecessor_plan.event_target,
            _parse_hex(self.rejection.event_target_hex, "event_target_hex"),
        ):
            _fail("predecessor plan event target differs from the IMP1 rejection")
        if not _same_binary64(
            self.predecessor_plan.macro_width,
            _parse_hex(self.rejection.attempted_width_hex, "attempted_width_hex"),
        ):
            _fail("predecessor executed width differs from the IMP1 rejection")
        expected = temporal_successor_plan(self.predecessor_plan)
        if self.successor_plan != expected:
            _fail("temporal successor is not the exact IMP1 half-cap plan")
        if not _same_binary64(
            self.successor_plan.requested_cap,
            _parse_hex(self.rejection.next_cap_hex, "next_cap_hex"),
        ):
            _fail("IMP1 half-cap identity differs from the successor requested_cap")
        append_imp1_rejection(self.predecessor_ledger, self.rejection.as_mapping())

    @property
    def updated_ledger(self) -> TDG11IMP1Ledger:
        return append_imp1_rejection(self.predecessor_ledger, self.rejection.as_mapping())

    @property
    def imp1_half_cap(self) -> float:
        return self.successor_plan.requested_cap

    def as_mapping(self) -> dict[str, object]:
        predecessor_checkpoint = imp1_checkpoint_extension(self.predecessor_ledger)
        mapping = {
            "rejection": self.rejection.as_mapping(),
            "predecessor_plan": encode_tdg7_plan(self.predecessor_plan),
            "successor_plan": encode_tdg7_plan(self.successor_plan),
            "predecessor_imp1_ledger": predecessor_checkpoint,
            "predecessor_imp1_ledger_sha256": ledger_sha256(self.predecessor_ledger),
            "preparation_sha256": self.preparation_sha256,
            "assessment_sha256": self.assessment_sha256,
            "rejection_sha256": _digest(
                self.rejection.canonical_text().encode("ascii")
            ),
            "imp1_half_cap_hex": self.imp1_half_cap.hex(),
        }
        if set(mapping) != _TEMPORAL_KEY_SET:
            _fail("temporal identity is incomplete or broadened")
        return mapping

    @classmethod
    def from_mapping(cls, value: object) -> "HLT17TemporalIdentity":
        item = _exact_mapping(value, _TEMPORAL_KEY_SET, "temporal identity")
        identity = cls(
            rejection=_rejection_from_mapping(item["rejection"]),
            predecessor_plan=decode_tdg7_plan(item["predecessor_plan"]),
            successor_plan=decode_tdg7_plan(item["successor_plan"]),
            predecessor_ledger=_ledger_from_checkpoint(
                item["predecessor_imp1_ledger"], "predecessor_imp1_ledger"
            ),
            preparation_sha256=item["preparation_sha256"],
            assessment_sha256=item["assessment_sha256"],
        )
        if identity.as_mapping() != item:
            _fail("temporal identity is not a redundant round trip")
        if identity.as_mapping()["predecessor_imp1_ledger_sha256"] != _digest_text(
            item["predecessor_imp1_ledger_sha256"], "predecessor_imp1_ledger_sha256"
        ):
            _fail("predecessor IMP1 ledger hash differs")
        if identity.as_mapping()["rejection_sha256"] != _digest_text(
            item["rejection_sha256"], "rejection_sha256"
        ):
            _fail("IMP1 rejection hash differs")
        if identity.imp1_half_cap.hex() != _text(
            item["imp1_half_cap_hex"], "imp1_half_cap_hex"
        ):
            _fail("IMP1 half-cap hex differs")
        return identity


@dataclass(frozen=True, slots=True)
class HLT17RetryOverlay:
    """Source/CFL reduction that does not replace the IMP1 half-cap successor."""

    owner: str
    attempted_plan: TDG7Binary64SubdivisionPlan
    next_plan: TDG7Binary64SubdivisionPlan | None
    evidence: dict[str, object]
    owner_current_count: int
    exhausted: bool

    def __post_init__(self) -> None:
        if self.owner not in {SOURCE_OWNER, CFL_OWNER}:
            _fail("overlay owner must be source or CFL, never temporal")
        validate_tdg7_binary64_subdivision_plan(self.attempted_plan)
        object.__setattr__(
            self,
            "owner_current_count",
            _nonnegative_int(self.owner_current_count, "owner_current_count"),
        )
        object.__setattr__(
            self, "exhausted", _bool(self.exhausted, None, "exhausted")
        )
        if self.owner_current_count < 1:
            _fail("overlay owner current count must be at least one")
        evidence = _overlay_evidence(self.evidence, self.owner)
        object.__setattr__(self, "evidence", evidence)
        if not _same_binary64(
            self.attempted_plan.requested_cap,
            _parse_hex(
                evidence["attempted_requested_cap_hex"],
                "attempted_requested_cap_hex",
            ),
        ) or not _same_binary64(
            self.attempted_plan.macro_width,
            _parse_hex(
                evidence["attempted_macro_width_hex"], "attempted_macro_width_hex"
            ),
        ):
            _fail("overlay evidence does not match the actually attempted plan")
        if self.exhausted:
            if self.next_plan is not None:
                _fail("exhausted overlay must not advertise a next lattice plan")
            return
        if self.next_plan is None:
            _fail("non-exhausted overlay must retain a strictly smaller next plan")
        validate_tdg7_binary64_subdivision_plan(self.next_plan)
        expected = strictly_smaller_tdg7_plan(self.attempted_plan, owner=self.owner)
        if self.next_plan != expected:
            _fail("overlay next plan is not the certified strictly smaller lattice")

    def as_mapping(self) -> dict[str, object]:
        mapping = {
            "owner": self.owner,
            "attempted_plan": encode_tdg7_plan(self.attempted_plan),
            "next_plan": (
                None if self.next_plan is None else encode_tdg7_plan(self.next_plan)
            ),
            "evidence": dict(self.evidence),
            "owner_current_count": self.owner_current_count,
            "exhausted": self.exhausted,
        }
        if set(mapping) != _OVERLAY_KEY_SET:
            _fail("overlay mapping is incomplete or broadened")
        return mapping

    @classmethod
    def from_mapping(cls, value: object) -> "HLT17RetryOverlay":
        item = _exact_mapping(value, _OVERLAY_KEY_SET, "retry overlay")
        overlay = cls(
            owner=_text(item["owner"], "owner"),
            attempted_plan=decode_tdg7_plan(item["attempted_plan"]),
            next_plan=_plan_or_none(item["next_plan"]),
            evidence=item["evidence"],
            owner_current_count=item["owner_current_count"],
            exhausted=item["exhausted"],
        )
        if overlay.as_mapping() != item:
            _fail("retry overlay is not a redundant round trip")
        return overlay


@dataclass(frozen=True, slots=True)
class HLT17IMP1Cursor:
    """Closed immutable HLT17 cursor; descriptor and u/p/q hashes are distinct.

    ``kernel_transition_parent_sha256`` / ``kernel_transition_digest`` hash
    this cursor transition. They are not a store-journal parent/tip.
    ``checkpoint_generation`` is a kernel accepted-boundary revision, not a
    FUTURE store checkpoint revision.
    """

    member_key: str
    identity_scope: str
    implementation_id: str
    mathematical_object: str
    reference_wire_artifact_id: str
    reference_wire_evaluator_id: str
    admission_runtime_id: str
    mode: str
    event_target: float
    accepted_time: float
    accepted_step_index: int
    accepted_transaction_serial: int
    descriptor_sha256: str
    physical_state_sha256: str
    coordinates_sha256: str
    synthetic_callable_binding_sha256: str
    monitor_sha256: str
    causal_state_sha256: str
    tracer_sha256: str
    guard_configuration_sha256: str
    ledger: TDG11IMP1Ledger
    temporal: HLT17TemporalIdentity | None
    overlay: HLT17RetryOverlay | None
    source_current: int
    source_total: int
    cfl_current: int
    cfl_total: int
    cursor_generation: int
    attempt_serial: int
    checkpoint_generation: int
    kernel_transition_parent_sha256: str
    kernel_transition_digest: str
    # kernel_transition_* hash this cursor transition. They are not a store
    # journal parent/tip and do not bind a durable checkpoint revision.

    def __post_init__(self) -> None:
        member = hlt17_member_spec(self.member_key)
        if self.identity_scope != member.identity_scope:
            _fail("cursor identity_scope differs from the member specification")
        try:
            require_hlt17_implementation_identity(
                {
                    "implementation_id": self.implementation_id,
                    "mathematical_object": self.mathematical_object,
                    "reference_wire_artifact_id": self.reference_wire_artifact_id,
                    "reference_wire_evaluator_id": self.reference_wire_evaluator_id,
                    "admission_runtime_id": self.admission_runtime_id,
                },
                label="HLT17 cursor",
            )
        except (TypeError, ValueError) as error:
            raise HLT17IMP1Error(str(error)) from error
        if self.mode not in {FRESH_READY, RETRY_PENDING}:
            _fail("cursor mode must be FRESH_READY or RETRY_PENDING")
        object.__setattr__(
            self,
            "event_target",
            _parse_hex(_canonical_hex(self.event_target, "event_target"), "event_target"),
        )
        object.__setattr__(
            self,
            "accepted_time",
            _parse_hex(
                _canonical_hex(self.accepted_time, "accepted_time"), "accepted_time"
            ),
        )
        object.__setattr__(
            self,
            "accepted_step_index",
            _nonnegative_int(self.accepted_step_index, "accepted_step_index"),
        )
        object.__setattr__(
            self,
            "accepted_transaction_serial",
            _nonnegative_int(
                self.accepted_transaction_serial, "accepted_transaction_serial"
            ),
        )
        for name in (
            "descriptor_sha256",
            "physical_state_sha256",
            "coordinates_sha256",
            "synthetic_callable_binding_sha256",
            "monitor_sha256",
            "causal_state_sha256",
            "tracer_sha256",
            "guard_configuration_sha256",
            "kernel_transition_parent_sha256",
            "kernel_transition_digest",
        ):
            object.__setattr__(self, name, _digest_text(getattr(self, name), name))
        if self.descriptor_sha256 == self.physical_state_sha256:
            _fail("descriptor and physical u/p/q identities are not interchangeable")
        if not isinstance(self.ledger, TDG11IMP1Ledger):
            raise TypeError("ledger must be TDG11IMP1Ledger")
        if self.ledger.method != member.method:
            _fail("IMP1 ledger method differs from the HLT17 member")
        if IMP1_METHOD_LABEL[self.ledger.method] != member.method_label:
            _fail("IMP1 method label differs from the HLT17 member")
        if not _same_binary64(self.ledger.last_accepted_time, self.accepted_time):
            _fail("accepted time differs from the IMP1 ledger")
        if self.ledger.current_step_index != self.accepted_step_index:
            _fail("accepted step index differs from the IMP1 ledger")
        if self.ledger.current_transaction_serial != self.accepted_transaction_serial:
            _fail("accepted transaction serial differs from the IMP1 ledger")
        if self.ledger.current_state_sha256 != self.physical_state_sha256:
            _fail("physical u/p/q hash differs from the IMP1 ledger state identity")
        if not (
            _same_binary64(self.event_target, self.accepted_time)
            or self.event_target > self.accepted_time
        ):
            _fail("event target must not precede the accepted time")
        for name in (
            "source_current",
            "source_total",
            "cfl_current",
            "cfl_total",
            "cursor_generation",
            "attempt_serial",
            "checkpoint_generation",
        ):
            object.__setattr__(self, name, _nonnegative_int(getattr(self, name), name))
        if self.source_current > self.source_total or self.cfl_current > self.cfl_total:
            _fail("source or CFL current count exceeds its cumulative total")
        if (
            self.source_current > HLT17_SOURCE_RETRY_CAP + 1
            or self.cfl_current > HLT17_CFL_RETRY_CAP + 1
        ):
            _fail("source or CFL current count exceeds the exhaustion edge")
        if self.temporal is not None:
            if type(self.temporal) is not HLT17TemporalIdentity:
                raise TypeError("temporal must be HLT17TemporalIdentity or None")
            self.temporal.__post_init__()
        if self.overlay is not None:
            if type(self.overlay) is not HLT17RetryOverlay:
                raise TypeError("overlay must be HLT17RetryOverlay or None")
            self.overlay.__post_init__()
        if self.mode == FRESH_READY:
            if self.temporal is not None:
                _fail("FRESH_READY must not fabricate temporal rejection evidence")
            if self.ledger.current_macro_step_temporal_retry_count != 0:
                _fail("FRESH_READY cannot carry an active IMP1 temporal retry")
        else:
            if self.temporal is None:
                _fail("RETRY_PENDING must retain the exact IMP1 temporal identity")
            if ledger_sha256(self.temporal.updated_ledger) != ledger_sha256(self.ledger):
                _fail("RETRY_PENDING ledger is not the updated IMP1 rejection ledger")
            if self.temporal.rejection.method != self.ledger.method:
                _fail("temporal rejection method differs")
            if self.temporal.rejection.accepted_state_sha256 != self.physical_state_sha256:
                _fail("temporal rejection accepted state differs from physical u/p/q")
            if self.temporal.rejection.accepted_step_index != self.accepted_step_index:
                _fail("temporal rejection step index differs")
            if (
                self.temporal.rejection.accepted_transaction_serial
                != self.accepted_transaction_serial
            ):
                _fail("temporal rejection transaction serial differs")
            if not _same_binary64(
                self.temporal.predecessor_plan.event_target, self.event_target
            ):
                _fail("temporal predecessor target differs from the cursor target")
            if self.ledger.current_macro_step_temporal_retry_count < 1:
                _fail("RETRY_PENDING IMP1 current temporal retry count is zero")
        if self.overlay is not None:
            attempted = self.overlay.attempted_plan
            if not _same_binary64(attempted.current, self.accepted_time):
                _fail("overlay attempted plan is not bound to the accepted time")
            if not _same_binary64(attempted.event_target, self.event_target):
                _fail("overlay attempted plan is not bound to the event target")
            if self.mode == RETRY_PENDING:
                assert self.temporal is not None
                successor_width = self.temporal.successor_plan.macro_width
                attempted_width = self.overlay.attempted_plan.macro_width
                if self.overlay.attempted_plan == self.temporal.predecessor_plan:
                    _fail("source/CFL overlay must not replace the IMP1 successor plan")
                if not (
                    _same_binary64(attempted_width, successor_width)
                    or attempted_width < successor_width
                ):
                    _fail("overlay attempted width exceeds the IMP1 half-cap successor")
            owner_current = (
                self.source_current
                if self.overlay.owner == SOURCE_OWNER
                else self.cfl_current
            )
            if self.overlay.owner_current_count != owner_current:
                _fail("overlay independent counter differs from the owner current count")
            cap = (
                HLT17_SOURCE_RETRY_CAP
                if self.overlay.owner == SOURCE_OWNER
                else HLT17_CFL_RETRY_CAP
            )
            if owner_current > cap and not self.overlay.exhausted:
                _fail("owner current count exceeded its cap without exhaustion")

    @property
    def member(self) -> HLT17MemberSpec:
        return hlt17_member_spec(self.member_key)

    @property
    def method(self) -> str:
        return self.member.method

    @property
    def latest_owner(self) -> str | None:
        if self.overlay is not None:
            return self.overlay.owner
        if self.mode == RETRY_PENDING:
            return TEMPORAL_OWNER
        return None

    @property
    def active_temporal_retry_depth(self) -> int:
        return self.ledger.current_macro_step_temporal_retry_count

    def next_attempt_plan(self) -> TDG7Binary64SubdivisionPlan | None:
        """Plan the next attempt should execute; IMP1 successor stays untouched."""

        if self.overlay is not None:
            return self.overlay.next_plan
        if self.temporal is not None:
            return self.temporal.successor_plan
        return None

    def public_fresh_permitted(self) -> bool:
        return (
            self.mode == FRESH_READY
            and self.overlay is None
            and self.active_temporal_retry_depth == 0
        )


def encode_hlt17_cursor(
    cursor: HLT17IMP1Cursor, *, recompute_tip: bool = True
) -> dict[str, object]:
    if not isinstance(cursor, HLT17IMP1Cursor):
        raise TypeError("cursor must be HLT17IMP1Cursor")
    member = cursor.member
    checkpoint = imp1_checkpoint_extension(cursor.ledger)
    mapping: dict[str, object] = {
        "schema": HLT17_CURSOR_SCHEMA,
        "schema_version": HLT17_CURSOR_SCHEMA_VERSION,
        **_implementation_fields(cursor),
        "member_key": member.member_key,
        "identity_scope": member.identity_scope,
        "method": member.method,
        "method_label": member.method_label,
        "point_count": member.point_count,
        "owned_rows": member.owned_rows,
        "fallback_d01": member.fallback_d01,
        "fallback_d12": member.fallback_d12,
        "mode": cursor.mode,
        "event_target_hex": cursor.event_target.hex(),
        "accepted_time_hex": cursor.accepted_time.hex(),
        "accepted_step_index": cursor.accepted_step_index,
        "accepted_transaction_serial": cursor.accepted_transaction_serial,
        "descriptor_sha256": cursor.descriptor_sha256,
        "physical_state_sha256": cursor.physical_state_sha256,
        "coordinates_sha256": cursor.coordinates_sha256,
        "synthetic_callable_binding_sha256": cursor.synthetic_callable_binding_sha256,
        "monitor_sha256": cursor.monitor_sha256,
        "causal_state_sha256": cursor.causal_state_sha256,
        "tracer_sha256": cursor.tracer_sha256,
        "guard_configuration_sha256": cursor.guard_configuration_sha256,
        "imp1_ledger": checkpoint,
        "imp1_ledger_sha256": ledger_sha256(cursor.ledger),
        "inherited_snapshot_sha256": checkpoint["inherited_snapshot_sha256"],
        "temporal": None if cursor.temporal is None else cursor.temporal.as_mapping(),
        "overlay": None if cursor.overlay is None else cursor.overlay.as_mapping(),
        "source_current": cursor.source_current,
        "source_total": cursor.source_total,
        "cfl_current": cursor.cfl_current,
        "cfl_total": cursor.cfl_total,
        "source_cap": HLT17_SOURCE_RETRY_CAP,
        "cfl_cap": HLT17_CFL_RETRY_CAP,
        "cursor_generation": cursor.cursor_generation,
        "attempt_serial": cursor.attempt_serial,
        "checkpoint_generation": cursor.checkpoint_generation,
        "kernel_transition_parent_sha256": cursor.kernel_transition_parent_sha256,
        "kernel_transition_digest": cursor.kernel_transition_digest,
        "historical_origin_authenticated": False,
        "physical_classification": False,
        "campaign_execution_authorized": False,
        "origin_authenticated_by_this_cursor": False,
        "source_authenticated_by_this_cursor": False,
        "store_authenticated_by_this_cursor": False,
    }
    if set(mapping) != _CURSOR_KEY_SET:
        _fail("HLT17 cursor mapping is incomplete or broadened")
    if "pending_owner" in mapping:
        _fail("HLT17 cursor must not use exclusive pending_owner")
    if recompute_tip:
        expected = _digest(
            _canonical_bytes(
                {**mapping, "kernel_transition_digest": HLT17_GENESIS_DIGEST},
                label="HLT17 cursor",
            )
        )
        mapping["kernel_transition_digest"] = expected
    _canonical_bytes(mapping, label="HLT17 cursor")
    return mapping


def restore_hlt17_cursor(value: object) -> HLT17IMP1Cursor:
    item = _exact_mapping(value, _CURSOR_KEY_SET, "HLT17 cursor")
    if item["schema"] == HLT17_LEGACY_IMP1_CURSOR_SCHEMA:
        _fail("legacy IMP1 HLT17 cursor schema cannot be restored as C1R1")
    if item["schema"] != HLT17_CURSOR_SCHEMA:
        _fail("HLT17 cursor schema differs")
    if item["schema_version"] != HLT17_CURSOR_SCHEMA_VERSION:
        _fail("HLT17 cursor schema version differs")
    try:
        require_hlt17_implementation_identity(
            {name: item[name] for name in IMPLEMENTATION_IDENTITY_KEYS},
            label="HLT17 cursor",
        )
    except (TypeError, ValueError) as error:
        raise HLT17IMP1Error(str(error)) from error
    member = hlt17_member_spec(item["member_key"])
    if item["identity_scope"] != member.identity_scope:
        _fail("cursor identity_scope differs from the member specification")
    if item["method"] != member.method or item["method_label"] != member.method_label:
        _fail("cursor method identity differs from the HLT17 member")
    if (
        item["point_count"] != member.point_count
        or item["owned_rows"] != member.owned_rows
        or item["fallback_d01"] != member.fallback_d01
        or item["fallback_d12"] != member.fallback_d12
    ):
        _fail("cursor member sizes or 8N/16N ceilings differ")
    if member.identity_scope == IDENTITY_SCOPE_LIVE and member.member_key not in HLT17_MEMBER_KEYS:
        _fail("live cursor is outside the frozen six-member production table")
    if member.identity_scope == IDENTITY_SCOPE_SYNTHETIC and member.member_key in HLT17_MEMBER_KEYS:
        _fail("synthetic cursor cannot use a production member key")
    _bool(item["historical_origin_authenticated"], False, "historical_origin_authenticated")
    _bool(item["physical_classification"], False, "physical_classification")
    _bool(
        item["campaign_execution_authorized"], False, "campaign_execution_authorized"
    )
    _bool(
        item["origin_authenticated_by_this_cursor"],
        False,
        "origin_authenticated_by_this_cursor",
    )
    _bool(
        item["source_authenticated_by_this_cursor"],
        False,
        "source_authenticated_by_this_cursor",
    )
    _bool(
        item["store_authenticated_by_this_cursor"],
        False,
        "store_authenticated_by_this_cursor",
    )
    if item["source_cap"] != HLT17_SOURCE_RETRY_CAP or item["cfl_cap"] != HLT17_CFL_RETRY_CAP:
        _fail("source/CFL retry caps differ from the frozen HLT17 contract")
    temporal_value = item["temporal"]
    overlay_value = item["overlay"]
    cursor = HLT17IMP1Cursor(
        member_key=member.member_key,
        identity_scope=member.identity_scope,
        **_implementation_fields(),
        mode=_text(item["mode"], "mode"),
        event_target=_parse_hex(item["event_target_hex"], "event_target_hex"),
        accepted_time=_parse_hex(item["accepted_time_hex"], "accepted_time_hex"),
        accepted_step_index=item["accepted_step_index"],
        accepted_transaction_serial=item["accepted_transaction_serial"],
        descriptor_sha256=item["descriptor_sha256"],
        physical_state_sha256=item["physical_state_sha256"],
        coordinates_sha256=item["coordinates_sha256"],
        synthetic_callable_binding_sha256=item["synthetic_callable_binding_sha256"],
        monitor_sha256=item["monitor_sha256"],
        causal_state_sha256=item["causal_state_sha256"],
        tracer_sha256=item["tracer_sha256"],
        guard_configuration_sha256=item["guard_configuration_sha256"],
        ledger=_ledger_from_checkpoint(item["imp1_ledger"], "imp1_ledger"),
        temporal=(
            None if temporal_value is None else HLT17TemporalIdentity.from_mapping(temporal_value)
        ),
        overlay=(
            None if overlay_value is None else HLT17RetryOverlay.from_mapping(overlay_value)
        ),
        source_current=item["source_current"],
        source_total=item["source_total"],
        cfl_current=item["cfl_current"],
        cfl_total=item["cfl_total"],
        cursor_generation=item["cursor_generation"],
        attempt_serial=item["attempt_serial"],
        checkpoint_generation=item["checkpoint_generation"],
        kernel_transition_parent_sha256=item["kernel_transition_parent_sha256"],
        kernel_transition_digest=item["kernel_transition_digest"],
    )
    encoded = encode_hlt17_cursor(cursor)
    if encoded != item:
        _fail("HLT17 cursor identities are not redundant")
    if encoded["imp1_ledger_sha256"] != _digest_text(
        item["imp1_ledger_sha256"], "imp1_ledger_sha256"
    ):
        _fail("IMP1 ledger hash differs")
    if encoded["inherited_snapshot_sha256"] != _digest_text(
        item["inherited_snapshot_sha256"], "inherited_snapshot_sha256"
    ):
        _fail("inherited snapshot hash differs")
    return cursor


def cursor_sha256(cursor: HLT17IMP1Cursor) -> str:
    return _digest(_canonical_bytes(encode_hlt17_cursor(cursor), label="HLT17 cursor"))


def revalidate_hlt17_cursor(cursor: HLT17IMP1Cursor) -> HLT17IMP1Cursor:
    """Re-check the complete cursor, including nested temporal/overlay evidence."""

    if type(cursor) is not HLT17IMP1Cursor:
        raise TypeError("cursor must be HLT17IMP1Cursor")
    cursor.__post_init__()
    if cursor.temporal is not None:
        if type(cursor.temporal) is not HLT17TemporalIdentity:
            _fail("temporal identity has the wrong type")
        cursor.temporal.__post_init__()
    if cursor.overlay is not None:
        if type(cursor.overlay) is not HLT17RetryOverlay:
            _fail("overlay has the wrong type")
        cursor.overlay.__post_init__()
        mapped = cursor.overlay.as_mapping()
        if HLT17RetryOverlay.from_mapping(mapped).as_mapping() != mapped:
            _fail("overlay evidence is not a redundant round trip")
    restored = restore_hlt17_cursor(encode_hlt17_cursor(cursor))
    if encode_hlt17_cursor(restored) != encode_hlt17_cursor(cursor):
        _fail("HLT17 cursor failed complete nested revalidation")
    return restored


def _with_tip(
    cursor: HLT17IMP1Cursor, *, parent: str, generation: int, serial: int, checkpoint: int
) -> HLT17IMP1Cursor:
    provisional = HLT17IMP1Cursor(
        member_key=cursor.member_key,
        identity_scope=cursor.identity_scope,
        **_implementation_fields(cursor),
        mode=cursor.mode,
        event_target=cursor.event_target,
        accepted_time=cursor.accepted_time,
        accepted_step_index=cursor.accepted_step_index,
        accepted_transaction_serial=cursor.accepted_transaction_serial,
        descriptor_sha256=cursor.descriptor_sha256,
        physical_state_sha256=cursor.physical_state_sha256,
        coordinates_sha256=cursor.coordinates_sha256,
        synthetic_callable_binding_sha256=cursor.synthetic_callable_binding_sha256,
        monitor_sha256=cursor.monitor_sha256,
        causal_state_sha256=cursor.causal_state_sha256,
        tracer_sha256=cursor.tracer_sha256,
        guard_configuration_sha256=cursor.guard_configuration_sha256,
        ledger=cursor.ledger,
        temporal=cursor.temporal,
        overlay=cursor.overlay,
        source_current=cursor.source_current,
        source_total=cursor.source_total,
        cfl_current=cursor.cfl_current,
        cfl_total=cursor.cfl_total,
        cursor_generation=generation,
        attempt_serial=serial,
        checkpoint_generation=checkpoint,
        kernel_transition_parent_sha256=parent,
        kernel_transition_digest=HLT17_GENESIS_DIGEST,
    )
    mapping = encode_hlt17_cursor(provisional, recompute_tip=False)
    tip = _digest(
        _canonical_bytes(
            {**mapping, "kernel_transition_digest": HLT17_GENESIS_DIGEST},
            label="HLT17 cursor tip",
        )
    )
    return HLT17IMP1Cursor(
        member_key=provisional.member_key,
        identity_scope=provisional.identity_scope,
        **_implementation_fields(provisional),
        mode=provisional.mode,
        event_target=provisional.event_target,
        accepted_time=provisional.accepted_time,
        accepted_step_index=provisional.accepted_step_index,
        accepted_transaction_serial=provisional.accepted_transaction_serial,
        descriptor_sha256=provisional.descriptor_sha256,
        physical_state_sha256=provisional.physical_state_sha256,
        coordinates_sha256=provisional.coordinates_sha256,
        synthetic_callable_binding_sha256=provisional.synthetic_callable_binding_sha256,
        monitor_sha256=provisional.monitor_sha256,
        causal_state_sha256=provisional.causal_state_sha256,
        tracer_sha256=provisional.tracer_sha256,
        guard_configuration_sha256=provisional.guard_configuration_sha256,
        ledger=provisional.ledger,
        temporal=provisional.temporal,
        overlay=provisional.overlay,
        source_current=provisional.source_current,
        source_total=provisional.source_total,
        cfl_current=provisional.cfl_current,
        cfl_total=provisional.cfl_total,
        cursor_generation=generation,
        attempt_serial=serial,
        checkpoint_generation=checkpoint,
        kernel_transition_parent_sha256=parent,
        kernel_transition_digest=tip,
    )


def seed_hlt17_cursor(
    *,
    member_key: str,
    event_target: object,
    ledger: TDG11IMP1Ledger,
    descriptor_sha256: str,
    physical_state_sha256: str,
    coordinates_sha256: str,
    synthetic_callable_binding_sha256: str,
    monitor_sha256: str,
    causal_state_sha256: str,
    tracer_sha256: str,
    guard_configuration_sha256: str,
    source_total: int = 0,
    cfl_total: int = 0,
) -> HLT17IMP1Cursor:
    """Fresh cursor with no temporal rejection and no source/CFL overlay.

    ``source_total`` and ``cfl_total`` are inherited cumulative baselines.
    They are not zeroed at seed. Per-macro current counts start at zero.
    ``descriptor_sha256`` must be the actual accepted-encoding digest.
    """

    if not isinstance(ledger, TDG11IMP1Ledger):
        raise TypeError("ledger must be TDG11IMP1Ledger")
    if ledger.current_macro_step_temporal_retry_count != 0:
        _fail("fresh HLT17 seed cannot consume an active IMP1 retry ledger")
    spec = hlt17_member_spec(member_key)
    inherited_source = _nonnegative_int(source_total, "source_total")
    inherited_cfl = _nonnegative_int(cfl_total, "cfl_total")
    draft = HLT17IMP1Cursor(
        member_key=spec.member_key,
        identity_scope=spec.identity_scope,
        **_implementation_fields(),
        mode=FRESH_READY,
        event_target=_parse_hex(_canonical_hex(event_target, "event_target"), "event_target"),
        accepted_time=ledger.last_accepted_time,
        accepted_step_index=ledger.current_step_index,
        accepted_transaction_serial=ledger.current_transaction_serial,
        descriptor_sha256=descriptor_sha256,
        physical_state_sha256=physical_state_sha256,
        coordinates_sha256=coordinates_sha256,
        synthetic_callable_binding_sha256=synthetic_callable_binding_sha256,
        monitor_sha256=monitor_sha256,
        causal_state_sha256=causal_state_sha256,
        tracer_sha256=tracer_sha256,
        guard_configuration_sha256=guard_configuration_sha256,
        ledger=ledger,
        temporal=None,
        overlay=None,
        source_current=0,
        source_total=inherited_source,
        cfl_current=0,
        cfl_total=inherited_cfl,
        cursor_generation=0,
        attempt_serial=0,
        checkpoint_generation=0,
        kernel_transition_parent_sha256=HLT17_GENESIS_DIGEST,
        kernel_transition_digest=HLT17_GENESIS_DIGEST,
    )
    return _with_tip(
        draft,
        parent=HLT17_GENESIS_DIGEST,
        generation=0,
        serial=0,
        checkpoint=0,
    )


def close_temporal_rejection(
    cursor: HLT17IMP1Cursor,
    *,
    predecessor_plan: TDG7Binary64SubdivisionPlan,
    predecessor_ledger: TDG11IMP1Ledger,
    updated_ledger: TDG11IMP1Ledger,
    rejection: TDG11IMP1TemporalRejection | Mapping[str, object],
    preparation_sha256: str,
    assessment_sha256: str,
) -> HLT17IMP1Cursor:
    """Bind one durable IMP1 nonadmission as RETRY_PENDING; keep source/CFL counts."""

    if not isinstance(cursor, HLT17IMP1Cursor):
        raise TypeError("cursor must be HLT17IMP1Cursor")
    record = (
        rejection
        if isinstance(rejection, TDG11IMP1TemporalRejection)
        else _rejection_from_mapping(rejection)
    )
    temporal = HLT17TemporalIdentity(
        rejection=record,
        predecessor_plan=predecessor_plan,
        successor_plan=temporal_successor_plan(predecessor_plan),
        predecessor_ledger=predecessor_ledger,
        preparation_sha256=preparation_sha256,
        assessment_sha256=assessment_sha256,
    )
    if ledger_sha256(temporal.updated_ledger) != ledger_sha256(updated_ledger):
        _fail("updated IMP1 ledger differs from the appended rejection ledger")
    if cursor.mode == RETRY_PENDING:
        if cursor.temporal is None:
            _fail("pending cursor omitted its temporal identity")
        previous_next = cursor.next_attempt_plan()
        if previous_next is None:
            _fail("pending cursor has no executable next plan to reject")
        if predecessor_plan != previous_next:
            _fail(
                "new temporal rejection did not execute the exact next overlay or "
                "IMP1 successor plan"
            )
    elif cursor.overlay is not None:
        previous_next = cursor.next_attempt_plan()
        if previous_next is None or predecessor_plan != previous_next:
            _fail("fresh overlay temporal rejection did not execute the overlay plan")
    draft = HLT17IMP1Cursor(
        member_key=cursor.member_key,
        identity_scope=cursor.identity_scope,
        **_implementation_fields(cursor),
        mode=RETRY_PENDING,
        event_target=cursor.event_target,
        accepted_time=cursor.accepted_time,
        accepted_step_index=cursor.accepted_step_index,
        accepted_transaction_serial=cursor.accepted_transaction_serial,
        descriptor_sha256=cursor.descriptor_sha256,
        physical_state_sha256=cursor.physical_state_sha256,
        coordinates_sha256=cursor.coordinates_sha256,
        synthetic_callable_binding_sha256=cursor.synthetic_callable_binding_sha256,
        monitor_sha256=cursor.monitor_sha256,
        causal_state_sha256=cursor.causal_state_sha256,
        tracer_sha256=cursor.tracer_sha256,
        guard_configuration_sha256=cursor.guard_configuration_sha256,
        ledger=updated_ledger,
        temporal=temporal,
        overlay=None,
        source_current=cursor.source_current,
        source_total=cursor.source_total,
        cfl_current=cursor.cfl_current,
        cfl_total=cursor.cfl_total,
        cursor_generation=cursor.cursor_generation,
        attempt_serial=cursor.attempt_serial,
        checkpoint_generation=cursor.checkpoint_generation,
        kernel_transition_parent_sha256=cursor.kernel_transition_parent_sha256,
        kernel_transition_digest=HLT17_GENESIS_DIGEST,
    )
    return _with_tip(
        draft,
        parent=cursor.kernel_transition_digest,
        generation=cursor.cursor_generation + 1,
        serial=cursor.attempt_serial + 1,
        checkpoint=cursor.checkpoint_generation,
    )


def close_source_cfl_overlay(
    cursor: HLT17IMP1Cursor,
    *,
    owner: str,
    attempted_plan: TDG7Binary64SubdivisionPlan,
    evidence: Mapping[str, object],
) -> HLT17IMP1Cursor:
    """Record a source/CFL overlay without dropping temporal identity or IMP1 history."""

    if not isinstance(cursor, HLT17IMP1Cursor):
        raise TypeError("cursor must be HLT17IMP1Cursor")
    if owner not in {SOURCE_OWNER, CFL_OWNER}:
        _fail("overlay owner must be source or CFL")
    expected_plan = cursor.next_attempt_plan()
    if cursor.mode == FRESH_READY and cursor.overlay is None:
        expected_plan = attempted_plan
    elif expected_plan is None or attempted_plan != expected_plan:
        _fail("source/CFL overlay did not attempt the exact next durable plan")
    if cursor.mode == RETRY_PENDING:
        if cursor.temporal is None:
            _fail("pending overlay omitted temporal identity")
        if attempted_plan == cursor.temporal.successor_plan:
            pass
        elif encode_tdg7_plan(attempted_plan) == encode_tdg7_plan(
            cursor.temporal.predecessor_plan
        ):
            _fail("source/CFL overlay must not overwrite the IMP1 successor plan")
    if owner == SOURCE_OWNER:
        source_current = cursor.source_current + 1
        source_total = cursor.source_total + 1
        cfl_current = cursor.cfl_current
        cfl_total = cursor.cfl_total
        cap = HLT17_SOURCE_RETRY_CAP
        owner_current = source_current
    else:
        source_current = cursor.source_current
        source_total = cursor.source_total
        cfl_current = cursor.cfl_current + 1
        cfl_total = cursor.cfl_total + 1
        cap = HLT17_CFL_RETRY_CAP
        owner_current = cfl_current
    exhausted = owner_current > cap
    next_plan: TDG7Binary64SubdivisionPlan | None = None
    if not exhausted:
        try:
            next_plan = strictly_smaller_tdg7_plan(attempted_plan, owner=owner)
        except HLT17OwnerRetryExhausted:
            exhausted = True
            next_plan = None
    overlay = HLT17RetryOverlay(
        owner=owner,
        attempted_plan=attempted_plan,
        next_plan=next_plan,
        evidence=dict(evidence),
        owner_current_count=owner_current,
        exhausted=exhausted,
    )
    draft = HLT17IMP1Cursor(
        member_key=cursor.member_key,
        identity_scope=cursor.identity_scope,
        **_implementation_fields(cursor),
        mode=cursor.mode,
        event_target=cursor.event_target,
        accepted_time=cursor.accepted_time,
        accepted_step_index=cursor.accepted_step_index,
        accepted_transaction_serial=cursor.accepted_transaction_serial,
        descriptor_sha256=cursor.descriptor_sha256,
        physical_state_sha256=cursor.physical_state_sha256,
        coordinates_sha256=cursor.coordinates_sha256,
        synthetic_callable_binding_sha256=cursor.synthetic_callable_binding_sha256,
        monitor_sha256=cursor.monitor_sha256,
        causal_state_sha256=cursor.causal_state_sha256,
        tracer_sha256=cursor.tracer_sha256,
        guard_configuration_sha256=cursor.guard_configuration_sha256,
        ledger=cursor.ledger,
        temporal=cursor.temporal,
        overlay=overlay,
        source_current=source_current,
        source_total=source_total,
        cfl_current=cfl_current,
        cfl_total=cfl_total,
        cursor_generation=cursor.cursor_generation,
        attempt_serial=cursor.attempt_serial,
        checkpoint_generation=cursor.checkpoint_generation,
        kernel_transition_parent_sha256=cursor.kernel_transition_parent_sha256,
        kernel_transition_digest=HLT17_GENESIS_DIGEST,
    )
    return _with_tip(
        draft,
        parent=cursor.kernel_transition_digest,
        generation=cursor.cursor_generation + 1,
        serial=cursor.attempt_serial + 1,
        checkpoint=cursor.checkpoint_generation,
    )


def close_fine_adoption(
    cursor: HLT17IMP1Cursor,
    *,
    committed_ledger: TDG11IMP1Ledger,
    descriptor_sha256: str,
    physical_state_sha256: str,
    coordinates_sha256: str,
    monitor_sha256: str,
    causal_state_sha256: str,
    tracer_sha256: str,
    guard_configuration_sha256: str,
) -> HLT17IMP1Cursor:
    """Reset per-macro counters only; retain cumulative counts, history, and debit.

    ``descriptor_sha256`` is the actual new accepted-encoding digest. The new
    cursor names that digest. It must not retain the predecessor descriptor
    after the accepted state advances, and it is not a hash of this cursor.
    """

    if not isinstance(cursor, HLT17IMP1Cursor):
        raise TypeError("cursor must be HLT17IMP1Cursor")
    if not isinstance(committed_ledger, TDG11IMP1Ledger):
        raise TypeError("committed_ledger must be TDG11IMP1Ledger")
    new_descriptor = _digest_text(descriptor_sha256, "descriptor_sha256")
    new_physical = _digest_text(physical_state_sha256, "physical_state_sha256")
    if new_descriptor == new_physical:
        _fail("descriptor and physical u/p/q identities are not interchangeable")
    if new_descriptor == cursor.descriptor_sha256:
        _fail("fine adoption must not retain the old descriptor after state advance")
    if committed_ledger.current_macro_step_temporal_retry_count != 0:
        _fail("fine adoption must reset the IMP1 per-macro temporal retry count")
    if (
        committed_ledger.cumulative_temporal_retry_count
        != cursor.ledger.cumulative_temporal_retry_count
    ):
        _fail("fine adoption must retain cumulative IMP1 temporal history")
    if (
        committed_ledger.serialized_temporal_rejections
        != cursor.ledger.serialized_temporal_rejections
    ):
        _fail("fine adoption must retain canonical IMP1 rejection history")
    if committed_ledger.inherited_snapshot != cursor.ledger.inherited_snapshot:
        _fail("fine adoption must retain the frozen inherited TDG6 snapshot")
    if committed_ledger.accepted_macro_step_count != (
        cursor.ledger.accepted_macro_step_count + 1
    ):
        _fail("fine adoption must advance exactly one accepted IMP1 macro step")
    draft = HLT17IMP1Cursor(
        member_key=cursor.member_key,
        identity_scope=cursor.identity_scope,
        **_implementation_fields(cursor),
        mode=FRESH_READY,
        event_target=cursor.event_target,
        accepted_time=committed_ledger.last_accepted_time,
        accepted_step_index=committed_ledger.current_step_index,
        accepted_transaction_serial=committed_ledger.current_transaction_serial,
        descriptor_sha256=new_descriptor,
        physical_state_sha256=new_physical,
        coordinates_sha256=coordinates_sha256,
        synthetic_callable_binding_sha256=cursor.synthetic_callable_binding_sha256,
        monitor_sha256=monitor_sha256,
        causal_state_sha256=causal_state_sha256,
        tracer_sha256=tracer_sha256,
        guard_configuration_sha256=guard_configuration_sha256,
        ledger=committed_ledger,
        temporal=None,
        overlay=None,
        source_current=0,
        source_total=cursor.source_total,
        cfl_current=0,
        cfl_total=cursor.cfl_total,
        cursor_generation=cursor.cursor_generation,
        attempt_serial=cursor.attempt_serial,
        checkpoint_generation=cursor.checkpoint_generation,
        kernel_transition_parent_sha256=cursor.kernel_transition_parent_sha256,
        kernel_transition_digest=HLT17_GENESIS_DIGEST,
    )
    return _with_tip(
        draft,
        parent=cursor.kernel_transition_digest,
        generation=cursor.cursor_generation + 1,
        serial=cursor.attempt_serial + 1,
        checkpoint=cursor.checkpoint_generation + 1,
    )


__all__ = [
    "CFL_EVIDENCE_KEYS",
    "CFL_OWNER",
    "CURSOR_KEYS",
    "FRESH_READY",
    "FUTURE_AUTHENTICATED_ORIGIN_SEAM",
    "FUTURE_CAMPAIGN_TRANSACTION_SEAM",
    "FUTURE_CHECKPOINT_JOURNAL_SEAM",
    "FUTURE_CONTENT_ADDRESSED_STORE_SEAM",
    "FUTURE_MEMBER_CODEC_SEAM",
    "FUTURE_SOURCE_MANIFEST_SEAM",
    "HLT17IMP1Cursor",
    "HLT17IMP1Error",
    "HLT17InvalidPremise",
    "HLT17MemberSpec",
    "HLT17OwnerRetryExhausted",
    "HLT17ProductionMember",
    "HLT17RetryOverlay",
    "HLT17TemporalIdentity",
    "HLT17_ADMISSION_RUNTIME_ID",
    "HLT17_CFL_RETRY_CAP",
    "HLT17_CURSOR_SCHEMA",
    "HLT17_CURSOR_SCHEMA_VERSION",
    "HLT17_IMPLEMENTATION_ID",
    "HLT17_LEGACY_IMP1_CURSOR_SCHEMA",
    "HLT17_MATHEMATICAL_OBJECT",
    "HLT17_MEMBER_KEYS",
    "HLT17_REFERENCE_WIRE_ARTIFACT_ID",
    "HLT17_REFERENCE_WIRE_EVALUATOR_ID",
    "HLT17_PRODUCTION_MEMBERS",
    "HLT17_SOURCE_RETRY_CAP",
    "IMPLEMENTATION_IDENTITY_KEYS",
    "HLT17_SYNTHETIC_MEMBERS",
    "HLT17_SYNTHETIC_MEMBER_KEYS",
    "IDENTITY_SCOPE_LIVE",
    "IDENTITY_SCOPE_SYNTHETIC",
    "PLAN_KEYS",
    "RETRY_PENDING",
    "SOURCE_EVIDENCE_KEYS",
    "SOURCE_OWNER",
    "TEMPORAL_OWNER",
    "close_fine_adoption",
    "close_source_cfl_overlay",
    "close_temporal_rejection",
    "cursor_sha256",
    "decode_tdg7_plan",
    "encode_hlt17_cursor",
    "encode_tdg7_plan",
    "hlt17_member_spec",
    "imp1_enclosure_limits_for_member",
    "ledger_sha256",
    "plan_imp1_lattice",
    "production_member",
    "restore_hlt17_cursor",
    "revalidate_hlt17_cursor",
    "seed_hlt17_cursor",
    "strictly_smaller_tdg7_plan",
    "temporal_successor_plan",
]
