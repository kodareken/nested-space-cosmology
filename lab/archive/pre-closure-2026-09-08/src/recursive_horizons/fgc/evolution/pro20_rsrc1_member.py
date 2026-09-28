"""One-member physical reconstruction for prospective PRO20 RSRC1.

The closed PRO20 runner retained six live members in one process.  RSRC1
instead constructs exactly the scheduled member in an isolated child.  This
module performs no store read or write, publishes no endpoint, and grants no
execution authority.  A later exact freeze must bind every dependency before
the physical path may be measured.
"""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
from typing import Mapping

from .hlt16_member_codec import _digest, _template_identity, restore_member
from .hlt17_member_checkpoint import HLT17InMemoryCheckpoint
from .hlt17_member_codec import HLT17PredecessorIdentity, ORIGIN_PROTOCOL
from .hlt17_runtime_member import (
    HLT17ExternalReference,
    HLT17MemberIdentity,
    HLT17RuntimeMember,
    IDENTITY_SCOPE_LIVE,
)
from .numerical_engine import array_content_sha256
from .pro20_origin import (
    MEMBER_KEYS,
    Pro20OriginCapture,
    revalidate_pro20_origin_capture,
)
from .pro20_source_binding import Pro20GR0SourceBinding
from .proto19_gr0_static_factory import (
    AMPLITUDE,
    _PERSISTED_TEMPLATE_SHA256,
    _build_member,
    _decode_inputs,
    _mapping,
    _validated_sources,
)
from .tdg11_imp1_ledger import seed_imp1_ledger


ARTIFACT_ID = "FGC-1-PRO20-EV1-RSRC1"
PROTOCOL = "FGC-2-SF1-PROTO19"
CAMPAIGN_ID = "FGC-2-SF1-PRO20-RSRC1-EVENT1"
EVENT_TARGET = 3.0 / 2.0
PINNED_PRIVATE_DEPENDENCIES = (
    "proto19_gr0_static_factory._PERSISTED_TEMPLATE_SHA256",
    "proto19_gr0_static_factory._build_member",
    "proto19_gr0_static_factory._decode_inputs",
    "proto19_gr0_static_factory._mapping",
    "proto19_gr0_static_factory._validated_sources",
)


class Pro20RSRC1MemberError(ValueError):
    """The one-member reconstruction differs from the frozen physical owner."""


def _require_member_key(member_key: object) -> str:
    if type(member_key) is not str or member_key not in MEMBER_KEYS:
        raise Pro20RSRC1MemberError("RSRC1 member key is outside the live cohort")
    return member_key


def _require_sha256(value: object, label: str) -> str:
    if type(value) is not str or len(value) != 64 or value.lower() != value:
        raise Pro20RSRC1MemberError(f"{label} is not lowercase SHA-256")
    try:
        int(value, 16)
    except ValueError as error:
        raise Pro20RSRC1MemberError(f"{label} is not lowercase SHA-256") from error
    return value


def _origin_receipt(capture: Pro20OriginCapture) -> str:
    receipts = [
        leaf
        for leaf in capture.original_leaves
        if leaf.path.endswith("/receipts/generation-zero.json")
    ]
    if len(receipts) != 1:
        raise Pro20RSRC1MemberError("captured HLT15 origin receipt differs")
    return _require_sha256(receipts[0].content_id, "origin receipt")


def _captured_member(capture: Pro20OriginCapture, member_key: str):
    matches = [member for member in capture.members if member.member_key == member_key]
    if len(matches) != 1:
        raise Pro20RSRC1MemberError(f"captured origin member differs: {member_key}")
    return matches[0]


def build_rsrc1_static_shell(
    member_key: str,
    *,
    static_input_bytes: Mapping[str, bytes],
):
    """Construct exactly one authenticated amplitude-three GR-0 shell."""

    key = _require_member_key(member_key)
    try:
        cal9, rsp2, hlt10, rsp2_freeze = _decode_inputs(None, static_input_bytes)
        cal9, cal9_inputs, rsp2_input = _validated_sources(
            cal9, rsp2, hlt10, rsp2_freeze
        )
        cal9_physical = _mapping(cal9.get("physical_inputs"), "CAL9 physical inputs")
        cal9_numerics = _mapping(cal9.get("numerics"), "CAL9 numerics")
        cal9_thresholds = _mapping(
            cal9.get("universal_thresholds"), "CAL9 universal thresholds"
        )
        method_label, point_text = key.split("-", 1)
        point_count = int(point_text)
        if key == "SSPRK3-16385":
            physical = _mapping(rsp2.get("physical_inputs"), "RSP2 physical inputs")
            method = _mapping(rsp2.get("method"), "RSP2 method")
            numerics = _mapping(rsp2.get("numerics"), "RSP2 numerics")
            thresholds = _mapping(
                rsp2.get("universal_thresholds"), "RSP2 universal thresholds"
            )
            frozen_input = rsp2_input
        else:
            physical = cal9_physical
            method = _mapping(cal9.get(
                "primary_method" if method_label == "RK4" else "comparator_method"
            ), f"CAL9 {method_label} method")
            numerics = cal9_numerics
            thresholds = cal9_thresholds
            frozen_input = cal9_inputs[(method_label, point_count)]
        member = _build_member(
            physical=physical,
            method=method,
            numerics=numerics,
            thresholds=thresholds,
            frozen_input=frozen_input,
            point_count=point_count,
        )
    except Pro20RSRC1MemberError:
        raise
    except Exception as error:
        raise Pro20RSRC1MemberError(
            f"RSRC1 static construction failed: {key}"
        ) from error
    if _digest(_template_identity(member)) != _PERSISTED_TEMPLATE_SHA256[key]:
        raise Pro20RSRC1MemberError(
            f"RSRC1 static template differs from persisted HLT16: {key}"
        )
    return member


def _predecessor(captured, runtime: HLT17RuntimeMember) -> HLT17PredecessorIdentity:
    identity = {
        "protocol_artifact_id": ORIGIN_PROTOCOL,
        "campaign_id": CAMPAIGN_ID,
        "member_key": runtime.key,
        "method": runtime.method_label,
        "point_count": runtime.point_count,
        "accepted_boundary_time_hex": runtime.time.hex(),
        "accepted_state_sha256": runtime.physical_state_sha256(),
        "previous_step_index": runtime.step_index,
        "previous_transaction_serial": runtime.transaction_serial,
        "accepted_generation": 0,
    }
    return HLT17PredecessorIdentity(
        descriptor_sha256=captured.descriptor_content_id,
        cursor_identity=identity,
        cursor_sha256=captured.full_cursor_sha256,
    )


@dataclass(frozen=True, slots=True)
class Pro20RSRC1RebasedMember:
    member_key: str
    captured_origin_sha256: str
    generation1_bridge_content_id: str
    source_closure_sha256: str
    source_configuration_sha256: str
    source_binding: Pro20GR0SourceBinding
    member: HLT17RuntimeMember
    predecessor: HLT17PredecessorIdentity
    checkpoint: HLT17InMemoryCheckpoint

    def __post_init__(self) -> None:
        _require_member_key(self.member_key)
        _require_sha256(self.captured_origin_sha256, "captured origin")
        _require_sha256(self.generation1_bridge_content_id, "generation-one bridge")
        _require_sha256(self.source_closure_sha256, "source closure")
        _require_sha256(self.source_configuration_sha256, "source configuration")
        if self.member.key != self.member_key or self.checkpoint.member.key != self.member_key:
            raise Pro20RSRC1MemberError("RSRC1 rebased member identity differs")
        if self.source_binding.configuration_sha256 != self.source_configuration_sha256:
            raise Pro20RSRC1MemberError("RSRC1 source configuration differs")
        if self.member.source_operator_closure.sha256 != self.source_closure_sha256:
            raise Pro20RSRC1MemberError("RSRC1 source closure reference differs")
        self.source_binding.validate()
        self.member.__post_init__()
        self.checkpoint.agree()

    @property
    def bundle(self):
        return self.checkpoint.generation_bundle()

    @property
    def construction_sha256(self) -> str:
        payload = (
            f"{CAMPAIGN_ID}\n{self.member_key}\n{self.captured_origin_sha256}\n"
            f"{self.generation1_bridge_content_id}\n{self.source_closure_sha256}\n"
            f"{self.source_configuration_sha256}\n{self.bundle.descriptor_sha256}\n"
            f"{self.bundle.payload_sha256}\n{self.bundle.cursor_sha256}\n"
        ).encode("ascii")
        return sha256(payload).hexdigest()


def build_rsrc1_member(
    origin: Pro20OriginCapture,
    *,
    member_key: str,
    static_input_bytes: Mapping[str, bytes],
    source_closure_sha256: str,
) -> Pro20RSRC1RebasedMember:
    """Build and HLT15-restore one scheduled live member; write nothing."""

    key = _require_member_key(member_key)
    revalidate_pro20_origin_capture(origin)
    source_closure = _require_sha256(source_closure_sha256, "source closure")
    captured = _captured_member(origin, key)
    shell = build_rsrc1_static_shell(key, static_input_bytes=static_input_bytes)
    try:
        restore_member(shell, origin.fresh_hlt16_snapshot(key))
    except Exception as error:
        raise Pro20RSRC1MemberError(f"RSRC1 HLT15 restore failed: {key}") from error
    source_binding = Pro20GR0SourceBinding(shell.operator, shell.projector)
    state_sha256 = array_content_sha256(shell.state.u, shell.state.p, shell.state.q)
    receipt = _origin_receipt(origin)
    ledger = seed_imp1_ledger(
        shell.temporal_ledger,
        method=shell.integrator_id,
        state_sha256=state_sha256,
        step_index=shell.step_index,
        transaction_serial=shell.transaction_serial,
        origin_receipt_sha256=receipt,
    )
    identity = HLT17MemberIdentity(
        scope=IDENTITY_SCOPE_LIVE,
        method_label=shell.method_label,
        point_count=shell.point_count,
        integrator_id=shell.integrator_id,
        spatial_order=shell.spatial_order,
        campaign_id=CAMPAIGN_ID,
        amplitude=AMPLITUDE,
    )
    member = HLT17RuntimeMember(
        identity=identity,
        initial=shell.initial,
        operator=source_binding.rhs,
        projector=source_binding.project,
        transaction=shell.transaction,
        tracers=shell.tracers,
        state=shell.state,
        temporal_ledger=ledger,
        source_operator_closure=HLT17ExternalReference(
            kind="source_operator_closure", sha256=source_closure
        ),
        origin_reference=HLT17ExternalReference(
            kind="origin_receipt", sha256=receipt
        ),
        time=shell.time,
        step_index=shell.step_index,
        transaction_serial=shell.transaction_serial,
        CFL_retry_count=shell.CFL_retry_count,
        source_retry_count=shell.source_retry_count,
        executed_plan=None,
    )
    predecessor = _predecessor(captured, member)
    checkpoint = HLT17InMemoryCheckpoint.seed(
        member,
        event_target=EVENT_TARGET,
        origin_predecessor=predecessor,
    )
    return Pro20RSRC1RebasedMember(
        member_key=key,
        captured_origin_sha256=origin.capture_sha256,
        generation1_bridge_content_id=origin.generation1_bridge_content_id,
        source_closure_sha256=source_closure,
        source_configuration_sha256=source_binding.configuration_sha256,
        source_binding=source_binding,
        member=member,
        predecessor=predecessor,
        checkpoint=checkpoint,
    )


__all__ = [
    "ARTIFACT_ID",
    "CAMPAIGN_ID",
    "EVENT_TARGET",
    "PINNED_PRIVATE_DEPENDENCIES",
    "PROTOCOL",
    "Pro20RSRC1MemberError",
    "Pro20RSRC1RebasedMember",
    "build_rsrc1_member",
    "build_rsrc1_static_shell",
]
