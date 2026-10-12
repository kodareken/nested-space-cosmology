"""Physical-source construction and HLT17 rebase for prospective PRO20.

The public builder deliberately performs the preserved static GR-0
construction, whose historical factory evaluates the physical source at time
zero. It must therefore be called only by a separately committed source/
environment/resource qualification or later authority. It performs no file
write, store publication, campaign step, or endpoint adoption.
"""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
from types import MappingProxyType
from typing import Mapping

from recursive_horizons.evidence_io import canonical_json_bytes

from .hlt16_member_codec import restore_member
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
from .proto19_gr0_static_factory import build_static_gr0_shells
from .tdg11_imp1_ledger import seed_imp1_ledger


PROTOCOL = "FGC-2-SF1-PROTO19"
CAMPAIGN_ID = "FGC-2-SF1-PRO20-EVENT1"
EVENT_TARGET = 3.0 / 2.0


class Pro20SourceFactoryError(ValueError):
    """The captured origin, fresh shell, or HLT17 rebase differs."""


def _sha(value: object, label: str) -> str:
    if type(value) is not str or len(value) != 64 or value.lower() != value:
        raise Pro20SourceFactoryError(f"{label} must be a lowercase SHA-256 digest")
    try:
        int(value, 16)
    except ValueError as error:
        raise Pro20SourceFactoryError(f"{label} must be a lowercase SHA-256 digest") from error
    return value


def _origin_receipt(capture: Pro20OriginCapture) -> str:
    receipts = [
        leaf for leaf in capture.original_leaves
        if leaf.path.endswith("/receipts/generation-zero.json")
    ]
    if len(receipts) != 1:
        raise Pro20SourceFactoryError("captured HLT15 origin receipt differs")
    return _sha(receipts[0].content_id, "origin receipt")


def _captured_member(capture: Pro20OriginCapture, member_key: str):
    values = [member for member in capture.members if member.member_key == member_key]
    if len(values) != 1:
        raise Pro20SourceFactoryError(f"captured origin member differs: {member_key}")
    return values[0]


def _predecessor(
    captured,
    runtime: HLT17RuntimeMember,
) -> HLT17PredecessorIdentity:
    """Explicit PROTO17-to-PROTO19 rebase projection, not old cursor metadata."""

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
class Pro20RebasedMember:
    member_key: str
    captured_origin_sha256: str
    source_closure_sha256: str
    source_configuration_sha256: str
    source_binding: Pro20GR0SourceBinding
    member: HLT17RuntimeMember
    predecessor: HLT17PredecessorIdentity
    checkpoint: HLT17InMemoryCheckpoint

    def __post_init__(self) -> None:
        _sha(self.captured_origin_sha256, "captured origin")
        _sha(self.source_closure_sha256, "source closure")
        _sha(self.source_configuration_sha256, "source configuration")
        if self.member_key != self.member.key or self.member_key != self.checkpoint.member.key:
            raise Pro20SourceFactoryError("rebased member identity differs")
        if self.source_configuration_sha256 != self.source_binding.configuration_sha256:
            raise Pro20SourceFactoryError("rebased source configuration differs")
        if self.source_closure_sha256 != self.member.source_operator_closure.sha256:
            raise Pro20SourceFactoryError("rebased source closure reference differs")
        if self.predecessor.descriptor_sha256 == self.member.origin_reference.sha256:
            raise Pro20SourceFactoryError(
                "origin descriptor address collapsed into the receipt reference"
            )
        self.source_binding.validate()
        self.member.__post_init__()
        self.checkpoint.agree()

    @property
    def bundle(self):
        return self.checkpoint.generation_bundle()


@dataclass(frozen=True, slots=True)
class Pro20RuntimeOrigin:
    captured_origin_sha256: str
    source_closure_sha256: str
    members: Mapping[str, Pro20RebasedMember]

    def __post_init__(self) -> None:
        _sha(self.captured_origin_sha256, "captured origin")
        _sha(self.source_closure_sha256, "source closure")
        members = MappingProxyType(dict(self.members))
        object.__setattr__(self, "members", members)
        if tuple(members) != MEMBER_KEYS:
            raise Pro20SourceFactoryError("runtime origin member order differs")
        for key, value in members.items():
            if type(value) is not Pro20RebasedMember or value.member_key != key:
                raise Pro20SourceFactoryError("runtime origin member record differs")
            if value.captured_origin_sha256 != self.captured_origin_sha256:
                raise Pro20SourceFactoryError("runtime origin capture identity differs")
            if value.source_closure_sha256 != self.source_closure_sha256:
                raise Pro20SourceFactoryError("runtime origin source identity differs")

    @property
    def construction_sha256(self) -> str:
        return sha256(canonical_json_bytes({
            "campaign_id": CAMPAIGN_ID,
            "captured_origin_sha256": self.captured_origin_sha256,
            "member_bundles": {
                key: {
                    "configuration_sha256": value.source_configuration_sha256,
                    "descriptor_sha256": value.bundle.descriptor_sha256,
                    "physical_state_sha256": value.bundle.physical_state_sha256,
                    "cursor_sha256": value.bundle.cursor_sha256,
                }
                for key, value in self.members.items()
            },
            "protocol": PROTOCOL,
            "source_closure_sha256": self.source_closure_sha256,
        })).hexdigest()


def build_pro20_runtime_origin(
    origin: Pro20OriginCapture,
    *,
    static_input_bytes: Mapping[str, bytes],
    source_closure_sha256: str,
) -> Pro20RuntimeOrigin:
    """Build six fresh physical shells, restore HLT15, and seed HLT17.

    Construction invokes the existing physical source only inside
    ``build_static_gr0_shells``. Returned external references remain
    unauthenticated locally; later qualification/authority owns their proof.
    """

    revalidate_pro20_origin_capture(origin)
    source_closure = _sha(source_closure_sha256, "source closure")
    try:
        shells = build_static_gr0_shells(static_input_bytes=static_input_bytes)
    except Exception as error:
        raise Pro20SourceFactoryError("static GR-0 physical construction failed") from error
    if tuple(shells) != MEMBER_KEYS:
        raise Pro20SourceFactoryError("static GR-0 shell order differs")
    receipt = _origin_receipt(origin)
    rebased: dict[str, Pro20RebasedMember] = {}
    for key in MEMBER_KEYS:
        captured = _captured_member(origin, key)
        shell = shells[key]
        try:
            restore_member(shell, origin.fresh_hlt16_snapshot(key))
        except Exception as error:
            raise Pro20SourceFactoryError(f"{key} HLT15 restore failed") from error
        source_binding = Pro20GR0SourceBinding(shell.operator, shell.projector)
        state_sha256 = array_content_sha256(shell.state.u, shell.state.p, shell.state.q)
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
            amplitude="3",
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
        rebased[key] = Pro20RebasedMember(
            member_key=key,
            captured_origin_sha256=origin.capture_sha256,
            source_closure_sha256=source_closure,
            source_configuration_sha256=source_binding.configuration_sha256,
            source_binding=source_binding,
            member=member,
            predecessor=predecessor,
            checkpoint=checkpoint,
        )
    return Pro20RuntimeOrigin(
        captured_origin_sha256=origin.capture_sha256,
        source_closure_sha256=source_closure,
        members=rebased,
    )


__all__ = [
    "CAMPAIGN_ID",
    "EVENT_TARGET",
    "PROTOCOL",
    "Pro20RebasedMember",
    "Pro20RuntimeOrigin",
    "Pro20SourceFactoryError",
    "build_pro20_runtime_origin",
]
