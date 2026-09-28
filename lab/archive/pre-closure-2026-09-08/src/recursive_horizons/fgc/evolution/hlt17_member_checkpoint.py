"""Join HLT17 cursor/bridge and member/codec into one no-I/O checkpoint layer.

This is still a synthetic in-memory integration. It does not write a store,
read raw origin data, emit a compact certificate, or authorize a campaign.
A later exclusive writer may publish the returned generation-bundle parts
atomically; this module only returns those parts.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from hashlib import sha256
from typing import Any, Callable, Mapping

import numpy as np

from recursive_horizons.evidence_io import CanonicalJSONError, canonical_json_bytes, load_canonical_json

from .hlt17_admission_runtime import (
    FUTURE_STORE_IMPLEMENTATION_SEAM,
    TDG11C1R1CommittedRuntime,
    TDG11C1R1PreparedRuntime,
    hlt17_implementation_identity,
    require_hlt17_committed,
    require_hlt17_prepared,
)
from .hlt17_imp1_bridge import (
    FUTURE_INTEGRATION_SEAMS,
    HLT17BridgeResult,
    HLT17LiveBoundary,
    authenticate_live_boundary,
    commit_hlt17_kernel,
    finalize_hlt17_fine,
    live_boundary_identities,
    prepare_hlt17_fresh,
    prepare_hlt17_next,
    require_hlt17_admission,
)
from .hlt17_imp1_cursor import (
    FUTURE_AUTHENTICATED_ORIGIN_SEAM,
    FUTURE_CAMPAIGN_TRANSACTION_SEAM,
    FUTURE_CHECKPOINT_JOURNAL_SEAM,
    FUTURE_CONTENT_ADDRESSED_STORE_SEAM,
    FUTURE_SOURCE_MANIFEST_SEAM,
    HLT17IMP1Cursor,
    HLT17IMP1Error,
    HLT17InvalidPremise,
    cursor_sha256,
    encode_hlt17_cursor,
    imp1_enclosure_limits_for_member,
    ledger_sha256,
    restore_hlt17_cursor,
    revalidate_hlt17_cursor,
    seed_hlt17_cursor,
)
from .hlt17_member_codec import (
    ORIGIN_PROTOCOL,
    RUNTIME_PROTOCOL,
    SNAPSHOT_ROLE_ACCEPTED,
    HLT17MemberCodecError,
    HLT17MemberEncoding,
    HLT17PredecessorIdentity,
    encode_member,
    restore_member,
)
from .hlt17_runtime_member import (
    IDENTITY_SCOPE_LIVE,
    IDENTITY_SCOPE_SYNTHETIC,
    HLT17RollbackError,
    HLT17RuntimeMember,
    HLT17RuntimeMemberError,
)
from .tdg11_imp1_ledger import TDG11IMP1Ledger
from .proto5_runtime import GR0RuntimeStageTransaction

from .tdg6_temporal_admission_runtime import _clone_tracers
from .tdg7_binary64_subdivision_lattice import TDG7Binary64SubdivisionPlan


CHECKPOINT_PROTOCOL = "FGC-2-SF1-PROTO19"
PARKED_AT_EVENT_TARGET = "parked_at_event_target"
SEEDED = "seeded"
ACCEPTED_FINE = "accepted_fine"
FUTURE_ATOMIC_GENERATION_BUNDLE_SEAM = (
    "an exclusive writer must validate this in-memory checkpoint transition "
    "before publishing descriptor+payload+cursor as one atomic generation "
    "bundle; this layer returns the parts and does not write a store, bind a "
    "journal tip, or authorize a campaign"
)
PINNED_PRIVATE_DEPENDENCIES = (
    "recursive_horizons.fgc.evolution.tdg6_temporal_admission_runtime._clone_tracers",
)
IDENTITY_LABELS = (
    "physical_u_p_q",
    "accepted_descriptor",
    "full_kernel_cursor",
    "kernel_transition_digest",
    "future_store_journal",
)
SYNTHETIC_ORIGIN_CURSOR_PLACEHOLDER_DOMAIN = (
    b"hlt17-synthetic-origin-cursor-placeholder-v1\n"
)
LIVE_BOUNDARY_IDENTITY_NAMES = (
    "physical_state_sha256",
    "coordinates_sha256",
    "synthetic_callable_binding_sha256",
    "monitor_sha256",
    "causal_state_sha256",
    "tracer_sha256",
    "guard_configuration_sha256",
)


class HLT17CheckpointError(ValueError):
    """Malformed joint HLT17 checkpoint, restore, or adoption data."""


class HLT17CheckpointRollbackError(RuntimeError):
    """A late joint failure could not restore the complete prior boundary."""


def _fail(message: str) -> None:
    raise HLT17CheckpointError(message)


def _same_bits(left: float, right: float) -> bool:
    return float(left).hex() == float(right).hex()


def _clone_stage_transaction(
    transaction: GR0RuntimeStageTransaction,
) -> GR0RuntimeStageTransaction:
    if not isinstance(transaction, GR0RuntimeStageTransaction):
        raise TypeError("transaction must be GR0RuntimeStageTransaction")
    clone = GR0RuntimeStageTransaction(
        thresholds=transaction.thresholds,
        causal_state=transaction.causal_state,
        boundary_geometry=transaction.boundary_geometry,
        grid_spacing=transaction.grid_spacing,
        cfl_maximum=transaction.cfl_maximum,
        hat_normal_factor=transaction.hat_normal_factor,
    )
    clone.state = transaction.state
    return clone


def _isolated_member(member: HLT17RuntimeMember) -> HLT17RuntimeMember:
    member.__post_init__()
    return HLT17RuntimeMember(
        identity=member.identity,
        initial=member.initial,
        operator=member.operator,
        projector=member.projector,
        transaction=_clone_stage_transaction(member.transaction),
        tracers=_clone_tracers(member.tracers),
        state=member.state,
        temporal_ledger=member.temporal_ledger,
        source_operator_closure=member.source_operator_closure,
        origin_reference=member.origin_reference,
        time=member.time,
        step_index=member.step_index,
        transaction_serial=member.transaction_serial,
        CFL_retry_count=member.CFL_retry_count,
        source_retry_count=member.source_retry_count,
        executed_plan=member.executed_plan,
    )


def _origin_cursor_identity(member: HLT17RuntimeMember) -> dict[str, object]:
    return {
        "protocol_artifact_id": ORIGIN_PROTOCOL,
        "campaign_id": member.identity.campaign_id,
        "member_key": member.key,
        "method": member.identity.method_label,
        "point_count": member.identity.point_count,
        "accepted_boundary_time_hex": member.time.hex(),
        "accepted_state_sha256": member.physical_state_sha256(),
        "previous_step_index": member.step_index,
        "previous_transaction_serial": member.transaction_serial,
        "accepted_generation": 0,
    }


def synthetic_origin_predecessor_identity(
    member: HLT17RuntimeMember,
) -> HLT17PredecessorIdentity:
    """Synthetic-only PROTO17 placeholder. Not a real origin cursor digest."""

    if not isinstance(member, HLT17RuntimeMember):
        raise TypeError("member must be HLT17RuntimeMember")
    member.__post_init__()
    if member.identity.scope != IDENTITY_SCOPE_SYNTHETIC:
        _fail("synthetic origin cursor placeholder cannot be used for a live member")
    identity = _origin_cursor_identity(member)
    full = sha256(
        SYNTHETIC_ORIGIN_CURSOR_PLACEHOLDER_DOMAIN + canonical_json_bytes(identity)
    ).hexdigest()
    return HLT17PredecessorIdentity(
        descriptor_sha256=member.origin_reference.sha256,
        cursor_identity=identity,
        cursor_sha256=full,
    )


def origin_predecessor_identity(
    member: HLT17RuntimeMember,
    *,
    origin_cursor_sha256: str | None = None,
) -> HLT17PredecessorIdentity:
    """PROTO17 origin reference. Matching the digest is not origin authentication.

    Live members require an explicitly supplied origin cursor digest. Synthetic
    members may use a clearly labeled placeholder; that placeholder is not a
    manufactured real PROTO17 cursor.
    """

    if not isinstance(member, HLT17RuntimeMember):
        raise TypeError("member must be HLT17RuntimeMember")
    member.__post_init__()
    if member.identity.scope == IDENTITY_SCOPE_LIVE:
        if origin_cursor_sha256 is None:
            _fail("live members require an explicitly supplied origin cursor reference")
        identity = _origin_cursor_identity(member)
        return HLT17PredecessorIdentity(
            descriptor_sha256=member.origin_reference.sha256,
            cursor_identity=identity,
            cursor_sha256=origin_cursor_sha256,
        )
    if origin_cursor_sha256 is not None:
        identity = _origin_cursor_identity(member)
        return HLT17PredecessorIdentity(
            descriptor_sha256=member.origin_reference.sha256,
            cursor_identity=identity,
            cursor_sha256=origin_cursor_sha256,
        )
    return synthetic_origin_predecessor_identity(member)


def predecessor_from_accepted_cursor(
    cursor: HLT17IMP1Cursor,
    *,
    campaign_id: str,
    accepted_generation: int,
    descriptor_sha256: str,
) -> HLT17PredecessorIdentity:
    """Name the old full kernel cursor from a new accepted descriptor.

    ``identity_sha256`` hashes the small projection. ``cursor_sha256`` is the
    actual full-cursor digest. Neither matching is authentication.
    """

    revalidate_hlt17_cursor(cursor)
    if type(campaign_id) is not str or not campaign_id:
        _fail("campaign_id must be nonempty text")
    if type(accepted_generation) is not int or accepted_generation < 0:
        raise TypeError("accepted_generation must be a nonnegative built-in integer")
    identity = {
        "protocol_artifact_id": RUNTIME_PROTOCOL,
        "campaign_id": campaign_id,
        "member_key": cursor.member_key,
        "method": cursor.member.method_label,
        "point_count": cursor.member.point_count,
        "accepted_boundary_time_hex": cursor.accepted_time.hex(),
        "accepted_state_sha256": cursor.physical_state_sha256,
        "previous_step_index": cursor.accepted_step_index,
        "previous_transaction_serial": cursor.accepted_transaction_serial,
        "accepted_generation": accepted_generation,
    }
    full = cursor_sha256(cursor)
    projected = sha256(canonical_json_bytes(identity)).hexdigest()
    if projected == full:
        _fail("projected cursor-identity digest collided with the full-cursor digest")
    return HLT17PredecessorIdentity(
        descriptor_sha256=descriptor_sha256,
        cursor_identity=identity,
        cursor_sha256=full,
    )


def _descriptor_mapping(encoding: HLT17MemberEncoding) -> dict[str, Any]:
    if not isinstance(encoding, HLT17MemberEncoding):
        raise TypeError("encoding must be HLT17MemberEncoding")
    try:
        metadata = load_canonical_json(encoding.descriptor)
    except (CanonicalJSONError, TypeError, ValueError) as error:
        raise HLT17CheckpointError("accepted descriptor is not canonical JSON") from error
    if type(metadata) is not dict:
        _fail("accepted descriptor is not a JSON object")
    return metadata


def _require_no_hash_cycle(
    *,
    old_cursor: HLT17IMP1Cursor,
    new_cursor: HLT17IMP1Cursor,
    encoding: HLT17MemberEncoding,
) -> None:
    metadata = _descriptor_mapping(encoding)
    if metadata.get("snapshot_role") != SNAPSHOT_ROLE_ACCEPTED:
        _fail("fine adoption must produce an accepted physical-state descriptor")
    if any(
        key in metadata
        for key in ("successor_cursor", "successor_cursor_sha256", "cursor_chain_sha256")
    ):
        _fail("new descriptor embeds a successor cursor")
    old_full = cursor_sha256(old_cursor)
    new_full = cursor_sha256(new_cursor)
    if metadata["predecessor_cursor_sha256"] != old_full:
        _fail("new descriptor must name the actual old full kernel cursor")
    if metadata["predecessor_cursor_identity_sha256"] == metadata["predecessor_cursor_sha256"]:
        _fail("projected cursor-identity digest is not the full-cursor digest")
    if metadata["predecessor_descriptor_sha256"] != old_cursor.descriptor_sha256:
        _fail("new descriptor must name the old accepted descriptor")
    if new_cursor.descriptor_sha256 != encoding.descriptor_sha256:
        _fail("new cursor must name the actual new descriptor")
    if new_cursor.descriptor_sha256 == old_cursor.descriptor_sha256:
        _fail("new cursor retained the old descriptor after state advance")
    if new_full == encoding.descriptor_sha256 or old_full == encoding.descriptor_sha256:
        _fail("descriptor digest collided with a full-cursor digest")
    if new_full in metadata.values() or encoding.descriptor_sha256 == old_full:
        _fail("descriptor/cursor digests formed a circular definition")
    text = encoding.descriptor.decode("ascii")
    if new_full in text:
        _fail("new descriptor includes the new full-cursor digest")
    if new_cursor.physical_state_sha256 != encoding.physical_state_sha256:
        _fail("new cursor physical hash differs from the new encoding")
    if new_cursor.physical_state_sha256 == new_cursor.descriptor_sha256:
        _fail("descriptor and physical u/p/q identities are not interchangeable")


def _as_encoding(encoding: object) -> HLT17MemberEncoding:
    if isinstance(encoding, HLT17MemberEncoding):
        encoding.__post_init__()
        return encoding
    raise TypeError("encoding must be HLT17MemberEncoding")


def _full_cursor_bytes(cursor: HLT17IMP1Cursor) -> bytes:
    return canonical_json_bytes(encode_hlt17_cursor(revalidate_hlt17_cursor(cursor)))


def _restore_cursor_bytes(raw: object) -> tuple[bytes, HLT17IMP1Cursor]:
    if type(raw) is not bytes:
        raise TypeError("cursor_bytes must be bytes")
    try:
        mapping = load_canonical_json(raw)
    except (CanonicalJSONError, TypeError, ValueError) as error:
        raise HLT17CheckpointError("full-cursor bytes are not canonical JSON") from error
    try:
        cursor = restore_hlt17_cursor(mapping)
    except (HLT17IMP1Error, TypeError, ValueError) as error:
        raise HLT17CheckpointError("full-cursor bytes are not a closed HLT17 cursor") from error
    canonical = _full_cursor_bytes(cursor)
    if canonical != raw:
        _fail("full-cursor bytes are not the canonical cursor encoding")
    return canonical, cursor


def _require_revision_relationship(
    *, accepted_generation: int, cursor: HLT17IMP1Cursor
) -> None:
    if type(accepted_generation) is not int or accepted_generation < 0:
        _fail("accepted_generation must be a nonnegative built-in integer")
    if cursor.checkpoint_generation != accepted_generation:
        _fail("kernel checkpoint generation differs from accepted physical generation")
    if cursor.cursor_generation < accepted_generation:
        _fail("cursor generation cannot precede accepted physical generation")
    if cursor.attempt_serial < accepted_generation:
        _fail("attempt serial cannot precede accepted physical generation")
    if cursor.ledger.accepted_macro_step_count != accepted_generation:
        _fail("IMP1 accepted macro-step count differs from accepted generation")


def _require_legal_retry_suffix(
    accepted: TDG11IMP1Ledger, live: TDG11IMP1Ledger
) -> None:
    """Live IMP1 ledger must be the accepted ledger or its current-step suffix."""

    if not isinstance(accepted, TDG11IMP1Ledger) or not isinstance(live, TDG11IMP1Ledger):
        raise TypeError("retry suffix comparison requires TDG11IMP1Ledger")
    if accepted.current_macro_step_temporal_retry_count != 0:
        _fail("accepted encoding cannot carry an active retry ledger")
    if ledger_sha256(accepted) == ledger_sha256(live):
        if live.current_macro_step_temporal_retry_count != 0:
            _fail("identical accepted ledger cannot carry an active temporal retry")
        return
    if live.inherited_snapshot != accepted.inherited_snapshot:
        _fail("retry ledger cannot rewrite the inherited TDG6 sibling")
    if live.origin_receipt_sha256 != accepted.origin_receipt_sha256:
        _fail("retry ledger origin receipt differs")
    if live.origin_state_sha256 != accepted.origin_state_sha256:
        _fail("retry ledger origin state differs")
    if live.origin_step_index != accepted.origin_step_index:
        _fail("retry ledger origin step differs")
    if live.origin_transaction_serial != accepted.origin_transaction_serial:
        _fail("retry ledger origin serial differs")
    if live.method != accepted.method:
        _fail("retry ledger method differs")
    if not _same_bits(live.last_accepted_time, accepted.last_accepted_time):
        _fail("retry ledger cannot advance accepted time")
    if live.current_state_sha256 != accepted.current_state_sha256:
        _fail("retry ledger cannot replace the accepted physical state")
    if live.current_step_index != accepted.current_step_index:
        _fail("retry ledger cannot advance the step index")
    if live.current_transaction_serial != accepted.current_transaction_serial:
        _fail("retry ledger cannot advance the transaction serial")
    if live.accepted_macro_step_count != accepted.accepted_macro_step_count:
        _fail("retry ledger cannot change accepted macro-step count")
    if live.accumulated_debit_vector != accepted.accumulated_debit_vector:
        _fail("retry ledger cannot rewrite accumulated debit")
    if (
        live.last_accepted_macro_step_temporal_retry_count
        != accepted.last_accepted_macro_step_temporal_retry_count
    ):
        _fail("retry ledger cannot rewrite last-accepted temporal retry history")
    accepted_rejections = accepted.serialized_temporal_rejections
    live_rejections = live.serialized_temporal_rejections
    if live_rejections[: len(accepted_rejections)] != accepted_rejections:
        _fail("retry ledger is not a prefix-preserving suffix of the accepted ledger")
    appended = live_rejections[len(accepted_rejections) :]
    if not appended:
        _fail("retry ledger changed without appending a current-step rejection")
    if live.current_macro_step_temporal_retry_count != len(appended):
        _fail("appended rejection count differs from the current-step temporal counter")
    if live.cumulative_temporal_retry_count != (
        accepted.cumulative_temporal_retry_count + len(appended)
    ):
        _fail("cumulative temporal retry count is not the accepted count plus suffix")


@dataclass(frozen=True, slots=True)
class HLT17GenerationBundle:
    """Validated descriptor, payload, and full-cursor bytes for later publication.

    Hashes, revisions, projections, and authorization flags are derived from
    those bytes after codec and cursor reauthentication. They cannot be forged
    through construction or dataclass replace.
    """

    descriptor: bytes
    payload: bytes
    cursor_bytes: bytes
    _encoding: HLT17MemberEncoding = field(init=False, repr=False, compare=False)
    _cursor: HLT17IMP1Cursor = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        try:
            encoding = HLT17MemberEncoding(self.descriptor, self.payload)
        except (HLT17MemberCodecError, TypeError, ValueError) as error:
            raise HLT17CheckpointError(
                "generation bundle encoding failed codec reauthentication"
            ) from error
        try:
            cursor_bytes, cursor = _restore_cursor_bytes(self.cursor_bytes)
        except (HLT17CheckpointError, TypeError, ValueError):
            raise
        if cursor.descriptor_sha256 != encoding.descriptor_sha256:
            _fail("full cursor does not name the accepted descriptor")
        if cursor.physical_state_sha256 != encoding.physical_state_sha256:
            _fail("full cursor physical hash differs from the accepted encoding")
        _require_legal_retry_suffix(encoding.decoded.ledger, cursor.ledger)
        metadata = encoding.metadata
        accepted_generation = int(metadata["runtime_identity"]["accepted_generation"])
        _require_revision_relationship(
            accepted_generation=accepted_generation, cursor=cursor
        )
        if metadata["snapshot_role"] != SNAPSHOT_ROLE_ACCEPTED:
            _fail("generation bundle must export an accepted physical-state descriptor")
        if bool(metadata["campaign_execution_authorized"]) is not False:
            _fail("generation bundle cannot authorize a campaign")
        if bool(metadata["historical_origin_authenticated"]) is not False:
            _fail("generation bundle cannot authenticate origin")
        if bool(metadata["durable_store_bound"]) is not False:
            _fail("generation bundle cannot bind a durable store")
        expected_identity = hlt17_implementation_identity()
        runtime_identity = metadata["runtime_identity"]
        for name, value in expected_identity.items():
            if runtime_identity.get(name) != value:
                _fail(f"generation bundle {name} is not the fixed C1R1 identity")
            if getattr(cursor, name) != value:
                _fail(f"generation bundle cursor {name} is not the fixed C1R1 identity")
        object.__setattr__(self, "descriptor", bytes(encoding.descriptor))
        object.__setattr__(self, "payload", bytes(encoding.payload))
        object.__setattr__(self, "cursor_bytes", cursor_bytes)
        object.__setattr__(self, "_encoding", encoding)
        object.__setattr__(self, "_cursor", cursor)

    @property
    def descriptor_sha256(self) -> str:
        return self._encoding.descriptor_sha256

    @property
    def payload_sha256(self) -> str:
        return self._encoding.payload_sha256

    @property
    def physical_state_sha256(self) -> str:
        return self._encoding.physical_state_sha256

    @property
    def cursor_mapping(self) -> dict[str, object]:
        return encode_hlt17_cursor(self._cursor)

    @property
    def cursor_sha256(self) -> str:
        return cursor_sha256(self._cursor)

    @property
    def accepted_generation(self) -> int:
        return int(self._encoding.metadata["runtime_identity"]["accepted_generation"])

    @property
    def implementation_id(self) -> str:
        return str(self._cursor.implementation_id)

    @property
    def admission_runtime_id(self) -> str:
        return str(self._cursor.admission_runtime_id)

    @property
    def cursor_generation(self) -> int:
        return self._cursor.cursor_generation

    @property
    def attempt_serial(self) -> int:
        return self._cursor.attempt_serial

    @property
    def kernel_checkpoint_generation(self) -> int:
        return self._cursor.checkpoint_generation

    @property
    def kernel_transition_parent_sha256(self) -> str:
        return self._cursor.kernel_transition_parent_sha256

    @property
    def kernel_transition_digest(self) -> str:
        return self._cursor.kernel_transition_digest

    @property
    def predecessor_cursor_identity_sha256(self) -> str:
        return str(self._encoding.metadata["predecessor_cursor_identity_sha256"])

    @property
    def predecessor_cursor_sha256(self) -> str:
        return str(self._encoding.metadata["predecessor_cursor_sha256"])

    @property
    def store_journal_parent_sha256(self) -> None:
        return None

    @property
    def store_journal_tip_sha256(self) -> None:
        return None

    @property
    def store_checkpoint_revision(self) -> None:
        return None

    @property
    def durable_store_bound(self) -> bool:
        return False

    @property
    def campaign_execution_authorized(self) -> bool:
        return False

    @property
    def historical_origin_authenticated(self) -> bool:
        return False

    @property
    def publication_authorized(self) -> bool:
        return False

    @property
    def atomic_publication_seam(self) -> str:
        return FUTURE_ATOMIC_GENERATION_BUNDLE_SEAM


@dataclass(frozen=True, slots=True)
class HLT17CheckpointResult:
    disposition: str
    cursor: HLT17IMP1Cursor
    encoding: HLT17MemberEncoding
    bundle: HLT17GenerationBundle
    prepared: TDG11C1R1PreparedRuntime | None = None
    committed: TDG11C1R1CommittedRuntime | None = None
    bridge: HLT17BridgeResult | None = None
    last_accepted_executed_plan: TDG7Binary64SubdivisionPlan | None = None
    parked: bool = False
    accepted_state_advanced: bool = False
    evidence: Mapping[str, object] = field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(self, "evidence", dict(self.evidence))
        encoding = _as_encoding(self.encoding)
        object.__setattr__(self, "encoding", encoding)
        cursor = revalidate_hlt17_cursor(self.cursor)
        object.__setattr__(self, "cursor", cursor)
        expected = _bundle_from(encoding=encoding, cursor=cursor)
        if (
            self.bundle.descriptor != expected.descriptor
            or self.bundle.payload != expected.payload
            or self.bundle.cursor_bytes != expected.cursor_bytes
        ):
            _fail("result bundle does not match the encoding and full-cursor bytes")
        if self.last_accepted_executed_plan != encoding.executed_plan:
            _fail("result last accepted plan differs from the accepted descriptor")
        if type(self.parked) is not bool or type(self.accepted_state_advanced) is not bool:
            raise TypeError("parked and accepted_state_advanced must be built-in bools")
        if self.parked and self.accepted_state_advanced:
            _fail("parking at the event target cannot advance accepted state")
        if self.disposition == PARKED_AT_EVENT_TARGET and self.parked is not True:
            _fail("parked_at_event_target requires parked=True")
        if self.disposition == ACCEPTED_FINE:
            if self.committed is None or self.accepted_state_advanced is not True:
                _fail("accepted_fine requires the committed runtime and the advance flag")


def _bundle_from(
    *,
    encoding: HLT17MemberEncoding,
    cursor: HLT17IMP1Cursor,
) -> HLT17GenerationBundle:
    encoding = _as_encoding(encoding)
    return HLT17GenerationBundle(
        descriptor=encoding.descriptor,
        payload=encoding.payload,
        cursor_bytes=_full_cursor_bytes(cursor),
    )


class HLT17InMemoryCheckpoint:
    """Joined no-I/O member+cursor checkpoint. Not a store, origin, or campaign."""

    def __init__(
        self,
        member: HLT17RuntimeMember,
        cursor: HLT17IMP1Cursor,
        encoding: HLT17MemberEncoding,
    ) -> None:
        if not isinstance(member, HLT17RuntimeMember):
            raise TypeError("member must be HLT17RuntimeMember")
        if not isinstance(cursor, HLT17IMP1Cursor):
            raise TypeError("cursor must be HLT17IMP1Cursor")
        encoding = _as_encoding(encoding)
        self.member = member
        self._cursor = revalidate_hlt17_cursor(cursor)
        self._descriptor = encoding.descriptor
        self._payload = encoding.payload
        self.agree()

    @property
    def cursor(self) -> HLT17IMP1Cursor:
        return revalidate_hlt17_cursor(self._cursor)

    @cursor.setter
    def cursor(self, cursor: HLT17IMP1Cursor) -> None:
        self._cursor = revalidate_hlt17_cursor(cursor)

    @property
    def accepted_encoding(self) -> HLT17MemberEncoding:
        return HLT17MemberEncoding(self._descriptor, self._payload)

    @accepted_encoding.setter
    def accepted_encoding(self, encoding: HLT17MemberEncoding) -> None:
        encoding = _as_encoding(encoding)
        self._descriptor = encoding.descriptor
        self._payload = encoding.payload

    @property
    def accepted_generation(self) -> int:
        return int(self.accepted_encoding.metadata["runtime_identity"]["accepted_generation"])

    @property
    def last_accepted_executed_plan(self) -> TDG7Binary64SubdivisionPlan | None:
        return self.accepted_encoding.executed_plan

    @classmethod
    def seed(
        cls,
        member: HLT17RuntimeMember,
        *,
        event_target: object,
        origin_predecessor: HLT17PredecessorIdentity | None = None,
        fault_hook: Callable[[str], None] | None = None,
    ) -> "HLT17InMemoryCheckpoint":
        """Encode the accepted boundary, then seed the cursor with that digest."""

        if not isinstance(member, HLT17RuntimeMember):
            raise TypeError("member must be HLT17RuntimeMember")
        member.__post_init__()
        if origin_predecessor is None:
            if member.identity.scope == IDENTITY_SCOPE_LIVE:
                _fail("live seeding requires an explicitly supplied origin predecessor")
            predecessor = synthetic_origin_predecessor_identity(member)
        else:
            predecessor = origin_predecessor
        if not isinstance(predecessor, HLT17PredecessorIdentity):
            raise TypeError("origin_predecessor must be HLT17PredecessorIdentity")
        if predecessor.cursor_identity["protocol_artifact_id"] != ORIGIN_PROTOCOL:
            _fail("initial encoding must name the sealed PROTO17 origin identity")
        encoding = encode_member(member, predecessor, generation=0)
        if fault_hook is not None:
            fault_hook("after_initial_encoding")
        if encoding.descriptor_sha256 != sha256(encoding.descriptor).hexdigest():
            _fail("initial encoding descriptor hash is not redundant")
        identities = live_boundary_identities(
            state=member.state,
            rhs=member.operator,
            projector=member.projector,
            transaction=member.transaction,
            tracers=member.tracers,
            coordinates=member.coordinates,
        )
        if identities["physical_state_sha256"] != encoding.physical_state_sha256:
            _fail("initial encoding physical hash differs from the live member")
        cursor = seed_hlt17_cursor(
            member_key=member.key,
            event_target=event_target,
            ledger=member.temporal_ledger,
            descriptor_sha256=encoding.descriptor_sha256,
            physical_state_sha256=identities["physical_state_sha256"],
            coordinates_sha256=identities["coordinates_sha256"],
            synthetic_callable_binding_sha256=identities[
                "synthetic_callable_binding_sha256"
            ],
            monitor_sha256=identities["monitor_sha256"],
            causal_state_sha256=identities["causal_state_sha256"],
            tracer_sha256=identities["tracer_sha256"],
            guard_configuration_sha256=identities["guard_configuration_sha256"],
            source_total=member.source_retry_count,
            cfl_total=member.CFL_retry_count,
        )
        if cursor.descriptor_sha256 != encoding.descriptor_sha256:
            _fail("seeded cursor did not receive the actual descriptor digest")
        if cursor.source_total != member.source_retry_count:
            _fail("seeded cursor zeroed the inherited source baseline")
        if cursor.cfl_total != member.CFL_retry_count:
            _fail("seeded cursor zeroed the inherited CFL baseline")
        return cls(member, cursor, encoding)

    def live_boundary(self) -> HLT17LiveBoundary:
        self.member.__post_init__()
        revalidate_hlt17_cursor(self.cursor)
        return HLT17LiveBoundary(
            member_key=self.member.key,
            state=self.member.state,
            rhs=self.member.operator,
            projector=self.member.projector,
            transaction=self.member.transaction,
            tracers=self.member.tracers,
            coordinates=self.member.coordinates,
            temporal_ledger=self.member.temporal_ledger,
            previous_step_index=self.member.step_index,
            previous_transaction_serial=self.member.transaction_serial,
            descriptor_sha256=self.cursor.descriptor_sha256,
            event_target=self.cursor.event_target,
            limits=imp1_enclosure_limits_for_member(self.member.key),
        )

    def at_event_target(self) -> bool:
        return _same_bits(self.member.time, self.cursor.event_target)

    def generation_bundle(self) -> HLT17GenerationBundle:
        self.agree()
        return _bundle_from(encoding=self.accepted_encoding, cursor=self.cursor)

    def remaining_seams(self) -> dict[str, str]:
        return {
            "content_addressed_store": FUTURE_CONTENT_ADDRESSED_STORE_SEAM,
            "checkpoint_journal": FUTURE_CHECKPOINT_JOURNAL_SEAM,
            "authenticated_origin": FUTURE_AUTHENTICATED_ORIGIN_SEAM,
            "campaign_transaction": FUTURE_CAMPAIGN_TRANSACTION_SEAM,
            "source_manifest": FUTURE_SOURCE_MANIFEST_SEAM,
            "atomic_generation_bundle": FUTURE_ATOMIC_GENERATION_BUNDLE_SEAM,
            "bridge_store": FUTURE_INTEGRATION_SEAMS["content_addressed_store"],
            "store_implementation_identity": FUTURE_STORE_IMPLEMENTATION_SEAM,
        }

    def agree(self) -> None:
        """Reauthenticate encoding and cross-check the live boundary with the cursor."""

        member = self.member
        try:
            member.__post_init__()
        except HLT17RuntimeMemberError as error:
            raise HLT17CheckpointError(str(error)) from error
        try:
            encoding = HLT17MemberEncoding(self._descriptor, self._payload)
        except (HLT17MemberCodecError, TypeError, ValueError) as error:
            raise HLT17CheckpointError(
                "accepted encoding failed codec reauthentication"
            ) from error
        cursor = revalidate_hlt17_cursor(self._cursor)
        try:
            authenticate_live_boundary(cursor, self.live_boundary())
        except (HLT17InvalidPremise, HLT17IMP1Error) as error:
            raise HLT17CheckpointError(str(error)) from error
        decoded = encoding.decoded
        metadata = dict(decoded.metadata)
        if member.key != cursor.member_key:
            _fail("member/cursor identity differs")
        if member.source_retry_count != cursor.source_total:
            _fail("member source total differs from the cursor source total")
        if member.CFL_retry_count != cursor.cfl_total:
            _fail("member CFL total differs from the cursor CFL total")
        encoding_source = int(metadata["counters"]["source_retry_count"])
        encoding_cfl = int(metadata["counters"]["CFL_retry_count"])
        if cursor.source_total < encoding_source:
            _fail("cursor source total lost the accepted source baseline")
        if cursor.cfl_total < encoding_cfl:
            _fail("cursor CFL total lost the accepted CFL baseline")
        if not _same_bits(member.time, cursor.accepted_time):
            _fail("member accepted time differs from the cursor")
        if member.step_index != cursor.accepted_step_index:
            _fail("member step index differs from the cursor")
        if member.transaction_serial != cursor.accepted_transaction_serial:
            _fail("member transaction serial differs from the cursor")
        if member.physical_state_sha256() != cursor.physical_state_sha256:
            _fail("member physical hash differs from the cursor")
        if member.physical_state_sha256() != encoding.physical_state_sha256:
            _fail("member physical hash differs from the accepted encoding")
        if decoded.physical_state_sha256 != encoding.physical_state_sha256:
            _fail("decoded physical hash differs from the encoding")
        if cursor.descriptor_sha256 != encoding.descriptor_sha256:
            _fail("cursor does not name the accepted descriptor")
        if encoding.descriptor_sha256 == encoding.physical_state_sha256:
            _fail("descriptor and physical identities are not interchangeable")
        identities = live_boundary_identities(
            state=member.state,
            rhs=member.operator,
            projector=member.projector,
            transaction=member.transaction,
            tracers=member.tracers,
            coordinates=member.coordinates,
        )
        for name in LIVE_BOUNDARY_IDENTITY_NAMES:
            if identities[name] != getattr(cursor, name):
                _fail(f"live {name} differs from the HLT17 cursor")
        if not np.array_equal(member.state.u, decoded.arrays["u"]):
            _fail("live u differs from the accepted payload")
        if not np.array_equal(member.state.p, decoded.arrays["p"]):
            _fail("live p differs from the accepted payload")
        if not np.array_equal(member.state.q, decoded.arrays["q"]):
            _fail("live q differs from the accepted payload")
        if not np.array_equal(member.coordinates, decoded.arrays["grid_coordinates"]):
            _fail("live coordinates differ from the accepted payload")
        if not np.array_equal(member.tracers.labels, decoded.arrays["tracer_labels"]):
            _fail("live tracer labels differ from the accepted payload")
        if not np.array_equal(member.tracers.positions, decoded.arrays["tracer_positions"]):
            _fail("live tracer positions differ from the accepted payload")
        if not np.array_equal(
            member.tracers.proper_times, decoded.arrays["tracer_proper_times"]
        ):
            _fail("live tracer proper times differ from the accepted payload")
        if not np.array_equal(
            np.stack(tuple(member.tracers.event_proper_times)),
            decoded.arrays["event_proper_times"],
        ):
            _fail("live tracer event proper times differ from the accepted payload")
        if not np.array_equal(
            np.stack(tuple(member.tracers.event_fields)),
            decoded.arrays["event_fields"],
        ):
            _fail("live tracer event fields differ from the accepted payload")
        accepted_generation = int(metadata["runtime_identity"]["accepted_generation"])
        _require_revision_relationship(
            accepted_generation=accepted_generation, cursor=cursor
        )
        if ledger_sha256(member.temporal_ledger) != ledger_sha256(cursor.ledger):
            _fail("member ledger differs from the cursor ledger")
        _require_legal_retry_suffix(decoded.ledger, cursor.ledger)
        if member.temporal_ledger.inherited_snapshot != cursor.ledger.inherited_snapshot:
            _fail("member/cursor inherited TDG6 sibling differs")
        if member.executed_plan != encoding.executed_plan:
            _fail("live executed plan differs from the accepted descriptor")
        if metadata["source_operator_closure"] != member.source_operator_closure.as_mapping():
            _fail("live source/operator closure differs from the accepted descriptor")
        if metadata["origin_reference"] != member.origin_reference.as_mapping():
            _fail("live origin reference differs from the accepted descriptor")
        expected_identity = hlt17_implementation_identity()
        runtime_identity = metadata["runtime_identity"]
        for name, value in expected_identity.items():
            if runtime_identity.get(name) != value:
                _fail(f"accepted encoding {name} is not the fixed C1R1 identity")
            if getattr(cursor, name) != value:
                _fail(f"cursor {name} is not the fixed C1R1 identity")

    def _parked_result(self) -> HLT17CheckpointResult:
        self.agree()
        return HLT17CheckpointResult(
            disposition=PARKED_AT_EVENT_TARGET,
            cursor=self.cursor,
            encoding=self.accepted_encoding,
            bundle=self.generation_bundle(),
            last_accepted_executed_plan=self.last_accepted_executed_plan,
            parked=True,
            accepted_state_advanced=False,
            evidence={
                "kind": PARKED_AT_EVENT_TARGET,
                "accepted_time_hex": self.member.time.hex(),
                "event_target_hex": self.cursor.event_target.hex(),
                "source_invoked": False,
                "retarget_authorized": False,
                "campaign_execution_authorized": False,
            },
        )

    def prepare(
        self,
        *,
        requested_cap: object,
        replay: object | None = None,
    ) -> HLT17BridgeResult | HLT17CheckpointResult:
        """Prepare the next attempt, or park at the current event target."""

        self.agree()
        if self.at_event_target():
            return self._parked_result()
        boundary = self.live_boundary()
        authenticate_live_boundary(self.cursor, boundary)
        if self.cursor.public_fresh_permitted():
            return prepare_hlt17_fresh(
                self.cursor, boundary, requested_cap=requested_cap
            )
        return prepare_hlt17_next(self.cursor, boundary, replay=replay)

    def require_admission(
        self, prepared: TDG11C1R1PreparedRuntime
    ) -> HLT17BridgeResult | HLT17CheckpointResult:
        self.agree()
        if self.at_event_target():
            return self._parked_result()
        boundary = self.live_boundary()
        return require_hlt17_admission(self.cursor, boundary, prepared)

    def absorb_nonfine(
        self,
        result: HLT17BridgeResult,
        *,
        fault_hook: Callable[[str], None] | None = None,
    ) -> HLT17CheckpointResult:
        """Install a source/CFL/temporal cursor overlay without advancing state."""

        if not isinstance(result, HLT17BridgeResult):
            raise TypeError("result must be HLT17BridgeResult")
        if result.disposition == ACCEPTED_FINE or result.accepted_state_advanced:
            _fail("fine adoption must use adopt_fine, not absorb_nonfine")
        prior = self.member.capture()
        prior_cursor = self.cursor
        try:
            self.member.install_retry_overlay_state(
                result.cursor.ledger,
                source_retry_count=result.cursor.source_total,
                CFL_retry_count=result.cursor.cfl_total,
                executed_plan=self.last_accepted_executed_plan,
                fault_hook=fault_hook,
            )
            self.cursor = revalidate_hlt17_cursor(result.cursor)
            if self.member.executed_plan != self.last_accepted_executed_plan:
                _fail("retry overlay dropped the last accepted executed plan")
            self.agree()
            if fault_hook is not None:
                fault_hook("after_overlay_result")
            return HLT17CheckpointResult(
                disposition=result.disposition,
                cursor=self.cursor,
                encoding=self.accepted_encoding,
                bundle=self.generation_bundle(),
                prepared=result.prepared,
                bridge=result,
                last_accepted_executed_plan=self.last_accepted_executed_plan,
                parked=False,
                accepted_state_advanced=False,
                evidence=dict(result.evidence),
            )
        except BaseException as error:
            try:
                self.member.restore(prior)
                self.cursor = prior_cursor
                self.member.__post_init__()
            except BaseException as rollback_error:
                raise HLT17CheckpointRollbackError(
                    "HLT17 complete overlay rollback failed"
                ) from rollback_error
            if isinstance(error, (HLT17RuntimeMemberError, HLT17RollbackError)):
                raise HLT17CheckpointError(str(error)) from error
            raise

    def adopt_fine(
        self,
        prepared: TDG11C1R1PreparedRuntime,
        *,
        fault_hook: Callable[[str], None] | None = None,
    ) -> HLT17CheckpointResult:
        """Commit once on an isolated copy, encode, finalize, then adopt live."""

        self.agree()
        if self.at_event_target():
            _fail("cannot adopt a fine step after the member is parked at the event target")
        prepared = require_hlt17_prepared(prepared)
        prior = self.member.capture()
        prior_cursor = self.cursor
        prior_encoding = self.accepted_encoding
        old_cursor = self.cursor
        try:
            provisional = _isolated_member(self.member)
            boundary = HLT17LiveBoundary(
                member_key=provisional.key,
                state=provisional.state,
                rhs=provisional.operator,
                projector=provisional.projector,
                transaction=provisional.transaction,
                tracers=provisional.tracers,
                coordinates=provisional.coordinates,
                temporal_ledger=provisional.temporal_ledger,
                previous_step_index=provisional.step_index,
                previous_transaction_serial=provisional.transaction_serial,
                descriptor_sha256=old_cursor.descriptor_sha256,
                event_target=old_cursor.event_target,
                limits=imp1_enclosure_limits_for_member(provisional.key),
            )
            authenticate_live_boundary(old_cursor, boundary)
            provisional_prior = provisional.capture()
            committed = commit_hlt17_kernel(old_cursor, boundary, prepared)
            if fault_hook is not None:
                fault_hook("after_kernel_adoption")
            committed = require_hlt17_committed(committed)
            provisional.adopt_committed(
                committed,
                executed_plan=prepared.plan,
                source_retry_count=old_cursor.source_total,
                CFL_retry_count=old_cursor.cfl_total,
                rollback_to=provisional_prior,
            )
            predecessor = predecessor_from_accepted_cursor(
                old_cursor,
                campaign_id=self.member.identity.campaign_id,
                accepted_generation=self.accepted_generation,
                descriptor_sha256=old_cursor.descriptor_sha256,
            )
            encoding = encode_member(
                provisional,
                predecessor,
                generation=self.accepted_generation + 1,
            )
            if fault_hook is not None:
                fault_hook("after_provisional_encoding")
            finalized = finalize_hlt17_fine(
                old_cursor,
                boundary,
                prepared,
                committed,
                descriptor_sha256=encoding.descriptor_sha256,
            )
            if fault_hook is not None:
                fault_hook("after_cursor_finalization")
            new_cursor = finalized.cursor
            _require_no_hash_cycle(
                old_cursor=old_cursor, new_cursor=new_cursor, encoding=encoding
            )
            provisional.source_retry_count = new_cursor.source_total
            provisional.CFL_retry_count = new_cursor.cfl_total
            provisional.__post_init__()
            bundle = _bundle_from(encoding=encoding, cursor=new_cursor)
            result = HLT17CheckpointResult(
                disposition=ACCEPTED_FINE,
                cursor=new_cursor,
                encoding=encoding,
                bundle=bundle,
                prepared=prepared,
                committed=committed,
                bridge=finalized,
                last_accepted_executed_plan=encoding.executed_plan,
                parked=False,
                accepted_state_advanced=True,
                evidence=dict(finalized.evidence),
            )
            if fault_hook is not None:
                fault_hook("after_result_construction")
            self.member.restore(provisional.capture())
            self.cursor = new_cursor
            self.accepted_encoding = encoding
            self.agree()
            return result
        except BaseException as error:
            try:
                self.member.restore(prior)
                self.cursor = prior_cursor
                self.accepted_encoding = prior_encoding
                self.member.__post_init__()
            except BaseException as rollback_error:
                raise HLT17CheckpointRollbackError(
                    "HLT17 complete fine-adoption rollback failed"
                ) from rollback_error
            if isinstance(error, (HLT17RuntimeMemberError, HLT17RollbackError, HLT17IMP1Error)):
                raise HLT17CheckpointError(str(error)) from error
            raise

    def restore(
        self,
        encoding: HLT17MemberEncoding,
        cursor: HLT17IMP1Cursor | Mapping[str, object],
        *,
        last_accepted_executed_plan: TDG7Binary64SubdivisionPlan | None = None,
        fault_hook: Callable[[str], None] | None = None,
    ) -> HLT17MemberEncoding:
        """Restore accepted descriptor and full cursor, then apply the overlay."""

        return restore_hlt17_checkpoint(
            self,
            encoding,
            cursor,
            last_accepted_executed_plan=last_accepted_executed_plan,
            fault_hook=fault_hook,
        )

    def advance_once(
        self,
        *,
        requested_cap: object,
        fault_hook: Callable[[str], None] | None = None,
    ) -> HLT17CheckpointResult:
        """Single prepare -> admit/overlay/fine path. No duplicate admission policy."""

        prepared_result = self.prepare(requested_cap=requested_cap)
        if isinstance(prepared_result, HLT17CheckpointResult):
            return prepared_result
        if prepared_result.disposition in {
            "source_retry_required",
            "cfl_retry_required",
            "source_retry_exhausted",
            "cfl_retry_exhausted",
            "invalid_premise",
        }:
            return self.absorb_nonfine(prepared_result, fault_hook=fault_hook)
        if prepared_result.prepared is None:
            _fail(f"{prepared_result.disposition} omitted a prepared family")
        admitted = self.require_admission(prepared_result.prepared)
        if isinstance(admitted, HLT17CheckpointResult):
            return admitted
        if admitted.disposition in {
            "temporal_retry_required",
            "temporal_retry_exhausted",
            "invalid_premise",
            "source_retry_required",
            "cfl_retry_required",
            "source_retry_exhausted",
            "cfl_retry_exhausted",
        }:
            return self.absorb_nonfine(admitted, fault_hook=fault_hook)
        if admitted.prepared is None:
            _fail("admitted result omitted a prepared family")
        return self.adopt_fine(admitted.prepared, fault_hook=fault_hook)


def restore_hlt17_checkpoint(
    checkpoint: HLT17InMemoryCheckpoint,
    encoding: HLT17MemberEncoding,
    cursor: HLT17IMP1Cursor | Mapping[str, object],
    *,
    last_accepted_executed_plan: TDG7Binary64SubdivisionPlan | None = None,
    fault_hook: Callable[[str], None] | None = None,
) -> HLT17MemberEncoding:
    """Authenticate accepted bytes and the full cursor, then apply overlay."""

    if not isinstance(checkpoint, HLT17InMemoryCheckpoint):
        raise TypeError("checkpoint must be HLT17InMemoryCheckpoint")
    if not isinstance(encoding, HLT17MemberEncoding):
        raise TypeError("encoding must be HLT17MemberEncoding")
    prior = checkpoint.member.capture()
    prior_cursor = checkpoint.cursor
    prior_encoding = checkpoint.accepted_encoding
    try:
        encoding = _as_encoding(encoding)
        decoded = encoding.decoded
        if last_accepted_executed_plan is not None and (
            last_accepted_executed_plan != decoded.executed_plan
        ):
            _fail(
                "caller-supplied last accepted plan differs from the accepted descriptor"
            )
        accepted_plan = decoded.executed_plan
        restored_cursor = (
            cursor
            if isinstance(cursor, HLT17IMP1Cursor)
            else restore_hlt17_cursor(cursor)
        )
        restored_cursor = revalidate_hlt17_cursor(restored_cursor)
        metadata = dict(decoded.metadata)
        if restored_cursor.descriptor_sha256 != encoding.descriptor_sha256:
            _fail("restored cursor does not name the accepted descriptor")
        if restored_cursor.physical_state_sha256 != encoding.physical_state_sha256:
            _fail("restored cursor physical hash differs from the accepted encoding")
        if restored_cursor.member_key != checkpoint.member.key:
            _fail("restored cursor member differs")
        if metadata["source_operator_closure"] != (
            checkpoint.member.source_operator_closure.as_mapping()
        ):
            _fail("restored encoding source/operator closure differs")
        if metadata["origin_reference"] != checkpoint.member.origin_reference.as_mapping():
            _fail("restored encoding origin reference differs")
        _require_legal_retry_suffix(decoded.ledger, restored_cursor.ledger)
        restore_member(checkpoint.member, encoding, fault_hook=fault_hook)
        encoding_source = int(metadata["counters"]["source_retry_count"])
        encoding_cfl = int(metadata["counters"]["CFL_retry_count"])
        if restored_cursor.source_total < encoding_source:
            _fail("restored cursor lost the inherited source baseline")
        if restored_cursor.cfl_total < encoding_cfl:
            _fail("restored cursor lost the inherited CFL baseline")
        checkpoint.member.install_retry_overlay_state(
            restored_cursor.ledger,
            source_retry_count=restored_cursor.source_total,
            CFL_retry_count=restored_cursor.cfl_total,
            executed_plan=accepted_plan,
            fault_hook=fault_hook,
        )
        checkpoint.cursor = restored_cursor
        checkpoint.accepted_encoding = encoding
        if checkpoint.member.executed_plan != accepted_plan:
            _fail("restore dropped the last accepted executed plan")
        checkpoint.agree()
        return encoding
    except BaseException as error:
        try:
            checkpoint.member.restore(prior)
            checkpoint.cursor = prior_cursor
            checkpoint.accepted_encoding = prior_encoding
            checkpoint.member.__post_init__()
        except BaseException as rollback_error:
            raise HLT17CheckpointRollbackError(
                "HLT17 complete checkpoint restore rollback failed"
            ) from rollback_error
        if isinstance(
            error, (HLT17RuntimeMemberError, HLT17RollbackError, HLT17IMP1Error)
        ):
            raise HLT17CheckpointError(str(error)) from error
        raise


__all__ = [
    "ACCEPTED_FINE",
    "CHECKPOINT_PROTOCOL",
    "FUTURE_ATOMIC_GENERATION_BUNDLE_SEAM",
    "HLT17CheckpointError",
    "HLT17CheckpointResult",
    "HLT17CheckpointRollbackError",
    "HLT17GenerationBundle",
    "HLT17InMemoryCheckpoint",
    "IDENTITY_LABELS",
    "PARKED_AT_EVENT_TARGET",
    "PINNED_PRIVATE_DEPENDENCIES",
    "SEEDED",
    "SYNTHETIC_ORIGIN_CURSOR_PLACEHOLDER_DOMAIN",
    "origin_predecessor_identity",
    "predecessor_from_accepted_cursor",
    "restore_hlt17_checkpoint",
    "synthetic_origin_predecessor_identity",
]
