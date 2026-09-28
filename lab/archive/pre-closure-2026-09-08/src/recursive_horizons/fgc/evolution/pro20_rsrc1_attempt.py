"""Store-blind one-attempt reducer for prospective PRO20 RSRC1.

This module is intentionally independent of every campaign-store module.  It
restores one accepted HLT17 bundle into one reconstructed live member, executes
one scheduled prepare/admit/adopt path, and derives the structural result from
bytes.  Publication, scheduling across members, authority, and endpoint use
belong to later parent layers.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Mapping

from recursive_horizons.evidence_io import load_canonical_json

from .hlt17_admission_runtime import (
    C1R1CoordinateLatticeStop,
    C1R1RefinementPathStop,
    C1R1RuntimeResourceStop,
    C1R1TemporalRetryExhausted,
)
from .hlt17_imp1_cursor import (
    CFL_OWNER,
    FRESH_READY,
    HLT17_CFL_RETRY_CAP,
    HLT17InvalidPremise,
    HLT17OwnerRetryExhausted,
    HLT17_SOURCE_RETRY_CAP,
    HLT17IMP1Cursor,
    RETRY_PENDING,
    SOURCE_OWNER,
    ledger_sha256,
    revalidate_hlt17_cursor,
    restore_hlt17_cursor,
)
from .hlt17_member_checkpoint import (
    ACCEPTED_FINE,
    HLT17BridgeResult,
    HLT17CheckpointResult,
    HLT17GenerationBundle,
    HLT17InMemoryCheckpoint,
)
from .hlt17_member_codec import HLT17MemberEncoding
from .pro20_origin import MEMBER_KEYS
from .protocol_v19 import (
    BLOB_UNCHANGED_KINDS,
    CURSOR_CHANGED_RETRY_KINDS,
    RECORD_KINDS,
    TEMPORAL_EXHAUSTION_REASONS,
    TERMINAL_KINDS,
)
from .tdg11_imp1_ledger import imp1_checkpoint_extension


ARTIFACT_ID = "FGC-1-PRO20-EV1-RSRC1"
SOURCE_RETRY_CAP = HLT17_SOURCE_RETRY_CAP
CFL_RETRY_CAP = HLT17_CFL_RETRY_CAP


class RSRC1AttemptError(ValueError):
    """A scheduled attempt or its immutable byte relation differs."""


class RSRC1AttemptPremiseError(RSRC1AttemptError):
    """A declared attempt premise is absent or inconsistent."""


class RSRC1AttemptUnexpectedError(RuntimeError):
    """The child produced an unclassifiable or state-unsafe relation."""


def _fail(message: str, cls: type[Exception] = RSRC1AttemptError) -> None:
    raise cls(message)


def _same_bits(left: float, right: float) -> bool:
    return float(left).hex() == float(right).hex()


def _cursor_from_bundle(bundle: HLT17GenerationBundle) -> HLT17IMP1Cursor:
    mapping = load_canonical_json(bundle.cursor_bytes)
    if type(mapping) is not dict:
        _fail("bundle cursor is not a JSON object", RSRC1AttemptPremiseError)
    return restore_hlt17_cursor(mapping)


def _bundles_identical(
    left: HLT17GenerationBundle,
    right: HLT17GenerationBundle,
) -> bool:
    return (
        left.descriptor == right.descriptor
        and left.payload == right.payload
        and left.cursor_bytes == right.cursor_bytes
    )


def _encoding_changed(
    left: HLT17GenerationBundle,
    right: HLT17GenerationBundle,
) -> bool:
    return left.descriptor != right.descriptor or left.payload != right.payload


def fresh_requested_cap(checkpoint: HLT17InMemoryCheckpoint) -> float:
    """Return the unchanged physical CFL cap formula for a fresh cursor."""

    if not isinstance(checkpoint, HLT17InMemoryCheckpoint):
        raise TypeError("checkpoint must be HLT17InMemoryCheckpoint")
    checkpoint.agree()
    cursor = checkpoint.cursor
    if not cursor.public_fresh_permitted():
        _fail("fresh cap is not available on a pending cursor", RSRC1AttemptPremiseError)
    remaining = float(cursor.event_target) - float(cursor.accepted_time)
    if not math.isfinite(remaining) or remaining <= 0.0:
        _fail("fresh cap requires positive remaining time", RSRC1AttemptPremiseError)
    transaction = checkpoint.member.transaction
    spacing = float(transaction.grid_spacing)
    cfl_maximum = float(transaction.cfl_maximum)
    speed = float(transaction.causal_state.previous_speed_upper)
    if not all(
        math.isfinite(value) and value > 0.0
        for value in (spacing, cfl_maximum, speed)
    ):
        _fail("fresh cap inputs are not positive", RSRC1AttemptPremiseError)
    cfl_cap = cfl_maximum * spacing / speed
    if not math.isfinite(cfl_cap) or cfl_cap <= 0.0:
        _fail("CFL cap is not a positive finite width", RSRC1AttemptPremiseError)
    return remaining if remaining < cfl_cap else cfl_cap


def pending_plan_from_cursor(cursor: HLT17IMP1Cursor):
    """Return only the exact cursor-owned pending plan."""

    if not isinstance(cursor, HLT17IMP1Cursor):
        raise TypeError("cursor must be HLT17IMP1Cursor")
    cursor = revalidate_hlt17_cursor(cursor)
    if cursor.public_fresh_permitted():
        return None
    plan = cursor.next_attempt_plan()
    if plan is None:
        _fail(
            "pending retry cursor has no executable next plan",
            RSRC1AttemptPremiseError,
        )
    return plan


@dataclass(frozen=True, slots=True)
class AttemptClassification:
    """Protocol kind derived from bytes and typed stops, never a child label."""

    kind: str
    member_key: str
    successor_bundle: HLT17GenerationBundle
    terminal_evidence: Mapping[str, object] | None
    accepted_state_advanced: bool
    executable_retry_cursor: bool
    physical_state_preserved: bool
    structural_class: str

    def __post_init__(self) -> None:
        if self.kind not in RECORD_KINDS or self.kind == "seed":
            _fail("classified kind is not a production attempt kind", RSRC1AttemptPremiseError)
        if self.member_key not in MEMBER_KEYS:
            _fail("classified member is outside the live cohort", RSRC1AttemptPremiseError)
        if type(self.accepted_state_advanced) is not bool:
            raise TypeError("accepted_state_advanced must be a built-in bool")
        if type(self.executable_retry_cursor) is not bool:
            raise TypeError("executable_retry_cursor must be a built-in bool")
        if type(self.physical_state_preserved) is not bool:
            raise TypeError("physical_state_preserved must be a built-in bool")
        if self.kind in BLOB_UNCHANGED_KINDS:
            if self.accepted_state_advanced or self.executable_retry_cursor:
                _fail(
                    f"{self.kind} cannot advance state or advertise an executable cursor",
                    RSRC1AttemptPremiseError,
                )
            if type(self.terminal_evidence) is not dict:
                _fail(f"{self.kind} omitted terminal evidence", RSRC1AttemptPremiseError)
            if self.terminal_evidence.get("kind") != self.kind:
                _fail(f"{self.kind} evidence kind differs", RSRC1AttemptPremiseError)
            if self.terminal_evidence.get("executable_retry_cursor") is not False:
                _fail(
                    f"{self.kind} evidence advertised an executable retry cursor",
                    RSRC1AttemptPremiseError,
                )
        if self.kind == "accepted_fine" and self.accepted_state_advanced is not True:
            _fail("accepted_fine requires accepted-state advance", RSRC1AttemptPremiseError)
        if self.kind in TERMINAL_KINDS and self.executable_retry_cursor:
            _fail("typed terminal cannot remain executable", RSRC1AttemptPremiseError)
        object.__setattr__(
            self,
            "terminal_evidence",
            None if self.terminal_evidence is None else dict(self.terminal_evidence),
        )


def _overlay_exhausted(cursor: HLT17IMP1Cursor) -> bool:
    overlay = cursor.overlay
    if overlay is None:
        return False
    if overlay.exhausted is True or overlay.next_plan is None:
        return True
    cap = SOURCE_RETRY_CAP if overlay.owner == SOURCE_OWNER else CFL_RETRY_CAP
    return overlay.owner_current_count > cap


def classify_attempt(
    *,
    member_key: str,
    predecessor_bundle: HLT17GenerationBundle,
    successor_bundle: HLT17GenerationBundle,
    accepted_state_advanced: bool,
    resource_stop: Mapping[str, object] | None = None,
    premise_stop: Mapping[str, object] | None = None,
    temporal_exhaustion_evidence: Mapping[str, object] | None = None,
) -> AttemptClassification:
    """Independently reduce one immutable predecessor/successor relation."""

    if member_key not in MEMBER_KEYS:
        _fail("member_key is outside the live cohort", RSRC1AttemptPremiseError)
    if type(accepted_state_advanced) is not bool:
        raise TypeError("accepted_state_advanced must be a built-in bool")
    identical = _bundles_identical(predecessor_bundle, successor_bundle)
    encoding_changed = _encoding_changed(predecessor_bundle, successor_bundle)
    predecessor_cursor = _cursor_from_bundle(predecessor_bundle)
    successor_cursor = _cursor_from_bundle(successor_bundle)
    physical_preserved = (
        predecessor_cursor.physical_state_sha256
        == successor_cursor.physical_state_sha256
        and _same_bits(predecessor_cursor.accepted_time, successor_cursor.accepted_time)
    )
    if resource_stop is not None:
        if accepted_state_advanced:
            _fail(
                "resource stop cannot advertise accepted-state advance",
                RSRC1AttemptPremiseError,
            )
        if not identical:
            _fail(
                "resource stop must return the unchanged accepted bundle",
                RSRC1AttemptPremiseError,
            )
        extras = {
            key: value
            for key, value in dict(resource_stop).items()
            if key not in {"kind", "reason", "executable_retry_cursor", "spends_retry"}
        }
        evidence = {
            **extras,
            "kind": "resource_exhausted",
            "reason": str(resource_stop.get("reason") or "resource_ceiling"),
            "executable_retry_cursor": False,
            "spends_retry": False,
        }
        return AttemptClassification(
            kind="resource_exhausted",
            member_key=member_key,
            successor_bundle=predecessor_bundle,
            terminal_evidence=evidence,
            accepted_state_advanced=False,
            executable_retry_cursor=False,
            physical_state_preserved=True,
            structural_class="typed_resource_stop",
        )
    if premise_stop is not None:
        if accepted_state_advanced:
            _fail(
                "premise stop cannot advertise accepted-state advance",
                RSRC1AttemptPremiseError,
            )
        if not identical:
            _fail(
                "premise stop must return the unchanged accepted bundle",
                RSRC1AttemptPremiseError,
            )
        evidence = {
            "kind": "invalid_premise",
            "reason": str(premise_stop.get("reason") or "invalid_premise"),
            "executable_retry_cursor": False,
            "spends_retry": False,
        }
        return AttemptClassification(
            kind="invalid_premise",
            member_key=member_key,
            successor_bundle=predecessor_bundle,
            terminal_evidence=evidence,
            accepted_state_advanced=False,
            executable_retry_cursor=False,
            physical_state_preserved=True,
            structural_class="typed_premise_stop",
        )
    if accepted_state_advanced:
        if identical or not encoding_changed or physical_preserved:
            _fail(
                "accepted advance did not change the accepted physical encoding",
                RSRC1AttemptPremiseError,
            )
        if successor_cursor.mode != FRESH_READY:
            _fail("accepted fine did not return a fresh cursor", RSRC1AttemptPremiseError)
        return AttemptClassification(
            kind="accepted_fine",
            member_key=member_key,
            successor_bundle=successor_bundle,
            terminal_evidence=None,
            accepted_state_advanced=True,
            executable_retry_cursor=False,
            physical_state_preserved=False,
            structural_class="accepted_physical_advance",
        )
    if encoding_changed:
        _fail(
            "non-fine attempt changed accepted encoding without an advance flag",
            RSRC1AttemptUnexpectedError,
        )
    if identical:
        if temporal_exhaustion_evidence is None:
            _fail(
                "unchanged bundle is not a typed terminal without evidence",
                RSRC1AttemptUnexpectedError,
            )
        if accepted_state_advanced:
            _fail(
                "temporal exhaustion cannot advertise accepted-state advance",
                RSRC1AttemptPremiseError,
            )
        evidence = dict(temporal_exhaustion_evidence)
        evidence["kind"] = "temporal_retry_exhausted"
        evidence["executable_retry_cursor"] = False
        if evidence.get("reason") not in TEMPORAL_EXHAUSTION_REASONS:
            _fail("temporal exhaustion reason differs", RSRC1AttemptPremiseError)
        return AttemptClassification(
            kind="temporal_retry_exhausted",
            member_key=member_key,
            successor_bundle=predecessor_bundle,
            terminal_evidence=evidence,
            accepted_state_advanced=False,
            executable_retry_cursor=False,
            physical_state_preserved=True,
            structural_class="temporal_exhausted_unchanged",
        )
    overlay = successor_cursor.overlay
    if overlay is not None:
        if not physical_preserved:
            _fail(
                "source/CFL retry changed accepted physical state",
                RSRC1AttemptUnexpectedError,
            )
        exhausted = _overlay_exhausted(successor_cursor)
        if overlay.owner == SOURCE_OWNER:
            kind = "source_retry_exhausted" if exhausted else "source_retry"
            structural = "source_overlay_exhausted" if exhausted else "source_overlay_retry"
        elif overlay.owner == CFL_OWNER:
            kind = "cfl_retry_exhausted" if exhausted else "cfl_retry"
            structural = "cfl_overlay_exhausted" if exhausted else "cfl_overlay_retry"
        else:
            _fail("overlay owner is not source or CFL", RSRC1AttemptUnexpectedError)
            raise AssertionError("unreachable")
        evidence = None
        if kind in TERMINAL_KINDS:
            evidence = {
                "kind": kind,
                "exhausted": True,
                "executable_retry_cursor": False,
                "next_plan": None,
            }
        return AttemptClassification(
            kind=kind,
            member_key=member_key,
            successor_bundle=successor_bundle,
            terminal_evidence=evidence,
            accepted_state_advanced=False,
            executable_retry_cursor=(
                kind in CURSOR_CHANGED_RETRY_KINDS and kind not in TERMINAL_KINDS
            ),
            physical_state_preserved=True,
            structural_class=structural,
        )
    if (
        successor_cursor.mode == RETRY_PENDING
        and successor_cursor.temporal is not None
        and successor_cursor.next_attempt_plan() is not None
    ):
        if not physical_preserved:
            _fail("temporal retry changed accepted physical state", RSRC1AttemptUnexpectedError)
        return AttemptClassification(
            kind="temporal_retry",
            member_key=member_key,
            successor_bundle=successor_bundle,
            terminal_evidence=None,
            accepted_state_advanced=False,
            executable_retry_cursor=True,
            physical_state_preserved=True,
            structural_class="temporal_pending_successor",
        )
    _fail(
        "cursor changed without an executable retry plan; temporal exhaustion "
        "must not be absorbed into an executable cursor",
        RSRC1AttemptUnexpectedError,
    )
    raise AssertionError("unreachable")


def restore_checkpoint_from_bundle(
    checkpoint: HLT17InMemoryCheckpoint,
    bundle: HLT17GenerationBundle,
) -> None:
    encoding = HLT17MemberEncoding(bundle.descriptor, bundle.payload)
    cursor = _cursor_from_bundle(bundle)
    checkpoint.restore(
        encoding,
        cursor,
        last_accepted_executed_plan=encoding.executed_plan,
    )


def _resource_evidence(reason: str, **fields: object) -> dict[str, object]:
    return {"reason": reason, **fields}


def execute_scheduled_attempt(
    checkpoint: HLT17InMemoryCheckpoint,
    *,
    requested_cap: float | None = None,
) -> AttemptClassification:
    """Execute one prepare/admit/adopt path without importing a store."""

    if not isinstance(checkpoint, HLT17InMemoryCheckpoint):
        raise TypeError("checkpoint must be HLT17InMemoryCheckpoint")
    checkpoint.agree()
    predecessor = checkpoint.generation_bundle()
    member_key = checkpoint.member.key
    cursor = checkpoint.cursor
    if checkpoint.at_event_target():
        _fail("parked member cannot be advanced", RSRC1AttemptPremiseError)
    try:
        if cursor.public_fresh_permitted():
            cap = fresh_requested_cap(checkpoint) if requested_cap is None else float(
                requested_cap
            )
            prepared = checkpoint.prepare(requested_cap=cap)
        else:
            if pending_plan_from_cursor(cursor) is None:
                _fail("pending cursor omitted its next plan", RSRC1AttemptPremiseError)
            prepared = checkpoint.prepare(requested_cap=requested_cap)
    except HLT17OwnerRetryExhausted as error:
        return classify_attempt(
            member_key=member_key,
            predecessor_bundle=predecessor,
            successor_bundle=predecessor,
            accepted_state_advanced=False,
            premise_stop={
                "reason": (
                    f"persisted_{error.owner}_exhaustion_without_terminal:"
                    f"{error.reason}"
                )
            },
        )
    except HLT17InvalidPremise as error:
        return classify_attempt(
            member_key=member_key,
            predecessor_bundle=predecessor,
            successor_bundle=predecessor,
            accepted_state_advanced=False,
            premise_stop={"reason": str(error) or "invalid_premise"},
        )
    except C1R1RuntimeResourceStop as error:
        return classify_attempt(
            member_key=member_key,
            predecessor_bundle=predecessor,
            successor_bundle=predecessor,
            accepted_state_advanced=False,
            resource_stop=_resource_evidence(str(error) or "c1r1_runtime_resource"),
        )
    except (C1R1CoordinateLatticeStop, C1R1RefinementPathStop) as error:
        return classify_attempt(
            member_key=member_key,
            predecessor_bundle=predecessor,
            successor_bundle=predecessor,
            accepted_state_advanced=False,
            premise_stop={"reason": str(error) or type(error).__name__},
        )
    if isinstance(prepared, HLT17CheckpointResult):
        _fail("prepare returned a checkpoint result", RSRC1AttemptPremiseError)
    if not isinstance(prepared, HLT17BridgeResult):
        raise TypeError("prepare must return HLT17BridgeResult")
    if prepared.accepted_state_advanced:
        _fail("prepare advanced accepted state", RSRC1AttemptUnexpectedError)
    overlays = {
        "source_retry_required",
        "cfl_retry_required",
        "source_retry_exhausted",
        "cfl_retry_exhausted",
    }
    if prepared.disposition in overlays | {"temporal_retry_required"}:
        absorbed = checkpoint.absorb_nonfine(prepared)
        return classify_attempt(
            member_key=member_key,
            predecessor_bundle=predecessor,
            successor_bundle=absorbed.bundle,
            accepted_state_advanced=False,
        )
    if prepared.disposition == "temporal_retry_exhausted":
        evidence = dict(prepared.evidence)
        evidence["kind"] = "temporal_retry_exhausted"
        evidence["executable_retry_cursor"] = False
        return classify_attempt(
            member_key=member_key,
            predecessor_bundle=predecessor,
            successor_bundle=predecessor,
            accepted_state_advanced=False,
            temporal_exhaustion_evidence=evidence,
        )
    if prepared.disposition == "invalid_premise":
        return classify_attempt(
            member_key=member_key,
            predecessor_bundle=predecessor,
            successor_bundle=predecessor,
            accepted_state_advanced=False,
            premise_stop={"reason": str(prepared.evidence.get("reason") or "invalid_premise")},
        )
    if prepared.prepared is None:
        _fail("prepare omitted an admissible family", RSRC1AttemptPremiseError)
    try:
        admitted = checkpoint.require_admission(prepared.prepared)
    except C1R1RuntimeResourceStop as error:
        return classify_attempt(
            member_key=member_key,
            predecessor_bundle=predecessor,
            successor_bundle=predecessor,
            accepted_state_advanced=False,
            resource_stop=_resource_evidence(str(error) or "c1r1_runtime_resource"),
        )
    except C1R1TemporalRetryExhausted as error:
        evidence = {
            "kind": "temporal_retry_exhausted",
            "reason": error.reason,
            "updated_imp1_ledger": imp1_checkpoint_extension(error.updated_ledger),
            "updated_ledger_sha256": ledger_sha256(error.updated_ledger),
            "executable_retry_cursor": False,
        }
        return classify_attempt(
            member_key=member_key,
            predecessor_bundle=predecessor,
            successor_bundle=predecessor,
            accepted_state_advanced=False,
            temporal_exhaustion_evidence=evidence,
        )
    except (HLT17InvalidPremise, HLT17OwnerRetryExhausted) as error:
        return classify_attempt(
            member_key=member_key,
            predecessor_bundle=predecessor,
            successor_bundle=predecessor,
            accepted_state_advanced=False,
            premise_stop={"reason": str(error) or "invalid_premise"},
        )
    if isinstance(admitted, HLT17CheckpointResult):
        _fail("admission returned a checkpoint result", RSRC1AttemptPremiseError)
    if admitted.disposition in overlays | {"temporal_retry_required"}:
        absorbed = checkpoint.absorb_nonfine(admitted)
        return classify_attempt(
            member_key=member_key,
            predecessor_bundle=predecessor,
            successor_bundle=absorbed.bundle,
            accepted_state_advanced=False,
        )
    if admitted.disposition == "temporal_retry_exhausted":
        evidence = dict(admitted.evidence)
        evidence["kind"] = "temporal_retry_exhausted"
        evidence["executable_retry_cursor"] = False
        reason = str(evidence.get("reason") or "maximum_temporal_retries")
        evidence["reason"] = (
            reason if reason in TEMPORAL_EXHAUSTION_REASONS else "maximum_temporal_retries"
        )
        return classify_attempt(
            member_key=member_key,
            predecessor_bundle=predecessor,
            successor_bundle=predecessor,
            accepted_state_advanced=False,
            temporal_exhaustion_evidence=evidence,
        )
    if admitted.disposition == "invalid_premise":
        return classify_attempt(
            member_key=member_key,
            predecessor_bundle=predecessor,
            successor_bundle=predecessor,
            accepted_state_advanced=False,
            premise_stop={"reason": str(admitted.evidence.get("reason") or "invalid_premise")},
        )
    adopted = checkpoint.adopt_fine(admitted.prepared)
    if adopted.disposition != ACCEPTED_FINE or adopted.accepted_state_advanced is not True:
        _fail("fine adoption did not advance accepted state", RSRC1AttemptUnexpectedError)
    return classify_attempt(
        member_key=member_key,
        predecessor_bundle=predecessor,
        successor_bundle=adopted.bundle,
        accepted_state_advanced=True,
    )


__all__ = [
    "ARTIFACT_ID",
    "AttemptClassification",
    "RSRC1AttemptError",
    "RSRC1AttemptPremiseError",
    "RSRC1AttemptUnexpectedError",
    "classify_attempt",
    "execute_scheduled_attempt",
    "fresh_requested_cap",
    "pending_plan_from_cursor",
    "restore_checkpoint_from_bundle",
]
