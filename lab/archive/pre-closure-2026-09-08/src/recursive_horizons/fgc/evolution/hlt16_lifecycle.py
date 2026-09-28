"""Pure, fail-closed PROTO18 cursor transitions for HLT16.

This owner does not read stores, advance PDEs, or select physical branches.
PROTO15 remains the authority for cursor, TDG6, and TDG7 representations.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping

from . import proto15_runtime as p15
from . import tdg5_stage_complete_refinement_runtime as tdg5
from . import tdg6_temporal_admission_runtime as tdg6


PROTO17_ARTIFACT_ID = "FGC-2-SF1-PROTO17"
PROTO18_ARTIFACT_ID = "FGC-2-SF1-PROTO18"
GEN0_EVENT_INDEX = 23
GEN0_TIME_HEX = float(23 / 16).hex()
FIRST_EVENT_TARGET_HEX = float(3 / 2).hex()


class HLT16LifecycleError(ValueError):
    """A prospective cursor transition violates the frozen lifecycle."""


@dataclass(frozen=True, slots=True)
class RetryExhaustion:
    """Evidence for a TDG6 numerical stop; it is never physical evidence."""

    reason: str
    member_key: str
    accepted_state_sha256: str
    accepted_time: Mapping[str, str]
    updated_ledger: p15.TDG6TemporalLedger
    durable_rejection_record_sha256: str
    predecessor_cursor_sha256: str


def _wrap(function: Any, *args: Any, **kwargs: Any) -> Any:
    try:
        return function(*args, **kwargs)
    except (TypeError, ValueError, KeyError, AttributeError) as error:
        raise HLT16LifecycleError(str(error)) from error


def _sha(name: str, value: object, *, non_root: bool = False) -> str:
    digest = _wrap(p15._require_sha, name, value)
    if non_root and digest == p15.ROOT_SHA256:
        raise HLT16LifecycleError(f"{name} must not be the root digest")
    return digest


def _time(value: object) -> dict[str, str]:
    return _wrap(p15._time_identity, value)


def _checked_cursor(cursor: p15.Proto15Cursor) -> p15.Proto15Cursor:
    if not isinstance(cursor, p15.Proto15Cursor):
        raise HLT16LifecycleError("cursor must be a Proto15Cursor")
    _wrap(cursor.validate)
    return cursor


def _checked_ledger(ledger: p15.TDG6TemporalLedger) -> p15.TDG6TemporalLedger:
    """Round-trip through PROTO15's strict persisted representation."""
    if not isinstance(ledger, p15.TDG6TemporalLedger):
        raise HLT16LifecycleError("ledger must be TDG6TemporalLedger")
    mapping = _wrap(p15._ledger_mapping, ledger)
    rebuilt = _wrap(p15._ledger_from_mapping, mapping)
    if _wrap(p15._ledger_mapping, rebuilt) != mapping:
        raise HLT16LifecycleError("TDG6 ledger does not round-trip canonically")
    return rebuilt


_CUBIC_ENVELOPE_FIELDS = frozenset(
    {
        "polynomial_count",
        "endpoint_candidate_count",
        "real_interior_root_count",
        "ambiguous_discriminant_count",
        "raw_candidate_maximum",
        "coefficient_construction_debit",
        "outward_arithmetic_debit",
        "bernstein_certification_slack",
        "certified_continuous_upper_bound",
    }
)
_MAGNITUDE_INTERVAL_FIELDS = frozenset(
    {
        "lower_bound",
        "upper_bound",
        "subinterval_count",
        "owned_row_count",
        "envelope",
    }
)


def _exact_mapping(
    name: str, value: object, expected_fields: frozenset[str]
) -> dict[str, object]:
    if not isinstance(value, Mapping):
        raise HLT16LifecycleError(f"{name} must be a mapping")
    mapping = dict(value)
    if frozenset(mapping) != expected_fields:
        raise HLT16LifecycleError(f"{name} fields differ")
    return mapping


def _decode_magnitude_interval(
    name: str, value: object
) -> tdg6.TDG6Binary64MagnitudeInterval:
    """Strictly inflate one persisted TDG6 interval and nested envelope."""
    interval = _exact_mapping(name, value, _MAGNITUDE_INTERVAL_FIELDS)
    envelope_mapping = _exact_mapping(
        f"{name} envelope", interval["envelope"], _CUBIC_ENVELOPE_FIELDS
    )
    envelope = _wrap(tdg5.Binary64CubicEnvelope, **envelope_mapping)
    return _wrap(
        tdg6.TDG6Binary64MagnitudeInterval,
        **{**interval, "envelope": envelope},
    )


def _ledger_hash(ledger: p15.TDG6TemporalLedger) -> str:
    return _wrap(p15._hash, _wrap(p15._ledger_mapping, _checked_ledger(ledger)))


def _require_cursor_ledger(cursor: p15.Proto15Cursor, ledger: p15.TDG6TemporalLedger) -> None:
    if cursor.payload["TDG6_ledger_sha256"] != _ledger_hash(ledger):
        raise HLT16LifecycleError("cursor TDG6 ledger hash differs from supplied ledger")
    if cursor.payload["accepted_boundary_time"]["binary64_hex"] != ledger.last_accepted_time.hex():
        raise HLT16LifecycleError("cursor accepted time differs from TDG6 ledger")


def _require_proto18(cursor: p15.Proto15Cursor, ledger: p15.TDG6TemporalLedger) -> None:
    _checked_cursor(cursor)
    if cursor.payload["protocol_artifact_id"] != PROTO18_ARTIFACT_ID:
        raise HLT16LifecycleError("transition requires a PROTO18 cursor")
    _require_cursor_ledger(cursor, ledger)


def _next_identity(cursor: p15.Proto15Cursor) -> tuple[int, int, str]:
    generation = int(cursor.payload["cursor_generation"]) + 1
    serial = int(cursor.payload["attempt_serial"]) + 1
    return (
        generation,
        serial,
        f"{cursor.payload['campaign_id']}:{cursor.payload['member_key']}:"
        f"{cursor.payload['committed_common_event_index']}:{serial}",
    )


def _new_cursor(payload: Mapping[str, object]) -> p15.Proto15Cursor:
    return _wrap(p15._cursor, payload)


def _observed_cursor(payload: Mapping[str, object]) -> p15.Proto15Cursor:
    """Validate an observed cursor without normalizing its supplied digest."""
    if not isinstance(payload, Mapping) or "cursor_chain_sha256" not in payload:
        raise HLT16LifecycleError("observed cursor digest is absent")
    return _checked_cursor(_wrap(p15.Proto15Cursor, dict(payload)))


def canonical_mapping(cursor: p15.Proto15Cursor) -> dict[str, object]:
    return dict(_checked_cursor(cursor).payload)


def canonical_hash(cursor: p15.Proto15Cursor) -> str:
    return _checked_cursor(cursor).sha256


def bridge_gen0(
    gen0: Mapping[str, object] | p15.Proto15Cursor,
    ledger: p15.TDG6TemporalLedger,
    descriptor_hash: str,
    *,
    successor_campaign_id: str | None = None,
) -> p15.Proto15Cursor:
    """Bridge the one authenticated PROTO17 GEN0 state into PROTO18.

    GEN0 contains source-array identities, not restartable evolved payloads.
    This bridge is therefore valid only at the frozen initial event/time with a
    zero TDG6 progression ledger and root journal.
    """
    old = _checked_cursor(
        gen0 if isinstance(gen0, p15.Proto15Cursor) else _observed_cursor(gen0)
    )
    checked = _checked_ledger(ledger)
    descriptor = _sha("authenticated GEN0 descriptor", descriptor_hash, non_root=True)
    p = old.payload
    if (
        p["protocol_artifact_id"] != PROTO17_ARTIFACT_ID
        or p["mode"] != "FRESH_READY"
        or p["committed_common_event_index"] != GEN0_EVENT_INDEX
        or p["cursor_generation"] != 0
        or p["attempt_serial"] != 0
        or p["journal_tip_sha256"] != p15.ROOT_SHA256
        or p["cursor_chain_parent_sha256"] != p15.ROOT_SHA256
        or p["accepted_boundary_time"]["binary64_hex"] != GEN0_TIME_HEX
        or p["event_target_time"]["binary64_hex"] != FIRST_EVENT_TARGET_HEX
    ):
        raise HLT16LifecycleError("GEN0 cursor differs from the frozen first-event boundary")
    _require_cursor_ledger(old, checked)
    if (
        checked.accepted_macro_step_count != 0
        or checked.cumulative_temporal_retry_count != 0
        or checked.current_macro_step_temporal_retry_count != 0
        or checked.last_accepted_macro_step_temporal_retry_count != 0
        or checked.serialized_temporal_rejections
        or any(value != 0.0 for value in checked.accumulated_debit_vector)
    ):
        raise HLT16LifecycleError("GEN0 TDG6 ledger is not the authenticated zero ledger")
    campaign_id = p["campaign_id"] if successor_campaign_id is None else successor_campaign_id
    if (
        not isinstance(campaign_id, str)
        or not campaign_id
        or ":" in campaign_id
    ):
        raise HLT16LifecycleError("successor campaign identity differs")
    successor = dict(p)
    successor.update(
        {
            "protocol_artifact_id": PROTO18_ARTIFACT_ID,
            "campaign_id": campaign_id,
            "accepted_state_sha256": descriptor,
            "TDG6_ledger_sha256": _ledger_hash(checked),
            "cursor_generation": 0,
            "attempt_serial": 0,
            "attempt_id": (
                f"{campaign_id}:{p['member_key']}:"
                f"{p['committed_common_event_index']}:0"
            ),
            "mode": "FRESH_READY",
            "retry_successor_payload_or_none": None,
            "journal_tip_sha256": p15.ROOT_SHA256,
            "cursor_chain_parent_sha256": old.sha256,
        }
    )
    successor.pop("cursor_chain_sha256", None)
    return _new_cursor(successor)


def validate_retry_evidence(
    cursor: p15.Proto15Cursor,
    prior_ledger: p15.TDG6TemporalLedger,
    evidence_mapping: Mapping[str, object],
    *,
    evolution_state_sha256: str,
) -> dict[str, object]:
    """Validate TDG6 evidence against the actual ``u,p,q`` identity.

    PROTO18 deliberately uses ``cursor.accepted_state_sha256`` for the wider
    persisted descriptor identity.  TDG6's ``initial_state_sha256`` instead
    identifies only the evolution fields ``u,p,q``.  Both are truthful and
    must remain distinct; callers bind the latter independently to the arrays
    owned by the descriptor.
    """

    prior = _checked_ledger(prior_ledger)
    _require_proto18(cursor, prior)
    physical = _sha("evolution u,p,q state", evolution_state_sha256, non_root=True)
    return dict(
        _wrap(
            p15._validate_retry_evidence_mapping,
            dict(evidence_mapping),
            method=cursor.payload["method"],
            state_sha256=physical,
            accepted_time=cursor.payload["accepted_boundary_time"],
            step_index=cursor.payload["previous_step_index"],
            transaction_serial=cursor.payload["previous_transaction_serial"],
            prior_ledger=prior,
        )
    )


def validate_retry_payload(
    cursor: p15.Proto15Cursor,
    ledger: p15.TDG6TemporalLedger,
    *,
    evolution_state_sha256: str,
) -> None:
    """Validate one PROTO18 pending cursor across both state identities.

    The legacy PROTO15 validator correctly assumes its cursor state address and
    TDG6 field-state hash are the same object.  HLT16 persists a broader
    descriptor, so this adapter first proves the real descriptor linkage, then
    runs the complete legacy relation check on a validation-only shadow whose
    state identity is the independently reconstructed ``u,p,q`` hash.  The
    actual cursor and durable evidence are never rewritten.
    """

    checked = _checked_ledger(ledger)
    _require_proto18(cursor, checked)
    physical = _sha("evolution u,p,q state", evolution_state_sha256, non_root=True)
    pending = cursor.payload["retry_successor_payload_or_none"]
    if not isinstance(pending, Mapping):
        raise HLT16LifecycleError("pending cursor omits retry payload")
    if pending.get("accepted_state_sha256") != cursor.payload["accepted_state_sha256"]:
        raise HLT16LifecycleError("pending cursor descriptor identity differs")
    shadow_payload = dict(cursor.payload)
    shadow_pending = dict(pending)
    shadow_payload["accepted_state_sha256"] = physical
    shadow_pending["accepted_state_sha256"] = physical
    shadow_payload["retry_successor_payload_or_none"] = shadow_pending
    shadow = p15.Proto15Cursor(shadow_payload)
    _wrap(p15._validate_retry_payload, shadow, checked)


def _validated_retry_outcome(
    cursor: p15.Proto15Cursor,
    prior: p15.TDG6TemporalLedger,
    outcome: object,
    *,
    evolution_state_sha256: str,
) -> tuple[dict[str, object], p15.TDG6TemporalLedger]:
    if not isinstance(outcome, (p15.TDG6TemporalRetryRequired, p15.TDG6TemporalRetryExhausted)):
        raise HLT16LifecycleError("outcome must be a typed TDG6 retry result")
    if not isinstance(outcome.evidence, p15.TDG6TemporalRetryEvidence):
        raise HLT16LifecycleError("typed retry outcome omits canonical evidence")
    evidence = validate_retry_evidence(
        cursor,
        prior,
        _wrap(p15._json_safe, outcome.evidence.as_mapping()),
        evolution_state_sha256=evolution_state_sha256,
    )
    updated = _checked_ledger(_wrap(p15._updated_retry_ledger, prior, evidence))
    if _wrap(p15._ledger_mapping, outcome.updated_ledger) != _wrap(p15._ledger_mapping, updated):
        raise HLT16LifecycleError("typed retry outcome ledger differs from its evidence")
    exhausted = (
        int(evidence["retry_count_for_current_macro_step"]) > p15.TDG6_MAXIMUM_TEMPORAL_RETRIES_PER_STEP
        or float(evidence["retry_step_size"]) < float(p15.TDG6_MINIMUM_MACRO_STEP)
    )
    if isinstance(outcome, p15.TDG6TemporalRetryRequired) and exhausted:
        raise HLT16LifecycleError("retry-required outcome is already exhausted")
    if isinstance(outcome, p15.TDG6TemporalRetryExhausted):
        expected = "maximum_temporal_retries" if int(evidence["retry_count_for_current_macro_step"]) > p15.TDG6_MAXIMUM_TEMPORAL_RETRIES_PER_STEP else "minimum_macro_step"
        if not exhausted or outcome.reason != expected:
            raise HLT16LifecycleError("retry exhaustion reason differs from its evidence")
    return dict(evidence), updated


def typed_retry_from_evidence(
    cursor: p15.Proto15Cursor,
    prior_ledger: p15.TDG6TemporalLedger,
    evidence_mapping: Mapping[str, object],
    *,
    evolution_state_sha256: str,
) -> p15.TDG6TemporalRetryRequired | p15.TDG6TemporalRetryExhausted:
    """Rebuild the one typed TDG6 result proved by persisted rejection evidence.

    This is a recovery-only decoder: it validates the frozen evidence against
    the accepted cursor and ledger before rebuilding a typed outcome.  It does
    not evaluate a PDE path or create a new numerical choice.
    """
    prior = _checked_ledger(prior_ledger)
    _require_proto18(cursor, prior)
    try:
        checked = validate_retry_evidence(
            cursor,
            prior,
            dict(evidence_mapping),
            evolution_state_sha256=evolution_state_sha256,
        )
        admissions = tuple(
            tdg6.TDG6RuntimeChannelAdmission(
                **{
                    **dict(item),
                    "outer_difference": _decode_magnitude_interval(
                        "outer_difference", item["outer_difference"]
                    ),
                    "finest_difference": _decode_magnitude_interval(
                        "finest_difference", item["finest_difference"]
                    ),
                }
            )
            for item in checked["channel_admissions"]
        )
        evidence = tdg6.TDG6TemporalRetryEvidence(
            **{
                **dict(checked),
                "failed_channels": tuple(checked["failed_channels"]),
                "channel_admissions": admissions,
            }
        )
        updated = p15._updated_retry_ledger(prior, checked)
    except (TypeError, ValueError, KeyError, AttributeError) as exc:
        raise HLT16LifecycleError("persisted TDG6 rejection evidence differs") from exc
    exhausted = (
        int(checked["retry_count_for_current_macro_step"])
        > p15.TDG6_MAXIMUM_TEMPORAL_RETRIES_PER_STEP
        or float(checked["retry_step_size"]) < float(p15.TDG6_MINIMUM_MACRO_STEP)
    )
    if exhausted:
        reason = "maximum_temporal_retries" if int(checked["retry_count_for_current_macro_step"]) > p15.TDG6_MAXIMUM_TEMPORAL_RETRIES_PER_STEP else "minimum_macro_step"
        return tdg6.TDG6TemporalRetryExhausted(reason=reason, evidence=evidence, updated_ledger=updated)
    return tdg6.TDG6TemporalRetryRequired(evidence, updated)


def close_retry(
    cursor: p15.Proto15Cursor,
    prior_ledger: p15.TDG6TemporalLedger,
    outcome: p15.TDG6TemporalRetryRequired,
    predecessor_plan: Mapping[str, object],
    successor_plan: Mapping[str, object],
    durable_rejection_record_sha256: str,
    predecessor_journal_tip_sha256: str,
    *,
    evolution_state_sha256: str,
) -> p15.Proto15Cursor:
    """Close one durable TDG6 rejection into exactly one RETRY_PENDING cursor."""
    prior = _checked_ledger(prior_ledger)
    _require_proto18(cursor, prior)
    if cursor.mode not in {"FRESH_READY", "RETRY_PENDING"}:
        raise HLT16LifecycleError("TDG6 rejection requires an executable cursor")
    # An accepted-fine cursor embeds the journal *parent* of its own accepted
    # record to avoid a cursor/record self-hash.  The complete campaign
    # checkpoint therefore owns the current journal tip.  Validate the caller's
    # checkpoint tip as a digest, but do not incorrectly equate it with the
    # cursor's older embedded ancestor.
    predecessor_journal_tip_sha256 = _sha(
        "predecessor checkpoint journal tip", predecessor_journal_tip_sha256
    )
    rejection = _sha("durable rejection record", durable_rejection_record_sha256, non_root=True)
    evidence, updated = _validated_retry_outcome(
        cursor,
        prior,
        outcome,
        evolution_state_sha256=evolution_state_sha256,
    )
    if not isinstance(outcome, p15.TDG6TemporalRetryRequired):
        raise HLT16LifecycleError("retry-pending closure requires TDG6TemporalRetryRequired")
    predecessor = _wrap(p15._plan_mapping, predecessor_plan)
    successor_plan_mapping = _wrap(p15._plan_mapping, successor_plan)
    if cursor.mode == "RETRY_PENDING":
        validate_retry_payload(
            cursor,
            prior,
            evolution_state_sha256=evolution_state_sha256,
        )
        pending = cursor.payload["retry_successor_payload_or_none"]
        assert isinstance(pending, Mapping)
        durable_predecessor = _wrap(p15._plan_mapping, pending["successor_plan"])
        if predecessor != durable_predecessor:
            raise HLT16LifecycleError(
                "pending TDG6 rejection did not execute the exact durable successor plan"
            )
    if (
        predecessor["current_hex"] != cursor.payload["accepted_boundary_time"]["binary64_hex"]
        or predecessor["event_target_hex"] != cursor.payload["event_target_time"]["binary64_hex"]
        or not p15._same_binary64(float.fromhex(predecessor["macro_width_hex"]), evidence["attempted_macro_step_size"])
    ):
        raise HLT16LifecycleError("TDG6 rejection is not bound to the current exact predecessor plan")
    expected = _wrap(
        p15._plan_mapping,
        p15.plan_forward_proto14_subdivision(
            float(evidence["initial_time"]),
            float.fromhex(str(cursor.payload["event_target_time"]["binary64_hex"])),
            float(evidence["retry_step_size"]),
            minimum_width=float(p15.TDG6_MINIMUM_MACRO_STEP),
        ),
    )
    if successor_plan_mapping != expected:
        raise HLT16LifecycleError("TDG6 retry successor is not the exact TDG7 half-step plan")
    generation, serial, attempt_id = _next_identity(cursor)
    pending: dict[str, object] = {
        "predecessor_plan": predecessor,
        "successor_plan": successor_plan_mapping,
        "TDG6_rejection_evidence": evidence,
        "complete_rejection_prefix": list(updated.serialized_temporal_rejections),
        "durable_rejection_record_sha256": rejection,
        "predecessor_journal_tip_sha256": predecessor_journal_tip_sha256,
        "predecessor_attempt_id": cursor.payload["attempt_id"],
        "successor_attempt_id": attempt_id,
        "member_key": cursor.payload["member_key"],
        "method": cursor.payload["method"],
        "point_count": cursor.payload["point_count"],
        "accepted_state_sha256": cursor.payload["accepted_state_sha256"],
        "accepted_boundary_time": dict(cursor.payload["accepted_boundary_time"]),
        "previous_step_index": cursor.payload["previous_step_index"],
        "previous_transaction_serial": cursor.payload["previous_transaction_serial"],
        "retry_count": updated.current_macro_step_temporal_retry_count,
        "half_cap": outcome.retry_step_size.hex(),
        "minimum_macro_step": float(p15.TDG6_MINIMUM_MACRO_STEP).hex(),
    }
    payload = dict(cursor.payload)
    payload.update(
        {
            "TDG6_ledger_sha256": _ledger_hash(updated),
            "cursor_generation": generation,
            "attempt_serial": serial,
            "attempt_id": attempt_id,
            "mode": "RETRY_PENDING",
            "retry_successor_payload_or_none": pending,
            "journal_tip_sha256": rejection,
            "cursor_chain_parent_sha256": cursor.sha256,
        }
    )
    payload.pop("cursor_chain_sha256", None)
    successor_cursor = _new_cursor(payload)
    validate_retry_payload(
        successor_cursor,
        updated,
        evolution_state_sha256=evolution_state_sha256,
    )
    return successor_cursor


def retry_exhaustion(
    cursor: p15.Proto15Cursor,
    prior_ledger: p15.TDG6TemporalLedger,
    outcome: p15.TDG6TemporalRetryExhausted,
    durable_rejection_record_sha256: str,
    predecessor_journal_tip_sha256: str,
    *,
    evolution_state_sha256: str,
) -> RetryExhaustion:
    """Validate a TDG6 numerical stop without creating a successor cursor.

    A stop can occur on the original fresh proposal or on the exact TDG7
    successor of one already durable rejection.  The latter is deliberately
    *not* a second chance to invent a pending history: the pending cursor must
    first validate against PROTO15's persisted rejection payload and the new
    exhaustion must be for that successor plan.  In both cases this function
    only returns terminal evidence; it never advances the accepted state or
    manufactures another cursor.
    """
    prior = _checked_ledger(prior_ledger)
    _require_proto18(cursor, prior)
    if cursor.mode not in {"FRESH_READY", "RETRY_PENDING"}:
        raise HLT16LifecycleError("retry exhaustion requires an executable cursor")
    predecessor_journal_tip_sha256 = _sha(
        "predecessor checkpoint journal tip", predecessor_journal_tip_sha256
    )
    if cursor.mode == "RETRY_PENDING":
        # This checks the entire first rejection chain: exact predecessor and
        # successor TDG7 plans, durable rejection identity, ledger prefix,
        # accepted state/time/serials, and the frozen half-cap relation.
        validate_retry_payload(
            cursor,
            prior,
            evolution_state_sha256=evolution_state_sha256,
        )
    evidence, updated = _validated_retry_outcome(
        cursor,
        prior,
        outcome,
        evolution_state_sha256=evolution_state_sha256,
    )
    if not isinstance(outcome, p15.TDG6TemporalRetryExhausted):
        raise HLT16LifecycleError("retry exhaustion requires TDG6TemporalRetryExhausted")
    if cursor.mode == "RETRY_PENDING":
        pending = cursor.payload["retry_successor_payload_or_none"]
        assert isinstance(pending, Mapping)  # checked by _validate_retry_payload
        successor = _wrap(p15._plan_mapping, pending["successor_plan"])
        if (
            not p15._same_binary64(
                evidence["initial_time"], float.fromhex(str(successor["current_hex"]))
            )
            or not p15._same_binary64(
                evidence["attempted_macro_step_size"],
                float.fromhex(str(successor["macro_width_hex"])),
            )
        ):
            raise HLT16LifecycleError(
                "pending retry exhaustion did not execute the exact durable TDG7 successor"
            )
        # _validated_retry_outcome already checks the exact increment against
        # ``prior``.  State the no-commit properties explicitly here so a
        # future ledger representation cannot accidentally promote a rejected
        # pending proposal into an accepted state transition.
        if (
            updated.last_accepted_time != prior.last_accepted_time
            or updated.accepted_macro_step_count != prior.accepted_macro_step_count
            or updated.accumulated_debit_vector != prior.accumulated_debit_vector
        ):
            raise HLT16LifecycleError(
                "pending retry exhaustion changed the accepted TDG6 boundary"
            )
    if evidence["initial_state_sha256"] != evolution_state_sha256:
        raise HLT16LifecycleError("retry exhaustion mutated the accepted state")
    return RetryExhaustion(
        reason=outcome.reason,
        member_key=str(cursor.payload["member_key"]),
        accepted_state_sha256=str(cursor.payload["accepted_state_sha256"]),
        accepted_time=_time(cursor.payload["accepted_boundary_time"]),
        updated_ledger=updated,
        durable_rejection_record_sha256=_sha("durable rejection record", durable_rejection_record_sha256, non_root=True),
        predecessor_cursor_sha256=cursor.sha256,
    )


def close_fine(
    cursor: p15.Proto15Cursor,
    prior_ledger: p15.TDG6TemporalLedger,
    ledger: p15.TDG6TemporalLedger,
    descriptor_hash: str,
    executed_plan: Mapping[str, object],
    accepted_time: Mapping[str, object],
    step_index: int,
    transaction_serial: int,
    predecessor_journal_tip_sha256: str,
) -> p15.Proto15Cursor:
    """Commit the sole accepted TDG6 fine path into a fresh PROTO18 cursor.

    The cursor intentionally retains its predecessor journal tip: the accepted
    journal record contains this cursor, so embedding that record hash here
    would create an unresolvable self-reference.  The checkpoint binds it.
    """
    prior = _checked_ledger(prior_ledger)
    updated = _checked_ledger(ledger)
    _require_proto18(cursor, prior)
    checkpoint_tip = _sha(
        "predecessor checkpoint journal tip",
        predecessor_journal_tip_sha256,
        non_root=True,
    )
    accepted = _time(accepted_time)
    plan = _wrap(p15._validated_accepted_plan, cursor, executed_plan, accepted)
    descriptor = _sha("accepted evolved descriptor", descriptor_hash, non_root=True)
    if descriptor == cursor.payload["accepted_state_sha256"]:
        raise HLT16LifecycleError("accepted fine must persist a new state descriptor")
    old_time = float.fromhex(str(cursor.payload["accepted_boundary_time"]["binary64_hex"]))
    new_time = float.fromhex(accepted["binary64_hex"])
    target = float.fromhex(str(cursor.payload["event_target_time"]["binary64_hex"]))
    if not old_time < new_time <= target:
        raise HLT16LifecycleError("accepted fine time does not advance within the event target")
    if updated.last_accepted_time.hex() != accepted["binary64_hex"]:
        raise HLT16LifecycleError("accepted fine ledger time differs from accepted state time")
    if (
        updated.current_macro_step_temporal_retry_count != 0
        or updated.accepted_macro_step_count != prior.accepted_macro_step_count + 1
        or updated.cumulative_temporal_retry_count != prior.cumulative_temporal_retry_count
        or updated.last_accepted_macro_step_temporal_retry_count != prior.current_macro_step_temporal_retry_count
        or updated.serialized_temporal_rejections != prior.serialized_temporal_rejections
        or any(new < old for new, old in zip(updated.accumulated_debit_vector, prior.accumulated_debit_vector))
    ):
        raise HLT16LifecycleError("accepted fine TDG6 ledger is not the exact monotone commit successor")
    if isinstance(step_index, bool) or not isinstance(step_index, int) or step_index <= cursor.payload["previous_step_index"]:
        raise HLT16LifecycleError("accepted fine step index does not advance")
    if isinstance(transaction_serial, bool) or not isinstance(transaction_serial, int) or transaction_serial <= cursor.payload["previous_transaction_serial"]:
        raise HLT16LifecycleError("accepted fine transaction serial does not advance")
    if plan["boundaries_hex"][-1] != accepted["binary64_hex"]:
        raise HLT16LifecycleError("accepted fine plan endpoint differs after validation")
    generation, serial, attempt_id = _next_identity(cursor)
    payload = dict(cursor.payload)
    payload.update(
        {
            "accepted_boundary_time": accepted,
            "accepted_state_sha256": descriptor,
            "previous_step_index": step_index,
            "previous_transaction_serial": transaction_serial,
            "TDG6_ledger_sha256": _ledger_hash(updated),
            "cursor_generation": generation,
            "attempt_serial": serial,
            "attempt_id": attempt_id,
            "mode": "FRESH_READY",
            "retry_successor_payload_or_none": None,
            # The accepted record embeds this successor cursor.  Therefore the
            # cursor can bind the already-known predecessor checkpoint tip but
            # cannot bind the hash of its own accepted record without a cycle.
            "journal_tip_sha256": checkpoint_tip,
            "cursor_chain_parent_sha256": cursor.sha256,
        }
    )
    payload.pop("cursor_chain_sha256", None)
    return _new_cursor(payload)


__all__ = [
    "FIRST_EVENT_TARGET_HEX", "GEN0_EVENT_INDEX", "GEN0_TIME_HEX", "HLT16LifecycleError",
    "PROTO17_ARTIFACT_ID", "PROTO18_ARTIFACT_ID", "RetryExhaustion", "bridge_gen0",
    "canonical_hash", "canonical_mapping", "close_fine", "close_retry", "retry_exhaustion",
    "typed_retry_from_evidence", "validate_retry_evidence", "validate_retry_payload",
]
