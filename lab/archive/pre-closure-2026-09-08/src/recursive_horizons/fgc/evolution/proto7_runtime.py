"""Fail-closed PROTO7 ownership for one complete unaccepted proposal.

PROTO7 changes no equation, source solver, threshold, integrator, or
non-source stop inherited from PROTO6.  It removes one numerical category
error: retryability is owned by acceptance state, not by an integrator's
stage label.  A source miss on the bitwise last accepted state is terminal.
A source-only miss anywhere in a proposal that has not committed -- including
the independently evaluated candidate endpoint -- may be rolled back and
retried.  Any non-source failure in that proposal vetoes retry and the
original chronological/within-stage stop priority remains authoritative.

The complete proposal is previewed before a retry decision.  Public evidence
therefore retains every failed evaluation, its raw source residual, its full
failure set, the accepted-state content hash, and exact no-commit facts.  The
caller owns the frozen retry budget, event-ledger durability, tracer/external
state, and cross-member rollback.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite
from numbers import Real
from typing import Any, Callable, Mapping

import numpy as np

from .boundary_domain import assess_causal_step
from .numerical_engine import (
    AcceptedStep,
    EvolutionRHS,
    EvolutionState,
    StageRecord,
    StepProposal,
    accept_step,
    array_content_sha256,
    propose_step,
)
from .proto5_runtime import (
    CFLRetryRequired,
    GR0RuntimeStageTransaction,
    GR0RuntimeStop,
    failed_gr0_universal_premises,
    gr0_universal_stage_snapshot,
)
from .proto6_runtime import (
    Proto6GR0EvolutionOperator,
    _latch_accepted_source_failure,
)


PROTO7_PROPOSAL_STAGE_NAMES = frozenset(
    {
        "rk4_k1",
        "rk4_k2",
        "rk4_k3",
        "rk4_k4",
        "ssprk3_s0",
        "ssprk3_s1",
        "ssprk3_s2",
        "candidate_endpoint",
    }
)


# The numerical source is exactly the PROTO6/CAL2 affine-root evaluator.
Proto7GR0EvolutionOperator = Proto6GR0EvolutionOperator


def _finite(name: str, value: object, *, positive: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a finite real scalar")
    answer = float(value)
    if not isfinite(answer):
        raise ValueError(f"{name} must be finite")
    if positive and answer <= 0.0:
        raise ValueError(f"{name} must be positive")
    return answer


def _nonnegative_integer(name: str, value: object) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return value


def _source_residual(record: StageRecord) -> float:
    try:
        value = record.rhs.diagnostics["source_residual_infinity"]
    except KeyError as error:
        raise ValueError("PROTO7 source evidence omits source_residual_infinity") from error
    return _finite("source_residual_infinity", value)


@dataclass(frozen=True, slots=True)
class FailedProposalEvaluation:
    """All failure information for one evaluation in an unaccepted proposal."""

    stage_name: str
    stage_time: float
    transaction_serial: int
    source_residual_infinity: float
    source_residual_maximum: float
    failures: tuple[str, ...]

    def __post_init__(self) -> None:
        if self.stage_name not in PROTO7_PROPOSAL_STAGE_NAMES:
            raise ValueError("failure evidence identifies no PROTO7 proposal stage")
        object.__setattr__(self, "stage_time", _finite("stage_time", self.stage_time))
        object.__setattr__(
            self,
            "transaction_serial",
            _nonnegative_integer("transaction_serial", self.transaction_serial),
        )
        object.__setattr__(
            self,
            "source_residual_infinity",
            _finite("source_residual_infinity", self.source_residual_infinity),
        )
        object.__setattr__(
            self,
            "source_residual_maximum",
            _finite(
                "source_residual_maximum",
                self.source_residual_maximum,
                positive=True,
            ),
        )
        if not self.failures or any(not isinstance(item, str) for item in self.failures):
            raise ValueError("failed evaluation must retain its complete failure tuple")


@dataclass(frozen=True, slots=True)
class ProposalFailureEvidence:
    """Complete decision evidence for one rejected or terminal proposal."""

    method: str
    initial_time: float
    proposed_step_size: float
    accepted_state_sha256: str
    accepted_step_index: int
    accepted_transaction_serial: int
    failed_evaluations: tuple[FailedProposalEvaluation, ...]
    complete_failure_set: tuple[str, ...]
    retryable_source_only: bool
    non_source_failure_vetoed_retry: bool
    selected_terminal_reason: str | None
    fields_bitwise_preserved: bool = True
    time_advanced: bool = False
    accepted_stage_count_advanced: bool = False
    causal_debit_advanced: bool = False
    preaccept_called: bool = False

    def __post_init__(self) -> None:
        object.__setattr__(self, "initial_time", _finite("initial_time", self.initial_time))
        object.__setattr__(
            self,
            "proposed_step_size",
            _finite("proposed_step_size", self.proposed_step_size, positive=True),
        )
        object.__setattr__(
            self,
            "accepted_step_index",
            _nonnegative_integer("accepted_step_index", self.accepted_step_index),
        )
        object.__setattr__(
            self,
            "accepted_transaction_serial",
            _nonnegative_integer(
                "accepted_transaction_serial", self.accepted_transaction_serial
            ),
        )
        if not isinstance(self.method, str) or not self.method:
            raise ValueError("method must be nonempty")
        if (
            not isinstance(self.accepted_state_sha256, str)
            or len(self.accepted_state_sha256) != 64
        ):
            raise ValueError("accepted_state_sha256 must be a SHA-256 digest")
        if not self.failed_evaluations or not self.complete_failure_set:
            raise ValueError("proposal failure evidence must be complete")
        observed = {
            failure
            for evaluation in self.failed_evaluations
            for failure in evaluation.failures
        }
        if observed != set(self.complete_failure_set):
            raise ValueError("complete_failure_set differs from failed evaluations")
        source_only = observed == {"newton_residual_limit"}
        if self.retryable_source_only is not source_only:
            raise ValueError("retry classification differs from complete proposal failures")
        if self.non_source_failure_vetoed_retry is not (not source_only):
            raise ValueError("non-source retry veto differs from complete failures")
        if source_only and self.selected_terminal_reason is not None:
            raise ValueError("retryable source-only proposal cannot select a terminal reason")
        if not source_only and self.selected_terminal_reason is None:
            raise ValueError("terminal proposal must retain its selected reason")
        if (
            self.fields_bitwise_preserved is not True
            or self.time_advanced is not False
            or self.accepted_stage_count_advanced is not False
            or self.causal_debit_advanced is not False
            or self.preaccept_called is not False
        ):
            raise ValueError("failed PROTO7 proposal violated no-commit semantics")


@dataclass(frozen=True, slots=True)
class AcceptedStateSourceFailureEvidence:
    """Terminal source miss observed before constructing any proposal."""

    accepted_time: float
    accepted_state_sha256: str
    accepted_step_index: int
    accepted_transaction_serial: int
    source_residual_infinity: float
    source_residual_maximum: float
    selected_terminal_reason: str = "newton_residual_limit"
    proposal_constructed: bool = False
    fields_bitwise_preserved: bool = True
    time_advanced: bool = False
    accepted_stage_count_advanced: bool = False
    causal_debit_advanced: bool = False
    preaccept_called: bool = False

    def __post_init__(self) -> None:
        object.__setattr__(self, "accepted_time", _finite("accepted_time", self.accepted_time))
        object.__setattr__(
            self,
            "accepted_step_index",
            _nonnegative_integer("accepted_step_index", self.accepted_step_index),
        )
        object.__setattr__(
            self,
            "accepted_transaction_serial",
            _nonnegative_integer(
                "accepted_transaction_serial", self.accepted_transaction_serial
            ),
        )
        object.__setattr__(
            self,
            "source_residual_infinity",
            _finite("source_residual_infinity", self.source_residual_infinity),
        )
        object.__setattr__(
            self,
            "source_residual_maximum",
            _finite(
                "source_residual_maximum",
                self.source_residual_maximum,
                positive=True,
            ),
        )
        if len(self.accepted_state_sha256) != 64:
            raise ValueError("accepted_state_sha256 must be a SHA-256 digest")
        if self.source_residual_infinity < self.source_residual_maximum:
            raise ValueError("accepted-state evidence does not cross the raw source gate")
        if (
            self.selected_terminal_reason != "newton_residual_limit"
            or self.proposal_constructed is not False
            or self.fields_bitwise_preserved is not True
            or self.time_advanced is not False
            or self.accepted_stage_count_advanced is not False
            or self.causal_debit_advanced is not False
            or self.preaccept_called is not False
        ):
            raise ValueError("accepted-state source stop semantics differ")


class Proto7TerminalStop(GR0RuntimeStop):
    """Terminal inherited scientific stop with complete PROTO7 evidence."""

    def __init__(
        self,
        stop: GR0RuntimeStop,
        evidence: ProposalFailureEvidence | AcceptedStateSourceFailureEvidence,
    ) -> None:
        super().__init__(stop.state, stop.failures)
        self.evidence = evidence


@dataclass(frozen=True, slots=True)
class Proto7StepAttempt:
    """Exactly one committed step or one retryable unaccepted proposal."""

    proposal: StepProposal
    accepted: AcceptedStep | None
    retry: ProposalFailureEvidence | None
    preaccept_payload: Any = None

    def __post_init__(self) -> None:
        if not isinstance(self.proposal, StepProposal):
            raise TypeError("proposal must be a StepProposal")
        if (self.accepted is None) == (self.retry is None):
            raise ValueError("PROTO7 attempt must be exactly accepted or retryable")
        if self.accepted is not None and not isinstance(self.accepted, AcceptedStep):
            raise TypeError("accepted must be an AcceptedStep")
        if self.retry is not None and self.retry.retryable_source_only is not True:
            raise ValueError("PROTO7 retry must be source-only")


def _complete_proposal_failure_evidence(
    transaction: GR0RuntimeStageTransaction,
    proposal: StepProposal,
    *,
    accepted_state_sha256: str,
    accepted_step_index: int,
    accepted_transaction_serial: int,
) -> ProposalFailureEvidence | None:
    if not proposal.stages or proposal.stages[-1].stage_name != "candidate_endpoint":
        raise ValueError("PROTO7 transaction requires an explicit candidate endpoint")
    speeds = tuple(transaction._stage_speed(record) for record in proposal.stages)
    boundary = assess_causal_step(
        transaction.causal_state,
        transaction.boundary_geometry,
        trial_time=proposal.final_time,
        stage_coordinate_speed_uppers=speeds[:-1],
        candidate_endpoint_coordinate_speed_upper=speeds[-1],
    )
    step_size = proposal.final_time - proposal.initial_time
    courant = step_size * boundary.interval_speed_upper / transaction.grid_spacing
    tolerance = 4096.0 * np.finfo(np.float64).eps * max(1.0, transaction.cfl_maximum)
    if courant > transaction.cfl_maximum + tolerance:
        raise CFLRetryRequired(
            observed_ratio=courant,
            maximum_ratio=transaction.cfl_maximum,
        )

    first_serial = transaction.state.last_transaction_serial + 1
    failed: list[FailedProposalEvaluation] = []
    for offset, record in enumerate(proposal.stages):
        serial = first_serial + offset
        snapshot = gr0_universal_stage_snapshot(
            record,
            transaction_serial=serial,
            boundary_margin_over_required_buffer=boundary.strict_margin_over_required_buffer,
            hat_normal_factor=transaction.hat_normal_factor,
        )
        failures = failed_gr0_universal_premises(snapshot, transaction.thresholds)
        if failures:
            failed.append(
                FailedProposalEvaluation(
                    stage_name=record.stage_name,
                    stage_time=record.time,
                    transaction_serial=serial,
                    source_residual_infinity=_source_residual(record),
                    source_residual_maximum=transaction.thresholds.source_residual_maximum,
                    failures=failures,
                )
            )
    if not failed:
        if not boundary.admissible:
            raise RuntimeError("PROTO7 boundary preview escaped the universal stop")
        return None

    complete_failure_set = tuple(
        name
        for name in transaction_stop_priority(transaction)
        if any(name in item.failures for item in failed)
    )
    source_only = complete_failure_set == ("newton_residual_limit",)
    selected_reason: str | None = None
    if not source_only:
        # Re-run the unchanged real transaction.  It owns chronological and
        # within-stage priority and latches the selected scientific stop.
        try:
            transaction(proposal)
        except GR0RuntimeStop as stop:
            selected_reason = stop.reason
        else:
            raise RuntimeError("real inherited transaction accepted a preview stop")

    return ProposalFailureEvidence(
        method=proposal.method,
        initial_time=proposal.initial_time,
        proposed_step_size=step_size,
        accepted_state_sha256=accepted_state_sha256,
        accepted_step_index=accepted_step_index,
        accepted_transaction_serial=accepted_transaction_serial,
        failed_evaluations=tuple(failed),
        complete_failure_set=complete_failure_set,
        retryable_source_only=source_only,
        non_source_failure_vetoed_retry=not source_only,
        selected_terminal_reason=selected_reason,
    )


def transaction_stop_priority(transaction: GR0RuntimeStageTransaction) -> tuple[str, ...]:
    """Recover the inherited ordering without defining a second stop list."""

    # A healthy synthetic snapshot is unnecessary: failed premise tuples are
    # already emitted in inherited order.  Collect their union by the stable
    # canonical order exported by proto5_runtime.
    from .proto5_runtime import GR0_RUNTIME_STOP_PRIORITY

    if not isinstance(transaction, GR0RuntimeStageTransaction):
        raise TypeError("transaction must be a GR0RuntimeStageTransaction")
    return tuple(GR0_RUNTIME_STOP_PRIORITY)


def attempt_proto7_step(
    *,
    method: str,
    time: Real,
    step_size: Real,
    state: EvolutionState,
    rhs: Callable[[float, EvolutionState], EvolutionRHS],
    projector: Callable[[float, EvolutionState], EvolutionState] | None,
    transaction: GR0RuntimeStageTransaction,
    previous_step_index: int,
    previous_transaction_serial: int,
    preaccept: Callable[[StepProposal], Any] | None = None,
) -> Proto7StepAttempt:
    """Evaluate one complete PROTO7 transaction without weakening a stop."""

    if not isinstance(state, EvolutionState):
        raise TypeError("state must be EvolutionState")
    accepted_time = _finite("time", time)
    proposed_size = _finite("step_size", step_size, positive=True)
    accepted_index = _nonnegative_integer("previous_step_index", previous_step_index)
    accepted_serial = _nonnegative_integer(
        "previous_transaction_serial", previous_transaction_serial
    )
    if transaction.state.first_failed_premise is not None:
        raise GR0RuntimeStop(
            transaction.state,
            (transaction.state.first_failed_premise,),
        )
    accepted_hash = array_content_sha256(state.u, state.p, state.q)
    accepted_rhs = rhs(accepted_time, state)
    accepted_record = StageRecord(
        "accepted_state_source_precheck", accepted_time, state, accepted_rhs
    )
    try:
        _latch_accepted_source_failure(transaction, accepted_record)
    except GR0RuntimeStop as stop:
        evidence = AcceptedStateSourceFailureEvidence(
            accepted_time=accepted_time,
            accepted_state_sha256=accepted_hash,
            accepted_step_index=accepted_index,
            accepted_transaction_serial=accepted_serial,
            source_residual_infinity=_source_residual(accepted_record),
            source_residual_maximum=transaction.thresholds.source_residual_maximum,
        )
        raise Proto7TerminalStop(stop, evidence) from stop

    proposal = propose_step(
        method=method,
        time=accepted_time,
        step_size=proposed_size,
        state=state,
        rhs=rhs,
        projector=projector,
    )
    evidence = _complete_proposal_failure_evidence(
        transaction,
        proposal,
        accepted_state_sha256=accepted_hash,
        accepted_step_index=accepted_index,
        accepted_transaction_serial=accepted_serial,
    )
    if evidence is not None:
        if evidence.retryable_source_only:
            return Proto7StepAttempt(proposal=proposal, accepted=None, retry=evidence)
        stop = GR0RuntimeStop(transaction.state, evidence.failed_evaluations[0].failures)
        # Use the selected reason and complete failure tuple from the real
        # transaction, not the first preview record when later stages vetoed.
        if stop.reason != evidence.selected_terminal_reason:
            raise RuntimeError("PROTO7 terminal evidence differs from inherited latch")
        raise Proto7TerminalStop(stop, evidence)

    payload = None if preaccept is None else preaccept(proposal)
    accepted = accept_step(
        proposal,
        previous_step_index=accepted_index,
        previous_transaction_serial=accepted_serial,
        guards=(transaction,),
    )
    return Proto7StepAttempt(
        proposal=proposal,
        accepted=accepted,
        retry=None,
        preaccept_payload=payload,
    )


__all__ = [
    "AcceptedStateSourceFailureEvidence",
    "FailedProposalEvaluation",
    "PROTO7_PROPOSAL_STAGE_NAMES",
    "ProposalFailureEvidence",
    "Proto7GR0EvolutionOperator",
    "Proto7StepAttempt",
    "Proto7TerminalStop",
    "attempt_proto7_step",
]
