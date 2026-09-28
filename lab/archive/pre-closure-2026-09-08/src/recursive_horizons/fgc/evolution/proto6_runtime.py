"""Fail-closed PROTO6 ownership for accepted and trial source failures.

PROTO5 correctly evaluated every Runge--Kutta stage, but its source evaluator
raised before the runtime transaction could distinguish an accepted-state
failure from a source-only failure in a still-unaccepted trial.  PROTO6 keeps
the physical equations, affine root, raw tolerance, and every non-source stop
unchanged.  It changes only the transaction boundary:

* the raw source gate is evaluated on the bitwise last accepted state before
  every proposal and remains a terminal ``newton_residual_limit``;
* a shadow copy of the complete PROTO5 transaction evaluates the proposal;
* only a source-only miss at an internal RK stage is returned as retryable;
* the real transaction, causal debit, tracer preview, and accepted state are
  untouched by a rejected trial; and
* the explicit candidate endpoint is not an internal stage and therefore
  retains PROTO5's terminal ownership.

The caller owns the frozen retry count/minimum-step policy and must serialize
every :class:`TrialSourceRetryEvidence` before making another proposal.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite
from numbers import Real
from typing import Any, Callable, Mapping

import numpy as np

from .boundary_domain import assess_causal_step
from .cal2_source_diagnosis import DiagnosticGR0EvolutionOperator
from .numerical_engine import (
    AcceptedStep,
    EvolutionRHS,
    EvolutionState,
    StageRecord,
    StepProposal,
    accept_step,
    propose_step,
)
from .proto5_runtime import (
    CFLRetryRequired,
    GR0RuntimeMonitorState,
    GR0RuntimeStageTransaction,
    GR0RuntimeStop,
    failed_gr0_universal_premises,
    gr0_universal_stage_snapshot,
)


PROTO6_INTERNAL_STAGE_NAMES = frozenset(
    {"rk4_k1", "rk4_k2", "rk4_k3", "rk4_k4", "ssprk3_s0", "ssprk3_s1", "ssprk3_s2"}
)


class Proto6GR0EvolutionOperator:
    """PROTO6-owned exposure of the unchanged affine GR-0 source result.

    CAL2's diagnostic evaluator deliberately was not an authorized evolution
    backend.  HLT4 promotes only its same-equation, same-refinement output into
    this transaction-aware adapter; :func:`attempt_proto6_step` still prevents
    a failed raw gate from entering an accepted state.
    """

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        self._diagnostic = DiagnosticGR0EvolutionOperator(*args, **kwargs)

    def __call__(self, time: float, state: EvolutionState) -> EvolutionRHS:
        return self._diagnostic(time, state)


def _finite(name: str, value: object, *, positive: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a finite real scalar")
    answer = float(value)
    if not isfinite(answer):
        raise ValueError(f"{name} must be finite")
    if positive and answer <= 0.0:
        raise ValueError(f"{name} must be positive")
    return answer


def _source_residual(record: StageRecord) -> float:
    try:
        value = record.rhs.diagnostics["source_residual_infinity"]
    except KeyError as error:
        raise ValueError("PROTO6 source evidence omits source_residual_infinity") from error
    return _finite("source_residual_infinity", value)


@dataclass(frozen=True, slots=True)
class TrialSourceStageEvidence:
    """One failed source-only stage inside a rejected proposal."""

    stage_name: str
    stage_time: float
    source_residual_infinity: float
    transaction_serial: int

    def __post_init__(self) -> None:
        if self.stage_name not in PROTO6_INTERNAL_STAGE_NAMES:
            raise ValueError("retry stage evidence must identify an internal RK stage")
        object.__setattr__(self, "stage_time", _finite("stage_time", self.stage_time))
        object.__setattr__(
            self,
            "source_residual_infinity",
            _finite("source_residual_infinity", self.source_residual_infinity),
        )
        if (
            isinstance(self.transaction_serial, bool)
            or not isinstance(self.transaction_serial, int)
            or self.transaction_serial < 0
        ):
            raise ValueError("transaction_serial must be nonnegative")


@dataclass(frozen=True, slots=True)
class TrialSourceRetryEvidence:
    """Complete public evidence for one rejected, unaccepted proposal."""

    initial_time: float
    proposed_step_size: float
    failed_stage_name: str
    failed_stage_time: float
    source_residual_infinity: float
    source_residual_maximum: float
    failed_transaction_serial: int
    failed_stages: tuple[TrialSourceStageEvidence, ...]
    accepted_state_bitwise_preserved: bool = True
    time_advanced: bool = False
    transaction_advanced: bool = False
    causal_debit_advanced: bool = False

    def __post_init__(self) -> None:
        for name in (
            "initial_time",
            "proposed_step_size",
            "failed_stage_time",
            "source_residual_infinity",
            "source_residual_maximum",
        ):
            object.__setattr__(
                self,
                name,
                _finite(name, getattr(self, name), positive=name in {"proposed_step_size", "source_residual_maximum"}),
            )
        if self.failed_stage_name not in PROTO6_INTERNAL_STAGE_NAMES:
            raise ValueError("retry evidence must identify an internal RK stage")
        if (
            isinstance(self.failed_transaction_serial, bool)
            or not isinstance(self.failed_transaction_serial, int)
            or self.failed_transaction_serial < 0
        ):
            raise ValueError("failed_transaction_serial must be nonnegative")
        if self.source_residual_infinity < self.source_residual_maximum:
            raise ValueError("retry evidence does not cross the raw source threshold")
        if not self.failed_stages:
            raise ValueError("retry evidence must expose every failed internal stage")
        first = self.failed_stages[0]
        if (
            first.stage_name != self.failed_stage_name
            or first.stage_time != self.failed_stage_time
            or first.source_residual_infinity != self.source_residual_infinity
            or first.transaction_serial != self.failed_transaction_serial
            or any(
                item.source_residual_infinity < self.source_residual_maximum
                for item in self.failed_stages
            )
        ):
            raise ValueError("primary retry fields differ from complete stage evidence")
        if (
            self.accepted_state_bitwise_preserved is not True
            or self.time_advanced is not False
            or self.transaction_advanced is not False
            or self.causal_debit_advanced is not False
        ):
            raise ValueError("rejected PROTO6 trials must have exact rollback semantics")


@dataclass(frozen=True, slots=True)
class Proto6StepAttempt:
    """Exactly one accepted step or one retryable rejected trial."""

    proposal: StepProposal
    accepted: AcceptedStep | None
    retry: TrialSourceRetryEvidence | None
    preaccept_payload: Any = None

    def __post_init__(self) -> None:
        if not isinstance(self.proposal, StepProposal):
            raise TypeError("proposal must be a StepProposal")
        if (self.accepted is None) == (self.retry is None):
            raise ValueError("PROTO6 attempt must be exactly accepted or retryable")
        if self.accepted is not None and not isinstance(self.accepted, AcceptedStep):
            raise TypeError("accepted must be an AcceptedStep")


def _latch_accepted_source_failure(
    transaction: GR0RuntimeStageTransaction,
    record: StageRecord,
) -> None:
    """Keep an accepted-state raw source miss terminal and immutable."""

    if transaction.state.first_failed_premise is not None:
        raise GR0RuntimeStop(
            transaction.state,
            (transaction.state.first_failed_premise,),
        )
    residual = _source_residual(record)
    if residual < transaction.thresholds.source_residual_maximum:
        return
    prior = transaction.state
    serial = max(prior.last_transaction_serial, 0)
    transaction.state = GR0RuntimeMonitorState(
        accepted_stage_count=prior.accepted_stage_count,
        last_transaction_serial=prior.last_transaction_serial,
        last_accepted_time=prior.last_accepted_time,
        first_failed_premise="newton_residual_limit",
        first_failed_transaction_serial=serial,
        first_failed_time=record.time,
    )
    raise GR0RuntimeStop(transaction.state, ("newton_residual_limit",))


def _trial_retry_evidence(
    transaction: GR0RuntimeStageTransaction,
    proposal: StepProposal,
) -> TrialSourceRetryEvidence | None:
    """Preview the complete transaction and isolate only a source-only miss."""

    if not proposal.stages or proposal.stages[-1].stage_name != "candidate_endpoint":
        raise ValueError("PROTO6 transaction requires an explicit candidate endpoint")
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
    tolerance = 4096.0 * np.finfo(np.float64).eps * max(
        1.0, transaction.cfl_maximum
    )
    if courant > transaction.cfl_maximum + tolerance:
        raise CFLRetryRequired(
            observed_ratio=courant,
            maximum_ratio=transaction.cfl_maximum,
        )

    first_serial = transaction.state.last_transaction_serial + 1
    failures_by_stage: list[tuple[StageRecord, int, tuple[str, ...]]] = []
    for offset, record in enumerate(proposal.stages):
        serial = first_serial + offset
        snapshot = gr0_universal_stage_snapshot(
            record,
            transaction_serial=serial,
            boundary_margin_over_required_buffer=(
                boundary.strict_margin_over_required_buffer
            ),
            hat_normal_factor=transaction.hat_normal_factor,
        )
        failures = failed_gr0_universal_premises(
            snapshot,
            transaction.thresholds,
        )
        if failures:
            failures_by_stage.append((record, serial, failures))
    if not failures_by_stage:
        if not boundary.admissible:
            raise RuntimeError("PROTO6 boundary preview escaped the universal stop")
        return None

    retryable = all(
        failures == ("newton_residual_limit",)
        and record.stage_name in PROTO6_INTERNAL_STAGE_NAMES
        for record, _serial, failures in failures_by_stage
    )
    if retryable:
        stage_evidence = tuple(
            TrialSourceStageEvidence(
                stage_name=record.stage_name,
                stage_time=record.time,
                source_residual_infinity=_source_residual(record),
                transaction_serial=serial,
            )
            for record, serial, _failures in failures_by_stage
        )
        first = stage_evidence[0]
        return TrialSourceRetryEvidence(
            initial_time=proposal.initial_time,
            proposed_step_size=step_size,
            failed_stage_name=first.stage_name,
            failed_stage_time=first.stage_time,
            source_residual_infinity=first.source_residual_infinity,
            source_residual_maximum=transaction.thresholds.source_residual_maximum,
            failed_transaction_serial=first.transaction_serial,
            failed_stages=stage_evidence,
        )

    # The complete preview found an endpoint or non-source failure. Re-run the
    # unchanged real transaction so its original chronological/within-stage
    # priority owns and latches the first scientific failure.
    transaction(proposal)
    raise RuntimeError("real PROTO5 transaction unexpectedly accepted a preview stop")


def attempt_proto6_step(
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
) -> Proto6StepAttempt:
    """Evaluate one PROTO6 proposal without hiding or weakening any stop.

    ``rhs`` must expose the unchanged affine-root diagnostics even when a raw
    residual misses its tolerance.  The accepted-state precheck is performed
    before the proposal.  A caller-supplied ``preaccept`` (the production use
    is the normal-flow tracer preview) runs only after the shadow transaction
    has proved that the proposal is neither retryable nor terminal.
    """

    if not isinstance(state, EvolutionState):
        raise TypeError("state must be EvolutionState")
    accepted_time = _finite("time", time)
    proposed_size = _finite("step_size", step_size, positive=True)
    accepted_rhs = rhs(accepted_time, state)
    accepted_record = StageRecord(
        "accepted_state_source_precheck",
        accepted_time,
        state,
        accepted_rhs,
    )
    _latch_accepted_source_failure(transaction, accepted_record)

    proposal = propose_step(
        method=method,
        time=accepted_time,
        step_size=proposed_size,
        state=state,
        rhs=rhs,
        projector=projector,
    )
    retry = _trial_retry_evidence(transaction, proposal)
    if retry is not None:
        return Proto6StepAttempt(proposal=proposal, accepted=None, retry=retry)

    payload = None if preaccept is None else preaccept(proposal)
    accepted = accept_step(
        proposal,
        previous_step_index=previous_step_index,
        previous_transaction_serial=previous_transaction_serial,
        guards=(transaction,),
    )
    return Proto6StepAttempt(
        proposal=proposal,
        accepted=accepted,
        retry=None,
        preaccept_payload=payload,
    )


def restore_common_event_snapshots(
    members: Mapping[str, Any],
    snapshots: Mapping[str, Mapping[str, Any]],
) -> None:
    """Restore every member or reject an incomplete cross-member rollback."""

    if set(members) != set(snapshots):
        raise ValueError("PROTO6 common-event rollback member keys differ")
    for key, member in members.items():
        restore = getattr(member, "restore", None)
        if not callable(restore):
            raise TypeError("PROTO6 rollback members must expose restore(snapshot)")
        restore(snapshots[key])


__all__ = [
    "PROTO6_INTERNAL_STAGE_NAMES",
    "Proto6GR0EvolutionOperator",
    "Proto6StepAttempt",
    "TrialSourceStageEvidence",
    "TrialSourceRetryEvidence",
    "attempt_proto6_step",
    "restore_common_event_snapshots",
]
