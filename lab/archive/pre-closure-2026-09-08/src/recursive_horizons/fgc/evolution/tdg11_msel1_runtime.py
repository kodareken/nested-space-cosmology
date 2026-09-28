"""Guarded SHADOW family builder for prospective TDG11 selection.

This module constructs one outer, two medium, and four fine PROTO7-guarded
shadow proposals from a bitwise accepted state.  It records the family and
returns it.  It does not classify channels, select a candidate, mutate a
production temporal ledger, admit or commit an accepted successor, touch a
store, or publish.

The only adaptation from the sealed PROTO7 attempt is the local proposer:
``numeric_engine.propose_step`` for legacy recorded arithmetic and
``propose_compensated_step`` for compensated arithmetic.  The proposer is
passed locally.  No historical module global is assigned.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from math import isfinite
from numbers import Real
from typing import Callable

import numpy as np

from . import numerical_engine as numeric_engine
from . import tdg5_stage_complete_refinement_runtime as tdg5
from . import tdg6_temporal_admission_runtime as tdg6
from .numerical_engine import (
    COMPARATOR_METHOD,
    PRIMARY_METHOD,
    EvolutionRHS,
    EvolutionState,
    StageRecord,
    StepProposal,
    accept_step,
    array_content_sha256,
)
from .proto5_runtime import GR0RuntimeStageTransaction, GR0RuntimeStop
from .proto6_runtime import _latch_accepted_source_failure
from .proto7_runtime import (
    AcceptedStateSourceFailureEvidence,
    ProposalFailureEvidence,
    Proto7StepAttempt,
    Proto7TerminalStop,
    _complete_proposal_failure_evidence,
    _source_residual,
)
from .tdg11_compensated_rk import (
    ARITHMETIC_ID as COMPENSATED_ARITHMETIC_ID,
    propose_compensated_step,
)
from .tdg11_msel1_reconstruction import (
    ORIGINAL_ARITHMETIC_ID,
    ValidatedRecordedFamily,
    validate_recorded_family,
)


# Historical private helpers.  No public clone API exists for GR-0
# transactions or normal-flow tracers.  These aliases isolate the
# underscore helpers used by the sealed TDG5/TDG6 shadow builders.
# This module never assigns into historical module globals.
_clone_gr0_transaction = tdg5._clone_gr0_transaction
_clone_tracers = tdg6._clone_tracers
_tracer_snapshot = tdg6._tracer_snapshot
_state_hash = tdg5._state_hash

_LEVEL_COUNTS: tuple[tuple[str, int], ...] = (("outer", 1), ("medium", 2), ("fine", 4))
_STAGES_PER_PROPOSAL = {PRIMARY_METHOD: 5, COMPARATOR_METHOD: 4}
_RECOGNIZED_ARITHMETIC = frozenset({ORIGINAL_ARITHMETIC_ID, COMPENSATED_ARITHMETIC_ID})
_SUCCESS_PRECHECKS = 7


class _CountingRHS:
    """Instrument actual RHS invocations; counts are never fabricated."""

    __slots__ = ("inner", "calls")

    def __init__(self, inner: Callable[[float, EvolutionState], EvolutionRHS]) -> None:
        self.inner = inner
        self.calls = 0

    def __call__(self, time: float, state: EvolutionState) -> EvolutionRHS:
        self.calls += 1
        return self.inner(time, state)


class _Progress:
    __slots__ = ("prechecks", "stages")

    def __init__(self) -> None:
        self.prechecks = 0
        self.stages = 0


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


def _same_binary64(left: float, right: float) -> bool:
    return np.float64(left).tobytes() == np.float64(right).tobytes()


def _immutable_array(name: str, value: object, *, ndim: int = 1) -> np.ndarray:
    array = np.asarray(value, dtype=np.float64)
    if array.ndim != ndim or np.any(~np.isfinite(array)):
        raise ValueError(f"{name} must be a finite {ndim}-dimensional array")
    answer = np.ascontiguousarray(array).copy()
    answer.setflags(write=False)
    return answer


def _array_identity(array: object) -> tuple[object, ...]:
    data = np.asarray(array)
    return (
        id(array) if isinstance(array, np.ndarray) else id(data),
        data.dtype.str,
        tuple(int(size) for size in data.shape),
        tuple(int(stride) for stride in data.strides),
        data.tobytes(),
        bool(data.flags.writeable),
        bool(data.flags.c_contiguous),
        bool(data.flags.f_contiguous),
        bool(data.flags.owndata),
        bool(data.flags.aligned),
    )


def _tracer_fingerprint(tracers: object) -> tuple[object, ...]:
    event_times = tuple(
        _array_identity(item) for item in getattr(tracers, "event_proper_times")
    )

    event_fields = tuple(
        _array_identity(item) for item in getattr(tracers, "event_fields")
    )
    return (
        id(tracers),
        _array_identity(getattr(tracers, "labels")),
        _array_identity(getattr(tracers, "positions")),
        _array_identity(getattr(tracers, "proper_times")),
        event_times,
        event_fields,
        float(getattr(tracers, "cutoff")).hex(),
        float(getattr(tracers, "outer_radius")).hex(),
    )


def _metadata_snapshot(value: object) -> object:
    """Copy scalar metadata values, including signed-zero float identity."""
    if type(value) in (float, np.float64):
        return ("binary64", float(value).hex())
    if type(value) in (str, int, bool) or value is None:
        return (type(value).__name__, value)
    if type(value) is dict:
        return tuple(
            (key, _metadata_snapshot(item)) for key, item in sorted(value.items())
        )
    if type(value) in (tuple, list):
        return (type(value).__name__, tuple(_metadata_snapshot(item) for item in value))
    raise TypeError("unsupported transaction metadata")


def _input_fingerprint(
    *,
    state: EvolutionState,
    transaction: GR0RuntimeStageTransaction,
    tracers: object,
    coordinates: object,
) -> tuple[object, ...]:
    return (
        id(state),
        _array_identity(state.u),
        _array_identity(state.p),
        _array_identity(state.q),
        id(transaction),
        id(transaction.state),
        id(transaction.causal_state),
        id(transaction.thresholds),
        id(transaction.boundary_geometry),
        _metadata_snapshot(asdict(transaction.state)),
        _metadata_snapshot(asdict(transaction.causal_state)),
        _metadata_snapshot(asdict(transaction.thresholds)),
        _metadata_snapshot(asdict(transaction.boundary_geometry)),
        float(transaction.grid_spacing).hex(),
        float(transaction.cfl_maximum).hex(),
        float(transaction.hat_normal_factor).hex(),
        _tracer_fingerprint(tracers),
        _array_identity(coordinates),
    )


def _select_proposer(arithmetic_id: str):
    if arithmetic_id == ORIGINAL_ARITHMETIC_ID:
        return numeric_engine.propose_step
    if arithmetic_id == COMPENSATED_ARITHMETIC_ID:
        return propose_compensated_step
    raise ValueError("unrecognized arithmetic")


def _subdivision_boundaries(
    start: float, width: float, count: int
) -> tuple[float, ...]:
    step = width / count
    if not isfinite(step) or step <= 0.0:
        raise ValueError("shadow subdivision step must be finite and positive")
    boundaries = tuple(start + index * step for index in range(count + 1))
    expected_final = start + width
    if not isfinite(expected_final) or not _same_binary64(
        boundaries[-1], expected_final
    ):
        raise ValueError("h, h/2, h/4 time lattice is not bitwise uniform")
    for left, right in zip(boundaries, boundaries[1:]):
        if not _same_binary64(right - left, step):
            raise ValueError("h, h/2, h/4 time lattice is not bitwise uniform")
        midpoint = left + (right - left) / 2.0
        if not isfinite(midpoint) or not (left < midpoint < right):
            raise ValueError("stage midpoint is not representable on the time lattice")
    return boundaries


def _family_boundaries(start: float, width: float) -> tuple[tuple[float, ...], ...]:
    half = width / 2.0
    quarter = width / 4.0
    if any(not isfinite(value) or value <= 0.0 for value in (width, half, quarter)):
        raise ValueError("h, h/2, and h/4 must be finite and positive")
    outer = _subdivision_boundaries(start, width, 1)
    medium = _subdivision_boundaries(start, width, 2)
    fine = _subdivision_boundaries(start, width, 4)
    if not (
        _same_binary64(outer[0], medium[0])
        and _same_binary64(outer[0], fine[0])
        and _same_binary64(outer[-1], medium[-1])
        and _same_binary64(outer[-1], fine[-1])
    ):
        raise ValueError("h, h/2, h/4 time lattice is not bitwise uniform")
    return outer, medium, fine


@dataclass(frozen=True, slots=True)
class GuardedShadowFamily:
    """One recorded 1/2/4 shadow family with instrumented RHS accounting."""

    paths: tuple[tdg6.TDG6ShadowPath, ...]
    recorded: ValidatedRecordedFamily
    accepted_source_precheck_count: int
    stage_record_count: int
    rhs_call_count: int

    def __post_init__(self) -> None:
        for name in (
            "accepted_source_precheck_count",
            "stage_record_count",
            "rhs_call_count",
        ):
            _nonnegative_integer(name, getattr(self, name))
        if len(self.paths) != 3:
            raise ValueError(
                "guarded shadow family must contain outer/medium/fine paths"
            )
        levels = tuple(path.level for path in self.paths)
        if levels != ("outer", "medium", "fine"):
            raise ValueError("guarded shadow paths are not in outer/medium/fine order")
        expected = tuple(count for _level, count in _LEVEL_COUNTS)
        observed = tuple(len(path.attempts) for path in self.paths)
        if observed != expected:
            raise ValueError("guarded shadow family does not have 1/2/4 proposals")
        if not isinstance(self.recorded, ValidatedRecordedFamily):
            raise TypeError("recorded must be ValidatedRecordedFamily")
        method = self.paths[0].method
        if not all(path.method == method for path in self.paths):
            raise ValueError("shadow paths use different methods")
        if self.recorded.method != method:
            raise ValueError("recorded family method differs from shadow paths")
        if self.accepted_source_precheck_count != _SUCCESS_PRECHECKS:
            raise ValueError(
                "successful family must retain seven accepted-state prechecks"
            )
        stages = sum(path.stage_record_count for path in self.paths)
        if self.stage_record_count != stages:
            raise ValueError("stage_record_count differs from retained shadow stages")
        expected_stages = _SUCCESS_PRECHECKS * _STAGES_PER_PROPOSAL[method]
        if stages != expected_stages:
            raise ValueError("retained stage/endpoint count differs from the method")
        if self.rhs_call_count != self.accepted_source_precheck_count + stages:
            raise ValueError(
                "RHS call count differs from prechecks plus retained stages"
            )
        recorded_paths = tuple(
            tuple(attempt.proposal for attempt in path.attempts) for path in self.paths
        )
        for left, right in zip(recorded_paths, self.recorded.paths, strict=True):
            if len(left) != len(right) or any(
                proposal is not recorded for proposal, recorded in zip(left, right)
            ):
                raise ValueError("recorded family does not retain the shadow proposals")


class TDG11ShadowPremiseStop(RuntimeError):
    """Diagnostic family stop. Not a scientific nonpass or informal retry."""

    def __init__(
        self,
        *,
        level: str,
        proposal_index: int,
        cause: BaseException | None = None,
        source_retry: ProposalFailureEvidence | None = None,
        accepted_source_precheck_count: int,
        stage_record_count: int,
        rhs_call_count: int,
    ) -> None:
        if level not in {name for name, _count in _LEVEL_COUNTS}:
            raise ValueError("unknown TDG11 shadow level")
        index = _nonnegative_integer("proposal_index", proposal_index)
        expected = dict(_LEVEL_COUNTS)[level]
        if index >= expected:
            raise ValueError("proposal_index exceeds the shadow level")
        if (cause is None) == (source_retry is None):
            raise ValueError(
                "TDG11 shadow stop requires exactly one cause or source retry"
            )
        if source_retry is not None:
            if not isinstance(source_retry, ProposalFailureEvidence):
                raise TypeError("source_retry must be ProposalFailureEvidence")
            if source_retry.retryable_source_only is not True:
                raise ValueError(
                    "source retry must remain PROTO7-owned and source-only"
                )
        for name, value in (
            ("accepted_source_precheck_count", accepted_source_precheck_count),
            ("stage_record_count", stage_record_count),
            ("rhs_call_count", rhs_call_count),
        ):
            _nonnegative_integer(name, value)
        self.level = level
        self.proposal_index = index
        self.cause = cause
        self.source_retry = source_retry
        self.accepted_source_precheck_count = accepted_source_precheck_count
        self.stage_record_count = stage_record_count
        self.rhs_call_count = rhs_call_count
        self.source_retry_owned_by_PROTO7 = source_retry is not None
        self.temporal_retry = False
        self.informal_retry = False
        self.real_transaction_advanced = False
        self.accepted_state_advanced = False
        self.real_tracer_advanced = False
        self.scientific_nonpass = False
        kind = "source_retry" if source_retry is not None else type(cause).__name__
        super().__init__(f"TDG11 {level}_{index} shadow path stopped: {kind}")


class TDG11ShadowIntegrityError(RuntimeError):
    """Accepted input objects mutated during diagnostic shadow construction."""


def _raise_premise_stop(
    *,
    level: str,
    proposal_index: int,
    progress: _Progress,
    rhs_calls: int,
    cause: BaseException | None = None,
    source_retry: ProposalFailureEvidence | None = None,
) -> None:
    stop = TDG11ShadowPremiseStop(
        level=level,
        proposal_index=proposal_index,
        cause=cause,
        source_retry=source_retry,
        accepted_source_precheck_count=progress.prechecks,
        stage_record_count=progress.stages,
        rhs_call_count=rhs_calls,
    )
    if cause is not None:
        raise stop from cause
    raise stop


def _guarded_shadow_attempt(
    *,
    proposer,
    method: str,
    time: float,
    step_size: float,
    state: EvolutionState,
    rhs: Callable[[float, EvolutionState], EvolutionRHS],
    projector: Callable[[float, EvolutionState], EvolutionState] | None,
    transaction: GR0RuntimeStageTransaction,
    previous_step_index: int,
    previous_transaction_serial: int,
    preaccept: Callable[[StepProposal], object],
    progress: _Progress,
) -> Proto7StepAttempt:
    """PROTO7 order with a locally supplied proposer and no module-global write."""

    if not isinstance(state, EvolutionState):
        raise TypeError("state must be EvolutionState")
    if transaction.state.first_failed_premise is not None:
        raise GR0RuntimeStop(
            transaction.state,
            (transaction.state.first_failed_premise,),
        )
    accepted_hash = array_content_sha256(state.u, state.p, state.q)
    accepted_rhs = rhs(time, state)
    progress.prechecks += 1
    accepted_record = StageRecord(
        "accepted_state_source_precheck", time, state, accepted_rhs
    )
    try:
        _latch_accepted_source_failure(transaction, accepted_record)
    except GR0RuntimeStop as stop:
        evidence = AcceptedStateSourceFailureEvidence(
            accepted_time=time,
            accepted_state_sha256=accepted_hash,
            accepted_step_index=previous_step_index,
            accepted_transaction_serial=previous_transaction_serial,
            source_residual_infinity=_source_residual(accepted_record),
            source_residual_maximum=transaction.thresholds.source_residual_maximum,
        )
        raise Proto7TerminalStop(stop, evidence) from stop

    proposal = proposer(
        method=method,
        time=time,
        step_size=step_size,
        state=state,
        rhs=rhs,
        projector=projector,
    )
    progress.stages += len(proposal.stages)
    evidence = _complete_proposal_failure_evidence(
        transaction,
        proposal,
        accepted_state_sha256=accepted_hash,
        accepted_step_index=previous_step_index,
        accepted_transaction_serial=previous_transaction_serial,
    )
    if evidence is not None:
        if evidence.retryable_source_only:
            return Proto7StepAttempt(proposal=proposal, accepted=None, retry=evidence)
        stop = GR0RuntimeStop(
            transaction.state, evidence.failed_evaluations[0].failures
        )
        if stop.reason != evidence.selected_terminal_reason:
            raise RuntimeError("PROTO7 terminal evidence differs from inherited latch")
        raise Proto7TerminalStop(stop, evidence)

    payload = preaccept(proposal)
    accepted = accept_step(
        proposal,
        previous_step_index=previous_step_index,
        previous_transaction_serial=previous_transaction_serial,
        guards=(transaction,),
    )
    return Proto7StepAttempt(
        proposal=proposal,
        accepted=accepted,
        retry=None,
        preaccept_payload=payload,
    )


def _shadow_path(
    *,
    level: str,
    method: str,
    boundaries: tuple[float, ...],
    state: EvolutionState,
    rhs: _CountingRHS,
    projector: Callable[[float, EvolutionState], EvolutionState] | None,
    transaction: GR0RuntimeStageTransaction,
    tracers: object,
    coordinates: np.ndarray,
    previous_step_index: int,
    previous_transaction_serial: int,
    proposer,
    progress: _Progress,
) -> tdg6.TDG6ShadowPath:
    expected = dict(_LEVEL_COUNTS)[level]
    if len(boundaries) != expected + 1:
        raise ValueError("TDG11 shadow boundaries differ from the level")
    shadow_transaction = _clone_gr0_transaction(transaction)
    shadow_tracers = _clone_tracers(tracers)
    current_state = state
    current_index = previous_step_index
    current_serial = previous_transaction_serial
    attempts: list[Proto7StepAttempt] = []
    for index, (left, right) in enumerate(zip(boundaries, boundaries[1:])):

        def preview(
            proposal: StepProposal, current=current_state, start=left
        ) -> object:
            return shadow_tracers.preview_advance(
                old_state=current,
                new_state=proposal.candidate_state,
                coordinates=coordinates,
                step_size=proposal.final_time - start,
            )

        try:
            attempt = _guarded_shadow_attempt(
                proposer=proposer,
                method=method,
                time=left,
                step_size=right - left,
                state=current_state,
                rhs=rhs,
                projector=projector,
                transaction=shadow_transaction,
                previous_step_index=current_index,
                previous_transaction_serial=current_serial,
                preaccept=preview,
                progress=progress,
            )
        except BaseException as error:
            if isinstance(error, (KeyboardInterrupt, SystemExit)):
                raise
            _raise_premise_stop(
                level=level,
                proposal_index=index,
                progress=progress,
                rhs_calls=rhs.calls,
                cause=error,
            )
        if attempt.retry is not None:
            _raise_premise_stop(
                level=level,
                proposal_index=index,
                progress=progress,
                rhs_calls=rhs.calls,
                source_retry=attempt.retry,
            )
        accepted = attempt.accepted
        if accepted is None or attempt.preaccept_payload is None:
            raise RuntimeError(
                "TDG11 shadow attempt omitted acceptance or tracer preview"
            )
        try:
            shadow_tracers.commit_advance(*attempt.preaccept_payload)
        except BaseException as error:
            if isinstance(error, (KeyboardInterrupt, SystemExit)):
                raise
            _raise_premise_stop(
                level=level,
                proposal_index=index,
                progress=progress,
                rhs_calls=rhs.calls,
                cause=error,
            )
        attempts.append(attempt)
        current_state = accepted.state
        current_index = accepted.step_index
        current_serial = accepted.transaction_serial
    return tdg6.TDG6ShadowPath(
        level=level,
        method=method,
        initial_time=boundaries[0],
        final_time=boundaries[-1],
        initial_state_sha256=_state_hash(state),
        attempts=tuple(attempts),
        final_monitor_state=shadow_transaction.state,
        final_causal_state=shadow_transaction.causal_state,
        final_tracer_positions=np.asarray(shadow_tracers.positions),
        final_tracer_proper_times=np.asarray(shadow_tracers.proper_times),
        final_tracer_snapshot=_tracer_snapshot(shadow_tracers),
    )


def _construct_family(
    *,
    method: str,
    arithmetic_id: str,
    start: float,
    width: float,
    state: EvolutionState,
    rhs: _CountingRHS,
    projector: Callable[[float, EvolutionState], EvolutionState] | None,
    transaction: GR0RuntimeStageTransaction,
    tracers: object,
    coordinates: np.ndarray,
    previous_step_index: int,
    previous_transaction_serial: int,
    owned_row_count: int,
) -> GuardedShadowFamily:
    proposer = _select_proposer(arithmetic_id)
    outer_bounds, medium_bounds, fine_bounds = _family_boundaries(start, width)
    progress = _Progress()
    paths = []
    for level, boundaries in (
        ("outer", outer_bounds),
        ("medium", medium_bounds),
        ("fine", fine_bounds),
    ):
        paths.append(
            _shadow_path(
                level=level,
                method=method,
                boundaries=boundaries,
                state=state,
                rhs=rhs,
                projector=projector,
                transaction=transaction,
                tracers=tracers,
                coordinates=coordinates,
                previous_step_index=previous_step_index,
                previous_transaction_serial=previous_transaction_serial,
                proposer=proposer,
                progress=progress,
            )
        )
    outer, medium, fine = paths
    recorded = validate_recorded_family(
        tuple(
            tuple(attempt.proposal for attempt in path.attempts)
            for path in (outer, medium, fine)
        ),
        method=method,
        arithmetic_id=arithmetic_id,
        expected_owned_row_count=owned_row_count,
    )
    if rhs.calls != progress.prechecks + progress.stages:
        raise RuntimeError(
            "instrumented RHS calls differ from precheck and stage records"
        )
    return GuardedShadowFamily(
        paths=(outer, medium, fine),
        recorded=recorded,
        accepted_source_precheck_count=progress.prechecks,
        stage_record_count=progress.stages,
        rhs_call_count=rhs.calls,
    )


def build_guarded_shadow_family(
    *,
    method: str,
    arithmetic_id: str,
    time: Real,
    step_size: Real,
    state: EvolutionState,
    rhs: Callable[[float, EvolutionState], EvolutionRHS],
    projector: Callable[[float, EvolutionState], EvolutionState] | None,
    transaction: GR0RuntimeStageTransaction,
    tracers: object,
    coordinates: object,
    previous_step_index: int,
    previous_transaction_serial: int,
) -> GuardedShadowFamily:
    """Build one guarded 1/2/4 shadow family without production mutation."""

    if method not in (PRIMARY_METHOD, COMPARATOR_METHOD):
        raise ValueError("unrecognized method")
    if arithmetic_id not in _RECOGNIZED_ARITHMETIC:
        raise ValueError("unrecognized arithmetic")
    if not isinstance(state, EvolutionState) or state.shape[1] != 6:
        raise ValueError("state must contain the six frozen fields")
    if not callable(rhs):
        raise TypeError("rhs must be callable")
    if projector is not None and not callable(projector):
        raise TypeError("projector must be callable or None")
    if not isinstance(transaction, GR0RuntimeStageTransaction):
        raise TypeError("transaction must be GR0RuntimeStageTransaction")
    start = _finite("time", time)
    width = _finite("step_size", step_size, positive=True)
    step_index = _nonnegative_integer("previous_step_index", previous_step_index)
    serial = _nonnegative_integer(
        "previous_transaction_serial", previous_transaction_serial
    )
    if serial != transaction.state.last_transaction_serial + 1:
        raise ValueError("member serial does not follow the transaction ledger")
    if not _same_binary64(transaction.causal_state.accepted_time, start):
        raise ValueError("causal ledger time differs from accepted time")
    if type(coordinates) is not np.ndarray or coordinates.dtype != np.dtype(np.float64):
        raise TypeError("coordinates must be a native binary64 array")
    coordinate_array = _immutable_array("coordinates", coordinates)
    if coordinate_array.size != state.shape[0]:
        raise ValueError("coordinate array differs from state grid")
    _tracer_snapshot(tracers)
    _family_boundaries(start, width)

    before = _input_fingerprint(
        state=state,
        transaction=transaction,
        tracers=tracers,
        coordinates=coordinates,
    )
    counted = _CountingRHS(rhs)
    try:
        family = _construct_family(
            method=method,
            arithmetic_id=arithmetic_id,
            start=start,
            width=width,
            state=state,
            rhs=counted,
            projector=projector,
            transaction=transaction,
            tracers=tracers,
            coordinates=coordinate_array,
            previous_step_index=step_index,
            previous_transaction_serial=serial,
            owned_row_count=state.shape[0] - 5,
        )
    except BaseException as error:
        after = _input_fingerprint(
            state=state,
            transaction=transaction,
            tracers=tracers,
            coordinates=coordinates,
        )
        if after != before:
            raise TDG11ShadowIntegrityError(
                "accepted diagnostic inputs mutated during guarded shadow family construction"
            ) from error
        raise
    after = _input_fingerprint(
        state=state,
        transaction=transaction,
        tracers=tracers,
        coordinates=coordinates,
    )
    if after != before:
        raise TDG11ShadowIntegrityError(
            "accepted diagnostic inputs mutated during guarded shadow family construction"
        )
    return family


__all__ = [
    "GuardedShadowFamily",
    "TDG11ShadowIntegrityError",
    "TDG11ShadowPremiseStop",
    "build_guarded_shadow_family",
]
