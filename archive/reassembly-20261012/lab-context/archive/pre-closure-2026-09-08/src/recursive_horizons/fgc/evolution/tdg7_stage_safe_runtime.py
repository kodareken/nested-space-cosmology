"""TDG7 stage-safe adapter for TDG6's immutable temporal compositor.

TDG7 owns only coordinate selection and its proof that the *actual* RK4 and
SSPRK3 abscissae are exact binary64 points.  TDG6 remains the sole owner of
shadows, admission, temporal retry evidence, and atomic fine-path commit.
This module deliberately does not reproduce numerical-engine stages.
"""

from __future__ import annotations

from dataclasses import dataclass
import json
from math import isfinite
from numbers import Real
from typing import Any, Callable

import numpy as np

from . import tdg6_temporal_admission_runtime as tdg6
from .numerical_engine import (
    COMPARATOR_METHOD,
    PRIMARY_METHOD,
    EvolutionRHS,
    EvolutionState,
    array_content_sha256,
)
from .proto5_runtime import GR0RuntimeStageTransaction
from .tdg6_temporal_admission_runtime import _shadow_path
from .tdg6_temporal_admission_design import (
    TDG6_MAXIMUM_TEMPORAL_RETRIES_PER_STEP,
    TDG6_MINIMUM_MACRO_STEP,
)
from .tdg7_binary64_subdivision_lattice import (
    CoordinateLatticeLimitReached,
    TDG7Binary64SubdivisionPlan,
    plan_forward_proto14_subdivision,
    validate_tdg7_binary64_subdivision_plan,
)


__all__ = (
    "FRESH_VS_RETRY_PROTOCOL_DISCRIMINATOR_IMPLEMENTED",
    "STATELESS_ADAPTER_CANNOT_DETECT_DISCARDED_EXTERNAL_RETRY",
    "TDG7CoordinateLatticeStop",
    "TDG7PreparedStageSafeRuntime",
    "TDG7RetryLatticeExhausted",
    "TDG7RetrySuccessor",
    "commit_tdg7_stage_safe_runtime",
    "prepare_tdg7_initial_stage_safe_runtime",
    "prepare_tdg7_retry_successor",
    "replan_tdg7_after_tdg6_retry",
    "require_tdg7_stage_safe_admission",
    "validate_tdg7_stage_safe_prepared",
)


# IMP3 is deliberately a stateless synthetic primitive.  A future frozen
# protocol/runtime cursor must decide whether an externally durable boundary is
# fresh or retrying; identical rewound inputs contain no information about a
# discarded prior rejection.  The disjoint public entrypoints below prevent an
# accidental optional-argument bypass without pretending to solve that later
# lifecycle problem here.
FRESH_VS_RETRY_PROTOCOL_DISCRIMINATOR_IMPLEMENTED = False
STATELESS_ADAPTER_CANNOT_DETECT_DISCARDED_EXTERNAL_RETRY = True


def _same_binary64(left: float, right: float) -> bool:
    return np.float64(left).tobytes() == np.float64(right).tobytes()


def _finite(name: str, value: object) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a finite real scalar")
    result = float(value)
    if not isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


class TDG7CoordinateLatticeStop(CoordinateLatticeLimitReached):
    """A nonphysical TDG7 stop raised before any runtime hook is touched."""

    physical_classification = False

    def __init__(self, reason: str) -> None:
        super().__init__(reason)


def _plan_forward_with_minimum(
    time: Real, event_target: Real, cap: Real, minimum_width: Real
) -> TDG7Binary64SubdivisionPlan:
    try:
        return plan_forward_proto14_subdivision(
            time, event_target, cap, minimum_width=minimum_width
        )
    except CoordinateLatticeLimitReached as error:
        raise TDG7CoordinateLatticeStop(error.reason) from error


def _tracer_snapshot(tracers: object) -> tdg6.TDG6TracerSnapshot:
    """Construct TDG6's public immutable boundary receipt without mutation."""

    required = (
        "labels",
        "positions",
        "proper_times",
        "event_proper_times",
        "event_fields",
        "cutoff",
        "outer_radius",
        "preview_advance",
        "commit_advance",
    )
    if any(not hasattr(tracers, item) for item in required):
        raise TypeError("tracers do not expose the normal-flow runtime contract")

    def array(name: str, value: object, ndim: int = 1) -> np.ndarray:
        answer = np.ascontiguousarray(np.asarray(value, dtype=np.float64))
        if answer.ndim != ndim or np.any(~np.isfinite(answer)):
            raise ValueError(f"{name} must be a finite {ndim}-dimensional array")
        return answer

    labels = array("tracer labels", tracers.labels)
    positions = array("tracer positions", tracers.positions)
    proper = array("tracer proper times", tracers.proper_times)
    if not (labels.shape == positions.shape == proper.shape):
        raise ValueError("tracer arrays have different shapes")
    event_times = tuple(
        array("event proper times", item) for item in tracers.event_proper_times
    )
    event_fields = tuple(
        array("event fields", item, 2) for item in tracers.event_fields
    )
    if len(event_times) != len(event_fields):
        raise ValueError("tracer event histories have different lengths")
    history = event_times + event_fields or (np.zeros((0,), dtype=np.float64),)
    return tdg6.TDG6TracerSnapshot(
        labels_sha256=array_content_sha256(labels),
        positions_sha256=array_content_sha256(positions),
        proper_times_sha256=array_content_sha256(proper),
        event_history_sha256=array_content_sha256(*history),
        event_count=len(event_times),
        cutoff=_finite("tracer cutoff", tracers.cutoff),
        outer_radius=_finite("tracer outer_radius", tracers.outer_radius),
    )


def _expected_stages(
    method: str, left: float, right: float
) -> tuple[tuple[str, float], ...]:
    midpoint = left + (right - left) / 2.0
    if method == PRIMARY_METHOD:
        return (
            ("rk4_k1", left),
            ("rk4_k2", midpoint),
            ("rk4_k3", midpoint),
            ("rk4_k4", right),
            ("candidate_endpoint", right),
        )
    if method == COMPARATOR_METHOD:
        return (
            ("ssprk3_s0", left),
            ("ssprk3_s1", right),
            ("ssprk3_s2", midpoint),
            ("candidate_endpoint", right),
        )
    raise ValueError("unknown TDG7 runtime method")


def _validate_path(
    plan: TDG7Binary64SubdivisionPlan, path: tdg6.TDG6ShadowPath
) -> None:
    expected_boundaries = {
        "outer": plan.one_full,
        "medium": plan.two_half,
        "fine": plan.four_quarter,
    }[path.level]
    if len(path.attempts) != len(expected_boundaries) - 1:
        raise RuntimeError("TDG7 proposal count differs from the shared plan")
    for attempt, (left, right) in zip(
        path.attempts, zip(expected_boundaries, expected_boundaries[1:]), strict=True
    ):
        proposal = attempt.proposal
        if not (
            _same_binary64(proposal.initial_time, left)
            and _same_binary64(proposal.final_time, right)
        ):
            raise RuntimeError("TDG7 proposal boundary differs from the shared plan")
        records = proposal.stages
        expected = _expected_stages(path.method, left, right)
        if len(records) != len(expected) or any(
            record.stage_name != name or not _same_binary64(record.time, time)
            for record, (name, time) in zip(records, expected, strict=True)
        ):
            raise RuntimeError("TDG7 stage record differs from the exact shared plan")
        if attempt.accepted is None or not _same_binary64(attempt.accepted.time, right):
            raise RuntimeError(
                "TDG7 accepted proposal endpoint differs from the shared plan"
            )
        if len(attempt.accepted.stage_records) != len(expected) or any(
            record.stage_name != name or not _same_binary64(record.time, time)
            for record, (name, time) in zip(
                attempt.accepted.stage_records, expected, strict=True
            )
        ):
            raise RuntimeError(
                "TDG7 accepted stage record differs from the exact shared plan"
            )


def validate_tdg7_stage_safe_prepared(prepared: "TDG7PreparedStageSafeRuntime") -> None:
    """Recheck the frozen plan against all retained TDG6 proposal records."""

    if not isinstance(prepared, TDG7PreparedStageSafeRuntime):
        raise TypeError("prepared must be TDG7PreparedStageSafeRuntime")
    validate_tdg7_binary64_subdivision_plan(prepared.plan)
    if prepared.plan.minimum_width is None or not _same_binary64(
        prepared.plan.minimum_width, float(TDG6_MINIMUM_MACRO_STEP)
    ):
        raise RuntimeError("TDG7 runtime plan changed the frozen TDG6 minimum width")
    inner = prepared.prepared_tdg6
    if not isinstance(inner, tdg6.TDG6PreparedGR0Compositor):
        raise TypeError("TDG7 wrapper omits TDG6 prepared compositor")
    if not (
        _same_binary64(inner.initial_time, prepared.plan.current)
        and _same_binary64(inner.final_time, prepared.plan.endpoint)
    ):
        raise RuntimeError("TDG7 prepared TDG6 bounds differ from the plan")
    for path in (inner.outer, inner.medium, inner.fine):
        _validate_path(prepared.plan, path)


@dataclass(frozen=True, slots=True)
class TDG7PreparedStageSafeRuntime:
    """The single immutable plan and the TDG6 preparation derived from it."""

    plan: TDG7Binary64SubdivisionPlan
    prepared_tdg6: tdg6.TDG6PreparedGR0Compositor

    def __post_init__(self) -> None:
        validate_tdg7_stage_safe_prepared(self)


class TDG7RetryLatticeExhausted(RuntimeError):
    """A TDG6 numerical retry cap cannot make a smaller TDG7-safe plan."""

    physical_classification = False

    def __init__(
        self, *, retry: tdg6.TDG6TemporalRetryRequired, lattice_reason: str
    ) -> None:
        self.evidence = retry.evidence
        self.updated_ledger = retry.updated_ledger
        self.requested_half_cap = retry.retry_step_size
        self.lattice_reason = lattice_reason
        super().__init__(f"TDG7 retry lattice exhausted: {lattice_reason}")


def _canonical_rejection(evidence: tdg6.TDG6TemporalRetryEvidence) -> str:
    return json.dumps(
        evidence.as_mapping(), sort_keys=True, separators=(",", ":"), allow_nan=False
    )


@dataclass(frozen=True, slots=True)
class TDG7RetrySuccessor:
    """One immutable TDG6-rejection lineage and its TDG7-safe next plan."""

    predecessor: "TDG7PreparedStageSafeRuntime"
    plan: TDG7Binary64SubdivisionPlan
    retry: tdg6.TDG6TemporalRetryRequired

    def __post_init__(self) -> None:
        if not isinstance(self.predecessor, TDG7PreparedStageSafeRuntime):
            raise TypeError("predecessor must be TDG7PreparedStageSafeRuntime")
        validate_tdg7_stage_safe_prepared(self.predecessor)
        predecessor_plan = self.predecessor.plan
        inner = self.predecessor.prepared_tdg6
        validate_tdg7_binary64_subdivision_plan(self.plan)
        if not isinstance(self.retry, tdg6.TDG6TemporalRetryRequired):
            raise TypeError("retry must be TDG6TemporalRetryRequired")
        evidence = self.retry.evidence
        ledger = self.retry.updated_ledger
        if (
            evidence.retry_count_for_current_macro_step
            > TDG6_MAXIMUM_TEMPORAL_RETRIES_PER_STEP
        ):
            raise ValueError("TDG7 retry lineage exceeds TDG6 retry limit")
        if not (
            evidence.method == inner.method
            and evidence.initial_state_sha256 == inner.initial_state_sha256
            and evidence.previous_step_index == inner.previous_step_index
            and evidence.previous_transaction_serial
            == inner.previous_transaction_serial
            and _same_binary64(evidence.initial_time, predecessor_plan.current)
            and _same_binary64(
                evidence.attempted_macro_step_size, predecessor_plan.macro_width
            )
            and _same_binary64(
                evidence.retry_step_size, predecessor_plan.macro_width / 2.0
            )
            and _same_binary64(self.plan.current, predecessor_plan.current)
            and _same_binary64(self.plan.requested_cap, evidence.retry_step_size)
            and self.plan.macro_width < predecessor_plan.macro_width
        ):
            raise ValueError("TDG7 retry successor differs from TDG6 retry evidence")
        original = inner.original_temporal_ledger
        expected_serialized = (
            *original.serialized_temporal_rejections,
            _canonical_rejection(evidence),
        )
        if (
            ledger.current_macro_step_temporal_retry_count
            != evidence.retry_count_for_current_macro_step
            or ledger.cumulative_temporal_retry_count
            != evidence.cumulative_temporal_retry_count
            or not _same_binary64(
                ledger.last_accepted_time,
                original.last_accepted_time,
            )
            or ledger.accepted_macro_step_count != original.accepted_macro_step_count
            or ledger.accumulated_debit_vector != original.accumulated_debit_vector
            or ledger.last_accepted_macro_step_temporal_retry_count
            != original.last_accepted_macro_step_temporal_retry_count
            or ledger.cumulative_temporal_retry_count
            != original.cumulative_temporal_retry_count + 1
            or ledger.current_macro_step_temporal_retry_count
            != original.current_macro_step_temporal_retry_count + 1
            or ledger.serialized_temporal_rejections != expected_serialized
        ):
            raise ValueError(
                "TDG7 retry successor omits canonical durable TDG6 evidence"
            )

    @property
    def temporal_ledger(self) -> tdg6.TDG6TemporalLedger:
        return self.retry.updated_ledger


def _prepare_tdg7_stage_safe_runtime(
    *,
    method: str,
    time: Real,
    event_target: Real,
    requested_cap: Real,
    state: EvolutionState,
    rhs: Callable[[float, EvolutionState], EvolutionRHS],
    projector: Callable[[float, EvolutionState], EvolutionState] | None,
    transaction: GR0RuntimeStageTransaction,
    tracers: object,
    coordinates: object,
    temporal_ledger: tdg6.TDG6TemporalLedger,
    previous_step_index: int,
    previous_transaction_serial: int,
    retry_successor: TDG7RetrySuccessor | None = None,
) -> TDG7PreparedStageSafeRuntime:
    """Private shared implementation for the disjoint initial/retry entrypoints."""

    # This must be the first operation: a lattice stop cannot invoke a source,
    # projector, tracer, transaction, or state hook.
    plan = _plan_forward_with_minimum(
        time, event_target, requested_cap, float(TDG6_MINIMUM_MACRO_STEP)
    )
    if retry_successor is not None:
        if not isinstance(retry_successor, TDG7RetrySuccessor):
            raise TypeError("retry_successor must be TDG7RetrySuccessor")
        if (
            plan != retry_successor.plan
            or temporal_ledger != retry_successor.temporal_ledger
        ):
            raise ValueError("TDG7 retry successor plan or durable ledger is stale")
    if method not in (PRIMARY_METHOD, COMPARATOR_METHOD):
        raise ValueError("unknown TDG7 runtime method")
    if not isinstance(state, EvolutionState) or state.shape[1] != 6:
        raise ValueError("TDG7 state must contain the six frozen fields")
    if not callable(rhs) or (projector is not None and not callable(projector)):
        raise TypeError("TDG7 rhs/projector contract is invalid")
    if not isinstance(transaction, GR0RuntimeStageTransaction):
        raise TypeError("transaction must be GR0RuntimeStageTransaction")
    if not isinstance(temporal_ledger, tdg6.TDG6TemporalLedger):
        raise TypeError("temporal_ledger must be TDG6TemporalLedger")
    if (
        retry_successor is None
        and temporal_ledger.current_macro_step_temporal_retry_count != 0
    ):
        raise ValueError(
            "TDG7 initial entrypoint cannot consume an active retry ledger"
        )
    coordinate_array = np.ascontiguousarray(np.asarray(coordinates, dtype=np.float64))
    if (
        coordinate_array.ndim != 1
        or coordinate_array.size != state.shape[0]
        or np.any(~np.isfinite(coordinate_array))
    ):
        raise ValueError("TDG7 coordinate array differs from state grid")
    if previous_step_index < 0 or previous_transaction_serial < 0:
        raise ValueError("TDG7 previous indices must be nonnegative")
    if previous_transaction_serial != transaction.state.last_transaction_serial + 1:
        raise ValueError("TDG7 member serial does not follow the transaction ledger")
    if not (
        _same_binary64(transaction.causal_state.accepted_time, plan.current)
        and _same_binary64(temporal_ledger.last_accepted_time, plan.current)
    ):
        raise ValueError("TDG7 accepted ledgers differ from planned time")
    if retry_successor is not None:
        evidence = retry_successor.retry.evidence
        if (
            evidence.method != method
            or evidence.initial_state_sha256
            != array_content_sha256(state.u, state.p, state.q)
            or evidence.previous_step_index != previous_step_index
            or evidence.previous_transaction_serial != previous_transaction_serial
        ):
            raise ValueError(
                "TDG7 retry successor differs from accepted runtime boundary"
            )

    original_monitor, original_causal = transaction.state, transaction.causal_state
    original_tracer, original_ledger = _tracer_snapshot(tracers), temporal_ledger
    state_hash = array_content_sha256(state.u, state.p, state.q)
    outer = _shadow_path(
        level="outer",
        method=method,
        boundaries=plan.one_full,
        state=state,
        rhs=rhs,
        projector=projector,
        transaction=transaction,
        tracers=tracers,
        coordinates=coordinate_array,
        previous_step_index=previous_step_index,
        previous_transaction_serial=previous_transaction_serial,
    )
    medium = _shadow_path(
        level="medium",
        method=method,
        boundaries=plan.two_half,
        state=state,
        rhs=rhs,
        projector=projector,
        transaction=transaction,
        tracers=tracers,
        coordinates=coordinate_array,
        previous_step_index=previous_step_index,
        previous_transaction_serial=previous_transaction_serial,
    )
    fine = _shadow_path(
        level="fine",
        method=method,
        boundaries=plan.four_quarter,
        state=state,
        rhs=rhs,
        projector=projector,
        transaction=transaction,
        tracers=tracers,
        coordinates=coordinate_array,
        previous_step_index=previous_step_index,
        previous_transaction_serial=previous_transaction_serial,
    )
    if (
        transaction.state != original_monitor
        or transaction.causal_state != original_causal
        or _tracer_snapshot(tracers) != original_tracer
        or array_content_sha256(state.u, state.p, state.q) != state_hash
        or temporal_ledger != original_ledger
    ):
        raise RuntimeError("TDG7 shadow preparation mutated a real boundary")
    inner = tdg6.TDG6PreparedGR0Compositor(
        method=method,
        initial_time=plan.current,
        final_time=plan.endpoint,
        initial_state_sha256=state_hash,
        previous_step_index=previous_step_index,
        previous_transaction_serial=previous_transaction_serial,
        original_monitor_state=original_monitor,
        original_causal_state=original_causal,
        original_tracer_snapshot=original_tracer,
        original_temporal_ledger=original_ledger,
        outer=outer,
        medium=medium,
        fine=fine,
        continuous_admission=tdg6.assess_tdg6_continuous_admission(
            outer=outer, medium=medium, fine=fine
        ),
    )
    return TDG7PreparedStageSafeRuntime(plan, inner)


def prepare_tdg7_initial_stage_safe_runtime(
    *,
    method: str,
    time: Real,
    event_target: Real,
    requested_cap: Real,
    state: EvolutionState,
    rhs: Callable[[float, EvolutionState], EvolutionRHS],
    projector: Callable[[float, EvolutionState], EvolutionState] | None,
    transaction: GR0RuntimeStageTransaction,
    tracers: object,
    coordinates: object,
    temporal_ledger: tdg6.TDG6TemporalLedger,
    previous_step_index: int,
    previous_transaction_serial: int,
) -> TDG7PreparedStageSafeRuntime:
    """Prepare a declared fresh macro attempt under the frozen TDG6 minimum.

    Choosing whether a durable external cursor is fresh or retrying belongs to
    a future protocol/runtime gate.  This entrypoint rejects a ledger which
    already declares an active retry; the disjoint retry entrypoint below is
    the only public route which consumes :class:`TDG7RetrySuccessor`.
    """

    return _prepare_tdg7_stage_safe_runtime(
        method=method,
        time=time,
        event_target=event_target,
        requested_cap=requested_cap,
        state=state,
        rhs=rhs,
        projector=projector,
        transaction=transaction,
        tracers=tracers,
        coordinates=coordinates,
        temporal_ledger=temporal_ledger,
        previous_step_index=previous_step_index,
        previous_transaction_serial=previous_transaction_serial,
        retry_successor=None,
    )


def require_tdg7_stage_safe_admission(
    prepared: TDG7PreparedStageSafeRuntime, **kwargs: Any
) -> None:
    validate_tdg7_stage_safe_prepared(prepared)
    tdg6.require_tdg6_temporal_admission(prepared.prepared_tdg6, **kwargs)


def commit_tdg7_stage_safe_runtime(
    prepared: TDG7PreparedStageSafeRuntime, **kwargs: Any
) -> tdg6.TDG6CommittedGR0Compositor:
    validate_tdg7_stage_safe_prepared(prepared)
    return tdg6.commit_tdg6_gr0_compositor(prepared.prepared_tdg6, **kwargs)


def replan_tdg7_after_tdg6_retry(
    prepared: TDG7PreparedStageSafeRuntime,
    retry: tdg6.TDG6TemporalRetryRequired,
) -> TDG7RetrySuccessor:
    """Consume TDG6's halved cap through TDG7; never return an unsafe raw cap."""

    validate_tdg7_stage_safe_prepared(prepared)
    plan = prepared.plan
    if not isinstance(retry, tdg6.TDG6TemporalRetryRequired):
        raise TypeError("retry must be TDG6TemporalRetryRequired")
    evidence = retry.evidence
    ledger = retry.updated_ledger
    inner = prepared.prepared_tdg6
    if (
        evidence.retry_count_for_current_macro_step
        > TDG6_MAXIMUM_TEMPORAL_RETRIES_PER_STEP
    ):
        raise ValueError("TDG7 retry lineage exceeds TDG6 retry limit")
    if (
        evidence.method != inner.method
        or evidence.initial_state_sha256 != inner.initial_state_sha256
        or evidence.previous_step_index != inner.previous_step_index
        or evidence.previous_transaction_serial != inner.previous_transaction_serial
        or ledger.current_macro_step_temporal_retry_count
        != evidence.retry_count_for_current_macro_step
        or ledger.cumulative_temporal_retry_count
        != evidence.cumulative_temporal_retry_count
        or not ledger.serialized_temporal_rejections
        or ledger.serialized_temporal_rejections[-1] != _canonical_rejection(evidence)
    ):
        raise ValueError("TDG7 retry lineage differs from durable TDG6 evidence")
    if not _same_binary64(retry.evidence.initial_time, plan.current):
        raise TDG7RetryLatticeExhausted(
            retry=retry, lattice_reason="retry_initial_time_differs_from_plan"
        )
    if not _same_binary64(retry.retry_step_size, plan.macro_width / 2.0):
        raise TDG7RetryLatticeExhausted(
            retry=retry, lattice_reason="retry_cap_differs_from_TDG6_half_step"
        )
    try:
        next_plan = _plan_forward_with_minimum(
            plan.current,
            plan.event_target,
            retry.retry_step_size,
            float(TDG6_MINIMUM_MACRO_STEP),
        )
    except TDG7CoordinateLatticeStop as error:
        raise TDG7RetryLatticeExhausted(
            retry=retry, lattice_reason=error.reason
        ) from error
    if not next_plan.macro_width < plan.macro_width:
        raise TDG7RetryLatticeExhausted(
            retry=retry, lattice_reason="retry_did_not_strictly_reduce_macro_width"
        )
    return TDG7RetrySuccessor(
        predecessor=prepared,
        plan=next_plan,
        retry=retry,
    )


def prepare_tdg7_retry_successor(
    *,
    successor: TDG7RetrySuccessor,
    method: str,
    state: EvolutionState,
    rhs: Callable[[float, EvolutionState], EvolutionRHS],
    projector: Callable[[float, EvolutionState], EvolutionState] | None,
    transaction: GR0RuntimeStageTransaction,
    tracers: object,
    coordinates: object,
    previous_step_index: int,
    previous_transaction_serial: int,
) -> TDG7PreparedStageSafeRuntime:
    """Prepare only with the durable ledger and exact cap bound to *successor*."""

    if not isinstance(successor, TDG7RetrySuccessor):
        raise TypeError("successor must be TDG7RetrySuccessor")
    return _prepare_tdg7_stage_safe_runtime(
        method=method,
        time=successor.plan.current,
        event_target=successor.plan.event_target,
        requested_cap=successor.plan.requested_cap,
        state=state,
        rhs=rhs,
        projector=projector,
        transaction=transaction,
        tracers=tracers,
        coordinates=coordinates,
        temporal_ledger=successor.temporal_ledger,
        previous_step_index=previous_step_index,
        previous_transaction_serial=previous_transaction_serial,
        retry_successor=successor,
    )
