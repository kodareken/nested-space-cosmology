"""Production-shaped TDG6 three-level temporal compositor.

The compositor implements the numerical object authorized by
``FGC-1-TDG6-PREF21`` without advancing a production campaign.  One full,
two half, and four quarter steps are evaluated from the same bitwise accepted
state through independent PROTO7 monitor, causal, and normal-flow-tracer
shadows.  The outer and medium paths are evidence only.  Only a complete,
admitted four-quarter path can be committed.

For every owned ``(u,p,q) x (alpha,v,lambda,R,phi,chi)`` channel, continuous
cubic-Hermite differences are enclosed for ``D01`` and ``D12``.  The frozen
exact ``8 U12^2 <= L01^2`` rule classifies the channel, and the complete
``U12`` vector is accumulated by outward nonnegative addition.  Temporal
retries have their own durable ledger and never become physical outcomes.

This is a production-shaped instrument plus synthetic qualification surface.
It does not authorize a trajectory, PROTO14, an eligible GR-0 case, SGB-L,
FGC-QR, DEF1, retained-EFT validity, a transition, or any physical claim.
"""

from __future__ import annotations

from copy import copy
from dataclasses import asdict, dataclass
from fractions import Fraction
import json
from math import isfinite, nextafter
from numbers import Real
from typing import Any, Callable, Mapping, Sequence

import numpy as np

from .boundary_domain import CausalBudgetState
from .numerical_engine import (
    EvolutionRHS,
    EvolutionState,
    array_content_sha256,
)
from .proto5_runtime import GR0RuntimeMonitorState, GR0RuntimeStageTransaction
from .proto7_runtime import ProposalFailureEvidence, Proto7StepAttempt, attempt_proto7_step
from . import tdg5_stage_complete_refinement_runtime as tdg5
from .tdg6_temporal_admission_design import (
    CertifiedMagnitudeInterval,
    TDG6_COMPLETE_STATE_CHANNELS,
    TDG6_MAXIMUM_TEMPORAL_RETRIES_PER_STEP,
    TDG6_MINIMUM_MACRO_STEP,
    TDG6_PASSING_CLASSES,
    TDG6_RETRY_FACTOR,
    classify_tdg6_channel,
)


TDG6_RUNTIME_SCHEMA_VERSION = 1
TDG6_LEVEL_STEP_COUNTS = {"outer": 1, "medium": 2, "fine": 4}
TDG6_PATH_NAMES = tuple(
    ["outer_0"]
    + [f"medium_{index}" for index in range(2)]
    + [f"fine_{index}" for index in range(4)]
)
_FIELD_NAMES = ("alpha", "v", "lambda", "R", "phi", "chi")


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


def _canonical_json(value: Mapping[str, object]) -> str:
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    )


def _outward_nonnegative_sum(left: float, right: float) -> float:
    """Add two binary64 upper bounds without rounding the result downward."""

    a = _finite("left debit", left)
    b = _finite("right debit", right)
    if a < 0.0 or b < 0.0:
        raise ValueError("temporal debits must be nonnegative")
    exact = Fraction.from_float(a) + Fraction.from_float(b)
    nominal = float(exact)
    if not isfinite(nominal):
        raise FloatingPointError("temporal debit accumulation overflowed")
    if Fraction.from_float(nominal) < exact:
        nominal = nextafter(nominal, float("inf"))
    return nominal


@dataclass(frozen=True, slots=True)
class TDG6Binary64MagnitudeInterval:
    """Outward binary64 interval for one channel's continuous difference."""

    lower_bound: float
    upper_bound: float
    subinterval_count: int
    owned_row_count: int
    envelope: tdg5.Binary64CubicEnvelope

    def __post_init__(self) -> None:
        lower = _finite("lower_bound", self.lower_bound)
        upper = _finite("upper_bound", self.upper_bound)
        if lower < 0.0 or upper < lower:
            raise ValueError("TDG6 magnitude interval is not ordered and nonnegative")
        if (
            isinstance(self.subinterval_count, bool)
            or self.subinterval_count not in (2, 4)
        ):
            raise ValueError("TDG6 interval must own two or four subintervals")
        if self.owned_row_count < 1:
            raise ValueError("TDG6 interval must own at least one radial row")
        if not isinstance(self.envelope, tdg5.Binary64CubicEnvelope):
            raise TypeError("envelope must be Binary64CubicEnvelope")
        if self.envelope.polynomial_count != self.subinterval_count * self.owned_row_count:
            raise ValueError("TDG6 envelope omits a row or subinterval")
        if upper != self.envelope.certified_continuous_upper_bound:
            raise ValueError("TDG6 upper bound differs from the certified envelope")

    def exact_interval(self) -> CertifiedMagnitudeInterval:
        return CertifiedMagnitudeInterval(
            Fraction.from_float(self.lower_bound),
            Fraction.from_float(self.upper_bound),
        )


def _magnitude_interval(
    coefficients: np.ndarray,
    radii: np.ndarray,
    *,
    subinterval_count: int,
    owned_row_count: int,
) -> TDG6Binary64MagnitudeInterval:
    values = np.asarray(coefficients, dtype=np.float64).reshape((-1, 4))
    uncertainty = np.asarray(radii, dtype=np.float64).reshape((-1, 4))
    envelope = tdg5.binary64_cubic_absolute_envelope(values, uncertainty)

    # The nominal value at every evaluated candidate is a legitimate point in
    # [0,1], even when the candidate is only an approximate stationary point.
    # Subtracting the complete coefficient and evaluation debits therefore
    # gives a conservative lower bound on the true supremum.  Fraction
    # arithmetic prevents the subtraction itself from rounding upward.
    exact_lower = max(
        Fraction(0),
        Fraction.from_float(envelope.raw_candidate_maximum)
        - Fraction.from_float(envelope.coefficient_construction_debit)
        - Fraction.from_float(envelope.outward_arithmetic_debit),
    )
    lower = float(exact_lower)
    if Fraction.from_float(lower) > exact_lower:
        lower = nextafter(lower, 0.0)
    return TDG6Binary64MagnitudeInterval(
        lower_bound=lower,
        upper_bound=envelope.certified_continuous_upper_bound,
        subinterval_count=subinterval_count,
        owned_row_count=owned_row_count,
        envelope=envelope,
    )


@dataclass(frozen=True, slots=True)
class TDG6RuntimeChannelAdmission:
    channel: str
    outer_difference: TDG6Binary64MagnitudeInterval
    finest_difference: TDG6Binary64MagnitudeInterval
    classification: str
    admission_passed: bool
    temporal_retry_permitted: bool
    order_threshold_resolved: bool
    order_threshold_passed: bool | None
    finest_pair_debit: float

    def __post_init__(self) -> None:
        if self.channel not in TDG6_COMPLETE_STATE_CHANNELS:
            raise ValueError("unknown TDG6 channel")
        if self.classification not in TDG6_PASSING_CLASSES | {
            "resolved_order_failure",
            "order_inconclusive",
        }:
            raise ValueError("unknown TDG6 runtime classification")
        if self.admission_passed is not (self.classification in TDG6_PASSING_CLASSES):
            raise ValueError("TDG6 runtime pass differs from classification")
        if self.temporal_retry_permitted is not (not self.admission_passed):
            raise ValueError("TDG6 runtime retry differs from admission")
        if self.finest_pair_debit != self.finest_difference.upper_bound:
            raise ValueError("TDG6 runtime debit differs from D12 upper bound")
        decision = classify_tdg6_channel(
            self.outer_difference.exact_interval(),
            self.finest_difference.exact_interval(),
        )
        if (
            decision.classification != self.classification
            or decision.admission_passed is not self.admission_passed
            or decision.temporal_retry_permitted is not self.temporal_retry_permitted
            or decision.order_threshold_resolved is not self.order_threshold_resolved
            or decision.order_threshold_passed is not self.order_threshold_passed
        ):
            raise ValueError("TDG6 runtime channel differs from exact frozen classifier")


@dataclass(frozen=True, slots=True)
class TDG6ContinuousAdmissionEvidence:
    method: str
    owned_row_count: int
    channel_admissions: tuple[TDG6RuntimeChannelAdmission, ...]
    complete_admission_passed: bool
    failed_channels: tuple[str, ...]
    finest_pair_debit_vector: tuple[float, ...]
    outer_subinterval_count: int = 2
    finest_subinterval_count: int = 4
    admission_is_all_of: bool = True
    debit_cancellation_permitted: bool = False
    absolute_state_tolerance_used: bool = False
    physical_signal_used_for_normalization: bool = False
    declared_cubic_is_exact_PDE_history: bool = False
    trajectory_asymptotic_regime_proved: bool = False
    rigorous_global_PDE_error_bound_proved: bool = False

    def __post_init__(self) -> None:
        if self.method not in tdg5.TDG5_RUNTIME_METHODS:
            raise ValueError("unknown TDG6 method")
        if self.owned_row_count < 1:
            raise ValueError("owned_row_count must be positive")
        channels = tuple(item.channel for item in self.channel_admissions)
        if channels != TDG6_COMPLETE_STATE_CHANNELS:
            raise ValueError("TDG6 channel evidence is incomplete or reordered")
        expected_failed = tuple(
            item.channel for item in self.channel_admissions if not item.admission_passed
        )
        if expected_failed != self.failed_channels:
            raise ValueError("TDG6 failed-channel ledger differs from channel evidence")
        if self.complete_admission_passed is not (not expected_failed):
            raise ValueError("TDG6 all-of admission differs from channel evidence")
        expected_debit = tuple(item.finest_pair_debit for item in self.channel_admissions)
        if expected_debit != self.finest_pair_debit_vector:
            raise ValueError("TDG6 debit vector differs from channel evidence")
        if (
            self.outer_subinterval_count != 2
            or self.finest_subinterval_count != 4
            or self.admission_is_all_of is not True
            or self.debit_cancellation_permitted is not False
            or self.absolute_state_tolerance_used is not False
            or self.physical_signal_used_for_normalization is not False
            or self.declared_cubic_is_exact_PDE_history is not False
            or self.trajectory_asymptotic_regime_proved is not False
            or self.rigorous_global_PDE_error_bound_proved is not False
        ):
            raise ValueError("TDG6 continuous evidence crossed its scientific boundary")


def classify_tdg6_runtime_intervals(
    *,
    method: str,
    outer_intervals: Mapping[str, TDG6Binary64MagnitudeInterval],
    finest_intervals: Mapping[str, TDG6Binary64MagnitudeInterval],
) -> TDG6ContinuousAdmissionEvidence:
    """Apply the exact frozen classifier to a complete ordered runtime map."""

    if method not in tdg5.TDG5_RUNTIME_METHODS:
        raise ValueError("unknown TDG6 method")
    if tuple(outer_intervals) != TDG6_COMPLETE_STATE_CHANNELS:
        raise ValueError("outer interval map is incomplete or reordered")
    if tuple(finest_intervals) != TDG6_COMPLETE_STATE_CHANNELS:
        raise ValueError("finest interval map is incomplete or reordered")
    admissions: list[TDG6RuntimeChannelAdmission] = []
    for channel in TDG6_COMPLETE_STATE_CHANNELS:
        outer = outer_intervals[channel]
        finest = finest_intervals[channel]
        if outer.owned_row_count != finest.owned_row_count:
            raise ValueError("TDG6 D01 and D12 owned-row counts differ")
        decision = classify_tdg6_channel(
            outer.exact_interval(), finest.exact_interval()
        )
        admissions.append(
            TDG6RuntimeChannelAdmission(
                channel=channel,
                outer_difference=outer,
                finest_difference=finest,
                classification=decision.classification,
                admission_passed=decision.admission_passed,
                temporal_retry_permitted=decision.temporal_retry_permitted,
                order_threshold_resolved=decision.order_threshold_resolved,
                order_threshold_passed=decision.order_threshold_passed,
                finest_pair_debit=finest.upper_bound,
            )
        )
    failed = tuple(item.channel for item in admissions if not item.admission_passed)
    return TDG6ContinuousAdmissionEvidence(
        method=method,
        owned_row_count=admissions[0].outer_difference.owned_row_count,
        channel_admissions=tuple(admissions),
        complete_admission_passed=not failed,
        failed_channels=failed,
        finest_pair_debit_vector=tuple(item.finest_pair_debit for item in admissions),
    )


def _proposal_coefficients(
    proposal: Any,
) -> tuple[np.ndarray, np.ndarray]:
    start_rhs, end_rhs = tdg5._proposal_endpoint_records(proposal)
    return tdg5._hermite_coefficients_with_radius(
        tdg5._owned_state(proposal.initial_state),
        tdg5._owned_rhs(start_rhs),
        tdg5._owned_state(proposal.candidate_state),
        tdg5._owned_rhs(end_rhs),
        width=proposal.final_time - proposal.initial_time,
    )


def _path_proposals(path: "TDG6ShadowPath") -> tuple[Any, ...]:
    return tuple(attempt.proposal for attempt in path.attempts)


def assess_tdg6_continuous_admission(
    *,
    outer: "TDG6ShadowPath",
    medium: "TDG6ShadowPath",
    fine: "TDG6ShadowPath",
) -> TDG6ContinuousAdmissionEvidence:
    """Construct complete channelwise D01/D12 intervals from seven proposals."""

    if not all(isinstance(item, TDG6ShadowPath) for item in (outer, medium, fine)):
        raise TypeError("TDG6 assessment requires three shadow paths")
    if (outer.level, medium.level, fine.level) != ("outer", "medium", "fine"):
        raise ValueError("TDG6 paths are not in outer/medium/fine order")
    if not (outer.method == medium.method == fine.method):
        raise ValueError("TDG6 paths use different methods")
    if not (
        outer.initial_state_sha256
        == medium.initial_state_sha256
        == fine.initial_state_sha256
    ):
        raise ValueError("TDG6 paths do not share one bitwise initial state")
    if not (
        _same_binary64(outer.initial_time, medium.initial_time)
        and _same_binary64(outer.initial_time, fine.initial_time)
        and _same_binary64(outer.final_time, medium.final_time)
        and _same_binary64(outer.final_time, fine.final_time)
    ):
        raise ValueError("TDG6 paths do not span one macro interval")

    outer_proposals = _path_proposals(outer)
    medium_proposals = _path_proposals(medium)
    fine_proposals = _path_proposals(fine)
    outer_coefficients, outer_radii = _proposal_coefficients(outer_proposals[0])
    medium_coefficients = tuple(_proposal_coefficients(item) for item in medium_proposals)
    fine_coefficients = tuple(_proposal_coefficients(item) for item in fine_proposals)

    d01_values: list[np.ndarray] = []
    d01_radii: list[np.ndarray] = []
    for half, (medium_value, medium_radius) in enumerate(medium_coefficients):
        restricted, restricted_radius = tdg5._restrict_to_half_with_radius(
            outer_coefficients, outer_radii, half=half
        )
        difference, difference_radius = tdg5._difference_with_radius(
            restricted, restricted_radius, medium_value, medium_radius
        )
        d01_values.append(difference)
        d01_radii.append(difference_radius)

    d12_values: list[np.ndarray] = []
    d12_radii: list[np.ndarray] = []
    for medium_index, (medium_value, medium_radius) in enumerate(medium_coefficients):
        for local_half in (0, 1):
            fine_index = 2 * medium_index + local_half
            restricted, restricted_radius = tdg5._restrict_to_half_with_radius(
                medium_value, medium_radius, half=local_half
            )
            fine_value, fine_radius = fine_coefficients[fine_index]
            difference, difference_radius = tdg5._difference_with_radius(
                restricted, restricted_radius, fine_value, fine_radius
            )
            d12_values.append(difference)
            d12_radii.append(difference_radius)

    outer_stack = np.stack(d01_values, axis=0)
    outer_radius_stack = np.stack(d01_radii, axis=0)
    finest_stack = np.stack(d12_values, axis=0)
    finest_radius_stack = np.stack(d12_radii, axis=0)
    # (subinterval, block, owned_row, field, coefficient)
    if outer_stack.ndim != 5 or outer_stack.shape[1] != 3 or outer_stack.shape[3] != 6:
        raise ValueError("TDG6 complete-state coefficient shape differs")
    owned_rows = outer_stack.shape[2]
    outer_intervals: dict[str, TDG6Binary64MagnitudeInterval] = {}
    finest_intervals: dict[str, TDG6Binary64MagnitudeInterval] = {}
    for block_index, block in enumerate(("u", "p", "q")):
        for field_index, field in enumerate(_FIELD_NAMES):
            channel = f"{block}:{field}"
            outer_intervals[channel] = _magnitude_interval(
                outer_stack[:, block_index, :, field_index, :],
                outer_radius_stack[:, block_index, :, field_index, :],
                subinterval_count=2,
                owned_row_count=owned_rows,
            )
            finest_intervals[channel] = _magnitude_interval(
                finest_stack[:, block_index, :, field_index, :],
                finest_radius_stack[:, block_index, :, field_index, :],
                subinterval_count=4,
                owned_row_count=owned_rows,
            )
    return classify_tdg6_runtime_intervals(
        method=outer.method,
        outer_intervals=outer_intervals,
        finest_intervals=finest_intervals,
    )


@dataclass(frozen=True, slots=True)
class TDG6TracerSnapshot:
    labels_sha256: str
    positions_sha256: str
    proper_times_sha256: str
    event_history_sha256: str
    event_count: int
    cutoff: float
    outer_radius: float

    def __post_init__(self) -> None:
        for name in (
            "labels_sha256",
            "positions_sha256",
            "proper_times_sha256",
            "event_history_sha256",
        ):
            value = getattr(self, name)
            if not isinstance(value, str) or len(value) != 64:
                raise ValueError(f"{name} must be a SHA-256 digest")
        _nonnegative_integer("event_count", self.event_count)
        _finite("cutoff", self.cutoff, positive=True)
        _finite("outer_radius", self.outer_radius, positive=True)


def _tracer_snapshot(tracers: object) -> TDG6TracerSnapshot:
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
    if any(not hasattr(tracers, name) for name in required):
        raise TypeError("tracers do not expose the normal-flow runtime contract")
    labels = _immutable_array("tracer labels", getattr(tracers, "labels"))
    positions = _immutable_array("tracer positions", getattr(tracers, "positions"))
    proper = _immutable_array("tracer proper times", getattr(tracers, "proper_times"))
    if not (labels.shape == positions.shape == proper.shape):
        raise ValueError("tracer arrays have different shapes")
    event_times = tuple(
        _immutable_array("event proper times", item)
        for item in getattr(tracers, "event_proper_times")
    )
    event_fields = tuple(
        _immutable_array("event fields", item, ndim=2)
        for item in getattr(tracers, "event_fields")
    )
    if len(event_times) != len(event_fields):
        raise ValueError("tracer event histories have different lengths")
    history_arrays: tuple[np.ndarray, ...] = event_times + event_fields
    if not history_arrays:
        history_arrays = (np.zeros((0,), dtype=np.float64),)
    return TDG6TracerSnapshot(
        labels_sha256=array_content_sha256(labels),
        positions_sha256=array_content_sha256(positions),
        proper_times_sha256=array_content_sha256(proper),
        event_history_sha256=array_content_sha256(*history_arrays),
        event_count=len(event_times),
        cutoff=_finite("tracer cutoff", getattr(tracers, "cutoff"), positive=True),
        outer_radius=_finite(
            "tracer outer_radius", getattr(tracers, "outer_radius"), positive=True
        ),
    )


def _clone_tracers(tracers: object) -> object:
    _tracer_snapshot(tracers)
    clone = copy(tracers)
    clone.labels = np.asarray(tracers.labels, dtype=np.float64).copy()
    clone.positions = np.asarray(tracers.positions, dtype=np.float64).copy()
    clone.proper_times = np.asarray(tracers.proper_times, dtype=np.float64).copy()
    clone.event_proper_times = [
        np.asarray(item, dtype=np.float64).copy() for item in tracers.event_proper_times
    ]
    clone.event_fields = [
        np.asarray(item, dtype=np.float64).copy() for item in tracers.event_fields
    ]
    return clone


def _restore_tracers(target: object, source: object) -> None:
    """Restore the complete mutable normal-flow tracer boundary."""

    target.labels = np.asarray(source.labels, dtype=np.float64).copy()
    target.positions = np.asarray(source.positions, dtype=np.float64).copy()
    target.proper_times = np.asarray(source.proper_times, dtype=np.float64).copy()
    target.event_proper_times = [
        np.asarray(item, dtype=np.float64).copy()
        for item in source.event_proper_times
    ]
    target.event_fields = [
        np.asarray(item, dtype=np.float64).copy() for item in source.event_fields
    ]
    target.cutoff = source.cutoff
    target.outer_radius = source.outer_radius


@dataclass(frozen=True, slots=True)
class TDG6ShadowPath:
    level: str
    method: str
    initial_time: float
    final_time: float
    initial_state_sha256: str
    attempts: tuple[Proto7StepAttempt, ...]
    final_monitor_state: GR0RuntimeMonitorState
    final_causal_state: CausalBudgetState
    final_tracer_positions: np.ndarray
    final_tracer_proper_times: np.ndarray
    final_tracer_snapshot: TDG6TracerSnapshot

    def __post_init__(self) -> None:
        if self.level not in TDG6_LEVEL_STEP_COUNTS:
            raise ValueError("unknown TDG6 level")
        if self.method not in tdg5.TDG5_RUNTIME_METHODS:
            raise ValueError("unknown TDG6 path method")
        expected = TDG6_LEVEL_STEP_COUNTS[self.level]
        if len(self.attempts) != expected:
            raise ValueError("TDG6 shadow path has the wrong proposal count")
        if len(self.initial_state_sha256) != 64:
            raise ValueError("initial_state_sha256 must be a SHA-256 digest")
        if any(item.accepted is None or item.retry is not None for item in self.attempts):
            raise ValueError("TDG6 shadow path retained an unaccepted proposal")
        if not _same_binary64(self.attempts[0].proposal.initial_time, self.initial_time):
            raise ValueError("TDG6 shadow path starts at the wrong time")
        if not _same_binary64(self.attempts[-1].proposal.final_time, self.final_time):
            raise ValueError("TDG6 shadow path ends at the wrong time")
        for left, right in zip(self.attempts, self.attempts[1:]):
            if not _same_binary64(left.proposal.final_time, right.proposal.initial_time):
                raise ValueError("TDG6 shadow proposals are not contiguous")
            if tdg5._state_hash(left.proposal.candidate_state) != tdg5._state_hash(
                right.proposal.initial_state
            ):
                raise ValueError("TDG6 shadow proposal states are not contiguous")
        positions = _immutable_array("final_tracer_positions", self.final_tracer_positions)
        proper = _immutable_array("final_tracer_proper_times", self.final_tracer_proper_times)
        if positions.shape != proper.shape:
            raise ValueError("TDG6 final tracer arrays have different shapes")
        object.__setattr__(self, "final_tracer_positions", positions)
        object.__setattr__(self, "final_tracer_proper_times", proper)

    @property
    def final_accepted(self) -> Any:
        accepted = self.attempts[-1].accepted
        if accepted is None:
            raise RuntimeError("TDG6 shadow path has no accepted endpoint")
        return accepted

    @property
    def stage_record_count(self) -> int:
        return sum(len(item.proposal.stages) for item in self.attempts)


class TDG6RefinementPathStop(RuntimeError):
    """A source retry or inherited non-source stop on one shadow proposal."""

    def __init__(
        self,
        *,
        path: str,
        cause: BaseException | None = None,
        source_retry: ProposalFailureEvidence | None = None,
    ) -> None:
        if path not in TDG6_PATH_NAMES:
            raise ValueError("unknown TDG6 refinement path")
        if (cause is None) == (source_retry is None):
            raise ValueError("TDG6 path stop requires exactly one cause or source retry")
        self.path = path
        self.cause = cause
        self.source_retry = source_retry
        self.source_retry_owned_by_PROTO7 = source_retry is not None
        self.temporal_retry = False
        self.real_transaction_advanced = False
        self.accepted_state_advanced = False
        self.real_tracer_advanced = False
        kind = "source_retry" if source_retry is not None else type(cause).__name__
        super().__init__(f"TDG6 {path} shadow path stopped: {kind}")


def _shadow_path(
    *,
    level: str,
    method: str,
    boundaries: Sequence[float],
    state: EvolutionState,
    rhs: Callable[[float, EvolutionState], EvolutionRHS],
    projector: Callable[[float, EvolutionState], EvolutionState] | None,
    transaction: GR0RuntimeStageTransaction,
    tracers: object,
    coordinates: np.ndarray,
    previous_step_index: int,
    previous_transaction_serial: int,
) -> TDG6ShadowPath:
    if level not in TDG6_LEVEL_STEP_COUNTS:
        raise ValueError("unknown TDG6 shadow level")
    if len(boundaries) != TDG6_LEVEL_STEP_COUNTS[level] + 1:
        raise ValueError("TDG6 shadow boundaries differ from the level")
    shadow_transaction = tdg5._clone_gr0_transaction(transaction)
    shadow_tracers = _clone_tracers(tracers)
    current_state = state
    current_index = previous_step_index
    current_serial = previous_transaction_serial
    attempts: list[Proto7StepAttempt] = []
    for index, (left, right) in enumerate(zip(boundaries, boundaries[1:])):
        path_name = f"{level}_{index}"

        def preview(proposal: Any) -> Any:
            return shadow_tracers.preview_advance(
                old_state=current_state,
                new_state=proposal.candidate_state,
                coordinates=coordinates,
                step_size=proposal.final_time - left,
            )

        try:
            attempt = attempt_proto7_step(
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
            )
        except BaseException as error:
            if isinstance(error, (KeyboardInterrupt, SystemExit)):
                raise
            raise TDG6RefinementPathStop(path=path_name, cause=error) from error
        if attempt.retry is not None:
            raise TDG6RefinementPathStop(
                path=path_name, source_retry=attempt.retry
            )
        accepted = attempt.accepted
        if accepted is None or attempt.preaccept_payload is None:
            raise RuntimeError("TDG6 shadow attempt omitted acceptance or tracer preview")
        try:
            shadow_tracers.commit_advance(*attempt.preaccept_payload)
        except BaseException as error:
            if isinstance(error, (KeyboardInterrupt, SystemExit)):
                raise
            raise TDG6RefinementPathStop(path=path_name, cause=error) from error
        attempts.append(attempt)
        current_state = accepted.state
        current_index = accepted.step_index
        current_serial = accepted.transaction_serial
    return TDG6ShadowPath(
        level=level,
        method=method,
        initial_time=boundaries[0],
        final_time=boundaries[-1],
        initial_state_sha256=tdg5._state_hash(state),
        attempts=tuple(attempts),
        final_monitor_state=shadow_transaction.state,
        final_causal_state=shadow_transaction.causal_state,
        final_tracer_positions=np.asarray(shadow_tracers.positions),
        final_tracer_proper_times=np.asarray(shadow_tracers.proper_times),
        final_tracer_snapshot=_tracer_snapshot(shadow_tracers),
    )


@dataclass(frozen=True, slots=True)
class TDG6TemporalLedger:
    last_accepted_time: float
    accepted_macro_step_count: int
    cumulative_temporal_retry_count: int
    current_macro_step_temporal_retry_count: int
    last_accepted_macro_step_temporal_retry_count: int
    accumulated_debit_vector: tuple[float, ...]
    serialized_temporal_rejections: tuple[str, ...]

    def __post_init__(self) -> None:
        _finite("last_accepted_time", self.last_accepted_time)
        for name in (
            "accepted_macro_step_count",
            "cumulative_temporal_retry_count",
            "current_macro_step_temporal_retry_count",
            "last_accepted_macro_step_temporal_retry_count",
        ):
            _nonnegative_integer(name, getattr(self, name))
        if (
            self.current_macro_step_temporal_retry_count
            > self.cumulative_temporal_retry_count
            or self.last_accepted_macro_step_temporal_retry_count
            > self.cumulative_temporal_retry_count
        ):
            raise ValueError("TDG6 per-step retry count exceeds its cumulative ledger")
        if len(self.accumulated_debit_vector) != len(TDG6_COMPLETE_STATE_CHANNELS):
            raise ValueError("TDG6 ledger debit vector must have 18 channels")
        for value in self.accumulated_debit_vector:
            if _finite("accumulated debit", value) < 0.0:
                raise ValueError("TDG6 accumulated debit must be nonnegative")
        if len(self.serialized_temporal_rejections) != self.cumulative_temporal_retry_count:
            raise ValueError("TDG6 retry count differs from serialized evidence")
        for record in self.serialized_temporal_rejections:
            if not isinstance(record, str):
                raise TypeError("TDG6 serialized rejection must be text")
            decoded = json.loads(record)
            if _canonical_json(decoded) != record:
                raise ValueError("TDG6 rejection evidence is not canonical JSON")
            if decoded.get("event_type") != "rejected_TDG6_temporal_admission":
                raise ValueError("TDG6 rejection evidence has the wrong event type")

    @classmethod
    def zero(cls, *, initial_time: Real) -> "TDG6TemporalLedger":
        return cls(
            last_accepted_time=_finite("initial_time", initial_time),
            accepted_macro_step_count=0,
            cumulative_temporal_retry_count=0,
            current_macro_step_temporal_retry_count=0,
            last_accepted_macro_step_temporal_retry_count=0,
            accumulated_debit_vector=(0.0,) * len(TDG6_COMPLETE_STATE_CHANNELS),
            serialized_temporal_rejections=(),
        )


@dataclass(frozen=True, slots=True)
class TDG6PreparedGR0Compositor:
    method: str
    initial_time: float
    final_time: float
    initial_state_sha256: str
    previous_step_index: int
    previous_transaction_serial: int
    original_monitor_state: GR0RuntimeMonitorState
    original_causal_state: CausalBudgetState
    original_tracer_snapshot: TDG6TracerSnapshot
    original_temporal_ledger: TDG6TemporalLedger
    outer: TDG6ShadowPath
    medium: TDG6ShadowPath
    fine: TDG6ShadowPath
    continuous_admission: TDG6ContinuousAdmissionEvidence
    outer_and_medium_are_evidence_only: bool = True
    only_complete_fine_path_is_committable: bool = True
    real_transaction_advanced_during_prepare: bool = False
    accepted_state_advanced_during_prepare: bool = False
    real_tracer_advanced_during_prepare: bool = False
    real_temporal_ledger_advanced_during_prepare: bool = False

    def __post_init__(self) -> None:
        if self.method not in tdg5.TDG5_RUNTIME_METHODS:
            raise ValueError("unknown TDG6 prepared method")
        if len(self.initial_state_sha256) != 64:
            raise ValueError("initial_state_sha256 must be a SHA-256 digest")
        _nonnegative_integer("previous_step_index", self.previous_step_index)
        _nonnegative_integer(
            "previous_transaction_serial", self.previous_transaction_serial
        )
        if (self.outer.level, self.medium.level, self.fine.level) != (
            "outer",
            "medium",
            "fine",
        ):
            raise ValueError("TDG6 prepared paths are reordered")
        if not all(
            path.method == self.method for path in (self.outer, self.medium, self.fine)
        ):
            raise ValueError("TDG6 prepared paths use different methods")
        if not _same_binary64(self.outer.initial_time, self.initial_time):
            raise ValueError("TDG6 prepared compositor starts at the wrong time")
        if not _same_binary64(self.fine.final_time, self.final_time):
            raise ValueError("TDG6 prepared compositor ends at the wrong time")
        if self.continuous_admission.method != self.method:
            raise ValueError("TDG6 prepared admission uses a different method")
        if (
            self.outer_and_medium_are_evidence_only is not True
            or self.only_complete_fine_path_is_committable is not True
            or self.real_transaction_advanced_during_prepare is not False
            or self.accepted_state_advanced_during_prepare is not False
            or self.real_tracer_advanced_during_prepare is not False
            or self.real_temporal_ledger_advanced_during_prepare is not False
        ):
            raise ValueError("TDG6 prepared compositor violates atomic ownership")


def _revalidate_prepared_boundary(
    prepared: TDG6PreparedGR0Compositor,
    *,
    transaction: GR0RuntimeStageTransaction,
    tracers: object,
    temporal_ledger: TDG6TemporalLedger,
    current_time: Real,
    current_state: EvolutionState,
    current_step_index: int,
    current_transaction_serial: int,
) -> None:
    """Prove that no accepted external boundary moved after preparation."""

    if not isinstance(prepared, TDG6PreparedGR0Compositor):
        raise TypeError("prepared must be TDG6PreparedGR0Compositor")
    if not isinstance(transaction, GR0RuntimeStageTransaction):
        raise TypeError("transaction must be GR0RuntimeStageTransaction")
    if not isinstance(temporal_ledger, TDG6TemporalLedger):
        raise TypeError("temporal_ledger must be TDG6TemporalLedger")
    if not isinstance(current_state, EvolutionState):
        raise TypeError("current_state must be EvolutionState")
    time = _finite("current_time", current_time)
    step_index = _nonnegative_integer("current_step_index", current_step_index)
    serial = _nonnegative_integer(
        "current_transaction_serial", current_transaction_serial
    )
    if (
        not _same_binary64(time, prepared.initial_time)
        or tdg5._state_hash(current_state) != prepared.initial_state_sha256
        or step_index != prepared.previous_step_index
        or serial != prepared.previous_transaction_serial
        or transaction.state != prepared.original_monitor_state
        or transaction.causal_state != prepared.original_causal_state
        or _tracer_snapshot(tracers) != prepared.original_tracer_snapshot
        or temporal_ledger != prepared.original_temporal_ledger
    ):
        raise RuntimeError("TDG6 external boundary drifted after shadow preparation")


def _subdivision_boundaries(start: float, width: float, count: int) -> tuple[float, ...]:
    step = width / count
    boundaries = tuple(start + index * step for index in range(count + 1))
    expected_final = start + width
    if not _same_binary64(boundaries[-1], expected_final):
        raise ValueError("TDG6 subdivision does not share the macro endpoint")
    for left, right in zip(boundaries, boundaries[1:]):
        if not _same_binary64(right - left, step):
            raise ValueError("TDG6 subdivision is not bitwise uniform")
    return boundaries


def prepare_tdg6_gr0_compositor(
    *,
    method: str,
    time: Real,
    step_size: Real,
    state: EvolutionState,
    rhs: Callable[[float, EvolutionState], EvolutionRHS],
    projector: Callable[[float, EvolutionState], EvolutionState] | None,
    transaction: GR0RuntimeStageTransaction,
    tracers: object,
    coordinates: object,
    temporal_ledger: TDG6TemporalLedger,
    previous_step_index: int,
    previous_transaction_serial: int,
) -> TDG6PreparedGR0Compositor:
    """Prepare all seven shadows and complete admission without real mutation."""

    if method not in tdg5.TDG5_RUNTIME_METHODS:
        raise ValueError("unknown TDG6 runtime method")
    if not isinstance(state, EvolutionState) or state.shape[1] != 6:
        raise ValueError("TDG6 state must contain the six frozen fields")
    if not callable(rhs):
        raise TypeError("rhs must be callable")
    if projector is not None and not callable(projector):
        raise TypeError("projector must be callable or None")
    if not isinstance(transaction, GR0RuntimeStageTransaction):
        raise TypeError("transaction must be GR0RuntimeStageTransaction")
    if not isinstance(temporal_ledger, TDG6TemporalLedger):
        raise TypeError("temporal_ledger must be TDG6TemporalLedger")
    start = _finite("time", time)
    width = _finite("step_size", step_size, positive=True)
    step_index = _nonnegative_integer("previous_step_index", previous_step_index)
    serial = _nonnegative_integer(
        "previous_transaction_serial", previous_transaction_serial
    )
    coordinate_array = _immutable_array("coordinates", coordinates)
    if coordinate_array.size != state.shape[0]:
        raise ValueError("TDG6 coordinate array differs from state grid")
    if serial != transaction.state.last_transaction_serial + 1:
        raise ValueError("TDG6 member serial does not follow the transaction ledger")
    if not _same_binary64(transaction.causal_state.accepted_time, start):
        raise ValueError("TDG6 causal ledger time differs from accepted time")
    if not _same_binary64(temporal_ledger.last_accepted_time, start):
        raise ValueError("TDG6 temporal ledger time differs from accepted time")

    original_monitor = transaction.state
    original_causal = transaction.causal_state
    original_tracer = _tracer_snapshot(tracers)
    original_ledger = temporal_ledger
    initial_hash = tdg5._state_hash(state)
    outer_boundaries = _subdivision_boundaries(start, width, 1)
    medium_boundaries = _subdivision_boundaries(start, width, 2)
    fine_boundaries = _subdivision_boundaries(start, width, 4)

    outer = _shadow_path(
        level="outer",
        method=method,
        boundaries=outer_boundaries,
        state=state,
        rhs=rhs,
        projector=projector,
        transaction=transaction,
        tracers=tracers,
        coordinates=coordinate_array,
        previous_step_index=step_index,
        previous_transaction_serial=serial,
    )
    medium = _shadow_path(
        level="medium",
        method=method,
        boundaries=medium_boundaries,
        state=state,
        rhs=rhs,
        projector=projector,
        transaction=transaction,
        tracers=tracers,
        coordinates=coordinate_array,
        previous_step_index=step_index,
        previous_transaction_serial=serial,
    )
    fine = _shadow_path(
        level="fine",
        method=method,
        boundaries=fine_boundaries,
        state=state,
        rhs=rhs,
        projector=projector,
        transaction=transaction,
        tracers=tracers,
        coordinates=coordinate_array,
        previous_step_index=step_index,
        previous_transaction_serial=serial,
    )
    if transaction.state != original_monitor or transaction.causal_state != original_causal:
        raise RuntimeError("TDG6 shadow preparation mutated the real transaction")
    if _tracer_snapshot(tracers) != original_tracer:
        raise RuntimeError("TDG6 shadow preparation mutated the real tracers")
    if tdg5._state_hash(state) != initial_hash:
        raise RuntimeError("TDG6 shadow preparation mutated the accepted state")
    if temporal_ledger != original_ledger:
        raise RuntimeError("TDG6 shadow preparation mutated the temporal ledger")

    admission = assess_tdg6_continuous_admission(
        outer=outer, medium=medium, fine=fine
    )
    return TDG6PreparedGR0Compositor(
        method=method,
        initial_time=start,
        final_time=start + width,
        initial_state_sha256=initial_hash,
        previous_step_index=step_index,
        previous_transaction_serial=serial,
        original_monitor_state=original_monitor,
        original_causal_state=original_causal,
        original_tracer_snapshot=original_tracer,
        original_temporal_ledger=original_ledger,
        outer=outer,
        medium=medium,
        fine=fine,
        continuous_admission=admission,
    )


@dataclass(frozen=True, slots=True)
class TDG6TemporalRetryEvidence:
    event_type: str
    method: str
    initial_time: float
    attempted_macro_step_size: float
    retry_step_size: float
    initial_state_sha256: str
    previous_step_index: int
    previous_transaction_serial: int
    retry_count_for_current_macro_step: int
    cumulative_temporal_retry_count: int
    failed_channels: tuple[str, ...]
    channel_admissions: tuple[TDG6RuntimeChannelAdmission, ...]
    accepted_state_bitwise_preserved: bool = True
    real_monitor_and_causal_ledgers_preserved: bool = True
    real_tracer_ledger_preserved: bool = True
    accumulated_debit_unchanged: bool = True
    classification_is_numerical_not_physical: bool = True

    def __post_init__(self) -> None:
        if self.event_type != "rejected_TDG6_temporal_admission":
            raise ValueError("TDG6 retry evidence has the wrong event type")
        if self.method not in tdg5.TDG5_RUNTIME_METHODS:
            raise ValueError("TDG6 retry evidence has an unknown method")
        start = _finite("initial_time", self.initial_time)
        attempted = _finite(
            "attempted_macro_step_size", self.attempted_macro_step_size, positive=True
        )
        retry = _finite("retry_step_size", self.retry_step_size, positive=True)
        if not _same_binary64(retry, attempted * float(TDG6_RETRY_FACTOR)):
            raise ValueError("TDG6 retry size differs from the frozen factor")
        if len(self.initial_state_sha256) != 64:
            raise ValueError("initial_state_sha256 must be a SHA-256 digest")
        _nonnegative_integer("previous_step_index", self.previous_step_index)
        _nonnegative_integer(
            "previous_transaction_serial", self.previous_transaction_serial
        )
        if self.retry_count_for_current_macro_step < 1:
            raise ValueError("TDG6 retry evidence must increment the step counter")
        if self.cumulative_temporal_retry_count < self.retry_count_for_current_macro_step:
            raise ValueError("TDG6 cumulative retry count is inconsistent")
        observed_failed = tuple(
            item.channel for item in self.channel_admissions if not item.admission_passed
        )
        if not observed_failed or observed_failed != self.failed_channels:
            raise ValueError("TDG6 retry evidence omits failed channels")
        if tuple(item.channel for item in self.channel_admissions) != TDG6_COMPLETE_STATE_CHANNELS:
            raise ValueError("TDG6 retry evidence omits channel decisions")
        if (
            self.accepted_state_bitwise_preserved is not True
            or self.real_monitor_and_causal_ledgers_preserved is not True
            or self.real_tracer_ledger_preserved is not True
            or self.accumulated_debit_unchanged is not True
            or self.classification_is_numerical_not_physical is not True
        ):
            raise ValueError("TDG6 retry evidence crossed its no-commit boundary")
        _ = start

    def as_mapping(self) -> dict[str, object]:
        return asdict(self)


class TDG6TemporalRetryRequired(RuntimeError):
    """An admitted numerical recovery: halve the macro step and retry."""

    def __init__(
        self,
        evidence: TDG6TemporalRetryEvidence,
        updated_ledger: TDG6TemporalLedger,
    ) -> None:
        self.evidence = evidence
        self.updated_ledger = updated_ledger
        self.retry_step_size = evidence.retry_step_size
        self.physical_classification = False
        super().__init__("TDG6 temporal admission requested a half-step retry")


class TDG6TemporalRetryExhausted(RuntimeError):
    """Frozen numerical retry count or minimum macro step was exhausted."""

    def __init__(
        self,
        *,
        reason: str,
        evidence: TDG6TemporalRetryEvidence,
        updated_ledger: TDG6TemporalLedger,
    ) -> None:
        if reason not in {"maximum_temporal_retries", "minimum_macro_step"}:
            raise ValueError("unknown TDG6 retry exhaustion reason")
        self.reason = reason
        self.evidence = evidence
        self.updated_ledger = updated_ledger
        self.physical_classification = False
        super().__init__(f"TDG6 numerical temporal retry exhausted: {reason}")


def require_tdg6_temporal_admission(
    prepared: TDG6PreparedGR0Compositor,
    *,
    transaction: GR0RuntimeStageTransaction,
    tracers: object,
    temporal_ledger: TDG6TemporalLedger,
    current_time: Real,
    current_state: EvolutionState,
    current_step_index: int,
    current_transaction_serial: int,
    durable_rejection_sink: Callable[[Mapping[str, object]], None],
) -> None:
    """Pass silently or durably record and raise one typed numerical retry."""

    if not callable(durable_rejection_sink):
        raise TypeError("durable_rejection_sink must be callable")
    _revalidate_prepared_boundary(
        prepared,
        transaction=transaction,
        tracers=tracers,
        temporal_ledger=temporal_ledger,
        current_time=current_time,
        current_state=current_state,
        current_step_index=current_step_index,
        current_transaction_serial=current_transaction_serial,
    )
    admission = prepared.continuous_admission
    if admission.complete_admission_passed:
        return

    next_step = (prepared.final_time - prepared.initial_time) * float(TDG6_RETRY_FACTOR)
    current_count = temporal_ledger.current_macro_step_temporal_retry_count + 1
    cumulative = temporal_ledger.cumulative_temporal_retry_count + 1
    evidence = TDG6TemporalRetryEvidence(
        event_type="rejected_TDG6_temporal_admission",
        method=prepared.method,
        initial_time=prepared.initial_time,
        attempted_macro_step_size=prepared.final_time - prepared.initial_time,
        retry_step_size=next_step,
        initial_state_sha256=prepared.initial_state_sha256,
        previous_step_index=prepared.previous_step_index,
        previous_transaction_serial=prepared.previous_transaction_serial,
        retry_count_for_current_macro_step=current_count,
        cumulative_temporal_retry_count=cumulative,
        failed_channels=admission.failed_channels,
        channel_admissions=admission.channel_admissions,
    )
    mapping = evidence.as_mapping()
    # The sink is intentionally called before any next proposal can be
    # represented.  An I/O failure leaves the caller's immutable ledger intact.
    durable_rejection_sink(mapping)
    serialized = _canonical_json(mapping)
    updated = TDG6TemporalLedger(
        last_accepted_time=temporal_ledger.last_accepted_time,
        accepted_macro_step_count=temporal_ledger.accepted_macro_step_count,
        cumulative_temporal_retry_count=cumulative,
        current_macro_step_temporal_retry_count=current_count,
        last_accepted_macro_step_temporal_retry_count=(
            temporal_ledger.last_accepted_macro_step_temporal_retry_count
        ),
        accumulated_debit_vector=temporal_ledger.accumulated_debit_vector,
        serialized_temporal_rejections=(
            *temporal_ledger.serialized_temporal_rejections,
            serialized,
        ),
    )
    if current_count > TDG6_MAXIMUM_TEMPORAL_RETRIES_PER_STEP:
        raise TDG6TemporalRetryExhausted(
            reason="maximum_temporal_retries",
            evidence=evidence,
            updated_ledger=updated,
        )
    if next_step < float(TDG6_MINIMUM_MACRO_STEP):
        raise TDG6TemporalRetryExhausted(
            reason="minimum_macro_step", evidence=evidence, updated_ledger=updated
        )
    raise TDG6TemporalRetryRequired(evidence, updated)


@dataclass(frozen=True, slots=True)
class TDG6CommittedGR0Compositor:
    method: str
    time: float
    step_index: int
    transaction_serial: int
    state: EvolutionState
    temporal_ledger: TDG6TemporalLedger
    continuous_admission: TDG6ContinuousAdmissionEvidence
    committed_path: str = "fine_four_quarter_steps"
    outer_path_committed: bool = False
    medium_path_committed: bool = False
    fine_quarter_steps_committed_atomically: bool = True

    def __post_init__(self) -> None:
        if self.method not in tdg5.TDG5_RUNTIME_METHODS:
            raise ValueError("unknown TDG6 committed method")
        _nonnegative_integer("step_index", self.step_index)
        _nonnegative_integer("transaction_serial", self.transaction_serial)
        if not isinstance(self.state, EvolutionState):
            raise TypeError("TDG6 committed state must be EvolutionState")
        if not self.continuous_admission.complete_admission_passed:
            raise ValueError("TDG6 committed a failed temporal admission")
        if (
            self.committed_path != "fine_four_quarter_steps"
            or self.outer_path_committed is not False
            or self.medium_path_committed is not False
            or self.fine_quarter_steps_committed_atomically is not True
        ):
            raise ValueError("TDG6 commit exposed the wrong path")


def commit_tdg6_gr0_compositor(
    prepared: TDG6PreparedGR0Compositor,
    *,
    transaction: GR0RuntimeStageTransaction,
    tracers: object,
    temporal_ledger: TDG6TemporalLedger,
    current_time: Real,
    current_state: EvolutionState,
    current_step_index: int,
    current_transaction_serial: int,
) -> TDG6CommittedGR0Compositor:
    """Atomically adopt only admitted fine state/monitor/causal/tracer/debit."""

    if not prepared.continuous_admission.complete_admission_passed:
        raise RuntimeError("TDG6 failed admission cannot commit")
    _revalidate_prepared_boundary(
        prepared,
        transaction=transaction,
        tracers=tracers,
        temporal_ledger=temporal_ledger,
        current_time=current_time,
        current_state=current_state,
        current_step_index=current_step_index,
        current_transaction_serial=current_transaction_serial,
    )

    debit = tuple(
        _outward_nonnegative_sum(old, new)
        for old, new in zip(
            temporal_ledger.accumulated_debit_vector,
            prepared.continuous_admission.finest_pair_debit_vector,
            strict=True,
        )
    )
    new_ledger = TDG6TemporalLedger(
        last_accepted_time=prepared.final_time,
        accepted_macro_step_count=temporal_ledger.accepted_macro_step_count + 1,
        cumulative_temporal_retry_count=temporal_ledger.cumulative_temporal_retry_count,
        current_macro_step_temporal_retry_count=0,
        last_accepted_macro_step_temporal_retry_count=(
            temporal_ledger.current_macro_step_temporal_retry_count
        ),
        accumulated_debit_vector=debit,
        serialized_temporal_rejections=temporal_ledger.serialized_temporal_rejections,
    )
    fine = prepared.fine.final_accepted
    committed = TDG6CommittedGR0Compositor(
        method=prepared.method,
        time=fine.time,
        step_index=fine.step_index,
        transaction_serial=fine.transaction_serial,
        state=fine.state,
        temporal_ledger=new_ledger,
        continuous_admission=prepared.continuous_admission,
    )
    prior_monitor = transaction.state
    prior_causal = transaction.causal_state
    prior_tracers = _clone_tracers(tracers)
    try:
        transaction.state = prepared.fine.final_monitor_state
        transaction.causal_state = prepared.fine.final_causal_state
        tracers.commit_advance(
            prepared.fine.final_tracer_positions,
            prepared.fine.final_tracer_proper_times,
        )
        if _tracer_snapshot(tracers) != prepared.fine.final_tracer_snapshot:
            raise RuntimeError("TDG6 fine tracer commit differs from its shadow")
    except BaseException:
        transaction.state = prior_monitor
        transaction.causal_state = prior_causal
        try:
            _restore_tracers(tracers, prior_tracers)
        except BaseException as rollback_error:
            raise RuntimeError("TDG6 tracer rollback failed") from rollback_error
        raise
    return committed


def tdg6_checkpoint_extension(
    temporal_ledger: TDG6TemporalLedger,
) -> tuple[dict[str, object], np.ndarray]:
    """Encode every TDG6 ledger field for an enclosing atomic checkpoint."""

    if not isinstance(temporal_ledger, TDG6TemporalLedger):
        raise TypeError("temporal_ledger must be TDG6TemporalLedger")
    debit = np.asarray(temporal_ledger.accumulated_debit_vector, dtype=np.float64)
    debit_hash = array_content_sha256(debit)
    metadata: dict[str, object] = {
        "schema_version": TDG6_RUNTIME_SCHEMA_VERSION,
        "channel_order": list(TDG6_COMPLETE_STATE_CHANNELS),
        "last_accepted_time": temporal_ledger.last_accepted_time,
        "accepted_macro_step_count": temporal_ledger.accepted_macro_step_count,
        "cumulative_temporal_retry_count": (
            temporal_ledger.cumulative_temporal_retry_count
        ),
        "current_macro_step_temporal_retry_count": (
            temporal_ledger.current_macro_step_temporal_retry_count
        ),
        "last_accepted_macro_step_temporal_retry_count": (
            temporal_ledger.last_accepted_macro_step_temporal_retry_count
        ),
        "accumulated_debit_sha256": debit_hash,
        "serialized_temporal_rejections": list(
            temporal_ledger.serialized_temporal_rejections
        ),
        "complete_rejection_evidence_present": True,
        "source_and_CFL_retry_counters_owned_elsewhere": True,
        "physical_classification": False,
    }
    return metadata, debit.copy()


def restore_tdg6_checkpoint_extension(
    metadata: Mapping[str, object],
    accumulated_debit: object,
) -> TDG6TemporalLedger:
    """Restore and validate the exact TDG6 checkpoint extension."""

    required = {
        "schema_version",
        "channel_order",
        "last_accepted_time",
        "accepted_macro_step_count",
        "cumulative_temporal_retry_count",
        "current_macro_step_temporal_retry_count",
        "last_accepted_macro_step_temporal_retry_count",
        "accumulated_debit_sha256",
        "serialized_temporal_rejections",
        "complete_rejection_evidence_present",
        "source_and_CFL_retry_counters_owned_elsewhere",
        "physical_classification",
    }
    if not isinstance(metadata, Mapping) or set(metadata) != required:
        raise ValueError("TDG6 checkpoint metadata is incomplete or broadened")
    if metadata["schema_version"] != TDG6_RUNTIME_SCHEMA_VERSION:
        raise ValueError("TDG6 checkpoint schema differs")
    if tuple(metadata["channel_order"]) != TDG6_COMPLETE_STATE_CHANNELS:
        raise ValueError("TDG6 checkpoint channel order differs")
    debit = _immutable_array("accumulated_debit", accumulated_debit)
    if debit.shape != (len(TDG6_COMPLETE_STATE_CHANNELS),):
        raise ValueError("TDG6 checkpoint debit shape differs")
    if array_content_sha256(debit) != metadata["accumulated_debit_sha256"]:
        raise ValueError("TDG6 checkpoint debit hash differs")
    if (
        metadata["complete_rejection_evidence_present"] is not True
        or metadata["source_and_CFL_retry_counters_owned_elsewhere"] is not True
        or metadata["physical_classification"] is not False
    ):
        raise ValueError("TDG6 checkpoint crossed its numerical scope")
    serialized = metadata["serialized_temporal_rejections"]
    if not isinstance(serialized, list) or any(not isinstance(item, str) for item in serialized):
        raise ValueError("TDG6 checkpoint rejection ledger is malformed")
    return TDG6TemporalLedger(
        last_accepted_time=metadata["last_accepted_time"],
        accepted_macro_step_count=metadata["accepted_macro_step_count"],
        cumulative_temporal_retry_count=metadata["cumulative_temporal_retry_count"],
        current_macro_step_temporal_retry_count=(
            metadata["current_macro_step_temporal_retry_count"]
        ),
        last_accepted_macro_step_temporal_retry_count=(
            metadata["last_accepted_macro_step_temporal_retry_count"]
        ),
        accumulated_debit_vector=tuple(float(item) for item in debit),
        serialized_temporal_rejections=tuple(serialized),
    )


__all__ = [
    "TDG6Binary64MagnitudeInterval",
    "TDG6CommittedGR0Compositor",
    "TDG6ContinuousAdmissionEvidence",
    "TDG6PreparedGR0Compositor",
    "TDG6RefinementPathStop",
    "TDG6RuntimeChannelAdmission",
    "TDG6ShadowPath",
    "TDG6TemporalLedger",
    "TDG6TemporalRetryEvidence",
    "TDG6TemporalRetryExhausted",
    "TDG6TemporalRetryRequired",
    "assess_tdg6_continuous_admission",
    "classify_tdg6_runtime_intervals",
    "commit_tdg6_gr0_compositor",
    "prepare_tdg6_gr0_compositor",
    "require_tdg6_temporal_admission",
    "restore_tdg6_checkpoint_extension",
    "tdg6_checkpoint_extension",
]
