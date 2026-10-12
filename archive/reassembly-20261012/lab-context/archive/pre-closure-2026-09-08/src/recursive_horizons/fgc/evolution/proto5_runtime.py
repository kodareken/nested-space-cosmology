"""Fail-closed PROTO5 runtime composition for the GR-0 calibration branch.

HLT2 partitions the historical HLT1 decisions into nine universal runtime
stops, eleven candidate-action-only stops, and six numerical decisions whose
successor is evaluated at common events.  GR-0 has no FGC-QR nonlinear branch,
effective-Planck deformation, or candidate characteristic frame.  Filling
those fields with convenient constants would therefore be scientifically
wrong.  This module instead implements the exact GR-0 partition:

* :class:`GR0RuntimeStageTransaction` evaluates the nine universal stops at
  every proposed Runge--Kutta stage and the candidate endpoint, composes them
  atomically with BND2's causal ledger, and treats a CFL miss as a retryable
  numerical proposal rather than a physical obstruction;
* :func:`proto5_gr0_common_event` evaluates the three-grid PROTO5 constraint
  rule and the inherited HLT2 spatial weighted spectrum at one bitwise-aligned
  coordinate time; and
* :func:`proto5_temporal_spectral_admission` evaluates a causal-past proper-
  time spectrum after at least 64 samples, including the interpolation error
  incurred when a nonuniform proper-time history is put on a uniform FFT grid.

Nothing here creates a run directory, advances a calibration by itself, or
classifies a trapped sphere.  The authorized runner must supply actual stage
records and common-event states.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite
from numbers import Real
from typing import Mapping, Sequence

import numpy as np

from .boundary_domain import (
    BoundaryGeometry,
    CausalBudgetState,
    assess_causal_step,
)
from .calibration_runtime import (
    GR0_UNIVERSAL_RUNTIME_STOP_IDS,
    SemidiscreteConstraintAdmission,
    gr0_semidiscrete_constraint_snapshot,
    proto5_semidiscrete_constraint_admission,
)
from .gr0_calibration import radial_null_observables
from .numerical_engine import EvolutionState, StageRecord, StepProposal, UniformRadialGrid
from .proto4_admission import (
    PROTO4_SPECTRAL_FIELD_ORDER,
    NestedSpectralAdmission,
    SpectralPowerBudget,
    SpectralThresholds,
    nested_spectral_admission,
    proper_radial_profile,
    spectral_field_budgets,
    windowed_spectral_power_budget,
)
from .health_monitor import compact_vacuum_buffer_window


GR0_RUNTIME_STOP_PRIORITY = tuple(GR0_UNIVERSAL_RUNTIME_STOP_IDS)
PROTO5_METHOD_CONTRACT: Mapping[str, Mapping[str, float | int]] = {
    "RK4": {
        "spatial_order": 4,
        "coarsest_constraint_guard": 1.0 / 50.0,
        "finest_constraint_guard": 1.0 / 1000.0,
    },
    "SSPRK3": {
        "spatial_order": 2,
        "coarsest_constraint_guard": 1.0 / 10.0,
        "finest_constraint_guard": 1.0 / 200.0,
    },
}


@dataclass(frozen=True)
class EventLogRecoveryPlan:
    action: str
    restored_payload: bytes
    orphaned_tail: bytes | None


def proto5_event_log_recovery_plan(
    current_payload: bytes | None,
    checkpoint_payload: bytes,
) -> EventLogRecoveryPlan:
    """Return the only admissible recovery from a two-file commit interruption.

    The atomic checkpoint owns accepted member state and embeds the exact log
    bytes that correspond to it.  A live log may be absent, identical, behind
    that checkpoint, or ahead by complete uncommitted JSONL records.  Only the
    last case yields an orphan tail; a divergent history is never merged.
    """

    if not isinstance(checkpoint_payload, bytes) or (
        checkpoint_payload and not checkpoint_payload.endswith(b"\n")
    ):
        raise ValueError("checkpoint event log must end at a complete line")
    if current_payload is None:
        return EventLogRecoveryPlan(
            "restored_missing_log_from_checkpoint",
            checkpoint_payload,
            None,
        )
    if not isinstance(current_payload, bytes) or (
        current_payload and not current_payload.endswith(b"\n")
    ):
        raise ValueError("current event log must end at a complete line")
    if current_payload == checkpoint_payload:
        return EventLogRecoveryPlan(
            "event_log_already_matches_checkpoint",
            checkpoint_payload,
            None,
        )
    if current_payload.startswith(checkpoint_payload):
        tail = current_payload[len(checkpoint_payload) :]
        if tail and not tail.endswith(b"\n"):
            raise ValueError("uncommitted event-log tail is incomplete")
        return EventLogRecoveryPlan(
            "archived_uncommitted_tail_and_restored_checkpoint",
            checkpoint_payload,
            tail,
        )
    if checkpoint_payload.startswith(current_payload):
        return EventLogRecoveryPlan(
            "completed_checkpoint_owned_log_write",
            checkpoint_payload,
            None,
        )
    raise ValueError("event log diverges from the atomic checkpoint history")


def _finite(name: str, value: Real, *, positive: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a real scalar")
    answer = float(value)
    if not isfinite(answer):
        raise ValueError(f"{name} must be finite")
    if positive and answer <= 0.0:
        raise ValueError(f"{name} must be positive")
    return answer


def _integer(name: str, value: object, *, minimum: int = 0) -> int:
    if isinstance(value, bool) or not isinstance(value, (int, np.integer)):
        raise TypeError(f"{name} must be an integer")
    answer = int(value)
    if answer < minimum:
        raise ValueError(f"{name} must be at least {minimum}")
    return answer


def _diagnostic(record: StageRecord, key: str) -> object:
    if key not in record.rhs.diagnostics:
        raise ValueError(f"GR-0 stage diagnostics omit {key}")
    return record.rhs.diagnostics[key]


def gr0_hat_lorentzian_margin(
    state: EvolutionState,
    *,
    hat_normal_factor: Real = 9.0,
) -> float:
    """Return an explicit positive margin for the spherical hat inverse.

    For the ADM base metric and ``q>1``, the auxiliary inverse has
    ``-hat_g^tt=q/alpha^2`` and negative radial-base determinant
    ``-det(hat_g_base)=q/(alpha^2 lambda^2)``.  Positivity of these quantities
    together with ``1/lambda^2`` is equivalent to a Lorentzian radial base
    with the coordinate-time slices spacelike.  Angular positivity follows
    from ``R>0`` away from the separately handled regular centre.
    """

    if not isinstance(state, EvolutionState) or state.shape[1] != 6:
        raise TypeError("state must contain the six canonical ADM fields")
    factor = _finite("hat_normal_factor", hat_normal_factor, positive=True)
    if factor <= 1.0:
        raise ValueError("hat_normal_factor must exceed one")
    alpha = state.u[:, 0]
    radial_metric = state.u[:, 2]
    if np.any(~np.isfinite(alpha)) or np.any(~np.isfinite(radial_metric)):
        raise ValueError("hat-cone metric factors became nonfinite")
    with np.errstate(divide="ignore", invalid="ignore"):
        time_margin = factor / alpha**2
        determinant_margin = factor / (alpha**2 * radial_metric**2)
        radial_spatial_margin = 1.0 / radial_metric**2
    margins = np.minimum.reduce(
        (time_margin, determinant_margin, radial_spatial_margin)
    )
    if np.any(~np.isfinite(margins)):
        return 0.0
    return float(np.min(margins, initial=np.inf))


@dataclass(frozen=True, slots=True)
class GR0UniversalThresholds:
    source_residual_maximum: float = 1.0e-12
    source_iteration_maximum: int = 16
    kinetic_condition_maximum: float = 1.0e10

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "source_residual_maximum",
            _finite(
                "source_residual_maximum",
                self.source_residual_maximum,
                positive=True,
            ),
        )
        object.__setattr__(
            self,
            "source_iteration_maximum",
            _integer(
                "source_iteration_maximum",
                self.source_iteration_maximum,
                minimum=0,
            ),
        )
        object.__setattr__(
            self,
            "kinetic_condition_maximum",
            _finite(
                "kinetic_condition_maximum",
                self.kinetic_condition_maximum,
                positive=True,
            ),
        )


@dataclass(frozen=True, slots=True)
class GR0UniversalStageSnapshot:
    time: float
    transaction_serial: int
    minimum_lapse: float
    minimum_radial_metric: float
    minimum_areal_radius_away_from_center: float
    minimum_hat_lorentzian_margin: float
    source_residual_infinity: float
    source_iterations: int
    source_residual_decreased_monotonically: bool
    kinetic_condition_infinity: float
    boundary_margin_over_required_buffer: float

    def __post_init__(self) -> None:
        object.__setattr__(self, "time", _finite("time", self.time))
        object.__setattr__(
            self,
            "transaction_serial",
            _integer("transaction_serial", self.transaction_serial),
        )
        object.__setattr__(
            self,
            "source_iterations",
            _integer("source_iterations", self.source_iterations),
        )
        if not isinstance(self.source_residual_decreased_monotonically, bool):
            raise TypeError("source_residual_decreased_monotonically must be bool")
        for name in (
            "minimum_lapse",
            "minimum_radial_metric",
            "minimum_areal_radius_away_from_center",
            "minimum_hat_lorentzian_margin",
            "source_residual_infinity",
            "kinetic_condition_infinity",
            "boundary_margin_over_required_buffer",
        ):
            object.__setattr__(self, name, _finite(name, getattr(self, name)))


def gr0_universal_stage_snapshot(
    record: StageRecord,
    *,
    transaction_serial: int,
    boundary_margin_over_required_buffer: Real,
    hat_normal_factor: Real = 9.0,
) -> GR0UniversalStageSnapshot:
    """Translate only actual GR-0 stage diagnostics into universal fields."""

    if not isinstance(record, StageRecord):
        raise TypeError("record must be a StageRecord")
    iterations = _diagnostic(record, "source_refinement_iterations")
    monotonic = _diagnostic(record, "source_residual_decreased_monotonically")
    if not isinstance(monotonic, (bool, np.bool_)):
        raise TypeError("source residual monotonic diagnostic must be bool")
    return GR0UniversalStageSnapshot(
        time=record.time,
        transaction_serial=transaction_serial,
        minimum_lapse=_diagnostic(record, "minimum_lapse"),
        minimum_radial_metric=_diagnostic(record, "minimum_radial_metric"),
        minimum_areal_radius_away_from_center=_diagnostic(
            record, "minimum_areal_radius_away_from_center"
        ),
        minimum_hat_lorentzian_margin=gr0_hat_lorentzian_margin(
            record.state,
            hat_normal_factor=hat_normal_factor,
        ),
        source_residual_infinity=_diagnostic(
            record, "source_residual_infinity"
        ),
        source_iterations=iterations,
        source_residual_decreased_monotonically=bool(monotonic),
        kinetic_condition_infinity=_diagnostic(
            record, "kinetic_condition_infinity"
        ),
        boundary_margin_over_required_buffer=boundary_margin_over_required_buffer,
    )


def failed_gr0_universal_premises(
    snapshot: GR0UniversalStageSnapshot,
    thresholds: GR0UniversalThresholds,
) -> tuple[str, ...]:
    """Return every failed universal premise in the inherited priority."""

    if not isinstance(snapshot, GR0UniversalStageSnapshot):
        raise TypeError("snapshot must be a GR0UniversalStageSnapshot")
    if not isinstance(thresholds, GR0UniversalThresholds):
        raise TypeError("thresholds must be GR0UniversalThresholds")
    checks = {
        "nonpositive_lapse": snapshot.minimum_lapse <= 0.0,
        "nonpositive_radial_metric": snapshot.minimum_radial_metric <= 0.0,
        "nonpositive_areal_radius_away_from_center": (
            snapshot.minimum_areal_radius_away_from_center <= 0.0
        ),
        "hat_cone_not_Lorentzian": snapshot.minimum_hat_lorentzian_margin <= 0.0,
        "newton_residual_limit": (
            snapshot.source_residual_infinity >= thresholds.source_residual_maximum
        ),
        "newton_iteration_limit": (
            snapshot.source_iterations > thresholds.source_iteration_maximum
        ),
        "newton_residual_not_monotonic": (
            not snapshot.source_residual_decreased_monotonically
        ),
        "kinetic_condition_limit": (
            snapshot.kinetic_condition_infinity
            >= thresholds.kinetic_condition_maximum
        ),
        "boundary_causal_buffer": (
            snapshot.boundary_margin_over_required_buffer <= 0.0
        ),
    }
    if tuple(checks) != GR0_RUNTIME_STOP_PRIORITY:
        raise RuntimeError("GR-0 stop priority differs from PROTO5 inheritance")
    return tuple(name for name in GR0_RUNTIME_STOP_PRIORITY if checks[name])


@dataclass(frozen=True, slots=True)
class GR0RuntimeMonitorState:
    accepted_stage_count: int = 0
    last_transaction_serial: int = -1
    last_accepted_time: float = 0.0
    first_failed_premise: str | None = None
    first_failed_transaction_serial: int | None = None
    first_failed_time: float | None = None


class GR0RuntimeStop(RuntimeError):
    """Typed scientific stop preserving the first universal failure."""

    def __init__(
        self,
        state: GR0RuntimeMonitorState,
        failures: Sequence[str],
    ) -> None:
        self.state = state
        self.failures = tuple(failures)
        self.reason = state.first_failed_premise
        super().__init__(f"GR-0 runtime stopped at {self.reason}")


class CFLRetryRequired(RuntimeError):
    """Retryable numerical refusal; never a physical or scientific stop."""

    def __init__(self, *, observed_ratio: float, maximum_ratio: float) -> None:
        self.observed_ratio = observed_ratio
        self.maximum_ratio = maximum_ratio
        super().__init__(
            "proposed GR-0 step exceeds the frozen CFL ratio and must be retried"
        )


@dataclass(frozen=True, slots=True)
class GR0RuntimeTransactionReceipt:
    accepted_stage_count: int
    first_transaction_serial: int
    last_transaction_serial: int
    interval_speed_upper: float
    courant_ratio: float
    remaining_causal_buffer: float
    strict_margin_over_required_buffer: float

    def as_mapping(self) -> Mapping[str, object]:
        return {
            "accepted_stage_count": self.accepted_stage_count,
            "first_transaction_serial": self.first_transaction_serial,
            "last_transaction_serial": self.last_transaction_serial,
            "interval_speed_upper": self.interval_speed_upper,
            "courant_ratio": self.courant_ratio,
            "remaining_causal_buffer": self.remaining_causal_buffer,
            "strict_margin_over_required_buffer": self.strict_margin_over_required_buffer,
        }


class GR0RuntimeStageTransaction:
    """Atomic nine-stop, CFL, and BND2 guard for one GR-0 grid."""

    def __init__(
        self,
        *,
        thresholds: GR0UniversalThresholds,
        causal_state: CausalBudgetState,
        boundary_geometry: BoundaryGeometry,
        grid_spacing: Real,
        cfl_maximum: Real,
        hat_normal_factor: Real = 9.0,
    ) -> None:
        if not isinstance(thresholds, GR0UniversalThresholds):
            raise TypeError("thresholds must be GR0UniversalThresholds")
        if not isinstance(causal_state, CausalBudgetState):
            raise TypeError("causal_state must be CausalBudgetState")
        if not isinstance(boundary_geometry, BoundaryGeometry):
            raise TypeError("boundary_geometry must be BoundaryGeometry")
        self.thresholds = thresholds
        self.causal_state = causal_state
        self.boundary_geometry = boundary_geometry
        self.grid_spacing = _finite("grid_spacing", grid_spacing, positive=True)
        self.cfl_maximum = _finite("cfl_maximum", cfl_maximum, positive=True)
        self.hat_normal_factor = _finite(
            "hat_normal_factor", hat_normal_factor, positive=True
        )
        if self.hat_normal_factor <= 1.0:
            raise ValueError("hat_normal_factor must exceed one")
        self.state = GR0RuntimeMonitorState()

    @staticmethod
    def _stage_speed(record: StageRecord) -> float:
        return _finite(
            "coordinate_speed_upper",
            _diagnostic(record, "coordinate_speed_upper"),
        )

    def __call__(self, proposal: StepProposal) -> Mapping[str, object]:
        if not isinstance(proposal, StepProposal):
            raise TypeError("proposal must be a StepProposal")
        if not proposal.stages or proposal.stages[-1].stage_name != "candidate_endpoint":
            raise ValueError("GR-0 transaction requires an explicit candidate endpoint")
        if self.state.first_failed_premise is not None:
            raise GR0RuntimeStop(self.state, (self.state.first_failed_premise,))

        speeds = tuple(self._stage_speed(record) for record in proposal.stages)
        boundary = assess_causal_step(
            self.causal_state,
            self.boundary_geometry,
            trial_time=proposal.final_time,
            stage_coordinate_speed_uppers=speeds[:-1],
            candidate_endpoint_coordinate_speed_upper=speeds[-1],
        )
        step_size = proposal.final_time - proposal.initial_time
        courant = step_size * boundary.interval_speed_upper / self.grid_spacing
        tolerance = 4096.0 * np.finfo(np.float64).eps * max(
            1.0, self.cfl_maximum
        )
        if courant > self.cfl_maximum + tolerance:
            raise CFLRetryRequired(
                observed_ratio=courant,
                maximum_ratio=self.cfl_maximum,
            )

        first_serial = self.state.last_transaction_serial + 1
        for offset, record in enumerate(proposal.stages):
            serial = first_serial + offset
            snapshot = gr0_universal_stage_snapshot(
                record,
                transaction_serial=serial,
                boundary_margin_over_required_buffer=(
                    boundary.strict_margin_over_required_buffer
                ),
                hat_normal_factor=self.hat_normal_factor,
            )
            failures = failed_gr0_universal_premises(snapshot, self.thresholds)
            if failures:
                self.state = GR0RuntimeMonitorState(
                    accepted_stage_count=self.state.accepted_stage_count,
                    last_transaction_serial=self.state.last_transaction_serial,
                    last_accepted_time=self.state.last_accepted_time,
                    first_failed_premise=failures[0],
                    first_failed_transaction_serial=serial,
                    first_failed_time=record.time,
                )
                raise GR0RuntimeStop(self.state, failures)

        if not boundary.admissible:
            raise RuntimeError(
                "BND2 rejected a proposal not caught by the universal boundary stop"
            )
        last_serial = first_serial + len(proposal.stages) - 1
        self.state = GR0RuntimeMonitorState(
            accepted_stage_count=self.state.accepted_stage_count
            + len(proposal.stages),
            last_transaction_serial=last_serial,
            last_accepted_time=proposal.final_time,
        )
        self.causal_state = boundary.candidate_state
        return GR0RuntimeTransactionReceipt(
            accepted_stage_count=len(proposal.stages),
            first_transaction_serial=first_serial,
            last_transaction_serial=last_serial,
            interval_speed_upper=boundary.interval_speed_upper,
            courant_ratio=courant,
            remaining_causal_buffer=boundary.remaining_causal_buffer,
            strict_margin_over_required_buffer=(
                boundary.strict_margin_over_required_buffer
            ),
        ).as_mapping()


@dataclass(frozen=True, slots=True)
class GR0CommonEventAssessment:
    method: str
    coordinate_time: float
    point_counts: tuple[int, ...]
    constraint_admission: SemidiscreteConstraintAdmission
    spatial_spectral_admission: NestedSpectralAdmission
    maximum_spatial_round_trip_interpolation_infinity: float
    admission_passed: bool


def proto5_gr0_common_event(
    states: Sequence[EvolutionState],
    grids: Sequence[UniformRadialGrid],
    *,
    method: str,
    coordinate_time: Real,
    cutoff: Real = 16.0,
    measurement_radius_maximum: Real = 24.0,
    taper_fraction: Real = 1.0 / 8.0,
    spectral_thresholds: SpectralThresholds | None = None,
) -> GR0CommonEventAssessment:
    """Compose PROTO5 constraints and spatial spectra at one common event."""

    if method not in PROTO5_METHOD_CONTRACT:
        raise ValueError("unknown PROTO5 GR-0 method")
    records = tuple(states)
    meshes = tuple(grids)
    if len(records) != 3 or len(meshes) != 3:
        raise ValueError("PROTO5 common events require exactly three grids")
    if any(
        not isinstance(state, EvolutionState)
        or not isinstance(grid, UniformRadialGrid)
        or state.shape != (grid.point_count, 6)
        for state, grid in zip(records, meshes, strict=True)
    ):
        raise ValueError("PROTO5 common-event state/grid pairs differ")
    counts = tuple(grid.point_count for grid in meshes)
    if counts != tuple(sorted(counts)) or len(set(counts)) != 3:
        raise ValueError("PROTO5 common-event grids must increase strictly")
    contract = PROTO5_METHOD_CONTRACT[method]
    order = int(contract["spatial_order"])
    time = _finite("coordinate_time", coordinate_time)
    constraint_snapshots = tuple(
        gr0_semidiscrete_constraint_snapshot(
            state,
            grid,
            coordinate_time=time,
            diagnostic_spatial_order=order,
        )
        for state, grid in zip(records, meshes, strict=True)
    )
    constraint = proto5_semidiscrete_constraint_admission(
        constraint_snapshots,
        method=method,
        coarsest_guard_maximum=contract["coarsest_constraint_guard"],
        finest_guard_maximum=contract["finest_constraint_guard"],
        minimum_finest_pair_order=1.5,
    )

    fraction = _finite("taper_fraction", taper_fraction, positive=True)
    if fraction >= 0.5:
        raise ValueError("taper_fraction must be below one half")
    budgets: list[Mapping[str, SpectralPowerBudget]] = []
    interpolation_maximum = 0.0
    for state, grid in zip(records, meshes, strict=True):
        profile = proper_radial_profile(
            state.u,
            state.q,
            grid.coordinates,
            cutoff=cutoff,
            measurement_radius_maximum=measurement_radius_maximum,
        )
        proper = profile.proper_coordinates
        window = compact_vacuum_buffer_window(
            proper,
            support_minimum=float(proper[0]),
            support_maximum=float(proper[-1]),
            taper_width=float((proper[-1] - proper[0]) * fraction),
        )
        budgets.append(
            spectral_field_budgets(
                profile.deviations,
                proper,
                cutoff=cutoff,
                window=window,
                thresholds=spectral_thresholds,
            )
        )
        interpolation_maximum = max(
            interpolation_maximum,
            float(
                np.max(
                    profile.round_trip_interpolation_infinity,
                    initial=0.0,
                )
            ),
        )
    spatial = nested_spectral_admission(
        counts,
        budgets,
        thresholds=spectral_thresholds,
    )
    return GR0CommonEventAssessment(
        method=method,
        coordinate_time=time,
        point_counts=counts,
        constraint_admission=constraint,
        spatial_spectral_admission=spatial,
        maximum_spatial_round_trip_interpolation_infinity=interpolation_maximum,
        admission_passed=(constraint.admission_passed and spatial.admission_passed),
    )


@dataclass(frozen=True, slots=True)
class TemporalSpectralAdmission:
    point_counts: tuple[int, ...]
    sample_count: int
    tracer_count: int
    field_names: tuple[str, ...]
    maximum_round_trip_interpolation_infinity: float
    nested_admission: NestedSpectralAdmission
    admission_passed: bool


@dataclass(frozen=True, slots=True)
class CalibrationTrappedAssessment:
    """Outcome-neutral robust trapped-sphere sign test at one common event."""

    coordinate_time: float
    point_counts: tuple[int, ...]
    primary_trapped_scores: tuple[float, ...]
    comparator_trapped_scores: tuple[float, ...]
    primary_Richardson_error: float
    comparator_Richardson_error: float
    fine_cross_method_error: float
    roundoff_error: float
    combined_error: float
    positive_margin_over_combined_error_factor: float | None
    both_common_event_admissions_passed: bool
    trapped_sign_passed: bool


def _trapped_scores_on_common_nodes(
    states: Sequence[EvolutionState],
    grids: Sequence[UniformRadialGrid],
    *,
    measurement_radius_maximum: float,
) -> tuple[float, ...]:
    counts = tuple(grid.point_count for grid in grids)
    coarse_intervals = counts[0] - 1
    scores: list[float] = []
    for state, grid in zip(states, grids, strict=True):
        intervals = grid.point_count - 1
        if intervals % coarse_intervals != 0:
            raise ValueError("trapped-score grids are not exactly nested")
        stride = intervals // coarse_intervals
        indices = np.arange(0, grid.point_count, stride, dtype=np.int64)
        radii = grid.coordinates[indices]
        selected = indices[(radii > 0.0) & (radii <= measurement_radius_maximum)]
        if selected.size < 1:
            raise ValueError("trapped-score measurement interval is empty")
        observed = radial_null_observables(state)
        # A positive value is the minimum distance of both expansions below
        # zero at the best common-grid sphere.  This retains the frozen null
        # orientation and never substitutes compactness alone for both signs.
        point_scores = np.minimum(
            -observed.theta_plus[selected],
            -observed.theta_minus[selected],
        )
        scores.append(float(np.max(point_scores, initial=-np.inf)))
    return tuple(scores)


def proto5_calibration_trapped_assessment(
    *,
    primary_states: Sequence[EvolutionState],
    comparator_states: Sequence[EvolutionState],
    grids: Sequence[UniformRadialGrid],
    primary_common_event: GR0CommonEventAssessment,
    comparator_common_event: GR0CommonEventAssessment,
    measurement_radius_maximum: Real = 24.0,
    minimum_observed_order: Real = 1.5,
    positive_margin_factor: Real = 4.0,
    roundoff_operation_budget: int = 4096,
) -> CalibrationTrappedAssessment:
    """Require two methods and three grids to support both trapped signs.

    Constraint and spectral admission are hard vetoes.  They are deliberately
    not converted into expansion units by an arbitrary multiplier.  The sign
    error is instead the additive sum of primary and comparator finest-pair
    Richardson estimates, their fine-grid method disagreement, and an
    explicit binary64 enclosure.  This calibration-only error contract is not
    the still-missing constraint-to-Raychaudhuri stability map for DEF1.
    """

    p_states = tuple(primary_states)
    c_states = tuple(comparator_states)
    meshes = tuple(grids)
    if len(p_states) != 3 or len(c_states) != 3 or len(meshes) != 3:
        raise ValueError("trapped assessment requires two methods on three grids")
    if (
        primary_common_event.method != "RK4"
        or comparator_common_event.method != "SSPRK3"
    ):
        raise ValueError("trapped assessment method ownership differs")
    time = primary_common_event.coordinate_time
    if np.float64(time).tobytes() != np.float64(
        comparator_common_event.coordinate_time
    ).tobytes():
        raise ValueError("trapped assessment events are not bitwise time aligned")
    counts = tuple(grid.point_count for grid in meshes)
    if (
        primary_common_event.point_counts != counts
        or comparator_common_event.point_counts != counts
    ):
        raise ValueError("trapped assessment grid lineage differs")
    maximum = _finite(
        "measurement_radius_maximum",
        measurement_radius_maximum,
        positive=True,
    )
    order = _finite("minimum_observed_order", minimum_observed_order, positive=True)
    factor = _finite("positive_margin_factor", positive_margin_factor, positive=True)
    budget = _integer(
        "roundoff_operation_budget", roundoff_operation_budget, minimum=1
    )
    primary_scores = _trapped_scores_on_common_nodes(
        p_states,
        meshes,
        measurement_radius_maximum=maximum,
    )
    comparator_scores = _trapped_scores_on_common_nodes(
        c_states,
        meshes,
        measurement_radius_maximum=maximum,
    )
    refinement = (counts[-1] - 1) / (counts[-2] - 1)
    if refinement <= 1.0:
        raise ValueError("trapped assessment refinement ratio must exceed one")
    denominator = refinement**order - 1.0
    primary_error = abs(primary_scores[-1] - primary_scores[-2]) / denominator
    comparator_error = (
        abs(comparator_scores[-1] - comparator_scores[-2]) / denominator
    )
    method_error = abs(primary_scores[-1] - comparator_scores[-1])
    scale = max(
        1.0,
        *(abs(value) for value in (*primary_scores, *comparator_scores)),
    )
    roundoff = budget * np.finfo(np.float64).eps * scale
    combined = primary_error + comparator_error + method_error + roundoff
    fine_margin = min(primary_scores[-1], comparator_scores[-1])
    admissions = (
        primary_common_event.admission_passed
        and comparator_common_event.admission_passed
    )
    passed = admissions and fine_margin > factor * combined
    return CalibrationTrappedAssessment(
        coordinate_time=time,
        point_counts=counts,
        primary_trapped_scores=primary_scores,
        comparator_trapped_scores=comparator_scores,
        primary_Richardson_error=primary_error,
        comparator_Richardson_error=comparator_error,
        fine_cross_method_error=method_error,
        roundoff_error=roundoff,
        combined_error=combined,
        positive_margin_over_combined_error_factor=(
            fine_margin / combined if combined > 0.0 else None
        ),
        both_common_event_admissions_passed=admissions,
        trapped_sign_passed=passed,
    )


def proto5_temporal_spectral_admission(
    *,
    point_counts: Sequence[int],
    proper_times_by_resolution: Sequence[object],
    field_histories_by_resolution: Sequence[object],
    field_names: Sequence[str] = PROTO4_SPECTRAL_FIELD_ORDER,
    cutoff: Real = 16.0,
    taper_fraction: Real = 1.0 / 8.0,
    thresholds: SpectralThresholds | None = None,
) -> TemporalSpectralAdmission:
    """Evaluate nested causal-past spectra on matched normal-flow tracers.

    ``proper_times`` has shape ``(samples, tracers)`` and ``field_histories``
    has shape ``(samples, tracers, fields)``.  The runner owns tracer
    transport and must provide identical tracer labels/order at every
    resolution.  Each nonuniform proper-time column is deterministically
    resampled before the FFT, and the reverse interpolation mismatch is
    retained as a public nonnegative error component.
    """

    counts = tuple(point_counts)
    times_records = tuple(proper_times_by_resolution)
    field_records = tuple(field_histories_by_resolution)
    names = tuple(field_names)
    if (
        len(counts) != 3
        or len(times_records) != 3
        or len(field_records) != 3
        or counts != tuple(sorted(counts))
        or len(set(counts)) != 3
    ):
        raise ValueError("temporal admission requires three nested resolutions")
    if names != PROTO4_SPECTRAL_FIELD_ORDER:
        raise ValueError("temporal field order differs from PROTO5 inheritance")
    fraction = _finite("taper_fraction", taper_fraction, positive=True)
    if fraction >= 0.5:
        raise ValueError("taper_fraction must be below one half")

    flattened_budgets: list[dict[str, SpectralPowerBudget]] = []
    reference_shape: tuple[int, int, int] | None = None
    maximum_round_trip = 0.0
    for resolution, (times_value, fields_value) in enumerate(
        zip(times_records, field_records, strict=True)
    ):
        times = np.asarray(times_value, dtype=np.float64)
        fields = np.asarray(fields_value, dtype=np.float64)
        if (
            times.ndim != 2
            or fields.ndim != 3
            or fields.shape[:2] != times.shape
            or fields.shape[2] != len(names)
            or not np.all(np.isfinite(times))
            or not np.all(np.isfinite(fields))
        ):
            raise ValueError("temporal history arrays have incompatible shapes")
        if times.shape[0] < 64 or np.any(np.diff(times, axis=0) <= 0.0):
            raise ValueError("temporal spectra require 64 increasing past samples")
        shape = fields.shape
        if reference_shape is None:
            reference_shape = shape
        elif shape != reference_shape:
            raise ValueError("temporal histories must share samples and tracers")

        budgets: dict[str, SpectralPowerBudget] = {}
        for tracer in range(times.shape[1]):
            proper = times[:, tracer]
            uniform = np.linspace(
                proper[0], proper[-1], proper.size, dtype=np.float64
            )
            width = fraction * (uniform[-1] - uniform[0])
            window = compact_vacuum_buffer_window(
                uniform,
                support_minimum=float(uniform[0]),
                support_maximum=float(uniform[-1]),
                taper_width=float(width),
            )
            for field, name in enumerate(names):
                values = fields[:, tracer, field]
                resampled = np.interp(uniform, proper, values)
                reversed_values = np.interp(proper, uniform, resampled)
                maximum_round_trip = max(
                    maximum_round_trip,
                    float(np.max(np.abs(reversed_values - values), initial=0.0)),
                )
                key = f"tracer_{tracer:04d}:{name}"
                budgets[key] = windowed_spectral_power_budget(
                    resampled,
                    uniform,
                    cutoff=cutoff,
                    window=window,
                    thresholds=thresholds,
                )
        flattened_budgets.append(budgets)

    nested = nested_spectral_admission(
        counts,
        flattened_budgets,
        thresholds=thresholds,
    )
    assert reference_shape is not None
    return TemporalSpectralAdmission(
        point_counts=counts,
        sample_count=reference_shape[0],
        tracer_count=reference_shape[1],
        field_names=names,
        maximum_round_trip_interpolation_infinity=maximum_round_trip,
        nested_admission=nested,
        admission_passed=nested.admission_passed,
    )


__all__ = [
    "CalibrationTrappedAssessment",
    "CFLRetryRequired",
    "EventLogRecoveryPlan",
    "GR0CommonEventAssessment",
    "GR0RuntimeMonitorState",
    "GR0RuntimeStageTransaction",
    "GR0RuntimeStop",
    "GR0RuntimeTransactionReceipt",
    "GR0UniversalStageSnapshot",
    "GR0UniversalThresholds",
    "GR0_RUNTIME_STOP_PRIORITY",
    "PROTO5_METHOD_CONTRACT",
    "TemporalSpectralAdmission",
    "failed_gr0_universal_premises",
    "gr0_hat_lorentzian_margin",
    "gr0_universal_stage_snapshot",
    "proto5_gr0_common_event",
    "proto5_calibration_trapped_assessment",
    "proto5_event_log_recovery_plan",
    "proto5_temporal_spectral_admission",
]
