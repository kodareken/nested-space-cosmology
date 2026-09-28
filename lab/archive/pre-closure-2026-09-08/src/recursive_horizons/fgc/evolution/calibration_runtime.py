"""Prospective semidiscrete controls for fresh GR-0 calibration.

This module is a successor adapter.  It does not alter the hash-bound ID2,
HLT2, GR-0 source, or numerical-engine implementations.  Its purpose is to
make the continuum-to-grid map and the floating-point interpretation of the
PROTO4 common-event constraint test explicit before a fresh trajectory is
opened.

The physical data ``u`` and ``p`` remain untouched.  Only the auxiliary
first-order variable ``q`` is projected to the derivative operator used by
the selected numerical method.  Physical and gauge constraints are then
evaluated from the unredefined GR-0 metric residual, while reduction
constraints remain ``q-D_h u``.  Raw normalized residuals always decide the
magnitude guards.  A declared floating-operation enclosure may classify a
component as numerically indistinguishable from exact zero for convergence
order only; its nonnegative size remains public and is never subtracted from
an observable margin.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite, log
from numbers import Real
from typing import Mapping, Sequence

import numpy as np

from .gr0_calibration import GR0GridInitialData, Q_CENTER_PARITIES
from .gr0_direct_source import gr0_ref1_residual_batch
from .numerical_engine import EvolutionState, SBPFirstDerivative, UniformRadialGrid
from .proto4_admission import (
    PROTO4_CONSTRAINT_ORDER,
    ConstraintNorms,
    normalized_constraint_norms,
)
from .vectorized_source import ADM_CENTER_PARITIES


GR0_UNIVERSAL_RUNTIME_STOP_IDS = (
    "nonpositive_lapse",
    "nonpositive_radial_metric",
    "nonpositive_areal_radius_away_from_center",
    "hat_cone_not_Lorentzian",
    "newton_residual_limit",
    "newton_iteration_limit",
    "newton_residual_not_monotonic",
    "kinetic_condition_limit",
    "boundary_causal_buffer",
)
CONSTRAINT_ROUNDOFF_OPERATION_BUDGET = 4096


def _finite(name: str, value: Real, *, positive: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a real scalar")
    answer = float(value)
    if not isfinite(answer):
        raise ValueError(f"{name} must be finite")
    if positive and answer <= 0.0:
        raise ValueError(f"{name} must be positive")
    return answer


def _freeze(name: str, value: object, *, ndim: int) -> np.ndarray:
    answer = np.asarray(value, dtype=np.float64)
    if answer.ndim != ndim or not np.all(np.isfinite(answer)):
        raise ValueError(f"{name} must be a finite rank-{ndim} binary64 array")
    answer = np.ascontiguousarray(answer).copy()
    answer.setflags(write=False)
    return answer


def project_gr0_semidiscrete_state(
    initial: GR0GridInitialData,
    *,
    spatial_order: int,
) -> EvolutionState:
    """Preserve physical initial data and impose ``q=D_h u`` exactly."""

    if not isinstance(initial, GR0GridInitialData):
        raise TypeError("initial must be GR0GridInitialData")
    derivative = SBPFirstDerivative(initial.grid, spatial_order)
    q = derivative.differentiate(
        initial.state.u,
        center_parities=ADM_CENTER_PARITIES,
    )
    return EvolutionState(initial.state.u, initial.state.p, q)


@dataclass(frozen=True, slots=True)
class SemidiscreteConstraintSnapshot:
    """Raw and roundoff-enclosed constraint evidence at one grid event."""

    point_count: int
    coordinate_time: float
    diagnostic_spatial_order: int
    raw_norms: ConstraintNorms
    effective_norms_for_order_only: ConstraintNorms
    normalized_roundoff_zero_enclosure: np.ndarray
    residuals: np.ndarray
    unredefined_term_contributions: np.ndarray

    def __post_init__(self) -> None:
        if (
            isinstance(self.point_count, bool)
            or not isinstance(self.point_count, int)
            or self.point_count < 9
        ):
            raise ValueError("point_count must be an integer of at least nine")
        object.__setattr__(
            self,
            "coordinate_time",
            _finite("coordinate_time", self.coordinate_time),
        )
        if self.diagnostic_spatial_order not in {2, 4}:
            raise ValueError("diagnostic_spatial_order must be two or four")
        for name in ("raw_norms", "effective_norms_for_order_only"):
            value = getattr(self, name)
            if not isinstance(value, ConstraintNorms):
                raise TypeError(f"{name} must be ConstraintNorms")
            if value.sample_count != self.point_count:
                raise ValueError(f"{name} sample count differs")
            if value.component_names != PROTO4_CONSTRAINT_ORDER:
                raise ValueError(f"{name} component order differs")
        enclosure = _freeze(
            "normalized_roundoff_zero_enclosure",
            self.normalized_roundoff_zero_enclosure,
            ndim=1,
        )
        if enclosure.shape != (len(PROTO4_CONSTRAINT_ORDER),) or np.any(enclosure < 0.0):
            raise ValueError("roundoff enclosure shape or sign differs")
        residuals = _freeze("residuals", self.residuals, ndim=2)
        terms = _freeze(
            "unredefined_term_contributions",
            self.unredefined_term_contributions,
            ndim=3,
        )
        if residuals.shape != (self.point_count, len(PROTO4_CONSTRAINT_ORDER)):
            raise ValueError("constraint residual shape differs")
        if terms.shape[:2] != residuals.shape or terms.shape[2] < 1:
            raise ValueError("constraint term-contribution shape differs")
        object.__setattr__(self, "normalized_roundoff_zero_enclosure", enclosure)
        object.__setattr__(self, "residuals", residuals)
        object.__setattr__(self, "unredefined_term_contributions", terms)


def gr0_semidiscrete_constraint_snapshot(
    state: EvolutionState,
    grid: UniformRadialGrid,
    *,
    coordinate_time: Real,
    diagnostic_spatial_order: int,
    planck_mass: Real = 2.0,
    scalar_mass: Real = 3.0,
    quartic_coupling: Real = 0.5,
    length_unit: Real = 4.0,
    roundoff_operation_budget: int = CONSTRAINT_ROUNDOFF_OPERATION_BUDGET,
) -> SemidiscreteConstraintSnapshot:
    """Evaluate ten constraints and their declared normalization at one event."""

    if not isinstance(state, EvolutionState) or state.shape != (grid.point_count, 6):
        raise ValueError("state and grid must define one six-field radial slice")
    if (
        isinstance(roundoff_operation_budget, bool)
        or not isinstance(roundoff_operation_budget, int)
        or roundoff_operation_budget < 1
    ):
        raise ValueError("roundoff_operation_budget must be a positive integer")
    derivative = SBPFirstDerivative(grid, diagnostic_spatial_order)
    differentiated_u = derivative.differentiate(
        state.u,
        center_parities=ADM_CENTER_PARITIES,
    )
    p_r = derivative.differentiate(
        state.p,
        center_parities=ADM_CENTER_PARITIES,
    )
    q_r = derivative.differentiate(
        state.q,
        center_parities=Q_CENTER_PARITIES,
    )
    source = gr0_ref1_residual_batch(
        state.u[1:],
        state.p[1:],
        state.q[1:],
        np.zeros((1, grid.point_count - 1, 6), dtype=np.float64),
        p_r[1:],
        q_r[1:],
        grid.coordinates[1:],
        planck_mass=planck_mass,
        scalar_mass=scalar_mass,
        quartic_coupling=quartic_coupling,
    )

    residuals = np.zeros((grid.point_count, len(PROTO4_CONSTRAINT_ORDER)))
    terms = np.zeros((grid.point_count, len(PROTO4_CONSTRAINT_ORDER), 3))
    residuals[1:, 0] = source.hamiltonian_constraint[0]
    residuals[1:, 1] = source.momentum_constraint[0]
    residuals[1:, 2:4] = source.gauge_constraint[0, :, :2]
    residuals[:, 4:] = state.q - differentiated_u

    metric = source.unredefined_metric_residual[0]
    shift = state.u[1:, 1]
    terms[1:, 0, 0] = metric[:, 0, 0]
    terms[1:, 0, 1] = -2.0 * shift * metric[:, 0, 1]
    terms[1:, 0, 2] = shift**2 * metric[:, 1, 1]
    terms[1:, 1, 0] = metric[:, 0, 1]
    terms[1:, 1, 1] = -shift * metric[:, 1, 1]
    terms[1:, 2:4, 0] = source.gauge_constraint[0, :, :2]
    terms[:, 4:, 0] = state.q
    terms[:, 4:, 1] = -differentiated_u

    tolerance = 2048.0 * np.finfo(np.float64).eps
    reconstructed = np.sum(terms, axis=2)
    if not np.allclose(reconstructed, residuals, rtol=tolerance, atol=tolerance):
        raise RuntimeError("constraint term decomposition does not reconstruct residuals")
    raw = normalized_constraint_norms(
        residuals,
        terms,
        planck_mass=planck_mass,
        length_unit=length_unit,
    )
    normalized_enclosure = np.full(
        len(PROTO4_CONSTRAINT_ORDER),
        roundoff_operation_budget * np.finfo(np.float64).eps,
        dtype=np.float64,
    )
    effective_components = np.maximum(
        raw.component_infinity - normalized_enclosure,
        0.0,
    )
    effective = ConstraintNorms(
        sample_count=raw.sample_count,
        component_names=raw.component_names,
        component_infinity=effective_components,
        global_infinity=float(np.max(effective_components, initial=0.0)),
        minimum_denominator=raw.minimum_denominator,
        maximum_denominator=raw.maximum_denominator,
    )
    return SemidiscreteConstraintSnapshot(
        point_count=grid.point_count,
        coordinate_time=coordinate_time,
        diagnostic_spatial_order=diagnostic_spatial_order,
        raw_norms=raw,
        effective_norms_for_order_only=effective,
        normalized_roundoff_zero_enclosure=normalized_enclosure,
        residuals=residuals,
        unredefined_term_contributions=terms,
    )


@dataclass(frozen=True, slots=True)
class SemidiscreteConstraintAdmission:
    method: str
    point_counts: tuple[int, ...]
    coordinate_time: float
    raw_global_norms: tuple[float, ...]
    component_finest_pair_orders: Mapping[str, float | None]
    component_status: Mapping[str, str]
    minimum_finite_finest_pair_order: float | None
    coarsest_guard_passed: bool
    finest_guard_passed: bool
    monotone_refinement_passed: bool
    finest_pair_order_passed: bool
    common_event_alignment_passed: bool
    admission_passed: bool


def _bitwise_equal_binary64(values: Sequence[float]) -> bool:
    if not values:
        return False
    first = np.float64(values[0]).tobytes()
    return all(np.float64(value).tobytes() == first for value in values[1:])


def proto5_semidiscrete_constraint_admission(
    samples: Sequence[SemidiscreteConstraintSnapshot],
    *,
    method: str,
    coarsest_guard_maximum: Real,
    finest_guard_maximum: Real,
    minimum_finest_pair_order: Real = 1.5,
) -> SemidiscreteConstraintAdmission:
    """Admit raw magnitude plus monotone, finest-pair asymptotic convergence."""

    records = tuple(samples)
    if len(records) != 3 or any(
        not isinstance(item, SemidiscreteConstraintSnapshot) for item in records
    ):
        raise ValueError("semidiscrete admission requires exactly three snapshots")
    counts = tuple(item.point_count for item in records)
    if tuple(sorted(counts)) != counts or len(set(counts)) != 3:
        raise ValueError("semidiscrete point counts must increase strictly")
    intervals = tuple(count - 1 for count in counts)
    refinement = tuple(
        right / left
        for left, right in zip(intervals[:-1], intervals[1:], strict=True)
    )
    if len(set(refinement)) != 1 or refinement[0] <= 1.0:
        raise ValueError("semidiscrete grids must share one refinement ratio")
    coarse_guard = _finite(
        "coarsest_guard_maximum", coarsest_guard_maximum, positive=True
    )
    fine_guard = _finite(
        "finest_guard_maximum", finest_guard_maximum, positive=True
    )
    order_minimum = _finite(
        "minimum_finest_pair_order", minimum_finest_pair_order, positive=True
    )
    times = tuple(item.coordinate_time for item in records)
    aligned = _bitwise_equal_binary64(times)
    component_orders: dict[str, float | None] = {}
    statuses: dict[str, str] = {}
    monotone = True
    order_pass = True
    finite_orders: list[float] = []
    for component, name in enumerate(PROTO4_CONSTRAINT_ORDER):
        values = tuple(
            float(item.effective_norms_for_order_only.component_infinity[component])
            for item in records
        )
        if any(
            right > left
            for left, right in zip(values[:-1], values[1:], strict=True)
        ):
            monotone = False
        medium, fine = values[-2:]
        if medium == 0.0 and fine == 0.0:
            component_orders[name] = None
            statuses[name] = "roundoff_enclosed_zero_pair"
        elif medium > 0.0 and fine == 0.0:
            component_orders[name] = None
            statuses[name] = "resolved_to_roundoff_enclosure"
        elif medium == 0.0 and fine > 0.0:
            component_orders[name] = None
            statuses[name] = "nonzero_reappeared_after_roundoff_enclosure"
            order_pass = False
        else:
            order = log(medium / fine) / log(refinement[-1])
            component_orders[name] = order
            statuses[name] = "finite_order"
            finite_orders.append(order)
            if order < order_minimum:
                order_pass = False
    raw_globals = tuple(item.raw_norms.global_infinity for item in records)
    coarse_pass = raw_globals[0] < coarse_guard
    fine_pass = raw_globals[-1] < fine_guard
    passed = aligned and coarse_pass and fine_pass and monotone and order_pass
    return SemidiscreteConstraintAdmission(
        method=method,
        point_counts=counts,
        coordinate_time=times[0],
        raw_global_norms=raw_globals,
        component_finest_pair_orders=component_orders,
        component_status=statuses,
        minimum_finite_finest_pair_order=(min(finite_orders) if finite_orders else None),
        coarsest_guard_passed=coarse_pass,
        finest_guard_passed=fine_pass,
        monotone_refinement_passed=monotone,
        finest_pair_order_passed=order_pass,
        common_event_alignment_passed=aligned,
        admission_passed=passed,
    )


__all__ = [
    "CONSTRAINT_ROUNDOFF_OPERATION_BUDGET",
    "GR0_UNIVERSAL_RUNTIME_STOP_IDS",
    "SemidiscreteConstraintAdmission",
    "SemidiscreteConstraintSnapshot",
    "gr0_semidiscrete_constraint_snapshot",
    "project_gr0_semidiscrete_state",
    "proto5_semidiscrete_constraint_admission",
]
