"""Outcome-neutral diagnosis of the PROTO7 common-event admission stop.

This module does not alter a trajectory, source equation, numerical method,
threshold, or physical input.  It supplies two prospective admission
contracts that can be evaluated on an already committed common-event slice:

* constraint convergence is evaluated on the rows owned by the evolution
  operator, excluding the explicitly frozen outer rows and every row whose
  derivative stencil touches them; accumulated binary64 roundoff is enclosed
  for convergence-order classification only;
* the coarsest spectral grid remains a convergence witness, while the two
  finest grids must each satisfy every unchanged absolute spectral budget and
  all adjacent tail ratios must still contract by the unchanged factor.

Raw residuals and every per-grid spectrum remain public.  Neither adapter can
authorize a run or promote a physical claim.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite, log
from numbers import Real
from typing import Mapping, Sequence

import numpy as np

from .calibration_runtime import (
    CONSTRAINT_ROUNDOFF_OPERATION_BUDGET,
    gr0_semidiscrete_constraint_snapshot,
)
from .numerical_engine import EvolutionState, SBPFirstDerivative, UniformRadialGrid
from .proto4_admission import (
    PROTO4_CONSTRAINT_ORDER,
    ConstraintNorms,
    NestedSpectralAdmission,
    SpectralPowerBudget,
    SpectralThresholds,
    nested_spectral_admission,
    normalized_constraint_norms,
)


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


@dataclass(frozen=True, slots=True)
class OwnedConstraintSnapshot:
    """One raw constraint slice plus its explicit diagnostic ownership."""

    point_count: int
    coordinate_time: float
    diagnostic_spatial_order: int
    accepted_stage_count: int
    fixed_outer_rows: int
    derivative_stencil_reach_rows: int
    excluded_outer_rows: int
    last_owned_radius: float
    full_domain_norms: ConstraintNorms
    owned_domain_norms: ConstraintNorms
    effective_owned_norms_for_order_only: ConstraintNorms
    accumulated_roundoff_enclosure: float
    full_component_maximum_radii: np.ndarray
    owned_component_maximum_radii: np.ndarray

    def __post_init__(self) -> None:
        if (
            isinstance(self.point_count, bool)
            or not isinstance(self.point_count, int)
            or self.point_count < 9
        ):
            raise ValueError("point_count must be an integer of at least nine")
        object.__setattr__(
            self, "coordinate_time", _finite("coordinate_time", self.coordinate_time)
        )
        if self.diagnostic_spatial_order not in {2, 4}:
            raise ValueError("diagnostic_spatial_order must be two or four")
        for name in (
            "accepted_stage_count",
            "fixed_outer_rows",
            "derivative_stencil_reach_rows",
            "excluded_outer_rows",
        ):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value < 0:
                raise ValueError(f"{name} must be a nonnegative integer")
        if self.fixed_outer_rows < 1 or self.derivative_stencil_reach_rows < 1:
            raise ValueError("outer ownership requires positive fixed and stencil rows")
        if self.excluded_outer_rows != (
            self.fixed_outer_rows + self.derivative_stencil_reach_rows
        ):
            raise ValueError("excluded rows must equal fixed rows plus stencil reach")
        if self.excluded_outer_rows >= self.point_count // 2:
            raise ValueError("outer ownership removes too much of the grid")
        object.__setattr__(
            self, "last_owned_radius", _finite("last_owned_radius", self.last_owned_radius)
        )
        enclosure = _finite(
            "accumulated_roundoff_enclosure",
            self.accumulated_roundoff_enclosure,
            positive=True,
        )
        object.__setattr__(self, "accumulated_roundoff_enclosure", enclosure)
        for name in (
            "full_domain_norms",
            "owned_domain_norms",
            "effective_owned_norms_for_order_only",
        ):
            value = getattr(self, name)
            if not isinstance(value, ConstraintNorms):
                raise TypeError(f"{name} must be ConstraintNorms")
            if value.component_names != PROTO4_CONSTRAINT_ORDER:
                raise ValueError(f"{name} component order differs")
        if self.full_domain_norms.sample_count != self.point_count:
            raise ValueError("full-domain sample count differs")
        owned_count = self.point_count - self.excluded_outer_rows
        if (
            self.owned_domain_norms.sample_count != owned_count
            or self.effective_owned_norms_for_order_only.sample_count != owned_count
        ):
            raise ValueError("owned-domain sample count differs")
        for name in (
            "full_component_maximum_radii",
            "owned_component_maximum_radii",
        ):
            value = _freeze(name, getattr(self, name), ndim=1)
            if value.shape != (len(PROTO4_CONSTRAINT_ORDER),) or np.any(value < 0.0):
                raise ValueError(f"{name} shape or sign differs")
            object.__setattr__(self, name, value)


def owned_semidiscrete_constraint_snapshot(
    state: EvolutionState,
    grid: UniformRadialGrid,
    *,
    coordinate_time: Real,
    diagnostic_spatial_order: int,
    fixed_outer_rows: int,
    accepted_stage_count: int,
    planck_mass: Real = 2.0,
    scalar_mass: Real = 3.0,
    quartic_coupling: Real = 0.5,
    length_unit: Real = 4.0,
) -> OwnedConstraintSnapshot:
    """Evaluate raw constraints on full and evolution-owned grid domains."""

    if isinstance(fixed_outer_rows, bool) or not isinstance(fixed_outer_rows, int):
        raise TypeError("fixed_outer_rows must be an integer")
    if fixed_outer_rows < 1:
        raise ValueError("fixed_outer_rows must be positive")
    if isinstance(accepted_stage_count, bool) or not isinstance(accepted_stage_count, int):
        raise TypeError("accepted_stage_count must be an integer")
    if accepted_stage_count < 0:
        raise ValueError("accepted_stage_count must be nonnegative")
    derivative = SBPFirstDerivative(grid, diagnostic_spatial_order)
    reach = derivative.stencil_reach_intervals
    excluded = fixed_outer_rows + reach
    owned_stop = grid.point_count - excluded
    if owned_stop < 9:
        raise ValueError("owned constraint domain is too small")

    full = gr0_semidiscrete_constraint_snapshot(
        state,
        grid,
        coordinate_time=coordinate_time,
        diagnostic_spatial_order=diagnostic_spatial_order,
        planck_mass=planck_mass,
        scalar_mass=scalar_mass,
        quartic_coupling=quartic_coupling,
        length_unit=length_unit,
    )
    owned_residuals = full.residuals[:owned_stop]
    owned_terms = full.unredefined_term_contributions[:owned_stop]
    owned = normalized_constraint_norms(
        owned_residuals,
        owned_terms,
        planck_mass=planck_mass,
        length_unit=length_unit,
    )

    # The original 4096-operation enclosure covers one diagnostic evaluation.
    # q and D_h u have instead accumulated independently through every accepted
    # stage.  Multiplying by (1 + accepted_stage_count) is a deterministic upper
    # budget, used only to decide whether an observed order is distinguishable
    # from accumulated binary64 cancellation.  Raw magnitudes are untouched.
    enclosure = float(
        CONSTRAINT_ROUNDOFF_OPERATION_BUDGET
        * np.finfo(np.float64).eps
        * (1 + accepted_stage_count)
    )
    effective_components = np.maximum(
        owned.component_infinity - enclosure,
        0.0,
    )
    effective = ConstraintNorms(
        sample_count=owned.sample_count,
        component_names=owned.component_names,
        component_infinity=effective_components,
        global_infinity=float(np.max(effective_components, initial=0.0)),
        minimum_denominator=owned.minimum_denominator,
        maximum_denominator=owned.maximum_denominator,
    )

    mass = _finite("planck_mass", planck_mass, positive=True)
    length = _finite("length_unit", length_unit, positive=True)
    denominator = np.maximum(
        np.sum(np.abs(full.unredefined_term_contributions), axis=2),
        (mass * mass) / (length * length),
    )
    pointwise = np.abs(full.residuals) / denominator
    full_indices = np.argmax(pointwise, axis=0)
    owned_indices = np.argmax(pointwise[:owned_stop], axis=0)
    coordinates = grid.coordinates
    return OwnedConstraintSnapshot(
        point_count=grid.point_count,
        coordinate_time=float(coordinate_time),
        diagnostic_spatial_order=diagnostic_spatial_order,
        accepted_stage_count=accepted_stage_count,
        fixed_outer_rows=fixed_outer_rows,
        derivative_stencil_reach_rows=reach,
        excluded_outer_rows=excluded,
        last_owned_radius=float(coordinates[owned_stop - 1]),
        full_domain_norms=full.raw_norms,
        owned_domain_norms=owned,
        effective_owned_norms_for_order_only=effective,
        accumulated_roundoff_enclosure=enclosure,
        full_component_maximum_radii=coordinates[full_indices],
        owned_component_maximum_radii=coordinates[owned_indices],
    )


@dataclass(frozen=True, slots=True)
class OwnedConstraintAdmission:
    method: str
    point_counts: tuple[int, ...]
    coordinate_time: float
    excluded_outer_rows: tuple[int, ...]
    last_owned_radii: tuple[float, ...]
    accepted_stage_counts: tuple[int, ...]
    accumulated_roundoff_enclosures: tuple[float, ...]
    full_domain_component_infinity: Mapping[str, tuple[float, ...]]
    owned_domain_component_infinity: Mapping[str, tuple[float, ...]]
    full_component_maximum_radii: Mapping[str, tuple[float, ...]]
    owned_component_maximum_radii: Mapping[str, tuple[float, ...]]
    owned_raw_global_norms: tuple[float, ...]
    component_finest_pair_orders: Mapping[str, float | None]
    component_status: Mapping[str, str]
    minimum_finite_finest_pair_order: float | None
    coarsest_guard_passed: bool
    finest_guard_passed: bool
    monotone_refinement_passed: bool
    finest_pair_order_passed: bool
    common_event_alignment_passed: bool
    admission_passed: bool


def owned_constraint_admission(
    samples: Sequence[OwnedConstraintSnapshot],
    *,
    method: str,
    coarsest_guard_maximum: Real,
    finest_guard_maximum: Real,
    minimum_finest_pair_order: Real = 1.5,
) -> OwnedConstraintAdmission:
    """Apply unchanged guards/order to the explicitly owned grid domain."""

    records = tuple(samples)
    if len(records) != 3 or any(
        not isinstance(item, OwnedConstraintSnapshot) for item in records
    ):
        raise ValueError("owned admission requires exactly three snapshots")
    counts = tuple(item.point_count for item in records)
    if tuple(sorted(counts)) != counts or len(set(counts)) != 3:
        raise ValueError("owned point counts must increase strictly")
    intervals = tuple(count - 1 for count in counts)
    ratios = tuple(
        right / left for left, right in zip(intervals[:-1], intervals[1:], strict=True)
    )
    if len(set(ratios)) != 1 or ratios[0] <= 1.0:
        raise ValueError("owned grids must share one refinement ratio")
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
    aligned = all(
        np.float64(value).tobytes() == np.float64(times[0]).tobytes()
        for value in times[1:]
    )

    full_values = {
        name: tuple(
            float(item.full_domain_norms.component_infinity[index])
            for item in records
        )
        for index, name in enumerate(PROTO4_CONSTRAINT_ORDER)
    }
    owned_values = {
        name: tuple(
            float(item.owned_domain_norms.component_infinity[index])
            for item in records
        )
        for index, name in enumerate(PROTO4_CONSTRAINT_ORDER)
    }
    full_radii = {
        name: tuple(float(item.full_component_maximum_radii[index]) for item in records)
        for index, name in enumerate(PROTO4_CONSTRAINT_ORDER)
    }
    owned_radii = {
        name: tuple(float(item.owned_component_maximum_radii[index]) for item in records)
        for index, name in enumerate(PROTO4_CONSTRAINT_ORDER)
    }

    component_orders: dict[str, float | None] = {}
    statuses: dict[str, str] = {}
    monotone = True
    order_pass = True
    finite_orders: list[float] = []
    for component, name in enumerate(PROTO4_CONSTRAINT_ORDER):
        values = tuple(
            float(item.effective_owned_norms_for_order_only.component_infinity[component])
            for item in records
        )
        if any(
            right > left for left, right in zip(values[:-1], values[1:], strict=True)
        ):
            monotone = False
        medium, fine = values[-2:]
        if medium == 0.0 and fine == 0.0:
            component_orders[name] = None
            statuses[name] = "accumulated_roundoff_enclosed_zero_pair"
        elif medium > 0.0 and fine == 0.0:
            component_orders[name] = None
            statuses[name] = "resolved_to_accumulated_roundoff_enclosure"
        elif medium == 0.0 and fine > 0.0:
            component_orders[name] = None
            statuses[name] = (
                "nonzero_reappeared_after_accumulated_roundoff_enclosure"
            )
            order_pass = False
        else:
            order = log(medium / fine) / log(ratios[-1])
            component_orders[name] = order
            statuses[name] = "finite_order"
            finite_orders.append(order)
            if order < order_minimum:
                order_pass = False
    owned_globals = tuple(item.owned_domain_norms.global_infinity for item in records)
    coarse_pass = owned_globals[0] < coarse_guard
    fine_pass = owned_globals[-1] < fine_guard
    passed = aligned and coarse_pass and fine_pass and monotone and order_pass
    return OwnedConstraintAdmission(
        method=method,
        point_counts=counts,
        coordinate_time=times[0],
        excluded_outer_rows=tuple(item.excluded_outer_rows for item in records),
        last_owned_radii=tuple(item.last_owned_radius for item in records),
        accepted_stage_counts=tuple(item.accepted_stage_count for item in records),
        accumulated_roundoff_enclosures=tuple(
            item.accumulated_roundoff_enclosure for item in records
        ),
        full_domain_component_infinity=full_values,
        owned_domain_component_infinity=owned_values,
        full_component_maximum_radii=full_radii,
        owned_component_maximum_radii=owned_radii,
        owned_raw_global_norms=owned_globals,
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


@dataclass(frozen=True, slots=True)
class FinestPairNestedSpectralAdmission:
    grid_point_counts: tuple[int, ...]
    individual_budget_passed_by_grid: Mapping[int, Mapping[str, bool]]
    coarsest_individual_budget_passed: bool
    finest_pair_individual_budgets_passed: bool
    every_nested_tail_ratio_passed: bool
    field_power_tail_ratios: Mapping[str, tuple[float | None, ...]]
    derivative_power_tail_ratios: Mapping[str, tuple[float | None, ...]]
    admission_passed: bool


def finest_pair_nested_spectral_admission(
    grid_point_counts: Sequence[int],
    budgets: Sequence[Mapping[str, SpectralPowerBudget]],
    *,
    thresholds: SpectralThresholds | None = None,
) -> FinestPairNestedSpectralAdmission:
    """Require unchanged absolute budgets on the finest pair plus all tails."""

    base: NestedSpectralAdmission = nested_spectral_admission(
        grid_point_counts,
        budgets,
        thresholds=thresholds,
    )
    counts = tuple(grid_point_counts)
    records = tuple(budgets)
    names = tuple(records[0])
    by_grid = {
        count: {
            name: bool(record[name].individual_admission_passed) for name in names
        }
        for count, record in zip(counts, records, strict=True)
    }
    coarse_pass = all(by_grid[counts[0]].values())
    fine_pair_pass = all(
        by_grid[count][name] for count in counts[-2:] for name in names
    )
    passed = fine_pair_pass and base.every_nested_tail_ratio_passed
    return FinestPairNestedSpectralAdmission(
        grid_point_counts=counts,
        individual_budget_passed_by_grid=by_grid,
        coarsest_individual_budget_passed=coarse_pass,
        finest_pair_individual_budgets_passed=fine_pair_pass,
        every_nested_tail_ratio_passed=base.every_nested_tail_ratio_passed,
        field_power_tail_ratios=base.field_power_tail_ratios,
        derivative_power_tail_ratios=base.derivative_power_tail_ratios,
        admission_passed=passed,
    )


__all__ = [
    "FinestPairNestedSpectralAdmission",
    "OwnedConstraintAdmission",
    "OwnedConstraintSnapshot",
    "finest_pair_nested_spectral_admission",
    "owned_constraint_admission",
    "owned_semidiscrete_constraint_snapshot",
]
