"""Pure PROTO13 method-owned calibration assessments.

PROTO13 keeps the PROTO12 equations, source, constraints, spectral gates, and
future-null orientation.  Its only runtime-level numerical change is that RK4
and SSPRK3 own different three-grid ladders.  This module therefore evaluates
the trapped-sphere sign on one exact physical-node set, gives each method its
own finest-pair Richardson enclosure, and requires those two enclosures to
overlap.  It performs no I/O and grants no run or physical authorization.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite
from numbers import Integral, Real
from typing import Any, Mapping, Sequence

import numpy as np

from .cal4_common_event_diagnosis import (
    OwnedConstraintAdmission,
    owned_constraint_admission,
)
from .gr0_calibration import radial_null_observables
from .health_monitor import compact_vacuum_buffer_window
from .numerical_engine import EvolutionState, UniformRadialGrid
from .proto11_runtime import reference_balanced_owned_constraint_snapshot
from .proto4_admission import (
    PROTO4_SPECTRAL_FIELD_ORDER,
    SpectralPowerBudget,
    SpectralThresholds,
    proper_radial_profile,
    spectral_field_budgets,
)
from .proto5_runtime import PROTO5_METHOD_CONTRACT
from .spectral_sensitivity import (
    SpectralTailSensitivity,
    ThreeGridProfileConvergence,
    spectral_tail_sensitivity,
    three_grid_profile_convergence,
)
from .spectral_sensitivity_v12 import (
    PairwiseResolvedOrSaturatedAdmission,
    pairwise_resolved_or_saturated_admission,
)


PROTO13_PRIMARY_POINT_COUNTS = (2049, 4097, 8193)
PROTO13_COMPARATOR_POINT_COUNTS = (4097, 8193, 16385)
PROTO13_COMMON_PHYSICAL_POINT_COUNT = 2049


def _finite(name: str, value: object, *, positive: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a finite real scalar")
    answer = float(value)
    if not isfinite(answer):
        raise ValueError(f"{name} must be finite")
    if positive and answer <= 0.0:
        raise ValueError(f"{name} must be positive")
    return answer


def _integer(name: str, value: object, *, minimum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral):
        raise TypeError(f"{name} must be an integer")
    answer = int(value)
    if answer < minimum:
        raise ValueError(f"{name} must be at least {minimum}")
    return answer


@dataclass(frozen=True, slots=True)
class Proto13MethodOwnedTrappedAssessment:
    """Two-method trapped-sign evidence on one exact physical-node set."""

    coordinate_time: float
    primary_point_counts: tuple[int, int, int]
    comparator_point_counts: tuple[int, int, int]
    common_physical_point_count: int
    primary_trapped_scores: tuple[float, float, float]
    comparator_trapped_scores: tuple[float, float, float]
    primary_Richardson_error: float
    comparator_Richardson_error: float
    primary_finest_interval: tuple[float, float]
    comparator_finest_interval: tuple[float, float]
    cross_method_interval_gap: float
    cross_method_interval_overlap: bool
    fine_cross_method_disagreement: float
    roundoff_error: float
    combined_error: float
    positive_margin_over_combined_error_factor: float | None
    both_method_owned_common_event_admissions_passed: bool
    trapped_sign_passed: bool


@dataclass(frozen=True, slots=True)
class Proto13GR0CommonEventAssessment:
    """Reference-balanced constraints plus the PROTO12 pairwise classifier."""

    method: str
    coordinate_time: float
    point_counts: tuple[int, int, int]
    accepted_stage_counts: tuple[int, int, int]
    constraint_admission: OwnedConstraintAdmission
    spectral_budgets: tuple[
        Mapping[str, SpectralPowerBudget],
        Mapping[str, SpectralPowerBudget],
        Mapping[str, SpectralPowerBudget],
    ]
    spectral_tail_sensitivities: tuple[
        Mapping[str, SpectralTailSensitivity],
        Mapping[str, SpectralTailSensitivity],
        Mapping[str, SpectralTailSensitivity],
    ]
    profile_convergence: ThreeGridProfileConvergence
    spatial_spectral_admission: PairwiseResolvedOrSaturatedAdmission
    maximum_spatial_round_trip_interpolation_infinity: float
    reference_state_map_applied: bool
    SRC4_source_backend_bound: bool
    interior_q_reprojected: bool
    admission_passed: bool


def proto13_gr0_common_event(
    states: Sequence[EvolutionState],
    grids: Sequence[UniformRadialGrid],
    *,
    accepted_stage_counts: Sequence[int],
    method: str,
    coordinate_time: Real,
    cutoff: Real = 16.0,
    measurement_radius_maximum: Real = 24.0,
    taper_fraction: Real = 1.0 / 8.0,
    fixed_outer_rows: int = 4,
    spectral_thresholds: SpectralThresholds | None = None,
) -> Proto13GR0CommonEventAssessment:
    """Evaluate one committed event on the method's exact owned ladder."""

    if method not in {"RK4", "SSPRK3"}:
        raise ValueError("unknown PROTO13 GR-0 method")
    expected = (
        PROTO13_PRIMARY_POINT_COUNTS
        if method == "RK4"
        else PROTO13_COMPARATOR_POINT_COUNTS
    )
    records, meshes = _validate_ladder(
        f"PROTO13 {method} common event", states, grids, expected
    )
    stages = tuple(accepted_stage_counts)
    if len(stages) != 3 or any(
        isinstance(value, bool) or not isinstance(value, Integral) or value < 0
        for value in stages
    ):
        raise ValueError("accepted_stage_counts must be three nonnegative integers")
    if fixed_outer_rows != 4:
        raise ValueError("PROTO13 fixed outer-row ownership differs")
    time = _finite("coordinate_time", coordinate_time)
    fraction = _finite("taper_fraction", taper_fraction, positive=True)
    if fraction >= 0.5:
        raise ValueError("taper_fraction must be below one half")
    contract = PROTO5_METHOD_CONTRACT[method]
    order = int(contract["spatial_order"])
    snapshots = tuple(
        reference_balanced_owned_constraint_snapshot(
            state,
            grid,
            coordinate_time=time,
            diagnostic_spatial_order=order,
            fixed_outer_rows=fixed_outer_rows,
            accepted_stage_count=int(stage),
        )
        for state, grid, stage in zip(records, meshes, stages, strict=True)
    )
    constraint = owned_constraint_admission(
        snapshots,
        method=method,
        coarsest_guard_maximum=contract["coarsest_constraint_guard"],
        finest_guard_maximum=contract["finest_constraint_guard"],
        minimum_finest_pair_order=1.5,
    )
    profiles = []
    budgets: list[dict[str, SpectralPowerBudget]] = []
    sensitivities: list[dict[str, SpectralTailSensitivity]] = []
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
        budget = dict(
            spectral_field_budgets(
                profile.deviations,
                proper,
                cutoff=cutoff,
                window=window,
                thresholds=spectral_thresholds,
            )
        )
        witness = {
            name: spectral_tail_sensitivity(
                profile.deviations[:, index],
                proper,
                window=window,
                round_trip_interpolation_infinity=(
                    profile.round_trip_interpolation_infinity[index]
                ),
                budget=budget[name],
            )
            for index, name in enumerate(PROTO4_SPECTRAL_FIELD_ORDER)
        }
        profiles.append(profile)
        budgets.append(budget)
        sensitivities.append(witness)
        interpolation_maximum = max(
            interpolation_maximum,
            float(np.max(profile.round_trip_interpolation_infinity, initial=0.0)),
        )
    convergence = three_grid_profile_convergence(expected, profiles)
    spatial = pairwise_resolved_or_saturated_admission(
        expected,
        budgets,
        sensitivities,
        convergence,
        thresholds=spectral_thresholds,
    )
    return Proto13GR0CommonEventAssessment(
        method=method,
        coordinate_time=time,
        point_counts=expected,
        accepted_stage_counts=(int(stages[0]), int(stages[1]), int(stages[2])),
        constraint_admission=constraint,
        spectral_budgets=(budgets[0], budgets[1], budgets[2]),
        spectral_tail_sensitivities=(
            sensitivities[0], sensitivities[1], sensitivities[2]
        ),
        profile_convergence=convergence,
        spatial_spectral_admission=spatial,
        maximum_spatial_round_trip_interpolation_infinity=interpolation_maximum,
        reference_state_map_applied=True,
        SRC4_source_backend_bound=True,
        interior_q_reprojected=False,
        admission_passed=(constraint.admission_passed and spatial.admission_passed),
    )


def _validate_ladder(
    name: str,
    states: Sequence[EvolutionState],
    grids: Sequence[UniformRadialGrid],
    expected: tuple[int, int, int],
) -> tuple[tuple[EvolutionState, ...], tuple[UniformRadialGrid, ...]]:
    state_tuple = tuple(states)
    grid_tuple = tuple(grids)
    if len(state_tuple) != 3 or len(grid_tuple) != 3:
        raise ValueError(f"{name} requires exactly three states and grids")
    counts = tuple(grid.point_count for grid in grid_tuple)
    if counts != expected:
        raise ValueError(f"{name} point counts differ from the frozen ladder")
    for state, grid in zip(state_tuple, grid_tuple, strict=True):
        if not isinstance(state, EvolutionState) or state.shape != (grid.point_count, 6):
            raise ValueError(f"{name} state/grid shape differs")
    return state_tuple, grid_tuple


def _shared_geometry(
    primary: Sequence[UniformRadialGrid],
    comparator: Sequence[UniformRadialGrid],
) -> tuple[float, float]:
    grids = (*primary, *comparator)
    lower = grids[0].minimum
    upper = grids[0].maximum
    for grid in grids:
        if (
            np.float64(grid.minimum).tobytes() != np.float64(lower).tobytes()
            or np.float64(grid.maximum).tobytes() != np.float64(upper).tobytes()
        ):
            raise ValueError("PROTO13 ladders do not share one exact radial domain")
        if (grid.point_count - 1) % (PROTO13_COMMON_PHYSICAL_POINT_COUNT - 1):
            raise ValueError("PROTO13 grid is not nested on the common physical mesh")
    return lower, upper


def _scores_on_common_physical_nodes(
    states: Sequence[EvolutionState],
    grids: Sequence[UniformRadialGrid],
    *,
    measurement_radius_maximum: float,
) -> tuple[float, float, float]:
    scores: list[float] = []
    common_intervals = PROTO13_COMMON_PHYSICAL_POINT_COUNT - 1
    for state, grid in zip(states, grids, strict=True):
        stride = (grid.point_count - 1) // common_intervals
        indices = np.arange(0, grid.point_count, stride, dtype=np.int64)
        if indices.size != PROTO13_COMMON_PHYSICAL_POINT_COUNT:
            raise RuntimeError("PROTO13 common-node restriction changed size")
        radii = grid.coordinates[indices]
        selected = indices[(radii > 0.0) & (radii <= measurement_radius_maximum)]
        if selected.size < 1:
            raise ValueError("PROTO13 trapped-score measurement interval is empty")
        observed = radial_null_observables(state)
        point_scores = np.minimum(
            -observed.theta_plus[selected],
            -observed.theta_minus[selected],
        )
        scores.append(float(np.max(point_scores, initial=-np.inf)))
    return tuple(scores)  # type: ignore[return-value]


def proto13_method_owned_trapped_assessment(
    *,
    primary_states: Sequence[EvolutionState],
    primary_grids: Sequence[UniformRadialGrid],
    comparator_states: Sequence[EvolutionState],
    comparator_grids: Sequence[UniformRadialGrid],
    primary_common_event: Any,
    comparator_common_event: Any,
    measurement_radius_maximum: Real = 24.0,
    minimum_observed_order: Real = 1.5,
    positive_margin_factor: Real = 4.0,
    roundoff_operation_budget: int = 4096,
) -> Proto13MethodOwnedTrappedAssessment:
    """Assess the fixed trapped sign without pretending unequal grids coincide.

    The scalar scores are evaluated on the same exact 2,049 physical nodes.
    Each method then supplies its own finest-pair Richardson interval.  The
    candidate sign can pass only when both common-event admissions pass, the
    two intervals overlap, and the smaller finest score exceeds the frozen
    multiple of the complete conservative error ledger.
    """

    p_states, p_grids = _validate_ladder(
        "PROTO13 primary ladder",
        primary_states,
        primary_grids,
        PROTO13_PRIMARY_POINT_COUNTS,
    )
    c_states, c_grids = _validate_ladder(
        "PROTO13 comparator ladder",
        comparator_states,
        comparator_grids,
        PROTO13_COMPARATOR_POINT_COUNTS,
    )
    _shared_geometry(p_grids, c_grids)
    if (
        getattr(primary_common_event, "method", None) != "RK4"
        or getattr(comparator_common_event, "method", None) != "SSPRK3"
        or tuple(getattr(primary_common_event, "point_counts", ()))
        != PROTO13_PRIMARY_POINT_COUNTS
        or tuple(getattr(comparator_common_event, "point_counts", ()))
        != PROTO13_COMPARATOR_POINT_COUNTS
    ):
        raise ValueError("PROTO13 common-event method or ladder ownership differs")
    time = _finite("coordinate_time", primary_common_event.coordinate_time)
    if np.float64(time).tobytes() != np.float64(
        comparator_common_event.coordinate_time
    ).tobytes():
        raise ValueError("PROTO13 common events are not bitwise time aligned")
    maximum = _finite(
        "measurement_radius_maximum", measurement_radius_maximum, positive=True
    )
    order = _finite("minimum_observed_order", minimum_observed_order, positive=True)
    factor = _finite("positive_margin_factor", positive_margin_factor, positive=True)
    budget = _integer(
        "roundoff_operation_budget", roundoff_operation_budget, minimum=1
    )
    primary_scores = _scores_on_common_physical_nodes(
        p_states, p_grids, measurement_radius_maximum=maximum
    )
    comparator_scores = _scores_on_common_physical_nodes(
        c_states, c_grids, measurement_radius_maximum=maximum
    )

    def richardson(scores: tuple[float, float, float], counts: tuple[int, int, int]) -> float:
        refinement = (counts[-1] - 1) / (counts[-2] - 1)
        if refinement <= 1.0:
            raise ValueError("PROTO13 refinement ratio must exceed one")
        return abs(scores[-1] - scores[-2]) / (refinement**order - 1.0)

    primary_error = richardson(primary_scores, PROTO13_PRIMARY_POINT_COUNTS)
    comparator_error = richardson(
        comparator_scores, PROTO13_COMPARATOR_POINT_COUNTS
    )
    primary_interval = (
        primary_scores[-1] - primary_error,
        primary_scores[-1] + primary_error,
    )
    comparator_interval = (
        comparator_scores[-1] - comparator_error,
        comparator_scores[-1] + comparator_error,
    )
    interval_gap = max(
        0.0,
        max(primary_interval[0], comparator_interval[0])
        - min(primary_interval[1], comparator_interval[1]),
    )
    interval_overlap = interval_gap == 0.0
    disagreement = abs(primary_scores[-1] - comparator_scores[-1])
    scale = max(
        1.0,
        *(abs(value) for value in (*primary_scores, *comparator_scores)),
    )
    roundoff = budget * np.finfo(np.float64).eps * scale
    combined = primary_error + comparator_error + interval_gap + roundoff
    fine_margin = min(primary_scores[-1], comparator_scores[-1])
    admissions = bool(
        primary_common_event.admission_passed
        and comparator_common_event.admission_passed
    )
    passed = (
        admissions
        and interval_overlap
        and fine_margin > factor * combined
    )
    return Proto13MethodOwnedTrappedAssessment(
        coordinate_time=time,
        primary_point_counts=PROTO13_PRIMARY_POINT_COUNTS,
        comparator_point_counts=PROTO13_COMPARATOR_POINT_COUNTS,
        common_physical_point_count=PROTO13_COMMON_PHYSICAL_POINT_COUNT,
        primary_trapped_scores=primary_scores,
        comparator_trapped_scores=comparator_scores,
        primary_Richardson_error=primary_error,
        comparator_Richardson_error=comparator_error,
        primary_finest_interval=primary_interval,
        comparator_finest_interval=comparator_interval,
        cross_method_interval_gap=interval_gap,
        cross_method_interval_overlap=interval_overlap,
        fine_cross_method_disagreement=disagreement,
        roundoff_error=roundoff,
        combined_error=combined,
        positive_margin_over_combined_error_factor=(
            fine_margin / combined if combined > 0.0 else None
        ),
        both_method_owned_common_event_admissions_passed=admissions,
        trapped_sign_passed=passed,
    )


__all__ = [
    "PROTO13_COMMON_PHYSICAL_POINT_COUNT",
    "PROTO13_COMPARATOR_POINT_COUNTS",
    "PROTO13_PRIMARY_POINT_COUNTS",
    "Proto13GR0CommonEventAssessment",
    "Proto13MethodOwnedTrappedAssessment",
    "proto13_gr0_common_event",
    "proto13_method_owned_trapped_assessment",
]
