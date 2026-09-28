"""Causal-domain bookkeeping for the scoped spherical FGC experiment.

HYP2 bounds characteristic roots in a physical orthonormal frame.  This
module performs the missing ADM conversion and accumulates a conservative
radial coordinate-distance budget.  It deliberately does *not* claim that the
outer boundary conditions form a nonlinear constraint-preserving IBVP.  The
contract is instead that no accepted measurement may lie in the continuum
domain of dependence of that boundary, after padding by the declared SBP
stencil reach.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import isfinite
from numbers import Real
from typing import Iterable, Sequence


Q = Fraction

# HYP2's real Riesz clusters in a physical orthonormal frame.  The radius is
# included because a separated contour localizes, but does not pin, each root.
LOCAL_CONE_SPEED_BOUNDS = {
    "physical": Q(6, 5),       # |+/-1| + 1/5
    "tilde": Q(11, 20),        # |+/-1/2| + 1/20
    "hat": Q(23, 60),          # |+/-1/3| + 1/20
}
ALL_CONE_LOCAL_SPEED_BOUND = max(LOCAL_CONE_SPEED_BOUNDS.values())


def _finite(name: str, value: Real) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a real number")
    answer = float(value)
    if not isfinite(answer):
        raise ValueError(f"{name} must be finite")
    return answer


def coordinate_speed_upper_bound(
    *,
    lapse: Real,
    radial_metric: Real,
    shift: Real,
    local_speed_bound: Real = ALL_CONE_LOCAL_SPEED_BOUND,
) -> float:
    """Return ``|shift| + lapse/local_metric * |z|`` for one ADM state.

    With ``ds^2=-N^2 dt^2+Lambda^2(dr+s dt)^2+R^2 dOmega^2``, a
    characteristic whose orthonormal radial speed is ``z`` has coordinate
    speed ``dr/dt=-s+(N/Lambda)z``.  The returned triangle bound covers both
    orientations and all HYP2 clusters when the default is used.
    """

    n = _finite("lapse", lapse)
    lam = _finite("radial_metric", radial_metric)
    beta = _finite("shift", shift)
    z = _finite("local_speed_bound", local_speed_bound)
    if n <= 0.0 or lam <= 0.0 or z < 0.0:
        raise ValueError("lapse and radial_metric must be positive and speed nonnegative")
    return abs(beta) + n * z / lam


def exact_coordinate_speed_upper_bound(
    *,
    lapse: Fraction,
    radial_metric: Fraction,
    shift: Fraction,
    local_speed_bound: Fraction = ALL_CONE_LOCAL_SPEED_BOUND,
) -> Fraction:
    """Exact rational counterpart used by the canonical certificate."""

    values = (lapse, radial_metric, shift, local_speed_bound)
    if any(not isinstance(value, Fraction) for value in values):
        raise TypeError("exact coordinate-speed inputs must be Fraction values")
    if lapse <= 0 or radial_metric <= 0 or local_speed_bound < 0:
        raise ValueError("lapse and radial_metric must be positive and speed nonnegative")
    return abs(shift) + lapse * local_speed_bound / radial_metric


@dataclass(frozen=True, slots=True)
class BoundaryGeometry:
    """Fixed extraction/boundary geometry for one numerical grid."""

    outer_radius: float
    measurement_radius: float
    minimum_causal_buffer: float
    sbp_stencil_reach: float

    def __post_init__(self) -> None:
        values = {
            "outer_radius": self.outer_radius,
            "measurement_radius": self.measurement_radius,
            "minimum_causal_buffer": self.minimum_causal_buffer,
            "sbp_stencil_reach": self.sbp_stencil_reach,
        }
        normalized = {name: _finite(name, value) for name, value in values.items()}
        if normalized["measurement_radius"] < 0.0:
            raise ValueError("measurement_radius must be nonnegative")
        if normalized["minimum_causal_buffer"] <= 0.0:
            raise ValueError("minimum_causal_buffer must be positive")
        if normalized["sbp_stencil_reach"] < 0.0:
            raise ValueError("sbp_stencil_reach must be nonnegative")
        if normalized["outer_radius"] <= normalized["measurement_radius"]:
            raise ValueError("outer_radius must exceed measurement_radius")
        for name, value in normalized.items():
            object.__setattr__(self, name, value)


@dataclass(frozen=True, slots=True)
class CausalBudgetState:
    """Immutable causal ledger at the last accepted Runge--Kutta stage."""

    accepted_time: float = 0.0
    accumulated_characteristic_distance: float = 0.0
    previous_speed_upper: float = 0.0

    def __post_init__(self) -> None:
        for name in (
            "accepted_time",
            "accumulated_characteristic_distance",
            "previous_speed_upper",
        ):
            value = _finite(name, getattr(self, name))
            if value < 0.0:
                raise ValueError(f"{name} must be nonnegative")
            object.__setattr__(self, name, value)


@dataclass(frozen=True, slots=True)
class CausalStepAssessment:
    """Pure preview of a trial stage interval; the caller chooses acceptance."""

    candidate_state: CausalBudgetState
    interval_speed_upper: float
    remaining_causal_buffer: float
    strict_margin_over_required_buffer: float
    admissible: bool


class BoundaryControlStop(RuntimeError):
    """Typed refusal raised before a boundary-contaminated stage is accepted."""

    def __init__(self, assessment: CausalStepAssessment) -> None:
        self.reason = "boundary_causal_buffer"
        self.assessment = assessment
        super().__init__(
            "boundary causal buffer would be exhausted before accepting trial stage"
        )


def assess_causal_step(
    state: CausalBudgetState,
    geometry: BoundaryGeometry,
    *,
    trial_time: Real,
    stage_coordinate_speed_uppers: Iterable[Real],
    candidate_endpoint_coordinate_speed_upper: Real,
) -> CausalStepAssessment:
    """Conservatively debit one proposed accepted-stage interval.

    The maximum includes the preceding accepted endpoint, every supplied
    internal trial-stage envelope, and the speed evaluated on the proposed
    final Runge--Kutta combination.  The endpoint is explicit because the last
    internal stage is not generally the accepted state.  This is the runtime
    contract NUM1 must call before committing a step.  A higher-fidelity
    integrator may supply additional dense-output samples; omitting the actual
    Runge--Kutta stages or final candidate endpoint is invalid.
    """

    if not isinstance(state, CausalBudgetState):
        raise TypeError("state must be CausalBudgetState")
    if not isinstance(geometry, BoundaryGeometry):
        raise TypeError("geometry must be BoundaryGeometry")
    end = _finite("trial_time", trial_time)
    if end <= state.accepted_time:
        raise ValueError("trial_time must be strictly later than accepted_time")
    samples = tuple(
        _finite(f"stage_coordinate_speed_uppers[{index}]", value)
        for index, value in enumerate(stage_coordinate_speed_uppers)
    )
    if not samples:
        raise ValueError("at least one trial-stage speed envelope is required")
    if any(value < 0.0 for value in samples):
        raise ValueError("stage coordinate-speed envelopes must be nonnegative")
    endpoint_speed = _finite(
        "candidate_endpoint_coordinate_speed_upper",
        candidate_endpoint_coordinate_speed_upper,
    )
    if endpoint_speed < 0.0:
        raise ValueError("candidate endpoint coordinate-speed envelope must be nonnegative")
    interval_speed = max((state.previous_speed_upper, *samples, endpoint_speed))
    distance = state.accumulated_characteristic_distance + interval_speed * (
        end - state.accepted_time
    )
    remaining = (
        geometry.outer_radius
        - geometry.measurement_radius
        - geometry.sbp_stencil_reach
        - distance
    )
    margin = remaining - geometry.minimum_causal_buffer
    candidate = CausalBudgetState(end, distance, endpoint_speed)
    return CausalStepAssessment(
        candidate_state=candidate,
        interval_speed_upper=interval_speed,
        remaining_causal_buffer=remaining,
        strict_margin_over_required_buffer=margin,
        admissible=margin > 0.0,
    )


def accept_causal_step(
    state: CausalBudgetState,
    geometry: BoundaryGeometry,
    *,
    trial_time: Real,
    stage_coordinate_speed_uppers: Iterable[Real],
    candidate_endpoint_coordinate_speed_upper: Real,
) -> CausalBudgetState:
    """Return a new immutable state or raise without altering ``state``."""

    assessment = assess_causal_step(
        state,
        geometry,
        trial_time=trial_time,
        stage_coordinate_speed_uppers=stage_coordinate_speed_uppers,
        candidate_endpoint_coordinate_speed_upper=candidate_endpoint_coordinate_speed_upper,
    )
    if not assessment.admissible:
        raise BoundaryControlStop(assessment)
    return assessment.candidate_state


@dataclass(frozen=True, slots=True)
class AlignedGrid:
    outer_radius: Fraction
    interval_count: int
    point_count: int
    spacing: Fraction
    measurement_index: int
    stencil_reach: Fraction


def aligned_outer_boundary_grids(
    *,
    nominal_outer_radius: Fraction,
    nominal_point_counts: Sequence[int],
    outer_radii: Sequence[Fraction],
    measurement_radius: Fraction,
    stencil_reach_intervals: int,
) -> tuple[AlignedGrid, ...]:
    """Move the boundary while preserving every nominal interior spacing."""

    if not isinstance(nominal_outer_radius, Fraction) or nominal_outer_radius <= 0:
        raise ValueError("nominal_outer_radius must be a positive Fraction")
    if not isinstance(measurement_radius, Fraction) or measurement_radius < 0:
        raise ValueError("measurement_radius must be a nonnegative Fraction")
    if isinstance(stencil_reach_intervals, bool) or stencil_reach_intervals < 0:
        raise ValueError("stencil_reach_intervals must be a nonnegative integer")
    if not nominal_point_counts or not outer_radii:
        raise ValueError("point counts and outer radii must be nonempty")
    answer: list[AlignedGrid] = []
    for point_count in nominal_point_counts:
        if isinstance(point_count, bool) or not isinstance(point_count, int) or point_count < 2:
            raise ValueError("nominal point counts must be integers at least two")
        spacing = nominal_outer_radius / (point_count - 1)
        measurement_intervals = measurement_radius / spacing
        if measurement_intervals.denominator != 1:
            raise ValueError("measurement radius is not aligned to the nominal grid")
        for outer in outer_radii:
            if not isinstance(outer, Fraction) or outer <= measurement_radius:
                raise ValueError("each outer radius must be a Fraction beyond measurement")
            intervals = outer / spacing
            if intervals.denominator != 1:
                raise ValueError("moved outer radius is not aligned at fixed spacing")
            answer.append(
                AlignedGrid(
                    outer_radius=outer,
                    interval_count=intervals.numerator,
                    point_count=intervals.numerator + 1,
                    spacing=spacing,
                    measurement_index=measurement_intervals.numerator,
                    stencil_reach=stencil_reach_intervals * spacing,
                )
            )
    return tuple(answer)


def exact_reference_budget(
    grid: AlignedGrid,
    *,
    final_time: Fraction,
    minimum_causal_buffer: Fraction,
    coordinate_speed_bound: Fraction = ALL_CONE_LOCAL_SPEED_BOUND,
) -> dict[str, Fraction | bool]:
    """Exact flat-reference margin for one aligned grid."""

    values = (final_time, minimum_causal_buffer, coordinate_speed_bound)
    if any(not isinstance(value, Fraction) for value in values):
        raise TypeError("reference-budget scalars must be Fractions")
    if final_time <= 0 or minimum_causal_buffer <= 0 or coordinate_speed_bound < 0:
        raise ValueError("reference-budget time/buffer must be positive and speed nonnegative")
    remaining = (
        grid.outer_radius
        - grid.measurement_index * grid.spacing
        - grid.stencil_reach
        - coordinate_speed_bound * final_time
    )
    margin = remaining - minimum_causal_buffer
    return {
        "characteristic_travel": coordinate_speed_bound * final_time,
        "remaining_causal_buffer": remaining,
        "strict_margin_over_required_buffer": margin,
        "passed": margin > 0,
    }
