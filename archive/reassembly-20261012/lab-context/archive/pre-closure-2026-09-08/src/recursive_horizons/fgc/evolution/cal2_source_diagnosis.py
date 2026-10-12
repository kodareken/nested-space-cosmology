"""Outcome-neutral diagnostics for the failed PROTO5 GR-0 source gate.

The PROTO5 campaign stopped when the direct affine GR-0 source solver missed
an absolute residual tolerance on an *unaccepted* Runge--Kutta stage.  This
module deliberately does not alter that solver.  It reconstructs the same
affine system and returns the residual/refinement evidence even when the raw
gate is missed, allowing CAL2-PREF4 to distinguish an accepted-state
continuum obstruction from a rejectable trial-step failure.

Nothing here authorizes a successor run.  In particular, the diagnostic
operator is not an evolution backend: it exists only to replay proposed
stages from an immutable checkpoint under declared step-size factors.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite
from numbers import Real

import numpy as np

from .gr0_calibration import (
    ADM_CENTER_PARITIES,
    FIELD_COUNT,
    Q_CENTER_PARITIES,
)
from .gr0_direct_source import gr0_ref1_residual_batch
from .numerical_engine import EvolutionRHS, EvolutionState, SBPFirstDerivative
from .vectorized_source import regular_center_acceleration_limit


def _finite(name: str, value: object, *, positive: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a finite real scalar")
    answer = float(value)
    if not isfinite(answer):
        raise ValueError(f"{name} must be finite")
    if positive and answer <= 0.0:
        raise ValueError(f"{name} must be positive")
    return answer


def _array(name: str, value: object, *, ndim: int) -> np.ndarray:
    answer = np.asarray(value, dtype=np.float64)
    if answer.ndim != ndim or not np.all(np.isfinite(answer)):
        raise ValueError(f"{name} must be a finite rank-{ndim} binary64 array")
    return np.ascontiguousarray(answer)


@dataclass(frozen=True, slots=True)
class ResidualIteration:
    iteration: int
    residual_infinity: float
    acceleration_correction_infinity: float
    relative_acceleration_correction: float
    maximum_residual_point_index: int
    maximum_residual_row_index: int
    maximum_residual_radius: float

    def __post_init__(self) -> None:
        if isinstance(self.iteration, bool) or not isinstance(self.iteration, int):
            raise TypeError("iteration must be int")
        if self.iteration < 0:
            raise ValueError("iteration must be nonnegative")
        for name in (
            "residual_infinity",
            "acceleration_correction_infinity",
            "relative_acceleration_correction",
            "maximum_residual_radius",
        ):
            value = _finite(name, getattr(self, name))
            if value < 0.0:
                raise ValueError(f"{name} must be nonnegative")
            object.__setattr__(self, name, value)
        for name, upper in (
            ("maximum_residual_point_index", None),
            ("maximum_residual_row_index", FIELD_COUNT),
        ):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value < 0:
                raise ValueError(f"{name} must be a nonnegative int")
            if upper is not None and value >= upper:
                raise ValueError(f"{name} lies outside the residual rows")


@dataclass(frozen=True, slots=True)
class AffineRootDiagnosis:
    accelerations: np.ndarray
    residuals: np.ndarray
    residual_infinity: float
    raw_tolerance: float
    raw_gate_passed: bool
    kinetic_condition_infinity_maximum: float
    residual_decreased_until_floor: bool
    roundoff_floor_reached: bool
    iterations: tuple[ResidualIteration, ...]

    def __post_init__(self) -> None:
        acceleration = _array("accelerations", self.accelerations, ndim=2).copy()
        residual = _array("residuals", self.residuals, ndim=2).copy()
        if acceleration.shape != residual.shape or acceleration.shape[1] != FIELD_COUNT:
            raise ValueError("diagnostic acceleration and residual shapes differ")
        acceleration.setflags(write=False)
        residual.setflags(write=False)
        object.__setattr__(self, "accelerations", acceleration)
        object.__setattr__(self, "residuals", residual)
        for name in (
            "residual_infinity",
            "raw_tolerance",
            "kinetic_condition_infinity_maximum",
        ):
            value = _finite(name, getattr(self, name), positive=name == "raw_tolerance")
            if value < 0.0:
                raise ValueError(f"{name} must be nonnegative")
            object.__setattr__(self, name, value)
        if not isinstance(self.raw_gate_passed, bool):
            raise TypeError("raw_gate_passed must be bool")
        if not isinstance(self.residual_decreased_until_floor, bool):
            raise TypeError("residual_decreased_until_floor must be bool")
        if not isinstance(self.roundoff_floor_reached, bool):
            raise TypeError("roundoff_floor_reached must be bool")
        if not self.iterations:
            raise ValueError("at least one residual iteration is required")


def diagnose_gr0_grid_accelerations(
    u: object,
    p: object,
    q: object,
    p_r: object,
    q_r: object,
    radii: object,
    *,
    raw_tolerance: Real = 1.0e-12,
    condition_number_maximum: Real = 1.0e10,
    maximum_refinement_iterations: int = 16,
) -> AffineRootDiagnosis:
    """Return the affine-root evidence without promoting a missed raw gate."""

    base = _array("u", u, ndim=2)
    velocity = _array("p", p, ndim=2)
    gradient = _array("q", q, ndim=2)
    mixed = _array("p_r", p_r, ndim=2)
    radial_second = _array("q_r", q_r, ndim=2)
    radius = _array("radii", radii, ndim=1)
    points, fields = base.shape
    if fields != FIELD_COUNT:
        raise ValueError("u must contain the six canonical ADM fields")
    if (
        velocity.shape != base.shape
        or gradient.shape != base.shape
        or mixed.shape != base.shape
        or radial_second.shape != base.shape
        or radius.shape != (points,)
    ):
        raise ValueError("diagnostic GR-0 source shapes differ")
    if np.any(radius <= 0.0):
        raise ValueError("diagnostic radii must be positive")
    tolerance = _finite("raw_tolerance", raw_tolerance, positive=True)
    condition_limit = _finite(
        "condition_number_maximum", condition_number_maximum, positive=True
    )
    if (
        isinstance(maximum_refinement_iterations, bool)
        or not isinstance(maximum_refinement_iterations, int)
        or not 0 <= maximum_refinement_iterations <= 16
    ):
        raise ValueError("maximum_refinement_iterations must lie in [0,16]")

    seeds = np.zeros((FIELD_COUNT + 1, points, FIELD_COUNT), dtype=np.float64)
    for field in range(FIELD_COUNT):
        seeds[field + 1, :, field] = 1.0
    evaluated = gr0_ref1_residual_batch(
        base, velocity, gradient, seeds, mixed, radial_second, radius
    )
    constant = evaluated.full_residual[0]
    jacobian = np.empty((points, FIELD_COUNT, FIELD_COUNT), dtype=np.float64)
    for field in range(FIELD_COUNT):
        jacobian[:, :, field] = evaluated.full_residual[field + 1] - constant
    try:
        inverse = np.linalg.inv(jacobian)
        acceleration = np.linalg.solve(jacobian, -constant[..., None])[..., 0]
    except np.linalg.LinAlgError as error:
        raise ValueError("diagnostic GR-0 kinetic block is singular") from error
    condition = np.max(np.sum(np.abs(jacobian), axis=2), axis=1) * np.max(
        np.sum(np.abs(inverse), axis=2), axis=1
    )
    condition_maximum = float(np.max(condition, initial=0.0))
    if not np.all(np.isfinite(condition)) or condition_maximum >= condition_limit:
        raise ValueError("diagnostic GR-0 kinetic condition limit reached")

    iterations: list[ResidualIteration] = []
    residual: np.ndarray | None = None
    residual_norm: float | None = None
    prior_norm: float | None = None
    decreased = True
    floor = False
    for iteration in range(maximum_refinement_iterations + 1):
        verified = gr0_ref1_residual_batch(
            base,
            velocity,
            gradient,
            acceleration[None, :, :],
            mixed,
            radial_second,
            radius,
        )
        residual = verified.full_residual[0]
        residual_norm = float(np.max(np.abs(residual), initial=0.0))
        maximum_index = np.unravel_index(
            int(np.argmax(np.abs(residual))), residual.shape
        )
        correction = np.linalg.solve(jacobian, -residual[..., None])[..., 0]
        correction_norm = float(np.max(np.abs(correction), initial=0.0))
        acceleration_scale = max(
            1.0, float(np.max(np.abs(acceleration), initial=0.0))
        )
        iterations.append(
            ResidualIteration(
                iteration=iteration,
                residual_infinity=residual_norm,
                acceleration_correction_infinity=correction_norm,
                relative_acceleration_correction=correction_norm
                / acceleration_scale,
                maximum_residual_point_index=int(maximum_index[0]),
                maximum_residual_row_index=int(maximum_index[1]),
                maximum_residual_radius=float(radius[maximum_index[0]]),
            )
        )
        if residual_norm <= tolerance or iteration == maximum_refinement_iterations:
            break
        candidate = acceleration + correction
        candidate_residual = gr0_ref1_residual_batch(
            base,
            velocity,
            gradient,
            candidate[None, :, :],
            mixed,
            radial_second,
            radius,
        ).full_residual[0]
        candidate_norm = float(
            np.max(np.abs(candidate_residual), initial=0.0)
        )
        if prior_norm is not None and residual_norm >= prior_norm:
            decreased = False
        if candidate_norm >= residual_norm:
            floor = True
            break
        prior_norm = residual_norm
        acceleration = candidate

    if residual is None or residual_norm is None:
        raise RuntimeError("diagnostic root produced no residual")
    return AffineRootDiagnosis(
        accelerations=acceleration,
        residuals=residual,
        residual_infinity=residual_norm,
        raw_tolerance=tolerance,
        raw_gate_passed=residual_norm <= tolerance,
        kinetic_condition_infinity_maximum=condition_maximum,
        residual_decreased_until_floor=decreased,
        roundoff_floor_reached=floor,
        iterations=tuple(iterations),
    )


class DiagnosticGR0EvolutionOperator:
    """Replay operator that exposes failed source roots without accepting them."""

    def __init__(
        self,
        grid: object,
        *,
        spatial_order: int,
        ko_dissipation: Real,
        raw_tolerance: Real,
        kinetic_condition_maximum: Real,
    ) -> None:
        self.grid = grid
        self.derivative = SBPFirstDerivative(grid, spatial_order)
        self.ko_dissipation = _finite("ko_dissipation", ko_dissipation)
        if self.ko_dissipation < 0.0:
            raise ValueError("ko_dissipation must be nonnegative")
        self.raw_tolerance = _finite("raw_tolerance", raw_tolerance, positive=True)
        self.kinetic_condition_maximum = _finite(
            "kinetic_condition_maximum", kinetic_condition_maximum, positive=True
        )

    def __call__(self, _time: float, state: EvolutionState) -> EvolutionRHS:
        if not isinstance(state, EvolutionState) or state.shape != (
            self.grid.point_count,
            FIELD_COUNT,
        ):
            raise ValueError("diagnostic GR-0 RHS state shape differs")
        p_r = self.derivative.differentiate(
            state.p, center_parities=ADM_CENTER_PARITIES
        )
        q_r = self.derivative.differentiate(
            state.q, center_parities=Q_CENTER_PARITIES
        )
        source = diagnose_gr0_grid_accelerations(
            state.u[1:],
            state.p[1:],
            state.q[1:],
            p_r[1:],
            q_r[1:],
            self.grid.coordinates[1:],
            raw_tolerance=self.raw_tolerance,
            condition_number_maximum=self.kinetic_condition_maximum,
        )
        acceleration = np.empty_like(state.p)
        acceleration[1:] = source.accelerations
        center = regular_center_acceleration_limit(acceleration[1:5])
        acceleration[0] = center.acceleration
        dissipation = self.derivative.kreiss_oliger_dissipation(
            state.p,
            coefficient=self.ko_dissipation,
            center_parities=ADM_CENTER_PARITIES,
        )
        reduction = state.q - self.derivative.differentiate(
            state.u, center_parities=ADM_CENTER_PARITIES
        )
        coordinate_speed = np.abs(state.u[:, 1]) + 1.2 * state.u[:, 0] / state.u[:, 2]
        last = source.iterations[-1]
        return EvolutionRHS(
            state.p,
            acceleration + dissipation,
            p_r,
            {
                "source_residual_infinity": source.residual_infinity,
                "source_raw_gate_passed": source.raw_gate_passed,
                "source_refinement_iterations": len(source.iterations) - 1,
                "source_residual_decreased_monotonically": (
                    source.residual_decreased_until_floor
                ),
                "source_roundoff_floor_reached": source.roundoff_floor_reached,
                "source_relative_acceleration_correction": (
                    last.relative_acceleration_correction
                ),
                "source_maximum_residual_point_index": (
                    last.maximum_residual_point_index
                ),
                "source_maximum_residual_row_index": last.maximum_residual_row_index,
                "source_maximum_residual_radius": last.maximum_residual_radius,
                "kinetic_condition_infinity": (
                    source.kinetic_condition_infinity_maximum
                ),
                "acceleration_infinity": float(
                    np.max(np.abs(acceleration), initial=0.0)
                ),
                "center_acceleration_estimator_infinity": center.estimator_infinity,
                "reduction_constraint_infinity": float(
                    np.max(np.abs(reduction), initial=0.0)
                ),
                "coordinate_speed_upper": float(
                    np.max(coordinate_speed, initial=0.0)
                ),
                "minimum_lapse": float(np.min(state.u[:, 0])),
                "minimum_radial_metric": float(np.min(state.u[:, 2])),
                "minimum_areal_radius_away_from_center": float(
                    np.min(state.u[1:, 3])
                ),
            },
        )


__all__ = [
    "AffineRootDiagnosis",
    "DiagnosticGR0EvolutionOperator",
    "ResidualIteration",
    "diagnose_gr0_grid_accelerations",
]
