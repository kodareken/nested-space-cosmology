"""Constraint-solved GR-0 grids and the frozen REF1 calibration operator.

The preserved PROTO3 contract would have selected the future FGC-QR matter
amplitude from GR-0 outcomes only.  CAL0 closes that executable contract, but
its outcome-neutral amplitude order remains the inherited premise for a valid
successor.  This module supplies the GR-0 calibration substrate without
opening an FGC-QR trajectory.  It embeds the already certified maximal
polar--areal GR-0 constraint solution into a complete centre-regular grid,
then evolves the unchanged GR-0 specialization of the REF1 equations through
the independently cross-checked direct-array source.

The outer rows are held at their analytic initial vacuum data.  This is not
promoted to a nonlinear constraint-preserving boundary condition: BND2's
causal ledger and the later outer-boundary move comparison must exclude those
rows from every retained conclusion.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite
from numbers import Real
from typing import Callable

import numpy as np

from ..initial_data_family import gauge_compatible_metric_time_derivatives
from ..initial_data_preflight import (
    PulseParameters,
    gr0_constraint_rhs,
    pulse_fields,
    solve_gr0_initial_slice,
)
from .boundary_domain import ALL_CONE_LOCAL_SPEED_BOUND
from .gr0_direct_source import solve_gr0_grid_accelerations
from .numerical_engine import (
    EvolutionRHS,
    EvolutionState,
    SBPFirstDerivative,
    UniformRadialGrid,
)
from .vectorized_source import ADM_CENTER_PARITIES, regular_center_acceleration_limit


FIELD_COUNT = ADM_CENTER_PARITIES.size
Q_CENTER_PARITIES = -ADM_CENTER_PARITIES


def _finite(name: str, value: object, *, positive: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a finite real scalar")
    answer = float(value)
    if not isfinite(answer):
        raise ValueError(f"{name} must be finite")
    if positive and answer <= 0.0:
        raise ValueError(f"{name} must be positive")
    return answer


def _freeze_array(name: str, value: object, *, ndim: int) -> np.ndarray:
    answer = np.asarray(value, dtype=np.float64)
    if answer.ndim != ndim or not np.all(np.isfinite(answer)):
        raise ValueError(f"{name} must be a finite rank-{ndim} binary64 array")
    answer = np.ascontiguousarray(answer).copy()
    answer.setflags(write=False)
    return answer


@dataclass(frozen=True, slots=True)
class GR0GridInitialData:
    """One complete regular-centre GR-0 calibration slice."""

    grid: UniformRadialGrid
    parameters: PulseParameters
    constraint_method: str
    state: EvolutionState
    support_minimum_index: int
    support_maximum_index: int
    peak_compactness: float
    peak_radius: float
    outer_mass: float
    vacuum_momentum_constant: float
    reduction_constraint_infinity: float
    no_initial_trapped_sphere: bool

    def __post_init__(self) -> None:
        if not isinstance(self.grid, UniformRadialGrid):
            raise TypeError("grid must be UniformRadialGrid")
        if not isinstance(self.parameters, PulseParameters):
            raise TypeError("parameters must be PulseParameters")
        if not isinstance(self.state, EvolutionState) or self.state.shape != (
            self.grid.point_count,
            FIELD_COUNT,
        ):
            raise ValueError("initial state does not match its six-field grid")
        if self.constraint_method not in {"RK4", "SSPRK3"}:
            raise ValueError("unknown initial constraint method")
        if not (
            0 < self.support_minimum_index < self.support_maximum_index
            < self.grid.point_count - 1
        ):
            raise ValueError("support indices do not define an interior shell")
        for name in (
            "peak_compactness",
            "peak_radius",
            "outer_mass",
            "vacuum_momentum_constant",
            "reduction_constraint_infinity",
        ):
            value = _finite(name, getattr(self, name))
            if name in {"peak_radius", "outer_mass"} and value <= 0.0:
                raise ValueError(f"{name} must be positive")
            if name == "reduction_constraint_infinity" and value < 0.0:
                raise ValueError("reduction constraint norm must be nonnegative")
            object.__setattr__(self, name, value)
        if not isinstance(self.no_initial_trapped_sphere, bool):
            raise TypeError("no_initial_trapped_sphere must be bool")


@dataclass(frozen=True, slots=True)
class GR0SliceObservables:
    """Metric-null expansions and Misner--Sharp diagnostics on one slice."""

    theta_plus: np.ndarray
    theta_minus: np.ndarray
    compactness: np.ndarray
    misner_sharp_mass: np.ndarray
    normal_areal_derivative: np.ndarray

    def __post_init__(self) -> None:
        arrays = []
        for name in (
            "theta_plus",
            "theta_minus",
            "compactness",
            "misner_sharp_mass",
            "normal_areal_derivative",
        ):
            value = _freeze_array(name, getattr(self, name), ndim=1)
            arrays.append(value)
            object.__setattr__(self, name, value)
        if len({value.shape for value in arrays}) != 1:
            raise ValueError("slice observable shapes differ")


def radial_null_observables(state: EvolutionState) -> GR0SliceObservables:
    """Return fixed-orientation physical radial-null expansion data.

    For ``ds^2=-alpha^2 dt^2+lambda^2(dr+shift dt)^2+R^2dOmega^2``,
    ``N(R)=(R_t-shift*R_r)/alpha`` and the two future radial directions
    have areal derivatives ``N(R)+/-R_r/lambda``.  The common positive
    normalization is irrelevant to the trapped-sphere sign test.
    """

    if not isinstance(state, EvolutionState) or state.shape[1] != FIELD_COUNT:
        raise TypeError("state must contain the six canonical ADM fields")
    alpha = state.u[:, 0]
    shift = state.u[:, 1]
    radial_metric = state.u[:, 2]
    areal = state.u[:, 3]
    if np.any(alpha <= 0.0) or np.any(radial_metric <= 0.0):
        raise ValueError("null observables require positive lapse and radial metric")
    normal = (state.p[:, 3] - shift * state.q[:, 3]) / alpha
    spatial = state.q[:, 3] / radial_metric
    plus = np.zeros_like(areal)
    minus = np.zeros_like(areal)
    positive = areal > 0.0
    plus[positive] = 2.0 * (normal[positive] + spatial[positive]) / areal[positive]
    minus[positive] = 2.0 * (normal[positive] - spatial[positive]) / areal[positive]
    compactness = 1.0 + normal**2 - spatial**2
    compactness[~positive] = 0.0
    mass = areal * compactness / 2.0
    return GR0SliceObservables(plus, minus, compactness, mass, normal)


def construct_gr0_grid_initial_data(
    parameters: PulseParameters,
    *,
    point_count: int,
    outer_radius: Real = 128.0,
    constraint_method: str = "RK4",
    diagnostic_spatial_order: int = 4,
) -> GR0GridInitialData:
    """Embed the GR-0 constraint solution and its exact vacuum continuations."""

    if not isinstance(parameters, PulseParameters):
        raise TypeError("parameters must be PulseParameters")
    outer = _finite("outer_radius", outer_radius, positive=True)
    if outer <= parameters.support_maximum:
        raise ValueError("outer radius must retain a vacuum buffer")
    grid = UniformRadialGrid(0.0, outer, point_count)
    spacing = grid.spacing

    def aligned_index(name: str, radius: float) -> int:
        raw = radius / spacing
        index = int(round(raw))
        tolerance = 8192.0 * np.finfo(np.float64).eps * max(1.0, abs(raw))
        if abs(raw - index) > tolerance:
            raise ValueError(f"{name} is not aligned to the requested grid")
        return index

    start_index = aligned_index("support minimum", parameters.support_minimum)
    end_index = aligned_index("support maximum", parameters.support_maximum)
    step_count = end_index - start_index
    constraint = solve_gr0_initial_slice(
        parameters,
        step_count=step_count,
        method=constraint_method,
    )
    radii = grid.coordinates
    radial_metric = np.ones(point_count, dtype=np.float64)
    angular_k = np.zeros(point_count, dtype=np.float64)
    radial_metric_r = np.zeros(point_count, dtype=np.float64)

    for offset, point in enumerate(constraint.points):
        index = start_index + offset
        radius = radii[index]
        tolerance = 8192.0 * np.finfo(np.float64).eps * max(1.0, radius)
        if abs(point.radius - radius) > tolerance:
            raise ValueError("constraint solution lost grid alignment")
        radial_metric[index] = point.radial_metric
        angular_k[index] = point.angular_extrinsic_curvature
        derivative, _ = gr0_constraint_rhs(
            radius,
            (point.radial_metric, point.angular_extrinsic_curvature),
            parameters,
        )
        radial_metric_r[index] = derivative

    exterior = np.arange(point_count) > end_index
    exterior_r = radii[exterior]
    mass = constraint.outer_mass
    momentum_constant = constraint.outer_k_times_r_cubed
    denominator = (
        1.0
        - 2.0 * mass / exterior_r
        + momentum_constant**2 / exterior_r**4
    )
    if np.any(denominator <= 0.0) or not np.all(np.isfinite(denominator)):
        raise ValueError("GR-0 analytic vacuum exterior lost its positive branch")
    radial_metric[exterior] = denominator**-0.5
    angular_k[exterior] = momentum_constant / exterior_r**3
    denominator_r = (
        2.0 * mass / exterior_r**2
        - 4.0 * momentum_constant**2 / exterior_r**5
    )
    radial_metric_r[exterior] = -0.5 * radial_metric[exterior] ** 3 * denominator_r

    u = np.zeros((point_count, FIELD_COUNT), dtype=np.float64)
    p = np.zeros_like(u)
    q = np.zeros_like(u)
    u[:, 0] = 1.0
    u[:, 2] = radial_metric
    u[:, 3] = radii
    q[:, 2] = radial_metric_r
    q[:, 3] = 1.0

    for index in range(start_index, end_index + 1):
        fields = pulse_fields(radii[index], parameters)
        u[index, 4] = fields["phi"]
        u[index, 5] = fields["chi"]
        p[index, 4] = fields["phi_pi"]
        p[index, 5] = fields["chi_pi"]
        q[index, 4] = fields["phi_r"]
        q[index, 5] = fields["chi_r"]

    positive = radii > 0.0
    h_tr_t = np.zeros(point_count, dtype=np.float64)
    for index in np.flatnonzero(positive):
        h_tr_t[index] = gauge_compatible_metric_time_derivatives(
            radii[index], radial_metric[index], radial_metric_r[index]
        )[1]
    p[:, 1] = h_tr_t / radial_metric**2
    p[:, 2] = 2.0 * radial_metric * angular_k
    p[:, 3] = -radii * angular_k

    # Exact parity and elementary flatness at the centre.
    u[0, ADM_CENTER_PARITIES == -1] = 0.0
    p[0, ADM_CENTER_PARITIES == -1] = 0.0
    q[0, Q_CENTER_PARITIES == -1] = 0.0
    q[0, 3] = u[0, 2]
    state = EvolutionState(u, p, q)
    derivative = SBPFirstDerivative(grid, diagnostic_spatial_order)
    reduction = q - derivative.differentiate(
        u, center_parities=ADM_CENTER_PARITIES
    )
    observables = radial_null_observables(state)
    initially_trapped = (observables.theta_plus[1:] < 0.0) & (
        observables.theta_minus[1:] < 0.0
    )
    return GR0GridInitialData(
        grid=grid,
        parameters=parameters,
        constraint_method=constraint_method,
        state=state,
        support_minimum_index=start_index,
        support_maximum_index=end_index,
        peak_compactness=constraint.peak_compactness,
        peak_radius=constraint.peak_radius,
        outer_mass=mass,
        vacuum_momentum_constant=momentum_constant,
        reduction_constraint_infinity=float(
            np.max(np.abs(reduction), initial=0.0)
        ),
        no_initial_trapped_sphere=not bool(np.any(initially_trapped)),
    )


def make_gr0_center_boundary_projector(
    initial: GR0GridInitialData,
    *,
    fixed_outer_rows: int,
) -> Callable[[float, EvolutionState], EvolutionState]:
    """Return the centre regularity and causally excluded outer-row projector."""

    if not isinstance(initial, GR0GridInitialData):
        raise TypeError("initial must be GR0GridInitialData")
    if (
        isinstance(fixed_outer_rows, bool)
        or not isinstance(fixed_outer_rows, int)
        or fixed_outer_rows < 1
        or fixed_outer_rows >= initial.grid.point_count // 2
    ):
        raise ValueError("fixed_outer_rows must be a small positive integer")
    reference = initial.state

    def projector(_time: float, state: EvolutionState) -> EvolutionState:
        if not isinstance(state, EvolutionState) or state.shape != reference.shape:
            raise ValueError("projected GR-0 state shape differs")
        u = state.u.copy()
        p = state.p.copy()
        q = state.q.copy()
        u[0, ADM_CENTER_PARITIES == -1] = 0.0
        p[0, ADM_CENTER_PARITIES == -1] = 0.0
        q[0, Q_CENTER_PARITIES == -1] = 0.0
        q[0, 3] = u[0, 2]
        u[-fixed_outer_rows:] = reference.u[-fixed_outer_rows:]
        p[-fixed_outer_rows:] = reference.p[-fixed_outer_rows:]
        q[-fixed_outer_rows:] = reference.q[-fixed_outer_rows:]
        return EvolutionState(u, p, q)

    return projector


class GR0EvolutionOperator:
    """Callable method-of-lines right-hand side for GR-0 REF1 calibration."""

    def __init__(
        self,
        grid: UniformRadialGrid,
        *,
        spatial_order: int,
        ko_dissipation: Real = 1.0 / 64.0,
        residual_tolerance: Real = 1.0e-12,
        kinetic_condition_maximum: Real = 1.0e10,
    ) -> None:
        if not isinstance(grid, UniformRadialGrid) or grid.minimum != 0.0:
            raise ValueError("GR-0 calibration requires a centre grid")
        self.grid = grid
        self.derivative = SBPFirstDerivative(grid, spatial_order)
        self.ko_dissipation = _finite("ko_dissipation", ko_dissipation)
        if self.ko_dissipation < 0.0:
            raise ValueError("ko_dissipation must be nonnegative")
        self.residual_tolerance = _finite(
            "residual_tolerance", residual_tolerance, positive=True
        )
        self.kinetic_condition_maximum = _finite(
            "kinetic_condition_maximum",
            kinetic_condition_maximum,
            positive=True,
        )

    def __call__(self, _time: float, state: EvolutionState) -> EvolutionRHS:
        if not isinstance(state, EvolutionState) or state.shape != (
            self.grid.point_count,
            FIELD_COUNT,
        ):
            raise ValueError("GR-0 RHS state shape differs")
        p_r = self.derivative.differentiate(
            state.p, center_parities=ADM_CENTER_PARITIES
        )
        q_r = self.derivative.differentiate(
            state.q, center_parities=Q_CENTER_PARITIES
        )
        source = solve_gr0_grid_accelerations(
            state.u[1:],
            state.p[1:],
            state.q[1:],
            p_r[1:],
            q_r[1:],
            self.grid.coordinates[1:],
            residual_tolerance=self.residual_tolerance,
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
        dp = acceleration + dissipation
        dq = p_r
        reduction = state.q - self.derivative.differentiate(
            state.u, center_parities=ADM_CENTER_PARITIES
        )
        coordinate_speed = np.abs(state.u[:, 1]) + (
            float(ALL_CONE_LOCAL_SPEED_BOUND)
            * state.u[:, 0]
            / state.u[:, 2]
        )
        diagnostics = {
            "source_residual_infinity": source.residual_infinity,
            "kinetic_condition_infinity": source.kinetic_condition_infinity_maximum,
            "source_refinement_iterations": source.refinement_iterations,
            "source_residual_decreased_monotonically": (
                source.residual_decreased_monotonically
            ),
            "acceleration_infinity": float(
                np.max(np.abs(acceleration), initial=0.0)
            ),
            "center_acceleration_estimator_infinity": center.estimator_infinity,
            "hamiltonian_constraint_infinity": float(
                np.max(np.abs(source.hamiltonian_constraint), initial=0.0)
            ),
            "momentum_constraint_infinity": float(
                np.max(np.abs(source.momentum_constraint), initial=0.0)
            ),
            "gauge_constraint_infinity": float(
                np.max(np.abs(source.gauge_constraint), initial=0.0)
            ),
            "reduction_constraint_infinity": float(
                np.max(np.abs(reduction), initial=0.0)
            ),
            "ricci_scalar_infinity": float(
                np.max(np.abs(source.ricci_scalar), initial=0.0)
            ),
            "ricci_squared_infinity": float(
                np.max(np.abs(source.ricci_squared), initial=0.0)
            ),
            "coordinate_speed_upper": float(
                np.max(coordinate_speed, initial=0.0)
            ),
            "minimum_lapse": float(np.min(state.u[:, 0])),
            "minimum_radial_metric": float(np.min(state.u[:, 2])),
            "minimum_areal_radius_away_from_center": float(
                np.min(state.u[1:, 3])
            ),
        }
        return EvolutionRHS(state.p, dp, dq, diagnostics)


__all__ = [
    "GR0EvolutionOperator",
    "GR0GridInitialData",
    "GR0SliceObservables",
    "Q_CENTER_PARITIES",
    "construct_gr0_grid_initial_data",
    "make_gr0_center_boundary_projector",
    "radial_null_observables",
]
