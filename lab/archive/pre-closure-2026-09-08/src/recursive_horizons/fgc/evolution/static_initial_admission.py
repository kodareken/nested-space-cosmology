"""Static PROTO4 admission for complete FGC-QR initial grids.

The first protocol versions described five FGC-QR holdout perturbations whose
matter amplitude is inherited from a later GR-0 calibration.  PROTO3 requires
their initial-data premises to be checked before an amplitude becomes eligible
for that calibration.  Consequently the pre-calibration ledger must evaluate
the Cartesian product of every declared amplitude and every held-out modifier;
checking only the five cases after selecting an amplitude would leak a dynamic
outcome back into input admission.

This module embeds the already certified radial FGC-QR constraint solution in
the complete regular-centre grid used by the numerical engine.  It supplies
the exact inner and outer vacuum continuations, frozen gauge-compatible time
derivatives, analytic radial derivatives, physical-constraint residuals, and a
discrete reduction-constraint diagnostic.  It never advances a time step,
reads a run namespace, or classifies collapse, activation, or defocusing.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite, sqrt
from numbers import Real

import numpy as np

from ..initial_data_family import (
    InitialDataParameters,
    InitialDataSolution,
    compact_family_fields,
    gauge_compatible_metric_time_derivatives,
    polar_areal_constraint_coefficients,
    polar_areal_constraint_residuals,
    solve_initial_data,
)
from .numerical_engine import EvolutionState, SBPFirstDerivative, UniformRadialGrid
from .vectorized_source import ADM_CENTER_PARITIES


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


@dataclass(frozen=True, slots=True)
class FGCQRGridInitialData:
    """One complete regular-centre FGC-QR initial hypersurface."""

    grid: UniformRadialGrid
    parameters: InitialDataParameters
    constraint_method: str
    constraint_solution: InitialDataSolution
    state: EvolutionState
    support_minimum_index: int
    support_maximum_index: int
    peak_compactness: float
    peak_radius: float
    outer_mass: float
    vacuum_momentum_constant: float
    physical_constraint_residual_infinity: float
    reduction_constraint_infinity: float
    minimum_radial_metric: float
    minimum_effective_planck_coefficient: float
    minimum_vacuum_metric_denominator: float
    exact_inner_vacuum_buffer: bool
    exact_outer_vacuum_buffer: bool
    finite_mass: bool
    regular_center: bool
    no_initial_trapped_sphere: bool
    compactness_inside_protocol_window: bool

    def __post_init__(self) -> None:
        if not isinstance(self.grid, UniformRadialGrid):
            raise TypeError("grid must be UniformRadialGrid")
        if not isinstance(self.parameters, InitialDataParameters):
            raise TypeError("parameters must be InitialDataParameters")
        if not isinstance(self.constraint_solution, InitialDataSolution):
            raise TypeError("constraint_solution must be InitialDataSolution")
        if not isinstance(self.state, EvolutionState) or self.state.shape != (
            self.grid.point_count,
            FIELD_COUNT,
        ):
            raise ValueError("FGC-QR initial state does not match its six-field grid")
        if self.constraint_method not in {"RK4", "SSPRK3"}:
            raise ValueError("unknown FGC-QR constraint method")
        if self.constraint_solution.method != self.constraint_method:
            raise ValueError("constraint method and embedded solution differ")
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
            "physical_constraint_residual_infinity",
            "reduction_constraint_infinity",
            "minimum_radial_metric",
            "minimum_effective_planck_coefficient",
            "minimum_vacuum_metric_denominator",
        ):
            value = _finite(name, getattr(self, name))
            if name in {
                "peak_radius",
                "outer_mass",
                "minimum_radial_metric",
                "minimum_effective_planck_coefficient",
                "minimum_vacuum_metric_denominator",
            } and value <= 0.0:
                raise ValueError(f"{name} must be positive")
            if name in {
                "physical_constraint_residual_infinity",
                "reduction_constraint_infinity",
            } and value < 0.0:
                raise ValueError(f"{name} must be nonnegative")
            object.__setattr__(self, name, value)
        for name in (
            "exact_inner_vacuum_buffer",
            "exact_outer_vacuum_buffer",
            "finite_mass",
            "regular_center",
            "no_initial_trapped_sphere",
            "compactness_inside_protocol_window",
        ):
            if not isinstance(getattr(self, name), bool):
                raise TypeError(f"{name} must be bool")


def _aligned_index(grid: UniformRadialGrid, name: str, radius: float) -> int:
    raw = radius / grid.spacing
    index = int(round(raw))
    tolerance = 8192.0 * np.finfo(np.float64).eps * max(1.0, abs(raw))
    if abs(raw - index) > tolerance:
        raise ValueError(f"{name} is not aligned to the requested grid")
    return index


def _null_trapping_exists(state: EvolutionState) -> bool:
    alpha = state.u[:, 0]
    shift = state.u[:, 1]
    radial_metric = state.u[:, 2]
    areal = state.u[:, 3]
    normal = (state.p[:, 3] - shift * state.q[:, 3]) / alpha
    spatial = state.q[:, 3] / radial_metric
    positive = areal > 0.0
    theta_plus = np.zeros_like(areal)
    theta_minus = np.zeros_like(areal)
    theta_plus[positive] = 2.0 * (normal[positive] + spatial[positive]) / areal[positive]
    theta_minus[positive] = 2.0 * (normal[positive] - spatial[positive]) / areal[positive]
    return bool(np.any((theta_plus[1:] < 0.0) & (theta_minus[1:] < 0.0)))


def construct_fgcqr_grid_initial_data(
    parameters: InitialDataParameters,
    *,
    point_count: int,
    constraint_method: str = "RK4",
    diagnostic_spatial_order: int = 4,
    maximum_constraint_residual: Real = 1.0e-10,
) -> FGCQRGridInitialData:
    """Embed one solved FGC-QR slice in a complete uniform radial grid."""

    if not isinstance(parameters, InitialDataParameters):
        raise TypeError("parameters must be InitialDataParameters")
    residual_limit = _finite(
        "maximum_constraint_residual", maximum_constraint_residual, positive=True
    )
    grid = UniformRadialGrid(0.0, parameters.outer_radius, point_count)
    start_index = _aligned_index(grid, "combined support minimum", parameters.support_minimum)
    end_index = _aligned_index(grid, "combined support maximum", parameters.support_maximum)
    step_count = end_index - start_index
    solution = solve_initial_data(
        parameters,
        step_count=step_count,
        method=constraint_method,
        maximum_constraint_residual=residual_limit,
    )
    radii = grid.coordinates
    radial_metric = np.ones(point_count, dtype=np.float64)
    angular_k = np.zeros(point_count, dtype=np.float64)
    radial_metric_r = np.zeros(point_count, dtype=np.float64)
    angular_k_r = np.zeros(point_count, dtype=np.float64)

    for offset, point in enumerate(solution.points):
        index = start_index + offset
        tolerance = 8192.0 * np.finfo(np.float64).eps * max(1.0, radii[index])
        if abs(point.radius - radii[index]) > tolerance:
            raise ValueError("FGC-QR constraint solution lost full-grid alignment")
        radial_metric[index] = point.radial_metric
        angular_k[index] = point.angular_extrinsic_curvature
        radial_metric_r[index] = point.radial_metric_derivative
        angular_k_r[index] = point.angular_extrinsic_curvature_derivative

    exterior = np.arange(point_count) > end_index
    exterior_r = radii[exterior]
    mass = solution.outer_mass
    momentum_constant = solution.vacuum_momentum_constant
    denominator = (
        1.0
        - 2.0 * mass / exterior_r
        + momentum_constant**2 / exterior_r**4
    )
    if np.any(denominator <= 0.0) or not np.all(np.isfinite(denominator)):
        raise ValueError("FGC-QR analytic vacuum exterior lost its positive branch")
    radial_metric[exterior] = denominator**-0.5
    angular_k[exterior] = momentum_constant / exterior_r**3
    denominator_r = (
        2.0 * mass / exterior_r**2
        - 4.0 * momentum_constant**2 / exterior_r**5
    )
    radial_metric_r[exterior] = -0.5 * radial_metric[exterior] ** 3 * denominator_r
    angular_k_r[exterior] = -3.0 * momentum_constant / exterior_r**4

    u = np.zeros((point_count, FIELD_COUNT), dtype=np.float64)
    p = np.zeros_like(u)
    q = np.zeros_like(u)
    u[:, 0] = 1.0
    u[:, 2] = radial_metric
    u[:, 3] = radii
    q[:, 2] = radial_metric_r
    q[:, 3] = 1.0

    for index in range(start_index, end_index + 1):
        fields = compact_family_fields(radii[index], parameters)
        u[index, 4] = fields["phi"]
        u[index, 5] = fields["chi"]
        p[index, 4] = fields["phi_pi"]
        p[index, 5] = fields["chi_pi"]
        q[index, 4] = fields["phi_r"]
        q[index, 5] = fields["chi_r"]

    for index in range(1, point_count):
        h_tr_t = gauge_compatible_metric_time_derivatives(
            radii[index], radial_metric[index], radial_metric_r[index]
        )[1]
        p[index, 1] = h_tr_t / radial_metric[index] ** 2
    p[:, 2] = 2.0 * radial_metric * angular_k
    p[:, 3] = -radii * angular_k

    # Exact parity and elementary flatness at the centre.
    u[0, ADM_CENTER_PARITIES == -1] = 0.0
    p[0, ADM_CENTER_PARITIES == -1] = 0.0
    q[0, Q_CENTER_PARITIES == -1] = 0.0
    q[0, 3] = u[0, 2]
    state = EvolutionState(u, p, q)

    maximum_physical_residual = 0.0
    action = parameters.action
    for index in range(1, point_count):
        fields = compact_family_fields(radii[index], parameters)
        coefficients = polar_areal_constraint_coefficients(
            radius=radii[index],
            radial_metric=radial_metric[index],
            angular_extrinsic_curvature=angular_k[index],
            phi=fields["phi"],
            phi_r=fields["phi_r"],
            phi_rr=fields["phi_rr"],
            phi_pi=fields["phi_pi"],
            phi_pi_r=fields["phi_pi_r"],
            chi_r=fields["chi_r"],
            chi_pi=fields["chi_pi"],
            planck_mass=action.planck_mass,
            beta=action.beta,
            scalar_mass=action.scalar_mass,
            quartic_coupling=action.quartic_coupling,
            eta=action.eta,
        )
        residuals = polar_areal_constraint_residuals(
            coefficients,
            radial_metric_derivative=radial_metric_r[index],
            angular_extrinsic_curvature_derivative=angular_k_r[index],
        )
        maximum_physical_residual = max(
            maximum_physical_residual,
            *(abs(float(value)) for value in residuals),
        )
    if maximum_physical_residual > residual_limit:
        raise ValueError("complete FGC-QR grid violates the physical-constraint limit")

    derivative = SBPFirstDerivative(grid, diagnostic_spatial_order)
    reduction = q - derivative.differentiate(u, center_parities=ADM_CENTER_PARITIES)
    # The static conclusion uses r<=24.  The analytic outer boundary closure is
    # deliberately excluded here, exactly as it will be by BND2 dynamically.
    measurement = radii <= 24.0
    reduction_infinity = float(np.max(np.abs(reduction[measurement]), initial=0.0))
    inner = np.arange(point_count) < start_index
    outer = np.arange(point_count) > end_index
    inner_exact = bool(
        np.all(u[inner, 2] == 1.0)
        and np.all(u[inner, 4:] == 0.0)
        and np.all(p[inner, 2:] == 0.0)
        and np.all(q[inner, 2] == 0.0)
        and np.all(q[inner, 4:] == 0.0)
    )
    outer_exact = bool(
        np.all(u[outer, 4:] == 0.0)
        and np.all(p[outer, 4:] == 0.0)
        and np.all(q[outer, 4:] == 0.0)
    )
    return FGCQRGridInitialData(
        grid=grid,
        parameters=parameters,
        constraint_method=constraint_method,
        constraint_solution=solution,
        state=state,
        support_minimum_index=start_index,
        support_maximum_index=end_index,
        peak_compactness=solution.peak_compactness,
        peak_radius=solution.peak_radius,
        outer_mass=mass,
        vacuum_momentum_constant=momentum_constant,
        physical_constraint_residual_infinity=maximum_physical_residual,
        reduction_constraint_infinity=reduction_infinity,
        minimum_radial_metric=float(np.min(radial_metric)),
        minimum_effective_planck_coefficient=solution.minimum_effective_planck_coefficient,
        minimum_vacuum_metric_denominator=solution.minimum_vacuum_metric_denominator,
        exact_inner_vacuum_buffer=inner_exact and solution.exact_inner_vacuum_buffer,
        exact_outer_vacuum_buffer=outer_exact and solution.exact_outer_vacuum_buffer,
        finite_mass=solution.finite_mass and isfinite(mass),
        regular_center=bool(
            solution.regular_center
            and u[0, 3] == 0.0
            and u[0, 2] == q[0, 3]
            and abs(radial_metric[0] - 1.0) == 0.0
        ),
        no_initial_trapped_sphere=(
            solution.no_initial_trapped_sphere and not _null_trapping_exists(state)
        ),
        compactness_inside_protocol_window=solution.compactness_inside_protocol_window,
    )


__all__ = ["FGCQRGridInitialData", "construct_fgcqr_grid_initial_data"]
