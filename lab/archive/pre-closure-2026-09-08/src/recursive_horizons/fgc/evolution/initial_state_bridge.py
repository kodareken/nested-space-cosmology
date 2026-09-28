"""Bridge solved ID1 hypersurfaces to complete REF1 second jets.

ID1 solves only the hypersurface constraints.  A principal-health calculation
also needs the spatial second derivatives of the metric and the six time
accelerations selected by the complete REF1 equations.  This module supplies
that missing, explicitly numerical bridge without changing either equation
set.

The radial constraint solution already stores ``lambda_r`` and ``k_r`` at
every point.  Their derivatives are reconstructed with a centred five-point
operator; the three-point result is retained as an error indicator.  All
other second-jet entries follow by differentiating the frozen polar--areal,
maximal-slice identities.  Finally the established same-evaluator nonlinear
source solver determines ``partial_t^2 u``.

This is an initial-slice adapter.  It is not a time evolution, a continuum
interpolation theorem, or a multidirectional hyperbolicity certificate.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite
from numbers import Real
from typing import Sequence

import numpy as np

from ..initial_data_family import (
    ConstraintPoint,
    InitialDataSolution,
    compact_family_fields,
    gauge_compatible_metric_time_derivatives,
)
from ..reference_connection import flat_spherical_annulus_reference
from ..spherical_reduction import SphericalState
from .floating_jet import FloatJet2, scalar_primal
from .nonlinear_source import (
    AccelerationSolveResult,
    NonlinearSolverConfig,
    acceleration_jacobian,
    parameter_vector,
    solve_accelerations,
    state_from_parameter_vector,
)


@dataclass(frozen=True, slots=True)
class SpatialSecondDerivativeDiagnostics:
    """Two independent derivative reconstructions at one support point."""

    radial_step: float
    lambda_rr_five_point: float
    lambda_rr_three_point: float
    k_rr_five_point: float
    k_rr_three_point: float
    lambda_rr_estimator: float
    k_rr_estimator: float


@dataclass(frozen=True, slots=True)
class InitialSecondJet:
    """One complete initial REF1 jet and its reconstruction diagnostics."""

    point_index: int
    coordinate_radius: float
    unaccelerated_state: SphericalState
    accelerated_state: SphericalState
    spatial_diagnostics: SpatialSecondDerivativeDiagnostics
    acceleration_solve: AccelerationSolveResult


def _finite(name: str, value: object, *, positive: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a finite real scalar")
    answer = float(value)
    if not isfinite(answer):
        raise ValueError(f"{name} must be finite")
    if positive and answer <= 0.0:
        raise ValueError(f"{name} must be positive")
    return answer


def _uniform_step(points: Sequence[ConstraintPoint]) -> float:
    if len(points) < 5:
        raise ValueError("initial solution needs at least five support points")
    step = _finite(
        "radial step", points[1].radius - points[0].radius, positive=True
    )
    tolerance = 8192.0 * np.finfo(np.float64).eps * max(
        1.0, abs(points[-1].radius)
    )
    if any(
        abs((right.radius - left.radius) - step) > tolerance
        for left, right in zip(points, points[1:], strict=False)
    ):
        raise ValueError("initial solution points are not uniformly spaced")
    return step


def _derivative_of_stored_first_derivative(
    points: Sequence[ConstraintPoint],
    index: int,
    attribute: str,
    step: float,
) -> tuple[float, float, float]:
    values = [
        _finite(attribute, getattr(points[index + offset], attribute))
        for offset in (-2, -1, 0, 1, 2)
    ]
    five = (values[0] - 8.0 * values[1] + 8.0 * values[3] - values[4]) / (
        12.0 * step
    )
    three = (values[3] - values[1]) / (2.0 * step)
    estimator = abs(five - three)
    if not all(isfinite(value) for value in (five, three, estimator)):
        raise ValueError("spatial second-derivative reconstruction is nonfinite")
    return five, three, estimator


def unaccelerated_initial_second_jet(
    solution: InitialDataSolution,
    point_index: int,
) -> tuple[SphericalState, SpatialSecondDerivativeDiagnostics]:
    """Construct all initial second-jet slots except ``partial_t^2 u``.

    The selected point must have two stored neighbors on each side.  Centre
    limits and the analytic exterior belong to their dedicated formulations;
    this bridge deliberately handles only the smooth compact-support chart.
    """

    if not isinstance(solution, InitialDataSolution):
        raise TypeError("solution must be an InitialDataSolution")
    if isinstance(point_index, bool) or not isinstance(point_index, int):
        raise TypeError("point_index must be an integer")
    points = solution.points
    if point_index < 2 or point_index >= len(points) - 2:
        raise ValueError("point_index needs two support neighbors on each side")
    step = _uniform_step(points)
    point = points[point_index]
    lambda_rr, lambda_rr_three, lambda_error = (
        _derivative_of_stored_first_derivative(
            points,
            point_index,
            "radial_metric_derivative",
            step,
        )
    )
    k_rr, k_rr_three, k_error = _derivative_of_stored_first_derivative(
        points,
        point_index,
        "angular_extrinsic_curvature_derivative",
        step,
    )

    radius = _finite("coordinate radius", point.radius, positive=True)
    radial_metric = _finite("radial metric", point.radial_metric, positive=True)
    lambda_r = _finite(
        "radial metric derivative", point.radial_metric_derivative
    )
    angular_k = _finite(
        "angular extrinsic curvature", point.angular_extrinsic_curvature
    )
    k_r = _finite(
        "angular extrinsic curvature derivative",
        point.angular_extrinsic_curvature_derivative,
    )
    fields = compact_family_fields(radius, solution.parameters)

    h_tt_t, h_tr_t = gauge_compatible_metric_time_derivatives(
        radius, radial_metric, lambda_r
    )
    recorded_tolerance = 16384.0 * np.finfo(np.float64).eps * max(
        1.0, abs(h_tr_t)
    )
    if (
        abs(float(h_tt_t) - point.h_tt_time_derivative) > recorded_tolerance
        or abs(float(h_tr_t) - point.h_tr_time_derivative) > recorded_tolerance
    ):
        raise ValueError("stored ID1 point violates its gauge-time-derivative identity")

    # h_tr,t = (lambda^2-1)/(2r) + lambda_r/(4lambda).
    h_tr_tr = (
        radial_metric * lambda_r / radius
        - (radial_metric**2 - 1.0) / (2.0 * radius**2)
        + (
            lambda_rr * radial_metric - lambda_r**2
        )
        / (4.0 * radial_metric**2)
    )
    action = solution.parameters.action
    state = SphericalState(
        h_tt=FloatJet2(-1.0, dt=h_tt_t),
        h_tr=FloatJet2(0.0, dt=h_tr_t, dtr=h_tr_tr),
        h_rr=FloatJet2(
            radial_metric**2,
            dt=4.0 * radial_metric**2 * angular_k,
            dr=2.0 * radial_metric * lambda_r,
            dtr=(
                8.0 * radial_metric * lambda_r * angular_k
                + 4.0 * radial_metric**2 * k_r
            ),
            drr=2.0 * lambda_r**2 + 2.0 * radial_metric * lambda_rr,
        ),
        areal_radius=FloatJet2(
            radius,
            dt=-radius * angular_k,
            dr=1.0,
            dtr=-angular_k - radius * k_r,
        ),
        phi=FloatJet2(
            fields["phi"],
            dt=fields["phi_pi"],
            dr=fields["phi_r"],
            dtr=fields["phi_pi_r"],
            drr=fields["phi_rr"],
        ),
        chi=FloatJet2(
            fields["chi"],
            dt=fields["chi_pi"],
            dr=fields["chi_r"],
            dtr=fields["chi_pi_r"],
            drr=fields["chi_rr"],
        ),
        planck_mass=action.planck_mass,
        beta=action.beta,
        mu=action.scalar_mass,
        g4=action.quartic_coupling,
        eta=action.eta,
        branch="FGC-QR",
    )
    diagnostics = SpatialSecondDerivativeDiagnostics(
        radial_step=step,
        lambda_rr_five_point=lambda_rr,
        lambda_rr_three_point=lambda_rr_three,
        k_rr_five_point=k_rr,
        k_rr_three_point=k_rr_three,
        lambda_rr_estimator=lambda_error,
        k_rr_estimator=k_error,
    )
    return state, diagnostics


def solve_initial_second_jet(
    solution: InitialDataSolution,
    point_index: int,
    *,
    warm_start: Sequence[Real] | None = None,
    solver_config: NonlinearSolverConfig | None = None,
) -> InitialSecondJet:
    """Complete one ID1 support jet with the six REF1 accelerations."""

    state, spatial = unaccelerated_initial_second_jet(solution, point_index)
    settings = solver_config or NonlinearSolverConfig(
        residual_tolerance=1.0e-12,
        maximum_iterations=16,
        condition_number_maximum=1.0e10,
        normalized_branch_displacement_maximum=1.0 / 16.0,
        parameter_half_width=1.0,
        acceleration_half_width=4.0,
    )
    radius = _finite(
        "coordinate radius", scalar_primal(state.areal_radius.value), positive=True
    )
    reference = flat_spherical_annulus_reference(radial_domain_minimum=0.5)
    if warm_start is None:
        # The ACT1/REF1 equations are linear in derivatives twice along any
        # one coordinate.  Use the same-evaluator acceleration Jacobian to
        # form the initial linear predictor.  This avoids treating the
        # arbitrary zero vector as a branch-continuation reference while
        # leaving the safeguarded solver responsible for acceptance.
        residual, jacobian, _condition = acceleration_jacobian(
            state,
            (0.0,) * 6,
            reference=reference,
            coordinate_radius=radius,
            tilde_normal_factor=4,
            hat_normal_factor=9,
        )
        try:
            initial_predictor = np.linalg.solve(jacobian, -residual)
        except np.linalg.LinAlgError as exc:
            raise ValueError("initial REF1 acceleration predictor is singular") from exc
        if not np.all(np.isfinite(initial_predictor)):
            raise ValueError("initial REF1 acceleration predictor is nonfinite")
        warm_start = tuple(float(value) for value in initial_predictor)

    solved = solve_accelerations(
        state,
        reference=reference,
        coordinate_radius=radius,
        tilde_normal_factor=4,
        hat_normal_factor=9,
        branch_center=state,
        warm_start=warm_start,
        config=settings,
    )
    accelerated = state_from_parameter_vector(
        state,
        parameter_vector(state),
        accelerations=solved.accelerations,
    )
    return InitialSecondJet(
        point_index=point_index,
        coordinate_radius=radius,
        unaccelerated_state=state,
        accelerated_state=accelerated,
        spatial_diagnostics=spatial,
        acceleration_solve=solved,
    )
