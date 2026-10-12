"""Constraint-compatible finite-mass spherical FGC-QR initial data.

The declared SF1 slice is maximal and polar--areal,

``alpha=1, shift=0, R=r, K^r_r=-2*k, K^theta_theta=K^phi_phi=k``.

On this slice the complete unredefined ACT1/VAR1 Hamiltonian constraint is
affine in ``partial_r lambda`` and the radial momentum constraint is affine in
``partial_r k``.  The small kernel below is the resulting specialization of
the full tensor evaluator.  It is not a replacement set of field equations:
exact controls compare it directly with
``physical_constraint_projections`` before any generated family is trusted.

The regular centre is an exact Minkowski buffer.  Compact matter and regulator
profiles are integrated only across their combined support.  Outside that
support the exact vacuum continuation

``k=J/r^3`` and ``lambda^-2=1-2*M/r+J^2/r^4``

provides an asymptotically flat finite-Misner--Sharp-mass end.  Two independent
explicit integrators are supplied for convergence checks.  This module
constructs initial hypersurfaces; it does not evolve them, establish a healthy
run domain, or authorize a holdout calculation.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import isfinite, sqrt
from numbers import Real
from typing import Any, Callable, Sequence

import numpy as np

from .initial_data_preflight import compact_bump_with_derivatives
from .modified_harmonic_constraints import physical_constraint_projections
from .modified_harmonic_reference import modified_harmonic_gauge_constraint
from .reference_connection import flat_spherical_annulus_reference
from .spherical_reduction import Jet2, SphericalState


Q = Fraction


STOP_OUTCOME_LABELS = {
    "invalid_input": "invalid_implementation_or_nonconverged_run",
    "singular_constraint_jacobian": "stopped_branch_loss",
    "constraint_residual_loss": "stopped_constraint_loss",
    "negative_metric_factor": "stopped_branch_loss",
    "nonfinite_continuation": "invalid_implementation_or_nonconverged_run",
    "lost_asymptotic_condition": "stopped_constraint_loss",
    "compactness_window_loss": "stopped_constraint_loss",
}


def _finite(name: str, value: object, *, positive: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a finite real scalar")
    answer = float(value)
    if not isfinite(answer):
        raise ValueError(f"{name} must be finite")
    if positive and answer <= 0.0:
        raise ValueError(f"{name} must be positive")
    return answer


def _finite_pair(name: str, value: Sequence[Real]) -> tuple[float, float]:
    if isinstance(value, (str, bytes)) or len(value) != 2:
        raise ValueError(f"{name} must contain two entries")
    return _finite(f"{name}[0]", value[0]), _finite(f"{name}[1]", value[1])


@dataclass(frozen=True, slots=True)
class FGCQRActionParameters:
    """The dimensionful code-unit parameters consumed by the ID1 kernel."""

    planck_mass: float = 2.0
    beta: float = -0.25
    scalar_mass: float = 3.0
    quartic_coupling: float = 0.5
    eta: float = 0.5

    def __post_init__(self) -> None:
        for name in (
            "planck_mass",
            "beta",
            "scalar_mass",
            "quartic_coupling",
            "eta",
        ):
            value = _finite(name, getattr(self, name))
            object.__setattr__(self, name, value)
        if (
            self.planck_mass <= 0.0
            or self.scalar_mass <= 0.0
            or self.quartic_coupling <= 0.0
            or self.beta == 0.0
            or self.eta <= 0.0
        ):
            raise ValueError("FGC-QR action parameters lie outside the declared branch")


@dataclass(frozen=True, slots=True)
class InitialDataParameters:
    """One member of the compact SF1 initial-data family."""

    chi_amplitude: float
    chi_half_width: float
    phi_amplitude: float
    center: float = 12.0
    phi_half_width: float = 2.0
    outer_radius: float = 128.0
    action: FGCQRActionParameters = FGCQRActionParameters()

    def __post_init__(self) -> None:
        for name in (
            "chi_amplitude",
            "chi_half_width",
            "center",
            "phi_half_width",
            "outer_radius",
        ):
            value = _finite(name, getattr(self, name), positive=True)
            object.__setattr__(self, name, value)
        phi_amplitude = _finite("phi_amplitude", self.phi_amplitude)
        if phi_amplitude < 0.0:
            raise ValueError("phi_amplitude must be nonnegative")
        object.__setattr__(self, "phi_amplitude", phi_amplitude)
        if not isinstance(self.action, FGCQRActionParameters):
            raise TypeError("action must be FGCQRActionParameters")
        if self.center <= max(self.chi_half_width, self.phi_half_width):
            raise ValueError("compact profiles must retain a strict centre buffer")
        if self.outer_radius <= self.support_maximum:
            raise ValueError("outer radius must retain a strict vacuum buffer")

    @property
    def support_minimum(self) -> float:
        return self.center - max(self.chi_half_width, self.phi_half_width)

    @property
    def support_maximum(self) -> float:
        return self.center + max(self.chi_half_width, self.phi_half_width)

    @property
    def center_buffer(self) -> float:
        return self.support_minimum

    @property
    def outer_vacuum_buffer(self) -> float:
        return self.outer_radius - self.support_maximum


@dataclass(frozen=True, slots=True)
class ConstraintCoefficients:
    """Affine physical-constraint coefficients in ``lambda_r`` and ``k_r``."""

    hamiltonian_constant: Any
    hamiltonian_lambda_r: Any
    momentum_constant: Any
    momentum_k_r: Any


@dataclass(frozen=True, slots=True)
class ConstraintPoint:
    radius: float
    radial_metric: float
    angular_extrinsic_curvature: float
    radial_metric_derivative: float
    angular_extrinsic_curvature_derivative: float
    hamiltonian_residual: float
    momentum_residual: float
    jacobian_determinant: float
    jacobian_condition_infinity: float
    compactness: float
    misner_sharp_mass: float
    h_tt_time_derivative: float
    h_tr_time_derivative: float


@dataclass(frozen=True, slots=True)
class InitialDataSolution:
    method: str
    step_count: int
    parameters: InitialDataParameters
    points: tuple[ConstraintPoint, ...]
    peak_compactness: float
    peak_radius: float
    outer_mass: float
    vacuum_momentum_constant: float
    outer_radial_metric: float
    outer_angular_extrinsic_curvature: float
    minimum_abs_constraint_jacobian_determinant: float
    maximum_constraint_jacobian_condition_infinity: float
    maximum_constraint_residual_infinity: float
    minimum_effective_planck_coefficient: float
    minimum_vacuum_metric_denominator: float
    finite_mass: bool
    regular_center: bool
    exact_inner_vacuum_buffer: bool
    exact_outer_vacuum_buffer: bool
    no_initial_trapped_sphere: bool
    compactness_inside_protocol_window: bool


class InitialDataSolveStop(RuntimeError):
    """Typed fail-closed termination of the ID1 constraint continuation."""

    def __init__(
        self,
        reason: str,
        message: str,
        diagnostics: dict[str, Any] | None = None,
    ) -> None:
        if reason not in STOP_OUTCOME_LABELS:
            raise ValueError("unknown ID1 stop reason")
        super().__init__(message)
        self.reason = reason
        self.outcome_label = STOP_OUTCOME_LABELS[reason]
        self.diagnostics = dict(diagnostics or {})


def compact_family_fields(
    radius: Real,
    parameters: InitialDataParameters,
) -> dict[str, float]:
    """Return the two compact scalar profiles and derivatives on the slice."""

    r = _finite("radius", radius)
    if r < 0.0:
        raise ValueError("radius must be nonnegative")
    chi_bump, chi_bump_r, chi_bump_rr = compact_bump_with_derivatives(
        r,
        center=parameters.center,
        half_width=parameters.chi_half_width,
    )
    phi_bump, phi_bump_r, phi_bump_rr = compact_bump_with_derivatives(
        r,
        center=parameters.center,
        half_width=parameters.phi_half_width,
    )
    if r == 0.0:
        if any(
            value != 0.0
            for value in (
                chi_bump,
                chi_bump_r,
                chi_bump_rr,
                phi_bump,
                phi_bump_r,
                phi_bump_rr,
            )
        ):
            raise ValueError("declared compact profiles reached the centre")
        chi = chi_r = chi_rr = chi_pi = chi_pi_r = 0.0
    else:
        amplitude = parameters.chi_amplitude
        chi = amplitude * chi_bump / r
        chi_r = amplitude * (chi_bump_r / r - chi_bump / r**2)
        chi_rr = amplitude * (
            chi_bump_rr / r
            - 2.0 * chi_bump_r / r**2
            + 2.0 * chi_bump / r**3
        )
        # Frozen future-ingoing orientation: d_t(r*chi)=d_r(r*chi).
        chi_pi = amplitude * chi_bump_r / r
        chi_pi_r = amplitude * (chi_bump_rr / r - chi_bump_r / r**2)
    phi = parameters.phi_amplitude * phi_bump
    phi_r = parameters.phi_amplitude * phi_bump_r
    phi_rr = parameters.phi_amplitude * phi_bump_rr
    return {
        "phi": phi,
        "phi_r": phi_r,
        "phi_rr": phi_rr,
        "phi_pi": 0.0,
        "phi_pi_r": 0.0,
        "chi": chi,
        "chi_r": chi_r,
        "chi_rr": chi_rr,
        "chi_pi": chi_pi,
        "chi_pi_r": chi_pi_r,
    }


def polar_areal_constraint_coefficients(
    *,
    radius: Any,
    radial_metric: Any,
    angular_extrinsic_curvature: Any,
    phi: Any,
    phi_r: Any,
    phi_rr: Any,
    phi_pi: Any,
    phi_pi_r: Any,
    chi_r: Any,
    chi_pi: Any,
    planck_mass: Any,
    beta: Any,
    scalar_mass: Any,
    quartic_coupling: Any,
    eta: Any,
) -> ConstraintCoefficients:
    """Return the exact affine coefficients of the two FGC-QR constraints.

    Arithmetic is deliberately generic: ordinary binary64 inputs support the
    numerical family solver, while ``Fraction`` inputs provide exact rational
    cross-checks against the complete tensor evaluator.
    """

    r = radius
    L = radial_metric
    K = angular_extrinsic_curvature
    p = phi
    pr = phi_r
    prr = phi_rr
    pp = phi_pi
    ppr = phi_pi_r
    xr = chi_r
    xp = chi_pi
    MP = planck_mass
    b = beta
    mu = scalar_mass
    g = quartic_coupling
    e = eta
    if r <= 0 or L <= 0 or MP <= 0 or mu <= 0 or g <= 0 or e <= 0 or b == 0:
        raise ValueError("polar-areal FGC-QR constraint point is outside its domain")
    if MP**2 + b * p**2 <= 0:
        raise ValueError("effective Planck coefficient is not positive")

    h0_numerator = (
        -48 * K**3 * L**5 * e * p * pp * r**2
        + 12 * K**2 * L**5 * MP**2 * r**2
        + 12 * K**2 * L**5 * b * p**2 * r**2
        - 32 * K**2 * L**3 * e * p * pr * r
        + 8 * K**2 * L**3 * e * p * prr * r**2
        + 8 * K**2 * L**3 * e * pr**2 * r**2
        - 16 * K * L**5 * e * p * pp
        + 16 * K * L**3 * e * p * pp
        - 4 * L**5 * MP**2
        - 4 * L**5 * b * p**2
        + L**5 * g * p**4 * r**2
        + 2 * L**5 * mu**2 * p**2 * r**2
        + 2 * L**5 * pp**2 * r**2
        + 2 * L**5 * r**2 * xp**2
        + 4 * L**3 * MP**2
        + 4 * L**3 * b * p**2
        + 16 * L**3 * b * p * pr * r
        + 8 * L**3 * b * p * prr * r**2
        + 8 * L**3 * b * pr**2 * r**2
        + 8 * L**3 * e * p * prr
        + 8 * L**3 * e * pr**2
        + 2 * L**3 * pr**2 * r**2
        + 2 * L**3 * r**2 * xr**2
        - 8 * L * e * p * prr
        - 8 * L * e * pr**2
    )
    h_lambda_numerator = (
        -8 * K**2 * L**2 * e * p * pr * r**2
        + 16 * K * L**2 * e * p * pp * r
        - 8 * L**2 * MP**2 * r
        - 8 * L**2 * b * p**2 * r
        - 8 * L**2 * b * p * pr * r**2
        - 8 * L**2 * e * p * pr
        + 24 * e * p * pr
    )
    h_scale = 4 * L**5 * r**2
    hamiltonian_constant = -h0_numerator / h_scale
    hamiltonian_lambda_r = -h_lambda_numerator / h_scale

    m0_numerator = (
        4 * K**3 * L**2 * e * p * pr * r**2
        - 12 * K**2 * L**2 * e * p * pp * r
        - 2 * K**2 * L**2 * e * p * ppr * r**2
        - 2 * K**2 * L**2 * e * pp * pr * r**2
        + 6 * K * L**2 * MP**2 * r
        + 6 * K * L**2 * b * p**2 * r
        + 4 * K * L**2 * b * p * pr * r**2
        + 4 * K * L**2 * e * p * pr
        - 16 * K * e * p * pr
        - 2 * L**2 * b * p * ppr * r**2
        - 2 * L**2 * b * pp * pr * r**2
        - 2 * L**2 * e * p * ppr
        - 2 * L**2 * e * pp * pr
        - L**2 * pp * pr * r**2
        - L**2 * r**2 * xp * xr
        + 2 * e * p * ppr
        + 2 * e * pp * pr
    )
    m_k_numerator = (
        -4 * K * L**2 * e * p * pp * r**2
        + 2 * L**2 * MP**2 * r**2
        + 2 * L**2 * b * p**2 * r**2
        - 4 * e * p * pr * r
    )
    m_scale = L**2 * r**2
    return ConstraintCoefficients(
        hamiltonian_constant=hamiltonian_constant,
        hamiltonian_lambda_r=hamiltonian_lambda_r,
        momentum_constant=m0_numerator / m_scale,
        momentum_k_r=m_k_numerator / m_scale,
    )


def polar_areal_constraint_residuals(
    coefficients: ConstraintCoefficients,
    *,
    radial_metric_derivative: Any,
    angular_extrinsic_curvature_derivative: Any,
) -> tuple[Any, Any]:
    if not isinstance(coefficients, ConstraintCoefficients):
        raise TypeError("coefficients must be ConstraintCoefficients")
    return (
        coefficients.hamiltonian_constant
        + coefficients.hamiltonian_lambda_r * radial_metric_derivative,
        coefficients.momentum_constant
        + coefficients.momentum_k_r * angular_extrinsic_curvature_derivative,
    )


def polar_areal_constraint_rhs(
    radius: Real,
    state: Sequence[Real],
    parameters: InitialDataParameters,
    *,
    minimum_abs_jacobian_entry: Real = 1.0e-10,
) -> tuple[tuple[float, float], ConstraintCoefficients, dict[str, float]]:
    """Solve the two affine physical constraints for ``lambda_r`` and ``k_r``."""

    r = _finite("radius", radius, positive=True)
    L, K = _finite_pair("state", state)
    if L <= 0.0:
        raise InitialDataSolveStop(
            "negative_metric_factor",
            "radial metric left the positive polar-areal branch",
            {"radius": r, "radial_metric": L},
        )
    threshold = _finite(
        "minimum_abs_jacobian_entry", minimum_abs_jacobian_entry, positive=True
    )
    fields = compact_family_fields(r, parameters)
    action = parameters.action
    try:
        coefficients = polar_areal_constraint_coefficients(
            radius=r,
            radial_metric=L,
            angular_extrinsic_curvature=K,
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
    except (ArithmeticError, TypeError, ValueError) as exc:
        raise InitialDataSolveStop(
            "nonfinite_continuation",
            "FGC-QR constraint kernel could not be evaluated",
            {"radius": r, "error": str(exc)},
        ) from exc
    diagonal = (
        float(coefficients.hamiltonian_lambda_r),
        float(coefficients.momentum_k_r),
    )
    if not all(isfinite(value) for value in diagonal) or min(map(abs, diagonal)) <= threshold:
        raise InitialDataSolveStop(
            "singular_constraint_jacobian",
            "FGC-QR physical-constraint Jacobian lost its declared margin",
            {"radius": r, "diagonal": diagonal, "threshold": threshold},
        )
    derivative = (
        -float(coefficients.hamiltonian_constant) / diagonal[0],
        -float(coefficients.momentum_constant) / diagonal[1],
    )
    if not all(isfinite(value) for value in derivative):
        raise InitialDataSolveStop(
            "nonfinite_continuation",
            "FGC-QR constraint derivative became nonfinite",
            {"radius": r},
        )
    residual = polar_areal_constraint_residuals(
        coefficients,
        radial_metric_derivative=derivative[0],
        angular_extrinsic_curvature_derivative=derivative[1],
    )
    residual_infinity = max(abs(float(value)) for value in residual)
    determinant = diagonal[0] * diagonal[1]
    condition = max(map(abs, diagonal)) / min(map(abs, diagonal))
    return derivative, coefficients, {
        "residual_infinity": residual_infinity,
        "jacobian_determinant": determinant,
        "jacobian_condition_infinity": condition,
    }


def gauge_compatible_metric_time_derivatives(
    radius: Any,
    radial_metric: Any,
    radial_metric_derivative: Any,
) -> tuple[Any, Any]:
    """Solve ``C^t=C^r=0`` for the two free metric time derivatives.

    The expression is the exact flat-reference MHG result for the frozen
    ``tilde_normal_factor=4`` and the declared unit-lapse, zero-shift slice.
    """

    r = radius
    L = radial_metric
    Lr = radial_metric_derivative
    if r <= 0 or L <= 0:
        raise ValueError("gauge-compatible annular point requires r>0 and lambda>0")
    return 0 * L, (2 * L**3 - 2 * L + Lr * r) / (4 * L * r)


def misner_sharp_compactness(
    radius: Real,
    radial_metric: Real,
    angular_extrinsic_curvature: Real,
) -> float:
    r = _finite("radius", radius, positive=True)
    L = _finite("radial_metric", radial_metric, positive=True)
    K = _finite("angular_extrinsic_curvature", angular_extrinsic_curvature)
    answer = 1.0 + r**2 * K**2 - 1.0 / L**2
    if not isfinite(answer):
        raise InitialDataSolveStop(
            "nonfinite_continuation", "Misner-Sharp compactness became nonfinite"
        )
    return answer


def _rk4_step(
    rhs: Callable[[float, Sequence[float]], np.ndarray],
    radius: float,
    state: np.ndarray,
    step: float,
) -> np.ndarray:
    first = rhs(radius, state)
    second = rhs(radius + step / 2.0, state + step * first / 2.0)
    third = rhs(radius + step / 2.0, state + step * second / 2.0)
    fourth = rhs(radius + step, state + step * third)
    return state + step * (first + 2.0 * second + 2.0 * third + fourth) / 6.0


def _ssprk3_step(
    rhs: Callable[[float, Sequence[float]], np.ndarray],
    radius: float,
    state: np.ndarray,
    step: float,
) -> np.ndarray:
    first = state + step * rhs(radius, state)
    second = 0.75 * state + 0.25 * (first + step * rhs(radius + step, first))
    return state / 3.0 + 2.0 * (
        second + step * rhs(radius + step / 2.0, second)
    ) / 3.0


def _vacuum_denominator_minimum(
    support_radius: float,
    outer_radius: float,
    mass: float,
    momentum_constant: float,
) -> float:
    candidates = [support_radius, outer_radius]
    if mass > 0.0 and momentum_constant != 0.0:
        critical = (2.0 * momentum_constant**2 / mass) ** (1.0 / 3.0)
        if support_radius < critical < outer_radius:
            candidates.append(critical)
    values = [
        1.0 - 2.0 * mass / radius + momentum_constant**2 / radius**4
        for radius in candidates
    ]
    return min(values)


def solve_initial_data(
    parameters: InitialDataParameters,
    *,
    step_count: int,
    method: str = "RK4",
    minimum_abs_jacobian_entry: Real = 1.0e-10,
    maximum_constraint_residual: Real = 1.0e-10,
) -> InitialDataSolution:
    """Construct one regular-centre, finite-mass FGC-QR constraint slice."""

    if not isinstance(parameters, InitialDataParameters):
        raise TypeError("parameters must be InitialDataParameters")
    if isinstance(step_count, bool) or not isinstance(step_count, int) or step_count <= 0:
        raise ValueError("step_count must be a positive integer")
    steppers = {"RK4": _rk4_step, "SSPRK3": _ssprk3_step}
    if method not in steppers:
        raise ValueError("unknown ID1 integration method")
    residual_limit = _finite(
        "maximum_constraint_residual", maximum_constraint_residual, positive=True
    )
    start = parameters.support_minimum
    end = parameters.support_maximum
    step = (end - start) / step_count
    state = np.asarray((1.0, 0.0), dtype=np.float64)

    def rhs(radius: float, values: Sequence[float]) -> np.ndarray:
        derivative, _coefficients, _diagnostics = polar_areal_constraint_rhs(
            radius,
            values,
            parameters,
            minimum_abs_jacobian_entry=minimum_abs_jacobian_entry,
        )
        return np.asarray(derivative, dtype=np.float64)

    raw_points: list[tuple[float, float, float]] = [(start, state[0], state[1])]
    radius = start
    for _ in range(step_count):
        try:
            state = steppers[method](rhs, radius, state, step)
        except InitialDataSolveStop:
            raise
        except (ArithmeticError, FloatingPointError, ValueError) as exc:
            raise InitialDataSolveStop(
                "nonfinite_continuation",
                "FGC-QR radial constraint continuation failed",
                {"radius": radius, "error": str(exc)},
            ) from exc
        radius += step
        if state[0] <= 0.0:
            raise InitialDataSolveStop(
                "negative_metric_factor",
                "radial metric became nonpositive during continuation",
                {"radius": radius, "radial_metric": float(state[0])},
            )
        if not np.all(np.isfinite(state)):
            raise InitialDataSolveStop(
                "nonfinite_continuation",
                "constraint continuation produced a nonfinite state",
                {"radius": radius},
            )
        raw_points.append((radius, float(state[0]), float(state[1])))

    points: list[ConstraintPoint] = []
    minimum_determinant = float("inf")
    maximum_condition = 0.0
    maximum_residual = 0.0
    minimum_effective_planck = float("inf")
    for radius, radial_metric, angular_k in raw_points:
        derivative, _coefficients, diagnostics = polar_areal_constraint_rhs(
            radius,
            (radial_metric, angular_k),
            parameters,
            minimum_abs_jacobian_entry=minimum_abs_jacobian_entry,
        )
        compactness = misner_sharp_compactness(radius, radial_metric, angular_k)
        mass = radius * compactness / 2.0
        gauge_dt = gauge_compatible_metric_time_derivatives(
            radius, radial_metric, derivative[0]
        )
        fields = compact_family_fields(radius, parameters)
        effective_planck = (
            parameters.action.planck_mass**2
            + parameters.action.beta * fields["phi"] ** 2
        )
        minimum_effective_planck = min(minimum_effective_planck, effective_planck)
        minimum_determinant = min(
            minimum_determinant, abs(diagnostics["jacobian_determinant"])
        )
        maximum_condition = max(
            maximum_condition, diagnostics["jacobian_condition_infinity"]
        )
        maximum_residual = max(maximum_residual, diagnostics["residual_infinity"])
        points.append(
            ConstraintPoint(
                radius=radius,
                radial_metric=radial_metric,
                angular_extrinsic_curvature=angular_k,
                radial_metric_derivative=derivative[0],
                angular_extrinsic_curvature_derivative=derivative[1],
                hamiltonian_residual=0.0
                if diagnostics["residual_infinity"] == 0.0
                else float(
                    polar_areal_constraint_residuals(
                        _coefficients,
                        radial_metric_derivative=derivative[0],
                        angular_extrinsic_curvature_derivative=derivative[1],
                    )[0]
                ),
                momentum_residual=0.0
                if diagnostics["residual_infinity"] == 0.0
                else float(
                    polar_areal_constraint_residuals(
                        _coefficients,
                        radial_metric_derivative=derivative[0],
                        angular_extrinsic_curvature_derivative=derivative[1],
                    )[1]
                ),
                jacobian_determinant=diagnostics["jacobian_determinant"],
                jacobian_condition_infinity=diagnostics[
                    "jacobian_condition_infinity"
                ],
                compactness=compactness,
                misner_sharp_mass=mass,
                h_tt_time_derivative=float(gauge_dt[0]),
                h_tr_time_derivative=float(gauge_dt[1]),
            )
        )
    if maximum_residual > residual_limit:
        raise InitialDataSolveStop(
            "constraint_residual_loss",
            "solved physical-constraint residual exceeds its tolerance",
            {"observed": maximum_residual, "limit": residual_limit},
        )

    support_point = points[-1]
    mass = support_point.misner_sharp_mass
    momentum_constant = (
        support_point.angular_extrinsic_curvature * support_point.radius**3
    )
    vacuum_minimum = _vacuum_denominator_minimum(
        support_point.radius,
        parameters.outer_radius,
        mass,
        momentum_constant,
    )
    if not isfinite(vacuum_minimum) or vacuum_minimum <= 0.0:
        raise InitialDataSolveStop(
            "lost_asymptotic_condition",
            "analytic exterior vacuum continuation lost its positive metric factor",
            {"minimum_denominator": vacuum_minimum},
        )
    outer_k = momentum_constant / parameters.outer_radius**3
    outer_denominator = (
        1.0
        - 2.0 * mass / parameters.outer_radius
        + momentum_constant**2 / parameters.outer_radius**4
    )
    outer_lambda = 1.0 / sqrt(outer_denominator)
    peak = max(points, key=lambda item: item.compactness)
    inside_window = 0.1 <= peak.compactness <= 0.75
    return InitialDataSolution(
        method=method,
        step_count=step_count,
        parameters=parameters,
        points=tuple(points),
        peak_compactness=peak.compactness,
        peak_radius=peak.radius,
        outer_mass=mass,
        vacuum_momentum_constant=momentum_constant,
        outer_radial_metric=outer_lambda,
        outer_angular_extrinsic_curvature=outer_k,
        minimum_abs_constraint_jacobian_determinant=minimum_determinant,
        maximum_constraint_jacobian_condition_infinity=maximum_condition,
        maximum_constraint_residual_infinity=maximum_residual,
        minimum_effective_planck_coefficient=minimum_effective_planck,
        minimum_vacuum_metric_denominator=vacuum_minimum,
        finite_mass=isfinite(mass),
        regular_center=True,
        exact_inner_vacuum_buffer=True,
        exact_outer_vacuum_buffer=True,
        no_initial_trapped_sphere=peak.compactness < 1.0,
        compactness_inside_protocol_window=inside_window,
    )


def exact_constraint_and_gauge_crosscheck(
    *,
    radius: int | Fraction,
    radial_metric: int | Fraction,
    angular_extrinsic_curvature: int | Fraction,
    radial_metric_derivative: int | Fraction,
    angular_extrinsic_curvature_derivative: int | Fraction,
    phi: int | Fraction,
    phi_r: int | Fraction,
    phi_rr: int | Fraction,
    phi_pi: int | Fraction,
    phi_pi_r: int | Fraction,
    chi: int | Fraction,
    chi_r: int | Fraction,
    chi_rr: int | Fraction,
    chi_pi: int | Fraction,
    chi_pi_r: int | Fraction,
) -> dict[str, Any]:
    """Cross-check the fast kernel and gauge solve against the full evaluator."""

    raw = {
        name: value if isinstance(value, Fraction) else Q(value)
        for name, value in locals().items()
    }
    r = raw["radius"]
    L = raw["radial_metric"]
    K = raw["angular_extrinsic_curvature"]
    Lr = raw["radial_metric_derivative"]
    Kr = raw["angular_extrinsic_curvature_derivative"]
    if r <= 0 or L <= 0:
        raise ValueError("exact cross-check requires r>0 and lambda>0")
    h_tt_dt, h_tr_dt = gauge_compatible_metric_time_derivatives(r, L, Lr)
    state = SphericalState(
        h_tt=Jet2(-1, dt=h_tt_dt),
        h_tr=Jet2(0, dt=h_tr_dt),
        h_rr=Jet2(
            L**2,
            dt=4 * L**2 * K,
            dr=2 * L * Lr,
            dtr=8 * L * Lr * K + 4 * L**2 * Kr,
            drr=2 * Lr**2,
        ),
        areal_radius=Jet2(r, dt=-r * K, dr=1, dtr=-K - r * Kr),
        phi=Jet2(
            raw["phi"],
            dt=raw["phi_pi"],
            dr=raw["phi_r"],
            dtr=raw["phi_pi_r"],
            drr=raw["phi_rr"],
        ),
        chi=Jet2(
            raw["chi"],
            dt=raw["chi_pi"],
            dr=raw["chi_r"],
            dtr=raw["chi_pi_r"],
            drr=raw["chi_rr"],
        ),
        planck_mass=2,
        beta=-Q(1, 4),
        mu=3,
        g4=Q(1, 2),
        eta=Q(1, 2),
        branch="FGC-QR",
    )
    coefficients = polar_areal_constraint_coefficients(
        radius=r,
        radial_metric=L,
        angular_extrinsic_curvature=K,
        phi=raw["phi"],
        phi_r=raw["phi_r"],
        phi_rr=raw["phi_rr"],
        phi_pi=raw["phi_pi"],
        phi_pi_r=raw["phi_pi_r"],
        chi_r=raw["chi_r"],
        chi_pi=raw["chi_pi"],
        planck_mass=Q(2),
        beta=-Q(1, 4),
        scalar_mass=Q(3),
        quartic_coupling=Q(1, 2),
        eta=Q(1, 2),
    )
    specialized = polar_areal_constraint_residuals(
        coefficients,
        radial_metric_derivative=Lr,
        angular_extrinsic_curvature_derivative=Kr,
    )
    projection = physical_constraint_projections(state)
    full = (projection["H"], projection["M"])
    if specialized != full:
        raise ValueError("polar-areal constraint kernel differs from ACT1/VAR1")
    reference = flat_spherical_annulus_reference(radial_domain_minimum=r / 2)
    gauge = modified_harmonic_gauge_constraint(
        state,
        reference=reference,
        coordinate_radius=r,
        tilde_normal_factor=Q(4),
    )
    if gauge.constraint_up[:2] != (Q(0), Q(0)):
        raise ValueError("polar-areal gauge solve differs from complete REF1 gauge")
    return {
        "specialized_constraints": specialized,
        "full_unredefined_constraints": full,
        "exact_constraint_pair_equality": True,
        "gauge_constraint_up_t_r": gauge.constraint_up[:2],
        "exact_metric_defined_gauge_zero": True,
        "constraint_jacobian_diagonal": (
            coefficients.hamiltonian_lambda_r,
            coefficients.momentum_k_r,
        ),
        "source_is_unredefined_ACT1_VAR1": projection[
            "source_is_unredefined_residual"
        ],
    }
