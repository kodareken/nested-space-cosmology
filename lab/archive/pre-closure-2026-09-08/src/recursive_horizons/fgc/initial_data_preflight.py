"""GR-0 initial-data preflight for the frozen SF1 pulse normalization.

The first SF1 protocol declared a compact pulse of ``r*chi`` and a target
initial Misner--Sharp compactness interval.  Those two choices must be checked
against the constraints before any FGC-QR holdout can be opened.

This module specializes the unredefined GR-0 Hamiltonian and momentum
constraints to a regular maximal polar-areal slice,

``alpha=1, shift=0, R=r, K^r_r=-2*k, K^theta_theta=K^phi_phi=k``.

The resulting two first-order equations are integrated independently with
RK4 and SSPRK3.  Exact rational point controls in the reproducer compare the
specialized formula with the full ACT1/VAR1 evaluator.  The module also gives
a conservative analytic compactness bound for the entire original amplitude
box.  It is a protocol preflight, not an FGC-QR mechanism result and not the
later constraint-compatible ID1 family.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import exp, isfinite, log
from numbers import Real
from typing import Any, Callable, Sequence

from .modified_harmonic_constraints import physical_constraint_projections
from .spherical_reduction import Jet2, SphericalState


Q = Fraction


def _finite(name: str, value: object, *, positive: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a finite real scalar")
    result = float(value)
    if not isfinite(result):
        raise ValueError(f"{name} must be finite")
    if positive and result <= 0.0:
        raise ValueError(f"{name} must be positive")
    return result


@dataclass(frozen=True, slots=True)
class PulseParameters:
    chi_amplitude: float
    center: float = 12.0
    half_width: float = 2.0
    phi_amplitude: float = 1.0 / 131072.0
    planck_mass: float = 2.0
    scalar_mass: float = 3.0
    quartic_coupling: float = 0.5

    def __post_init__(self) -> None:
        for name in (
            "chi_amplitude",
            "center",
            "half_width",
            "phi_amplitude",
            "planck_mass",
            "scalar_mass",
            "quartic_coupling",
        ):
            object.__setattr__(self, name, _finite(name, getattr(self, name), positive=True))
        if self.center <= self.half_width:
            raise ValueError("pulse support must retain a strict centre buffer")

    @property
    def support_minimum(self) -> float:
        return self.center - self.half_width

    @property
    def support_maximum(self) -> float:
        return self.center + self.half_width

    @property
    def planck_mass_squared(self) -> float:
        return self.planck_mass**2


@dataclass(frozen=True, slots=True)
class SlicePoint:
    radius: float
    radial_metric: float
    angular_extrinsic_curvature: float
    compactness: float
    misner_sharp_mass: float


@dataclass(frozen=True, slots=True)
class SliceSolution:
    method: str
    step_count: int
    parameters: PulseParameters
    points: tuple[SlicePoint, ...]
    peak_compactness: float
    peak_radius: float
    outer_mass: float
    outer_k_times_r_cubed: float
    finite_mass: bool
    no_initial_trapped_sphere: bool


def compact_bump_with_derivatives(
    radius: Real,
    *,
    center: Real,
    half_width: Real,
) -> tuple[float, float, float]:
    """Return the declared compact ``C-infinity`` bump and two derivatives."""

    r = _finite("radius", radius)
    c = _finite("center", center)
    width = _finite("half_width", half_width, positive=True)
    x = (r - c) / width
    if abs(x) >= 1.0:
        return 0.0, 0.0, 0.0
    denominator = 1.0 - x * x
    exponent = 1.0 - 1.0 / denominator
    first_exponent = -2.0 * x / denominator**2
    second_exponent = -2.0 / denominator**2 - 8.0 * x * x / denominator**3
    value = exp(exponent)
    first = value * first_exponent / width
    second = value * (second_exponent + first_exponent**2) / width**2
    return value, first, second


def pulse_fields(radius: Real, parameters: PulseParameters) -> dict[str, float]:
    """Evaluate the frozen ``r*chi`` pulse and the zero-momentum ``phi`` seed."""

    r = _finite("radius", radius, positive=True)
    bump, bump_r, bump_rr = compact_bump_with_derivatives(
        r,
        center=parameters.center,
        half_width=parameters.half_width,
    )
    amplitude = parameters.chi_amplitude
    chi = amplitude * bump / r
    chi_r = amplitude * (bump_r / r - bump / r**2)
    chi_rr = amplitude * (bump_rr / r - 2.0 * bump_r / r**2 + 2.0 * bump / r**3)
    # The frozen future-ingoing convention is d_t(r*chi)=d_r(r*chi).
    chi_pi = amplitude * bump_r / r
    chi_pi_r = amplitude * (bump_rr / r - bump_r / r**2)
    phi = parameters.phi_amplitude * bump
    phi_r = parameters.phi_amplitude * bump_r
    phi_rr = parameters.phi_amplitude * bump_rr
    return {
        "chi": chi,
        "chi_r": chi_r,
        "chi_rr": chi_rr,
        "chi_pi": chi_pi,
        "chi_pi_r": chi_pi_r,
        "phi": phi,
        "phi_r": phi_r,
        "phi_rr": phi_rr,
        "phi_pi": 0.0,
        "phi_pi_r": 0.0,
    }


def gr0_energy_and_momentum_density(
    radius: Real,
    radial_metric: Real,
    parameters: PulseParameters,
) -> tuple[float, float]:
    r = _finite("radius", radius, positive=True)
    lam = _finite("radial_metric", radial_metric, positive=True)
    fields = pulse_fields(r, parameters)
    phi = fields["phi"]
    potential = (
        0.5 * parameters.scalar_mass**2 * phi**2
        + 0.25 * parameters.quartic_coupling * phi**4
    )
    rho = 0.5 * (
        fields["phi_pi"] ** 2
        + fields["chi_pi"] ** 2
        + (fields["phi_r"] ** 2 + fields["chi_r"] ** 2) / lam**2
    ) + potential
    momentum = (
        fields["phi_pi"] * fields["phi_r"]
        + fields["chi_pi"] * fields["chi_r"]
    )
    return rho, momentum


def gr0_constraint_residuals(
    *,
    radius: int | Fraction,
    radial_metric: int | Fraction,
    angular_extrinsic_curvature: int | Fraction,
    radial_metric_derivative: int | Fraction,
    angular_extrinsic_curvature_derivative: int | Fraction,
    phi: int | Fraction,
    phi_r: int | Fraction,
    phi_pi: int | Fraction,
    chi: int | Fraction,
    chi_r: int | Fraction,
    chi_pi: int | Fraction,
    planck_mass: int | Fraction = 2,
    scalar_mass: int | Fraction = 3,
    quartic_coupling: int | Fraction = Q(1, 2),
) -> tuple[Fraction, Fraction]:
    """Exact specialized GR-0 Hamiltonian and momentum residuals."""

    values = tuple(
        value if isinstance(value, Fraction) else Q(value)
        for value in (
            radius,
            radial_metric,
            angular_extrinsic_curvature,
            radial_metric_derivative,
            angular_extrinsic_curvature_derivative,
            phi,
            phi_r,
            phi_pi,
            chi,
            chi_r,
            chi_pi,
            planck_mass,
            scalar_mass,
            quartic_coupling,
        )
    )
    (
        r,
        lam,
        k,
        lam_r,
        k_r,
        phi_value,
        phi_gradient,
        phi_momentum,
        _chi_value,
        chi_gradient,
        chi_momentum,
        m_planck,
        mu,
        g4,
    ) = values
    if r <= 0 or lam <= 0 or m_planck <= 0 or mu <= 0 or g4 <= 0:
        raise ValueError("exact GR-0 constraint point lies outside its positive domain")
    potential = mu**2 * phi_value**2 / 2 + g4 * phi_value**4 / 4
    rho = (
        phi_momentum**2
        + chi_momentum**2
        + (phi_gradient**2 + chi_gradient**2) / lam**2
    ) / 2 + potential
    momentum_density = phi_momentum * phi_gradient + chi_momentum * chi_gradient
    planck_squared = m_planck**2
    hamiltonian = planck_squared * (
        (1 - 1 / lam**2) / r**2
        + 2 * lam_r / (r * lam**3)
        - 3 * k**2
    ) - rho
    momentum = 2 * planck_squared * (k_r + 3 * k / r) - momentum_density
    return hamiltonian, momentum


def exact_gr0_constraint_crosscheck(
    *,
    radius: int | Fraction,
    radial_metric: int | Fraction,
    angular_extrinsic_curvature: int | Fraction,
    radial_metric_derivative: int | Fraction,
    angular_extrinsic_curvature_derivative: int | Fraction,
    phi: int | Fraction,
    phi_r: int | Fraction,
    phi_pi: int | Fraction,
    chi: int | Fraction,
    chi_r: int | Fraction,
    chi_pi: int | Fraction,
    planck_mass: int | Fraction = 2,
    scalar_mass: int | Fraction = 3,
    quartic_coupling: int | Fraction = Q(1, 2),
) -> dict[str, Any]:
    """Compare the specialized constraints with the full unredefined evaluator.

    The exact ADM data use unit lapse, zero shift, areal radius ``R=r``, and
    the maximal extrinsic-curvature eigenvalues

    ``K^r_r=-2*k`` and ``K^theta_theta=K^phi_phi=k``.

    With the repository's sign convention this fixes
    ``lambda_t=2*lambda*k`` and ``R_t=-r*k``.  Their radial derivatives are
    included explicitly so the full four-dimensional evaluator sees the same
    momentum data as the specialized constraint formula.  Coordinate-time
    accelerations are set to zero; CON4 independently proves that the physical
    normal projections do not depend on them.
    """

    names = (
        "radius",
        "radial_metric",
        "angular_extrinsic_curvature",
        "radial_metric_derivative",
        "angular_extrinsic_curvature_derivative",
        "phi",
        "phi_r",
        "phi_pi",
        "chi",
        "chi_r",
        "chi_pi",
        "planck_mass",
        "scalar_mass",
        "quartic_coupling",
    )
    raw_values = (
        radius,
        radial_metric,
        angular_extrinsic_curvature,
        radial_metric_derivative,
        angular_extrinsic_curvature_derivative,
        phi,
        phi_r,
        phi_pi,
        chi,
        chi_r,
        chi_pi,
        planck_mass,
        scalar_mass,
        quartic_coupling,
    )
    values: dict[str, Fraction] = {}
    for name, value in zip(names, raw_values, strict=True):
        if isinstance(value, bool) or not isinstance(value, (int, Fraction)):
            raise TypeError(f"{name} must be an exact rational")
        values[name] = value if isinstance(value, Fraction) else Q(value)
    r = values["radius"]
    lam = values["radial_metric"]
    k = values["angular_extrinsic_curvature"]
    lam_r = values["radial_metric_derivative"]
    k_r = values["angular_extrinsic_curvature_derivative"]
    if r <= 0 or lam <= 0:
        raise ValueError("exact GR-0 cross-check requires positive r and lambda")

    state = SphericalState(
        h_tt=Jet2.constant(-1),
        h_tr=Jet2.constant(0),
        h_rr=Jet2(
            lam**2,
            dt=4 * lam**2 * k,
            dr=2 * lam * lam_r,
            dtr=8 * lam * lam_r * k + 4 * lam**2 * k_r,
        ),
        areal_radius=Jet2(
            r,
            dt=-r * k,
            dr=Q(1),
            dtr=-k - r * k_r,
        ),
        phi=Jet2(
            values["phi"],
            dt=values["phi_pi"],
            dr=values["phi_r"],
        ),
        chi=Jet2(
            values["chi"],
            dt=values["chi_pi"],
            dr=values["chi_r"],
        ),
        planck_mass=values["planck_mass"],
        mu=values["scalar_mass"],
        g4=values["quartic_coupling"],
        branch="GR-0",
    )
    specialized = gr0_constraint_residuals(
        radius=r,
        radial_metric=lam,
        angular_extrinsic_curvature=k,
        radial_metric_derivative=lam_r,
        angular_extrinsic_curvature_derivative=k_r,
        phi=values["phi"],
        phi_r=values["phi_r"],
        phi_pi=values["phi_pi"],
        chi=values["chi"],
        chi_r=values["chi_r"],
        chi_pi=values["chi_pi"],
        planck_mass=values["planck_mass"],
        scalar_mass=values["scalar_mass"],
        quartic_coupling=values["quartic_coupling"],
    )
    projection = physical_constraint_projections(state)
    full = (projection["H"], projection["M"])
    if full != specialized:
        raise ValueError("specialized GR-0 constraints differ from unredefined ACT1/VAR1")
    return {
        "state": state,
        "specialized_H": specialized[0],
        "specialized_M": specialized[1],
        "full_unredefined_H": full[0],
        "full_unredefined_M": full[1],
        "exact_pair_equality": True,
        "source_is_unredefined_ACT1_VAR1": projection["source_is_unredefined_residual"],
        "coordinate_time_accelerations_set_to_zero_only_for_control": True,
    }


def gr0_constraint_rhs(
    radius: Real,
    state: Sequence[Real],
    parameters: PulseParameters,
) -> tuple[float, float]:
    """Solve the two specialized GR-0 constraints for ``lambda_r`` and ``k_r``."""

    if isinstance(state, (str, bytes)) or len(state) != 2:
        raise ValueError("GR-0 initial-data state must contain lambda and k")
    r = _finite("radius", radius, positive=True)
    lam = _finite("radial_metric", state[0], positive=True)
    k = _finite("angular_extrinsic_curvature", state[1])
    rho, momentum = gr0_energy_and_momentum_density(r, lam, parameters)
    planck_squared = parameters.planck_mass_squared
    lam_r = r * lam**3 / 2.0 * (
        rho / planck_squared
        - (1.0 - 1.0 / lam**2) / r**2
        + 3.0 * k**2
    )
    k_r = momentum / (2.0 * planck_squared) - 3.0 * k / r
    if not isfinite(lam_r) or not isfinite(k_r):
        raise ArithmeticError("GR-0 constraint right-hand side became nonfinite")
    return lam_r, k_r


def misner_sharp_compactness(
    radius: Real,
    radial_metric: Real,
    angular_extrinsic_curvature: Real,
) -> float:
    r = _finite("radius", radius, positive=True)
    lam = _finite("radial_metric", radial_metric, positive=True)
    k = _finite("angular_extrinsic_curvature", angular_extrinsic_curvature)
    value = 1.0 + r**2 * k**2 - 1.0 / lam**2
    if not isfinite(value):
        raise ArithmeticError("Misner--Sharp compactness became nonfinite")
    return value


def _rk4_step(
    rhs: Callable[[float, Sequence[float]], tuple[float, float]],
    radius: float,
    state: tuple[float, float],
    step: float,
) -> tuple[float, float]:
    def shifted(base: Sequence[float], tangent: Sequence[float], scale: float) -> tuple[float, float]:
        return base[0] + scale * tangent[0], base[1] + scale * tangent[1]

    first = rhs(radius, state)
    second = rhs(radius + step / 2.0, shifted(state, first, step / 2.0))
    third = rhs(radius + step / 2.0, shifted(state, second, step / 2.0))
    fourth = rhs(radius + step, shifted(state, third, step))
    return (
        state[0] + step * (first[0] + 2.0 * second[0] + 2.0 * third[0] + fourth[0]) / 6.0,
        state[1] + step * (first[1] + 2.0 * second[1] + 2.0 * third[1] + fourth[1]) / 6.0,
    )


def _ssprk3_step(
    rhs: Callable[[float, Sequence[float]], tuple[float, float]],
    radius: float,
    state: tuple[float, float],
    step: float,
) -> tuple[float, float]:
    first_rhs = rhs(radius, state)
    first = tuple(state[index] + step * first_rhs[index] for index in range(2))
    second_rhs = rhs(radius + step, first)
    second = tuple(
        0.75 * state[index]
        + 0.25 * (first[index] + step * second_rhs[index])
        for index in range(2)
    )
    third_rhs = rhs(radius + step / 2.0, second)
    return tuple(
        state[index] / 3.0
        + 2.0 * (second[index] + step * third_rhs[index]) / 3.0
        for index in range(2)
    )


def solve_gr0_initial_slice(
    parameters: PulseParameters,
    *,
    step_count: int,
    method: str = "RK4",
) -> SliceSolution:
    """Integrate the regular-centre GR-0 constraint slice across the pulse."""

    if isinstance(step_count, bool) or not isinstance(step_count, int) or step_count <= 0:
        raise ValueError("step_count must be a positive integer")
    steppers = {"RK4": _rk4_step, "SSPRK3": _ssprk3_step}
    if method not in steppers:
        raise ValueError("unknown GR-0 preflight integration method")
    start = parameters.support_minimum
    end = parameters.support_maximum
    step = (end - start) / step_count
    rhs = lambda radius, state: gr0_constraint_rhs(radius, state, parameters)
    state = (1.0, 0.0)
    points: list[SlicePoint] = []

    def point(radius: float, values: tuple[float, float]) -> SlicePoint:
        compactness = misner_sharp_compactness(radius, values[0], values[1])
        return SlicePoint(
            radius=radius,
            radial_metric=values[0],
            angular_extrinsic_curvature=values[1],
            compactness=compactness,
            misner_sharp_mass=radius * compactness / 2.0,
        )

    points.append(point(start, state))
    radius = start
    for _ in range(step_count):
        state = steppers[method](rhs, radius, state, step)
        radius += step
        if state[0] <= 0.0 or not all(isfinite(value) for value in state):
            raise ArithmeticError("GR-0 initial-data integration left its regular branch")
        points.append(point(radius, state))

    peak = max(points, key=lambda item: item.compactness)
    outer = points[-1]
    outer_k_r3 = outer.angular_extrinsic_curvature * outer.radius**3
    return SliceSolution(
        method=method,
        step_count=step_count,
        parameters=parameters,
        points=tuple(points),
        peak_compactness=peak.compactness,
        peak_radius=peak.radius,
        outer_mass=outer.misner_sharp_mass,
        outer_k_times_r_cubed=outer_k_r3,
        finite_mass=isfinite(outer.misner_sharp_mass),
        no_initial_trapped_sphere=peak.compactness < 1.0,
    )


def observed_convergence_order(coarse: float, medium: float, fine: float) -> float | None:
    numerator = abs(coarse - medium)
    denominator = abs(medium - fine)
    if denominator == 0.0:
        return None if numerator else float("inf")
    if numerator == 0.0:
        return 0.0
    return log(numerator / denominator, 2.0)


def proto1_compactness_upper_bound() -> dict[str, Fraction | bool | str]:
    """Conservative exact bootstrap bound for the complete PROTO1 amplitude box.

    The declared bump obeys ``B<=1``.  Writing
    ``y=1/(1-x^2)>=1`` gives
    ``|B_x| <= 2*y^2*exp(1-y) <= 8/e < 3``; at half-width two,
    ``|B_r|<3/2``.  The calculation assumes ``lambda^-2<=4`` in the
    matter-energy bound and returns a strictly stronger value, closing the
    bootstrap.
    """

    amplitude = Q(1, 8)
    phi_amplitude = Q(1, 131072)
    support_minimum = Q(10)
    support_maximum = Q(14)
    planck_squared = Q(4)
    inverse_metric_bootstrap = Q(4)
    bump_r_bound = Q(3, 2)
    chi_pi_bound = amplitude * bump_r_bound / support_minimum
    chi_r_bound = amplitude * (
        bump_r_bound / support_minimum + 1 / support_minimum**2
    )
    momentum_density_bound = chi_pi_bound * chi_r_bound
    phi_r_bound = phi_amplitude * bump_r_bound
    potential_bound = Q(9, 2) * phi_amplitude**2 + Q(1, 8) * phi_amplitude**4
    energy_density_bound = (
        chi_pi_bound**2
        + inverse_metric_bootstrap * chi_r_bound**2
        + inverse_metric_bootstrap * phi_r_bound**2
    ) / 2 + potential_bound
    k_r3_bound = (
        momentum_density_bound
        * (support_maximum**4 - support_minimum**4)
        / (8 * planck_squared)
    )
    integrated_energy_bound = (
        energy_density_bound
        * (support_maximum**3 - support_minimum**3)
        / 3
    )
    integrated_flux_bound = (
        k_r3_bound
        * momentum_density_bound
        * (support_maximum - support_minimum)
    )
    mass_absolute_bound = (
        integrated_energy_bound + integrated_flux_bound
    ) / (2 * planck_squared)
    compactness_absolute_bound = 2 * mass_absolute_bound / support_minimum
    r2_k2_bound = k_r3_bound**2 / support_minimum**4
    improved_inverse_metric_bound = 1 + r2_k2_bound + compactness_absolute_bound
    if not compactness_absolute_bound < Q(1, 64) < Q(1, 10):
        raise ValueError("PROTO1 analytic compactness separation did not close")
    if not improved_inverse_metric_bound < inverse_metric_bootstrap:
        raise ValueError("PROTO1 inverse-metric bootstrap did not close")
    return {
        "classification": "conditional_exact_maximal_polar_areal_GR0_compactness_bound",
        "bump_absolute_bound": Q(1),
        "bump_x_derivative_absolute_bound": Q(3),
        "bump_r_derivative_absolute_bound": bump_r_bound,
        "chi_pi_absolute_bound": chi_pi_bound,
        "chi_r_absolute_bound": chi_r_bound,
        "momentum_density_absolute_bound": momentum_density_bound,
        "energy_density_absolute_bound": energy_density_bound,
        "k_times_r_cubed_absolute_bound": k_r3_bound,
        "mass_absolute_bound": mass_absolute_bound,
        "compactness_absolute_bound": compactness_absolute_bound,
        "compactness_bound_below_one_over_64": True,
        "one_over_64_below_protocol_minimum_one_over_10": True,
        "assumed_inverse_radial_metric_squared_upper_bound": inverse_metric_bootstrap,
        "improved_inverse_radial_metric_squared_upper_bound": improved_inverse_metric_bound,
        "inverse_metric_bootstrap_closed": True,
        "regular_center_removes_vacuum_mass_and_r_minus_3_momentum_constants": True,
        "outer_vacuum_compactness_decreases_as_two_m_over_r": True,
    }
