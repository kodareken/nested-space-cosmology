"""Reduced homogeneous Einstein--Cartan spin-fluid collapse benchmark.

This module implements only the published closed-FLRW random-spin-fluid
reduction ``adot^2=-1+A/a^2-B/a^4``.  It is action-motivated, but is not a
full tetrad/Dirac evolution, a resolved vacuum boundary, or GMF-1B's global
parent-to-child solution.  In particular, it derives neither external ``Q``,
late dark energy, nor a varying causal speed.

Units are ``c=hbar=1`` and ``kappa=8*pi*G``.  Every public numerical API
rejects bools, non-real values, and non-finite values.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import cos, isfinite, pi, sin, sqrt, ulp
from numbers import Integral, Real
from typing import Literal


Branch = Literal["collapse", "expansion"]
Region = Literal["normal", "trapped", "anti_trapped", "marginal"]


def _real(name: str, value: Real) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a finite real number")
    result = float(value)
    if not isfinite(result):
        raise ValueError(f"{name} must be a finite real number")
    return result


def _positive(name: str, value: Real) -> float:
    result = _real(name, value)
    if result <= 0.0:
        raise ValueError(f"{name} must be positive")
    return result


def _nonnegative(name: str, value: Real) -> float:
    result = _real(name, value)
    if result < 0.0:
        raise ValueError(f"{name} must be non-negative")
    return result


def _finite_mapping(values: dict[str, object]) -> dict[str, object]:
    for value in values.values():
        if isinstance(value, float) and not isfinite(value):
            raise ValueError("derived quantity is outside the finite evaluation range")
    return values


@dataclass(frozen=True, slots=True)
class EinsteinCartanSpec:
    """Parameters of the closed-FLRW random-spin-fluid reduction.

    ``A`` has length squared and ``B`` length to the fourth power, so the
    scale factor has length.  ``A**2 > 4B`` gives two distinct positive
    turning points.  ``chi_boundary`` selects an outward-monotone finite cap;
    it is not a whole boundaryless closed ``S^3`` child.
    """

    A: float
    B: float
    chi_boundary: float = pi / 3.0
    gravitational_constant: float = 1.0

    def __post_init__(self) -> None:
        A = _positive("A", self.A)
        B = _positive("B", self.B)
        chi = _real("chi_boundary", self.chi_boundary)
        G = _positive("gravitational_constant", self.gravitational_constant)
        discriminant = A * A - 4.0 * B
        if not isfinite(discriminant):
            raise ValueError("A and B must give a finite turning-point discriminant")
        if discriminant <= 0.0:
            raise ValueError("A**2 must be strictly greater than 4*B")
        if not 0.0 < chi < pi / 2.0:
            raise ValueError("chi_boundary must satisfy 0 < chi_boundary < pi/2")
        object.__setattr__(self, "A", A)
        object.__setattr__(self, "B", B)
        object.__setattr__(self, "chi_boundary", chi)
        object.__setattr__(self, "gravitational_constant", G)

    @property
    def kappa(self) -> float:
        return 8.0 * pi * self.gravitational_constant


def turning_points(spec: EinsteinCartanSpec) -> dict[str, float]:
    """Return exact ``x=a^2`` and scale-factor turning points.

    The lower root is the reduced-model bounce and the upper root is the
    maximum closed-FLRW radius.  The allowed branch is ``a_min <= a <= a_max``.
    """

    discriminant = spec.A * spec.A - 4.0 * spec.B
    root = sqrt(discriminant)
    x_max = 0.5 * (spec.A + root)
    # The conjugate expression loses every significant digit when B << A^2.
    # Since x_min*x_max=B exactly in the polynomial, this is stable instead.
    x_min = spec.B / x_max
    return _finite_mapping({
        "discriminant": discriminant,
        "x_min": x_min,
        "x_max": x_max,
        "a_min": sqrt(x_min),
        "a_max": sqrt(x_max),
    })  # type: ignore[return-value]


def _allowed_scale_factor(spec: EinsteinCartanSpec, scale_factor: Real) -> float:
    a = _positive("scale_factor", scale_factor)
    roots = turning_points(spec)
    a_min = float(roots["a_min"])
    a_max = float(roots["a_max"])
    if a < a_min and abs(a - a_min) <= 32.0 * ulp(a_min):
        return a_min
    if a > a_max and abs(a - a_max) <= 32.0 * ulp(a_max):
        return a_max
    if a < a_min or a > a_max:
        raise ValueError("scale_factor is outside the real closed-FLRW branch")
    return a


def friedmann_F(spec: EinsteinCartanSpec, scale_factor: Real) -> float:
    """Return ``F(a)=adot^2=-1+A/a^2-B/a^4`` on the real branch."""

    a = _allowed_scale_factor(spec, scale_factor)
    roots = turning_points(spec)
    a_min = float(roots["a_min"])
    a_max = float(roots["a_max"])
    if a == a_min or a == a_max:
        return 0.0
    x = a * a
    result = (float(roots["x_max"]) - x) * (x - float(roots["x_min"])) / (x * x)
    if result < 0.0 or not isfinite(result):
        raise ValueError("Friedmann value is outside the real finite branch")
    return result


def scale_factor_acceleration(spec: EinsteinCartanSpec, scale_factor: Real) -> float:
    """Return ``addot=-A/a^3+2B/a^5`` from ``F'(a)/2``."""

    a = _allowed_scale_factor(spec, scale_factor)
    result = -spec.A / a**3 + 2.0 * spec.B / a**5
    if not isfinite(result):
        raise ValueError("scale-factor acceleration is outside the finite range")
    return result


def scale_factor_velocity(spec: EinsteinCartanSpec, scale_factor: Real, branch: Branch) -> float:
    """Return the signed velocity on the collapse or expansion branch."""

    if branch not in ("collapse", "expansion"):
        raise ValueError("branch must be 'collapse' or 'expansion'")
    magnitude = sqrt(friedmann_F(spec, scale_factor))
    return -magnitude if branch == "collapse" else magnitude


def fluid_state(spec: EinsteinCartanSpec, scale_factor: Real) -> dict[str, float]:
    """Return the total and component densities/pressures of the reduction.

    The positive term is radiation, ``w=1/3`` and scales as ``a^-4``.  The
    torsion term is negative, has ``w=1``, and scales as ``a^-6``.  These are
    a homogeneous effective-fluid closure; they are not dark components.
    """

    a = _allowed_scale_factor(spec, scale_factor)
    radiation_density = 3.0 * spec.A / (spec.kappa * a**4)
    radiation_pressure = radiation_density / 3.0
    torsion_density = -3.0 * spec.B / (spec.kappa * a**6)
    torsion_pressure = torsion_density
    density = radiation_density + torsion_density
    pressure = radiation_pressure + torsion_pressure
    return _finite_mapping({
        "rho_radiation": radiation_density,
        "p_radiation": radiation_pressure,
        "w_radiation": 1.0 / 3.0,
        "rho_torsion": torsion_density,
        "p_torsion": torsion_pressure,
        "w_torsion": 1.0,
        "rho_total": density,
        "p_total": pressure,
        "w_total": pressure / density,
        "rho_plus_p": (4.0 * spec.A / a**4 - 6.0 * spec.B / a**6) / spec.kappa,
        "rho_plus_3p": (6.0 * spec.A / a**4 - 12.0 * spec.B / a**6) / spec.kappa,
        "q_internal": 0.0,
    })  # type: ignore[return-value]


def continuity_residual(spec: EinsteinCartanSpec, scale_factor: Real, branch: Branch) -> float:
    """Return ``rhodot+3H(rho+p)`` for the total effective fluid."""

    a = _allowed_scale_factor(spec, scale_factor)
    velocity = scale_factor_velocity(spec, a, branch)
    state = fluid_state(spec, a)
    density_derivative = 3.0 * (-4.0 * spec.A / a**5 + 6.0 * spec.B / a**7) / spec.kappa
    hubble = velocity / a
    result = density_derivative * velocity + 3.0 * hubble * (
        float(state["rho_total"]) + float(state["p_total"])
    )
    return 0.0 if abs(result) <= 1.0e-13 else result


def curvature_invariants(spec: EinsteinCartanSpec, scale_factor: Real) -> dict[str, float]:
    """Return finite FLRW curvature invariants on the allowed branch."""

    a = _allowed_scale_factor(spec, scale_factor)
    acceleration = scale_factor_acceleration(spec, a)
    F = friedmann_F(spec, a)
    state = fluid_state(spec, a)
    ricci = 6.0 * spec.B / a**6
    ricci_squared = spec.kappa**2 * (
        float(state["rho_total"]) ** 2 + 3.0 * float(state["p_total"]) ** 2
    )
    kretschmann = 12.0 * ((acceleration / a) ** 2 + ((F + 1.0) / a**2) ** 2)
    return _finite_mapping({
        "ricci_scalar": ricci,
        "ricci_tensor_squared": ricci_squared,
        "kretschmann_scalar": kretschmann,
        "weyl_tensor_squared": 0.0,
    })  # type: ignore[return-value]


def null_expansions(
    spec: EinsteinCartanSpec, scale_factor: Real, chi: Real, branch: Branch
) -> dict[str, float]:
    """Return future radial expansions with ``k_plus dot k_minus=-2``."""

    a = _allowed_scale_factor(spec, scale_factor)
    chi = _real("chi", chi)
    if not 0.0 < chi < pi:
        raise ValueError("chi must satisfy 0 < chi < pi")
    velocity = scale_factor_velocity(spec, a, branch)
    hubble = velocity / a
    cotangent = cos(chi) / sin(chi)
    plus = 2.0 * (hubble + cotangent / a)
    minus = 2.0 * (hubble - cotangent / a)
    return _finite_mapping({
        "hubble": hubble,
        "theta_plus": plus,
        "theta_minus": minus,
        "normalization_k_plus_dot_k_minus": -2.0,
    })  # type: ignore[return-value]


def classify_round_sphere(
    spec: EinsteinCartanSpec,
    scale_factor: Real,
    chi: Real,
    branch: Branch,
    *,
    tolerance: Real = 1.0e-12,
) -> Region:
    """Classify a sphere as normal, trapped, anti-trapped, or marginal."""

    tolerance = _nonnegative("tolerance", tolerance)
    expansions = null_expansions(spec, scale_factor, chi, branch)
    plus = float(expansions["theta_plus"])
    minus = float(expansions["theta_minus"])
    if abs(plus) <= tolerance or abs(minus) <= tolerance:
        return "marginal"
    if plus < 0.0 and minus < 0.0:
        return "trapped"
    if plus > 0.0 and minus > 0.0:
        return "anti_trapped"
    return "normal"


def marginal_scale_factor_roots(spec: EinsteinCartanSpec, chi: Real) -> tuple[float, float]:
    """Return analytic marginal roots from ``F(a)=cot(chi)^2``.

    Both roots must lie on the real collapse branch.  A requested angle with
    no real marginal spheres raises ``ValueError`` rather than returning NaN.
    """

    chi = _real("chi", chi)
    if not 0.0 < chi < pi:
        raise ValueError("chi must satisfy 0 < chi < pi")
    cotangent_squared = (cos(chi) / sin(chi)) ** 2
    coefficient = 1.0 + cotangent_squared
    if not isfinite(coefficient):
        raise ValueError("chi is outside the finite marginal-root evaluation range")
    discriminant = spec.A * spec.A - 4.0 * spec.B * coefficient
    if not isfinite(discriminant) or discriminant <= 0.0:
        raise ValueError("no distinct real marginal roots exist for this chi")
    root = sqrt(discriminant)
    x_high = (spec.A + root) / (2.0 * coefficient)
    # The conjugate subtraction loses the small root when B*C << A**2.
    # The quadratic-root product is B/C, giving this stable equivalent.
    x_low = 2.0 * spec.B / (spec.A + root)
    if x_low <= 0.0 or x_high <= 0.0:
        raise ValueError("marginal roots are not positive")
    a_low, a_high = sqrt(x_low), sqrt(x_high)
    _allowed_scale_factor(spec, a_low)
    _allowed_scale_factor(spec, a_high)
    return (a_low, a_high)


def causal_sequence(spec: EinsteinCartanSpec) -> dict[str, object]:
    """Return a complete collapse--bounce--expansion sign sequence for the cap."""

    roots = turning_points(spec)
    a_min = float(roots["a_min"])
    a_max = float(roots["a_max"])
    marginal_low, marginal_high = marginal_scale_factor_roots(spec, spec.chi_boundary)
    probe = sqrt((marginal_low**2 + marginal_high**2) / 2.0)
    sequence = [
        {"label": "maximum_collapse_endpoint", "scale_factor": a_max, "branch": "collapse", "region": classify_round_sphere(spec, a_max, spec.chi_boundary, "collapse")},
        {"label": "collapse_outer_marginal", "scale_factor": marginal_high, "branch": "collapse", "region": classify_round_sphere(spec, marginal_high, spec.chi_boundary, "collapse", tolerance=1.0e-9)},
        {"label": "collapsing_trapped_probe", "scale_factor": probe, "branch": "collapse", "region": classify_round_sphere(spec, probe, spec.chi_boundary, "collapse")},
        {"label": "collapse_inner_marginal", "scale_factor": marginal_low, "branch": "collapse", "region": classify_round_sphere(spec, marginal_low, spec.chi_boundary, "collapse", tolerance=1.0e-9)},
        {"label": "bounce", "scale_factor": a_min, "branch": "expansion", "region": classify_round_sphere(spec, a_min, spec.chi_boundary, "expansion")},
        {"label": "expansion_inner_marginal", "scale_factor": marginal_low, "branch": "expansion", "region": classify_round_sphere(spec, marginal_low, spec.chi_boundary, "expansion", tolerance=1.0e-9)},
        {"label": "expanding_anti_trapped_probe", "scale_factor": probe, "branch": "expansion", "region": classify_round_sphere(spec, probe, spec.chi_boundary, "expansion")},
        {"label": "expansion_outer_marginal", "scale_factor": marginal_high, "branch": "expansion", "region": classify_round_sphere(spec, marginal_high, spec.chi_boundary, "expansion", tolerance=1.0e-9)},
        {"label": "maximum_expansion_endpoint", "scale_factor": a_max, "branch": "expansion", "region": classify_round_sphere(spec, a_max, spec.chi_boundary, "expansion")},
    ]
    expected = ["normal", "marginal", "trapped", "marginal", "normal", "marginal", "anti_trapped", "marginal", "normal"]
    return {
        "sequence": sequence,
        "marginal_scale_factors": [marginal_low, marginal_high],
        "single_interior_metric_causal_sequence": [item["region"] for item in sequence] == expected,
    }


def boundary_areal_radius(spec: EinsteinCartanSpec, scale_factor: Real) -> float:
    """Return the finite-cap boundary radius ``R_b=a sin(chi_b)``."""

    return _allowed_scale_factor(spec, scale_factor) * sin(spec.chi_boundary)


def boundary_misner_sharp_mass(spec: EinsteinCartanSpec, scale_factor: Real) -> float:
    """Return ``M_MS=sin^3(chi_b)(A/a-B/a^3)/(2G)`` at the cap boundary."""

    a = _allowed_scale_factor(spec, scale_factor)
    result = sin(spec.chi_boundary) ** 3 * (spec.A / a - spec.B / a**3) / (2.0 * spec.gravitational_constant)
    if not isfinite(result):
        raise ValueError("boundary Misner-Sharp mass is outside the finite range")
    return result


def boundary_work_residual(spec: EinsteinCartanSpec, scale_factor: Real, branch: Branch) -> float:
    """Return ``dM/dt + 4 pi p R^2 Rdot`` at the material cap boundary."""

    a = _allowed_scale_factor(spec, scale_factor)
    velocity = scale_factor_velocity(spec, a, branch)
    sine_cubed = sin(spec.chi_boundary) ** 3
    dM_da = sine_cubed * (-spec.A / a**2 + 3.0 * spec.B / a**4) / (2.0 * spec.gravitational_constant)
    state = fluid_state(spec, a)
    radius = boundary_areal_radius(spec, a)
    radius_velocity = velocity * sin(spec.chi_boundary)
    result = dM_da * velocity + 4.0 * pi * float(state["p_total"]) * radius**2 * radius_velocity
    return 0.0 if abs(result) <= 1.0e-13 else result


def static_parent_mass_drift(spec: EinsteinCartanSpec) -> dict[str, float | bool]:
    """Compare the upper-turning mass to the bounce mass.

    A fixed Schwarzschild/Kottler vacuum mass cannot equal both values when
    they differ.  The result therefore diagnoses the need for flux and/or a
    resolved layer; it does not construct either one.
    """

    roots = turning_points(spec)
    mass_at_maximum = boundary_misner_sharp_mass(spec, float(roots["a_max"]))
    mass_at_bounce = boundary_misner_sharp_mass(spec, float(roots["a_min"]))
    drift = mass_at_bounce - mass_at_maximum
    return _finite_mapping({
        "static_parent_mass_at_a_max": mass_at_maximum,
        "boundary_mass_at_bounce": mass_at_bounce,
        "boundary_mass_drift_to_bounce": drift,
        "bounce_to_initial_parent_mass_ratio": mass_at_bounce / mass_at_maximum,
        "fixed_vacuum_parent_can_remain_smooth_without_flux_or_layer": abs(drift) <= 1.0e-12,
    })  # type: ignore[return-value]


def half_cycle_proper_time(spec: EinsteinCartanSpec, panels: Integral = 512) -> float:
    """Integrate the bounce-to-maximum proper time by smooth Simpson quadrature.

    With ``x=a^2=x_min+(x_max-x_min)sin^2(theta)``, the endpoint singularity
    cancels exactly and ``dt=sqrt(x(theta)) dtheta`` on ``[0,pi/2]``.
    """

    if isinstance(panels, bool) or not isinstance(panels, Integral):
        raise ValueError("panels must be a positive even integer")
    panels = int(panels)
    if panels <= 0 or panels % 2:
        raise ValueError("panels must be a positive even integer")
    roots = turning_points(spec)
    x_min = float(roots["x_min"])
    difference = float(roots["x_max"]) - x_min
    step = (pi / 2.0) / panels

    def integrand(theta: float) -> float:
        return sqrt(x_min + difference * sin(theta) ** 2)

    total = integrand(0.0) + integrand(pi / 2.0)
    for index in range(1, panels):
        total += (4.0 if index % 2 else 2.0) * integrand(index * step)
    result = total * step / 3.0
    if not isfinite(result):
        raise ValueError("half-cycle proper time is outside the finite range")
    return result


def half_cycle_convergence(spec: EinsteinCartanSpec, coarse_panels: Integral = 128) -> dict[str, float]:
    """Return two Simpson refinements and their finite deterministic difference."""

    if isinstance(coarse_panels, bool) or not isinstance(coarse_panels, Integral):
        raise ValueError("coarse_panels must be a positive even integer")
    coarse_panels = int(coarse_panels)
    coarse = half_cycle_proper_time(spec, coarse_panels)
    fine = half_cycle_proper_time(spec, 2 * coarse_panels)
    finer = half_cycle_proper_time(spec, 4 * coarse_panels)
    return _finite_mapping({
        "coarse": coarse,
        "fine": fine,
        "finer": finer,
        "coarse_to_fine_difference": abs(fine - coarse),
        "fine_to_finer_difference": abs(finer - fine),
    })  # type: ignore[return-value]
