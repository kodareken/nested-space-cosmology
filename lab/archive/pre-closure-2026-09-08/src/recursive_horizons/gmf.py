"""GMF-1 spherical matching obstruction and local shell diagnostics.

This module is deliberately not a global black-hole-to-child solution.  It
checks exact spherical GR identities: Oppenheimer--Snyder matching, the
Misner--Sharp mass obstruction between a finite Kottler mass and the CCT-1
pure de-Sitter core, parent null expansions, and a prescribed Israel shell.
The shell is distributional and does not supply smooth microphysics or Q.
Units are c = hbar = 1; G is passed explicitly where dimensional.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite, pi, sqrt
from numbers import Integral, Real
from typing import Literal

from .transition import classify_surface

Region = Literal["normal", "trapped", "marginal"]


def _real(name: str, value: Real) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a finite real number")
    result = float(value)
    if not isfinite(result):
        raise ValueError(f"{name} must be a finite real number")
    return result


def _positive(name: str, value: float) -> float:
    value = _real(name, value)
    if value <= 0.0:
        raise ValueError(f"{name} must be positive")
    return value


def _nonnegative(name: str, value: float) -> float:
    value = _real(name, value)
    if value < 0.0:
        raise ValueError(f"{name} must be non-negative")
    return value


def _orientation(name: str, value: int) -> int:
    if (
        isinstance(value, bool)
        or not isinstance(value, Integral)
        or value not in (-1, 1)
    ):
        raise ValueError(f"{name} must be either -1 or +1")
    return int(value)


@dataclass(frozen=True, slots=True)
class KottlerRegion:
    """Spherical vacuum region with finite mass M >= 0 and Lambda >= 0."""

    mass: float
    cosmological_constant: float = 0.0
    gravitational_constant: float = 1.0

    def __post_init__(self) -> None:
        object.__setattr__(self, "mass", _nonnegative("mass", self.mass))
        object.__setattr__(self, "cosmological_constant", _nonnegative("cosmological_constant", self.cosmological_constant))
        object.__setattr__(self, "gravitational_constant", _positive("gravitational_constant", self.gravitational_constant))


@dataclass(frozen=True, slots=True)
class ClosedDeSitterChild:
    """The CCT-1 pure de-Sitter child benchmark with radius L > 0."""

    radius: float
    gravitational_constant: float = 1.0

    def __post_init__(self) -> None:
        object.__setattr__(self, "radius", _positive("radius", self.radius))
        object.__setattr__(self, "gravitational_constant", _positive("gravitational_constant", self.gravitational_constant))

    @property
    def cosmological_constant(self) -> float:
        return 3.0 / self.radius**2


@dataclass(frozen=True, slots=True)
class TimelikeShell:
    """Prescribed timelike Israel shell; this is distributional, not smooth."""

    surface_density: float
    epsilon_child: int = 1
    epsilon_parent: int = 1

    def __post_init__(self) -> None:
        object.__setattr__(self, "surface_density", _positive("surface_density", self.surface_density))
        object.__setattr__(self, "epsilon_child", _orientation("epsilon_child", self.epsilon_child))
        object.__setattr__(self, "epsilon_parent", _orientation("epsilon_parent", self.epsilon_parent))


def kottler_f(region: KottlerRegion, areal_radius: float) -> float:
    """Return f=1-2GM/R-Lambda R^2/3 in a Kottler region."""

    radius = _positive("areal_radius", areal_radius)
    return 1.0 - 2.0 * region.gravitational_constant * region.mass / radius - region.cosmological_constant * radius**2 / 3.0


def kottler_misner_sharp_mass(region: KottlerRegion, areal_radius: float) -> float:
    """Return geometric Misner--Sharp mass M + Lambda R^3/(6G)."""

    radius = _positive("areal_radius", areal_radius)
    return region.mass + region.cosmological_constant * radius**3 / (6.0 * region.gravitational_constant)


def child_misner_sharp_mass(child: ClosedDeSitterChild, areal_radius: float) -> float:
    """Return de-Sitter Misner--Sharp mass Lambda_child R^3/(6G)."""

    radius = _positive("areal_radius", areal_radius)
    return child.cosmological_constant * radius**3 / (6.0 * child.gravitational_constant)


def darmois_mass_jump(parent: KottlerRegion, child: ClosedDeSitterChild, areal_radius: float) -> float:
    """Return parent minus child Misner--Sharp mass at a proposed smooth join."""

    if parent.gravitational_constant != child.gravitational_constant:
        raise ValueError("parent and child must use the same gravitational_constant")
    return kottler_misner_sharp_mass(parent, areal_radius) - child_misner_sharp_mass(child, areal_radius)


def darmois_dynamic_gate(parent: KottlerRegion, child: ClosedDeSitterChild, radius_a: float, radius_b: float, *, tolerance: float = 1e-12) -> dict[str, object]:
    """Evaluate the necessary Misner--Sharp mass gate at two radii.

    With constant parent M and Lambdas, a nonzero mass jump cannot vanish at
    two distinct radii unless M=0 and the vacuum constants agree.  Passing
    this necessary gate is not, for generic matter regions, a proof that all
    induced-metric and extrinsic-curvature junction conditions hold.
    """

    radius_a = _positive("radius_a", radius_a)
    radius_b = _positive("radius_b", radius_b)
    tolerance = _nonnegative("tolerance", tolerance)
    if radius_a == radius_b:
        raise ValueError("radius_a and radius_b must be distinct")
    jump_a = darmois_mass_jump(parent, child, radius_a)
    jump_b = darmois_mass_jump(parent, child, radius_b)
    mass_gate_both = abs(jump_a) <= tolerance and abs(jump_b) <= tolerance
    trivial = (
        parent.mass == 0.0
        and parent.cosmological_constant == child.cosmological_constant
    )
    return {
        "mass_jump_radius_a": jump_a,
        "mass_jump_radius_b": jump_b,
        "mass_gate_passes_at_both_radii": mass_gate_both,
        "input_is_trivial_same_geometry_control": trivial,
        "full_darmois_match_proven": False,
        "classification": "darmois_mass_gate_not_global_solution",
    }


def os_boundary_radius(scale_factor: float, chi_boundary: float) -> float:
    """Return ``R_b=a sin(chi_b)`` for an outward-monotone closed-FLRW cap.

    The controlled Oppenheimer--Snyder branch keeps ``0 < chi_b < pi/2`` so
    the areal radius increases toward the Schwarzschild exterior.  Other
    orientation branches need a separate junction analysis.
    """

    from math import sin

    scale_factor = _positive("scale_factor", scale_factor)
    chi_boundary = _real("chi_boundary", chi_boundary)
    if not 0.0 < chi_boundary < pi / 2.0:
        raise ValueError("chi_boundary must satisfy 0 < chi_boundary < pi/2")
    return scale_factor * sin(chi_boundary)


def os_mass(density: float, scale_factor: float, chi_boundary: float, gravitational_constant: float = 1.0) -> float:
    """Return M=4 pi rho R_b^3/3 for the Oppenheimer--Snyder dust ball."""

    density = _positive("density", density)
    gravitational_constant = _positive("gravitational_constant", gravitational_constant)
    return 4.0 * pi * density * os_boundary_radius(scale_factor, chi_boundary) ** 3 / 3.0


def os_friedmann_residual(density: float, scale_factor: float, scale_factor_velocity: float, gravitational_constant: float = 1.0) -> float:
    """Return adot^2+1-8 pi G rho a^2/3 for closed dust FLRW."""

    density = _positive("density", density)
    scale_factor = _positive("scale_factor", scale_factor)
    velocity = _real("scale_factor_velocity", scale_factor_velocity)
    gravitational_constant = _positive("gravitational_constant", gravitational_constant)
    return velocity**2 + 1.0 - 8.0 * pi * gravitational_constant * density * scale_factor**2 / 3.0


def os_boundary_geodesic_residual(density: float, scale_factor: float, scale_factor_velocity: float, chi_boundary: float, gravitational_constant: float = 1.0) -> float:
    """Return Rdot^2-[2GM/R-sin(chi_b)^2] for the OS Schwarzschild boundary."""

    from math import sin

    gravitational_constant = _positive("gravitational_constant", gravitational_constant)
    radius = os_boundary_radius(scale_factor, chi_boundary)
    mass = os_mass(density, scale_factor, chi_boundary, gravitational_constant)
    velocity = _real("scale_factor_velocity", scale_factor_velocity)
    return (velocity * sin(chi_boundary)) ** 2 - (2.0 * gravitational_constant * mass / radius - sin(chi_boundary) ** 2)


def os_dust_density(
    reference_density: float,
    reference_scale_factor: float,
    scale_factor: float,
) -> float:
    """Return the exact conserved-dust scaling ``rho=rho_ref(a_ref/a)^3``."""

    density = _positive("reference_density", reference_density)
    reference_scale = _positive("reference_scale_factor", reference_scale_factor)
    scale = _positive("scale_factor", scale_factor)
    try:
        result = density * (reference_scale / scale) ** 3
    except OverflowError as error:
        raise ValueError("dust density is outside the finite evaluation range") from error
    if not isfinite(result):
        raise ValueError("dust density is outside the finite evaluation range")
    return result


def os_dust_ricci_scalar(
    reference_density: float,
    reference_scale_factor: float,
    scale_factor: float,
    gravitational_constant: float = 1.0,
) -> float:
    """Return ``R=8 pi G rho`` for the zero-Lambda pressureless OS interior."""

    gravitational_constant = _positive("gravitational_constant", gravitational_constant)
    result = 8.0 * pi * gravitational_constant * os_dust_density(
        reference_density,
        reference_scale_factor,
        scale_factor,
    )
    if not isfinite(result):
        raise ValueError("dust Ricci scalar is outside the finite evaluation range")
    return result


def os_endpoint_status() -> dict[str, object]:
    """State the exact physical limit of the OS benchmark."""

    return {
        "finite_parent_mass_matching": True,
        "smooth_match": True,
        "matching_status": "established_exact_OS_solution_reduced_identity_regression",
        "induced_metric_and_extrinsic_curvature_residuals_implemented": False,
        "regular_transition": False,
        "singularity_present_at_a_zero": True,
        "global_child_solution": False,
    }


def parent_ef_null_expansions(parent: KottlerRegion, areal_radius: float) -> tuple[float, float]:
    """Return future expansions in ingoing Eddington--Finkelstein normalization.

    The chosen null normals satisfy k_plus dot k_minus=-2, yielding
    theta_plus=sqrt(2) f/R and theta_minus=-2 sqrt(2)/R.
    """

    radius = _positive("areal_radius", areal_radius)
    return (sqrt(2.0) * kottler_f(parent, radius) / radius, -2.0 * sqrt(2.0) / radius)


def parent_expansion_product_residual(parent: KottlerRegion, areal_radius: float) -> float:
    """Check theta_plus theta_minus = -4 f/R^2."""

    radius = _positive("areal_radius", areal_radius)
    plus, minus = parent_ef_null_expansions(parent, radius)
    return plus * minus + 4.0 * kottler_f(parent, radius) / radius**2


def classify_parent_round_sphere(parent: KottlerRegion, areal_radius: float, *, tolerance: float = 1e-12) -> Region:
    """Classify normal/marginal/trapped parent spheres using EF expansions."""

    tolerance = _nonnegative("tolerance", tolerance)
    plus, minus = parent_ef_null_expansions(parent, areal_radius)
    if abs(plus) <= tolerance or abs(minus) <= tolerance:
        return "marginal"
    if plus < 0.0 and minus < 0.0:
        return "trapped"
    return "normal"


def child_region_sample(time: float, chi: float, child: ClosedDeSitterChild) -> str:
    """Classify a CCT-1 child sphere without asserting a join to the parent."""

    return classify_surface(time, chi, child.radius)


def shell_kappa(shell: TimelikeShell, gravitational_constant: float = 1.0) -> float:
    return 4.0 * pi * _positive("gravitational_constant", gravitational_constant) * shell.surface_density


def _child_f(child: ClosedDeSitterChild, radius: float) -> float:
    return 1.0 - child.cosmological_constant * radius**2 / 3.0


def shell_junction_residual(parent: KottlerRegion, child: ClosedDeSitterChild, shell: TimelikeShell, areal_radius: float, radial_velocity: float) -> float:
    """Return the unsquared timelike Israel residual epsilon_- beta_- - epsilon_+ beta_+ - kappa R."""

    if parent.gravitational_constant != child.gravitational_constant:
        raise ValueError("parent and child must use the same gravitational_constant")
    radius = _positive("areal_radius", areal_radius)
    velocity = _real("radial_velocity", radial_velocity)
    child_argument = velocity**2 + _child_f(child, radius)
    parent_argument = velocity**2 + kottler_f(parent, radius)
    if child_argument < 0.0 or parent_argument < 0.0:
        raise ValueError("timelike shell has negative beta squared")
    return shell.epsilon_child * sqrt(child_argument) - shell.epsilon_parent * sqrt(parent_argument) - shell_kappa(shell, parent.gravitational_constant) * radius


def shell_effective_potential(parent: KottlerRegion, child: ClosedDeSitterChild, shell: TimelikeShell, areal_radius: float) -> float:
    """Return the squared Israel potential; callers must also check the unsquared equation."""

    if parent.gravitational_constant != child.gravitational_constant:
        raise ValueError("parent and child must use the same gravitational_constant")
    radius = _positive("areal_radius", areal_radius)
    kappa = shell_kappa(shell, parent.gravitational_constant)
    numerator = _child_f(child, radius) - kottler_f(parent, radius) - kappa**2 * radius**2
    return kottler_f(parent, radius) - numerator**2 / (4.0 * kappa**2 * radius**2)


def shell_orientation_gate(parent: KottlerRegion, child: ClosedDeSitterChild, shell: TimelikeShell, areal_radius: float, radial_velocity: float, *, tolerance: float = 1e-10) -> dict[str, object]:
    """Accept a local shell state only if its unsquared and squared equations agree."""

    tolerance = _nonnegative("tolerance", tolerance)
    residual = shell_junction_residual(parent, child, shell, areal_radius, radial_velocity)
    potential = shell_effective_potential(parent, child, shell, areal_radius)
    return {"unsquared_residual": residual, "effective_potential": potential, "orientation_checked": True, "accepted": abs(residual) <= tolerance and abs(radial_velocity**2 + potential) <= tolerance, "thin_shell_distributional": True}


def shell_conservation_residual(surface_density: float, surface_pressure: float, areal_radius: float, radial_velocity: float, density_derivative: float = 0.0, flux_jump: float = 0.0) -> float:
    """Return d(sigma A)/dlambda+p dA/dlambda+A[T_un] for a spherical shell."""

    sigma = _positive("surface_density", surface_density)
    pressure = _real("surface_pressure", surface_pressure)
    radius = _positive("areal_radius", areal_radius)
    velocity = _real("radial_velocity", radial_velocity)
    sigma_dot = _real("density_derivative", density_derivative)
    flux = _real("flux_jump", flux_jump)
    area_dot = 8.0 * pi * radius * velocity
    return 4.0 * pi * radius**2 * sigma_dot + (sigma + pressure) * area_dot + 4.0 * pi * radius**2 * flux


def shell_static_gate(parent: KottlerRegion, child: ClosedDeSitterChild, shell: TimelikeShell, areal_radius: float, *, step: float = 1e-5, tolerance: float = 1e-8) -> dict[str, object]:
    """Return an oriented fixed-kappa potential diagnostic for a shell.

    The effective potential follows from a squared equation.  A static state
    is therefore accepted only when the original oriented Israel equation at
    zero radial velocity also passes.  The finite differences hold the shell
    surface density (and therefore kappa) fixed.  They are not a general shell
    stability calculation with an independently derived equation of state,
    flux law, or neighborhood admissibility proof.
    """

    radius = _positive("areal_radius", areal_radius)
    step = _positive("step", step)
    tolerance = _nonnegative("tolerance", tolerance)
    if step >= radius:
        raise ValueError("step must be smaller than areal_radius")
    center = shell_effective_potential(parent, child, shell, radius)
    left = shell_effective_potential(parent, child, shell, radius - step)
    right = shell_effective_potential(parent, child, shell, radius + step)
    first = (right - left) / (2.0 * step)
    second = (right - 2.0 * center + left) / step**2
    oriented_residual = shell_junction_residual(
        parent,
        child,
        shell,
        radius,
        0.0,
    )
    orientation_accepted = abs(oriented_residual) <= tolerance
    actually_static = (
        orientation_accepted
        and abs(center) <= tolerance
        and abs(first) <= tolerance
    )
    return {
        "potential": center,
        "potential_derivative": first,
        "potential_second_derivative": second,
        "oriented_junction_residual_at_rest": oriented_residual,
        "orientation_accepted_at_rest": orientation_accepted,
        "actually_static": actually_static,
        "fixed_kappa_potential_minimum_if_static": actually_static and second > 0.0,
        "fixed_kappa_during_derivative": True,
        "derivative_step": step,
        "neighborhood_orientation_admissibility_checked": False,
        "general_shell_stability_claimed": False,
        "thin_shell_distributional": True,
    }
