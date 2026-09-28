"""Matched SGB-L initial profiles and radial constraint-family continuation.

The linear branch is ``f(phi)=alpha_gb*phi`` with ``beta=eta=0``.  Profiles
are the compact SF1 bumps evaluated by the action-free helper
``compact_bump_with_derivatives``.  On the shared unit-lapse, zero-shift
slice, ``Pi_chi=chi_t=A_chi B_r/r`` is fixed before any geometry is solved
and does not depend on ``lambda``.  The two radial ODEs are the already
proved affine, diagonal ``H=H0+H_L L_r`` and ``M=M0+M_k k_r`` constraints.

RK4 and SSPRK3 here are radial *constraint integrators* on a small fixed
mesh.  They are not a time-evolution method, a trajectory, or a production
grid.  The interior is the exact Minkowski buffer; the exterior is the
initial GR vacuum continuation, which is not a static full SGB-L solution.
This slice does not set ``SGBL_branch_owned_and_healthy``.
"""

from __future__ import annotations

from dataclasses import dataclass, fields
from fractions import Fraction
from math import isfinite, sqrt
from numbers import Real
from typing import Any, Callable, Sequence

from .initial_data_family import gauge_compatible_metric_time_derivatives
from .initial_data_preflight import compact_bump_with_derivatives
from .modified_harmonic_constraints import physical_constraint_projections
from .sgb1_ctl1_constraints import SGBLAnnularInputs, sgbl_constraint_coefficients
from .spherical_reduction import Jet2, SphericalState, residuals


Q = Fraction

DECLARED_CHI_AMPLITUDE = 3.0
DECLARED_PHI_AMPLITUDE = 1.0 / 131072.0
DECLARED_CENTER = 12.0
DECLARED_HALF_WIDTH = 2.0
DECLARED_OUTER_RADIUS = 128.0
DECLARED_PLANCK_MASS = 2.0
DECLARED_SCALAR_MASS = 3.0
DECLARED_QUARTIC_COUPLING = 0.5
DECLARED_ALPHA_GB = -0.25
MAX_RADIAL_CONSTRAINT_STEPS = 64
CONSTRAINT_INTEGRATOR_NAMES = ("RK4", "SSPRK3")
FAMILY_STOP_REASONS = frozenset(
    {
        "singular_constraint_jacobian",
        "uncertified_floating_jacobian",
        "nonpositive_metric",
        "resource_limit",
        "nonfinite_continuation",
        "lost_asymptotic_condition",
    }
)


def _finite(name: str, value: object, *, positive: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a finite real scalar")
    result = float(value)
    if not isfinite(result):
        raise ValueError(f"{name} must be finite")
    if positive and result <= 0.0:
        raise ValueError(f"{name} must be positive")
    return result


def _finite_pair(name: str, value: Sequence[Real]) -> tuple[float, float]:
    if isinstance(value, (str, bytes)) or len(value) != 2:
        raise ValueError(f"{name} must contain two entries")
    return _finite(f"{name}[0]", value[0]), _finite(f"{name}[1]", value[1])


def _fraction(name: str, value: object) -> Fraction:
    if type(value) not in (int, Fraction):
        raise TypeError(f"{name} must be a Fraction or built-in int")
    return Fraction(value)


@dataclass(frozen=True, slots=True)
class SGBLFamilyParameters:
    """Declared SGB-L compact-family inputs.

    ``chi_amplitude`` is the one family member input (nominal 3).  The
    remaining numbers are frozen fixtures, not a scan or a fit.
    """

    chi_amplitude: float = DECLARED_CHI_AMPLITUDE
    phi_amplitude: float = DECLARED_PHI_AMPLITUDE
    center: float = DECLARED_CENTER
    chi_half_width: float = DECLARED_HALF_WIDTH
    phi_half_width: float = DECLARED_HALF_WIDTH
    outer_radius: float = DECLARED_OUTER_RADIUS
    planck_mass: float = DECLARED_PLANCK_MASS
    scalar_mass: float = DECLARED_SCALAR_MASS
    quartic_coupling: float = DECLARED_QUARTIC_COUPLING
    alpha_gb: float = DECLARED_ALPHA_GB

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "chi_amplitude", _finite("chi_amplitude", self.chi_amplitude, positive=True)
        )
        for name in (
            "phi_amplitude",
            "center",
            "chi_half_width",
            "phi_half_width",
            "outer_radius",
            "planck_mass",
            "scalar_mass",
            "quartic_coupling",
            "alpha_gb",
        ):
            object.__setattr__(self, name, _finite(name, getattr(self, name)))
        if self.phi_amplitude < 0.0:
            raise ValueError("phi_amplitude must be nonnegative")
        frozen = (
            ("phi_amplitude", DECLARED_PHI_AMPLITUDE),
            ("center", DECLARED_CENTER),
            ("chi_half_width", DECLARED_HALF_WIDTH),
            ("phi_half_width", DECLARED_HALF_WIDTH),
            ("outer_radius", DECLARED_OUTER_RADIUS),
            ("planck_mass", DECLARED_PLANCK_MASS),
            ("scalar_mass", DECLARED_SCALAR_MASS),
            ("quartic_coupling", DECLARED_QUARTIC_COUPLING),
            ("alpha_gb", DECLARED_ALPHA_GB),
        )
        for name, expected in frozen:
            if getattr(self, name) != expected:
                raise ValueError(
                    f"{name} is the declared SGB-L fixture {expected}, not a selected parameter"
                )
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
    def beta(self) -> float:
        return 0.0

    @property
    def eta(self) -> float:
        return 0.0


@dataclass(frozen=True, slots=True)
class SGBLAffineConstraintCoefficients:
    """Generic arithmetic record of ``H=H0+H_L L_r`` and ``M=M0+M_k k_r``."""

    hamiltonian_constant: Any
    hamiltonian_lambda_r: Any
    momentum_constant: Any
    momentum_k_r: Any

    @property
    def jacobian_diagonal(self) -> tuple[Any, Any]:
        return self.hamiltonian_lambda_r, self.momentum_k_r

    @property
    def jacobian_determinant(self) -> Any:
        return self.hamiltonian_lambda_r * self.momentum_k_r

    def evaluate(
        self,
        *,
        radial_metric_derivative: Any,
        angular_extrinsic_curvature_derivative: Any,
    ) -> tuple[Any, Any]:
        return (
            self.hamiltonian_constant + self.hamiltonian_lambda_r * radial_metric_derivative,
            self.momentum_constant
            + self.momentum_k_r * angular_extrinsic_curvature_derivative,
        )

    def derivatives(self) -> tuple[Any, Any]:
        """Return ``(L_r, k_r)`` or refuse a zero diagonal.

        An exact ``Fraction``/``int`` zero is the annular kernel singularity.
        A floating ``0.0`` is not that certificate: it raises
        ``uncertified_floating_jacobian``.  There is no conditioning floor.
        """

        diagonal = self.jacobian_diagonal
        exact = all(type(value) in (int, Fraction) for value in diagonal)
        if exact:
            if diagonal[0] == 0 or diagonal[1] == 0:
                raise SGBLConstraintFamilyStop(
                    "singular_constraint_jacobian",
                    "exact SGB-L physical-constraint Jacobian is singular",
                    {"diagonal": diagonal, "exact_kernel": True},
                )
        else:
            floats = tuple(float(value) for value in diagonal)
            if not all(isfinite(value) for value in floats):
                raise SGBLConstraintFamilyStop(
                    "nonfinite_continuation",
                    "SGB-L constraint Jacobian became nonfinite",
                    {"diagonal": floats},
                )
            if floats[0] == 0.0 or floats[1] == 0.0:
                raise SGBLConstraintFamilyStop(
                    "uncertified_floating_jacobian",
                    "a floating zero diagonal is not a certificate that the exact SGB-L Jacobian is singular",
                    {"diagonal": floats, "exact_kernel": False},
                )
        return (
            -self.hamiltonian_constant / self.hamiltonian_lambda_r,
            -self.momentum_constant / self.momentum_k_r,
        )


@dataclass(frozen=True, slots=True)
class SGBLFamilyPoint:
    """One stored radial node.  Compactness and mass must match the ADM formula."""

    radius: float
    radial_metric: float
    angular_extrinsic_curvature: float
    radial_metric_derivative: float
    angular_extrinsic_curvature_derivative: float
    hamiltonian_residual: float
    momentum_residual: float
    compactness: float
    misner_sharp_mass: float
    chi: float
    chi_pi: float
    phi: float
    phi_pi: float

    def __post_init__(self) -> None:
        for field in fields(self):
            object.__setattr__(self, field.name, _finite(field.name, getattr(self, field.name)))
        if self.radius <= 0.0:
            raise ValueError("radius must be positive")
        if self.radial_metric <= 0.0:
            raise ValueError("radial_metric must be positive")
        expected_compactness = (
            1.0
            + self.radius**2 * self.angular_extrinsic_curvature**2
            - 1.0 / self.radial_metric**2
        )
        if self.compactness != expected_compactness:
            raise ValueError("compactness does not match the polar-areal Misner-Sharp formula")
        expected_mass = self.radius * self.compactness / 2.0
        if self.misner_sharp_mass != expected_mass:
            raise ValueError("misner_sharp_mass does not match radius*compactness/2")


@dataclass(frozen=True, slots=True)
class SGBLFamilySolution:
    """Stored constraint-integrator nodes and declared parameters.

    Scope claims are derived properties and cannot be attached by
    construction or ``replace``.  Sampled-node untrappedness is not a
    continuous no-initial-trap certificate.  Exterior numbers are a
    floating evaluation of the analytic GR vacuum formula.
    """

    method: str
    step_count: int
    parameters: SGBLFamilyParameters
    points: tuple[SGBLFamilyPoint, ...]

    def __post_init__(self) -> None:
        if self.method not in CONSTRAINT_INTEGRATOR_NAMES:
            raise ValueError("unknown SGB-L radial constraint integrator")
        if isinstance(self.step_count, bool) or not isinstance(self.step_count, int) or self.step_count <= 0:
            raise ValueError("step_count must be a positive integer")
        if self.step_count > MAX_RADIAL_CONSTRAINT_STEPS:
            raise SGBLConstraintFamilyStop(
                "resource_limit",
                "SGB-L family mesh exceeds the small fixed constraint-integrator budget",
                {"step_count": self.step_count, "maximum": MAX_RADIAL_CONSTRAINT_STEPS},
            )
        if not isinstance(self.parameters, SGBLFamilyParameters):
            raise TypeError("parameters must be SGBLFamilyParameters")
        if not isinstance(self.points, tuple) or any(
            type(point) is not SGBLFamilyPoint for point in self.points
        ):
            raise TypeError("points must be a tuple of SGBLFamilyPoint")
        if len(self.points) != self.step_count + 1:
            raise ValueError("point count must equal step_count + 1")
        start = self.parameters.support_minimum
        end = self.parameters.support_maximum
        if self.points[0].radius != start or self.points[-1].radius != end:
            raise ValueError("mesh must run from the compact-support minimum to its maximum")
        widths = tuple(
            self.points[index + 1].radius - self.points[index].radius
            for index in range(len(self.points) - 1)
        )
        if any(width <= 0.0 for width in widths) or any(width != widths[0] for width in widths):
            raise ValueError("constraint-integrator mesh must be strictly increasing and uniform")
        first = self.points[0]
        if first.radial_metric != 1.0 or first.angular_extrinsic_curvature != 0.0:
            raise ValueError("family continuation must start at the Minkowski buffer")
        for point in self.points:
            expected = _point_from_state(
                point.radius,
                (point.radial_metric, point.angular_extrinsic_curvature),
                self.parameters,
            )
            if point != expected:
                raise ValueError("stored node does not match its declared family and defining RHS")
        vacuum_minimum = _vacuum_denominator_minimum(
            self.points[-1].radius,
            self.parameters.outer_radius,
            self.points[-1].misner_sharp_mass,
            self.points[-1].angular_extrinsic_curvature * self.points[-1].radius**3,
        )
        if not isfinite(vacuum_minimum) or vacuum_minimum <= 0.0:
            raise SGBLConstraintFamilyStop(
                "lost_asymptotic_condition",
                "analytic exterior vacuum continuation lost its positive metric factor",
                {"minimum_denominator": vacuum_minimum},
            )

    @property
    def peak_point(self) -> SGBLFamilyPoint:
        return max(self.points, key=lambda item: item.compactness)

    @property
    def peak_compactness(self) -> float:
        return self.peak_point.compactness

    @property
    def peak_radius(self) -> float:
        return self.peak_point.radius

    @property
    def outer_mass(self) -> float:
        return self.points[-1].misner_sharp_mass

    @property
    def vacuum_momentum_constant(self) -> float:
        last = self.points[-1]
        return last.angular_extrinsic_curvature * last.radius**3

    @property
    def floating_outer_radial_metric(self) -> float:
        return self._floating_exterior()[0]

    @property
    def floating_outer_angular_extrinsic_curvature(self) -> float:
        return self._floating_exterior()[1]

    @property
    def floating_minimum_vacuum_metric_denominator(self) -> float:
        return _vacuum_denominator_minimum(
            self.points[-1].radius,
            self.parameters.outer_radius,
            self.outer_mass,
            self.vacuum_momentum_constant,
        )

    @property
    def outer_radial_metric(self) -> float:
        return self.floating_outer_radial_metric

    @property
    def outer_angular_extrinsic_curvature(self) -> float:
        return self.floating_outer_angular_extrinsic_curvature

    @property
    def minimum_vacuum_metric_denominator(self) -> float:
        return self.floating_minimum_vacuum_metric_denominator

    @property
    def maximum_defining_rhs_residual_infinity(self) -> float:
        return max(
            max(abs(item.hamiltonian_residual), abs(item.momentum_residual))
            for item in self.points
        )

    @property
    def finite_mass(self) -> bool:
        return isfinite(self.outer_mass)

    @property
    def exact_inner_minkowski_buffer(self) -> bool:
        first = self.points[0]
        return (
            first.radius == self.parameters.support_minimum
            and first.radial_metric == 1.0
            and first.angular_extrinsic_curvature == 0.0
        )

    @property
    def analytic_exterior_vacuum_formula(self) -> str:
        return "k=J/r^3, lambda^{-2}=1-2M/r+J^2/r^4"

    @property
    def exact_outer_gr_vacuum_continuation(self) -> bool:
        """The analytic formula is used; the stored numbers are not exact rationals."""

        return False

    @property
    def sampled_nodes_untrapped(self) -> bool:
        return all(item.compactness < 1.0 for item in self.points)

    @property
    def continuous_no_initial_trapped_sphere(self) -> bool:
        return False

    @property
    def is_full_static_sgbl_solution(self) -> bool:
        return False

    @property
    def not_a_time_evolution_method(self) -> bool:
        return True

    @property
    def branch_owned_and_healthy(self) -> bool:
        return False

    def _floating_exterior(self) -> tuple[float, float, float]:
        return sgbl_floating_exterior_vacuum_evaluation(
            self.parameters.outer_radius,
            mass=self.outer_mass,
            momentum_constant=self.vacuum_momentum_constant,
        )


class SGBLConstraintFamilyStop(ArithmeticError):
    """Typed fail-closed stop of the SGB-L radial constraint continuation."""

    def __init__(
        self,
        reason: str,
        message: str,
        diagnostics: dict[str, Any] | None = None,
    ) -> None:
        if reason not in FAMILY_STOP_REASONS:
            raise ValueError("unknown SGB-L constraint-family stop reason")
        super().__init__(message)
        self.reason = reason
        self.diagnostics = dict(diagnostics or {})


def sgbl_compact_family_fields(
    radius: Real,
    parameters: SGBLFamilyParameters,
) -> dict[str, float]:
    """Return the compact scalar profiles.  Geometry, including ``lambda``, is not an input."""

    if not isinstance(parameters, SGBLFamilyParameters):
        raise TypeError("parameters must be SGBLFamilyParameters")
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
            chi_bump_rr / r - 2.0 * chi_bump_r / r**2 + 2.0 * chi_bump / r**3
        )
        # Unit-lapse, zero-shift, future-ingoing: d_t(r*chi)=d_r(r*chi).
        chi_pi = amplitude * chi_bump_r / r
        chi_pi_r = amplitude * (chi_bump_rr / r - chi_bump_r / r**2)
    return {
        "phi": parameters.phi_amplitude * phi_bump,
        "phi_r": parameters.phi_amplitude * phi_bump_r,
        "phi_rr": parameters.phi_amplitude * phi_bump_rr,
        "phi_pi": 0.0,
        "phi_pi_r": 0.0,
        "chi": chi,
        "chi_r": chi_r,
        "chi_rr": chi_rr,
        "chi_pi": chi_pi,
        "chi_pi_r": chi_pi_r,
    }


def sgbl_future_normal_chi_momentum(
    radius: Real,
    parameters: SGBLFamilyParameters,
) -> float:
    """``Pi_chi=chi_t=A_chi B_r/r`` on the unit-lapse, zero-shift slice."""

    r = _finite("radius", radius)
    if r < 0.0:
        raise ValueError("radius must be nonnegative")
    _bump, bump_r, _bump_rr = compact_bump_with_derivatives(
        r,
        center=parameters.center,
        half_width=parameters.chi_half_width,
    )
    if r == 0.0:
        return 0.0
    return parameters.chi_amplitude * bump_r / r


def _affine_from_values(
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
    scalar_mass: Any,
    quartic_coupling: Any,
    alpha_gb: Any,
) -> SGBLAffineConstraintCoefficients:
    r = radius
    L = radial_metric
    k = angular_extrinsic_curvature
    pr = phi_r
    pi = phi_pi
    mpl2 = planck_mass**2
    alpha = alpha_gb
    a0 = -2 * k**2
    a_l = 1 / (L**3 * r)
    b = k**2 + (1 - 1 / L**2) / r**2
    c0 = 3 * k / (L * r)
    c_k = 1 / L
    x0 = phi_rr / L**2 - 2 * k * pi
    x_l = -pr / L**3
    y = k * pi + pr / (L**2 * r)
    z = (phi_pi_r - 2 * k * pr) / L
    rho = (
        (pi**2 + chi_pi**2) / 2
        + (pr**2 + chi_r**2) / (2 * L**2)
        + scalar_mass**2 * phi**2 / 2
        + quartic_coupling * phi**4 / 4
    )
    return SGBLAffineConstraintCoefficients(
        hamiltonian_constant=mpl2 * (2 * a0 + b) - rho - 8 * alpha * (b * x0 + 2 * a0 * y),
        hamiltonian_lambda_r=2 * mpl2 * a_l - 8 * alpha * (b * x_l + 2 * a_l * y),
        momentum_constant=(
            2 * mpl2 * L * c0
            - pi * pr
            - chi_pi * chi_r
            - 8 * alpha * L * (b * z + 2 * c0 * y)
        ),
        momentum_k_r=2 * mpl2 * L * c_k - 16 * alpha * L * c_k * y,
    )


def sgbl_affine_constraint_coefficients(
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
    scalar_mass: Any,
    quartic_coupling: Any,
    alpha_gb: Any,
) -> SGBLAffineConstraintCoefficients:
    """Return the branch H/M coefficients in generic arithmetic.

    Exact ``Fraction`` or built-in ``int`` inputs are routed through the
    annular kernel.  Floating inputs use the same algebra and are not a
    second physical theory.
    """

    values = (
        radius,
        radial_metric,
        angular_extrinsic_curvature,
        phi,
        phi_r,
        phi_rr,
        phi_pi,
        phi_pi_r,
        chi_r,
        chi_pi,
        planck_mass,
        scalar_mass,
        quartic_coupling,
        alpha_gb,
    )
    names = (
        "radius",
        "radial_metric",
        "angular_extrinsic_curvature",
        "phi",
        "phi_r",
        "phi_rr",
        "phi_pi",
        "phi_pi_r",
        "chi_r",
        "chi_pi",
        "planck_mass",
        "scalar_mass",
        "quartic_coupling",
        "alpha_gb",
    )
    if all(type(value) in (int, Fraction) for value in values):
        exact = {name: _fraction(name, value) for name, value in zip(names, values, strict=True)}
        coefficients = sgbl_constraint_coefficients(
            SGBLAnnularInputs(
                radius=exact["radius"],
                radial_metric=exact["radial_metric"],
                angular_extrinsic_curvature=exact["angular_extrinsic_curvature"],
                phi=exact["phi"],
                phi_r=exact["phi_r"],
                phi_rr=exact["phi_rr"],
                phi_pi=exact["phi_pi"],
                phi_pi_r=exact["phi_pi_r"],
                chi_r=exact["chi_r"],
                chi_pi=exact["chi_pi"],
                planck_mass=exact["planck_mass"],
                scalar_mass=exact["scalar_mass"],
                quartic_coupling=exact["quartic_coupling"],
                alpha_gb=exact["alpha_gb"],
            )
        )
        return SGBLAffineConstraintCoefficients(
            hamiltonian_constant=coefficients.hamiltonian_constant,
            hamiltonian_lambda_r=coefficients.hamiltonian_lambda_r,
            momentum_constant=coefficients.momentum_constant,
            momentum_k_r=coefficients.momentum_k_r,
        )
    return _affine_from_values(
        radius=radius,
        radial_metric=radial_metric,
        angular_extrinsic_curvature=angular_extrinsic_curvature,
        phi=phi,
        phi_r=phi_r,
        phi_rr=phi_rr,
        phi_pi=phi_pi,
        phi_pi_r=phi_pi_r,
        chi_r=chi_r,
        chi_pi=chi_pi,
        planck_mass=planck_mass,
        scalar_mass=scalar_mass,
        quartic_coupling=quartic_coupling,
        alpha_gb=alpha_gb,
    )


def sgbl_polar_areal_state(
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
    radial_metric_second_derivative: int | Fraction = 0,
    planck_mass: int | Fraction = 2,
    scalar_mass: int | Fraction = 3,
    quartic_coupling: int | Fraction = Q(1, 2),
    alpha_gb: int | Fraction = -Q(1, 4),
) -> SphericalState:
    """Exact unit-lapse, zero-shift polar--areal two-jet for tensor checks."""

    raw = {
        name: _fraction(name, value)
        for name, value in {
            "radius": radius,
            "radial_metric": radial_metric,
            "angular_extrinsic_curvature": angular_extrinsic_curvature,
            "radial_metric_derivative": radial_metric_derivative,
            "angular_extrinsic_curvature_derivative": angular_extrinsic_curvature_derivative,
            "phi": phi,
            "phi_r": phi_r,
            "phi_rr": phi_rr,
            "phi_pi": phi_pi,
            "phi_pi_r": phi_pi_r,
            "chi": chi,
            "chi_r": chi_r,
            "chi_rr": chi_rr,
            "chi_pi": chi_pi,
            "chi_pi_r": chi_pi_r,
            "radial_metric_second_derivative": radial_metric_second_derivative,
            "planck_mass": planck_mass,
            "scalar_mass": scalar_mass,
            "quartic_coupling": quartic_coupling,
            "alpha_gb": alpha_gb,
        }.items()
    }
    r = raw["radius"]
    L = raw["radial_metric"]
    k = raw["angular_extrinsic_curvature"]
    Lr = raw["radial_metric_derivative"]
    kr = raw["angular_extrinsic_curvature_derivative"]
    Lrr = raw["radial_metric_second_derivative"]
    if r <= 0 or L <= 0:
        raise ValueError("polar-areal tensor state requires r>0 and lambda>0")
    h_tt_dt, h_tr_dt = gauge_compatible_metric_time_derivatives(r, L, Lr)
    return SphericalState(
        h_tt=Jet2(-1, dt=h_tt_dt),
        h_tr=Jet2(0, dt=h_tr_dt),
        h_rr=Jet2(
            L**2,
            dt=4 * L**2 * k,
            dr=2 * L * Lr,
            dtr=8 * L * Lr * k + 4 * L**2 * kr,
            drr=2 * Lr**2 + 2 * L * Lrr,
        ),
        areal_radius=Jet2(r, dt=-r * k, dr=1, dtr=-k - r * kr),
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
        planck_mass=raw["planck_mass"],
        mu=raw["scalar_mass"],
        g4=raw["quartic_coupling"],
        alpha=raw["alpha_gb"],
        beta=0,
        eta=0,
        branch="SGB-L" if raw["alpha_gb"] else "GR-0",
    )


def sgbl_constraint_rhs(
    radius: Real,
    state: Sequence[Real],
    parameters: SGBLFamilyParameters,
) -> tuple[tuple[float, float], SGBLAffineConstraintCoefficients, dict[str, float]]:
    """Solve the two affine SGB-L constraints for ``lambda_r`` and ``k_r``.

    A floating zero diagonal is rechecked on the same finite binary-rational
    inputs before it can be called an exact kernel singularity. No
    quantitative Jacobian floor is introduced.
    """

    if not isinstance(parameters, SGBLFamilyParameters):
        raise TypeError("parameters must be SGBLFamilyParameters")
    r = _finite("radius", radius, positive=True)
    L, k = _finite_pair("state", state)
    if L <= 0.0:
        raise SGBLConstraintFamilyStop(
            "nonpositive_metric",
            "radial metric left the positive polar-areal branch",
            {"radius": r, "radial_metric": L},
        )
    fields = sgbl_compact_family_fields(r, parameters)
    try:
        coefficients = sgbl_affine_constraint_coefficients(
            radius=r,
            radial_metric=L,
            angular_extrinsic_curvature=k,
            phi=fields["phi"],
            phi_r=fields["phi_r"],
            phi_rr=fields["phi_rr"],
            phi_pi=fields["phi_pi"],
            phi_pi_r=fields["phi_pi_r"],
            chi_r=fields["chi_r"],
            chi_pi=fields["chi_pi"],
            planck_mass=parameters.planck_mass,
            scalar_mass=parameters.scalar_mass,
            quartic_coupling=parameters.quartic_coupling,
            alpha_gb=parameters.alpha_gb,
        )
    except (ArithmeticError, TypeError, ValueError, ZeroDivisionError) as exc:
        raise SGBLConstraintFamilyStop(
            "nonfinite_continuation",
            "SGB-L constraint kernel could not be evaluated",
            {"radius": r, "error": str(exc)},
        ) from exc
    diagonal = (
        float(coefficients.hamiltonian_lambda_r),
        float(coefficients.momentum_k_r),
    )
    if not all(isfinite(value) for value in diagonal):
        raise SGBLConstraintFamilyStop(
            "nonfinite_continuation",
            "SGB-L constraint Jacobian became nonfinite",
            {"radius": r, "diagonal": diagonal},
        )
    try:
        solved = coefficients.derivatives()
    except SGBLConstraintFamilyStop as exc:
        if exc.reason == "uncertified_floating_jacobian":
            exact = _exact_affine_coefficients_from_finite_inputs(
                radius=r,
                radial_metric=L,
                angular_extrinsic_curvature=k,
                fields=fields,
                parameters=parameters,
            )
            if exact is not None:
                try:
                    exact.derivatives()
                except SGBLConstraintFamilyStop as certified:
                    if certified.reason == "singular_constraint_jacobian":
                        certified.diagnostics = {
                            "radius": r,
                            "certified_on_finite_binary_rationals": True,
                            **certified.diagnostics,
                        }
                        raise certified from exc
            exc.diagnostics = {"radius": r, "exact_kernel": False, **exc.diagnostics}
        else:
            exc.diagnostics = {"radius": r, **exc.diagnostics}
        raise
    derivative = (float(solved[0]), float(solved[1]))
    if not all(isfinite(value) for value in derivative):
        raise SGBLConstraintFamilyStop(
            "nonfinite_continuation",
            "SGB-L constraint derivative became nonfinite",
            {"radius": r},
        )
    residual = coefficients.evaluate(
        radial_metric_derivative=derivative[0],
        angular_extrinsic_curvature_derivative=derivative[1],
    )
    return derivative, coefficients, {
        "residual_infinity": max(abs(float(value)) for value in residual),
        "jacobian_determinant": diagonal[0] * diagonal[1],
        "chi_pi": fields["chi_pi"],
        "chi": fields["chi"],
        "phi": fields["phi"],
        "phi_pi": fields["phi_pi"],
    }


def sgbl_misner_sharp_compactness(
    radius: Real,
    radial_metric: Real,
    angular_extrinsic_curvature: Real,
) -> float:
    r = _finite("radius", radius, positive=True)
    L = _finite("radial_metric", radial_metric, positive=True)
    k = _finite("angular_extrinsic_curvature", angular_extrinsic_curvature)
    answer = 1.0 + r**2 * k**2 - 1.0 / L**2
    if not isfinite(answer):
        raise SGBLConstraintFamilyStop(
            "nonfinite_continuation",
            "Misner-Sharp compactness became nonfinite",
        )
    return answer


def sgbl_analytic_exterior_vacuum_denominator(
    radius: Any,
    mass: Any,
    momentum_constant: Any,
) -> Any:
    """Analytic GR vacuum factor ``1-2M/r+J^2/r^4``.  Not a floating certificate."""

    return 1 - 2 * mass / radius + momentum_constant**2 / radius**4


def sgbl_floating_exterior_vacuum_evaluation(
    radius: Real,
    *,
    mass: Real,
    momentum_constant: Real,
) -> tuple[float, float, float]:
    """Floating evaluation of the analytic exterior formula at one radius.

    The formula is exact; the returned ``(lambda, k, denominator)`` are
    binary64 values.  This is not a static full SGB-L solution.
    """

    r = _finite("radius", radius, positive=True)
    mass_value = _finite("mass", mass)
    momentum = _finite("momentum_constant", momentum_constant)
    denominator = float(sgbl_analytic_exterior_vacuum_denominator(r, mass_value, momentum))
    if not isfinite(denominator):
        raise SGBLConstraintFamilyStop(
            "nonfinite_continuation",
            "outer vacuum denominator became nonfinite",
            {"radius": r},
        )
    if denominator <= 0.0:
        raise SGBLConstraintFamilyStop(
            "lost_asymptotic_condition",
            "analytic exterior vacuum continuation lost its positive metric factor",
            {"radius": r, "denominator": denominator},
        )
    return 1.0 / sqrt(denominator), momentum / r**3, denominator


def sgbl_initial_vacuum_outer_continuation(
    radius: Real,
    *,
    mass: Real,
    momentum_constant: Real,
) -> tuple[float, float, float]:
    """Floating evaluation of the analytic GR vacuum continuation."""

    return sgbl_floating_exterior_vacuum_evaluation(
        radius, mass=mass, momentum_constant=momentum_constant
    )


def _exact_affine_coefficients_from_finite_inputs(
    *,
    radius: float,
    radial_metric: float,
    angular_extrinsic_curvature: float,
    fields: dict[str, float],
    parameters: SGBLFamilyParameters,
) -> SGBLAffineConstraintCoefficients | None:
    """Exact annular kernel on the binary-rational images of finite inputs."""

    try:
        return sgbl_affine_constraint_coefficients(
            radius=Fraction(radius),
            radial_metric=Fraction(radial_metric),
            angular_extrinsic_curvature=Fraction(angular_extrinsic_curvature),
            phi=Fraction(fields["phi"]),
            phi_r=Fraction(fields["phi_r"]),
            phi_rr=Fraction(fields["phi_rr"]),
            phi_pi=Fraction(fields["phi_pi"]),
            phi_pi_r=Fraction(fields["phi_pi_r"]),
            chi_r=Fraction(fields["chi_r"]),
            chi_pi=Fraction(fields["chi_pi"]),
            planck_mass=Fraction(parameters.planck_mass),
            scalar_mass=Fraction(parameters.scalar_mass),
            quartic_coupling=Fraction(parameters.quartic_coupling),
            alpha_gb=Fraction(parameters.alpha_gb),
        )
    except (ArithmeticError, TypeError, ValueError, ZeroDivisionError):
        return None


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
    return min(
        1.0 - 2.0 * mass / radius + momentum_constant**2 / radius**4
        for radius in candidates
    )


def sgbl_radial_constraint_rk4_step(
    rhs: Callable[[float, Sequence[float]], tuple[float, float]],
    radius: float,
    state: tuple[float, float],
    step: float,
) -> tuple[float, float]:
    """One radial RK4 *constraint-integrator* step in ``(lambda, k)``."""

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


def sgbl_radial_constraint_ssprk3_step(
    rhs: Callable[[float, Sequence[float]], tuple[float, float]],
    radius: float,
    state: tuple[float, float],
    step: float,
) -> tuple[float, float]:
    """One radial SSPRK3 *constraint-integrator* step in ``(lambda, k)``."""

    first_rhs = rhs(radius, state)
    first = (state[0] + step * first_rhs[0], state[1] + step * first_rhs[1])
    second_rhs = rhs(radius + step, first)
    second = (
        0.75 * state[0] + 0.25 * (first[0] + step * second_rhs[0]),
        0.75 * state[1] + 0.25 * (first[1] + step * second_rhs[1]),
    )
    third_rhs = rhs(radius + step / 2.0, second)
    return (
        state[0] / 3.0 + 2.0 * (second[0] + step * third_rhs[0]) / 3.0,
        state[1] / 3.0 + 2.0 * (second[1] + step * third_rhs[1]) / 3.0,
    )


CONSTRAINT_INTEGRATORS = {
    "RK4": sgbl_radial_constraint_rk4_step,
    "SSPRK3": sgbl_radial_constraint_ssprk3_step,
}


@dataclass(frozen=True, slots=True)
class SGBLVacuumConstraintRefinement:
    """Known GR-vacuum constraint ODE versus the analytic exterior formula.

    This checks radial constraint-integrator order on a nontrivial exact
    vacuum solution.  It is not a time method, a family sweep, or a static
    SGB-L spacetime.
    """

    method: str
    step_count: int
    start_radius: float
    end_radius: float
    mass: float
    momentum_constant: float
    analytic_end_radial_metric: float
    analytic_end_k: float
    integrated_end_radial_metric: float
    integrated_end_k: float

    def __post_init__(self) -> None:
        if self.method not in CONSTRAINT_INTEGRATOR_NAMES:
            raise ValueError("unknown SGB-L radial constraint integrator")
        for name in (
            "start_radius",
            "end_radius",
            "mass",
            "momentum_constant",
            "analytic_end_radial_metric",
            "analytic_end_k",
            "integrated_end_radial_metric",
            "integrated_end_k",
        ):
            object.__setattr__(self, name, _finite(name, getattr(self, name)))
        if isinstance(self.step_count, bool) or not isinstance(self.step_count, int) or self.step_count <= 0:
            raise ValueError("step_count must be a positive integer")
        if self.step_count > MAX_RADIAL_CONSTRAINT_STEPS:
            raise SGBLConstraintFamilyStop(
                "resource_limit",
                "SGB-L family mesh exceeds the small fixed constraint-integrator budget",
                {"step_count": self.step_count, "maximum": MAX_RADIAL_CONSTRAINT_STEPS},
            )
        if self.start_radius <= 0.0 or self.end_radius <= self.start_radius:
            raise ValueError("vacuum refinement interval must be positive and increasing")
        expected_lambda, expected_k, _denominator = sgbl_floating_exterior_vacuum_evaluation(
            self.end_radius,
            mass=self.mass,
            momentum_constant=self.momentum_constant,
        )
        if self.analytic_end_radial_metric != expected_lambda or self.analytic_end_k != expected_k:
            raise ValueError("analytic end state does not match the exterior vacuum formula")
        if self.analytic_end_radial_metric <= 0.0 or self.integrated_end_radial_metric <= 0.0:
            raise ValueError("vacuum refinement radial metric must stay positive")

    @property
    def radial_metric_error(self) -> float:
        return abs(self.integrated_end_radial_metric - self.analytic_end_radial_metric)

    @property
    def k_error(self) -> float:
        return abs(self.integrated_end_k - self.analytic_end_k)

    @property
    def is_full_static_sgbl_solution(self) -> bool:
        return False

    @property
    def not_a_time_evolution_method(self) -> bool:
        return True


def sgbl_vacuum_constraint_rhs(
    radius: Real,
    state: Sequence[Real],
) -> tuple[float, float]:
    """Affine H/M right-hand side with identically vanishing scalars."""

    r = _finite("radius", radius, positive=True)
    L, k = _finite_pair("state", state)
    if L <= 0.0:
        raise SGBLConstraintFamilyStop(
            "nonpositive_metric",
            "radial metric left the positive polar-areal branch",
            {"radius": r, "radial_metric": L},
        )
    coefficients = sgbl_affine_constraint_coefficients(
        radius=r,
        radial_metric=L,
        angular_extrinsic_curvature=k,
        phi=0.0,
        phi_r=0.0,
        phi_rr=0.0,
        phi_pi=0.0,
        phi_pi_r=0.0,
        chi_r=0.0,
        chi_pi=0.0,
        planck_mass=DECLARED_PLANCK_MASS,
        scalar_mass=DECLARED_SCALAR_MASS,
        quartic_coupling=DECLARED_QUARTIC_COUPLING,
        alpha_gb=DECLARED_ALPHA_GB,
    )
    derivative = coefficients.derivatives()
    return float(derivative[0]), float(derivative[1])


def integrate_sgbl_vacuum_constraint_refinement(
    *,
    mass: Real,
    momentum_constant: Real,
    start_radius: Real,
    end_radius: Real,
    step_count: int,
    method: str = "RK4",
) -> SGBLVacuumConstraintRefinement:
    """Integrate a known vacuum constraint solution and compare to the analytic formula."""

    start = _finite("start_radius", start_radius, positive=True)
    end = _finite("end_radius", end_radius, positive=True)
    mass_value = _finite("mass", mass)
    momentum = _finite("momentum_constant", momentum_constant)
    if end <= start:
        raise ValueError("vacuum refinement interval must be increasing")
    if isinstance(step_count, bool) or not isinstance(step_count, int) or step_count <= 0:
        raise ValueError("step_count must be a positive integer")
    if step_count > MAX_RADIAL_CONSTRAINT_STEPS:
        raise SGBLConstraintFamilyStop(
            "resource_limit",
            "SGB-L family mesh exceeds the small fixed constraint-integrator budget",
            {"step_count": step_count, "maximum": MAX_RADIAL_CONSTRAINT_STEPS},
        )
    if method not in CONSTRAINT_INTEGRATORS:
        raise ValueError("unknown SGB-L radial constraint integrator")
    start_lambda, start_k, _start_den = sgbl_floating_exterior_vacuum_evaluation(
        start, mass=mass_value, momentum_constant=momentum
    )
    end_lambda, end_k, _end_den = sgbl_floating_exterior_vacuum_evaluation(
        end, mass=mass_value, momentum_constant=momentum
    )
    step = (end - start) / step_count
    state = (start_lambda, start_k)
    radius = start
    for _ in range(step_count):
        state = CONSTRAINT_INTEGRATORS[method](sgbl_vacuum_constraint_rhs, radius, state, step)
        radius += step
        if state[0] <= 0.0 or not all(isfinite(value) for value in state):
            raise SGBLConstraintFamilyStop(
                "nonpositive_metric" if state[0] <= 0.0 else "nonfinite_continuation",
                "vacuum constraint continuation left the regular branch",
                {"radius": radius},
            )
    return SGBLVacuumConstraintRefinement(
        method=method,
        step_count=step_count,
        start_radius=start,
        end_radius=end,
        mass=mass_value,
        momentum_constant=momentum,
        analytic_end_radial_metric=end_lambda,
        analytic_end_k=end_k,
        integrated_end_radial_metric=state[0],
        integrated_end_k=state[1],
    )


def _point_from_state(
    radius: float,
    state: tuple[float, float],
    parameters: SGBLFamilyParameters,
) -> SGBLFamilyPoint:
    derivative, coefficients, diagnostics = sgbl_constraint_rhs(radius, state, parameters)
    compactness = sgbl_misner_sharp_compactness(radius, state[0], state[1])
    residual = coefficients.evaluate(
        radial_metric_derivative=derivative[0],
        angular_extrinsic_curvature_derivative=derivative[1],
    )
    return SGBLFamilyPoint(
        radius=radius,
        radial_metric=state[0],
        angular_extrinsic_curvature=state[1],
        radial_metric_derivative=derivative[0],
        angular_extrinsic_curvature_derivative=derivative[1],
        hamiltonian_residual=float(residual[0]),
        momentum_residual=float(residual[1]),
        compactness=compactness,
        misner_sharp_mass=radius * compactness / 2.0,
        chi=diagnostics["chi"],
        chi_pi=diagnostics["chi_pi"],
        phi=diagnostics["phi"],
        phi_pi=diagnostics["phi_pi"],
    )


def integrate_sgbl_radial_constraints(
    parameters: SGBLFamilyParameters,
    *,
    step_count: int,
    method: str = "RK4",
) -> SGBLFamilySolution:
    """Integrate the two SGB-L constraint ODEs across the compact support.

    The stepping methods are radial constraint integrators.  They are not a
    time method.  Meshes larger than ``MAX_RADIAL_CONSTRAINT_STEPS`` are a
    typed resource stop: this slice does not run a production-grid sweep.
    """

    if not isinstance(parameters, SGBLFamilyParameters):
        raise TypeError("parameters must be SGBLFamilyParameters")
    if isinstance(step_count, bool) or not isinstance(step_count, int) or step_count <= 0:
        raise ValueError("step_count must be a positive integer")
    if step_count > MAX_RADIAL_CONSTRAINT_STEPS:
        raise SGBLConstraintFamilyStop(
            "resource_limit",
            "SGB-L family mesh exceeds the small fixed constraint-integrator budget",
            {"step_count": step_count, "maximum": MAX_RADIAL_CONSTRAINT_STEPS},
        )
    if method not in CONSTRAINT_INTEGRATORS:
        raise ValueError("unknown SGB-L radial constraint integrator")
    start = parameters.support_minimum
    end = parameters.support_maximum
    step = (end - start) / step_count
    state = (1.0, 0.0)

    def rhs(radius: float, values: Sequence[float]) -> tuple[float, float]:
        derivative, _coefficients, _diagnostics = sgbl_constraint_rhs(
            radius, values, parameters
        )
        return derivative

    points = [_point_from_state(start, state, parameters)]
    radius = start
    for _ in range(step_count):
        try:
            state = CONSTRAINT_INTEGRATORS[method](rhs, radius, state, step)
        except SGBLConstraintFamilyStop:
            raise
        except (ArithmeticError, FloatingPointError, ValueError, ZeroDivisionError) as exc:
            raise SGBLConstraintFamilyStop(
                "nonfinite_continuation",
                "SGB-L radial constraint continuation failed",
                {"radius": radius, "error": str(exc)},
            ) from exc
        radius += step
        if state[0] <= 0.0:
            raise SGBLConstraintFamilyStop(
                "nonpositive_metric",
                "radial metric became nonpositive during continuation",
                {"radius": radius, "radial_metric": state[0]},
            )
        if not all(isfinite(value) for value in state):
            raise SGBLConstraintFamilyStop(
                "nonfinite_continuation",
                "constraint continuation produced a nonfinite state",
                {"radius": radius},
            )
        points.append(_point_from_state(radius, state, parameters))

    return SGBLFamilySolution(
        method=method,
        step_count=step_count,
        parameters=parameters,
        points=tuple(points),
    )


def sgbl_reconstruct_constraints_by_finite_difference(
    solution: SGBLFamilySolution,
) -> tuple[dict[str, float], ...]:
    """Reconstruct ``(H, M)`` from stored ``lambda(r), k(r)`` by central differences.

    This is independent of the defining constraint-integrator right-hand
    side.  Vanishing defining-RHS residuals are not this check.
    """

    if not isinstance(solution, SGBLFamilySolution):
        raise TypeError("solution must be SGBLFamilySolution")
    points = solution.points
    if len(points) < 3:
        raise ValueError("finite-difference reconstruction needs at least three radial nodes")
    records: list[dict[str, float]] = []
    for index in range(1, len(points) - 1):
        previous = points[index - 1]
        current = points[index]
        nxt = points[index + 1]
        step_left = current.radius - previous.radius
        step_right = nxt.radius - current.radius
        if step_left <= 0.0 or step_right <= 0.0:
            raise ValueError("finite-difference reconstruction requires increasing radii")
        if step_left != step_right:
            raise ValueError("finite-difference reconstruction requires a uniform mesh")
        lambda_r = (nxt.radial_metric - previous.radial_metric) / (2.0 * step_left)
        k_r = (
            nxt.angular_extrinsic_curvature - previous.angular_extrinsic_curvature
        ) / (2.0 * step_left)
        fields = sgbl_compact_family_fields(current.radius, solution.parameters)
        coefficients = sgbl_affine_constraint_coefficients(
            radius=current.radius,
            radial_metric=current.radial_metric,
            angular_extrinsic_curvature=current.angular_extrinsic_curvature,
            phi=fields["phi"],
            phi_r=fields["phi_r"],
            phi_rr=fields["phi_rr"],
            phi_pi=fields["phi_pi"],
            phi_pi_r=fields["phi_pi_r"],
            chi_r=fields["chi_r"],
            chi_pi=fields["chi_pi"],
            planck_mass=solution.parameters.planck_mass,
            scalar_mass=solution.parameters.scalar_mass,
            quartic_coupling=solution.parameters.quartic_coupling,
            alpha_gb=solution.parameters.alpha_gb,
        )
        hamiltonian, momentum = coefficients.evaluate(
            radial_metric_derivative=lambda_r,
            angular_extrinsic_curvature_derivative=k_r,
        )
        records.append(
            {
                "radius": current.radius,
                "finite_difference_lambda_r": lambda_r,
                "finite_difference_k_r": k_r,
                "defining_lambda_r": current.radial_metric_derivative,
                "defining_k_r": current.angular_extrinsic_curvature_derivative,
                "finite_difference_H": float(hamiltonian),
                "finite_difference_M": float(momentum),
                "defining_H": current.hamiltonian_residual,
                "defining_M": current.momentum_residual,
                "chi_pi": fields["chi_pi"],
                "chi": fields["chi"],
            }
        )
    return tuple(records)


def sgbl_branch_health_gate() -> dict[str, Any]:
    """Holdout health remains closed; that closure is not a hyperbolicity proof."""

    return {
        "SGBL_branch_owned_and_healthy": False,
        "missing_health_closes_the_gate": True,
        "missing_health_proves_nonhyperbolicity": False,
        "interval_invertibility": "unqualified",
        "source_runtime": "unqualified",
        "kinetic_cone_health": "unqualified",
        "evolving_center_after_matter_arrives": "unqualified",
        "continuous_no_initial_trapped_sphere": "unqualified",
        "sampled_nodes_untrapped_is_not_continuous_no_initial_trap": True,
    }


def sgbl_tensor_constraint_pair(state: SphericalState) -> tuple[Any, Any]:
    """Independent RED1 normal projections ``(H, M)``."""

    if not isinstance(state, SphericalState):
        raise TypeError("state must be a SphericalState")
    projection = physical_constraint_projections(state)
    return projection["H"], projection["M"]


def sgbl_scalar_source_on_vacuum_geometry(state: SphericalState) -> Any:
    """Unredefined scalar residual; nonzero GB shows the outer continuation is not static SGB-L."""

    if not isinstance(state, SphericalState):
        raise TypeError("state must be a SphericalState")
    return residuals(state)["phi"]
