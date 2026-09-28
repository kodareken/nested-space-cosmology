"""Validated continuous initial-slice compactness and evolving-centre monitor.

The matched SGB-L family is the continuum object with an exact Minkowski
buffer, compact-support scalar cells, and the analytic GR vacuum exterior
``k=J/r^3``, ``lambda^{-2}=1-2M/r+J^2/r^4``.  Sampled constraint-integrator
nodes are not a continuous no-trap certificate.

Compact-support geometry is enclosed in the Misner-Sharp chart ``(C,k)``
with ``C=1+r^2 k^2-lambda^{-2}`` and ``D=1+r^2 k^2-C=lambda^{-2}``.
Reconstruction uses the positive square-root branch ``lambda=D^{-1/2}``
with ``D`` strictly positive.  The affine H/M solve for ``lambda_r`` and
``k_r`` is unchanged.  The state ODEs are ``k_r`` and
``C_r=2 r k^2 + 2 r^2 k k_r + 2 lambda_r/lambda^3``, started from exact
``(C,k)=(0,0)``.  Picard self-maps in ``(C,k)``.  The product-box
``(lambda,k,C)`` graph and the ``(C,J)`` chart are regression only.  The
exterior identity ``C=2M/r`` is read from the enclosed support-boundary
state.

The evolving-centre monitor is a fail-closed contract grounded in
:mod:`sgb1_ctl1_center` and :mod:`regular_center`.  An unexecuted monitor is
not a pre-holdout pass.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import isqrt
from numbers import Integral
from types import MappingProxyType
from typing import Any, Mapping

from .exact_interval import Interval, interval
from .interval_tangent import IntervalFirstTangent, primal_and_tangent
from .sgb1_ctl1_center import (
    sgbl_empty_minkowski_center_profile,
    sgbl_initial_center_series,
    validate_sgbl_initial_center_profile,
)
from .sgb1_ctl1_constraints import SGBLAnnularInputs, sgbl_constraint_coefficients
from .sgb1_ctl1_family import (
    CONSTRAINT_INTEGRATOR_NAMES,
    DECLARED_ALPHA_GB,
    DECLARED_CENTER,
    DECLARED_CHI_AMPLITUDE,
    DECLARED_HALF_WIDTH,
    DECLARED_OUTER_RADIUS,
    DECLARED_PHI_AMPLITUDE,
    DECLARED_PLANCK_MASS,
    DECLARED_QUARTIC_COUPLING,
    DECLARED_SCALAR_MASS,
    SGBLFamilyParameters,
    integrate_sgbl_radial_constraints,
)


Q = Fraction
BUMP_X_DERIVATIVE_BOUND = Q(3)
BUMP_XX_DERIVATIVE_BOUND = Q(80)
FACTORIAL_EIGHT = 40320
DECLARED_LAMBDA_DOMAIN = Interval(Q(1, 8), Q(16))
DECLARED_K_DOMAIN = Interval(Q(-4), Q(4))
DECLARED_C_DOMAIN = Interval(Q(-16), Q(16))
DECLARED_J_DOMAIN = Interval(-Q(16384), Q(16384))
CERTIFICATE_CHART_KIND = "c_k"
DECLARED_BASE_CELLS = 16
DECLARED_MAX_BISECTION_DEPTH = 8
DECLARED_MAX_CELLS = 1024
DECLARED_MAX_RHS_EVALUATIONS = 8192
DECLARED_MAX_PICARD_ITERATIONS = 12
DECLARED_MAX_RATIONAL_BITS = 16384
DECLARED_EULER_INFLATION = 2
DECLARED_DYADIC_DENOMINATOR = 2 ** 64
MONITOR_REQUIRED_CHECKS = (
    "even_parity",
    "elementary_flatness_A_equals_lambda_value_dt_dtt",
    "certified_negative_powers_absent",
    "emptiness_not_assumed_after_matter_arrives",
)
HEALTH_STOP_REASONS = frozenset(
    {
        "domain_error",
        "resource_limit",
        "interval_inconclusive",
        "monitor_not_executed",
        "missing_trajectory_owner",
    }
)
ODE_INCONCLUSIVE_REASONS = frozenset(
    {
        "picard_strict_self_map_failed",
        "jacobian_diagonal_contains_zero",
        "physical_domain_escape",
        "compactness_enclosure_not_below_one",
        "residual_enclosure_misses_origin",
        "exterior_compactness_not_below_one",
        "zero_in_interval_reciprocal",
        "mean_value_inconsistent",
        "compactness_mean_value_inconsistent",
        "algebraic_compactness_mean_value_inconsistent",
        "compactness_invariant_misses_origin",
        "chart_denominator_not_strictly_positive",
        "sqrt_branch_not_positive",
        "chart_roundtrip_misses_origin",
    }
)


def _fraction(name: str, value: object, *, nonnegative: bool = False) -> Fraction:
    if type(value) not in (int, Fraction):
        raise TypeError(f"{name} must be a Fraction or built-in int")
    result = Fraction(value)
    if nonnegative and result < 0:
        raise ValueError(f"{name} must be nonnegative")
    return result


def _positive_int(name: str, value: object, *, allow_zero: bool = False) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral):
        raise ValueError(f"{name} must be an integer")
    result = int(value)
    if result < 0 or (result == 0 and not allow_zero):
        raise ValueError(f"{name} must be a positive integer")
    return result


def _float_to_exact(name: str, value: float) -> Fraction:
    result = Fraction(value)
    if float(result) != float(value):
        raise SGBLInitialHealthStop(
            "domain_error",
            f"{name} is not an exact dyadic family parameter",
            {"value": value},
        )
    return result


def _as_interval(value: object) -> Interval:
    if type(value) is Interval:
        return value
    if type(value) in (int, Fraction):
        return Interval.singleton(Fraction(value))
    raise TypeError("value must be an Interval or exact rational")


def _intersect(left: Interval, right: Interval) -> Interval | None:
    lower = max(left.lower, right.lower)
    upper = min(left.upper, right.upper)
    if lower > upper:
        return None
    return Interval(lower, upper)


def _interval_bits(value: Interval) -> int:
    return max(
        abs(value.lower.numerator).bit_length(),
        value.lower.denominator.bit_length(),
        abs(value.upper.numerator).bit_length(),
        value.upper.denominator.bit_length(),
    )


def _floor_dyadic(value: Fraction, denominator: int) -> Fraction:
    return Fraction((value.numerator * denominator) // value.denominator, denominator)


def _ceil_dyadic(value: Fraction, denominator: int) -> Fraction:
    return -_floor_dyadic(-value, denominator)


def _outward(value: Interval, denominator: int = DECLARED_DYADIC_DENOMINATOR) -> Interval:
    """Valid outward rounding onto a declared dyadic grid.  This is a resource."""

    return Interval(_floor_dyadic(value.lower, denominator), _ceil_dyadic(value.upper, denominator))


def _perfect_square_sqrt(value: Fraction) -> Fraction | None:
    if value < 0:
        return None
    numerator_root = isqrt(value.numerator)
    denominator_root = isqrt(value.denominator)
    if (
        numerator_root * numerator_root == value.numerator
        and denominator_root * denominator_root == value.denominator
    ):
        return Fraction(numerator_root, denominator_root)
    return None


def _sqrt_floor_dyadic(
    value: Fraction, denominator: int = DECLARED_DYADIC_DENOMINATOR
) -> Fraction:
    if value < 0:
        raise SGBLInitialHealthStop("domain_error", "square-root lower bound requires a nonnegative argument")
    if value == 0:
        return Q(0)
    scaled = value.numerator * denominator * denominator
    return Fraction(isqrt(scaled // value.denominator), denominator)


def _sqrt_ceil_dyadic(
    value: Fraction, denominator: int = DECLARED_DYADIC_DENOMINATOR
) -> Fraction:
    if value < 0:
        raise SGBLInitialHealthStop("domain_error", "square-root upper bound requires a nonnegative argument")
    if value == 0:
        return Q(0)
    scaled = value.numerator * denominator * denominator
    root = isqrt(scaled // value.denominator)
    if root * root * value.denominator < scaled:
        root += 1
    return Fraction(root, denominator)


def sgbl_interval_sqrt(value: Interval) -> Interval:
    """Positive outward square-root enclosure.  Zero or sign change is fail-closed."""

    value = _as_interval(value)
    if not value.strictly_positive():
        raise SGBLInitialHealthStop(
            "interval_inconclusive",
            "square-root chart coordinate is not strictly positive",
            {"obstruction": "chart_denominator_not_strictly_positive"},
        )
    exact_lower = _perfect_square_sqrt(value.lower)
    exact_upper = _perfect_square_sqrt(value.upper)
    lower = exact_lower if exact_lower is not None else _sqrt_floor_dyadic(value.lower)
    upper = exact_upper if exact_upper is not None else _sqrt_ceil_dyadic(value.upper)
    if lower > upper:
        raise SGBLInitialHealthStop(
            "interval_inconclusive",
            "outward square-root enclosure inverted",
            {"obstruction": "chart_denominator_not_strictly_positive"},
        )
    return Interval(lower, upper)


class SGBLInitialHealthStop(ArithmeticError):
    """Typed fail-closed stop of the continuous initial-health owner."""

    def __init__(
        self,
        reason: str,
        message: str,
        payload: Mapping[str, Any] | None = None,
    ) -> None:
        if reason not in HEALTH_STOP_REASONS:
            raise ValueError("unknown SGB-L initial-health stop reason")
        self.reason = reason
        self.payload = dict(payload or {})
        super().__init__(message)


@dataclass
class _ODEBudget:
    rhs_evaluations: int = 0
    compactness_rhs_evaluations: int = 0
    cells: int = 0
    max_bits: int = 0
    mean_value_strictly_tighter: bool = False
    compactness_mean_value_strictly_tighter: bool = False
    algebraic_mean_value_strictly_tighter: bool = False
    correlated_strictly_tighter: bool = False


@dataclass(frozen=True, slots=True, kw_only=True)
class SGBLValidatedODEPolicy:
    """Declared Picard/domain/refinement policy.  Not fitted after an outcome."""

    lambda_domain: Interval = DECLARED_LAMBDA_DOMAIN
    k_domain: Interval = DECLARED_K_DOMAIN
    compactness_domain: Interval = DECLARED_C_DOMAIN
    angular_momentum_domain: Interval = DECLARED_J_DOMAIN
    base_cells: int = DECLARED_BASE_CELLS
    max_bisection_depth: int = DECLARED_MAX_BISECTION_DEPTH
    max_cells: int = DECLARED_MAX_CELLS
    max_rhs_evaluations: int = DECLARED_MAX_RHS_EVALUATIONS
    max_picard_iterations: int = DECLARED_MAX_PICARD_ITERATIONS
    max_rational_bit_length: int = DECLARED_MAX_RATIONAL_BITS

    def __post_init__(self) -> None:
        if type(self.lambda_domain) is not Interval or type(self.k_domain) is not Interval:
            raise TypeError("Picard domains must be exact intervals")
        if type(self.compactness_domain) is not Interval:
            raise TypeError("Picard compactness domain must be an exact interval")
        if type(self.angular_momentum_domain) is not Interval:
            raise TypeError("Picard angular-momentum domain must be an exact interval")
        if not self.lambda_domain.strictly_positive():
            raise SGBLInitialHealthStop(
                "domain_error",
                "declared lambda domain must be strictly positive",
            )
        if self.compactness_domain.upper <= 1:
            raise SGBLInitialHealthStop(
                "domain_error",
                "declared compactness resource domain must not assume C<=1",
                {
                    "compactness_domain": (
                        self.compactness_domain.lower,
                        self.compactness_domain.upper,
                    )
                },
            )
        if not (self.compactness_domain.lower < 0 < self.compactness_domain.upper):
            raise SGBLInitialHealthStop(
                "domain_error",
                "declared compactness resource domain must contain a neighborhood of C=0",
                {
                    "compactness_domain": (
                        self.compactness_domain.lower,
                        self.compactness_domain.upper,
                    )
                },
            )
        if not (
            self.angular_momentum_domain.lower < 0 < self.angular_momentum_domain.upper
        ):
            raise SGBLInitialHealthStop(
                "domain_error",
                "declared angular-momentum resource domain must contain a neighborhood of J=0",
                {
                    "angular_momentum_domain": (
                        self.angular_momentum_domain.lower,
                        self.angular_momentum_domain.upper,
                    )
                },
            )
        object.__setattr__(self, "base_cells", _positive_int("base_cells", self.base_cells))
        object.__setattr__(
            self,
            "max_bisection_depth",
            _positive_int("max_bisection_depth", self.max_bisection_depth, allow_zero=True),
        )
        object.__setattr__(self, "max_cells", _positive_int("max_cells", self.max_cells))
        object.__setattr__(
            self,
            "max_rhs_evaluations",
            _positive_int("max_rhs_evaluations", self.max_rhs_evaluations),
        )
        object.__setattr__(
            self,
            "max_picard_iterations",
            _positive_int("max_picard_iterations", self.max_picard_iterations),
        )
        object.__setattr__(
            self,
            "max_rational_bit_length",
            _positive_int("max_rational_bit_length", self.max_rational_bit_length),
        )


DECLARED_ODE_POLICY = SGBLValidatedODEPolicy()


@dataclass(frozen=True, slots=True, kw_only=True)
class SGBLExactInitialSlice:
    """Exact declared matched-family geometry for the continuum certificate."""

    chi_amplitude: Fraction | int = 3
    phi_amplitude: Fraction | int = Q(1, 131072)
    center: Fraction | int = 12
    half_width: Fraction | int = 2
    outer_radius: Fraction | int = 128
    planck_mass: Fraction | int = 2
    scalar_mass: Fraction | int = 3
    quartic_coupling: Fraction | int = Q(1, 2)
    alpha_gb: Fraction | int = -Q(1, 4)

    def __post_init__(self) -> None:
        for name in (
            "chi_amplitude",
            "phi_amplitude",
            "center",
            "half_width",
            "outer_radius",
            "planck_mass",
            "scalar_mass",
            "quartic_coupling",
            "alpha_gb",
        ):
            object.__setattr__(
                self,
                name,
                _fraction(
                    name,
                    getattr(self, name),
                    nonnegative=name not in {"alpha_gb"},
                ),
            )
        if self.chi_amplitude <= 0 or self.half_width <= 0 or self.planck_mass <= 0:
            raise ValueError("chi amplitude, half-width and Planck mass must be positive")
        if self.center <= self.half_width:
            raise ValueError("compact profiles must retain a strict centre buffer")
        if self.outer_radius <= self.support_maximum:
            raise ValueError("outer radius must retain a strict vacuum buffer")

    @property
    def support_minimum(self) -> Fraction:
        return self.center - self.half_width

    @property
    def support_maximum(self) -> Fraction:
        return self.center + self.half_width

    @property
    def planck_mass_squared(self) -> Fraction:
        return self.planck_mass**2

    @classmethod
    def from_family_parameters(cls, parameters: SGBLFamilyParameters) -> "SGBLExactInitialSlice":
        if type(parameters) is not SGBLFamilyParameters:
            raise TypeError("parameters must be SGBLFamilyParameters")
        return cls(
            chi_amplitude=_float_to_exact("chi_amplitude", parameters.chi_amplitude),
            phi_amplitude=_float_to_exact("phi_amplitude", parameters.phi_amplitude),
            center=_float_to_exact("center", parameters.center),
            half_width=_float_to_exact("chi_half_width", parameters.chi_half_width),
            outer_radius=_float_to_exact("outer_radius", parameters.outer_radius),
            planck_mass=_float_to_exact("planck_mass", parameters.planck_mass),
            scalar_mass=_float_to_exact("scalar_mass", parameters.scalar_mass),
            quartic_coupling=_float_to_exact(
                "quartic_coupling", parameters.quartic_coupling
            ),
            alpha_gb=_float_to_exact("alpha_gb", parameters.alpha_gb),
        )


def _exp_neg_upper(value: Fraction) -> Fraction:
    """Return an exact rational upper bound for ``e^{-t}`` on ``t>=0``."""

    if value < 0:
        raise SGBLInitialHealthStop("domain_error", "exp(-t) upper bound requires t>=0")
    if value == 0:
        return Q(1)
    if value >= 8:
        return Fraction(FACTORIAL_EIGHT) / value**8
    total = Q(1)
    term = Q(1)
    for index in range(1, 17):
        term *= value / index
        total += term
    return 1 / total


def _exp_neg_lower(value: Fraction) -> Fraction:
    """Return a nonnegative exact rational lower bound for ``e^{-t}``."""

    if value < 0:
        raise SGBLInitialHealthStop("domain_error", "exp(-t) lower bound requires t>=0")
    if value == 0:
        return Q(1)
    if value >= 8:
        return Q(0)
    total = Q(1)
    term = Q(1)
    for index in range(1, 17):
        term *= value / index
        total += term
    exponential_upper = total + term * value / (17 - value)
    return 1 / exponential_upper


def _exp_neg_enclosure(value: Interval) -> Interval:
    if value.lower < 0:
        raise SGBLInitialHealthStop(
            "domain_error",
            "compact-bump exponent enclosure requires a nonnegative argument",
        )
    return Interval(_exp_neg_lower(value.upper), _exp_neg_upper(value.lower))


def _bump_global_enclosure(half_width: Fraction) -> tuple[Interval, Interval, Interval]:
    first = interval(-BUMP_X_DERIVATIVE_BOUND, BUMP_X_DERIVATIVE_BOUND) / half_width
    second = interval(-BUMP_XX_DERIVATIVE_BOUND, BUMP_XX_DERIVATIVE_BOUND) / (
        half_width * half_width
    )
    return interval(0, 1), first, second


def sgbl_compact_bump_enclosure(
    radius: Interval,
    *,
    center: Fraction,
    half_width: Fraction,
) -> tuple[Interval, Interval, Interval]:
    """Interval enclosure of the declared compact bump and two radial derivatives.

    Cells wholly inside ``(-1,1)`` in the scaled coordinate use that full
    interval, including ``(-1,-7/8)`` and ``(7/8,1)``.  Cells that touch or
    cross the support boundary use the global enclosure.  The bump is not
    clipped away from the edge.
    """

    if half_width <= 0:
        raise SGBLInitialHealthStop("domain_error", "bump half-width must be positive")
    if radius.upper < radius.lower:
        raise SGBLInitialHealthStop("domain_error", "radius interval is empty")
    scaled = (radius - center) / half_width
    if scaled.upper <= -1 or scaled.lower >= 1:
        zero = Interval.singleton(0)
        return zero, zero, zero
    if scaled.lower <= -1 or scaled.upper >= 1:
        return _bump_global_enclosure(half_width)
    square = scaled ** 2
    one_minus_square = Interval.singleton(1) - square
    if not one_minus_square.strictly_positive():
        return _bump_global_enclosure(half_width)
    peak = one_minus_square.reciprocal()
    exponent_shift = peak - Interval.singleton(1)
    bump = _exp_neg_enclosure(exponent_shift)
    first_exponent = (interval(-2, -2) * scaled) / (one_minus_square * one_minus_square)
    second_exponent = interval(-2, -2) / (one_minus_square * one_minus_square) + (
        interval(-8, -8) * square
    ) / (one_minus_square * one_minus_square * one_minus_square)
    first = bump * first_exponent / half_width
    second = bump * (second_exponent + first_exponent ** 2) / (half_width * half_width)
    global_value, global_first, global_second = _bump_global_enclosure(half_width)
    value = Interval(
        max(bump.lower, global_value.lower), min(bump.upper, global_value.upper)
    )
    first = Interval(max(first.lower, global_first.lower), min(first.upper, global_first.upper))
    second = Interval(
        max(second.lower, global_second.lower), min(second.upper, global_second.upper)
    )
    return value, first, second


def sgbl_interval_family_fields(
    radius: Interval,
    spec: SGBLExactInitialSlice,
) -> dict[str, Interval]:
    """Enclose the lambda-free compact profiles on a radius interval."""

    if not radius.strictly_positive():
        raise SGBLInitialHealthStop(
            "domain_error",
            "compact-support cells must lie at positive radius",
        )
    bump, bump_r, bump_rr = sgbl_compact_bump_enclosure(
        radius, center=spec.center, half_width=spec.half_width
    )
    inv_r = radius.reciprocal()
    chi_amp = Interval.singleton(spec.chi_amplitude)
    phi_amp = Interval.singleton(spec.phi_amplitude)
    return {
        "phi": _outward(phi_amp * bump),
        "phi_r": _outward(phi_amp * bump_r),
        "phi_rr": _outward(phi_amp * bump_rr),
        "phi_pi": Interval.singleton(0),
        "phi_pi_r": Interval.singleton(0),
        "chi_r": _outward(chi_amp * (bump_r * inv_r - bump * inv_r * inv_r)),
        "chi_pi": _outward(chi_amp * bump_r * inv_r),
    }


def _primal(value: object) -> Interval:
    if isinstance(value, IntervalFirstTangent):
        return value.primal
    return _as_interval(value)


def sgbl_generic_affine_coefficients(
    *,
    radius: Interval,
    radial_metric: Interval | IntervalFirstTangent,
    angular_extrinsic_curvature: Interval | IntervalFirstTangent,
    spec: SGBLExactInitialSlice,
) -> dict[str, Interval | IntervalFirstTangent]:
    """Affine H/M coefficients on Interval or IntervalFirstTangent state.

    The algebra is the generic annular kernel.  Seeding ``lambda`` or ``k``
    with an interval first tangent yields an enclosure of the corresponding
    directional derivative.
    """

    if type(radius) is not Interval:
        raise TypeError("radius must be an Interval")
    if not radius.strictly_positive() or not _primal(radial_metric).strictly_positive():
        raise SGBLInitialHealthStop(
            "domain_error",
            "interval affine coefficients require positive r and lambda",
        )
    fields = sgbl_interval_family_fields(radius, spec)
    mpl2 = Interval.singleton(spec.planck_mass_squared)
    alpha = Interval.singleton(spec.alpha_gb)
    two = Interval.singleton(2)
    eight = Interval.singleton(8)
    sixteen = Interval.singleton(16)
    lam = radial_metric
    k = angular_extrinsic_curvature
    a0 = Interval.singleton(-2) * (k ** 2)
    a_l = (lam ** 3 * radius).reciprocal()
    b = (k ** 2) + (Interval.singleton(1) - (lam ** 2).reciprocal()) / (radius ** 2)
    c0 = Interval.singleton(3) * k / (lam * radius)
    c_k = lam.reciprocal()
    x0 = fields["phi_rr"] / (lam ** 2) - two * k * fields["phi_pi"]
    x_l = -fields["phi_r"] / (lam ** 3)
    y = k * fields["phi_pi"] + fields["phi_r"] / (lam ** 2 * radius)
    z = (fields["phi_pi_r"] - two * k * fields["phi_r"]) / lam
    rho = (
        (fields["phi_pi"] ** 2 + fields["chi_pi"] ** 2) / two
        + (fields["phi_r"] ** 2 + fields["chi_r"] ** 2) / (two * (lam ** 2))
        + Interval.singleton(spec.scalar_mass ** 2) * (fields["phi"] ** 2) / two
        + Interval.singleton(spec.quartic_coupling) * (fields["phi"] ** 4) / Interval.singleton(4)
    )
    return {
        "hamiltonian_constant": mpl2 * (two * a0 + b) - rho - eight * alpha * (b * x0 + two * a0 * y),
        "hamiltonian_lambda_r": two * mpl2 * a_l - eight * alpha * (b * x_l + two * a_l * y),
        "momentum_constant": (
            two * mpl2 * lam * c0
            - fields["phi_pi"] * fields["phi_r"]
            - fields["chi_pi"] * fields["chi_r"]
            - eight * alpha * lam * (b * z + two * c0 * y)
        ),
        "momentum_k_r": two * mpl2 * lam * c_k - sixteen * alpha * lam * c_k * y,
    }


def sgbl_interval_affine_coefficients(
    *,
    radius: Interval,
    radial_metric: Interval,
    angular_extrinsic_curvature: Interval,
    spec: SGBLExactInitialSlice,
) -> dict[str, Interval]:
    """Natural interval extension of the generic affine H/M coefficients."""

    coefficients = sgbl_generic_affine_coefficients(
        radius=radius,
        radial_metric=radial_metric,
        angular_extrinsic_curvature=angular_extrinsic_curvature,
        spec=spec,
    )
    return {name: _primal(value) for name, value in coefficients.items()}


def sgbl_misner_sharp_compactness_interval(
    radius: Interval,
    radial_metric: Interval,
    angular_extrinsic_curvature: Interval,
) -> Interval:
    """Direct interval compactness ``C=1+r^2 k^2-lambda^{-2}``."""

    if not radius.strictly_positive() or not radial_metric.strictly_positive():
        raise SGBLInitialHealthStop(
            "domain_error",
            "compactness enclosure requires positive r and lambda",
        )
    return (
        Interval.singleton(1)
        + (radius ** 2) * (angular_extrinsic_curvature ** 2)
        - (radial_metric ** 2).reciprocal()
    )


def sgbl_misner_sharp_compactness_state_gradient(
    radius: Interval,
    radial_metric: Interval,
    angular_extrinsic_curvature: Interval,
) -> dict[str, Interval]:
    """Exact partials of algebraic ``C(r,lambda,k)=1+r^2 k^2-lambda^{-2}``.

    Holding the other two arguments fixed:

    ``dC/dr=2 r k^2``, ``dC/dlambda=2/lambda^3``, ``dC/dk=2 r^2 k``.
    These are not the ODE chain-rule ``C_r``, which also includes
    ``lambda_r`` and ``k_r``.
    """

    if not radius.strictly_positive() or not radial_metric.strictly_positive():
        raise SGBLInitialHealthStop(
            "domain_error",
            "compactness gradient requires positive r and lambda",
        )
    return {
        "d_c_d_r": 2 * radius * (angular_extrinsic_curvature ** 2),
        "d_c_d_lambda": 2 * (radial_metric ** 3).reciprocal(),
        "d_c_d_k": 2 * (radius ** 2) * angular_extrinsic_curvature,
    }


def sgbl_centered_mean_value_algebraic_compactness(
    radius: Interval,
    radial_metric: Interval,
    angular_extrinsic_curvature: Interval,
) -> Interval:
    """Centered multivariable mean-value enclosure of algebraic ``C``.

    ``C(Y) subset C(mid)+dC/dr(Y)(r-r_mid)+dC/dlambda(Y)(lambda-lambda_mid)
    +dC/dk(Y)(k-k_mid)`` with the exact partials evaluated on the full box.
    """

    radius = _as_interval(radius)
    radial_metric = _as_interval(radial_metric)
    angular_extrinsic_curvature = _as_interval(angular_extrinsic_curvature)
    r_mid = Interval.singleton(radius.midpoint())
    lambda_mid = Interval.singleton(radial_metric.midpoint())
    k_mid = Interval.singleton(angular_extrinsic_curvature.midpoint())
    center = sgbl_misner_sharp_compactness_interval(
        r_mid, lambda_mid, k_mid
    )
    gradient = sgbl_misner_sharp_compactness_state_gradient(
        radius, radial_metric, angular_extrinsic_curvature
    )
    return (
        center
        + gradient["d_c_d_r"] * (radius - r_mid)
        + gradient["d_c_d_lambda"] * (radial_metric - lambda_mid)
        + gradient["d_c_d_k"] * (angular_extrinsic_curvature - k_mid)
    )


def sgbl_intersect_algebraic_compactness_enclosures(
    radius: Interval,
    radial_metric: Interval,
    angular_extrinsic_curvature: Interval,
    *,
    propagated: Interval | None = None,
) -> dict[str, Interval]:
    """Intersect natural algebraic ``C``, centered mean-value ``C``, and optional propagated ``C``."""

    natural = sgbl_misner_sharp_compactness_interval(
        radius, radial_metric, angular_extrinsic_curvature
    )
    centered = sgbl_centered_mean_value_algebraic_compactness(
        radius, radial_metric, angular_extrinsic_curvature
    )
    algebraic = _intersect(natural, centered)
    if algebraic is None:
        raise SGBLInitialHealthStop(
            "interval_inconclusive",
            "natural and centered algebraic compactness enclosures are disjoint",
            {"obstruction": "algebraic_compactness_mean_value_inconsistent"},
        )
    correlated = algebraic
    if propagated is not None:
        if type(propagated) is not Interval:
            raise TypeError("propagated compactness must be an Interval")
        correlated = _intersect(algebraic, propagated)
        if correlated is None:
            raise SGBLInitialHealthStop(
                "interval_inconclusive",
                "algebraic and propagated compactness enclosures are disjoint",
                {"obstruction": "compactness_invariant_misses_origin"},
            )
    return {
        "natural": natural,
        "centered": centered,
        "algebraic": algebraic,
        "correlated": correlated,
    }


def sgbl_misner_sharp_compactness_radial_derivative(
    radius: Interval | IntervalFirstTangent,
    radial_metric: Interval | IntervalFirstTangent,
    angular_extrinsic_curvature: Interval | IntervalFirstTangent,
    lambda_r: Interval | IntervalFirstTangent,
    k_r: Interval | IntervalFirstTangent,
) -> Interval | IntervalFirstTangent:
    """Exact chain-rule ``C_r=2 r k^2 + 2 r^2 k k_r + 2 lambda_r/lambda^3``.

    From ``C=1+r^2 k^2-lambda^{-2}``:

    ``d/dr(r^2 k^2)=2 r k^2 + 2 r^2 k k_r``,
    ``d/dr(-lambda^{-2})=2 lambda_r / lambda^3``.
    """

    if not _primal(radius).strictly_positive() or not _primal(radial_metric).strictly_positive():
        raise SGBLInitialHealthStop(
            "domain_error",
            "compactness derivative requires positive r and lambda",
        )
    return (
        2 * radius * (angular_extrinsic_curvature ** 2)
        + 2 * (radius ** 2) * angular_extrinsic_curvature * k_r
        + 2 * lambda_r / (radial_metric ** 3)
    )


def sgbl_misner_sharp_compactness_invariant_residual(
    radius: Interval,
    radial_metric: Interval,
    angular_extrinsic_curvature: Interval,
    compactness: Interval,
) -> Interval:
    """Enclosure of ``C-(1+r^2 k^2-lambda^{-2})``.  The identity requires zero."""

    return compactness - sgbl_misner_sharp_compactness_interval(
        radius, radial_metric, angular_extrinsic_curvature
    )


def sgbl_misner_sharp_compactness_chain_rule_identity(
    radius: Interval,
    radial_metric: Interval,
    angular_extrinsic_curvature: Interval,
    lambda_r: Interval,
    k_r: Interval,
) -> dict[str, Interval]:
    """Checked first-order identity: the r-jet of ``C`` matches ``C_r``.

    Seeding ``r`` with tangent 1 and ``(lambda,k)`` with tangents
    ``(lambda_r,k_r)`` differentiates ``C=1+r^2 k^2-lambda^{-2}`` by the
    interval first-tangent product rule.  The residual against the closed
    form must contain zero; at singletons it is exactly zero.
    """

    radius = _as_interval(radius)
    radial_metric = _as_interval(radial_metric)
    angular_extrinsic_curvature = _as_interval(angular_extrinsic_curvature)
    lambda_r = _as_interval(lambda_r)
    k_r = _as_interval(k_r)
    r_jet = IntervalFirstTangent.seed(radius)
    lam_jet = IntervalFirstTangent(radial_metric, lambda_r)
    k_jet = IntervalFirstTangent(angular_extrinsic_curvature, k_r)
    compactness_jet = (
        IntervalFirstTangent.constant(1)
        + (r_jet ** 2) * (k_jet ** 2)
        - (lam_jet ** 2).reciprocal()
    )
    formula = sgbl_misner_sharp_compactness_radial_derivative(
        radius, radial_metric, angular_extrinsic_curvature, lambda_r, k_r
    )
    formula_primal = _primal(formula)
    algebraic = sgbl_misner_sharp_compactness_interval(
        radius, radial_metric, angular_extrinsic_curvature
    )
    return {
        "jet_primal": compactness_jet.primal,
        "jet_tangent": compactness_jet.tangent,
        "formula": formula_primal,
        "algebraic_compactness": algebraic,
        "invariant_residual": compactness_jet.primal - algebraic,
        "chain_rule_residual": compactness_jet.tangent - formula_primal,
    }


def sgbl_misner_sharp_angular_momentum_interval(
    radius: Interval,
    angular_extrinsic_curvature: Interval,
) -> Interval:
    """Chart coordinate ``J=r^3 k``."""

    if not radius.strictly_positive():
        raise SGBLInitialHealthStop(
            "domain_error",
            "angular-momentum chart requires positive radius",
        )
    return (radius ** 3) * angular_extrinsic_curvature


def sgbl_misner_sharp_angular_momentum_radial_derivative(
    radius: Interval,
    angular_momentum: Interval,
    k_r: Interval,
) -> Interval:
    """Exact ``J_r=3 J/r + r^3 k_r`` from ``J=r^3 k``."""

    if not radius.strictly_positive():
        raise SGBLInitialHealthStop(
            "domain_error",
            "angular-momentum derivative requires positive radius",
        )
    return 3 * angular_momentum / radius + (radius ** 3) * k_r


def sgbl_misner_sharp_chart_state(
    radius: Interval,
    compactness: Interval,
    angular_momentum: Interval,
) -> dict[str, Interval]:
    """Reconstruct ``k=J/r^3`` and ``lambda=D^{-1/2}`` from ``D=1+J^2/r^4-C``.

    ``D`` must be strictly positive.  The square root is the positive branch
    so ``lambda>0``.  A zero or sign-changing ``D`` is fail-closed.
    """

    radius = _as_interval(radius)
    compactness = _as_interval(compactness)
    angular_momentum = _as_interval(angular_momentum)
    if not radius.strictly_positive():
        raise SGBLInitialHealthStop(
            "domain_error",
            "Misner-Sharp chart reconstruction requires positive radius",
        )
    angular_extrinsic_curvature = angular_momentum / (radius ** 3)
    denominator = (
        Interval.singleton(1)
        + (angular_momentum ** 2) / (radius ** 4)
        - compactness
    )
    if not denominator.strictly_positive():
        raise SGBLInitialHealthStop(
            "interval_inconclusive",
            "chart denominator D=1+J^2/r^4-C is not strictly positive",
            {
                "obstruction": "chart_denominator_not_strictly_positive",
                "denominator": (denominator.lower, denominator.upper),
            },
        )
    sqrt_denominator = sgbl_interval_sqrt(denominator)
    radial_metric = sqrt_denominator.reciprocal()
    if not radial_metric.strictly_positive():
        raise SGBLInitialHealthStop(
            "interval_inconclusive",
            "positive square-root branch did not produce lambda>0",
            {"obstruction": "sqrt_branch_not_positive"},
        )
    return {
        "angular_extrinsic_curvature": angular_extrinsic_curvature,
        "denominator": denominator,
        "sqrt_denominator": sqrt_denominator,
        "radial_metric": radial_metric,
        "negative_sqrt_branch": -sqrt_denominator,
    }


def sgbl_misner_sharp_chart_roundtrip(
    radius: Interval,
    radial_metric: Interval,
    angular_extrinsic_curvature: Interval,
) -> dict[str, Interval]:
    """Exact chart equivalence: ``(lambda,k) -> (C,J) -> (lambda,k)``.

    At singletons with a perfect-square denominator the residuals are
    ``{0}``.  In general the residual enclosures must contain zero.
    """

    compactness = sgbl_misner_sharp_compactness_interval(
        radius, radial_metric, angular_extrinsic_curvature
    )
    angular_momentum = sgbl_misner_sharp_angular_momentum_interval(
        radius, angular_extrinsic_curvature
    )
    reconstructed = sgbl_misner_sharp_chart_state(
        radius, compactness, angular_momentum
    )
    lambda_inv_sq = (radial_metric ** 2).reciprocal()
    return {
        "compactness": compactness,
        "angular_momentum": angular_momentum,
        "denominator": reconstructed["denominator"],
        "lambda_residual": reconstructed["radial_metric"] - radial_metric,
        "k_residual": reconstructed["angular_extrinsic_curvature"]
        - angular_extrinsic_curvature,
        "denominator_residual": reconstructed["denominator"] - lambda_inv_sq,
        "reconstructed_lambda": reconstructed["radial_metric"],
        "reconstructed_k": reconstructed["angular_extrinsic_curvature"],
    }


def sgbl_misner_sharp_ck_chart_state(
    radius: Interval,
    compactness: Interval,
    angular_extrinsic_curvature: Interval,
) -> dict[str, Interval]:
    """Reconstruct ``lambda=D^{-1/2}`` from ``D=1+r^2 k^2-C``.

    ``k`` is a chart coordinate, not reconstructed from ``J``.  ``D`` must
    be strictly positive.  The square root is the positive branch.
    """

    radius = _as_interval(radius)
    compactness = _as_interval(compactness)
    angular_extrinsic_curvature = _as_interval(angular_extrinsic_curvature)
    if not radius.strictly_positive():
        raise SGBLInitialHealthStop(
            "domain_error",
            "Misner-Sharp (C,k) chart reconstruction requires positive radius",
        )
    denominator = (
        Interval.singleton(1)
        + (radius ** 2) * (angular_extrinsic_curvature ** 2)
        - compactness
    )
    if not denominator.strictly_positive():
        raise SGBLInitialHealthStop(
            "interval_inconclusive",
            "chart denominator D=1+r^2 k^2-C is not strictly positive",
            {
                "obstruction": "chart_denominator_not_strictly_positive",
                "denominator": (denominator.lower, denominator.upper),
            },
        )
    sqrt_denominator = sgbl_interval_sqrt(denominator)
    radial_metric = sqrt_denominator.reciprocal()
    if not radial_metric.strictly_positive():
        raise SGBLInitialHealthStop(
            "interval_inconclusive",
            "positive square-root branch did not produce lambda>0",
            {"obstruction": "sqrt_branch_not_positive"},
        )
    return {
        "angular_extrinsic_curvature": angular_extrinsic_curvature,
        "denominator": denominator,
        "sqrt_denominator": sqrt_denominator,
        "radial_metric": radial_metric,
        "negative_sqrt_branch": -sqrt_denominator,
    }


def sgbl_misner_sharp_ck_chart_roundtrip(
    radius: Interval,
    radial_metric: Interval,
    angular_extrinsic_curvature: Interval,
) -> dict[str, Interval]:
    """Exact ``(C,k)`` chart equivalence: ``(lambda,k) -> (C,k) -> lambda``.

    At singletons with a perfect-square denominator the lambda residual is
    ``{0}``.  The ``k`` coordinate is copied, not reconstructed.
    """

    compactness = sgbl_misner_sharp_compactness_interval(
        radius, radial_metric, angular_extrinsic_curvature
    )
    reconstructed = sgbl_misner_sharp_ck_chart_state(
        radius, compactness, angular_extrinsic_curvature
    )
    lambda_inv_sq = (radial_metric ** 2).reciprocal()
    return {
        "compactness": compactness,
        "denominator": reconstructed["denominator"],
        "lambda_residual": reconstructed["radial_metric"] - radial_metric,
        "k_residual": reconstructed["angular_extrinsic_curvature"]
        - angular_extrinsic_curvature,
        "denominator_residual": reconstructed["denominator"] - lambda_inv_sq,
        "reconstructed_lambda": reconstructed["radial_metric"],
        "reconstructed_k": reconstructed["angular_extrinsic_curvature"],
    }


def _observe(budget: _ODEBudget, policy: SGBLValidatedODEPolicy, *values: Interval) -> None:
    for value in values:
        budget.max_bits = max(budget.max_bits, _interval_bits(value))
        if budget.max_bits > policy.max_rational_bit_length:
            raise SGBLInitialHealthStop(
                "resource_limit",
                "interval endpoints exceeded the declared rational bit cap",
                {"observed": budget.max_bits, "limit": policy.max_rational_bit_length},
            )


def _count_eval(budget: _ODEBudget, policy: SGBLValidatedODEPolicy) -> None:
    budget.rhs_evaluations += 1
    if budget.rhs_evaluations > policy.max_rhs_evaluations:
        raise SGBLInitialHealthStop(
            "resource_limit",
            "validated ODE exceeded the residual-evaluation budget",
            {
                "evaluations": budget.rhs_evaluations,
                "limit": policy.max_rhs_evaluations,
            },
        )


def _solve_affine_rhs(coefficients: Mapping[str, object]) -> tuple[object, object]:
    h0 = coefficients["hamiltonian_constant"]
    h_l = coefficients["hamiltonian_lambda_r"]
    m0 = coefficients["momentum_constant"]
    m_k = coefficients["momentum_k_r"]
    if _primal(h_l).contains_zero() or _primal(m_k).contains_zero():
        raise SGBLInitialHealthStop(
            "interval_inconclusive",
            "affine constraint Jacobian diagonal contains zero",
            {
                "obstruction": "jacobian_diagonal_contains_zero",
                "hamiltonian_lambda_r": (_primal(h_l).lower, _primal(h_l).upper),
                "momentum_k_r": (_primal(m_k).lower, _primal(m_k).upper),
            },
        )
    return -h0 / h_l, -m0 / m_k


def _affine_eval(
    radius: Interval,
    radial_metric: Interval | IntervalFirstTangent,
    angular_extrinsic_curvature: Interval | IntervalFirstTangent,
    spec: SGBLExactInitialSlice,
    policy: SGBLValidatedODEPolicy,
    budget: _ODEBudget,
) -> dict[str, Interval | IntervalFirstTangent]:
    _count_eval(budget, policy)
    try:
        return sgbl_generic_affine_coefficients(
            radius=radius,
            radial_metric=radial_metric,
            angular_extrinsic_curvature=angular_extrinsic_curvature,
            spec=spec,
        )
    except ZeroDivisionError as exc:
        raise SGBLInitialHealthStop(
            "interval_inconclusive",
            "interval reciprocal encountered a zero-containing denominator",
            {"obstruction": "zero_in_interval_reciprocal"},
        ) from exc


def sgbl_constraint_rhs_state_jacobian(
    radius: Interval,
    radial_metric: Interval,
    angular_extrinsic_curvature: Interval,
    spec: SGBLExactInitialSlice,
    *,
    policy: SGBLValidatedODEPolicy | None = None,
    budget: _ODEBudget | None = None,
) -> dict[str, Interval]:
    """Enclose ``d(lambda_r,k_r)/d(lambda,k)`` by interval first tangents."""

    policy = policy or DECLARED_ODE_POLICY
    budget = budget or _ODEBudget()
    lambda_seed = IntervalFirstTangent.seed(radial_metric)
    k_held = IntervalFirstTangent.constant(angular_extrinsic_curvature)
    k_seed = IntervalFirstTangent.seed(angular_extrinsic_curvature)
    lambda_held = IntervalFirstTangent.constant(radial_metric)
    lambda_column = _solve_affine_rhs(
        _affine_eval(radius, lambda_seed, k_held, spec, policy, budget)
    )
    k_column = _solve_affine_rhs(
        _affine_eval(radius, lambda_held, k_seed, spec, policy, budget)
    )
    compactness_lambda = sgbl_misner_sharp_compactness_radial_derivative(
        radius, lambda_seed, k_held, lambda_column[0], lambda_column[1]
    )
    compactness_k = sgbl_misner_sharp_compactness_radial_derivative(
        radius, lambda_held, k_seed, k_column[0], k_column[1]
    )
    if not isinstance(compactness_lambda, IntervalFirstTangent) or not isinstance(
        compactness_k, IntervalFirstTangent
    ):
        raise SGBLInitialHealthStop(
            "interval_inconclusive",
            "compactness-derivative jet failed to carry a first tangent",
            {"obstruction": "mean_value_inconsistent"},
        )
    return {
        "d_lambda_r_d_lambda": _outward(primal_and_tangent(lambda_column[0])[1]),
        "d_k_r_d_lambda": _outward(primal_and_tangent(lambda_column[1])[1]),
        "d_lambda_r_d_k": _outward(primal_and_tangent(k_column[0])[1]),
        "d_k_r_d_k": _outward(primal_and_tangent(k_column[1])[1]),
        "d_compactness_r_d_lambda": _outward(compactness_lambda.tangent),
        "d_compactness_r_d_k": _outward(compactness_k.tangent),
    }


def _centered_state_rhs(
    radius: Interval,
    radial_metric: Interval,
    angular_extrinsic_curvature: Interval,
    spec: SGBLExactInitialSlice,
    policy: SGBLValidatedODEPolicy,
    budget: _ODEBudget,
) -> dict[str, Interval]:
    """Centered mean-value of ``(lambda_r,k_r,C_r)`` from one Jacobian evaluation."""

    lambda_mid = Interval.singleton(radial_metric.midpoint())
    k_mid = Interval.singleton(angular_extrinsic_curvature.midpoint())
    center_rhs = _solve_affine_rhs(
        _affine_eval(radius, lambda_mid, k_mid, spec, policy, budget)
    )
    center_lambda_r = _outward(_primal(center_rhs[0]))
    center_k_r = _outward(_primal(center_rhs[1]))
    jacobian = sgbl_constraint_rhs_state_jacobian(
        radius,
        radial_metric,
        angular_extrinsic_curvature,
        spec,
        policy=policy,
        budget=budget,
    )
    deviation_lambda = radial_metric - lambda_mid
    deviation_k = angular_extrinsic_curvature - k_mid
    lambda_r = _outward(
        center_lambda_r
        + jacobian["d_lambda_r_d_lambda"] * deviation_lambda
        + jacobian["d_lambda_r_d_k"] * deviation_k
    )
    k_r = _outward(
        center_k_r
        + jacobian["d_k_r_d_lambda"] * deviation_lambda
        + jacobian["d_k_r_d_k"] * deviation_k
    )
    center_compactness_r = _outward(
        _primal(
            sgbl_misner_sharp_compactness_radial_derivative(
                radius, lambda_mid, k_mid, center_lambda_r, center_k_r
            )
        )
    )
    compactness_r = _outward(
        center_compactness_r
        + jacobian["d_compactness_r_d_lambda"] * deviation_lambda
        + jacobian["d_compactness_r_d_k"] * deviation_k
    )
    return {
        "centered_lambda_r": lambda_r,
        "centered_k_r": k_r,
        "centered_compactness_r": compactness_r,
        "center_lambda_r": center_lambda_r,
        "center_k_r": center_k_r,
        "center_compactness_r": center_compactness_r,
        "d_lambda_r_d_lambda": jacobian["d_lambda_r_d_lambda"],
        "d_k_r_d_lambda": jacobian["d_k_r_d_lambda"],
        "d_lambda_r_d_k": jacobian["d_lambda_r_d_k"],
        "d_k_r_d_k": jacobian["d_k_r_d_k"],
        "d_compactness_r_d_lambda": jacobian["d_compactness_r_d_lambda"],
        "d_compactness_r_d_k": jacobian["d_compactness_r_d_k"],
    }


def sgbl_centered_mean_value_rhs(
    radius: Interval,
    radial_metric: Interval,
    angular_extrinsic_curvature: Interval,
    spec: SGBLExactInitialSlice,
    *,
    policy: SGBLValidatedODEPolicy | None = None,
    budget: _ODEBudget | None = None,
) -> tuple[Interval, Interval]:
    """Centered mean-value enclosure of ``F`` in ``(lambda,k)`` at fixed radius box."""

    policy = policy or DECLARED_ODE_POLICY
    budget = budget or _ODEBudget()
    payload = _centered_state_rhs(
        radius, radial_metric, angular_extrinsic_curvature, spec, policy, budget
    )
    return payload["centered_lambda_r"], payload["centered_k_r"]


def sgbl_centered_mean_value_compactness_rhs(
    radius: Interval,
    radial_metric: Interval,
    angular_extrinsic_curvature: Interval,
    spec: SGBLExactInitialSlice,
    *,
    policy: SGBLValidatedODEPolicy | None = None,
    budget: _ODEBudget | None = None,
) -> Interval:
    """Centered mean-value enclosure of ``C_r`` in ``(lambda,k)`` at fixed radius box."""

    policy = policy or DECLARED_ODE_POLICY
    budget = budget or _ODEBudget()
    payload = _centered_state_rhs(
        radius, radial_metric, angular_extrinsic_curvature, spec, policy, budget
    )
    return payload["centered_compactness_r"]


def _rhs(
    radius: Interval,
    radial_metric: Interval,
    angular_extrinsic_curvature: Interval,
    spec: SGBLExactInitialSlice,
    policy: SGBLValidatedODEPolicy,
    budget: _ODEBudget,
) -> tuple[Interval, Interval, dict[str, Interval]]:
    coefficients = _affine_eval(
        radius, radial_metric, angular_extrinsic_curvature, spec, policy, budget
    )
    h0 = _outward(_primal(coefficients["hamiltonian_constant"]))
    h_l = _outward(_primal(coefficients["hamiltonian_lambda_r"]))
    m0 = _outward(_primal(coefficients["momentum_constant"]))
    m_k = _outward(_primal(coefficients["momentum_k_r"]))
    _observe(budget, policy, h0, h_l, m0, m_k)
    natural_lambda_r, natural_k_r = _solve_affine_rhs(
        {
            "hamiltonian_constant": h0,
            "hamiltonian_lambda_r": h_l,
            "momentum_constant": m0,
            "momentum_k_r": m_k,
        }
    )
    natural_lambda_r = _outward(_primal(natural_lambda_r))
    natural_k_r = _outward(_primal(natural_k_r))
    centered = _centered_state_rhs(
        radius,
        radial_metric,
        angular_extrinsic_curvature,
        spec,
        policy,
        budget,
    )
    centered_lambda_r = centered["centered_lambda_r"]
    centered_k_r = centered["centered_k_r"]
    lambda_r = _intersect(natural_lambda_r, centered_lambda_r)
    k_r = _intersect(natural_k_r, centered_k_r)
    if lambda_r is None or k_r is None:
        raise SGBLInitialHealthStop(
            "interval_inconclusive",
            "natural and centered mean-value RHS enclosures are disjoint",
            {"obstruction": "mean_value_inconsistent"},
        )
    lambda_r = _outward(lambda_r)
    k_r = _outward(k_r)
    if (
        lambda_r.width() < natural_lambda_r.width()
        or k_r.width() < natural_k_r.width()
    ):
        budget.mean_value_strictly_tighter = True
    natural_compactness_r = _outward(
        _primal(
            sgbl_misner_sharp_compactness_radial_derivative(
                radius,
                radial_metric,
                angular_extrinsic_curvature,
                natural_lambda_r,
                natural_k_r,
            )
        )
    )
    composed_compactness_r = _outward(
        _primal(
            sgbl_misner_sharp_compactness_radial_derivative(
                radius, radial_metric, angular_extrinsic_curvature, lambda_r, k_r
            )
        )
    )
    centered_compactness_r = centered["centered_compactness_r"]
    compactness_r = _intersect(natural_compactness_r, centered_compactness_r)
    if compactness_r is None:
        raise SGBLInitialHealthStop(
            "interval_inconclusive",
            "natural and centered mean-value compactness RHS enclosures are disjoint",
            {"obstruction": "compactness_mean_value_inconsistent"},
        )
    compactness_r = _intersect(compactness_r, composed_compactness_r)
    if compactness_r is None:
        raise SGBLInitialHealthStop(
            "interval_inconclusive",
            "composed compactness RHS is disjoint from the mean-value intersection",
            {"obstruction": "compactness_mean_value_inconsistent"},
        )
    compactness_r = _outward(compactness_r)
    budget.compactness_rhs_evaluations += 1
    if compactness_r.width() < natural_compactness_r.width():
        budget.compactness_mean_value_strictly_tighter = True
    hamiltonian = _outward(h0 + h_l * natural_lambda_r)
    momentum = _outward(m0 + m_k * natural_k_r)
    _observe(
        budget,
        policy,
        lambda_r,
        k_r,
        compactness_r,
        hamiltonian,
        momentum,
        natural_compactness_r,
        centered_compactness_r,
    )
    return lambda_r, k_r, {
        "hamiltonian_constant": h0,
        "hamiltonian_lambda_r": h_l,
        "momentum_constant": m0,
        "momentum_k_r": m_k,
        "lambda_r": lambda_r,
        "k_r": k_r,
        "compactness_r": compactness_r,
        "natural_lambda_r": natural_lambda_r,
        "natural_k_r": natural_k_r,
        "centered_lambda_r": centered_lambda_r,
        "centered_k_r": centered_k_r,
        "natural_compactness_r": natural_compactness_r,
        "centered_compactness_r": centered_compactness_r,
        "composed_compactness_r": composed_compactness_r,
        "hamiltonian_residual": hamiltonian,
        "momentum_residual": momentum,
    }


def _chart_rhs(
    radius: Interval,
    compactness: Interval,
    angular_momentum: Interval,
    spec: SGBLExactInitialSlice,
    policy: SGBLValidatedODEPolicy,
    budget: _ODEBudget,
) -> tuple[Interval, Interval, dict[str, Interval]]:
    state = sgbl_misner_sharp_chart_state(radius, compactness, angular_momentum)
    radial_metric = _intersect(state["radial_metric"], policy.lambda_domain)
    angular_extrinsic_curvature = _intersect(
        state["angular_extrinsic_curvature"], policy.k_domain
    )
    if radial_metric is None or angular_extrinsic_curvature is None:
        raise SGBLInitialHealthStop(
            "interval_inconclusive",
            "reconstructed lambda/k escaped the declared physical domain",
            {"obstruction": "physical_domain_escape"},
        )
    radial_metric = _outward(radial_metric)
    angular_extrinsic_curvature = _outward(angular_extrinsic_curvature)
    lambda_r, k_r, payload = _rhs(
        radius, radial_metric, angular_extrinsic_curvature, spec, policy, budget
    )
    angular_momentum_r = _outward(
        sgbl_misner_sharp_angular_momentum_radial_derivative(
            radius, angular_momentum, k_r
        )
    )
    composed_angular_momentum_r = _outward(
        3 * (radius ** 2) * angular_extrinsic_curvature + (radius ** 3) * k_r
    )
    angular_momentum_rhs = _intersect(angular_momentum_r, composed_angular_momentum_r)
    if angular_momentum_rhs is None:
        raise SGBLInitialHealthStop(
            "interval_inconclusive",
            "equivalent J_r formulae are disjoint",
            {"obstruction": "chart_roundtrip_misses_origin"},
        )
    angular_momentum_r = _outward(angular_momentum_rhs)
    _observe(budget, policy, angular_momentum_r, state["denominator"], radial_metric)
    payload = dict(payload)
    payload.update(
        {
            "angular_momentum_r": angular_momentum_r,
            "denominator": _outward(state["denominator"]),
            "sqrt_denominator": _outward(state["sqrt_denominator"]),
            "reconstructed_lambda": radial_metric,
            "reconstructed_k": angular_extrinsic_curvature,
        }
    )
    return payload["compactness_r"], angular_momentum_r, payload


def _ck_chart_rhs(
    radius: Interval,
    compactness: Interval,
    angular_extrinsic_curvature: Interval,
    spec: SGBLExactInitialSlice,
    policy: SGBLValidatedODEPolicy,
    budget: _ODEBudget,
) -> tuple[Interval, Interval, dict[str, Interval]]:
    state = sgbl_misner_sharp_ck_chart_state(
        radius, compactness, angular_extrinsic_curvature
    )
    radial_metric = _intersect(state["radial_metric"], policy.lambda_domain)
    k_phys = _intersect(angular_extrinsic_curvature, policy.k_domain)
    if radial_metric is None or k_phys is None:
        raise SGBLInitialHealthStop(
            "interval_inconclusive",
            "reconstructed lambda or chart k escaped the declared physical domain",
            {"obstruction": "physical_domain_escape"},
        )
    radial_metric = _outward(radial_metric)
    k_phys = _outward(k_phys)
    lambda_r, k_r, payload = _rhs(radius, radial_metric, k_phys, spec, policy, budget)
    _observe(budget, policy, state["denominator"], radial_metric, k_phys)
    payload = dict(payload)
    payload.update(
        {
            "denominator": _outward(state["denominator"]),
            "sqrt_denominator": _outward(state["sqrt_denominator"]),
            "reconstructed_lambda": radial_metric,
            "reconstructed_k": k_phys,
        }
    )
    return payload["compactness_r"], k_r, payload


@dataclass(frozen=True, slots=True, kw_only=True)
class SGBLODECellEnclosure:
    """Validated graph enclosure of a constraint chart on one radial cell."""

    radius: Interval
    lambda_box: Interval
    k_box: Interval
    compactness_box: Interval
    j_box: Interval
    lambda_left: Interval
    k_left: Interval
    compactness_left: Interval
    j_left: Interval
    lambda_right: Interval
    k_right: Interval
    compactness_right: Interval
    j_right: Interval
    compactness: Interval
    direct_compactness: Interval
    centered_compactness: Interval
    propagated_compactness: Interval
    invariant_residual: Interval
    denominator: Interval
    denominator_margin: Fraction
    chart_coordinates: bool
    chart_kind: str
    depth: int
    picard_strict_self_map: bool
    residual_contains_origin: bool
    compactness_invariant_contains_zero: bool
    correlated_strictly_tighter: bool

    @property
    def compactness_upper(self) -> Fraction:
        return self.compactness.upper


def _try_picard(
    r0: Fraction,
    r1: Fraction,
    lambda_left: Interval,
    k_left: Interval,
    compactness_left: Interval,
    spec: SGBLExactInitialSlice,
    policy: SGBLValidatedODEPolicy,
    budget: _ODEBudget,
    depth: int,
) -> SGBLODECellEnclosure | str:
    radius = interval(r0, r1)
    width = r1 - r0
    if width <= 0:
        return "domain_error"
    step = interval(0, width)
    try:
        seed_lambda_r, seed_k_r, seed_payload = _rhs(
            radius, lambda_left, k_left, spec, policy, budget
        )
    except SGBLInitialHealthStop as exc:
        if exc.reason in {"resource_limit", "domain_error"}:
            raise
        return str(exc.payload.get("obstruction") or "zero_in_interval_reciprocal")
    seed_compactness_r = seed_payload["compactness_r"]
    lambda_radius = DECLARED_EULER_INFLATION * width * seed_lambda_r.abs_upper()
    k_radius = DECLARED_EULER_INFLATION * width * seed_k_r.abs_upper()
    compactness_radius = DECLARED_EULER_INFLATION * width * seed_compactness_r.abs_upper()
    seed_lambda = _outward(
        Interval(lambda_left.lower - lambda_radius, lambda_left.upper + lambda_radius)
    )
    seed_k = _outward(Interval(k_left.lower - k_radius, k_left.upper + k_radius))
    seed_compactness = _outward(
        Interval(
            compactness_left.lower - compactness_radius,
            compactness_left.upper + compactness_radius,
        )
    )
    graph_lambda = _intersect(seed_lambda, policy.lambda_domain)
    graph_k = _intersect(seed_k, policy.k_domain)
    graph_compactness = _intersect(seed_compactness, policy.compactness_domain)
    if graph_lambda is None or graph_k is None or graph_compactness is None:
        return "physical_domain_escape"
    graph_lambda = _outward(graph_lambda)
    graph_k = _outward(graph_k)
    graph_compactness = _outward(graph_compactness)
    image_lambda = seed_lambda
    image_k = seed_k
    image_compactness = seed_compactness
    payload: dict[str, Interval] | None = None
    for _ in range(policy.max_picard_iterations):
        try:
            lambda_r, k_r, payload = _rhs(
                radius, graph_lambda, graph_k, spec, policy, budget
            )
        except SGBLInitialHealthStop as exc:
            if exc.reason in {"resource_limit", "domain_error"}:
                raise
            return str(exc.payload.get("obstruction") or "zero_in_interval_reciprocal")
        compactness_r = payload["compactness_r"]
        image_lambda = _outward(lambda_left + step * lambda_r)
        image_k = _outward(k_left + step * k_r)
        image_compactness = _outward(compactness_left + step * compactness_r)
        if (
            image_lambda.strictly_inside(graph_lambda)
            and image_k.strictly_inside(graph_k)
            and image_compactness.strictly_inside(graph_compactness)
        ):
            right_lambda = _outward(lambda_left + Interval.singleton(width) * lambda_r)
            right_k = _outward(k_left + Interval.singleton(width) * k_r)
            right_compactness = _outward(
                compactness_left + Interval.singleton(width) * compactness_r
            )
            propagated = _outward(image_compactness)
            try:
                cell_enclosures = sgbl_intersect_algebraic_compactness_enclosures(
                    radius, graph_lambda, graph_k, propagated=propagated
                )
                right_enclosures = sgbl_intersect_algebraic_compactness_enclosures(
                    Interval.singleton(r1),
                    right_lambda,
                    right_k,
                    propagated=right_compactness,
                )
            except SGBLInitialHealthStop as exc:
                if exc.reason in {"resource_limit", "domain_error"}:
                    raise
                return str(exc.payload.get("obstruction") or "zero_in_interval_reciprocal")
            direct = _outward(cell_enclosures["natural"])
            centered = _outward(cell_enclosures["centered"])
            correlated = _outward(cell_enclosures["correlated"])
            right_compactness = _outward(right_enclosures["correlated"])
            invariant = _outward(
                sgbl_misner_sharp_compactness_invariant_residual(
                    radius, graph_lambda, graph_k, correlated
                )
            )
            if not invariant.contains_zero():
                return "compactness_invariant_misses_origin"
            _observe(
                budget,
                policy,
                correlated,
                direct,
                centered,
                propagated,
                invariant,
                right_lambda,
                right_k,
                right_compactness,
            )
            residual_ok = payload["hamiltonian_residual"].contains_zero() and payload[
                "momentum_residual"
            ].contains_zero()
            if not residual_ok:
                return "residual_enclosure_misses_origin"
            if centered.width() < direct.width():
                budget.algebraic_mean_value_strictly_tighter = True
            tighter = correlated.width() < direct.width()
            if tighter:
                budget.correlated_strictly_tighter = True
            j_left = _outward(Interval.singleton(r0) ** 3 * k_left)
            j_right = _outward(Interval.singleton(r1) ** 3 * right_k)
            j_box = _outward((radius ** 3) * graph_k)
            denominator = _outward((graph_lambda ** 2).reciprocal())
            return SGBLODECellEnclosure(
                radius=radius,
                lambda_box=graph_lambda,
                k_box=graph_k,
                compactness_box=graph_compactness,
                j_box=j_box,
                lambda_left=lambda_left,
                k_left=k_left,
                compactness_left=compactness_left,
                j_left=j_left,
                lambda_right=right_lambda,
                k_right=right_k,
                compactness_right=right_compactness,
                j_right=j_right,
                compactness=correlated,
                direct_compactness=direct,
                centered_compactness=centered,
                propagated_compactness=propagated,
                invariant_residual=invariant,
                denominator=denominator,
                denominator_margin=denominator.lower,
                chart_coordinates=False,
                chart_kind="lambda_k",
                depth=depth,
                picard_strict_self_map=True,
                residual_contains_origin=True,
                compactness_invariant_contains_zero=True,
                correlated_strictly_tighter=tighter,
            )
        next_lambda = _intersect(image_lambda, graph_lambda)
        next_k = _intersect(image_k, graph_k)
        next_compactness = _intersect(image_compactness, graph_compactness)
        if next_lambda is None or next_k is None or next_compactness is None:
            return "physical_domain_escape"
        next_lambda = _outward(next_lambda)
        next_k = _outward(next_k)
        next_compactness = _outward(next_compactness)
        if (
            next_lambda == graph_lambda
            and next_k == graph_k
            and next_compactness == graph_compactness
        ):
            break
        graph_lambda, graph_k, graph_compactness = next_lambda, next_k, next_compactness
    return "picard_strict_self_map_failed"


def sgbl_validate_ode_cell_inventory(
    cells: tuple[SGBLODECellEnclosure, ...],
    *,
    origin: Fraction,
    terminus: Fraction | None = None,
) -> tuple[Fraction, Fraction]:
    """Require a strictly ordered, nonoverlapping, gap-free radial inventory.

    Coverage runs from ``origin`` through the last returned right boundary.
    If ``terminus`` is supplied, the inventory must end there exactly.
    An empty inventory covers the degenerate interval ``[origin, origin]``.
    """

    if not cells:
        if terminus is not None and terminus != origin:
            raise ValueError("empty ODE inventory cannot cover a nonempty interval")
        return origin, origin
    if cells[0].radius.lower != origin:
        raise ValueError("ODE inventory does not start at the requested left boundary")
    previous_right = origin
    for cell in cells:
        if cell.radius.upper <= cell.radius.lower:
            raise ValueError("ODE cell radius is not a positive-length interval")
        if cell.radius.lower < previous_right:
            raise ValueError("ODE inventory contains overlapping cells")
        if cell.radius.lower > previous_right:
            raise ValueError("ODE inventory contains a radial gap")
        previous_right = cell.radius.upper
    if terminus is not None and previous_right != terminus:
        raise ValueError("ODE inventory does not reach the required right boundary")
    return origin, previous_right


def _combine_partition(
    origin: Fraction,
    terminus: Fraction,
    left_cells: tuple[SGBLODECellEnclosure, ...],
    right_cells: tuple[SGBLODECellEnclosure, ...],
    right_reason: str | None,
) -> tuple[tuple[SGBLODECellEnclosure, ...], str | None]:
    """Keep left siblings and append the right recursion in radial order."""

    combined = left_cells + right_cells
    if right_reason is None:
        sgbl_validate_ode_cell_inventory(combined, origin=origin, terminus=terminus)
    else:
        sgbl_validate_ode_cell_inventory(combined, origin=origin)
    return combined, right_reason


def _partition_cell(
    r0: Fraction,
    r1: Fraction,
    lambda_left: Interval,
    k_left: Interval,
    compactness_left: Interval,
    spec: SGBLExactInitialSlice,
    policy: SGBLValidatedODEPolicy,
    budget: _ODEBudget,
    depth: int,
) -> tuple[tuple[SGBLODECellEnclosure, ...], str | None]:
    budget.cells += 1
    if budget.cells > policy.max_cells:
        raise SGBLInitialHealthStop(
            "resource_limit",
            "validated ODE exceeded the declared cell cap",
            {"cells": budget.cells, "limit": policy.max_cells},
        )
    outcome = _try_picard(
        r0, r1, lambda_left, k_left, compactness_left, spec, policy, budget, depth
    )
    if type(outcome) is SGBLODECellEnclosure:
        if outcome.compactness.upper < 1:
            return (outcome,), None
        if depth >= policy.max_bisection_depth:
            return (outcome,), "compactness_enclosure_not_below_one"
        reason = "compactness_enclosure_not_below_one"
    else:
        reason = outcome
        if depth >= policy.max_bisection_depth:
            return (), reason
    mid = (r0 + r1) / 2
    left_cells, left_reason = _partition_cell(
        r0,
        mid,
        lambda_left,
        k_left,
        compactness_left,
        spec,
        policy,
        budget,
        depth + 1,
    )
    if left_reason is not None:
        sgbl_validate_ode_cell_inventory(left_cells, origin=r0)
        return left_cells, left_reason
    sgbl_validate_ode_cell_inventory(left_cells, origin=r0, terminus=mid)
    last = left_cells[-1]
    right_cells, right_reason = _partition_cell(
        mid,
        r1,
        last.lambda_right,
        last.k_right,
        last.compactness_right,
        spec,
        policy,
        budget,
        depth + 1,
    )
    return _combine_partition(r0, r1, left_cells, right_cells, right_reason)


def _try_picard_chart(
    r0: Fraction,
    r1: Fraction,
    compactness_left: Interval,
    angular_momentum_left: Interval,
    spec: SGBLExactInitialSlice,
    policy: SGBLValidatedODEPolicy,
    budget: _ODEBudget,
    depth: int,
) -> SGBLODECellEnclosure | str:
    radius = interval(r0, r1)
    width = r1 - r0
    if width <= 0:
        return "domain_error"
    step = interval(0, width)
    try:
        seed_compactness_r, seed_angular_momentum_r, _seed_payload = _chart_rhs(
            radius, compactness_left, angular_momentum_left, spec, policy, budget
        )
    except SGBLInitialHealthStop as exc:
        if exc.reason in {"resource_limit", "domain_error"}:
            raise
        return str(exc.payload.get("obstruction") or "zero_in_interval_reciprocal")
    compactness_radius = DECLARED_EULER_INFLATION * width * seed_compactness_r.abs_upper()
    angular_momentum_radius = (
        DECLARED_EULER_INFLATION * width * seed_angular_momentum_r.abs_upper()
    )
    seed_compactness = _outward(
        Interval(
            compactness_left.lower - compactness_radius,
            compactness_left.upper + compactness_radius,
        )
    )
    seed_angular_momentum = _outward(
        Interval(
            angular_momentum_left.lower - angular_momentum_radius,
            angular_momentum_left.upper + angular_momentum_radius,
        )
    )
    graph_compactness = _intersect(seed_compactness, policy.compactness_domain)
    graph_angular_momentum = _intersect(
        seed_angular_momentum, policy.angular_momentum_domain
    )
    if graph_compactness is None or graph_angular_momentum is None:
        return "physical_domain_escape"
    graph_compactness = _outward(graph_compactness)
    graph_angular_momentum = _outward(graph_angular_momentum)
    payload: dict[str, Interval] | None = None
    for _ in range(policy.max_picard_iterations):
        try:
            compactness_r, angular_momentum_r, payload = _chart_rhs(
                radius,
                graph_compactness,
                graph_angular_momentum,
                spec,
                policy,
                budget,
            )
        except SGBLInitialHealthStop as exc:
            if exc.reason in {"resource_limit", "domain_error"}:
                raise
            return str(exc.payload.get("obstruction") or "zero_in_interval_reciprocal")
        image_compactness = _outward(compactness_left + step * compactness_r)
        image_angular_momentum = _outward(
            angular_momentum_left + step * angular_momentum_r
        )
        if image_compactness.strictly_inside(
            graph_compactness
        ) and image_angular_momentum.strictly_inside(graph_angular_momentum):
            right_compactness = _outward(
                compactness_left + Interval.singleton(width) * compactness_r
            )
            right_angular_momentum = _outward(
                angular_momentum_left + Interval.singleton(width) * angular_momentum_r
            )
            propagated = _outward(image_compactness)
            reconstructed_lambda = payload["reconstructed_lambda"]
            reconstructed_k = payload["reconstructed_k"]
            denominator = payload["denominator"]
            try:
                cell_enclosures = sgbl_intersect_algebraic_compactness_enclosures(
                    radius,
                    reconstructed_lambda,
                    reconstructed_k,
                    propagated=propagated,
                )
                right_state = sgbl_misner_sharp_chart_state(
                    Interval.singleton(r1), right_compactness, right_angular_momentum
                )
                right_enclosures = sgbl_intersect_algebraic_compactness_enclosures(
                    Interval.singleton(r1),
                    right_state["radial_metric"],
                    right_state["angular_extrinsic_curvature"],
                    propagated=right_compactness,
                )
            except SGBLInitialHealthStop as exc:
                if exc.reason in {"resource_limit", "domain_error"}:
                    raise
                return str(exc.payload.get("obstruction") or "zero_in_interval_reciprocal")
            direct = _outward(cell_enclosures["natural"])
            centered = _outward(cell_enclosures["centered"])
            correlated = _outward(cell_enclosures["correlated"])
            right_compactness = _outward(right_enclosures["correlated"])
            right_lambda = _outward(right_state["radial_metric"])
            right_k = _outward(right_state["angular_extrinsic_curvature"])
            invariant = _outward(
                sgbl_misner_sharp_compactness_invariant_residual(
                    radius, reconstructed_lambda, reconstructed_k, correlated
                )
            )
            if not invariant.contains_zero():
                return "chart_roundtrip_misses_origin"
            _observe(
                budget,
                policy,
                correlated,
                direct,
                centered,
                propagated,
                invariant,
                reconstructed_lambda,
                reconstructed_k,
                right_compactness,
                right_angular_momentum,
            )
            residual_ok = payload["hamiltonian_residual"].contains_zero() and payload[
                "momentum_residual"
            ].contains_zero()
            if not residual_ok:
                return "residual_enclosure_misses_origin"
            if centered.width() < direct.width():
                budget.algebraic_mean_value_strictly_tighter = True
            tighter = correlated.width() < direct.width()
            if tighter:
                budget.correlated_strictly_tighter = True
            try:
                left_state = sgbl_misner_sharp_chart_state(
                    Interval.singleton(r0), compactness_left, angular_momentum_left
                )
            except SGBLInitialHealthStop as exc:
                if exc.reason in {"resource_limit", "domain_error"}:
                    raise
                return str(exc.payload.get("obstruction") or "zero_in_interval_reciprocal")
            return SGBLODECellEnclosure(
                radius=radius,
                lambda_box=reconstructed_lambda,
                k_box=reconstructed_k,
                compactness_box=graph_compactness,
                j_box=graph_angular_momentum,
                lambda_left=_outward(left_state["radial_metric"]),
                k_left=_outward(left_state["angular_extrinsic_curvature"]),
                compactness_left=compactness_left,
                j_left=angular_momentum_left,
                lambda_right=right_lambda,
                k_right=right_k,
                compactness_right=right_compactness,
                j_right=right_angular_momentum,
                compactness=correlated,
                direct_compactness=direct,
                centered_compactness=centered,
                propagated_compactness=propagated,
                invariant_residual=invariant,
                denominator=denominator,
                denominator_margin=denominator.lower,
                chart_coordinates=True,
                chart_kind="c_j",
                depth=depth,
                picard_strict_self_map=True,
                residual_contains_origin=True,
                compactness_invariant_contains_zero=True,
                correlated_strictly_tighter=tighter,
            )
        next_compactness = _intersect(image_compactness, graph_compactness)
        next_angular_momentum = _intersect(
            image_angular_momentum, graph_angular_momentum
        )
        if next_compactness is None or next_angular_momentum is None:
            return "physical_domain_escape"
        next_compactness = _outward(next_compactness)
        next_angular_momentum = _outward(next_angular_momentum)
        if (
            next_compactness == graph_compactness
            and next_angular_momentum == graph_angular_momentum
        ):
            break
        graph_compactness, graph_angular_momentum = (
            next_compactness,
            next_angular_momentum,
        )
    return "picard_strict_self_map_failed"


def _partition_chart_cell(
    r0: Fraction,
    r1: Fraction,
    compactness_left: Interval,
    angular_momentum_left: Interval,
    spec: SGBLExactInitialSlice,
    policy: SGBLValidatedODEPolicy,
    budget: _ODEBudget,
    depth: int,
) -> tuple[tuple[SGBLODECellEnclosure, ...], str | None]:
    budget.cells += 1
    if budget.cells > policy.max_cells:
        raise SGBLInitialHealthStop(
            "resource_limit",
            "validated ODE exceeded the declared cell cap",
            {"cells": budget.cells, "limit": policy.max_cells},
        )
    outcome = _try_picard_chart(
        r0,
        r1,
        compactness_left,
        angular_momentum_left,
        spec,
        policy,
        budget,
        depth,
    )
    if type(outcome) is SGBLODECellEnclosure:
        if outcome.compactness.upper < 1:
            return (outcome,), None
        if depth >= policy.max_bisection_depth:
            return (outcome,), "compactness_enclosure_not_below_one"
        reason = "compactness_enclosure_not_below_one"
    else:
        reason = outcome
        if depth >= policy.max_bisection_depth:
            return (), reason
    mid = (r0 + r1) / 2
    left_cells, left_reason = _partition_chart_cell(
        r0,
        mid,
        compactness_left,
        angular_momentum_left,
        spec,
        policy,
        budget,
        depth + 1,
    )
    if left_reason is not None:
        sgbl_validate_ode_cell_inventory(left_cells, origin=r0)
        return left_cells, left_reason
    sgbl_validate_ode_cell_inventory(left_cells, origin=r0, terminus=mid)
    last = left_cells[-1]
    right_cells, right_reason = _partition_chart_cell(
        mid,
        r1,
        last.compactness_right,
        last.j_right,
        spec,
        policy,
        budget,
        depth + 1,
    )
    return _combine_partition(r0, r1, left_cells, right_cells, right_reason)


def _try_picard_ck(
    r0: Fraction,
    r1: Fraction,
    compactness_left: Interval,
    k_left: Interval,
    spec: SGBLExactInitialSlice,
    policy: SGBLValidatedODEPolicy,
    budget: _ODEBudget,
    depth: int,
) -> SGBLODECellEnclosure | str:
    radius = interval(r0, r1)
    width = r1 - r0
    if width <= 0:
        return "domain_error"
    step = interval(0, width)
    try:
        seed_compactness_r, seed_k_r, _seed_payload = _ck_chart_rhs(
            radius, compactness_left, k_left, spec, policy, budget
        )
    except SGBLInitialHealthStop as exc:
        if exc.reason in {"resource_limit", "domain_error"}:
            raise
        return str(exc.payload.get("obstruction") or "zero_in_interval_reciprocal")
    compactness_radius = DECLARED_EULER_INFLATION * width * seed_compactness_r.abs_upper()
    k_radius = DECLARED_EULER_INFLATION * width * seed_k_r.abs_upper()
    seed_compactness = _outward(
        Interval(
            compactness_left.lower - compactness_radius,
            compactness_left.upper + compactness_radius,
        )
    )
    seed_k = _outward(Interval(k_left.lower - k_radius, k_left.upper + k_radius))
    graph_compactness = _intersect(seed_compactness, policy.compactness_domain)
    graph_k = _intersect(seed_k, policy.k_domain)
    if graph_compactness is None or graph_k is None:
        return "physical_domain_escape"
    graph_compactness = _outward(graph_compactness)
    graph_k = _outward(graph_k)
    payload: dict[str, Interval] | None = None
    for _ in range(policy.max_picard_iterations):
        try:
            compactness_r, k_r, payload = _ck_chart_rhs(
                radius, graph_compactness, graph_k, spec, policy, budget
            )
        except SGBLInitialHealthStop as exc:
            if exc.reason in {"resource_limit", "domain_error"}:
                raise
            return str(exc.payload.get("obstruction") or "zero_in_interval_reciprocal")
        image_compactness = _outward(compactness_left + step * compactness_r)
        image_k = _outward(k_left + step * k_r)
        if image_compactness.strictly_inside(graph_compactness) and image_k.strictly_inside(
            graph_k
        ):
            right_compactness = _outward(
                compactness_left + Interval.singleton(width) * compactness_r
            )
            right_k = _outward(k_left + Interval.singleton(width) * k_r)
            propagated = _outward(image_compactness)
            reconstructed_lambda = payload["reconstructed_lambda"]
            denominator = payload["denominator"]
            try:
                cell_enclosures = sgbl_intersect_algebraic_compactness_enclosures(
                    radius,
                    reconstructed_lambda,
                    graph_k,
                    propagated=propagated,
                )
                right_state = sgbl_misner_sharp_ck_chart_state(
                    Interval.singleton(r1), right_compactness, right_k
                )
                right_enclosures = sgbl_intersect_algebraic_compactness_enclosures(
                    Interval.singleton(r1),
                    right_state["radial_metric"],
                    right_k,
                    propagated=right_compactness,
                )
            except SGBLInitialHealthStop as exc:
                if exc.reason in {"resource_limit", "domain_error"}:
                    raise
                return str(exc.payload.get("obstruction") or "zero_in_interval_reciprocal")
            direct = _outward(cell_enclosures["natural"])
            centered = _outward(cell_enclosures["centered"])
            correlated = _outward(cell_enclosures["correlated"])
            right_compactness = _outward(right_enclosures["correlated"])
            right_lambda = _outward(right_state["radial_metric"])
            invariant = _outward(
                sgbl_misner_sharp_compactness_invariant_residual(
                    radius, reconstructed_lambda, graph_k, correlated
                )
            )
            if not invariant.contains_zero():
                return "chart_roundtrip_misses_origin"
            _observe(
                budget,
                policy,
                correlated,
                direct,
                centered,
                propagated,
                invariant,
                reconstructed_lambda,
                graph_k,
                right_compactness,
                right_k,
            )
            residual_ok = payload["hamiltonian_residual"].contains_zero() and payload[
                "momentum_residual"
            ].contains_zero()
            if not residual_ok:
                return "residual_enclosure_misses_origin"
            if centered.width() < direct.width():
                budget.algebraic_mean_value_strictly_tighter = True
            tighter = correlated.width() < direct.width()
            if tighter:
                budget.correlated_strictly_tighter = True
            try:
                left_state = sgbl_misner_sharp_ck_chart_state(
                    Interval.singleton(r0), compactness_left, k_left
                )
            except SGBLInitialHealthStop as exc:
                if exc.reason in {"resource_limit", "domain_error"}:
                    raise
                return str(exc.payload.get("obstruction") or "zero_in_interval_reciprocal")
            j_left = _outward(Interval.singleton(r0) ** 3 * k_left)
            j_right = _outward(Interval.singleton(r1) ** 3 * right_k)
            j_box = _outward((radius ** 3) * graph_k)
            return SGBLODECellEnclosure(
                radius=radius,
                lambda_box=reconstructed_lambda,
                k_box=graph_k,
                compactness_box=graph_compactness,
                j_box=j_box,
                lambda_left=_outward(left_state["radial_metric"]),
                k_left=k_left,
                compactness_left=compactness_left,
                j_left=j_left,
                lambda_right=right_lambda,
                k_right=right_k,
                compactness_right=right_compactness,
                j_right=j_right,
                compactness=correlated,
                direct_compactness=direct,
                centered_compactness=centered,
                propagated_compactness=propagated,
                invariant_residual=invariant,
                denominator=denominator,
                denominator_margin=denominator.lower,
                chart_coordinates=True,
                chart_kind=CERTIFICATE_CHART_KIND,
                depth=depth,
                picard_strict_self_map=True,
                residual_contains_origin=True,
                compactness_invariant_contains_zero=True,
                correlated_strictly_tighter=tighter,
            )
        next_compactness = _intersect(image_compactness, graph_compactness)
        next_k = _intersect(image_k, graph_k)
        if next_compactness is None or next_k is None:
            return "physical_domain_escape"
        next_compactness = _outward(next_compactness)
        next_k = _outward(next_k)
        if next_compactness == graph_compactness and next_k == graph_k:
            break
        graph_compactness, graph_k = next_compactness, next_k
    return "picard_strict_self_map_failed"


def _partition_ck_cell(
    r0: Fraction,
    r1: Fraction,
    compactness_left: Interval,
    k_left: Interval,
    spec: SGBLExactInitialSlice,
    policy: SGBLValidatedODEPolicy,
    budget: _ODEBudget,
    depth: int,
) -> tuple[tuple[SGBLODECellEnclosure, ...], str | None]:
    budget.cells += 1
    if budget.cells > policy.max_cells:
        raise SGBLInitialHealthStop(
            "resource_limit",
            "validated ODE exceeded the declared cell cap",
            {"cells": budget.cells, "limit": policy.max_cells},
        )
    outcome = _try_picard_ck(
        r0, r1, compactness_left, k_left, spec, policy, budget, depth
    )
    if type(outcome) is SGBLODECellEnclosure:
        if outcome.compactness.upper < 1:
            return (outcome,), None
        if depth >= policy.max_bisection_depth:
            return (outcome,), "compactness_enclosure_not_below_one"
        reason = "compactness_enclosure_not_below_one"
    else:
        reason = outcome
        if depth >= policy.max_bisection_depth:
            return (), reason
    mid = (r0 + r1) / 2
    left_cells, left_reason = _partition_ck_cell(
        r0, mid, compactness_left, k_left, spec, policy, budget, depth + 1
    )
    if left_reason is not None:
        sgbl_validate_ode_cell_inventory(left_cells, origin=r0)
        return left_cells, left_reason
    sgbl_validate_ode_cell_inventory(left_cells, origin=r0, terminus=mid)
    last = left_cells[-1]
    right_cells, right_reason = _partition_ck_cell(
        mid,
        r1,
        last.compactness_right,
        last.k_right,
        spec,
        policy,
        budget,
        depth + 1,
    )
    return _combine_partition(r0, r1, left_cells, right_cells, right_reason)


@dataclass(frozen=True, slots=True, kw_only=True)
class SGBLValidatedODERecord:
    """Piecewise validated enclosure of the compact-support constraint ODE."""

    spec: SGBLExactInitialSlice
    policy: SGBLValidatedODEPolicy
    cells: tuple[SGBLODECellEnclosure, ...]
    rhs_evaluations: int
    compactness_rhs_evaluations: int
    max_rational_bit_length_observed: int
    mean_value_strictly_tighter: bool
    compactness_mean_value_strictly_tighter: bool
    algebraic_mean_value_strictly_tighter: bool
    correlated_strictly_tighter_than_direct: bool
    chart_coordinates: bool
    chart_kind: str
    obstruction: str | None

    def __post_init__(self) -> None:
        _, coverage_right = sgbl_validate_ode_cell_inventory(
            self.cells, origin=self.spec.support_minimum
        )
        if self.obstruction is None:
            if not self.cells:
                raise ValueError("a complete ODE record requires a nonempty cell inventory")
            if coverage_right != self.spec.support_maximum:
                raise ValueError("a complete ODE record must tile the entire compact support")

    @property
    def coverage_left(self) -> Fraction:
        left, _right = sgbl_validate_ode_cell_inventory(
            self.cells, origin=self.spec.support_minimum
        )
        return left

    @property
    def coverage_right(self) -> Fraction:
        _left, right = sgbl_validate_ode_cell_inventory(
            self.cells, origin=self.spec.support_minimum
        )
        return right

    @property
    def tiles_compact_support(self) -> bool:
        return (
            bool(self.cells)
            and self.coverage_left == self.spec.support_minimum
            and self.coverage_right == self.spec.support_maximum
        )

    @property
    def covers_complete_support(self) -> bool:
        return self.obstruction is None and self.tiles_compact_support

    @property
    def lambda_end(self) -> Interval:
        return self.cells[-1].lambda_right if self.cells else DECLARED_LAMBDA_DOMAIN

    @property
    def k_end(self) -> Interval:
        return self.cells[-1].k_right if self.cells else DECLARED_K_DOMAIN

    @property
    def compactness_end(self) -> Interval:
        return self.cells[-1].compactness_right if self.cells else DECLARED_C_DOMAIN

    @property
    def angular_momentum_end(self) -> Interval:
        return self.cells[-1].j_right if self.cells else DECLARED_J_DOMAIN

    @property
    def denominator_margin(self) -> Fraction | None:
        if not self.cells:
            return None
        return min(cell.denominator_margin for cell in self.cells)

    @property
    def support_compactness(self) -> Interval:
        if not self.cells:
            return Interval.singleton(2)
        lower = min(cell.compactness.lower for cell in self.cells)
        upper = max(cell.compactness.upper for cell in self.cells)
        return Interval(lower, upper)

    @property
    def direct_support_compactness(self) -> Interval:
        if not self.cells:
            return Interval.singleton(2)
        lower = min(cell.direct_compactness.lower for cell in self.cells)
        upper = max(cell.direct_compactness.upper for cell in self.cells)
        return Interval(lower, upper)

    @property
    def invariant_residual(self) -> Interval | None:
        if not self.cells:
            return None
        lower = min(cell.invariant_residual.lower for cell in self.cells)
        upper = max(cell.invariant_residual.upper for cell in self.cells)
        return Interval(lower, upper)


def sgbl_validated_lambda_k_constraint_ode(
    spec: SGBLExactInitialSlice,
    *,
    policy: SGBLValidatedODEPolicy | None = None,
) -> SGBLValidatedODERecord:
    """Regression enclosure of ``(lambda,k,C)`` on the compact support."""

    if type(spec) is not SGBLExactInitialSlice:
        raise TypeError("spec must be SGBLExactInitialSlice")
    policy = policy or DECLARED_ODE_POLICY
    if type(policy) is not SGBLValidatedODEPolicy:
        raise TypeError("policy must be SGBLValidatedODEPolicy")
    start = spec.support_minimum
    end = spec.support_maximum
    width = (end - start) / policy.base_cells
    budget = _ODEBudget()
    cells: list[SGBLODECellEnclosure] = []
    lambda_left = Interval.singleton(1)
    k_left = Interval.singleton(0)
    compactness_left = Interval.singleton(0)
    obstruction: str | None = None
    for index in range(policy.base_cells):
        left = start + index * width
        right = end if index == policy.base_cells - 1 else start + (index + 1) * width
        piece, obstruction = _partition_cell(
            left,
            right,
            lambda_left,
            k_left,
            compactness_left,
            spec,
            policy,
            budget,
            0,
        )
        cells.extend(piece)
        if obstruction is not None:
            break
        lambda_left = cells[-1].lambda_right
        k_left = cells[-1].k_right
        compactness_left = cells[-1].compactness_right
    return SGBLValidatedODERecord(
        spec=spec,
        policy=policy,
        cells=tuple(cells),
        rhs_evaluations=budget.rhs_evaluations,
        compactness_rhs_evaluations=budget.compactness_rhs_evaluations,
        max_rational_bit_length_observed=budget.max_bits,
        mean_value_strictly_tighter=budget.mean_value_strictly_tighter,
        compactness_mean_value_strictly_tighter=budget.compactness_mean_value_strictly_tighter,
        algebraic_mean_value_strictly_tighter=budget.algebraic_mean_value_strictly_tighter,
        correlated_strictly_tighter_than_direct=budget.correlated_strictly_tighter,
        chart_coordinates=False,
        chart_kind="lambda_k",
        obstruction=obstruction,
    )


def sgbl_validated_cj_constraint_ode(
    spec: SGBLExactInitialSlice,
    *,
    policy: SGBLValidatedODEPolicy | None = None,
) -> SGBLValidatedODERecord:
    """Regression enclosure of the Misner-Sharp ``(C,J)`` chart."""

    if type(spec) is not SGBLExactInitialSlice:
        raise TypeError("spec must be SGBLExactInitialSlice")
    policy = policy or DECLARED_ODE_POLICY
    if type(policy) is not SGBLValidatedODEPolicy:
        raise TypeError("policy must be SGBLValidatedODEPolicy")
    start = spec.support_minimum
    end = spec.support_maximum
    width = (end - start) / policy.base_cells
    budget = _ODEBudget()
    cells: list[SGBLODECellEnclosure] = []
    compactness_left = Interval.singleton(0)
    angular_momentum_left = Interval.singleton(0)
    obstruction: str | None = None
    for index in range(policy.base_cells):
        left = start + index * width
        right = end if index == policy.base_cells - 1 else start + (index + 1) * width
        piece, obstruction = _partition_chart_cell(
            left,
            right,
            compactness_left,
            angular_momentum_left,
            spec,
            policy,
            budget,
            0,
        )
        cells.extend(piece)
        if obstruction is not None:
            break
        compactness_left = cells[-1].compactness_right
        angular_momentum_left = cells[-1].j_right
    return SGBLValidatedODERecord(
        spec=spec,
        policy=policy,
        cells=tuple(cells),
        rhs_evaluations=budget.rhs_evaluations,
        compactness_rhs_evaluations=budget.compactness_rhs_evaluations,
        max_rational_bit_length_observed=budget.max_bits,
        mean_value_strictly_tighter=budget.mean_value_strictly_tighter,
        compactness_mean_value_strictly_tighter=budget.compactness_mean_value_strictly_tighter,
        algebraic_mean_value_strictly_tighter=budget.algebraic_mean_value_strictly_tighter,
        correlated_strictly_tighter_than_direct=budget.correlated_strictly_tighter,
        chart_coordinates=True,
        chart_kind="c_j",
        obstruction=obstruction,
    )


def sgbl_validated_constraint_ode(
    spec: SGBLExactInitialSlice,
    *,
    policy: SGBLValidatedODEPolicy | None = None,
) -> SGBLValidatedODERecord:
    """Enclose compact-support geometry in the Misner-Sharp ``(C,k)`` chart."""

    if type(spec) is not SGBLExactInitialSlice:
        raise TypeError("spec must be SGBLExactInitialSlice")
    policy = policy or DECLARED_ODE_POLICY
    if type(policy) is not SGBLValidatedODEPolicy:
        raise TypeError("policy must be SGBLValidatedODEPolicy")
    start = spec.support_minimum
    end = spec.support_maximum
    width = (end - start) / policy.base_cells
    budget = _ODEBudget()
    cells: list[SGBLODECellEnclosure] = []
    compactness_left = Interval.singleton(0)
    k_left = Interval.singleton(0)
    obstruction: str | None = None
    for index in range(policy.base_cells):
        left = start + index * width
        right = end if index == policy.base_cells - 1 else start + (index + 1) * width
        piece, obstruction = _partition_ck_cell(
            left,
            right,
            compactness_left,
            k_left,
            spec,
            policy,
            budget,
            0,
        )
        cells.extend(piece)
        if obstruction is not None:
            break
        compactness_left = cells[-1].compactness_right
        k_left = cells[-1].k_right
    return SGBLValidatedODERecord(
        spec=spec,
        policy=policy,
        cells=tuple(cells),
        rhs_evaluations=budget.rhs_evaluations,
        compactness_rhs_evaluations=budget.compactness_rhs_evaluations,
        max_rational_bit_length_observed=budget.max_bits,
        mean_value_strictly_tighter=budget.mean_value_strictly_tighter,
        compactness_mean_value_strictly_tighter=budget.compactness_mean_value_strictly_tighter,
        algebraic_mean_value_strictly_tighter=budget.algebraic_mean_value_strictly_tighter,
        correlated_strictly_tighter_than_direct=budget.correlated_strictly_tighter,
        chart_coordinates=True,
        chart_kind=CERTIFICATE_CHART_KIND,
        obstruction=obstruction,
    )


def sgbl_exact_minkowski_rhs_reduction(spec: SGBLExactInitialSlice) -> dict[str, Fraction]:
    """Exact affine RHS at the buffer endpoint must match the interval extension."""

    point = SGBLAnnularInputs(
        radius=spec.support_minimum,
        radial_metric=1,
        angular_extrinsic_curvature=0,
        phi=0,
        phi_r=0,
        phi_rr=0,
        phi_pi=0,
        phi_pi_r=0,
        chi_r=0,
        chi_pi=0,
        planck_mass=spec.planck_mass,
        scalar_mass=spec.scalar_mass,
        quartic_coupling=spec.quartic_coupling,
        alpha_gb=spec.alpha_gb,
    )
    exact = sgbl_constraint_coefficients(point)
    interval_coefficients = sgbl_interval_affine_coefficients(
        radius=Interval.singleton(spec.support_minimum),
        radial_metric=Interval.singleton(1),
        angular_extrinsic_curvature=Interval.singleton(0),
        spec=spec,
    )
    for name in (
        "hamiltonian_constant",
        "hamiltonian_lambda_r",
        "momentum_constant",
        "momentum_k_r",
    ):
        boxed = interval_coefficients[name]
        value = getattr(exact, name)
        if not boxed.is_singleton() or boxed.lower != value:
            raise ValueError(f"{name} interval reduction differs from the exact kernel")
    lambda_r = -exact.hamiltonian_constant / exact.hamiltonian_lambda_r
    k_r = -exact.momentum_constant / exact.momentum_k_r
    radius = Interval.singleton(spec.support_minimum)
    radial_metric = Interval.singleton(1)
    angular_extrinsic_curvature = Interval.singleton(0)
    compactness = sgbl_misner_sharp_compactness_interval(
        radius, radial_metric, angular_extrinsic_curvature
    )
    compactness_r = sgbl_misner_sharp_compactness_radial_derivative(
        radius,
        radial_metric,
        angular_extrinsic_curvature,
        Interval.singleton(lambda_r),
        Interval.singleton(k_r),
    )
    identity = sgbl_misner_sharp_compactness_chain_rule_identity(
        radius,
        radial_metric,
        angular_extrinsic_curvature,
        Interval.singleton(lambda_r),
        Interval.singleton(k_r),
    )
    angular_momentum = sgbl_misner_sharp_angular_momentum_interval(
        radius, angular_extrinsic_curvature
    )
    angular_momentum_r = sgbl_misner_sharp_angular_momentum_radial_derivative(
        radius, angular_momentum, Interval.singleton(k_r)
    )
    roundtrip = sgbl_misner_sharp_chart_roundtrip(
        radius, radial_metric, angular_extrinsic_curvature
    )
    ck_roundtrip = sgbl_misner_sharp_ck_chart_roundtrip(
        radius, radial_metric, angular_extrinsic_curvature
    )
    if (
        not compactness.is_singleton()
        or compactness.lower != 0
        or not _primal(compactness_r).is_singleton()
        or _primal(compactness_r).lower != 0
        or not identity["chain_rule_residual"].is_singleton()
        or identity["chain_rule_residual"].lower != 0
        or not identity["invariant_residual"].is_singleton()
        or identity["invariant_residual"].lower != 0
        or not angular_momentum.is_singleton()
        or angular_momentum.lower != 0
        or not angular_momentum_r.is_singleton()
        or angular_momentum_r.lower != 0
        or not roundtrip["lambda_residual"].is_singleton()
        or roundtrip["lambda_residual"].lower != 0
        or not roundtrip["k_residual"].is_singleton()
        or roundtrip["k_residual"].lower != 0
        or not roundtrip["denominator"].is_singleton()
        or roundtrip["denominator"].lower != 1
        or not ck_roundtrip["lambda_residual"].is_singleton()
        or ck_roundtrip["lambda_residual"].lower != 0
        or not ck_roundtrip["k_residual"].is_singleton()
        or ck_roundtrip["k_residual"].lower != 0
        or not ck_roundtrip["denominator"].is_singleton()
        or ck_roundtrip["denominator"].lower != 1
    ):
        raise ValueError("Minkowski-buffer compactness jet must vanish")
    return {
        "lambda_r": lambda_r,
        "k_r": k_r,
        "compactness": compactness.lower,
        "compactness_r": _primal(compactness_r).lower,
        "angular_momentum": angular_momentum.lower,
        "angular_momentum_r": angular_momentum_r.lower,
        "denominator": roundtrip["denominator"].lower,
        "hamiltonian_lambda_r": exact.hamiltonian_lambda_r,
        "momentum_k_r": exact.momentum_k_r,
    }


@dataclass(frozen=True, slots=True, kw_only=True)
class SGBLContinuousCompactnessRecord:
    """Continuum compactness certificate.  Sampled nodes are not stored."""

    slice: SGBLExactInitialSlice
    policy: SGBLValidatedODEPolicy
    minkowski_buffer_compactness: Fraction
    compact_support_compactness_upper: Fraction
    exterior_compactness_upper: Fraction
    compactness_margin: Fraction
    ode_cell_count: int
    picard_strict_inclusions: int
    rhs_evaluations: int
    compactness_rhs_evaluations: int
    classification: str
    inconclusive_reason: str | None
    missing_theorem: str | None
    theorem: str
    direct_compact_support_compactness_upper: Fraction
    correlated_strictly_tighter_than_direct: bool
    compactness_invariant_residual: Interval | None = None
    denominator_margin: Fraction | None = None
    chart_coordinates: bool = True
    chart_kind: str = CERTIFICATE_CHART_KIND
    ode: SGBLValidatedODERecord | None = None

    def __post_init__(self) -> None:
        if type(self.slice) is not SGBLExactInitialSlice:
            raise TypeError("slice must be SGBLExactInitialSlice")
        if type(self.policy) is not SGBLValidatedODEPolicy:
            raise TypeError("policy must be SGBLValidatedODEPolicy")
        if self.minkowski_buffer_compactness != 0:
            raise ValueError("exact Minkowski buffer compactness must be zero")
        if self.classification not in {
            "continuous_no_initial_trapped_sphere",
            "interval_inconclusive",
        }:
            raise ValueError("unknown compactness classification")
        if self.ode is not None and self.ode_cell_count != len(self.ode.cells):
            raise ValueError("reported cell count does not match the ODE inventory")
        if self.classification == "continuous_no_initial_trapped_sphere":
            if self.ode is None or not self.ode.covers_complete_support:
                raise ValueError("proved continuum no-trap requires a complete support tiling")
        if self.continuous_no_initial_trapped_sphere:
            if self.classification != "continuous_no_initial_trapped_sphere":
                raise ValueError("proved continuum no-trap requires its classification")
            if not (
                self.compact_support_compactness_upper < 1
                and self.exterior_compactness_upper < 1
                and self.compactness_margin > 0
            ):
                raise ValueError("proved continuum no-trap requires C<1 throughout")
            if self.inconclusive_reason is not None or self.missing_theorem is not None:
                raise ValueError("proved continuum no-trap cannot retain an obstruction")
            if (
                self.compactness_invariant_residual is not None
                and not self.compactness_invariant_residual.contains_zero()
            ):
                raise ValueError(
                    "proved continuum no-trap requires the compactness invariant residual to contain zero"
                )
            if (
                self.ode is not None
                and not all(cell.compactness_invariant_contains_zero for cell in self.ode.cells)
            ):
                raise ValueError("proved continuum no-trap requires a per-cell compactness invariant")
            if not self.chart_coordinates or self.chart_kind != CERTIFICATE_CHART_KIND:
                raise ValueError("proved continuum no-trap requires the (C,k) chart Picard")
            if self.ode is not None and (
                not self.ode.chart_coordinates or self.ode.chart_kind != CERTIFICATE_CHART_KIND
            ):
                raise ValueError("proved continuum no-trap requires the (C,k) chart Picard")
            if self.ode is not None and not all(
                cell.chart_coordinates and cell.chart_kind == CERTIFICATE_CHART_KIND
                for cell in self.ode.cells
            ):
                raise ValueError("proved continuum no-trap requires per-cell (C,k) chart coordinates")
        elif self.inconclusive_reason is None:
            raise ValueError("failed continuum no-trap must be typed inconclusive")

    @property
    def compactness_upper_bound(self) -> Fraction:
        return max(
            self.minkowski_buffer_compactness,
            self.compact_support_compactness_upper,
            self.exterior_compactness_upper,
        )

    @property
    def coverage_left(self) -> Fraction | None:
        return None if self.ode is None else self.ode.coverage_left

    @property
    def coverage_right(self) -> Fraction | None:
        return None if self.ode is None else self.ode.coverage_right

    @property
    def covers_complete_support(self) -> bool:
        return self.ode is not None and self.ode.covers_complete_support

    @property
    def continuous_no_initial_trapped_sphere(self) -> bool:
        return (
            self.classification == "continuous_no_initial_trapped_sphere"
            and self.compactness_upper_bound < 1
            and self.compactness_margin > 0
            and self.covers_complete_support
        )

    @property
    def sampled_nodes_are_not_the_certificate(self) -> bool:
        return True

    @property
    def SGBL_branch_owned_and_healthy(self) -> bool:
        return False

    @property
    def copied_gr0_or_fgcqr_health_evidence(self) -> bool:
        return False


def _obstruction_theorem(reason: str | None) -> str | None:
    return {
        None: None,
        "picard_strict_self_map_failed": (
            "strict_interval_Picard_self_map_of_the_affine_constraint_ODE_on_the_declared_refinement_tree"
        ),
        "jacobian_diagonal_contains_zero": (
            "affine_H_L_and_M_k_diagonals_excluding_zero_on_every_declared_cell_box"
        ),
        "physical_domain_escape": (
            "Picard_image_remaining_inside_the_declared_physical_state_domain"
        ),
        "compactness_enclosure_not_below_one": (
            "correlated_interval_compactness_C_strictly_below_one_on_the_validated_augmented_graph"
        ),
        "residual_enclosure_misses_origin": (
            "affine_H_and_M_residual_enclosures_containing_the_origin_on_every_cell"
        ),
        "exterior_compactness_not_below_one": (
            "analytic_exterior_C_equals_two_M_over_r_with_maximum_strictly_below_one"
        ),
        "zero_in_interval_reciprocal": (
            "interval_reciprocals_of_r_and_lambda_remaining_defined_on_every_cell_box"
        ),
        "mean_value_inconsistent": (
            "centered_mean_value_RHS_intersecting_the_natural_interval_extension"
        ),
        "compactness_mean_value_inconsistent": (
            "centered_mean_value_compactness_RHS_intersecting_the_natural_chain_rule_extension"
        ),
        "algebraic_compactness_mean_value_inconsistent": (
            "centered_mean_value_algebraic_C_intersecting_the_natural_interval_extension"
        ),
        "compactness_invariant_misses_origin": (
            "algebraic_and_propagated_compactness_enclosures_overlapping_with_invariant_residual_zero"
        ),
        "chart_denominator_not_strictly_positive": (
            "chart_denominator_D_equals_one_plus_J_squared_over_r_to_the_fourth_minus_C_strictly_positive"
        ),
        "sqrt_branch_not_positive": (
            "positive_square_root_branch_of_D_producing_strictly_positive_lambda"
        ),
        "chart_roundtrip_misses_origin": (
            "Misner_Sharp_chart_roundtrip_residuals_containing_the_origin"
        ),
    }.get(reason, reason)


def _ode_with_obstruction(ode: SGBLValidatedODERecord, reason: str) -> SGBLValidatedODERecord:
    return SGBLValidatedODERecord(
        spec=ode.spec,
        policy=ode.policy,
        cells=ode.cells,
        rhs_evaluations=ode.rhs_evaluations,
        compactness_rhs_evaluations=ode.compactness_rhs_evaluations,
        max_rational_bit_length_observed=ode.max_rational_bit_length_observed,
        mean_value_strictly_tighter=ode.mean_value_strictly_tighter,
        compactness_mean_value_strictly_tighter=ode.compactness_mean_value_strictly_tighter,
        algebraic_mean_value_strictly_tighter=ode.algebraic_mean_value_strictly_tighter,
        correlated_strictly_tighter_than_direct=ode.correlated_strictly_tighter_than_direct,
        chart_coordinates=ode.chart_coordinates,
        chart_kind=ode.chart_kind,
        obstruction=reason,
    )


def sgbl_continuous_initial_compactness(
    spec: SGBLExactInitialSlice | SGBLFamilyParameters,
    *,
    policy: SGBLValidatedODEPolicy | None = None,
) -> SGBLContinuousCompactnessRecord:
    """Prove or refuse ``C=2m/R<1`` by a validated constraint-ODE enclosure."""

    if type(spec) is SGBLFamilyParameters:
        slice_spec = SGBLExactInitialSlice.from_family_parameters(spec)
    elif type(spec) is SGBLExactInitialSlice:
        slice_spec = spec
    else:
        raise TypeError("spec must be SGBLExactInitialSlice or SGBLFamilyParameters")
    policy = policy or DECLARED_ODE_POLICY
    minkowski = sgbl_exact_minkowski_rhs_reduction(slice_spec)
    if (
        minkowski["lambda_r"] != 0
        or minkowski["k_r"] != 0
        or minkowski["compactness"] != 0
        or minkowski["compactness_r"] != 0
        or minkowski["angular_momentum"] != 0
        or minkowski["angular_momentum_r"] != 0
        or minkowski["denominator"] != 1
    ):
        raise ValueError("Minkowski-buffer constraint RHS and compactness jet must vanish")
    ode = sgbl_validated_constraint_ode(slice_spec, policy=policy)
    support_upper = ode.support_compactness.upper if ode.cells else Q(2)
    direct_upper = ode.direct_support_compactness.upper if ode.cells else Q(2)
    exterior_upper = support_upper
    if ode.cells and ode.obstruction is None:
        end_radius = Interval.singleton(slice_spec.support_maximum)
        try:
            end_enclosures = sgbl_intersect_algebraic_compactness_enclosures(
                end_radius, ode.lambda_end, ode.k_end, propagated=ode.compactness_end
            )
        except SGBLInitialHealthStop as exc:
            if exc.reason in {"resource_limit", "domain_error"}:
                raise
            ode = _ode_with_obstruction(
                ode, str(exc.payload.get("obstruction") or "compactness_invariant_misses_origin")
            )
            end_direct = sgbl_misner_sharp_compactness_interval(
                end_radius, ode.lambda_end, ode.k_end
            )
            end_correlated = end_direct
        else:
            end_direct = end_enclosures["natural"]
            end_correlated = end_enclosures["correlated"]
        # Exterior C=2M/r with M=r_end C_end/2, so max_{r>=r_end} C = C_end when C_end>=0.
        exterior_upper = max(end_correlated.upper, Q(0))
        support_upper = max(support_upper, end_correlated.upper)
        direct_upper = max(direct_upper, end_direct.upper)
        if ode.obstruction is None and exterior_upper >= 1:
            ode = _ode_with_obstruction(ode, "exterior_compactness_not_below_one")
    proved = (
        ode.obstruction is None
        and ode.covers_complete_support
        and support_upper < 1
        and exterior_upper < 1
        and all(
            cell.picard_strict_self_map
            and cell.residual_contains_origin
            and cell.compactness_invariant_contains_zero
            for cell in ode.cells
        )
    )
    margin = 1 - max(support_upper, exterior_upper)
    tighter = ode.correlated_strictly_tighter_than_direct or (
        ode.cells and support_upper < direct_upper
    )
    theorem_core = (
        "On the exact Minkowski buffer (lambda,k,C)=(1,0,0) the constraint "
        "RHS vanishes and D=lambda^{-2}=1.  Compact-support geometry is enclosed "
        "in the Misner-Sharp chart (C,k) with C=1+r^2 k^2-lambda^{-2} and "
        "D=1+r^2 k^2-C=lambda^{-2}.  Reconstruction uses a strictly positive D "
        "and the positive square-root branch lambda=D^{-1/2}.  The affine H/M "
        "solve for lambda_r and k_r is unchanged.  The state ODEs are k_r and "
        "C_r=2 r k^2 + 2 r^2 k k_r + 2 lambda_r/lambda^3, started from "
        "(C,k)=(0,0).  Interval Picard self-maps in (C,k) on the declared tree.  "
        "The partition concatenates left and right children so the cell "
        "inventory is a strictly ordered, gap-free covering from the support "
        "minimum through the last proved right boundary; a complete record "
        "tiles the entire compact support.  "
        "The product-box (lambda,k,C) graph and the (C,J) chart are regression, "
        "not the certificate.  The analytic exterior satisfies C=2M/r with M "
        "taken from the enclosed support-boundary state.  Sampled RK4/SSPRK3 "
        "nodes are not the certificate."
    )
    if proved:
        return SGBLContinuousCompactnessRecord(
            slice=slice_spec,
            policy=policy,
            minkowski_buffer_compactness=Q(0),
            compact_support_compactness_upper=support_upper,
            exterior_compactness_upper=exterior_upper,
            compactness_margin=margin,
            ode_cell_count=len(ode.cells),
            picard_strict_inclusions=len(ode.cells),
            rhs_evaluations=ode.rhs_evaluations,
            compactness_rhs_evaluations=ode.compactness_rhs_evaluations,
            classification="continuous_no_initial_trapped_sphere",
            inconclusive_reason=None,
            missing_theorem=None,
            theorem=theorem_core,
            direct_compact_support_compactness_upper=direct_upper,
            correlated_strictly_tighter_than_direct=tighter,
            compactness_invariant_residual=ode.invariant_residual,
            denominator_margin=ode.denominator_margin,
            chart_coordinates=ode.chart_coordinates,
            chart_kind=ode.chart_kind,
            ode=ode,
        )
    reason = ode.obstruction or "compactness_enclosure_not_below_one"
    return SGBLContinuousCompactnessRecord(
        slice=slice_spec,
        policy=policy,
        minkowski_buffer_compactness=Q(0),
        compact_support_compactness_upper=support_upper,
        exterior_compactness_upper=exterior_upper,
        compactness_margin=margin,
        ode_cell_count=len(ode.cells),
        picard_strict_inclusions=sum(1 for cell in ode.cells if cell.picard_strict_self_map),
        rhs_evaluations=ode.rhs_evaluations,
        compactness_rhs_evaluations=ode.compactness_rhs_evaluations,
        classification="interval_inconclusive",
        inconclusive_reason=reason,
        missing_theorem=_obstruction_theorem(reason),
        theorem=(
            "The Minkowski buffer still has C=0.  The declared interval-Picard "
            "tree for the Misner-Sharp (C,k) chart ODE did not produce a complete "
            "validated graph with C<1.  Sampled integrator nodes cannot complete "
            f"the continuum proof.  Obstruction: {reason}."
        ),
        direct_compact_support_compactness_upper=direct_upper,
        correlated_strictly_tighter_than_direct=tighter,
        compactness_invariant_residual=ode.invariant_residual,
        denominator_margin=ode.denominator_margin,
        chart_coordinates=ode.chart_coordinates,
        chart_kind=ode.chart_kind,
        ode=ode,
    )


def sgbl_integrator_containment_regression(
    record: SGBLContinuousCompactnessRecord,
    *,
    method: str,
    step_count: int = 16,
) -> dict[str, Any]:
    """Regression-only: floating RK4/SSPRK3 nodes must lie in the interval graph.

    This is not a continuum proof.  The family integrators are sampled
    constraint ODE steps; the certificate is the Picard enclosure.
    """

    if method not in CONSTRAINT_INTEGRATOR_NAMES:
        raise ValueError("unknown SGB-L radial constraint integrator")
    if record.ode is None or not record.ode.cells:
        raise SGBLInitialHealthStop(
            "interval_inconclusive",
            "containment regression requires a nonempty validated graph",
        )
    parameters = SGBLFamilyParameters(chi_amplitude=float(record.slice.chi_amplitude))
    solution = integrate_sgbl_radial_constraints(
        parameters, step_count=step_count, method=method
    )
    contained = 0
    eligible = 0
    for point in solution.points:
        radius = Fraction(point.radius)
        if float(radius) != point.radius:
            raise SGBLInitialHealthStop(
                "domain_error",
                "regression node radius is not an exact dyadic",
            )
        cell = next(
            (
                item
                for item in record.ode.cells
                if item.radius.lower <= radius <= item.radius.upper
            ),
            None,
        )
        if cell is None:
            continue
        eligible += 1
        lambda_ok = float(cell.lambda_box.lower) <= point.radial_metric <= float(
            cell.lambda_box.upper
        )
        k_ok = float(cell.k_box.lower) <= point.angular_extrinsic_curvature <= float(
            cell.k_box.upper
        )
        if lambda_ok and k_ok:
            contained += 1
    return {
        "method": method,
        "step_count": step_count,
        "nodes": len(solution.points),
        "eligible_nodes": eligible,
        "contained_nodes": contained,
        "all_eligible_nodes_contained": eligible > 0 and contained == eligible,
        "covers_complete_support": record.continuous_no_initial_trapped_sphere,
        "not_the_continuum_certificate": True,
    }


def _chart_overlap_stats(
    chart_cells: tuple[SGBLODECellEnclosure, ...],
    other_cells: tuple[SGBLODECellEnclosure, ...],
) -> dict[str, int | bool]:
    overlapping = 0
    lambda_overlap = 0
    k_overlap = 0
    compactness_overlap = 0
    for chart_cell in chart_cells:
        for other_cell in other_cells:
            if _intersect(chart_cell.radius, other_cell.radius) is None:
                continue
            overlapping += 1
            if _intersect(chart_cell.lambda_box, other_cell.lambda_box) is not None:
                lambda_overlap += 1
            if _intersect(chart_cell.k_box, other_cell.k_box) is not None:
                k_overlap += 1
            if _intersect(chart_cell.compactness, other_cell.compactness) is not None:
                compactness_overlap += 1
            break
    return {
        "overlapping_radial_cells": overlapping,
        "lambda_enclosures_overlap": lambda_overlap,
        "k_enclosures_overlap": k_overlap,
        "compactness_enclosures_overlap": compactness_overlap,
        "all_overlapping_cells_consistent": (
            overlapping > 0
            and lambda_overlap == overlapping
            and k_overlap == overlapping
            and compactness_overlap == overlapping
        ),
    }


def sgbl_chart_versus_lambda_k_regression(
    record: SGBLContinuousCompactnessRecord,
) -> dict[str, Any]:
    """Regression-only comparison of the ``(C,k)`` chart against prior graphs."""

    if record.ode is None or not record.ode.cells:
        raise SGBLInitialHealthStop(
            "interval_inconclusive",
            "chart versus prior-chart regression requires a nonempty chart graph",
        )
    lambda_k = sgbl_validated_lambda_k_constraint_ode(record.slice, policy=record.policy)
    cj = sgbl_validated_cj_constraint_ode(record.slice, policy=record.policy)
    lambda_k_stats = _chart_overlap_stats(record.ode.cells, lambda_k.cells)
    cj_stats = _chart_overlap_stats(record.ode.cells, cj.cells)
    return {
        "not_the_continuum_certificate": True,
        "chart_kind": record.ode.chart_kind,
        "chart_cells": len(record.ode.cells),
        "chart_compactness_upper": record.ode.support_compactness.upper,
        "lambda_k_cells": len(lambda_k.cells),
        "lambda_k_obstruction": lambda_k.obstruction,
        "lambda_k_compactness_upper": lambda_k.support_compactness.upper,
        "lambda_k_overlapping_radial_cells": lambda_k_stats["overlapping_radial_cells"],
        "lambda_k_consistent": lambda_k_stats["all_overlapping_cells_consistent"],
        "cj_cells": len(cj.cells),
        "cj_obstruction": cj.obstruction,
        "cj_compactness_upper": cj.support_compactness.upper,
        "cj_overlapping_radial_cells": cj_stats["overlapping_radial_cells"],
        "cj_consistent": cj_stats["all_overlapping_cells_consistent"],
        "overlapping_radial_cells": lambda_k_stats["overlapping_radial_cells"],
        "all_overlapping_cells_consistent": lambda_k_stats["all_overlapping_cells_consistent"],
    }


@dataclass(frozen=True, slots=True, kw_only=True)
class SGBLEvolvingCenterMonitorContract:
    """Fail-closed evolving-centre monitor contract for later trajectories."""

    executed: bool
    classification: str
    grounded_in: tuple[str, ...]
    required_checks: tuple[str, ...]
    initial_empty_center_series: Mapping[str, Any]
    inconclusive_reason: str | None
    missing_theorem: str | None

    def __post_init__(self) -> None:
        if self.executed:
            raise ValueError("this slice has no trajectory owner that can execute the monitor")
        if self.classification != "monitor_not_executed":
            raise ValueError("unexecuted monitor classification differs")
        if self.pre_holdout_pass:
            raise ValueError("an unexecuted monitor cannot be a pre-holdout pass")
        object.__setattr__(
            self, "initial_empty_center_series", MappingProxyType(dict(self.initial_empty_center_series))
        )

    @property
    def pre_holdout_pass(self) -> bool:
        return False

    @property
    def unexecuted_monitor_is_not_a_pre_holdout_pass(self) -> bool:
        return True

    @property
    def evolving_center_after_matter_arrives(self) -> bool:
        return False

    @property
    def SGBL_branch_owned_and_healthy(self) -> bool:
        return False


def sgbl_evolving_center_monitor_contract() -> SGBLEvolvingCenterMonitorContract:
    """Return the unexecuted monitor contract grounded in the series owner."""

    profile = sgbl_empty_minkowski_center_profile()
    parsed = validate_sgbl_initial_center_profile(profile)
    series = sgbl_initial_center_series(profile)
    if not parsed["elementary_flatness_value_dt_dtt_valid"] or not series["initial_empty_center"]:
        raise SGBLInitialHealthStop(
            "domain_error",
            "evolving-centre contract requires the exact empty Minkowski series owner",
        )
    return SGBLEvolvingCenterMonitorContract(
        executed=False,
        classification="monitor_not_executed",
        grounded_in=(
            "sgb1_ctl1_center.sgbl_initial_center_series",
            "sgb1_ctl1_center.validate_sgbl_initial_center_profile",
            "regular_center.LaurentSeries",
            "regular_center.SeriesJet2",
        ),
        required_checks=MONITOR_REQUIRED_CHECKS,
        initial_empty_center_series={
            "profile_id": series["profile_id"],
            "center_limits": series["center_limits"],
            "initial_empty_center": series["initial_empty_center"],
            "evolving_center_after_matter_arrives": series["evolving_center_after_matter_arrives"],
            "certified_negative_powers_absent": series["certified_negative_powers_absent"],
        },
        inconclusive_reason="monitor_not_executed",
        missing_theorem="evolving_center_monitor_on_a_later_branch_owned_trajectory",
    )


def sgbl_execute_evolving_center_monitor(payload: object = None) -> SGBLEvolvingCenterMonitorContract:
    """Refuse execution: this slice has no trajectory owner."""

    raise SGBLInitialHealthStop(
        "missing_trajectory_owner" if payload is not None else "monitor_not_executed",
        "the evolving-centre monitor cannot run without a later trajectory owner",
        {"payload_supplied": payload is not None},
    )


def sgbl_initial_health_gate(
    compactness: SGBLContinuousCompactnessRecord | None = None,
    monitor: SGBLEvolvingCenterMonitorContract | None = None,
) -> dict[str, Any]:
    """Aggregate health remains closed.  Local continuum facts stay local."""

    return {
        "SGBL_branch_owned_and_healthy": False,
        "continuous_no_initial_trapped_sphere": (
            False if compactness is None else compactness.continuous_no_initial_trapped_sphere
        ),
        "sampled_nodes_are_not_the_certificate": True,
        "evolving_center_after_matter_arrives": False,
        "unexecuted_monitor_is_not_a_pre_holdout_pass": True,
        "monitor_executed": False if monitor is None else monitor.executed,
        "pre_holdout_pass": False,
        "FRZ1": False,
        "PREF1": False,
        "execution_authorized": False,
        "holdout_authorized": False,
        "copied_gr0_or_fgcqr_health_evidence": False,
        "declared_family_floats": {
            "chi_amplitude": DECLARED_CHI_AMPLITUDE,
            "phi_amplitude": DECLARED_PHI_AMPLITUDE,
            "center": DECLARED_CENTER,
            "half_width": DECLARED_HALF_WIDTH,
            "outer_radius": DECLARED_OUTER_RADIUS,
            "planck_mass": DECLARED_PLANCK_MASS,
            "scalar_mass": DECLARED_SCALAR_MASS,
            "quartic_coupling": DECLARED_QUARTIC_COUPLING,
            "alpha_gb": DECLARED_ALPHA_GB,
        },
        "missing_health_closes_the_gate": True,
    }


__all__ = [
    "BUMP_X_DERIVATIVE_BOUND",
    "BUMP_XX_DERIVATIVE_BOUND",
    "DECLARED_BASE_CELLS",
    "DECLARED_C_DOMAIN",
    "CERTIFICATE_CHART_KIND",
    "DECLARED_J_DOMAIN",
    "DECLARED_ODE_POLICY",
    "SGBLContinuousCompactnessRecord",
    "SGBLEvolvingCenterMonitorContract",
    "SGBLExactInitialSlice",
    "SGBLInitialHealthStop",
    "SGBLODECellEnclosure",
    "SGBLValidatedODEPolicy",
    "SGBLValidatedODERecord",
    "sgbl_chart_versus_lambda_k_regression",
    "sgbl_compact_bump_enclosure",
    "sgbl_continuous_initial_compactness",
    "sgbl_evolving_center_monitor_contract",
    "sgbl_exact_minkowski_rhs_reduction",
    "sgbl_execute_evolving_center_monitor",
    "sgbl_initial_health_gate",
    "sgbl_integrator_containment_regression",
    "sgbl_centered_mean_value_algebraic_compactness",
    "sgbl_centered_mean_value_compactness_rhs",
    "sgbl_centered_mean_value_rhs",
    "sgbl_constraint_rhs_state_jacobian",
    "sgbl_generic_affine_coefficients",
    "sgbl_interval_affine_coefficients",
    "sgbl_interval_family_fields",
    "sgbl_interval_sqrt",
    "sgbl_intersect_algebraic_compactness_enclosures",
    "sgbl_misner_sharp_angular_momentum_interval",
    "sgbl_misner_sharp_angular_momentum_radial_derivative",
    "sgbl_misner_sharp_chart_roundtrip",
    "sgbl_misner_sharp_chart_state",
    "sgbl_misner_sharp_ck_chart_roundtrip",
    "sgbl_misner_sharp_ck_chart_state",
    "sgbl_misner_sharp_compactness_chain_rule_identity",
    "sgbl_misner_sharp_compactness_interval",
    "sgbl_misner_sharp_compactness_invariant_residual",
    "sgbl_misner_sharp_compactness_radial_derivative",
    "sgbl_misner_sharp_compactness_state_gradient",
    "sgbl_validate_ode_cell_inventory",
    "sgbl_validated_cj_constraint_ode",
    "sgbl_validated_constraint_ode",
    "sgbl_validated_lambda_k_constraint_ode",
]
