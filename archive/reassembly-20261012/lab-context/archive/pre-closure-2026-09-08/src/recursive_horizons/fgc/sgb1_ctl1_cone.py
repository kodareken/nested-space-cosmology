"""Branch-owned SGB-L continuum cone, kinetic, and symmetrizer enclosure.

The coefficient map is the ACT1 principal symbol specialized to the linear
branch ``F=Mpl^2`` (constant), ``F'=0``, ``f'=alpha_gb``, and
``Hess(f)=alpha_gb Hess(phi)``.  Complete ``A``, ``B_i``, and ``C_ij``
tensors are enclosed over a caller-declared orthonormal-frame box by exact
interval arithmetic, not by an FGC-QR or HYP2 pass bit.

All-direction control uses the n-independent operator-norm sum of those
coefficient tensors together with an exact unnormalized radial Einstein--two-
scalar eigenframe: six four-dimensional rational kernels, a rank-24 frame,
and a Bauer--Fike resolvent bound from spectral gaps and a rational
``||V||_2||V^{-1}||_2`` upper bound.  Physical action-energy restrictions
are proved by exact Sylvester minors and a Gershgorin lower bound.  Arc
samples and binary64 spectra are regressions, never the proof.  Structural
Hessian symmetry and the pure-gauge identity are algebraic statements on a
complete exact basis: all ten symmetric Hessian generators and twenty
linearly independent algebraic-curvature generators (Kulkarni--Nomizu
products with exact rank-20 selection).  A single Hess(00) or R(0101)
representative is not sufficient.  Aggregate ``SGBL_branch_owned_and_healthy``
stays false.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from itertools import product
from math import isqrt
from numbers import Integral
from types import MappingProxyType
from typing import Any, Mapping, Sequence

import numpy as np

from .evolution.covariant_principal_health import (
    DIMENSION,
    FIELD_COUNT,
    FIELD_ORDER,
    FIRST_ORDER_FIELD_COUNT,
    METRIC_COMPONENTS,
    CovariantPrincipalBackground,
    esf_background,
    principal_coefficient_tensors,
)
from .exact_interval import Interval, interval, interval_matrix_infinity_row_sum_bound
from .exact_interval_linear_algebra import center_preconditioned_neumann_inverse
from .exact_linear_algebra import determinant, inverse, matrix_multiply, rank, transpose
from .sgb1_ctl1_principal import (
    ANGULAR_DIRECTION,
    DECLARED_HAT_NORMAL_FACTOR,
    DECLARED_TILDE_NORMAL_FACTOR,
    MIXED_DIRECTION,
    RADIAL_DIRECTION,
    sgbl_covariant_principal_background,
)
from .sgb1_ctl1_source import (
    HAT_NORMAL_FACTOR,
    SGBLSourceInputs,
    TILDE_NORMAL_FACTOR,
    sgbl_source_solve,
    sgbl_source_state,
)
from .spherical_reduction import SphericalState


Q = Fraction
ORTHONORMAL_MINKOWSKI_CHART = "orthonormal_minkowski_frame"
ZERO = Interval.singleton(0)
ONE = Interval.singleton(1)
HALF = Interval.singleton(Q(1, 2))
TWO = Interval.singleton(2)
FOUR = Interval.singleton(4)
ETA_SIGNS = (-1, 1, 1, 1)
BINARY64_UNIT = Q(1, 1 << 52)
DEFAULT_ARC_COUNT = 192
DEFAULT_SQRT2_BITS = 80
DEFAULT_ROUND_UNITS = 65536
DEFAULT_MAX_BITS = 4096
DEFAULT_MAX_EVALUATIONS = 128
DEFAULT_MAX_HALF_WIDTH = Q(1)
PHYSICAL_CONTOUR_RADIUS = Q(1, 8)
AUXILIARY_CONTOUR_RADIUS = Q(1, 16)
ESF_CLUSTER_CENTERS = (
    ("physical_minus", Q(-1)),
    ("tilde_minus", Q(-1, 2)),
    ("hat_minus", Q(-1, 3)),
    ("hat_plus", Q(1, 3)),
    ("tilde_plus", Q(1, 2)),
    ("physical_plus", Q(1)),
)
REGRESSION_DIRECTIONS = (RADIAL_DIRECTION, ANGULAR_DIRECTION, MIXED_DIRECTION)
FORBIDDEN_HEALTH_IMPORTS = (
    "WeakCouplingThresholds",
    "weak_coupling_health_certificate",
    "esf_reference_cluster_certificate",
    "canonical_health_monitor_values",
    "background_from_spherical_state",
)
SYMMETRIC_INDEX_PAIRS = tuple(
    (first, second)
    for first in range(DIMENSION)
    for second in range(first, DIMENSION)
)
HESSIAN_GENERATOR_COUNT = len(SYMMETRIC_INDEX_PAIRS)
ALGEBRAIC_CURVATURE_GENERATOR_COUNT = 20
KULKARNI_NOMIZU_CANDIDATE_COUNT = (
    HESSIAN_GENERATOR_COUNT * (HESSIAN_GENERATOR_COUNT + 1) // 2
)


def _fraction(name: str, value: object, *, nonnegative: bool = False) -> Fraction:
    if type(value) not in (int, Fraction):
        raise TypeError(f"{name} must be a Fraction or built-in int")
    result = Fraction(value)
    if nonnegative and result < 0:
        raise ValueError(f"{name} must be nonnegative")
    return result


def _positive_int(name: str, value: object) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or int(value) <= 0:
        raise ValueError(f"{name} must be a positive integer")
    return int(value)


def _interval(name: str, value: object) -> Interval:
    if type(value) is Interval:
        return value
    return Interval.singleton(_fraction(name, value))


def _interval_bits(value: Interval) -> int:
    return max(
        abs(value.lower.numerator).bit_length(),
        value.lower.denominator.bit_length(),
        abs(value.upper.numerator).bit_length(),
        value.upper.denominator.bit_length(),
    )


def _rational_sqrt_enclosure(value: Fraction, bits: int) -> Interval:
    """Return a closed rational interval containing ``sqrt(value)``."""

    if value < 0:
        raise ValueError("square root of a negative rational is not real")
    if value == 0:
        return Interval.singleton(0)
    scale = 1 << bits
    target = value * scale * scale
    floor_int = target.numerator // target.denominator
    root = isqrt(floor_int)
    while (root + 1) ** 2 <= target:
        root += 1
    while root > 0 and root ** 2 > target:
        root -= 1
    lower = Q(root, scale)
    if Q(root * root, scale * scale) == value:
        return Interval.singleton(lower)
    return Interval(lower, Q(root + 1, scale))


def _sqrt2_enclosure(bits: int) -> Interval:
    return _rational_sqrt_enclosure(Q(2), bits)


def _spectral_norm_upper(matrix: Sequence[Sequence[Interval]], bits: int) -> Fraction:
    """Return a rational upper bound of the matrix 2-norm over the interval."""

    inf_norm = interval_matrix_infinity_row_sum_bound(matrix)
    one_norm = max(
        sum(matrix[row][column].abs_upper() for row in range(len(matrix)))
        for column in range(len(matrix[0]))
    )
    return _rational_sqrt_enclosure(inf_norm * one_norm, bits).upper


def _matrix_sha256(matrices: Sequence[Sequence[Sequence[Interval]]]) -> str:
    parts: list[str] = []
    for matrix in matrices:
        for row in matrix:
            for entry in row:
                parts.append(
                    f"{entry.lower.numerator}/{entry.lower.denominator}:"
                    f"{entry.upper.numerator}/{entry.upper.denominator}"
                )
    return sha256(",".join(parts).encode()).hexdigest()


class SGBLConeStop(ArithmeticError):
    """Typed cone stop; wrapping/inconclusive records are not this type."""

    def __init__(
        self,
        reason: str,
        message: str,
        payload: Mapping[str, Any] | None = None,
    ) -> None:
        self.reason = reason
        self.payload = dict(payload or {})
        super().__init__(message)


@dataclass(frozen=True, slots=True, kw_only=True)
class SGBLConeLimits:
    """Prospective arithmetic, arc, and bit caps. Not fitted physical widths."""

    max_residual_evaluations: int = DEFAULT_MAX_EVALUATIONS
    max_rational_bit_length: int = DEFAULT_MAX_BITS
    max_parameter_half_width: Fraction | int = DEFAULT_MAX_HALF_WIDTH
    arc_count: int = DEFAULT_ARC_COUNT
    sqrt2_bits: int = DEFAULT_SQRT2_BITS
    roundoff_unit_count: int = DEFAULT_ROUND_UNITS
    physical_contour_radius: Fraction | int = PHYSICAL_CONTOUR_RADIUS
    auxiliary_contour_radius: Fraction | int = AUXILIARY_CONTOUR_RADIUS

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "max_residual_evaluations",
            _positive_int("max_residual_evaluations", self.max_residual_evaluations),
        )
        object.__setattr__(
            self,
            "max_rational_bit_length",
            _positive_int("max_rational_bit_length", self.max_rational_bit_length),
        )
        object.__setattr__(
            self,
            "arc_count",
            _positive_int("arc_count", self.arc_count),
        )
        object.__setattr__(
            self,
            "sqrt2_bits",
            _positive_int("sqrt2_bits", self.sqrt2_bits),
        )
        object.__setattr__(
            self,
            "roundoff_unit_count",
            _positive_int("roundoff_unit_count", self.roundoff_unit_count),
        )
        object.__setattr__(
            self,
            "max_parameter_half_width",
            _fraction(
                "max_parameter_half_width",
                self.max_parameter_half_width,
                nonnegative=True,
            ),
        )
        object.__setattr__(
            self,
            "physical_contour_radius",
            _fraction("physical_contour_radius", self.physical_contour_radius),
        )
        object.__setattr__(
            self,
            "auxiliary_contour_radius",
            _fraction("auxiliary_contour_radius", self.auxiliary_contour_radius),
        )
        if self.physical_contour_radius <= 0 or self.auxiliary_contour_radius <= 0:
            raise ValueError("contour radii must be positive")
        if self.arc_count == 512:
            raise SGBLConeStop(
                "resource_limit",
                "arc_count 512 is the HYP2 sample count and is not this policy",
            )
        _assert_contours_disjoint(
            self.physical_contour_radius,
            self.auxiliary_contour_radius,
        )


def _assert_contours_disjoint(physical: Fraction, auxiliary: Fraction) -> None:
    radii = {
        "physical_minus": physical,
        "physical_plus": physical,
        "tilde_minus": auxiliary,
        "tilde_plus": auxiliary,
        "hat_minus": auxiliary,
        "hat_plus": auxiliary,
    }
    centers = {name: center for name, center in ESF_CLUSTER_CENTERS}
    names = tuple(centers)
    for left, right in product(names, names):
        if left >= right:
            continue
        gap = abs(centers[left] - centers[right])
        if gap <= radii[left] + radii[right]:
            raise SGBLConeStop(
                "resource_limit",
                "declared ESF contours are not strictly disjoint",
                {"left": left, "right": right, "gap": gap},
            )


def _zero_tensor4() -> tuple[tuple[tuple[tuple[Interval, ...], ...], ...], ...]:
    return tuple(
        tuple(tuple(tuple(ZERO for _ in range(DIMENSION)) for _ in range(DIMENSION))
              for _ in range(DIMENSION))
        for _ in range(DIMENSION)
    )


def _zero_matrix4() -> tuple[tuple[Interval, ...], ...]:
    return tuple(tuple(ZERO for _ in range(DIMENSION)) for _ in range(DIMENSION))


def _inflate(entry: Interval, width: Fraction) -> Interval:
    if width == 0:
        return entry
    return Interval(entry.lower - width, entry.upper + width)


def _interval_tensor4(
    name: str,
    value: object,
    *,
    half_width: Fraction | None = None,
) -> tuple[tuple[tuple[tuple[Interval, ...], ...], ...], ...]:
    width = Q(0) if half_width is None else half_width
    if value is None:
        if width == 0:
            return _zero_tensor4()
        bump = interval(-width, width)
        return tuple(
            tuple(
                tuple(tuple(bump for _ in range(DIMENSION)) for _ in range(DIMENSION))
                for _ in range(DIMENSION)
            )
            for _ in range(DIMENSION)
        )
    if type(value) is not tuple or len(value) != DIMENSION:
        raise TypeError(f"{name} must be a 4x4x4x4 nested tuple of intervals")
    rows = []
    for first in value:
        if type(first) is not tuple or len(first) != DIMENSION:
            raise TypeError(f"{name} must be a 4x4x4x4 nested tuple of intervals")
        plane = []
        for second in first:
            if type(second) is not tuple or len(second) != DIMENSION:
                raise TypeError(f"{name} must be a 4x4x4x4 nested tuple of intervals")
            line = []
            for third in second:
                if type(third) is not tuple or len(third) != DIMENSION:
                    raise TypeError(f"{name} must be a 4x4x4x4 nested tuple")
                line.append(
                    tuple(_inflate(_interval(name, entry), width) for entry in third)
                )
            plane.append(tuple(line))
        rows.append(tuple(plane))
    return tuple(rows)


def _interval_matrix4(
    name: str,
    value: object,
    *,
    half_width: Fraction | None = None,
) -> tuple[tuple[Interval, ...], ...]:
    width = Q(0) if half_width is None else half_width
    if value is None:
        if width == 0:
            return _zero_matrix4()
        bump = interval(-width, width)
        matrix = tuple(tuple(bump for _ in range(DIMENSION)) for _ in range(DIMENSION))
    else:
        if type(value) is not tuple or len(value) != DIMENSION:
            raise TypeError(f"{name} must be a 4x4 nested tuple of intervals")
        matrix = tuple(
            tuple(_inflate(_interval(name, entry), width) for entry in row)
            for row in value
        )
        if any(len(row) != DIMENSION for row in matrix):
            raise TypeError(f"{name} must be a 4x4 nested tuple of intervals")
    for first in range(DIMENSION):
        for second in range(first + 1, DIMENSION):
            if matrix[first][second] != matrix[second][first]:
                raise ValueError(f"{name} must be symmetric as an interval matrix")
    return matrix


def _all_zero_tensor4(value: Sequence[Sequence[Sequence[Sequence[Interval]]]]) -> bool:
    return all(
        entry == ZERO
        for first in value
        for second in first
        for third in second
        for entry in third
    )


def _all_zero_matrix4(value: Sequence[Sequence[Interval]]) -> bool:
    return all(entry == ZERO for row in value for entry in row)


@dataclass(frozen=True, slots=True, kw_only=True)
class SGBLPrincipalBackgroundBox:
    """Caller-declared orthonormal principal-background product box.

    The linear-branch constraints are part of the contract: ``F=Mpl^2``,
    ``F'=0``, ``f'=alpha_gb``, and ``Hess(f)=alpha_gb Hess(phi)``.  Widths are
    declared by the caller; they are not fitted from a later margin.
    """

    planck_mass: Interval | Fraction | int
    alpha_gb: Interval | Fraction | int
    riemann_lower: tuple | None = None
    hessian_phi_lower: tuple | None = None
    riemann_half_width: Fraction | int = 0
    hessian_phi_half_width: Fraction | int = 0
    tilde_normal_factor: Fraction | int = TILDE_NORMAL_FACTOR
    hat_normal_factor: Fraction | int = HAT_NORMAL_FACTOR
    limits: SGBLConeLimits = SGBLConeLimits()
    allow_zero_coupling_control: bool = False
    chart: str = ORTHONORMAL_MINKOWSKI_CHART

    def __post_init__(self) -> None:
        if type(self.limits) is not SGBLConeLimits:
            raise TypeError("limits must be SGBLConeLimits")
        if type(self.allow_zero_coupling_control) is not bool:
            raise TypeError("allow_zero_coupling_control must be boolean")
        if self.chart != ORTHONORMAL_MINKOWSKI_CHART:
            raise SGBLConeStop(
                "cone_chart_error",
                "SGB-L cone enclosure is the orthonormal Minkowski principal frame",
            )
        object.__setattr__(self, "planck_mass", _interval("planck_mass", self.planck_mass))
        object.__setattr__(self, "alpha_gb", _interval("alpha_gb", self.alpha_gb))
        object.__setattr__(
            self,
            "riemann_half_width",
            _fraction("riemann_half_width", self.riemann_half_width, nonnegative=True),
        )
        object.__setattr__(
            self,
            "hessian_phi_half_width",
            _fraction(
                "hessian_phi_half_width",
                self.hessian_phi_half_width,
                nonnegative=True,
            ),
        )
        object.__setattr__(
            self,
            "tilde_normal_factor",
            _fraction("tilde_normal_factor", self.tilde_normal_factor),
        )
        object.__setattr__(
            self,
            "hat_normal_factor",
            _fraction("hat_normal_factor", self.hat_normal_factor),
        )
        if not self.planck_mass.strictly_positive():
            raise SGBLConeStop(
                "sign_chart_error",
                "Planck mass interval must be strictly positive",
            )
        if self.alpha_gb.contains_zero() and not self.allow_zero_coupling_control:
            raise ValueError(
                "production SGB-L cone box requires nonzero alpha_gb; "
                "zero coupling is only a named control"
            )
        if not 1 < self.tilde_normal_factor < self.hat_normal_factor:
            raise SGBLConeStop(
                "sign_chart_error",
                "auxiliary normal factors must satisfy 1<tilde<hat",
            )
        if (
            self.tilde_normal_factor != DECLARED_TILDE_NORMAL_FACTOR
            or self.hat_normal_factor != DECLARED_HAT_NORMAL_FACTOR
        ):
            raise SGBLConeStop(
                "sign_chart_error",
                "SGB-L cone uses the locked tilde=4, hat=9 auxiliary cones",
            )
        width = max(self.riemann_half_width, self.hessian_phi_half_width)
        if width > self.limits.max_parameter_half_width:
            raise SGBLConeStop(
                "resource_limit",
                "declared principal-background half-width exceeds the arithmetic cap",
            )
        object.__setattr__(
            self,
            "riemann_lower",
            _interval_tensor4(
                "riemann_lower",
                self.riemann_lower,
                half_width=self.riemann_half_width,
            ),
        )
        object.__setattr__(
            self,
            "hessian_phi_lower",
            _interval_matrix4(
                "hessian_phi_lower",
                self.hessian_phi_lower,
                half_width=self.hessian_phi_half_width,
            ),
        )

    @property
    def effective_planck_squared(self) -> Interval:
        return self.planck_mass * self.planck_mass

    @property
    def effective_planck_prime(self) -> Interval:
        return ZERO

    @property
    def gb_coupling_prime(self) -> Interval:
        return self.alpha_gb

    @property
    def hessian_gb_lower(self) -> tuple[tuple[Interval, ...], ...]:
        return tuple(
            tuple(self.alpha_gb * entry for entry in row)
            for row in self.hessian_phi_lower
        )

    @property
    def is_singleton(self) -> bool:
        return (
            self.planck_mass.is_singleton()
            and self.alpha_gb.is_singleton()
            and self.riemann_half_width == 0
            and self.hessian_phi_half_width == 0
            and all(
                entry.is_singleton()
                for plane in self.riemann_lower
                for row in plane
                for line in row
                for entry in line
            )
            and all(
                entry.is_singleton()
                for row in self.hessian_phi_lower
                for entry in row
            )
        )

    @property
    def curvature_free(self) -> bool:
        return _all_zero_tensor4(self.riemann_lower) and _all_zero_matrix4(
            self.hessian_phi_lower
        )


def sgbl_esf_principal_box(
    *,
    planck_mass: Interval | Fraction | int = 2,
    alpha_gb: Interval | Fraction | int = Q(-1, 4),
    riemann_half_width: Fraction | int = 0,
    hessian_phi_half_width: Fraction | int = 0,
    limits: SGBLConeLimits | None = None,
    allow_zero_coupling_control: bool = False,
) -> SGBLPrincipalBackgroundBox:
    """Return a Minkowski-frame ESF box with the SGB-L derivative contract."""

    return SGBLPrincipalBackgroundBox(
        planck_mass=planck_mass,
        alpha_gb=alpha_gb,
        riemann_half_width=riemann_half_width,
        hessian_phi_half_width=hessian_phi_half_width,
        limits=limits or SGBLConeLimits(),
        allow_zero_coupling_control=allow_zero_coupling_control,
    )


def sgbl_principal_box_from_state(
    state: SphericalState,
    *,
    riemann_half_width: Fraction | int = 0,
    hessian_phi_half_width: Fraction | int = 0,
    limits: SGBLConeLimits | None = None,
) -> SGBLPrincipalBackgroundBox:
    """Ground a box in the SGB-L principal adapter, not in an FGC-QR builder.

    Float64 orthonormal curvature is converted to exact dyadics.  That is a
    declared-box contract on the adapter image, not a claim that the binary64
    curvature equals an exact spherical tensor.
    """

    background = sgbl_covariant_principal_background(state)
    alpha = Q(*float(background.gb_coupling_prime).as_integer_ratio())
    if alpha == 0:
        raise ValueError("adapter box requires nonzero alpha_gb")
    hessian_gb = tuple(
        tuple(
            Interval.singleton(Q(*float(background.hessian_gb_lower[a, b]).as_integer_ratio()))
            for b in range(DIMENSION)
        )
        for a in range(DIMENSION)
    )
    hessian_phi = tuple(
        tuple(hessian_gb[a][b] / Interval.singleton(alpha) for b in range(DIMENSION))
        for a in range(DIMENSION)
    )
    riemann = tuple(
        tuple(
            tuple(
                tuple(
                    Interval.singleton(
                        Q(*float(background.riemann_lower[a, b, c, d]).as_integer_ratio())
                    )
                    for d in range(DIMENSION)
                )
                for c in range(DIMENSION)
            )
            for b in range(DIMENSION)
        )
        for a in range(DIMENSION)
    )
    mass = Q(*float(state.planck_mass).as_integer_ratio()) if not isinstance(
        state.planck_mass, (int, Fraction)
    ) else Fraction(state.planck_mass)
    box = SGBLPrincipalBackgroundBox(
        planck_mass=mass,
        alpha_gb=alpha,
        riemann_lower=riemann,
        hessian_phi_lower=hessian_phi,
        riemann_half_width=riemann_half_width,
        hessian_phi_half_width=hessian_phi_half_width,
        limits=limits or SGBLConeLimits(),
    )
    if box.effective_planck_prime != ZERO:
        raise SGBLConeStop("action_identity_error", "SGB-L box must keep F'=0")
    return box


def sgbl_nonflat_source_principal_box(
    *,
    riemann_half_width: Fraction | int = 0,
    hessian_phi_half_width: Fraction | int = 0,
    limits: SGBLConeLimits | None = None,
) -> SGBLPrincipalBackgroundBox:
    """Adapter image of the locked nonzero-shift two-scalar SGB-L source fixture."""

    point = SGBLSourceInputs(
        coordinate_radius=Q(5, 2),
        alpha=(2, Q(1, 3), Q(-1, 2), Q(1, 5), Q(-1, 7)),
        shift=(Q(1, 3), Q(-2, 5), Q(1, 4), Q(-1, 6), Q(2, 9)),
        radial_metric=(Q(3, 2), Q(1, 8), Q(-1, 3), Q(1, 7), Q(-1, 4)),
        areal_radius=(4, Q(-1, 5), Q(3, 2), Q(1, 9), Q(-2, 5)),
        phi=(Q(1, 2), Q(2, 3), Q(-4, 5), Q(1, 4), Q(-1, 8)),
        chi=(Q(-3, 4), Q(1, 2), Q(3, 5), Q(-2, 7), Q(1, 6)),
        planck_mass=2,
        scalar_mass=3,
        quartic_coupling=Q(1, 2),
        alpha_gb=Q(-1, 4),
    )
    return sgbl_principal_box_from_state(
        sgbl_source_state(point, sgbl_source_solve(point)),
        riemann_half_width=riemann_half_width,
        hessian_phi_half_width=hessian_phi_half_width,
        limits=limits,
    )


def _symmetric_basis(component: tuple[int, int]) -> list[list[Interval]]:
    first, second = component
    answer = [[ZERO for _ in range(DIMENSION)] for _ in range(DIMENSION)]
    answer[first][second] = ONE
    answer[second][first] = ONE
    return answer


def _variational(value: Sequence[Sequence[Interval]], component: tuple[int, int]) -> Interval:
    first, second = component
    entry = value[first][second]
    return entry if first == second else TWO * entry


def _raise_matrix(matrix: Sequence[Sequence[Interval]]) -> list[list[Interval]]:
    return [
        [Interval.singleton(ETA_SIGNS[a] * ETA_SIGNS[b]) * matrix[a][b] for b in range(DIMENSION)]
        for a in range(DIMENSION)
    ]


def _principal_riemann(
    perturbation: Sequence[Sequence[Interval]],
    xi: Sequence[Interval],
) -> list[list[list[list[Interval]]]]:
    answer = [
        [[[ZERO for _ in range(DIMENSION)] for _ in range(DIMENSION)] for _ in range(DIMENSION)]
        for _ in range(DIMENSION)
    ]
    for a, b, c, d in product(range(DIMENSION), repeat=4):
        answer[a][b][c][d] = HALF * (
            xi[b] * xi[c] * perturbation[a][d]
            + xi[a] * xi[d] * perturbation[b][c]
            - xi[b] * xi[d] * perturbation[a][c]
            - xi[a] * xi[c] * perturbation[b][d]
        )
    return answer


def _ricci_and_scalar(
    riemann: Sequence[Sequence[Sequence[Sequence[Interval]]]],
) -> tuple[list[list[Interval]], Interval]:
    ricci = [[ZERO for _ in range(DIMENSION)] for _ in range(DIMENSION)]
    for b, d in product(range(DIMENSION), repeat=2):
        total = ZERO
        for a in range(DIMENSION):
            total = total + Interval.singleton(ETA_SIGNS[a]) * riemann[a][b][a][d]
        ricci[b][d] = total
    scalar = ZERO
    for b in range(DIMENSION):
        scalar = scalar + Interval.singleton(ETA_SIGNS[b]) * ricci[b][b]
    return ricci, scalar


def _double_dual(
    riemann: Sequence[Sequence[Sequence[Sequence[Interval]]]],
    ricci: Sequence[Sequence[Interval]],
    scalar: Interval,
) -> list[list[list[list[Interval]]]]:
    metric = tuple(
        tuple(Interval.singleton(ETA_SIGNS[a] if a == b else 0) for b in range(DIMENSION))
        for a in range(DIMENSION)
    )
    answer = [
        [[[ZERO for _ in range(DIMENSION)] for _ in range(DIMENSION)] for _ in range(DIMENSION)]
        for _ in range(DIMENSION)
    ]
    for a, b, c, d in product(range(DIMENSION), repeat=4):
        answer[a][b][c][d] = (
            riemann[a][b][c][d]
            - metric[a][c] * ricci[d][b]
            + metric[a][d] * ricci[c][b]
            + metric[b][c] * ricci[d][a]
            - metric[b][d] * ricci[c][a]
            + HALF * scalar * (metric[a][c] * metric[d][b] - metric[a][d] * metric[c][b])
        )
    return answer


def _gauss_bonnet_variation(
    background_riemann: Sequence[Sequence[Sequence[Sequence[Interval]]]],
    principal_riemann: Sequence[Sequence[Sequence[Sequence[Interval]]]],
) -> Interval:
    background_ricci, background_scalar = _ricci_and_scalar(background_riemann)
    principal_ricci, principal_scalar = _ricci_and_scalar(principal_riemann)
    background_ricci_up = _raise_matrix(background_ricci)
    total = TWO * background_scalar * principal_scalar
    ricci_pair = ZERO
    for a, b in product(range(DIMENSION), repeat=2):
        ricci_pair = ricci_pair + background_ricci_up[a][b] * principal_ricci[a][b]
    riemann_pair = ZERO
    for a, b, c, d in product(range(DIMENSION), repeat=4):
        sign = ETA_SIGNS[a] * ETA_SIGNS[b] * ETA_SIGNS[c] * ETA_SIGNS[d]
        riemann_pair = riemann_pair + Interval.singleton(sign) * background_riemann[a][b][c][d] * principal_riemann[a][b][c][d]
    return total - FOUR * TWO * ricci_pair + TWO * riemann_pair


def _ungauged_interval_symbol(
    box: SGBLPrincipalBackgroundBox,
    xi: Sequence[Interval],
) -> list[list[Interval]]:
    xi_up = tuple(Interval.singleton(ETA_SIGNS[a]) * xi[a] for a in range(DIMENSION))
    xi_squared = sum((xi[a] * xi_up[a] for a in range(DIMENSION)), ZERO)
    riemann = box.riemann_lower
    hessian_gb = box.hessian_gb_lower
    f_prime = box.gb_coupling_prime
    planck = box.effective_planck_squared
    riemann_zero = _all_zero_tensor4(riemann)
    hessian_zero = _all_zero_matrix4(hessian_gb)
    if riemann_zero:
        ricci = [[ZERO for _ in range(DIMENSION)] for _ in range(DIMENSION)]
        scalar = ZERO
        double_dual = None
    else:
        ricci, scalar = _ricci_and_scalar(riemann)
        double_dual = _double_dual(riemann, ricci, scalar)
    hessian_gb_up = (
        [[ZERO for _ in range(DIMENSION)] for _ in range(DIMENSION)]
        if hessian_zero
        else _raise_matrix(hessian_gb)
    )
    answer = [[ZERO for _ in range(FIELD_COUNT)] for _ in range(FIELD_COUNT)]
    metric = tuple(
        tuple(Interval.singleton(ETA_SIGNS[a] if a == b else 0) for b in range(DIMENSION))
        for a in range(DIMENSION)
    )
    metric_scalar_cross = [[ZERO for _ in range(DIMENSION)] for _ in range(DIMENSION)]
    if not riemann_zero and double_dual is not None:
        for a, b in product(range(DIMENSION), repeat=2):
            contraction = ZERO
            for c, d in product(range(DIMENSION), repeat=2):
                contraction = contraction + double_dual[a][c][b][d] * xi_up[c] * xi_up[d]
            metric_scalar_cross[a][b] = FOUR * f_prime * contraction

    for column, component in enumerate(METRIC_COMPONENTS):
        perturbation = _symmetric_basis(component)
        delta_riemann = _principal_riemann(perturbation, xi)
        delta_ricci, delta_scalar = _ricci_and_scalar(delta_riemann)
        delta_einstein = [
            [
                delta_ricci[a][b] - HALF * metric[a][b] * delta_scalar
                for b in range(DIMENSION)
            ]
            for a in range(DIMENSION)
        ]
        hess_term = [[ZERO for _ in range(DIMENSION)] for _ in range(DIMENSION)]
        if not hessian_zero:
            delta_double_dual = _double_dual(delta_riemann, delta_ricci, delta_scalar)
            for a, b in product(range(DIMENSION), repeat=2):
                total = ZERO
                for c, d in product(range(DIMENSION), repeat=2):
                    total = total + delta_double_dual[a][c][b][d] * hessian_gb_up[c][d]
                hess_term[a][b] = FOUR * total
        response_lower = [
            [HALF * planck * delta_einstein[a][b] + hess_term[a][b] for b in range(DIMENSION)]
            for a in range(DIMENSION)
        ]
        response = _raise_matrix(response_lower)
        response = [[-response[a][b] for b in range(DIMENSION)] for a in range(DIMENSION)]
        for row, equation in enumerate(METRIC_COMPONENTS):
            answer[row][column] = _variational(response, equation)
        scalar_metric = ZERO
        if not riemann_zero:
            scalar_metric = f_prime * _gauss_bonnet_variation(riemann, delta_riemann)
        answer[10][column] = scalar_metric

    raised_cross = _raise_matrix(metric_scalar_cross)
    raised_cross = [[-raised_cross[a][b] for b in range(DIMENSION)] for a in range(DIMENSION)]
    for row, equation in enumerate(METRIC_COMPONENTS):
        answer[row][10] = _variational(raised_cross, equation)
    answer[10][10] = xi_squared
    answer[11][11] = xi_squared
    return answer


def _diag_projector_xi(
    diag: Sequence[Fraction],
    alpha: int,
    mu: int,
    nu: int,
    xi: Sequence[Interval],
) -> Interval:
    """Contract a diagonal trace-reversal projector against ``xi``."""

    hat_nu = Interval.singleton(diag[nu])
    hat_mu = Interval.singleton(diag[mu])
    term = ZERO
    if alpha == mu:
        term = term + hat_nu * xi[nu]
    if alpha == nu:
        term = term + hat_mu * xi[mu]
    if mu == nu:
        term = term - hat_mu * xi[alpha]
    return HALF * term


def _diag_projector_h(
    diag: Sequence[Fraction],
    beta: int,
    xi: Sequence[Interval],
    perturbation: Sequence[Sequence[Interval]],
) -> Interval:
    """Contract a diagonal trace-reversal projector against ``xi`` and ``h``."""

    total = ZERO
    for rho, sigma in product(range(DIMENSION), repeat=2):
        entry = perturbation[rho][sigma]
        if entry == ZERO:
            continue
        for delta in range(DIMENSION):
            term = ZERO
            if beta == rho:
                term = term + Interval.singleton(diag[sigma] if sigma == delta else 0) * xi[delta]
            if beta == sigma:
                term = term + Interval.singleton(diag[rho] if rho == delta else 0) * xi[delta]
            if beta == delta:
                term = term - Interval.singleton(diag[rho] if rho == sigma else 0) * xi[delta]
            total = total + HALF * term * entry
    return total


def _gauge_interval_symbol(
    box: SGBLPrincipalBackgroundBox,
    xi: Sequence[Interval],
) -> list[list[Interval]]:
    tilde_diag = (-box.tilde_normal_factor, Q(1), Q(1), Q(1))
    hat_diag = (-box.hat_normal_factor, Q(1), Q(1), Q(1))
    scale = HALF * box.effective_planck_squared
    answer = [[ZERO for _ in range(FIELD_COUNT)] for _ in range(FIELD_COUNT)]
    for column, component in enumerate(METRIC_COMPONENTS):
        perturbation = _symmetric_basis(component)
        response_up = [[ZERO for _ in range(DIMENSION)] for _ in range(DIMENSION)]
        tilde_by_beta = [
            _diag_projector_h(tilde_diag, beta, xi, perturbation)
            for beta in range(DIMENSION)
        ]
        for mu, nu in product(range(DIMENSION), repeat=2):
            total = ZERO
            for alpha in range(DIMENSION):
                hat_xi = _diag_projector_xi(hat_diag, alpha, mu, nu, xi)
                eta = Interval.singleton(ETA_SIGNS[alpha])
                total = total + hat_xi * eta * tilde_by_beta[alpha]
            response_up[mu][nu] = scale * total
        for row, equation in enumerate(METRIC_COMPONENTS):
            answer[row][column] = _variational(response_up, equation)
    return answer


def _add_matrices(
    left: Sequence[Sequence[Interval]],
    right: Sequence[Sequence[Interval]],
) -> list[list[Interval]]:
    return [
        [left[row][column] + right[row][column] for column in range(len(left[0]))]
        for row in range(len(left))
    ]


def _frobenius_normalize(
    matrix: Sequence[Sequence[Interval]],
    sqrt2: Interval,
) -> list[list[Interval]]:
    scales = []
    for first, second in METRIC_COMPONENTS:
        scales.append(sqrt2 if first != second else ONE)
    scales.extend((ONE, ONE))
    return [
        [matrix[row][column] / (scales[row] * scales[column]) for column in range(FIELD_COUNT)]
        for row in range(FIELD_COUNT)
    ]


def _covector(*components: int | Fraction) -> tuple[Interval, ...]:
    return tuple(Interval.singleton(Q(component)) for component in components)


class _EvaluationCounter:
    def __init__(self, limit: int) -> None:
        self.limit = limit
        self.count = 0
        self.bits = 0

    def bump(self, matrix: Sequence[Sequence[Interval]], *, bit_limit: int) -> None:
        self.count += 1
        if self.count > self.limit:
            raise SGBLConeStop(
                "resource_limit",
                "principal-symbol evaluations exceeded the declared cap",
                {"count": self.count, "limit": self.limit},
            )
        observed = max(_interval_bits(entry) for row in matrix for entry in row)
        self.bits = max(self.bits, observed)
        if observed > bit_limit:
            raise SGBLConeStop(
                "resource_limit",
                "interval endpoints exceeded the declared rational bit cap",
                {"observed": observed, "limit": bit_limit},
            )


def _complete_symbol(
    box: SGBLPrincipalBackgroundBox,
    xi: Sequence[Interval],
    *,
    include_gauge: bool,
    sqrt2: Interval,
    counter: _EvaluationCounter,
    normalize: bool = True,
) -> list[list[Interval]]:
    raw = _ungauged_interval_symbol(box, xi)
    if include_gauge:
        raw = _add_matrices(raw, _gauge_interval_symbol(box, xi))
    matrix = _frobenius_normalize(raw, sqrt2) if normalize else raw
    counter.bump(matrix, bit_limit=box.limits.max_rational_bit_length)
    return matrix


def _clone_matrix(matrix: Sequence[Sequence[Interval]]) -> list[list[Interval]]:
    return [[entry for entry in row] for row in matrix]


def _scale_metric_block(
    matrix: Sequence[Sequence[Interval]],
    factor: Interval,
) -> list[list[Interval]]:
    answer = _clone_matrix(matrix)
    for row in range(10):
        for column in range(10):
            answer[row][column] = factor * matrix[row][column]
    return answer


@lru_cache(maxsize=8)
def _unit_esf_coefficient_cache(
    tilde: Fraction,
    hat: Fraction,
    include_gauge: bool,
    sqrt2_bits: int,
    normalize: bool,
) -> tuple:
    unit_box = SGBLPrincipalBackgroundBox(
        planck_mass=1,
        alpha_gb=0,
        limits=SGBLConeLimits(
            sqrt2_bits=sqrt2_bits,
            max_residual_evaluations=256,
            arc_count=DEFAULT_ARC_COUNT,
        ),
        allow_zero_coupling_control=True,
    )
    counter = _EvaluationCounter(256)
    extracted = _extract_coefficient_tensors_general(
        unit_box,
        include_gauge=include_gauge,
        sqrt2=_sqrt2_enclosure(sqrt2_bits),
        counter=counter,
        normalize=normalize,
    )
    frozen = (
        _freeze_matrix(extracted["time_time"]),
        tuple(_freeze_matrix(matrix) for matrix in extracted["time_space"]),
        tuple(
            tuple(
                _freeze_matrix(extracted["space_space"][first][second])
                for second in range(3)
            )
            for first in range(3)
        ),
    )
    return frozen


def _curvature_free_coefficients(
    box: SGBLPrincipalBackgroundBox,
    *,
    include_gauge: bool,
    sqrt2: Interval,
    counter: _EvaluationCounter,
    normalize: bool,
) -> dict[str, Any]:
    _ = sqrt2
    unit_time, unit_space, unit_ss = _unit_esf_coefficient_cache(
        box.tilde_normal_factor,
        box.hat_normal_factor,
        include_gauge,
        box.limits.sqrt2_bits,
        normalize,
    )
    factor = box.effective_planck_squared
    time_time = _scale_metric_block(unit_time, factor)
    time_space = [_scale_metric_block(matrix, factor) for matrix in unit_space]
    space_space = [
        [_scale_metric_block(unit_ss[first][second], factor) for second in range(3)]
        for first in range(3)
    ]
    for matrix in (time_time, *time_space, *(space_space[i][j] for i in range(3) for j in range(3))):
        counter.bump(matrix, bit_limit=box.limits.max_rational_bit_length)
    return {
        "time_time": time_time,
        "time_space": time_space,
        "space_space": space_space,
        "origin_unused": _covector(0, 0, 0, 0),
    }


def _extract_coefficient_tensors(
    box: SGBLPrincipalBackgroundBox,
    *,
    include_gauge: bool,
    sqrt2: Interval,
    counter: _EvaluationCounter,
    normalize: bool = True,
) -> dict[str, Any]:
    if box.curvature_free:
        return _curvature_free_coefficients(
            box,
            include_gauge=include_gauge,
            sqrt2=sqrt2,
            counter=counter,
            normalize=normalize,
        )
    return _extract_coefficient_tensors_general(
        box,
        include_gauge=include_gauge,
        sqrt2=sqrt2,
        counter=counter,
        normalize=normalize,
    )


def _extract_coefficient_tensors_general(
    box: SGBLPrincipalBackgroundBox,
    *,
    include_gauge: bool,
    sqrt2: Interval,
    counter: _EvaluationCounter,
    normalize: bool = True,
) -> dict[str, Any]:
    origin = _covector(0, 0, 0, 0)
    time = _covector(1, 0, 0, 0)
    time_time = _complete_symbol(
        box,
        time,
        include_gauge=include_gauge,
        sqrt2=sqrt2,
        counter=counter,
        normalize=normalize,
    )
    time_space: list[list[list[Interval]]] = []
    space_space = [
        [[[ZERO for _ in range(FIELD_COUNT)] for _ in range(FIELD_COUNT)] for _ in range(3)]
        for _ in range(3)
    ]
    spatial = []
    for index in range(3):
        components = [0, 0, 0, 0]
        components[index + 1] = 1
        vector = _covector(*components)
        spatial.append(vector)
        space_space[index][index] = _complete_symbol(
            box,
            vector,
            include_gauge=include_gauge,
            sqrt2=sqrt2,
            counter=counter,
            normalize=normalize,
        )
        plus = _covector(*tuple(1 if slot == 0 else components[slot] for slot in range(4)))
        minus = _covector(*tuple(-1 if slot == 0 else components[slot] for slot in range(4)))
        plus_symbol = _complete_symbol(
            box,
            plus,
            include_gauge=include_gauge,
            sqrt2=sqrt2,
            counter=counter,
            normalize=normalize,
        )
        minus_symbol = _complete_symbol(
            box,
            minus,
            include_gauge=include_gauge,
            sqrt2=sqrt2,
            counter=counter,
            normalize=normalize,
        )
        time_space.append(
            [
                [HALF * (plus_symbol[row][column] - minus_symbol[row][column])
                 for column in range(FIELD_COUNT)]
                for row in range(FIELD_COUNT)
            ]
        )
    for first in range(3):
        for second in range(first + 1, 3):
            plus_components = [0, 0, 0, 0]
            minus_components = [0, 0, 0, 0]
            plus_components[first + 1] = 1
            plus_components[second + 1] = 1
            minus_components[first + 1] = 1
            minus_components[second + 1] = -1
            plus_symbol = _complete_symbol(
                box,
                _covector(*plus_components),
                include_gauge=include_gauge,
                sqrt2=sqrt2,
                counter=counter,
                normalize=normalize,
            )
            minus_symbol = _complete_symbol(
                box,
                _covector(*minus_components),
                include_gauge=include_gauge,
                sqrt2=sqrt2,
                counter=counter,
                normalize=normalize,
            )
            mixed = [
                [
                    Interval.singleton(Q(1, 4))
                    * (plus_symbol[row][column] - minus_symbol[row][column])
                    for column in range(FIELD_COUNT)
                ]
                for row in range(FIELD_COUNT)
            ]
            space_space[first][second] = mixed
            space_space[second][first] = mixed
    return {
        "time_time": time_time,
        "time_space": time_space,
        "space_space": space_space,
        "origin_unused": origin,
    }


def _freeze_matrix(matrix: Sequence[Sequence[Interval]]) -> tuple[tuple[Interval, ...], ...]:
    return tuple(tuple(entry for entry in row) for row in matrix)


def _midpoint_matrix(matrix: Sequence[Sequence[Interval]]) -> tuple[tuple[Fraction, ...], ...]:
    return tuple(tuple(entry.midpoint() for entry in row) for row in matrix)


def _contains_float_matrix(
    enclosure: Sequence[Sequence[Interval]],
    array: np.ndarray,
) -> bool:
    for row in range(FIELD_COUNT):
        for column in range(FIELD_COUNT):
            dyadic = Q(*float(array[row, column]).as_integer_ratio())
            entry = enclosure[row][column]
            scale = max(Q(1), abs(entry.lower), abs(entry.upper), abs(dyadic))
            slack = DEFAULT_ROUND_UNITS * BINARY64_UNIT * scale
            if dyadic < entry.lower - slack or dyadic > entry.upper + slack:
                return False
    return True


def _symmetry_defect(matrix: Sequence[Sequence[Interval]]) -> Interval:
    lower = min(
        (matrix[row][column] - matrix[column][row]).lower
        for row in range(len(matrix))
        for column in range(len(matrix))
    )
    upper = max(
        (matrix[row][column] - matrix[column][row]).upper
        for row in range(len(matrix))
        for column in range(len(matrix))
    )
    return Interval(lower, upper)


def _quadratic_symbol_tensor(
    coefficients: Mapping[str, Any],
) -> list:
    quadratic = [
        [
            [[ZERO for _ in range(FIELD_COUNT)] for _ in range(FIELD_COUNT)]
            for _ in range(DIMENSION)
        ]
        for _ in range(DIMENSION)
    ]
    quadratic[0][0] = coefficients["time_time"]
    for index in range(3):
        half = [
            [
                HALF * coefficients["time_space"][index][row][column]
                for column in range(FIELD_COUNT)
            ]
            for row in range(FIELD_COUNT)
        ]
        quadratic[0][index + 1] = half
        quadratic[index + 1][0] = half
    for first in range(3):
        for second in range(3):
            quadratic[first + 1][second + 1] = coefficients["space_space"][first][second]
    return quadratic


def _pure_gauge_field_coefficients(*, mutate_sign: int = 1) -> list:
    if mutate_sign not in (-1, 1):
        raise ValueError("gauge mutation sign must be +1 or -1")
    gauge = [
        [[ZERO for _ in range(DIMENSION)] for _ in range(FIELD_COUNT)]
        for _ in range(DIMENSION)
    ]
    for mu in range(DIMENSION):
        for field, (first, second) in enumerate(METRIC_COMPONENTS):
            for vector in range(DIMENSION):
                gauge[mu][field][vector] = Interval.singleton(
                    int(first == mu and second == vector)
                    + mutate_sign * int(second == mu and first == vector)
                )
    return gauge


def _pure_gauge_operator_sum(
    coefficients: Mapping[str, Any],
    bits: int,
    *,
    mutate_sign: int = 1,
) -> Interval:
    """n-independent cubic bound of ``P(xi) h(xi,X)``.

    Unnormalized field coordinates keep the identity free of ``sqrt(2)``
    wrapping.  The triangle inequality on the symmetrized 12-vectors bounds
    every unit covector.
    """

    if bits <= 0:
        raise ValueError("sqrt2_bits must be positive")
    quadratic = _quadratic_symbol_tensor(coefficients)
    gauge = _pure_gauge_field_coefficients(mutate_sign=mutate_sign)
    permutations = (
        (0, 1, 2),
        (0, 2, 1),
        (1, 0, 2),
        (1, 2, 0),
        (2, 0, 1),
        (2, 1, 0),
    )
    total_upper = Q(0)
    identity_lower = Q(0)
    six = Interval.singleton(6)
    for first, second, third in product(range(DIMENSION), repeat=3):
        slots = (first, second, third)
        for vector in range(DIMENSION):
            accum = [ZERO for _ in range(FIELD_COUNT)]
            for perm in permutations:
                left, mid, right = slots[perm[0]], slots[perm[1]], slots[perm[2]]
                for row in range(FIELD_COUNT):
                    entry = ZERO
                    for inner in range(FIELD_COUNT):
                        entry = (
                            entry
                            + quadratic[left][mid][row][inner]
                            * gauge[right][inner][vector]
                        )
                    accum[row] = accum[row] + entry
            for row in range(FIELD_COUNT):
                entry = accum[row] / six
                identity_lower = min(identity_lower, entry.lower)
                total_upper += entry.abs_upper()
    if identity_lower > 0:
        return Interval(identity_lower, total_upper)
    return Interval(min(Q(0), identity_lower), total_upper)


def _action_energy_deformation(
    action: Mapping[str, Any],
    reference: Mapping[str, Any],
    bits: int,
) -> Fraction:
    delta_a = [
        [
            action["time_time"][row][column] - reference["time_time"][row][column]
            for column in range(FIELD_COUNT)
        ]
        for row in range(FIELD_COUNT)
    ]
    bound = _spectral_norm_upper(delta_a, bits)
    for index in range(3):
        delta_b = [
            [
                action["time_space"][index][row][column]
                - reference["time_space"][index][row][column]
                for column in range(FIELD_COUNT)
            ]
            for row in range(FIELD_COUNT)
        ]
        bound += _spectral_norm_upper(delta_b, bits)
    return bound


def _companion_blocks(
    coefficients: Mapping[str, Any],
    inverse_a: Sequence[Sequence[Interval]],
) -> tuple[list[list[list[Interval]]], list[list[list[list[Interval]]]]]:
    time_space = []
    for index in range(3):
        applied = [
            [
                sum(
                    (
                        inverse_a[row][inner] * coefficients["time_space"][index][inner][column]
                        for inner in range(FIELD_COUNT)
                    ),
                    ZERO,
                )
                for column in range(FIELD_COUNT)
            ]
            for row in range(FIELD_COUNT)
        ]
        time_space.append([[-entry for entry in row] for row in applied])
    space_space: list[list[list[list[Interval]]]] = [
        [[[ZERO for _ in range(FIELD_COUNT)] for _ in range(FIELD_COUNT)] for _ in range(3)]
        for _ in range(3)
    ]
    for first, second in product(range(3), repeat=2):
        applied = [
            [
                sum(
                    (
                        inverse_a[row][inner]
                        * coefficients["space_space"][first][second][inner][column]
                        for inner in range(FIELD_COUNT)
                    ),
                    ZERO,
                )
                for column in range(FIELD_COUNT)
            ]
            for row in range(FIELD_COUNT)
        ]
        space_space[first][second] = [[-entry for entry in row] for row in applied]
    return time_space, space_space


def _companion_deformation(
    time_space: Sequence[Sequence[Sequence[Interval]]],
    space_space: Sequence[Sequence[Sequence[Sequence[Interval]]]],
    reference_time_space: Sequence[Sequence[Sequence[Interval]]],
    reference_space_space: Sequence[Sequence[Sequence[Sequence[Interval]]]],
    bits: int,
) -> Fraction:
    bound = Q(0)
    for index in range(3):
        delta = [
            [
                time_space[index][row][column] - reference_time_space[index][row][column]
                for column in range(FIELD_COUNT)
            ]
            for row in range(FIELD_COUNT)
        ]
        bound += _spectral_norm_upper(delta, bits)
    for first, second in product(range(3), repeat=2):
        delta = [
            [
                space_space[first][second][row][column]
                - reference_space_space[first][second][row][column]
                for column in range(FIELD_COUNT)
            ]
            for row in range(FIELD_COUNT)
        ]
        bound += _spectral_norm_upper(delta, bits)
    return bound


def _interval_inverse(
    matrix: Sequence[Sequence[Interval]],
) -> Mapping[str, Any] | None:
    center = _midpoint_matrix(matrix)
    try:
        center_inverse = inverse(center)
    except ValueError as exc:
        raise SGBLConeStop(
            "lost_kinetic",
            "center kinetic block A is exactly singular",
        ) from exc
    try:
        return center_preconditioned_neumann_inverse(matrix, center_inverse)
    except ValueError:
        return None


def _as_rational_matrix(
    matrix: Sequence[Sequence[Interval]],
    *,
    name: str,
) -> tuple[tuple[Fraction, ...], ...]:
    rows = []
    for row in matrix:
        converted = []
        for entry in row:
            if not entry.is_singleton():
                raise SGBLConeStop(
                    "interval_inconclusive",
                    f"{name} is not an exact rational singleton",
                )
            converted.append(entry.lower)
        rows.append(tuple(converted))
    return tuple(rows)


def _one_norm_rational(matrix: Sequence[Sequence[Fraction]]) -> Fraction:
    size = len(matrix)
    return max(sum(abs(matrix[row][column]) for row in range(size)) for column in range(size))


def _inf_norm_rational(matrix: Sequence[Sequence[Fraction]]) -> Fraction:
    return max(sum(abs(entry) for entry in row) for row in matrix)


def _spectral_norm_upper_rational(
    matrix: Sequence[Sequence[Fraction]],
    bits: int,
) -> Fraction:
    return _rational_sqrt_enclosure(
        _one_norm_rational(matrix) * _inf_norm_rational(matrix),
        bits,
    ).upper


def _negate_rational(
    matrix: Sequence[Sequence[Fraction]],
) -> tuple[tuple[Fraction, ...], ...]:
    return tuple(tuple(-entry for entry in row) for row in matrix)


def _block_companion(
    time_space: Sequence[Sequence[Fraction]],
    space_space: Sequence[Sequence[Fraction]],
) -> tuple[tuple[Fraction, ...], ...]:
    identity = tuple(
        tuple(Q(int(row == column)) for column in range(FIELD_COUNT))
        for row in range(FIELD_COUNT)
    )
    zero = tuple(tuple(Q(0) for _ in range(FIELD_COUNT)) for _ in range(FIELD_COUNT))
    top = tuple(zero[row] + identity[row] for row in range(FIELD_COUNT))
    bottom = tuple(
        _negate_rational(space_space)[row] + _negate_rational(time_space)[row]
        for row in range(FIELD_COUNT)
    )
    return top + bottom


def _exact_radial_companion(
    time_time: Sequence[Sequence[Fraction]],
    time_space: Sequence[Sequence[Fraction]],
    space_space: Sequence[Sequence[Fraction]],
) -> tuple[tuple[Fraction, ...], ...]:
    try:
        inverse_a = inverse(time_time)
    except ValueError as exc:
        raise SGBLConeStop("lost_kinetic", "exact ESF kinetic block A is singular") from exc
    return _block_companion(
        matrix_multiply(inverse_a, time_space),
        matrix_multiply(inverse_a, space_space),
    )


def _exact_kernel(
    matrix: Sequence[Sequence[Fraction]],
    eigenvalue: Fraction,
) -> tuple[tuple[Fraction, ...], ...]:
    """Return an exact rational basis of ``ker(M - lambda I)``."""

    size = len(matrix)
    rows = [list(row) for row in matrix]
    for index in range(size):
        rows[index][index] -= eigenvalue
    pivots: list[int] = []
    rank_index = 0
    for column in range(size):
        pivot = next(
            (row for row in range(rank_index, size) if rows[row][column] != 0),
            None,
        )
        if pivot is None:
            continue
        rows[rank_index], rows[pivot] = rows[pivot], rows[rank_index]
        pivot_value = rows[rank_index][column]
        rows[rank_index] = [entry / pivot_value for entry in rows[rank_index]]
        for row in range(size):
            if row == rank_index:
                continue
            factor = rows[row][column]
            if factor:
                rows[row] = [
                    left - factor * right
                    for left, right in zip(rows[row], rows[rank_index])
                ]
        pivots.append(column)
        rank_index += 1
        if rank_index == size:
            break
    free = [column for column in range(size) if column not in set(pivots)]
    basis = []
    for free_column in free:
        vector = [Q(0)] * size
        vector[free_column] = Q(1)
        for pivot_row, pivot_column in enumerate(pivots):
            vector[pivot_column] = -rows[pivot_row][free_column]
        basis.append(tuple(vector))
    return tuple(basis)


def _columns_to_matrix(
    columns: Sequence[Sequence[Fraction]],
) -> tuple[tuple[Fraction, ...], ...]:
    height = len(columns[0])
    return tuple(tuple(column[row] for column in columns) for row in range(height))


def _leading_minors(
    matrix: Sequence[Sequence[Fraction]],
) -> tuple[Fraction, ...]:
    return tuple(
        determinant(tuple(row[:size] for row in matrix[:size]))
        for size in range(1, len(matrix) + 1)
    )


def _gershgorin_lower(matrix: Sequence[Sequence[Fraction]]) -> Fraction:
    return min(
        matrix[row][row]
        - sum(abs(matrix[row][column]) for column in range(len(matrix)) if column != row)
        for row in range(len(matrix))
    )


def _restrict_form(
    form: Sequence[Sequence[Fraction]],
    columns: Sequence[Sequence[Fraction]],
) -> tuple[tuple[Fraction, ...], ...]:
    frame = _columns_to_matrix(columns)
    return matrix_multiply(transpose(frame), matrix_multiply(form, frame))


def _unnormalized_esf_coefficients(
    *,
    planck_mass: Fraction,
    include_gauge: bool,
    limits: SGBLConeLimits,
) -> dict[str, Any]:
    box = sgbl_esf_principal_box(
        planck_mass=planck_mass,
        alpha_gb=0,
        limits=limits,
        allow_zero_coupling_control=True,
    )
    counter = _EvaluationCounter(limits.max_residual_evaluations)
    return _extract_coefficient_tensors(
        box,
        include_gauge=include_gauge,
        sqrt2=_sqrt2_enclosure(limits.sqrt2_bits),
        counter=counter,
        normalize=False,
    )


def _esf_exact_reference_certificate(
    *,
    planck_mass: Fraction,
    limits: SGBLConeLimits,
) -> dict[str, Any]:
    """Exact unnormalized radial ESF eigenframe, resolvent, and energy form.

    No floating SVD or fitted ulp allowance enters the certified bounds.
    """

    complete = _unnormalized_esf_coefficients(
        planck_mass=planck_mass, include_gauge=True, limits=limits
    )
    action = _unnormalized_esf_coefficients(
        planck_mass=planck_mass, include_gauge=False, limits=limits
    )
    time_time = _as_rational_matrix(complete["time_time"], name="ESF A")
    time_space = _as_rational_matrix(complete["time_space"][0], name="ESF B_r")
    space_space = _as_rational_matrix(
        complete["space_space"][0][0], name="ESF C_rr"
    )
    action_a = _as_rational_matrix(action["time_time"], name="ESF action A")
    action_b = _as_rational_matrix(action["time_space"][0], name="ESF action B_r")
    companion = _exact_radial_companion(time_time, time_space, space_space)
    columns: list[tuple[Fraction, ...]] = []
    clusters: dict[str, Any] = {}
    speeds = tuple(center for _name, center in ESF_CLUSTER_CENTERS)
    for name, center in ESF_CLUSTER_CENTERS:
        kernel = _exact_kernel(companion, center)
        if len(kernel) != 4:
            raise SGBLConeStop(
                "lost_kinetic",
                f"ESF kernel at speed {center} has dimension {len(kernel)}, not four",
            )
        for vector in kernel:
            image = matrix_multiply(companion, tuple((entry,) for entry in vector))
            expected = tuple((center * entry,) for entry in vector)
            if image != expected:
                raise SGBLConeStop(
                    "lost_kinetic",
                    f"ESF eigenvector identity failed at speed {center}",
                )
        columns.extend(kernel)
        radius = (
            limits.physical_contour_radius
            if name.startswith("physical")
            else limits.auxiliary_contour_radius
        )
        other_gap = min(abs(center - speed) for speed in speeds if speed != center)
        if radius * 2 >= other_gap:
            raise SGBLConeStop(
                "resource_limit",
                "declared contour radius is not strictly less than half the spectral gap",
                {"center": center, "radius": radius, "gap": other_gap},
            )
        clusters[name] = {
            "center": center,
            "radius": radius,
            "eigenvalue_count": 4,
            "exact_nullity": 4,
            "nearest_spectral_gap": other_gap,
            "contour_distance_lower": radius,
        }
    frame = _columns_to_matrix(columns)
    frame_rank = rank(frame)
    if frame_rank != FIRST_ORDER_FIELD_COUNT:
        raise SGBLConeStop(
            "lost_kinetic",
            f"ESF eigenframe rank is {frame_rank}, not 24",
        )
    try:
        frame_inverse = inverse(frame)
    except ValueError as exc:
        raise SGBLConeStop("lost_kinetic", "ESF eigenframe is singular") from exc
    kappa = _spectral_norm_upper_rational(
        frame, limits.sqrt2_bits
    ) * _spectral_norm_upper_rational(frame_inverse, limits.sqrt2_bits)
    if kappa <= 0:
        raise SGBLConeStop("resource_limit", "ESF eigenframe condition upper bound is nonpositive")
    for record in clusters.values():
        record["condition_number_2_upper"] = kappa
        record["certified_lower_bound"] = record["contour_distance_lower"] / kappa
    physical_energy: dict[str, Any] = {}
    energy_min = None
    form_norm_upper = Q(0)
    for name, target, slice_start in (
        ("physical_minus", Q(-1), 0),
        ("physical_plus", Q(1), 20),
    ):
        form = _block_rational_energy(
            action_b, action_a, sign=Q(1) if target < 0 else Q(-1)
        )
        restricted = _restrict_form(form, columns[slice_start:slice_start + 4])
        minors = _leading_minors(restricted)
        if any(minor <= 0 for minor in minors):
            raise SGBLConeStop(
                "lost_kinetic",
                f"ESF physical energy form at {name} fails Sylvester's criterion",
            )
        gershgorin = _gershgorin_lower(restricted)
        physical_energy[name] = {
            "dimension": 4,
            "sylvester_leading_minors": minors,
            "gershgorin_lower": gershgorin,
            "determinant": minors[-1],
        }
        energy_min = gershgorin if energy_min is None else min(energy_min, gershgorin)
        form_norm_upper = max(form_norm_upper, _spectral_norm_upper_rational(form, limits.sqrt2_bits))
    if energy_min is None or energy_min <= 0:
        raise SGBLConeStop(
            "lost_kinetic",
            "ESF physical energy Gershgorin lower bound is not strictly positive",
        )
    return {
        "classification": "exact_unnormalized_radial_eigenframe_and_Bauer_Fike_resolvent",
        "basis": "unnormalized_symmetric_tensor_components",
        "direction_representative": (1, 0, 0),
        "eigenframe_rank": frame_rank,
        "condition_number_2_upper": kappa,
        "clusters": clusters,
        "physical_energy": physical_energy,
        "physical_energy_coercivity_lower": energy_min,
        "physical_energy_operator_norm_upper": form_norm_upper,
        "six_contours_partition_all_24_characteristics": True,
        "proof_uses_floating_svd": False,
        "arc_samples_are_regression_only": True,
        "hyp2_sample_count_not_used": True,
        "ulp_allowance_not_used": True,
    }


def _block_rational_energy(
    action_b: Sequence[Sequence[Fraction]],
    action_a: Sequence[Sequence[Fraction]],
    *,
    sign: Fraction,
) -> tuple[tuple[Fraction, ...], ...]:
    """Return ``sign * [[B, A], [A, 0]]`` as a 24-by-24 rational matrix."""

    scaled_b = tuple(tuple(sign * entry for entry in row) for row in action_b)
    scaled_a = tuple(tuple(sign * entry for entry in row) for row in action_a)
    zero = tuple(tuple(Q(0) for _ in range(FIELD_COUNT)) for _ in range(FIELD_COUNT))
    top = tuple(scaled_b[row] + scaled_a[row] for row in range(FIELD_COUNT))
    bottom = tuple(scaled_a[row] + zero[row] for row in range(FIELD_COUNT))
    return top + bottom


def _exact_zero_interval(value: Interval) -> bool:
    return value.lower == 0 and value.upper == 0


def _extracted_field_symmetry(action: Mapping[str, Any]) -> Interval:
    symmetry = _symmetry_defect(action["time_time"])
    for index in range(3):
        other = _symmetry_defect(action["time_space"][index])
        symmetry = Interval(
            min(symmetry.lower, other.lower),
            max(symmetry.upper, other.upper),
        )
    for first, second in product(range(3), repeat=2):
        other = _symmetry_defect(action["space_space"][first][second])
        symmetry = Interval(
            min(symmetry.lower, other.lower),
            max(symmetry.upper, other.upper),
        )
    return symmetry


def _identities_on_action(
    action: Mapping[str, Any],
    bits: int,
) -> tuple[bool, bool, Interval, Interval]:
    symmetry = _extracted_field_symmetry(action)
    gauge = _pure_gauge_operator_sum(action, bits)
    return (
        _exact_zero_interval(symmetry),
        _exact_zero_interval(gauge),
        symmetry,
        gauge,
    )


def _rational_symmetric_matrix(component: tuple[int, int]) -> tuple[tuple[Fraction, ...], ...]:
    first, second = component
    data = [[Q(0) for _ in range(DIMENSION)] for _ in range(DIMENSION)]
    data[first][second] = Q(1)
    data[second][first] = Q(1)
    return tuple(tuple(row) for row in data)


def _interval_matrix_from_rational(
    matrix: Sequence[Sequence[Fraction]],
) -> tuple[tuple[Interval, ...], ...]:
    return tuple(tuple(Interval.singleton(entry) for entry in row) for row in matrix)


def _interval_tensor_from_rational(
    tensor: Sequence[Sequence[Sequence[Sequence[Fraction]]]],
) -> tuple[tuple[tuple[tuple[Interval, ...], ...], ...], ...]:
    return tuple(
        tuple(
            tuple(
                tuple(Interval.singleton(entry) for entry in line)
                for line in plane
            )
            for plane in volume
        )
        for volume in tensor
    )


def _flatten_matrix4(matrix: Sequence[Sequence[Fraction]]) -> tuple[Fraction, ...]:
    return tuple(matrix[row][column] for row in range(DIMENSION) for column in range(DIMENSION))


def _flatten_tensor4(
    tensor: Sequence[Sequence[Sequence[Sequence[Fraction]]]],
) -> tuple[Fraction, ...]:
    return tuple(
        tensor[first][second][third][fourth]
        for first, second, third, fourth in product(range(DIMENSION), repeat=4)
    )


def _kulkarni_nomizu(
    left: Sequence[Sequence[Fraction]],
    right: Sequence[Sequence[Fraction]],
) -> tuple[tuple[tuple[tuple[Fraction, ...], ...], ...], ...]:
    """Return the Kulkarni--Nomizu product of two symmetric 4-tensors."""

    data = [
        [[[Q(0) for _ in range(DIMENSION)] for _ in range(DIMENSION)] for _ in range(DIMENSION)]
        for _ in range(DIMENSION)
    ]
    for first, second, third, fourth in product(range(DIMENSION), repeat=4):
        data[first][second][third][fourth] = (
            left[first][third] * right[second][fourth]
            + left[second][fourth] * right[first][third]
            - left[first][fourth] * right[second][third]
            - left[second][third] * right[first][fourth]
        )
    return tuple(
        tuple(tuple(tuple(entry for entry in line) for line in plane) for plane in volume)
        for volume in data
    )


def _hessian_is_symmetric(matrix: Sequence[Sequence[Fraction]]) -> bool:
    return all(
        matrix[row][column] == matrix[column][row]
        for row, column in product(range(DIMENSION), repeat=2)
    )


def _algebraic_curvature_identities(
    tensor: Sequence[Sequence[Sequence[Sequence[Fraction]]]],
) -> dict[str, bool]:
    pair_antisymmetric = True
    pair_exchange = True
    first_bianchi = True
    for first, second, third, fourth in product(range(DIMENSION), repeat=4):
        value = tensor[first][second][third][fourth]
        if value != -tensor[second][first][third][fourth]:
            pair_antisymmetric = False
        if value != -tensor[first][second][fourth][third]:
            pair_antisymmetric = False
        if value != tensor[third][fourth][first][second]:
            pair_exchange = False
        if (
            value
            + tensor[first][third][fourth][second]
            + tensor[first][fourth][second][third]
            != 0
        ):
            first_bianchi = False
    return {
        "pair_antisymmetric": pair_antisymmetric,
        "pair_exchange": pair_exchange,
        "first_bianchi": first_bianchi,
    }


def _riemann_0101_rational() -> tuple[tuple[tuple[tuple[Fraction, ...], ...], ...], ...]:
    return _kulkarni_nomizu(
        _rational_symmetric_matrix((0, 0)),
        _rational_symmetric_matrix((1, 1)),
    )


def _mutate_tensor_break_bianchi(
    tensor: Sequence[Sequence[Sequence[Sequence[Fraction]]]],
) -> tuple[tuple[tuple[tuple[Fraction, ...], ...], ...], ...]:
    data = [
        [[list(line) for line in plane] for plane in volume]
        for volume in tensor
    ]
    data[0][1][2][3] = data[0][1][2][3] + 1
    return tuple(
        tuple(tuple(tuple(entry for entry in line) for line in plane) for plane in volume)
        for volume in data
    )


def _mutate_matrix_break_symmetry(
    matrix: Sequence[Sequence[Fraction]],
) -> tuple[tuple[Fraction, ...], ...]:
    data = [list(row) for row in matrix]
    data[0][1] = data[0][1] + 1
    return tuple(tuple(row) for row in data)


def _build_structural_generator_basis(
    *,
    omit_hessian_index: int | None = None,
    omit_curvature_index: int | None = None,
    mutate_curvature_bianchi: bool = False,
    mutate_hessian_symmetry: bool = False,
) -> dict[str, Any]:
    """Construct and independently certify the complete exact generator basis."""

    hessian_rationals = tuple(
        _rational_symmetric_matrix(component) for component in SYMMETRIC_INDEX_PAIRS
    )
    if mutate_hessian_symmetry:
        broken = list(hessian_rationals)
        broken[0] = _mutate_matrix_break_symmetry(broken[0])
        hessian_rationals = tuple(broken)
    hessian_labels = list(SYMMETRIC_INDEX_PAIRS)
    if omit_hessian_index is not None:
        if not 0 <= omit_hessian_index < len(hessian_rationals):
            raise ValueError("omit_hessian_index is out of range")
        hessian_rationals = tuple(
            matrix
            for index, matrix in enumerate(hessian_rationals)
            if index != omit_hessian_index
        )
        hessian_labels = [
            label
            for index, label in enumerate(hessian_labels)
            if index != omit_hessian_index
        ]
    hessian_symmetric = all(_hessian_is_symmetric(matrix) for matrix in hessian_rationals)
    hessian_rank = (
        rank(tuple(_flatten_matrix4(matrix) for matrix in hessian_rationals))
        if hessian_rationals
        else 0
    )

    candidates: list[tuple[tuple[tuple[int, int], tuple[int, int]], tuple]] = []
    for first_index, left_label in enumerate(SYMMETRIC_INDEX_PAIRS):
        left = _rational_symmetric_matrix(left_label)
        for right_label in SYMMETRIC_INDEX_PAIRS[first_index:]:
            right = _rational_symmetric_matrix(right_label)
            candidates.append(((left_label, right_label), _kulkarni_nomizu(left, right)))
    if len(candidates) != KULKARNI_NOMIZU_CANDIDATE_COUNT:
        raise SGBLConeStop(
            "broken_gauge",
            "Kulkarni-Nomizu candidate enumeration is incomplete",
            {"count": len(candidates)},
        )
    candidate_identities = [
        _algebraic_curvature_identities(tensor) for _, tensor in candidates
    ]
    candidates_all_algebraic = all(
        item["pair_antisymmetric"] and item["pair_exchange"] and item["first_bianchi"]
        for item in candidate_identities
    )
    candidate_rank = rank(tuple(_flatten_tensor4(tensor) for _, tensor in candidates))

    selected: list[tuple[tuple[tuple[int, int], tuple[int, int]], tuple]] = []
    selected_flat: list[tuple[Fraction, ...]] = []
    for label, tensor in candidates:
        trial = tuple(selected_flat + [_flatten_tensor4(tensor)])
        if rank(trial) == len(selected) + 1:
            selected.append((label, tensor))
            selected_flat.append(_flatten_tensor4(tensor))
        if len(selected) == ALGEBRAIC_CURVATURE_GENERATOR_COUNT:
            break
    curvature_labels = [label for label, _ in selected]
    curvature_rationals = tuple(tensor for _, tensor in selected)
    if mutate_curvature_bianchi and curvature_rationals:
        broken_curvature = list(curvature_rationals)
        broken_curvature[0] = _mutate_tensor_break_bianchi(broken_curvature[0])
        curvature_rationals = tuple(broken_curvature)
    if omit_curvature_index is not None:
        if not 0 <= omit_curvature_index < len(curvature_rationals):
            raise ValueError("omit_curvature_index is out of range")
        curvature_rationals = tuple(
            tensor
            for index, tensor in enumerate(curvature_rationals)
            if index != omit_curvature_index
        )
        curvature_labels = [
            label
            for index, label in enumerate(curvature_labels)
            if index != omit_curvature_index
        ]
    curvature_identity_records = [
        _algebraic_curvature_identities(tensor) for tensor in curvature_rationals
    ]
    curvature_pair_antisymmetric = all(
        item["pair_antisymmetric"] for item in curvature_identity_records
    )
    curvature_pair_exchange = all(
        item["pair_exchange"] for item in curvature_identity_records
    )
    curvature_first_bianchi = all(
        item["first_bianchi"] for item in curvature_identity_records
    )
    curvature_rank = (
        rank(tuple(_flatten_tensor4(tensor) for tensor in curvature_rationals))
        if curvature_rationals
        else 0
    )
    representative_atom = _riemann_0101_rational()
    representative_in_span = (
        rank(
            tuple(_flatten_tensor4(tensor) for tensor in curvature_rationals)
            + (_flatten_tensor4(representative_atom),)
        )
        == curvature_rank
        if curvature_rationals
        else False
    )
    complete = (
        len(hessian_rationals) == HESSIAN_GENERATOR_COUNT
        and hessian_rank == HESSIAN_GENERATOR_COUNT
        and hessian_symmetric
        and len(curvature_rationals) == ALGEBRAIC_CURVATURE_GENERATOR_COUNT
        and curvature_rank == ALGEBRAIC_CURVATURE_GENERATOR_COUNT
        and candidate_rank == ALGEBRAIC_CURVATURE_GENERATOR_COUNT
        and curvature_pair_antisymmetric
        and curvature_pair_exchange
        and curvature_first_bianchi
        and candidates_all_algebraic
    )
    if not complete:
        raise SGBLConeStop(
            "broken_gauge",
            "complete Hessian/curvature generator basis failed independent rank or identity checks",
            {
                "hessian_generator_count": len(hessian_rationals),
                "hessian_generator_rank": hessian_rank,
                "hessian_symmetric": hessian_symmetric,
                "algebraic_curvature_generator_count": len(curvature_rationals),
                "algebraic_curvature_generator_rank": curvature_rank,
                "kulkarni_nomizu_span_rank": candidate_rank,
                "pair_antisymmetric": curvature_pair_antisymmetric,
                "pair_exchange": curvature_pair_exchange,
                "first_bianchi": curvature_first_bianchi,
            },
        )
    return {
        "hessian_generator_count": len(hessian_rationals),
        "hessian_generator_rank": hessian_rank,
        "hessian_labels": tuple(hessian_labels),
        "hessian_generators": tuple(
            _interval_matrix_from_rational(matrix) for matrix in hessian_rationals
        ),
        "hessian_symmetric": hessian_symmetric,
        "algebraic_curvature_generator_count": len(curvature_rationals),
        "algebraic_curvature_generator_rank": curvature_rank,
        "algebraic_curvature_labels": tuple(curvature_labels),
        "algebraic_curvature_generators": tuple(
            _interval_tensor_from_rational(tensor) for tensor in curvature_rationals
        ),
        "kulkarni_nomizu_candidate_count": len(candidates),
        "kulkarni_nomizu_span_rank": candidate_rank,
        "pair_antisymmetric": curvature_pair_antisymmetric,
        "pair_exchange": curvature_pair_exchange,
        "first_bianchi": curvature_first_bianchi,
        "candidates_all_algebraic": candidates_all_algebraic,
        "representative_r0101_in_span": representative_in_span,
        "representative_generator_sufficient": False,
        "complete_exact_basis": True,
    }


@lru_cache(maxsize=1)
def _cached_structural_generator_basis() -> tuple:
    payload = _build_structural_generator_basis()
    return (
        payload["hessian_generators"],
        payload["algebraic_curvature_generators"],
        payload["hessian_labels"],
        payload["algebraic_curvature_labels"],
        MappingProxyType(
            {
                key: value
                for key, value in payload.items()
                if key
                not in {
                    "hessian_generators",
                    "algebraic_curvature_generators",
                    "hessian_labels",
                    "algebraic_curvature_labels",
                }
            }
        ),
    )


def sgbl_structural_generator_basis(
    *,
    omit_hessian_index: int | None = None,
    omit_curvature_index: int | None = None,
    mutate_curvature_bianchi: bool = False,
    mutate_hessian_symmetry: bool = False,
) -> dict[str, Any]:
    """Return the complete exact Hessian and algebraic-curvature generator basis."""

    mutated = (
        omit_hessian_index is not None
        or omit_curvature_index is not None
        or mutate_curvature_bianchi
        or mutate_hessian_symmetry
    )
    if mutated:
        return _build_structural_generator_basis(
            omit_hessian_index=omit_hessian_index,
            omit_curvature_index=omit_curvature_index,
            mutate_curvature_bianchi=mutate_curvature_bianchi,
            mutate_hessian_symmetry=mutate_hessian_symmetry,
        )
    hessian, curvature, hess_labels, curv_labels, meta = _cached_structural_generator_basis()
    payload = dict(meta)
    payload["hessian_generators"] = hessian
    payload["algebraic_curvature_generators"] = curvature
    payload["hessian_labels"] = hess_labels
    payload["algebraic_curvature_labels"] = curv_labels
    return payload


def _action_identities_for_box(
    box: SGBLPrincipalBackgroundBox,
) -> tuple[bool, bool, Interval, Interval]:
    counter = _EvaluationCounter(box.limits.max_residual_evaluations)
    action = _extract_coefficient_tensors(
        box,
        include_gauge=False,
        sqrt2=_sqrt2_enclosure(box.limits.sqrt2_bits),
        counter=counter,
        normalize=False,
    )
    return _identities_on_action(action, box.limits.sqrt2_bits)


@lru_cache(maxsize=4)
def _cached_esf_structural_identities(sqrt2_bits: int) -> tuple[bool, bool, Interval, Interval]:
    box = sgbl_esf_principal_box(
        planck_mass=1,
        alpha_gb=0,
        limits=SGBLConeLimits(sqrt2_bits=sqrt2_bits, max_residual_evaluations=256),
        allow_zero_coupling_control=True,
    )
    return _action_identities_for_box(box)


@lru_cache(maxsize=4)
def _cached_complete_generator_identities(sqrt2_bits: int) -> dict[str, Any]:
    basis = sgbl_structural_generator_basis()
    limits = SGBLConeLimits(sqrt2_bits=sqrt2_bits, max_residual_evaluations=1024)
    symmetric = True
    gauge = True
    hessian_checked = 0
    for hessian in basis["hessian_generators"]:
        box = SGBLPrincipalBackgroundBox(
            planck_mass=1,
            alpha_gb=1,
            hessian_phi_lower=hessian,
            limits=limits,
        )
        identities = _action_identities_for_box(box)
        symmetric = symmetric and identities[0]
        gauge = gauge and identities[1]
        hessian_checked += 1
    curvature_checked = 0
    for riemann in basis["algebraic_curvature_generators"]:
        box = SGBLPrincipalBackgroundBox(
            planck_mass=1,
            alpha_gb=1,
            riemann_lower=riemann,
            limits=limits,
        )
        identities = _action_identities_for_box(box)
        symmetric = symmetric and identities[0]
        gauge = gauge and identities[1]
        curvature_checked += 1
    if (
        hessian_checked != HESSIAN_GENERATOR_COUNT
        or curvature_checked != ALGEBRAIC_CURVATURE_GENERATOR_COUNT
    ):
        raise SGBLConeStop(
            "broken_gauge",
            "structural identities omitted a required generator",
            {
                "checked_hessian_generators": hessian_checked,
                "checked_curvature_generators": curvature_checked,
            },
        )
    return {
        "structural_hessian_symmetric": symmetric,
        "structural_pure_gauge_identity": gauge,
        "checked_hessian_generators": hessian_checked,
        "checked_curvature_generators": curvature_checked,
        "hessian_generator_rank": basis["hessian_generator_rank"],
        "algebraic_curvature_generator_rank": basis["algebraic_curvature_generator_rank"],
        "kulkarni_nomizu_candidate_count": basis["kulkarni_nomizu_candidate_count"],
        "kulkarni_nomizu_span_rank": basis["kulkarni_nomizu_span_rank"],
        "pair_antisymmetric": basis["pair_antisymmetric"],
        "pair_exchange": basis["pair_exchange"],
        "first_bianchi": basis["first_bianchi"],
        "representative_generator_sufficient": False,
        "complete_exact_basis": True,
    }


def sgbl_structural_identities(
    box: SGBLPrincipalBackgroundBox,
) -> dict[str, Any]:
    """Algebraic Hessian symmetry and pure-gauge identity on the complete basis.

    Interval wrapping of a candidate box is not this certificate.  One Hess(00)
    or R(0101) representative is not sufficient.
    """

    if type(box) is not SGBLPrincipalBackgroundBox:
        raise TypeError("box must be SGBLPrincipalBackgroundBox")
    esf = _cached_esf_structural_identities(box.limits.sqrt2_bits)
    generators = _cached_complete_generator_identities(box.limits.sqrt2_bits)
    symmetric = esf[0] and generators["structural_hessian_symmetric"]
    gauge = esf[1] and generators["structural_pure_gauge_identity"]
    return {
        "structural_hessian_symmetric": symmetric,
        "structural_pure_gauge_identity": gauge,
        "checked_esf_generator": True,
        "checked_hessian_generators": generators["checked_hessian_generators"],
        "checked_curvature_generators": generators["checked_curvature_generators"],
        "hessian_generator_count": HESSIAN_GENERATOR_COUNT,
        "hessian_generator_rank": generators["hessian_generator_rank"],
        "algebraic_curvature_generator_count": ALGEBRAIC_CURVATURE_GENERATOR_COUNT,
        "algebraic_curvature_generator_rank": generators["algebraic_curvature_generator_rank"],
        "kulkarni_nomizu_candidate_count": generators["kulkarni_nomizu_candidate_count"],
        "kulkarni_nomizu_span_rank": generators["kulkarni_nomizu_span_rank"],
        "pair_antisymmetric": generators["pair_antisymmetric"],
        "pair_exchange": generators["pair_exchange"],
        "first_bianchi": generators["first_bianchi"],
        "representative_generator_sufficient": False,
        "complete_exact_basis": True,
    }


def _kinetic_singular_value_lower(
    matrix: Sequence[Sequence[Interval]],
    inverse_payload: Mapping[str, Any],
    bits: int,
) -> Fraction:
    """Lower-bound ``sigma_min(A)`` from a Neumann inverse 2-norm upper bound."""

    inverse_norm = _spectral_norm_upper(inverse_payload["inverse_enclosure"], bits)
    if inverse_norm <= 0:
        raise SGBLConeStop("lost_kinetic", "kinetic inverse-norm upper bound is nonpositive")
    return Q(1) / inverse_norm


@dataclass(frozen=True, slots=True, kw_only=True)
class SGBLCoefficientEnclosure:
    """Outward interval enclosure of the complete principal coefficient tensors."""

    time_time: tuple[tuple[Interval, ...], ...]
    time_space: tuple[tuple[tuple[Interval, ...], ...], ...]
    space_space: tuple[tuple[tuple[tuple[Interval, ...], ...], ...], ...]
    ungauged_time_time: tuple[tuple[Interval, ...], ...]
    ungauged_time_space: tuple[tuple[tuple[Interval, ...], ...], ...]
    ungauged_space_space: tuple[tuple[tuple[tuple[Interval, ...], ...], ...], ...]
    sha256: str
    evaluations: int
    max_rational_bit_length_observed: int

    def __post_init__(self) -> None:
        object.__setattr__(self, "time_time", _freeze_matrix(self.time_time))
        object.__setattr__(
            self,
            "time_space",
            tuple(_freeze_matrix(matrix) for matrix in self.time_space),
        )
        object.__setattr__(
            self,
            "space_space",
            tuple(
                tuple(_freeze_matrix(self.space_space[first][second]) for second in range(3))
                for first in range(3)
            ),
        )
        expected = _matrix_sha256(
            (self.time_time,)
            + self.time_space
            + tuple(self.space_space[first][second] for first in range(3) for second in range(3))
        )
        if self.sha256 != expected:
            raise ValueError("coefficient enclosure hash does not match the tensors")


@dataclass(frozen=True, slots=True, kw_only=True)
class SGBLConeRecord:
    """Continuum cone/symmetrizer record. Aggregate branch health stays false."""

    box: SGBLPrincipalBackgroundBox
    branch: str
    action: Mapping[str, Fraction | str]
    coefficient_enclosure: SGBLCoefficientEnclosure
    esf_reference: Mapping[str, Any]
    kinetic_rho_infinity: Fraction | None
    kinetic_singular_value_lower: Fraction | None
    companion_deformation_upper: Fraction | None
    action_energy_deformation_upper: Fraction | None
    cluster_resolvent_margins: Mapping[str, Fraction] | None
    physical_projector_displacement_upper: Fraction | None
    physical_energy_coercivity_lower: Fraction | None
    gauge_identity_defect: Interval
    action_hessian_symmetry_defect: Interval
    structural_hessian_symmetric: bool
    structural_pure_gauge_identity: bool
    exact_eigenframe_rank: int | None
    continuum_direction_bound_independent_of_n: bool
    real_complete_basis: bool
    positive_symmetrizer: bool
    classification: str
    inconclusive_reason: str | None
    residual_evaluations: int
    max_rational_bit_length_observed: int
    numpy_esf_tensors_contained: bool
    angular_regression_directions: tuple[tuple[float, float, float], ...]

    def __post_init__(self) -> None:
        if type(self.box) is not SGBLPrincipalBackgroundBox:
            raise TypeError("box must be SGBLPrincipalBackgroundBox")
        if self.branch != "SGB-L":
            raise ValueError("cone record is owned by the linear branch")
        if self.classification not in {
            "continuum_cone_proved",
            "interval_inconclusive",
        }:
            raise ValueError("unknown cone classification")
        if self.classification == "continuum_cone_proved":
            if not (
                self.real_complete_basis
                and self.positive_symmetrizer
                and self.continuum_direction_bound_independent_of_n
                and self.inconclusive_reason is None
            ):
                raise ValueError("proved continuum cone requires retained strict margins")
            if self.kinetic_singular_value_lower is None or self.kinetic_singular_value_lower <= 0:
                raise ValueError("proved cone requires a strictly positive kinetic margin")
            if self.physical_energy_coercivity_lower is None or self.physical_energy_coercivity_lower <= 0:
                raise ValueError("proved cone requires a strictly positive energy margin")
            if self.cluster_resolvent_margins is None or any(
                margin <= 0 for margin in self.cluster_resolvent_margins.values()
            ):
                raise ValueError("proved cone requires strictly positive resolvent margins")
            if not self.structural_hessian_symmetric or not self.structural_pure_gauge_identity:
                raise ValueError("proved cone requires algebraic Hessian and pure-gauge identities")
            if self.exact_eigenframe_rank != FIRST_ORDER_FIELD_COUNT:
                raise ValueError("proved cone requires an exact rank-24 ESF eigenframe")
        else:
            if self.real_complete_basis or self.positive_symmetrizer:
                raise ValueError("inconclusive cone cannot retain a real basis or symmetrizer")
        object.__setattr__(self, "action", MappingProxyType(dict(self.action)))
        if self.cluster_resolvent_margins is not None:
            object.__setattr__(
                self,
                "cluster_resolvent_margins",
                MappingProxyType(dict(self.cluster_resolvent_margins)),
            )
        object.__setattr__(self, "esf_reference", MappingProxyType(dict(self.esf_reference)))

    @property
    def strongly_hyperbolic(self) -> bool:
        return self.classification == "continuum_cone_proved"

    @property
    def SGBL_branch_owned_and_healthy(self) -> bool:
        return False

    @property
    def quantitative_all_covector_weak_coupling_health_envelope_passed(self) -> bool:
        return False

    @property
    def cone_certificate(self) -> str:
        return self.classification

    @property
    def representative_generator_sufficient(self) -> bool:
        return False


def sgbl_enclose_principal_coefficients(
    box: SGBLPrincipalBackgroundBox,
) -> SGBLCoefficientEnclosure:
    """Derive outward enclosures of complete A, B_i, C_ij and the ungauged form."""

    if type(box) is not SGBLPrincipalBackgroundBox:
        raise TypeError("box must be SGBLPrincipalBackgroundBox")
    sqrt2 = _sqrt2_enclosure(box.limits.sqrt2_bits)
    counter = _EvaluationCounter(box.limits.max_residual_evaluations)
    complete = _extract_coefficient_tensors(
        box, include_gauge=True, sqrt2=sqrt2, counter=counter
    )
    action = _extract_coefficient_tensors(
        box, include_gauge=False, sqrt2=sqrt2, counter=counter
    )
    frozen_time = _freeze_matrix(complete["time_time"])
    frozen_space = tuple(_freeze_matrix(matrix) for matrix in complete["time_space"])
    frozen_ss = tuple(
        tuple(_freeze_matrix(complete["space_space"][first][second]) for second in range(3))
        for first in range(3)
    )
    digest = _matrix_sha256(
        (frozen_time,)
        + frozen_space
        + tuple(frozen_ss[first][second] for first in range(3) for second in range(3))
    )
    return SGBLCoefficientEnclosure(
        time_time=frozen_time,
        time_space=frozen_space,
        space_space=frozen_ss,
        ungauged_time_time=_freeze_matrix(action["time_time"]),
        ungauged_time_space=tuple(_freeze_matrix(matrix) for matrix in action["time_space"]),
        ungauged_space_space=tuple(
            tuple(_freeze_matrix(action["space_space"][first][second]) for second in range(3))
            for first in range(3)
        ),
        sha256=digest,
        evaluations=counter.count,
        max_rational_bit_length_observed=counter.bits,
    )


def sgbl_enclose_principal_cone(box: SGBLPrincipalBackgroundBox) -> SGBLConeRecord:
    """Prove or fail-closed the continuum all-direction SGB-L cone on one box."""

    enclosure = sgbl_enclose_principal_coefficients(box)
    bits = box.limits.sqrt2_bits
    reference_box = sgbl_esf_principal_box(
        planck_mass=box.planck_mass.midpoint(),
        alpha_gb=0,
        limits=box.limits,
        allow_zero_coupling_control=True,
    )
    reference_enclosure = sgbl_enclose_principal_coefficients(reference_box)
    esf_reference = _esf_exact_reference_certificate(
        planck_mass=box.planck_mass.midpoint(),
        limits=box.limits,
    )
    numpy_tensors = principal_coefficient_tensors(
        esf_background(float(box.planck_mass.midpoint()))
    )
    numpy_contained = _contains_float_matrix(
        reference_enclosure.time_time, numpy_tensors.time_time
    ) and all(
        _contains_float_matrix(
            reference_enclosure.time_space[index], numpy_tensors.time_space[index]
        )
        for index in range(3)
    )
    identity_counter = _EvaluationCounter(box.limits.max_residual_evaluations)
    unnormalized_complete = _extract_coefficient_tensors(
        box,
        include_gauge=True,
        sqrt2=_sqrt2_enclosure(bits),
        counter=identity_counter,
        normalize=False,
    )
    unnormalized_action = _extract_coefficient_tensors(
        box,
        include_gauge=False,
        sqrt2=_sqrt2_enclosure(bits),
        counter=identity_counter,
        normalize=False,
    )
    symmetry = _extracted_field_symmetry(unnormalized_action)
    gauge_defect = _pure_gauge_operator_sum(unnormalized_action, bits)
    structural = sgbl_structural_identities(box)
    if (
        not structural["structural_hessian_symmetric"]
        or not structural["structural_pure_gauge_identity"]
    ):
        raise SGBLConeStop(
            "broken_gauge",
            "algebraic Hessian symmetry or pure-gauge identity failed on formula generators",
            structural,
        )
    action_map = {
        "planck_mass": box.planck_mass.midpoint(),
        "alpha_gb": box.alpha_gb.midpoint(),
        "beta": Q(0),
        "eta": Q(0),
        "F_prime": Q(0),
        "chart": ORTHONORMAL_MINKOWSKI_CHART,
    }

    def _inconclusive(reason: str, **fields: Any) -> SGBLConeRecord:
        return SGBLConeRecord(
            box=box,
            branch="SGB-L",
            action=action_map,
            coefficient_enclosure=enclosure,
            esf_reference=esf_reference,
            kinetic_rho_infinity=fields.get("kinetic_rho_infinity"),
            kinetic_singular_value_lower=fields.get("kinetic_singular_value_lower"),
            companion_deformation_upper=fields.get("companion_deformation_upper"),
            action_energy_deformation_upper=fields.get("action_energy_deformation_upper"),
            cluster_resolvent_margins=fields.get("cluster_resolvent_margins"),
            physical_projector_displacement_upper=fields.get(
                "physical_projector_displacement_upper"
            ),
            physical_energy_coercivity_lower=fields.get("physical_energy_coercivity_lower"),
            gauge_identity_defect=gauge_defect,
            action_hessian_symmetry_defect=symmetry,
            structural_hessian_symmetric=structural["structural_hessian_symmetric"],
            structural_pure_gauge_identity=structural["structural_pure_gauge_identity"],
            exact_eigenframe_rank=esf_reference.get("eigenframe_rank"),
            continuum_direction_bound_independent_of_n=True,
            real_complete_basis=False,
            positive_symmetrizer=False,
            classification="interval_inconclusive",
            inconclusive_reason=reason,
            residual_evaluations=(
                enclosure.evaluations
                + reference_enclosure.evaluations
                + identity_counter.count
            ),
            max_rational_bit_length_observed=max(
                enclosure.max_rational_bit_length_observed,
                reference_enclosure.max_rational_bit_length_observed,
            ),
            numpy_esf_tensors_contained=numpy_contained,
            angular_regression_directions=REGRESSION_DIRECTIONS,
        )

    unnorm_inverse = _interval_inverse(unnormalized_complete["time_time"])
    if unnorm_inverse is None:
        return _inconclusive("neumann_rho_not_below_one")
    kinetic_lower = _kinetic_singular_value_lower(
        unnormalized_complete["time_time"], unnorm_inverse, bits
    )
    if kinetic_lower <= 0:
        raise SGBLConeStop("lost_kinetic", "kinetic singular-value lower bound is not positive")
    reference_unnormalized = _extract_coefficient_tensors(
        reference_box,
        include_gauge=True,
        sqrt2=_sqrt2_enclosure(bits),
        counter=identity_counter,
        normalize=False,
    )
    reference_unnormalized_action = _extract_coefficient_tensors(
        reference_box,
        include_gauge=False,
        sqrt2=_sqrt2_enclosure(bits),
        counter=identity_counter,
        normalize=False,
    )
    reference_unnorm_inverse = _interval_inverse(reference_unnormalized["time_time"])
    if reference_unnorm_inverse is None:
        raise SGBLConeStop("lost_kinetic", "ESF unnormalized kinetic Neumann certificate failed")
    time_space, space_space = _companion_blocks(
        unnormalized_complete, unnorm_inverse["inverse_enclosure"]
    )
    reference_time_space, reference_space_space = _companion_blocks(
        reference_unnormalized, reference_unnorm_inverse["inverse_enclosure"]
    )
    companion_delta = _companion_deformation(
        time_space,
        space_space,
        reference_time_space,
        reference_space_space,
        bits,
    )
    energy_delta = _action_energy_deformation(
        unnormalized_action, reference_unnormalized_action, bits
    )
    cluster_margins = {
        name: record["certified_lower_bound"] - companion_delta
        for name, record in esf_reference["clusters"].items()
    }
    if any(margin <= 0 for margin in cluster_margins.values()):
        return _inconclusive(
            "resolvent_margin_not_strictly_positive",
            kinetic_rho_infinity=unnorm_inverse["rho_infinity"],
            kinetic_singular_value_lower=kinetic_lower,
            companion_deformation_upper=companion_delta,
            action_energy_deformation_upper=energy_delta,
            cluster_resolvent_margins=cluster_margins,
        )
    physical_lower = min(
        esf_reference["clusters"][name]["certified_lower_bound"]
        for name in ("physical_minus", "physical_plus")
    )
    physical_radius = box.limits.physical_contour_radius
    resolvent_neumann = companion_delta / physical_lower
    if resolvent_neumann >= 1:
        return _inconclusive(
            "physical_resolvent_neumann_not_below_one",
            kinetic_rho_infinity=unnorm_inverse["rho_infinity"],
            kinetic_singular_value_lower=kinetic_lower,
            companion_deformation_upper=companion_delta,
            action_energy_deformation_upper=energy_delta,
            cluster_resolvent_margins=cluster_margins,
        )
    projector = (
        physical_radius
        * (Q(1) / physical_lower) ** 2
        * companion_delta
        / (1 - resolvent_neumann)
    )
    if projector >= 1:
        return _inconclusive(
            "physical_projector_displacement_not_below_one",
            kinetic_rho_infinity=unnorm_inverse["rho_infinity"],
            kinetic_singular_value_lower=kinetic_lower,
            companion_deformation_upper=companion_delta,
            action_energy_deformation_upper=energy_delta,
            cluster_resolvent_margins=cluster_margins,
            physical_projector_displacement_upper=projector,
        )
    reference_energy = esf_reference["physical_energy_coercivity_lower"]
    reference_norm = esf_reference["physical_energy_operator_norm_upper"]
    transported = (
        reference_energy * (1 - projector) ** 2
        - reference_norm * (2 * (1 + projector) * projector + projector ** 2)
        - energy_delta
    )
    if transported <= 0:
        return _inconclusive(
            "physical_energy_coercivity_not_strictly_positive",
            kinetic_rho_infinity=unnorm_inverse["rho_infinity"],
            kinetic_singular_value_lower=kinetic_lower,
            companion_deformation_upper=companion_delta,
            action_energy_deformation_upper=energy_delta,
            cluster_resolvent_margins=cluster_margins,
            physical_projector_displacement_upper=projector,
            physical_energy_coercivity_lower=transported,
        )
    return SGBLConeRecord(
        box=box,
        branch="SGB-L",
        action=action_map,
        coefficient_enclosure=enclosure,
        esf_reference=esf_reference,
        kinetic_rho_infinity=unnorm_inverse["rho_infinity"],
        kinetic_singular_value_lower=kinetic_lower,
        companion_deformation_upper=companion_delta,
        action_energy_deformation_upper=energy_delta,
        cluster_resolvent_margins=cluster_margins,
        physical_projector_displacement_upper=projector,
        physical_energy_coercivity_lower=transported,
        gauge_identity_defect=gauge_defect,
        action_hessian_symmetry_defect=symmetry,
        structural_hessian_symmetric=structural["structural_hessian_symmetric"],
        structural_pure_gauge_identity=structural["structural_pure_gauge_identity"],
        exact_eigenframe_rank=esf_reference["eigenframe_rank"],
        continuum_direction_bound_independent_of_n=True,
        real_complete_basis=True,
        positive_symmetrizer=True,
        classification="continuum_cone_proved",
        inconclusive_reason=None,
        residual_evaluations=(
            enclosure.evaluations
            + reference_enclosure.evaluations
            + identity_counter.count
        ),
        max_rational_bit_length_observed=max(
            enclosure.max_rational_bit_length_observed,
            reference_enclosure.max_rational_bit_length_observed,
            identity_counter.bits,
        ),
        numpy_esf_tensors_contained=numpy_contained,
        angular_regression_directions=REGRESSION_DIRECTIONS,
    )


def sgbl_esf_exact_reference(
    *,
    planck_mass: Fraction | int = 2,
    limits: SGBLConeLimits | None = None,
) -> dict[str, Any]:
    """Public exact unnormalized ESF radial eigenframe/resolvent/energy certificate."""

    return _esf_exact_reference_certificate(
        planck_mass=_fraction("planck_mass", planck_mass),
        limits=limits or SGBLConeLimits(),
    )


def sgbl_pure_gauge_identity_defect(
    box: SGBLPrincipalBackgroundBox,
    *,
    mutate_sign: int = 1,
) -> Interval:
    """Return the unnormalized cubic gauge-identity bound, with an optional mutation."""

    counter = _EvaluationCounter(box.limits.max_residual_evaluations)
    action = _extract_coefficient_tensors(
        box,
        include_gauge=False,
        sqrt2=_sqrt2_enclosure(box.limits.sqrt2_bits),
        counter=counter,
        normalize=False,
    )
    return _pure_gauge_operator_sum(
        action, box.limits.sqrt2_bits, mutate_sign=mutate_sign
    )


def sgbl_angular_regression_imaginary_parts(
    box: SGBLPrincipalBackgroundBox,
) -> dict[tuple[float, float, float], float]:
    """Weaker directional spectra. These never certify the continuum cone."""

    background = CovariantPrincipalBackground(
        effective_planck_squared=float(box.effective_planck_squared.midpoint()),
        effective_planck_prime=0.0,
        gb_coupling_prime=float(box.alpha_gb.midpoint()),
        riemann_lower=np.array(
            [
                [
                    [
                        [float(box.riemann_lower[a][b][c][d].midpoint()) for d in range(4)]
                        for c in range(4)
                    ]
                    for b in range(4)
                ]
                for a in range(4)
            ],
            dtype=np.float64,
        ),
        hessian_gb_lower=np.array(
            [
                [float(box.hessian_gb_lower[a][b].midpoint()) for b in range(4)]
                for a in range(4)
            ],
            dtype=np.float64,
        ),
    )
    tensors = principal_coefficient_tensors(
        background,
        tilde_normal_factor=float(box.tilde_normal_factor),
        hat_normal_factor=float(box.hat_normal_factor),
    )
    answers: dict[tuple[float, float, float], float] = {}
    for direction in REGRESSION_DIRECTIONS:
        spatial = np.asarray(direction, dtype=np.float64)
        a_matrix = tensors.time_time
        b_matrix = sum(
            tensors.time_space[index] * spatial[index] for index in range(3)
        )
        c_matrix = np.einsum("i,ijab,j->ab", spatial, tensors.space_space, spatial)
        companion = np.block(
            [
                [np.zeros((FIELD_COUNT, FIELD_COUNT)), np.eye(FIELD_COUNT)],
                [-np.linalg.solve(a_matrix, c_matrix), -np.linalg.solve(a_matrix, b_matrix)],
            ]
        )
        answers[direction] = float(np.max(np.abs(np.linalg.eigvals(companion).imag)))
    return answers


__all__ = [
    "ALGEBRAIC_CURVATURE_GENERATOR_COUNT",
    "ANGULAR_DIRECTION",
    "DEFAULT_ARC_COUNT",
    "FIELD_ORDER",
    "FORBIDDEN_HEALTH_IMPORTS",
    "HESSIAN_GENERATOR_COUNT",
    "KULKARNI_NOMIZU_CANDIDATE_COUNT",
    "MIXED_DIRECTION",
    "ORTHONORMAL_MINKOWSKI_CHART",
    "RADIAL_DIRECTION",
    "REGRESSION_DIRECTIONS",
    "SGBLCoefficientEnclosure",
    "SGBLConeLimits",
    "SGBLConeRecord",
    "SGBLConeStop",
    "SGBLPrincipalBackgroundBox",
    "SYMMETRIC_INDEX_PAIRS",
    "sgbl_angular_regression_imaginary_parts",
    "sgbl_enclose_principal_coefficients",
    "sgbl_enclose_principal_cone",
    "sgbl_esf_exact_reference",
    "sgbl_esf_principal_box",
    "sgbl_nonflat_source_principal_box",
    "sgbl_principal_box_from_state",
    "sgbl_pure_gauge_identity_defect",
    "sgbl_structural_generator_basis",
    "sgbl_structural_identities",
]
