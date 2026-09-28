"""Matched-family SGB-L principal-background feeder (prospective).

This owner consumes an authenticated completed initial-health ODE cell or
record together with the exact family specification and builds interval
lower/full two-jets without treating binary64 values as enclosures.  Compact
profiles are extended to the chi two-jet, including ``chi_pi_r = chi_tr``.
``lambda_r`` and ``k_r`` come from the affine H/M constraints.
``lambda_rr`` and ``k_rr`` are enclosed by exact interval differentiation of
that RHS, or the third-profile derivative owner is named and the slot is
left unfilled.  Unit-lapse, zero-shift polar-areal ADM two-jets receive
metric/scalar ``dtt`` only from a parametric source-admission certificate.
Direct interval 4D Riemann and Hess(phi) are pulled back by a checked
positive orthonormal frame with outward rational square-root factors, and
the image is a :class:`SGBLPrincipalBackgroundBox` for
``sgbl_enclose_principal_cone``.

The exact Minkowski buffer, completed support cells, and analytic exterior
are emitted only from authenticated inputs.  A nominal ``A_chi=3``
incomplete graph is one typed ``family_geometry_incomplete`` record, not a
fabricated exterior or health pass.  Aggregate
``SGBL_branch_owned_and_healthy`` stays false.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from fractions import Fraction
from numbers import Integral
from types import MappingProxyType
from typing import Any, Mapping, Sequence

from .exact_interval import Interval, interval
from .interval_tangent import IntervalFirstTangent, primal_and_tangent
from .sgb1_ctl1_admission import (
    SGBLAdmissionInconclusive,
    SGBLSourceAdmissionRecord,
    sgbl_parametric_source_admission,
)
from .sgb1_ctl1_cone import (
    ORTHONORMAL_MINKOWSKI_CHART,
    SGBLConeLimits,
    SGBLPrincipalBackgroundBox,
)
from .sgb1_ctl1_family import sgbl_affine_constraint_coefficients
from .sgb1_ctl1_initial_health import (
    SGBLContinuousCompactnessRecord,
    SGBLExactInitialSlice,
    SGBLInitialHealthStop,
    SGBLODECellEnclosure,
    sgbl_compact_bump_enclosure,
    sgbl_interval_sqrt,
)
from .sgb1_ctl1_interval_health import SGBLIntervalBox, SGBLIntervalLimits
from .sgb1_ctl1_source import (
    REFERENCE_RADIAL_MINIMUM,
    SGBLSourceInputs,
    SGBLSourceJacobianSolveStop,
    sgbl_source_solve,
    sgbl_source_state,
)
from .spherical_reduction import (
    Jet2,
    SphericalState,
    _state_from_adm_pg_jets,
    direct_4d_curvature,
)


Q = Fraction
DIMENSION = 4
ZERO = Interval.singleton(0)
ONE = Interval.singleton(1)
ETA_SIGNS = (-1, 1, 1, 1)
MINKOWSKI = tuple(
    tuple(Interval.singleton(ETA_SIGNS[index] if index == column else 0) for column in range(DIMENSION))
    for index in range(DIMENSION)
)
DECLARED_DYADIC_DENOMINATOR = 2 ** 64
BUMP_XXX_X_DERIVATIVE_BOUND = 39312
FAMILY_ACCELERATION_HALF_WIDTH = Q(1, 1 << 10)
DEFAULT_MAX_FAMILY_PRINCIPAL_CELLS = 16
FORBIDDEN_HEALTH_IMPORTS = (
    "WeakCouplingThresholds",
    "weak_coupling_health_certificate",
    "esf_reference_cluster_certificate",
    "canonical_health_monitor_values",
    "background_from_spherical_state",
)
FAMILY_PRINCIPAL_STOP_REASONS = frozenset(
    {
        "family_geometry_incomplete",
        "missing_third_profile_derivative_owner",
        "source_admission_incomplete",
        "cone_chart_error",
        "sign_chart_error",
        "resource_limit",
        "domain_error",
        "frame_not_orthonormal",
        "float_enclosure_refused",
    }
)
SLOT_OWNERSHIP = (
    {
        "slot": "r",
        "owner": "SGBLODECellEnclosure.radius or exact buffer/exterior radius",
        "formula": "authenticated radial interval or singleton",
    },
    {
        "slot": "lambda",
        "owner": "SGBLODECellEnclosure.lambda_box",
        "formula": "Picard graph of the affine constraint ODE",
    },
    {
        "slot": "k",
        "owner": "SGBLODECellEnclosure.k_box",
        "formula": "Picard graph; Minkowski buffer k=0",
    },
    {
        "slot": "phi, phi_r, phi_rr",
        "owner": "sgbl_compact_bump_enclosure",
        "formula": "A_phi * (B, B_r, B_rr)",
    },
    {
        "slot": "phi_rrr",
        "owner": "sgbl_compact_bump_third_derivative_enclosure",
        "formula": "A_phi * B_rrr; B_rrr = B (f' S + S') / w^3",
    },
    {
        "slot": "phi_pi, phi_pi_r, phi_pi_rr",
        "owner": "declared compact family",
        "formula": "0 on the unit-lapse zero-shift slice",
    },
    {
        "slot": "chi",
        "owner": "sgbl_complete_interval_family_fields",
        "formula": "A_chi * B / r",
    },
    {
        "slot": "chi_r",
        "owner": "sgbl_complete_interval_family_fields",
        "formula": "A_chi * (B_r/r - B/r^2)",
    },
    {
        "slot": "chi_rr",
        "owner": "sgbl_complete_interval_family_fields",
        "formula": "A_chi * (B_rr/r - 2 B_r/r^2 + 2 B/r^3)",
    },
    {
        "slot": "chi_pi",
        "owner": "sgbl_complete_interval_family_fields",
        "formula": "A_chi * B_r / r = chi_t",
    },
    {
        "slot": "chi_pi_r / chi_tr",
        "owner": "sgbl_complete_interval_family_fields",
        "formula": "A_chi * (B_rr/r - B_r/r^2)",
    },
    {
        "slot": "lambda_r",
        "owner": "affine H = H0 + H_L lambda_r",
        "formula": "lambda_r = -H0 / H_L",
    },
    {
        "slot": "k_r",
        "owner": "affine M = M0 + M_k k_r",
        "formula": "k_r = -M0 / M_k",
    },
    {
        "slot": "lambda_rr, k_rr",
        "owner": "interval differentiation of the affine RHS",
        "formula": "d/dr of (-H0/H_L, -M0/M_k) with IntervalFirstTangent; needs phi_rrr",
    },
    {
        "slot": "alpha",
        "owner": "unit-lapse polar-areal gauge",
        "formula": "alpha = 1, spatial derivatives 0",
    },
    {
        "slot": "shift, shift_t, shift_tr",
        "owner": "zero-shift plus C^t=C^r=0",
        "formula": "shift=0, shift_t=(2 L^3-2 L+L_r r)/(4 L^3 r)",
    },
    {
        "slot": "lambda_t, lambda_tr",
        "owner": "polar-areal K^r_r=-2k",
        "formula": "lambda_t=2 L k, lambda_tr=2(L_r k + L k_r)",
    },
    {
        "slot": "R, R_t, R_r, R_tr, R_rr",
        "owner": "polar-areal identity R=r",
        "formula": "R_t=-r k, R_r=1, R_tr=-k-r k_r, R_rr=0",
    },
    {
        "slot": "alpha_tt, shift_tt, lambda_tt, R_tt, phi_tt, chi_tt",
        "owner": "sgbl_parametric_source_admission",
        "formula": "Krawczyk image of R(a;z)=R0(z)+J(z)a; never a hidden zero",
    },
    {
        "slot": "Riemann_coord",
        "owner": "spherical_reduction.direct_4d_curvature",
        "formula": "Christoffel/Riemann from the complete two-jet",
    },
    {
        "slot": "Hess(phi)_coord",
        "owner": "covariant Hessian on the spherical chart",
        "formula": "partial_a partial_b phi - Gamma^c_ab partial_c phi",
    },
    {
        "slot": "orthonormal frame",
        "owner": "checked positive polar frame",
        "formula": "n=-lapse g^{-1} dt, e_r=(0, 1/sqrt(h_rr)), e_theta=e_phi=1/R",
    },
    {
        "slot": "Riemann_orth, Hess(phi)_orth",
        "owner": "frame pullback",
        "formula": "R_abcd = e^mu_a e^nu_b e^rho_c e^sigma_d R_munurhosigma",
    },
    {
        "slot": "SGBLPrincipalBackgroundBox",
        "owner": "sgb1_ctl1_cone",
        "formula": "F=Mpl^2, F'=0, f'=alpha_gb, Hess(f)=alpha_gb Hess(phi)",
    },
)


def _fraction(name: str, value: object) -> Fraction:
    if type(value) not in (int, Fraction):
        raise TypeError(f"{name} must be a Fraction or built-in int")
    return Fraction(value)


def _positive_int(name: str, value: object) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or int(value) <= 0:
        raise ValueError(f"{name} must be a positive integer")
    return int(value)


def _as_interval(value: object) -> Interval:
    if type(value) is Interval:
        return value
    if type(value) is IntervalFirstTangent:
        return value.primal
    if type(value) in (int, Fraction):
        return Interval.singleton(Fraction(value))
    if isinstance(value, float):
        raise SGBLFamilyPrincipalStop(
            "float_enclosure_refused",
            "family principal jets refuse binary64-as-enclosure inputs",
        )
    raise TypeError("value must be an Interval or exact rational")


def _floor_dyadic(value: Fraction, denominator: int) -> Fraction:
    return Fraction((value.numerator * denominator) // value.denominator, denominator)


def _ceil_dyadic(value: Fraction, denominator: int) -> Fraction:
    return -_floor_dyadic(-value, denominator)


def _outward(value: Interval, denominator: int = DECLARED_DYADIC_DENOMINATOR) -> Interval:
    """Valid outward rounding onto the declared dyadic grid."""

    return Interval(
        _floor_dyadic(value.lower, denominator),
        _ceil_dyadic(value.upper, denominator),
    )


def _witness(box: Interval) -> Fraction:
    """Return an exact rational inside ``box``.  Zero is preferred when present."""

    if box.is_singleton():
        return box.lower
    if box.contains_zero():
        return Q(0)
    return box.midpoint()


def _interval_matrix4(
    value: Sequence[Sequence[object]],
) -> tuple[tuple[Interval, ...], ...]:
    matrix = tuple(tuple(_as_interval(entry) for entry in row) for row in value)
    if len(matrix) != DIMENSION or any(len(row) != DIMENSION for row in matrix):
        raise TypeError("expected a 4x4 interval matrix")
    hulled = [[matrix[row][column] for column in range(DIMENSION)] for row in range(DIMENSION)]
    for first in range(DIMENSION):
        for second in range(first + 1, DIMENSION):
            lower = min(hulled[first][second].lower, hulled[second][first].lower)
            upper = max(hulled[first][second].upper, hulled[second][first].upper)
            merged = Interval(lower, upper)
            hulled[first][second] = merged
            hulled[second][first] = merged
    return tuple(tuple(row) for row in hulled)


def _interval_tensor4(
    value: Sequence[Sequence[Sequence[Sequence[object]]]],
) -> tuple[tuple[tuple[tuple[Interval, ...], ...], ...], ...]:
    tensor = tuple(
        tuple(
            tuple(tuple(_as_interval(entry) for entry in line) for line in plane)
            for plane in row
        )
        for row in value
    )
    if len(tensor) != DIMENSION:
        raise TypeError("expected a 4x4x4x4 interval tensor")
    return tensor


class SGBLFamilyPrincipalStop(ArithmeticError):
    """Typed feeder stop.  Incomplete family graphs are records, not this type."""

    def __init__(
        self,
        reason: str,
        message: str,
        payload: Mapping[str, Any] | None = None,
    ) -> None:
        if reason not in FAMILY_PRINCIPAL_STOP_REASONS:
            raise ValueError("unknown SGB-L family-principal stop reason")
        self.reason = reason
        self.payload = dict(payload or {})
        super().__init__(message)


@dataclass(frozen=True, slots=True, kw_only=True)
class SGBLFamilyPrincipalMissingOwner:
    """Immutable record of a derivative or source-box theorem that is not filled."""

    slot: str
    owner: str
    reason: str
    theorem: str

    def __post_init__(self) -> None:
        if not self.slot or not self.owner or not self.reason or not self.theorem:
            raise ValueError("missing-owner records must name slot, owner, reason and theorem")
        if self.reason in {"substituted_zero", "binary64_placeholder"}:
            raise ValueError("missing owners must not be filled by zero or binary64")


def sgbl_family_principal_slot_table() -> tuple[Mapping[str, str], ...]:
    """Return the complete slot-ownership table for this feeder."""

    return tuple(MappingProxyType(dict(entry)) for entry in SLOT_OWNERSHIP)


def sgbl_compact_bump_third_derivative_enclosure(
    radius: Interval,
    *,
    center: Fraction,
    half_width: Fraction,
) -> Interval:
    """Interval enclosure of the compact-bump third radial derivative.

    The scaled third derivative of ``B(x)=exp(1-1/(1-x^2))`` obeys
    ``|B'''(x)| < 39312`` on the real line, from ``e<3`` and
    ``e^{-u} u^n <= n!``.  Cells wholly outside the support are exactly
    zero.  Interior cells use the natural interval extension intersected
    with that global bound.  Boundary-touching cells use the global bound.
    """

    if half_width <= 0:
        raise SGBLFamilyPrincipalStop("domain_error", "bump half-width must be positive")
    global_third = interval(
        -BUMP_XXX_X_DERIVATIVE_BOUND, BUMP_XXX_X_DERIVATIVE_BOUND
    ) / (half_width ** 3)
    scaled = (radius - center) / half_width
    if scaled.upper <= -1 or scaled.lower >= 1:
        return ZERO
    if scaled.lower <= -1 or scaled.upper >= 1:
        return global_third
    square = scaled ** 2
    one_minus_square = ONE - square
    if not one_minus_square.strictly_positive():
        return global_third
    value, first, second = sgbl_compact_bump_enclosure(
        radius, center=center, half_width=half_width
    )
    first_exponent = (interval(-2, -2) * scaled) / (one_minus_square * one_minus_square)
    second_exponent = interval(-2, -2) / (one_minus_square * one_minus_square) + (
        interval(-8, -8) * square
    ) / (one_minus_square * one_minus_square * one_minus_square)
    third_exponent = (interval(-24, -24) * scaled) / (
        one_minus_square * one_minus_square * one_minus_square
    ) + (interval(-48, -48) * (scaled ** 3)) / (
        one_minus_square * one_minus_square * one_minus_square * one_minus_square
    )
    curvature_s = second_exponent + first_exponent ** 2
    curvature_s_prime = third_exponent + interval(2, 2) * first_exponent * second_exponent
    third = (
        value * (first_exponent * curvature_s + curvature_s_prime) / (half_width ** 3)
    )
    return Interval(
        max(third.lower, global_third.lower),
        min(third.upper, global_third.upper),
    )


def sgbl_complete_interval_family_fields(
    radius: Interval,
    spec: SGBLExactInitialSlice,
) -> dict[str, Interval]:
    """Enclose the lambda-free compact profiles including chi two-jet slots."""

    if type(spec) is not SGBLExactInitialSlice:
        raise TypeError("spec must be SGBLExactInitialSlice")
    if not radius.strictly_positive():
        raise SGBLFamilyPrincipalStop(
            "domain_error",
            "compact-support cells must lie at positive radius",
        )
    bump, bump_r, bump_rr = sgbl_compact_bump_enclosure(
        radius, center=spec.center, half_width=spec.half_width
    )
    bump_rrr = sgbl_compact_bump_third_derivative_enclosure(
        radius, center=spec.center, half_width=spec.half_width
    )
    inv_r = radius.reciprocal()
    inv_r2 = inv_r * inv_r
    inv_r3 = inv_r2 * inv_r
    chi_amp = Interval.singleton(spec.chi_amplitude)
    phi_amp = Interval.singleton(spec.phi_amplitude)
    chi_pi_r = _outward(chi_amp * (bump_rr * inv_r - bump_r * inv_r2))
    return {
        "phi": _outward(phi_amp * bump),
        "phi_r": _outward(phi_amp * bump_r),
        "phi_rr": _outward(phi_amp * bump_rr),
        "phi_rrr": _outward(phi_amp * bump_rrr),
        "phi_pi": ZERO,
        "phi_pi_r": ZERO,
        "phi_pi_rr": ZERO,
        "chi": _outward(chi_amp * bump * inv_r),
        "chi_r": _outward(chi_amp * (bump_r * inv_r - bump * inv_r2)),
        "chi_rr": _outward(
            chi_amp * (bump_rr * inv_r - interval(2, 2) * bump_r * inv_r2 + interval(2, 2) * bump * inv_r3)
        ),
        "chi_pi": _outward(chi_amp * bump_r * inv_r),
        "chi_pi_r": chi_pi_r,
        "chi_tr": chi_pi_r,
    }


def _affine_from_fields(
    radius: object,
    radial_metric: object,
    angular_extrinsic_curvature: object,
    fields: Mapping[str, object],
    spec: SGBLExactInitialSlice,
) -> dict[str, object]:
    """Generic H/M coefficients.  Accepts Interval or IntervalFirstTangent."""

    r = radius
    lam = radial_metric
    k = angular_extrinsic_curvature
    mpl2 = spec.planck_mass_squared
    alpha = spec.alpha_gb
    two = 2
    eight = 8
    sixteen = 16
    a0 = -2 * (k ** 2)
    a_l = 1 / (lam ** 3 * r)
    b = (k ** 2) + (1 - 1 / (lam ** 2)) / (r ** 2)
    c0 = 3 * k / (lam * r)
    c_k = 1 / lam
    x0 = fields["phi_rr"] / (lam ** 2) - two * k * fields["phi_pi"]
    x_l = -fields["phi_r"] / (lam ** 3)
    y = k * fields["phi_pi"] + fields["phi_r"] / (lam ** 2 * r)
    z = (fields["phi_pi_r"] - two * k * fields["phi_r"]) / lam
    rho = (
        (fields["phi_pi"] ** 2 + fields["chi_pi"] ** 2) / two
        + (fields["phi_r"] ** 2 + fields["chi_r"] ** 2) / (two * (lam ** 2))
        + spec.scalar_mass ** 2 * (fields["phi"] ** 2) / two
        + spec.quartic_coupling * (fields["phi"] ** 4) / 4
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


def _primal_interval(value: object) -> Interval:
    return primal_and_tangent(value)[0] if type(value) is IntervalFirstTangent else _as_interval(value)


def _solve_affine(coefficients: Mapping[str, object]) -> tuple[object, object]:
    h_l = coefficients["hamiltonian_lambda_r"]
    m_k = coefficients["momentum_k_r"]
    if _primal_interval(h_l).contains_zero() or _primal_interval(m_k).contains_zero():
        raise SGBLFamilyPrincipalStop(
            "domain_error",
            "affine constraint Jacobian diagonal contains zero",
            {"obstruction": "jacobian_diagonal_contains_zero"},
        )
    return (
        -coefficients["hamiltonian_constant"] / h_l,
        -coefficients["momentum_constant"] / m_k,
    )


def sgbl_family_affine_jets(
    radius: Interval,
    radial_metric: Interval,
    angular_extrinsic_curvature: Interval,
    spec: SGBLExactInitialSlice,
) -> dict[str, Interval]:
    """Enclose ``lambda_r, k_r`` and their radial derivatives on one cell box."""

    fields = sgbl_complete_interval_family_fields(radius, spec)
    if "phi_rrr" not in fields:
        raise SGBLFamilyPrincipalStop(
            "missing_third_profile_derivative_owner",
            "lambda_rr requires the compact-bump third radial derivative",
            {
                "missing_owner": "sgbl_compact_bump_third_derivative_enclosure",
                "slot": "phi_rrr",
            },
        )
    coefficients = _affine_from_fields(
        radius, radial_metric, angular_extrinsic_curvature, fields, spec
    )
    lambda_r, k_r = _solve_affine(coefficients)
    lambda_r_i = _outward(_primal_interval(lambda_r))
    k_r_i = _outward(_primal_interval(k_r))
    seeded_fields = {
        "phi": IntervalFirstTangent.seed(fields["phi"], fields["phi_r"]),
        "phi_r": IntervalFirstTangent.seed(fields["phi_r"], fields["phi_rr"]),
        "phi_rr": IntervalFirstTangent.seed(fields["phi_rr"], fields["phi_rrr"]),
        "phi_pi": IntervalFirstTangent.constant(fields["phi_pi"]),
        "phi_pi_r": IntervalFirstTangent.constant(fields["phi_pi_r"]),
        "chi_r": IntervalFirstTangent.seed(fields["chi_r"], fields["chi_rr"]),
        "chi_pi": IntervalFirstTangent.seed(fields["chi_pi"], fields["chi_pi_r"]),
    }
    seeded_coefficients = _affine_from_fields(
        IntervalFirstTangent.seed(radius, ONE),
        IntervalFirstTangent.seed(radial_metric, lambda_r_i),
        IntervalFirstTangent.seed(angular_extrinsic_curvature, k_r_i),
        seeded_fields,
        spec,
    )
    lambda_r_jet, k_r_jet = _solve_affine(seeded_coefficients)
    lambda_rr = _outward(primal_and_tangent(lambda_r_jet)[1])
    k_rr = _outward(primal_and_tangent(k_r_jet)[1])
    payload = dict(fields)
    payload.update(
        {
            "lambda_r": lambda_r_i,
            "k_r": k_r_i,
            "lambda_rr": lambda_rr,
            "k_rr": k_rr,
            "hamiltonian_lambda_r": _outward(_primal_interval(coefficients["hamiltonian_lambda_r"])),
            "momentum_k_r": _outward(_primal_interval(coefficients["momentum_k_r"])),
        }
    )
    return payload


def _shift_time_derivative(radius: object, radial_metric: object, lambda_r: object) -> object:
    """``shift_t = (2 L^3 - 2 L + L_r r) / (4 L^3 r)`` from C^t=C^r=0."""

    lam = radial_metric
    return (2 * lam ** 3 - 2 * lam + lambda_r * radius) / (4 * lam ** 3 * radius)


def sgbl_family_adm_lower_jets(
    *,
    radius: Interval,
    radial_metric: Interval,
    angular_extrinsic_curvature: Interval,
    spec: SGBLExactInitialSlice,
) -> dict[str, Jet2]:
    """Unit-lapse, zero-shift polar-areal lower two-jets (``dtt=0``)."""

    affine = sgbl_family_affine_jets(
        radius, radial_metric, angular_extrinsic_curvature, spec
    )
    lam = radial_metric
    k = angular_extrinsic_curvature
    lambda_r = affine["lambda_r"]
    k_r = affine["k_r"]
    lambda_rr = affine["lambda_rr"]
    shift_t = _shift_time_derivative(radius, lam, lambda_r)
    shift_tr = primal_and_tangent(
        _shift_time_derivative(
            IntervalFirstTangent.seed(radius, ONE),
            IntervalFirstTangent.seed(lam, lambda_r),
            IntervalFirstTangent.seed(lambda_r, lambda_rr),
        )
    )[1]
    return {
        "alpha": Jet2(ONE),
        "shift": Jet2(ZERO, dt=_outward(_as_interval(shift_t)), dtr=_outward(_as_interval(shift_tr))),
        "lambda": Jet2(
            lam,
            dt=_outward(interval(2, 2) * lam * k),
            dr=lambda_r,
            dtr=_outward(interval(2, 2) * (lambda_r * k + lam * k_r)),
            drr=affine["lambda_rr"],
        ),
        "areal_radius": Jet2(
            radius,
            dt=_outward(-(radius * k)),
            dr=ONE,
            dtr=_outward(-(k + radius * k_r)),
        ),
        "phi": Jet2(
            affine["phi"],
            dt=affine["phi_pi"],
            dr=affine["phi_r"],
            dtr=affine["phi_pi_r"],
            drr=affine["phi_rr"],
        ),
        "chi": Jet2(
            affine["chi"],
            dt=affine["chi_pi"],
            dr=affine["chi_r"],
            dtr=affine["chi_pi_r"],
            drr=affine["chi_rr"],
        ),
    }


def _action_mapping(spec: SGBLExactInitialSlice) -> dict[str, Fraction]:
    return {
        "planck_mass": spec.planck_mass,
        "scalar_mass": spec.scalar_mass,
        "quartic_coupling": spec.quartic_coupling,
        "alpha_gb": spec.alpha_gb,
        "beta": Q(0),
        "eta": Q(0),
    }


def _source_inputs_from_witness(
    *,
    radius: Fraction,
    radial_metric: Fraction,
    angular_extrinsic_curvature: Fraction,
    spec: SGBLExactInitialSlice,
) -> tuple[SGBLSourceInputs, dict[str, Interval]]:
    """Declared exact center inside the cell product box."""

    radius_i = Interval.singleton(radius)
    lam_i = Interval.singleton(radial_metric)
    k_i = Interval.singleton(angular_extrinsic_curvature)
    affine = sgbl_family_affine_jets(radius_i, lam_i, k_i, spec)
    fields = {
        name: _witness(affine[name])
        for name in (
            "phi",
            "phi_r",
            "phi_rr",
            "phi_pi",
            "phi_pi_r",
            "chi",
            "chi_r",
            "chi_rr",
            "chi_pi",
            "chi_pi_r",
        )
    }
    exact = sgbl_affine_constraint_coefficients(
        radius=radius,
        radial_metric=radial_metric,
        angular_extrinsic_curvature=angular_extrinsic_curvature,
        phi=fields["phi"],
        phi_r=fields["phi_r"],
        phi_rr=fields["phi_rr"],
        phi_pi=fields["phi_pi"],
        phi_pi_r=fields["phi_pi_r"],
        chi_r=fields["chi_r"],
        chi_pi=fields["chi_pi"],
        planck_mass=spec.planck_mass,
        scalar_mass=spec.scalar_mass,
        quartic_coupling=spec.quartic_coupling,
        alpha_gb=spec.alpha_gb,
    )
    lambda_r, k_r = exact.derivatives()
    lambda_rr = _witness(affine["lambda_rr"])
    shift_t = _shift_time_derivative(radius, radial_metric, lambda_r)
    shift_tr = _witness(
        _outward(
            _as_interval(
                primal_and_tangent(
                    _shift_time_derivative(
                        IntervalFirstTangent.seed(radius_i, ONE),
                        IntervalFirstTangent.seed(lam_i, Interval.singleton(lambda_r)),
                        IntervalFirstTangent.seed(
                            Interval.singleton(lambda_r), affine["lambda_rr"]
                        ),
                    )
                )[1]
            )
        )
    )
    point = SGBLSourceInputs(
        coordinate_radius=radius,
        alpha=(1, 0, 0, 0, 0),
        shift=(0, shift_t, 0, shift_tr, 0),
        radial_metric=(
            radial_metric,
            2 * radial_metric * angular_extrinsic_curvature,
            lambda_r,
            2 * (lambda_r * angular_extrinsic_curvature + radial_metric * k_r),
            lambda_rr,
        ),
        areal_radius=(
            radius,
            -radius * angular_extrinsic_curvature,
            1,
            -angular_extrinsic_curvature - radius * k_r,
            0,
        ),
        phi=(fields["phi"], fields["phi_pi"], fields["phi_r"], fields["phi_pi_r"], fields["phi_rr"]),
        chi=(fields["chi"], fields["chi_pi"], fields["chi_r"], fields["chi_pi_r"], fields["chi_rr"]),
        planck_mass=spec.planck_mass,
        scalar_mass=spec.scalar_mass,
        quartic_coupling=spec.quartic_coupling,
        alpha_gb=spec.alpha_gb,
    )
    return point, affine


def sgbl_family_cell_source_admission(
    point: SGBLSourceInputs,
    *,
    parameter_half_width: Fraction | int = 0,
    acceleration_half_width: Fraction | int = FAMILY_ACCELERATION_HALF_WIDTH,
    limits: SGBLIntervalLimits | None = None,
) -> SGBLSourceAdmissionRecord:
    """Parametric uniqueness certificate at a declared family source point."""

    box = SGBLIntervalBox(
        center=point,
        parameter_half_width=parameter_half_width,
        acceleration_half_width=acceleration_half_width,
        limits=limits or SGBLIntervalLimits(),
    )
    return sgbl_parametric_source_admission(box)


def _full_two_jet_state(
    point: SGBLSourceInputs,
    accelerations: Sequence[object],
    spec: SGBLExactInitialSlice,
) -> SphericalState:
    if len(accelerations) != 6:
        raise ValueError("ADM accelerations must have six slots")
    jets = {
        "alpha": replace(point.alpha, dtt=_as_interval(accelerations[0])),
        "shift": replace(point.shift, dtt=_as_interval(accelerations[1])),
        "lambda": replace(point.radial_metric, dtt=_as_interval(accelerations[2])),
        "areal_radius": replace(point.areal_radius, dtt=_as_interval(accelerations[3])),
        "phi": replace(point.phi, dtt=_as_interval(accelerations[4])),
        "chi": replace(point.chi, dtt=_as_interval(accelerations[5])),
    }
    return _state_from_adm_pg_jets(
        {"model_id": "SGB-L", "action_parameters": _action_mapping(spec)},
        jets,
    )


def sgbl_interval_scalar_hessian(
    state: SphericalState,
    connection: Sequence[Sequence[Sequence[object]]],
) -> tuple[tuple[Interval, ...], ...]:
    """Coordinate Hessian ``partial_a partial_b phi - Gamma^c_ab partial_c phi``."""

    gradient = (
        _as_interval(state.phi.dt),
        _as_interval(state.phi.dr),
        ZERO,
        ZERO,
    )
    second = [[ZERO for _ in range(DIMENSION)] for _ in range(DIMENSION)]
    second[0][0] = _as_interval(state.phi.dtt)
    second[0][1] = second[1][0] = _as_interval(state.phi.dtr)
    second[1][1] = _as_interval(state.phi.drr)
    hessian = [[ZERO for _ in range(DIMENSION)] for _ in range(DIMENSION)]
    for first in range(DIMENSION):
        for second_index in range(DIMENSION):
            total = second[first][second_index]
            for covector in range(DIMENSION):
                total = total - _as_interval(connection[covector][first][second_index]) * gradient[
                    covector
                ]
            hessian[first][second_index] = total
    return _interval_matrix4(hessian)


def sgbl_interval_orthonormal_frame(
    state: SphericalState,
) -> tuple[tuple[Interval, ...], ...]:
    """Positive orthonormal frame with outward rational square-root factors."""

    h_tt = _as_interval(state.h_tt.value)
    h_tr = _as_interval(state.h_tr.value)
    h_rr = _as_interval(state.h_rr.value)
    radius = _as_interval(state.areal_radius.value)
    if not h_rr.strictly_positive() or not radius.strictly_positive():
        raise SGBLFamilyPrincipalStop(
            "sign_chart_error",
            "orthonormal frame requires strictly positive h_rr and areal radius",
        )
    determinant = h_tt * h_rr - h_tr * h_tr
    if not determinant.strictly_negative():
        raise SGBLFamilyPrincipalStop(
            "sign_chart_error",
            "two-dimensional base must be strictly Lorentzian",
        )
    inverse_00 = h_rr / determinant
    inverse_01 = -h_tr / determinant
    if not inverse_00.strictly_negative():
        raise SGBLFamilyPrincipalStop(
            "sign_chart_error",
            "coordinate-time covector must be strictly timelike",
        )
    try:
        lapse = ONE / sgbl_interval_sqrt(-inverse_00)
        radial = ONE / sgbl_interval_sqrt(h_rr)
    except (SGBLInitialHealthStop, ZeroDivisionError, ValueError) as exc:
        raise SGBLFamilyPrincipalStop(
            "sign_chart_error",
            "positive square-root branch of the orthonormal frame failed",
        ) from exc
    angular = radius.reciprocal()
    frame = [[ZERO for _ in range(DIMENSION)] for _ in range(DIMENSION)]
    frame[0][0] = _outward(-lapse * inverse_00)
    frame[1][0] = _outward(-lapse * inverse_01)
    frame[1][1] = _outward(radial)
    frame[2][2] = _outward(angular)
    frame[3][3] = _outward(angular)
    metric = [[ZERO for _ in range(DIMENSION)] for _ in range(DIMENSION)]
    metric[0][0] = h_tt
    metric[0][1] = metric[1][0] = h_tr
    metric[1][1] = h_rr
    metric[2][2] = radius * radius
    metric[3][3] = radius * radius
    for first in range(DIMENSION):
        for second_index in range(DIMENSION):
            pulled = ZERO
            for mu in range(DIMENSION):
                for nu in range(DIMENSION):
                    pulled = pulled + frame[mu][first] * metric[mu][nu] * frame[nu][second_index]
            expected = MINKOWSKI[first][second_index]
            if pulled.lower > expected.lower or pulled.upper < expected.upper:
                raise SGBLFamilyPrincipalStop(
                    "frame_not_orthonormal",
                    "pulled spherical frame does not contain the Minkowski metric",
                    {
                        "indices": (first, second_index),
                        "pulled": (pulled.lower, pulled.upper),
                    },
                )
    return tuple(tuple(row) for row in frame)


def sgbl_checked_orthonormal_pullback(
    riemann_coordinate: Sequence[Sequence[Sequence[Sequence[object]]]],
    hessian_coordinate: Sequence[Sequence[object]],
    frame: Sequence[Sequence[Interval]],
) -> tuple[
    tuple[tuple[tuple[tuple[Interval, ...], ...], ...], ...],
    tuple[tuple[Interval, ...], ...],
]:
    """Pull Riemann and Hess(phi) to the orthonormal Minkowski frame."""

    riemann = [
        [[[ZERO for _ in range(DIMENSION)] for _ in range(DIMENSION)] for _ in range(DIMENSION)]
        for _ in range(DIMENSION)
    ]
    for a_index in range(DIMENSION):
        for b_index in range(DIMENSION):
            for c_index in range(DIMENSION):
                for d_index in range(DIMENSION):
                    total = ZERO
                    for mu in range(DIMENSION):
                        if frame[mu][a_index] == ZERO:
                            continue
                        for nu in range(DIMENSION):
                            if frame[nu][b_index] == ZERO:
                                continue
                            for rho in range(DIMENSION):
                                if frame[rho][c_index] == ZERO:
                                    continue
                                for sigma in range(DIMENSION):
                                    if frame[sigma][d_index] == ZERO:
                                        continue
                                    total = total + (
                                        frame[mu][a_index]
                                        * frame[nu][b_index]
                                        * frame[rho][c_index]
                                        * frame[sigma][d_index]
                                        * _as_interval(
                                            riemann_coordinate[mu][nu][rho][sigma]
                                        )
                                    )
                    riemann[a_index][b_index][c_index][d_index] = _outward(total)
    hessian = [[ZERO for _ in range(DIMENSION)] for _ in range(DIMENSION)]
    for a_index in range(DIMENSION):
        for b_index in range(DIMENSION):
            total = ZERO
            for mu in range(DIMENSION):
                for nu in range(DIMENSION):
                    total = total + frame[mu][a_index] * _as_interval(
                        hessian_coordinate[mu][nu]
                    ) * frame[nu][b_index]
            hessian[a_index][b_index] = _outward(total)
    return _interval_tensor4(riemann), _interval_matrix4(hessian)


def sgbl_principal_box_from_spherical_state(
    state: SphericalState,
    *,
    limits: SGBLConeLimits | None = None,
) -> SGBLPrincipalBackgroundBox:
    """Build an orthonormal interval box from a complete spherical two-jet."""

    if not isinstance(state, SphericalState):
        raise TypeError("state must be a SphericalState")
    if state.branch != "SGB-L":
        raise ValueError("family principal boxes require the linear branch")
    riemann_coordinate, _inverse, connection = direct_4d_curvature(state)
    hessian_coordinate = sgbl_interval_scalar_hessian(state, connection)
    frame = sgbl_interval_orthonormal_frame(state)
    riemann, hessian = sgbl_checked_orthonormal_pullback(
        riemann_coordinate, hessian_coordinate, frame
    )
    mass = _as_interval(state.planck_mass)
    alpha = _as_interval(state.alpha)
    if not mass.is_singleton() or not alpha.is_singleton():
        raise SGBLFamilyPrincipalStop(
            "domain_error",
            "family principal boxes keep exact singleton action parameters",
        )
    return SGBLPrincipalBackgroundBox(
        planck_mass=mass.lower,
        alpha_gb=alpha.lower,
        riemann_lower=riemann,
        hessian_phi_lower=hessian,
        limits=limits or SGBLConeLimits(),
        chart=ORTHONORMAL_MINKOWSKI_CHART,
    )


def sgbl_minkowski_buffer_principal_box(
    spec: SGBLExactInitialSlice | None = None,
    *,
    radius: Fraction | int | None = None,
    limits: SGBLConeLimits | None = None,
) -> SGBLPrincipalBackgroundBox:
    """Exact curvature-free box on the authenticated Minkowski buffer."""

    spec = spec or SGBLExactInitialSlice()
    if type(spec) is not SGBLExactInitialSlice:
        raise TypeError("spec must be SGBLExactInitialSlice")
    buffer_radius = _fraction("radius", spec.support_minimum if radius is None else radius)
    if buffer_radius <= REFERENCE_RADIAL_MINIMUM:
        raise SGBLFamilyPrincipalStop(
            "domain_error",
            "buffer radius must lie strictly above the frozen MHG annulus minimum",
        )
    if buffer_radius > spec.support_minimum:
        raise SGBLFamilyPrincipalStop(
            "domain_error",
            "Minkowski-buffer boxes cannot be taken from the compact support",
        )
    point = SGBLSourceInputs(
        coordinate_radius=buffer_radius,
        alpha=(1, 0, 0, 0, 0),
        shift=(0, 0, 0, 0, 0),
        radial_metric=(1, 0, 0, 0, 0),
        areal_radius=(buffer_radius, 0, 1, 0, 0),
        phi=(0, 0, 0, 0, 0),
        chi=(0, 0, 0, 0, 0),
        planck_mass=spec.planck_mass,
        scalar_mass=spec.scalar_mass,
        quartic_coupling=spec.quartic_coupling,
        alpha_gb=spec.alpha_gb,
    )
    state = sgbl_source_state(point, (0, 0, 0, 0, 0, 0))
    box = sgbl_principal_box_from_spherical_state(state, limits=limits)
    if not box.curvature_free:
        raise SGBLFamilyPrincipalStop(
            "domain_error",
            "exact Minkowski buffer must produce a curvature-free principal box",
        )
    return box


def _fixture_a_inputs() -> SGBLSourceInputs:
    return SGBLSourceInputs(
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


def sgbl_nonflat_fixture_a_interval_principal_box(
    *,
    limits: SGBLConeLimits | None = None,
) -> SGBLPrincipalBackgroundBox:
    """Exact-interval orthonormal box of locked source fixture A.

    Coordinate Riemann and Hess(phi) are the existing exact two-jet owners.
    The orthonormal image uses checked rational square-root factors, not
    binary64-as-exact conversion.
    """

    point = _fixture_a_inputs()
    state = sgbl_source_state(point, sgbl_source_solve(point))
    return sgbl_principal_box_from_spherical_state(state, limits=limits)


def _analytic_exterior_available(record: SGBLContinuousCompactnessRecord) -> bool:
    return (
        record.continuous_no_initial_trapped_sphere
        and record.ode is not None
        and record.ode.obstruction is None
        and bool(record.ode.cells)
        and record.ode.cells[-1].radius.upper == record.slice.support_maximum
    )


def _graph_covers_support(record: SGBLContinuousCompactnessRecord) -> bool:
    if (
        record.ode is None
        or record.ode.obstruction is not None
        or not record.ode.cells
        or not record.continuous_no_initial_trapped_sphere
    ):
        return False
    cells = record.ode.cells
    if cells[0].radius.lower != record.slice.support_minimum:
        return False
    if cells[-1].radius.upper != record.slice.support_maximum:
        return False
    return all(
        cell.picard_strict_self_map
        and cell.residual_contains_origin
        and cell.compactness_invariant_contains_zero
        for cell in cells
    )


def _cell_strictly_interior(cell: SGBLODECellEnclosure, spec: SGBLExactInitialSlice) -> bool:
    scaled = (cell.radius - spec.center) / spec.half_width
    return scaled.lower > -1 and scaled.upper < 1


@dataclass(frozen=True, slots=True, kw_only=True)
class SGBLFamilyIntervalJets:
    """Complete interval lower two-jet slots on one authenticated region."""

    radius: Interval
    radial_metric: Interval
    angular_extrinsic_curvature: Interval
    fields: Mapping[str, Interval]
    coverage: str

    def __post_init__(self) -> None:
        object.__setattr__(self, "fields", MappingProxyType(dict(self.fields)))
        required = {
            "chi",
            "chi_r",
            "chi_rr",
            "chi_pi",
            "chi_pi_r",
            "chi_tr",
            "lambda_r",
            "k_r",
            "lambda_rr",
            "k_rr",
            "phi_rrr",
        }
        missing = required - set(self.fields)
        if missing:
            raise ValueError(f"family interval jets missing slots: {sorted(missing)}")


@dataclass(frozen=True, slots=True, kw_only=True)
class SGBLFamilyPrincipalCellRecord:
    """Principal-box attempt on one authenticated ODE cell."""

    cell: SGBLODECellEnclosure
    jets: SGBLFamilyIntervalJets | None
    source_point: SGBLSourceInputs | None
    admission: SGBLSourceAdmissionRecord | None
    box: SGBLPrincipalBackgroundBox | None
    coverage: str
    classification: str
    missing_owner: SGBLFamilyPrincipalMissingOwner | None

    def __post_init__(self) -> None:
        if type(self.cell) is not SGBLODECellEnclosure:
            raise TypeError("cell must be SGBLODECellEnclosure")
        if self.classification not in {
            "principal_box",
            "missing_owner",
            "source_admission_incomplete",
        }:
            raise ValueError("unknown family-cell principal classification")
        if self.classification == "principal_box" and type(self.box) is not SGBLPrincipalBackgroundBox:
            raise ValueError("principal_box classification requires a cone box")
        if self.classification != "principal_box" and self.missing_owner is None:
            raise ValueError("failed cell records must retain the missing owner")


@dataclass(frozen=True, slots=True, kw_only=True)
class SGBLFamilyPrincipalRecord:
    """Aggregate matched-family principal feeder record."""

    slice: SGBLExactInitialSlice
    compactness: SGBLContinuousCompactnessRecord
    buffer_box: SGBLPrincipalBackgroundBox | None
    cell_records: tuple[SGBLFamilyPrincipalCellRecord, ...]
    exterior_box: SGBLPrincipalBackgroundBox | None
    classification: str
    inconclusive_reason: str | None
    missing_theorem: str | None
    missing_owners: tuple[SGBLFamilyPrincipalMissingOwner, ...]
    theorem: str

    def __post_init__(self) -> None:
        if type(self.slice) is not SGBLExactInitialSlice:
            raise TypeError("slice must be SGBLExactInitialSlice")
        if type(self.compactness) is not SGBLContinuousCompactnessRecord:
            raise TypeError("compactness must be SGBLContinuousCompactnessRecord")
        if self.classification not in {
            "family_principal_boxes",
            "family_geometry_incomplete",
        }:
            raise ValueError("unknown family-principal classification")
        if self.classification == "family_geometry_incomplete":
            if self.exterior_box is not None:
                raise ValueError("incomplete family graphs cannot fabricate an exterior box")
            if self.inconclusive_reason != "family_geometry_incomplete":
                raise ValueError("incomplete family graphs retain family_geometry_incomplete")
        object.__setattr__(self, "cell_records", tuple(self.cell_records))
        object.__setattr__(self, "missing_owners", tuple(self.missing_owners))

    @property
    def SGBL_branch_owned_and_healthy(self) -> bool:
        return False

    @property
    def FRZ1(self) -> bool:
        return False

    @property
    def PREF1(self) -> bool:
        return False

    @property
    def execution_authorized(self) -> bool:
        return False

    @property
    def holdout_authorized(self) -> bool:
        return False


def sgbl_family_cell_principal_box(
    cell: SGBLODECellEnclosure,
    spec: SGBLExactInitialSlice,
    *,
    acceleration_half_width: Fraction | int = FAMILY_ACCELERATION_HALF_WIDTH,
    limits: SGBLConeLimits | None = None,
    interval_limits: SGBLIntervalLimits | None = None,
) -> SGBLFamilyPrincipalCellRecord:
    """Build a principal box from one authenticated Picard cell."""

    if type(cell) is not SGBLODECellEnclosure:
        raise TypeError("cell must be SGBLODECellEnclosure")
    if type(spec) is not SGBLExactInitialSlice:
        raise TypeError("spec must be SGBLExactInitialSlice")
    if not (
        cell.picard_strict_self_map
        and cell.residual_contains_origin
        and cell.compactness_invariant_contains_zero
    ):
        missing = SGBLFamilyPrincipalMissingOwner(
            slot="validated_ODE_cell",
            owner="sgbl_validated_constraint_ode",
            reason="unauthenticated_cell",
            theorem="strict_Picard_self_map_with_affine_residuals_and_compactness_invariant",
        )
        return SGBLFamilyPrincipalCellRecord(
            cell=cell,
            jets=None,
            source_point=None,
            admission=None,
            box=None,
            coverage="none",
            classification="missing_owner",
            missing_owner=missing,
        )
    affine = sgbl_family_affine_jets(
        cell.radius, cell.lambda_box, cell.k_box, spec
    )
    jets = SGBLFamilyIntervalJets(
        radius=cell.radius,
        radial_metric=cell.lambda_box,
        angular_extrinsic_curvature=cell.k_box,
        fields=affine,
        coverage="cell_product_lower_jets",
    )
    radius = _witness(cell.radius)
    radial_metric = _witness(cell.lambda_box)
    angular_k = _witness(cell.k_box)
    point, _affine_at_center = _source_inputs_from_witness(
        radius=radius,
        radial_metric=radial_metric,
        angular_extrinsic_curvature=angular_k,
        spec=spec,
    )
    covering_owner = SGBLFamilyPrincipalMissingOwner(
        slot="dtt_over_cell_product_box",
        owner="sgbl_parametric_source_admission",
        reason="uniform_covering_box_not_proved",
        theorem=(
            "parametric_Krawczyk_uniqueness_on_the_non_uniform_family_cell_product_box"
        ),
    )
    try:
        admission = sgbl_family_cell_source_admission(
            point,
            parameter_half_width=0,
            acceleration_half_width=acceleration_half_width,
            limits=interval_limits,
        )
    except (SGBLAdmissionInconclusive, SGBLSourceJacobianSolveStop) as exc:
        reason = getattr(exc, "reason", "source_admission_incomplete")
        missing = SGBLFamilyPrincipalMissingOwner(
            slot="dtt",
            owner="sgbl_parametric_source_admission",
            reason=str(reason),
            theorem="parametric_Krawczyk_unique_acceleration_root_at_the_declared_family_witness",
        )
        return SGBLFamilyPrincipalCellRecord(
            cell=cell,
            jets=jets,
            source_point=point,
            admission=None,
            box=None,
            coverage="declared_witness_in_cell_product_box",
            classification="source_admission_incomplete",
            missing_owner=missing,
        )
    if not admission.unique_acceleration_root_for_every_declared_parameter_point:
        missing = SGBLFamilyPrincipalMissingOwner(
            slot="dtt",
            owner="sgbl_parametric_source_admission",
            reason=str(admission.inconclusive_reason or "source_admission_incomplete"),
            theorem="parametric_Krawczyk_unique_acceleration_root_at_the_declared_family_witness",
        )
        return SGBLFamilyPrincipalCellRecord(
            cell=cell,
            jets=jets,
            source_point=point,
            admission=admission,
            box=None,
            coverage="declared_witness_in_cell_product_box",
            classification="source_admission_incomplete",
            missing_owner=missing,
        )
    state = _full_two_jet_state(point, admission.acceleration_box, spec)
    box = sgbl_principal_box_from_spherical_state(state, limits=limits)
    return SGBLFamilyPrincipalCellRecord(
        cell=cell,
        jets=jets,
        source_point=point,
        admission=admission,
        box=box,
        coverage="declared_witness_in_cell_product_box",
        classification="principal_box",
        missing_owner=covering_owner,
    )


def _exterior_boundary_box(
    record: SGBLContinuousCompactnessRecord,
    *,
    limits: SGBLConeLimits | None = None,
    interval_limits: SGBLIntervalLimits | None = None,
) -> SGBLPrincipalBackgroundBox | SGBLFamilyPrincipalMissingOwner:
    spec = record.slice
    end = record.ode.cells[-1]
    try:
        cell_record = sgbl_family_cell_principal_box(
            end,
            spec,
            limits=limits,
            interval_limits=interval_limits,
        )
    except SGBLFamilyPrincipalStop as exc:
        return SGBLFamilyPrincipalMissingOwner(
            slot="exterior_boundary_box",
            owner="analytic_GR_vacuum_continuation",
            reason=exc.reason,
            theorem="authenticated_support_boundary_state_with_C_equals_two_M_over_r",
        )
    if cell_record.box is None:
        return cell_record.missing_owner or SGBLFamilyPrincipalMissingOwner(
            slot="exterior_boundary_box",
            owner="analytic_GR_vacuum_continuation",
            reason="source_admission_incomplete",
            theorem="authenticated_support_boundary_state_with_C_equals_two_M_over_r",
        )
    return cell_record.box


def sgbl_family_principal_from_health(
    record: SGBLContinuousCompactnessRecord,
    *,
    spec: SGBLExactInitialSlice | None = None,
    max_cells: int = DEFAULT_MAX_FAMILY_PRINCIPAL_CELLS,
    limits: SGBLConeLimits | None = None,
    interval_limits: SGBLIntervalLimits | None = None,
    include_buffer: bool = True,
) -> SGBLFamilyPrincipalRecord:
    """Feed principal boxes from an authenticated compactness record.

    Incomplete graphs return one typed ``family_geometry_incomplete``
    record and do not fabricate an exterior or set health true.
    """

    if type(record) is not SGBLContinuousCompactnessRecord:
        raise TypeError("record must be SGBLContinuousCompactnessRecord")
    slice_spec = record.slice
    if spec is not None:
        if type(spec) is not SGBLExactInitialSlice:
            raise TypeError("spec must be SGBLExactInitialSlice")
        if spec != slice_spec:
            raise SGBLFamilyPrincipalStop(
                "domain_error",
                "family specification does not match the authenticated health record",
            )
    max_cells = _positive_int("max_cells", max_cells)
    theorem = (
        "Authenticated Minkowski-buffer, compact-support, and exterior principal "
        "boxes are built from exact interval two-jets, affine constraint "
        "differentiation, and a parametric source-admission certificate.  "
        "Binary64 values are not treated as enclosures.  Missing derivative or "
        "source-box owners are named rather than filled by zero."
    )
    if not _graph_covers_support(record):
        return SGBLFamilyPrincipalRecord(
            slice=slice_spec,
            compactness=record,
            buffer_box=None,
            cell_records=(),
            exterior_box=None,
            classification="family_geometry_incomplete",
            inconclusive_reason="family_geometry_incomplete",
            missing_theorem=(
                "complete_validated_constraint_ODE_graph_covering_the_compact_support"
            ),
            missing_owners=(
                SGBLFamilyPrincipalMissingOwner(
                    slot="compact_support_graph",
                    owner="sgbl_continuous_initial_compactness",
                    reason="family_geometry_incomplete",
                    theorem="complete_validated_constraint_ODE_graph_covering_the_compact_support",
                ),
            ),
            theorem=(
                "The compactness record does not cover the compact support with a "
                "validated Picard graph.  Exterior vacuum continuation and branch "
                "health are not fabricated from the incomplete cells.  "
                f"Obstruction: {record.inconclusive_reason}."
            ),
        )
    buffer_box = (
        sgbl_minkowski_buffer_principal_box(slice_spec, limits=limits)
        if include_buffer
        else None
    )
    ordered = tuple(
        cell
        for cell in record.ode.cells
        if _cell_strictly_interior(cell, slice_spec)
    ) + tuple(
        cell
        for cell in record.ode.cells
        if not _cell_strictly_interior(cell, slice_spec)
    )
    cell_records: list[SGBLFamilyPrincipalCellRecord] = []
    missing_owners: list[SGBLFamilyPrincipalMissingOwner] = []
    for cell in ordered:
        if len(cell_records) >= max_cells:
            break
        cell_record = sgbl_family_cell_principal_box(
            cell,
            slice_spec,
            limits=limits,
            interval_limits=interval_limits,
        )
        cell_records.append(cell_record)
        if cell_record.missing_owner is not None:
            missing_owners.append(cell_record.missing_owner)
    exterior_box: SGBLPrincipalBackgroundBox | None = None
    if _analytic_exterior_available(record):
        exterior = _exterior_boundary_box(
            record, limits=limits, interval_limits=interval_limits
        )
        if type(exterior) is SGBLPrincipalBackgroundBox:
            exterior_box = exterior
        elif type(exterior) is SGBLFamilyPrincipalMissingOwner:
            missing_owners.append(exterior)
    return SGBLFamilyPrincipalRecord(
        slice=slice_spec,
        compactness=record,
        buffer_box=buffer_box,
        cell_records=tuple(cell_records),
        exterior_box=exterior_box,
        classification="family_principal_boxes",
        inconclusive_reason=None,
        missing_theorem=None if not missing_owners else missing_owners[0].theorem,
        missing_owners=tuple(missing_owners),
        theorem=theorem,
    )


def sgbl_family_principal_health_gate(
    record: SGBLFamilyPrincipalRecord,
) -> dict[str, Any]:
    """Aggregate production flags remain closed after a local principal box."""

    if type(record) is not SGBLFamilyPrincipalRecord:
        raise TypeError("record must be SGBLFamilyPrincipalRecord")
    return {
        "SGBL_branch_owned_and_healthy": False,
        "FRZ1": False,
        "PREF1": False,
        "execution_authorized": False,
        "holdout_authorized": False,
        "family_geometry_incomplete": record.classification == "family_geometry_incomplete",
        "copied_gr0_or_fgcqr_health_evidence": False,
        "float_used_as_enclosure": False,
        "missing_health_closes_the_gate": True,
        "buffer_box_emitted": record.buffer_box is not None,
        "exterior_box_emitted": record.exterior_box is not None,
        "principal_cell_boxes": sum(
            1 for cell in record.cell_records if cell.box is not None
        ),
    }


__all__ = [
    "BUMP_XXX_X_DERIVATIVE_BOUND",
    "FAMILY_ACCELERATION_HALF_WIDTH",
    "FORBIDDEN_HEALTH_IMPORTS",
    "SLOT_OWNERSHIP",
    "SGBLFamilyIntervalJets",
    "SGBLFamilyPrincipalCellRecord",
    "SGBLFamilyPrincipalMissingOwner",
    "SGBLFamilyPrincipalRecord",
    "SGBLFamilyPrincipalStop",
    "sgbl_checked_orthonormal_pullback",
    "sgbl_compact_bump_third_derivative_enclosure",
    "sgbl_complete_interval_family_fields",
    "sgbl_family_adm_lower_jets",
    "sgbl_family_affine_jets",
    "sgbl_family_cell_principal_box",
    "sgbl_family_cell_source_admission",
    "sgbl_family_principal_from_health",
    "sgbl_family_principal_health_gate",
    "sgbl_family_principal_slot_table",
    "sgbl_interval_orthonormal_frame",
    "sgbl_interval_scalar_hessian",
    "sgbl_minkowski_buffer_principal_box",
    "sgbl_nonflat_fixture_a_interval_principal_box",
    "sgbl_principal_box_from_spherical_state",
]
