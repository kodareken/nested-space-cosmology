"""Exact pointwise SGB-L source residual and acceleration Jacobian.

At fixed lower jets the complete six-row MHG residual of the linear branch
``f(phi)=alpha_gb*phi``, ``beta=eta=0``, is affine in the six ADM
coordinate-time accelerations

    a = (alpha_tt, shift_tt, lambda_tt, R_tt, phi_tt, chi_tt).

That identity follows from the highest-time curvature cancellations already
proved: it is not a sample interpolation, a claim that the map from lower
jets to accelerations is linear, or a proof that ``J`` is invertible.  The
reference instrument is RED1's full MHG tensor residual, not a substituted
evolution system.  ``R0`` and ``J`` are read from seven exact constant and
basis evaluations.  The optional algebraic solve re-evaluates that same
residual at the root and requires the exact zero vector.

This is finite pointwise algebra only.  It does not construct a matched
family, continue through the centre, certify interval invertibility, supply
a source runtime, or prove branch kinetic or cone health.  The six-by-six
source Jacobian is not the two-dimensional annular constraint Jacobian.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from fractions import Fraction
from typing import Sequence

from .exact_linear_algebra import determinant, inverse, matrix_multiply
from .modified_harmonic_reference import (
    MHG2_FULL_EQUATION_ORDER,
    modified_harmonic_full_residuals,
)
from .reference_connection import flat_spherical_annulus_reference
from .spherical_reduction import Jet2, SphericalState, _state_from_adm_pg_jets


Q = Fraction

SOURCE_ACCELERATION_ORDER = (
    "alpha_tt",
    "shift_tt",
    "lambda_tt",
    "R_tt",
    "phi_tt",
    "chi_tt",
)
SOURCE_EQUATION_ORDER = MHG2_FULL_EQUATION_ORDER
SOURCE_FIELD_ORDER = ("alpha", "shift", "lambda", "R", "phi", "chi")
LOWER_JET_COMPONENTS = ("value", "dt", "dr", "dtr", "drr")
REFERENCE_RADIAL_MINIMUM = Q(1, 2)
TILDE_NORMAL_FACTOR = 4
HAT_NORMAL_FACTOR = 9
_ZERO_ACCELERATION = (Q(0), Q(0), Q(0), Q(0), Q(0), Q(0))
_SOURCE_REFERENCE = flat_spherical_annulus_reference(
    radial_domain_minimum=REFERENCE_RADIAL_MINIMUM,
)


def _fraction(name: str, value: object) -> Fraction:
    if type(value) not in (int, Fraction):
        raise TypeError(f"{name} must be a Fraction or built-in int")
    return Fraction(value)


def _lower_jet(name: str, value: object) -> Jet2:
    """Accept an exact Jet2 with ``dtt=0`` or the five-tuple excluding ``dtt``."""

    if type(value) is Jet2:
        for component in ("value", "dt", "dr", "dtt", "dtr", "drr"):
            if type(getattr(value, component)) is not Fraction:
                raise TypeError(f"{name}.{component} must be a Fraction")
        if value.dtt != 0:
            raise ValueError(
                f"{name} lower jet must not include a nonzero dtt; "
                "accelerations are supplied separately"
            )
        return value
    if isinstance(value, (tuple, list)):
        if len(value) == 6:
            raise ValueError(
                f"{name} lower jet excludes dtt; use the five-tuple "
                f"{LOWER_JET_COMPONENTS} or a Jet2 with dtt=0"
            )
        if len(value) != 5:
            raise TypeError(
                f"{name} must be Jet2 or a 5-tuple {LOWER_JET_COMPONENTS}"
            )
        parts = tuple(
            _fraction(f"{name}.{component}", item)
            for component, item in zip(LOWER_JET_COMPONENTS, value, strict=True)
        )
        return Jet2(parts[0], dt=parts[1], dr=parts[2], dtt=0, dtr=parts[3], drr=parts[4])
    raise TypeError(f"{name} must be Jet2 or a 5-tuple {LOWER_JET_COMPONENTS}")


def _six(name: str, value: object) -> tuple[Fraction, ...]:
    if not isinstance(value, (tuple, list)) or len(value) != 6:
        raise TypeError(f"{name} must be a 6-tuple of Fraction or built-in int")
    return tuple(_fraction(f"{name}[{index}]", item) for index, item in enumerate(value))


def _jacobian(name: str, value: object) -> tuple[tuple[Fraction, ...], ...]:
    if not isinstance(value, (tuple, list)) or len(value) != 6:
        raise TypeError(f"{name} must be a 6-by-6 tuple of tuples")
    return tuple(_six(f"{name}[{index}]", row) for index, row in enumerate(value))


def _apply_jacobian(
    jacobian: Sequence[Sequence[Fraction]],
    accelerations: Sequence[Fraction],
) -> tuple[Fraction, ...]:
    return tuple(
        sum(row[column] * accelerations[column] for column in range(6))
        for row in jacobian
    )


@dataclass(frozen=True, slots=True, kw_only=True)
class SGBLSourceInputs:
    """Exact annular source-point data with coordinate-time accelerations excluded.

    Jet fields accept only :class:`Jet2` with ``dtt=0`` or the five-tuple
    ``(value, dt, dr, dtr, drr)``.  A supplied nonzero ``dtt`` is refused
    rather than dropped.  Lapse ``alpha``, ADM ``lambda`` (``radial_metric``),
    areal radius ``R``, coordinate radius, ``planck_mass``, ``scalar_mass``
    (``mu``), and ``quartic_coupling`` (``g4``) are positive.
    ``planck_mass`` is the mass, not its square.  ``alpha_gb=0`` is allowed
    only as an algebraic coupling-limit control, not as a production SGB-L
    branch.  The frozen MHG reference is the flat spherical annulus with
    minimum ``1/2``, ``tilde=4``, ``hat=9``, and ``beta=eta=0``.
    """

    coordinate_radius: Fraction | int
    alpha: Jet2 | tuple[object, ...]
    shift: Jet2 | tuple[object, ...]
    radial_metric: Jet2 | tuple[object, ...]
    areal_radius: Jet2 | tuple[object, ...]
    phi: Jet2 | tuple[object, ...]
    chi: Jet2 | tuple[object, ...]
    planck_mass: Fraction | int
    scalar_mass: Fraction | int
    quartic_coupling: Fraction | int
    alpha_gb: Fraction | int

    def __post_init__(self) -> None:
        for name in (
            "coordinate_radius",
            "planck_mass",
            "scalar_mass",
            "quartic_coupling",
            "alpha_gb",
        ):
            object.__setattr__(self, name, _fraction(name, getattr(self, name)))
        for name in (
            "alpha",
            "shift",
            "radial_metric",
            "areal_radius",
            "phi",
            "chi",
        ):
            object.__setattr__(self, name, _lower_jet(name, getattr(self, name)))
        for name, value in (
            ("coordinate_radius", self.coordinate_radius),
            ("alpha", self.alpha.value),
            ("radial_metric", self.radial_metric.value),
            ("areal_radius", self.areal_radius.value),
            ("planck_mass", self.planck_mass),
            ("scalar_mass", self.scalar_mass),
            ("quartic_coupling", self.quartic_coupling),
        ):
            if value <= 0:
                raise ValueError(f"{name} must be positive")


@dataclass(frozen=True, slots=True, kw_only=True)
class SGBLSourceCoefficients:
    """``R0`` and the six-by-six source Jacobian ``J`` in ``R(a)=R0+J a``.

    The record stores only those exact arrays.  Coefficient identities and
    solved-residual claims are not fields and cannot be attached by
    construction or ``replace``.
    """

    constant: tuple[Fraction | int, ...]
    jacobian: tuple[tuple[Fraction | int, ...], ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "constant", _six("constant", self.constant))
        object.__setattr__(self, "jacobian", _jacobian("jacobian", self.jacobian))

    @property
    def equation_order(self) -> tuple[str, ...]:
        return SOURCE_EQUATION_ORDER

    @property
    def acceleration_order(self) -> tuple[str, ...]:
        return SOURCE_ACCELERATION_ORDER

    @property
    def jacobian_determinant(self) -> Fraction:
        return determinant(self.jacobian)

    def evaluate(
        self,
        accelerations: Sequence[Fraction | int],
    ) -> tuple[Fraction, ...]:
        """Return the affine reconstruction ``R0+J a`` at exact accelerations."""

        acceleration = _six("accelerations", accelerations)
        applied = _apply_jacobian(self.jacobian, acceleration)
        return tuple(
            constant + linear
            for constant, linear in zip(self.constant, applied, strict=True)
        )

    def solve(self) -> tuple[Fraction, ...]:
        """Return the exact algebraic root of ``R0+J a=0``, or refuse ``det J=0``.

        Nonzero determinant is an exact algebraic condition, not a
        quantitative Jacobian floor, a neighborhood or existence theorem, or
        a kinetic-health bound.  This linear solve is not a re-evaluation of
        the tensor residual; :func:`sgbl_source_solve` performs that check.
        """

        if self.jacobian_determinant == 0:
            raise SGBLSourceJacobianSolveStop(self)
        inverse_jacobian = inverse(self.jacobian)
        rhs = tuple((-entry,) for entry in self.constant)
        product = matrix_multiply(inverse_jacobian, rhs)
        return tuple(row[0] for row in product)


class SGBLSourceJacobianSolveStop(ArithmeticError):
    """An exactly singular six-by-six source Jacobian; no state is advanced."""

    reason = "singular_source_jacobian"

    def __init__(self, coefficients: SGBLSourceCoefficients) -> None:
        self.coefficients = coefficients
        super().__init__("singular source Jacobian")


def sgbl_source_state(
    point: SGBLSourceInputs,
    accelerations: Sequence[Fraction | int] = _ZERO_ACCELERATION,
) -> SphericalState:
    """Map an exact source point and ADM accelerations to a spherical two-jet."""

    if type(point) is not SGBLSourceInputs:
        raise TypeError("point must be SGBLSourceInputs")
    acceleration = _six("accelerations", accelerations)
    jets = {
        "alpha": replace(point.alpha, dtt=acceleration[0]),
        "shift": replace(point.shift, dtt=acceleration[1]),
        "lambda": replace(point.radial_metric, dtt=acceleration[2]),
        "areal_radius": replace(point.areal_radius, dtt=acceleration[3]),
        "phi": replace(point.phi, dtt=acceleration[4]),
        "chi": replace(point.chi, dtt=acceleration[5]),
    }
    return _state_from_adm_pg_jets(
        {
            "model_id": "SGB-L" if point.alpha_gb else "GR-0",
            "action_parameters": {
                "planck_mass": point.planck_mass,
                "scalar_mass": point.scalar_mass,
                "quartic_coupling": point.quartic_coupling,
                "alpha_gb": point.alpha_gb,
            },
        },
        jets,
    )


def sgbl_source_residual(
    point: SGBLSourceInputs,
    accelerations: Sequence[Fraction | int] = _ZERO_ACCELERATION,
) -> tuple[Fraction, ...]:
    """Evaluate the complete six-row MHG residual at exact accelerations."""

    state = sgbl_source_state(point, accelerations)
    return tuple(
        modified_harmonic_full_residuals(
            state,
            reference=_SOURCE_REFERENCE,
            coordinate_radius=point.coordinate_radius,
            tilde_normal_factor=TILDE_NORMAL_FACTOR,
            hat_normal_factor=HAT_NORMAL_FACTOR,
        )["full_residual_vector"]
    )


def sgbl_source_coefficients(point: SGBLSourceInputs) -> SGBLSourceCoefficients:
    """Read ``R0`` and ``J`` from seven exact constant and unit-basis residuals.

    The highest-time identity makes those seven evaluations exact for the
    affine map; they are not a numerical fit.  Scalar and chi momenta in the
    lower jet are retained.  Invertibility of ``J`` is a separate question
    from affinity.
    """

    if type(point) is not SGBLSourceInputs:
        raise TypeError("point must be SGBLSourceInputs")
    constant = sgbl_source_residual(point, _ZERO_ACCELERATION)
    columns: list[tuple[Fraction, ...]] = []
    for index in range(6):
        seed = list(_ZERO_ACCELERATION)
        seed[index] = Q(1)
        evaluated = sgbl_source_residual(point, seed)
        columns.append(
            tuple(
                entry - origin
                for entry, origin in zip(evaluated, constant, strict=True)
            )
        )
    jacobian = tuple(
        tuple(columns[column][row] for column in range(6))
        for row in range(6)
    )
    return SGBLSourceCoefficients(constant=constant, jacobian=jacobian)


def sgbl_source_solve(point: SGBLSourceInputs) -> tuple[Fraction, ...]:
    """Solve ``R0+J a=0`` and require the actual full residual to vanish.

    The tensor residual is re-evaluated at the algebraic root.  Exact zero
    is required; there is no conditioning floor.
    """

    accelerations = sgbl_source_coefficients(point).solve()
    residual = sgbl_source_residual(point, accelerations)
    if residual != _ZERO_ACCELERATION:
        raise RuntimeError(
            "full MHG residual is not exact zero at the algebraic source root"
        )
    return accelerations


__all__ = [
    "HAT_NORMAL_FACTOR",
    "LOWER_JET_COMPONENTS",
    "REFERENCE_RADIAL_MINIMUM",
    "SGBLSourceCoefficients",
    "SGBLSourceInputs",
    "SGBLSourceJacobianSolveStop",
    "SOURCE_ACCELERATION_ORDER",
    "SOURCE_EQUATION_ORDER",
    "SOURCE_FIELD_ORDER",
    "TILDE_NORMAL_FACTOR",
    "sgbl_source_coefficients",
    "sgbl_source_residual",
    "sgbl_source_solve",
    "sgbl_source_state",
]
