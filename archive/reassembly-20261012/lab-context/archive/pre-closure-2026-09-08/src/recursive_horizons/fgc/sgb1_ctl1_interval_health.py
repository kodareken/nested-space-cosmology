"""Declared-box invertibility of the SGB-L ADM source Jacobian.

At fixed lower jets the complete six-row MHG residual is affine in the ADM
accelerations ``(alpha_tt, shift_tt, lambda_tt, R_tt, phi_tt, chi_tt)``.
This module encloses that six-by-six map over a caller-declared lower-jet
box.  A Neumann inverse is emitted only when a predeclared infinity-norm
condition proves ``rho<1``.  Wide wrapping with ``rho>=1`` is typed
interval-inconclusive; it is not a fitted floor and not
``singular_source_jacobian``.

The instrument does not freeze a physical production width, claim
all-direction cone health, or run a trajectory.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from hashlib import sha256
from numbers import Integral
from types import MappingProxyType
from typing import Any, Mapping, Sequence

from .exact_interval import Interval, interval, interval_matrix_infinity_row_sum_bound
from .exact_interval_krawczyk import parametric_krawczyk_inclusion
from .exact_interval_linear_algebra import center_preconditioned_neumann_inverse
from .exact_linear_algebra import determinant, inverse
from .interval_tangent import IntervalFirstTangent, primal_and_tangent
from .modified_harmonic_reference import modified_harmonic_full_residuals
from .reference_connection import (
    ReferenceConnection,
    flat_spherical_annulus_reference,
)
from .sgb1_ctl1_source import (
    HAT_NORMAL_FACTOR,
    REFERENCE_RADIAL_MINIMUM,
    SGBLSourceInputs,
    SGBLSourceJacobianSolveStop,
    SOURCE_ACCELERATION_ORDER,
    SOURCE_EQUATION_ORDER,
    TILDE_NORMAL_FACTOR,
    sgbl_source_coefficients,
)
from .spherical_reduction import ADM_FIELD_ORDER, Jet2, _state_from_adm_pg_jets


Q = Fraction
ADM_SOURCE_CHART = "adm_accelerations"
LOWER_JET_BOX_SLOTS = ("value", "dt", "dr", "dtr", "drr")
ADM_JET_ATTR = {
    "alpha": "alpha",
    "shift": "shift",
    "lambda": "radial_metric",
    "areal_radius": "areal_radius",
    "phi": "phi",
    "chi": "chi",
}
DEFAULT_MAX_RESIDUAL_EVALUATIONS = 7
DEFAULT_MAX_RATIONAL_BIT_LENGTH = 16384
DEFAULT_MAX_PARAMETER_HALF_WIDTH = Q(1)
_CHART_FAILURE_MARKERS = (
    "two-dimensional base must be Lorentzian",
    "areal_radius must be positive",
    "singular metric",
    "singular jet reciprocal",
    "ADM/PG alpha and lambda must be positive",
    "effective Planck coefficient F must be positive",
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


def _interval_bit_length(value: Interval) -> int:
    return max(
        abs(value.lower.numerator).bit_length(),
        value.lower.denominator.bit_length(),
        abs(value.upper.numerator).bit_length(),
        value.upper.denominator.bit_length(),
    )


def _matrix_sha256(matrix: Sequence[Sequence[Fraction]]) -> str:
    payload = ",".join(
        f"{entry.numerator}/{entry.denominator}"
        for row in matrix
        for entry in row
    ).encode()
    return sha256(payload).hexdigest()


@dataclass(frozen=True, slots=True, kw_only=True)
class SGBLIntervalLimits:
    """Arithmetic and resource caps. These are not physical box widths."""

    max_residual_evaluations: int = DEFAULT_MAX_RESIDUAL_EVALUATIONS
    max_rational_bit_length: int = DEFAULT_MAX_RATIONAL_BIT_LENGTH
    max_parameter_half_width: Fraction | int = DEFAULT_MAX_PARAMETER_HALF_WIDTH

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "max_residual_evaluations",
            _positive_int(
                "max_residual_evaluations", self.max_residual_evaluations
            ),
        )
        object.__setattr__(
            self,
            "max_rational_bit_length",
            _positive_int(
                "max_rational_bit_length", self.max_rational_bit_length
            ),
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


@dataclass(frozen=True, slots=True, kw_only=True)
class SGBLIntervalBox:
    """Declared lower-jet product box around an exact SGB-L source point.

    ``parameter_half_width`` is applied to every non-acceleration jet slot.
    The caller declares it; this owner does not fit it from a condition
    number.  Accelerations stay exact unless a positive
    ``acceleration_half_width`` is supplied for an optional unique-root
    Krawczyk test.
    """

    center: SGBLSourceInputs
    parameter_half_width: Fraction | int
    acceleration_half_width: Fraction | int = 0
    limits: SGBLIntervalLimits = SGBLIntervalLimits()
    allow_zero_coupling_control: bool = False
    chart: str = ADM_SOURCE_CHART

    def __post_init__(self) -> None:
        if type(self.center) is not SGBLSourceInputs:
            raise TypeError("center must be SGBLSourceInputs")
        object.__setattr__(
            self,
            "parameter_half_width",
            _fraction(
                "parameter_half_width",
                self.parameter_half_width,
                nonnegative=True,
            ),
        )
        object.__setattr__(
            self,
            "acceleration_half_width",
            _fraction(
                "acceleration_half_width",
                self.acceleration_half_width,
                nonnegative=True,
            ),
        )
        if type(self.limits) is not SGBLIntervalLimits:
            raise TypeError("limits must be SGBLIntervalLimits")
        if type(self.allow_zero_coupling_control) is not bool:
            raise TypeError("allow_zero_coupling_control must be boolean")
        if self.chart != ADM_SOURCE_CHART:
            raise SGBLIntervalInconclusive(
                "interval_chart_error",
                "SGB-L source Jacobian is d(MHG rows)/d(ADM accelerations)",
            )
        if self.center.alpha_gb == 0 and not self.allow_zero_coupling_control:
            raise ValueError(
                "production SGB-L box requires nonzero alpha_gb; "
                "zero coupling is only a named control"
            )
        if self.parameter_half_width > self.limits.max_parameter_half_width:
            raise SGBLIntervalInconclusive(
                "resource_limit",
                "declared parameter half-width exceeds the arithmetic cap",
            )

    @property
    def coordinate_radius(self) -> Fraction:
        return self.center.coordinate_radius

    @property
    def is_singleton(self) -> bool:
        return self.parameter_half_width == 0


@dataclass(frozen=True, slots=True, kw_only=True)
class SGBLIntervalHealthRecord:
    """Enclosure of the ADM source Jacobian on one declared box.

    Invertibility, unique-root and cone/health bits are independent.  Health
    and pass fields are derived and false.
    """

    box: SGBLIntervalBox
    branch: str
    action: Mapping[str, Fraction]
    reference: ReferenceConnection
    tilde_normal_factor: Fraction
    hat_normal_factor: Fraction
    exact_center_jacobian: tuple[tuple[Fraction, ...], ...]
    exact_center_det: Fraction
    exact_center_inverse: tuple[tuple[Fraction, ...], ...] | None
    exact_center_jacobian_sha256: str
    jacobian_box: tuple[tuple[Interval, ...], ...]
    remainder_bound: Fraction
    rho_infinity: Fraction | None
    inverse_tail_entrywise_bound: Fraction | None
    source_inverse: Mapping[str, Any] | None
    source_invertible_over_box: bool
    unique_acceleration_root_over_box: bool
    independent_exact_singleton_reduction: bool
    residual_evaluations: int
    max_rational_bit_length_observed: int
    classification: str
    inconclusive_reason: str | None
    spherical_ref1_blocks: None = None
    principal_cone_margins: None = None
    constraint_center_valid: None = None

    def __post_init__(self) -> None:
        if type(self.box) is not SGBLIntervalBox:
            raise TypeError("box must be SGBLIntervalBox")
        if self.branch != "SGB-L":
            raise ValueError("interval record is owned by the linear branch")
        if self.source_invertible_over_box and self.source_inverse is None:
            raise ValueError("proved invertibility must retain the inverse payload")
        if self.classification not in {
            "neumann_inverse_enclosure",
            "interval_inconclusive",
            "exact_singleton_inverse",
        }:
            raise ValueError("unknown interval classification")
        exact = tuple(tuple(Fraction(value) for value in row) for row in self.exact_center_jacobian)
        if len(exact) != 6 or any(len(row) != 6 for row in exact):
            raise ValueError("exact center Jacobian must be6x6")
        object.__setattr__(self, "exact_center_jacobian", exact)
        if determinant(exact) != self.exact_center_det:
            raise ValueError("exact center determinant does not match the Jacobian")
        if _matrix_sha256(exact) != self.exact_center_jacobian_sha256:
            raise ValueError("exact center Jacobian hash differs")
        if self.exact_center_det == 0 or self.exact_center_inverse is None:
            raise ValueError("nonsingular interval record must retain the center inverse")
        if dict(self.action) != _action_mapping(self.box.center):
            raise ValueError("interval action differs from the box center")
        object.__setattr__(self, "action", MappingProxyType(dict(self.action)))
        jacobian_box = tuple(tuple(value for value in row) for row in self.jacobian_box)
        if len(jacobian_box) != 6 or any(
            len(row) != 6 or any(type(value) is not Interval for value in row)
            for row in jacobian_box
        ):
            raise ValueError("interval Jacobian must be a6x6 Interval matrix")
        for row in range(6):
            for column in range(6):
                interval_value = jacobian_box[row][column]
                if not interval_value.lower <= exact[row][column] <= interval_value.upper:
                    raise ValueError("interval Jacobian does not contain its exact center")
        object.__setattr__(self, "jacobian_box", jacobian_box)
        deviation = tuple(
            tuple(
                jacobian_box[row][column] - Interval.singleton(exact[row][column])
                for column in range(6)
            )
            for row in range(6)
        )
        if interval_matrix_infinity_row_sum_bound(deviation) != self.remainder_bound:
            raise ValueError("interval remainder bound differs")
        if self.classification == "interval_inconclusive":
            if self.source_invertible_over_box or self.source_inverse is not None:
                raise ValueError("inconclusive interval cannot retain a proved inverse")
        else:
            if (
                not self.source_invertible_over_box
                or self.source_inverse is None
                or self.rho_infinity is None
                or not self.rho_infinity < 1
            ):
                raise ValueError("proved interval inverse requires retained rho<1 evidence")
        if self.classification == "exact_singleton_inverse" and (
            not self.box.is_singleton or not self.independent_exact_singleton_reduction
        ):
            raise ValueError("exact singleton classification differs from the box")
        if self.unique_acceleration_root_over_box and (
            not self.source_invertible_over_box or self.box.acceleration_half_width <= 0
        ):
            raise ValueError("unique-root claim lacks a proved inverse/acceleration box")
        if self.source_inverse is not None:
            object.__setattr__(
                self, "source_inverse", MappingProxyType(dict(self.source_inverse))
            )

    @property
    def acceleration_order(self) -> tuple[str, ...]:
        return SOURCE_ACCELERATION_ORDER

    @property
    def equation_order(self) -> tuple[str, ...]:
        return SOURCE_EQUATION_ORDER

    @property
    def strongly_hyperbolic(self) -> bool:
        return False

    @property
    def SGBL_branch_owned_and_healthy(self) -> bool:
        return False

    @property
    def cone_certificate(self) -> str:
        return "unqualified"


class SGBLIntervalInconclusive(ArithmeticError):
    """Typed stop when the box cannot support a proved inverse.

    This is not ``singular_source_jacobian`` and not a scientific nonpass.
    """

    def __init__(self, reason: str, message: str, payload: Mapping[str, Any] | None = None) -> None:
        self.reason = reason
        self.payload = dict(payload or {})
        super().__init__(message)


def _action_mapping(point: SGBLSourceInputs) -> dict[str, Fraction]:
    return {
        "planck_mass": point.planck_mass,
        "scalar_mass": point.scalar_mass,
        "quartic_coupling": point.quartic_coupling,
        "alpha_gb": point.alpha_gb,
        "beta": Q(0),
        "eta": Q(0),
    }


def _boxed_adm_jets(
    point: SGBLSourceInputs,
    *,
    half_width: Fraction,
    accelerations: Sequence[Fraction],
    seed_field: str | None,
) -> dict[str, Jet2]:
    fields: dict[str, Jet2] = {}
    for index, adm_name in enumerate(ADM_FIELD_ORDER):
        source = getattr(point, ADM_JET_ATTR[adm_name])
        entries: dict[str, object] = {}
        for slot in LOWER_JET_BOX_SLOTS:
            value = getattr(source, slot)
            entries[slot] = interval(value - half_width, value + half_width)
        dtt: object = Interval.singleton(accelerations[index])
        if seed_field == adm_name:
            dtt = IntervalFirstTangent.seed(dtt)
        entries["dtt"] = dtt
        fields[adm_name] = Jet2(**entries)
    return fields


def _interval_state(
    point: SGBLSourceInputs,
    *,
    half_width: Fraction,
    accelerations: Sequence[Fraction],
    seed_field: str | None = None,
):
    return _state_from_adm_pg_jets(
        {
            "model_id": "SGB-L",
            "action_parameters": _action_mapping(point),
        },
        _boxed_adm_jets(
            point,
            half_width=half_width,
            accelerations=accelerations,
            seed_field=seed_field,
        ),
    )


def _full_residual(point: SGBLSourceInputs, state: object) -> tuple[object, ...]:
    return tuple(
        modified_harmonic_full_residuals(
            state,
            reference=flat_spherical_annulus_reference(
                radial_domain_minimum=REFERENCE_RADIAL_MINIMUM,
            ),
            coordinate_radius=point.coordinate_radius,
            tilde_normal_factor=TILDE_NORMAL_FACTOR,
            hat_normal_factor=HAT_NORMAL_FACTOR,
        )["full_residual_vector"]
    )


def _raise_if_chart_or_reciprocal(exc: BaseException) -> None:
    message = str(exc)
    if isinstance(exc, ZeroDivisionError) or "contains zero" in message:
        raise SGBLIntervalInconclusive(
            "zero_in_interval_reciprocal",
            "interval reciprocal encountered a zero-containing denominator",
        ) from exc
    if any(marker in message for marker in _CHART_FAILURE_MARKERS):
        raise SGBLIntervalInconclusive(
            "interval_chart_domain",
            "declared box is not provably inside the ADM source chart",
        ) from exc


def _enforce_bits(values: Sequence[Interval], *, limit: int) -> int:
    observed = max(_interval_bit_length(entry) for entry in values)
    if observed > limit:
        raise SGBLIntervalInconclusive(
            "resource_limit",
            "interval endpoints exceeded the declared rational bit cap",
            {"observed": observed, "limit": limit},
        )
    return observed


def _singleton_equals_exact(
    jacobian_box: Sequence[Sequence[Interval]],
    exact: Sequence[Sequence[Fraction]],
) -> bool:
    return all(
        jacobian_box[row][column].is_singleton()
        and jacobian_box[row][column].lower == exact[row][column]
        for row in range(6)
        for column in range(6)
    )


def sgbl_enclose_adm_source_jacobian(
    box: SGBLIntervalBox,
) -> SGBLIntervalHealthRecord:
    """Enclose ``J(Z)`` and, when ``rho<1``, a Neumann inverse.

    The exact center Jacobian is the existing seven-evaluation ADM map.
    Interval columns are the first tangents of the same residual with respect
    to the six ADM accelerations, with lower jets boxed.  That is not the
    metric-``dtt`` Jacobian.
    """

    if type(box) is not SGBLIntervalBox:
        raise TypeError("box must be SGBLIntervalBox")
    point = box.center
    coefficients = sgbl_source_coefficients(point)
    if coefficients.jacobian_determinant == 0:
        raise SGBLSourceJacobianSolveStop(coefficients)
    center_inverse = inverse(coefficients.jacobian)
    limits = box.limits
    evaluations = 0
    observed_bits = 0
    zero_acc = tuple(Q(0) for _ in range(6))
    columns: list[tuple[Interval, ...]] = []
    primals: list[tuple[Interval, ...]] = []
    try:
        for field in ADM_FIELD_ORDER:
            evaluations += 1
            if evaluations > limits.max_residual_evaluations:
                raise SGBLIntervalInconclusive(
                    "resource_limit",
                    "ADM Jacobian enclosure exceeded the residual-evaluation budget",
                )
            parts = tuple(
                primal_and_tangent(value)
                for value in _full_residual(
                    point,
                    _interval_state(
                        point,
                        half_width=box.parameter_half_width,
                        accelerations=zero_acc,
                        seed_field=field,
                    ),
                )
            )
            primal = tuple(item[0] for item in parts)
            tangent = tuple(item[1] for item in parts)
            observed_bits = max(
                observed_bits,
                _enforce_bits(primal + tangent, limit=limits.max_rational_bit_length),
            )
            primals.append(primal)
            columns.append(tangent)
    except SGBLIntervalInconclusive:
        raise
    except (ValueError, ZeroDivisionError, TypeError) as exc:
        _raise_if_chart_or_reciprocal(exc)
        raise
    if any(item != primals[0] for item in primals):
        raise ValueError("interval tangent primal residual changed across ADM seeds")
    jacobian_box = tuple(
        tuple(columns[column][row] for column in range(6))
        for row in range(6)
    )
    deviation = tuple(
        tuple(
            jacobian_box[row][column]
            - Interval.singleton(coefficients.jacobian[row][column])
            for column in range(6)
        )
        for row in range(6)
    )
    remainder = interval_matrix_infinity_row_sum_bound(deviation)
    singleton = box.is_singleton and _singleton_equals_exact(
        jacobian_box, coefficients.jacobian
    )
    inverse_payload: Mapping[str, Any] | None = None
    rho: Fraction | None = None
    tail: Fraction | None = None
    invertible = False
    reason: str | None = None
    classification = "interval_inconclusive"
    try:
        inverse_payload = center_preconditioned_neumann_inverse(
            jacobian_box, center_inverse
        )
        rho = inverse_payload["rho_infinity"]
        tail = inverse_payload["inverse_tail_entrywise_bound"]
        invertible = bool(inverse_payload["rho_strictly_below_one"])
        classification = (
            "exact_singleton_inverse" if singleton else "neumann_inverse_enclosure"
        )
    except ValueError as exc:
        if "rho < 1" not in str(exc) and "does not contract" not in str(exc):
            raise
        reason = "neumann_rho_not_below_one"
        classification = "interval_inconclusive"
        invertible = False
        inverse_payload = None
    unique_root = False
    if invertible and box.acceleration_half_width > 0:
        unique_root = _unique_acceleration_root(
            primals[0],
            jacobian_box,
            center_inverse,
            box.acceleration_half_width,
        )
    return SGBLIntervalHealthRecord(
        box=box,
        branch="SGB-L",
        action=_action_mapping(point),
        reference=flat_spherical_annulus_reference(
            radial_domain_minimum=REFERENCE_RADIAL_MINIMUM,
        ),
        tilde_normal_factor=Q(TILDE_NORMAL_FACTOR),
        hat_normal_factor=Q(HAT_NORMAL_FACTOR),
        exact_center_jacobian=coefficients.jacobian,
        exact_center_det=coefficients.jacobian_determinant,
        exact_center_inverse=center_inverse,
        exact_center_jacobian_sha256=_matrix_sha256(coefficients.jacobian),
        jacobian_box=jacobian_box,
        remainder_bound=remainder,
        rho_infinity=rho,
        inverse_tail_entrywise_bound=tail,
        source_inverse=inverse_payload,
        source_invertible_over_box=invertible,
        unique_acceleration_root_over_box=unique_root,
        independent_exact_singleton_reduction=singleton,
        residual_evaluations=evaluations,
        max_rational_bit_length_observed=observed_bits,
        classification=classification,
        inconclusive_reason=reason,
        spherical_ref1_blocks=None,
        principal_cone_margins=None,
        constraint_center_valid=None,
    )


def _unique_acceleration_root(
    residual_at_zero: Sequence[Interval],
    jacobian_box: Sequence[Sequence[Interval]],
    center_inverse: Sequence[Sequence[Fraction]],
    half_width: Fraction,
) -> bool:
    displacement = tuple(
        interval(-half_width, half_width) for _ in range(6)
    )
    try:
        payload = parametric_krawczyk_inclusion(
            center_inverse=center_inverse,
            residual_at_center=residual_at_zero,
            jacobian_box=jacobian_box,
            displacement_box=displacement,
        )
    except ValueError:
        return False
    return bool(payload["krawczyk_image_strictly_inside_displacement_box"])


def jacobian_contains_point(
    jacobian_box: Sequence[Sequence[Interval]],
    point: Sequence[Sequence[Fraction]],
) -> bool:
    """Return whether an exact Jacobian lies in an interval enclosure."""

    return all(
        jacobian_box[row][column].lower
        <= point[row][column]
        <= jacobian_box[row][column].upper
        for row in range(6)
        for column in range(6)
    )


__all__ = [
    "ADM_SOURCE_CHART",
    "DEFAULT_MAX_PARAMETER_HALF_WIDTH",
    "DEFAULT_MAX_RATIONAL_BIT_LENGTH",
    "DEFAULT_MAX_RESIDUAL_EVALUATIONS",
    "LOWER_JET_BOX_SLOTS",
    "SGBLIntervalBox",
    "SGBLIntervalHealthRecord",
    "SGBLIntervalInconclusive",
    "SGBLIntervalLimits",
    "jacobian_contains_point",
    "sgbl_enclose_adm_source_jacobian",
]
