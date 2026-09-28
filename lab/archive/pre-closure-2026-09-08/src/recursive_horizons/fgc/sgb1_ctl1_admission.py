"""Parametric interval-Newton/Krawczyk admission for the affine SGB-L source.

At fixed lower jets ``z`` the complete six-row MHG residual is the affine map

    R(a; z) = R0(z) + J(z) a,

already owned by :mod:`sgb1_ctl1_source`.  This module takes a *declared*
lower-jet box and a *declared* acceleration box whose interior contains the
exact centre root.  It does not search for a later box, cite PROTO4
``newton_residual_limit`` or ``1e-12``, or treat a merely finite residual as
admission.

Conditional on valid interval enclosures of ``R0(Z)`` and ``J(Z)``, a strict
Krawczyk contraction ``rho_∞<1`` together with strict interior inclusion of
the Krawczyk image, and a residual enclosure of that image that contains the
origin, prove one acceleration root for every lower-jet point in the box.
That local uniqueness bit is not ``SGBL_branch_owned_and_healthy`` and not a
production ``source_solve_admission_qualified`` gate.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from types import MappingProxyType
from typing import Any, Mapping, Sequence

from .exact_interval import Interval, interval, interval_matrix_vector
from .exact_interval_krawczyk import parametric_krawczyk_inclusion
from .interval_tangent import primal_and_tangent
from .sgb1_ctl1_interval_health import (
    ADM_SOURCE_CHART,
    SGBLIntervalBox,
    SGBLIntervalHealthRecord,
    SGBLIntervalInconclusive,
    _full_residual,
    _interval_state,
    sgbl_enclose_adm_source_jacobian,
)
from .sgb1_ctl1_source import (
    SGBLSourceJacobianSolveStop,
    SOURCE_ACCELERATION_ORDER,
    SOURCE_EQUATION_ORDER,
    sgbl_source_solve,
)


Q = Fraction
_ZERO6 = (Q(0), Q(0), Q(0), Q(0), Q(0), Q(0))
ADMISSION_CLASSIFICATIONS = frozenset(
    {
        "parametric_krawczyk_unique_root",
        "interval_inconclusive",
    }
)
INCONCLUSIVE_REASONS = frozenset(
    {
        "acceleration_box_not_interior",
        "krawczyk_contraction_not_below_one",
        "krawczyk_image_not_strictly_inside",
        "residual_enclosure_misses_origin",
        "neumann_rho_not_below_one",
        "interval_chart_error",
        "interval_chart_domain",
        "zero_in_interval_reciprocal",
        "resource_limit",
    }
)


def _fraction(name: str, value: object) -> Fraction:
    if type(value) not in (int, Fraction):
        raise TypeError(f"{name} must be a Fraction or built-in int")
    return Fraction(value)


def _six_fraction(name: str, value: object) -> tuple[Fraction, ...]:
    if not isinstance(value, (tuple, list)) or len(value) != 6:
        raise TypeError(f"{name} must be a 6-tuple of Fraction or built-in int")
    return tuple(_fraction(f"{name}[{index}]", item) for index, item in enumerate(value))


def _as_interval(value: object) -> Interval:
    if type(value) is Interval:
        return value
    primal, _tangent = primal_and_tangent(value)
    if type(primal) is Interval:
        return primal
    if type(primal) in (int, Fraction):
        return Interval.singleton(Fraction(primal))
    raise TypeError("residual enclosure must be an exact interval")


def _vector_contains_origin(box: Sequence[Interval]) -> bool:
    return all(entry.contains_zero() for entry in box)


def _affine_image(
    constant: Sequence[Interval],
    jacobian: Sequence[Sequence[Interval]],
    accelerations: Sequence[object],
) -> tuple[Interval, ...]:
    return tuple(
        constant_entry + linear_entry
        for constant_entry, linear_entry in zip(
            constant,
            interval_matrix_vector(jacobian, accelerations),
            strict=True,
        )
    )


class SGBLAdmissionInconclusive(ArithmeticError):
    """Typed stop when a declared box cannot support a uniqueness certificate.

    This is not ``singular_source_jacobian``, not a PROTO4 residual floor, and
    not a scientific nonpass of the holdout.
    """

    def __init__(
        self,
        reason: str,
        message: str,
        payload: Mapping[str, Any] | None = None,
    ) -> None:
        if reason not in INCONCLUSIVE_REASONS:
            raise ValueError("unknown SGB-L admission inconclusive reason")
        self.reason = reason
        self.payload = dict(payload or {})
        super().__init__(message)


@dataclass(frozen=True, slots=True, kw_only=True)
class SGBLSourceAdmissionRecord:
    """Immutable parametric uniqueness record on one declared source box.

    Local uniqueness may be true.  Aggregate health, production admission,
    execution and holdout flags remain false or absent.
    """

    box: SGBLIntervalBox
    acceleration_center: tuple[Fraction, ...]
    acceleration_box: tuple[Interval, ...]
    jacobian_record: SGBLIntervalHealthRecord
    residual_at_acceleration_center: tuple[Interval, ...]
    residual_on_krawczyk_image: tuple[Interval, ...]
    krawczyk: Mapping[str, Any] | None
    interval_newton_image: tuple[Interval, ...] | None
    interval_newton_strictly_inside_acceleration_box: bool
    unique_acceleration_root_for_every_declared_parameter_point: bool
    classification: str
    inconclusive_reason: str | None
    residual_evaluations: int
    theorem: str
    missing_theorem: str | None

    def __post_init__(self) -> None:
        if type(self.box) is not SGBLIntervalBox:
            raise TypeError("box must be SGBLIntervalBox")
        if type(self.jacobian_record) is not SGBLIntervalHealthRecord:
            raise TypeError("jacobian_record must be SGBLIntervalHealthRecord")
        object.__setattr__(
            self, "acceleration_center", _six_fraction("acceleration_center", self.acceleration_center)
        )
        acceleration_box = tuple(entry for entry in self.acceleration_box)
        if len(acceleration_box) != 6 or any(type(entry) is not Interval for entry in acceleration_box):
            raise ValueError("acceleration box must be six exact intervals")
        object.__setattr__(self, "acceleration_box", acceleration_box)
        residual_center = tuple(_as_interval(entry) for entry in self.residual_at_acceleration_center)
        residual_image = tuple(_as_interval(entry) for entry in self.residual_on_krawczyk_image)
        if len(residual_center) != 6 or len(residual_image) != 6:
            raise ValueError("residual enclosures must have six rows")
        object.__setattr__(self, "residual_at_acceleration_center", residual_center)
        object.__setattr__(self, "residual_on_krawczyk_image", residual_image)
        if self.classification not in ADMISSION_CLASSIFICATIONS:
            raise ValueError("unknown admission classification")
        if self.unique_acceleration_root_for_every_declared_parameter_point:
            if self.classification != "parametric_krawczyk_unique_root":
                raise ValueError("proved uniqueness requires the Krawczyk classification")
            if self.krawczyk is None or self.inconclusive_reason is not None:
                raise ValueError("proved uniqueness must retain a Krawczyk payload")
            required = {
                "rho_infinity_upper_bound",
                "rho_strictly_below_one",
                "krawczyk_image_strictly_inside_displacement_box",
                "minimum_strict_componentwise_inclusion_margin",
                "conditional_theorem",
            }
            if not required <= set(self.krawczyk):
                raise ValueError("proved uniqueness must retain a Krawczyk payload")
            if (
                not self.krawczyk["rho_strictly_below_one"]
                or self.krawczyk["rho_infinity_upper_bound"] >= 1
                or not self.krawczyk["krawczyk_image_strictly_inside_displacement_box"]
                or self.krawczyk["minimum_strict_componentwise_inclusion_margin"] <= 0
            ):
                raise ValueError("proved uniqueness requires a strict Krawczyk contraction and inclusion")
            if not _vector_contains_origin(residual_image):
                raise ValueError("proved uniqueness requires a residual enclosure of the origin")
        elif self.classification != "interval_inconclusive" or self.inconclusive_reason is None:
            raise ValueError("failed uniqueness must be typed interval-inconclusive")
        if self.krawczyk is not None:
            object.__setattr__(self, "krawczyk", MappingProxyType(dict(self.krawczyk)))
        if self.interval_newton_image is not None:
            newton = tuple(_as_interval(entry) for entry in self.interval_newton_image)
            if len(newton) != 6:
                raise ValueError("interval-Newton image must have six rows")
            object.__setattr__(self, "interval_newton_image", newton)

    @property
    def acceleration_order(self) -> tuple[str, ...]:
        return SOURCE_ACCELERATION_ORDER

    @property
    def equation_order(self) -> tuple[str, ...]:
        return SOURCE_EQUATION_ORDER

    @property
    def chart(self) -> str:
        return ADM_SOURCE_CHART

    @property
    def source_solve_admission_qualified(self) -> bool:
        return False

    @property
    def SGBL_branch_owned_and_healthy(self) -> bool:
        return False

    @property
    def uses_proto4_newton_residual_limit(self) -> bool:
        return False

    @property
    def finite_residual_alone_is_admission(self) -> bool:
        return False

    @property
    def box_was_fitted_after_the_outcome(self) -> bool:
        return False


def _declared_acceleration_box(
    center: Sequence[Fraction],
    half_width: Fraction,
) -> tuple[Interval, ...]:
    if half_width <= 0:
        raise SGBLAdmissionInconclusive(
            "acceleration_box_not_interior",
            "parametric uniqueness requires a declared acceleration box with nonempty interior",
            {"acceleration_half_width": half_width},
        )
    return tuple(interval(value - half_width, value + half_width) for value in center)


def _enclose_residual_at_accelerations(
    box: SGBLIntervalBox,
    accelerations: Sequence[Fraction],
    *,
    evaluations_used: int,
) -> tuple[tuple[Interval, ...], int]:
    limits = box.limits
    used = evaluations_used + 1
    if used > limits.max_residual_evaluations:
        raise SGBLAdmissionInconclusive(
            "resource_limit",
            "parametric residual enclosure exceeded the residual-evaluation budget",
            {"evaluations": used, "limit": limits.max_residual_evaluations},
        )
    try:
        residual = tuple(
            _as_interval(entry)
            for entry in _full_residual(
                box.center,
                _interval_state(
                    box.center,
                    half_width=box.parameter_half_width,
                    accelerations=tuple(accelerations),
                ),
            )
        )
    except SGBLIntervalInconclusive as exc:
        raise SGBLAdmissionInconclusive(exc.reason, str(exc), exc.payload) from exc
    except (ValueError, ZeroDivisionError, TypeError) as exc:
        message = str(exc)
        if isinstance(exc, ZeroDivisionError) or "contains zero" in message:
            raise SGBLAdmissionInconclusive(
                "zero_in_interval_reciprocal",
                "interval reciprocal encountered a zero-containing denominator",
            ) from exc
        raise
    observed = max(
        max(
            abs(entry.lower.numerator).bit_length(),
            entry.lower.denominator.bit_length(),
            abs(entry.upper.numerator).bit_length(),
            entry.upper.denominator.bit_length(),
        )
        for entry in residual
    )
    if observed > limits.max_rational_bit_length:
        raise SGBLAdmissionInconclusive(
            "resource_limit",
            "interval endpoints exceeded the declared rational bit cap",
            {"observed": observed, "limit": limits.max_rational_bit_length},
        )
    return residual, used


def _krawczyk_reason(exc: ValueError) -> str:
    message = str(exc)
    if "contraction" in message or "below one" in message:
        return "krawczyk_contraction_not_below_one"
    if "strictly inside" in message or "not strictly" in message:
        return "krawczyk_image_not_strictly_inside"
    if "interior" in message:
        return "acceleration_box_not_interior"
    raise exc


def sgbl_parametric_source_admission(
    box: SGBLIntervalBox,
    *,
    acceleration_center: Sequence[Fraction | int] | None = None,
) -> SGBLSourceAdmissionRecord:
    """Prove or refuse a unique acceleration root on a declared source box.

    The acceleration box is the declared ``acceleration_half_width`` about
    ``acceleration_center``.  The default centre is the exact algebraic root
    of ``R0(z0)+J(z0)a=0``.  A failed inclusion does not enlarge, shrink, or
    otherwise fit that box.
    """

    if type(box) is not SGBLIntervalBox:
        raise TypeError("box must be SGBLIntervalBox")
    try:
        jacobian_record = sgbl_enclose_adm_source_jacobian(box)
    except SGBLSourceJacobianSolveStop:
        raise
    except SGBLIntervalInconclusive as exc:
        raise SGBLAdmissionInconclusive(exc.reason, str(exc), exc.payload) from exc
    if acceleration_center is None:
        center = sgbl_source_solve(box.center)
    else:
        center = _six_fraction("acceleration_center", acceleration_center)
    acceleration_box = _declared_acceleration_box(center, box.acceleration_half_width)
    displacement = tuple(
        interval(-box.acceleration_half_width, box.acceleration_half_width) for _ in range(6)
    )
    residual_center, evaluations = _enclose_residual_at_accelerations(
        box,
        center,
        evaluations_used=jacobian_record.residual_evaluations,
    )
    newton_image: tuple[Interval, ...] | None = None
    newton_inside = False
    if jacobian_record.source_inverse is not None:
        constant = _affine_image(
            residual_center,
            jacobian_record.jacobian_box,
            tuple(-entry for entry in center),
        )
        newton_image = tuple(
            -entry
            for entry in interval_matrix_vector(
                jacobian_record.source_inverse["inverse_enclosure"],
                constant,
            )
        )
        newton_inside = all(
            newton_image[index].strictly_inside(acceleration_box[index])
            for index in range(6)
        )
    krawczyk: Mapping[str, Any] | None = None
    reason: str | None = None
    unique = False
    residual_image = residual_center
    if not jacobian_record.source_invertible_over_box:
        reason = jacobian_record.inconclusive_reason or "neumann_rho_not_below_one"
    else:
        try:
            payload = parametric_krawczyk_inclusion(
                center_inverse=jacobian_record.exact_center_inverse,
                residual_at_center=residual_center,
                jacobian_box=jacobian_record.jacobian_box,
                displacement_box=displacement,
            )
        except ValueError as exc:
            reason = _krawczyk_reason(exc)
        else:
            # Affine reconstruction about the Krawczyk centre: F(a0+d)=F(a0)+J d.
            residual_image = _affine_image(
                residual_center,
                jacobian_record.jacobian_box,
                payload["krawczyk_displacement_image_box"],
            )
            if not _vector_contains_origin(residual_image):
                reason = "residual_enclosure_misses_origin"
            else:
                krawczyk = payload
                unique = True
    classification = (
        "parametric_krawczyk_unique_root" if unique else "interval_inconclusive"
    )
    return SGBLSourceAdmissionRecord(
        box=box,
        acceleration_center=center,
        acceleration_box=acceleration_box,
        jacobian_record=jacobian_record,
        residual_at_acceleration_center=residual_center,
        residual_on_krawczyk_image=residual_image,
        krawczyk=krawczyk,
        interval_newton_image=newton_image,
        interval_newton_strictly_inside_acceleration_box=newton_inside,
        unique_acceleration_root_for_every_declared_parameter_point=unique,
        classification=classification,
        inconclusive_reason=None if unique else reason,
        residual_evaluations=evaluations,
        theorem=(
            "For every lower-jet point z in the declared box, if the supplied "
            "R0 and J enclosures are valid for the affine C1 map "
            "R(a;z)=R0(z)+J(z)a, a strict infinity-norm Krawczyk contraction "
            "together with strict interior inclusion of the Krawczyk image and "
            "a residual enclosure of that image containing the origin prove a "
            "unique acceleration root in the declared acceleration box.  The "
            "declared box is an input, not a fitted output."
        ),
        missing_theorem=None if unique else (
            None
            if reason
            in {
                "acceleration_box_not_interior",
                "resource_limit",
                "interval_chart_error",
                "interval_chart_domain",
                "zero_in_interval_reciprocal",
            }
            else "stricter_parameter_or_acceleration_box_with_the_same_declared_inputs"
        ),
    )


def sgbl_source_admission_health_gate(record: SGBLSourceAdmissionRecord) -> dict[str, Any]:
    """Aggregate production flags remain closed even after a local uniqueness proof."""

    if type(record) is not SGBLSourceAdmissionRecord:
        raise TypeError("record must be SGBLSourceAdmissionRecord")
    return {
        "SGBL_branch_owned_and_healthy": False,
        "source_solve_admission_qualified": False,
        "parametric_unique_root_certified": (
            record.unique_acceleration_root_for_every_declared_parameter_point
        ),
        "FRZ1": False,
        "PREF1": False,
        "execution_authorized": False,
        "holdout_authorized": False,
        "uses_proto4_newton_residual_limit": False,
        "finite_residual_alone_is_admission": False,
        "box_was_fitted_after_the_outcome": False,
        "missing_health_closes_the_gate": True,
    }


__all__ = [
    "SGBLAdmissionInconclusive",
    "SGBLSourceAdmissionRecord",
    "sgbl_parametric_source_admission",
    "sgbl_source_admission_health_gate",
]
