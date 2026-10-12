"""Fail-closed whole-cell SGB-L source-acceleration admission (prospective).

The existing parametric owner boxes every non-acceleration jet slot by one
uniform half-width.  The matched-family feeder therefore only certifies a
declared witness and names whole-cell covering as a missing owner.  This
module accepts an *explicit nonuniform* lower-jet product box with exact
slot provenance, evaluates interval ``R0(Z)`` and ``J(Z)`` on those slots,
and runs a predeclared acceleration box through strict parametric
Krawczyk / interval-Newton inclusion.

One six-acceleration root for every point of ``Z`` is proved only when
that inclusion closes.  Nonuniform radii are never replaced by one uniform
half-width.  A singleton witness is never promoted to a covering theorem.
Aggregate ``SGBL_branch_owned_and_healthy``, ``FRZ1``, and ``PREF1`` stay
false.  This owner does not edit the existing feeder or uniform admission.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from numbers import Integral
from types import MappingProxyType
from typing import Any, Mapping, Sequence

from .exact_interval import Interval, interval, interval_matrix_vector
from .exact_interval_krawczyk import parametric_krawczyk_inclusion
from .exact_interval_linear_algebra import center_preconditioned_neumann_inverse
from .exact_linear_algebra import inverse
from .interval_tangent import IntervalFirstTangent, primal_and_tangent
from .sgb1_ctl1_family_principal import (
    FAMILY_ACCELERATION_HALF_WIDTH,
    SGBLFamilyPrincipalStop,
    sgbl_family_adm_lower_jets,
    _source_inputs_from_witness,
    _witness,
)
from .sgb1_ctl1_initial_health import (
    SGBLContinuousCompactnessRecord,
    SGBLExactInitialSlice,
    SGBLODECellEnclosure,
)
from .sgb1_ctl1_interval_health import (
    ADM_SOURCE_CHART,
    LOWER_JET_BOX_SLOTS,
    SGBLIntervalInconclusive,
    SGBLIntervalLimits,
    _full_residual,
)
from .sgb1_ctl1_source import (
    SGBLSourceInputs,
    SGBLSourceJacobianSolveStop,
    SOURCE_ACCELERATION_ORDER,
    SOURCE_EQUATION_ORDER,
    sgbl_source_coefficients,
    sgbl_source_solve,
)
from .spherical_reduction import ADM_FIELD_ORDER, Jet2, _state_from_adm_pg_jets


Q = Fraction
DECLARED_DYADIC_DENOMINATOR = 2**64
_ZERO6 = (Q(0), Q(0), Q(0), Q(0), Q(0), Q(0))
ADM_JET_ATTR = {
    "alpha": "alpha",
    "shift": "shift",
    "lambda": "radial_metric",
    "areal_radius": "areal_radius",
    "phi": "phi",
    "chi": "chi",
}
CANONICAL_SLOT_NAMES = tuple(
    f"{field}.{component}"
    for field in ADM_FIELD_ORDER
    for component in LOWER_JET_BOX_SLOTS
)
ADMISSION_CLASSIFICATIONS = frozenset(
    {
        "parametric_krawczyk_unique_root",
        "interval_inconclusive",
        "family_geometry_incomplete",
        "unauthenticated_cell",
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
        "omitted_slot",
        "uniformization_alias",
        "correlated_uniformization",
        "promoted_singleton_witness",
        "family_geometry_incomplete",
        "unauthenticated_cell",
        "slot_center_outside_interval",
        "tampered_slot_map",
        "float_enclosure_refused",
        "domain_error",
    }
)
WRAPPING_OBSTRUCTIONS = frozenset(
    {
        "neumann_rho_not_below_one",
        "krawczyk_contraction_not_below_one",
        "krawczyk_image_not_strictly_inside",
        "residual_enclosure_misses_origin",
    }
)
COVERAGE_KINDS = frozenset(
    {
        "explicit_product",
        "flat_product",
        "exact_singleton",
        "family_cell_product",
    }
)
FORBIDDEN_HEALTH_IMPORTS = (
    "WeakCouplingThresholds",
    "weak_coupling_health_certificate",
    "esf_reference_cluster_certificate",
    "canonical_health_monitor_values",
    "background_from_spherical_state",
    "newton_residual_limit",
    "solve_accelerations",
    "FGCQRActionParameters",
    "solve_initial_data",
)
_CHART_FAILURE_MARKERS = (
    "two-dimensional base must be Lorentzian",
    "areal_radius must be positive",
    "singular metric",
    "singular jet reciprocal",
    "ADM/PG alpha and lambda must be positive",
    "effective Planck coefficient F must be positive",
)
FAMILY_SLOT_PROVENANCE: Mapping[str, tuple[str, str]] = MappingProxyType(
    {
        "alpha.value": (
            "unit-lapse polar-areal gauge",
            "alpha = 1",
        ),
        "alpha.dt": (
            "unit-lapse polar-areal gauge",
            "alpha_t = 0",
        ),
        "alpha.dr": (
            "unit-lapse polar-areal gauge",
            "alpha_r = 0",
        ),
        "alpha.dtr": (
            "unit-lapse polar-areal gauge",
            "alpha_tr = 0",
        ),
        "alpha.drr": (
            "unit-lapse polar-areal gauge",
            "alpha_rr = 0",
        ),
        "shift.value": (
            "zero-shift plus C^t=C^r=0",
            "shift = 0",
        ),
        "shift.dt": (
            "zero-shift plus C^t=C^r=0",
            "shift_t = (2 L^3-2 L+L_r r)/(4 L^3 r)",
        ),
        "shift.dr": (
            "zero-shift plus C^t=C^r=0",
            "shift_r = 0",
        ),
        "shift.dtr": (
            "zero-shift plus C^t=C^r=0",
            "shift_tr from interval d/dr of shift_t",
        ),
        "shift.drr": (
            "zero-shift plus C^t=C^r=0",
            "shift_rr = 0 on the initial slice",
        ),
        "lambda.value": (
            "SGBLODECellEnclosure.lambda_box",
            "Picard graph of the affine constraint ODE",
        ),
        "lambda.dt": (
            "polar-areal K^r_r=-2k",
            "lambda_t = 2 L k",
        ),
        "lambda.dr": (
            "affine H = H0 + H_L lambda_r",
            "lambda_r = -H0 / H_L",
        ),
        "lambda.dtr": (
            "polar-areal K^r_r=-2k",
            "lambda_tr = 2(L_r k + L k_r)",
        ),
        "lambda.drr": (
            "interval differentiation of the affine RHS",
            "d/dr(-H0/H_L) with IntervalFirstTangent; needs phi_rrr",
        ),
        "areal_radius.value": (
            "SGBLODECellEnclosure.radius",
            "authenticated radial interval; polar-areal R=r",
        ),
        "areal_radius.dt": (
            "polar-areal identity R=r",
            "R_t = -r k",
        ),
        "areal_radius.dr": (
            "polar-areal identity R=r",
            "R_r = 1",
        ),
        "areal_radius.dtr": (
            "polar-areal identity R=r",
            "R_tr = -k - r k_r",
        ),
        "areal_radius.drr": (
            "polar-areal identity R=r",
            "R_rr = 0",
        ),
        "phi.value": (
            "sgbl_compact_bump_enclosure",
            "A_phi B",
        ),
        "phi.dt": (
            "declared compact family",
            "phi_pi = 0 on the unit-lapse zero-shift slice",
        ),
        "phi.dr": (
            "sgbl_compact_bump_enclosure",
            "A_phi B_r",
        ),
        "phi.dtr": (
            "declared compact family",
            "phi_pi_r = 0",
        ),
        "phi.drr": (
            "sgbl_compact_bump_enclosure",
            "A_phi B_rr",
        ),
        "chi.value": (
            "sgbl_complete_interval_family_fields",
            "A_chi B / r",
        ),
        "chi.dt": (
            "sgbl_complete_interval_family_fields",
            "chi_t = A_chi B_r / r",
        ),
        "chi.dr": (
            "sgbl_complete_interval_family_fields",
            "A_chi (B_r/r - B/r^2)",
        ),
        "chi.dtr": (
            "sgbl_complete_interval_family_fields",
            "chi_tr = A_chi (B_rr/r - B_r/r^2)",
        ),
        "chi.drr": (
            "sgbl_complete_interval_family_fields",
            "A_chi (B_rr/r - 2 B_r/r^2 + 2 B/r^3)",
        ),
    }
)
EXPLICIT_SLOT_OWNER = "caller-declared nonuniform product slot"
EXPLICIT_SLOT_FORMULA = (
    "exact interval for this jet component only; never a shared uniform "
    "half-width and never a hidden zero"
)
ACCELERATION_SLOT_OWNER = "predeclared acceleration box"
ACCELERATION_SLOT_FORMULA = (
    "Krawczyk/interval-Newton image of R(a;z)=R0(z)+J(z)a; never a hidden "
    "zero, never a singleton-witness promotion, never a fitted residual floor"
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


def _six_fraction(name: str, value: object) -> tuple[Fraction, ...]:
    if not isinstance(value, (tuple, list)) or len(value) != 6:
        raise TypeError(f"{name} must be a 6-tuple of Fraction or built-in int")
    return tuple(_fraction(f"{name}[{index}]", item) for index, item in enumerate(value))


def _as_interval(value: object) -> Interval:
    if type(value) is Interval:
        return value
    if type(value) is IntervalFirstTangent:
        return value.primal
    if type(value) in (int, Fraction):
        return Interval.singleton(Fraction(value))
    if isinstance(value, float):
        raise SGBLCellAdmissionInconclusive(
            "float_enclosure_refused",
            "whole-cell product boxes refuse binary64-as-enclosure inputs",
        )
    raise TypeError("product-box slots must be exact intervals or rationals")


def _contains_rational(box: Interval, value: Fraction) -> bool:
    return box.lower <= value <= box.upper


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


def _interval_bit_length(value: Interval) -> int:
    return max(
        abs(value.lower.numerator).bit_length(),
        value.lower.denominator.bit_length(),
        abs(value.upper.numerator).bit_length(),
        value.upper.denominator.bit_length(),
    )


def _action_mapping(point: SGBLSourceInputs) -> dict[str, Fraction]:
    return {
        "planck_mass": point.planck_mass,
        "scalar_mass": point.scalar_mass,
        "quartic_coupling": point.quartic_coupling,
        "alpha_gb": point.alpha_gb,
        "beta": Q(0),
        "eta": Q(0),
    }


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


def _point_component(point: SGBLSourceInputs, slot_name: str) -> Fraction:
    field, component = slot_name.split(".")
    return getattr(getattr(point, ADM_JET_ATTR[field]), component)


class SGBLCellAdmissionInconclusive(ArithmeticError):
    """Typed stop when a declared product box cannot support uniqueness.

    This is not ``singular_source_jacobian``, not a PROTO4 residual floor,
    and not a scientific nonpass of the holdout.
    """

    def __init__(
        self,
        reason: str,
        message: str,
        payload: Mapping[str, Any] | None = None,
    ) -> None:
        if reason not in INCONCLUSIVE_REASONS:
            raise ValueError("unknown SGB-L whole-cell admission inconclusive reason")
        self.reason = reason
        self.payload = dict(payload or {})
        super().__init__(message)


@dataclass(frozen=True, slots=True, kw_only=True)
class SGBLCellSlot:
    """One lower-jet product-box slot with exact provenance."""

    name: str
    interval: Interval
    owner: str
    formula: str

    def __post_init__(self) -> None:
        if self.name not in CANONICAL_SLOT_NAMES:
            raise ValueError(f"unknown whole-cell product slot {self.name}")
        if type(self.interval) is not Interval:
            raise TypeError("slot interval must be an exact Interval")
        if not self.owner or not self.formula:
            raise ValueError("every product-box slot must retain owner and formula")
        if self.owner in {"substituted_zero", "binary64_placeholder"}:
            raise ValueError("slot provenance must not be a hidden zero or binary64")

    @property
    def radius(self) -> Fraction:
        return self.interval.radius()

    @property
    def is_singleton(self) -> bool:
        return self.interval.is_singleton()


def sgbl_cell_admission_slot_table() -> tuple[Mapping[str, str], ...]:
    """Return the static slot-ownership table, including acceleration dtt."""

    rows = [
        {
            "slot": name,
            "owner": FAMILY_SLOT_PROVENANCE[name][0],
            "formula": FAMILY_SLOT_PROVENANCE[name][1],
        }
        for name in CANONICAL_SLOT_NAMES
    ]
    rows.append(
        {
            "slot": "coordinate_radius",
            "owner": "authenticated cell radius or exact source-point radius",
            "formula": (
                "MHG reference connection is evaluated at the declared exact "
                "centre; polar-areal R is the areal_radius.value product slot"
            ),
        }
    )
    rows.append(
        {
            "slot": "alpha_tt, shift_tt, lambda_tt, R_tt, phi_tt, chi_tt",
            "owner": ACCELERATION_SLOT_OWNER,
            "formula": ACCELERATION_SLOT_FORMULA,
        }
    )
    return tuple(MappingProxyType(dict(row)) for row in rows)


def sgbl_declared_acceleration_box(
    center: Sequence[Fraction | int],
    half_width: Fraction | int,
) -> tuple[Interval, ...]:
    """Return a predeclared interior acceleration box.  It is not fitted."""

    values = _six_fraction("acceleration_center", center)
    width = _fraction("acceleration_half_width", half_width)
    if width <= 0:
        raise SGBLCellAdmissionInconclusive(
            "acceleration_box_not_interior",
            "whole-cell uniqueness requires a declared acceleration box with nonempty interior",
            {"acceleration_half_width": width},
        )
    return tuple(interval(value - width, value + width) for value in values)


def _canonicalize_slots(slots: Sequence[SGBLCellSlot]) -> tuple[SGBLCellSlot, ...]:
    if not isinstance(slots, (tuple, list)):
        raise TypeError("product-box slots must be a sequence of SGBLCellSlot")
    converted = tuple(slot for slot in slots)
    if any(type(slot) is not SGBLCellSlot for slot in converted):
        raise TypeError("product-box slots must be SGBLCellSlot records")
    names = tuple(slot.name for slot in converted)
    missing = tuple(name for name in CANONICAL_SLOT_NAMES if name not in names)
    extra = tuple(name for name in names if name not in CANONICAL_SLOT_NAMES)
    if missing:
        raise SGBLCellAdmissionInconclusive(
            "omitted_slot",
            "whole-cell product box omitted required lower-jet slots",
            {"omitted": missing, "extra": extra},
        )
    if extra or len(names) != len(CANONICAL_SLOT_NAMES) or len(set(names)) != len(names):
        raise ValueError("product-box slots must be the 30 canonical lower-jet names without aliases")
    by_name = {slot.name: slot for slot in converted}
    return tuple(by_name[name] for name in CANONICAL_SLOT_NAMES)


def _slots_from_jets(
    jets: Mapping[str, Jet2],
    *,
    provenance: Mapping[str, tuple[str, str]],
) -> tuple[SGBLCellSlot, ...]:
    slots: list[SGBLCellSlot] = []
    for name in CANONICAL_SLOT_NAMES:
        field, component = name.split(".")
        owner, formula = provenance[name]
        slots.append(
            SGBLCellSlot(
                name=name,
                interval=_as_interval(getattr(jets[field], component)),
                owner=owner,
                formula=formula,
            )
        )
    return tuple(slots)


def _explicit_provenance() -> dict[str, tuple[str, str]]:
    return {name: (EXPLICIT_SLOT_OWNER, EXPLICIT_SLOT_FORMULA) for name in CANONICAL_SLOT_NAMES}


def _jets_from_point(
    point: SGBLSourceInputs,
    slot_half_widths: Mapping[str, Fraction],
) -> dict[str, Jet2]:
    fields: dict[str, Jet2] = {}
    for adm_name in ADM_FIELD_ORDER:
        source = getattr(point, ADM_JET_ATTR[adm_name])
        entries: dict[str, object] = {}
        for component in LOWER_JET_BOX_SLOTS:
            name = f"{adm_name}.{component}"
            width = slot_half_widths[name]
            value = getattr(source, component)
            entries[component] = interval(value - width, value + width)
        entries["dtt"] = Interval.singleton(0)
        fields[adm_name] = Jet2(**entries)
    return fields


def _family_jets(cell: SGBLODECellEnclosure, spec: SGBLExactInitialSlice) -> dict[str, Jet2]:
    return sgbl_family_adm_lower_jets(
        radius=cell.radius,
        radial_metric=cell.lambda_box,
        angular_extrinsic_curvature=cell.k_box,
        spec=spec,
    )


def _family_slots(cell: SGBLODECellEnclosure, spec: SGBLExactInitialSlice) -> tuple[SGBLCellSlot, ...]:
    return _slots_from_jets(_family_jets(cell, spec), provenance=FAMILY_SLOT_PROVENANCE)


def _cell_authenticated(cell: SGBLODECellEnclosure) -> bool:
    return bool(
        cell.picard_strict_self_map
        and cell.residual_contains_origin
        and cell.compactness_invariant_contains_zero
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
    return all(_cell_authenticated(cell) for cell in cells)


@dataclass(frozen=True, slots=True, kw_only=True)
class SGBLCellProductBox:
    """Declared nonuniform ADM lower-jet product box with slot provenance.

    There is no ``parameter_half_width``.  Each of the thirty lower-jet
    slots carries its own interval.  A family cell box is reconstructed
    from the authenticated Picard/profile owners and is refused if those
    radii are later replaced by one shared half-width.
    """

    center: SGBLSourceInputs
    slots: tuple[SGBLCellSlot, ...]
    acceleration_box: tuple[Interval, ...]
    coverage: str
    limits: SGBLIntervalLimits = SGBLIntervalLimits()
    allow_zero_coupling_control: bool = False
    chart: str = ADM_SOURCE_CHART
    family_cell: SGBLODECellEnclosure | None = None
    family_spec: SGBLExactInitialSlice | None = None
    radius_box: Interval | None = None

    def __post_init__(self) -> None:
        if type(self.center) is not SGBLSourceInputs:
            raise TypeError("center must be SGBLSourceInputs")
        if self.coverage not in COVERAGE_KINDS:
            raise ValueError("unknown whole-cell product-box coverage")
        if type(self.limits) is not SGBLIntervalLimits:
            raise TypeError("limits must be SGBLIntervalLimits")
        if type(self.allow_zero_coupling_control) is not bool:
            raise TypeError("allow_zero_coupling_control must be boolean")
        if self.chart != ADM_SOURCE_CHART:
            raise SGBLCellAdmissionInconclusive(
                "interval_chart_error",
                "SGB-L whole-cell source Jacobian is d(MHG rows)/d(ADM accelerations)",
            )
        if self.center.alpha_gb == 0 and not self.allow_zero_coupling_control:
            raise ValueError(
                "production SGB-L product box requires nonzero alpha_gb; "
                "zero coupling is only a named control"
            )
        object.__setattr__(self, "slots", _canonicalize_slots(self.slots))
        acceleration = tuple(entry for entry in self.acceleration_box)
        if len(acceleration) != 6 or any(type(entry) is not Interval for entry in acceleration):
            raise ValueError("acceleration box must be six exact intervals")
        object.__setattr__(self, "acceleration_box", acceleration)
        if self.radius_box is not None and type(self.radius_box) is not Interval:
            raise TypeError("radius_box must be an Interval or None")
        for slot in self.slots:
            value = _point_component(self.center, slot.name)
            if not _contains_rational(slot.interval, value):
                raise SGBLCellAdmissionInconclusive(
                    "slot_center_outside_interval",
                    "declared source centre lies outside a product-box slot",
                    {"slot": slot.name},
                )
        if self.coverage == "family_cell_product":
            self._validate_family_product()
        elif self.family_cell is not None or self.family_spec is not None:
            raise ValueError("family cell/spec provenance is only for family_cell_product coverage")
        if self.coverage == "exact_singleton" and not self.is_singleton:
            raise ValueError("exact_singleton coverage requires every lower-jet slot to be a singleton")

    def _validate_family_product(self) -> None:
        if type(self.family_cell) is not SGBLODECellEnclosure:
            raise TypeError("family_cell_product coverage requires SGBLODECellEnclosure")
        if type(self.family_spec) is not SGBLExactInitialSlice:
            raise TypeError("family_cell_product coverage requires SGBLExactInitialSlice")
        if not _cell_authenticated(self.family_cell):
            raise SGBLCellAdmissionInconclusive(
                "unauthenticated_cell",
                "whole-cell family admission requires a validated Picard cell",
            )
        try:
            expected = _family_slots(self.family_cell, self.family_spec)
        except SGBLFamilyPrincipalStop as exc:
            raise SGBLCellAdmissionInconclusive(
                "domain_error" if exc.reason == "domain_error" else "zero_in_interval_reciprocal",
                str(exc),
                exc.payload,
            ) from exc
        if self.is_singleton and not self.family_cell.radius.is_singleton():
            raise SGBLCellAdmissionInconclusive(
                "promoted_singleton_witness",
                "a non-singleton family cell cannot be admitted as a singleton witness box",
            )
        if self.slots != expected:
            expected_radii = tuple(slot.radius for slot in expected)
            observed_radii = tuple(slot.radius for slot in self.slots)
            positive = tuple(radius for radius in observed_radii if radius > 0)
            correlated = (
                bool(positive) and len(set(positive)) == 1 and expected_radii != observed_radii
            )
            reason = "correlated_uniformization" if correlated else "tampered_slot_map"
            raise SGBLCellAdmissionInconclusive(
                reason,
                "family product-box slots must retain the reconstructed nonuniform radii",
                {
                    "expected_radii": expected_radii,
                    "observed_radii": observed_radii,
                },
            )
        if self.parameter_radii_are_uniform and not self.is_singleton:
            raise SGBLCellAdmissionInconclusive(
                "uniformization_alias",
                "family cell product boxes do not replace nonuniform radii by one uniform half-width",
            )
        if self.radius_box != self.family_cell.radius:
            raise SGBLCellAdmissionInconclusive(
                "tampered_slot_map",
                "family product-box radius provenance must equal the authenticated cell radius",
            )

    def slot_map(self) -> tuple[Mapping[str, object], ...]:
        """Return the exact per-slot interval and provenance."""

        rows: list[Mapping[str, object]] = []
        for slot in self.slots:
            rows.append(
                MappingProxyType(
                    {
                        "slot": slot.name,
                        "owner": slot.owner,
                        "formula": slot.formula,
                        "lower": slot.interval.lower,
                        "upper": slot.interval.upper,
                        "radius": slot.radius,
                        "singleton": slot.is_singleton,
                    }
                )
            )
        return tuple(rows)

    def slot_interval(self, field: str, component: str) -> Interval:
        name = f"{field}.{component}"
        for slot in self.slots:
            if slot.name == name:
                return slot.interval
        raise KeyError(name)

    @property
    def is_singleton(self) -> bool:
        return all(slot.is_singleton for slot in self.slots)

    @property
    def parameter_radii_are_uniform(self) -> bool:
        radii = {slot.radius for slot in self.slots}
        return len(radii) == 1

    @property
    def max_slot_radius(self) -> Fraction:
        return max(slot.radius for slot in self.slots)

    @property
    def coordinate_radius(self) -> Fraction:
        return self.center.coordinate_radius


def sgbl_cell_product_box_from_point(
    point: SGBLSourceInputs,
    *,
    slot_half_widths: Mapping[str, Fraction | int] | None = None,
    acceleration_half_width: Fraction | int,
    acceleration_center: Sequence[Fraction | int] | None = None,
    acceleration_box: Sequence[Interval] | None = None,
    limits: SGBLIntervalLimits | None = None,
    coverage: str = "explicit_product",
    allow_zero_coupling_control: bool = False,
    uniform_half_width: Fraction | int | None = None,
) -> SGBLCellProductBox:
    """Build an explicit product box.  A single uniform width is refused."""

    if type(point) is not SGBLSourceInputs:
        raise TypeError("point must be SGBLSourceInputs")
    if uniform_half_width is not None:
        raise SGBLCellAdmissionInconclusive(
            "uniformization_alias",
            "whole-cell product boxes do not replace slot radii by one uniform half-width",
            {"uniform_half_width": uniform_half_width},
        )
    widths = {name: Q(0) for name in CANONICAL_SLOT_NAMES}
    if slot_half_widths is not None:
        if not isinstance(slot_half_widths, Mapping):
            raise TypeError("slot_half_widths must be a mapping from slot name to half-width")
        unknown = tuple(name for name in slot_half_widths if name not in CANONICAL_SLOT_NAMES)
        if unknown:
            raise ValueError(f"unknown product-box slot names: {unknown}")
        for name, width in slot_half_widths.items():
            widths[name] = _fraction(f"slot_half_widths[{name}]", width, nonnegative=True)
    if acceleration_center is None:
        try:
            center_acc = sgbl_source_solve(point)
        except SGBLSourceJacobianSolveStop:
            center_acc = _ZERO6
    else:
        center_acc = _six_fraction("acceleration_center", acceleration_center)
    if acceleration_box is None:
        declared = sgbl_declared_acceleration_box(center_acc, acceleration_half_width)
    else:
        declared = tuple(_as_interval(entry) for entry in acceleration_box)
        if len(declared) != 6:
            raise ValueError("acceleration box must have six slots")
    jets = _jets_from_point(point, widths)
    coverage_kind = coverage
    if coverage_kind == "explicit_product" and all(width == 0 for width in widths.values()):
        coverage_kind = "exact_singleton"
    return SGBLCellProductBox(
        center=point,
        slots=_slots_from_jets(jets, provenance=_explicit_provenance()),
        acceleration_box=declared,
        coverage=coverage_kind if coverage == "explicit_product" else coverage,
        limits=limits or SGBLIntervalLimits(),
        allow_zero_coupling_control=allow_zero_coupling_control,
        radius_box=Interval.singleton(point.coordinate_radius),
    )


def sgbl_flat_product_box(
    *,
    slot_half_widths: Mapping[str, Fraction | int] | None = None,
    acceleration_half_width: Fraction | int = Q(1, 1 << 10),
    limits: SGBLIntervalLimits | None = None,
    alpha_gb: Fraction | int = Q(-1, 4),
) -> SGBLCellProductBox:
    """Minkowski-vacuum product box with caller-declared nonuniform radii."""

    point = SGBLSourceInputs(
        coordinate_radius=Q(5, 2),
        alpha=(1, 0, 0, 0, 0),
        shift=(0, 0, 0, 0, 0),
        radial_metric=(1, 0, 0, 0, 0),
        areal_radius=(Q(5, 2), 0, 1, 0, 0),
        phi=(0, 0, 0, 0, 0),
        chi=(0, 0, 0, 0, 0),
        planck_mass=2,
        scalar_mass=3,
        quartic_coupling=Q(1, 2),
        alpha_gb=alpha_gb,
    )
    declared_widths = dict(slot_half_widths or {})
    if not declared_widths:
        declared_widths = {
            "lambda.value": Q(1, 1 << 22),
            "areal_radius.value": Q(1, 1 << 24),
        }
    return sgbl_cell_product_box_from_point(
        point,
        slot_half_widths=declared_widths,
        acceleration_half_width=acceleration_half_width,
        limits=limits,
        coverage="flat_product",
    )


def sgbl_cell_product_box_from_family_cell(
    cell: SGBLODECellEnclosure,
    spec: SGBLExactInitialSlice,
    *,
    acceleration_half_width: Fraction | int = FAMILY_ACCELERATION_HALF_WIDTH,
    limits: SGBLIntervalLimits | None = None,
) -> SGBLCellProductBox:
    """Build the nonuniform family lower-jet product box of one Picard cell."""

    if type(cell) is not SGBLODECellEnclosure:
        raise TypeError("cell must be SGBLODECellEnclosure")
    if type(spec) is not SGBLExactInitialSlice:
        raise TypeError("spec must be SGBLExactInitialSlice")
    if not _cell_authenticated(cell):
        raise SGBLCellAdmissionInconclusive(
            "unauthenticated_cell",
            "whole-cell family admission requires a validated Picard cell",
        )
    try:
        slots = _family_slots(cell, spec)
        point, _affine = _source_inputs_from_witness(
            radius=_witness(cell.radius),
            radial_metric=_witness(cell.lambda_box),
            angular_extrinsic_curvature=_witness(cell.k_box),
            spec=spec,
        )
    except SGBLFamilyPrincipalStop as exc:
        raise SGBLCellAdmissionInconclusive(
            "domain_error" if exc.reason == "domain_error" else "zero_in_interval_reciprocal",
            str(exc),
            exc.payload,
        ) from exc
    center_acc = sgbl_source_solve(point)
    return SGBLCellProductBox(
        center=point,
        slots=slots,
        acceleration_box=sgbl_declared_acceleration_box(center_acc, acceleration_half_width),
        coverage="family_cell_product",
        limits=limits or SGBLIntervalLimits(),
        family_cell=cell,
        family_spec=spec,
        radius_box=cell.radius,
    )


@dataclass(frozen=True, slots=True, kw_only=True)
class SGBLCellAdmissionRecord:
    """Immutable whole-cell uniqueness record on one declared product box.

    Local uniqueness may be true.  Aggregate health, production admission,
    execution and holdout flags remain false.
    """

    box: SGBLCellProductBox
    acceleration_center: tuple[Fraction, ...]
    residual_at_acceleration_center: tuple[Interval, ...]
    residual_constant_box: tuple[Interval, ...]
    jacobian_box: tuple[tuple[Interval, ...], ...]
    exact_center_jacobian: tuple[tuple[Fraction, ...], ...]
    exact_center_det: Fraction
    residual_on_krawczyk_image: tuple[Interval, ...]
    krawczyk: Mapping[str, Any] | None
    interval_newton_image: tuple[Interval, ...] | None
    interval_newton_strictly_inside_acceleration_box: bool
    unique_acceleration_root_for_every_declared_parameter_point: bool
    classification: str
    inconclusive_reason: str | None
    residual_evaluations: int
    max_rational_bit_length_observed: int
    independent_exact_singleton_reduction: bool
    theorem: str
    missing_theorem: str | None

    def __post_init__(self) -> None:
        if type(self.box) is not SGBLCellProductBox:
            raise TypeError("box must be SGBLCellProductBox")
        object.__setattr__(
            self, "acceleration_center", _six_fraction("acceleration_center", self.acceleration_center)
        )
        residual_center = tuple(_as_interval(entry) for entry in self.residual_at_acceleration_center)
        residual_constant = tuple(_as_interval(entry) for entry in self.residual_constant_box)
        residual_image = tuple(_as_interval(entry) for entry in self.residual_on_krawczyk_image)
        if len(residual_center) != 6 or len(residual_constant) != 6 or len(residual_image) != 6:
            raise ValueError("residual enclosures must have six rows")
        object.__setattr__(self, "residual_at_acceleration_center", residual_center)
        object.__setattr__(self, "residual_constant_box", residual_constant)
        object.__setattr__(self, "residual_on_krawczyk_image", residual_image)
        jacobian = tuple(tuple(entry for entry in row) for row in self.jacobian_box)
        if len(jacobian) != 6 or any(
            len(row) != 6 or any(type(entry) is not Interval for entry in row) for row in jacobian
        ):
            raise ValueError("interval Jacobian must be a 6x6 Interval matrix")
        object.__setattr__(self, "jacobian_box", jacobian)
        if self.classification not in ADMISSION_CLASSIFICATIONS:
            raise ValueError("unknown whole-cell admission classification")
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
            if self.box.coverage == "family_cell_product" and self.box.is_singleton:
                if self.box.family_cell is not None and not self.box.family_cell.radius.is_singleton():
                    raise ValueError("proved family uniqueness cannot be a promoted singleton witness")
        elif self.classification == "parametric_krawczyk_unique_root" or (
            self.classification
            not in {
                "interval_inconclusive",
                "family_geometry_incomplete",
                "unauthenticated_cell",
            }
            or self.inconclusive_reason is None
        ):
            raise ValueError("failed uniqueness must be a typed incomplete or interval-inconclusive record")
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
    def FRZ1(self) -> bool:
        return False

    @property
    def PREF1(self) -> bool:
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

    @property
    def promoted_singleton_witness(self) -> bool:
        return False

    @property
    def uniformized_nonuniform_radii(self) -> bool:
        return False

    @property
    def wrapping_obstruction(self) -> bool:
        return self.inconclusive_reason in WRAPPING_OBSTRUCTIONS


@dataclass(frozen=True, slots=True, kw_only=True)
class SGBLCellAdmissionGraphRecord:
    """Aggregate whole-cell attempts on an authenticated compactness record."""

    compactness: SGBLContinuousCompactnessRecord
    cell_records: tuple[SGBLCellAdmissionRecord, ...]
    classification: str
    inconclusive_reason: str | None
    missing_theorem: str | None
    theorem: str

    def __post_init__(self) -> None:
        if type(self.compactness) is not SGBLContinuousCompactnessRecord:
            raise TypeError("compactness must be SGBLContinuousCompactnessRecord")
        object.__setattr__(self, "cell_records", tuple(self.cell_records))
        if self.classification not in {
            "whole_cell_family_records",
            "family_geometry_incomplete",
        }:
            raise ValueError("unknown whole-cell graph classification")
        if self.classification == "family_geometry_incomplete":
            if self.cell_records:
                raise ValueError("incomplete family graphs cannot fabricate whole-cell records")
            if self.inconclusive_reason != "family_geometry_incomplete":
                raise ValueError("incomplete family graphs retain family_geometry_incomplete")

    @property
    def SGBL_branch_owned_and_healthy(self) -> bool:
        return False

    @property
    def FRZ1(self) -> bool:
        return False

    @property
    def PREF1(self) -> bool:
        return False


def _product_interval_state(
    box: SGBLCellProductBox,
    accelerations: Sequence[Fraction],
    *,
    seed_field: str | None = None,
):
    fields: dict[str, Jet2] = {}
    for index, adm_name in enumerate(ADM_FIELD_ORDER):
        entries: dict[str, object] = {}
        for component in LOWER_JET_BOX_SLOTS:
            entries[component] = box.slot_interval(adm_name, component)
        dtt: object = Interval.singleton(accelerations[index])
        if seed_field == adm_name:
            dtt = IntervalFirstTangent.seed(dtt)
        entries["dtt"] = dtt
        fields[adm_name] = Jet2(**entries)
    return _state_from_adm_pg_jets(
        {
            "model_id": "SGB-L",
            "action_parameters": _action_mapping(box.center),
        },
        fields,
    )


def _maybe_outward(value: Interval, *, singleton_box: bool) -> Interval:
    return value if singleton_box else _outward(value)


def _raise_if_chart_or_reciprocal(exc: BaseException) -> None:
    message = str(exc)
    if isinstance(exc, ZeroDivisionError) or "contains zero" in message:
        raise SGBLCellAdmissionInconclusive(
            "zero_in_interval_reciprocal",
            "interval reciprocal encountered a zero-containing denominator",
        ) from exc
    if any(marker in message for marker in _CHART_FAILURE_MARKERS):
        raise SGBLCellAdmissionInconclusive(
            "interval_chart_domain",
            "declared product box is not provably inside the ADM source chart",
        ) from exc


def _enforce_bits(values: Sequence[Interval], *, limit: int) -> int:
    observed = max(_interval_bit_length(entry) for entry in values)
    if observed > limit:
        raise SGBLCellAdmissionInconclusive(
            "resource_limit",
            "interval endpoints exceeded the declared rational bit cap",
            {"observed": observed, "limit": limit},
        )
    return observed


def _enclose_residual(
    box: SGBLCellProductBox,
    accelerations: Sequence[Fraction],
    *,
    evaluations_used: int,
    seed_field: str | None = None,
) -> tuple[tuple[Interval, ...], tuple[Interval, ...] | None, int, int]:
    used = evaluations_used + 1
    if used > box.limits.max_residual_evaluations:
        raise SGBLCellAdmissionInconclusive(
            "resource_limit",
            "whole-cell residual enclosure exceeded the residual-evaluation budget",
            {"evaluations": used, "limit": box.limits.max_residual_evaluations},
        )
    try:
        parts = tuple(
            primal_and_tangent(value)
            for value in _full_residual(
                box.center,
                _product_interval_state(box, accelerations, seed_field=seed_field),
            )
        )
    except SGBLIntervalInconclusive as exc:
        reason = exc.reason if exc.reason in INCONCLUSIVE_REASONS else "interval_chart_domain"
        raise SGBLCellAdmissionInconclusive(reason, str(exc), exc.payload) from exc
    except (ValueError, ZeroDivisionError, TypeError) as exc:
        _raise_if_chart_or_reciprocal(exc)
        raise
    primal = tuple(_maybe_outward(item[0], singleton_box=box.is_singleton) for item in parts)
    tangent = None
    if seed_field is not None:
        tangent = tuple(_maybe_outward(item[1], singleton_box=box.is_singleton) for item in parts)
    observed = _enforce_bits(
        primal if tangent is None else primal + tangent,
        limit=box.limits.max_rational_bit_length,
    )
    return primal, tangent, used, observed


def _krawczyk_reason(exc: ValueError) -> str:
    message = str(exc)
    if "contraction" in message or "below one" in message:
        return "krawczyk_contraction_not_below_one"
    if "strictly inside" in message or "not strictly" in message:
        return "krawczyk_image_not_strictly_inside"
    if "interior" in message:
        return "acceleration_box_not_interior"
    raise exc


def _displacement_box(
    acceleration_box: Sequence[Interval],
    center: Sequence[Fraction],
) -> tuple[Interval, ...]:
    displacement = tuple(
        interval(box.lower - centre, box.upper - centre)
        for box, centre in zip(acceleration_box, center, strict=True)
    )
    if any(not (entry.lower < 0 < entry.upper) for entry in displacement):
        raise SGBLCellAdmissionInconclusive(
            "acceleration_box_not_interior",
            "every declared acceleration interval must contain the centre strictly in its interior",
        )
    return displacement


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


def sgbl_whole_cell_source_admission(
    box: SGBLCellProductBox,
    *,
    acceleration_center: Sequence[Fraction | int] | None = None,
) -> SGBLCellAdmissionRecord:
    """Prove or refuse a unique acceleration root on a nonuniform product box.

    A failed inclusion does not enlarge, shrink, uniformize, or otherwise
    fit the declared slots or acceleration box.  A singleton witness is
    not a covering theorem.
    """

    if type(box) is not SGBLCellProductBox:
        raise TypeError("box must be SGBLCellProductBox")
    coefficients = sgbl_source_coefficients(box.center)
    if coefficients.jacobian_determinant == 0:
        raise SGBLSourceJacobianSolveStop(coefficients)
    center_inverse = inverse(coefficients.jacobian)
    if acceleration_center is None:
        center = sgbl_source_solve(box.center)
    else:
        center = _six_fraction("acceleration_center", acceleration_center)
    displacement = _displacement_box(box.acceleration_box, center)
    evaluations = 0
    observed_bits = 0
    columns: list[tuple[Interval, ...]] = []
    primals: list[tuple[Interval, ...]] = []
    zero_acc = _ZERO6
    for field in ADM_FIELD_ORDER:
        primal, tangent, evaluations, bits = _enclose_residual(
            box,
            zero_acc,
            evaluations_used=evaluations,
            seed_field=field,
        )
        observed_bits = max(observed_bits, bits)
        if tangent is None:
            raise RuntimeError("ADM Jacobian enclosure lost its interval tangent")
        primals.append(primal)
        columns.append(tangent)
    if any(item != primals[0] for item in primals):
        raise ValueError("interval tangent primal residual changed across ADM seeds")
    jacobian_box = tuple(tuple(columns[column][row] for column in range(6)) for row in range(6))
    residual_center, _unused_tangent, evaluations, bits = _enclose_residual(
        box,
        center,
        evaluations_used=evaluations,
    )
    observed_bits = max(observed_bits, bits)
    residual_constant = _affine_image(
        residual_center,
        jacobian_box,
        tuple(-entry for entry in center),
    )
    newton_image: tuple[Interval, ...] | None = None
    newton_inside = False
    inverse_payload: Mapping[str, Any] | None = None
    reason: str | None = None
    try:
        inverse_payload = center_preconditioned_neumann_inverse(jacobian_box, center_inverse)
    except ValueError as exc:
        if "rho < 1" not in str(exc) and "does not contract" not in str(exc):
            raise
        reason = "neumann_rho_not_below_one"
    else:
        newton_image = tuple(
            -entry
            for entry in interval_matrix_vector(
                inverse_payload["inverse_enclosure"],
                residual_constant,
            )
        )
        newton_inside = all(
            newton_image[index].strictly_inside(box.acceleration_box[index]) for index in range(6)
        )
    krawczyk: Mapping[str, Any] | None = None
    unique = False
    residual_image = residual_center
    if reason is None:
        try:
            payload = parametric_krawczyk_inclusion(
                center_inverse=center_inverse,
                residual_at_center=residual_center,
                jacobian_box=jacobian_box,
                displacement_box=displacement,
            )
        except ValueError as exc:
            reason = _krawczyk_reason(exc)
        else:
            residual_image = _affine_image(
                residual_center,
                jacobian_box,
                payload["krawczyk_displacement_image_box"],
            )
            if not _vector_contains_origin(residual_image):
                reason = "residual_enclosure_misses_origin"
            else:
                krawczyk = payload
                unique = True
    classification = "parametric_krawczyk_unique_root" if unique else "interval_inconclusive"
    singleton = box.is_singleton and _singleton_equals_exact(jacobian_box, coefficients.jacobian)
    return SGBLCellAdmissionRecord(
        box=box,
        acceleration_center=center,
        residual_at_acceleration_center=residual_center,
        residual_constant_box=residual_constant,
        jacobian_box=jacobian_box,
        exact_center_jacobian=coefficients.jacobian,
        exact_center_det=coefficients.jacobian_determinant,
        residual_on_krawczyk_image=residual_image,
        krawczyk=krawczyk,
        interval_newton_image=newton_image,
        interval_newton_strictly_inside_acceleration_box=newton_inside,
        unique_acceleration_root_for_every_declared_parameter_point=unique,
        classification=classification,
        inconclusive_reason=None if unique else reason,
        residual_evaluations=evaluations,
        max_rational_bit_length_observed=observed_bits,
        independent_exact_singleton_reduction=singleton,
        theorem=(
            "For every lower-jet point z in the declared nonuniform product box, "
            "if the supplied R0 and J enclosures are valid for the affine C1 map "
            "R(a;z)=R0(z)+J(z)a, a strict infinity-norm Krawczyk contraction "
            "together with strict interior inclusion of the Krawczyk image and "
            "a residual enclosure of that image containing the origin prove a "
            "unique acceleration root in the predeclared acceleration box.  "
            "Slot radii are inputs with provenance, not a fitted uniform "
            "half-width, and a singleton witness is not a covering theorem."
        ),
        missing_theorem=None
        if unique
        else (
            None
            if reason
            in {
                "acceleration_box_not_interior",
                "resource_limit",
                "interval_chart_error",
                "interval_chart_domain",
                "zero_in_interval_reciprocal",
                "omitted_slot",
                "uniformization_alias",
                "correlated_uniformization",
                "promoted_singleton_witness",
                "tampered_slot_map",
                "unauthenticated_cell",
                "family_geometry_incomplete",
                "float_enclosure_refused",
                "slot_center_outside_interval",
                "domain_error",
            }
            else "stricter_nonuniform_product_or_acceleration_box_with_the_same_declared_inputs"
        ),
    )


def sgbl_family_cell_whole_cell_admission(
    cell: SGBLODECellEnclosure,
    spec: SGBLExactInitialSlice,
    *,
    acceleration_half_width: Fraction | int = FAMILY_ACCELERATION_HALF_WIDTH,
    limits: SGBLIntervalLimits | None = None,
) -> SGBLCellAdmissionRecord:
    """Admit one authenticated family cell, or return the exact obstruction."""

    if not _cell_authenticated(cell):
        raise SGBLCellAdmissionInconclusive(
            "unauthenticated_cell",
            "whole-cell family admission requires a validated Picard cell",
        )
    box = sgbl_cell_product_box_from_family_cell(
        cell,
        spec,
        acceleration_half_width=acceleration_half_width,
        limits=limits,
    )
    return sgbl_whole_cell_source_admission(box)


def sgbl_whole_cell_admission_from_health(
    record: SGBLContinuousCompactnessRecord,
    *,
    spec: SGBLExactInitialSlice | None = None,
    max_cells: int = 2,
    acceleration_half_width: Fraction | int = FAMILY_ACCELERATION_HALF_WIDTH,
    limits: SGBLIntervalLimits | None = None,
) -> SGBLCellAdmissionGraphRecord:
    """Admit completed family cells, or refuse an incomplete compactness graph."""

    if type(record) is not SGBLContinuousCompactnessRecord:
        raise TypeError("record must be SGBLContinuousCompactnessRecord")
    slice_spec = record.slice
    if spec is not None:
        if type(spec) is not SGBLExactInitialSlice:
            raise TypeError("spec must be SGBLExactInitialSlice")
        if spec != slice_spec:
            raise SGBLCellAdmissionInconclusive(
                "domain_error",
                "family specification does not match the authenticated health record",
            )
    max_cells = _positive_int("max_cells", max_cells)
    theorem = (
        "Whole-cell unique acceleration roots are proved only by strict "
        "parametric Krawczyk inclusion on the reconstructed nonuniform family "
        "product box.  Incomplete graphs, omitted slots, uniformized radii, "
        "and singleton-witness promotion are typed refusals.  Aggregate "
        "branch health stays false."
    )
    if not _graph_covers_support(record):
        return SGBLCellAdmissionGraphRecord(
            compactness=record,
            cell_records=(),
            classification="family_geometry_incomplete",
            inconclusive_reason="family_geometry_incomplete",
            missing_theorem=(
                "complete_validated_constraint_ODE_graph_covering_the_compact_support"
            ),
            theorem=(
                "The compactness record does not cover the compact support with a "
                "validated Picard graph.  Whole-cell source uniqueness is not "
                "fabricated from incomplete cells or from a singleton witness.  "
                f"Obstruction: {record.inconclusive_reason}."
            ),
        )
    interior = tuple(
        cell
        for cell in record.ode.cells
        if cell.radius.lower > slice_spec.support_minimum
        and cell.radius.upper < slice_spec.support_maximum
    )
    selected = interior[:max_cells] if interior else tuple(record.ode.cells[:max_cells])
    cell_records = tuple(
        sgbl_family_cell_whole_cell_admission(
            cell,
            slice_spec,
            acceleration_half_width=acceleration_half_width,
            limits=limits,
        )
        for cell in selected
    )
    return SGBLCellAdmissionGraphRecord(
        compactness=record,
        cell_records=cell_records,
        classification="whole_cell_family_records",
        inconclusive_reason=None,
        missing_theorem=None
        if all(item.unique_acceleration_root_for_every_declared_parameter_point for item in cell_records)
        else (
            cell_records[0].missing_theorem
            if cell_records
            else "stricter_nonuniform_product_or_acceleration_box_with_the_same_declared_inputs"
        ),
        theorem=theorem,
    )


def sgbl_cell_admission_health_gate(record: SGBLCellAdmissionRecord | SGBLCellAdmissionGraphRecord) -> dict[str, Any]:
    """Aggregate production flags remain closed after a local uniqueness proof."""

    if type(record) not in (SGBLCellAdmissionRecord, SGBLCellAdmissionGraphRecord):
        raise TypeError("record must be a whole-cell admission record")
    unique = False
    wrapping = False
    if type(record) is SGBLCellAdmissionRecord:
        unique = record.unique_acceleration_root_for_every_declared_parameter_point
        wrapping = record.wrapping_obstruction
        incomplete = record.classification == "family_geometry_incomplete"
    else:
        unique = bool(record.cell_records) and all(
            item.unique_acceleration_root_for_every_declared_parameter_point for item in record.cell_records
        )
        wrapping = any(item.wrapping_obstruction for item in record.cell_records)
        incomplete = record.classification == "family_geometry_incomplete"
    return {
        "SGBL_branch_owned_and_healthy": False,
        "source_solve_admission_qualified": False,
        "parametric_unique_root_certified": unique,
        "FRZ1": False,
        "PREF1": False,
        "execution_authorized": False,
        "holdout_authorized": False,
        "uses_proto4_newton_residual_limit": False,
        "finite_residual_alone_is_admission": False,
        "box_was_fitted_after_the_outcome": False,
        "promoted_singleton_witness": False,
        "uniformized_nonuniform_radii": False,
        "family_geometry_incomplete": incomplete,
        "wrapping_obstruction": wrapping,
        "missing_health_closes_the_gate": True,
    }


__all__ = [
    "ACCELERATION_SLOT_FORMULA",
    "CANONICAL_SLOT_NAMES",
    "FAMILY_SLOT_PROVENANCE",
    "FORBIDDEN_HEALTH_IMPORTS",
    "INCONCLUSIVE_REASONS",
    "SGBLCellAdmissionGraphRecord",
    "SGBLCellAdmissionInconclusive",
    "SGBLCellAdmissionRecord",
    "SGBLCellProductBox",
    "SGBLCellSlot",
    "WRAPPING_OBSTRUCTIONS",
    "sgbl_cell_admission_health_gate",
    "sgbl_cell_admission_slot_table",
    "sgbl_cell_product_box_from_family_cell",
    "sgbl_cell_product_box_from_point",
    "sgbl_declared_acceleration_box",
    "sgbl_family_cell_whole_cell_admission",
    "sgbl_flat_product_box",
    "sgbl_whole_cell_admission_from_health",
    "sgbl_whole_cell_source_admission",
]
