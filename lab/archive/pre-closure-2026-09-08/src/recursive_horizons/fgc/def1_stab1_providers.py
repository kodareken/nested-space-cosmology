"""Outcome-blind DEF1-STAB1 provider conversions and joint-box assembly.

Each of the fourteen frozen ``Q``-error components is converted into explicit
26-input geometry radii and/or a justified additive ``Q`` debit. Geometric
sensitivities are evaluated once on the joint box
``nominal +/- sum_c radii_c``. Component contributions are then
``sum_i Lip_i * eps_{c,i} + additive_c``. Independent smaller boxes are not
used.

This slice does not read COL1 data, does not certify a global PDE error, and
does not set ``DEF1_error_map_passed``. Richardson and IMP1 remain
conditional local estimators under frozen hypotheses. Missing C-to-jet,
transport, interpolation or extraction derivative bounds refuse rather than
default to zero. A universal nonlinear PDE theorem is not a prerequisite.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from numbers import Integral
from typing import Final, Mapping, Sequence

from .def1_geometry_error import (
    INPUT_NAMES,
    QGeometryPremiseError,
    QGeometryResourceError,
    q_sensitivity_enclosures,
)
from .def1_stab1 import (
    CONDITIONAL_SOURCES,
    ERROR_BUDGET_COMPONENTS,
    PREMISE_SOURCES,
    PREMISE_STATUS_CONDITIONAL,
    PREMISE_STATUS_PROVEN,
    PREMISE_STATUSES,
    PROTOCOL_MASS_FLUX_INCONSISTENCY_TOKEN,
    PROVEN_SOURCES,
    QComponentPremise,
    SOURCE_IMP1_ADMISSION_DEBIT,
    SOURCE_RICHARDSON,
    AssembledQErrorBudget,
    Def1BooleanRecord,
    Def1Stab1Error,
    MassFluxAssessment,
    assemble_q_error_budget,
)
from .exact_interval import Interval
from .spherical_reduction import INDEPENDENT_EQUATION_ORDER


Q = Fraction
TENSOR_DIM: Final[int] = 4
FULL_RESIDUAL_LENGTH: Final[int] = len(INDEPENDENT_EQUATION_ORDER)
METRIC_ACCELERATION_INPUTS: Final[tuple[str, ...]] = (
    "alpha.dtt",
    "shift.dtt",
    "lambda.dtt",
    "R.dtt",
)
METRIC_FIRST_DERIVATIVE_INPUTS: Final[tuple[str, ...]] = tuple(
    f"{field}.{part}"
    for field in ("alpha", "shift", "lambda", "R")
    for part in ("dt", "dr")
)
TANGENT_INPUTS: Final[tuple[str, ...]] = ("k.t", "k.r")

CONVERSION_RICHARDSON: Final[str] = (
    "conditional local Richardson debit or radii; not a global PDE-error theorem"
)
CONVERSION_IMP1: Final[str] = (
    "conditional IMP1 local admission arithmetic debit; "
    "not a global PDE-error theorem"
)
CONVERSION_PHYSICAL_RESIDUAL: Final[str] = (
    "componentwise |E_ab||k^a||k^b|/F residual Q-error debit, enclosed over "
    "the joint k box; not a signed Ekk identity and not a global "
    "solution-error theorem"
)
CONVERSION_GAUGE_EXTENSION: Final[str] = (
    "componentwise complete gauge-extension |E_ab||k^a||k^b|/F residual "
    "Q-error debit, enclosed over the joint k box; not a C-to-jet inverse, "
    "not a signed Ekk identity, and not a global solution-error theorem"
)
CONVERSION_GAUGE_C_TO_JET: Final[str] = (
    "certified C-to-jet inverse bound times the full gauge vector, mapped to "
    "geometry radii; C alone is not used"
)
CONVERSION_SOURCE_INVERSION: Final[str] = (
    "certified inverse-J bound times the full residual, assigned to metric "
    "acceleration radii"
)
CONVERSION_REDUCTION: Final[str] = (
    "explicit first-derivative discrepancies times a supplied derivative/"
    "operator bound"
)
CONVERSION_SUPPLIED_RADII: Final[str] = (
    "explicit complete supplied 26-input radii; missing slots are not filled"
)
CONVERSION_AFFINE_GRONWALL: Final[str] = (
    "affine transport defect times a supplied finite-interval Gronwall factor; "
    "not residual times step by assertion"
)
CONVERSION_NULLNESS: Final[str] = (
    "null residual times a supplied null-to-tangent bound, assigned to k radii"
)
CONVERSION_INTERPOLATION_LIPSCHITZ: Final[str] = (
    "Lipschitz interpolation remainder from enclosed first derivatives times "
    "sample spacing; not a second-order remainder"
)
CONVERSION_INTERPOLATION_SECOND_ORDER: Final[str] = (
    "Lagrange second-order interpolation remainder from enclosed second "
    "derivatives"
)
CONVERSION_EXTRACTION: Final[str] = (
    "unsampled affine-gap remainder: supplied dQ/dlambda bound times max gap; "
    "positive samples are not a continuous-Q certificate"
)
CONVERSION_BOUNDARY_DEBIT: Final[str] = (
    "supplied additional boundary-to-Q debit; physical causality is not an "
    "exact zero numerical-boundary theorem"
)
CONVERSION_BOUNDARY_ZERO: Final[str] = (
    "zero additional boundary debit under explicit no-influence, coverage, "
    "and independent numerical-boundary guards"
)
CONVERSION_CONSERVATION_DEBIT: Final[str] = (
    "supplied additional conservation-to-Q debit inside a valid mass-flux ledger"
)
CONVERSION_CONSERVATION_ZERO: Final[str] = (
    "zero additional conservation debit under an explicit coverage premise "
    "and a valid mass-flux ledger"
)
CONVERSION_ARITHMETIC: Final[str] = (
    "explicit arithmetic or IMP1 local admission debit in complete-Q units"
)
CONVERSION_EXACT_ZERO: Final[str] = (
    "explicit exact-algebra zero radii and additive debit for a synthetic control"
)

REMAINING_GATE_WORK: Final[tuple[str, ...]] = (
    "spatial_temporal Richardson and IMP1 remain declared conditional local "
    "estimators under frozen hypotheses, not independently qualified "
    "continuum enclosures",
    "gauge C-to-jet inverse is not derived internally; C-only evidence refuses",
    "affine Gronwall factor is a supplied finite-interval stability premise",
    "interpolation second-order remainder requires enclosed second derivatives "
    "covering the claimed slots",
    "extraction continuous-Q remainder requires a derivative bound on "
    "unsampled affine gaps; positive samples are not that bound",
    "boundary numerical error is not implied by physical causality",
    "conservation-to-Q conversion is not a global flux theorem; mass-flux "
    "veto remains independent",
    "COL1 trajectory binding and DEF1-PREF1 classification are later; this "
    "map gate is not routed after holdout",
    "DEF1_error_map_passed remains false until independent qualification of "
    "these conversions; a universal nonlinear PDE theorem is not a prerequisite",
)

ZERO_BOUNDARY_GUARDS: Final[frozenset[str]] = frozenset(
    {
        "no_influence_premise",
        "coverage_premise",
        "independent_guard",
    }
)
ZERO_CONSERVATION_GUARDS: Final[frozenset[str]] = frozenset(
    {
        "coverage_premise",
        "mass_flux_ledger_valid",
    }
)
ALLOWED_BOUNDARY_EXTRA_GUARDS: Final[frozenset[str]] = frozenset(
    {"physical_causality_passed"}
)


class Def1Stab1ProviderError(Def1Stab1Error):
    """Fail-closed provider conversion or joint-box assembly error."""


def _rational(name: str, value: object) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, (Fraction, Integral)):
        raise TypeError(f"{name} must be an exact rational value")
    return value if isinstance(value, Fraction) else Q(int(value))


def _nonnegative_rational(name: str, value: object) -> Fraction:
    result = _rational(name, value)
    if result < 0:
        raise Def1Stab1ProviderError(f"{name} must be nonnegative")
    return result


def _positive_rational(name: str, value: object) -> Fraction:
    result = _rational(name, value)
    if result <= 0:
        raise Def1Stab1ProviderError(f"{name} must be strictly positive")
    return result


def _boolean(name: str, value: object) -> bool:
    if not isinstance(value, bool):
        raise TypeError(f"{name} must be a bool")
    return value


def _text(name: str, value: object) -> str:
    if not isinstance(value, str) or not value.strip():
        raise Def1Stab1ProviderError(f"{name} must be a nonempty string")
    return value


def _optional_kwargs(component: str, kwargs: Mapping[str, object]) -> None:
    if kwargs:
        unexpected = ", ".join(sorted(kwargs))
        raise Def1Stab1ProviderError(
            f"{component} refuses unexpected arguments ({unexpected}); "
            "H/M-only, C-only, or measured Q cannot replace a complete conversion"
        )


def _radii_from_mapping(
    values: Mapping[str, object], *, label: str
) -> tuple[tuple[str, Fraction], ...]:
    if set(values) != set(INPUT_NAMES):
        missing = [name for name in INPUT_NAMES if name not in values]
        extra = sorted(name for name in values if name not in INPUT_NAMES)
        raise Def1Stab1ProviderError(
            f"{label} must be exactly the 26 named inputs; "
            f"missing={missing} extra={extra}"
        )
    return tuple(
        (name, _nonnegative_rational(f"{label}[{name}]", values[name]))
        for name in INPUT_NAMES
    )


def _radii_pairs(value: object, *, label: str) -> tuple[tuple[str, Fraction], ...]:
    """Reject duplicate names, aliases, and out-of-order inventories."""

    if isinstance(value, Mapping):
        return _radii_from_mapping(value, label=label)
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise TypeError(f"{label} must be a mapping or a sequence of (name, radius) pairs")
    if len(value) != len(INPUT_NAMES):
        raise Def1Stab1ProviderError(
            f"{label} must list exactly the 26 geometry inputs in INPUT_NAMES order"
        )
    names: list[str] = []
    ordered: list[tuple[str, Fraction]] = []
    for index, item in enumerate(value):
        if not isinstance(item, tuple) or len(item) != 2:
            raise TypeError(f"{label}[{index}] must be a (name, radius) pair")
        name, raw = item
        if not isinstance(name, str):
            raise TypeError(f"{label}[{index}] name must be a string")
        names.append(name)
        ordered.append(
            (name, _nonnegative_rational(f"{label}[{name}]", raw))
        )
    if len(set(names)) != len(names):
        raise Def1Stab1ProviderError(
            f"{label} duplicate radius names are refused rather than collapsed"
        )
    if tuple(names) != INPUT_NAMES:
        raise Def1Stab1ProviderError(
            f"{label} must follow INPUT_NAMES order without aliases or reordering"
        )
    return tuple(ordered)


def _radii_tuple(values: Mapping[str, Fraction]) -> tuple[tuple[str, Fraction], ...]:
    return _radii_from_mapping(values, label="geometry radii")


def complete_zero_radii() -> dict[str, Fraction]:
    """Explicit zeros for slots a conversion does not occupy.

    This is a conversion output, not a missing-evidence default.
    """

    return {name: Q(0) for name in INPUT_NAMES}


def radii_from_slots(slots: Mapping[str, object]) -> dict[str, Fraction]:
    """Fill unmentioned slots with explicit conversion zeros."""

    if not isinstance(slots, Mapping):
        raise TypeError("radii slots must be a mapping")
    extra = sorted(name for name in slots if name not in INPUT_NAMES)
    if extra:
        raise Def1Stab1ProviderError(
            f"unknown geometry-radius slots {extra}; they are not filled or ignored"
        )
    values = complete_zero_radii()
    for name, raw in slots.items():
        values[name] = _nonnegative_rational(f"radii[{name}]", raw)
    return values


def require_complete_radii(values: object, *, label: str) -> dict[str, Fraction]:
    """Require every 26-input radius explicitly; missing slots refuse."""

    if not isinstance(values, Mapping):
        raise TypeError(f"{label} must be a mapping of the 26 geometry inputs")
    missing = [name for name in INPUT_NAMES if name not in values]
    extra = sorted(name for name in values if name not in INPUT_NAMES)
    if missing or extra:
        raise Def1Stab1ProviderError(
            f"{label} requires exactly the 26 named inputs; "
            f"missing={missing} extra={extra}"
        )
    return {
        name: _nonnegative_rational(f"{label}[{name}]", values[name])
        for name in INPUT_NAMES
    }


def _tensor4(name: str, value: object) -> tuple[tuple[Fraction, ...], ...]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise TypeError(f"{name} must be a 4-by-4 exact tensor")
    if len(value) != TENSOR_DIM:
        raise Def1Stab1ProviderError(
            f"{name} must be the complete 4-by-4 residual; H/M-only is insufficient"
        )
    rows: list[tuple[Fraction, ...]] = []
    for row_index, row in enumerate(value):
        if not isinstance(row, Sequence) or isinstance(row, (str, bytes)):
            raise TypeError(f"{name}[{row_index}] must be a length-4 exact row")
        if len(row) != TENSOR_DIM:
            raise Def1Stab1ProviderError(
                f"{name}[{row_index}] must have four complete components"
            )
        rows.append(
            tuple(
                _rational(f"{name}[{row_index},{column}]", entry)
                for column, entry in enumerate(row)
            )
        )
    return rows[0], rows[1], rows[2], rows[3]


def _tangent4(name: str, value: object) -> tuple[Fraction, Fraction, Fraction, Fraction]:
    """Exact four-vector for the signed ``E_ab k^a k^b`` diagnostic."""

    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise TypeError(f"{name} must be a length-2 or length-4 exact tangent")
    if len(value) not in {2, 4}:
        raise Def1Stab1ProviderError(f"{name} must be a radial (k.t, k.r) tangent")
    components = [_rational(f"{name}[{index}]", entry) for index, entry in enumerate(value)]
    if len(components) == 2:
        components.extend([Q(0), Q(0)])
    return components[0], components[1], components[2], components[3]


def _radial_future_tangent(
    name: str, value: object
) -> tuple[Fraction, Fraction]:
    """Future-directed radial tangent used by Q-error residual debits."""

    vector = _tangent4(name, value)
    if vector[2] != 0 or vector[3] != 0:
        raise Def1Stab1ProviderError(
            f"{name} must be a radial tangent; nonzero angular components refuse"
        )
    if vector[0] <= 0:
        raise Def1Stab1ProviderError(
            f"{name} must be future-directed with strictly positive k.t; "
            "geometry refuses k.t=0 and so do residual Q-error debits"
        )
    return vector[0], vector[1]


def _sequence(
    name: str, value: object, *, length: int
) -> tuple[Fraction, ...]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise TypeError(f"{name} must be an exact sequence of length {length}")
    if len(value) != length:
        raise Def1Stab1ProviderError(
            f"{name} must have length {length}; partial residuals refuse"
        )
    return tuple(_rational(f"{name}[{index}]", entry) for index, entry in enumerate(value))


def _inf_norm(values: Sequence[Fraction]) -> Fraction:
    return max((abs(item) for item in values), default=Q(0))


def residual_kk(
    residual: Sequence[Sequence[object]], tangent: Sequence[object]
) -> Fraction:
    """Exact signed ``E_ab k^a k^b``. Cancellation here is not a debit defect."""

    tensor = _tensor4("residual E_ab", residual)
    vector = _tangent4("tangent k", tangent)
    total = Q(0)
    for a in range(TENSOR_DIM):
        for b in range(TENSOR_DIM):
            total += tensor[a][b] * vector[a] * vector[b]
    return total


def conservative_residual_q_debit(
    residual: Sequence[Sequence[object]],
    tangent: Sequence[object],
    coupling_F: object,
) -> Fraction:
    """Entrywise ``sum |E_ab| |k^a| |k^b| / F`` Q-error debit at one tangent.

    This is not ``|E_ab k^a k^b|/F``: signed cancellation is refused as an
    error bound. The tangent must be future-directed and radial. The result
    is a point debit, not a whole-box enclosure.
    """

    tensor = _tensor4("residual E_ab", residual)
    k_t, k_r = _radial_future_tangent("tangent k", tangent)
    coupling = _positive_rational("coupling F", coupling_F)
    vector = (k_t, k_r, Q(0), Q(0))
    total = Q(0)
    for a in range(TENSOR_DIM):
        for b in range(TENSOR_DIM):
            total += abs(tensor[a][b]) * abs(vector[a]) * abs(vector[b])
    return total / coupling


def residual_q_debit(
    residual: Sequence[Sequence[object]],
    tangent: Sequence[object],
    coupling_F: object,
) -> Fraction:
    """Conservative pointwise Q-error debit; see ``conservative_residual_q_debit``."""

    return conservative_residual_q_debit(residual, tangent, coupling_F)


def enclose_residual_q_debit(
    residual: Sequence[Sequence[object]],
    k_t: Interval,
    k_r: Interval,
    coupling_F_lower: object,
) -> Fraction:
    """Enclose the componentwise residual debit over a joint ``k`` box."""

    if type(k_t) is not Interval or type(k_r) is not Interval:
        raise TypeError("joint k box components must be Intervals")
    if not k_t.strictly_positive():
        raise Def1Stab1ProviderError(
            "joint k.t box must be strictly future-directed; k.t=0 refuses"
        )
    coupling = _positive_rational("coupling F lower bound", coupling_F_lower)
    tensor = _tensor4("residual E_ab", residual)
    magnitudes = (k_t.abs_upper(), k_r.abs_upper(), Q(0), Q(0))
    total = Q(0)
    for a in range(TENSOR_DIM):
        for b in range(TENSOR_DIM):
            total += abs(tensor[a][b]) * magnitudes[a] * magnitudes[b]
    return total / coupling


def _source_status(component: str, status: str, source: str) -> tuple[str, str]:
    if status not in PREMISE_STATUSES:
        raise Def1Stab1ProviderError(
            f"{component} status must be proven_enclosure or conditional_premise"
        )
    if source not in PREMISE_SOURCES:
        raise Def1Stab1ProviderError(f"{component} source is not an allowed premise")
    if source in CONDITIONAL_SOURCES and status != PREMISE_STATUS_CONDITIONAL:
        raise Def1Stab1ProviderError(
            f"{component} source {source} is a conditional premise, not a proven enclosure"
        )
    if source in PROVEN_SOURCES and status != PREMISE_STATUS_PROVEN:
        raise Def1Stab1ProviderError(
            f"{component} source {source} must be labeled proven_enclosure"
        )
    if source == SOURCE_RICHARDSON and component != "spatial_temporal":
        raise Def1Stab1ProviderError(
            "Richardson may only be declared on spatial_temporal, and never "
            "as a global PDE error"
        )
    if source == SOURCE_IMP1_ADMISSION_DEBIT and component not in {
        "spatial_temporal",
        "arithmetic",
    }:
        raise Def1Stab1ProviderError(
            "IMP1 admission debit may only be declared on spatial_temporal "
            "or arithmetic, and never as a global PDE error"
        )
    return status, source


def _declaration(
    *,
    component: str,
    status: str,
    source: str,
    context: object,
    unit: object,
    provenance: object,
    conversion: str,
    guards: Sequence[str] = (),
) -> dict[str, object]:
    status, source = _source_status(component, status, source)
    if not isinstance(guards, Sequence) or isinstance(guards, (str, bytes)):
        raise TypeError(f"{component} guards must be a sequence of strings")
    frozen_guards = tuple(
        _text(f"{component}.guards[{index}]", item)
        for index, item in enumerate(guards)
    )
    if len(frozen_guards) != len(set(frozen_guards)):
        raise Def1Stab1ProviderError(
            f"{component} duplicate guards are refused"
        )
    return {
        "component": component,
        "status": status,
        "source": source,
        "context": _text(f"{component}.context", context),
        "unit": _text(f"{component}.unit", unit),
        "provenance": _text(f"{component}.provenance", provenance),
        "conversion": _text(f"{component}.conversion", conversion),
        "guards": frozen_guards,
    }


def _require_guard_membership(
    component: str,
    guards: Sequence[str],
    required: frozenset[str],
    allowed_extra: frozenset[str] = frozenset(),
) -> None:
    """Require every named guard; extras and duplicates refuse."""

    present = set(guards)
    unknown = [name for name in guards if name not in required | allowed_extra]
    if unknown:
        raise Def1Stab1ProviderError(
            f"{component} zero-debit guards contain unrelated names {unknown}; "
            "required-subset membership is not a set-proper-subset check"
        )
    missing = [name for name in sorted(required) if name not in present]
    if missing:
        if component == "boundary" and "physical_causality_passed" in present:
            raise Def1Stab1ProviderError(
                "physical causality alone does not prove exact zero "
                "numerical boundary error"
            )
        raise Def1Stab1ProviderError(
            f"{component} zero additional debit requires {sorted(required)}; "
            f"missing={missing}"
        )


@dataclass(frozen=True, slots=True)
class ProviderRecord:
    """One converted component: 26 radii, additive debit, and declarations."""

    component: str
    radii: tuple[tuple[str, Fraction], ...]
    additive_q: Fraction
    status: str
    source: str
    context: str
    unit: str
    provenance: str
    conversion: str
    guards: tuple[str, ...] = ()
    error_residual: tuple[tuple[Fraction, ...], ...] | None = None
    coupling_F_lower: Fraction | None = None
    declared_tangent: tuple[Fraction, Fraction] | None = None

    def __post_init__(self) -> None:
        if self.component not in ERROR_BUDGET_COMPONENTS:
            raise Def1Stab1ProviderError(
                f"unknown Q-error component {self.component!r}"
            )
        radii = _radii_pairs(self.radii, label=f"{self.component}.radii")
        object.__setattr__(self, "radii", radii)
        object.__setattr__(
            self,
            "additive_q",
            _nonnegative_rational(f"{self.component}.additive_q", self.additive_q),
        )
        declared = _declaration(
            component=self.component,
            status=self.status,
            source=self.source,
            context=self.context,
            unit=self.unit,
            provenance=self.provenance,
            conversion=self.conversion,
            guards=self.guards,
        )
        for key, value in declared.items():
            object.__setattr__(self, key, value)
        residual_fields = (
            self.error_residual,
            self.coupling_F_lower,
            self.declared_tangent,
        )
        if any(field is None for field in residual_fields) and any(
            field is not None for field in residual_fields
        ):
            raise Def1Stab1ProviderError(
                f"{self.component} residual debit requires residual, F lower "
                "bound, and declared tangent together"
            )
        residual_owner = self.component in {"physical_constraint", "gauge_constraint"}
        if self.error_residual is not None and not residual_owner:
            raise Def1Stab1ProviderError(
                f"{self.component} cannot relabel a physical/gauge residual debit"
            )
        physical_conversion = self.conversion == CONVERSION_PHYSICAL_RESIDUAL
        gauge_conversions = {
            CONVERSION_GAUGE_EXTENSION,
            CONVERSION_GAUGE_EXTENSION + "; plus " + CONVERSION_GAUGE_C_TO_JET,
        }
        gauge_extension = self.component == "gauge_constraint" and self.conversion in (
            gauge_conversions
        )
        if self.component == "physical_constraint" and not physical_conversion:
            raise Def1Stab1ProviderError(
                "physical_constraint requires the canonical complete-residual conversion"
            )
        if physical_conversion or gauge_extension:
            if self.error_residual is None:
                raise Def1Stab1ProviderError(
                    f"{self.component} complete-residual conversion omitted its bound inputs"
                )
        elif self.error_residual is not None:
            raise Def1Stab1ProviderError(
                f"{self.component} residual inputs disagree with its conversion"
            )
        if self.error_residual is not None:
            tangent = _radial_future_tangent(
                f"{self.component}.declared_tangent", self.declared_tangent
            )
            coupling = _positive_rational(
                f"{self.component}.coupling_F_lower", self.coupling_F_lower
            )
            tensor = _tensor4(f"{self.component}.error_residual", self.error_residual)
            object.__setattr__(self, "error_residual", tensor)
            object.__setattr__(self, "coupling_F_lower", coupling)
            object.__setattr__(self, "declared_tangent", tangent)
            expected = conservative_residual_q_debit(tensor, tangent, coupling)
            if self.additive_q != expected:
                raise Def1Stab1ProviderError(
                    f"{self.component} additive_q must equal the componentwise "
                    "point residual debit; signed Ekk cancellation is not an "
                    "error bound"
                )
        if self.source == SOURCE_RICHARDSON and self.conversion != CONVERSION_RICHARDSON:
            raise Def1Stab1ProviderError(
                "Richardson conversion must remain the canonical conditional label"
            )
        if self.source == SOURCE_IMP1_ADMISSION_DEBIT and self.conversion != CONVERSION_IMP1:
            raise Def1Stab1ProviderError(
                "IMP1 conversion must remain the canonical conditional label"
            )
        if self.component == "boundary" and self.zero_contribution:
            _require_guard_membership(
                self.component,
                self.guards,
                ZERO_BOUNDARY_GUARDS,
                ALLOWED_BOUNDARY_EXTRA_GUARDS,
            )
        if self.component == "conservation" and self.zero_contribution:
            _require_guard_membership(
                self.component,
                self.guards,
                ZERO_CONSERVATION_GUARDS,
            )

    @property
    def radii_mapping(self) -> dict[str, Fraction]:
        return {name: value for name, value in self.radii}

    @property
    def zero_contribution(self) -> bool:
        return self.additive_q == 0 and all(value == 0 for _, value in self.radii)


def _make_record(
    *,
    component: str,
    radii: Mapping[str, Fraction],
    additive_q: Fraction,
    status: str,
    source: str,
    context: object,
    unit: object,
    provenance: object,
    conversion: str,
    guards: Sequence[str] = (),
    error_residual: Sequence[Sequence[object]] | None = None,
    coupling_F_lower: object | None = None,
    declared_tangent: Sequence[object] | None = None,
) -> ProviderRecord:
    return ProviderRecord(
        component=component,
        radii=_radii_pairs(radii, label=f"{component}.radii"),
        additive_q=additive_q,
        status=status,
        source=source,
        context=context,
        unit=unit,
        provenance=provenance,
        conversion=conversion,
        guards=tuple(guards),
        error_residual=(
            None
            if error_residual is None
            else _tensor4(f"{component}.error_residual", error_residual)
        ),
        coupling_F_lower=coupling_F_lower,
        declared_tangent=(
            None
            if declared_tangent is None
            else _radial_future_tangent(
                f"{component}.declared_tangent", declared_tangent
            )
        ),
    )


def spatial_temporal_provider(
    *,
    source: str,
    status: str,
    context: object,
    unit: object,
    provenance: object,
    radii: Mapping[str, object] | None = None,
    additive_q: object = 0,
    **kwargs: object,
) -> ProviderRecord:
    """Convert a declared spatial/temporal debit. Richardson/IMP1 stay conditional."""

    _optional_kwargs("spatial_temporal", kwargs)
    if source == SOURCE_RICHARDSON:
        conversion = CONVERSION_RICHARDSON
    elif source == SOURCE_IMP1_ADMISSION_DEBIT:
        conversion = CONVERSION_IMP1
    else:
        conversion = CONVERSION_EXACT_ZERO if radii is None and additive_q == 0 else (
            "declared spatial-temporal radii or additive debit; not a global "
            "PDE-error theorem"
        )
    additive = _nonnegative_rational("spatial_temporal.additive_q", additive_q)
    mapping = (
        complete_zero_radii()
        if radii is None
        else require_complete_radii(radii, label="spatial_temporal.radii")
    )
    zero_debit = additive == 0 and all(value == 0 for value in mapping.values())
    if zero_debit and source not in PROVEN_SOURCES:
        raise Def1Stab1ProviderError(
            "spatial_temporal refuses a missing conversion; Richardson/IMP1 "
            "zero-only stubs are not a continuum debit"
        )
    return _make_record(
        component="spatial_temporal",
        radii=mapping,
        additive_q=additive,
        status=status,
        source=source,
        context=context,
        unit=unit,
        provenance=provenance,
        conversion=conversion,
    )


def physical_constraint_provider(
    *,
    residual_E: Sequence[Sequence[object]],
    tangent: Sequence[object],
    coupling_F: object,
    status: str,
    source: str,
    context: object,
    unit: object,
    provenance: object,
    **kwargs: object,
) -> ProviderRecord:
    """Map the complete physical residual to a pointwise additive Q debit."""

    _optional_kwargs("physical_constraint", kwargs)
    debit = conservative_residual_q_debit(residual_E, tangent, coupling_F)
    return _make_record(
        component="physical_constraint",
        radii=complete_zero_radii(),
        additive_q=debit,
        status=status,
        source=source,
        context=context,
        unit=unit,
        provenance=provenance,
        conversion=CONVERSION_PHYSICAL_RESIDUAL,
        error_residual=residual_E,
        coupling_F_lower=coupling_F,
        declared_tangent=tangent,
    )


def gauge_constraint_provider(
    *,
    status: str,
    source: str,
    context: object,
    unit: object,
    provenance: object,
    constraint_C: Sequence[object] | None = None,
    c_to_jet_inverse_bound: object | None = None,
    extension_E: Sequence[Sequence[object]] | None = None,
    tangent: Sequence[object] | None = None,
    coupling_F: object | None = None,
    **kwargs: object,
) -> ProviderRecord:
    """Convert complete gauge-extension residual or a certified C-to-jet inverse.

    ``C`` alone refuses. H/M-only evidence is not accepted.
    """

    _optional_kwargs("gauge_constraint", kwargs)
    has_extension = extension_E is not None
    has_inverse = c_to_jet_inverse_bound is not None
    if constraint_C is not None and not has_inverse and not has_extension:
        raise Def1Stab1ProviderError(
            "gauge C-to-jet inverse or complete extension residual is missing; "
            "C-only evidence is not replaced by zero"
        )
    if not has_extension and not has_inverse:
        raise Def1Stab1ProviderError(
            "gauge_constraint requires a complete extension residual or a "
            "certified C-to-jet inverse bound"
        )
    additive = Q(0)
    radii = complete_zero_radii()
    conversion = CONVERSION_GAUGE_EXTENSION
    if has_extension:
        if tangent is None or coupling_F is None:
            raise Def1Stab1ProviderError(
                "gauge-extension residual debit requires tangent k and coupling F"
            )
        additive = conservative_residual_q_debit(extension_E, tangent, coupling_F)
        conversion = CONVERSION_GAUGE_EXTENSION
    if has_inverse:
        if constraint_C is None:
            raise Def1Stab1ProviderError(
                "C-to-jet inverse bound requires the complete gauge vector C"
            )
        if len(constraint_C) not in {2, 4}:
            raise Def1Stab1ProviderError(
                "gauge C must be the complete spherical vector; C-only scalars refuse"
            )
        bound = _nonnegative_rational(
            "c_to_jet_inverse_bound", c_to_jet_inverse_bound
        )
        vector = _sequence(
            "constraint_C",
            constraint_C,
            length=len(constraint_C),
        )
        scale = bound * _inf_norm(vector)
        if scale == 0 and _inf_norm(vector) > 0:
            raise Def1Stab1ProviderError(
                "nonzero gauge C cannot default to zero jet radii"
            )
        radii = {name: scale for name in INPUT_NAMES}
        conversion = CONVERSION_GAUGE_C_TO_JET if not has_extension else (
            CONVERSION_GAUGE_EXTENSION + "; plus " + CONVERSION_GAUGE_C_TO_JET
        )
    residual_fields: dict[str, object] = {}
    if has_extension:
        residual_fields = {
            "error_residual": extension_E,
            "coupling_F_lower": coupling_F,
            "declared_tangent": tangent,
        }
    return _make_record(
        component="gauge_constraint",
        radii=radii,
        additive_q=additive,
        status=status,
        source=source,
        context=context,
        unit=unit,
        provenance=provenance,
        conversion=conversion,
        **residual_fields,
    )


def reduction_constraint_provider(
    *,
    first_derivative_discrepancies: Mapping[str, object],
    derivative_operator_bound: object,
    status: str,
    source: str,
    context: object,
    unit: object,
    provenance: object,
    **kwargs: object,
) -> ProviderRecord:
    """Map explicit first-derivative discrepancies through an operator bound."""

    _optional_kwargs("reduction_constraint", kwargs)
    if not isinstance(first_derivative_discrepancies, Mapping):
        raise TypeError("first_derivative_discrepancies must be a mapping")
    missing = [
        name
        for name in METRIC_FIRST_DERIVATIVE_INPUTS
        if name not in first_derivative_discrepancies
    ]
    extra = sorted(
        name
        for name in first_derivative_discrepancies
        if name not in METRIC_FIRST_DERIVATIVE_INPUTS
    )
    if missing or extra:
        raise Def1Stab1ProviderError(
            "reduction_constraint requires explicit first-derivative "
            f"discrepancies on {METRIC_FIRST_DERIVATIVE_INPUTS}; "
            f"missing={missing} extra={extra}"
        )
    operator = _nonnegative_rational(
        "derivative_operator_bound", derivative_operator_bound
    )
    slots: dict[str, Fraction] = {}
    for name in METRIC_FIRST_DERIVATIVE_INPUTS:
        discrepancy = abs(
            _rational(
                f"first_derivative_discrepancies[{name}]",
                first_derivative_discrepancies[name],
            )
        )
        if discrepancy > 0 and operator == 0:
            raise Def1Stab1ProviderError(
                "nonzero reduction discrepancies cannot default to zero radii"
            )
        slots[name] = discrepancy * operator
    return _make_record(
        component="reduction_constraint",
        radii=radii_from_slots(slots),
        additive_q=Q(0),
        status=status,
        source=source,
        context=context,
        unit=unit,
        provenance=provenance,
        conversion=CONVERSION_REDUCTION,
    )


def initial_data_provider(
    *,
    radii: Mapping[str, object],
    status: str,
    source: str,
    context: object,
    unit: object,
    provenance: object,
    additive_q: object = 0,
    **kwargs: object,
) -> ProviderRecord:
    """Require explicit complete 26-input initial-data radii."""

    _optional_kwargs("initial_data", kwargs)
    return _make_record(
        component="initial_data",
        radii=require_complete_radii(radii, label="initial_data.radii"),
        additive_q=_nonnegative_rational("initial_data.additive_q", additive_q),
        status=status,
        source=source,
        context=context,
        unit=unit,
        provenance=provenance,
        conversion=CONVERSION_SUPPLIED_RADII,
    )


def nonlinear_source_provider(
    *,
    full_residual: Sequence[object],
    inverse_jacobian_bound: object,
    status: str,
    source: str,
    context: object,
    unit: object,
    provenance: object,
    **kwargs: object,
) -> ProviderRecord:
    """Map certified inverse-J bound times the full residual to metric dtt radii."""

    _optional_kwargs("nonlinear_source", kwargs)
    residual = _sequence(
        "full_residual", full_residual, length=FULL_RESIDUAL_LENGTH
    )
    bound = _nonnegative_rational("inverse_jacobian_bound", inverse_jacobian_bound)
    scale = bound * _inf_norm(residual)
    if _inf_norm(residual) > 0 and bound == 0:
        raise Def1Stab1ProviderError(
            "nonzero source residual cannot default to zero acceleration radii"
        )
    slots = {name: scale for name in METRIC_ACCELERATION_INPUTS}
    return _make_record(
        component="nonlinear_source",
        radii=radii_from_slots(slots),
        additive_q=Q(0),
        status=status,
        source=source,
        context=context,
        unit=unit,
        provenance=provenance,
        conversion=CONVERSION_SOURCE_INVERSION,
    )


def affine_provider(
    *,
    transport_defect: object,
    affine_interval: object,
    gronwall_factor: object,
    status: str,
    source: str,
    context: object,
    unit: object,
    provenance: object,
    **kwargs: object,
) -> ProviderRecord:
    """Bound tangent error by a supplied finite-interval Gronwall factor."""

    _optional_kwargs("affine", kwargs)
    defect = _nonnegative_rational("transport_defect", transport_defect)
    _positive_rational("affine_interval", affine_interval)
    factor = _nonnegative_rational("gronwall_factor", gronwall_factor)
    if defect > 0 and factor == 0:
        raise Def1Stab1ProviderError(
            "nonzero affine defect cannot default to zero tangent radii; "
            "residual times step without a Gronwall factor refuses"
        )
    scale = factor * defect
    slots = {name: scale for name in TANGENT_INPUTS}
    return _make_record(
        component="affine",
        radii=radii_from_slots(slots),
        additive_q=Q(0),
        status=status,
        source=source,
        context=context,
        unit=unit,
        provenance=provenance,
        conversion=CONVERSION_AFFINE_GRONWALL,
    )


def affine_residual_times_step_refused(
    *,
    transport_defect: object,
    affine_interval: object,
) -> None:
    """Document that residual times step is not an affine tangent bound."""

    defect = _nonnegative_rational("transport_defect", transport_defect)
    interval = _positive_rational("affine_interval", affine_interval)
    raise Def1Stab1ProviderError(
        "affine defect-to-tangent conversion refuses residual times step "
        f"({defect}*{interval}) without an explicit finite-interval Gronwall factor"
    )


def nullness_provider(
    *,
    null_residual: object,
    null_to_tangent_bound: object,
    status: str,
    source: str,
    context: object,
    unit: object,
    provenance: object,
    **kwargs: object,
) -> ProviderRecord:
    """Map ``g(k,k)`` residual through a supplied tangent conversion bound."""

    _optional_kwargs("nullness", kwargs)
    residual = abs(_rational("null_residual", null_residual))
    bound = _nonnegative_rational("null_to_tangent_bound", null_to_tangent_bound)
    if residual > 0 and bound == 0:
        raise Def1Stab1ProviderError(
            "nonzero null residual cannot default to zero tangent radii"
        )
    scale = bound * residual
    slots = {name: scale for name in TANGENT_INPUTS}
    return _make_record(
        component="nullness",
        radii=radii_from_slots(slots),
        additive_q=Q(0),
        status=status,
        source=source,
        context=context,
        unit=unit,
        provenance=provenance,
        conversion=CONVERSION_NULLNESS,
    )


def trajectory_alignment_provider(
    *,
    radii: Mapping[str, object],
    status: str,
    source: str,
    context: object,
    unit: object,
    provenance: object,
    additive_q: object = 0,
    **kwargs: object,
) -> ProviderRecord:
    """Require explicit complete alignment radii; this is not a trajectory reader."""

    _optional_kwargs("trajectory_alignment", kwargs)
    return _make_record(
        component="trajectory_alignment",
        radii=require_complete_radii(radii, label="trajectory_alignment.radii"),
        additive_q=_nonnegative_rational(
            "trajectory_alignment.additive_q", additive_q
        ),
        status=status,
        source=source,
        context=context,
        unit=unit,
        provenance=provenance,
        conversion=CONVERSION_SUPPLIED_RADII,
    )


def interpolation_provider(
    *,
    sample_spacing: object,
    first_derivative_enclosures: Mapping[str, object],
    status: str,
    source: str,
    context: object,
    unit: object,
    provenance: object,
    remainder_radii: Mapping[str, object] | None = None,
    second_derivative_enclosures: Mapping[str, object] | None = None,
    claim_second_order_remainder: bool = False,
    **kwargs: object,
) -> ProviderRecord:
    """Lipschitz interpolation from first derivatives; second-order claims refuse."""

    _optional_kwargs("interpolation", kwargs)
    spacing = _nonnegative_rational("sample_spacing", sample_spacing)
    if not isinstance(first_derivative_enclosures, Mapping):
        raise TypeError("first_derivative_enclosures must be a mapping")
    extra = sorted(
        name for name in first_derivative_enclosures if name not in INPUT_NAMES
    )
    if extra:
        raise Def1Stab1ProviderError(
            f"interpolation derivative enclosures have unknown slots {extra}"
        )
    if not first_derivative_enclosures:
        raise Def1Stab1ProviderError(
            "interpolation derivative bounds are missing; they are not replaced by zero"
        )
    if _boolean("claim_second_order_remainder", claim_second_order_remainder):
        if not isinstance(second_derivative_enclosures, Mapping):
            raise Def1Stab1ProviderError(
                "second-order interpolation remainder requires enclosed second "
                "derivatives; a Lipschitz first-derivative bound is not that remainder"
            )
        if not second_derivative_enclosures:
            raise Def1Stab1ProviderError(
                "second-order interpolation remainder refuses an empty derivative "
                "enclosure; remainder zeros are not a second-order proof"
            )
        conversion = CONVERSION_INTERPOLATION_SECOND_ORDER
        slots: dict[str, Fraction] = {}
        for name, raw in second_derivative_enclosures.items():
            if name not in INPUT_NAMES:
                raise Def1Stab1ProviderError(
                    f"unknown interpolation second-derivative slot {name}"
                )
            enclosure = _interval_enclosure(
                f"second_derivative_enclosures[{name}]", raw
            )
            slots[name] = enclosure.abs_upper() * spacing * spacing / 2
    else:
        conversion = CONVERSION_INTERPOLATION_LIPSCHITZ
        slots = {}
        for name, raw in first_derivative_enclosures.items():
            enclosure = _interval_enclosure(
                f"first_derivative_enclosures[{name}]", raw
            )
            slots[name] = enclosure.abs_upper() * spacing
    covered = set(slots)
    remainder = remainder_radii or {}
    if not isinstance(remainder, Mapping):
        raise TypeError("remainder_radii must be a mapping")
    overlap = sorted(name for name in remainder if name in covered)
    if overlap:
        raise Def1Stab1ProviderError(
            f"interpolation remainder radii overlap converted slots {overlap}"
        )
    for name, raw in remainder.items():
        if name not in INPUT_NAMES:
            raise Def1Stab1ProviderError(
                f"unknown interpolation remainder slot {name}"
            )
        slots[name] = _nonnegative_rational(f"remainder_radii[{name}]", raw)
        covered.add(name)
    missing = [name for name in INPUT_NAMES if name not in covered]
    if missing:
        raise Def1Stab1ProviderError(
            "interpolation requires explicit remainder radii for jet slots "
            f"without derivative enclosures; missing={missing}"
        )
    return _make_record(
        component="interpolation",
        radii=radii_from_slots(slots),
        additive_q=Q(0),
        status=status,
        source=source,
        context=context,
        unit=unit,
        provenance=provenance,
        conversion=conversion,
    )


def _interval_enclosure(name: str, value: object) -> Interval:
    if type(value) is Interval:
        return Interval(
            _rational(f"{name}.lower", value.lower),
            _rational(f"{name}.upper", value.upper),
        )
    number = _rational(name, value)
    return Interval.singleton(number)


def extraction_provider(
    *,
    affine_samples: Sequence[object],
    derivative_bound: object | None = None,
    status: str,
    source: str,
    context: object,
    unit: object,
    provenance: object,
    sample_q_values: Sequence[object] | None = None,
    **kwargs: object,
) -> ProviderRecord:
    """Cover unsampled affine gaps; positive samples are not a continuous certificate."""

    _optional_kwargs("extraction", kwargs)
    if not isinstance(affine_samples, Sequence) or isinstance(
        affine_samples, (str, bytes)
    ):
        raise TypeError("affine_samples must be a sequence")
    if len(affine_samples) < 2:
        raise Def1Stab1ProviderError(
            "extraction between affine samples requires at least two samples "
            "and a derivative bound on the unsampled gaps"
        )
    samples = tuple(
        _rational(f"affine_samples[{index}]", value)
        for index, value in enumerate(affine_samples)
    )
    for index in range(1, len(samples)):
        if samples[index] <= samples[index - 1]:
            raise Def1Stab1ProviderError(
                "extraction samples must be strictly increasing"
            )
    if derivative_bound is None:
        raise Def1Stab1ProviderError(
            "extraction derivative bound is missing; positive sample values "
            "are not a continuous-Q certificate"
        )
    bound = _nonnegative_rational("extraction.derivative_bound", derivative_bound)
    gaps = tuple(samples[index] - samples[index - 1] for index in range(1, len(samples)))
    max_gap = max(gaps)
    if sample_q_values is not None:
        if not isinstance(sample_q_values, Sequence) or isinstance(
            sample_q_values, (str, bytes)
        ):
            raise TypeError("sample_q_values must be a sequence")
        if len(sample_q_values) != len(samples):
            raise Def1Stab1ProviderError(
                "sample_q_values must match the affine sample list"
            )
        for index, raw in enumerate(sample_q_values):
            _rational(f"sample_q_values[{index}]", raw)
    debit = bound * max_gap
    return _make_record(
        component="extraction",
        radii=complete_zero_radii(),
        additive_q=debit,
        status=status,
        source=source,
        context=context,
        unit=unit,
        provenance=provenance,
        conversion=CONVERSION_EXTRACTION,
    )


def boundary_provider(
    *,
    additional_debit: object,
    status: str,
    source: str,
    context: object,
    unit: object,
    provenance: object,
    physical_causality_passed: object = False,
    no_influence_premise: object = False,
    coverage_premise: object = False,
    independent_guard: object = False,
    **kwargs: object,
) -> ProviderRecord:
    """Zero additional boundary debit needs coverage/no-influence plus a guard."""

    _optional_kwargs("boundary", kwargs)
    debit = _nonnegative_rational("boundary.additional_debit", additional_debit)
    causality = _boolean("physical_causality_passed", physical_causality_passed)
    no_influence = _boolean("no_influence_premise", no_influence_premise)
    coverage = _boolean("coverage_premise", coverage_premise)
    independent = _boolean("independent_guard", independent_guard)
    guards: list[str] = []
    if causality:
        guards.append("physical_causality_passed")
    if no_influence:
        guards.append("no_influence_premise")
    if coverage:
        guards.append("coverage_premise")
    if independent:
        guards.append("independent_guard")
    if debit == 0:
        if causality and not (no_influence and coverage and independent):
            raise Def1Stab1ProviderError(
                "physical causality alone does not prove exact zero numerical "
                "boundary error"
            )
        if not (no_influence and coverage and independent):
            raise Def1Stab1ProviderError(
                "zero additional boundary debit requires no-influence, coverage, "
                "and a separate numerical-boundary guard"
            )
        conversion = CONVERSION_BOUNDARY_ZERO
    else:
        conversion = CONVERSION_BOUNDARY_DEBIT
    return _make_record(
        component="boundary",
        radii=complete_zero_radii(),
        additive_q=debit,
        status=status,
        source=source,
        context=context,
        unit=unit,
        provenance=provenance,
        conversion=conversion,
        guards=tuple(guards),
    )


def conservation_provider(
    *,
    additional_debit: object,
    mass_flux: MassFluxAssessment,
    status: str,
    source: str,
    context: object,
    unit: object,
    provenance: object,
    coverage_premise: object = False,
    **kwargs: object,
) -> ProviderRecord:
    """Zero additional conservation debit needs coverage and a valid ledger."""

    _optional_kwargs("conservation", kwargs)
    if not isinstance(mass_flux, MassFluxAssessment):
        raise TypeError("mass_flux must be a MassFluxAssessment")
    if not mass_flux.ledger_valid:
        token = mass_flux.protocol_token or PROTOCOL_MASS_FLUX_INCONSISTENCY_TOKEN
        raise Def1Stab1ProviderError(
            "conservation provider refuses a mass-flux veto; "
            f"protocol_token={token}"
        )
    debit = _nonnegative_rational("conservation.additional_debit", additional_debit)
    coverage = _boolean("coverage_premise", coverage_premise)
    guards = ["mass_flux_ledger_valid"]
    if coverage:
        guards.append("coverage_premise")
    if debit == 0 and not coverage:
        raise Def1Stab1ProviderError(
            "zero additional conservation debit requires an explicit coverage premise"
        )
    conversion = (
        CONVERSION_CONSERVATION_ZERO if debit == 0 else CONVERSION_CONSERVATION_DEBIT
    )
    return _make_record(
        component="conservation",
        radii=complete_zero_radii(),
        additive_q=debit,
        status=status,
        source=source,
        context=context,
        unit=unit,
        provenance=provenance,
        conversion=conversion,
        guards=tuple(guards),
    )


def arithmetic_provider(
    *,
    additive_q: object,
    source: str,
    status: str,
    context: object,
    unit: object,
    provenance: object,
    radii: Mapping[str, object] | None = None,
    **kwargs: object,
) -> ProviderRecord:
    """IMP1 local admission arithmetic or an explicit arithmetic debit."""

    _optional_kwargs("arithmetic", kwargs)
    conversion = CONVERSION_IMP1 if source == SOURCE_IMP1_ADMISSION_DEBIT else (
        CONVERSION_ARITHMETIC
    )
    mapping = (
        complete_zero_radii()
        if radii is None
        else require_complete_radii(radii, label="arithmetic.radii")
    )
    additive = _nonnegative_rational("arithmetic.additive_q", additive_q)
    if (
        source == SOURCE_IMP1_ADMISSION_DEBIT
        and additive == 0
        and all(value == 0 for value in mapping.values())
    ):
        raise Def1Stab1ProviderError(
            "IMP1 arithmetic controls require a nonzero local admission debit, "
            "not a registered zero stub"
        )
    return _make_record(
        component="arithmetic",
        radii=mapping,
        additive_q=additive,
        status=status,
        source=source,
        context=context,
        unit=unit,
        provenance=provenance,
        conversion=conversion,
    )


@dataclass(frozen=True, slots=True)
class ComponentContribution:
    """Joint-box contribution of one frozen component."""

    component: str
    radii_term: Fraction
    additive_q: Fraction
    total: Fraction

    def __post_init__(self) -> None:
        if self.component not in ERROR_BUDGET_COMPONENTS:
            raise Def1Stab1ProviderError(
                f"unknown Q-error component {self.component!r}"
            )
        radii_term = _nonnegative_rational(
            f"{self.component}.radii_term", self.radii_term
        )
        additive = _nonnegative_rational(
            f"{self.component}.additive_q", self.additive_q
        )
        total = _nonnegative_rational(f"{self.component}.total", self.total)
        if total != radii_term + additive:
            raise Def1Stab1ProviderError(
                f"{self.component} total must equal radii_term plus additive_q"
            )
        object.__setattr__(self, "radii_term", radii_term)
        object.__setattr__(self, "additive_q", additive)
        object.__setattr__(self, "total", total)


@dataclass(frozen=True, slots=True)
class AssembledProviderMap:
    """Checked provider DTO wrapping the frozen error-budget assembly."""

    nominal: tuple[tuple[str, Fraction], ...]
    records: tuple[ProviderRecord, ...]
    contributions: tuple[ComponentContribution, ...]
    joint_radii: tuple[tuple[str, Fraction], ...]
    joint_box: tuple[tuple[str, Interval], ...]
    lipschitz: tuple[tuple[str, Interval], ...]
    assembled: AssembledQErrorBudget
    mass_flux: MassFluxAssessment
    remaining_gate_work: tuple[str, ...]
    evaluation_context: str
    unit: str
    def1_error_map_passed: bool
    used_measured_q: bool
    global_pde_error_certified: bool
    richardson_treated_as_global_pde_error: bool
    imp1_admission_debit_treated_as_global_pde_error: bool
    def1_booleans: Def1BooleanRecord | None = None

    def __post_init__(self) -> None:
        expected = _derive_assembled_map(
            nominal={name: value for name, value in self.nominal},
            records=self.records,
            mass_flux=self.mass_flux,
            def1_booleans=self.def1_booleans,
        )
        for name, value in expected.items():
            if getattr(self, name) != value:
                raise Def1Stab1ProviderError(
                    f"assembled provider map {name} does not match the joint-box "
                    "rederivation from stored nominal, records, and mass-flux; "
                    "cached Jacobian, inner budget, and flags are not trusted"
                )


def _lipschitz_term(
    lipschitz: Mapping[str, Interval], radii: Mapping[str, Fraction]
) -> Fraction:
    return sum(
        (lipschitz[name].abs_upper() * radii[name] for name in INPUT_NAMES),
        Q(0),
    )


def _nominal_point(nominal: Mapping[str, object]) -> dict[str, Fraction]:
    if not isinstance(nominal, Mapping):
        raise TypeError("nominal geometry must be a mapping")
    if set(nominal) != set(INPUT_NAMES):
        missing = [name for name in INPUT_NAMES if name not in nominal]
        extra = sorted(name for name in nominal if name not in INPUT_NAMES)
        raise Def1Stab1ProviderError(
            "nominal geometry requires exactly the 26 named inputs; "
            f"missing={missing} extra={extra}"
        )
    return {name: _rational(name, nominal[name]) for name in INPUT_NAMES}


def _nominal_pairs(nominal: Mapping[str, object]) -> tuple[tuple[str, Fraction], ...]:
    center = _nominal_point(nominal)
    return tuple((name, center[name]) for name in INPUT_NAMES)


def _shared_evaluation_identity(
    records: Sequence[ProviderRecord],
) -> tuple[str, str]:
    contexts = {record.context for record in records}
    units = {record.unit for record in records}
    if len(contexts) != 1 or len(units) != 1:
        raise Def1Stab1ProviderError(
            "provider assembly requires one shared evaluation-context and unit "
            "identity; component provenance may vary, and these strings remain "
            "declarations rather than authenticated production provenance"
        )
    context = next(iter(contexts))
    unit = next(iter(units))
    return context, unit


def _ordered_records(records: Sequence[ProviderRecord]) -> tuple[ProviderRecord, ...]:
    if not isinstance(records, Sequence) or isinstance(records, (str, bytes)):
        raise TypeError("records must be a sequence of ProviderRecord")
    seen: dict[str, ProviderRecord] = {}
    for index, record in enumerate(records):
        if not isinstance(record, ProviderRecord):
            raise TypeError(f"records[{index}] must be a ProviderRecord")
        if record.component in seen:
            raise Def1Stab1ProviderError(
                f"duplicate Q-error component ownership for {record.component}"
            )
        seen[record.component] = record
    missing = [name for name in ERROR_BUDGET_COMPONENTS if name not in seen]
    extra = sorted(name for name in seen if name not in ERROR_BUDGET_COMPONENTS)
    if missing or extra:
        raise Def1Stab1ProviderError(
            "provider records must be exactly the frozen component list; "
            f"missing={missing} extra={extra}"
        )
    return tuple(seen[name] for name in ERROR_BUDGET_COMPONENTS)


def _record_additive_on_box(
    record: ProviderRecord,
    center: Mapping[str, Fraction],
    box: Mapping[str, Interval],
) -> Fraction:
    if record.error_residual is None:
        return record.additive_q
    declared = record.declared_tangent
    if declared != (center["k.t"], center["k.r"]):
        raise Def1Stab1ProviderError(
            f"{record.component} declared tangent is foreign to the nominal "
            "joint evaluation context; a point contraction is not a whole-box bound"
        )
    return enclose_residual_q_debit(
        record.error_residual,
        box["k.t"],
        box["k.r"],
        record.coupling_F_lower,
    )


def joint_input_radii(records: Sequence[ProviderRecord]) -> dict[str, Fraction]:
    """Return the componentwise sum of all provider radii."""

    summed = complete_zero_radii()
    for record in records:
        if not isinstance(record, ProviderRecord):
            raise TypeError("records must contain ProviderRecord values")
        for name, value in record.radii:
            summed[name] += value
    return summed


def _derive_assembled_map(
    *,
    nominal: Mapping[str, object],
    records: Sequence[ProviderRecord],
    mass_flux: MassFluxAssessment,
    def1_booleans: Def1BooleanRecord | None,
) -> dict[str, object]:
    ordered = _ordered_records(records)
    context, unit = _shared_evaluation_identity(ordered)
    if not isinstance(mass_flux, MassFluxAssessment):
        raise TypeError("mass_flux must be a MassFluxAssessment")
    if not mass_flux.ledger_valid:
        token = mass_flux.protocol_token or PROTOCOL_MASS_FLUX_INCONSISTENCY_TOKEN
        raise Def1Stab1ProviderError(
            "assembled provider map refuses a mass-flux veto; "
            f"protocol_token={token}"
        )
    if def1_booleans is not None and not isinstance(def1_booleans, Def1BooleanRecord):
        raise TypeError("def1_booleans must be an independent Def1BooleanRecord")
    center = _nominal_point(nominal)
    summed = joint_input_radii(ordered)
    box = {
        name: Interval(center[name] - summed[name], center[name] + summed[name])
        for name in INPUT_NAMES
    }
    try:
        derivatives = q_sensitivity_enclosures(box)
    except QGeometryResourceError:
        raise
    except QGeometryPremiseError as error:
        raise Def1Stab1ProviderError(
            "joint geometry box is not a certified annular Lorentzian domain"
        ) from error
    contributions = []
    premises = []
    for record in ordered:
        radii_term = _lipschitz_term(derivatives, record.radii_mapping)
        additive = _record_additive_on_box(record, center, box)
        total = radii_term + additive
        contributions.append(
            ComponentContribution(
                component=record.component,
                radii_term=radii_term,
                additive_q=additive,
                total=total,
            )
        )
        premises.append(
            QComponentPremise(
                component=record.component,
                sensitivity=Interval.singleton(1),
                input_error=Interval.singleton(total),
                status=record.status,
                source=record.source,
            )
        )
    assembled = assemble_q_error_budget(premises)
    for record, (name, status), (_, source) in zip(
        ordered, assembled.component_statuses, assembled.component_sources, strict=True
    ):
        if record.component != name or record.status != status or record.source != source:
            raise Def1Stab1ProviderError(
                f"{record.component} outer source/status must match the inner "
                "AssembledQErrorBudget inventory"
            )
    if (
        assembled.used_measured_q
        or assembled.global_pde_error_certified
        or assembled.richardson_treated_as_global_pde_error
        or assembled.imp1_admission_debit_treated_as_global_pde_error
    ):
        raise Def1Stab1ProviderError(
            "assembled budget cannot certify measured Q or a global PDE error"
        )
    return {
        "nominal": _nominal_pairs(center),
        "records": ordered,
        "contributions": tuple(contributions),
        "joint_radii": _radii_tuple(summed),
        "joint_box": tuple((name, box[name]) for name in INPUT_NAMES),
        "lipschitz": tuple((name, derivatives[name]) for name in INPUT_NAMES),
        "assembled": assembled,
        "mass_flux": mass_flux,
        "remaining_gate_work": REMAINING_GATE_WORK,
        "evaluation_context": context,
        "unit": unit,
        "def1_error_map_passed": False,
        "used_measured_q": False,
        "global_pde_error_certified": False,
        "richardson_treated_as_global_pde_error": False,
        "imp1_admission_debit_treated_as_global_pde_error": False,
        "def1_booleans": def1_booleans,
    }


def _trusted_assembled_map(fields: Mapping[str, object]) -> AssembledProviderMap:
    """Factory construction: fields already derived, so skip a second Jacobian."""

    inst = object.__new__(AssembledProviderMap)
    for name in AssembledProviderMap.__dataclass_fields__:
        object.__setattr__(inst, name, fields[name])
    return inst


def assemble_provider_error_map(
    nominal: Mapping[str, object],
    records: Sequence[ProviderRecord],
    *,
    mass_flux: MassFluxAssessment,
    def1_booleans: Def1BooleanRecord | None = None,
    **kwargs: object,
) -> AssembledProviderMap:
    """Assemble the fourteen conversions on the single joint geometry box."""

    if kwargs:
        unexpected = ", ".join(sorted(kwargs))
        raise Def1Stab1ProviderError(
            "provider assembly refuses unexpected arguments "
            f"({unexpected}); measured Q cannot enter the bound"
        )
    return _trusted_assembled_map(
        _derive_assembled_map(
            nominal=nominal,
            records=records,
            mass_flux=mass_flux,
            def1_booleans=def1_booleans,
        )
    )


__all__ = [
    "AssembledProviderMap",
    "CONVERSION_AFFINE_GRONWALL",
    "CONVERSION_EXTRACTION",
    "CONVERSION_GAUGE_EXTENSION",
    "CONVERSION_IMP1",
    "CONVERSION_INTERPOLATION_LIPSCHITZ",
    "CONVERSION_PHYSICAL_RESIDUAL",
    "CONVERSION_RICHARDSON",
    "CONVERSION_SOURCE_INVERSION",
    "CONVERSION_SUPPLIED_RADII",
    "ComponentContribution",
    "Def1Stab1ProviderError",
    "conservative_residual_q_debit",
    "enclose_residual_q_debit",
    "FULL_RESIDUAL_LENGTH",
    "METRIC_ACCELERATION_INPUTS",
    "METRIC_FIRST_DERIVATIVE_INPUTS",
    "ProviderRecord",
    "REMAINING_GATE_WORK",
    "affine_provider",
    "affine_residual_times_step_refused",
    "arithmetic_provider",
    "assemble_provider_error_map",
    "boundary_provider",
    "complete_zero_radii",
    "conservation_provider",
    "extraction_provider",
    "gauge_constraint_provider",
    "initial_data_provider",
    "interpolation_provider",
    "joint_input_radii",
    "nonlinear_source_provider",
    "nullness_provider",
    "physical_constraint_provider",
    "radii_from_slots",
    "reduction_constraint_provider",
    "require_complete_radii",
    "residual_kk",
    "residual_q_debit",
    "spatial_temporal_provider",
    "trajectory_alignment_provider",
]
