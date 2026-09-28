"""Bounded outcome-blind DEF1-STAB1 core.

This module owns the frozen ``RaychaudhuriErrorBudget`` and
``ActivationAssessment`` contracts, fail-closed componentwise ``Q``-error
assembly from supplied sensitivity and input-error enclosures, exact
covariant Misner--Sharp algebra, typed ``mass_flux_inconsistency`` mapping,
activation arithmetic with an explicit case reference, and independent
margin helpers.  It does not read a trajectory, certify a global PDE error,
or emit a complete physical DEF1 gate.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from fractions import Fraction
from math import inf, isfinite, nextafter
from numbers import Integral
from typing import Final, Mapping, Sequence

from .exact_interval import Interval, coerce_interval


Q = Fraction
BASE_DIM: Final[int] = 2

ERROR_BUDGET_COMPONENTS: Final[tuple[str, ...]] = (
    "spatial_temporal",
    "physical_constraint",
    "gauge_constraint",
    "reduction_constraint",
    "initial_data",
    "nonlinear_source",
    "affine",
    "nullness",
    "trajectory_alignment",
    "interpolation",
    "extraction",
    "boundary",
    "conservation",
    "arithmetic",
)

DEF1_BOOLEAN_NAMES: Final[tuple[str, ...]] = (
    "resolved_activation",
    "control_dominance",
    "resolved_trapped_interval",
    "resolved_complete_defocusing",
    "direct_Raychaudhuri_agreement",
    "all_health_constraints_scales_valid",
    "mass_flux_ledger_valid",
    "three_resolution_convergence",
    "two_method_agreement",
)

PROTOCOL_MASS_FLUX_INCONSISTENCY_TOKEN: Final[str] = "mass_flux_inconsistency"
STOPPED_MASS_FLUX_INCONSISTENCY: Final[str] = "stopped_mass_flux_inconsistency"

GR0_PHI_POLICY_CANONICAL_ZERO: Final[str] = "canonical_gr0_phi_zero"
GR0_PHI_POLICY_HISTORICAL_PLANTED_SEED: Final[str] = (
    "historical_planted_calibration_seed"
)
GR0_PHI_POLICIES: Final[frozenset[str]] = frozenset(
    {
        GR0_PHI_POLICY_CANONICAL_ZERO,
        GR0_PHI_POLICY_HISTORICAL_PLANTED_SEED,
    }
)

# Declaration labels supplied by the caller, not independently
# authenticated enclosure certificates.
PREMISE_STATUS_PROVEN: Final[str] = "proven_enclosure"
PREMISE_STATUS_CONDITIONAL: Final[str] = "conditional_premise"
PREMISE_STATUSES: Final[frozenset[str]] = frozenset(
    {PREMISE_STATUS_PROVEN, PREMISE_STATUS_CONDITIONAL}
)

SOURCE_SUPPLIED_CERTIFIED: Final[str] = "supplied_certified_enclosure"
SOURCE_EXACT_ALGEBRA: Final[str] = "exact_algebra"
SOURCE_RICHARDSON: Final[str] = "richardson"
SOURCE_IMP1_ADMISSION_DEBIT: Final[str] = "imp1_admission_debit"
SOURCE_DECLARED_CONDITIONAL: Final[str] = "declared_conditional_enclosure"
PREMISE_SOURCES: Final[frozenset[str]] = frozenset(
    {
        SOURCE_SUPPLIED_CERTIFIED,
        SOURCE_EXACT_ALGEBRA,
        SOURCE_RICHARDSON,
        SOURCE_IMP1_ADMISSION_DEBIT,
        SOURCE_DECLARED_CONDITIONAL,
    }
)
CONDITIONAL_SOURCES: Final[frozenset[str]] = frozenset(
    {
        SOURCE_RICHARDSON,
        SOURCE_IMP1_ADMISSION_DEBIT,
        SOURCE_DECLARED_CONDITIONAL,
    }
)
PROVEN_SOURCES: Final[frozenset[str]] = frozenset(
    {SOURCE_SUPPLIED_CERTIFIED, SOURCE_EXACT_ALGEBRA}
)
RICHARDSON_OR_IMP1_SOURCES: Final[frozenset[str]] = frozenset(
    {SOURCE_RICHARDSON, SOURCE_IMP1_ADMISSION_DEBIT}
)

REQUIRED_MATCHED_CONTROL_NAMES: Final[tuple[str, ...]] = ("GR-0", "SGB-L")
MASS_FLUX_INCONSISTENCY_REASONS: Final[frozenset[str]] = frozenset(
    {
        "conservation_residual_exceeds_enclosure",
        "missing_trace_adjusted_projector",
        "hamiltonian_momentum_only",
        "nonpositive_coupling_F",
        "malformed_mass_inputs",
        "incomplete_base_tensor",
    }
)

MIN_ACTIVATION_LOWER_BOUND: Final[Fraction] = Q(2)
TRAPPEDNESS_ERROR_FACTOR: Final[int] = 4
COMPLETE_Q_ERROR_FACTOR: Final[int] = 4
MIN_AFFINE_SAMPLES: Final[int] = 8
MIN_AFFINE_INTERVAL_OVER_L0: Final[Fraction] = Q(1, 64)
DIRECT_RAYCHAUDHURI_AGREEMENT_MAX: Final[Fraction] = Q(1, 100_000_000)
AFFINE_NORMALIZATION_RESIDUAL_MAX: Final[Fraction] = Q(1, 10_000_000_000)
MIN_NESTED_RESOLUTIONS: Final[int] = 3
MIN_INDEPENDENT_METHODS: Final[int] = 2
MIN_OBSERVED_ORDER: Final[Fraction] = Q(3, 2)

Matrix2 = tuple[tuple[Fraction, Fraction], tuple[Fraction, Fraction]]
Vector2 = tuple[Fraction, Fraction]


class Def1Stab1Error(ValueError):
    """Fail-closed DEF1-STAB1 contract, algebra, or premise error."""


def _rational(name: str, value: object) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, (Fraction, Integral)):
        raise TypeError(f"{name} must be an exact rational value")
    return value if isinstance(value, Fraction) else Q(int(value))


def _nonnegative_rational(name: str, value: object) -> Fraction:
    result = _rational(name, value)
    if result < 0:
        raise Def1Stab1Error(f"{name} must be nonnegative")
    return result


def _positive_rational(name: str, value: object) -> Fraction:
    result = _rational(name, value)
    if result <= 0:
        raise Def1Stab1Error(f"{name} must be strictly positive")
    return result


def _boolean(name: str, value: object) -> bool:
    if not isinstance(value, bool):
        raise TypeError(f"{name} must be a bool")
    return value


def _interval(name: str, value: object) -> Interval:
    try:
        interval = coerce_interval(value)
    except TypeError as exc:
        raise TypeError(f"{name} must be an exact interval or rational") from exc
    except ValueError as exc:
        raise Def1Stab1Error(f"{name} is a malformed or overlapping interval") from exc
    return interval


def _nonnegative_interval(name: str, value: object) -> Interval:
    interval = _interval(name, value)
    if interval.lower < 0:
        raise Def1Stab1Error(f"{name} enclosure must be nonnegative")
    return interval


def _finite_nonnegative_float(name: str, value: object) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise TypeError(f"{name} must be a real float")
    try:
        number = float(value)
    except OverflowError as exc:
        raise Def1Stab1Error(f"{name} is outside the finite binary64 range") from exc
    if not isfinite(number) or number < 0.0:
        raise Def1Stab1Error(f"{name} must be finite and nonnegative")
    if isinstance(value, int) and Q(number) != value:
        raise Def1Stab1Error(f"{name} cannot be losslessly stored as binary64")
    return number


def _outward_upper_float(value: Fraction) -> float:
    if value < 0:
        raise Def1Stab1Error("upper-bound conversion requires a nonnegative rational")
    try:
        approx = float(value)
    except OverflowError as exc:
        raise Def1Stab1Error("upper-bound conversion left the finite range") from exc
    if not isfinite(approx) or approx < 0.0:
        raise Def1Stab1Error("upper-bound conversion left the finite range")
    if Q(approx) < value:
        approx = nextafter(approx, inf)
    if not isfinite(approx) or approx < 0.0 or Q(approx) < value:
        raise Def1Stab1Error("upper-bound conversion failed to round outward")
    return approx


def _outward_lower_float(value: Fraction) -> float:
    try:
        approx = float(value)
    except OverflowError as exc:
        raise Def1Stab1Error("lower-bound conversion left the finite range") from exc
    if not isfinite(approx):
        raise Def1Stab1Error("lower-bound conversion left the finite range")
    if Q(approx) > value:
        approx = nextafter(approx, -inf)
    if not isfinite(approx) or Q(approx) > value:
        raise Def1Stab1Error("lower-bound conversion failed to round outward")
    return approx


def _fraction_from_nonnegative_float(name: str, value: float) -> Fraction:
    if not isfinite(value) or value < 0.0:
        raise Def1Stab1Error(f"{name} must be finite and nonnegative")
    return Q(value)


def _stored_component_sum(component_floats: Sequence[float]) -> Fraction:
    total = Q(0)
    for index, value in enumerate(component_floats):
        total += _fraction_from_nonnegative_float(f"component[{index}]", value)
    return total


def _dominating_total(
    component_floats: Sequence[float], exact_total: Fraction
) -> float:
    stored_sum = _stored_component_sum(component_floats)
    required = exact_total if exact_total >= stored_sum else stored_sum
    total = _outward_upper_float(required)
    while Q(total) < stored_sum or Q(total) < exact_total:
        nxt = nextafter(total, inf)
        if not isfinite(nxt) or nxt <= total:
            raise Def1Stab1Error("error-budget total overflowed")
        total = nxt
    return total


def _finite_float(name: str, value: object) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise TypeError(f"{name} must be a real float")
    try:
        number = float(value)
    except OverflowError as exc:
        raise Def1Stab1Error(f"{name} is outside the finite binary64 range") from exc
    if not isfinite(number):
        raise Def1Stab1Error(f"{name} must be finite")
    if isinstance(value, int) and Q(number) != value:
        raise Def1Stab1Error(f"{name} cannot be losslessly stored as binary64")
    return number


def _nonnegative_integer(name: str, value: object) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral):
        raise TypeError(f"{name} must be an integer")
    number = int(value)
    if number < 0:
        raise Def1Stab1Error(f"{name} must be nonnegative")
    return number


def _optional_rational(name: str, value: object) -> Fraction | None:
    if value is None:
        return None
    return _rational(name, value)


def _optional_vector2(name: str, value: object) -> Vector2 | None:
    if value is None:
        return None
    return _vector2(name, value)


def _matrix2(name: str, value: object) -> Matrix2:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise TypeError(f"{name} must be a 2-by-2 exact matrix")
    if len(value) != BASE_DIM:
        raise Def1Stab1Error(f"{name} must be a 2-by-2 exact matrix")
    rows: list[tuple[Fraction, Fraction]] = []
    for row_index, row in enumerate(value):
        if not isinstance(row, Sequence) or isinstance(row, (str, bytes)):
            raise TypeError(f"{name}[{row_index}] must be a length-2 exact row")
        if len(row) != BASE_DIM:
            raise Def1Stab1Error(f"{name}[{row_index}] must be a length-2 exact row")
        rows.append(
            (
                _rational(f"{name}[{row_index},0]", row[0]),
                _rational(f"{name}[{row_index},1]", row[1]),
            )
        )
    return rows[0], rows[1]


def _vector2(name: str, value: object) -> Vector2:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise TypeError(f"{name} must be a length-2 exact vector")
    if len(value) != BASE_DIM:
        raise Def1Stab1Error(f"{name} must be a length-2 exact vector")
    return (
        _rational(f"{name}[0]", value[0]),
        _rational(f"{name}[1]", value[1]),
    )


def _scale_matrix(scalar: Fraction, matrix: Matrix2) -> Matrix2:
    return (
        (scalar * matrix[0][0], scalar * matrix[0][1]),
        (scalar * matrix[1][0], scalar * matrix[1][1]),
    )


def _subtract_matrix(left: Matrix2, right: Matrix2) -> Matrix2:
    return (
        (left[0][0] - right[0][0], left[0][1] - right[0][1]),
        (left[1][0] - right[1][0], left[1][1] - right[1][1]),
    )


def _trace(matrix: Matrix2) -> Fraction:
    return matrix[0][0] + matrix[1][1]


@dataclass(frozen=True, slots=True)
class RaychaudhuriErrorBudget:
    """Frozen componentwise upper bound on error in complete ``Q``."""

    spatial_temporal: float
    physical_constraint: float
    gauge_constraint: float
    reduction_constraint: float
    initial_data: float
    nonlinear_source: float
    affine: float
    nullness: float
    trajectory_alignment: float
    interpolation: float
    extraction: float
    boundary: float
    conservation: float
    arithmetic: float
    total_upper_bound: float

    def __post_init__(self) -> None:
        values = []
        for name in ERROR_BUDGET_COMPONENTS:
            values.append(_finite_nonnegative_float(name, getattr(self, name)))
        total = _finite_nonnegative_float(
            "total_upper_bound", self.total_upper_bound
        )
        stored_sum = _stored_component_sum(values)
        if Q(total) < stored_sum:
            raise Def1Stab1Error(
                "total_upper_bound must dominate the exact sum of the stored "
                "binary64 component upper bounds"
            )
        for name, value in zip(ERROR_BUDGET_COMPONENTS, values, strict=True):
            object.__setattr__(self, name, value)
        object.__setattr__(self, "total_upper_bound", total)


@dataclass(frozen=True, slots=True)
class ActivationAssessment:
    """Resolved activation factor against the fixed threshold and controls."""

    raw_factor: float
    lower_bound: float
    matched_control_upper_bound: float
    threshold_passed: bool
    control_dominance_passed: bool

    def __post_init__(self) -> None:
        raw = _finite_nonnegative_float("raw_factor", self.raw_factor)
        lower = _finite_float("lower_bound", self.lower_bound)
        control = _finite_nonnegative_float(
            "matched_control_upper_bound", self.matched_control_upper_bound
        )
        threshold = _boolean("threshold_passed", self.threshold_passed)
        dominance = _boolean(
            "control_dominance_passed", self.control_dominance_passed
        )
        stored_lower = Q(lower)
        stored_control = Q(control)
        if lower > raw:
            raise Def1Stab1Error("activation lower bound cannot exceed its raw factor")
        if threshold is not (stored_lower >= MIN_ACTIVATION_LOWER_BOUND):
            raise Def1Stab1Error(
                "threshold_passed must match the saved lower bound against 2"
            )
        if dominance is not (stored_lower > stored_control):
            raise Def1Stab1Error(
                "control_dominance_passed must be a strict comparison of the "
                "saved outward bounds; overlapping or equal saved bounds fail"
            )
        object.__setattr__(self, "raw_factor", raw)
        object.__setattr__(self, "lower_bound", lower)
        object.__setattr__(self, "matched_control_upper_bound", control)


@dataclass(frozen=True, slots=True)
class QComponentPremise:
    """One named ``Q``-error component with an explicit enclosure status."""

    component: str
    sensitivity: Interval
    input_error: Interval
    status: str
    source: str

    def __post_init__(self) -> None:
        if self.component not in ERROR_BUDGET_COMPONENTS:
            raise Def1Stab1Error(
                f"unknown Q-error component {self.component!r}"
            )
        object.__setattr__(
            self,
            "sensitivity",
            _interval(f"{self.component}.sensitivity", self.sensitivity),
        )
        object.__setattr__(
            self,
            "input_error",
            _nonnegative_interval(f"{self.component}.input_error", self.input_error),
        )
        if self.status not in PREMISE_STATUSES:
            raise Def1Stab1Error(
                f"{self.component} status must be proven_enclosure or "
                "conditional_premise"
            )
        if self.source not in PREMISE_SOURCES:
            raise Def1Stab1Error(f"{self.component} source is not an allowed premise")
        if self.source in CONDITIONAL_SOURCES and self.status != PREMISE_STATUS_CONDITIONAL:
            raise Def1Stab1Error(
                f"{self.component} source {self.source} is a conditional premise, "
                "not a proven enclosure"
            )
        if self.source in PROVEN_SOURCES and self.status != PREMISE_STATUS_PROVEN:
            raise Def1Stab1Error(
                f"{self.component} source {self.source} must be labeled "
                "proven_enclosure"
            )
        if self.source == SOURCE_RICHARDSON and self.component != "spatial_temporal":
            raise Def1Stab1Error(
                "Richardson may only be declared on spatial_temporal, and never "
                "as a global PDE error"
            )
        if self.source == SOURCE_IMP1_ADMISSION_DEBIT and self.component not in {
            "spatial_temporal",
            "arithmetic",
        }:
            raise Def1Stab1Error(
                "IMP1 admission debit may only be declared on spatial_temporal "
                "or arithmetic, and never as a global PDE error"
            )


def _ordered_component_inventory(
    name: str,
    inventory: object,
    allowed: frozenset[str],
) -> tuple[tuple[str, str], ...]:
    if not isinstance(inventory, tuple):
        raise TypeError(f"{name} must be a tuple")
    if len(inventory) != len(ERROR_BUDGET_COMPONENTS):
        raise Def1Stab1Error(f"{name} must follow the frozen component order")
    ordered: list[tuple[str, str]] = []
    for index, item in enumerate(inventory):
        if not isinstance(item, tuple) or len(item) != 2:
            raise TypeError(f"{name}[{index}] must be a (component, label) pair")
        component, label = item
        if component != ERROR_BUDGET_COMPONENTS[index]:
            raise Def1Stab1Error(f"{name} must follow the frozen component order")
        if not isinstance(label, str) or label not in allowed:
            raise Def1Stab1Error(f"{name}[{component}] is a malformed declaration")
        ordered.append((component, label))
    return tuple(ordered)


def _exact_component_inventory(
    inventory: object,
) -> tuple[tuple[str, Fraction], ...]:
    if not isinstance(inventory, tuple):
        raise TypeError("exact_components must be a tuple")
    if len(inventory) != len(ERROR_BUDGET_COMPONENTS):
        raise Def1Stab1Error("exact_components must follow the frozen component order")
    ordered: list[tuple[str, Fraction]] = []
    for index, item in enumerate(inventory):
        if not isinstance(item, tuple) or len(item) != 2:
            raise TypeError("exact_components entries must be (component, Fraction) pairs")
        component, value = item
        if component != ERROR_BUDGET_COMPONENTS[index]:
            raise Def1Stab1Error("exact_components must follow the frozen component order")
        ordered.append(
            (component, _nonnegative_rational(f"exact_components[{component}]", value))
        )
    return tuple(ordered)


@dataclass(frozen=True, slots=True)
class AssembledQErrorBudget:
    """Exact assembly record wrapping the frozen float error budget."""

    budget: RaychaudhuriErrorBudget
    exact_components: tuple[tuple[str, Fraction], ...]
    exact_total: Fraction
    component_statuses: tuple[tuple[str, str], ...]
    component_sources: tuple[tuple[str, str], ...]
    all_components_declared_proven_enclosures: bool
    global_pde_error_certified: bool
    used_measured_q: bool
    richardson_treated_as_global_pde_error: bool
    imp1_admission_debit_treated_as_global_pde_error: bool

    def __post_init__(self) -> None:
        if not isinstance(self.budget, RaychaudhuriErrorBudget):
            raise TypeError("budget must be a RaychaudhuriErrorBudget")
        exact_components = _exact_component_inventory(self.exact_components)
        object.__setattr__(self, "exact_components", exact_components)
        exact_total = _nonnegative_rational("exact_total", self.exact_total)
        object.__setattr__(self, "exact_total", exact_total)
        if exact_total != sum((value for _, value in exact_components), Q(0)):
            raise Def1Stab1Error("exact_total must equal the exact component sum")
        for name, exact in exact_components:
            stored = Q(getattr(self.budget, name))
            if stored < exact:
                raise Def1Stab1Error(
                    f"stored {name} bound drifted below its exact upper bound"
                )
        if Q(self.budget.total_upper_bound) < exact_total:
            raise Def1Stab1Error(
                "stored total_upper_bound drifted below the exact component sum"
            )
        statuses = _ordered_component_inventory(
            "component_statuses", self.component_statuses, PREMISE_STATUSES
        )
        sources = _ordered_component_inventory(
            "component_sources", self.component_sources, PREMISE_SOURCES
        )
        object.__setattr__(self, "component_statuses", statuses)
        object.__setattr__(self, "component_sources", sources)
        for (component, status), (_, source) in zip(statuses, sources, strict=True):
            if source in CONDITIONAL_SOURCES and status != PREMISE_STATUS_CONDITIONAL:
                raise Def1Stab1Error(
                    f"{component} source {source} is a conditional declaration"
                )
            if source in PROVEN_SOURCES and status != PREMISE_STATUS_PROVEN:
                raise Def1Stab1Error(
                    f"{component} source {source} must be declared proven_enclosure"
                )
            if source == SOURCE_RICHARDSON and component != "spatial_temporal":
                raise Def1Stab1Error(
                    "Richardson may only be declared on spatial_temporal"
                )
            if source == SOURCE_IMP1_ADMISSION_DEBIT and component not in {
                "spatial_temporal",
                "arithmetic",
            }:
                raise Def1Stab1Error(
                    "IMP1 admission debit may only be declared on spatial_temporal "
                    "or arithmetic"
                )
        declared_all_proven = all(
            status == PREMISE_STATUS_PROVEN for _, status in statuses
        )
        declared = _boolean(
            "all_components_declared_proven_enclosures",
            self.all_components_declared_proven_enclosures,
        )
        if declared is not declared_all_proven:
            raise Def1Stab1Error(
                "all_components_declared_proven_enclosures must match the "
                "supplied status declarations"
            )
        if _boolean("used_measured_q", self.used_measured_q):
            raise Def1Stab1Error("Q-error assembly must remain independent of measured Q")
        if _boolean(
            "global_pde_error_certified", self.global_pde_error_certified
        ) or _boolean(
            "richardson_treated_as_global_pde_error",
            self.richardson_treated_as_global_pde_error,
        ) or _boolean(
            "imp1_admission_debit_treated_as_global_pde_error",
            self.imp1_admission_debit_treated_as_global_pde_error,
        ):
            raise Def1Stab1Error(
                "Richardson and IMP1 admission debit are not global PDE error"
            )


@dataclass(frozen=True, slots=True)
class MatchedControlSample:
    """Matched-control peak and activation-error enclosure."""

    name: str
    max_abs_phi: Fraction
    activation_error: Fraction

    def __post_init__(self) -> None:
        if self.name not in REQUIRED_MATCHED_CONTROL_NAMES:
            raise Def1Stab1Error(
                "matched controls must be the required GR-0 and SGB-L pair"
            )
        object.__setattr__(
            self,
            "max_abs_phi",
            _nonnegative_rational(f"{self.name}.max_abs_phi", self.max_abs_phi),
        )
        object.__setattr__(
            self,
            "activation_error",
            _nonnegative_rational(f"{self.name}.activation_error", self.activation_error),
        )


@dataclass(frozen=True, slots=True)
class AffineIntervalAssessment:
    """Independent affine-sample and interval helper, not a trajectory reader."""

    sample_count: int
    interval_over_L0: Fraction
    samples_passed: bool
    interval_passed: bool
    passed: bool

    def __post_init__(self) -> None:
        sample_count = _nonnegative_integer("sample_count", self.sample_count)
        object.__setattr__(self, "sample_count", sample_count)
        interval = _nonnegative_rational("interval_over_L0", self.interval_over_L0)
        object.__setattr__(self, "interval_over_L0", interval)
        samples_passed = _boolean("samples_passed", self.samples_passed)
        interval_passed = _boolean("interval_passed", self.interval_passed)
        passed = _boolean("passed", self.passed)
        expected_samples = sample_count >= MIN_AFFINE_SAMPLES
        expected_interval = interval >= MIN_AFFINE_INTERVAL_OVER_L0
        if samples_passed is not expected_samples:
            raise Def1Stab1Error(
                "samples_passed must match sample_count >= 8"
            )
        if interval_passed is not expected_interval:
            raise Def1Stab1Error(
                "interval_passed must match Delta lambda / L0 >= 1/64"
            )
        if passed is not (expected_samples and expected_interval):
            raise Def1Stab1Error(
                "passed must be the conjunction of the sample and interval flags"
            )


@dataclass(frozen=True, slots=True)
class MassFluxAssessment:
    """Misner--Sharp ledger assessment with typed protocol-token mapping."""

    mass: Fraction | None
    mass_gradient: Vector2 | None
    equation_defect: Vector2 | None
    residual: Fraction
    enclosure: Fraction
    ledger_valid: bool
    protocol_token: str | None
    veto_class: str | None
    typed_reason: str | None
    used_hamiltonian_momentum_only: bool

    def __post_init__(self) -> None:
        ledger_valid = _boolean("ledger_valid", self.ledger_valid)
        if _boolean(
            "used_hamiltonian_momentum_only", self.used_hamiltonian_momentum_only
        ):
            raise Def1Stab1Error(
                "Hamiltonian/momentum residuals alone are insufficient for the "
                "full spacetime Misner-Sharp mass-flux vector"
            )
        residual = _nonnegative_rational("residual", self.residual)
        enclosure = _nonnegative_rational("enclosure", self.enclosure)
        object.__setattr__(self, "residual", residual)
        object.__setattr__(self, "enclosure", enclosure)
        object.__setattr__(self, "mass", _optional_rational("mass", self.mass))
        object.__setattr__(
            self,
            "mass_gradient",
            _optional_vector2("mass_gradient", self.mass_gradient),
        )
        object.__setattr__(
            self,
            "equation_defect",
            _optional_vector2("equation_defect", self.equation_defect),
        )
        expected_valid = residual <= enclosure
        if ledger_valid is not expected_valid:
            raise Def1Stab1Error(
                "ledger_valid must match residual <= enclosure"
            )
        if ledger_valid:
            if (
                self.protocol_token is not None
                or self.veto_class is not None
                or self.typed_reason is not None
            ):
                raise Def1Stab1Error(
                    "a valid mass-flux ledger cannot carry a veto token"
                )
            return
        if self.protocol_token != PROTOCOL_MASS_FLUX_INCONSISTENCY_TOKEN:
            raise Def1Stab1Error(
                "mass-flux veto must map to the existing protocol token "
                "mass_flux_inconsistency"
            )
        if self.veto_class != STOPPED_MASS_FLUX_INCONSISTENCY:
            raise Def1Stab1Error(
                "mass-flux veto class must remain stopped_mass_flux_inconsistency"
            )
        if self.typed_reason not in MASS_FLUX_INCONSISTENCY_REASONS:
            raise Def1Stab1Error("mass-flux veto reason is not in the typed map")
        if self.typed_reason != "conservation_residual_exceeds_enclosure":
            raise Def1Stab1Error(
                "a residual that exceeds its enclosure must use "
                "conservation_residual_exceeds_enclosure"
            )


@dataclass(frozen=True, slots=True)
class Def1BooleanRecord:
    """Nine independently supplied DEF1 booleans; never a collapsed label."""

    resolved_activation: bool
    control_dominance: bool
    resolved_trapped_interval: bool
    resolved_complete_defocusing: bool
    direct_Raychaudhuri_agreement: bool
    all_health_constraints_scales_valid: bool
    mass_flux_ledger_valid: bool
    three_resolution_convergence: bool
    two_method_agreement: bool

    def __post_init__(self) -> None:
        for name in DEF1_BOOLEAN_NAMES:
            _boolean(name, getattr(self, name))

    def as_mapping(self) -> dict[str, bool]:
        return {name: bool(getattr(self, name)) for name in DEF1_BOOLEAN_NAMES}


def _contribution(premise: QComponentPremise) -> Fraction:
    return premise.sensitivity.abs_upper() * premise.input_error.upper


def assemble_q_error_budget(
    premises: Sequence[QComponentPremise],
    **kwargs: object,
) -> AssembledQErrorBudget:
    """Assemble a conservative ``Q``-error upper bound from supplied premises.

    The bound is independent of any measured positive ``Q``.  Missing,
    extra, duplicate-owned, negative, or nonfinite premises refuse.
    Correlated sources across different components remain allowed under the
    conservative sum.  Richardson
    and IMP1 admission debit remain conditional and are never relabeled as
    a global PDE error.
    """

    if kwargs:
        unexpected = ", ".join(sorted(kwargs))
        raise Def1Stab1Error(
            "Q-error assembly refuses unexpected arguments "
            f"({unexpected}); measured Q cannot enter the bound"
        )
    if not isinstance(premises, Sequence) or isinstance(premises, (str, bytes)):
        raise TypeError("premises must be a sequence of QComponentPremise")
    seen: dict[str, QComponentPremise] = {}
    for index, premise in enumerate(premises):
        if not isinstance(premise, QComponentPremise):
            raise TypeError(f"premises[{index}] must be a QComponentPremise")
        if premise.component in seen:
            raise Def1Stab1Error(
                f"duplicate Q-error component ownership for {premise.component}; "
                "correlated sources remain allowed under a conservative sum"
            )
        seen[premise.component] = premise
    missing = [name for name in ERROR_BUDGET_COMPONENTS if name not in seen]
    extra = sorted(name for name in seen if name not in ERROR_BUDGET_COMPONENTS)
    if missing or extra:
        raise Def1Stab1Error(
            "Q-error premises must be exactly the frozen component list; "
            f"missing={missing} extra={extra}"
        )
    exact_components = tuple(
        (name, _contribution(seen[name])) for name in ERROR_BUDGET_COMPONENTS
    )
    exact_total = sum((value for _, value in exact_components), Q(0))
    component_floats = tuple(
        _outward_upper_float(value) for _, value in exact_components
    )
    budget = RaychaudhuriErrorBudget(
        **{
            name: component_floats[index]
            for index, name in enumerate(ERROR_BUDGET_COMPONENTS)
        },
        total_upper_bound=_dominating_total(component_floats, exact_total),
    )
    statuses = tuple((name, seen[name].status) for name in ERROR_BUDGET_COMPONENTS)
    sources = tuple((name, seen[name].source) for name in ERROR_BUDGET_COMPONENTS)
    return AssembledQErrorBudget(
        budget=budget,
        exact_components=exact_components,
        exact_total=exact_total,
        component_statuses=statuses,
        component_sources=sources,
        all_components_declared_proven_enclosures=all(
            status == PREMISE_STATUS_PROVEN for _, status in statuses
        ),
        global_pde_error_certified=False,
        used_measured_q=False,
        richardson_treated_as_global_pde_error=False,
        imp1_admission_debit_treated_as_global_pde_error=False,
    )


def explicit_zero_input_error() -> Interval:
    """Explicit zero enclosure for exact algebraic controls.

    This is not a certified scientific input and must not fill missing evidence.
    """

    return Interval.singleton(0)


def inject_component_input_error(
    premises: Sequence[QComponentPremise],
    component: str,
    extra_error: object,
) -> tuple[QComponentPremise, ...]:
    """Return premises with an injected nonnegative extra input error."""

    extra = _nonnegative_rational("injected extra error", extra_error)
    if component not in ERROR_BUDGET_COMPONENTS:
        raise Def1Stab1Error(f"unknown Q-error component {component!r}")
    updated: list[QComponentPremise] = []
    found = False
    for premise in premises:
        if not isinstance(premise, QComponentPremise):
            raise TypeError("premises must contain QComponentPremise values")
        if premise.component != component:
            updated.append(premise)
            continue
        found = True
        new_error = Interval(
            premise.input_error.lower, premise.input_error.upper + extra
        )
        updated.append(replace(premise, input_error=new_error))
    if not found:
        raise Def1Stab1Error(
            f"cannot inject into missing component {component}; missing "
            "premises refuse rather than default to zero"
        )
    return tuple(updated)


def _lorentzian_symmetric_matrix(name: str, value: object) -> Matrix2:
    matrix = _matrix2(name, value)
    if matrix[0][1] != matrix[1][0]:
        raise Def1Stab1Error(f"{name} must be symmetric")
    determinant = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]
    if determinant >= 0:
        raise Def1Stab1Error(
            f"{name} must be Lorentzian with strictly negative determinant"
        )
    return matrix


def inverse_base_metric(metric: object) -> Matrix2:
    """Return the inverse of a Lorentzian 2D base metric."""

    h = _lorentzian_symmetric_matrix("base metric", metric)
    determinant = h[0][0] * h[1][1] - h[0][1] * h[1][0]
    return (
        (h[1][1] / determinant, -h[0][1] / determinant),
        (-h[1][0] / determinant, h[0][0] / determinant),
    )


def inverse_base_metric_from_adm(
    lapse: object, shift: object, radial_scale: object
) -> Matrix2:
    """Exact inverse 2D metric from positive ADM lapse, shift, and ``lambda``."""

    alpha = _positive_rational("lapse", lapse)
    v = _rational("shift", shift)
    lam = _positive_rational("radial scale", radial_scale)
    lam2 = lam * lam
    metric: Matrix2 = (
        (-(alpha * alpha) + lam2 * v * v, lam2 * v),
        (lam2 * v, lam2),
    )
    return inverse_base_metric(metric)


def mix_base_tensor(covariant: object, inverse_metric: object) -> Matrix2:
    """Raise the second index of a 2D covariant tensor with ``h^{cb}``.

    The inverse metric must be symmetric and Lorentzian.  The covariant
    tensor is not required to be symmetric.
    """

    lowered = _matrix2("covariant base tensor", covariant)
    inverse = _lorentzian_symmetric_matrix("inverse base metric", inverse_metric)
    return (
        (
            lowered[0][0] * inverse[0][0] + lowered[0][1] * inverse[1][0],
            lowered[0][0] * inverse[0][1] + lowered[0][1] * inverse[1][1],
        ),
        (
            lowered[1][0] * inverse[0][0] + lowered[1][1] * inverse[1][0],
            lowered[1][0] * inverse[0][1] + lowered[1][1] * inverse[1][1],
        ),
    )


def misner_sharp_mass(
    radius: object,
    inverse_metric: object,
    radius_derivatives: object,
) -> Fraction:
    """Exact ``m = R (1 - h^{ab} R_a R_b) / 2`` on the 2D base."""

    areal = _positive_rational("areal radius", radius)
    inverse = _lorentzian_symmetric_matrix("inverse base metric", inverse_metric)
    gradient = _vector2("radius derivatives", radius_derivatives)
    contraction = (
        inverse[0][0] * gradient[0] * gradient[0]
        + (inverse[0][1] + inverse[1][0]) * gradient[0] * gradient[1]
        + inverse[1][1] * gradient[1] * gradient[1]
    )
    return areal * (1 - contraction) / 2


def trace_adjusted_radial_projector(
    mixed_tensor: object, radius_derivatives: object
) -> Vector2:
    """Return ``(T_a{}^b - delta_a{}^b T_c{}^c) R_b`` with 2D trace."""

    mixed = _matrix2("mixed base tensor", mixed_tensor)
    gradient = _vector2("radius derivatives", radius_derivatives)
    trace = _trace(mixed)
    projector = []
    for index in range(BASE_DIM):
        projector.append(
            (mixed[index][0] - (trace if index == 0 else Q(0))) * gradient[0]
            + (mixed[index][1] - (trace if index == 1 else Q(0))) * gradient[1]
        )
    return projector[0], projector[1]


def misner_sharp_mass_gradient(
    radius: object,
    mixed_einstein: object,
    radius_derivatives: object,
) -> Vector2:
    """Exact ``d_a m = (R^2 / 2) (G_a{}^b - delta_a{}^b G_c{}^c) R_b``."""

    areal = _positive_rational("areal radius", radius)
    projector = trace_adjusted_radial_projector(mixed_einstein, radius_derivatives)
    factor = (areal * areal) / 2
    return factor * projector[0], factor * projector[1]


def misner_sharp_equation_defect(
    radius: object,
    coupling_F: object,
    mixed_equation_residual: object,
    radius_derivatives: object,
) -> Vector2:
    """Defect of ``E = F G - T_eff``: ``R^2 / (2 F)`` times the projector of ``E``.

    ``F`` must be strictly positive.  Hamiltonian and momentum residuals
    alone are insufficient for the full spacetime mass-flux vector: in an
    orthonormal frame, ``P[E]_n = -E_ss N(R) + E_ns S(R)`` uses ``E_ss``
    for the normal/time component, while
    ``P[E]_s = -E_ns N(R) + E_nn S(R)``.  Mixed tensors are not required
    to be symmetric.
    """

    areal = _positive_rational("areal radius", radius)
    coupling = _rational("coupling F", coupling_F)
    if coupling <= 0:
        raise Def1Stab1Error("coupling F must be strictly positive")
    projector = trace_adjusted_radial_projector(
        mixed_equation_residual, radius_derivatives
    )
    factor = (areal * areal) / (2 * coupling)
    return factor * projector[0], factor * projector[1]


def misner_sharp_identity_holds(
    *,
    radius: object,
    coupling_F: object,
    mixed_einstein: object,
    mixed_effective_source: object,
    mixed_equation_residual: object,
    radius_derivatives: object,
) -> bool:
    """Algebraic ``E = F G - T_eff`` split on supplied mixed 2D tensors.

    This checks that the supplied residual equals ``F G - T_eff`` and that
    the projector identity splits as ``dm = flux + defect``.  It is not an
    independent geometric, trajectory, or conservation-ledger check.
    """

    coupling = _rational("coupling F", coupling_F)
    if coupling <= 0:
        raise Def1Stab1Error("coupling F must be strictly positive")
    einstein = _matrix2("mixed Einstein", mixed_einstein)
    source = _matrix2("mixed effective source", mixed_effective_source)
    residual = _matrix2("mixed equation residual", mixed_equation_residual)
    reconstructed = _subtract_matrix(_scale_matrix(coupling, einstein), source)
    if reconstructed != residual:
        raise Def1Stab1Error(
            "supplied E does not equal F G - T_eff on the 2D base"
        )
    gradient = misner_sharp_mass_gradient(radius, einstein, radius_derivatives)
    flux = misner_sharp_equation_defect(
        radius, coupling, source, radius_derivatives
    )
    defect = misner_sharp_equation_defect(
        radius, coupling, residual, radius_derivatives
    )
    return gradient == (flux[0] + defect[0], flux[1] + defect[1])


def map_mass_flux_inconsistency_token(reason: str) -> str:
    """Map a typed mass-flux reason onto the existing protocol token."""

    if reason not in MASS_FLUX_INCONSISTENCY_REASONS:
        raise Def1Stab1Error("mass-flux reason is outside the typed veto map")
    return PROTOCOL_MASS_FLUX_INCONSISTENCY_TOKEN


def assess_mass_flux_ledger(
    *,
    radius: object,
    inverse_metric: object,
    radius_derivatives: object,
    coupling_F: object,
    mixed_equation_residual: object,
    delta_mass: object,
    integrated_flux: object,
    residual_enclosure: object,
    hamiltonian_momentum_only: bool = False,
    mixed_einstein: object | None = None,
) -> MassFluxAssessment:
    """Assess ``|Delta m - F_integrated| <= E_M`` with the 2D projector.

    Hamiltonian/momentum-only evidence is refused because it does not
    determine the full spacetime mass-flux vector.  Protocol tokens are
    not renamed.  The inverse metric must be symmetric and Lorentzian.
    """

    if hamiltonian_momentum_only:
        raise Def1Stab1Error(
            "Hamiltonian/momentum residuals alone are insufficient for the "
            "full spacetime Misner-Sharp mass-flux vector"
        )
    mixed = _matrix2("mixed equation residual", mixed_equation_residual)
    mass = misner_sharp_mass(radius, inverse_metric, radius_derivatives)
    defect = misner_sharp_equation_defect(
        radius, coupling_F, mixed, radius_derivatives
    )
    gradient = None
    if mixed_einstein is not None:
        gradient = misner_sharp_mass_gradient(
            radius, mixed_einstein, radius_derivatives
        )
    residual = abs(
        _rational("delta mass", delta_mass)
        - _rational("integrated flux", integrated_flux)
    )
    enclosure = _nonnegative_rational("mass-flux enclosure", residual_enclosure)
    if residual <= enclosure:
        return MassFluxAssessment(
            mass=mass,
            mass_gradient=gradient,
            equation_defect=defect,
            residual=residual,
            enclosure=enclosure,
            ledger_valid=True,
            protocol_token=None,
            veto_class=None,
            typed_reason=None,
            used_hamiltonian_momentum_only=False,
        )
    return MassFluxAssessment(
        mass=mass,
        mass_gradient=gradient,
        equation_defect=defect,
        residual=residual,
        enclosure=enclosure,
        ledger_valid=False,
        protocol_token=map_mass_flux_inconsistency_token(
            "conservation_residual_exceeds_enclosure"
        ),
        veto_class=STOPPED_MASS_FLUX_INCONSISTENCY,
        typed_reason="conservation_residual_exceeds_enclosure",
        used_hamiltonian_momentum_only=False,
    )


def assess_activation(
    *,
    max_abs_phi: object,
    activation_error: object,
    s_ref: object,
    gr0_phi_policy: str,
    matched_controls: Sequence[MatchedControlSample],
    historical_planted_seed: object | None = None,
    **kwargs: object,
) -> ActivationAssessment:
    """Assess ``A = max|phi| / S_ref`` with a declared GR-0 phi policy.

    Canonical GR-0 has ``phi = 0``.  The historical planted calibration seed
    is a different declared policy.  ``0/0`` is refused, and the policy is
    never selected silently.
    """

    if kwargs:
        raise Def1Stab1Error(
            "activation assessment refuses unexpected arguments; GR-0 policy "
            "cannot be silently changed or selected"
        )
    if gr0_phi_policy not in GR0_PHI_POLICIES:
        raise Def1Stab1Error("GR-0 phi policy must be declared and recognized")
    if gr0_phi_policy == GR0_PHI_POLICY_CANONICAL_ZERO:
        if historical_planted_seed is not None:
            raise Def1Stab1Error(
                "canonical GR-0 phi=0 cannot carry a historical planted seed"
            )
    else:
        if historical_planted_seed is None:
            raise Def1Stab1Error(
                "historical planted calibration seed is missing; it is not "
                "selected silently and 0/0 is refused"
            )
        _positive_rational(
            "historical planted calibration seed", historical_planted_seed
        )
    reference = _positive_rational("case reference S_ref", s_ref)
    peak = _nonnegative_rational("max_abs_phi", max_abs_phi)
    error = _nonnegative_rational("activation_error", activation_error)
    if not isinstance(matched_controls, Sequence) or isinstance(
        matched_controls, (str, bytes)
    ):
        raise TypeError("matched_controls must be a sequence")
    if len(matched_controls) != len(REQUIRED_MATCHED_CONTROL_NAMES):
        raise Def1Stab1Error("matched GR-0 and SGB-L controls are both required")
    seen: dict[str, MatchedControlSample] = {}
    for control in matched_controls:
        if not isinstance(control, MatchedControlSample):
            raise TypeError("matched controls must be MatchedControlSample values")
        if control.name in seen:
            raise Def1Stab1Error(f"overlapping matched-control sample {control.name}")
        seen[control.name] = control
    if tuple(sorted(seen)) != tuple(sorted(REQUIRED_MATCHED_CONTROL_NAMES)):
        raise Def1Stab1Error("matched controls must be exactly GR-0 and SGB-L")
    gr0 = seen["GR-0"]
    if gr0_phi_policy == GR0_PHI_POLICY_CANONICAL_ZERO and gr0.max_abs_phi != 0:
        raise Def1Stab1Error(
            "canonical GR-0 phi=0 forbids a nonzero GR-0 numerator; the "
            "historical planted seed is a different declared policy"
        )
    raw = peak / reference
    lower = raw - error
    control_uppers = []
    for name in REQUIRED_MATCHED_CONTROL_NAMES:
        control = seen[name]
        control_uppers.append(
            control.max_abs_phi / reference + control.activation_error
        )
    matched_upper = max(control_uppers)
    raw_float = float(raw)
    if not isfinite(raw_float) or raw_float < 0.0:
        raise Def1Stab1Error("raw activation factor left the finite range")
    stored_lower = _outward_lower_float(lower)
    stored_control = _outward_upper_float(matched_upper)
    # Flags are decided from saved outward bounds so the five-field DTO
    # cannot advertise a strict comparison after rounding overlap.
    # Exact inequalities are never weakened: a True flag still requires
    # the exact comparison and the saved-bound comparison.
    threshold = (
        lower >= MIN_ACTIVATION_LOWER_BOUND
        and Q(stored_lower) >= MIN_ACTIVATION_LOWER_BOUND
    )
    dominance = lower > matched_upper and Q(stored_lower) > Q(stored_control)
    return ActivationAssessment(
        raw_factor=raw_float,
        lower_bound=stored_lower,
        matched_control_upper_bound=stored_control,
        threshold_passed=threshold,
        control_dominance_passed=dominance,
    )


def trappedness_margin_passed(
    theta_plus: object,
    theta_minus: object,
    error_plus: object,
    error_minus: object,
) -> bool:
    """Return True iff both expansions satisfy ``theta < -4 E_theta`` strictly."""

    plus = _rational("theta_plus", theta_plus)
    minus = _rational("theta_minus", theta_minus)
    err_plus = _nonnegative_rational("E_theta_plus", error_plus)
    err_minus = _nonnegative_rational("E_theta_minus", error_minus)
    return plus < -TRAPPEDNESS_ERROR_FACTOR * err_plus and minus < (
        -TRAPPEDNESS_ERROR_FACTOR * err_minus
    )


def complete_q_margin_passed(min_q: object, q_error: object) -> bool:
    """Return True iff ``min Q > 4 E_Q > 0``."""

    minimum = _rational("min Q", min_q)
    error = _nonnegative_rational("E_Q", q_error)
    return error > 0 and minimum > COMPLETE_Q_ERROR_FACTOR * error


def affine_interval_and_sample_gate(
    affine_parameters: Sequence[object],
    length_scale_L0: object,
) -> AffineIntervalAssessment:
    """Require eight strictly increasing samples and ``Delta lambda / L0 >= 1/64``."""

    if not isinstance(affine_parameters, Sequence) or isinstance(
        affine_parameters, (str, bytes)
    ):
        raise TypeError("affine_parameters must be a sequence")
    if len(affine_parameters) == 0:
        raise Def1Stab1Error("affine sample list is missing")
    samples = tuple(
        _rational(f"affine_parameters[{index}]", value)
        for index, value in enumerate(affine_parameters)
    )
    for index in range(1, len(samples)):
        if samples[index] <= samples[index - 1]:
            raise Def1Stab1Error(
                "affine samples must be strictly increasing; overlapping or "
                "unordered samples refuse"
            )
    scale = _positive_rational("L0", length_scale_L0)
    interval = (samples[-1] - samples[0]) / scale
    samples_passed = len(samples) >= MIN_AFFINE_SAMPLES
    interval_passed = interval >= MIN_AFFINE_INTERVAL_OVER_L0
    return AffineIntervalAssessment(
        sample_count=len(samples),
        interval_over_L0=interval,
        samples_passed=samples_passed,
        interval_passed=interval_passed,
        passed=samples_passed and interval_passed,
    )


def direct_raychaudhuri_routes_agree(direct: object, assembled: object) -> bool:
    """Return True iff ``|direct - assembled| <= 1/100000000``."""

    residual = abs(_rational("direct Q", direct) - _rational("assembled Q", assembled))
    return residual <= DIRECT_RAYCHAUDHURI_AGREEMENT_MAX


def affine_normalization_residual_passed(residual: object) -> bool:
    """Return True iff the affine residual is at most ``1/10000000000``."""

    return abs(_rational("affine residual", residual)) <= AFFINE_NORMALIZATION_RESIDUAL_MAX


def declared_resolution_count_meets_minimum(count: object) -> bool:
    """Count helper only; it does not certify three-resolution convergence."""

    if isinstance(count, bool) or not isinstance(count, Integral):
        raise TypeError("resolution count must be an integer")
    if int(count) < 0:
        raise Def1Stab1Error("resolution count must be nonnegative")
    return int(count) >= MIN_NESTED_RESOLUTIONS


def declared_method_count_meets_minimum(count: object) -> bool:
    """Count helper only; it does not certify two-method agreement."""

    if isinstance(count, bool) or not isinstance(count, Integral):
        raise TypeError("method count must be an integer")
    if int(count) < 0:
        raise Def1Stab1Error("method count must be nonnegative")
    return int(count) >= MIN_INDEPENDENT_METHODS


def observed_order_meets_minimum(order: object) -> bool:
    """Order helper only; it does not certify three-resolution convergence."""

    return _rational("observed order", order) >= MIN_OBSERVED_ORDER


def activation_threshold_passed(lower_bound: object) -> bool:
    """Return True iff ``A - E_A >= 2``."""

    return _rational("activation lower bound", lower_bound) >= MIN_ACTIVATION_LOWER_BOUND


def control_dominance_passed(
    lower_bound: object, matched_control_upper_bound: object
) -> bool:
    """Return True iff the candidate lower bound strictly exceeds the control upper."""

    return _rational("activation lower bound", lower_bound) > _nonnegative_rational(
        "matched control upper bound", matched_control_upper_bound
    )


def record_def1_booleans(
    *,
    resolved_activation: bool,
    control_dominance: bool,
    resolved_trapped_interval: bool,
    resolved_complete_defocusing: bool,
    direct_Raychaudhuri_agreement: bool,
    all_health_constraints_scales_valid: bool,
    mass_flux_ledger_valid: bool,
    three_resolution_convergence: bool,
    two_method_agreement: bool,
) -> Def1BooleanRecord:
    """Record all nine DEF1 booleans; missing entries cannot be inferred."""

    return Def1BooleanRecord(
        resolved_activation=resolved_activation,
        control_dominance=control_dominance,
        resolved_trapped_interval=resolved_trapped_interval,
        resolved_complete_defocusing=resolved_complete_defocusing,
        direct_Raychaudhuri_agreement=direct_Raychaudhuri_agreement,
        all_health_constraints_scales_valid=all_health_constraints_scales_valid,
        mass_flux_ledger_valid=mass_flux_ledger_valid,
        three_resolution_convergence=three_resolution_convergence,
        two_method_agreement=two_method_agreement,
    )


def record_def1_booleans_from_mapping(
    values: Mapping[str, object],
) -> Def1BooleanRecord:
    """Build the nine-boolean record from an explicit complete mapping."""

    if not isinstance(values, Mapping):
        raise TypeError("DEF1 booleans must be supplied as a mapping")
    names = set(values)
    required = set(DEF1_BOOLEAN_NAMES)
    if names != required:
        raise Def1Stab1Error(
            "DEF1 booleans must be supplied independently and completely; "
            f"missing={sorted(required - names)} extra={sorted(names - required)}. "
            "Health and convergence cannot be inferred from margins."
        )
    return record_def1_booleans(
        **{name: _boolean(name, values[name]) for name in DEF1_BOOLEAN_NAMES}
    )


def zero_error_premises(
    *,
    status: str = PREMISE_STATUS_PROVEN,
    source: str = SOURCE_EXACT_ALGEBRA,
) -> tuple[QComponentPremise, ...]:
    """Explicit zero-error premises for exact algebraic controls.

    This helper is not a certified scientific enclosure and must not fill
    missing physical evidence.
    """

    if status == PREMISE_STATUS_CONDITIONAL and source in PROVEN_SOURCES:
        raise Def1Stab1Error("proven sources cannot be labeled conditional")
    if status == PREMISE_STATUS_PROVEN and source in CONDITIONAL_SOURCES:
        raise Def1Stab1Error("conditional sources cannot be labeled proven")
    return tuple(
        QComponentPremise(
            component=name,
            sensitivity=Interval.singleton(1),
            input_error=explicit_zero_input_error(),
            status=status,
            source=source,
        )
        for name in ERROR_BUDGET_COMPONENTS
    )


__all__ = [
    "AFFINE_NORMALIZATION_RESIDUAL_MAX",
    "AssembledQErrorBudget",
    "ActivationAssessment",
    "AffineIntervalAssessment",
    "COMPLETE_Q_ERROR_FACTOR",
    "CONDITIONAL_SOURCES",
    "DEF1_BOOLEAN_NAMES",
    "DIRECT_RAYCHAUDHURI_AGREEMENT_MAX",
    "Def1BooleanRecord",
    "Def1Stab1Error",
    "ERROR_BUDGET_COMPONENTS",
    "GR0_PHI_POLICIES",
    "GR0_PHI_POLICY_CANONICAL_ZERO",
    "GR0_PHI_POLICY_HISTORICAL_PLANTED_SEED",
    "MIN_ACTIVATION_LOWER_BOUND",
    "MIN_AFFINE_INTERVAL_OVER_L0",
    "MIN_AFFINE_SAMPLES",
    "MIN_INDEPENDENT_METHODS",
    "MIN_NESTED_RESOLUTIONS",
    "MIN_OBSERVED_ORDER",
    "MassFluxAssessment",
    "MatchedControlSample",
    "PREMISE_STATUS_CONDITIONAL",
    "PREMISE_STATUS_PROVEN",
    "PROTOCOL_MASS_FLUX_INCONSISTENCY_TOKEN",
    "PROVEN_SOURCES",
    "QComponentPremise",
    "RaychaudhuriErrorBudget",
    "SOURCE_DECLARED_CONDITIONAL",
    "SOURCE_EXACT_ALGEBRA",
    "SOURCE_IMP1_ADMISSION_DEBIT",
    "SOURCE_RICHARDSON",
    "SOURCE_SUPPLIED_CERTIFIED",
    "STOPPED_MASS_FLUX_INCONSISTENCY",
    "TRAPPEDNESS_ERROR_FACTOR",
    "activation_threshold_passed",
    "affine_interval_and_sample_gate",
    "affine_normalization_residual_passed",
    "assemble_q_error_budget",
    "assess_activation",
    "assess_mass_flux_ledger",
    "explicit_zero_input_error",
    "complete_q_margin_passed",
    "control_dominance_passed",
    "declared_method_count_meets_minimum",
    "declared_resolution_count_meets_minimum",
    "direct_raychaudhuri_routes_agree",
    "inject_component_input_error",
    "inverse_base_metric",
    "inverse_base_metric_from_adm",
    "map_mass_flux_inconsistency_token",
    "misner_sharp_equation_defect",
    "misner_sharp_identity_holds",
    "misner_sharp_mass",
    "misner_sharp_mass_gradient",
    "mix_base_tensor",
    "observed_order_meets_minimum",
    "record_def1_booleans",
    "record_def1_booleans_from_mapping",
    "trace_adjusted_radial_projector",
    "trappedness_margin_passed",
    "zero_error_premises",
]
