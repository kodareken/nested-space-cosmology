"""Prospective Taylor-model enclosure of the nominal SGB-L constraint ODE.

This owner does not edit the Picard or chart instruments.  It freezes a
finite Taylor order and a subdivision/resource/bit policy, then evaluates
that whole policy on the unchanged ``A_chi=3`` family, support ``[10,14]``,
16 base cells, ``(lambda,k,C)`` graph, physical domains, affine H/M
equations, and strict ``C<1`` gate.

Every required ``r``/state derivative is produced by exact interval
automatic differentiation (univariate jets).  Compact-bump derivatives
use the same jets, with a Lagrange remainder; sampled finite differences
are not a proof.  Predecessor depth-ladder hashes are regression only.

Aggregate health, ``FRZ1``, ``PREF1``, execution, and holdout remain false.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from hashlib import sha256
from json import dumps
from numbers import Integral
from time import perf_counter
from types import MappingProxyType
from typing import Any, Mapping, Sequence

from .exact_interval import Interval, interval
from .sgb1_ctl1_initial_health import (
    BUMP_X_DERIVATIVE_BOUND,
    BUMP_XX_DERIVATIVE_BOUND,
    CERTIFICATE_CHART_KIND,
    DECLARED_BASE_CELLS,
    DECLARED_C_DOMAIN,
    DECLARED_J_DOMAIN,
    DECLARED_K_DOMAIN,
    DECLARED_LAMBDA_DOMAIN,
    SGBLExactInitialSlice,
)
from .sgb1_ctl1_trap_refinement2 import (
    PREDECESSOR_NOMINAL_CONTRACT_SHA256 as PREDECESSOR_TRAP_REFINEMENT_CONTRACT_SHA256,
)


Q = Fraction
ZERO = Interval.singleton(0)
ONE = Interval.singleton(1)
RationalLike = Fraction | int
INSTRUMENT_ID = "FGC-1-SGB1-CTL1-TRAP-TAYLOR"
PRODUCT_BOX_CHART_KIND = "lambda_k"
DECLARED_TAYLOR_ORDER = 4
DECLARED_TAYLOR_LEVEL_DEPTHS = (2, 4)
DECLARED_MAX_TAYLOR_PICARD_ITERATIONS = 8
DECLARED_EULER_INFLATION = 2
DECLARED_DYADIC_DENOMINATOR = 2 ** 64
DECLARED_MAX_RATIONAL_BITS = 16384
DECLARED_COMPACTNESS_STRICT_UPPER = 1
DECLARED_JET_CALLS_PER_PARTITION = 2 + DECLARED_MAX_TAYLOR_PICARD_ITERATIONS
NOMINAL_CHI_AMPLITUDE = Q(3)
SMALL_AMPLITUDE_CONTROL = Q(1, 8)
FACTORIAL_EIGHT = 40320
PREDECESSOR_TRAP_REFINEMENT2_CONTRACT_SHA256 = (
    "2cc8ba39a756afc433f1b04ce6fa3e26f45d795f2ba4c617ec080e7342cc396d"
)
LEVEL_CLASSIFICATIONS = frozenset(
    {
        "continuous_no_initial_trapped_sphere",
        "interval_inconclusive",
        "resource_limit",
        "domain_error",
    }
)
RESOURCE_REASONS = frozenset(
    {
        "max_cells",
        "max_jet_evaluations",
        "max_rational_bit_length",
    }
)
CONTRACT_STOP_REASONS = frozenset(
    {
        "changed_policy",
        "reordered_levels",
        "skipped_level",
        "tampered_contract",
        "gapped_inventory",
        "wrong_predecessor",
        "wrong_remainder",
        "omitted_derivative",
    }
)
EVALUATION_STOP_REASONS = frozenset(
    {
        "domain_error",
        "resource_limit",
        "interval_inconclusive",
    }
)
ODE_INCONCLUSIVE_REASONS = frozenset(
    {
        "taylor_self_map_failed",
        "jacobian_diagonal_contains_zero",
        "physical_domain_escape",
        "compactness_enclosure_not_below_one",
        "residual_enclosure_misses_origin",
        "exterior_compactness_not_below_one",
        "zero_in_interval_reciprocal",
        "taylor_remainder_inconsistent",
        "compactness_invariant_misses_origin",
        "bump_remainder_not_enclosed",
        "incomplete_support_coverage",
    }
)
MISSING_THEOREMS = MappingProxyType(
    {
        "taylor_self_map_failed": (
            "strict_Taylor_model_self_map_of_the_affine_constraint_ODE_on_the_declared_jet_tree"
        ),
        "jacobian_diagonal_contains_zero": (
            "affine_H_L_and_M_k_diagonals_excluding_zero_on_every_declared_cell_box"
        ),
        "physical_domain_escape": (
            "Taylor_image_remaining_inside_the_declared_physical_state_domain"
        ),
        "compactness_enclosure_not_below_one": (
            "correlated_Taylor_model_compactness_C_strictly_below_one_on_the_validated_graph"
        ),
        "residual_enclosure_misses_origin": (
            "affine_H_and_M_residual_enclosures_containing_the_origin_on_every_cell"
        ),
        "exterior_compactness_not_below_one": (
            "analytic_exterior_C_equals_two_M_over_r_with_maximum_strictly_below_one"
        ),
        "zero_in_interval_reciprocal": (
            "interval_reciprocals_of_r_and_lambda_remaining_defined_on_every_cell_jet"
        ),
        "taylor_remainder_inconsistent": (
            "Lagrange_remainder_Taylor_model_intersecting_the_first_order_existence_box"
        ),
        "compactness_invariant_misses_origin": (
            "algebraic_and_propagated_compactness_enclosures_overlapping_with_invariant_residual_zero"
        ),
        "bump_remainder_not_enclosed": (
            "compact_bump_Lagrange_remainder_from_exact_interval_jets_or_proved_global_derivative_bounds"
        ),
        "incomplete_support_coverage": (
            "complete_compact_support_covering_by_Taylor_model_cells"
        ),
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


def _as_interval(value: object) -> Interval:
    if type(value) is Interval:
        return value
    if type(value) in (int, Fraction):
        return Interval.singleton(Fraction(value))
    raise TypeError("value must be an Interval or exact rational")


def _fraction_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def _interval_text(value: Interval) -> tuple[str, str]:
    return (_fraction_text(value.lower), _fraction_text(value.upper))


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
    """Valid outward rounding onto the declared dyadic grid.  This is a resource."""

    return Interval(_floor_dyadic(value.lower, denominator), _ceil_dyadic(value.upper, denominator))


def _factorials(limit: int) -> tuple[int, ...]:
    values = [1]
    acc = 1
    for index in range(1, limit + 1):
        acc *= index
        values.append(acc)
    return tuple(values)


FACTORIAL = _factorials(24)


def _binom(n: int, k: int) -> int:
    if k < 0 or k > n:
        return 0
    return FACTORIAL[n] // (FACTORIAL[k] * FACTORIAL[n - k])


def _eval_poly(coefficients: Sequence[Interval], argument: Interval) -> Interval:
    result = coefficients[-1]
    for coefficient in reversed(coefficients[:-1]):
        result = result * argument + coefficient
    return result


class SGBLTaylorTrapStop(ValueError):
    """Typed contract stop of the prospective Taylor-model instrument."""

    def __init__(
        self,
        reason: str,
        message: str,
        payload: Mapping[str, Any] | None = None,
    ) -> None:
        if reason not in CONTRACT_STOP_REASONS:
            raise ValueError("unknown Taylor-trap contract stop reason")
        self.reason = reason
        self.payload = dict(payload or {})
        super().__init__(message)


class SGBLTaylorEvaluationStop(ArithmeticError):
    """Typed fail-closed stop of a Taylor-model evaluation."""

    def __init__(
        self,
        reason: str,
        message: str,
        payload: Mapping[str, Any] | None = None,
    ) -> None:
        if reason not in EVALUATION_STOP_REASONS:
            raise ValueError("unknown Taylor-trap evaluation stop reason")
        self.reason = reason
        self.payload = dict(payload or {})
        super().__init__(message)


@dataclass
class _TaylorBudget:
    jet_evaluations: int = 0
    compactness_jet_evaluations: int = 0
    cells: int = 0
    max_bits: int = 0
    correlated_strictly_tighter: bool = False


def _observe(budget: _TaylorBudget, limit: int, *values: Interval) -> None:
    for value in values:
        budget.max_bits = max(budget.max_bits, _interval_bits(value))
        if budget.max_bits > limit:
            raise SGBLTaylorEvaluationStop(
                "resource_limit",
                "Taylor-model endpoints exceeded the declared rational bit cap",
                {
                    "observed": budget.max_bits,
                    "limit": limit,
                    "resource_reason": "max_rational_bit_length",
                },
            )


def _outward_all(budget: _TaylorBudget, limit: int, *values: Interval) -> tuple[Interval, ...]:
    rounded = tuple(_outward(value) for value in values)
    _observe(budget, limit, *rounded)
    return rounded


@dataclass(frozen=True, slots=True, eq=False)
class IntervalJet:
    """Univariate exact-interval derivative jet.

    ``derivatives[k]`` encloses ``f^{(k)}`` at the expansion point, or on a
    domain when the primal is a non-singleton.  Higher derivatives are never
    invented as silent zeros: that is ``omitted_derivative``.
    """

    derivatives: tuple[Interval, ...]

    def __post_init__(self) -> None:
        if not self.derivatives:
            raise ValueError("interval jet must contain the primal")
        object.__setattr__(
            self,
            "derivatives",
            tuple(_as_interval(value) for value in self.derivatives),
        )

    @property
    def order(self) -> int:
        return len(self.derivatives) - 1

    @property
    def primal(self) -> Interval:
        return self.derivatives[0]

    @classmethod
    def constant(cls, value: Interval | RationalLike, order: int) -> "IntervalJet":
        degree = _positive_int("order", order, allow_zero=True)
        return cls((_as_interval(value),) + (ZERO,) * degree)

    @classmethod
    def identity(cls, value: Interval | RationalLike, order: int) -> "IntervalJet":
        degree = _positive_int("order", order, allow_zero=True)
        if degree == 0:
            return cls((_as_interval(value),))
        return cls((_as_interval(value), ONE) + (ZERO,) * (degree - 1))

    def truncate(self, order: int) -> "IntervalJet":
        degree = _positive_int("order", order, allow_zero=True)
        if degree > self.order:
            raise SGBLTaylorTrapStop(
                "omitted_derivative",
                "jet truncation cannot invent missing derivatives",
                {"have": self.order, "need": degree},
            )
        return IntervalJet(self.derivatives[: degree + 1])

    def __bool__(self) -> bool:
        raise TypeError("interval jet truth value is ambiguous; inspect its primal")

    def _coerce(self, other: object) -> "IntervalJet":
        if isinstance(other, IntervalJet):
            if other.order != self.order:
                raise SGBLTaylorTrapStop(
                    "omitted_derivative",
                    "jet arithmetic requires matching derivative orders",
                    {"left": self.order, "right": other.order},
                )
            return other
        return IntervalJet.constant(other, self.order)  # type: ignore[arg-type]

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, IntervalJet):
            return False
        return self.derivatives == other.derivatives

    def __neg__(self) -> "IntervalJet":
        return IntervalJet(tuple(-value for value in self.derivatives))

    def __add__(self, other: object) -> "IntervalJet":
        value = self._coerce(other)
        return IntervalJet(
            tuple(left + right for left, right in zip(self.derivatives, value.derivatives))
        )

    __radd__ = __add__

    def __sub__(self, other: object) -> "IntervalJet":
        return self + (-self._coerce(other))

    def __rsub__(self, other: object) -> "IntervalJet":
        return self._coerce(other) - self

    def __mul__(self, other: object) -> "IntervalJet":
        value = self._coerce(other)
        derivatives = []
        for degree in range(self.order + 1):
            acc = ZERO
            for index in range(degree + 1):
                acc = acc + (
                    _binom(degree, index)
                    * self.derivatives[index]
                    * value.derivatives[degree - index]
                )
            derivatives.append(acc)
        return IntervalJet(tuple(derivatives))

    __rmul__ = __mul__

    def reciprocal(self) -> "IntervalJet":
        if self.primal.contains_zero():
            raise ZeroDivisionError("interval jet reciprocal denominator contains zero")
        inverse = [self.primal.reciprocal()]
        for degree in range(1, self.order + 1):
            acc = ZERO
            for index in range(1, degree + 1):
                acc = acc + (
                    _binom(degree, index)
                    * self.derivatives[index]
                    * inverse[degree - index]
                )
            inverse.append(-acc * inverse[0])
        return IntervalJet(tuple(inverse))

    def __truediv__(self, other: object) -> "IntervalJet":
        return self * self._coerce(other).reciprocal()

    def __rtruediv__(self, other: object) -> "IntervalJet":
        return self._coerce(other) * self.reciprocal()

    def __pow__(self, exponent: object) -> "IntervalJet":
        if isinstance(exponent, bool) or not isinstance(exponent, Integral):
            return NotImplemented
        power = int(exponent)
        if power == 0:
            return IntervalJet.constant(1, self.order)
        if power < 0:
            return self.reciprocal() ** (-power)
        result = self
        for _ in range(power - 1):
            result = result * self
        return result


def sgbl_interval_jet_shift(jet: IntervalJet, count: int, order: int) -> IntervalJet:
    """Return the jet of ``f^{(count)}`` as a function, of the requested order."""

    if type(jet) is not IntervalJet:
        raise TypeError("jet must be IntervalJet")
    shift = _positive_int("count", count, allow_zero=True)
    degree = _positive_int("order", order, allow_zero=True)
    need = degree + shift
    if jet.order < need:
        raise SGBLTaylorTrapStop(
            "omitted_derivative",
            "missing higher derivative for a jet shift",
            {"have": jet.order, "need": need, "shift": shift},
        )
    return IntervalJet(jet.derivatives[shift : shift + degree + 1])


def sgbl_lagrange_remainder(
    remainder_derivative: Interval,
    displacement: Interval,
    order: int,
) -> Interval:
    """Taylor theorem: ``f^{(n+1)}(ξ)/(n+1)! * (x-c)^{n+1}``.

    ``remainder_derivative`` must enclose ``f^{(n+1)}`` on a convex set
    containing ``c`` and ``x``.  This is the Lagrange remainder, not a
    sampled stencil.
    """

    degree = _positive_int("order", order, allow_zero=True)
    if degree + 1 >= len(FACTORIAL):
        raise SGBLTaylorEvaluationStop(
            "resource_limit",
            "Taylor remainder order exceeded the frozen factorial table",
            {"order": degree, "resource_reason": "max_jet_evaluations"},
        )
    derivative = _as_interval(remainder_derivative)
    step = _as_interval(displacement)
    return derivative / FACTORIAL[degree + 1] * (step ** (degree + 1))


@dataclass(frozen=True, slots=True)
class TaylorModel:
    """Polynomial of degree ``n`` plus a Lagrange remainder interval.

    Coefficients are Taylor coefficients ``f^{(k)}(c)/k!``.  The stored
    ``remainder_derivative`` is an enclosure of ``f^{(n+1)}`` on ``domain``.
    """

    coefficients: tuple[Interval, ...]
    remainder_derivative: Interval
    center: Fraction
    domain: Interval

    def __post_init__(self) -> None:
        if not self.coefficients:
            raise ValueError("Taylor model must contain the degree-0 coefficient")
        object.__setattr__(
            self,
            "coefficients",
            tuple(_as_interval(value) for value in self.coefficients),
        )
        object.__setattr__(self, "remainder_derivative", _as_interval(self.remainder_derivative))
        object.__setattr__(self, "center", _fraction("center", self.center))
        object.__setattr__(self, "domain", _as_interval(self.domain))
        if not (
            self.domain.lower <= self.center <= self.domain.upper
        ):
            raise SGBLTaylorEvaluationStop(
                "domain_error",
                "Taylor expansion point must lie in the remainder domain",
            )

    @property
    def order(self) -> int:
        return len(self.coefficients) - 1

    def displacement(self) -> Interval:
        return self.domain - Interval.singleton(self.center)

    def remainder_on_domain(self) -> Interval:
        return _outward(
            sgbl_lagrange_remainder(
                self.remainder_derivative, self.displacement(), self.order
            )
        )

    def polynomial_range(self) -> Interval:
        return _outward(_eval_poly(self.coefficients, self.displacement()))

    def range(self) -> Interval:
        return _outward(self.polynomial_range() + self.remainder_on_domain())

    def at(self, point: Fraction) -> Interval:
        value = _fraction("point", point)
        if value < self.domain.lower or value > self.domain.upper:
            raise SGBLTaylorEvaluationStop(
                "domain_error",
                "Taylor evaluation point escaped the remainder domain",
            )
        step = Interval.singleton(value - self.center)
        polynomial = _eval_poly(self.coefficients, step)
        remainder = sgbl_lagrange_remainder(
            self.remainder_derivative, step, self.order
        )
        return _outward(polynomial + remainder)


def sgbl_taylor_model_from_jet(
    jet: IntervalJet,
    remainder_derivative: Interval,
    center: Fraction,
    domain: Interval,
) -> TaylorModel:
    """Convert a derivative jet into a Taylor model by Taylor's theorem."""

    if type(jet) is not IntervalJet:
        raise TypeError("jet must be IntervalJet")
    coefficients = tuple(
        _outward(jet.derivatives[index] / FACTORIAL[index]) for index in range(jet.order + 1)
    )
    return TaylorModel(
        coefficients=coefficients,
        remainder_derivative=_as_interval(remainder_derivative),
        center=center,
        domain=domain,
    )


def sgbl_require_honest_remainder(
    claimed: Interval,
    remainder_derivative: Interval,
    displacement: Interval,
    order: int,
) -> Interval:
    """Refuse a remainder that does not contain the Lagrange enclosure."""

    independent = _outward(
        sgbl_lagrange_remainder(remainder_derivative, displacement, order)
    )
    claimed_interval = _as_interval(claimed)
    if not independent.subset_of(claimed_interval):
        raise SGBLTaylorTrapStop(
            "wrong_remainder",
            "claimed remainder does not contain the Lagrange remainder enclosure",
            {
                "claimed": _interval_text(claimed_interval),
                "independent": _interval_text(independent),
            },
        )
    return independent


def sgbl_require_complete_jet(jet: IntervalJet, order: int) -> IntervalJet:
    """Refuse a jet that silently dropped a required derivative."""

    if type(jet) is not IntervalJet:
        raise TypeError("jet must be IntervalJet")
    degree = _positive_int("order", order, allow_zero=True)
    if jet.order != degree:
        raise SGBLTaylorTrapStop(
            "omitted_derivative",
            "jet does not carry every derivative through the frozen order",
            {"have": jet.order, "need": degree},
        )
    return jet


@dataclass(frozen=True)
class _BumpTerm:
    coeff: Fraction
    x_power: int
    t_power: int


def _collect_terms(terms: Sequence[_BumpTerm]) -> tuple[_BumpTerm, ...]:
    buckets: dict[tuple[int, int], Fraction] = {}
    for term in terms:
        key = (term.x_power, term.t_power)
        buckets[key] = buckets.get(key, Q(0)) + term.coeff
    return tuple(
        _BumpTerm(coeff, x_power, t_power)
        for (x_power, t_power), coeff in sorted(buckets.items())
        if coeff != 0
    )


def _differentiate_terms(terms: Sequence[_BumpTerm]) -> tuple[_BumpTerm, ...]:
    out: list[_BumpTerm] = []
    for term in terms:
        if term.x_power:
            out.append(_BumpTerm(term.coeff * term.x_power, term.x_power - 1, term.t_power))
        if term.t_power:
            out.append(
                _BumpTerm(
                    term.coeff * 2 * term.t_power,
                    term.x_power + 1,
                    term.t_power + 1,
                )
            )
    return _collect_terms(out)


def _multiply_terms(
    left: Sequence[_BumpTerm], right: Sequence[_BumpTerm]
) -> tuple[_BumpTerm, ...]:
    out = [
        _BumpTerm(
            first.coeff * second.coeff,
            first.x_power + second.x_power,
            first.t_power + second.t_power,
        )
        for first in left
        for second in right
    ]
    return _collect_terms(out)


def _e_neg_u_times_u_plus_one_bound(power: int) -> Fraction:
    """``e^{-u}(u+1)^m <= sum_k binom(m,k) k!`` for ``u>=0``.

    Uses ``e^{-u} u^k <= k!``, which follows from ``e^u >= u^k/k!``.
    """

    degree = _positive_int("power", power, allow_zero=True)
    total = Q(0)
    for index in range(degree + 1):
        total += _binom(degree, index) * FACTORIAL[index]
    return total


def _bump_scaled_derivative_bounds(max_order: int) -> tuple[Fraction, ...]:
    """Global bounds for ``|d^n B/dx^n|`` of ``B=exp(1-1/(1-x^2))`` on ``|x|<=1``.

    Logarithmic derivatives satisfy ``B^{(n)}=Q_n B`` with the recurrence
    ``Q_n = Q_{n-1}' + Q_{n-1} Q_1``, ``Q_1=-2x(1-x^2)^{-2}``.  Each monomial
    ``B x^p (1-x^2)^{-m}`` is bounded by ``e^{-u}(u+1)^m`` with ``|x|<=1``.
    """

    degree = _positive_int("max_order", max_order, allow_zero=True)
    q0 = (_BumpTerm(Q(1), 0, 0),)
    q1 = (_BumpTerm(Q(-2), 1, 2),)
    qs: list[tuple[_BumpTerm, ...]] = [q0, q1]
    while len(qs) <= degree:
        previous = qs[-1]
        qs.append(
            _collect_terms(tuple(_differentiate_terms(previous) + _multiply_terms(previous, q1)))
        )
    bounds = []
    for polynomials in qs[: degree + 1]:
        total = Q(0)
        for term in polynomials:
            total += abs(term.coeff) * _e_neg_u_times_u_plus_one_bound(term.t_power)
        bounds.append(total)
    return tuple(bounds)


DECLARED_BUMP_JET_ORDER = DECLARED_TAYLOR_ORDER + 3
_RAW_BUMP_SCALED_DERIVATIVE_BOUNDS = _bump_scaled_derivative_bounds(DECLARED_BUMP_JET_ORDER)
_ESTABLISHED_BUMP_SCALED_BOUNDS = (
    Q(1),
    Q(BUMP_X_DERIVATIVE_BOUND),
    Q(BUMP_XX_DERIVATIVE_BOUND),
)
BUMP_SCALED_DERIVATIVE_BOUNDS = tuple(
    min(raw, _ESTABLISHED_BUMP_SCALED_BOUNDS[index])
    if index < len(_ESTABLISHED_BUMP_SCALED_BOUNDS)
    else raw
    for index, raw in enumerate(_RAW_BUMP_SCALED_DERIVATIVE_BOUNDS)
)


def _exp_neg_upper(value: Fraction) -> Fraction:
    if value < 0:
        raise SGBLTaylorEvaluationStop("domain_error", "exp(-t) upper bound requires t>=0")
    if value == 0:
        return Q(1)
    if value >= 8:
        return Fraction(FACTORIAL_EIGHT) / value ** 8
    total = Q(1)
    term = Q(1)
    for index in range(1, 17):
        term *= value / index
        total += term
    return 1 / total


def _exp_neg_lower(value: Fraction) -> Fraction:
    if value < 0:
        raise SGBLTaylorEvaluationStop("domain_error", "exp(-t) lower bound requires t>=0")
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
        raise SGBLTaylorEvaluationStop(
            "domain_error",
            "compact-bump exponent enclosure requires a nonnegative argument",
        )
    return Interval(_exp_neg_lower(value.upper), _exp_neg_upper(value.lower))


def _jet_exp_neg(exponent: IntervalJet) -> IntervalJet:
    primal = exponent.primal
    if primal.lower < 0:
        if primal.upper < 0:
            raise SGBLTaylorEvaluationStop(
                "domain_error",
                "compact-bump exponent jet is strictly negative",
            )
        primal = Interval(Q(0), primal.upper)
        exponent = IntervalJet((primal,) + exponent.derivatives[1:])
    values = [_exp_neg_enclosure(exponent.primal)]
    for degree in range(1, exponent.order + 1):
        acc = ZERO
        for index in range(degree):
            acc = acc + (
                _binom(degree - 1, index)
                * exponent.derivatives[index + 1]
                * values[degree - 1 - index]
            )
        values.append(-acc)
    return IntervalJet(tuple(values))


def _global_bump_jet(half_width: Fraction, order: int) -> IntervalJet:
    degree = _positive_int("order", order, allow_zero=True)
    if degree > DECLARED_BUMP_JET_ORDER:
        raise SGBLTaylorTrapStop(
            "omitted_derivative",
            "bump jet order exceeds the frozen global-bound table",
            {"have": DECLARED_BUMP_JET_ORDER, "need": degree},
        )
    width = Interval.singleton(half_width)
    derivatives = [interval(0, 1)]
    for index in range(1, degree + 1):
        bound = BUMP_SCALED_DERIVATIVE_BOUNDS[index]
        derivatives.append(interval(-bound, bound) / (width ** index))
    return IntervalJet(tuple(derivatives))


def _intersect_jets(left: IntervalJet, right: IntervalJet) -> IntervalJet:
    if left.order != right.order:
        raise SGBLTaylorTrapStop(
            "omitted_derivative",
            "jet intersection requires matching orders",
        )
    derivatives = []
    for first, second in zip(left.derivatives, right.derivatives):
        overlap = _intersect(first, second)
        if overlap is None:
            raise SGBLTaylorEvaluationStop(
                "interval_inconclusive",
                "formula jet and global bump bound are disjoint",
                {"obstruction": "bump_remainder_not_enclosed"},
            )
        derivatives.append(overlap)
    return IntervalJet(tuple(derivatives))


def sgbl_compact_bump_interval_jet(
    radius: Interval | IntervalJet,
    *,
    center: Fraction,
    half_width: Fraction,
    order: int,
) -> IntervalJet:
    """Exact interval jet of the compact bump ``B=exp(1-1/(1-x^2))``.

    Interior cells use automatic differentiation of the closed form.
    Cells that touch the support edge use the proved global derivative
    bounds.  At a singleton on or outside the support the jet is zero.
    """

    degree = _positive_int("order", order, allow_zero=True)
    if half_width <= 0:
        raise SGBLTaylorEvaluationStop("domain_error", "bump half-width must be positive")
    if degree > DECLARED_BUMP_JET_ORDER:
        raise SGBLTaylorTrapStop(
            "omitted_derivative",
            "requested bump jet exceeds the frozen derivative table",
            {"have": DECLARED_BUMP_JET_ORDER, "need": degree},
        )
    primal = radius.primal if isinstance(radius, IntervalJet) else _as_interval(radius)
    r_jet = IntervalJet.identity(primal, degree)
    scaled = (primal - center) / half_width
    global_jet = _global_bump_jet(half_width, degree)
    if scaled.upper <= -1 or scaled.lower >= 1:
        return IntervalJet.constant(0, degree)
    if scaled.is_singleton() and (scaled.lower <= -1 or scaled.lower >= 1):
        return IntervalJet.constant(0, degree)
    if scaled.lower <= -1 or scaled.upper >= 1:
        return global_jet
    one_minus = Interval.singleton(1) - scaled ** 2
    if not one_minus.strictly_positive():
        return global_jet
    try:
        x_jet = (r_jet - center) / half_width
        one_minus_jet = IntervalJet.constant(1, degree) - x_jet * x_jet
        if not one_minus_jet.primal.strictly_positive():
            return global_jet
        peak = one_minus_jet.reciprocal()
        exponent = peak - IntervalJet.constant(1, degree)
        bump = _jet_exp_neg(exponent)
    except ZeroDivisionError as exc:
        raise SGBLTaylorEvaluationStop(
            "interval_inconclusive",
            "compact-bump jet reciprocal encountered a zero-containing denominator",
            {"obstruction": "zero_in_interval_reciprocal"},
        ) from exc
    return _intersect_jets(bump, global_jet)


def _family_field_jets(
    radius_jet: IntervalJet,
    spec: SGBLExactInitialSlice,
    field_order: int,
) -> dict[str, IntervalJet]:
    degree = _positive_int("field_order", field_order, allow_zero=True)
    bump_order = degree + 2
    radius_jet = IntervalJet.identity(radius_jet.primal, bump_order)
    bump = sgbl_compact_bump_interval_jet(
        radius_jet,
        center=spec.center,
        half_width=spec.half_width,
        order=bump_order,
    )
    sgbl_require_complete_jet(bump, bump_order)
    phi_amp = IntervalJet.constant(spec.phi_amplitude, degree)
    chi_amp = IntervalJet.constant(spec.chi_amplitude, degree)
    r_field = radius_jet.truncate(degree)
    bump_value = bump.truncate(degree)
    bump_r = sgbl_interval_jet_shift(bump, 1, degree)
    bump_rr = sgbl_interval_jet_shift(bump, 2, degree)
    try:
        inv_r = r_field.reciprocal()
    except ZeroDivisionError as exc:
        raise SGBLTaylorEvaluationStop(
            "interval_inconclusive",
            "family-field jet reciprocal encountered a zero-containing radius",
            {"obstruction": "zero_in_interval_reciprocal"},
        ) from exc
    return {
        "phi": phi_amp * bump_value,
        "phi_r": phi_amp * bump_r,
        "phi_rr": phi_amp * bump_rr,
        "phi_pi": IntervalJet.constant(0, degree),
        "phi_pi_r": IntervalJet.constant(0, degree),
        "chi_r": chi_amp * (bump_r * inv_r - bump_value * inv_r * inv_r),
        "chi_pi": chi_amp * bump_r * inv_r,
    }


def _count_jet(budget: _TaylorBudget, limit: int) -> None:
    budget.jet_evaluations += 1
    if budget.jet_evaluations > limit:
        raise SGBLTaylorEvaluationStop(
            "resource_limit",
            "Taylor-model exceeded the declared jet-evaluation budget",
            {
                "evaluations": budget.jet_evaluations,
                "limit": limit,
                "resource_reason": "max_jet_evaluations",
            },
        )


def _outward_jet(jet: IntervalJet, budget: _TaylorBudget, bit_limit: int) -> IntervalJet:
    rounded = tuple(_outward(value) for value in jet.derivatives)
    _observe(budget, bit_limit, *rounded)
    return IntervalJet(rounded)


def sgbl_affine_constraint_jets(
    radius: IntervalJet,
    radial_metric: IntervalJet,
    angular_extrinsic_curvature: IntervalJet,
    spec: SGBLExactInitialSlice,
    *,
    budget: _TaylorBudget,
    jet_limit: int,
    bit_limit: int,
    family: Mapping[str, IntervalJet] | None = None,
) -> dict[str, IntervalJet]:
    """Affine H/M jets and ``(lambda_r,k_r,C_r)`` by exact interval AD."""

    if type(spec) is not SGBLExactInitialSlice:
        raise TypeError("spec must be SGBLExactInitialSlice")
    _count_jet(budget, jet_limit)
    order = radial_metric.order
    if angular_extrinsic_curvature.order != order or radius.order < order:
        raise SGBLTaylorTrapStop(
            "omitted_derivative",
            "affine jets require matching lambda/k orders and a long enough r jet",
        )
    if not radius.primal.strictly_positive() or not radial_metric.primal.strictly_positive():
        raise SGBLTaylorEvaluationStop(
            "domain_error",
            "affine jets require positive r and lambda",
        )
    r_jet = radius.truncate(order)
    fields = family or _family_field_jets(radius, spec, order)
    try:
        mpl2 = IntervalJet.constant(spec.planck_mass_squared, order)
        alpha = IntervalJet.constant(spec.alpha_gb, order)
        two = IntervalJet.constant(2, order)
        eight = IntervalJet.constant(8, order)
        sixteen = IntervalJet.constant(16, order)
        lam = radial_metric
        k = angular_extrinsic_curvature
        a0 = IntervalJet.constant(-2, order) * (k ** 2)
        a_l = (lam ** 3 * r_jet).reciprocal()
        b = (k ** 2) + (IntervalJet.constant(1, order) - (lam ** 2).reciprocal()) / (
            r_jet ** 2
        )
        c0 = IntervalJet.constant(3, order) * k / (lam * r_jet)
        c_k = lam.reciprocal()
        x0 = fields["phi_rr"] / (lam ** 2) - two * k * fields["phi_pi"]
        x_l = -fields["phi_r"] / (lam ** 3)
        y = k * fields["phi_pi"] + fields["phi_r"] / (lam ** 2 * r_jet)
        z = (fields["phi_pi_r"] - two * k * fields["phi_r"]) / lam
        rho = (
            (fields["phi_pi"] ** 2 + fields["chi_pi"] ** 2) / two
            + (fields["phi_r"] ** 2 + fields["chi_r"] ** 2) / (two * (lam ** 2))
            + IntervalJet.constant(spec.scalar_mass ** 2, order)
            * (fields["phi"] ** 2)
            / two
            + IntervalJet.constant(spec.quartic_coupling, order)
            * (fields["phi"] ** 4)
            / IntervalJet.constant(4, order)
        )
        hamiltonian_constant = mpl2 * (two * a0 + b) - rho - eight * alpha * (
            b * x0 + two * a0 * y
        )
        hamiltonian_lambda_r = two * mpl2 * a_l - eight * alpha * (b * x_l + two * a_l * y)
        momentum_constant = (
            two * mpl2 * lam * c0
            - fields["phi_pi"] * fields["phi_r"]
            - fields["chi_pi"] * fields["chi_r"]
            - eight * alpha * lam * (b * z + two * c0 * y)
        )
        momentum_k_r = two * mpl2 * lam * c_k - sixteen * alpha * lam * c_k * y
    except ZeroDivisionError as exc:
        raise SGBLTaylorEvaluationStop(
            "interval_inconclusive",
            "affine jet reciprocal encountered a zero-containing denominator",
            {"obstruction": "zero_in_interval_reciprocal"},
        ) from exc
    if hamiltonian_lambda_r.primal.contains_zero() or momentum_k_r.primal.contains_zero():
        raise SGBLTaylorEvaluationStop(
            "interval_inconclusive",
            "affine constraint Jacobian diagonal contains zero",
            {
                "obstruction": "jacobian_diagonal_contains_zero",
                "hamiltonian_lambda_r": _interval_text(hamiltonian_lambda_r.primal),
                "momentum_k_r": _interval_text(momentum_k_r.primal),
            },
        )
    try:
        lambda_r = -hamiltonian_constant / hamiltonian_lambda_r
        k_r = -momentum_constant / momentum_k_r
        compactness_r = (
            two * r_jet * (k ** 2)
            + two * (r_jet ** 2) * k * k_r
            + two * lambda_r / (lam ** 3)
        )
    except ZeroDivisionError as exc:
        raise SGBLTaylorEvaluationStop(
            "interval_inconclusive",
            "constraint RHS jet reciprocal encountered a zero-containing denominator",
            {"obstruction": "zero_in_interval_reciprocal"},
        ) from exc
    hamiltonian_residual = hamiltonian_constant + hamiltonian_lambda_r * lambda_r
    momentum_residual = momentum_constant + momentum_k_r * k_r
    payload = {
        "hamiltonian_constant": hamiltonian_constant,
        "hamiltonian_lambda_r": hamiltonian_lambda_r,
        "momentum_constant": momentum_constant,
        "momentum_k_r": momentum_k_r,
        "lambda_r": lambda_r,
        "k_r": k_r,
        "compactness_r": compactness_r,
        "hamiltonian_residual": hamiltonian_residual,
        "momentum_residual": momentum_residual,
    }
    return {
        name: _outward_jet(jet, budget, bit_limit) for name, jet in payload.items()
    }


def _fill_state_jets(
    radius_primal: Interval,
    lambda0: Interval,
    k0: Interval,
    spec: SGBLExactInitialSlice,
    order: int,
    budget: _TaylorBudget,
    jet_limit: int,
    bit_limit: int,
) -> dict[str, IntervalJet]:
    """Triangular ODE jets: ``y^{(m+1)} = f^{(m)}`` along the constraint flow."""

    degree = _positive_int("order", order, allow_zero=True)
    bump_order = degree + 2
    r_jet = IntervalJet.identity(radius_primal, bump_order)
    family = _family_field_jets(r_jet, spec, degree)
    lambda_derivatives = [lambda0] + [ZERO] * (degree + 1)
    k_derivatives = [k0] + [ZERO] * (degree + 1)
    compactness_derivatives = [None] * (degree + 1)
    last: dict[str, IntervalJet] | None = None
    for index in range(degree + 1):
        rhs = sgbl_affine_constraint_jets(
            r_jet.truncate(index),
            IntervalJet(tuple(lambda_derivatives[: index + 1])),
            IntervalJet(tuple(k_derivatives[: index + 1])),
            spec,
            budget=budget,
            jet_limit=jet_limit,
            bit_limit=bit_limit,
            family={name: jet.truncate(index) for name, jet in family.items()},
        )
        lambda_derivatives[index + 1] = rhs["lambda_r"].derivatives[index]
        k_derivatives[index + 1] = rhs["k_r"].derivatives[index]
        compactness_derivatives[index] = rhs["compactness_r"].derivatives[index]
        last = rhs
    assert last is not None
    compactness_left = (
        Interval.singleton(1)
        + (radius_primal ** 2) * (k0 ** 2)
        - (lambda0 ** 2).reciprocal()
    )
    compactness_jet = IntervalJet(
        (_outward(compactness_left),) + tuple(compactness_derivatives)
    )
    return {
        "lambda": IntervalJet(tuple(lambda_derivatives[: degree + 1])),
        "k": IntervalJet(tuple(k_derivatives[: degree + 1])),
        "compactness": compactness_jet.truncate(degree),
        "lambda_remainder": IntervalJet(tuple(lambda_derivatives)),
        "k_remainder": IntervalJet(tuple(k_derivatives)),
        "compactness_remainder": compactness_jet,
        "hamiltonian_residual": last["hamiltonian_residual"],
        "momentum_residual": last["momentum_residual"],
        "hamiltonian_lambda_r": last["hamiltonian_lambda_r"],
        "momentum_k_r": last["momentum_k_r"],
    }


def sgbl_algebraic_compactness_interval(
    radius: Interval,
    radial_metric: Interval,
    angular_extrinsic_curvature: Interval,
) -> Interval:
    if not radius.strictly_positive() or not radial_metric.strictly_positive():
        raise SGBLTaylorEvaluationStop(
            "domain_error",
            "algebraic compactness requires positive r and lambda",
        )
    return (
        Interval.singleton(1)
        + (radius ** 2) * (angular_extrinsic_curvature ** 2)
        - (radial_metric ** 2).reciprocal()
    )


def _explicit_rhs_jet(
    radius: IntervalJet,
    state: IntervalJet,
    coefficients: Sequence[Interval],
) -> IntervalJet:
    """Jet of a polynomial right-hand side ``sum a_k r^k`` independent of state."""

    del state
    result = IntervalJet.constant(coefficients[0], radius.order)
    power = IntervalJet.constant(1, radius.order)
    for coefficient in coefficients[1:]:
        power = power * radius
        result = result + IntervalJet.constant(coefficient, radius.order) * power
    return result


def sgbl_taylor_polynomial_ode_control(
    *,
    order: int = DECLARED_TAYLOR_ORDER,
) -> tuple[Mapping[str, Any], ...]:
    """Exact remainder-zero reconstructions of polynomial IVPs ``y'=p'(r)``.

    For ``deg p <= order`` the Lagrange remainder derivative is identically
    zero, so the Taylor model reconstructs ``p`` exactly on any cell.
    """

    degree = _positive_int("order", order, allow_zero=True)
    records = []
    left = Q(0)
    right = Q(1, 2)
    domain = interval(left, right)
    for poly_degree in range(degree + 1):
        # y = r^d, y' = d r^{d-1}, y(0)=0 for d>=1 and y=1 for d=0.
        y0 = Interval.singleton(1 if poly_degree == 0 else 0)
        r_jet = IntervalJet.identity(Interval.singleton(left), degree + 1)
        rhs_coefficients = [ZERO] * max(poly_degree, 1)
        if poly_degree == 0:
            rhs_coefficients = [ZERO]
        else:
            rhs_coefficients = [ZERO] * poly_degree
            rhs_coefficients[poly_degree - 1] = Interval.singleton(poly_degree)
        derivatives = [y0]
        for index in range(degree + 1):
            state = IntervalJet(tuple(derivatives[: index + 1]))
            rhs = _explicit_rhs_jet(r_jet.truncate(index), state, rhs_coefficients)
            derivatives.append(rhs.derivatives[index])
        jet = IntervalJet(tuple(derivatives[: degree + 1]))
        remainder_derivative = derivatives[degree + 1]
        model = sgbl_taylor_model_from_jet(jet, remainder_derivative, left, domain)
        honest = sgbl_require_honest_remainder(
            model.remainder_on_domain(),
            remainder_derivative,
            model.displacement(),
            degree,
        )
        exact_right = Interval.singleton(right ** poly_degree if poly_degree else 1)
        reconstructed = model.at(right)
        records.append(
            MappingProxyType(
                {
                    "polynomial_degree": poly_degree,
                    "taylor_order": degree,
                    "remainder_derivative_is_zero": remainder_derivative == ZERO,
                    "honest_remainder_is_zero": honest == ZERO,
                    "contains_exact_right_endpoint": (
                        reconstructed.lower <= exact_right.lower
                        and exact_right.upper <= reconstructed.upper
                    ),
                    "exact_at_right_endpoint": reconstructed == exact_right,
                    "range": _interval_text(model.range()),
                }
            )
        )
        if poly_degree <= degree and remainder_derivative != ZERO:
            raise SGBLTaylorTrapStop(
                "wrong_remainder",
                "polynomial ODE remainder derivative must vanish at sufficient order",
                {"degree": poly_degree, "order": degree},
            )
        if poly_degree <= degree and reconstructed != exact_right:
            raise SGBLTaylorTrapStop(
                "wrong_remainder",
                "polynomial ODE Taylor model must reconstruct the exact endpoint",
                {"degree": poly_degree, "got": _interval_text(reconstructed)},
            )
    return tuple(records)


def sgbl_taylor_known_linear_ode_control(
    *,
    order: int = DECLARED_TAYLOR_ORDER,
) -> Mapping[str, Any]:
    """Validated Taylor model of ``y'=2r``, ``y(0)=0``, whose solution is ``r^2``."""

    degree = _positive_int("order", order, allow_zero=True)
    left = Q(0)
    right = Q(1, 4)
    domain = interval(left, right)
    r_jet = IntervalJet.identity(Interval.singleton(left), degree + 1)
    derivatives = [ZERO]
    for index in range(degree + 1):
        rhs = IntervalJet.constant(2, index) * r_jet.truncate(index)
        derivatives.append(rhs.derivatives[index])
    jet = IntervalJet(tuple(derivatives[: degree + 1]))
    remainder_derivative = derivatives[degree + 1]
    model = sgbl_taylor_model_from_jet(jet, remainder_derivative, left, domain)
    sgbl_require_honest_remainder(
        model.remainder_on_domain(),
        remainder_derivative,
        model.displacement(),
        degree,
    )
    reconstructed = model.at(right)
    exact = Interval.singleton(right ** 2)
    return MappingProxyType(
        {
            "ode": "y_prime_equals_two_r",
            "solution": "r_squared",
            "remainder_derivative_is_zero": remainder_derivative == ZERO,
            "contains_exact_right_endpoint": (
                reconstructed.lower <= exact.lower and exact.upper <= reconstructed.upper
            ),
            "exact_at_right_endpoint": reconstructed == exact,
            "reconstructed": _interval_text(reconstructed),
            "exact": _interval_text(exact),
        }
    )


@dataclass(frozen=True, slots=True, kw_only=True)
class SGBLTaylorResourceLevel:
    """One prospectively declared Taylor-model resource level."""

    name: str
    taylor_order: int
    max_bisection_depth: int
    max_cells: int
    max_jet_evaluations: int
    max_rational_bit_length: int
    base_cells: int = DECLARED_BASE_CELLS
    max_picard_iterations: int = DECLARED_MAX_TAYLOR_PICARD_ITERATIONS
    lambda_domain: Interval = DECLARED_LAMBDA_DOMAIN
    k_domain: Interval = DECLARED_K_DOMAIN
    compactness_domain: Interval = DECLARED_C_DOMAIN
    angular_momentum_domain: Interval = DECLARED_J_DOMAIN

    def __post_init__(self) -> None:
        object.__setattr__(self, "name", str(self.name))
        object.__setattr__(
            self, "taylor_order", _positive_int("taylor_order", self.taylor_order)
        )
        object.__setattr__(
            self,
            "max_bisection_depth",
            _positive_int("max_bisection_depth", self.max_bisection_depth, allow_zero=True),
        )
        object.__setattr__(self, "max_cells", _positive_int("max_cells", self.max_cells))
        object.__setattr__(
            self,
            "max_jet_evaluations",
            _positive_int("max_jet_evaluations", self.max_jet_evaluations),
        )
        object.__setattr__(
            self,
            "max_rational_bit_length",
            _positive_int("max_rational_bit_length", self.max_rational_bit_length),
        )
        object.__setattr__(self, "base_cells", _positive_int("base_cells", self.base_cells))
        object.__setattr__(
            self,
            "max_picard_iterations",
            _positive_int("max_picard_iterations", self.max_picard_iterations),
        )
        if self.taylor_order != DECLARED_TAYLOR_ORDER:
            raise ValueError("Taylor-model ladder must keep the frozen Taylor order")
        if self.base_cells != DECLARED_BASE_CELLS:
            raise ValueError("Taylor-model ladder must keep the declared 16 base cells")
        if type(self.lambda_domain) is not Interval or type(self.k_domain) is not Interval:
            raise TypeError("ladder domains must be exact intervals")
        if type(self.compactness_domain) is not Interval:
            raise TypeError("ladder compactness domain must be an exact interval")
        if (
            self.lambda_domain != DECLARED_LAMBDA_DOMAIN
            or self.k_domain != DECLARED_K_DOMAIN
            or self.compactness_domain != DECLARED_C_DOMAIN
        ):
            raise ValueError("Taylor-model ladder must keep the declared lambda/k/C domains")
        if self.compactness_domain.upper <= DECLARED_COMPACTNESS_STRICT_UPPER:
            raise ValueError("compactness resource domain must not assume C<=1")
        if not (self.compactness_domain.lower < 0 < self.compactness_domain.upper):
            raise ValueError("compactness resource domain must contain a neighborhood of C=0")

    def as_contract_mapping(self) -> dict[str, int | str | tuple[str, str]]:
        return {
            "name": self.name,
            "taylor_order": self.taylor_order,
            "max_bisection_depth": self.max_bisection_depth,
            "max_cells": self.max_cells,
            "max_jet_evaluations": self.max_jet_evaluations,
            "max_rational_bit_length": self.max_rational_bit_length,
            "base_cells": self.base_cells,
            "max_picard_iterations": self.max_picard_iterations,
            "lambda_domain": _interval_text(self.lambda_domain),
            "k_domain": _interval_text(self.k_domain),
            "compactness_domain": _interval_text(self.compactness_domain),
        }


def sgbl_declared_taylor_cell_cap(max_bisection_depth: int) -> int:
    """Worst-case partition calls on the 16-base-cell binary tree."""

    depth = _positive_int("max_bisection_depth", max_bisection_depth, allow_zero=True)
    return DECLARED_BASE_CELLS * (2 ** (depth + 1) - 1)


def sgbl_declared_taylor_jet_cap(max_bisection_depth: int) -> int:
    """Worst-case affine jet evaluations: seed, Picard iterates, remainder."""

    return sgbl_declared_taylor_cell_cap(max_bisection_depth) * DECLARED_JET_CALLS_PER_PARTITION * (
        DECLARED_TAYLOR_ORDER + 2
    )


def _declare_taylor_resource_ladder() -> tuple[SGBLTaylorResourceLevel, ...]:
    levels = []
    for depth in DECLARED_TAYLOR_LEVEL_DEPTHS:
        levels.append(
            SGBLTaylorResourceLevel(
                name=f"order_{DECLARED_TAYLOR_ORDER}_depth_{depth}",
                taylor_order=DECLARED_TAYLOR_ORDER,
                max_bisection_depth=depth,
                max_cells=sgbl_declared_taylor_cell_cap(depth),
                max_jet_evaluations=sgbl_declared_taylor_jet_cap(depth),
                max_rational_bit_length=DECLARED_MAX_RATIONAL_BITS,
            )
        )
    return tuple(levels)


DECLARED_TAYLOR_RESOURCE_LADDER = _declare_taylor_resource_ladder()


def _require_declared_ladder(ladder: object) -> tuple[SGBLTaylorResourceLevel, ...]:
    if type(ladder) is not tuple:
        raise SGBLTaylorTrapStop(
            "changed_policy",
            "Taylor resource ladder must be the frozen declared tuple",
        )
    if not ladder:
        raise SGBLTaylorTrapStop(
            "skipped_level",
            "declared Taylor ladder levels were skipped",
        )
    if not all(type(level) is SGBLTaylorResourceLevel for level in ladder):
        raise SGBLTaylorTrapStop(
            "changed_policy",
            "Taylor ladder levels must be SGBLTaylorResourceLevel",
        )
    depths = tuple(level.max_bisection_depth for level in ladder)
    declared_set = set(DECLARED_TAYLOR_LEVEL_DEPTHS)
    if len(ladder) < len(DECLARED_TAYLOR_LEVEL_DEPTHS) and set(depths) <= declared_set:
        raise SGBLTaylorTrapStop(
            "skipped_level",
            "a declared Taylor ladder level was skipped",
            {"depths": depths},
        )
    if (
        len(ladder) == len(DECLARED_TAYLOR_LEVEL_DEPTHS)
        and set(depths) == declared_set
        and depths != DECLARED_TAYLOR_LEVEL_DEPTHS
    ):
        raise SGBLTaylorTrapStop(
            "reordered_levels",
            "declared Taylor ladder levels were reordered",
            {"depths": depths},
        )
    if ladder != DECLARED_TAYLOR_RESOURCE_LADDER:
        raise SGBLTaylorTrapStop(
            "changed_policy",
            "Taylor resource ladder is not the frozen declaration",
            {"depths": depths},
        )
    return ladder


def _require_predecessor_hashes(
    trap_refinement: object,
    trap_refinement2: object,
) -> tuple[str, str]:
    if type(trap_refinement) is not str or type(trap_refinement2) is not str:
        raise SGBLTaylorTrapStop(
            "wrong_predecessor",
            "predecessor contract hashes must be the frozen depth-ladder digests",
        )
    if trap_refinement != PREDECESSOR_TRAP_REFINEMENT_CONTRACT_SHA256:
        raise SGBLTaylorTrapStop(
            "wrong_predecessor",
            "Taylor instrument does not bind the preserved depth-8/10/12 hash",
            {
                "supplied": trap_refinement,
                "required": PREDECESSOR_TRAP_REFINEMENT_CONTRACT_SHA256,
            },
        )
    if trap_refinement2 != PREDECESSOR_TRAP_REFINEMENT2_CONTRACT_SHA256:
        raise SGBLTaylorTrapStop(
            "wrong_predecessor",
            "Taylor instrument does not bind the preserved depth-14/16/18 hash",
            {
                "supplied": trap_refinement2,
                "required": PREDECESSOR_TRAP_REFINEMENT2_CONTRACT_SHA256,
            },
        )
    return trap_refinement, trap_refinement2


@dataclass(frozen=True, slots=True, kw_only=True)
class SGBLTaylorCellEnclosure:
    """Validated Taylor-model graph enclosure on one radial cell."""

    radius: Interval
    lambda_box: Interval
    k_box: Interval
    compactness_box: Interval
    lambda_left: Interval
    k_left: Interval
    compactness_left: Interval
    lambda_right: Interval
    k_right: Interval
    compactness_right: Interval
    compactness: Interval
    direct_compactness: Interval
    propagated_compactness: Interval
    invariant_residual: Interval
    taylor_order: int
    depth: int
    taylor_strict_self_map: bool
    residual_contains_origin: bool
    compactness_invariant_contains_zero: bool
    correlated_strictly_tighter: bool
    chart_kind: str
    chart_coordinates: bool

    @property
    def compactness_upper(self) -> Fraction:
        return self.compactness.upper


def sgbl_validate_taylor_cell_inventory(
    cells: tuple[SGBLTaylorCellEnclosure, ...],
    *,
    origin: Fraction,
    terminus: Fraction | None = None,
) -> tuple[Fraction, Fraction]:
    """Require a strictly ordered, nonoverlapping, gap-free radial inventory."""

    if not cells:
        if terminus is not None and terminus != origin:
            raise ValueError("empty Taylor inventory cannot cover a nonempty interval")
        return origin, origin
    if cells[0].radius.lower != origin:
        raise ValueError("Taylor inventory does not start at the requested left boundary")
    previous_right = origin
    for cell in cells:
        if cell.radius.upper <= cell.radius.lower:
            raise ValueError("Taylor cell radius is not a positive-length interval")
        if cell.radius.lower < previous_right:
            raise ValueError("Taylor inventory contains overlapping cells")
        if cell.radius.lower > previous_right:
            raise ValueError("Taylor inventory contains a radial gap")
        previous_right = cell.radius.upper
    if terminus is not None and previous_right != terminus:
        raise ValueError("Taylor inventory does not reach the required right boundary")
    return origin, previous_right


def sgbl_taylor_inventory_span(
    cells: tuple[SGBLTaylorCellEnclosure, ...],
    spec: SGBLExactInitialSlice,
) -> tuple[Fraction, Fraction]:
    if type(spec) is not SGBLExactInitialSlice:
        raise TypeError("spec must be SGBLExactInitialSlice")
    if any(type(cell) is not SGBLTaylorCellEnclosure for cell in cells):
        raise TypeError("cells must be SGBLTaylorCellEnclosure")
    try:
        return sgbl_validate_taylor_cell_inventory(cells, origin=spec.support_minimum)
    except ValueError as exc:
        raise SGBLTaylorTrapStop(
            "gapped_inventory",
            "Taylor ladder refuses a dropped or gapped cell inventory",
            {"detail": str(exc)},
        ) from exc


def sgbl_taylor_cells_cover_support(
    cells: tuple[SGBLTaylorCellEnclosure, ...],
    spec: SGBLExactInitialSlice,
) -> bool:
    if not cells:
        return False
    _left, right = sgbl_taylor_inventory_span(cells, spec)
    return right == spec.support_maximum


def sgbl_taylor_shared_coverage_monotone(
    coarser: tuple[SGBLTaylorCellEnclosure, ...],
    finer: tuple[SGBLTaylorCellEnclosure, ...],
) -> bool:
    for left in coarser:
        for right in finer:
            if _intersect(left.radius, right.radius) is None:
                continue
            if (
                _intersect(left.lambda_box, right.lambda_box) is None
                or _intersect(left.k_box, right.k_box) is None
                or _intersect(left.compactness, right.compactness) is None
            ):
                return False
            if right.radius.subset_of(left.radius) and not (
                right.lambda_box.subset_of(left.lambda_box)
                and right.k_box.subset_of(left.k_box)
                and right.compactness.subset_of(left.compactness)
            ):
                return False
    return True


def _order_zero_rhs(
    radius: Interval,
    radial_metric: Interval,
    angular_extrinsic_curvature: Interval,
    spec: SGBLExactInitialSlice,
    budget: _TaylorBudget,
    jet_limit: int,
    bit_limit: int,
) -> dict[str, Interval]:
    payload = sgbl_affine_constraint_jets(
        IntervalJet.identity(radius, 0),
        IntervalJet.constant(radial_metric, 0),
        IntervalJet.constant(angular_extrinsic_curvature, 0),
        spec,
        budget=budget,
        jet_limit=jet_limit,
        bit_limit=bit_limit,
    )
    budget.compactness_jet_evaluations += 1
    return {name: jet.primal for name, jet in payload.items()}


def _try_taylor_cell(
    r0: Fraction,
    r1: Fraction,
    lambda_left: Interval,
    k_left: Interval,
    compactness_left: Interval,
    spec: SGBLExactInitialSlice,
    level: SGBLTaylorResourceLevel,
    budget: _TaylorBudget,
    depth: int,
) -> SGBLTaylorCellEnclosure | str:
    radius = interval(r0, r1)
    width = r1 - r0
    if width <= 0:
        return "domain_error"
    step = interval(0, width)
    try:
        seed = _order_zero_rhs(
            radius,
            lambda_left,
            k_left,
            spec,
            budget,
            level.max_jet_evaluations,
            level.max_rational_bit_length,
        )
    except SGBLTaylorEvaluationStop as exc:
        if exc.reason in {"resource_limit", "domain_error"}:
            raise
        return str(exc.payload.get("obstruction") or "zero_in_interval_reciprocal")
    lambda_radius = DECLARED_EULER_INFLATION * width * seed["lambda_r"].abs_upper()
    k_radius = DECLARED_EULER_INFLATION * width * seed["k_r"].abs_upper()
    compactness_radius = DECLARED_EULER_INFLATION * width * seed["compactness_r"].abs_upper()
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
    graph_lambda = _intersect(seed_lambda, level.lambda_domain)
    graph_k = _intersect(seed_k, level.k_domain)
    graph_compactness = _intersect(seed_compactness, level.compactness_domain)
    if graph_lambda is None or graph_k is None or graph_compactness is None:
        return "physical_domain_escape"
    graph_lambda, graph_k, graph_compactness = _outward_all(
        budget,
        level.max_rational_bit_length,
        graph_lambda,
        graph_k,
        graph_compactness,
    )
    payload: dict[str, Interval] | None = None
    image_lambda = seed_lambda
    image_k = seed_k
    image_compactness = seed_compactness
    for _ in range(level.max_picard_iterations):
        try:
            payload = _order_zero_rhs(
                radius,
                graph_lambda,
                graph_k,
                spec,
                budget,
                level.max_jet_evaluations,
                level.max_rational_bit_length,
            )
        except SGBLTaylorEvaluationStop as exc:
            if exc.reason in {"resource_limit", "domain_error"}:
                raise
            return str(exc.payload.get("obstruction") or "zero_in_interval_reciprocal")
        image_lambda = _outward(lambda_left + step * payload["lambda_r"])
        image_k = _outward(k_left + step * payload["k_r"])
        image_compactness = _outward(compactness_left + step * payload["compactness_r"])
        if (
            image_lambda.strictly_inside(graph_lambda)
            and image_k.strictly_inside(graph_k)
            and image_compactness.strictly_inside(graph_compactness)
        ):
            break
        next_lambda = _intersect(image_lambda, graph_lambda)
        next_k = _intersect(image_k, graph_k)
        next_compactness = _intersect(image_compactness, graph_compactness)
        if next_lambda is None or next_k is None or next_compactness is None:
            return "physical_domain_escape"
        next_lambda, next_k, next_compactness = _outward_all(
            budget,
            level.max_rational_bit_length,
            next_lambda,
            next_k,
            next_compactness,
        )
        if (
            next_lambda == graph_lambda
            and next_k == graph_k
            and next_compactness == graph_compactness
        ):
            return "taylor_self_map_failed"
        graph_lambda, graph_k, graph_compactness = next_lambda, next_k, next_compactness
    else:
        return "taylor_self_map_failed"
    if payload is None:
        return "taylor_self_map_failed"
    if not (
        image_lambda.strictly_inside(graph_lambda)
        and image_k.strictly_inside(graph_k)
        and image_compactness.strictly_inside(graph_compactness)
    ):
        return "taylor_self_map_failed"
    try:
        point_jets = _fill_state_jets(
            Interval.singleton(r0),
            lambda_left,
            k_left,
            spec,
            level.taylor_order,
            budget,
            level.max_jet_evaluations,
            level.max_rational_bit_length,
        )
        box_jets = _fill_state_jets(
            radius,
            graph_lambda,
            graph_k,
            spec,
            level.taylor_order,
            budget,
            level.max_jet_evaluations,
            level.max_rational_bit_length,
        )
    except SGBLTaylorEvaluationStop as exc:
        if exc.reason in {"resource_limit", "domain_error"}:
            raise
        return str(exc.payload.get("obstruction") or "zero_in_interval_reciprocal")
    except ZeroDivisionError:
        return "zero_in_interval_reciprocal"
    lambda_model = sgbl_taylor_model_from_jet(
        point_jets["lambda"],
        box_jets["lambda_remainder"].derivatives[level.taylor_order + 1],
        r0,
        radius,
    )
    k_model = sgbl_taylor_model_from_jet(
        point_jets["k"],
        box_jets["k_remainder"].derivatives[level.taylor_order + 1],
        r0,
        radius,
    )
    compactness_model = sgbl_taylor_model_from_jet(
        point_jets["compactness"],
        box_jets["compactness_remainder"].derivatives[level.taylor_order + 1],
        r0,
        radius,
    )
    sgbl_require_honest_remainder(
        lambda_model.remainder_on_domain(),
        lambda_model.remainder_derivative,
        lambda_model.displacement(),
        level.taylor_order,
    )
    sgbl_require_honest_remainder(
        k_model.remainder_on_domain(),
        k_model.remainder_derivative,
        k_model.displacement(),
        level.taylor_order,
    )
    sgbl_require_honest_remainder(
        compactness_model.remainder_on_domain(),
        compactness_model.remainder_derivative,
        compactness_model.displacement(),
        level.taylor_order,
    )
    lambda_tm = lambda_model.range()
    k_tm = k_model.range()
    compactness_tm = compactness_model.range()
    correlated_lambda = _intersect(lambda_tm, graph_lambda)
    correlated_k = _intersect(k_tm, graph_k)
    if correlated_lambda is None or correlated_k is None:
        return "taylor_remainder_inconsistent"
    try:
        direct = _outward(
            sgbl_algebraic_compactness_interval(radius, correlated_lambda, correlated_k)
        )
    except SGBLTaylorEvaluationStop as exc:
        if exc.reason in {"resource_limit", "domain_error"}:
            raise
        return str(exc.payload.get("obstruction") or "zero_in_interval_reciprocal")
    except ZeroDivisionError:
        return "zero_in_interval_reciprocal"
    propagated = _outward(image_compactness)
    correlated = _intersect(direct, propagated)
    if correlated is None:
        return "compactness_invariant_misses_origin"
    correlated_tm = _intersect(correlated, compactness_tm)
    if correlated_tm is None:
        return "taylor_remainder_inconsistent"
    correlated = _outward(correlated_tm)
    invariant = _outward(correlated - direct)
    if not invariant.contains_zero():
        return "compactness_invariant_misses_origin"
    residual_ok = payload["hamiltonian_residual"].contains_zero() and payload[
        "momentum_residual"
    ].contains_zero()
    if not residual_ok:
        return "residual_enclosure_misses_origin"
    right_lambda = _outward(lambda_model.at(r1))
    right_k = _outward(k_model.at(r1))
    right_from_euler = _outward(lambda_left + Interval.singleton(width) * payload["lambda_r"])
    right_k_euler = _outward(k_left + Interval.singleton(width) * payload["k_r"])
    right_lambda_join = _intersect(right_lambda, right_from_euler)
    right_k_join = _intersect(right_k, right_k_euler)
    if right_lambda_join is None or right_k_join is None:
        return "taylor_remainder_inconsistent"
    right_lambda = _outward(right_lambda_join)
    right_k = _outward(right_k_join)
    right_compactness_tm = _outward(compactness_model.at(r1))
    right_compactness_euler = _outward(
        compactness_left + Interval.singleton(width) * payload["compactness_r"]
    )
    try:
        right_direct = _outward(
            sgbl_algebraic_compactness_interval(Interval.singleton(r1), right_lambda, right_k)
        )
    except (SGBLTaylorEvaluationStop, ZeroDivisionError):
        return "zero_in_interval_reciprocal"
    right_compactness = _intersect(right_compactness_tm, right_compactness_euler)
    if right_compactness is None:
        return "taylor_remainder_inconsistent"
    right_compactness = _intersect(right_compactness, right_direct)
    if right_compactness is None:
        return "compactness_invariant_misses_origin"
    right_compactness = _outward(right_compactness)
    tighter = correlated.width() < direct.width()
    if tighter:
        budget.correlated_strictly_tighter = True
    _observe(
        budget,
        level.max_rational_bit_length,
        correlated,
        direct,
        propagated,
        invariant,
        right_lambda,
        right_k,
        right_compactness,
        correlated_lambda,
        correlated_k,
    )
    return SGBLTaylorCellEnclosure(
        radius=radius,
        lambda_box=correlated_lambda,
        k_box=correlated_k,
        compactness_box=graph_compactness,
        lambda_left=lambda_left,
        k_left=k_left,
        compactness_left=compactness_left,
        lambda_right=right_lambda,
        k_right=right_k,
        compactness_right=right_compactness,
        compactness=correlated,
        direct_compactness=direct,
        propagated_compactness=propagated,
        invariant_residual=invariant,
        taylor_order=level.taylor_order,
        depth=depth,
        taylor_strict_self_map=True,
        residual_contains_origin=True,
        compactness_invariant_contains_zero=True,
        correlated_strictly_tighter=tighter,
        chart_kind=PRODUCT_BOX_CHART_KIND,
        chart_coordinates=False,
    )


def _combine_partition(
    origin: Fraction,
    terminus: Fraction,
    left_cells: tuple[SGBLTaylorCellEnclosure, ...],
    right_cells: tuple[SGBLTaylorCellEnclosure, ...],
    right_reason: str | None,
) -> tuple[tuple[SGBLTaylorCellEnclosure, ...], str | None]:
    combined = left_cells + right_cells
    if right_reason is None:
        sgbl_validate_taylor_cell_inventory(combined, origin=origin, terminus=terminus)
    else:
        sgbl_validate_taylor_cell_inventory(combined, origin=origin)
    return combined, right_reason


def _partition_cell(
    r0: Fraction,
    r1: Fraction,
    lambda_left: Interval,
    k_left: Interval,
    compactness_left: Interval,
    spec: SGBLExactInitialSlice,
    level: SGBLTaylorResourceLevel,
    budget: _TaylorBudget,
    depth: int,
) -> tuple[tuple[SGBLTaylorCellEnclosure, ...], str | None]:
    budget.cells += 1
    if budget.cells > level.max_cells:
        raise SGBLTaylorEvaluationStop(
            "resource_limit",
            "Taylor-model exceeded the declared cell cap",
            {
                "cells": budget.cells,
                "limit": level.max_cells,
                "resource_reason": "max_cells",
            },
        )
    outcome = _try_taylor_cell(
        r0,
        r1,
        lambda_left,
        k_left,
        compactness_left,
        spec,
        level,
        budget,
        depth,
    )
    if type(outcome) is SGBLTaylorCellEnclosure:
        if outcome.compactness.upper < DECLARED_COMPACTNESS_STRICT_UPPER:
            return (outcome,), None
        if depth >= level.max_bisection_depth:
            return (outcome,), "compactness_enclosure_not_below_one"
        reason: str | None = "compactness_enclosure_not_below_one"
    else:
        reason = outcome
        if depth >= level.max_bisection_depth:
            return (), reason
    mid = (r0 + r1) / 2
    left_cells, left_reason = _partition_cell(
        r0,
        mid,
        lambda_left,
        k_left,
        compactness_left,
        spec,
        level,
        budget,
        depth + 1,
    )
    if left_reason is not None:
        sgbl_validate_taylor_cell_inventory(left_cells, origin=r0)
        return left_cells, left_reason
    sgbl_validate_taylor_cell_inventory(left_cells, origin=r0, terminus=mid)
    last = left_cells[-1]
    right_cells, right_reason = _partition_cell(
        mid,
        r1,
        last.lambda_right,
        last.k_right,
        last.compactness_right,
        spec,
        level,
        budget,
        depth + 1,
    )
    return _combine_partition(r0, r1, left_cells, right_cells, right_reason)


@dataclass(frozen=True, slots=True, kw_only=True)
class SGBLTaylorODERecord:
    """Taylor-model inventory of the ``(lambda,k,C)`` constraint graph."""

    spec: SGBLExactInitialSlice
    level: SGBLTaylorResourceLevel
    cells: tuple[SGBLTaylorCellEnclosure, ...]
    jet_evaluations: int
    compactness_jet_evaluations: int
    max_rational_bit_length_observed: int
    correlated_strictly_tighter_than_direct: bool
    chart_kind: str
    chart_coordinates: bool
    obstruction: str | None

    @property
    def tiles_compact_support(self) -> bool:
        if not self.cells:
            return False
        return (
            self.cells[0].radius.lower == self.spec.support_minimum
            and self.cells[-1].radius.upper == self.spec.support_maximum
            and self.obstruction is None
        )

    @property
    def covers_complete_support(self) -> bool:
        return self.tiles_compact_support

    @property
    def support_compactness(self) -> Interval:
        if not self.cells:
            return interval(0, 2)
        lower = min(cell.compactness.lower for cell in self.cells)
        upper = max(cell.compactness.upper for cell in self.cells)
        return Interval(lower, upper)

    @property
    def lambda_end(self) -> Interval:
        if not self.cells:
            return Interval.singleton(1)
        return self.cells[-1].lambda_right

    @property
    def k_end(self) -> Interval:
        if not self.cells:
            return Interval.singleton(0)
        return self.cells[-1].k_right

    @property
    def compactness_end(self) -> Interval:
        if not self.cells:
            return Interval.singleton(0)
        return self.cells[-1].compactness_right


def sgbl_validated_taylor_lambda_k_constraint_ode(
    spec: SGBLExactInitialSlice,
    level: SGBLTaylorResourceLevel,
) -> SGBLTaylorODERecord:
    """Run the frozen Taylor-model policy on one resource level."""

    if type(spec) is not SGBLExactInitialSlice:
        raise TypeError("spec must be SGBLExactInitialSlice")
    if type(level) is not SGBLTaylorResourceLevel:
        raise TypeError("level must be SGBLTaylorResourceLevel")
    start = spec.support_minimum
    end = spec.support_maximum
    width = (end - start) / level.base_cells
    budget = _TaylorBudget()
    cells: list[SGBLTaylorCellEnclosure] = []
    lambda_left = Interval.singleton(1)
    k_left = Interval.singleton(0)
    compactness_left = Interval.singleton(0)
    obstruction: str | None = None
    for index in range(level.base_cells):
        left = start + index * width
        right = end if index == level.base_cells - 1 else start + (index + 1) * width
        piece, obstruction = _partition_cell(
            left,
            right,
            lambda_left,
            k_left,
            compactness_left,
            spec,
            level,
            budget,
            0,
        )
        cells.extend(piece)
        if obstruction is not None:
            break
        lambda_left = cells[-1].lambda_right
        k_left = cells[-1].k_right
        compactness_left = cells[-1].compactness_right
    return SGBLTaylorODERecord(
        spec=spec,
        level=level,
        cells=tuple(cells),
        jet_evaluations=budget.jet_evaluations,
        compactness_jet_evaluations=budget.compactness_jet_evaluations,
        max_rational_bit_length_observed=budget.max_bits,
        correlated_strictly_tighter_than_direct=budget.correlated_strictly_tighter,
        chart_kind=PRODUCT_BOX_CHART_KIND,
        chart_coordinates=False,
        obstruction=obstruction,
    )


def _resource_reason_from_payload(payload: Mapping[str, Any]) -> str:
    if payload.get("resource_reason") in RESOURCE_REASONS:
        return str(payload["resource_reason"])
    if "evaluations" in payload:
        return "max_jet_evaluations"
    if "cells" in payload:
        return "max_cells"
    if "observed" in payload:
        return "max_rational_bit_length"
    return "max_jet_evaluations"


def _support_compactness_upper(ode: SGBLTaylorODERecord | None) -> Fraction:
    if ode is None or not ode.cells:
        return Q(2)
    return ode.support_compactness.upper


def _level_proved(
    *,
    classification: str,
    obstruction: str | None,
    complete_support_coverage: bool,
    compactness_upper: Fraction,
    exterior_compactness_upper: Fraction,
    compactness_margin: Fraction,
    ode: SGBLTaylorODERecord | None,
    resource_reason: str | None,
) -> bool:
    if classification != "continuous_no_initial_trapped_sphere":
        return False
    if (
        obstruction is not None
        or resource_reason is not None
        or not complete_support_coverage
        or compactness_upper >= DECLARED_COMPACTNESS_STRICT_UPPER
        or exterior_compactness_upper >= DECLARED_COMPACTNESS_STRICT_UPPER
        or compactness_margin <= 0
        or ode is None
        or not ode.cells
        or ode.chart_kind != PRODUCT_BOX_CHART_KIND
        or ode.chart_coordinates
        or ode.obstruction is not None
        or not ode.tiles_compact_support
        or not ode.covers_complete_support
    ):
        return False
    return all(
        cell.taylor_strict_self_map
        and cell.residual_contains_origin
        and cell.compactness_invariant_contains_zero
        and cell.chart_kind == PRODUCT_BOX_CHART_KIND
        and not cell.chart_coordinates
        and cell.compactness.upper < DECLARED_COMPACTNESS_STRICT_UPPER
        for cell in ode.cells
    )


def _exterior_upper(
    spec: SGBLExactInitialSlice,
    ode: SGBLTaylorODERecord,
) -> tuple[Fraction, str | None]:
    try:
        end = sgbl_algebraic_compactness_interval(
            Interval.singleton(spec.support_maximum),
            ode.lambda_end,
            ode.k_end,
        )
        correlated = _intersect(end, ode.compactness_end)
        if correlated is None:
            return max(ode.compactness_end.upper, Q(0)), "compactness_invariant_misses_origin"
        return max(correlated.upper, Q(0)), None
    except SGBLTaylorEvaluationStop as exc:
        if exc.reason in {"resource_limit", "domain_error"}:
            raise
        return max(ode.compactness_end.upper, Q(0)), str(
            exc.payload.get("obstruction") or "compactness_invariant_misses_origin"
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class SGBLTaylorLevelRecord:
    """Outcome of one prospectively declared Taylor-model resource level."""

    level: SGBLTaylorResourceLevel
    classification: str
    obstruction: str | None
    complete_support_coverage: bool
    compactness_upper: Fraction
    exterior_compactness_upper: Fraction
    compactness_margin: Fraction
    coverage_left: Fraction
    coverage_right: Fraction
    ode_cell_count: int
    jet_evaluations: int
    compactness_jet_evaluations: int
    max_rational_bit_length_observed: int
    taylor_strict_inclusions: int
    continuous_no_initial_trapped_sphere: bool
    resource_reason: str | None
    missing_theorem: str | None
    diagnostic_wall_seconds: float
    ode: SGBLTaylorODERecord | None = None

    def __post_init__(self) -> None:
        if type(self.level) is not SGBLTaylorResourceLevel:
            raise TypeError("level must be SGBLTaylorResourceLevel")
        if self.classification not in LEVEL_CLASSIFICATIONS:
            raise ValueError("unknown Taylor ladder classification")
        if type(self.compactness_upper) is not Fraction:
            raise TypeError("compactness_upper must be a Fraction")
        if type(self.exterior_compactness_upper) is not Fraction:
            raise TypeError("exterior_compactness_upper must be a Fraction")
        if type(self.compactness_margin) is not Fraction:
            raise TypeError("compactness_margin must be a Fraction")
        if type(self.coverage_left) is not Fraction or type(self.coverage_right) is not Fraction:
            raise TypeError("coverage endpoints must be Fractions")
        if self.compactness_margin != (
            1 - max(self.compactness_upper, self.exterior_compactness_upper)
        ):
            raise ValueError("compactness_margin must be the exact C<1 remainder")
        if self.resource_reason is not None and self.resource_reason not in RESOURCE_REASONS:
            raise ValueError("unknown resource reason")
        if self.classification == "resource_limit" and self.resource_reason is None:
            raise ValueError("resource_limit requires a typed resource reason")
        if self.classification != "resource_limit" and self.resource_reason is not None:
            raise ValueError("resource reasons cannot label a non-resource classification")
        if self.obstruction is not None and self.obstruction not in ODE_INCONCLUSIVE_REASONS:
            raise ValueError("unknown Taylor obstruction")
        proved = _level_proved(
            classification=self.classification,
            obstruction=self.obstruction,
            complete_support_coverage=self.complete_support_coverage,
            compactness_upper=self.compactness_upper,
            exterior_compactness_upper=self.exterior_compactness_upper,
            compactness_margin=self.compactness_margin,
            ode=self.ode,
            resource_reason=self.resource_reason,
        )
        if self.continuous_no_initial_trapped_sphere != proved:
            raise ValueError("continuous-no-trap bit does not match the Taylor evidence")
        if proved and self.missing_theorem is not None:
            raise ValueError("a proved Taylor level cannot retain a missing theorem")
        if not proved and self.classification == "continuous_no_initial_trapped_sphere":
            raise ValueError("unproved Taylor level cannot use the pass classification")
        if self.ode is not None and self.ode.chart_kind != PRODUCT_BOX_CHART_KIND:
            raise ValueError("Taylor ladder must keep chart_kind lambda_k")
        if self.ode is not None and self.ode.chart_coordinates:
            raise ValueError("Taylor ladder is not the (C,k) certificate chart")
        if self.ode is None:
            if self.complete_support_coverage:
                raise ValueError("a resource or domain stop cannot claim complete support tiling")
            if self.ode_cell_count != 0:
                raise ValueError("missing ODE record cannot report produced cells")
        else:
            left, right = sgbl_validate_taylor_cell_inventory(
                self.ode.cells, origin=self.ode.spec.support_minimum
            )
            if self.coverage_left != left or self.coverage_right != right:
                raise ValueError("coverage span does not match the gap-free inventory")
            if self.ode_cell_count != len(self.ode.cells):
                raise ValueError("ode_cell_count does not match the returned inventory")
            tiles = self.ode.tiles_compact_support
            if self.complete_support_coverage != tiles:
                raise ValueError("complete_support_coverage does not match compact-support tiling")
            if self.continuous_no_initial_trapped_sphere and not tiles:
                raise ValueError("a local pass requires complete support tiling")

    @property
    def sampled_nodes_are_not_the_certificate(self) -> bool:
        return True

    @property
    def chart_kind(self) -> str:
        return PRODUCT_BOX_CHART_KIND

    @property
    def product_box_is_not_the_certificate_chart(self) -> bool:
        return True

    def as_contract_mapping(self) -> dict[str, Any]:
        return {
            "name": self.level.name,
            "taylor_order": self.level.taylor_order,
            "max_bisection_depth": self.level.max_bisection_depth,
            "classification": self.classification,
            "obstruction": self.obstruction,
            "complete_support_coverage": self.complete_support_coverage,
            "coverage_left": _fraction_text(self.coverage_left),
            "coverage_right": _fraction_text(self.coverage_right),
            "compactness_upper": _fraction_text(self.compactness_upper),
            "exterior_compactness_upper": _fraction_text(self.exterior_compactness_upper),
            "compactness_margin": _fraction_text(self.compactness_margin),
            "ode_cell_count": self.ode_cell_count,
            "jet_evaluations": self.jet_evaluations,
            "compactness_jet_evaluations": self.compactness_jet_evaluations,
            "max_rational_bit_length_observed": self.max_rational_bit_length_observed,
            "taylor_strict_inclusions": self.taylor_strict_inclusions,
            "continuous_no_initial_trapped_sphere": self.continuous_no_initial_trapped_sphere,
            "resource_reason": self.resource_reason,
            "missing_theorem": self.missing_theorem,
            "sampled_nodes_are_not_the_certificate": True,
            "chart_kind": PRODUCT_BOX_CHART_KIND,
        }


def sgbl_evaluate_taylor_ladder_level(
    spec: SGBLExactInitialSlice,
    level: SGBLTaylorResourceLevel,
) -> SGBLTaylorLevelRecord:
    """Evaluate one declared or attack resource level.  Caps are not raised."""

    if type(spec) is not SGBLExactInitialSlice:
        raise TypeError("spec must be SGBLExactInitialSlice")
    if type(level) is not SGBLTaylorResourceLevel:
        raise TypeError("level must be SGBLTaylorResourceLevel")
    started = perf_counter()
    ode: SGBLTaylorODERecord | None = None
    obstruction: str | None = None
    resource_reason: str | None = None
    classification = "interval_inconclusive"
    jet_evaluations = 0
    compactness_jet_evaluations = 0
    observed_bits = 0
    try:
        ode = sgbl_validated_taylor_lambda_k_constraint_ode(spec, level)
    except SGBLTaylorEvaluationStop as exc:
        if exc.reason == "resource_limit":
            classification = "resource_limit"
            resource_reason = _resource_reason_from_payload(exc.payload)
            jet_evaluations = int(exc.payload.get("evaluations") or 0)
            observed_bits = int(exc.payload.get("observed") or 0)
        elif exc.reason == "domain_error":
            classification = "domain_error"
        else:
            obstruction = str(
                exc.payload.get("obstruction") or "zero_in_interval_reciprocal"
            )
    wall = perf_counter() - started
    if ode is not None:
        obstruction = ode.obstruction
        jet_evaluations = ode.jet_evaluations
        compactness_jet_evaluations = ode.compactness_jet_evaluations
        observed_bits = ode.max_rational_bit_length_observed
        if ode.chart_kind != PRODUCT_BOX_CHART_KIND or ode.chart_coordinates:
            raise SGBLTaylorTrapStop(
                "tampered_contract",
                "Taylor ladder received a non-lambda_k ODE record",
            )
    cells = ode.cells if ode is not None else ()
    if ode is not None:
        coverage_left, coverage_right = sgbl_taylor_inventory_span(cells, spec)
        complete = sgbl_taylor_cells_cover_support(cells, spec)
        if complete != ode.tiles_compact_support:
            raise SGBLTaylorTrapStop(
                "tampered_contract",
                "Taylor tiling bit does not match the validated inventory",
            )
    else:
        coverage_left = spec.support_minimum
        coverage_right = spec.support_minimum
        complete = False
    compactness_upper = _support_compactness_upper(ode)
    exterior_upper = compactness_upper
    if (
        complete
        and ode is not None
        and obstruction is None
        and resource_reason is None
        and classification not in {"resource_limit", "domain_error"}
    ):
        exterior_upper, exterior_reason = _exterior_upper(spec, ode)
        if exterior_reason is not None:
            obstruction = exterior_reason
        elif exterior_upper >= DECLARED_COMPACTNESS_STRICT_UPPER:
            obstruction = "exterior_compactness_not_below_one"
    margin = 1 - max(compactness_upper, exterior_upper)
    if classification not in {"resource_limit", "domain_error"}:
        classification = "interval_inconclusive"
    proved = _level_proved(
        classification="continuous_no_initial_trapped_sphere",
        obstruction=obstruction,
        complete_support_coverage=complete,
        compactness_upper=compactness_upper,
        exterior_compactness_upper=exterior_upper,
        compactness_margin=margin,
        ode=ode,
        resource_reason=resource_reason,
    )
    if proved:
        classification = "continuous_no_initial_trapped_sphere"
    if classification == "resource_limit":
        missing = "declared_taylor_model_resource_budget"
    elif classification == "domain_error":
        missing = "declared_taylor_model_physical_domain"
    elif obstruction is not None:
        missing = MISSING_THEOREMS.get(obstruction, obstruction)
    elif not complete:
        missing = MISSING_THEOREMS["incomplete_support_coverage"]
    else:
        missing = None
    return SGBLTaylorLevelRecord(
        level=level,
        classification=classification,
        obstruction=obstruction,
        complete_support_coverage=complete,
        compactness_upper=compactness_upper,
        exterior_compactness_upper=exterior_upper,
        compactness_margin=margin,
        coverage_left=coverage_left,
        coverage_right=coverage_right,
        ode_cell_count=len(cells),
        jet_evaluations=jet_evaluations,
        compactness_jet_evaluations=compactness_jet_evaluations,
        max_rational_bit_length_observed=observed_bits,
        taylor_strict_inclusions=sum(1 for cell in cells if cell.taylor_strict_self_map),
        continuous_no_initial_trapped_sphere=proved,
        resource_reason=resource_reason,
        missing_theorem=missing,
        diagnostic_wall_seconds=wall,
        ode=ode,
    )


def _shared_coverage_monotone(levels: tuple[SGBLTaylorLevelRecord, ...]) -> bool:
    for earlier, later in zip(levels, levels[1:]):
        earlier_cells = earlier.ode.cells if earlier.ode is not None else ()
        later_cells = later.ode.cells if later.ode is not None else ()
        if not sgbl_taylor_shared_coverage_monotone(earlier_cells, later_cells):
            return False
    return True


def _compactness_upper_nonincreasing(levels: tuple[SGBLTaylorLevelRecord, ...]) -> bool:
    previous: Fraction | None = None
    for record in levels:
        if record.ode is None or not record.ode.cells:
            continue
        if previous is not None and record.compactness_upper > previous:
            return False
        previous = record.compactness_upper
    return True


def _contract_payload(
    spec: SGBLExactInitialSlice,
    ladder: tuple[SGBLTaylorResourceLevel, ...],
    levels: tuple[SGBLTaylorLevelRecord, ...],
    *,
    shared_coverage_monotone: bool,
    compactness_upper_nonincreasing: bool,
    predecessor_trap_refinement_contract_sha256: str,
    predecessor_trap_refinement2_contract_sha256: str,
) -> dict[str, Any]:
    proved_names = tuple(
        record.level.name
        for record in levels
        if record.continuous_no_initial_trapped_sphere
    )
    return {
        "INSTRUMENT_ID": INSTRUMENT_ID,
        "chart_kind": PRODUCT_BOX_CHART_KIND,
        "certificate_chart_kind": CERTIFICATE_CHART_KIND,
        "product_box_is_not_the_certificate_chart": True,
        "taylor_order": DECLARED_TAYLOR_ORDER,
        "chi_amplitude": _fraction_text(spec.chi_amplitude),
        "phi_amplitude": _fraction_text(spec.phi_amplitude),
        "center": _fraction_text(spec.center),
        "half_width": _fraction_text(spec.half_width),
        "compactness_gate": "strict_C_lt_1",
        "compactness_strict_upper": _fraction_text(Q(DECLARED_COMPACTNESS_STRICT_UPPER)),
        "base_cells": DECLARED_BASE_CELLS,
        "declared_depths": list(DECLARED_TAYLOR_LEVEL_DEPTHS),
        "ladder": [level.as_contract_mapping() for level in ladder],
        "levels": [record.as_contract_mapping() for record in levels],
        "evaluated_every_declared_level": len(levels) == len(DECLARED_TAYLOR_LEVEL_DEPTHS),
        "shared_coverage_monotone": shared_coverage_monotone,
        "compactness_upper_nonincreasing": compactness_upper_nonincreasing,
        "proved_level_names": list(proved_names),
        "any_declared_level_proves_continuous_no_initial_trap": bool(proved_names),
        "predecessor_trap_refinement_contract_sha256": predecessor_trap_refinement_contract_sha256,
        "predecessor_trap_refinement2_contract_sha256": predecessor_trap_refinement2_contract_sha256,
        "predecessors_are_regression_only": True,
        "sampled_finite_differences_are_not_the_proof": True,
        "SGBL_branch_owned_and_healthy": False,
        "FRZ1": False,
        "PREF1": False,
        "holdout_authorized": False,
        "execution_authorized": False,
        "sampled_nodes_are_not_the_certificate": True,
    }


def sgbl_taylor_contract_sha256(payload: Mapping[str, Any]) -> str:
    raw = dumps(dict(payload), sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return sha256(raw.encode("ascii")).hexdigest()


@dataclass(frozen=True, slots=True, kw_only=True)
class SGBLTaylorLadderRecord:
    """Complete evaluation of the frozen Taylor-model resource ladder.

    A true local no-trap bit is not ``SGBL_branch_owned_and_healthy``,
    ``FRZ1``, ``PREF1``, or holdout.
    """

    slice: SGBLExactInitialSlice
    ladder: tuple[SGBLTaylorResourceLevel, ...]
    levels: tuple[SGBLTaylorLevelRecord, ...]
    predecessor_trap_refinement_contract_sha256: str
    predecessor_trap_refinement2_contract_sha256: str
    shared_coverage_monotone: bool
    compactness_upper_nonincreasing: bool
    proved_level_names: tuple[str, ...]
    any_declared_level_proves_continuous_no_initial_trap: bool
    contract_payload: Mapping[str, Any]
    contract_sha256: str

    def __post_init__(self) -> None:
        if type(self.slice) is not SGBLExactInitialSlice:
            raise TypeError("slice must be SGBLExactInitialSlice")
        _require_predecessor_hashes(
            self.predecessor_trap_refinement_contract_sha256,
            self.predecessor_trap_refinement2_contract_sha256,
        )
        if self.ladder != DECLARED_TAYLOR_RESOURCE_LADDER:
            raise SGBLTaylorTrapStop(
                "changed_policy",
                "aggregate record must keep the frozen Taylor resource ladder",
            )
        if len(self.levels) != len(DECLARED_TAYLOR_LEVEL_DEPTHS):
            raise SGBLTaylorTrapStop(
                "skipped_level",
                "aggregate record must evaluate every declared Taylor level",
            )
        if tuple(record.level for record in self.levels) != DECLARED_TAYLOR_RESOURCE_LADDER:
            raise SGBLTaylorTrapStop(
                "reordered_levels",
                "aggregate record levels must follow the frozen Taylor ladder order",
            )
        expected_monotone = _shared_coverage_monotone(self.levels)
        expected_nonincreasing = _compactness_upper_nonincreasing(self.levels)
        if self.shared_coverage_monotone != expected_monotone:
            raise SGBLTaylorTrapStop(
                "tampered_contract",
                "shared-coverage monotone bit does not match the returned cells",
            )
        if self.compactness_upper_nonincreasing != expected_nonincreasing:
            raise SGBLTaylorTrapStop(
                "tampered_contract",
                "compactness-upper monotone bit does not match the returned levels",
            )
        proved = tuple(
            record.level.name
            for record in self.levels
            if record.continuous_no_initial_trapped_sphere
        )
        if self.proved_level_names != proved:
            raise SGBLTaylorTrapStop(
                "tampered_contract",
                "proved level names do not match the level records",
            )
        if self.any_declared_level_proves_continuous_no_initial_trap != bool(proved):
            raise SGBLTaylorTrapStop(
                "tampered_contract",
                "aggregate pass bit does not match proved levels",
            )
        if proved and not self.shared_coverage_monotone:
            raise SGBLTaylorTrapStop(
                "tampered_contract",
                "a Taylor pass cannot be claimed on nonmonotone shared coverage",
            )
        expected_payload = _contract_payload(
            self.slice,
            self.ladder,
            self.levels,
            shared_coverage_monotone=self.shared_coverage_monotone,
            compactness_upper_nonincreasing=self.compactness_upper_nonincreasing,
            predecessor_trap_refinement_contract_sha256=(
                self.predecessor_trap_refinement_contract_sha256
            ),
            predecessor_trap_refinement2_contract_sha256=(
                self.predecessor_trap_refinement2_contract_sha256
            ),
        )
        if dict(self.contract_payload) != expected_payload:
            raise SGBLTaylorTrapStop(
                "tampered_contract",
                "contract payload does not match the frozen Taylor ladder record",
            )
        expected_hash = sgbl_taylor_contract_sha256(expected_payload)
        if self.contract_sha256 != expected_hash:
            raise SGBLTaylorTrapStop(
                "tampered_contract",
                "contract hash does not match the nonpromoting payload",
            )
        if any(
            dict(self.contract_payload)[name] is not False
            for name in (
                "SGBL_branch_owned_and_healthy",
                "FRZ1",
                "PREF1",
                "holdout_authorized",
                "execution_authorized",
            )
        ):
            raise SGBLTaylorTrapStop(
                "tampered_contract",
                "nonpromoting contract payload cannot set aggregate health true",
            )

    @property
    def sampled_nodes_are_not_the_certificate(self) -> bool:
        return True

    @property
    def product_box_is_not_the_certificate_chart(self) -> bool:
        return True

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
    def holdout_authorized(self) -> bool:
        return False

    @property
    def execution_authorized(self) -> bool:
        return False

    @property
    def evaluated_every_declared_level(self) -> bool:
        return len(self.levels) == len(DECLARED_TAYLOR_LEVEL_DEPTHS)

    @property
    def diagnostic_wall_seconds(self) -> tuple[float, ...]:
        return tuple(record.diagnostic_wall_seconds for record in self.levels)


def sgbl_evaluate_taylor_resource_ladder(
    spec: SGBLExactInitialSlice | None = None,
    *,
    ladder: tuple[SGBLTaylorResourceLevel, ...] | None = None,
    predecessor_trap_refinement_contract_sha256: str | None = None,
    predecessor_trap_refinement2_contract_sha256: str | None = None,
) -> SGBLTaylorLadderRecord:
    """Evaluate every frozen Taylor level.  A pass does not stop later levels."""

    if spec is None:
        spec = SGBLExactInitialSlice()
    elif type(spec) is not SGBLExactInitialSlice:
        raise TypeError("spec must be SGBLExactInitialSlice")
    if spec.chi_amplitude != NOMINAL_CHI_AMPLITUDE and spec.chi_amplitude != SMALL_AMPLITUDE_CONTROL:
        raise ValueError(
            "Taylor ladder slice must be nominal A_chi=3 or the named small-amplitude control"
        )
    frozen = _require_declared_ladder(
        DECLARED_TAYLOR_RESOURCE_LADDER if ladder is None else ladder
    )
    predecessor_one, predecessor_two = _require_predecessor_hashes(
        PREDECESSOR_TRAP_REFINEMENT_CONTRACT_SHA256
        if predecessor_trap_refinement_contract_sha256 is None
        else predecessor_trap_refinement_contract_sha256,
        PREDECESSOR_TRAP_REFINEMENT2_CONTRACT_SHA256
        if predecessor_trap_refinement2_contract_sha256 is None
        else predecessor_trap_refinement2_contract_sha256,
    )
    records = []
    for level in frozen:
        records.append(sgbl_evaluate_taylor_ladder_level(spec, level))
    levels = tuple(records)
    monotone = _shared_coverage_monotone(levels)
    nonincreasing = _compactness_upper_nonincreasing(levels)
    proved = tuple(
        record.level.name for record in levels if record.continuous_no_initial_trapped_sphere
    )
    payload = _contract_payload(
        spec,
        frozen,
        levels,
        shared_coverage_monotone=monotone,
        compactness_upper_nonincreasing=nonincreasing,
        predecessor_trap_refinement_contract_sha256=predecessor_one,
        predecessor_trap_refinement2_contract_sha256=predecessor_two,
    )
    return SGBLTaylorLadderRecord(
        slice=spec,
        ladder=frozen,
        levels=levels,
        predecessor_trap_refinement_contract_sha256=predecessor_one,
        predecessor_trap_refinement2_contract_sha256=predecessor_two,
        shared_coverage_monotone=monotone,
        compactness_upper_nonincreasing=nonincreasing,
        proved_level_names=proved,
        any_declared_level_proves_continuous_no_initial_trap=bool(proved),
        contract_payload=MappingProxyType(payload),
        contract_sha256=sgbl_taylor_contract_sha256(payload),
    )


def sgbl_taylor_health_gate(record: SGBLTaylorLadderRecord) -> dict[str, Any]:
    """Aggregate flags stay false even if a declared Taylor level locally passes."""

    if type(record) is not SGBLTaylorLadderRecord:
        raise TypeError("record must be SGBLTaylorLadderRecord")
    return {
        "any_declared_level_proves_continuous_no_initial_trap": (
            record.any_declared_level_proves_continuous_no_initial_trap
        ),
        "proved_level_names": record.proved_level_names,
        "evaluated_every_declared_level": record.evaluated_every_declared_level,
        "shared_coverage_monotone": record.shared_coverage_monotone,
        "sampled_nodes_are_not_the_certificate": True,
        "product_box_is_not_the_certificate_chart": True,
        "predecessors_are_regression_only": True,
        "sampled_finite_differences_are_not_the_proof": True,
        "SGBL_branch_owned_and_healthy": False,
        "FRZ1": False,
        "PREF1": False,
        "holdout_authorized": False,
        "execution_authorized": False,
        "copied_gr0_or_fgcqr_health_evidence": False,
    }


__all__ = [
    "BUMP_SCALED_DERIVATIVE_BOUNDS",
    "DECLARED_BUMP_JET_ORDER",
    "DECLARED_COMPACTNESS_STRICT_UPPER",
    "DECLARED_JET_CALLS_PER_PARTITION",
    "DECLARED_TAYLOR_LEVEL_DEPTHS",
    "DECLARED_TAYLOR_ORDER",
    "DECLARED_TAYLOR_RESOURCE_LADDER",
    "INSTRUMENT_ID",
    "IntervalJet",
    "NOMINAL_CHI_AMPLITUDE",
    "PREDECESSOR_TRAP_REFINEMENT2_CONTRACT_SHA256",
    "PREDECESSOR_TRAP_REFINEMENT_CONTRACT_SHA256",
    "PRODUCT_BOX_CHART_KIND",
    "SGBLTaylorCellEnclosure",
    "SGBLTaylorEvaluationStop",
    "SGBLTaylorLadderRecord",
    "SGBLTaylorLevelRecord",
    "SGBLTaylorODERecord",
    "SGBLTaylorResourceLevel",
    "SGBLTaylorTrapStop",
    "SMALL_AMPLITUDE_CONTROL",
    "TaylorModel",
    "sgbl_affine_constraint_jets",
    "sgbl_algebraic_compactness_interval",
    "sgbl_compact_bump_interval_jet",
    "sgbl_declared_taylor_cell_cap",
    "sgbl_declared_taylor_jet_cap",
    "sgbl_evaluate_taylor_ladder_level",
    "sgbl_evaluate_taylor_resource_ladder",
    "sgbl_interval_jet_shift",
    "sgbl_lagrange_remainder",
    "sgbl_require_complete_jet",
    "sgbl_require_honest_remainder",
    "sgbl_taylor_cells_cover_support",
    "sgbl_taylor_contract_sha256",
    "sgbl_taylor_health_gate",
    "sgbl_taylor_inventory_span",
    "sgbl_taylor_known_linear_ode_control",
    "sgbl_taylor_model_from_jet",
    "sgbl_taylor_polynomial_ode_control",
    "sgbl_taylor_shared_coverage_monotone",
    "sgbl_validate_taylor_cell_inventory",
    "sgbl_validated_taylor_lambda_k_constraint_ode",
]
