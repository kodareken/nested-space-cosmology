"""Independent exact evaluator for the prospective TDG9 arithmetic gate.

This standard-library module deliberately does not import the primary TDG9
implementation or historical TDG5/TDG6 private arithmetic.  It independently
rebuilds binary64 Hermite coefficients, half restrictions, and difference
cubics.  Cubic extrema are certified by exact rational derivative-root
isolation, repeated bisection, and interval Horner evaluation.

The returned rational intervals enclose generally algebraic absolute
suprema.  They certify only the finite-dimensional Hermite replay supplied to
this function; they do not certify a historical proposal, a PDE trajectory,
or an asymptotic temporal order.  Work is bounded by explicit depth and node
limits, and ambiguity is a typed resource exhaustion.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from hashlib import sha256
from math import isfinite, isqrt
from numbers import Rational
from typing import Final, Iterable, Sequence


Q = Fraction
TDG9_INDEPENDENT_SCHEMA_VERSION: Final[int] = 1
TDG9_INDEPENDENT_EVALUATOR_ID: Final[str] = (
    "tdg9_exact_derivative_isolation_interval_horner_v1"
)
_HASH_DOMAIN: Final[bytes] = b"TDG9-EXACT-CUBIC-STREAM-v1\n"
_COMBINED_HASH_DOMAIN: Final[str] = "TDG9-EXACT-COMBINED-v1\n"
_RESOURCE_EXHAUSTION_REASONS: Final[frozenset[str]] = frozenset(
    {
        "max_nodes_exhausted_before_stream_completed",
        "max_nodes_exhausted_before_contraction_resolved",
        "max_depth_reached_before_contraction_resolved",
    }
)

Cubic = tuple[Fraction, Fraction, Fraction, Fraction]
Polynomial = tuple[Fraction, ...]


def _sha256_digest(value: object, *, name: str) -> str:
    if (
        not isinstance(value, str)
        or len(value) != 64
        or any(character not in "0123456789abcdef" for character in value)
    ):
        raise ValueError(f"{name} must be a lowercase SHA-256 digest")
    return value


def _combined_coefficient_hash(outer_digest: str, finest_digest: str) -> str:
    return sha256(
        (_COMBINED_HASH_DOMAIN + outer_digest + "\n" + finest_digest + "\n").encode(
            "ascii"
        )
    ).hexdigest()


def _nonnegative(value: Rational, *, name: str) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, Rational):
        raise TypeError(f"{name} must be rational")
    answer = Q(value)
    if answer < 0:
        raise ValueError(f"{name} must be nonnegative")
    return answer


def _integer(value: object, *, name: str, zero_allowed: bool = False) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{name} must be an integer")
    if value < (0 if zero_allowed else 1):
        raise ValueError(
            f"{name} must be {'nonnegative' if zero_allowed else 'positive'}"
        )
    return value


def _sequence(value: object, *, size: int, name: str) -> tuple[object, ...]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise TypeError(f"{name} must be a sequence")
    answer = tuple(value)
    if len(answer) != size:
        raise ValueError(f"{name} must contain exactly {size} entries")
    return answer


def _dyadic(value: object, *, name: str) -> Fraction:
    if type(value) is not float:
        raise TypeError(f"{name} must be a built-in binary64 float")
    if not isfinite(value):
        raise ValueError(f"{name} must be finite")
    numerator, denominator = value.as_integer_ratio()
    return Q(numerator, denominator)


def _segment(value: object, *, name: str) -> tuple[Fraction, ...]:
    raw = _sequence(value, size=5, name=name)
    answer = tuple(
        _dyadic(item, name=f"{name}[{index}]") for index, item in enumerate(raw)
    )
    if answer[4] <= 0:
        raise ValueError(f"{name} width must be positive")
    return answer


def independent_exact_hermite_coefficients(segment: object) -> Cubic:
    """Independently expand the four normalized cubic-Hermite basis terms."""

    y0, f0, y1, f1, width = _segment(segment, name="segment")
    s0 = width * f0
    s1 = width * f1
    # Basis vectors are expanded directly rather than copied from the primary.
    basis = (
        (Q(1), Q(0), Q(-3), Q(2)),
        (Q(0), Q(1), Q(-2), Q(1)),
        (Q(0), Q(0), Q(3), Q(-2)),
        (Q(0), Q(0), Q(-1), Q(1)),
    )
    data = (y0, s0, y1, s1)
    return tuple(
        sum(
            (data[index] * basis[index][power] for index in range(4)),
            Q(0),
        )
        for power in range(4)
    )  # type: ignore[return-value]


def _exact_cubic(value: object, *, name: str) -> Cubic:
    raw = _sequence(value, size=4, name=name)
    answer: list[Fraction] = []
    for index, item in enumerate(raw):
        if isinstance(item, bool) or not isinstance(item, Rational):
            raise TypeError(f"{name}[{index}] must be rational")
        answer.append(Q(item))
    return tuple(answer)  # type: ignore[return-value]


def independent_restrict_exact_cubic_to_half(
    coefficients: object, *, half: int
) -> Cubic:
    """Independently substitute ``t=(half+s)/2`` into one exact cubic."""

    a0, a1, a2, a3 = _exact_cubic(coefficients, name="coefficients")
    if isinstance(half, bool) or half not in (0, 1):
        raise ValueError("half must be zero or one")
    offset = Q(half, 2)
    scale = Q(1, 2)
    # Binomial expansion of sum a_k (offset + scale*s)^k.
    return (
        a0 + a1 * offset + a2 * offset**2 + a3 * offset**3,
        scale * (a1 + 2 * a2 * offset + 3 * a3 * offset**2),
        scale**2 * (a2 + 3 * a3 * offset),
        scale**3 * a3,
    )


def independent_subtract_exact_cubics(left: object, right: object) -> Cubic:
    left_values = _exact_cubic(left, name="left")
    right_values = _exact_cubic(right, name="right")
    return tuple(left_values[index] - right_values[index] for index in range(4))  # type: ignore[return-value]


def independent_three_halves_order_passes_exact(
    *, outer_lower_squared: Rational, finest_upper_squared: Rational
) -> bool:
    """Independently apply the sufficient ``8 U12^2 <= L01^2`` test."""

    outer = _nonnegative(outer_lower_squared, name="outer_lower_squared")
    finest = _nonnegative(finest_upper_squared, name="finest_upper_squared")
    return 8 * finest <= outer


@dataclass(frozen=True, slots=True)
class IndependentExactRationalInterval:
    lower: Fraction
    upper: Fraction

    def __post_init__(self) -> None:
        lower = _nonnegative(self.lower, name="lower")
        upper = _nonnegative(self.upper, name="upper")
        if lower > upper:
            raise ValueError("interval lower bound exceeds upper bound")
        object.__setattr__(self, "lower", lower)
        object.__setattr__(self, "upper", upper)


@dataclass(frozen=True, slots=True)
class IndependentExactSupremumEvidence:
    interval: IndependentExactRationalInterval
    polynomial_count: int
    coefficient_stream_sha256: str
    nodes_visited: int
    maximum_depth_reached: int
    evaluator_id: str = TDG9_INDEPENDENT_EVALUATOR_ID

    def __post_init__(self) -> None:
        if not isinstance(self.interval, IndependentExactRationalInterval):
            raise TypeError("interval must be IndependentExactRationalInterval")
        _integer(self.polynomial_count, name="polynomial_count")
        _integer(self.nodes_visited, name="nodes_visited")
        _integer(
            self.maximum_depth_reached,
            name="maximum_depth_reached",
            zero_allowed=True,
        )
        _sha256_digest(
            self.coefficient_stream_sha256,
            name="coefficient stream hash",
        )
        if self.evaluator_id != TDG9_INDEPENDENT_EVALUATOR_ID:
            raise ValueError("unexpected independent evaluator identifier")


@dataclass(frozen=True, slots=True)
class IndependentExactTemporalArithmeticEvidence:
    outer_difference: IndependentExactSupremumEvidence
    finest_difference: IndependentExactSupremumEvidence
    row_count: int
    classification: str
    sufficient_pass_certified: bool
    contraction_threshold_resolved: bool
    contraction_threshold_passed: bool | None
    sufficient_pass_left: Fraction
    sufficient_pass_right: Fraction
    max_depth: int
    max_nodes: int
    nodes_visited: int
    combined_coefficient_stream_sha256: str
    schema_version: int = TDG9_INDEPENDENT_SCHEMA_VERSION
    evaluator_id: str = TDG9_INDEPENDENT_EVALUATOR_ID
    absolute_tolerance_used: bool = False
    physical_signal_used_for_normalization: bool = False
    replay_identity_certified: bool = False
    historical_proposal_certified: bool = False
    pde_trajectory_order_certified: bool = False

    def __post_init__(self) -> None:
        row_count = _integer(self.row_count, name="row_count")
        max_depth = _integer(self.max_depth, name="max_depth", zero_allowed=True)
        max_nodes = _integer(self.max_nodes, name="max_nodes")
        nodes_visited = _integer(self.nodes_visited, name="nodes_visited")
        if not isinstance(
            self.outer_difference, IndependentExactSupremumEvidence
        ) or not isinstance(self.finest_difference, IndependentExactSupremumEvidence):
            raise TypeError(
                "difference evidence must use IndependentExactSupremumEvidence"
            )
        if self.outer_difference.polynomial_count != 2 * row_count:
            raise ValueError("outer polynomial count must equal twice row_count")
        if self.finest_difference.polynomial_count != 4 * row_count:
            raise ValueError("finest polynomial count must equal four times row_count")
        nested_nodes = (
            self.outer_difference.nodes_visited + self.finest_difference.nodes_visited
        )
        if nodes_visited != nested_nodes:
            raise ValueError("total nodes differ from nested supremum evidence")
        if nodes_visited > max_nodes:
            raise ValueError("visited nodes exceed the independent node budget")
        if (
            self.outer_difference.maximum_depth_reached > max_depth
            or self.finest_difference.maximum_depth_reached > max_depth
        ):
            raise ValueError("nested supremum depth exceeds the declared maximum")
        if self.classification not in {
            "exact_zero",
            "sufficient_contraction_pass",
            "sufficient_contraction_failure",
        }:
            raise ValueError("unknown independent temporal classification")
        if self.sufficient_pass_certified is not (
            self.classification == "sufficient_contraction_pass"
        ):
            raise ValueError("sufficient-pass flag differs from classification")
        resolved = self.classification != "exact_zero"
        if self.contraction_threshold_resolved is not resolved:
            raise ValueError("threshold-resolution flag differs from classification")
        expected = (
            True
            if self.classification == "sufficient_contraction_pass"
            else False
            if self.classification == "sufficient_contraction_failure"
            else None
        )
        if self.contraction_threshold_passed is not expected:
            raise ValueError("threshold-pass flag differs from classification")
        if self.sufficient_pass_left != 8 * self.finest_difference.interval.upper**2:
            raise ValueError("sufficient-pass left side differs from 8 U12^2")
        if self.sufficient_pass_right != self.outer_difference.interval.lower**2:
            raise ValueError("sufficient-pass right side differs from L01^2")
        if (
            self.schema_version != TDG9_INDEPENDENT_SCHEMA_VERSION
            or self.evaluator_id != TDG9_INDEPENDENT_EVALUATOR_ID
            or self.absolute_tolerance_used
            or self.physical_signal_used_for_normalization
            or self.replay_identity_certified
            or self.historical_proposal_certified
            or self.pde_trajectory_order_certified
        ):
            raise ValueError("independent evidence crossed its scope")
        combined_hash = _sha256_digest(
            self.combined_coefficient_stream_sha256,
            name="combined coefficient hash",
        )
        expected_combined_hash = _combined_coefficient_hash(
            self.outer_difference.coefficient_stream_sha256,
            self.finest_difference.coefficient_stream_sha256,
        )
        if combined_hash != expected_combined_hash:
            raise ValueError("combined coefficient hash differs from component hashes")


@dataclass(frozen=True, slots=True)
class IndependentExactResourceExhaustionEvidence:
    reason: str
    max_depth: int
    max_nodes: int
    nodes_visited: int
    completed_rows: int
    outer_interval: IndependentExactRationalInterval | None
    finest_interval: IndependentExactRationalInterval | None
    outer_coefficient_prefix_sha256: str
    finest_coefficient_prefix_sha256: str
    evaluator_id: str = TDG9_INDEPENDENT_EVALUATOR_ID

    def __post_init__(self) -> None:
        _integer(self.max_depth, name="max_depth", zero_allowed=True)
        max_nodes = _integer(self.max_nodes, name="max_nodes")
        nodes_visited = _integer(
            self.nodes_visited, name="nodes_visited", zero_allowed=True
        )
        completed_rows = _integer(
            self.completed_rows, name="completed_rows", zero_allowed=True
        )
        if self.reason not in _RESOURCE_EXHAUSTION_REASONS:
            raise ValueError("unknown independent exhaustion reason")
        if nodes_visited > max_nodes:
            raise ValueError("visited nodes exceed the independent node budget")
        for interval, name in (
            (self.outer_interval, "outer_interval"),
            (self.finest_interval, "finest_interval"),
        ):
            if interval is not None and not isinstance(
                interval, IndependentExactRationalInterval
            ):
                raise TypeError(
                    f"{name} must be IndependentExactRationalInterval or None"
                )
        if self.outer_interval is None and self.finest_interval is not None:
            raise ValueError("finest interval cannot precede an outer interval")
        if completed_rows and (
            self.outer_interval is None or self.finest_interval is None
        ):
            raise ValueError("completed rows require both partial intervals")
        if self.reason == "max_nodes_exhausted_before_stream_completed":
            if nodes_visited != max_nodes:
                raise ValueError("stream node exhaustion must consume the node budget")
        else:
            if completed_rows == 0:
                raise ValueError("contraction exhaustion requires a completed row")
            if self.outer_interval is None or self.finest_interval is None:
                raise ValueError("contraction exhaustion requires both intervals")
            if (
                self.reason == "max_nodes_exhausted_before_contraction_resolved"
                and nodes_visited < max_nodes - 1
            ):
                raise ValueError(
                    "contraction node exhaustion is inconsistent with budget"
                )
        _sha256_digest(
            self.outer_coefficient_prefix_sha256,
            name="outer coefficient prefix hash",
        )
        _sha256_digest(
            self.finest_coefficient_prefix_sha256,
            name="finest coefficient prefix hash",
        )
        if self.evaluator_id != TDG9_INDEPENDENT_EVALUATOR_ID:
            raise ValueError("unexpected independent evaluator identifier")


class IndependentExactTemporalArithmeticResourceExhausted(RuntimeError):
    def __init__(self, evidence: IndependentExactResourceExhaustionEvidence) -> None:
        self.evidence = evidence
        super().__init__(f"tdg9_independent_resource_exhausted:{evidence.reason}")


class _NodeBudgetExhausted(Exception):
    """Private control-flow sentinel owned only by this evaluator."""


@dataclass(slots=True)
class _NodeBudget:
    maximum: int
    used: int = 0

    def claim(self, count: int) -> bool:
        if self.used + count > self.maximum:
            return False
        self.used += count
        return True


def _evaluate(coefficients: Cubic, point: Fraction) -> Fraction:
    a0, a1, a2, a3 = coefficients
    return a0 + point * (a1 + point * (a2 + point * a3))


def _trim(polynomial: Polynomial) -> Polynomial:
    values = list(polynomial)
    while len(values) > 1 and values[-1] == 0:
        values.pop()
    return tuple(values)


def _poly_evaluate(polynomial: Polynomial, point: Fraction) -> Fraction:
    answer = Q(0)
    for coefficient in reversed(polynomial):
        answer = coefficient + point * answer
    return answer


def _poly_derivative(polynomial: Polynomial) -> Polynomial:
    if len(polynomial) <= 1:
        return (Q(0),)
    return _trim(
        tuple(index * polynomial[index] for index in range(1, len(polynomial)))
    )


def _poly_remainder(dividend: Polynomial, divisor: Polynomial) -> Polynomial:
    numerator = list(_trim(dividend))
    denominator = _trim(divisor)
    if denominator == (Q(0),):
        raise ZeroDivisionError("zero polynomial divisor")
    while len(numerator) >= len(denominator) and any(numerator):
        scale = numerator[-1] / denominator[-1]
        shift = len(numerator) - len(denominator)
        for index, coefficient in enumerate(denominator):
            numerator[index + shift] -= scale * coefficient
        while len(numerator) > 1 and numerator[-1] == 0:
            numerator.pop()
    return _trim(tuple(numerator))


def _sturm_sequence(polynomial: Polynomial) -> tuple[Polynomial, ...]:
    first = _trim(polynomial)
    second = _poly_derivative(first)
    sequence = [first, second]
    while second != (Q(0),):
        remainder = _poly_remainder(first, second)
        if remainder == (Q(0),):
            break
        remainder = tuple(-item for item in remainder)
        sequence.append(remainder)
        first, second = second, remainder
    return tuple(sequence)


def _sign_variations(sequence: tuple[Polynomial, ...], point: Fraction) -> int:
    signs: list[int] = []
    for polynomial in sequence:
        value = _poly_evaluate(polynomial, point)
        if value:
            signs.append(1 if value > 0 else -1)
    return sum(left != right for left, right in zip(signs, signs[1:]))


def _rational_square_root(value: Fraction) -> Fraction | None:
    if value < 0:
        return None
    numerator = isqrt(value.numerator)
    denominator = isqrt(value.denominator)
    if numerator * numerator != value.numerator:
        return None
    if denominator * denominator != value.denominator:
        return None
    return Q(numerator, denominator)


def _interval_product(
    left: tuple[Fraction, Fraction], right: tuple[Fraction, Fraction]
) -> tuple[Fraction, Fraction]:
    products = tuple(a * b for a in left for b in right)
    return min(products), max(products)


def _interval_horner(
    coefficients: Cubic, lower: Fraction, upper: Fraction
) -> tuple[Fraction, Fraction]:
    variable = (lower, upper)
    interval = (coefficients[3], coefficients[3])
    for coefficient in (coefficients[2], coefficients[1], coefficients[0]):
        product = _interval_product(interval, variable)
        interval = (product[0] + coefficient, product[1] + coefficient)
    return interval


def _absolute_upper(interval: tuple[Fraction, Fraction]) -> Fraction:
    return max(abs(interval[0]), abs(interval[1]))


def _derivative_candidates(
    coefficients: Cubic,
    *,
    max_depth: int,
    budget: _NodeBudget,
) -> tuple[
    tuple[Fraction, ...],
    tuple[tuple[Fraction, Fraction], ...],
    int,
    bool,
]:
    c = coefficients[1]
    b = 2 * coefficients[2]
    a = 3 * coefficients[3]
    if a == 0:
        if b == 0:
            return (), (), 0, False
        root = -c / b
        return ((root,) if 0 < root < 1 else ()), (), 0, False
    discriminant = b * b - 4 * a * c
    if discriminant < 0:
        return (), (), 0, False
    square_root = _rational_square_root(discriminant)
    if square_root is not None:
        roots = tuple(
            sorted(
                {
                    root
                    for root in (
                        (-b - square_root) / (2 * a),
                        (-b + square_root) / (2 * a),
                    )
                    if 0 < root < 1
                }
            )
        )
        return roots, (), 0, False

    derivative = (c, b, a)
    sturm = _sturm_sequence(derivative)
    root_count = _sign_variations(sturm, Q(0)) - _sign_variations(sturm, Q(1))
    if root_count <= 0:
        return (), (), 0, False
    stack: list[tuple[Fraction, Fraction, int, int]] = [(Q(0), Q(1), 0, root_count)]
    leaves: list[tuple[Fraction, Fraction]] = []
    maximum_depth_reached = 0
    limited = False
    while stack:
        lower, upper, depth, count = stack.pop()
        maximum_depth_reached = max(maximum_depth_reached, depth)
        if depth >= max_depth:
            leaves.append((lower, upper))
            continue
        if not budget.claim(2):
            leaves.append((lower, upper))
            limited = True
            continue
        midpoint = (lower + upper) / 2
        left_count = _sign_variations(sturm, lower) - _sign_variations(sturm, midpoint)
        right_count = count - left_count
        if right_count:
            stack.append((midpoint, upper, depth + 1, right_count))
        if left_count:
            stack.append((lower, midpoint, depth + 1, left_count))
    return (), tuple(sorted(leaves)), maximum_depth_reached, limited


def _certify_one(
    coefficients: Cubic, *, max_depth: int, budget: _NodeBudget
) -> tuple[IndependentExactRationalInterval, int, int, bool]:
    start_nodes = budget.used
    if not budget.claim(1):
        raise _NodeBudgetExhausted
    endpoint_values = (
        abs(_evaluate(coefficients, Q(0))),
        abs(_evaluate(coefficients, Q(1))),
    )
    lower = max(endpoint_values)
    upper = lower
    roots, root_intervals, reached, limited = _derivative_candidates(
        coefficients, max_depth=max_depth, budget=budget
    )
    for root in roots:
        value = abs(_evaluate(coefficients, root))
        lower = max(lower, value)
        upper = max(upper, value)
    for left, right in root_intervals:
        midpoint = (left + right) / 2
        lower = max(lower, abs(_evaluate(coefficients, midpoint)))
        upper = max(upper, _absolute_upper(_interval_horner(coefficients, left, right)))
    return (
        IndependentExactRationalInterval(lower, upper),
        reached,
        budget.used - start_nodes,
        limited,
    )


def _hash_coefficients(hasher: object, ordinal: int, coefficients: Cubic) -> None:
    line = (
        f"{ordinal}|"
        + "|".join(f"{item.numerator}/{item.denominator}" for item in coefficients)
        + "\n"
    )
    hasher.update(line.encode("ascii"))  # type: ignore[attr-defined]


def _optional_interval(
    lower: Fraction | None, upper: Fraction | None
) -> IndependentExactRationalInterval | None:
    if lower is None or upper is None:
        return None
    return IndependentExactRationalInterval(lower, upper)


def assess_exact_temporal_refinement_independently(
    rows: Iterable[object], *, max_depth: int, max_nodes: int
) -> IndependentExactTemporalArithmeticEvidence:
    """Independently certify the streamed Hermite contraction comparison."""

    depth_limit = _integer(max_depth, name="max_depth", zero_allowed=True)
    node_limit = _integer(max_nodes, name="max_nodes")
    budget = _NodeBudget(node_limit)
    outer_hasher = sha256(_HASH_DOMAIN)
    fine_hasher = sha256(_HASH_DOMAIN)
    outer_lower: Fraction | None = None
    outer_upper: Fraction | None = None
    fine_lower: Fraction | None = None
    fine_upper: Fraction | None = None
    outer_count = fine_count = row_count = 0
    outer_nodes = fine_nodes = 0
    outer_depth = fine_depth = 0
    limited = False

    try:
        for row_index, row in enumerate(rows):
            raw_outer, raw_medium, raw_fine = _sequence(
                row, size=3, name=f"rows[{row_index}]"
            )
            medium = _sequence(raw_medium, size=2, name=f"rows[{row_index}].medium")
            fine = _sequence(raw_fine, size=4, name=f"rows[{row_index}].fine")
            outer_segment = _segment(raw_outer, name=f"rows[{row_index}].outer")
            medium_segments = tuple(
                _segment(item, name=f"rows[{row_index}].medium[{index}]")
                for index, item in enumerate(medium)
            )
            fine_segments = tuple(
                _segment(item, name=f"rows[{row_index}].fine[{index}]")
                for index, item in enumerate(fine)
            )
            width = outer_segment[4]
            if any(item[4] != width / 2 for item in medium_segments):
                raise ValueError("medium widths must equal half the outer width")
            if any(item[4] != width / 4 for item in fine_segments):
                raise ValueError(
                    "fine widths must equal one quarter of the outer width"
                )

            outer_coefficients = independent_exact_hermite_coefficients(raw_outer)
            medium_coefficients = tuple(
                independent_exact_hermite_coefficients(item) for item in medium
            )
            fine_coefficients = tuple(
                independent_exact_hermite_coefficients(item) for item in fine
            )
            for half, child in enumerate(medium_coefficients):
                difference = independent_subtract_exact_cubics(
                    independent_restrict_exact_cubic_to_half(
                        outer_coefficients, half=half
                    ),
                    child,
                )
                _hash_coefficients(outer_hasher, outer_count, difference)
                outer_count += 1
                interval, reached, nodes, hit_limit = _certify_one(
                    difference, max_depth=depth_limit, budget=budget
                )
                outer_nodes += nodes
                outer_depth = max(outer_depth, reached)
                limited = limited or hit_limit
                outer_lower = (
                    interval.lower
                    if outer_lower is None
                    else max(outer_lower, interval.lower)
                )
                outer_upper = (
                    interval.upper
                    if outer_upper is None
                    else max(outer_upper, interval.upper)
                )
            for medium_index, parent in enumerate(medium_coefficients):
                for half in (0, 1):
                    child_index = 2 * medium_index + half
                    difference = independent_subtract_exact_cubics(
                        independent_restrict_exact_cubic_to_half(parent, half=half),
                        fine_coefficients[child_index],
                    )
                    _hash_coefficients(fine_hasher, fine_count, difference)
                    fine_count += 1
                    interval, reached, nodes, hit_limit = _certify_one(
                        difference, max_depth=depth_limit, budget=budget
                    )
                    fine_nodes += nodes
                    fine_depth = max(fine_depth, reached)
                    limited = limited or hit_limit
                    fine_lower = (
                        interval.lower
                        if fine_lower is None
                        else max(fine_lower, interval.lower)
                    )
                    fine_upper = (
                        interval.upper
                        if fine_upper is None
                        else max(fine_upper, interval.upper)
                    )
            row_count += 1
    except _NodeBudgetExhausted as exc:
        raise IndependentExactTemporalArithmeticResourceExhausted(
            IndependentExactResourceExhaustionEvidence(
                reason="max_nodes_exhausted_before_stream_completed",
                max_depth=depth_limit,
                max_nodes=node_limit,
                nodes_visited=budget.used,
                completed_rows=row_count,
                outer_interval=_optional_interval(outer_lower, outer_upper),
                finest_interval=_optional_interval(fine_lower, fine_upper),
                outer_coefficient_prefix_sha256=outer_hasher.hexdigest(),
                finest_coefficient_prefix_sha256=fine_hasher.hexdigest(),
            )
        ) from exc

    if row_count == 0:
        raise ValueError("rows must contain at least one refinement row")
    assert outer_lower is not None and outer_upper is not None
    assert fine_lower is not None and fine_upper is not None
    outer_interval = IndependentExactRationalInterval(outer_lower, outer_upper)
    fine_interval = IndependentExactRationalInterval(fine_lower, fine_upper)
    pass_left = 8 * fine_interval.upper**2
    pass_right = outer_interval.lower**2
    if outer_interval.upper == 0 and fine_interval.upper == 0:
        classification = "exact_zero"
        threshold_resolved = False
        threshold_passed: bool | None = None
    elif pass_left <= pass_right:
        classification = "sufficient_contraction_pass"
        threshold_resolved = True
        threshold_passed = True
    elif 8 * fine_interval.lower**2 > outer_interval.upper**2:
        classification = "sufficient_contraction_failure"
        threshold_resolved = True
        threshold_passed = False
    else:
        raise IndependentExactTemporalArithmeticResourceExhausted(
            IndependentExactResourceExhaustionEvidence(
                reason=(
                    "max_nodes_exhausted_before_contraction_resolved"
                    if limited
                    else "max_depth_reached_before_contraction_resolved"
                ),
                max_depth=depth_limit,
                max_nodes=node_limit,
                nodes_visited=budget.used,
                completed_rows=row_count,
                outer_interval=outer_interval,
                finest_interval=fine_interval,
                outer_coefficient_prefix_sha256=outer_hasher.hexdigest(),
                finest_coefficient_prefix_sha256=fine_hasher.hexdigest(),
            )
        )

    outer_evidence = IndependentExactSupremumEvidence(
        interval=outer_interval,
        polynomial_count=outer_count,
        coefficient_stream_sha256=outer_hasher.hexdigest(),
        nodes_visited=outer_nodes,
        maximum_depth_reached=outer_depth,
    )
    fine_evidence = IndependentExactSupremumEvidence(
        interval=fine_interval,
        polynomial_count=fine_count,
        coefficient_stream_sha256=fine_hasher.hexdigest(),
        nodes_visited=fine_nodes,
        maximum_depth_reached=fine_depth,
    )
    combined = _combined_coefficient_hash(
        outer_hasher.hexdigest(), fine_hasher.hexdigest()
    )
    return IndependentExactTemporalArithmeticEvidence(
        outer_difference=outer_evidence,
        finest_difference=fine_evidence,
        row_count=row_count,
        classification=classification,
        sufficient_pass_certified=classification == "sufficient_contraction_pass",
        contraction_threshold_resolved=threshold_resolved,
        contraction_threshold_passed=threshold_passed,
        sufficient_pass_left=pass_left,
        sufficient_pass_right=pass_right,
        max_depth=depth_limit,
        max_nodes=node_limit,
        nodes_visited=budget.used,
        combined_coefficient_stream_sha256=combined,
    )


__all__ = [
    "IndependentExactRationalInterval",
    "IndependentExactResourceExhaustionEvidence",
    "IndependentExactSupremumEvidence",
    "IndependentExactTemporalArithmeticEvidence",
    "IndependentExactTemporalArithmeticResourceExhausted",
    "TDG9_INDEPENDENT_EVALUATOR_ID",
    "assess_exact_temporal_refinement_independently",
    "independent_exact_hermite_coefficients",
    "independent_restrict_exact_cubic_to_half",
    "independent_subtract_exact_cubics",
    "independent_three_halves_order_passes_exact",
]
