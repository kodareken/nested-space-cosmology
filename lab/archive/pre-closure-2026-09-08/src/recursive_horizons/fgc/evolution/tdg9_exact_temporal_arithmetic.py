"""Exact dyadic temporal-refinement arithmetic for the prospective TDG9 gate.

This module is an implementation instrument, not an authorization or a live
diagnostic.  It consumes a stream of supplied replay binary64 cubic-Hermite
data and reconstructs every input with :meth:`float.as_integer_ratio`.  A
later artifact must own the provenance and bitwise identity of that replay.
This module does not import the historical TDG5/TDG6 arithmetic and never
advances a PDE state.

One input row is ``(outer, medium, fine)``.  ``outer`` contains one segment,
``medium`` two segments, and ``fine`` four segments.  A segment is the plain
five-tuple ``(left_value, left_rhs, right_value, right_rhs, width)``.  The
coarse polynomial is restricted exactly to each child interval before the
child polynomial is subtracted.  Hence a row contributes two ``D01`` cubics
and four ``D12`` cubics.

The primary evaluator converts each exact cubic to Bernstein form and uses
de Casteljau subdivision to enclose its absolute supremum.  Bounds are exact
rationals enclosing the generally algebraic supremum; the code does not call
an irrational extremum "an exact rational supremum".  Work is bounded by the
explicit ``max_depth`` and ``max_nodes`` arguments.  An unresolved replay
contraction comparison raises typed resource-exhaustion evidence instead of
weakening the unchanged ``8 U12^2 <= L01^2`` sufficient test.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from hashlib import sha256
import heapq
import math
from numbers import Rational
from typing import Final, Iterable, Sequence


Q = Fraction
TDG9_EXACT_ARITHMETIC_SCHEMA_VERSION: Final[int] = 1
TDG9_EXACT_ARITHMETIC_EVALUATOR_ID: Final[str] = "tdg9_exact_bernstein_de_casteljau_v1"
TDG9_ORDER_SQUARED_MULTIPLIER: Final[int] = 8
_COEFFICIENT_HASH_DOMAIN: Final[bytes] = b"TDG9-EXACT-CUBIC-STREAM-v1\n"
_COMBINED_HASH_DOMAIN: Final[str] = "TDG9-EXACT-COMBINED-v1\n"
_RESOURCE_EXHAUSTION_REASONS: Final[frozenset[str]] = frozenset(
    {
        "max_nodes_exhausted_before_stream_completed",
        "max_nodes_exhausted_before_contraction_resolved",
        "max_depth_reached_before_contraction_resolved",
    }
)

Cubic = tuple[Fraction, Fraction, Fraction, Fraction]


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


def _exact_nonnegative(value: Rational, *, name: str) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, Rational):
        raise TypeError(f"{name} must be rational")
    answer = Q(value)
    if answer < 0:
        raise ValueError(f"{name} must be nonnegative")
    return answer


def _positive_integer(value: object, *, name: str, allow_zero: bool = False) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{name} must be an integer")
    if value < (0 if allow_zero else 1):
        qualifier = "nonnegative" if allow_zero else "positive"
        raise ValueError(f"{name} must be {qualifier}")
    return value


def _binary64_fraction(value: object, *, name: str) -> Fraction:
    """Return the exact rational encoded by one built-in binary64 value."""

    if type(value) is not float:  # Deliberately reject bool, int, Decimal, numpy.
        raise TypeError(f"{name} must be a built-in binary64 float")
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite")
    numerator, denominator = value.as_integer_ratio()
    return Q(numerator, denominator)


def _fixed_sequence(value: object, *, length: int, name: str) -> tuple[object, ...]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise TypeError(f"{name} must be a sequence")
    answer = tuple(value)
    if len(answer) != length:
        raise ValueError(f"{name} must contain exactly {length} entries")
    return answer


def _parse_segment(value: object, *, name: str) -> tuple[Fraction, ...]:
    raw = _fixed_sequence(value, length=5, name=name)
    exact = tuple(
        _binary64_fraction(item, name=f"{name}[{index}]")
        for index, item in enumerate(raw)
    )
    if exact[4] <= 0:
        raise ValueError(f"{name} width must be positive")
    return exact


def exact_hermite_coefficients(segment: object) -> Cubic:
    """Reconstruct one normalized cubic from exact binary64 Hermite data."""

    left, left_rhs, right, right_rhs, width = _parse_segment(segment, name="segment")
    left_slope = width * left_rhs
    right_slope = width * right_rhs
    return (
        left,
        left_slope,
        -3 * left - 2 * left_slope + 3 * right - right_slope,
        2 * left + left_slope - 2 * right + right_slope,
    )


def _cubic(value: object, *, name: str) -> Cubic:
    raw = _fixed_sequence(value, length=4, name=name)
    answer: list[Fraction] = []
    for index, item in enumerate(raw):
        if isinstance(item, bool) or not isinstance(item, Rational):
            raise TypeError(f"{name}[{index}] must be rational")
        answer.append(Q(item))
    return tuple(answer)  # type: ignore[return-value]


def restrict_exact_cubic_to_half(coefficients: object, *, half: int) -> Cubic:
    """Restrict a monomial cubic to one half and renormalize to ``[0, 1]``."""

    values = _cubic(coefficients, name="coefficients")
    if isinstance(half, bool) or half not in (0, 1):
        raise ValueError("half must be zero or one")
    a0, a1, a2, a3 = values
    if half == 0:
        return (a0, a1 / 2, a2 / 4, a3 / 8)
    return (
        a0 + a1 / 2 + a2 / 4 + a3 / 8,
        a1 / 2 + a2 / 2 + 3 * a3 / 8,
        a2 / 4 + 3 * a3 / 8,
        a3 / 8,
    )


def subtract_exact_cubics(left: object, right: object) -> Cubic:
    """Subtract two exact monomial cubics coefficient by coefficient."""

    left_values = _cubic(left, name="left")
    right_values = _cubic(right, name="right")
    return tuple(a - b for a, b in zip(left_values, right_values, strict=True))  # type: ignore[return-value]


def three_halves_order_passes_exact(
    *, outer_lower_squared: Rational, finest_upper_squared: Rational
) -> bool:
    """Apply the unchanged exact sufficient test ``8 U12^2 <= L01^2``."""

    outer = _exact_nonnegative(outer_lower_squared, name="outer_lower_squared")
    finest = _exact_nonnegative(finest_upper_squared, name="finest_upper_squared")
    return TDG9_ORDER_SQUARED_MULTIPLIER * finest <= outer


@dataclass(frozen=True, slots=True)
class ExactRationalInterval:
    """Rational isolating enclosure of a nonnegative algebraic magnitude."""

    lower: Fraction
    upper: Fraction

    def __post_init__(self) -> None:
        lower = _exact_nonnegative(self.lower, name="lower")
        upper = _exact_nonnegative(self.upper, name="upper")
        if lower > upper:
            raise ValueError("interval lower bound exceeds upper bound")
        object.__setattr__(self, "lower", lower)
        object.__setattr__(self, "upper", upper)

    @property
    def width(self) -> Fraction:
        return self.upper - self.lower


@dataclass(frozen=True, slots=True)
class ExactSupremumEvidence:
    """Stream-wide rational enclosure and its exact provenance summary."""

    interval: ExactRationalInterval
    polynomial_count: int
    coefficient_stream_sha256: str
    nodes_visited: int
    maximum_depth_reached: int
    evaluator_id: str = TDG9_EXACT_ARITHMETIC_EVALUATOR_ID

    def __post_init__(self) -> None:
        if not isinstance(self.interval, ExactRationalInterval):
            raise TypeError("interval must be ExactRationalInterval")
        _positive_integer(self.polynomial_count, name="polynomial_count")
        _positive_integer(self.nodes_visited, name="nodes_visited")
        _positive_integer(
            self.maximum_depth_reached,
            name="maximum_depth_reached",
            allow_zero=True,
        )
        _sha256_digest(
            self.coefficient_stream_sha256,
            name="coefficient stream hash",
        )
        if self.evaluator_id != TDG9_EXACT_ARITHMETIC_EVALUATOR_ID:
            raise ValueError("unexpected exact-arithmetic evaluator identifier")


@dataclass(frozen=True, slots=True)
class ExactTemporalArithmeticEvidence:
    """Completed exact D01/D12 Hermite-replay assessment for one channel.

    The classification concerns only the finite-dimensional replay cubics.  It
    is not a certification of historical proposal order, PDE convergence, or
    trajectory asymptotics.
    """

    outer_difference: ExactSupremumEvidence
    finest_difference: ExactSupremumEvidence
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
    schema_version: int = TDG9_EXACT_ARITHMETIC_SCHEMA_VERSION
    evaluator_id: str = TDG9_EXACT_ARITHMETIC_EVALUATOR_ID
    absolute_tolerance_used: bool = False
    physical_signal_used_for_normalization: bool = False
    replay_identity_certified: bool = False
    historical_proposal_certified: bool = False
    pde_trajectory_order_certified: bool = False

    def __post_init__(self) -> None:
        row_count = _positive_integer(self.row_count, name="row_count")
        max_depth = _positive_integer(self.max_depth, name="max_depth", allow_zero=True)
        max_nodes = _positive_integer(self.max_nodes, name="max_nodes")
        nodes_visited = _positive_integer(self.nodes_visited, name="nodes_visited")
        if not isinstance(
            self.outer_difference, ExactSupremumEvidence
        ) or not isinstance(self.finest_difference, ExactSupremumEvidence):
            raise TypeError("difference evidence must use ExactSupremumEvidence")
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
            raise ValueError("visited nodes exceed the exact-arithmetic budget")
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
            raise ValueError("unknown exact temporal classification")
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
            self.schema_version != TDG9_EXACT_ARITHMETIC_SCHEMA_VERSION
            or self.evaluator_id != TDG9_EXACT_ARITHMETIC_EVALUATOR_ID
            or self.absolute_tolerance_used
            or self.physical_signal_used_for_normalization
            or self.replay_identity_certified
            or self.historical_proposal_certified
            or self.pde_trajectory_order_certified
        ):
            raise ValueError("exact temporal evidence crossed its scope")
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
class ExactResourceExhaustionEvidence:
    """Typed fail-closed record for bounded exact work that did not resolve."""

    reason: str
    max_depth: int
    max_nodes: int
    nodes_visited: int
    completed_rows: int
    outer_interval: ExactRationalInterval | None
    finest_interval: ExactRationalInterval | None
    outer_coefficient_prefix_sha256: str
    finest_coefficient_prefix_sha256: str
    evaluator_id: str = TDG9_EXACT_ARITHMETIC_EVALUATOR_ID

    def __post_init__(self) -> None:
        _positive_integer(self.max_depth, name="max_depth", allow_zero=True)
        max_nodes = _positive_integer(self.max_nodes, name="max_nodes")
        nodes_visited = _positive_integer(
            self.nodes_visited, name="nodes_visited", allow_zero=True
        )
        completed_rows = _positive_integer(
            self.completed_rows, name="completed_rows", allow_zero=True
        )
        if self.reason not in _RESOURCE_EXHAUSTION_REASONS:
            raise ValueError("unknown exact-arithmetic exhaustion reason")
        if nodes_visited > max_nodes:
            raise ValueError("visited nodes exceed the exact-arithmetic budget")
        for interval, name in (
            (self.outer_interval, "outer_interval"),
            (self.finest_interval, "finest_interval"),
        ):
            if interval is not None and not isinstance(interval, ExactRationalInterval):
                raise TypeError(f"{name} must be ExactRationalInterval or None")
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
        if self.evaluator_id != TDG9_EXACT_ARITHMETIC_EVALUATOR_ID:
            raise ValueError("unexpected exact-arithmetic evaluator identifier")


class ExactTemporalArithmeticResourceExhausted(RuntimeError):
    """Raised when exact bounded work cannot certify pass or failure."""

    def __init__(self, evidence: ExactResourceExhaustionEvidence) -> None:
        self.evidence = evidence
        super().__init__(f"tdg9_exact_resource_exhausted:{evidence.reason}")


class _NodeBudgetExhausted(Exception):
    """Private control-flow sentinel owned only by the Bernstein evaluator."""


def _bernstein_controls(coefficients: Cubic) -> Cubic:
    a0, a1, a2, a3 = coefficients
    return (
        a0,
        a0 + a1 / 3,
        a0 + 2 * a1 / 3 + a2 / 3,
        a0 + a1 + a2 + a3,
    )


def _split_bernstein(controls: Cubic) -> tuple[Cubic, Cubic]:
    b0, b1, b2, b3 = controls
    b01 = (b0 + b1) / 2
    b12 = (b1 + b2) / 2
    b23 = (b2 + b3) / 2
    b012 = (b01 + b12) / 2
    b123 = (b12 + b23) / 2
    midpoint = (b012 + b123) / 2
    return (b0, b01, b012, midpoint), (midpoint, b123, b23, b3)


def _control_interval(controls: Cubic) -> ExactRationalInterval:
    return ExactRationalInterval(
        max(abs(controls[0]), abs(controls[3])),
        max(abs(item) for item in controls),
    )


@dataclass(slots=True)
class _Budget:
    maximum: int
    used: int = 0

    def claim(self, count: int) -> bool:
        if self.used + count > self.maximum:
            return False
        self.used += count
        return True


def _certify_cubic(
    coefficients: Cubic, *, max_depth: int, budget: _Budget
) -> tuple[ExactRationalInterval, int, int, bool]:
    starting_nodes = budget.used
    if not budget.claim(1):
        raise _NodeBudgetExhausted
    controls = _bernstein_controls(coefficients)
    initial = _control_interval(controls)
    lower = initial.lower
    serial = 0
    # Heap key is exact and deterministic: largest upper first, then depth/order.
    leaves: list[tuple[Fraction, int, int, Cubic]] = [
        (-initial.upper, 0, serial, controls)
    ]
    maximum_depth_reached = 0
    budget_limited = False
    while leaves:
        upper = -leaves[0][0]
        if upper == lower:
            break
        _, depth, leaf_serial, current = heapq.heappop(leaves)
        maximum_depth_reached = max(maximum_depth_reached, depth)
        if depth >= max_depth:
            heapq.heappush(leaves, (-upper, depth, leaf_serial, current))
            break
        if not budget.claim(2):
            heapq.heappush(leaves, (-upper, depth, leaf_serial, current))
            budget_limited = True
            break
        left, right = _split_bernstein(current)
        for child in (left, right):
            serial += 1
            child_interval = _control_interval(child)
            lower = max(lower, child_interval.lower)
            heapq.heappush(
                leaves,
                (-child_interval.upper, depth + 1, serial, child),
            )
        maximum_depth_reached = max(maximum_depth_reached, depth + 1)
    upper = max((-item[0] for item in leaves), default=lower)
    return (
        ExactRationalInterval(lower, upper),
        maximum_depth_reached,
        budget.used - starting_nodes,
        budget_limited,
    )


def _feed_coefficient_hash(hasher: object, ordinal: int, coefficients: Cubic) -> None:
    line = (
        f"{ordinal}|"
        + "|".join(f"{item.numerator}/{item.denominator}" for item in coefficients)
        + "\n"
    )
    hasher.update(line.encode("ascii"))  # type: ignore[attr-defined]


def _partial_interval(
    lower: Fraction | None, upper: Fraction | None
) -> ExactRationalInterval | None:
    if lower is None or upper is None:
        return None
    return ExactRationalInterval(lower, upper)


def assess_exact_temporal_refinement(
    rows: Iterable[object], *, max_depth: int, max_nodes: int
) -> ExactTemporalArithmeticEvidence:
    """Certify a streamed D01/D12 replay contraction or fail closed.

    The input iterable is consumed once and retained only one row at a time.
    ``max_nodes`` counts Bernstein tree nodes across the complete stream.
    """

    depth_limit = _positive_integer(max_depth, name="max_depth", allow_zero=True)
    node_limit = _positive_integer(max_nodes, name="max_nodes")
    budget = _Budget(node_limit)
    outer_hasher = sha256(_COEFFICIENT_HASH_DOMAIN)
    finest_hasher = sha256(_COEFFICIENT_HASH_DOMAIN)
    outer_lower: Fraction | None = None
    outer_upper: Fraction | None = None
    finest_lower: Fraction | None = None
    finest_upper: Fraction | None = None
    outer_count = 0
    finest_count = 0
    row_count = 0
    maximum_depth_reached_outer = 0
    maximum_depth_reached_finest = 0
    outer_nodes = 0
    finest_nodes = 0
    budget_limited = False

    try:
        for row_index, row in enumerate(rows):
            raw_outer, raw_medium, raw_fine = _fixed_sequence(
                row, length=3, name=f"rows[{row_index}]"
            )
            medium = _fixed_sequence(
                raw_medium, length=2, name=f"rows[{row_index}].medium"
            )
            fine = _fixed_sequence(raw_fine, length=4, name=f"rows[{row_index}].fine")
            outer_segment = _parse_segment(raw_outer, name=f"rows[{row_index}].outer")
            medium_segments = tuple(
                _parse_segment(item, name=f"rows[{row_index}].medium[{index}]")
                for index, item in enumerate(medium)
            )
            fine_segments = tuple(
                _parse_segment(item, name=f"rows[{row_index}].fine[{index}]")
                for index, item in enumerate(fine)
            )
            width = outer_segment[4]
            if any(item[4] != width / 2 for item in medium_segments):
                raise ValueError("medium widths must equal half the outer width")
            if any(item[4] != width / 4 for item in fine_segments):
                raise ValueError(
                    "fine widths must equal one quarter of the outer width"
                )

            # Reuse only exact coefficients computed in this module.  No TDG5/6
            # coefficient or rounded envelope enters this evaluator.
            outer_coefficients = exact_hermite_coefficients(raw_outer)
            medium_coefficients = tuple(
                exact_hermite_coefficients(item) for item in medium
            )
            fine_coefficients = tuple(exact_hermite_coefficients(item) for item in fine)
            for half, child in enumerate(medium_coefficients):
                difference = subtract_exact_cubics(
                    restrict_exact_cubic_to_half(outer_coefficients, half=half),
                    child,
                )
                _feed_coefficient_hash(outer_hasher, outer_count, difference)
                outer_count += 1
                interval, reached, nodes, limited = _certify_cubic(
                    difference, max_depth=depth_limit, budget=budget
                )
                outer_nodes += nodes
                maximum_depth_reached_outer = max(maximum_depth_reached_outer, reached)
                budget_limited = budget_limited or limited
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
                    difference = subtract_exact_cubics(
                        restrict_exact_cubic_to_half(parent, half=half),
                        fine_coefficients[child_index],
                    )
                    _feed_coefficient_hash(finest_hasher, finest_count, difference)
                    finest_count += 1
                    interval, reached, nodes, limited = _certify_cubic(
                        difference, max_depth=depth_limit, budget=budget
                    )
                    finest_nodes += nodes
                    maximum_depth_reached_finest = max(
                        maximum_depth_reached_finest, reached
                    )
                    budget_limited = budget_limited or limited
                    finest_lower = (
                        interval.lower
                        if finest_lower is None
                        else max(finest_lower, interval.lower)
                    )
                    finest_upper = (
                        interval.upper
                        if finest_upper is None
                        else max(finest_upper, interval.upper)
                    )
            row_count += 1
    except _NodeBudgetExhausted as exc:
        evidence = ExactResourceExhaustionEvidence(
            reason="max_nodes_exhausted_before_stream_completed",
            max_depth=depth_limit,
            max_nodes=node_limit,
            nodes_visited=budget.used,
            completed_rows=row_count,
            outer_interval=_partial_interval(outer_lower, outer_upper),
            finest_interval=_partial_interval(finest_lower, finest_upper),
            outer_coefficient_prefix_sha256=outer_hasher.hexdigest(),
            finest_coefficient_prefix_sha256=finest_hasher.hexdigest(),
        )
        raise ExactTemporalArithmeticResourceExhausted(evidence) from exc

    if row_count == 0:
        raise ValueError("rows must contain at least one refinement row")
    assert outer_lower is not None and outer_upper is not None
    assert finest_lower is not None and finest_upper is not None
    outer_interval = ExactRationalInterval(outer_lower, outer_upper)
    finest_interval = ExactRationalInterval(finest_lower, finest_upper)
    outer_evidence = ExactSupremumEvidence(
        interval=outer_interval,
        polynomial_count=outer_count,
        coefficient_stream_sha256=outer_hasher.hexdigest(),
        nodes_visited=outer_nodes,
        maximum_depth_reached=maximum_depth_reached_outer,
    )
    finest_evidence = ExactSupremumEvidence(
        interval=finest_interval,
        polynomial_count=finest_count,
        coefficient_stream_sha256=finest_hasher.hexdigest(),
        nodes_visited=finest_nodes,
        maximum_depth_reached=maximum_depth_reached_finest,
    )
    pass_left = 8 * finest_interval.upper**2
    pass_right = outer_interval.lower**2
    if outer_interval.upper == 0 and finest_interval.upper == 0:
        classification = "exact_zero"
        threshold_resolved = False
        threshold_passed: bool | None = None
    elif pass_left <= pass_right:
        classification = "sufficient_contraction_pass"
        threshold_resolved = True
        threshold_passed = True
    elif 8 * finest_interval.lower**2 > outer_interval.upper**2:
        # This strict reverse enclosure is a sufficient failure.  Merely
        # negating the sufficient-pass inequality would not be certified.
        classification = "sufficient_contraction_failure"
        threshold_resolved = True
        threshold_passed = False
    else:
        reason = (
            "max_nodes_exhausted_before_contraction_resolved"
            if budget_limited
            else "max_depth_reached_before_contraction_resolved"
        )
        raise ExactTemporalArithmeticResourceExhausted(
            ExactResourceExhaustionEvidence(
                reason=reason,
                max_depth=depth_limit,
                max_nodes=node_limit,
                nodes_visited=budget.used,
                completed_rows=row_count,
                outer_interval=outer_interval,
                finest_interval=finest_interval,
                outer_coefficient_prefix_sha256=outer_hasher.hexdigest(),
                finest_coefficient_prefix_sha256=finest_hasher.hexdigest(),
            )
        )

    combined_hash = _combined_coefficient_hash(
        outer_hasher.hexdigest(), finest_hasher.hexdigest()
    )
    return ExactTemporalArithmeticEvidence(
        outer_difference=outer_evidence,
        finest_difference=finest_evidence,
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
        combined_coefficient_stream_sha256=combined_hash,
    )


__all__ = [
    "ExactRationalInterval",
    "ExactResourceExhaustionEvidence",
    "ExactSupremumEvidence",
    "ExactTemporalArithmeticEvidence",
    "ExactTemporalArithmeticResourceExhausted",
    "TDG9_EXACT_ARITHMETIC_EVALUATOR_ID",
    "TDG9_EXACT_ARITHMETIC_SCHEMA_VERSION",
    "assess_exact_temporal_refinement",
    "exact_hermite_coefficients",
    "restrict_exact_cubic_to_half",
    "subtract_exact_cubics",
    "three_halves_order_passes_exact",
]
