"""Exact, bounded localization of absolute extrema of rational cubics.

This module is deliberately standard-library-only.  It accepts already-built
rational cubics with caller-owned metadata and localizes the global absolute
maximum on ``[0, 1]`` to endpoints or derivative-root intervals.  It preserves
every co-maximizer whose exact enclosures cannot be separated within the
frozen refinement budget.  Binary64 ULPs are not used here as tolerances.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from numbers import Rational
from typing import Final, Iterable, Mapping, Sequence


Q = Fraction
EVALUATOR_ID: Final[str] = "tdg9_loc1_derivative_monotone_bisection_v1"


def _q(value: object, *, name: str) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, Rational):
        raise TypeError(f"{name} must be rational")
    return Q(value)


def _cubic(value: object) -> tuple[Fraction, Fraction, Fraction, Fraction]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise TypeError("coefficients must be a sequence")
    raw = tuple(value)
    if len(raw) != 4:
        raise ValueError("coefficients must contain four entries")
    return tuple(_q(item, name=f"coefficients[{index}]") for index, item in enumerate(raw))  # type: ignore[return-value]


def _poly(coefficients: tuple[Fraction, ...], x: Fraction) -> Fraction:
    answer = Q(0)
    for coefficient in reversed(coefficients):
        answer = answer * x + coefficient
    return answer


def _value_interval(
    coefficients: tuple[Fraction, Fraction, Fraction, Fraction],
    lower: Fraction,
    upper: Fraction,
) -> tuple[Fraction, Fraction]:
    """Natural rational interval extension on nonnegative ``[lower, upper]``."""

    lo = hi = Q(0)
    for power, coefficient in enumerate(coefficients):
        left = coefficient * lower**power
        right = coefficient * upper**power
        lo += min(left, right)
        hi += max(left, right)
    return lo, hi


def _absolute_interval(lower: Fraction, upper: Fraction) -> tuple[Fraction, Fraction]:
    if lower <= 0 <= upper:
        return Q(0), max(-lower, upper)
    return min(abs(lower), abs(upper)), max(abs(lower), abs(upper))


def _derivative(coefficients: tuple[Fraction, Fraction, Fraction, Fraction]) -> tuple[Fraction, Fraction, Fraction]:
    _, a1, a2, a3 = coefficients
    return a1, 2 * a2, 3 * a3


def _stationary_intervals(
    coefficients: tuple[Fraction, Fraction, Fraction, Fraction],
    *,
    depth: int,
) -> tuple[tuple[Fraction, Fraction], ...]:
    """Isolate all derivative roots in ``(0,1)`` by monotone pieces.

    The derivative is quadratic.  Splitting at its rational vertex makes each
    piece monotone, so endpoint signs isolate every simple root.  A repeated
    root is exactly the rational vertex and is retained as a point interval.
    """

    if isinstance(depth, bool) or not isinstance(depth, int) or depth < 0:
        raise ValueError("depth must be a nonnegative integer")
    derivative = _derivative(coefficients)
    b0, b1, b2 = derivative
    if b2 == 0:
        if b1 == 0:
            return ()
        root = -b0 / b1
        return ((root, root),) if 0 < root < 1 else ()
    vertex = -b1 / (2 * b2)
    cuts = [Q(0)]
    if 0 < vertex < 1:
        cuts.append(vertex)
    cuts.append(Q(1))
    roots: list[tuple[Fraction, Fraction]] = []
    if 0 < vertex < 1 and _poly(derivative, vertex) == 0:
        roots.append((vertex, vertex))
    for left, right in zip(cuts, cuts[1:]):
        f_left, f_right = _poly(derivative, left), _poly(derivative, right)
        if f_left == 0 and 0 < left < 1 and (left, left) not in roots:
            roots.append((left, left))
        if f_right == 0 and 0 < right < 1 and (right, right) not in roots:
            roots.append((right, right))
        if f_left == 0 or f_right == 0 or (f_left < 0) == (f_right < 0):
            continue
        lo, hi = left, right
        sign = f_left < 0
        for _ in range(depth):
            middle = (lo + hi) / 2
            value = _poly(derivative, middle)
            if value == 0:
                lo = hi = middle
                break
            if (value < 0) == sign:
                lo = middle
            else:
                hi = middle
        roots.append((lo, hi))
    return tuple(sorted(set(roots)))


@dataclass(frozen=True, slots=True)
class LocalCubic:
    coefficients: tuple[Fraction, Fraction, Fraction, Fraction]
    metadata: tuple[tuple[str, object], ...]

    def __init__(self, coefficients: object, metadata: Mapping[str, object] | Sequence[tuple[str, object]]) -> None:
        object.__setattr__(self, "coefficients", _cubic(coefficients))
        pairs = tuple(metadata.items()) if isinstance(metadata, Mapping) else tuple(metadata)
        if len({key for key, _ in pairs}) != len(pairs) or any(not isinstance(key, str) for key, _ in pairs):
            raise ValueError("metadata keys must be unique strings")
        object.__setattr__(self, "metadata", pairs)


@dataclass(frozen=True, slots=True)
class LocalizedCandidate:
    polynomial_ordinal: int
    location: str
    location_ordinal: int
    parameter_lower: Fraction
    parameter_upper: Fraction
    absolute_lower: Fraction
    absolute_upper: Fraction
    metadata: tuple[tuple[str, object], ...]


@dataclass(frozen=True, slots=True)
class LocalizationEvidence:
    classification: str
    polynomial_count: int
    candidate_count: int
    candidates: tuple[LocalizedCandidate, ...]
    global_absolute_lower: Fraction
    global_absolute_upper: Fraction
    maximum_candidates: int
    refinement_depth: int
    evaluator_id: str = EVALUATOR_ID
    tolerance_used: bool = False


def localize_absolute_maximum(
    cubics: Iterable[LocalCubic],
    *,
    maximum_candidates: int,
    refinement_depth: int = 160,
) -> LocalizationEvidence:
    """Return the exact-separated winner set, or every unresolved co-winner."""

    if isinstance(maximum_candidates, bool) or not isinstance(maximum_candidates, int) or maximum_candidates < 1:
        raise ValueError("maximum_candidates must be positive")
    all_candidates: list[LocalizedCandidate] = []
    polynomial_count = 0
    for ordinal, item in enumerate(cubics):
        if not isinstance(item, LocalCubic):
            raise TypeError("cubics must contain LocalCubic values")
        polynomial_count += 1
        locations = [("left_endpoint", 0, Q(0), Q(0)), ("right_endpoint", 0, Q(1), Q(1))]
        locations.extend(("interior_stationary", root_ordinal, lo, hi) for root_ordinal, (lo, hi) in enumerate(_stationary_intervals(item.coefficients, depth=refinement_depth)))
        for location, location_ordinal, lo, hi in locations:
            value_lo, value_hi = _value_interval(item.coefficients, lo, hi)
            absolute_lo, absolute_hi = _absolute_interval(value_lo, value_hi)
            all_candidates.append(LocalizedCandidate(ordinal, location, location_ordinal, lo, hi, absolute_lo, absolute_hi, item.metadata))
            if len(all_candidates) > maximum_candidates:
                raise RuntimeError("candidate_ceiling_exhausted")
    if polynomial_count == 0:
        raise ValueError("at least one cubic is required")
    global_lower = max(item.absolute_lower for item in all_candidates)
    survivors = tuple(item for item in all_candidates if item.absolute_upper >= global_lower)
    global_upper = max(item.absolute_upper for item in survivors)
    classification = "unique_maximum" if len(survivors) == 1 else "nonunique_or_interval_inconclusive"
    return LocalizationEvidence(classification, polynomial_count, len(all_candidates), survivors, global_lower, global_upper, maximum_candidates, refinement_depth)
