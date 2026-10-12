"""Independent exact localization route for TDG9 LOC1 rational cubics.

Unlike the primary monotone-piece route, this module enumerates derivative
roots with an exact quadratic discriminant and rational square-root brackets.
It imports no TDG9 evaluator or primary localization implementation.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import isqrt
from numbers import Rational
from typing import Final, Iterable, Mapping, Sequence


Q = Fraction
EVALUATOR_ID: Final[str] = "tdg9_loc1_quadratic_discriminant_v1"


def _q(value: object) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, Rational):
        raise TypeError("coefficient must be rational")
    return Q(value)


def _cubic(value: object) -> tuple[Fraction, Fraction, Fraction, Fraction]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence) or len(value) != 4:
        raise ValueError("coefficients must contain four rational entries")
    return tuple(_q(item) for item in value)  # type: ignore[return-value]


def _sqrt_bounds(value: Fraction, bits: int) -> tuple[Fraction, Fraction]:
    if value < 0:
        raise ValueError("square root requires nonnegative input")
    scale = 1 << bits
    floor = isqrt((value.numerator * scale * scale) // value.denominator)
    while Q((floor + 1) ** 2, scale * scale) <= value:
        floor += 1
    while Q(floor * floor, scale * scale) > value:
        floor -= 1
    lower = Q(floor, scale)
    return (lower, lower) if lower * lower == value else (lower, Q(floor + 1, scale))


def _roots(coefficients: tuple[Fraction, Fraction, Fraction, Fraction], bits: int) -> tuple[tuple[Fraction, Fraction], ...]:
    _, a1, a2, a3 = coefficients
    linear, quadratic = 2 * a2, 3 * a3
    if quadratic == 0:
        if linear == 0:
            return ()
        root = -a1 / linear
        return ((root, root),) if 0 < root < 1 else ()
    discriminant = linear * linear - 4 * quadratic * a1
    if discriminant < 0:
        return ()
    low_sqrt, high_sqrt = _sqrt_bounds(discriminant, bits)
    intervals: list[tuple[Fraction, Fraction]] = []
    for sign in (-1, 1):
        numerators = (-linear + sign * low_sqrt, -linear + sign * high_sqrt)
        values = tuple(item / (2 * quadratic) for item in numerators)
        lo, hi = min(values), max(values)
        if hi <= 0 or lo >= 1:
            continue
        intervals.append((max(Q(0), lo), min(Q(1), hi)))
    return tuple(sorted(set(intervals)))


def _value_interval(coefficients: tuple[Fraction, ...], lo: Fraction, hi: Fraction) -> tuple[Fraction, Fraction]:
    lower = upper = Q(0)
    for power, coefficient in enumerate(coefficients):
        values = coefficient * lo**power, coefficient * hi**power
        lower += min(values)
        upper += max(values)
    return lower, upper


def _absolute(lower: Fraction, upper: Fraction) -> tuple[Fraction, Fraction]:
    if lower <= 0 <= upper:
        return Q(0), max(-lower, upper)
    return min(abs(lower), abs(upper)), max(abs(lower), abs(upper))


@dataclass(frozen=True, slots=True)
class IndependentLocalCubic:
    coefficients: tuple[Fraction, Fraction, Fraction, Fraction]
    metadata: tuple[tuple[str, object], ...]

    def __init__(self, coefficients: object, metadata: Mapping[str, object] | Sequence[tuple[str, object]]) -> None:
        object.__setattr__(self, "coefficients", _cubic(coefficients))
        pairs = tuple(metadata.items()) if isinstance(metadata, Mapping) else tuple(metadata)
        if len({key for key, _ in pairs}) != len(pairs):
            raise ValueError("metadata keys must be unique")
        object.__setattr__(self, "metadata", pairs)


@dataclass(frozen=True, slots=True)
class IndependentCandidate:
    polynomial_ordinal: int
    location: str
    location_ordinal: int
    parameter_lower: Fraction
    parameter_upper: Fraction
    absolute_lower: Fraction
    absolute_upper: Fraction
    metadata: tuple[tuple[str, object], ...]


@dataclass(frozen=True, slots=True)
class IndependentLocalizationEvidence:
    classification: str
    polynomial_count: int
    candidate_count: int
    candidates: tuple[IndependentCandidate, ...]
    global_absolute_lower: Fraction
    global_absolute_upper: Fraction
    maximum_candidates: int
    refinement_bits: int
    evaluator_id: str = EVALUATOR_ID
    tolerance_used: bool = False


def localize_absolute_maximum_independently(
    cubics: Iterable[IndependentLocalCubic],
    *,
    maximum_candidates: int,
    refinement_bits: int = 192,
) -> IndependentLocalizationEvidence:
    if isinstance(maximum_candidates, bool) or not isinstance(maximum_candidates, int) or maximum_candidates < 1:
        raise ValueError("maximum_candidates must be positive")
    candidates: list[IndependentCandidate] = []
    count = 0
    for ordinal, item in enumerate(cubics):
        if not isinstance(item, IndependentLocalCubic):
            raise TypeError("cubics must contain IndependentLocalCubic values")
        count += 1
        stationary = tuple((lo, hi) for lo, hi in _roots(item.coefficients, refinement_bits) if lo < hi or 0 < lo < 1)
        locations = [("left_endpoint", 0, Q(0), Q(0)), ("right_endpoint", 0, Q(1), Q(1))]
        locations.extend(("interior_stationary", root_ordinal, lo, hi) for root_ordinal, (lo, hi) in enumerate(stationary))
        for location, location_ordinal, lo, hi in locations:
            value = _value_interval(item.coefficients, lo, hi)
            absolute = _absolute(*value)
            candidates.append(IndependentCandidate(ordinal, location, location_ordinal, lo, hi, *absolute, item.metadata))
            if len(candidates) > maximum_candidates:
                raise RuntimeError("candidate_ceiling_exhausted")
    if count == 0:
        raise ValueError("at least one cubic is required")
    lower = max(item.absolute_lower for item in candidates)
    survivors = tuple(item for item in candidates if item.absolute_upper >= lower)
    upper = max(item.absolute_upper for item in survivors)
    classification = "unique_maximum" if len(survivors) == 1 else "nonunique_or_interval_inconclusive"
    return IndependentLocalizationEvidence(classification, count, len(candidates), survivors, lower, upper, maximum_candidates, refinement_bits)
