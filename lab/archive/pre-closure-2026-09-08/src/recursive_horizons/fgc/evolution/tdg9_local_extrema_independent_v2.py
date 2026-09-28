"""Independent exact LOC2 localization for rational cubic extrema.

This route never imports the primary localizer.  It first removes exact
derivative roots at the unit-interval endpoints, then uses exact rational
square detection or a bounded adaptive discriminant enclosure.  A nonsquare
root enclosure is never clipped into ``(0, 1)``: unresolved location or root
separation stops with typed, non-scientific evidence.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from hashlib import sha256
from math import gcd, isqrt, lcm
from numbers import Rational
from typing import Final, Iterable, Mapping, Sequence


Q = Fraction
EVALUATOR_ID: Final[str] = "tdg9_loc2_endpoint_deflated_adaptive_discriminant_v2"
INITIAL_REFINEMENT_BITS: Final[int] = 192
REFINEMENT_SCHEDULE: Final[str] = "doubling_plus_exact_proof_cap"
GLOBAL_PROOF_BIT_CEILING: Final[int] = 9216


class RootIsolationInconclusive(RuntimeError):
    """A typed diagnostic stop, never a physical or numerical classification."""

    reason: str

    def __init__(self, reason: str, detail: object) -> None:
        if reason not in {
            "root_proof_cap_exceeded",
            "root_location_inconclusive",
            "root_separation_inconclusive",
            "root_count_invariant_failure",
        }:
            raise ValueError("unknown root-isolation reason")
        self.reason = reason
        self.detail = str(detail)
        super().__init__(f"{reason}: {self.detail}")


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
    value = Q(0)
    for coefficient in reversed(coefficients):
        value = value * x + coefficient
    return value


def _value_interval(
    coefficients: tuple[Fraction, Fraction, Fraction, Fraction],
    lower: Fraction,
    upper: Fraction,
) -> tuple[Fraction, Fraction]:
    lo = hi = Q(0)
    for power, coefficient in enumerate(coefficients):
        endpoints = coefficient * lower**power, coefficient * upper**power
        lo += min(endpoints)
        hi += max(endpoints)
    return lo, hi


def _absolute_interval(lower: Fraction, upper: Fraction) -> tuple[Fraction, Fraction]:
    if lower <= 0 <= upper:
        return Q(0), max(-lower, upper)
    return min(abs(lower), abs(upper)), max(abs(lower), abs(upper))


def _perfect_rational_square(value: Fraction) -> Fraction | None:
    if value < 0:
        return None
    numerator = isqrt(value.numerator)
    denominator = isqrt(value.denominator)
    if numerator * numerator == value.numerator and denominator * denominator == value.denominator:
        return Q(numerator, denominator)
    return None


def _sqrt_bounds_integer_absolute(value: int, bits: int) -> tuple[Fraction, Fraction]:
    """Absolute dyadic enclosure with exact width ``2^-bits``."""

    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError("square root requires a nonnegative integer")
    if isinstance(bits, bool) or not isinstance(bits, int) or bits < 1:
        raise ValueError("bits must be positive")
    if value == 0:
        return Q(0), Q(0)
    scale = 1 << bits
    floor = isqrt(value << (2 * bits))
    lower = Q(floor, scale)
    if floor * floor == value << (2 * bits):
        return lower, lower
    return lower, Q(floor + 1, scale)


def _primitive_integer_quadratic(
    b0: Fraction,
    b1: Fraction,
    b2: Fraction,
) -> tuple[int, int, int]:
    """Return primitive ``A*x^2+B*x+C`` with a deterministic positive lead."""

    denominator = lcm(b0.denominator, b1.denominator, b2.denominator)
    raw = (
        b2.numerator * (denominator // b2.denominator),
        b1.numerator * (denominator // b1.denominator),
        b0.numerator * (denominator // b0.denominator),
    )
    common = gcd(gcd(abs(raw[0]), abs(raw[1])), abs(raw[2]))
    if common:
        raw = tuple(item // common for item in raw)
    sign_owner = next((item for item in raw if item != 0), 1)
    if sign_owner < 0:
        raw = tuple(-item for item in raw)
    return raw  # type: ignore[return-value]


def _proof_bits(A: int, B: int, C: int) -> int:
    height_bits = max(abs(A), abs(B), abs(C), 1).bit_length()
    return 2 * height_bits + 8


def _refinement_schedule(proof_bits: int) -> tuple[int, ...]:
    if proof_bits > GLOBAL_PROOF_BIT_CEILING:
        raise RootIsolationInconclusive(
            "root_proof_cap_exceeded",
            f"proof bits {proof_bits} exceed authorization ceiling {GLOBAL_PROOF_BIT_CEILING}",
        )
    terminal = max(INITIAL_REFINEMENT_BITS, proof_bits)
    values = [INITIAL_REFINEMENT_BITS]
    while values[-1] < terminal:
        doubled = 2 * values[-1]
        if doubled >= terminal:
            break
        values.append(doubled)
    if values[-1] != terminal:
        values.append(terminal)
    return tuple(values)


def _sign_variations(values: tuple[Fraction, ...]) -> int:
    signs = [1 if value > 0 else -1 for value in values if value != 0]
    return sum(left != right for left, right in zip(signs, signs[1:], strict=False))


def _sturm_open_unit_root_count(A: int, B: int, C: int) -> int:
    """Exact distinct-root count in ``(0,1)`` for a nondegenerate quadratic."""

    discriminant = B * B - 4 * A * C
    if discriminant <= 0:
        return 0 if discriminant < 0 else int(0 < Q(-B, 2 * A) < 1)
    tail = Q(discriminant, 4 * A)
    at_zero = (Q(C), Q(B), tail)
    at_one = (Q(A + B + C), Q(2 * A + B), tail)
    return _sign_variations(at_zero) - _sign_variations(at_one)


@dataclass(frozen=True, slots=True)
class RootInterval:
    lower: Fraction
    upper: Fraction
    method: str
    refinement_bits: int

    def __post_init__(self) -> None:
        if self.lower > self.upper:
            raise ValueError("root interval is reversed")
        if self.method not in {
            "exact_endpoint_deflation",
            "exact_linear",
            "exact_repeated",
            "exact_rational_square",
            "adaptive_nonsquare",
        }:
            raise ValueError("unknown root method")


def _exact_open_root(value: Fraction, method: str) -> tuple[RootInterval, ...]:
    if 0 < value < 1:
        return (RootInterval(value, value, method, 0),)
    return ()


def _endpoint_deflated_roots(
    C: int,
    B: int,
    A: int,
) -> tuple[RootInterval, ...] | None:
    """Return exact interior roots when q has x=0 or x=1 as an exact factor."""

    if C == 0:
        if B == 0:
            return ()
        if A == 0:
            return ()
        return _exact_open_root(Q(-B, A), "exact_endpoint_deflation")
    if A + B + C == 0:
        if A == 0:
            return ()
        return _exact_open_root(Q(C, A), "exact_endpoint_deflation")
    return None


def _quadratic_intervals(
    C: int,
    B: int,
    A: int,
    bits: int,
) -> tuple[RootInterval, ...]:
    discriminant = B * B - 4 * A * C
    if discriminant < 0:
        return ()
    if discriminant == 0:
        root = Q(-B, 2 * A)
        return (RootInterval(root, root, "exact_repeated", 0),)
    exact = _perfect_rational_square(discriminant)
    if exact is not None:
        exact_integer = exact.numerator
        roots = tuple(sorted((Q(-B - exact_integer, 2 * A), Q(-B + exact_integer, 2 * A))))
        return tuple(RootInterval(root, root, "exact_rational_square", 0) for root in roots)
    sqrt_lower, sqrt_upper = _sqrt_bounds_integer_absolute(discriminant, bits)
    intervals: list[RootInterval] = []
    for sign in (-1, 1):
        numerators = (Q(-B) + sign * sqrt_lower, Q(-B) + sign * sqrt_upper)
        values = tuple(item / (2 * A) for item in numerators)
        intervals.append(
            RootInterval(min(values), max(values), "adaptive_nonsquare", bits)
        )
    return tuple(sorted(intervals, key=lambda item: (item.lower, item.upper)))


def isolate_open_unit_roots(
    b0: object,
    b1: object,
    b2: object,
) -> tuple[RootInterval, ...]:
    """Isolate derivative roots wholly inside ``(0,1)`` or fail closed."""

    coefficients = (
        _q(b0, name="b0"),
        _q(b1, name="b1"),
        _q(b2, name="b2"),
    )
    b0_q, b1_q, b2_q = coefficients
    A, B, C = _primitive_integer_quadratic(b0_q, b1_q, b2_q)
    if A == 0:
        if B == 0:
            return ()
        return _exact_open_root(Q(-C, B), "exact_linear")
    deflated = _endpoint_deflated_roots(C, B, A)
    if deflated is not None:
        return deflated
    expected_count = _sturm_open_unit_root_count(A, B, C)
    last_location: tuple[RootInterval, ...] = ()
    schedule = _refinement_schedule(_proof_bits(A, B, C))
    for index, bits in enumerate(schedule):
        roots = _quadratic_intervals(C, B, A, bits)
        if not roots:
            return ()
        location_resolved = all(
            (root.upper <= 0)
            or (root.lower >= 1)
            or (0 < root.lower and root.upper < 1)
            for root in roots
        )
        if not location_resolved:
            last_location = roots
            continue
        if len(roots) == 2 and roots[0].upper >= roots[1].lower:
            if index + 1 < len(schedule):
                continue
            raise RootIsolationInconclusive("root_separation_inconclusive", roots)
        interior = tuple(root for root in roots if 0 < root.lower and root.upper < 1)
        if len(interior) != expected_count:
            raise RootIsolationInconclusive(
                "root_count_invariant_failure",
                {"expected": expected_count, "observed": len(interior), "roots": roots},
            )
        return interior
    raise RootIsolationInconclusive("root_location_inconclusive", last_location)


@dataclass(frozen=True, slots=True)
class IndependentLocalCubicV2:
    coefficients: tuple[Fraction, Fraction, Fraction, Fraction]
    metadata: tuple[tuple[str, object], ...]

    def __init__(self, coefficients: object, metadata: Mapping[str, object] | Sequence[tuple[str, object]]) -> None:
        object.__setattr__(self, "coefficients", _cubic(coefficients))
        pairs = tuple(metadata.items()) if isinstance(metadata, Mapping) else tuple(metadata)
        if len({key for key, _ in pairs}) != len(pairs) or any(not isinstance(key, str) for key, _ in pairs):
            raise ValueError("metadata keys must be unique strings")
        object.__setattr__(self, "metadata", pairs)


@dataclass(frozen=True, slots=True)
class IndependentCandidateV2:
    polynomial_ordinal: int
    location: str
    location_ordinal: int
    parameter_lower: Fraction
    parameter_upper: Fraction
    absolute_lower: Fraction
    absolute_upper: Fraction
    metadata: tuple[tuple[str, object], ...]
    root_method: str
    refinement_bits: int


@dataclass(frozen=True, slots=True)
class IndependentLocalizationEvidenceV2:
    classification: str
    polynomial_count: int
    candidate_count: int
    candidates: tuple[IndependentCandidateV2, ...]
    global_absolute_lower: Fraction
    global_absolute_upper: Fraction
    maximum_candidates: int
    initial_refinement_bits: int
    refinement_schedule: str
    global_proof_bit_ceiling: int
    stationary_count_stream_sha256: str
    evaluator_id: str = EVALUATOR_ID
    tolerance_used: bool = False
    boundary_clipping_used: bool = False


def localize_absolute_maximum_independently_v2(
    cubics: Iterable[IndependentLocalCubicV2],
    *,
    maximum_candidates: int,
) -> IndependentLocalizationEvidenceV2:
    if isinstance(maximum_candidates, bool) or not isinstance(maximum_candidates, int) or maximum_candidates < 1:
        raise ValueError("maximum_candidates must be positive")
    candidates: list[IndependentCandidateV2] = []
    stationary_digest = sha256(b"TDG9-LOC2-STATIONARY-COUNT-STREAM-v1\n")
    polynomial_count = 0
    for ordinal, item in enumerate(cubics):
        if not isinstance(item, IndependentLocalCubicV2):
            raise TypeError("cubics must contain IndependentLocalCubicV2 values")
        polynomial_count += 1
        _, a1, a2, a3 = item.coefficients
        roots = isolate_open_unit_roots(
            a1,
            2 * a2,
            3 * a3,
        )
        stationary_digest.update(f"{ordinal}|{len(roots)}\n".encode())
        locations: list[tuple[str, int, Fraction, Fraction, str, int]] = [
            ("left_endpoint", 0, Q(0), Q(0), "exact_endpoint", 0),
            ("right_endpoint", 0, Q(1), Q(1), "exact_endpoint", 0),
        ]
        locations.extend(
            (
                "interior_stationary",
                root_ordinal,
                root.lower,
                root.upper,
                root.method,
                root.refinement_bits,
            )
            for root_ordinal, root in enumerate(roots)
        )
        for location, location_ordinal, lower, upper, method, bits in locations:
            value_lower, value_upper = _value_interval(item.coefficients, lower, upper)
            absolute_lower, absolute_upper = _absolute_interval(value_lower, value_upper)
            candidates.append(
                IndependentCandidateV2(
                    ordinal,
                    location,
                    location_ordinal,
                    lower,
                    upper,
                    absolute_lower,
                    absolute_upper,
                    item.metadata,
                    method,
                    bits,
                )
            )
            if len(candidates) > maximum_candidates:
                raise RuntimeError("candidate_ceiling_exhausted")
    if polynomial_count == 0:
        raise ValueError("at least one cubic is required")
    global_lower = max(item.absolute_lower for item in candidates)
    survivors = tuple(item for item in candidates if item.absolute_upper >= global_lower)
    global_upper = max(item.absolute_upper for item in survivors)
    classification = "unique_maximum" if len(survivors) == 1 else "nonunique_or_interval_inconclusive"
    return IndependentLocalizationEvidenceV2(
        classification,
        polynomial_count,
        len(candidates),
        survivors,
        global_lower,
        global_upper,
        maximum_candidates,
        INITIAL_REFINEMENT_BITS,
        REFINEMENT_SCHEDULE,
        GLOBAL_PROOF_BIT_CEILING,
        stationary_digest.hexdigest(),
    )


__all__ = [
    "EVALUATOR_ID",
    "INITIAL_REFINEMENT_BITS",
    "REFINEMENT_SCHEDULE",
    "GLOBAL_PROOF_BIT_CEILING",
    "RootIsolationInconclusive",
    "RootInterval",
    "IndependentLocalCubicV2",
    "IndependentCandidateV2",
    "IndependentLocalizationEvidenceV2",
    "isolate_open_unit_roots",
    "localize_absolute_maximum_independently_v2",
]
