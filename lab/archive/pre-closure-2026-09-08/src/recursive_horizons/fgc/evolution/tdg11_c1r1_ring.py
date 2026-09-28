"""Standard-library exact integer/exponent arithmetic on (1/3)Z[1/2].

This module is the C1R1 representation owner.  It does not classify a
channel, hash an IMP1 receipt, or admit a production step.

``(1/3)Z[1/2]`` is the IMP1 direct Bernstein domain: an additive dyadic
module, not a ring closed under arbitrary products.  In particular
``(1/3)*(1/3) = 1/9`` leaves the domain.  The identifier ``RingElement``
names this implementation's representative; it does not claim ring
closure.  Values are

``mantissa * 2^{exp2} / 3^{exp3}``

with ``exp3 in {0, 1}`` for the direct domain and ``exp3 = 2`` only for
squared witnesses of ``8 U12^2 <= L01^2`` in the square-extended module
``(1/9)Z[1/2]``.  Constructors cancel factors of two and cancel threes
against the denominator so the triple is a unique representative.
Instances are immutable so aliases and the shared zero/one constants
cannot be mutated in place.

Dyadic linear operations, Hermite-to-Bernstein conversion, five-sample
evaluation at ``0, 1/4, 1/2, 3/4, 1`` and de Casteljau bisection stay
inside this representation.  They do not call ``fractions.Fraction``
on the hot path.  Conversion to a reduced numerator/denominator pair is
deferred until a caller needs IMP1 reference-wire text or a Fraction.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Final, Iterable, Sequence


class DirectRingConversionError(ValueError):
    """Exact rational outside the requested dyadic submodule or (1/3)Z[1/2]."""


class RingElement:
    """Immutable unique integer-exponent representative of one domain value."""

    __slots__ = ("mantissa", "exp2", "exp3")

    def __init__(self, mantissa: int, exp2: int = 0, exp3: int = 0) -> None:
        if type(mantissa) is not int or type(exp2) is not int or type(exp3) is not int:
            raise TypeError("domain components must be built-in integers")
        if exp3 < 0 or exp3 > 2:
            raise ValueError("3-adic denominator exponent must be 0, 1, or 2")
        if mantissa == 0:
            object.__setattr__(self, "mantissa", 0)
            object.__setattr__(self, "exp2", 0)
            object.__setattr__(self, "exp3", 0)
            return
        while exp3 > 0 and mantissa % 3 == 0:
            mantissa //= 3
            exp3 -= 1
        trailing = (mantissa & -mantissa).bit_length() - 1
        if trailing:
            mantissa >>= trailing
            exp2 += trailing
        object.__setattr__(self, "mantissa", mantissa)
        object.__setattr__(self, "exp2", exp2)
        object.__setattr__(self, "exp3", exp3)

    def __setattr__(self, name: str, value: object) -> None:
        raise AttributeError("RingElement representatives are immutable")

    def __delattr__(self, name: str) -> None:
        raise AttributeError("RingElement representatives are immutable")

    def is_zero(self) -> bool:
        return self.mantissa == 0

    def is_nonnegative(self) -> bool:
        return self.mantissa >= 0

    def is_dyadic(self) -> bool:
        return self.exp3 == 0

    def in_direct_ring(self) -> bool:
        return self.exp3 in (0, 1)

    def __eq__(self, other: object) -> bool:
        if type(other) is not RingElement:
            return NotImplemented
        return (
            self.mantissa == other.mantissa
            and self.exp2 == other.exp2
            and self.exp3 == other.exp3
        )

    def __hash__(self) -> int:
        return hash((self.mantissa, self.exp2, self.exp3))

    def __repr__(self) -> str:
        return (
            f"RingElement({self.mantissa!r}, {self.exp2!r}, {self.exp3!r})"
        )

    def reduced_pair(self) -> tuple[int, int]:
        """Return the reduced numerator and positive denominator."""

        mantissa = self.mantissa
        exp2 = self.exp2
        exp3 = self.exp3
        if mantissa == 0:
            return 0, 1
        three = 1 if exp3 == 0 else 3 if exp3 == 1 else 9
        if exp2 >= 0:
            return mantissa << exp2, three
        return mantissa, (1 << -exp2) * three

    def bit_size(self) -> int:
        numerator, denominator = self.reduced_pair()
        return max(abs(numerator).bit_length(), denominator.bit_length())

    def as_fraction(self) -> Fraction:
        numerator, denominator = self.reduced_pair()
        return Fraction(numerator, denominator)

    def canonical_text(self) -> str:
        numerator, denominator = self.reduced_pair()
        return int_to_decimal_text(numerator) + "/" + int_to_decimal_text(denominator)

    def abs(self) -> RingElement:
        if self.mantissa >= 0:
            return self
        return RingElement(-self.mantissa, self.exp2, self.exp3)

    def neg(self) -> RingElement:
        if self.mantissa == 0:
            return RING_ZERO
        return RingElement(-self.mantissa, self.exp2, self.exp3)

    def add(self, other: RingElement) -> RingElement:
        return add(self, other)

    def sub(self, other: RingElement) -> RingElement:
        return add(self, other.neg())

    def mul(self, other: RingElement) -> RingElement:
        return mul(self, other)

    def square(self) -> RingElement:
        if self.exp3 > 1:
            raise ValueError(
                "squaring a 1/9-witness leaves the square-extended module"
            )
        return RingElement(
            self.mantissa * self.mantissa,
            self.exp2 + self.exp2,
            self.exp3 + self.exp3,
        )


RING_ZERO: Final[RingElement] = RingElement(0)
RING_ONE: Final[RingElement] = RingElement(1)
_DECIMAL_CHUNK_BASE: Final[int] = 1_000_000_000
_THREE_POWERS: Final[tuple[int, int, int]] = (1, 3, 9)
Cubic = tuple[RingElement, RingElement, RingElement, RingElement]


def int_to_decimal_text(value: int) -> str:
    """Format an int without Python's global decimal-digit conversion limit."""

    if type(value) is not int:
        raise TypeError("decimal conversion requires a built-in integer")
    if value == 0:
        return "0"
    sign = "-" if value < 0 else ""
    remaining = abs(value)
    chunks: list[int] = []
    base = _DECIMAL_CHUNK_BASE
    while remaining:
        remaining, group = divmod(remaining, base)
        chunks.append(group)
    digits = str(chunks.pop())
    digits += "".join(format(group, "09d") for group in reversed(chunks))
    return sign + digits


def add(left: RingElement, right: RingElement) -> RingElement:
    if left.mantissa == 0:
        return right
    if right.mantissa == 0:
        return left
    exp3 = left.exp3 if left.exp3 >= right.exp3 else right.exp3
    mantissa_left = left.mantissa * _THREE_POWERS[exp3 - left.exp3]
    mantissa_right = right.mantissa * _THREE_POWERS[exp3 - right.exp3]
    if left.exp2 == right.exp2:
        return RingElement(mantissa_left + mantissa_right, left.exp2, exp3)
    if left.exp2 > right.exp2:
        return RingElement(
            (mantissa_left << (left.exp2 - right.exp2)) + mantissa_right,
            right.exp2,
            exp3,
        )
    return RingElement(
        mantissa_left + (mantissa_right << (right.exp2 - left.exp2)),
        left.exp2,
        exp3,
    )


def mul(left: RingElement, right: RingElement) -> RingElement:
    exp3 = left.exp3 + right.exp3
    if exp3 > 2:
        raise ValueError("product would leave the square-extended module (1/9)Z[1/2]")
    return RingElement(
        left.mantissa * right.mantissa,
        left.exp2 + right.exp2,
        exp3,
    )


def mul_int(value: RingElement, coefficient: int) -> RingElement:
    if type(coefficient) is not int:
        raise TypeError("integer coefficient must be a built-in integer")
    if coefficient == 0 or value.mantissa == 0:
        return RING_ZERO
    return RingElement(value.mantissa * coefficient, value.exp2, value.exp3)


def div_pow2(value: RingElement, power: int) -> RingElement:
    if type(power) is not int or power < 0:
        raise ValueError("binary denominator power must be a nonnegative integer")
    if value.mantissa == 0 or power == 0:
        return value
    return RingElement(value.mantissa, value.exp2 - power, value.exp3)


def div3(value: RingElement) -> RingElement:
    if value.exp3 >= 2:
        raise ValueError(
            "dividing by 3 would leave the square-extended module (1/9)Z[1/2]"
        )
    return RingElement(value.mantissa, value.exp2, value.exp3 + 1)


def cmp(left: RingElement, right: RingElement) -> int:
    """Return -1, 0, or 1 as ``left ? right``, matching Fraction order."""

    if left.mantissa == 0 and right.mantissa == 0:
        return 0
    if left.mantissa == 0:
        return -1 if right.mantissa > 0 else 1
    if right.mantissa == 0:
        return 1 if left.mantissa > 0 else -1
    left_scale = left.mantissa * _THREE_POWERS[right.exp3]
    right_scale = right.mantissa * _THREE_POWERS[left.exp3]
    if left.exp2 == right.exp2:
        if left_scale < right_scale:
            return -1
        if left_scale > right_scale:
            return 1
        return 0
    if left.exp2 > right.exp2:
        shifted = left_scale << (left.exp2 - right.exp2)
        if shifted < right_scale:
            return -1
        if shifted > right_scale:
            return 1
        return 0
    shifted = right_scale << (right.exp2 - left.exp2)
    if left_scale < shifted:
        return -1
    if left_scale > shifted:
        return 1
    return 0


def max_ring(values: Iterable[RingElement]) -> RingElement:
    iterator = iter(values)
    try:
        best = next(iterator)
    except StopIteration as exc:
        raise ValueError("max_ring requires at least one value") from exc
    for item in iterator:
        if cmp(item, best) > 0:
            best = item
    return best


def _intern(value: RingElement) -> RingElement:
    if value.is_zero():
        return RING_ZERO
    if value == RING_ONE:
        return RING_ONE
    return value


def from_int(value: int) -> RingElement:
    if type(value) is not int:
        raise TypeError("direct-domain integers must be built-in int")
    return _intern(RingElement(value))


def from_fraction(
    value: Fraction,
    *,
    allow_three: bool,
    allow_nine: bool = False,
) -> RingElement:
    if type(value) is not Fraction:
        raise TypeError("from_fraction requires an exact Fraction")
    if value.denominator <= 0:
        raise ValueError("denominator must be positive")
    numerator = value.numerator
    denominator = value.denominator
    trailing = (denominator & -denominator).bit_length() - 1
    odd = denominator >> trailing
    if odd == 1:
        exp3 = 0
    elif odd == 3:
        if not allow_three:
            raise DirectRingConversionError("not an exact dyadic rational")
        exp3 = 1
    elif odd == 9:
        if not allow_nine:
            raise DirectRingConversionError("not in (1/3)Z[1/2]")
        exp3 = 2
    else:
        raise DirectRingConversionError("not in (1/3)Z[1/2]")
    return _intern(RingElement(numerator, -trailing, exp3))


def from_input(
    value: object,
    *,
    name: str,
    allow_three: bool,
    allow_nine: bool = False,
) -> RingElement:
    """Convert an IMP1-admissible scalar, or signal an unsupported fast domain."""

    if type(value) is int:
        return from_int(value)
    if type(value) is Fraction:
        if value.denominator <= 0:
            raise ValueError(f"{name} denominator must be positive")
        try:
            return from_fraction(
                value, allow_three=allow_three, allow_nine=allow_nine
            )
        except DirectRingConversionError as exc:
            raise DirectRingConversionError(f"{name}: {exc}") from exc
    raise TypeError(f"{name} must be an exact Fraction or built-in int")


def hermite_bernstein(
    y0: RingElement,
    f0: RingElement,
    y1: RingElement,
    f1: RingElement,
    width: RingElement,
) -> Cubic:
    """Exact Hermite cubic ``(y0, y0 + h f0 / 3, y1 - h f1 / 3, y1)``."""

    increment_left = div3(mul(width, f0))
    increment_right = div3(mul(width, f1))
    return (
        y0,
        add(y0, increment_left),
        add(y1, increment_right.neg()),
        y1,
    )


def split_bernstein(controls: Cubic) -> tuple[Cubic, Cubic]:
    """Exact dyadic de Casteljau restriction of one cubic to [0, 1/2] and [1/2, 1]."""

    b0, b1, b2, b3 = controls
    c01 = div_pow2(add(b0, b1), 1)
    c12 = div_pow2(add(b1, b2), 1)
    c23 = div_pow2(add(b2, b3), 1)
    c012 = div_pow2(add(c01, c12), 1)
    c123 = div_pow2(add(c12, c23), 1)
    midpoint = div_pow2(add(c012, c123), 1)
    return (b0, c01, c012, midpoint), (midpoint, c123, c23, b3)


def subtract_bernstein(left: Cubic, right: Cubic) -> Cubic:
    return (
        add(left[0], right[0].neg()),
        add(left[1], right[1].neg()),
        add(left[2], right[2].neg()),
        add(left[3], right[3].neg()),
    )


def identically_zero(controls: Cubic) -> bool:
    return (
        controls[0].is_zero()
        and controls[1].is_zero()
        and controls[2].is_zero()
        and controls[3].is_zero()
    )


def hull_upper(controls: Cubic) -> RingElement:
    return max_ring(item.abs() for item in controls)


def bernstein_samples(controls: Cubic) -> tuple[RingElement, ...]:
    """Exact values at the IMP1 five-sample grid ``0, 1/4, 1/2, 3/4, 1``."""

    b0, b1, b2, b3 = controls
    at_half = div_pow2(
        add(add(b0, b3), add(mul_int(b1, 3), mul_int(b2, 3))),
        3,
    )
    at_quarter = div_pow2(
        add(add(mul_int(b0, 27), mul_int(b1, 27)), add(mul_int(b2, 9), b3)),
        6,
    )
    at_three_quarter = div_pow2(
        add(add(b0, mul_int(b1, 9)), add(mul_int(b2, 27), mul_int(b3, 27))),
        6,
    )
    return (b0, at_quarter, at_half, at_three_quarter, b3)


def sample_lower(controls: Cubic) -> RingElement:
    samples = bernstein_samples(controls)
    if identically_zero(controls):
        if any(not item.is_zero() for item in samples):
            raise ValueError("zero Bernstein cubic produced a nonzero sample")
        return RING_ZERO
    magnitude = max_ring(item.abs() for item in samples)
    if magnitude.is_zero():
        raise ValueError("nonzero cubic vanished at all five dyadic samples")
    return magnitude


def cubic_from_inputs(
    value: object,
    *,
    name: str,
    allow_three: bool = True,
) -> Cubic:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise TypeError(f"{name} must be a sequence")
    items = tuple(value)
    if len(items) != 4:
        raise ValueError(f"{name} must contain exactly 4 entries")
    return tuple(
        from_input(
            item,
            name=f"{name}[{index}]",
            allow_three=allow_three,
        )
        for index, item in enumerate(items)
    )  # type: ignore[return-value]


__all__ = [
    "Cubic",
    "DirectRingConversionError",
    "RING_ONE",
    "RING_ZERO",
    "RingElement",
    "add",
    "bernstein_samples",
    "cmp",
    "cubic_from_inputs",
    "div3",
    "div_pow2",
    "from_fraction",
    "from_input",
    "from_int",
    "hermite_bernstein",
    "hull_upper",
    "identically_zero",
    "int_to_decimal_text",
    "max_ring",
    "mul",
    "mul_int",
    "sample_lower",
    "split_bernstein",
    "subtract_bernstein",
]
