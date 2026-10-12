"""Exact interval first tangents for verified local derivative enclosures."""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from numbers import Integral
from typing import Any

from .exact_interval import Interval, Rational, coerce_interval


@dataclass(frozen=True, slots=True, eq=False)
class IntervalFirstTangent:
    """A primal interval and an interval enclosing one directional derivative."""

    primal: Interval
    tangent: Interval = Interval(Fraction(0), Fraction(0))

    def __post_init__(self) -> None:
        object.__setattr__(self, "primal", coerce_interval(self.primal))
        object.__setattr__(self, "tangent", coerce_interval(self.tangent))

    @classmethod
    def constant(cls, value: Interval | Rational) -> "IntervalFirstTangent":
        return cls(coerce_interval(value))

    @classmethod
    def seed(
        cls,
        value: Interval | Rational,
        derivative: Interval | Rational = Fraction(1),
    ) -> "IntervalFirstTangent":
        return cls(coerce_interval(value), coerce_interval(derivative))

    @staticmethod
    def _coerce(value: object) -> "IntervalFirstTangent" | Any:
        if isinstance(value, IntervalFirstTangent):
            return value
        try:
            return IntervalFirstTangent.constant(value)  # type: ignore[arg-type]
        except TypeError:
            return NotImplemented

    def __bool__(self) -> bool:
        raise TypeError("interval tangent truth value is ambiguous; inspect its primal interval")

    def _comparison(self, _other: object) -> Any:
        raise TypeError("interval tangent ordering is ambiguous; inspect its primal interval")

    __lt__ = _comparison
    __le__ = _comparison
    __gt__ = _comparison
    __ge__ = _comparison

    def __eq__(self, other: object) -> bool:
        value = self._coerce(other)
        return False if value is NotImplemented else self.primal == value.primal and self.tangent == value.tangent

    def __neg__(self) -> "IntervalFirstTangent":
        return IntervalFirstTangent(-self.primal, -self.tangent)

    def __add__(self, other: object) -> "IntervalFirstTangent" | Any:
        value = self._coerce(other)
        if value is NotImplemented:
            return NotImplemented
        return IntervalFirstTangent(self.primal + value.primal, self.tangent + value.tangent)

    __radd__ = __add__

    def __sub__(self, other: object) -> "IntervalFirstTangent" | Any:
        value = self._coerce(other)
        if value is NotImplemented:
            return NotImplemented
        return IntervalFirstTangent(self.primal - value.primal, self.tangent - value.tangent)

    def __rsub__(self, other: object) -> "IntervalFirstTangent" | Any:
        value = self._coerce(other)
        return NotImplemented if value is NotImplemented else value - self

    def __mul__(self, other: object) -> "IntervalFirstTangent" | Any:
        value = self._coerce(other)
        if value is NotImplemented:
            return NotImplemented
        return IntervalFirstTangent(
            self.primal * value.primal,
            self.tangent * value.primal + self.primal * value.tangent,
        )

    __rmul__ = __mul__

    def reciprocal(self) -> "IntervalFirstTangent":
        reciprocal = self.primal.reciprocal()
        return IntervalFirstTangent(reciprocal, -(self.tangent * (reciprocal**2)))

    def __truediv__(self, other: object) -> "IntervalFirstTangent" | Any:
        value = self._coerce(other)
        if value is NotImplemented:
            return NotImplemented
        return self * value.reciprocal()

    def __rtruediv__(self, other: object) -> "IntervalFirstTangent" | Any:
        value = self._coerce(other)
        return NotImplemented if value is NotImplemented else value / self

    def __pow__(self, exponent: object) -> "IntervalFirstTangent" | Any:
        if isinstance(exponent, bool) or not isinstance(exponent, Integral):
            return NotImplemented
        power = int(exponent)
        if power == 0:
            return IntervalFirstTangent.constant(1)
        if power < 0:
            return self.reciprocal() ** (-power)
        return IntervalFirstTangent(
            self.primal**power,
            power * (self.primal ** (power - 1)) * self.tangent,
        )


def primal_and_tangent(value: object) -> tuple[Interval, Interval]:
    """Return enclosed primal/tangent parts, treating a rational as zero tangent."""

    if isinstance(value, IntervalFirstTangent):
        return value.primal, value.tangent
    return coerce_interval(value), Interval.singleton(0)
