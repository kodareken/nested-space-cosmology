"""Exact first-tangent scalars for finite-dimensional residual derivatives.

``FirstTangent(x, dx)`` represents ``x + epsilon*dx`` with
``epsilon**2 = 0``.  Arithmetic therefore propagates one exact directional
derivative without floating point, finite differencing, or a symbolic-algebra
dependency.  IMP1 uses the type only as an algebraic tangent in local two-jet
space; it is not an additional spacetime derivative.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from numbers import Integral
from typing import Any


def _fraction(name: str, value: object) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, (Fraction, Integral)):
        raise TypeError(f"{name} must be an exact rational value")
    return value if isinstance(value, Fraction) else Fraction(int(value))


@dataclass(frozen=True, slots=True, eq=False)
class FirstTangent:
    """An exact rational primal value and first directional derivative."""

    primal: Fraction
    tangent: Fraction = Fraction(0)

    def __post_init__(self) -> None:
        object.__setattr__(self, "primal", _fraction("primal", self.primal))
        object.__setattr__(self, "tangent", _fraction("tangent", self.tangent))

    @classmethod
    def seed(
        cls, value: int | Fraction, derivative: int | Fraction = 1
    ) -> "FirstTangent":
        return cls(_fraction("seed value", value), _fraction("seed derivative", derivative))

    @staticmethod
    def _coerce(value: object) -> "FirstTangent" | Any:
        if isinstance(value, FirstTangent):
            return value
        if isinstance(value, bool) or not isinstance(value, (Fraction, Integral)):
            return NotImplemented
        return FirstTangent(_fraction("operand", value))

    def __add__(self, other: object) -> "FirstTangent" | Any:
        value = self._coerce(other)
        if value is NotImplemented:
            return NotImplemented
        return FirstTangent(
            self.primal + value.primal,
            self.tangent + value.tangent,
        )

    __radd__ = __add__

    def __neg__(self) -> "FirstTangent":
        return FirstTangent(-self.primal, -self.tangent)

    def __sub__(self, other: object) -> "FirstTangent" | Any:
        value = self._coerce(other)
        if value is NotImplemented:
            return NotImplemented
        return FirstTangent(
            self.primal - value.primal,
            self.tangent - value.tangent,
        )

    def __rsub__(self, other: object) -> "FirstTangent" | Any:
        value = self._coerce(other)
        if value is NotImplemented:
            return NotImplemented
        return value - self

    def __mul__(self, other: object) -> "FirstTangent" | Any:
        value = self._coerce(other)
        if value is NotImplemented:
            return NotImplemented
        return FirstTangent(
            self.primal * value.primal,
            self.tangent * value.primal + self.primal * value.tangent,
        )

    __rmul__ = __mul__

    def reciprocal(self) -> "FirstTangent":
        if self.primal == 0:
            raise ZeroDivisionError("first-tangent reciprocal has zero primal")
        return FirstTangent(
            1 / self.primal,
            -self.tangent / self.primal**2,
        )

    def __truediv__(self, other: object) -> "FirstTangent" | Any:
        value = self._coerce(other)
        if value is NotImplemented:
            return NotImplemented
        return self * value.reciprocal()

    def __rtruediv__(self, other: object) -> "FirstTangent" | Any:
        value = self._coerce(other)
        if value is NotImplemented:
            return NotImplemented
        return value / self

    def __pow__(self, exponent: object) -> "FirstTangent" | Any:
        if isinstance(exponent, bool) or not isinstance(exponent, Integral):
            return NotImplemented
        power = int(exponent)
        if power == 0:
            return FirstTangent(Fraction(1))
        if power < 0:
            return (self.reciprocal()) ** (-power)
        return FirstTangent(
            self.primal**power,
            power * self.primal ** (power - 1) * self.tangent,
        )

    def _comparison_operand(self, other: object) -> Fraction | Any:
        value = self._coerce(other)
        return NotImplemented if value is NotImplemented else value.primal

    def __eq__(self, other: object) -> bool:
        value = self._coerce(other)
        if value is NotImplemented:
            return False
        return self.primal == value.primal and self.tangent == value.tangent

    def __lt__(self, other: object) -> bool | Any:
        value = self._comparison_operand(other)
        return NotImplemented if value is NotImplemented else self.primal < value

    def __le__(self, other: object) -> bool | Any:
        value = self._comparison_operand(other)
        return NotImplemented if value is NotImplemented else self.primal <= value

    def __gt__(self, other: object) -> bool | Any:
        value = self._comparison_operand(other)
        return NotImplemented if value is NotImplemented else self.primal > value

    def __ge__(self, other: object) -> bool | Any:
        value = self._comparison_operand(other)
        return NotImplemented if value is NotImplemented else self.primal >= value


def primal_and_tangent(value: object) -> tuple[Fraction, Fraction]:
    """Return exact primal/tangent parts, treating a rational as zero tangent."""

    if isinstance(value, FirstTangent):
        return value.primal, value.tangent
    return _fraction("tangent output", value), Fraction(0)
