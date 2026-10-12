"""Binary64 first-tangent scalars and two-jets for the exact REF1 code path.

``FloatJet2`` deliberately subclasses the established :class:`Jet2`.  The
existing tensor evaluator therefore sees the same six field jets and executes
the same algebra, while this adapter preserves finite binary64 scalars instead
of coercing them to ``Fraction``.  No field equation is copied here.

The adapter is an implementation tool, not an interval proof.  Every use must
remain paired with exact rational fixtures and independent finite-difference
checks at the artifact boundary.
"""

from __future__ import annotations

from math import isfinite
from numbers import Integral, Real
import sys
from typing import Any

from ..exact_tangent import FirstTangent
from ..spherical_reduction import Jet2


def _finite_float(name: str, value: object) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a finite real scalar")
    result = float(value)
    if not isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


class FloatFirstTangent(FirstTangent):
    """A finite binary64 primal and one forward derivative."""

    __slots__ = ()

    def __post_init__(self) -> None:
        object.__setattr__(self, "primal", _finite_float("primal", self.primal))
        object.__setattr__(self, "tangent", _finite_float("tangent", self.tangent))

    @classmethod
    def seed(cls, value: Real, derivative: Real = 1.0) -> "FloatFirstTangent":
        return cls(_finite_float("seed value", value), _finite_float("seed derivative", derivative))

    @staticmethod
    def _coerce(value: object) -> "FloatFirstTangent" | Any:
        if isinstance(value, FloatFirstTangent):
            return value
        if isinstance(value, bool) or not isinstance(value, Real):
            return NotImplemented
        return FloatFirstTangent(_finite_float("operand", value))

    def __add__(self, other: object) -> "FloatFirstTangent" | Any:
        value = self._coerce(other)
        if value is NotImplemented:
            return NotImplemented
        return FloatFirstTangent(
            self.primal + value.primal,
            self.tangent + value.tangent,
        )

    __radd__ = __add__

    def __neg__(self) -> "FloatFirstTangent":
        return FloatFirstTangent(-self.primal, -self.tangent)

    def __sub__(self, other: object) -> "FloatFirstTangent" | Any:
        value = self._coerce(other)
        return NotImplemented if value is NotImplemented else self + (-value)

    def __rsub__(self, other: object) -> "FloatFirstTangent" | Any:
        value = self._coerce(other)
        return NotImplemented if value is NotImplemented else value - self

    def __mul__(self, other: object) -> "FloatFirstTangent" | Any:
        value = self._coerce(other)
        if value is NotImplemented:
            return NotImplemented
        return FloatFirstTangent(
            self.primal * value.primal,
            self.tangent * value.primal + self.primal * value.tangent,
        )

    __rmul__ = __mul__

    def reciprocal(self) -> "FloatFirstTangent":
        if self.primal == 0.0:
            raise ZeroDivisionError("floating first-tangent reciprocal has zero primal")
        return FloatFirstTangent(
            1.0 / self.primal,
            -self.tangent / self.primal**2,
        )

    def __truediv__(self, other: object) -> "FloatFirstTangent" | Any:
        value = self._coerce(other)
        return NotImplemented if value is NotImplemented else self * value.reciprocal()

    def __rtruediv__(self, other: object) -> "FloatFirstTangent" | Any:
        value = self._coerce(other)
        return NotImplemented if value is NotImplemented else value / self

    def __pow__(self, exponent: object) -> "FloatFirstTangent" | Any:
        if isinstance(exponent, bool) or not isinstance(exponent, Integral):
            return NotImplemented
        power = int(exponent)
        if power == 0:
            return FloatFirstTangent(1.0)
        if power < 0:
            return self.reciprocal() ** (-power)
        return FloatFirstTangent(
            self.primal**power,
            power * self.primal ** (power - 1) * self.tangent,
        )

    def _comparison_operand(self, other: object) -> float | Any:
        value = self._coerce(other)
        return NotImplemented if value is NotImplemented else value.primal

    def __eq__(self, other: object) -> bool:
        value = self._coerce(other)
        if value is NotImplemented:
            return False
        tolerance = 64.0 * sys.float_info.epsilon
        primal_scale = max(1.0, abs(self.primal), abs(value.primal))
        tangent_scale = max(1.0, abs(self.tangent), abs(value.tangent))
        return (
            abs(self.primal - value.primal) <= tolerance * primal_scale
            and abs(self.tangent - value.tangent) <= tolerance * tangent_scale
        )

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

    def __abs__(self) -> float:
        return abs(self.primal)


FloatScalar = FloatFirstTangent


def _float_scalar(name: str, value: object) -> FloatScalar:
    if isinstance(value, FloatFirstTangent):
        return value
    return FloatFirstTangent(_finite_float(name, value))


def scalar_primal(value: object) -> float:
    """Return a finite primal without silently discarding a derivative."""

    if isinstance(value, FloatFirstTangent):
        return value.primal
    return _finite_float("scalar", value)


def scalar_primal_and_tangent(value: object) -> tuple[float, float]:
    if isinstance(value, FloatFirstTangent):
        return value.primal, value.tangent
    return _finite_float("scalar", value), 0.0


class FloatJet2(Jet2):
    """A finite binary64/tangent implementation of the existing ``Jet2`` API."""

    def __post_init__(self) -> None:
        for name in ("value", "dt", "dr", "dtt", "dtr", "drr"):
            object.__setattr__(self, name, _float_scalar(name, getattr(self, name)))

    @classmethod
    def constant(cls, value: Real | FloatFirstTangent) -> "FloatJet2":
        return cls(_float_scalar("constant", value))

    @staticmethod
    def coerce(value: object) -> "FloatJet2":
        if isinstance(value, FloatJet2):
            return value
        if isinstance(value, Jet2):
            return FloatJet2(
                value.value,
                value.dt,
                value.dr,
                value.dtt,
                value.dtr,
                value.drr,
            )
        if isinstance(value, bool) or not isinstance(value, (Real, FloatFirstTangent)):
            raise TypeError("floating jet operand must be a real scalar or Jet2")
        return FloatJet2.constant(value)

    def __add__(self, other: object) -> "FloatJet2":
        value = self.coerce(other)
        return FloatJet2(
            self.value + value.value,
            self.dt + value.dt,
            self.dr + value.dr,
            self.dtt + value.dtt,
            self.dtr + value.dtr,
            self.drr + value.drr,
        )

    __radd__ = __add__

    def __neg__(self) -> "FloatJet2":
        return FloatJet2(
            -self.value,
            -self.dt,
            -self.dr,
            -self.dtt,
            -self.dtr,
            -self.drr,
        )

    def __sub__(self, other: object) -> "FloatJet2":
        return self + (-self.coerce(other))

    def __rsub__(self, other: object) -> "FloatJet2":
        return self.coerce(other) - self

    def __mul__(self, other: object) -> "FloatJet2":
        value = self.coerce(other)
        return FloatJet2(
            self.value * value.value,
            self.dt * value.value + self.value * value.dt,
            self.dr * value.value + self.value * value.dr,
            self.dtt * value.value
            + 2 * self.dt * value.dt
            + self.value * value.dtt,
            self.dtr * value.value
            + self.dt * value.dr
            + self.dr * value.dt
            + self.value * value.dtr,
            self.drr * value.value
            + 2 * self.dr * value.dr
            + self.value * value.drr,
        )

    __rmul__ = __mul__

    def compose(
        self,
        value: Real | FloatFirstTangent,
        first: Real | FloatFirstTangent,
        second: Real | FloatFirstTangent,
    ) -> "FloatJet2":
        q0 = _float_scalar("composition value", value)
        q1 = _float_scalar("composition first derivative", first)
        q2 = _float_scalar("composition second derivative", second)
        return FloatJet2(
            q0,
            q1 * self.dt,
            q1 * self.dr,
            q1 * self.dtt + q2 * self.dt * self.dt,
            q1 * self.dtr + q2 * self.dt * self.dr,
            q1 * self.drr + q2 * self.dr * self.dr,
        )

    def reciprocal(self) -> "FloatJet2":
        value = scalar_primal(self.value)
        if value == 0.0:
            raise ZeroDivisionError("floating jet reciprocal has zero value")
        return self.compose(
            1.0 / self.value,
            -1.0 / (self.value * self.value),
            2.0 / (self.value**3),
        )

    def __truediv__(self, other: object) -> "FloatJet2":
        return self * self.coerce(other).reciprocal()

    def __rtruediv__(self, other: object) -> "FloatJet2":
        return self.coerce(other) / self


def float_jet(value: Jet2, *, dtt: FloatScalar | None = None) -> FloatJet2:
    """Copy one exact or floating jet, optionally replacing its acceleration."""

    if not isinstance(value, Jet2):
        raise TypeError("value must implement the established Jet2 interface")
    return FloatJet2(
        value.value,
        value.dt,
        value.dr,
        value.dtt if dtt is None else dtt,
        value.dtr,
        value.drr,
    )
