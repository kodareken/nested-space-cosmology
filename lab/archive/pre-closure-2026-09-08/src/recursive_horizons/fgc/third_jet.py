"""Exact scalar two-coordinate third jets for future metric-derived checks.

``Jet3`` is a data/calculus carrier only.  It stores an arbitrary exact local
third jet in the frozen ``(t,r)`` coordinate chart and exposes its two-jet
truncation plus the two exact first directional two-jets.  It deliberately
does not evaluate the REF1 residual, differentiate a connection, solve an
implicit branch, or assert constraint propagation.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from numbers import Integral
from typing import Any, Mapping

from .spherical_reduction import Jet2


Q = Fraction
JET3_COMPONENT_ORDER = (
    "value",
    "dt",
    "dr",
    "dtt",
    "dtr",
    "drr",
    "dttt",
    "dttr",
    "dtrr",
    "drrr",
)
THIRD_DERIVATIVE_ORDER = ("dttt", "dttr", "dtrr", "drrr")


def _fraction(name: str, value: object) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, (Fraction, Integral)):
        raise TypeError(f"{name} must be an exact rational")
    return value if isinstance(value, Fraction) else Q(int(value))


@dataclass(frozen=True, slots=True)
class Jet3:
    """One exact scalar third jet in the declared symmetric ``t,r`` chart."""

    value: Fraction
    dt: Fraction = Q(0)
    dr: Fraction = Q(0)
    dtt: Fraction = Q(0)
    dtr: Fraction = Q(0)
    drr: Fraction = Q(0)
    dttt: Fraction = Q(0)
    dttr: Fraction = Q(0)
    dtrr: Fraction = Q(0)
    drrr: Fraction = Q(0)

    def __post_init__(self) -> None:
        for name in JET3_COMPONENT_ORDER:
            object.__setattr__(self, name, _fraction(name, getattr(self, name)))

    @classmethod
    def constant(cls, value: Fraction | int) -> "Jet3":
        """Return the exact constant third jet."""

        return cls(_fraction("value", value))

    @classmethod
    def from_mapping(cls, value: Mapping[str, object]) -> "Jet3":
        """Parse an exact mapping, permitting omitted derivative slots as zero."""

        if not isinstance(value, Mapping) or "value" not in value:
            raise ValueError("third jet mapping requires value")
        unknown = set(value) - set(JET3_COMPONENT_ORDER)
        if unknown:
            raise ValueError(f"third jet mapping has unknown keys: {sorted(unknown)}")
        return cls(**{name: _fraction(name, value.get(name, Q(0))) for name in JET3_COMPONENT_ORDER})

    def as_mapping(self) -> dict[str, Fraction]:
        """Return the full canonical exact component mapping."""

        return {name: getattr(self, name) for name in JET3_COMPONENT_ORDER}

    def as_jet2(self) -> Jet2:
        """Forget third derivatives while preserving the exact two-jet."""

        return Jet2(self.value, self.dt, self.dr, self.dtt, self.dtr, self.drr)

    def directional_jet2(self, direction: str) -> Jet2:
        """Return the exact two-jet of ``partial_direction`` of this scalar."""

        if direction == "t":
            return Jet2(self.dt, self.dtt, self.dtr, self.dttt, self.dttr, self.dtrr)
        if direction == "r":
            return Jet2(self.dr, self.dtr, self.drr, self.dttr, self.dtrr, self.drrr)
        raise ValueError("third-jet direction must be t or r")

    def partial_t_jet2(self) -> Jet2:
        """Return the two-jet of the exact ``t`` derivative."""

        return self.directional_jet2("t")

    def partial_r_jet2(self) -> Jet2:
        """Return the two-jet of the exact ``r`` derivative."""

        return self.directional_jet2("r")
