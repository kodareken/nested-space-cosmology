"""Prospective exact binary64 lattice for TDG7 1/2/4 subdivisions.

This module is deliberately pure arithmetic.  It neither prepares TDG6
shadows nor mutates a PROTO14 member.  It is a prospective diagnosis/design
aid for the forward PROTO14 coordinate-time envelope ``[23/16, 32]``.

The central construction is conservative: a binary64 cap is an upper bound,
not an instruction to use an unrepresentable width.  For the current/event
interval we choose ``Q = max(ulp(current), ulp(event_target))`` and floor the
cap to a positive multiple of ``G = 8 Q``.  The five boundaries are made once
from exact dyadics and then sliced into the one-, two-, and four-step paths.

``G = 4 Q`` would certify only those path *endpoints*.  Its fine width may be
an odd multiple of ``Q``, making the ``c=1/2`` time used by both RK4 and
SSPRK3 a half-``Q`` value that binary64 rounds (or aliases to an endpoint).
``G = 8 Q`` makes every fine interval a multiple of ``2 Q``.  This certifies
the actual stage abscissae ``c in {0, 1/2, 1}``, not merely their endpoints.

``G`` alone is not enough.  Both endpoints must be on the ``Q`` lattice.  An
off-lattice phase can remain harmless below a binade yet become unrepresentable
after crossing into one with a larger ULP; reducing the width cannot remove
that phase.  Such input is a typed coordinate-lattice limit, not a rounded
landing opportunity.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import math
import struct
from typing import Final


PROTO14_TIME_LOWER: Final[float] = 23.0 / 16.0
PROTO14_TIME_UPPER: Final[float] = 32.0


class CoordinateLatticeLimitReached(ValueError):
    """A requested binary64 coordinate step has no certified TDG7 plan."""

    def __init__(self, reason: str, *, detail: str = "") -> None:
        self.reason = reason
        self.classification = "coordinate_lattice_limit_reached"
        super().__init__(
            f"coordinate_lattice_limit_reached:{reason}"
            + (f": {detail}" if detail else "")
        )


def _bits(value: float) -> bytes:
    return struct.pack(">d", value)


def _finite_fraction(name: str, value: object) -> tuple[float, Fraction]:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise CoordinateLatticeLimitReached("nonfinite_coordinate", detail=name)
    try:
        answer = float(value)
    except (OverflowError, ValueError) as exc:
        raise CoordinateLatticeLimitReached(
            "nonfinite_coordinate", detail=name
        ) from exc
    if not math.isfinite(answer):
        raise CoordinateLatticeLimitReached("nonfinite_coordinate", detail=name)
    return answer, Fraction.from_float(answer)


def _exact_float(value: Fraction) -> float:
    answer = float(value)
    if not math.isfinite(answer) or Fraction.from_float(answer) != value:
        raise CoordinateLatticeLimitReached("internal_exactness_failure")
    return answer


@dataclass(frozen=True, slots=True)
class TDG7Binary64SubdivisionPlan:
    """A certified shared 1/2/4, stage-safe boundary plan with positive widths."""

    current: float
    event_target: float
    requested_cap: float
    quantum: float
    selection_budget: float
    macro_width: float
    fine_width: float
    conservative_reduction: float
    minimum_width: float | None
    boundaries: tuple[float, float, float, float, float]
    target_limited: bool

    def __post_init__(self) -> None:
        """Make ``dataclasses.replace`` obey the same exact proof contract."""

        validate_tdg7_binary64_subdivision_plan(self)

    @property
    def one_full(self) -> tuple[float, float]:
        return (self.boundaries[0], self.boundaries[4])

    @property
    def two_half(self) -> tuple[float, float, float]:
        return (self.boundaries[0], self.boundaries[2], self.boundaries[4])

    @property
    def four_quarter(self) -> tuple[float, float, float, float, float]:
        return self.boundaries

    @staticmethod
    def _stage_times_for(
        boundaries: tuple[float, ...],
    ) -> tuple[tuple[float, float, float], ...]:
        """Return ``(c=0, c=1/2, c=1)`` times for each exact substep.

        These are shared time coordinates, not a statement that an RK method
        has only three evaluations.  RK4 and SSPRK3 can revisit a stage time;
        TDG7 requires only that the distinct geometric midpoint is exact and
        strictly interior rather than silently rounded to an endpoint.
        """

        return tuple(
            (left, left + (right - left) / 2.0, right)
            for left, right in zip(boundaries, boundaries[1:])
        )

    @property
    def one_full_stage_times(self) -> tuple[tuple[float, float, float], ...]:
        return self._stage_times_for(self.one_full)

    @property
    def two_half_stage_times(self) -> tuple[tuple[float, float, float], ...]:
        return self._stage_times_for(self.two_half)

    @property
    def four_quarter_stage_times(self) -> tuple[tuple[float, float, float], ...]:
        return self._stage_times_for(self.four_quarter)

    @property
    def endpoint(self) -> float:
        return self.boundaries[4]


@dataclass(frozen=True, slots=True)
class HistoricalTDG6GuardReconstruction:
    """Read-only reconstruction of the historical independent 1/2/4 guards."""

    start: float
    requested_width: float
    expected_final: float
    nominal_steps: tuple[float, float, float]
    rounded_boundaries: tuple[
        tuple[float, float],
        tuple[float, float, float],
        tuple[float, float, float, float, float],
    ]
    rounded_differences: tuple[tuple[float, ...], tuple[float, ...], tuple[float, ...]]
    failed_counts: tuple[int, ...]


def reconstruct_historical_tdg6_guard(
    start: object, requested_width: object
) -> HistoricalTDG6GuardReconstruction:
    """Reproduce TDG6's old independent 1/2/4 boundary guard exactly.

    This intentionally uses the historical rounded expression
    ``start + index * (width / count)`` separately for each count.  It does
    not execute a shadow or touch campaign state.
    """

    start_float, _ = _finite_fraction("current", start)
    width_float, width = _finite_fraction("requested_cap", requested_width)
    if width <= 0:
        raise CoordinateLatticeLimitReached("nonpositive_requested_cap")
    expected_final = start_float + width_float
    all_boundaries: list[tuple[float, ...]] = []
    all_differences: list[tuple[float, ...]] = []
    nominal_steps: list[float] = []
    failed: list[int] = []
    for count in (1, 2, 4):
        nominal = width_float / count
        boundaries = tuple(start_float + index * nominal for index in range(count + 1))
        differences = tuple(right - left for left, right in zip(boundaries, boundaries[1:]))
        nominal_steps.append(nominal)
        all_boundaries.append(boundaries)
        all_differences.append(differences)
        if _bits(boundaries[-1]) != _bits(expected_final) or any(
            _bits(delta) != _bits(nominal) for delta in differences
        ):
            failed.append(count)
    return HistoricalTDG6GuardReconstruction(
        start=start_float,
        requested_width=width_float,
        expected_final=expected_final,
        nominal_steps=(nominal_steps[0], nominal_steps[1], nominal_steps[2]),
        rounded_boundaries=(all_boundaries[0], all_boundaries[1], all_boundaries[2]),  # type: ignore[arg-type]
        rounded_differences=(all_differences[0], all_differences[1], all_differences[2]),
        failed_counts=tuple(failed),
    )


def validate_tdg7_binary64_subdivision_plan(plan: TDG7Binary64SubdivisionPlan) -> None:
    """Fail closed unless *plan* still satisfies every exact/bitwise invariant."""

    if not isinstance(plan, TDG7Binary64SubdivisionPlan):
        raise TypeError("plan must be TDG7Binary64SubdivisionPlan")
    x_float, x = _finite_fraction("current", plan.current)
    target_float, target = _finite_fraction("event_target", plan.event_target)
    cap_float, cap = _finite_fraction("requested_cap", plan.requested_cap)
    budget_float, budget = _finite_fraction("selection_budget", plan.selection_budget)
    q_float, q = _finite_fraction("quantum", plan.quantum)
    macro_float, macro = _finite_fraction("macro_width", plan.macro_width)
    fine_float, fine = _finite_fraction("fine_width", plan.fine_width)
    reduction_float, reduction = _finite_fraction(
        "conservative_reduction", plan.conservative_reduction
    )
    if not (Fraction.from_float(PROTO14_TIME_LOWER) <= x < target <= Fraction.from_float(PROTO14_TIME_UPPER)):
        raise CoordinateLatticeLimitReached("internal_exactness_failure")
    if cap <= 0 or budget <= 0 or macro <= 0 or fine <= 0 or reduction < 0:
        raise CoordinateLatticeLimitReached("internal_exactness_failure")
    if q != Fraction.from_float(max(math.ulp(x_float), math.ulp(target_float))):
        raise CoordinateLatticeLimitReached("internal_exactness_failure")
    if x % q != 0 or target % q != 0:
        raise CoordinateLatticeLimitReached("internal_exactness_failure")
    # Eight Q ticks per macro step make every quarter-step two Q ticks, so the
    # stage midpoint of each fine step is an exact one-Q lattice point.
    grid = 8 * q
    remaining = target - x
    if remaining % grid != 0 or budget != min(cap, remaining):
        raise CoordinateLatticeLimitReached("internal_exactness_failure")
    if (
        macro != (budget // grid) * grid
        or macro % grid != 0
        or macro > budget
        or fine != macro / 4
        or fine % (2 * q) != 0
        or reduction != budget - macro
    ):
        raise CoordinateLatticeLimitReached("internal_exactness_failure")
    if plan.minimum_width is not None:
        _, minimum = _finite_fraction("minimum_width", plan.minimum_width)
        if minimum <= 0 or macro < minimum:
            raise CoordinateLatticeLimitReached("internal_exactness_failure")
    if plan.target_limited != (macro == remaining):
        raise CoordinateLatticeLimitReached("internal_exactness_failure")
    if len(plan.boundaries) != 5:
        raise CoordinateLatticeLimitReached("internal_exactness_failure")
    exact = tuple(x + index * fine for index in range(5))
    if any(
        Fraction.from_float(actual) != expected
        for actual, expected in zip(plan.boundaries, exact, strict=True)
    ) or _bits(plan.boundaries[0]) != _bits(x_float):
        raise CoordinateLatticeLimitReached("internal_exactness_failure")
    expected_stage_paths = (
        ((x, x + macro / 2, x + macro),),
        (
            (x, x + macro / 4, x + macro / 2),
            (x + macro / 2, x + 3 * macro / 4, x + macro),
        ),
        tuple(
            (
                x + index * fine,
                x + (2 * index + 1) * fine / 2,
                x + (index + 1) * fine,
            )
            for index in range(4)
        ),
    )
    stage_paths = (
        plan.one_full_stage_times,
        plan.two_half_stage_times,
        plan.four_quarter_stage_times,
    )
    if (
        Fraction.from_float(plan.boundaries[4] - plan.boundaries[0]) != macro
        or any(
            Fraction.from_float(right - left) != fine
            for left, right in zip(plan.boundaries, plan.boundaries[1:])
        )
        or _bits(plan.one_full[-1]) != _bits(plan.two_half[-1])
        or _bits(plan.two_half[-1]) != _bits(plan.four_quarter[-1])
        or any(
            len(actual_path) != len(expected_path)
            or any(
                Fraction.from_float(actual_time) != expected_time
                for actual_step, expected_step in zip(
                    actual_path, expected_path, strict=True
                )
                for actual_time, expected_time in zip(
                    actual_step, expected_step, strict=True
                )
            )
            or any(
                not (left < midpoint < right)
                for left, midpoint, right in actual_path
            )
            for actual_path, expected_path in zip(
                stage_paths, expected_stage_paths, strict=True
            )
        )
        or not all(math.isfinite(item) for item in (budget_float, q_float, macro_float, fine_float, reduction_float))
    ):
        raise CoordinateLatticeLimitReached("internal_exactness_failure")


def plan_forward_proto14_subdivision(
    current: object, event_target: object, requested_cap: object, *, minimum_width: object | None = None
) -> TDG7Binary64SubdivisionPlan:
    """Return the largest forward 1/2/4-uniform width no greater than *cap*.

    All comparisons and flooring are over the exact dyadic values represented
    by the supplied binary64 inputs.  The resulting plan is intentionally
    restricted to the forward PROTO14 envelope; this makes the proof simple
    and prevents this design helper from silently becoming a generic runtime.
    """

    x_float, x = _finite_fraction("current", current)
    t_float, target = _finite_fraction("event_target", event_target)
    cap_float, cap = _finite_fraction("requested_cap", requested_cap)
    lower = Fraction.from_float(PROTO14_TIME_LOWER)
    upper = Fraction.from_float(PROTO14_TIME_UPPER)
    if not (lower <= x < target <= upper):
        reason = "target_not_ahead" if target <= x else "outside_frozen_time_envelope"
        raise CoordinateLatticeLimitReached(reason)
    if cap <= 0:
        raise CoordinateLatticeLimitReached("nonpositive_requested_cap")
    minimum_float: float | None = None
    minimum: Fraction | None = None
    if minimum_width is not None:
        minimum_float, minimum = _finite_fraction("minimum_width", minimum_width)
        if minimum <= 0:
            raise CoordinateLatticeLimitReached("below_minimum_aligned_width")

    q_float = max(math.ulp(x_float), math.ulp(t_float))
    q = Fraction.from_float(q_float)
    # ``8 Q`` is required rather than ``4 Q`` because the actual RK4/SSPRK3
    # ``c=1/2`` stages must be exact and interior on every fine substep.
    grid = 8 * q
    # This predicate, rather than decrementing a candidate width, is what
    # rules out an unrepresentable phase at an upward binade crossing.
    if x % q != 0 or target % q != 0:
        raise CoordinateLatticeLimitReached("endpoint_not_q_aligned")
    remaining = target - x
    if remaining % grid != 0:
        raise CoordinateLatticeLimitReached("target_not_stage_lattice_aligned")

    budget = min(cap, remaining)
    ticks = budget // grid
    if ticks < 1:
        raise CoordinateLatticeLimitReached("no_positive_aligned_width")
    macro = ticks * grid
    fine = macro / 4
    if minimum is not None and macro < minimum:
        raise CoordinateLatticeLimitReached("below_minimum_aligned_width")
    exact_boundaries = tuple(x + index * fine for index in range(5))
    boundaries = tuple(_exact_float(item) for item in exact_boundaries)
    if len(boundaries) != 5:  # placates type checkers and protects the shape.
        raise CoordinateLatticeLimitReached("internal_exactness_failure")
    b0, b1, b2, b3, b4 = boundaries

    stage_paths = (
        ((b0, b0 + (b4 - b0) / 2.0, b4),),
        (
            (b0, b0 + (b2 - b0) / 2.0, b2),
            (b2, b2 + (b4 - b2) / 2.0, b4),
        ),
        tuple(
            (left, left + (right - left) / 2.0, right)
            for left, right in zip(boundaries, boundaries[1:])
        ),
    )
    expected_stage_paths = (
        ((x, x + macro / 2, x + macro),),
        (
            (x, x + macro / 4, x + macro / 2),
            (x + macro / 2, x + 3 * macro / 4, x + macro),
        ),
        tuple(
            (
                x + index * fine,
                x + (2 * index + 1) * fine / 2,
                x + (index + 1) * fine,
            )
            for index in range(4)
        ),
    )

    # Exact oracle obligations first, then the corresponding binary64 values.
    if (
        exact_boundaries[4] - exact_boundaries[0] != macro
        or any(
            exact_boundaries[index + 1] - exact_boundaries[index] != fine
            for index in range(4)
        )
        or exact_boundaries[4] > target
        or any(Fraction.from_float(item) != exact for item, exact in zip(boundaries, exact_boundaries, strict=True))
        or Fraction.from_float(b4 - b0) != macro
        or any(
            Fraction.from_float(right - left) != fine
            for left, right in zip(boundaries, boundaries[1:])
        )
        or _bits(b0) != _bits(x_float)
        or any(
            len(actual_path) != len(expected_path)
            or any(
                Fraction.from_float(actual_time) != expected_time
                for actual_step, expected_step in zip(
                    actual_path, expected_path, strict=True
                )
                for actual_time, expected_time in zip(
                    actual_step, expected_step, strict=True
                )
            )
            or any(not (left < midpoint < right) for left, midpoint, right in actual_path)
            for actual_path, expected_path in zip(
                stage_paths, expected_stage_paths, strict=True
            )
        )
    ):
        raise CoordinateLatticeLimitReached("internal_exactness_failure")
    # The endpoint object is literally shared by all three path slices.
    if not (_bits((b0, b4)[1]) == _bits((b0, b2, b4)[2]) == _bits((b0, b1, b2, b3, b4)[4])):
        raise CoordinateLatticeLimitReached("internal_exactness_failure")

    return TDG7Binary64SubdivisionPlan(
        current=b0,
        event_target=t_float,
        requested_cap=cap_float,
        quantum=q_float,
        selection_budget=_exact_float(budget),
        macro_width=b4 - b0,
        fine_width=b1 - b0,
        conservative_reduction=_exact_float(budget - macro),
        minimum_width=minimum_float,
        boundaries=(b0, b1, b2, b3, b4),
        target_limited=(macro == remaining),
    )


def cal11_event24_witness() -> TDG7Binary64SubdivisionPlan:
    """Return the exact RK4-2049 event-24 cap witness, without running a campaign."""

    return plan_forward_proto14_subdivision(
        23.0 / 16.0,
        24.0 / 16.0,
        float.fromhex("0x1.aaa90b0fb5c26p-8"),
    )
