"""Prospective TDG6 temporal-admission threshold design.

This module contains no campaign I/O and advances no numerical state.  It
encodes the exact, outcome-independent discriminator frozen after TDG5-IMP1:

* three same-grid refinement levels (one full, two half, four quarter steps),
* a conservative interval test for at least three-halves contraction order,
* an enclosure-dominated class whose entire uncertainty remains a public
  debit rather than being rounded into an order claim, and
* a nonnegative finest-pair debit with no stand-alone absolute tolerance.

The future production compositor and the independently bound TDG6 theorem are
deliberately outside this module.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from numbers import Rational
from typing import Final


TDG6_MINIMUM_OBSERVED_ORDER: Final[Fraction] = Fraction(3, 2)
TDG6_LEVEL_STEP_COUNTS: Final[tuple[int, int, int]] = (1, 2, 4)
TDG6_ORDER_SQUARED_MULTIPLIER: Final[int] = 8
TDG6_RETRY_FACTOR: Final[Fraction] = Fraction(1, 2)
TDG6_MAXIMUM_TEMPORAL_RETRIES_PER_STEP: Final[int] = 32
TDG6_MINIMUM_MACRO_STEP: Final[Fraction] = Fraction(1, 1_073_741_824)
TDG6_METHOD_STAGE_RECORDS: Final[dict[str, int]] = {
    "RK4": 5,
    "SSPRK3": 4,
}
TDG6_COMPLETE_STATE_CHANNELS: Final[tuple[str, ...]] = tuple(
    f"{block}:{field}"
    for block in ("u", "p", "q")
    for field in ("alpha", "v", "lambda", "R", "phi", "chi")
)
TDG6_PASSING_CLASSES: Final[frozenset[str]] = frozenset(
    {
        "exact_zero",
        "enclosure_dominated_debit_only",
        "resolved_order_pass",
    }
)


def _fraction(name: str, value: Rational) -> Fraction:
    if not isinstance(value, Rational) or isinstance(value, bool):
        raise TypeError(f"{name} must be rational")
    answer = Fraction(value)
    if answer < 0:
        raise ValueError(f"{name} must be nonnegative")
    return answer


@dataclass(frozen=True, slots=True)
class CertifiedMagnitudeInterval:
    """One nonnegative outward interval for a continuous path difference."""

    lower: Fraction
    upper: Fraction

    def __post_init__(self) -> None:
        lower = _fraction("lower", self.lower)
        upper = _fraction("upper", self.upper)
        if lower > upper:
            raise ValueError("magnitude interval lower bound exceeds upper bound")
        object.__setattr__(self, "lower", lower)
        object.__setattr__(self, "upper", upper)


@dataclass(frozen=True, slots=True)
class TDG6ChannelAdmission:
    """Exact prospective classification of one complete-state channel."""

    classification: str
    admission_passed: bool
    temporal_retry_permitted: bool
    minimum_observed_order: Fraction
    outer_difference: CertifiedMagnitudeInterval
    finest_difference: CertifiedMagnitudeInterval
    finest_pair_debit: Fraction
    order_threshold_resolved: bool
    order_threshold_passed: bool | None
    absolute_state_tolerance_used: bool = False
    physical_signal_used_for_normalization: bool = False

    def __post_init__(self) -> None:
        if self.classification not in TDG6_PASSING_CLASSES | {
            "resolved_order_failure",
            "order_inconclusive",
        }:
            raise ValueError("unknown TDG6 channel classification")
        if self.admission_passed is not (
            self.classification in TDG6_PASSING_CLASSES
        ):
            raise ValueError("TDG6 pass flag differs from classification")
        if self.temporal_retry_permitted is not (not self.admission_passed):
            raise ValueError("only a failed numerical admission may request retry")
        if Fraction(self.minimum_observed_order) != TDG6_MINIMUM_OBSERVED_ORDER:
            raise ValueError("TDG6 minimum observed order differs from the freeze")
        if self.finest_pair_debit != self.finest_difference.upper:
            raise ValueError("TDG6 debit must retain the full finest-pair upper bound")
        if self.absolute_state_tolerance_used or self.physical_signal_used_for_normalization:
            raise ValueError("TDG6 may not hide an absolute or outcome-scaled tolerance")
        resolved = self.classification in {
            "resolved_order_pass",
            "resolved_order_failure",
        }
        if self.order_threshold_resolved is not resolved:
            raise ValueError("TDG6 order-resolution flag differs from classification")
        expected_order_flag = (
            True
            if self.classification == "resolved_order_pass"
            else False
            if self.classification == "resolved_order_failure"
            else None
        )
        if self.order_threshold_passed is not expected_order_flag:
            raise ValueError("TDG6 order-pass flag differs from classification")


def three_halves_order_passes_squared(
    *,
    outer_lower_squared: Rational,
    finest_upper_squared: Rational,
) -> bool:
    """Return the exact ``p >= 3/2`` contraction decision.

    For nonnegative differences ``D_01`` and ``D_12``, order at least three
    halves means ``D_12 <= D_01 / 2^(3/2)``.  Squaring removes the irrational
    comparison and gives the exact rational inequality

    ``8 D_12^2 <= D_01^2``.
    """

    outer = _fraction("outer_lower_squared", outer_lower_squared)
    finest = _fraction("finest_upper_squared", finest_upper_squared)
    return TDG6_ORDER_SQUARED_MULTIPLIER * finest <= outer


def classify_tdg6_channel(
    outer_difference: CertifiedMagnitudeInterval,
    finest_difference: CertifiedMagnitudeInterval,
) -> TDG6ChannelAdmission:
    """Classify one channel without an absolute or outcome-scaled tolerance."""

    if not isinstance(outer_difference, CertifiedMagnitudeInterval):
        raise TypeError("outer_difference must be CertifiedMagnitudeInterval")
    if not isinstance(finest_difference, CertifiedMagnitudeInterval):
        raise TypeError("finest_difference must be CertifiedMagnitudeInterval")

    if outer_difference.upper == 0 and finest_difference.upper == 0:
        classification = "exact_zero"
        resolved = False
        passed: bool | None = None
    elif (
        outer_difference.lower == 0
        and finest_difference.lower == 0
        and finest_difference.upper <= outer_difference.upper
    ):
        # No convergence order is invented when both certified intervals touch
        # zero.  The complete fine uncertainty is retained as a debit, and the
        # class passes only because refinement did not enlarge that enclosure.
        classification = "enclosure_dominated_debit_only"
        resolved = False
        passed = None
    elif outer_difference.lower > 0:
        passed = three_halves_order_passes_squared(
            outer_lower_squared=outer_difference.lower**2,
            finest_upper_squared=finest_difference.upper**2,
        )
        classification = (
            "resolved_order_pass" if passed else "resolved_order_failure"
        )
        resolved = True
    else:
        classification = "order_inconclusive"
        resolved = False
        passed = None

    return TDG6ChannelAdmission(
        classification=classification,
        admission_passed=classification in TDG6_PASSING_CLASSES,
        temporal_retry_permitted=classification not in TDG6_PASSING_CLASSES,
        minimum_observed_order=TDG6_MINIMUM_OBSERVED_ORDER,
        outer_difference=outer_difference,
        finest_difference=finest_difference,
        finest_pair_debit=finest_difference.upper,
        order_threshold_resolved=resolved,
        order_threshold_passed=passed,
    )


def tdg6_path_accounting(method: str) -> dict[str, int]:
    """Return exact proposal/stage ownership for the frozen three levels."""

    if method not in TDG6_METHOD_STAGE_RECORDS:
        raise ValueError("unknown TDG6 method")
    records = TDG6_METHOD_STAGE_RECORDS[method]
    total_proposals = sum(TDG6_LEVEL_STEP_COUNTS)
    committed_proposals = TDG6_LEVEL_STEP_COUNTS[-1]
    return {
        "coarse_proposals": TDG6_LEVEL_STEP_COUNTS[0],
        "medium_proposals": TDG6_LEVEL_STEP_COUNTS[1],
        "fine_proposals": committed_proposals,
        "total_shadow_proposals": total_proposals,
        "total_shadow_stage_records": total_proposals * records,
        "committed_fine_proposals": committed_proposals,
        "committed_fine_stage_records": committed_proposals * records,
    }


def tdg6_design_preflight() -> dict[str, object]:
    """Execute exact no-trajectory controls for the prospective freeze."""

    exact_zero = classify_tdg6_channel(
        CertifiedMagnitudeInterval(Fraction(0), Fraction(0)),
        CertifiedMagnitudeInterval(Fraction(0), Fraction(0)),
    )
    floor = classify_tdg6_channel(
        CertifiedMagnitudeInterval(Fraction(0), Fraction(1, 2**40)),
        CertifiedMagnitudeInterval(Fraction(0), Fraction(1, 2**41)),
    )
    rk4 = classify_tdg6_channel(
        CertifiedMagnitudeInterval(Fraction(16), Fraction(16)),
        CertifiedMagnitudeInterval(Fraction(1), Fraction(1)),
    )
    ssprk3 = classify_tdg6_channel(
        CertifiedMagnitudeInterval(Fraction(8), Fraction(8)),
        CertifiedMagnitudeInterval(Fraction(1), Fraction(1)),
    )
    small_nonconvergent = classify_tdg6_channel(
        CertifiedMagnitudeInterval(Fraction(1, 10**12), Fraction(1, 10**12)),
        CertifiedMagnitudeInterval(Fraction(9, 10**13), Fraction(9, 10**13)),
    )
    large_convergent = classify_tdg6_channel(
        CertifiedMagnitudeInterval(Fraction(80), Fraction(80)),
        CertifiedMagnitudeInterval(Fraction(10), Fraction(10)),
    )
    inconclusive = classify_tdg6_channel(
        CertifiedMagnitudeInterval(Fraction(0), Fraction(1, 2**40)),
        CertifiedMagnitudeInterval(Fraction(1, 2**42), Fraction(1, 2**39)),
    )

    boundary_pass = three_halves_order_passes_squared(
        outer_lower_squared=Fraction(8), finest_upper_squared=Fraction(1)
    )
    boundary_miss = three_halves_order_passes_squared(
        outer_lower_squared=Fraction(8) - Fraction(1, 2**52),
        finest_upper_squared=Fraction(1),
    )
    controls_pass = (
        exact_zero.classification == "exact_zero"
        and exact_zero.finest_pair_debit == 0
        and floor.classification == "enclosure_dominated_debit_only"
        and floor.finest_pair_debit == Fraction(1, 2**41)
        and rk4.classification == "resolved_order_pass"
        and ssprk3.classification == "resolved_order_pass"
        and small_nonconvergent.classification == "resolved_order_failure"
        and not small_nonconvergent.admission_passed
        and large_convergent.classification == "resolved_order_pass"
        and large_convergent.finest_pair_debit == 10
        and inconclusive.classification == "order_inconclusive"
        and inconclusive.temporal_retry_permitted
        and boundary_pass
        and not boundary_miss
        and tdg6_path_accounting("RK4")["committed_fine_stage_records"] == 20
        and tdg6_path_accounting("SSPRK3")["committed_fine_stage_records"] == 16
    )
    if not controls_pass:
        raise AssertionError("TDG6 exact design preflight failed")

    def admission(value: TDG6ChannelAdmission) -> dict[str, object]:
        return {
            "classification": value.classification,
            "admission_passed": value.admission_passed,
            "temporal_retry_permitted": value.temporal_retry_permitted,
            "finest_pair_debit": str(value.finest_pair_debit),
            "order_threshold_resolved": value.order_threshold_resolved,
            "order_threshold_passed": value.order_threshold_passed,
        }

    return {
        "minimum_observed_order": str(TDG6_MINIMUM_OBSERVED_ORDER),
        "exact_squared_order_inequality": "8*D12_upper_squared<=D01_lower_squared",
        "complete_state_channel_count": len(TDG6_COMPLETE_STATE_CHANNELS),
        "path_accounting": {
            method: tdg6_path_accounting(method)
            for method in TDG6_METHOD_STAGE_RECORDS
        },
        "controls": {
            "exact_zero": admission(exact_zero),
            "enclosure_dominated": admission(floor),
            "RK4_formal_contraction": admission(rk4),
            "SSPRK3_formal_contraction": admission(ssprk3),
            "small_nonconvergent": admission(small_nonconvergent),
            "large_convergent": admission(large_convergent),
            "order_inconclusive": admission(inconclusive),
            "exact_boundary_passed": boundary_pass,
            "one_rational_unit_below_boundary_rejected": not boundary_miss,
        },
        "no_absolute_state_tolerance": True,
        "physical_signal_used_for_normalization": False,
        "controls_passed": controls_pass,
    }


__all__ = [
    "CertifiedMagnitudeInterval",
    "TDG6ChannelAdmission",
    "TDG6_COMPLETE_STATE_CHANNELS",
    "TDG6_LEVEL_STEP_COUNTS",
    "TDG6_MAXIMUM_TEMPORAL_RETRIES_PER_STEP",
    "TDG6_MINIMUM_MACRO_STEP",
    "TDG6_MINIMUM_OBSERVED_ORDER",
    "TDG6_ORDER_SQUARED_MULTIPLIER",
    "TDG6_RETRY_FACTOR",
    "classify_tdg6_channel",
    "tdg6_design_preflight",
    "tdg6_path_accounting",
    "three_halves_order_passes_squared",
]
