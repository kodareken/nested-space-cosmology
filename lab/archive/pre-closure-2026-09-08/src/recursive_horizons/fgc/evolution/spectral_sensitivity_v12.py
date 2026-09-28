"""PROTO12's versioned pairwise diagnostic-conditioning classifier.

The pre-PROTO12 :mod:`spectral_sensitivity` module is immutable historical
implementation evidence for CAL6, HLT8, and HLT9.  This module extends that
API without changing its bytes.  Direct adjacent-tail contraction remains the
primary route, every raw ratio remains public, and diagnostic saturation is
neither a continuum-error bound nor physical-resolution evidence.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

from .proto4_admission import (
    PROTO4_SPECTRAL_FIELD_ORDER,
    SpectralPowerBudget,
    SpectralThresholds,
)
from .spectral_sensitivity import (
    SpectralTailSensitivity,
    ThreeGridProfileConvergence,
)


@dataclass(frozen=True, slots=True)
class PairwiseResolvedOrSaturatedAdmission:
    """Direct-or-conditioned classification on both adjacent grid pairs."""

    point_counts: tuple[int, int, int]
    field_power_tail_ratios: Mapping[str, tuple[float | None, float | None]]
    derivative_power_tail_ratios: Mapping[
        str, tuple[float | None, float | None]
    ]
    field_pair_classification: Mapping[str, tuple[str, str]]
    derivative_pair_classification: Mapping[str, tuple[str, str]]
    pair_saturation_available: Mapping[str, tuple[bool, bool]]
    medium_and_fine_individual_budgets_passed: bool
    every_profile_contracted: bool
    saturation_used: bool
    every_tail_resolved_or_saturated: bool
    admission_passed: bool


def pairwise_resolved_or_saturated_admission(
    point_counts: Sequence[int],
    budgets: Sequence[Mapping[str, SpectralPowerBudget]],
    sensitivities: Sequence[Mapping[str, SpectralTailSensitivity]],
    profile_convergence: ThreeGridProfileConvergence,
    *,
    thresholds: SpectralThresholds | None = None,
) -> PairwiseResolvedOrSaturatedAdmission:
    """Classify both adjacent pairs without hiding any failed raw ratio.

    Direct ``ratio < maximum_nested_tail_ratio`` remains the primary route.
    A failed ratio may be labelled ``diagnostically_saturated`` only when both
    grids in that pair pass their unchanged absolute budgets, both map-scale
    witnesses pass, and the complete three-grid profile/round-trip contraction
    passes.  No measured map residual is used as a continuum-error estimate.
    """

    limits = SpectralThresholds() if thresholds is None else thresholds
    if not isinstance(limits, SpectralThresholds):
        raise TypeError("thresholds must be SpectralThresholds")
    counts = tuple(point_counts)
    records = tuple(budgets)
    witnesses = tuple(sensitivities)
    names = tuple(PROTO4_SPECTRAL_FIELD_ORDER)
    if (
        len(counts) != 3
        or len(records) != 3
        or len(witnesses) != 3
        or tuple(profile_convergence.point_counts) != counts
        or any(tuple(record) != names for record in records)
        or any(tuple(record) != names for record in witnesses)
    ):
        raise ValueError(
            "pairwise resolved-or-saturated inputs differ from the three-grid contract"
        )

    def ratio(numerator: float, denominator: float) -> float | None:
        if denominator > 0.0:
            return numerator / denominator
        if numerator == 0.0:
            return 0.0
        return None

    def classify(value: float | None, saturation_available: bool) -> str:
        if value is not None and value < limits.maximum_nested_tail_ratio:
            return "directly_resolved"
        if saturation_available:
            return "diagnostically_saturated"
        return "failed"

    field_ratios: dict[str, tuple[float | None, float | None]] = {}
    derivative_ratios: dict[str, tuple[float | None, float | None]] = {}
    field_classification: dict[str, tuple[str, str]] = {}
    derivative_classification: dict[str, tuple[str, str]] = {}
    saturation_by_pair: dict[str, tuple[bool, bool]] = {}
    saturation_used = False
    for name in names:
        field_values = tuple(
            record[name].top_band_field_power_fraction for record in records
        )
        derivative_values = tuple(
            record[name].top_band_derivative_weighted_power_fraction
            for record in records
        )
        f_ratios = (
            ratio(field_values[1], field_values[0]),
            ratio(field_values[2], field_values[1]),
        )
        d_ratios = (
            ratio(derivative_values[1], derivative_values[0]),
            ratio(derivative_values[2], derivative_values[1]),
        )
        available = tuple(
            records[pair][name].individual_admission_passed
            and records[pair + 1][name].individual_admission_passed
            and witnesses[pair][name].diagnostically_saturated
            and witnesses[pair + 1][name].diagnostically_saturated
            and profile_convergence.complete_profile_contraction_by_field[name]
            for pair in range(2)
        )
        f_class = tuple(
            classify(f_ratios[pair], available[pair]) for pair in range(2)
        )
        d_class = tuple(
            classify(d_ratios[pair], available[pair]) for pair in range(2)
        )
        field_ratios[name] = f_ratios
        derivative_ratios[name] = d_ratios
        field_classification[name] = (f_class[0], f_class[1])
        derivative_classification[name] = (d_class[0], d_class[1])
        saturation_by_pair[name] = (available[0], available[1])
        saturation_used = saturation_used or any(
            value == "diagnostically_saturated" for value in (*f_class, *d_class)
        )

    every_tail = all(
        value != "failed"
        for classifications in (
            *field_classification.values(),
            *derivative_classification.values(),
        )
        for value in classifications
    )
    medium_and_fine_budgets = all(
        records[index][name].individual_admission_passed
        for index in (1, 2)
        for name in names
    )
    admitted = (
        medium_and_fine_budgets
        and profile_convergence.every_field_contracted
        and every_tail
    )
    return PairwiseResolvedOrSaturatedAdmission(
        point_counts=(int(counts[0]), int(counts[1]), int(counts[2])),
        field_power_tail_ratios=field_ratios,
        derivative_power_tail_ratios=derivative_ratios,
        field_pair_classification=field_classification,
        derivative_pair_classification=derivative_classification,
        pair_saturation_available=saturation_by_pair,
        medium_and_fine_individual_budgets_passed=medium_and_fine_budgets,
        every_profile_contracted=profile_convergence.every_field_contracted,
        saturation_used=saturation_used,
        every_tail_resolved_or_saturated=every_tail,
        admission_passed=admitted,
    )


__all__ = [
    "PairwiseResolvedOrSaturatedAdmission",
    "pairwise_resolved_or_saturated_admission",
]
