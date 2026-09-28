"""Outcome-blind conditioning checks for the proper-grid spectrum.

The PROTO4/PROTO8 spatial diagnostic maps fields from the native radial grid
to a uniform proper-distance grid before applying an FFT.  Its round-trip
residual is deliberately *not* a continuum interpolation-error bound.  It is,
however, a directly measured non-idempotence scale of that diagnostic map.

This module keeps those two statements separate.  A nested tail ratio retains
its original direct ``< 1/4`` route.  A failed ratio may only be labelled
``diagnostically_saturated`` when all of the following are independently true:

* the unchanged absolute budgets pass on the medium and fine grids;
* the coarse-to-medium ratio still passes directly;
* erasing the fine-grid top band changes the diagnostic input by no more than
  the measured round-trip non-idempotence on both medium and fine grids;
* the round-trip residual contracts strictly on both refinement steps; and
* the complete resampled profile contracts strictly in both RMS and infinity
  norm from the first grid pair to the second.

The result is a conditioning classification, not permission to call the
round-trip residual a physical uncertainty or a continuum error estimate.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite, pi, sqrt
from numbers import Integral, Real
from typing import Mapping, Sequence

import numpy as np

from .proto4_admission import (
    PROTO4_SPECTRAL_FIELD_ORDER,
    ProperRadialProfile,
    SpectralPowerBudget,
    SpectralThresholds,
)


def _finite(name: str, value: Real, *, nonnegative: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{name} must be real")
    answer = float(value)
    if not isfinite(answer):
        raise ValueError(f"{name} must be finite")
    if nonnegative and answer < 0.0:
        raise ValueError(f"{name} must be nonnegative")
    return answer


def _array(name: str, value: object, *, ndim: int) -> np.ndarray:
    answer = np.asarray(value, dtype=np.float64)
    if answer.ndim != ndim or not np.all(np.isfinite(answer)):
        raise ValueError(f"{name} must be a finite {ndim}-dimensional array")
    return answer


@dataclass(frozen=True, slots=True)
class SpectralTailSensitivity:
    """One field's top-band response relative to map non-idempotence."""

    sample_count: int
    top_band_first_bin: int
    nyquist_angular_frequency: float
    round_trip_interpolation_infinity: float
    top_band_field_rms_amplitude: float
    top_band_derivative_rms_amplitude: float
    top_band_erasure_perturbation_infinity: float
    erasure_over_round_trip: float | None
    field_tail_within_round_trip_scale: bool
    derivative_tail_within_nyquist_round_trip_scale: bool
    erasure_within_round_trip_scale: bool
    diagnostically_saturated: bool

    def __post_init__(self) -> None:
        for name in ("sample_count", "top_band_first_bin"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, Integral) or value < 0:
                raise ValueError(f"{name} must be a nonnegative integer")
        if self.sample_count < 16 or self.top_band_first_bin > self.sample_count // 2:
            raise ValueError("spectral sensitivity indices are inconsistent")
        for name in (
            "nyquist_angular_frequency",
            "round_trip_interpolation_infinity",
            "top_band_field_rms_amplitude",
            "top_band_derivative_rms_amplitude",
            "top_band_erasure_perturbation_infinity",
        ):
            object.__setattr__(
                self,
                name,
                _finite(name, getattr(self, name), nonnegative=True),
            )
        ratio = self.erasure_over_round_trip
        if ratio is not None:
            object.__setattr__(
                self,
                "erasure_over_round_trip",
                _finite("erasure_over_round_trip", ratio, nonnegative=True),
            )
        for name in (
            "field_tail_within_round_trip_scale",
            "derivative_tail_within_nyquist_round_trip_scale",
            "erasure_within_round_trip_scale",
            "diagnostically_saturated",
        ):
            if not isinstance(getattr(self, name), bool):
                raise TypeError(f"{name} must be bool")


def spectral_tail_sensitivity(
    samples: object,
    proper_coordinates: object,
    *,
    window: object,
    round_trip_interpolation_infinity: Real,
    budget: SpectralPowerBudget,
) -> SpectralTailSensitivity:
    """Measure whether a top band lies below the map's observed sensitivity.

    The top-band erasure is the orthogonal real-FFT projection obtained by
    setting the diagnostic's top-band coefficients to zero.  Its infinity
    norm is an explicit pass/fail-flipping witness in diagnostic-input space.
    No claim is made that the round-trip value bounds the continuum field.
    """

    data = _array("samples", samples, ndim=1)
    proper = _array("proper_coordinates", proper_coordinates, ndim=1)
    weights = _array("window", window, ndim=1)
    if data.shape != proper.shape or data.shape != weights.shape:
        raise ValueError("spectral sensitivity arrays must share one shape")
    if data.size != budget.sample_count:
        raise ValueError("spectral sensitivity budget sample count differs")
    if np.any(weights < 0.0) or np.any(weights > 1.0):
        raise ValueError("spectral sensitivity window must lie in [0,1]")
    spacing = np.diff(proper)
    if spacing.size == 0 or np.any(spacing <= 0.0):
        raise ValueError("proper coordinates must increase strictly")
    tolerance = 4096.0 * np.finfo(np.float64).eps * max(1.0, abs(spacing[0]))
    if not np.allclose(spacing, spacing[0], rtol=0.0, atol=tolerance):
        raise ValueError("spectral sensitivity requires a uniform proper grid")
    round_trip = _finite(
        "round_trip_interpolation_infinity",
        round_trip_interpolation_infinity,
        nonnegative=True,
    )

    transformed = np.fft.rfft(data * weights)
    filtered = transformed.copy()
    filtered[budget.top_band_first_bin :] = 0.0
    erased = np.fft.irfft(filtered, n=data.size)
    erasure = float(np.max(np.abs(erased - data * weights), initial=0.0))

    top_field_power = (
        budget.total_field_power * budget.top_band_field_power_fraction
    )
    top_derivative_power = (
        budget.total_derivative_weighted_power
        * budget.top_band_derivative_weighted_power_fraction
    )
    field_rms = sqrt(max(0.0, top_field_power))
    derivative_rms = sqrt(max(0.0, top_derivative_power))
    nyquist = pi / float(spacing[0])
    if round_trip > 0.0:
        ratio: float | None = erasure / round_trip
        field_within = field_rms <= round_trip
        derivative_within = derivative_rms <= nyquist * round_trip
        erasure_within = erasure <= round_trip
        saturated = field_within and derivative_within and erasure_within
    else:
        ratio = 0.0 if erasure == 0.0 else None
        field_within = field_rms == 0.0
        derivative_within = derivative_rms == 0.0
        erasure_within = erasure == 0.0
        saturated = False
    return SpectralTailSensitivity(
        sample_count=data.size,
        top_band_first_bin=budget.top_band_first_bin,
        nyquist_angular_frequency=nyquist,
        round_trip_interpolation_infinity=round_trip,
        top_band_field_rms_amplitude=field_rms,
        top_band_derivative_rms_amplitude=derivative_rms,
        top_band_erasure_perturbation_infinity=erasure,
        erasure_over_round_trip=ratio,
        field_tail_within_round_trip_scale=field_within,
        derivative_tail_within_nyquist_round_trip_scale=derivative_within,
        erasure_within_round_trip_scale=erasure_within,
        diagnostically_saturated=saturated,
    )


@dataclass(frozen=True, slots=True)
class AdjacentProfileDifference:
    coarse_point_count: int
    fine_point_count: int
    comparison_sample_count: int
    common_proper_radius_maximum: float
    infinity_by_field: Mapping[str, float]
    rms_by_field: Mapping[str, float]


@dataclass(frozen=True, slots=True)
class ThreeGridProfileConvergence:
    point_counts: tuple[int, int, int]
    adjacent_differences: tuple[AdjacentProfileDifference, AdjacentProfileDifference]
    exact_zero_by_field: Mapping[str, bool]
    infinity_contracted_by_field: Mapping[str, bool]
    rms_contracted_by_field: Mapping[str, bool]
    round_trip_contracted_by_field: Mapping[str, bool]
    complete_profile_contraction_by_field: Mapping[str, bool]
    every_field_contracted: bool


def _adjacent_profile_difference(
    coarse_count: int,
    fine_count: int,
    coarse: ProperRadialProfile,
    fine: ProperRadialProfile,
) -> AdjacentProfileDifference:
    common_maximum = min(
        float(coarse.proper_coordinates[-1]),
        float(fine.proper_coordinates[-1]),
    )
    sample_count = min(coarse.proper_coordinates.size, fine.proper_coordinates.size)
    common = np.linspace(0.0, common_maximum, sample_count, dtype=np.float64)
    coarse_values = np.column_stack(
        [
            np.interp(common, coarse.proper_coordinates, coarse.deviations[:, index])
            for index in range(coarse.deviations.shape[1])
        ]
    )
    fine_values = np.column_stack(
        [
            np.interp(common, fine.proper_coordinates, fine.deviations[:, index])
            for index in range(fine.deviations.shape[1])
        ]
    )
    difference = coarse_values - fine_values
    infinity = np.max(np.abs(difference), axis=0)
    rms = np.sqrt(np.mean(difference * difference, axis=0))
    return AdjacentProfileDifference(
        coarse_point_count=coarse_count,
        fine_point_count=fine_count,
        comparison_sample_count=sample_count,
        common_proper_radius_maximum=common_maximum,
        infinity_by_field={
            name: float(infinity[index])
            for index, name in enumerate(PROTO4_SPECTRAL_FIELD_ORDER)
        },
        rms_by_field={
            name: float(rms[index])
            for index, name in enumerate(PROTO4_SPECTRAL_FIELD_ORDER)
        },
    )


def three_grid_profile_convergence(
    point_counts: Sequence[int],
    profiles: Sequence[ProperRadialProfile],
) -> ThreeGridProfileConvergence:
    """Require strict whole-profile and map-residual contraction on three grids."""

    counts = tuple(point_counts)
    records = tuple(profiles)
    if (
        len(counts) != 3
        or len(records) != 3
        or any(isinstance(value, bool) or not isinstance(value, Integral) for value in counts)
        or not counts[0] < counts[1] < counts[2]
    ):
        raise ValueError("profile convergence requires three increasing point counts")
    if any(
        profile.deviations.shape[1] != len(PROTO4_SPECTRAL_FIELD_ORDER)
        for profile in records
    ):
        raise ValueError("profile convergence field order differs")
    adjacent = (
        _adjacent_profile_difference(counts[0], counts[1], records[0], records[1]),
        _adjacent_profile_difference(counts[1], counts[2], records[1], records[2]),
    )
    exact_zero: dict[str, bool] = {}
    infinity_contracted: dict[str, bool] = {}
    rms_contracted: dict[str, bool] = {}
    round_trip_contracted: dict[str, bool] = {}
    complete: dict[str, bool] = {}
    for index, name in enumerate(PROTO4_SPECTRAL_FIELD_ORDER):
        infinity_values = tuple(item.infinity_by_field[name] for item in adjacent)
        rms_values = tuple(item.rms_by_field[name] for item in adjacent)
        round_trip_values = tuple(
            float(profile.round_trip_interpolation_infinity[index])
            for profile in records
        )
        zero = (
            infinity_values == (0.0, 0.0)
            and rms_values == (0.0, 0.0)
            and round_trip_values == (0.0, 0.0, 0.0)
        )
        inf_pass = zero or (
            infinity_values[0] > 0.0
            and 0.0 <= infinity_values[1] < infinity_values[0]
        )
        rms_pass = zero or (
            rms_values[0] > 0.0 and 0.0 <= rms_values[1] < rms_values[0]
        )
        round_trip_pass = zero or (
            round_trip_values[0] > round_trip_values[1] > round_trip_values[2] >= 0.0
        )
        exact_zero[name] = zero
        infinity_contracted[name] = inf_pass
        rms_contracted[name] = rms_pass
        round_trip_contracted[name] = round_trip_pass
        complete[name] = inf_pass and rms_pass and round_trip_pass
    return ThreeGridProfileConvergence(
        point_counts=(int(counts[0]), int(counts[1]), int(counts[2])),
        adjacent_differences=adjacent,
        exact_zero_by_field=exact_zero,
        infinity_contracted_by_field=infinity_contracted,
        rms_contracted_by_field=rms_contracted,
        round_trip_contracted_by_field=round_trip_contracted,
        complete_profile_contraction_by_field=complete,
        every_field_contracted=all(complete.values()),
    )


@dataclass(frozen=True, slots=True)
class ResolvedOrSaturatedAdmission:
    point_counts: tuple[int, int, int]
    field_power_tail_ratios: Mapping[str, tuple[float | None, float | None]]
    derivative_power_tail_ratios: Mapping[str, tuple[float | None, float | None]]
    coarse_to_medium_direct_passed_by_field: Mapping[str, bool]
    coarse_to_medium_direct_passed_by_derivative: Mapping[str, bool]
    finest_pair_field_classification: Mapping[str, str]
    finest_pair_derivative_classification: Mapping[str, str]
    finest_pair_individual_budgets_passed: bool
    every_profile_contracted: bool
    saturation_used: bool
    every_tail_resolved_or_saturated: bool
    admission_passed: bool


def resolved_or_saturated_admission(
    point_counts: Sequence[int],
    budgets: Sequence[Mapping[str, SpectralPowerBudget]],
    sensitivities: Sequence[Mapping[str, SpectralTailSensitivity]],
    profile_convergence: ThreeGridProfileConvergence,
    *,
    thresholds: SpectralThresholds | None = None,
) -> ResolvedOrSaturatedAdmission:
    """Apply the direct coarse ratio and narrow finest saturation contract."""

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
        raise ValueError("resolved-or-saturated inputs differ from the three-grid contract")

    def ratio(numerator: float, denominator: float) -> float | None:
        if denominator > 0.0:
            return numerator / denominator
        if numerator == 0.0:
            return 0.0
        return None

    field_ratios: dict[str, tuple[float | None, float | None]] = {}
    derivative_ratios: dict[str, tuple[float | None, float | None]] = {}
    coarse_field: dict[str, bool] = {}
    coarse_derivative: dict[str, bool] = {}
    fine_field: dict[str, str] = {}
    fine_derivative: dict[str, str] = {}
    fine_budgets = all(
        records[index][name].individual_admission_passed
        for index in (1, 2)
        for name in names
    )
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
        field_ratios[name] = f_ratios
        derivative_ratios[name] = d_ratios
        coarse_field[name] = (
            f_ratios[0] is not None
            and f_ratios[0] < limits.maximum_nested_tail_ratio
        )
        coarse_derivative[name] = (
            d_ratios[0] is not None
            and d_ratios[0] < limits.maximum_nested_tail_ratio
        )
        saturation_available = (
            records[1][name].individual_admission_passed
            and records[2][name].individual_admission_passed
            and witnesses[1][name].diagnostically_saturated
            and witnesses[2][name].diagnostically_saturated
            and profile_convergence.complete_profile_contraction_by_field[name]
        )

        if f_ratios[1] is not None and f_ratios[1] < limits.maximum_nested_tail_ratio:
            fine_field[name] = "directly_resolved"
        elif saturation_available:
            fine_field[name] = "diagnostically_saturated"
            saturation_used = True
        else:
            fine_field[name] = "failed"
        if d_ratios[1] is not None and d_ratios[1] < limits.maximum_nested_tail_ratio:
            fine_derivative[name] = "directly_resolved"
        elif saturation_available:
            fine_derivative[name] = "diagnostically_saturated"
            saturation_used = True
        else:
            fine_derivative[name] = "failed"

    every_tail = (
        all(coarse_field.values())
        and all(coarse_derivative.values())
        and all(value != "failed" for value in fine_field.values())
        and all(value != "failed" for value in fine_derivative.values())
    )
    admitted = (
        fine_budgets
        and profile_convergence.every_field_contracted
        and every_tail
    )
    return ResolvedOrSaturatedAdmission(
        point_counts=(int(counts[0]), int(counts[1]), int(counts[2])),
        field_power_tail_ratios=field_ratios,
        derivative_power_tail_ratios=derivative_ratios,
        coarse_to_medium_direct_passed_by_field=coarse_field,
        coarse_to_medium_direct_passed_by_derivative=coarse_derivative,
        finest_pair_field_classification=fine_field,
        finest_pair_derivative_classification=fine_derivative,
        finest_pair_individual_budgets_passed=fine_budgets,
        every_profile_contracted=profile_convergence.every_field_contracted,
        saturation_used=saturation_used,
        every_tail_resolved_or_saturated=every_tail,
        admission_passed=admitted,
    )


__all__ = [
    "AdjacentProfileDifference",
    "ResolvedOrSaturatedAdmission",
    "SpectralTailSensitivity",
    "ThreeGridProfileConvergence",
    "resolved_or_saturated_admission",
    "spectral_tail_sensitivity",
    "three_grid_profile_convergence",
]
