"""Prospective absolute-tail discriminator for the TDG2 diagnosis.

TDG1 proved that normalizing every temporal tail by its own total power
destroys amplitude-convergence information.  This module does not modify that
historical result and does not define a replacement admission gate.  It
provides an outcome-neutral discriminator for the next read-only analysis:

* align each three-resolution history on its common proper-time interval;
* measure *absolute*, compact-windowed top-band field and derivative power;
* independently recompute the five top-band DFT coefficients with Decimal;
* expose binary64/reference disagreement and a conditional interpolation
  perturbation enclosure; and
* distinguish power contracting toward zero, nonzero resolution-independent
  or convergent power, enclosure/floor domination, and unresolved behavior.

The interpolation interval is conditional on the public round-trip debit being
used as a pointwise perturbation radius.  It is deliberately not described as
a theorem about the unknown continuum history.  No classification produced
here can retroactively pass PROTO13, authorize PROTO14, or open a candidate
trajectory.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from decimal import Decimal, localcontext
from functools import lru_cache
from math import ceil, hypot, isfinite, log, pi
from numbers import Real
from typing import Any, Mapping, Sequence

import numpy as np

from .health_monitor import (
    CONSTRAINT_NORMALIZATION_FLOOR,
    compact_vacuum_buffer_window,
)
from .proto4_admission import PROTO4_SPECTRAL_FIELD_ORDER


TDG2_SAMPLE_COUNT = 64
TDG2_TRACER_COUNT = 48
TDG2_FIELD_NAMES = PROTO4_SPECTRAL_FIELD_ORDER
TDG2_TAPER_FRACTION = 1.0 / 8.0
TDG2_UNRESOLVED_TOP_FRACTION = 1.0 / 8.0
TDG2_ABSOLUTE_AMPLITUDE_FLOOR = CONSTRAINT_NORMALIZATION_FLOOR
TDG2_MINIMUM_HISTORY_ORDER = 3.0 / 2.0
TDG2_MINIMUM_ZERO_POWER_ORDER = 3.0
TDG2_DECIMAL_PRECISION_LOW = 72
TDG2_DECIMAL_PRECISION_HIGH = 96
TDG2_DECIMAL_PI = (
    "3.141592653589793238462643383279502884197169399375105820974944592307"
    "8164062862089986280348253421170679821480865132823066470938446"
)
TDG2_METHOD_POINT_COUNTS: Mapping[str, tuple[int, int, int]] = {
    "RK4": (2049, 4097, 8193),
    "SSPRK3": (4097, 8193, 16385),
}
TDG2_MEASURE_NAMES = (
    "field_tail_power",
    "derivative_tail_power",
)
TDG2_CLASS_NAMES = (
    "contracts_to_zero_at_required_power_order",
    "nonzero_resolution_independent_within_enclosure",
    "converges_to_nonzero_resolved_power",
    "instrument_floor_or_interpolation_enclosure_dominated",
    "unresolved_absolute_tail_behavior",
)
TDG2_SYNTHETIC_CONTROL_NAMES = (
    "exact_zero",
    "subfloor_nonzero",
    "power_contracts_to_zero",
    "nonzero_resolution_independent",
    "nonzero_convergent",
    "unresolved_nonmonotone",
    "interpolation_enclosure_dominated",
)


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _finite(name: str, value: object, *, positive: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a real scalar")
    answer = float(value)
    if not isfinite(answer):
        raise ValueError(f"{name} must be finite")
    if positive and answer <= 0.0:
        raise ValueError(f"{name} must be positive")
    return answer


def _array(name: str, value: object, *, ndim: int) -> np.ndarray:
    answer = np.asarray(value, dtype=np.float64)
    if answer.ndim != ndim:
        raise ValueError(f"{name} must have rank {ndim}")
    if not np.all(np.isfinite(answer)):
        raise ValueError(f"{name} must contain only finite values")
    return answer


def _outward_interval(value: float, debit: float) -> tuple[float, float]:
    centre = _finite("interval value", value)
    radius = _finite("interval debit", debit)
    if centre < 0.0 or radius < 0.0:
        raise ValueError("power intervals require nonnegative values and debits")
    lower = max(0.0, centre - radius)
    upper = centre + radius
    return (
        float(np.nextafter(lower, -np.inf)) if lower > 0.0 else 0.0,
        float(np.nextafter(upper, np.inf)),
    )


@dataclass(frozen=True, slots=True)
class AbsolutePowerInterval:
    """One absolute power measurement and its public conditional enclosure."""

    observed: float
    binary64_fft_value: float
    arithmetic_debit: float
    interpolation_debit: float
    total_debit: float
    lower: float
    upper: float

    def __post_init__(self) -> None:
        for name in (
            "observed",
            "binary64_fft_value",
            "arithmetic_debit",
            "interpolation_debit",
            "total_debit",
            "lower",
            "upper",
        ):
            value = _finite(name, getattr(self, name))
            if value < 0.0:
                raise ValueError(f"{name} must be nonnegative")
            object.__setattr__(self, name, value)
        if self.total_debit < max(self.arithmetic_debit, self.interpolation_debit):
            raise ValueError("total power debit omits an owning contribution")
        if not self.lower <= self.observed <= self.upper:
            raise ValueError("power observation lies outside its interval")


@dataclass(frozen=True, slots=True)
class AbsoluteTailObservation:
    """Absolute top-band evidence for one aligned 64-sample signal."""

    sample_count: int
    top_band_first_bin: int
    nyquist_bin: int
    proper_time_start: float
    proper_time_stop: float
    proper_time_spacing: float
    round_trip_interpolation_infinity: float
    window_L1: float
    peak_top_band_coefficient_amplitude: float
    binary64_peak_top_band_coefficient_amplitude: float
    coefficient_arithmetic_debit: float
    coefficient_interpolation_debit: float
    declared_absolute_amplitude_floor: float
    top_band_is_below_declared_amplitude_floor: bool
    field_tail_power: AbsolutePowerInterval
    derivative_tail_power: AbsolutePowerInterval

    def __post_init__(self) -> None:
        if self.sample_count != TDG2_SAMPLE_COUNT:
            raise ValueError("TDG2 observation sample count differs")
        if self.nyquist_bin != self.sample_count // 2:
            raise ValueError("TDG2 Nyquist bin differs")
        expected_first = int(
            ceil((1.0 - TDG2_UNRESOLVED_TOP_FRACTION) * self.nyquist_bin)
        )
        if self.top_band_first_bin != expected_first:
            raise ValueError("TDG2 top-band boundary differs")
        for name in (
            "proper_time_start",
            "proper_time_stop",
            "proper_time_spacing",
            "round_trip_interpolation_infinity",
            "window_L1",
            "peak_top_band_coefficient_amplitude",
            "binary64_peak_top_band_coefficient_amplitude",
            "coefficient_arithmetic_debit",
            "coefficient_interpolation_debit",
            "declared_absolute_amplitude_floor",
        ):
            value = _finite(name, getattr(self, name))
            object.__setattr__(self, name, value)
        if (
            self.proper_time_stop <= self.proper_time_start
            or self.proper_time_spacing <= 0.0
            or min(
                self.round_trip_interpolation_infinity,
                self.window_L1,
                self.peak_top_band_coefficient_amplitude,
                self.binary64_peak_top_band_coefficient_amplitude,
                self.coefficient_arithmetic_debit,
                self.coefficient_interpolation_debit,
                self.declared_absolute_amplitude_floor,
            )
            < 0.0
        ):
            raise ValueError("TDG2 observation geometry or debit differs")
        if not isinstance(self.top_band_is_below_declared_amplitude_floor, bool):
            raise TypeError("TDG2 floor decision must be bool")


@dataclass(frozen=True, slots=True)
class AbsolutePowerClassification:
    """One outcome-neutral classification of a three-resolution power ladder."""

    measure: str
    classification: str
    point_counts: tuple[int, int, int]
    refinement_ratio: float
    observed_powers: tuple[float, float, float]
    lower_bounds: tuple[float, float, float]
    upper_bounds: tuple[float, float, float]
    zero_power_order_lower_bounds: tuple[float | None, float | None]
    difference_order_lower_bound: float | None
    conservative_continuum_lower_bound: float | None
    every_resolution_below_declared_amplitude_floor: bool
    enclosure_controls_classification: bool
    replacement_admission_authorized: bool = False

    def __post_init__(self) -> None:
        if self.measure not in TDG2_MEASURE_NAMES:
            raise ValueError("unknown TDG2 power measure")
        if self.classification not in TDG2_CLASS_NAMES:
            raise ValueError("unknown TDG2 classification")
        if len(self.point_counts) != 3 or len(set(self.point_counts)) != 3:
            raise ValueError("TDG2 classification requires three resolutions")
        if self.refinement_ratio <= 1.0:
            raise ValueError("TDG2 refinement ratio must exceed one")
        if any(len(record) != 3 for record in (self.observed_powers, self.lower_bounds, self.upper_bounds)):
            raise ValueError("TDG2 power vectors must contain three values")
        if self.replacement_admission_authorized:
            raise ValueError("TDG2 may not authorize a replacement admission")


def _compact_window(coordinate: np.ndarray) -> np.ndarray:
    width = TDG2_TAPER_FRACTION * float(coordinate[-1] - coordinate[0])
    return compact_vacuum_buffer_window(
        coordinate,
        support_minimum=float(coordinate[0]),
        support_maximum=float(coordinate[-1]),
        taper_width=width,
    )


def _decimal_sin_cos(angle: Decimal, *, precision: int) -> tuple[Decimal, Decimal]:
    """Evaluate sine and cosine deterministically after reduction to [-pi, pi]."""

    with localcontext() as context:
        context.prec = precision + 20
        decimal_pi = Decimal(TDG2_DECIMAL_PI)
        two_pi = Decimal(2) * decimal_pi
        reduced = angle % two_pi
        if reduced > decimal_pi:
            reduced -= two_pi
        x2 = reduced * reduced
        tolerance = Decimal(10) ** Decimal(-(precision + 8))

        sine_term = reduced
        sine_sum = reduced
        cosine_term = Decimal(1)
        cosine_sum = Decimal(1)
        index = 1
        while True:
            sine_term *= -x2 / Decimal((2 * index) * (2 * index + 1))
            cosine_term *= -x2 / Decimal((2 * index - 1) * (2 * index))
            sine_sum += sine_term
            cosine_sum += cosine_term
            if abs(sine_term) < tolerance and abs(cosine_term) < tolerance:
                break
            index += 1
            if index > 256:
                raise ArithmeticError("TDG2 Decimal trigonometric series did not converge")
        context.prec = precision
        return +sine_sum, +cosine_sum


@lru_cache(maxsize=None)
def _decimal_top_band_twiddles(
    precision: int,
) -> tuple[tuple[tuple[Decimal, Decimal], ...], ...]:
    if precision not in (TDG2_DECIMAL_PRECISION_LOW, TDG2_DECIMAL_PRECISION_HIGH):
        raise ValueError("TDG2 Decimal precision differs")
    nyquist = TDG2_SAMPLE_COUNT // 2
    first_top = int(ceil((1.0 - TDG2_UNRESOLVED_TOP_FRACTION) * nyquist))
    with localcontext() as context:
        context.prec = precision + 20
        two_pi = Decimal(2) * Decimal(TDG2_DECIMAL_PI)
        records: list[tuple[tuple[Decimal, Decimal], ...]] = []
        for frequency_bin in range(first_top, nyquist + 1):
            row: list[tuple[Decimal, Decimal]] = []
            for sample_index in range(TDG2_SAMPLE_COUNT):
                angle = (
                    two_pi
                    * Decimal(frequency_bin)
                    * Decimal(sample_index)
                    / Decimal(TDG2_SAMPLE_COUNT)
                )
                sine, cosine = _decimal_sin_cos(angle, precision=precision)
                row.append((cosine, -sine))
            records.append(tuple(row))
        return tuple(records)


def _decimal_tail_reference(
    weighted_samples: np.ndarray,
    *,
    spacing: float,
    precision: int,
) -> dict[str, object]:
    twiddles = _decimal_top_band_twiddles(precision)
    with localcontext() as context:
        context.prec = precision
        values = tuple(Decimal.from_float(float(value)) for value in weighted_samples)
        n_decimal = Decimal(TDG2_SAMPLE_COUNT)
        two_pi = Decimal(2) * Decimal(TDG2_DECIMAL_PI)
        spacing_decimal = Decimal.from_float(spacing)
        first_top = int(
            ceil(
                (1.0 - TDG2_UNRESOLVED_TOP_FRACTION)
                * (TDG2_SAMPLE_COUNT // 2)
            )
        )
        field_power = Decimal(0)
        derivative_power = Decimal(0)
        coefficients: list[tuple[Decimal, Decimal]] = []
        for offset, row in enumerate(twiddles):
            frequency_bin = first_top + offset
            real = sum(
                (value * cosine for value, (cosine, _) in zip(values, row, strict=True)),
                Decimal(0),
            )
            imaginary = sum(
                (value * sine for value, (_, sine) in zip(values, row, strict=True)),
                Decimal(0),
            )
            coefficients.append((real, imaginary))
            amplitude_squared = real * real + imaginary * imaginary
            multiplicity = Decimal(1 if frequency_bin == TDG2_SAMPLE_COUNT // 2 else 2)
            bin_power = multiplicity * amplitude_squared / (n_decimal * n_decimal)
            angular = (
                two_pi
                * Decimal(frequency_bin)
                / (n_decimal * spacing_decimal)
            )
            field_power += bin_power
            derivative_power += angular * angular * bin_power
        peak = max(
            ((real * real + imaginary * imaginary).sqrt() for real, imaginary in coefficients),
            default=Decimal(0),
        )
        return {
            "field": +field_power,
            "derivative": +derivative_power,
            "peak": +peak,
            "coefficients": tuple(coefficients),
        }


def _binary64_tail_reference(
    weighted_samples: np.ndarray,
    *,
    spacing: float,
) -> dict[str, object]:
    transformed = np.fft.rfft(weighted_samples)
    nyquist = TDG2_SAMPLE_COUNT // 2
    first_top = int(ceil((1.0 - TDG2_UNRESOLVED_TOP_FRACTION) * nyquist))
    coefficients = transformed[first_top : nyquist + 1]
    amplitudes = np.abs(coefficients)
    multiplicities = np.full(coefficients.size, 2.0, dtype=np.float64)
    multiplicities[-1] = 1.0
    powers = multiplicities * amplitudes**2 / float(TDG2_SAMPLE_COUNT**2)
    bins = np.arange(first_top, nyquist + 1, dtype=np.float64)
    angular = 2.0 * pi * bins / float(TDG2_SAMPLE_COUNT * spacing)
    return {
        "field": float(np.sum(powers)),
        "derivative": float(np.sum(angular**2 * powers)),
        "peak": float(np.max(amplitudes, initial=0.0)),
        "coefficients": coefficients,
    }


def _reference_debit(high: Decimal, low: Decimal, binary64: float) -> float:
    high_float = float(high)
    low_float = float(low)
    disagreement = max(abs(high_float - low_float), abs(high_float - binary64))
    outward = abs(float(np.nextafter(high_float, np.inf)) - high_float)
    return float(np.nextafter(disagreement + outward, np.inf))


def absolute_tail_observation(
    proper_times: object,
    values: object,
    *,
    round_trip_interpolation_infinity: Real,
) -> AbsoluteTailObservation:
    """Measure one aligned signal with independent arithmetic cross-checks."""

    proper = _array("proper_times", proper_times, ndim=1)
    data = _array("values", values, ndim=1)
    _require(
        proper.shape == (TDG2_SAMPLE_COUNT,)
        and data.shape == (TDG2_SAMPLE_COUNT,),
        "TDG2 signals must contain exactly 64 samples",
    )
    spacings = np.diff(proper)
    _require(np.all(spacings > 0.0), "TDG2 proper times must increase")
    spacing = float(spacings[0])
    tolerance = 4096.0 * np.finfo(np.float64).eps * max(1.0, abs(spacing))
    _require(
        np.allclose(spacings, spacing, rtol=0.0, atol=tolerance),
        "TDG2 aligned proper-time grid must be uniform",
    )
    interpolation_radius = _finite(
        "round_trip_interpolation_infinity",
        round_trip_interpolation_infinity,
    )
    _require(interpolation_radius >= 0.0, "TDG2 interpolation debit is negative")

    window = _compact_window(proper)
    weighted = data * window
    low = _decimal_tail_reference(
        weighted,
        spacing=spacing,
        precision=TDG2_DECIMAL_PRECISION_LOW,
    )
    high = _decimal_tail_reference(
        weighted,
        spacing=spacing,
        precision=TDG2_DECIMAL_PRECISION_HIGH,
    )
    binary = _binary64_tail_reference(weighted, spacing=spacing)

    coefficient_disagreement = 0.0
    for high_pair, low_pair, binary_value in zip(
        high["coefficients"],
        low["coefficients"],
        binary["coefficients"],
        strict=True,
    ):
        high_real, high_imaginary = (float(value) for value in high_pair)
        low_real, low_imaginary = (float(value) for value in low_pair)
        coefficient_disagreement = max(
            coefficient_disagreement,
            hypot(high_real - low_real, high_imaginary - low_imaginary),
            hypot(
                high_real - float(np.real(binary_value)),
                high_imaginary - float(np.imag(binary_value)),
            ),
        )
    peak = float(high["peak"])
    peak_outward = abs(float(np.nextafter(peak, np.inf)) - peak)
    coefficient_arithmetic_debit = float(
        np.nextafter(coefficient_disagreement + peak_outward, np.inf)
    )
    window_l1 = float(np.sum(np.abs(window)))
    coefficient_interpolation_debit = window_l1 * interpolation_radius

    first_top = int(
        ceil(
            (1.0 - TDG2_UNRESOLVED_TOP_FRACTION)
            * (TDG2_SAMPLE_COUNT // 2)
        )
    )
    high_coefficients = high["coefficients"]
    field_interpolation = 0.0
    derivative_interpolation = 0.0
    for offset, pair in enumerate(high_coefficients):
        frequency_bin = first_top + offset
        amplitude = hypot(float(pair[0]), float(pair[1]))
        multiplicity = 1.0 if frequency_bin == TDG2_SAMPLE_COUNT // 2 else 2.0
        perturbation = (
            multiplicity
            * (
                2.0 * amplitude * coefficient_interpolation_debit
                + coefficient_interpolation_debit**2
            )
            / float(TDG2_SAMPLE_COUNT**2)
        )
        field_interpolation += perturbation
        angular = 2.0 * pi * frequency_bin / float(TDG2_SAMPLE_COUNT * spacing)
        derivative_interpolation += angular**2 * perturbation

    field_observed = float(high["field"])
    derivative_observed = float(high["derivative"])
    field_arithmetic = _reference_debit(
        high["field"], low["field"], float(binary["field"])
    )
    derivative_arithmetic = _reference_debit(
        high["derivative"], low["derivative"], float(binary["derivative"])
    )
    field_total = field_arithmetic + field_interpolation
    derivative_total = derivative_arithmetic + derivative_interpolation
    field_lower, field_upper = _outward_interval(field_observed, field_total)
    derivative_lower, derivative_upper = _outward_interval(
        derivative_observed, derivative_total
    )
    floor_upper = (
        peak
        + coefficient_arithmetic_debit
        + coefficient_interpolation_debit
    )
    return AbsoluteTailObservation(
        sample_count=TDG2_SAMPLE_COUNT,
        top_band_first_bin=first_top,
        nyquist_bin=TDG2_SAMPLE_COUNT // 2,
        proper_time_start=float(proper[0]),
        proper_time_stop=float(proper[-1]),
        proper_time_spacing=spacing,
        round_trip_interpolation_infinity=interpolation_radius,
        window_L1=window_l1,
        peak_top_band_coefficient_amplitude=peak,
        binary64_peak_top_band_coefficient_amplitude=float(binary["peak"]),
        coefficient_arithmetic_debit=coefficient_arithmetic_debit,
        coefficient_interpolation_debit=coefficient_interpolation_debit,
        declared_absolute_amplitude_floor=TDG2_ABSOLUTE_AMPLITUDE_FLOOR,
        top_band_is_below_declared_amplitude_floor=(
            floor_upper <= TDG2_ABSOLUTE_AMPLITUDE_FLOOR
        ),
        field_tail_power=AbsolutePowerInterval(
            observed=field_observed,
            binary64_fft_value=float(binary["field"]),
            arithmetic_debit=field_arithmetic,
            interpolation_debit=field_interpolation,
            total_debit=field_total,
            lower=field_lower,
            upper=field_upper,
        ),
        derivative_tail_power=AbsolutePowerInterval(
            observed=derivative_observed,
            binary64_fft_value=float(binary["derivative"]),
            arithmetic_debit=derivative_arithmetic,
            interpolation_debit=derivative_interpolation,
            total_debit=derivative_total,
            lower=derivative_lower,
            upper=derivative_upper,
        ),
    )


def _common_overlap_histories(
    proper_times_by_resolution: Sequence[np.ndarray],
    values_by_resolution: Sequence[np.ndarray],
) -> tuple[np.ndarray, tuple[np.ndarray, ...], tuple[float, ...]]:
    start = max(float(record[0]) for record in proper_times_by_resolution)
    stop = min(float(record[-1]) for record in proper_times_by_resolution)
    _require(stop > start, "TDG2 histories have no common proper-time span")
    common = np.linspace(start, stop, TDG2_SAMPLE_COUNT, dtype=np.float64)
    aligned: list[np.ndarray] = []
    debits: list[float] = []
    for proper, values in zip(
        proper_times_by_resolution, values_by_resolution, strict=True
    ):
        aligned_values = np.interp(common, proper, values)
        mask = (proper >= start) & (proper <= stop)
        _require(np.count_nonzero(mask) >= 16, "TDG2 common overlap is too short")
        reverse = np.interp(proper[mask], common, aligned_values)
        debits.append(
            float(np.max(np.abs(reverse - values[mask]), initial=0.0))
        )
        aligned.append(aligned_values)
    return common, tuple(aligned), tuple(debits)


def _method_counts(method: str, point_counts: Sequence[int]) -> tuple[int, int, int]:
    _require(method in TDG2_METHOD_POINT_COUNTS, "TDG2 method differs")
    counts = tuple(point_counts)
    _require(counts == TDG2_METHOD_POINT_COUNTS[method], "TDG2 method ladder differs")
    return counts


def _power_interval(
    observation: AbsoluteTailObservation,
    measure: str,
) -> AbsolutePowerInterval:
    if measure not in TDG2_MEASURE_NAMES:
        raise ValueError("unknown TDG2 measure")
    return getattr(observation, measure)


def _order_lower(
    coarse_lower: float,
    coarse_upper: float,
    fine_lower: float,
    fine_upper: float,
    *,
    refinement_ratio: float,
) -> float | None:
    if coarse_lower <= 0.0 or fine_upper <= 0.0:
        return None
    if coarse_lower <= fine_upper:
        return None
    return log(coarse_lower / fine_upper) / log(refinement_ratio)


def classify_absolute_power_ladder(
    observations: Sequence[AbsoluteTailObservation],
    *,
    method: str,
    point_counts: Sequence[int],
    measure: str,
) -> AbsolutePowerClassification:
    """Classify one absolute power ladder without defining a pass/fail gate."""

    counts = _method_counts(method, point_counts)
    _require(len(observations) == 3, "TDG2 requires three observations")
    records = tuple(_power_interval(value, measure) for value in observations)
    ratio_a = (counts[1] - 1) / (counts[0] - 1)
    ratio_b = (counts[2] - 1) / (counts[1] - 1)
    _require(ratio_a == ratio_b and ratio_a > 1.0, "TDG2 refinement ratios differ")
    observed = tuple(value.observed for value in records)
    lower = tuple(value.lower for value in records)
    upper = tuple(value.upper for value in records)
    every_floor = all(
        value.top_band_is_below_declared_amplitude_floor
        for value in observations
    )
    any_conditional_debit = any(
        value.interpolation_debit > 0.0 or value.arithmetic_debit > 0.0
        for value in records
    )

    zero_orders = (
        _order_lower(lower[0], upper[0], lower[1], upper[1], refinement_ratio=ratio_a),
        _order_lower(lower[1], upper[1], lower[2], upper[2], refinement_ratio=ratio_a),
    )
    difference_order: float | None = None
    continuum_lower: float | None = None
    enclosure_controls = False

    if every_floor or lower[2] == 0.0:
        classification = "instrument_floor_or_interpolation_enclosure_dominated"
        enclosure_controls = True
    elif all(
        value is not None and value >= TDG2_MINIMUM_ZERO_POWER_ORDER
        for value in zero_orders
    ):
        classification = "contracts_to_zero_at_required_power_order"
    else:
        common_lower = max(lower)
        common_upper = min(upper)
        if common_lower <= common_upper and common_lower > 0.0:
            classification = "nonzero_resolution_independent_within_enclosure"
            continuum_lower = common_lower
        else:
            decreasing_coarse = (lower[0] - upper[1], upper[0] - lower[1])
            decreasing_fine = (lower[1] - upper[2], upper[1] - lower[2])
            increasing_coarse = (lower[1] - upper[0], upper[1] - lower[0])
            increasing_fine = (lower[2] - upper[1], upper[2] - lower[1])
            direction: str | None = None
            coarse_difference: tuple[float, float] | None = None
            fine_difference: tuple[float, float] | None = None
            if decreasing_coarse[0] > 0.0 and decreasing_fine[0] > 0.0:
                direction = "decreasing"
                coarse_difference = decreasing_coarse
                fine_difference = decreasing_fine
            elif increasing_coarse[0] > 0.0 and increasing_fine[0] > 0.0:
                direction = "increasing"
                coarse_difference = increasing_coarse
                fine_difference = increasing_fine

            if direction is not None and coarse_difference and fine_difference:
                if coarse_difference[0] > fine_difference[1] > 0.0:
                    difference_order = log(
                        coarse_difference[0] / fine_difference[1]
                    ) / log(ratio_a)
                if (
                    difference_order is not None
                    and difference_order >= TDG2_MINIMUM_HISTORY_ORDER
                ):
                    if direction == "decreasing":
                        contraction = ratio_a ** (-TDG2_MINIMUM_HISTORY_ORDER)
                        remaining_upper = (
                            contraction
                            / (1.0 - contraction)
                            * fine_difference[1]
                        )
                        continuum_lower = max(0.0, lower[2] - remaining_upper)
                    else:
                        continuum_lower = lower[2]
                    if continuum_lower > 0.0:
                        classification = "converges_to_nonzero_resolved_power"
                    else:
                        classification = "unresolved_absolute_tail_behavior"
                else:
                    classification = "unresolved_absolute_tail_behavior"
            elif any_conditional_debit and (
                max(lower[0], lower[1]) <= min(upper[0], upper[1])
                or max(lower[1], lower[2]) <= min(upper[1], upper[2])
            ):
                classification = "instrument_floor_or_interpolation_enclosure_dominated"
                enclosure_controls = True
            else:
                classification = "unresolved_absolute_tail_behavior"

    return AbsolutePowerClassification(
        measure=measure,
        classification=classification,
        point_counts=counts,
        refinement_ratio=ratio_a,
        observed_powers=observed,
        lower_bounds=lower,
        upper_bounds=upper,
        zero_power_order_lower_bounds=zero_orders,
        difference_order_lower_bound=difference_order,
        conservative_continuum_lower_bound=continuum_lower,
        every_resolution_below_declared_amplitude_floor=every_floor,
        enclosure_controls_classification=enclosure_controls,
    )


def diagnose_absolute_tail_ladder(
    proper_times_by_resolution: Sequence[object],
    values_by_resolution: Sequence[object],
    *,
    method: str,
    point_counts: Sequence[int],
) -> dict[str, Any]:
    """Apply the frozen absolute-tail discriminator to one scalar history ladder."""

    counts = _method_counts(method, point_counts)
    times = tuple(
        _array("proper_times", value, ndim=1)
        for value in proper_times_by_resolution
    )
    values = tuple(_array("values", value, ndim=1) for value in values_by_resolution)
    _require(
        len(times) == 3
        and len(values) == 3
        and all(item.shape == (TDG2_SAMPLE_COUNT,) for item in (*times, *values))
        and all(np.all(np.diff(item) > 0.0) for item in times),
        "TDG2 ladder inputs must be three increasing 64-sample records",
    )
    common, aligned, debits = _common_overlap_histories(times, values)
    observations = tuple(
        absolute_tail_observation(
            common,
            data,
            round_trip_interpolation_infinity=debit,
        )
        for data, debit in zip(aligned, debits, strict=True)
    )
    classifications = {
        measure: asdict(
            classify_absolute_power_ladder(
                observations,
                method=method,
                point_counts=counts,
                measure=measure,
            )
        )
        for measure in TDG2_MEASURE_NAMES
    }
    return {
        "method": method,
        "point_counts": list(counts),
        "common_proper_time_start": float(common[0]),
        "common_proper_time_stop": float(common[-1]),
        "observations": [asdict(value) for value in observations],
        "classifications": classifications,
        "round_trip_interpolation_debit_is_a_conditional_perturbation_radius": True,
        "classification_is_not_a_replacement_admission": True,
    }


def _synthetic_times(*, jitter: float = 0.0) -> tuple[np.ndarray, ...]:
    base = np.linspace(0.0, 63.0 / 16.0, TDG2_SAMPLE_COUNT, dtype=np.float64)
    if jitter == 0.0:
        return (base.copy(), base.copy(), base.copy())
    normalized = np.linspace(0.0, 1.0, TDG2_SAMPLE_COUNT, dtype=np.float64)
    records = [base]
    for phase in (0.0, 0.37):
        perturbation = jitter * np.sin(2.0 * pi * (normalized + phase))
        perturbation[0] = 0.0
        perturbation[-1] = 0.0
        candidate = base + perturbation
        _require(np.all(np.diff(candidate) > 0.0), "synthetic TDG2 time jitter folded")
        records.append(candidate)
    return tuple(records)


def synthetic_absolute_tail_controls() -> dict[str, Any]:
    """Prove the frozen classifier separates its intended classes without raw data."""

    times = _synthetic_times()
    coordinate = times[0]
    sample_index = np.arange(TDG2_SAMPLE_COUNT, dtype=np.float64)
    carrier = np.sin(2.0 * pi * 29.0 * sample_index / TDG2_SAMPLE_COUNT)
    counts = TDG2_METHOD_POINT_COUNTS["RK4"]

    controls: dict[str, tuple[tuple[np.ndarray, ...], tuple[np.ndarray, ...]]] = {
        "exact_zero": (times, tuple(np.zeros_like(coordinate) for _ in range(3))),
        "subfloor_nonzero": (
            times,
            tuple(1.0e-40 * carrier for _ in range(3)),
        ),
        "power_contracts_to_zero": (
            times,
            tuple(amplitude * carrier for amplitude in (1.0e-3, 1.25e-4, 1.5625e-5)),
        ),
        "nonzero_resolution_independent": (
            times,
            tuple(1.0e-3 * carrier for _ in range(3)),
        ),
        "nonzero_convergent": (
            times,
            tuple(
                1.0e-3 * amplitude * carrier
                for amplitude in (1.25, 1.0625, 1.015625)
            ),
        ),
        "unresolved_nonmonotone": (
            times,
            tuple(
                1.0e-3 * amplitude * carrier
                for amplitude in (1.0, 0.5, 0.75)
            ),
        ),
    }
    jittered_times = _synthetic_times(jitter=2.0e-2)
    jittered_values = tuple(
        1.0e-8 * np.sin(
            2.0
            * pi
            * 29.0
            * (record - record[0])
            / (record[-1] - record[0])
        )
        for record in jittered_times
    )
    controls["interpolation_enclosure_dominated"] = (
        jittered_times,
        jittered_values,
    )
    _require(
        tuple(controls) == TDG2_SYNTHETIC_CONTROL_NAMES,
        "TDG2 synthetic control order differs",
    )

    expected = {
        "exact_zero": "instrument_floor_or_interpolation_enclosure_dominated",
        "subfloor_nonzero": "instrument_floor_or_interpolation_enclosure_dominated",
        "power_contracts_to_zero": "contracts_to_zero_at_required_power_order",
        "nonzero_resolution_independent": "nonzero_resolution_independent_within_enclosure",
        "nonzero_convergent": "converges_to_nonzero_resolved_power",
        "unresolved_nonmonotone": "unresolved_absolute_tail_behavior",
        "interpolation_enclosure_dominated": "instrument_floor_or_interpolation_enclosure_dominated",
    }
    records: dict[str, Any] = {}
    all_matched = True
    for name, (control_times, control_values) in controls.items():
        diagnosis = diagnose_absolute_tail_ladder(
            control_times,
            control_values,
            method="RK4",
            point_counts=counts,
        )
        observed_classes = {
            measure: diagnosis["classifications"][measure]["classification"]
            for measure in TDG2_MEASURE_NAMES
        }
        matched = all(value == expected[name] for value in observed_classes.values())
        all_matched = all_matched and matched
        records[name] = {
            "expected_classification": expected[name],
            "observed_classifications": observed_classes,
            "classification_matched": matched,
            "diagnosis": diagnosis,
        }
    return {
        "controls": records,
        "all_control_classes_separated": all_matched,
        "actual_terminal_histories_consumed": False,
        "replacement_temporal_admission_defined": False,
        "mathematical_consequence": (
            "absolute_unnormalized_tail_power_preserves_amplitude_convergence_"
            "and_separates_zero_nonzero_enclosure_and_unresolved_controls"
            if all_matched
            else "TDG2_synthetic_discriminator_not_validated"
        ),
    }


__all__ = [
    "AbsolutePowerClassification",
    "AbsolutePowerInterval",
    "AbsoluteTailObservation",
    "TDG2_ABSOLUTE_AMPLITUDE_FLOOR",
    "TDG2_CLASS_NAMES",
    "TDG2_DECIMAL_PRECISION_HIGH",
    "TDG2_DECIMAL_PRECISION_LOW",
    "TDG2_FIELD_NAMES",
    "TDG2_MEASURE_NAMES",
    "TDG2_METHOD_POINT_COUNTS",
    "TDG2_MINIMUM_HISTORY_ORDER",
    "TDG2_MINIMUM_ZERO_POWER_ORDER",
    "TDG2_SAMPLE_COUNT",
    "TDG2_SYNTHETIC_CONTROL_NAMES",
    "TDG2_TAPER_FRACTION",
    "TDG2_TRACER_COUNT",
    "TDG2_UNRESOLVED_TOP_FRACTION",
    "absolute_tail_observation",
    "classify_absolute_power_ladder",
    "diagnose_absolute_tail_ladder",
    "synthetic_absolute_tail_controls",
]
