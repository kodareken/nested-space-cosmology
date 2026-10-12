"""Prospective native-grid discriminator for TDG3.

TDG2 showed that its common-grid round-trip interpolation enclosure dominates
most unresolved terminal-history classifications.  This module asks a narrower
question without loading those histories: can the same absolute top-band shape
be measured directly on each native nonuniform proper-time grid?

Two estimators are frozen.  The primary estimator analytically integrates the
piecewise-linear interpolant of the native, compact-windowed samples.  The
comparator applies native Voronoi/trapezoid weights directly to those samples.
Neither estimator resamples a value onto a common time grid.  Decimal
calculations at two precisions independently audit the binary64 coefficients
and powers.  Agreement between the two estimators is descriptive evidence
only: neither the piecewise-linear surrogate nor the quadrature rule is claimed
to enclose the unknown continuum history, and no result can retroactively pass
PROTO13 or define a replacement temporal admission.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from decimal import Decimal, localcontext
from math import fsum, hypot, isfinite, log, pi
from numbers import Real
from typing import Any, Mapping, Sequence

import numpy as np

from .health_monitor import compact_vacuum_buffer_window
from .tdg2_absolute_tail_diagnosis import (
    TDG2_ABSOLUTE_AMPLITUDE_FLOOR,
    TDG2_DECIMAL_PI,
    TDG2_DECIMAL_PRECISION_HIGH,
    TDG2_DECIMAL_PRECISION_LOW,
    TDG2_FIELD_NAMES,
    TDG2_METHOD_POINT_COUNTS,
    TDG2_MINIMUM_HISTORY_ORDER,
    TDG2_MINIMUM_ZERO_POWER_ORDER,
    TDG2_SAMPLE_COUNT,
    TDG2_TAPER_FRACTION,
    _decimal_sin_cos,
    diagnose_absolute_tail_ladder,
)


TDG3_SAMPLE_COUNT = TDG2_SAMPLE_COUNT
TDG3_FIELD_NAMES = TDG2_FIELD_NAMES
TDG3_METHOD_POINT_COUNTS = TDG2_METHOD_POINT_COUNTS
TDG3_TOP_BINS = (28, 29, 30, 31, 32)
TDG3_ESTIMATOR_NAMES = (
    "piecewise_linear_exact_integral",
    "native_trapezoid_quadrature",
)
TDG3_MEASURE_NAMES = (
    "phase_field_tail_power",
    "phase_derivative_tail_power",
)
TDG3_POWER_CLASS_NAMES = (
    "contracts_to_zero_at_required_power_order",
    "nonzero_resolution_independent_within_arithmetic_enclosure",
    "converges_to_nonzero_native_power",
    "arithmetic_floor_dominated",
    "unresolved_native_power_behavior",
)
TDG3_COMBINED_CLASS_NAMES = (
    "native_estimators_agree_zero_contracting",
    "native_estimators_agree_nonzero",
    "native_estimators_agree_floor_dominated",
    "native_estimators_agree_unresolved",
    "native_estimators_disagree",
)
TDG3_NORMALIZED_COEFFICIENT_FLOOR = (
    TDG2_ABSOLUTE_AMPLITUDE_FLOOR / TDG3_SAMPLE_COUNT
)
TDG3_SYNTHETIC_CONTROL_NAMES = (
    "exact_zero",
    "subfloor_nonzero",
    "power_contracts_to_zero",
    "nonzero_resolution_independent",
    "nonzero_convergent",
    "unresolved_nonmonotone",
    "common_grid_interpolation_owner",
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


def _array(name: str, value: object) -> np.ndarray:
    answer = np.asarray(value, dtype=np.float64)
    if answer.ndim != 1 or answer.shape != (TDG3_SAMPLE_COUNT,):
        raise ValueError(f"{name} must contain exactly 64 scalars")
    if not np.all(np.isfinite(answer)):
        raise ValueError(f"{name} must contain only finite values")
    return answer


def _outward_interval(value: float, debit: float) -> tuple[float, float]:
    centre = _finite("interval value", value)
    radius = _finite("interval debit", debit)
    if centre < 0.0 or radius < 0.0:
        raise ValueError("TDG3 power intervals must be nonnegative")
    lower = max(0.0, centre - radius)
    upper = centre + radius
    return (
        float(np.nextafter(lower, -np.inf)) if lower > 0.0 else 0.0,
        float(np.nextafter(upper, np.inf)),
    )


def _complex_multiply(
    left: tuple[Decimal, Decimal], right: tuple[Decimal, Decimal]
) -> tuple[Decimal, Decimal]:
    a, b = left
    c, d = right
    return (a * c - b * d, a * d + b * c)


def _decimal_exponential(
    angle: Decimal, *, precision: int
) -> tuple[Decimal, Decimal]:
    sine, cosine = _decimal_sin_cos(angle, precision=precision)
    return (cosine, -sine)


def _native_inputs(
    proper_times: object, values: object
) -> tuple[np.ndarray, np.ndarray, float, np.ndarray]:
    proper = _array("proper_times", proper_times)
    data = _array("values", values)
    spacings = np.diff(proper)
    _require(np.all(spacings > 0.0), "TDG3 proper times must increase")
    span = float(proper[-1] - proper[0])
    _require(span > 0.0, "TDG3 proper-time span must be positive")
    phase = (proper - proper[0]) / span
    _require(phase[0] == 0.0 and phase[-1] == 1.0, "TDG3 phase map differs")
    window = compact_vacuum_buffer_window(
        phase,
        support_minimum=0.0,
        support_maximum=1.0,
        taper_width=TDG2_TAPER_FRACTION,
    )
    weighted = data * window
    return phase, weighted, span, window


def _binary_piecewise_coefficients(
    phase: np.ndarray, weighted: np.ndarray
) -> tuple[complex, ...]:
    coefficients: list[complex] = []
    for frequency_bin in TDG3_TOP_BINS:
        angular = 2.0 * pi * frequency_bin
        real_terms: list[float] = []
        imaginary_terms: list[float] = []
        for index in range(TDG3_SAMPLE_COUNT - 1):
            x0 = float(phase[index])
            width = float(phase[index + 1] - phase[index])
            y0 = float(weighted[index])
            slope = float(weighted[index + 1] - weighted[index]) / width
            e0 = np.exp(-1j * angular * x0)
            eh = np.exp(-1j * angular * width)
            i0 = e0 * (1.0 - eh) / (1j * angular)
            i1 = e0 * (eh * (1.0 + 1j * angular * width) - 1.0) / (
                angular * angular
            )
            contribution = y0 * i0 + slope * i1
            real_terms.append(float(contribution.real))
            imaginary_terms.append(float(contribution.imag))
        coefficients.append(complex(fsum(real_terms), fsum(imaginary_terms)))
    return tuple(coefficients)


def _decimal_piecewise_coefficients(
    phase: np.ndarray, weighted: np.ndarray, *, precision: int
) -> tuple[tuple[Decimal, Decimal], ...]:
    if precision not in (TDG2_DECIMAL_PRECISION_LOW, TDG2_DECIMAL_PRECISION_HIGH):
        raise ValueError("TDG3 Decimal precision differs")
    with localcontext() as context:
        context.prec = precision
        two_pi = Decimal(2) * Decimal(TDG2_DECIMAL_PI)
        x = tuple(Decimal.from_float(float(value)) for value in phase)
        y = tuple(Decimal.from_float(float(value)) for value in weighted)
        coefficients: list[tuple[Decimal, Decimal]] = []
        for frequency_bin in TDG3_TOP_BINS:
            angular = two_pi * Decimal(frequency_bin)
            angular_squared = angular * angular
            real = Decimal(0)
            imaginary = Decimal(0)
            for index in range(TDG3_SAMPLE_COUNT - 1):
                x0 = x[index]
                width = x[index + 1] - x[index]
                slope = (y[index + 1] - y[index]) / width
                e0 = _decimal_exponential(
                    angular * x0, precision=precision
                )
                eh = _decimal_exponential(
                    angular * width, precision=precision
                )
                one_minus_eh = (Decimal(1) - eh[0], -eh[1])
                i0_local = (
                    one_minus_eh[1] / angular,
                    -one_minus_eh[0] / angular,
                )
                i0 = _complex_multiply(e0, i0_local)
                phase_width = angular * width
                eh_linear = (
                    eh[0] - eh[1] * phase_width - Decimal(1),
                    eh[0] * phase_width + eh[1],
                )
                i1 = _complex_multiply(
                    e0,
                    (
                        eh_linear[0] / angular_squared,
                        eh_linear[1] / angular_squared,
                    ),
                )
                real += y[index] * i0[0] + slope * i1[0]
                imaginary += y[index] * i0[1] + slope * i1[1]
            coefficients.append((+real, +imaginary))
        return tuple(coefficients)


def _native_trapezoid_weights(phase: np.ndarray) -> np.ndarray:
    weights = np.empty_like(phase)
    weights[0] = 0.5 * (phase[1] - phase[0])
    weights[-1] = 0.5 * (phase[-1] - phase[-2])
    weights[1:-1] = 0.5 * (phase[2:] - phase[:-2])
    _require(np.all(weights > 0.0), "TDG3 native quadrature weights differ")
    tolerance = 4096.0 * np.finfo(np.float64).eps
    _require(abs(float(np.sum(weights)) - 1.0) <= tolerance, "TDG3 weights do not close")
    return weights


def _binary_trapezoid_coefficients(
    phase: np.ndarray, weighted: np.ndarray
) -> tuple[complex, ...]:
    weights = _native_trapezoid_weights(phase)
    coefficients: list[complex] = []
    for frequency_bin in TDG3_TOP_BINS:
        angular = 2.0 * pi * frequency_bin
        values = weights * weighted
        real = fsum(
            float(value) * float(np.cos(angular * float(position)))
            for position, value in zip(phase, values, strict=True)
        )
        imaginary = fsum(
            -float(value) * float(np.sin(angular * float(position)))
            for position, value in zip(phase, values, strict=True)
        )
        coefficients.append(complex(real, imaginary))
    return tuple(coefficients)


def _decimal_trapezoid_coefficients(
    phase: np.ndarray, weighted: np.ndarray, *, precision: int
) -> tuple[tuple[Decimal, Decimal], ...]:
    if precision not in (TDG2_DECIMAL_PRECISION_LOW, TDG2_DECIMAL_PRECISION_HIGH):
        raise ValueError("TDG3 Decimal precision differs")
    weights = _native_trapezoid_weights(phase)
    with localcontext() as context:
        context.prec = precision
        two_pi = Decimal(2) * Decimal(TDG2_DECIMAL_PI)
        x = tuple(Decimal.from_float(float(value)) for value in phase)
        y = tuple(Decimal.from_float(float(value)) for value in weighted)
        q = tuple(Decimal.from_float(float(value)) for value in weights)
        coefficients: list[tuple[Decimal, Decimal]] = []
        for frequency_bin in TDG3_TOP_BINS:
            angular = two_pi * Decimal(frequency_bin)
            real = Decimal(0)
            imaginary = Decimal(0)
            for position, value, weight in zip(x, y, q, strict=True):
                exponential = _decimal_exponential(
                    angular * position, precision=precision
                )
                real += weight * value * exponential[0]
                imaginary += weight * value * exponential[1]
            coefficients.append((+real, +imaginary))
        return tuple(coefficients)


def _coefficient_bundle(
    phase: np.ndarray, weighted: np.ndarray, *, estimator: str
) -> tuple[
    tuple[complex, ...],
    tuple[tuple[Decimal, Decimal], ...],
    tuple[tuple[Decimal, Decimal], ...],
]:
    if estimator == "piecewise_linear_exact_integral":
        binary = _binary_piecewise_coefficients(phase, weighted)
        low = _decimal_piecewise_coefficients(
            phase, weighted, precision=TDG2_DECIMAL_PRECISION_LOW
        )
        high = _decimal_piecewise_coefficients(
            phase, weighted, precision=TDG2_DECIMAL_PRECISION_HIGH
        )
    elif estimator == "native_trapezoid_quadrature":
        binary = _binary_trapezoid_coefficients(phase, weighted)
        low = _decimal_trapezoid_coefficients(
            phase, weighted, precision=TDG2_DECIMAL_PRECISION_LOW
        )
        high = _decimal_trapezoid_coefficients(
            phase, weighted, precision=TDG2_DECIMAL_PRECISION_HIGH
        )
    else:
        raise ValueError("unknown TDG3 native estimator")
    return binary, low, high


def _decimal_powers(
    coefficients: Sequence[tuple[Decimal, Decimal]]
) -> tuple[Decimal, Decimal, Decimal]:
    with localcontext() as context:
        context.prec = TDG2_DECIMAL_PRECISION_HIGH
        two_pi = Decimal(2) * Decimal(TDG2_DECIMAL_PI)
        field = Decimal(0)
        derivative = Decimal(0)
        peak = Decimal(0)
        for frequency_bin, (real, imaginary) in zip(
            TDG3_TOP_BINS, coefficients, strict=True
        ):
            amplitude_squared = real * real + imaginary * imaginary
            amplitude = amplitude_squared.sqrt()
            peak = max(peak, amplitude)
            field += Decimal(2) * amplitude_squared
            angular = two_pi * Decimal(frequency_bin)
            derivative += Decimal(2) * angular * angular * amplitude_squared
        return (+field, +derivative, +peak)


def _binary_powers(coefficients: Sequence[complex]) -> tuple[float, float, float]:
    amplitudes = [abs(value) for value in coefficients]
    field = fsum(2.0 * value * value for value in amplitudes)
    derivative = fsum(
        2.0 * (2.0 * pi * frequency_bin) ** 2 * value * value
        for frequency_bin, value in zip(TDG3_TOP_BINS, amplitudes, strict=True)
    )
    return field, derivative, max(amplitudes, default=0.0)


def _reference_debit(high: Decimal, low: Decimal, binary64: float) -> float:
    high_float = float(high)
    disagreement = max(abs(high_float - float(low)), abs(high_float - binary64))
    outward = abs(float(np.nextafter(high_float, np.inf)) - high_float)
    return float(np.nextafter(disagreement + outward, np.inf))


@dataclass(frozen=True, slots=True)
class NativePowerInterval:
    observed: float
    binary64_value: float
    arithmetic_debit: float
    lower: float
    upper: float

    def __post_init__(self) -> None:
        for name in self.__dataclass_fields__:
            value = _finite(name, getattr(self, name))
            if value < 0.0:
                raise ValueError("TDG3 power evidence must be nonnegative")
            object.__setattr__(self, name, value)
        if not self.lower <= self.observed <= self.upper:
            raise ValueError("TDG3 power observation lies outside its interval")


@dataclass(frozen=True, slots=True)
class NativeTailObservation:
    estimator: str
    sample_count: int
    top_bins: tuple[int, ...]
    proper_time_start: float
    proper_time_stop: float
    proper_time_span: float
    minimum_normalized_spacing: float
    maximum_normalized_spacing: float
    normalized_spacing_ratio: float
    peak_coefficient_amplitude: float
    binary64_peak_coefficient_amplitude: float
    coefficient_arithmetic_debit: float
    declared_normalized_coefficient_floor: float
    top_band_is_below_declared_amplitude_floor: bool
    phase_field_tail_power: NativePowerInterval
    phase_derivative_tail_power: NativePowerInterval
    common_grid_resampling_performed: bool = False
    continuum_surrogate_enclosure_claimed: bool = False

    def __post_init__(self) -> None:
        if self.estimator not in TDG3_ESTIMATOR_NAMES:
            raise ValueError("unknown TDG3 estimator")
        if self.sample_count != TDG3_SAMPLE_COUNT or self.top_bins != TDG3_TOP_BINS:
            raise ValueError("TDG3 observation shape differs")
        for name in (
            "proper_time_start",
            "proper_time_stop",
            "proper_time_span",
            "minimum_normalized_spacing",
            "maximum_normalized_spacing",
            "normalized_spacing_ratio",
            "peak_coefficient_amplitude",
            "binary64_peak_coefficient_amplitude",
            "coefficient_arithmetic_debit",
            "declared_normalized_coefficient_floor",
        ):
            object.__setattr__(self, name, _finite(name, getattr(self, name)))
        if (
            self.proper_time_stop <= self.proper_time_start
            or self.proper_time_span <= 0.0
            or self.minimum_normalized_spacing <= 0.0
            or self.maximum_normalized_spacing < self.minimum_normalized_spacing
            or self.normalized_spacing_ratio < 1.0
        ):
            raise ValueError("TDG3 observation geometry differs")
        if self.common_grid_resampling_performed or self.continuum_surrogate_enclosure_claimed:
            raise ValueError("TDG3 may not claim resampling or continuum enclosure")


def native_tail_observation(
    proper_times: object, values: object, *, estimator: str
) -> NativeTailObservation:
    """Measure native-grid top-band phase power with no common-grid resampling."""

    phase, weighted, span, _ = _native_inputs(proper_times, values)
    binary, low_coefficients, high_coefficients = _coefficient_bundle(
        phase, weighted, estimator=estimator
    )
    low_field, low_derivative, low_peak = _decimal_powers(low_coefficients)
    high_field, high_derivative, high_peak = _decimal_powers(high_coefficients)
    binary_field, binary_derivative, binary_peak = _binary_powers(binary)
    coefficient_disagreement = 0.0
    for binary_value, low_pair, high_pair in zip(
        binary, low_coefficients, high_coefficients, strict=True
    ):
        coefficient_disagreement = max(
            coefficient_disagreement,
            hypot(
                float(high_pair[0] - low_pair[0]),
                float(high_pair[1] - low_pair[1]),
            ),
            hypot(
                float(high_pair[0]) - binary_value.real,
                float(high_pair[1]) - binary_value.imag,
            ),
        )
    peak = float(high_peak)
    coefficient_debit = float(
        np.nextafter(
            coefficient_disagreement
            + abs(float(np.nextafter(peak, np.inf)) - peak),
            np.inf,
        )
    )
    field_debit = _reference_debit(high_field, low_field, binary_field)
    derivative_debit = _reference_debit(
        high_derivative, low_derivative, binary_derivative
    )
    field_lower, field_upper = _outward_interval(float(high_field), field_debit)
    derivative_lower, derivative_upper = _outward_interval(
        float(high_derivative), derivative_debit
    )
    spacings = np.diff(phase)
    floor_upper = peak + coefficient_debit
    return NativeTailObservation(
        estimator=estimator,
        sample_count=TDG3_SAMPLE_COUNT,
        top_bins=TDG3_TOP_BINS,
        proper_time_start=float(np.asarray(proper_times)[0]),
        proper_time_stop=float(np.asarray(proper_times)[-1]),
        proper_time_span=span,
        minimum_normalized_spacing=float(np.min(spacings)),
        maximum_normalized_spacing=float(np.max(spacings)),
        normalized_spacing_ratio=float(np.max(spacings) / np.min(spacings)),
        peak_coefficient_amplitude=peak,
        binary64_peak_coefficient_amplitude=binary_peak,
        coefficient_arithmetic_debit=coefficient_debit,
        declared_normalized_coefficient_floor=TDG3_NORMALIZED_COEFFICIENT_FLOOR,
        top_band_is_below_declared_amplitude_floor=(
            floor_upper <= TDG3_NORMALIZED_COEFFICIENT_FLOOR
        ),
        phase_field_tail_power=NativePowerInterval(
            observed=float(high_field),
            binary64_value=binary_field,
            arithmetic_debit=field_debit,
            lower=field_lower,
            upper=field_upper,
        ),
        phase_derivative_tail_power=NativePowerInterval(
            observed=float(high_derivative),
            binary64_value=binary_derivative,
            arithmetic_debit=derivative_debit,
            lower=derivative_lower,
            upper=derivative_upper,
        ),
    )


def _method_counts(method: str, point_counts: Sequence[int]) -> tuple[int, int, int]:
    _require(method in TDG3_METHOD_POINT_COUNTS, "TDG3 method differs")
    counts = tuple(point_counts)
    _require(counts == TDG3_METHOD_POINT_COUNTS[method], "TDG3 method ladder differs")
    return counts


def _order_lower(
    coarse_lower: float,
    fine_upper: float,
    *,
    refinement_ratio: float,
) -> float | None:
    if coarse_lower <= 0.0 or fine_upper <= 0.0 or coarse_lower <= fine_upper:
        return None
    return log(coarse_lower / fine_upper) / log(refinement_ratio)


def classify_native_power_ladder(
    observations: Sequence[NativeTailObservation],
    *,
    method: str,
    point_counts: Sequence[int],
    measure: str,
) -> dict[str, Any]:
    counts = _method_counts(method, point_counts)
    _require(len(observations) == 3, "TDG3 requires three observations")
    _require(measure in TDG3_MEASURE_NAMES, "unknown TDG3 measure")
    estimators = {record.estimator for record in observations}
    _require(len(estimators) == 1, "TDG3 ladder mixes estimators")
    intervals = tuple(getattr(record, measure) for record in observations)
    ratio_a = (counts[1] - 1) / (counts[0] - 1)
    ratio_b = (counts[2] - 1) / (counts[1] - 1)
    _require(ratio_a == ratio_b and ratio_a > 1.0, "TDG3 refinement ratios differ")
    observed = tuple(record.observed for record in intervals)
    lower = tuple(record.lower for record in intervals)
    upper = tuple(record.upper for record in intervals)
    zero_orders = (
        _order_lower(lower[0], upper[1], refinement_ratio=ratio_a),
        _order_lower(lower[1], upper[2], refinement_ratio=ratio_a),
    )
    difference_order: float | None = None
    continuum_lower: float | None = None
    every_floor = all(
        record.top_band_is_below_declared_amplitude_floor
        for record in observations
    )
    if every_floor or lower[2] == 0.0:
        classification = "arithmetic_floor_dominated"
    elif all(
        order is not None and order >= TDG2_MINIMUM_ZERO_POWER_ORDER
        for order in zero_orders
    ):
        classification = "contracts_to_zero_at_required_power_order"
    else:
        common_lower = max(lower)
        common_upper = min(upper)
        if common_lower <= common_upper and common_lower > 0.0:
            classification = "nonzero_resolution_independent_within_arithmetic_enclosure"
            continuum_lower = common_lower
        else:
            directions = []
            for sign in (1.0, -1.0):
                coarse = (sign * lower[1] - sign * upper[0], sign * upper[1] - sign * lower[0])
                fine = (sign * lower[2] - sign * upper[1], sign * upper[2] - sign * lower[1])
                if coarse[0] > 0.0 and fine[0] > 0.0:
                    directions.append((sign, coarse, fine))
            if len(directions) == 1:
                sign, coarse, fine = directions[0]
                if coarse[0] > fine[1] > 0.0:
                    difference_order = log(coarse[0] / fine[1]) / log(ratio_a)
                if (
                    difference_order is not None
                    and difference_order >= TDG2_MINIMUM_HISTORY_ORDER
                ):
                    if sign < 0.0:
                        contraction = ratio_a ** (-TDG2_MINIMUM_HISTORY_ORDER)
                        remainder = contraction / (1.0 - contraction) * fine[1]
                        continuum_lower = max(0.0, lower[2] - remainder)
                    else:
                        continuum_lower = lower[2]
                    classification = (
                        "converges_to_nonzero_native_power"
                        if continuum_lower > 0.0
                        else "unresolved_native_power_behavior"
                    )
                else:
                    classification = "unresolved_native_power_behavior"
            else:
                classification = "unresolved_native_power_behavior"
    _require(classification in TDG3_POWER_CLASS_NAMES, "TDG3 class differs")
    return {
        "estimator": next(iter(estimators)),
        "measure": measure,
        "classification": classification,
        "point_counts": list(counts),
        "refinement_ratio": ratio_a,
        "observed_powers": list(observed),
        "lower_bounds": list(lower),
        "upper_bounds": list(upper),
        "zero_power_order_lower_bounds": list(zero_orders),
        "difference_order_lower_bound": difference_order,
        "conservative_native_surrogate_lower_bound": continuum_lower,
        "every_resolution_below_declared_amplitude_floor": every_floor,
        "replacement_admission_authorized": False,
    }


def _combined_class(left: str, right: str) -> str:
    zero = "contracts_to_zero_at_required_power_order"
    nonzero = {
        "nonzero_resolution_independent_within_arithmetic_enclosure",
        "converges_to_nonzero_native_power",
    }
    floor = "arithmetic_floor_dominated"
    unresolved = "unresolved_native_power_behavior"
    if left == right == zero:
        return "native_estimators_agree_zero_contracting"
    if left in nonzero and right in nonzero:
        return "native_estimators_agree_nonzero"
    if left == right == floor:
        return "native_estimators_agree_floor_dominated"
    if left == right == unresolved:
        return "native_estimators_agree_unresolved"
    return "native_estimators_disagree"


def diagnose_native_tail_ladder(
    proper_times_by_resolution: Sequence[object],
    values_by_resolution: Sequence[object],
    *,
    method: str,
    point_counts: Sequence[int],
) -> dict[str, Any]:
    counts = _method_counts(method, point_counts)
    _require(
        len(proper_times_by_resolution) == 3
        and len(values_by_resolution) == 3,
        "TDG3 requires three native histories",
    )
    observations: dict[str, tuple[NativeTailObservation, ...]] = {}
    classifications: dict[str, dict[str, Any]] = {}
    for estimator in TDG3_ESTIMATOR_NAMES:
        records = tuple(
            native_tail_observation(proper, values, estimator=estimator)
            for proper, values in zip(
                proper_times_by_resolution, values_by_resolution, strict=True
            )
        )
        observations[estimator] = records
        classifications[estimator] = {
            measure: classify_native_power_ladder(
                records,
                method=method,
                point_counts=counts,
                measure=measure,
            )
            for measure in TDG3_MEASURE_NAMES
        }
    combined = {
        measure: _combined_class(
            classifications[TDG3_ESTIMATOR_NAMES[0]][measure]["classification"],
            classifications[TDG3_ESTIMATOR_NAMES[1]][measure]["classification"],
        )
        for measure in TDG3_MEASURE_NAMES
    }
    _require(
        set(combined.values()) <= set(TDG3_COMBINED_CLASS_NAMES),
        "TDG3 combined class differs",
    )
    return {
        "method": method,
        "point_counts": list(counts),
        "observations": {
            estimator: [asdict(record) for record in records]
            for estimator, records in observations.items()
        },
        "classifications": classifications,
        "combined_classifications": combined,
        "common_grid_resampling_performed": False,
        "piecewise_linear_or_quadrature_surrogate_is_not_a_continuum_enclosure": True,
        "classification_is_not_a_replacement_admission": True,
    }


def _synthetic_times(jitter: float = 0.0) -> np.ndarray:
    base = np.linspace(0.0, 63.0 / 16.0, TDG3_SAMPLE_COUNT, dtype=np.float64)
    if jitter == 0.0:
        return base
    phase = np.linspace(0.0, 1.0, TDG3_SAMPLE_COUNT, dtype=np.float64)
    perturbation = jitter * np.sin(2.0 * pi * phase) * np.sin(pi * phase)
    candidate = base + perturbation
    _require(np.all(np.diff(candidate) > 0.0), "TDG3 synthetic grid folded")
    return candidate


def synthetic_native_tail_controls() -> dict[str, Any]:
    """Exercise the frozen native-grid discriminator without campaign arrays."""

    base = _synthetic_times()
    phase = (base - base[0]) / (base[-1] - base[0])
    carrier = np.sin(2.0 * pi * 29.0 * phase)
    times = (base.copy(), base.copy(), base.copy())
    controls: dict[str, tuple[tuple[np.ndarray, ...], tuple[np.ndarray, ...]]] = {
        "exact_zero": (times, tuple(np.zeros_like(base) for _ in range(3))),
        "subfloor_nonzero": (times, tuple(1.0e-40 * carrier for _ in range(3))),
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
            tuple(1.0e-3 * amplitude * carrier for amplitude in (1.25, 1.0625, 1.015625)),
        ),
        "unresolved_nonmonotone": (
            times,
            tuple(1.0e-3 * amplitude * carrier for amplitude in (1.0, 0.5, 0.75)),
        ),
    }
    jittered_times = tuple(_synthetic_times(value) for value in (0.02, 0.01, 0.005))
    jittered_values = tuple(
        1.0e-8
        * np.sin(
            2.0
            * pi
            * 29.0
            * (record - record[0])
            / (record[-1] - record[0])
        )
        for record in jittered_times
    )
    controls["common_grid_interpolation_owner"] = (jittered_times, jittered_values)
    _require(
        tuple(controls) == TDG3_SYNTHETIC_CONTROL_NAMES,
        "TDG3 synthetic control order differs",
    )
    expected = {
        "exact_zero": "native_estimators_agree_floor_dominated",
        "subfloor_nonzero": "native_estimators_agree_floor_dominated",
        "power_contracts_to_zero": "native_estimators_agree_zero_contracting",
        "nonzero_resolution_independent": "native_estimators_agree_nonzero",
        "nonzero_convergent": "native_estimators_agree_nonzero",
        "unresolved_nonmonotone": "native_estimators_agree_unresolved",
        "common_grid_interpolation_owner": "native_estimators_agree_nonzero",
    }
    records: dict[str, Any] = {}
    all_matched = True
    for name, (control_times, control_values) in controls.items():
        diagnosis = diagnose_native_tail_ladder(
            control_times,
            control_values,
            method="RK4",
            point_counts=TDG3_METHOD_POINT_COUNTS["RK4"],
        )
        observed = diagnosis["combined_classifications"]
        matched = all(value == expected[name] for value in observed.values())
        all_matched = all_matched and matched
        records[name] = {
            "expected_combined_classification": expected[name],
            "observed_combined_classifications": observed,
            "classification_matched": matched,
            "diagnosis": diagnosis,
        }
    interpolation_control = diagnose_absolute_tail_ladder(
        jittered_times,
        jittered_values,
        method="RK4",
        point_counts=TDG3_METHOD_POINT_COUNTS["RK4"],
    )
    tdg2_classes = {
        measure: interpolation_control["classifications"][measure]["classification"]
        for measure in ("field_tail_power", "derivative_tail_power")
    }
    tdg2_identifies_interpolation_owner = all(
        value == "instrument_floor_or_interpolation_enclosure_dominated"
        for value in tdg2_classes.values()
    )
    return {
        "controls": records,
        "all_control_classes_separated": (
            all_matched and tdg2_identifies_interpolation_owner
        ),
        "common_grid_interpolation_owner_control": {
            "TDG2_classifications": tdg2_classes,
            "TDG2_identifies_interpolation_owner": tdg2_identifies_interpolation_owner,
            "TDG3_combined_classifications": records[
                "common_grid_interpolation_owner"
            ]["observed_combined_classifications"],
            "TDG3_native_estimators_agree_nonzero": all(
                value == "native_estimators_agree_nonzero"
                for value in records["common_grid_interpolation_owner"][
                    "observed_combined_classifications"
                ].values()
            ),
        },
        "actual_terminal_histories_consumed": False,
        "common_grid_resampling_performed": False,
        "continuum_surrogate_enclosure_claimed": False,
        "replacement_temporal_admission_defined": False,
    }
