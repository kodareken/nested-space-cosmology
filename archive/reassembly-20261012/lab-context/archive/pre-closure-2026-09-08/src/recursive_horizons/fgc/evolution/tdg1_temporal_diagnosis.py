"""Prospective diagnostics for PROTO13's terminal temporal-spectrum stop.

The inherited PROTO5 temporal gate applies ``nested_spectral_admission`` to
*normalized* temporal tail fractions at three spatial resolutions.  TDG1 does
not change that historical rule.  It supplies two narrower instruments:

* synthetic controls that ask whether the rule is even a valid general test
  of spatial convergence for a nonzero temporal continuum signal; and
* outcome-neutral descriptions of the immutable 64-sample terminal histories
  that separate absolute signal scale, taper/detrending sensitivity,
  resampling sensitivity, interpolation debit, and raw history-difference
  convergence.

No routine in this module advances a state, changes the historical ``1/4``
threshold, authorizes a replacement gate, or classifies FGC-QR physics.
"""

from __future__ import annotations

from dataclasses import asdict
from math import isfinite, log
from numbers import Real
from typing import Any, Mapping, Sequence

import numpy as np

from .health_monitor import (
    CONSTRAINT_NORMALIZATION_FLOOR,
    compact_vacuum_buffer_window,
)
from .proto4_admission import (
    PROTO4_SPECTRAL_FIELD_ORDER,
    SpectralPowerBudget,
    SpectralThresholds,
    windowed_spectral_power_budget,
)
from .proto5_runtime import proto5_temporal_spectral_admission


TDG1_SAMPLE_COUNT = 64
TDG1_TRACER_COUNT = 48
TDG1_FIELD_NAMES = PROTO4_SPECTRAL_FIELD_ORDER
TDG1_CUTOFF = 16.0
TDG1_TAPER_FRACTION = 1.0 / 8.0
TDG1_MAXIMUM_NESTED_TAIL_RATIO = 1.0 / 4.0
TDG1_MINIMUM_HISTORY_DIFFERENCE_ORDER = 3.0 / 2.0
TDG1_ABSOLUTE_AMPLITUDE_FLOOR = CONSTRAINT_NORMALIZATION_FLOOR
TDG1_DOWNSAMPLE_COUNTS = (16, 32, 64)
TDG1_SYNTHETIC_CONTROL_NAMES = (
    "zero",
    "constant",
    "affine",
    "low_sine",
    "smooth_pulse",
    "amplitude_scaled_low_sine",
)
TDG1_VIEW_NAMES = (
    "inherited_compact",
    "rectangular",
    "mean_centered_compact",
    "affine_detrended_compact",
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


def _serial_budget(value: SpectralPowerBudget) -> dict[str, Any]:
    return asdict(value)


def _uniform_resample(
    proper_times: np.ndarray,
    values: np.ndarray,
    *,
    sample_count: int = TDG1_SAMPLE_COUNT,
) -> tuple[np.ndarray, np.ndarray, float]:
    _require(proper_times.shape == values.shape, "time/value shapes differ")
    _require(
        proper_times.size >= 16 and np.all(np.diff(proper_times) > 0.0),
        "proper times must increase strictly",
    )
    _require(
        isinstance(sample_count, int) and not isinstance(sample_count, bool),
        "sample_count must be an integer",
    )
    _require(sample_count >= 16, "sample_count must be at least 16")
    uniform = np.linspace(
        proper_times[0], proper_times[-1], sample_count, dtype=np.float64
    )
    resampled = np.interp(uniform, proper_times, values)
    reverse = np.interp(proper_times, uniform, resampled)
    debit = float(np.max(np.abs(reverse - values), initial=0.0))
    return uniform, resampled, debit


def _compact_window(coordinate: np.ndarray) -> np.ndarray:
    width = TDG1_TAPER_FRACTION * float(coordinate[-1] - coordinate[0])
    return compact_vacuum_buffer_window(
        coordinate,
        support_minimum=float(coordinate[0]),
        support_maximum=float(coordinate[-1]),
        taper_width=width,
    )


def _affine_residual(coordinate: np.ndarray, values: np.ndarray) -> np.ndarray:
    centered_coordinate = coordinate - float(np.mean(coordinate))
    centered_values = values - float(np.mean(values))
    denominator = float(np.dot(centered_coordinate, centered_coordinate))
    _require(denominator > 0.0, "affine detrending coordinate is degenerate")
    slope = float(np.dot(centered_coordinate, centered_values)) / denominator
    return centered_values - slope * centered_coordinate


def _failure_reasons(
    budget: SpectralPowerBudget,
    thresholds: SpectralThresholds,
) -> tuple[str, ...]:
    reasons: list[str] = []
    if (
        budget.top_band_field_power_fraction
        >= thresholds.maximum_top_band_field_power_fraction
    ):
        reasons.append("field_tail")
    if (
        budget.top_band_derivative_weighted_power_fraction
        >= thresholds.maximum_top_band_derivative_power_fraction
    ):
        reasons.append("derivative_tail")
    if budget.rms_scale_over_cutoff >= thresholds.maximum_rms_scale_over_cutoff:
        reasons.append("RMS_scale")
    return tuple(reasons)


def _power_bin(total_power: float) -> str:
    if total_power == 0.0:
        return "exact_zero"
    finfo = np.finfo(np.float64)
    if total_power < finfo.tiny:
        return "subnormal_power"
    if total_power < finfo.eps**2:
        return "below_epsilon_squared"
    if total_power < finfo.eps:
        return "between_epsilon_squared_and_epsilon"
    return "at_or_above_epsilon"


def temporal_signal_views(
    proper_times: object,
    values: object,
    *,
    cutoff: Real = TDG1_CUTOFF,
    thresholds: SpectralThresholds | None = None,
) -> dict[str, Any]:
    """Return fixed diagnostic views of one immutable temporal signal.

    The views are diagnostic only.  A favorable view cannot replace the
    inherited admission rule or authorize a trajectory.
    """

    proper = _array("proper_times", proper_times, ndim=1)
    data = _array("values", values, ndim=1)
    _require(
        proper.shape == (TDG1_SAMPLE_COUNT,)
        and data.shape == (TDG1_SAMPLE_COUNT,),
        "TDG1 signals must contain exactly 64 samples",
    )
    scale = _finite("cutoff", cutoff, positive=True)
    limits = SpectralThresholds() if thresholds is None else thresholds
    if not isinstance(limits, SpectralThresholds):
        raise TypeError("thresholds must be SpectralThresholds")

    uniform, resampled, interpolation_debit = _uniform_resample(proper, data)
    compact = _compact_window(uniform)
    centered = resampled - float(np.mean(resampled))
    detrended = _affine_residual(uniform, resampled)
    view_inputs = {
        "inherited_compact": (resampled, compact),
        "rectangular": (resampled, None),
        "mean_centered_compact": (centered, compact),
        "affine_detrended_compact": (detrended, compact),
    }
    _require(tuple(view_inputs) == TDG1_VIEW_NAMES, "TDG1 view order differs")
    views: dict[str, dict[str, Any]] = {}
    for name, (signal, window) in view_inputs.items():
        budget = windowed_spectral_power_budget(
            signal,
            uniform,
            cutoff=scale,
            window=window,
            thresholds=limits,
        )
        weighted = signal if window is None else signal * window
        transformed = np.fft.rfft(weighted)
        peak = float(np.max(np.abs(transformed), initial=0.0))
        views[name] = {
            "budget": _serial_budget(budget),
            "failure_reasons": list(_failure_reasons(budget, limits)),
            "windowed_FFT_peak_amplitude": peak,
            "below_declared_absolute_amplitude_floor": (
                peak <= limits.absolute_last_bin_amplitude_floor
            ),
            "top_band_field_power": (
                budget.total_field_power
                * budget.top_band_field_power_fraction
            ),
            "top_band_derivative_weighted_power": (
                budget.total_derivative_weighted_power
                * budget.top_band_derivative_weighted_power_fraction
            ),
        }

    downsampled: dict[str, dict[str, Any]] = {}
    for count in TDG1_DOWNSAMPLE_COUNTS:
        coordinate = np.linspace(uniform[0], uniform[-1], count, dtype=np.float64)
        signal = np.interp(coordinate, uniform, resampled)
        budget = windowed_spectral_power_budget(
            signal,
            coordinate,
            cutoff=scale,
            window=_compact_window(coordinate),
            thresholds=limits,
        )
        downsampled[str(count)] = {
            "budget": _serial_budget(budget),
            "failure_reasons": list(_failure_reasons(budget, limits)),
        }

    inherited = views["inherited_compact"]
    return {
        "raw_maximum_absolute": float(np.max(np.abs(data), initial=0.0)),
        "raw_variation_infinity": float(
            np.max(np.abs(data - float(np.mean(data))), initial=0.0)
        ),
        "round_trip_interpolation_infinity": interpolation_debit,
        "inherited_total_power_bin": _power_bin(
            float(inherited["budget"]["total_field_power"])
        ),
        "views": views,
        "resampling_sensitivity_views": downsampled,
        "resampling_views_are_not_independent_temporal_convergence_evidence": True,
    }


def _all_ratios_equal(
    admission: object,
    expected: float,
) -> bool:
    nested = admission.nested_admission
    values = (
        *nested.field_power_tail_ratios.values(),
        *nested.derivative_power_tail_ratios.values(),
    )
    return all(
        ratio is not None and ratio == expected
        for record in values
        for ratio in record
    )


def synthetic_temporal_gate_controls() -> dict[str, Any]:
    """Exercise the inherited gate without consuming campaign histories."""

    x = np.linspace(0.0, 63.0 / 16.0, TDG1_SAMPLE_COUNT, dtype=np.float64)
    tracer_count = 2
    times = np.repeat(x[:, None], tracer_count, axis=1)
    normalized = (x - x[0]) / (x[-1] - x[0])
    signals: dict[str, np.ndarray] = {
        "zero": np.zeros_like(x),
        "constant": np.full_like(x, 1.0 / 8.0),
        "affine": 1.0 / 8.0 + normalized / 64.0,
        "low_sine": np.sin(2.0 * np.pi * normalized) / 8.0,
        "smooth_pulse": np.exp(-((normalized - 0.5) / 0.2) ** 2) / 8.0,
    }
    records: dict[str, Any] = {}
    counts = (2049, 4097, 8193)
    for name, signal in signals.items():
        history = np.broadcast_to(
            signal[:, None, None],
            (TDG1_SAMPLE_COUNT, tracer_count, len(TDG1_FIELD_NAMES)),
        ).copy()
        admission = proto5_temporal_spectral_admission(
            point_counts=counts,
            proper_times_by_resolution=(times, times, times),
            field_histories_by_resolution=(history, history, history),
        )
        records[name] = {
            "admission_passed": admission.admission_passed,
            "every_individual_budget_passed": (
                admission.nested_admission.every_individual_budget_passed
            ),
            "every_nested_tail_ratio_passed": (
                admission.nested_admission.every_nested_tail_ratio_passed
            ),
            "all_field_and_derivative_ratios_equal": (
                0.0 if name == "zero" else 1.0
            ),
            "ratio_identity_verified": _all_ratios_equal(
                admission, 0.0 if name == "zero" else 1.0
            ),
        }

    base = signals["low_sine"]
    scaled_histories = tuple(
        np.broadcast_to(
            (amplitude * base)[:, None, None],
            (TDG1_SAMPLE_COUNT, tracer_count, len(TDG1_FIELD_NAMES)),
        ).copy()
        for amplitude in (1.0, 1.0 / 4.0, 1.0 / 16.0)
    )
    scaled = proto5_temporal_spectral_admission(
        point_counts=counts,
        proper_times_by_resolution=(times, times, times),
        field_histories_by_resolution=scaled_histories,
    )
    records["amplitude_scaled_low_sine"] = {
        "amplitude_sequence": [1.0, 1.0 / 4.0, 1.0 / 16.0],
        "amplitude_error_contracts_by_four_per_resolution": True,
        "power_error_contracts_by_sixteen_per_resolution": True,
        "admission_passed": scaled.admission_passed,
        "every_individual_budget_passed": (
            scaled.nested_admission.every_individual_budget_passed
        ),
        "every_nested_tail_ratio_passed": (
            scaled.nested_admission.every_nested_tail_ratio_passed
        ),
        "all_field_and_derivative_ratios_equal": 1.0,
        "ratio_identity_verified": _all_ratios_equal(scaled, 1.0),
    }
    _require(
        tuple(records) == TDG1_SYNTHETIC_CONTROL_NAMES,
        "TDG1 synthetic control order differs",
    )
    zero = records["zero"]
    nonzero = [records[name] for name in TDG1_SYNTHETIC_CONTROL_NAMES[1:]]
    confirmed = (
        zero["admission_passed"] is True
        and zero["ratio_identity_verified"] is True
        and all(item["every_individual_budget_passed"] is True for item in nonzero)
        and all(item["every_nested_tail_ratio_passed"] is False for item in nonzero)
        and all(item["admission_passed"] is False for item in nonzero)
        and all(item["ratio_identity_verified"] is True for item in nonzero)
    )
    return {
        "controls": records,
        "structural_counterexample_confirmed": confirmed,
        "mathematical_consequence": (
            "normalized_raw_temporal_tail_fractions_are_not_a_general_"
            "spatial_convergence_observable_for_nonzero_continuum_histories"
            if confirmed
            else "counterexample_not_confirmed"
        ),
        "historical_threshold_unchanged": TDG1_MAXIMUM_NESTED_TAIL_RATIO,
    }


def _common_overlap_histories(
    proper_times_by_resolution: Sequence[np.ndarray],
    values_by_resolution: Sequence[np.ndarray],
) -> tuple[np.ndarray, tuple[np.ndarray, ...], tuple[float, ...]]:
    start = max(float(record[0]) for record in proper_times_by_resolution)
    stop = min(float(record[-1]) for record in proper_times_by_resolution)
    _require(stop > start, "resolution histories have no common proper-time span")
    common = np.linspace(start, stop, TDG1_SAMPLE_COUNT, dtype=np.float64)
    aligned: list[np.ndarray] = []
    debits: list[float] = []
    for proper, values in zip(
        proper_times_by_resolution, values_by_resolution, strict=True
    ):
        aligned_values = np.interp(common, proper, values)
        mask = (proper >= start) & (proper <= stop)
        _require(np.count_nonzero(mask) >= 16, "common proper-time overlap is too short")
        reverse = np.interp(proper[mask], common, aligned_values)
        debit = float(
            np.max(np.abs(reverse - values[mask]), initial=0.0)
        )
        aligned.append(aligned_values)
        debits.append(debit)
    return common, tuple(aligned), tuple(debits)


def temporal_history_difference_diagnosis(
    proper_times_by_resolution: Sequence[object],
    values_by_resolution: Sequence[object],
    *,
    point_counts: Sequence[int],
    minimum_order: Real = TDG1_MINIMUM_HISTORY_DIFFERENCE_ORDER,
) -> dict[str, Any]:
    """Describe raw three-resolution history differences on common proper time."""

    counts = tuple(point_counts)
    _require(
        len(counts) == 3
        and all(
            isinstance(value, int) and not isinstance(value, bool)
            for value in counts
        )
        and all(b > a for a, b in zip(counts, counts[1:])),
        "history-difference point counts differ",
    )
    ratios = tuple(
        (b - 1) / (a - 1) for a, b in zip(counts, counts[1:])
    )
    _require(
        ratios[0] == ratios[1] and ratios[0] > 1.0,
        "history-difference refinement ratios differ",
    )
    required = _finite("minimum_order", minimum_order, positive=True)
    times = tuple(_array("proper_times", value, ndim=1) for value in proper_times_by_resolution)
    values = tuple(_array("values", value, ndim=1) for value in values_by_resolution)
    _require(
        len(times) == 3
        and len(values) == 3
        and all(item.shape == (TDG1_SAMPLE_COUNT,) for item in (*times, *values))
        and all(np.all(np.diff(item) > 0.0) for item in times),
        "history-difference inputs must be three increasing 64-sample records",
    )
    common, aligned, debits = _common_overlap_histories(times, values)
    coarse_difference = aligned[0] - aligned[1]
    fine_difference = aligned[1] - aligned[2]
    coarse_linf = float(np.max(np.abs(coarse_difference), initial=0.0))
    fine_linf = float(np.max(np.abs(fine_difference), initial=0.0))
    coarse_l2 = float(np.sqrt(np.mean(coarse_difference**2)))
    fine_l2 = float(np.sqrt(np.mean(fine_difference**2)))

    if coarse_linf == 0.0 and fine_linf == 0.0:
        classification = "exactly_resolution_identical"
        order = None
        raw_order_passed = True
    elif fine_linf == 0.0 and coarse_linf > 0.0:
        classification = "fine_pair_exact_after_nonzero_coarse_difference"
        order = None
        raw_order_passed = True
    elif coarse_linf > 0.0 and fine_linf > 0.0:
        order = log(coarse_linf / fine_linf) / log(ratios[0])
        if fine_linf < coarse_linf and order >= required:
            classification = "raw_order_pass"
            raw_order_passed = True
        elif fine_linf < coarse_linf:
            classification = "raw_order_below_required"
            raw_order_passed = False
        else:
            classification = "raw_difference_nonmonotone"
            raw_order_passed = False
    else:
        classification = "raw_difference_nonmonotone"
        order = None
        raw_order_passed = False

    pair_debits = (debits[0] + debits[1], debits[1] + debits[2])
    interpolation_limited = (
        coarse_linf <= pair_debits[0] or fine_linf <= pair_debits[1]
    )
    window = _compact_window(common)
    coarse_budget = windowed_spectral_power_budget(
        coarse_difference, common, cutoff=TDG1_CUTOFF, window=window
    )
    fine_budget = windowed_spectral_power_budget(
        fine_difference, common, cutoff=TDG1_CUTOFF, window=window
    )
    return {
        "point_counts": list(counts),
        "refinement_ratio": ratios[0],
        "common_proper_time_start": float(common[0]),
        "common_proper_time_stop": float(common[-1]),
        "raw_difference_infinity": [coarse_linf, fine_linf],
        "raw_difference_L2": [coarse_l2, fine_l2],
        "raw_observed_order": order,
        "minimum_required_order": required,
        "raw_order_passed": raw_order_passed,
        "classification": classification,
        "round_trip_interpolation_debits": list(debits),
        "pairwise_interpolation_debits": list(pair_debits),
        "interpolation_limited": interpolation_limited,
        "coarse_difference_spectral_budget": _serial_budget(coarse_budget),
        "fine_difference_spectral_budget": _serial_budget(fine_budget),
        "diagnostic_is_not_a_replacement_admission": True,
    }


def _empty_field_counts() -> dict[str, int]:
    return {name: 0 for name in TDG1_FIELD_NAMES}


def _aggregate_resolution_signals(
    proper_times: np.ndarray,
    fields: np.ndarray,
) -> dict[str, Any]:
    raw_failed = 0
    subfloor_failed = 0
    above_floor_failed = 0
    persistent_all_views = 0
    rescued = {name: 0 for name in TDG1_VIEW_NAMES[1:]}
    downsample_changes = {"16": 0, "32": 0}
    failures_by_field = _empty_field_counts()
    subfloor_by_field = _empty_field_counts()
    power_bins = {
        "exact_zero": 0,
        "subnormal_power": 0,
        "below_epsilon_squared": 0,
        "between_epsilon_squared_and_epsilon": 0,
        "at_or_above_epsilon": 0,
    }
    maximum_interpolation = 0.0
    worst: list[dict[str, Any]] = []
    for tracer in range(TDG1_TRACER_COUNT):
        for field_index, field_name in enumerate(TDG1_FIELD_NAMES):
            record = temporal_signal_views(
                proper_times[:, tracer], fields[:, tracer, field_index]
            )
            maximum_interpolation = max(
                maximum_interpolation,
                float(record["round_trip_interpolation_infinity"]),
            )
            power_bins[record["inherited_total_power_bin"]] += 1
            inherited = record["views"]["inherited_compact"]
            inherited_pass = bool(inherited["budget"]["individual_admission_passed"])
            if not inherited_pass:
                raw_failed += 1
                failures_by_field[field_name] += 1
                if inherited["below_declared_absolute_amplitude_floor"]:
                    subfloor_failed += 1
                    subfloor_by_field[field_name] += 1
                else:
                    above_floor_failed += 1
                if all(
                    not record["views"][name]["budget"]["individual_admission_passed"]
                    for name in TDG1_VIEW_NAMES
                ):
                    persistent_all_views += 1
                for name in TDG1_VIEW_NAMES[1:]:
                    if record["views"][name]["budget"]["individual_admission_passed"]:
                        rescued[name] += 1
                for count in ("16", "32"):
                    if (
                        record["resampling_sensitivity_views"][count]["budget"]
                        ["individual_admission_passed"]
                        != inherited_pass
                    ):
                        downsample_changes[count] += 1
                worst.append(
                    {
                        "tracer": tracer,
                        "field": field_name,
                        "total_field_power": inherited["budget"]["total_field_power"],
                        "FFT_peak_amplitude": inherited[
                            "windowed_FFT_peak_amplitude"
                        ],
                        "below_declared_absolute_amplitude_floor": inherited[
                            "below_declared_absolute_amplitude_floor"
                        ],
                        "failure_reasons": inherited["failure_reasons"],
                        "field_tail_fraction": inherited["budget"]
                        ["top_band_field_power_fraction"],
                        "derivative_tail_fraction": inherited["budget"]
                        ["top_band_derivative_weighted_power_fraction"],
                        "RMS_scale_over_cutoff": inherited["budget"]
                        ["rms_scale_over_cutoff"],
                        "round_trip_interpolation_infinity": record[
                            "round_trip_interpolation_infinity"
                        ],
                    }
                )
    worst.sort(
        key=lambda value: (
            float(value["RMS_scale_over_cutoff"]),
            float(value["field_tail_fraction"]),
            -int(value["tracer"]),
            str(value["field"]),
        ),
        reverse=True,
    )
    return {
        "signal_count": TDG1_TRACER_COUNT * len(TDG1_FIELD_NAMES),
        "inherited_individual_failure_count": raw_failed,
        "failed_below_declared_absolute_amplitude_floor": subfloor_failed,
        "failed_above_declared_absolute_amplitude_floor": above_floor_failed,
        "failed_persisting_under_all_fixed_views": persistent_all_views,
        "failures_rescued_by_fixed_view": rescued,
        "downsampled_view_classification_changes": downsample_changes,
        "failures_by_field": failures_by_field,
        "subfloor_failures_by_field": subfloor_by_field,
        "total_power_bins": power_bins,
        "maximum_round_trip_interpolation_infinity": maximum_interpolation,
        "worst_twelve_inherited_failures": worst[:12],
    }


def diagnose_terminal_method_histories(
    *,
    method: str,
    point_counts: Sequence[int],
    proper_times_by_resolution: Sequence[object],
    field_histories_by_resolution: Sequence[object],
) -> dict[str, Any]:
    """Apply the prospectively frozen TDG1 diagnostics to one method ladder."""

    counts = tuple(point_counts)
    _require(method in ("RK4", "SSPRK3"), "TDG1 method differs")
    expected = (2049, 4097, 8193) if method == "RK4" else (4097, 8193, 16385)
    _require(counts == expected, "TDG1 method-owned ladder differs")
    times = tuple(_array("proper_times", value, ndim=2) for value in proper_times_by_resolution)
    fields = tuple(_array("field_histories", value, ndim=3) for value in field_histories_by_resolution)
    _require(
        len(times) == 3
        and len(fields) == 3
        and all(item.shape == (TDG1_SAMPLE_COUNT, TDG1_TRACER_COUNT) for item in times)
        and all(
            item.shape
            == (TDG1_SAMPLE_COUNT, TDG1_TRACER_COUNT, len(TDG1_FIELD_NAMES))
            for item in fields
        )
        and all(np.all(np.diff(item, axis=0) > 0.0) for item in times),
        "TDG1 method histories have incompatible shapes",
    )
    by_resolution = {
        str(count): _aggregate_resolution_signals(time, field)
        for count, time, field in zip(counts, times, fields, strict=True)
    }

    classifications: dict[str, int] = {
        "exactly_resolution_identical": 0,
        "fine_pair_exact_after_nonzero_coarse_difference": 0,
        "raw_order_pass": 0,
        "raw_order_below_required": 0,
        "raw_difference_nonmonotone": 0,
    }
    interpolation_limited_count = 0
    raw_orders: list[float] = []
    by_field = {
        name: {
            "raw_order_pass": 0,
            "raw_order_below_required": 0,
            "raw_difference_nonmonotone": 0,
            "interpolation_limited": 0,
        }
        for name in TDG1_FIELD_NAMES
    }
    worst: list[dict[str, Any]] = []
    for tracer in range(TDG1_TRACER_COUNT):
        for field_index, field_name in enumerate(TDG1_FIELD_NAMES):
            record = temporal_history_difference_diagnosis(
                [item[:, tracer] for item in times],
                [item[:, tracer, field_index] for item in fields],
                point_counts=counts,
            )
            classification = str(record["classification"])
            classifications[classification] += 1
            if record["interpolation_limited"]:
                interpolation_limited_count += 1
                by_field[field_name]["interpolation_limited"] += 1
            if classification in by_field[field_name]:
                by_field[field_name][classification] += 1
            if record["raw_observed_order"] is not None:
                raw_orders.append(float(record["raw_observed_order"]))
            if not record["raw_order_passed"]:
                worst.append(
                    {
                        "tracer": tracer,
                        "field": field_name,
                        "classification": classification,
                        "raw_difference_infinity": record[
                            "raw_difference_infinity"
                        ],
                        "raw_observed_order": record["raw_observed_order"],
                        "interpolation_limited": record[
                            "interpolation_limited"
                        ],
                        "pairwise_interpolation_debits": record[
                            "pairwise_interpolation_debits"
                        ],
                    }
                )
    worst.sort(
        key=lambda value: (
            float("inf")
            if value["raw_observed_order"] is None
            else -float(value["raw_observed_order"]),
            -float(value["raw_difference_infinity"][1]),
            -int(value["tracer"]),
            str(value["field"]),
        )
    )
    return {
        "method": method,
        "point_counts": list(counts),
        "sample_count": TDG1_SAMPLE_COUNT,
        "tracer_count": TDG1_TRACER_COUNT,
        "field_names": list(TDG1_FIELD_NAMES),
        "individual_signal_views_by_resolution": by_resolution,
        "history_difference_classification_counts": classifications,
        "history_difference_interpolation_limited_count": interpolation_limited_count,
        "history_difference_raw_order_summary": {
            "finite_count": len(raw_orders),
            "minimum": min(raw_orders) if raw_orders else None,
            "median": float(np.median(raw_orders)) if raw_orders else None,
            "maximum": max(raw_orders) if raw_orders else None,
        },
        "history_difference_counts_by_field": by_field,
        "worst_twelve_nonpassing_history_differences": worst[:12],
        "diagnostic_does_not_authorize_a_replacement_gate": True,
    }


def diagnose_terminal_temporal_histories(
    *,
    proper_times_by_method: Mapping[str, Sequence[object]],
    field_histories_by_method: Mapping[str, Sequence[object]],
) -> dict[str, Any]:
    """Compose the synthetic theorem-control and immutable history diagnosis."""

    _require(
        tuple(proper_times_by_method) == ("RK4", "SSPRK3")
        and tuple(field_histories_by_method) == ("RK4", "SSPRK3"),
        "TDG1 method order differs",
    )
    synthetic = synthetic_temporal_gate_controls()
    methods = {
        "RK4": diagnose_terminal_method_histories(
            method="RK4",
            point_counts=(2049, 4097, 8193),
            proper_times_by_resolution=proper_times_by_method["RK4"],
            field_histories_by_resolution=field_histories_by_method["RK4"],
        ),
        "SSPRK3": diagnose_terminal_method_histories(
            method="SSPRK3",
            point_counts=(4097, 8193, 16385),
            proper_times_by_resolution=proper_times_by_method["SSPRK3"],
            field_histories_by_resolution=field_histories_by_method["SSPRK3"],
        ),
    }
    subfloor = sum(
        record["failed_below_declared_absolute_amplitude_floor"]
        for method in methods.values()
        for record in method["individual_signal_views_by_resolution"].values()
    )
    above_floor = sum(
        record["failed_above_declared_absolute_amplitude_floor"]
        for method in methods.values()
        for record in method["individual_signal_views_by_resolution"].values()
    )
    return {
        "synthetic_controls": synthetic,
        "methods": methods,
        "aggregate": {
            "inherited_raw_tail_ratio_gate_structurally_invalid": synthetic[
                "structural_counterexample_confirmed"
            ],
            "subfloor_inherited_individual_failure_count": subfloor,
            "above_floor_inherited_individual_failure_count": above_floor,
            "individual_budget_failure_cause_fully_derived": False,
            "replacement_temporal_admission_defined": False,
            "new_trajectory_authorized": False,
        },
    }


__all__ = [
    "TDG1_ABSOLUTE_AMPLITUDE_FLOOR",
    "TDG1_CUTOFF",
    "TDG1_DOWNSAMPLE_COUNTS",
    "TDG1_FIELD_NAMES",
    "TDG1_MAXIMUM_NESTED_TAIL_RATIO",
    "TDG1_MINIMUM_HISTORY_DIFFERENCE_ORDER",
    "TDG1_SAMPLE_COUNT",
    "TDG1_SYNTHETIC_CONTROL_NAMES",
    "TDG1_TAPER_FRACTION",
    "TDG1_TRACER_COUNT",
    "TDG1_VIEW_NAMES",
    "diagnose_terminal_method_histories",
    "diagnose_terminal_temporal_histories",
    "synthetic_temporal_gate_controls",
    "temporal_history_difference_diagnosis",
    "temporal_signal_views",
]
