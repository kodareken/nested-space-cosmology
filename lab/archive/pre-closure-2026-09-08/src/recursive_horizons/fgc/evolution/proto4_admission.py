"""Outcome-neutral PROTO4 numerical-admission machinery.

This module succeeds the six HLT1 decisions that CAL0 showed were not valid
physical admission rules for the declared compact input.  It does *not*
rewrite HLT1 history.  Instead it supplies three independent contracts:

* an exact per-branch applicability map for the 26 historical stop IDs;
* weighted, resolution-aware spectra in proper time or proper distance; and
* common-event, three-resolution constraint convergence with deterministic
  term-sum normalization.

No routine in this module evolves a candidate branch or classifies a collapse
outcome.  Candidate-only stops activate only after a branch-owned definition
has been supplied.  In particular, the existing HYP2 definition owns FGC-QR
only and cannot be reused silently for SGB-L.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import ceil, isfinite, log, pi, sqrt
from numbers import Real
from typing import Mapping, Sequence

import numpy as np

from ..scoped_run_authorization import (
    PROTO4_CANDIDATE_BRANCH_STOP_IDS,
    PROTO4_REPLACED_STOP_IDS,
    PROTO4_UNIVERSAL_STOP_IDS,
)
from .health_monitor import (
    CONSTRAINT_NORMALIZATION_FLOOR,
    FIELD_ORDER,
    STOP_PRIORITY,
)


PROTO4_BRANCH_ORDER = ("GR-0", "SGB-L", "FGC-QR")
PROTO4_SPECTRAL_FIELD_ORDER = (
    "alpha_minus_1",
    "shift_over_r",
    "lambda_minus_1",
    "R_over_r_minus_1",
    "phi_over_Lambda",
    "chi_over_Lambda",
)
PROTO4_CONSTRAINT_ORDER = (
    "Hamiltonian",
    "radial_momentum",
    "gauge_t",
    "gauge_r",
    "reduction_alpha",
    "reduction_shift",
    "reduction_lambda",
    "reduction_R",
    "reduction_phi",
    "reduction_chi",
)

# These are evidence owners, not aliases.  HYP2, DOM4, and SRC1 all state an
# FGC-QR scope.  No corresponding SGB-L bundle exists yet, so SGB-L remains a
# typed missing premise rather than inheriting the FGC-QR numbers.
FGCQR_BRANCH_DEFINITION_OWNER = (
    "FGC-1-SRC1-NL1+FGC-1-DOM4-RUN1+FGC-1-HYP2-MD1"
)
BRANCH_DEFINITION_OWNERS: Mapping[str, str | None] = {
    "GR-0": None,
    "SGB-L": None,
    "FGC-QR": FGCQR_BRANCH_DEFINITION_OWNER,
}


def _finite(name: str, value: Real, *, positive: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a real scalar")
    answer = float(value)
    if not isfinite(answer):
        raise ValueError(f"{name} must be finite")
    if positive and answer <= 0.0:
        raise ValueError(f"{name} must be positive")
    return answer


def _array(name: str, value: object, *, ndim: int | None = None) -> np.ndarray:
    answer = np.asarray(value, dtype=np.float64)
    if ndim is not None and answer.ndim != ndim:
        raise ValueError(f"{name} must have rank {ndim}")
    if not np.all(np.isfinite(answer)):
        raise ValueError(f"{name} must contain only finite values")
    return answer


@dataclass(frozen=True, slots=True)
class BranchStopApplicability:
    """Exact PROTO4 interpretation of the historical HLT1 stop IDs."""

    branch: str
    universal_stop_ids: tuple[str, ...]
    candidate_branch_stop_ids: tuple[str, ...]
    replaced_stop_ids: tuple[str, ...]
    candidate_definition_owner: str | None
    candidate_definition_available: bool
    runtime_monitor_complete: bool
    missing_premise: str | None


def branch_stop_applicability(
    branch: str,
    *,
    candidate_definition_owner: str | None = None,
) -> BranchStopApplicability:
    """Resolve stop ownership without borrowing another branch's evidence.

    GR-0 uses the nine universal runtime stops and must not be rejected by any
    target-action stop.  Candidate branches require their exact evidence
    owner before the eleven candidate stops become active.  The six replaced
    decisions are never active under PROTO4; their successor admission is
    evaluated separately.
    """

    if branch not in PROTO4_BRANCH_ORDER:
        raise ValueError("unknown PROTO4 branch")
    partition = (
        *PROTO4_UNIVERSAL_STOP_IDS,
        *PROTO4_CANDIDATE_BRANCH_STOP_IDS,
        *PROTO4_REPLACED_STOP_IDS,
    )
    # Category order intentionally differs from HLT1's scientific priority.
    if (
        len(partition) != len(STOP_PRIORITY)
        or len(set(partition)) != len(partition)
        or set(partition) != set(STOP_PRIORITY)
    ):
        raise RuntimeError("PROTO4 no longer partitions every HLT1 stop exactly once")

    expected_owner = BRANCH_DEFINITION_OWNERS[branch]
    if branch == "GR-0":
        if candidate_definition_owner is not None:
            raise ValueError("GR-0 cannot accept a candidate-branch definition owner")
        return BranchStopApplicability(
            branch=branch,
            universal_stop_ids=tuple(PROTO4_UNIVERSAL_STOP_IDS),
            candidate_branch_stop_ids=(),
            replaced_stop_ids=tuple(PROTO4_REPLACED_STOP_IDS),
            candidate_definition_owner=None,
            candidate_definition_available=True,
            runtime_monitor_complete=True,
            missing_premise=None,
        )

    if candidate_definition_owner is not None and candidate_definition_owner != expected_owner:
        if expected_owner is None:
            raise ValueError(f"{branch} has no certified candidate-branch definition owner")
        raise ValueError(f"{branch} candidate-branch definition owner differs")
    available = expected_owner is not None and candidate_definition_owner == expected_owner
    return BranchStopApplicability(
        branch=branch,
        universal_stop_ids=tuple(PROTO4_UNIVERSAL_STOP_IDS),
        candidate_branch_stop_ids=(
            tuple(PROTO4_CANDIDATE_BRANCH_STOP_IDS) if available else ()
        ),
        replaced_stop_ids=tuple(PROTO4_REPLACED_STOP_IDS),
        candidate_definition_owner=candidate_definition_owner,
        candidate_definition_available=available,
        runtime_monitor_complete=available,
        missing_premise=(
            None if available else f"{branch}_candidate_branch_health_definition_missing"
        ),
    )


@dataclass(frozen=True, slots=True)
class ProperRadialProfile:
    """Six dimensionless reference deviations on a uniform proper grid."""

    proper_coordinates: np.ndarray
    deviations: np.ndarray
    round_trip_interpolation_infinity: np.ndarray

    def __post_init__(self) -> None:
        proper = _array("proper_coordinates", self.proper_coordinates, ndim=1)
        deviations = _array("deviations", self.deviations, ndim=2)
        error = _array(
            "round_trip_interpolation_infinity",
            self.round_trip_interpolation_infinity,
            ndim=1,
        )
        if deviations.shape != (proper.size, len(PROTO4_SPECTRAL_FIELD_ORDER)):
            raise ValueError("proper-radial deviations have the wrong shape")
        if error.shape != (len(PROTO4_SPECTRAL_FIELD_ORDER),) or np.any(error < 0.0):
            raise ValueError("interpolation diagnostics have the wrong shape")
        spacing = np.diff(proper)
        if proper.size < 16 or np.any(spacing <= 0.0):
            raise ValueError("proper grid must contain at least 16 increasing points")
        tolerance = 4096.0 * np.finfo(np.float64).eps * max(1.0, abs(spacing[0]))
        if not np.allclose(spacing, spacing[0], rtol=0.0, atol=tolerance):
            raise ValueError("proper grid must be uniform")
        proper = np.ascontiguousarray(proper).copy()
        deviations = np.ascontiguousarray(deviations).copy()
        error = np.ascontiguousarray(error).copy()
        for array in (proper, deviations, error):
            array.setflags(write=False)
        object.__setattr__(self, "proper_coordinates", proper)
        object.__setattr__(self, "deviations", deviations)
        object.__setattr__(self, "round_trip_interpolation_infinity", error)


def minkowski_reference_deviations(
    u: object,
    q: object,
    radii: object,
    *,
    cutoff: Real,
) -> np.ndarray:
    """Return the six frozen PROTO4 dimensionless spectral fields.

    At a regular centre, ``shift/r`` and ``R/r`` are evaluated from their
    analytic first radial derivatives.  Numerical epsilon division is never
    used as a substitute for the centre limit.
    """

    fields = _array("u", u, ndim=2)
    radial_derivatives = _array("q", q, ndim=2)
    radius = _array("radii", radii, ndim=1)
    scale = _finite("cutoff", cutoff, positive=True)
    if fields.shape != radial_derivatives.shape or fields.shape != (
        radius.size,
        len(FIELD_ORDER),
    ):
        raise ValueError("u, q, and radii must share the canonical six-field grid")
    if np.any(np.diff(radius) <= 0.0) or radius[0] < 0.0:
        raise ValueError("radii must be nonnegative and strictly increasing")

    positive = radius > 0.0
    shift_over_r = np.empty_like(radius)
    areal_over_r = np.empty_like(radius)
    shift_over_r[positive] = fields[positive, 1] / radius[positive]
    areal_over_r[positive] = fields[positive, 3] / radius[positive]
    if not positive[0]:
        if radius[0] != 0.0:
            raise ValueError("only an exact r=0 point may use analytic centre limits")
        shift_over_r[0] = radial_derivatives[0, 1]
        areal_over_r[0] = radial_derivatives[0, 3]
    answer = np.column_stack(
        (
            fields[:, 0] - 1.0,
            shift_over_r,
            fields[:, 2] - 1.0,
            areal_over_r - 1.0,
            fields[:, 4] / scale,
            fields[:, 5] / scale,
        )
    )
    if not np.all(np.isfinite(answer)):
        raise ValueError("reference deviations became nonfinite")
    return answer


def proper_radial_profile(
    u: object,
    q: object,
    radii: object,
    *,
    cutoff: Real,
    measurement_radius_maximum: Real,
) -> ProperRadialProfile:
    """Resample the frozen reference deviations onto proper radial distance.

    Proper distance is accumulated by a deterministic trapezoidal integral
    of the positive radial metric ``lambda``.  Linear interpolation is used
    only to put the data on a uniform Fourier grid.  The round-trip mismatch
    is returned explicitly; it is a diagnostic/error input, not a proof of a
    continuum interpolation bound.
    """

    fields = _array("u", u, ndim=2)
    radial_derivatives = _array("q", q, ndim=2)
    radius = _array("radii", radii, ndim=1)
    maximum = _finite(
        "measurement_radius_maximum", measurement_radius_maximum, positive=True
    )
    if fields.shape != radial_derivatives.shape or fields.shape[0] != radius.size:
        raise ValueError("proper-radial inputs have inconsistent shapes")
    if fields.shape[1] != len(FIELD_ORDER):
        raise ValueError("proper-radial fields are not in canonical order")
    if radius[0] != 0.0 or np.any(np.diff(radius) <= 0.0):
        raise ValueError("proper-radial coordinates require a regular r=0 grid")
    radial_metric = fields[:, 2]
    if np.any(radial_metric <= 0.0):
        raise ValueError("proper distance requires a positive radial metric")
    selected = np.flatnonzero(radius <= maximum)
    if selected.size < 16 or selected[-1] != selected.size - 1:
        raise ValueError("measurement interval must be a contiguous prefix with at least 16 points")
    last = int(selected[-1])
    r = radius[: last + 1]
    lam = radial_metric[: last + 1]
    proper = np.zeros_like(r)
    proper[1:] = np.cumsum(0.5 * (lam[:-1] + lam[1:]) * np.diff(r))
    if proper[-1] <= 0.0 or np.any(np.diff(proper) <= 0.0):
        raise ValueError("proper radial integral is not strictly increasing")

    deviations = minkowski_reference_deviations(
        fields[: last + 1],
        radial_derivatives[: last + 1],
        r,
        cutoff=cutoff,
    )
    uniform = np.linspace(proper[0], proper[-1], proper.size, dtype=np.float64)
    resampled = np.column_stack(
        [np.interp(uniform, proper, deviations[:, index]) for index in range(deviations.shape[1])]
    )
    reverse = np.column_stack(
        [np.interp(proper, uniform, resampled[:, index]) for index in range(deviations.shape[1])]
    )
    round_trip = np.max(np.abs(reverse - deviations), axis=0)
    return ProperRadialProfile(uniform, resampled, round_trip)


def causal_past_window(proper_times: object, *, taper_fraction: Real = 1.0 / 8.0) -> np.ndarray:
    """Return a two-edge smooth window built only from an available past history."""

    times = _array("proper_times", proper_times, ndim=1)
    fraction = _finite("taper_fraction", taper_fraction, positive=True)
    if times.size < 64 or np.any(np.diff(times) <= 0.0):
        raise ValueError("causal temporal spectrum requires 64 increasing past samples")
    if fraction >= 0.5:
        raise ValueError("temporal taper fraction must be below one half")
    span = times[-1] - times[0]
    width = fraction * span
    from .health_monitor import compact_vacuum_buffer_window

    return compact_vacuum_buffer_window(
        times,
        support_minimum=times[0],
        support_maximum=times[-1],
        taper_width=width,
    )


@dataclass(frozen=True, slots=True)
class SpectralThresholds:
    unresolved_top_fraction: float = 1.0 / 8.0
    maximum_top_band_field_power_fraction: float = 2.0**-20
    maximum_top_band_derivative_power_fraction: float = 2.0**-10
    maximum_rms_scale_over_cutoff: float = 1.0 / 4.0
    maximum_nested_tail_ratio: float = 1.0 / 4.0
    relative_last_bin_amplitude_floor: float = 2.0**-40
    absolute_last_bin_amplitude_floor: float = CONSTRAINT_NORMALIZATION_FLOOR

    def __post_init__(self) -> None:
        for name in (
            "unresolved_top_fraction",
            "maximum_top_band_field_power_fraction",
            "maximum_top_band_derivative_power_fraction",
            "maximum_rms_scale_over_cutoff",
            "maximum_nested_tail_ratio",
            "relative_last_bin_amplitude_floor",
        ):
            value = _finite(name, getattr(self, name), positive=True)
            if value >= 1.0:
                raise ValueError(f"{name} must lie below one")
            object.__setattr__(self, name, value)
        absolute = _finite(
            "absolute_last_bin_amplitude_floor",
            self.absolute_last_bin_amplitude_floor,
        )
        if absolute < 0.0:
            raise ValueError("absolute last-bin floor must be nonnegative")
        object.__setattr__(self, "absolute_last_bin_amplitude_floor", absolute)


@dataclass(frozen=True, slots=True)
class SpectralPowerBudget:
    sample_count: int
    top_band_first_bin: int
    nyquist_bin: int
    total_field_power: float
    top_band_field_power_fraction: float
    total_derivative_weighted_power: float
    top_band_derivative_weighted_power_fraction: float
    rms_angular_scale: float
    rms_scale_over_cutoff: float
    last_occupied_bin: int
    last_bin_alias_diagnostic: bool
    individual_admission_passed: bool

    def __post_init__(self) -> None:
        for name in (
            "sample_count",
            "top_band_first_bin",
            "nyquist_bin",
            "last_occupied_bin",
        ):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value < 0:
                raise ValueError(f"{name} must be a nonnegative integer")
        for name in (
            "total_field_power",
            "top_band_field_power_fraction",
            "total_derivative_weighted_power",
            "top_band_derivative_weighted_power_fraction",
            "rms_angular_scale",
            "rms_scale_over_cutoff",
        ):
            value = _finite(name, getattr(self, name))
            if value < 0.0:
                raise ValueError(f"{name} must be nonnegative")
            object.__setattr__(self, name, value)
        for name in ("last_bin_alias_diagnostic", "individual_admission_passed"):
            if not isinstance(getattr(self, name), bool):
                raise TypeError(f"{name} must be bool")
        if self.sample_count < 16:
            raise ValueError("sample_count must be at least 16")
        if self.nyquist_bin != self.sample_count // 2:
            raise ValueError("nyquist_bin does not match sample_count")
        if not 0 <= self.top_band_first_bin <= self.nyquist_bin:
            raise ValueError("top-band index lies outside the resolved spectrum")
        if not 0 <= self.last_occupied_bin <= self.nyquist_bin:
            raise ValueError("last occupied bin lies outside the resolved spectrum")
        if self.top_band_field_power_fraction > 1.0:
            raise ValueError("field-power tail fraction cannot exceed one")
        if self.top_band_derivative_weighted_power_fraction > 1.0:
            raise ValueError("derivative-power tail fraction cannot exceed one")


def windowed_spectral_power_budget(
    samples: object,
    proper_coordinates: object,
    *,
    cutoff: Real,
    window: object | None = None,
    thresholds: SpectralThresholds | None = None,
) -> SpectralPowerBudget:
    """Measure weighted field/derivative tails and RMS physical scale.

    One-sided real-FFT multiplicities are included, so the power sums obey the
    real Parseval convention.  The largest occupied-bin Boolean is retained
    solely as a diagnostic and never participates in ``individual_admission``.
    """

    data = _array("samples", samples, ndim=1)
    coordinate = _array("proper_coordinates", proper_coordinates, ndim=1)
    scale = _finite("cutoff", cutoff, positive=True)
    limits = SpectralThresholds() if thresholds is None else thresholds
    if not isinstance(limits, SpectralThresholds):
        raise TypeError("thresholds must be SpectralThresholds")
    if data.size != coordinate.size or data.size < 16:
        raise ValueError("spectral budget requires at least 16 matched samples")
    spacings = np.diff(coordinate)
    if np.any(spacings <= 0.0):
        raise ValueError("proper coordinates must be strictly increasing")
    spacing = float(spacings[0])
    tolerance = 4096.0 * np.finfo(np.float64).eps * max(1.0, abs(spacing))
    if not np.allclose(spacings, spacing, rtol=0.0, atol=tolerance):
        raise ValueError("spectral budget requires a uniform proper grid")
    if window is None:
        weights = np.ones_like(data)
    else:
        weights = _array("window", window, ndim=1)
        if weights.shape != data.shape or np.any(weights < 0.0) or np.any(weights > 1.0):
            raise ValueError("spectral window must match samples and lie in [0,1]")

    transformed = np.fft.rfft(data * weights)
    amplitudes = np.abs(transformed)
    one_sided = np.full(amplitudes.size, 2.0, dtype=np.float64)
    one_sided[0] = 1.0
    if data.size % 2 == 0:
        one_sided[-1] = 1.0
    power = one_sided * amplitudes**2 / float(data.size * data.size)
    angular = 2.0 * pi * np.fft.rfftfreq(data.size, d=spacing)
    derivative_power = angular**2 * power
    total = float(np.sum(power))
    derivative_total = float(np.sum(derivative_power))
    nyquist_bin = amplitudes.size - 1
    first_top = int(ceil((1.0 - limits.unresolved_top_fraction) * nyquist_bin))
    top_total = float(np.sum(power[first_top:]))
    top_derivative = float(np.sum(derivative_power[first_top:]))
    field_fraction = top_total / total if total > 0.0 else 0.0
    derivative_fraction = (
        top_derivative / derivative_total if derivative_total > 0.0 else 0.0
    )
    rms = sqrt(derivative_total / total) if total > 0.0 else 0.0

    peak = float(np.max(amplitudes, initial=0.0))
    amplitude_floor = max(
        limits.absolute_last_bin_amplitude_floor,
        limits.relative_last_bin_amplitude_floor * peak,
    )
    occupied = np.flatnonzero(amplitudes > amplitude_floor)
    last_bin = int(occupied[-1]) if occupied.size else 0
    diagnostic = last_bin >= first_top
    admitted = (
        field_fraction < limits.maximum_top_band_field_power_fraction
        and derivative_fraction < limits.maximum_top_band_derivative_power_fraction
        and rms / scale < limits.maximum_rms_scale_over_cutoff
    )
    return SpectralPowerBudget(
        sample_count=data.size,
        top_band_first_bin=first_top,
        nyquist_bin=nyquist_bin,
        total_field_power=total,
        top_band_field_power_fraction=field_fraction,
        total_derivative_weighted_power=derivative_total,
        top_band_derivative_weighted_power_fraction=derivative_fraction,
        rms_angular_scale=rms,
        rms_scale_over_cutoff=rms / scale,
        last_occupied_bin=last_bin,
        last_bin_alias_diagnostic=diagnostic,
        individual_admission_passed=admitted,
    )


def spectral_field_budgets(
    samples: object,
    proper_coordinates: object,
    *,
    cutoff: Real,
    window: object | None = None,
    field_names: Sequence[str] = PROTO4_SPECTRAL_FIELD_ORDER,
    thresholds: SpectralThresholds | None = None,
) -> dict[str, SpectralPowerBudget]:
    """Evaluate one matched multi-field proper spectrum."""

    data = _array("samples", samples, ndim=2)
    names = tuple(field_names)
    if len(set(names)) != len(names) or data.shape[1] != len(names):
        raise ValueError("spectral field names must uniquely match the sample columns")
    return {
        name: windowed_spectral_power_budget(
            data[:, index],
            proper_coordinates,
            cutoff=cutoff,
            window=window,
            thresholds=thresholds,
        )
        for index, name in enumerate(names)
    }


def _tail_ratio(numerator: float, denominator: float) -> float | None:
    if denominator > 0.0:
        return numerator / denominator
    if numerator == 0.0:
        return 0.0
    return None


@dataclass(frozen=True, slots=True)
class NestedSpectralAdmission:
    grid_point_counts: tuple[int, ...]
    field_power_tail_ratios: Mapping[str, tuple[float | None, ...]]
    derivative_power_tail_ratios: Mapping[str, tuple[float | None, ...]]
    every_individual_budget_passed: bool
    every_nested_tail_ratio_passed: bool
    admission_passed: bool


def nested_spectral_admission(
    grid_point_counts: Sequence[int],
    budgets: Sequence[Mapping[str, SpectralPowerBudget]],
    *,
    thresholds: SpectralThresholds | None = None,
) -> NestedSpectralAdmission:
    """Require weighted-tail decay on at least three nested resolutions."""

    limits = SpectralThresholds() if thresholds is None else thresholds
    if not isinstance(limits, SpectralThresholds):
        raise TypeError("thresholds must be SpectralThresholds")
    counts = tuple(grid_point_counts)
    if len(counts) < 3 or len(counts) != len(budgets):
        raise ValueError("nested spectrum requires at least three matched resolutions")
    if any(isinstance(n, bool) or not isinstance(n, int) or n < 16 for n in counts):
        raise ValueError("grid point counts must be integers of at least 16")
    if any(later <= earlier for earlier, later in zip(counts, counts[1:])):
        raise ValueError("grid point counts must increase strictly")
    names = tuple(budgets[0])
    if not names or len(set(names)) != len(names):
        raise ValueError("spectral budgets require unique field names")
    if any(tuple(record) != names for record in budgets):
        raise ValueError("all nested spectra must retain one field order")
    if any(not isinstance(value, SpectralPowerBudget) for record in budgets for value in record.values()):
        raise TypeError("nested spectra must contain SpectralPowerBudget values")

    field_ratios: dict[str, tuple[float | None, ...]] = {}
    derivative_ratios: dict[str, tuple[float | None, ...]] = {}
    every_ratio = True
    every_individual = True
    for name in names:
        field_values = [record[name].top_band_field_power_fraction for record in budgets]
        derivative_values = [
            record[name].top_band_derivative_weighted_power_fraction for record in budgets
        ]
        field = tuple(_tail_ratio(b, a) for a, b in zip(field_values, field_values[1:]))
        derivative = tuple(
            _tail_ratio(b, a) for a, b in zip(derivative_values, derivative_values[1:])
        )
        field_ratios[name] = field
        derivative_ratios[name] = derivative
        every_ratio = every_ratio and all(
            value is not None and value < limits.maximum_nested_tail_ratio
            for value in (*field, *derivative)
        )
        every_individual = every_individual and all(
            record[name].individual_admission_passed for record in budgets
        )
    return NestedSpectralAdmission(
        grid_point_counts=counts,
        field_power_tail_ratios=field_ratios,
        derivative_power_tail_ratios=derivative_ratios,
        every_individual_budget_passed=every_individual,
        every_nested_tail_ratio_passed=every_ratio,
        admission_passed=every_individual and every_ratio,
    )


@dataclass(frozen=True, slots=True)
class ConstraintNorms:
    sample_count: int
    component_names: tuple[str, ...]
    component_infinity: np.ndarray
    global_infinity: float
    minimum_denominator: float
    maximum_denominator: float

    def __post_init__(self) -> None:
        if (
            isinstance(self.sample_count, bool)
            or not isinstance(self.sample_count, int)
            or self.sample_count < 1
        ):
            raise ValueError("constraint sample_count must be a positive integer")
        names = tuple(self.component_names)
        if not names or len(set(names)) != len(names):
            raise ValueError("constraint names must be nonempty and unique")
        components = _array("component_infinity", self.component_infinity, ndim=1)
        if components.shape != (len(names),) or np.any(components < 0.0):
            raise ValueError("constraint component norms have the wrong shape")
        global_value = _finite("global_infinity", self.global_infinity)
        minimum = _finite("minimum_denominator", self.minimum_denominator, positive=True)
        maximum = _finite("maximum_denominator", self.maximum_denominator, positive=True)
        if maximum < minimum:
            raise ValueError("constraint denominator extrema are reversed")
        tolerance = 128.0 * np.finfo(np.float64).eps * max(1.0, global_value)
        if abs(global_value - float(np.max(components, initial=0.0))) > tolerance:
            raise ValueError("global constraint norm does not match component maxima")
        components = np.ascontiguousarray(components).copy()
        components.setflags(write=False)
        object.__setattr__(self, "component_names", names)
        object.__setattr__(self, "component_infinity", components)
        object.__setattr__(self, "global_infinity", global_value)
        object.__setattr__(self, "minimum_denominator", minimum)
        object.__setattr__(self, "maximum_denominator", maximum)


def normalized_constraint_norms(
    residuals: object,
    unredefined_term_contributions: object,
    *,
    planck_mass: Real,
    length_unit: Real,
    component_names: Sequence[str] = PROTO4_CONSTRAINT_ORDER,
) -> ConstraintNorms:
    """Apply PROTO4's deterministic componentwise term-sum normalization."""

    residual = _array("residuals", residuals, ndim=2)
    terms = _array("unredefined_term_contributions", unredefined_term_contributions, ndim=3)
    mass = _finite("planck_mass", planck_mass, positive=True)
    length = _finite("length_unit", length_unit, positive=True)
    names = tuple(component_names)
    if residual.shape[:2] != terms.shape[:2] or residual.shape[1] != len(names):
        raise ValueError("constraint residuals, terms, and names must agree")
    if terms.shape[2] < 1 or len(set(names)) != len(names):
        raise ValueError("each uniquely named constraint requires at least one unredefined term")
    physical_floor = (mass * mass) / (length * length)
    denominator = np.maximum(
        np.sum(np.abs(terms), axis=2),
        max(physical_floor, CONSTRAINT_NORMALIZATION_FLOOR),
    )
    if np.any(denominator <= 0.0) or not np.all(np.isfinite(denominator)):
        raise ValueError("constraint normalization denominator must be finite and positive")
    ratio = np.abs(residual) / denominator
    component = np.max(ratio, axis=0, initial=0.0)
    return ConstraintNorms(
        sample_count=residual.shape[0],
        component_names=names,
        component_infinity=component,
        global_infinity=float(np.max(component, initial=0.0)),
        minimum_denominator=float(np.min(denominator, initial=np.inf)),
        maximum_denominator=float(np.max(denominator, initial=0.0)),
    )


@dataclass(frozen=True, slots=True)
class ConstraintEventSample:
    point_count: int
    coordinate_time: float
    norms: ConstraintNorms
    affine_label: float | None = None
    affine_origin: float | None = None

    def __post_init__(self) -> None:
        if isinstance(self.point_count, bool) or not isinstance(self.point_count, int) or self.point_count < 3:
            raise ValueError("point_count must be an integer of at least three")
        object.__setattr__(self, "coordinate_time", _finite("coordinate_time", self.coordinate_time))
        if not isinstance(self.norms, ConstraintNorms):
            raise TypeError("norms must be ConstraintNorms")
        if self.point_count != self.norms.sample_count:
            raise ValueError("point_count must match the normalized constraint grid")
        if (self.affine_label is None) != (self.affine_origin is None):
            raise ValueError("affine label and origin must be supplied together")
        if self.affine_label is not None:
            object.__setattr__(self, "affine_label", _finite("affine_label", self.affine_label))
            object.__setattr__(self, "affine_origin", _finite("affine_origin", self.affine_origin))


@dataclass(frozen=True, slots=True)
class CommonEventConstraintAdmission:
    point_counts: tuple[int, ...]
    component_names: tuple[str, ...]
    event_alignment_contract: str
    adjacent_refinement_ratios: tuple[float, ...]
    component_observed_orders: Mapping[str, tuple[float | None, ...]]
    component_order_status: Mapping[str, tuple[str, ...]]
    minimum_finite_observed_order: float | None
    coarsest_guard_passed: bool
    finest_guard_passed: bool
    convergence_order_passed: bool
    common_event_alignment_passed: bool
    admission_passed: bool


def _order_pair(coarse: float, fine: float, refinement_ratio: float) -> tuple[float | None, str]:
    if coarse == 0.0 and fine == 0.0:
        return None, "exact_zero_pair"
    if coarse > 0.0 and fine == 0.0:
        return None, "resolved_to_exact_zero"
    if coarse == 0.0 and fine > 0.0:
        return None, "nonzero_reappeared_after_exact_zero"
    return log(coarse / fine) / log(refinement_ratio), "finite_order"


def _binary64_equal(values: Sequence[float]) -> bool:
    """Return exact IEEE-754 binary64 equality, including the sign of zero."""

    if not values:
        return False
    reference = np.float64(values[0]).tobytes()
    return all(np.float64(value).tobytes() == reference for value in values[1:])


def common_event_constraint_admission(
    samples: Sequence[ConstraintEventSample],
    *,
    event_alignment_contract: str = "GR0_coordinate_time",
    coarsest_guard_maximum: Real = 1.0 / 50.0,
    finest_guard_maximum: Real = 1.0 / 1000.0,
    minimum_observed_order: Real = 3.0 / 2.0,
) -> CommonEventConstraintAdmission:
    """Require small, convergent constraints at one predeclared common event."""

    records = tuple(samples)
    if len(records) < 3 or any(not isinstance(item, ConstraintEventSample) for item in records):
        raise ValueError("constraint admission requires at least three event samples")
    coarse_guard = _finite("coarsest_guard_maximum", coarsest_guard_maximum, positive=True)
    fine_guard = _finite("finest_guard_maximum", finest_guard_maximum, positive=True)
    order_minimum = _finite("minimum_observed_order", minimum_observed_order, positive=True)
    counts = tuple(item.point_count for item in records)
    if any(later <= earlier for earlier, later in zip(counts, counts[1:])):
        raise ValueError("constraint resolutions must increase strictly")
    names = records[0].norms.component_names
    if any(item.norms.component_names != names for item in records):
        raise ValueError("constraint component order differs across resolutions")
    if event_alignment_contract not in {
        "GR0_coordinate_time",
        "DEF1_coordinate_time_affine_label_and_origin",
    }:
        raise ValueError("unknown common-event alignment contract")
    times = tuple(item.coordinate_time for item in records)
    affine = tuple((item.affine_label, item.affine_origin) for item in records)
    if event_alignment_contract == "GR0_coordinate_time":
        aligned = _binary64_equal(times) and all(
            value == (None, None) for value in affine
        )
    else:
        aligned = (
            affine[0] != (None, None)
            and all(value == times[0] for value in times)
            and all(value == affine[0] for value in affine)
        )

    intervals = tuple(count - 1 for count in counts)
    refinement = tuple(
        float(fine / coarse) for coarse, fine in zip(intervals, intervals[1:])
    )
    if any(value <= 1.0 for value in refinement):
        raise ValueError("constraint refinement ratios must exceed one")

    orders: dict[str, tuple[float | None, ...]] = {}
    statuses: dict[str, tuple[str, ...]] = {}
    all_order_pass = True
    finite_orders: list[float] = []
    for index, name in enumerate(names):
        values = [float(item.norms.component_infinity[index]) for item in records]
        pairs = tuple(
            _order_pair(coarse, fine, ratio)
            for coarse, fine, ratio in zip(values, values[1:], refinement)
        )
        component_orders = tuple(value for value, _status in pairs)
        component_status = tuple(status for _value, status in pairs)
        orders[name] = component_orders
        statuses[name] = component_status
        for value, status in pairs:
            if status in {"exact_zero_pair", "resolved_to_exact_zero"}:
                continue
            if status != "finite_order" or value is None:
                all_order_pass = False
                continue
            finite_orders.append(value)
            if value < order_minimum:
                all_order_pass = False

    coarse_pass = records[0].norms.global_infinity < coarse_guard
    fine_pass = records[-1].norms.global_infinity < fine_guard
    passed = aligned and coarse_pass and fine_pass and all_order_pass
    return CommonEventConstraintAdmission(
        point_counts=counts,
        component_names=names,
        event_alignment_contract=event_alignment_contract,
        adjacent_refinement_ratios=refinement,
        component_observed_orders=orders,
        component_order_status=statuses,
        minimum_finite_observed_order=(min(finite_orders) if finite_orders else None),
        coarsest_guard_passed=coarse_pass,
        finest_guard_passed=fine_pass,
        convergence_order_passed=all_order_pass,
        common_event_alignment_passed=aligned,
        admission_passed=passed,
    )


def conservative_richardson_error(
    fine_values: object,
    coarse_values_restricted_to_fine_event: object,
    *,
    refinement_ratio: Real,
    assumed_order: Real,
    event_alignment: CommonEventConstraintAdmission,
) -> float:
    """Return a nonnegative Richardson estimate that may only be added.

    The caller must first express both arrays in the same observable units and
    align them at the same event/trajectory label.  This routine deliberately
    has no subtraction mode and does not turn an asymptotic estimate into a
    rigorous continuum theorem.
    """

    if not isinstance(event_alignment, CommonEventConstraintAdmission):
        raise TypeError("event_alignment must be a common-event admission record")
    if not event_alignment.common_event_alignment_passed:
        raise ValueError("Richardson inputs require a passed common-event alignment")
    fine = _array("fine_values", fine_values)
    coarse = _array(
        "coarse_values_restricted_to_fine_event",
        coarse_values_restricted_to_fine_event,
    )
    ratio = _finite("refinement_ratio", refinement_ratio, positive=True)
    order = _finite("assumed_order", assumed_order, positive=True)
    if fine.shape != coarse.shape or ratio <= 1.0:
        raise ValueError("Richardson inputs must be matched and refinement must exceed one")
    denominator = ratio**order - 1.0
    if denominator <= 0.0 or not isfinite(denominator):
        raise ValueError("Richardson denominator must be finite and positive")
    return float(np.max(np.abs(fine - coarse), initial=0.0) / denominator)


def conservative_observable_error_sum(
    components_in_common_observable_units: Mapping[str, Real],
) -> float:
    """Add nonnegative named error bounds; never cancel one against another."""

    if not components_in_common_observable_units or "constraint" not in components_in_common_observable_units:
        raise ValueError("observable error budget must include a constraint contribution")
    values = []
    for name, value in components_in_common_observable_units.items():
        number = _finite(f"error component {name}", value)
        if number < 0.0:
            raise ValueError("observable error components must be nonnegative")
        values.append(number)
    return float(sum(values))
