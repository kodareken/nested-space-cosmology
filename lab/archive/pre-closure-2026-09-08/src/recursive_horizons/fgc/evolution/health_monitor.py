"""Executable, transactional health monitors for the scoped FGC run.

The monitor is intentionally stricter than a diagnostics logger: a trial stage
is either checked against every declared premise and accepted, or it is
rejected while the last accepted state and the first failure remain intact.
The scale checks are classical stop proxies.  They are never promoted here to
a Wilsonian EFT remainder estimate.
"""

from __future__ import annotations

from dataclasses import dataclass, fields
from math import ceil, exp, isfinite, pi, sqrt
from numbers import Real
from typing import Mapping, Sequence

import numpy as np


FIELD_ORDER = ("alpha", "shift", "lambda", "R", "phi", "chi")
# Natural-unit mass dimensions.  R has length dimension +1, hence mass -1.
FIELD_MASS_DIMENSIONS = np.asarray((0, 0, 0, -1, 1, 1), dtype=np.int64)
CONSTRAINT_NORMALIZATION_FLOOR = 1.0e-30


def _finite(name: str, value: Real) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a real number")
    answer = float(value)
    if not isfinite(answer):
        raise ValueError(f"{name} must be finite")
    return answer


def _array(name: str, value: object, *, ndim: int | None = None) -> np.ndarray:
    answer = np.asarray(value, dtype=np.float64)
    if ndim is not None and answer.ndim != ndim:
        raise ValueError(f"{name} must have {ndim} dimensions")
    if not np.all(np.isfinite(answer)):
        raise ValueError(f"{name} must contain only finite values")
    return answer


def dimensionless_state_components(
    u: object,
    p: object,
    q: object,
    *,
    length_unit: Real,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Scale ``(u, partial_t u, partial_r u)`` by canonical mass dimensions."""

    scale = _finite("length_unit", length_unit)
    if scale <= 0.0:
        raise ValueError("length_unit must be positive")
    arrays = tuple(_array(name, value) for name, value in (("u", u), ("p", p), ("q", q)))
    if any(array.shape[-1:] != (len(FIELD_ORDER),) for array in arrays):
        raise ValueError("u, p, and q must have final axis in canonical six-field order")
    u_scale = np.power(scale, FIELD_MASS_DIMENSIONS, dtype=np.float64)
    derivative_scale = np.power(scale, FIELD_MASS_DIMENSIONS + 1, dtype=np.float64)
    return arrays[0] * u_scale, arrays[1] * derivative_scale, arrays[2] * derivative_scale


def induced_scaled_infinity_norm(
    matrix: object,
    *,
    row_scales: object,
    column_scales: object,
) -> float:
    """Induced infinity norm after declared row/column nondimensionalization."""

    value = _array("matrix", matrix, ndim=2)
    rows = _array("row_scales", row_scales, ndim=1)
    columns = _array("column_scales", column_scales, ndim=1)
    if value.shape != (rows.size, columns.size):
        raise ValueError("matrix shape must match row and column scales")
    if np.any(rows <= 0.0) or np.any(columns <= 0.0):
        raise ValueError("row and column scales must be positive")
    scaled = rows[:, None] * value / columns[None, :]
    return float(np.max(np.sum(np.abs(scaled), axis=1), initial=0.0))


def normalized_componentwise_infinity(
    residuals: object,
    source_scales: object,
    *,
    floor: Real = CONSTRAINT_NORMALIZATION_FLOOR,
) -> float:
    """Return ``max_i |C_i|/max(|source_i|, floor)`` fail-closed."""

    residual = _array("residuals", residuals)
    source = _array("source_scales", source_scales)
    if residual.shape != source.shape:
        raise ValueError("residual and source-scale shapes must agree")
    minimum = _finite("floor", floor)
    if minimum < 0.0:
        raise ValueError("normalization floor must be nonnegative")
    denominator = np.maximum(np.abs(source), minimum)
    zero = denominator == 0.0
    if np.any(zero & (residual != 0.0)):
        raise ValueError("nonzero residual has zero normalization denominator")
    ratios = np.zeros_like(residual)
    np.divide(np.abs(residual), denominator, out=ratios, where=~zero)
    return float(np.max(ratios, initial=0.0))


def normalized_branch_displacement(
    candidate: object,
    previous_accepted: object,
    *,
    component_scales: object,
) -> float:
    """Infinity distance from the preceding accepted continuation root."""

    current = _array("candidate", candidate, ndim=1)
    previous = _array("previous_accepted", previous_accepted, ndim=1)
    scales = _array("component_scales", component_scales, ndim=1)
    if current.shape != previous.shape or current.shape != scales.shape:
        raise ValueError("branch vectors and component scales must have equal shape")
    if np.any(scales <= 0.0):
        raise ValueError("branch component scales must be positive")
    return float(np.max(np.abs(current - previous) / scales, initial=0.0))


def curvature_scale_and_ratio(
    kretschmann: object,
    ricci_eigenvalues: object,
    *,
    cutoff: Real,
) -> tuple[float, float]:
    """Return ``kappa`` and the dimensionless monitor ``(kappa/Lambda)^2``.

    ``kappa=max(|K|^(1/4), sqrt(|lambda_Ricci|))`` has mass dimension one.
    Squaring its cutoff ratio is the unique dimensionless reading of the
    protocol key ``curvature_scale_over_Lambda_squared``.
    """

    k = _array("kretschmann", kretschmann)
    ricci = _array("ricci_eigenvalues", ricci_eigenvalues)
    scale = _finite("cutoff", cutoff)
    if scale <= 0.0:
        raise ValueError("cutoff must be positive")
    kappa_k = float(np.max(np.power(np.abs(k), 0.25), initial=0.0))
    kappa_r = float(np.max(np.sqrt(np.abs(ricci)), initial=0.0))
    kappa = max(kappa_k, kappa_r)
    return kappa, (kappa / scale) ** 2


def _smooth_step(value: np.ndarray) -> np.ndarray:
    answer = np.zeros_like(value)
    answer[value >= 1.0] = 1.0
    mask = (value > 0.0) & (value < 1.0)
    x = value[mask]
    left = np.exp(-1.0 / x)
    right = np.exp(-1.0 / (1.0 - x))
    answer[mask] = left / (left + right)
    return answer


def compact_vacuum_buffer_window(
    coordinates: object,
    *,
    support_minimum: Real,
    support_maximum: Real,
    taper_width: Real,
) -> np.ndarray:
    """C-infinity compact window equal to one away from both vacuum edges."""

    x = _array("coordinates", coordinates, ndim=1)
    lower = _finite("support_minimum", support_minimum)
    upper = _finite("support_maximum", support_maximum)
    width = _finite("taper_width", taper_width)
    if width <= 0.0 or upper - lower <= 2.0 * width:
        raise ValueError("window support must contain two positive tapers and a plateau")
    return _smooth_step((x - lower) / width) * _smooth_step((upper - x) / width)


@dataclass(frozen=True, slots=True)
class SpectralSupport:
    angular_support: float
    support_bin: int
    nyquist_bin: int
    first_unresolved_bin: int
    alias_band_occupied: bool
    amplitude_threshold: float


def windowed_spectral_support(
    samples: object,
    proper_coordinates: object,
    *,
    window: object | None = None,
    relative_amplitude_floor: Real = 2.0**-40,
    absolute_amplitude_floor: Real = CONSTRAINT_NORMALIZATION_FLOOR,
    unresolved_top_fraction: Real = 1.0 / 8.0,
) -> SpectralSupport:
    """Resolve one causal temporal or proper-radial spectral support.

    Coordinates must be uniformly spaced proper time or proper distance.  The
    caller supplies a compact window for spatial data; a causal past-time
    window may be supplied for frequency data.  Occupation of the top eighth
    is a separate alias stop even when the physical-cutoff ratio is small.
    """

    data = _array("samples", samples, ndim=1)
    coordinate = _array("proper_coordinates", proper_coordinates, ndim=1)
    if data.size != coordinate.size or data.size < 16:
        raise ValueError("spectral support requires at least 16 matched samples")
    spacings = np.diff(coordinate)
    if np.any(spacings <= 0.0):
        raise ValueError("proper coordinates must be strictly increasing")
    spacing = float(spacings[0])
    tolerance = 4096.0 * np.finfo(np.float64).eps * max(1.0, abs(spacing))
    if not np.allclose(spacings, spacing, rtol=0.0, atol=tolerance):
        raise ValueError("proper coordinates must be uniformly spaced or resampled first")
    relative = _finite("relative_amplitude_floor", relative_amplitude_floor)
    absolute = _finite("absolute_amplitude_floor", absolute_amplitude_floor)
    top = _finite("unresolved_top_fraction", unresolved_top_fraction)
    if not 0.0 < relative < 1.0 or absolute < 0.0 or not 0.0 < top < 1.0:
        raise ValueError("spectral floor/fraction contract differs")
    if window is None:
        weights = np.ones(data.size, dtype=np.float64)
    else:
        weights = _array("window", window, ndim=1)
        if weights.shape != data.shape or np.any(weights < 0.0) or np.any(weights > 1.0):
            raise ValueError("spectral window must match samples and lie in [0,1]")
    transformed = np.fft.rfft(data * weights)
    amplitudes = np.abs(transformed)
    peak = float(np.max(amplitudes, initial=0.0))
    threshold = max(absolute, relative * peak)
    occupied = np.flatnonzero(amplitudes > threshold)
    support_bin = int(occupied[-1]) if occupied.size else 0
    nyquist_bin = amplitudes.size - 1
    first_unresolved = int(ceil((1.0 - top) * nyquist_bin))
    angular = 2.0 * pi * float(np.fft.rfftfreq(data.size, d=spacing)[support_bin])
    return SpectralSupport(
        angular_support=angular,
        support_bin=support_bin,
        nyquist_bin=nyquist_bin,
        first_unresolved_bin=first_unresolved,
        alias_band_occupied=support_bin >= first_unresolved,
        amplitude_threshold=threshold,
    )


@dataclass(frozen=True, slots=True)
class HealthThresholds:
    newton_residual_max: float
    newton_iteration_max: int
    kinetic_condition_max: float
    branch_displacement_max: float
    acceleration_max: float
    companion_deformation_max: float
    action_energy_deformation_max: float
    kinetic_deformation_max: float
    effective_planck_ratio_min: float
    characteristic_imaginary_max: float
    eigenframe_condition_max: float
    symmetrizer_coercivity_min: float
    normalized_constraint_max: float
    constraint_convergence_order_min: float
    phi_over_cutoff_max: float
    proper_frequency_over_cutoff_max: float
    proper_wavenumber_over_cutoff_max: float
    curvature_scale_squared_over_cutoff_squared_max: float

    def __post_init__(self) -> None:
        for field in fields(self):
            value = getattr(self, field.name)
            if field.name == "newton_iteration_max":
                if isinstance(value, bool) or not isinstance(value, int) or value < 0:
                    raise ValueError("newton_iteration_max must be a nonnegative integer")
            else:
                normalized = _finite(field.name, value)
                if normalized <= 0.0:
                    raise ValueError(f"{field.name} must be positive")
                object.__setattr__(self, field.name, normalized)


@dataclass(frozen=True, slots=True)
class StageHealthSnapshot:
    time: float
    # Global monitor-transaction serial, not a method-local Butcher stage.
    # SSPRK3 evaluates internal states at non-monotone physical abscissae, so
    # transaction order must be carried independently of ``time``.
    stage_index: int
    minimum_lapse: float
    minimum_radial_metric: float
    minimum_areal_radius_away_from_center: float
    minimum_hat_lorentzian_margin: float
    newton_residual_infinity: float
    newton_iterations: int
    newton_residual_monotonic: bool
    kinetic_condition_infinity: float
    normalized_branch_displacement: float
    acceleration_infinity: float
    companion_deformation_operator_2: float
    action_energy_deformation_operator_2: float
    kinetic_deformation_operator_2: float
    effective_planck_ratio: float
    characteristic_imaginary_part: float
    eigenframe_condition_number: float
    symmetrizer_coercivity: float
    normalized_constraint_infinity: float
    observed_constraint_convergence_order: float
    phi_over_cutoff: float
    proper_frequency_over_cutoff: float
    proper_wavenumber_over_cutoff: float
    curvature_scale_squared_over_cutoff_squared: float
    temporal_alias_band_occupied: bool
    radial_alias_band_occupied: bool
    boundary_margin_over_required_buffer: float

    def __post_init__(self) -> None:
        if isinstance(self.stage_index, bool) or not isinstance(self.stage_index, int) or self.stage_index < 0:
            raise ValueError("stage_index must be a nonnegative integer")
        if isinstance(self.newton_iterations, bool) or not isinstance(self.newton_iterations, int) or self.newton_iterations < 0:
            raise ValueError("newton_iterations must be a nonnegative integer")
        for name in ("newton_residual_monotonic", "temporal_alias_band_occupied", "radial_alias_band_occupied"):
            if not isinstance(getattr(self, name), bool):
                raise TypeError(f"{name} must be bool")
        for field in fields(self):
            if field.name in {"stage_index", "newton_iterations", "newton_residual_monotonic", "temporal_alias_band_occupied", "radial_alias_band_occupied"}:
                continue
            object.__setattr__(self, field.name, _finite(field.name, getattr(self, field.name)))


STOP_PRIORITY = (
    "nonpositive_lapse",
    "nonpositive_radial_metric",
    "nonpositive_areal_radius_away_from_center",
    "hat_cone_not_Lorentzian",
    "newton_residual_limit",
    "newton_iteration_limit",
    "newton_residual_not_monotonic",
    "kinetic_condition_limit",
    "branch_displacement_limit",
    "acceleration_limit",
    "all_covector_companion_deformation_limit",
    "action_energy_deformation_limit",
    "kinetic_deformation_limit",
    "effective_planck_ratio_limit",
    "characteristic_imaginary_part_limit",
    "eigenframe_condition_limit",
    "symmetrizer_coercivity_limit",
    "normalized_constraint_limit",
    "constraint_convergence_order_limit",
    "phi_cutoff_proxy_limit",
    "proper_frequency_cutoff_proxy_limit",
    "proper_wavenumber_cutoff_proxy_limit",
    "temporal_alias_band_occupied",
    "radial_alias_band_occupied",
    "curvature_cutoff_proxy_limit",
    "boundary_causal_buffer",
)


def failed_health_premises(
    snapshot: StageHealthSnapshot,
    thresholds: HealthThresholds,
) -> tuple[str, ...]:
    """Return every failed premise in deterministic scientific priority."""

    checks = {
        "nonpositive_lapse": snapshot.minimum_lapse <= 0.0,
        "nonpositive_radial_metric": snapshot.minimum_radial_metric <= 0.0,
        "nonpositive_areal_radius_away_from_center": snapshot.minimum_areal_radius_away_from_center <= 0.0,
        "hat_cone_not_Lorentzian": snapshot.minimum_hat_lorentzian_margin <= 0.0,
        "newton_residual_limit": snapshot.newton_residual_infinity >= thresholds.newton_residual_max,
        "newton_iteration_limit": snapshot.newton_iterations > thresholds.newton_iteration_max,
        "newton_residual_not_monotonic": not snapshot.newton_residual_monotonic,
        "kinetic_condition_limit": snapshot.kinetic_condition_infinity >= thresholds.kinetic_condition_max,
        "branch_displacement_limit": snapshot.normalized_branch_displacement >= thresholds.branch_displacement_max,
        "acceleration_limit": snapshot.acceleration_infinity >= thresholds.acceleration_max,
        "all_covector_companion_deformation_limit": snapshot.companion_deformation_operator_2 >= thresholds.companion_deformation_max,
        "action_energy_deformation_limit": snapshot.action_energy_deformation_operator_2 >= thresholds.action_energy_deformation_max,
        "kinetic_deformation_limit": snapshot.kinetic_deformation_operator_2 >= thresholds.kinetic_deformation_max,
        "effective_planck_ratio_limit": snapshot.effective_planck_ratio <= thresholds.effective_planck_ratio_min,
        "characteristic_imaginary_part_limit": snapshot.characteristic_imaginary_part >= thresholds.characteristic_imaginary_max,
        "eigenframe_condition_limit": snapshot.eigenframe_condition_number >= thresholds.eigenframe_condition_max,
        "symmetrizer_coercivity_limit": snapshot.symmetrizer_coercivity <= thresholds.symmetrizer_coercivity_min,
        "normalized_constraint_limit": snapshot.normalized_constraint_infinity >= thresholds.normalized_constraint_max,
        "constraint_convergence_order_limit": snapshot.observed_constraint_convergence_order <= thresholds.constraint_convergence_order_min,
        "phi_cutoff_proxy_limit": snapshot.phi_over_cutoff >= thresholds.phi_over_cutoff_max,
        "proper_frequency_cutoff_proxy_limit": snapshot.proper_frequency_over_cutoff >= thresholds.proper_frequency_over_cutoff_max,
        "proper_wavenumber_cutoff_proxy_limit": snapshot.proper_wavenumber_over_cutoff >= thresholds.proper_wavenumber_over_cutoff_max,
        "temporal_alias_band_occupied": snapshot.temporal_alias_band_occupied,
        "radial_alias_band_occupied": snapshot.radial_alias_band_occupied,
        "curvature_cutoff_proxy_limit": snapshot.curvature_scale_squared_over_cutoff_squared >= thresholds.curvature_scale_squared_over_cutoff_squared_max,
        "boundary_causal_buffer": snapshot.boundary_margin_over_required_buffer <= 0.0,
    }
    if tuple(checks) != STOP_PRIORITY:
        raise RuntimeError("health-stop priority and implementation order diverged")
    return tuple(name for name in STOP_PRIORITY if checks[name])


@dataclass(frozen=True, slots=True)
class HealthMonitorState:
    accepted_stage_count: int = 0
    last_accepted_time: float = 0.0
    last_accepted_stage_index: int = -1
    first_failed_premise: str | None = None
    first_failed_stage_index: int | None = None
    first_failed_time: float | None = None


class ClassicalHealthStop(RuntimeError):
    def __init__(self, state: HealthMonitorState, failures: Sequence[str]) -> None:
        self.reason = state.first_failed_premise
        self.state = state
        self.failures = tuple(failures)
        super().__init__(f"classical health monitor stopped at {self.reason}")


class ClassicalHealthMonitor:
    """Stateful transaction preserving the first failed scientific premise."""

    def __init__(self, thresholds: HealthThresholds) -> None:
        if not isinstance(thresholds, HealthThresholds):
            raise TypeError("thresholds must be HealthThresholds")
        self.thresholds = thresholds
        self.state = HealthMonitorState()

    def assess_and_accept(self, snapshot: StageHealthSnapshot) -> HealthMonitorState:
        if not isinstance(snapshot, StageHealthSnapshot):
            raise TypeError("snapshot must be StageHealthSnapshot")
        if self.state.first_failed_premise is not None:
            raise ClassicalHealthStop(self.state, (self.state.first_failed_premise,))
        if snapshot.stage_index <= self.state.last_accepted_stage_index:
            raise ValueError(
                "stage_index must be a strictly increasing global transaction serial"
            )
        failures = failed_health_premises(snapshot, self.thresholds)
        if failures:
            self.state = HealthMonitorState(
                accepted_stage_count=self.state.accepted_stage_count,
                last_accepted_time=self.state.last_accepted_time,
                last_accepted_stage_index=self.state.last_accepted_stage_index,
                first_failed_premise=failures[0],
                first_failed_stage_index=snapshot.stage_index,
                first_failed_time=snapshot.time,
            )
            raise ClassicalHealthStop(self.state, failures)
        self.state = HealthMonitorState(
            accepted_stage_count=self.state.accepted_stage_count + 1,
            last_accepted_time=snapshot.time,
            last_accepted_stage_index=snapshot.stage_index,
        )
        return self.state
