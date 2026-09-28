"""Outcome-neutral PROTO12 GR-0 runtime components.

PROTO12 keeps PROTO11's independently evolved ``(u,p,q)`` state, reference-
balanced radial derivative map, constraints, physical inputs, time
integrators, thresholds, and stop rules.  It changes only two numerical
premises:

* the unchanged GR-0 ACT1/VAR1/REF1 affine source is evaluated through the
  SRC4 tensor-contracted backend in fixed, bounded batches; and
* both adjacent spatial-spectral pairs are classified by the frozen PROTO12
  direct-or-fully-guarded diagnostic rule.

Every raw source residual, raw spectral ratio, and predecessor assessment
remains visible.  Diagnostic saturation is neither a continuum-error bound
nor evidence of physical resolution.  This module performs no I/O and grants
no run or physical authorization by itself.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite
from numbers import Integral, Real
from typing import Mapping, Sequence

import numpy as np

from .cal4_common_event_diagnosis import OwnedConstraintAdmission
from .gr0_calibration import FIELD_COUNT
from .numerical_engine import (
    EvolutionRHS,
    EvolutionState,
    SBPFirstDerivative,
    UniformRadialGrid,
)
from .proto11_runtime import (
    Proto11GR0CommonEventAssessment,
    exact_spherical_minkowski_reference,
    proto11_gr0_common_event,
    reference_balanced_spatial_derivatives,
)
from .proto4_admission import SpectralPowerBudget, SpectralThresholds
from .proto8_runtime import Proto8GR0CommonEventAssessment
from .proto9_runtime import PROTO9_POINT_COUNTS
from .spectral_sensitivity import (
    SpectralTailSensitivity,
    ThreeGridProfileConvergence,
)
from .spectral_sensitivity_v12 import (
    PairwiseResolvedOrSaturatedAdmission,
    pairwise_resolved_or_saturated_admission,
)
from .src4_vectorized_reference_source import (
    diagnose_gr0_vectorized_reference_accelerations,
)
from .vectorized_source import (
    ADM_CENTER_PARITIES,
    regular_center_acceleration_limit,
)


PROTO12_POINT_COUNTS = PROTO9_POINT_COUNTS
PROTO12_SOURCE_POINT_BATCH_SIZE = 2048


def _finite(name: str, value: object, *, positive: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a finite real scalar")
    answer = float(value)
    if not isfinite(answer):
        raise ValueError(f"{name} must be finite")
    if positive and answer <= 0.0:
        raise ValueError(f"{name} must be positive")
    return answer


def _point_batch_size(value: object) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral):
        raise TypeError("point_batch_size must be an integer")
    answer = int(value)
    if not 1 <= answer <= PROTO12_SOURCE_POINT_BATCH_SIZE:
        raise ValueError(
            "point_batch_size must lie in "
            f"[1,{PROTO12_SOURCE_POINT_BATCH_SIZE}]"
        )
    return answer


def _source_arrays(
    u: object,
    p: object,
    q: object,
    p_r: object,
    q_r: object,
    radii: object,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    fields = tuple(
        np.asarray(value, dtype=np.float64) for value in (u, p, q, p_r, q_r)
    )
    radius = np.asarray(radii, dtype=np.float64)
    if (
        any(value.ndim != 2 for value in fields)
        or len({value.shape for value in fields}) != 1
        or fields[0].shape[1] != FIELD_COUNT
        or radius.ndim != 1
        or radius.shape[0] != fields[0].shape[0]
        or not all(np.all(np.isfinite(value)) for value in fields)
        or not np.all(np.isfinite(radius))
    ):
        raise ValueError(
            "PROTO12 source arrays must share one finite six-field grid"
        )
    if radius.size == 0 or np.any(radius <= 0.0):
        raise ValueError("PROTO12 source radii must be positive and nonempty")
    return (*fields, radius)


@dataclass(frozen=True, slots=True)
class Proto12ChunkedSourceDiagnosis:
    """Complete SRC4 pointwise affine-root evidence after fixed batching."""

    accelerations: np.ndarray
    residuals: np.ndarray
    residual_infinity: float
    raw_tolerance: float
    raw_gate_passed: bool
    kinetic_condition_infinity_maximum: float
    maximum_refinement_iterations_used: int
    residual_decreased_monotonically: bool
    roundoff_floor_reached: bool
    maximum_relative_acceleration_correction: float
    maximum_residual_point_index: int
    maximum_residual_row_index: int
    maximum_residual_radius: float
    point_batch_size: int
    point_batch_count: int

    def __post_init__(self) -> None:
        acceleration = np.asarray(self.accelerations, dtype=np.float64)
        residual = np.asarray(self.residuals, dtype=np.float64)
        if (
            acceleration.ndim != 2
            or acceleration.shape != residual.shape
            or acceleration.shape[1] != FIELD_COUNT
            or acceleration.shape[0] < 1
            or not np.all(np.isfinite(acceleration))
            or not np.all(np.isfinite(residual))
        ):
            raise ValueError("PROTO12 chunked source arrays differ")
        for name in (
            "residual_infinity",
            "raw_tolerance",
            "kinetic_condition_infinity_maximum",
            "maximum_relative_acceleration_correction",
            "maximum_residual_radius",
        ):
            value = _finite(name, getattr(self, name))
            if name in {"raw_tolerance", "kinetic_condition_infinity_maximum", "maximum_residual_radius"} and value <= 0.0:
                raise ValueError(f"{name} must be positive")
            if name in {"residual_infinity", "maximum_relative_acceleration_correction"} and value < 0.0:
                raise ValueError(f"{name} must be nonnegative")
            object.__setattr__(self, name, value)
        for name in (
            "maximum_refinement_iterations_used",
            "maximum_residual_point_index",
            "maximum_residual_row_index",
            "point_batch_size",
            "point_batch_count",
        ):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, Integral) or value < 0:
                raise ValueError(f"{name} must be a nonnegative integer")
            object.__setattr__(self, name, int(value))
        if self.maximum_residual_point_index >= acceleration.shape[0]:
            raise ValueError("PROTO12 maximum residual point lies outside the grid")
        if self.maximum_residual_row_index >= FIELD_COUNT:
            raise ValueError("PROTO12 maximum residual row lies outside REF1")
        if self.point_batch_size < 1 or self.point_batch_count < 1:
            raise ValueError("PROTO12 batching evidence must be nonempty")
        for name in (
            "raw_gate_passed",
            "residual_decreased_monotonically",
            "roundoff_floor_reached",
        ):
            if not isinstance(getattr(self, name), bool):
                raise TypeError(f"{name} must be bool")
        acceleration = np.ascontiguousarray(acceleration).copy()
        residual = np.ascontiguousarray(residual).copy()
        acceleration.setflags(write=False)
        residual.setflags(write=False)
        object.__setattr__(self, "accelerations", acceleration)
        object.__setattr__(self, "residuals", residual)


def diagnose_gr0_proto12_accelerations_chunked(
    u: object,
    p: object,
    q: object,
    p_r: object,
    q_r: object,
    radii: object,
    *,
    raw_tolerance: Real = 1.0e-12,
    condition_number_maximum: Real = 1.0e10,
    maximum_refinement_iterations: int = 16,
    point_batch_size: int = PROTO12_SOURCE_POINT_BATCH_SIZE,
) -> Proto12ChunkedSourceDiagnosis:
    """Evaluate the unchanged SRC4 systems with bounded peak memory."""

    base, velocity, gradient, mixed, radial_second, radius = _source_arrays(
        u, p, q, p_r, q_r, radii
    )
    tolerance = _finite("raw_tolerance", raw_tolerance, positive=True)
    condition_limit = _finite(
        "condition_number_maximum", condition_number_maximum, positive=True
    )
    batch_size = _point_batch_size(point_batch_size)
    if (
        isinstance(maximum_refinement_iterations, bool)
        or not isinstance(maximum_refinement_iterations, int)
        or not 0 <= maximum_refinement_iterations <= 16
    ):
        raise ValueError("maximum_refinement_iterations must lie in [0,16]")
    diagnoses = []
    accelerations = []
    residuals = []
    for start in range(0, radius.size, batch_size):
        stop = min(start + batch_size, radius.size)
        diagnosis = diagnose_gr0_vectorized_reference_accelerations(
            base[start:stop],
            velocity[start:stop],
            gradient[start:stop],
            mixed[start:stop],
            radial_second[start:stop],
            radius[start:stop],
            raw_tolerance=tolerance,
            condition_number_maximum=condition_limit,
            maximum_refinement_iterations=maximum_refinement_iterations,
        )
        diagnoses.append(diagnosis)
        accelerations.append(diagnosis.accelerations)
        residuals.append(diagnosis.residuals)
    acceleration = np.concatenate(accelerations, axis=0)
    residual = np.concatenate(residuals, axis=0)
    maximum = np.unravel_index(int(np.argmax(np.abs(residual))), residual.shape)
    return Proto12ChunkedSourceDiagnosis(
        accelerations=acceleration,
        residuals=residual,
        residual_infinity=float(np.max(np.abs(residual), initial=0.0)),
        raw_tolerance=tolerance,
        raw_gate_passed=all(item.raw_gate_passed for item in diagnoses),
        kinetic_condition_infinity_maximum=max(
            item.kinetic_condition_infinity_maximum for item in diagnoses
        ),
        maximum_refinement_iterations_used=max(
            len(item.iterations) - 1 for item in diagnoses
        ),
        residual_decreased_monotonically=all(
            item.residual_decreased_until_floor for item in diagnoses
        ),
        roundoff_floor_reached=any(
            item.roundoff_floor_reached for item in diagnoses
        ),
        maximum_relative_acceleration_correction=max(
            item.iterations[-1].relative_acceleration_correction
            for item in diagnoses
        ),
        maximum_residual_point_index=int(maximum[0]),
        maximum_residual_row_index=int(maximum[1]),
        maximum_residual_radius=float(radius[maximum[0]]),
        point_batch_size=batch_size,
        point_batch_count=len(diagnoses),
    )


class Proto12GR0EvolutionOperator:
    """PROTO11 derivative map plus SRC4's unchanged-equation source backend."""

    def __init__(
        self,
        grid: UniformRadialGrid,
        *,
        spatial_order: int,
        ko_dissipation: Real,
        raw_tolerance: Real,
        kinetic_condition_maximum: Real,
        maximum_refinement_iterations: int = 16,
        point_batch_size: int = PROTO12_SOURCE_POINT_BATCH_SIZE,
    ) -> None:
        if not isinstance(grid, UniformRadialGrid) or grid.minimum != 0.0:
            raise ValueError("PROTO12 GR-0 operator requires a centre grid")
        self.grid = grid
        self.derivative = SBPFirstDerivative(grid, spatial_order)
        self.ko_dissipation = _finite("ko_dissipation", ko_dissipation)
        if self.ko_dissipation < 0.0:
            raise ValueError("ko_dissipation must be nonnegative")
        self.raw_tolerance = _finite(
            "raw_tolerance", raw_tolerance, positive=True
        )
        self.kinetic_condition_maximum = _finite(
            "kinetic_condition_maximum",
            kinetic_condition_maximum,
            positive=True,
        )
        if (
            isinstance(maximum_refinement_iterations, bool)
            or not isinstance(maximum_refinement_iterations, int)
            or not 0 <= maximum_refinement_iterations <= 16
        ):
            raise ValueError("maximum_refinement_iterations must lie in [0,16]")
        self.maximum_refinement_iterations = maximum_refinement_iterations
        self.point_batch_size = _point_batch_size(point_batch_size)

    def __call__(self, _time: float, state: EvolutionState) -> EvolutionRHS:
        if not isinstance(state, EvolutionState) or state.shape != (
            self.grid.point_count,
            FIELD_COUNT,
        ):
            raise ValueError("PROTO12 GR-0 RHS state shape differs")
        p_r, q_r, differentiated_u = reference_balanced_spatial_derivatives(
            state, self.derivative
        )
        source = diagnose_gr0_proto12_accelerations_chunked(
            state.u[1:],
            state.p[1:],
            state.q[1:],
            p_r[1:],
            q_r[1:],
            self.grid.coordinates[1:],
            raw_tolerance=self.raw_tolerance,
            condition_number_maximum=self.kinetic_condition_maximum,
            maximum_refinement_iterations=self.maximum_refinement_iterations,
            point_batch_size=self.point_batch_size,
        )
        reference = exact_spherical_minkowski_reference(self.grid)
        exact_reference = (
            np.array_equal(state.u, reference.u)
            and np.array_equal(state.p, reference.p)
            and np.array_equal(state.q, reference.q)
        )
        acceleration = np.empty_like(state.p)
        acceleration[1:] = 0.0 if exact_reference else source.accelerations
        center = regular_center_acceleration_limit(acceleration[1:5])
        acceleration[0] = center.acceleration
        dissipation = self.derivative.kreiss_oliger_dissipation(
            state.p,
            coefficient=self.ko_dissipation,
            center_parities=ADM_CENTER_PARITIES,
        )
        reduction = state.q - differentiated_u
        coordinate_speed = np.abs(state.u[:, 1]) + (
            1.2 * state.u[:, 0] / state.u[:, 2]
        )
        return EvolutionRHS(
            state.p,
            acceleration + dissipation,
            p_r,
            {
                "source_residual_infinity": source.residual_infinity,
                "source_raw_gate_passed": source.raw_gate_passed,
                "source_refinement_iterations": (
                    source.maximum_refinement_iterations_used
                ),
                "source_residual_decreased_monotonically": (
                    source.residual_decreased_monotonically
                ),
                "source_roundoff_floor_reached": source.roundoff_floor_reached,
                "source_relative_acceleration_correction": (
                    source.maximum_relative_acceleration_correction
                ),
                "source_maximum_residual_point_index": (
                    source.maximum_residual_point_index
                ),
                "source_maximum_residual_row_index": (
                    source.maximum_residual_row_index
                ),
                "source_maximum_residual_radius": source.maximum_residual_radius,
                "kinetic_condition_infinity": (
                    source.kinetic_condition_infinity_maximum
                ),
                "acceleration_infinity": float(
                    np.max(np.abs(acceleration), initial=0.0)
                ),
                "center_acceleration_estimator_infinity": (
                    center.estimator_infinity
                ),
                "reduction_constraint_infinity": float(
                    np.max(np.abs(reduction), initial=0.0)
                ),
                "coordinate_speed_upper": float(
                    np.max(coordinate_speed, initial=0.0)
                ),
                "minimum_lapse": float(np.min(state.u[:, 0])),
                "minimum_radial_metric": float(np.min(state.u[:, 2])),
                "minimum_areal_radius_away_from_center": float(
                    np.min(state.u[1:, 3])
                ),
                "SRC4_tensor_contracted_reference_source": True,
                "PROTO12_source_point_batch_size": self.point_batch_size,
                "PROTO12_source_point_batch_count": source.point_batch_count,
                "PROTO11_reference_state_map_retained": True,
                "PROTO11_exact_reference_equilibrium_applied": exact_reference,
                "PROTO11_interior_q_reprojected": False,
            },
        )


@dataclass(frozen=True, slots=True)
class Proto12GR0CommonEventAssessment:
    """PROTO11 constraints/raw evidence plus PROTO12 pairwise classification."""

    method: str
    coordinate_time: float
    point_counts: tuple[int, int, int]
    accepted_stage_counts: tuple[int, int, int]
    constraint_admission: OwnedConstraintAdmission
    legacy_PROTO11_common_event: Proto11GR0CommonEventAssessment
    raw_PROTO9_common_event: Proto8GR0CommonEventAssessment
    raw_spatial_spectral_admission: object
    spectral_budgets: tuple[
        Mapping[str, SpectralPowerBudget],
        Mapping[str, SpectralPowerBudget],
        Mapping[str, SpectralPowerBudget],
    ]
    spectral_tail_sensitivities: tuple[
        Mapping[str, SpectralTailSensitivity],
        Mapping[str, SpectralTailSensitivity],
        Mapping[str, SpectralTailSensitivity],
    ]
    profile_convergence: ThreeGridProfileConvergence
    spatial_spectral_admission: PairwiseResolvedOrSaturatedAdmission
    maximum_spatial_round_trip_interpolation_infinity: float
    raw_PROTO9_admission_passed: bool
    legacy_PROTO10_admission_passed: bool
    legacy_PROTO11_admission_passed: bool
    reference_balanced_constraint_admission_passed: bool
    reference_state_map_applied: bool
    SRC4_source_backend_bound: bool
    interior_q_reprojected: bool
    admission_passed: bool


def proto12_gr0_common_event(
    states: Sequence[EvolutionState],
    grids: Sequence[UniformRadialGrid],
    *,
    accepted_stage_counts: Sequence[int],
    method: str,
    coordinate_time: Real,
    cutoff: Real = 16.0,
    measurement_radius_maximum: Real = 24.0,
    taper_fraction: Real = 1.0 / 8.0,
    fixed_outer_rows: int = 4,
    spectral_thresholds: SpectralThresholds | None = None,
) -> Proto12GR0CommonEventAssessment:
    """Recompute PROTO11, then apply PROTO12 to both adjacent pairs."""

    legacy = proto11_gr0_common_event(
        states,
        grids,
        accepted_stage_counts=accepted_stage_counts,
        method=method,
        coordinate_time=coordinate_time,
        cutoff=cutoff,
        measurement_radius_maximum=measurement_radius_maximum,
        taper_fraction=taper_fraction,
        fixed_outer_rows=fixed_outer_rows,
        spectral_thresholds=spectral_thresholds,
    )
    if legacy.point_counts != PROTO12_POINT_COUNTS:
        raise RuntimeError("PROTO12 lost the frozen calibration ladder")
    pairwise = pairwise_resolved_or_saturated_admission(
        PROTO12_POINT_COUNTS,
        legacy.spectral_budgets,
        legacy.spectral_tail_sensitivities,
        legacy.profile_convergence,
        thresholds=spectral_thresholds,
    )
    raw = legacy.raw_spatial_spectral_admission
    if (
        dict(raw.field_power_tail_ratios)
        != dict(pairwise.field_power_tail_ratios)
        or dict(raw.derivative_power_tail_ratios)
        != dict(pairwise.derivative_power_tail_ratios)
    ):
        raise RuntimeError("PROTO12 recomputation changed a raw nested-tail ratio")
    if pairwise.medium_and_fine_individual_budgets_passed is not all(
        legacy.spectral_budgets[index][name].individual_admission_passed
        for index in (1, 2)
        for name in legacy.spectral_budgets[index]
    ):
        raise RuntimeError("PROTO12 changed an absolute spectral budget")
    admitted = (
        legacy.constraint_admission.admission_passed
        and pairwise.admission_passed
    )
    stages = tuple(int(value) for value in accepted_stage_counts)
    if len(stages) != 3:
        raise ValueError("PROTO12 requires three accepted-stage counts")
    return Proto12GR0CommonEventAssessment(
        method=legacy.method,
        coordinate_time=legacy.coordinate_time,
        point_counts=PROTO12_POINT_COUNTS,
        accepted_stage_counts=(stages[0], stages[1], stages[2]),
        constraint_admission=legacy.constraint_admission,
        legacy_PROTO11_common_event=legacy,
        raw_PROTO9_common_event=legacy.raw_PROTO9_common_event,
        raw_spatial_spectral_admission=raw,
        spectral_budgets=legacy.spectral_budgets,
        spectral_tail_sensitivities=legacy.spectral_tail_sensitivities,
        profile_convergence=legacy.profile_convergence,
        spatial_spectral_admission=pairwise,
        maximum_spatial_round_trip_interpolation_infinity=(
            legacy.maximum_spatial_round_trip_interpolation_infinity
        ),
        raw_PROTO9_admission_passed=legacy.raw_PROTO9_admission_passed,
        legacy_PROTO10_admission_passed=(
            legacy.legacy_PROTO10_admission_passed
        ),
        legacy_PROTO11_admission_passed=legacy.admission_passed,
        reference_balanced_constraint_admission_passed=(
            legacy.reference_balanced_constraint_admission_passed
        ),
        reference_state_map_applied=True,
        SRC4_source_backend_bound=True,
        interior_q_reprojected=False,
        admission_passed=admitted,
    )


__all__ = [
    "PROTO12_POINT_COUNTS",
    "PROTO12_SOURCE_POINT_BATCH_SIZE",
    "Proto12ChunkedSourceDiagnosis",
    "Proto12GR0CommonEventAssessment",
    "Proto12GR0EvolutionOperator",
    "diagnose_gr0_proto12_accelerations_chunked",
    "proto12_gr0_common_event",
]
