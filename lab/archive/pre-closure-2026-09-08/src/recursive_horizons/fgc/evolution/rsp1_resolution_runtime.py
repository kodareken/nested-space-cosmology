"""Bounded SRC3 runtime for the amplitude-three resolution-spectrum study.

RSP1 is not a successor collapse protocol.  It asks one numerical-premise
question left open by CAL8/PREF10: does the evolved ``phi/Lambda`` derivative
tail enter the unchanged direct ``< 1/4`` contraction regime when one genuinely
new nested grid is added?

The module therefore keeps the PROTO11 physical state, reference-balanced
radial derivative map, methods, constraints, and thresholds, but evaluates the
unchanged GR-0 REF1 rows through SRC3's reference-covariant arithmetic.  The
pointwise affine systems are processed in fixed batches so the new 16,385-point
grid has a bounded memory footprint.  Batching changes neither an equation nor
the per-point solve.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite
from numbers import Integral, Real
from typing import Mapping, Sequence

import numpy as np

from .cal4_common_event_diagnosis import (
    FinestPairNestedSpectralAdmission,
    OwnedConstraintAdmission,
    OwnedConstraintSnapshot,
    finest_pair_nested_spectral_admission,
    owned_constraint_admission,
)
from .calibration_runtime import (
    CONSTRAINT_ROUNDOFF_OPERATION_BUDGET,
    SemidiscreteConstraintSnapshot,
)
from .gr0_calibration import FIELD_COUNT, Q_CENTER_PARITIES
from .health_monitor import compact_vacuum_buffer_window
from .numerical_engine import (
    EvolutionRHS,
    EvolutionState,
    SBPFirstDerivative,
    UniformRadialGrid,
)
from .proto4_admission import (
    PROTO4_CONSTRAINT_ORDER,
    PROTO4_SPECTRAL_FIELD_ORDER,
    ConstraintNorms,
    SpectralPowerBudget,
    SpectralThresholds,
    normalized_constraint_norms,
    proper_radial_profile,
    spectral_field_budgets,
)
from .proto5_runtime import PROTO5_METHOD_CONTRACT
from .proto8_runtime import PROTO8_CONSTRAINT_OWNERSHIP
from .proto11_runtime import (
    exact_spherical_minkowski_reference,
    reference_balanced_spatial_derivatives,
)
from .spectral_sensitivity import (
    ThreeGridProfileConvergence,
    three_grid_profile_convergence,
)
from .src3_reference_balanced_source import (
    diagnose_gr0_reference_balanced_accelerations,
    gr0_reference_balanced_ref1_residual_batch,
)
from .vectorized_source import (
    ADM_CENTER_PARITIES,
    regular_center_acceleration_limit,
)


RSP1_POINT_COUNTS = (4097, 8193, 16385)
RSP1_SOURCE_POINT_BATCH_SIZE = 2048
RSP1_TARGET_FIELD = "phi_over_Lambda"


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
    if answer < 1 or answer > RSP1_SOURCE_POINT_BATCH_SIZE:
        raise ValueError(
            f"point_batch_size must lie in [1,{RSP1_SOURCE_POINT_BATCH_SIZE}]"
        )
    return answer


def _state_arrays(
    u: object,
    p: object,
    q: object,
    p_r: object,
    q_r: object,
    radii: object,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    values = tuple(
        np.asarray(value, dtype=np.float64) for value in (u, p, q, p_r, q_r)
    )
    radius = np.asarray(radii, dtype=np.float64)
    if (
        any(value.ndim != 2 for value in values)
        or len({value.shape for value in values}) != 1
        or values[0].shape[1] != FIELD_COUNT
        or radius.ndim != 1
        or radius.shape[0] != values[0].shape[0]
        or not all(np.all(np.isfinite(value)) for value in values)
        or not np.all(np.isfinite(radius))
    ):
        raise ValueError("RSP1 source arrays must share one finite six-field grid")
    if radius.size == 0 or np.any(radius <= 0.0):
        raise ValueError("RSP1 source radii must be positive and nonempty")
    return (*values, radius)


@dataclass(frozen=True, slots=True)
class RSP1ChunkedSourceDiagnosis:
    """Complete pointwise affine-root evidence after bounded batching."""

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
            raise ValueError("RSP1 chunked source arrays differ")
        for name in (
            "residual_infinity",
            "raw_tolerance",
            "kinetic_condition_infinity_maximum",
            "maximum_relative_acceleration_correction",
            "maximum_residual_radius",
        ):
            value = _finite(name, getattr(self, name))
            if name != "maximum_relative_acceleration_correction" and value <= 0.0:
                if name != "residual_infinity" or value != 0.0:
                    raise ValueError(f"{name} has an invalid sign")
            if name == "maximum_relative_acceleration_correction" and value < 0.0:
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
            raise ValueError("RSP1 maximum residual point lies outside the grid")
        if self.maximum_residual_row_index >= FIELD_COUNT:
            raise ValueError("RSP1 maximum residual row lies outside REF1")
        if self.point_batch_size < 1 or self.point_batch_count < 1:
            raise ValueError("RSP1 batching evidence must be nonempty")
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


def diagnose_gr0_reference_balanced_accelerations_chunked(
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
    point_batch_size: int = RSP1_SOURCE_POINT_BATCH_SIZE,
) -> RSP1ChunkedSourceDiagnosis:
    """Solve all pointwise SRC3 affine systems with bounded peak memory."""

    base, velocity, gradient, mixed, radial_second, radius = _state_arrays(
        u, p, q, p_r, q_r, radii
    )
    tolerance = _finite("raw_tolerance", raw_tolerance, positive=True)
    condition_limit = _finite(
        "condition_number_maximum", condition_number_maximum, positive=True
    )
    batch_size = _point_batch_size(point_batch_size)
    accelerations: list[np.ndarray] = []
    residuals: list[np.ndarray] = []
    diagnoses = []
    for start in range(0, radius.size, batch_size):
        stop = min(start + batch_size, radius.size)
        diagnosis = diagnose_gr0_reference_balanced_accelerations(
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
    residual_infinity = float(np.max(np.abs(residual), initial=0.0))
    return RSP1ChunkedSourceDiagnosis(
        accelerations=acceleration,
        residuals=residual,
        residual_infinity=residual_infinity,
        raw_tolerance=tolerance,
        raw_gate_passed=(
            residual_infinity < tolerance
            and all(item.raw_gate_passed for item in diagnoses)
        ),
        kinetic_condition_infinity_maximum=max(
            item.kinetic_condition_infinity_maximum for item in diagnoses
        ),
        maximum_refinement_iterations_used=max(
            len(item.iterations) - 1 for item in diagnoses
        ),
        residual_decreased_monotonically=all(
            item.residual_decreased_until_floor for item in diagnoses
        ),
        roundoff_floor_reached=any(item.roundoff_floor_reached for item in diagnoses),
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


class RSP1GR0EvolutionOperator:
    """PROTO11 semidiscrete map with SRC3's unchanged-equation evaluator."""

    def __init__(
        self,
        grid: UniformRadialGrid,
        *,
        spatial_order: int,
        ko_dissipation: Real,
        raw_tolerance: Real,
        kinetic_condition_maximum: Real,
        point_batch_size: int = RSP1_SOURCE_POINT_BATCH_SIZE,
    ) -> None:
        if not isinstance(grid, UniformRadialGrid) or grid.minimum != 0.0:
            raise ValueError("RSP1 GR-0 operator requires a regular-centre grid")
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
        self.point_batch_size = _point_batch_size(point_batch_size)

    def __call__(self, _time: float, state: EvolutionState) -> EvolutionRHS:
        if not isinstance(state, EvolutionState) or state.shape != (
            self.grid.point_count,
            FIELD_COUNT,
        ):
            raise ValueError("RSP1 GR-0 RHS state shape differs")
        p_r, q_r, differentiated_u = reference_balanced_spatial_derivatives(
            state, self.derivative
        )
        source = diagnose_gr0_reference_balanced_accelerations_chunked(
            state.u[1:],
            state.p[1:],
            state.q[1:],
            p_r[1:],
            q_r[1:],
            self.grid.coordinates[1:],
            raw_tolerance=self.raw_tolerance,
            condition_number_maximum=self.kinetic_condition_maximum,
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
                "SRC3_reference_covariant_source": True,
                "RSP1_source_point_batch_size": self.point_batch_size,
                "RSP1_source_point_batch_count": source.point_batch_count,
                "PROTO11_reference_state_map_retained": True,
                "PROTO11_exact_reference_equilibrium_applied": exact_reference,
                "PROTO11_interior_q_reprojected": False,
            },
        )


def _src3_constraint_arrays(
    state: EvolutionState,
    grid: UniformRadialGrid,
    *,
    diagnostic_spatial_order: int,
    point_batch_size: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    derivative = SBPFirstDerivative(grid, diagnostic_spatial_order)
    p_r, q_r, differentiated_u = reference_balanced_spatial_derivatives(
        state, derivative
    )
    residuals = np.zeros((grid.point_count, len(PROTO4_CONSTRAINT_ORDER)))
    terms = np.zeros((grid.point_count, len(PROTO4_CONSTRAINT_ORDER), 3))
    batch_size = _point_batch_size(point_batch_size)
    for start in range(1, grid.point_count, batch_size):
        stop = min(start + batch_size, grid.point_count)
        source = gr0_reference_balanced_ref1_residual_batch(
            state.u[start:stop],
            state.p[start:stop],
            state.q[start:stop],
            np.zeros((1, stop - start, FIELD_COUNT), dtype=np.float64),
            p_r[start:stop],
            q_r[start:stop],
            grid.coordinates[start:stop],
        )
        residuals[start:stop, 0] = source.hamiltonian_constraint[0]
        residuals[start:stop, 1] = source.momentum_constraint[0]
        residuals[start:stop, 2:4] = source.gauge_constraint[0, :, :2]
        metric = source.unredefined_metric_residual[0]
        shift = state.u[start:stop, 1]
        terms[start:stop, 0, 0] = metric[:, 0, 0]
        terms[start:stop, 0, 1] = -2.0 * shift * metric[:, 0, 1]
        terms[start:stop, 0, 2] = shift**2 * metric[:, 1, 1]
        terms[start:stop, 1, 0] = metric[:, 0, 1]
        terms[start:stop, 1, 1] = -shift * metric[:, 1, 1]
        terms[start:stop, 2:4, 0] = source.gauge_constraint[0, :, :2]
    residuals[:, 4:] = state.q - differentiated_u
    terms[:, 4:, 0] = state.q
    terms[:, 4:, 1] = -differentiated_u
    tolerance = 2048.0 * np.finfo(np.float64).eps
    if not np.allclose(
        np.sum(terms, axis=2), residuals, rtol=tolerance, atol=tolerance
    ):
        raise RuntimeError("RSP1 SRC3 constraint terms do not reconstruct")
    return residuals, terms, differentiated_u


def rsp1_semidiscrete_constraint_snapshot(
    state: EvolutionState,
    grid: UniformRadialGrid,
    *,
    coordinate_time: Real,
    diagnostic_spatial_order: int,
    planck_mass: Real = 2.0,
    length_unit: Real = 4.0,
    point_batch_size: int = RSP1_SOURCE_POINT_BATCH_SIZE,
    roundoff_operation_budget: int = CONSTRAINT_ROUNDOFF_OPERATION_BUDGET,
) -> SemidiscreteConstraintSnapshot:
    """Evaluate all GR-0 constraints through SRC3's arithmetic graph."""

    if not isinstance(state, EvolutionState) or state.shape != (
        grid.point_count,
        FIELD_COUNT,
    ):
        raise ValueError("RSP1 constraint state and grid differ")
    if (
        isinstance(roundoff_operation_budget, bool)
        or not isinstance(roundoff_operation_budget, int)
        or roundoff_operation_budget < 1
    ):
        raise ValueError("roundoff_operation_budget must be a positive integer")
    residuals, terms, _differentiated_u = _src3_constraint_arrays(
        state,
        grid,
        diagnostic_spatial_order=diagnostic_spatial_order,
        point_batch_size=point_batch_size,
    )
    raw = normalized_constraint_norms(
        residuals,
        terms,
        planck_mass=planck_mass,
        length_unit=length_unit,
    )
    enclosure = np.full(
        len(PROTO4_CONSTRAINT_ORDER),
        roundoff_operation_budget * np.finfo(np.float64).eps,
        dtype=np.float64,
    )
    effective_components = np.maximum(raw.component_infinity - enclosure, 0.0)
    effective = ConstraintNorms(
        sample_count=raw.sample_count,
        component_names=raw.component_names,
        component_infinity=effective_components,
        global_infinity=float(np.max(effective_components, initial=0.0)),
        minimum_denominator=raw.minimum_denominator,
        maximum_denominator=raw.maximum_denominator,
    )
    return SemidiscreteConstraintSnapshot(
        point_count=grid.point_count,
        coordinate_time=_finite("coordinate_time", coordinate_time),
        diagnostic_spatial_order=diagnostic_spatial_order,
        raw_norms=raw,
        effective_norms_for_order_only=effective,
        normalized_roundoff_zero_enclosure=enclosure,
        residuals=residuals,
        unredefined_term_contributions=terms,
    )


def rsp1_owned_constraint_snapshot(
    state: EvolutionState,
    grid: UniformRadialGrid,
    *,
    coordinate_time: Real,
    diagnostic_spatial_order: int,
    fixed_outer_rows: int,
    accepted_stage_count: int,
    planck_mass: Real = 2.0,
    length_unit: Real = 4.0,
    point_batch_size: int = RSP1_SOURCE_POINT_BATCH_SIZE,
) -> OwnedConstraintSnapshot:
    """Apply unchanged PROTO8 ownership to the SRC3 constraint surface."""

    if (
        isinstance(fixed_outer_rows, bool)
        or not isinstance(fixed_outer_rows, int)
        or fixed_outer_rows < 1
    ):
        raise ValueError("fixed_outer_rows must be a positive integer")
    if (
        isinstance(accepted_stage_count, bool)
        or not isinstance(accepted_stage_count, int)
        or accepted_stage_count < 0
    ):
        raise ValueError("accepted_stage_count must be a nonnegative integer")
    derivative = SBPFirstDerivative(grid, diagnostic_spatial_order)
    reach = derivative.stencil_reach_intervals
    excluded = fixed_outer_rows + reach
    owned_stop = grid.point_count - excluded
    if owned_stop < 9:
        raise ValueError("RSP1 owned constraint domain is too small")
    full = rsp1_semidiscrete_constraint_snapshot(
        state,
        grid,
        coordinate_time=coordinate_time,
        diagnostic_spatial_order=diagnostic_spatial_order,
        planck_mass=planck_mass,
        length_unit=length_unit,
        point_batch_size=point_batch_size,
    )
    owned_residuals = full.residuals[:owned_stop]
    owned_terms = full.unredefined_term_contributions[:owned_stop]
    owned = normalized_constraint_norms(
        owned_residuals,
        owned_terms,
        planck_mass=planck_mass,
        length_unit=length_unit,
    )
    enclosure = float(
        CONSTRAINT_ROUNDOFF_OPERATION_BUDGET
        * np.finfo(np.float64).eps
        * (1 + accepted_stage_count)
    )
    effective_components = np.maximum(
        owned.component_infinity - enclosure, 0.0
    )
    effective = ConstraintNorms(
        sample_count=owned.sample_count,
        component_names=owned.component_names,
        component_infinity=effective_components,
        global_infinity=float(np.max(effective_components, initial=0.0)),
        minimum_denominator=owned.minimum_denominator,
        maximum_denominator=owned.maximum_denominator,
    )
    mass = _finite("planck_mass", planck_mass, positive=True)
    length = _finite("length_unit", length_unit, positive=True)
    denominator = np.maximum(
        np.sum(np.abs(full.unredefined_term_contributions), axis=2),
        (mass * mass) / (length * length),
    )
    pointwise = np.abs(full.residuals) / denominator
    full_indices = np.argmax(pointwise, axis=0)
    owned_indices = np.argmax(pointwise[:owned_stop], axis=0)
    coordinates = grid.coordinates
    return OwnedConstraintSnapshot(
        point_count=grid.point_count,
        coordinate_time=float(coordinate_time),
        diagnostic_spatial_order=diagnostic_spatial_order,
        accepted_stage_count=accepted_stage_count,
        fixed_outer_rows=fixed_outer_rows,
        derivative_stencil_reach_rows=reach,
        excluded_outer_rows=excluded,
        last_owned_radius=float(coordinates[owned_stop - 1]),
        full_domain_norms=full.raw_norms,
        owned_domain_norms=owned,
        effective_owned_norms_for_order_only=effective,
        accumulated_roundoff_enclosure=enclosure,
        full_component_maximum_radii=coordinates[full_indices],
        owned_component_maximum_radii=coordinates[owned_indices],
    )


@dataclass(frozen=True, slots=True)
class RSP1CommonEventAssessment:
    """One method's complete evidence for the targeted RSP1 question."""

    method: str
    coordinate_time: float
    point_counts: tuple[int, int, int]
    accepted_stage_counts: tuple[int, int, int]
    constraint_admission: OwnedConstraintAdmission
    raw_spatial_spectral_admission: FinestPairNestedSpectralAdmission
    spectral_budgets: tuple[
        Mapping[str, SpectralPowerBudget],
        Mapping[str, SpectralPowerBudget],
        Mapping[str, SpectralPowerBudget],
    ]
    profile_convergence: ThreeGridProfileConvergence
    target_field: str
    target_derivative_tail_ratios: tuple[float | None, float | None]
    target_direct_pair_passed: tuple[bool, bool]
    target_medium_and_fine_absolute_budgets_passed: bool
    target_complete_profile_contraction_passed: bool
    all_medium_and_fine_absolute_budgets_passed: bool
    maximum_spatial_round_trip_interpolation_infinity: float
    study_admission_passed: bool


def rsp1_gr0_common_event(
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
    point_batch_size: int = RSP1_SOURCE_POINT_BATCH_SIZE,
    spectral_thresholds: SpectralThresholds | None = None,
) -> RSP1CommonEventAssessment:
    """Assess only the predeclared amplitude-three spectrum question."""

    if method not in PROTO5_METHOD_CONTRACT:
        raise ValueError("unknown RSP1 GR-0 method")
    records = tuple(states)
    meshes = tuple(grids)
    stages = tuple(accepted_stage_counts)
    if (
        len(records) != 3
        or len(meshes) != 3
        or len(stages) != 3
        or tuple(grid.point_count for grid in meshes) != RSP1_POINT_COUNTS
        or any(
            isinstance(value, bool) or not isinstance(value, int) or value < 0
            for value in stages
        )
    ):
        raise ValueError("RSP1 common event differs from the frozen grid contract")
    if any(
        not isinstance(state, EvolutionState)
        or state.shape != (grid.point_count, FIELD_COUNT)
        for state, grid in zip(records, meshes, strict=True)
    ):
        raise ValueError("RSP1 common-event state/grid pairs differ")
    if fixed_outer_rows != PROTO8_CONSTRAINT_OWNERSHIP[
        "projector_fixed_outer_rows"
    ]:
        raise ValueError("RSP1 fixed outer-row ownership differs")
    batch_size = _point_batch_size(point_batch_size)
    limits = SpectralThresholds() if spectral_thresholds is None else spectral_thresholds
    if not isinstance(limits, SpectralThresholds):
        raise TypeError("spectral_thresholds must be SpectralThresholds")
    fraction = _finite("taper_fraction", taper_fraction, positive=True)
    if fraction >= 0.5:
        raise ValueError("taper_fraction must be below one half")
    time = _finite("coordinate_time", coordinate_time)
    contract = PROTO5_METHOD_CONTRACT[method]
    order = int(contract["spatial_order"])
    snapshots = tuple(
        rsp1_owned_constraint_snapshot(
            state,
            grid,
            coordinate_time=time,
            diagnostic_spatial_order=order,
            fixed_outer_rows=fixed_outer_rows,
            accepted_stage_count=accepted,
            point_batch_size=batch_size,
        )
        for state, grid, accepted in zip(records, meshes, stages, strict=True)
    )
    constraint = owned_constraint_admission(
        snapshots,
        method=method,
        coarsest_guard_maximum=contract["coarsest_constraint_guard"],
        finest_guard_maximum=contract["finest_constraint_guard"],
        minimum_finest_pair_order=1.5,
    )
    profiles = []
    budgets: list[dict[str, SpectralPowerBudget]] = []
    interpolation_maximum = 0.0
    for state, grid in zip(records, meshes, strict=True):
        profile = proper_radial_profile(
            state.u,
            state.q,
            grid.coordinates,
            cutoff=cutoff,
            measurement_radius_maximum=measurement_radius_maximum,
        )
        proper = profile.proper_coordinates
        window = compact_vacuum_buffer_window(
            proper,
            support_minimum=float(proper[0]),
            support_maximum=float(proper[-1]),
            taper_width=float((proper[-1] - proper[0]) * fraction),
        )
        profiles.append(profile)
        budgets.append(
            dict(
                spectral_field_budgets(
                    profile.deviations,
                    proper,
                    cutoff=cutoff,
                    window=window,
                    thresholds=limits,
                )
            )
        )
        interpolation_maximum = max(
            interpolation_maximum,
            float(np.max(profile.round_trip_interpolation_infinity, initial=0.0)),
        )
    raw = finest_pair_nested_spectral_admission(
        RSP1_POINT_COUNTS,
        budgets,
        thresholds=limits,
    )
    convergence = three_grid_profile_convergence(RSP1_POINT_COUNTS, profiles)
    ratios = raw.derivative_power_tail_ratios[RSP1_TARGET_FIELD]
    pair_passed = tuple(
        value is not None and value < limits.maximum_nested_tail_ratio
        for value in ratios
    )
    target_budgets = all(
        budgets[index][RSP1_TARGET_FIELD].individual_admission_passed
        for index in (1, 2)
    )
    all_budgets = all(
        budgets[index][name].individual_admission_passed
        for index in (1, 2)
        for name in PROTO4_SPECTRAL_FIELD_ORDER
    )
    target_profile = convergence.complete_profile_contraction_by_field[
        RSP1_TARGET_FIELD
    ]
    admitted = (
        constraint.admission_passed
        and all_budgets
        and target_budgets
        and target_profile
        and all(pair_passed)
    )
    return RSP1CommonEventAssessment(
        method=method,
        coordinate_time=time,
        point_counts=RSP1_POINT_COUNTS,
        accepted_stage_counts=(int(stages[0]), int(stages[1]), int(stages[2])),
        constraint_admission=constraint,
        raw_spatial_spectral_admission=raw,
        spectral_budgets=(budgets[0], budgets[1], budgets[2]),
        profile_convergence=convergence,
        target_field=RSP1_TARGET_FIELD,
        target_derivative_tail_ratios=ratios,
        target_direct_pair_passed=(bool(pair_passed[0]), bool(pair_passed[1])),
        target_medium_and_fine_absolute_budgets_passed=target_budgets,
        target_complete_profile_contraction_passed=target_profile,
        all_medium_and_fine_absolute_budgets_passed=all_budgets,
        maximum_spatial_round_trip_interpolation_infinity=interpolation_maximum,
        study_admission_passed=admitted,
    )


__all__ = [
    "RSP1ChunkedSourceDiagnosis",
    "RSP1CommonEventAssessment",
    "RSP1GR0EvolutionOperator",
    "RSP1_POINT_COUNTS",
    "RSP1_SOURCE_POINT_BATCH_SIZE",
    "RSP1_TARGET_FIELD",
    "diagnose_gr0_reference_balanced_accelerations_chunked",
    "rsp1_gr0_common_event",
    "rsp1_owned_constraint_snapshot",
    "rsp1_semidiscrete_constraint_snapshot",
]
