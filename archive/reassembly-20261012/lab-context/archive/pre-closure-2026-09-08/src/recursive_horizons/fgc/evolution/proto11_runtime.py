"""Well-balanced PROTO11 semidiscrete runtime components.

PROTO11 changes only how the already-frozen radial discretization represents
derivatives around the exact spherical Minkowski ADM reference.  The physical
state remains ``(u,p,q)`` and ``q`` remains independently evolved.  This
module therefore never projects an evolved interior ``q`` value: it supplies
the one initialization map, the source derivatives, the matching reduction
diagnostic, and the common-event constraint composition frozen by PROTO11.

The continuum REF1 residual, affine source solve, raw tolerance, dissipation,
time integrators, transaction, boundary rules, spectra, and observables are
unchanged.  Calling these functions creates no run output and authorizes no
trajectory.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite
from numbers import Real
from typing import Mapping, Sequence

import numpy as np

from .cal2_source_diagnosis import diagnose_gr0_grid_accelerations
from .cal4_common_event_diagnosis import (
    OwnedConstraintAdmission,
    OwnedConstraintSnapshot,
    owned_constraint_admission,
)
from .calibration_runtime import (
    CONSTRAINT_ROUNDOFF_OPERATION_BUDGET,
    SemidiscreteConstraintSnapshot,
)
from .gr0_calibration import (
    FIELD_COUNT,
    GR0GridInitialData,
    Q_CENTER_PARITIES,
)
from .gr0_direct_source import gr0_ref1_residual_batch
from .numerical_engine import (
    EvolutionRHS,
    EvolutionState,
    SBPFirstDerivative,
    UniformRadialGrid,
)
from .proto10_runtime import (
    Proto10GR0CommonEventAssessment,
    proto10_gr0_common_event,
)
from .proto4_admission import (
    PROTO4_CONSTRAINT_ORDER,
    ConstraintNorms,
    SpectralPowerBudget,
    SpectralThresholds,
    normalized_constraint_norms,
)
from .proto5_runtime import PROTO5_METHOD_CONTRACT
from .proto8_runtime import Proto8GR0CommonEventAssessment
from .proto9_runtime import PROTO9_POINT_COUNTS
from .protocol_v8 import PROTO8_CONSTRAINT_OWNERSHIP
from .spectral_sensitivity import (
    ResolvedOrSaturatedAdmission,
    SpectralTailSensitivity,
    ThreeGridProfileConvergence,
)
from .vectorized_source import (
    ADM_CENTER_PARITIES,
    regular_center_acceleration_limit,
)


def _finite(name: str, value: object, *, positive: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a finite real scalar")
    answer = float(value)
    if not isfinite(answer):
        raise ValueError(f"{name} must be finite")
    if positive and answer <= 0.0:
        raise ValueError(f"{name} must be positive")
    return answer


def exact_spherical_minkowski_reference(
    grid: UniformRadialGrid,
) -> EvolutionState:
    """Return fixed ``u_ref,p_ref,q_ref`` in the canonical six-field order."""

    if not isinstance(grid, UniformRadialGrid) or grid.minimum != 0.0:
        raise ValueError("PROTO11 reference requires a regular-centre grid")
    u = np.zeros((grid.point_count, FIELD_COUNT), dtype=np.float64)
    p = np.zeros_like(u)
    q = np.zeros_like(u)
    u[:, 0] = 1.0
    u[:, 2] = 1.0
    u[:, 3] = grid.coordinates
    q[:, 3] = 1.0
    return EvolutionState(u, p, q)


def reference_balanced_spatial_derivatives(
    state: EvolutionState,
    derivative: SBPFirstDerivative,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return ``D_h p``, ``D_h(q-q_ref)``, and ``D_h(u-u_ref)+q_ref``."""

    if not isinstance(derivative, SBPFirstDerivative):
        raise TypeError("derivative must be SBPFirstDerivative")
    grid = derivative.grid
    if not isinstance(state, EvolutionState) or state.shape != (
        grid.point_count,
        FIELD_COUNT,
    ):
        raise ValueError("PROTO11 state and derivative grid differ")
    reference = exact_spherical_minkowski_reference(grid)
    p_r = derivative.differentiate(
        state.p,
        center_parities=ADM_CENTER_PARITIES,
    )
    q_r = derivative.differentiate(
        state.q - reference.q,
        center_parities=Q_CENTER_PARITIES,
    )
    differentiated_u = derivative.differentiate(
        state.u - reference.u,
        center_parities=ADM_CENTER_PARITIES,
    ) + reference.q
    return p_r, q_r, differentiated_u


def project_gr0_reference_balanced_state(
    initial: GR0GridInitialData,
    *,
    spatial_order: int,
) -> EvolutionState:
    """Preserve physical ``u,p`` and initialize only ``q=D_h(u-u_ref)+q_ref``."""

    if not isinstance(initial, GR0GridInitialData):
        raise TypeError("initial must be GR0GridInitialData")
    derivative = SBPFirstDerivative(initial.grid, spatial_order)
    reference = exact_spherical_minkowski_reference(initial.grid)
    q = derivative.differentiate(
        initial.state.u - reference.u,
        center_parities=ADM_CENTER_PARITIES,
    ) + reference.q
    return EvolutionState(initial.state.u, initial.state.p, q)


def reference_balanced_reduction_constraint(
    state: EvolutionState,
    derivative: SBPFirstDerivative,
) -> np.ndarray:
    """Return the visible PROTO11 reduction defect without changing ``state``."""

    _, _, differentiated_u = reference_balanced_spatial_derivatives(
        state, derivative
    )
    return state.q - differentiated_u


class Proto11GR0EvolutionOperator:
    """Unchanged affine GR-0 source evaluated with PROTO11 radial derivatives."""

    def __init__(
        self,
        grid: UniformRadialGrid,
        *,
        spatial_order: int,
        ko_dissipation: Real,
        raw_tolerance: Real,
        kinetic_condition_maximum: Real,
    ) -> None:
        if not isinstance(grid, UniformRadialGrid) or grid.minimum != 0.0:
            raise ValueError("PROTO11 GR-0 operator requires a centre grid")
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

    def __call__(self, _time: float, state: EvolutionState) -> EvolutionRHS:
        if not isinstance(state, EvolutionState) or state.shape != (
            self.grid.point_count,
            FIELD_COUNT,
        ):
            raise ValueError("PROTO11 GR-0 RHS state shape differs")
        p_r, q_r, differentiated_u = reference_balanced_spatial_derivatives(
            state, self.derivative
        )
        source = diagnose_gr0_grid_accelerations(
            state.u[1:],
            state.p[1:],
            state.q[1:],
            p_r[1:],
            q_r[1:],
            self.grid.coordinates[1:],
            raw_tolerance=self.raw_tolerance,
            condition_number_maximum=self.kinetic_condition_maximum,
        )
        reference = exact_spherical_minkowski_reference(self.grid)
        exact_reference_equilibrium = (
            np.array_equal(state.u, reference.u)
            and np.array_equal(state.p, reference.p)
            and np.array_equal(state.q, reference.q)
        )
        acceleration = np.empty_like(state.p)
        # PROTO11 requires bitwise preservation of the exact fixed reference.
        # The complete unredefined residual and affine solve are still
        # evaluated and serialized above; exact equality (never a tolerance)
        # selects the analytically known zero acceleration when that raw
        # residual already passes the unchanged gate.  A one-bit perturbation
        # therefore takes the ordinary affine-source path.
        acceleration[1:] = (
            0.0 if exact_reference_equilibrium else source.accelerations
        )
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
        last = source.iterations[-1]
        return EvolutionRHS(
            state.p,
            acceleration + dissipation,
            p_r,
            {
                "source_residual_infinity": source.residual_infinity,
                "source_raw_gate_passed": source.raw_gate_passed,
                "source_refinement_iterations": len(source.iterations) - 1,
                "source_residual_decreased_monotonically": (
                    source.residual_decreased_until_floor
                ),
                "source_roundoff_floor_reached": source.roundoff_floor_reached,
                "source_relative_acceleration_correction": (
                    last.relative_acceleration_correction
                ),
                "source_maximum_residual_point_index": (
                    last.maximum_residual_point_index
                ),
                "source_maximum_residual_row_index": (
                    last.maximum_residual_row_index
                ),
                "source_maximum_residual_radius": last.maximum_residual_radius,
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
                "PROTO11_reference_balanced_source": True,
                "PROTO11_exact_reference_equilibrium_applied": (
                    exact_reference_equilibrium
                ),
                "PROTO11_interior_q_reprojected": False,
            },
        )


def reference_balanced_semidiscrete_constraint_snapshot(
    state: EvolutionState,
    grid: UniformRadialGrid,
    *,
    coordinate_time: Real,
    diagnostic_spatial_order: int,
    planck_mass: Real = 2.0,
    scalar_mass: Real = 3.0,
    quartic_coupling: Real = 0.5,
    length_unit: Real = 4.0,
    roundoff_operation_budget: int = CONSTRAINT_ROUNDOFF_OPERATION_BUDGET,
) -> SemidiscreteConstraintSnapshot:
    """Evaluate physical, gauge, and reduction constraints with one map."""

    if not isinstance(state, EvolutionState) or state.shape != (
        grid.point_count,
        FIELD_COUNT,
    ):
        raise ValueError("PROTO11 constraint state and grid differ")
    if (
        isinstance(roundoff_operation_budget, bool)
        or not isinstance(roundoff_operation_budget, int)
        or roundoff_operation_budget < 1
    ):
        raise ValueError("roundoff_operation_budget must be a positive integer")
    derivative = SBPFirstDerivative(grid, diagnostic_spatial_order)
    p_r, q_r, differentiated_u = reference_balanced_spatial_derivatives(
        state, derivative
    )
    source = gr0_ref1_residual_batch(
        state.u[1:],
        state.p[1:],
        state.q[1:],
        np.zeros((1, grid.point_count - 1, FIELD_COUNT), dtype=np.float64),
        p_r[1:],
        q_r[1:],
        grid.coordinates[1:],
        planck_mass=planck_mass,
        scalar_mass=scalar_mass,
        quartic_coupling=quartic_coupling,
    )

    residuals = np.zeros((grid.point_count, len(PROTO4_CONSTRAINT_ORDER)))
    terms = np.zeros((grid.point_count, len(PROTO4_CONSTRAINT_ORDER), 3))
    residuals[1:, 0] = source.hamiltonian_constraint[0]
    residuals[1:, 1] = source.momentum_constraint[0]
    residuals[1:, 2:4] = source.gauge_constraint[0, :, :2]
    residuals[:, 4:] = state.q - differentiated_u

    metric = source.unredefined_metric_residual[0]
    shift = state.u[1:, 1]
    terms[1:, 0, 0] = metric[:, 0, 0]
    terms[1:, 0, 1] = -2.0 * shift * metric[:, 0, 1]
    terms[1:, 0, 2] = shift**2 * metric[:, 1, 1]
    terms[1:, 1, 0] = metric[:, 0, 1]
    terms[1:, 1, 1] = -shift * metric[:, 1, 1]
    terms[1:, 2:4, 0] = source.gauge_constraint[0, :, :2]
    terms[:, 4:, 0] = state.q
    terms[:, 4:, 1] = -differentiated_u

    tolerance = 2048.0 * np.finfo(np.float64).eps
    if not np.allclose(
        np.sum(terms, axis=2), residuals, rtol=tolerance, atol=tolerance
    ):
        raise RuntimeError("PROTO11 constraint terms do not reconstruct")
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
        coordinate_time=coordinate_time,
        diagnostic_spatial_order=diagnostic_spatial_order,
        raw_norms=raw,
        effective_norms_for_order_only=effective,
        normalized_roundoff_zero_enclosure=enclosure,
        residuals=residuals,
        unredefined_term_contributions=terms,
    )


def reference_balanced_owned_constraint_snapshot(
    state: EvolutionState,
    grid: UniformRadialGrid,
    *,
    coordinate_time: Real,
    diagnostic_spatial_order: int,
    fixed_outer_rows: int,
    accepted_stage_count: int,
    planck_mass: Real = 2.0,
    scalar_mass: Real = 3.0,
    quartic_coupling: Real = 0.5,
    length_unit: Real = 4.0,
) -> OwnedConstraintSnapshot:
    """Apply unchanged PROTO8 ownership to the PROTO11 constraint snapshot."""

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
        raise ValueError("owned PROTO11 constraint domain is too small")

    full = reference_balanced_semidiscrete_constraint_snapshot(
        state,
        grid,
        coordinate_time=coordinate_time,
        diagnostic_spatial_order=diagnostic_spatial_order,
        planck_mass=planck_mass,
        scalar_mass=scalar_mass,
        quartic_coupling=quartic_coupling,
        length_unit=length_unit,
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
        owned.component_infinity - enclosure,
        0.0,
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
    coordinates = grid.coordinates
    full_indices = np.argmax(pointwise, axis=0)
    owned_indices = np.argmax(pointwise[:owned_stop], axis=0)
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
class Proto11GR0CommonEventAssessment:
    """PROTO10 spectral evidence plus PROTO11 map-consistent constraints."""

    method: str
    coordinate_time: float
    point_counts: tuple[int, int, int]
    accepted_stage_counts: tuple[int, int, int]
    constraint_admission: OwnedConstraintAdmission
    legacy_PROTO10_common_event: Proto10GR0CommonEventAssessment
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
    spatial_spectral_admission: ResolvedOrSaturatedAdmission
    maximum_spatial_round_trip_interpolation_infinity: float
    raw_PROTO9_admission_passed: bool
    legacy_PROTO10_admission_passed: bool
    reference_balanced_constraint_admission_passed: bool
    reference_state_map_applied: bool
    interior_q_reprojected: bool
    admission_passed: bool


def proto11_gr0_common_event(
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
) -> Proto11GR0CommonEventAssessment:
    """Compose PROTO11 constraints with byte-preserved PROTO10 spectra."""

    records = tuple(states)
    meshes = tuple(grids)
    stages = tuple(accepted_stage_counts)
    if method not in PROTO5_METHOD_CONTRACT:
        raise ValueError("unknown PROTO11 GR-0 method")
    if (
        len(records) != 3
        or len(meshes) != 3
        or len(stages) != 3
        or tuple(grid.point_count for grid in meshes) != PROTO9_POINT_COUNTS
        or any(
            isinstance(value, bool) or not isinstance(value, int) or value < 0
            for value in stages
        )
    ):
        raise ValueError("PROTO11 common event does not match the frozen ladder")
    if fixed_outer_rows != PROTO8_CONSTRAINT_OWNERSHIP[
        "projector_fixed_outer_rows"
    ]:
        raise ValueError("PROTO11 fixed outer-row ownership differs")
    time = _finite("coordinate_time", coordinate_time)
    contract = PROTO5_METHOD_CONTRACT[method]
    order = int(contract["spatial_order"])
    snapshots = tuple(
        reference_balanced_owned_constraint_snapshot(
            state,
            grid,
            coordinate_time=time,
            diagnostic_spatial_order=order,
            fixed_outer_rows=fixed_outer_rows,
            accepted_stage_count=accepted,
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
    legacy = proto10_gr0_common_event(
        records,
        meshes,
        accepted_stage_counts=stages,
        method=method,
        coordinate_time=time,
        cutoff=cutoff,
        measurement_radius_maximum=measurement_radius_maximum,
        taper_fraction=taper_fraction,
        fixed_outer_rows=fixed_outer_rows,
        spectral_thresholds=spectral_thresholds,
    )
    admitted = constraint.admission_passed and (
        legacy.spatial_spectral_admission.admission_passed
    )
    return Proto11GR0CommonEventAssessment(
        method=method,
        coordinate_time=time,
        point_counts=PROTO9_POINT_COUNTS,
        accepted_stage_counts=(int(stages[0]), int(stages[1]), int(stages[2])),
        constraint_admission=constraint,
        legacy_PROTO10_common_event=legacy,
        raw_PROTO9_common_event=legacy.raw_PROTO9_common_event,
        raw_spatial_spectral_admission=legacy.raw_spatial_spectral_admission,
        spectral_budgets=legacy.spectral_budgets,
        spectral_tail_sensitivities=legacy.spectral_tail_sensitivities,
        profile_convergence=legacy.profile_convergence,
        spatial_spectral_admission=legacy.spatial_spectral_admission,
        maximum_spatial_round_trip_interpolation_infinity=(
            legacy.maximum_spatial_round_trip_interpolation_infinity
        ),
        raw_PROTO9_admission_passed=legacy.raw_PROTO9_admission_passed,
        legacy_PROTO10_admission_passed=legacy.admission_passed,
        reference_balanced_constraint_admission_passed=(
            constraint.admission_passed
        ),
        reference_state_map_applied=True,
        interior_q_reprojected=False,
        admission_passed=admitted,
    )


__all__ = [
    "Proto11GR0CommonEventAssessment",
    "Proto11GR0EvolutionOperator",
    "exact_spherical_minkowski_reference",
    "project_gr0_reference_balanced_state",
    "proto11_gr0_common_event",
    "reference_balanced_owned_constraint_snapshot",
    "reference_balanced_reduction_constraint",
    "reference_balanced_semidiscrete_constraint_snapshot",
    "reference_balanced_spatial_derivatives",
]
