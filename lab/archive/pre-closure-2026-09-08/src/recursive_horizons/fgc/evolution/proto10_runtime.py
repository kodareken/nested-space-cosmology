"""Fail-closed PROTO10 common-event spectral composition.

PROTO10 preserves the complete PROTO9 common-event assessment and adds one
conditioning-guarded interpretation of the finest nested spectral pair.  The
raw ratios, absolute budgets, constraint assessment, grid ladder, and direct
``< 1/4`` route remain public.  A failed finest-pair ratio can be classified
as diagnostically saturated only by :mod:`spectral_sensitivity`'s complete
tail-erasure, map-residual, and whole-profile guard set.

This module performs no evolution, reads no run output, and does not interpret
diagnostic saturation as a continuum-error bound or physical resolution.
"""

from __future__ import annotations

from dataclasses import dataclass
from numbers import Real
from typing import Mapping, Sequence

import numpy as np

from .cal4_common_event_diagnosis import (
    FinestPairNestedSpectralAdmission,
    OwnedConstraintAdmission,
)
from .health_monitor import compact_vacuum_buffer_window
from .numerical_engine import EvolutionState, UniformRadialGrid
from .proto4_admission import (
    PROTO4_SPECTRAL_FIELD_ORDER,
    SpectralPowerBudget,
    SpectralThresholds,
    proper_radial_profile,
    spectral_field_budgets,
)
from .proto8_runtime import Proto8GR0CommonEventAssessment
from .proto9_runtime import PROTO9_POINT_COUNTS, proto9_gr0_common_event
from .spectral_sensitivity import (
    ResolvedOrSaturatedAdmission,
    SpectralTailSensitivity,
    ThreeGridProfileConvergence,
    resolved_or_saturated_admission,
    spectral_tail_sensitivity,
    three_grid_profile_convergence,
)


@dataclass(frozen=True, slots=True)
class Proto10GR0CommonEventAssessment:
    """Raw PROTO9 evidence plus PROTO10's guarded spectral interpretation."""

    method: str
    coordinate_time: float
    point_counts: tuple[int, int, int]
    accepted_stage_counts: tuple[int, int, int]
    constraint_admission: OwnedConstraintAdmission
    raw_PROTO9_common_event: Proto8GR0CommonEventAssessment
    raw_spatial_spectral_admission: FinestPairNestedSpectralAdmission
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
    admission_passed: bool


def _raw_ratios_equal(
    raw: FinestPairNestedSpectralAdmission,
    guarded: ResolvedOrSaturatedAdmission,
) -> bool:
    return (
        dict(raw.field_power_tail_ratios)
        == dict(guarded.field_power_tail_ratios)
        and dict(raw.derivative_power_tail_ratios)
        == dict(guarded.derivative_power_tail_ratios)
    )


def proto10_gr0_common_event(
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
) -> Proto10GR0CommonEventAssessment:
    """Recompute and compose every PROTO10 guard on one committed event.

    The immutable PROTO9 compositor is evaluated first.  PROTO10 then rebuilds
    the same proper-distance profiles and budgets, measures the diagnostic-map
    sensitivity, and applies the guarded finest-pair rule.  Internal equality
    checks make any change to the raw ratios or absolute-budget decisions a
    runtime error rather than an alternative interpretation.
    """

    records = tuple(states)
    meshes = tuple(grids)
    raw = proto9_gr0_common_event(
        records,
        meshes,
        accepted_stage_counts=accepted_stage_counts,
        method=method,
        coordinate_time=coordinate_time,
        cutoff=cutoff,
        measurement_radius_maximum=measurement_radius_maximum,
        taper_fraction=taper_fraction,
        fixed_outer_rows=fixed_outer_rows,
        spectral_thresholds=spectral_thresholds,
    )
    if raw.point_counts != PROTO9_POINT_COUNTS:
        raise RuntimeError("PROTO10 lost the exact PROTO9 resolution ladder")

    fraction = float(taper_fraction)
    if not np.isfinite(fraction) or fraction <= 0.0 or fraction >= 0.5:
        raise ValueError("taper_fraction must be finite and lie in (0, 1/2)")

    profiles = []
    budgets: list[dict[str, SpectralPowerBudget]] = []
    sensitivities: list[dict[str, SpectralTailSensitivity]] = []
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
        budget = dict(
            spectral_field_budgets(
                profile.deviations,
                proper,
                cutoff=cutoff,
                window=window,
                thresholds=spectral_thresholds,
            )
        )
        witness = {
            name: spectral_tail_sensitivity(
                profile.deviations[:, index],
                proper,
                window=window,
                round_trip_interpolation_infinity=(
                    profile.round_trip_interpolation_infinity[index]
                ),
                budget=budget[name],
            )
            for index, name in enumerate(PROTO4_SPECTRAL_FIELD_ORDER)
        }
        profiles.append(profile)
        budgets.append(budget)
        sensitivities.append(witness)

    convergence = three_grid_profile_convergence(PROTO9_POINT_COUNTS, profiles)
    guarded = resolved_or_saturated_admission(
        PROTO9_POINT_COUNTS,
        budgets,
        sensitivities,
        convergence,
        thresholds=spectral_thresholds,
    )

    recomputed_by_grid = {
        point_count: {
            name: bool(budget[name].individual_admission_passed)
            for name in PROTO4_SPECTRAL_FIELD_ORDER
        }
        for point_count, budget in zip(PROTO9_POINT_COUNTS, budgets, strict=True)
    }
    if raw.spatial_spectral_admission.individual_budget_passed_by_grid != recomputed_by_grid:
        raise RuntimeError("PROTO10 recomputation changed a raw absolute budget decision")
    if not _raw_ratios_equal(raw.spatial_spectral_admission, guarded):
        raise RuntimeError("PROTO10 recomputation changed a raw nested-tail ratio")
    if (
        raw.spatial_spectral_admission.finest_pair_individual_budgets_passed
        is not guarded.finest_pair_individual_budgets_passed
    ):
        raise RuntimeError("PROTO10 changed the finest-pair absolute-budget gate")

    admitted = raw.constraint_admission.admission_passed and guarded.admission_passed
    return Proto10GR0CommonEventAssessment(
        method=raw.method,
        coordinate_time=raw.coordinate_time,
        point_counts=(
            int(raw.point_counts[0]),
            int(raw.point_counts[1]),
            int(raw.point_counts[2]),
        ),
        accepted_stage_counts=(
            int(raw.accepted_stage_counts[0]),
            int(raw.accepted_stage_counts[1]),
            int(raw.accepted_stage_counts[2]),
        ),
        constraint_admission=raw.constraint_admission,
        raw_PROTO9_common_event=raw,
        raw_spatial_spectral_admission=raw.spatial_spectral_admission,
        spectral_budgets=(budgets[0], budgets[1], budgets[2]),
        spectral_tail_sensitivities=(
            sensitivities[0],
            sensitivities[1],
            sensitivities[2],
        ),
        profile_convergence=convergence,
        spatial_spectral_admission=guarded,
        maximum_spatial_round_trip_interpolation_infinity=(
            raw.maximum_spatial_round_trip_interpolation_infinity
        ),
        raw_PROTO9_admission_passed=raw.admission_passed,
        admission_passed=admitted,
    )


__all__ = [
    "Proto10GR0CommonEventAssessment",
    "proto10_gr0_common_event",
]
