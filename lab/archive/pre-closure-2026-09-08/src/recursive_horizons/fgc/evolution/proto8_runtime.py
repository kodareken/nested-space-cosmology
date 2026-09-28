"""Fail-closed PROTO8 common-event composition for GR-0 calibration.

PROTO8 changes no evolution equation, source solve, integrator, grid, physical
input, numerical threshold, or trajectory transaction inherited from PROTO7.
It owns only two diagnostics evaluated after all six members have committed to
one synchronized common event:

* constraint magnitude and convergence are assessed on rows owned by the
  evolution operator, while full-domain residuals remain public and
  accumulated binary64 roundoff may affect order classification only;
* the coarsest spatial spectrum remains a public convergence witness, the
  medium and finest spectra must each pass every unchanged absolute budget,
  and every adjacent nested tail must contract below the unchanged ceiling.

This module performs no evolution and creates no run output.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite
from numbers import Real
from typing import Sequence

import numpy as np

from .cal4_common_event_diagnosis import (
    FinestPairNestedSpectralAdmission,
    OwnedConstraintAdmission,
    finest_pair_nested_spectral_admission,
    owned_constraint_admission,
    owned_semidiscrete_constraint_snapshot,
)
from .health_monitor import compact_vacuum_buffer_window
from .numerical_engine import EvolutionState, UniformRadialGrid
from .proto4_admission import (
    SpectralPowerBudget,
    SpectralThresholds,
    proper_radial_profile,
    spectral_field_budgets,
)
from .proto5_runtime import PROTO5_METHOD_CONTRACT
from .protocol_v8 import PROTO8_CONSTRAINT_OWNERSHIP


@dataclass(frozen=True, slots=True)
class Proto8GR0CommonEventAssessment:
    """Complete PROTO8 admission evidence at one committed common event."""

    method: str
    coordinate_time: float
    point_counts: tuple[int, ...]
    accepted_stage_counts: tuple[int, ...]
    constraint_admission: OwnedConstraintAdmission
    spatial_spectral_admission: FinestPairNestedSpectralAdmission
    maximum_spatial_round_trip_interpolation_infinity: float
    admission_passed: bool


def _finite(name: str, value: Real, *, positive: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a real scalar")
    answer = float(value)
    if not isfinite(answer):
        raise ValueError(f"{name} must be finite")
    if positive and answer <= 0.0:
        raise ValueError(f"{name} must be positive")
    return answer


def _stage_counts(values: Sequence[int]) -> tuple[int, ...]:
    answer = tuple(values)
    if len(answer) != 3 or any(
        isinstance(item, bool) or not isinstance(item, int) or item < 0
        for item in answer
    ):
        raise ValueError("accepted_stage_counts must contain three nonnegative integers")
    return answer


def proto8_gr0_common_event(
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
) -> Proto8GR0CommonEventAssessment:
    """Compose the frozen PROTO8 constraints and spatial spectra.

    The caller must provide three already committed, synchronized grid states.
    Accepted-stage counts are resolution-local and enter only the declared
    roundoff enclosure used for constraint-order classification.
    """

    if method not in PROTO5_METHOD_CONTRACT:
        raise ValueError("unknown PROTO8 GR-0 method")
    records = tuple(states)
    meshes = tuple(grids)
    stages = _stage_counts(accepted_stage_counts)
    if len(records) != 3 or len(meshes) != 3:
        raise ValueError("PROTO8 common events require exactly three grids")
    if any(
        not isinstance(state, EvolutionState)
        or not isinstance(grid, UniformRadialGrid)
        or state.shape != (grid.point_count, 6)
        for state, grid in zip(records, meshes, strict=True)
    ):
        raise ValueError("PROTO8 common-event state/grid pairs differ")
    counts = tuple(grid.point_count for grid in meshes)
    if counts != tuple(sorted(counts)) or len(set(counts)) != 3:
        raise ValueError("PROTO8 common-event grids must increase strictly")
    frozen_fixed = PROTO8_CONSTRAINT_OWNERSHIP["projector_fixed_outer_rows"]
    if fixed_outer_rows != frozen_fixed:
        raise ValueError("PROTO8 fixed outer-row ownership differs")

    contract = PROTO5_METHOD_CONTRACT[method]
    order = int(contract["spatial_order"])
    time = _finite("coordinate_time", coordinate_time)
    snapshots = tuple(
        owned_semidiscrete_constraint_snapshot(
            state,
            grid,
            coordinate_time=time,
            diagnostic_spatial_order=order,
            fixed_outer_rows=fixed_outer_rows,
            accepted_stage_count=accepted,
        )
        for state, grid, accepted in zip(records, meshes, stages, strict=True)
    )
    expected_reach = (
        PROTO8_CONSTRAINT_OWNERSHIP["RK4_derivative_stencil_reach_rows"]
        if method == "RK4"
        else PROTO8_CONSTRAINT_OWNERSHIP["SSPRK3_derivative_stencil_reach_rows"]
    )
    if any(item.derivative_stencil_reach_rows != expected_reach for item in snapshots):
        raise RuntimeError("PROTO8 derivative stencil ownership differs")
    constraint = owned_constraint_admission(
        snapshots,
        method=method,
        coarsest_guard_maximum=contract["coarsest_constraint_guard"],
        finest_guard_maximum=contract["finest_constraint_guard"],
        minimum_finest_pair_order=1.5,
    )

    fraction = _finite("taper_fraction", taper_fraction, positive=True)
    if fraction >= 0.5:
        raise ValueError("taper_fraction must be below one half")
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
        budgets.append(
            dict(
                spectral_field_budgets(
                    profile.deviations,
                    proper,
                    cutoff=cutoff,
                    window=window,
                    thresholds=spectral_thresholds,
                )
            )
        )
        interpolation_maximum = max(
            interpolation_maximum,
            float(np.max(profile.round_trip_interpolation_infinity, initial=0.0)),
        )
    spatial = finest_pair_nested_spectral_admission(
        counts,
        budgets,
        thresholds=spectral_thresholds,
    )
    return Proto8GR0CommonEventAssessment(
        method=method,
        coordinate_time=time,
        point_counts=counts,
        accepted_stage_counts=stages,
        constraint_admission=constraint,
        spatial_spectral_admission=spatial,
        maximum_spatial_round_trip_interpolation_infinity=interpolation_maximum,
        admission_passed=(constraint.admission_passed and spatial.admission_passed),
    )


__all__ = ["Proto8GR0CommonEventAssessment", "proto8_gr0_common_event"]
