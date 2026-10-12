"""Resolution-only PROTO9 adapter for committed GR-0 common events.

PROTO9 inherits the complete PROTO8 diagnostic compositor and every evolution,
transaction, physical, and numerical rule.  Its only runtime responsibility is
to fail closed unless the synchronized nested ladder is exactly
``2049 -> 4097 -> 8193``.  In particular, this adapter cannot admit a
conditional projection or relabel a PROTO8 state as new finest-grid evidence.
"""

from __future__ import annotations

from numbers import Real
from typing import Sequence

from .numerical_engine import EvolutionState, UniformRadialGrid
from .proto4_admission import SpectralThresholds
from .proto8_runtime import (
    Proto8GR0CommonEventAssessment,
    proto8_gr0_common_event,
)


PROTO9_POINT_COUNTS = (2049, 4097, 8193)


def validate_proto9_common_event_assessment(
    assessment: Proto8GR0CommonEventAssessment,
) -> Proto8GR0CommonEventAssessment:
    """Reject an inherited assessment unless it carries the frozen ladder."""

    if not isinstance(assessment, Proto8GR0CommonEventAssessment):
        raise TypeError("PROTO9 requires a PROTO8 common-event assessment")
    if assessment.point_counts != PROTO9_POINT_COUNTS:
        raise ValueError("PROTO9 assessment does not carry new 8193-point evidence")
    return assessment


def proto9_gr0_common_event(
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
    """Apply unchanged PROTO8 admission on the exact PROTO9 grid ladder."""

    records = tuple(states)
    meshes = tuple(grids)
    counts = tuple(
        grid.point_count if isinstance(grid, UniformRadialGrid) else None
        for grid in meshes
    )
    if len(records) != 3 or len(meshes) != 3 or counts != PROTO9_POINT_COUNTS:
        raise ValueError(
            "PROTO9 common events require exact point counts "
            f"{PROTO9_POINT_COUNTS!r}"
        )
    assessment = proto8_gr0_common_event(
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
    return validate_proto9_common_event_assessment(assessment)


__all__ = [
    "PROTO9_POINT_COUNTS",
    "proto9_gr0_common_event",
    "validate_proto9_common_event_assessment",
]
