from __future__ import annotations

from dataclasses import replace
from pathlib import Path
import sys
import unittest

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from recursive_horizons.fgc.evolution.cal4_common_event_diagnosis import (  # noqa: E402
    finest_pair_nested_spectral_admission,
    owned_semidiscrete_constraint_snapshot,
)
from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    EvolutionState,
    SBPFirstDerivative,
    UniformRadialGrid,
)
from recursive_horizons.fgc.evolution.proto4_admission import (  # noqa: E402
    PROTO4_CONSTRAINT_ORDER,
    SpectralPowerBudget,
)
from recursive_horizons.fgc.evolution.vectorized_source import (  # noqa: E402
    ADM_CENTER_PARITIES,
)


def _flat_state(grid: UniformRadialGrid, order: int) -> EvolutionState:
    u = np.zeros((grid.point_count, 6), dtype=np.float64)
    u[:, 0] = 1.0
    u[:, 2] = 1.0
    u[:, 3] = grid.coordinates
    p = np.zeros_like(u)
    q = SBPFirstDerivative(grid, order).differentiate(
        u,
        center_parities=ADM_CENTER_PARITIES,
    )
    return EvolutionState(u, p, q)


def _budget(
    sample_count: int,
    *,
    field_fraction: float,
    derivative_fraction: float,
    admitted: bool,
) -> SpectralPowerBudget:
    nyquist = sample_count // 2
    return SpectralPowerBudget(
        sample_count=sample_count,
        top_band_first_bin=max(0, nyquist - 2),
        nyquist_bin=nyquist,
        total_field_power=1.0,
        top_band_field_power_fraction=field_fraction,
        total_derivative_weighted_power=1.0,
        top_band_derivative_weighted_power_fraction=derivative_fraction,
        rms_angular_scale=1.0,
        rms_scale_over_cutoff=0.1,
        last_occupied_bin=nyquist,
        last_bin_alias_diagnostic=True,
        individual_admission_passed=admitted,
    )


class FGCCAL4CommonEventDiagnosisTests(unittest.TestCase):
    def test_projector_owned_outer_defect_remains_public_but_is_not_owned(self) -> None:
        grid = UniformRadialGrid(0.0, 8.0, 65)
        state = _flat_state(grid, 4)
        q = state.q.copy()
        q[-5, 0] += 1.0e-4
        diagnosed = owned_semidiscrete_constraint_snapshot(
            EvolutionState(state.u, state.p, q),
            grid,
            coordinate_time=0.25,
            diagnostic_spatial_order=4,
            fixed_outer_rows=4,
            accepted_stage_count=20,
        )
        component = PROTO4_CONSTRAINT_ORDER.index("reduction_alpha")
        self.assertEqual(diagnosed.excluded_outer_rows, 7)
        self.assertGreater(
            diagnosed.full_domain_norms.component_infinity[component],
            1.0e-6,
        )
        self.assertLess(
            diagnosed.owned_domain_norms.component_infinity[component],
            diagnosed.accumulated_roundoff_enclosure,
        )
        self.assertGreater(
            diagnosed.full_component_maximum_radii[component],
            diagnosed.last_owned_radius,
        )

    def test_interior_reduction_defect_cannot_hide_in_roundoff_enclosure(self) -> None:
        grid = UniformRadialGrid(0.0, 8.0, 65)
        state = _flat_state(grid, 4)
        q = state.q.copy()
        q[20, 0] += 1.0e-6
        diagnosed = owned_semidiscrete_constraint_snapshot(
            EvolutionState(state.u, state.p, q),
            grid,
            coordinate_time=0.25,
            diagnostic_spatial_order=4,
            fixed_outer_rows=4,
            accepted_stage_count=20,
        )
        component = PROTO4_CONSTRAINT_ORDER.index("reduction_alpha")
        self.assertGreater(
            diagnosed.owned_domain_norms.component_infinity[component],
            1000.0 * diagnosed.accumulated_roundoff_enclosure,
        )
        self.assertLessEqual(
            diagnosed.owned_component_maximum_radii[component],
            diagnosed.last_owned_radius,
        )

    def test_coarse_only_spectral_miss_passes_only_with_resolved_finest_pair(self) -> None:
        counts = (33, 65, 129)
        budgets = (
            {"field": _budget(33, field_fraction=1.0e-4, derivative_fraction=1.0e-2, admitted=False)},
            {"field": _budget(65, field_fraction=1.0e-6, derivative_fraction=1.0e-4, admitted=True)},
            {"field": _budget(129, field_fraction=1.0e-8, derivative_fraction=1.0e-6, admitted=True)},
        )
        result = finest_pair_nested_spectral_admission(counts, budgets)
        self.assertFalse(result.coarsest_individual_budget_passed)
        self.assertTrue(result.finest_pair_individual_budgets_passed)
        self.assertTrue(result.every_nested_tail_ratio_passed)
        self.assertTrue(result.admission_passed)

        unresolved_medium = list(budgets)
        unresolved_medium[1] = {
            "field": replace(unresolved_medium[1]["field"], individual_admission_passed=False)
        }
        self.assertFalse(
            finest_pair_nested_spectral_admission(counts, unresolved_medium).admission_passed
        )

    def test_noncontracting_nested_tail_remains_a_veto(self) -> None:
        counts = (33, 65, 129)
        budgets = (
            {"field": _budget(33, field_fraction=1.0e-4, derivative_fraction=1.0e-2, admitted=False)},
            {"field": _budget(65, field_fraction=1.0e-6, derivative_fraction=1.0e-4, admitted=True)},
            {"field": _budget(129, field_fraction=5.0e-7, derivative_fraction=5.0e-5, admitted=True)},
        )
        result = finest_pair_nested_spectral_admission(counts, budgets)
        self.assertTrue(result.finest_pair_individual_budgets_passed)
        self.assertFalse(result.every_nested_tail_ratio_passed)
        self.assertFalse(result.admission_passed)


if __name__ == "__main__":
    unittest.main()
