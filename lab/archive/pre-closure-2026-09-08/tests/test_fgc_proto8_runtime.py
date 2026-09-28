from __future__ import annotations

from pathlib import Path
import sys
import unittest

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    EvolutionState,
    SBPFirstDerivative,
    UniformRadialGrid,
)
from recursive_horizons.fgc.evolution.proto8_runtime import (  # noqa: E402
    proto8_gr0_common_event,
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
        u, center_parities=ADM_CENTER_PARITIES
    )
    return EvolutionState(u, p, q)


class FGCProto8RuntimeTests(unittest.TestCase):
    def test_both_methods_apply_the_frozen_owned_rows_and_nested_roles(self) -> None:
        grids = tuple(UniformRadialGrid(0.0, 8.0, count) for count in (33, 65, 129))
        for method, order, excluded in (("RK4", 4, 7), ("SSPRK3", 2, 5)):
            with self.subTest(method=method):
                states = tuple(_flat_state(grid, order) for grid in grids)
                assessed = proto8_gr0_common_event(
                    states,
                    grids,
                    accepted_stage_counts=(4, 8, 16),
                    method=method,
                    coordinate_time=0.25,
                    measurement_radius_maximum=8.0,
                )
                self.assertEqual(
                    assessed.constraint_admission.excluded_outer_rows,
                    (excluded, excluded, excluded),
                )
                self.assertEqual(assessed.accepted_stage_counts, (4, 8, 16))
                self.assertTrue(assessed.constraint_admission.admission_passed)
                self.assertTrue(
                    assessed.spatial_spectral_admission.finest_pair_individual_budgets_passed
                )
                self.assertTrue(
                    assessed.spatial_spectral_admission.every_nested_tail_ratio_passed
                )
                self.assertTrue(assessed.admission_passed)

    def test_outer_defect_is_public_while_interior_defect_remains_a_veto(self) -> None:
        grids = tuple(UniformRadialGrid(0.0, 8.0, count) for count in (33, 65, 129))
        states = [_flat_state(grid, 4) for grid in grids]
        outer_states = []
        for state in states:
            q = state.q.copy()
            q[-5, 0] += 1.0e-4
            outer_states.append(EvolutionState(state.u, state.p, q))
        outer = proto8_gr0_common_event(
            outer_states,
            grids,
            accepted_stage_counts=(4, 8, 16),
            method="RK4",
            coordinate_time=0.25,
            measurement_radius_maximum=8.0,
        )
        full = outer.constraint_admission.full_domain_component_infinity[
            "reduction_alpha"
        ]
        owned = outer.constraint_admission.owned_domain_component_infinity[
            "reduction_alpha"
        ]
        self.assertGreater(max(full), 1.0e-6)
        self.assertLess(max(owned), max(outer.constraint_admission.accumulated_roundoff_enclosures))
        self.assertTrue(outer.constraint_admission.admission_passed)

        interior_states = []
        for state in states:
            q = state.q.copy()
            q[state.shape[0] // 3, 0] += 1.0e-6
            interior_states.append(EvolutionState(state.u, state.p, q))
        interior = proto8_gr0_common_event(
            interior_states,
            grids,
            accepted_stage_counts=(4, 8, 16),
            method="RK4",
            coordinate_time=0.25,
            measurement_radius_maximum=8.0,
        )
        self.assertFalse(interior.constraint_admission.admission_passed)
        self.assertFalse(interior.admission_passed)

    def test_frozen_ownership_and_input_contract_fail_closed(self) -> None:
        grids = tuple(UniformRadialGrid(0.0, 8.0, count) for count in (33, 65, 129))
        states = tuple(_flat_state(grid, 4) for grid in grids)
        with self.assertRaisesRegex(ValueError, "fixed outer-row"):
            proto8_gr0_common_event(
                states,
                grids,
                accepted_stage_counts=(0, 0, 0),
                method="RK4",
                coordinate_time=0.0,
                measurement_radius_maximum=8.0,
                fixed_outer_rows=3,
            )
        with self.assertRaisesRegex(ValueError, "accepted_stage_counts"):
            proto8_gr0_common_event(
                states,
                grids,
                accepted_stage_counts=(0, 0),
                method="RK4",
                coordinate_time=0.0,
                measurement_radius_maximum=8.0,
            )
        with self.assertRaisesRegex(ValueError, "increase strictly"):
            proto8_gr0_common_event(
                tuple(reversed(states)),
                tuple(reversed(grids)),
                accepted_stage_counts=(0, 0, 0),
                method="RK4",
                coordinate_time=0.0,
                measurement_radius_maximum=8.0,
            )


if __name__ == "__main__":
    unittest.main()
