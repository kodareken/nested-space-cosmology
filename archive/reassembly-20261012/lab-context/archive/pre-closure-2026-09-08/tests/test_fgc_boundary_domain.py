from __future__ import annotations

from fractions import Fraction
from pathlib import Path
import sys
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from recursive_horizons.fgc.evolution.boundary_domain import (  # noqa: E402
    ALL_CONE_LOCAL_SPEED_BOUND,
    BoundaryControlStop,
    BoundaryGeometry,
    CausalBudgetState,
    accept_causal_step,
    aligned_outer_boundary_grids,
    coordinate_speed_upper_bound,
    exact_coordinate_speed_upper_bound,
    exact_reference_budget,
)


Q = Fraction


class FGCBoundaryDomainTests(unittest.TestCase):
    def test_ADM_frame_conversion_includes_shift_lapse_metric_and_all_cones(self) -> None:
        self.assertEqual(ALL_CONE_LOCAL_SPEED_BOUND, Q(6, 5))
        self.assertEqual(
            exact_coordinate_speed_upper_bound(
                lapse=Q(3, 2), radial_metric=Q(5, 4), shift=Q(-1, 10)
            ),
            Q(77, 50),
        )
        self.assertAlmostEqual(
            coordinate_speed_upper_bound(lapse=1.5, radial_metric=1.25, shift=-0.1),
            1.54,
        )
        with self.assertRaisesRegex(ValueError, "positive"):
            coordinate_speed_upper_bound(lapse=1.0, radial_metric=0.0, shift=0.0)

    def test_moved_boundaries_preserve_spacing_and_exact_reference_margin(self) -> None:
        grids = aligned_outer_boundary_grids(
            nominal_outer_radius=Q(128),
            nominal_point_counts=(1025, 2049, 4097),
            outer_radii=(Q(96), Q(128), Q(160)),
            measurement_radius=Q(24),
            stencil_reach_intervals=3,
        )
        self.assertEqual(len(grids), 9)
        self.assertEqual([grid.point_count for grid in grids[:3]], [769, 1025, 1281])
        smallest = exact_reference_budget(
            grids[0], final_time=Q(32), minimum_causal_buffer=Q(16)
        )
        self.assertEqual(smallest["remaining_causal_buffer"], Q(1329, 40))
        self.assertEqual(smallest["strict_margin_over_required_buffer"], Q(689, 40))
        self.assertTrue(all(exact_reference_budget(grid, final_time=Q(32), minimum_causal_buffer=Q(16))["passed"] for grid in grids))

    def test_trial_crossing_is_rejected_without_mutating_last_accepted_state(self) -> None:
        geometry = BoundaryGeometry(96, 24, 16, 3 / 8)
        initial = CausalBudgetState(previous_speed_upper=1.2)
        accepted = accept_causal_step(
            initial,
            geometry,
            trial_time=1,
            stage_coordinate_speed_uppers=(1.1, 1.2, 1.15, 1.2),
            candidate_endpoint_coordinate_speed_upper=1.3,
        )
        with self.assertRaises(BoundaryControlStop) as caught:
            accept_causal_step(
                accepted,
                geometry,
                trial_time=32,
                stage_coordinate_speed_uppers=(4, 4, 4, 4),
                candidate_endpoint_coordinate_speed_upper=5,
            )
        self.assertEqual(caught.exception.reason, "boundary_causal_buffer")
        self.assertLess(caught.exception.assessment.strict_margin_over_required_buffer, 0)
        self.assertEqual(accepted.accepted_time, 1.0)
        self.assertEqual(accepted.accumulated_characteristic_distance, 1.3)
        self.assertEqual(accepted.previous_speed_upper, 1.3)

    def test_candidate_endpoint_speed_is_not_substituted_by_last_RK_stage(self) -> None:
        geometry = BoundaryGeometry(96, 24, 16, 3 / 8)
        accepted = accept_causal_step(
            CausalBudgetState(previous_speed_upper=1.0),
            geometry,
            trial_time=1,
            stage_coordinate_speed_uppers=(1.0, 1.1, 1.2, 1.1),
            candidate_endpoint_coordinate_speed_upper=2.0,
        )
        self.assertEqual(accepted.accumulated_characteristic_distance, 2.0)
        self.assertEqual(accepted.previous_speed_upper, 2.0)


if __name__ == "__main__":
    unittest.main()
