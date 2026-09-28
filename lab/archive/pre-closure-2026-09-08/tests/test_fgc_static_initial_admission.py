from __future__ import annotations

from pathlib import Path
import sys
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "src"))

from recursive_horizons.fgc.evolution.static_initial_admission import (  # noqa: E402
    construct_fgcqr_grid_initial_data,
)
from recursive_horizons.fgc.initial_data_family import (  # noqa: E402
    InitialDataParameters,
)


PHI_SEED = 1.0 / 131072.0


class FGCStaticInitialAdmissionTests(unittest.TestCase):
    def test_complete_grid_preserves_constraints_buffers_and_center(self) -> None:
        initial = construct_fgcqr_grid_initial_data(
            InitialDataParameters(2.5, 2.0, PHI_SEED),
            point_count=1025,
            constraint_method="RK4",
            diagnostic_spatial_order=4,
        )
        self.assertEqual(initial.state.shape, (1025, 6))
        self.assertEqual(initial.support_minimum_index, 80)
        self.assertEqual(initial.support_maximum_index, 112)
        self.assertTrue(initial.regular_center)
        self.assertTrue(initial.exact_inner_vacuum_buffer)
        self.assertTrue(initial.exact_outer_vacuum_buffer)
        self.assertTrue(initial.finite_mass)
        self.assertTrue(initial.no_initial_trapped_sphere)
        self.assertTrue(initial.compactness_inside_protocol_window)
        self.assertLess(initial.physical_constraint_residual_infinity, 1.0e-14)
        self.assertGreater(initial.minimum_effective_planck_coefficient, 3.9)
        self.assertGreater(initial.minimum_vacuum_metric_denominator, 0.8)
        self.assertEqual(initial.state.u[0, 3], 0.0)
        self.assertEqual(initial.state.q[0, 3], initial.state.u[0, 2])

    def test_both_constraint_integrators_construct_the_widest_case(self) -> None:
        parameters = InitialDataParameters(3.0, 2.25, PHI_SEED)
        primary = construct_fgcqr_grid_initial_data(
            parameters,
            point_count=2049,
            constraint_method="RK4",
            diagnostic_spatial_order=4,
        )
        comparator = construct_fgcqr_grid_initial_data(
            parameters,
            point_count=2049,
            constraint_method="SSPRK3",
            diagnostic_spatial_order=2,
        )
        self.assertEqual(primary.constraint_solution.step_count, 72)
        self.assertEqual(comparator.constraint_solution.step_count, 72)
        self.assertLess(abs(primary.outer_mass - comparator.outer_mass), 1.0e-3)
        self.assertLess(
            abs(primary.peak_compactness - comparator.peak_compactness), 1.0e-3
        )
        self.assertTrue(primary.no_initial_trapped_sphere)
        self.assertTrue(comparator.no_initial_trapped_sphere)

    def test_compactness_window_failures_are_preserved_not_clipped(self) -> None:
        below = construct_fgcqr_grid_initial_data(
            InitialDataParameters(2.0, 2.25, PHI_SEED),
            point_count=4097,
        )
        above = construct_fgcqr_grid_initial_data(
            InitialDataParameters(5.0, 1.75, PHI_SEED),
            point_count=4097,
        )
        self.assertLess(below.peak_compactness, 0.1)
        self.assertFalse(below.compactness_inside_protocol_window)
        self.assertGreater(above.peak_compactness, 0.75)
        self.assertFalse(above.compactness_inside_protocol_window)
        self.assertTrue(below.no_initial_trapped_sphere)
        self.assertTrue(above.no_initial_trapped_sphere)

    def test_grid_alignment_and_constraint_limit_fail_closed(self) -> None:
        parameters = InitialDataParameters(2.5, 2.25, PHI_SEED)
        with self.assertRaisesRegex(ValueError, "aligned"):
            construct_fgcqr_grid_initial_data(parameters, point_count=1026)
        with self.assertRaisesRegex(ValueError, "positive"):
            construct_fgcqr_grid_initial_data(
                parameters,
                point_count=1025,
                maximum_constraint_residual=0.0,
            )


if __name__ == "__main__":
    unittest.main()
