import unittest

from recursive_horizons.nsc_scale_binding import (
    ScaleBindingConventions,
    ScaleBindingSolver,
)


class ScaleBindingTests(unittest.TestCase):
    def test_fixed_q_radius_root_and_physical_cutoff_order(self):
        solver = ScaleBindingSolver(ScaleBindingConventions(
            magnetic_flux=4,
            compact_modes=18,
            quadrature_points=72,
            quadrature_tolerance=8e-8,
        ))
        result = solver.solve((3.0, 5.0), recursive_transfer=0.0446327233759311)
        self.assertLess(abs(result.radius_residual), 2e-8)
        self.assertGreater(result.cutoff_over_magnetic_scale, 1)
        self.assertLess(result.child_H_over_cutoff, 1)
        self.assertGreater(result.conditional_gap_squared_over_cutoff, 0)
        self.assertLess(result.recursive_endpoint_product, 1)

    def test_cutoff_must_exceed_physical_matching_scale(self):
        solver = ScaleBindingSolver(ScaleBindingConventions(magnetic_flux=4))
        with self.assertRaises(ValueError):
            solver.coefficients(2.0)


if __name__ == "__main__":
    unittest.main()

