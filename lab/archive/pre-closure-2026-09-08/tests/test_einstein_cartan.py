from __future__ import annotations

import json
import math
import sys
import unittest
from pathlib import Path

REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "src"))

from recursive_horizons.einstein_cartan import (  # noqa: E402
    EinsteinCartanSpec,
    boundary_misner_sharp_mass,
    boundary_work_residual,
    causal_sequence,
    continuity_residual,
    curvature_invariants,
    fluid_state,
    friedmann_F,
    half_cycle_convergence,
    marginal_scale_factor_roots,
    scale_factor_acceleration,
    static_parent_mass_drift,
    turning_points,
)


class EinsteinCartanReductionTests(unittest.TestCase):
    def setUp(self) -> None:
        self.spec = EinsteinCartanSpec(A=10.0, B=1.0)

    def test_exact_turning_points_and_acceleration(self) -> None:
        roots = turning_points(self.spec)
        self.assertAlmostEqual(roots["x_min"] + roots["x_max"], 10.0)
        self.assertAlmostEqual(roots["x_min"] * roots["x_max"], 1.0)
        self.assertAlmostEqual(friedmann_F(self.spec, roots["a_min"]), 0.0)
        self.assertAlmostEqual(friedmann_F(self.spec, roots["a_max"]), 0.0)
        self.assertGreater(scale_factor_acceleration(self.spec, roots["a_min"]), 0.0)
        self.assertLess(scale_factor_acceleration(self.spec, roots["a_max"]), 0.0)

    def test_fluid_continuity_and_finite_bounce_invariants(self) -> None:
        roots = turning_points(self.spec)
        a_min = roots["a_min"]
        state = fluid_state(self.spec, a_min)
        self.assertEqual(state["w_radiation"], 1.0 / 3.0)
        self.assertEqual(state["w_torsion"], 1.0)
        self.assertEqual(state["q_internal"], 0.0)
        self.assertGreater(state["rho_radiation"], 0.0)
        self.assertLess(state["rho_torsion"], 0.0)
        self.assertLess(state["rho_plus_p"], 0.0)
        self.assertLess(state["rho_plus_3p"], 0.0)
        self.assertAlmostEqual(continuity_residual(self.spec, a_min, "expansion"), 0.0)
        invariants = curvature_invariants(self.spec, a_min)
        self.assertEqual(invariants["weyl_tensor_squared"], 0.0)
        self.assertTrue(all(math.isfinite(value) for value in invariants.values()))
        self.assertGreater(invariants["ricci_scalar"], 0.0)

    def test_marginal_roots_and_one_metric_causal_sequence(self) -> None:
        roots = turning_points(self.spec)
        low, high = marginal_scale_factor_roots(self.spec, self.spec.chi_boundary)
        self.assertLess(roots["a_min"], low)
        self.assertLess(low, high)
        self.assertLess(high, roots["a_max"])
        sequence = causal_sequence(self.spec)
        self.assertTrue(sequence["single_interior_metric_causal_sequence"])
        self.assertEqual(
            [item["region"] for item in sequence["sequence"]],
            [
                "normal",
                "marginal",
                "trapped",
                "marginal",
                "normal",
                "marginal",
                "anti_trapped",
                "marginal",
                "normal",
            ],
        )

    def test_high_dynamic_range_roots_and_factorized_friedmann_value(self) -> None:
        high_range = EinsteinCartanSpec(A=1.0e16, B=1.0)
        roots = turning_points(high_range)
        self.assertGreater(roots["x_min"], 0.0)
        self.assertAlmostEqual(roots["x_min"] * roots["x_max"], 1.0, places=12)
        self.assertEqual(friedmann_F(high_range, roots["a_min"]), 0.0)
        self.assertEqual(friedmann_F(high_range, roots["a_max"]), 0.0)
        marginal_low, marginal_high = marginal_scale_factor_roots(
            high_range, high_range.chi_boundary
        )
        self.assertGreater(marginal_low, 0.0)
        self.assertGreater(marginal_high, marginal_low)
        self.assertLess(marginal_high, roots["a_max"])

    def test_boundary_work_identity_and_static_parent_mass_drift(self) -> None:
        roots = turning_points(self.spec)
        probe = (roots["a_min"] + roots["a_max"]) / 2.0
        self.assertAlmostEqual(boundary_work_residual(self.spec, probe, "collapse"), 0.0)
        self.assertAlmostEqual(boundary_work_residual(self.spec, probe, "expansion"), 0.0)
        drift = static_parent_mass_drift(self.spec)
        self.assertNotEqual(drift["boundary_mass_drift_to_bounce"], 0.0)
        self.assertFalse(drift["fixed_vacuum_parent_can_remain_smooth_without_flux_or_layer"])
        self.assertAlmostEqual(drift["bounce_to_initial_parent_mass_ratio"], 0.1010205144, places=8)
        self.assertNotAlmostEqual(
            boundary_misner_sharp_mass(self.spec, roots["a_min"]),
            boundary_misner_sharp_mass(self.spec, roots["a_max"]),
        )

    def test_smooth_simpson_refinement(self) -> None:
        convergence = half_cycle_convergence(self.spec, 32)
        self.assertGreater(convergence["coarse"], 0.0)
        self.assertLess(convergence["fine_to_finer_difference"], convergence["coarse_to_fine_difference"])

    def test_invalid_inputs_and_json_safe_outputs(self) -> None:
        for bad in (True, "10", math.inf, math.nan, 0.0, -1.0):
            with self.subTest(bad=repr(bad)):
                with self.assertRaises(ValueError):
                    EinsteinCartanSpec(A=bad, B=1.0)  # type: ignore[arg-type]
        with self.assertRaises(ValueError):
            EinsteinCartanSpec(A=2.0, B=1.0)
        with self.assertRaises(ValueError):
            EinsteinCartanSpec(A=1.0e308, B=1.0)
        with self.assertRaises(ValueError):
            EinsteinCartanSpec(A=10.0, B=1.0, chi_boundary=math.pi / 2.0)
        with self.assertRaises(ValueError):
            friedmann_F(self.spec, 0.1)
        with self.assertRaises(ValueError):
            continuity_residual(self.spec, 1.0, "neither")  # type: ignore[arg-type]
        with self.assertRaises(ValueError):
            marginal_scale_factor_roots(self.spec, 0.0)
        record = {
            "turning": turning_points(self.spec),
            "fluid": fluid_state(self.spec, turning_points(self.spec)["a_min"]),
            "curvature": curvature_invariants(self.spec, turning_points(self.spec)["a_min"]),
            "sequence": causal_sequence(self.spec),
            "time": half_cycle_convergence(self.spec),
        }
        json.dumps(record, allow_nan=False)


if __name__ == "__main__":
    unittest.main()
