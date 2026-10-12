from __future__ import annotations

import json
import math
import sys
import unittest
from pathlib import Path

REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "src"))

from recursive_horizons.nested_spectrum import (  # noqa: E402
    eventual_occurrence_probability,
    finite_band_selection,
    required_nonclaims,
    two_level_recurrence_ratio,
)


class NestedSpectrumControlTests(unittest.TestCase):
    def test_finite_band_selector_has_declared_nonzero_scale(self) -> None:
        control = finite_band_selection(alpha=2.0, beta=1.0)
        self.assertEqual(control["k_star_squared"], 1.0)
        self.assertEqual(control["k_star"], 1.0)
        self.assertEqual(control["ell_star"], 1.0)
        self.assertEqual(control["energy_coefficient_at_k_star"], -1.0)
        self.assertAlmostEqual(control["stationary_derivative_residual"], 0.0)
        self.assertTrue(control["nonzero_finite_scale_selected"])
        self.assertIn("not_covariant_gravity", control["classification"])

    def test_golden_ratio_is_only_the_one_one_recurrence_fixture(self) -> None:
        golden = two_level_recurrence_ratio(a=1.0, b=1.0)
        expected = 0.5 * (1.0 + math.sqrt(5.0))
        self.assertAlmostEqual(golden["positive_ratio"], expected)
        self.assertAlmostEqual(golden["scaled_characteristic_residual"], 0.0)
        self.assertTrue(golden["golden_ratio_fixture"])
        self.assertFalse(golden["coefficients_derived_from_physics"])

        other = two_level_recurrence_ratio(a=2.0, b=1.0)
        self.assertAlmostEqual(other["positive_ratio"], 1.0 + math.sqrt(2.0))
        self.assertFalse(other["golden_ratio_fixture"])

        large_b = two_level_recurrence_ratio(a=0.0, b=1.0e308)
        self.assertAlmostEqual(large_b["positive_ratio"], 1.0e154)
        self.assertAlmostEqual(large_b["scaled_characteristic_residual"], 0.0)

        cancellation_control = two_level_recurrence_ratio(a=1.0e100, b=1.0)
        self.assertEqual(cancellation_control["positive_ratio"], 1.0e100)
        self.assertLess(
            abs(cancellation_control["scaled_characteristic_residual"]),
            1.0e-15,
        )

    def test_eternal_opportunity_control_is_conditional_and_stable(self) -> None:
        control = eventual_occurrence_probability(1.0e-6, 1_000_000)
        expected = -math.expm1(1_000_000 * math.log1p(-1.0e-6))
        self.assertAlmostEqual(control["probability_at_least_one"], expected)
        self.assertAlmostEqual(control["normalization_residual"], 0.0)
        self.assertEqual(control["infinite_independent_trial_limit"], 1.0)
        self.assertTrue(control["independence_assumed"])
        self.assertFalse(control["ergodicity_derived"])
        self.assertFalse(control["logical_necessity_proven"])

        self.assertEqual(eventual_occurrence_probability(0.0, 10)["probability_at_least_one"], 0.0)
        self.assertEqual(eventual_occurrence_probability(1.0, 1)["probability_at_least_one"], 1.0)
        self.assertEqual(eventual_occurrence_probability(0.5, 0)["probability_at_least_one"], 0.0)

    def test_invalid_inputs_are_rejected(self) -> None:
        for values in ((0.0, 1.0), (1.0, 0.0), (math.inf, 1.0), (True, 1.0)):
            with self.subTest(selector=values):
                with self.assertRaises(ValueError):
                    finite_band_selection(*values)  # type: ignore[arg-type]
        for values in ((-1.0, 1.0), (1.0, 0.0), (math.nan, 1.0)):
            with self.subTest(recurrence=values):
                with self.assertRaises(ValueError):
                    two_level_recurrence_ratio(*values)
        for values in ((-0.1, 1), (1.1, 1), (0.5, -1), (0.5, True), (0.5, 10**15 + 1)):
            with self.subTest(probability=values):
                with self.assertRaises(ValueError):
                    eventual_occurrence_probability(*values)  # type: ignore[arg-type]

    def test_record_is_json_safe_and_retains_every_nonclaim(self) -> None:
        from scripts.reproduce_nested_spectrum import record

        payload = record()
        self.assertEqual(payload["schema_version"], 1)
        self.assertEqual(payload["project_version"], "0.11.0")
        self.assertEqual(payload["model_id"], "NGS-0")
        self.assertEqual(payload["nonclaims"], required_nonclaims())
        self.assertTrue(all(value is False for value in payload["nonclaims"].values()))
        json.dumps(payload, allow_nan=False)


if __name__ == "__main__":
    unittest.main()
