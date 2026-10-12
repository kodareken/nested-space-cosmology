from __future__ import annotations

import json
import math
import sys
import unittest
from pathlib import Path

REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "src"))

from recursive_horizons import (  # noqa: E402
    C,
    CosmologyInputs,
    compactness,
    critical_density,
    flat_lcdm_age,
    hubble_mass,
    hubble_parameter_si,
    hubble_radius,
    lambda_from_h0_omega_lambda,
    parent_inference_from_lambda,
    reproduce_core,
)


class StandardCosmologyTests(unittest.TestCase):
    def setUp(self) -> None:
        self.inputs = CosmologyInputs()
        self.hubble_si = hubble_parameter_si(self.inputs.h0_km_s_mpc)

    def test_hubble_conversion(self) -> None:
        self.assertAlmostEqual(self.hubble_si, 2.1842852410855023e-18, places=30)

    def test_critical_density_scale(self) -> None:
        expected = 8.532855163739318e-27
        self.assertAlmostEqual(
            critical_density(self.hubble_si), expected, delta=expected * 1e-15
        )

    def test_hubble_sphere_compactness_is_identity(self) -> None:
        ratio = compactness(
            hubble_mass(self.hubble_si), hubble_radius(self.hubble_si)
        )
        self.assertAlmostEqual(ratio, 1.0, places=14)

    def test_lambda_conversion(self) -> None:
        value = lambda_from_h0_omega_lambda(67.4, 0.685)
        expected = 1.090910502838093e-52
        self.assertAlmostEqual(value, expected, delta=expected * 1e-15)

    def test_age_benchmark(self) -> None:
        age_gyr = flat_lcdm_age(self.inputs) / (1e9 * 31_557_600.0)
        self.assertAlmostEqual(age_gyr, 13.796234644007173, places=11)


class ConditionalInferenceTests(unittest.TestCase):
    def test_hsat_is_labeled_as_inference(self) -> None:
        lambda_m2 = lambda_from_h0_omega_lambda(67.4, 0.685)
        result = parent_inference_from_lambda(lambda_m2)
        self.assertEqual(
            result["classification"],
            "conditional_inference_not_independent_prediction",
        )
        self.assertLessEqual(float(result["relative_closure_error"]), 5e-16)

    def test_record_is_finite_and_json_serializable(self) -> None:
        record = reproduce_core()
        encoded = json.dumps(record, allow_nan=False)
        self.assertIn("algebraic_resubstitution_not_prediction", encoded)


class ValidationTests(unittest.TestCase):
    def test_nonpositive_values_are_rejected(self) -> None:
        for bad in (0.0, -1.0, math.inf, math.nan):
            with self.subTest(bad=bad):
                with self.assertRaises(ValueError):
                    hubble_parameter_si(bad)

    def test_age_requires_flat_matter_plus_lambda_inputs(self) -> None:
        with self.assertRaises(ValueError):
            flat_lcdm_age(CosmologyInputs(67.4, 0.3, 0.6))

    def test_compactness_has_no_hidden_c_override(self) -> None:
        radius = 2.0
        mass_for_unity = C**2 * radius / (2.0 * 6.67430e-11)
        self.assertAlmostEqual(compactness(mass_for_unity, radius), 1.0)


if __name__ == "__main__":
    unittest.main()
