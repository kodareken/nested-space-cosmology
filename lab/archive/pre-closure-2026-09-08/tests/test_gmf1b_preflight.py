from __future__ import annotations

import json
import math
import sys
import unittest
from pathlib import Path

REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "src"))

from recursive_horizons.gmf1b_preflight import (  # noqa: E402
    EC1BounceSpec,
    canonical_scalar_focusing_gate,
    canonical_scalar_null_energy,
    canonical_scalar_raychaudhuri_derivative,
    ec1_bounce_thermal_gate,
    thermal_ec_state,
)


class GMF1BPreflightTests(unittest.TestCase):
    def test_canonical_scalar_null_focusing_is_nonpositive(self) -> None:
        state = canonical_scalar_raychaudhuri_derivative(1.0, -1.0, 0.25, 0.5)
        self.assertEqual(state["T_kk"], 0.25)
        self.assertLess(state["raychaudhuri_derivative"], 0.0)
        gate = canonical_scalar_focusing_gate(1.0, -1.0, 0.25, 0.5)
        self.assertTrue(gate["same_generator_defocusing_to_zero_excluded_under_stated_assumptions"])
        self.assertEqual(canonical_scalar_null_energy(-3.0), 9.0)

    def test_derived_overflow_is_rejected(self) -> None:
        with self.assertRaises(ValueError):
            canonical_scalar_raychaudhuri_derivative(1.0e308, 1.0e308, 0.0, 0.0)
        with self.assertRaises(ValueError):
            thermal_ec_state(1.0e308)

    def test_thermal_regime_boundaries(self) -> None:
        self.assertEqual(thermal_ec_state(0.1)["classification"], "causal_barotrope")
        self.assertEqual(thermal_ec_state(2.0 / 9.0)["classification"], "degenerate_zero_sound_speed")
        self.assertEqual(thermal_ec_state(0.4)["classification"], "gradient_unstable")
        singular = thermal_ec_state(2.0 / 3.0)
        self.assertEqual(singular["classification"], "singular_enthalpy_and_density_derivative")
        self.assertIsNone(singular["cs2"])
        superluminal = thermal_ec_state(0.8)
        self.assertEqual(superluminal["classification"], "superluminal_relative_to_metric_cone")
        self.assertGreater(superluminal["cs2"], 1.0)
        immediately_below = thermal_ec_state(
            math.nextafter(2.0 / 3.0, -math.inf)
        )
        immediately_above = thermal_ec_state(
            math.nextafter(2.0 / 3.0, math.inf)
        )
        self.assertEqual(immediately_below["classification"], "gradient_unstable")
        self.assertLess(immediately_below["cs2"], 0.0)
        self.assertEqual(
            immediately_above["classification"],
            "superluminal_relative_to_metric_cone",
        )
        self.assertGreater(immediately_above["cs2"], 1.0)

    def test_ec1_bounce_rejects_only_naive_thermal_closure(self) -> None:
        gate = ec1_bounce_thermal_gate(EC1BounceSpec(10.0, 1.0))
        self.assertTrue(gate["z_bounce_strictly_between_one_half_and_one"])
        self.assertFalse(gate["healthy_causal_barotrope_at_bounce"])
        self.assertTrue(gate["naive_thermal_averaged_perfect_fluid_rejected_as_gmf1b_material_law"])
        self.assertFalse(gate["full_einstein_cartan_dirac_rejected"])
        self.assertFalse(gate["ec1_background_algebra_rejected"])

        near_degenerate = ec1_bounce_thermal_gate(EC1BounceSpec(2.0001, 1.0))
        self.assertEqual(
            near_degenerate["thermal_classification_at_bounce"],
            "gradient_unstable",
        )
        self.assertLess(near_degenerate["cs2_at_bounce"], 0.0)
        self.assertFalse(near_degenerate["healthy_causal_barotrope_at_bounce"])

    def test_high_dynamic_range_ec1_bounce_gate(self) -> None:
        gate = ec1_bounce_thermal_gate(EC1BounceSpec(1.0e16, 1.0))
        self.assertGreater(gate["x_min"], 0.0)
        self.assertGreater(gate["one_minus_z_bounce"], 0.0)
        self.assertTrue(gate["z_bounce_strictly_between_one_half_and_one"])
        self.assertFalse(gate["healthy_causal_barotrope_at_bounce"])
        with self.assertRaises(ValueError):
            ec1_bounce_thermal_gate(EC1BounceSpec(1.0e154, 5.0e-324))

    def test_invalid_inputs_and_json_safety(self) -> None:
        for bad in (True, "0.1", complex(1.0, 0.0), math.inf, math.nan, -1.0):
            with self.subTest(bad=repr(bad)):
                with self.assertRaises(ValueError):
                    thermal_ec_state(bad)  # type: ignore[arg-type]
        for values in ((True, 1.0), (2.0, 1.0), (1.0e308, 1.0)):
            with self.subTest(values=values):
                with self.assertRaises(ValueError):
                    EC1BounceSpec(*values)  # type: ignore[arg-type]
        with self.assertRaises(ValueError):
            canonical_scalar_raychaudhuri_derivative(1.0, -1.0, -0.1, 0.0)
        outputs = [
            canonical_scalar_focusing_gate(1.0, -1.0, 0.0, 0.0),
            thermal_ec_state(2.0 / 3.0),
            ec1_bounce_thermal_gate(EC1BounceSpec(10.0, 1.0)),
        ]
        json.dumps(outputs, allow_nan=False)

    def test_record_has_explicit_non_no_go_boundaries(self) -> None:
        from scripts.reproduce_gmf1b_preflight import record

        payload = record()
        self.assertFalse(payload["general_recursive_horizons_no_go_proven"])
        self.assertFalse(payload["canonical_scalar_regular_bridge_supported_under_stated_assumptions"])


if __name__ == "__main__":
    unittest.main()
