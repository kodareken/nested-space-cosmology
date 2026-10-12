from __future__ import annotations

import json
import math
import sys
import unittest
from pathlib import Path

REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "src"))

from recursive_horizons.falsification import (  # noqa: E402
    DESI_DR2_FLAT_WCDM_SUMMARIES,
    BoundaryScalingSpec,
    PublishedGaussianSummary,
    approximate_gaussian_w_gate,
    arrival_delay_fraction,
    boundary_acceleration_gate,
    boundary_density_ratio,
    boundary_equation_of_state,
    disformal_photon_cone,
    gw170817_cone_gate,
    tensor_photon_delta,
)


class RelativeConeTests(unittest.TestCase):
    def test_delta_delay_and_gw170817_inclusive_gate(self) -> None:
        self.assertEqual(tensor_photon_delta(1.0, 1.0), 0.0)
        self.assertEqual(arrival_delay_fraction(0.0), 0.0)
        self.assertTrue(gw170817_cone_gate(-3.0e-15)["within_inclusive_bound"])
        self.assertTrue(gw170817_cone_gate(7.0e-16)["within_inclusive_bound"])
        self.assertFalse(gw170817_cone_gate(1.0e-14)["within_inclusive_bound"])
        self.assertFalse(gw170817_cone_gate(-1.0e-14)["within_inclusive_bound"])

    def test_disformal_cone_has_stable_null_limit(self) -> None:
        null = disformal_photon_cone(1.0, 4.0, 0.0)
        self.assertEqual(null["c_gamma_over_c_tensor"], 1.0)
        self.assertEqual(null["delta_tensor_over_photon"], 0.0)
        shifted = disformal_photon_cone(2.0, 0.5, 0.5)
        self.assertAlmostEqual(shifted["c_gamma_over_c_tensor"], math.sqrt(0.75))
        self.assertAlmostEqual(shifted["delta_tensor_over_photon"], 1.0 / math.sqrt(0.75) - 1.0)

    def test_invalid_relative_cone_inputs_are_rejected(self) -> None:
        for bad in (True, "1", complex(1.0, 0.0), math.inf, math.nan, 0.0, -1.0):
            with self.subTest(bad=repr(bad)):
                with self.assertRaises(ValueError):
                    tensor_photon_delta(bad, 1.0)
        with self.assertRaises(ValueError):
            gw170817_cone_gate(-1.0)
        with self.assertRaises(ValueError):
            disformal_photon_cone(1.0, 1.0, -0.1)
        with self.assertRaises(ValueError):
            disformal_photon_cone(1.0, 2.0, 0.5)


class BoundaryScalingTests(unittest.TestCase):
    def test_literal_scalings_and_acceleration_gate(self) -> None:
        vacuum = BoundaryScalingSpec(q=0.0, s=1.0, separately_conserved=True)
        surface = BoundaryScalingSpec(q=1.0, s=1.0, separately_conserved=True)
        inverse_area = BoundaryScalingSpec(q=2.0, s=1.0, separately_conserved=True)
        self.assertEqual(boundary_equation_of_state(vacuum), -1.0)
        self.assertAlmostEqual(boundary_equation_of_state(surface), -2.0 / 3.0)
        self.assertAlmostEqual(boundary_equation_of_state(inverse_area), -1.0 / 3.0)
        self.assertTrue(boundary_acceleration_gate(surface)["accelerates"])
        self.assertFalse(boundary_acceleration_gate(inverse_area)["accelerates"])
        self.assertAlmostEqual(boundary_density_ratio(surface, 2.0), 0.5)

    def test_nonconserved_scaling_has_no_w_or_acceleration_gate(self) -> None:
        spec = BoundaryScalingSpec(q=1.0, s=1.0, separately_conserved=False)
        self.assertAlmostEqual(boundary_density_ratio(spec, 2.0), 0.5)
        with self.assertRaises(ValueError):
            boundary_equation_of_state(spec)
        with self.assertRaises(ValueError):
            boundary_acceleration_gate(spec)

    def test_derived_overflow_and_density_underflow_are_rejected(self) -> None:
        product_overflow = BoundaryScalingSpec(
            q=1.0e308,
            s=1.0e308,
            separately_conserved=True,
        )
        with self.assertRaises(ValueError):
            boundary_equation_of_state(product_overflow)
        with self.assertRaises(ValueError):
            boundary_acceleration_gate(product_overflow)
        with self.assertRaises(ValueError):
            boundary_density_ratio(product_overflow, 2.0)

        density_underflow = BoundaryScalingSpec(
            q=1.0e308,
            s=1.0,
            separately_conserved=True,
        )
        with self.assertRaises(ValueError):
            boundary_density_ratio(density_underflow, 2.0)


class GaussianSummaryTests(unittest.TestCase):
    def test_constant_tension_wall_is_rejected_and_lambda_is_inside_all_summaries(self) -> None:
        for summary in DESI_DR2_FLAT_WCDM_SUMMARIES:
            with self.subTest(summary=summary.label):
                self.assertFalse(approximate_gaussian_w_gate(-2.0 / 3.0, summary)["within_approximate_interval"])
                self.assertTrue(approximate_gaussian_w_gate(-1.0, summary)["within_approximate_interval"])

    def test_validation_and_json_safe_records(self) -> None:
        with self.assertRaises(ValueError):
            BoundaryScalingSpec(q=True, s=1.0, separately_conserved=True)
        with self.assertRaises(ValueError):
            BoundaryScalingSpec(q=1.0, s=1.0, separately_conserved=1)  # type: ignore[arg-type]
        with self.assertRaises(ValueError):
            PublishedGaussianSummary("", -1.0, 0.1)
        with self.assertRaises(ValueError):
            PublishedGaussianSummary("bad", -1.0, 0.0)
        records = [
            gw170817_cone_gate(0.0),
            disformal_photon_cone(1.0, 0.0, 0.0),
            boundary_acceleration_gate(BoundaryScalingSpec(1.0, 1.0, True)),
            approximate_gaussian_w_gate(-1.0, DESI_DR2_FLAT_WCDM_SUMMARIES[0]),
        ]
        json.dumps(records, allow_nan=False)
        for record in records:
            self.assertTrue(all(not isinstance(value, float) or math.isfinite(value) for value in record.values()))


if __name__ == "__main__":
    unittest.main()
