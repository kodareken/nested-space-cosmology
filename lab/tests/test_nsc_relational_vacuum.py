import unittest

from recursive_horizons.nsc_mmp_embedding import CommonCoefficients
from recursive_horizons.nsc_relational_vacuum import (
    RelationalVacuumNormalization,
    uniqueness_residual,
)


class RelationalVacuumTests(unittest.TestCase):
    def setUp(self):
        self.raw = CommonCoefficients(0.0055, 0.0048, 0.1725)
        self.law = RelationalVacuumNormalization()

    def test_normalization_is_idempotent_and_preserves_gradient_coefficients(self):
        first = self.law.normalize(self.raw).normalized
        second = self.law.normalize(first).normalized
        self.assertEqual(first, second)
        self.assertEqual(first.vacuum, 0.0)
        self.assertEqual(first.einstein, self.raw.einstein)
        self.assertEqual(first.gauge, self.raw.gauge)

    def test_direct_sum_scale_and_ctp_normalization(self):
        left = self.law.counterterm(self.raw.vacuum, 2.0)
        right = self.law.counterterm(self.raw.vacuum, 3.0)
        combined = self.law.counterterm(self.raw.vacuum, 5.0)
        self.assertAlmostEqual(left + right, combined)
        self.assertEqual(self.law.scale_covariance_defect(0.3, 7.0, 4.5), 0.0)
        self.assertEqual(self.law.ctp_counterterm(0.3, 7.0, 7.0), 0.0)

    def test_normalized_coefficients_pass_the_mmp_vacuum_gate(self):
        result = self.law.mmp_gate(self.raw)
        self.assertEqual(result["vacuum_seed_Xi"], 0.0)
        self.assertTrue(result["charged_seed_bound_passes"])
        self.assertTrue(result["asymptotically_flat_bulk_vacuum_passes"])
        self.assertTrue(result["Einstein_coefficient_unchanged"])
        self.assertTrue(result["gauge_coefficient_unchanged"])

    def test_projector_requirements_vanish(self):
        residuals = uniqueness_residual()
        for name, value in residuals.items():
            if name.endswith("preserved") or name.endswith("annihilated") or name.endswith("idempotence"):
                self.assertEqual(value, 0.0, name)


if __name__ == "__main__":
    unittest.main()
