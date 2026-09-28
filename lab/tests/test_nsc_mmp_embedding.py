import unittest

from recursive_horizons.nsc_mmp_embedding import CommonCoefficients


class MMPEmbeddingTests(unittest.TestCase):
    def test_common_action_dictionary_preserves_charge_radius(self):
        coefficients = CommonCoefficients(einstein=2.5, gauge=0.2, vacuum=0.0)
        one = coefficients.mmp_map(flux=1)
        three = coefficients.mmp_map(flux=3)
        self.assertAlmostEqual(three["r_e_squared"], 9 * one["r_e_squared"])
        self.assertAlmostEqual(three["ell_MMP"], 9 * one["ell_MMP"])
        self.assertEqual(one["lambda_4"], 0.0)
        self.assertEqual(one["Xi"], 0.0)

    def test_flux_is_a_nonzero_integer_sector(self):
        coefficients = CommonCoefficients(einstein=1.0, gauge=1.0, vacuum=0.0)
        for value in (0, 1.5, True):
            with self.assertRaises(ValueError):
                coefficients.mmp_map(flux=value)


if __name__ == "__main__":
    unittest.main()
