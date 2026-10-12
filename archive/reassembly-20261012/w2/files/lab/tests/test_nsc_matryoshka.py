"""Bounded checks for the Matryoshka comparison.

The two exact recurrences must reproduce the published prefixes, the
growth ordering must separate the factorial polygon count from the
Schroeder chain count, and the resolvent-level structure must carry one
inverse factor and one extra 1/Omega per level.
"""
import unittest

from recursive_horizons import nsc_matryoshka as matryoshka


class SequenceTests(unittest.TestCase):
    def test_matryoshka_prefix_matches_paper(self):
        computed = matryoshka.matryoshka_numbers(len(matryoshka.MATRYOSHKA_PREFIX) + 2)
        self.assertEqual(tuple(computed[: len(matryoshka.MATRYOSHKA_PREFIX)]),
                         matryoshka.MATRYOSHKA_PREFIX)

    def test_matryoshka_rejects_below_x2(self):
        with self.assertRaises(ValueError):
            matryoshka.matryoshka_numbers(1)

    def test_schroder_prefix_matches_reference(self):
        computed = matryoshka.schroder_numbers(len(matryoshka.SCHRODER_PREFIX))
        self.assertEqual(tuple(computed), matryoshka.SCHRODER_PREFIX)

    def test_schroder_rejects_empty_chain(self):
        with self.assertRaises(ValueError):
            matryoshka.schroder_numbers(0)


class GrowthTests(unittest.TestCase):
    def test_matryoshka_growth_diverges_beyond_chain(self):
        matryoshka_ratios = matryoshka.growth_ratio(matryoshka.matryoshka_numbers(7))
        schroder_ratios = matryoshka.growth_ratio(matryoshka.schroder_numbers(6))
        self.assertTrue(all(r > s for r, s in zip(matryoshka_ratios, schroder_ratios[:len(matryoshka_ratios)])))
        # The polygon count is factorial: its ratio grows with n, while the
        # chain count is exponential with an approximately constant ratio.
        self.assertTrue(all(matryoshka_ratios[i] > matryoshka_ratios[i - 1]
                            for i in range(1, len(matryoshka_ratios))))


class ResolventStructureTests(unittest.TestCase):
    def test_each_level_has_one_inverse_and_one_rescale(self):
        levels = matryoshka.resolvent_levels(4.0, 3)
        self.assertEqual(len(levels), 3)
        for index, level in enumerate(levels, start=1):
            self.assertEqual(level["level"], index)
            self.assertEqual(level["own_inverse_factors"], 1)
            self.assertAlmostEqual(level["argument_rescale"], 1.0 / (4.0 ** index), places=12)
            self.assertAlmostEqual(level["coefficient_scale_per_link"], 0.25, places=12)

    def test_resolvent_levels_reject_bad_scale(self):
        with self.assertRaises(ValueError):
            matryoshka.resolvent_levels(0.0, 2)
        with self.assertRaises(ValueError):
            matryoshka.resolvent_levels(-2.0, 2)

    def test_comparison_names_the_divergence(self):
        body = matryoshka.comparison()
        self.assertTrue(body["matryoshka_prefix_matches_paper"])
        self.assertTrue(body["schroder_prefix_matches_reference"])
        self.assertTrue(body["weight_comparison"]["skeleton_matches"])
        self.assertFalse(body["weight_comparison"]["weight_type_matches"])
        self.assertTrue(body["weight_comparison"]["no_cosmohedron_equivalence_identified"])


if __name__ == "__main__":
    unittest.main()