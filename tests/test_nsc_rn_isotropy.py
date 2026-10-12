"""Bounded checks for the isotropy trajectory and its registered test."""

import unittest
from pathlib import Path

import numpy as np

from recursive_horizons import nsc_rn_isotropy as isotropy

LAB = Path(__file__).resolve().parents[1]
FAMILY = LAB / "results/development/nsc-rn-neutral-family-v1"


class DepartureTests(unittest.TestCase):
    def test_isotropic_profile_has_zero_departure(self):
        balance = np.full(16, 0.5)
        self.assertAlmostEqual(isotropy.isotropy_index(balance), 0.0, places=12)

    def test_known_departure_is_the_mean(self):
        balance = np.array([0.5, 0.7, 0.3, 1.0])
        expected = float(np.mean([0.0, 0.2, 0.2, 0.5]))
        self.assertAlmostEqual(isotropy.isotropy_index(balance), expected, places=12)

    def test_profile_is_symmetric_about_the_isotropic_point(self):
        below = isotropy.departure_profile(np.array([0.4, 0.45]))
        above = isotropy.departure_profile(np.array([0.6, 0.55]))
        self.assertTrue(np.allclose(below, above, atol=1e-12))

    def test_nonfinite_cells_are_excluded(self):
        balance = np.array([0.5, np.nan, 1.0])
        self.assertAlmostEqual(isotropy.isotropy_index(balance), 0.25, places=12)

    def test_empty_profile_is_rejected(self):
        with self.assertRaises(ValueError):
            isotropy.isotropy_index(np.array([]))


class RegisteredVerdictTests(unittest.TestCase):
    def _rows(self, indices, times=None, refinement=None):
        times = times or list(range(len(indices)))
        return [{"station": f"s{i}", "t_over_rm": float(t), "isotropy_index": float(value),
                 "refinement_indicator": refinement}
                for i, (t, value) in enumerate(zip(times, indices))]

    def test_decay_is_observed(self):
        rows = self._rows([0.10, 0.40, 0.05], refinement=0.001)
        verdict = isotropy.relaxation_verdict(rows)
        self.assertEqual(verdict["verdict"], "relaxation_observed")

    def test_growth_is_not_relaxation(self):
        rows = self._rows([0.10, 0.20, 0.45], refinement=0.001)
        verdict = isotropy.relaxation_verdict(rows)
        self.assertEqual(verdict["verdict"], "no_relaxation_observed")

    def test_small_fall_is_inconclusive_not_a_pass(self):
        # The tail does fall past the registered fraction, but the fall is
        # inside the refinement multiple: that is inconclusive, not a pass.
        rows = self._rows([0.40, 0.199], refinement=0.1)
        verdict = isotropy.relaxation_verdict(rows)
        self.assertEqual(verdict["verdict"], "inconclusive_on_this_domain")

    def test_fall_short_of_the_fraction_is_no_relaxation(self):
        rows = self._rows([0.40, 0.399], refinement=0.001)
        verdict = isotropy.relaxation_verdict(rows)
        self.assertEqual(verdict["verdict"], "no_relaxation_observed")

    def test_single_station_is_inconclusive(self):
        rows = self._rows([0.30])
        verdict = isotropy.relaxation_verdict(rows)
        self.assertEqual(verdict["verdict"], "no_relaxation_observed")

    def test_registration_is_written_into_the_verdict(self):
        verdict = isotropy.relaxation_verdict(self._rows([0.2, 0.1], refinement=0.0005))
        self.assertIn("before the relaxation run", verdict["registration"])
        self.assertEqual(verdict["refinement_multiple_registered"], isotropy.REFINEMENT_MULTIPLE)


class BalanceMirrorTests(unittest.TestCase):
    def test_balance_matches_the_pointwise_definition(self):
        radial = np.array([1.0, -2.0, 0.0, 3.0])
        angular = np.array([1.0, 2.0, 0.0, 1.0])
        expected = np.array([0.5, 0.5, 0.0, 0.75])
        self.assertTrue(np.allclose(isotropy.balance_from_pressures(radial, angular),
                                    expected, atol=1e-12))

    @unittest.skipUnless((FAMILY / "mu-1.010-eps-0.001/station-0.npz").is_file(),
                         "recorded RN family payloads are not present")
    def test_balance_mirrors_the_owning_module_exactly(self):
        from recursive_horizons import nsc_rn_gradient_bubble as bubble
        with np.load(FAMILY / "mu-1.010-eps-0.001/station-0.npz") as payload:
            owned = np.asarray(bubble.gradient_bubble_profile(payload)["balance"])
            mirrored = isotropy.balance_from_pressures(
                np.asarray(payload["radial_pressure"]), np.asarray(payload["angular_pressure"]))
        self.assertTrue(np.array_equal(owned, mirrored))


@unittest.skipUnless((FAMILY / "mu-1.010-eps-0.001/station-0.npz").is_file(),
                     "recorded RN family payloads are not present")
class RecordedStationTests(unittest.TestCase):
    def test_real_station_reduces(self):
        from recursive_horizons import nsc_rn_gradient_bubble as bubble
        with np.load(FAMILY / "mu-1.010-eps-0.001/station-0.npz") as payload:
            profile = bubble.gradient_bubble_profile(payload)
        index = isotropy.isotropy_index(np.asarray(profile["balance"]))
        self.assertTrue(np.isfinite(index))
        self.assertGreaterEqual(index, 0.0)
        self.assertLessEqual(index, 0.5)


if __name__ == "__main__":
    unittest.main()