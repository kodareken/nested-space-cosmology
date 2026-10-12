"""Bounded checks for the effective reach of the two RN radial forces.

The synthetic cases pin the definitions on inputs whose answer is known
exactly, in particular that reach is taken at the LAST crossing and so
is immune to a near-throat dip.
"""
import json
import unittest
from pathlib import Path

import numpy as np

from recursive_horizons import nsc_rn_force_reach as reach

LAB = Path(__file__).resolve().parents[1]
RECORD = LAB / "results/development/nsc-rn-force-reach-v1.json"


def _payload(inner, outer, radius=None):
    inner = np.asarray(inner, dtype=float)
    if radius is None:
        radius = np.linspace(1.04, 32.0, inner.size)
    return {"radius": np.asarray(radius, dtype=float),
            "radial_pressure": inner,
            "angular_pressure": np.asarray(outer, dtype=float)}


class LastCrossingTests(unittest.TestCase):
    def test_reach_ignores_a_near_throat_dip(self):
        """The definition must not be fooled by the ringing."""
        radius = np.linspace(1.0, 10.0, 1000)
        inner = np.ones_like(radius)
        inner[1:20] = 1.0e-12          # a deep dip just inside the edge
        outer = 1.0e-9 * np.ones_like(radius)
        profile = reach.force_reach(_payload(inner, outer, radius))
        self.assertAlmostEqual(profile["reach_inner"], 10.0, delta=0.02)

    def test_reach_stops_where_the_force_actually_dies(self):
        radius = np.linspace(1.0, 10.0, 1000)
        inner = np.ones_like(radius)
        inner[500:] = 1.0e-12
        outer = 1.0e-9 * np.ones_like(radius)
        profile = reach.force_reach(_payload(inner, outer, radius))
        self.assertAlmostEqual(profile["reach_inner"], 5.5, delta=0.05)

    def test_supported_prefix_ignores_a_later_return(self):
        radius = np.linspace(1.0, 10.0, 1000)
        inner = np.ones_like(radius)
        inner[500:520] = 1.0e-12
        inner[520:] = 1.0
        outer = 1.0e-9 * np.ones_like(radius)
        profile = reach.force_reach(_payload(inner, outer, radius))
        self.assertLess(profile["supported_prefix_frac_inner"], 0.6)
        self.assertAlmostEqual(profile["reach_inner"], 10.0, delta=0.02)


class VerdictTests(unittest.TestCase):
    def test_larger_force_reaching_further(self):
        radius = np.linspace(1.0, 10.0, 500)
        inner = np.exp(-radius / 3.0)
        outer = np.exp(-radius / 1.0)
        profile = reach.force_reach(_payload(inner, outer, radius))
        self.assertEqual(reach.force_reach_verdict(profile), "larger_reaches_further")

    def test_larger_force_reaching_less_far(self):
        """A steep dominant force can out-weigh but under-reach."""
        radius = np.linspace(1.0, 10.0, 500)
        inner = 100.0 * np.exp(-radius / 0.5)
        outer = np.exp(-radius / 4.0)
        profile = reach.force_reach(_payload(inner, outer, radius))
        self.assertGreater(profile["peak_inner"], profile["peak_outer"])
        self.assertEqual(profile["reach_longer"], "angular")
        self.assertEqual(reach.force_reach_verdict(profile), "larger_reaches_less_far")

    def test_equal_reach(self):
        radius = np.linspace(1.0, 10.0, 500)
        inner = np.exp(-radius / 2.0)
        outer = np.exp(-radius / 2.0)
        profile = reach.force_reach(_payload(inner, outer, radius))
        self.assertEqual(profile["reach_longer"], "equal")
        self.assertEqual(reach.force_reach_verdict(profile), "equal_reach")

    def test_dominant_fraction_is_one_when_one_force_always_wins(self):
        radius = np.linspace(1.0, 10.0, 200)
        profile = reach.force_reach(_payload(np.ones_like(radius),
                                             0.5 * np.ones_like(radius), radius))
        self.assertAlmostEqual(profile["dominant_fraction_inner"], 1.0)


class CrossoverTests(unittest.TestCase):
    def test_single_swap_is_found(self):
        radius = np.linspace(1.0, 10.0, 500)
        inner = np.where(radius < 5.0, 2.0, 0.5)
        outer = np.ones_like(radius)
        profile = reach.force_reach(_payload(inner, outer, radius))
        self.assertEqual(profile["crossover_count"], 1)
        self.assertAlmostEqual(profile["crossover_radii"][0], 5.0, delta=0.05)


class InputGuardTests(unittest.TestCase):
    def test_non_increasing_grid_refused(self):
        with self.assertRaises(ValueError):
            reach.force_reach(_payload(np.ones(3), np.ones(3), np.array([1.0, 0.9, 1.1])))

    def test_shape_mismatch_refused(self):
        with self.assertRaises(ValueError):
            reach.force_reach({"radius": np.ones(5),
                               "radial_pressure": np.ones(4),
                               "angular_pressure": np.ones(5)})

    def test_short_grid_refused(self):
        with self.assertRaises(ValueError):
            reach.force_reach(_payload(np.ones(2), np.ones(2)))


@unittest.skipUnless(RECORD.is_file(), "force-reach record not written yet")
class RecordTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record = json.loads(RECORD.read_text(encoding="utf-8"))

    def test_record_is_creation_only(self):
        self.assertTrue(self.record["creation_only"])

    def test_record_covers_the_whole_family(self):
        self.assertEqual(self.record["coverage"]["cases"], 6)
        self.assertEqual(self.record["coverage"]["total_stations"], 177)

    def test_scope_note_disclaims_a_force_law(self):
        self.assertIn("not a force law", self.record["scope_note"])

    def test_the_larger_force_dominates_the_collar(self):
        """Amplitude dominance is unambiguous even where reach is not."""
        self.assertGreater(self.record["summary"]["dominant_fraction_inner_min"], 0.9)

    def test_reach_is_not_uniformly_longer_for_the_larger_force(self):
        """The honest answer: mostly yes, and it flips."""
        counts = self.record["reach_longer_counts"]
        self.assertGreater(counts.get("radial", 0), counts.get("angular", 0))
        self.assertGreater(counts.get("angular", 0), 0)


if __name__ == "__main__":
    unittest.main()