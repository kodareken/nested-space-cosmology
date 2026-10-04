"""Bounded checks for the two-force gradient balance on RN radial states.

The synthetic cases pin the numerics on inputs whose answer is known in
closed form, and the recorded-case checks confirm the provenance guard
that keeps the analytic source-free reference out of the force ratio.
"""
import json
import unittest
from pathlib import Path

import numpy as np

from recursive_horizons import nsc_rn_gradient_bubble as gb

LAB = Path(__file__).resolve().parents[1]
RECORD = LAB / "results/development/nsc-rn-gradient-bubble-v2.json"


def _payload(inner, outer, radius=None):
    inner = np.asarray(inner, dtype=float)
    if radius is None:
        radius = np.linspace(1.04, 32.0, inner.size)
    return {"radius": np.asarray(radius, dtype=float),
            "radial_pressure": inner,
            "angular_pressure": np.asarray(outer, dtype=float)}


class LogSlopeTests(unittest.TestCase):
    def test_power_law_slope_is_exact(self):
        radius = np.linspace(1.04, 32.0, 401)
        for exponent in (-1.0, -3.5, -7.0):
            with self.subTest(exponent=exponent):
                force = 2.0 * radius ** exponent
                slope = gb._log_slope(radius, force)
                interior = slope[1:-1]
                self.assertTrue(np.allclose(interior, exponent, rtol=1e-8, atol=1e-8))

    def test_slope_is_nan_when_a_neighbour_vanishes(self):
        # The centred stencil at i reads force[i-1] and force[i+1]. A zero at
        # the centre does not make the slope undefined; a zero in the
        # stencil does.
        radius = np.linspace(1.0, 2.0, 5)
        force = np.array([1.0, 0.0, 1.0, 1.0, 1.0])
        slope = gb._log_slope(radius, force)
        self.assertFalse(np.isnan(slope[1]))
        self.assertAlmostEqual(slope[1], 0.0)

        hollow = np.array([1.0, 1.0, 0.0, 1.0, 1.0])
        hollow_slope = gb._log_slope(radius, hollow)
        self.assertTrue(np.isnan(hollow_slope[1]))
        self.assertFalse(np.isnan(hollow_slope[2]))
        self.assertTrue(np.isnan(hollow_slope[3]))

    def test_endpoints_are_undefined(self):
        radius = np.linspace(1.0, 2.0, 9)
        slope = gb._log_slope(radius, radius ** -2.0)
        self.assertTrue(np.isnan(slope[0]))
        self.assertTrue(np.isnan(slope[-1]))


class BalanceTests(unittest.TestCase):
    def test_equal_forces_give_one_half(self):
        force = np.linspace(1e-6, 1e-8, 64)
        profile = gb.gradient_bubble_profile(_payload(force, force.copy()))
        self.assertTrue(np.allclose(profile["balance"], 0.5))

    def test_absent_outer_force_gives_one(self):
        force = np.linspace(1e-6, 1e-8, 64)
        profile = gb.gradient_bubble_profile(_payload(force, np.zeros_like(force)))
        self.assertTrue(np.allclose(profile["balance"], 1.0))

    def test_integrated_force_matches_analytic_sphere(self):
        # 4 pi int_1^2 r^2 dr = 4 pi (8-1)/3
        radius = np.linspace(1.0, 2.0, 200001)
        value = gb._integrated_force(radius, np.ones_like(radius))
        self.assertAlmostEqual(value, 4.0 * np.pi * 7.0 / 3.0, places=4)

    def test_sign_change_counter(self):
        self.assertEqual(gb._sign_changes(np.array([1.0, 2.0, -1.0, -2.0, 5.0])), 2)
        self.assertEqual(gb._sign_changes(np.array([1.0, 2.0, 3.0])), 0)


class WallTests(unittest.TestCase):
    """The verdict basis: a measured wall, not a chosen slope limit."""

    def _profile(self, balance_values, radius=None):
        balance = np.asarray(balance_values, dtype=float)
        if radius is None:
            radius = np.linspace(1.04, 32.0, balance.size)
        return gb._balance_walls(radius, balance, np.ones(balance.size, dtype=bool))

    def test_single_wall_is_located(self):
        # Radial above one half for the first third, angular below after.
        radius = np.linspace(1.04, 32.0, 300)
        balance = np.where(radius < 10.0, 0.8, 0.2)
        walls = gb._balance_walls(radius, balance, np.ones(radius.size, dtype=bool))
        self.assertEqual(walls.size, 1)
        self.assertAlmostEqual(float(walls[0]), 10.0, delta=0.2)

    def test_single_cell_wobble_is_not_a_wall(self):
        radius = np.linspace(1.04, 32.0, 300)
        balance = np.full(radius.size, 0.8)
        balance[100] = 0.2
        walls = gb._balance_walls(radius, balance, np.ones(radius.size, dtype=bool))
        self.assertEqual(walls.size, 0)

    def test_crossing_outside_the_live_region_is_ignored(self):
        radius = np.linspace(1.04, 32.0, 300)
        balance = np.where(radius < 10.0, 0.8, 0.2)
        # Live only below the crossing, so the wall cannot be believed.
        live = np.zeros(radius.size, dtype=bool)
        live[:25] = True
        walls = gb._balance_walls(radius, balance, live)
        self.assertEqual(walls.size, 0)

    def test_verdict_is_two_force_wall_when_a_wall_exists(self):
        radius = np.linspace(1.04, 32.0, 300)
        # Constant magnitudes keep both forces live out to the crossing.
        inner = np.where(radius < 10.0, 2.0, 0.2)
        outer = np.full(radius.size, 0.5)
        verdict = gb.gradient_bubble_assessment(_payload(inner, outer, radius))
        self.assertEqual(verdict["verdict"], "two_force_wall")
        self.assertEqual(verdict["profile"]["wall_count"], 1)
        self.assertAlmostEqual(verdict["profile"]["r_first_wall"], 10.0, delta=0.2)

    def test_verdict_is_single_force_collar_without_a_wall(self):
        radius = np.linspace(1.04, 32.0, 300)
        # The outer force is present but never catches the radial one, so
        # there is no wall and yet this is not a missing-force dead end.
        inner = np.full(radius.size, 2.0)
        outer = np.full(radius.size, 1.0)
        verdict = gb.gradient_bubble_assessment(_payload(inner, outer, radius))
        self.assertEqual(verdict["profile"]["wall_count"], 0)
        self.assertEqual(verdict["dead_end_flags"], [])
        self.assertEqual(verdict["verdict"], "single_force_collar")

    def test_negligible_outer_force_is_a_hard_dead_end(self):
        radius = np.linspace(1.04, 32.0, 300)
        inner = np.full(radius.size, 2.0)
        outer = np.full(radius.size, 0.01)
        verdict = gb.gradient_bubble_assessment(_payload(inner, outer, radius))
        self.assertIn("single_force", verdict["dead_end_flags"])
        self.assertEqual(verdict["verdict"], "dead_end")

    def test_steep_throat_is_reported_but_does_not_decide(self):
        """The v1 mistake: a steep slope used to override the structure."""
        radius = np.linspace(1.04, 32.0, 300)
        steep = radius ** -30.0
        paired = np.where(radius < 10.0, radius ** -3.0, radius ** -9.0)
        # A steep radial force with a genuine wall crossing it.
        verdict = gb.gradient_bubble_assessment(_payload(steep * 0 + paired, paired, radius))
        steep_only = gb.gradient_bubble_assessment(
            _payload(steep, 0.6 * steep, radius))
        self.assertIn("throat_steep_reported", steep_only["dead_end_flags"])
        self.assertNotIn("throat_steep_reported", verdict["dead_end_flags"])


class DeadEndTests(unittest.TestCase):
    def test_single_force_flag_when_outer_is_negligible(self):
        radius = np.linspace(1.04, 32.0, 401)
        inner = radius ** -4.0
        outer = 1e-9 * radius ** -4.0
        verdict = gb.gradient_bubble_assessment(_payload(inner, outer, radius))
        self.assertEqual(verdict["verdict"], "dead_end")
        self.assertIn("single_force", verdict["dead_end_flags"])

    def test_throat_runaway_flag_on_steep_profile(self):
        radius = np.linspace(1.04, 32.0, 401)
        inner = radius ** -30.0
        outer = 0.5 * radius ** -30.0
        verdict = gb.gradient_bubble_assessment(_payload(inner, outer, radius))
        self.assertIn("throat_steep_reported", verdict["dead_end_flags"])
        self.assertNotEqual(verdict["verdict"], "two_force_wall")

    def test_balanced_verdict_when_both_forces_act(self):
        radius = np.linspace(1.04, 32.0, 401)
        inner = np.where(radius < 10.0, 2.0, 0.2)
        outer = np.full(radius.size, 0.5)
        verdict = gb.gradient_bubble_assessment(_payload(inner, outer, radius))
        self.assertEqual(verdict["dead_end_flags"], [])
        self.assertEqual(verdict["verdict"], "two_force_wall")

    def test_curvature_dominates_uses_same_provenance(self):
        radius = np.linspace(1.04, 32.0, 401)
        inner = radius ** -5.0
        outer = 0.6 * radius ** -6.0
        payload = _payload(inner, outer, radius)
        benign = gb.gradient_bubble_assessment(payload, {"R4": 1e-12 * np.ones_like(radius)})
        loud = gb.gradient_bubble_assessment(payload, {"R4": 1e4 * np.ones_like(radius)})
        self.assertNotIn("curvature_dominates", benign["dead_end_flags"])
        self.assertIn("curvature_dominates", loud["dead_end_flags"])

    def test_reference_ricci2_is_rejected_by_name(self):
        """The analytic source-free reference must not enter the force ratio."""
        radius = np.linspace(1.04, 32.0, 401)
        inner = radius ** -5.0
        outer = 0.6 * radius ** -6.0
        payload = _payload(inner, outer, radius)
        reference = 4.0 * radius ** -8.0
        with_ref = gb.gradient_bubble_assessment(payload, {"R4": np.zeros_like(radius),
                                                           "Ricci2": reference})
        without = gb.gradient_bubble_assessment(payload)
        self.assertEqual(with_ref["curvature_force_ratio"], without["curvature_force_ratio"])
        self.assertIn("Ricci2", with_ref["curvature_keys_rejected"])


class InputGuardTests(unittest.TestCase):
    def test_non_increasing_grid_refused(self):
        radius = np.array([1.04, 1.02, 1.0])
        with self.assertRaises(ValueError):
            gb.gradient_bubble_profile(_payload(np.ones(3), np.ones(3), radius))

    def test_shape_mismatch_refused(self):
        with self.assertRaises(ValueError):
            gb.gradient_bubble_profile({"radius": np.ones(5),
                                        "radial_pressure": np.ones(4),
                                        "angular_pressure": np.ones(5)})

    def test_short_grid_refused(self):
        with self.assertRaises(ValueError):
            gb.gradient_bubble_profile(_payload(np.ones(2), np.ones(2)))


@unittest.skipUnless(RECORD.is_file(), "gradient-bubble record not written yet")
class RecordTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record = json.loads(RECORD.read_text(encoding="utf-8"))

    def test_record_is_creation_only(self):
        self.assertTrue(self.record["creation_only"])

    def test_record_reflects_the_incomplete_family(self):
        coverage = self.record["coverage"]
        self.assertEqual(coverage["cases"], 6)
        self.assertGreater(coverage["total_stations"], 0)
        self.assertFalse(coverage["family_complete"])

    def test_scope_note_disclaims_a_cosmological_claim(self):
        self.assertIn("does not measure", self.record["scope_note"])

    def test_outer_force_is_never_negligible_in_any_case(self):
        ratios = [c["stations"][-1]["force_ratio_integrated"] for c in self.record["cases"]]
        self.assertTrue(all(r > 0.1 for r in ratios), ratios)

    def test_record_supersedes_the_threshold_verdict(self):
        self.assertIsNotNone(self.record["supersedes"])
        self.assertIn("straddled", self.record["supersedes"]["superseded_reason"])

    def test_wall_persistence_falls_as_the_exterior_strengthens(self):
        """The derived relation: a stronger exterior dissolves the wall."""
        by_mu = {}
        for case in self.record["cases"]:
            by_mu.setdefault(case["mass_over_rm"], []).append(
                case["wall_presence_fraction"])
        weak = max(by_mu[1.002])
        strong = max(by_mu[1.04])
        self.assertGreater(weak, strong)
        self.assertGreater(weak, 0.5)
        self.assertLess(strong, 0.5)

    def test_wall_appears_in_every_case_at_some_station(self):
        self.assertEqual(self.record["coverage"]["cases_showing_a_wall_at_any_station"], 6)

    def test_wall_is_never_read_only_at_the_final_station(self):
        """The wall oscillates, so a single snapshot would misreport it."""
        coverage = self.record["coverage"]
        self.assertLess(coverage["cases_with_a_wall_at_the_final_station"], 6)
        self.assertGreater(coverage["station_fraction_with_a_wall"], 0.0)

    def test_no_hard_dead_end_flag_fires(self):
        self.assertEqual(
            sorted({f for f in self.record["dead_end_flags_observed"] if f != "throat_steep_reported"}),
            [])


if __name__ == "__main__":
    unittest.main()