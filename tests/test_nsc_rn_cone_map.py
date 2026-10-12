"""Bounded checks for the cone map on RN radial states.

Synthetic cases pin the cone facts on inputs whose answer is known in
closed form. Recorded-case checks confirm the reconstruction guard: the
module must refuse a station whose mass does not match the stored
coupled mass, and it must reproduce the march's route on a real payload
when one exists.
"""
import unittest
from pathlib import Path

import numpy as np

from recursive_horizons import nsc_rn_cone_map as cone

LAB = Path(__file__).resolve().parents[1]
FAMILY = LAB / "results/development/nsc-rn-neutral-family-v1"


def _profile(radius, lapse, shift):
    return cone.cone_profile(
        {"radius": np.asarray(radius, dtype=float)},
        {"lapse": np.asarray(lapse, dtype=float),
         "shift": np.asarray(shift, dtype=float),
         "shift_factor": np.asarray(np.array(shift) / np.array(lapse), dtype=float)},
    )


class ConeIdentityTests(unittest.TestCase):
    def test_speed_identities(self):
        radius = np.linspace(1.0, 2.0, 33)
        lapse = 1.0 - 0.15 * ((radius - 1.0) ** 2)
        shift = 0.2 * radius
        profile = _profile(radius, lapse, shift)
        self.assertTrue(np.allclose(profile["v_out"] + profile["v_in"], -2.0 * shift, atol=1e-12))
        self.assertTrue(np.allclose(profile["v_out"] - profile["v_in"], 2.0 * lapse, atol=1e-12))

    def test_local_speed_statement_is_one(self):
        profile = _profile(np.linspace(1.0, 2.0, 9), np.full(9, 0.7), np.linspace(0.1, 0.4, 9))
        self.assertEqual(profile["local_speed_everywhere"], 1.0)

    def test_linear_trap_crossing_is_exact(self):
        radius = np.linspace(1.0, 2.0, 65)
        lapse = np.ones_like(radius)
        shift = 2.4 - 0.8 * radius  # v_out = 0.8 r - 1.4: negative inside, zero at r = 1.75
        profile = _profile(radius, lapse, shift)
        self.assertAlmostEqual(profile["trap_radius"], 1.75, places=12)
        self.assertGreater(profile["trapped_below_index"], 0)

    def test_no_crossing_raises_unresolved_flag(self):
        radius = np.linspace(1.0, 2.0, 9)
        profile = _profile(radius, np.ones_like(radius), np.full_like(radius, 0.2))
        self.assertTrue(np.isnan(profile["trap_radius"]))

    def test_cumulative_crossing_matches_trapezoid(self):
        radius = np.linspace(1.0, 4.0, 129)
        lapse = np.ones_like(radius)
        shift = np.full_like(radius, 0.3)
        profile = _profile(radius, lapse, shift)
        expected = np.trapezoid(1.0 / (-profile["v_in"]), radius)
        self.assertTrue(np.isclose(profile["crossing_coordinate_time"][0], expected, rtol=1e-12))
        self.assertTrue(np.isclose(profile["light_delay_vs_flat"][0],
                                   expected - (radius[-1] - radius[0]), rtol=1e-12))

    def test_redshift_is_inverse_lapse(self):
        radius = np.linspace(1.0, 2.0, 17)
        lapse = np.exp(-0.25 * (radius - 1.0))
        profile = _profile(radius, lapse, np.zeros_like(radius))
        self.assertTrue(np.allclose(profile["ray_frequency_ratio"], 1.0 / lapse, atol=1e-12))
        self.assertAlmostEqual(profile["redshift_at_excision"],
                               1.0 / lapse[0], places=12)


class GridAndReferenceTests(unittest.TestCase):
    def test_station_grid_is_production_grid(self):
        radius = np.linspace(1.04, 32.0, 97)
        grid = cone.station_grid({"radius": radius})
        self.assertEqual(grid["points"], 97)
        self.assertTrue(np.all(np.asarray(grid["weights"]) > 0.0))
        self.assertEqual(float(grid["charge"]), 1.0)
        self.assertAlmostEqual(grid["spacing"], radius[1] - radius[0], places=12)

    def test_station_grid_rejects_nonuniform_radius(self):
        with self.assertRaises(ValueError):
            cone.station_grid({"radius": np.array([1.0, 1.1, 1.25])})

    def test_reference_r_plus(self):
        for mu in (1.002, 1.01, 1.04):
            self.assertAlmostEqual(
                cone.reference_r_plus(mu), mu + np.sqrt(mu * mu - 1.0), places=12)

    def test_reference_r_plus_rejects_unphysical(self):
        with self.assertRaises(ValueError):
            cone.reference_r_plus(0.9)


class ReconstructionGuardTests(unittest.TestCase):
    def _vacuum_payload(self, mismatched=False):
        radius = np.linspace(1.8, 2.2, 12)
        phi = np.zeros((2, 12, 1), dtype=complex)
        mass = np.full(12, 1.5)
        if mismatched:
            mass = mass + 1.0e-6
        return {
            "radius": radius,
            "coupled_phi": phi,
            "occupations": np.array([0.5]),
            "coupled_inner_mass": np.array([1.5]),
            "coupled_mass": mass,
            "coupled_t": np.array([0.0]),
        }

    def test_vacuum_mass_residual_is_roundoff(self):
        payload = self._vacuum_payload()
        grid = cone.station_grid(payload)
        metric = cone.coupled_metric(payload, grid)
        self.assertTrue(metric["reconstruction_is_the_state_metric"])
        self.assertLess(metric["mass_residual"], 1e-12)

    def test_mismatched_mass_is_refused(self):
        payload = self._vacuum_payload(mismatched=True)
        with self.assertRaises(ValueError):
            cone.station_assessment(payload)

    def test_vacuum_redshift_dilution_absent(self):
        payload = self._vacuum_payload()
        assessment = cone.station_assessment(payload)
        self.assertEqual(assessment["verdict"], "cone_map_measured")
        self.assertTrue(np.isclose(assessment["profile"]["redshift_at_excision"], 1.0, atol=1e-12))


@unittest.skipUnless((FAMILY / "prepare.json").is_file() and (FAMILY / "mu-1.010-eps-0.001/station-0.npz").is_file(),
                     "recorded RN family payloads are not present")
class RecordedPayloadTests(unittest.TestCase):
    def test_real_station_reconstructs_the_stored_mass(self):
        with np.load(FAMILY / "mu-1.010-eps-0.001/station-0.npz") as payload:
            assessment = cone.station_assessment(payload, mass_over_rm=1.01)
        self.assertLess(assessment["mass_residual"], cone.MASS_RECONSTRUCTION_TOLERANCE)
        self.assertEqual(assessment["verdict"], "cone_map_measured")
        profile = assessment["profile"]
        self.assertTrue(np.isfinite(profile["trap_radius"]))
        self.assertLess(profile["chart_speeds_at_inner"][0], 0.0)
        self.assertGreater(profile["chart_speeds_at_outer"][0], 0.0)
        self.assertTrue(np.isfinite(profile["redshift_at_excision"]))
        self.assertGreater(profile["redshift_at_excision"], 0.0)


if __name__ == "__main__":
    unittest.main()