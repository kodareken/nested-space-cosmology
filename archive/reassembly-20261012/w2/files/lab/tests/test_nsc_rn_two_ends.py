"""Bounded checks for the two-ends room-quality measure."""

import unittest
from pathlib import Path

import numpy as np

from recursive_horizons import nsc_rn_two_ends as ends

LAB = Path(__file__).resolve().parents[1]
FAMILY = LAB / "results/development/nsc-rn-neutral-family-v1"


class ClosedFormTests(unittest.TestCase):
    def test_uniform_room_quality_is_inverse_radius(self):
        radius = np.linspace(1.0, 32.0, 129)
        lapse = np.full_like(radius, 1.2)
        shift = np.zeros_like(radius)
        scales = ends.local_room_scales(radius, lapse, shift, spacing=radius[1] - radius[0])
        expected = 1.2 / radius / (1.2 / (radius[1] - radius[0]))
        self.assertTrue(np.allclose(scales["room_quality"], expected, rtol=1e-12))

    def test_quality_is_scale_invariant(self):
        radius = np.linspace(1.0, 16.0, 65)
        spacing = radius[1] - radius[0]
        lapse = 1.0 - 0.1 * ((radius - 1.0) / 15.0)
        shift = 0.2 * radius
        # Self-similar rescale: coordinates double, spacing doubles, and the
        # dimensionless lapse and shift keep the same profile as functions of
        # the coordinate ratio, which here means the identical arrays.
        base = ends.local_room_scales(radius, lapse, shift, spacing)["room_quality"]
        scaled = ends.local_room_scales(2.0 * radius, lapse, shift, 2.0 * spacing)["room_quality"]
        self.assertTrue(np.allclose(base, scaled, rtol=1e-12))

    def test_shift_raises_the_cutoff_end(self):
        radius = np.linspace(1.0, 16.0, 65)
        spacing = radius[1] - radius[0]
        lapse = np.ones_like(radius)
        still = ends.local_room_scales(radius, lapse, np.zeros_like(radius), spacing)
        moving = ends.local_room_scales(radius, lapse, np.full_like(radius, 0.5), spacing)
        self.assertGreater(moving["omega_cutoff"][0], still["omega_cutoff"][0])
        self.assertLess(moving["room_quality"][0], still["room_quality"][0])

    def test_rejects_nonpositive_inputs(self):
        radius = np.linspace(1.0, 2.0, 9)
        with self.assertRaises(ValueError):
            ends.local_room_scales(radius, np.zeros_like(radius), np.zeros_like(radius), spacing=0.1)


@unittest.skipUnless((FAMILY / "prepare.json").is_file() and (FAMILY / "mu-1.010-eps-0.001/station-0.npz").is_file(),
                     "recorded RN family payloads are not present")
class RecordedRoomTests(unittest.TestCase):
    def test_recorded_room_has_two_ends(self):
        with np.load(FAMILY / "mu-1.010-eps-0.001/station-0.npz") as payload:
            assessment = ends.station_two_ends(payload, mass_over_rm=1.01)
        self.assertTrue(assessment["scale_invariant_definition"].startswith("N^2 spacing"))
        self.assertTrue(assessment["dimensionless_invariant_is_per_room"])
        self.assertGreater(assessment["room_quality_inner"], assessment["room_quality_outer"])


if __name__ == "__main__":
    unittest.main()