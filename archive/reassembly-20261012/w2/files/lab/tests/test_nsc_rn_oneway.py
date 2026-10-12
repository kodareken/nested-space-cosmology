"""Bounded checks for the one-way enclosure measure.

The operator matrix must reproduce ``nsc_rn_pg.dirac_rates`` exactly, the
direction statistic is bounded in [-1, 1] and agrees with hand-set
vectors, and the interior weight partitions the field. A recorded-station
check runs only when the family payloads exist.
"""
import unittest
from pathlib import Path

import numpy as np

from recursive_horizons import nsc_rn_cone_map as cone
from recursive_horizons import nsc_rn_oneway as oneway
from recursive_horizons import nsc_rn_pg as pg

LAB = Path(__file__).resolve().parents[1]
FAMILY = LAB / "results/development/nsc-rn-neutral-family-v1"


def _small_grid(points=24):
    radius = np.linspace(1.04, 3.0, points)
    return cone.station_grid({"radius": radius}), radius


def _geometry(points):
    radius = np.linspace(1.04, 3.0, points)
    lapse = 1.0 - 0.05 * ((radius - 1.04) ** 2)
    shift = 0.3 * radius
    return radius, lapse, shift


class OperatorMatrixTests(unittest.TestCase):
    def test_matrix_reproduces_dirac_rates(self):
        grid, radius = _small_grid()
        lapse = np.full_like(radius, 0.9)
        shift = 0.25 * radius
        matrix = oneway.operator_matrix(grid, lapse, shift)
        rng = np.random.default_rng(3)
        plus = rng.normal(size=radius.shape[0])
        minus = rng.normal(size=radius.shape[0])
        vector = np.concatenate([plus, minus])
        phi = np.zeros((2, radius.shape[0], 1))
        phi[0, :, 0] = plus
        phi[1, :, 0] = minus
        rates, _sat = pg.dirac_rates(grid, phi, lapse, shift, kappa=1)
        expected = np.concatenate([rates[0, :, 0], rates[1, :, 0]])
        self.assertTrue(np.allclose(matrix @ vector, expected, atol=1e-12, rtol=1e-12))


class DirectionStatisticTests(unittest.TestCase):
    def test_outward_vector_is_positive(self):
        radius = np.linspace(1.04, 3.0, 12)
        v_out = np.full_like(radius, 0.7)
        v_in = np.full_like(radius, -1.7)
        plus = np.ones_like(radius)  # all weight in the outward channel
        minus = np.zeros_like(radius)
        value = oneway._direction_statistic(v_out, v_in, plus, minus)
        self.assertAlmostEqual(value, 1.0, places=12)

    def test_inward_vector_is_negative(self):
        radius = np.linspace(1.04, 3.0, 12)
        v_out = np.full_like(radius, 0.7)
        v_in = np.full_like(radius, -1.7)
        plus = np.zeros_like(radius)
        minus = np.ones_like(radius)
        value = oneway._direction_statistic(v_out, v_in, plus, minus)
        self.assertAlmostEqual(value, -1.0, places=12)

    def test_trapped_channel_content_counts_as_inward(self):
        # Inside the trap both characteristics are negative; the plus
        # channel is inward there, so a pure trapped plus-mode has a
        # negative signed current. The v1 magnitude statistic read this
        # same configuration as outward; v2 pins the correction.
        radius = np.linspace(1.04, 3.0, 12)
        v_out = np.full_like(radius, -0.5)
        v_in = np.full_like(radius, -1.5)
        plus = np.ones_like(radius)
        minus = np.zeros_like(radius)
        value = oneway._direction_statistic(v_out, v_in, plus, minus)
        self.assertAlmostEqual(value, -1.0, places=12)

    def test_balanced_vector_is_zero(self):
        radius = np.linspace(1.04, 3.0, 12)
        v_out = np.full_like(radius, 0.7)
        v_in = np.full_like(radius, -0.7)
        plus = np.ones_like(radius)
        minus = np.ones_like(radius)
        value = oneway._direction_statistic(v_out, v_in, plus, minus)
        self.assertAlmostEqual(value, 0.0, places=12)

    def test_vanishing_field_is_zero_not_nan(self):
        radius = np.linspace(1.04, 3.0, 8)
        value = oneway._direction_statistic(
            np.full_like(radius, 0.7), np.full_like(radius, -0.7),
            np.zeros_like(radius), np.zeros_like(radius))
        self.assertEqual(value, 0.0)

    def test_interior_weight_partitions(self):
        radius = np.linspace(1.0, 2.0, 16)
        plus = np.ones_like(radius)
        minus = np.zeros_like(radius)
        weight = oneway._interior_weight(radius, 1.5, plus, minus)
        self.assertAlmostEqual(weight, 8.0 / 16.0, places=9)


@unittest.skipUnless((FAMILY / "prepare.json").is_file() and (FAMILY / "mu-1.010-eps-0.001/station-0.npz").is_file(),
                     "recorded RN family payloads are not present")
class RecordedOnewayTests(unittest.TestCase):
    def test_real_station_measures(self):
        with np.load(FAMILY / "mu-1.010-eps-0.001/station-0.npz") as payload:
            assessment = oneway.station_oneway(payload, mass_over_rm=1.01)
        self.assertEqual(assessment["verdict"], "oneway_measured")
        self.assertGreater(assessment["mode_count"], 0)
        self.assertTrue(np.isfinite(assessment["interior_max_direction"]))
        self.assertTrue(np.isfinite(assessment["exterior_min_direction"]))
        self.assertLessEqual(assessment["interior_max_direction"], 1.0 + 1e-9)
        self.assertGreaterEqual(assessment["exterior_min_direction"], -1.0 - 1e-9)


if __name__ == "__main__":
    unittest.main()