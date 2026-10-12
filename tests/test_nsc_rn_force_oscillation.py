"""Closed-form crossing and correlation tests; independent numpy checks."""
import unittest
import numpy as np
from recursive_horizons.nsc_rn_force_oscillation import (
    crossing_kernel, association_kernel, profile_summary, spacing_summary)


class OscillationTests(unittest.TestCase):
    def test_dispatchers(self):
        for fn in (crossing_kernel, association_kernel):
            self.assertEqual(type(fn).__name__, 'CPUDispatcher')

    def test_alternating_two_cell_half_wave(self):
        # Samples of cos(pi*j/2), represented exactly (zeros included).
        y = np.tile([1., 0., -1., 0.], 5)[None, :]
        p, c = crossing_kernel(y, 0.)
        self.assertEqual(c[0], 9)
        np.testing.assert_array_equal(p[0, :9], np.arange(1., 18., 2.))
        s = spacing_summary(np.diff(p[0, :9]), .25)
        self.assertEqual(s['mean_cells'], 2.)
        self.assertEqual(s['mean_r_m'], .5)
        self.assertEqual(s['le3_cells'], 8)

    def test_resolved_closed_form_wave(self):
        y = np.cos(np.pi * (np.arange(81) + .25) / 10)[None, :]
        p, c = crossing_kernel(y, 0.)
        self.assertEqual(c[0], 8)
        np.testing.assert_allclose(np.diff(p[0, :8]), 10., atol=1e-13)

    def test_plateau_and_touch(self):
        y = np.array([[1., 0., 0., -1., 0., -1., 0., 1.],
                      [1., 0., 0., 1., 0., 1., 0., 1.]])
        p, c = crossing_kernel(y, 0.)
        np.testing.assert_array_equal(c, [2, 0])
        np.testing.assert_array_equal(p[0, :2], [1.5, 6.])

    def test_linear_interpolation(self):
        p, c = crossing_kernel(np.array([[3., -1., -1., 1.]]), 0.)
        np.testing.assert_allclose(p[0, :c[0]], [.75, 2.5])

    def test_floor_removes_tiny_ringing(self):
        y = np.array([[1., 1e-14, -1e-14, 1e-14, 1.]])
        _, c = crossing_kernel(y, 0.)
        self.assertEqual(c[0], 2)
        _, c = crossing_kernel(y, 1e-12)
        self.assertEqual(c[0], 0)

    def test_all_zero(self):
        _, c = crossing_kernel(np.zeros((2, 10)), 0.)
        np.testing.assert_array_equal(c, [0, 0])
        self.assertIsNone(spacing_summary(np.array([]), .1)['mean_cells'])

    def test_radial_bins(self):
        r = np.linspace(1., 32., 125)
        x = np.array([1., 2., 6., 14., 30., 62., 100.])
        s = profile_summary(x, r)
        self.assertEqual(s['sign_changes'], 7)
        self.assertEqual([a['crossings'] for a in s['regions']], [2, 1, 1, 1, 2])
        self.assertEqual(sum(a['intervals'] for a in s['regions']), 6)
        self.assertEqual(s['mean_cells'], 16.5)
        self.assertEqual(s['mean_r_m'], 4.125)

    def test_correlation_closed_form(self):
        m = np.zeros((6, 20), dtype=np.int64)
        m[0, [2, 10]] = 1
        m[2, [2, 10]] = 1
        m[3, [5, 15]] = 1
        a = association_kernel(m)
        self.assertEqual(a[0, 0, 0], 2)
        self.assertAlmostEqual(a[0, 0, 1], 1.)
        self.assertEqual(a[0, 0, 2], 1.)
        self.assertEqual(a[0, 0, 3], .5)
        self.assertEqual(a[0, 1, 0], 0)
        self.assertAlmostEqual(a[0, 1, 1], -1/9)
        self.assertEqual(a[0, 1, 2], 0)
        self.assertTrue(np.isnan(a[1, 0, 1]))

    def test_independent_numpy_correlations(self):
        rng = np.random.default_rng(2026)
        m = (rng.random((6, 100)) < .2).astype(np.int64)
        a = association_kernel(m)
        for f in range(2):
            for s in range(4):
                self.assertAlmostEqual(a[f, s, 1], np.corrcoef(m[f], m[s+2])[0, 1])


if __name__ == '__main__':
    unittest.main()
