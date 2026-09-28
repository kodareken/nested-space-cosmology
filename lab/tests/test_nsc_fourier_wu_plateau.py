import unittest

import numpy as np

from recursive_horizons.nsc_evolved_incoming_constraints import compatible_history_slots
from recursive_horizons.nsc_local_incoming_family import LocalIncomingFamily, plateau_derivatives
from recursive_horizons.nsc_local_incoming_fourier_family import (
    ALLOWED_FOURIER_COUNTS, LocalFourierAxialFunction, LocalFourierIncomingFamily,
    fourier_mode_kinds, fourier_wave_derivatives, project_callable_on_fourier,
)


class FourierWuPlateauTests(unittest.TestCase):
    def test_known_cosine_recovers_itself_on_I(self):
        center, inner = 0.0, 0.03
        kinds = fourier_mode_kinds(8)
        index = kinds.index(('cos', 1))
        coeff = np.zeros(8)
        coeff[index] = 1.3
        function = LocalFourierAxialFunction(coeff, center)
        for z in np.linspace(center - inner, center + inner, 11):
            expected = fourier_wave_derivatives(z, center, inner, 'cos', 1)
            for order in range(5):
                self.assertAlmostEqual(function(z, order), 1.3 * expected[order], delta=2e-12)
        for z in (-0.07, -0.06, 0.06, 0.07):
            self.assertEqual([function(z, n) for n in range(5)], [0.0] * 5)

    def test_leibniz_rule_in_the_collar(self):
        center, inner, outer = 0.0, 0.03, 0.06
        kinds = fourier_mode_kinds(8)
        index = kinds.index(('sin', 2))
        coeff = np.zeros(8)
        coeff[index] = -0.4
        function = LocalFourierAxialFunction(coeff, center)
        z = 0.043
        window = plateau_derivatives(z, center, inner, outer)
        series = fourier_wave_derivatives(z, center, inner, 'sin', 2) * -0.4
        from math import comb
        for n in range(5):
            expected = sum(comb(n, k) * window[k] * series[n - k] for k in range(n + 1))
            self.assertAlmostEqual(function(z, n), expected, delta=2e-12)

    def test_metric_slots_are_derivatives_of_same_two_functions(self):
        coeff = np.zeros((2, 8))
        coeff[0, 0] = 0.02
        coeff[0, 1] = -0.01
        coeff[1, 2] = 1.5
        family = LocalFourierIncomingFamily(coeff)
        z = family.collocation_nodes()
        slots, tangent = compatible_history_slots(family.metric(), z, 16)
        w, U = family.functions
        expected = np.array([[w(v, 0), w(v, 1), w(v, 2), w(v, 3), U(v, 0), U(v, 1)] for v in z])
        np.testing.assert_allclose(slots, expected, atol=3e-11, rtol=2e-15)
        direction, step = 3, 1e-5
        altered = coeff.copy()
        altered.ravel()[direction] += step
        plus, _ = compatible_history_slots(LocalFourierIncomingFamily(altered).metric(), z, 16)
        np.testing.assert_allclose((plus - slots) / step, tangent[direction], atol=2e-7, rtol=2e-10)
        self.assertGreater(family.radius_lower_bound(), 0.)
        self.assertGreater(family.interval[1] - family.interval[0], 0.)
        self.assertEqual(ALLOWED_FOURIER_COUNTS, (8, 16))

    def test_projection_recovers_a_known_pair_and_rejects_chebyshev_width(self):
        center, inner = 0.0, 0.03
        kinds = fourier_mode_kinds(8)
        target = np.zeros(8)
        target[kinds.index(('const', 0))] = 0.2
        target[kinds.index(('sin', 1))] = -0.15
        recovered = project_callable_on_fourier(
            lambda z: LocalFourierAxialFunction(target, center)(z, 0),
            center, inner, 8)
        np.testing.assert_allclose(recovered, target, atol=2e-12)
        with self.assertRaises(ValueError):
            LocalFourierIncomingFamily(np.zeros((2, 7)))
        with self.assertRaises(ValueError):
            LocalFourierIncomingFamily(np.zeros((2, 32)))
        self.assertIn('Chebyshev', LocalIncomingFamily(np.zeros((2, 8))).description()['coefficient_order'])


if __name__ == '__main__':
    unittest.main()
