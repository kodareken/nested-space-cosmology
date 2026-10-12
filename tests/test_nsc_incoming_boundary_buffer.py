"""Controls of the new extension/frame/phase path, without archived field runs."""
import unittest

import numpy as np
from scipy.linalg import expm
from scipy.sparse import bmat, csr_matrix, diags

from recursive_horizons.nsc_incoming_boundary_buffer import (
    append_extension, characteristic_from_ks_path, continuum_incoming_buffer,
    exact_harmonic_reference, extended_grid, ks_from_characteristic,
    original_grid, prefixes_equal,
)
from recursive_horizons.nsc_transmitting_history_jets import (
    FourthOrderModePropagator, incoming_reference_columns,
)


class IncomingBoundaryBufferTests(unittest.TestCase):
    def test_prefix_bytes_and_spin_layout(self):
        for old_points, points in ((801, 951), (3201, 3801)):
            old, new = original_grid(old_points), extended_grid(points)
            phi = np.arange(2*old_points*12).reshape(2*old_points, 12).astype(complex)
            extension = -np.arange(2*(points-old_points)*12).reshape(2*(points-old_points),12)
            joined = append_extension(old, phi, new, extension)
            self.assertTrue(prefixes_equal(old, phi, new, joined))
            np.testing.assert_array_equal(joined.reshape(2, points, 12)[:, old_points:],
                                          extension.reshape(2, points-old_points, 12))
            changed = new.copy(); changed[2] += 1e-12
            with self.assertRaises(ValueError):
                append_extension(old, phi, changed, extension)

    def test_inverse_restriction_clock_and_source_order(self):
        rng = np.random.default_rng(52)
        E = np.array([-.97, -.55, .55, .97])
        rho = np.array([1.2, 1.8])
        A = rng.normal(size=(2,4,2,3))+1j*rng.normal(size=(2,4,2,3))
        phi = characteristic_from_ks_path(rho, E, A)
        for index in range(2):
            actual, _ = ks_from_characteristic(rho, phi, E, index)
            np.testing.assert_allclose(actual, A[index], atol=3e-14, rtol=0)

    def test_harmonic_block_and_nodal_pde_against_dense_exponential(self):
        x = np.linspace(-2., 1.2, 17)
        owner = FourthOrderModePropagator(x, .7, .8)
        rng = np.random.default_rng(11)
        phi = rng.normal(size=(34,12))+1j*rng.normal(size=(34,12))
        E = np.array([-1.3, -.4, .4, 1.3])
        times = np.linspace(0., .0002, 3)
        got = exact_harmonic_reference(owner, phi, E, times)
        L, _ = owner.generator(owner.reference)
        B = incoming_reference_columns(owner, phi)
        M = bmat([[L, csr_matrix(B)], [None, diags(-1j*np.repeat(E,3))]]).toarray()
        initial = np.vstack((phi, np.eye(12)))
        take = got['sample_rows']; R = got['restriction']
        for k, t in enumerate(times):
            field = expm(t*M)@initial
            np.testing.assert_allclose(got['columns'][k], R@field[take], atol=2e-13, rtol=0)
            np.testing.assert_allclose(got['pde_axial_columns'][k], R@(M@field)[take], atol=3e-11, rtol=0)
        np.testing.assert_allclose(got['phase'], got['expected_phase'], atol=2e-14, rtol=0)
        with self.assertRaises(ValueError):
            exact_harmonic_reference(owner, phi, E, np.array([0., .001, .003]))

    def test_continuum_buffer_scope(self):
        result = continuum_incoming_buffer()
        self.assertGreater(result['causal_margin_lower_bound'], 0.)
        self.assertGreater(result['travel_time_quadrature'], result['travel_time_lower_bound'])
        self.assertFalse(result['finite_grid_exact_causality'])


if __name__ == '__main__':
    unittest.main()
