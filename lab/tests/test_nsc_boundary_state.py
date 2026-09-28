"""One new integration slice: correlated boundary reconstruction.

The existing smooth-cell Hamiltonian is an input. No old source, spectrum
scan, regulator comparison or publication generator is rerun.
"""
import unittest

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import expm

from recursive_horizons.nsc_boundary_state import GaussianBoundaryState
from recursive_horizons.nsc_covariant_operator import smooth_metric
from recursive_horizons.nsc_influence import (
    canonical_hamiltonian, ground_covariance, influence,
)


class BoundaryStateIntegration(unittest.TestCase):
    def test_correlated_history_and_action_on_existing_smooth_cell(self):
        metric = smooth_metric(12, general=True)
        h = canonical_hamiltonian(metric)
        _, _, vacuum = ground_covariance(h)
        p = np.tile(metric.x >= 0, 2)
        # A unitary relative-phase preparation retains the pure covariance
        # and its cross-room correlations; it does not select a new spectrum.
        phase = np.exp(.37j * p)
        c = phase[:, None] * vacuum * phase.conj()[None, :]
        state = GaussianBoundaryState(c, p)
        b, hc = h[np.ix_(p, ~p)], h[np.ix_(~p, ~p)]
        t, s = .63, .27
        ut, us = expm(-1j * h * t), expm(-1j * h * s)
        # Independently assemble V(t) from its causal convolution rather
        # than taking the off-diagonal block of full U(t).
        nodes, weights = leggauss(32)
        def convolution(time):
            result = np.zeros_like(b)
            for x, w in zip((nodes + 1) * time / 2, weights * time / 2):
                a = expm(-1j * h * (time - x))[np.ix_(p, p)]
                result += -1j * w * a @ b @ expm(-1j * hc * x)
            return result
        vt, vs = convolution(t), convolution(s)
        at, ass = ut[np.ix_(p, p)], us[np.ix_(p, p)]
        reconstructed = state.reconstruct(at, vt, ass, vs)
        direct = (ut @ c @ us.conj().T)[np.ix_(p, p)]
        covariance_error = np.max(abs(reconstructed['total'] - direct))
        self.assertLess(covariance_error, 2e-12)
        missing_preparation = reconstructed['parent'] + reconstructed['child']
        omission = np.linalg.norm(missing_preparation - direct)
        self.assertGreater(omission, .01)

        k = state.kernels(t, s, link_t=b, child_t=expm(-1j * hc * t),
                          link_s=b, child_s=expm(-1j * hc * s))
        np.testing.assert_allclose(k['retarded'], k['greater'] - k['lesser'], atol=2e-12)
        ft = b @ expm(-1j * hc * t)
        np.testing.assert_allclose(k['initial_occupied'], ft @ c[np.ix_(~p, p)])

        # Endpoint Grassmann elimination retains correlated preparation
        # even for a pure state, without inverting C or I-C.
        probe = np.diag(np.cos(np.tile(metric.x, 2)))
        plus = expm(-.4j * (h + .17 * probe))
        minus = expm(-.4j * (h - .11 * probe))
        factored = state.overlap_schur(plus, minus)
        canonical = influence(c, plus, minus)
        action_error = abs(factored['amplitude'] - canonical['amplitude'])
        self.assertLess(action_error, 2e-12)
        self.assertLess(abs(factored['principal_action'] - canonical['principal_action']), 2e-12)
        print(f'Boundary reconstruction max error: {covariance_error:.9g}; '
              f'initial-correlation omission norm: {omission:.9g}; '
              f'normalized-action factorization error: {action_error:.9g}')


if __name__ == '__main__':
    unittest.main()
