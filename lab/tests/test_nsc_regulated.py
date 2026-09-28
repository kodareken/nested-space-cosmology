import unittest
import numpy as np
from scipy.linalg import expm
from recursive_horizons.nsc_regulated import (
    OperatorConventions, RegulatedOperator, recursive_response, direct_chain, pg_principal_witness,
)


class RegulatedTests(unittest.TestCase):
    def test_degenerate_spectrum_hessian_and_basis_covariance(self):
        d = np.diag([.8, .8, -1.2])
        perturbation = np.array([[.1, .2j, .4], [-.2j, -.2, .1], [.4, .1, .3]])
        op = RegulatedOperator(d)
        gradient, hessian = op.variation(perturbation)
        eps = .0002
        finite = (RegulatedOperator(d+eps*perturbation).action()-2*op.action()+RegulatedOperator(d-eps*perturbation).action())/eps**2
        self.assertAlmostEqual(finite, hessian, delta=2e-7)
        u = expm(1j*perturbation)
        rotated = RegulatedOperator(u@d@u.conj().T).variation(u@perturbation@u.conj().T)
        np.testing.assert_allclose(rotated, [gradient, hessian], atol=1e-12)

    def test_zero_mode_and_nonhermitian_domain_rejected(self):
        with self.assertRaisesRegex(ValueError, 'zero'):
            RegulatedOperator([[0., 0.], [0., 1.]])
        with self.assertRaisesRegex(ValueError, 'Hermitian'):
            RegulatedOperator([[1., 1.], [0., 1.]])

    def test_raw_determinant_limit_keeps_normalization(self):
        d = np.diag([.5, -1.2])
        op = RegulatedOperator(d, OperatorConventions(cutoff=1e5, normalization=.7))
        self.assertAlmostEqual(op.action(), -np.log(abs(np.linalg.det(d/.7))), places=8)

    def test_chain_schur_order_and_normalization(self):
        h = np.array([[.7, .2j], [-.2j, 1.1]])
        b = np.array([[.2j, .1], [.05, .3]])
        for omega in (1., 1.7):
            np.testing.assert_allclose(recursive_response(h,b,.4+.1j,omega,8),
                                       direct_chain(h,b,.4+.1j,omega,8),atol=1e-12)

    def test_trapped_symbol_has_nonzero_characteristic_spatial_covector(self):
        result = pg_principal_witness()
        self.assertLess(result['A'], 0)
        self.assertGreater(np.linalg.norm(result['covector']), 0)
        self.assertAlmostEqual(result['principal_eigenvalues'][1], 0., places=12)
        self.assertFalse(result['covariant_Lorentzian_Dirac_ill_posed_inferred'])


if __name__ == '__main__':
    unittest.main()
