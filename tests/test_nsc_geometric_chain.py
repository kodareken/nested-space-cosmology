import unittest
import numpy as np
from recursive_horizons.nsc_geometric_chain import geometric_motif, continuum_discriminant


class GeometricChainTests(unittest.TestCase):
    def test_removing_angular_potential_closes_zero_phase_gap(self):
        free=geometric_motif(2.,32,kappa=0.)
        confined=geometric_motif(2.,32,kappa=1.)
        def gap(m): return np.min(np.abs(np.linalg.eigvalsh(m['H']+m['B']+m['B'].T)))
        self.assertLess(gap(free),1e-12)
        self.assertGreater(gap(confined),.6)

    def test_monodromy_at_zero_matches_analytic_integral(self):
        for radius in (2.,4.):
            trace,det=continuum_discriminant(0.,radius)
            self.assertAlmostEqual(trace,2*np.cosh(2*np.arcsinh(radius)),delta=1e-9)
            self.assertAlmostEqual(det,1.,delta=1e-9)


if __name__=='__main__':unittest.main()
