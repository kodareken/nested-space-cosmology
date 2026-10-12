"""Cross-check independent representations of the same throat equations."""
import unittest
import numpy as np
from recursive_horizons.nsc_boundary import continuum_half_map, staggered_forward
from recursive_horizons.nsc_geometric_chain import geometric_motif


class ThroatIntegrationTests(unittest.TestCase):
    def test_nonreal_first_order_maps_approach_exact_squared_threshold_jump(self):
        radius=4.
        t=np.arcsinh(radius)
        ip=np.exp(3*t)/6+np.exp(t)/2-2/3
        im=2/3-np.exp(-t)/2-np.exp(-3*t)/6
        expected=1/ip+1/im
        errors=[]
        for energy in (1e-3j,1e-4j):
            parent=continuum_half_map(radius,energy,side='parent')['weyl_m_throat']
            child=continuum_half_map(radius,energy,side='child')['weyl_m_throat']
            jump=energy*(child-parent)
            errors.append(abs(jump-expected))
        self.assertLess(errors[1],1e-6)
        self.assertLess(errors[1],errors[0]/50)

    def test_geometric_chain_links_are_the_boundary_discretization_entries(self):
        radius=4.; count=32
        forward=staggered_forward(np.linspace(-radius,radius,count+1),1.)
        motif=geometric_motif(radius,count)
        np.testing.assert_allclose(motif['H'][count:,:count],forward[:,:count],atol=1e-14)
        self.assertAlmostEqual(motif['B'][-1,0],forward[-1,-1],places=13)
        self.assertEqual(np.count_nonzero(motif['B']),1)


if __name__=='__main__':unittest.main()
