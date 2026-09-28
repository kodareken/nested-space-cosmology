import unittest
import numpy as np
from recursive_horizons.nsc_lorentzian import sbp_operator,evolve,geometry,horizon,characteristic


class LorentzianTests(unittest.TestCase):
    def test_discrete_green_identity_includes_only_physical_inflow(self):
        op=sbp_operator(101)
        self.assertLess(op['boundary_form_residual'],1e-12)
        self.assertEqual([(ch,node) for ch,node,_ in op['incoming']],[(1,-1)])
        self.assertTrue(np.all(op['speeds'][:,0]<0))

    def test_outward_local_light_direction_moves_toward_child_at_throat(self):
        self.assertLess(1-float(geometry(0.)[0]),0)
        self.assertLess(characteristic(.4,.4,0),.4)
        self.assertAlmostEqual(float(geometry(horizon())[2]),0.,places=11)

    def test_number_plus_outflow_ledger(self):
        row=evolve(401,final_time=.8)
        self.assertLess(abs(row['probability_balance_residual']),1e-8)
        self.assertLess(abs(row['child_balance_residual']),1e-8)
        self.assertGreater(row['child_probability'],.98)


if __name__=='__main__':unittest.main()
