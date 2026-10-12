"""New normalization and regional-state checks, not a cosmological campaign."""
from pathlib import Path
import sys
import unittest
import sympy as sp

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from recursive_horizons.nsc_nested_qualities import finite_window,exact_controls


class NestedQualitiesTests(unittest.TestCase):
    def test_exact_new_controls_and_mutations(self):
        r=exact_controls()
        self.assertTrue(all(r['checks'].values()))
        self.assertTrue(r['checks']['omitted_omega_factor_rejected'])
        self.assertTrue(r['checks']['same_law_nonstationary_state'])

    def test_shift_holds_at_different_depths_and_scales(self):
        h=sp.Matrix([[2,sp.I],[-sp.I,3]]);b=sp.Matrix([[1,2],[sp.I,0]])
        for depth in [1,2,4]:
            for omega in [sp.Rational(2,3),sp.Integer(2)]:
                self.assertEqual(finite_window(h,b,omega,-1,depth)/omega,
                                 finite_window(h,b,omega,-2,depth))

    def test_adding_a_region_is_not_a_shift_of_the_same_window(self):
        h=sp.Matrix([[1]]);b=sp.Matrix([[1]]);z=sp.I
        responses=[]
        for depth in [2,3]:
            j=finite_window(h,b,2,0,depth)
            responses.append((z*sp.eye(depth)-j).inv()[0,0])
        self.assertNotEqual(sp.simplify(responses[0]-responses[1]),0)

    def test_invalid_operators_and_scale_rejected(self):
        for h,b,o,n in [(sp.Matrix([[sp.I]]),sp.ones(1),2,2),
                        (sp.eye(2),sp.ones(1),2,2),(sp.eye(1),sp.ones(1),0,2),
                        (sp.eye(1),sp.ones(1),2,0)]:
            with self.assertRaises(ValueError):finite_window(h,b,o,0,n)

if __name__=='__main__':unittest.main()
