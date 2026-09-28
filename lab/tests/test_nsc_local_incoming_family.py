import unittest

import mpmath as mp
import numpy as np

from recursive_horizons.nsc_local_incoming_family import (
    LocalAxialFunction, LocalIncomingFamily, plateau_derivatives,
)
from recursive_horizons.nsc_evolved_incoming_constraints import compatible_history_slots


class LocalIncomingFamilyTests(unittest.TestCase):
    def test_window_derivatives_against_independent_high_precision_differentiation(self):
        with mp.workdps(60):
            for point in (-.059,-.043,-.031,.031,.043,.059):
                def reference(z):
                    u=(abs(z)-mp.mpf(.03))/(mp.mpf(.06)-mp.mpf(.03))
                    return 1/(1+mp.exp(1/(1-u)-1/u))
                actual=plateau_derivatives(point,0.)
                for n in range(5):
                    expected=float(mp.diff(reference,mp.mpf(point),n))
                    self.assertLessEqual(abs(actual[n]-expected),2e-11*max(1.,abs(expected)))

    def test_exact_plateau_support_and_polynomial_jets(self):
        f=LocalAxialFunction((.2,.1,-.3,.4),0.)
        poly=np.polynomial.Chebyshev(f.coefficients,domain=[-.06,.06])
        for z in (-.03,0.,.03):
            for n in range(5):
                self.assertAlmostEqual(f(z,n),poly.deriv(n)(z),delta=2e-10)
        for z in (-.07,-.06,.06,.07):
            self.assertEqual([f(z,n) for n in range(5)],[0.]*5)

    def test_metric_slots_are_derivatives_of_same_two_functions(self):
        coeff=np.zeros((2,8));coeff[0,:4]=[.03,-.02,.01,.005];coeff[1,:3]=[2.,-1.,.5]
        family=LocalIncomingFamily(coeff)
        metric=family.metric();z=family.collocation_nodes()
        slots,tangent=compatible_history_slots(metric,z,16)
        w,U=family.functions
        expected=np.array([[w(v,0),w(v,1),w(v,2),w(v,3),U(v,0),U(v,1)] for v in z])
        np.testing.assert_allclose(slots,expected,atol=3e-11,rtol=2e-15)
        direction=5;step=1e-5
        altered=coeff.copy();altered.ravel()[direction]+=step
        plus=LocalIncomingFamily(altered).metric()
        value,_=compatible_history_slots(plus,z,16)
        np.testing.assert_allclose((value-slots)/step,tangent[direction],atol=2e-7,rtol=2e-10)
        self.assertGreater(family.radius_lower_bound(),0.)
        self.assertGreater(family.interval[1]-family.interval[0],0.)

    def test_unsafe_radius_and_invalid_inputs_are_visible(self):
        coeff=np.zeros((2,8));coeff[0,0]=100.
        self.assertLess(LocalIncomingFamily(coeff).radius_lower_bound(),0.)
        for value in (np.zeros((2,7)),np.zeros((1,8)),np.full((2,8),np.nan)):
            with self.assertRaises(ValueError):LocalIncomingFamily(value)
        with self.assertRaises(ValueError):LocalIncomingFamily(np.zeros((2,8)),normal_inner=.04)


if __name__=='__main__':unittest.main()
