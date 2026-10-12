"""Same-functional Weyl identities, metric variation and discretization controls."""
from dataclasses import replace
import unittest
import numpy as np
from scipy.linalg import eigh
from recursive_horizons.nsc_covariant_operator import smooth_metric,euclidean_operator,SIGMA2
from recursive_horizons.nsc_influence import canonical_hamiltonian
from recursive_horizons.nsc_weyl_cocycle import (
    conformal_metric,ratio_hamiltonian,ratio_operator,cutoff_functional,fiber,
    heat_path_cocycle,operator_variations,compensated_lift,continuum_identities)
from recursive_horizons.nsc_regulated import RegulatedOperator,OperatorConventions


class WeylCompletionTests(unittest.TestCase):
    def setUp(self):
        self.m=smooth_metric(12,general=True)
        self.s=.04+.07*np.cos(2*np.pi*self.m.x/self.m.length)

    def test_continuum_ratios_and_connection_are_invariant(self):
        self.assertTrue(all(value=='0' for value in continuum_identities().values()))

    def test_finite_ratio_operator_preserves_the_actual_sandwich(self):
        changed=conformal_metric(self.m,self.s)
        h=ratio_hamiltonian(self.m);hc=ratio_hamiltonian(changed)
        np.testing.assert_allclose(h,hc,atol=1e-13,rtol=1e-13)
        w=np.exp(-np.r_[self.s,self.s]/2)
        d=ratio_operator(self.m,.7)
        np.testing.assert_allclose(ratio_operator(changed,.7),w[:,None]*d*w[None,:],atol=1e-13)
        np.testing.assert_allclose(d,d.conj().T,atol=1e-13)

    def test_ratio_and_prior_discretization_agree_on_ultrastatic_cell(self):
        m=smooth_metric(16)
        np.testing.assert_allclose(ratio_operator(m,.7),euclidean_operator(m,.7),atol=1e-14)
        np.testing.assert_allclose(ratio_hamiltonian(m),canonical_hamiltonian(m),atol=1e-14)

    def test_paired_physical_chirality_is_preserved(self):
        m=conformal_metric(self.m,self.s);n=2*m.points;zero=np.zeros((n,n),complex)
        h=np.block([[ratio_hamiltonian(m,1),zero],[zero,ratio_hamiltonian(m,-1)]])
        d=np.block([[ratio_operator(m,.7,1),zero],[zero,ratio_operator(m,.7,-1)]])
        off=-np.kron(SIGMA2,np.eye(m.points));gamma5=np.block([[zero,off],[off,zero]])
        np.testing.assert_allclose(h@gamma5,gamma5@h,atol=1e-12)
        np.testing.assert_allclose(d@gamma5,-gamma5@d,atol=1e-12)

    def test_finite_metric_gradients_match_independent_energy_changes(self):
        row=fiber(self.m,.7,2,2.)
        for i,name in enumerate(('lapse','radial_scale','sphere_radius')):
            direction=.2+.1*np.cos(2*np.pi*self.m.x/self.m.length);step=1e-5
            original=getattr(self.m,name)
            plus=replace(self.m,**{name:original*np.exp(step*direction)})
            minus=replace(self.m,**{name:original*np.exp(-step*direction)})
            ep=fiber(plus,.7,2,2.,False)['energy'];em=fiber(minus,.7,2,2.,False)['energy']
            self.assertAlmostEqual((ep-em)/(2*step),np.dot(row['log_gradients'][i],direction),delta=2e-9)

    def test_first_and_second_operator_variations_are_independent_of_difference_step(self):
        direction=np.array([self.s,.2*self.s,-.4*self.s])
        d,d1,d2=operator_variations(self.m,.7,1,direction);step=1e-3
        def changed(t):return replace(self.m,**{name:getattr(self.m,name)*np.exp(t*f)
            for name,f in zip(('lapse','radial_scale','sphere_radius'),direction)})
        ep=ratio_operator(changed(step),.7);em=ratio_operator(changed(-step),.7)
        np.testing.assert_allclose((ep-em)/(2*step),d1,atol=3e-9,rtol=1e-7)
        np.testing.assert_allclose((ep-2*d+em)/step**2,d2,atol=2e-8,rtol=1e-6)

    def test_local_Weyl_source_is_heat_insertion_and_has_a_nonzero_hessian(self):
        row=cutoff_functional(self.m,2.,4,48)
        self.assertLess(row['local_Weyl_Ward_residual'],1e-12)
        direction=np.array([self.s]*3)
        d,d1,d2=operator_variations(self.m,.7,1,direction)
        operator=RegulatedOperator(d,OperatorConventions(cutoff=2.))
        gradient,hessian=operator.variation(d1,d2)
        one=fiber(self.m,.7,1,2.)
        self.assertAlmostEqual(gradient,np.dot(self.s,one['heat_diagonal']),delta=1e-12)
        self.assertGreater(abs(hessian),1e-6)

    def test_integrated_heat_path_matches_difference_of_the_same_action(self):
        base=cutoff_functional(self.m,2.,4,48,False)['energy']
        end=cutoff_functional(conformal_metric(self.m,self.s),2.,4,48,False)['energy']
        path=heat_path_cocycle(self.m,self.s,2.,4,48,8)
        self.assertAlmostEqual(end-base,path,delta=3e-10)

    def test_conditional_lift_has_gauge_redundancy_without_stationarity(self):
        row=compensated_lift(self.m,self.s,2.,4,48)
        shift=.03*np.sin(2*np.pi*self.m.x/self.m.length)
        other=compensated_lift(conformal_metric(self.m,shift),self.s-shift,2.,4,48)
        self.assertAlmostEqual(row['lifted_energy'],other['lifted_energy'],delta=1e-12)
        self.assertLess(row['Ward_residual'],1e-12)
        self.assertGreater(np.linalg.norm(row['phi_gradient']),1.)
        self.assertAlmostEqual(row['base_energy']+row['compensating_action'],row['lifted_energy'],delta=1e-12)

    def test_positive_Weyl_congruence_preserves_signature(self):
        before=eigh(ratio_operator(self.m,.7),eigvals_only=True)
        after=eigh(ratio_operator(conformal_metric(self.m,3*self.s),.7),eigvals_only=True)
        self.assertEqual(sum(before<0),sum(after<0))
        self.assertGreater(min(abs(after)),.1)


if __name__=='__main__':unittest.main()
