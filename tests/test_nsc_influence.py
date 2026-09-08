from pathlib import Path
import sys
import unittest

import numpy as np
from scipy.linalg import expm

from recursive_horizons.nsc_covariant_operator import smooth_metric
from recursive_horizons.nsc_influence import (
    canonical_hamiltonian,geometric_vertex,ground_covariance,thermal_covariance,
    influence,connected,gaussian_smeared_vertex,noise_kernel,response_spectrum,
    mixed_response_from_influence,retarded_response,pulse_unitary,
    with_real_branch_action,fock_annihilators,second_quantize,
)
from recursive_horizons.nsc_vacuum_work import leading_response
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from check_nsc_influence import fock_controls,compare_record


class InfluenceTests(unittest.TestCase):
    def setUp(self):
        m=smooth_metric(12);self.h=canonical_hamiltonian(m);self.v=geometric_vertex(m)
        _,_,self.c=ground_covariance(self.h)

    def test_small_Fock_oracle_matches_noncommuting_mixed_and_pure_states(self):
        r=fock_controls();self.assertLess(r['mixed_state_error'],1e-12)
        self.assertLess(abs(r['pure_orbital_determinant']-r['pure_orbital_Fock_trace']),1e-12)

    def test_Fock_operators_have_the_required_CAR(self):
        aa=fock_annihilators(3);identity=np.eye(8)
        for i in range(3):
            for j in range(3):
                np.testing.assert_allclose(aa[i]@aa[j].conj().T+aa[j].conj().T@aa[i],int(i==j)*identity)

    def test_equal_histories_and_swapped_histories(self):
        plus=expm(-.4j*(self.h+.2*self.v));minus=expm(-.4j*(self.h-.1*self.v))
        self.assertLess(abs(influence(self.c,plus,plus)['amplitude']-1),1e-12)
        a=influence(self.c,plus,minus)['amplitude'];b=influence(self.c,minus,plus)['amplitude']
        self.assertLess(abs(a.conjugate()-b),1e-12);self.assertLessEqual(abs(a),1+1e-12)

    def test_invalid_covariance_and_nonunitary_history_are_rejected(self):
        i=np.eye(len(self.h))
        with self.assertRaises(ValueError):influence(1.1*i,i,i)
        with self.assertRaises(ValueError):influence(self.c,2*i,i)

    def test_zero_overlap_is_valid_but_its_log_is_not_a_branch(self):
        r=influence(np.array([[.5]]),np.array([[-1.]]),np.array([[1.]]))
        self.assertEqual(r['amplitude'],0.);self.assertIsNone(r['principal_action'])

    def test_small_nonzero_overlap_retains_log_when_amplitude_underflows(self):
        i=np.eye(64);phase=np.exp(1j*(np.pi-1e-8))
        r=influence(.5*i,phase*i,i)
        self.assertTrue(r['amplitude_underflow'])
        self.assertIsNotNone(r['principal_action'])
        self.assertTrue(np.isfinite(r['log_modulus']))

    def test_smeared_variance_equals_independent_pair_creation_weight(self):
        a=gaussian_smeared_vertex(self.h,self.v,.5)
        self.assertAlmostEqual(connected(self.c,a,a).real,leading_response(self.h,self.v)['pair_number_coefficient'],delta=1e-12)

    def test_geometric_noise_matrix_is_positive(self):
        matrix=noise_kernel(self.h,self.c,self.v,[0.,.2,.6,1.])
        self.assertGreater(np.min(np.linalg.eigvalsh(matrix)),-1e-12)

    def test_contour_ordering_gives_retarded_response_and_future_cancellation(self):
        for source in (.2,1.1):
            observed=mixed_response_from_influence(self.h,self.c,self.v,self.v,.7,source)
            expected=-retarded_response(self.h,self.c,self.v,self.v,.7,source)
            self.assertAlmostEqual(observed,expected,delta=2e-7)

    def test_thermal_positive_frequency_fluctuation_dissipation(self):
        c=thermal_covariance(self.h,1.7)
        for row in response_spectrum(self.h,c,self.v):
            self.assertAlmostEqual(row['noise_plus'],row['minus_im_chi_plus']/np.tanh(.85*row['frequency']),delta=1e-12)

    def test_temperature_depends_on_energy_differences_and_has_ground_limit(self):
        a=thermal_covariance(self.h,1.7)
        b=thermal_covariance(self.h+2*np.eye(len(self.h)),1.7,chemical_potential=2.)
        np.testing.assert_allclose(a,b,atol=1e-12)
        np.testing.assert_allclose(thermal_covariance(self.h,40.),self.c,atol=1e-12)

    def test_time_ordered_geometric_pulse_preserves_fermionic_occupation(self):
        unitary=pulse_unitary(self.h,self.v,.02,steps=128)
        c=unitary@self.c@unitary.conj().T
        np.testing.assert_allclose(c@c,c,atol=1e-11)
        self.assertAlmostEqual(np.trace(c).real,np.trace(self.c).real,delta=1e-10)

    def test_local_real_phase_does_not_fix_absolute_source(self):
        z=.8*np.exp(.2j)
        shifted=with_real_branch_action(z,.7,.1)
        self.assertAlmostEqual(abs(shifted),abs(z),delta=1e-14)
        self.assertEqual(with_real_branch_action(1.,.7,.7),1.)
        self.assertGreater(abs(shifted-z),.1)

    def test_corrupted_scope_or_noise_sign_is_rejected(self):
        with self.assertRaises(RuntimeError):compare_record({'gravity_closed':False},{'gravity_closed':True})
        with self.assertRaises(RuntimeError):compare_record({'noise':.035},{'noise':-.035})


if __name__=='__main__':unittest.main()
