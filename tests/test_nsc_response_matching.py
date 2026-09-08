"""Independent finite-frequency and metric-coordinate response controls."""
from dataclasses import replace
from pathlib import Path
import sys
import unittest
import numpy as np
from recursive_horizons.nsc_covariant_operator import smooth_metric
from recursive_horizons.nsc_response_matching import (
    canonical_data, static_proper_time, radius_path, retarded_susceptibility,
    euclidean_kernel, rational_frequency_bubble, abel_extrapolation)
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from check_nsc_response_matching import compare_record


class MatchingTests(unittest.TestCase):
    def test_matrix_frequency_loop_matches_spectral_retarded_function(self):
        for eta,n in ((0.,13),(.5,12)):
            m=smooth_metric(n,eta=eta);d=canonical_data(m)
            for nu in (0.,.6,1.4):
                self.assertAlmostEqual(rational_frequency_bubble(m,nu),
                                       retarded_susceptibility(d,1j*nu).real,delta=2e-11)

    def test_static_proper_time_bubble_has_independent_frequency_owner(self):
        m=smooth_metric(12);d=canonical_data(m)
        self.assertAlmostEqual(euclidean_kernel(m,0.,2.,frequency_points=64),
                               static_proper_time(d,2.)['inverse_radius_hessian'],delta=2e-10)

    def test_metric_coordinate_contact_matches_independent_energy_derivative(self):
        m=smooth_metric(12);d=canonical_data(m);s=static_proper_time(d,2.)
        step=.004
        for coordinate,key in (('inverse_radius','inverse_radius_hessian'),('log_radius','log_radius_hessian')):
            ep=static_proper_time(canonical_data(radius_path(m,step,coordinate)),2.)['energy']
            em=static_proper_time(canonical_data(radius_path(m,-step,coordinate)),2.)['energy']
            self.assertAlmostEqual((ep-2*s['energy']+em)/step**2,s[key],delta=2e-6)
        self.assertGreater(abs(s['log_radius_contact']),.01)

    def test_retarded_schwarz_and_positive_absorption(self):
        d=canonical_data(smooth_metric(12))
        z=.8+.3j
        self.assertAlmostEqual(abs(retarded_susceptibility(d,-z.conjugate())
                                    -retarded_susceptibility(d,z).conjugate()),0.,delta=1e-13)
        self.assertLess(retarded_susceptibility(d,z).imag,0.)

    def test_finite_proper_time_is_not_unregulated_vacuum_response(self):
        d=canonical_data(smooth_metric(12))
        p=static_proper_time(d,2.)['inverse_radius_hessian']
        k=retarded_susceptibility(d,0.).real
        self.assertGreater(abs(p-k),.1)
        self.assertLess(abs(static_proper_time(d,1e7)['inverse_radius_hessian']-k),1e-6)

    def test_Abel_extrapolation_is_quadratic_interpolation_at_zero(self):
        p=canonical_data(smooth_metric(33,eta=0.));a=canonical_data(smooth_metric(32))
        row=abel_extrapolation(p,a,.6,8.)
        coefficients=np.polyfit(1/np.array(row['cutoffs'])**2,row['differences'],2)
        self.assertAlmostEqual(coefficients[-1],row['extrapolated'],delta=1e-13)
        self.assertFalse(row['remainder_bound_proved'])

    def test_domains_and_zero_mode_prescriptions_are_not_silently_extended(self):
        m=smooth_metric(12)
        with self.assertRaises(ValueError):canonical_data(replace(m,lapse=np.full(12,1.1)))
        with self.assertRaises(ValueError):radius_path(m,.1,'unknown')
        with self.assertRaises(ValueError):retarded_susceptibility(canonical_data(m),-.2j)
        with self.assertRaises(ValueError):euclidean_kernel(m,-.1,2.)

    def test_all_field_comparator_rejects_scope_hash_and_numeric_changes(self):
        for left,right in [({'a':1.},{'a':1.001}),({'hash':'a'},{'hash':'b'}),
                           ({'self_sourced':False},{'self_sourced':True}),({'i':1},{'i':True})]:
            with self.assertRaises(RuntimeError):compare_record(left,right)


if __name__=='__main__':unittest.main()
