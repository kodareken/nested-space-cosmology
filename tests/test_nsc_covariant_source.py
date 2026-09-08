from dataclasses import replace
from pathlib import Path
import sys
import unittest

import numpy as np

from recursive_horizons.nsc_covariant_operator import (
    CovariantStaticMetric,cutoff_response,cylinder_metric,cylinder_zero_winding_energy,
    directional_variation,euclidean_operator,frequency_density,heat_coefficients,
    scalar_derivative,smooth_metric,spatial_cutoff_proxy,ultraviolet_subtraction,
    ultrastatic_reference,ward_summary,
)
from recursive_horizons.nsc_covariant_measure import finite_cutoff_cylinder_stress
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from check_nsc_covariant_source import compare_record,schur_sign_controls,spectrum_controls


class CovariantSourceTests(unittest.TestCase):
    def test_fourier_domains_keep_symmetric_signed_momenta(self):
        for eta,n in ((0.,17),(.5,16)):
            m=smooth_metric(n,eta=eta)
            np.testing.assert_allclose(m.momenta,-m.momenta[::-1],atol=1e-14)
            np.testing.assert_allclose(m.momentum_matrix,m.momentum_matrix.conj().T,atol=1e-14)
            np.testing.assert_allclose(m.momentum_matrix.real,0.,atol=1e-13)
        with self.assertRaises(ValueError):smooth_metric(16,eta=0.)

    def test_scalar_metric_derivative_keeps_odd_highest_frequency(self):
        for n in (11,12):
            length=4.;x=np.arange(n)*length/n;mode=5
            values=np.cos(2*np.pi*mode*x/length)
            expected=-2*np.pi*mode/length*np.sin(2*np.pi*mode*x/length)
            np.testing.assert_allclose(scalar_derivative(values,length),expected,atol=1e-13)

    def test_physical_chirality_square_and_sign_reductions(self):
        for record in spectrum_controls():
            self.assertLess(record['full_paired_chirality_residual'],1e-11)
            self.assertLess(record['finite_commutator_square_residual'],1e-10)
            self.assertGreater(record['finite_metric_vs_Weyl_sandwich_operator_difference_norm'],1e-3)

    def test_covariant_frequency_integral_matches_ultrastatic_reference(self):
        m=smooth_metric(16)
        value=cutoff_response(m,angular_max=8)['energy']
        self.assertAlmostEqual(value,ultrastatic_reference(m,angular_max=8),delta=1e-9)

    def test_full_cylinder_stress_matches_independent_winding_response(self):
        for eta,n in ((0.,25),(.5,24)):
            m=cylinder_metric(n,eta=eta);r=cutoff_response(m,angular_max=16)
            ref=finite_cutoff_cylinder_stress(4.,cutoff=2.,eta=eta)
            self.assertAlmostEqual(r['energy']-cylinder_zero_winding_energy(angular_max=16),ref['energy'],delta=2e-10)
            np.testing.assert_allclose(r['radial_null'],ref['axial_null_stress_integral'],atol=2e-10)

    def test_covariant_lapse_scaling_distinguishes_spatial_cutoff(self):
        m=smooth_metric(16);base=cutoff_response(m,angular_max=8)['energy']
        changed=replace(m,lapse=1.2*m.lapse)
        correct=cutoff_response(changed,angular_max=8)['energy']
        proxy=spatial_cutoff_proxy(changed,angular_max=8)
        self.assertAlmostEqual(correct,1.2*base,delta=1e-9)
        self.assertGreater(abs(proxy-correct),1.)

    def test_metric_HF_gradient_matches_independent_energy_change(self):
        m=smooth_metric(16,general=True);r=cutoff_response(m,angular_max=8)
        shape=np.sin(2*np.pi*m.x/m.length)+.3
        for j,name in enumerate(('lapse','radial_scale','sphere_radius')):
            delta=getattr(m,name)*shape;step=1e-4
            plus=replace(m,**{name:getattr(m,name)+step*delta})
            minus=replace(m,**{name:getattr(m,name)-step*delta})
            fd=(cutoff_response(plus,angular_max=8,gradients=False)['energy']-cutoff_response(minus,angular_max=8,gradients=False)['energy'])/(2*step)
            self.assertAlmostEqual(fd,float(np.dot(r['metric_gradients'][j],delta)),delta=2e-7)

    def test_actual_varying_metric_stress_conserves_under_refinement(self):
        residual=[]
        for n in (16,24,32):
            m=smooth_metric(n,general=True);r=cutoff_response(m,angular_max=12)
            residual.append(ward_summary(m,r)['maximum_radial_conservation_residual'])
        self.assertLess(residual[1],residual[0]/20)
        self.assertLess(residual[2],residual[1]/20)

    def test_normalization_derivative_has_opposite_local_and_remainder_signs(self):
        m=smooth_metric(24,general=True);step=1e-4
        r=ultraviolet_subtraction(m,2.,1.)
        fd=(ultraviolet_subtraction(m,2.,np.exp(step))['energy']-ultraviolet_subtraction(m,2.,np.exp(-step))['energy'])/(2*step)
        self.assertAlmostEqual(fd,r['local_normalization_log_derivative'],delta=1e-10)
        self.assertEqual(r['local_normalization_log_derivative'],-r['remainder_normalization_log_derivative'])

    def test_product_heat_coefficients_match_exact_geometric_values(self):
        m=cylinder_metric(24);a=heat_coefficients(m)
        self.assertAlmostEqual(a['a0'],64*np.pi,delta=1e-12)
        self.assertAlmostEqual(a['a2'],-32*np.pi/3,delta=1e-12)
        self.assertAlmostEqual(a['a4'],-16*np.pi/15,delta=1e-12)

    def test_schur_reciprocal_and_frame_change_do_not_create_vacuum_drain(self):
        r=schur_sign_controls()
        self.assertGreater(r['self_energy_norm'],.01)
        self.assertGreater(r['current_operator_norm'],.01)
        self.assertLess(abs(r['vacuum_current']),1e-12)
        self.assertLess(r['consistent_child_frame_change_residual'],1e-12)
        self.assertTrue(np.all(np.array(r['real_energy_sign_controls'][0]['Schur_correction_eigenvalues'])>0))
        self.assertTrue(np.all(np.array(r['real_energy_sign_controls'][1]['Schur_correction_eigenvalues'])<0))

    def test_invalid_frequency_angular_sector_and_complex_metric_are_rejected(self):
        m=smooth_metric(16)
        with self.assertRaises(ValueError):euclidean_operator(m,np.nan)
        with self.assertRaises(ValueError):euclidean_operator(m,0.,0)
        with self.assertRaises(ValueError):replace(m,lapse=m.lapse+1j)
        with self.assertRaises(ValueError):cutoff_response(m,angular_max=True)

    def test_reproduction_rejects_changed_sign_claim_and_dimensions(self):
        with self.assertRaises(RuntimeError):compare_record({'source':-.009},{'source':.009})
        with self.assertRaises(RuntimeError):compare_record({'physical_closure':False},{'physical_closure':True})
        with self.assertRaises(RuntimeError):compare_record({'points':24},{'points':25})


if __name__=='__main__':unittest.main()
