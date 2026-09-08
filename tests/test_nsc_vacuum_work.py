from dataclasses import replace
from pathlib import Path
import sys
import unittest

import numpy as np

from recursive_horizons.nsc_shape_response import smooth_metric, StaticAxialMetric
from recursive_horizons.nsc_energy_transfer import dirac_hamiltonian
from recursive_horizons.nsc_vacuum_work import radius_vertex, pulse_profile, leading_response, evolve_vacuum_pulse

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from check_nsc_vacuum_work import compare_record


class VacuumWorkTests(unittest.TestCase):
    def test_vertex_is_exact_finite_geometry_intervention(self):
        m=smooth_metric(2.,32)
        for eta in (0.,.5):
            h=dirac_hamiltonian(m,2,eta)
            v=radius_vertex(m,2,eta)
            for epsilon in (-.2,.3):
                changed=replace(m,sphere_radius=m.sphere_radius/(1+epsilon*pulse_profile(m)))
                np.testing.assert_allclose(dirac_hamiltonian(changed,2,eta),h+epsilon*v,atol=1e-12)

    def test_uniform_cylinder_matches_independent_momentum_formula(self):
        n=32
        m=StaticAxialMetric(4.,np.ones(n),np.ones(n),np.ones(n))
        phase=2*np.pi*(np.arange(n)+.5)/n
        p=2*np.sin(phase/2)/m.spacing
        mass=np.cos(phase/2)
        e=np.hypot(p,mass)
        probability=(mass*p/e)**2*(np.pi/2)*np.exp(-e*e)
        response=leading_response(dirac_hamiltonian(m),radius_vertex(m,profile=np.ones(n)))
        self.assertAlmostEqual(response['pair_number_coefficient'],sum(probability),delta=1e-12)
        self.assertAlmostEqual(response['energy_coefficient'],sum(2*e*probability),delta=1e-12)

    def test_clock_reparametrization_has_no_pair_production(self):
        h=dirac_hamiltonian(smooth_metric(2.,32))
        response=leading_response(h,h)
        self.assertLess(response['energy_coefficient'],1e-25)

    def test_vacuum_production_energy_includes_particle_and_hole(self):
        m=smooth_metric(2.,32)
        r=leading_response(dirac_hamiltonian(m),radius_vertex(m))
        self.assertGreater(r['energy_coefficient'],0.)
        self.assertAlmostEqual(r['energy_coefficient'],r['particle_energy_coefficient']+r['hole_energy_coefficient'],delta=1e-14)
        slow=leading_response(dirac_hamiltonian(m),radius_vertex(m),2.)
        self.assertLess(slow['energy_coefficient'],r['energy_coefficient']*1e-5)

    def test_finite_vacuum_evolution_accounts_for_work_and_fermion_pairs(self):
        r=evolve_vacuum_pulse(smooth_metric(2.,32),.02,rtol=2e-13)
        self.assertGreater(r['excitation_energy'],0.)
        self.assertLess(r['work_balance_residual'],1e-11)
        self.assertLess(r['regional_balance_residual'],1e-11)
        self.assertGreater(abs(r['integrated_regional_inflow']),1e-6)
        self.assertLess(abs(r['pair_number']-r['hole_number']),1e-12)
        self.assertLess(r['orthonormality_residual'],1e-13)

    def test_zero_drive_and_invalid_radius_are_separate_controls(self):
        m=smooth_metric(2.,32)
        r=evolve_vacuum_pulse(m,0.)
        self.assertLess(abs(r['excitation_energy']),1e-12)
        self.assertEqual(r['pair_number'],0.)
        with self.assertRaises(ValueError): evolve_vacuum_pulse(m,-1.)
        with self.assertRaises(ValueError): leading_response(dirac_hamiltonian(m),radius_vertex(m),0.)

    def test_reproduction_rejects_changed_small_signal_or_physical_claim(self):
        with self.assertRaises(RuntimeError): compare_record({'energy':3.4e-5},{'energy':3.5e-5})
        with self.assertRaises(RuntimeError): compare_record({'Q_derived':False},{'Q_derived':True})
        with self.assertRaises(RuntimeError): compare_record({'points':32},{'points':33})


if __name__=='__main__': unittest.main()
