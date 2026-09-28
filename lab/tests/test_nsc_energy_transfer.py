import unittest

import numpy as np
from scipy.linalg import expm

from recursive_horizons.nsc_energy_transfer import (
    dirac_hamiltonian, positions, regional_operators, split_cut_currents,
    vacuum, expectation, positive_packet, evolve_packet, integrated_current,
    schur_state_data,
)
from recursive_horizons.nsc_shape_response import smooth_metric


class EnergyTransferTests(unittest.TestCase):
    def setUp(self):
        self.metric = smooth_metric(2., 32)
        self.h = dirac_hamiltonian(self.metric)
        self.mask = positions(self.metric) >= 0
        self.e, self.v, self.c = vacuum(self.h)
        self.ops = regional_operators(self.h, self.mask)

    def test_half_link_and_bare_energy_accounting_agree(self):
        o = self.ops
        np.testing.assert_allclose(o['energy_a']+o['energy_b'], self.h)
        np.testing.assert_allclose(o['bare_a']+o['bare_b']+o['link'], self.h)
        np.testing.assert_allclose(o['energy_a'], o['bare_a']+o['link']/2)
        p = np.diag(self.mask.astype(float))
        np.testing.assert_allclose(o['current_a'], .5j*(self.h@self.h@p-p@self.h@self.h), atol=1e-12)

    def test_ground_state_stationary_after_complex_basis_change(self):
        phase = np.exp(1j*np.arange(len(self.h))*.713)
        changed = phase[:, None]*self.h*phase.conj()[None, :]
        _, _, c = vacuum(changed)
        current = regional_operators(changed, self.mask)['current_a']
        self.assertGreater(np.max(np.abs(current.real)), 1.)
        self.assertLess(abs(expectation(c, current)), 1e-10)

    def test_occupation_changes_flow_without_changing_retarded_response(self):
        psi = positive_packet(self.metric, self.e, self.v)
        psi = evolve_packet(self.e, self.v, psi, [.6])[:, 0]
        excited = self.c+np.outer(psi, psi.conj())
        a = schur_state_data(self.h, self.mask, .3+.4j, self.c)
        b = schur_state_data(self.h, self.mask, .3+.4j, excited)
        self.assertEqual(a['self_energy_norm'], b['self_energy_norm'])
        self.assertLess(a['direct_reduced_error'], 1e-12)
        self.assertGreater(abs(expectation(excited, self.ops['current_a'])), 1.)
        self.assertLess(abs(expectation(self.c, self.ops['current_a'])), 1e-12)

    def test_excitation_satisfies_pauli_constraint(self):
        psi = positive_packet(self.metric, self.e, self.v)
        c = self.c+np.outer(psi,psi.conj())
        np.testing.assert_allclose(c@c, c, atol=1e-12)
        self.assertLess(np.linalg.norm(self.c@psi), 1e-12)

    def test_energy_rate_matches_independent_time_derivative(self):
        psi = positive_packet(self.metric, self.e, self.v)
        psi = expm(-.51j*self.h)@psi
        step = 1e-5
        energy = lambda t: np.vdot(expm(-1j*t*self.h)@psi,
                                   self.ops['energy_a']@(expm(-1j*t*self.h)@psi)).real
        fd = (energy(step)-energy(-step))/(2*step)
        exact = np.vdot(psi, self.ops['current_a']@psi).real
        self.assertAlmostEqual(fd, exact, delta=1e-7)

    def test_cut_integrals_sum_to_regional_energy_change(self):
        psi = positive_packet(self.metric, self.e, self.v)
        cuts = split_cut_currents(self.ops['current_a'], positions(self.metric), self.metric.length)
        self.assertGreater(np.linalg.norm(cuts['seam']), 0.)
        np.testing.assert_allclose(cuts['seam']+cuts['throat'], self.ops['current_a'])
        final = expm(-1.2j*self.h)@psi
        delta = np.vdot(final,self.ops['energy_a']@final).real-np.vdot(psi,self.ops['energy_a']@psi).real
        total = sum(integrated_current(self.e,self.v,psi,op,1.2) for op in cuts.values())
        self.assertAlmostEqual(delta,total,delta=1e-11)

    def test_unit_rescaling_distinguishes_energy_from_transfer_rate(self):
        scale = 2.3
        current = regional_operators(self.h/scale,self.mask)['current_a']
        np.testing.assert_allclose(current*scale**2,self.ops['current_a'],atol=1e-11)

    def test_time_dependent_energy_work_term_is_required(self):
        # A(t)=H0+tV, evolving a fixed initial vector for an infinitesimal step.
        # The derivative at zero separates state evolution and explicit work.
        psi = positive_packet(self.metric,self.e,self.v)
        perturbation = np.diag(np.sin(positions(self.metric)))
        work_a = regional_operators(perturbation,self.mask)['energy_a']
        dt = 1e-6
        values = []
        for t in (-dt,dt):
            state = expm(-1j*self.h*t)@psi
            operator = self.ops['energy_a']+t*work_a
            values.append(np.vdot(state,operator@state).real)
        exact = np.vdot(psi,(self.ops['current_a']+work_a)@psi).real
        self.assertAlmostEqual((values[1]-values[0])/(2*dt),exact,delta=1e-8)

    def test_initial_child_amplitude_is_not_removed_by_schur_elimination(self):
        # Laplace resolvent identity includes the source from the eliminated side.
        z = .3+.4j
        rng = np.random.default_rng(73)
        initial = rng.normal(size=len(self.h))+1j*rng.normal(size=len(self.h))
        p = self.mask
        aa,bb,b = self.h[np.ix_(p,p)],self.h[np.ix_(~p,~p)],self.h[np.ix_(p,~p)]
        rb = np.linalg.inv(z*np.eye(len(bb))-bb)
        effective = z*np.eye(len(aa))-aa-b@rb@b.T.conj()
        correct = np.linalg.solve(effective,initial[p]+b@rb@initial[~p])
        direct = np.linalg.solve(z*np.eye(len(self.h))-self.h,initial)[p]
        np.testing.assert_allclose(correct,direct,atol=1e-12)
        omitted = np.linalg.solve(effective,initial[p])
        self.assertGreater(np.linalg.norm(omitted-direct),.1)

    def test_exact_two_level_control_has_transfer_and_recurrence(self):
        # H=e I+b sigma1, occupied left site. E_A=e cos²(bt).
        energy,coupling = 2.,.7
        h = np.array([[energy,coupling],[coupling,energy]])
        o = regional_operators(h,[1,0])
        for time in (0.,.3,1.,np.pi/coupling):
            psi = expm(-1j*time*h)@np.array([1.,0.])
            self.assertAlmostEqual(np.vdot(psi,o['energy_a']@psi).real,
                                   energy*np.cos(coupling*time)**2,delta=1e-12)
            self.assertAlmostEqual(np.vdot(psi,o['current_a']@psi).real,
                                   -energy*coupling*np.sin(2*coupling*time),delta=1e-12)


if __name__ == '__main__':
    unittest.main()
