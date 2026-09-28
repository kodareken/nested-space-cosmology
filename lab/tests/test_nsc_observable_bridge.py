from pathlib import Path
import sys
import unittest

import numpy as np
import sympy as sp

from recursive_horizons.nsc_spinor_bridge import pauli, weyl_matrices, exact_checks, chiral_evolution
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from check_nsc_observable_bridge import conservation_checks


class ObservableBridgeTests(unittest.TestCase):
    def test_exact_clifford_charge_and_radial_identities(self):
        result=exact_checks()
        self.assertEqual(set(result['identities'].values()),{'0'})
        self.assertEqual(set(result['joint_projector_ranks'].values()),{2})
        self.assertFalse(result['candidate_restriction']['selected_by_actual_throat_or_action'])

    def test_restricted_embedding_matches_direct_evolution(self):
        w=weyl_matrices();beta=np.array(w['gamma'][0],complex)
        kinetic=sum(p*np.array(a,complex) for p,a in zip((.3,-.7,1.1),w['alpha']))
        h=np.kron(np.eye(2),kinetic)+.6*np.kron(np.array(pauli()[0],complex),beta)
        embedding=np.eye(8)[:,[0,1,6,7]]
        initial=np.array([1,2j,-.3,.7j]);initial/=np.linalg.norm(initial)
        from scipy.linalg import expm
        full=expm(-.8j*h)@embedding@initial
        reduced=embedding@expm(-.8j*(kinetic+.6*beta))@initial
        np.testing.assert_allclose(full,reduced,atol=1e-12)

    def test_identity_link_closes_gap_at_finite_momentum(self):
        w=weyl_matrices();alpha=np.array(w['alpha'][2],complex)
        h=np.kron(np.eye(2),alpha)+np.kron(np.array(pauli()[0],complex),np.eye(4))
        self.assertLess(np.min(abs(np.linalg.eigvalsh(h))),1e-12)

    def test_massive_energy_eigenstate_has_constant_chirality(self):
        result=chiral_evolution()
        self.assertLess(np.ptp([r['positive_energy_right_probability'] for r in result['controls']]),1e-12)
        self.assertGreater(result['pure_left_negative_energy_weight'],0.)
        self.assertLess(result['pure_left_negative_energy_weight'],1.)

    def test_massless_positive_energy_left_state_has_no_flip(self):
        result=chiral_evolution(momentum=1.,mass=0.,helicity=-1)
        self.assertEqual(result['pure_left_negative_energy_weight'],0.)
        self.assertTrue(all(r['pure_left_to_right']==0. for r in result['controls']))

    def test_charge_conjugation_conjugates_U1_phase(self):
        w=weyl_matrices();b=np.array(w['charge_conjugation_matrix'],complex)
        psi=np.array([1.,.3j,-.2,.5j])
        phase=np.exp(.37j)
        np.testing.assert_allclose(b@(phase*psi).conj(),phase.conjugate()*(b@psi.conj()),atol=1e-12)

    def test_source_partition_formula_matches_numeric_derivative(self):
        # Nonzero pressures and both internal/external sources prevent an
        # accidentally specialized dust/closed formula passing this control.
        b,d,h,wb,wd,q,jb,jd=2.,5.,.7,.1,-.8,.3,.2,-.1
        bp=-3*(1+wb)*b+(q+jb)/h
        dp=-3*(1+wd)*d+(-q+jd)/h
        step=1e-6
        frac=lambda t:(d+t*dp)/(b+t*bp+d+t*dp)
        fd=(frac(step)-frac(-step))/(2*step)
        f=d/(b+d)
        expected=3*f*(1-f)*(wb-wd)-q/(h*(b+d))+jd/(h*(b+d))-f*(jb+jd)/(h*(b+d))
        self.assertAlmostEqual(fd,expected,delta=1e-9)

    def test_zero_exchange_does_not_force_deceleration(self):
        result=conservation_checks()
        control=result['zero_Q_counterexamples'][0]
        self.assertEqual(control['internal_Q'],0.)
        self.assertLess(control['rows'][-1]['q_dec'],-.99)
        self.assertGreater(control['rows'][-1]['dark_fraction'],.999)
        other=result['zero_Q_counterexamples'][1]
        self.assertNotEqual(control['rows'][2]['H_over_H0'],other['rows'][2]['H_over_H0'])


if __name__=='__main__':unittest.main()
