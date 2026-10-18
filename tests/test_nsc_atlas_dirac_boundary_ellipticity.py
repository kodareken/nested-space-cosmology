"""Real-frequency boundary and full-angular-direction ellipticity controls."""
import numpy as np
import pytest

from recursive_horizons.nsc_atlas_dirac_evolution import real_frequency_boundary_symbol


def test_actual_real_boundary_gram_and_glancing_are_separate():
    beta=np.array([[0.,1.],[1.,0.]])
    for W in (-2.3,2.3):
        for fraction in (-.9,0.,.9):
            tau=fraction*abs(W)
            row=real_frequency_boundary_symbol(tau,W)
            M=row['matrix']
            expected=2*(W*W*np.eye(2)-tau*W*beta)
            assert M.conj().T@M==pytest.approx(expected,abs=2e-14)
            assert np.linalg.svd(M,compute_uv=False)[-1]**2==pytest.approx(row['squared_lower_bound'],abs=2e-14)
    for tau,W in ((1.,1.),(-1.,1.),(2.,1.),(0.,0.)):
        row=real_frequency_boundary_symbol(tau,W)
        assert row['status']=='characteristic_or_glancing'
        assert row['matrix'] is None


def test_full_tangent_symbol_has_the_same_uniform_complementing_bound():
    # Full product spinor, before choosing an angular eigenvector. This
    # verifies direction-uniformity rather than assuming a scalar k covers it.
    s1=np.array([[0.,1.],[1.,0.]])
    s2=np.array([[0.,-1j],[1j,0.]])
    normal=np.kron(np.diag([1.,-1.]),np.kron(s2,np.eye(2)))
    J=np.vstack((np.eye(4),np.eye(4)))/np.sqrt(2)
    for eta in ((1.3,.7),(-1.3,.7),(.2,-2.1)):
        K=eta[0]*s1+eta[1]*s2
        tangent=np.kron(np.eye(2),np.kron(s1,K))
        angular_norm=np.hypot(*eta)
        tau=.6*angular_norm;decay=np.sqrt(angular_norm**2-tau**2)
        M=(1j*decay*normal+tangent-tau*np.eye(8))@J
        lower=2*angular_norm*(angular_norm-abs(tau))
        assert np.linalg.svd(M,compute_uv=False)[-1]**2==pytest.approx(lower,abs=2e-14)
