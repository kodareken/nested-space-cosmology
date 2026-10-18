"""Independent bounded-band endpoint and exact continuum resolvent controls."""
import numpy as np
import pytest
from scipy.linalg import solve_sylvester

from recursive_horizons.nsc_atlas_dirac_bounds import (
    balanced_band_bound,double_band_bounds,resolvent_separation_bound,
    frozen_projector_difference_bound)


def power(A,n):
    e,U=np.linalg.eigh(A)
    return (U*e**n)@U.conj().T


def test_balanced_endpoints_and_full_double_graph_bounds_against_independent_solves():
    rng=np.random.default_rng(704)
    for trial in range(5):
        v=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3))
        w=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4))
        A=np.eye(3)+v@v.conj().T;B=2*np.eye(4)+w@w.conj().T
        Y=rng.normal(size=(3,4))+1j*rng.normal(size=(3,4))
        norm=np.linalg.norm(Y,2);X=Y.copy()
        for n in range(1,6):
            X=solve_sylvester(A,B,X)
            balanced=power(A,n/2)@X@power(B,n/2)
            assert np.linalg.norm(balanced,2)<=balanced_band_bound(norm,n)+2e-12
            if n==2:
                bounds=double_band_bounds(norm,np.linalg.eigvalsh(A)[0],np.linalg.eigvalsh(B)[0])
                assert np.linalg.norm(A@X,2)<=bounds['left_energy']+2e-13
                assert np.linalg.norm(X@B,2)<=bounds['right_energy']+2e-13
                assert np.linalg.norm(X,2)<=bounds['operator_norm']+2e-13
    # Equal scalar bands attain the constants; the balanced endpoint
    # cannot be dismissed by a divergent integral of separate sup norms.
    A=3*np.eye(3);Y=np.eye(3);X=Y.copy()
    for n in range(1,6):
        X=solve_sylvester(A,A,X)
        assert np.linalg.norm(power(A,n/2)@X@power(A,n/2),2)==pytest.approx(2**(-n),abs=2e-14)
    second=solve_sylvester(A,A,solve_sylvester(A,A,Y))
    assert np.linalg.norm(A@second,2)==pytest.approx(double_band_bounds(1.,3.,3.)['left_energy'])


def test_separated_bound_against_exact_constant_mass_whole_line_resolvent():
    alpha=np.array([[0.,-1j],[1j,0.]])
    beta=np.array([[0.,1.],[1.,0.]])
    for mass,zeta,d in ((1.7,0.,.3),(1.7,2.1,.4),(3.2,-.8,.7)):
        rho=np.hypot(mass,zeta)
        # Exact Green kernel between x<=-d and y>=0 factorizes into
        # e^(rho*x)*e^(-rho*y)*M/(2rho). Its half-line norm is explicit.
        M=-1j*rho*alpha+mass*beta+1j*zeta*np.eye(2)
        exact=np.exp(-rho*d)*np.linalg.norm(M,2)/(4*rho**2)
        assert exact==pytest.approx(np.exp(-rho*d)/(2*rho),abs=2e-15)
        assert exact<=resolvent_separation_bound(mass,zeta,d)


def test_frozen_projector_bound_has_the_correct_units_and_angular_dependence():
    base=frozen_projector_difference_bound(.8,.2,.3,2)
    assert frozen_projector_difference_bound(.8,.2,.3,-2)==base
    assert frozen_projector_difference_bound(.8/1000,.2/1000,.3*1000,2)==pytest.approx(base)
    assert frozen_projector_difference_bound(.8,.2,.3,4)<base
    assert frozen_projector_difference_bound(.8,0.,.3,2)==0
    assert 0<frozen_projector_difference_bound(.8,.2,.3,10000)<base
    assert 0<resolvent_separation_bound(10000.,1.,1.)
    assert frozen_projector_difference_bound(1.,100.,0.,1)==1.
