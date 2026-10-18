"""Independent tetrad, norm, squared-Hamiltonian and interface controls."""
import math

import numpy as np
import pytest
import sympy as s

from recursive_horizons.nsc_atlas_dirac_probe import (
    cap_radial_dirac,quiet_reference_curvatures,reference_dirac_interface,
)


def test_direct_static_tetrad_operator_after_spinor_rescaling():
    r,h=s.symbols('r h',positive=True)
    A=1-h*h*r*r
    rescale=A**(-s.Rational(1,4))/r
    u=1+r+2*r**3
    # The independent orthonormal-tetrad radial term, multiplied by sqrt(A).
    tetrad=A*(s.diff(rescale*u,r)+(1/r+s.diff(A,r)/(4*A))*rescale*u)
    assert s.simplify(tetrad/rescale-A*s.diff(u,r))==0
    # Original Dirac norm r^2/sqrt(A) dr becomes dz=dr/A.
    assert s.simplify(r*r/s.sqrt(A)*rescale**2-1/A)==0


def test_squared_first_order_hamiltonian_has_both_partner_potentials():
    z=s.symbols('z',positive=True)
    W=1/s.sinh(z)
    sigma1=s.Matrix([[0,1],[1,0]])
    sigma2=s.Matrix([[0,-s.I],[s.I,0]])
    u=s.Matrix([z+z**3,1+2*z*z])
    def H(v): return -s.I*sigma2*v.diff(z)+sigma1*W*v
    expected=-u.diff(z,2)+s.diag(W*W-s.diff(W,z),W*W+s.diff(W,z))*u
    assert all(s.simplify(v)==0 for v in H(H(u))-expected)
    r=math.tanh(.6)
    row=cap_radial_dirac(r)
    assert row['W']==pytest.approx(1/math.sinh(.6))
    assert row['dW_dz']==pytest.approx(-math.cosh(.6)/math.sinh(.6)**2)


def test_flat_and_angular_sign_limits():
    for k in (-3,-1,1,2):
        row=cap_radial_dirac(.5,hubble=0,angular_k=k)
        assert row['V_plus']==pytest.approx(k*(k-1)/.5**2)
        assert row['V_minus']==pytest.approx(k*(k+1)/.5**2)
        curved=cap_radial_dirac(.5,angular_k=k)
        opposite=cap_radial_dirac(.5,angular_k=-k)
        assert curved['V_plus']==opposite['V_minus']


def test_wall_domain_from_direct_warped_tetrad_not_scalar_boundary_law():
    eta=s.symbols('eta',real=True)
    # In Gaussian normal form h_ab=e^(2a eta) h0_ab, K=3a.
    # A tangentially constant spinor exp(-3a eta/2) obeys
    # (d_eta+K/2)psi=0 separately on both sides.
    left,right=s.Rational(1,7),s.Rational(1,7)-2
    derivatives=[]
    for a in (left,right):
        psi=s.exp(-3*a*eta/2)
        assert s.diff(psi,eta)+3*a*psi/2==0
        derivatives.append(s.diff(psi,eta).subs(eta,0))
    assert derivatives[1]-derivatives[0]==3
    for r in (.45,.5,.55):
        row=reference_dirac_interface(r)
        assert row['trace_K_jump']==pytest.approx(-6)
        assert row['squared_Dirac_normal_derivative_jump_coefficient']==pytest.approx(3)
        assert row['angular_inverse_radius_normal_jump']==pytest.approx(2/r)
        assert not row['exterior_static_time_is_timelike']


def test_curvatures_match_einstein_cap_and_vacuum_exterior():
    row=quiet_reference_curvatures(.5)
    cap=row['cap']
    assert cap['Ricci_squared']==cap['R']**2/4
    assert cap['Riemann_squared']-2*cap['Ricci_squared']+cap['R']**2/3==0
    assert row['exterior']['Weyl_squared']==pytest.approx(491.52)


@pytest.mark.parametrize('kwargs',[{'radius':0},{'radius':1},{'radius':.5,'angular_k':0},
    {'radius':.5,'angular_k':True},{'radius':.5,'hubble':np.inf}])
def test_invalid_radial_domain(kwargs):
    with pytest.raises(ValueError): cap_radial_dirac(**kwargs)
