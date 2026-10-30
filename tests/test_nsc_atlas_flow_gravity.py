import math
import numpy as np
import pytest
import sympy as s
from recursive_horizons.nsc_atlas_flow_gravity import boundary_coefficients,horizon_flux,connections


def test_boundary_momentum_constraint_has_the_ns_normalization():
    rc=s.symbols('rc',positive=True);dv,dp,conv,lap=s.symbols('dv dp conv lap')
    # gamma^tautau=-1/rc; T_tau_i=-v_i/sqrt(rc).
    divergence=dv/rc**s.Rational(3,2)+(conv+dp)/rc**s.Rational(3,2)-lap/s.sqrt(rc)
    assert s.simplify(rc**s.Rational(3,2)*divergence-(dv+conv+dp-rc*lap))==0


def test_three_dimensional_antisymmetric_gradient_has_two_curl_norms():
    a,b,c=s.symbols('a b c',real=True)
    F=s.Matrix([[0,a,b],[-a,0,c],[-b,-c,0]])
    assert sum(x*x for x in F)==2*(a*a+b*b+c*c)
    row=boundary_coefficients(5,2,.001)
    assert row['integrated_Kretschmann_coefficient']==pytest.approx(-6e9)


def test_horizon_flux_follows_direct_scalar_derivatives_and_null_generator():
    w,m,Omega,A=s.symbols('w m Omega A',real=True)
    dt=-s.I*w*A;dphi=s.I*m*A;dk=dt+Omega*dphi
    E=2*s.re(s.conjugate(dk)*dt);J=-2*s.re(s.conjugate(dk)*dphi)
    assert s.simplify(E-2*A*A*w*(w-m*Omega))==0
    assert s.simplify(J-2*A*A*m*(w-m*Omega))==0
    assert s.simplify(E-Omega*J-2*A*A*(w-m*Omega)**2)==0
    for x in (.5,1.,1.5):
        q=horizon_flux(x,1,1)
        assert q['corotating_energy']==pytest.approx(q['null_stress_integral'])
    assert horizon_flux(.5,1,1)['energy']<0


def test_sphere_isoperimetric_equality_and_scope():
    R=s.symbols('R',positive=True)
    A=4*s.pi*R**2;V=4*s.pi*R**3/3
    assert s.simplify(A**3-36*s.pi*V**2)==0
    c=connections()
    assert c['hydrodynamic_lift']['bulk_spacetime_dimensions']==5
    assert c['hydrodynamic_lift']['inherited_4D_parent_matched'] is False
    assert c['rotating_horizon']['child_spacetime_created'] is False
