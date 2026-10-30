import math
import numpy as np
import pytest
import sympy as s
from recursive_horizons.nsc_atlas_parke_nucleation import thin_wall
from recursive_horizons.nsc_atlas_smooth_horizon import (
    calculation,metric_f,horizons,infall_at,horizon_response)


def test_unsquared_junction_and_zero_tension_asymptote():
    for sigma in (.2,.02,.002,.0002):
        row=thin_wall(sigma=sigma)
        assert abs(row['angular_residual'])<3e-16
        assert abs(row['inside_constraint'])<5e-16
    t=s.symbols('t',positive=True)
    inv=(s.Rational(7,16)/(3*t))**2+s.Rational(7,48)+t*t/16
    assert s.limit(1/s.sqrt(inv),t,0,dir='+')==0
    assert s.limit(1/s.sqrt(inv)/t,t,0,dir='+')==s.Rational(48,7)


def test_common_chart_matching_requires_equal_mass_and_vacuum():
    r,dm,dl=s.symbols('r dm dl',nonzero=True)
    jump=-2*dm/r-dl*r*r/3
    required=-dl*r**3/6
    assert s.simplify(jump.subs(dm,required))==0
    assert s.simplify(s.diff(jump,r).subs(dm,required)+dl*r)==0


def test_full_einstein_tensor_of_moving_horizon_comparison():
    # Independent four-dimensional tensor derivation: detects sign/normalization errors.
    v,r,theta,phi=s.symbols('v r theta phi',real=True)
    lam=s.symbols('L',real=True);m=s.Function('m')(v)
    coords=(v,r,theta,phi);f=1-2*m/r-lam*r*r/3
    g=s.Matrix([[-f,1,0,0],[1,0,0,0],[0,0,r*r,0],[0,0,0,r*r*s.sin(theta)**2]])
    inv=g.inv();n=4
    gamma=[[[s.simplify(sum(inv[a,d]*(s.diff(g[d,c],coords[b])+s.diff(g[d,b],coords[c])
        -s.diff(g[b,c],coords[d])) for d in range(n))/2) for c in range(n)]
        for b in range(n)] for a in range(n)]
    ric=s.zeros(n)
    for a in range(n):
        for b in range(n):
            ric[a,b]=s.simplify(sum(s.diff(gamma[c][a][b],coords[c])-s.diff(gamma[c][a][c],coords[b])
                +sum(gamma[c][c][d]*gamma[d][a][b]-gamma[c][b][d]*gamma[d][a][c] for d in range(n))
                for c in range(n)))
    scalar=s.simplify(sum(inv[a,b]*ric[a,b] for a in range(n) for b in range(n)))
    assert scalar==4*lam
    stress=(ric-scalar*g/2+lam*g).applyfunc(s.simplify)
    expected=s.zeros(n);expected[0,0]=2*s.diff(m,v)/r**2
    assert stress==expected
    assert s.simplify(g.det()+r**4*s.sin(theta)**2)==0


def test_horizon_crossing_and_drain_response_against_independent_roots():
    lam=.22;m=.5
    for r in horizons(m,lam).values():
        for shift in (-1e-5,0,1e-5):
            p=infall_at(r+shift,m,lam)
            assert p['norm']==pytest.approx(-1,abs=2e-15)
        p=infall_at(r,m,lam)
        assert p['dr_dtau']==pytest.approx(-1)
        assert p['dv_dtau']==pytest.approx(.5)
    plus=horizons(m+1e-6,lam);minus=horizons(m-1e-6,lam)
    for key,r in horizons(m,lam).items():
        assert (plus[key]-minus[key])/2e-6==pytest.approx(horizon_response(r,lam),rel=3e-9)


def test_c_symmetric_motion_and_global_neutrality_are_different():
    # Arbitrary nonzero motion-dependent, C-even mixing preserves a C-even state.
    from scipy.linalg import expm
    C=np.array([[0.,1.],[1.,0.]]);Q=np.diag([1.,-1.]);rho=np.array([[.5,.1],[.1,.5]],complex)
    for angle in (.17,-.8,2.3):
        U=expm(-1j*angle*C);rho=U@rho@U.conj().T
        assert np.trace(rho@Q)==pytest.approx(0,abs=1e-15)
    # Global charge zero alone permits separated charge on two sides.
    psi=np.array([0.,1.,0.,0.]);ident=np.eye(2)
    assert psi@(np.kron(Q,ident)+np.kron(ident,Q))@psi==0
    assert psi@np.kron(Q,ident)@psi==1


def test_spherical_current_balance_distinguishes_transfer_from_creation():
    # Finite smooth conserved current in Minkowski spherical coordinates.
    r,t,R=s.symbols('r t R',positive=True)
    density=s.exp(-t);radial=r*s.exp(-t)/3
    assert s.simplify(s.diff(density,t)+s.diff(r*r*radial,r)/r**2)==0
    charge=4*s.pi*R**3*density/3
    assert s.simplify(s.diff(charge,t)+4*s.pi*R**2*radial.subs(r,R))==0


def test_offset_is_oriented_charge_and_expansion_preserves_comoving_charge():
    t,A,B,w,delta,H,nu=s.symbols('t A B w delta H nu',real=True)
    x=A*s.cos(w*t);y=B*s.cos(w*t-delta)
    Q=x*s.diff(y,t)-y*s.diff(x,t)
    assert s.trigsimp(Q-A*B*w*s.sin(delta))==0
    energy=(s.diff(x,t)**2+s.diff(y,t)**2+w*w*(x*x+y*y))/2
    assert s.trigsimp(energy-w*w*(A*A+B*B)/2)==0
    xd=s.exp(-3*H*t/2)*x.subs(w,nu);yd=s.exp(-3*H*t/2)*y.subs(w,nu)
    assert s.trigsimp(s.exp(3*H*t)*(xd*s.diff(yd,t)-yd*s.diff(xd,t))-A*B*nu*s.sin(delta))==0
    # Exact ensemble averages at fixed time, not one-cycle averaging of a changing envelope.
    kinetic_average=A*A*s.exp(-3*H*t)*(9*H*H/4+nu*nu)/4
    potential_average=A*A*s.exp(-3*H*t)*w*w/4
    assert s.simplify((kinetic_average-potential_average).subs(nu*nu,w*w-9*H*H/4))==0
    assert s.simplify((kinetic_average+potential_average).subs(nu*nu,w*w-9*H*H/4)
                      -A*A*w*w*s.exp(-3*H*t)/2)==0


def test_report_preserves_conditional_values_and_missing_inputs():
    row=calculation();v=row['admitted_vacuum']
    assert v['Lambda_Rc_squared']==pytest.approx(.3249667388661,rel=3e-12)
    assert v['horizon_m0']>row['shell_limit']['reference_radius']
    assert row['charge']['physical_eta'] is None
    assert row['radius_inversion_control']['mass_selected'] is False
    assert row['drain_comparison']['vacuum_added_by_null_flux']==0
    assert row['drain_comparison']['flux_rate_selected'] is False
    assert row['physical_eta_or_Omega_prediction'] is False
