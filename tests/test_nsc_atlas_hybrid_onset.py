"""Source-to-hybrid onset and an incoming cap; no physical source replay."""
import numpy as np
import pytest
import sympy as s

from control_budget import ControlBudget
from recursive_horizons.nsc_atlas_nonlinear_boundary import hybrid_wall_boundary

CONTROLS=ControlBudget(maximum_evaluations=40000,cpu_seconds=20.)


@pytest.fixture
def budget():
    with CONTROLS.active() as value: yield value


@pytest.fixture(scope='module',autouse=True)
def report():
    yield
    print('hybrid onset controls:',CONTROLS.budget.evaluations,
          'charged values;',CONTROLS.active_cpu_seconds,'active CPU seconds')


def test_actual_source_and_scalar_prefix_leading_jet(budget):
    budget.check(12000)
    t,k,l,m=s.symbols('t k l m')
    zeta=-1+k*t+l*t**2+m*t**3
    g=s.Poly((1-zeta**2)**4,t)
    assert g.nth(4)==16*k**4
    assert all(g.nth(j)==0 for j in range(4))
    assert g.diff().nth(3)==64*k**4
    x,a,U,dV=s.symbols('x a U dV',nonzero=True)
    V=(1-x/a)*s.exp(x/a)/U
    assert s.simplify(4*V.diff(x)/dV+4*x*s.exp(x/a)/(a*a*U*dV))==0
    # Exact scalar prefix at the wall, with a fixed bounded spatial map.
    # Free fourth-order prefix and geometry coefficients must not change t^3.
    r0,rx0,rxw,R,alpha,i4,c5,g5=s.symbols('r0 rx0 rxw R alpha i4 c5 g5',nonzero=True)
    source=16*k**4*t**4+g5*t**5
    u0=source/r0;uw=c5*t**5/R
    D0=1+alpha*u0**2/6;Dw=1+alpha*uw**2/6
    cV=uw*rxw+s.sqrt(Dw)*((source.diff(t)-u0*rx0)/s.sqrt(D0)-i4*t**4)
    jet=s.Poly(s.series(cV,t,0,5).removeO(),t)
    assert all(jet.nth(j)==0 for j in range(3))
    assert jet.nth(3)==64*k**4
    assert s.simplify(jet.nth(4)-(5*g5-16*k**4*rx0/r0-i4))==0


def test_actual_hybrid_drive_and_scalar_gap_leading_coefficient(budget):
    budget.check(12000)
    t,R,A,H,Vp,alpha,c5,w4,g3,g4=s.symbols(
        't R A H Vp alpha c5 w4 g3 g4',positive=True)
    c=c5*t**5;u=c/R;D=1+alpha*u*u/6
    cV=g3*t**3+g4*t**4
    cm=A/H*(R*s.sqrt(D)*w4*t**4-u)
    wplus=(2*Vp*cV/H-u)/(R*s.sqrt(D))
    drive=A*R*wplus
    target=cV+cm/(2*Vp)
    assert s.series(drive,t,0,4).removeO()==2*A*Vp*g3*t**3/H
    gap=s.series(c.diff(t)-target,t,0,5).removeO().expand()
    assert gap.coeff(t,3)==-g3
    assert s.simplify(gap.coeff(t,4)-(5*c5-g4-A*R*w4/(2*H*Vp)))==0
    # Independent conversion of the strong boundary lag gives the same factor.
    tx,rx=s.symbols('tx rx')
    converted=-(rx+tx*A)/(2*A)*(2*A*Vp*g3/H)
    assert s.simplify(converted.subs(rx,H/Vp-tx*A)+g3)==0


def test_incoming_cap_satisfies_wave_wall_and_integrated_mass(budget):
    budget.check(14000)
    T,z,Z,front,a,b,bz=s.symbols('T z Z front a b bz',real=True)
    primitive=a*(T+z-front)**4  # Nonzero side; zero below the C3 front.
    assert s.diff(primitive,T,2)==s.diff(primitive,z,2)
    F=primitive.subs(z,Z);k=s.diff(F,T)
    lower=front-T  # This control has 0 < lower < Z: centre still quiet.
    N=s.integrate((s.diff(primitive,T)**2+s.diff(primitive,z)**2)/4,
                  (z,lower,Z))-b*F**2/4
    assert s.simplify(N-(8*a*a*(T+Z-front)**7/7-b*F**2/4))==0
    C=2*k-b*F
    assert s.simplify(N.diff(T)-k*(C-k)/2)==0
    shape=N.diff(Z)+N.diff(b)*bz
    assert s.simplify(shape-(k*k/2-b*F*k/2-bz*F*F/4))==0
    # Exercise the production hybrid independently of the symbolic formulas.
    R=.5;A=1-R*R;value=1/64;cV=1/32
    out=hybrid_wall_boundary(
        dict(r=R,c=value,eta_plus=0.),
        dict(r_V=-1/8,c_V=cV,log_C_V=0.),
        dict(A=A,delta=0.,w_minus=value/R**2),
        dict(C_wall=1.,log_C_s=0.),alpha=0.)
    tx=float(out['T_s']/out['Vdot_plus']);rx=float(out['r_s']/out['Vdot_plus'])
    label_rate=tx+rx/A;derivative=cV/label_rate
    p=derivative/A;h=p-value/R
    assert out['c_U_plus']==0
    assert out['c_s']/out['Vdot_plus']==cV
    np.testing.assert_allclose([R*out['Pi'],R*out['Phi']],[p,h],rtol=2e-15,atol=0.)
    np.testing.assert_allclose(derivative*label_rate,cV,rtol=2e-15,atol=0.)
