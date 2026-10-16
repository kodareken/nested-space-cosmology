"""Independent polar PDE, incoming wave and common wall port controls."""
import mpmath as mp
import numpy as np
import pytest
import sympy as s

from control_budget import ControlBudget
from recursive_horizons.nsc_atlas_cap_characteristics import (
    characteristic_point_rhs,characteristic_wall_port,_frame_gap)
from recursive_horizons.nsc_atlas_characteristic_initialization import Budget

CONTROLS=ControlBudget(maximum_evaluations=100000,cpu_seconds=30.)


@pytest.fixture
def budget():
    with CONTROLS.active() as value: yield value


@pytest.fixture(scope='module',autouse=True)
def report():
    yield
    print('cap characteristic controls:',CONTROLS.budget.evaluations,
          'charged values;',CONTROLS.active_cpu_seconds,'active CPU seconds')


def test_conservative_transform_and_primitive_constraint(budget):
    budget.check(18000)
    r=s.symbols('r');p,q,B,A0,S=s.symbols('p q B A0 S',cls=s.Function)
    p,q,B,A0,S=(fn(r) for fn in (p,q,B,A0,S))
    pT=s.diff(B*q,r)-S;qT=s.diff(B*p,r)
    beta=B/A0
    assert s.simplify(A0*(pT+qT)-A0*s.diff(beta*A0*(p+q),r)+A0*S)==0
    assert s.simplify(A0*(pT-qT)+A0*s.diff(beta*A0*(p-q),r)+A0*S)==0
    # d/dz = A0 d/dr. The two identical sources cancel exactly.
    plusT=A0*s.diff(B*(p+q),r)-A0*S
    minusT=-A0*s.diff(B*(p-q),r)-A0*S
    assert s.simplify(A0*s.diff(B*p,r)-(plusT-minusT)/2)==0


def test_exact_incoming_wave_has_no_outgoing_source_at_zero_coupling(budget):
    z=np.array([.4,.5,.6]);R=np.tanh(z);d=.1+z-.35
    F=d**4;Fp=4*d**3;Fpp=12*d*d
    out=characteristic_point_rhs(R,F,2*Fp,0.,Fp,2*Fpp,0.,mass=0.,lapse=0.,
                                 alpha=0.,clock_rate=1.,budget=budget)
    np.testing.assert_allclose(out['chi_rate'],Fp,rtol=2e-15,atol=0.)
    np.testing.assert_allclose(out['E_plus_rate'],2*Fpp,rtol=2e-15,atol=0.)
    assert np.all(out['E_minus_rate']==0.)
    assert np.all(out['primitive_residual']==0.)
    assert np.all(out['beta_z']==0.)
    assert np.all(out['source']==0.)


def test_finite_coupling_matches_independent_polar_equations(budget):
    budget.check(6000)
    with mp.workdps(70):
        r,alpha,n,delta=map(mp.mpf,('.37','.04','.002','.03'))
        f=mp.mpf('.1')+mp.mpf('.2')*r*r
        p=mp.mpf('.03')*r+mp.mpf('.01')*r**3
        q=mp.mpf('.1')+mp.mpf('.6')*r*r
        pr=mp.mpf('.03')+mp.mpf('.03')*r*r;qr=mp.mpf('1.2')*r
        A0=1-r*r;A=A0-2*alpha*n/r;mu=mp.exp(delta);B=A*mu
        af=mp.sqrt(alpha/6)*f
        P1=(3*mp.cosh(af)**4-3)/alpha
        U=2*f*(mp.sinh(af)/af)*mp.cosh(af)**3
        h=q-f;rho=r*r*P1/2+A*(p*p+h*h)/4
        Ar=-2*r-2*alpha*(rho/r-n/(r*r))
        deltar=alpha*(p*p+h*h)/(2*r)
        Br=mu*(Ar+A*deltar)
        pT=Br*q+B*qr-Br*f-r*mu*U
        qT=Br*p+B*pr
        ep=A0*(p+q);em=A0*(p-q)
        epz=A0*(-2*r*(p+q)+A0*(pr+qr))
        emz=A0*(-2*r*(p-q)+A0*(pr-qr))
        expected=list(map(float,[B*p,A0*(pT+qT),A0*(pT-qT),rho,A*B*p*h/2]))
        values=list(map(float,[r,r*f,ep,em,A0*q,epz,emz]))
    out=characteristic_point_rhs(*values,mass=float(n),lapse=float(delta),alpha=float(alpha),
                                 clock_rate=1.,budget=budget)
    keys=('chi_rate','E_plus_rate','E_minus_rate','mass_density','temporal_mass_target')
    np.testing.assert_allclose([out[k] for k in keys],expected,rtol=4e-13,atol=2e-16)
    assert out['source']!=0.
    assert out['numerical_error_enclosure'] is False
    shifted=characteristic_point_rhs(*values[:4],values[4]+.01,*values[5:],
        mass=float(n),lapse=float(delta),alpha=float(alpha),clock_rate=2.,coordinate_rate=.1,budget=budget)
    np.testing.assert_allclose(shifted['primitive_residual'],.01,rtol=0.,atol=2e-17)
    np.testing.assert_allclose([shifted[k] for k in keys[:3]],
        2*np.array(expected[:3])+.1*np.array([values[4]+.01,values[5],values[6]]),
        rtol=4e-13,atol=2e-16)


def test_wall_port_shares_value_and_reconstructs_both_hybrid_traces(budget):
    R=.48;chi=.013;em=-.004;alpha=.02;n=.001;delta=.01;calls=[]
    def exterior(common):
        calls.append(dict(common))
        with pytest.raises(TypeError): common['c']=0.
        return dict(r_V=-.11,c_V=.007,log_C_V=.014)
    out=characteristic_wall_port(R,chi,em,mass=n,lapse=delta,alpha=alpha,
        exterior_wall=exterior,C_wall=.9,eta_plus=.03,eta_plus_x=-.01,budget=budget)
    assert len(calls)==1
    with mp.workdps(60):
        c=float(mp.mpf(R)*mp.sinh(mp.sqrt(mp.mpf(alpha)/6)*mp.mpf(chi)/R)/mp.sqrt(mp.mpf(alpha)/6))
    np.testing.assert_allclose(calls[0]['c'],c,rtol=2e-15)
    hy=out['hybrid'];rates=out['wall_rates'];A=1-R*R-2*alpha*n/R
    np.testing.assert_allclose([out['p'],out['h']],R*np.array([hy['Pi'],hy['Phi']]),rtol=4e-14,atol=1e-17)
    np.testing.assert_allclose(out['cap_c_x'],rates['c'],rtol=4e-14,atol=1e-17)
    np.testing.assert_allclose(out['temporal_mass_target'],A*A*np.exp(delta)*out['p']*out['h']/2,rtol=3e-15)
    np.testing.assert_allclose(out['moving_mass_target'],rates['T']*out['temporal_mass_target']+
                              rates['r']*out['mass_density'],rtol=3e-15)


def test_quiet_wall_and_invalid_domains_and_budget(budget):
    kwargs=dict(mass=0.,lapse=0.,alpha=0.,clock_rate=1.,budget=budget)
    for r in (0.,1.,True,np.nan):
        with pytest.raises(ValueError): characteristic_point_rhs(r,0.,0.,0.,0.,0.,0.,**kwargs)
    for update in (dict(mass=10.,alpha=1.),dict(alpha=True),dict(lapse=1000.),dict(clock_rate=0.)):
        with pytest.raises(ValueError): characteristic_point_rhs(.4,0.,0.,0.,0.,0.,0.,**(kwargs|update))
    with pytest.raises(RuntimeError,match='budget'):
        characteristic_point_rhs(.4,0.,0.,0.,0.,0.,0.,**(kwargs|dict(budget=Budget(maximum_evaluations=1,cpu_seconds=1.))))
    wall=dict(mass=0.,lapse=0.,alpha=0.,exterior_wall=lambda common:dict(r_V=-.1,c_V=0.,log_C_V=0.),
              C_wall=1.,eta_plus=0.,eta_plus_x=0.,budget=budget)
    quiet=characteristic_wall_port(.5,0.,0.,**wall)
    assert quiet['E_plus']==quiet['wall_rates']['chi']==quiet['cap_c_x']==0.
    assert quiet['temporal_mass_target']==quiet['moving_mass_target']==0.
    with pytest.raises(ValueError,match='three outgoing'):
        characteristic_wall_port(.5,0.,0.,**(wall|dict(exterior_wall=lambda common:dict(c_V=0.))))
    with pytest.raises(ValueError,match='positive A'):
        characteristic_wall_port(.5,0.,0.,**(wall|dict(mass=100.,alpha=1.)))


def test_small_frame_correction_survives_leading_term_cancellation(budget):
    budget.check(2000)
    for f,alpha in ((.2,1e-12),(-.2,.02),(.2,6.)):
        with mp.workdps(80):
            ff,al=mp.mpf(f),mp.mpf(alpha)
            q=mp.sqrt(al/6)*ff
            expected=float(ff*(mp.sinh(q)/q-mp.cosh(q)))
        np.testing.assert_allclose(_frame_gap(f,alpha),expected,rtol=4e-14,atol=0.)
    assert _frame_gap(.2,0.)==0.
