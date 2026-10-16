"""Manufactured smooth pieces, independent integrals, no retained cap."""
import mpmath as mp
import numpy as np
from numpy.polynomial import Polynomial as Poly
import pytest

from control_budget import ControlBudget
from recursive_horizons.nsc_atlas_characteristic_initialization import Budget
from recursive_horizons.nsc_atlas_cap_transport_reference import transport_radial_piece,STATE_KEYS

CONTROLS=ControlBudget(maximum_evaluations=200000,cpu_seconds=20.)


@pytest.fixture
def budget():
    with CONTROLS.active() as value: yield value


@pytest.fixture(scope='module',autouse=True)
def report():
    yield
    print('cap transport reference controls:',CONTROLS.budget.evaluations,
          'charged values;',CONTROLS.active_cpu_seconds,'active CPU seconds')


def initial(**changes):
    return dict.fromkeys(STATE_KEYS,0.)|changes


def zero_fields(r):
    return dict(f_x=0.,Pi=0.,Phi=0.,Pi_x=0.,Phi_x=0.,Pi_r=0.,Phi_r=0.,Pdiff=0.,Pdiff_f=0.)


def test_constant_field_mass_and_potential_channel_have_exact_cubic_integrals(budget):
    def fields(r):
        assert r>0
        budget.check(9)
        return zero_fields(r)|dict(f_x=.3,Pdiff=.2,Pdiff_f=.7)
    R=.4
    result=transport_radial_piece(fields,left=0.,right=R,initial=initial(delta=-.1),
        alpha=.01,clock_rate=.9,budget=budget)
    np.testing.assert_allclose(result['n'],.2*R**3/6,rtol=2e-11,atol=2e-14)
    np.testing.assert_allclose(result['Z_f'],.7*.3*R**3/6,rtol=2e-11,atol=2e-14)
    assert result['delta']==-.1
    assert [result[k] for k in ('Z_Pi','Z_Phi','Z_jump')]==[0.,0.,0.]
    assert result['coefficient_calls']<result['rhs_calls']
    assert result['sampled_A_positive'] and not result['uniform_error_enclosure']


def test_constant_attenuation_mass_lapse_and_incoming_channels_against_high_precision(budget):
    low,high,p,alpha=.1,.4,.7,.8
    a=alpha*p*p/2
    seed=initial(n=.0003,delta=-.02,Z_f=.1,Z_Pi=-.2,Z_Phi=.3,Z_jump=-.4)
    def fields(r):
        budget.check(9)
        return zero_fields(r)|dict(Pi=p/np.sqrt(r),Pi_r=-p/(2*r**1.5))
    result=transport_radial_piece(fields,left=low,right=high,initial=seed,
        alpha=alpha,clock_rate=0.,budget=budget)
    with mp.workdps(60):
        l,h,am,pm=map(mp.mpf,(low,high,a,p))
        weight=mp.exp(-am*(h-l))
        mass=weight*mp.mpf(seed['n'])+mp.quad(lambda r:mp.exp(-am*(h-r))*pm**2*r*(1-r*r)/4,[l,h])
        expected=[float(mass),seed['delta']+a*(high-low),
                  *[float(weight*seed[k]) for k in STATE_KEYS[2:]]]
    np.testing.assert_allclose([result[k] for k in STATE_KEYS],expected,rtol=2e-10,atol=2e-13)
    assert seed==initial(n=.0003,delta=-.02,Z_f=.1,Z_Pi=-.2,Z_Phi=.3,Z_jump=-.4)


def test_polynomial_bulk_channels_match_separate_exact_integrals(budget):
    r=Poly([0,1]); A=1-r*r
    f,p,phi=r*r,1+r*r,2*r
    fx,px,phix=.3+r,r+r*r,-.2*r
    delta,tx=.17,1.7
    B=A*np.exp(delta); U=2*f
    sources=[r*r*U*(fx-tx*B*p)/2,
        A*p*(r*r*px-tx*(r*r*B*phi).deriv()+tx*r*r*np.exp(delta)*U)/2,
        r*r*A*phi*(phix-tx*(B*p).deriv())/2]
    density=r*r*(f*f/2+A*(p*p+phi*phi)/4)
    density_x=r*r*U*fx/2+r*r*A*(p*px+phi*phix)/2
    flux=r*r*A*B*p*phi/2
    def fields(at):
        budget.check(9)
        return dict(f_x=fx(at),Pi=p(at),Phi=phi(at),Pi_x=px(at),Phi_x=phix(at),
            Pi_r=p.deriv()(at),Phi_r=phi.deriv()(at),Pdiff=f(at)**2,Pdiff_f=U(at))
    R=.4
    result=transport_radial_piece(fields,left=0.,right=R,initial=initial(delta=delta),
        alpha=0.,clock_rate=tx,budget=budget)
    integral=lambda poly:poly.integ()(R)-poly.integ()(0.)
    expected=[integral(poly) for poly in sources]
    np.testing.assert_allclose([result[k] for k in ('Z_f','Z_Pi','Z_Phi')],expected,rtol=2e-11,atol=2e-13)
    np.testing.assert_allclose(result['n'],integral(density),rtol=2e-11,atol=2e-13)
    np.testing.assert_allclose(sum(expected),integral(density_x)-tx*flux(R),rtol=2e-13,atol=2e-14)
    assert result['delta']==delta and result['Z_jump']==0


def test_born_piece_keeps_channels_without_querying_fields(budget):
    def no_fields(r): raise AssertionError('a born piece has no radial flow')
    seed=initial(n=.001,delta=.2,Z_jump=.3)
    result=transport_radial_piece(no_fields,left=.2,right=.2,initial=seed,alpha=.1,
                                 clock_rate=2.,budget=budget)
    assert {k:result[k] for k in STATE_KEYS}==seed
    assert result['born'] and result['rhs_calls']==result['coefficient_calls']==0
    assert result['evaluations']==8


def test_invalid_state_domain_jets_and_exhaustion_stop_locally(budget):
    options=dict(left=0.,right=.2,initial=initial(),alpha=0.,clock_rate=1.,budget=budget)
    for change,message in ((dict(initial=initial(n=1.)),'regular centre'),
                           (dict(right=1.),'0 <= left'),(dict(rtol=1e-16),'tolerances')):
        with pytest.raises(ValueError,match=message):
            transport_radial_piece(zero_fields,**(options|change))
    with pytest.raises(ValueError,match='lacks required'):
        transport_radial_piece(lambda r:{},**options)
    with pytest.raises(ValueError,match='finite'):
        transport_radial_piece(lambda r:zero_fields(r)|dict(Pi=np.nan),**options)
    with pytest.raises(RuntimeError,match='budget'):
        transport_radial_piece(zero_fields,**(options|dict(budget=Budget(maximum_evaluations=1,cpu_seconds=1.))))
