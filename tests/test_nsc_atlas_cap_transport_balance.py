"""Full manufactured cap balances and independent piecewise mass integrals."""
import copy
import math

import numpy as np
import pytest

from control_budget import ControlBudget
from test_nsc_atlas_cap_radial_comparison import setup
from recursive_horizons.nsc_atlas_cap_transport_balance import (
    cap_transport_balance,reference_transport_cap,_integrate_pieces,
)
from recursive_horizons.nsc_atlas_background_cap import background_cap_rhs,background_cap_profile
from recursive_horizons.nsc_atlas_cap_mass_tangent import mass_tangent,mass_balance_ledger
from recursive_horizons.nsc_atlas_polar_cap import _potential

CONTROLS=ControlBudget(maximum_evaluations=2000000,cpu_seconds=60.)


@pytest.fixture
def budget():
    with CONTROLS.active() as value: yield value


@pytest.fixture(scope='module',autouse=True)
def report():
    yield
    print('full cap transport controls:',CONTROLS.budget.evaluations,
          'charged values;',CONTROLS.active_cpu_seconds,'active CPU seconds')


@pytest.mark.parametrize('R',[1/256,1/2])
def test_constant_zero_coupling_cap_matches_hand_mass_balance(R,budget):
    # Dyadic data make the constant-field limiter ties exactly representable.
    # Decimal moment/rate rounding can select different one-sided ties in
    # component-first and direction-first arithmetic; that domain guard stays.
    p,q,rx,pdot,qdot=1/4,1/8,R/8,1/16,-1/32
    faces=R*np.linspace(0,1,5)
    state=dict(Ip=p*np.diff(faces**2)/2,Iq=q*np.diff(faces))
    rates=dict(R=rx,T=1/4,Ip=(pdot/p+2*rx/R)*state['Ip'],Iq=(qdot/q+rx/R)*state['Iq'])
    metric,rates,tangent=setup(budget,radius=R,alpha=0.,state=state,rates=rates)
    result=cap_transport_balance(metric,rates,tangent,budget=budget)
    expected=q*qdot*R**3/3+p*pdot*(R**3/6-R**5/10)
    np.testing.assert_allclose([result['total'],result['native_balance']],expected,rtol=2e-10,atol=2e-14)
    assert abs(result['reference_rate_defect'])<2e-13
    assert abs(result['reference_mass_defect'])<2e-14
    assert result['actual_balance_defect']==0
    assert not result['uniform_error_enclosure']


@pytest.mark.parametrize('alpha',[0.,.002])
def test_nonconstant_full_cap_channels_recover_independent_native_balance(alpha,budget):
    metric,rates,tangent=setup(budget,cubic=True,alpha=alpha,order=4)
    before=copy.deepcopy((metric['state'],rates))
    result=cap_transport_balance(metric,rates,tangent,budget=budget)
    assert abs(result['balance_defect'])<2e-10
    assert abs(result['reference_rate_defect'])<2e-10
    assert abs(result['reference_mass_defect'])<2e-11
    assert abs(result['reference']['jump_sum_defect'])<2e-11
    assert abs(result['channels']['interfaces'])>1e-5
    assert abs(result['total']-result['channels']['interfaces']-result['native_balance'])>1e-5
    for a,b in zip(before,(metric['state'],rates)):
        for key in a: np.testing.assert_array_equal(a[key],b[key])
    assert result['reference']['sampled_A_positive']


def test_existing_coupled_rates_and_flux_ledger_match_reference_channels(budget):
    metric,_,_=setup(budget,cubic=True,alpha=.002,order=8)
    rx,tx=.02,.3
    rhs=background_cap_rhs(metric,T_x=tx,radius_x=rx,incoming_w_plus=.7,budget=budget)
    rates=dict(R=rx,T=tx,**rhs['rates'])
    tangent=mass_tangent(metric,rates,budget=budget)
    ledger=mass_balance_ledger(metric,rhs,tangent,R_x=rx,T_x=tx,budget=budget)
    result=cap_transport_balance(metric,rates,tangent,budget=budget)
    R=metric['radius']; w=metric['wall']
    row=background_cap_profile(metric,np.array([R]),budget=budget)
    pi,phi=float(row['Pi'][0]),float(row['Phi'][0])
    P=float(_potential(np.array([w['f']]),metric['alpha'])[1][0])
    rho=R*R*(P/2+w['A']*(pi*pi+phi*phi)/4)
    flux=R*R*w['A']**2*math.exp(w['delta'])*pi*phi/2
    expected=ledger['total']-rx*rho-tx*flux
    np.testing.assert_allclose(result['total'],expected,rtol=2e-10,atol=2e-11)
    assert abs(result['reference_rate_defect'])<2e-10
    assert abs(result['reference_mass_defect'])<2e-11
    assert abs(result['reference']['jump_sum_defect'])<2e-11
    print('manufactured coupled balance:',expected,'assembly defect:',result['balance_defect'])


def piece(left,right,vl,vr,owner,side='right'):
    return dict(left=left,right=right,left_rate=vl,right_rate=vr,owner=owner,
                small_radius_side=side,born=left==right)


class PolynomialView:
    def __init__(self,coefficients): self.coefficients=coefficients

    def reference_sample(self,owner,radii,*,budget,small_radius_side):
        p,c=self.coefficients[owner]
        r=np.asarray(radii)
        row=dict(Pi=np.full_like(r,p),Phi=c*r,Pi_r=np.zeros_like(r),Phi_r=np.full_like(r,c),
            f_x=np.zeros_like(r),Pi_x=np.zeros_like(r),Phi_x=np.zeros_like(r),
            Pdiff=np.zeros_like(r),Pdiff_f=np.zeros_like(r))
        budget.check(9*r.size)
        return row


def analytic_density(r,coefficients):
    p,c=coefficients
    return r*r*(1-r*r)*(p*p+c*c*r*r)/4


def analytic_flux(r,coefficients):
    p,c=coefficients
    return r**3*(1-r*r)**2*p*c/2


def analytic_mass(r,coefficients):
    p,c=coefficients
    return p*p*r**3/12+(c*c-p*p)*r**5/20-c*c*r**7/28


@pytest.mark.parametrize('where',['switch','wall','axis'])
def test_born_interface_chain_matches_piecewise_integrals(where,budget):
    coefficients=[(.2,.3),(.4,-.1),(.5,.2),(.6,-.2)]
    R=.4; tx=1.3
    if where=='switch':
        panels=[piece(0.,.01,0.,.7,0,'left'),piece(.01,.01,.7,0.,1,'left'),
                piece(.01,.01,0.,-.2,2),piece(.01,R,-.2,.05,3)]
    elif where=='wall':
        panels=[piece(0.,R,0.,.1,0),piece(R,R,.1,.05,1),piece(R,R,.05,.02,2)]
    else:
        panels=[piece(0.,0.,0.,-.1,0),piece(0.,R,-.1,.02,1)]
    result=_integrate_pieces(PolynomialView(coefficients),panels,center_delta=0.,alpha=0.,
                             clock_rate=tx,budget=budget,rtol=2e-12,atol=2e-14)
    mass=sum(analytic_mass(p['right'],coefficients[p['owner']])-
             analytic_mass(p['left'],coefficients[p['owner']]) for p in panels)
    moving=sum(p['right_rate']*analytic_density(p['right'],coefficients[p['owner']])-
               p['left_rate']*analytic_density(p['left'],coefficients[p['owner']]) for p in panels)
    last=panels[-1]; last_fields=coefficients[last['owner']]
    expected=moving-last['right_rate']*analytic_density(R,last_fields)-tx*analytic_flux(R,last_fields)
    summed=math.fsum(result[k] for k in ('Z_f','Z_Pi','Z_Phi','Z_jump'))
    np.testing.assert_allclose(result['n'],mass,rtol=2e-11,atol=2e-14)
    np.testing.assert_allclose(summed,expected,rtol=2e-10,atol=2e-13)
    assert len(result['jumps'])==len(panels)-1
    assert abs(result['jump_sum_defect'])<2e-13
    if where!='axis':
        assert result['stages'][0]['state']['n']==result['stages'][1]['state']['n']
    if where=='wall': assert all(j['r']==R for j in result['jumps'])
    if where=='axis': assert result['jumps'][0]['mass_motion']==result['jumps'][0]['clock_flux']==0
    else: assert abs(result['Z_jump'])>1e-10


def test_center_wall_switch_and_broken_partition_rejected(budget):
    metric,rates,tangent=setup(budget)
    for bad,message in ((tangent|dict(reverse=np.nan),'native reverse'),
                        (tangent|dict(reverse=tangent['reverse']+.01),'contractions differ')):
        with pytest.raises(ValueError,match=message):
            cap_transport_balance(metric,rates,bad,budget=budget)
    metric['center']['delta']+=.01
    with pytest.raises(ValueError,match='centre differs'):
        reference_transport_cap(metric,rates,budget=budget)
    metric,rates,tangent=setup(budget,radius=.01)
    before=budget.evaluations
    with pytest.raises(ValueError,match='wall on the small-radius'):
        cap_transport_balance(metric,rates,tangent,budget=budget)
    assert budget.evaluations==before
    panels=[piece(0.,.1,0.,.2,0),piece(.1,.4,.3,.1,1)]
    with pytest.raises(ValueError,match='coordinates/rates'):
        _integrate_pieces(PolynomialView([(.2,.3),(.4,-.1)]),panels,center_delta=0.,alpha=0.,
            clock_rate=1.,budget=budget,rtol=2e-11,atol=2e-13)
