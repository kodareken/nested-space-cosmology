"""Independent forced equations, CK closure and normalization of critical data."""
import math

import numpy as np
import pytest
from scipy.special import eval_gegenbauer,roots_chebyu,roots_jacobi
from threadpoolctl import threadpool_limits

from recursive_horizons.nsc_atlas_shape_transition import critical_parent_response
from recursive_horizons.nsc_atlas_closed_pocket import solve_critical_pocket_data,PocketValve


@pytest.fixture(scope='module')
def response():
    with threadpool_limits(limits=1):return critical_parent_response(budget_seconds=30.)


def test_actual_BVP_polynomial_derivatives_solve_independently_projected_forcing(response):
    r=response;T=np.array(r['audit']['T']);q=np.array(r['audit']['carrier_q']);h=np.array(r['audit']['unit_q0'])
    u,w=roots_chebyu(8);w/=math.pi/2
    U=np.array([eval_gegenbauer(l,1,u) for l in range(5)])
    chi=q@U[:3];pref=6*r['model']['mu2']/(math.sqrt(84)*r['H2'])
    sources=[pref*(chi**2*w)@U.T,2*pref*(h[:,None]*chi*w)@U[:3].T,
             pref*(h[:,None]**2*np.ones_like(chi)*w)@U[:1].T]
    for key,source in zip(('base','first_response','second_response'),sources):
        block=r['controls'][key];a=block['audit'];labels=np.array(block['channels'])
        z=np.array(a['z']);zp=np.array(a['z_T']);zpp=np.array(a['z_TT'])
        np.testing.assert_allclose(a['z_T_from_polynomial'],zp,rtol=0,atol=5e-9)
        residual=zpp+3*np.tanh(T)[:,None]*zp+(
            labels*(labels+2)/np.cosh(T)[:,None]**2-r['model']['alpha'])*z-source
        np.testing.assert_allclose(residual/(1+abs(source)),0.,rtol=0,atol=3e-8)
        for i,l in enumerate(labels):
            assert (zp[0,i] if l%2==0 else z[0,i])==pytest.approx(0.,abs=2e-11)
        np.testing.assert_allclose(zp[-1]+(3+r['model']['p'])*z[-1],0.,rtol=0,atol=2e-11)


def test_critical_coefficients_keep_units_symmetry_and_controlled_reference_range(response):
    r=response;assert r['admitted']
    initial=r['audit']['carrier_initial'];A=r['carrier_A'];B=r['carrier_B'];H=r['H']
    np.testing.assert_allclose(initial,[0.,A/2,0.,-B*math.sqrt(.6)/(12*H),0.,-B*math.sqrt(.6)/(4*H)],rtol=0,atol=2e-17)
    S=np.array(r['canonical_S_critical_U']);Pi=np.array(r['canonical_PiS_critical_U'])
    np.testing.assert_array_equal(S[[1,3]],0.);np.testing.assert_array_equal(Pi[[0,2,4]],0.)
    np.testing.assert_allclose(S,[.872748522264,0.,.118690949915,0.,-.000180420861],rtol=0,atol=2e-7)
    np.testing.assert_allclose(Pi,[0.,3.93448553253e-5,0.,-.000460536214302,0.],rtol=0,atol=2e-8)
    assert r['delta_PiChi']==pytest.approx(-.0063847845374,abs=2e-8)
    assert r['controls']['cutoff_initial_S_max_spread']<1e-7
    assert r['controls']['cutoff_initial_PiS_max_spread']<1e-7


def test_quadratic_matches_direct_CK_and_the_weak_zero_source_branch(response):
    r=response;basis=r['response_basis'];c=r['CK_coefficients']
    u,w=roots_jacobi(12,1.5,1.5);w/=sum(w)
    U=np.array([eval_gegenbauer(l,1,u) for l in range(5)])
    Du=np.array([np.zeros_like(u) if l==0 else 2*eval_gegenbauer(l-1,2,u) for l in range(5)])
    S0=np.array(basis['S0_U']);P0=np.array(basis['Pi0_U']);S1=np.array(basis['S1_U']);P1=np.array(basis['Pi1_U']);S2=np.array(basis['S2_U'])
    independent=[w@((P0@U)*(S0@Du)),w@((P1@U)*(S0@Du)+(P0@U)*(S1@Du)),w@((P1@U)*(S1@Du))]
    np.testing.assert_allclose([c['C0'],c['C1'],c['C2']],independent,rtol=2e-13,atol=3e-18)
    for delta in (r['delta_PiChi'],.003):
        S=S0+delta*S1+delta**2*S2;Pi=P0+delta*P1
        PiChi=np.array([-r['carrier_B']*math.sqrt(.6)/12+delta,0.,-r['carrier_B']*math.sqrt(.6)/4])
        direct=w@((Pi@U)*(S@Du))+r['carrier_A']*(w@(PiChi@U[:3]))
        polynomial=c['C2']*delta**2+(r['carrier_A']+c['C1'])*delta+c['C0']
        assert direct==pytest.approx(polynomial,abs=4e-18)
    assert abs(r['CK_direct_weighted_residual'])<2e-16
    assert c['minimum_reduced_discriminant_along_carrier_scaling']>0
    lam=.01;A=r['carrier_A'];C0,C1,C2=c['C0'],c['C1'],c['C2']
    D=(A+lam**2*C1)**2-4*lam**4*C2*C0
    weak=-2*lam**3*C0/(A+lam**2*C1+math.sqrt(D))
    assert weak/lam**3==pytest.approx(-C0/A,rel=1e-4)
    assert abs(weak)<abs(r['delta_PiChi'])/100000


def test_invalid_input_fails_before_a_response_is_computed():
    for value in (-.1,float('nan'),True):
        with pytest.raises(ValueError,match='carrier scale'):critical_parent_response(value)
    with pytest.raises(ValueError,match='allowance'):critical_parent_response(budget_seconds=0.)


@pytest.fixture(scope='module')
def constrained(response):
    with threadpool_limits(limits=1):
        return solve_critical_pocket_data(40,critical_record=response,compression=.01)


def test_joint_factory_uses_closed_history_verbatim_and_rejects_another_carrier_correction(response,constrained):
    data=constrained;normal=math.sqrt(2*math.pi**2)
    expected=np.array(response['canonical_PiS_critical_U']);expected[2]-=math.sqrt(84)*.01/3
    np.testing.assert_array_equal(data.momentum_S_coefficients,normal*expected)
    np.testing.assert_array_equal(data.momentum_carrier_coefficients,
        normal*np.array(response['updated_PiChi_conformal_U']))
    assert response['canonical_PiS_critical_U'][2]==0.  # Supplied shared record was not changed.
    bad=dict(response);bad['updated_PiChi_conformal_U']=list(response['updated_PiChi_conformal_U'])
    bad['updated_PiChi_conformal_U'][0]+=1e-4
    with pytest.raises(ValueError,match='no second correction'):
        solve_critical_pocket_data(40,critical_record=bad,compression=.01)


def test_joint_variable_momentum_York_and_dense_physical_constraints(constrained):
    d=constrained;u,w=roots_jacobi(401,.5,.5);w*=4*math.pi;v=1-u*u
    H=d.parent.H;normal=math.sqrt(2*math.pi**2);psi=d.values(u);pu=d.values(u,1);puu=d.values(u,2)
    state=d.state_values(u);S=d.S_values(u);Su=d.S_values(u,1);PiS,PiChi=state[6:8]
    labels=d.York_primitive_labels;coef=d.York_primitive_coefficients/normal
    cover=d.parent.kappa*sum(c*8*eval_gegenbauer(int(l)-2,3,u)/(l*(l+2)-3)
        for c,l in zip(coef,labels) if l>=2)
    cover_u=d.parent.kappa*sum(c*48*eval_gegenbauer(int(l)-3,4,u)/(l*(l+2)-3)
        for c,l in zip(coef,labels) if l>=3)
    Ar=v*cover;Aru=-2*u*cover+v*cover_u
    J=(PiS*psi**6)*Su+(PiChi*psi**6)*d.carrier_A
    divergence=H*(-np.sqrt(v)*Aru+3*u/np.sqrt(v)*Ar)
    np.testing.assert_allclose(divergence,d.parent.kappa*H*np.sqrt(v)*J,rtol=3e-11,atol=5e-13)
    # The primitive derivative tests BOTH variable momenta directly; it does
    # not reuse the constant-PiS formula of the older preparation.
    F_u=sum(c*2*eval_gegenbauer(int(l)-1,2,u) for c,l in zip(coef,labels) if l>=1)
    np.testing.assert_allclose(F_u,J,rtol=3e-11,atol=5e-13)
    Kr=state[2]+v*state[3];Kt=state[2]
    Ktu=-.5*Aru*psi**-6+3*Ar*psi**-7*pu
    s=d.sstar+S/math.sqrt(84);chi=d.carrier_A*u
    potential=7*np.exp(-12*s)-7.5*np.exp(-14*s)+25/32*np.exp(-24*s)+200*np.exp(-12*s)*chi**2
    grad2=H*H*v*psi**-4*(Su**2+d.carrier_A**2)
    rho=.5*(PiS**2+PiChi**2+grad2)+potential
    R3=H*H*psi**-5*(-8*(v*puu-3*u*pu)+6*psi)
    ham=R3+(Kr+2*Kt)**2-Kr**2-2*Kt**2-2*d.parent.kappa*rho
    q=psi**2/H;j=np.sqrt(v)/q*(PiS*Su+PiChi*d.carrier_A)
    momentum=2*((u/np.sqrt(v)-2*np.sqrt(v)*pu/psi)*(Kr-Kt)+np.sqrt(v)*Ktu)/q-d.parent.kappa*j
    measure=w*psi**6/H**3;rms=lambda x:math.sqrt(float(measure@(x*x)/sum(measure)))/H**2
    assert rms(ham)<2e-8 and rms(momentum)<5e-10
    lower=(d.coefficients[0]-np.dot(d.labels[1:]+1,abs(d.coefficients[1:])))/normal
    assert lower>0 and min(psi)>=lower-1e-14
    system=PocketValve(64,initial_data=d)
    assert not system.reflection_symmetry and system.calls==0
    np.testing.assert_array_equal(system.initial,d.state_values(system.u))
