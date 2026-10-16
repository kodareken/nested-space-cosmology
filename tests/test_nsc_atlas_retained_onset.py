"""Independent polynomial, moment and constrained-mass onset controls."""
from fractions import Fraction as Q
from math import comb

import mpmath as mp
import numpy as np
from numpy.polynomial import Polynomial
from flint import ctx
import pytest

from control_budget import ControlBudget
from recursive_horizons.nsc_atlas_cap_background_jets import RetainedCapBackground
from recursive_horizons.nsc_atlas_cap_background_bound import background_metric_bound,_background_controls
from recursive_horizons.nsc_atlas_cap_fixed_bound import _ball
from recursive_horizons.nsc_atlas_cap_fixed_basis import FixedCapBasis
from recursive_horizons.nsc_atlas_cap_fixed_reference import fixed_cap_reference
from recursive_horizons.nsc_atlas_cap_fixed_transfer import transfer_initial_background
from recursive_horizons.nsc_atlas_cap_mass_derivatives import third_moments

CONTROLS=ControlBudget(maximum_evaluations=2000000,cpu_seconds=30.)
P=Polynomial([0,0,0,0,1/64,1/128,-1/256,1/512])


@pytest.fixture
def budget():
    with CONTROLS.active() as value: yield value


@pytest.fixture(scope='module',autouse=True)
def report():
    yield
    print('retained onset controls:',CONTROLS.budget.evaluations,
          'charged values;',CONTROLS.active_cpu_seconds,'active CPU seconds')


def background(front=0.,width=.5):
    labels=front+np.array([-width,0.,width,2*width])
    x=np.maximum(labels-front,0.)
    return RetainedCapBackground(labels,P(x),P.deriv()(x),P.deriv(2)(x),
                                 onset=dict(index=1,quartic_coefficient=P.coef[4]))


def test_one_sided_front_and_right_join_reproduce_independent_polynomial(budget):
    bg=background();profile=bg.profile
    query=np.array([0.,1e-10,.03,.125,.5])
    jets=profile.jet(query,budget=budget)
    for order,key in enumerate(('F','F_prime','F_second','F_third')):
        # At the next knot the public third jet selects the following C2
        # quintic; equality of third derivatives there is not required.
        keep=slice(None,-1) if order==3 else slice(None)
        np.testing.assert_allclose(jets[key][keep],P.deriv(order)(query[keep]),rtol=2e-12,atol=1e-40)
    quiet=profile.jet(np.array([-.3,-1e-10]),budget=budget)
    assert all(np.all(v==0.) for v in quiet.values())
    # The piece to the left stays identically zero, including its fourth jet.
    assert np.all(profile._spline.c[:,0]==0.)
    assert profile._spline.derivative(4)(0.)==pytest.approx(24*P.coef[4],rel=2e-15)
    for order in range(3):
        left=np.polyval(np.polyder(profile._spline.c[:,1],order),.5)
        right=np.polyval(np.polyder(profile._spline.c[:,2],order),0.)
        assert left==pytest.approx(right,abs=2e-15)
    assert np.polyval(np.polyder(profile._spline.c[:,1],3),.5)==pytest.approx(P.deriv(3)(.5),abs=2e-15)
    onset=bg.onset;onset['quartic_coefficient']=99.
    data=bg.inputs;data['values'][:]=99.
    assert bg.onset['quartic_coefficient']==P.coef[4]
    assert bg.inputs['values'][0]==0.
    with pytest.raises(ValueError,match='must remain in the background'):
        transfer_initial_background(bg,0.,.25,budget=budget)


def test_degree_seven_moments_cross_quiet_front_exactly(budget):
    bg=background();T=.01;distance=.03
    # Integrate independent physical polynomials after substitution x=T+d*u.
    variable=Polynomial([T,distance]);cut=-T/distance
    def integral(poly):
        primitive=poly.integ();return primitive(1.)-primitive(cut)
    first=P.deriv()(variable);second=P.deriv(2)(variable);third=P.deriv(3)(variable)
    expected=[integral(first),integral(second),integral(second*Polynomial([0,1]))]
    observed=bg.profile._split_integrals(T,distance,budget)
    np.testing.assert_allclose(observed,expected,rtol=2e-12,atol=1e-22)
    expected=[integral(third*Polynomial([0]*k+[1])) for k in range(3)]
    np.testing.assert_allclose(third_moments(bg.profile,T,distance,budget=budget),expected,rtol=2e-12,atol=1e-22)


@pytest.mark.parametrize('z',[1e-10,.01,np.nextafter(.01,np.inf)])
def test_centre_and_switch_jets_against_80_digit_polynomial(z,budget):
    bg=background();row=bg.reduced_jets(.125,z,budget=budget)
    with mp.workdps(80):
        F=lambda x:sum(mp.mpf(float(c))*x**i for i,c in enumerate(P.coef))
        def fields(t,z):
            r=mp.tanh(z);A=1/mp.cosh(z)**2
            f=(F(t+z)-F(t-z))/r
            return dict(f=f,p=(mp.diff(F,t+z)-mp.diff(F,t-z))/A,
                        h=(mp.diff(F,t+z)+mp.diff(F,t-z))/A-f)
        t,zz=mp.mpf('.125'),mp.mpf(float(z))
        for key in ('f','p','h'):
            expected={key:fields(t,zz)[key],key+'_T':mp.diff(lambda q:fields(q,zz)[key],t),
                      key+'_z':mp.diff(lambda q:fields(t,q)[key],zz)}
            for name,value in expected.items():
                assert abs(mp.mpf(float(row[name]))-value)<=mp.mpf('3e-10')*abs(value),(name,z,value,row[name])
    budget.check(1000)


def test_original_input_bernstein_hulls_reconstruct_the_same_seventh_degree_piece(budget):
    bg=background();data=bg.inputs
    with ctx.workprec(160):
        labels,values,first,second=([_ball(v) for v in data[k]] for k in ('labels','values','first','second'))
        controls=_background_controls(bg,labels,values,first,second)[1]
        for x in (Q(0),Q(1,16),Q(1,4),Q(1,2)):
            t=x/Q(1,2)
            for order,coefs in enumerate(controls):
                n=len(coefs)-1
                value=sum(comb(n,i)*_ball(t)**i*_ball(1-t)**(n-i)*c for i,c in enumerate(coefs))
                derivative=P.deriv(order)
                exact=sum(Q(float(c))*x**i for i,c in enumerate(derivative.coef))
                assert value.contains(_ball(exact))
    budget.check(5000)


def test_actual_mass_time_derivative_and_both_metric_certificates(budget):
    bg=background(.5,.25);basis=FixedCapBasis([0,1],budget=budget)
    zero=np.zeros(basis.dimension);T=.125;Z=Q(1,2);d=T+float(Z)-.5
    R=np.tanh(float(Z));A0=1-R*R;b=A0/R;F=P(d);k=P.deriv()(d)
    opts=dict(T_x=1.,radius_rate=0.,incoming_w_plus=(2*k-b*F)/(A0*R),budget=budget,background=bg)
    row=fixed_cap_reference(basis,zero,zero,Z,alpha=0.,background_time=T,**opts)
    primitive=(P.deriv()*P.deriv()).integ()
    exact=(primitive(d)-primitive(0))/2-b*F*F/4
    assert row['mass']==pytest.approx(exact,abs=2e-14)
    assert row['gradient']['T']==pytest.approx(k*(k-b*F)/2,abs=2e-13)
    assert abs(row['explicit_time_remainder'])<2e-13
    finite=fixed_cap_reference(basis,zero,zero,Z,alpha=.002,background_time=T,**opts)
    h=1e-4
    plus=fixed_cap_reference(basis,zero,zero,Z,alpha=.002,background_time=T+h,**opts)
    minus=fixed_cap_reference(basis,zero,zero,Z,alpha=.002,background_time=T-h,**opts)
    assert (plus['mass']-minus['mass'])/(2*h)==pytest.approx(finite['gradient']['T'],abs=2e-12)
    assert abs(finite['remainder_identity_defect'])<2e-13
    other=background_metric_bound(bg,T,float(Z),[0.],[0.],alpha=.002,budget=budget)
    for certificate in (finite['metric_certificate'],other):
        assert certificate['mass_upper']>=finite['mass']>0.
        assert finite['sampled_minimum_A']>=certificate['A_lower_float_down']>0.


def test_invalid_onset_contract_is_rejected(budget):
    labels=np.array([0.,.5,1.]);jets=[P.deriv(i)(labels) for i in range(3)]
    for onset in ({},{'index':True,'quartic_coefficient':1.},
                  {'index':2,'quartic_coefficient':1.},{'index':0,'quartic_coefficient':0.},
                  {'index':0,'quartic_coefficient':np.inf},{'index':1,'quartic_coefficient':1.}):
        with pytest.raises(ValueError): RetainedCapBackground(labels,*jets,onset=onset)
    with pytest.raises(ValueError,match='unresolved'):
        RetainedCapBackground([0.,1e-50],[0.,0.],[0.,0.],[0.,0.],
                              onset=dict(index=0,quartic_coefficient=1e-300))
    with pytest.raises(ValueError,match='supplied finite'):
        background().profile.jet(1.01,budget=budget)
