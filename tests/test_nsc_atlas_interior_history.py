"""Degree-seven history reconstruction and original-input certificate controls."""
from fractions import Fraction as Q
from math import comb

import numpy as np
from numpy.polynomial import Polynomial
from flint import ctx
import pytest

from control_budget import ControlBudget
from recursive_horizons.nsc_atlas_characteristic_background import _enriched_bernstein
from recursive_horizons.nsc_atlas_cap_background_jets import RetainedCapBackground
from recursive_horizons.nsc_atlas_cap_background_bound import _background_controls,background_metric_bound
from recursive_horizons.nsc_atlas_cap_fixed_bound import _ball,fixed_background_metric_bound
from recursive_horizons.nsc_atlas_cap_fixed_basis import FixedCapBasis
from recursive_horizons.nsc_atlas_cap_fixed_reference import fixed_cap_reference
from recursive_horizons.nsc_atlas_cap_fixed_transfer import transfer_initial_background
from recursive_horizons.nsc_atlas_cap_mass_derivatives import third_moments

CONTROLS=ControlBudget(maximum_evaluations=2000000,cpu_seconds=60.)
P=Polynomial([0,0,0,0,1/64,1/128,-1/256,1/512])


@pytest.fixture
def budget():
    with CONTROLS.active() as value: yield value


@pytest.fixture(scope='module',autouse=True)
def report():
    yield
    print('interior history controls:',CONTROLS.budget.evaluations,'values;',CONTROLS.active_cpu_seconds,'CPU seconds')


def background(enriched=True):
    labels=np.array([-.5,0.,.5,1.]);x=np.maximum(labels,0.)
    interior=dict(indices=[2],labels=[.6875],values=P([.6875]),first=P.deriv()([.6875]))
    return RetainedCapBackground(labels,P(x),P.deriv()(x),P.deriv(2)(x),
        onset=dict(index=1,quartic_coefficient=P.coef[4]),interior_jets=interior if enriched else None)


def test_rational_controls_reproduce_independent_seventh_degree(budget):
    coefficients=[Q(2,-7),Q(1,3),Q(-2,5),Q(4,9),Q(-1,8),Q(1,11),Q(3,13),Q(-5,17)]
    def jet(x,order=0):
        from math import factorial
        return sum(c*Q(factorial(i),factorial(i-order))*x**(i-order)
                   for i,c in enumerate(coefficients) if i>=order)
    eta=Q(3,7);h=Q(2,3)
    controls=_enriched_bernstein(jet(0),jet(0,1)/h,jet(0,2)/h**2,
        jet(1),jet(1,1)/h,jet(1,2)/h**2,h,eta,jet(eta),jet(eta,1)/h)
    for x in (Q(0),Q(1,8),eta,Q(5,8),Q(1)):
        assert sum(comb(7,i)*x**i*(1-x)**(7-i)*c for i,c in enumerate(controls))==jet(x)
    budget.check(5000)


def test_runtime_joins_original_onset_and_copied_inputs(budget):
    bg=background();old=background(False);queries=np.array([.51,.625,.6875,.875,1.])
    actual=bg.profile.jet(queries,budget=budget)
    for order,key in enumerate(('F','F_prime','F_second','F_third')):
        np.testing.assert_allclose(actual[key],P.deriv(order)(queries),rtol=3e-12,atol=2e-15)
    np.testing.assert_array_equal(bg.profile._spline.c[:,:2],old.profile._spline.c[:,:2])
    for order in range(3):
        assert np.polyval(np.polyder(bg.profile._spline.c[:,1],order),.5)==pytest.approx(
            np.polyval(np.polyder(bg.profile._spline.c[:,2],order),0.),abs=2e-14)
    inputs=bg.inputs;inputs['interior_values'][0]=99.
    interior=bg.interior_jets;interior['values'][0]=99.
    assert bg.interior_jets['values'][0]==P(.6875)
    assert set(old.inputs)=={'labels','values','first','second'}
    with pytest.raises(ValueError,match='must remain in the background'):
        transfer_initial_background(bg,.125,.5,budget=budget)


def test_arb_controls_and_moments_use_interior_inputs(budget):
    bg=background();data=bg.inputs
    with ctx.workprec(160):
        arrays=[list(map(_ball,data[k])) for k in ('labels','values','first','second')]
        controls=_background_controls(bg,*arrays)[2]
        for x in (Q(5,8),Q(3,4),Q(7,8)):
            t=2*x-1
            for order,row in enumerate(controls):
                n=len(row)-1
                observed=sum(comb(n,i)*_ball(t)**i*_ball(1-t)**(n-i)*c for i,c in enumerate(row))
                expected=sum(Q(float(c))*x**i for i,c in enumerate(P.deriv(order).coef))
                assert observed.contains(_ball(expected))
    T=.75;z=.0625;variable=Polynomial([T,z]);u=Polynomial([0,1])
    expected=[]
    for n in range(3):
        antiderivative=(P.deriv(3)(variable)*u**n).integ()
        expected.append(antiderivative(1)-antiderivative(-1))
    np.testing.assert_allclose(third_moments(bg.profile,T,z,budget=budget),expected,rtol=3e-12,atol=2e-15)
    budget.check(5000)


def test_mass_and_both_certificates_include_the_same_enrichment(budget):
    bg=background();basis=FixedCapBasis([0,1],budget=budget);zero=np.zeros(basis.dimension)
    T=.125;Z=Q(1,2);R=np.tanh(.5);A0=1-R*R;F=P(.625);k=P.deriv()(.625)
    row=fixed_cap_reference(basis,zero,zero,Z,alpha=0.,T_x=1.,radius_rate=0.,
        incoming_w_plus=(2*k-A0*F/R)/(A0*R),background=bg,background_time=T,budget=budget)
    integral=(P.deriv()*P.deriv()).integ()
    expected=(integral(.625)-integral(0))/2-A0*F*F/(4*R)
    assert row['mass']==pytest.approx(expected,abs=3e-14)
    assert abs(row['wall']['incoming_discrepancy'])<3e-13
    # The forward ODE has its own absolute integration error. Check the
    # independent weighted quadrature against the exact zero identity.
    adjoint_velocity=row['wall']['velocity_discrepancy']-row['projected_gradient_discrepancy']['wall_velocity']
    assert abs(adjoint_velocity)<3e-13
    # An interior bubble with every endpoint jet zero would be completely
    # invisible to an incorrect quintic-only certificate.
    bubble=RetainedCapBackground([0.,.5,1.],[0.]*3,[0.]*3,[0.]*3,
        interior_jets=dict(indices=[0],labels=[.25],values=[1/1024],first=[0.]))
    assert background_metric_bound(bubble,.25,.125,[0.],[0.],alpha=.002,budget=budget)['mass_upper']>0.
    bound=fixed_background_metric_bound(basis,zero,zero,Z,background=bubble,T=.25,alpha=.002,budget=budget)
    assert bound['mass_upper']>0.
    with pytest.raises(ValueError,match='enriched history'):
        transfer_initial_background(bubble,.25,.5,budget=budget)


def test_invalid_interior_queries_and_quiet_onset_replacement_fail(budget):
    labels=[-.5,0.,.5,1.];values=[0.,0.,.1,.2]
    valid=dict(indices=[2],labels=[.75],values=[.15],first=[.2])
    def build(interior):
        return RetainedCapBackground(labels,values,[0.]*4,[0.]*4,
            onset=dict(index=1,quartic_coefficient=1.),interior_jets=interior)
    for change in ({'indices':[1]},{'indices':[True]},{'indices':[2.5]},
                   {'labels':[.5]},{'labels':[1.]},{'values':[float('nan')]},
                   {'first':[]},{'indices':[2,2]}):
        with pytest.raises(ValueError): build(dict(valid,**change))
    budget.check(1000)
