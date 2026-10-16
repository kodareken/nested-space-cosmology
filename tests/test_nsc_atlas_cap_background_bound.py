"""Exact polynomial reconstruction and radial-envelope controls."""
import math

import numpy as np
from numpy.polynomial import Polynomial
from flint import arb, ctx, fmpq
import pytest
import sympy as s

from control_budget import ControlBudget
from recursive_horizons.nsc_atlas_cap_background_jets import RetainedCapBackground
from recursive_horizons.nsc_atlas_cap_background_bound import (
    background_metric_bound, _hermite_controls, _envelope_step,
)


CONTROLS=ControlBudget(maximum_evaluations=500000,cpu_seconds=30.)


@pytest.fixture
def budget():
    with CONTROLS.active() as value:
        yield value


@pytest.fixture(scope='module',autouse=True)
def receipt():
    yield
    print('background bound controls:',CONTROLS.budget.evaluations,'charged values;',
          CONTROLS.active_cpu_seconds,'active CPU seconds')


def background(poly,knots=(-1.,0.,.5,1.)):
    knots=np.array(knots)
    return RetainedCapBackground(knots,poly(knots),poly.deriv()(knots),poly.deriv(2)(knots))


def test_bernstein_controls_reconstruct_exact_quintic_and_derivatives(budget):
    x=s.symbols('x')
    P=s.Rational(1,4)+x/8-2*x*x+x**3/2+x**4/8-x**5/4
    def ball(v):
        n,d=s.fraction(v)
        return arb(fmpq(int(n),int(d)))
    with ctx.workprec(160):
        left,right=s.Rational(-1,4),s.Rational(1,2)
        inputs=[ball(s.diff(P,x,order).subs(x,point)) for point in (left,right) for order in range(3)]
        controls=_hermite_controls(*inputs,ball(right-left))
        for t in [s.Rational(0),s.Rational(1,7),s.Rational(1,2),s.Rational(1)]:
            for order,coefficients in enumerate(controls):
                degree=len(coefficients)-1
                value=sum(math.comb(degree,j)*ball(t)**j*ball(1-t)**(degree-j)*c for j,c in enumerate(coefficients))
                expected=ball(s.diff(P,x,order).subs(x,left+(right-left)*t))
                assert value.contains(expected)
    budget.check(4000)


def test_affine_envelope_maximum_is_at_an_endpoint(budget):
    r,U,b,left=s.symbols('r U b left',positive=True)
    expression=(U+b*(r-left))/r
    assert s.simplify(s.diff(expression,r)-(b*left-U)/r**2)==0
    with ctx.workprec(160):
        following,ratio=_envelope_step(arb(2),arb(3),arb(1),arb(2))
        assert following>=5 and ratio>=arb(fmpq(5,2))
        assert float(ratio)-2.5<1e-14
        following,ratio=_envelope_step(arb(0),arb(3),arb(0),arb(2))
        assert following>=6 and ratio>=3
    budget.check(500)


def test_quiet_floor_and_mass_are_exact_geometric_limits(budget):
    bg=background(Polynomial([0,0,0,1]),(0.,.5,1.))
    row=background_metric_bound(bg,-1.,.5,[0.],[0.],alpha=.002,budget=budget)
    expected=1/np.cosh(.5)**2
    assert row['mass_upper']==0.
    assert row['A_lower_float_down']<=expected
    assert expected-row['A_lower_float_down']<2e-15
    assert row['precision_bits']==160


def test_nonzero_hermite_source_and_departure_bound_geometric_mass(budget):
    Z,c=.5,.03
    bg=background(Polynomial([0,c]))
    row=background_metric_bound(bg,0.,Z,[.001],[.002],alpha=0.,budget=budget)
    # chi=(2c+.001/Z)z and k=.002 z/Z, integrated independently.
    slope=2*c+.001/Z
    beta=(1/np.cosh(Z)**2)/np.tanh(Z)
    mass=((.002/Z)**2*Z**3/3+slope*slope*Z-beta*(slope*Z)**2)/4
    assert row['mass_upper']>=mass>0.
    finite=background_metric_bound(bg,0.,Z,[.001],[.002],alpha=.002,budget=budget)
    assert 0<finite['A_lower_float_down']<row['A_lower_float_down']


def test_failed_sufficient_bound_does_not_claim_a_physical_zero(budget):
    with pytest.raises(ValueError,match='physical metric sign is undetermined'):
        background_metric_bound(background(Polynomial([0,1.])),0.,.5,[0.],[0.],alpha=1.,budget=budget)


def test_exact_domain_rejects_rounded_endpoint_escape(budget):
    bg=background(Polynomial([0,0,0,1]),(0.,.5,1.))
    # Floating sum rounds to1 although the exact sum of these binaries is >1.
    T,Z=.1,.9
    assert T+Z==1.
    with pytest.raises(ValueError,match='exact slice exits'):
        background_metric_bound(bg,T,Z,[0.],[0.],alpha=0.,budget=budget)


def test_count_domain_and_budget_controls(budget):
    bg=background(Polynomial([0.,1.]))
    for T,Z,c,k,al in [(0.,0.,[0.],[0.],0.),(0.,.5,[],[],0.),
                       (0.,.5,[0.]*17,[0.]*17,0.),(0.,.5,[0.],[0.,0.],0.),
                       (0.,.5,[0.],[0.],-1.)]:
        with pytest.raises(ValueError):
            background_metric_bound(bg,T,Z,c,k,alpha=al,budget=budget)
    maximum=budget.maximum_evaluations
    budget.maximum_evaluations=budget.evaluations+1
    try:
        with pytest.raises(RuntimeError,match='budget exhausted'):
            background_metric_bound(bg,0.,.5,[0.],[0.],alpha=0.,budget=budget)
    finally:
        budget.maximum_evaluations=maximum
