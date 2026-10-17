"""Independent algebra and mixed-equation controls; no original field run."""
from fractions import Fraction as Q

import pytest
import sympy as s

from control_budget import ControlBudget
from recursive_horizons import nsc_atlas_scalar_departure_bound as m

CONTROLS=ControlBudget(maximum_evaluations=12000,cpu_seconds=20.)


@pytest.fixture
def budget():
    with CONTROLS.active() as value:yield value


@pytest.fixture(scope='module',autouse=True)
def report():
    yield
    print('scalar departure controls:',CONTROLS.budget.evaluations,
          'charged values;',CONTROLS.active_cpu_seconds,'active CPU seconds')


def inputs():
    return dict(duration=Q(1,4),width=2,kernel=Q(1,2),source_X=1,
                source_Q=3,wall_P=2,kernel_pointwise=Q(1,3),
                kernel_integrated=Q(2,3),reference_X=6)


def test_feedback_against_eliminated_two_by_two_system(budget):
    out=m.exterior_comparison(**inputs(),budget=budget)
    # X=1+2P, P=2+(X/2+2)/4 -> 3P/4=21/8.
    assert out==dict(loop=Q(1,4),P=Q(7,2),X=8,Q=15)


def test_zero_departure_forcing(budget):
    args=inputs()
    for key in ('source_X','source_Q','wall_P','kernel_pointwise','kernel_integrated'):
        args[key]=0
    out=m.exterior_comparison(**args,budget=budget)
    assert all(out[k]==0 for k in ('X','P','Q'))


@pytest.mark.parametrize('change,match',[
    ({'duration':1},'reaches one'),({'source_X':-1},'nonnegative'),
    ({'duration':.25},'exact rational')])
def test_invalid_comparison_stops(change,match,budget):
    args=inputs();args.update(change)
    with pytest.raises(ValueError,match=match):m.exterior_comparison(**args,budget=budget)


def test_reference_kernel_log_derivative(budget):
    budget.check(1024)
    x,a=s.symbols('x a',positive=True);r=s.Function('r')(x)
    C=x/r*s.exp((x-r)/a);K=C*a/(4*r**3)
    assert s.simplify(s.diff(K,x)/K-(1/x+1/a-s.diff(r,x)*(1/a+4/r)))==0


def test_arrival_inverse_second_derivative(budget):
    budget.check(1024)
    y=s.symbols('y');f=s.Function('f');q=s.Function('q')(y)
    w=f(q)
    recovered=(s.diff(w,y,2)*s.diff(q,y)-s.diff(w,y)*s.diff(q,y,2))/s.diff(q,y)**3
    assert s.simplify(recovered-s.diff(f(q),q,2))==0


def test_common_z_collar_chain_rule(budget):
    budget.check(1024)
    x,eta,zstar=s.symbols('x eta zstar');L=1+x+x*x;T=x+x*x/3
    z=zstar+eta*L
    # Exact reference wave F(T+z) with cubic F: mixed and transport
    # derivatives independently expose both clock and changing collar.
    chi=(T+z)**3;ep=6*(T+z)**2
    chi_speed=eta*s.diff(L,x)/L
    plus_speed=(eta*s.diff(L,x)+s.diff(T,x))/L
    assert s.simplify(s.diff(chi,x)-chi_speed*s.diff(chi,eta)-s.diff(T,x)*ep/2)==0
    assert s.simplify(s.diff(ep,x)-plus_speed*s.diff(ep,eta))==0


def test_mixed_equation_primitives_keep_both_variables(budget):
    budget.check(1024)
    u,x,U=s.symbols('u x U');X=1+x+u*x+u*u*x*x
    P=s.diff(X,u);Qx=s.diff(X,x);forcing=-s.diff(X,u,x)
    assert s.simplify(Qx-Qx.subs(u,0)+s.integrate(forcing,(u,0,u)))==0
    assert s.simplify(X-X.subs(u,0)-s.integrate(P,(u,0,u)))==0
    assert s.simplify(P-P.subs(x,0)+s.integrate(forcing,(x,0,x)))==0
    assert s.simplify(P+forcing*X)!=0  # a first-order replacement fails


def test_budget_stops_before_certificate_access(budget):
    class Spent:
        def check(self,count=0):raise RuntimeError('spent')
    with pytest.raises(RuntimeError,match='spent'):
        m.scalar_departure_bound({}, {}, {}, {}, {},budget=Spent())
