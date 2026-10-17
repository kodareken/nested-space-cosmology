"""Source Taylor identities and homogeneous controls, without source replay."""
from fractions import Fraction as Q

import pytest
import sympy as s

from control_budget import ControlBudget
from recursive_horizons import nsc_atlas_source_remainder_bound as m

CONTROLS=ControlBudget(maximum_evaluations=12000,cpu_seconds=20.)


@pytest.fixture
def budget():
    with CONTROLS.active() as value:yield value


@pytest.fixture(scope='module',autouse=True)
def report():
    yield
    print('source remainder controls:',CONTROLS.budget.evaluations,
          'charged values;',CONTROLS.active_cpu_seconds,'active CPU seconds')


def inputs():
    return dict(g=Q(1,8),g_x=Q(1,4),g_xx=Q(1,2),h=Q(1,32),h_x=Q(1,16),
                duration=Q(1,16),right=Q(1,2),r_min=Q(1,4),F_min=Q(3,4))


def test_quiet_source(budget):
    v=inputs()
    for key in ('g','g_x','g_xx','h','h_x'):v[key]=0
    result=m.source_remainder_comparison(**v,budget=budget)
    assert all(value==0 for value in result['coefficients'].values())


def test_quartic_homogeneity_at_fixed_tube(budget):
    v=inputs();k=Q(3,2)
    original=m.source_remainder_comparison(**v,budget=budget)
    for key in ('g','g_x','g_xx'):v[key]*=k
    for key in ('h','h_x'):v[key]*=k*k
    scaled=m.source_remainder_comparison(**v,budget=budget)
    assert scaled['coefficients']=={key:k**4*value for key,value in original['coefficients'].items()}


@pytest.mark.parametrize('change,match',[
    ({'r_min':0},'positive Jordan'),({'F_min':Q(1,2)},'positive Jordan'),
    ({'g':-1},'nonnegative'),({'h':.1},'exact integer')])
def test_domain_guards(change,match,budget):
    v=inputs();v.update(change)
    with pytest.raises(ValueError,match=match):m.source_remainder_comparison(**v,budget=budget)


def test_frame_taylor_coefficients(budget):
    budget.check(1024)
    al,x,g,H,K=s.symbols('al x g H K',positive=True)
    rj=x+al*H+al*al*K
    r=s.sqrt(rj*rj-al*g*g/6)
    logF=s.log(1-al*g*g/(6*rj*rj))
    assert s.simplify(s.diff(r,al).subs(al,0)-(H-g*g/(12*x)))==0
    assert s.simplify(s.diff(r,al,2).subs(al,0)/2-(K+g*g*H/(12*x*x)-g**4/(288*x**3)))==0
    assert s.simplify(s.diff(logF,al).subs(al,0)+g*g/(6*x*x))==0
    assert s.simplify(s.diff(logF,al,2).subs(al,0)/2-(g*g*H/(3*x**3)-g**4/(72*x**4)))==0


def test_product_and_lapse_slope_identities(budget):
    budget.check(1024)
    x,al=s.symbols('x al');rj=s.Function('rj')(x);g=s.Function('g')(x)
    rad2=rj*rj-al*g*g/6;F=rad2/(rj*rj)
    assert s.simplify(s.diff(rad2,x)/2-(rj*s.diff(rj,x)-al*g*s.diff(g,x)/6))==0
    expected=-al*g/(3*F)*(s.diff(g,x)/rj**2-g*s.diff(rj,x)/rj**3)
    assert s.simplify(s.diff(F,x)/F-expected)==0


def test_radius_quotient_derivative(budget):
    budget.check(1024)
    u=s.symbols('u');root=s.sqrt(1-u)
    f=1/(1+root)
    assert s.simplify(s.diff(f,u)-1/(2*root*(1+root)**2))==0
    assert s.simplify(root-(1-u*f))==0


def test_budget_stops_before_source_access(budget):
    class Spent:
        def check(self,count=0):raise RuntimeError('spent')
    with pytest.raises(RuntimeError,match='spent'):
        m.source_remainder_bound({}, {},budget=Spent())
