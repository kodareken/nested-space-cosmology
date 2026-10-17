"""Frame-identity and finite-box controls without original-data evaluation."""
from fractions import Fraction as Q

from flint import arb,ctx
import pytest
import sympy as s

from control_budget import ControlBudget
from recursive_horizons import nsc_atlas_jordan_radius_bound as m
from recursive_horizons.nsc_atlas_geometric_remainder_bound import Jet,even_product,jet_box

CONTROLS=ControlBudget(maximum_evaluations=16384,cpu_seconds=20.)


@pytest.fixture
def budget():
    with CONTROLS.active() as value:
        with ctx.workprec(160):yield value


@pytest.fixture(scope='module',autouse=True)
def report():
    yield
    print('Jordan radius controls:',CONTROLS.budget.evaluations,'charged values;',CONTROLS.active_cpu_seconds,'active CPU seconds')


def test_canonical_invariant_frame_identity(budget):
    budget.check(512)
    q=s.symbols('q',real=True)
    assert s.simplify(s.cosh(q)**2-1-s.sinh(q)**2)==0
    # c=R*sinh(a*f)/a, hence R_J^2=R^2+alpha*c^2/6.
    # Using c/R as canonical f would reverse this quartic coefficient.
    t=s.symbols('t')
    assert s.diff(1+t/2+t*t/24,t,2)==s.Rational(1,12)
    assert s.diff(s.sqrt(1+t),t,2).subs(t,0)==s.Rational(-1,4)


def test_symbolic_zero_amplitude_jet(budget):
    budget.check(2048)
    a,R,X,cr,cx=s.symbols('a R X cr cx',real=True)
    r=R+a*cr;x=X+a*cx;q=a*(x/r)**2/6
    frame=r*(1+q/2+q*q/24)
    assert s.simplify(s.diff(frame,a).subs(a,0)-cr-X**2/(12*R))==0
    expected=X*cx/(6*R)-X**2*cr/(12*R**2)+X**4/(864*R**3)
    assert s.simplify(s.diff(frame,a,2).subs(a,0)/2-expected)==0
    actual=Jet(2,Q(1,7));field=Jet(Q(1,3),Q(-1,11))
    _,ch,_,_=even_product(Q(0),field/actual)
    exact=Q(str(expected.subs({R:2,X:s.Rational(1,3),cr:s.Rational(1,7),cx:s.Rational(-1,11)})))
    value=(actual*ch).second/2
    assert m.dyadic(value.lower())<=exact<=m.dyadic(value.upper())


def test_positive_bound_contains_finite_alpha_box_hessian(budget):
    args=dict(radius_lower=Q(2),radius_upper=Q(3),X_upper=Q(1,3),
              radius_difference=Q(1,7),X_difference=Q(1,11),alpha_max=Q(1,10))
    bound=m.frame_half_hessian(**args,budget=budget)
    # Interior value boxes leave room for Arb's outward radius rounding;
    # both first-direction boxes still cross zero.
    R=jet_box(Q(9,4),Q(11,4),Q(1,7));X=jet_box(Q(-1,4),Q(1,4),Q(1,11))
    _,ch,_,_=even_product(Q(1,10),X/R)
    actual=(R*ch).second/2
    assert abs(actual).upper()<=m._ball(bound['half_hessian'])
    assert 0<bound['conformal_factor_lower']<=1


def test_zero_source_and_first_difference_limit(budget):
    out=m.frame_half_hessian(radius_lower=1,radius_upper=2,X_upper=0,
        radius_difference=Q(1,7),X_difference=0,alpha_max=1,budget=budget)
    assert out['half_hessian']==0 and out['conformal_factor_lower']==1


def test_signed_event_shift_and_invalid_inputs(budget):
    budget.check(1024)
    first=dict(constants={k:dict(lower=1,upper=2) for k in ('a','v','d')},
        outputs=dict(eta=dict(lower=-3,upper=1),rho_proper_time=dict(lower=-2,upper=-1)))
    assert m.first_jordan_interval(first)==dict(lower=-50,upper=15)
    first['constants']['d']['lower']=0
    with pytest.raises(ValueError,match='positive reference'):
        m.first_jordan_interval(first)
    class Spent:
        def check(self,count=0):raise RuntimeError('spent')
    with pytest.raises(RuntimeError,match='spent'):
        m.jordan_radius_bound({}, {}, {}, {}, {},budget=Spent())


def test_frame_rejects_lost_radius_floor(budget):
    with pytest.raises(ValueError,match='positive radius'):
        m.frame_half_hessian(radius_lower=0,radius_upper=1,X_upper=1,
            radius_difference=0,X_difference=0,alpha_max=1,budget=budget)
