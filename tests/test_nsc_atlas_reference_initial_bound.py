"""Initial inverse, quiet-history and common-coordinate controls."""
from fractions import Fraction as Q

from flint import ctx
import numpy as np
import pytest
import sympy as s

from control_budget import ControlBudget
from recursive_horizons.nsc_atlas_characteristic_background import CharacteristicProfile
from recursive_horizons.nsc_atlas_reference_gauge_bound import _ball
from recursive_horizons.nsc_atlas_reference_initial_bound import inverse_seed_enclosure

CONTROLS=ControlBudget(maximum_evaluations=20000,cpu_seconds=15.)


@pytest.fixture
def budget():
    with CONTROLS.active() as value:
        yield value


@pytest.fixture(scope='module',autouse=True)
def report():
    yield
    print('reference initial controls:',CONTROLS.budget.evaluations,
          'charged values;',CONTROLS.active_cpu_seconds,'active CPU seconds')


def test_inverse_residual_encloses_a_known_nonlinear_root(budget):
    budget.check(1000)
    with ctx.workprec(160):
        seed,root=_ball(Q(1,3)),_ball(Q(1,4))
        residual=(-seed-seed**2)-(-root-root**2)
        enclosure,error=inverse_seed_enclosure(seed,residual,_ball(Q(7,5)))
        assert enclosure.contains(root)
        assert error>abs(root-seed)
        with pytest.raises(ValueError,match='positive derivative'):
            inverse_seed_enclosure(seed,residual,_ball(0))


def test_zero_scalar_closes_both_wall_curvatures_and_rates(budget):
    budget.check(8000)
    R,x,v=s.symbols('R x v',positive=True);a=s.Rational(4,5)
    bm=1/(5*R**2)+3*R/4;bp=1/(5*R**2)-5*R/4
    q=25*R**2/16-1+3/(10*R)+1/(25*R**4)
    assert s.simplify(bm-bp-2*R)==0
    assert s.simplify((s.diff(q,R)/2-R)/bm-s.diff(bm,R))==0
    assert s.simplify((s.diff(q,R)/2+a/(2*R**2))/bp-s.diff(bp,R))==0
    assert s.simplify(s.diff(bm-bp,R)-2)==0
    sx=2*x*(v-bp)/(a-x)
    rx=x*(a-R)/(R*(a-x))
    def reduce(expression):
        return s.factor(s.rem(s.together(expression).as_numer_denom()[0],v*v-q,v))
    assert reduce(2*rx/sx+2*R-(v+bm))==0  # H at zero field
    assert reduce((1-R**2)/(v+bm)-(bm-v))==0
    # The wall's two exterior null derivatives give the same proper speed.
    assert reduce(sx*(a-x)/(4*x)+rx/sx-v)==0


def test_quiet_prefix_is_exact_and_rejects_a_nonzero_input(budget):
    labels=np.array([-1.,0.,.5,1.]);values=np.array([0.,0.,0.,1.])
    first=np.array([0.,0.,0.,2.]);second=np.array([0.,0.,0.,3.])
    profile=CharacteristicProfile(labels,values,first,second,
                                  onset=dict(index=2,quartic_coefficient=1.))
    budget.check(1000)
    jets=profile.jet(np.array([-2.,-1.,0.,.5]),budget=budget)
    assert all(np.all(jets[k]==0) for k in ('F','F_prime','F_second'))
    bad=values.copy();bad[1]=1e-100
    with pytest.raises(ValueError,match='exactly quiet'):
        CharacteristicProfile(labels,bad,first,second,onset=dict(index=2,quartic_coefficient=1.))


def test_common_reference_pullback_retains_both_moving_terms(budget):
    budget.check(5000)
    x,u,eta=s.symbols('x u eta',positive=True)
    W=1+x;F=x*x*eta**3+s.sin(x)*eta
    physical=F.subs(eta,u/W)
    fx=(s.diff(F,x)-eta*s.diff(W,x)/W*s.diff(F,eta)).subs(eta,u/W)
    mixed=((s.diff(F,x,eta)-s.diff(W,x)/W*(s.diff(F,eta)+eta*s.diff(F,eta,2)))/W).subs(eta,u/W)
    assert s.simplify(fx-s.diff(physical,x))==0
    assert s.simplify(mixed-s.diff(physical,x,u))==0
