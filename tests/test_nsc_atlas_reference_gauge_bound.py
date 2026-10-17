"""Independent identities behind the exact reference interval certificate."""
from fractions import Fraction as Q

from flint import acb, ctx
import mpmath as mp
import pytest
import sympy as s

from control_budget import ControlBudget
from recursive_horizons.nsc_atlas_reference_gauge_bound import _advanced, _ball, _record

CONTROLS=ControlBudget(maximum_evaluations=20000,cpu_seconds=15.)


@pytest.fixture
def budget():
    with CONTROLS.active() as value:
        yield value


@pytest.fixture(scope='module',autouse=True)
def report():
    yield
    print('reference gauge controls:',CONTROLS.budget.evaluations,
          'charged values;',CONTROLS.active_cpu_seconds,'active CPU seconds')


def test_junction_rationalization_and_wall_normalization(budget):
    budget.check(10000)
    R,x,v=s.symbols('R x v',positive=True); a=s.Rational(4,5)
    bm=1/(5*R**2)+3*R/4; bp=1/(5*R**2)-5*R/4
    q=25*R**2/16-1+3/(10*R)+1/(25*R**4)
    assert s.factor(q-bp**2-(a-R)/R)==0
    assert s.factor(bm**2-q-(1-R**2))==0
    def reduce(expression):
        numerator=s.together(expression).as_numer_denom()[0]
        return s.factor(s.rem(numerator,v*v-q,v))
    assert reduce(-1/(v*(v-bp))+R/(a-R)*(1+bp/v))==0
    Rx=2*x*v*(v-bp)/(a-x); sx=Rx/v
    # Differentiate W=4a(1-F(R)/F(x)) along the implicit inverse;
    # C*F(R)/F(x)=x(a-R)/(R(a-x)), so exponentials cancel exactly.
    CWx=4*x*(a-R)/(R*(a-x))*(R*Rx/(a-R)-x/(a-x))
    assert reduce(CWx-sx*sx)==0
    assert reduce(bm*sx/(1-R*R)-Rx/(1-R*R)-sx/(bm+v))==0


def test_rapidity_formula_matches_original_second_jet(budget):
    budget.check(5000)
    with mp.workdps(70):
        x,R=mp.mpf('0.43'),mp.mpf('0.49'); a=mp.mpf(4)/5
        bp=1/(5*R**2)-5*R/4
        v=mp.sqrt(25*R**2/16-1+3/(10*R)+1/(25*R**4))
        vR=(25*R/8-3/(10*R**2)-4/(25*R**5))/(2*v)
        bpR=-2/(5*R**3)-mp.mpf(5)/4
        lam=-x/(a*(a-x)); lamx=-1/(a-x)**2
        advance=-1/(v*(v-bp)); logadvance=-vR/v-(vR-bpR)/(v-bp)
        Rx=2*a*lam/advance; Rxx=Rx*(lamx/lam-logadvance*Rx)
        p=-R/(a*(a-R)); px=-1/(a-R)**2
        Z=p*Rx-lam; Zx=px*Rx**2+p*Rxx-lamx
        original=-(Z*Z+Zx)/(2*Z)
        candidate=-(1/x+1/(a-x)+(vR-bpR)*Rx/(v-bp))
        candidate+=(1/x+1/a-(1/R+1/a)*Rx)/2
        assert abs(original-candidate)<mp.mpf('1e-65')


def test_complex_integrand_and_signed_dyadic_records(budget):
    budget.check(1000)
    with ctx.workprec(160):
        R=acb(_ball(Q(1,2)))
        q=25*R**2/16-1+3/(10*R)+1/(25*R**4)
        bp=1/(5*R**2)-5*R/4
        direct=-1/(q.sqrt(analytic=True)*(q.sqrt(analytic=True)-bp))
        assert (_advanced(R,True)-direct).contains(0)
        row=_record(_advanced(R,True).real)
        def exact(key):
            v=row[key]; return Q(int(v['mantissa']))*Q(2)**v['exponent']
        assert Q(row['lower_float'])<=exact('lower')<=exact('upper')<=Q(row['upper_float'])<0
        with pytest.raises(ValueError,match='exact rational'):
            _ball(.5)
