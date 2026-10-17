"""Exact controls for integrated exterior feedback and cap metric variation."""
import pytest
import sympy as s
import mpmath as mp
import numpy as np

from control_budget import ControlBudget
from recursive_horizons.nsc_atlas_nonlinear_bulk import bulk_operator


CONTROLS=ControlBudget(maximum_evaluations=30000,cpu_seconds=15.)


@pytest.fixture
def budget():
    with CONTROLS.active() as value:
        yield value


@pytest.fixture(scope='module',autouse=True)
def report():
    yield
    print('integrated metric feedback controls:',CONTROLS.budget.evaluations,
          'charged values;',CONTROLS.active_cpu_seconds,'active CPU seconds')


def test_cap_constraint_difference_keeps_its_metric_coefficient(budget):
    budget.check(8000)
    r,alpha,A0,N2,p1,p2,h1,h2,P1,P2=s.symbols('r alpha A0 N2 p1 p2 h1 h2 P1 P2')
    a1=alpha*(p1*p1+h1*h1)/(2*r); a2=alpha*(p2*p2+h2*h2)/(2*r)
    b1=r*r*P1/2+A0*(p1*p1+h1*h1)/4
    b2=r*r*P2/2+A0*(p2*p2+h2*h2)/4
    A2=A0-2*alpha*N2/r
    target=r*r*(P1-P2)/2+A2*((p1*p1+h1*h1)-(p2*p2+h2*h2))/4
    assert s.expand(b1-b2-N2*(a1-a2)-target)==0
    # The wall metric input is nonzero; the integration norm supplies its delay.
    H,U,S=s.symbols('H U S',nonzero=True)
    wall=(S*r-A2/H)/(2*U)
    assert s.simplify(s.diff(wall,N2)-alpha/(U*H*r))==0
    pr,vr,C=s.symbols('pr vr C')
    mixed=-(pr*vr+C/4)/r
    assert s.expand(r*mixed+pr*vr+C/4)==0


def test_moving_interval_retains_the_whole_boundary_history(budget):
    budget.check(8000)
    x,tau,W0,v,F=s.symbols('x tau W0 v F',positive=True)
    boundary=tau**3
    # Old points start with zero jet. Born points inherit b(tau), then p_x=F.
    integral=W0*F*x+s.integrate(v*(boundary+F*(x-tau)),(tau,0,x))
    expected=W0*F*x+v*x**4/4+v*F*x*x/2
    assert s.expand(integral-expected)==0
    assert s.expand(s.diff(integral,x)-(v*x**3+F*(W0+v*x)))==0
    # A pointwise current wall value is neither the spatial history nor zero.
    born=integral.subs(F,0)
    assert born!=0
    assert s.expand(born-(W0+v*x)*x**3)!=0
    assert s.limit(born/x**3,x,0)==0


def test_coupled_instantaneous_graph_and_time_feedback_are_distinct(budget):
    budget.check(10000)
    a,b,c=s.symbols('a b c',nonnegative=True)
    # Order: Jr, JL, JX, clock, incoming, outgoing, cap value.
    B=s.zeros(7); B[4,0]=a; B[4,1]=b; B[4,2]=c; B[6,2]=1
    assert B**2==s.zeros(7)
    assert (s.eye(7)-B)*(s.eye(7)+B)==s.eye(7)
    B=B.subs({a:2,b:3,c:1})
    A=s.zeros(7)
    A[0,4]=4; A[1,5]=2; A[2,4]=3; A[3,6]=1
    A[4,5]=2; A[5,6]=3; A[6,4]=4
    assert (A*B)[0,0]!=0  # Cap-to-metric feedback is kept in the time block.
    G=s.eye(7)+B; K=G*A
    duration=s.Rational(1,100)
    rate=max(sum(K[i,j] for j in range(7)) for i in range(7))
    assert duration*rate<1
    d=s.Matrix([1,0,0,0,0,0,0])
    exact=(s.eye(7)-B-duration*A).inv()*d
    upper=max(G*d)/(1-duration*rate)
    assert all(0<=value<=upper for value in exact)
    assert exact[0]>1  # Omitting the time feedback would give the wrong result.


def test_canonical_lapse_source_has_no_extra_conformal_metric_factor(budget):
    budget.check(3000)
    with mp.workdps(70):
        r,C,alpha,f,pr,vr,P,Q=map(mp.mpf,('.5','.7','.4','.2','.25','1.1','.3','.4'))
        a=mp.sqrt(alpha/6); u=mp.sinh(a*f)/a; root=mp.cosh(a*f)
        fU=(P-f*pr)/r; fx=(Q-f*vr)/r
        c=r*u; cU=u*pr+r*root*fU; cx=u*vr+r*root*fx
        expected=2*(pr*vr+C/4)/(r*r)-alpha*fU*fx
        wrong=2*(pr*vr+C/4)/(r*r)-alpha*C*fU*fx
        assert abs(expected-wrong)>mp.mpf('1e-3')
        state=dict(r=float(r),C=float(C),c=float(c))
        jets=dict(r_U=float(pr),r_V=float(vr),c_U=float(cU),c_V=float(cx))
    actual=bulk_operator(state,jets,alpha=float(alpha),region='exterior')
    np.testing.assert_allclose(actual['log_C_UV'],float(expected),rtol=2e-15)
