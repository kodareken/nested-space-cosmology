"""Exact onset order and coefficient before transferring the prepared history."""
import sympy as s
import pytest

from control_budget import ControlBudget

CONTROLS=ControlBudget(maximum_evaluations=1,cpu_seconds=20.)


@pytest.fixture(autouse=True)
def budget():
    with CONTROLS.active(): yield


@pytest.fixture(scope='module',autouse=True)
def report():
    yield
    print('history onset controls:',CONTROLS.budget.evaluations,
          'field samples;',CONTROLS.active_cpu_seconds,'active CPU seconds')


def test_quintic_and_sextic_do_not_preserve_both_order_and_coefficient():
    u=s.symbols('u');y,v,w=s.symbols('y v w')
    a,b,c=s.symbols('a b c')
    quintic=a*u**3+b*u**4+c*u**5
    eq=[quintic.subs(u,1)-y,s.diff(quintic,u).subs(u,1)-v,s.diff(quintic,u,2).subs(u,1)-w]
    solution=s.solve(eq,(a,b,c))
    assert s.expand(solution[a]-(10*y-4*v+w/2))==0
    A,B,h=s.symbols('A B h',nonzero=True)
    x=s.symbols('x');exact=A*x**4+B*x**6
    jets={y:exact.subs(x,h),v:h*s.diff(exact,x).subs(x,h),w:h*h*s.diff(exact,x,2).subs(x,h)}
    assert s.expand(solution[a].subs(jets))==B*h**6
    sextic=a*u**4+b*u**5+c*u**6
    eq=[sextic.subs(u,1)-y,s.diff(sextic,u).subs(u,1)-v,s.diff(sextic,u,2).subs(u,1)-w]
    solution=s.solve(eq,(a,b,c))
    exact=A*x**4+B*x**7
    jets={y:exact.subs(x,h),v:h*s.diff(exact,x).subs(x,h),w:h*h*s.diff(exact,x,2).subs(x,h)}
    assert s.expand(solution[a].subs(jets))==A*h**4+B*h**7


def test_seventh_degree_first_piece_preserves_given_quartic_and_next_jets():
    u=s.symbols('u');A,h,y,v,w=s.symbols('A h y v w',nonzero=True)
    C4=A*h**4;Y=y-C4;V=v-4*C4;W=w-12*C4
    expected=(C4,21*Y-6*V+W/2,11*V-W-35*Y,(W-10*V+30*Y)/2)
    coefficients=s.symbols('c0:8');poly=sum(c*u**i for i,c in enumerate(coefficients))
    eq=[s.diff(poly,u,i).subs(u,0) for i in range(4)]
    eq += [s.diff(poly,u,4).subs(u,0)-24*C4,poly.subs(u,1)-y,
           s.diff(poly,u).subs(u,1)-v,s.diff(poly,u,2).subs(u,1)-w]
    matrix,rhs=s.linear_eq_to_matrix(eq,coefficients)
    assert matrix.det()!=0
    solved=matrix.inv()*rhs
    assert all(s.expand(solved[i])==0 for i in range(4))
    assert all(s.expand(solved[i+4]-value)==0 for i,value in enumerate(expected))
    exact=sum(s.symbols('a4:8')[i-4]*(h*u)**i for i in range(4,8))
    jets={A:s.symbols('a4'),y:exact.subs(u,1),v:s.diff(exact,u).subs(u,1),w:s.diff(exact,u,2).subs(u,1)}
    result=sum(value*u**(i+4) for i,value in enumerate(expected)).subs(jets,simultaneous=True)
    assert s.expand(result-exact)==0


def test_compact_source_coefficient_under_arrival_clock():
    dv,qv=s.symbols('dv qv',positive=True)
    x=s.symbols('x');q2,q3,e5=s.symbols('q2 q3 e5')
    zeta=-1+4*x/dv
    wall=(1-zeta*zeta)**4+e5*x**5
    arrival=qv*x+q2*x*x+q3*x**3
    assert s.limit(wall/arrival**4,x,0)==16*(4/(dv*qv))**4


def test_quintic_midpoint_derivative_can_hide_the_leading_error():
    x=s.symbols('x');remainder=x**3*(x-1)**3
    assert all(s.diff(remainder,x,j).subs(x,edge)==0 for j in range(3) for edge in (0,1))
    assert s.diff(remainder,x).subs(x,s.Rational(1,2))==0
    assert s.diff(remainder,x).subs(x,s.Rational(1,4))!=0
    # x^6 minus this remainder is exactly its endpoint-jet quintic.
    assert s.degree(s.expand(x**6-remainder),x)==5


def test_interior_value_and_first_jet_enrichment_preserves_endpoint_jets():
    x=s.symbols('x');psi=x**3*(1-x)**3
    eta=s.Rational(3,7);a,b=s.symbols('a b')
    # Two extra degrees remove the full degree-six/seven remainder.
    error=psi*(a+b*x)
    dv=error.subs(x,eta);dx=s.diff(error,x).subs(x,eta)
    ca=dv/psi.subs(x,eta)
    cb=(dx-s.diff(psi,x).subs(x,eta)*ca)/psi.subs(x,eta)
    correction=s.expand(psi*(ca+cb*(x-eta)))
    assert s.expand(correction-error)==0
    assert all(s.diff(correction,x,j).subs(x,edge)==0 for j in range(3) for edge in (0,1))
    assert correction.subs(x,eta)==dv
    assert s.diff(correction,x).subs(x,eta)==dx
    dm,dp,h=s.symbols('dm dp h')
    midpoint=64*psi*(dm+h*dp*(x-s.Rational(1,2)))
    assert midpoint.subs(x,s.Rational(1,2))==dm
    assert s.diff(midpoint,x).subs(x,s.Rational(1,2))==h*dp
