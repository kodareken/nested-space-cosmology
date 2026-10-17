"""Canonical equation/port checks and continuum comparison controls."""
from fractions import Fraction as Q

import mpmath as mp
import numpy as np
import pytest
import sympy as s

from control_budget import ControlBudget
from recursive_horizons.nsc_atlas_nonlinear_boundary import hybrid_wall_boundary
from recursive_horizons.nsc_atlas_transmission_bound import (
    COMPONENTS, COEFFICIENTS, scalar_inverse_bound, canonical_source_error_bound)


CONTROLS = ControlBudget(maximum_evaluations=100000, cpu_seconds=30.)


@pytest.fixture
def budget():
    with CONTROLS.active() as value:
        yield value


@pytest.fixture(scope='module', autouse=True)
def report():
    yield
    print('connected scalar inverse controls:', CONTROLS.budget.evaluations,
          'charged values;', CONTROLS.active_cpu_seconds, 'active CPU seconds')


def coefficients(**changes):
    return dict.fromkeys(COEFFICIENTS, Q(0)) | changes


def data(**changes):
    return dict.fromkeys(COMPONENTS, Q(0)) | changes


def test_canonical_exterior_equation_and_nonlinear_wall_derivative(budget):
    budget.check(15000)
    r, rU, rV, C, f, fU, fV = s.symbols('r rU rV C f fU fV', nonzero=True)
    rUV = -(rU*rV+C/4)/r
    fUV = -(rU*fV+rV*fU)/r
    XUV = rUV*f+rU*fV+rV*fU+r*fUV
    assert s.simplify(XUV+(rU*rV+C/4)/r**2*(r*f)) == 0
    y, a, R = s.symbols('y a R', positive=True)
    J = 2*R*s.cosh(y)**3*(y-s.tanh(y))/a
    derivative = s.diff(J, y)*a/R
    expected = 2*s.cosh(y)*s.sinh(y)*(3*y*s.cosh(y)-2*s.sinh(y))
    assert s.simplify(s.expand_trig(derivative-expected)) == 0
    assert s.limit(J, y, 0) == 0
    assert s.limit(expected/y**2, y, 0) == 2


def test_canonical_port_agrees_with_existing_invariant_hybrid(budget):
    budget.check(7000)
    # A finite, deliberately nonzero nonlinear jump; no physical input or solve.
    for alpha in (0., .7):
        R, f, em, q, A, C, eta, rv = .48, .6, -.013, .021, .7, .9, .03, -.11
        with mp.workdps(75):
            rr, ff, al = map(mp.mpf, (R, f, alpha))
            aa = mp.sqrt(al/6)
            uu = mp.sinh(aa*ff)/aa if aa else ff
            root = mp.cosh(aa*ff)
            jj = ff-mp.tanh(aa*ff)/aa if aa else mp.mpf(0)
            U = mp.exp(-mp.mpf(eta))/mp.sqrt(mp.mpf(C))
            V = mp.exp(mp.mpf(eta))/mp.sqrt(mp.mpf(C))
            H = 2*V*mp.mpf(rv)+2*root**3*rr
            J = 2*root**3*rr*jj
            A0 = 1-rr*rr
            c = rr*uu
            cv = root*(mp.mpf(q)-jj*mp.mpf(rv))
            wm = (mp.mpf(em)/A0+ff)/rr
            expected_ep = A0/H*(2*V*mp.mpf(q)+J)
            expected_p = (mp.mpf(A)*mp.mpf(em)/(H*A0)+J)/(2*U)
        out = hybrid_wall_boundary(dict(r=R,c=float(c),eta_plus=eta),
            dict(r_V=rv,c_V=float(cv),log_C_V=.04),
            dict(A=A,delta=.02,w_minus=float(wm)),dict(C_wall=C,log_C_s=.01),alpha=alpha)
        ep = (1-R*R)*(R*float(out['w_plus'])+f)
        canon_p = float(out['c_U_plus'])/float(root)+float(jj)*float(out['r_U_plus'])
        np.testing.assert_allclose([ep,canon_p],[float(expected_ep),float(expected_p)],
                                   rtol=2e-14,atol=2e-16)
        canon_s = float(out['c_s'])/float(root)+float(jj)*float(out['r_s'])
        nplus = float(out['n_c'])/float(root)+float(jj)*float(out['beta_plus'])
        nminus = float(out['n_c'])/float(root)+float(jj)*float(out['beta_minus'])
        np.testing.assert_allclose(nplus-nminus,-float(J),rtol=2e-13,atol=1e-16)
        np.testing.assert_allclose(canon_s/float(V),q+(float(U)/float(V))*canon_p,rtol=2e-14)
        if alpha:
            assert abs(nplus-nminus) > 1e-3  # Omitting the wall term must fail.


def test_large_cross_gains_need_no_instantaneous_small_gain_assumption(budget):
    c = coefficients(width=2,exterior_potential=3,cap_reaction=1,cap_drift=1,
        cap_speed=3,wall_speed=2,incoming_gain=7,outgoing_gain=11,
        incoming_value_gain=2,outgoing_value_gain=3)
    d = data(cap_outgoing=1,wall_value=Q(1,100))
    out = scalar_inverse_bound(c,d,duration=Q(1,10000),budget=budget)
    budget.check(10000)
    B,A,G = (s.Matrix(out[key]) for key in ('instantaneous','evolution','inverse'))
    assert B**5 == s.zeros(7)
    assert B**4 != s.zeros(7)
    assert (s.eye(7)-B)*G == s.eye(7)
    assert G == (s.eye(7)-B).inv()
    # An independently solved saturated sup inequality must lie inside the bound.
    saturated = (s.eye(7)-B-s.Rational(1,10000)*A).inv()*s.Matrix([d[k] for k in COMPONENTS])
    assert all(0 <= saturated[i] <= out['bounds'][key] for i,key in enumerate(COMPONENTS))
    assert out['theta'] < 1
    assert c['incoming_gain']*c['outgoing_gain'] > 1


def test_exact_incoming_quartic_on_expanding_cap_obeys_connected_inequalities(budget):
    budget.check(13000)
    t = s.symbols('t',nonnegative=True)
    z = s.symbols('z',real=True)
    speed = s.Rational(1,5)
    # Active branch; its positive-part extension is C3 with zero quiet data.
    field = (t+z-1)**4/(1+speed)**4
    ep = 2*s.diff(field,t)
    assert s.diff(ep,t)-s.diff(ep,z) == 0
    assert s.diff(field,t)-ep/2 == 0
    assert s.simplify(field.subs(z,1+speed*t)-t**4) == 0
    a = 2/(1+speed)
    assert s.simplify(ep.subs(z,1+speed*t)-a*4*t**3) == 0
    c = coefficients(width=2,cap_speed=1,wall_speed=Q(3,10),
                     incoming_gain=Q(a),outgoing_gain=3)
    duration = Q(1,20)
    d = data(exterior_value=duration**4,exterior_longitudinal=4*duration**3)
    out = scalar_inverse_bound(c,d,duration=duration,budget=budget)
    u = s.Matrix([0,t**4,0,t**4,4*t**3,a*4*t**3,t**4])
    dd = s.Matrix([0,0,0,t**4,4*t**3,0,0])
    B,A = s.Matrix(out['instantaneous']),s.Matrix(out['evolution'])
    gap = dd+B*u+(A*u).applyfunc(lambda v:s.integrate(v,(t,0,t)))-u
    assert all(s.simplify(v).is_nonnegative for v in gap)
    assert all(u[i].subs(t,s.Rational(duration)) <= out['bounds'][key]
               for i,key in enumerate(COMPONENTS))


def test_quiet_core_radial_difference_uses_damping_not_exponential_growth(budget):
    budget.check(6000)
    # Exact constant-coefficient constraints on a collar, independently solved.
    with mp.workdps(70):
        length = mp.mpf('.1'); a1,a2,b1,b2 = map(mp.mpf,('.3','.2','.07','.05'))
        n1 = b1*(-mp.expm1(-a1*length))/a1
        n2 = b2*(-mp.expm1(-a2*length))/a2
        bound = length*(abs(b1-b2)+n2*abs(a1-a2))
        assert abs(n1-n2) < bound
        assert abs((a1-a2)*length) == length*abs(a1-a2)
    n1,n2,a1,a2,b1,b2 = s.symbols('n1 n2 a1 a2 b1 b2')
    residual = (b1-a1*n1)-(b2-a2*n2)+a1*(n1-n2)-(b1-b2-n2*(a1-a2))
    assert s.expand(residual) == 0


def test_exact_input_guards_zero_forcing_weights_and_failed_margin(budget):
    c = coefficients(width=1,exterior_potential=1,cap_speed=2,incoming_gain=2,outgoing_gain=2)
    out = scalar_inverse_bound(c,data(),duration=Q(1,100),budget=budget,
                              weights=dict.fromkeys(COMPONENTS,Q(2)))
    assert set(out['bounds'].values()) == {0}
    with pytest.raises(ValueError,match='margin'):
        scalar_inverse_bound(c,data(),duration=10,budget=budget)
    with pytest.raises(ValueError,match='exact integer'):
        scalar_inverse_bound(c,data(),duration=.01,budget=budget)
    with pytest.raises(ValueError,match='positive domain'):
        scalar_inverse_bound(c,data(),duration=-1,budget=budget)


def test_canonical_source_jacobian_and_finite_difference_enclosure(budget):
    budget.check(10000)
    r,c,rx,cx,a = s.symbols('r c rx cx a',positive=True)
    u=c/r; h=1/s.sqrt(1+a*a*u*u)
    field=r*s.asinh(a*u)/a; j=s.asinh(a*u)/a-u*h
    first=h*cx+j*rx
    derivative=-a*a*u*h**3*(cx-u*rx)/r
    assert s.simplify(s.diff(field,c)-h)==0
    assert s.simplify(s.diff(field,r)-j)==0
    for variable,expected in ((c,derivative),(r,-u*derivative),(cx,h),(rx,j)):
        assert s.simplify(s.diff(first,variable)-expected)==0
    errors=dict(r=Q(1,10000),c=Q(1,20000),r_x=Q(1,1000),c_x=Q(1,2000))
    args=dict(alpha=Q(3,5),radius_floor=Q(2,5),u_bound=Q(1,2),c_x_bound=1,
              r_x_bound=2,errors=errors,budget=budget)
    bound=canonical_source_error_bound(**args)
    with mp.workdps(75):
        def convert(row):
            rr,cc,rrx,ccx=row; aa=mp.sqrt(mp.mpf(3)/30)
            uu=cc/rr; hh=1/mp.sqrt(1+aa*aa*uu*uu)
            ff=mp.asinh(aa*uu)/aa
            return rr*ff, hh*ccx+(ff-uu*hh)*rrx
        left=list(map(mp.mpf,('.45','.1','1.2','.3')))
        right=[v+mp.mpf(e.numerator)/e.denominator for v,e in zip(left,errors.values())]
        actual=[abs(v-w) for v,w in zip(convert(left),convert(right))]
        assert all(value <= mp.mpf(bound[key].numerator)/bound[key].denominator
                   for value,key in zip(actual,('value','first')))
    zero=canonical_source_error_bound(**(args|dict(alpha=0)))
    assert zero['value']==errors['c'] and zero['first']==errors['c_x']
