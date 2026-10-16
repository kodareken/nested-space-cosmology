"""Exact reconstruction and formula-switch identities; no field sampling."""
import pytest
import sympy as s

from control_budget import ControlBudget

CONTROLS = ControlBudget(maximum_evaluations=1, cpu_seconds=20.)


@pytest.fixture(autouse=True)
def exact_budget():
    with CONTROLS.active():
        yield


@pytest.fixture(scope="module", autouse=True)
def report():
    yield
    print("cap interface algebra:", CONTROLS.budget.evaluations,
          "field samples;", CONTROLS.active_cpu_seconds, "active CPU seconds")


def test_scalar_primitive_matches_cell_moments_including_the_regular_first_cell():
    a,b,c = s.symbols("a b c", positive=True)
    chi,iq0,iq1,slope = s.symbols("chi iq0 iq1 slope")
    # An ordinary cell has arbitrary limited slope, but its integral is
    # exactly its cell moment. Thus either side gives the same chi/r.
    qc = iq0/(b-a)-slope*(a+b)/2
    end = chi+(b-a)*(qc+slope*(b+a)/2)
    assert s.expand(end-chi-iq0) == 0
    # The first cell instead fits an even quadratic to its first two
    # cell means. Its endpoint primitive still equals Iq[0].
    qm0,qm1 = b*b/3,(b*b+b*c+c*c)/3
    q2 = (iq1/(c-b)-iq0/b)/(qm1-qm0)
    q0 = iq0/b-q2*qm0
    first = b*q0+b**3*q2/3
    assert s.simplify(first-iq0) == 0
    assert s.simplify(first/b-iq0/b) == 0
    # A constant first-cell q0 alone would lose this continuity.
    assert s.simplify(b*q0-iq0) != 0


def test_potential_formula_switch_has_a_nonzero_value_and_derivative_jump():
    q,f = s.symbols("q f", positive=True)
    exact = s.sinh(q)/q
    truncated = 1+q*q/6+q**4/120
    factor = (s.cosh(q)**2+1)/2
    difference = factor*(exact**2-truncated**2)
    # P_exact-P_small = f^2*difference, with q proportional to f.
    # The derivative difference is f*(2*difference+q*difference').
    assert s.series(difference,q,0,8).removeO().expand() == q**6/2520
    derivative = 2*difference+q*s.diff(difference,q)
    assert s.series(derivative,q,0,8).removeO().expand() == q**6/315
    # For a simple level-set crossing, retaining only the bulk derivative
    # misses its moving jump. Check against an independently integrated mass.
    r,x,h,jump = s.symbols("r x h jump", positive=True)
    wall = s.Rational(3,4)
    level = h-x
    left = (r+x)**2
    right = left+jump
    mass = s.integrate(r*r*left/2,(r,0,level))+s.integrate(r*r*right/2,(r,level,wall))
    bulk = s.integrate(r*r*(r+x),(r,0,wall))
    interface = -s.diff(level,x)*level**2*jump/2
    assert s.expand(s.diff(mass,x)-bulk-interface) == 0
    assert s.expand(s.diff(mass,x)-bulk) != 0


def test_small_radius_formula_jump_has_its_declared_leading_order():
    r = s.symbols("r", positive=True)
    exact = (r/(1-r*r)-s.atanh(r))/(r*r)
    polynomial = sum(s.Rational(2*k,2*k+1)*r**(2*k-1) for k in range(1,6))
    difference = exact-polynomial
    assert s.series(difference,r,0,15).removeO().expand() == (
        s.Rational(12,13)*r**11+s.Rational(14,15)*r**13)
    # The background scalar and Pi agree between formulas; Phi has this
    # coefficient times the common first moment. Its mass/lapse jump cannot
    # be erased merely because the coefficient is small.
    alpha,pi,phi,j,n = s.symbols("alpha pi phi j n")
    a_left = alpha*r*(pi*pi+phi*phi)/2
    a_right = alpha*r*(pi*pi+(phi+j)**2)/2
    b_jump = r*r*(1-r*r)*((phi+j)**2-phi*phi)/4
    rho_jump = b_jump-(a_right-a_left)*n
    A = 1-r*r-2*alpha*n/r
    assert s.expand(rho_jump-r*r*A*(2*phi*j+j*j)/4) == 0
