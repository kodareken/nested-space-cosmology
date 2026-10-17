"""Exact source/profile and finite-slab stability controls, without field runs."""
import pytest
import sympy as s

from control_budget import ControlBudget


CONTROLS = ControlBudget(maximum_evaluations=0, cpu_seconds=20.)


@pytest.fixture
def budget():
    with CONTROLS.active() as value:
        yield value


@pytest.fixture(scope="module", autouse=True)
def report():
    yield
    print("profile stability controls:", CONTROLS.budget.evaluations,
          "physical samples;", CONTROLS.active_cpu_seconds, "active CPU seconds")


def test_source_equation_from_full_conformal_null_stress(budget):
    x, alpha = s.symbols("x alpha", real=True)
    r, g = s.Function("r")(x), s.Function("g")(x)
    field = g/r
    frame = 1-alpha*field**2/6
    stress_numerator = s.diff(field, x)**2-s.diff(field**2, x, 2)/6
    # Affine Jordan Ric_kk=-2*r''/r and F*Tfull=Akk, before solving for r''.
    constraint = -2*frame*s.diff(r, x, 2)/r-alpha*stress_numerator
    expected = alpha*((g*s.diff(g, x, 2)-2*s.diff(g, x)**2)/(6*r)
                      +g*s.diff(g, x)*s.diff(r, x)/(3*r**2))
    assert s.simplify(constraint.subs(s.diff(r, x, 2), expected)) == 0
    assert expected.subs(alpha, 0) == 0


def test_profile_derivatives_and_source_gronwall_coefficients(budget):
    r, v, g, p, z = s.symbols("r v g p z", real=True)
    A = (g*z-2*p**2)/(6*r)+g*p*v/(3*r**2)
    expected = (z/(6*r)+p*v/(3*r**2),
                -2*p/(3*r)+g*v/(3*r**2), g/(6*r))
    for variable, derivative in zip((g, p, z), expected):
        assert s.simplify(s.diff(A, variable)-derivative) == 0
    dr, dv, dg, dp, dz, epsilon = s.symbols("dr dv dg dp dz epsilon", real=True)
    varied = A.subs({r:r+epsilon*dr, v:v+epsilon*dv, g:g+epsilon*dg,
                    p:p+epsilon*dp, z:z+epsilon*dz}, simultaneous=True)
    linear = (s.diff(A, r)*dr+s.diff(A, v)*dv
              +expected[0]*dg+expected[1]*dp+expected[2]*dz)
    assert s.simplify(s.diff(varied, epsilon).subs(epsilon, 0)-linear) == 0
    G0, G1, G2, floor, vmax = s.symbols("G0 G1 G2 floor vmax", positive=True)
    term_bounds = (G2/(6*floor)+G1*vmax/(3*floor**2),
                   2*G1/(3*floor)+G0*vmax/(3*floor**2), G0/(6*floor))
    Bg = (G2+4*G1+G0)/(6*floor)+(G1+G0)*vmax/(3*floor**2)
    assert s.simplify(sum(term_bounds)-Bg) == 0
    t, L, forcing = s.symbols("t L forcing", positive=True)
    exact = forcing*(s.exp(L*t)-1)/L
    upper = forcing*t*s.exp(L*t)
    assert s.simplify(s.diff(exact, t)-L*exact-forcing) == 0
    assert exact.subs(t, 0) == upper.subs(t, 0) == 0
    assert s.simplify(s.diff(upper-exact, t)-forcing*L*t*s.exp(L*t)) == 0


def test_data_are_propagated_through_each_slab(budget):
    b, d, delta = s.symbols("b d delta", positive=True)
    q = s.Rational(1, 2)
    A, B = b/(1-q), d/(1-q)
    unknown = s.symbols("unknown")
    previous = s.Integer(0)
    for count in range(1, 5):
        # An affine manufactured fixed point attains the comparison recurrence.
        value = s.solve(s.Eq(unknown, q*unknown+b*previous+d*delta), unknown)[0]
        expected = B*delta*sum(A**i for i in range(count))
        assert s.simplify(value-expected) == 0
        assert value.subs(delta, 0) == 0
        previous = value
    # Omitting inherited data fails after the first slab.
    assert s.simplify(previous-B*delta) != 0


def test_changing_characteristic_speed_keeps_other_solution_gradient(budget):
    x, eta, a, b = s.symbols("x eta a b", real=True)
    profile = lambda z: z**3+2*z
    w1, w2 = profile(eta+a*x), profile(eta+b*x)
    difference = w2-w1
    left = s.diff(difference, x)-a*s.diff(difference, eta)
    right = (b-a)*s.diff(w2, eta)
    assert s.simplify(left-right) == 0
    wrong = (b-a)*s.diff(w1, eta)
    assert s.simplify(left-wrong) != 0
