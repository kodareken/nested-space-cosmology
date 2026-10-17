"""Exact full-stress normalization; no original pulse integration or evolution."""
import pytest
import sympy as s

from control_budget import ControlBudget


CONTROLS = ControlBudget(maximum_evaluations=0, cpu_seconds=10.)


@pytest.fixture
def budget():
    with CONTROLS.active() as value:
        yield value


@pytest.fixture(scope="module", autouse=True)
def report():
    yield
    print("null normalization controls:", CONTROLS.budget.evaluations,
          "field samples;", CONTROLS.active_cpu_seconds, "active CPU seconds")


def test_full_stress_elimination_boundary_identity_and_weight(budget):
    x = s.symbols("x", real=True)
    kappa = s.symbols("kappa", positive=True)
    phi = s.Function("phi")(x)
    F = 1 - kappa*phi**2/6
    A = s.diff(phi, x)**2 - s.diff(phi**2, x, 2)/6
    T = s.symbols("T")
    # Solve T=A+phi^2*G/6 together with G=kappa*T, retaining backreaction.
    full = s.solve(T - A - kappa*phi**2*T/6, T)[0]
    assert s.cancel(full - A/F) == 0
    W = (2 + F)/(3*F**2)
    primitive = phi*s.diff(phi, x)/(3*F)
    assert s.cancel(full - W*s.diff(phi, x)**2
                    + s.diff(primitive, x)) == 0
    assert s.cancel(W - (1-kappa*phi**2/18)/F**2) == 0
    z = s.symbols("z", positive=True)
    weight = (2+z)/(3*z**2)
    assert s.cancel(weight-1-(1-z)*(3*z+2)/(3*z**2)) == 0
    assert s.cancel(s.diff(weight, z)+(z+4)/(3*z**3)) == 0


def test_compact_polynomial_weak_limit_and_nonquiet_endpoint(budget):
    x, epsilon = s.symbols("x epsilon", real=True)
    kappa = s.symbols("kappa", positive=True)
    p = (1-x*x)**4
    dp = s.diff(p, x)
    for endpoint in (-1, 1):
        assert p.subs(x, endpoint) == dp.subs(x, endpoint) == 0
    phi = epsilon*p
    A = s.diff(phi, x)**2 - s.diff(phi**2, x, 2)/6
    full = A/(1-kappa*phi**2/6)
    weak = s.diff(full, epsilon, 2).subs(epsilon, 0)/2
    expected = s.Rational(2, 3)*dp**2-p*s.diff(p, x, 2)/3
    assert s.expand(weak-expected) == 0
    J0 = s.integrate(s.expand(dp**2), (x, -1, 1))
    assert J0 > 0
    assert s.integrate(s.expand(weak), (x, -1, 1)) == J0
    assert full.subs(epsilon, 0) == 0
    # A nonquiet interval must retain the endpoint term, already at kappa=0.
    nonquiet = 1+x
    numerator = s.diff(nonquiet, x)**2-s.diff(nonquiet**2, x, 2)/6
    primitive = nonquiet*s.diff(nonquiet, x)/3
    integral = s.integrate(numerator, (x, 0, 1))
    positive = s.integrate(s.diff(nonquiet, x)**2, (x, 0, 1))
    boundary = primitive.subs(x, 1)-primitive.subs(x, 0)
    assert integral == positive-boundary == s.Rational(2, 3)
    assert integral != positive


def test_affine_anchor_and_dimensional_scaling(budget):
    p, v, a = s.symbols("p v a", real=True)
    kappa, c, omega, ell = s.symbols("kappa c omega ell", positive=True)
    epsilon = s.symbols("epsilon", real=True)

    def stress(phi, first, second):
        return (2*first**2-phi*second)/(3*(1-kappa*phi**2/6))

    original = stress(p, v, a)
    transformed = stress(p, v/c, a/c**2)
    assert s.cancel(c*transformed-original/c) == 0
    assert s.cancel(c*transformed/(omega/c)-original/omega) == 0
    alpha = kappa*epsilon**2/ell**2
    F = 1-alpha*p**2/6
    physical_positive = (epsilon*v/ell**2)**2*(2+F)/(3*F**2)
    dimensionless_density = v**2*(2+F)/(3*F**2)
    assert s.cancel(kappa*ell*ell*physical_positive
                    -alpha*dimensionless_density) == 0
    # Contracting the tensor with dx differs from the tangent d/dlambda.
    assert s.cancel(ell**2*stress(epsilon*p/ell, epsilon*v/ell**2,
                                epsilon*a/ell**3)
                    -stress(epsilon*p/ell, epsilon*v/ell,
                            epsilon*a/ell)) == 0


def test_positive_weight_does_not_order_changing_profiles(budget):
    x, alpha = s.symbols("x alpha", real=True)
    p = (1-x*x)**4
    # Logical control only: this rescaling is not a sourced radius solution.
    q = p/(1+alpha)
    F = 1-alpha*q*q/6
    density = s.diff(q, x)**2*(2+F)/(3*F**2)
    derivative = s.diff(density, alpha).subs(alpha, 0)
    expected = s.diff(p, x)**2*(-2+s.Rational(5, 18)*p*p)
    assert s.expand(derivative-expected) == 0
    # On [-1,1], 0<=p<=1: the bracket is at most -31/18.
    assert s.integrate(s.expand(derivative), (x, -1, 1)) < 0
