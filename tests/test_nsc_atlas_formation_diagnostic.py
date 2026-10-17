"""Exact formation identities and scope controls; no physical field solve."""
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
    print("formation diagnostic controls:", CONTROLS.budget.evaluations,
          "physical samples;", CONTROLS.active_cpu_seconds, "active CPU seconds")


def test_quiet_sphere_expansions_from_metric_area_and_mass(budget):
    ell, d, boost = s.symbols("ell d boost", positive=True)
    r = s.Rational(2, 5)
    gradient = ell*s.Matrix([s.Rational(1, 4), 1])
    metric = s.Matrix([[0, -ell**2/2], [-ell**2/2, 0]])
    plus, minus = s.Matrix([0, 2/ell]), s.Matrix([2/ell, 0])
    R = ell*r
    area_derivative = lambda k: 8*s.pi*R*(k.dot(gradient))
    expansion = lambda k: s.cancel(area_derivative(k)/(4*s.pi*R**2))
    assert (plus.T*metric*minus)[0] == -2
    assert expansion(plus) == 10/ell
    assert expansion(minus) == s.Rational(5, 2)/ell
    assert s.simplify(expansion(boost*plus)*expansion(minus/boost)
                      -25/ell**2) == 0
    norm = (gradient.T*metric.inv()*gradient)[0]
    mass = s.cancel(R*(1-norm)/2)
    assert mass == R
    assert s.simplify(expansion(plus)*expansion(minus)
                      +4*(1-2*mass/R)/R**2) == 0
    normal = s.Matrix([d, 1])/(ell*s.sqrt(d))
    tangent = s.Matrix([-d, 1])/(ell*s.sqrt(d))
    assert s.simplify(expansion(normal+tangent)-10/(ell*s.sqrt(d))) == 0
    assert s.simplify(expansion(normal-tangent)-5*s.sqrt(d)/(2*ell)) == 0


def test_affine_constraint_focusing_and_regular_centre_control(budget):
    x = s.symbols("x", real=True)
    R, C, psi = [s.Function(name)(x) for name in ("R", "C", "psi")]
    # A coordinate null ray has affine k=C^-1 d/dx (constant scale omitted).
    affine = lambda f: s.diff(f, x)/C
    constraint = s.diff(C, x)*s.diff(R, x)/C-R*s.diff(psi, x)**2/2
    theta = 2*affine(R)/R
    residual = affine(theta)+theta**2/2+affine(psi)**2
    assert s.simplify(residual.subs(s.diff(R, x, 2), constraint)) == 0
    R0, ell = s.symbols("R0 ell", positive=True)
    # Flat Cartesian null curve (t,z)=(x,R0-x) reaches a regular axis.
    tangent = s.Matrix([1, -1])
    assert (tangent.T*s.diag(-1, 1)*tangent)[0] == 0
    flat_radius = R0-x
    flat_theta = 2*s.diff(flat_radius, x)/flat_radius
    assert s.diff(flat_radius, x, 2) == 0
    assert s.simplify(s.diff(flat_theta, x)+flat_theta**2/2) == 0
    assert flat_radius.subs(x, R0) == 0
    assert 2/(10/ell) == ell/5
    assert 2/(s.Rational(5, 2)/ell) == 4*ell/5


def test_surface_stress_and_israel_give_same_crossing_impulse(budget):
    q, sigma, kappa, radius = s.symbols("q sigma kappa radius", positive=True)
    # Local orthonormal wall frame: N=(0,1), radial null k=(q,+/-q).
    metric = s.diag(-1, 1)
    N = s.Matrix([0, 1])
    projector = metric-(metric*N)*(metric*N).T
    stress_coefficient = -sigma*projector
    jump_normal_radius = -kappa*sigma*radius/2
    for orientation in (-1, 1):
        k = s.Matrix([q, orientation*q])
        normal_speed = (N.T*metric*k)[0]
        assert (k.T*metric*k)[0] == 0
        Tkk = (k.T*stress_coefficient*k)[0]
        assert Tkk == sigma*normal_speed**2
        # delta(eta(lambda)) contributes 1/abs(N.k); crossing order reverses.
        raychaudhuri_jump = -kappa*Tkk/s.Abs(normal_speed)
        israel_jump = 2*orientation*normal_speed*jump_normal_radius/radius
        assert s.simplify(raychaudhuri_jump-israel_jump) == 0
        assert raychaudhuri_jump == -kappa*sigma*q


def test_conformal_connection_affine_scale_and_global_bound_control(budget):
    t, z = s.symbols("t z", real=True)
    F = s.Function("F")(t, z)
    metric = s.diag(-1, 1)
    coordinates = (t, z)
    inverse = metric.inv()
    # Compute the conformal connection from its metric, not a recalled rule.
    conformal = F*metric
    inverse_conformal = inverse/F
    connection = [[[sum(inverse_conformal[a, e]*(
        s.diff(conformal[e, b], coordinates[c])
        +s.diff(conformal[e, c], coordinates[b])
        -s.diff(conformal[b, c], coordinates[e])) for e in range(2))/2
        for c in range(2)] for b in range(2)] for a in range(2)]
    kJ = s.Matrix([1, 1])
    kE = kJ/F
    for a in range(2):
        acceleration = sum(kE[b]*s.diff(kE[a], coordinates[b])
            +sum(connection[a][b][c]*kE[b]*kE[c] for c in range(2))
            for b in range(2))
        assert s.simplify(acceleration) == 0
    u = s.symbols("u", nonnegative=True)
    parameter = s.integrate(1/(1+u)**2, (u, 0, u))
    assert s.simplify(parameter-u/(1+u)) == 0
    assert s.limit(parameter, u, s.oo) == 1
    assert s.limit(1/(1+u)**2, u, s.oo) == 0
