"""Exact chart/commutation controls; no physical evolution or existence assertion."""
import pytest
import sympy as s

from control_budget import ControlBudget


CONTROLS = ControlBudget(maximum_evaluations=24000, cpu_seconds=30.)


@pytest.fixture
def budget():
    with CONTROLS.active() as value:
        yield value


@pytest.fixture(scope='module', autouse=True)
def report():
    yield
    print('null continuation controls:', CONTROLS.budget.evaluations,
          'charged values;', CONTROLS.active_cpu_seconds, 'active CPU seconds')


def test_second_jet_hierarchy_and_nonstraight_initial_curve(budget):
    budget.check(8000)
    u, v = s.symbols('u v', real=True)
    y = s.exp(u + v + u*v)
    p, q = s.diff(y, u), s.diff(y, v)
    P, Q = s.diff(p, u), s.diff(q, v)
    R, S = s.diff(q, v, 2), s.diff(p, u, 2)
    Y, I, O = s.symbols('Y I O', nonzero=True)
    variables = (Y, I, O)
    nonlinear = Y + I*O/Y
    values = {Y: y, I: p, O: q}
    F = nonlinear.subs(values)
    assert s.simplify(F - s.diff(y, u, v)) == 0
    gradient = [s.diff(nonlinear, z).subs(values) for z in variables]
    hessian = s.hessian(nonlinear, variables).subs(values)
    Fu, Fv = s.diff(F, u), s.diff(F, v)
    Zu, Zv = s.Matrix((p, P, F)), s.Matrix((q, F, Q))
    Fy, Fp, Fq = gradient
    outgoing = (Zv.T*hessian*Zv)[0] + Fy*Q + Fp*Fv + Fq*R
    incoming = (Zu.T*hessian*Zu)[0] + Fy*P + Fp*S + Fq*Fu
    assert s.simplify(s.diff(R, u) - outgoing) == 0
    assert s.simplify(s.diff(S, v) - incoming) == 0

    # A curved initial line exercises ell'' and both mixed endpoint terms.
    lower = 1 - v - v*v
    q0 = q.subs(u, lower)
    endpoint = (s.diff(q0, v, 2)
                - s.diff(lower, v, 2)*F.subs(u, lower)
                - s.diff(lower, v)**2*Fu.subs(u, lower)
                - 2*s.diff(lower, v)*Fv.subs(u, lower))
    assert s.simplify(endpoint - R.subs(u, lower)) == 0


def test_second_incoming_entry_uses_outgoing_second_jet(budget):
    budget.check(8000)
    u, v, t = s.symbols('u v t', real=True)
    y = s.exp(u + v + u*v)
    p, q = s.diff(y, u), s.diff(y, v)
    F = s.diff(y, u, v)
    P, Q = s.diff(p, u), s.diff(q, v)
    R, S = s.diff(q, v, 2), s.diff(p, u, 2)
    wall_u, wall_v = 1 + t + t*t, 2 + t + t*t/2
    a, b = s.diff(wall_u, t), s.diff(wall_v, t)
    at_wall = lambda expression: expression.subs({u: wall_u, v: wall_v})
    g = at_wall(p)
    entry_P = (s.diff(g, t) - b*at_wall(F))/a
    entry_S = s.diff(entry_P, t)/a - b*at_wall(s.diff(F, u))/a
    assert s.simplify(entry_P - at_wall(P)) == 0
    assert s.simplify(entry_S - at_wall(S)) == 0
    outgoing_second = (s.diff(a, t)*at_wall(F) + s.diff(b, t)*at_wall(Q)
                       + a*a*at_wall(s.diff(F, u))
                       + 2*a*b*at_wall(s.diff(F, v)) + b*b*at_wall(R))
    assert s.simplify(s.diff(at_wall(q), t, 2) - outgoing_second) == 0


def test_polar_null_metric_scalar_and_wall_normalization(budget):
    budget.check(8000)
    A0, A, beta, jp, H = s.symbols('A0 A beta jp H', positive=True)
    jm = s.symbols('jm', negative=True)
    Ep, Em, root, gap = s.symbols('Ep Em root gap', real=True)
    jacobian = s.Matrix(((beta*jp, jp), (-beta*jm, jm)))
    inverse = jacobian.inv()
    metric = s.diag(-A0*A0*beta*beta/A, A0*A0/A)
    native_metric = s.simplify(inverse.T*metric*inverse)
    C = -A0*A0/(A*jp*jm)
    assert native_metric == s.Matrix(((0, -C/2), (-C/2, 0)))
    dr = s.Matrix(((0, A0),))*inverse
    assert s.simplify(dr - s.Matrix(((A0/(2*jp), A0/(2*jm)),))) == s.zeros(1, 2)
    assert s.simplify(-4*dr[0]*dr[1]/C - A) == 0
    dchi = s.Matrix(((beta*(Ep+Em)/2, (Ep-Em)/2),))
    dc = (root*dchi + gap*s.Matrix(((0, A0),)))*inverse
    expected = s.Matrix((((root*Ep+A0*gap)/(2*jp),
                          (-root*Em+A0*gap)/(2*jm)),))
    assert s.simplify(dc - expected) == s.zeros(1, 2)
    radial_speed, radial_beta = (H-A/H)/2, (H+A/H)/2
    rates = s.simplify(jacobian*s.Matrix((radial_beta/(beta*A0), radial_speed/A0)))
    assert s.simplify(rates - s.Matrix((jp*H/A0, -jm*A/(H*A0)))) == s.zeros(2, 1)
    assert s.simplify(C*rates[0]*rates[1] - 1) == 0
