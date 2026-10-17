"""Exact past-collar and Cauchy-data controls; no physical field calculation."""
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
    print("past preparation controls:", CONTROLS.budget.evaluations,
          "physical samples;", CONTROLS.active_cpu_seconds, "active CPU seconds")


def test_oriented_primitive_keeps_source_transverse_jet_and_bottom(budget):
    u, v, t, w = s.symbols("u v t w", real=True)
    # Manufactured H=Y+P+Q, Y=exp(2u+3v), including a nonzero bottom.
    Y = s.exp(2*u+3*v)
    S, B0 = s.exp(3*v), s.exp(2*u)
    Q = s.diff(Y, v)
    P = s.diff(B0, u)*s.exp(v)+s.integrate(
        s.exp(v-t)*(Y+Q).subs(v, t), (t, 0, v))
    assert s.simplify(P-s.diff(Y, u)) == 0
    H = Y+P+Q
    newY = S+s.integrate(P.subs(u, w), (w, 0, u))
    newQ = s.diff(S, v)+s.integrate(H.subs(u, w), (w, 0, u))
    assert s.simplify(newY-Y) == 0
    assert s.simplify(newQ-Q) == 0
    assert s.simplify(newY.subs(v, 0)-B0) == 0
    assert s.simplify(s.diff(P, v)-H) == 0
    assert s.simplify(P.subs(u, 0)-s.diff(B0, u).subs(u, 0)) != 0
    # Under U=-u the old mixed equation has the opposite sign.
    U = s.symbols("U", real=True)
    forward = Y.subs(u, -U)
    assert s.simplify(s.diff(forward, U, v)+H.subs(u, -U)) == 0
    # A conserved r*constraint remains conserved in the reversed coordinate.
    r = 2+u+v
    E = 7/r
    assert s.diff(r*E, u) == 0
    assert s.diff((r*E).subs(u, -U), U) == 0


def test_transverse_ode_lipschitz_constants_are_attained(budget):
    u, v, t, epsilon = s.symbols("u v t epsilon", positive=True)
    # H=Y+P+Q: equal perturbations epsilon in Y,Q and fixed bottom P.
    deltaP = s.integrate(s.exp(v-t)*2*epsilon, (t, 0, v))
    assert s.simplify(deltaP-2*epsilon*(s.exp(v)-1)) == 0
    deltaYnew = u*deltaP
    deltaQnew = u*(2*epsilon+deltaP)
    assert s.simplify(deltaYnew/epsilon-2*u*(s.exp(v)-1)) == 0
    assert s.simplify(deltaQnew/epsilon-2*u*s.exp(v)) == 0
    # Both contractions need the small u-width; a long v-leg alone is not small.
    assert deltaYnew.subs(u, 0) == deltaQnew.subs(u, 0) == 0
    E, B, rho, delta, eta = s.symbols("E B rho delta eta", positive=True)
    displacement = E*B*delta+2*(E-1)*eta
    # On the E>=2 branch, the declared choices leave half the P tube unused.
    saturated = displacement.subs({delta: rho/(4*E*B), eta: rho/(8*(E-1))})
    assert s.cancel(saturated-rho/2) == 0
    # On 1<=E<=2 the eta=rho/8 branch is no larger than rho/2.
    other = displacement.subs({delta: rho/(4*E*B), eta: rho/8})
    assert s.cancel(rho/2-other-rho*(2-E)/4) == 0


def test_jordan_spacelike_frame_and_scalar_jet_reconstruction(budget):
    C, d, ell, F = s.symbols("C d ell F", positive=True)
    metric = s.Matrix([[0, -ell**2*C/2], [-ell**2*C/2, 0]])
    tangent = s.Matrix([-d, 1])
    normal = s.Matrix([d, 1])/(ell*s.sqrt(C*d))
    assert s.simplify((normal.T*metric*normal)[0]) == -1
    assert s.simplify((normal.T*metric*tangent)[0]) == 0
    assert s.simplify((tangent.T*metric*tangent)[0]-ell**2*C*d) == 0
    U, x = s.symbols("U x", real=True)
    c, r = s.Function("c")(U, x), s.Function("r")(U, x)
    epsilon = s.symbols("epsilon", real=True)

    def n(value):
        return normal[0]*s.diff(value, U)+normal[1]*s.diff(value, x)

    phi = epsilon*c/(ell*r)
    assert s.simplify(n(phi)-epsilon/ell*(n(c)/r-c*n(r)/r**2)) == 0
    tangential = -d*s.diff(phi, U)+s.diff(phi, x)
    assert s.simplify((ell*s.sqrt(C*d)*n(phi)-tangential)/(2*d)
                      -s.diff(phi, U)) == 0
    assert s.simplify((ell*s.sqrt(C*d)*n(phi)+tangential)/2
                      -s.diff(phi, x)) == 0
    einstein_normal = normal/s.sqrt(F)
    assert s.simplify((einstein_normal.T*(F*metric)*einstein_normal)[0]) == -1


def test_extrinsic_curvature_from_metric_connection(budget):
    U, x = s.symbols("U x", real=True)
    ell, d = s.symbols("ell d", positive=True)
    C, r = s.Function("C")(U, x), s.Function("r")(U, x)
    coordinates = (U, x)
    metric = s.Matrix([[0, -ell**2*C/2], [-ell**2*C/2, 0]])
    inverse = metric.inv()
    normal = s.Matrix([d, 1])/(ell*s.sqrt(C*d))
    tangent = s.Matrix([-d, 1])/(ell*s.sqrt(C*d))
    connection = [[[s.simplify(sum(inverse[a, e]*(
        s.diff(metric[e, b], coordinates[c])
        +s.diff(metric[e, c], coordinates[b])
        -s.diff(metric[b, c], coordinates[e])) for e in range(2))/2)
        for c in range(2)] for b in range(2)] for a in range(2)]
    derivative = s.Matrix([sum(tangent[b]*(s.diff(normal[a], coordinates[b])
        +sum(connection[a][b][c]*normal[c] for c in range(2)))
        for b in range(2)) for a in range(2)])
    radial_K = -(tangent.T*metric*derivative)[0]
    nC = sum(normal[a]*s.diff(C, coordinates[a]) for a in range(2))
    nr = sum(normal[a]*s.diff(r, coordinates[a]) for a in range(2))
    assert s.simplify(radial_K+nC/(2*C)) == 0
    # Compute Gamma^a_theta_theta directly from the angular metric ell^2*r^2.
    angular_connection = -inverse*s.Matrix(
        [s.diff(ell**2*r**2, q) for q in coordinates])/2
    angular_K = ((metric*normal).T*angular_connection)[0]
    assert s.simplify(angular_K/(ell**2*r**2)+nr/r) == 0
