"""Exact cap transport/edge-map algebra; no retained state or numerical solve."""
import sympy as s
import pytest

from control_budget import ControlBudget

CONTROLS = ControlBudget(maximum_evaluations=1, cpu_seconds=20.)


@pytest.fixture(autouse=True)
def exact_budget():
    with CONTROLS.active():
        yield


@pytest.fixture(scope="module", autouse=True)
def report():
    yield
    print("cap transport algebra:", CONTROLS.budget.evaluations,
          "field samples;", CONTROLS.active_cpu_seconds, "active CPU seconds")


def test_full_parameter_direction_identity_without_imposing_field_equations():
    r, x, alpha, tx = s.symbols("r x alpha T_x", positive=True)
    n, delta, f, pi, phi = (s.Function(k)(r,x) for k in ("n","delta","f","pi","phi"))
    potential = s.Function("P")(f)
    U = s.diff(potential, f)
    A = 1-r*r-2*alpha*n/r
    B = A*s.exp(delta)
    a = alpha*r*(pi*pi+phi*phi)/2
    b = r*r*potential/2+r*r*(1-r*r)*(pi*pi+phi*phi)/4
    J = r*r*A*B*pi*phi/2
    er = s.diff(n,r)+a*n-b
    ed = s.diff(delta,r)-a
    M = s.diff(n,x)-tx*J
    rf = s.diff(f,x)-tx*B*pi
    rp = s.diff(pi,x)-tx*s.diff(r*r*B*phi,r)/(r*r)+tx*s.exp(delta)*U
    rq = s.diff(phi,x)-tx*s.diff(B*pi,r)
    field = r*r*(U*rf+A*(pi*rp+phi*rq))/2
    difference = s.diff(M,r)+a*M-s.diff(er,x)-tx*J*ed-field
    assert s.expand(difference) == 0
    # An arbitrary reconstruction need not satisfy any of these equations.
    # Dropping the lapse or potential-reaction channel is not an identity.
    assert s.expand(difference+tx*J*ed) != 0
    assert s.expand(difference-r*r*A*pi*tx*s.exp(delta)*U/2) != 0


def test_shadow_wall_conversion_retains_radius_and_flux_differences():
    a,b,ne,error,dne,derror,rx,tx,je,jq = s.symbols("a b ne e dne de Rx Tx JE JQ")
    nq = ne+error
    rhoe, rhoq = b-a*ne, b-a*nq
    # dne denotes the fixed-radius derivative; the edge has speed Rx.
    direct = dne+rx*rhoe+derror-rx*rhoq-tx*jq
    channels = derror+rx*a*error+tx*(je-jq)+(dne-tx*je)
    assert s.expand(direct-channels) == 0
    assert s.expand(direct-(channels-rx*a*error)) == rx*a*error


def test_panel_error_sum_and_derivative_include_entering_mass_and_attenuation():
    x = s.symbols("x")
    exact_a = [s.Function("a"+str(j))(x) for j in range(3)]
    exact_b = [s.Function("b"+str(j))(x) for j in range(3)]
    numeric_a = [s.Function("c"+str(j))(x) for j in range(3)]
    numeric_b = [s.Function("d"+str(j))(x) for j in range(3)]
    exact, numeric = [s.S.Zero], [s.S.Zero]
    local, after = [], []
    for j in range(3):
        exact.append(exact_a[j]*exact[-1]+exact_b[j])
        numeric.append(numeric_a[j]*numeric[-1]+numeric_b[j])
        local.append(numeric[-1]-(exact_a[j]*numeric[j]+exact_b[j]))
        after.append(s.prod(exact_a[j+1:]))
    error = sum(w*d for w,d in zip(after,local))
    assert s.expand(numeric[-1]-exact[-1]-error) == 0
    full = sum(s.diff(w,x)*d+w*s.diff(d,x) for w,d in zip(after,local))
    assert s.expand(s.diff(numeric[-1]-exact[-1],x)-full) == 0
    omitted_attenuation = sum(w*s.diff(d,x) for w,d in zip(after,local))
    assert s.expand(full-omitted_attenuation) != 0
    # The local error is evaluated at the NUMERICAL entering mass. Its
    # derivative includes (c_j-a_j) times the derivative of that mass.
    omitted_entry = full-sum(after[j]*(numeric_a[j]-exact_a[j])*s.diff(numeric[j],x)
                            for j in range(3))
    assert s.expand(full-omitted_entry) != 0


def test_moving_interface_balance_against_an_independent_exact_radial_solution():
    r,x,u = s.symbols("r x u", real=True)
    face = s.Rational(1,5)+x/25
    wall = s.Rational(3,5)+3*x/100
    tx = 2+x
    A = 1-r*r  # alpha=0; normalized mass is still a nontrivial functional.
    f = (1+x)*r*r
    pairs = [((1+x)*r,(2-x)*r), ((3+x)*r,(1+2*x)*r)]
    densities,fluxes,fields = [],[],[]
    for pi,phi in pairs:
        density = r*r*f*f/2+r*r*A*(pi*pi+phi*phi)/4
        flux = r*r*A*A*pi*phi/2
        rf = s.diff(f,x)-tx*A*pi
        rp = s.diff(pi,x)-tx*s.diff(r*r*A*phi,r)/(r*r)+tx*2*f
        rq = s.diff(phi,x)-tx*s.diff(A*pi,r)
        field = r*r*(2*f*rf+A*(pi*rp+phi*rq))/2
        densities.append(density); fluxes.append(flux); fields.append(s.expand(field))
    first = s.integrate(densities[0].subs(r,u),(u,0,face))
    mass = first+s.integrate(densities[1].subs(r,u),(u,face,wall))
    direct = s.diff(mass,x)-s.diff(wall,x)*densities[1].subs(r,wall)-tx*fluxes[1].subs(r,wall)
    bulk = s.integrate(fields[0],(r,0,face))+s.integrate(fields[1],(r,face,wall))
    jump = (s.diff(face,x)*(densities[1]-densities[0])+tx*(fluxes[1]-fluxes[0])).subs(r,face)
    assert s.expand(direct-(bulk-jump)) == 0
    assert s.expand(direct-bulk).subs(x,0) != 0
    assert s.expand(direct-(bulk+jump)).subs(x,0) != 0
