"""Exact and manufactured controls for the finite boundary jet and lift."""
import time
from types import SimpleNamespace

import numpy as np
import pytest
import sympy as sp
from scipy.interpolate import CubicHermiteSpline

from recursive_horizons import nsc_atlas_scalar_lift as model
from recursive_horizons.nsc_atlas_bulk_response import BulkBudget
from recursive_horizons.nsc_atlas_cap_wall_response import reference_coefficients
from recursive_horizons.nsc_atlas_prepared_response import PreparedFieldBudget
from recursive_horizons.nsc_atlas_scalar_boundary import ScalarBoundary


@pytest.fixture(autouse=True)
def analytic_budget():
    started = time.process_time()
    yield
    assert time.process_time()-started < 30


def test_mass_null_transport_and_potential_derivative_from_independent_constraints():
    r, nu, f, fp = sp.symbols("r nu f fp", positive=True)
    A = 1-r*r
    p, pT, pr = f/r, fp/r, fp/(r*A)-f/r**2
    nuT = r*r*A*pT*pr/2
    nur = r*r*(pT*pT/A+A*pr*pr)/4+r*r*p*p/2
    nuu = sp.factor((nuT-A*nur)/2)
    assert sp.factor(nuu+A*(1+r*r)*f*f/(8*r*r)) == 0
    K = A*nu/(2*r**3)-A*f*fp/(12*r**3)+A*(1-3*r*r)*f*f/(24*r**4)
    derivative = sp.diff(K, nu)*nuu-sp.diff(K, r)*A/2
    expected = (A*(3-r*r)*nu/(4*r**4)-A*(3-r*r)*f*fp/(24*r**4)
                +A*(1-8*r*r+3*r**4)*f*f/(48*r**5))
    assert sp.factor(derivative-expected) == 0
    values = model.cap_potential_u(.5, .2, .1/.5, .3/.5, .3/(.5*.75)-.1/.5**2, .1, .3)
    assert abs(float(values["potential_u_identity_defect"])) <= 4*np.finfo(float).eps*abs(float(values["K1_u"]))
    assert float(values["mass_null_identity_defect"]) == pytest.approx(0., abs=2e-17)


def test_prearrival_profile_moment_is_the_regular_kinetic_and_potential_mass():
    T, r = sp.symbols("T r", real=True)
    q = T+sp.atanh(r)
    f = sp.Function("f")
    C = sp.Function("C")
    A = 1-r*r
    p = f(q)/r
    mass = C(q)-A*f(q)**2/(4*r)
    pT, pr = sp.diff(p, T), sp.diff(p, r)
    # Replace C'(q) structurally, without assuming any field equation beyond
    # the regular no-return profile p=f(q)/r.
    derivative = sp.diff(mass, r)
    Cprime = next(term for term in derivative.atoms(sp.Subs) if term.expr.has(C))
    density = r*r*(pT*pT/A+A*pr*pr)/4+r*r*p*p/2
    assert sp.simplify(derivative.subs(Cprime, sp.diff(f(q), T)**2/2)-density) == 0
    time_derivative = sp.diff(mass, T).subs(Cprime, sp.diff(f(q), T)**2/2)
    assert sp.simplify(time_derivative-r*r*A*pT*pr/2) == 0
    # The quartic onset has f(0)=0 and C(0)=0, fixing the integration constant.
    explicit = 8*q**7/7-A*q**8/(4*r)
    assert sp.simplify(explicit.subs(T, -sp.atanh(r))) == 0


class AnalyticCap:
    def __init__(self):
        self.budget = PreparedFieldBudget()
        self.domain = dict(qa=0., qb=.33, Ta=-np.arctanh(.45),
                           Tb=.33-np.arctanh(.55), qmin=-2*np.arctanh(.45))
        self.edges = np.array([0., .11, .22, .33])

    def radial_limit(self, T):
        return np.full_like(np.asarray(T), .55, dtype=float)

    def profile(self, q):
        q = np.asarray(q)
        return {"f": q**4, "f_prime": 4*q**3}


def test_nonconstant_manufactured_outgoing_ray_derivative_by_differentiating_fixed_integral(monkeypatch):
    cap = AnalyticCap()
    budget = BulkBudget(started=cap.budget.started, external=cap.budget)

    def analytic_mass(field, T, r, *, radial_order, budget):
        T, r = np.broadcast_arrays(T, r)
        q, A = T+np.arctanh(r), 1-r*r
        f, fp = q**4, 4*q**3
        budget.check(metric=T.size)
        return {"nu_E": 8*q**7/7-A*f*f/(4*r), "p": f/r,
                "p_T": fp/r, "p_r": fp/(r*A)-f/r**2}
    monkeypatch.setattr(model, "radial_metric", analytic_mass)
    s, z = sp.symbols("s z", real=True)
    u, upper = -sp.Rational(9, 10)+s*sp.Rational(3, 100), sp.Rational(1, 5)+s*sp.Rational(1, 25)
    q = upper*z
    r = sp.tanh((q-u)/2)
    A, f, fp, C = 1-r*r, q**4, 4*q**3, 8*q**7/7
    K = A*C/(2*r**3)-A*f*fp/(12*r**3)-A*f*f/(12*r**4)
    fixed_interval_integrand = -K*f*upper
    derivative = sp.lambdify((s, z), sp.diff(fixed_interval_integrand, s), "numpy")
    value = sp.lambdify((s, z), fixed_interval_integrand, "numpy")
    node, weight = np.polynomial.legendre.leggauss(64)
    points, weights = (node+1)/2, weight/2
    location = .4
    result = model.outgoing_cap_jet(cap, float(u.subs(s, location)), float(upper.subs(s, location)),
                                    budget=budget, radial_order=4, ray_order=16)
    computed = -.04*result["K1_wall"]*result["wall_profile"]-.03*result["integral_Ku_f"]
    assert result["J_minus"] == pytest.approx(float(weights @ value(location, points)), rel=5e-13, abs=1e-20)
    assert computed == pytest.approx(float(weights @ derivative(location, points)), rel=8e-13, abs=1e-20)
    assert abs(computed) > 1e-12


def test_wall_boundary_slope_matches_total_derivative_of_unexpanded_expression():
    x, V, w, p, rho, V1, J = sp.symbols("x V w p rho V1 J", real=True)
    m = sp.Rational(2, 5)
    bp = m/(2*x*x)-sp.Rational(5, 4)*x
    bm = bp+2*x
    v = sp.sqrt(bp*bp-1+2*m/x)
    d, A0 = v-bp, 1-x*x
    Vs = -V/(4*m*d)
    c1r = 2*(bm-v)*J/(A0*d)-(rho-x*p*p/12-v*V1/Vs)*2*m*w/(x**3*d*d)
    rU = -4*m*m*V*sp.exp(-x/(2*m))/x
    point = {x: .51, V: -.7, w: .102, p: .2, rho: .04, V1: .03, J: .07}
    reference = reference_coefficients(.51)
    speed, clock = float(reference["v"]), float(Vs.subs(point))
    rates = {x: speed, V: clock, w: .3, p: (.3-speed*.2)/.51,
             rho: -.02, V1: -.05, J: -.13}
    expected = sum(float(sp.diff(rU*c1r, variable).subs(point))*rate for variable, rate in rates.items())
    row = {**{key: float(value) for key, value in reference.items()},
           "V": -.7, "w": .102, "p": .2, "rho": .04, "V1": .03, "Vdot": clock,
           "j": .3, "A": rates[p], "rho_s": -.02, "V1_s": -.05}
    result = model.wall_boundary_derivative(row, .07, -.13)
    assert result["b_s"] == pytest.approx(expected, rel=3e-14, abs=3e-15)
    assert result["b_V"] == pytest.approx(expected/clock, rel=3e-14, abs=3e-15)


class PolynomialJets(model.BoundaryJet):
    """Analytic cubic b(U) with a nonlinear, monotone reference wall map."""
    def __init__(self):
        external = PreparedFieldBudget()
        self.budget = BulkBudget(started=external.started, external=external)
        self.cap_field = SimpleNamespace(budget=external)
        self.domain = dict(U0=-3., Ua=-2., Ub=0., Va=-2., Vb=-1.)
        self.geometry = {"wall": CubicHermiteSpline([-2., -1.], [-2., 0.], [1., 3.])}
        self.settings = dict(radial_order=4, ray_order=4)
        self.cache, self.ray_nodes = {}, 0

    def sample(self, t):
        t = float(t)
        if t not in self.cache:
            V = -2+t
            U, Uv = float(self.geometry["wall"](V)), float(self.geometry["wall"](V, 1))
            z = U+2
            b, bU = z*z+z**3, 2*z+3*z*z
            self.cache[t] = dict(t=t, U=U, V=V, b=b, b_U=bU, b_V=bU*Uv, U_V=Uv, J_minus_s=0.)
        return dict(self.cache[t])


def test_quiet_adapter_boundary_jet_and_lift_are_exactly_zero_and_share_cache():
    cap = AnalyticCap()
    cap.domain.update(U0=-3., Ua=-2., Ub=-1., Va=-2., Vb=-1., xa=.45, xb=.55)
    cap.geometry = {"wall": CubicHermiteSpline([-2., -1.], [-2., -1.], [1., 1.]),
                    "arrival": CubicHermiteSpline([-2., -1.], [0., .33], [.33, .33])}
    cap.prepared_field = object()
    cap.quiet_certificate = {"quiet_cap": True}
    cap.profile = lambda q: {"f": np.zeros_like(np.asarray(q)), "f_prime": np.zeros_like(np.asarray(q))}
    cap.jet = lambda T, r: {key: np.zeros_like(np.broadcast_arrays(T, r)[0])
                            for key in ("p", "p_T", "p_r", "p_center")}
    budget = BulkBudget(started=cap.budget.started, external=cap.budget)
    def wall_sample(t):
        x, V = .45+.1*t, -2+t
        row = {key: float(value) for key, value in reference_coefficients(x).items()}
        return {**row, "V": V, "Vdot": -V/(1.6*(row["v"]-row["beta_plus"])),
                **dict.fromkeys(("p", "w", "j", "A", "rho", "rho_s", "V1", "V1_s"), 0.)}
    wall = SimpleNamespace(budget=budget, sample=wall_sample, geometry=cap.domain,
                           wedge=SimpleNamespace(prepared_field=cap.prepared_field))
    boundary = ScalarBoundary(cap, wall, radial_order=4, ray_order=4)
    boundary.sample(.5)
    jets = model.BoundaryJet(boundary)
    lift = model.BoundaryLift(jets, knots=5, validation_times=(.125, .375, .625, .875))
    assert len(jets.cache) == 9
    for row in jets.cache.values():
        for key in ("J_minus", "J_minus_s", "b", "b_s", "b_V", "b_U", "c1_r", "c1_r_s"):
            assert row[key] == 0
    assert jets.cache[.5]["boundary_cache_checked"]
    assert jets.cache[.5]["boundary_cache_value_defect"] == 0
    previous = budget.metric_values
    jets.sample(.5)
    assert budget.metric_values == previous
    assert lift.summary()["terminal_lift"] == 0
    assert np.count_nonzero(lift.B(np.linspace(-3., -1., 31))) == 0
    wall.budget = BulkBudget(started=cap.budget.started, external=cap.budget)
    with pytest.raises(ValueError, match="cumulative wall"):
        model.BoundaryJet(boundary)


def test_cubic_boundary_reproduction_exact_lift_and_nonlinear_wall_transport():
    jets = PolynomialJets()
    lift = model.BoundaryLift(jets, knots=5, validation_times=(.1, .4, .6, .9))
    U = np.linspace(-3., 0., 97)
    z = np.maximum(U+2, 0.)
    assert lift.b_at_U(U) == pytest.approx(z*z+z**3, abs=3e-15)
    assert lift.b_U(U) == pytest.approx(2*z+3*z*z, abs=8e-15)
    assert lift.B(U) == pytest.approx(z**3/3+z**4/4, abs=3e-15)
    V = np.linspace(-2., -1., 31)
    mapped = jets.geometry["wall"](V)
    assert lift.B_at_wall(V) == pytest.approx(lift.B(mapped), abs=2e-15)
    assert lift.b(V) == pytest.approx(lift.b_at_U(mapped), abs=2e-15)
    report = lift.summary()
    assert max(report["maximum_absolute_defects"].values()) < 1e-14
    assert report["terminal_lift"] == pytest.approx(20/3, abs=3e-15)
    assert lift.B(-2.) == lift.B(-3.) == lift.b_at_U(-2.) == 0


def test_nested_fits_share_queries_and_withheld_points_are_not_fitted():
    jets = PolynomialJets()
    coarse = model.BoundaryLift(jets, knots=9)
    fine = model.BoundaryLift(jets, knots=17, validation_times=np.linspace(0., 1., 33)[1::2])
    assert len(jets.cache) == 33
    assert fine.knot_count == 17
    assert len(fine.validation_rows) == 16
    assert coarse.B(0.) == pytest.approx(fine.B(0.), abs=3e-15)
    with pytest.raises(ValueError, match="fitted knot"):
        fine.validate([.5])
    count = len(jets.cache)
    fine.validate([1/32])
    assert len(jets.cache) == count


@pytest.mark.parametrize("bad", [True, 9., 3, 65])
def test_unadmitted_knot_counts_fail(bad):
    with pytest.raises(ValueError, match="knot count"):
        model.BoundaryLift(PolynomialJets(), knots=bad)


def test_invalid_domains_shapes_and_nonquiet_anchors_fail():
    lift = model.BoundaryLift(PolynomialJets(), knots=5)
    for invalid in (True, np.nan, [-3.1], [.01], np.zeros((2, 2))):
        with pytest.raises(ValueError):
            lift.B(invalid)
    for invalid in (-2.01, -.99, np.inf):
        with pytest.raises(ValueError):
            lift.b(invalid)
    old_wall = lift.geometry["wall"]
    lift.geometry["wall"] = lambda V: np.asarray(V)*0+1.
    with pytest.raises(ValueError, match="wall map"):
        lift.B_at_wall(-1.5)
    lift.geometry["wall"] = old_wall
    bad = PolynomialJets()
    bad.sample(0.)
    bad.cache[0.]["b_U"] = 1.
    with pytest.raises(ValueError, match="quiet"):
        model.BoundaryLift(bad, knots=5)
    bad = PolynomialJets()
    bad.domain["Ub"] = .1
    with pytest.raises(ValueError, match="full finite"):
        model.BoundaryLift(bad, knots=5)


def test_expired_shared_clock_rejects_even_polynomial_queries():
    jets = PolynomialJets()
    lift = model.BoundaryLift(jets, knots=5)
    jets.budget.started -= 121
    jets.cap_field.budget.started = jets.budget.started
    with pytest.raises(RuntimeError, match="budget"):
        lift.B(-1.)
