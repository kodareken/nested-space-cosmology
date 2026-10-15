"""Independent manufactured cap transport and fixed-clock reconstruction."""
import math
from types import SimpleNamespace

import numpy as np
import pytest
from scipy.integrate import quad

from recursive_horizons import nsc_atlas_central_scalar as model
from recursive_horizons.nsc_atlas_bulk_response import BulkBudget
from recursive_horizons.nsc_atlas_cap_wall_response import reference_coefficients
from recursive_horizons.nsc_atlas_prepared_response import PreparedFieldBudget


class PolynomialProfile:
    """Exact free reduced wave f(q)=scale*q**degree, not a physical source."""
    def __init__(self, degree=1, scale=1.):
        self.degree, self.scale = degree, scale
        self.budget = PreparedFieldBudget(cpu_seconds=30.)
        self.prepared_field = object()
        self.domain = dict(qmin=-1., qa=.7, qb=1.5, Ta=-1., Tb=1., xa=.45, xb=.55)
        self.edges = np.array([.7, .8, 1., 1.5])
        self.profile_sizes, self.jet_sizes = [], []

    def profile(self, q, *, second=False):
        q = np.asarray(q)
        assert q.size <= 513
        if np.any(q < -1.) or np.any(q > 1.5):
            raise ValueError("manufactured profile exits its supplied interval")
        self.profile_sizes.append(q.size)
        n, s = self.degree, self.scale
        result = dict(f=s*q**n, f_prime=s*n*q**(n-1))
        if second:
            result["f_second"] = np.zeros_like(q) if n == 1 else np.full_like(q, 2*s)
        return result

    def jet(self, T, r):
        T, r = np.broadcast_arrays(T, r)
        assert T.size <= 513
        if np.any(r < 0) or np.any(r > .55) or np.any(T-np.arctanh(r) < -1.) or np.any(T+np.arctanh(r) > 1.5):
            raise ValueError("manufactured cap domain exceeded")
        self.jet_sizes.append(T.size)
        ratio = np.divide(np.arctanh(r), r, out=np.ones_like(r), where=r != 0)
        # derivative of atanh(r)/r, with an independent convergent central series
        derivative = np.zeros_like(r)
        small = r <= 1e-4
        derivative[small] = sum(2*k*r[small]**(2*k-1)/(2*k+1) for k in range(1, 6))
        rr = r[~small]
        derivative[~small] = (rr/(1-rr*rr)-np.arctanh(rr))/rr**2
        if self.degree == 1:
            factor, factor_T = np.full_like(T, 2*self.scale), np.zeros_like(T)
        else:
            factor, factor_T = 4*self.scale*T, np.full_like(T, 4*self.scale)
        return dict(p=factor*ratio, p_T=factor_T*ratio, p_r=factor*derivative, p_center=factor)


def shared(field):
    return BulkBudget(cpu_seconds=30., started=field.budget.started, external=field.budget)


def reference_wall(x=.5):
    state = {k: float(v) for k, v in reference_coefficients(x).items()}
    V = -2.
    Vs = -V/(1.6*(state["v"]-state["beta_plus"]))
    return {**state, "V": V, "Vdot": Vs, "V1": 0., "V1_s": 0.}


def zero_inputs():
    wall = reference_wall()
    cap = dict(q0=.9, x=.5, e=0., e_s=0., T1=0., T1_s=0., nu_E_wall=0., delta_E_wall=0., ray_integral=0.)
    profile = dict(f=0., f_prime=0., f_second=0.)
    integrals = dict(q=.9, x=.5, forcing_integral=0., optical_time_integral=0., lambda_plus_integral=0.)
    return wall, cap, profile, integrals


def reconstruct(wall, cap, profile, integrals, **overrides):
    Q = (wall["beta_minus"]+wall["v"])/(1-wall["x"]**2)
    options = dict(y_wall=0., y_wall_s=0., J_minus=0., arrival_spline_derivative=Q/wall["Vdot"])
    return model.reconstruct_central_scalar(wall, cap, profile, integrals, **(options | overrides))


def test_optical_time_derivative_matches_independent_quadratic_metric_coefficient():
    field = PolynomialProfile(degree=2)
    T, r = .63, .47
    result = model.optical_time_derivative(field, T, r, radial_order=32, budget=shared(field))
    # For f=q^2, p=T*P(r). Extract the exact T^2 coefficients of the
    # original radial Einstein constraints, rather than differencing a clock.
    def P(z):
        return 4*math.atanh(z)/z if z else 4.
    def Pr(z):
        return 4*(z/(1-z*z)-math.atanh(z))/z**2 if z else 0.
    nu2 = quad(lambda z: z*z*(1-z*z)*Pr(z)**2/4+z*z*P(z)**2/2, 0., r, epsabs=1e-13)[0]
    delta2 = -16/12+quad(lambda z: z*Pr(z)**2/2, 0., r, epsabs=1e-13)[0]
    exact = 2*T*(delta2-2*nu2/(r*(1-r*r)))
    assert float(result["D_T"]) == pytest.approx(exact, abs=2e-13)
    assert float(result["central_anchor"]) == pytest.approx(-16*T/6)
    assert float(result["central_profile_defect"]) == 0.
    centre = model.optical_time_derivative(field, T, 0., budget=shared(field))
    assert centre["radial_nodes"] == 0
    assert float(centre["D_T"]) == pytest.approx(-16*T/6)


def test_incoming_integral_matches_closed_static_mass_and_independent_null_rule():
    field = PolynomialProfile()
    q, x = .9, .5
    result = model.incoming_response_integrals(field, q, x, radial_order=32, ray_order=32, budget=shared(field))
    # f=q gives p=2*atanh(r)/r. The exact regular mass is
    # nu=atanh(r)-(1-r^2)*atanh(r)^2/r; no radial_metric call in this oracle.
    def Kc(u):
        z = (q-u)/2
        r = math.tanh(z)
        if r < 1e-5:
            K = 4*r*r/135+4*r**4/189
        else:
            A = 1-r*r
            mass = z-A*z*z/r
            p, pr = 2*z/r, 2*(r/A-z)/r**2
            K = A*(2*mass/r**3-p*p/3+A*pr*pr/6)/4
        return K*2*z
    expected = 2*quad(Kc, q-2*math.atanh(x), q, epsabs=2e-13)[0]
    assert result["forcing_integral"] == pytest.approx(expected, abs=2e-12)
    assert result["forcing_integral"] > 0
    assert result["forcing_negative"] == pytest.approx(0., abs=1e-12)
    assert result["optical_time_integral"] == 0
    # Dropping f(u) would give c0=q rather than 2*atanh(r), including at centre.
    assert result["reflected_profile_defect"] < 5e-16
    assert max(field.profile_sizes+field.jet_sizes) <= 513


def test_forced_null_polynomial_cancels_at_regular_centre_with_moving_wall():
    wall, cap, profile, integrals = zero_inputs()
    x, q, sigma = wall["x"], cap["q0"], .73
    A = 1-x*x
    Q = (wall["beta_minus"]+wall["v"])/A
    a_prime = (wall["beta_minus"]-wall["v"])/(A*Q)
    separation = 2*math.atanh(x)
    # c1=g(v)-g(u)-sigma*(v-u)^5/20, S=sigma*(v-u)^3.
    # Choose g(u)=0 on the wall's negative outgoing label, g(v)=v^3.
    H = q**3-sigma*separation**5/20
    Hq = 3*q*q-sigma*separation**4*(1-a_prime)/4
    J = sigma*separation**4/4
    integrals["forcing_integral"] = -sigma*separation**4/2
    result = reconstruct(wall, cap, profile, integrals, y_wall=H, y_wall_s=Q*Hq, J_minus=J)
    assert result["p1_center"] == pytest.approx(6*q*q, abs=2e-15)
    assert abs(result["outgoing_wall_term"]) > .01
    assert result["incoming_forcing_term"] > 0
    assert result["p1_center"] == pytest.approx(sum(result[k] for k in (
        "exterior_trace_term", "arrival_rate_term", "arrival_shift_term", "outgoing_wall_term", "incoming_forcing_term")))


def test_closed_optical_clock_derivative_and_current_exterior_clock_pullback():
    wall, cap, profile, integrals = zero_inputs()
    wall.update(V1=.23, V1_s=-.17)
    cap.update(e=.04, e_s=-.02, T1=.11, T1_s=.07)
    profile.update(f=.4, f_prime=.8, f_second=.7)
    q, x, a = cap["q0"], wall["x"], .21
    xi, A, v, Vs = math.atanh(x), 1-x*x, wall["v"], wall["Vdot"]
    # Manufactured D(T,r)=a*T has a closed incoming null integral.
    cap["delta_E_wall"] = a*(q-xi)
    integrals["lambda_plus_integral"] = a*(xi*xi/2-q*xi)
    integrals["optical_time_integral"] = a*xi
    result = reconstruct(wall, cap, profile, integrals, y_wall=.3, y_wall_s=.2)
    Q, Qs, Vss = result["q_plus_s"], result["q_plus_ss"], result["V0_ss"]
    lambda_s = a*((xi-q)*v/A-Q*xi)
    Delta_s = .07-.02/A+.04*2*x*v/A**2+lambda_s
    ratio_s = (Qs*Vs-Q*Vss)/Vs**2
    Cs = Delta_s-ratio_s*.23-(Q/Vs)*(-.17)
    assert result["Delta_C_fixed_s_s"] == pytest.approx(Delta_s, abs=2e-15)
    assert result["Delta_C_fixed_V_s"] == pytest.approx(Cs, abs=2e-15)
    assert result["H_q"] == pytest.approx(.2/Q-Cs*.8/Q-result["Delta_C_fixed_V"]*.7)
    # No extra -C*p0' is appended after this fixed-clock reconstruction.
    assert result["p1_center"] == pytest.approx(2*result["H_q"])
    assert result["central_frame_algebra_defect"] == pytest.approx(0., abs=2e-15)
    changed = reconstruct({**wall, "V1": -.1}, cap, profile, integrals, y_wall=.3, y_wall_s=.2)
    assert changed["Delta_C_fixed_V"]-result["Delta_C_fixed_V"] == pytest.approx(Q/Vs*.33)
    mismatch = reconstruct(wall, cap, profile, integrals, arrival_spline_derivative=Q/Vs+.001)
    assert mismatch["arrival_rate_defect"] == pytest.approx(.001*Vs)
    assert mismatch["arrival_continuum_derivative"] == result["arrival_continuum_derivative"]


def test_zero_source_and_cubic_field_scaling():
    zero = PolynomialProfile(scale=0.)
    integrals = model.incoming_response_integrals(zero, .9, .5, radial_order=8, ray_order=8, budget=shared(zero))
    assert integrals["forcing_integral"] == integrals["optical_time_integral"] == integrals["lambda_plus_integral"] == 0
    wall, cap, profile, _ = zero_inputs()
    assert reconstruct(wall, cap, profile, integrals)["p1_center"] == 0
    first, second = PolynomialProfile(degree=2), PolynomialProfile(degree=2, scale=-.6)
    a = model.incoming_response_integrals(first, .9, .5, budget=shared(first))
    b = model.incoming_response_integrals(second, .9, .5, budget=shared(second))
    assert b["forcing_integral"] == pytest.approx(-.6**3*a["forcing_integral"], rel=2e-11)
    for key in ("optical_time_integral", "lambda_plus_integral"):
        assert b[key] == pytest.approx(.6**2*a[key], rel=2e-12)


def test_incoming_rejects_outside_domain_and_optical_queries_are_bounded():
    field = PolynomialProfile()
    for q, x in ((.69, .5), (1.51, .5), (.9, .6), (True, .5)):
        with pytest.raises(ValueError):
            model.incoming_response_integrals(field, q, x, budget=shared(field))
    for T, r in ((True, .2), (.9, -.1), ([.9]*514, [.2]*514), (.9, math.nan)):
        with pytest.raises(ValueError):
            model.optical_time_derivative(field, T, r, budget=shared(field))
    with pytest.raises(ValueError, match="shared"):
        model.optical_time_derivative(field, .9, .2, budget=BulkBudget())
    wall, cap, profile, integrals = zero_inputs()
    with pytest.raises(ValueError, match="coordinates"):
        reconstruct(wall, cap, profile, {**integrals, "q": .91})


def test_live_adapter_reuses_geometry_but_not_supplied_scalar_coefficients():
    field = PolynomialProfile(scale=0.)
    budget = shared(field)
    wall, cap, _, _ = zero_inputs()
    Q = (wall["beta_minus"]+wall["v"])/(1-wall["x"]**2)
    field.geometry = {"arrival": lambda V, derivative=0: Q/wall["Vdot"] if derivative else .9}
    calls = []
    wall_solution = SimpleNamespace(budget=budget, wedge=SimpleNamespace(prepared_field=field.prepared_field),
        sample=lambda t: calls.append(("wall", t)) or dict(wall))
    cap_solution = SimpleNamespace(cap_field=field, wall_solution=wall_solution, budget=budget,
        domain=field.domain, settings=dict(radial_order=8, ray_order=8, time_panels=8),
        sample=lambda t: calls.append(("cap", t)) or dict(cap))
    response = model.CentralScalarResponse(cap_solution)
    quiet = response.sample(.5, y_wall=0., y_wall_s=0., J_minus=0.)
    metric_count = budget.metric_values
    changed = response.sample(.5, y_wall=0., y_wall_s=Q, J_minus=0.)
    assert quiet["p1_center"] == 0 and changed["p1_center"] == 2
    assert calls == [("cap", .5), ("wall", .5)]
    assert budget.metric_values == metric_count
    assert quiet["endpoint_is_one_sided_formal_limit"] is False
    assert response.budget is budget
