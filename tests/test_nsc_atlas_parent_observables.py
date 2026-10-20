"""Independent geometric and observer controls for the parent atlas."""
import math

import numpy as np
import pytest
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq

from recursive_horizons.nsc_atlas_parent_observables import (
    DustParent, de_sitter_bulk_slice, population_balance,
    projected_brane_source, quadratic_brane_stress, shadow_lensing,
    critical_tension_wall, tension_potential,
)


@pytest.mark.parametrize("tau", [0., .01, 1., 4.])
def test_darmois_matching_by_direct_contraction(tau):
    p = DustParent()
    s = p.state(tau)
    R, f, v, uT = (s[k] for k in ("radius", "f", "radial_velocity", "T_rate"))
    g = np.diag([-f, 1/f])
    u, n = np.array([uT, v]), np.array([v/f, 1.])
    assert u @ g @ u == pytest.approx(-1., abs=3e-15)
    assert n @ g @ n == pytest.approx(1., abs=3e-15)
    assert u @ g @ n == pytest.approx(0., abs=3e-15)
    # Gamma^R_TT=f f'/2 and Gamma^R_RR=-f'/(2f).
    fp = 2*p.mass/R**2
    ar = s["radial_acceleration"]+f*fp*uT*uT/2-fp*v*v/(2*f)
    at = -fp*v/f**2+fp/f*uT*v
    assert abs(ar)+abs(at) < 2e-15
    # K_theta^theta=n(R)/R, K_tau^tau=n.a, on both sides.
    assert n[1]/R == s["K_angular"]
    assert n @ g @ np.array([at, ar]) == pytest.approx(s["K_time"], abs=2e-15)
    assert s["mass_removed"] == pytest.approx(p.mass, abs=4e-16)
    assert s["H"]**2*(1-2*s["deceleration"]) == s["parent_Lambda"]


def test_retained_collar_has_a_spacelike_constraint_satisfying_extension():
    p = DustParent()
    r = np.array([.9, 1., 1.7, 2.8, 3.5, 4.])
    s = p.initial_slice(r)
    assert np.all(s["radial_metric"] > 0)
    assert s["T_radius"][0] == 0.
    assert s["radial_metric"][-1] == pytest.approx(1., abs=3e-16)
    assert s["normal_radius_rate"][-1] == pytest.approx(p.state(0.)["radial_velocity"])
    # Direct intrinsic curvature plus extrinsic terms in the Schwarzschild
    # constraint. Finite differences are independent of the matching code.
    for r0 in (1.4, 2., 3.1):
        h = 1e-5
        def EV(r1):
            s1 = p.initial_slice(r1)
            return 1/math.sqrt(float(s1["radial_metric"])), float(s1["normal_radius_rate"])
        E, V = EV(r0)
        Em, Vm = EV(r0-h)
        Ep, Vp = EV(r0+h)
        ER, VR = (Ep-Em)/(2*h), (Vp-Vm)/(2*h)
        C = -4*E*ER/r0+2*(1-E*E)/r0**2
        assert abs(C+2*V*V/r0**2+4*V*VR/r0) < 2e-10


def test_null_clock_against_ode_and_arrival_time_derivative():
    p = DustParent()
    result = p.arrival()
    def ode(t, y):
        return [(1-2*p.mass/y[0])/p.state(t)["f"]]
    def contact(t, y):
        return y[0]-p.state(t)["radius"]
    contact.terminal = True
    contact.direction = 1.
    sol = solve_ivp(ode, (0., 10.), [2.], events=contact, method="DOP853",
                    rtol=2e-12, atol=2e-13)
    assert sol.success and len(sol.t_events[0]) == 1
    assert sol.t_events[0][0] == pytest.approx(result["reception_time"], abs=3e-10)
    assert abs(result["null_equation_residual"]) < 3e-12
    tau = result["reception_time"]
    h = 1e-4
    def receive(emission_T):
        return brentq(lambda t: p.static_time(t)-emission_T
                      -p.tortoise(p.state(t)["radius"])+p.tortoise(2.),
                      tau-.1, tau+.1, xtol=2e-13)
    clock_ratio = (receive(h)-receive(-h))/(2*h*math.sqrt(.6))
    assert clock_ratio == pytest.approx(result["one_plus_redshift"], abs=1e-9)


def test_mass_and_clock_scaling():
    a, b = DustParent(), DustParent(2., 20.)
    r1, r2 = a.arrival(), b.arrival(10.)
    assert r2["reception_time"] == pytest.approx(5*r1["reception_time"], abs=2e-12)
    assert r2["one_plus_redshift"] == pytest.approx(r1["one_plus_redshift"], abs=2e-14)
    assert b.H0 == pytest.approx(a.H0/5)


@pytest.mark.parametrize("rho,p", [(0., 0.), (2., 0.), (2., 2/3), (2., -2.)])
def test_quadratic_stress_independent_perfect_fluid_limit(rho, p):
    T = np.diag([rho, p, p, p])
    expected = np.diag([rho*rho/12]+[rho*(rho+2*p)/12]*3)
    assert quadratic_brane_stress(T) == pytest.approx(expected, abs=3e-15)
    r = projected_brane_source(T, tension=6.6, kappa5_squared=1., L5=1.)
    assert r["projection_residual"] < 2e-14


def test_full_gauss_projection_includes_momentum_and_anisotropic_stress():
    T = np.array([[.7, .1, -.2, 0.], [.1, .2, .03, 0.],
                  [-.2, .03, .3, .01], [0., 0., .01, -.1]])
    E = np.diag([.6, .1, .2, .3])
    r = projected_brane_source(T, tension=6.6, kappa5_squared=1., L5=1., weyl=E)
    assert r["projection_residual"] < 3e-15
    assert r["Lambda4"] == pytest.approx(.63, abs=2e-15)
    r0 = projected_brane_source(T, tension=6., kappa5_squared=1., L5=1.)
    assert r0["Lambda4"] == 0.
    with pytest.raises(ValueError, match="traceless"):
        projected_brane_source(T, tension=6., kappa5_squared=1., L5=1., weyl=np.eye(4))


def test_positive_induced_lambda_has_an_actual_constant_curvature_bulk():
    y = np.array([0., .3, .7, 1.])
    s = de_sitter_bulk_slice(y)
    assert s["Lambda4_L5_squared"] == pytest.approx(.63, abs=2e-15)
    # Both independent sectional curvatures are those of AdS5.
    assert s["normal_bulk_curvature"] == pytest.approx(-np.ones(4), abs=3e-15)
    assert s["tangential_bulk_curvature"] == pytest.approx(-np.ones(4), abs=3e-14)
    assert s["derivative"][0]/s["warp"][0] == -1.1
    with pytest.raises(ValueError, match="coordinate horizon"):
        de_sitter_bulk_slice(s["coordinate_horizon"])


def test_shadow_light_deflection_from_the_null_path_integral():
    for f in (0., .5, 1.):
        r = shadow_lensing(f)
        gamma = r["Psi_over_Phi"]
        # b=1 and GM_dyn=1. z=tan(theta) regularizes the whole straight ray.
        alpha = quad(lambda theta: (1+gamma)*math.cos(theta),
                     -math.pi/2, math.pi/2, epsabs=2e-13)[0]
        assert alpha/4 == pytest.approx(r["lensing_over_GR_same_dynamical_mass"], abs=1e-14)
    assert shadow_lensing(1.)["lensing_over_GR_same_dynamical_mass"] == .75
    assert shadow_lensing(0., separation_over_L5=100.)["Psi_over_Phi"] == 1.


def test_constant_bh_density_does_not_fix_the_total_physical_pressure():
    lo, hi = population_balance(1., exponent=3.), population_balance(2., exponent=3.)
    assert lo["black_hole_density"] == hi["black_hole_density"] == .1
    assert hi["reservoir_density"] == pytest.approx(.025)
    assert hi["closed_effective_w_if_separately_conserved"] == -1.
    assert hi["total_deceleration"] == .5
    # Derivatives w.r.t. log a: both open component balances cancel.
    a, h = 1.5, 1e-5
    m = population_balance(a*math.exp(-h), exponent=3.)
    p = population_balance(a*math.exp(h), exponent=3.)
    c = population_balance(a, exponent=3.)
    reservoir_rate = (p["reservoir_density"]-m["reservoir_density"])/(2*h)
    assert reservoir_rate+3*c["reservoir_density"] == pytest.approx(
        -c["transfer_Q_over_H"], abs=2e-10)
    with pytest.raises(ValueError, match="exhausted"):
        population_balance(3., exponent=3.)


def test_domains_are_explicit():
    with pytest.raises(ValueError):
        DustParent(initial_boundary=.7)
    with pytest.raises(ValueError):
        DustParent().state(-10.)
    with pytest.raises(ValueError):
        DustParent().arrival(5.)
    with pytest.raises(ValueError):
        shadow_lensing(1.1)


def test_potential_is_the_existing_signed_wall_motion():
    from recursive_horizons.nsc_atlas_nested_slice import tension_wall
    for R in (.5, 1.5, 2.):
        old = tension_wall(R, 0., 1., .4, 0., 2.)
        p = tension_potential(2*R, mu=.8, u=.25)
        assert float(p["V"]) == pytest.approx(-old["Rdot"]**2, abs=2e-14)
        assert float(p["beta_minus"]) == pytest.approx(old["beta_minus"], abs=2e-14)
        assert float(p["beta_plus"]) == pytest.approx(old["beta_plus"], abs=2e-14)
        nested = tension_wall(R, .01, 1., .4, 0., 2.)
        comparison = tension_potential(2*R, mu=.78, u=.25)
        assert float(comparison["V"])-.04/(2*R) == pytest.approx(
            -nested["Rdot"]**2, abs=2e-14)


@pytest.mark.parametrize("u", [0., .25, 1., 4.])
def test_critical_wall_obeys_both_signed_israel_conditions_and_is_unstable(u):
    s = critical_tension_wall(u)
    x, mu, bm, bp = (s[k] for k in ("x", "mu", "beta_minus", "beta_plus"))
    assert abs(s["V"])+abs(s["derivative"]) < 2e-14
    assert bm > 0 > bp
    assert bm*bm == pytest.approx(1-u*x*x, abs=2e-15)
    assert bp*bp == pytest.approx(1-2*mu/x, abs=2e-15)
    assert bm-bp == pytest.approx(x, abs=2e-15)
    assert mu/(x*x*bp)+u*x/bm == pytest.approx(-1., abs=2e-14)
    assert s["second_derivative"] < 0
    assert s["second_derivative"] == pytest.approx(
        -12*mu*mu/x**6-1.5*(1+u)**2, abs=2e-13)


def test_no_stable_equilibrium_is_an_identity_not_a_grid_inference():
    import sympy as sp
    r, A, B, C = sp.symbols("r A B C", positive=True)
    V = 1-A/r**4-B/r-C*r*r
    B_stationary = sp.solve(sp.diff(V, r), B)[0]
    reduced = sp.factor(sp.diff(V, r, 2).subs(B, B_stationary))
    assert sp.simplify(reduced+12*A/r**6+6*C) == 0
    assert critical_tension_wall(0.)["mu"] == pytest.approx(16/27)
    # mu_critical'(u) = -V_u/V_mu; signed factors prove it is negative.
    for u in (0., .25, 1., 4.):
        s = critical_tension_wall(u)
        Vu = s["x"]*s["beta_plus"]
        Vmu = -2*s["beta_minus"]/s["x"]**2
        assert -Vu/Vmu < 0
        # The sole positive stationary radius scales as mu**(1/3).
        x_max = s["x"]*(.78/s["mu"])**(1/3)
        state = tension_potential(x_max, mu=.78, u=u)
        assert abs(float(state["derivative"])) < 2e-13
        assert float(state["V"]) == pytest.approx(1-(.78/s["mu"])**(2/3), abs=2e-14)
        assert state["V"] < 0
    assert .4*2 > 16/27


def test_scale_inheritance_preserves_the_free_dimensionless_family():
    M, k, H, R = .4, 2., 1., 1.5
    first = tension_potential(k*R, mu=k*M, u=(H/k)**2)
    scale = 7.
    M, k, H, R = M*scale, k/scale, H/scale, R*scale
    second = tension_potential(k*R, mu=k*M, u=(H/k)**2)
    for key in first:
        assert first[key] == pytest.approx(second[key], abs=2e-15)
