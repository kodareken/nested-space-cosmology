"""Physical junction, independent wave and causal controls for the nested body."""
import math

import numpy as np
import pytest
from scipy.integrate import solve_ivp
import sympy as sp

from recursive_horizons.nsc_atlas_gravitational_response import (
    Budget, Region, analytic_cap_state, axial_metric_data, cap_state,
    dirichlet_kernel, interface_transfer, nested_body, response, stiff_inside, transfer,
)


def test_stiff_scalar_radial_law_and_exact_one_shell_stability():
    r = sp.symbols("r", positive=True)
    R = sp.Rational(5, 2)
    a = (sp.sqrt(5)+sp.sqrt(21))/8
    H2 = (1-a*a)/R**2
    k0 = (a-1/sp.sqrt(5))/R
    # Canonical homogeneous surface scalar: d(R^2*theta_dot)/d tau=0.
    sigma = sp.symbols("sigma", positive=True)*(R/r)**4
    assert sp.simplify(sp.diff(sigma*r*r, r)+sigma*2*r) == 0
    w = k0*R*(R/r)**3
    V = 1-H2*r*r-((2/r-H2*r*r+w*w)/(2*w))**2
    assert sp.simplify(V.subs(r, R)) == 0
    assert sp.simplify(sp.diff(V, r).subs(r, R)) == 0
    assert sp.simplify(sp.diff(V, r, 2).subs(r, R)
                       -(4*sp.sqrt(105)-28)/25) == 0
    # The first compact trial fails rather than being accepted by appearance.
    trial = stiff_inside(Region(1., 0.), inside_mass=0., radius=2.2)
    assert trial.radial_potential(2.2)[2] < 0


def test_two_level_background_solves_all_junctions_and_force_balance():
    body = nested_body()
    for shell in (body.inner, body.outer):
        s = shell.record()
        assert s["sigma"] > 0
        assert s["pressure"] == pytest.approx(s["sigma"], abs=1e-15)
        for key in ("V", "Vprime", "induced_clock_gap", "israel_gap", "force_gap"):
            assert abs(s[key]) < 2e-13
        assert s["Vsecond"] > 0
        assert s["q"] == pytest.approx(4*math.pi*s["sigma"]*s["clock_lapse"], abs=1e-15)
    assert body.core.f(body.inner.radius) > 0
    # The Kottler lapse is concave, so positive endpoint values suffice.
    assert body.parent.f(body.inner.radius) > 0
    assert body.parent.f(body.outer.radius) > 0
    assert body.parent.fpp(1.) < 0


@pytest.mark.parametrize("shell_index", [0, 1])
def test_all_axial_metric_and_surface_scalar_junctions(shell_index):
    body = nested_body()
    shell = (body.inner, body.outer)[shell_index]
    R, omega = shell.radius, 1.3
    minus_state = np.array([.7+.2j, -.4+.8j])
    plus_state = interface_transfer(shell) @ minus_state
    minus = axial_metric_data(shell.inside, R, minus_state, omega)
    plus = axial_metric_data(shell.outside, R, plus_state, omega)
    a, b = math.sqrt(shell.inside.f(R)), math.sqrt(shell.outside.f(R))
    assert plus["h0"] == pytest.approx(minus["h0"], abs=2e-15)
    assert b*plus["h1"] == pytest.approx(a*minus["h1"], abs=2e-15)
    # The tA Israel condition for the irrotational surface scalar has S_tA=p h_tA.
    # K_tA = sqrt(f)*(h0'+i*omega*h1)/2; the remaining condition is independent
    # of the two matching assertions and uses the bulk Einstein wave equation.
    curvature_jump = (b*(plus["h0p"]+1j*omega*plus["h1"])
                      -a*(minus["h0p"]+1j*omega*minus["h1"]))/2
    assert curvature_jump == pytest.approx((b-a)*minus["h0"]/R, abs=3e-14)


@pytest.mark.parametrize("omega", [.4, 2.])
def test_regular_cap_against_independent_hypergeometric_solution(omega):
    body = nested_body()
    numerical = cap_state(body, omega, 2, Budget(), rtol=1e-11, fraction=.001)
    exact = analytic_cap_state(body.core, body.inner.radius, omega, 2)
    assert numerical[1]/numerical[0] == pytest.approx(exact[1]/exact[0], rel=2e-11)


def test_composition_is_symplectic_and_agrees_with_an_unsplit_bulk_propagation():
    body, omega = nested_body(), .8
    R0, R1, R2 = body.inner.radius, 1., body.outer.radius
    budget = Budget()
    # The low-frequency barrier makes this fundamental matrix ill conditioned;
    # resolve its determinant more tightly than a normalized response requires.
    one = transfer(body.parent, R0, R2, omega, 2, budget, rtol=2e-13)
    left = transfer(body.parent, R0, R1, omega, 2, budget, rtol=2e-13)
    right = transfer(body.parent, R1, R2, omega, 2, budget, rtol=2e-13)
    assert right @ left == pytest.approx(one, rel=3e-11, abs=2e-10)
    assert np.linalg.det(one) == pytest.approx(1., abs=1e-10)
    J = np.array([[0., 1.], [-1., 0.]])
    for shell in (body.inner, body.outer):
        S = interface_transfer(shell)
        assert S.T @ J @ S == pytest.approx(J, abs=1e-15)


@pytest.mark.parametrize("omega", [.4, 2.])
def test_regular_and_ingoing_responses_preserve_physical_flux(omega):
    body = nested_body()
    cap = response(body, omega)
    bh = response(body, omega, termination="horizon")
    assert cap["reflected_power"] == pytest.approx(1., abs=2e-14)
    assert bh["absorbed_power"] > 0
    assert abs(bh["energy_balance_gap"]) < 3e-10
    assert bh["flux_relative_gap"] < 2e-8
    assert abs(cap["reflection"]-bh["reflection"]) > 1e-7
    if omega == 2.:
        assert bh["absorbed_power"] > .4


def test_common_clock_mass_scaling_preserves_dimensionless_response():
    first, second = nested_body(), nested_body(3.)
    for t in ("cap", "horizon"):
        a = response(first, 2., termination=t)
        b = response(second, 2/3, termination=t)
        assert a["reflection"] == pytest.approx(b["reflection"], abs=2e-10)
    assert second.causal_time() == pytest.approx(3*first.causal_time(), rel=2e-14)


def test_actual_gravitational_region_gluing_is_the_schur_complement():
    body, omega, budget = nested_body(), 1.1, Budget()
    left = transfer(body.parent, body.inner.radius, body.outer.radius,
                    omega, 2, budget, rtol=2e-13)
    right = transfer(body.exterior, body.outer.radius, body.observer,
                     omega, 2, budget, rtol=2e-13)
    K1, K2 = dirichlet_kernel(left), dirichlet_kernel(right)
    full = np.zeros((3, 3), dtype=complex)
    full[:2, :2] += K1
    full[1:, 1:] += K2
    full[1, 1] += body.outer.jump_strength
    keep = [0, 2]
    reduced = full[np.ix_(keep, keep)]-np.outer(full[keep, 1], full[1, keep])/full[1, 1]
    composed = right @ interface_transfer(body.outer) @ left
    assert reduced == pytest.approx(dirichlet_kernel(composed), rel=3e-11, abs=2e-10)


def test_return_time_against_an_independent_null_ode():
    body = nested_body()
    # Integrate coordinate time along an ingoing principal characteristic.
    elapsed = 0.
    for region, start, stop in [(body.exterior, body.observer, body.outer.radius),
                                (body.parent, body.outer.radius, body.inner.radius)]:
        sol = solve_ivp(lambda r, y: [-1/region.speed(r)], (start, stop), [0.],
                        rtol=2e-12, atol=2e-13, method="DOP853")
        assert sol.success
        elapsed += sol.y[0, -1]
    proper = 2*math.sqrt(body.exterior.f(body.observer))*elapsed
    assert proper == pytest.approx(body.causal_time(), abs=2e-10)
    assert body.causal_time(center=True) > proper


def test_domain_and_budget_failures_are_local():
    with pytest.raises(ValueError):
        nested_body(0.)
    with pytest.raises(ValueError):
        stiff_inside(Region(1., 0.), inside_mass=0., radius=1.9)
    with pytest.raises(ValueError):
        response(nested_body(), .4, termination="arbitrary")
    with pytest.raises(RuntimeError, match="allowance"):
        transfer(Region(1., 0.), 3., 4., .4, 2, Budget(limit=1))
