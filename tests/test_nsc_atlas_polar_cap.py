"""Independent radial constraints, Reynolds balances and regular-axis controls."""
import time

import numpy as np
import pytest
from scipy.integrate import quad

from recursive_horizons.nsc_atlas_characteristic_initialization import Budget
from recursive_horizons.nsc_atlas_polar_cap import cap_metric, cap_rhs, cap_diagnostics


@pytest.fixture(scope="module")
def budget():
    # One aggregate allowance for this complete manufactured suite.
    result = Budget(maximum_evaluations=20000, cpu_seconds=30)
    yield result
    result.check()
    print("polar-cap manufactured evaluations:", result.evaluations,
          "process CPU:", time.process_time()-result.started)


def data(n=8, *, f=0., Pi=0., Q=0.):
    return np.linspace(0, 1, n+1), {k: np.full(n, v) for k, v in (("f", f), ("Pi", Pi), ("Q", Q))}


@pytest.mark.parametrize("alpha", [0., .2, 1e-22])
def test_exact_vacuum_geometry_and_zero_field_motion(alpha, budget):
    y, state = data()
    metric = cap_metric(y, .55, state, alpha=alpha, budget=budget)
    r = metric["r_faces"]
    np.testing.assert_array_equal(metric["m_faces"], r**3/2)
    np.testing.assert_array_equal(metric["A_faces"], 1-r*r)
    np.testing.assert_array_equal(metric["delta_faces"], np.zeros_like(r))
    np.testing.assert_array_equal(metric["n_faces"], np.zeros_like(r))
    for velocity in (-.1, 0., .1):
        rhs = cap_rhs(metric, T_x=.8, radius_x=velocity,
                      incoming_w_plus=0., wall_f=0.)
        for rate in rhs["rates"].values():
            np.testing.assert_array_equal(rate, np.zeros(8))
        assert rhs["Phi_faces"][0] == 0
        assert rhs["speed_plus"][-1] < 0 < rhs["speed_minus"][-1]


def test_zero_amplitude_retains_scalar_and_nonzero_metric_sensitivity(budget):
    y, state = data(f=.4, Pi=.3, Q=.2)
    metric = cap_metric(y, .5, state, alpha=0., budget=budget)
    np.testing.assert_array_equal(metric["potential_derivative_normalized"], 2*state["f"])
    np.testing.assert_array_equal(metric["potential_difference_over_alpha"], state["f"]**2)
    assert metric["n_faces"][-1] > 0
    assert metric["center"]["p_J"] == .4
    np.testing.assert_array_equal(metric["A_faces"], 1-metric["r_faces"]**2)
    assert metric["wall"]["w_minus"] == pytest.approx(.2)


def test_tiny_amplitude_keeps_lapse_and_mass_response_without_subtraction(budget):
    y, state = data(f=.4, Pi=.2, Q=.3)
    zero = cap_metric(y, .5, state, alpha=0., budget=budget)
    alpha = 1e-20
    tiny = cap_metric(y, .5, state, alpha=alpha, budget=budget)
    r = tiny["r_faces"]
    expected = -.4**2/12+.2**2*r*r/4+.3**2*r**4/8
    np.testing.assert_allclose(tiny["delta_faces"]/alpha, expected, rtol=2e-14, atol=1e-16)
    np.testing.assert_allclose(tiny["n_faces"], zero["n_faces"], rtol=2e-14, atol=1e-17)
    assert tiny["delta_faces"][0] != 0
    assert tiny["center"]["p_J"] == pytest.approx(.4)


def test_mass_matches_independent_full_constraint_integral(budget):
    alpha, field, momentum, radius = .2, .7, .4, .55
    y, state = data(6, f=field, Pi=momentum)
    metric = cap_metric(y, radius, state, alpha=alpha, quadrature_order=4, budget=budget)
    potential = 3*np.cosh(np.sqrt(alpha/6)*field)**4
    exponent = alpha*momentum**2/4
    # Directly integrate the original m equation, with its vacuum term;
    # this oracle never uses the implementation's subtracted n equation.
    expected = [np.exp(-exponent*r*r)*quad(
        lambda s: np.exp(exponent*s*s)*s*s*(potential/2+alpha*momentum**2/4),
        0, r, epsabs=1e-14, epsrel=1e-13)[0] for r in metric["r_faces"]]
    np.testing.assert_allclose(metric["m_faces"], expected, rtol=3e-13, atol=1e-15)
    np.testing.assert_allclose(metric["delta_faces"],
                               -np.log(np.cosh(np.sqrt(alpha/6)*field))+exponent*metric["r_faces"]**2,
                               rtol=3e-13, atol=1e-15)
    assert metric["center"]["mu"] == pytest.approx(potential/6+alpha*momentum**2/12)


def test_nonconstant_reconstructed_momentum_and_Q_match_radial_integral(budget):
    # An arbitrary manufactured constraint slice; f and Phi need not yet
    # satisfy derivative compatibility for this independent integral check.
    alpha, R, f, P0, P2, Q = .12, .55, .2, .3, .7, .4
    y, state = data(6, f=f, Q=Q)
    r = R*y
    state["Pi"] = P0+P2*3*np.diff(r**5)/(5*np.diff(r**3))
    metric = cap_metric(y, R, state, alpha=alpha, budget=budget)
    potential = 3*np.cosh(np.sqrt(alpha/6)*f)**4
    def I(z):
        return alpha*(P0*P0*z*z/4+(2*P0*P2+Q*Q)*z**4/8+P2*P2*z**6/12)
    def density(z):
        return np.exp(I(z))*z*z*(potential/2+alpha*((P0+P2*z*z)**2+Q*Q*z*z)/4)
    expected = [np.exp(-I(z))*quad(density, 0, z, epsabs=1e-14, epsrel=1e-13)[0] for z in r]
    np.testing.assert_allclose(metric["m_faces"], expected, rtol=2e-13, atol=1e-15)
    np.testing.assert_allclose(metric["delta_faces"]-metric["center"]["delta"], I(r), atol=2e-17)


@pytest.mark.parametrize("velocity", [-.08, 0., .08])
def test_constant_Q_and_Pi_have_correct_axis_and_moving_volume_rates(velocity, budget):
    R, P0, Q0 = .5, .3, .4
    y, state = data(f=0., Pi=P0, Q=Q0)
    metric = cap_metric(y, R, state, alpha=0., budget=budget)
    row = cap_rhs(metric, T_x=1., radius_x=velocity,
                  incoming_w_plus=P0+R*Q0, wall_f=0.)
    r = metric["r_faces"]
    rho2 = 3*np.diff(r**5)/(5*np.diff(r**3))
    np.testing.assert_allclose(row["rates"]["Q"], -2*P0, atol=3e-14)
    np.testing.assert_allclose(row["rates"]["Pi"], Q0*(3-5*rho2), atol=2e-14)
    assert row["Phi_flux"][0] == pytest.approx(P0)
    assert row["Pi_flux"][0] == 0
    assert row["axis_Pi_rate_limit"] == pytest.approx(3*Q0)


@pytest.mark.parametrize("n", [4, 8, 16])
def test_even_quadratic_momentum_is_consistent_in_every_axis_cell(n, budget):
    R, P0, P2 = .5, .3, .6
    y, state = data(n)
    r = R*y
    # These are actual r²-weighted averages of the analytic quadratic.
    rho2 = 3*np.diff(r**5)/(5*np.diff(r**3))
    state["Pi"] = P0+P2*rho2
    metric = cap_metric(y, R, state, alpha=0., budget=budget)
    row = cap_rhs(metric, T_x=1., radius_x=.03,
                  incoming_w_plus=P0+P2*R*R, wall_f=0.)
    np.testing.assert_allclose(metric["Pi_intercept"], P0, atol=2e-15)
    np.testing.assert_allclose(metric["Pi_rho_slope"], P2, atol=2e-13)
    # Average (B Pi)_r/r with the Q cell's r weight.  This includes the
    # first cell: an error independent of dr there would fail this test.
    rho1 = np.diff(r**4)/(2*np.diff(r*r))
    expected = 2*(P2-P0)-4*P2*rho1
    np.testing.assert_allclose(row["rates"]["Q"], expected, rtol=1e-11, atol=5e-13)
    assert row["Pi_faces"][0] == pytest.approx(P0)


def test_global_reynolds_balances_use_the_actual_boundary_fluxes(budget):
    y, state = data(7)
    state.update(f=np.linspace(.1, .3, 7), Pi=np.linspace(.2, -.1, 7), Q=np.linspace(-.2, .4, 7))
    R, Tx, Rx = .53, .7, .04
    metric = cap_metric(y, R, state, alpha=.03, budget=budget)
    row = cap_rhs(metric, T_x=Tx, radius_x=Rx, incoming_w_plus=.35, wall_f=.32)
    dy, d1, d2 = np.diff(y), np.diff(y*y)/2, np.diff(y**3)/3
    rate = row["rates"]
    pi_storage = np.dot(R**3*d2, rate["Pi"]+3*Rx/R*state["Pi"]+Tx*row["source_Pi"])
    q_storage = np.dot(R*R*d1, rate["Q"]+2*Rx/R*state["Q"])
    f_storage = np.dot(R*dy, rate["f"]+Rx/R*state["f"]-Tx*row["source_f"])
    assert pi_storage == pytest.approx(R*R*(row["Pi_flux"][-1]-row["Pi_flux"][0]), abs=2e-16)
    assert q_storage == pytest.approx(row["Phi_flux"][-1]-row["Phi_flux"][0], abs=2e-16)
    assert f_storage == pytest.approx(row["f_flux"][-1]-row["f_flux"][0], abs=2e-16)
    assert row["Pi_faces"][-1]+row["Phi_faces"][-1] == pytest.approx(.35)
    assert row["Pi_faces"][-1]-row["Phi_faces"][-1] == pytest.approx(metric["wall"]["w_minus"])
    other = cap_rhs(metric, T_x=Tx/2.3, radius_x=Rx/2.3, incoming_w_plus=.35, wall_f=.32)
    for key in ("f", "Pi", "Q"):
        np.testing.assert_allclose(other["rates"][key], rate[key]/2.3, rtol=1e-13, atol=1e-14)


@pytest.mark.parametrize("velocity", [-.1, .1])
def test_constant_value_mesh_transport_and_single_incoming_wave(velocity, budget):
    y, state = data(f=.3)
    metric = cap_metric(y, .5, state, alpha=0., budget=budget)
    row = cap_rhs(metric, T_x=1., radius_x=velocity, incoming_w_plus=0., wall_f=.3)
    np.testing.assert_allclose(row["rates"]["f"], 0., atol=5e-16)
    np.testing.assert_allclose(row["rates"]["Pi"], -.6, atol=5e-15)
    np.testing.assert_allclose(row["rates"]["Q"], 0., atol=5e-16)
    changed = cap_rhs(metric, T_x=1., radius_x=velocity, incoming_w_plus=0., wall_f=.7)
    np.testing.assert_array_equal(changed["rates"]["Pi"], row["rates"]["Pi"])
    np.testing.assert_array_equal(changed["rates"]["Q"], row["rates"]["Q"])
    if velocity < 0:
        np.testing.assert_array_equal(changed["rates"]["f"], row["rates"]["f"])
    else:
        assert changed["rates"]["f"][-1] != row["rates"]["f"][-1]


def test_independent_constraint_differences_and_derivative_compatibility(budget):
    errors = []
    for n in (4, 8):
        y, state = data(n, Q=1.)
        r = .5*y
        state["f"] = (r[:-1]**2+r[:-1]*r[1:]+r[1:]**2)/6
        metric = cap_metric(y, .5, state, alpha=0., budget=budget)
        diag = cap_diagnostics(metric)
        np.testing.assert_allclose(diag["Phi_minus_fd_f_r"], 0, atol=3e-15)
        # Even exact vacuum nodes have a nonzero trapezoidal derivative
        # diagnostic.  It is not a disguised check of the solved constraint.
        np.testing.assert_allclose(diag["mass_trapezoid_defect"], -np.diff(r)**2/4, atol=2e-16)
        errors.append(diag["max_abs_mass_trapezoid_defect"])
    assert errors[0]/errors[1] == pytest.approx(4.)
    metric["state"]["f"] += np.arange(8)*.01
    assert cap_diagnostics(metric)["max_abs_Phi_minus_fd_f_r"] > .1


@pytest.mark.parametrize("change", ["negative_alpha", "bad_radius", "nonuniform", "missing_field", "bad_shape", "nan", "order", "boolean"])
def test_invalid_metric_inputs_are_rejected(change, budget):
    y, state = data(4)
    kwargs = dict(alpha=0., budget=budget)
    R = .5
    if change == "negative_alpha": kwargs["alpha"] = -.01
    if change == "bad_radius": R = False
    if change == "nonuniform": y[1] += .01
    if change == "missing_field": del state["f"]
    if change == "bad_shape": state["Q"] = np.zeros((2, 2))
    if change == "nan": state["Pi"][1] = np.nan
    if change == "order": kwargs["quadrature_order"] = True
    if change == "boolean": state["Q"] = [True, 0., 0., 0.]
    with pytest.raises(ValueError):
        cap_metric(y, R, state, **kwargs)


def test_trapped_slice_timelike_guard_and_budget_failures(budget):
    y, state = data(4)
    with pytest.raises(ValueError, match="A,B>0"):
        cap_metric(y, 1.1, state, alpha=0., budget=budget)
    metric = cap_metric(y, .5, state, alpha=0., budget=budget)
    with pytest.raises(ValueError, match="timelike"):
        cap_rhs(metric, T_x=1., radius_x=.75, incoming_w_plus=0., wall_f=0.)
    with pytest.raises(ValueError, match="same cumulative"):
        cap_rhs(metric, T_x=1., radius_x=0., incoming_w_plus=0., wall_f=0., budget=Budget())
    with pytest.raises(RuntimeError, match="budget exhausted"):
        cap_metric(y, .5, state, alpha=0., budget=Budget(maximum_evaluations=1, cpu_seconds=30))
    expired = Budget(maximum_evaluations=20000, cpu_seconds=30, started=time.process_time()-31)
    with pytest.raises(RuntimeError, match="budget exhausted"):
        cap_metric(y, .5, state, alpha=0., budget=expired)


def test_exact_charge_counts_and_caller_inputs_are_preserved(budget):
    y, state = data(4, f=.2, Pi=.3)
    saved = {key: value.copy() for key, value in state.items()}
    before = budget.evaluations
    metric = cap_metric(y, .5, state, alpha=.1, quadrature_order=2, budget=budget)
    assert budget.evaluations-before == metric["charged_evaluations"] == 36
    before = budget.evaluations
    row = cap_rhs(metric, T_x=1., radius_x=0., incoming_w_plus=.3, wall_f=.2)
    assert budget.evaluations-before == row["charged_evaluations"] == 13
    for key in saved:
        np.testing.assert_array_equal(state[key], saved[key])
    metric["state"]["f"][0] = 99
    assert state["f"][0] == .2
