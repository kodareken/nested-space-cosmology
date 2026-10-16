"""Regular reduced variables, independent constraints and incoming-wave controls."""
import time

import numpy as np
import pytest
from scipy.integrate import quad

from recursive_horizons.nsc_atlas_characteristic_initialization import Budget
from recursive_horizons.nsc_atlas_nonlinear_boundary import hybrid_wall_boundary
from recursive_horizons.nsc_atlas_reduced_cap import (
    _fields, _gamma, _potential_reaction_quotient,
    reduced_cap_metric, reduced_cap_rhs, reduced_cap_profile, reduced_cap_diagnostics,
)


@pytest.fixture(scope="module")
def budget():
    result = Budget(maximum_evaluations=20000, cpu_seconds=30)
    yield result
    result.check()
    print("reduced-cap manufactured evaluations:", result.evaluations,
          "process CPU:", time.process_time()-result.started)


def polynomial(n=6, R=.5, f0=0., f2=0., P0=0., P2=0.):
    y = np.linspace(0, 1, n+1)
    r = R*y
    return y, {"Ip": P0*np.diff(r*r)/2+P2*np.diff(r**4)/4,
               "Iq": f0*np.diff(r)+f2*np.diff(r**3)}


@pytest.mark.parametrize("alpha", [0., .2, 1e-22])
def test_exact_vacuum_and_zero_motion_field(alpha, budget):
    y, state = polynomial()
    metric = reduced_cap_metric(y, .5, state, alpha=alpha, budget=budget)
    r = metric["r_faces"]
    np.testing.assert_array_equal(metric["m_faces"], r**3/2)
    np.testing.assert_array_equal(metric["A_faces"], 1-r*r)
    np.testing.assert_array_equal(metric["delta_faces"], np.zeros_like(r))
    assert metric["wall"]["c"] == metric["wall"]["chi"] == metric["center"]["f"] == 0
    for speed in (-.15, 0., .15):
        row = reduced_cap_rhs(metric, T_x=.8, radius_x=speed, incoming_w_plus=0.)
        for rates in row["rates"].values():
            np.testing.assert_array_equal(rates, np.zeros(6))
        assert row["q_flux"][0] == row["wall_chi_x"] == row["wall_c_x"] == 0


def test_axis_parity_true_centre_and_integral_preservation(budget):
    R, f0, f2, P0, P2 = .5, .7, .4, .3, -.2
    y, state = polynomial(R=R, f0=f0, f2=f2, P0=P0, P2=P2)
    metric = reduced_cap_metric(y, R, state, alpha=.04, budget=budget)
    query = np.array([0., 1e-8, R/12])
    row = reduced_cap_profile(metric, query)
    np.testing.assert_allclose(row["f"], f0+f2*query**2, atol=1e-15)
    np.testing.assert_allclose(row["Pi"], P0+P2*query**2, atol=1e-15)
    np.testing.assert_allclose(row["Phi"], 2*f2*query, atol=1e-15)
    np.testing.assert_allclose(row["q"], f0+3*f2*query**2, atol=1e-15)
    assert row["p"][0] == row["chi"][0] == row["Phi"][0] == 0
    assert metric["center"]["f"] == pytest.approx(f0, abs=1e-15)
    assert metric["f_cell_average"][0] == pytest.approx(f0+f2*(R/6)**2/3, abs=1e-15)
    assert metric["center"]["delta"] == pytest.approx(-np.log(np.cosh(np.sqrt(.04/6)*f0)), abs=1e-16)
    diagnostics = reduced_cap_diagnostics(metric)
    assert diagnostics["max_abs_constructed_p_integral_defect"] < 1e-17
    assert diagnostics["max_abs_constructed_q_integral_defect"] < 1e-16
    assert metric["wall"]["chi"] == pytest.approx(R*(f0+f2*R*R), abs=1e-15)


def test_profile_primitive_and_plain_cell_average_are_same_function(budget):
    y, state = polynomial(f0=.2, f2=.9, P0=.3)
    metric = reduced_cap_metric(y, .5, state, alpha=0., budget=budget)
    rec = metric["reconstruction"]
    for i, (a, b) in enumerate(zip(metric["r_faces"][:-1], metric["r_faces"][1:])):
        expected = quad(lambda r: float(_fields(rec, i, r)["f"]), a, b, epsabs=1e-13)[0]/(b-a)
        assert metric["f_cell_average"][i] == pytest.approx(expected, abs=2e-15)
        qi = quad(lambda r: float(_fields(rec, i, r)["q"]), a, b, epsabs=1e-13)[0]
        assert qi == pytest.approx(state["Iq"][i], abs=2e-16)
        # Analytic derivative of the independently written primitive.
        if i:
            r = (a+b)/2
            data = _fields(rec, i, r)
            assert (data["q"]-data["f"])/r == pytest.approx(data["Phi"], abs=2e-14)


def test_gamma_primitive_matches_independent_positive_integrand(budget):
    y, state = polynomial(n=5, f0=.4, f2=-.6, P0=.3, P2=.8)
    metric = reduced_cap_metric(y, .5, state, alpha=.12, budget=budget)
    rec = metric["reconstruction"]
    for i, (a, b) in enumerate(zip(metric["r_faces"][:-1], metric["r_faces"][1:])):
        def density(r):
            data = _fields(rec, i, r)
            return .12*r*(data["Pi"]**2+data["Phi"]**2)/2
        for end in (a+(b-a)*1e-7, (a+b)/2, b):
            expected = quad(density, a, end, epsabs=1e-14)[0]
            assert _gamma(rec, i, end, .12) == pytest.approx(expected, rel=3e-11, abs=5e-18)


def test_original_mass_constraint_independent_full_integral(budget):
    alpha, f, Pi, R = .12, .6, .4, .5
    y, state = polynomial(n=6, R=R, f0=f, P0=Pi)
    metric = reduced_cap_metric(y, R, state, alpha=alpha, budget=budget)
    potential = 3*np.cosh(np.sqrt(alpha/6)*f)**4
    a = alpha*Pi*Pi/4
    expected = [np.exp(-a*r*r)*quad(
        lambda z: np.exp(a*z*z)*z*z*(potential/2+alpha*Pi*Pi/4), 0, r,
        epsabs=1e-14)[0] for r in metric["r_faces"]]
    np.testing.assert_allclose(metric["m_faces"], expected, rtol=3e-13, atol=1e-15)
    np.testing.assert_allclose(metric["delta_faces"],
                               -np.log(np.cosh(np.sqrt(alpha/6)*f))+a*metric["r_faces"]**2,
                               rtol=2e-13, atol=1e-16)


def test_tiny_alpha_response_and_stable_reaction_limit(budget):
    y, state = polynomial(f0=.4, P0=.3)
    zero = reduced_cap_metric(y, .5, state, alpha=0., budget=budget)
    tiny = reduced_cap_metric(y, .5, state, alpha=1e-22, budget=budget)
    expected = -.4**2/12+.3**2*tiny["r_faces"]**2/4
    np.testing.assert_allclose(tiny["delta_faces"]/1e-22, expected, atol=1e-16)
    np.testing.assert_allclose(tiny["n_faces"], zero["n_faces"], atol=1e-17)
    values = np.array([0., -.7, .4, 1.2])
    # sinhc(z)cosh(z)^3=1+(5/3)z^2+O(z^4).
    np.testing.assert_allclose(_potential_reaction_quotient(values, 0.), 5*values**3/9, atol=2e-16)
    np.testing.assert_allclose(_potential_reaction_quotient(values, 1e-22), 5*values**3/9, atol=2e-16)
    row = reduced_cap_rhs(tiny, T_x=1., radius_x=.1, incoming_w_plus=.2)
    assert np.any(row["reaction_nodes"] != 0)
    assert np.all(np.isfinite(row["reaction_nodes"]/1e-22))


def test_finite_alpha_reaction_matches_original_canonical_operator(budget):
    alpha = .15
    y, state = polynomial(f0=.4, f2=.2, P0=.3, P2=.1)
    metric = reduced_cap_metric(y, .5, state, alpha=alpha, budget=budget)
    row = reduced_cap_rhs(metric, T_x=1., radius_x=.1, incoming_w_plus=.2)
    r, f = metric["r_nodes"], metric["f_nodes"]
    z = np.sqrt(alpha/6)*f
    potential = 3*np.cosh(z)**4
    U = 2*f*np.sinh(z)/z*np.cosh(z)**3
    # B_r follows the original full mass/lapse constraints.  This oracle
    # retains their vacuum terms instead of using the factored reaction.
    expected = np.exp(metric["delta_nodes"])*(
        (2*metric["m_nodes"]/r**2-r*potential)*f+r*U)
    np.testing.assert_allclose(row["reaction_nodes"], expected, rtol=2e-11, atol=4e-16)
    assert np.max(abs(row["reaction_nodes"])) > 1e-5


@pytest.mark.parametrize("velocity", [-.1, 0., .1])
def test_regular_linear_p_constant_q_exact_moving_flux_and_center_rate(velocity, budget):
    R, f0, P0 = .5, .4, .3
    y, state = polynomial(R=R, f0=f0, P0=P0)
    metric = reduced_cap_metric(y, R, state, alpha=0., budget=budget)
    row = reduced_cap_rhs(metric, T_x=1., radius_x=velocity, incoming_w_plus=P0)
    r, A = R*y, 1-(R*y)**2
    np.testing.assert_allclose(row["p_flux"], A*f0+velocity*y*r*P0, atol=2e-16)
    np.testing.assert_allclose(row["q_flux"], A*r*P0+velocity*y*f0, atol=2e-16)
    expected_q_rate = np.diff(A*r*P0+velocity*y*f0)
    np.testing.assert_allclose(row["rates"]["Iq"], expected_q_rate, atol=2e-16)
    # The true centre Q0=(7*qbar0-qbar1)/6, and qbar includes moving volume.
    qbar_rate = row["rates"]["Iq"]/np.diff(r)-velocity/R*state["Iq"]/np.diff(r)
    center_rate = (7*qbar_rate[0]-qbar_rate[1])/6
    assert center_rate == pytest.approx(P0, abs=2e-15)


@pytest.mark.parametrize("alpha", [0., .03])
def test_hybrid_wall_chain_uses_same_conserved_scalar(alpha, budget):
    R, velocity = .5, .2
    y, state = polynomial(R=R, f0=.3, f2=.2, P0=.4)
    metric = reduced_cap_metric(y, R, state, alpha=alpha, budget=budget)
    c, A = metric["wall"]["c"], metric["wall"]["A"]
    D = 1+alpha*(c/R)**2/6
    H = velocity+np.sqrt(velocity*velocity+A)
    Rp = H-2*D**1.5*R
    wall = hybrid_wall_boundary({"r": R, "c": c, "eta_plus": 0.},
        {"r_V": Rp/2, "c_V": .15, "log_C_V": -.1},
        {key: metric["wall"][key] for key in ("A", "delta", "w_minus")},
        {"C_wall": 1., "log_C_s": 0.}, alpha=alpha)
    row = reduced_cap_rhs(metric, T_x=float(wall["T_s"]), radius_x=float(wall["r_s"]),
                          incoming_w_plus=float(wall["w_plus"]))
    assert row["wall_chi_x"] == pytest.approx(velocity*metric["wall"]["f"]+R*wall["f_s"], abs=2e-16)
    assert row["wall_c_x"] == pytest.approx(wall["c_s"], abs=3e-16)
    assert abs(row["constructed_wall_flux_defect"]) < 2e-16
    assert row["q_flux"][0] == 0
    assert row["speed_plus"][-1] < 0 < row["speed_minus"][-1]


def bump(zeta):
    z = (np.asarray(zeta)-.65)/.10
    active = abs(z) < 1
    safe = np.where(active, z, 0.)
    return np.where(active, .01*(1-safe*safe)**6, 0.)


def bump_prime(zeta):
    z = (np.asarray(zeta)-.65)/.10
    active = abs(z) < 1
    safe = np.where(active, z, 0.)
    return np.where(active, -1.2*safe*(1-safe*safe)**5, 0.)


def incoming_state(n, T):
    R = .45+.5*T
    y = np.linspace(0, 1, n+1)
    chi = bump(T+np.arctanh(R*y))
    integral = np.diff(chi)
    return y, R, {"Ip": integral.copy(), "Iq": integral.copy()}


def test_nonconstant_incoming_stage_has_zero_outgoing_and_correct_moving_rate(budget):
    errors = []
    for n in (16, 32):
        T = .15
        y, R, state = incoming_state(n, T)
        metric = reduced_cap_metric(y, R, state, alpha=0., quadrature_order=2, budget=budget)
        incoming = 2*bump_prime(T+np.arctanh(R))/(1-R*R)
        row = reduced_cap_rhs(metric, T_x=1., radius_x=.5,
                              incoming_w_plus=(incoming-metric["wall"]["f"])/R)
        np.testing.assert_array_equal(row["rates"]["Ip"], row["rates"]["Iq"])
        np.testing.assert_array_equal(row["reaction_nodes"], np.zeros_like(row["reaction_nodes"]))
        np.testing.assert_array_equal(row["a_minus_faces"], np.zeros(n+1))
        r = R*y
        exact_flux = bump_prime(T+np.arctanh(r))*(1+.5*y/(1-r*r))
        errors.append(float(np.sum(abs(row["rates"]["Iq"]-np.diff(exact_flux)))))
        assert row["wall_chi_x"] == pytest.approx(exact_flux[-1], abs=1e-16)
    assert errors[1] < errors[0]
    print("incoming stage L1 integral-rate errors:", errors)


def test_short_incoming_propagation_without_resetting_wall_value(budget):
    # Actual SSPRK2 evolution of a supplied analytic manufactured wave,
    # not the original prepared source or a finite-amplitude physical run.
    # Only a_plus=2*F'/A is prescribed.  Its w_plus conversion uses the
    # actual evolved wall f, so no exact wall value is imposed as extra data.
    errors = []
    for n in (16, 32):
        start, stop, steps = .14, .16, 8
        y, _, state = incoming_state(n, start)
        dt = (stop-start)/steps
        def rate(T, values):
            R = .45+.5*T
            metric = reduced_cap_metric(y, R, values, alpha=0., quadrature_order=2, budget=budget)
            incoming = 2*bump_prime(T+np.arctanh(R))/(1-R*R)
            result = reduced_cap_rhs(metric, T_x=1., radius_x=.5,
                incoming_w_plus=(incoming-metric["wall"]["f"])/R)
            assert dt*result["characteristic_courant_rate"] < .4
            assert result["q_flux"][0] == 0
            assert abs(result["constructed_wall_flux_defect"]) < 1e-16
            return result["rates"]
        for step in range(steps):
            T = start+step*dt
            first = rate(T, state)
            trial = {key: state[key]+dt*first[key] for key in state}
            second = rate(T+dt, trial)
            state = {key: (state[key]+trial[key]+dt*second[key])/2 for key in state}
        # The first axis cell enforces regular parity, rather than forcing
        # p=q there.  A tiny upwind precursor can reach that cell before
        # the exact compact wave; retain its actual outgoing defect.
        outgoing_defect = float(np.max(abs(state["Ip"]-state["Iq"])))
        assert outgoing_defect < 1e-16
        R = .45+.5*stop
        chi = np.r_[0., np.cumsum(state["Iq"])]
        expected = bump(stop+np.arctanh(R*y))
        errors.append(float(np.max(abs(chi-expected))))
        times = np.linspace(start, stop, steps+1)
        radii = .45+.5*times
        exact_wall_flux = bump_prime(times+np.arctanh(radii))*(1+.5/(1-radii*radii))
        budget.check(len(times))
        # SSPRK2 applies the trapezoidal time rule to this prescribed flux.
        # Compare that independent quadrature, without resetting chi_W.
        trapezoidal_wall = bump(start+np.arctanh(radii[0]))+np.trapezoid(exact_wall_flux, times)
        assert chi[-1] == pytest.approx(trapezoidal_wall, abs=2e-16)
        print("incoming propagated outgoing integral defect/wall time error:",
              n, outgoing_defect, float(abs(chi[-1]-expected[-1])))
    assert errors[1] < errors[0]
    print("incoming propagated primitive max errors:", errors)


@pytest.mark.parametrize("change", ["state", "negative_alpha", "grid", "radius", "order", "boolean"])
def test_invalid_inputs(change, budget):
    y, state = polynomial(n=3)
    args = {"alpha": .1, "budget": budget}
    R = .5
    if change == "state":
        state["f"] = np.zeros(3)
    elif change == "negative_alpha":
        args["alpha"] = -.1
    elif change == "grid":
        y[1] += .01
    elif change == "radius":
        R = 0.
    elif change == "order":
        args["quadrature_order"] = 3
    else:
        state["Ip"] = [False, 0., 0.]
    with pytest.raises(ValueError):
        reduced_cap_metric(y, R, state, **args)


def test_stage_and_query_guards_and_budget(budget):
    y, state = polynomial(n=3)
    metric = reduced_cap_metric(y, .5, state, alpha=0., budget=budget)
    with pytest.raises(ValueError, match="timelike"):
        reduced_cap_rhs(metric, T_x=1., radius_x=.75, incoming_w_plus=0.)
    with pytest.raises(ValueError, match="share"):
        reduced_cap_rhs(metric, T_x=1., radius_x=0., incoming_w_plus=0., budget=Budget())
    with pytest.raises(ValueError, match="finite cap"):
        reduced_cap_profile(metric, [-1e-15])
    with pytest.raises(ValueError, match="finite cap"):
        reduced_cap_profile(metric, np.zeros(514))
    with pytest.raises(RuntimeError, match="budget"):
        reduced_cap_metric(y, .5, state, alpha=0., budget=Budget(maximum_evaluations=0))
    with pytest.raises(ValueError, match="A,B"):
        reduced_cap_metric(y, 1.1, state, alpha=0., budget=budget)
