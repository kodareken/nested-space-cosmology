"""Bounded derivative/inflow controls, not a coupled physical evolution."""
import time

import numpy as np
import pytest

from recursive_horizons.nsc_atlas_characteristic_initialization import Budget
from recursive_horizons.nsc_atlas_departure_gradient import (
    first_order_departure_gradient as first_order,
    quadratic_departure_gradient as quadratic,
)


@pytest.fixture(scope="module")
def budget():
    shared = Budget(maximum_evaluations=1_500_000, cpu_seconds=45.)
    yield shared
    shared.check()
    print("departure-gradient aggregate", shared.evaluations,
          "points", time.process_time()-shared.started, "process CPU seconds")


def test_polynomials_all_rows_channels_and_input_preservation(budget):
    eta = np.linspace(0., 1., 17)
    values = np.stack((np.full_like(eta, .3), -.2+2*eta, .7-eta+3*eta**2))
    before = values.copy(); before_eta = eta.copy()
    start = budget.evaluations
    result = quadratic(eta, values, budget=budget)
    assert budget.evaluations-start == values.size
    np.testing.assert_allclose(result, np.stack((0*eta, 0*eta+2, -1+6*eta)), atol=3e-14)
    np.testing.assert_array_equal(values, before)
    np.testing.assert_array_equal(eta, before_eta)
    assert not np.shares_memory(result, values)
    np.testing.assert_array_equal(quadratic(eta, 0*eta, budget=budget), 0*eta)
    np.testing.assert_array_equal(quadratic(eta, values[2], budget=budget), result[2])


def test_inflow_orientation_and_independent_wall_derivative(budget):
    eta = np.linspace(0., 1., 9); h = eta[1]
    # One changed incoming wall value reaches the penultimate and prior row;
    # there is no replacement by a requested/known wall derivative.
    wall = np.zeros(9); wall[-1] = 1.
    expected = np.zeros(9); expected[-3:] = (-1/(2*h), 1/(2*h), 3/(2*h))
    np.testing.assert_array_equal(quadratic(eta, wall, budget=budget), expected)
    # Outflow/source data remain data: the helper does not overwrite them.
    source = np.zeros(9); source[0] = 2.
    expected = np.zeros(9); expected[0] = -3/h
    np.testing.assert_array_equal(quadratic(eta, source, budget=budget), expected)
    assert source[0] == 2. and wall[-1] == 1.


def test_first_order_is_the_frozen_secant_comparator(budget):
    eta = np.linspace(0., 1., 17)
    values = np.stack((np.sin(2*eta), eta**2, np.exp(eta)))
    expected = np.stack([np.r_[np.diff(z), z[-1]-z[-2]]/np.diff(eta)[0] for z in values])
    np.testing.assert_array_equal(first_order(eta, values, budget=budget), expected)
    np.testing.assert_array_equal(first_order([0., .5, 1.], [1., 2., 3.], budget=budget), [2., 2., 2.])
    with pytest.raises(ValueError, match="4..512"):
        quadratic([0., .5, 1.], [1., 2., 3.], budget=budget)


@pytest.mark.parametrize("eta,values", [
    ([0., .2, .5, .75, 1.], [0.]*5),
    ([.01, .25, .5, .75, 1.], [0.]*5),
    ([0., .25, .5, .75, .99], [0.]*5),
    ([0., .25, .5, .5, 1.], [0.]*5),
    ([0., .25, .5, .75, 1.], [0.]*4),
    ([0., .25, .5, .75, 1.], [[0.]*5]*4),
    ([0., .25, .5, .75, 1.], np.zeros((1, 1, 5))),
    ([0., .25, .5, .75, 1.], [0., True, 0., 0., 0.]),
    ([False, .25, .5, .75, 1.], [0.]*5),
    ([0., .25, .5, .75, 1.], [0., np.nan, 0., 0., 0.]),
    ([0., .25, .5, .75, 1.], [0., 1j, 0., 0., 0.]),
    (np.linspace(0., 1., 514), np.zeros(514)),
])
def test_malformed_inputs(eta, values, budget):
    before = budget.evaluations
    with pytest.raises(ValueError):
        quadratic(eta, values, budget=budget)
    assert budget.evaluations == before


def test_budget_live_shared_limits_and_overflow(budget):
    eta = np.linspace(0., 1., 5)
    for overrides in ({"evaluations": -1}, {"maximum_evaluations": 0},
                      {"cpu_seconds": 0.}, {"started": time.process_time()+100.},
                      {"cpu_seconds": complex(1., 1.)}):
        bad = Budget()
        for key, value in overrides.items():
            setattr(bad, key, value)
        with pytest.raises(ValueError):
            quadratic(eta, eta, budget=bad)
    with pytest.raises(ValueError, match="budget"):
        quadratic(eta, eta, budget=None)
    # Charge the attempted tiny-budget work to the enclosing control allowance.
    budget.check(5)
    tiny = Budget(maximum_evaluations=4, cpu_seconds=45.)
    with pytest.raises(RuntimeError, match="exhausted"):
        quadratic(eta, eta, budget=tiny)
    assert tiny.evaluations == 5
    with pytest.raises(ValueError, match="nonfinite departure"):
        quadratic(eta, np.full(5, np.finfo(float).max), budget=budget)


def _ssprk3(state, t, dt, rhs, impose, budget):
    """Independent nonautonomous Shu--Osher updates at t,t+dt,t+dt/2."""
    state = impose(state.copy(), t)
    k0 = rhs(t, state)
    budget.check(state.size)
    predictor = impose(state+dt*k0, t+dt)
    k1 = rhs(t+dt, predictor)
    budget.check(state.size)
    midpoint = impose(.75*state+.25*(predictor+dt*k1), t+dt/2)
    k2 = rhs(t+dt/2, midpoint)
    budget.check(state.size)
    accepted = impose(state/3+2*(midpoint+dt*k2)/3, t+dt)
    return accepted, (predictor, midpoint, accepted)


def _shape(s):
    return np.sin(2*np.pi*s)+.2*np.cos(3*np.pi*s)


def _shape_prime(s):
    return 2*np.pi*np.cos(2*np.pi*s)-.6*np.pi*np.sin(3*np.pi*s)


def _smooth_run(cells, method, kind, budget):
    eta = np.linspace(0., 1., cells+1)
    final = .25; steps = 4*cells; dt = final/steps
    def exact(t, points):
        budget.check(np.size(points))
        if kind == "constant":
            return _shape(points+.6*t)-_shape(.6*t)
        return _shape(points*np.exp(.6*t+.1*t*t))-_shape(0.)
    def impose(z, t):
        z[[0, -1]] = exact(t, np.array([0., 1.]))
        assert z[0] == 0.
        return z
    maximum_courant = 0.
    def rhs(t, z):
        nonlocal maximum_courant
        gradient = method(eta, z, budget=budget)
        budget.check(z.size)
        if kind == "constant":
            speed = np.full_like(eta, .6)
            result = speed*gradient-.6*_shape_prime(.6*t)
        else:
            speed = eta*(.6+.2*t)
            result = speed*gradient
        maximum_courant = max(maximum_courant, float(np.max(speed))*dt*cells)
        result[[0, -1]] = 0.  # boundary values, not an evolved boundary ODE
        return result
    state = exact(0., eta)
    for step in range(steps):
        state, _ = _ssprk3(state, step*dt, dt, rhs, impose, budget)
        assert np.all(np.isfinite(state))
    truth = exact(final, eta)
    budget.check(3*state.size)
    difference = state-truth
    return {"linf": float(np.max(abs(difference[1:-1]))),
            "l2": float(np.sqrt(np.trapezoid(difference**2, eta))),
            "wall_jet_error": abs(float(method(eta, state, budget=budget)[-1])-
                float((.6*0+1.)*_shape_prime(1.+.6*final) if kind == "constant" else
                      np.exp(.6*final+.1*final*final)*_shape_prime(np.exp(.6*final+.1*final*final)))),
            "courant": maximum_courant}


@pytest.mark.parametrize("kind", ["constant", "eta_flow"])
def test_manufactured_inflow_convergence_and_first_order_comparison(kind, budget):
    results = {name: [_smooth_run(n, method, kind, budget) for n in (16, 32, 64)]
               for name, method in (("first", first_order), ("quadratic", quadratic))}
    print("smooth inflow", kind, results)
    for rows in results.values():
        assert max(row["courant"] for row in rows) < .05
        assert all(rows[i+1]["linf"] < rows[i]["linf"] for i in range(2))
    second_errors = np.array([row["linf"] for row in results["quadratic"]])
    ratios = second_errors[:-1]/second_errors[1:]
    print("quadratic observed spatial ratios", kind, ratios.tolist())
    assert np.all(ratios > 3.) and np.all(ratios < 5.)
    assert results["quadratic"][-1]["linf"] < .2*results["first"][-1]["linf"]
    assert results["quadratic"][-1]["wall_jet_error"] < results["first"][-1]["wall_jet_error"]


@pytest.mark.parametrize("kind", ["constant", "eta_flow"])
def test_finite_homogeneous_amplification_with_fixed_inflow(kind, budget):
    observations = []
    for cells in (16, 32, 64):
        eta = np.linspace(0., 1., cells+1)
        # Smooth, oscillatory, and near-wall localized initial data.  The last
        # two deliberately stress a nonnormal closure; no spectral argument
        # is substituted for the measured finite propagation.
        initial = np.stack((np.sin(np.pi*eta), np.sin(.8*np.pi*np.arange(cells+1)),
                            eta*(1-eta)*np.exp(-80*(1-eta)**2)))
        initial[:, [0, -1]] = 0.
        budget.check(2*initial.size)
        norms = np.sqrt(np.trapezoid(initial**2, eta, axis=-1))
        initial /= norms[:, None]
        final = .5; speed_bound = 1. if kind == "constant" else .8
        steps = int(np.ceil(final*speed_bound*cells/.3)); dt = final/steps
        max_gain = 1.; maximum_courant = 0.
        def impose(z, t):
            z[:, [0, -1]] = 0.
            return z
        def rhs(t, z):
            nonlocal maximum_courant
            gradient = quadratic(eta, z, budget=budget)
            speed = np.ones_like(eta) if kind == "constant" else eta*(.6+.4*t)
            budget.check(z.size)
            out = speed*gradient
            out[:, [0, -1]] = 0.
            maximum_courant = max(maximum_courant, dt*cells*float(np.max(speed)))
            return out
        state = initial.copy()
        for step in range(steps):
            state, stages = _ssprk3(state, step*dt, dt, rhs, impose, budget)
            for stage in stages:
                budget.check(stage.size)
                assert np.all(np.isfinite(stage)) and np.all(stage[:, [0, -1]] == 0.)
                max_gain = max(max_gain, float(np.max(np.sqrt(np.trapezoid(stage**2, eta, axis=-1)))))
        observations.append({"cells": cells, "maximum_l2_gain": max_gain,
                             "maximum_courant": maximum_courant})
    print("finite homogeneous amplification", kind, observations)
    assert max(row["maximum_courant"] for row in observations) <= .30000000000001
    assert max(row["maximum_l2_gain"] for row in observations) < 1.5
    # This bounded three-grid experiment is not a grid-uniform norm estimate.


def test_closure_matrix_signs_are_only_a_local_spectral_check(budget):
    for a, b in ((1., 1.), (.8, .9)):
        h = 1/32
        block = np.array([[-3*a/(2*h), 2*a/h], [-b/(2*h), 0.]])
        budget.check(4)
        assert np.trace(block) < 0.
        assert np.linalg.det(block) == pytest.approx(a*b/h**2)
        assert np.all(np.real(np.linalg.eigvals(block)) < 0.)
        # The symmetric part is not negative definite: eigenvalues alone
        # cannot supply a uniform energy/stability claim for this closure.
        assert np.max(np.linalg.eigvalsh((block+block.T)/2)) > 0.
