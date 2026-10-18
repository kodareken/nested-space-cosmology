"""Independent small-grid controls of the full conformal Hamiltonian derivative."""
import time

import numpy as np
from numpy.testing import assert_allclose
import pytest

from recursive_horizons.nsc_conformal_tangent_evolution import (
    FIELDS, GEOMETRY, _constraint_mixed, conformal_jvp, conformal_tangent_constraints,
    variational_rk4_step,
)
from recursive_horizons.nsc_spherical_coupling import CauchyState
from recursive_horizons.nsc_spherical_galerkin_coupling import (
    _constraint_directional, blank_state, build_grid, compose_fine_hamiltonian,
    conformal_constraint_transport, energy, manufactured_columns, prolong_geometry,
    prolong_state, pull_geometry, rk4_step,
)

ACTIVE_CPU = [0.]


@pytest.fixture(autouse=True)
def control_budget():
    started = time.process_time()
    yield
    ACTIVE_CPU[0] += time.process_time()-started
    assert ACTIVE_CPU[0] < 30., "causal tangent manufactured controls exhausted 30 CPU seconds"


@pytest.fixture(scope="module")
def example():
    grid = build_grid(10, quadrature=40, gauge="conformal")
    state = blank_state(grid, *manufactured_columns(grid.nf))
    x = 2*np.pi*grid.xi_g/grid.length
    state.Q = .8 + .09*np.sin(x) + .01*np.cos(3*x)
    state.r = 2.1 + .13*np.cos(x) + .01*np.sin(4*x)
    state.chi = .04*np.sin(2*x)
    state.p_Q = .02*np.cos(2*x)
    state.p_r = .03*np.sin(x)
    state.p_chi = .01*np.cos(3*x)
    return grid, state


def direction(grid, seed):
    rng = np.random.default_rng(seed)
    values = {name: .1*rng.normal(size=grid.ng) for name in GEOMETRY}
    values.update({name: .1*(rng.normal(size=(grid.nf, 6)) +
                            1j*rng.normal(size=(grid.nf, 6))) for name in ("phi0", "phi1")})
    return CauchyState(**values)


def offset(state, tangent, step):
    return CauchyState(**{name: getattr(state, name)+step*getattr(tangent, name) for name in FIELDS})


def symplectic(grid, first, second):
    value = grid.dx_g*sum(np.sum(
        getattr(first, q)*getattr(second, p)-getattr(first, p)*getattr(second, q)
    ) for q, p in (("Q", "p_Q"), ("r", "p_r"), ("chi", "p_chi")))
    for name in ("phi0", "phi1"):
        value += 2*grid.fine.multiplicity*np.sum(grid.fine.occupations[None, :]*np.imag(
            getattr(first, name).conj()*getattr(second, name)))
    return float(value)


def test_all_eight_rates_match_an_independent_directional_difference(example):
    grid, state = example
    tangent = direction(grid, 31)
    step = 3e-6
    plus, _ = compose_fine_hamiltonian(grid, offset(state, tangent, step))
    minus, _ = compose_fine_hamiltonian(grid, offset(state, tangent, -step))
    actual = conformal_jvp(grid, state, tangent)
    for name in FIELDS:
        difference = (getattr(plus, name)-getattr(minus, name))/(2*step)
        assert_allclose(getattr(actual, name), difference, rtol=2e-6, atol=2e-7,
                        err_msg="full real directional derivative: "+name)


def test_full_geometry_and_weighted_column_symplectic_identity(example):
    grid, state = example
    first, second = direction(grid, 12), direction(grid, 18)
    df = conformal_jvp(grid, state, first)
    ds = conformal_jvp(grid, state, second)
    left, right = symplectic(grid, df, second), symplectic(grid, first, ds)
    assert abs(left+right) < 2e-10*max(1., abs(left), abs(right))


def test_zero_direction_and_reused_baseline_bundle(example):
    grid, state = example
    tangent = direction(grid, 2)
    zero = CauchyState(**{name: np.zeros_like(getattr(tangent, name)) for name in FIELDS})
    _, bundle = compose_fine_hamiltonian(grid, state)
    actual = conformal_jvp(grid, state, zero, bundle=bundle)
    for name in FIELDS:
        assert_allclose(getattr(actual, name), 0., atol=1e-15)
    direct = conformal_jvp(grid, state, tangent)
    reused = conformal_jvp(grid, state, tangent, bundle=bundle)
    for name in FIELDS:
        assert_allclose(getattr(direct, name), getattr(reused, name), rtol=1e-13, atol=1e-13)


def test_zero_momentum_conformal_radius_has_wave_acceleration(example):
    grid, original = example
    state = original.copy()
    state.Q = np.full(grid.ng, .7)
    for name in ("chi", "p_Q", "p_r", "p_chi"):
        setattr(state, name, np.zeros(grid.ng))
    rate, bundle = compose_fine_hamiltonian(grid, state)
    acceleration = conformal_jvp(grid, state, rate, bundle=bundle)
    exact = pull_geometry(grid, grid.derivative @ (
        grid.derivative @ prolong_geometry(grid, state.r)))
    assert_allclose(acceleration.r, exact, rtol=1e-10, atol=1e-10)


def test_wrong_gauge_and_nonreal_or_nonfinite_directions_are_rejected(example):
    grid, state = example
    tangent = direction(grid, 4)
    with pytest.raises(ValueError, match="conformal"):
        conformal_jvp(build_grid(10), state, tangent)
    tangent.Q = tangent.Q.astype(complex)
    with pytest.raises(ValueError, match="real"):
        conformal_jvp(grid, state, tangent)
    tangent.Q = np.full(grid.ng, np.nan)
    with pytest.raises(ValueError, match="finite"):
        conformal_jvp(grid, state, tangent)


def test_mixed_constraints_and_forcing_match_independent_directional_differences(example):
    grid, state = example
    first, second = direction(grid, 51), direction(grid, 23)
    rate, bundle = compose_fine_hamiltonian(grid, state)
    fine_first, fine_second = prolong_state(grid, first), prolong_state(grid, second)
    fields = lambda value: tuple(getattr(value, name) for name in FIELDS)
    actual = _constraint_mixed(bundle["fine_system"], bundle["fine_state"], fine_first, fine_second)
    step = 3e-6
    projected = []
    forcing = []
    current = []
    for sign in (1, -1):
        trial = offset(state, first, sign*step)
        trial_rate, trial_bundle = compose_fine_hamiltonian(grid, trial)
        _, _, h, d = _constraint_directional(
            trial_bundle["fine_system"], trial_bundle["fine_state"], fields(fine_second),
            trial_bundle["source"],
        )
        projected.append((h, d))
        forcing.append(conformal_constraint_transport(
            grid, trial, trial_rate, trial_bundle, include_vectors=True))
        current.append(trial_bundle["source"]["force_beta"]/grid.dx_q)
    for index in range(2):
        difference = (projected[0][index]-projected[1][index])/(2*step)
        assert_allclose(actual[index], difference, rtol=3e-6, atol=2e-6)
    report = conformal_tangent_constraints(grid, state, first, rate=rate,
                                          bundle=bundle, include_vectors=True)
    for name in ("forcing_h", "forcing_d"):
        difference = (forcing[0][name]-forcing[1][name])/(2*step)
        assert_allclose(report[name], difference, rtol=3e-6, atol=2e-6,
                        err_msg="variation of the actual constraint forcing: "+name)
    assert_allclose(report["shift_source"], (current[0]-current[1])/(2*step),
                    rtol=3e-6, atol=2e-6)
    assert report["constraint_tangent_norm"] > 1e-3
    assert report["continuum_constraint_certified"] is False
    assert report["sbp_identity_error"] < 1e-9*max(
        1., abs(report["constraint_tangent_energy_rate"]))


def test_variational_rk4_is_the_derivative_of_the_existing_discrete_step(example):
    grid, state = example
    tangent = direction(grid, 7)
    dt, epsilon = .001, 3e-6
    evolved, response = variational_rk4_step(grid, state, tangent, dt)
    reference = rk4_step(grid, state, dt)
    plus = rk4_step(grid, offset(state, tangent, epsilon), dt)
    minus = rk4_step(grid, offset(state, tangent, -epsilon), dt)
    for name in FIELDS:
        assert_allclose(getattr(evolved, name), getattr(reference, name), rtol=1e-13, atol=1e-13)
        numerical = (getattr(plus, name)-getattr(minus, name))/(2*epsilon)
        assert_allclose(getattr(response, name), numerical, rtol=2e-6, atol=2e-8)
    rate, _ = compose_fine_hamiltonian(grid, evolved)
    variation = symplectic(grid, rate, response)
    independent = (energy(grid, plus)-energy(grid, minus))/(2*epsilon)
    assert_allclose(variation, independent, rtol=2e-6, atol=2e-7)
