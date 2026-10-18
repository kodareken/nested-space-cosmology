"""Fixed-environment Hamiltonian check independent of the projected rate formula."""
import numpy as np
import pytest

from recursive_horizons.nsc_spherical_galerkin_coupling import (
    blank_state, build_grid, compose_fine_hamiltonian, energy, manufactured_columns,
)
from recursive_horizons.nsc_conformal_resolution import (
    from_log_state, log_state, log_rates, project_geometry, project_rates,
)
from recursive_horizons.nsc_conformal_vacuum_generator import lambda_generator
from measure_nsc_conformal_resolution import saved_rate_replay_errors, GEOMETRY_FIELDS


@pytest.fixture(scope="module")
def example():
    grid = build_grid(16, gauge="conformal")
    state = blank_state(grid, *manufactured_columns(grid.nf))
    x = 2*np.pi*grid.xi_g/grid.length
    state.r = 2.*np.exp(.02*np.cos(x)+.012*np.cos(3*x))
    state.Q = np.exp(.01*np.sin(x)+.008*np.sin(4*x))
    state.chi = .01*np.cos(3*x)
    state.p_r = .004*np.sin(3*x)
    state.p_Q = .003*np.cos(2*x)
    state.p_chi = .002*np.cos(4*x)
    basis = np.column_stack((np.ones(grid.ng),np.cos(x),np.sin(x)))
    U, _ = np.linalg.qr(basis)
    P = U@U.T
    return grid,state,P


def test_correlation_gradient_at_fixed_environment(example):
    grid,state,P = example
    coarse = project_geometry(state,P)
    full_variables, coarse_variables = log_state(state),log_state(coarse)
    environment = {key:full_variables[key]-coarse_variables[key] for key in full_variables}
    full_rate, _ = compose_fine_hamiltonian(grid,state)
    coarse_rate, _ = compose_fine_hamiltonian(grid,coarse)
    xf,xc = log_rates(state,full_rate),log_rates(coarse,coarse_rate)
    delta = project_rates({key:xf[key]-xc[key] for key in xf},P)
    rng = np.random.default_rng(46)
    direction = {key:P@rng.normal(size=grid.ng) for key in full_variables}
    exact = grid.dx_g*sum(np.dot(delta[q],direction[p])-np.dot(delta[p],direction[q])
                          for q,p in (("a","pi_a"),("b","pi_b"),("chi","p_chi")))

    def correlation(step):
        retained = {key:coarse_variables[key]+step*direction[key] for key in direction}
        restored = {key:retained[key]+environment[key] for key in direction}
        return (energy(grid,from_log_state(restored,state))-
                energy(grid,from_log_state(retained,state)))

    step = 2e-6
    numerical = (correlation(step)-correlation(-step))/(2*step)
    assert numerical == pytest.approx(exact,rel=2e-5,abs=2e-6)


def test_energy_equality_does_not_imply_vacuum_force(example):
    grid,state,P = example
    state = state.copy()
    # Constant pi_a with unresolved radius gives a retained kinetic
    # contribution through <1/r^2>. A volume potential has no a rate.
    state.p_r = .1/state.r
    state.p_chi = np.zeros(grid.ng)
    coarse = project_geometry(state,P)
    full_rate, _ = compose_fine_hamiltonian(grid,state)
    coarse_rate, _ = compose_fine_hamiltonian(grid,coarse)
    xf,xc = log_rates(state,full_rate),log_rates(coarse,coarse_rate)
    actual = project_rates({key:xf[key]-xc[key] for key in xf},P)
    unit = lambda_generator(grid,coarse)
    coefficient = (energy(grid,state)-energy(grid,coarse))/unit["energy"]
    assert coefficient*unit["energy"] == pytest.approx(energy(grid,state)-energy(grid,coarse))
    candidate = project_rates(unit["canonical_rates"],P)
    assert np.max(abs(actual["a"]))>1e-8
    assert np.array_equal(candidate["a"],np.zeros(grid.ng))


def test_archive_rate_interface_uses_prolonged_projected_values(example):
    grid,state,P = example
    rate, _ = compose_fine_hamiltonian(grid,state)
    saved = {name:grid.A_g@getattr(rate,name) for name in GEOMETRY_FIELDS}
    assert max(saved_rate_replay_errors(grid,rate,saved).values()) == 0
    saved["r"] = np.zeros(grid.ng)
    with pytest.raises(ValueError,match="quadrature grid"):
        saved_rate_replay_errors(grid,rate,saved)
