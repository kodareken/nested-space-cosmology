"""Independent volume and directional-variation controls, no archive evaluation."""
import numpy as np
import pytest

from recursive_horizons.nsc_spherical_coupling import CauchyState
from recursive_horizons.nsc_spherical_galerkin_coupling import build_grid
from recursive_horizons.nsc_conformal_vacuum_generator import (
    lambda_generator, ungauged_volume_jet,
)


@pytest.fixture(scope="module")
def grid():
    return build_grid(16, gauge="conformal")


def state_for(grid):
    x = grid.xi_g * 2*np.pi/grid.length
    zeros = np.zeros(grid.ng)
    phi = np.zeros((grid.nf, grid.fine.occupations.size), dtype=complex)
    return CauchyState(1.1+.03*np.cos(x), 2.+.04*np.sin(2*x), zeros.copy(),
                       zeros.copy(), zeros.copy(), zeros.copy(), phi, phi.copy())


def test_ungauged_volume_is_four_dimensional():
    r, L, Q, A = 2., .7, 1.3, .2
    jet = ungauged_volume_jet(r, L, Q, A)
    # The independent ADM volume is 4*pi*N*q*r^2, rho=2*A*Lambda.
    assert jet["density"] == pytest.approx(4*np.pi*(r*L)*(r*Q)*r*r*(2*A))
    for name, value in (("dr", r), ("dL", L), ("dQ", Q)):
        args = [r, L, Q]
        i = {"dr": 0, "dL": 1, "dQ": 2}[name]
        step = 1e-5*value
        args[i] += step
        plus = ungauged_volume_jet(*args, A)["density"]
        args[i] -= 2*step
        minus = ungauged_volume_jet(*args, A)["density"]
        assert (plus-minus)/(2*step) == pytest.approx(jet[name], rel=3e-9)


def test_unit_generator_only_changes_momenta(grid):
    result = lambda_generator(grid, state_for(grid))
    assert result["energy"] > 0
    for name in ("a", "b", "chi", "p_chi"):
        assert np.array_equal(result["canonical_rates"][name], np.zeros(grid.ng))
    assert np.max(result["canonical_rates"]["pi_a"]) < 0
    assert np.max(result["canonical_rates"]["pi_b"]) < 0
    assert result["matter_rate_is_zero"]


def test_canonical_gradient_matches_independent_volume_difference(grid):
    state = state_for(grid)
    directions = (np.cos(2*np.pi*grid.xi_g/grid.length),
                  np.sin(4*np.pi*grid.xi_g/grid.length))
    base = lambda_generator(grid, state)
    step = 1e-5
    for field, momentum, direction in zip(("r", "Q"), ("pi_a", "pi_b"), directions):
        plus, minus = state.copy(), state.copy()
        setattr(plus, field, getattr(state, field)*np.exp(step*direction))
        setattr(minus, field, getattr(state, field)*np.exp(-step*direction))
        numerical = (lambda_generator(grid, plus)["energy"] -
                     lambda_generator(grid, minus)["energy"])/(2*step)
        canonical = -grid.dx_g*np.dot(base["canonical_rates"][momentum], direction)
        assert numerical == pytest.approx(canonical, rel=2e-7, abs=2e-7)


@pytest.mark.parametrize("r,L,Q,A", [(0,1,1,1), (1,-1,1,1), (1,1,np.nan,1), (1,1,1,0)])
def test_invalid_volume_is_rejected(r,L,Q,A):
    with pytest.raises(ValueError):
        ungauged_volume_jet(r,L,Q,A)
