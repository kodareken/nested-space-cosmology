"""Focused ownership/sign checks of the dynamic LLL geometric allocation."""
import numpy as np
import pytest

from recursive_horizons.nsc_horizon_source import conformal_stress
from recursive_horizons.nsc_lorentzian import geometry
from recursive_horizons.nsc_lll_geometric_history import LLLGeometricHistoryAllocation
from recursive_horizons.nsc_spherical_local_history import SphericalHistoryGrid, OwnedCompactKSHistory


def stationary_control(history):
    x = history.grid.radius_coordinate
    beta, bp, A = geometry(x)
    Ap = -2*beta*bp
    App = -6*(np.pi/2-np.arctan(x))+6*x/(1+x*x)
    # Existing allocation with zero state constants, not the total LLL tensor.
    uu, uv, vv = conformal_stress(A, Ap, App, 0., 0., central_charge=4)
    expected = np.stack(((uu-2*uv+vv)/(-A), (-uu+vv)/(-A),
                         (uu+2*uv+vv)/(-A), np.zeros_like(A)), axis=-1)
    got = LLLGeometricHistoryAllocation(4).evaluate(history, np.zeros(4))
    return float(np.max(abs(got['stress_2D_allocation']-expected[None, ...])))


def test_inherits_stationary_allocation_and_factor_once():
    h = OwnedCompactKSHistory(SphericalHistoryGrid.gaussian(12, 12))
    assert stationary_control(h) < 3e-11
    a = [.002, .001, .003, .001]
    got = LLLGeometricHistoryAllocation(4).evaluate(h, a)
    r = got['geometry']['r']
    assert np.max(abs(got['stress_4D_allocation']*4*np.pi*r[..., None]**2
                      -got['stress_2D_allocation'])) < 3e-11
    assert np.max(abs(got['stress_2D_allocation'][..., 3])) == 0
    assert got['action_gradient_from_stress'][3] == 0


def test_compact_history_action_gradient_matches_tensor_pairing():
    h = OwnedCompactKSHistory(SphericalHistoryGrid.gaussian(64, 64))
    owner = LLLGeometricHistoryAllocation(4); a = [.002, .001, .003, .001]
    got = owner.evaluate(h, a)
    assert np.max(abs(owner.action_gradient(h, a)-got['action_gradient_from_stress'])) < 3e-11
    # Local correction has a genuine momentum channel on spatially varying g.
    assert np.max(abs(got['stress_2D_allocation'][..., 1])) > 1e-3
    with pytest.raises(ValueError): LLLGeometricHistoryAllocation(0)
    with pytest.raises(ValueError): LLLGeometricHistoryAllocation(True)
    with pytest.raises(ValueError): LLLGeometricHistoryAllocation(4.5)
