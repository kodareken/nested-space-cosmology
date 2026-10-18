"""Manufactured controls for a fixed reference resolution. No archived state is read."""
import numpy as np
import pytest

from recursive_horizons.nsc_conformal_resolution import (
    from_log_state, log_rates, log_state, project_geometry, project_rates,
    reference_projector,
)
from recursive_horizons.nsc_spherical_coupling import CauchyRate, CauchyState
from recursive_horizons.nsc_spherical_galerkin_coupling import build_grid


def _grid():
    return build_grid(10)


def _state(grid, radius=1.7, ratio=0.6):
    count = grid.ng
    coordinate = np.linspace(-0.2, 0.3, count)
    phi0 = np.zeros((grid.nf, 6), dtype=complex)
    phi1 = np.zeros((grid.nf, 6), dtype=complex)
    phi0[0, 0] = 1.0 + 0.2j
    return CauchyState(
        Q=np.full(count, ratio) + 0.05 * np.sin(coordinate),
        r=np.full(count, radius) + 0.03 * np.cos(coordinate),
        chi=0.1 * coordinate,
        p_Q=0.02 * coordinate,
        p_r=-0.04 * coordinate,
        p_chi=0.01 * np.cos(coordinate),
        phi0=phi0,
        phi1=phi1,
    )


def _rate(state):
    count = state.r.shape[0]
    coordinate = np.linspace(0.01, -0.02, count)
    return CauchyRate(
        Q=0.03 * coordinate, r=0.02 * np.cos(coordinate), chi=-0.01 * coordinate,
        p_Q=0.04 * np.sin(coordinate), p_r=0.05 * coordinate, p_chi=np.full(count, 0.006),
        phi0=np.zeros_like(state.phi0), phi1=np.zeros_like(state.phi1),
        force_L=np.zeros(count), force_Q=np.zeros(count), force_beta=np.zeros(count),
        fieldwork_power=0.0,
    )


def _mock_grid(prolongation, derivative):
    class MockGrid:
        pass

    grid = MockGrid()
    grid.A_g = np.asarray(prolongation, dtype=float)
    grid.derivative = np.asarray(derivative, dtype=float)
    grid.nq, grid.ng = grid.A_g.shape
    grid.dx_q = 0.5
    return grid


def test_constant_metric_has_analytic_spectrum_and_scales_lengths():
    grid = _grid()
    radius, ratio = 2.0, 0.5
    reference = np.full(grid.ng, radius)
    lapse_ratio = np.full(grid.ng, ratio)
    ruler = 1.0
    projector = reference_projector(grid, reference, lapse_ratio, ruler)
    wave = 2 * np.pi * grid.modes_g / grid.length
    expected = np.sort((wave / (radius * ratio)) ** 2)
    assert projector.eigenvalues == pytest.approx(expected, abs=1e-8, rel=1e-8)
    assert projector.retained == int(np.count_nonzero(expected * ruler**2 <= 1))
    assert projector.diagnostics["eigen_residual"] < 1e-8
    scale = 1e6
    scaled = reference_projector(grid, scale * reference, lapse_ratio, scale * ruler)
    assert scaled.eigenvalues*scale**2 == pytest.approx(projector.eigenvalues, abs=1e-8, rel=1e-8)
    assert scaled.retained == projector.retained
    assert scaled.matrix == pytest.approx(projector.matrix, abs=1e-8)
    wide = reference_projector(grid, reference, lapse_ratio, 1e6)
    assert wide.retained >= 1
    assert wide.matrix @ np.ones(grid.ng) == pytest.approx(np.ones(grid.ng), abs=1e-8)
    assert wide.diagnostics["smallest_excluded_eigenvalue"] is not None


def test_projector_is_idempotent_symmetric_and_canonically_orthogonal():
    grid = _grid()
    projector = reference_projector(grid, np.full(grid.ng, 1.4), np.full(grid.ng, 0.8), 2.5)
    matrix = projector.matrix
    assert matrix @ matrix == pytest.approx(matrix, abs=1e-8)
    assert matrix.T == pytest.approx(matrix, abs=1e-8)
    generator = np.random.default_rng(3)
    left, right = generator.normal(size=grid.ng), generator.normal(size=grid.ng)
    pairing = grid.dx_g * np.dot(matrix @ left, right - matrix @ right)
    assert pairing == pytest.approx(0.0, abs=1e-8)
    state = _state(grid)
    once = project_geometry(state, matrix)
    twice = project_geometry(once, matrix)
    assert log_state(twice)["a"] == pytest.approx(log_state(once)["a"], abs=1e-10)
    assert log_state(twice)["pi_b"] == pytest.approx(log_state(once)["pi_b"], abs=1e-10)


def test_identity_resolution_round_trips_geometry_and_keeps_fermions():
    grid = _grid()
    state = _state(grid)
    radius = state.r.copy()
    phi0 = state.phi0.copy()
    projector = reference_projector(grid, np.full(grid.ng, 1.0), np.full(grid.ng, 1.0), 1e-3)
    assert projector.retained == grid.ng
    assert projector.matrix == pytest.approx(np.eye(grid.ng), abs=1e-8)
    assert projector.diagnostics["smallest_excluded_eigenvalue"] is None
    projected = project_geometry(state, projector.matrix)
    recovered = from_log_state(log_state(state), state)
    assert projected.r == pytest.approx(state.r, abs=1e-10)
    assert projected.Q == pytest.approx(state.Q, abs=1e-10)
    assert projected.p_r == pytest.approx(state.p_r, abs=1e-10)
    assert projected.p_Q == pytest.approx(state.p_Q, abs=1e-10)
    assert recovered.chi == pytest.approx(state.chi, abs=1e-12)
    assert np.array_equal(projected.phi0, phi0)
    assert np.array_equal(projected.phi1, state.phi1)
    projected.phi0[0, 0] = 9.0
    assert np.array_equal(state.phi0, phi0)
    assert state.r == pytest.approx(radius, abs=0.0)
    rates = project_rates(log_rates(state, _rate(state)), projector.matrix)
    assert rates["a"] == pytest.approx(log_rates(state, _rate(state))["a"], abs=1e-12)


def test_log_rate_matches_a_central_finite_difference():
    grid = _grid()
    state = _state(grid)
    rate = _rate(state)
    step = 1e-6

    def displaced(sign):
        return CauchyState(
            Q=state.Q + sign * step * rate.Q,
            r=state.r + sign * step * rate.r,
            chi=state.chi + sign * step * rate.chi,
            p_Q=state.p_Q + sign * step * rate.p_Q,
            p_r=state.p_r + sign * step * rate.p_r,
            p_chi=state.p_chi + sign * step * rate.p_chi,
            phi0=state.phi0.copy(), phi1=state.phi1.copy(),
        )

    derivative = log_rates(state, rate)
    for name in derivative:
        difference = (log_state(displaced(1))[name] - log_state(displaced(-1))[name]) / (2 * step)
        assert derivative[name] == pytest.approx(difference, abs=1e-8, rel=1e-7)


def test_bad_shapes_and_prolonged_sign_fail_without_nodal_failure():
    with pytest.raises(ValueError, match="derivative"):
        reference_projector(_mock_grid(np.eye(4, 3), np.eye(3)), np.ones(3), np.ones(3), 1.0)
    prolongation = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [-2.0, 0.4, 0.4]])
    assert np.min(prolongation @ np.ones(3)) < 0
    with pytest.raises(ValueError, match="prolonged"):
        reference_projector(_mock_grid(prolongation, np.zeros((4, 4))), np.ones(3), np.ones(3), 1.0)
    grid = _grid()
    state = _state(grid)
    values = log_state(state)
    values["a"] = np.full(grid.ng, np.inf)
    with pytest.raises(ValueError):
        from_log_state(values, state)
    values["a"] = np.full(grid.ng,-1000.)
    with pytest.raises(ValueError):
        from_log_state(values,state)
