"""Analytic Fourier control of the preparation's source/radius normalization.

The manufactured r=1 background is an operator control, not the sourced
initial radius root. All first variations are single resolved harmonics.
"""
import numpy as np
from numpy.testing import assert_allclose

from recursive_horizons.nsc_spherical_galerkin_coupling import (
    _constraint_directional, _projected_radius_jacobian, blank_state,
    build_grid, prolong_columns, prolong_state, pull_geometry,
)


def test_paired_plane_wave_phase_has_the_analytic_radius_response():
    grid = build_grid(16, gauge="conformal")
    columns0 = np.empty((grid.nf, 6), dtype=complex)
    columns1 = np.empty_like(columns0)
    for pair in range(3):
        carrier = np.exp(1j * (2*pair+1)*np.pi*grid.xi_f/grid.length)
        columns0[:, 2*pair] = carrier / np.sqrt(2*grid.nf)
        columns1[:, 2*pair] = 1j*columns0[:, 2*pair]
        columns0[:, 2*pair+1] = columns0[:, 2*pair].conj()
        columns1[:, 2*pair+1] = columns1[:, 2*pair].conj()
    state = blank_state(grid, columns0, columns1)
    fine = prolong_state(grid, state)
    q = state.Q[0]
    wave = 2*np.pi/grid.length
    generator = np.cos(wave*grid.xi_f)
    dphi0 = generator[:, None]*columns1
    dphi1 = -generator[:, None]*columns0
    zero = np.zeros(grid.nq)
    direction = (zero,)*6 + (
        prolong_columns(grid, dphi0), prolong_columns(grid, dphi1),
    )
    _, _, h_dot, current_dot = _constraint_directional(grid.fine, fine, direction)

    # Every normalized chiral column contributes f'/nq to delta_K.
    # The weights sum to 3; multiply by M once and divide by dx_q.
    factor = np.sum(grid.fine.occupations)*grid.fine.multiplicity/grid.length
    derivative = -wave*np.sin(wave*grid.xi_q)
    density_dot = factor*derivative/q
    assert_allclose(h_dot/q, density_dot, rtol=1e-11, atol=1e-11)
    assert_allclose(current_dot, 0., atol=1e-11)

    # At constant r=1, J=16*pi*A*(D^2/q-q), hence each mode has
    # eigenvalue -16*pi*A*(q+k^2/q). This is derived without J itself.
    radius_dot = (
        factor*(-wave*np.sin(wave*grid.xi_g))
        / (16*np.pi*grid.fine.A*(q*q+wave*wave))
    )
    jacobian = _projected_radius_jacobian(grid, state.r)
    assert_allclose(jacobian@radius_dot, -pull_geometry(grid, density_dot),
                    rtol=1e-10, atol=1e-10)
    combined = list(direction)
    combined[1] = grid.A_g@radius_dot
    _, _, total_h_dot, total_d_dot = _constraint_directional(
        grid.fine, fine, tuple(combined),
    )
    assert_allclose(total_h_dot/q, 0., atol=1e-10)
    assert_allclose(total_d_dot, 0., atol=1e-11)
