"""Holomorphic frequency extension of the published real PG/Jost owner."""
import numpy as np
import pytest

from recursive_horizons.nsc_complex_horizon_modes import (
    analytic_horizon_frame, complex_outgoing_ratio, complex_reflection,
    finite_pg_horizon_basis, outgoing_momentum, wirtinger_residual,
)
from recursive_horizons.nsc_massive_jost_modes import outgoing_ratio, solve_jost
from recursive_horizons.nsc_paired_horizon_preparation import PairedHorizonSeedMap
from recursive_horizons.nsc_pg_massive_modes import MassivePGModeResolution
from recursive_horizons.nsc_unruh_state import horizon_frame

P = PairedHorizonSeedMap(1.9006916054701435, .23832579963401956, 3.973074368754331, 1e-10, 1e-9)
MASS = np.pi / 2
ANGULAR = 0.
KAPPA = P.surface_gravity
BACKGROUND = P.background
COLLAR = 1e-10
DISTANCE = np.sqrt(2 * COLLAR / KAPPA)
TORTOISE = BACKGROUND.near_tortoise(COLLAR)
EFFECTIVE = MASS * P.horizon_radius
Z_OFF = 1.2 + 0.04j


@pytest.fixture(scope='module')
def above_gap_section():
    energy = 2.5
    jost = solve_jost(BACKGROUND, energy, MASS, ANGULAR)
    section = MassivePGModeResolution(P).section(energy, MASS, ANGULAR, [0.], [3.], jost=jost)
    return energy, jost, section


def test_analytic_frame_matches_real_horizon_frame():
    for timelike in (False, True):
        for frequency in (0.2, 0.7, 1.2, MASS - 0.05, 2.5, 3.2):
            owned = horizon_frame(frequency, KAPPA, EFFECTIVE, DISTANCE, TORTOISE, timelike=timelike)
            analytic = analytic_horizon_frame(
                frequency, KAPPA, EFFECTIVE, DISTANCE, TORTOISE, timelike=timelike)
            assert np.linalg.norm(owned - analytic) < 3e-14


def test_analytic_frame_cauchy_riemann():
    for timelike, point in ((False, Z_OFF), (True, 2.5 + 0.04j)):
        residual = wirtinger_residual(
            lambda z: analytic_horizon_frame(z, KAPPA, EFFECTIVE, DISTANCE, TORTOISE, timelike=timelike),
            point, 1e-6)
        assert residual['relative'] < 1e-8
        assert residual['cauchy_riemann'] < 1e-6


def test_second_column_is_conjugate_analytic_not_z_conjugation():
    analytic = analytic_horizon_frame(Z_OFF, KAPPA, EFFECTIVE, DISTANCE, TORTOISE, timelike=True)
    a, b = analytic[0, 0], analytic[1, 0]
    conjugated_z = np.array([[a, -b.conjugate()], [b, a.conjugate()]], complex)
    assert np.linalg.norm(analytic - conjugated_z) > 1
    at_conjugate = analytic_horizon_frame(
        Z_OFF.conjugate(), KAPPA, EFFECTIVE, DISTANCE, TORTOISE, timelike=True)
    assert np.linalg.norm(at_conjugate - analytic.conj()) > 1
    real = analytic_horizon_frame(1.2, KAPPA, EFFECTIVE, DISTANCE, TORTOISE, timelike=True)
    a, b = real[0, 0], real[1, 0]
    assert np.linalg.norm(real - np.array([[a, -b.conjugate()], [b, a.conjugate()]], complex)) < 3e-14


def test_complex_series_matches_real_outgoing_ratio():
    for energy, radius in ((0.7, 80.), (1.2, 200.), (2.5, 200.), (3.2, 120.)):
        owned = outgoing_ratio(energy, MASS, ANGULAR, radius)
        holomorphic = complex_outgoing_ratio(energy, MASS, ANGULAR, radius)
        assert abs(owned[0] - holomorphic[0]) < 3e-14
        assert abs(owned[1] - holomorphic[1]) < 3e-14
    below = outgoing_momentum(1.0 + 0.05j, MASS)
    assert below.real >= 0 and below.imag > 0
    above = outgoing_momentum(2.5, MASS)
    assert above.imag == 0 and above.real > 0


def test_complex_reflection_matches_real_jost(above_gap_section):
    energy, jost, _ = above_gap_section
    holomorphic = complex_reflection(energy, MASS, ANGULAR, BACKGROUND)
    assert abs(holomorphic.reflection - jost.reflection) < 3e-11
    assert holomorphic.residuals['outer_Riccati_residual'] < 3e-13
    for frequency in (0.7, 1.2):
        owned = solve_jost(BACKGROUND, frequency, MASS, ANGULAR)
        value = complex_reflection(frequency, MASS, ANGULAR, BACKGROUND)
        assert abs(value.reflection - owned.reflection) < 3e-11


def test_complex_reflection_cauchy_riemann():
    below = wirtinger_residual(
        lambda z: complex_reflection(z, MASS, ANGULAR, BACKGROUND).reflection,
        1.0 + 0.06j, 2e-3, order=4)
    above = wirtinger_residual(
        lambda z: complex_reflection(z, MASS, ANGULAR, BACKGROUND).reflection,
        2.2 + 0.05j, 1e-3)
    assert below['relative'] < 3e-4
    assert above['relative'] < 3e-6


def test_finite_basis_matches_real_pg_columns(above_gap_section):
    energy, jost, section = above_gap_section
    basis = finite_pg_horizon_basis(energy, MASS, ANGULAR, P)
    # Propagated fields use the declared numerical 3e-11 tolerance;
    # the algebraic Frobenius-frame identity above retains 3e-14.
    assert np.linalg.norm(basis.interior - section.interior_fundamental[0]) < 3e-11
    assert np.linalg.norm(basis.exterior[:, 1] - section.exterior_ingoing[0]) < 3e-9
    assert np.linalg.norm(basis.interior - basis.interior_conjugate) < 3e-14
    assert basis.residuals['pg_clock_is_geometric'] == 0
    assert basis.residuals['frequency']['interior_ode'] < 3e-9
    assert basis.residuals['frequency']['exterior_ode'] < 3e-9
    assert jost.mass == MASS


def test_finite_basis_conjugate_analytic_on_both_charts():
    basis = finite_pg_horizon_basis(Z_OFF, MASS, ANGULAR, P)
    assert basis.interior.shape == basis.exterior.shape == (2, 2)
    assert basis.conjugate_frequency == Z_OFF.conjugate()
    assert basis.residuals['naive_interior_conjugation'] > 1e-2
    assert basis.residuals['naive_exterior_conjugation'] > 1e-2
    assert basis.residuals['pg_clock_is_geometric'] < 3e-14
    for key in ('frequency', 'conjugate_frequency'):
        assert basis.residuals[key]['interior_ode'] < 3e-9
        assert basis.residuals[key]['exterior_ode'] < 3e-9
        assert basis.residuals[key]['interior_Dirac_frame'] < 3e-11
        assert basis.residuals[key]['exterior_Dirac_frame'] < 3e-11
    cr = wirtinger_residual(
        lambda z: finite_pg_horizon_basis(z, MASS, ANGULAR, P).interior, 1.2 + 0.05j, 2e-3)
    assert cr['relative'] < 3e-4


def test_domain_rejects_pole_state_and_lower_half():
    with pytest.raises(ValueError, match='kappa/2'):
        analytic_horizon_frame(0.2 + 1j * (KAPPA / 2), KAPPA, EFFECTIVE, DISTANCE, TORTOISE)
    with pytest.raises(ValueError, match='history duration'):
        finite_pg_horizon_basis(1.2 + 1j * KAPPA, MASS, ANGULAR, P)
    with pytest.raises(ValueError, match='Im z'):
        complex_reflection(1.2 - 0.01j, MASS, ANGULAR, BACKGROUND)
    with pytest.raises(TypeError):
        complex_reflection(2.5, MASS, ANGULAR, BACKGROUND, seed_covariance=np.eye(2))
    with pytest.raises(TypeError):
        finite_pg_horizon_basis(2.5, MASS, ANGULAR, P, seed_covariance=np.eye(2))
    with pytest.raises(TypeError):
        finite_pg_horizon_basis(2.5, MASS, ANGULAR, BACKGROUND)
