"""Equation and continuation controls for one bounded, fixed-charge saddle.

The module-scoped pair of BVP solutions is generated here, not loaded from an
old cache. These tests do not establish the full physical fluctuation index.
"""
import math

import numpy as np
import pytest
from scipy.integrate import simpson
from threadpoolctl import threadpool_limits

from recursive_horizons import nsc_atlas_shape_transition as transition


@pytest.fixture(scope="module")
def bounces():
    with threadpool_limits(limits=1):
        coarse = transition.solve_bounce()
        fine = transition.solve_bounce(
            tolerance=5e-8, cutoff=5e-5, mesh=321, previous=coarse)
    return coarse, fine


def _direct_potential(fields):
    a, b = fields
    return math.exp(-3*(a+3*b)) * (
        7 - 15/8*(math.exp(-2*a)+3*math.exp(-2*b))
        + 25/64*(math.exp(-6*(a+b))+math.exp(-12*b)))


def test_restricted_kinetic_metric_and_potential_derivatives():
    np.testing.assert_array_equal(transition.METRIC, [[7.5, 13.5], [13.5, 49.5]])
    velocity = np.array([.31, -.12])
    full_velocity = np.array([velocity[0], *([velocity[1]]*3)])
    direct_norm = 3*(full_velocity@full_velocity) + 4.5*sum(full_velocity)**2
    assert velocity@transition.METRIC@velocity == pytest.approx(direct_norm, abs=2e-16)
    point = np.array([-.12, .047])
    assert transition.potential(point) == pytest.approx(_direct_potential(point), abs=2e-15)
    step = 2e-5
    basis = np.eye(2)
    direct_gradient = np.array([
        (_direct_potential(point+step*e)-_direct_potential(point-step*e))/(2*step)
        for e in basis])
    direct_hessian = np.empty((2, 2))
    zero = _direct_potential(point)
    for i, e in enumerate(basis):
        direct_hessian[i, i] = (
            _direct_potential(point+step*e)-2*zero+_direct_potential(point-step*e))/step**2
        for j in range(i):
            f = basis[j]
            direct_hessian[i, j] = direct_hessian[j, i] = (
                _direct_potential(point+step*(e+f))-_direct_potential(point+step*(e-f))
                -_direct_potential(point+step*(-e+f))+_direct_potential(point-step*(e+f))
            )/(4*step**2)
    np.testing.assert_allclose(transition.gradient(point), direct_gradient, rtol=3e-7, atol=2e-8)
    np.testing.assert_allclose(transition.hessian(point), direct_hessian, rtol=3e-6, atol=2e-5)


def test_exact_constant_field_four_sphere_has_the_correct_euclidean_signs():
    false, _, _ = transition.reference_points()
    fields = np.asarray(false["log_radii"])[[0, 1]]
    H = math.sqrt(false["V"]/3)
    length = math.pi/H
    x = np.array([.08, .23, .61, .94])
    radius = np.sin(math.pi*x)/H
    slope = np.cos(math.pi*x)
    state = np.vstack([np.repeat(fields[:, None], len(x), axis=1),
                       np.zeros((2, len(x))), radius, slope])
    derivative = transition.euclidean_rhs(x, state, length)
    expected = np.vstack([np.zeros((4, len(x))), length*slope, -length*H**2*radius])
    np.testing.assert_allclose(derivative, expected, rtol=2e-13, atol=2e-11)
    np.testing.assert_allclose(slope**2-1+radius**2*false["V"]/3, 0., atol=3e-16)


def test_nonconstant_profile_rejects_the_hawking_moss_first_seed(bounces):
    coarse, _ = bounces
    assert len(coarse.attempts) == 2
    assert coarse.attempts[0]["seed_amplitude"] == .4
    assert not coarse.attempts[0]["nonconstant_profile"]
    assert coarse.attempts[1]["seed_amplitude"] == .8
    assert coarse.attempts[1]["nonconstant_profile"]
    assert np.max(np.ptp(coarse.solution.y[:2], axis=1)) > .3
    assert np.min(coarse.solution.y[4]) > 0
    assert coarse.solution.y[5, 0] > .999999 and coarse.solution.y[5, -1] < -.999999


def test_refined_exit_and_action_converge_without_a_new_exit_state(bounces):
    coarse, fine = bounces
    a = transition.bounce_diagnostics(coarse)
    b = transition.bounce_diagnostics(fine)
    np.testing.assert_allclose(coarse.pole_fields(), fine.pole_fields(), rtol=0., atol=2e-6)
    assert abs(coarse.length-fine.length) < 2e-5
    assert abs(a["action_trace"]-b["action_trace"]) < 2e-4
    assert b["max_Euclidean_constraint"] < a["max_Euclidean_constraint"]
    assert b["max_Euclidean_constraint"] < 5e-8
    assert not b["full_physical_negative_mode_count_proved"]


def test_constraint_and_full_action_identity_keep_the_pole_amplification_visible(bounces):
    _, bounce = bounces
    diagnostic = transition.bounce_diagnostics(bounce)
    assert abs(diagnostic["action_identity_difference"]) < 2e-5
    x = np.linspace(bounce.cutoff, 1-bounce.cutoff, 4001)
    y = bounce.solution.sol(x)
    dy = bounce.solution.sol(x, 1)
    radius, slope = y[4:]
    velocity = y[2:4].T
    kinetic = np.einsum("ni,ij,nj->n", velocity, transition.METRIC, velocity)
    potential = transition.potential(y[:2].T)
    actual_radius_second = dy[5]/bounce.length
    constraint = slope**2-1-radius**2*(kinetic/2-potential)/3
    hamiltonian = -4*actual_radius_second/radius + 2*(1-slope**2)/radius**2-kinetic-2*potential
    acceleration_residual = actual_radius_second + radius*(kinetic+potential)/3
    np.testing.assert_allclose(
        hamiltonian, -4*acceleration_residual/radius-2*constraint/radius**2,
        rtol=0., atol=2e-9)
    weighted_rms = math.sqrt(simpson(radius**2*hamiltonian**2, x=x)/simpson(radius**2, x=x))
    raw_maximum = float(np.max(abs(hamiltonian)))
    assert np.isfinite(raw_maximum)
    assert raw_maximum == pytest.approx(diagnostic["max_common_Cauchy_Hamiltonian"], abs=2e-14)
    # The actual maximum is retained above. A small first-integral residual
    # need not give a small pointwise Hamiltonian near the shrinking pole.
    assert weighted_rms < 2e-5, f"weighted RMS={weighted_rms}, unweighted max={raw_maximum}"
    assert max(abs(constraint)) < 5e-8


def test_open_pole_is_regular_and_reaches_its_first_turn_without_material_inflow(bounces):
    _, bounce = bounces
    with threadpool_limits(limits=1):
        solution, diagnostic = transition.continue_pole(bounce)
    assert diagnostic["endpoint_reason"] == "first_scale_factor_turn"
    assert solution.t_events[0].size == 1
    assert abs(solution.y[5, -1]) < 2e-12
    assert diagnostic["scale_at_endpoint"] > 0
    assert diagnostic["incoming_material_flux"] == 0
    assert diagnostic["maximum_Friedmann_constraint"] < 2e-8
    assert abs(diagnostic["density_work_balance_residual"]) < 2e-8
    # Compute near-pole invariants from the actual scale-factor history,
    # rather than declaring a=0 to be a curvature singularity.
    t = .002
    step = 1e-4
    middle = solution.sol(t)
    a, slope = middle[4:]
    acceleration = (solution.sol(t+step)[4]-2*a+solution.sol(t-step)[4])/step**2
    radial_curvature = acceleration/a
    tangential_curvature = (slope**2-1)/a**2
    ricci = 6*(radial_curvature+tangential_curvature)
    kretschmann = 12*(radial_curvature**2+tangential_curvature**2)
    assert ricci == pytest.approx(diagnostic["initial_Ricci_scalar_limit"], rel=3e-4, abs=2e-6)
    assert kretschmann == pytest.approx(diagnostic["initial_Kretschmann_limit"], rel=6e-4, abs=1e-6)


def test_shared_cauchy_equator_is_a_spatial_maximum_with_covariant_screen_area(bounces):
    _, bounce = bounces
    geometry = transition.connecting_geometry(bounce)
    assert geometry["T0_Kij"] == geometry["T0_normal_field_momenta"] == 0
    assert geometry["normal_material_flux"] == 0
    assert geometry["T0_common_Cauchy_volume"] > 0
    assert geometry["T0_common_Cauchy_normal_energy"] > 0
    x = geometry["equator_x"]
    y = bounce.solution.sol(x)
    assert abs(y[5]) < 2e-12
    assert geometry["equator_spatial_second_derivative"] < 0
    assert geometry["equator_geodesic_radius_second_proper_derivative"] > 0
    expansions = geometry["equator_null_expansions"]
    assert all(value < 0 for value in expansions["-1.0"])
    assert expansions["0.0"] == [0., 0.]
    assert all(value > 0 for value in expansions["1.0"])
    assert geometry["dimensional_screen_area_ratio"] == pytest.approx(1., abs=5e-16)
