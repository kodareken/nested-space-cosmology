"""Analytic and mutation controls for finite-band Jost phase transport."""
import math

import numpy as np
import pytest
from flint import arb, ctx
from scipy.integrate._ivp.common import OdeSolution

from recursive_horizons.nsc_ks_ball_trajectory import restored_upper
from recursive_horizons.nsc_ks_trajectory import TrajectorySegment
from recursive_horizons.nsc_massive_jost_modes import solve_jost
from recursive_horizons.nsc_massive_jost_transport_bound import (
    MISSING_STEP,
    PhaseTransportSegment,
    _composed_metric_tail,
    _hairer_series,
    _metric_truncated,
    _series_context,
    _variable,
    enclose_phase_transport_cell,
    metric_series,
    original_dense_solution,
    phase_transport_segments,
    subgap_phase_transport_bound,
)
from recursive_horizons.nsc_paired_horizon_preparation import PairedHorizonSeedMap


ENERGY = 0.0013248831260437577
MASS = np.pi / 2
ANGULAR = np.sqrt(5.)
PREP = PairedHorizonSeedMap(
    1.9006916054701435, .23832579963401956, 3.973074368754331, 1e-10, 1e-9)
HORIZON = PREP.horizon_rho


@pytest.fixture(scope='module')
def group14_mode():
    return solve_jost(PREP.background, ENERGY, MASS, ANGULAR)


def _outer_cell(y_start=None, width=0.01):
    if y_start is None:
        y_start = math.log(60. - HORIZON)
    return PhaseTransportSegment(
        y_start, y_start - width, 0., 0., np.zeros(6))


def test_original_dense_is_odesolution_named_dense(group14_mode):
    dense = original_dense_solution(group14_mode.run)
    assert group14_mode.run.sol.__code__.co_freevars == ('dense',)
    assert isinstance(dense, OdeSolution)
    assert dense is group14_mode.run.sol.__closure__[0].cell_contents
    assert len(dense.interpolants) == len(group14_mode.run.t) - 1


def test_next_interpolant_y_old_is_the_exact_numerical_endpoint(group14_mode):
    dense = original_dense_solution(group14_mode.run)
    for current, following in zip(dense.interpolants, dense.interpolants[1:]):
        np.testing.assert_array_equal(current.y_old + current.F[0], following.y_old)


def test_trajectory_segment_basis_matches_scipy_interpolant(group14_mode):
    dense = original_dense_solution(group14_mode.run)
    interpolant = dense.interpolants[10]
    following = dense.interpolants[11]
    segment = TrajectorySegment(
        interpolant.t_old, interpolant.t, interpolant.y_old, following.y_old,
        interpolant.F[1:])
    for fraction in (0., .2, .7, 1.):
        node = interpolant.t_old + fraction * interpolant.h
        phase = np.array([segment.evaluate(fraction)[0], interpolant(node)[0]])
        np.testing.assert_allclose(phase[0], phase[1], atol=2e-15, rtol=0)
        if fraction in (0., 1.):
            np.testing.assert_array_equal(
                segment.evaluate(fraction), interpolant(node) if fraction == 0.
                else following.y_old)


def test_hairer_series_matches_dense_phase(group14_mode):
    dense = original_dense_solution(group14_mode.run)
    interpolant = dense.interpolants[3]
    following = dense.interpolants[4]
    segment = TrajectorySegment(
        interpolant.t_old, interpolant.t, interpolant.y_old, following.y_old,
        interpolant.F[1:])
    with _series_context(128, 7):
        start = float(interpolant.y_old[0].real)
        end = float(following.y_old[0].real)
        series0 = _hairer_series(start, end, interpolant.F[1:, 0].real,
                                 _variable(arb(0), 7), 7)
        series1 = _hairer_series(start, end, interpolant.F[1:, 0].real,
                                 _variable(arb(1), 7), 7)
        assert series0[0].contains(arb(start))
        assert series1[0].contains(arb(end))
        for fraction in (0., .25, 1.):
            series = _hairer_series(
                start, end, interpolant.F[1:, 0].real, _variable(arb(fraction), 7), 7)
            numeric = float(segment.evaluate(fraction)[0].real)
            assert abs(series[0] - arb(numeric)) < arb('1e-15')


def test_phase_rhs_times_offset_is_f_y():
    for sign in (-1, 1):
        radius, energy, mass, theta = 10., .4, 1.2, -.31
        angular = sign * 2.1
        metric = 1 - 3 * ((1 + radius**2) * np.arctan(1 / radius) - radius)
        v1 = angular * np.sqrt(metric) / np.sqrt(1 + radius**2)
        v2 = mass * np.sqrt(metric)
        phase_r = 2 / metric * (v1 * np.cos(theta) + v2 * np.sin(theta) - energy)
        offset = radius - HORIZON
        y = np.log(offset)
        f_y = (2 * np.exp(y) / metric * (
            np.sqrt(metric) * (angular / np.sqrt(1 + radius**2) * np.cos(theta)
                               + mass * np.sin(theta)) - energy))
        assert abs(offset * phase_r - f_y) < 1e-14
        assert abs(np.exp(y) * phase_r - f_y) < 1e-14


def test_metric_tail_composed_with_inv_r_contains_closed_form_y_derivatives():
    with _series_context(256, 8):
        y0 = (arb(10) - arb(HORIZON)).log()
        step = arb('-0.05')
        x = _variable(arb(0), 8)
        radius = arb(HORIZON) + (y0 + step * x).exp()
        inv_r = radius.inv()
        enclosed = metric_series(inv_r, 8, 8)
        exact = 1 - 3 * ((1 + radius * radius) * radius.inv().atan() - radius)
        for k in range(9):
            assert enclosed[k].contains(exact[k])


def test_metric_u_derivatives_are_not_used_as_y_derivatives():
    # Tail in 1/r starts at u^3 for one metric term. Its y-Taylor coefficients
    # are the Faà-di-Bruno composition, not the u-Taylor coefficients themselves.
    with _series_context(192, 6):
        y0 = (arb(8) - arb(HORIZON)).log()
        x = _variable(arb(0), 6)
        radius = arb(HORIZON) + (y0 + x).exp()
        inv_r = radius.inv()
        truncated = _metric_truncated(inv_r, 1, 6)
        composed = _composed_metric_tail(inv_r, 1, 6)
        exact = 1 - 3 * ((1 + radius * radius) * radius.inv().atan() - radius)
        u_max = inv_r[0].abs_upper()
        majorant = arb(6) * _variable(u_max, 6)**3 / (1 - _variable(u_max, 6))
        tail = exact - truncated
        assert abs(tail[1]) < arb(majorant[1]).abs_lower() / 10
        assert not composed[1].contains(majorant[1])
        for k in range(7):
            assert (truncated[k] + composed[k]).contains(exact[k])


def test_exact_zero_phase_cell_has_zero_defect():
    cell = enclose_phase_transport_cell(
        _outer_cell(), 0., 1., 0., HORIZON)
    assert cell['defect_upper'].is_zero()
    assert cell['phase_derivative_lower'] > 0
    assert cell['step'] < 0


def test_mutated_interpolant_is_rejected_or_has_larger_defect():
    exact = enclose_phase_transport_cell(_outer_cell(), 0., 1., 0., HORIZON)
    cell = _outer_cell()
    mutated = PhaseTransportSegment(
        cell.y_start, cell.y_end, 0., 0., np.array([1e-3, 0, 0, 0, 0, 0], float))
    try:
        result = enclose_phase_transport_cell(mutated, 0., 1., 0., HORIZON)
    except ArithmeticError:
        return
    assert result['defect_upper'] > exact['defect_upper']
    assert result['defect_upper'] > arb('1e-8')


def test_representative_group14_finite_band_bound(group14_mode):
    mass = arb(MASS).union(arb.pi() / 2)
    angular = arb(ANGULAR).union(arb(5).sqrt())
    result = subgap_phase_transport_bound(
        ENERGY, mass, angular, PREP.background, mode=group14_mode)
    error = restored_upper(result['phase_error_inner_upper'])
    inner = restored_upper({
        'mantissa': result['inner_radius_certified_lower']['mantissa'],
        'exponent': result['inner_radius_certified_lower']['exponent'],
    })
    offset = restored_upper({
        'mantissa': result['inner_offset_certified_lower']['mantissa'],
        'exponent': result['inner_offset_certified_lower']['exponent'],
    })
    assert result['certified_finite_band'] is True
    assert result['all_exterior_radii'] is False
    assert result['horizon_sewing_error'] is None
    assert result['physical_rho1_source_error'] is None
    assert result['missing_step'] == MISSING_STEP
    assert result['outer_radius'] == group14_mode.outer_radius
    assert result['cells'] > 100
    assert error < arb('1e-8')
    assert inner >= 8
    assert offset >= arb('1e-4')
    assert restored_upper(result['initializer_phase_error_upper']) <= error


def test_negative_angular_sign_has_a_finite_outer_cell_bound():
    mode = solve_jost(PREP.background, ENERGY, MASS, -ANGULAR)
    dense = original_dense_solution(mode.run)
    segments = phase_transport_segments(dense, HORIZON, 8.)
    cell = enclose_phase_transport_cell(
        segments[0], ENERGY, MASS, -ANGULAR, HORIZON)
    assert cell['defect_upper'] < arb('1e-8')
    assert cell['phase_derivative_lower'] > 0


def test_inadmissible_inner_radius_cannot_emit_bound(group14_mode):
    with pytest.raises(ValueError, match='inner radius'):
        subgap_phase_transport_bound(
            ENERGY, MASS, ANGULAR, PREP.background, inner_radius=1.,
            mode=group14_mode)
    with pytest.raises(ValueError, match='subgap'):
        subgap_phase_transport_bound(2., MASS, ANGULAR, PREP.background)


def test_flint_precision_and_series_cap_are_restored(group14_mode):
    before = (ctx.prec, ctx.cap)
    enclose_phase_transport_cell(_outer_cell(), 0., 1., 0., HORIZON)
    assert (ctx.prec, ctx.cap) == before
    with pytest.raises(ValueError):
        subgap_phase_transport_bound(
            ENERGY, MASS, ANGULAR, PREP.background, inner_radius=1.,
            mode=group14_mode)
    assert (ctx.prec, ctx.cap) == before


def test_transport_rejects_degree_below_full_phase_polynomial():
    with pytest.raises(ValueError, match='Taylor degree'):
        enclose_phase_transport_cell(_outer_cell(), 0., 1., 0., HORIZON, degree=6)


def test_cell_width_keeps_exact_binary_endpoint_difference():
    cell = PhaseTransportSegment(2.**-4, -2.**-60, 0., 0., np.zeros(6))
    enclosed = enclose_phase_transport_cell(cell, 0., 1., 0., HORIZON)
    with ctx.workprec(192):
        exact = arb(cell.y_end)-arb(cell.y_start)
        assert enclosed['step'] == exact
        assert not enclosed['step'].contains(arb(cell.y_end-cell.y_start))


def test_actual_initializer_is_bound_even_with_different_seed_order():
    mode = solve_jost(PREP.background, ENERGY, MASS, ANGULAR, outer_radius=60., order=4)
    proof = subgap_phase_transport_bound(ENERGY,MASS,ANGULAR,PREP.background,mode=mode)
    actual = float(original_dense_solution(mode.run).interpolants[0].y_old[0].real)
    assert proof['actual_initializer_phase_hex'] == actual.hex()
    assert restored_upper(proof['initializer_phase_error_upper']) > arb('1e-14')


def test_wrong_radial_coordinate_origin_is_rejected(group14_mode):
    from recursive_horizons.nsc_unruh_state import ParentDirac
    changed = ParentDirac(HORIZON+1e-5, PREP.surface_gravity)
    with pytest.raises(ValueError, match='coordinate origin'):
        subgap_phase_transport_bound(ENERGY,MASS,ANGULAR,changed,mode=group14_mode)


def test_subdivision_bounds_the_same_cell_without_width_scaling(group14_mode):
    dense = original_dense_solution(group14_mode.run)
    cells = phase_transport_segments(dense,HORIZON,HORIZON+1.01e-4)
    last = cells[-1]
    coarse = enclose_phase_transport_cell(last,ENERGY,MASS,ANGULAR,HORIZON,
                                         degree=16,metric_terms=48,defect_subdivisions=1)
    fine = enclose_phase_transport_cell(last,ENERGY,MASS,ANGULAR,HORIZON,
                                       degree=16,metric_terms=48,defect_subdivisions=4)
    assert fine['defect_upper'] <= coarse['defect_upper']
    assert fine['step'] == coarse['step']
    assert fine['phase_derivative_lower'] > 0
