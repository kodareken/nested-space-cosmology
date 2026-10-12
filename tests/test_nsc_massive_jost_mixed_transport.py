"""Signed growth, contraction, and tube controls for mixed phase transport."""
import math

import numpy as np
import pytest
from flint import arb, ctx

from recursive_horizons.nsc_ks_ball_trajectory import restored_upper
from recursive_horizons.nsc_massive_jost_mixed_transport import (
    MISSING_STEP,
    advance_mixed_tube,
    comparison_value,
    enclose_mixed_phase_transport_cell,
    phi_weight,
    signed_transport_image,
    subgap_mixed_phase_transport_bound,
)
from recursive_horizons.nsc_massive_jost_modes import solve_jost
from recursive_horizons.nsc_massive_jost_transport_bound import (
    PhaseTransportSegment,
    _contraction,
    enclose_phase_transport_cell,
    original_dense_solution,
    phase_transport_segments,
)
from recursive_horizons.nsc_paired_horizon_preparation import PairedHorizonSeedMap


ENERGY = 0.24867511687395624
MASS = np.pi / 2
ANGULAR = np.sqrt(5.)
PREP = PairedHorizonSeedMap(
    1.9006916054701435, .23832579963401956, 3.973074368754331, 1e-10, 1e-9)
HORIZON = PREP.horizon_rho
NEAR = HORIZON + 1.01e-4


@pytest.fixture(scope='module')
def row15_mode():
    return solve_jost(
        PREP.background, ENERGY, MASS, ANGULAR, order=8, radial_collar=1e-10,
        rtol=2e-13, atol=2e-15)


def _outer_cell(width=0.01, theta=0.):
    y_start = math.log(60. - HORIZON)
    return PhaseTransportSegment(
        y_start, y_start - width, theta, theta, np.zeros(6))


def _panel_integral(rate, panels=32):
    total = arb(0)
    rate = arb(rate)
    for index in range(panels):
        left, right = arb(index) / panels, arb(index + 1) / panels
        total += (-rate * left.union(right)).exp() * (right - left)
    return total


def test_phi_zero_near_zero_and_straddle_do_not_divide():
    with ctx.workprec(128):
        zero = phi_weight(arb(0))
        assert zero.is_finite() and zero.contains(arb(1))
        predecessor = _contraction(arb('1e-6'), arb('1e-8'), arb(0))
        assert not predecessor.is_finite()
        ein, eps = arb('1e-6'), arb('1e-8')
        image = signed_transport_image(ein, eps, arb(0))
        assert image.is_finite()
        assert (image - (ein + eps)).abs_upper() < arb('1e-24')
        for rate in (arb('1e-12'), arb('-1e-12'), arb('-1e-9').union(arb('2e-9'))):
            weight = phi_weight(rate)
            assert weight.is_finite() and weight > 0
            assert (weight - 1).abs_upper() < arb('1e-8')
        straddling = phi_weight(arb(-1).union(arb('0.5')))
        assert straddling.is_finite() and straddling.contains(arb(1))
        assert straddling.overlaps(phi_weight(arb(-1)))
        assert straddling.overlaps(phi_weight(arb('0.5')))


def test_phi_matches_enclosed_integral_and_signed_formula():
    with ctx.workprec(128):
        for rate in (arb('-0.4'), arb('1e-8'), arb('2.5'), arb(0)):
            assert _panel_integral(rate, 48).contains(phi_weight(rate))
        ein, eps, rate = arb('1e-6'), arb('2e-8'), arb('-0.3')
        image = signed_transport_image(ein, eps, rate)
        assert image.overlaps(_contraction(ein, eps, rate))
        assert image.lower() > ein
        contracted = signed_transport_image(ein, arb(0), arb('1.25'))
        assert contracted.upper() < ein
        assert contracted.overlaps(_contraction(ein, arb(0), arb('1.25')))


def test_constant_rate_peak_lies_at_an_endpoint():
    with ctx.workprec(128):
        samples = (
            (arb('1e-4'), arb('1e-6'), arb('0.8')),
            (arb('1e-4'), arb('1e-3'), arb('0.2')),
            (arb('1e-4'), arb(0), arb('-0.5')),
            (arb('1e-4'), arb('1e-5'), arb(0)),
            (arb('1e-4'), arb('1e-5'), arb('-1e-8').union(arb('2e-8'))),
        )
        for ein, eps, rate in samples:
            values = [comparison_value(ein, eps, rate, x) for x in (0, .25, .5, .8, 1)]
            peak = values[0].union(values[-1])
            for value in values:
                assert peak.contains(value) or value.upper() <= peak.upper()
            k = float(rate.lower())
            slope = float(eps) - k * float(ein)
            if slope > 0:
                assert values[-1].lower() >= values[0].upper() or values[-1].overlaps(values[0])
            elif slope < 0:
                assert values[0].lower() >= values[-1].upper() or values[0].overlaps(values[-1])


def test_continuous_tube_exit_is_rejected():
    cell = enclose_mixed_phase_transport_cell(
        _outer_cell(width=0.05), 0., -1., 0., HORIZON, tube='1e-4')
    assert cell['transport_rate_lower'] < 0
    assert cell['defect_upper'].is_zero()
    ein, tube = arb('2e-6'), arb('1e-4')
    assert comparison_value(ein, cell['defect_upper'], cell['transport_rate_lower'], 0) < tube
    assert comparison_value(ein, cell['defect_upper'], cell['transport_rate_lower'], 1) > tube
    with pytest.raises(ArithmeticError, match='tube'):
        advance_mixed_tube(ein, cell, tube)


def test_manufactured_cells_cover_growth_contraction_and_zero_rate():
    growing = enclose_mixed_phase_transport_cell(_outer_cell(), 0., -1., 0., HORIZON)
    shrinking = enclose_mixed_phase_transport_cell(_outer_cell(), 0., 1., 0., HORIZON)
    flat = enclose_mixed_phase_transport_cell(_outer_cell(), 0., 0., 0., HORIZON)
    near = enclose_mixed_phase_transport_cell(_outer_cell(), 0., 1e-8, 0., HORIZON)
    assert growing['defect_upper'].is_zero() and growing['phase_derivative_lower'] < 0
    assert shrinking['defect_upper'].is_zero() and shrinking['phase_derivative_lower'] > 0
    assert flat['defect_upper'].is_zero()
    assert flat['transport_rate_lower'].abs_upper() < arb('1e-8')
    assert near['transport_rate_lower'].abs_upper() < arb('1e-4')
    assert signed_transport_image(arb('1e-8'), arb(0), growing['transport_rate_lower']).lower() > arb('1e-8')
    assert signed_transport_image(arb('1e-8'), arb(0), shrinking['transport_rate_lower']).upper() < arb('1e-8')
    flat_image = signed_transport_image(arb('1e-6'), arb('1e-9'), flat['transport_rate_lower'])
    assert (flat_image - arb('1.001e-6')).abs_upper() < arb('1e-10')


def test_predecessor_still_rejects_nonpositive_derivative():
    cell = _outer_cell()
    with pytest.raises(ArithmeticError, match='monotonicity'):
        enclose_phase_transport_cell(cell, 0., -1., 0., HORIZON)
    signed = enclose_mixed_phase_transport_cell(cell, 0., -1., 0., HORIZON)
    assert signed['phase_derivative_lower'] < 0


def test_unresolved_metric_is_not_returned_as_a_bound():
    y_start = math.log(0.2)
    cell = PhaseTransportSegment(y_start, y_start - 0.01, 0., 0., np.zeros(6))
    before = (ctx.prec, ctx.cap)
    with pytest.raises(ArithmeticError, match='metric'):
        enclose_mixed_phase_transport_cell(cell, 0., 1., 0., 0.5)
    assert (ctx.prec, ctx.cap) == before


def test_cell_width_keeps_exact_binary_endpoint_difference():
    cell = PhaseTransportSegment(2.**-4, -2.**-60, 0., 0., np.zeros(6))
    enclosed = enclose_mixed_phase_transport_cell(cell, 0., 1., 0., HORIZON)
    with ctx.workprec(192):
        exact = arb(cell.y_end) - arb(cell.y_start)
        assert enclosed['step'] == exact
        assert not enclosed['step'].contains(arb(cell.y_end - cell.y_start))
        assert enclosed['transport_rate_lower'] == (
            (-enclosed['step']) * enclosed['phase_derivative_lower']).lower()


def test_subdivision_does_not_insert_a_second_width(row15_mode):
    dense = original_dense_solution(row15_mode.run)
    segments = phase_transport_segments(dense, HORIZON, NEAR)

    def lower_estimate(segment):
        y = 0.5 * (segment.y_start + segment.y_end)
        theta = 0.5 * (segment.theta_start + segment.theta_end)
        radius = HORIZON + math.exp(y)
        metric = 1 - 3 * ((1 + radius * radius) * math.atan(1 / radius) - radius)
        return (2 * math.exp(y) / math.sqrt(metric) * (
            MASS * math.cos(theta) - ANGULAR * math.sin(theta) / math.sqrt(1 + radius * radius)))
    segment = min(segments, key=lower_estimate)
    assert lower_estimate(segment) < 0
    kwargs = dict(degree=12, metric_terms=48)
    coarse = enclose_mixed_phase_transport_cell(
        segment, ENERGY, MASS, ANGULAR, HORIZON, defect_subdivisions=1, **kwargs)
    fine = enclose_mixed_phase_transport_cell(
        segment, ENERGY, MASS, ANGULAR, HORIZON, defect_subdivisions=4, **kwargs)
    assert fine['step'] == coarse['step']
    assert fine['defect_upper'] <= coarse['defect_upper']
    assert fine['phase_derivative_lower'] < 0 and coarse['phase_derivative_lower'] < 0
    assert fine['transport_rate_lower'].abs_upper() < 2 * coarse['transport_rate_lower'].abs_upper()
    assert coarse['transport_rate_lower'].abs_upper() < 2 * fine['transport_rate_lower'].abs_upper()


def test_original_row15_expanding_cell_is_enclosed(row15_mode):
    dense = original_dense_solution(row15_mode.run)
    segments = phase_transport_segments(dense, HORIZON, NEAR)
    def lower_estimate(segment):
        radius = HORIZON + math.exp(segment.y_end)
        metric = 1 - 3 * ((1 + radius * radius) * math.atan(1 / radius) - radius)
        theta = segment.theta_end
        return (2 * math.exp(segment.y_end) / math.sqrt(metric) * (
            MASS * math.cos(theta) - ANGULAR * math.sin(theta) / math.sqrt(1 + radius * radius)))
    segment = min(segments, key=lower_estimate)
    with pytest.raises(ArithmeticError, match='monotonicity'):
        enclose_phase_transport_cell(
            segment, ENERGY, MASS, ANGULAR, HORIZON, degree=12, metric_terms=48,
            defect_subdivisions=4)
    signed = enclose_mixed_phase_transport_cell(
        segment, ENERGY, MASS, ANGULAR, HORIZON, degree=12, metric_terms=48,
        defect_subdivisions=4)
    assert signed['phase_derivative_lower'] < 0
    assert signed['transport_rate_lower'] < 0
    assert signed['metric_lower'] > 0
    assert signed['defect_upper'] < arb('1e-5')
    assert signed['normalized_defect'] == 'theta_x-h*f_y'


def test_outer_initializer_transport_stays_in_tube(row15_mode):
    proof = subgap_mixed_phase_transport_bound(
        ENERGY, arb(MASS).union(arb.pi() / 2), arb(ANGULAR).union(arb(5).sqrt()),
        PREP.background, mode=row15_mode, inner_radius=30., degree=12, metric_terms=12)
    error = restored_upper(proof['phase_error_inner_upper'])
    assert proof['certified_finite_band'] is True
    assert proof['source_columns_replaced'] is False
    assert proof['horizon_sewing_error'] is None
    assert proof['physical_rho1_source_error'] is None
    assert proof['missing_step'] == MISSING_STEP
    assert proof['rate_signs']['expanding'] == 0
    assert proof['cells'] >= 2
    assert error < arb('1e-8')
    assert restored_upper(proof['initializer_phase_error_upper']) < arb('1e-8')


def test_precision_is_restored_and_invalid_input_raises():
    before = (ctx.prec, ctx.cap)
    enclose_mixed_phase_transport_cell(_outer_cell(), 0., 1., 0., HORIZON)
    assert (ctx.prec, ctx.cap) == before
    with pytest.raises(ValueError, match='inner radius'):
        subgap_mixed_phase_transport_bound(
            ENERGY, MASS, ANGULAR, PREP.background, inner_radius=1.)
    with pytest.raises(ValueError, match='subgap'):
        subgap_mixed_phase_transport_bound(2., MASS, ANGULAR, PREP.background)
    with pytest.raises(ValueError, match='Taylor degree'):
        enclose_mixed_phase_transport_cell(_outer_cell(), 0., 1., 0., HORIZON, degree=6)
    with pytest.raises(ArithmeticError, match='tube'):
        advance_mixed_tube(arb('1e-4'), {
            'defect_upper': arb(0), 'transport_rate_lower': arb('-2'),
        }, '1e-5')
    assert (ctx.prec, ctx.cap) == before


@pytest.mark.parametrize('missing', [None, False, True])
def test_missing_or_boolean_proof_inputs_are_not_zero(missing):
    with pytest.raises(ValueError, match='explicit'):phi_weight(missing)
    for values in ((missing,0,0),(0,missing,0),(0,0,missing)):
        with pytest.raises(ValueError, match='explicit'):signed_transport_image(*values)
    with pytest.raises(ValueError, match='explicit'):comparison_value(1,0,missing,.5)
    with pytest.raises(ValueError, match='explicit'):comparison_value(1,0,0,missing)
    with pytest.raises(ValueError, match='explicit'):
        advance_mixed_tube(0,{'defect_upper':arb(0),'transport_rate_lower':arb(0)},missing)
