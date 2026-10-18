"""Bounded controls for the quiet wall-pinned Dirac time chart.

These checks stay on manufactured points and one local jet. They do not call
the nine-cut reference report or search for the coordinate fold.
"""
import math

import numpy as np
import pytest
import sympy as s

from recursive_horizons.nsc_atlas_dirac_slice import (
    constant_rate_bridge, outgoing_slice, source_anchor_slice,
)
from recursive_horizons.nsc_atlas_dirac_time_chart import (
    ChartDomainError, cap_chart_fields, chart_fields, fixed_wall_transmission,
    initial_bridge_points, initial_outgoing_points, make_time_chart, wall_time_data,
)


@pytest.fixture(scope='module')
def chart():
    return make_time_chart()


def test_initial_slice_velocity_and_signed_descending_slope(chart):
    q = np.array([0.5 * chart['q_wall'], 0.0, chart['q_horizon']])
    assert q[0] < 0.0 < q[2]
    base = initial_bridge_points(q, chart)
    reference = constant_rate_bridge(q, source_anchor_slice())
    assert base['R0'] == pytest.approx(reference['R'])
    assert base['R0_y'] == pytest.approx(reference['Rprime'])
    assert base['V0'] == pytest.approx(reference['V'])
    assert base['R0_y'][0] < 0.0
    assert base['R0_y'][1] == pytest.approx(0.0, abs=1e-12)
    assert base['R0'][2] == pytest.approx(0.8)
    fields = chart_fields(0.45, base, chart)
    assert fields['a'] == pytest.approx(1.0)
    assert fields['N'] == pytest.approx(base['N0'])
    assert fields['shift'] == pytest.approx(base['b0'])
    assert fields['U_t'] == pytest.approx(base['Xdot0_U'])
    assert fields['R_t'] == pytest.approx(base['Xdot0_R'])
    assert fields['Kll'] == pytest.approx(-base['V_R'], abs=1e-9)
    assert fields['mass'] == pytest.approx(0.4)
    assert fields['hamiltonian_residual'] == pytest.approx(0.0, abs=1e-9)
    assert fields['momentum_residual'] == pytest.approx(0.0, abs=1e-9)
    assert fields['jacobian'] == pytest.approx(fields['N'] * fields['a'])
    assert np.min(fields['spatial_norm']) > 0.0
    assert np.all(np.isfinite(fields['U_y']))
    assert max(abs(chart['de_t0'][0]), abs(chart['de_t0'][1])) < 1e-9
    outgoing = initial_outgoing_points(np.array([1.0, 1.1, 2.0]), chart)
    vacuum = outgoing_slice(np.array([1.1, 2.0]), source_anchor_slice())
    assert outgoing['V0'][1:] == pytest.approx(vacuum['V'])
    assert outgoing['R0_y'][1:] == pytest.approx(vacuum['Rprime'])
    assert outgoing['V0'][2] == pytest.approx(0.0)
    assert outgoing['R0_y'][2] == pytest.approx(math.sqrt(1.0 - 0.8 / 2.0))
    end = initial_bridge_points(chart['q_taper'], chart)
    assert outgoing['R0'][0] == pytest.approx(end['R0'])
    assert outgoing['R0_y'][0] == pytest.approx(end['R0_y'])
    assert outgoing['y'][0] == pytest.approx(end['y'])
    assert outgoing['U0'][0] == pytest.approx(end['U0'], abs=1e-12)
    outside = chart_fields(0.45, outgoing, chart)
    assert outside['mass'] == pytest.approx(0.4)
    assert outside['hamiltonian_residual'] == pytest.approx(0.0, abs=1e-9)
    assert outside['dirac']['speed'] == pytest.approx(outside['c'])
    assert outside['dirac']['derivative_companion'] == pytest.approx(0.5 * outside['c_y'])
    assert outside['dirac']['shift'] == pytest.approx(outside['shift'])
    assert outside['dirac']['shift_companion'] == pytest.approx(0.5 * outside['shift_y'])


def test_one_local_time_and_y_finite_difference(chart):
    q, radius, step_q, step_r = -0.1, 0.48, 1e-5, 1e-5

    def sample(parameter, wall):
        return chart_fields(wall, initial_bridge_points(parameter, chart), chart)

    centre, ahead, behind = sample(q, radius), sample(q + step_q, radius), sample(q - step_q, radius)
    dy = (ahead['y'] - behind['y']) / (2.0 * step_q)
    assert (ahead['U'] - behind['U']) / (2.0 * step_q) == pytest.approx(centre['U_y'] * dy, abs=1e-7)
    assert (ahead['R'] - behind['R']) / (2.0 * step_q) == pytest.approx(centre['R_y'] * dy, abs=1e-7)
    assert (ahead['N'] - behind['N']) / (2.0 * step_q) == pytest.approx(centre['N_y'] * dy, abs=1e-6)
    later, earlier = sample(q, radius + step_r), sample(q, radius - step_r)
    dt = chart['q0'] * 2.0 * step_r
    assert (later['U'] - earlier['U']) / dt == pytest.approx(centre['U_t'], abs=1e-6)
    assert (later['R'] - earlier['R']) / dt == pytest.approx(centre['R_t'], abs=1e-6)
    assert centre['jacobian'].item() > 0.0
    assert centre['spatial_norm'].item() > 0.0
    assert centre['N'].item() > 0.0


@pytest.mark.parametrize('radius', [0.45, 0.5, 0.55])
def test_future_wall_metric_israel_and_spinor(chart, radius):
    row = fixed_wall_transmission(radius, chart, np.array([0.3 - 0.2j, 0.1 + 0.4j]))
    for key in ('radial_length_error', 'lapse_error', 'shift_error', 'cap_lapse_error',
                'cap_shift_error', 'cap_radius_error', 'tangent_U_error', 'tangent_R_error',
                'position_U_error', 'position_R_error', 'israel_Ktau_error', 'beta_plus_error',
                'spinor_parallel_error'):
        assert row[key] == pytest.approx(0.0, abs=1e-9)
    assert row['surface_spinor_force'] == 0.0
    assert row['nested_Psi_over_Phi'] is None
    assert row['cap_Kll'] == 0.0 and row['cap_a_y'] == 0.0
    if radius == 0.45:
        assert row['exterior_Kll'] == pytest.approx(0.0, abs=1e-9)
    else:
        assert abs(row['exterior_Kll']) > 0.1


def test_shifted_half_density_green_identity():
    coordinate = s.symbols('y', real=True)
    speed = 1 + coordinate**2 / 5
    shift = s.Rational(1, 5) * s.cos(coordinate)
    angular = (1 + coordinate) / 3
    alpha = s.Matrix([[0, -s.I], [s.I, 0]])
    sigma1 = s.Matrix([[0, 1], [1, 0]])
    left = s.Matrix([1 + s.I * coordinate, coordinate**2])
    right = s.Matrix([coordinate**3, 2 - s.I * coordinate])

    def generator(spinor):
        return (-s.I * alpha * (speed * spinor.diff(coordinate) + s.diff(speed, coordinate) * spinor / 2)
            + s.I * (shift * spinor.diff(coordinate) + s.diff(shift, coordinate) * spinor / 2)
            + sigma1 * angular * spinor)

    density = (s.conjugate(left).T * generator(right) - s.conjugate(generator(left)).T * right)[0]
    current = -s.I * (s.conjugate(left).T * (speed * alpha - shift * s.eye(2)) * right)[0]
    assert s.simplify(density - s.diff(current, coordinate)) == 0


def test_manufactured_reversed_jacobian_is_rejected(chart):
    base = initial_outgoing_points(2.0, chart)
    assert base['y'].item() > 1.0
    assert base['S'].item() == 0.0
    reversed_velocity = dict(base)
    reversed_velocity['Xdot0_U'] = -base['Xdot0_U']
    reversed_velocity['Xdot0_R'] = -base['Xdot0_R']
    with pytest.raises(ChartDomainError) as caught:
        chart_fields(0.45, reversed_velocity, chart)
    assert caught.value.details['min_jacobian'] < 0.0
    assert caught.value.details['min_spatial_norm'] > 0.0
    assert caught.value.details['radius'] == pytest.approx(0.45)


def test_domain_and_nonfinite_inputs_are_rejected(chart):
    base = initial_bridge_points(0.0, chart)
    with pytest.raises(ValueError):
        initial_bridge_points(np.nan, chart)
    with pytest.raises(ValueError):
        initial_bridge_points(chart['q_wall'] - 0.05, chart)
    with pytest.raises(ValueError):
        initial_outgoing_points(np.array([np.nan]), chart)
    with pytest.raises(ValueError):
        initial_outgoing_points(0.9, chart)
    with pytest.raises(ValueError):
        initial_outgoing_points(4.1, chart)
    with pytest.raises(ValueError):
        wall_time_data(np.inf, chart)
    with pytest.raises(ValueError):
        wall_time_data(0.4, chart)
    with pytest.raises(ValueError):
        chart_fields(0.6, base, chart)
    with pytest.raises(ValueError):
        initial_bridge_points(0.0, chart, order=16)
    with pytest.raises(ValueError):
        cap_chart_fields(0.5, 0.05, chart)
    with pytest.raises(ValueError):
        cap_chart_fields(0.5, -2.0, chart)
    with pytest.raises(ValueError):
        fixed_wall_transmission(0.45, chart, np.ones(3))
    centre = cap_chart_fields(0.5, -math.asin(0.5), chart)
    assert centre['R'].item() == pytest.approx(0.0, abs=1e-15)
    assert centre['angular_coefficient_evaluated'] is False
    assert 'dirac' not in centre
    cap_wall = cap_chart_fields(0.5, 0.0, chart)
    assert cap_wall['R'].item() == pytest.approx(0.5)
    assert cap_wall['a'].item() == 1.0
    assert 'dirac' in cap_wall


def test_one_coordinate_clock_quadrature_comparison(chart):
    bridge_32 = initial_bridge_points(chart['q_taper'], chart, order=32)['U0'].item()
    bridge_64 = initial_bridge_points(chart['q_taper'], chart, order=64)['U0'].item()
    wall_32 = wall_time_data(0.5, chart, order=32)
    wall_64 = wall_time_data(0.5, chart, order=64)
    assert abs(bridge_32 - bridge_64) < 1e-10
    assert abs(wall_32['U'].item() - wall_64['U'].item()) < 1e-10
    assert abs(wall_32['Tcap'].item() - wall_64['Tcap'].item()) < 1e-10
    assert abs(wall_32['tau'].item() - wall_64['tau'].item()) < 1e-10
    assert wall_32['N_wall'].item() > 0.0
