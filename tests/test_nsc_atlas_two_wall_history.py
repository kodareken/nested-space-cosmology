"""Independent host embedding, clock orientation and characteristic controls."""
import math

import pytest
from scipy.integrate import quad

from recursive_horizons.nsc_atlas_nested_spacetime import TwoWallHostHistory


def test_tortoise_partial_fractions_against_direct_static_integral():
    model = TwoWallHostHistory()
    for radius in (.08, .3, .8):
        direct = quad(lambda r: 1/(1-.02/r-r*r), .5, radius,
                      epsabs=2e-13, epsrel=2e-13)[0]
        assert model.tortoise(radius) == pytest.approx(direct, abs=2e-13)
    assert model.tortoise(.5) == 0.


def test_shared_initial_embedding_recovers_metric_and_normal_on_both_sides():
    model = TwoWallHostHistory()
    for branch, radius in (('inner', .1), ('inner', .76), ('inner', .8),
                           ('outer', .5), ('outer', .72), ('outer', model.rc),
                           ('outer', 1.4), ('outer', 1.5)):
        row = model.initial_event(branch, radius)
        lake = row['lake']
        for key in ('spatial_norm_residual','normal_norm_residual',
                    'orthogonality_residual','normal_radius_residual'):
            assert abs(lake[key]) < 3e-13
        assert lake['N_u'] < 0 and lake['N_w'] > 0
        if row['T'] is not None:
            assert row['T']-model.tortoise(radius) == pytest.approx(
                -math.log(abs(lake['u'])/model.ell)/model.kappa, abs=5e-15)


def test_bifurcation_chart_is_smooth_and_agrees_with_static_overlaps():
    model = TwoWallHostHistory()
    for y in (-math.sqrt(.1-model.rb), -.001, 0., .001,
              math.sqrt(.5-model.rb)):
        row = model.bridge_event(y)
        assert row['induced_dy_squared'] == pytest.approx(row['dl_dy']**2, rel=4e-13)
        assert abs(row['spatial_norm_residual']) < 4e-13
        assert abs(row['normal_radius_residual']) < 4e-13
        assert row['future_N_u'] < 0 and row['future_N_w'] > 0
        if y:
            side = 'inner' if y < 0 else 'outer'
            if abs(y) > .01:
                reference = model.initial_event(side, .1 if y < 0 else .5)
                assert row['u'] == pytest.approx(reference['lake']['u'], rel=3e-13)
                assert row['w'] == pytest.approx(reference['lake']['w'], rel=3e-13)
    throat = model.bridge_event(0.)
    assert throat['u'] == throat['w'] == 0.
    assert throat['induced_dy_squared'] == pytest.approx(2/model.kappa, rel=4e-13)


def test_initial_source_derivatives_and_arrival_map_from_independent_differences():
    model = TwoWallHostHistory()
    for name, radius in (('inner', .77), ('outer', 1.46)):
        h = 2e-5
        row = model.initial_event(name, radius)
        derivative = (model._initial_null(name,radius-2*h)
                      -8*model._initial_null(name,radius-h)
                      +8*model._initial_null(name,radius+h)
                      -model._initial_null(name,radius+2*h))/(12*h)
        Rp, V = row['Rprime'], row['V']
        expected = 1/(Rp*(Rp-V)) if name == 'inner' else -1/(Rp*(Rp+V))
        assert derivative == pytest.approx(expected, rel=2e-10)
        event = model.wall_event(name,.02)
        origin = model.characteristic_origin(name,event)
        dt = 2e-5
        left = model.characteristic_origin(name,model.wall_event(name,.02-dt))
        right = model.characteristic_origin(name,model.wall_event(name,.02+dt))
        numerical = 2*dt/(right['proper_label']-left['proper_label'])
        assert origin['arrival_d_tau_d_label'] == pytest.approx(numerical, rel=2e-8)


def test_shared_clock_gauge_and_length_scaling():
    unit = TwoWallHostHistory()
    large = TwoWallHostHistory(length_scale=7.)
    for name in ('inner','outer'):
        row = unit.wall_event(name,.04)
        scaled = large.wall_event(name,.28)
        for key in ('radius','T','q','null_coordinate'):
            assert scaled[key] == pytest.approx(7*row[key], rel=2e-12,abs=2e-13)
        for key in ('Rdot','host_beta'):
            assert scaled[key] == pytest.approx(row[key], rel=2e-12)
        for key in ('u','w'):
            assert scaled['lake'][key] == pytest.approx(7*row['lake'][key],rel=2e-12)
        origin = unit.characteristic_origin(name,row)
        big_origin = large.characteristic_origin(name,scaled)
        assert big_origin['proper_label'] == pytest.approx(7*origin['proper_label'],rel=2e-12)
        assert big_origin['arrival_d_tau_d_label'] == pytest.approx(
            origin['arrival_d_tau_d_label'],rel=2e-12)
    assert unit.initial_event('inner',.1)['T'] == 0.
    assert unit.initial_event('outer',.5)['T'] == 0.
    assert unit.initial_wall('inner')['T'] != 0.
    assert unit.initial_wall('outer')['T'] != 0.


def test_domain_rejections_preserve_original_inputs():
    model = TwoWallHostHistory()
    with pytest.raises(ValueError,match='admitted'):
        model.wall_event('inner',.041)
    with pytest.raises(ValueError,match='left taper'):
        model.initial_event('inner',.81)
    with pytest.raises(ValueError,match='right host'):
        model.initial_event('outer',1.51)
    with pytest.raises(ValueError,match='order'):
        TwoWallHostHistory(order=48)


def test_finite_shared_history_against_independent_evolution_and_causal_orientation():
    from check_nsc_atlas_two_wall_history import calculation
    row = calculation()
    assert row['counts']['scalar_geometry_or_RHS_evaluations'] <= 20000
    assert row['process_cpu_seconds'] < 60
    assert row['counts']['new_field_grids'] == 0
    assert not row['causal_access']['inner_segment_to_outer_segment']
    assert not row['causal_access']['outer_segment_to_inner_segment']
    for name, case in row['cases'].items():
        assert max(case['comparisons'].values()) < 2e-9
        assert max(case['residuals'].values()) < 2e-10
        assert case['endpoint_source']['radius'] < case['initial_event']['radius']
        arrivals = case['arrival_derivatives_at_initial_midpoint_endpoint']
        assert all(value > 0 for value in arrivals) if name == 'inner' else all(value < 0 for value in arrivals)
        endpoint = case['endpoint']
        assert endpoint['lake']['u_tau'] < 0 < endpoint['lake']['w_tau']
        certificate = case['whole_interval_certificate']
        assert certificate['source_endpoint_bracket_product'] < 0
        assert certificate['speed_squared_positive_lower_bound'] > 0
        assert certificate['subtracted_growth_remainder_upper_bound'] < 1
        assert endpoint['radius'] < certificate['analytic_radius_upper_bound']
        if name == 'inner':
            assert certificate['host_f_positive_lower_bound'] > 0
