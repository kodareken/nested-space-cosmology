"""Numerical comparison controls: circular phase and unchanged preparation."""
import copy
import math

import pytest

from compare_nsc_atlas_gh_evolution import ERRORS, compare


def record(degree=8, phase=math.pi-.1, *, rtol=3e-8, ds=-.001):
    point = dict(observer_proper_duration=0, metric_valid=True,
        observer_volume_mode=-.035, observer_phase_wrapped=0,
        geodesic_observer_total_density_over_initial_vacuum=1.25,
        slice_proper_volume_over_initial_round_reference=1,
        **{name:.01 for name in ERRORS})
    final = dict(point, observer_proper_duration=4., observer_phase_wrapped=phase,
                 observer_volume_mode=point['observer_volume_mode']+ds)
    return dict(epsilon=.25, degree=degree, prescribed_GH_gauge='lapse_controlled',
        gauge_constraint_damping=.27, lapse_feedback_coefficient=.27,
        formulation='GH', symmetry='paired', source='same two modes',
        requested_observer_proper_endpoint=4., model_inputs={'sstar':-.035},
        status='endpoint_reached', snapshots=[point,final],
        temporal_tolerances=dict(rtol=rtol, atol=rtol/100, max_step=.05))


def test_principal_phase_and_exact_local_coupling():
    result = compare([record(), record(12, -math.pi+.1, ds=-.002)])
    change = result['spatial_changes'][0]
    assert change['signed_change']['observer_phase_principal_difference'] == pytest.approx(.2)
    assert result['cases'][1]['final']['g_KK_over_initial_local_value'] == pytest.approx(math.exp(.014))
    assert result['cases'][1]['final']['g_KK_fractional_change_from_initial_local_value'] == pytest.approx(math.expm1(.014))
    assert not result['physical_prediction_admitted']
    assert result['cases'][0]['work']['residual_over_initial_excitation'] is None


def test_time_control_is_separate_from_spatial_change_and_keeps_work():
    loose, fine = record(12), record(12, rtol=3e-9)
    fine['snapshots'][0]['excitation_normal_energy_over_initial_round_vacuum_energy'] = .05
    fine['snapshots'][-1]['excitation_work_balance_residual'] = .0001
    fine['passive_matter_work'] = dict(enabled=True, valid=True)
    result = compare([fine,loose])
    assert result['spatial_changes'] == []
    assert len(result['temporal_changes']) == 1
    assert result['temporal_changes'][0]['from_temporal_tolerances']['rtol'] == 3e-8
    assert result['cases'][0]['work']['residual_over_initial_excitation'] == pytest.approx(.002)


@pytest.mark.parametrize('change', ['epsilon', 'prescribed_GH_gauge', 'model_inputs', 'endpoint'])
def test_incompatible_families_are_not_merged(change):
    a,b = record(),record(12)
    if change == 'endpoint':
        b['requested_observer_proper_endpoint'] = 5
    elif change == 'model_inputs':
        b['model_inputs']['sstar'] = -.04
    elif change == 'epsilon':
        b[change] = 1
    else:
        b[change] = 'harmonic'
    with pytest.raises(ValueError, match='different preparations'):
        compare([a,b])


def test_incomplete_endpoints_and_duplicate_controls_are_rejected():
    a = record()
    with pytest.raises(ValueError, match='duplicate'):
        compare([a,copy.deepcopy(a)])
    a['status'] = 'incomplete'
    with pytest.raises(ValueError, match='reached proper endpoint'):
        compare([a])
