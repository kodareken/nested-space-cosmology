"""Synthetic KS history-evaluator controls. Not a source campaign or gate claim."""
import time
from dataclasses import replace

import numpy as np
import pytest

from recursive_horizons.nsc_compatible_history_geometry import CompatibleIncomingMetric
from recursive_horizons.nsc_ks_batched_constraints import map_negative_batch
from recursive_horizons.nsc_ks_history_evaluator import (
    KSHistoryEvaluator, fixed_target_union, group_signed_operator_families,
    history_operator_cache_key, signed_symmetry_applies, subset_node_indices,
)
from recursive_horizons.nsc_ks_profile_identity import profile_identity
from recursive_horizons.nsc_ks_source_envelope import (
    computational_z_grid, physical_incoming_interval, usual_axial_support,
)
from recursive_horizons.nsc_local_history_newton import (
    HistoryCollocation, bind_history_evaluator_receipt,
)
from recursive_horizons.nsc_local_incoming_constraints import local_error_budget
from recursive_horizons.nsc_local_incoming_family import LocalAxialFunction, LocalIncomingFamily
from test_nsc_ks_batched_constraints import (
    ANGULAR, MASS, RHO_UP, bind, channel, coefficients_for, make_panel,
    negative_source,
)
from test_nsc_ks_difference_envelope import inputs as difference_inputs


CPU = {'used': 0.0}
SOLVER = dict(rtol=2e-10, atol=2e-12, max_step=0.002)
ENERGY_INTERVAL = (0.3, 1.2)
DEGREE = 1
BASELINE = np.array([0.31, -0.17])


@pytest.fixture(autouse=True)
def _cpu_budget():
    start = time.process_time()
    yield
    CPU['used'] += time.process_time() - start


def star_coefficients(scale=1.0):
    coeff = np.zeros((2, 8))
    coeff[0, :4] = scale * np.array([0.012, -0.006, 0.003, -0.001])
    coeff[1, :4] = scale * np.array([0.020, 0.008, -0.004, 0.002])
    return coeff


def make_family(scale=1.0):
    return LocalIncomingFamily(star_coefficients(scale))


def nodes_for(family=None):
    family = make_family() if family is None else family
    collocation = HistoryCollocation(8)
    return collocation.solve_nodes(family), collocation.verification_nodes(family)


def signed_batches(energy=0.7, weight=0.2, seed=33):
    panel = make_panel(name='synthetic-history', energies=(energy,), weights=(weight,), seed=seed)
    positive = bind(panel)
    negative = map_negative_batch(positive, negative_source(panel))
    ch_p = channel(14, MASS, ANGULAR, 1)
    ch_m = channel(14, MASS, ANGULAR, -1)
    return panel, positive, negative, ch_p, ch_m


def make_evaluator(entries, family=None, **kwargs):
    family = make_family() if family is None else family
    solve, verify = nodes_for(family)
    settings = dict(
        baseline_gradient=BASELINE, coefficients=coefficients_for(),
        energy_interval=ENERGY_INTERVAL, interpolant_degree=DEGREE,
        z_grid=computational_z_grid(16), solve_nodes=solve, verification_nodes=verify,
        axial_support=usual_axial_support(), ell0_channels=(channel(13, MASS, 0., 1, copies=1, degeneracy=4),),
        error_budget=local_error_budget(), **SOLVER)
    settings.update(kwargs)
    return KSHistoryEvaluator(entries, **settings)


def test_signed_grouping_reuses_operator_only_for_actual_s3_partners():
    _, positive, negative, ch_p, ch_m = signed_batches()
    other = bind(make_panel(name='opposite-angular-positive', sign=-1, energies=(0.7,),
                            weights=(0.2,), seed=8))
    ch_other = channel(14, MASS, ANGULAR, -1)
    paired = group_signed_operator_families((
        ('p', positive, ch_p), ('m', negative, ch_m),
    ))
    assert len(paired) == 1
    assert len(paired[0]['applies']) == 1
    assert len(paired[0]['maps']) == 1
    assert paired[0]['maps'][0][1] == 'm'
    assert signed_symmetry_applies(positive, negative)
    assert not signed_symmetry_applies(positive, other)
    split = group_signed_operator_families((
        ('p', positive, ch_p), ('o', other, ch_other),
    ))
    assert len(split) == 2
    assert all(len(group['maps']) == 0 for group in split)
    same_channel = group_signed_operator_families((
        ('a', bind(make_panel(name='a', energies=(0.4,), weights=(0.1,), seed=1)), ch_p),
        ('b', bind(make_panel(name='b', energies=(0.8,), weights=(0.2,), seed=2)), ch_p),
    ))
    assert len(same_channel) == 1
    assert len(same_channel[0]['applies']) == 2


def test_cache_key_binds_profile_directions_grid_and_channel_not_amplitudes_alone():
    family = difference_inputs(0.0)[2]
    direction = family.directions[0]
    changed = replace(
        direction, w=LocalAxialFunction((0., 1., 0.2), direction.w.center))
    other = CompatibleIncomingMetric(family.amplitudes, (changed,))
    assert tuple(float(v) for v in family.amplitudes) == tuple(float(v) for v in other.amplitudes)
    kwargs = dict(
        mass=MASS, angular=ANGULAR, rho_up=RHO_UP, z_grid=computational_z_grid(16),
        target_z=np.linspace(*physical_incoming_interval(), 5), energy_interval=ENERGY_INTERVAL,
        degree=DEGREE, axial_support=usual_axial_support(), **SOLVER)
    key = history_operator_cache_key(family, **kwargs)
    assert key != history_operator_cache_key(other, **kwargs)
    assert key[0] == profile_identity(family)
    shifted_grid = dict(kwargs)
    shifted_grid['z_grid'] = computational_z_grid(32)
    assert key != history_operator_cache_key(family, **shifted_grid)
    shifted_rho = dict(kwargs)
    shifted_rho['rho_up'] = 1.04
    assert key != history_operator_cache_key(family, **shifted_rho)
    shifted_mass = dict(kwargs)
    shifted_mass['mass'] = MASS + 0.1
    assert key != history_operator_cache_key(family, **shifted_mass)
    assert key != history_operator_cache_key(family, **kwargs, integrator='cf4')
    assert key != history_operator_cache_key(family, **kwargs, step_control='joint')
    assert key != history_operator_cache_key(
        family, **kwargs, tangent_rtol=SOLVER['rtol'] / 2,
        tangent_atol=SOLVER['atol'])


def test_absent_nodes_rejected_and_union_is_exact_subset_owner():
    solve, verify = nodes_for()
    union = fixed_target_union(solve, verify)
    assert np.all(np.diff(union) > 0)
    assert subset_node_indices(union, solve).shape == (len(solve),)
    assert subset_node_indices(union, verify).shape == (len(verify),)
    shifted = np.array(solve, dtype=float, copy=True)
    shifted[0] += 3e-12
    with pytest.raises(ValueError, match='absent'):
        subset_node_indices(union, shifted)
    with pytest.raises(ValueError, match='absent'):
        subset_node_indices(union, np.array([union[0] + 1.0]))


def test_node_subsetting_cache_reuse_and_receipt_contract():
    _, positive, negative, ch_p, ch_m = signed_batches()
    family = make_family()
    evaluator = make_evaluator((('p', positive, ch_p), ('m', negative, ch_m)), family=family)
    solve, verify = nodes_for(family)
    first = evaluator.evaluate(family, z=solve)
    second = evaluator.evaluate(family, z=solve)
    verify_raw = evaluator.evaluate(family, z=verify)
    union_raw = evaluator.evaluate(family)
    assert evaluator.operator_evolution_count == 1
    assert evaluator.signed_operator_reuse[0]['signed_symmetry_maps'] == 1
    assert np.array_equal(first['z'], solve)
    assert np.array_equal(verify_raw['z'], verify)
    assert np.array_equal(union_raw['z'], evaluator.target_nodes)
    assert first['history_jacobian'].shape == (16, len(solve), 2)
    assert np.max(np.abs(first['history_jacobian'])) > 0
    np.testing.assert_array_equal(first['action_gradient'], second['action_gradient'])
    receipt = bind_history_evaluator_receipt(first, family, expected_nodes=solve)
    assert receipt.full_retarded_state_derivative is True
    assert receipt.physical_constraint_status == 'OPEN'
    assert receipt.missing_error_budget is None
    assert first['profile_identity'] == profile_identity(family)
    assert first['history_identity'] == profile_identity(family)
    cached = np.array(second['action_gradient'], copy=True)
    first['action_gradient'] = np.full_like(first['action_gradient'], 123.0)
    np.testing.assert_array_equal(
        evaluator.evaluate(family, z=solve)['action_gradient'], cached)
    with pytest.raises(ValueError, match='absent'):
        evaluator.evaluate(family, z=solve + 3e-12)
    other = make_family(1.2)
    evaluator.evaluate(other, z=solve)
    assert evaluator.operator_evolution_count == 2


def test_source_mutation_changes_identity_and_is_rejected_across_bindings():
    family = make_family()
    solve, _ = nodes_for(family)
    _, a_pos, a_neg, ch_p, ch_m = signed_batches(seed=33)
    _, b_pos, b_neg, _, _ = signed_batches(seed=99)
    assert a_pos.source.digest == b_pos.source.digest
    assert a_pos.preparation_digest != b_pos.preparation_digest
    first = make_evaluator((('p', a_pos, ch_p), ('m', a_neg, ch_m)), family=family)
    second = make_evaluator((('p', b_pos, ch_p), ('m', b_neg, ch_m)), family=family)
    raw_a = first.evaluate(family, z=solve)
    raw_b = second.evaluate(family, z=solve)
    assert raw_a['source_identity'] != raw_b['source_identity']
    assert np.max(np.abs(raw_a['action_gradient'] - raw_b['action_gradient'])) > 0
    bound = bind_history_evaluator_receipt(raw_a, family, expected_nodes=solve)
    with pytest.raises(ValueError, match='source identity'):
        bind_history_evaluator_receipt(
            raw_b, family, expected_nodes=solve,
            expected_source_identity=bound.source_identity)
    _, c_pos, c_neg, _, _ = signed_batches(energy=0.9, seed=33)
    assert c_pos.source.digest != a_pos.source.digest
    raw_c = make_evaluator((('p', c_pos, ch_p), ('m', c_neg, ch_m)), family=family).evaluate(
        family, z=solve)
    with pytest.raises(ValueError, match='source identity'):
        bind_history_evaluator_receipt(
            raw_c, family, expected_nodes=solve,
            expected_source_identity=bound.source_identity)
    assert 13 in raw_a['scope']['analytic_zero_response_groups']
    np.testing.assert_array_equal(raw_a['baseline_action_gradient'][0], BASELINE)


def test_explicit_signed_source_weights_and_no_factor_two():
    _, positive, negative, ch_p, ch_m = signed_batches()
    family = make_family()
    solve, _ = nodes_for(family)
    both = make_evaluator((('p', positive, ch_p), ('m', negative, ch_m)), family=family)
    pos_only = make_evaluator((('p', positive, ch_p),), family=family,
                              ell0_channels=())
    raw_both = both.evaluate(family, z=solve)
    raw_pos = pos_only.evaluate(family, z=solve)
    assert np.array_equal(negative.source.column_weights, positive.source.column_weights)
    assert abs(negative.source.covariance[0, 1]) > 0
    complement = np.eye(negative.source.covariance.shape[0]) - positive.source.covariance.conj()
    assert np.max(np.abs(negative.source.covariance - complement)) < 3e-13
    pos_id = next(key for key, rec in raw_both['batch_records'].items() if rec['energy_sign'] > 0)
    neg_id = next(key for key, rec in raw_both['batch_records'].items() if rec['energy_sign'] < 0)
    assert raw_both['batch_records'][neg_id]['source_digest'] == negative.source.digest
    assert raw_both['batch_records'][pos_id]['source_digest'] == positive.source.digest
    both_corr = sum(raw_both['batch_corrections'].values())
    pos_corr = sum(raw_pos['batch_corrections'].values())
    assert np.max(np.abs(both_corr - 2 * pos_corr)) > 1e-8
    assert raw_both['scope']['positive_result_doubled'] is False
    assert raw_both['scope']['energy_folding_factor'] == 1
    assert raw_both['full_retarded_state_derivative'] is True
    assert raw_both['physical_EXISTENCE_certificate'] is False
    assert raw_both['error_budget']['changed_history_finite_energy'] is None
    assert raw_both['missing_error_budget'] is None
    assert raw_both['scope']['integrator'] == 'dop853'
    assert raw_both['scope']['step_control'] == 'primal'
    assert raw_both['scope']['tangent_rtol'] == SOLVER['rtol']
    assert raw_both['scope']['tangent_atol'] == SOLVER['atol']


def test_integrator_and_step_control_are_validated_before_evolution():
    _, positive, negative, ch_p, ch_m = signed_batches()
    entries = (('p', positive, ch_p), ('m', negative, ch_m))
    with pytest.raises(ValueError, match='integrator'):
        make_evaluator(entries, integrator='euler')
    with pytest.raises(ValueError, match='step_control'):
        make_evaluator(entries, step_control='shared')
    with pytest.raises(ValueError, match='separate tangent'):
        make_evaluator(entries, step_control='joint', tangent_rtol=1e-8)


def test_retarded_jacobian_matches_adapter_finite_difference():
    _, positive, negative, ch_p, ch_m = signed_batches()
    family = make_family()
    solve, _ = nodes_for(family)
    evaluator = make_evaluator((('p', positive, ch_p), ('m', negative, ch_m)), family=family)
    center = evaluator.evaluate(family, z=solve)
    h = 1e-4
    plus = make_family()
    minus = make_family()
    plus_coeff = np.array(plus.coefficients, float)
    minus_coeff = np.array(minus.coefficients, float)
    plus_coeff[0, 0] += h
    minus_coeff[0, 0] -= h
    plus = LocalIncomingFamily(plus_coeff)
    minus = LocalIncomingFamily(minus_coeff)
    fd = (evaluator.evaluate(plus, z=solve)['action_gradient']
          - evaluator.evaluate(minus, z=solve)['action_gradient']) / (2 * h)
    np.testing.assert_allclose(fd, center['history_jacobian'][0], atol=3e-6, rtol=0)
    geom = center['local_reference_gradient_tangent']
    assert geom.shape == center['history_jacobian'].shape
    assert np.max(np.abs(center['history_jacobian'] - geom)) > 0
    assert evaluator.operator_evolution_count == 3
    assert center['scope']['baseline_and_geometry_counted_once'] is True
    assert center['scope']['source_cutoff_coincidence_bridge'] == 'OPEN'


def test_frozen_c0_and_ell0_numerical_batch_rejected():
    _, positive, _, ch_p, _ = signed_batches()
    family = make_family()
    with pytest.raises(TypeError, match='KSUpstreamBatch'):
        make_evaluator(((object(), ch_p),), family=family)
    with pytest.raises(ValueError, match='Gamma_rest'):
        make_evaluator((('p', positive, ch_p),), family=family, Gamma_rest=1.0)
    ell0 = bind(make_panel(name='ell0', group=13, angular=0., energies=(0.7,), weights=(0.2,)))
    with pytest.raises(ValueError, match='analytically'):
        make_evaluator((('z', ell0, channel(13, MASS, 0., 1)),), family=family)


def test_zzz_cpu_budget_is_control_only():
    assert 0.0 < CPU['used'] < 30.0, CPU['used']


def test_caller_buffers_and_nested_receipts_cannot_change_cached_science():
    _, positive, negative, ch_p, ch_m = signed_batches()
    baseline = BASELINE.copy()
    coefficients = coefficients_for()
    expected_d = np.array(coefficients['D'], copy=True)
    evaluator = make_evaluator(((positive, ch_p), (negative, ch_m)),
                               baseline_gradient=baseline, coefficients=coefficients)
    owned = evaluator._entries[0][1]
    assert not np.shares_memory(owned.initial_columns, positive.initial_columns)
    assert not np.shares_memory(owned.source.covariance, positive.source.covariance)
    baseline[:] = 999
    coefficients['D'][:] = 999
    np.testing.assert_array_equal(evaluator._baseline, BASELINE)
    np.testing.assert_array_equal(evaluator._coefficients['D'], expected_d)
    raw = evaluator.evaluate(make_family())
    first_key = next(iter(raw['batch_records']))
    with pytest.raises(TypeError):
        raw['batch_records'][first_key]['preparation_digest'] = 'changed'
    with pytest.raises(TypeError):
        raw['scope']['analytic_zero_response_groups'][13]['baseline_retained'] = False
    np.testing.assert_array_equal(raw['baseline_action_gradient'][0], BASELINE)
