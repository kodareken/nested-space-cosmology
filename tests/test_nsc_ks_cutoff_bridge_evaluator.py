"""Source-cutoff assembly controls; no complete spectrum or physical root."""
import numpy as np
import pytest

from recursive_horizons.nsc_ks_cutoff_bridge_evaluator import (
    KSCutoffBridgeEvaluator, signed_family_angular_square_weights)
from recursive_horizons.nsc_local_history_newton import bind_history_evaluator_receipt
from recursive_horizons.nsc_local_incoming_family import LocalIncomingFamily
from test_nsc_ks_history_evaluator import make_evaluator, make_family, signed_batches


def setup(family=None):
    _, positive, negative, cp, cn = signed_batches()
    raw = make_evaluator(((positive, cp), (negative, cn)), family=family)
    return raw, KSCutoffBridgeEvaluator(raw, phase_gauss_nodes=8)


def test_same_state_plus_edge_and_retarded_derivative_once_per_signed_family():
    family = make_family(.1)
    raw, owned = setup(family)
    base = raw.evaluate(family)
    result = owned.evaluate(family)
    np.testing.assert_array_equal(result['raw_action_gradient'], base['action_gradient'])
    np.testing.assert_allclose(result['action_gradient']-base['action_gradient'],
        result['source_cutoff_edge_gradient'], atol=1e-14, rtol=0)
    np.testing.assert_allclose(result['history_jacobian']-base['history_jacobian'],
        result['source_cutoff_edge_history_tangent'], atol=1e-10, rtol=0)
    np.testing.assert_allclose(result['source_cutoff_angular_square_weights'], [60., 60.], atol=2e-14, rtol=0)
    assert np.max(abs(result['source_cutoff_edge_gradient'])) > 1e-8
    assert np.max(abs(result['source_cutoff_edge_history_tangent'])) > 1e-8
    assert raw.operator_evolution_count == 1
    first = bind_history_evaluator_receipt(result, family, expected_nodes=owned.target_nodes)
    small = owned.evaluate(family, owned.target_nodes[::2])
    np.testing.assert_array_equal(small['source_cutoff_edge_gradient'], result['source_cutoff_edge_gradient'][::2])
    bind_history_evaluator_receipt(small, family,
        expected_nodes=owned.target_nodes[::2], expected_source_identity=first.source_identity)
    assert len(owned._phase_cache) == 1
    assert result['physical_EXISTENCE_certificate'] is False
    assert result['error_budget']['changed_history_tail'] is None
    assert result['scope']['new_action_term'] is False


def test_full_joined_jacobian_finite_difference_with_w_and_U_changes():
    family = make_family(.1)
    _, owned = setup(family)
    direction = np.zeros((2, 8))
    direction[0, :2] = (.25, 1.)
    direction[1, :2] = (.75, -.5)
    step = 1e-5
    coefficient = np.array(family.coefficients)
    center = owned.evaluate(family)
    plus = owned.evaluate(LocalIncomingFamily(coefficient+step*direction))
    minus = owned.evaluate(LocalIncomingFamily(coefficient-step*direction))
    fd = (plus['action_gradient']-minus['action_gradient'])/(2*step)
    tangent = np.tensordot(direction.ravel(), center['history_jacobian'], axes=(0, 0))
    np.testing.assert_allclose(fd, tangent, atol=3e-8, rtol=0)


def test_zero_history_keeps_edge_tangent_and_requires_evolved_owner():
    family = LocalIncomingFamily(np.zeros((2, 8)))
    raw, owned = setup(family)
    result = owned.evaluate(family)
    np.testing.assert_array_equal(result['source_cutoff_edge_gradient'], 0.)
    assert np.max(abs(result['source_cutoff_edge_history_tangent'])) > 0
    with pytest.raises(TypeError, match='KSHistoryEvaluator'):
        KSCutoffBridgeEvaluator({'matter': 'arbitrary'})
    record = dict(raw.evaluate(family)['family_records'])
    record['duplicate'] = next(iter(record.values()))
    with pytest.raises(ValueError, match='unique'):
        signed_family_angular_square_weights(record)
