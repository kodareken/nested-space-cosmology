"""Rectangular retarded trust-region and refresh policy controls."""
import numpy as np

from recursive_horizons.nsc_ks_retarded_trust_region import (
    RefreshPolicy, acceptance, broyden_update, flatten_constraints,
    levenberg_step, metric_norm, update_radius)


def test_rectangular_flattening_and_trust_step_reduce_linear_residual():
    gradient = np.array([[2., -1.], [1., .5], [-.5, .25]])
    tangent = np.zeros((2, 3, 2))
    tangent[0, :, 0] = (1., .5, -.2)
    tangent[0, :, 1] = (-.3, .4, .2)
    tangent[1, :, 0] = (.2, -.7, .6)
    tangent[1, :, 1] = (.8, .1, -.5)
    residual, jacobian = flatten_constraints(gradient, tangent)
    assert jacobian.shape == (6, 2)
    metric = np.diag([2., 3.])
    proposal = levenberg_step(residual, jacobian, metric, radius=.5)
    assert proposal['metric_norm'] <= .5 * (1 + 1e-12)
    assert proposal['predicted_merit'] < proposal['parent_merit']


def test_broyden_secant_and_acceptance_use_fresh_measured_values():
    jacobian = np.array([[1., 0.], [0., 1.], [1., 1.]])
    step = np.array([.1, -.2])
    measured = np.array([.11, -.18, -.06])
    updated = broyden_update(jacobian, step, measured)
    np.testing.assert_allclose(updated @ step, measured)
    parent = np.array([[1., .5], [.4, .2]])
    candidate = parent * .5
    predicted = np.concatenate((candidate[:, 0], candidate[:, 1])) * .8
    result = acceptance(parent, candidate, predicted,
                        radius_positive=True, preparation_supported=True)
    assert result['accepted']
    assert result['candidate_merit'] == .5
    rejected = acceptance(parent, candidate, predicted,
                          radius_positive=False, preparation_supported=True)
    assert not rejected['accepted']


def test_radius_and_refresh_rules_match_declared_policy():
    assert update_radius(1., .1, True) == .25
    assert update_radius(1., .9, True, maximum=1.5) == 1.5
    policy = RefreshPolicy()
    for _ in range(4):
        policy, refresh = policy.update(accepted=True, prediction_ratio=.8)
        assert not refresh
    policy, refresh = policy.update(accepted=True, prediction_ratio=.8)
    assert refresh and policy == RefreshPolicy()
    policy, refresh = RefreshPolicy().update(accepted=False, prediction_ratio=.1)
    assert not refresh
    policy, refresh = policy.update(accepted=False, prediction_ratio=.1)
    assert refresh


def test_metric_norm_rejects_shape_mismatch():
    try:
        metric_norm(np.ones(2), np.eye(3))
    except ValueError as error:
        assert 'shape' in str(error)
    else:
        raise AssertionError('metric shape mismatch must fail')
