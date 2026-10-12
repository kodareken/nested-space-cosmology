"""Manufactured checks for the local covariance-error insertion."""
import numpy as np
import pytest
from flint import arb, ctx

pytest.importorskip('flint')

from recursive_horizons.nsc_evolved_incoming_constraints import source_column_matter
from recursive_horizons.nsc_evolved_incoming_state import FixedSourcePreparation
from recursive_horizons.nsc_ks_ball_trajectory import restored_upper
from recursive_horizons.nsc_ks_batched_constraints import KSUpstreamBatch
from recursive_horizons.nsc_ks_covariance_insertion import (
    compare_kernel_to_source_column_matter, covariance_insertion_factor,
    energy_diagonal_eta, homogeneous_gram_lower, insert_from_endpoint_slots,
    insert_from_upstream_batch, numerical_action_operators,
    numerical_covariance_insertion, numerical_operator_trace,
    numerical_retarded_kernels, require_energy_diagonal_preparation,
    signed_covariance_factors, source_difference_insertion_bound,
)
from recursive_horizons.nsc_ks_source_envelope import RHO_SIGMA, _fixed_preparation_digest


def _matter(field, axial, covariance, **parameters):
    empty = np.zeros((0, *np.asarray(field).shape), complex)
    return source_column_matter(
        field, axial, covariance, empty, empty, **parameters)['action_gradient']


def _upper(packed):
    with ctx.workprec(160):
        return restored_upper(packed)


def _sharp_fiber():
    gamma = 2.0 ** -10
    root = 2.0 ** -5
    epsilon = 2.0 ** -12
    upstream = np.array([[root, 0.0, 0.0], [0.0, 1.0, 0.0]], complex)
    kernel = np.array([[1.0, 0.0], [1.0, 0.0]], complex)
    field = (kernel @ upstream)[None, :, :]
    axial = np.zeros_like(field)
    delta_C = np.zeros((3, 3), complex)
    delta_C[0, 0] = epsilon / gamma
    covariance = np.eye(3, dtype=complex) * 0.4
    parameters = dict(mass=1.25, angular=0.0, axial_scale=1.0, radius=1.0, multiplicity=2.0)
    return {
        'upstream': upstream, 'field': field, 'axial': axial, 'delta_C': delta_C,
        'covariance': covariance, 'epsilon': epsilon, 'gamma': gamma, 'parameters': parameters,
    }


def _norms(field):
    return float(np.nextafter(np.linalg.norm(field, axis=(1, 2)).max(), np.inf))


def test_kernel_trace_matches_source_column_matter_with_actual_derivative():
    upstream = np.array([[0.8, 0.2, 0.1j], [-0.1j, 0.7, -0.3]], complex)
    kernel = np.stack((
        np.array([[0.4 + 0.1j, -0.2], [0.05, 0.3 - 0.1j]], complex),
        np.array([[-0.1, 0.2j], [0.3, 0.15 + 0.05j]], complex),
    ))
    axial_kernel = np.stack((
        np.array([[0.07, -0.03j], [-0.02, 0.04]], complex),
        np.array([[0.01 + 0.02j, 0.05], [-0.04, -0.02]], complex),
    ))
    energy = 0.6
    assert not np.allclose(axial_kernel, -1j * energy * kernel)
    weight = np.sqrt(0.25 / (2 * np.pi))
    field = weight * np.einsum('zia,as->zis', kernel, upstream, optimize=True)
    axial = weight * np.einsum('zia,as->zis', axial_kernel, upstream, optimize=True)
    substituted = -1j * energy * field
    covariance = np.diag([0.55, 0.42, 0.31]).astype(complex)
    parameters = dict(mass=0.3, angular=-0.8, axial_scale=0.9, radius=1.4, multiplicity=2.0)
    upstream.setflags(write=False)
    field.setflags(write=False)
    axial.setflags(write=False)
    covariance.setflags(write=False)
    frozen = upstream.tobytes(), field.tobytes(), axial.tobytes(), covariance.tobytes()
    check = compare_kernel_to_source_column_matter(upstream, field, axial, covariance, **parameters)
    assert check['discrepancy'] < 1e-9
    assert check['certificate'] is False
    assert check['quadrature_applied_again'] is False
    assert check['evolved_momentum_replaced_by_minus_energy'] is False
    assert (upstream.tobytes(), field.tobytes(), axial.tobytes(), covariance.tobytes()) == frozen
    direct = _matter(field, axial, covariance, **parameters)
    wrong_derivative = _matter(field, substituted, covariance, **parameters)
    assert not np.allclose(direct, wrong_derivative)
    other_radius = dict(parameters, radius=2.0)
    assert not np.allclose(direct, _matter(field, axial, covariance, **other_radius))
    assert compare_kernel_to_source_column_matter(
        upstream, field, axial, covariance, **other_radius)['discrepancy'] < 1e-9
    twice = weight * field
    twice_axial = weight * axial
    assert not np.allclose(direct, _matter(twice, twice_axial, covariance, **parameters))


def test_history_difference_matches_the_inserted_kernel_trace():
    upstream = np.array([[0.9, -0.15, 0.05], [0.1j, 0.6, 0.2]], complex)
    weight = np.sqrt(0.2 / (2 * np.pi))
    reference_kernel = np.stack((np.eye(2), 0.25 * np.eye(2))).astype(complex)
    reference_axial = np.stack((
        0.05 * np.array([[0.0, -1j], [1j, 0.0]]), np.zeros((2, 2))))
    history_kernel = reference_kernel + np.array([[[0.02, -0.01j], [0.03, 0.0]],
                                                  [[-0.01, 0.02], [0.0, 0.04j]]])
    history_axial = reference_axial + np.array([[[0.01, 0.0], [-0.02j, 0.015]],
                                                [[0.0, 0.01], [0.02, -0.01]]])
    assert not np.allclose(history_axial, -1j * 0.4 * history_kernel)
    reference = weight * np.einsum('zia,as->zis', reference_kernel, upstream, optimize=True)
    reference_z = weight * np.einsum('zia,as->zis', reference_axial, upstream, optimize=True)
    history = weight * np.einsum('zia,as->zis', history_kernel, upstream, optimize=True)
    history_z = weight * np.einsum('zia,as->zis', history_axial, upstream, optimize=True)
    delta_Q = np.array([[1.5e-4, -2e-5j], [2e-5j, -1.0e-4]], complex)
    delta_C = numerical_covariance_insertion(upstream, delta_Q)['delta_C']
    covariance = np.diag([0.5, 0.45, 0.35]).astype(complex)
    assert np.linalg.eigvalsh(covariance + delta_C).max() <= 1
    parameters = dict(mass=0.35, angular=0.5, axial_scale=0.8, radius=1.2, multiplicity=1.5)
    base = _matter(history, history_z, covariance, **parameters) - _matter(
        reference, reference_z, covariance, **parameters)
    changed = _matter(history, history_z, covariance + delta_C, **parameters) - _matter(
        reference, reference_z, covariance + delta_C, **parameters)
    actual = changed - base
    history_kernels = numerical_retarded_kernels(upstream, history, history_z)
    reference_kernels = numerical_retarded_kernels(upstream, reference, reference_z)
    history_ops = numerical_action_operators(history_kernels['K'], history_kernels['L'], **parameters)
    reference_ops = numerical_action_operators(
        reference_kernels['K'], reference_kernels['L'], **parameters)
    difference = {
        'H_N': history_ops['H_N'] - reference_ops['H_N'],
        'H_beta': history_ops['H_beta'] - reference_ops['H_beta'],
    }
    reproduced = upstream @ delta_C @ upstream.conj().T
    traced = numerical_operator_trace(reproduced, difference)['action_trace'].real
    assert np.allclose(actual, traced, atol=1e-9, rtol=0)
    assert numerical_operator_trace(reproduced, difference)['certificate'] is False
    assert history_ops['evolved_momentum_replaced_by_minus_energy'] is False
    epsilon = float(np.nextafter(np.linalg.norm(reproduced, 2), np.inf))
    factor = covariance_insertion_factor(upstream, epsilon, bits=120)
    bound = source_difference_insertion_bound(
        _norms(reference), _norms(history - reference), _norms(reference_z),
        _norms(history_z - reference_z), 0.0, 0.0, 0.0, 0.0, factor['eta_upper'],
        mass=0.35, absolute_angular=0.5, axial_lower=0.8, radius_lower=1.2,
        multiplicity=1.5, axial_scale=0.8, radius=1.2, bits=120)
    limits = [_upper(bound['N']), _upper(bound['beta'])]
    assert np.all(np.max(np.abs(actual), axis=0) <= [float(item) for item in limits])
    assert bound['physical_source_coverage'] is False
    assert bound['full_upstream_gate'] == 'OPEN'
    assert bound['quadrature_applied_again'] is False


def test_zero_history_change_has_zero_source_difference_error():
    upstream = np.array([[0.7, 0.1, -0.2], [0.0, 0.5, 0.3j]], complex)
    kernel = np.array([[[0.2, -0.1j], [0.05, 0.4]]], complex)
    axial_kernel = np.array([[[0.03, 0.01], [-0.02, 0.05j]]], complex)
    field = np.einsum('zia,as->zis', kernel, upstream, optimize=True)
    axial = np.einsum('zia,as->zis', axial_kernel, upstream, optimize=True)
    covariance = np.diag([0.6, 0.3, 0.2]).astype(complex)
    delta_C = numerical_covariance_insertion(
        upstream, np.array([[2e-4, 0.0], [0.0, -1e-4]], complex))['delta_C']
    parameters = dict(mass=0.4, angular=-0.3, axial_scale=1.1, radius=1.3, multiplicity=1.0)
    absolute = _matter(field, axial, covariance + delta_C, **parameters) - _matter(
        field, axial, covariance, **parameters)
    difference = (_matter(field, axial, covariance + delta_C, **parameters)
                  - _matter(field, axial, covariance + delta_C, **parameters))
    assert np.max(np.abs(absolute)) > 1e-6
    assert np.max(np.abs(difference)) == 0
    bound = source_difference_insertion_bound(
        _norms(field), 0.0, _norms(axial), 0.0, 0.0, 0.0, 0.0, 0.0, 0.2,
        mass=0.4, absolute_angular=0.3, axial_lower=1.1, radius_lower=1.3,
        multiplicity=1.0, bits=80)
    assert _upper(bound['N']) == 0
    assert _upper(bound['beta']) == 0
    assert bound['history_minus_reference_cancellation'] is True


def test_near_singular_fiber_ignoring_gamma_underestimates():
    sample = _sharp_fiber()
    inserted = numerical_covariance_insertion(
        sample['upstream'], np.diag([sample['epsilon'], 0.0]).astype(complex))
    operator_norm = np.linalg.norm(inserted['delta_C'], 2)
    assert operator_norm > sample['epsilon'] * 100
    factor = covariance_insertion_factor(sample['upstream'], sample['epsilon'], bits=120)
    eta = _upper(factor['eta_upper'])
    gamma = _upper(factor['gamma_lower'])
    assert gamma <= sample['gamma']
    assert gamma > sample['gamma'] * (1 - 1e-8)
    assert eta >= operator_norm
    assert factor['gram_identity_assumed'] is False
    assert factor['column_phases_compared'] is False
    parameters = sample['parameters']
    actual = _matter(
        sample['field'], sample['axial'], sample['covariance'] + sample['delta_C'], **parameters)
    actual -= _matter(sample['field'], sample['axial'], sample['covariance'], **parameters)
    expected = parameters['multiplicity'] * parameters['mass'] * sample['epsilon'] * 2
    assert np.allclose(actual, [[expected, 0.0]], atol=1e-12, rtol=0)
    common = dict(
        reference_norm=0.0, reference_axial_norm=0.0, reference_error=0.0,
        reference_axial_error=0.0, difference_axial_norm=0.0, difference_axial_error=0.0,
        difference_error=0.0, mass=parameters['mass'], absolute_angular=0.0,
        axial_lower=1.0, radius_lower=1.0, multiplicity=parameters['multiplicity'], bits=120)
    enclosed = source_difference_insertion_bound(
        difference_norm=_norms(sample['field']), eta=factor['eta_upper'], **common)
    ignored = source_difference_insertion_bound(
        difference_norm=_norms(sample['field']), eta=sample['epsilon'], **common)
    change = float(np.max(np.abs(actual)))
    assert change <= float(_upper(enclosed['N']))
    assert change > float(_upper(ignored['N'])) * 100
    assert float(_upper(enclosed['N'])) < change * 2


def test_weight_is_applied_once():
    sample = _sharp_fiber()
    weight = np.sqrt(8.0 / (2 * np.pi))
    assert weight > 1
    weighted = weight * sample['field']
    parameters = sample['parameters']
    once = _matter(weighted, sample['axial'], sample['covariance'] + sample['delta_C'], **parameters)
    once -= _matter(weighted, sample['axial'], sample['covariance'], **parameters)
    twice = _matter(
        weight * weighted, sample['axial'], sample['covariance'] + sample['delta_C'], **parameters)
    twice -= _matter(weight * weighted, sample['axial'], sample['covariance'], **parameters)
    bound = source_difference_insertion_bound(
        0.0, _norms(weighted), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        sample['epsilon'] / sample['gamma'], mass=parameters['mass'], absolute_angular=0.0,
        axial_lower=1.0, radius_lower=1.0, multiplicity=parameters['multiplicity'], bits=120)
    # Float contraction can sit a step above the exact value. The certificate
    # has to cover that exact once-weighted contraction, not a second factor.
    limit = _upper(bound['N'])
    once_change = float(np.max(np.abs(once)))
    twice_change = float(np.max(np.abs(twice)))
    assert twice_change > once_change * weight
    with ctx.workprec(180):
        exact_once = arb(5) / arb(1024) / arb.pi()
        exact_twice = exact_once * arb(4) / arb.pi()
        assert exact_once <= limit
        assert exact_twice > limit * arb(weight)
    assert bound['quadrature_applied_again'] is False


def test_field_and_source_errors_jointly_remain_enclosed():
    sample = _sharp_fiber()
    true = 1.5 * sample['field']
    parameters = sample['parameters']
    actual = _matter(true, sample['axial'], sample['covariance'] + sample['delta_C'], **parameters)
    actual -= _matter(true, sample['axial'], sample['covariance'], **parameters)
    captured = _norms(sample['field'])
    missing_enlargement = source_difference_insertion_bound(
        0.0, captured, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, sample['epsilon'] / sample['gamma'],
        mass=parameters['mass'], absolute_angular=0.0, axial_lower=1.0, radius_lower=1.0,
        multiplicity=parameters['multiplicity'], bits=120)
    joint = source_difference_insertion_bound(
        0.0, captured, 0.0, 0.0, 0.0, _norms(true - sample['field']), 0.0, 0.0,
        sample['epsilon'] / sample['gamma'], mass=parameters['mass'], absolute_angular=0.0,
        axial_lower=1.0, radius_lower=1.0, multiplicity=parameters['multiplicity'], bits=120)
    change = float(np.max(np.abs(actual)))
    assert change > float(_upper(missing_enlargement['N']))
    assert change <= float(_upper(joint['N']))
    assert _upper(joint['enlarged_difference_norm']) >= _upper(
        missing_enlargement['enlarged_difference_norm']) * 1.5 * (1 - 1e-12)
    assert joint['field_errors_added_to_norm_uppers'] is True


def test_negative_bound_is_not_copied_and_gram_complement_is_required():
    positive = np.array([[1.0, 0.0, 0.0], [0.0, 2.0, 0.0]], complex)
    negative = np.diag([1.0, -1.0]) @ positive.conj()
    with pytest.raises(ValueError, match='own covariance'):
        signed_covariance_factors(positive, 0.0, negative, None)
    with pytest.raises(ValueError, match='own covariance'):
        signed_covariance_factors(positive, 0.0, negative, 3.0, copy_positive_bound=True)
    factors = signed_covariance_factors(positive, 0.0, negative, 3.0, bits=120)
    assert factors['negative_covariance_bound_copied_from_positive'] is False
    assert factors['physical_source_coverage'] is False
    assert float(_upper(factors['positive']['eta_upper'])) == 0
    assert float(_upper(factors['negative']['eta_upper'])) >= 3
    delta_Q = np.diag([0.0, 3.0]).astype(complex)
    delta_C = numerical_covariance_insertion(negative, delta_Q)['delta_C']
    covariance = np.eye(3, dtype=complex) * 0.2
    assert np.linalg.eigvalsh(covariance + delta_C).max() < 1
    field = negative[None, :, :]
    axial = -1j * field
    parameters = dict(mass=0.0, angular=0.0, axial_scale=1.0, radius=1.0, multiplicity=1.0)
    actual = _matter(field, axial, covariance + delta_C, **parameters) - _matter(
        field, axial, covariance, **parameters)
    assert np.allclose(actual, [[-3.0, -3.0]], atol=1e-10, rtol=0)
    copied = source_difference_insertion_bound(
        0.0, _norms(field), 0.0, _norms(axial), 0.0, 0.0, 0.0, 0.0,
        factors['positive']['eta_upper'], mass=0.0, absolute_angular=0.0,
        axial_lower=1.0, radius_lower=1.0, multiplicity=1.0, bits=120)
    accounted = source_difference_insertion_bound(
        0.0, _norms(field), 0.0, _norms(axial), 0.0, 0.0, 0.0, 0.0,
        factors['negative']['eta_upper'], mass=0.0, absolute_angular=0.0,
        axial_lower=1.0, radius_lower=1.0, multiplicity=1.0, bits=120)
    change = np.max(np.abs(actual), axis=0)
    assert np.any(change > [float(_upper(copied['N'])), float(_upper(copied['beta']))])
    assert np.all(change <= [float(_upper(accounted['N'])), float(_upper(accounted['beta']))])
    assert accounted['gram_identity_assumed'] is False
    assert accounted['upstream_records_loaded'] is False


def test_column_phase_is_not_required_for_the_gram_factor():
    upstream = np.array([[0.6, -0.2j, 0.1], [0.15, 0.5, -0.25j]], complex)
    phased = upstream * np.exp(1j * np.array([0.3, -1.1, 2.0]))
    with ctx.workprec(160):
        original = float(homogeneous_gram_lower(upstream, bits=160))
        rotated = float(homogeneous_gram_lower(phased, bits=160))
    assert original == pytest.approx(rotated, rel=1e-8)
    delta_Q = np.array([[1e-3, 2e-4j], [-2e-4j, -5e-4]], complex)
    inserted = numerical_covariance_insertion(phased, delta_Q)['delta_C']
    factor = covariance_insertion_factor(
        phased, float(np.nextafter(np.linalg.norm(delta_Q, 2), np.inf)), bits=120)
    assert np.linalg.norm(inserted, 2) <= float(_upper(factor['eta_upper'])) * (1 + 1e-9)


def test_absent_incomplete_and_invalid_sign_bounds_are_rejected():
    upstream = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], complex)
    good = dict(
        reference_norm=1.0, difference_norm=0.1, reference_axial_norm=1.0,
        difference_axial_norm=0.1, reference_error=0.0, difference_error=0.0,
        reference_axial_error=0.0, difference_axial_error=0.0, eta=0.01,
        mass=0.2, absolute_angular=0.3, axial_lower=0.8, radius_lower=1.1,
        multiplicity=1.0, bits=80)
    with pytest.raises(ValueError, match='missing'):
        source_difference_insertion_bound(**{**good, 'difference_error': None})
    with pytest.raises(ValueError, match='missing'):
        source_difference_insertion_bound(**{**good, 'reference_axial_error': None})
    with pytest.raises(ValueError, match='missing'):
        covariance_insertion_factor(upstream, None)
    with pytest.raises(ValueError, match='invalid sign'):
        source_difference_insertion_bound(**{**good, 'axial_lower': arb(-0.2).union(arb(0.5))})
    with pytest.raises(ValueError, match='invalid sign'):
        source_difference_insertion_bound(**{**good, 'mass': arb(-0.1).union(arb(0.2))})
    with pytest.raises(ValueError, match='declared slice'):
        source_difference_insertion_bound(**{**good, 'axial_scale': 0.7})
    with pytest.raises(ValueError, match='rank failure'):
        homogeneous_gram_lower(np.array([[1.0, 0.0, 0.0], [0.0, 0.0, 0.0]], complex))
    with pytest.raises(ValueError, match='incomplete|2x3'):
        homogeneous_gram_lower(np.eye(2, dtype=complex))
    labels = np.array([0.4, 0.4, -0.4])
    with pytest.raises(ValueError, match='not collapsed'):
        require_energy_diagonal_preparation(np.zeros((2, 3), complex), labels)
    cross = np.zeros((6, 6), complex)
    cross[0, 3] = 0.2
    with pytest.raises(ValueError, match='energy-diagonal'):
        require_energy_diagonal_preparation(
            np.zeros((2, 6), complex), np.array([0.2, 0.2, 0.2, 0.5, 0.5, 0.5]), cross)
    with pytest.raises(ValueError, match='own covariance'):
        energy_diagonal_eta([
            {'columns': upstream, 'epsilon': 0.01, 'energy': 0.4, 'energy_sign': 1},
            {'columns': upstream, 'epsilon': None, 'energy': -0.4, 'energy_sign': -1},
        ])
    mixed = energy_diagonal_eta([
        {'columns': upstream, 'epsilon': 1e-6, 'energy': 0.4, 'energy_sign': 1},
        {'columns': _sharp_fiber()['upstream'], 'epsilon': _sharp_fiber()['epsilon'],
         'energy': 1.7, 'energy_sign': 1},
    ], bits=120)
    assert float(_upper(mixed['eta_upper'])) >= 0.25 * (1 - 1e-8)
    assert mixed['negative_covariance_bound_copied_from_positive'] is False
    assert mixed['energies_collapsed'] is False


def _slots(**overrides):
    values = dict(
        reference_norm=0.2, reference_axial_norm=0.1, difference_norm=0.05,
        difference_axial_norm=0.04, reference_error=0.0, reference_axial_error=0.0,
        difference_error=0.01, difference_axial_error=0.0, mass=0.2,
        absolute_angular=0.4, axial_lower=0.9, radius_lower=1.2, multiplicity=1.0)
    values.update(overrides)
    return {name: {'value': item, 'owner': 'test slot', 'reason': None} for name, item in values.items()}


def test_endpoint_slots_require_errors_and_do_not_fill_a_budget():
    upstream = np.array([[1.0, 0.2, 0.0], [0.0, 0.8, 0.1]], complex)
    fibers = [{'columns': upstream, 'epsilon': 1e-3, 'energy': 0.3, 'energy_sign': 1}]
    matter = _slots()
    matter['source_norm'] = {'value': 50.0, 'owner': 'not an epsilon', 'reason': None}
    scope = 'manufactured one-fiber history difference; not a physical source'
    result = insert_from_endpoint_slots(matter, fibers, caller_scope=scope, bits=80)
    assert result['caller_scope'] == scope
    assert result['physical_upstream_budget_component'] is None
    assert result['full_physical_budget_populated'] is False
    assert result['physical_source_coverage'] is False
    assert result['source_bound_scope'].startswith('caller-supplied')
    assert result['full_upstream_gate'] == 'OPEN'
    assert result['endpoint_slot_source_norm_used_as_epsilon'] is False
    assert result['partial_source_scope'] is True
    assert float(_upper(result['eta_upper'])) < 1
    untouched = insert_from_endpoint_slots(
        _slots(), fibers, caller_scope=scope, bits=80)
    assert result['N'] == untouched['N']
    assert result['beta'] == untouched['beta']
    broken = _slots()
    broken['difference_error'] = {'value': None, 'owner': None, 'reason': 'absent'}
    with pytest.raises(ValueError, match='missing'):
        insert_from_endpoint_slots(broken, fibers, caller_scope=scope)
    with pytest.raises(ValueError, match='scope'):
        insert_from_endpoint_slots(matter, fibers, caller_scope='  ')
    assert matter['difference_norm']['value'] == 0.05


def _batch(upstream, covariance, energies, *, energy_sign, angular_sign, angular):
    weights = np.full(upstream.shape[1], np.sqrt(0.25 / (2 * np.pi)))
    source = FixedSourcePreparation(covariance, weights, energies)
    rho_up = 1.05
    mass = 0.2
    digest = _fixed_preparation_digest(source, upstream, mass, angular, rho_up, RHO_SIGMA)
    return KSUpstreamBatch(
        'synthetic/rows0-1', 'synthetic', (0, 1), 14, angular_sign, energy_sign,
        mass, angular, rho_up, source, upstream, digest,
        {'energy_folding_factor': 1, 'history_independent': True},
        {'rows': [0, 1]})


def test_upstream_batch_consumes_its_own_epsilon_without_physical_coverage():
    upstream = np.array([[0.9, 0.1, 0.0], [-0.2, 0.7, 0.05]], complex)
    covariance = np.diag([0.5, 0.4, 0.3]).astype(complex)
    batch = _batch(upstream, covariance, np.full(3, 0.5), energy_sign=1, angular_sign=1, angular=0.4)
    norms = dict(
        reference_norm=0.3, difference_norm=0.02, reference_axial_norm=0.2,
        difference_axial_norm=0.01, reference_error=0.0, difference_error=1e-4,
        reference_axial_error=0.0, difference_axial_error=0.0, mass=0.2,
        absolute_angular=0.4, axial_lower=0.8, radius_lower=1.2, multiplicity=1.0)
    before = batch.initial_columns.tobytes(), batch.source.covariance.tobytes()
    result = insert_from_upstream_batch(
        batch, [1e-3], norms, caller_scope='synthetic positive fiber, not the proved row', bits=80)
    assert (batch.initial_columns.tobytes(), batch.source.covariance.tobytes()) == before
    assert result['physical_source_coverage'] is False
    assert result['source_evidence_authenticated'] is False
    assert result['preparation_ode_rerun'] is False
    assert result['archived_columns_modified'] is False
    assert result['physical_upstream_budget_component'] is None
    assert result['full_upstream_gate'] == 'OPEN'
    assert result['energy_sign'] == 1
    negative = _batch(
        np.diag([1.0, -1.0]) @ upstream.conj(), covariance, np.full(3, -0.5),
        energy_sign=-1, angular_sign=-1, angular=-0.4)
    with pytest.raises(ValueError, match='own covariance'):
        insert_from_upstream_batch(
            negative, None, norms, caller_scope='negative fiber without its own epsilon')
    relabeled = KSUpstreamBatch(
        'group14/low16_1/rows0-1', 'group14/low16_1', (0, 1), 14, 1, 1, 0.2, 0.4, 1.05,
        batch.source, upstream, batch.preparation_digest,
        {'energy_folding_factor': 1}, {'rows': [0, 1]})
    labeled = insert_from_upstream_batch(
        relabeled, [1e-3], norms, caller_scope='label matches a known row but is not evidence',
        bits=80)
    assert labeled['source_evidence_authenticated'] is False
    assert labeled['physical_source_coverage'] is False
    assert labeled['upstream_records_loaded'] is False
    for name in ('mass','absolute_angular'):
        with pytest.raises(ValueError,match='channel bounds'):
            insert_from_upstream_batch(batch,[1e-3],{**norms,name:0},caller_scope='mismatched channel')


@pytest.mark.parametrize('denominator',('radius','axial'))
def test_geometry_denominators_use_the_lower_endpoint(denominator):
    from recursive_horizons.nsc_transmitting_dirac_domain import S1,S2,S3
    angle=.1;epsilon=.1
    U=np.cos(angle)*np.eye(2)+1j*np.sin(angle)*S1
    A=np.array([[1,0,0],[0,1,0]],complex)
    reference=A[None];history=(U@A)[None]
    energy=1. if denominator=='axial' else 0.
    reference_z=-1j*energy*reference;history_z=-1j*energy*history
    source=.4*np.eye(3,dtype=complex);delta=np.zeros((3,3),complex)
    delta[:2,:2]=epsilon*(S2 if denominator=='axial' else S3)
    parameters=dict(mass=0.,angular=0. if denominator=='axial' else 1.,
                    axial_scale=1.,radius=1.,multiplicity=1.)
    actual=(_matter(history,history_z,source+delta,**parameters)-_matter(reference,reference_z,source+delta,**parameters)
            -_matter(history,history_z,source,**parameters)+_matter(reference,reference_z,source,**parameters))
    geometry={'axial_lower':1.,'radius_lower':1.}
    geometry['axial_lower' if denominator=='axial' else 'radius_lower']=arb('1.5','.5')
    bound=source_difference_insertion_bound(_norms(reference),_norms(history-reference),
        _norms(reference_z),_norms(history_z-reference_z),0,0,0,0,epsilon,
        mass=0,absolute_angular=parameters['angular'],multiplicity=1,**geometry,bits=160)
    assert max(abs(actual[:,0]))<=float(_upper(bound['N']))


def test_integer_bounds_and_explicit_energy_labels_are_preserved():
    upstream=np.array([[1,0,0],[0,1,0]],complex)
    exact=2**100+1
    factor=covariance_insertion_factor(upstream,exact,bits=160)
    with ctx.workprec(180):assert _upper(factor['eta_upper'])>=arb(exact)
    with pytest.raises(ValueError,match='energy sign'):
        energy_diagonal_eta([{'columns':upstream,'epsilon':1e-4,'energy':1.,'energy_sign':1.5}])
    with pytest.raises(ValueError,match='source energy'):
        energy_diagonal_eta([{'columns':upstream,'epsilon':1e-4,'energy_sign':1}])


def test_relative_source_phase_is_not_erased_with_coherence():
    A=np.array([[1,0,0],[0,1,0]],complex)
    C=np.array([[.5,.2j,0],[-.2j,.5,0],[0,0,.1]],complex)
    moved=A*np.exp(1j*np.array([.3,-.7,0.]))
    assert abs(float(homogeneous_gram_lower(A))-float(homogeneous_gram_lower(moved)))<1e-14
    assert np.linalg.norm(moved@C@moved.conj().T-A@C@A.conj().T,2)>.1
