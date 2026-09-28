"""Signed source identity and fixed reference continuation controls."""
import numpy as np
import pytest

from recursive_horizons.nsc_ks_source_inventory import ReferenceSourcePanel
from recursive_horizons.nsc_ks_difference_envelope import evolve_ks_difference_envelope
from recursive_horizons.nsc_ks_source_envelope import computational_z_grid, physical_incoming_interval, usual_axial_support
from recursive_horizons.nsc_local_incoming_family import LocalIncomingFamily
from recursive_horizons.nsc_paired_horizon_preparation import source_covariance
from recursive_horizons.nsc_retarded_radial_response import ks_generator
from recursive_horizons.nsc_transmitting_dirac_domain import S3


CONFIG = {'surface_gravity': .24, 'omega': 3.9}


def panel():
    energies = np.array([.4, .8, 2.])
    rng = np.random.default_rng(33)
    columns = rng.normal(size=(3, 2, 3)) + 1j * rng.normal(size=(3, 2, 3))
    C = np.array([source_covariance(e, CONFIG['surface_gravity'], CONFIG['omega'], 1.5) for e in energies])
    return ReferenceSourcePanel('synthetic', 14, 1, 1.5, np.sqrt(5), energies, np.array([.1, .2, .3]),
                                C, columns, {'fixture': True})


def test_signed_partner_preserves_the_owned_operator_and_source_complement():
    positive = panel()
    negative = positive.negative_partner(CONFIG)
    assert negative.angular_sign == -1
    np.testing.assert_allclose(negative.covariance, np.eye(3) - positive.covariance.conj(), rtol=0, atol=2e-16)
    for rho in (1., 1.015, 1.03):
        Gp = ks_generator(rho, positive.energies, positive.mass, positive.angular)
        Gn = ks_generator(rho, negative.energies, negative.mass, negative.angular)
        np.testing.assert_array_equal(Gn, S3 @ Gp.conj() @ S3)
    np.testing.assert_array_equal(positive.weights, negative.weights)


def test_batches_neither_drop_coherence_nor_double_weights():
    source = panel()
    pieces = list(source.batches(2))
    for name in ('energies', 'weights', 'covariance', 'amplitudes_at_one'):
        np.testing.assert_array_equal(np.concatenate([getattr(p, name) for p in pieces]), getattr(source, name))
    assert abs(pieces[0].covariance[0, 0, 1]) > 0
    with pytest.raises(ValueError):
        list(source.batches(0))


def test_reference_continuation_and_zero_history_return_the_original_columns():
    original = next(panel().batches(2))
    options = dict(rtol=2e-12, atol=2e-14, max_step=.001)
    source, initial, report = original.fixed_upstream(rho_up=1.03, **options)
    target = np.linspace(*physical_incoming_interval(), 3)
    result = evolve_ks_difference_envelope(source, initial, LocalIncomingFamily(np.zeros((2, 8))),
        computational_z_grid(16), target, original.mass, original.angular, 1.03,
        axial_support=usual_axial_support(), tangents='zero', **options)
    expected = original.amplitudes_at_one.transpose(1, 0, 2).reshape(2, -1)
    np.testing.assert_allclose(result.reference_amplitudes, expected, atol=2e-13, rtol=0)
    assert report['physical_source_error_bound'] is None
    assert report['incoming_state_frozen'] is False
    np.testing.assert_array_equal(source.column_weights, np.repeat(np.sqrt(original.weights / (2 * np.pi)), 3))
