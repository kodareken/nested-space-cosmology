"""Same-slice moments of actual global modes; no seed or history evolution."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from recursive_horizons.nsc_common_ks_trace import restrict_resolved_modes
from recursive_horizons.nsc_incoming_cauchy_jets import reference_incoming_metric_jets
from recursive_horizons.nsc_incoming_state_moments import (
    incoming_adiabatic_reference, incoming_state_moments, incoming_group_factor,
    RADIUS, AXIAL,
)
from recursive_horizons.nsc_spatial_reference_symbol import reference_projector
from recursive_horizons.nsc_transmitting_dirac_domain import S1, S2, S3


ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope='module')
def data():
    record = json.loads((ROOT/'results/development/nsc-pg-massive-mode-resolution.json').read_text())
    path = ROOT/record['payload']['path']
    assert hashlib.sha256(path.read_bytes()).hexdigest() == record['payload']['sha256']
    with np.load(path, allow_pickle=False) as fields:
        yield fields


def test_homogeneous_bloch_equals_existing_proper_KS_reference():
    E = np.array([.2, 2., 9.])
    cases = [(0., np.sqrt(5)), (np.pi/2, 0.), (np.pi/2, -np.sqrt(5))]
    errors = []
    for mass, angular in cases:
        actual, _ = incoming_adiabatic_reference(E, mass, angular)
        independent = reference_projector(*reference_incoming_metric_jets(), -E,
                                          np.full(3, mass), np.full(3, angular))
        expected = independent['orders'].sum(axis=0)
        errors.append(float(np.max(abs(actual-expected))))
    assert max(errors) < 3e-11
    print(json.dumps({'proper_KS_reference_errors': errors}))


@pytest.mark.parametrize('group', [1, 13, 14])
def test_global_mode_moments_and_horizon_coherence(data, group):
    label = data['label']
    candidates = np.flatnonzero((label[:, 0] == group) & (label[:, 1] == 1))
    indices = candidates[[np.argmin(abs(label[candidates, 5]-target)) for target in (.2, 2., 9.)]]
    E = label[indices, 5]
    mass, angular = label[indices[0], 3:5]
    phi = data['interior'][indices, 2]
    source, projector = data['source_covariance'][indices], data['source_projector'][indices]
    before = source.copy()
    result = incoming_state_moments(E, phi, source, projector, mass, angular)
    assert np.array_equal(source, before)
    assert np.isfinite(result['kernels']).all()
    assert max(result['diagnostics'].values()) < 3e-9
    canonical, _ = restrict_resolved_modes(1., E, phi)
    diagonal = np.zeros_like(source)
    diagonal[:, np.arange(3), np.arange(3)] = source[:, np.arange(3), np.arange(3)]
    no_coherence = canonical@diagonal@canonical.swapaxes(-1, -2).conj()
    assert np.linalg.norm(result['covariance'][0]-no_coherence[0]) > 1e-6
    trace_expected = data['transmission'][indices, 0]*(source[:, 2, 2]-source[:, 0, 0])
    trace_error = float(np.max(abs(np.trace(result['difference'], axis1=1, axis2=2)-trace_expected)))
    assert trace_error < 3e-9
    # The momentum channel is covariant T01=-k Tr(C-P)/a, including its sign.
    flux_error = float(np.max(abs(result['kernels'][:, 2]-E/AXIAL*trace_expected.real)))
    assert flux_error < 3e-8
    assert np.array_equal(result['kernels'][:, [0, 1, 3]], result['raw_kernels'][:, [0, 1, 3]])
    assert np.max(abs(result['transmission_norm']-data['transmission'][indices, 0])) < 3e-11
    raw = np.einsum('nij,nji->n', result['difference'],
                    (E/AXIAL)[:, None, None]*np.eye(2)).real
    assert np.max(abs(result['raw_kernels'][:, 2]-raw)) < 3e-15
    np.testing.assert_allclose(result['kernels'][:, 2],
        E/AXIAL*result['transmission_norm']*(source[:, 2, 2]-source[:, 0, 0]).real,
        rtol=4*np.finfo(float).eps, atol=0.)
    print(json.dumps({'group': group, 'diagnostics': result['diagnostics'],
                      'source_trace_identity_error': trace_error,
                      'covariant_flux_error': flux_error,
                      'raw_flux_discrepancy': float(result['per_energy_diagnostics']['raw_flux_discrepancy'].max()),
                      'largest_roundoff_indicator': float(result['roundoff_indicator'].max())}))


def test_signed_energy_fold_and_group_multiplicity(data):
    label = data['label']
    group, energy_target = 14, 2.
    results = {}
    for sign in (-1, 1):
        indices = np.flatnonzero((label[:, 0] == group) & (label[:, 1] == sign))
        index = indices[np.argmin(abs(label[indices, 5]-energy_target))]
        row = label[index]
        results[sign] = incoming_state_moments([row[5]], data['interior'][index:index+1, 2],
            data['source_covariance'][index:index+1], data['source_projector'][index:index+1], row[3], row[4])
    E = results[1]['energies'][0]
    mass, angular = np.pi/2, np.sqrt(5)
    minus_difference = -S3@results[-1]['difference'][0].conj()@S3
    positive_momentum = E/AXIAL
    vertices = np.array([-mass*S1+angular/RADIUS*S2+positive_momentum*S3,
                         positive_momentum*S3, -positive_momentum*np.eye(2), angular/(2*RADIUS)*S2])
    negative_kernels = np.einsum('ij,vji->v', minus_difference, vertices).real
    assert np.max(abs(negative_kernels-results[-1]['kernels'][0])) < 3e-11
    channel = {'angular_eigenvalue': angular, 'copy_count': 2, 'degeneracy': 12}
    factor = incoming_group_factor(channel, (-1, 1))
    unfolded_factor = 2*(12/2)/(4*np.pi*RADIUS**2)/(2*np.pi*AXIAL)
    assert abs(factor-2*unfolded_factor) < 3e-15
    with pytest.raises(ValueError, match='complete actual angular-sign inventory'):
        incoming_group_factor(channel, (1,))


def test_seed_and_massless_LLL_are_not_accepted_as_massive_input():
    with pytest.raises(ValueError, match='massless LLL'):
        incoming_adiabatic_reference([1.], 0., 0.)
    with pytest.raises(ValueError, match='source mode columns'):
        incoming_state_moments([1.], np.eye(2)[None], np.eye(3)[None], np.eye(3)[None], 1., 0.)


def test_coisometry_without_physical_source_sewing_is_rejected(data):
    label = data['label']
    indices = np.flatnonzero((label[:, 0] == 13) & (label[:, 1] == 1) & (label[:, 5] > 2.))
    index = indices[0]
    E = np.array([label[index, 5]])
    # Independent source-space rotation produces a canonical coisometry but
    # destroys the horizon-partner unit norm and column orthogonality.
    wrong = np.array([[1/np.sqrt(2), 1/np.sqrt(2), 0.], [0., 0., 1.]], complex)
    assert np.max(abs(wrong@wrong.conj().T-np.eye(2))) < 3e-15
    probe = np.column_stack((np.eye(2), np.zeros(2)))[None]
    conversion, _ = restrict_resolved_modes(1., E, probe)
    pg = np.linalg.solve(conversion[0, :, :2], wrong)[None]
    with pytest.raises(ValueError, match='physical horizon/incoming sewing'):
        incoming_state_moments(E, pg, data['source_covariance'][index:index+1],
                               data['source_projector'][index:index+1], np.pi/2, 0.)
