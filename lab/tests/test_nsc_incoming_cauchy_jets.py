"""Same-surface normal-jet freedom and a bounded inherited-channel UV screen."""
import numpy as np
import pytest

from recursive_horizons.nsc_incoming_cauchy_jets import (
    IncomingNormalJetChange, incoming_cauchy_jets, incoming_reference_uv_difference,
    NORMAL_ENTRIES,
)
from recursive_horizons.nsc_spatial_reference_symbol import INDEX, INDICES, I


def derivative_probe():
    # Numerical derivative probes only. These values are not candidate
    # physical extrinsic curvature or a selected history/duration.
    return incoming_cauchy_jets((
        IncomingNormalJetChange('a', 1, 0, .02),
        IncomingNormalJetChange('r', 1, 0, -.01),
        IncomingNormalJetChange('a', 1, 1, .003),
        IncomingNormalJetChange('r', 1, 1, .005),
        IncomingNormalJetChange('a', 2, 1, -.002),
        IncomingNormalJetChange('r', 1, 2, .004),
        IncomingNormalJetChange('r', 3, 0, -.003),
        IncomingNormalJetChange('a', 4, 0, .001),
    ))


def test_normal_jets_preserve_same_surface_and_commuting_derivatives():
    domain = derivative_probe(); domain.validate()
    assert len(NORMAL_ENTRIES) == 10
    data = domain.normal_geometry()
    assert abs(data['intrinsic_a']**2-(3*np.pi/2-4)) < 3e-14
    assert abs(data['intrinsic_r']**2-2) < 3e-14
    for base, changed in zip(domain.baseline, domain.fields):
        assert np.max(abs(changed.data[[i for i, (t, z, k) in enumerate(INDICES) if t == 0]]
                          -base.data[[i for i, (t, z, k) in enumerate(INDICES) if t == 0]])) == 0
    for field in domain.fields[2:]:
        assert np.max(abs(field.derivative(t=1).derivative(z=1).data
                          -field.derivative(z=1).derivative(t=1).data)) == 0
    rows = {(r['field'], r['normal_order'], r['spatial_order']): r for r in domain.changed_normal_entries()}
    assert len(rows) == len(domain.changes)
    for change in domain.changes:
        assert abs(rows[change.field, change.normal_order, change.spatial_order]['delta']-change.delta) < 3e-14
    assert abs(data['dz_k_parallel']-.003/data['intrinsic_a']) < 3e-14
    assert abs(data['dz_k_perp']-.005/data['intrinsic_r']) < 3e-14
    with pytest.raises(ValueError, match='normal order'):
        IncomingNormalJetChange('r', 0, 1, .1)
    with pytest.raises(ValueError, match='only normal'):
        IncomingNormalJetChange('N', 1, 0, .1)
    with pytest.raises(ValueError, match='normal order'):
        IncomingNormalJetChange('a', 3, 2, .1)
    # The owner also rejects a hidden direct edit outside its declared slots.
    altered = incoming_cauchy_jets()
    altered.fields[0].data[INDEX[1, 0, 0]] += .01*I
    with pytest.raises(ValueError, match='normal/frame'):
        altered.validate()
    # Matching two shifted copies is not proof of the same physical slice.
    wrong_slice = incoming_cauchy_jets()
    wrong_slice.baseline[2].data[INDEX[0, 0, 0]] += .01*I
    wrong_slice.fields[2].data[INDEX[0, 0, 0]] += .01*I
    with pytest.raises(ValueError, match='existing rho=1'):
        wrong_slice.validate()


def test_inherited_channel_uv_difference_and_exact_lll_separation():
    channels = (
        {'name': 'LLL', 'mass': 0., 'angular': 0.},
        {'name': 'angular', 'mass': 0., 'angular': np.sqrt(5)},
        {'name': 'compact_mass_only', 'mass': np.pi/2, 'angular': 0.},
        {'name': 'compact_mixed', 'mass': np.pi/2, 'angular': np.sqrt(5)},
    )
    row = incoming_reference_uv_difference(derivative_probe(), channels, [16., 32., 64., 128.])
    assert row['initial_hamiltonian_difference'] == 0
    assert row['initial_raw_vertex_difference'] == 0
    assert row['initial_projector_zero_order_difference'] == 0
    assert row['maximum_formal_reference_residual'] < 3e-11
    assert row['maximum_symmetric_vertex_imaginary_part'] < 3e-14
    assert np.max(abs(row['projector_difference_orders'][:, 0])) == 0
    assert np.max(abs(row['reference_vertex_difference_orders'][:, 0])) == 0
    assert np.max(abs(row['projector_difference_orders'][1:, 1:])) > 1e-8
    # This is a measured sample-tail screen, not an infinite-UV certificate.
    envelope = np.max(row['unsigned_vertex_envelope'][1:], axis=-1)
    assert np.all(envelope[:, 1:] < envelope[:, :-1])
    assert np.all(envelope[:, -1] < envelope[:, 0]/20)
    assert row['scope']['constraints_solved'] is False
    assert row['scope']['full_UV_or_Hadamard_certificate'] is False
