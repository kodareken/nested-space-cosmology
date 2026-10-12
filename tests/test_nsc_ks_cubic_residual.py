"""Cubic leading term added into the saved finite-cutoff value residual."""
from pathlib import Path

import numpy as np
import pytest
from flint import arb

from recursive_horizons.nsc_ks_cubic_inventory import channel_moments, load_inventory
from recursive_horizons.nsc_ks_cubic_residual import (
    VALUE_LABEL, apply_cubic_leading, certify, load_saved_finite_cutoff_value,
    require_moments)
from recursive_horizons.nsc_ks_finite_history_uv_coefficients import ARCHIVE_RHO_UP
from recursive_horizons.nsc_ks_massive_cubic_uv import (
    CURRENT_HISTORY_IDENTITY, require_original_history)
from recursive_horizons.nsc_ks_value_evaluator import target_nodes

ROOT = Path(__file__).resolve().parents[1]
CHANNELS = ('N2', 'N4', 'NM', 'J2', 'J4', 'JM')


@pytest.fixture(scope='module')
def saved():
    return load_saved_finite_cutoff_value(ROOT)


def _bases(nodes, identity, *, j_sign=1, nm_shift=0, omit=None, scalar=False,
           repeat_center=False, surface=False):
    nodes = np.asarray(nodes, dtype=float)
    count = nodes.size
    center = arb(float(nodes[count // 2]))
    bases = {}
    for sign in (1, -1):
        channels = {}
        for offset, name in enumerate(CHANNELS):
            if scalar:
                channels[name] = arb(offset + 1)
                continue
            values = []
            for index in range(count):
                value = arb(index + 1) * arb(offset + 1) + arb(sign)
                if name in ('J2', 'J4', 'JM'):
                    value *= arb(j_sign)
                if name == 'NM':
                    value += arb(nm_shift)
                values.append(value)
            channels[name] = values
        if omit is not None:
            channels.pop(omit, None)
        if surface:
            channels['surface_mass_correction'] = [arb(1) for _ in range(count)]
        domains = [center for _ in range(count)] if repeat_center else [
            arb(float(node)) for node in nodes]
        if scalar:
            domains = center
        bases[sign] = {
            'sign': sign, 'massless_A': True, 'profile_identity': identity,
            'C_M': None, 'rho_up': arb(float(ARCHIVE_RHO_UP)),
            'z_domain': domains, 'channels': channels,
        }
    return bases


def test_saved_replay_binds_history_source_and_nodes(saved):
    family, identity = require_original_history(ROOT)
    assert identity == CURRENT_HISTORY_IDENTITY == saved['profile_identity']
    assert identity.startswith('0b0e4ced')
    assert saved['label'] == VALUE_LABEL
    assert saved['node_count'] == 129
    assert saved['family_count'] == 60 and saved['signed_family_count'] == 120
    assert saved['merit'] == saved['register_merit']
    assert np.array_equal(saved['nodes'], target_nodes(family, 129))
    assert float(saved['nodes'][saved['center_index']]) == float(family.center)
    assert saved['rho_up'] == float(ARCHIVE_RHO_UP)
    assert list(saved['weights']) == [107952.0, 107952.0]
    rows = load_inventory(ROOT)
    assert {int(row['group']) for row in rows if float(row['cutoff']) == 320.0} == {10, 11, 12, 31, 32}
    moments = require_moments(channel_moments(rows))
    assert moments['energy_folding_factor'] == 1 and moments['C_M'] is None


def test_real_addition_changes_every_saved_node(saved):
    bases = _bases(saved['nodes'], saved['profile_identity'])
    report = apply_cubic_leading(
        saved['residual'], bases, ROOT, saved['nodes'],
        coefficient_source='test-nodal-basis')
    assert report['value_residual_array_modified'] is True
    assert report['cubic_leading_applications'] == 1
    assert report['used_raw_coefficient'] is False
    assert report['surface_mass_added_again'] is False
    assert report['center_was_broadcast'] is False
    assert report['all_exact_target_nodes'] is True
    assert report['approximate_is_enclosure'] is False
    assert report['higher_uv_remainder_C_M'] is None
    assert report['physical_local_gate'] == 'OPEN'
    assert np.array_equal(report['residual'], saved['residual'] + report['approximate_leading'])
    assert not np.array_equal(report['residual'], saved['residual'])
    assert not np.allclose(report['approximate_leading'][0], report['approximate_leading'][-1])
    assert not np.array_equal(report['approximate_leading'], np.broadcast_to(
        report['approximate_leading'][saved['center_index']], report['approximate_leading'].shape))
    twice = apply_cubic_leading(
        report['residual'], bases, ROOT, saved['nodes'])
    assert not np.array_equal(twice['residual'], report['residual'])
    with pytest.raises(ValueError, match='higher remainder'):
        certify(report)


def test_omission_sign_double_count_and_node_mismatch(saved):
    identity = saved['profile_identity']
    nodes = saved['nodes']
    honest = _bases(nodes, identity)
    report = apply_cubic_leading(saved['residual'], honest, ROOT, nodes)
    with pytest.raises(ValueError, match='omitted cubic channel'):
        apply_cubic_leading(saved['residual'], _bases(nodes, identity, omit='N4'), ROOT, nodes)
    with pytest.raises(ValueError, match='both characteristic'):
        apply_cubic_leading(saved['residual'], {1: honest[1]}, ROOT, nodes)
    flipped = apply_cubic_leading(
        saved['residual'], _bases(nodes, identity, j_sign=-1), ROOT, nodes)
    for honest_beta, flipped_beta in zip(report['leading_beta'], flipped['leading_beta']):
        assert (honest_beta + flipped_beta).contains(0)
        assert not (honest_beta - flipped_beta).contains(0)
    for honest_n, flipped_n in zip(report['leading_N'], flipped['leading_N']):
        assert (honest_n - flipped_n).contains(0)
    shifted = apply_cubic_leading(
        saved['residual'], _bases(nodes, identity, nm_shift=1), ROOT, nodes)
    assert not np.array_equal(shifted['approximate_leading'], report['approximate_leading'])
    with pytest.raises(ValueError, match='surface mass'):
        apply_cubic_leading(saved['residual'], _bases(nodes, identity, surface=True), ROOT, nodes)
    with pytest.raises(ValueError, match='node-profile mismatch'):
        apply_cubic_leading(saved['residual'][:-1], _bases(nodes[:-1], identity), ROOT, nodes)
    reversed_nodes = nodes[::-1]
    with pytest.raises(ValueError, match='node-profile mismatch'):
        apply_cubic_leading(
            saved['residual'][::-1], _bases(reversed_nodes, identity), ROOT, nodes)
    with pytest.raises(ValueError, match='center scalar must not be broadcast'):
        apply_cubic_leading(saved['residual'], _bases(nodes, identity, scalar=True), ROOT, nodes)
    with pytest.raises(ValueError, match='center scalar must not be broadcast'):
        apply_cubic_leading(
            saved['residual'], _bases(nodes, identity, repeat_center=True), ROOT, nodes)
    with pytest.raises(ValueError, match='center scalar must not be broadcast'):
        apply_cubic_leading(saved['residual'][0], honest, ROOT, nodes)


def test_partial_profile_is_not_filled_and_moments_reject_folding(saved):
    identity = saved['profile_identity']
    subset = [saved['center_index'], 0, saved['node_count'] - 1]
    bases = _bases(saved['nodes'][subset], identity)
    report = apply_cubic_leading(
        saved['residual'][subset], bases, ROOT, saved['nodes'], node_indices=subset)
    assert report['all_exact_target_nodes'] is False
    assert report['completed_node_count'] == 3
    assert report['node_indices'] == subset
    assert 'residual_in_node_order' not in report
    with pytest.raises(ValueError, match='higher remainder'):
        certify(report)
    moments = channel_moments(load_inventory(ROOT))
    bad = dict(moments)
    bad['energy_folding_factor'] = 2
    with pytest.raises(ValueError, match='energy folding'):
        require_moments(bad)
    bad = dict(moments)
    bad['angular_families'] = 59
    with pytest.raises(ValueError, match='omitted angular family'):
        require_moments(bad)


def test_supplied_remainder_is_not_a_constructed_certificate(saved):
    bases = _bases(saved['nodes'][:1], saved['profile_identity'])
    bases[1]['C_M'] = arb(0)
    with pytest.raises(ValueError, match='higher remainder was not constructed'):
        apply_cubic_leading(
            saved['residual'][:1], bases, ROOT, saved['nodes'], node_indices=[0])
