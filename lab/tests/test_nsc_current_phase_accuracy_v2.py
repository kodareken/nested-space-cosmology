"""CURRENT-HISTORY phase-value v2: bindings, true-ball distance, no overwrite."""
import json
import os
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np
import pytest
pytest.importorskip('flint')
from flint import arb, ctx

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import derive_nsc_ks_current_phase_accuracy_v2 as P

from recursive_horizons.nsc_dirac_source_phase_bound import pack_upper, unpack_upper
from recursive_horizons.nsc_ks_cutoff_bridge_evaluator import (
    signed_family_angular_square_weights)
from recursive_horizons.nsc_ks_value_evaluator import (
    edge_change, one_direction_metric)
from recursive_horizons.nsc_local_incoming_family import LocalIncomingFamily


def test_write_once_refuses_overwrite_and_foreign_staging(tmp_path):
    path = tmp_path / 'record.json'
    P.write_once(path, {'verdict': 'OPEN'})
    original = path.read_bytes()
    with pytest.raises(FileExistsError):
        P.write_once(path, {'verdict': 'EXISTENCE'})
    assert path.read_bytes() == original
    stage = path.with_name(path.name + f'.{os.getpid()}.tmp')
    stage.write_bytes(b'foreign checkpoint')
    with pytest.raises(FileExistsError):
        P.write_once(path, {'verdict': 'OPEN'})
    assert stage.read_bytes() == b'foreign checkpoint'


def test_true_ball_distance_includes_midpoint_displacement_not_only_radius():
    with ctx.workprec(120):
        ball = arb(1) + arb(0, arb(1) / 10)
        radius_only = ball.rad().upper()
        on_center = P.directed_distance_upper(1.0, ball)
        displaced = P.directed_distance_upper(1.2, ball)
        assert on_center <= radius_only * arb(1) + arb(1) / 10**18
        assert displaced > radius_only
        assert displaced >= arb(3) / 10 - arb(1) / 10**12


def test_pilot_indices_are_first_middle_last_and_all_is_blocked():
    assert P.select_indices(129, 'pilot') == (0, 64, 128)
    assert P.select_indices(257, 'pilot') == (0, 128, 256)
    assert P.select_indices(129, 'all') == tuple(range(129))
    assert P.select_indices(129, 'pilot', (0, 64, 128)) == (0, 64, 128)
    with pytest.raises(ValueError, match='refusing all declared nodes'):
        P.record_run(
            P.DEFAULT_HISTORY, 129, 'all', None, 'blocked-all', 120.0, False)


def test_current_history_identity_and_n32_collocation_on_I():
    path, parent, family, identity = P.load_history()
    assert identity == parent['profile_identity']
    assert identity.startswith('0b0e4ced')
    assert np.asarray(family.coefficients).shape == (2, 32)
    nodes = P.collocation(family, 129)
    assert len(nodes) == 129
    lo, hi = family.interval
    assert float(nodes[0]) == pytest.approx(lo, abs=0)
    assert float(nodes[-1]) == pytest.approx(hi, abs=0)
    assert float(nodes[64]) == family.center
    metric = P.one_direction(family)
    assert len(metric.directions) == 1
    assert len(family.metric().directions) == 64
    with pytest.raises(ValueError, match='129 or 257'):
        P.collocation(family, 47)


def test_history_identity_mismatch_is_rejected(tmp_path, monkeypatch):
    monkeypatch.setattr(P, 'ROOT', tmp_path)
    original = json.loads((ROOT / P.DEFAULT_HISTORY).read_text())
    original['profile_identity'] = '0' * 64
    dest = tmp_path / 'results' / 'development'
    dest.mkdir(parents=True)
    path = dest / 'nsc-ks-gate-history-lm-broyden.json'
    path.write_text(json.dumps(original))
    with pytest.raises(ValueError, match='identity'):
        P.load_history(path)


def test_original_inventory_weights_match_production_ledger():
    from recursive_horizons.nsc_ks_retained_upstream_archive import RetainedUpstreamArchive
    archive = RetainedUpstreamArchive(ROOT)
    weights, records = P.production_angular_weights(archive)
    assert len(records) == 120
    np.testing.assert_allclose(weights, 107952.0, rtol=0, atol=0)
    np.testing.assert_array_equal(
        signed_family_angular_square_weights(records), weights)
    assert float(P.original_archive_rho_up()) == P.ORIGINAL_RHO_UP


def test_true_factor_uses_exact_weight_and_numeric_weights_are_separate():
    with ctx.workprec(120):
        plus, minus = arb(1) / 10**6, arb(2) / 10**6
        a_true = arb(4) / 5
        truth = P.true_edge_gradient(plus, minus, a_true)
        other = P.true_edge_gradient(plus, minus, a_true, weight=107952)
        for left, right in zip(truth, other):
            assert abs(left - right).upper() <= arb(1) / 10**30
        with pytest.raises(ValueError, match='exact original'):
            P.true_edge_gradient(plus, minus, a_true, weight=107953)
        numeric = P.numeric_contraction(
            float(plus.mid()), float(minus.mid()),
            [107952.0, 107952.0], arb(0.844031385900246))
        assert all(value.is_finite() for value in numeric)


def test_checkpoint_binding_mismatch_and_replay_strip_runtime(tmp_path):
    directory = tmp_path / 'artifacts'
    binding = P.evaluation_binding(
        'a' * 64, ROOT / P.DEFAULT_HISTORY, np.array([1.2, 1.21, 1.22]),
        (0, 1, 2), P.ORIGINAL_RHO_UP, np.array([107952.0, 107952.0]),
        0.844031385900246, {'scripts/x.py': '0' * 64})
    P.persist_binding(directory, binding)
    other = dict(binding)
    other['bits'] = 80
    with pytest.raises(ValueError, match='another source/history/numerics'):
        P.persist_binding(directory, other)
    row = {'node': 0, 'z': 1.2, 'cpu_seconds': 3.5, 'complete': True}
    stripped = P.without_runtime(row)
    assert 'cpu_seconds' not in stripped
    assert stripped['node'] == 0
    left = {'coverage': {'mean_cpu_seconds_per_completed_node': 0.7,
                         'estimated_cpu_seconds_for_declared_nodes': 90.0,
                         'requested_nodes': [0, 64, 128]},
            'cpu_seconds': 2.0, 'status': 'OPEN'}
    right = {'coverage': {'mean_cpu_seconds_per_completed_node': 0.9,
                          'estimated_cpu_seconds_for_declared_nodes': 110.0,
                          'requested_nodes': [0, 64, 128]},
             'cpu_seconds': 3.0, 'status': 'OPEN'}
    assert P.canonical(left) == P.canonical(right)


def test_resume_skips_existing_node_and_does_not_rewrite_it(tmp_path, monkeypatch):
    calls = []

    def fake_context(history, node_count):
        family = LocalIncomingFamily(np.zeros((2, 32)))
        nodes = np.linspace(1.1768780671724792, 1.2368780671724793, node_count)
        identity = 'b' * 64
        hashes = {'scripts/derive_nsc_ks_current_phase_accuracy_v2.py': 'c' * 64}
        return {
            'parent_path': ROOT / P.DEFAULT_HISTORY, 'parent': {'profile_identity': identity},
            'family': family, 'identity': identity, 'nodes': nodes,
            'rho_up': P.ORIGINAL_RHO_UP, 'weights': np.array([107952.0, 107952.0]),
            'records': {}, 'a_numeric': 0.844031385900246,
            'a_interval': (0.8440313859002458, 0.8440313859002461),
            'hashes': hashes, 'archive': SimpleNamespace(input_hashes={}),
        }

    def fake_compute(*args, **kwargs):
        calls.append(args[3])
        index = args[3]
        point = float(args[2])
        return {
            'node': int(index), 'z': point, 'z_hex': point.hex(),
            'history_identity': 'b' * 64, 'complete': True, 'obstruction': None,
            'budget_exceeded': False, 'cpu_seconds': 0.01,
            'quadrature': [], 'source_hashes': fake_context(None, 129)['hashes'],
            'rho_up': P.ORIGINAL_RHO_UP,
        }

    monkeypatch.setattr(P, 'load_context', fake_context)
    monkeypatch.setattr(P, 'compute_node', fake_compute)
    monkeypatch.setattr(P, 'paths', lambda prefix: (tmp_path / 'art', tmp_path / 'out.json'))
    ctx = fake_context(None, 129)
    binding = P.evaluation_binding(
        ctx['identity'], ctx['parent_path'], ctx['nodes'], (0, 64, 128),
        ctx['rho_up'], ctx['weights'], ctx['a_numeric'], ctx['hashes'])
    art = tmp_path / 'art'
    P.write_once(P.binding_file(art), binding)
    point = float(ctx['nodes'][0])
    first = {
        'node': 0, 'z': point, 'z_hex': point.hex(),
        'history_identity': ctx['identity'], 'complete': True, 'obstruction': None,
        'budget_exceeded': False, 'cpu_seconds': 0.01, 'quadrature': [],
        'source_hashes': binding['source_hashes'], 'rho_up': P.ORIGINAL_RHO_UP,
    }
    P.write_once(P.node_file(art, 0), first)
    original = P.node_file(art, 0).read_bytes()
    result = P.record_run(
        P.DEFAULT_HISTORY, 129, 'pilot', None, 'resume-unit', 120.0, False)
    assert calls == [64, 128]
    assert P.node_file(art, 0).read_bytes() == original
    assert result['requested_node_indices'] == [0, 64, 128]
    assert result['scope']['between_node_phase_error_bound'] is None
    assert result['scope']['physical_local_gate'] == 'OPEN'
    assert result['scope']['physical_EXISTENCE_certificate'] is False
    assert not str(result['status']).startswith('PASS')
    assert (tmp_path / 'out.json').exists()


def test_edge_change_matches_one_direction_contraction_on_zero_history():
    family = LocalIncomingFamily(np.zeros((2, 32)))
    z = np.array([family.center], float)
    weights = np.array([107952.0, 107952.0])
    a_numeric, _interval = 0.844031385900246, None
    metric = one_direction_metric(family)
    assert len(metric.directions) == 1
    production = edge_change(family, z, weights, a_numeric, P.ORIGINAL_RHO_UP, gauss_nodes=8)
    np.testing.assert_array_equal(production, np.zeros((1, 2)))
    assert P.PHASE_COMPONENT_TOLERANCE == 1e-12
    assert P.GAUSS_COUNTS == (192, 384)


def test_summarize_keeps_between_node_none_and_open_gate():
    binding = {
        'profile_identity': '0b0e4ced' + '0' * 56,
        'history_path': P.DEFAULT_HISTORY,
        'node_count': 129,
        'rho_up': P.ORIGINAL_RHO_UP,
        'signed_angular_square_weights': [107952.0, 107952.0],
        'source_hashes': {},
    }
    with ctx.workprec(120):
        error = [pack_upper(arb(1) / 10**13), pack_upper(arb(1) / 10**13)]
    rows = [{
        'node': 0, 'complete': True, 'obstruction': None, 'cpu_seconds': 1.5,
        'quadrature': [
            {'gauss_nodes_per_normal_panel': 192, 'value_error_upper_N_beta': error,
             'phase_value_control_pass': True},
            {'gauss_nodes_per_normal_panel': 384, 'value_error_upper_N_beta': error,
             'phase_value_control_pass': True},
        ],
    }]
    result = P.summarize(binding, rows, (0,), 129, 1.5, 120.0, True)
    assert result['scope']['between_node_phase_error_bound'] is None
    assert result['coverage']['all_declared_nodes_covered'] is False
    assert result['coverage']['estimated_cpu_seconds_for_declared_nodes'] == pytest.approx(1.5 * 129)
    assert result['status'].startswith('OPEN')
    assert result['quadrature_attempts']['192']['phase_value_control_pass'] is True
    assert unpack_upper(result['quadrature_attempts']['192']['value_error_upper_N_beta'][0]) <= arb(1) / 10**12
