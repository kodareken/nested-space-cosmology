"""CURRENT-HISTORY phase-value v3: tighter radius, descriptors, failed-node count."""
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
import derive_nsc_ks_current_phase_accuracy_v2 as V2
import derive_nsc_ks_current_phase_accuracy_v3 as P

from recursive_horizons.nsc_dirac_source_phase_bound import (
    directed_source_phase_z_enclosure, pack_upper, unit_angular_target,
    unpack_upper)
from recursive_horizons.nsc_ks_value_evaluator import one_direction_metric
from recursive_horizons.nsc_local_incoming_family import LocalIncomingFamily


def test_v2_globals_are_not_patched_by_the_private_enclosure_namespace():
    assert V2.SCHEMA == 'NSC-KS-CURRENT-PHASE-ACCURACY-v2'
    assert P.SCHEMA == 'NSC-KS-CURRENT-PHASE-ACCURACY-v3'
    assert V2.compute_node.__globals__['directed_source_phase_z_enclosure'] is (
        directed_source_phase_z_enclosure)
    assert P._compute_node_core.__globals__ is not V2.compute_node.__globals__
    wrapper = P._compute_node_core.__globals__['directed_source_phase_z_enclosure']
    assert wrapper is P._TIGHTER_ENCLOSURE
    assert wrapper is not directed_source_phase_z_enclosure
    assert wrapper._inner is directed_source_phase_z_enclosure
    assert P.TARGET_RADIUS_DIVISOR == 16


def test_declared_target_is_unit_angular_over_sixteen():
    with ctx.workprec(120):
        expected = unit_angular_target() / arb(16)
        got = P.declared_target_radius()
        assert abs(got - expected).upper() <= arb(1) / 10**40
        packed = P.packed_declared_target()
        assert unpack_upper(packed) >= expected
        assert unpack_upper(packed) <= expected * arb(1) + arb(1) / 10**30


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


def test_true_gradient_uses_exact_107952_and_a_and_distance_includes_midpoint():
    with ctx.workprec(120):
        plus = arb(1) / 10**6 + arb(0, arb(1) / 10**16)
        minus = arb(2) / 10**6 + arb(0, arb(2) / 10**16)
        a_true = arb(4) / 5
        truth = P.true_edge_gradient(plus, minus, a_true)
        same = P.true_edge_gradient(plus, minus, a_true, weight=107952)
        weight = arb(107952)
        two_pi = 2 * arb.pi()
        plus_term, minus_term = weight * plus, weight * minus
        expected = (
            (plus_term - minus_term) / (a_true * two_pi),
            -(plus_term + minus_term) / two_pi,
        )
        for left, right, extra in zip(truth, expected, same):
            assert left.overlaps(right)
            assert extra.overlaps(left)
            assert abs(left.mid() - right.mid()).upper() <= arb(1) / 10**30
        with pytest.raises(ValueError, match='exact original'):
            P.true_edge_gradient(plus, minus, a_true, weight=107953)
        radius = truth[0].rad().upper()
        mid = float(truth[0].mid())
        on_center = P.directed_distance_upper(mid, truth[0])
        displaced = P.directed_distance_upper(mid + 1.0e-10, truth[0])
        assert on_center <= radius * arb(1) + arb(1) / 10**16
        assert displaced > radius
        production = (float(truth[0].mid()), float(truth[1].mid()))
        report = P.contraction_report(
            (float(plus.mid()), float(minus.mid())),
            (107952.0, 107952.0), float(a_true.mid()),
            (float(a_true.mid()), float(a_true.mid())), a_true, truth, production)
        assert report['value_error_is_direct_distance_to_true_ball'] is True
        assert P.GAUSS_COUNTS == (192, 384)
        assert V2.GAUSS_COUNTS == (192, 384)


def test_pilot_indices_are_first_middle_last_and_all_is_blocked():
    assert P.select_indices(129, 'pilot') == (0, 64, 128)
    assert P.select_indices(257, 'pilot') == (0, 128, 256)
    assert P.select_indices(129, 'all') == tuple(range(129))
    with pytest.raises(ValueError, match='refusing all declared nodes'):
        P.record_run(
            P.DEFAULT_HISTORY, 129, 'all', None, 'blocked-all', 60.0, False)


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
    metric = P.one_direction(family)
    assert len(metric.directions) == 1
    assert len(family.metric().directions) == 64
    assert len(one_direction_metric(family).directions) == 1


def test_binding_includes_tighter_radius_and_source_code_hashes():
    hashes = P.source_hashes()
    assert 'scripts/derive_nsc_ks_current_phase_accuracy_v3.py' in hashes
    assert 'docs/nsc-ks-current-phase-accuracy-v3.md' in hashes
    assert 'scripts/derive_nsc_ks_current_phase_accuracy_v2.py' in hashes
    binding = P.evaluation_binding(
        'a' * 64, ROOT / P.DEFAULT_HISTORY, np.array([1.2, 1.21, 1.22]),
        (0, 1, 2), P.ORIGINAL_RHO_UP, np.array([107952.0, 107952.0]),
        0.844031385900246, hashes)
    assert binding['schema'] == P.SCHEMA
    assert binding['parent_schema'] == V2.SCHEMA
    assert binding['target_radius_divisor'] == 16
    assert binding['source_hashes'] == hashes
    with ctx.workprec(120):
        assert unpack_upper(binding['target_radius_upper']) <= (
            unit_angular_target() / arb(16)) * arb(1) + arb(1) / 10**30
        assert unpack_upper(binding['inherited_unit_angular_target_upper']) >= (
            unit_angular_target())


def test_status_names_actual_failed_node_count_not_sample_size():
    hashes = {'scripts/x.py': '0' * 64}
    binding = P.evaluation_binding(
        '0b0e4ced' + '0' * 56, ROOT / P.DEFAULT_HISTORY,
        np.array([1.2, 1.21, 1.22]), (0, 1, 2), P.ORIGINAL_RHO_UP,
        np.array([107952.0, 107952.0]), 0.844031385900246, hashes)
    with ctx.workprec(120):
        miss = [pack_upper(arb(2) / 10**12), pack_upper(arb(1) / 10**13)]
        hit = [pack_upper(arb(1) / 10**13), pack_upper(arb(1) / 10**13)]
    def row(index, errors, passed):
        return {
            'node': index, 'complete': True, 'obstruction': None, 'cpu_seconds': 0.4,
            'quadrature': [
                {'gauss_nodes_per_normal_panel': 192, 'value_error_upper_N_beta': errors,
                 'phase_value_control_pass': passed},
                {'gauss_nodes_per_normal_panel': 384, 'value_error_upper_N_beta': errors,
                 'phase_value_control_pass': passed},
            ],
        }
    rows = [row(0, miss, False), row(1, hit, True), row(2, miss, False)]
    result = P.summarize(binding, rows, (0, 1, 2), 129, 1.2, 60.0, True)
    assert result['quadrature_attempts']['192']['failed_node_count'] == 2
    assert result['quadrature_attempts']['192']['failed_node_indices'] == [0, 2]
    assert 'misses 1e-12 on 2 of 3 requested nodes (failed nodes: 0, 2)' in result['status']
    assert 'misses 1e-12 on 3 nodes' not in result['status']
    assert result['scope']['between_node_phase_error_bound'] is None
    assert result['scope']['physical_local_gate'] == 'OPEN'
    assert result['coverage']['all_declared_nodes_covered'] is False
    assert result['scope']['all_declared_nodes_covered'] is False
    assert 'rows' not in result


def test_all_declared_nodes_covered_requires_every_index():
    hashes = {}
    binding = P.evaluation_binding(
        'b' * 64, ROOT / P.DEFAULT_HISTORY, np.linspace(1.2, 1.22, 3),
        (0, 1, 2), P.ORIGINAL_RHO_UP, np.array([107952.0, 107952.0]),
        0.844031385900246, hashes)
    with ctx.workprec(120):
        hit = [pack_upper(arb(1) / 10**13), pack_upper(arb(1) / 10**13)]
    rows = [{
        'node': index, 'complete': True, 'obstruction': None, 'cpu_seconds': 0.1,
        'quadrature': [
            {'gauss_nodes_per_normal_panel': 192, 'value_error_upper_N_beta': hit,
             'phase_value_control_pass': True},
            {'gauss_nodes_per_normal_panel': 384, 'value_error_upper_N_beta': hit,
             'phase_value_control_pass': True},
        ],
    } for index in range(129)]
    partial = P.summarize(binding, rows[:3], (0, 64, 128), 129, 0.3, 60.0, True)
    assert partial['coverage']['all_declared_nodes_covered'] is False
    full = P.summarize(
        binding, rows, tuple(range(129)), 129, 12.9, 300.0, True)
    assert full['coverage']['all_declared_nodes_covered'] is True
    assert full['scope']['all_declared_nodes_covered'] is True
    assert P.nodal_component_pass(full) is True


def test_checkpoint_rejects_radius_schema_and_hash_changes(tmp_path):
    hashes = {'scripts/derive_nsc_ks_current_phase_accuracy_v3.py': 'c' * 64}
    nodes = np.array([1.2, 1.21, 1.22])
    binding = P.evaluation_binding(
        'b' * 64, ROOT / P.DEFAULT_HISTORY, nodes, (0, 1, 2),
        P.ORIGINAL_RHO_UP, np.array([107952.0, 107952.0]),
        0.844031385900246, hashes)
    directory = tmp_path / 'art'
    P.persist_binding(directory, binding)
    point = float(nodes[0])
    with ctx.workprec(120):
        hit = [pack_upper(arb(1) / 10**13), pack_upper(arb(1) / 10**13)]
    row = {
        'schema': P.SCHEMA, 'node': 0, 'z': point, 'z_hex': point.hex(),
        'history_identity': 'b' * 64, 'complete': True, 'obstruction': None,
        'budget_exceeded': False, 'cpu_seconds': 0.01,
        'bits': P.BITS, 'phase_component_tolerance': P.PHASE_COMPONENT_TOLERANCE,
        'gauss_nodes_per_normal_panel': list(P.GAUSS_COUNTS),
        'exact_angular_square_weight': 107952, 'one_direction_metric': True,
        'target_radius_divisor': 16,
        'target_radius_upper': binding['target_radius_upper'],
        'signed_angular_square_weights': [107952.0, 107952.0],
        'production_axial_scale': 0.844031385900246,
        'rho_up': P.ORIGINAL_RHO_UP, 'source_hashes': hashes,
        'edge_gradient_balls_N_beta': [{'binary64_mid': 0.0}],
        'scope': {'value_error_is_direct_distance_to_true_ball': True},
        'quadrature': [
            {'gauss_nodes_per_normal_panel': 192, 'value_error_upper_N_beta': hit,
             'phase_value_control_pass': True},
            {'gauss_nodes_per_normal_panel': 384, 'value_error_upper_N_beta': hit,
             'phase_value_control_pass': True},
        ],
    }
    P.write_once(P.node_file(directory, 0), row)
    loaded = P.load_node_checkpoint(directory, 0, binding)
    assert loaded['node'] == 0
    mutated = dict(row)
    mutated['target_radius_divisor'] = 1
    P.node_file(directory, 0).write_text(json.dumps(mutated, sort_keys=True) + '\n')
    with pytest.raises(ValueError, match='target_radius_divisor'):
        P.load_node_checkpoint(directory, 0, binding)
    mutated = dict(row)
    mutated['schema'] = V2.SCHEMA
    P.node_file(directory, 0).write_text(json.dumps(mutated, sort_keys=True) + '\n')
    with pytest.raises(ValueError, match='schema'):
        P.load_node_checkpoint(directory, 0, binding)
    mutated = dict(row)
    mutated['source_hashes'] = {'scripts/x.py': 'd' * 64}
    P.node_file(directory, 0).write_text(json.dumps(mutated, sort_keys=True) + '\n')
    with pytest.raises(ValueError, match='code or ledger hashes'):
        P.load_node_checkpoint(directory, 0, binding)


def test_recursive_payload_descriptor_rejects_mutated_node_and_hash(tmp_path):
    directory = tmp_path / 'results' / 'development' / 'artifacts' / 'tree'
    directory.mkdir(parents=True)
    P.write_once(P.binding_file(directory), {'schema': P.SCHEMA})
    P.write_once(P.node_file(directory, 0), {'node': 0})
    P.write_once(P.node_file(directory, 1), {'node': 1})
    payloads = P.artifact_descriptor(directory, (0, 1))
    assert payloads['schema'] == P.PAYLOAD_SCHEMA
    assert payloads['sha256'] == payloads['listing_sha256']
    assert payloads['binding']['path'].endswith('binding.json')
    P.authenticate_payloads(payloads, directory, (0, 1))
    original = P.node_file(directory, 0).read_bytes()
    P.node_file(directory, 0).write_bytes(original + b' ')
    with pytest.raises(ValueError, match='payload hash descriptor changed'):
        P.authenticate_payloads(payloads, directory, (0, 1))
    P.node_file(directory, 0).write_bytes(original)
    broken = json.loads(json.dumps(payloads))
    broken['nodes'][0]['sha256'] = '0' * 64
    with pytest.raises(ValueError, match='payload hash descriptor changed'):
        P.authenticate_payloads(broken, directory, (0, 1))
    broken = json.loads(json.dumps(payloads))
    broken['sha256'] = '0' * 64
    with pytest.raises(ValueError, match='recursive payload hash descriptor changed'):
        P.authenticate_payloads(broken, directory, (0, 1))


def _fake_context(node_count=129):
    family = LocalIncomingFamily(np.zeros((2, 32)))
    nodes = np.linspace(1.1768780671724792, 1.2368780671724793, node_count)
    identity = 'b' * 64
    hashes = P.source_hashes({'results/development/nsc-ks-gate-history-lm-broyden.json': 'e' * 64})
    hashes['scripts/derive_nsc_ks_current_phase_accuracy_v3.py'] = 'c' * 64
    return {
        'parent_path': ROOT / P.DEFAULT_HISTORY, 'parent': {'profile_identity': identity},
        'family': family, 'identity': identity, 'nodes': nodes,
        'rho_up': P.ORIGINAL_RHO_UP, 'weights': np.array([107952.0, 107952.0]),
        'records': {}, 'a_numeric': 0.844031385900246,
        'a_interval': (0.8440313859002458, 0.8440313859002461),
        'hashes': hashes, 'archive': SimpleNamespace(input_hashes={}),
    }


def _fake_row(ctx, index, hashes):
    point = float(ctx['nodes'][index])
    with ctx_work():
        hit = [pack_upper(arb(1) / 10**13), pack_upper(arb(1) / 10**13)]
    return {
        'schema': P.SCHEMA, 'node': int(index), 'z': point, 'z_hex': point.hex(),
        'history_identity': ctx['identity'], 'complete': True, 'obstruction': None,
        'budget_exceeded': False, 'cpu_seconds': 0.01,
        'bits': P.BITS, 'phase_component_tolerance': P.PHASE_COMPONENT_TOLERANCE,
        'gauss_nodes_per_normal_panel': list(P.GAUSS_COUNTS),
        'exact_angular_square_weight': 107952, 'one_direction_metric': True,
        'target_radius_divisor': 16,
        'target_radius_upper': P.packed_declared_target(),
        'signed_angular_square_weights': [107952.0, 107952.0],
        'production_axial_scale': ctx['a_numeric'],
        'rho_up': P.ORIGINAL_RHO_UP, 'source_hashes': hashes,
        'edge_gradient_balls_N_beta': [{'binary64_mid': 0.0}],
        'scope': {
            'between_node_phase_error_bound': None,
            'value_error_is_direct_distance_to_true_ball': True,
            'physical_local_gate': 'OPEN',
        },
        'quadrature': [
            {'gauss_nodes_per_normal_panel': 192, 'value_error_upper_N_beta': hit,
             'phase_value_control_pass': True},
            {'gauss_nodes_per_normal_panel': 384, 'value_error_upper_N_beta': hit,
             'phase_value_control_pass': True},
        ],
    }


def ctx_work():
    return ctx.workprec(120)


def test_resume_skips_existing_node_and_does_not_rewrite_it(tmp_path, monkeypatch):
    calls = []
    ctx_data = _fake_context()

    def fake_context(history, node_count):
        return ctx_data

    def fake_compute(*args, **kwargs):
        calls.append(args[3])
        return _fake_row(ctx_data, args[3], ctx_data['hashes'])

    monkeypatch.setattr(P, 'load_context', fake_context)
    monkeypatch.setattr(P, 'compute_node', fake_compute)
    monkeypatch.setattr(P, 'paths', lambda prefix: (tmp_path / 'art', tmp_path / 'out.json'))
    binding = P.evaluation_binding(
        ctx_data['identity'], ctx_data['parent_path'], ctx_data['nodes'], (0, 64, 128),
        ctx_data['rho_up'], ctx_data['weights'], ctx_data['a_numeric'], ctx_data['hashes'])
    art = tmp_path / 'art'
    first = _fake_row(ctx_data, 0, binding['source_hashes'])
    P.write_once(P.binding_file(art), binding)
    P.write_once(P.node_file(art, 0), first)
    original = P.node_file(art, 0).read_bytes()
    result = P.record_run(
        P.DEFAULT_HISTORY, 129, 'pilot', None, 'resume-unit', 60.0, False)
    assert calls == [64, 128]
    assert P.node_file(art, 0).read_bytes() == original
    assert result['requested_node_indices'] == [0, 64, 128]
    assert result['scope']['between_node_phase_error_bound'] is None
    assert result['scope']['physical_local_gate'] == 'OPEN'
    assert result['payloads']['nodes'][0]['node'] == 0
    assert result['hash_descriptor'] == result['payloads']['sha256']
    assert 'rows' not in result
    assert (tmp_path / 'out.json').exists()


def test_check_replays_and_rejects_mutated_node_payload_and_hash(tmp_path, monkeypatch):
    ctx_data = _fake_context()

    def fake_context(history, node_count):
        return ctx_data

    def fake_compute(*args, **kwargs):
        return _fake_row(ctx_data, args[3], ctx_data['hashes'])

    monkeypatch.setattr(P, 'load_context', fake_context)
    monkeypatch.setattr(P, 'compute_node', fake_compute)
    monkeypatch.setattr(P, 'paths', lambda prefix: (tmp_path / 'art', tmp_path / 'out.json'))
    recorded = P.record_run(
        P.DEFAULT_HISTORY, 129, 'pilot', None, 'check-unit', 60.0, False)
    replayed = P.check_run(
        P.DEFAULT_HISTORY, 129, 'pilot', None, 'check-unit', 60.0, False)
    assert P.canonical(recorded) == P.canonical(replayed)
    assert replayed['scope']['physical_local_gate'] == 'OPEN'

    node_path = P.node_file(tmp_path / 'art', 0)
    original = node_path.read_bytes()
    node_path.write_bytes(original + b' ')
    with pytest.raises(ValueError, match='payload hash descriptor changed'):
        P.check_run(P.DEFAULT_HISTORY, 129, 'pilot', None, 'check-unit', 60.0, False)
    node_path.write_bytes(original)

    saved_path = tmp_path / 'out.json'
    saved = json.loads(saved_path.read_text())
    saved['payloads']['nodes'][1]['sha256'] = '0' * 64
    saved_path.write_text(json.dumps(saved, sort_keys=True, indent=2) + '\n')
    with pytest.raises(ValueError, match='payload hash descriptor changed'):
        P.check_run(P.DEFAULT_HISTORY, 129, 'pilot', None, 'check-unit', 60.0, False)
    saved_path.write_text(json.dumps(recorded, sort_keys=True, indent=2) + '\n')

    mutated = json.loads(original.decode())
    mutated['node'] = 7
    node_path.write_text(json.dumps(mutated, sort_keys=True) + '\n')
    with pytest.raises(ValueError):
        P.check_run(P.DEFAULT_HISTORY, 129, 'pilot', None, 'check-unit', 60.0, False)


def test_recorded_result_keeps_gate_open_and_between_node_none():
    for prefix in ('pilot', 'n129'):
        _directory, output = P.paths(prefix)
        if not output.exists():
            continue
        saved = json.loads(output.read_text())
        assert saved['schema'] == P.SCHEMA
        assert saved['history_identity'].startswith('0b0e4ced')
        assert saved['scope']['physical_local_gate'] == 'OPEN'
        assert saved['scope']['between_node_phase_error_bound'] is None
        assert saved['scope']['physical_EXISTENCE_certificate'] is False
        assert saved['target_radius_divisor'] == 16
        assert 'rows' not in saved
        assert saved['hash_descriptor'] == saved['payloads']['sha256']
        P.authenticate_payloads(
            saved['payloads'], P.paths(prefix)[0], saved['requested_node_indices'])
        if prefix == 'pilot':
            assert saved['coverage']['all_declared_nodes_covered'] is False
        if prefix == 'n129':
            covered = saved['coverage']['all_declared_nodes_covered']
            if covered:
                assert saved['requested_node_indices'] == list(range(129))
                assert saved['completed_node_indices'] == list(range(129))
