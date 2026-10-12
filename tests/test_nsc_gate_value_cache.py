#!/usr/bin/env python3
"""Value-only family cache identity: v2 bindings, v1 immutability, budgeted dispatch."""
import copy
from hashlib import sha256
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
sys.path.insert(0, str(ROOT / 'scripts'))
import derive_nsc_ks_gate_value as V
from recursive_horizons.nsc_ks_evaluation_binding import (
    BINDING_DIMENSIONS, IMPLEMENTATION_OWNERS, LEGACY_REPLAY_MODE,
    PRODUCTION_MODE, SCHEMA_V1, SCHEMA_V2, BindingMismatch,
    authenticate_source_hashes, binding_digest, compare_bindings,
    evaluation_binding, family_identity, family_key_text,
    family_source_binding, require_complete_binding, source_identity_of,
    structured_digest, write_json_atomic)
from recursive_horizons.nsc_ks_profile_identity import profile_identity
from recursive_horizons.nsc_ks_source_envelope import computational_z_grid
from recursive_horizons.nsc_local_incoming_family import LocalIncomingFamily
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes


@pytest.fixture
def value_root(tmp_path, monkeypatch):
    (tmp_path / 'results' / 'development' / 'artifacts').mkdir(parents=True)
    monkeypatch.setattr(V, 'ROOT', tmp_path)
    return tmp_path


def hex64(seed):
    return sha256(str(seed).encode()).hexdigest()


def dict_batch(**overrides):
    row = {
        'batch_id': 'synthetic:+1:0:1',
        'panel_name': 'synthetic',
        'original_panel': 'synthetic',
        'rows': [0, 1],
        'group': 14,
        'angular_sign': 1,
        'energy_sign': 1,
        'mass': 1.5,
        'angular': 2.236067977,
        'rho_up': 1.03,
        'source_digest': hex64('source'),
        'preparation_digest': hex64('prep'),
        'covariance_digest': hex64('cov'),
        'covariance_shape': [3, 3],
        'weights_digest': hex64('weights'),
        'energies_digest': hex64('energies'),
        'initial_columns_digest': hex64('initial'),
        'initial_columns_shape': [2, 3],
        'copy_count': 2,
        'degeneracy': 12,
        'channel': {
            'group': 14, 'angular_sign': 1, 'energy_sign': 1,
            'compact_mass': 1.5, 'angular_eigenvalue': 2.236067977,
            'copy_count': 2, 'degeneracy': 12,
        },
    }
    row.update(overrides)
    return row


def synthetic_ctx(keys=((14, 1),)):
    family = LocalIncomingFamily(np.zeros((2, 8)))
    solver = dict(V.BASE_SOLVER)
    target = np.linspace(0.12, 0.18, 47)
    grid = computational_z_grid(solver['grid_nodes'], solver['length'])
    hashes = {name: hex64(name) for name in IMPLEMENTATION_OWNERS}
    hashes['results/development/nsc-incoming-source-update-v5.json'] = hex64('v5')
    hashes['results/development/history.json'] = hex64('history')
    families = {}
    for key in keys:
        batch = dict_batch(
            batch_id=f'synthetic-{key[0]}:{key[1]:+d}',
            group=key[0], angular_sign=key[1],
            source_digest=hex64(('source', key)),
            preparation_digest=hex64(('prep', key)),
            covariance_digest=hex64(('cov', key)),
            weights_digest=hex64(('w', key)),
            energies_digest=hex64(('e', key)),
            initial_columns_digest=hex64(('a', key)))
        batch['channel'] = {
            **batch['channel'], 'group': key[0], 'angular_sign': key[1],
        }
        families[family_key_text(key)] = family_source_binding([batch])
    return {
        'family': family, 'target': target, 'grid': grid, 'solver': solver,
        'identity': profile_identity(family, include_normal_window=True),
        'history_path': 'results/development/history.json',
        'hashes': hashes,
        'coeff': {
            'a': 1.0, 'r': 1.0, 'A0': 1.0, 'C': 1.0,
            'D': np.zeros(5), 'F': np.ones(3),
        },
        'baseline': np.array([0.1, -0.2]),
        'archived': {key: {} for key in keys},
        'node_count': 47, 'binding': None, 'rho_up': 1.03,
        'source_families': families,
    }


def make_binding(ctx=None):
    ctx = synthetic_ctx() if ctx is None else ctx
    return evaluation_binding(
        family=ctx['family'], target=ctx['target'], grid=ctx['grid'],
        solver=ctx['solver'], hashes=ctx['hashes'],
        history_path=ctx['history_path'],
        profile_identity_value=ctx['identity'],
        coefficients=ctx['coeff'], baseline=ctx['baseline'],
        source_families=ctx['source_families'],
        node_count=ctx['node_count'])


def write_payload(path, matter, metadata, z=None):
    arrays = {
        'matter_change': np.asarray(matter, float),
        'metadata_json': np.frombuffer(
            json.dumps(metadata, sort_keys=True).encode(), np.uint8),
    }
    if z is not None:
        arrays['z'] = np.asarray(z, float)
    path.write_bytes(deterministic_npz_bytes(arrays))
    return sha256(path.read_bytes()).hexdigest()


def write_family(directory, key, record, matter, metadata, z=None, root=None):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    json_path, npz_path = V.family_paths(directory, key)
    if record.get('schema') == SCHEMA_V2:
        metadata = {**{name: record[name] for name in (
            'family_identity', 'evaluation_identity', 'source_identity', 'binding_digest')},
            **metadata}
    digest = write_payload(npz_path, matter, metadata, z=z)
    record = dict(record)
    record['payload'] = {
        'path': str(npz_path.relative_to(root or V.ROOT)),
        'sha256': digest, 'bytes': npz_path.stat().st_size,
    }
    write_json_atomic(json_path, record, exclusive=True)
    return json_path, npz_path, record


def v2_record(binding, key, **overrides):
    record = {
        'schema': SCHEMA_V2,
        'positive_family': [int(key[0]), int(key[1])],
        'family_identity': family_identity(binding, key),
        'evaluation_identity': binding_digest(binding),
        'source_identity': binding['source_families'][family_key_text(key)]['source_identity'],
        'binding_digest': binding_digest(binding),
        'node_count': binding['node_count'],
        'matter_shape': [int(binding['node_count']), 2],
        'solver': binding['solver'],
        'tangents': 'zero',
        'profile_identity': binding['profile_identity'],
        'source_hashes': binding['source_hashes'],
    }
    record.update(overrides)
    return record


def v1_record(key, identity, hashes, **overrides):
    record = {
        'schema': SCHEMA_V1,
        'positive_family': [int(key[0]), int(key[1])],
        'tangents': 'zero',
        'profile_identity': identity,
        'source_hashes': hashes,
    }
    record.update(overrides)
    return record


def test_evaluation_binding_is_complete_and_digest_stable():
    binding = make_binding()
    require_complete_binding(binding)
    assert binding['schema'].endswith('-v2')
    assert binding_digest(binding) == binding_digest(copy.deepcopy(binding))
    assert family_identity(binding, (14, 1)) != binding_digest(binding)


def test_missing_implementation_hash_is_rejected_before_complete():
    ctx = synthetic_ctx()
    ctx['hashes'] = {
        name: ctx['hashes'][name] for name in ctx['hashes']
        if name != IMPLEMENTATION_OWNERS[-1]
    }
    with pytest.raises(ValueError, match='missing binding field'):
        make_binding(ctx)


def _refresh_source_identities(binding):
    families = binding['source_families']
    for record in families.values():
        record['source_identity'] = source_identity_of(record['batches'])
    binding['source_identity'] = structured_digest({
        key: families[key]['source_identity'] for key in sorted(families)
    })


def mutate_dimension(binding, dimension):
    mutated = copy.deepcopy(binding)
    key = next(iter(mutated['source_families']))
    if dimension == 'source':
        mutated['source_families'][key]['batches'][0]['source_digest'] = hex64('tamper-source')
        mutated['source_families'][key]['batches'][0]['covariance_digest'] = hex64('tamper-cov')
        _refresh_source_identities(mutated)
    elif dimension == 'history':
        mutated['profile_identity'] = hex64('tamper-history')
        mutated['geometry']['profile_identity'] = mutated['profile_identity']
        mutated['geometry']['w_digest'] = hex64('tamper-w')
        mutated['history_path'] = 'results/development/other-history.json'
    elif dimension == 'channel':
        mutated['source_families'][key]['batches'][0]['channel']['copy_count'] = 99
        mutated['source_families'][key]['batches'][0]['channel']['degeneracy'] = 3
    elif dimension == 'settings':
        mutated['solver']['integrator'] = 'cf4'
        mutated['solver']['max_step'] = (0.001).hex()
    elif dimension == 'nodes':
        mutated['node_count'] = 129
        mutated['target']['shape'] = [129]
        mutated['target']['digest'] = hex64('tamper-target')
        mutated['grid']['nodes'] = 128
    else:
        mutated['implementation_hashes'][IMPLEMENTATION_OWNERS[0]] = hex64('tamper-code')
        mutated['source_hashes'][IMPLEMENTATION_OWNERS[0]] = hex64('tamper-code')
    return mutated


@pytest.mark.parametrize('dimension', BINDING_DIMENSIONS)
def test_each_binding_dimension_invalidates_a_matching_cache(value_root, dimension):
    ctx = synthetic_ctx()
    binding = make_binding(ctx)
    mutated = mutate_dimension(binding, dimension)
    with pytest.raises(BindingMismatch, match=dimension):
        compare_bindings(binding, mutated)
    key = (14, 1)
    write_family(
        value_root, key, v2_record(binding, key), np.zeros((47, 2)),
        {'family_records': {'14_1_E+1': {}}}, z=ctx['target'], root=value_root)
    ctx['binding'] = mutated
    with pytest.raises((BindingMismatch, ValueError)):
        V.inspect_existing_families(ctx, value_root, mode=PRODUCTION_MODE)


def test_tampered_payload_and_descendant_metadata_are_rejected(value_root):
    ctx = synthetic_ctx()
    binding = make_binding(ctx)
    ctx['binding'] = binding
    key = (14, 1)
    json_path, npz_path, record = write_family(
        value_root, key, v2_record(binding, key), np.zeros((47, 2)),
        {'family_records': {'14_1_E+1': {}}}, z=ctx['target'], root=value_root)
    npz_path.write_bytes(npz_path.read_bytes() + b'\x00')
    with pytest.raises(BindingMismatch, match='payload'):
        V.inspect_existing_families(ctx, value_root, mode=PRODUCTION_MODE)
    write_payload(
        npz_path, np.zeros((47, 2)), {'family_records': {'14_1_E+1': {}}},
        z=ctx['target'])
    record['payload']['sha256'] = hex64('forged')
    write_json_atomic(json_path, record)
    with pytest.raises(BindingMismatch, match='payload'):
        V.inspect_existing_families(ctx, value_root, mode=PRODUCTION_MODE)


def test_mixed_v1_cache_is_rejected_before_evolution(value_root, monkeypatch):
    ctx = synthetic_ctx(keys=((14, 1), (15, 1)))
    binding = make_binding(ctx)
    ctx['binding'] = binding
    write_family(
        value_root, (14, 1),
        v1_record((14, 1), ctx['identity'], ctx['hashes']),
        np.zeros((47, 2)), {'family_records': {}}, z=ctx['target'], root=value_root)
    evolved = []

    def forbidden(*_args, **_kwargs):
        evolved.append(True)
        raise AssertionError('evolution must not start on a mixed cache')

    monkeypatch.setattr(V, 'evolve_family', forbidden)
    with pytest.raises(ValueError, match='v1 family cache'):
        V.inspect_existing_families(ctx, value_root, mode=PRODUCTION_MODE)
    assert evolved == []


def test_v1_top_register_is_immutable(tmp_path, monkeypatch):
    development = tmp_path / 'results' / 'development'
    development.mkdir(parents=True)
    (development / 'artifacts').mkdir()
    monkeypatch.setattr(V, 'ROOT', tmp_path)
    output = development / 'nsc-ks-gate-value-iterate6-47.json'
    output.write_text(json.dumps({
        'schema': SCHEMA_V1, 'merit': 0.001, 'tangents': 'zero',
        'profile_identity': 'x' * 64, 'source_hashes': {},
    }, sort_keys=True) + '\n')
    original = output.read_bytes()
    with pytest.raises(ValueError, match='immutable'):
        V.run(47, 1, 10.0, label='iterate6-47')
    assert output.read_bytes() == original


def test_existing_v2_register_is_not_overwritten(tmp_path, monkeypatch):
    development = tmp_path / 'results' / 'development'
    development.mkdir(parents=True)
    (development / 'artifacts').mkdir()
    monkeypatch.setattr(V, 'ROOT', tmp_path)
    output = development / 'nsc-ks-gate-value-v2-47.json'
    output.write_text(json.dumps({
        'schema': SCHEMA_V2, 'merit': 0.001, 'tangents': 'zero',
    }, sort_keys=True) + '\n')
    original = output.read_bytes()
    with pytest.raises(FileExistsError, match='already exists'):
        V.run(47, 1, 10.0)
    assert output.read_bytes() == original


def test_legacy_replay_resolves_git_bytes_and_production_stays_strict(tmp_path):
    relative = 'scripts/derive_nsc_ks_gate_value.py'
    historical = b'historical-value-owner\n'
    expected = sha256(historical).hexdigest()
    (tmp_path / 'scripts').mkdir()
    (tmp_path / relative).write_bytes(b'changed-working-tree\n')
    hashes = {relative: expected}
    with pytest.raises(ValueError, match='source hash mismatch'):
        authenticate_source_hashes(hashes, tmp_path, mode=PRODUCTION_MODE)

    def git_log(_relative, _max_commits):
        return ['deadbeefdeadbeefdeadbeefdeadbeefdeadbeef']

    def git_show(commit, path):
        assert path == relative
        assert commit.startswith('deadbeef')
        return historical

    origins = authenticate_source_hashes(
        hashes, tmp_path, mode=LEGACY_REPLAY_MODE,
        git_log=git_log, git_show=git_show)
    assert origins[relative] == 'deadbeefdeadbeefdeadbeefdeadbeefdeadbeef'


def test_legacy_replay_does_not_walk_result_artifacts_in_git(tmp_path):
    relative = 'results/development/nsc-incoming-source-update-v5.json'
    (tmp_path / 'results' / 'development').mkdir(parents=True)
    (tmp_path / relative).write_bytes(b'current-archive\n')
    hashes = {relative: hex64('old-archive')}

    def git_log(*_args, **_kwargs):
        raise AssertionError('numerical result history must not be walked')

    with pytest.raises(ValueError, match='legacy recorded data hash mismatch'):
        authenticate_source_hashes(
            hashes, tmp_path, mode=LEGACY_REPLAY_MODE, git_log=git_log)


def test_explicit_legacy_replay_accepts_v1_family_and_production_does_not(value_root):
    ctx = synthetic_ctx()
    key = (14, 1)
    top = {
        'schema': SCHEMA_V1, 'profile_identity': ctx['identity'],
        'tangents': 'zero', 'source_hashes': ctx['hashes'],
        'node_count': 47, 'merit': 0.1,
    }
    write_family(
        value_root, key, v1_record(key, ctx['identity'], ctx['hashes']),
        np.zeros((47, 2)), {'family_records': {}}, z=ctx['target'], root=value_root)
    inspected = V.inspect_existing_families(
        ctx, value_root, mode=LEGACY_REPLAY_MODE, recorded=top)
    assert inspected['pending'] == []
    assert len(inspected['done']) == 1
    with pytest.raises(ValueError, match='v1 family cache'):
        V.inspect_existing_families(ctx, value_root, mode=PRODUCTION_MODE)


def test_cpu_budget_stops_dispatch_and_keeps_checkpoints(tmp_path, monkeypatch):
    ctx = synthetic_ctx(keys=((14, 1), (15, 1), (16, 1)))
    binding = make_binding(ctx)
    ctx['binding'] = binding
    development = tmp_path / 'results' / 'development'
    artifacts = development / 'artifacts' / 'nsc-ks-gate-value-budget'
    development.mkdir(parents=True)
    artifacts.mkdir(parents=True)
    monkeypatch.setattr(V, 'ROOT', tmp_path)
    monkeypatch.setattr(V, 'LOCK', tmp_path / 'dirac.lock')
    monkeypatch.setattr(V, 'context', lambda *args, **kwargs: ctx)
    launched = []

    def fake_run(payload):
        key = tuple(payload[1])
        launched.append(key)
        json_path, npz_path = V.family_paths(artifacts, key)
        write_family(
            artifacts, key, v2_record(binding, key), np.zeros((47, 2)),
            {'family_records': {}}, z=ctx['target'], root=tmp_path)
        assert json_path.exists() and npz_path.exists()
        return {'positive_family': list(key), 'CPU_seconds': 10.0}

    monkeypatch.setattr(V, '_run_one', fake_run)
    monkeypatch.setattr(V, 'output_paths', lambda label: (artifacts, development / f'nsc-ks-gate-value-{label}.json'))
    result = V.run(47, 2, 15.0, label='budget', inline=True, lock_lines=[])
    assert launched == [(14, 1), (15, 1)]
    assert result['families_evolved_this_call'] == 2
    assert result['completed_positive_families'] == 2
    assert not (development / 'nsc-ks-gate-value-budget.json').exists()
    leftover = artifacts / 'family-16-+1.npz'
    leftover.write_bytes(b'orphan-checkpoint')
    with pytest.raises(ValueError, match='orphan family payload'):
        V.inspect_existing_families(ctx, artifacts, mode=PRODUCTION_MODE)
    assert leftover.exists()


def test_family_json_write_is_exclusive(tmp_path):
    path = tmp_path / 'family-14-+1.json'
    write_json_atomic(path, {'schema': SCHEMA_V2}, exclusive=True)
    original = path.read_bytes()
    with pytest.raises(FileExistsError):
        write_json_atomic(path, {'schema': SCHEMA_V2, 'tampered': True}, exclusive=True)
    assert path.read_bytes() == original


def test_default_production_label_is_v2_and_not_a_legacy_name():
    assert V.default_label(47) == 'v2-47'


def test_bound_context_rejects_changed_runtime_settings():
    ctx = synthetic_ctx()
    V.bind_context(ctx)
    ctx['solver']['rtol'] *= 0.1
    with pytest.raises(BindingMismatch, match='settings'):
        V.bind_context(ctx)


def test_payload_identity_is_checked_even_with_a_new_payload_hash(value_root):
    ctx = synthetic_ctx()
    binding = V.bind_context(ctx)
    key = (14, 1)
    write_family(
        value_root, key, v2_record(binding, key), np.zeros((47, 2)),
        {'family_records': {}, 'family_identity': hex64('different-family')},
        z=ctx['target'], root=value_root)
    with pytest.raises(BindingMismatch, match='payload identity'):
        V.inspect_existing_families(ctx, value_root, mode=PRODUCTION_MODE)


def test_implementation_hashes_follow_relative_imports(tmp_path):
    from recursive_horizons.nsc_ks_evaluation_binding import implementation_hashes

    package = tmp_path/'src/recursive_horizons'
    package.mkdir(parents=True)
    scripts = tmp_path/'scripts'
    scripts.mkdir()
    (scripts/'owner.py').write_text('from recursive_horizons.foo import value\n')
    (package/'__init__.py').write_text('__version__ = "0.1"\n')
    (package/'foo.py').write_text('from . import bar\nvalue = bar.value\n')
    (package/'bar.py').write_text('value = 1\n')
    (package/'unrelated.py').write_text('unused = 42\n')
    owners = ('scripts/owner.py',)
    before = implementation_hashes(tmp_path, owners)
    assert set(before) == {'scripts/owner.py', 'src/recursive_horizons/foo.py',
                           'src/recursive_horizons/bar.py'}
    (package/'bar.py').write_text('value = 2\n')
    after = implementation_hashes(tmp_path, owners)
    assert before['src/recursive_horizons/bar.py'] != after['src/recursive_horizons/bar.py']
    (package/'__init__.py').write_text('from .bar import value\n')
    assert 'src/recursive_horizons/__init__.py' in implementation_hashes(tmp_path, owners)
    assert V.default_label(47) != 'iterate6-47'


def test_actual_family_writer_output_roundtrips_strict_cache_validation(value_root):
    ctx = synthetic_ctx()
    binding = make_binding(ctx)
    ctx['binding'] = binding
    directory, _ = V.output_paths('writer-roundtrip')
    directory.mkdir(parents=True)
    V.persist_binding(directory, binding)
    # Exercise the real serializer. A fixture previously supplied the identity
    # omitted by this writer; scientific operator behavior is tested elsewhere.
    operator = SimpleNamespace(**{name: np.zeros(1) for name in V.KEYS})
    operator.z = ctx['target']
    operator.digest = hex64('operator')
    operator.diagnostics = {}
    V.save_family(ctx, directory, (14, 1), np.zeros((47, 2)), {}, operator)
    result = V.inspect_existing_families(ctx, directory, mode=PRODUCTION_MODE)
    assert len(result['done']) == 1
    assert result['pending'] == []
