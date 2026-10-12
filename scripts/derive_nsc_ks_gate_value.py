#!/usr/bin/env python3
"""Value-only residual of one declared history. Tangents are not evolved.

Production writes NSC-KS-GATE-VALUE-v2 with a structured evaluation binding.
v1 registers are immutable and are reconstructed only through explicit
legacy replay. Default labels are v2-{nodes}.

Owner API for calibration and tests: SCHEMA_V2, DEFAULT_LABEL_PREFIX,
default_label, bind_context, persist_binding, load_persisted_binding,
inspect_existing_families, assemble(mode=...), reconstruct_gradient(replay=...),
dispatch_with_cpu_budget, run(..., inline=, lock_lines=).
"""
import os
import re
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
import argparse
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from hashlib import sha256
import json
from pathlib import Path
import signal
import subprocess
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
sys.path.insert(0, str(ROOT/'scripts'))
import derive_nsc_ks_coupled_retained_control as R
import recursive_horizons.nsc_ks_history_evaluator as H
from recursive_horizons.nsc_ks_cutoff_bridge_evaluator import (
    signed_family_angular_square_weights)
from recursive_horizons.nsc_ks_energy_propagator import (
    KSEnergyPropagator, evolve_energy_propagator)
from recursive_horizons.nsc_ks_evaluation_binding import (
    BINDING_SCHEMA, IMPLEMENTATION_OWNERS, LEGACY_REPLAY_MODE, PRODUCTION_MODE,
    SCHEMA_V1, SCHEMA_V2, authenticate_payload, authenticate_source_hashes,
    binding_digest, classify_family_record, classify_top_record,
    compare_family_cache, compare_bindings, evaluation_binding, family_identity, family_key_text,
    family_source_binding, legacy_family_consistent, payload_descriptor,
    require_complete_binding, write_bytes_atomic, write_json_atomic,
    implementation_hashes)
from recursive_horizons.nsc_ks_family_pool import (
    acquire_dirac_lock, orphan_family_payloads, release_dirac_lock,
    _worker_entry)
from recursive_horizons.nsc_ks_profile_identity import profile_identity
from recursive_horizons.nsc_ks_signed_state import negative_angular_partner
from recursive_horizons.nsc_ks_source_envelope import (
    KSEnvelopeBinding, computational_z_grid, usual_axial_support, _sample_axial_profiles)
from recursive_horizons.nsc_ks_value_evaluator import (
    MERIT_DEFINITION, TARGET_COUNTS, assemble_residual, edge_change,
    geometry_change, matter_change, merit, target_nodes)
from recursive_horizons.nsc_local_incoming_constraints import _channel_record
from recursive_horizons.nsc_local_incoming_family import LocalIncomingFamily
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes

HISTORY = ROOT/'results/development/nsc-ks-coupled-newton-n32-truncated-iterate6.json'
LOCK = ROOT/'results/development/nsc-ks-gate-dirac.lock'
FAMILY_CAP = 7200.0
PHASE_NODES = 192
LABEL_RE = re.compile(r'^[a-z0-9-]+$')
BASE_SOLVER = {
    'grid_nodes': R.GRID_NODES,
    'degree': R.DEGREE,
    'rtol': R.OPTIONS['rtol'],
    'atol': R.OPTIONS['atol'],
    'max_step': R.OPTIONS['max_step'],
    'phase_nodes': PHASE_NODES,
    'length': 0.4,
    'integrator': 'dop853',
    'step_control': 'primal',
}
KEYS = ('energies', 'z', 'reference', 'difference', 'difference_z', 'tangent', 'tangent_z')
OWNERS = IMPLEMENTATION_OWNERS
DEFAULT_LABEL_PREFIX = 'v2'
LEGACY_REPLAY_MATH = (
    'src/recursive_horizons/nsc_ks_value_evaluator.py',
    'src/recursive_horizons/nsc_evolved_incoming_constraints.py',
    'src/recursive_horizons/nsc_compatible_history_geometry.py',
    'src/recursive_horizons/nsc_local_incoming_family.py',
    'src/recursive_horizons/nsc_dirac_source_phase.py',
    'src/recursive_horizons/nsc_ks_cutoff_bridge_evaluator.py',
    'src/recursive_horizons/nsc_ks_source_envelope.py',
    'src/recursive_horizons/nsc_ks_spacetime_variation.py',
    'src/recursive_horizons/nsc_lorentzian.py',
    'scripts/derive_nsc_evolved_incoming_constraints.py',
)


def digest(path):
    path = Path(path)
    return sha256((path if path.is_absolute() else ROOT/path).read_bytes()).hexdigest()


def require_label(label):
    if not isinstance(label, str) or LABEL_RE.fullmatch(label) is None:
        raise ValueError('label must match [a-z0-9-]+')
    return label


def output_paths(label):
    safe = require_label(label)
    artifacts = (ROOT/'results/development/artifacts').resolve()
    development = (ROOT/'results/development').resolve()
    directory = (artifacts/f'nsc-ks-gate-value-{safe}').resolve()
    output = (development/f'nsc-ks-gate-value-{safe}.json').resolve()
    if directory.parent != artifacts or output.parent != development:
        raise ValueError('label escapes the development directories')
    return directory, output


def solver_settings(overrides=None):
    settings = dict(BASE_SOLVER)
    for key, value in (overrides or {}).items():
        if key not in settings:
            raise ValueError('unknown solver setting ' + key)
        settings[key] = value
    if int(settings['grid_nodes']) not in (64, 128, 256, 512, 1024, 2048, 4096):
        raise ValueError('grid must be 64, 128, 256, 512, 1024, 2048 or 4096')
    if int(settings['degree']) not in (32, 48, 64):
        raise ValueError('energy degree must be 32, 48 or 64')
    if int(settings['phase_nodes']) not in (192, 384):
        raise ValueError('phase quadrature must be 192 or 384')
    for name in ('rtol', 'atol', 'max_step'):
        number = float(settings[name])
        if not np.isfinite(number) or number <= 0:
            raise ValueError(name + ' must be positive and finite')
        settings[name] = number
    settings['grid_nodes'] = int(settings['grid_nodes'])
    settings['degree'] = int(settings['degree'])
    settings['phase_nodes'] = int(settings['phase_nodes'])
    length = float(settings['length'])
    if not np.isfinite(length) or length <= 0:
        raise ValueError('computational length must be positive')
    settings['length'] = length
    if settings['integrator'] not in ('dop853', 'cf4'):
        raise ValueError("integrator must be 'dop853' or 'cf4'")
    if settings['step_control'] not in ('primal', 'joint'):
        raise ValueError("step control must be primal or joint")
    return settings


def family_paths(directory, key):
    group, sign = int(key[0]), int(key[1])
    if group < 0 or group > 99 or sign not in (-1, 1):
        raise ValueError('family key is not a retained angular family')
    stem = f'family-{group:02d}-{sign:+d}'
    return directory/(stem+'.json'), directory/(stem+'.npz')


def history_path(path=None):
    chosen = HISTORY if path is None else Path(path)
    if not chosen.is_absolute():
        chosen = ROOT/chosen
    resolved = chosen.resolve()
    development = (ROOT/'results/development').resolve()
    if resolved.parent != development or resolved.suffix != '.json':
        raise ValueError('history must be one json register in results/development')
    return resolved


def context(node_count, settings=None, history=None):
    base = R.context()
    parent_path = history_path(history)
    parent = json.loads(parent_path.read_text())
    family = LocalIncomingFamily(np.asarray(parent['history']['coefficients'], float))
    if profile_identity(family, include_normal_window=True) != parent['profile_identity']:
        raise ValueError('history identity does not match its coefficients')
    solver = solver_settings(settings)
    target = np.asarray(target_nodes(family, node_count), float)
    return {
        **base, 'family': family, 'target': target, 'node_count': int(node_count),
        'solver': solver,
        'grid': computational_z_grid(solver['grid_nodes'], solver['length']),
        'identity': parent['profile_identity'],
        'history_path': str(parent_path.relative_to(ROOT)),
        'hashes': {**base['hashes'], **base['archive'].input_hashes,
                   **implementation_hashes(ROOT),
                   str(parent_path.relative_to(ROOT)): digest(parent_path)},
        'binding': None,
    }


def source_families_from_ctx(ctx):
    families = {}
    for key in sorted(ctx['archived']):
        entries = []
        for batch, channel in ctx['archive'].family_entries(key):
            entries.append((batch, _channel_record('_', {'_': channel})))
        families[family_key_text(key)] = family_source_binding(entries)
    return families


def bind_context(ctx):
    """Structured production binding. Missing fields fail before any dispatch."""
    previous = ctx.get('binding')
    if (previous is not None and '_archive_identity' in ctx
            and ctx['_archive_identity'] != id(ctx.get('archive'))):
        raise ValueError('source archive changed inside an evaluation context')
    families = ctx.get('source_families') or source_families_from_ctx(ctx)
    binding = evaluation_binding(
        family=ctx['family'], target=ctx['target'], grid=ctx['grid'],
        solver=ctx['solver'], hashes=ctx['hashes'],
        history_path=ctx['history_path'],
        profile_identity_value=ctx['identity'],
        coefficients=ctx['coeff'], baseline=ctx['baseline'],
        source_families=families, node_count=ctx['node_count'])
    if previous is not None:
        compare_bindings(previous, binding)
    ctx['binding'] = binding
    ctx['_archive_identity'] = id(ctx.get('archive'))
    ctx['source_families'] = families
    ctx['evaluation_identity'] = binding_digest(binding)
    return binding


def binding_file(directory):
    return Path(directory)/'binding.json'


def persist_binding(directory, binding):
    path = binding_file(directory)
    raw = json.dumps(binding, sort_keys=True, indent=2, allow_nan=False)+'\n'
    if path.exists():
        if path.read_text() != raw:
            raise ValueError('checkpoint belongs to another source/history/numerics')
        return path
    write_bytes_atomic(path, raw.encode(), exclusive=True)
    return path


def load_persisted_binding(directory):
    path = binding_file(directory)
    if not path.exists():
        raise ValueError('missing binding field: evaluation binding file')
    return require_complete_binding(json.loads(path.read_text()))


def contract_family(ctx, key, operator):
    archived = ctx['archived'][key]
    entries = ctx['archive'].family_entries(key)
    tagged = []
    for batch, channel in entries:
        normalized = _channel_record('_', {'_': channel})
        tagged.append((
            f'{batch.group}_{batch.angular_sign}_E{batch.energy_sign:+d}:{batch.panel_name}',
            batch, normalized))
    groups = H.group_signed_operator_families(tagged)
    if len(groups) != 1:
        raise ValueError('one operator channel per positive angular family required')
    group = groups[0]
    metric = ctx['family'].metric()
    if operator.tangent.shape[1] != 0:
        raise ValueError('value operator must carry zero tangent directions')
    corrections, family_records, states = {}, {}, {}

    def record_signed(signed_id, change, channel, energy_sign, batch_ids):
        corrections[signed_id] = corrections.get(signed_id, 0) + change
        if signed_id in family_records:
            family_records[signed_id]['batch_ids'].extend(batch_ids)
            return
        family_records[signed_id] = {
            'group': channel['group'], 'angular_sign': channel['angular_sign'],
            'energy_sign': energy_sign, 'compact_mass': channel['compact_mass'],
            'angular_eigenvalue': channel['angular_eigenvalue'],
            'signed_angular': channel['angular_sign'] * channel['angular_eigenvalue'],
            'copy_count': channel['copy_count'], 'degeneracy': channel['degeneracy'],
            'multiplicity': channel['multiplicity'],
            'ell0_pure_radius_identity': channel['ell0_pure_radius_identity'],
            'batch_ids': list(batch_ids),
        }

    for batch_id, batch, channel in group['applies']:
        prepared = operator.apply(batch.source, batch.initial_columns)
        prepared.require_history(metric)
        change = matter_change(prepared, channel, ctx['coeff'])
        signed_id = f'{channel["group"]}_{channel["angular_sign"]}_E{batch.energy_sign:+d}'
        record_signed(signed_id, change, channel, batch.energy_sign, [batch_id])
        states[batch_id] = prepared
    for pos_id, neg_id, neg_batch, neg_channel in group['maps']:
        partner = negative_angular_partner(states[pos_id], neg_batch.source)
        partner.require_history(metric)
        change = matter_change(partner, neg_channel, ctx['coeff'])
        signed_id = f'{neg_channel["group"]}_{neg_channel["angular_sign"]}_E{neg_batch.energy_sign:+d}'
        record_signed(signed_id, change, neg_channel, neg_batch.energy_sign, [neg_id])
    total = sum(corrections.values())
    return total, family_records, operator


def evolve_family(ctx, key):
    archived = ctx['archived'][key]
    interval = R.interpolation_interval(archived)
    entries = ctx['archive'].family_entries(key)
    tagged = []
    for batch, channel in entries:
        normalized = _channel_record('_', {'_': channel})
        tagged.append((
            f'{batch.group}_{batch.angular_sign}_E{batch.energy_sign:+d}:{batch.panel_name}',
            batch, normalized))
    group = H.group_signed_operator_families(tagged)[0]
    first = next(batch for _, batch, _ in group['applies'])
    solver = ctx['solver']
    operator = evolve_energy_propagator(
        interval, solver['degree'], ctx['family'], ctx['grid'], ctx['target'],
        first.mass, first.angular, first.rho_up,
        axial_support=usual_axial_support(), tangents='zero',
        rtol=solver['rtol'], atol=solver['atol'], max_step=solver['max_step'],
        integrator=solver.get('integrator', 'dop853'),
        step_control=solver.get('step_control', 'primal'))
    return contract_family(ctx, key, operator)


def save_family(ctx, directory, key, total, family_records, operator, cpu_seconds=None,
                inherited_from=None):
    target_json, target_npz = family_paths(directory, key)
    if target_json.exists() or target_npz.exists():
        raise FileExistsError('value family output already exists')
    binding = ctx.get('binding')
    if binding is None and binding_file(directory).exists():
        binding = load_persisted_binding(directory)
        ctx['binding'] = binding
    binding = bind_context(ctx) if binding is None else require_complete_binding(binding)
    arrays = {name: np.asarray(getattr(operator, name)) for name in KEYS}
    arrays['matter_change'] = np.asarray(total, float)
    meta = {'operator': {'digest': operator.digest, 'diagnostics': dict(operator.diagnostics)},
            'family_records': family_records,
            'family_identity': family_identity(binding, key),
            'evaluation_identity': binding_digest(binding),
            'source_identity': binding['source_families'][family_key_text(key)]['source_identity'],
            'binding_digest': binding_digest(binding)}
    arrays['metadata_json'] = np.frombuffer(json.dumps(meta, sort_keys=True).encode(), np.uint8)
    raw = deterministic_npz_bytes(arrays)
    write_bytes_atomic(target_npz, raw, exclusive=True)
    record = {
        'schema': SCHEMA_V2,
        'positive_family': [int(key[0]), int(key[1])],
        'family_identity': family_identity(binding, key),
        'evaluation_identity': binding_digest(binding),
        'source_identity': binding['source_families'][family_key_text(key)]['source_identity'],
        'binding_digest': binding_digest(binding),
        'node_count': ctx['node_count'],
        'matter_shape': [int(ctx['node_count']), 2],
        'solver': binding['solver'],
        'operator_digest': operator.digest,
        'matter_change_maxima': np.max(np.abs(total), axis=0).tolist(),
        'CPU_seconds': cpu_seconds,
        'new_operator_solves': 1 if inherited_from is None else 0,
        'tangents': 'zero',
        'source_hashes': ctx['hashes'],
        'profile_identity': ctx['identity'],
        'payload': payload_descriptor(target_npz, ROOT),
    }
    if inherited_from is not None:
        record['inherited_from'] = inherited_from
    write_json_atomic(target_json, record, exclusive=True)
    return record


def _run_one(payload):
    started = time.process_time()
    node_count, key, directory, settings, history = payload
    ctx = context(node_count, settings, history)
    directory = Path(directory)
    # A worker must validate its own current inputs against the parent's
    # checkpoint before spending scientific CPU, not copy the old label.
    compare_bindings(load_persisted_binding(directory), bind_context(ctx))

    def stop(*_args):
        raise TimeoutError('one value-only family exceeded its CPU allocation')

    previous = signal.signal(signal.SIGPROF, stop)
    signal.setitimer(signal.ITIMER_PROF, FAMILY_CAP)
    try:
        total, family_records, operator = evolve_family(ctx, tuple(key))
    finally:
        signal.setitimer(signal.ITIMER_PROF, 0)
        signal.signal(signal.SIGPROF, previous)
    cpu_seconds = time.process_time() - started
    record = save_family(
        ctx, directory, tuple(key), total, family_records, operator, cpu_seconds)
    return {'positive_family': record['positive_family'], 'CPU_seconds': record['CPU_seconds']}


def _rho_up(ctx):
    if 'rho_up' in ctx:
        return float(ctx['rho_up'])
    first = ctx['archive'].family_entries(next(iter(ctx['archived'])))[0][0]
    return float(first.rho_up)


def inspect_existing_families(ctx, directory, *, mode, recorded=None):
    """Authenticate stored families. Never deletes orphan payloads."""
    directory = Path(directory)
    orphan_family_payloads(directory)
    expected = bind_context(ctx) if mode == PRODUCTION_MODE else None
    pending = []
    done = []
    schemas = set()
    matter = np.zeros((len(ctx['target']), 2), float)
    family_records = {}
    for key in sorted(ctx['archived']):
        path, payload = family_paths(directory, key)
        if not path.exists():
            if payload.exists():
                raise ValueError('orphan family payload: ' + payload.name)
            pending.append(key)
            continue
        record = json.loads(path.read_text())
        schema = classify_family_record(record)
        schemas.add(schema)
        arrays = authenticate_payload(record, payload, ROOT, target=ctx['target'])
        compare_family_cache(
            record, expected, key, mode=mode, payload_arrays=arrays)
        if mode == LEGACY_REPLAY_MODE:
            if recorded is None:
                raise ValueError('explicit legacy replay requires the top record')
            if record.get('profile_identity') != ctx['identity']:
                raise ValueError('legacy family belongs to a different evaluation')
            legacy_family_consistent(recorded, record)
        stored = arrays.get('metadata') or {}
        matter = matter + arrays['matter_change']
        done.append(record)
        overlap = set(family_records) & set(stored.get('family_records', ()))
        if overlap:
            raise ValueError('signed family counted twice')
        family_records.update(stored.get('family_records', {}))
    if mode == PRODUCTION_MODE and SCHEMA_V1 in schemas:
        raise ValueError(
            'v1 family cache is not a production evaluation; choose a new v2 label')
    if mode == PRODUCTION_MODE and schemas and schemas != {SCHEMA_V2}:
        raise ValueError('mixed value-family cache is rejected before evolution')
    if mode == LEGACY_REPLAY_MODE and SCHEMA_V2 in schemas:
        raise ValueError('explicit legacy replay cannot mix v2 production caches')
    return {
        'pending': pending, 'done': done, 'matter': matter,
        'family_records': family_records, 'schemas': schemas,
        'binding': expected,
    }


def assemble(ctx, directory, *, mode=PRODUCTION_MODE, recorded=None):
    inspected = inspect_existing_families(
        ctx, directory, mode=mode, recorded=recorded)
    done = inspected['done']
    family_records = inspected['family_records']
    matter = inspected['matter']
    complete = len(done) == len(ctx['archived'])
    if complete and inspected['pending']:
        raise ValueError('complete value assembly still has pending families')
    geometry = edge = gradient = None
    if complete:
        if mode == PRODUCTION_MODE and inspected['binding'] is None:
            raise ValueError('missing binding field: evaluation binding')
        if len(family_records) != 120:
            raise ValueError('one hundred twenty explicit signed families required')
        weights = signed_family_angular_square_weights(family_records)
        if not np.allclose(weights, R.EXPECTED_WEIGHT, rtol=0, atol=1e-6):
            raise ValueError('signed angular-square weights must be 107952 for each energy sign')
        geometry = geometry_change(ctx['family'], ctx['target'], ctx['coeff'])
        edge = edge_change(
            ctx['family'], ctx['target'], weights, ctx['coeff']['a'], _rho_up(ctx),
            ctx['solver']['phase_nodes'])
        gradient = assemble_residual(ctx['baseline'], geometry, matter, edge)
    schema = SCHEMA_V2 if mode == PRODUCTION_MODE else SCHEMA_V1
    result = {
        'schema': schema,
        'accountable_author': 'Douglas Ek',
        'status': (
            'OPEN: value-only residual; uncertified physical gate'
            if complete else 'OPEN: partial value-only residual'),
        'merit_definition': MERIT_DEFINITION,
        'profile_identity': ctx['identity'],
        'history_path': ctx['history_path'],
        'node_count': ctx['node_count'],
        'solver': ctx['solver'],
        'tangents': 'zero',
        'completed_positive_families': len(done),
        'required_positive_families': len(ctx['archived']),
        'gradient': gradient,
        'constraint_maxima': None if gradient is None else np.max(np.abs(gradient), axis=0).tolist(),
        'merit': None if gradient is None else merit(gradient),
        'families_evolved_this_call': None,
        'physical_EXISTENCE_certificate': False,
        'physical_NONEXISTENCE_certificate': False,
        'prediction_is_not_a_residual': True,
        'source_hashes': ctx['hashes'] if mode == PRODUCTION_MODE else (
            recorded or {}).get('source_hashes', ctx['hashes']),
    }
    if mode == PRODUCTION_MODE:
        binding = inspected['binding'] or bind_context(ctx)
        result['evaluation_identity'] = binding_digest(binding)
        result['binding_schema'] = BINDING_SCHEMA
        result['source_identity'] = binding['source_identity']
        result['input_hashes'] = {
            str(binding_file(directory).relative_to(ROOT)): digest(binding_file(directory)),
            **{str(family_paths(directory, tuple(row['positive_family']))[0].relative_to(ROOT)):
               digest(family_paths(directory, tuple(row['positive_family']))[0])
               for row in done},
        }
        result['auxiliary_payloads'] = [row['payload'] for row in done]
    return result


def history_named_by(result):
    matched = []
    for key in result['source_hashes']:
        if not key.endswith('.json') or '/artifacts/' in key:
            continue
        path = ROOT/key
        if not path.is_file():
            continue
        payload = json.loads(path.read_text())
        if payload.get('profile_identity') == result['profile_identity'] and 'history' in payload:
            matched.append(path)
    if len(matched) != 1:
        raise ValueError('value register does not name one history')
    return matched[0]


def default_label(node_count):
    return f'{DEFAULT_LABEL_PREFIX}-{int(node_count)}'


def reconstruct_gradient(label, *, replay=None, git_log=None, git_show=None):
    """Nodal residual from a finished value-only run. No Dirac.

    v1 registers take the explicit historical replay path unless replay=False.
    Production reconstruction stays strict on current working-tree hashes.
    """
    directory, output = output_paths(label)
    result = json.loads(output.read_text())
    schema = classify_top_record(result)
    if replay is None:
        replay = schema == SCHEMA_V1
    if schema == SCHEMA_V1 and not replay:
        raise ValueError('v1 value register is not a production cache')
    if replay:
        return _reconstruct_legacy_gradient(
            result, directory, record_path=output, git_log=git_log, git_show=git_show)
    if schema != SCHEMA_V2:
        raise ValueError('unknown value register schema')
    history = result.get('history_path') or history_named_by(result)
    ctx = context(result['node_count'], result['solver'], history)
    authenticate_source_hashes(result['source_hashes'], ROOT, mode=PRODUCTION_MODE)
    authenticate_source_hashes(result['input_hashes'], ROOT, mode=PRODUCTION_MODE)
    assembled = assemble(ctx, directory, mode=PRODUCTION_MODE)
    for key in ('evaluation_identity', 'source_identity', 'profile_identity',
                'input_hashes', 'auxiliary_payloads'):
        if assembled[key] != result[key]:
            raise ValueError('production replay differs: ' + key)
    return _finish_reconstructed_gradient(assembled, result)


def _finish_reconstructed_gradient(assembled, result):
    gradient = assembled.get('gradient')
    if gradient is None:
        raise ValueError('value run is not complete')
    if merit(gradient) != result['merit']:
        raise ValueError('reconstructed merit does not match the register')
    return np.asarray(gradient, float)


def _reconstruct_legacy_gradient(result, directory, *, record_path, git_log=None, git_show=None):
    # Recorded fields can be contracted without rerunning the old propagator.
    # The mathematics that IS executed must still be the recorded version.
    relative_record = str(Path(record_path).relative_to(ROOT))
    commit = subprocess.check_output(
        ['git', '-C', str(ROOT), 'log', '-1', '--format=%H', '--', relative_record],
        text=True, timeout=30).strip()
    if not commit:
        raise ValueError('legacy replay requires a committed historical register')
    frozen_record = subprocess.check_output(
        ['git', '-C', str(ROOT), 'show', f'{commit}:{relative_record}'], timeout=30)
    if frozen_record != Path(record_path).read_bytes():
        raise ValueError('legacy value register differs from its pinned Git artifact')
    for name in LEGACY_REPLAY_MATH:
        frozen = subprocess.check_output(
            ['git', '-C', str(ROOT), 'show', f'{commit}:{name}'], timeout=30)
        if sha256(frozen).hexdigest() != digest(name):
            raise ValueError('legacy replay mathematics changed; use the pinned checkout: ' + name)
    authenticate_source_hashes(
        result['source_hashes'], ROOT, mode=LEGACY_REPLAY_MODE,
        git_log=git_log, git_show=git_show)
    if result.get('tangents') != 'zero':
        raise ValueError('legacy value register tangents are not zero')
    history = history_named_by(result)
    ctx = context(result['node_count'], result.get('solver'), history)
    # This path contracts already recorded v1 fields; it launches no solver.
    # Retain the old settings exactly instead of relabeling them as v2.
    if result.get('solver') is not None:
        ctx['solver'] = dict(result['solver'])
    if ctx['identity'] != result['profile_identity']:
        raise ValueError('legacy history identity does not match the top record')
    if int(ctx['node_count']) != int(result['node_count']):
        raise ValueError('legacy node count does not match the top record')
    if result.get('solver') is not None and ctx['solver'] != result['solver']:
        raise ValueError('legacy solver settings do not match the top record')
    assembled = assemble(
        ctx, directory, mode=LEGACY_REPLAY_MODE, recorded=result)
    return _finish_reconstructed_gradient(assembled, result)


def dispatch_with_cpu_budget(function, items, workers, cpu_budget, *, inline=False):
    """Launch family jobs until completed work meets the CPU budget.

    In-flight workers may finish after the budget is crossed. Checkpoints are
    kept. Orphan payloads are not deleted.
    """
    items = list(items)
    if not items:
        return ()
    budget = float(cpu_budget)
    if not np.isfinite(budget) or budget <= 0:
        raise ValueError('CPU budget must be positive')
    if inline:
        done = []
        used = 0.0
        for item in items:
            if used >= budget:
                break
            result = function(item)
            done.append(result)
            used += float(result['CPU_seconds'])
        return tuple(done)
    count = min(int(workers), len(items))
    if count < 1:
        raise ValueError('at least one family worker required')
    queue = list(items)
    in_flight = {}
    done = []
    used = 0.0

    def launch(pool):
        while queue and len(in_flight) < count and used < budget:
            item = queue.pop(0)
            in_flight[pool.submit(_worker_entry, function, item)] = item

    with ProcessPoolExecutor(max_workers=count) as pool:
        launch(pool)
        while in_flight:
            finished, _ = wait(tuple(in_flight), return_when=FIRST_COMPLETED)
            for fut in finished:
                in_flight.pop(fut)
                result = fut.result()
                done.append(result)
                used += float(result['CPU_seconds'])
            if used < budget:
                launch(pool)
    return tuple(done)


def refuse_existing_register(output):
    output = Path(output)
    if not output.exists():
        return
    existing = json.loads(output.read_text())
    schema = classify_top_record(existing)
    if schema == SCHEMA_V1:
        raise ValueError(
            'v1 value register is immutable; explicit legacy replay only')
    raise FileExistsError('v2 value register already exists; use --check or a new label')


def run(node_count, workers, cpu_budget, label=None, settings=None, history=None,
        *, inline=False, lock_lines=None):
    label = default_label(node_count) if label is None else require_label(label)
    directory, output = output_paths(label)
    refuse_existing_register(output)
    directory.mkdir(parents=True, exist_ok=True)
    ctx = context(node_count, settings, history)
    binding = bind_context(ctx)
    persist_binding(directory, binding)
    inspected = inspect_existing_families(ctx, directory, mode=PRODUCTION_MODE)
    pending = [
        (node_count, list(key), str(directory), ctx['solver'], ctx['history_path'])
        for key in inspected['pending']
    ]
    acquire_dirac_lock(LOCK, 'derive_nsc_ks_gate_value.py', os.getpid(), lines=lock_lines)
    started = time.perf_counter()
    cpu = time.process_time()
    launched = ()
    try:
        if pending:
            launched = dispatch_with_cpu_budget(
                _run_one, pending, workers, cpu_budget, inline=inline)
    finally:
        release_dirac_lock(LOCK, os.getpid())
    result = assemble(ctx, directory, mode=PRODUCTION_MODE)
    result.pop('gradient', None)
    result['wall_seconds'] = time.perf_counter() - started
    result['parent_cpu_seconds'] = time.process_time() - cpu
    result['families_evolved_this_call'] = len(launched)
    result['cpu_budget'] = float(cpu_budget)
    if result['completed_positive_families'] == result['required_positive_families']:
        if output.exists():
            raise FileExistsError('v2 value register already exists')
        write_json_atomic(output, result, exclusive=True, indent=2)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--run', action='store_true')
    mode.add_argument('--check', action='store_true')
    parser.add_argument('--nodes', type=int, default=47)
    parser.add_argument('--label', default=None)
    parser.add_argument('--workers', type=int, default=9)
    parser.add_argument('--cpu-budget', type=float, default=14400.0)
    parser.add_argument('--grid', type=int)
    parser.add_argument('--degree', type=int)
    parser.add_argument('--rtol', type=float)
    parser.add_argument('--atol', type=float)
    parser.add_argument('--max-step', type=float)
    parser.add_argument('--phase-nodes', type=int)
    parser.add_argument('--length', type=float)
    parser.add_argument('--integrator', choices=('dop853', 'cf4'))
    parser.add_argument('--history', default=None)
    args = parser.parse_args()
    overrides = {}
    for key, value in (
            ('grid_nodes', args.grid), ('degree', args.degree), ('rtol', args.rtol),
            ('atol', args.atol), ('max_step', args.max_step),
            ('phase_nodes', args.phase_nodes), ('length', args.length),
            ('integrator', args.integrator)):
        if value is not None:
            overrides[key] = value
    label = args.label or default_label(args.nodes)
    directory, output = output_paths(label)
    if args.run:
        result = run(args.nodes, args.workers, args.cpu_budget, label, overrides or None, args.history)
    else:
        saved = json.loads(output.read_text())
        schema = classify_top_record(saved)
        if schema == SCHEMA_V1:
            reconstruct_gradient(label, replay=True)
            result = saved
        else:
            ctx = context(args.nodes, overrides or None, args.history)
            result = assemble(ctx, directory, mode=PRODUCTION_MODE)
            result.pop('gradient', None)
            if saved['merit'] != result['merit'] or saved['constraint_maxima'] != result['constraint_maxima']:
                raise ValueError('value residual replay changed')
            if saved.get('evaluation_identity') != result.get('evaluation_identity'):
                raise ValueError('value residual replay changed')
    print(json.dumps({
        'status': result['status'],
        'merit': result['merit'],
        'constraint_maxima': result['constraint_maxima'],
        'completed_positive_families': result['completed_positive_families'],
        'wall_seconds': result.get('wall_seconds'),
        'merit_definition': result['merit_definition'],
    }, indent=2))
