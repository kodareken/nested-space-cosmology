#!/usr/bin/env python3
"""CURRENT-HISTORY nodal phase-value enclosure v3. No Dirac or source evolution.

Successor of scripts/derive_nsc_ks_current_phase_accuracy_v2.py. Reuses that
owner's helpers, context, true-gradient contraction and production 192/384
path. The directed enclosure is the same function with a declared tighter
target_radius = unit_angular_target()/16. V2 module globals are not patched.
Between-node phase remainder stays None. The physical local gate stays OPEN.
"""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path
import signal
import sys
import time
import types

from flint import arb, ctx

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
sys.path.insert(0, str(ROOT / 'scripts'))

import derive_nsc_ks_current_phase_accuracy_v2 as V2
from recursive_horizons.nsc_dirac_source_phase_bound import (
    directed_source_phase_z_enclosure, pack_upper, unit_angular_target,
    unpack_upper)

SCHEMA = 'NSC-KS-CURRENT-PHASE-ACCURACY-v3'
PAYLOAD_SCHEMA = 'NSC-KS-CURRENT-PHASE-ACCURACY-PAYLOADS-v3'
TARGET_RADIUS_DIVISOR = 16
DEFAULT_HISTORY = V2.DEFAULT_HISTORY
BITS = V2.BITS
PHASE_COMPONENT_TOLERANCE = V2.PHASE_COMPONENT_TOLERANCE
GAUSS_COUNTS = V2.GAUSS_COUNTS
ANGULAR_SQUARE_WEIGHT = V2.ANGULAR_SQUARE_WEIGHT
ORIGINAL_RHO_UP = V2.ORIGINAL_RHO_UP
OWNERS = (
    'scripts/derive_nsc_ks_current_phase_accuracy_v3.py',
    'docs/nsc-ks-current-phase-accuracy-v3.md',
    'scripts/derive_nsc_ks_current_phase_accuracy_v2.py',
    'src/recursive_horizons/nsc_dirac_source_phase_bound.py',
    'src/recursive_horizons/nsc_dirac_source_phase.py',
    'src/recursive_horizons/nsc_ks_cutoff_bridge_evaluator.py',
    'src/recursive_horizons/nsc_ks_value_evaluator.py',
    'src/recursive_horizons/nsc_ks_ball_geometry.py',
    'src/recursive_horizons/nsc_ks_ball_operator.py',
    'src/recursive_horizons/nsc_ks_profile_identity.py',
    'src/recursive_horizons/nsc_local_incoming_family.py',
)
RUNTIME_KEYS = frozenset(V2.RUNTIME_KEYS) | frozenset({
    'max_cpu_seconds_per_completed_node',
})
CHECKPOINT_KEYS = (
    'schema', 'bits', 'phase_component_tolerance',
    'gauss_nodes_per_normal_panel', 'exact_angular_square_weight',
    'one_direction_metric', 'target_radius_divisor', 'target_radius_upper',
    'signed_angular_square_weights',
)

digest = V2.digest
write_once = V2.write_once
select_indices = V2.select_indices
node_file = V2.node_file
binding_file = V2.binding_file
persist_binding = V2.persist_binding
true_edge_gradient = V2.true_edge_gradient
directed_distance_upper = V2.directed_distance_upper
numeric_contraction = V2.numeric_contraction
contraction_report = V2.contraction_report
one_direction = V2.one_direction
load_history = V2.load_history
collocation = V2.collocation
original_archive_rho_up = V2.original_archive_rho_up
production_angular_weights = V2.production_angular_weights
parse_indices = V2.parse_indices


def declared_target_radius():
    """unit_angular_target()/16 at the enclosure working precision."""
    return unit_angular_target() / arb(int(TARGET_RADIUS_DIVISOR))


def packed_declared_target():
    with ctx.workprec(BITS):
        return pack_upper(declared_target_radius())


def packed_inherited_target():
    with ctx.workprec(BITS):
        return pack_upper(unit_angular_target())


class _TighterDirectedEnclosure:
    """Private enclosure callable. Does not assign V2 module globals."""

    def __init__(self):
        self._inner = directed_source_phase_z_enclosure
        self.divisor = TARGET_RADIUS_DIVISOR

    def __call__(self, *args, **kwargs):
        if kwargs.get('target_radius') is None:
            kwargs['target_radius'] = unit_angular_target() / arb(int(self.divisor))
        return self._inner(*args, **kwargs)


def _bind_v2_compute_node(enclosure):
    """Copy v2 compute_node into a private globals dict with our enclosure."""
    copied = dict(V2.compute_node.__globals__)
    copied['directed_source_phase_z_enclosure'] = enclosure
    rebound = types.FunctionType(
        V2.compute_node.__code__, copied, '_compute_node_core',
        V2.compute_node.__defaults__, V2.compute_node.__closure__)
    rebound.__kwdefaults__ = V2.compute_node.__kwdefaults__
    rebound.__module__ = __name__
    rebound.__qualname__ = '_compute_node_core'
    return rebound, copied


_TIGHTER_ENCLOSURE = _TighterDirectedEnclosure()
_compute_node_core, _COMPUTE_NODE_GLOBALS = _bind_v2_compute_node(_TIGHTER_ENCLOSURE)


def without_runtime(value):
    if isinstance(value, dict):
        return {key: without_runtime(item) for key, item in value.items()
                if key not in RUNTIME_KEYS}
    if isinstance(value, list):
        return [without_runtime(item) for item in value]
    return value


def canonical(value):
    return json.loads(json.dumps(without_runtime(value), sort_keys=True, allow_nan=False))


def paths(prefix):
    safe = V2.require_prefix(prefix)
    artifacts = (ROOT / 'results/development/artifacts').resolve()
    development = (ROOT / 'results/development').resolve()
    directory = (artifacts / f'nsc-ks-current-phase-accuracy-v3-{safe}').resolve()
    output = (development / f'nsc-ks-current-phase-accuracy-v3-{safe}.json').resolve()
    if directory.parent != artifacts or output.parent != development:
        raise ValueError('prefix escapes the development directories')
    return directory, output


def repository_relative(path):
    path = Path(path).resolve()
    try:
        text = path.relative_to(ROOT.resolve()).as_posix()
    except ValueError:
        text = path.as_posix()
    if not text or text == '.':
        raise ValueError('payload path must be repository-relative')
    return text


def resolve_declared_path(declared):
    candidate = Path(declared)
    return candidate if candidate.is_absolute() else (ROOT / candidate).resolve()


def file_descriptor(path):
    path = Path(path)
    raw = path.read_bytes()
    return {
        'path': repository_relative(path),
        'sha256': sha256(raw).hexdigest(),
        'bytes': len(raw),
    }


def artifact_descriptor(directory, indices):
    """Recursive payload tree: directory digest is the hash of child hashes."""
    binding = file_descriptor(binding_file(directory))
    nodes = []
    for index in indices:
        item = file_descriptor(node_file(directory, index))
        item['node'] = int(index)
        nodes.append(item)
    listing = json.dumps({
        'binding': binding['sha256'],
        'nodes': [[item['node'], item['sha256'], item['bytes']] for item in nodes],
    }, sort_keys=True, separators=(',', ':')).encode()
    return {
        'schema': PAYLOAD_SCHEMA,
        'directory': repository_relative(directory),
        'sha256': sha256(listing).hexdigest(),
        'bytes': len(listing),
        'listing_sha256': sha256(listing).hexdigest(),
        'binding': binding,
        'nodes': nodes,
    }


def visit_hash_descriptors(value):
    """Every nested {path, sha256} object must match live bytes."""
    if isinstance(value, dict):
        if 'path' in value and 'sha256' in value:
            path = resolve_declared_path(value['path'])
            raw = path.read_bytes()
            digest_value = sha256(raw).hexdigest()
            if digest_value != value['sha256']:
                raise ValueError('payload hash descriptor changed: ' + value['path'])
            if value.get('bytes') is not None and int(value['bytes']) != len(raw):
                raise ValueError('payload size descriptor changed: ' + value['path'])
        for item in value.values():
            visit_hash_descriptors(item)
    elif isinstance(value, list):
        for item in value:
            visit_hash_descriptors(item)


def authenticate_payloads(payloads, directory, indices):
    recomputed = artifact_descriptor(directory, indices)
    if canonical(payloads) != canonical(recomputed):
        raise ValueError('recursive payload hash descriptor changed')
    visit_hash_descriptors(payloads)
    return recomputed


def source_hashes(extra=None):
    hashes = {name: digest(name) for name in OWNERS}
    hashes.update({name: digest(name) for name in (
        V2.BRANCH, V2.INVENTORY, V2.RETAINED, V2.RHO_UP_LEDGER)})
    if extra:
        hashes.update(extra)
    return hashes


def load_context(history, node_count):
    ctx_data = V2.load_context(history, node_count)
    extra = {
        str(ctx_data['parent_path'].relative_to(ROOT)): digest(ctx_data['parent_path']),
        **ctx_data['archive'].input_hashes,
    }
    ctx_data['hashes'] = source_hashes(extra)
    return ctx_data


def evaluation_binding(identity, history, nodes, indices, rho_up, weights, a_numeric,
                       hashes):
    binding = V2.evaluation_binding(
        identity, history, nodes, indices, rho_up, weights, a_numeric, hashes)
    binding['schema'] = SCHEMA
    binding['parent_schema'] = V2.SCHEMA
    binding['target_radius_divisor'] = TARGET_RADIUS_DIVISOR
    binding['target_radius_upper'] = packed_declared_target()
    binding['inherited_unit_angular_target_upper'] = packed_inherited_target()
    binding['production_axial_scale_hex'] = float(a_numeric).hex()
    binding['source_hashes'] = hashes
    return binding


def annotate_node(row):
    packed = packed_declared_target()
    row = dict(row)
    row['schema'] = SCHEMA
    row['target_radius_divisor'] = TARGET_RADIUS_DIVISOR
    row['target_radius_upper'] = packed
    enclosure = row.get('enclosure')
    if isinstance(enclosure, dict):
        enclosure = dict(enclosure)
        enclosure['target_radius_divisor'] = TARGET_RADIUS_DIVISOR
        enclosure['target_radius_upper'] = packed
        row['enclosure'] = enclosure
    return row


def compute_node(*args, **kwargs):
    return annotate_node(_compute_node_core(*args, **kwargs))


def require_node_matches_binding(row, binding, index):
    if row.get('budget_exceeded'):
        raise ValueError('incomplete budget-exceeded node file must not be reused')
    for key in CHECKPOINT_KEYS:
        if key not in binding:
            continue
        if row.get(key) != binding[key]:
            raise ValueError(f'node checkpoint {key} changed')
    axial = row.get('production_axial_scale')
    if axial is not None:
        expected = binding.get('production_axial_scale_hex')
        if expected is None:
            expected = float(binding['production_axial_scale']).hex()
        if float(axial).hex() != expected:
            raise ValueError('node checkpoint production axial scale changed')
    if row.get('complete'):
        quadrature = row.get('quadrature')
        if not isinstance(quadrature, list) or not quadrature:
            raise ValueError('complete node checkpoint missing quadrature')
        counts = [item.get('gauss_nodes_per_normal_panel') for item in quadrature]
        if counts != list(GAUSS_COUNTS):
            raise ValueError('node checkpoint gauss quadrature changed')
        if row.get('edge_gradient_balls_N_beta') is None:
            raise ValueError('complete node checkpoint missing true-gradient balls')
        if 'value_error_is_direct_distance_to_true_ball' not in (
                row.get('scope') or {}):
            raise ValueError('complete node checkpoint missing true-ball distance flag')


def load_node_checkpoint(directory, index, binding):
    row = V2.load_node_checkpoint(directory, index, binding)
    if row is None:
        return None
    require_node_matches_binding(row, binding, index)
    return row


def node_summary(row):
    summary = {
        'node': int(row['node']),
        'complete': bool(row.get('complete')),
        'obstruction': row.get('obstruction'),
        'budget_exceeded': bool(row.get('budget_exceeded')),
        'cpu_seconds': row.get('cpu_seconds'),
        'quadrature': [],
    }
    for count in GAUSS_COUNTS:
        match = next((item for item in row.get('quadrature', [])
                      if item.get('gauss_nodes_per_normal_panel') == count), None)
        if match is None:
            continue
        summary['quadrature'].append({
            'gauss_nodes_per_normal_panel': count,
            'phase_value_control_pass': bool(match.get('phase_value_control_pass')),
            'value_error_upper_N_beta': match.get('value_error_upper_N_beta'),
        })
    return summary


def failed_nodes_for_count(rows, count):
    failed = []
    compared = 0
    for row in rows:
        if not row.get('complete'):
            continue
        match = next((item for item in row.get('quadrature', [])
                      if item.get('gauss_nodes_per_normal_panel') == count), None)
        if match is None:
            continue
        compared += 1
        if not match.get('phase_value_control_pass'):
            failed.append(int(row['node']))
    return failed, compared


def all_declared_covered(rows, requested, declared_count, complete_run):
    declared = set(range(int(declared_count)))
    requested_set = {int(v) for v in requested}
    completed = {int(row['node']) for row in rows if row.get('complete')}
    return bool(
        complete_run
        and requested_set == declared
        and completed == declared
        and len(rows) == int(declared_count)
        and all(row.get('complete') for row in rows))


def nodal_component_pass(result):
    """True when every compared production gradient meets 1e-12 in N and beta."""
    attempts = result.get('quadrature_attempts') or {}
    if not attempts:
        return False
    for count in GAUSS_COUNTS:
        item = attempts.get(str(count))
        if item is None or not item.get('phase_value_control_pass'):
            return False
        if int(item.get('nodes_compared') or 0) <= 0:
            return False
        if int(item.get('failed_node_count') or 0) != 0:
            return False
    return True


def summarize(binding, rows, requested, declared_count, cpu, cap, complete_run,
              payloads=None):
    finished = [row for row in rows if row.get('complete')]
    attempts = {}
    for count in GAUSS_COUNTS:
        errors = []
        passed = []
        failed, compared = failed_nodes_for_count(rows, count)
        for row in finished:
            match = next((item for item in row.get('quadrature', [])
                          if item['gauss_nodes_per_normal_panel'] == count), None)
            if match is None:
                continue
            errors.append(match['value_error_upper_N_beta'])
            passed.append(match['phase_value_control_pass'])
        maximum = None
        if errors:
            with ctx.workprec(BITS):
                maximum = []
                for component in (0, 1):
                    value = arb(0)
                    for row in errors:
                        value = value.max(unpack_upper(row[component])).upper()
                    maximum.append(pack_upper(value))
        attempts[str(count)] = {
            'gauss_nodes_per_normal_panel': count,
            'nodes_compared': compared,
            'failed_node_indices': failed,
            'failed_node_count': len(failed),
            'phase_value_control_pass': bool(passed) and all(passed) and compared > 0,
            'value_error_upper_N_beta': maximum,
        }
    mean = None
    estimate = None
    maximum_cpu = None
    if rows:
        costs = [float(row['cpu_seconds']) for row in rows]
        mean = float(sum(costs) / len(costs))
        estimate = mean * int(declared_count)
        maximum_cpu = float(max(costs))
    obstructions = [row['obstruction'] for row in rows if row.get('obstruction')]
    covered = all_declared_covered(rows, requested, declared_count, complete_run)
    if not complete_run:
        status = (
            'OPEN: CPU budget exceeded after '
            f'{len(rows)}/{len(requested)} requested current-history nodes; '
            'resumable; physical local gate OPEN')
    elif obstructions:
        status = ('OPEN: directed enclosure obstruction '
                  f'{obstructions[0]}; physical local gate OPEN')
    else:
        gauss = []
        for count in GAUSS_COUNTS:
            item = attempts[str(count)]
            compared = item['nodes_compared']
            failed = item['failed_node_indices']
            if compared == 0:
                gauss.append(f'{count} not compared')
            elif item['phase_value_control_pass']:
                gauss.append(
                    f'{count} meets 1e-12 on {compared} of {compared} requested nodes')
            else:
                listed = ', '.join(str(v) for v in failed)
                gauss.append(
                    f'{count} misses 1e-12 on {len(failed)} of {compared} '
                    f'requested nodes (failed nodes: {listed})')
        status = (
            'OPEN: nodal phase-value on requested current-history nodes ('
            + '; '.join(gauss) + '); between-node None; physical local gate OPEN')
    result = {
        'schema': SCHEMA,
        'accountable_author': 'Douglas Ek',
        'status': status,
        'history_identity': binding['profile_identity'],
        'history_path': binding['history_path'],
        'node_count': binding['node_count'],
        'requested_node_indices': list(requested),
        'completed_node_indices': [row['node'] for row in rows],
        'rho_up': binding['rho_up'],
        'bits': BITS,
        'phase_component_tolerance': PHASE_COMPONENT_TOLERANCE,
        'target_radius_divisor': TARGET_RADIUS_DIVISOR,
        'target_radius_upper': binding['target_radius_upper'],
        'inherited_unit_angular_target_upper': binding[
            'inherited_unit_angular_target_upper'],
        'signed_angular_square_weights': binding['signed_angular_square_weights'],
        'exact_angular_square_weight': ANGULAR_SQUARE_WEIGHT,
        'gauss_nodes_per_normal_panel': list(GAUSS_COUNTS),
        'quadrature_attempts': attempts,
        'node_summaries': [node_summary(row) for row in rows],
        'coverage': {
            'declared_family_nodes': int(declared_count),
            'requested_nodes': [int(v) for v in requested],
            'completed_nodes': [row['node'] for row in rows],
            'all_declared_nodes_covered': covered,
            'mean_cpu_seconds_per_completed_node': mean,
            'max_cpu_seconds_per_completed_node': maximum_cpu,
            'estimated_cpu_seconds_for_declared_nodes': estimate,
        },
        'scope': {
            'physical_local_gate': 'OPEN',
            'all_declared_nodes_covered': covered,
            'between_node_phase_error_bound': None,
            'phase_tangent_error_bound': None,
            'finite_source_tail_bound': None,
            'source_and_field_error_bound': None,
            'value_error_is_direct_distance_to_true_ball': True,
            'new_field_or_source_evolutions': 0,
            'physical_optimizer_run': False,
            'metric_timestep': False,
            'physical_EXISTENCE_certificate': False,
        },
        'source_hashes': binding['source_hashes'],
        'binding': binding,
        'runtime': {'CPU_seconds': float(cpu), 'CPU_cap': float(cap),
                    'complete': bool(complete_run)},
    }
    if payloads is not None:
        result['payloads'] = payloads
        result['hash_descriptor'] = payloads['sha256']
    return result


def record_run(history, node_count, selection, node_indices, prefix, cpu_budget,
               allow_all):
    if not math.isfinite(cpu_budget) or cpu_budget <= 0:
        raise ValueError('CPU budget must be positive')
    if selection == 'all' and node_indices is None and not allow_all:
        raise ValueError(
            'refusing all declared nodes before a completed pilot; '
            'record selection=pilot first, then pass --allow-all')
    directory, output = paths(prefix)
    if output.exists():
        raise FileExistsError('current phase-accuracy record exists; use --check')
    ctx_data = load_context(history, node_count)
    indices = select_indices(node_count, selection, node_indices)
    binding = evaluation_binding(
        ctx_data['identity'], ctx_data['parent_path'], ctx_data['nodes'], indices,
        ctx_data['rho_up'], ctx_data['weights'], ctx_data['a_numeric'],
        ctx_data['hashes'])
    persist_binding(directory, binding)
    started = time.process_time()
    rows = []
    finished_all = False

    def exhausted(*_):
        raise TimeoutError('current phase enclosure exceeded its CPU budget')

    previous = signal.signal(signal.SIGPROF, exhausted)
    try:
        for index in indices:
            existing = load_node_checkpoint(directory, index, binding)
            if existing is not None:
                rows.append(existing)
                continue
            remaining = float(cpu_budget) - (time.process_time() - started)
            if remaining <= 0:
                break
            signal.setitimer(signal.ITIMER_PROF, remaining)
            try:
                row = compute_node(
                    ctx_data['family'], ctx_data['identity'],
                    ctx_data['nodes'][index], index, ctx_data['rho_up'],
                    ctx_data['weights'], ctx_data['a_numeric'],
                    ctx_data['a_interval'], ctx_data['hashes'],
                    cpu_limit=remaining)
            except TimeoutError:
                break
            finally:
                signal.setitimer(signal.ITIMER_PROF, 0)
            if row.get('budget_exceeded') and not row.get('complete'):
                break
            write_once(node_file(directory, index), row)
            rows.append(row)
        else:
            finished_all = len(rows) == len(indices)
    finally:
        signal.setitimer(signal.ITIMER_PROF, 0)
        signal.signal(signal.SIGPROF, previous)
    cpu = time.process_time() - started
    payloads = None
    if finished_all:
        payloads = artifact_descriptor(directory, indices)
    result = summarize(
        binding, rows, indices, node_count, cpu, cpu_budget, finished_all,
        payloads=payloads)
    if finished_all:
        write_once(output, result)
    return result


def replay_node(ctx_data, binding, index, cpu_limit):
    return compute_node(
        ctx_data['family'], ctx_data['identity'], ctx_data['nodes'][index],
        index, ctx_data['rho_up'], ctx_data['weights'], ctx_data['a_numeric'],
        ctx_data['a_interval'], ctx_data['hashes'], cpu_limit=cpu_limit)


def check_run(history, node_count, selection, node_indices, prefix, cpu_budget,
              allow_all):
    directory, output = paths(prefix)
    if not output.exists():
        raise FileNotFoundError('current phase-accuracy record missing; use --record')
    saved = json.loads(output.read_text())
    ctx_data = load_context(history, node_count)
    indices = select_indices(node_count, selection, node_indices)
    binding = evaluation_binding(
        ctx_data['identity'], ctx_data['parent_path'], ctx_data['nodes'], indices,
        ctx_data['rho_up'], ctx_data['weights'], ctx_data['a_numeric'],
        ctx_data['hashes'])
    persisted = json.loads(binding_file(directory).read_text())
    if persisted != binding or saved['binding'] != binding:
        raise ValueError('checkpoint belongs to another source/history/numerics')
    if saved.get('schema') != SCHEMA:
        raise ValueError('record schema is not the v3 phase-accuracy owner')
    authenticate_payloads(saved['payloads'], directory, indices)
    started = time.process_time()
    rows = []
    remaining = float(cpu_budget)
    for index in indices:
        path = node_file(directory, index)
        previous = json.loads(path.read_text())
        require_node_matches_binding(previous, binding, index)
        row = replay_node(ctx_data, binding, index, remaining)
        remaining = float(cpu_budget) - (time.process_time() - started)
        if canonical(row) != canonical(previous):
            raise ValueError(f'current phase-value replay differs at node {index}')
        rows.append(row)
        if remaining <= 0 and index != indices[-1]:
            raise TimeoutError('current phase-accuracy --check exceeded its CPU budget')
    cpu = time.process_time() - started
    payloads = artifact_descriptor(directory, indices)
    result = summarize(
        binding, rows, indices, node_count, cpu, cpu_budget, True, payloads=payloads)
    if canonical(result) != canonical(saved):
        raise ValueError('current phase-value replay differs')
    return result


def report(result):
    attempts = result['quadrature_attempts']
    compact = {}
    for key, item in attempts.items():
        compact[key] = {
            'failed_node_count': item['failed_node_count'],
            'failed_node_indices': item['failed_node_indices'],
            'nodes_compared': item['nodes_compared'],
            'phase_value_control_pass': item['phase_value_control_pass'],
            'value_error_upper_N_beta': item['value_error_upper_N_beta'],
        }
    return {
        'status': result['status'],
        'schema': result['schema'],
        'history_identity': result['history_identity'],
        'requested_node_indices': result['requested_node_indices'],
        'completed_node_indices': result['completed_node_indices'],
        'target_radius_divisor': result.get('target_radius_divisor'),
        'quadrature_attempts': compact,
        'coverage': result['coverage'],
        'hash_descriptor': result.get('hash_descriptor'),
        'CPU_seconds': result['runtime']['CPU_seconds'],
        'complete': result['runtime']['complete'],
        'physical_local_gate': result['scope']['physical_local_gate'],
        'between_node_phase_error_bound': result['scope']['between_node_phase_error_bound'],
        'all_declared_nodes_covered': result['scope']['all_declared_nodes_covered'],
        'nodal_component_pass': nodal_component_pass(result),
    }


def main(argv=None):
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--record', action='store_true')
    mode.add_argument('--check', action='store_true')
    parser.add_argument('--history', default=DEFAULT_HISTORY)
    parser.add_argument('--nodes', type=int, required=True, choices=(129, 257))
    parser.add_argument('--selection', choices=('pilot', 'all'), default='pilot')
    parser.add_argument('--node-indices', default=None)
    parser.add_argument('--output-prefix', default='pilot')
    parser.add_argument('--cpu-budget', type=float, default=60.0)
    parser.add_argument('--allow-all', action='store_true')
    args = parser.parse_args(argv)
    indices = parse_indices(args.node_indices)
    if args.record:
        result = record_run(
            args.history, args.nodes, args.selection, indices, args.output_prefix,
            args.cpu_budget, args.allow_all)
    else:
        result = check_run(
            args.history, args.nodes, args.selection, indices, args.output_prefix,
            args.cpu_budget, args.allow_all)
    print(json.dumps(report(result), indent=2, sort_keys=True))
    return result


if __name__ == '__main__':
    main()
