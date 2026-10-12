#!/usr/bin/env python3
"""CURRENT-HISTORY nodal phase-value enclosure. No Dirac or source evolution.

Successor of scripts/derive_nsc_ks_coupled_phase_accuracy.py for the live
n=32 history. Each node proof binds that history, the actual collocation
point, the original angular ledger, the numeric phase settings and the
code. Production 192/384 unit-angular gradients are compared to the
directed true-gradient ball, not only a midpoint/radius. Between-node
phase remainder stays None. The physical local gate stays OPEN.
"""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'

import argparse
from hashlib import sha256
import json
from pathlib import Path
import re
import signal
import sys
import time

import numpy as np
from flint import arb, ctx

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))

from recursive_horizons.nsc_dirac_source_phase import formal_source_phase_coefficient
from recursive_horizons.nsc_dirac_source_phase_bound import (
    ANGULAR_SQUARE_WEIGHT, BITS as ENCLOSURE_BITS, SOURCE_SIGNS,
    directed_source_phase_z_enclosure, pack_real_ball, pack_upper, unpack_upper)
from recursive_horizons.nsc_ks_ball_geometry import background_series
from recursive_horizons.nsc_ks_cutoff_bridge_evaluator import (
    signed_family_angular_square_weights, source_edge_from_phase)
from recursive_horizons.nsc_ks_profile_identity import profile_identity
from recursive_horizons.nsc_ks_retained_upstream_archive import RetainedUpstreamArchive
from recursive_horizons.nsc_ks_value_evaluator import (
    TARGET_COUNTS, one_direction_metric, target_nodes)
from recursive_horizons.nsc_local_incoming_constraints import _channel_record
from recursive_horizons.nsc_local_incoming_family import LocalIncomingFamily

SCHEMA = 'NSC-KS-CURRENT-PHASE-ACCURACY-v2'
DEFAULT_HISTORY = 'results/development/nsc-ks-gate-history-lm-broyden.json'
BRANCH = 'results/development/nsc-incoming-surface-regular-branch.json'
INVENTORY = 'results/development/nsc-ks-source-inventory.json'
RETAINED = 'results/development/nsc-ks-retained-response-sum.json'
RHO_UP_LEDGER = 'results/development/nsc-ks-cutoff-bridge-control.json'
ORIGINAL_RHO_UP = 1.0300000000000002
EXPECTED_WEIGHT = float(ANGULAR_SQUARE_WEIGHT)
EXPECTED_SIGNED_FAMILIES = 120
BITS = ENCLOSURE_BITS
PHASE_COMPONENT_TOLERANCE = 1e-12
GAUSS_COUNTS = (192, 384)
PREFIX_RE = re.compile(r'^[a-z0-9-]+$')
OWNERS = (
    'scripts/derive_nsc_ks_current_phase_accuracy_v2.py',
    'docs/nsc-ks-current-phase-accuracy-v2.md',
    'src/recursive_horizons/nsc_dirac_source_phase_bound.py',
    'src/recursive_horizons/nsc_dirac_source_phase.py',
    'src/recursive_horizons/nsc_ks_cutoff_bridge_evaluator.py',
    'src/recursive_horizons/nsc_ks_value_evaluator.py',
    'src/recursive_horizons/nsc_ks_ball_geometry.py',
    'src/recursive_horizons/nsc_ks_ball_operator.py',
    'src/recursive_horizons/nsc_ks_profile_identity.py',
    'src/recursive_horizons/nsc_local_incoming_family.py',
)
RUNTIME_KEYS = frozenset({
    'cpu_seconds', 'CPU_seconds', 'runtime',
    'mean_cpu_seconds_per_completed_node',
    'estimated_cpu_seconds_for_declared_nodes',
})


def digest(path):
    path = Path(path)
    return sha256((path if path.is_absolute() else ROOT / path).read_bytes()).hexdigest()


def write_once(path, data):
    """Publish complete bytes without replacing an existing scientific file."""
    raw = data if isinstance(data, bytes) else (
        json.dumps(data, sort_keys=True, indent=2, allow_nan=False) + '\n').encode()
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f'.{os.getpid()}.tmp')
    created = False
    try:
        with temporary.open('xb') as stream:
            created = True
            stream.write(raw)
            stream.flush()
            os.fsync(stream.fileno())
        os.link(temporary, path)
    finally:
        if created and temporary.exists():
            temporary.unlink()


def without_runtime(value):
    if isinstance(value, dict):
        return {key: without_runtime(item) for key, item in value.items()
                if key not in RUNTIME_KEYS}
    if isinstance(value, list):
        return [without_runtime(item) for item in value]
    return value


def canonical(value):
    """JSON round-trip so --check compares scientific bytes, not CPU or live types."""
    return json.loads(json.dumps(without_runtime(value), sort_keys=True, allow_nan=False))


def require_prefix(prefix):
    if not isinstance(prefix, str) or PREFIX_RE.fullmatch(prefix) is None:
        raise ValueError('output prefix must match [a-z0-9-]+')
    return prefix


def history_path(path=None):
    chosen = ROOT / DEFAULT_HISTORY if path is None else Path(path)
    if not chosen.is_absolute():
        chosen = ROOT / chosen
    resolved = chosen.resolve()
    development = (ROOT / 'results/development').resolve()
    if resolved.parent != development or resolved.suffix != '.json':
        raise ValueError('history must be one json register in results/development')
    return resolved


def paths(prefix):
    safe = require_prefix(prefix)
    artifacts = (ROOT / 'results/development/artifacts').resolve()
    development = (ROOT / 'results/development').resolve()
    directory = (artifacts / f'nsc-ks-current-phase-accuracy-v2-{safe}').resolve()
    output = (development / f'nsc-ks-current-phase-accuracy-v2-{safe}.json').resolve()
    if directory.parent != artifacts or output.parent != development:
        raise ValueError('prefix escapes the development directories')
    return directory, output


def node_file(directory, index):
    return Path(directory) / f'node-{int(index):03d}.json'


def binding_file(directory):
    return Path(directory) / 'binding.json'


def load_history(path=None):
    parent_path = history_path(path)
    parent = json.loads(parent_path.read_text())
    coefficients = np.asarray(parent['history']['coefficients'], float)
    if coefficients.shape != (2, 32):
        raise ValueError('current history must be n=32 w/U')
    family = LocalIncomingFamily(coefficients)
    identity = profile_identity(family, include_normal_window=True)
    if identity != parent['profile_identity']:
        raise ValueError('history identity does not match its coefficients')
    if family.radius_lower_bound() <= 0:
        raise ValueError('positive-radius current history required')
    recorded_lower = parent['history'].get('radius_lower_bound')
    if recorded_lower is not None and float(recorded_lower) != family.radius_lower_bound():
        raise ValueError('recorded radius lower bound does not match the coefficients')
    return parent_path, parent, family, identity


def original_archive_rho_up():
    """Bind the original upstream slice. This does not restore source columns."""
    ledger = json.loads((ROOT / RHO_UP_LEDGER).read_text())
    rho_up = float(ledger['rho_up'])
    if rho_up != ORIGINAL_RHO_UP:
        raise ValueError('original archive rho_up is not the bound binary64 slice')
    if not 1.03 <= rho_up <= 33 / 32:
        raise ValueError('original rho_up must lie in the owned preparation slab')
    return rho_up


def original_signed_family_records(archive):
    """Original 120 signed families from the inventory ledger. No family_entries."""
    records = {}
    channels = archive.meta['channels']
    for group, sign in archive.family_keys:
        channel = channels[int(group)]
        positive = _channel_record('_', {'_': {**channel, 'angular_sign': int(sign)}})
        negative = _channel_record('_', {'_': {**channel, 'angular_sign': -int(sign)}})
        for energy_sign, item in ((1, positive), (-1, negative)):
            signed_id = f'{item["group"]}_{item["angular_sign"]}_E{energy_sign:+d}'
            if signed_id in records:
                raise ValueError('duplicate original signed family')
            records[signed_id] = {
                'group': item['group'],
                'angular_sign': item['angular_sign'],
                'energy_sign': energy_sign,
                'compact_mass': item['compact_mass'],
                'angular_eigenvalue': item['angular_eigenvalue'],
                'signed_angular': item['angular_sign'] * item['angular_eigenvalue'],
                'copy_count': item['copy_count'],
                'degeneracy': item['degeneracy'],
                'multiplicity': item['multiplicity'],
                'ell0_pure_radius_identity': item['ell0_pure_radius_identity'],
            }
    if len(records) != EXPECTED_SIGNED_FAMILIES:
        raise ValueError('one hundred twenty explicit signed families required')
    return records


def production_angular_weights(archive):
    records = original_signed_family_records(archive)
    weights = signed_family_angular_square_weights(records)
    if not np.allclose(weights, EXPECTED_WEIGHT, rtol=0, atol=1e-6):
        raise ValueError('signed angular-square weights must be 107952 for each energy sign')
    return weights, records


def axial_scale(branch=None):
    record = branch if branch is not None else json.loads((ROOT / BRANCH).read_text())
    interval = record['branch']['intrinsic_intervals']['a']
    if len(interval) != 2 or not np.isfinite(interval).all() or not interval[0] <= interval[1]:
        raise ValueError('certified incoming axial-scale interval required')
    numeric = float(np.mean(interval))
    if not np.isfinite(numeric) or numeric <= 0:
        raise ValueError('positive production axial scale required')
    return numeric, (float(interval[0]), float(interval[1]))


def collocation(family, count):
    if count not in TARGET_COUNTS or count not in (129, 257):
        raise ValueError('current-history phase nodes must be 129 or 257')
    nodes = np.asarray(target_nodes(family, count), float)
    if len(nodes) != count or not np.isfinite(nodes).all():
        raise ValueError('finite current-family collocation required')
    lo, hi = family.interval
    if not (nodes[0] >= lo and nodes[-1] <= hi):
        raise ValueError('collocation nodes must lie on the declared interval I')
    return nodes


def select_indices(count, selection='pilot', node_indices=None):
    if node_indices is not None:
        indices = tuple(int(v) for v in node_indices)
        if not indices or len(set(indices)) != len(indices):
            raise ValueError('distinct node indices required')
        if any(index < 0 or index >= count for index in indices):
            raise ValueError('node index outside the declared collocation')
        return indices
    if selection == 'pilot':
        return (0, count // 2, count - 1)
    if selection == 'all':
        return tuple(range(count))
    raise ValueError("selection must be 'pilot' or 'all'")


def true_edge_gradient(plus, minus, a_true, weight=ANGULAR_SQUARE_WEIGHT):
    """Exact original-ledger contraction of unit-angular f_z balls."""
    if int(weight) != ANGULAR_SQUARE_WEIGHT:
        raise ValueError('true factor uses the exact original angular-square weight')
    weight = arb(int(weight))
    two_pi = 2 * arb.pi()
    plus_term, minus_term = weight * plus, weight * minus
    return ((plus_term - minus_term) / (a_true * two_pi),
            -(plus_term + minus_term) / two_pi)


def directed_distance_upper(value, ball):
    """Distance of a production number to the whole true ball, not only its radius."""
    return abs(arb(float(value)) - ball).upper()


def numeric_contraction(plus, minus, weights, axial):
    plus_term = arb(float(weights[0])) * arb(float(plus))
    minus_term = arb(float(weights[1])) * arb(float(minus))
    two_pi = 2 * arb.pi()
    return ((plus_term - minus_term) / (axial * two_pi),
            -(plus_term + minus_term) / two_pi)


def contraction_report(phase_fz, weights, a_numeric, a_interval, a_true, true_gradient,
                       production):
    plus, minus = (float(phase_fz[0]), float(phase_fz[1]))
    a_ball = arb(float(a_interval[0])).union(arb(float(a_interval[1])))
    exact_weight = (ANGULAR_SQUARE_WEIGHT, ANGULAR_SQUARE_WEIGHT)
    phase_only = numeric_contraction(plus, minus, exact_weight, a_true)
    numeric_a = numeric_contraction(plus, minus, weights, arb(float(a_numeric)))
    interval_a = numeric_contraction(plus, minus, weights, a_ball)
    produced = (arb(float(production[0])), arb(float(production[1])))
    pack = lambda pair, truth: [pack_upper(abs(pair[k] - truth[k]).upper()) for k in (0, 1)]
    return {
        'phase_only_to_true_ball_upper_N_beta': pack(phase_only, true_gradient),
        'numeric_a_weight_to_true_ball_upper_N_beta': pack(numeric_a, true_gradient),
        'a_interval_to_true_ball_upper_N_beta': pack(interval_a, true_gradient),
        'float_contraction_rounding_upper_N_beta': [
            pack_upper(abs(produced[k] - numeric_a[k]).upper()) for k in (0, 1)],
        'production_to_true_ball_upper_N_beta': pack(produced, true_gradient),
        'value_error_is_direct_distance_to_true_ball': True,
    }


def one_direction(family):
    metric = one_direction_metric(family)
    if len(metric.directions) != 1 or len(metric.amplitudes) != 1:
        raise ValueError('unit-angular physical one-direction metric required')
    if len(family.metric().directions) != 2 * len(family.coefficients[0]):
        raise ValueError('basis metric still has one direction per coefficient')
    w, u = family.functions
    direction = metric.directions[0]
    if (direction.w.coefficients != w.coefficients or direction.U.coefficients != u.coefficients
            or float(direction.w.center) != float(w.center)
            or float(direction.U.center) != float(u.center)):
        raise ValueError('one-direction metric must carry the actual w/U profiles')
    if (float(direction.inner_radius) != family.normal_inner
            or float(direction.outer_radius) != family.normal_outer):
        raise ValueError('one-direction metric must keep the owned normal window')
    return metric


def source_hashes(extra=None):
    hashes = {name: digest(name) for name in OWNERS}
    hashes.update({name: digest(name) for name in (BRANCH, INVENTORY, RETAINED, RHO_UP_LEDGER)})
    if extra:
        hashes.update(extra)
    return hashes


def evaluation_binding(identity, history, nodes, indices, rho_up, weights, a_numeric,
                       hashes):
    return {
        'schema': SCHEMA,
        'profile_identity': identity,
        'history_path': str(Path(history).relative_to(ROOT)),
        'history_sha256': digest(history),
        'node_count': int(len(nodes)),
        'requested_node_indices': [int(v) for v in indices],
        'target_nodes': [float(nodes[i]) for i in indices],
        'target_node_hex': [float(nodes[i]).hex() for i in indices],
        'rho_up': float(rho_up),
        'rho_up_hex': float(rho_up).hex(),
        'signed_angular_square_weights': [float(v) for v in weights],
        'exact_angular_square_weight': ANGULAR_SQUARE_WEIGHT,
        'production_axial_scale': float(a_numeric),
        'gauss_nodes_per_normal_panel': list(GAUSS_COUNTS),
        'bits': BITS,
        'phase_component_tolerance': PHASE_COMPONENT_TOLERANCE,
        'one_direction_metric': True,
        'new_field_or_source_evolutions': 0,
        'source_hashes': hashes,
    }


def pack_quadrature(count, phase, gradient, weights, a_numeric, a_interval, a_true,
                    true_gradient):
    row = np.asarray(gradient, float)
    if row.shape != (2,) or not np.isfinite(row).all():
        raise ValueError('finite production (N, beta) edge required')
    report = contraction_report(
        np.asarray(phase.f_z)[:, 0], weights, a_numeric, a_interval, a_true,
        true_gradient, row)
    errors = report['production_to_true_ball_upper_N_beta']
    passed = all(unpack_upper(item) <= arb(1) / 10**12 for item in errors)
    return {
        'gauss_nodes_per_normal_panel': int(count),
        'phase_binding': phase.binding,
        'numeric_f_z_unit_angular': [float(v) for v in np.asarray(phase.f_z)[:, 0]],
        'numeric_f_z_hex': [float(v).hex() for v in np.asarray(phase.f_z)[:, 0]],
        'edge_gradient': row.tolist(),
        'value_error_upper_N_beta': errors,
        'phase_value_control_pass': bool(passed),
        'contraction': report,
        'one_direction_amplitudes': 1,
        'tangent_directions': int(np.asarray(phase.delta_f_z).shape[1]),
    }


def compute_node(family, identity, point, index, rho_up, weights, a_numeric, a_interval,
                 hashes, *, cpu_limit):
    started = time.process_time()
    metric = one_direction(family)
    if profile_identity(family, include_normal_window=True) != identity:
        raise ValueError('history identity changed before the node proof')
    z = np.asarray([float(point)], float)
    with ctx.workprec(BITS):
        a_true = background_series(arb(1), 0).a[0]
        remaining = float(cpu_limit) - (time.process_time() - started)
        if remaining <= 0:
            return {'node': int(index), 'z': float(point), 'budget_exceeded': True,
                    'obstruction': 'cpu_budget', 'complete': False,
                    'cpu_seconds': time.process_time() - started}
        enclosure = directed_source_phase_z_enclosure(
            metric, float(point), rho_up=rho_up, bits=BITS, cpu_limit=remaining)
        packed_true = None
        true_gradient = None
        quadrature = []
        if enclosure.obstruction == 'cpu_budget':
            return {'node': int(index), 'z': float(point), 'budget_exceeded': True,
                    'obstruction': 'cpu_budget', 'complete': False,
                    'cpu_seconds': time.process_time() - started}
        if enclosure.obstruction is None:
            if len(enclosure.balls) != 2 or enclosure.signs != SOURCE_SIGNS:
                raise ArithmeticError('both original source signs required')
            plus, minus = enclosure.balls
            true_gradient = true_edge_gradient(plus, minus, a_true)
            packed_true = [pack_real_ball(v) for v in true_gradient]
            for count in GAUSS_COUNTS:
                remaining = float(cpu_limit) - (time.process_time() - started)
                if remaining <= 0:
                    return {'node': int(index), 'z': float(point),
                            'budget_exceeded': True, 'obstruction': 'cpu_budget',
                            'complete': False,
                            'cpu_seconds': time.process_time() - started}
                phase = formal_source_phase_coefficient(
                    metric, z, angular=1.0, rho_up=float(rho_up), gauss_nodes=count)
                gradient, tangent = source_edge_from_phase(phase, weights, a_numeric)
                if int(np.asarray(tangent).shape[0]) != 1:
                    raise ValueError('one-direction edge tangent has the wrong count')
                quadrature.append(pack_quadrature(
                    count, phase, np.asarray(gradient, float)[0], weights, a_numeric,
                    a_interval, a_true, true_gradient))
        row = {
            'node': int(index),
            'z': float(point),
            'z_hex': float(point).hex(),
            'history_identity': identity,
            'enclosed_profile_identity': enclosure.profile_identity,
            'one_direction_metric': True,
            'rho_up': float(rho_up),
            'signed_angular_square_weights': [float(v) for v in weights],
            'exact_angular_square_weight': ANGULAR_SQUARE_WEIGHT,
            'bits': BITS,
            'gauss_nodes_per_normal_panel': list(GAUSS_COUNTS),
            'phase_component_tolerance': PHASE_COMPONENT_TOLERANCE,
            'enclosure': {
                'status': enclosure.status,
                'obstruction': enclosure.obstruction,
                'f_z_unit_angular': None if enclosure.obstruction else [
                    pack_real_ball(v) for v in enclosure.balls],
                'cells': list(enclosure.cells),
                'max_depth': list(enclosure.max_depth),
                'splits': list(enclosure.splits),
                'cover_cells': enclosure.cover_cells,
                'dyadic_scale': enclosure.dyadic_scale,
                'taylor_degree': enclosure.taylor_degree,
                'remainder_order': enclosure.remainder_order,
            },
            'a_true': pack_real_ball(a_true),
            'production_axial_scale': float(a_numeric),
            'production_axial_interval': [float(v) for v in a_interval],
            'edge_gradient_balls_N_beta': packed_true,
            'quadrature': quadrature,
            'source_hashes': hashes,
            'complete': enclosure.obstruction is None,
            'budget_exceeded': False,
            'obstruction': enclosure.obstruction,
            'scope': {
                'between_node_phase_error_bound': None,
                'phase_tangent_error_bound': None,
                'finite_source_tail_bound': None,
                'source_and_field_error_bound': None,
                'physical_local_gate': 'OPEN',
                'new_field_or_source_evolutions': 0,
                'physical_optimizer_run': False,
                'value_error_is_direct_distance_to_true_ball': True,
            },
            'cpu_seconds': time.process_time() - started,
        }
    return row


def summarize(binding, rows, requested, declared_count, cpu, cap, complete_run):
    finished = [row for row in rows if row.get('complete')]
    attempts = {}
    for count in GAUSS_COUNTS:
        errors = []
        passed = []
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
            'nodes_compared': len(errors),
            'phase_value_control_pass': bool(passed) and all(passed),
            'value_error_upper_N_beta': maximum,
        }
    mean = None
    estimate = None
    if rows:
        mean = float(sum(float(row['cpu_seconds']) for row in rows) / len(rows))
        estimate = mean * int(declared_count)
    obstructions = [row['obstruction'] for row in rows if row.get('obstruction')]
    if not complete_run:
        status = ('OPEN: CPU budget exceeded after '
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
            if compared == 0:
                gauss.append(f'{count} not compared')
            elif item['phase_value_control_pass']:
                gauss.append(f'{count} meets 1e-12 on {compared} nodes')
            else:
                gauss.append(f'{count} misses 1e-12 on {compared} nodes')
        status = (
            'OPEN: nodal phase-value on requested current-history nodes ('
            + '; '.join(gauss) + '); between-node None; physical local gate OPEN')
    return {
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
        'signed_angular_square_weights': binding['signed_angular_square_weights'],
        'exact_angular_square_weight': ANGULAR_SQUARE_WEIGHT,
        'gauss_nodes_per_normal_panel': list(GAUSS_COUNTS),
        'quadrature_attempts': attempts,
        'rows': rows,
        'coverage': {
            'declared_family_nodes': int(declared_count),
            'requested_nodes': [int(v) for v in requested],
            'completed_nodes': [row['node'] for row in rows],
            'all_declared_nodes_covered': (
                complete_run and len(requested) == int(declared_count)
                and all(row.get('complete') for row in rows)),
            'mean_cpu_seconds_per_completed_node': mean,
            'estimated_cpu_seconds_for_declared_nodes': estimate,
        },
        'scope': {
            'physical_local_gate': 'OPEN',
            'all_declared_nodes_covered': (
                complete_run and len(requested) == int(declared_count)),
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


def load_context(history, node_count):
    parent_path, parent, family, identity = load_history(history)
    nodes = collocation(family, node_count)
    rho_up = original_archive_rho_up()
    archive = RetainedUpstreamArchive(ROOT)
    weights, records = production_angular_weights(archive)
    a_numeric, a_interval = axial_scale()
    hashes = source_hashes({
        str(parent_path.relative_to(ROOT)): digest(parent_path),
        **archive.input_hashes,
    })
    return {
        'parent_path': parent_path, 'parent': parent, 'family': family,
        'identity': identity, 'nodes': nodes, 'rho_up': rho_up,
        'weights': weights, 'records': records, 'a_numeric': a_numeric,
        'a_interval': a_interval, 'hashes': hashes, 'archive': archive,
    }


def persist_binding(directory, binding):
    path = binding_file(directory)
    if path.exists():
        if json.loads(path.read_text()) != binding:
            raise ValueError('checkpoint belongs to another source/history/numerics')
        return path
    write_once(path, binding)
    return path


def load_node_checkpoint(directory, index, binding):
    path = node_file(directory, index)
    if not path.exists():
        return None
    row = json.loads(path.read_text())
    if row.get('budget_exceeded'):
        raise ValueError('incomplete budget-exceeded node file must not be reused')
    if row.get('history_identity') != binding['profile_identity']:
        raise ValueError('node checkpoint history identity changed')
    if row.get('node') != int(index):
        raise ValueError('node checkpoint index changed')
    expected_hex = None
    for requested, hex_value in zip(binding['requested_node_indices'], binding['target_node_hex']):
        if requested == int(index):
            expected_hex = hex_value
            break
    if expected_hex is None or row.get('z_hex') != expected_hex:
        raise ValueError('node checkpoint z changed')
    if float(row['rho_up']).hex() != binding['rho_up_hex']:
        raise ValueError('node checkpoint rho_up changed')
    if row.get('source_hashes') != binding['source_hashes']:
        raise ValueError('node checkpoint code or ledger hashes changed')
    return row


def record_run(history, node_count, selection, node_indices, prefix, cpu_budget,
               allow_all):
    if not np.isfinite(cpu_budget) or cpu_budget <= 0:
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
    result = summarize(
        binding, rows, indices, node_count, cpu, cpu_budget, finished_all)
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
    started = time.process_time()
    rows = []
    remaining = float(cpu_budget)
    for index in indices:
        path = node_file(directory, index)
        previous = json.loads(path.read_text())
        row = replay_node(ctx_data, binding, index, remaining)
        remaining = float(cpu_budget) - (time.process_time() - started)
        if canonical(row) != canonical(previous):
            raise ValueError(f'current phase-value replay differs at node {index}')
        rows.append(row)
        if remaining <= 0 and index != indices[-1]:
            raise TimeoutError('current phase-accuracy --check exceeded its CPU budget')
    cpu = time.process_time() - started
    result = summarize(
        binding, rows, indices, node_count, cpu, cpu_budget, True)
    if canonical(result) != canonical(saved):
        raise ValueError('current phase-value replay differs')
    return result


def parse_indices(text):
    if text is None:
        return None
    return tuple(int(part) for part in text.split(',') if part != '')


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
    parser.add_argument('--cpu-budget', type=float, default=120.0)
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
    print(json.dumps({
        'status': result['status'],
        'history_identity': result['history_identity'],
        'requested_node_indices': result['requested_node_indices'],
        'completed_node_indices': result['completed_node_indices'],
        'quadrature_attempts': result['quadrature_attempts'],
        'coverage': result['coverage'],
        'CPU_seconds': result['runtime']['CPU_seconds'],
        'complete': result['runtime']['complete'],
        'physical_local_gate': result['scope']['physical_local_gate'],
    }, indent=2, sort_keys=True))
    return result


if __name__ == '__main__':
    main()
