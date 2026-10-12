#!/usr/bin/env python3
"""Matter columns of Chebyshev modes 32..63 at the iterate6 history."""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
import argparse
from hashlib import sha256
import json
from pathlib import Path
import signal
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
sys.path.insert(0, str(ROOT/'scripts'))
import derive_nsc_ks_coupled_newton_n16_damped as D
import derive_nsc_ks_coupled_retained_control as R
import recursive_horizons.nsc_ks_history_evaluator as H
from recursive_horizons.nsc_compatible_history_geometry import (
    CompatibleIncomingMetric, CompatibleRadiusDirection)
from recursive_horizons.nsc_dirac_source_phase import formal_source_phase_coefficient
from recursive_horizons.nsc_evolved_incoming_constraints import (
    compatible_history_slots, surface_geometry_response)
from recursive_horizons.nsc_ks_cutoff_bridge_evaluator import (
    signed_family_angular_square_weights, source_edge_from_phase)
from recursive_horizons.nsc_ks_energy_propagator import (
    KSEnergyPropagator, evolve_energy_propagator)
from recursive_horizons.nsc_ks_local_constraints import _family_matter
from recursive_horizons.nsc_ks_signed_state import negative_angular_partner
from recursive_horizons.nsc_ks_source_envelope import (
    KSEnvelopeBinding, usual_axial_support, _sample_axial_profiles)
from recursive_horizons.nsc_local_incoming_constraints import _channel_record
from recursive_horizons.nsc_local_incoming_family import (
    LocalAxialFunction, LocalIncomingFamily)
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes

DIRECTORY = ROOT/'results/development/artifacts/nsc-ks-coupled-newton-n64-high-modes'
PARENT_ARTIFACTS = ROOT/'results/development/artifacts/nsc-ks-coupled-newton-n32-truncated-iterate6'
PARENT = ROOT/'results/development/nsc-ks-coupled-newton-n32-truncated-iterate6.json'
OUTPUT = ROOT/'results/development/nsc-ks-coupled-newton-n64-high-modes.json'
HIGH = 64
FROZEN_PLUS = 65
FAMILY_CAP = 180.0
EXISTENCE = 3e-11
OWNERS = (
    'scripts/derive_nsc_ks_coupled_newton_n64_high_modes.py',
    'docs/nsc-ks-coupled-newton-n64-high-modes.md',
)


def digest(path):
    path = Path(path)
    return sha256((path if path.is_absolute() else ROOT/path).read_bytes()).hexdigest()


def family_paths(key):
    group = int(key[0])
    sign = int(key[1])
    if group < 0 or group > 99 or sign not in (-1, 1):
        raise ValueError('family key is not a retained angular family')
    stem = f'family-{group:02d}-{sign:+d}'
    return DIRECTORY/(stem+'.json'), DIRECTORY/(stem+'.npz')


def parent_family_paths(key):
    stem = f'family-{key[0]:02d}-{key[1]:+d}'
    return PARENT_ARTIFACTS/(stem+'.json'), PARENT_ARTIFACTS/(stem+'.npz')


def next_metric(family):
    coeff = np.asarray(family.coefficients, float)
    center, ai, ao = family.center, family.axial_inner, family.axial_outer
    ni, no = family.normal_inner, family.normal_outer
    zero = LocalAxialFunction((0.,), center, ai, ao)
    frozen = CompatibleRadiusDirection(
        LocalAxialFunction(tuple(coeff[0]), center, ai, ao),
        LocalAxialFunction(tuple(coeff[1]), center, ai, ao), ni, no)
    eye = np.eye(64)
    basis = [LocalAxialFunction(tuple(eye[i]), center, ai, ao) for i in range(32, 64)]
    high_w = tuple(CompatibleRadiusDirection(basis[i], zero, ni, no) for i in range(32))
    high_u = tuple(CompatibleRadiusDirection(zero, basis[i], ni, no) for i in range(32))
    return CompatibleIncomingMetric((1.,) + (0.,) * HIGH, (frozen,) + high_w + high_u)


def context():
    base = R.context()
    parent = json.loads(PARENT.read_text())
    family = LocalIncomingFamily(np.asarray(parent['history']['coefficients'], float))
    if family.radius_lower_bound() <= 0:
        raise ValueError('iterate6 history must keep a positive radius')
    metric = next_metric(family)
    solve = family.collocation_nodes(16)
    verify = family.collocation_nodes(33)
    target = H.fixed_target_union(solve, verify)
    if target.tolist() != parent['z']:
        raise ValueError('n=64 modes must keep the iterate6 nodes')
    return {
        **base, 'family': family, 'metric': metric, 'target': target,
        'identity': parent['profile_identity'], 'parent': parent,
        'hashes': {**base['archive'].input_hashes, **{p: digest(p) for p in OWNERS},
                   str(PARENT.relative_to(ROOT)): digest(PARENT)},
    }


def evaluate_family(ctx, key, saved=None, allow_new=True):
    archived = ctx['archived'][key]
    interval = R.interpolation_interval(archived)
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
    first = next(batch for _, batch, _ in group['applies'])
    metric = ctx['metric']
    if saved is None:
        if not allow_new:
            raise RuntimeError('replay must not call Dirac evolution')
        operator = evolve_energy_propagator(
            interval, R.DEGREE, metric, ctx['grid'], ctx['target'],
            first.mass, first.angular, first.rho_up,
            axial_support=usual_axial_support(), tangents='all', **R.OPTIONS)
    else:
        arrays, meta = saved
        w, u = _sample_axial_profiles(metric.directions, ctx['grid'])
        binding = KSEnvelopeBinding(
            ctx['grid'], w, u, tuple(metric.amplitudes),
            tuple(d.inner_radius for d in metric.directions),
            tuple(d.outer_radius for d in metric.directions),
            usual_axial_support(), first.rho_up, 1., **R.OPTIONS)
        operator = KSEnergyPropagator(
            interval, *(np.asarray(arrays['operator/'+name]) for name in R.KEYS),
            first.mass, first.angular, first.rho_up, binding, meta['operator']['diagnostics'])
        if operator.digest != meta['operator']['digest']:
            raise ValueError('saved n=64 operator digest changed')
    if operator.tangent.shape[1] != FROZEN_PLUS:
        raise ValueError('frozen-plus-next-mode tangents required')
    corrections, tangents, family_records, states = {}, {}, {}, {}

    def record_signed(signed_id, matter, channel, energy_sign, batch_ids):
        if signed_id in corrections:
            corrections[signed_id] = corrections[signed_id] + matter['correction']
            tangents[signed_id] = tangents[signed_id] + matter['tangent'][1:]
            family_records[signed_id]['batch_ids'].extend(batch_ids)
            return
        corrections[signed_id] = matter['correction']
        tangents[signed_id] = matter['tangent'][1:]
        family_records[signed_id] = {
            'group': channel['group'], 'angular_sign': channel['angular_sign'],
            'energy_sign': energy_sign, 'compact_mass': channel['compact_mass'],
            'angular_eigenvalue': channel['angular_eigenvalue'],
            'signed_angular': matter['signed_angular'],
            'copy_count': channel['copy_count'], 'degeneracy': channel['degeneracy'],
            'n_actual_angular_signs': channel['n_actual_angular_signs'],
            'multiplicity': matter['multiplicity'],
            'ell0_pure_radius_identity': matter['ell0_pure_radius_identity'],
            'batch_ids': list(batch_ids),
        }

    for batch_id, batch, channel in group['applies']:
        prepared = operator.apply(batch.source, batch.initial_columns)
        prepared.require_history(metric)
        matter = _family_matter(prepared, channel, ctx['coeff'], FROZEN_PLUS, ctx['target'])
        signed_id = f'{channel["group"]}_{channel["angular_sign"]}_E{batch.energy_sign:+d}'
        record_signed(signed_id, matter, channel, batch.energy_sign, [batch_id])
        states[batch_id] = prepared
    for pos_id, neg_id, neg_batch, neg_channel in group['maps']:
        partner = negative_angular_partner(states[pos_id], neg_batch.source)
        partner.require_history(metric)
        matter = _family_matter(partner, neg_channel, ctx['coeff'], FROZEN_PLUS, ctx['target'])
        signed_id = f'{neg_channel["group"]}_{neg_channel["angular_sign"]}_E{neg_batch.energy_sign:+d}'
        record_signed(signed_id, matter, neg_channel, neg_batch.energy_sign, [neg_id])
    arrays = {f'operator/{name}': np.asarray(getattr(operator, name)) for name in R.KEYS}
    arrays['paired_matter_change'] = sum(corrections.values())
    arrays['paired_high_matter_tangent'] = sum(tangents.values())
    meta = {'operator': {'digest': operator.digest, 'diagnostics': dict(operator.diagnostics)},
            'family_records': R.jsonable(family_records)}
    arrays['metadata_json'] = np.frombuffer(json.dumps(meta, sort_keys=True).encode(), np.uint8)
    return arrays, meta, sorted(corrections), operator


def new_family(ctx, key, cpu_limit):
    target_json, target_npz = family_paths(key)
    if target_json.exists() or target_npz.exists():
        raise FileExistsError('n=64 family output already exists')
    cpu = time.process_time()
    limit = min(FAMILY_CAP, float(cpu_limit))

    def stop(*_):
        raise TimeoutError('one n=64 mode operator exceeded its CPU allocation')

    previous = signal.signal(signal.SIGPROF, stop)
    signal.setitimer(signal.ITIMER_PROF, limit)
    try:
        arrays, meta, signed, operator = evaluate_family(ctx, key, saved=None, allow_new=True)
    finally:
        signal.setitimer(signal.ITIMER_PROF, 0)
        signal.signal(signal.SIGPROF, previous)
    parent_json, parent_npz = parent_family_paths(key)
    with np.load(parent_npz, allow_pickle=False) as handle:
        parent_change = np.asarray(handle['paired_matter_change'], float)
    evolved = np.asarray(arrays['paired_matter_change'], float)
    rel = float(np.max(np.abs(evolved - parent_change)) / max(float(np.max(np.abs(parent_change))), 1e-30))
    if rel > 0.05:
        raise ValueError(f'frozen iterate6 matter change drifted: relative {rel}')
    raw = deterministic_npz_bytes(arrays)
    record = {
        'positive_family': list(key),
        'operator_digest': operator.digest,
        'signed_families': signed,
        'family_records': meta['family_records'],
        'matter_relative_drift': rel,
        'CPU_seconds': time.process_time() - cpu,
        'new_operator_solves': 1,
        'source_hashes': ctx['hashes'],
        'profile_identity': ctx['identity'],
        'payload': {'path': str(target_npz.relative_to(ROOT)),
                    'sha256': sha256(raw).hexdigest(), 'bytes': len(raw)},
        'parent_payload_sha256': digest(parent_npz),
    }
    tmp = target_npz.with_suffix('.npz.tmp')
    tmp.write_bytes(raw)
    tmp.replace(target_npz)
    target_json.write_text(json.dumps(record, sort_keys=True)+'\n')
    return record


def load_tangent(ctx, key):
    path, payload = family_paths(key)
    record = json.loads(path.read_text())
    if digest(payload) != record['payload']['sha256']:
        raise ValueError('n=64 family payload changed')
    if record['source_hashes'] != ctx['hashes'] or record['profile_identity'] != ctx['identity']:
        raise ValueError('n=64 family belongs to a different history')
    with np.load(payload, allow_pickle=False) as handle:
        tangent = np.asarray(handle['paired_high_matter_tangent'], float)
    if tangent.shape != (HIGH, len(ctx['target']), 2):
        raise ValueError('high-mode tangent shape changed')
    return tangent, record


def assemble(ctx):
    records = []
    tangent = np.zeros((HIGH, len(ctx['target']), 2), float)
    family_records = {}
    for key in sorted(ctx['archived']):
        path, _payload = family_paths(key)
        if not path.exists():
            continue
        column, record = load_tangent(ctx, key)
        tangent += column
        records.append(record)
        family_records.update(record['family_records'])
    complete = len(records) == len(ctx['archived'])
    gradient = np.asarray(ctx['parent']['action_gradient'], float)
    base = np.asarray(ctx['parent']['history_jacobian'], float)
    report = {
        'completed_positive_families': len(records),
        'required_positive_families': len(ctx['archived']),
        'complete': complete,
    }
    if not complete:
        return report
    slots, slot_tangents = compatible_history_slots(ctx['metric'], ctx['target'], FROZEN_PLUS)
    geometry = surface_geometry_response(slots, slot_tangents, ctx['coeff'])
    weights = signed_family_angular_square_weights(family_records)
    first = ctx['archive'].family_entries(next(iter(ctx['archived'])))[0][0]
    phase = formal_source_phase_coefficient(
        ctx['metric'], ctx['target'], angular=1., rho_up=first.rho_up, gauss_nodes=192)
    _edge, edge_tangent = source_edge_from_phase(phase, weights, ctx['coeff']['a'])
    extra = geometry['action_gradient_tangent'][1:] + tangent + edge_tangent[1:]
    full = np.concatenate([base, extra], axis=0)
    residual = np.concatenate([gradient[:, 0], gradient[:, 1]])
    matrix = np.vstack([full[:, :, 0].T, full[:, :, 1].T])
    row = np.linalg.norm(matrix, axis=1)
    coln = np.linalg.norm(matrix, axis=0)
    row = np.where(row > 0, row, 1.0)
    coln = np.where(coln > 0, coln, 1.0)
    left, singular, right = np.linalg.svd(matrix / row[:, None] / coln[None, :], full_matrices=False)
    scaled_residual = residual / row
    n = len(gradient)
    parent = np.max(np.abs(gradient), axis=0)
    rows = []
    for count in (64, 80, 96, 112, 128):
        coef = np.zeros_like(singular)
        coef[:count] = (left[:, :count].T @ (-scaled_residual)) / singular[:count]
        delta = (right.T @ coef) / coln
        prediction = residual + matrix @ delta
        maxima = np.array([
            float(np.max(np.abs(prediction[:n]))),
            float(np.max(np.abs(prediction[n:]))),
        ])
        norm = float(np.linalg.norm(delta))
        used = delta if norm <= 1 else delta / norm
        clipped = residual + matrix @ used
        clipped_max = np.array([
            float(np.max(np.abs(clipped[:n]))),
            float(np.max(np.abs(clipped[n:]))),
        ])
        rows.append({
            'kept_modes': count,
            'step_norm': norm,
            'unclipped_maxima': maxima.tolist(),
            'unit_clipped_maxima': clipped_max.tolist(),
            'unclipped_reaches_existence': bool(np.all(maxima <= EXISTENCE)),
            'clipped_beats_parent': bool(np.all(clipped_max < parent)),
        })
    return {
        'schema': 'NSC-KS-COUPLED-NEWTON-N64-HIGH-MODES-v1',
        'accountable_author': 'Douglas Ek',
        'status': 'OPEN: n=64 mode columns measured at the iterate6 history',
        'profile_identity': ctx['identity'],
        'parent_constraint_maxima': parent.tolist(),
        'ranks': rows,
        'completed_positive_families': len(records),
        'required_positive_families': len(ctx['archived']),
        'physical_EXISTENCE_certificate': False,
        'physical_NONEXISTENCE_certificate': False,
        'source_hashes': ctx['hashes'],
    }


def run(cpu_budget, max_new):
    ctx = context()
    DIRECTORY.mkdir(parents=True, exist_ok=True)
    start = time.process_time()
    new = 0
    for key in sorted(ctx['archived']):
        path, _payload = family_paths(key)
        if path.exists():
            continue
        if new >= max_new:
            continue
        remaining = cpu_budget - (time.process_time() - start)
        if remaining <= 1:
            break
        record = new_family(ctx, key, remaining - 1)
        new += 1
        print(json.dumps({
            'positive_family': record['positive_family'],
            'CPU_seconds': record['CPU_seconds'],
            'matter_relative_drift': record['matter_relative_drift'],
        }), flush=True)
    result = assemble(ctx)
    if result.get('completed_positive_families') == result.get('required_positive_families'):
        OUTPUT.write_text(json.dumps(result, sort_keys=True, indent=2)+'\n')
    return result


def stored_family_hashes():
    ctx = context()
    stored = None
    keys = []
    for key in sorted(ctx['archived']):
        path, payload = family_paths(key)
        if not path.exists():
            raise ValueError('n=64 family is missing: '+path.name)
        record = json.loads(path.read_text())
        if digest(payload) != record['payload']['sha256']:
            raise ValueError('n=64 family payload changed: '+path.name)
        if stored is None:
            stored = record['source_hashes']
        elif record['source_hashes'] != stored:
            raise ValueError('n=64 family hash sets disagree')
        keys.append(key)
    if len(keys) != len(ctx['archived']):
        raise ValueError('n=64 columns are incomplete')
    return ctx, stored, keys


def replay_family(key, ctx=None):
    if ctx is None:
        ctx, _stored, _keys = stored_family_hashes()
    path, payload = family_paths(key)
    record = json.loads(path.read_text())
    with np.load(payload, allow_pickle=False) as handle:
        arrays = {name: np.array(handle[name], copy=True) for name in handle.files}
    meta = json.loads(arrays['metadata_json'].tobytes())
    current, _meta, signed, operator = evaluate_family(
        ctx, key, saved=(arrays, meta), allow_new=False)
    if signed != record['signed_families']:
        raise ValueError('n=64 replay signed coverage changed')
    if operator.digest != record['operator_digest']:
        raise ValueError('n=64 replay operator digest changed')
    if not np.array_equal(current['paired_high_matter_tangent'], arrays['paired_high_matter_tangent']):
        raise ValueError('n=64 replay tangent changed')
    if not np.array_equal(current['paired_matter_change'], arrays['paired_matter_change']):
        raise ValueError('n=64 replay matter change changed')
    return {'positive_family': list(key), 'operator_digest': operator.digest}


def replay(family_key=(1, -1)):
    ctx, stored, _keys = stored_family_hashes()
    replayed = replay_family(family_key, ctx=ctx)
    ctx['hashes'] = stored
    result = assemble(ctx)
    saved = json.loads(OUTPUT.read_text())
    if result.get('ranks') != saved['ranks'] or result.get('profile_identity') != saved['profile_identity']:
        raise ValueError('n=64 column replay ranks changed')
    result['replayed_family'] = replayed['positive_family']
    result['replay_dirac_evolution'] = False
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--run', action='store_true')
    mode.add_argument('--check', action='store_true')
    parser.add_argument('--cpu-budget', type=float, default=7200.0)
    parser.add_argument('--max-new', type=int, default=60)
    args = parser.parse_args()
    result = run(args.cpu_budget, args.max_new) if args.run else replay()
    print(json.dumps({
        'completed_positive_families': result.get('completed_positive_families'),
        'required_positive_families': result.get('required_positive_families'),
        'ranks': result.get('ranks'),
        'status': result.get('status'),
        'replayed_family': result.get('replayed_family'),
        'replay_dirac_evolution': result.get('replay_dirac_evolution'),
    }, indent=2))
