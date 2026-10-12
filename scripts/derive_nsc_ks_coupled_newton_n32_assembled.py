#!/usr/bin/env python3
"""Assembled Chebyshev n=32 Jacobian at the accepted iterate4 history.

Zero-pads e15982e2… from 16 to 32 coefficients per function. Reuses the
already-evolved n=16 matter columns by index remap. Evolves ONLY the new
high-mode directions (w16..31, U16..31) through the retained-source
evaluator path, with the low-mode history frozen as one amplitude-1
direction so the base radius stays the same physical g. No slot-transfer,
no geometry-only map, no invented matter columns.
"""
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
import derive_nsc_ks_coupled_newton_n16 as N16
import derive_nsc_ks_coupled_newton_n16_damped as D
import derive_nsc_ks_coupled_newton_n16_damped_iterate4 as I4
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
from recursive_horizons.nsc_ks_profile_identity import profile_identity
from recursive_horizons.nsc_ks_signed_state import negative_angular_partner
from recursive_horizons.nsc_ks_source_envelope import (
    KSEnvelopeBinding, _sample_axial_profiles, usual_axial_support)
from recursive_horizons.nsc_local_history_newton import (
    LocalHistoryNewtonSettings, _clip_step, _linear_step)
from recursive_horizons.nsc_local_incoming_constraints import _channel_record
from recursive_horizons.nsc_local_incoming_family import (
    LocalAxialFunction, LocalIncomingFamily)
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes

DIRECTORY = ROOT/'results/development/artifacts/nsc-ks-coupled-newton-n32-assembled'
OUTPUT = ROOT/'results/development/nsc-ks-coupled-newton-n32-assembled.json'
I4_ARTIFACTS = I4.DIRECTORY
EXPECTED_PARENT = I4.EXPECTED_IDENTITY
EXPECTED_PADDED = '49223e613233a84ece27eb4efd544eb462b4343453d7dcbc8f2cf0cfec9dff68'
PARENT_RESIDUAL = (0.004107193193038788, 0.005790788097514318)
EXISTENCE_TOLERANCE = 3e-11
CONDITION_BLOWUP = 1e15
SOLVE_COUNT = 16
VERIFY_COUNT = 33
RETARDED = 64
HIGH_RETARDED = 32
FROZEN_PLUS_HIGH = 33
PHASE_NODES = 192
MAX_NEW = 60
CPU_CAP = 1200.
FAMILY_CAP = 40.
OWNERS = (
    'scripts/derive_nsc_ks_coupled_newton_n32_assembled.py',
    'docs/nsc-ks-coupled-newton-n32-assembled.md',
    'scripts/derive_nsc_ks_coupled_newton_n16_damped_iterate4.py',
    'scripts/derive_nsc_ks_coupled_newton_n16.py',
    'scripts/derive_nsc_ks_coupled_retained_control.py',
    'src/recursive_horizons/nsc_ks_history_evaluator.py',
    'src/recursive_horizons/nsc_ks_cutoff_bridge_evaluator.py',
    'src/recursive_horizons/nsc_dirac_source_phase.py',
    'src/recursive_horizons/nsc_ks_energy_propagator.py',
    'src/recursive_horizons/nsc_ks_difference_envelope.py',
    'src/recursive_horizons/nsc_ks_batched_constraints.py',
    'src/recursive_horizons/nsc_local_incoming_family.py',
)


def digest(path):
    path = Path(path)
    return sha256((path if path.is_absolute() else ROOT/path).read_bytes()).hexdigest()


def family_paths(key):
    stem = f'family-{key[0]:02d}-{key[1]:+d}'
    return DIRECTORY/(stem+'.json'), DIRECTORY/(stem+'.npz')


def i4_family_paths(key):
    stem = f'family-{key[0]:02d}-{key[1]:+d}'
    return I4_ARTIFACTS/(stem+'.json'), I4_ARTIFACTS/(stem+'.npz')


def zero_pad(coefficients):
    raw = np.asarray(coefficients, float)
    if raw.shape != (2, 16):
        raise ValueError('n=32 assembled start must pad the accepted 16-coefficient history')
    padded = np.zeros((2, 32), float)
    padded[:, :16] = raw
    family = LocalIncomingFamily(padded)
    if family.radius_lower_bound() <= 0:
        raise ValueError('zero-padded n=32 family must keep a positive radius bound')
    return family


def high_mode_metric(family32):
    """Frozen low-mode history plus the 32 new Chebyshev directions only."""
    coeff = np.asarray(family32.coefficients, float)
    low = coeff[:, :16]
    center, ai, ao = family32.center, family32.axial_inner, family32.axial_outer
    ni, no = family32.normal_inner, family32.normal_outer
    zero = LocalAxialFunction((0.,), center, ai, ao)
    frozen = CompatibleRadiusDirection(
        LocalAxialFunction(tuple(low[0]), center, ai, ao),
        LocalAxialFunction(tuple(low[1]), center, ai, ao),
        ni, no)
    basis = [LocalAxialFunction(tuple(row), center, ai, ao) for row in np.eye(32)]
    high_w = tuple(CompatibleRadiusDirection(basis[i], zero, ni, no) for i in range(16, 32))
    high_u = tuple(CompatibleRadiusDirection(zero, basis[i], ni, no) for i in range(16, 32))
    return CompatibleIncomingMetric((1.,) + (0.,) * HIGH_RETARDED, (frozen,) + high_w + high_u)


def remap_low_tangent(tangent16):
    """Map n=16 (32,nz,2) into n=32 (64,nz,2) low-mode slots; high slots stay zero."""
    src = np.asarray(tangent16, float)
    if src.ndim != 3 or src.shape[0] != 32:
        raise ValueError('n=16 matter tangent must have shape (32, nz, 2)')
    out = np.zeros((RETARDED, src.shape[1], src.shape[2]), float)
    out[0:16] = src[0:16]
    out[32:48] = src[16:32]
    return out


def splice_high_tangent(low64, high32):
    out = np.array(low64, float, copy=True)
    high = np.asarray(high32, float)
    if out.shape[0] != RETARDED or high.shape[0] != HIGH_RETARDED:
        raise ValueError('spliced tangent shapes must be (64,*,*) and (32,*,*)')
    if out.shape[1:] != high.shape[1:]:
        raise ValueError('high-mode tangent node shape mismatch')
    out[16:32] = high[0:16]
    out[48:64] = high[16:32]
    return out


def flatten_system(gradient, tangent):
    residual = np.concatenate([gradient[:, 0], gradient[:, 1]])
    jacobian = np.vstack([tangent[:, :, 0].T, tangent[:, :, 1].T])
    return residual, jacobian, len(gradient)


def leftover_of(gradient, tangent, settings):
    residual, jacobian, n_nodes = flatten_system(gradient, tangent)
    delta, diagnostics = _linear_step(jacobian, residual, jacobian, settings)
    if delta is None:
        pred = residual
        clipped_pred = residual
        clip_norm = None
        step_norm = None
    else:
        pred = residual + jacobian @ delta
        clipped, clip_norm = _clip_step(delta, settings.max_step_norm)
        clipped_pred = residual + jacobian @ clipped
        step_norm = float(np.linalg.norm(delta))
    return {
        'unclipped': {
            'N_max': float(np.max(abs(pred[:n_nodes]))),
            'beta_max': float(np.max(abs(pred[n_nodes:]))),
            'control_max': float(np.max(abs(pred))),
            'step_norm': step_norm,
        },
        'clipped': {
            'N_max': float(np.max(abs(clipped_pred[:n_nodes]))),
            'beta_max': float(np.max(abs(clipped_pred[n_nodes:]))),
            'control_max': float(np.max(abs(clipped_pred))),
            'step_norm': clip_norm,
        },
        'condition_number': float(diagnostics['condition_number']),
        'singular_values': np.asarray(diagnostics['singular_values'], float),
        'rank_or_conditioning': bool(diagnostics['rank_or_conditioning']),
        'delta': None if delta is None else np.asarray(delta, float),
    }


def reconstruct_operator(interval, metric, grid, target, mass, angular, rho_up, arrays, diagnostics):
    w, u = _sample_axial_profiles(metric.directions, grid)
    binding = KSEnvelopeBinding(
        grid, w, u, tuple(metric.amplitudes),
        tuple(d.inner_radius for d in metric.directions),
        tuple(d.outer_radius for d in metric.directions),
        usual_axial_support(), rho_up, 1., **R.OPTIONS)
    operator = KSEnergyPropagator(
        interval, *(np.asarray(arrays[key]) for key in R.KEYS),
        mass, angular, rho_up, binding, diagnostics)
    if not np.array_equal(operator.z, target):
        raise ValueError('saved high-mode operator target nodes changed')
    if operator.tangent.shape[1] != FROZEN_PLUS_HIGH:
        raise ValueError('frozen-plus-high retarded directions required')
    return operator


def context():
    base = R.context()
    trial = json.loads(I4.OUTPUT.read_text())
    for path, expected in trial['source_hashes'].items():
        if digest(path) != expected:
            raise ValueError('iterate4 input changed: '+path)
    if trial['profile_identity'] != EXPECTED_PARENT:
        raise ValueError('n=32 assembled must start from accepted iterate4 identity')
    if trial['completed_positive_families'] != trial['required_positive_families']:
        raise ValueError('incomplete iterate4 cannot seed n=32 assembled')
    measured = np.asarray(trial['constraint_maxima'], float)
    if not np.allclose(measured, PARENT_RESIDUAL, rtol=0, atol=1e-15):
        raise ValueError('iterate4 residual maxima changed')
    family = zero_pad(trial['history']['coefficients'])
    identity = profile_identity(family, include_normal_window=True)
    if identity != EXPECTED_PADDED:
        raise ValueError('zero-padded n=32 profile identity changed')
    if identity in D.FORBIDDEN_IDENTITIES or EXPECTED_PARENT in D.FORBIDDEN_IDENTITIES:
        raise ValueError('forbidden identity selected')
    solve = family.collocation_nodes(SOLVE_COUNT)
    verify = family.collocation_nodes(VERIFY_COUNT)
    target = H.fixed_target_union(solve, verify)
    if len(target) != 47:
        raise ValueError('n=32 assembled solve/verification union must have 47 nodes')
    if not np.array_equal(target, np.asarray(trial['z'], float)):
        raise ValueError('n=32 assembled must keep the iterate4 target union')
    metric = high_mode_metric(family)
    hashes = {
        **base['archive'].input_hashes,
        **{path: digest(path) for path in OWNERS},
        str(I4.OUTPUT.relative_to(ROOT)): digest(I4.OUTPUT),
    }
    return {
        **base,
        'family': family,
        'high_metric': metric,
        'identity': identity,
        'parent_identity': EXPECTED_PARENT,
        'solve': solve,
        'verify': verify,
        'target': target,
        'hashes': hashes,
        'parent': trial,
        'measured': measured,
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
    metric = ctx['high_metric']
    saved_arrays = saved_meta = None
    if saved is not None:
        saved_arrays, saved_meta = saved
    if saved_arrays is not None:
        operator = reconstruct_operator(
            interval, metric, ctx['grid'], ctx['target'], first.mass, first.angular,
            first.rho_up, {name: saved_arrays['operator/'+name] for name in R.KEYS},
            saved_meta['operator']['diagnostics'])
        if operator.digest != saved_meta['operator']['digest']:
            raise ValueError('saved high-mode operator digest changed')
    else:
        if not allow_new:
            raise RuntimeError('replay must not call Dirac evolution')
        operator = evolve_energy_propagator(
            interval, R.DEGREE, metric, ctx['grid'], ctx['target'],
            first.mass, first.angular, first.rho_up,
            axial_support=usual_axial_support(), tangents='all', **R.OPTIONS)
    if operator.tangent.shape[1] != FROZEN_PLUS_HIGH:
        raise ValueError('high-mode evolve must retain frozen-plus-high directions')

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

    corrections = {}
    tangents = {}
    family_records = {}
    states = {}
    for batch_id, batch, channel in group['applies']:
        prepared = operator.apply(batch.source, batch.initial_columns)
        prepared.require_history(metric)
        if not np.array_equal(prepared.z, ctx['target']):
            raise ValueError('operator targets must be the frozen solver/verification union')
        if prepared.column_tangents.shape[0] != FROZEN_PLUS_HIGH:
            raise ValueError('prepared high-mode tangents must include the frozen direction')
        matter = _family_matter(
            prepared, channel, ctx['coeff'], FROZEN_PLUS_HIGH, ctx['target'])
        signed_id = f'{channel["group"]}_{channel["angular_sign"]}_E{batch.energy_sign:+d}'
        record_signed(signed_id, matter, channel, batch.energy_sign, [batch_id])
        states[batch_id] = prepared
    for pos_id, neg_id, neg_batch, neg_channel in group['maps']:
        partner = negative_angular_partner(states[pos_id], neg_batch.source)
        partner.require_history(metric)
        matter = _family_matter(
            partner, neg_channel, ctx['coeff'], FROZEN_PLUS_HIGH, ctx['target'])
        signed_id = (
            f'{neg_channel["group"]}_{neg_channel["angular_sign"]}_E{neg_batch.energy_sign:+d}')
        record_signed(signed_id, matter, neg_channel, neg_batch.energy_sign, [neg_id])
    expected_ids = {f'{key[0]}_{key[1]}_E+1', f'{key[0]}_{-key[1]}_E-1'}
    signed = sorted(corrections)
    if set(signed) != expected_ids or len(signed) != 2:
        raise ValueError('explicit signed partners required for this family')
    for name, value in tangents.items():
        if value.shape != (HIGH_RETARDED, len(ctx['target']), 2):
            raise ValueError('high-mode matter tangent must have shape (32, nz, 2)')
    arrays = {f'operator/{name}': np.asarray(getattr(operator, name)) for name in R.KEYS}
    arrays['paired_matter_change'] = sum(corrections.values())
    arrays['paired_high_matter_tangent'] = sum(tangents.values())
    for name, value in corrections.items():
        arrays['correction/'+name] = value
    for name, value in tangents.items():
        arrays['tangent/'+name] = value
    operator_record = {
        'mass': float(operator.mass), 'angular': float(operator.angular),
        'rho_up': float(operator.rho_up), 'digest': operator.digest,
        'diagnostics': dict(operator.diagnostics),
        'reused_from_coupled_pilot': False,
        'high_mode_only': True,
        'frozen_plus_high_directions': FROZEN_PLUS_HIGH,
        'recorded_high_matter_directions': HIGH_RETARDED,
    }
    meta = {
        'operator': operator_record,
        'family_records': R.jsonable(family_records),
    }
    arrays['metadata_json'] = np.frombuffer(json.dumps(meta, sort_keys=True).encode(), np.uint8)
    return arrays, meta, signed, operator


def new_family(ctx, key, cpu_limit=FAMILY_CAP):
    target_json, target_npz = family_paths(key)
    if target_json.exists() or target_npz.exists():
        raise FileExistsError('family output already exists')
    cpu = time.process_time()
    limit = min(FAMILY_CAP, float(cpu_limit))
    if not np.isfinite(limit) or limit <= 0:
        raise ValueError('positive remaining family CPU allocation required')

    def stop(*_):
        raise TimeoutError('one n=32 high-mode operator exceeded its remaining CPU allocation')

    previous = signal.signal(signal.SIGPROF, stop)
    signal.setitimer(signal.ITIMER_PROF, limit)
    try:
        arrays, meta, signed, operator = evaluate_family(ctx, key, saved=None, allow_new=True)
    finally:
        signal.setitimer(signal.ITIMER_PROF, 0)
        signal.signal(signal.SIGPROF, previous)
    archived = ctx['archived'][key]
    entries = ctx['archive'].family_entries(key)
    pos_rows = sum(len(batch.source.energies)//3 for batch, _ in entries if batch.energy_sign > 0)
    if pos_rows != archived['source_rows']:
        raise ValueError('restored rows differ from the original family coverage')
    parent_json, parent_npz = i4_family_paths(key)
    parent = json.loads(parent_json.read_text())
    if digest(parent_npz) != parent['payload']['sha256']:
        raise ValueError('iterate4 family payload changed: '+str(parent_npz))
    with np.load(parent_npz, allow_pickle=False) as f:
        parent_change = np.asarray(f['paired_matter_change'], float)
    # Adaptive DOP853 couples tangent state into step selection, so the frozen
    # high-mode base field can drift by ~1e-2 relative. Assembly always reuses
    # the iterate4 matter change; this check only guards against order-of-
    # magnitude mistakes.
    evolved_change = np.asarray(arrays['paired_matter_change'], float)
    rel = float(np.max(abs(evolved_change - parent_change))
                / max(float(np.max(abs(parent_change))), 1e-30))
    if rel > 0.05:
        raise ValueError(
            'high-mode base matter change drifted too far from iterate4: '
            f'relative {rel}')
    arrays['paired_matter_change_high_mode_evolve'] = evolved_change
    arrays['paired_matter_change'] = parent_change
    raw = deterministic_npz_bytes(arrays)
    record = {
        'schema': 'NSC-KS-COUPLED-NEWTON-N32-ASSEMBLED-FAMILY-v1',
        'positive_family': list(key),
        'operator_energy_interval': list(R.interpolation_interval(archived)),
        'interpolant_degree': R.DEGREE,
        'operator_digest': operator.digest,
        'operator': meta['operator'],
        'source_rows': pos_rows,
        'source_panels': list(archived['source_panels']),
        'quadrature_measure': float(archived['quadrature_measure']),
        'signed_families': signed,
        'family_records': meta['family_records'],
        'matter_change_maxima': np.max(abs(parent_change), axis=0).tolist(),
        'high_mode_evolve_matter_change_maxima': np.max(abs(evolved_change), axis=0).tolist(),
        'high_mode_evolve_matter_relative_drift': rel,
        'high_matter_tangent_maxima': np.max(abs(arrays['paired_high_matter_tangent']), axis=(0, 1)).tolist(),
        'CPU_seconds': time.process_time()-cpu,
        'family_CPU_cap': FAMILY_CAP,
        'new_operator_solves': 1,
        'upstream_batches': len(entries)//2,
        'source_hashes': ctx['hashes'],
        'profile_identity': ctx['identity'],
        'parent_profile_identity': ctx['parent_identity'],
        'physical_local_gate': 'OPEN',
        'source_and_numerical_error_bound': None,
        'payload': {'path': str(target_npz.relative_to(ROOT)),
                    'sha256': sha256(raw).hexdigest(), 'bytes': len(raw)},
        'iterate4_family': {
            'path': str(parent_json.relative_to(ROOT)),
            'sha256': digest(parent_json),
            'payload_sha256': parent['payload']['sha256'],
        },
    }
    tmp = target_npz.with_suffix('.npz.tmp')
    tmp.write_bytes(raw)
    tmp.replace(target_npz)
    target_json.write_text(R.dump(record))
    return record


def read_family(ctx, key, replay=False):
    path, payload = family_paths(key)
    record = json.loads(path.read_text())
    archived = ctx['archived'][key]
    if record['source_hashes'] != ctx['hashes'] or record['positive_family'] != list(key):
        raise ValueError('family signature/identity changed')
    if record.get('profile_identity') != ctx['identity']:
        raise ValueError('saved operator belongs to a different history g')
    if record['source_rows'] != archived['source_rows'] or record['source_panels'] != archived['source_panels']:
        raise ValueError('family source coverage changed')
    if list(record['operator_energy_interval']) != list(R.interpolation_interval(archived)):
        raise ValueError('family interpolation interval changed')
    if record['interpolant_degree'] != R.DEGREE:
        raise ValueError('family interpolant degree changed')
    if digest(record['payload']['path']) != record['payload']['sha256']:
        raise ValueError('family payload changed')
    with np.load(ROOT/record['payload']['path'], allow_pickle=False) as f:
        arrays = {name: np.array(f[name], copy=True) for name in f.files}
    meta = json.loads(arrays['metadata_json'].tobytes())
    if meta['operator']['digest'] != record['operator_digest']:
        raise ValueError('family operator metadata changed')
    if sorted(name.split('/', 1)[1] for name in arrays if name.startswith('correction/')) != record['signed_families']:
        raise ValueError('signed family coverage changed')
    if not np.array_equal(sum(value for name, value in arrays.items() if name.startswith('correction/')),
                          arrays['paired_matter_change_high_mode_evolve']):
        raise ValueError('signed response sum changed')
    if not np.array_equal(sum(value for name, value in arrays.items() if name.startswith('tangent/')),
                          arrays['paired_high_matter_tangent']):
        raise ValueError('signed high-mode tangent sum changed')
    if not np.array_equal(arrays['paired_matter_change'], parent_matter_change(ctx, key)):
        raise ValueError('assembled matter change must equal the iterate4 family change')
    if replay:
        current, current_meta, signed, operator = evaluate_family(
            ctx, key, saved=(arrays, meta), allow_new=False)
        if signed != record['signed_families']:
            raise ValueError('family replay signed coverage changed')
        if operator.digest != record['operator_digest']:
            raise ValueError('family replay operator digest changed')
        for name, value in current.items():
            if name in ('metadata_json', 'paired_matter_change'):
                continue
            if name == 'paired_matter_change_high_mode_evolve':
                continue
            if not np.array_equal(value, arrays[name]):
                raise ValueError('family replay changed: '+name)
        if not np.array_equal(
                current['paired_matter_change'],
                arrays['paired_matter_change_high_mode_evolve']):
            raise ValueError('family replay evolved matter change changed')
        if R.jsonable(current_meta['family_records']) != record['family_records']:
            raise ValueError('family replay records changed')
    return (record, arrays['paired_matter_change'], arrays['paired_high_matter_tangent'],
            record['family_records'])


def load_remapped_low_tangent(ctx, key):
    path, payload = i4_family_paths(key)
    record = json.loads(path.read_text())
    if digest(payload) != record['payload']['sha256']:
        raise ValueError('iterate4 family payload changed during remap')
    if record.get('profile_identity') != ctx['parent_identity']:
        raise ValueError('iterate4 family belongs to a different parent g')
    with np.load(payload, allow_pickle=False) as f:
        tangent = np.asarray(f['paired_matter_tangent'], float)
        change = np.asarray(f['paired_matter_change'], float)
    return remap_low_tangent(tangent), change, record


def parent_matter_change(ctx, key):
    return load_remapped_low_tangent(ctx, key)[1]


def damped_proposals(family, gradient, jacobian, measured, settings):
    residual, matrix, _ = flatten_system(gradient, jacobian)
    delta, diagnostics = _linear_step(matrix, residual, matrix, settings)
    rows = []
    selected = None
    if delta is None:
        return rows, selected, diagnostics, None
    clipped, clip_norm = _clip_step(delta, settings.max_step_norm)
    scale = 1.0
    while scale >= settings.min_line_search_scale - 1e-15:
        step = scale * clipped
        candidate = LocalIncomingFamily(
            np.asarray(family.coefficients, float) + step.reshape(2, 32))
        radius = candidate.radius_lower_bound()
        identity = profile_identity(candidate, include_normal_window=True)
        predicted = D.predicted_maxima(gradient, jacobian, step)
        forbidden = identity in D.FORBIDDEN_IDENTITIES or identity == EXPECTED_PARENT
        improves = radius > 0 and not forbidden and D.beats_measured(predicted, measured)
        row = {
            'step_scale': scale,
            'step_norm': float(np.linalg.norm(step)),
            'clipped_direction_norm': float(clip_norm),
            'radius_lower_bound': radius,
            'candidate_profile_identity': identity,
            'predicted_all_node_residual_maxima': predicted.tolist(),
            'prediction_beats_measured': improves,
            'forbidden_identity': forbidden,
            'positive_radius': radius > 0,
            'prediction_is_not_a_nonlinear_residual': True,
            'layout': 'rectangular-n32-assembled',
            'direction': step.reshape(2, 32).tolist(),
            'unclipped_step_norm': float(np.linalg.norm(delta)),
        }
        rows.append(row)
        if selected is None and improves:
            selected = row
        scale *= settings.line_search_factor
    return rows, selected, diagnostics, delta


def assemble(ctx, replay=False):
    records = []
    matter = np.zeros((len(ctx['target']), 2), float)
    low_tangent = np.zeros((RETARDED, len(ctx['target']), 2), float)
    high_tangent = np.zeros((HIGH_RETARDED, len(ctx['target']), 2), float)
    signed = []
    family_records = {}
    for key in sorted(ctx['archived']):
        path, payload = family_paths(key)
        if not path.exists():
            if payload.exists():
                raise ValueError('orphan operator payload retained for recovery: '+str(payload))
            continue
        record, delta, dJ_high, local_records = read_family(ctx, key, replay)
        remapped, parent_change, parent_record = load_remapped_low_tangent(ctx, key)
        if not np.array_equal(delta, parent_change):
            raise ValueError('assembled matter change must match remapped iterate4 family')
        records.append(record)
        matter = matter + parent_change
        low_tangent = low_tangent + remapped
        high_tangent = high_tangent + dJ_high
        signed.extend(record['signed_families'])
        overlap = set(family_records) & set(local_records)
        if overlap:
            raise ValueError('signed family counted twice: '+sorted(overlap)[0])
        family_records.update(local_records)
        _ = parent_record
    if len(signed) != len(set(signed)):
        raise ValueError('signed family counted twice')
    tangent = splice_high_tangent(low_tangent, high_tangent)
    slots, slot_tangents = compatible_history_slots(ctx['family'].metric(), ctx['target'], RETARDED)
    geometry = surface_geometry_response(slots, slot_tangents, ctx['coeff'])
    baseline = np.broadcast_to(ctx['baseline'], (len(ctx['target']), 2)).copy()
    raw_gradient = baseline + geometry['action_gradient_change'] + matter
    raw_jacobian = geometry['action_gradient_tangent'] + tangent
    complete = len(records) == len(ctx['archived'])
    weights = edge = edge_tangent = None
    phase_binding = None
    if complete:
        if len(family_records) != 120:
            raise ValueError('one hundred twenty explicit signed families required')
        weights = signed_family_angular_square_weights(family_records)
        if not np.allclose(weights, R.EXPECTED_WEIGHT, rtol=0, atol=1e-6):
            raise ValueError('signed angular-square weights must be 107952 for each energy sign')
        first = ctx['archive'].family_entries(next(iter(ctx['archived'])))[0][0]
        phase = formal_source_phase_coefficient(
            ctx['family'], ctx['target'], angular=1., rho_up=first.rho_up, gauss_nodes=PHASE_NODES)
        if phase.profile_identity != ctx['identity'] or not np.array_equal(phase.z, ctx['target']):
            raise ValueError('unit-angular phase must use the complete target union and the same g')
        if phase.delta_f_z.shape[1] != RETARDED:
            raise ValueError('phase must retain all64 retarded directions')
        edge, edge_tangent = source_edge_from_phase(phase, weights, ctx['coeff']['a'])
        phase_binding = phase.binding
        gradient = raw_gradient + edge
        jacobian = raw_jacobian + edge_tangent
        joined_matter = matter + edge
        joined_tangent = tangent + edge_tangent
    else:
        gradient = raw_gradient
        jacobian = raw_jacobian
        joined_matter = matter
        joined_tangent = tangent
    parent_gradient = np.asarray(ctx['parent']['action_gradient'], float)
    measured = np.max(abs(gradient), axis=0)
    settings = LocalHistoryNewtonSettings()
    leftover = leftover_of(gradient, jacobian, settings) if complete else None
    proposals = selected = diagnostics = delta = None
    evolve_authorized = False
    named_stall = None
    if complete:
        if not np.allclose(gradient, parent_gradient, rtol=1e-9, atol=1e-12):
            raise ValueError('zero-padded n=32 residual must match iterate4 at the same physical g')
        if leftover['condition_number'] >= CONDITION_BLOWUP:
            named_stall = 'n32_assembled_condition_blowup_finite_image'
        elif leftover['unclipped']['N_max'] > EXISTENCE_TOLERANCE or leftover['unclipped']['beta_max'] > EXISTENCE_TOLERANCE:
            named_stall = 'finite_wu_n32_assembled_image_cannot_reach_existence'
        proposals, selected, diagnostics, delta = damped_proposals(
            ctx['family'], gradient, jacobian, ctx['measured'], settings)
        # Authorize a trial only when clipped prediction beats the parent residual
        # and the kept-mode leftover is not stuck above the gate with a blown condition.
        condition_ok = leftover['condition_number'] < CONDITION_BLOWUP
        leftover_blocks = (
            leftover['unclipped']['N_max'] > EXISTENCE_TOLERANCE
            or leftover['unclipped']['beta_max'] > EXISTENCE_TOLERANCE)
        wild_step = (
            leftover['unclipped']['step_norm'] is not None
            and leftover['unclipped']['step_norm'] > 10.0
            and leftover_blocks)
        evolve_authorized = bool(
            selected is not None and condition_ok and not wild_step
            and leftover['clipped']['N_max'] < PARENT_RESIDUAL[0]
            and leftover['clipped']['beta_max'] < PARENT_RESIDUAL[1])
        if selected is None and named_stall is None:
            named_stall = 'n32_assembled_clipped_prediction_does_not_improve'
        if wild_step and named_stall is None:
            named_stall = 'finite_wu_n32_assembled_image_cannot_reach_existence'
    status = (
        'OPEN: complete n=32 assembled residual of the accepted iterate4 history; '
        'uncertified physical gate' if complete
        else 'OPEN: partial n=32 assembled residual of the accepted iterate4 history')
    if complete and named_stall is not None:
        status = (
            'OPEN: assembled Chebyshev n=32 leftover cannot reach EXISTENCE; '
            'finite-image stall; no wild full-scale step')
    return {
        'schema': 'NSC-KS-COUPLED-NEWTON-N32-ASSEMBLED-v1',
        'accountable_author': 'Douglas Ek',
        'status': status,
        'history': ctx['family'].description(),
        'profile_identity': ctx['identity'],
        'parent_profile_identity': ctx['parent_identity'],
        'layout': 'rectangular-n32-assembled-high-mode',
        'state_law': 'C_Sigma[g]=U_g C_up U_g†; original fixed source',
        'completed_positive_families': len(records),
        'required_positive_families': len(ctx['archived']),
        'all_retained_nonzero_angular_families_included': complete,
        'positive_energy_rows': sum(item['source_rows'] for item in records),
        'signed_family_ids': sorted(signed),
        'grid_nodes': R.GRID_NODES,
        'energy_interpolation_degree': R.DEGREE,
        'retarded_directions': RETARDED,
        'high_mode_directions_evolved': HIGH_RETARDED,
        'low_mode_directions_reused': 32,
        'options': R.OPTIONS,
        'target_nodes': ctx['target'].tolist(),
        'z': ctx['target'].tolist(),
        'summed_matter_change': matter.tolist(),
        'summed_matter_tangent': tangent.tolist(),
        'summed_high_matter_tangent': high_tangent.tolist() if complete else None,
        'raw_action_gradient': raw_gradient.tolist(),
        'raw_history_jacobian': raw_jacobian.tolist(),
        'source_cutoff_edge_gradient': None if edge is None else edge.tolist(),
        'source_cutoff_edge_history_tangent': None if edge_tangent is None else edge_tangent.tolist(),
        'action_gradient': gradient.tolist(),
        'history_jacobian': jacobian.tolist(),
        'constraint_maxima': measured.tolist(),
        'raw_constraint_maxima': np.max(abs(raw_gradient), axis=0).tolist(),
        'matter_change_maxima': np.max(abs(matter), axis=0).tolist(),
        'joined_matter_change_maxima': np.max(abs(joined_matter), axis=0).tolist(),
        'retarded_matter_plus_edge_jacobian_maxima': np.max(abs(joined_tangent), axis=(0, 1)).tolist(),
        'source_cutoff_angular_square_weights': None if weights is None else np.asarray(weights).tolist(),
        'source_cutoff_phase_binding': phase_binding,
        'source_cutoff_phase_gauss_nodes': PHASE_NODES,
        'previous_constraint_maxima': ctx['measured'].tolist(),
        'parent_constraint_maxima': ctx['measured'].tolist(),
        'prediction_is_not_a_nonlinear_residual': True,
        'measured_all_node_residual_maxima': measured.tolist() if complete else None,
        'radius_lower_bound': ctx['family'].radius_lower_bound(),
        'assembled_unclipped_leftover': None if leftover is None else [
            leftover['unclipped']['N_max'], leftover['unclipped']['beta_max']],
        'assembled_clipped_leftover': None if leftover is None else [
            leftover['clipped']['N_max'], leftover['clipped']['beta_max']],
        'assembled_condition_number': None if leftover is None else leftover['condition_number'],
        'assembled_unclipped_step_norm': None if leftover is None else leftover['unclipped']['step_norm'],
        'unclipped_leftover_reaches_existence_tolerance': None if leftover is None else bool(
            leftover['unclipped']['N_max'] <= EXISTENCE_TOLERANCE
            and leftover['unclipped']['beta_max'] <= EXISTENCE_TOLERANCE),
        'proposals': proposals,
        'selected_proposal': selected,
        'linear_step_diagnostics': None if diagnostics is None else {
            'condition_number': float(diagnostics['condition_number']),
            'rank_or_conditioning': bool(diagnostics['rank_or_conditioning']),
        },
        'evolve_authorized': evolve_authorized,
        'families_evolved': len(records),
        'named_search_stall': named_stall,
        'slot_transfer': False,
        'geometry_only_map': False,
        'invented_matter_columns': False,
        'forbidden_profiles_evolved': False,
        'family_records': [{'path': str(family_paths(tuple(item['positive_family']))[0].relative_to(ROOT)),
                            'sha256': digest(family_paths(tuple(item['positive_family']))[0]),
                            'operator_digest': item['operator_digest'],
                            'new_operator_solves': item['new_operator_solves']}
                           for item in records],
        'CPU_seconds': sum(item['CPU_seconds'] for item in records),
        'CPU_cap': CPU_CAP,
        'family_CPU_cap': FAMILY_CAP,
        'new_operator_solves': sum(item['new_operator_solves'] for item in records),
        'reused_operator_solves': 0,
        'reused_n16_matter_columns': True,
        'source_hashes': ctx['hashes'],
        'scope': {
            'local_interval': 'S(1)+[.12,.18]',
            'physical_local_gate': 'OPEN',
            'physical_root_claimed': False,
            'physical_EXISTENCE_certificate': False,
            'physical_NONEXISTENCE_certificate': False,
            'full_source_coverage': complete,
            'source_cutoff_edge_and_tangent_included': complete,
            'baseline_and_geometry_counted_once': True,
            'changed_history_UV_bound': None,
            'between_node_error_bound': None,
            'field_error_bound': None,
            'low_subgap_error_bound': None,
            'metric_timestep': False,
            'new_action_term': False,
            'Gamma_rest_assigned': False,
        },
        'reproducer': 'python3 scripts/derive_nsc_ks_coupled_newton_n32_assembled.py --check',
    }


def run(cpu_budget, max_new):
    if not np.isfinite(cpu_budget) or cpu_budget <= 0 or max_new < 0:
        raise ValueError('nonnegative work bounds required')
    ctx = context()
    DIRECTORY.mkdir(parents=True, exist_ok=True)
    start = time.process_time()
    new = 0
    spent = 0.
    for key in sorted(ctx['archived']):
        path, payload = family_paths(key)
        if path.exists():
            record, _, _, _ = read_family(ctx, key)
            spent += record['CPU_seconds']
            continue
        if payload.exists():
            raise ValueError('orphan operator payload retained for recovery: '+str(payload))
        if new >= max_new:
            continue
        remaining = min(cpu_budget - (time.process_time()-start), CPU_CAP - spent)
        if remaining <= 1.:
            break
        try:
            record = new_family(ctx, key, cpu_limit=remaining-1.)
        except TimeoutError as error:
            print(json.dumps({'status': 'OPEN', 'timeout': str(error),
                              'positive_family': list(key),
                              'physical_NONEXISTENCE_certificate': False}), flush=True)
            break
        new += record['new_operator_solves']
        spent += record['CPU_seconds']
        print(json.dumps({k: record[k] for k in
                          ('positive_family', 'source_rows', 'matter_change_maxima',
                           'high_matter_tangent_maxima', 'CPU_seconds',
                           'new_operator_solves')}), flush=True)
        progress = assemble(ctx)
        (DIRECTORY/'progress.json').write_text(R.dump(progress))
    result = assemble(ctx)
    if result['all_retained_nonzero_angular_families_included']:
        progress = DIRECTORY/'progress.json'
        if progress.exists():
            progress.unlink()
        OUTPUT.write_text(R.dump(result))
    else:
        (DIRECTORY/'progress.json').write_text(R.dump(result))
    result['runtime'] = {
        'CPU_seconds': time.process_time()-start,
        'CPU_cap': CPU_CAP,
        'family_CPU_seconds': result['CPU_seconds'],
        'new_operator_solves': result['new_operator_solves'],
        'reused_operator_solves': 0,
    }
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--run', action='store_true')
    mode.add_argument('--check', action='store_true')
    parser.add_argument('--cpu-budget', type=float, default=CPU_CAP)
    parser.add_argument('--max-new', type=int, default=MAX_NEW)
    args = parser.parse_args()
    if args.run:
        result = run(args.cpu_budget, args.max_new)
        runtime = result.get('runtime')
    else:
        ctx = context()
        result = assemble(ctx, replay=True)
        path = OUTPUT if OUTPUT.exists() else DIRECTORY/'progress.json'
        previous = json.loads(path.read_text())
        previous.pop('runtime', None)
        if R.jsonable(result) != previous:
            raise ValueError('Newton n=32 assembled retained-source replay differs')
        runtime = {'CPU_seconds': None, 'new_operator_solves': 0,
                   'reused_operator_solves': 0, 'replay_dirac_evolution': False}
    print(json.dumps({key: result[key] for key in (
        'status', 'completed_positive_families', 'required_positive_families',
        'constraint_maxima', 'assembled_unclipped_leftover',
        'assembled_clipped_leftover', 'assembled_condition_number',
        'evolve_authorized', 'families_evolved', 'named_search_stall',
        'CPU_seconds', 'new_operator_solves') if key in result}, indent=2))
    print(json.dumps(runtime, indent=2))
