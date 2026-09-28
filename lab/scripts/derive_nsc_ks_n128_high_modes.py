#!/usr/bin/env python3
"""Chebyshev modes 64..127 of the declared (w, U) class at ad759424.

The n=64 Jacobian is reused for modes 0..63. This owner evolves only the
new modes, under primal-block DOP853 control, and concatenates them in
family order [w0:128, U0:128]. A direction is not evolved as a Newton step
unless its linear N gain can clear 1e-10 while beta moves. No s^5 term.
"""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
import argparse
from dataclasses import dataclass
from hashlib import sha256
import json
from pathlib import Path
import signal
import sys
import time

import numpy as np
from scipy.optimize import linprog

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
sys.path.insert(0, str(ROOT/'scripts'))
import derive_nsc_ks_coupled_newton_n16 as N16
import derive_nsc_ks_coupled_retained_control as R
import recursive_horizons.nsc_ks_energy_propagator as P
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
from recursive_horizons.nsc_ks_primal_step_control import evolve_primal_controlled
from recursive_horizons.nsc_ks_profile_identity import profile_identity
from recursive_horizons.nsc_ks_signed_state import negative_angular_partner
from recursive_horizons.nsc_ks_source_envelope import (
    KSEnvelopeBinding, usual_axial_support, _sample_axial_profiles)
from recursive_horizons.nsc_ks_spacetime_variation import chart_coordinates
from recursive_horizons.nsc_local_incoming_constraints import _channel_record
from recursive_horizons.nsc_local_incoming_family import (
    LocalAxialFunction, LocalIncomingFamily, plateau_derivatives)
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes

import derive_nsc_ks_coupled_newton_n16_damped as D

DIRECTORY = ROOT/'results/development/artifacts/nsc-ks-n128-high-modes'
PARENT_ARTIFACTS = ROOT/'results/development/artifacts/nsc-ks-coupled-newton-n64-primal-norm-trust03'
PARENT = ROOT/'results/development/nsc-ks-coupled-newton-n64-primal-norm-trust03.json'
OUTPUT = ROOT/'results/development/nsc-ks-n128-high-modes.json'
MODES = 128
LOW = 64
HIGH = 128
FROZEN_PLUS = 129
FAMILY_CAP = 240.0
DRIFT_ABS = 1e-12
EXISTENCE = 3e-11
BAR = 1e-10
PARENT_IDENTITY = 'ad75942438dad1de599ee260e5d9aa7f17d60af26b533bee2b0691f2a6dec777'
LIFTED_IDENTITY = '3af275b0ee0d0d0237799290e7a76639513d47cc5f1bc6da1d4d33593392f888'
SOLVER = 'primal-block DOP853 RMS on (A, D); tangent directions excluded'
OWNERS = (
    'scripts/derive_nsc_ks_n128_high_modes.py',
    'docs/nsc-ks-n128-high-modes.md',
    'src/recursive_horizons/nsc_ks_primal_step_control.py',
)
LP_OPTIONS = {
    'primal_feasibility_tolerance': 1e-10,
    'dual_feasibility_tolerance': 1e-10,
}


@dataclass(frozen=True)
class LocalIncomingFamily128(LocalIncomingFamily):
    """Same (w, U) class at 128 Chebyshev coefficients per function.

    The base allowlist stays 8, 16, 32, 64 so historical source hashes of
    that file remain valid. This subclass is the n=128 owner.
    """

    def __post_init__(self):
        raw = np.asarray(self.coefficients)
        if (raw.ndim != 2 or raw.shape != (2, MODES)
                or np.iscomplexobj(raw) or not np.isfinite(raw).all()):
            raise ValueError('two finite real vectors with 128 coefficients required')
        center = chart_coordinates(1.)[1]+.15 if self.center is None else float(self.center)
        plateau_derivatives(center, center, self.axial_inner, self.axial_outer)
        if (not np.isfinite([self.normal_inner, self.normal_outer]).all()
                or not 0 < self.normal_inner < self.normal_outer):
            raise ValueError('positive ordered normal support radii required')
        object.__setattr__(self, 'center', center)
        object.__setattr__(self, 'coefficients', tuple(tuple(map(float, row)) for row in raw))


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


def use_primal_control():
    P.evolve_ks_difference_envelope = evolve_primal_controlled


def lifted_coefficients(parent):
    coeff = np.asarray(parent['history']['coefficients'], float)
    if coeff.shape != (2, LOW):
        raise ValueError('parent history must be 64 coefficients per function')
    padded = np.zeros((2, MODES), float)
    padded[:, :LOW] = coeff
    if np.max(np.abs(padded[:, LOW:])) != 0.0:
        raise ValueError('zero pad must leave modes 64..127 unused')
    return padded


def next_metric(family):
    coeff = np.asarray(family.coefficients, float)
    center, ai, ao = family.center, family.axial_inner, family.axial_outer
    ni, no = family.normal_inner, family.normal_outer
    zero = LocalAxialFunction((0.,), center, ai, ao)
    frozen = CompatibleRadiusDirection(
        LocalAxialFunction(tuple(coeff[0]), center, ai, ao),
        LocalAxialFunction(tuple(coeff[1]), center, ai, ao), ni, no)
    eye = np.eye(MODES)
    basis = [LocalAxialFunction(tuple(eye[i]), center, ai, ao) for i in range(LOW, MODES)]
    high_w = tuple(CompatibleRadiusDirection(basis[i], zero, ni, no) for i in range(HIGH // 2))
    high_u = tuple(CompatibleRadiusDirection(zero, basis[i], ni, no) for i in range(HIGH // 2))
    return CompatibleIncomingMetric((1.,) + (0.,) * HIGH, (frozen,) + high_w + high_u)


def context():
    use_primal_control()
    base = R.context()
    parent = json.loads(PARENT.read_text())
    if parent['profile_identity'] != PARENT_IDENTITY or not parent['newton_step_accepted']:
        raise ValueError('n=128 parent must be the accepted ad759424 history')
    if parent['solver'] != SOLVER:
        raise ValueError('n=128 parent Jacobian must come from primal-block control')
    padded = lifted_coefficients(parent)
    family = LocalIncomingFamily128(padded)
    if family.radius_lower_bound() <= 0:
        raise ValueError('lifted history must keep a positive radius')
    identity = profile_identity(family, include_normal_window=True)
    if identity != LIFTED_IDENTITY:
        raise ValueError('lifted profile identity changed')
    if identity in D.FORBIDDEN_IDENTITIES or identity == PARENT_IDENTITY:
        raise ValueError('lifted identity is forbidden or is not a new owner')
    metric = next_metric(family)
    if len(metric.directions) != FROZEN_PLUS:
        raise ValueError('frozen history plus modes 64..127 required')
    solve = family.collocation_nodes(N16.SOLVE_COUNT)
    verify = family.collocation_nodes(N16.VERIFY_COUNT)
    target = H.fixed_target_union(solve, verify)
    if target.tolist() != parent['z']:
        raise ValueError('n=128 modes must keep the ad759424 nodes')
    slots, tangents = compatible_history_slots(metric, target, FROZEN_PLUS)
    owned, owned_tangents = compatible_history_slots(family.metric(), target, 2 * MODES)
    if owned_tangents.shape[0] != 2 * MODES:
        raise ValueError('n=128 family must carry 256 directions')
    if np.max(np.abs(owned_tangents[LOW:MODES] - tangents[1:1 + LOW])) != 0.0:
        raise ValueError('high w slot tangents left the n=128 family basis')
    if np.max(np.abs(owned_tangents[MODES + LOW:] - tangents[1 + LOW:])) != 0.0:
        raise ValueError('high U slot tangents left the n=128 family basis')
    if np.max(np.abs(slots - owned)) > 1e-11:
        raise ValueError('frozen history slots left the lifted family')
    return {
        **base, 'family': family, 'metric': metric, 'target': target,
        'identity': identity, 'parent': parent,
        'evolution_identity': profile_identity(metric, include_normal_window=True),
        'hashes': {**{path: digest(path) for path in OWNERS},
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
            raise ValueError('saved n=128 operator digest changed')
    if operator.tangent.shape[1] != FROZEN_PLUS:
        raise ValueError('frozen-plus-high-mode tangents required')
    corrections, tangents, family_records, states = {}, {}, {}, {}

    def record_signed(signed_id, matter, channel, energy_sign, batch_ids):
        column = np.asarray(matter['tangent'][1:], float)
        if not np.isfinite(column).all() or not np.isfinite(matter['correction']).all():
            raise ValueError('high-mode matter tangent is not finite')
        if signed_id in corrections:
            corrections[signed_id] = corrections[signed_id] + matter['correction']
            tangents[signed_id] = tangents[signed_id] + column
            family_records[signed_id]['batch_ids'].extend(batch_ids)
            return
        corrections[signed_id] = matter['correction']
        tangents[signed_id] = column
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
        raise FileExistsError('n=128 family output already exists')
    cpu = time.process_time()
    limit = min(FAMILY_CAP, float(cpu_limit))

    def stop(*_):
        raise TimeoutError('one n=128 mode operator exceeded its CPU allocation')

    previous = signal.signal(signal.SIGPROF, stop)
    signal.setitimer(signal.ITIMER_PROF, limit)
    try:
        arrays, meta, signed, operator = evaluate_family(ctx, key, saved=None, allow_new=True)
    finally:
        signal.setitimer(signal.ITIMER_PROF, 0)
        signal.signal(signal.SIGPROF, previous)
    _parent_json, parent_npz = parent_family_paths(key)
    with np.load(parent_npz, allow_pickle=False) as handle:
        parent_change = np.asarray(handle['paired_matter_change'], float)
    evolved = np.asarray(arrays['paired_matter_change'], float)
    drift = float(np.max(np.abs(evolved - parent_change)))
    if drift > DRIFT_ABS:
        raise ValueError(f'frozen ad759424 matter change drifted by {drift}')
    if arrays['paired_high_matter_tangent'].shape != (HIGH, len(ctx['target']), 2):
        raise ValueError('high-mode tangent shape changed')
    raw = deterministic_npz_bytes(arrays)
    record = {
        'positive_family': list(key),
        'operator_digest': operator.digest,
        'signed_families': signed,
        'family_records': meta['family_records'],
        'matter_absolute_drift': drift,
        'CPU_seconds': time.process_time() - cpu,
        'accepted_steps': int(operator.diagnostics['accepted_steps']),
        'new_operator_solves': 1,
        'solver': SOLVER,
        'source_hashes': ctx['hashes'],
        'profile_identity': ctx['identity'],
        'parent_profile_identity': PARENT_IDENTITY,
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
        raise ValueError('n=128 family payload changed')
    if record['source_hashes'] != ctx['hashes'] or record['profile_identity'] != ctx['identity']:
        raise ValueError('n=128 family belongs to a different history')
    if record['solver'] != SOLVER or record['parent_profile_identity'] != PARENT_IDENTITY:
        raise ValueError('n=128 family left primal-block control at ad759424')
    with np.load(payload, allow_pickle=False) as handle:
        tangent = np.asarray(handle['paired_high_matter_tangent'], float)
    if tangent.shape != (HIGH, len(ctx['target']), 2):
        raise ValueError('high-mode tangent shape changed')
    return tangent, record


def family_order(low, extra):
    """[w0:64, U0:64] plus [w64:128, U64:128] becomes [w0:128, U0:128]."""
    if low.shape[0] != 2 * LOW or extra.shape != (HIGH, low.shape[1], 2):
        raise ValueError('low and high blocks do not form an n=128 Jacobian')
    half = LOW
    w = np.concatenate([low[:half], extra[:half]], axis=0)
    u = np.concatenate([low[half:], extra[half:]], axis=0)
    return np.concatenate([w, u], axis=0)


def predicted_maxima(residual, matrix, step):
    prediction = np.asarray(residual, float) + np.asarray(matrix, float) @ np.asarray(step, float)
    nodes = residual.size // 2
    return np.array([
        float(np.max(np.abs(prediction[:nodes]))),
        float(np.max(np.abs(prediction[nodes:]))),
    ])


def feasible(response, residual, n_cap, beta_cap, h):
    nodes = residual.size // 2
    width = response.shape[1]
    scale = 1e6 / max(float(np.max(np.abs(residual))), 1e-30)
    rn = residual[:nodes] * scale
    rb = residual[nodes:] * scale
    an = response[:nodes] * (h * scale)
    ab = response[nodes:] * (h * scale)
    caps = np.array([n_cap, beta_cap], float) * scale
    cuts = []
    objective = np.concatenate([np.zeros(width), np.ones(width)])
    eye = np.eye(width)
    for _cut in range(16):
        gain_rows = np.vstack([an, -an, ab, -ab])
        gain_limits = np.concatenate([
            caps[0] - rn, caps[0] + rn, caps[1] - rb, caps[1] + rb,
        ])
        if cuts:
            gain_rows = np.vstack([gain_rows, np.vstack(cuts)])
            gain_limits = np.concatenate([gain_limits, np.full(len(cuts), 1.0)])
        pad = np.zeros((gain_rows.shape[0], width))
        augmented = np.vstack([
            np.hstack([gain_rows, pad]),
            np.hstack([eye, -eye]),
            np.hstack([-eye, -eye]),
        ])
        bounds = np.concatenate([gain_limits, np.zeros(2 * width)])
        solved = linprog(
            objective, A_ub=augmented, b_ub=bounds, bounds=(None, None),
            method='highs', options=LP_OPTIONS)
        if not solved.success:
            return False, None
        direction = solved.x[:width]
        if float(np.linalg.norm(direction)) <= 1.0 + 1e-7:
            return True, direction
        cuts.append(direction / np.linalg.norm(direction))
    return False, None


def max_n_gain(response, residual, measured, h, beta_drop):
    lo = 0.0
    hi = float(measured[0])
    best = None
    for _step in range(24):
        mid = 0.5 * (lo + hi)
        ok, direction = feasible(
            response, residual, measured[0] - mid, measured[1] - beta_drop, h)
        if ok:
            lo = mid
            best = direction
        else:
            hi = mid
    return lo, best


def whiten(right, col, start, end):
    raw = right[start:end].T / col[:, None]
    gram = raw.T @ raw
    evals, evecs = np.linalg.eigh(gram)
    keep = evals > evals[-1] * 1e-10
    return raw @ (evecs[:, keep] / np.sqrt(evals[keep]))


def ray_ceiling(step, matrix, residual, measured):
    nodes = residual.size // 2
    step = np.asarray(step, float)
    step = step / max(float(np.linalg.norm(step)), 1e-30)
    response = matrix @ step
    best_gain = 0.0
    best_h = 0.0
    best_beta = 0.0
    for sign in (1.0, -1.0):
        column = sign * response
        for h in np.geomspace(1e-12, 1e-3, 25):
            prediction = residual + float(h) * column
            maxima = np.array([
                float(np.max(np.abs(prediction[:nodes]))),
                float(np.max(np.abs(prediction[nodes:]))),
            ])
            gain = measured - maxima
            if gain[1] >= -1e-18 and gain[0] > best_gain:
                best_gain = float(gain[0])
                best_h = float(h)
                best_beta = float(gain[1])
    slope = best_gain / best_h if best_h else 0.0
    return best_gain, best_h, slope, best_beta


def window_row(name, basis, matrix, residual, measured, h):
    response = matrix @ basis
    gain, direction = max_n_gain(response, residual, measured, h, 0.0)
    beta_gain = 0.0
    if direction is not None:
        step = h * (basis @ direction)
        predicted = predicted_maxima(residual, matrix, step)
        beta_gain = float(measured[1] - predicted[1])
    slope = gain / h if h else 0.0
    return {
        'window': name,
        'width': int(basis.shape[1]),
        'h': float(h),
        'n_gain': float(gain),
        'beta_gain': beta_gain,
        'slope': float(slope),
        'quad_coeff': None,
        'crossover_n_gain_upper_bound': float(gain),
        'clears_bar_with_beta_moving': bool(gain >= BAR and beta_gain > 1e-16),
        'probed': False,
    }


def screen_new_columns(jacobian, gradient, measured):
    residual = np.concatenate([gradient[:, 0], gradient[:, 1]])
    matrix = np.vstack([jacobian[:, :, 0].T, jacobian[:, :, 1].T])
    high = np.concatenate([np.arange(LOW, MODES), np.arange(MODES + LOW, 2 * MODES)])
    block = matrix[:, high]
    row = np.linalg.norm(block, axis=1)
    col = np.linalg.norm(block, axis=0)
    row = np.where(row > 0, row, 1.0)
    col = np.where(col > 0, col, 1.0)
    _left, spectrum, right = np.linalg.svd(block / row[:, None] / col[None, :], full_matrices=False)
    rows = []
    spans = (
        ('new_w_64_80', 0, 16, 1e-3),
        ('new_w_80_96', 16, 32, 1e-3),
        ('new_w_96_112', 32, 48, 1e-3),
        ('new_w_112_128', 48, 64, 1e-3),
        ('new_U_64_80', 64, 80, 1e-3),
        ('new_U_80_96', 80, 96, 1e-3),
        ('new_U_96_112', 96, 112, 1e-3),
        ('new_U_112_128', 112, 128, 1e-3),
        ('all_new', 0, 128, 1e-6),
        ('all_new', 0, 128, 1e-3),
        ('new_svd_16', 0, 16, 1e-3),
        ('new_svd_32', 0, 32, 1e-3),
        ('new_svd_64', 0, 64, 1e-3),
    )
    for name, start, end, h in spans:
        if name.startswith('new_svd'):
            basis = whiten(right, col, start, end)
        else:
            basis = np.eye(HIGH)[:, start:end]
        rows.append(window_row(name, basis, block, residual, measured, h))
    best_mode = {'index': 0, 'n_gain': 0.0, 'h': 0.0, 'slope': 0.0, 'beta_gain': 0.0}
    for index in range(right.shape[0]):
        gain, h, slope, beta = ray_ceiling(right[index] / col, block, residual, measured)
        if gain > best_mode['n_gain']:
            best_mode = {
                'index': int(index + 1), 'n_gain': gain, 'h': h,
                'slope': slope, 'beta_gain': beta,
            }
    best_coordinate = {'index': 0, 'n_gain': 0.0, 'h': 0.0, 'slope': 0.0, 'beta_gain': 0.0}
    for index in range(HIGH):
        step = np.zeros(HIGH)
        step[index] = 1.0
        gain, h, slope, beta = ray_ceiling(step, block, residual, measured)
        if gain > best_coordinate['n_gain']:
            best_coordinate = {
                'index': int(high[index]), 'n_gain': gain, 'h': h,
                'slope': slope, 'beta_gain': beta,
            }
    small = next(row for row in rows if row['window'] == 'all_new' and row['h'] == 1e-6)
    large = next(row for row in rows if row['window'] == 'all_new' and row['h'] == 1e-3)
    plateau = abs(small['n_gain'] - large['n_gain']) <= 1e-14 * max(1.0, abs(large['n_gain']))
    clears = any(row['clears_bar_with_beta_moving'] for row in rows)
    clears = clears or best_mode['n_gain'] >= BAR or best_coordinate['n_gain'] >= BAR
    ratio = float(spectrum[0] / spectrum[-1]) if spectrum[-1] > 0 else None
    return {
        'directions_screened': rows,
        'best_new_svd_mode': best_mode,
        'best_new_coordinate_ray': best_coordinate,
        'new_block_singular_ratio': ratio,
        'new_block_plateau': bool(plateau),
        'probe_authorized': bool(clears),
        'quad_coeff_measured': False,
    }


def assemble(ctx):
    records = []
    tangent = np.zeros((HIGH, len(ctx['target']), 2), float)
    family_records = {}
    drifts = []
    for key in sorted(ctx['archived']):
        path, _payload = family_paths(key)
        if not path.exists():
            continue
        column, record = load_tangent(ctx, key)
        tangent += column
        records.append(record)
        drifts.append(record['matter_absolute_drift'])
        family_records.update(record['family_records'])
    complete = len(records) == len(ctx['archived'])
    report = {
        'schema': 'NSC-KS-N128-HIGH-MODES-v1',
        'accountable_author': 'Douglas Ek',
        'completed_positive_families': len(records),
        'required_positive_families': len(ctx['archived']),
        'complete': complete,
        'solver': SOLVER,
        'profile_identity': ctx['identity'],
        'parent_profile_identity': PARENT_IDENTITY,
        'history_class': 'zero pad of ad759424 inside declared (w, U); modes 64..127 added; no s^5',
        'physical_EXISTENCE_certificate': False,
        'physical_NONEXISTENCE_certificate': False,
        'enclosed_error': None,
        'uv_field_low_subgap_between_node': None,
        'families_evolved_as_newton_steps': 0,
    }
    if not complete:
        report['status'] = 'OPEN: partial n=128 high-mode columns'
        return report
    if len(family_records) != 120:
        raise ValueError('one hundred twenty explicit signed families required')
    slots, slot_tangents = compatible_history_slots(ctx['metric'], ctx['target'], FROZEN_PLUS)
    geometry = surface_geometry_response(slots, slot_tangents, ctx['coeff'])
    weights = signed_family_angular_square_weights(family_records)
    first = ctx['archive'].family_entries(next(iter(ctx['archived'])))[0][0]
    phase = formal_source_phase_coefficient(
        ctx['metric'], ctx['target'], angular=1., rho_up=first.rho_up, gauss_nodes=N16.PHASE_NODES)
    if phase.delta_f_z.shape[1] != FROZEN_PLUS:
        raise ValueError('phase must retain the frozen history and modes 64..127')
    _edge, edge_tangent = source_edge_from_phase(phase, weights, ctx['coeff']['a'])
    extra = geometry['action_gradient_tangent'][1:] + tangent + edge_tangent[1:]
    if not np.isfinite(extra).all():
        raise ValueError('assembled high-mode columns are not finite')
    low = np.asarray(ctx['parent']['history_jacobian'], float)
    gradient = np.asarray(ctx['parent']['action_gradient'], float)
    full = family_order(low, extra)
    measured = np.asarray(ctx['parent']['constraint_maxima'], float)
    screened = screen_new_columns(full, gradient, measured)
    n_gap = float(measured[0] - EXISTENCE)
    beta_gap = float(measured[1] - EXISTENCE)
    if screened['probe_authorized']:
        status = 'OPEN: an n=128 high-mode direction clears the 1e-10 linear N bar; quadratic probe not yet measured'
        stall = None
        next_owner = 'measure that direction\'s own quadratic N coefficient with a small primal-controlled probe before any accepted step'
    else:
        status = 'OPEN: n=128 high-mode linear N gain stays below 1e-10'
        stall = 'n128_high_mode_linear_gain_below_1e-10'
        next_owner = (
            'declared (w, U) class on I; n=128 columns do not clear 1e-10. '
            'A further doubling to 256 coefficients is the same subclass pattern, not a new residual class, and is not authorized as a Newton step from this table'
        )
    report.update({
        'status': status,
        'stall': stall,
        'best_constraint_maxima': measured.tolist(),
        'n_gap_to_existence': n_gap,
        'beta_gap_to_existence': beta_gap,
        'reused_parent_columns': LOW * 2,
        'new_columns': HIGH,
        'matter_absolute_drift_max': float(np.max(drifts)),
        'history_jacobian': full.tolist(),
        'radius_lower_bound': ctx['family'].radius_lower_bound(),
        'evolution_metric_identity': ctx['evolution_identity'],
        'newton_step_accepted': False,
        'probe_evolved': 0,
        'evolve_authorized': False,
        'nonexistence_not_claimed': 'joint linear gain on the reused n=64 block stays positive and the geometry cokernel stays the recorded trivial one',
        'higher_n_pathway': 'n=256 would be another LocalIncomingFamily subclass; the base allowlist remains 8, 16, 32, 64',
        'next_owner': next_owner,
        'source_hashes': ctx['hashes'],
        **screened,
    })
    return report


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
            'accepted_steps': record['accepted_steps'],
            'matter_absolute_drift': record['matter_absolute_drift'],
        }), flush=True)
    result = assemble(ctx)
    target = OUTPUT if result.get('complete') else DIRECTORY/'progress.json'
    target.write_text(json.dumps(R.jsonable(result), indent=2, sort_keys=True)+'\n')
    if result.get('complete'):
        progress = DIRECTORY/'progress.json'
        if progress.exists():
            progress.unlink()
    return result


def replay_family(ctx, key):
    path, payload = family_paths(key)
    record = json.loads(path.read_text())
    with np.load(payload, allow_pickle=False) as handle:
        arrays = {name: np.array(handle[name], copy=True) for name in handle.files}
    meta = json.loads(arrays['metadata_json'].tobytes())
    current, _meta, signed, operator = evaluate_family(
        ctx, key, saved=(arrays, meta), allow_new=False)
    if signed != record['signed_families']:
        raise ValueError('n=128 replay signed coverage changed')
    if operator.digest != record['operator_digest']:
        raise ValueError('n=128 replay operator digest changed')
    if not np.array_equal(current['paired_high_matter_tangent'], arrays['paired_high_matter_tangent']):
        raise ValueError('n=128 replay tangent changed')
    if not np.array_equal(current['paired_matter_change'], arrays['paired_matter_change']):
        raise ValueError('n=128 replay matter change changed')
    return record


def check():
    ctx = context()
    saved = json.loads(OUTPUT.read_text())
    replay_family(ctx, (1, -1))
    result = assemble(ctx)
    for key in (
        'stall', 'profile_identity', 'best_constraint_maxima',
        'n_gap_to_existence', 'beta_gap_to_existence',
        'directions_screened', 'best_new_svd_mode', 'best_new_coordinate_ray',
        'probe_authorized', 'matter_absolute_drift_max',
        'physical_NONEXISTENCE_certificate', 'history_jacobian',
    ):
        if saved.get(key) != result.get(key):
            raise ValueError('n=128 register changed: '+key)
    if saved['newton_step_accepted'] or saved['physical_EXISTENCE_certificate']:
        raise ValueError('n=128 column assembly is not a gate certificate')
    return saved


def preflight():
    ctx = context()
    parent = ctx['parent']
    return {
        'profile_identity': ctx['identity'],
        'parent_profile_identity': PARENT_IDENTITY,
        'radius_lower_bound': ctx['family'].radius_lower_bound(),
        'directions': len(ctx['metric'].directions),
        'nodes': len(ctx['target']),
        'parent_maxima': parent['constraint_maxima'],
        'evolution_metric_identity': ctx['evolution_identity'],
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--run', action='store_true')
    mode.add_argument('--check', action='store_true')
    mode.add_argument('--preflight', action='store_true')
    parser.add_argument('--cpu-budget', type=float, default=14400.0)
    parser.add_argument('--max-new', type=int, default=60)
    args = parser.parse_args()
    if args.preflight:
        print(json.dumps(preflight(), indent=2))
    elif args.check:
        saved = check()
        print(json.dumps({
            'stall': saved['stall'],
            'status': saved['status'],
            'profile_identity': saved['profile_identity'],
            'probe_authorized': saved['probe_authorized'],
            'matter_absolute_drift_max': saved['matter_absolute_drift_max'],
        }, indent=2))
    else:
        result = run(args.cpu_budget, args.max_new)
        print(json.dumps({
            'completed_positive_families': result.get('completed_positive_families'),
            'required_positive_families': result.get('required_positive_families'),
            'status': result.get('status'),
            'stall': result.get('stall'),
            'probe_authorized': result.get('probe_authorized'),
        }, indent=2))
