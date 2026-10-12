#!/usr/bin/env python3
"""Dirac-free n=32 geometry leftover at iterate4 g. No high-mode matter."""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
sys.path.insert(0, str(ROOT/'scripts'))
import derive_nsc_ks_coupled_newton_n16_damped as D
import derive_nsc_ks_coupled_newton_n16_damped_iterate4 as I4
import derive_nsc_ks_n16_geometry_killing_current as K
import derive_nsc_ks_n16_leftover_anatomy as A
from recursive_horizons.nsc_evolved_incoming_constraints import (
    compatible_history_slots, surface_geometry_response)
from recursive_horizons.nsc_ks_profile_identity import profile_identity
from recursive_horizons.nsc_local_history_newton import (
    HistoryCollocation, LocalHistoryNewtonSettings, _clip_step, _linear_step)
from recursive_horizons.nsc_local_incoming_family import LocalIncomingFamily

OUTPUT = ROOT/'results/development/nsc-ks-n32-geometry-leftover.json'
OWNERS = (
    'scripts/derive_nsc_ks_n32_geometry_leftover.py',
    'docs/nsc-ks-n32-geometry-leftover.md',
    'src/recursive_horizons/nsc_local_incoming_family.py',
)


def digest(path):
    path = Path(path)
    return sha256((path if path.is_absolute() else ROOT/path).read_bytes()).hexdigest()


def _round(value, digits=12):
    array = np.asarray(value, float)
    if array.shape == ():
        return None if not np.isfinite(array) else round(float(array), digits)
    return [None if not np.isfinite(item) else round(float(item), digits)
            for item in array.ravel()]


def zero_pad(family):
    coeff = np.asarray(family.coefficients, float)
    padded = np.zeros((2, 32), float)
    padded[:, :coeff.shape[1]] = coeff
    return LocalIncomingFamily(padded)


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
    n_n = n_nodes
    return {
        'unclipped': {
            'N_max': float(np.max(abs(pred[:n_n]))),
            'beta_max': float(np.max(abs(pred[n_n:]))),
            'control_max': float(np.max(abs(pred))),
            'step_norm': step_norm,
        },
        'clipped': {
            'N_max': float(np.max(abs(clipped_pred[:n_n]))),
            'beta_max': float(np.max(abs(clipped_pred[n_n:]))),
            'control_max': float(np.max(abs(clipped_pred))),
            'step_norm': clip_norm,
        },
        'condition_number': float(diagnostics['condition_number']),
        'singular': np.asarray(diagnostics['singular_values'], float),
        'rank_or_conditioning': bool(diagnostics['rank_or_conditioning']),
        'delta': delta,
    }


def truncated(residual, jacobian, singular, n_nodes, settings):
    factors = A.scaled_factors(jacobian, residual, jacobian, settings)
    rows = []
    for cutoff in A.CUTOFFS:
        keep = factors['singular'] > (float(factors['singular'][0]) / cutoff)
        delta = A.delta_from_keep(factors, keep)
        pred = residual if delta is None else residual + jacobian @ delta
        rows.append({
            'condition_cutoff': cutoff,
            'kept_modes': int(np.count_nonzero(keep)),
            'N_max': _round(np.max(abs(pred[:n_nodes]))),
            'beta_max': _round(np.max(abs(pred[n_nodes:]))),
            'control_max': _round(np.max(abs(pred))),
            'step_norm': None if delta is None else _round(np.linalg.norm(delta)),
        })
    return rows


def compute():
    killing = json.loads(K.OUTPUT.read_text())
    trial = json.loads(I4.OUTPUT.read_text())
    for path, expected in {**killing['source_hashes'], **killing['input_hashes'],
                           **trial['source_hashes']}.items():
        if digest(path) != expected:
            raise ValueError('n=32 geometry leftover input changed: '+path)
    if killing['physical_NONEXISTENCE_certificate']:
        raise ValueError('n=32 leftover is only for a failed NON-EXISTENCE branch')
    ctx = I4.context()
    family16 = ctx['family']
    family32 = zero_pad(family16)
    if family32.radius_lower_bound() <= 0:
        raise ValueError('zero-padded n=32 family must keep a positive radius bound')
    z = np.asarray(trial['z'], float)
    gradient = np.asarray(trial['action_gradient'], float)
    slots, ds = compatible_history_slots(family32.metric(), z, 64)
    geometry = surface_geometry_response(slots, ds, ctx['coeff'])
    settings = LocalHistoryNewtonSettings()
    all_node = leftover_of(gradient, geometry['action_gradient_tangent'], settings)
    collocation = HistoryCollocation(16, layout='rectangular')
    solve_nodes = collocation.solve_nodes(family16)
    indices = np.searchsorted(z, solve_nodes)
    if not np.array_equal(z[indices], solve_nodes):
        raise ValueError('exact owned n=16 solve-node subset required')
    solve = leftover_of(
        gradient[indices], geometry['action_gradient_tangent'][:, indices], settings)
    residual, jacobian, n_nodes = flatten_system(
        gradient, geometry['action_gradient_tangent'])
    truncated_rows = truncated(
        residual, jacobian, all_node['singular'], n_nodes, settings)
    leftover = np.array([
        all_node['unclipped']['N_max'], all_node['unclipped']['beta_max']], float)
    measured = np.asarray(trial['constraint_maxima'], float)
    anatomy = json.loads(A.OUTPUT.read_text())
    n16_leftover = np.asarray([
        anatomy['all_node_unclipped_leftover_from_solve_step']['N_max'],
        anatomy['all_node_unclipped_leftover_from_solve_step']['beta_max']], float)
    reaches = bool(np.all(leftover <= A.EXISTENCE_TOLERANCE))
    order_below = bool(np.all(leftover <= 0.1 * n16_leftover))
    clipped_beats = bool(np.all(np.array([
        all_node['clipped']['N_max'], all_node['clipped']['beta_max']]) < measured))
    clip_norm = all_node['clipped']['step_norm']
    evolve_authorized = bool(
        all_node['delta'] is not None
        and (reaches or order_below) and clipped_beats
        and clip_norm is not None and clip_norm <= 1.0 + 1e-12)
    return {
        'schema': 'NSC-KS-N32-GEOMETRY-LEFTOVER-v1',
        'accountable_author': 'Douglas Ek',
        'status': (
            'OPEN: n=32 geometry leftover cannot reach EXISTENCE; '
            'no family evolved'
            if not evolve_authorized else
            'n=32 geometry leftover authorizes one retarded evolve'),
        'current_profile_identity': profile_identity(family16),
        'padded_profile_identity': profile_identity(family32),
        'n32_coefficient_count': 32,
        'n32_direction_count': 64,
        'all_node_unclipped_leftover': _round(leftover),
        'all_node_clipped_leftover': _round([
            all_node['clipped']['N_max'], all_node['clipped']['beta_max']]),
        'all_node_unclipped_step_norm': _round(all_node['unclipped']['step_norm']),
        'all_node_clipped_step_norm': _round(all_node['clipped']['step_norm']),
        'all_node_condition_number': _round(all_node['condition_number']),
        'all_node_rank_or_conditioning': all_node['rank_or_conditioning'],
        'solve_unclipped_leftover': _round([
            solve['unclipped']['N_max'], solve['unclipped']['beta_max']]),
        'solve_unclipped_step_norm': _round(solve['unclipped']['step_norm']),
        'solve_condition_number': _round(solve['condition_number']),
        'truncated_svd_all_node': truncated_rows,
        'n16_all_node_leftover': _round(n16_leftover),
        'unclipped_reaches_existence': reaches,
        'unclipped_is_order_below_n16_leftover': order_below,
        'clipped_beats_measured': clipped_beats,
        'evolve_authorized': evolve_authorized,
        'existence_tolerance': A.EXISTENCE_TOLERANCE,
        'forbidden_identities': list(D.FORBIDDEN_IDENTITIES),
        'named_search_stall': (
            None if evolve_authorized else 'n32_geometry_leftover_cannot_reach_existence'),
        'scope': {
            'physical_local_gate': 'OPEN',
            'prediction_is_not_a_nonlinear_residual': True,
            'high_mode_matter_columns_invented': False,
            'new_field_or_source_evolutions': 0,
            'physical_EXISTENCE_certificate': False,
            'physical_NONEXISTENCE_claimed': False,
            'metric_timestep': False,
            'full_scale_n16_not_evolved': True,
            'n8_next_not_evolved': True,
        },
        'source_hashes': {p: digest(p) for p in OWNERS},
        'input_hashes': {
            str(K.OUTPUT.relative_to(ROOT)): digest(K.OUTPUT),
            str(I4.OUTPUT.relative_to(ROOT)): digest(I4.OUTPUT),
            str(A.OUTPUT.relative_to(ROOT)): digest(A.OUTPUT),
        },
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--record', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = compute()
    if args.record:
        if OUTPUT.exists():
            raise FileExistsError('n=32 geometry leftover exists; use --check')
        OUTPUT.write_text(json.dumps(result, sort_keys=True, indent=2, allow_nan=False)+'\n')
    elif result != json.loads(OUTPUT.read_text()):
        raise ValueError('n=32 geometry leftover replay differs')
    print(json.dumps({
        'status': result['status'],
        'evolve_authorized': result['evolve_authorized'],
        'all_node_unclipped_leftover': result['all_node_unclipped_leftover'],
        'solve_unclipped_leftover': result['solve_unclipped_leftover'],
        'all_node_condition_number': result['all_node_condition_number'],
        'named_search_stall': result['named_search_stall'],
    }, indent=2))
