#!/usr/bin/env python3
"""Refuse NON-EXISTENCE and measure the dropped-jet all-node leftover."""
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
import derive_nsc_ks_n16_leftover_anatomy as A
from recursive_horizons.nsc_evolved_incoming_constraints import (
    compatible_history_slots, surface_geometry_response)
from recursive_horizons.nsc_local_history_newton import (
    DeclaredJetBoundary, HistoryCollocation, LocalHistoryNewtonSettings,
    _clip_step, _control_system, _linear_step)
from recursive_horizons.nsc_local_incoming_family import LocalIncomingFamily

OUTPUT = ROOT/'results/development/nsc-ks-n16-leftover-identification.json'
OWNERS = (
    'scripts/derive_nsc_ks_n16_leftover_identification.py',
    'docs/nsc-ks-n16-leftover-identification.md',
    'src/recursive_horizons/nsc_local_history_newton.py',
)


def digest(path):
    path = Path(path)
    return sha256((path if path.is_absolute() else ROOT/path).read_bytes()).hexdigest()


def _round(value, digits=12):
    array = np.asarray(value, float)
    if array.shape == ():
        return None if not np.isfinite(array) else round(float(array), digits)
    return [None if not np.isfinite(item) else round(float(item), digits) for item in array.ravel()]


def component_maxima(values):
    values = np.asarray(values, float)
    return [float(np.max(abs(values[:, 0]))), float(np.max(abs(values[:, 1])))]


def compute():
    anatomy = json.loads(A.OUTPUT.read_text())
    trial = json.loads(I4.OUTPUT.read_text())
    for path, expected in {**anatomy['source_hashes'], **anatomy['input_hashes'],
                           **trial['source_hashes']}.items():
        if digest(path) != expected:
            raise ValueError('leftover identification input changed: '+path)
    if anatomy['walk_projection_sign_stable']:
        raise ValueError('sign-stable leftover would require a different owner')
    ctx = I4.context()
    family = ctx['family']
    z = np.asarray(trial['z'], float)
    gradient = np.asarray(trial['action_gradient'], float)
    jacobian = np.asarray(trial['history_jacobian'], float)
    matter = np.asarray(trial['summed_matter_change'], float)
    edge = np.asarray(trial['source_cutoff_edge_gradient'], float)
    raw = np.asarray(trial['raw_action_gradient'], float)
    slots, ds = compatible_history_slots(family.metric(), z, 32)
    geometry = surface_geometry_response(slots, ds, ctx['coeff'])['action_gradient_change']
    baseline = raw - geometry - matter
    leftover_dir = gradient / np.linalg.norm(gradient)
    pieces = {
        'baseline': component_maxima(baseline),
        'geometry': component_maxima(geometry),
        'matter': component_maxima(matter),
        'edge': component_maxima(edge),
        'raw': component_maxima(raw),
        'assembled': component_maxima(gradient),
        'leftover_direction_dot_baseline': float(np.sum(leftover_dir * baseline)),
        'leftover_direction_dot_geometry': float(np.sum(leftover_dir * geometry)),
        'leftover_direction_dot_matter': float(np.sum(leftover_dir * matter)),
        'leftover_direction_dot_edge': float(np.sum(leftover_dir * edge)),
    }
    source_identity = json.loads(D.OUTPUT.read_text())['source_identity']
    settings = LocalHistoryNewtonSettings()
    collocation = HistoryCollocation(16, layout='rectangular')
    solve_nodes = collocation.solve_nodes(family)
    indices = np.searchsorted(z, solve_nodes)
    solve_slots, solve_ds = compatible_history_slots(family.metric(), solve_nodes, 32)
    solve_receipt = A.receipt_on(
        family, solve_nodes, gradient[indices], jacobian[:, indices],
        source_identity,
        surface_geometry_response(solve_slots, solve_ds, ctx['coeff'])['action_gradient_tangent'])
    solve_system = _control_system(
        solve_receipt, collocation,
        DeclaredJetBoundary(family.center, 0., ctx['seed_slope'], 0.))
    drop_r, drop_j, drop_p, drop_n, drop_b = A.dropped_jet_system(solve_system)
    drop_delta, drop_diag = _linear_step(drop_j, drop_r, drop_p, settings)
    if drop_delta is None:
        raise ValueError('dropped-jet rectangular step must remain available')
    clipped, clip_norm = _clip_step(drop_delta, settings.max_step_norm)
    unclipped_all = gradient + np.tensordot(drop_delta, jacobian, axes=(0, 0))
    clipped_all = gradient + np.tensordot(clipped, jacobian, axes=(0, 0))
    unclipped_max = component_maxima(unclipped_all)
    clipped_max = component_maxima(clipped_all)
    measured = np.asarray(trial['constraint_maxima'], float)
    unclipped_reaches = bool(np.all(np.asarray(unclipped_max) <= A.EXISTENCE_TOLERANCE))
    clipped_beats = bool(np.all(np.asarray(clipped_max) < measured))
    order_below = bool(np.all(np.asarray(clipped_max) <= 0.1 * measured))
    evolve_authorized = bool(
        (unclipped_reaches or order_below) and clipped_beats and clip_norm <= 1.0 + 1e-12)
    return {
        'schema': 'NSC-KS-N16-LEFTOVER-IDENTIFICATION-v1',
        'accountable_author': 'Douglas Ek',
        'status': (
            'OPEN: leftover direction is not sign-stable; '
            'no scoped NON-EXISTENCE certificate'),
        'profile_identity': trial['profile_identity'],
        'walk_projection_sign_stable': False,
        'walk_projection_min': anatomy['walk_projection_min'],
        'walk_projection_max': anatomy['walk_projection_max'],
        'solve_leftover_is_orthogonal_leftover': True,
        'jet_rows_do_not_own_leftover': True,
        'residual_split_maxima': {key: _round(value) if not isinstance(value, float)
                                  else _round(value)
                                  for key, value in pieces.items()},
        'locked_current_identified': False,
        'physical_NONEXISTENCE_certificate': False,
        'dropped_jet_step_norm': _round(float(np.linalg.norm(drop_delta))),
        'dropped_jet_clipped_step_norm': _round(clip_norm),
        'dropped_jet_unclipped_all_node_leftover': _round(unclipped_max),
        'dropped_jet_clipped_all_node_leftover': _round(clipped_max),
        'dropped_jet_unclipped_reaches_existence': unclipped_reaches,
        'dropped_jet_clipped_beats_measured': clipped_beats,
        'dropped_jet_clipped_is_order_below_measured': order_below,
        'evolve_authorized': evolve_authorized,
        'winning_linearization': None if not evolve_authorized else 'freed_jet_clipped',
        'scope': {
            'physical_local_gate': 'OPEN',
            'failed_optimization_is_not_NONEXISTENCE': True,
            'sign_stability_is_not_a_class_proof': True,
            'prediction_is_not_a_nonlinear_residual': True,
            'new_field_or_source_evolutions': 0,
            'physical_EXISTENCE_certificate': False,
            'metric_timestep': False,
        },
        'source_hashes': {p: digest(p) for p in OWNERS},
        'input_hashes': {
            str(A.OUTPUT.relative_to(ROOT)): digest(A.OUTPUT),
            str(I4.OUTPUT.relative_to(ROOT)): digest(I4.OUTPUT),
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
            raise FileExistsError('leftover identification exists; use --check')
        OUTPUT.write_text(json.dumps(result, sort_keys=True, indent=2, allow_nan=False)+'\n')
    elif result != json.loads(OUTPUT.read_text()):
        raise ValueError('leftover identification replay differs')
    print(json.dumps({
        'status': result['status'],
        'physical_NONEXISTENCE_certificate': result['physical_NONEXISTENCE_certificate'],
        'locked_current_identified': result['locked_current_identified'],
        'dropped_jet_unclipped_all_node_leftover': result['dropped_jet_unclipped_all_node_leftover'],
        'dropped_jet_clipped_all_node_leftover': result['dropped_jet_clipped_all_node_leftover'],
        'dropped_jet_step_norm': result['dropped_jet_step_norm'],
        'evolve_authorized': result['evolve_authorized'],
        'winning_linearization': result['winning_linearization'],
        'residual_split_maxima': result['residual_split_maxima'],
    }, indent=2))
