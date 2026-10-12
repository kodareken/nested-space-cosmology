#!/usr/bin/env python3
"""Accept the fifth clipped n=16 trial and propose the next unit step."""
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
import derive_nsc_ks_coupled_newton_n16 as N16
import derive_nsc_ks_coupled_newton_n16_damped as D
import derive_nsc_ks_coupled_newton_n16_damped_iterate4 as I4
from recursive_horizons.nsc_local_history_newton import (
    DeclaredJetBoundary, LocalHistoryNewtonSettings)
from recursive_horizons.nsc_ks_profile_identity import profile_identity

OUTPUT = ROOT/'results/development/nsc-ks-coupled-newton-n16-damped-next5.json'
OWNERS = (
    'scripts/derive_nsc_ks_coupled_newton_n16_damped_next5.py',
    'docs/nsc-ks-coupled-newton-n16-damped-next5.md',
    'src/recursive_horizons/nsc_local_history_newton.py',
    'src/recursive_horizons/nsc_evolved_incoming_constraints.py',
    'src/recursive_horizons/nsc_local_incoming_family.py',
)
EXISTENCE_TOLERANCE = 3e-11


def digest(path):
    path = Path(path)
    return sha256((path if path.is_absolute() else ROOT/path).read_bytes()).hexdigest()


def unclipped_leftover(gradient, jacobian, direction, unclipped_norm):
    step = np.asarray(direction, float).reshape(-1) * float(unclipped_norm)
    return D.predicted_maxima(gradient, jacobian, step)


def inputs():
    trial = json.loads(I4.OUTPUT.read_text())
    for path, expected in trial['source_hashes'].items():
        if digest(path) != expected:
            raise ValueError('fifth clipped trial input changed: '+path)
    if trial['profile_identity'] != I4.EXPECTED_IDENTITY:
        raise ValueError('fifth clipped trial identity changed')
    if trial['completed_positive_families'] != trial['required_positive_families']:
        raise ValueError('incomplete fifth clipped trial cannot be accepted')
    if trial['new_operator_solves'] != 60 or trial['reused_operator_solves'] != 0:
        raise ValueError('fifth clipped trial must evolve all sixty families')
    if not trial['residual_improved_on_all_components'] or trial['radius_lower_bound'] <= 0:
        raise ValueError('fifth clipped trial must improve both residuals with positive radius')
    gradient = np.asarray(trial['action_gradient'], float)
    jacobian = np.asarray(trial['history_jacobian'], float)
    measured = np.asarray(trial['constraint_maxima'], float)
    if not np.all(measured < np.asarray(trial['previous_constraint_maxima'], float)):
        raise ValueError('fifth clipped trial residual did not improve on both components')
    ctx = I4.context()
    source_identity = json.loads(D.OUTPUT.read_text())['source_identity']
    return ctx, trial, source_identity, gradient, jacobian, measured


def compute():
    ctx, trial, source_identity, gradient, jacobian, measured = inputs()
    family = ctx['family']
    settings = LocalHistoryNewtonSettings()
    D.FORBIDDEN_IDENTITIES = D.FORBIDDEN_IDENTITIES + (trial['profile_identity'],)
    proposals = []
    selected_rect = None
    for layout in ('rectangular', 'tau'):
        row, selected = D.propose_layout(
            ctx, layout, gradient, jacobian, source_identity, measured, settings)
        proposals.append(row)
        if layout == 'rectangular':
            selected_rect = selected
    rect = proposals[0]
    leftover = unclipped_leftover(
        gradient, jacobian, rect['direction'], rect['unclipped_step_norm'])
    leftover_decides = bool(np.all(leftover <= EXISTENCE_TOLERANCE))
    clip_can_decide = leftover_decides
    evolve_authorized = selected_rect is not None and clip_can_decide
    if evolve_authorized:
        status = (
            'OPEN: fifth clipped n=16 trial accepted as a solver step; '
            'unclipped leftover can reach the gate and needs fresh evolution')
        stall = None
    elif selected_rect is not None:
        status = (
            'OPEN: fifth clipped n=16 trial accepted as a solver step; '
            'clipped prediction improves but unclipped leftover cannot decide the gate')
        stall = 'n16_damped_frozen_jacobian_leftover_cannot_reach_existence'
    else:
        status = (
            'OPEN: fifth clipped n=16 trial accepted as a solver step; '
            'next damped linear prediction does not improve')
        stall = 'n16_damped_next5_linear_prediction_does_not_improve'
    prediction = np.asarray(trial['linear_predicted_all_node_residual_maxima'], float)
    return {
        'schema': 'NSC-KS-COUPLED-NEWTON-N16-DAMPED-NEXT5-v1',
        'accountable_author': 'Douglas Ek',
        'status': status,
        'accepted_trial': {
            'path': str(I4.OUTPUT.relative_to(ROOT)),
            'profile_identity': trial['profile_identity'],
            'solver_step_accepted': True,
            'measured_constraint_maxima': measured.tolist(),
            'previous_constraint_maxima': trial['previous_constraint_maxima'],
            'linear_predicted_all_node_residual_maxima': prediction.tolist(),
            'prediction_minus_measured': (prediction - measured).tolist(),
            'radius_lower_bound': trial['radius_lower_bound'],
            'all_sixty_families_evolved_for_this_g': True,
            'physical_local_gate': 'OPEN',
            'physical_EXISTENCE_certificate': False,
        },
        'current_history': family.description(),
        'current_profile_identity': profile_identity(family),
        'parent_n8_profile_identity': N16.PARENT_IDENTITY,
        'source_identity': source_identity,
        'z': ctx['target'].tolist(),
        'action_gradient_with_192_node_phase': gradient.tolist(),
        'full_retarded_history_jacobian': jacobian.tolist(),
        'constraint_maxima_with_192_node_phase': measured.tolist(),
        'constraint_maxima_with_verified_phase_values': measured.tolist(),
        'proposals': proposals,
        'unclipped_linear_predicted_all_node_residual_maxima': leftover.tolist(),
        'unclipped_leftover_reaches_existence_tolerance': leftover_decides,
        'existence_tolerance': EXISTENCE_TOLERANCE,
        'selected_for_next_nonlinear_trial': 'rectangular' if evolve_authorized else None,
        'evolve_authorized': evolve_authorized,
        'named_search_stall': stall,
        'boundary_control': DeclaredJetBoundary(
            family.center, 0., ctx['seed_slope'], 0.).description(),
        'scope': {
            'physical_local_gate': 'OPEN',
            'missing_total_error_bound': None,
            'between_node_bound': None,
            'new_field_or_source_evolutions': 0,
            'metric_timestep': False,
            'solver_step_is_not_physical_EXISTENCE': True,
            'prediction_is_not_a_nonlinear_residual': True,
            'physical_NONEXISTENCE_claimed': False,
            'frozen_jacobian_walk_stopped': not leftover_decides,
        },
        'source_hashes': {p: digest(p) for p in OWNERS},
        'input_hashes': {str(p.relative_to(ROOT)): digest(p) for p in (I4.OUTPUT, D.OUTPUT)},
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
            raise FileExistsError('fifth next clipped proposal exists; use --check')
        OUTPUT.write_text(json.dumps(result, sort_keys=True, indent=2, allow_nan=False)+'\n')
    elif result != json.loads(OUTPUT.read_text()):
        raise ValueError('fifth next clipped proposal replay differs')
    print(json.dumps({
        'status': result['status'],
        'accepted_maxima': result['accepted_trial']['measured_constraint_maxima'],
        'evolve_authorized': result['evolve_authorized'],
        'named_search_stall': result['named_search_stall'],
        'unclipped_leftover': result['unclipped_linear_predicted_all_node_residual_maxima'],
        'unclipped_leftover_reaches_existence_tolerance': result[
            'unclipped_leftover_reaches_existence_tolerance'],
        'proposals': [{k: row[k] for k in row if k not in ('direction', 'candidate', 'scales')}
                      for row in result['proposals']],
    }, indent=2))
