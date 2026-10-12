#!/usr/bin/env python3
"""Accept the clipped n=16 trial as a solver step and propose the next clip.

The measured residual matched the linear prediction to about 1e-9. Acceptance
is a solver-step decision, not a physical local-gate certificate. The next
step reuses the saved trial Jacobian, clips to max_step_norm=1, and walks
the existing line-search scales. A linear prediction is not a residual.
"""
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
import derive_nsc_ks_coupled_newton_n16_damped_trial as T
from recursive_horizons.nsc_local_history_newton import (
    DeclaredJetBoundary, LocalHistoryNewtonSettings)
from recursive_horizons.nsc_ks_profile_identity import profile_identity

OUTPUT = ROOT/'results/development/nsc-ks-coupled-newton-n16-damped-next.json'
OWNERS = (
    'scripts/derive_nsc_ks_coupled_newton_n16_damped_next.py',
    'docs/nsc-ks-coupled-newton-n16-damped-next.md',
    'src/recursive_horizons/nsc_local_history_newton.py',
    'src/recursive_horizons/nsc_evolved_incoming_constraints.py',
    'src/recursive_horizons/nsc_local_incoming_family.py',
)


def digest(path):
    path = Path(path)
    return sha256((path if path.is_absolute() else ROOT/path).read_bytes()).hexdigest()


def inputs():
    trial = json.loads(T.OUTPUT.read_text())
    for path, expected in trial['source_hashes'].items():
        if digest(path) != expected:
            raise ValueError('clipped n=16 trial input changed: '+path)
    if trial['profile_identity'] != T.EXPECTED_IDENTITY:
        raise ValueError('accepted trial identity changed')
    if trial['completed_positive_families'] != trial['required_positive_families']:
        raise ValueError('incomplete clipped trial cannot be accepted')
    if trial['new_operator_solves'] != 60 or trial['reused_operator_solves'] != 0:
        raise ValueError('accepted clipped trial must evolve all sixty families')
    if trial['source_cutoff_phase_gauss_nodes'] != 192:
        raise ValueError('accepted clipped trial must use the 192-node phase')
    if trial['retarded_directions'] != 32:
        raise ValueError('accepted clipped trial must keep 32 retarded directions')
    if trial['scope'].get('group14_operators_reused_from_coupled_pilot'):
        raise ValueError('accepted trial reused a different-g operator')
    if not trial['residual_improved_on_all_components']:
        raise ValueError('clipped trial residual did not improve on both components')
    if trial['radius_lower_bound'] <= 0:
        raise ValueError('accepted clipped trial must keep a positive radius bound')
    gradient = np.asarray(trial['action_gradient'], float)
    jacobian = np.asarray(trial['history_jacobian'], float)
    measured = np.asarray(trial['constraint_maxima'], float)
    previous = np.asarray(trial['previous_constraint_maxima'], float)
    if not np.isfinite(gradient).all() or not np.isfinite(jacobian).all():
        raise ValueError('accepted clipped trial residual/Jacobian must be finite')
    if not np.all(measured < previous):
        raise ValueError('clipped trial residual did not improve on both components')
    ctx = T.context()
    if profile_identity(ctx['family']) != trial['profile_identity']:
        raise ValueError('next damped step must start from the measured clipped trial')
    if ctx['target'].tolist() != trial['z']:
        raise ValueError('next damped step must keep the measured n=16 node union')
    source_identity = json.loads(D.OUTPUT.read_text())['source_identity']
    if source_identity.get('source_preparation_changed'):
        raise ValueError('next damped step must keep the same upstream source')
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
    evolve_authorized = selected_rect is not None
    if evolve_authorized:
        status = (
            'OPEN: clipped n=16 trial accepted as a solver step; '
            'next clipped prediction improves and needs fresh evolution')
        stall = None
        selected_name = 'rectangular'
    else:
        status = (
            'OPEN: clipped n=16 trial accepted as a solver step; '
            'next damped linear prediction does not improve')
        stall = 'n16_damped_next_linear_prediction_does_not_improve'
        selected_name = None
    prediction = np.asarray(trial['linear_predicted_all_node_residual_maxima'], float)
    return {
        'schema': 'NSC-KS-COUPLED-NEWTON-N16-DAMPED-NEXT-v1',
        'accountable_author': 'Douglas Ek',
        'status': status,
        'accepted_trial': {
            'path': str(T.OUTPUT.relative_to(ROOT)),
            'profile_identity': trial['profile_identity'],
            'solver_step_accepted': True,
            'physical_trial_accepted_in_source_register': False,
            'measured_constraint_maxima': measured.tolist(),
            'previous_constraint_maxima': trial['previous_constraint_maxima'],
            'linear_predicted_all_node_residual_maxima': prediction.tolist(),
            'prediction_minus_measured': (prediction - measured).tolist(),
            'radius_lower_bound': trial['radius_lower_bound'],
            'all_nodes_finite': True,
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
        'phase_gauss_nodes': 192,
        'nodal_phase_value_error_upper_N_beta': None,
        'boundary_control': DeclaredJetBoundary(
            family.center, 0., ctx['seed_slope'], 0.).description(),
        'forbidden_identities': list(D.FORBIDDEN_IDENTITIES),
        'proposals': proposals,
        'selected_for_next_nonlinear_trial': selected_name,
        'evolve_authorized': evolve_authorized,
        'named_search_stall': stall,
        'scope': {
            'physical_local_gate': 'OPEN',
            'missing_total_error_bound': None,
            'between_node_bound': None,
            'source_cutoff_phase_quadrature_error_bound': None,
            'old_phase_bounds_transferred': False,
            'new_field_or_source_evolutions': 0,
            'metric_timestep': False,
            'used_existing_clip_and_line_search': True,
            'full_scale_n16_not_evolved': True,
            'n8_next_not_evolved': True,
            'solver_step_is_not_physical_EXISTENCE': True,
            'prediction_is_not_a_nonlinear_residual': True,
            'physical_NONEXISTENCE_claimed': False,
        },
        'source_hashes': {p: digest(p) for p in OWNERS},
        'input_hashes': {str(p.relative_to(ROOT)): digest(p) for p in (T.OUTPUT, D.OUTPUT)},
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
            raise FileExistsError('next clipped n=16 proposal exists; use --check')
        OUTPUT.write_text(json.dumps(result, sort_keys=True, indent=2, allow_nan=False)+'\n')
    elif result != json.loads(OUTPUT.read_text()):
        raise ValueError('next clipped n=16 proposal replay differs')
    print(json.dumps({
        'status': result['status'],
        'accepted_maxima': result['accepted_trial']['measured_constraint_maxima'],
        'evolve_authorized': result['evolve_authorized'],
        'named_search_stall': result['named_search_stall'],
        'selected_for_next_nonlinear_trial': result['selected_for_next_nonlinear_trial'],
        'proposals': [{k: row[k] for k in row if k not in ('direction', 'candidate', 'scales')}
                      for row in result['proposals']],
    }, indent=2))
