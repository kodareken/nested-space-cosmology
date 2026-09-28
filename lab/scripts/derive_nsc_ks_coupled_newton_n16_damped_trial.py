#!/usr/bin/env python3
"""Restartable full-source evaluation of the clipped n=16 rectangular candidate.

Reuses the n=16 evaluator, 32 retarded directions, 47-node union and 192-node
phase. Operators from any other g, including the zero-pad n=16 history, are
not reused. The linear prediction is not this residual.
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
import derive_nsc_ks_coupled_retained_control as R
import recursive_horizons.nsc_ks_history_evaluator as H
from recursive_horizons.nsc_ks_profile_identity import profile_identity
from recursive_horizons.nsc_local_incoming_family import LocalIncomingFamily

DIRECTORY = ROOT/'results/development/artifacts/nsc-ks-coupled-newton-n16-damped'
OUTPUT = ROOT/'results/development/nsc-ks-coupled-newton-n16-damped-trial.json'
EXPECTED_IDENTITY = 'ee033df12b7a5f91e1c97d2e0674ee6cba5a705e6a6a7aff5fed648e31cee645'
N16_ASSEMBLE = N16.assemble
OWNERS = (
    'scripts/derive_nsc_ks_coupled_newton_n16_damped_trial.py',
    'docs/nsc-ks-coupled-newton-n16-damped-trial.md',
    'scripts/derive_nsc_ks_coupled_newton_n16.py',
    'scripts/derive_nsc_ks_coupled_newton_n16_damped.py',
    'src/recursive_horizons/nsc_ks_retained_upstream_archive.py',
    'src/recursive_horizons/nsc_ks_history_evaluator.py',
    'src/recursive_horizons/nsc_ks_cutoff_bridge_evaluator.py',
    'src/recursive_horizons/nsc_dirac_source_phase.py',
    'src/recursive_horizons/nsc_ks_energy_propagator.py',
    'src/recursive_horizons/nsc_ks_difference_envelope.py',
    'src/recursive_horizons/nsc_ks_batched_constraints.py',
    'src/recursive_horizons/nsc_local_incoming_family.py',
)


def digest(path):
    return sha256((ROOT/path).read_bytes()).hexdigest()


def selected_rectangular(proposal):
    if not proposal.get('evolve_authorized'):
        raise ValueError('damped n=16 proposal did not authorize evolution')
    if proposal['selected_for_next_nonlinear_trial'] != 'rectangular':
        raise ValueError('damped nonlinear trial must remain the rectangular layout')
    matches = [row for row in proposal['proposals'] if row['layout'] == 'rectangular']
    if len(matches) != 1:
        raise ValueError('exactly one rectangular damped n=16 proposal required')
    row = matches[0]
    if row['candidate_profile_identity'] != EXPECTED_IDENTITY:
        raise ValueError('clipped rectangular candidate identity changed')
    if row['candidate_profile_identity'] in D.FORBIDDEN_IDENTITIES:
        raise ValueError('forbidden full-scale or stalled n=8 identity selected')
    if row.get('unclipped_full_scale_rejected') is not True:
        raise ValueError('unclipped full-scale n=16 candidate must stay rejected')
    if row.get('trial_accepted') is not False:
        raise ValueError('proposal must not pre-accept the trial')
    if not row.get('prediction_is_not_a_nonlinear_residual'):
        raise ValueError('linear prediction must remain marked as not a residual')
    if float(row['step_norm']) > 1.0 + 1e-12:
        raise ValueError('damped trial must keep the clipped unit step')
    return row


def context():
    base = R.context()
    proposal = json.loads(D.OUTPUT.read_text())
    for path, expected in {**proposal['source_hashes'], **proposal['input_hashes']}.items():
        if digest(path) != expected:
            raise ValueError('damped n=16 proposal input changed: '+path)
    selected = selected_rectangular(proposal)
    family = LocalIncomingFamily(np.asarray(selected['candidate']['coefficients'], float))
    identity = profile_identity(family, include_normal_window=True)
    if identity != selected['candidate_profile_identity'] or family.radius_lower_bound() <= 0:
        raise ValueError('same positive-radius clipped rectangular candidate required')
    solve = family.collocation_nodes(N16.SOLVE_COUNT)
    verify = family.collocation_nodes(N16.VERIFY_COUNT)
    target = H.fixed_target_union(solve, verify)
    if len(target) != 47:
        raise ValueError('damped n=16 solve/verification union must have 47 nodes')
    if target.tolist() != proposal['z']:
        raise ValueError('damped trial must keep the measured n=16 node union')
    hashes = {
        **base['archive'].input_hashes,
        **{path: digest(path) for path in OWNERS},
        str(D.OUTPUT.relative_to(ROOT)): digest(D.OUTPUT),
    }
    return {
        **base,
        'family': family,
        'identity': identity,
        'solve': solve,
        'verify': verify,
        'target': target,
        'hashes': hashes,
        'proposal': proposal,
        'selected': selected,
    }


def assemble(ctx, replay=False):
    result = N16_ASSEMBLE(ctx, replay=replay)
    complete = result['all_retained_nonzero_angular_families_included']
    predicted = np.asarray(ctx['selected']['predicted_all_node_residual_maxima'], float)
    previous = np.asarray(ctx['proposal']['constraint_maxima_with_192_node_phase'], float)
    measured = np.asarray(result['constraint_maxima'], float)
    improved = bool(complete and np.all(measured < previous) and ctx['family'].radius_lower_bound() > 0)
    result['schema'] = 'NSC-KS-COUPLED-NEWTON-N16-DAMPED-TRIAL-v1'
    result['status'] = (
        'OPEN: complete clipped n=16 rectangular residual; uncertified physical gate'
        if complete else
        'OPEN: partial clipped n=16 rectangular residual')
    result['parent_profile_identity'] = ctx['proposal']['current_profile_identity']
    result['parent_n8_profile_identity'] = N16.PARENT_IDENTITY
    result['layout'] = 'rectangular-n16-damped'
    result['linear_predicted_all_node_residual_maxima'] = predicted.tolist()
    result['prediction_minus_measured'] = (
        None if not complete else (predicted - measured).tolist())
    result['residual_improved_on_all_components'] = improved if complete else None
    result['n16_constraint_maxima'] = previous.tolist()
    result['reproducer'] = 'python3 scripts/derive_nsc_ks_coupled_newton_n16_damped_trial.py --check'
    return result


def family_paths(key):
    return N16.family_paths(key)


def bind():
    N16.DIRECTORY = DIRECTORY
    N16.OUTPUT = OUTPUT
    N16.OWNERS = OWNERS
    N16.context = context
    N16.assemble = assemble


if __name__ == '__main__':
    bind()
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--run', action='store_true')
    mode.add_argument('--check', action='store_true')
    parser.add_argument('--cpu-budget', type=float, default=R.CPU_CAP)
    parser.add_argument('--max-new', type=int, default=N16.MAX_NEW)
    args = parser.parse_args()
    if args.run:
        result = N16.run(args.cpu_budget, args.max_new)
        runtime = result.get('runtime')
    else:
        ctx = context()
        result = assemble(ctx, replay=True)
        path = OUTPUT if OUTPUT.exists() else DIRECTORY/'progress.json'
        previous = json.loads(path.read_text())
        previous.pop('runtime', None)
        if R.jsonable(result) != previous:
            raise ValueError('clipped n=16 trial retained-source replay differs')
        runtime = {'CPU_seconds': None, 'new_operator_solves': 0,
                   'reused_operator_solves': 0, 'replay_dirac_evolution': False}
    print(json.dumps({key: result[key] for key in (
        'status', 'completed_positive_families', 'required_positive_families',
        'positive_energy_rows', 'constraint_maxima', 'raw_constraint_maxima',
        'matter_change_maxima', 'linear_predicted_all_node_residual_maxima',
        'residual_improved_on_all_components', 'trial_accepted',
        'CPU_seconds', 'new_operator_solves', 'reused_operator_solves') if key in result}, indent=2))
    print(json.dumps(runtime, indent=2))
