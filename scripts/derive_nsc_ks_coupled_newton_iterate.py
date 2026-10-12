#!/usr/bin/env python3
"""Restartable full-source evaluation of the next rectangular Newton candidate.

This reuses the first-trial evaluator with a different history identity,
artifact directory and owner hashes. Operators from any previous g are not
reused. The linear prediction is not this residual.
"""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
sys.path.insert(0, str(ROOT/'scripts'))
import derive_nsc_ks_coupled_newton_next as N
import derive_nsc_ks_coupled_newton_trial as T
import derive_nsc_ks_coupled_retained_control as R

T.PROPOSAL = N.OUTPUT
T.DIRECTORY = ROOT/'results/development/artifacts/nsc-ks-coupled-newton-iterate'
T.OUTPUT = ROOT/'results/development/nsc-ks-coupled-newton-iterate.json'
T.OWNERS = (
    'scripts/derive_nsc_ks_coupled_newton_iterate.py',
    'docs/nsc-ks-coupled-newton-iterate.md',
    'scripts/derive_nsc_ks_coupled_newton_next.py',
    'scripts/derive_nsc_ks_coupled_newton_trial.py',
    'src/recursive_horizons/nsc_ks_retained_upstream_archive.py',
    'src/recursive_horizons/nsc_ks_history_evaluator.py',
    'src/recursive_horizons/nsc_ks_cutoff_bridge_evaluator.py',
    'src/recursive_horizons/nsc_dirac_source_phase.py',
    'src/recursive_horizons/nsc_ks_energy_propagator.py',
    'src/recursive_horizons/nsc_ks_difference_envelope.py',
    'src/recursive_horizons/nsc_ks_batched_constraints.py',
    'src/recursive_horizons/nsc_local_incoming_family.py',
)


def bind_next_candidate():
    proposal = json.loads(N.OUTPUT.read_text())
    if proposal['selected_for_next_nonlinear_trial'] != 'rectangular':
        raise ValueError('next nonlinear trial must remain the rectangular layout')
    matches = [row for row in proposal['proposals'] if row['layout'] == 'rectangular']
    if len(matches) != 1:
        raise ValueError('exactly one next rectangular Newton proposal required')
    T.EXPECTED_IDENTITY = matches[0]['candidate_profile_identity']
    return T.EXPECTED_IDENTITY


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--run', action='store_true')
    mode.add_argument('--check', action='store_true')
    parser.add_argument('--cpu-budget', type=float, default=R.CPU_CAP)
    parser.add_argument('--max-new', type=int, default=60)
    args = parser.parse_args()
    bind_next_candidate()
    if args.run:
        result = T.run(args.cpu_budget, args.max_new)
        runtime = result.get('runtime')
    else:
        ctx = T.context()
        result = T.assemble(ctx, replay=True)
        path = T.OUTPUT if T.OUTPUT.exists() else T.DIRECTORY/'progress.json'
        previous = json.loads(path.read_text())
        previous.pop('runtime', None)
        if R.jsonable(result) != previous:
            raise ValueError('Newton-iterate retained-source replay differs')
        runtime = {'CPU_seconds': None, 'new_operator_solves': 0,
                   'reused_operator_solves': 0, 'replay_dirac_evolution': False}
    print(json.dumps({key: result[key] for key in (
        'status', 'completed_positive_families', 'required_positive_families',
        'positive_energy_rows', 'constraint_maxima', 'raw_constraint_maxima',
        'matter_change_maxima', 'previous_constraint_maxima',
        'linear_predicted_all_node_residual_maxima',
        'residual_improved_on_all_components', 'trial_accepted',
        'CPU_seconds', 'new_operator_solves', 'reused_operator_solves') if key in result}, indent=2))
    print(json.dumps(runtime, indent=2))
