#!/usr/bin/env python3
"""Restartable full-source evaluation of the fourth clipped n=16 candidate."""
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
import derive_nsc_ks_coupled_newton_n16_damped_next3 as N3
import derive_nsc_ks_coupled_retained_control as R
import recursive_horizons.nsc_ks_history_evaluator as H
from recursive_horizons.nsc_ks_profile_identity import profile_identity
from recursive_horizons.nsc_local_incoming_family import LocalIncomingFamily

DIRECTORY = ROOT/'results/development/artifacts/nsc-ks-coupled-newton-n16-damped-iterate3'
OUTPUT = ROOT/'results/development/nsc-ks-coupled-newton-n16-damped-iterate3.json'
EXPECTED_IDENTITY = '045a3ab26eb7ff85ca40772ca1aaeb7e0ae6a02fa4b51f25386e036521c177d0'
OWNERS = (
    'scripts/derive_nsc_ks_coupled_newton_n16_damped_iterate3.py',
    'docs/nsc-ks-coupled-newton-n16-damped-iterate3.md',
    'scripts/derive_nsc_ks_coupled_newton_n16.py',
    'scripts/derive_nsc_ks_coupled_newton_n16_damped_next3.py',
    'src/recursive_horizons/nsc_ks_retained_upstream_archive.py',
    'src/recursive_horizons/nsc_ks_history_evaluator.py',
    'src/recursive_horizons/nsc_ks_cutoff_bridge_evaluator.py',
    'src/recursive_horizons/nsc_dirac_source_phase.py',
    'src/recursive_horizons/nsc_ks_energy_propagator.py',
    'src/recursive_horizons/nsc_ks_difference_envelope.py',
    'src/recursive_horizons/nsc_ks_batched_constraints.py',
    'src/recursive_horizons/nsc_local_incoming_family.py',
)
N16_ASSEMBLE = N16.assemble


def digest(path):
    return sha256((ROOT/path).read_bytes()).hexdigest()


def selected_rectangular(proposal):
    if not proposal.get('evolve_authorized'):
        raise ValueError('fourth clipped proposal did not authorize evolution')
    matches = [row for row in proposal['proposals'] if row['layout'] == 'rectangular']
    row = matches[0]
    if row['candidate_profile_identity'] != EXPECTED_IDENTITY:
        raise ValueError('fourth clipped candidate identity changed')
    if row['candidate_profile_identity'] in D.FORBIDDEN_IDENTITIES:
        raise ValueError('forbidden identity selected')
    if float(row['step_norm']) > 1.0 + 1e-12:
        raise ValueError('fourth clipped trial must keep the unit step')
    return row


def context():
    base = R.context()
    proposal = json.loads(N3.OUTPUT.read_text())
    for path, expected in {**proposal['source_hashes'], **proposal['input_hashes']}.items():
        if digest(path) != expected:
            raise ValueError('fourth clipped proposal input changed: '+path)
    selected = selected_rectangular(proposal)
    proposal = dict(proposal)
    proposal['constraint_maxima_with_verified_phase_values'] = proposal[
        'constraint_maxima_with_192_node_phase']
    family = LocalIncomingFamily(np.asarray(selected['candidate']['coefficients'], float))
    identity = profile_identity(family, include_normal_window=True)
    if identity != selected['candidate_profile_identity'] or family.radius_lower_bound() <= 0:
        raise ValueError('same positive-radius fourth clipped candidate required')
    solve = family.collocation_nodes(N16.SOLVE_COUNT)
    verify = family.collocation_nodes(N16.VERIFY_COUNT)
    target = H.fixed_target_union(solve, verify)
    hashes = {
        **base['archive'].input_hashes,
        **{path: digest(path) for path in OWNERS},
        str(N3.OUTPUT.relative_to(ROOT)): digest(N3.OUTPUT),
    }
    return {
        **base, 'family': family, 'identity': identity, 'solve': solve,
        'verify': verify, 'target': target, 'hashes': hashes,
        'proposal': proposal, 'selected': selected,
    }


def assemble(ctx, replay=False):
    result = N16_ASSEMBLE(ctx, replay=replay)
    complete = result['all_retained_nonzero_angular_families_included']
    predicted = np.asarray(ctx['selected']['predicted_all_node_residual_maxima'], float)
    previous = np.asarray(ctx['proposal']['constraint_maxima_with_192_node_phase'], float)
    measured = np.asarray(result['constraint_maxima'], float)
    improved = bool(complete and np.all(measured < previous) and ctx['family'].radius_lower_bound() > 0)
    result['schema'] = 'NSC-KS-COUPLED-NEWTON-N16-DAMPED-ITERATE3-v1'
    result['status'] = (
        'OPEN: complete fourth clipped n=16 residual; uncertified physical gate'
        if complete else 'OPEN: partial fourth clipped n=16 residual')
    result['parent_profile_identity'] = ctx['proposal']['current_profile_identity']
    result['layout'] = 'rectangular-n16-damped-iterate3'
    result['linear_predicted_all_node_residual_maxima'] = predicted.tolist()
    result['prediction_minus_measured'] = None if not complete else (predicted - measured).tolist()
    result['residual_improved_on_all_components'] = improved if complete else None
    result['reproducer'] = 'python3 scripts/derive_nsc_ks_coupled_newton_n16_damped_iterate3.py --check'
    return result


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
            raise ValueError('fourth clipped trial retained-source replay differs')
        runtime = {'CPU_seconds': None, 'new_operator_solves': 0,
                   'reused_operator_solves': 0, 'replay_dirac_evolution': False}
    print(json.dumps({key: result[key] for key in (
        'status', 'completed_positive_families', 'required_positive_families',
        'constraint_maxima', 'linear_predicted_all_node_residual_maxima',
        'residual_improved_on_all_components', 'new_operator_solves') if key in result}, indent=2))
    print(json.dumps(runtime, indent=2))
