#!/usr/bin/env python3
"""One rank-64 step from the primal-block remeasure of iterate6.

The parent residual is the zero-padded iterate6 history evolved with
primal-block DOP853 control. The step is the equilibrated rank-64 truncation
of that measured Jacobian, already inside the unit ball. It is evolved with
the same controller. A linear prediction is not a residual.
"""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
sys.path.insert(0, str(ROOT/'scripts'))
import derive_nsc_ks_coupled_newton_n16 as N16
import derive_nsc_ks_coupled_newton_n16_damped as D
import derive_nsc_ks_coupled_retained_control as R
import recursive_horizons.nsc_ks_energy_propagator as P
import recursive_horizons.nsc_ks_history_evaluator as H
from recursive_horizons.nsc_ks_primal_step_control import evolve_primal_controlled
from recursive_horizons.nsc_ks_profile_identity import profile_identity

import derive_nsc_ks_coupled_newton_n64_truncated_iterate as I64

PARENT = ROOT/'results/development/nsc-ks-coupled-newton-n64-primal-norm-zeropad.json'
ITERATE6 = ROOT/'results/development/nsc-ks-coupled-newton-n32-truncated-iterate6.json'
DIRECTORY = ROOT/'results/development/artifacts/nsc-ks-coupled-newton-n64-primal-norm-rank64'
OUTPUT = ROOT/'results/development/nsc-ks-coupled-newton-n64-primal-norm-rank64.json'
PROPOSAL = ROOT/'results/development/nsc-ks-coupled-newton-n64-primal-norm-rank64-proposal.json'
OWNERS = (
    'scripts/derive_nsc_ks_n64_primal_norm_rank64.py',
    'docs/nsc-ks-n64-primal-norm-rank64.md',
    'src/recursive_horizons/nsc_ks_primal_step_control.py',
)
KEPT = 64
EXISTENCE = 3e-11


def digest(path):
    path = Path(path)
    return sha256((path if path.is_absolute() else ROOT/path).read_bytes()).hexdigest()


def proposal():
    parent = json.loads(PARENT.read_text())
    iterate6 = json.loads(ITERATE6.read_text())
    if parent['schema'] != 'NSC-KS-N64-PRIMAL-NORM-ZEROPAD-v1':
        raise ValueError('rank-64 parent must be the primal-block zero pad')
    if parent['coefficient_step_evolved'] or parent['newton_step_accepted']:
        raise ValueError('the zero pad is not a previous coefficient step')
    gradient = np.asarray(parent['action_gradient'], float)
    jacobian = np.asarray(parent['history_jacobian'], float)
    if jacobian.shape[0] != I64.DIRECTIONS:
        raise ValueError('parent Jacobian must carry 128 family-ordered directions')
    coefficients = np.asarray(parent['history']['coefficients'], float)
    delta, residual, matrix, singular = I64.rank_truncated_delta(gradient, jacobian, KEPT)
    _clipped, step, raw_norm, step_norm = I64.clip_and_scale(delta, 1.0)
    if not np.array_equal(step, delta):
        raise ValueError('rank-64 primal-norm step must already lie inside the unit ball')
    predicted = I64.predicted_maxima(residual, matrix, step)
    parent_max = np.asarray(parent['constraint_maxima'], float)
    iterate6_max = np.asarray(iterate6['constraint_maxima'], float)
    if not np.all(predicted < parent_max) or not np.all(predicted < iterate6_max):
        raise ValueError('rank-64 prediction must beat both the primal-norm parent and iterate6')
    candidate = I64.LocalIncomingFamily64(coefficients + step.reshape(2, KEPT))
    identity = profile_identity(candidate, include_normal_window=True)
    if identity in D.FORBIDDEN_IDENTITIES or identity == parent['profile_identity']:
        raise ValueError('candidate identity is forbidden or unchanged')
    if candidate.radius_lower_bound() <= 0:
        raise ValueError('candidate radius bound must stay positive')
    solve = candidate.collocation_nodes(N16.SOLVE_COUNT)
    verify = candidate.collocation_nodes(N16.VERIFY_COUNT)
    target = H.fixed_target_union(solve, verify)
    if target.tolist() != parent['z']:
        raise ValueError('rank-64 step must stay on the iterate6 nodes')
    return {
        'schema': 'NSC-KS-N64-PRIMAL-NORM-RANK64-PROPOSAL-v1',
        'accountable_author': 'Douglas Ek',
        'status': 'OPEN: one rank-64 primal-norm step authorized for measurement',
        'kept_modes': KEPT,
        'step_scale': 1.0,
        'raw_step_norm': raw_norm,
        'step_norm': step_norm,
        'step': step.tolist(),
        'parent_profile_identity': parent['profile_identity'],
        'parent_constraint_maxima': parent_max.tolist(),
        'iterate6_constraint_maxima': iterate6_max.tolist(),
        'predicted_all_node_residual_maxima': predicted.tolist(),
        'prediction_beats_iterate6_on_both_components': True,
        'prediction_is_not_a_nonlinear_residual': True,
        'candidate_profile_identity': identity,
        'candidate_coefficients': candidate.description()['coefficients'],
        'radius_lower_bound': candidate.radius_lower_bound(),
        'z': parent['z'],
        'evolve_authorized': True,
        'physical_EXISTENCE_certificate': False,
        'physical_NONEXISTENCE_certificate': False,
        'source_hashes': {path: digest(path) for path in OWNERS},
        'input_hashes': {
            str(PARENT.relative_to(ROOT)): digest(PARENT),
            str(ITERATE6.relative_to(ROOT)): digest(ITERATE6),
        },
    }


def context(record):
    if not record.get('evolve_authorized'):
        raise ValueError('rank-64 evolution requires an authorized prediction')
    family = I64.LocalIncomingFamily64(np.asarray(record['candidate_coefficients'], float))
    identity = profile_identity(family, include_normal_window=True)
    if identity != record['candidate_profile_identity'] or identity in D.FORBIDDEN_IDENTITIES:
        raise ValueError('candidate identity changed or is forbidden')
    base = R.context()
    solve = family.collocation_nodes(N16.SOLVE_COUNT)
    verify = family.collocation_nodes(N16.VERIFY_COUNT)
    target = H.fixed_target_union(solve, verify)
    if target.tolist() != record['z']:
        raise ValueError('rank-64 nodes changed')
    return {
        **base,
        'family': family,
        'identity': identity,
        'solve': solve,
        'verify': verify,
        'target': target,
        'hashes': {
            **base['hashes'],
            **{path: digest(path) for path in OWNERS},
            **record['input_hashes'],
        },
        'proposal': {
            'constraint_maxima_with_verified_phase_values': record['parent_constraint_maxima'],
        },
        'selected': {'step_scale': record['step_scale'], 'step_norm': record['step_norm']},
        'rank64_proposal': record,
    }


def annotate(result, ctx):
    record = ctx['rank64_proposal']
    measured = np.asarray(result['constraint_maxima'], float)
    predicted = np.asarray(record['predicted_all_node_residual_maxima'], float)
    parent = np.asarray(record['parent_constraint_maxima'], float)
    iterate6 = np.asarray(record['iterate6_constraint_maxima'], float)
    complete = result['all_retained_nonzero_angular_families_included']
    improved_parent = bool(complete and np.all(measured < parent) and ctx['family'].radius_lower_bound() > 0)
    beats_iterate6 = bool(complete and np.all(measured < iterate6))
    result['schema'] = 'NSC-KS-N64-PRIMAL-NORM-RANK64-v1'
    result['status'] = (
        'OPEN: measured rank-64 primal-norm step; uncertified physical gate'
        if complete else 'OPEN: partial rank-64 primal-norm step')
    result['solver'] = 'primal-block DOP853 RMS on (A, D); tangent directions excluded'
    result['kept_modes'] = KEPT
    result['step_scale'] = record['step_scale']
    result['step_norm'] = record['step_norm']
    result['linear_predicted_all_node_residual_maxima'] = predicted.tolist()
    result['prediction_minus_measured'] = None if not complete else (predicted - measured).tolist()
    result['iterate6_constraint_maxima'] = iterate6.tolist()
    result['parent_constraint_maxima'] = parent.tolist()
    result['parent_profile_identity'] = record['parent_profile_identity']
    result['residual_improved_on_all_components'] = improved_parent if complete else None
    result['beats_iterate6_on_both_components'] = beats_iterate6 if complete else None
    result['prediction_is_not_a_nonlinear_residual'] = True
    result['coefficient_step_evolved'] = True
    result['newton_step_accepted'] = bool(beats_iterate6 and improved_parent)
    result['physical_EXISTENCE_certificate'] = False
    result['physical_NONEXISTENCE_certificate'] = False
    result['nodal_residual_reaches_existence_tolerance'] = bool(
        complete and np.all(measured <= EXISTENCE))
    result['enclosed_error'] = None
    result['reproducer'] = 'python3 scripts/derive_nsc_ks_n64_primal_norm_rank64.py --check'
    return result


def run(cpu_budget, max_new):
    P.evolve_ks_difference_envelope = evolve_primal_controlled
    N16.RETARDED = I64.DIRECTIONS
    N16.DIRECTORY = DIRECTORY
    N16.OUTPUT = OUTPUT
    R.FAMILY_CAP = 240.0
    R.CPU_CAP = 14400.0
    if PROPOSAL.exists():
        record = json.loads(PROPOSAL.read_text())
    else:
        record = proposal()
        PROPOSAL.write_text(json.dumps(record, indent=2, sort_keys=True)+'\n')
    ctx = context(record)
    DIRECTORY.mkdir(parents=True, exist_ok=True)
    start = time.process_time()
    new = 0
    for key in sorted(ctx['archived']):
        path, _payload = N16.family_paths(key)
        if path.exists():
            continue
        if new >= max_new:
            continue
        remaining = cpu_budget - (time.process_time() - start)
        if remaining <= 1.0:
            break
        saved = N16.new_family(ctx, key, cpu_limit=remaining - 1.0)
        new += saved['new_operator_solves']
        print(json.dumps({
            'positive_family': saved['positive_family'],
            'CPU_seconds': saved['CPU_seconds'],
            'accepted_steps': saved['operator']['diagnostics']['accepted_steps'],
            'matter_change_maxima': saved['matter_change_maxima'],
        }), flush=True)
    result = annotate(N16.assemble(ctx), ctx)
    target = OUTPUT if result['all_retained_nonzero_angular_families_included'] else DIRECTORY/'progress.json'
    if target == OUTPUT:
        progress = DIRECTORY/'progress.json'
        if progress.exists():
            progress.unlink()
    target.write_text(json.dumps(R.jsonable(result), indent=2, sort_keys=True)+'\n')
    return result


def check():
    record = json.loads(PROPOSAL.read_text())
    fresh = proposal()
    if fresh['candidate_profile_identity'] != record['candidate_profile_identity']:
        raise ValueError('rank-64 proposal identity changed')
    if fresh['predicted_all_node_residual_maxima'] != record['predicted_all_node_residual_maxima']:
        raise ValueError('rank-64 prediction changed')
    if fresh['step_norm'] != record['step_norm']:
        raise ValueError('rank-64 step norm changed')
    P.evolve_ks_difference_envelope = evolve_primal_controlled
    N16.RETARDED = I64.DIRECTIONS
    N16.DIRECTORY = DIRECTORY
    N16.OUTPUT = OUTPUT
    ctx = context(record)
    result = annotate(N16.assemble(ctx, replay=True), ctx)
    saved = json.loads(OUTPUT.read_text())
    if saved['constraint_maxima'] != result['constraint_maxima']:
        raise ValueError('rank-64 residual changed')
    if saved['prediction_minus_measured'] != result['prediction_minus_measured']:
        raise ValueError('rank-64 prediction error changed')
    if saved['physical_EXISTENCE_certificate'] or saved['physical_NONEXISTENCE_certificate']:
        raise ValueError('rank-64 measurement is not a gate certificate')
    return saved


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--propose', action='store_true')
    mode.add_argument('--run', action='store_true')
    mode.add_argument('--check', action='store_true')
    parser.add_argument('--cpu-budget', type=float, default=14400.0)
    parser.add_argument('--max-new', type=int, default=60)
    args = parser.parse_args()
    if args.propose:
        if PROPOSAL.exists():
            raise SystemExit('proposal already exists')
        record = proposal()
        PROPOSAL.write_text(json.dumps(record, indent=2, sort_keys=True)+'\n')
        print(json.dumps({
            'candidate_profile_identity': record['candidate_profile_identity'],
            'step_norm': record['step_norm'],
            'predicted_all_node_residual_maxima': record['predicted_all_node_residual_maxima'],
            'parent_constraint_maxima': record['parent_constraint_maxima'],
            'iterate6_constraint_maxima': record['iterate6_constraint_maxima'],
            'radius_lower_bound': record['radius_lower_bound'],
        }, indent=2))
    elif args.check:
        saved = check()
        print(json.dumps({
            'constraint_maxima': saved['constraint_maxima'],
            'prediction_minus_measured': saved['prediction_minus_measured'],
            'beats_iterate6_on_both_components': saved['beats_iterate6_on_both_components'],
        }, indent=2))
    else:
        result = run(args.cpu_budget, args.max_new)
        print(json.dumps({
            'status': result['status'],
            'completed_positive_families': result['completed_positive_families'],
            'constraint_maxima': result['constraint_maxima'],
            'prediction_minus_measured': result['prediction_minus_measured'],
            'beats_iterate6_on_both_components': result['beats_iterate6_on_both_components'],
            'profile_identity': result['profile_identity'],
        }, indent=2))
