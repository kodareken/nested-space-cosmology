#!/usr/bin/env python3
"""One more beta-primary step from the measured 1e-8 beta history.

The parent prediction matched to a few picounits. The same 48-mode block
can lower beta by another 3e-9 and N by five times that new N error.
Measurement uses primal-block DOP853 control.
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
from scipy.optimize import linprog

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

PARENT = ROOT/'results/development/nsc-ks-coupled-newton-n64-primal-norm-beta48.json'
DIRECTORY = ROOT/'results/development/artifacts/nsc-ks-coupled-newton-n64-primal-norm-beta48b'
OUTPUT = ROOT/'results/development/nsc-ks-coupled-newton-n64-primal-norm-beta48b.json'
PROPOSAL = ROOT/'results/development/nsc-ks-coupled-newton-n64-primal-norm-beta48b-proposal.json'
OWNERS = (
    'scripts/derive_nsc_ks_n64_primal_norm_beta48b.py',
    'docs/nsc-ks-n64-primal-norm-beta48b.md',
    'src/recursive_horizons/nsc_ks_primal_step_control.py',
)
PARENT_IDENTITY = '6534e55b1f3dc3b968c2dfcdf56805113a5f4a354199e10a16704455b2316478'
KEPT = 48
N_GAIN_FACTOR = 5.0
BETA_GAIN = 3e-9
MAX_STEP_NORM = 1e-6
MAX_RATIO = 1e4
SCALE = 1e12
EXISTENCE = 3e-11


def digest(path):
    path = Path(path)
    return sha256((path if path.is_absolute() else ROOT/path).read_bytes()).hexdigest()


def sup_step(gradient, jacobian, floor):
    _delta, residual, matrix, singular = I64.rank_truncated_delta(gradient, jacobian, 1)
    row = np.linalg.norm(matrix, axis=1)
    col = np.linalg.norm(matrix, axis=0)
    row = np.where(row > 0, row, 1.0)
    col = np.where(col > 0, col, 1.0)
    scaled = matrix / row[:, None] / col[None, :]
    _left, spectrum, right = np.linalg.svd(scaled, full_matrices=False)
    if spectrum.size < KEPT:
        raise ValueError('sup-norm subspace is larger than the singular spectrum')
    ratio = float(spectrum[0] / spectrum[KEPT - 1])
    if ratio > MAX_RATIO:
        raise ValueError('sup-norm subspace is past the well-conditioned block')
    measured = np.array([
        float(np.max(np.abs(residual[:gradient.shape[0]]))),
        float(np.max(np.abs(residual[gradient.shape[0]:]))),
    ])
    floor = np.asarray(floor, float)
    gain = np.array([N_GAIN_FACTOR * floor[0], BETA_GAIN])
    if np.any(measured <= gain):
        raise ValueError('requested sup-norm gain exceeds the measured residual')
    basis = right[:KEPT].T / col[:, None]
    response = matrix @ basis
    nodes = gradient.shape[0]
    bounds = np.concatenate([
        np.full(nodes, (measured[0] - gain[0]) * SCALE),
        np.full(nodes, (measured[1] - gain[1]) * SCALE),
    ])
    scaled_residual = residual * SCALE
    scaled_response = response * SCALE
    rows = scaled_response.shape[0]
    objective = np.concatenate([np.zeros(KEPT), np.ones(KEPT)])
    inequalities = np.zeros((2 * rows + 2 * KEPT, 2 * KEPT))
    inequalities[:rows, :KEPT] = scaled_response
    inequalities[rows:2 * rows, :KEPT] = -scaled_response
    inequalities[2 * rows:2 * rows + KEPT, :KEPT] = np.eye(KEPT)
    inequalities[2 * rows:2 * rows + KEPT, KEPT:] = -np.eye(KEPT)
    inequalities[2 * rows + KEPT:, :KEPT] = -np.eye(KEPT)
    inequalities[2 * rows + KEPT:, KEPT:] = -np.eye(KEPT)
    limits = np.concatenate([
        bounds - scaled_residual,
        bounds + scaled_residual,
        np.zeros(2 * KEPT),
    ])
    options = {
        'primal_feasibility_tolerance': 1e-10,
        'dual_feasibility_tolerance': 1e-10,
    }
    solved = [
        linprog(objective, A_ub=inequalities, b_ub=limits, bounds=(None, None),
                method='highs', options=options)
        for _repeat in range(2)
    ]
    if not all(item.success for item in solved):
        raise ValueError('sup-norm step is infeasible in the well-conditioned block')
    coefficients = [item.x[:KEPT] for item in solved]
    if not np.array_equal(coefficients[0], coefficients[1]):
        raise ValueError('sup-norm linear program is not repeatable')
    step = basis @ coefficients[0]
    if not np.isfinite(step).all():
        raise ValueError('sup-norm step is not finite')
    return step, residual, matrix, spectrum, ratio, measured


def proposal():
    parent = json.loads(PARENT.read_text())
    if parent['profile_identity'] != PARENT_IDENTITY:
        raise ValueError('beta-primary parent must be the measured 1e-8 beta history')
    if parent['solver'] != 'primal-block DOP853 RMS on (A, D); tangent directions excluded':
        raise ValueError('sup-norm parent Jacobian must come from primal-block control')
    if not parent['newton_step_accepted'] or not parent['residual_improved_on_all_components']:
        raise ValueError('sup-norm parent must be an accepted improvement on both components')
    gradient = np.asarray(parent['action_gradient'], float)
    jacobian = np.asarray(parent['history_jacobian'], float)
    if jacobian.shape[0] != I64.DIRECTIONS:
        raise ValueError('parent Jacobian must carry 128 family-ordered directions')
    floor = np.abs(np.asarray(parent['prediction_minus_measured'], float))
    step, residual, matrix, _spectrum, ratio, measured = sup_step(gradient, jacobian, floor)
    step_norm = float(np.linalg.norm(step))
    if step_norm > MAX_STEP_NORM:
        raise ValueError('sup-norm step leaves the previously linear regime')
    predicted = I64.predicted_maxima(residual, matrix, step)
    gain = measured - predicted
    requested = np.array([N_GAIN_FACTOR * floor[0], BETA_GAIN])
    if not np.all(gain >= requested * (1 - 1e-6)):
        raise ValueError('beta-primary prediction must clear the requested gains')
    coefficients = np.asarray(parent['history']['coefficients'], float)
    candidate = I64.LocalIncomingFamily64(coefficients + step.reshape(2, 64))
    identity = profile_identity(candidate, include_normal_window=True)
    if identity in D.FORBIDDEN_IDENTITIES or identity == parent['profile_identity']:
        raise ValueError('candidate identity is forbidden or unchanged')
    if candidate.radius_lower_bound() <= 0:
        raise ValueError('candidate radius bound must stay positive')
    solve = candidate.collocation_nodes(N16.SOLVE_COUNT)
    verify = candidate.collocation_nodes(N16.VERIFY_COUNT)
    target = H.fixed_target_union(solve, verify)
    if target.tolist() != parent['z']:
        raise ValueError('sup-norm step must stay on the measured nodes')
    return {
        'schema': 'NSC-KS-N64-PRIMAL-NORM-BETA48B-PROPOSAL-v1',
        'accountable_author': 'Douglas Ek',
        'status': 'OPEN: one sup-norm primal step authorized for measurement',
        'kept_modes': KEPT,
        'n_gain_factor_over_prediction_error': N_GAIN_FACTOR,
        'requested_beta_gain': BETA_GAIN,
        'prediction_error_floor': floor.tolist(),
        'singular_ratio': ratio,
        'step_norm': step_norm,
        'step': step.tolist(),
        'parent_profile_identity': parent['profile_identity'],
        'parent_constraint_maxima': measured.tolist(),
        'predicted_all_node_residual_maxima': predicted.tolist(),
        'predicted_gain': gain.tolist(),
        'prediction_is_not_a_nonlinear_residual': True,
        'ill_conditioned_tail_not_evolved': 'ranks from 74 have norms 0.30-0.79 and singular ratios from 6.7e8',
        'candidate_profile_identity': identity,
        'candidate_coefficients': candidate.description()['coefficients'],
        'radius_lower_bound': candidate.radius_lower_bound(),
        'z': parent['z'],
        'evolve_authorized': True,
        'physical_EXISTENCE_certificate': False,
        'physical_NONEXISTENCE_certificate': False,
        'source_hashes': {path: digest(path) for path in OWNERS},
        'input_hashes': {str(PARENT.relative_to(ROOT)): digest(PARENT)},
    }


def context(record):
    if not record.get('evolve_authorized'):
        raise ValueError('sup-norm evolution requires an authorized prediction')
    family = I64.LocalIncomingFamily64(np.asarray(record['candidate_coefficients'], float))
    identity = profile_identity(family, include_normal_window=True)
    if identity != record['candidate_profile_identity'] or identity in D.FORBIDDEN_IDENTITIES:
        raise ValueError('candidate identity changed or is forbidden')
    base = R.context()
    solve = family.collocation_nodes(N16.SOLVE_COUNT)
    verify = family.collocation_nodes(N16.VERIFY_COUNT)
    target = H.fixed_target_union(solve, verify)
    if target.tolist() != record['z']:
        raise ValueError('sup-norm nodes changed')
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
        'selected': {'step_scale': 1.0, 'step_norm': record['step_norm']},
        'beta48b_proposal': record,
    }


def annotate(result, ctx):
    record = ctx['beta48b_proposal']
    measured = np.asarray(result['constraint_maxima'], float)
    predicted = np.asarray(record['predicted_all_node_residual_maxima'], float)
    parent = np.asarray(record['parent_constraint_maxima'], float)
    complete = result['all_retained_nonzero_angular_families_included']
    gain = parent - predicted
    error = None if not complete else predicted - measured
    improved = bool(complete and np.all(measured < parent) and ctx['family'].radius_lower_bound() > 0)
    honest = bool(complete and np.all(np.abs(error) < gain))
    reaches = bool(complete and np.all(measured <= EXISTENCE))
    result['schema'] = 'NSC-KS-N64-PRIMAL-NORM-BETA48B-v1'
    result['status'] = (
        'OPEN: measured sup-norm primal step; uncertified physical gate'
        if complete else 'OPEN: partial sup-norm primal step')
    result['solver'] = 'primal-block DOP853 RMS on (A, D); tangent directions excluded'
    result['kept_modes'] = KEPT
    result['n_gain_factor_over_prediction_error'] = N_GAIN_FACTOR
    result['requested_beta_gain'] = BETA_GAIN
    result['step_norm'] = record['step_norm']
    result['singular_ratio'] = record['singular_ratio']
    result['linear_predicted_all_node_residual_maxima'] = predicted.tolist()
    result['prediction_minus_measured'] = None if error is None else error.tolist()
    result['prediction_minus_measured_is_honest'] = honest if complete else None
    result['parent_constraint_maxima'] = parent.tolist()
    result['parent_profile_identity'] = record['parent_profile_identity']
    result['residual_improved_on_all_components'] = improved if complete else None
    result['prediction_is_not_a_nonlinear_residual'] = True
    result['coefficient_step_evolved'] = True
    result['newton_step_accepted'] = bool(improved and honest)
    result['physical_EXISTENCE_certificate'] = False
    result['physical_NONEXISTENCE_certificate'] = False
    result['nodal_residual_reaches_existence_tolerance'] = reaches if complete else None
    result['enclosed_error'] = None
    result['reproducer'] = 'python3 scripts/derive_nsc_ks_n64_primal_norm_beta48b.py --check'
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
        raise ValueError('sup-norm proposal identity changed')
    if fresh['predicted_all_node_residual_maxima'] != record['predicted_all_node_residual_maxima']:
        raise ValueError('sup-norm prediction changed')
    if fresh['step_norm'] != record['step_norm']:
        raise ValueError('sup-norm step norm changed')
    P.evolve_ks_difference_envelope = evolve_primal_controlled
    N16.RETARDED = I64.DIRECTIONS
    N16.DIRECTORY = DIRECTORY
    N16.OUTPUT = OUTPUT
    ctx = context(record)
    result = annotate(N16.assemble(ctx, replay=True), ctx)
    saved = json.loads(OUTPUT.read_text())
    if saved['constraint_maxima'] != result['constraint_maxima']:
        raise ValueError('sup-norm residual changed')
    if saved['prediction_minus_measured'] != result['prediction_minus_measured']:
        raise ValueError('sup-norm prediction error changed')
    if saved['physical_EXISTENCE_certificate'] or saved['physical_NONEXISTENCE_certificate']:
        raise ValueError('sup-norm measurement is not a gate certificate')
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
            'singular_ratio': record['singular_ratio'],
            'predicted_all_node_residual_maxima': record['predicted_all_node_residual_maxima'],
            'predicted_gain': record['predicted_gain'],
            'parent_constraint_maxima': record['parent_constraint_maxima'],
            'radius_lower_bound': record['radius_lower_bound'],
        }, indent=2))
    elif args.check:
        saved = check()
        print(json.dumps({
            'constraint_maxima': saved['constraint_maxima'],
            'prediction_minus_measured': saved['prediction_minus_measured'],
            'residual_improved_on_all_components': saved['residual_improved_on_all_components'],
            'prediction_minus_measured_is_honest': saved['prediction_minus_measured_is_honest'],
            'newton_step_accepted': saved['newton_step_accepted'],
        }, indent=2))
    else:
        result = run(args.cpu_budget, args.max_new)
        print(json.dumps({
            'status': result['status'],
            'completed_positive_families': result['completed_positive_families'],
            'constraint_maxima': result['constraint_maxima'],
            'prediction_minus_measured': result['prediction_minus_measured'],
            'residual_improved_on_all_components': result['residual_improved_on_all_components'],
            'prediction_minus_measured_is_honest': result['prediction_minus_measured_is_honest'],
            'newton_step_accepted': result['newton_step_accepted'],
            'profile_identity': result['profile_identity'],
        }, indent=2))
