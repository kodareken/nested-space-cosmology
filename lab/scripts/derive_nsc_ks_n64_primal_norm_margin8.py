#!/usr/bin/env python3
"""One joint step inside the measured quadratic-N regime.

The rejected norm-1e-7 step lowered beta and raised N. Its N error and the
accepted norm-1.3e-8 step both scale as the square of the step norm. In the
first 48 equilibrated modes a minimum-L1 step, pulled into the Euclidean ball
of radius 6.5e-9, predicts a beta drop of 3.5e-9 and an N drop at least eight
times that quadratic coefficient at the actual step norm.
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
DIRECTORY = ROOT/'results/development/artifacts/nsc-ks-coupled-newton-n64-primal-norm-margin8'
OUTPUT = ROOT/'results/development/nsc-ks-coupled-newton-n64-primal-norm-margin8.json'
PROPOSAL = ROOT/'results/development/nsc-ks-coupled-newton-n64-primal-norm-margin8-proposal.json'
OWNERS = (
    'scripts/derive_nsc_ks_n64_primal_norm_margin8.py',
    'docs/nsc-ks-n64-primal-norm-margin8.md',
    'src/recursive_horizons/nsc_ks_primal_step_control.py',
)
PARENT_IDENTITY = '6534e55b1f3dc3b968c2dfcdf56805113a5f4a354199e10a16704455b2316478'
KEPT = 48
N_MARGIN = 8.0
BETA_GAIN = 3.5e-9
MAX_STEP_NORM = 6.5e-9
MAX_RATIO = 2e3
MAX_CUTS = 8
SCALE = 1e12
EXISTENCE = 3e-11
OPTIONS = {
    'primal_feasibility_tolerance': 1e-10,
    'dual_feasibility_tolerance': 1e-10,
}


def digest(path):
    path = Path(path)
    return sha256((path if path.is_absolute() else ROOT/path).read_bytes()).hexdigest()


def quadratic_n_coefficient(parent):
    error = abs(float(parent['prediction_minus_measured'][0]))
    norm = float(parent['step_norm'])
    if error <= 0.0 or norm <= 0.0:
        raise ValueError('parent quadratic coefficient needs a measured N error and a step norm')
    return error / norm**2


def subspace(jacobian):
    residual_width = jacobian.shape[1]
    matrix = np.vstack([jacobian[:, :, 0].T, jacobian[:, :, 1].T])
    row = np.linalg.norm(matrix, axis=1)
    col = np.linalg.norm(matrix, axis=0)
    row = np.where(row > 0, row, 1.0)
    col = np.where(col > 0, col, 1.0)
    scaled = matrix / row[:, None] / col[None, :]
    _left, spectrum, right = np.linalg.svd(scaled, full_matrices=False)
    if spectrum.size < KEPT:
        raise ValueError('margin subspace is larger than the singular spectrum')
    ratio = float(spectrum[0] / spectrum[KEPT - 1])
    if ratio > MAX_RATIO:
        raise ValueError('margin subspace is past the well-conditioned block')
    basis = right[:KEPT].T / col[:, None]
    response = matrix @ basis
    return residual_width, matrix, basis, response, ratio


def linear_program(inequalities, limits):
    objective = np.concatenate([np.zeros(KEPT), np.ones(KEPT)])
    eye = np.eye(KEPT)
    pad = np.zeros((inequalities.shape[0], KEPT))
    augmented = np.vstack([
        np.hstack([inequalities, pad]),
        np.hstack([eye, -eye]),
        np.hstack([-eye, -eye]),
    ])
    bounds = np.concatenate([limits, np.zeros(2 * KEPT)])
    solved = [
        linprog(objective, A_ub=augmented, b_ub=bounds, bounds=(None, None),
                method='highs', options=OPTIONS)
        for _repeat in range(2)
    ]
    if not all(item.success for item in solved):
        return None
    if not np.array_equal(solved[0].x, solved[1].x):
        raise ValueError('margin linear program is not repeatable')
    return solved[0].x[:KEPT]


def margin_step(gradient, jacobian, coefficient):
    nodes, matrix, basis, response, ratio = subspace(jacobian)
    residual = np.concatenate([gradient[:, 0], gradient[:, 1]])
    measured = np.array([
        float(np.max(np.abs(residual[:nodes]))),
        float(np.max(np.abs(residual[nodes:]))),
    ])
    n_gain = N_MARGIN * coefficient * MAX_STEP_NORM**2
    if np.any(measured <= np.array([n_gain, BETA_GAIN])):
        raise ValueError('requested margin gains exceed the measured residual')
    caps = np.concatenate([
        np.full(nodes, measured[0] - n_gain),
        np.full(nodes, measured[1] - BETA_GAIN),
    ])
    scaled_response = response * SCALE
    scaled_residual = residual * SCALE
    scaled_caps = caps * SCALE
    rows = scaled_response.shape[0]
    base = np.zeros((2 * rows, KEPT))
    base[:rows] = scaled_response
    base[rows:] = -scaled_response
    base_limits = np.concatenate([
        scaled_caps - scaled_residual,
        scaled_caps + scaled_residual,
    ])
    cuts = []
    cut_limits = []
    for _cut in range(MAX_CUTS + 1):
        if cuts:
            inequalities = np.vstack([base, np.vstack(cuts)])
            limits = np.concatenate([base_limits, np.asarray(cut_limits, float)])
        else:
            inequalities = base
            limits = base_limits
        coefficients = linear_program(inequalities, limits)
        if coefficients is None:
            raise ValueError('margin step is infeasible in the Euclidean ball')
        step = basis @ coefficients
        norm = float(np.linalg.norm(step))
        if not np.isfinite(step).all() or not np.isfinite(norm):
            raise ValueError('margin step is not finite')
        if norm <= MAX_STEP_NORM * (1 + 1e-9):
            return step, residual, matrix, ratio, measured, n_gain, len(cuts)
        if len(cuts) == MAX_CUTS:
            break
        # One supporting cut of the Euclidean ball. L1 minimization alone can
        # sit outside the norm where the measured N nonlinearity dominates.
        gradient_coefficients = basis.T @ step
        cuts.append(gradient_coefficients / norm)
        cut_limits.append(MAX_STEP_NORM * norm)
    raise ValueError('margin step did not enter the Euclidean ball')


def proposal():
    parent = json.loads(PARENT.read_text())
    if parent['profile_identity'] != PARENT_IDENTITY:
        raise ValueError('margin parent must be the accepted beta-primary history')
    if parent['solver'] != 'primal-block DOP853 RMS on (A, D); tangent directions excluded':
        raise ValueError('margin parent Jacobian must come from primal-block control')
    if not parent['newton_step_accepted'] or not parent['residual_improved_on_all_components']:
        raise ValueError('margin parent must be an accepted improvement on both components')
    gradient = np.asarray(parent['action_gradient'], float)
    jacobian = np.asarray(parent['history_jacobian'], float)
    if jacobian.shape[0] != I64.DIRECTIONS:
        raise ValueError('parent Jacobian must carry 128 family-ordered directions')
    coefficient = quadratic_n_coefficient(parent)
    step, residual, matrix, ratio, measured, n_gain, cuts = margin_step(
        gradient, jacobian, coefficient)
    step_norm = float(np.linalg.norm(step))
    if step_norm > MAX_STEP_NORM * (1 + 1e-9):
        raise ValueError('margin step leaves the quadratic-N ball')
    predicted = I64.predicted_maxima(residual, matrix, step)
    gain = measured - predicted
    n_margin = float(gain[0] / (coefficient * step_norm**2))
    if gain[0] < n_gain * (1 - 1e-6) or gain[1] < BETA_GAIN * (1 - 1e-6):
        raise ValueError('margin prediction missed a requested gain')
    if n_margin < N_MARGIN:
        raise ValueError('actual N margin is below eight quadratic coefficients')
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
        raise ValueError('margin step must stay on the measured nodes')
    return {
        'schema': 'NSC-KS-N64-PRIMAL-NORM-MARGIN8-PROPOSAL-v1',
        'accountable_author': 'Douglas Ek',
        'status': 'OPEN: one quadratic-margin primal step authorized for measurement',
        'kept_modes': KEPT,
        'quadratic_n_coefficient': coefficient,
        'quadratic_n_coefficient_source': 'abs(parent prediction_minus_measured N) / parent step_norm^2',
        'n_margin_factor': N_MARGIN,
        'requested_n_gain': n_gain,
        'requested_beta_gain': BETA_GAIN,
        'max_step_norm': MAX_STEP_NORM,
        'euclidean_ball_cuts': cuts,
        'actual_n_margin_over_quadratic': n_margin,
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
        raise ValueError('margin evolution requires an authorized prediction')
    family = I64.LocalIncomingFamily64(np.asarray(record['candidate_coefficients'], float))
    identity = profile_identity(family, include_normal_window=True)
    if identity != record['candidate_profile_identity'] or identity in D.FORBIDDEN_IDENTITIES:
        raise ValueError('candidate identity changed or is forbidden')
    base = R.context()
    solve = family.collocation_nodes(N16.SOLVE_COUNT)
    verify = family.collocation_nodes(N16.VERIFY_COUNT)
    target = H.fixed_target_union(solve, verify)
    if target.tolist() != record['z']:
        raise ValueError('margin nodes changed')
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
        'margin8_proposal': record,
    }


def annotate(result, ctx):
    record = ctx['margin8_proposal']
    measured = np.asarray(result['constraint_maxima'], float)
    predicted = np.asarray(record['predicted_all_node_residual_maxima'], float)
    parent = np.asarray(record['parent_constraint_maxima'], float)
    complete = result['all_retained_nonzero_angular_families_included']
    gain = parent - predicted
    error = None if not complete else predicted - measured
    improved = bool(complete and np.all(measured < parent) and ctx['family'].radius_lower_bound() > 0)
    honest = bool(complete and np.all(np.abs(error) < gain))
    reaches = bool(complete and np.all(measured <= EXISTENCE))
    result['schema'] = 'NSC-KS-N64-PRIMAL-NORM-MARGIN8-v1'
    result['status'] = (
        'OPEN: measured quadratic-margin primal step; uncertified physical gate'
        if complete else 'OPEN: partial quadratic-margin primal step')
    result['solver'] = 'primal-block DOP853 RMS on (A, D); tangent directions excluded'
    result['kept_modes'] = KEPT
    result['n_margin_factor'] = N_MARGIN
    result['requested_beta_gain'] = BETA_GAIN
    result['quadratic_n_coefficient'] = record['quadratic_n_coefficient']
    result['actual_n_margin_over_quadratic'] = record['actual_n_margin_over_quadratic']
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
    result['reproducer'] = 'python3 scripts/derive_nsc_ks_n64_primal_norm_margin8.py --check'
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
        raise ValueError('margin proposal identity changed')
    if fresh['predicted_all_node_residual_maxima'] != record['predicted_all_node_residual_maxima']:
        raise ValueError('margin prediction changed')
    if fresh['step_norm'] != record['step_norm']:
        raise ValueError('margin step norm changed')
    if fresh['actual_n_margin_over_quadratic'] != record['actual_n_margin_over_quadratic']:
        raise ValueError('margin factor changed')
    P.evolve_ks_difference_envelope = evolve_primal_controlled
    N16.RETARDED = I64.DIRECTIONS
    N16.DIRECTORY = DIRECTORY
    N16.OUTPUT = OUTPUT
    ctx = context(record)
    result = annotate(N16.assemble(ctx, replay=True), ctx)
    saved = json.loads(OUTPUT.read_text())
    if saved['constraint_maxima'] != result['constraint_maxima']:
        raise ValueError('margin residual changed')
    if saved['prediction_minus_measured'] != result['prediction_minus_measured']:
        raise ValueError('margin prediction error changed')
    if saved['physical_EXISTENCE_certificate'] or saved['physical_NONEXISTENCE_certificate']:
        raise ValueError('margin measurement is not a gate certificate')
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
        again = proposal()
        if again['candidate_profile_identity'] != record['candidate_profile_identity']:
            raise SystemExit('proposal identity is not repeatable')
        if again['step'] != record['step']:
            raise SystemExit('proposal step is not repeatable')
        PROPOSAL.write_text(json.dumps(record, indent=2, sort_keys=True)+'\n')
        print(json.dumps({
            'candidate_profile_identity': record['candidate_profile_identity'],
            'step_norm': record['step_norm'],
            'singular_ratio': record['singular_ratio'],
            'actual_n_margin_over_quadratic': record['actual_n_margin_over_quadratic'],
            'euclidean_ball_cuts': record['euclidean_ball_cuts'],
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
