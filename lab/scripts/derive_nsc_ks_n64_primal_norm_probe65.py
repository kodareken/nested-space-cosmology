#!/usr/bin/env python3
"""Measure the quadratic N coefficient of equilibrated modes 65-73.

The step is the min-L1 direction in that nine-mode block, inside a Euclidean
ball of radius 1e-6. Authorization does not reuse the five-times coefficient
or the spent rank-56/64/70/74/80/88/94 rays.
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

PARENT = ROOT/'results/development/nsc-ks-coupled-newton-n64-primal-norm-trust03.json'
DIRECTORY = ROOT/'results/development/artifacts/nsc-ks-coupled-newton-n64-primal-norm-probe65'
OUTPUT = ROOT/'results/development/nsc-ks-coupled-newton-n64-primal-norm-probe65.json'
PROPOSAL = ROOT/'results/development/nsc-ks-coupled-newton-n64-primal-norm-probe65-proposal.json'
OWNERS = (
    'scripts/derive_nsc_ks_n64_primal_norm_probe65.py',
    'docs/nsc-ks-n64-primal-norm-probe65.md',
    'src/recursive_horizons/nsc_ks_primal_step_control.py',
)
PARENT_IDENTITY = 'ad75942438dad1de599ee260e5d9aa7f17d60af26b533bee2b0691f2a6dec777'
H_CAP = 1e-6
N_GAIN = 1.5e-13
BETA_GAIN = 1e-16
MODE_START = 64
MODE_END = 73
EXISTENCE = 3e-11
OPTIONS = {
    'primal_feasibility_tolerance': 1e-10,
    'dual_feasibility_tolerance': 1e-10,
}


def digest(path):
    path = Path(path)
    return sha256((path if path.is_absolute() else ROOT/path).read_bytes()).hexdigest()


def subspace(jacobian):
    matrix = np.vstack([jacobian[:, :, 0].T, jacobian[:, :, 1].T])
    row = np.linalg.norm(matrix, axis=1)
    col = np.linalg.norm(matrix, axis=0)
    row = np.where(row > 0, row, 1.0)
    col = np.where(col > 0, col, 1.0)
    scaled = matrix / row[:, None] / col[None, :]
    _left, spectrum, right = np.linalg.svd(scaled, full_matrices=False)
    raw = right[MODE_START:MODE_END].T / col[:, None]
    gram = raw.T @ raw
    evals, evecs = np.linalg.eigh(gram)
    keep = evals > evals[-1] * 1e-10
    basis = raw @ (evecs[:, keep] / np.sqrt(evals[keep]))
    ratio = float(spectrum[MODE_START] / spectrum[MODE_END - 1])
    return matrix, basis, matrix @ basis, ratio


def coefficients_of(response, residual, measured):
    nodes = residual.size // 2
    width = response.shape[1]
    scale = 1e6 / max(float(np.max(np.abs(residual))), 1e-30)
    rn = residual[:nodes] * scale
    rb = residual[nodes:] * scale
    an = response[:nodes] * (H_CAP * scale)
    ab = response[nodes:] * (H_CAP * scale)
    caps = np.array([measured[0] - N_GAIN, measured[1] - BETA_GAIN]) * scale
    cuts = []
    objective = np.concatenate([np.zeros(width), np.ones(width)])
    eye = np.eye(width)
    for _cut in range(16):
        gain_rows = np.vstack([an, -an, ab, -ab])
        gain_limits = np.concatenate([
            caps[0] - rn, caps[0] + rn, caps[1] - rb, caps[1] + rb,
        ])
        if cuts:
            gain_rows = np.vstack([gain_rows, np.vstack(cuts)])
            gain_limits = np.concatenate([gain_limits, np.full(len(cuts), 1.0)])
        pad = np.zeros((gain_rows.shape[0], width))
        augmented = np.vstack([
            np.hstack([gain_rows, pad]),
            np.hstack([eye, -eye]),
            np.hstack([-eye, -eye]),
        ])
        bounds = np.concatenate([gain_limits, np.zeros(2 * width)])
        solved = [
            linprog(objective, A_ub=augmented, b_ub=bounds, bounds=(None, None),
                    method='highs', options=OPTIONS)
            for _repeat in range(2)
        ]
        if not all(item.success for item in solved):
            return None
        if not np.array_equal(solved[0].x, solved[1].x):
            raise ValueError('probe linear program is not repeatable')
        direction = solved[0].x[:width]
        if float(np.linalg.norm(direction)) <= 1.0 + 1e-7:
            return direction
        cuts.append(direction / np.linalg.norm(direction))
    return None


def proposal():
    parent = json.loads(PARENT.read_text())
    if parent['profile_identity'] != PARENT_IDENTITY or not parent['newton_step_accepted']:
        raise ValueError('probe parent must be the accepted shrink')
    if parent['solver'] != 'primal-block DOP853 RMS on (A, D); tangent directions excluded':
        raise ValueError('probe parent Jacobian must come from primal-block control')
    gradient = np.asarray(parent['action_gradient'], float)
    jacobian = np.asarray(parent['history_jacobian'], float)
    if jacobian.shape[0] != I64.DIRECTIONS:
        raise ValueError('parent Jacobian must carry 128 family-ordered directions')
    measured = np.asarray(parent['constraint_maxima'], float)
    residual = np.concatenate([gradient[:, 0], gradient[:, 1]])
    matrix, basis, response, ratio = subspace(jacobian)
    direction = coefficients_of(response, residual, measured)
    if direction is None:
        raise ValueError('modes 65-73 probe is infeasible')
    step = H_CAP * (basis @ direction)
    step_norm = float(np.linalg.norm(step))
    if not 0.0 < step_norm <= H_CAP * (1.0 + 1e-8):
        raise ValueError('probe step leaves the Euclidean ball')
    predicted = I64.predicted_maxima(residual, matrix, step)
    gain = measured - predicted
    if gain[0] < N_GAIN * (1 - 1e-6) or gain[1] < BETA_GAIN * (1 - 1e-6):
        raise ValueError('probe prediction missed a requested gain')
    coefficients = np.asarray(parent['history']['coefficients'], float)
    candidate = I64.LocalIncomingFamily64(coefficients + step.reshape(2, 64))
    identity = profile_identity(candidate, include_normal_window=True)
    if identity in D.FORBIDDEN_IDENTITIES or identity == PARENT_IDENTITY:
        raise ValueError('candidate identity is forbidden or already measured')
    if candidate.radius_lower_bound() <= 0:
        raise ValueError('candidate radius bound must stay positive')
    solve = candidate.collocation_nodes(N16.SOLVE_COUNT)
    verify = candidate.collocation_nodes(N16.VERIFY_COUNT)
    target = H.fixed_target_union(solve, verify)
    if target.tolist() != parent['z']:
        raise ValueError('probe step must stay on the measured nodes')
    return {
        'schema': 'NSC-KS-N64-PRIMAL-NORM-PROBE65-PROPOSAL-v1',
        'accountable_author': 'Douglas Ek',
        'status': 'OPEN: one modes 65-73 probe authorized for its own N coefficient',
        'mode_start': MODE_START + 1,
        'mode_end': MODE_END,
        'subspace_singular_ratio': ratio,
        'max_step_norm': H_CAP,
        'requested_n_gain': N_GAIN,
        'requested_beta_gain': BETA_GAIN,
        'spent_truncations_not_reused': [56, 64, 70, 74, 80, 88, 94],
        'quadratic_coefficient_not_recycled': True,
        'step_norm': step_norm,
        'step': step.tolist(),
        'parent_profile_identity': parent['profile_identity'],
        'parent_constraint_maxima': measured.tolist(),
        'predicted_all_node_residual_maxima': predicted.tolist(),
        'predicted_gain': gain.tolist(),
        'prediction_is_not_a_nonlinear_residual': True,
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
        raise ValueError('probe evolution requires an authorized prediction')
    family = I64.LocalIncomingFamily64(np.asarray(record['candidate_coefficients'], float))
    identity = profile_identity(family, include_normal_window=True)
    if identity != record['candidate_profile_identity'] or identity in D.FORBIDDEN_IDENTITIES:
        raise ValueError('candidate identity changed or is forbidden')
    base = R.context()
    solve = family.collocation_nodes(N16.SOLVE_COUNT)
    verify = family.collocation_nodes(N16.VERIFY_COUNT)
    target = H.fixed_target_union(solve, verify)
    if target.tolist() != record['z']:
        raise ValueError('probe nodes changed')
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
        'probe65_proposal': record,
    }


def annotate(result, ctx):
    record = ctx['probe65_proposal']
    measured = np.asarray(result['constraint_maxima'], float)
    predicted = np.asarray(record['predicted_all_node_residual_maxima'], float)
    parent = np.asarray(record['parent_constraint_maxima'], float)
    complete = result['all_retained_nonzero_angular_families_included']
    gain = parent - predicted
    error = None if not complete else predicted - measured
    improved = bool(complete and np.all(measured < parent) and ctx['family'].radius_lower_bound() > 0)
    honest = bool(complete and np.all(np.abs(error) < gain))
    reaches = bool(complete and np.all(measured <= EXISTENCE))
    result['schema'] = 'NSC-KS-N64-PRIMAL-NORM-PROBE65-v1'
    result['status'] = (
        'OPEN: modes 65-73 probe measured; uncertified physical gate'
        if complete else 'OPEN: partial modes 65-73 probe')
    result['solver'] = 'primal-block DOP853 RMS on (A, D); tangent directions excluded'
    result['step_norm'] = record['step_norm']
    result['subspace_singular_ratio'] = record['subspace_singular_ratio']
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
    result['reproducer'] = 'python3 scripts/derive_nsc_ks_n64_primal_norm_probe65.py --check'
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
        raise ValueError('probe proposal identity changed')
    if fresh['predicted_all_node_residual_maxima'] != record['predicted_all_node_residual_maxima']:
        raise ValueError('probe prediction changed')
    if fresh['step_norm'] != record['step_norm']:
        raise ValueError('probe step norm changed')
    P.evolve_ks_difference_envelope = evolve_primal_controlled
    N16.RETARDED = I64.DIRECTIONS
    N16.DIRECTORY = DIRECTORY
    N16.OUTPUT = OUTPUT
    ctx = context(record)
    result = annotate(N16.assemble(ctx, replay=True), ctx)
    saved = json.loads(OUTPUT.read_text())
    if saved['constraint_maxima'] != result['constraint_maxima']:
        raise ValueError('probe residual changed')
    if saved['prediction_minus_measured'] != result['prediction_minus_measured']:
        raise ValueError('probe prediction error changed')
    if saved['physical_EXISTENCE_certificate'] or saved['physical_NONEXISTENCE_certificate']:
        raise ValueError('probe measurement is not a gate certificate')
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
        if again['step'] != record['step'] or again['candidate_profile_identity'] != record['candidate_profile_identity']:
            raise SystemExit('proposal is not repeatable')
        PROPOSAL.write_text(json.dumps(record, indent=2, sort_keys=True)+'\n')
        print(json.dumps({
            'candidate_profile_identity': record['candidate_profile_identity'],
            'step_norm': record['step_norm'],
            'subspace_singular_ratio': record['subspace_singular_ratio'],
            'predicted_all_node_residual_maxima': record['predicted_all_node_residual_maxima'],
            'predicted_gain': record['predicted_gain'],
            'radius_lower_bound': record['radius_lower_bound'],
        }, indent=2))
    elif args.check:
        saved = check()
        print(json.dumps({
            'constraint_maxima': saved['constraint_maxima'],
            'prediction_minus_measured': saved['prediction_minus_measured'],
            'newton_step_accepted': saved['newton_step_accepted'],
            'profile_identity': saved['profile_identity'],
        }, indent=2))
    else:
        result = run(args.cpu_budget, args.max_new)
        print(json.dumps({
            'status': result['status'],
            'completed_positive_families': result['completed_positive_families'],
            'constraint_maxima': result['constraint_maxima'],
            'prediction_minus_measured': result['prediction_minus_measured'],
            'newton_step_accepted': result['newton_step_accepted'],
            'profile_identity': result['profile_identity'],
        }, indent=2))
