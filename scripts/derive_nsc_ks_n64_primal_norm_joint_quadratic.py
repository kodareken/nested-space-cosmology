#!/usr/bin/env python3
"""Stall: best joint linear N gain sits inside the measured quadratic term.

The accepted shrink at norm 2.586e-10 is the largest step whose N prediction
error stayed near 6e-16. Rank 74 and above, and the unused mid-band rays,
are tested on that Jacobian. None is evolved. The comparison coefficient is
the measured five-times step, not a coefficient of these unevolved directions.
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
from scipy.optimize import linprog

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
import derive_nsc_ks_coupled_newton_n64_truncated_iterate as I64

PARENT = ROOT/'results/development/nsc-ks-coupled-newton-n64-primal-norm-trust03.json'
QUADRATIC = ROOT/'results/development/nsc-ks-coupled-newton-n64-primal-norm-trust05-quadratic.json'
OUTPUT = ROOT/'results/development/nsc-ks-coupled-newton-n64-primal-norm-joint-quadratic.json'
OWNERS = (
    'scripts/derive_nsc_ks_n64_primal_norm_joint_quadratic.py',
    'docs/nsc-ks-n64-primal-norm-joint-quadratic.md',
)
PARENT_IDENTITY = 'ad75942438dad1de599ee260e5d9aa7f17d60af26b533bee2b0691f2a6dec777'
BETA_GAIN = 1e-10
EXISTENCE = 3e-11
DAMPED_RANKS = (56, 64, 70, 74, 80, 88, 94)
OPTIONS = {
    'primal_feasibility_tolerance': 1e-10,
    'dual_feasibility_tolerance': 1e-10,
}


def digest(path):
    path = Path(path)
    return sha256((path if path.is_absolute() else ROOT/path).read_bytes()).hexdigest()


def equilibrated(gradient, jacobian):
    residual = np.concatenate([gradient[:, 0], gradient[:, 1]])
    matrix = np.vstack([jacobian[:, :, 0].T, jacobian[:, :, 1].T])
    row = np.linalg.norm(matrix, axis=1)
    col = np.linalg.norm(matrix, axis=0)
    row = np.where(row > 0, row, 1.0)
    col = np.where(col > 0, col, 1.0)
    scaled = matrix / row[:, None] / col[None, :]
    _left, spectrum, right = np.linalg.svd(scaled, full_matrices=False)
    return residual, matrix, spectrum, right, col


def whiten(right, col, start, end):
    raw = right[start:end].T / col[:, None]
    gram = raw.T @ raw
    evals, evecs = np.linalg.eigh(gram)
    keep = evals > evals[-1] * 1e-10
    basis = raw @ (evecs[:, keep] / np.sqrt(evals[keep]))
    return basis


def feasible(response, residual, n_cap, beta_cap, h):
    nodes = residual.size // 2
    width = response.shape[1]
    scale = 1e6 / max(float(np.max(np.abs(residual))), 1e-30)
    rn = residual[:nodes] * scale
    rb = residual[nodes:] * scale
    an = response[:nodes] * (h * scale)
    ab = response[nodes:] * (h * scale)
    caps = np.array([n_cap, beta_cap], float) * scale
    cuts = []
    objective = np.concatenate([np.zeros(width), np.ones(width)])
    eye = np.eye(width)
    for _cut in range(20):
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
        solved = linprog(
            objective, A_ub=augmented, b_ub=bounds, bounds=(None, None),
            method='highs', options=OPTIONS)
        if not solved.success:
            return False
        direction = solved.x[:width]
        if float(np.linalg.norm(direction)) <= 1.0 + 1e-7:
            return True
        cuts.append(direction / np.linalg.norm(direction))
    return False


def max_n_gain(response, residual, measured, h, beta_gain):
    lo = 0.0
    hi = float(measured[0])
    for _step in range(32):
        mid = 0.5 * (lo + hi)
        if feasible(response, residual, measured[0] - mid, measured[1] - beta_gain, h):
            lo = mid
        else:
            hi = mid
    return lo


def damped_rows(gradient, jacobian, measured, h):
    rows = []
    for rank in DAMPED_RANKS:
        delta, residual, matrix, singular = I64.rank_truncated_delta(gradient, jacobian, rank)
        raw = float(np.linalg.norm(delta))
        step = delta * (h / raw)
        predicted = I64.predicted_maxima(residual, matrix, step)
        gain = measured - predicted
        rows.append({
            'kept_modes': rank,
            'raw_step_norm': raw,
            'singular_ratio': float(singular[0] / singular[rank - 1]),
            'damped_step_norm': h,
            'predicted_gain': [float(gain[0]), float(gain[1])],
            'both_components_improve': bool(gain[0] > 0 and gain[1] > 0),
        })
    return rows


def crossover(slope, coefficient, measured_n):
    h_star = slope / coefficient
    gain = slope * h_star
    return {
        'step_norm': h_star,
        'predicted_n_gain': gain,
        'residual_after_gain': measured_n - gain,
        'can_reach_existence': bool(measured_n - gain <= EXISTENCE),
    }


def compute():
    parent = json.loads(PARENT.read_text())
    quadratic = json.loads(QUADRATIC.read_text())
    if parent['profile_identity'] != PARENT_IDENTITY or not parent['newton_step_accepted']:
        raise ValueError('joint search parent must be the accepted shrink')
    if parent['solver'] != 'primal-block DOP853 RMS on (A, D); tangent directions excluded':
        raise ValueError('joint search parent must use primal-block control')
    if quadratic['stall'] != 'n64_primal_norm_trust05_n_rose':
        raise ValueError('quadratic reference must stay the rejected five-times step')
    if quadratic['best_profile_identity'] != PARENT_IDENTITY:
        raise ValueError('quadratic reference must still name the accepted shrink')
    gradient = np.asarray(parent['action_gradient'], float)
    jacobian = np.asarray(parent['history_jacobian'], float)
    measured = np.asarray(parent['constraint_maxima'], float)
    if jacobian.shape[0] != I64.DIRECTIONS:
        raise ValueError('parent Jacobian must carry 128 family-ordered directions')
    h = float(parent['step_norm'])
    coefficient = float(quadratic['direction_quadratic_n_coefficient'])
    spoilage = coefficient * h * h
    residual, matrix, spectrum, right, col = equilibrated(gradient, jacobian)
    nodes = gradient.shape[0]
    if not np.allclose(residual[:nodes], gradient[:, 0]):
        raise ValueError('residual packing does not match the N column')
    head = whiten(right, col, 0, 48)
    tail = whiten(right, col, 73, 94)
    head_response = matrix @ head
    tail_response = matrix @ tail
    held = max_n_gain(head_response, residual, measured, h, 0.0)
    joint = max_n_gain(head_response, residual, measured, h, BETA_GAIN)
    if not (0.0 < joint <= held):
        raise ValueError('joint N gain must be positive and no larger than the beta-held gain')
    if held >= spoilage or joint >= spoilage:
        raise ValueError('linear N gain clears the measured quadratic term; do not record this stall')
    tail_lip = float(np.max(np.linalg.norm(tail_response[:nodes], axis=1)))
    head_slope = joint / h
    head_cross = crossover(head_slope, coefficient, float(measured[0]))
    tail_cross = crossover(tail_lip, coefficient, float(measured[0]))
    if head_cross['can_reach_existence'] or tail_cross['can_reach_existence']:
        raise ValueError('crossover gain would require a certification review')
    half_crossover_net = head_cross['predicted_n_gain'] / 4.0
    if half_crossover_net >= 1e-11:
        raise ValueError('half-crossover net N gain is large enough to reconsider the stall')
    rows = damped_rows(gradient, jacobian, measured, h)
    uphill = [row for row in rows if row['kept_modes'] in (56, 64, 70)]
    if not all(row['predicted_gain'][0] < 0 and row['predicted_gain'][1] < 0 for row in uphill):
        raise ValueError('mid-rank damped rays must raise both residuals at the linear norm')
    tail_rays = [row for row in rows if row['kept_modes'] >= 74]
    if not all(row['predicted_gain'][0] < 1e-16 for row in tail_rays):
        raise ValueError('damped rank>=74 N gain is no longer negligible at the linear norm')
    if any(row['raw_step_norm'] < 0.1 for row in tail_rays):
        raise ValueError('rank>=74 raw step entered the measured linear norm')
    return {
        'schema': 'NSC-KS-N64-PRIMAL-NORM-JOINT-QUADRATIC-v1',
        'accountable_author': 'Douglas Ek',
        'status': 'OPEN: best joint linear N gain is inside the measured quadratic term',
        'stall': 'n64_primal_norm_joint_gain_inside_quadratic',
        'best_profile_identity': parent['profile_identity'],
        'best_constraint_maxima': [float(measured[0]), float(measured[1])],
        'best_prediction_minus_measured': parent['prediction_minus_measured'],
        'linear_norm': h,
        'comparison_direction': 'rejected five-times step baf90351; coefficient is not remeasured on the unevolved directions',
        'direction_quadratic_n_coefficient': coefficient,
        'quadratic_spoilage_at_linear_norm': spoilage,
        'modes_1_48_width': int(head.shape[1]),
        'modes_1_48_singular_ratio': float(spectrum[0] / spectrum[47]),
        'n_gain_with_beta_held': held,
        'required_beta_gain': BETA_GAIN,
        'n_gain_with_required_beta_gain': joint,
        'joint_n_slope': head_slope,
        'joint_crossover_against_measured_coefficient': head_cross,
        'half_crossover_net_n_gain': half_crossover_net,
        'modes_74_94_width': int(tail.shape[1]),
        'modes_74_94_n_lipschitz': tail_lip,
        'tail_crossover_against_measured_coefficient': tail_cross,
        'damped_truncations_at_linear_norm': rows,
        'families_evolved': 0,
        'evolve_authorized': False,
        'newton_step_accepted': False,
        'physical_EXISTENCE_certificate': False,
        'physical_NONEXISTENCE_certificate': False,
        'nonexistence_not_claimed': 'positive joint linear gain and the already recorded trivial geometry cokernel do not exclude the declared class',
        'enclosed_error': None,
        'uv_field_low_subgap_between_node': None,
        'next_owner': 'declared (w,U) class on I=S(1)+[0.12,0.18]',
    }


def record():
    payload = compute()
    payload['source_hashes'] = {path: digest(path) for path in OWNERS}
    payload['input_hashes'] = {
        str(PARENT.relative_to(ROOT)): digest(PARENT),
        str(QUADRATIC.relative_to(ROOT)): digest(QUADRATIC),
    }
    return payload


def check():
    saved = json.loads(OUTPUT.read_text())
    fresh = record()
    for key in (
        'stall', 'best_constraint_maxima', 'linear_norm',
        'direction_quadratic_n_coefficient', 'quadratic_spoilage_at_linear_norm',
        'n_gain_with_beta_held', 'n_gain_with_required_beta_gain',
        'joint_crossover_against_measured_coefficient',
        'half_crossover_net_n_gain',
        'tail_crossover_against_measured_coefficient',
        'damped_truncations_at_linear_norm',
        'evolve_authorized', 'physical_NONEXISTENCE_certificate',
    ):
        if saved[key] != fresh[key]:
            raise ValueError('joint quadratic register changed: '+key)
    return saved


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    if args.write:
        if OUTPUT.exists():
            raise SystemExit('joint quadratic register already exists')
        payload = record()
        OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+'\n')
        print(json.dumps({
            'stall': payload['stall'],
            'n_gain_with_beta_held': payload['n_gain_with_beta_held'],
            'n_gain_with_required_beta_gain': payload['n_gain_with_required_beta_gain'],
            'quadratic_spoilage_at_linear_norm': payload['quadratic_spoilage_at_linear_norm'],
            'joint_crossover_against_measured_coefficient': payload['joint_crossover_against_measured_coefficient'],
            'tail_crossover_against_measured_coefficient': payload['tail_crossover_against_measured_coefficient'],
            'damped_rank74_gain': payload['damped_truncations_at_linear_norm'][3]['predicted_gain'],
            'damped_rank74_raw_norm': payload['damped_truncations_at_linear_norm'][3]['raw_step_norm'],
        }, indent=2))
    else:
        saved = check()
        print(json.dumps({'stall': saved['stall'], 'status': saved['status']}, indent=2))
