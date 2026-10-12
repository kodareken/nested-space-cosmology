#!/usr/bin/env python3
"""Stall: fresh n=64 directions do not clear an N gain of 1e-10 with beta moving.

Spent truncations 56, 64, 70, 74, 80, 88 and 94 are not re-probed. No new
history is evolved. Quadratic coefficients stay unmeasured where the linear
ceiling is already below the prune bar for every nonnegative coefficient.
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
PROBE = ROOT/'results/development/nsc-ks-coupled-newton-n64-primal-norm-probe65.json'
OUTPUT = ROOT/'results/development/nsc-ks-coupled-newton-n64-primal-norm-fresh-bar.json'
OWNERS = (
    'scripts/derive_nsc_ks_n64_primal_norm_fresh_bar.py',
    'docs/nsc-ks-n64-primal-norm-fresh-bar.md',
)
PARENT_IDENTITY = 'ad75942438dad1de599ee260e5d9aa7f17d60af26b533bee2b0691f2a6dec777'
BAR = 1e-10
EXISTENCE = 3e-11
SPENT = (56, 64, 70, 74, 80, 88, 94)
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
    return raw @ (evecs[:, keep] / np.sqrt(evals[keep]))


def n_weighted_basis(matrix, col, count):
    weights = np.ones(matrix.shape[0])
    weights[:matrix.shape[0] // 2] = 1e3
    scaled = (matrix * weights[:, None]) / col[None, :]
    _left, _spectrum, right = np.linalg.svd(scaled, full_matrices=False)
    raw = right[:count].T / col[:, None]
    gram = raw.T @ raw
    evals, evecs = np.linalg.eigh(gram)
    keep = evals > evals[-1] * 1e-10
    return raw @ (evecs[:, keep] / np.sqrt(evals[keep]))


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
        solved = linprog(
            objective, A_ub=augmented, b_ub=bounds, bounds=(None, None),
            method='highs', options=OPTIONS)
        if not solved.success:
            return False, None
        direction = solved.x[:width]
        if float(np.linalg.norm(direction)) <= 1.0 + 1e-7:
            return True, direction
        cuts.append(direction / np.linalg.norm(direction))
    return False, None


def max_n_gain(response, residual, measured, h, beta_drop):
    lo = 0.0
    hi = float(measured[0])
    best = None
    for _step in range(24):
        mid = 0.5 * (lo + hi)
        ok, direction = feasible(
            response, residual, measured[0] - mid, measured[1] - beta_drop, h)
        if ok:
            lo = mid
            best = direction
        else:
            hi = mid
    return lo, best


def ray_ceiling(step, matrix, residual, measured):
    nodes = residual.size // 2
    step = np.asarray(step, float)
    step = step / max(float(np.linalg.norm(step)), 1e-30)
    response = matrix @ step
    best_gain = 0.0
    best_h = 0.0
    best_beta = 0.0
    for sign in (1.0, -1.0):
        column = sign * response
        for h in np.geomspace(1e-12, 1.0, 37):
            prediction = residual + float(h) * column
            maxima = np.array([
                float(np.max(np.abs(prediction[:nodes]))),
                float(np.max(np.abs(prediction[nodes:]))),
            ])
            gain = measured - maxima
            if gain[1] >= -1e-18 and gain[0] > best_gain:
                best_gain = float(gain[0])
                best_h = float(h)
                best_beta = float(gain[1])
    slope = best_gain / best_h if best_h else 0.0
    return best_gain, best_h, slope, best_beta


def window_row(name, basis, matrix, residual, measured, h, beta_drop):
    response = matrix @ basis
    gain, direction = max_n_gain(response, residual, measured, h, beta_drop)
    beta_gain = 0.0
    if direction is not None:
        step = h * (basis @ direction)
        predicted = I64.predicted_maxima(residual, matrix, step)
        beta_gain = float(measured[1] - predicted[1])
    slope = gain / h if h else 0.0
    return {
        'window': name,
        'width': int(basis.shape[1]),
        'h': h,
        'beta_drop_required': beta_drop,
        'n_gain': gain,
        'beta_gain': beta_gain,
        'slope': slope,
        'quad_coeff': None,
        'crossover_n_gain_upper_bound': gain,
        'c_required_for_bar': None if slope <= 0.0 else slope**2 / BAR,
        'clears_bar_with_beta_moving': bool(gain >= BAR and beta_gain > 1e-16),
        'probed': False,
    }


def measured_probe(parent):
    probe = json.loads(PROBE.read_text())
    if probe['profile_identity'] == parent['profile_identity'] or probe['newton_step_accepted']:
        raise ValueError('the modes 65-73 probe must stay a rejected new history')
    if probe['parent_profile_identity'] != PARENT_IDENTITY:
        raise ValueError('probe parent identity changed')
    if probe['solver'] != 'primal-block DOP853 RMS on (A, D); tangent directions excluded':
        raise ValueError('probe must use primal-block control')
    parent_max = parent['constraint_maxima']
    predicted = probe['linear_predicted_all_node_residual_maxima']
    measured = probe['constraint_maxima']
    if not (measured[0] > parent_max[0] and measured[1] < parent_max[1]):
        raise ValueError('the probe must raise N and lower beta')
    step_norm = float(probe['step_norm'])
    n_gain = float(parent_max[0] - predicted[0])
    beta_gain = float(parent_max[1] - predicted[1])
    n_error = abs(float(probe['prediction_minus_measured'][0]))
    if n_error <= n_gain or step_norm <= 0.0 or n_gain <= 0.0 or beta_gain <= 0.0:
        raise ValueError('the probe must have an N error above its predicted N gain')
    coefficient = n_error / step_norm**2
    slope = n_gain / step_norm
    crossover_norm = slope / coefficient
    crossover_gain = slope * crossover_norm
    if parent_max[0] - crossover_gain <= EXISTENCE:
        raise ValueError('probe crossover would require a certification review')
    return {
        'profile_identity': probe['profile_identity'],
        'constraint_maxima': measured,
        'prediction_minus_measured': probe['prediction_minus_measured'],
        'step_norm': step_norm,
        'subspace_singular_ratio': probe['subspace_singular_ratio'],
        'slope': slope,
        'quad_coeff': coefficient,
        'crossover_norm': crossover_norm,
        'crossover_n_gain': crossover_gain,
        'predicted_n_gain': n_gain,
        'predicted_beta_gain': beta_gain,
        'newton_step_accepted': False,
        'closes_beta_gap': False,
    }


def compute():
    parent = json.loads(PARENT.read_text())
    if parent['profile_identity'] != PARENT_IDENTITY or not parent['newton_step_accepted']:
        raise ValueError('fresh-bar parent must be the accepted shrink')
    if parent['solver'] != 'primal-block DOP853 RMS on (A, D); tangent directions excluded':
        raise ValueError('fresh-bar parent must use primal-block control')
    gradient = np.asarray(parent['action_gradient'], float)
    jacobian = np.asarray(parent['history_jacobian'], float)
    if jacobian.shape[0] != I64.DIRECTIONS:
        raise ValueError('parent Jacobian must carry 128 family-ordered directions')
    measured = np.asarray(parent['constraint_maxima'], float)
    residual, matrix, spectrum, right, col = equilibrated(gradient, jacobian)
    nodes = gradient.shape[0]
    if not np.allclose(residual[:nodes], gradient[:, 0]):
        raise ValueError('residual packing does not match the N column')
    windows = (
        ('modes_1_16', 0, 16, 1e-3),
        ('modes_17_32', 16, 32, 1e-3),
        ('modes_33_48', 32, 48, 1e-3),
        ('modes_1_40', 0, 40, 1e-3),
        ('modes_1_47', 0, 47, 1e-3),
        ('modes_1_47_small', 0, 47, 1e-6),
        ('modes_49_55', 48, 55, 1e-3),
        ('modes_1_73', 0, 73, 1e-3),
        ('modes_1_94', 0, 94, 1e-3),
        ('modes_65_73', 64, 73, 1e-2),
    )
    rows = []
    for name, start, end, h in windows:
        rows.append(window_row(
            name, whiten(right, col, start, end), matrix, residual, measured, h, 0.0))
    weighted = n_weighted_basis(matrix, col, 32)
    rows.append(window_row(
        'n_weighted_svd_32', weighted, matrix, residual, measured, 1e-3, 0.0))
    tail = next(row for row in rows if row['window'] == 'modes_65_73')
    tail_beta = window_row(
        'modes_65_73_beta_1e-12', whiten(right, col, 64, 73),
        matrix, residual, measured, 1e-2, 1e-12)
    rows.append(tail_beta)
    if tail['n_gain'] < BAR or tail_beta['n_gain'] < BAR:
        raise ValueError('modes 65-73 linear joint gain left the recorded regime')
    probe = measured_probe(parent)
    if probe['crossover_n_gain'] >= BAR:
        raise ValueError('measured modes 65-73 crossover clears the prune bar')
    plateau = next(row for row in rows if row['window'] == 'modes_1_47')
    plateau_small = next(row for row in rows if row['window'] == 'modes_1_47_small')
    if plateau['n_gain'] >= BAR or plateau_small['n_gain'] >= BAR:
        raise ValueError('the well-conditioned block reached the N bar')
    if plateau['n_gain'] > 2.0 * plateau_small['n_gain']:
        raise ValueError('well-conditioned N gain is still growing with step norm')
    best_mode = {'index': 0, 'n_gain': 0.0, 'h': 0.0, 'slope': 0.0, 'beta_gain': 0.0}
    for index in range(spectrum.size):
        if (index + 1) in SPENT:
            continue
        basis = whiten(right, col, index, index + 1)
        gain, h, slope, beta = ray_ceiling(basis[:, 0], matrix, residual, measured)
        if gain > best_mode['n_gain']:
            best_mode = {
                'index': index + 1, 'n_gain': gain, 'h': h,
                'slope': slope, 'beta_gain': beta,
            }
    best_coordinate = {'index': 0, 'n_gain': 0.0, 'h': 0.0, 'slope': 0.0, 'beta_gain': 0.0}
    for index in range(matrix.shape[1]):
        step = np.zeros(matrix.shape[1])
        step[index] = 1.0
        gain, h, slope, beta = ray_ceiling(step, matrix, residual, measured)
        if gain > best_coordinate['n_gain']:
            best_coordinate = {
                'index': index, 'n_gain': gain, 'h': h,
                'slope': slope, 'beta_gain': beta,
            }
    if best_mode['n_gain'] >= BAR or best_coordinate['n_gain'] >= BAR:
        raise ValueError('a single fresh ray reaches the N bar')
    n_gap = float(measured[0] - EXISTENCE)
    beta_gap = float(measured[1] - EXISTENCE)
    return {
        'schema': 'NSC-KS-N64-PRIMAL-NORM-FRESH-BAR-v1',
        'accountable_author': 'Douglas Ek',
        'status': 'OPEN: the one probed fresh direction has crossover N gain 1.71e-16',
        'stall': 'n64_primal_norm_fresh_directions_below_1e-10',
        'best_profile_identity': parent['profile_identity'],
        'best_constraint_maxima': [float(measured[0]), float(measured[1])],
        'n_gap_to_existence': n_gap,
        'beta_gap_to_existence': beta_gap,
        'stated_7p29e-8_drop_is_the_beta_gap': True,
        'prune_bar_on_n_gain': BAR,
        'spent_truncations_not_reprobed': list(SPENT),
        'directions_screened': rows,
        'best_nonspent_svd_mode': best_mode,
        'best_coordinate_ray': best_coordinate,
        'modes_65_73_linear_n_gain_at_h_0_01': tail_beta['n_gain'],
        'modes_65_73_linear_beta_gain_at_h_0_01': tail_beta['beta_gain'],
        'probe': probe,
        'probes_evolved': 1,
        'families_evolved': 60,
        'evolve_authorized': False,
        'multi_step_composition': False,
        'higher_n_pathway': 'LocalIncomingFamily allows 8, 16, 32 or 64 coefficients; no n>64 assembled owner exists',
        'physical_EXISTENCE_certificate': False,
        'physical_NONEXISTENCE_certificate': False,
        'nonexistence_not_claimed': 'joint linear gain remains positive and the geometry cokernel stays the recorded trivial one',
        'enclosed_error': None,
        'uv_field_low_subgap_between_node': None,
        'next_owner': 'declared (w,U) class on I=S(1)+[0.12,0.18]; next resolution is an owned Chebyshev block above 64 from this identity',
    }


def record():
    payload = compute()
    payload['source_hashes'] = {path: digest(path) for path in OWNERS}
    payload['input_hashes'] = {
        str(PARENT.relative_to(ROOT)): digest(PARENT),
        str(PROBE.relative_to(ROOT)): digest(PROBE),
    }
    return payload


def check():
    saved = json.loads(OUTPUT.read_text())
    fresh = record()
    for key in (
        'stall', 'best_constraint_maxima', 'n_gap_to_existence',
        'beta_gap_to_existence', 'directions_screened',
        'best_nonspent_svd_mode', 'best_coordinate_ray',
        'probe', 'probes_evolved',
        'physical_NONEXISTENCE_certificate',
    ):
        if saved[key] != fresh[key]:
            raise ValueError('fresh-bar register changed: '+key)
    return saved


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    if args.write:
        if OUTPUT.exists():
            raise SystemExit('fresh-bar register already exists')
        payload = record()
        OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+'\n')
        summary = [{
            'window': row['window'],
            'h': row['h'],
            'n_gain': row['n_gain'],
            'beta_gain': row['beta_gain'],
            'slope': row['slope'],
        } for row in payload['directions_screened']]
        print(json.dumps({
            'stall': payload['stall'],
            'n_gap_to_existence': payload['n_gap_to_existence'],
            'beta_gap_to_existence': payload['beta_gap_to_existence'],
            'rows': summary,
            'best_nonspent_svd_mode': payload['best_nonspent_svd_mode'],
            'probe_crossover_n_gain': payload['probe']['crossover_n_gain'],
            'probe_quad_coeff': payload['probe']['quad_coeff'],
            'probes_evolved': payload['probes_evolved'],
        }, indent=2))
    else:
        saved = check()
        print(json.dumps({'stall': saved['stall'], 'status': saved['status']}, indent=2))
