#!/usr/bin/env python3
"""One physical-metric trust-region step from nodal value-only residuals.

The Jacobian columns are measured differences. A prediction is not a residual.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
import argparse
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
sys.path.insert(0, str(ROOT/'scripts'))
import derive_nsc_ks_gate_step as Step
import derive_nsc_ks_gate_value as Value
from recursive_horizons.nsc_ks_profile_identity import profile_identity
from recursive_horizons.nsc_ks_value_evaluator import merit
from recursive_horizons.nsc_local_incoming_family import LocalAxialFunction, LocalIncomingFamily

PRINCIPAL = ROOT/'results/development/nsc-incoming-surface-principal.json'
PARENT_HISTORY = ROOT/'results/development/nsc-ks-gate-step-u0-trial.json'
PARENT_LABEL = 'step-u0-h5e-4'
ITERATE5 = ROOT/'results/development/nsc-ks-coupled-newton-n32-truncated-iterate5.json'
ITERATE6 = ROOT/'results/development/nsc-ks-coupled-newton-n32-truncated-iterate6.json'
REGISTER = ROOT/'results/development/nsc-ks-gate-lm-subspace.json'
# Time movement at this clock, from the recorded CF4 halving. Drops at or
# below ten times that movement are not accepted as search progress.
TIME_MOVEMENT = 5.53e-7
LOOP = {
    'grid_nodes': 256,
    'degree': 32,
    'rtol': 2e-11,
    'atol': 2e-13,
    'max_step': 0.0005,
    'phase_nodes': 192,
    'length': 1.6,
    'integrator': 'cf4',
}
# (function row, Chebyshev index, probe width, label)
# w1–w3 at 1e-5 left the linear range. The fine labels replace them.
PROBES = (
    (1, 0, 1e-4, 'lm-u0'), (1, 1, 1e-4, 'lm-u1'),
    (1, 2, 1e-4, 'lm-u2'), (1, 3, 1e-4, 'lm-u3'),
    (0, 0, 1e-5, 'lm-w0'),
    (0, 1, 1e-8, 'lm-w1-fine'), (0, 2, 1e-8, 'lm-w2-fine'),
    (0, 3, 1e-8, 'lm-w3-fine'),
)


def principal_weights():
    matrix = json.loads(PRINCIPAL.read_text())['matrix']['matrix_intervals']
    mid = np.array([[0.5 * (lo + hi) for lo, hi in row] for row in matrix], float)
    # Columns are U and w curvature. Weights are squared residual scales.
    return float(np.dot(mid[:, 1], mid[:, 1])), float(np.dot(mid[:, 0], mid[:, 0]))


def physical_metric(n, center, inner=0.03, outer=0.06, samples=2001):
    """H^3 Gram on w and H^1 Gram on U, scaled by the certified principal matrix."""
    w_weight, u_weight = principal_weights()
    left, right = center - outer, center + outer
    z = np.linspace(left, right, int(samples))
    dz = float(z[1] - z[0])
    eye = np.eye(n)
    values = np.empty((n, 4, len(z)))
    for index in range(n):
        function = LocalAxialFunction(tuple(eye[index]), center, inner, outer)
        for order in range(4):
            values[index, order] = [function(float(point), order) for point in z]
    gram = lambda orders: sum(values[:, order] @ values[:, order].T for order in orders) * dz
    metric = np.zeros((2 * n, 2 * n))
    metric[:n, :n] = w_weight * gram(range(4))
    metric[n:, n:] = u_weight * gram(range(2))
    if np.linalg.eigvalsh(metric).min() <= 0:
        raise ArithmeticError('physical metric is not positive definite')
    return metric


def trust_radius(metric):
    """Physical size of iterate6's accepted repair. Not a residual."""
    parent = np.asarray(json.loads(ITERATE5.read_text())['history']['coefficients'], float)
    child = np.asarray(json.loads(ITERATE6.read_text())['history']['coefficients'], float)
    step = (child - parent).ravel()
    if step.shape != (metric.shape[0],):
        raise ValueError('iterate6 repair is not an n=32 step')
    size = float(np.sqrt(step @ metric @ step))
    if not np.isfinite(size) or size <= 0:
        raise ValueError('iterate6 physical size is not positive')
    return size


def levenberg_step(residual, jacobian, metric, radius):
    """Minimise the linear model inside the physical trust region."""
    residual = np.asarray(residual, float).ravel()
    jacobian = np.asarray(jacobian, float)
    metric = np.asarray(metric, float)
    if jacobian.ndim != 2 or jacobian.shape[0] != residual.shape[0]:
        raise ValueError('Jacobian and residual shapes differ')
    if metric.shape != (jacobian.shape[1], jacobian.shape[1]):
        raise ValueError('metric must match the probed coordinates')
    normal = jacobian.T @ jacobian
    slope = jacobian.T @ residual

    def sized(mu):
        step = np.linalg.solve(normal + mu * metric, -slope)
        size = float(np.sqrt(step @ metric @ step))
        return step, size

    step, size = sized(0.0)
    if size > radius:
        lo, hi = 0.0, 1.0
        while sized(hi)[1] > radius:
            hi *= 4.0
        for _ in range(50):
            mid = 0.5 * (lo + hi)
            step, size = sized(mid)
            if size > radius:
                lo = mid
            else:
                hi = mid
        step, size = sized(hi)
    predicted = residual + jacobian @ step
    return step, predicted


def sampled_step(residual, jacobian, metric, radius, box):
    """Keep every coordinate inside the width where its column was measured."""
    step, predicted = levenberg_step(residual, jacobian, metric, radius)
    if np.all(np.abs(step) <= box):
        return step, predicted, radius
    lo, hi = 0.0, radius
    chosen = (np.zeros_like(step), residual, 0.0)
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        step, predicted = levenberg_step(residual, jacobian, metric, mid)
        if np.all(np.abs(step) <= box * (1 + 1e-8)):
            chosen = (step, predicted, mid)
            lo = mid
        else:
            hi = mid
    return chosen


def _coefficients(path):
    payload = json.loads(Path(path).read_text())
    return np.asarray(payload['history']['coefficients'], float)


def write_history(coefficients, name):
    if not Value.LABEL_RE.fullmatch(name):
        raise ValueError('history name must match [a-z0-9-]+')
    family = LocalIncomingFamily(coefficients)
    identity = profile_identity(family, include_normal_window=True)
    if identity in Step.FORBIDDEN:
        raise ValueError('forbidden history')
    if family.radius_lower_bound() <= 0:
        raise ValueError('history must keep a positive radius')
    path = (ROOT/'results/development'/f'nsc-ks-gate-history-{name}.json').resolve()
    if path.parent != (ROOT/'results/development').resolve():
        raise ValueError('history path escapes development')
    record = {
        'schema': 'NSC-KS-GATE-HISTORY-v1',
        'profile_identity': identity,
        'history': family.description(),
        'prediction_is_not_a_residual': True,
    }
    path.write_text(json.dumps(record, sort_keys=True, indent=2)+'\n')
    return path, identity


def probe_plan(parent):
    plan = []
    for row, index, width, name in PROBES:
        coefficients = parent.copy()
        coefficients[row, index] += width
        plan.append((name, coefficients, row, index, width))
    return plan


def measured_column(parent_gradient, label, width):
    gradient = Value.reconstruct_gradient(label)
    return (gradient - parent_gradient) / width


def subspace_metric(full, columns):
    index = [row * 32 + mode for row, mode in columns]
    return full[np.ix_(index, index)], index


def propose():
    parent = _coefficients(PARENT_HISTORY)
    parent_gradient = Value.reconstruct_gradient(PARENT_LABEL)
    columns = []
    jacobian_columns = []
    for name, _coefficients_probe, row, index, width in probe_plan(parent):
        label_path = Value.output_paths(name)[1]
        if not label_path.exists():
            raise FileNotFoundError('probe is not measured: ' + name)
        columns.append((row, index))
        jacobian_columns.append(measured_column(parent_gradient, name, width).ravel())
    jacobian = np.column_stack(jacobian_columns)
    full = physical_metric(32, LocalIncomingFamily(parent).center)
    metric, index = subspace_metric(full, columns)
    repair = trust_radius(full)
    # A coordinate may move by at most the width where its column was measured.
    box = np.zeros(len(columns))
    for slot, (_name, _coeff, _row, _mode, width) in enumerate(probe_plan(parent)):
        box[slot] = width
    sampled = float(np.sqrt(box @ metric @ box))
    radius = min(repair, sampled)
    step, predicted, radius = sampled_step(
        parent_gradient, jacobian, metric, radius, box)
    coefficients = parent.copy()
    for (row, mode), value in zip(columns, step):
        coefficients[row, mode] += value
    path, identity = write_history(coefficients, 'lm-candidate')
    proposal = Step.proposal(parent, coefficients - parent, 
        json.loads(PARENT_HISTORY.read_text())['profile_identity'],
        merit(parent_gradient))
    record = {
        'schema': 'NSC-KS-GATE-LM-SUBSPACE-v1',
        'accountable_author': 'Douglas Ek',
        'status': 'OPEN: linear model only; candidate is not a residual',
        'loop_numerics': LOOP,
        'time_movement_at_this_clock': TIME_MOVEMENT,
        'acceptance_requires_merit_drop_above': 10 * TIME_MOVEMENT,
        'parent_label': PARENT_LABEL,
        'parent_merit': merit(parent_gradient),
        'probed_coordinates': [{'row': row, 'mode': mode} for row, mode in columns],
        'iterate6_physical_size': repair,
        'sampled_box_radius': sampled,
        'trust_radius': radius,
        'step': step.tolist(),
        'predicted_merit': merit(predicted.reshape(parent_gradient.shape)),
        'prediction_is_not_a_residual': True,
        'candidate_history': str(path.relative_to(ROOT)),
        'candidate_profile_identity': identity,
        'proposal': proposal,
        'physical_EXISTENCE_certificate': False,
        'physical_NONEXISTENCE_certificate': False,
    }
    REGISTER.write_text(json.dumps(record, sort_keys=True, indent=2)+'\n')
    return record


def broyden_propose():
    """One rank-one update from the measured candidate. No new probes."""
    record = json.loads(REGISTER.read_text())
    parent = _coefficients(PARENT_HISTORY)
    parent_gradient = Value.reconstruct_gradient(PARENT_LABEL)
    child_gradient = Value.reconstruct_gradient('lm-candidate')
    columns = []
    jacobian_columns = []
    widths = []
    for name, _probe, row, index, width in probe_plan(parent):
        columns.append((row, index))
        widths.append(width)
        jacobian_columns.append(measured_column(parent_gradient, name, width).ravel())
    jacobian = np.column_stack(jacobian_columns)
    step = np.asarray(record['step'], float)
    change = (child_gradient - parent_gradient).ravel()
    jacobian = jacobian + np.outer(change - jacobian @ step, step) / np.dot(step, step)
    full = physical_metric(32, LocalIncomingFamily(parent).center)
    metric, _index = subspace_metric(full, columns)
    box = np.asarray(widths, float)
    radius = min(trust_radius(full), float(np.sqrt(box @ metric @ box)))
    next_step, predicted, radius = sampled_step(child_gradient, jacobian, metric, radius, box)
    coefficients = _coefficients(ROOT/record['candidate_history'])
    for (row, mode), value in zip(columns, next_step):
        coefficients[row, mode] += value
    path, identity = write_history(coefficients, 'lm-broyden')
    out = {
        'schema': 'NSC-KS-GATE-LM-BROYDEN-v1',
        'accountable_author': 'Douglas Ek',
        'status': 'OPEN: Broyden prediction is not a residual',
        'parent_profile_identity': record['candidate_profile_identity'],
        'parent_merit': merit(child_gradient),
        'predicted_merit': merit(predicted.reshape(child_gradient.shape)),
        'prediction_is_not_a_residual': True,
        'trust_radius': radius,
        'step': next_step.tolist(),
        'candidate_history': str(path.relative_to(ROOT)),
        'candidate_profile_identity': identity,
        'physical_EXISTENCE_certificate': False,
        'physical_NONEXISTENCE_certificate': False,
    }
    path_out = ROOT/'results/development/nsc-ks-gate-lm-broyden.json'
    path_out.write_text(json.dumps(out, sort_keys=True, indent=2)+'\n')
    return out


def next_probe():
    parent = _coefficients(PARENT_HISTORY)
    for name, coefficients, *_rest in probe_plan(parent):
        if not Value.output_paths(name)[1].exists():
            path, _identity = write_history(coefficients, name)
            return name, path
    return None, None


def run_next(workers):
    name, path = next_probe()
    if name is None:
        record = propose() if not REGISTER.exists() else json.loads(REGISTER.read_text())
        if not Value.output_paths('lm-candidate')[1].exists():
            return Value.run(
                47, workers, 14400.0, 'lm-candidate', LOOP,
                record['candidate_history'])
        return record
    return Value.run(47, workers, 14400.0, name, LOOP, str(path))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--run-next', action='store_true')
    mode.add_argument('--propose', action='store_true')
    mode.add_argument('--broyden', action='store_true')
    parser.add_argument('--workers', type=int, default=9)
    args = parser.parse_args()
    if args.broyden:
        result = broyden_propose()
        print(json.dumps({
            'predicted_merit': result['predicted_merit'],
            'parent_merit': result['parent_merit'],
            'candidate_profile_identity': result['candidate_profile_identity'],
        }, indent=2))
    elif args.propose:
        result = propose()
        print(json.dumps({
            'predicted_merit': result['predicted_merit'],
            'parent_merit': result['parent_merit'],
            'trust_radius': result['trust_radius'],
            'candidate_profile_identity': result['candidate_profile_identity'],
        }, indent=2))
    else:
        result = run_next(args.workers)
        print(json.dumps({
            'merit': result.get('merit'),
            'constraint_maxima': result.get('constraint_maxima'),
            'status': result.get('status'),
            'wall_seconds': result.get('wall_seconds'),
        }, indent=2))
