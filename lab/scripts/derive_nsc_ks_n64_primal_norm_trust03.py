#!/usr/bin/env python3
"""Shrink the measured ball direction to 0.3 of its N break-even norm.

The rejected step aa7e1982… fixed this direction's N quadratic coefficient.
At 0.3 times the norm where linear N gain meets that coefficient, the
predicted N gain is 10/3 of the quadratic term and beta still falls.
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
DIRECTION = ROOT/'results/development/nsc-ks-coupled-newton-n64-primal-norm-margin8-proposal.json'
MEASURED = ROOT/'results/development/nsc-ks-coupled-newton-n64-primal-norm-margin8.json'
DIRECTORY = ROOT/'results/development/artifacts/nsc-ks-coupled-newton-n64-primal-norm-trust03'
OUTPUT = ROOT/'results/development/nsc-ks-coupled-newton-n64-primal-norm-trust03.json'
PROPOSAL = ROOT/'results/development/nsc-ks-coupled-newton-n64-primal-norm-trust03-proposal.json'
OWNERS = (
    'scripts/derive_nsc_ks_n64_primal_norm_trust03.py',
    'docs/nsc-ks-n64-primal-norm-trust03.md',
    'src/recursive_horizons/nsc_ks_primal_step_control.py',
)
PARENT_IDENTITY = '6534e55b1f3dc3b968c2dfcdf56805113a5f4a354199e10a16704455b2316478'
REJECTED_IDENTITY = 'aa7e198206f59e794b1869b760d1c3216e65140df91ae57e5e33e317c9931a70'
FACTOR = 0.3
N_MARGIN = 3.0
EXISTENCE = 3e-11


def digest(path):
    path = Path(path)
    return sha256((path if path.is_absolute() else ROOT/path).read_bytes()).hexdigest()


def load_direction():
    parent = json.loads(PARENT.read_text())
    direction = json.loads(DIRECTION.read_text())
    measured = json.loads(MEASURED.read_text())
    if parent['profile_identity'] != PARENT_IDENTITY:
        raise ValueError('trust parent must be the accepted beta-primary history')
    if parent['solver'] != 'primal-block DOP853 RMS on (A, D); tangent directions excluded':
        raise ValueError('trust parent Jacobian must come from primal-block control')
    if not parent['newton_step_accepted']:
        raise ValueError('trust parent must be an accepted improvement on both components')
    if direction['candidate_profile_identity'] != REJECTED_IDENTITY:
        raise ValueError('trust direction must be the rejected ball step')
    if measured['profile_identity'] != REJECTED_IDENTITY or measured['newton_step_accepted']:
        raise ValueError('trust direction must stay a rejected measurement')
    if measured['solver'] != parent['solver']:
        raise ValueError('trust direction must use primal-block control')
    step0 = np.asarray(direction['step'], float)
    if step0.shape != (I64.DIRECTIONS,):
        raise ValueError('trust direction must be one 128-coefficient step')
    h0 = float(np.linalg.norm(step0))
    if not np.isfinite(h0) or h0 <= 0.0 or abs(h0 - float(measured['step_norm'])) > 0.0:
        raise ValueError('stored ball step does not match its measured norm')
    gradient = np.asarray(parent['action_gradient'], float)
    jacobian = np.asarray(parent['history_jacobian'], float)
    residual = np.concatenate([gradient[:, 0], gradient[:, 1]])
    matrix = np.vstack([jacobian[:, :, 0].T, jacobian[:, :, 1].T])
    full = residual + matrix @ step0
    nodes = gradient.shape[0]
    predicted_full = np.array([
        float(np.max(np.abs(full[:nodes]))),
        float(np.max(np.abs(full[nodes:]))),
    ])
    stored = np.asarray(measured['linear_predicted_all_node_residual_maxima'], float)
    if not np.array_equal(predicted_full, stored):
        raise ValueError('stored ball step does not reproduce its linear prediction')
    parent_max = np.asarray(parent['constraint_maxima'], float)
    n_gain = float(parent_max[0] - stored[0])
    n_error = abs(float(measured['prediction_minus_measured'][0]))
    if n_gain <= 0.0 or n_error <= n_gain:
        raise ValueError('ball step must have raised N above its linear gain')
    coefficient = n_error / h0**2
    break_even = (n_gain / h0) / coefficient
    return parent, step0, h0, residual, matrix, parent_max, coefficient, break_even


def proposal():
    parent, step0, h0, residual, matrix, parent_max, coefficient, break_even = load_direction()
    step = step0 * (FACTOR * break_even / h0)
    step_norm = float(np.linalg.norm(step))
    if not (0.0 < step_norm < break_even):
        raise ValueError('trust step must sit below the measured N crossing')
    predicted = I64.predicted_maxima(residual, matrix, step)
    gain = parent_max - predicted
    n_margin = float(gain[0] / (coefficient * step_norm**2))
    if gain[0] <= 0.0 or gain[1] <= 0.0:
        raise ValueError('trust prediction must lower both components')
    if n_margin < N_MARGIN:
        raise ValueError('trust N gain is inside the measured quadratic term')
    if gain[1] <= 1e-10:
        raise ValueError('trust beta drop is below 1e-10')
    coefficients = np.asarray(parent['history']['coefficients'], float)
    candidate = I64.LocalIncomingFamily64(coefficients + step.reshape(2, 64))
    identity = profile_identity(candidate, include_normal_window=True)
    if identity in D.FORBIDDEN_IDENTITIES or identity in (PARENT_IDENTITY, REJECTED_IDENTITY):
        raise ValueError('candidate identity is forbidden or already measured')
    if candidate.radius_lower_bound() <= 0:
        raise ValueError('candidate radius bound must stay positive')
    solve = candidate.collocation_nodes(N16.SOLVE_COUNT)
    verify = candidate.collocation_nodes(N16.VERIFY_COUNT)
    target = H.fixed_target_union(solve, verify)
    if target.tolist() != parent['z']:
        raise ValueError('trust step must stay on the measured nodes')
    return {
        'schema': 'NSC-KS-N64-PRIMAL-NORM-TRUST03-PROPOSAL-v1',
        'accountable_author': 'Douglas Ek',
        'status': 'OPEN: one trust-region shrink authorized for measurement',
        'shrink_factor_of_break_even': FACTOR,
        'break_even_step_norm': break_even,
        'quadratic_n_coefficient': coefficient,
        'quadratic_n_coefficient_source': 'rejected aa7e1982 N prediction error / step_norm^2',
        'n_margin_over_quadratic': n_margin,
        'step_norm': step_norm,
        'step': step.tolist(),
        'parent_profile_identity': parent['profile_identity'],
        'rejected_profile_identity': REJECTED_IDENTITY,
        'parent_constraint_maxima': parent_max.tolist(),
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
        'input_hashes': {
            str(PARENT.relative_to(ROOT)): digest(PARENT),
            str(DIRECTION.relative_to(ROOT)): digest(DIRECTION),
            str(MEASURED.relative_to(ROOT)): digest(MEASURED),
        },
    }


def context(record):
    if not record.get('evolve_authorized'):
        raise ValueError('trust evolution requires an authorized prediction')
    family = I64.LocalIncomingFamily64(np.asarray(record['candidate_coefficients'], float))
    identity = profile_identity(family, include_normal_window=True)
    if identity != record['candidate_profile_identity'] or identity in D.FORBIDDEN_IDENTITIES:
        raise ValueError('candidate identity changed or is forbidden')
    base = R.context()
    solve = family.collocation_nodes(N16.SOLVE_COUNT)
    verify = family.collocation_nodes(N16.VERIFY_COUNT)
    target = H.fixed_target_union(solve, verify)
    if target.tolist() != record['z']:
        raise ValueError('trust nodes changed')
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
        'trust03_proposal': record,
    }


def annotate(result, ctx):
    record = ctx['trust03_proposal']
    measured = np.asarray(result['constraint_maxima'], float)
    predicted = np.asarray(record['predicted_all_node_residual_maxima'], float)
    parent = np.asarray(record['parent_constraint_maxima'], float)
    complete = result['all_retained_nonzero_angular_families_included']
    gain = parent - predicted
    error = None if not complete else predicted - measured
    improved = bool(complete and np.all(measured < parent) and ctx['family'].radius_lower_bound() > 0)
    honest = bool(complete and np.all(np.abs(error) < gain))
    reaches = bool(complete and np.all(measured <= EXISTENCE))
    result['schema'] = 'NSC-KS-N64-PRIMAL-NORM-TRUST03-v1'
    result['status'] = (
        'OPEN: measured trust-region shrink; uncertified physical gate'
        if complete else 'OPEN: partial trust-region shrink')
    result['solver'] = 'primal-block DOP853 RMS on (A, D); tangent directions excluded'
    result['shrink_factor_of_break_even'] = FACTOR
    result['quadratic_n_coefficient'] = record['quadratic_n_coefficient']
    result['n_margin_over_quadratic'] = record['n_margin_over_quadratic']
    result['step_norm'] = record['step_norm']
    result['linear_predicted_all_node_residual_maxima'] = predicted.tolist()
    result['prediction_minus_measured'] = None if error is None else error.tolist()
    result['prediction_minus_measured_is_honest'] = honest if complete else None
    result['parent_constraint_maxima'] = parent.tolist()
    result['parent_profile_identity'] = record['parent_profile_identity']
    result['rejected_profile_identity'] = record['rejected_profile_identity']
    result['residual_improved_on_all_components'] = improved if complete else None
    result['prediction_is_not_a_nonlinear_residual'] = True
    result['coefficient_step_evolved'] = True
    result['newton_step_accepted'] = bool(improved and honest)
    result['physical_EXISTENCE_certificate'] = False
    result['physical_NONEXISTENCE_certificate'] = False
    result['nodal_residual_reaches_existence_tolerance'] = reaches if complete else None
    result['enclosed_error'] = None
    result['reproducer'] = 'python3 scripts/derive_nsc_ks_n64_primal_norm_trust03.py --check'
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
        raise ValueError('trust proposal identity changed')
    if fresh['predicted_all_node_residual_maxima'] != record['predicted_all_node_residual_maxima']:
        raise ValueError('trust prediction changed')
    if fresh['step_norm'] != record['step_norm']:
        raise ValueError('trust step norm changed')
    if fresh['n_margin_over_quadratic'] != record['n_margin_over_quadratic']:
        raise ValueError('trust N margin changed')
    P.evolve_ks_difference_envelope = evolve_primal_controlled
    N16.RETARDED = I64.DIRECTIONS
    N16.DIRECTORY = DIRECTORY
    N16.OUTPUT = OUTPUT
    ctx = context(record)
    result = annotate(N16.assemble(ctx, replay=True), ctx)
    saved = json.loads(OUTPUT.read_text())
    if saved['constraint_maxima'] != result['constraint_maxima']:
        raise ValueError('trust residual changed')
    if saved['prediction_minus_measured'] != result['prediction_minus_measured']:
        raise ValueError('trust prediction error changed')
    if saved['physical_EXISTENCE_certificate'] or saved['physical_NONEXISTENCE_certificate']:
        raise ValueError('trust measurement is not a gate certificate')
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
            'break_even_step_norm': record['break_even_step_norm'],
            'n_margin_over_quadratic': record['n_margin_over_quadratic'],
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
