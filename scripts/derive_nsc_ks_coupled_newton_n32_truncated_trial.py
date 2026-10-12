#!/usr/bin/env python3
"""Evolve one truncated-SVD n=32 step from the assembled iterate4 Jacobian."""
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
import derive_nsc_ks_coupled_newton_n32_assembled as N32
import derive_nsc_ks_coupled_retained_control as R
import recursive_horizons.nsc_ks_history_evaluator as H
from recursive_horizons.nsc_ks_profile_identity import profile_identity
from recursive_horizons.nsc_local_incoming_family import LocalIncomingFamily

DIRECTORY = ROOT/'results/development/artifacts/nsc-ks-coupled-newton-n32-truncated-trial'
OUTPUT = ROOT/'results/development/nsc-ks-coupled-newton-n32-truncated-trial.json'
PROPOSAL = ROOT/'results/development/nsc-ks-coupled-newton-n32-truncated-proposal.json'
CUTOFF = 1e6
PARENT_RESIDUAL = np.array([0.004107193193038788, 0.005790788097514318])
OWNERS = (
    'scripts/derive_nsc_ks_coupled_newton_n32_truncated_trial.py',
    'docs/nsc-ks-coupled-newton-n32-truncated-trial.md',
    'scripts/derive_nsc_ks_coupled_newton_n16.py',
    'scripts/derive_nsc_ks_coupled_newton_n32_assembled.py',
    'src/recursive_horizons/nsc_ks_history_evaluator.py',
    'src/recursive_horizons/nsc_ks_cutoff_bridge_evaluator.py',
    'src/recursive_horizons/nsc_dirac_source_phase.py',
    'src/recursive_horizons/nsc_ks_energy_propagator.py',
    'src/recursive_horizons/nsc_ks_difference_envelope.py',
    'src/recursive_horizons/nsc_ks_batched_constraints.py',
    'src/recursive_horizons/nsc_local_incoming_family.py',
)


def digest(path):
    path = Path(path)
    return sha256((path if path.is_absolute() else ROOT/path).read_bytes()).hexdigest()


def truncated_step(gradient, jacobian, cutoff):
    residual = np.concatenate([gradient[:, 0], gradient[:, 1]])
    matrix = np.vstack([jacobian[:, :, 0].T, jacobian[:, :, 1].T])
    row = np.linalg.norm(matrix, axis=1)
    col = np.linalg.norm(matrix, axis=0)
    row = np.where(row > 0, row, 1.0)
    col = np.where(col > 0, col, 1.0)
    scaled = matrix / row[:, None] / col[None, :]
    left, singular, right = np.linalg.svd(scaled, full_matrices=False)
    if singular[0] <= 0:
        raise ValueError('assembled Jacobian has no positive singular value')
    keep = singular > (singular[0] / cutoff)
    if not np.any(keep):
        raise ValueError('truncated SVD kept no modes')
    delta_scaled = right[keep].T @ ((left[:, keep].T @ (-residual / row)) / singular[keep])
    delta = delta_scaled / col
    if not np.isfinite(delta).all():
        raise ValueError('truncated step is not finite')
    return delta, int(np.count_nonzero(keep)), float(singular[0] / singular[-1])


def proposal_from_parent():
    parent = json.loads(N32.OUTPUT.read_text())
    for path, expected in parent['source_hashes'].items():
        if digest(path) != expected:
            raise ValueError('n=32 assembled input changed: '+path)
    if parent['profile_identity'] != N32.EXPECTED_PADDED:
        raise ValueError('truncated trial must start from the zero-padded iterate4 history')
    if not parent['all_retained_nonzero_angular_families_included']:
        raise ValueError('incomplete n=32 assembly cannot seed a trial')
    gradient = np.asarray(parent['action_gradient'], float)
    jacobian = np.asarray(parent['history_jacobian'], float)
    measured = np.asarray(parent['constraint_maxima'], float)
    if not np.allclose(measured, PARENT_RESIDUAL):
        raise ValueError('parent residual changed')
    coefficients = np.asarray(parent['history']['coefficients'], float)
    delta, kept, condition = truncated_step(gradient, jacobian, CUTOFF)
    step_norm = float(np.linalg.norm(delta))
    if step_norm > 1.0 + 1e-12:
        raise ValueError('this cutoff must stay inside the unit step')
    candidate = LocalIncomingFamily(coefficients + delta.reshape(2, 32))
    identity = profile_identity(candidate, include_normal_window=True)
    if identity in D.FORBIDDEN_IDENTITIES or identity in (N32.EXPECTED_PADDED, N32.EXPECTED_PARENT):
        raise ValueError('truncated candidate is a forbidden or unchanged identity')
    if candidate.radius_lower_bound() <= 0:
        raise ValueError('truncated candidate must keep a positive radius bound')
    predicted = D.predicted_maxima(gradient, jacobian, delta)
    if not D.beats_measured(predicted, measured):
        raise ValueError('truncated prediction must beat the parent residual on both components')
    return {
        'schema': 'NSC-KS-COUPLED-NEWTON-N32-TRUNCATED-PROPOSAL-v1',
        'accountable_author': 'Douglas Ek',
        'status': 'OPEN: truncated n=32 step authorized for one measured evolution',
        'condition_cutoff': CUTOFF,
        'kept_modes': kept,
        'parent_condition_number': condition,
        'parent_profile_identity': parent['profile_identity'],
        'parent_constraint_maxima': measured.tolist(),
        'step_norm': step_norm,
        'predicted_all_node_residual_maxima': predicted.tolist(),
        'prediction_is_not_a_nonlinear_residual': True,
        'candidate_profile_identity': identity,
        'candidate_coefficients': candidate.description()['coefficients'],
        'radius_lower_bound': candidate.radius_lower_bound(),
        'evolve_authorized': True,
        'z': parent['z'],
        'source_hashes': {p: digest(p) for p in OWNERS},
        'input_hashes': {str(N32.OUTPUT.relative_to(ROOT)): digest(N32.OUTPUT)},
    }


def context(proposal):
    base = R.context()
    family = LocalIncomingFamily(np.asarray(proposal['candidate_coefficients'], float))
    identity = profile_identity(family, include_normal_window=True)
    if identity != proposal['candidate_profile_identity']:
        raise ValueError('truncated candidate identity changed')
    solve = family.collocation_nodes(N16.SOLVE_COUNT)
    verify = family.collocation_nodes(N16.VERIFY_COUNT)
    target = H.fixed_target_union(solve, verify)
    if target.tolist() != proposal['z']:
        raise ValueError('truncated trial must keep the iterate4 node union')
    hashes = {
        **base['archive'].input_hashes,
        **{path: digest(path) for path in OWNERS},
        str(PROPOSAL.relative_to(ROOT)): digest(PROPOSAL),
    }
    return {
        **base,
        'family': family,
        'identity': identity,
        'solve': solve,
        'verify': verify,
        'target': target,
        'hashes': hashes,
        'proposal': {
            'constraint_maxima_with_verified_phase_values': proposal['parent_constraint_maxima'],
        },
        'selected': {
            'step_scale': 1.0,
            'step_norm': proposal['step_norm'],
        },
    }


def bind():
    N16.RETARDED = 64
    N16.DIRECTORY = DIRECTORY
    N16.OUTPUT = OUTPUT
    R.FAMILY_CAP = 120.0
    R.CPU_CAP = 7200.0


def run(cpu_budget, max_new):
    bind()
    if PROPOSAL.exists():
        proposal = json.loads(PROPOSAL.read_text())
    else:
        proposal = proposal_from_parent()
        PROPOSAL.write_text(json.dumps(proposal, sort_keys=True, indent=2)+'\n')
    ctx = context(proposal)
    DIRECTORY.mkdir(parents=True, exist_ok=True)
    start = time.process_time()
    new = 0
    spent = 0.0
    for key in sorted(ctx['archived']):
        path, _payload = N16.family_paths(key)
        if path.exists():
            record, _, _, _ = N16.read_family(ctx, key)
            spent += record['CPU_seconds']
            continue
        if new >= max_new:
            continue
        remaining = cpu_budget - (time.process_time()-start)
        if remaining <= 1.0:
            break
        record = N16.new_family(ctx, key, cpu_limit=remaining-1.0)
        new += record['new_operator_solves']
        spent += record['CPU_seconds']
        print(json.dumps({
            'positive_family': record['positive_family'],
            'CPU_seconds': record['CPU_seconds'],
            'matter_change_maxima': record['matter_change_maxima'],
        }), flush=True)
        N16.assemble(ctx)
    result = N16.assemble(ctx)
    measured = np.asarray(result['constraint_maxima'], float)
    predicted = np.asarray(proposal['predicted_all_node_residual_maxima'], float)
    complete = result['all_retained_nonzero_angular_families_included']
    improved = bool(complete and np.all(measured < PARENT_RESIDUAL) and ctx['family'].radius_lower_bound() > 0)
    result['schema'] = 'NSC-KS-COUPLED-NEWTON-N32-TRUNCATED-TRIAL-v1'
    result['status'] = (
        'OPEN: complete truncated n=32 measured residual; uncertified physical gate'
        if complete else 'OPEN: partial truncated n=32 measured residual')
    result['parent_profile_identity'] = proposal['parent_profile_identity']
    result['layout'] = 'rectangular-n32-truncated-svd-1e6'
    result['condition_cutoff'] = CUTOFF
    result['kept_modes'] = proposal['kept_modes']
    result['linear_predicted_all_node_residual_maxima'] = predicted.tolist()
    result['prediction_minus_measured'] = None if not complete else (predicted - measured).tolist()
    result['residual_improved_on_all_components'] = improved if complete else None
    result['prediction_is_not_a_nonlinear_residual'] = True
    result['physical_EXISTENCE_certificate'] = False
    result['physical_NONEXISTENCE_certificate'] = False
    result['named_search_stall'] = None
    result['reproducer'] = 'python3 scripts/derive_nsc_ks_coupled_newton_n32_truncated_trial.py --check'
    result['runtime'] = {
        'CPU_seconds': time.process_time()-start,
        'new_operator_solves': result['new_operator_solves'],
    }
    if complete:
        OUTPUT.write_text(json.dumps(R.jsonable(result), sort_keys=True, indent=2)+'\n')
    else:
        (DIRECTORY/'progress.json').write_text(json.dumps(R.jsonable(result), sort_keys=True, indent=2)+'\n')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--run', action='store_true')
    mode.add_argument('--check', action='store_true')
    parser.add_argument('--cpu-budget', type=float, default=7200.0)
    parser.add_argument('--max-new', type=int, default=60)
    args = parser.parse_args()
    if args.run:
        result = run(args.cpu_budget, args.max_new)
    else:
        bind()
        proposal = json.loads(PROPOSAL.read_text())
        ctx = context(proposal)
        result = N16.assemble(ctx, replay=True)
        saved = json.loads(OUTPUT.read_text())
        if not result['all_retained_nonzero_angular_families_included']:
            raise ValueError('truncated trial replay is incomplete')
        if np.asarray(result['constraint_maxima'], float).tolist() != saved['constraint_maxima']:
            raise ValueError('truncated trial replay residual changed')
    print(json.dumps({
        'status': result['status'],
        'completed_positive_families': result['completed_positive_families'],
        'required_positive_families': result['required_positive_families'],
        'constraint_maxima': result.get('constraint_maxima'),
        'linear_predicted_all_node_residual_maxima': result.get('linear_predicted_all_node_residual_maxima'),
        'prediction_minus_measured': result.get('prediction_minus_measured'),
        'residual_improved_on_all_components': result.get('residual_improved_on_all_components'),
        'physical_EXISTENCE_certificate': result.get('physical_EXISTENCE_certificate', False),
    }, indent=2))
