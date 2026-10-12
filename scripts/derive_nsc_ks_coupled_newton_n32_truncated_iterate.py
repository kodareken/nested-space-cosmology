#!/usr/bin/env python3
"""One more truncated-SVD n=32 step from a measured assembled Jacobian."""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
import argparse
from hashlib import sha256
import json
from pathlib import Path
import re
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
sys.path.insert(0, str(ROOT/'scripts'))
import derive_nsc_ks_coupled_newton_n16 as N16
import derive_nsc_ks_coupled_newton_n16_damped as D
import derive_nsc_ks_coupled_retained_control as R
import recursive_horizons.nsc_ks_history_evaluator as H
from recursive_horizons.nsc_ks_profile_identity import profile_identity
from recursive_horizons.nsc_local_incoming_family import LocalIncomingFamily

NAME_RE = re.compile(r'^[a-z0-9-]+$')
OWNERS = (
    'scripts/derive_nsc_ks_coupled_newton_n32_truncated_iterate.py',
    'docs/nsc-ks-coupled-newton-n32-truncated-iterate.md',
    'scripts/derive_nsc_ks_coupled_newton_n16.py',
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


def truncated_delta(gradient, jacobian, cutoff):
    residual = np.concatenate([gradient[:, 0], gradient[:, 1]])
    matrix = np.vstack([jacobian[:, :, 0].T, jacobian[:, :, 1].T])
    row = np.linalg.norm(matrix, axis=1)
    col = np.linalg.norm(matrix, axis=0)
    row = np.where(row > 0, row, 1.0)
    col = np.where(col > 0, col, 1.0)
    scaled = matrix / row[:, None] / col[None, :]
    left, singular, right = np.linalg.svd(scaled, full_matrices=False)
    keep = singular > (singular[0] / float(cutoff))
    if not np.any(keep):
        raise ValueError('truncated SVD kept no modes')
    delta = (right[keep].T @ ((left[:, keep].T @ (-residual / row)) / singular[keep])) / col
    if not np.isfinite(delta).all():
        raise ValueError('truncated step is not finite')
    return delta, int(np.count_nonzero(keep)), float(singular[0] / singular[-1])


def require_name(name):
    # The name is interpolated into artifact paths. Only the regex match is
    # used, so the path cannot carry a separator or an extra suffix.
    if not isinstance(name, str):
        raise ValueError('name must match [a-z0-9-]+')
    match = NAME_RE.fullmatch(name)
    if match is None:
        raise ValueError('name must match [a-z0-9-]+')
    return match.group(0)


def paths(name):
    safe = require_name(name)
    artifacts = (ROOT/'results/development/artifacts').resolve()
    development = (ROOT/'results/development').resolve()
    directory = (artifacts/f'nsc-ks-coupled-newton-n32-truncated-{safe}').resolve()
    output = (development/f'nsc-ks-coupled-newton-n32-truncated-{safe}.json').resolve()
    proposal = (development/f'nsc-ks-coupled-newton-n32-truncated-{safe}-proposal.json').resolve()
    if directory.parent != artifacts or output.parent != development or proposal.parent != development:
        raise ValueError('name escapes the development directories')
    return directory, output, proposal


def relative_input(path):
    resolved = Path(path).resolve()
    root = ROOT.resolve()
    if root not in resolved.parents and resolved != root:
        raise ValueError('proposal input must stay inside the lab')
    return str(resolved.relative_to(root))


def build_proposal(parent_path, cutoff, scale, campaign_best_path=None):
    parent = json.loads(Path(parent_path).read_text())
    gradient = np.asarray(parent['action_gradient'], float)
    jacobian = np.asarray(parent['history_jacobian'], float)
    immediate = np.asarray(parent['constraint_maxima'], float)
    acceptance = immediate
    best_rel = None
    if campaign_best_path is not None:
        best = json.loads(Path(campaign_best_path).read_text())
        acceptance = np.asarray(best['constraint_maxima'], float)
        best_rel = relative_input(campaign_best_path)
        if best['profile_identity'] == parent['profile_identity']:
            raise ValueError('campaign best must be a different measured history')
    coefficients = np.asarray(parent['history']['coefficients'], float)
    delta, kept, condition = truncated_delta(gradient, jacobian, cutoff)
    step = float(scale) * delta
    step_norm = float(np.linalg.norm(step))
    if step_norm > 1.0 + 1e-12:
        raise ValueError('scaled truncated step must stay inside the unit ball')
    candidate = LocalIncomingFamily(coefficients + step.reshape(2, 32))
    identity = profile_identity(candidate, include_normal_window=True)
    if identity in D.FORBIDDEN_IDENTITIES or identity == parent['profile_identity']:
        raise ValueError('iterate candidate is forbidden or unchanged')
    if candidate.radius_lower_bound() <= 0:
        raise ValueError('iterate candidate must keep a positive radius bound')
    predicted = D.predicted_maxima(gradient, jacobian, step)
    if not D.beats_measured(predicted, acceptance):
        raise ValueError('truncated prediction must beat the acceptance residual on both components')
    parent_rel = relative_input(parent_path)
    proposal = {
        'schema': 'NSC-KS-COUPLED-NEWTON-N32-TRUNCATED-ITERATE-PROPOSAL-v1',
        'accountable_author': 'Douglas Ek',
        'status': 'OPEN: truncated n=32 iterate authorized for one measured evolution',
        'parent_path': parent_rel,
        'condition_cutoff': float(cutoff),
        'step_scale': float(scale),
        'kept_modes': kept,
        'parent_condition_number': condition,
        'parent_profile_identity': parent['profile_identity'],
        'immediate_parent_constraint_maxima': immediate.tolist(),
        'parent_constraint_maxima': acceptance.tolist(),
        'step_norm': step_norm,
        'predicted_all_node_residual_maxima': predicted.tolist(),
        'prediction_is_not_a_nonlinear_residual': True,
        'candidate_profile_identity': identity,
        'candidate_coefficients': candidate.description()['coefficients'],
        'radius_lower_bound': candidate.radius_lower_bound(),
        'evolve_authorized': True,
        'z': parent['z'],
        'source_hashes': {p: digest(p) for p in OWNERS},
        'input_hashes': {parent_rel: digest(parent_path)},
    }
    if best_rel is not None:
        proposal['campaign_best_path'] = best_rel
        proposal['campaign_best_constraint_maxima'] = acceptance.tolist()
        proposal['status'] = (
            'OPEN: truncated step from the immediate parent, accepted only if it '
            'beats the campaign best on both components')
        proposal['input_hashes'][best_rel] = digest(campaign_best_path)
    return proposal


def context(proposal):
    base = R.context()
    family = LocalIncomingFamily(np.asarray(proposal['candidate_coefficients'], float))
    identity = profile_identity(family, include_normal_window=True)
    if identity != proposal['candidate_profile_identity']:
        raise ValueError('iterate candidate identity changed')
    solve = family.collocation_nodes(N16.SOLVE_COUNT)
    verify = family.collocation_nodes(N16.VERIFY_COUNT)
    target = H.fixed_target_union(solve, verify)
    if target.tolist() != proposal['z']:
        raise ValueError('iterate must keep the measured node union')
    return {
        **base,
        'family': family,
        'identity': identity,
        'solve': solve,
        'verify': verify,
        'target': target,
        'hashes': {
            **base['archive'].input_hashes,
            **{path: digest(path) for path in OWNERS},
            **proposal['input_hashes'],
        },
        'proposal': {
            'constraint_maxima_with_verified_phase_values': proposal['parent_constraint_maxima'],
        },
        'selected': {'step_scale': proposal['step_scale'], 'step_norm': proposal['step_norm']},
    }


def bind(directory, output):
    N16.RETARDED = 64
    N16.DIRECTORY = directory
    N16.OUTPUT = output
    R.FAMILY_CAP = 120.0
    R.CPU_CAP = 7200.0


def run(parent_path, cutoff, scale, name, cpu_budget, max_new):
    directory, output, proposal_path = paths(name)
    bind(directory, output)
    if proposal_path.exists():
        proposal = json.loads(proposal_path.read_text())
    else:
        proposal = build_proposal(parent_path, cutoff, scale)
        proposal_path.write_text(json.dumps(proposal, sort_keys=True, indent=2)+'\n')
    ctx = context(proposal)
    directory.mkdir(parents=True, exist_ok=True)
    start = time.process_time()
    new = 0
    for key in sorted(ctx['archived']):
        path, _payload = N16.family_paths(key)
        if path.exists():
            continue
        if new >= max_new:
            continue
        remaining = cpu_budget - (time.process_time()-start)
        if remaining <= 1.0:
            break
        record = N16.new_family(ctx, key, cpu_limit=remaining-1.0)
        new += record['new_operator_solves']
        print(json.dumps({
            'positive_family': record['positive_family'],
            'CPU_seconds': record['CPU_seconds'],
            'matter_change_maxima': record['matter_change_maxima'],
        }), flush=True)
    result = N16.assemble(ctx)
    measured = np.asarray(result['constraint_maxima'], float)
    parent_residual = np.asarray(proposal['parent_constraint_maxima'], float)
    predicted = np.asarray(proposal['predicted_all_node_residual_maxima'], float)
    complete = result['all_retained_nonzero_angular_families_included']
    improved = bool(complete and np.all(measured < parent_residual) and ctx['family'].radius_lower_bound() > 0)
    result['schema'] = 'NSC-KS-COUPLED-NEWTON-N32-TRUNCATED-ITERATE-v1'
    result['status'] = (
        'OPEN: complete truncated n=32 iterate residual; uncertified physical gate'
        if complete else 'OPEN: partial truncated n=32 iterate residual')
    result['parent_profile_identity'] = proposal['parent_profile_identity']
    result['layout'] = f'rectangular-n32-truncated-svd-{cutoff:g}-scale-{scale:g}'
    result['condition_cutoff'] = float(cutoff)
    result['step_scale'] = float(scale)
    result['kept_modes'] = proposal['kept_modes']
    result['linear_predicted_all_node_residual_maxima'] = predicted.tolist()
    result['prediction_minus_measured'] = None if not complete else (predicted - measured).tolist()
    result['residual_improved_on_all_components'] = improved if complete else None
    result['prediction_is_not_a_nonlinear_residual'] = True
    result['physical_EXISTENCE_certificate'] = False
    result['physical_NONEXISTENCE_certificate'] = False
    result['reproducer'] = 'python3 scripts/derive_nsc_ks_coupled_newton_n32_truncated_iterate.py --check'
    if complete:
        output.write_text(json.dumps(R.jsonable(result), sort_keys=True, indent=2)+'\n')
    else:
        (directory/'progress.json').write_text(json.dumps(R.jsonable(result), sort_keys=True, indent=2)+'\n')
    return result


def propose(parent_path, cutoff, scale, name, campaign_best_path=None):
    _directory, _output, proposal_path = paths(name)
    proposal = build_proposal(parent_path, cutoff, scale, campaign_best_path)
    if proposal_path.exists():
        old = json.loads(proposal_path.read_text())
        if old['candidate_profile_identity'] != proposal['candidate_profile_identity']:
            raise ValueError('refusing to replace a proposal with a different candidate')
        if not np.array_equal(
                np.asarray(old['candidate_coefficients'], float),
                np.asarray(proposal['candidate_coefficients'], float)):
            raise ValueError('refusing to replace a proposal with different coefficients')
    proposal_path.write_text(json.dumps(proposal, sort_keys=True, indent=2)+'\n')
    return proposal


def historical_context(name):
    """Replay against the hash set stored with the families.

    Owner-byte repairs change this script's digest. The saved operators stay
    the scientific record, so replay uses the hash set they share.
    """
    directory, output, proposal_path = paths(name)
    bind(directory, output)
    proposal = json.loads(proposal_path.read_text())
    ctx = context(proposal)
    stored = None
    for key in sorted(ctx['archived']):
        path, payload = N16.family_paths(key)
        if not path.exists():
            raise ValueError('truncated trail family is missing: '+path.name)
        record = json.loads(path.read_text())
        if digest(payload) != record['payload']['sha256']:
            raise ValueError('truncated trail payload changed: '+path.name)
        if stored is None:
            stored = record['source_hashes']
        elif record['source_hashes'] != stored:
            raise ValueError('truncated trail family hash sets disagree')
    ctx['hashes'] = stored
    return ctx, output


def replay_name(name, family_key=None):
    ctx, output = historical_context(name)
    saved = json.loads(output.read_text())
    key = family_key if family_key is not None else sorted(ctx['archived'])[0]
    _record, _change, _tangent, _local = N16.read_family(ctx, key, replay=True)
    result = N16.assemble(ctx, replay=False)
    if np.asarray(result['constraint_maxima'], float).tolist() != saved['constraint_maxima']:
        raise ValueError('truncated iterate replay residual changed')
    result['replayed_family'] = list(key)
    result['replay_dirac_evolution'] = False
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--run', action='store_true')
    mode.add_argument('--propose', action='store_true')
    mode.add_argument('--check', action='store_true')
    parser.add_argument('--parent')
    parser.add_argument('--cutoff', type=float)
    parser.add_argument('--scale', type=float, default=1.0)
    parser.add_argument('--campaign-best')
    parser.add_argument('--name', required=True)
    parser.add_argument('--cpu-budget', type=float, default=7200.0)
    parser.add_argument('--max-new', type=int, default=60)
    args = parser.parse_args()
    require_name(args.name)
    if args.run or args.propose:
        if args.parent is None or args.cutoff is None:
            parser.error('--parent and --cutoff are required to propose or run')
    if args.run:
        result = run(args.parent, args.cutoff, args.scale, args.name, args.cpu_budget, args.max_new)
    elif args.propose:
        result = propose(args.parent, args.cutoff, args.scale, args.name, args.campaign_best)
    else:
        result = replay_name(args.name)
    print(json.dumps({
        'status': result['status'],
        'candidate_profile_identity': result.get('candidate_profile_identity'),
        'completed_positive_families': result.get('completed_positive_families'),
        'constraint_maxima': result.get('constraint_maxima'),
        'linear_predicted_all_node_residual_maxima': result.get(
            'linear_predicted_all_node_residual_maxima',
            result.get('predicted_all_node_residual_maxima')),
        'prediction_minus_measured': result.get('prediction_minus_measured'),
        'residual_improved_on_all_components': result.get('residual_improved_on_all_components'),
        'replayed_family': result.get('replayed_family'),
    }, indent=2))
