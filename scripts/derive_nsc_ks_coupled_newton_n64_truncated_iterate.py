#!/usr/bin/env python3
"""One truncated, damped, unit-clipped n=64 step from measured columns.

Iterate6 owns a 64-column Jacobian. The measured Chebyshev modes 32..63 are
appended in the blocked basis [w0:32, U0:32, w32:64, U32:64]. A later
measured n=64 history already stores the family order [w0:64, U0:64]. The
step is an equilibrated truncated SVD solution, clipped to the unit ball and
scaled. It is evolved only when that linear prediction beats the measured
residual on both components. The prediction is not the residual.
"""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
import argparse
from dataclasses import dataclass
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
import derive_nsc_ks_coupled_newton_n64_high_modes as N64
import derive_nsc_ks_coupled_retained_control as R
import recursive_horizons.nsc_ks_history_evaluator as H
from recursive_horizons.nsc_dirac_source_phase import formal_source_phase_coefficient
from recursive_horizons.nsc_evolved_incoming_constraints import (
    compatible_history_slots, surface_geometry_response)
from recursive_horizons.nsc_ks_cutoff_bridge_evaluator import (
    signed_family_angular_square_weights, source_edge_from_phase)
from recursive_horizons.nsc_ks_profile_identity import profile_identity
from recursive_horizons.nsc_ks_spacetime_variation import chart_coordinates
from recursive_horizons.nsc_local_incoming_family import (
    LocalIncomingFamily, plateau_derivatives)

EXISTENCE = 3e-11
MODES = 64
DIRECTIONS = 128
OWNERS = (
    'scripts/derive_nsc_ks_coupled_newton_n64_truncated_iterate.py',
    'docs/nsc-ks-coupled-newton-n64-truncated-iterate.md',
    'scripts/derive_nsc_ks_coupled_newton_n16.py',
    'scripts/derive_nsc_ks_coupled_newton_n64_high_modes.py',
    'src/recursive_horizons/nsc_ks_history_evaluator.py',
    'src/recursive_horizons/nsc_ks_cutoff_bridge_evaluator.py',
    'src/recursive_horizons/nsc_dirac_source_phase.py',
    'src/recursive_horizons/nsc_ks_energy_propagator.py',
    'src/recursive_horizons/nsc_ks_difference_envelope.py',
    'src/recursive_horizons/nsc_ks_batched_constraints.py',
    'src/recursive_horizons/nsc_local_incoming_family.py',
)


@dataclass(frozen=True)
class LocalIncomingFamily64(LocalIncomingFamily):
    """Same (w, U) class at 64 Chebyshev coefficients per function."""

    def __post_init__(self):
        raw = np.asarray(self.coefficients)
        if (raw.ndim != 2 or raw.shape != (2, MODES)
                or np.iscomplexobj(raw) or not np.isfinite(raw).all()):
            raise ValueError('two finite real vectors with 64 coefficients required')
        center = chart_coordinates(1.)[1]+.15 if self.center is None else float(self.center)
        plateau_derivatives(center, center, self.axial_inner, self.axial_outer)
        if (not np.isfinite([self.normal_inner, self.normal_outer]).all()
                or not 0 < self.normal_inner < self.normal_outer):
            raise ValueError('positive ordered normal support radii required')
        object.__setattr__(self, 'center', center)
        object.__setattr__(self, 'coefficients', tuple(tuple(map(float, row)) for row in raw))


def digest(path):
    path = Path(path)
    return sha256((path if path.is_absolute() else ROOT/path).read_bytes()).hexdigest()


def paths(name):
    directory = ROOT/'results/development/artifacts'/f'nsc-ks-coupled-newton-n64-truncated-{name}'
    output = ROOT/'results/development'/f'nsc-ks-coupled-newton-n64-truncated-{name}.json'
    proposal = ROOT/'results/development'/f'nsc-ks-coupled-newton-n64-truncated-{name}-proposal.json'
    return directory, output, proposal


def blocked_to_family(delta):
    """Map [w0:32, U0:32, w32:64, U32:64] into family order [w0:64, U0:64]."""
    delta = np.asarray(delta, float)
    if delta.shape != (DIRECTIONS,):
        raise ValueError('blocked step must have 128 amplitudes')
    w = np.concatenate([delta[0:32], delta[64:96]])
    u = np.concatenate([delta[32:64], delta[96:128]])
    return np.vstack([w, u])


def rank_truncated_delta(gradient, jacobian, kept_modes):
    residual = np.concatenate([gradient[:, 0], gradient[:, 1]])
    matrix = np.vstack([jacobian[:, :, 0].T, jacobian[:, :, 1].T])
    row = np.linalg.norm(matrix, axis=1)
    col = np.linalg.norm(matrix, axis=0)
    row = np.where(row > 0, row, 1.0)
    col = np.where(col > 0, col, 1.0)
    scaled = matrix / row[:, None] / col[None, :]
    left, singular, right = np.linalg.svd(scaled, full_matrices=False)
    count = int(kept_modes)
    if isinstance(kept_modes, bool) or count < 1 or count > singular.size:
        raise ValueError('kept mode count is outside the singular spectrum')
    coef = np.zeros_like(singular)
    coef[:count] = (left[:, :count].T @ (-residual / row)) / singular[:count]
    delta = (right.T @ coef) / col
    if not np.isfinite(delta).all():
        raise ValueError('truncated step is not finite')
    return delta, residual, matrix, singular


def predicted_maxima(residual, matrix, step):
    prediction = residual + matrix @ step
    n = prediction.size // 2
    return np.array([
        float(np.max(np.abs(prediction[:n]))),
        float(np.max(np.abs(prediction[n:]))),
    ])


def clip_and_scale(delta, scale):
    norm = float(np.linalg.norm(delta))
    if not np.isfinite(norm) or norm == 0.0:
        raise ValueError('truncated step must be a finite nonzero direction')
    clipped = delta if norm <= 1.0 else delta / norm
    step = float(scale) * clipped
    step_norm = float(np.linalg.norm(step))
    if not 0.0 < float(scale) <= 1.0:
        raise ValueError('damping scale must lie in (0, 1]')
    if step_norm > 1.0 + 1e-12:
        raise ValueError('damped clipped step must stay inside the unit ball')
    return clipped, step, norm, step_norm


def full_blocked_jacobian(ctx):
    tangent = np.zeros((N64.HIGH, len(ctx['target']), 2), float)
    family_records = {}
    seen = 0
    for key in sorted(ctx['archived']):
        column, record = N64.load_tangent(ctx, key)
        tangent += column
        family_records.update(record['family_records'])
        seen += 1
    if seen != len(ctx['archived']):
        raise ValueError('n=64 high-mode columns are incomplete')
    slots, slot_tangents = compatible_history_slots(ctx['metric'], ctx['target'], N64.FROZEN_PLUS)
    geometry = surface_geometry_response(slots, slot_tangents, ctx['coeff'])
    weights = signed_family_angular_square_weights(family_records)
    first = ctx['archive'].family_entries(next(iter(ctx['archived'])))[0][0]
    phase = formal_source_phase_coefficient(
        ctx['metric'], ctx['target'], angular=1., rho_up=first.rho_up, gauss_nodes=192)
    _edge, edge_tangent = source_edge_from_phase(phase, weights, ctx['coeff']['a'])
    extra = geometry['action_gradient_tangent'][1:] + tangent + edge_tangent[1:]
    base = np.asarray(ctx['parent']['history_jacobian'], float)
    if base.shape[0] != 64 or extra.shape != (N64.HIGH, len(ctx['target']), 2):
        raise ValueError('blocked Jacobian blocks do not have 64 columns each')
    return np.asarray(ctx['parent']['action_gradient'], float), np.concatenate([base, extra], axis=0)


def system_from_parent(parent_path):
    parent = json.loads(Path(parent_path).read_text())
    coefficients = np.asarray(parent['history']['coefficients'], float)
    jacobian = np.asarray(parent['history_jacobian'], float)
    gradient = np.asarray(parent['action_gradient'], float)
    if coefficients.shape == (2, 32) and jacobian.shape[0] == 64:
        published = json.loads(N64.OUTPUT.read_text())
        if parent['profile_identity'] != published['profile_identity']:
            raise ValueError('blocked n=64 columns belong only to the iterate6 history')
        ctx = N64.context()
        gradient, jacobian = full_blocked_jacobian(ctx)
        return {
            'layout': 'blocked-n32-plus-n64-high',
            'apply': 'blocked',
            'gradient': gradient,
            'jacobian': jacobian,
            'coefficients': coefficients,
            'parent': parent,
            'published_ranks': published['ranks'],
        }
    if coefficients.shape == (2, MODES) and jacobian.shape[0] == DIRECTIONS:
        return {
            'layout': 'family-n64',
            'apply': 'family',
            'gradient': gradient,
            'jacobian': jacobian,
            'coefficients': coefficients,
            'parent': parent,
            'published_ranks': None,
        }
    raise ValueError('parent must be iterate6 or a measured 64-coefficient history')


def require_published_rank(system, kept_modes, raw_norm, clipped_maxima):
    ranks = system['published_ranks']
    if ranks is None:
        return
    row = next(item for item in ranks if item['kept_modes'] == int(kept_modes))
    if float(raw_norm) != float(row['step_norm']):
        raise ValueError('truncated step norm does not match the published n=64 rank')
    if clipped_maxima.tolist() != row['unit_clipped_maxima']:
        raise ValueError('unit-clipped prediction does not match the published n=64 rank')


def candidate_from_step(system, step):
    coefficients = np.asarray(system['coefficients'], float)
    if system['apply'] == 'blocked':
        update = blocked_to_family(step)
        padded = np.zeros((2, MODES), float)
        padded[:, :32] = coefficients
        values = padded + update
    else:
        values = coefficients + np.asarray(step, float).reshape(2, MODES)
    family = LocalIncomingFamily64(values)
    identity = profile_identity(family, include_normal_window=True)
    parent_identity = system['parent']['profile_identity']
    if identity in D.FORBIDDEN_IDENTITIES or identity == parent_identity:
        raise ValueError('n=64 candidate is forbidden or unchanged')
    if family.radius_lower_bound() <= 0:
        raise ValueError('n=64 candidate must keep a positive radius bound')
    return family, identity


def build_proposal(parent_path, kept_modes, scale, name):
    system = system_from_parent(parent_path)
    gradient = system['gradient']
    jacobian = system['jacobian']
    measured = np.max(np.abs(gradient), axis=0)
    parent_measured = np.asarray(system['parent']['constraint_maxima'], float)
    if not np.array_equal(measured, parent_measured):
        raise ValueError('parent residual does not match its action gradient')
    delta, residual, matrix, singular = rank_truncated_delta(gradient, jacobian, kept_modes)
    _clipped, step, raw_norm, step_norm = clip_and_scale(delta, scale)
    clipped_maxima = predicted_maxima(residual, matrix, _clipped)
    if float(scale) == 1.0:
        require_published_rank(system, kept_modes, raw_norm, clipped_maxima)
    predicted = predicted_maxima(residual, matrix, step)
    if not D.beats_measured(predicted, measured):
        raise ValueError('damped clipped prediction must beat the measured residual on both components')
    family, identity = candidate_from_step(system, step)
    solve = family.collocation_nodes(N16.SOLVE_COUNT)
    verify = family.collocation_nodes(N16.VERIFY_COUNT)
    target = H.fixed_target_union(solve, verify)
    if target.tolist() != system['parent']['z']:
        raise ValueError('n=64 step must keep the measured node union')
    kept = int(kept_modes)
    return {
        'schema': 'NSC-KS-COUPLED-NEWTON-N64-TRUNCATED-ITERATE-PROPOSAL-v1',
        'accountable_author': 'Douglas Ek',
        'status': 'OPEN: truncated n=64 iterate authorized for one measured evolution',
        'name': name,
        'parent_path': str(Path(parent_path).relative_to(ROOT)) if Path(parent_path).is_absolute() else parent_path,
        'column_layout': system['layout'],
        'kept_modes': kept,
        'step_scale': float(scale),
        'raw_step_norm': raw_norm,
        'step_norm': step_norm,
        'retained_singular_ratio': float(singular[0] / singular[kept-1]),
        'parent_profile_identity': system['parent']['profile_identity'],
        'parent_constraint_maxima': measured.tolist(),
        'predicted_all_node_residual_maxima': predicted.tolist(),
        'prediction_reaches_existence': bool(np.all(predicted <= EXISTENCE)),
        'prediction_is_not_a_nonlinear_residual': True,
        'candidate_profile_identity': identity,
        'candidate_coefficients': family.description()['coefficients'],
        'radius_lower_bound': family.radius_lower_bound(),
        'evolve_authorized': True,
        'z': system['parent']['z'],
        'physical_EXISTENCE_certificate': False,
        'physical_NONEXISTENCE_certificate': False,
        'source_hashes': {p: digest(p) for p in OWNERS},
        'input_hashes': {
            str(Path(parent_path).resolve().relative_to(ROOT)): digest(parent_path),
            str(N64.OUTPUT.relative_to(ROOT)): digest(N64.OUTPUT),
        },
    }


def context(proposal):
    if not proposal.get('evolve_authorized'):
        raise ValueError('n=64 evolution requires an authorized prediction')
    base = R.context()
    family = LocalIncomingFamily64(np.asarray(proposal['candidate_coefficients'], float))
    identity = profile_identity(family, include_normal_window=True)
    if identity != proposal['candidate_profile_identity']:
        raise ValueError('n=64 candidate identity changed')
    if identity in D.FORBIDDEN_IDENTITIES:
        raise ValueError('forbidden history must not be evolved')
    solve = family.collocation_nodes(N16.SOLVE_COUNT)
    verify = family.collocation_nodes(N16.VERIFY_COUNT)
    target = H.fixed_target_union(solve, verify)
    if target.tolist() != proposal['z']:
        raise ValueError('n=64 iterate must keep the measured node union')
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
    N16.RETARDED = DIRECTIONS
    N16.DIRECTORY = directory
    N16.OUTPUT = output
    R.FAMILY_CAP = 240.0
    R.CPU_CAP = 14400.0


def run(parent_path, kept_modes, scale, name, cpu_budget, max_new):
    directory, output, proposal_path = paths(name)
    bind(directory, output)
    if proposal_path.exists():
        proposal = json.loads(proposal_path.read_text())
    else:
        proposal = build_proposal(parent_path, kept_modes, scale, name)
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
    reaches = bool(complete and np.all(measured <= EXISTENCE))
    result['schema'] = 'NSC-KS-COUPLED-NEWTON-N64-TRUNCATED-ITERATE-v1'
    result['status'] = (
        'OPEN: complete truncated n=64 iterate residual; uncertified physical gate'
        if complete else 'OPEN: partial truncated n=64 iterate residual')
    result['parent_profile_identity'] = proposal['parent_profile_identity']
    result['column_layout'] = proposal['column_layout']
    result['layout'] = (
        f"rectangular-n64-truncated-rank-{proposal['kept_modes']}-scale-{proposal['step_scale']:g}")
    result['kept_modes'] = proposal['kept_modes']
    result['step_scale'] = proposal['step_scale']
    result['step_norm'] = proposal['step_norm']
    result['linear_predicted_all_node_residual_maxima'] = predicted.tolist()
    result['prediction_minus_measured'] = None if not complete else (predicted - measured).tolist()
    result['residual_improved_on_all_components'] = improved if complete else None
    result['nodal_residual_reaches_existence_tolerance'] = reaches if complete else None
    result['enclosed_error'] = None
    result['prediction_is_not_a_nonlinear_residual'] = True
    result['physical_EXISTENCE_certificate'] = False
    result['physical_NONEXISTENCE_certificate'] = False
    result['reproducer'] = 'python3 scripts/derive_nsc_ks_coupled_newton_n64_truncated_iterate.py --check --name '+name
    if complete:
        output.write_text(json.dumps(R.jsonable(result), sort_keys=True, indent=2)+'\n')
        progress = directory/'progress.json'
        if progress.exists():
            progress.unlink()
    else:
        (directory/'progress.json').write_text(json.dumps(R.jsonable(result), sort_keys=True, indent=2)+'\n')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--run', action='store_true')
    mode.add_argument('--propose', action='store_true')
    mode.add_argument('--check', action='store_true')
    parser.add_argument('--parent')
    parser.add_argument('--kept-modes', type=int, default=80)
    parser.add_argument('--scale', type=float, default=1.0)
    parser.add_argument('--name', required=True)
    parser.add_argument('--cpu-budget', type=float, default=7200.0)
    parser.add_argument('--max-new', type=int, default=60)
    args = parser.parse_args()
    if args.propose:
        if not args.parent:
            raise SystemExit('propose requires --parent')
        proposal = build_proposal(args.parent, args.kept_modes, args.scale, args.name)
        _directory, _output, proposal_path = paths(args.name)
        if proposal_path.exists():
            raise SystemExit('proposal already exists')
        proposal_path.write_text(json.dumps(proposal, sort_keys=True, indent=2)+'\n')
        print(json.dumps({
            'status': proposal['status'],
            'candidate_profile_identity': proposal['candidate_profile_identity'],
            'step_norm': proposal['step_norm'],
            'raw_step_norm': proposal['raw_step_norm'],
            'predicted_all_node_residual_maxima': proposal['predicted_all_node_residual_maxima'],
            'parent_constraint_maxima': proposal['parent_constraint_maxima'],
            'radius_lower_bound': proposal['radius_lower_bound'],
            'evolve_authorized': proposal['evolve_authorized'],
        }, indent=2))
    elif args.run:
        if not args.parent and not paths(args.name)[2].exists():
            raise SystemExit('run requires --parent or an existing proposal')
        result = run(
            args.parent, args.kept_modes, args.scale, args.name, args.cpu_budget, args.max_new)
        print(json.dumps({
            'status': result['status'],
            'completed_positive_families': result['completed_positive_families'],
            'constraint_maxima': result.get('constraint_maxima'),
            'linear_predicted_all_node_residual_maxima': result.get('linear_predicted_all_node_residual_maxima'),
            'prediction_minus_measured': result.get('prediction_minus_measured'),
            'residual_improved_on_all_components': result.get('residual_improved_on_all_components'),
            'physical_EXISTENCE_certificate': result.get('physical_EXISTENCE_certificate'),
        }, indent=2))
    else:
        directory, output, proposal_path = paths(args.name)
        bind(directory, output)
        proposal = json.loads(proposal_path.read_text())
        result = N16.assemble(context(proposal), replay=True)
        saved = json.loads(output.read_text())
        if np.asarray(result['constraint_maxima'], float).tolist() != saved['constraint_maxima']:
            raise ValueError('n=64 truncated iterate replay residual changed')
        print(json.dumps({
            'status': saved['status'],
            'constraint_maxima': saved['constraint_maxima'],
            'residual_improved_on_all_components': saved['residual_improved_on_all_components'],
            'replay_dirac_evolution': False,
        }, indent=2))
