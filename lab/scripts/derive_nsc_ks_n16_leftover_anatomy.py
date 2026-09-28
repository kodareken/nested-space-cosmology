#!/usr/bin/env python3
"""Split the n=16 leftover without evolving a new history."""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
sys.path.insert(0, str(ROOT/'scripts'))
import derive_nsc_ks_coupled_newton_n16_damped as D
import derive_nsc_ks_coupled_newton_n16_damped_iterate4 as I4
from recursive_horizons.nsc_evolved_incoming_constraints import (
    compatible_history_slots, surface_geometry_response)
from recursive_horizons.nsc_ks_profile_identity import profile_identity
from recursive_horizons.nsc_local_history_newton import (
    DeclaredJetBoundary, HistoryCollocation, LocalHistoryNewtonSettings,
    bind_history_evaluator_receipt, _boundary_system, _control_system,
    _linear_step)
from recursive_horizons.nsc_local_incoming_constraints import local_error_budget
from recursive_horizons.nsc_local_incoming_family import LocalIncomingFamily

OUTPUT = ROOT/'results/development/nsc-ks-n16-leftover-anatomy.json'
OWNERS = (
    'scripts/derive_nsc_ks_n16_leftover_anatomy.py',
    'docs/nsc-ks-n16-leftover-anatomy.md',
    'src/recursive_horizons/nsc_local_history_newton.py',
)
EXISTENCE_TOLERANCE = 3e-11
CUTOFFS = (1e4, 1e6, 1e8)
WALK = (
    ROOT/'results/development/nsc-ks-coupled-newton-n16-damped-trial.json',
    ROOT/'results/development/nsc-ks-coupled-newton-n16-damped-iterate.json',
    ROOT/'results/development/nsc-ks-coupled-newton-n16-damped-iterate2.json',
    ROOT/'results/development/nsc-ks-coupled-newton-n16-damped-iterate3.json',
    ROOT/'results/development/nsc-ks-coupled-newton-n16-damped-iterate4.json',
)


def digest(path):
    path = Path(path)
    return sha256((path if path.is_absolute() else ROOT/path).read_bytes()).hexdigest()


def _round(value, digits=12):
    array = np.asarray(value, float)
    if array.shape == ():
        return None if not np.isfinite(array) else round(float(array), digits)
    return [None if not np.isfinite(item) else round(float(item), digits) for item in array.ravel()]


def scaled_factors(jacobian, residual, preconditioner, settings):
    J = np.asarray(jacobian, float)
    r = np.asarray(residual, float)
    use_geom = bool(settings.geometric_preconditioner and preconditioner is not None)
    P = np.asarray(preconditioner, float) if use_geom else None
    if P is not None and P.shape != J.shape:
        P = None
    row = np.linalg.norm(J, axis=1)
    col = np.linalg.norm(J, axis=0)
    if P is not None:
        col_p = np.linalg.norm(P, axis=0)
        if np.all(col_p > 0) and np.all(np.isfinite(col_p)):
            col = np.where(col > 0, np.sqrt(col * col_p), col_p)
    row = np.where(row > 0, row, 1.0)
    col = np.where(col > 0, col, 1.0)
    Js = J / row[:, None] / col[None, :]
    rs = r / row
    u, singular, vt = np.linalg.svd(Js, full_matrices=False)
    return {
        'J': J, 'r': r, 'row': row, 'col': col, 'Js': Js, 'rs': rs,
        'u': u, 'singular': singular, 'vt': vt,
    }


def delta_from_keep(factors, keep):
    u, singular, vt, rs, col = (
        factors['u'], factors['singular'], factors['vt'], factors['rs'], factors['col'])
    if not np.any(keep):
        return None
    delta_scaled = vt[keep].T @ ((u[:, keep].T @ (-rs)) / singular[keep])
    delta = delta_scaled / col
    if not np.isfinite(delta).all():
        return None
    return delta


def leftover_blocks(residual, jacobian, delta, n_n, n_beta, has_jet):
    pred = residual if delta is None else residual + jacobian @ delta
    n_end = n_n
    beta_end = n_n + n_beta
    blocks = {
        'control_max': float(np.max(abs(pred))),
        'N_max': float(np.max(abs(pred[:n_end]))),
        'beta_max': float(np.max(abs(pred[n_end:beta_end]))),
        'step_norm': None if delta is None else float(np.linalg.norm(delta)),
    }
    if has_jet:
        blocks['jet_max'] = float(np.max(abs(pred[beta_end:])))
    return blocks


def orthogonal_max(factors):
    projector = factors['u'] @ factors['u'].T
    leftover = (np.eye(len(factors['rs'])) - projector) @ factors['rs']
    return float(np.max(abs(leftover * factors['row'])))


def truncated_rows(factors, n_n, n_beta, has_jet):
    smax = float(factors['singular'][0])
    rows = []
    for cutoff in CUTOFFS:
        keep = factors['singular'] > (smax / cutoff)
        delta = delta_from_keep(factors, keep)
        row = leftover_blocks(
            factors['r'], factors['J'], delta, n_n, n_beta, has_jet)
        row['condition_cutoff'] = cutoff
        row['kept_modes'] = int(np.count_nonzero(keep))
        rows.append(row)
    return rows


def receipt_on(family, z, gradient, jacobian, source_identity, geometry_tangent):
    raw = {
        'z': z, 'action_gradient': gradient, 'history_jacobian': jacobian,
        'source_identity': source_identity,
        'history_identity': profile_identity(family),
        'full_retarded_state_derivative': True,
        'local_reference_gradient_tangent': geometry_tangent,
        'physical_constraint_status': 'OPEN',
        'physical_EXISTENCE_certificate': False,
        'error_budget': local_error_budget(),
    }
    return bind_history_evaluator_receipt(raw, family, expected_nodes=z)


def all_node_system(family, z, gradient, jacobian, geometry_tangent, boundary):
    r_n = np.asarray(gradient, float)[:, 0]
    r_beta = np.asarray(gradient, float)[:, 1]
    j_n = np.asarray(jacobian, float)[:, :, 0].T
    j_beta = np.asarray(jacobian, float)[:, :, 1].T
    r_bc, j_bc, _ = _boundary_system(family, boundary)
    geom_n = np.asarray(geometry_tangent, float)[:, :, 0].T
    geom_beta = np.asarray(geometry_tangent, float)[:, :, 1].T
    residual = np.concatenate([r_n, r_beta, r_bc])
    matrix = np.vstack([j_n, j_beta, j_bc])
    preconditioner = np.vstack([geom_n, geom_beta, j_bc])
    return residual, matrix, preconditioner, len(r_n), len(r_beta)


def dropped_jet_system(system):
    n_n = len(system['n_residual'])
    n_beta = len(system['beta_residual'])
    residual = np.concatenate([system['n_residual'], system['beta_residual']])
    jacobian = system['jacobian'][:n_n + n_beta]
    preconditioner = None if system['preconditioner'] is None else system['preconditioner'][:n_n + n_beta]
    return residual, jacobian, preconditioner, n_n, n_beta


def classify(solve_leftover, all_node_leftover, jet_leftover, truncated, singular):
    tiny = 100 * EXISTENCE_TOLERANCE
    solve_small = solve_leftover <= tiny
    verify_large = all_node_leftover >= 1e-5
    jet_owns = jet_leftover >= 0.5 * max(solve_leftover, EXISTENCE_TOLERANCE)
    smax = float(singular[0])
    tiny_modes = int(np.count_nonzero(singular < smax / 1e6))
    truncated_win = any(
        row['control_max'] <= 0.1 * all_node_leftover and row['kept_modes'] < 32
        for row in truncated)
    if solve_small and verify_large:
        return 'all_node_collocation'
    if truncated_win or tiny_modes >= 2:
        return 'truncated_svd'
    if jet_owns:
        return 'freed_jet'
    if solve_leftover > tiny:
        return 'necessary_relation_candidate'
    return 'no_winning_linearization'


def compute():
    trial = json.loads(I4.OUTPUT.read_text())
    for path, expected in trial['source_hashes'].items():
        if digest(path) != expected:
            raise ValueError('leftover anatomy input changed: '+path)
    if trial['profile_identity'] != I4.EXPECTED_IDENTITY:
        raise ValueError('anatomy must use the fifth clipped history')
    ctx = I4.context()
    family = ctx['family']
    z = np.asarray(trial['z'], float)
    gradient = np.asarray(trial['action_gradient'], float)
    jacobian = np.asarray(trial['history_jacobian'], float)
    source_identity = json.loads(D.OUTPUT.read_text())['source_identity']
    slots, ds = compatible_history_slots(family.metric(), z, 32)
    geometry = surface_geometry_response(slots, ds, ctx['coeff'])
    boundary = DeclaredJetBoundary(family.center, 0., ctx['seed_slope'], 0.)
    settings = LocalHistoryNewtonSettings()
    collocation = HistoryCollocation(16, layout='rectangular')
    solve_nodes = collocation.solve_nodes(family)
    indices = np.searchsorted(z, solve_nodes)
    if not np.array_equal(z[indices], solve_nodes):
        raise ValueError('exact owned n=16 solve-node subset required')
    solve_receipt = receipt_on(
        family, solve_nodes, gradient[indices], jacobian[:, indices],
        source_identity, geometry['action_gradient_tangent'][:, indices])
    solve_system = _control_system(solve_receipt, collocation, boundary)
    delta, diagnostics = _linear_step(
        solve_system['jacobian'], solve_system['residual'],
        solve_system['preconditioner'], settings)
    if delta is None:
        raise ValueError('solve rectangular step must remain available for anatomy')
    solve_factors = scaled_factors(
        solve_system['jacobian'], solve_system['residual'],
        solve_system['preconditioner'], settings)
    n_n = len(solve_system['n_residual'])
    n_beta = len(solve_system['beta_residual'])
    unclipped = leftover_blocks(
        solve_system['residual'], solve_system['jacobian'], delta, n_n, n_beta, True)
    all_node_pred = gradient + np.tensordot(delta, jacobian, axes=(0, 0))
    all_node = {
        'N_max': float(np.max(abs(all_node_pred[:, 0]))),
        'beta_max': float(np.max(abs(all_node_pred[:, 1]))),
        'control_max': float(np.max(abs(all_node_pred))),
    }
    all_residual, all_jacobian, all_pre, all_n, all_beta = all_node_system(
        family, z, gradient, jacobian, geometry['action_gradient_tangent'], boundary)
    all_delta, all_diag = _linear_step(all_jacobian, all_residual, all_pre, settings)
    all_factors = scaled_factors(all_jacobian, all_residual, all_pre, settings)
    all_unclipped = leftover_blocks(
        all_residual, all_jacobian, all_delta, all_n, all_beta, True)
    drop_r, drop_j, drop_p, drop_n, drop_b = dropped_jet_system(solve_system)
    drop_delta, drop_diag = _linear_step(drop_j, drop_r, drop_p, settings)
    drop_unclipped = leftover_blocks(drop_r, drop_j, drop_delta, drop_n, drop_b, False)
    moved = DeclaredJetBoundary(float(family.interval[0]), 0., ctx['seed_slope'], 0.)
    moved_system = _control_system(solve_receipt, collocation, moved)
    moved_delta, moved_diag = _linear_step(
        moved_system['jacobian'], moved_system['residual'],
        moved_system['preconditioner'], settings)
    moved_unclipped = leftover_blocks(
        moved_system['residual'], moved_system['jacobian'], moved_delta,
        len(moved_system['n_residual']), len(moved_system['beta_residual']), True)
    leftover_dir = all_node_pred / np.linalg.norm(all_node_pred)
    walk = []
    for path in WALK:
        record = json.loads(path.read_text())
        residual = np.asarray(record['action_gradient'], float)
        if residual.shape != leftover_dir.shape:
            raise ValueError('walk residual shape changed: '+str(path))
        walk.append({
            'path': str(path.relative_to(ROOT)),
            'profile_identity': record['profile_identity'],
            'constraint_maxima': record['constraint_maxima'],
            'leftover_direction_dot_residual': float(np.sum(leftover_dir * residual)),
        })
    dots = np.array([row['leftover_direction_dot_residual'] for row in walk])
    branch = classify(
        unclipped['control_max'], all_node['control_max'], unclipped['jet_max'],
        truncated_rows(solve_factors, n_n, n_beta, True), solve_factors['singular'])
    return {
        'schema': 'NSC-KS-N16-LEFTOVER-ANATOMY-v1',
        'accountable_author': 'Douglas Ek',
        'status': 'OPEN: leftover anatomy recorded; linear prediction is not a residual',
        'profile_identity': trial['profile_identity'],
        'existence_tolerance': EXISTENCE_TOLERANCE,
        'solve_unclipped_leftover': {key: _round(value) for key, value in unclipped.items()},
        'all_node_unclipped_leftover_from_solve_step': {
            key: _round(value) for key, value in all_node.items()},
        'all_node_system_unclipped_leftover': {
            key: _round(value) for key, value in all_unclipped.items()},
        'dropped_jet_unclipped_leftover': {
            key: _round(value) for key, value in drop_unclipped.items()},
        'moved_jet_unclipped_leftover': {
            key: _round(value) for key, value in moved_unclipped.items()},
        'solve_orthogonal_leftover_max': _round(orthogonal_max(solve_factors)),
        'all_node_orthogonal_leftover_max': _round(orthogonal_max(all_factors)),
        'solve_singular_values': _round(solve_factors['singular']),
        'solve_condition_number': _round(diagnostics['condition_number']),
        'all_node_condition_number': _round(all_diag['condition_number']),
        'dropped_jet_condition_number': _round(drop_diag['condition_number']),
        'moved_jet_condition_number': _round(moved_diag['condition_number']),
        'truncated_svd_solve': [
            {key: _round(value) if key != 'condition_cutoff' and key != 'kept_modes' else value
             for key, value in row.items()}
            for row in truncated_rows(solve_factors, n_n, n_beta, True)],
        'walk_leftover_direction_projections': walk,
        'walk_projection_min': _round(np.min(dots)),
        'walk_projection_max': _round(np.max(dots)),
        'walk_projection_sign_stable': bool(np.all(dots > 0) or np.all(dots < 0)),
        'classified_branch': branch,
        'scope': {
            'physical_local_gate': 'OPEN',
            'prediction_is_not_a_nonlinear_residual': True,
            'new_field_or_source_evolutions': 0,
            'physical_EXISTENCE_certificate': False,
            'physical_NONEXISTENCE_claimed': False,
            'metric_timestep': False,
        },
        'source_hashes': {p: digest(p) for p in OWNERS},
        'input_hashes': {str(I4.OUTPUT.relative_to(ROOT)): digest(I4.OUTPUT)},
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--record', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = compute()
    if args.record:
        if OUTPUT.exists():
            raise FileExistsError('leftover anatomy register exists; use --check')
        OUTPUT.write_text(json.dumps(result, sort_keys=True, indent=2, allow_nan=False)+'\n')
    elif result != json.loads(OUTPUT.read_text()):
        raise ValueError('leftover anatomy replay differs')
    print(json.dumps({
        'status': result['status'],
        'classified_branch': result['classified_branch'],
        'solve_unclipped_leftover': result['solve_unclipped_leftover'],
        'all_node_unclipped_leftover_from_solve_step': result[
            'all_node_unclipped_leftover_from_solve_step'],
        'all_node_system_unclipped_leftover': result['all_node_system_unclipped_leftover'],
        'dropped_jet_unclipped_leftover': result['dropped_jet_unclipped_leftover'],
        'moved_jet_unclipped_leftover': result['moved_jet_unclipped_leftover'],
        'solve_orthogonal_leftover_max': result['solve_orthogonal_leftover_max'],
        'solve_condition_number': result['solve_condition_number'],
        'truncated_svd_solve': result['truncated_svd_solve'],
        'walk_projection_sign_stable': result['walk_projection_sign_stable'],
        'walk_projection_min': result['walk_projection_min'],
        'walk_projection_max': result['walk_projection_max'],
    }, indent=2))
