#!/usr/bin/env python3
"""Propose a clipped/damped n=16 Newton step from the saved residual.

The unclipped n=16 rectangular step has condition ~7e8 and step norm ~19.85.
This owner reuses that recorded residual/Jacobian, clips to the existing
max_step_norm, then walks the existing line-search scales. It does not evolve
families and does not treat a linear prediction as a residual.
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

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
sys.path.insert(0, str(ROOT/'scripts'))
import derive_nsc_ks_coupled_newton_n16 as N16
import derive_nsc_ks_coupled_newton_n16_next as NEXT
from recursive_horizons.nsc_evolved_incoming_constraints import (
    compatible_history_slots, surface_geometry_response)
from recursive_horizons.nsc_local_history_newton import (
    DeclaredJetBoundary, HistoryCollocation, LocalHistoryNewtonSettings,
    bind_history_evaluator_receipt, _clip_step, _control_system, _linear_step)
from recursive_horizons.nsc_local_incoming_family import LocalIncomingFamily
from recursive_horizons.nsc_local_incoming_constraints import local_error_budget
from recursive_horizons.nsc_ks_profile_identity import profile_identity

OUTPUT = ROOT/'results/development/nsc-ks-coupled-newton-n16-damped.json'
FORBIDDEN_IDENTITIES = (
    '2d2588c3c9b1c83fe1a3f291314b8030300c4f14d3b12c24eadd2664eb20231b',
    '949a188181bd9ba0df43bc8ea00db2c89ce35a458017052b7af928059f819680',
)
OWNERS = (
    'scripts/derive_nsc_ks_coupled_newton_n16_damped.py',
    'docs/nsc-ks-coupled-newton-n16-damped.md',
    'src/recursive_horizons/nsc_local_history_newton.py',
    'src/recursive_horizons/nsc_evolved_incoming_constraints.py',
    'src/recursive_horizons/nsc_local_incoming_family.py',
)


def digest(path):
    path = Path(path)
    return sha256((path if path.is_absolute() else ROOT/path).read_bytes()).hexdigest()


def line_search_scales(settings):
    scales = []
    scale = 1.0
    while scale >= settings.min_line_search_scale - 1e-15:
        scales.append(float(scale))
        scale *= settings.line_search_factor
    return tuple(scales)


def predicted_maxima(gradient, jacobian, step):
    prediction = gradient + np.tensordot(step, jacobian, axes=(0, 0))
    return np.max(abs(prediction), axis=0)


def beats_measured(predicted, measured):
    return bool(np.all(np.asarray(predicted, float) < np.asarray(measured, float)))


def damped_scale_rows(family, delta, gradient, jacobian, measured, settings):
    clipped, clip_norm = _clip_step(delta, settings.max_step_norm)
    rows = []
    selected = None
    for scale in line_search_scales(settings):
        step = scale * clipped
        candidate = LocalIncomingFamily(
            np.asarray(family.coefficients, float) + step.reshape(2, 16))
        radius = candidate.radius_lower_bound()
        identity = profile_identity(candidate)
        predicted = predicted_maxima(gradient, jacobian, step)
        forbidden = identity in FORBIDDEN_IDENTITIES
        improves = radius > 0 and not forbidden and beats_measured(predicted, measured)
        row = {
            'step_scale': scale,
            'step_norm': float(np.linalg.norm(step)),
            'clipped_direction_norm': float(clip_norm),
            'radius_lower_bound': radius,
            'candidate_profile_identity': identity,
            'predicted_all_node_residual_maxima': predicted.tolist(),
            'prediction_beats_measured': improves,
            'forbidden_identity': forbidden,
            'positive_radius': radius > 0,
            'prediction_is_not_a_nonlinear_residual': True,
        }
        rows.append(row)
        if selected is None and improves:
            selected = {
                **row,
                'direction': clipped.reshape(2, 16).tolist(),
                'candidate': candidate.description(),
            }
    return clipped, clip_norm, rows, selected


def inputs():
    proposal = json.loads(NEXT.OUTPUT.read_text())
    for path, expected in {**proposal['source_hashes'], **proposal['input_hashes']}.items():
        if digest(path) != expected:
            raise ValueError('n=16 Newton proposal input changed: '+path)
    ctx = N16.context()
    family = ctx['family']
    if profile_identity(family) != proposal['current_profile_identity']:
        raise ValueError('damped step must start from the measured n=16 history')
    if family.radius_lower_bound() <= 0:
        raise ValueError('n=16 family must keep a positive radius bound')
    if ctx['target'].tolist() != proposal['z']:
        raise ValueError('damped step must keep the measured n=16 node union')
    gradient = np.asarray(proposal['action_gradient_with_192_node_phase'], float)
    jacobian = np.asarray(proposal['full_retarded_history_jacobian'], float)
    measured = np.asarray(proposal['constraint_maxima_with_192_node_phase'], float)
    if not np.isfinite(gradient).all() or not np.isfinite(jacobian).all():
        raise ValueError('saved n=16 residual/Jacobian must be finite')
    if jacobian.shape[0] != 32:
        raise ValueError('n=16 damped step requires all 32 retarded directions')
    source_identity = proposal['source_identity']
    if source_identity.get('source_preparation_changed'):
        raise ValueError('damped step must keep the same upstream source')
    return ctx, proposal, source_identity, gradient, jacobian, measured


def propose_layout(ctx, layout, gradient, jacobian, source_identity, measured, settings):
    family, z = ctx['family'], ctx['target']
    slots, ds = compatible_history_slots(family.metric(), z, 32)
    geometry = surface_geometry_response(slots, ds, ctx['coeff'])
    boundary = DeclaredJetBoundary(family.center, 0., ctx['seed_slope'], 0.)
    collocation = HistoryCollocation(16, layout=layout)
    solve_nodes = collocation.solve_nodes(family)
    indices = np.searchsorted(z, solve_nodes)
    if not np.array_equal(z[indices], solve_nodes):
        raise ValueError('exact owned n=16 solve-node subset required')
    raw = {
        'z': solve_nodes, 'action_gradient': gradient[indices],
        'history_jacobian': jacobian[:, indices], 'source_identity': source_identity,
        'history_identity': profile_identity(family), 'full_retarded_state_derivative': True,
        'local_reference_gradient_tangent': geometry['action_gradient_tangent'][:, indices],
        'physical_constraint_status': 'OPEN', 'physical_EXISTENCE_certificate': False,
        'error_budget': local_error_budget(),
    }
    receipt = bind_history_evaluator_receipt(raw, family, expected_nodes=solve_nodes)
    system = _control_system(receipt, collocation, boundary)
    delta, diagnostics = _linear_step(
        system['jacobian'], system['residual'], system['preconditioner'], settings)
    row = {
        'layout': layout,
        'rank': diagnostics['rank'],
        'unknowns': diagnostics['n_unknowns'],
        'equations': diagnostics['n_equations'],
        'condition_number': diagnostics['condition_number'],
        'rank_or_conditioning': diagnostics['rank_or_conditioning'],
        'initial_control_residual_max': float(np.max(abs(system['residual']))),
        'full_minus_geometric_jacobian_max': system['full_minus_geometric_max'],
        'full_retarded_jacobian_consumed': True,
        'trial_accepted': False,
        'max_step_norm': settings.max_step_norm,
        'unclipped_step_norm': None,
        'clipped_step_norm': None,
        'scales': [],
        'selected_scale': None,
    }
    if delta is None:
        row['proposal_status'] = 'OPEN: rank or conditioning rejected the linear step'
        return row, None
    unclipped_norm = float(np.linalg.norm(delta))
    unclipped = LocalIncomingFamily(
        np.asarray(family.coefficients, float) + delta.reshape(2, 16))
    clipped, clip_norm, scales, selected = damped_scale_rows(
        family, delta, gradient, jacobian, measured, settings)
    row.update({
        'unclipped_step_norm': unclipped_norm,
        'clipped_step_norm': float(clip_norm),
        'unclipped_candidate_profile_identity': profile_identity(unclipped),
        'unclipped_full_scale_rejected': profile_identity(unclipped) in FORBIDDEN_IDENTITIES,
        'scales': scales,
        'selected_scale': None if selected is None else selected['step_scale'],
    })
    if selected is None:
        row['proposal_status'] = (
            'OPEN: no clipped n=16 scale predicts a residual below the measured maxima')
        return row, None
    control_pred = float(np.max(abs(
        system['residual'] + selected['step_scale'] * system['jacobian'] @ clipped)))
    row.update({
        'proposal_status': 'OPEN: fresh nonlinear source evolution required',
        'step_scale': selected['step_scale'],
        'step_norm': selected['step_norm'],
        'direction': selected['direction'],
        'candidate': selected['candidate'],
        'candidate_profile_identity': selected['candidate_profile_identity'],
        'predicted_control_residual_max': control_pred,
        'predicted_all_node_residual_maxima': selected['predicted_all_node_residual_maxima'],
        'prediction_is_not_a_nonlinear_residual': True,
        'radius_lower_bound': selected['radius_lower_bound'],
    })
    return row, selected


def compute():
    ctx, proposal, source_identity, gradient, jacobian, measured = inputs()
    family, z = ctx['family'], ctx['target']
    settings = LocalHistoryNewtonSettings()
    proposals = []
    selected_rect = None
    for layout in ('rectangular', 'tau'):
        row, selected = propose_layout(
            ctx, layout, gradient, jacobian, source_identity, measured, settings)
        proposals.append(row)
        if layout == 'rectangular':
            selected_rect = selected
    evolve_authorized = selected_rect is not None
    if evolve_authorized:
        status = (
            'OPEN: clipped n=16 rectangular prediction improves; '
            'fresh evolution required before acceptance')
        stall = None
        selected_name = 'rectangular'
    else:
        status = (
            'OPEN: n=16 damped linear prediction does not improve; '
            'named search stall, no trial evolved')
        stall = 'n16_damped_linear_prediction_does_not_improve'
        selected_name = None
    return {
        'schema': 'NSC-KS-COUPLED-NEWTON-N16-DAMPED-v1',
        'accountable_author': 'Douglas Ek',
        'status': status,
        'current_history': family.description(),
        'current_profile_identity': profile_identity(family),
        'parent_n8_profile_identity': N16.PARENT_IDENTITY,
        'parent_n16_zero_pad_identity': proposal['current_profile_identity'],
        'source_identity': source_identity,
        'z': z.tolist(),
        'action_gradient_with_192_node_phase': gradient.tolist(),
        'full_retarded_history_jacobian': jacobian.tolist(),
        'constraint_maxima_with_192_node_phase': measured.tolist(),
        'constraint_maxima_with_verified_phase_values': measured.tolist(),
        'phase_gauss_nodes': 192,
        'nodal_phase_value_error_upper_N_beta': None,
        'boundary_control': DeclaredJetBoundary(
            family.center, 0., ctx['seed_slope'], 0.).description(),
        'forbidden_identities': list(FORBIDDEN_IDENTITIES),
        'proposals': proposals,
        'selected_for_next_nonlinear_trial': selected_name,
        'evolve_authorized': evolve_authorized,
        'named_search_stall': stall,
        'scope': {
            'physical_local_gate': 'OPEN',
            'missing_total_error_bound': None,
            'between_node_bound': None,
            'source_cutoff_phase_quadrature_error_bound': None,
            'old_phase_bounds_transferred': False,
            'new_field_or_source_evolutions': 0,
            'metric_timestep': False,
            'used_existing_clip_and_line_search': True,
            'full_scale_n16_not_evolved': True,
            'n8_next_not_evolved': True,
            'prediction_is_not_a_nonlinear_residual': True,
            'physical_NONEXISTENCE_claimed': False,
        },
        'source_hashes': {p: digest(p) for p in OWNERS},
        'input_hashes': {str(p.relative_to(ROOT)): digest(p) for p in (NEXT.OUTPUT,)},
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
            raise FileExistsError('n=16 damped Newton proposal exists; use --check')
        OUTPUT.write_text(json.dumps(result, sort_keys=True, indent=2, allow_nan=False)+'\n')
    elif result != json.loads(OUTPUT.read_text()):
        raise ValueError('n=16 damped Newton proposal replay differs')
    print(json.dumps({
        'status': result['status'],
        'current_maxima': result['constraint_maxima_with_192_node_phase'],
        'evolve_authorized': result['evolve_authorized'],
        'named_search_stall': result['named_search_stall'],
        'selected_for_next_nonlinear_trial': result['selected_for_next_nonlinear_trial'],
        'proposals': [{k: row[k] for k in row if k not in ('direction', 'candidate', 'scales')}
                      for row in result['proposals']],
    }, indent=2))
