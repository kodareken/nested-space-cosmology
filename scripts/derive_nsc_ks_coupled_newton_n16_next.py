#!/usr/bin/env python3
"""Propose the first n=16 rectangular Newton step from the zero-pad residual.

This is a linear prediction. It is not a measured residual at the new g.
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
import derive_nsc_ks_coupled_retained_control as R
from recursive_horizons.nsc_evolved_incoming_constraints import (
    compatible_history_slots, surface_geometry_response)
from recursive_horizons.nsc_local_history_newton import (
    DeclaredJetBoundary, HistoryCollocation, LocalHistoryNewtonSettings,
    bind_history_evaluator_receipt, _control_system, _linear_step)
from recursive_horizons.nsc_local_incoming_family import LocalIncomingFamily
from recursive_horizons.nsc_local_incoming_constraints import local_error_budget
from recursive_horizons.nsc_ks_profile_identity import profile_identity

OUTPUT = ROOT/'results/development/nsc-ks-coupled-newton-n16-next.json'
OWNERS = (
    'scripts/derive_nsc_ks_coupled_newton_n16_next.py',
    'docs/nsc-ks-coupled-newton-n16-next.md',
    'src/recursive_horizons/nsc_local_history_newton.py',
    'src/recursive_horizons/nsc_evolved_incoming_constraints.py',
    'src/recursive_horizons/nsc_local_incoming_family.py',
)


def digest(path):
    path = Path(path)
    return sha256((path if path.is_absolute() else ROOT/path).read_bytes()).hexdigest()


def inputs():
    trial = json.loads(N16.OUTPUT.read_text())
    ctx = N16.context()
    if R.jsonable(N16.assemble(ctx, replay=False)) != trial:
        raise ValueError('n=16 residual assembly receipts or aggregate changed')
    if trial['completed_positive_families'] != trial['required_positive_families']:
        raise ValueError('incomplete n=16 residual cannot propose a step')
    if trial['new_operator_solves'] != 60 or trial['reused_operator_solves'] != 0:
        raise ValueError('n=16 residual must evolve all sixty families')
    if trial['retarded_directions'] != 32 or trial['source_cutoff_phase_gauss_nodes'] != 192:
        raise ValueError('n=16 residual must keep 32 retarded directions and 192-node phase')
    gradient = np.asarray(trial['action_gradient'], float)
    jacobian = np.asarray(trial['history_jacobian'], float)
    if not np.isfinite(gradient).all() or not np.isfinite(jacobian).all():
        raise ValueError('n=16 residual/Jacobian must be finite')
    if trial['radius_lower_bound'] <= 0 or ctx['family'].radius_lower_bound() <= 0:
        raise ValueError('n=16 family must keep a positive radius bound')
    if jacobian.shape[0] != 32:
        raise ValueError('n=16 Newton step requires all 32 retarded directions')
    source_ids = []
    for item in trial['family_records']:
        rec = json.loads((ROOT/item['path']).read_text())
        with np.load(ROOT/rec['payload']['path'], allow_pickle=False) as f:
            meta = json.loads(f['metadata_json'].tobytes())
        source_ids.extend(meta['source_identity'])
    source_ids.sort(key=lambda value: value[0])
    if len({value[0] for value in source_ids}) != len(source_ids):
        raise ValueError('unique original upstream batch identities required')
    source_identity = {
        'fixed_upstream_batch_sha256': sha256(json.dumps(
            source_ids, sort_keys=True, separators=(',', ':')).encode()).hexdigest(),
        'batch_count': len(source_ids),
        'source_preparation_changed': False,
    }
    return ctx, trial, source_identity, gradient, jacobian


def compute():
    ctx, trial, source_identity, gradient, jacobian = inputs()
    family, z = ctx['family'], ctx['target']
    slots, ds = compatible_history_slots(family.metric(), z, 32)
    geometry = surface_geometry_response(slots, ds, ctx['coeff'])
    boundary = DeclaredJetBoundary(family.center, 0., ctx['seed_slope'], 0.)
    settings = LocalHistoryNewtonSettings()
    proposals = []
    for layout in ('rectangular', 'tau'):
        collocation = HistoryCollocation(16, layout=layout)
        solve_nodes = collocation.solve_nodes(family)
        indices = np.searchsorted(z, solve_nodes)
        if not np.array_equal(z[indices], solve_nodes):
            raise ValueError('exact owned n=16 solve-node subset required')
        raw = {'z': solve_nodes, 'action_gradient': gradient[indices],
            'history_jacobian': jacobian[:, indices], 'source_identity': source_identity,
            'history_identity': profile_identity(family), 'full_retarded_state_derivative': True,
            'local_reference_gradient_tangent': geometry['action_gradient_tangent'][:, indices],
            'physical_constraint_status': 'OPEN', 'physical_EXISTENCE_certificate': False,
            'error_budget': local_error_budget()}
        receipt = bind_history_evaluator_receipt(raw, family, expected_nodes=solve_nodes)
        system = _control_system(receipt, collocation, boundary)
        delta, diagnostics = _linear_step(system['jacobian'], system['residual'],
                                          system['preconditioner'], settings)
        row = {'layout': layout, 'rank': diagnostics['rank'], 'unknowns': diagnostics['n_unknowns'],
            'equations': diagnostics['n_equations'], 'condition_number': diagnostics['condition_number'],
            'rank_or_conditioning': diagnostics['rank_or_conditioning'],
            'initial_control_residual_max': float(np.max(abs(system['residual']))),
            'full_minus_geometric_jacobian_max': system['full_minus_geometric_max'],
            'full_retarded_jacobian_consumed': True, 'trial_accepted': False}
        if delta is not None:
            scale = 1.
            while scale >= 1/64:
                candidate = LocalIncomingFamily(np.asarray(family.coefficients)+scale*delta.reshape(2, 16))
                if candidate.radius_lower_bound() > 0:
                    break
                scale /= 2
            if scale < 1/64:
                row['proposal_status'] = 'OPEN: no positive-radius proposed scale'
            else:
                prediction = gradient+np.tensordot(scale*delta, jacobian, axes=(0, 0))
                row.update({'proposal_status': 'OPEN: fresh nonlinear source evolution required',
                    'step_scale': scale, 'step_norm': float(np.linalg.norm(scale*delta)),
                    'direction': delta.reshape(2, 16).tolist(), 'candidate': candidate.description(),
                    'candidate_profile_identity': profile_identity(candidate),
                    'predicted_control_residual_max': float(np.max(abs(system['residual']+scale*system['jacobian']@delta))),
                    'predicted_all_node_residual_maxima': np.max(abs(prediction), axis=0).tolist(),
                    'prediction_is_not_a_nonlinear_residual': True})
        proposals.append(row)
    measured = np.max(abs(gradient), axis=0)
    return {
        'schema': 'NSC-KS-COUPLED-NEWTON-N16-NEXT-v1',
        'accountable_author': 'Douglas Ek',
        'status': 'OPEN: n=16 rectangular Newton proposal; no trial accepted or physical root claimed',
        'current_history': family.description(),
        'current_profile_identity': profile_identity(family),
        'parent_n8_profile_identity': N16.PARENT_IDENTITY,
        'source_identity': source_identity,
        'z': z.tolist(),
        'action_gradient_with_192_node_phase': gradient.tolist(),
        'full_retarded_history_jacobian': jacobian.tolist(),
        'constraint_maxima_with_192_node_phase': measured.tolist(),
        'constraint_maxima_with_verified_phase_values': measured.tolist(),
        'phase_gauss_nodes': 192,
        'nodal_phase_value_error_upper_N_beta': None,
        'boundary_control': boundary.description(),
        'proposals': proposals,
        'selected_for_next_nonlinear_trial': 'rectangular',
        'selected_for_first_nonlinear_trial': 'rectangular',
        'scope': {
            'physical_local_gate': 'OPEN',
            'missing_total_error_bound': None,
            'between_node_bound': None,
            'source_cutoff_phase_quadrature_error_bound': None,
            'old_phase_bounds_transferred': False,
            'new_field_or_source_evolutions': 0,
            'metric_timestep': False,
            'resolution_raised_because_n8_stalled': True,
            'physical_NONEXISTENCE_claimed': False,
        },
        'source_hashes': {p: digest(p) for p in OWNERS},
        'input_hashes': {str(p.relative_to(ROOT)): digest(p) for p in (N16.OUTPUT,)},
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
            raise FileExistsError('n=16 Newton proposal exists; use --check')
        OUTPUT.write_text(json.dumps(result, sort_keys=True, indent=2, allow_nan=False)+'\n')
    elif result != json.loads(OUTPUT.read_text()):
        raise ValueError('n=16 Newton proposal replay differs')
    print(json.dumps({
        'status': result['status'],
        'current_maxima': result['constraint_maxima_with_192_node_phase'],
        'proposals': [{k: row[k] for k in row if k not in ('direction', 'candidate')}
                      for row in result['proposals']],
    }, indent=2))
