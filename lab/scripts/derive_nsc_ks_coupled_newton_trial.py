#!/usr/bin/env python3
"""Restartable full-source evaluation of the selected rectangular Newton candidate.

Reuses the frozen upstream archive and existing evaluators. Operators from a
different g, including the committed group14 pilot, are not reused. The
linear Newton prediction is not this residual. No physical root is claimed.
"""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
import argparse
from hashlib import sha256
import json
from pathlib import Path
import signal
import sys
import time
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
sys.path.insert(0, str(ROOT/'scripts'))
import recursive_horizons.nsc_ks_history_evaluator as H
import derive_nsc_ks_coupled_retained_control as R
from recursive_horizons.nsc_dirac_source_phase import formal_source_phase_coefficient
from recursive_horizons.nsc_evolved_incoming_constraints import (
    compatible_history_slots, surface_geometry_response)
from recursive_horizons.nsc_ks_cutoff_bridge_evaluator import (
    signed_family_angular_square_weights, source_edge_from_phase)
from recursive_horizons.nsc_ks_profile_identity import profile_identity
from recursive_horizons.nsc_local_incoming_family import LocalIncomingFamily
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes

PROPOSAL = ROOT/'results/development/nsc-ks-coupled-newton-proposal.json'
DIRECTORY = ROOT/'results/development/artifacts/nsc-ks-coupled-newton-trial'
OUTPUT = ROOT/'results/development/nsc-ks-coupled-newton-trial.json'
EXPECTED_IDENTITY = '8baa9ec9642d234e43e2a9716824a73a16fe454032a9a928cfc3bc85bd3b681b'
PHASE_NODES = 192
MAX_NEW = 60
OWNERS = (
    'scripts/derive_nsc_ks_coupled_newton_trial.py',
    'docs/nsc-ks-coupled-newton-trial.md',
    'scripts/derive_nsc_ks_coupled_retained_control.py',
    'src/recursive_horizons/nsc_ks_retained_upstream_archive.py',
    'src/recursive_horizons/nsc_ks_history_evaluator.py',
    'src/recursive_horizons/nsc_ks_cutoff_bridge_evaluator.py',
    'src/recursive_horizons/nsc_dirac_source_phase.py',
    'src/recursive_horizons/nsc_ks_energy_propagator.py',
    'src/recursive_horizons/nsc_ks_difference_envelope.py',
    'src/recursive_horizons/nsc_ks_batched_constraints.py',
    'src/recursive_horizons/nsc_local_incoming_family.py',
)


def digest(path):
    return sha256((ROOT/path).read_bytes()).hexdigest()


def family_paths(key):
    stem = f'family-{key[0]:02d}-{key[1]:+d}'
    return DIRECTORY/(stem+'.json'), DIRECTORY/(stem+'.npz')


def selected_rectangular(proposal):
    if proposal['selected_for_first_nonlinear_trial'] != 'rectangular':
        raise ValueError('first nonlinear trial must remain the rectangular layout')
    matches = [row for row in proposal['proposals'] if row['layout'] == 'rectangular']
    if len(matches) != 1:
        raise ValueError('exactly one rectangular Newton proposal required')
    row = matches[0]
    if row['candidate_profile_identity'] != EXPECTED_IDENTITY:
        raise ValueError('rectangular candidate identity changed')
    if row.get('trial_accepted') is not False:
        raise ValueError('proposal must not pre-accept the trial')
    if not row.get('prediction_is_not_a_nonlinear_residual'):
        raise ValueError('linear prediction must remain marked as not a residual')
    return row


def context():
    base = R.context()
    proposal = json.loads(PROPOSAL.read_text())
    for path, expected in {**proposal['source_hashes'], **proposal['input_hashes']}.items():
        if digest(path) != expected:
            raise ValueError('Newton proposal input changed: '+path)
    selected = selected_rectangular(proposal)
    family = LocalIncomingFamily(np.asarray(selected['candidate']['coefficients'], float))
    identity = profile_identity(family, include_normal_window=True)
    if identity != selected['candidate_profile_identity'] or family.radius_lower_bound() <= 0:
        raise ValueError('same positive-radius rectangular candidate required')
    solve = family.collocation_nodes(R.SOLVE_COUNT)
    verify = family.collocation_nodes(R.VERIFY_COUNT)
    target = H.fixed_target_union(solve, verify)
    if target.tolist() != base['target'].tolist() or target.tolist() != proposal['z']:
        raise ValueError('trial must keep the same solve/verification union as the proposal')
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
        'proposal': proposal,
        'selected': selected,
    }


def evaluate_family(ctx, key, saved=None, allow_new=True):
    archived = ctx['archived'][key]
    interval = R.interpolation_interval(archived)
    entries = ctx['archive'].family_entries(key)
    first = next(batch for batch, _ in entries if batch.energy_sign > 0)
    if (first.group, first.angular_sign) != key:
        raise ValueError('restored family channel mismatch')
    saved_arrays = saved_meta = None
    if saved is not None:
        saved_arrays, saved_meta = saved
    evaluator = H.KSHistoryEvaluator(
        entries, baseline_gradient=ctx['baseline'], coefficients=ctx['coeff'],
        energy_interval=interval, interpolant_degree=R.DEGREE, z_grid=ctx['grid'],
        solve_nodes=ctx['solve'], verification_nodes=ctx['verify'], **R.OPTIONS)
    actual_evolve = H.evolve_energy_propagator

    def evolve(interval_e, degree, candidate, grid, target, mass, angular, rho_up, **options):
        if degree != R.DEGREE:
            raise ValueError('degree32 interpolant required')
        if tuple(float(v) for v in interval_e) != interval:
            raise ValueError('operator interval is not the archived family interval')
        if (options.get('rtol') != R.OPTIONS['rtol'] or options.get('atol') != R.OPTIONS['atol']
                or options.get('max_step') != R.OPTIONS['max_step'] or options.get('tangents') != 'all'):
            raise ValueError('solver options/grid directions changed')
        if not np.array_equal(grid, ctx['grid']) or not np.array_equal(target, ctx['target']):
            raise ValueError('computational grid or target union changed')
        if profile_identity(candidate, include_normal_window=True) != ctx['identity']:
            raise ValueError('operator history g changed')
        if (mass, angular, rho_up) != (first.mass, first.angular, first.rho_up):
            raise ValueError('operator channel mass/angular/rho_up changed')
        if saved_arrays is not None:
            operator = R.reconstruct_operator(
                interval, candidate, grid, target, mass, angular, rho_up,
                {name: saved_arrays['operator/'+name] for name in R.KEYS},
                saved_meta['operator']['diagnostics'])
            if operator.digest != saved_meta['operator']['digest']:
                raise ValueError('saved family operator digest changed')
            return operator
        if not allow_new:
            raise RuntimeError('replay must not call Dirac evolution')
        return actual_evolve(
            interval_e, degree, candidate, grid, target, mass, angular, rho_up, **options)

    with patch.object(H, 'evolve_energy_propagator', side_effect=evolve):
        result = evaluator.evaluate(ctx['family'])
    if evaluator.operator_evolution_count != 1 or len(evaluator.signed_operator_reuse) != 1:
        raise AssertionError('one operator channel per positive angular family required')
    if result['history_jacobian'].shape != (16, len(ctx['target']), 2):
        raise AssertionError('all16 retarded directions required')
    expected_ids = {f'{key[0]}_{key[1]}_E+1', f'{key[0]}_{-key[1]}_E-1'}
    signed = sorted(result['family_records'])
    if set(signed) != expected_ids or len(signed) != 2:
        raise ValueError('explicit signed partners required for this family')
    corrections = {name: np.asarray(value) for name, value in result['family_corrections'].items()}
    tangents = {name: np.asarray(value) for name, value in result['family_matter_tangents'].items()}
    if set(corrections) != expected_ids or set(tangents) != expected_ids:
        raise ValueError('family correction coverage changed')
    operator = next(iter(evaluator._operator_cache.values()))
    arrays = {f'operator/{name}': np.asarray(getattr(operator, name)) for name in R.KEYS}
    arrays['paired_matter_change'] = sum(corrections.values())
    arrays['paired_matter_tangent'] = sum(tangents.values())
    for name, value in corrections.items():
        arrays['correction/'+name] = value
    for name, value in tangents.items():
        arrays['tangent/'+name] = value
    operator_record = {
        'mass': float(operator.mass), 'angular': float(operator.angular),
        'rho_up': float(operator.rho_up), 'digest': operator.digest,
        'diagnostics': dict(operator.diagnostics),
        'reused_from_coupled_pilot': False,
    }
    meta = {
        'operator': operator_record,
        'family_records': R.jsonable(result['family_records']),
        'source_identity': R.jsonable(result['source_identity']),
    }
    arrays['metadata_json'] = np.frombuffer(json.dumps(meta, sort_keys=True).encode(), np.uint8)
    return arrays, meta, signed, operator


def new_family(ctx, key, cpu_limit=R.FAMILY_CAP):
    target_json, target_npz = family_paths(key)
    if target_json.exists() or target_npz.exists():
        raise FileExistsError('family output already exists')
    cpu = time.process_time()
    limit = min(R.FAMILY_CAP, float(cpu_limit))
    if not np.isfinite(limit) or limit <= 0:
        raise ValueError('positive remaining family CPU allocation required')
    def stop(*_):
        raise TimeoutError('one Newton-trial operator exceeded its remaining CPU allocation')
    previous = signal.signal(signal.SIGPROF, stop)
    signal.setitimer(signal.ITIMER_PROF, limit)
    try:
        arrays, meta, signed, operator = evaluate_family(ctx, key, saved=None, allow_new=True)
    finally:
        signal.setitimer(signal.ITIMER_PROF, 0)
        signal.signal(signal.SIGPROF, previous)
    archived = ctx['archived'][key]
    entries = ctx['archive'].family_entries(key)
    pos_rows = sum(len(batch.source.energies)//3 for batch, _ in entries if batch.energy_sign > 0)
    if pos_rows != archived['source_rows']:
        raise ValueError('restored rows differ from the original family coverage')
    raw = deterministic_npz_bytes(arrays)
    record = {
        'schema': 'NSC-KS-COUPLED-NEWTON-TRIAL-FAMILY-v1',
        'positive_family': list(key),
        'operator_energy_interval': list(R.interpolation_interval(archived)),
        'interpolant_degree': R.DEGREE,
        'operator_digest': operator.digest,
        'operator': meta['operator'],
        'source_rows': pos_rows,
        'source_panels': list(archived['source_panels']),
        'quadrature_measure': float(archived['quadrature_measure']),
        'signed_families': signed,
        'family_records': meta['family_records'],
        'matter_change_maxima': np.max(abs(arrays['paired_matter_change']), axis=0).tolist(),
        'matter_tangent_maxima': np.max(abs(arrays['paired_matter_tangent']), axis=(0, 1)).tolist(),
        'CPU_seconds': time.process_time()-cpu,
        'family_CPU_cap': R.FAMILY_CAP,
        'new_operator_solves': 1,
        'upstream_batches': len(entries)//2,
        'source_hashes': ctx['hashes'],
        'profile_identity': ctx['identity'],
        'physical_local_gate': 'OPEN',
        'source_and_numerical_error_bound': None,
        'payload': {'path': str(target_npz.relative_to(ROOT)),
                    'sha256': sha256(raw).hexdigest(), 'bytes': len(raw)},
    }
    tmp = target_npz.with_suffix('.npz.tmp')
    tmp.write_bytes(raw)
    tmp.replace(target_npz)
    target_json.write_text(R.dump(record))
    return record


def read_family(ctx, key, replay=False):
    path, payload = family_paths(key)
    record = json.loads(path.read_text())
    archived = ctx['archived'][key]
    if record['source_hashes'] != ctx['hashes'] or record['positive_family'] != list(key):
        raise ValueError('family signature/identity changed')
    if record.get('profile_identity') != ctx['identity']:
        raise ValueError('saved operator belongs to a different history g')
    if record['source_rows'] != archived['source_rows'] or record['source_panels'] != archived['source_panels']:
        raise ValueError('family source coverage changed')
    if list(record['operator_energy_interval']) != list(R.interpolation_interval(archived)):
        raise ValueError('family interpolation interval changed')
    if record['interpolant_degree'] != R.DEGREE:
        raise ValueError('family interpolant degree changed')
    if digest(record['payload']['path']) != record['payload']['sha256']:
        raise ValueError('family payload changed')
    if record['operator'].get('reused_from_coupled_pilot'):
        raise ValueError('Newton trial must not reuse a different-g group14 operator')
    with np.load(ROOT/record['payload']['path'], allow_pickle=False) as f:
        arrays = {name: np.array(f[name], copy=True) for name in f.files}
    meta = json.loads(arrays['metadata_json'].tobytes())
    if meta['operator']['digest'] != record['operator_digest']:
        raise ValueError('family operator metadata changed')
    if sorted(name.split('/', 1)[1] for name in arrays if name.startswith('correction/')) != record['signed_families']:
        raise ValueError('signed family coverage changed')
    if not np.array_equal(sum(value for name, value in arrays.items() if name.startswith('correction/')),
                          arrays['paired_matter_change']):
        raise ValueError('signed response sum changed')
    if not np.array_equal(sum(value for name, value in arrays.items() if name.startswith('tangent/')),
                          arrays['paired_matter_tangent']):
        raise ValueError('signed tangent sum changed')
    if replay:
        current, current_meta, signed, operator = evaluate_family(
            ctx, key, saved=(arrays, meta), allow_new=False)
        if signed != record['signed_families']:
            raise ValueError('family replay signed coverage changed')
        if operator.digest != record['operator_digest']:
            raise ValueError('family replay operator digest changed')
        for name, value in current.items():
            if name == 'metadata_json':
                continue
            if not np.array_equal(value, arrays[name]):
                raise ValueError('family replay changed: '+name)
        if R.jsonable(current_meta['family_records']) != record['family_records']:
            raise ValueError('family replay records changed')
    return record, arrays['paired_matter_change'], arrays['paired_matter_tangent'], record['family_records']


def assemble(ctx, replay=False):
    records = []
    matter = np.zeros((len(ctx['target']), 2), float)
    tangent = np.zeros((16, len(ctx['target']), 2), float)
    signed = []
    family_records = {}
    for key in sorted(ctx['archived']):
        path, payload = family_paths(key)
        if not path.exists():
            if payload.exists():
                raise ValueError('orphan operator payload retained for recovery: '+str(payload))
            continue
        record, delta, dJ, local_records = read_family(ctx, key, replay)
        records.append(record)
        matter = matter + delta
        tangent = tangent + dJ
        signed.extend(record['signed_families'])
        overlap = set(family_records) & set(local_records)
        if overlap:
            raise ValueError('signed family counted twice: '+sorted(overlap)[0])
        family_records.update(local_records)
    if len(signed) != len(set(signed)):
        raise ValueError('signed family counted twice')
    slots, slot_tangents = compatible_history_slots(ctx['family'].metric(), ctx['target'], 16)
    geometry = surface_geometry_response(slots, slot_tangents, ctx['coeff'])
    baseline = np.broadcast_to(ctx['baseline'], (len(ctx['target']), 2)).copy()
    raw_gradient = baseline + geometry['action_gradient_change'] + matter
    raw_jacobian = geometry['action_gradient_tangent'] + tangent
    complete = len(records) == len(ctx['archived'])
    weights = edge = edge_tangent = None
    phase_binding = None
    if complete:
        if len(family_records) != 120:
            raise ValueError('one hundred twenty explicit signed families required')
        weights = signed_family_angular_square_weights(family_records)
        if not np.allclose(weights, R.EXPECTED_WEIGHT, rtol=0, atol=1e-6):
            raise ValueError('signed angular-square weights must be 107952 for each energy sign')
        first = ctx['archive'].family_entries(next(iter(ctx['archived'])))[0][0]
        phase = formal_source_phase_coefficient(
            ctx['family'], ctx['target'], angular=1., rho_up=first.rho_up, gauss_nodes=PHASE_NODES)
        if phase.profile_identity != ctx['identity'] or not np.array_equal(phase.z, ctx['target']):
            raise ValueError('unit-angular phase must use the complete target union and the same g')
        if phase.delta_f_z.shape[1] != 16:
            raise ValueError('phase must retain all16 retarded directions')
        edge, edge_tangent = source_edge_from_phase(phase, weights, ctx['coeff']['a'])
        phase_binding = phase.binding
        gradient = raw_gradient + edge
        jacobian = raw_jacobian + edge_tangent
        joined_matter = matter + edge
        joined_tangent = tangent + edge_tangent
    else:
        gradient = raw_gradient
        jacobian = raw_jacobian
        joined_matter = matter
        joined_tangent = tangent
    previous = np.asarray(ctx['proposal']['constraint_maxima_with_verified_phase_values'], float)
    predicted = np.asarray(ctx['selected']['predicted_all_node_residual_maxima'], float)
    measured = np.max(abs(gradient), axis=0)
    node_maxima = np.max(abs(gradient), axis=1)
    improved = bool(complete and np.all(measured < previous) and ctx['family'].radius_lower_bound() > 0)
    status = ('OPEN: complete rectangular Newton-trial finite-source residual; '
              'uncertified physical gate' if complete
              else 'OPEN: partial rectangular Newton-trial finite-source residual')
    return {
        'schema': 'NSC-KS-COUPLED-NEWTON-TRIAL-v1',
        'accountable_author': 'Douglas Ek',
        'status': status,
        'history': ctx['family'].description(),
        'profile_identity': ctx['identity'],
        'parent_profile_identity': ctx['proposal']['current_profile_identity'],
        'layout': 'rectangular',
        'step_scale': ctx['selected']['step_scale'],
        'step_norm': ctx['selected']['step_norm'],
        'seed': {'purpose': 'geometric solver control only',
                 'w_at_center': 0., 'w_z_at_center': ctx['seed_slope'],
                 'w_zz_at_center': 0., 'U_constant': ctx['seed_U'],
                 'source_evolved_from_original_upstream': True,
                 'physical_root_claimed': False},
        'state_law': 'C_Sigma[g]=U_g C_up U_g†; original fixed source',
        'completed_positive_families': len(records),
        'required_positive_families': len(ctx['archived']),
        'all_retained_nonzero_angular_families_included': complete,
        'positive_energy_rows': sum(item['source_rows'] for item in records),
        'signed_family_ids': sorted(signed),
        'sampled_signed_families': sorted(signed),
        'grid_nodes': R.GRID_NODES,
        'energy_interpolation_degree': R.DEGREE,
        'retarded_directions': 16,
        'options': R.OPTIONS,
        'target_nodes': ctx['target'].tolist(),
        'z': ctx['target'].tolist(),
        'summed_matter_change': matter.tolist(),
        'summed_matter_tangent': tangent.tolist(),
        'raw_action_gradient': raw_gradient.tolist(),
        'raw_history_jacobian': raw_jacobian.tolist(),
        'source_cutoff_edge_gradient': None if edge is None else edge.tolist(),
        'source_cutoff_edge_history_tangent': None if edge_tangent is None else edge_tangent.tolist(),
        'action_gradient': gradient.tolist(),
        'history_jacobian': jacobian.tolist(),
        'constraint_maxima': measured.tolist(),
        'raw_constraint_maxima': np.max(abs(raw_gradient), axis=0).tolist(),
        'matter_change_maxima': np.max(abs(matter), axis=0).tolist(),
        'joined_matter_change_maxima': np.max(abs(joined_matter), axis=0).tolist(),
        'retarded_matter_plus_edge_jacobian_maxima': np.max(abs(joined_tangent), axis=(0, 1)).tolist(),
        'source_cutoff_angular_square_weights': None if weights is None else np.asarray(weights).tolist(),
        'source_cutoff_phase_binding': phase_binding,
        'source_cutoff_phase_gauss_nodes': PHASE_NODES,
        'previous_constraint_maxima': previous.tolist(),
        'linear_predicted_all_node_residual_maxima': predicted.tolist(),
        'prediction_is_not_a_nonlinear_residual': True,
        'measured_all_node_residual_maxima': measured.tolist() if complete else None,
        'measured_nodewise_residual_maxima': node_maxima.tolist() if complete else None,
        'radius_lower_bound': ctx['family'].radius_lower_bound(),
        'residual_improved_on_all_components': improved,
        'trial_accepted': False,
        'acceptance_rule': 'accept only after complete measured residual, positive radius, and explicit review; linear prediction is not acceptance',
        'family_records': [{'path': str(family_paths(tuple(item['positive_family']))[0].relative_to(ROOT)),
                            'sha256': digest(family_paths(tuple(item['positive_family']))[0]),
                            'operator_digest': item['operator_digest'],
                            'new_operator_solves': item['new_operator_solves']}
                           for item in records],
        'CPU_seconds': sum(item['CPU_seconds'] for item in records),
        'CPU_cap': R.CPU_CAP,
        'family_CPU_cap': R.FAMILY_CAP,
        'new_operator_solves': sum(item['new_operator_solves'] for item in records),
        'reused_operator_solves': 0,
        'source_hashes': ctx['hashes'],
        'scope': {
            'local_interval': 'S(1)+[.12,.18]',
            'physical_local_gate': 'OPEN',
            'physical_root_claimed': False,
            'physical_optimizer_run': False,
            'physical_EXISTENCE_certificate': False,
            'physical_NONEXISTENCE_certificate': False,
            'full_source_coverage': complete,
            'source_preparation_repeated': False,
            'source_cutoff_edge_and_tangent_included': complete,
            'source_cutoff_edge_computed_once': True,
            'baseline_and_geometry_counted_once': True,
            'energy_folding_factor': 1,
            'analytic_zero_radius_response_groups': [0, 13, 23],
            'ell0_baseline_retained': True,
            'ell0_matter_change_exact': [0., 0.],
            'group14_operators_reused_from_coupled_pilot': False,
            'old_phase_bounds_transferred': False,
            'source_cutoff_phase_quadrature_error_bound': None,
            'full_source_error_bound': None,
            'source_accuracy': None,
            'field_error_bound': None,
            'node_evolution_error_bound': None,
            'operator_interpolation_error_bound': None,
            'source_cutoff_finite_tail_error_bound': None,
            'changed_history_UV_bound': None,
            'between_node_error_bound': None,
            'finite_tail_and_full_error_bound': None,
            'metric_timestep': False,
            'new_action_term': False,
            'Gamma_rest_assigned': False,
        },
        'reproducer': 'python3 scripts/derive_nsc_ks_coupled_newton_trial.py --check',
    }


def run(cpu_budget, max_new):
    if not np.isfinite(cpu_budget) or cpu_budget <= 0 or max_new < 0:
        raise ValueError('nonnegative work bounds required')
    ctx = context()
    DIRECTORY.mkdir(parents=True, exist_ok=True)
    start = time.process_time()
    new = 0
    spent = 0.
    for key in sorted(ctx['archived']):
        path, payload = family_paths(key)
        if path.exists():
            record, _, _, _ = read_family(ctx, key)
            spent += record['CPU_seconds']
            continue
        if payload.exists():
            raise ValueError('orphan operator payload retained for recovery: '+str(payload))
        if new >= max_new:
            continue
        remaining = min(cpu_budget - (time.process_time()-start), R.CPU_CAP - spent)
        if remaining <= 1.:
            break
        try:
            record = new_family(ctx, key, cpu_limit=remaining-1.)
        except TimeoutError as error:
            print(json.dumps({'status': 'OPEN', 'timeout': str(error),
                              'positive_family': list(key),
                              'physical_NONEXISTENCE_certificate': False}), flush=True)
            break
        new += record['new_operator_solves']
        spent += record['CPU_seconds']
        print(json.dumps({k: record[k] for k in
                          ('positive_family', 'source_rows', 'matter_change_maxima',
                           'CPU_seconds', 'new_operator_solves')}), flush=True)
        progress = assemble(ctx)
        (DIRECTORY/'progress.json').write_text(R.dump(progress))
    result = assemble(ctx)
    if result['all_retained_nonzero_angular_families_included']:
        progress = DIRECTORY/'progress.json'
        if progress.exists():
            progress.unlink()
        OUTPUT.write_text(R.dump(result))
    else:
        (DIRECTORY/'progress.json').write_text(R.dump(result))
    result['runtime'] = {
        'CPU_seconds': time.process_time()-start,
        'CPU_cap': R.CPU_CAP,
        'family_CPU_seconds': result['CPU_seconds'],
        'new_operator_solves': result['new_operator_solves'],
        'reused_operator_solves': 0,
    }
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--run', action='store_true')
    mode.add_argument('--check', action='store_true')
    parser.add_argument('--cpu-budget', type=float, default=R.CPU_CAP)
    parser.add_argument('--max-new', type=int, default=MAX_NEW)
    args = parser.parse_args()
    if args.run:
        result = run(args.cpu_budget, args.max_new)
        runtime = result.get('runtime')
    else:
        ctx = context()
        result = assemble(ctx, replay=True)
        path = OUTPUT if OUTPUT.exists() else DIRECTORY/'progress.json'
        previous = json.loads(path.read_text())
        previous.pop('runtime', None)
        if R.jsonable(result) != previous:
            raise ValueError('Newton-trial retained-source replay differs')
        runtime = {'CPU_seconds': None, 'new_operator_solves': 0,
                   'reused_operator_solves': 0, 'replay_dirac_evolution': False}
    print(json.dumps({key: result[key] for key in (
        'status', 'completed_positive_families', 'required_positive_families',
        'positive_energy_rows', 'constraint_maxima', 'raw_constraint_maxima',
        'matter_change_maxima', 'previous_constraint_maxima',
        'linear_predicted_all_node_residual_maxima',
        'residual_improved_on_all_components', 'trial_accepted',
        'CPU_seconds', 'new_operator_solves', 'reused_operator_solves') if key in result}, indent=2))
    print(json.dumps(runtime, indent=2))
