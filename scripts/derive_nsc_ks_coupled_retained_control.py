#!/usr/bin/env python3
"""One full retained-source evaluation of the committed coupled w/U history."""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
import argparse
from collections.abc import Mapping
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
from recursive_horizons.nsc_dirac_source_phase import formal_source_phase_coefficient
from recursive_horizons.nsc_evolved_incoming_constraints import (
    compatible_history_slots, surface_geometry_response)
from recursive_horizons.nsc_ks_cutoff_bridge_evaluator import (
    signed_family_angular_square_weights, source_edge_from_phase)
from recursive_horizons.nsc_ks_energy_propagator import KSEnergyPropagator
from recursive_horizons.nsc_ks_profile_identity import profile_identity
from recursive_horizons.nsc_ks_retained_upstream_archive import RetainedUpstreamArchive
from recursive_horizons.nsc_ks_source_envelope import (
    KSEnvelopeBinding, computational_z_grid, _sample_axial_profiles, usual_axial_support)
from recursive_horizons.nsc_local_incoming_family import LocalIncomingFamily
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes
from derive_nsc_evolved_incoming_constraints import coefficients_from_records

DIRECTORY = ROOT/'results/development/artifacts/nsc-ks-coupled-retained-control'
OUTPUT = ROOT/'results/development/nsc-ks-coupled-retained-control.json'
PILOT_JSON = ROOT/'results/development/nsc-ks-coupled-history-pilot.json'
PILOT_NPZ = ROOT/'results/development/artifacts/nsc-ks-coupled-history-pilot.npz'
OPTIONS = {'rtol': 2e-11, 'atol': 2e-13, 'max_step': .002}
DEGREE = 32
GRID_NODES = 64
SOLVE_COUNT = 8
VERIFY_COUNT = 17
PHASE_NODES = 32
FAMILY_CAP = 40.
CPU_CAP = 600.
MAX_NEW = 58
EXPECTED_WEIGHT = 107952.
SEED_SLOPE = .003446010557459475
SEED_U = 12.12435557608503
KEYS = ('energies', 'z', 'reference', 'difference', 'difference_z', 'tangent', 'tangent_z')
OWNERS = (
    'scripts/derive_nsc_ks_coupled_retained_control.py',
    'docs/nsc-ks-coupled-retained-control.md',
    'src/recursive_horizons/nsc_ks_retained_upstream_archive.py',
    'src/recursive_horizons/nsc_ks_history_evaluator.py',
    'src/recursive_horizons/nsc_ks_cutoff_bridge_evaluator.py',
    'src/recursive_horizons/nsc_dirac_source_phase.py',
    'src/recursive_horizons/nsc_ks_energy_propagator.py',
    'src/recursive_horizons/nsc_ks_difference_envelope.py',
    'src/recursive_horizons/nsc_ks_batched_constraints.py',
)


def digest(path):
    return sha256((ROOT/path).read_bytes()).hexdigest()


def jsonable(value):
    if isinstance(value, Mapping):
        return {str(key): jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(item) for item in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    return value


def dump(record):
    return json.dumps(jsonable(record), sort_keys=True, indent=2, allow_nan=False)+'\n'


def family_paths(key):
    stem = f'family-{key[0]:02d}-{key[1]:+d}'
    return DIRECTORY/(stem+'.json'), DIRECTORY/(stem+'.npz')


def interpolation_interval(record):
    lo, hi = (float(v) for v in record['operator_energy_interval'])
    measure = float(record['quadrature_measure'])
    interval = (0., hi) if lo == 0. else (0., measure)
    if interval not in ((0., 160.), (0., 320.)):
        raise ValueError('original family interpolation interval must be [0,160] or [0,320]')
    if abs(measure-interval[1]) > 1e-11:
        raise ValueError('archived quadrature measure must match the positive interpolation endpoint')
    return interval


def reconstruct_operator(interval, family, grid, target, mass, angular, rho_up, arrays, diagnostics):
    metric = family.metric()
    w, u = _sample_axial_profiles(metric.directions, grid)
    binding = KSEnvelopeBinding(
        grid, w, u, tuple(metric.amplitudes),
        tuple(d.inner_radius for d in metric.directions),
        tuple(d.outer_radius for d in metric.directions),
        usual_axial_support(), rho_up, 1., **OPTIONS)
    operator = KSEnergyPropagator(
        interval, *(np.asarray(arrays[key]) for key in KEYS),
        mass, angular, rho_up, binding, diagnostics)
    if not np.array_equal(operator.z, target):
        raise ValueError('saved operator target nodes changed')
    if operator.tangent.shape[1] != 16:
        raise ValueError('all16 retarded directions required')
    return operator


def context():
    archive = RetainedUpstreamArchive(ROOT)
    pilot = json.loads(PILOT_JSON.read_text())
    if digest(PILOT_NPZ.relative_to(ROOT)) != pilot['payload']['sha256']:
        raise ValueError('coupled-history pilot payload changed')
    if (pilot['energy_interpolation_degree'] != DEGREE or pilot['grid_nodes'] != GRID_NODES
            or list(pilot['energy_interval']) != [0., 160.] or pilot['options'] != OPTIONS
            or pilot['retarded_directions'] != 16):
        raise ValueError('committed group14 pilot does not use the required n8/degree32/64-grid control')
    for name in ('nsc-incoming-surface-coefficients.json',
                 'nsc-incoming-surface-regular-branch.json',
                 'nsc-incoming-source-update-v5.json'):
        archive._record('results/development/'+name)
    coefficient = json.loads((ROOT/'results/development/nsc-incoming-surface-coefficients.json').read_text())
    branch = json.loads((ROOT/'results/development/nsc-incoming-surface-regular-branch.json').read_text())
    source = json.loads((ROOT/'results/development/nsc-incoming-source-update-v5.json').read_text())
    coeff = coefficients_from_records(coefficient, branch)
    baseline = np.asarray(source['baseline']['action_gradient_approximant'], float)
    slope = -baseline[1]/coeff['F'][0]
    U = -(baseline[0]+coeff['C']*slope*slope)/coeff['A0']
    if abs(slope-SEED_SLOPE) > 1e-15 or abs(U-SEED_U) > 1e-15:
        raise ValueError('geometric solver-control seed changed')
    if abs(pilot['seed']['w_z_at_center']-SEED_SLOPE) > 1e-15 or abs(pilot['seed']['U_constant']-SEED_U) > 1e-15:
        raise ValueError('committed pilot seed is not the declared geometric control')
    family = LocalIncomingFamily(np.asarray(pilot['history']['coefficients'], float))
    identity = profile_identity(family, include_normal_window=True)
    if identity != pilot['profile_identity']:
        raise ValueError('committed coupled-history g changed')
    grid = computational_z_grid(GRID_NODES)
    solve = family.collocation_nodes(SOLVE_COUNT)
    verify = family.collocation_nodes(VERIFY_COUNT)
    target = H.fixed_target_union(solve, verify)
    if target.tolist() != pilot['target_nodes']:
        raise ValueError('solve/verification union differs from the committed group14 pilot')
    with np.load(PILOT_NPZ, allow_pickle=False) as f:
        pilot_arrays = {key: np.array(f[key], copy=True) for key in f.files}
    pilot_meta = json.loads(pilot_arrays.pop('metadata_json').tobytes())
    if len(pilot_meta['operators']) != 2:
        raise ValueError('pilot must retain both original group14 angular operators')
    archived = {}
    for key in archive.family_keys:
        record = archive._families[key]
        archive._payload(record['payload'])
        archived[key] = record
    if len(archived) != 60:
        raise ValueError('sixty original nonzero-angular families required')
    hashes = {
        **archive.input_hashes,
        **{path: digest(path) for path in OWNERS},
        str(PILOT_JSON.relative_to(ROOT)): digest(PILOT_JSON.relative_to(ROOT)),
        str(PILOT_NPZ.relative_to(ROOT)): digest(PILOT_NPZ.relative_to(ROOT)),
    }
    ell0 = archive.analytic_zero_channels()
    if tuple(channel['index'] for channel in ell0) != (0, 13, 23):
        raise ValueError('ell=0 groups 0,13,23 required')
    return {
        'archive': archive, 'pilot': pilot, 'pilot_arrays': pilot_arrays,
        'pilot_meta': pilot_meta, 'family': family, 'identity': identity,
        'coeff': coeff, 'baseline': baseline, 'grid': grid, 'solve': solve,
        'verify': verify, 'target': target, 'archived': archived,
        'hashes': hashes, 'ell0': ell0, 'seed_slope': float(slope), 'seed_U': float(U),
    }


def operator_from_pilot(ctx, mass, angular, rho_up, interval):
    matches = [index for index, rec in enumerate(ctx['pilot_meta']['operators'])
               if rec['mass'] == mass and rec['angular'] == angular and rec['rho_up'] == rho_up]
    if len(matches) != 1:
        raise ValueError('exactly one committed group14 operator required for this channel')
    index = matches[0]
    arrays = {key: ctx['pilot_arrays'][f'operator{index}/{key}'] for key in KEYS}
    operator = reconstruct_operator(
        interval, ctx['family'], ctx['grid'], ctx['target'], mass, angular, rho_up,
        arrays, ctx['pilot_meta']['operators'][index]['diagnostics'])
    expected = ctx['pilot_meta']['operators'][index]['digest']
    if operator.digest != expected:
        raise ValueError('reconstructed group14 operator digest differs from the committed pilot')
    return operator, {'mass': mass, 'angular': angular, 'rho_up': rho_up,
                      'digest': operator.digest,
                      'diagnostics': dict(operator.diagnostics),
                      'reused_from_coupled_pilot': True,
                      'pilot_operator_index': index}


def evaluate_family(ctx, key, saved=None, allow_new=True):
    archived = ctx['archived'][key]
    interval = interpolation_interval(archived)
    entries = ctx['archive'].family_entries(key)
    first = next(batch for batch, _ in entries if batch.energy_sign > 0)
    if (first.group, first.angular_sign) != key:
        raise ValueError('restored family channel mismatch')
    reused = key[0] == 14
    if reused and interval != (0., 160.):
        raise ValueError('group14 must reuse the committed [0,160] pilot interpolant')
    saved_arrays = saved_meta = None
    if saved is not None:
        saved_arrays, saved_meta = saved
    evaluator = H.KSHistoryEvaluator(
        entries, baseline_gradient=ctx['baseline'], coefficients=ctx['coeff'],
        energy_interval=interval, interpolant_degree=DEGREE, z_grid=ctx['grid'],
        solve_nodes=ctx['solve'], verification_nodes=ctx['verify'], **OPTIONS)
    actual_evolve = H.evolve_energy_propagator

    def evolve(interval_e, degree, candidate, grid, target, mass, angular, rho_up, **options):
        if degree != DEGREE:
            raise ValueError('degree32 interpolant required')
        if tuple(float(v) for v in interval_e) != interval:
            raise ValueError('operator interval is not the archived family interval')
        if (options.get('rtol') != OPTIONS['rtol'] or options.get('atol') != OPTIONS['atol']
                or options.get('max_step') != OPTIONS['max_step'] or options.get('tangents') != 'all'):
            raise ValueError('solver options/grid directions changed')
        if not np.array_equal(grid, ctx['grid']) or not np.array_equal(target, ctx['target']):
            raise ValueError('computational grid or target union changed')
        if profile_identity(candidate, include_normal_window=True) != ctx['identity']:
            raise ValueError('operator history g changed')
        if (mass, angular, rho_up) != (first.mass, first.angular, first.rho_up):
            raise ValueError('operator channel mass/angular/rho_up changed')
        if saved_arrays is not None:
            operator = reconstruct_operator(
                interval, candidate, grid, target, mass, angular, rho_up,
                {key_name: saved_arrays['operator/'+key_name] for key_name in KEYS},
                saved_meta['operator']['diagnostics'])
            if operator.digest != saved_meta['operator']['digest']:
                raise ValueError('saved family operator digest changed')
            return operator
        if reused:
            operator, _ = operator_from_pilot(ctx, mass, angular, rho_up, interval)
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
    arrays = {f'operator/{name}': np.asarray(getattr(operator, name)) for name in KEYS}
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
        'reused_from_coupled_pilot': bool(reused),
    }
    if reused:
        operator_record['pilot_digest'] = operator.digest
    meta = {
        'operator': operator_record,
        'family_records': jsonable(result['family_records']),
        'source_identity': jsonable(result['source_identity']),
    }
    arrays['metadata_json'] = np.frombuffer(json.dumps(meta, sort_keys=True).encode(), np.uint8)
    return arrays, meta, signed, operator


def new_family(ctx, key, cpu_limit=FAMILY_CAP):
    target_json, target_npz = family_paths(key)
    if target_json.exists() or target_npz.exists():
        raise FileExistsError('family output already exists')
    cpu = time.process_time()
    limit = min(FAMILY_CAP, float(cpu_limit))
    if not np.isfinite(limit) or limit <= 0:
        raise ValueError('positive remaining family CPU allocation required')
    def stop(*_):
        raise TimeoutError('one coupled retained operator exceeded its remaining CPU allocation')
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
        'schema': 'NSC-KS-COUPLED-RETAINED-FAMILY-v1',
        'positive_family': list(key),
        'operator_energy_interval': list(interpolation_interval(archived)),
        'interpolant_degree': DEGREE,
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
        'family_CPU_cap': FAMILY_CAP,
        'new_operator_solves': 0 if key[0] == 14 else 1,
        'upstream_batches': len(entries)//2,
        'source_hashes': ctx['hashes'],
        'physical_local_gate': 'OPEN',
        'source_and_numerical_error_bound': None,
        'payload': {'path': str(target_npz.relative_to(ROOT)),
                    'sha256': sha256(raw).hexdigest(), 'bytes': len(raw)},
    }
    tmp = target_npz.with_suffix('.npz.tmp')
    tmp.write_bytes(raw)
    tmp.replace(target_npz)
    target_json.write_text(dump(record))
    return record


def read_family(ctx, key, replay=False):
    path, payload = family_paths(key)
    record = json.loads(path.read_text())
    archived = ctx['archived'][key]
    if record['source_hashes'] != ctx['hashes'] or record['positive_family'] != list(key):
        raise ValueError('family signature/identity changed')
    if record['source_rows'] != archived['source_rows'] or record['source_panels'] != archived['source_panels']:
        raise ValueError('family source coverage changed')
    if list(record['operator_energy_interval']) != list(interpolation_interval(archived)):
        raise ValueError('family interpolation interval changed')
    if record['interpolant_degree'] != DEGREE:
        raise ValueError('family interpolant degree changed')
    if digest(record['payload']['path']) != record['payload']['sha256']:
        raise ValueError('family payload changed')
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
        if jsonable(current_meta['family_records']) != record['family_records']:
            raise ValueError('family replay records changed')
        if key[0] == 14:
            reused, _ = operator_from_pilot(
                ctx, operator.mass, operator.angular, operator.rho_up,
                interpolation_interval(archived))
            if reused.digest != operator.digest:
                raise ValueError('replayed group14 operator is not the committed pilot operator')
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
        if not np.allclose(weights, EXPECTED_WEIGHT, rtol=0, atol=1e-6):
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
    status = ('OPEN: complete coupled-history retained finite-source control; '
              'uncertified physical gate' if complete
              else 'OPEN: partial coupled-history retained finite-source control')
    return {
        'schema': 'NSC-KS-COUPLED-RETAINED-CONTROL-v1',
        'accountable_author': 'Douglas Ek',
        'status': status,
        'history': ctx['family'].description(),
        'profile_identity': ctx['identity'],
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
        'grid_nodes': GRID_NODES,
        'energy_interpolation_degree': DEGREE,
        'retarded_directions': 16,
        'options': OPTIONS,
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
        'constraint_maxima': np.max(abs(gradient), axis=0).tolist(),
        'raw_constraint_maxima': np.max(abs(raw_gradient), axis=0).tolist(),
        'matter_change_maxima': np.max(abs(matter), axis=0).tolist(),
        'joined_matter_change_maxima': np.max(abs(joined_matter), axis=0).tolist(),
        'retarded_matter_plus_edge_jacobian_maxima': np.max(abs(joined_tangent), axis=(0, 1)).tolist(),
        'source_cutoff_angular_square_weights': None if weights is None else np.asarray(weights).tolist(),
        'source_cutoff_phase_binding': phase_binding,
        'source_cutoff_phase_gauss_nodes': PHASE_NODES,
        'family_records': [{'path': str(family_paths(tuple(item['positive_family']))[0].relative_to(ROOT)),
                            'sha256': digest(family_paths(tuple(item['positive_family']))[0]),
                            'operator_digest': item['operator_digest'],
                            'new_operator_solves': item['new_operator_solves']}
                           for item in records],
        'CPU_seconds': sum(item['CPU_seconds'] for item in records),
        'CPU_cap': CPU_CAP,
        'family_CPU_cap': FAMILY_CAP,
        'new_operator_solves': sum(item['new_operator_solves'] for item in records),
        'reused_operator_solves': sum(1 for item in records if item['new_operator_solves'] == 0),
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
            'group14_operators_reused_from_coupled_pilot': True,
            'full_source_error_bound': None,
            'source_accuracy': None,
            'field_error_bound': None,
            'node_evolution_error_bound': None,
            'operator_interpolation_error_bound': None,
            'source_cutoff_finite_tail_error_bound': None,
            'source_cutoff_phase_quadrature_error_bound': None,
            'changed_history_UV_bound': None,
            'between_node_error_bound': None,
            'finite_tail_and_full_error_bound': None,
            'metric_timestep': False,
            'new_action_term': False,
            'Gamma_rest_assigned': False,
        },
        'reproducer': 'python3 scripts/derive_nsc_ks_coupled_retained_control.py --check',
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
        needed = 0 if key[0] == 14 else 1
        if needed and new >= max_new:
            continue
        remaining = min(cpu_budget - (time.process_time()-start), CPU_CAP - spent)
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
        (DIRECTORY/'progress.json').write_text(dump(progress))
    result = assemble(ctx)
    if result['all_retained_nonzero_angular_families_included']:
        progress = DIRECTORY/'progress.json'
        if progress.exists():
            progress.unlink()
        OUTPUT.write_text(dump(result))
    else:
        (DIRECTORY/'progress.json').write_text(dump(result))
    result['runtime'] = {
        'CPU_seconds': time.process_time()-start,
        'CPU_cap': CPU_CAP,
        'family_CPU_seconds': result['CPU_seconds'],
        'new_operator_solves': result['new_operator_solves'],
        'reused_operator_solves': result['reused_operator_solves'],
    }
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--run', action='store_true')
    mode.add_argument('--check', action='store_true')
    parser.add_argument('--cpu-budget', type=float, default=CPU_CAP)
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
        if jsonable(result) != previous:
            raise ValueError('coupled retained-source replay differs')
        runtime = {'CPU_seconds': None, 'new_operator_solves': 0,
                   'reused_operator_solves': 0, 'replay_dirac_evolution': False}
    print(json.dumps({key: result[key] for key in (
        'status', 'completed_positive_families', 'required_positive_families',
        'positive_energy_rows', 'constraint_maxima', 'raw_constraint_maxima',
        'matter_change_maxima', 'CPU_seconds', 'new_operator_solves',
        'reused_operator_solves')}, indent=2))
    print(json.dumps(runtime, indent=2))
