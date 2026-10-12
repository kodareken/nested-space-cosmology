#!/usr/bin/env python3
"""Two bounded exact-phase buffer runs; --check replays saved arrays only."""
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

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
sys.path.insert(0, str(ROOT / 'scripts'))
from recursive_horizons.nsc_evolved_incoming_state import FixedSourcePreparation
from recursive_horizons.nsc_incoming_boundary_buffer import (
    ATOL, COARSE_POINTS, CPU_CAP, DRIFT_TARGET, EXTENDED_RIGHT, FINE_POINTS,
    FINE_SPACING, MAX_STEP, MULTIPLICITY, ORIGINAL_FINE_POINTS, ORIGINAL_RIGHT,
    RTOL, SAMPLE_TIMES, TIME_STOP, append_extension, characteristic_from_ks_path,
    continue_canonical_amplitudes, continuum_incoming_buffer, exact_harmonic_reference,
    extended_grid, fiber_gram_residuals, initial_stationary_defect, ks_from_characteristic,
    matter_drift, maxabs, nested_prefix, original_grid, pg_from_ks_amplitudes,
    pg_from_characteristic, prefixes_equal, sigma_parameters, slice_state_residuals,
)
from recursive_horizons.nsc_ks_spacetime_variation import chart_coordinates
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes
from recursive_horizons.nsc_retarded_radial_response import local_frame_checks
from recursive_horizons.nsc_transmitting_history_jets import FourthOrderModePropagator

OUTPUT = 'results/development/nsc-incoming-boundary-buffer.json'
ARTIFACT = 'results/development/artifacts/nsc-incoming-boundary-buffer.npz'
OWNED = (
    'src/recursive_horizons/nsc_incoming_boundary_buffer.py',
    'scripts/derive_nsc_incoming_boundary_buffer.py',
    'tests/test_nsc_incoming_boundary_buffer.py',
    'docs/nsc-incoming-boundary-buffer.md',
)
INPUTS = (
    'results/development/nsc-retarded-response-resolution.json',
    'results/development/nsc-exact-source-phase.json',
    'results/development/nsc-incoming-fixed-transfer.json',
    'src/recursive_horizons/nsc_retarded_radial_response.py',
    'src/recursive_horizons/nsc_common_ks_trace.py',
    'src/recursive_horizons/nsc_transmitting_dirac_domain.py',
    'src/recursive_horizons/nsc_ks_spacetime_variation.py',
    'src/recursive_horizons/nsc_transmitting_history_jets.py',
    'src/recursive_horizons/nsc_transmitting_history_modes.py',
    'src/recursive_horizons/nsc_evolved_incoming_state.py',
    'src/recursive_horizons/nsc_evolved_incoming_constraints.py',
    'src/recursive_horizons/nsc_lorentzian.py',
    'scripts/derive_nsc_exact_source_phase.py',
    'scripts/derive_nsc_retarded_response_resolution.py',
)
CASES = (('coarse', COARSE_POINTS), ('fine', FINE_POINTS))


class CPUExceeded(RuntimeError):
    pass


def digest(path):
    return sha256((ROOT / path).read_bytes()).hexdigest()


def signatures():
    return {path: digest(path) for path in (*OWNED, *INPUTS)}


def case_arrays(arrays, name):
    return {key[len(name) + 1:]: value for key, value in arrays.items() if key.startswith(name + '/')}


def inputs():
    resolution = json.loads((ROOT / INPUTS[0]).read_text())
    exact = json.loads((ROOT / INPUTS[1]).read_text())
    transfer = json.loads((ROOT / INPUTS[2]).read_text())
    payload = resolution['payload']
    if digest(payload['path']) != payload['sha256']:
        raise ValueError('resolution artifact changed')
    if digest(exact['payload']['path']) != exact['payload']['sha256']:
        raise ValueError('exact-source-phase artifact changed')
    channel = transfer['actual_source_bindings']['channel']
    multiplicity = channel['copy_count'] * channel['degeneracy'] / 2
    if multiplicity != MULTIPLICITY:
        raise ValueError('canonical multiplicity 12 required')
    if channel['compact_mass'] <= 0 or channel['angular_eigenvalue'] <= 0:
        raise ValueError('inherited signed family required')
    with np.load(ROOT / payload['path'], allow_pickle=False) as archive:
        x = np.array(archive['reference/x'], copy=True)
        phi = np.array(archive['reference/phi'], copy=True)
        data = {key: np.array(archive['source/' + key], copy=True)
                for key in ('energies', 'weights', 'source', 'projector', 'mass', 'angular')}
    if not np.array_equal(x, original_grid(ORIGINAL_FINE_POINTS)):
        raise ValueError('saved 3201 prefix coordinates changed')
    if phi.shape != (2 * ORIGINAL_FINE_POINTS, 12):
        raise ValueError('saved 3201 prefix columns changed')
    archived_x, archived_phi = nested_prefix(x, phi, 801)
    return {
        'x': x, 'phi': phi, 'data': data, 'archived_x': archived_x, 'archived_phi': archived_phi,
        'resolution_payload': payload, 'exact_payload': exact['payload'],
        'old_exact_phase_drift': exact['measurements']['PDE_matter_drift'],
        'old_midpoint_forcing_drift': exact['control']['old_midpoint_forcing_drift'],
        'channel': channel, 'multiplicity': float(multiplicity),
    }


def analyze(arrays, meta):
    if 'extension/amplitudes' not in arrays:
        return {'prefix': {}, 'per_case': {}, 'maximum_PDE_matter_drift': None,
                'maximum_phase_or_flux_residual': None, 'partial_control_target': DRIFT_TARGET,
                'numerical_partial_control_passed': False, 'field_runs': 0,
                'field_runs_attempted': int(meta.get('field_runs_attempted', 0)),
                'preparation_completed': False}
    times = arrays['times']
    C = arrays['source/covariance']
    weights = arrays['source/column_weights']
    projector = arrays['source/projector']
    original_x, original_phi = arrays['original/x'], arrays['original/phi']
    archived_x, archived_phi = nested_prefix(original_x, original_phi, 801)
    extension = fiber_gram_residuals(
        np.concatenate((arrays['extension/start_amplitudes'][None], arrays['extension/amplitudes'])),
        projector)
    rows = {}
    prefix = {
        'fine_prefix_equal': prefixes_equal(original_x, original_phi, arrays['fine/x'], arrays['fine/phi'])
        if 'fine/x' in arrays else False,
        'coarse_prefix_equal': False,
        'archived801_nested_residual': maxabs(archived_phi - arrays.get('archived/phi', archived_phi)),
        'rho1_source_columns_unchanged': True,
    }
    if 'coarse/x' in arrays:
        x801, phi801 = nested_prefix(original_x, original_phi, 801)
        prefix['coarse_prefix_equal'] = prefixes_equal(x801, phi801, arrays['coarse/x'], arrays['coarse/phi'])
        prefix['archived801_nested_residual'] = maxabs(phi801 - arrays['archived/phi']) if 'archived/phi' in arrays else maxabs(phi801 - archived_phi)
    for name, _points in CASES:
        if name not in meta.get('completed', []):
            continue
        case = case_arrays(arrays, name)
        parameters = dict(meta['parameters'])
        drift = matter_drift(case['columns'], case['pde_axial_columns'], case['reference_columns'],
                             case['reference_axial_columns'], C, weights, parameters, times)
        F0 = case['columns'][0]
        original_index = int(np.flatnonzero(original_x == 1.)[0])
        take = np.array([original_index, len(original_x) + original_index])
        restriction = case['restriction']
        expected0 = restriction @ original_phi[take]
        if name == 'coarse':
            x801, phi801 = nested_prefix(original_x, original_phi, 801)
            ix = int(np.flatnonzero(x801 == 1.)[0])
            expected0 = restriction @ phi801[np.array([ix, 801 + ix])]
        rho1_match = maxabs(F0 - expected0)
        prefix['rho1_source_columns_unchanged'] = bool(
            prefix['rho1_source_columns_unchanged'] and rho1_match == 0.)
        axial_drift = np.max(np.abs(case['pde_axial_columns'] - case['reference_axial_columns']), axis=(1, 2))
        peak = int(np.argmax(axial_drift))
        rows[name] = {
            'grid_points': int(len(case['x'])),
            'right_edge': float(case['x'][-1]),
            'spacing': float(np.diff(case['x'])[0]),
            'exact_harmonic_phase_residual': maxabs(case['phase'] - case['expected_phase']),
            'flux_algebra': float(meta['cases'][name]['flux_algebra']),
            'incoming_column_drift': maxabs(case['columns'] - case['reference_columns']),
            'pde_axial_column_drift': maxabs(case['pde_axial_columns'] - case['reference_axial_columns']),
            'incoming_derivative_peak_time': float(times[peak]),
            'incoming_derivative_peak': float(axial_drift[peak]),
            'rho1_initial_column_residual': float(rho1_match),
            'stationary_defect': meta['cases'][name]['stationary_defect'],
            **drift,
        }
    algebra = [row['exact_harmonic_phase_residual'] for row in rows.values()] + [
        row['flux_algebra'] for row in rows.values()]
    drifts = [component for row in rows.values() for component in row['PDE_matter_drift']]
    completed = list(meta.get('completed', []))
    numerical_pass = (
        completed == ['coarse', 'fine']
        and prefix['fine_prefix_equal'] and prefix['coarse_prefix_equal']
        and prefix['rho1_source_columns_unchanged']
        and (max(algebra) if algebra else np.inf) < DRIFT_TARGET
        and rows['fine']['meets_partial_control_target']
        and max(extension.values()) < DRIFT_TARGET
        and max(meta['local_frame'].values()) < DRIFT_TARGET
        and meta['causal_buffer']['causal_margin_lower_bound'] > 0
        and not meta.get('cpu_budget_exceeded', False)
    )
    return {
        'prefix': prefix,
        'extension_fiber_gram': extension,
        'local_frame': meta.get('local_frame', {}),
        'restriction_roundtrip': meta.get('restriction_roundtrip', {}),
        'slice_residuals': meta.get('slice_residuals', {}),
        'causal_buffer': meta.get('causal_buffer', {}),
        'per_case': rows,
        'maximum_phase_or_flux_residual': max(algebra) if algebra else None,
        'maximum_PDE_matter_drift': max(drifts) if drifts else None,
        'partial_control_target': DRIFT_TARGET,
        'numerical_partial_control_passed': numerical_pass,
        'field_runs': len(completed),
        'field_runs_attempted': int(meta.get('field_runs_attempted', len(completed))),
        'preparation_completed': True,
        'fine_grid_target_met': bool(rows.get('fine', {}).get('meets_partial_control_target', False)),
        'coarse_grid_role': 'comparison, not a second fine-grid accuracy requirement',
    }


def record(arrays, meta, path, digest_value, size):
    checks = analyze(arrays, meta)
    passed = checks['numerical_partial_control_passed']
    status = ('PASS: fine-grid numerical boundary-buffer drift below 3e-11 on this partial source control; '
              'physical/full-source OPEN' if passed else
              'OPEN: numerical boundary-buffer control completed; physical/full-source remain OPEN')
    if meta.get('cpu_budget_exceeded'):
        status = 'OPEN: CPU cap reached; no automatic extra field case'
    return {
        'schema': 'NSC-INCOMING-BOUNDARY-BUFFER-v1',
        'accountable_author': 'Douglas Ek',
        'status': status,
        'control': {
            'family': '14_1',
            'history_amplitude': 0.,
            'field_runs_authorized': 2,
            'field_runs': checks['field_runs'],
            'field_runs_attempted': checks['field_runs_attempted'],
            'cases': ['coarse 951/65 h=.004 x=[-2,1.8]', 'fine 3801/65 h=.001 x=[-2,1.8]'],
            'sample_times': SAMPLE_TIMES,
            'time_coordinates': 'inherited numerical [0,.3], not a selected physical duration',
            'method': 'shared dense DOP853 KS continuation of archived columns at 1.2, then exact harmonic phase block',
            'axial_derivative': 'actual semi-discrete PDE at fixed rho1; no cubic substitution',
            'integration': {'rtol': RTOL, 'atol': ATOL, 'max_step': MAX_STEP,
                            'interval': [ORIGINAL_RIGHT, EXTENDED_RIGHT], 'energies': 4, 'coupled_ode': 1},
            'CPU_budget_seconds': CPU_CAP,
            'partial_control_target': DRIFT_TARGET,
            'target_applies_to': 'fine 3801 grid; coarse 951 is the spatial comparison',
            'canonical_multiplicity': MULTIPLICITY,
            'old_exact_phase_drift': meta.get('old_exact_phase_drift'),
            'old_midpoint_forcing_drift': meta.get('old_midpoint_forcing_drift'),
        },
        'measurements': checks,
        'runtime': {key: meta[key] for key in
                    ('CPU_seconds', 'wall_seconds', 'cpu_budget_seconds', 'cpu_budget_exceeded',
                     'case_CPU_seconds', 'continuation_CPU_seconds', 'stop_reason') if key in meta},
        'scope': {
            'source_law_changed': False,
            'source_phase_fitted': False,
            'column_normalized': False,
            'new_source_covariance': False,
            'interior_counterforce': False,
            'stress_drift_subtracted': False,
            'incoming_C0_imposed': False,
            'metric_evolution': False,
            'physical_constraint_solution': False,
            'full_source_error_bound': None,
            'continuum_field_error_bound': None,
            'finite_grid_exact_causality': False,
            'physical_Cup_changed': False,
            'new_horizon_or_source_grid': False,
            'nonzero_history_accuracy_verified': False,
            'publication': False,
            'mode_accuracy_from_Gram': False,
        },
        'source_hashes': {path: meta['signature'][path] for path in OWNED},
        'input_hashes': {path: meta['signature'][path] for path in INPUTS},
        'reused_artifacts': meta.get('reused_artifacts', []),
        'payload': {'path': path, 'sha256': digest_value, 'bytes': size},
        'reproducer': 'python3 scripts/derive_nsc_incoming_boundary_buffer.py --check',
    }


def publish(arrays, meta):
    encoded = {**arrays, 'metadata_json': np.frombuffer(json.dumps(meta, sort_keys=True).encode(), np.uint8)}
    raw = deterministic_npz_bytes(encoded)
    (ROOT / ARTIFACT).write_bytes(raw)
    result = record(encoded, meta, ARTIFACT, sha256(raw).hexdigest(), len(raw))
    (ROOT / OUTPUT).write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + '\n')
    return result


def run():
    if (ROOT / OUTPUT).exists() or (ROOT / ARTIFACT).exists():
        raise FileExistsError('buffer output already exists; use --check, never repeat this control implicitly')
    before = signatures()
    loaded = inputs()
    start = time.process_time()
    wall = time.monotonic()
    arrays = {}
    meta = {
        'signature': before, 'completed': [], 'cases': {}, 'case_CPU_seconds': {},
        'cpu_budget_seconds': CPU_CAP, 'cpu_budget_exceeded': False, 'stop_reason': 'checkpoint',
        'old_exact_phase_drift': loaded['old_exact_phase_drift'],
        'old_midpoint_forcing_drift': loaded['old_midpoint_forcing_drift'],
        'reused_artifacts': [loaded['resolution_payload'], loaded['exact_payload']],
        'parameters': None, 'field_runs_attempted': 0,
    }

    def stop(*_):
        raise CPUExceeded('60 CPU-second incoming-boundary-buffer budget exhausted')

    previous = signal.signal(signal.SIGPROF, stop)
    signal.setitimer(signal.ITIMER_PROF, CPU_CAP)

    def checkpoint(reason):
        meta.update(CPU_seconds=time.process_time() - start, wall_seconds=time.monotonic() - wall,
                    stop_reason=reason)
        if signatures() != before:
            raise ValueError('buffer owners changed during execution')
        return publish(arrays, meta)

    try:
        data = loaded['data']
        energies = data['energies']
        mass = float(data['mass'].item())
        angular = float(data['angular'].item())
        source = FixedSourcePreparation.from_signed_blocks(energies, data['weights'], data['source'])
        original_x, original_phi = loaded['x'], loaded['phi']
        index = int(np.flatnonzero(original_x == ORIGINAL_RIGHT)[0])
        if index != len(original_x) - 1:
            raise ValueError('saved 3201 prefix must end at rho=1.2')
        start_amplitudes, chart = ks_from_characteristic(original_x, original_phi, energies, index)
        meta['restriction_roundtrip'] = {
            'clock_shift': float(chart['clock_shift']),
            'KS_time': float(chart['KS_time']),
            'start_pg_residual': maxabs(
                pg_from_ks_amplitudes(ORIGINAL_RIGHT, energies, start_amplitudes)
                - pg_from_characteristic(original_phi, index)),
        }
        meta['local_frame'] = local_frame_checks((1., ORIGINAL_RIGHT, 1.5, EXTENDED_RIGHT),
                                                 energies, mass, angular)
        meta['slice_residuals'] = {
            'rho_1.2': slice_state_residuals(energies, pg_from_ks_amplitudes(
                ORIGINAL_RIGHT, energies, start_amplitudes), data['source'], data['projector'], ORIGINAL_RIGHT),
        }
        x_fine = extended_grid(FINE_POINTS)
        new_rho = x_fine[x_fine > ORIGINAL_RIGHT]
        begin = time.process_time()
        sampled, solver = continue_canonical_amplitudes(
            energies, mass, angular, start_amplitudes, new_rho, rtol=RTOL, atol=ATOL, max_step=MAX_STEP)
        meta['continuation_CPU_seconds'] = time.process_time() - begin
        meta['solver'] = solver
        extension = characteristic_from_ks_path(new_rho, energies, sampled)
        phi_fine = append_extension(original_x, original_phi, x_fine, extension)
        if not prefixes_equal(original_x, original_phi, x_fine, phi_fine):
            raise ArithmeticError('fine prefix is not the saved 3201 bytes')
        x_coarse, phi_coarse = nested_prefix(x_fine, phi_fine, COARSE_POINTS)
        x801, phi801 = nested_prefix(original_x, original_phi, 801)
        if not prefixes_equal(x801, phi801, x_coarse, phi_coarse):
            raise ArithmeticError('coarse prefix is not the nested 801 subset of the saved 3201 bytes')
        arrays.update({
            'times': np.linspace(0., TIME_STOP, SAMPLE_TIMES),
            'z': np.linspace(0., TIME_STOP, SAMPLE_TIMES) + chart_coordinates(1.)[1],
            'original/x': original_x, 'original/phi': original_phi,
            'archived/x': loaded['archived_x'], 'archived/phi': loaded['archived_phi'],
            'source/energies': energies, 'source/weights': data['weights'],
            'source/source': data['source'], 'source/projector': data['projector'],
            'source/mass': data['mass'], 'source/angular': data['angular'],
            'source/covariance': source.covariance, 'source/column_weights': source.column_weights,
            'extension/rho': new_rho, 'extension/amplitudes': sampled,
            'extension/start_amplitudes': start_amplitudes,
            'coarse/x': x_coarse, 'coarse/phi': phi_coarse,
            'fine/x': x_fine, 'fine/phi': phi_fine,
        })
        meta['causal_buffer'] = continuum_incoming_buffer()
        meta['slice_residuals']['rho_1.8'] = slice_state_residuals(
            energies, pg_from_ks_amplitudes(EXTENDED_RIGHT, energies, sampled[-1]),
            data['source'], data['projector'], EXTENDED_RIGHT)
        grids = {'coarse': (x_coarse, phi_coarse), 'fine': (x_fine, phi_fine)}
        times = arrays['times']
        for name, (x, phi) in grids.items():
            begin = time.process_time()
            owner = FourthOrderModePropagator(x, mass, angular)
            if meta['parameters'] is None:
                parameters = sigma_parameters(owner, loaded['multiplicity'])
                if abs(parameters['multiplicity'] - MULTIPLICITY) > 0:
                    raise ValueError('canonical multiplicity 12 required')
                if parameters['mass'] != mass or parameters['angular'] != angular:
                    raise ValueError('same inherited signed family required')
                meta['parameters'] = parameters
            _, spins, flux = initial_stationary_defect(owner, phi, energies)
            meta['field_runs_attempted'] += 1
            result = exact_harmonic_reference(owner, phi, energies, times)
            meta['case_CPU_seconds'][name] = time.process_time() - begin
            meta['cases'][name] = {'flux_algebra': result['flux_algebra'],
                                   'stationary_defect': spins,
                                   'initial_flux_algebra': flux}
            arrays.update({
                name + '/columns': result['columns'],
                name + '/pde_axial_columns': result['pde_axial_columns'],
                name + '/reference_columns': result['reference_columns'],
                name + '/reference_axial_columns': result['reference_axial_columns'],
                name + '/phase': result['phase'],
                name + '/expected_phase': result['expected_phase'],
                name + '/restriction': result['restriction'],
                name + '/sample_rows': result['sample_rows'],
            })
            meta['completed'].append(name)
            checkpoint('completed ' + name)
            print(name + ' exact-phase CPU ' + str(meta['case_CPU_seconds'][name]), flush=True)
        meta['stop_reason'] = 'completed exactly two prescribed reference buffer runs'
    except CPUExceeded:
        meta['cpu_budget_exceeded'] = True
        meta['stop_reason'] = 'OPEN: CPU cap reached; no automatic extra field case'
    finally:
        signal.setitimer(signal.ITIMER_PROF, 0)
        signal.signal(signal.SIGPROF, previous)
    return checkpoint(meta['stop_reason'])


def check():
    saved = json.loads((ROOT / OUTPUT).read_text())
    loaded = inputs()
    if digest(ARTIFACT) != saved['payload']['sha256']:
        raise ValueError('boundary-buffer artifact changed')
    with np.load(ROOT / ARTIFACT, allow_pickle=False) as archive:
        arrays = {key: archive[key] for key in archive.files}
    meta = json.loads(arrays['metadata_json'].tobytes())
    if meta['signature'] != signatures():
        raise ValueError('boundary-buffer owners/inputs changed')
    if meta.get('completed') and (not np.array_equal(arrays['original/x'], loaded['x']) or not np.array_equal(arrays['original/phi'], loaded['phi'])):
        raise ValueError('saved 3201 prefix bytes changed')
    if meta.get('completed') and not np.array_equal(arrays['source/source'], loaded['data']['source']):
        raise ValueError('original source fibers changed')
    if meta.get('completed') and not np.array_equal(arrays['source/energies'], loaded['data']['energies']):
        raise ValueError('original signed energies changed')
    actual = record(arrays, meta, saved['payload']['path'], saved['payload']['sha256'],
                    (ROOT / ARTIFACT).stat().st_size)
    if saved != actual:
        raise ValueError('incoming-boundary-buffer replay differs')
    return actual


def main():
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--run', action='store_true')
    group.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = run() if args.run else check()
    print(json.dumps({
        'status': result['status'],
        'runtime': result['runtime'],
        'measurements': {
            'prefix': result['measurements']['prefix'],
            'per_case': {name: {key: value for key, value in row.items()
                                if key not in ('PDE_values', 'reference_values')}
                         for name, row in result['measurements']['per_case'].items()},
            'maximum_PDE_matter_drift': result['measurements']['maximum_PDE_matter_drift'],
            'numerical_partial_control_passed': result['measurements']['numerical_partial_control_passed'],
        },
    }, indent=2))


if __name__ == '__main__':
    main()
