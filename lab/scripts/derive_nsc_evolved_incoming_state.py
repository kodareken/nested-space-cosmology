#!/usr/bin/env python3
"""Bounded source-fixed incoming-state control; replay never propagates fields."""
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
import derive_nsc_retarded_compatible_response as RETARDED
from recursive_horizons.nsc_compatible_history_geometry import CompatibleRadiusDirection, CompatibleIncomingMetric
from recursive_horizons.nsc_evolved_incoming_state import (
    AmplitudeOnlyMetric, CachedCompatibleIncomingMetric, FixedSourcePreparation,
    evolve_incoming_state, restriction_covariance_tangent, rho1_restriction_map,
    rho1_node_index,
)
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes
from recursive_horizons.nsc_transmitting_history_jets import FourthOrderModePropagator
import recursive_horizons.nsc_evolved_incoming_state as STATE

OUTPUT = 'results/development/nsc-evolved-incoming-state.json'
ARTIFACT = 'results/development/artifacts/nsc-evolved-incoming-state.npz'
OWNED = (
    'src/recursive_horizons/nsc_evolved_incoming_state.py',
    'scripts/derive_nsc_evolved_incoming_state.py',
    'tests/test_nsc_evolved_incoming_state.py',
    'docs/nsc-evolved-incoming-state.md',
)
RECORDS = (
    'results/development/nsc-retarded-compatible-response.json',
    'results/development/nsc-transmitting-history-modes.json',
    'results/development/nsc-pg-ctp-mode-jets.json',
    'results/development/nsc-compatible-prepared-history.json',
)
OWNERS = tuple(
    'src/recursive_horizons/' + name + '.py' for name in (
        'nsc_prepared_history_jets', 'nsc_compatible_history_geometry',
        'nsc_transmitting_history_jets', 'nsc_transmitting_history_modes',
        'nsc_transmitting_dirac_domain', 'nsc_ks_spacetime_variation',
        'nsc_pg_ks_metric_pullback', 'nsc_transmitting_resolvent',
        'nsc_lorentzian', 'nsc_mode_resolved_cauchy_state',
    )
) + (
    'scripts/derive_nsc_retarded_compatible_response.py',
    'docs/nsc-retarded-compatible-preparation.md',
    'docs/nsc-preparation-domain-decision.md',
)
AMPLITUDE = 0.001
STEP = 1e-4
POINTS, STEPS = 801, 64
CPU_CAP = 120.
FIELD_TOL = 3e-8
ALGEBRA_TOL = 3e-11
MAXABS = lambda x: float(np.max(np.abs(x)))


def digest(path):
    return sha256((ROOT / path).read_bytes()).hexdigest()


def signature():
    return {p: digest(p) for p in (*OWNED, *RECORDS, *OWNERS)}


def load_retarded_inputs():
    x, phi, data, records = RETARDED.inputs()
    guard = RETARDED.support_guard()
    retarded = json.loads((ROOT / RECORDS[0]).read_text())
    if digest(retarded['payload']['path']) != retarded['payload']['sha256']:
        raise ValueError('retarded response artifact changed')
    return x, phi, data, records, guard, retarded


def family(amplitude, guard):
    direction = CompatibleRadiusDirection(
        RETARDED.axial_bump(guard['axial_center']), lambda z, n: 0., .007, .03)
    return CompatibleIncomingMetric((float(amplitude),), (direction,))


def incident(owner, phi, data, ndir):
    n = len(owner.x)
    energies = np.repeat(data['energies'], 3)
    tangents = np.zeros((ndir, *phi.shape), complex)
    right = np.array([n - 1, 2 * n - 1])
    N, beta, q, _ = owner.reference
    speed = np.array([N[-1] / q[-1] - beta[-1], -N[-1] / q[-1] - beta[-1]])
    sat = -speed / owner.weights[-1]

    def incoming(time):
        forcing = np.zeros_like(phi)
        forcing[right] = sat[:, None] * phi[right] * np.exp(-1j * energies * time)[None, :]
        return forcing, np.zeros_like(tangents)

    return phi, tangents, incoming


def kernel(columns, weights, covariance):
    nsrc = covariance.shape[0]
    F = (np.asarray(columns) * np.asarray(weights)).reshape(-1, nsrc)
    return F @ covariance @ F.conj().T


def probe_cache(guard):
    lower, upper = guard['radial_support']
    x = np.array([-2., lower, 0.99, 1., 1.01, upper, 1.2])
    metric = family(AMPLITUDE, guard)
    cache = CachedCompatibleIncomingMetric(x, metric)
    times = (0., 0.11, 0.15, 0.21, 0.3)
    result = cache.verify_against_owned(times)
    if result['maximum_radius_change_from_background'] <= 0:
        raise ArithmeticError('nonzero-amplitude cache did not change the actual radius')
    zero = CompatibleIncomingMetric((0.,), metric.directions)
    if MAXABS(cache.values(0.15, x)[3] - zero.values(0.15, x)[3]) <= 0:
        raise ArithmeticError('cache agrees with a frozen zero-radius background at nonzero amplitude')
    return result


def retarded_restriction_replay(retarded, data):
    path = ROOT / retarded['payload']['path']
    with np.load(path, allow_pickle=False) as arrays:
        field = arrays['801_64/field']
        tangent = arrays['801_64/field_tangent']
        trace = arrays['801_64/trace']
        trace_tangent = arrays['801_64/trace_tangent']
        sample_rows = arrays['801_64/sample_rows']
        x = arrays['801_64/x']
        saved_kernel = arrays['801_64/covariance_tangent']
        source = arrays['source/source']
        weights = arrays['source/weights']
    owner = FourthOrderModePropagator(x, data['mass'].item(), data['angular'].item())
    index = rho1_node_index(x)
    N, beta, q, radius = owner.reference[:, index]
    restriction = rho1_restriction_map(N, q, beta, radius)
    field_trace = restriction @ field[sample_rows]
    tangent_trace = restriction @ tangent[sample_rows]
    root = np.repeat(np.sqrt(weights / (2 * np.pi)), 3)
    C = FixedSourcePreparation.from_signed_blocks(
        data['energies'], weights, source).covariance
    F = (trace * root).reshape(-1, 12)
    dF = (trace_tangent * root).reshape(-1, 12)
    rebuilt = restriction_covariance_tangent(F, dF, C, np.zeros_like(C))
    return {
        'final_trace_residual': MAXABS(field_trace - trace[-1]),
        'final_trace_tangent_residual': MAXABS(tangent_trace - trace_tangent[-1]),
        'sampled_kernel_residual': MAXABS(rebuilt - saved_kernel),
        'restriction_parameter_derivative': 0.,
        'old_solver_rerun': False,
    }


def one_state(amplitude, ndir, owner, phi, data, times, source, guard):
    metric = family(amplitude, guard)
    provider = CachedCompatibleIncomingMetric(owner.x, metric)
    if ndir == 0:
        provider = AmplitudeOnlyMetric(provider)
    fields, tangents, incoming = incident(owner, phi, data, ndir)
    return evolve_incoming_state(
        owner, times, provider, fields, tangents, incoming, source)


class CPUExceeded(RuntimeError):
    pass


def analyze(arrays, meta):
    weights = arrays['source/column_weights']
    C = arrays['source/covariance']
    field_fd = (arrays['plus/field'] - arrays['minus/field']) / (2 * STEP)
    trace_fd = (arrays['plus/columns'] - arrays['minus/columns']) / (2 * STEP)
    cov_fd = (kernel(arrays['plus/columns'], weights, C)
              - kernel(arrays['minus/columns'], weights, C)) / (2 * STEP)
    analytic_cov = arrays['center/covariance_tangent'].reshape(
        arrays['center/columns'].shape[0] * 2, arrays['center/columns'].shape[0] * 2)
    residuals = {
        'complete_field_derivative': MAXABS(field_fd - arrays['center/field_tangent']),
        'incoming_trace_derivative': MAXABS(trace_fd - arrays['center/column_tangents']),
        'complete_coherent_covariance_derivative': MAXABS(cov_fd - analytic_cov),
        'norm_flux_residual': float(meta['algebra']['norm_flux_residual']),
        'metric_flux_tangent_residual': float(meta['algebra']['metric_flux_tangent_residual']),
        'auxiliary_identity_residual': float(meta['algebra']['auxiliary_identity_residual']),
        'source_covariance_hermiticity': MAXABS(C - C.conj().T),
        'retarded_final_trace_residual': float(meta['retarded_restriction']['final_trace_residual']),
        'retarded_kernel_residual': float(meta['retarded_restriction']['sampled_kernel_residual']),
        'cache_value_residual': float(meta['cache']['value_residual']),
        'cache_direction_residual': float(meta['cache']['direction_residual']),
    }
    F = (arrays['center/columns'] * weights).reshape(-1, C.shape[0])
    dF = (arrays['center/column_tangents'][None, ...] * weights).reshape(1, F.shape[0], C.shape[0])
    frozen = restriction_covariance_tangent(F, np.zeros_like(dF), C, np.zeros((1, *C.shape), complex))
    fake = np.zeros((1, *C.shape), complex)
    fake[0, 0, 1] = 0.05
    fake[0, 1, 0] = 0.05
    wrong_source = restriction_covariance_tangent(F, dF, C, fake)
    omitted = {
        'frozen_incoming_columns': MAXABS(frozen[0] - analytic_cov),
        'nonzero_source_covariance_tangent': MAXABS(wrong_source[0] - analytic_cov),
        'finite_difference_versus_frozen_columns': MAXABS(cov_fd - frozen[0]),
    }
    algebra = max(residuals[k] for k in (
        'norm_flux_residual', 'metric_flux_tangent_residual', 'auxiliary_identity_residual'))
    passed = (
        residuals['complete_field_derivative'] < FIELD_TOL
        and residuals['incoming_trace_derivative'] < FIELD_TOL
        and residuals['complete_coherent_covariance_derivative'] < FIELD_TOL
        and algebra < ALGEBRA_TOL
        and residuals['source_covariance_hermiticity'] < ALGEBRA_TOL
        and residuals['retarded_final_trace_residual'] < 3e-12
        and residuals['retarded_kernel_residual'] < 3e-12
        and residuals['cache_value_residual'] < 3e-12
        and min(omitted.values()) > 1e-16
        and meta['field_runs'] == 3
        and not meta['cpu_budget_exceeded']
    )
    return residuals, omitted, passed, algebra


def record_from_arrays(arrays, meta, artifact_hash, artifact_bytes):
    residuals, omitted, passed, algebra = analyze(arrays, meta)
    current = meta['verification_signature']
    return {
        'schema': 'NSC-EVOLVED-INCOMING-STATE-v1',
        'accountable_author': 'Douglas Ek',
        'status': (
            'PASS: source-fixed adapter control; sampled source and continuum error OPEN'
            if passed else
            'OPEN: source-fixed adapter control residual or coverage gate not met'
        ),
        'control': {
            'family': '14_1',
            'amplitude': AMPLITUDE,
            'finite_difference_step': STEP,
            'grid_points': POINTS,
            'midpoint_steps': STEPS,
            'PG_time_window': [0., 0.3],
            'selected_physical_duration': None,
            'numerical_normal_window_radii': [0.007, 0.03],
            'field_runs': meta['field_runs'],
            'field_run_cap': 3,
            'initial_and_incident_tangents': 'explicit 0 by unchanged entire past/upstream support',
            'source_law': 'unchanged original coherent source blocks; explicit dC_src=0',
            'weights': 'sqrt(w/(2*pi)) repeated on each source column; no multiplicity factor added',
            'restriction': 'actual rho1 MODE_TO_CURRENT/normal half-density; dR=0 for this radius family',
            'stationary_phase_reapplied': False,
            'outgoing_momentum_substitution': None,
            'field_derivative_tolerance': FIELD_TOL,
            'covariance_derivative_tolerance': FIELD_TOL,
            'flux_algebra_tolerance': ALGEBRA_TOL,
            'support': meta['support'],
        },
        'residuals': residuals,
        'omitted_term_effects': omitted,
        'algebra_maximum': algebra,
        'cache_verification': meta['cache'],
        'retarded_array_only_restriction': meta['retarded_restriction'],
        'runtime': {k: meta[k] for k in (
            'CPU_seconds', 'wall_seconds', 'cpu_budget_seconds', 'cpu_budget_exceeded',
            'case_CPU_seconds', 'field_runs')},
        'scope': {
            'source_fixed_incoming_state': True,
            'frozen_incoming_C0': False,
            'sampled_source': 'OPEN',
            'continuum_error_bound': None,
            'finite_source_error': 'OPEN',
            'discretization_error': 'OPEN',
            'full_CAR': 'OPEN',
            'physical_history_solution': 'OPEN',
            'cosmology_EXISTENCE': 'OPEN',
            'metric_timestep': False,
            'new_stress_or_Gamma': False,
            'A_q_Omega_zeta_Vfull_refit': False,
            'old_baseline_scattering_horizon_or_source_grids_regenerated': False,
            'publication': False,
        },
        'source_hashes': {p: current[p] for p in OWNED},
        'input_hashes': {p: current[p] for p in (*RECORDS, *OWNERS)},
        'reused_artifacts': meta['reused_artifacts'],
        'payload': {'path': ARTIFACT, 'sha256': artifact_hash, 'bytes': artifact_bytes},
        'reproducer': 'python3 scripts/derive_nsc_evolved_incoming_state.py --check',
    }


def publish(arrays, meta):
    current = signature()
    meta['verification_signature'] = dict(current)
    meta['signature'] = current
    payload = {**arrays, 'metadata_json': np.frombuffer(
        json.dumps({k: meta[k] for k in meta if k != 'signature'}, sort_keys=True).encode(), np.uint8)}
    raw = deterministic_npz_bytes(payload)
    path = ROOT / ARTIFACT
    path.write_bytes(raw)
    record = record_from_arrays(payload, meta, sha256(raw).hexdigest(), len(raw))
    (ROOT / OUTPUT).write_text(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + '\n')
    return record


def run():
    start_cpu = time.process_time()
    start_wall = time.monotonic()
    before = signature()
    meta = {
        'completed': [], 'case_CPU_seconds': {}, 'cpu_budget_seconds': CPU_CAP,
        'cpu_budget_exceeded': False, 'field_runs': 0,
    }

    def stop(signum, frame):
        raise CPUExceeded('120 second producer CPU budget exhausted')

    old_handler = signal.signal(signal.SIGPROF, stop)
    signal.setitimer(signal.ITIMER_PROF, CPU_CAP)
    inner = STATE.evolve_prepared_field_jets
    field_runs = [0]

    def counted(*args, **kwargs):
        field_runs[0] += 1
        if field_runs[0] > 3:
            raise RuntimeError('more than three 801-grid field runs were requested')
        return inner(*args, **kwargs)

    STATE.evolve_prepared_field_jets = counted
    arrays = {}
    try:
        x, phi, data, records, guard, retarded = load_retarded_inputs()
        meta['support'] = guard
        meta['reused_artifacts'] = [r['payload'] for r in records] + [retarded['payload']]
        meta['cache'] = probe_cache(guard)
        meta['retarded_restriction'] = retarded_restriction_replay(retarded, data)
        if x.shape != (POINTS,) or phi.shape != (1602, 12):
            raise ValueError('authenticated group14 801-grid columns required')
        owner = FourthOrderModePropagator(x, data['mass'].item(), data['angular'].item())
        source = FixedSourcePreparation.from_signed_blocks(
            data['energies'], data['weights'], data['source'])
        times = np.linspace(0., 0.3, STEPS + 1)
        print('support and cache checked; three prescribed 801/64 solves only', flush=True)
        jobs = (('center', AMPLITUDE, 1), ('plus', AMPLITUDE + STEP, 0),
                ('minus', AMPLITUDE - STEP, 0))
        states = {}
        for name, amplitude, ndir in jobs:
            begin = time.process_time()
            states[name] = one_state(amplitude, ndir, owner, phi, data, times, source, guard)
            meta['case_CPU_seconds'][name] = time.process_time() - begin
            meta['completed'].append(name)
            print(name + ' checkpoint, CPU ' + str(meta['case_CPU_seconds'][name]), flush=True)
        meta['field_runs'] = field_runs[0]
        center, plus, minus = states['center'], states['plus'], states['minus']
        center.require_history(family(AMPLITUDE, guard), times, x)
        meta['algebra'] = {
            'norm_flux_residual': center.diagnostics['norm_flux_residual'],
            'metric_flux_tangent_residual': center.diagnostics['metric_flux_tangent_residual'],
            'auxiliary_identity_residual': center.diagnostics['auxiliary_identity_residual'],
        }
        arrays.update({
            'source/covariance': np.asarray(source.covariance),
            'source/column_weights': np.asarray(source.column_weights),
            'source/energies': np.asarray(source.energies),
            'source/blocks': np.asarray(data['source']),
            'source/weights': np.asarray(data['weights']),
            'source/mass': np.asarray(data['mass']),
            'source/angular': np.asarray(data['angular']),
            'center/z': np.asarray(center.z),
            'center/columns': np.asarray(center.columns),
            'center/column_tangents': np.asarray(center.column_tangents[0]),
            'center/covariance': np.asarray(center.covariance()),
            'center/covariance_tangent': np.asarray(center.covariance_tangent()[0]),
            'center/field': np.asarray(center.diagnostics['evolved_field']),
            'center/field_tangent': np.asarray(center.diagnostics['tangent_fields'][0]),
            'center/restriction': np.asarray(center.history.restriction),
            'center/rho1_metric': np.asarray(center.history.rho1_metric),
            'plus/columns': np.asarray(plus.columns),
            'plus/field': np.asarray(plus.diagnostics['evolved_field']),
            'minus/columns': np.asarray(minus.columns),
            'minus/field': np.asarray(minus.diagnostics['evolved_field']),
            'times': times,
            'x': np.asarray(x),
        })
        meta.update(CPU_seconds=time.process_time() - start_cpu,
                    wall_seconds=time.monotonic() - start_wall)
        if signature() != before:
            raise ValueError('input/producer changed during bounded adapter control')
        if field_runs[0] != 3:
            raise RuntimeError('source-fixed control must execute exactly three field runs')
        return publish(arrays, meta)
    except CPUExceeded:
        meta['cpu_budget_exceeded'] = True
        meta.update(CPU_seconds=time.process_time() - start_cpu,
                    wall_seconds=time.monotonic() - start_wall)
        raise
    finally:
        signal.setitimer(signal.ITIMER_PROF, 0)
        signal.signal(signal.SIGPROF, old_handler)
        STATE.evolve_prepared_field_jets = inner


def check():
    old = json.loads((ROOT / OUTPUT).read_text())
    if digest(old['payload']['path']) != old['payload']['sha256']:
        raise ValueError('evolved incoming-state artifact changed')
    path = ROOT / old['payload']['path']
    with np.load(path, allow_pickle=False) as loaded:
        arrays = {k: loaded[k] for k in loaded.files}
    meta = json.loads(arrays['metadata_json'].tobytes())
    current = signature()
    if meta['verification_signature'] != current:
        raise ValueError('evolved incoming-state producer/input changed')
    meta['verification_signature'] = current
    record = record_from_arrays(arrays, meta, digest(old['payload']['path']), path.stat().st_size)
    if record != old:
        raise ValueError('evolved incoming-state receipt differs from replay')
    return record


def main():
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--run', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    if args.run:
        record = run()
    else:
        record = check()
    print(json.dumps({k: record[k] for k in (
        'status', 'residuals', 'omitted_term_effects', 'runtime', 'scope')}, indent=2))


if __name__ == '__main__':
    main()
