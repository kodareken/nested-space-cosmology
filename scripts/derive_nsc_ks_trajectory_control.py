#!/usr/bin/env python3
"""One actual-source trajectory capture, at most 30 CPU seconds; no sweep."""
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
import scipy

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
sys.path.insert(0, str(ROOT / 'scripts'))
import derive_nsc_ks_difference_control as C
from recursive_horizons import nsc_ks_source_envelope as E
from recursive_horizons.nsc_ks_difference_envelope import KSDifferenceIncoming
from recursive_horizons.nsc_ks_trajectory import capture_ks_trajectory, TrajectorySegment, TrajectoryResidualSampler
from recursive_horizons.nsc_evolved_incoming_state import FixedSourcePreparation
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes

OUTPUT = 'results/development/nsc-ks-trajectory-control.json'
PAYLOAD = 'results/development/artifacts/nsc-ks-trajectory-control.npz'
OWNED = ('scripts/derive_nsc_ks_trajectory_control.py', 'docs/nsc-ks-trajectory-control.md',
         'src/recursive_horizons/nsc_ks_trajectory.py', 'docs/nsc-ks-trajectory.md')
INPUTS = tuple(dict.fromkeys((*C.OWNED, *C.INPUTS, C.OUTPUT)))
POINTS, CAP, FRACTIONS = 64, 30., (.25, .5, .75)
FIELDS = ('columns', 'axial_columns', 'column_tangents', 'axial_tangents',
          'reference_amplitudes', 'envelope_difference', 'envelope_difference_z')


def digest(path):
    return sha256((ROOT / path).read_bytes()).hexdigest()


def signatures():
    return {p: digest(p) for p in (*OWNED, *INPUTS)}


def unpack(arrays, meta):
    family = C.R.C.P.family(C.R.C.P.ALPHA)
    mesh = E.computational_z_grid(POINTS)
    w, U = E._sample_axial_profiles(family.directions, mesh)
    binding = E.KSEnvelopeBinding(mesh, w, U, family.amplitudes,
        tuple(d.inner_radius for d in family.directions), tuple(d.outer_radius for d in family.directions),
        E.usual_axial_support(), meta['rho_up'], 1., **meta['options'])
    if binding.fingerprint != meta['binding']:
        raise ValueError('saved trajectory history binding changed')
    get = lambda key: arrays['prepared/' + key]
    prepared = KSDifferenceIncoming(arrays['target_z'], get('columns'), get('column_tangents'),
        get('axial_columns'), get('axial_tangents'), arrays['source/covariance'],
        arrays['source/column_weights'], arrays['source/energies'], meta['parameters']['mass'],
        meta['parameters']['angular'], arrays['initial_canonical_columns'], meta['rho_up'],
        binding, meta['preparation'], None, {}, get('reference_amplitudes'),
        get('envelope_difference'), get('envelope_difference_z'))
    nodes, coefficients, rho = arrays['state_nodes'], arrays['dense_corrections'], arrays['rho_nodes']
    segments = tuple(TrajectorySegment(rho[j], rho[j+1], nodes[j], nodes[j+1], coefficients[j])
                     for j in range(len(coefficients)))
    return prepared, family, segments


def analyze(arrays, meta):
    prepared, family, segments = unpack(arrays, meta)
    sampler = TrajectoryResidualSampler(prepared, family, oversampling=2)
    norms = np.array([[sampler.sample(part, x)['weighted_L2_indicators'] for x in FRACTIONS]
                      for part in segments])
    endpoint = sampler._envelope(segments[-1].end)
    Dtarget = np.moveaxis(E.trigonometric_polynomial(endpoint[1], sampler.native, prepared.z), -1, 0)
    endpoint_error = max(float(np.max(abs(endpoint[0] - prepared.reference_amplitudes))),
                         float(np.max(abs(Dtarget - prepared.envelope_difference))))
    widths = -np.diff(arrays['rho_nodes'])
    return norms, {
        'maximum_sampled_residual_norms_order_0_1_2': norms.max(axis=(0, 1)).tolist(),
        'sampled_rectangle_integral_indicators_not_bounds': (widths[:, None] * norms.max(axis=1)).sum(axis=0).tolist(),
        'endpoint_reconstruction_residual': endpoint_error,
        'endpoint_reconstruction_tolerance': 3e-14,
        'segment_count': len(segments), 'continuous_residual_integral_bounds': None,
        'physical_local_gate': 'OPEN',
    }


def record(arrays, meta, payload, measurements):
    return {'schema': 'NSC-KS-TRAJECTORY-CONTROL-v1', 'accountable_author': 'Douglas Ek',
        'status': 'OPEN: actual trajectory captured; continuous residual enclosures missing',
        'measurements': measurements,
        'control': {'maximum_field_runs': 1, 'envelope_points': POINTS, 'CPU_cap': CAP,
                    'source_family': '14_1', 'history_amplitude': C.R.C.P.ALPHA,
                    'local_interval': 'S(1)+[.12,.18]', 'options': meta['options']},
        'runtime': meta['runtime'], 'source_hashes': {p: meta['signature'][p] for p in OWNED},
        'input_hashes': {p: meta['signature'][p] for p in INPUTS}, 'reused_artifact': meta['reused_artifact'],
        'scope': {'sampled_norms_are_error_bounds': False, 'physical_local_gate': 'OPEN',
                  'source_changed': False, 'stress_drift_subtracted': False,
                  'source_error_bound': None, 'metric_timestep': False,
                  'value_only_capture': True, 'retarded_derivative_owner_unchanged': True},
        'payload': payload, 'reproducer': 'python3 scripts/derive_nsc_ks_trajectory_control.py --check'}


def run():
    if (ROOT / OUTPUT).exists() or (ROOT / PAYLOAD).exists():
        raise FileExistsError('trajectory capture exists; use --check')
    sig = signatures()
    old = C.check()
    with np.load(ROOT / old['payload']['path'], allow_pickle=False) as loaded:
        original = {k: loaded[k] for k in loaded.files}
    previous = json.loads(original['metadata_json'].tobytes())
    arrays = {k: v for k, v in original.items() if k.startswith('source/') or k in ('initial_canonical_columns', 'target_z')}
    source = FixedSourcePreparation(arrays['source/covariance'], arrays['source/column_weights'], arrays['source/energies'])
    meta = {'signature': sig, 'rho_up': previous['rho_up'], 'parameters': previous['parameters'],
            'options': C.R.C.TIGHT_OPTIONS, 'reused_artifact': old['payload'],
            'numpy_version': np.__version__, 'scipy_version': scipy.__version__}
    cpu, wall = time.process_time(), time.monotonic()
    def stop(*_):
        raise RuntimeError('one trajectory capture exceeded the 30 CPU-second cap')
    prior = signal.signal(signal.SIGPROF, stop)
    signal.setitimer(signal.ITIMER_PROF, CAP)
    try:
        prepared, segments = capture_ks_trajectory(source, arrays['initial_canonical_columns'],
            C.R.C.P.family(C.R.C.P.ALPHA), E.computational_z_grid(POINTS), arrays['target_z'],
            meta['parameters']['mass'], meta['parameters']['angular'], meta['rho_up'],
            axial_support=E.usual_axial_support(), tangents='zero', **meta['options'])
        if prepared.fixed_preparation_digest != previous['cases']['tighter512']['preparation']:
            raise ValueError('actual upstream preparation changed')
        meta.update(preparation=prepared.fixed_preparation_digest, binding=prepared.binding.fingerprint)
        arrays.update({'prepared/' + k: getattr(prepared, k) for k in FIELDS})
        arrays['rho_nodes'] = np.array([segments[0].rho_start] + [p.rho_end for p in segments])
        arrays['state_nodes'] = np.stack([segments[0].start] + [p.end for p in segments])
        arrays['dense_corrections'] = np.stack([p.coefficients for p in segments])
        norms, measurements = analyze(arrays, meta)
        arrays['sampled_residual_norms'] = norms
    finally:
        signal.setitimer(signal.ITIMER_PROF, 0)
        signal.signal(signal.SIGPROF, prior)
    meta['runtime'] = {'CPU_seconds': time.process_time() - cpu, 'wall_seconds': time.monotonic() - wall,
                       'field_runs': 1, 'CPU_cap': CAP}
    if signatures() != sig:
        raise ValueError('trajectory source/input changed')
    arrays['metadata_json'] = np.frombuffer(json.dumps(meta, sort_keys=True).encode(), np.uint8)
    raw = deterministic_npz_bytes(arrays)
    payload = {'path': PAYLOAD, 'sha256': sha256(raw).hexdigest(), 'bytes': len(raw)}
    result = record(arrays, meta, payload, measurements)
    (ROOT / PAYLOAD).write_bytes(raw)
    (ROOT / OUTPUT).write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + '\n')
    return result


def check():
    saved = json.loads((ROOT / OUTPUT).read_text())
    C.check()
    if digest(PAYLOAD) != saved['payload']['sha256']:
        raise ValueError('trajectory payload changed')
    with np.load(ROOT / PAYLOAD, allow_pickle=False) as loaded:
        arrays = {k: loaded[k] for k in loaded.files}
    meta = json.loads(arrays['metadata_json'].tobytes())
    if signatures() != meta['signature']:
        raise ValueError('trajectory owner/input changed')
    norms, measurements = analyze(arrays, meta)
    if not np.array_equal(norms, arrays['sampled_residual_norms']):
        raise ValueError('saved continuous-polynomial samples no longer replay')
    if record(arrays, meta, saved['payload'], measurements) != saved:
        raise ValueError('trajectory record no longer replays')
    return saved


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--run', action='store_true')
    group.add_argument('--check', action='store_true')
    result = run() if parser.parse_args().run else check()
    print(json.dumps({k: result[k] for k in ('status', 'measurements', 'runtime')}, indent=2))
