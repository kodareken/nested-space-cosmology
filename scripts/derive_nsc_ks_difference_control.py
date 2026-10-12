#!/usr/bin/env python3
"""Three bounded actual-source checks of the joint reference/difference solve."""
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
import derive_nsc_ks_envelope_resolution as R
from recursive_horizons.nsc_evolved_incoming_state import FixedSourcePreparation
from recursive_horizons.nsc_ks_difference_envelope import evolve_ks_difference_envelope
from recursive_horizons.nsc_ks_source_envelope import computational_z_grid, usual_axial_support
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes

OUTPUT = 'results/development/nsc-ks-difference-control.json'
ARTIFACT = 'results/development/artifacts/nsc-ks-difference-control.npz'
OWNED = ('scripts/derive_nsc_ks_difference_control.py', 'docs/nsc-ks-difference-control.md',
         'src/recursive_horizons/nsc_ks_difference_envelope.py', 'docs/nsc-ks-difference-envelope.md')
INPUTS = tuple(dict.fromkeys((*R.OWNED, *R.INPUTS, R.OUTPUT)))
CASES = (('n256', 256, R.C.TIGHT_OPTIONS), ('n512', 512, R.C.TIGHT_OPTIONS),
         ('tighter512', 512, R.TIGHTER))
CAP = 30.
TARGET = 1e-11


class BudgetReached(RuntimeError):
    pass


def digest(path):
    return sha256((ROOT / path).read_bytes()).hexdigest()


def signatures():
    return {p: digest(p) for p in (*OWNED, *INPUTS)}


def analyze(arrays, meta):
    values = {n: R.C.matter(arrays, n, meta['parameters']) for n in ['previous512', *meta['completed']]}
    checks = {}
    for first, second, label in (('n256', 'n512', 'space'), ('n512', 'tighter512', 'time'),
                                 ('previous512', 'tighter512', 'representation')):
        if first in values and second in values:
            delta = np.max(abs(values[second]['action_gradient'] - values[first]['action_gradient']), axis=0)
            tangent = np.max(abs(values[second]['action_gradient_tangent'] -
                                 values[first]['action_gradient_tangent']), axis=(0, 1))
            checks[label] = {'matter_change': delta.tolist(), 'tangent_change': tangent.tolist(),
                             'target': TARGET, 'passes_indicator': bool(np.max(delta) <= TARGET)}
    return checks


def record(arrays, meta, payload):
    checks = analyze(arrays, meta)
    passed = all(k in checks and checks[k]['passes_indicator'] for k in ('space', 'time'))
    return {
        'schema': 'NSC-KS-DIFFERENCE-CONTROL-v1', 'accountable_author': 'Douglas Ek',
        'status': ('PASS' if passed else 'OPEN') + ': local numerical indicators; physical local gate OPEN',
        'control': {'source_family': '14_1', 'history_amplitude': R.C.P.ALPHA,
                    'local_interval': 'S(1)+[.12,.18]', 'max_runs': 3, 'cases': CASES,
                    'source_fixed_across_all_cases': True},
        'measurements': checks, 'runtime': meta['runtime'],
        'scope': {'physical_local_gate': 'OPEN', 'full_source_error': None,
                  'continuum_error': None, 'between_node_error': None,
                  'indicators_are_error_bounds': False, 'new_PG_runs': 0,
                  'stress_drift_subtracted': False, 'metric_timestep': False,
                  'same_Dirac_law': True, 'frozen_C0_promoted': False},
        'source_hashes': {p: meta['signature'][p] for p in OWNED},
        'input_hashes': {p: meta['signature'][p] for p in INPUTS},
        'reused_artifacts': meta['reused_artifacts'], 'payload': payload,
        'reproducer': 'python3 scripts/derive_nsc_ks_difference_control.py --check',
    }


def publish(arrays, meta):
    raw = deterministic_npz_bytes({**arrays, 'metadata_json': np.frombuffer(
        json.dumps(meta, sort_keys=True).encode(), np.uint8)})
    (ROOT / ARTIFACT).write_bytes(raw)
    result = record(arrays, meta, {'path': ARTIFACT, 'sha256': sha256(raw).hexdigest(), 'bytes': len(raw)})
    (ROOT / OUTPUT).write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + '\n')
    return result


def run():
    if (ROOT / OUTPUT).exists() or (ROOT / ARTIFACT).exists():
        raise FileExistsError('existing difference control; use --check')
    sig = signatures()
    previous = R.check()
    with np.load(ROOT / previous['payload']['path'], allow_pickle=False) as f:
        old = {k: f[k] for k in f.files}
    oldmeta = json.loads(old['metadata_json'].tobytes())
    arrays = {k: v for k, v in old.items() if k.startswith('source/') or k in ('target_z', 'initial_canonical_columns')}
    for k in ('F', 'Fz', 'dF', 'dFz'):
        arrays['previous512/' + k] = old['tighter512/' + k]
    source = FixedSourcePreparation(arrays['source/covariance'], arrays['source/column_weights'], arrays['source/energies'])
    meta = {'signature': sig, 'parameters': oldmeta['parameters'], 'rho_up': oldmeta['rho_up'],
            'completed': [], 'attempted': [], 'cases': {}, 'reused_artifacts': [previous['payload']]}
    cpu, wall, capped = time.process_time(), time.monotonic(), False

    def stop(*_):
        raise BudgetReached('three difference solves: 30 CPU-second cap')

    def checkpoint():
        meta['runtime'] = {'CPU_seconds': time.process_time() - cpu, 'wall_seconds': time.monotonic() - wall,
                           'CPU_cap': CAP, 'completed': list(meta['completed']),
                           'attempted': list(meta['attempted']), 'cap_reached': capped}
        if signatures() != sig:
            raise ValueError('difference control owner/input changed')
        return publish(arrays, meta)

    prior = signal.signal(signal.SIGPROF, stop)
    signal.setitimer(signal.ITIMER_PROF, CAP)
    try:
        for name, n, options in CASES:
            begin = time.process_time()
            meta['attempted'].append(name)
            result = evolve_ks_difference_envelope(source, arrays['initial_canonical_columns'],
                R.C.P.family(R.C.P.ALPHA), computational_z_grid(n), arrays['target_z'],
                meta['parameters']['mass'], meta['parameters']['angular'], meta['rho_up'],
                axial_support=usual_axial_support(), tangents='all', **options)
            for key, value in (('F', result.columns), ('Fz', result.axial_columns),
                               ('dF', result.column_tangents), ('dFz', result.axial_tangents),
                               ('A', result.reference_amplitudes), ('D', result.envelope_difference),
                               ('Dz', result.envelope_difference_z)):
                arrays[name + '/' + key] = value
            if result.fixed_preparation_digest != oldmeta['cases']['tighter512']['preparation']:
                raise ValueError('fixed preparation changed')
            meta['cases'][name] = {'CPU_seconds': time.process_time() - begin,
                'preparation': result.fixed_preparation_digest, 'binding': result.binding.fingerprint,
                'function_evaluations': result.diagnostics['function_evaluations']}
            meta['completed'].append(name)
            checkpoint()
    except BudgetReached:
        capped = True
    finally:
        signal.setitimer(signal.ITIMER_PROF, 0)
        signal.signal(signal.SIGPROF, prior)
    return checkpoint()


def check():
    saved = json.loads((ROOT / OUTPUT).read_text())
    R.check()
    if digest(ARTIFACT) != saved['payload']['sha256']:
        raise ValueError('difference payload changed')
    with np.load(ROOT / ARTIFACT, allow_pickle=False) as f:
        arrays = {k: f[k] for k in f.files}
    meta = json.loads(arrays['metadata_json'].tobytes())
    if signatures() != meta['signature']:
        raise ValueError('difference source/input changed')
    rebuilt = json.loads(json.dumps(record(arrays, meta, saved['payload'])))
    if rebuilt != saved:
        raise ValueError('difference array replay differs')
    return saved


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--run', action='store_true')
    group.add_argument('--check', action='store_true')
    result = run() if parser.parse_args().run else check()
    print(json.dumps({k: result[k] for k in ('status', 'measurements', 'runtime')}, indent=2))
