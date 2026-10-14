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
from recursive_horizons.nsc_posix import install as _install_posix
_install_posix()
from recursive_horizons.nsc_ks_source_envelope import resolution_control as R
from recursive_horizons.nsc_evolved_incoming_state import FixedSourcePreparation
from recursive_horizons.nsc_ks_difference_envelope import evolve_ks_difference_envelope
from recursive_horizons.nsc_ks_source_envelope import computational_z_grid, usual_axial_support
from recursive_horizons.nsc_ks_source_envelope import _require_envelope_capture_signature
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes

OUTPUT = 'results/development/nsc-ks-difference-control.json'
ARTIFACT = 'results/development/artifacts/nsc-ks-difference-control.npz'
ORIGINAL_SOURCE_CONTEXT = '329ef809e466a86ed4a290ee32c19485f85add9d'
ORIGINAL_OWNED = ('scripts/derive_nsc_ks_difference_control.py', 'docs/nsc-ks-difference-control.md',
         'src/recursive_horizons/nsc_ks_difference_envelope.py', 'docs/nsc-ks-difference-envelope.md')
ORIGINAL_INPUTS = tuple(dict.fromkeys((*R.OWNED, *R.INPUTS, R.OUTPUT)))
OWNED = ('scripts/derive_nsc_ks_difference_control.py', 'docs/nsc-ks-source-envelope.md',
         'src/recursive_horizons/nsc_ks_difference_envelope.py',
         'src/recursive_horizons/nsc_ks_source_envelope.py',
         'src/recursive_horizons/nsc_ks_matter_receipt_replay.py')
_CURRENT_INPUT_OWNERS = {
    'scripts/derive_nsc_ks_envelope_resolution.py':'src/recursive_horizons/nsc_ks_source_envelope.py',
    'scripts/derive_nsc_ks_envelope_comparison.py':'src/recursive_horizons/nsc_ks_source_envelope.py',
    'docs/nsc-ks-envelope-resolution.md':'docs/nsc-ks-source-envelope.md',
    'docs/nsc-ks-envelope-comparison.md':'docs/nsc-ks-source-envelope.md',
    'scripts/derive_nsc_evolved_incoming_constraints.py':'src/recursive_horizons/nsc_evolved_incoming_constraints.py',
    'scripts/derive_nsc_local_prepared_response.py':'scripts/derive_nsc_cosmological_atlas_v8.py',
    'docs/nsc-local-prepared-response.md':'docs/nsc-atlas-prepared-response.md',
}
INPUTS = tuple(dict.fromkeys(_CURRENT_INPUT_OWNERS.get(path,path) for path in ORIGINAL_INPUTS))
CASES = (('n256', 256, R.C.TIGHT_OPTIONS), ('n512', 512, R.C.TIGHT_OPTIONS),
         ('tighter512', 512, R.TIGHTER))
CAP = 30.
TARGET = 1e-11


from recursive_horizons.nsc_frozen_sources import (
    frozen_source_context, selected_source_context, selected_source_sha256)

class BudgetReached(RuntimeError):
    pass


def digest(path):
    from recursive_horizons.nsc_frozen_sources import selected_source_context, authenticate_source
    pins={
        'src/recursive_horizons/nsc_local_incoming_family.py':('09f3faa364ca88b8ff482c413022ce38f40bef0d','9d2b9fc57e2057d6c34d8943384428341b314370021bcd4156af88320c3de04b'),
        'src/recursive_horizons/nsc_ks_difference_envelope.py':('e62d28f3230ae6ddb2712df71eb7c61a8ed65ad6','77c89c8e720ca97336ccdee441a645b5c73472cb6f15b09e8120f68ab53db196'),
    }
    if selected_source_context() is not None and path in pins:
        context,expected=pins[path]
        manifest=json.loads((ROOT/'archive/source-history/manifest.json').read_text(encoding='utf-8'))
        if not any(p.get('path')==path and p.get('sha256')==expected and p.get('commit')==context for p in manifest['pins']):
            raise ValueError('original difference source pin changed')
        authenticate_source(ROOT,path,expected,context=context)
        return expected
    return selected_source_sha256(ROOT, path)


def signatures():
    owned,inputs = _binding_paths()
    return {p: digest(p) for p in (*owned, *inputs)}


def _binding_paths():
    if selected_source_context() is not None:
        return ORIGINAL_OWNED,ORIGINAL_INPUTS
    return OWNED,INPUTS


def replay_binding_receipt(saved,meta,arrays,name,binding,family):
    """Keep a stored sampling fingerprint separate from the current binding.

    The old profile sample arrays were not saved. Its fingerprint is bound
    by the exact local NPZ and original producer pins; we cannot reconstruct
    those bits from a different numerical runtime. Solver, domain recipe,
    coefficients, target and source arrays remain unchanged. This does not
    supply a field/geometry error certificate.
    """
    from recursive_horizons.nsc_frozen_sources import selected_source_context, source_sha256
    if selected_source_context() is None:
        raise ValueError('explicit original source context required for binding replay')
    if digest(ARTIFACT)!=saved['payload']['sha256']:
        raise ValueError('original binding receipt bytes or sources changed')
    _require_envelope_capture_signature(ROOT,meta['signature'],signatures())
    with np.load(ROOT/ARTIFACT,allow_pickle=False) as original:
        if set(arrays)!=set(original.files) or any(
                arrays[k].dtype!=original[k].dtype or not np.array_equal(arrays[k],original[k])
                for k in arrays):
            raise ValueError('original target, source or field arrays changed')
    case=next(c for c in saved['control']['cases'] if c[0]==name)
    _name,points,options=case
    expected=next(c for c in CASES if c[0]==name)
    if list(expected)!=case or any(getattr(binding,k)!=v for k,v in options.items()):
        raise ValueError('original binding solver recipe changed')
    if (binding.rho_up!=meta['rho_up'] or binding.rho_sigma!=1.
            or binding.amplitudes!=(saved['control']['history_amplitude'],)
            or binding.inner_radii!=(.007,) or binding.outer_radii!=(.03,)
            or not np.array_equal(binding.computational_z,computational_z_grid(points))
            or binding.axial_support!=usual_axial_support()
            or family.directions[0].w.coefficients!=(0.,1.)
            or family.directions[0].U.coefficients!=(0.,)
            or family.directions[0].w.inner!=.03 or family.directions[0].w.outer!=.06
            or arrays['target_z'].shape!=(21,) or not np.all(np.diff(arrays['target_z'])>0)):
        raise ValueError('original coefficient, domain or target recipe changed')
    original=meta['cases'][name]['binding']
    if not isinstance(original,str) or len(original)!=64 or any(c not in '0123456789abcdef' for c in original):
        raise ValueError('original binding fingerprint malformed')
    return {'original_fingerprint':original,'current_fingerprint':binding.fingerprint,
            'original_samples_retained':False,'original_source_context':selected_source_context(),
            'current_source_hashes':{p:source_sha256(ROOT,p) for p in (
                'scripts/derive_nsc_ks_difference_control.py',
                'src/recursive_horizons/nsc_ks_difference_envelope.py',
                'src/recursive_horizons/nsc_local_incoming_family.py')},
            'physical_field_error_bound':None}


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
    owned,inputs = _binding_paths()
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
        'source_hashes': {p: meta['signature'][p] for p in owned},
        'input_hashes': {p: meta['signature'][p] for p in inputs},
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
    with frozen_source_context(ORIGINAL_SOURCE_CONTEXT):
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
    _require_envelope_capture_signature(ROOT,meta['signature'],signatures())
    rebuilt = json.loads(json.dumps(record(arrays, meta, saved['payload'])))
    from recursive_horizons.nsc_ks_matter_receipt_replay import require_envelope_record
    require_envelope_record(saved,rebuilt,arrays,meta)
    return saved


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--run', action='store_true')
    group.add_argument('--check', action='store_true')
    parser.add_argument('--source-context', help='Full original source commit for saved --check; numerical inputs stay local.')
    args=parser.parse_args()
    if args.source_context is not None:
        if args.run:
            parser.error('--source-context selects historical --check only')
        with frozen_source_context(args.source_context):
            result=check()
    else:
        result=run() if args.run else check()
    print(json.dumps({k: result[k] for k in ('status', 'measurements', 'runtime')}, indent=2))
