#!/usr/bin/env python3
"""Bind resolved global modes and Gaussian variations to a common KS trace."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT/'src'))
from recursive_horizons.nsc_common_ks_trace import incoming_KS_state, restrict_resolved_modes, reference_KS_vertices
from recursive_horizons.nsc_transmitting_dirac_domain import MODE_TO_CURRENT
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes

OUTPUT = 'results/development/nsc-common-ks-trace.json'
MODE = 'results/development/nsc-pg-ctp-mode-jets.json'
HISTORY = 'results/development/nsc-transmitting-history-modes.json'
VERTEX = 'results/development/nsc-ks-spacetime-variation.json'
SOURCES = ('src/recursive_horizons/nsc_common_ks_trace.py', 'scripts/derive_nsc_common_ks_trace.py',
           'tests/test_nsc_common_ks_trace.py', 'docs/nsc-common-ks-trace.md')
INPUTS = (MODE, HISTORY, VERTEX, 'src/recursive_horizons/nsc_transmitting_dirac_domain.py',
          'src/recursive_horizons/nsc_ks_spacetime_variation.py', 'src/recursive_horizons/nsc_lorentzian.py',
          'src/recursive_horizons/nsc_paired_horizon_preparation.py', 'src/recursive_horizons/nsc_pg_massive_modes.py',
          'docs/nsc-causal-common-functional.md', 'docs/nsc-compact-ctp-completion.md',
          'docs/nsc-declared-action-scope.md')


def read(p): return json.loads((ROOT/p).read_text())
def sha(p): return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def load(p):
    with np.load(p, allow_pickle=False) as a: return {k: a[k].copy() for k in a.files}
def payload(r):
    p = r['payload']
    if sha(p['path']) != p['sha256']: raise ValueError('input mode artifact changed')
    return load(ROOT/p['path'])
def family(a, name): return {k[len(name)+1:]: v for k, v in a.items() if k.startswith(name+'/')}


def prepare():
    modes, history, old = (payload(read(p)) for p in (MODE, HISTORY, VERTEX))
    names = sorted({k.split('/')[0] for k in modes if '/' in k})
    out = {}
    for name in names:
        d = family(modes, name); E = d['energies']; x = history[name+'/x']
        ix = int(np.argmin(abs(x-1.)))
        if abs(x[ix]-1.) > 3e-12: raise ValueError('authenticated source grid must contain rho=1')
        fields = history[name+'/reference_field'].reshape(2, len(x), len(E), 3).transpose(1, 2, 0, 3)
        current = np.einsum('ab,ebs->eas', MODE_TO_CURRENT, fields[ix])
        state = incoming_KS_state(E, current, d['source'], d['projector'], rho=float(x[ix]))
        reference, _ = restrict_resolved_modes(0., E, d['at_zero'])
        at_zero = reference@d['source']@reference.swapaxes(-1, -2).conj()
        vertex = reference_KS_vertices(d)
        diag = dict(state['residuals'])
        diag['PG_KS_full_vertex'] = float(np.max(abs(vertex-old[name+'/vertex'])))
        diag['covariance_spectrum_transport'] = float(np.max(abs(np.linalg.eigvalsh(at_zero)-state['eigenvalues'])))
        diag['state_change_from_unit_radius'] = float(np.max(abs(at_zero-state['covariance'])))
        diag['source_complement_correlation_norm'] = float(np.linalg.norm(state['source_complement_correlations']))
        for key in ('canonical_momenta', 'mode_map', 'covariance', 'source_complement_correlations', 'source_complement', 'eigenvalues'):
            out[name+'/'+key] = state[key]
        out[name+'/KS_vertex'] = vertex
        out[name+'/diagnostics_json'] = np.frombuffer(json.dumps(diag, sort_keys=True).encode(), np.uint8)
        out[name+'/chart_json'] = np.frombuffer(json.dumps(state['chart'], sort_keys=True).encode(), np.uint8)
    metadata = {'signature': {p: sha(p) for p in (*SOURCES, *INPUTS)},
                'input_payloads': {p: read(p)['payload'] for p in (MODE, HISTORY, VERTEX)}}
    out['metadata_json'] = np.frombuffer(json.dumps(metadata, sort_keys=True).encode(), np.uint8)
    raw = deterministic_npz_bytes(out); digest = hashlib.sha256(raw).hexdigest()
    path = ROOT/f'results/development/artifacts/nsc-common-ks-trace.{digest}.npz'
    path.write_bytes(raw); return path


def make_record(path):
    a = load(path); metadata = json.loads(a['metadata_json'].tobytes())
    for p, expected in metadata['signature'].items():
        if sha(p) != expected: raise ValueError('common trace dependency changed: '+p)
    names = sorted({k.split('/')[0] for k in a if '/' in k})
    rows, maxima = {}, {}
    checks = ('source_coisometry', 'closed_source_projection', 'source_complement_projector',
              'Hermiticity', 'CAR_lower', 'CAR_upper', 'PG_KS_full_vertex', 'covariance_spectrum_transport')
    for name in names:
        d = json.loads(a[name+'/diagnostics_json'].tobytes())
        rows[name] = {'diagnostics': d, 'momenta': a[name+'/canonical_momenta'].tolist(),
                      'covariance_eigenvalues': a[name+'/eigenvalues'].tolist()}
        for key in checks: maxima[key] = max(maxima.get(key, 0.), d[key])
    failures = {k: v for k, v in maxima.items() if v > 3e-11}
    return {
        'schema': 'NSC-COMMON-KS-TRACE-v1', 'accountable_author': 'Douglas Ek',
        'status': 'PASS: incoming interior KS Cauchy restriction and common reference trace' if not failures else 'OPEN',
        'source_hashes': {p: sha(p) for p in SOURCES}, 'input_hashes': {p: sha(p) for p in INPUTS},
        'payload': {'path': str(path.relative_to(ROOT)), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest(), 'bytes': path.stat().st_size},
        'locked_inputs': read(MODE)['locked_inputs'], 'maxima': maxima, 'tolerance': 3e-11, 'failures': failures,
        'per_family': rows,
        'domain': {'incoming_slice': {'rho': 1., **json.loads(a[names[0]+'/chart_json'].tobytes())},
            'region': 'future interior trapped collar, unwrapped z line; both radial characteristics point downstream',
            'spectrum': 'same global source fibers and their inherited derivative-control frequencies; continuum not replaced by a finite closed system',
            'state': 'restriction of resolved global modes with C_H and inherited incoming occupation, not copied seed covariance',
            'Fourier_normalization': 'k=-E, dk/(2*pi)=dE/(2*pi)',
            'spacetime_measure': 'dT dz=d tau d rho/a0',
            'weak_vertex': 'all four raw KS directions, midpoint k=-(Ei+Ej)/2, phase and canonical-normal factors retained',
            'covariance_scope': 'physical Fourier-symbol samples; source-complement correlations are not an exterior joint-state reconstruction',
            'old_generators_rerun': False, 'source_frames': 'same-event frame applied to global solution restriction; not U_L0'},
        'selected_action': {'inventory': 'Gamma_G[C0]+Gamma_sub,ad4+local locked branch difference',
            'full_raw_heat_equivalence': 'separate UV interpretation, not a required execution gate',
            'finite_trace_requirement': 'common regulated pairing with cutoff/boundary terms and convergent spectral tails remains required',
            'new_prescription_or_coefficient': False},
        'gate': {'incoming_KS_state_restriction': 'PASS' if not failures else 'OPEN',
            'reference_geometry_common_trace': 'PASS' if not failures else 'OPEN',
            'nonlinear_combined_source': 'OPEN: common evolution/subtraction trace and spectral convergence required',
            'physical_EndpointBranchJets': 'OPEN', 'B2': 'OPEN', 'V_c': None, 'extended_stationarity': 'OPEN',
            'absolute_stress': None, 'nulls': None, 'updated_constraints': None, 'metric_timestep': False,
            'history_selected': False, 'Gamma_rest_extra_independent_term': 'none introduced', 'Z3': 'OUT OF SCOPE',
            'homogeneous_nonexistence': 'preserved in its original scope, not rerun'},
        'comparison': {'fields': 'all', 'float_atol': 3e-12, 'float_rtol': 3e-10, 'exact': 'hashes, domains, structure'},
        'reproducer': 'python3 scripts/derive_nsc_common_ks_trace.py --check',
    }


def main():
    p = argparse.ArgumentParser(); p.add_argument('--prepare', action='store_true'); p.add_argument('--check', action='store_true')
    args = p.parse_args()
    if args.prepare == args.check: p.error('choose prepare or check')
    path = prepare() if args.prepare else ROOT/read(OUTPUT)['payload']['path']
    result = make_record(path)
    if args.check:
        from derive_nsc_transmitting_boundary_binding import compare
        compare(read(OUTPUT), result)
    else: (ROOT/OUTPUT).write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    print(json.dumps({k: result[k] for k in ('status', 'maxima', 'failures')}, indent=2), flush=True)
    if result['failures']: raise SystemExit(1)


if __name__ == '__main__': main()
