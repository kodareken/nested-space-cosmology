#!/usr/bin/env python3
"""Reference-band action and its explicit local variation decomposition."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT/'src'))
from recursive_horizons.nsc_spatial_reference_symbol import supplied_KS_metric_jets, reference_projector
from recursive_horizons.nsc_reference_band_action import (
    band_frame, direction_series, action_variation_identity, values,
    exact_massless_reference,
)
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes

OUTPUT = 'results/development/nsc-reference-band-action.json'
PREVIOUS = 'results/development/nsc-spatial-reference-symbol.json'
GENERIC_DIAGNOSTIC = 'results/development/artifacts/nsc-reference-band-action.98c734ac57cbb9834c2b446c0e95538caea7c112456b57f5c47cfd61d6945003.npz'
SOURCES = ('src/recursive_horizons/nsc_reference_band_action.py',
           'scripts/derive_nsc_reference_band_action.py', 'tests/test_nsc_reference_band_action.py',
           'docs/nsc-reference-band-action.md')
INPUTS = (PREVIOUS, GENERIC_DIAGNOSTIC, 'src/recursive_horizons/nsc_spatial_reference_symbol.py',
          'src/recursive_horizons/nsc_ks_spacetime_variation.py',
          'src/recursive_horizons/nsc_pg_ks_metric_pullback.py', 'src/recursive_horizons/nsc_lorentzian.py',
          'src/recursive_horizons/nsc_compact_matching.py', 'src/recursive_horizons/nsc_scale_binding.py',
          'docs/nsc-vacuum-matched-ctp.md', 'docs/nsc-adm-source-constraints.md',
          'docs/nsc-compact-ctp-completion.md', 'docs/nsc-compact-matching.md',
          'docs/nsc-declared-action-scope.md')
AMPLITUDES = np.array([.002, .001, .003, .001])
STEP = 1e-4


def sha(p): return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def read(p): return json.loads((ROOT/p).read_text())
def load(p):
    with np.load(p, allow_pickle=False) as f: return {k: f[k].copy() for k in f.files}


def variation(job):
    B, k, m, l, Pi, baseP, frame = job
    weights = np.array([1., -8., 8., -1.])/(12*STEP)
    rows = []
    for offset in (-2., -1., 1., 2.):
        amplitudes = AMPLITUDES+offset*STEP*np.eye(4)[B]
        P = reference_projector(*supplied_KS_metric_jets(.019, .23, amplitudes), k, m, l)
        P = exact_massless_reference(P, m, l, k)
        rows.append((P, band_frame(P, Pi)))
    dH = sum(w*p['hamiltonian'] for w, (p, f) in zip(weights, rows))
    dU = direction_series([f['U'] for p, f in rows], weights)
    dh = direction_series([f['effective'] for p, f in rows], weights)
    output = action_variation_identity(baseP, frame, dH, dU, dh)
    return B, output


def prepare(workers):
    old = read(PREVIOUS); p = old['payload']
    if sha(p['path']) != p['sha256']: raise ValueError('old reference symbol artifact changed')
    a = load(ROOT/p['path']); oldmeta = json.loads(a['metadata_json'].tobytes())
    k, m, l = a['momentum'], a['mass'], a['angular']
    # The full Taylor data are needed for the frame, beyond the old point values.
    baseP = reference_projector(*supplied_KS_metric_jets(.019, .23), k, m, l)
    previous = float(np.max(abs(baseP['orders']-a['spatial_orders'])))
    baseP = exact_massless_reference(baseP, m, l, k)
    Pi = np.tile(np.diag([1., 0.]), (len(k), 1, 1)).astype(complex)
    for j in np.flatnonzero((m == 0) & (l == 0) & (k > 0)): Pi[j] = np.diag([0., 1.])
    frame = band_frame(baseP, Pi)
    arrays = {'target': Pi, 'band_energy': frame['band_energy'],
              'effective_coefficients': values(frame['effective']),
              'connection_coefficients': values(frame['connection'])}
    for name in ('U', 'effective', 'connection'):
        for n, jet in enumerate(frame[name]): arrays[f'jets/{name}/{n}'] = jet.data
    for kind in ('raw_residuals', 'scaled_residuals'):
        for name, data in frame[kind].items(): arrays[kind+'/'+name] = data
    print('Reference frame and band action constructed on 126 inherited local symbols', flush=True)
    with ProcessPoolExecutor(max_workers=workers) as pool:
        jobs = [(B, k, m, l, Pi, baseP, frame) for B in range(4)]
        for B, result in pool.map(variation, jobs):
            arrays.update({f'variation/{B}/{name}': data for name, data in result.items()})
            print(f'Raw-KS action variation {B+1}/4 completed', flush=True)
    fields = supplied_KS_metric_jets(.019, .23)
    expected = -baseP['gap']-fields[1].value[0, 0, 0]*k
    metadata = {'signature': {p: sha(p) for p in (*SOURCES, *INPUTS)},
                'input_payload': old['payload'], 'labels': oldmeta['labels'],
                'previous_projector_residual': previous,
                'leading_overlap_identity': frame['leading_overlap_identity'],
                'LLL_generic_point_delta': baseP['LLL_generic_point_delta'],
                'LLL_chiral_operator_residual': baseP['LLL_chiral_operator_residual'],
                'zero_order_energy_residual': float(np.max(abs(frame['band_energy'][0]-expected)))}
    arrays['metadata_json'] = np.frombuffer(json.dumps(metadata, sort_keys=True).encode(), np.uint8)
    raw = deterministic_npz_bytes(arrays); digest = hashlib.sha256(raw).hexdigest()
    path = ROOT/f'results/development/artifacts/nsc-reference-band-action.{digest}.npz'
    path.write_bytes(raw); return path


def make_record(path):
    a = load(path); meta = json.loads(a['metadata_json'].tobytes())
    for p, expected in meta['signature'].items():
        if sha(p) != expected: raise ValueError('band-action dependency changed: '+p)
    maxima = {k: float(np.max(v)) for k, v in a.items() if k.startswith('scaled_residuals/')}
    maxima.update({k: meta[k] for k in ('previous_projector_residual', 'leading_overlap_identity',
        'zero_order_energy_residual', 'LLL_generic_point_delta', 'LLL_chiral_operator_residual')})
    raw = {k: float(np.max(v)) for k, v in a.items() if k.startswith('raw_residuals/')}
    for name in ('identity_scaled', 'conjugation_scaled'):
        maxima[name] = max(float(np.max(a[f'variation/{B}/{name}'])) for B in range(4))
    for name in ('identity_absolute', 'conjugation_absolute'):
        raw[name] = max(float(np.max(a[f'variation/{B}/{name}'])) for B in range(4))
    maxima['real_band_action'] = float(np.max(abs(a['band_energy'].imag)))
    tolerances = {k: 3e-8 if k in ('identity_scaled', 'conjugation_scaled') else 3e-11 for k in maxima}
    failures = {k: v for k, v in maxima.items() if v > tolerances[k]}
    generic = load(ROOT/GENERIC_DIAGNOSTIC)
    initial = {name: max(float(np.max(generic[f'variation/{B}/{name}'])) for B in range(4))
               for name in ('identity_scaled', 'conjugation_scaled')}
    families = {}
    for j, label in enumerate(meta['labels']):
        families[label] = {
            'band_energy_by_order': a['band_energy'][:, j].real.tolist(),
            'frame_maximum_scaled_residual': max(float(v[:, j].max()) for k, v in a.items() if k.startswith('scaled_residuals/')),
            'variation_maximum_scaled_residual': max(float(a[f'variation/{B}/identity_scaled'][:, j].max()) for B in range(4)),
            'time_connection_maximum': float(np.max(abs(a['connection_coefficients'][:, j]))),
            'local_star_trace_exchange_maximum': max(float(np.max(abs(a[f'variation/{B}/star_trace_defect_density'][:, j]))) for B in range(4)),
        }
    return {
        'schema': 'NSC-REFERENCE-BAND-ACTION-v1', 'accountable_author': 'Douglas Ek',
        'status': 'PASS: formal reference action density and variation; finite matching OPEN' if not failures else 'OPEN',
        'source_hashes': {p: sha(p) for p in SOURCES}, 'input_hashes': {p: sha(p) for p in INPUTS},
        'payload': {'path': str(path.relative_to(ROOT)), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest(), 'bytes': path.stat().st_size},
        'locked_inputs': read(PREVIOUS)['locked_inputs'], 'maxima': maxima, 'raw_residuals': raw,
        'tolerances': tolerances, 'failures': failures, 'per_symbol': families,
        'numerical_repair': {'initial_generic_reference_diagnostic': initial,
            'artifact': {'path': GENERIC_DIAGNOSTIC, 'sha256': sha(GENERIC_DIAGNOSTIC)},
            'cause': 'generic LLL square-root/normalization cancellation amplified under metric differentiation',
            'repair': 'same exact chiral operator projector; no threshold loosening, no new state',
            'massive_channels': 'full star-normalized frame retained'},
        'domain': {'formal_order': 4, 'point': {'PG_tau': .019, 'rho': .23},
            'history': 'same supplied metric amplitudes(.002,.001,.003,.001),omega.4; no history selection',
            'phase_space': 'actual KS T,z and canonical k; local gapped reference-band charts',
            'canonical_state': 'unchanged; reference frame is not Cauchy transport or a Gaussian covariance',
            'action_density': 'tr(Pi*(U star H star Udag + i*epsilon*dT(U) star Udag))',
            'CTP_subtraction': 'same reference phase on plus branch minus minus branch; positive reference-vertex sign',
            'trace_boundary_terms': 'explicit time derivative and local star-trace exchange retained, not set to zero',
            'finite_difference_amplitude_step': STEP, 'finite_difference_stencil': 'fourth-order centered',
            'local_star_exchange_interpretation': 'density before integration; neither an evaluated cutoff residual nor physical stress',
            'regulator': 'not assigned by this formal construction; use the declared covariant heat matching',
            'LLL': 'exact constant chiral reference on each gapped momentum-sign chart; not reused for any massive/angular channel',
            'imported_method': 'PST space-adiabatic band construction, applied to owned NSC symbols'},
        'finite_matching': {
            'identity': '[Gamma_G,vac^E-Gamma_ref,matched^E+Gamma_local,locked^E-Gamma_heat,full^E]_(<=4)=0',
            'locked_complement': 'H_mu(y)=h_Lambda(y)-E1(y/mu^2); mu_B=sqrt(abs(q))/r_e',
            'required_inventory': 'same light and retained positive compact canonical channels on both sides',
            'not_allowed': 'replace heat matching by a Gaussian of shifted coordinate energy, ignore cutoff trace defects, or refit coefficients',
            'status': 'OPEN: finite general-history regulator and measure conversion not evaluated'},
        'gate': {'formal_reference_action': 'PASS' if not failures else 'OPEN',
            'finite_matched_reference': 'OPEN', 'complete_source': 'OPEN', 'physical_EndpointBranchJets': 'OPEN',
            'B2': 'OPEN', 'V_c': None, 'extended_stationarity': 'OPEN', 'absolute_stress': None,
            'nulls': None, 'updated_constraints': None, 'physical_C0_changed': False, 'metric_timestep': False,
            'history_selected': False, 'Gamma_rest_extra_independent_term': 'none introduced',
            'Z3': 'OUT OF SCOPE', 'Weyl_time_node_diagnostic': 93.54264532195464,
            'homogeneous_nonexistence': 'preserved in original scope, not rerun'},
        'comparison': {'fields': 'all', 'float_atol': 3e-11, 'float_rtol': 3e-10, 'exact': 'hashes, scope, structure'},
        'reproducer': 'python3 scripts/derive_nsc_reference_band_action.py --check',
    }


def main():
    p = argparse.ArgumentParser(); p.add_argument('--prepare', action='store_true'); p.add_argument('--check', action='store_true')
    p.add_argument('--workers', type=int, default=3); args = p.parse_args()
    if args.prepare == args.check: p.error('choose prepare or check')
    path = prepare(args.workers) if args.prepare else ROOT/read(OUTPUT)['payload']['path']
    result = make_record(path)
    if args.check:
        from derive_nsc_transmitting_boundary_binding import compare
        compare(read(OUTPUT), result)
    else: (ROOT/OUTPUT).write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    print(json.dumps({k: result[k] for k in ('status', 'maxima', 'failures')}, indent=2), flush=True)
    if result['failures']: raise SystemExit(1)


if __name__ == '__main__': main()
