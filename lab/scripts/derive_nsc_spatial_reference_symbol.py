#!/usr/bin/env python3
"""NSC-specific application of the formal reference-symbol recursion."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT/'src'))
from recursive_horizons.nsc_spatial_reference_symbol import (
    reference_projector, homogeneous_seed_metric_jets, supplied_KS_metric_jets, SIGMA,
)
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes
from recursive_horizons.nsc_transmitting_history_modes import SuppliedKSHarmonicMetric
from recursive_horizons.nsc_pg_ks_metric_pullback import pg_to_ks, reference_chart

OUTPUT = 'results/development/nsc-spatial-reference-symbol.json'
SEED = 'results/development/nsc-mode-resolved-cauchy-state.json'
MODE = 'results/development/nsc-pg-ctp-mode-jets.json'
SOURCES = ('src/recursive_horizons/nsc_spatial_reference_symbol.py',
           'scripts/derive_nsc_spatial_reference_symbol.py', 'tests/test_nsc_spatial_reference_symbol.py',
           'docs/nsc-spatial-reference-symbol.md')
INPUTS = (SEED, MODE, 'results/development/nsc-spherical-local-history.json',
          'src/recursive_horizons/nsc_general_ks_reference.py',
          'src/recursive_horizons/nsc_compact_ctp_neck.py',
          'src/recursive_horizons/nsc_ks_spacetime_variation.py',
          'src/recursive_horizons/nsc_pg_ks_metric_pullback.py',
          'src/recursive_horizons/nsc_transmitting_history_modes.py',
          'docs/nsc-adm-source-constraints.md', 'docs/nsc-declared-action-scope.md')


def read(p): return json.loads((ROOT/p).read_text())
def sha(p): return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def load(p):
    with np.load(p, allow_pickle=False) as f: return {k: f[k].copy() for k in f.files}
def payload(record):
    p = record['payload']
    if sha(p['path']) != p['sha256']: raise ValueError('authenticated payload changed')
    return load(ROOT/p['path'])


def metric_check(fields):
    tau, rho = .019, .23; amplitudes = (.002, .001, .003, .001)
    owner = SuppliedKSHarmonicMetric(amplitudes, .4, (0., .05))
    def values(t, x): return pg_to_ks(x, *owner.values(t, np.array([x]))[:, 0])
    direct = values(tau, rho)
    value = np.array([f.value[0, 0, 0].real for f in fields])
    dt, dx = 8e-5, 1e-3
    Dtau = (values(tau-2*dt, rho)-8*values(tau-dt, rho)+8*values(tau+dt, rho)-values(tau+2*dt, rho))/(12*dt)
    Drho = (values(tau, rho-2*dx)-8*values(tau, rho-dx)+8*values(tau, rho+dx)-values(tau, rho+2*dx))/(12*dx)
    b, _, a, _, _ = reference_chart(rho)
    DT = b/a*Dtau-a*Drho
    gotT = np.array([f.derivative(t=1).value[0, 0, 0].real for f in fields])
    gotz = np.array([f.derivative(z=1).value[0, 0, 0].real for f in fields])
    return {'value': float(np.max(abs(value-direct))),
            'T_derivative': float(np.max(abs(gotT-DT))), 'z_derivative': float(np.max(abs(gotz-Dtau)))}


def prepare():
    seed = read(SEED); data = payload(seed); modes = payload(read(MODE))
    count = len(data['frequency']); mass = np.zeros(count); angular = mass.copy()
    for channel in seed['channels']:
        sl = slice(channel['sample_offset'], channel['sample_offset']+channel['sample_count'])
        mass[sl], angular[sl] = channel['compact_mass'], channel['angular_eigenvalue']
    homogeneous = []; fields = homogeneous_seed_metric_jets()
    for lo in range(0, count, 128):
        hi = min(lo+128, count)
        row = reference_projector(*fields, -data['frequency'][lo:hi], mass[lo:hi], angular[lo:hi])
        homogeneous.append(row['orders'])
    arrays = {'homogeneous_orders': np.concatenate(homogeneous, axis=1)}
    names = sorted({k.split('/')[0] for k in modes if '/' in k})
    labels, mm, ll, kk = [], [], [], []
    for name in names:
        for sign in (1, -1):
            labels.append(f'{name}/k{sign:+d}')
            mm.append(modes[name+'/mass'].item()); ll.append(modes[name+'/angular'].item())
            kk.append(sign*float(modes[name+'/energies'][2]))
    fields = supplied_KS_metric_jets(.019, .23)
    result = reference_projector(*fields, kk, mm, ll)
    arrays['spatial_orders'] = result['orders']; arrays['gap'] = result['gap']
    arrays['mass'], arrays['angular'], arrays['momentum'] = np.array(mm), np.array(ll), np.array(kk)
    for n, row in enumerate(result['residuals'], 1):
        for name, values in row.items():
            for kind, values in values.items(): arrays[f'residual/{n}/{name}/{kind}'] = values
    N, beta, axial, radius = fields
    av, rv, rz = axial.value[0, 0, 0], radius.value[0, 0, 0], radius.derivative(z=1).value[0, 0, 0]
    normal_gap = np.sqrt(np.array(mm)**2+(np.array(ll)/rv)**2+(np.array(kk)/av)**2)
    arrays['first_scalar_trace_expected'] = -.5*np.array(mm)*np.array(ll)*rz/(av*rv**2*normal_gap**3)
    metadata = {'signature': {p: sha(p) for p in (*SOURCES, *INPUTS)},
                'seed_payload': seed['payload'], 'mode_payload': read(MODE)['payload'],
                'labels': labels, 'metric_check': metric_check(fields), 'physical_point': {'PG_tau': .019, 'rho': .23}}
    arrays['metadata_json'] = np.frombuffer(json.dumps(metadata, sort_keys=True).encode(), np.uint8)
    raw = deterministic_npz_bytes(arrays); digest = hashlib.sha256(raw).hexdigest()
    path = ROOT/f'results/development/artifacts/nsc-spatial-reference-symbol.{digest}.npz'
    path.write_bytes(raw); return path


def make_record(path):
    a = load(path); meta = json.loads(a['metadata_json'].tobytes())
    for p, expected in meta['signature'].items():
        if sha(p) != expected: raise ValueError('reference-symbol dependency changed: '+p)
    data = payload({'payload': meta['seed_payload']}); seed = read(SEED)
    P = a['homogeneous_orders'].sum(axis=0); extra = a['homogeneous_orders'][1:].sum(axis=0)
    controls = {}; maximum = {'compact': 0., 'angular_energy': 0., 'angular_parallel': 0., 'LLL': 0.}
    for c in seed['channels']:
        lo, hi = c['sample_offset'], c['sample_offset']+c['sample_count']; sl = slice(lo, hi)
        if c['compact_mass'] > 0 or c['angular_eigenvalue'] == 0:
            key = 'compact' if c['compact_mass'] > 0 else 'LLL'
            row = {key: float(np.max(abs(P[sl]-data['covariance_reference_seed'][sl])))}
        else:
            h = data['hy_seed'][sl, None, None]*SIGMA[1]+data['hz_seed'][sl, None, None]*SIGMA[2]
            energy = np.einsum('fij,fji->f', extra[sl], h).real
            pressure = data['hz_seed'][sl]*np.einsum('fij,ji->f', extra[sl], SIGMA[2]).real
            row = {'angular_energy': float(np.max(abs(energy-data['reference_energy'][sl]))),
                   'angular_parallel': float(np.max(abs(pressure-data['reference_parallel'][sl])))}
        controls[str(c['index'])] = {'count': c['sample_count'], 'residuals': row}
        for key, value in row.items(): maximum[key] = max(maximum[key], value)
    per_family = {}; absolute = {}; scaled = {}
    trace = np.trace(a['spatial_orders'][1], axis1=-2, axis2=-1)
    for j, label in enumerate(meta['labels']):
        row = {'gap': float(a['gap'][j]), 'first_trace_real': float(trace[j].real),
               'first_trace_imaginary': float(trace[j].imag), 'orders': {}}
        for n in range(1, 5):
            r = {}
            for name in ('idempotence', 'transport', 'G_off_block', 'F_diagonal_block', 'Hermiticity'):
                values = {kind: float(a[f'residual/{n}/{name}/{kind}'][j]) for kind in ('absolute', 'scaled')}
                absolute[name] = max(absolute.get(name, 0.), values['absolute'])
                scaled[name] = max(scaled.get(name, 0.), values['scaled'])
                r[name] = values
            row['orders'][str(n)] = r
        per_family[label] = row
    maxima = {**{'homogeneous_'+k: v for k, v in maximum.items()},
              **{'formal_'+k: v for k, v in scaled.items()},
              **{'coordinate_'+k: v for k, v in meta['metric_check'].items()},
              'first_scalar_trace_formula': float(np.max(abs(trace-a['first_scalar_trace_expected'])))}
    tolerance = {k: 3e-8 if k in ('coordinate_T_derivative', 'coordinate_z_derivative') else 3e-11 for k in maxima}
    failures = {k: v for k, v in maxima.items() if v > tolerance[k]}
    return {
        'schema': 'NSC-SPATIAL-REFERENCE-SYMBOL-v1', 'accountable_author': 'Douglas Ek',
        'status': 'PASS: formal fourth-order spatial reference symbol; finite subtraction matching OPEN' if not failures else 'OPEN',
        'source_hashes': {p: sha(p) for p in SOURCES}, 'input_hashes': {p: sha(p) for p in INPUTS},
        'payload': {'path': str(path.relative_to(ROOT)), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest(), 'bytes': path.stat().st_size},
        'locked_inputs': read('results/development/nsc-spherical-local-history.json')['locked_inputs'],
        'homogeneous_controls': controls, 'spatial_families': per_family,
        'maxima': maxima, 'formal_absolute_residuals': absolute, 'tolerances': tolerance, 'failures': failures,
        'domain': {'order': 4, 'coordinates': ['KS T', 'KS z', 'canonical k'], 'physical_point': meta['physical_point'],
            'source_history': 'same owned amplitudes(.002,.001,.003,.001),omega.4,PG support(0,.05), fixed rho=0 cut',
            'homogeneous_reference_blocks': len(P), 'local_signed_symbols': len(meta['labels']),
            'sampling': 'one inherited control momentum per signed family plus its opposite sign; no spectral integral',
            'principal_symbol': '-beta*k*I+N*(-m*sigma1+lambda/r*sigma2+k/a*sigma3)',
            'ordering': 'canonical symmetric Dirac operator, Weyl symbol, fixed Pauli frame; no extra subprincipal imaginary term',
            'formal_parameter': 'derivative grading only, no physical epsilon chosen',
            'gap_condition': 'positive normal band gap; zero gap rejected',
            'residual_scaling': 'each coefficient residual divided by 1 plus the magnitude of its equation terms; absolute errors also retained',
            'pointwise_CAR': 'not applicable to formal Weyl symbol; physical C0 unchanged',
            'scalar_trace_maximum': float(np.max(abs(trace))),
            'imported_method': 'https://arxiv.org/abs/math-ph/0201055; projector recursion and time-dependent extension',
            'old_generators_rerun': False},
        'gate': {'formal_spatial_reference': 'PASS' if not failures else 'OPEN',
            'homogeneous_reference_recovery': 'PASS' if max(maximum.values()) <= 3e-11 else 'OPEN',
            'generating_reference_functional': 'OPEN: integrable subtraction variation not yet constructed',
            'finite_regulator_local_matching': 'OPEN: same-allocation identity required before source composition',
            'full_renormalized_source': 'OPEN', 'physical_EndpointBranchJets': 'OPEN', 'B2': 'OPEN',
            'V_c': None, 'extended_stationarity': 'OPEN', 'absolute_stress': None, 'nulls': None,
            'updated_constraints': None, 'physical_C0_changed': False, 'metric_timestep': False,
            'history_selected': False, 'Gamma_rest_extra_independent_term': 'none introduced',
            'Z3': 'OUT OF SCOPE', 'Weyl_time_node_diagnostic': 93.54264532195464,
            'homogeneous_nonexistence': 'preserved in original scope, not rerun'},
        'comparison': {'fields': 'all', 'float_atol': 3e-12, 'float_rtol': 3e-10, 'exact': 'hashes, scopes, structure'},
        'reproducer': 'python3 scripts/derive_nsc_spatial_reference_symbol.py --check',
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
