#!/usr/bin/env python3
"""New finite subgap/history source panel; historical inputs are read only."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import inspect
import json
from pathlib import Path
import sys
import tempfile

import numpy as np
from numpy.polynomial.legendre import leggauss

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from recursive_horizons.nsc_paired_horizon_preparation import PairedHorizonSeedMap
from recursive_horizons.nsc_complex_horizon_modes import complex_reflection
from recursive_horizons.nsc_subgap_history_response import (
    inner_horizon_basis, analytic_columns_on_grid, paired_history_response,
    AnalyticResponsePanel, integrate_subgap_response,
)
from recursive_horizons.nsc_transmitting_history_modes import SuppliedKSHarmonicMetric
from recursive_horizons.nsc_transmitting_history_jets import FourthOrderModePropagator
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes
from derive_nsc_transmitting_history_modes import TIME, AMPLITUDES, OMEGA, payload, load
import derive_nsc_pg_threshold_projection as previous
import derive_nsc_pg_retarded_threshold as retarded

OUTPUT = 'results/development/nsc-subgap-history-response.json'
RESTART = 'results/development/nsc-compact-matched-restart.json'
MODE = 'results/development/nsc-pg-ctp-mode-jets.json'
JETS = 'results/development/nsc-transmitting-history-jets.json'
SOURCES = ('src/recursive_horizons/nsc_subgap_history_response.py',
           'scripts/derive_nsc_subgap_history_response.py',
           'tests/test_nsc_subgap_history_response.py', 'docs/nsc-subgap-history-response.md')
INPUTS = (RESTART, MODE, JETS, 'results/development/nsc-pg-retained-subgap.json',
          'src/recursive_horizons/nsc_complex_horizon_modes.py',
          'src/recursive_horizons/nsc_paired_horizon_preparation.py',
          'src/recursive_horizons/nsc_pg_massive_modes.py',
          'src/recursive_horizons/nsc_transmitting_history_modes.py',
          'src/recursive_horizons/nsc_transmitting_history_jets.py',
          'src/recursive_horizons/nsc_transmitting_ctp_variation.py',
          'src/recursive_horizons/nsc_pg_threshold_projection.py',
          'src/recursive_horizons/nsc_ks_spacetime_variation.py',
          'src/recursive_horizons/nsc_pg_ks_metric_pullback.py',
          'src/recursive_horizons/nsc_transmitting_dirac_domain.py',
          'src/recursive_horizons/nsc_common_time_bulk_split.py',
          'src/recursive_horizons/nsc_lorentzian.py',
          'docs/nsc-declared-action-scope.md')
CONFIGS = ((401, 64), (801, 64), (801, 128))


def read(p): return json.loads((ROOT/p).read_text())
def sha(p): return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()


def preparation(config):
    return PairedHorizonSeedMap(*(config[k] for k in
        ('horizon_rho', 'surface_gravity', 'omega', 'horizon_offset', 'scattering_tolerance')))


def response_node(job):
    z, config, configs = job
    p = preparation(config); m = np.pi/2
    v = inner_horizon_basis(z, m, 0., p)
    vd = v if not complex(z).imag else inner_horizon_basis(complex(z).conjugate(), m, 0., p)
    arrays = {'at_zero': v, 'at_zero_dual': vd}
    grids = {}
    for points, steps in configs:
        if points not in grids:
            x = np.linspace(-2., 1.2, points)
            upper = analytic_columns_on_grid(z, m, 0., v, x)
            dual = upper if not complex(z).imag else analytic_columns_on_grid(complex(z).conjugate(), m, 0., vd, x)
            grids[points] = FourthOrderModePropagator(x, m, 0.), upper, dual
        owner, upper, dual = grids[points]
        row = paired_history_response(owner, z, upper, dual, np.linspace(*TIME, steps+1),
                                     SuppliedKSHarmonicMetric(tuple(AMPLITUDES), OMEGA, TIME))
        key = f'{points}_{steps}'
        arrays[key+'/response'] = row.pop('response')
        arrays[key+'/dual_response'] = row.pop('dual_response')
        arrays[key+'/diagnostics_json'] = np.frombuffer(json.dumps(row, sort_keys=True).encode(), np.uint8)
    return z, arrays


def reflection_nodes(config, channel):
    """Reuse only caches keyed by the exact old node owner and configuration."""
    old_signature = {'source_hashes': {q: previous.sha(q) for q in previous.DEPENDENCIES},
                     'node_function_sha256': hashlib.sha256(inspect.getsource(previous.evaluate_node).encode()).hexdigest(),
                     'config': config, 'channel': channel, 'angular_sign': 1}
    old_tag = hashlib.sha256(json.dumps(old_signature, sort_keys=True).encode()).hexdigest()[:16]
    signature = {**old_signature, 'retarded_sources': {q: previous.sha(q) for q in
        ('src/recursive_horizons/nsc_pg_retarded_packets.py', 'src/recursive_horizons/nsc_pg_retarded_threshold.py')},
        'retarded_node_function_sha256': hashlib.sha256(inspect.getsource(retarded.evaluate).encode()).hexdigest()}
    tag = hashlib.sha256(json.dumps(signature, sort_keys=True).encode()).hexdigest()[:16]
    caches = [Path(tempfile.gettempdir())/('nsc-threshold-contour-'+old_tag),
              Path(tempfile.gettempdir())/('nsc-retarded-threshold-'+tag)]
    p = preparation(config); m = np.pi/2
    nodes = []
    for points, fraction in ((24, 3), (48, 3), (48, 4)):
        x, _ = leggauss(points); h = p.surface_gravity/fraction
        for a, b in ((1.+0j, 1.+1j*h), (1.+1j*h, m+1j*h), (m+1j*h, m+0j)):
            nodes.extend(a+(xi+1)*(b-a)/2 for xi in x)
    values, reused = {}, 0
    for z in dict.fromkeys(nodes):
        hit = None
        for cache in caches:
            file = cache/(previous.node_key(z)+'.npz')
            if file.exists():
                with np.load(file, allow_pickle=False) as data:
                    if 'reflection' in data and data['frequency'].item() == z:
                        hit = data['reflection'].item(); break
        if hit is None:
            hit = complex_reflection(z, m, 0., p.background, radial_collar=p.collar_delta_q).reflection
        else:
            reused += 1
        values[z] = hit
    return np.array(list(values)), np.array(list(values.values())), reused


def prepare(workers):
    config = read(RESTART)['scattering_provenance']['config']
    m = np.pi/2; energies = (1+m)/2+(m-1)/2*np.cos(np.pi*np.arange(17)/16)
    zcheck = (1+m)/2+1j*config['surface_gravity']/3
    mode = payload(read(MODE)); oldjet = payload(read(JETS))
    Econtrol = float(mode['13_1/energies'][2])
    jobs = [(complex(E), config, CONFIGS) for E in energies]
    jobs += [(zcheck, config, (CONFIGS[-1],)), (complex(Econtrol), config, (CONFIGS[1],))]
    arrays = {'energies': energies, 'complex_control_frequency': np.array(zcheck),
              'normalization_control_frequency': np.array(Econtrol)}
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for j, (z, rows) in enumerate(pool.map(response_node, jobs), 1):
            prefix = f'node{j-1}' if j <= len(energies) else ('complex_control' if j == len(energies)+1 else 'normalization_control')
            arrays.update({prefix+'/'+k: v for k, v in rows.items()})
            print(f'{j}/{len(jobs)} new subgap response nodes', flush=True)
    channel = next(c for c in read('results/development/nsc-mode-resolved-cauchy-state.json')['channels'] if c['index'] == 13)
    frequency, reflection, reused = reflection_nodes(config, channel)
    arrays['reflection_frequency'], arrays['reflection'] = frequency, reflection
    print(f'{len(frequency)} contour reflection samples: {reused} reused from exact prior cache', flush=True)
    arrays['control/physical_at_zero'] = mode['13_1/at_zero'][2]
    arrays['control/old_overlap'] = oldjet['13_1/801_64/overlap'][:, 6:9, 6:9]
    metadata = {'signature': {p: sha(p) for p in (*SOURCES, *INPUTS)}, 'config': config,
                'mode_payload': read(MODE)['payload'], 'oldjet_payload': read(JETS)['payload'],
                'reused_reflections': reused, 'all_reflections': len(frequency)}
    arrays['metadata_json'] = np.frombuffer(json.dumps(metadata, sort_keys=True).encode(), np.uint8)
    raw = deterministic_npz_bytes(arrays); digest = hashlib.sha256(raw).hexdigest()
    path = ROOT/f'results/development/artifacts/nsc-subgap-history-response.{digest}.npz'
    path.write_bytes(raw)
    return path


def make_record(path):
    a = load(path); meta = json.loads(a['metadata_json'].tobytes())
    for p, expected in meta['signature'].items():
        if sha(p) != expected: raise ValueError('subgap history dependency changed: '+p)
    E = a['energies']; m = np.pi/2; kappa = meta['config']['surface_gravity']
    reflection = dict(zip(a['reflection_frequency'], a['reflection']))
    tables, maxima = {}, {}
    for points, steps in CONFIGS:
        key = f'{points}_{steps}'
        B = np.array([a[f'node{i}/{key}/response'] for i in range(len(E))])
        panel = AnalyticResponsePanel(E, B, 1., m)
        for quad, fraction in ((24, 3), (48, 3), (48, 4)):
            value = integrate_subgap_response(panel, kappa, reflection.__getitem__, points=quad, height=kappa/fraction)
            tables[f'{key}/q{quad}_h{fraction}'] = {k: (v.tolist() if not np.iscomplexobj(v) else
                {'real': v.real.tolist(), 'imag': v.imag.tolist()}) if isinstance(v, np.ndarray) else v for k, v in value.items()}
        for i in range(len(E)):
            diag = json.loads(a[f'node{i}/{key}/diagnostics_json'].tobytes())
            for name in ('analytic_adjoint_residual', 'metric_flux_tangent', 'norm_flux', 'tangent_boundary_trace'):
                maxima[name] = max(maxima.get(name, 0.), diag[name])
        if key == '801_128':
            low = AnalyticResponsePanel(E[::2], B[::2], 1., m)
            lowvalue = integrate_subgap_response(low, kappa, reflection.__getitem__)
            value = integrate_subgap_response(panel, kappa, reflection.__getitem__)
            maxima['spectral_interpolation'] = float(np.max(abs(value['signed_action_derivative']-lowvalue['signed_action_derivative'])))
            maxima['independent_complex_pair'] = float(np.max(abs(panel(a['complex_control_frequency'].item())-a['complex_control/801_128/response'])))
            diag = json.loads(a['complex_control/801_128/diagnostics_json'].tobytes())
            maxima['complex_pair_adjoint'] = diag['analytic_adjoint_residual']
    get = lambda key: np.array(tables[key]['signed_action_derivative'])
    maxima['spatial_source_refinement'] = float(np.max(abs(get('801_64/q48_h3')-get('401_64/q48_h3'))))
    maxima['time_source_refinement'] = float(np.max(abs(get('801_128/q48_h3')-get('801_64/q48_h3'))))
    maxima['contour_quadrature'] = float(np.max(abs(get('801_128/q48_h3')-get('801_128/q24_h3'))))
    maxima['contour_deformation'] = float(np.max(abs(get('801_128/q48_h3')-get('801_128/q48_h4'))))
    v = a['normalization_control/at_zero']; physical = a['control/physical_at_zero']
    R = np.vdot(v[:, 0], physical[:, 0])/np.vdot(v[:, 0], v[:, 0])
    D = np.diag([R, 1.]); reconstructed = np.column_stack((R*v[:, 0], v[:, 1], np.zeros(2)))
    maxima['owned_column_normalization'] = float(np.max(abs(reconstructed-physical)))
    maxima['owned_control_reflection_modulus'] = float(abs(abs(R)-1))
    B = a['normalization_control/801_64/response']
    mapped = np.array([D.conj().T@(-1j*b)@D for b in B])
    maxima['owned_Gaussian_overlap'] = float(np.max(abs(mapped-a['control/old_overlap'][:, :2, :2])))
    tolerances = {k: 3e-8 if k in ('spatial_source_refinement', 'time_source_refinement', 'tangent_boundary_trace') else 3e-11 for k in maxima}
    failures = {k: v for k, v in maxima.items() if v > tolerances[k]}
    return {
        'schema': 'NSC-SUBGAP-HISTORY-RESPONSE-v1', 'accountable_author': 'Douglas Ek',
        'status': 'PASS: group13 finite subgap Gaussian history variation' if not failures else 'OPEN',
        'source_hashes': {p: sha(p) for p in SOURCES}, 'input_hashes': {p: sha(p) for p in INPUTS},
        'payload': {'path': str(path.relative_to(ROOT)), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest(), 'bytes': path.stat().st_size},
        'locked_inputs': read(MODE)['locked_inputs'], 'maxima': maxima, 'tolerances': tolerances, 'failures': failures,
        'domain': {'group': 13, 'mass': m, 'angular': 0., 'signed_panels': [[-m, -1.], [1., m]],
            'normalization': 'bare Gaussian derivative per reduced radial family; no angular/copy/compact multiplier',
            'state': 'owned horizon pair and whole-line subgap reflection; seed covariance not input',
            'history': {'coordinates': list(TIME), 'amplitudes': AMPLITUDES.tolist(), 'omega': OMEGA, 'selected': False},
            'source_formula': '-(smooth+2*Re(reflected))/pi; exact signed pair, not a new subtraction',
            'complex_frequency': 'analytic bilinear using independently propagated conjugate-energy dual, never a Gaussian state',
            'causal_scope': 'compact perturbation in trapped patch; exterior u=0 in its initial domain of dependence',
            'interpolation': '17 unsewn real Chebyshev nodes; nested9 control; reflection never interpolated'},
        'source_integrals': tables,
        'numerics': {'configurations': [list(c) for c in CONFIGS], 'reflection_samples': meta['all_reflections'],
                     'reused_reflections': meta['reused_reflections'], 'old_generators_rerun': False},
        'gate': {'finite_group13_subgap_source': 'PASS' if not failures else 'OPEN',
                 'full_spectral_reference_source': 'OPEN: remaining panels/families and inhomogeneous reference/local variation',
                 'physical_EndpointBranchJets': 'OPEN', 'B2': 'OPEN', 'V_c': None, 'extended_stationarity': 'OPEN',
                 'absolute_stress': None, 'nulls': None, 'updated_constraints': None, 'metric_timestep': False,
                 'Z3': 'OUT OF SCOPE', 'Weyl_time_node_diagnostic': 93.54264532195464,
                 'Gamma_rest_extra_term': 'none introduced', 'homogeneous_nonexistence': 'preserved, not rerun'},
        'comparison': {'fields': 'all', 'float_atol': 3e-12, 'float_rtol': 3e-10, 'exact': 'hashes, scope, structure'},
        'reproducer': 'python3 scripts/derive_nsc_subgap_history_response.py --check',
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
    else:
        (ROOT/OUTPUT).write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    print(json.dumps({k: result[k] for k in ('status', 'maxima', 'failures')}, indent=2), flush=True)
    if result['failures']: raise SystemExit(1)


if __name__ == '__main__': main()
