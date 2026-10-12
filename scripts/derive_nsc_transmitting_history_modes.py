#!/usr/bin/env python3
"""New nonlinear spatial-mode history calculation, with scoped refinements."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from numpy.polynomial.legendre import leggauss

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from recursive_horizons.nsc_transmitting_history_modes import (
    RelativeModePropagator, SuppliedKSHarmonicMetric, reference_fields_on_grid,
)
from recursive_horizons.nsc_transmitting_resolvent import profile
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes

OUTPUT = 'results/development/nsc-transmitting-history-modes.json'
MODE = 'results/development/nsc-pg-ctp-mode-jets.json'
LINEAR = 'results/development/nsc-ks-spacetime-variation.json'
SOURCES = ('src/recursive_horizons/nsc_transmitting_history_modes.py',
           'scripts/derive_nsc_transmitting_history_modes.py',
           'tests/test_nsc_transmitting_history_modes.py',
           'docs/nsc-transmitting-history-modes.md')
INPUTS = (MODE, LINEAR,
          'src/recursive_horizons/nsc_ks_spacetime_variation.py',
          'src/recursive_horizons/nsc_pg_ks_metric_pullback.py',
          'src/recursive_horizons/nsc_common_time_bulk_split.py',
          'src/recursive_horizons/nsc_transmitting_dirac_domain.py',
          'src/recursive_horizons/nsc_transmitting_resolvent.py',
          'src/recursive_horizons/nsc_lorentzian.py',
          'docs/nsc-declared-action-scope.md')
AMPLITUDES = np.array([.002, .001, .003, .001])
TIME = (0., .05)
OMEGA = .4
CONFIGS = ((401, 64), (801, 64), (801, 128))


def sha(path):
    return hashlib.sha256((ROOT/path).read_bytes()).hexdigest()


def read(path):
    return json.loads((ROOT/path).read_text())


def load(path):
    with np.load(path, allow_pickle=False) as a:
        return {k: a[k].copy() for k in a.files}


def payload(record):
    p = record['payload']
    if sha(p['path']) != p['sha256']:
        raise ValueError('input payload hash changed: '+p['path'])
    return load(ROOT/p['path'])


def linear_prediction(data, vertex):
    """Independent time quadrature of the locked reference Duhamel kernel."""
    E = np.repeat(data['energies'], 3)
    M = np.einsum('a,aij->ij', AMPLITUDES, vertex)
    result = np.zeros_like(M)
    t, w = leggauss(64)
    ti, tf = TIME
    for tau, weight in zip((ti+tf)/2+(tf-ti)*t/2, (tf-ti)*w/2):
        bump = profile(2*(tau-ti)/(tf-ti)-1)[0]
        H1 = .5*(np.exp(-1j*OMEGA*tau)*M+np.exp(1j*OMEGA*tau)*M.conj().T)
        result += -1j*weight*bump*np.exp(-1j*E*(tf-tau))[:, None]*H1*np.exp(-1j*E*(tau-ti))[None, :]
    return result


def one_family(job):
    name, data, vertex = job
    E = data['energies']
    mass, angular = float(data['mass'].item()), float(data['angular'].item())
    arrays, diagnostics, grids = {}, {}, {}
    for points, steps in CONFIGS:
        x = np.linspace(-2., 1.2, points)
        if points not in grids:
            fields = reference_fields_on_grid(E, mass, angular, data['at_zero'], x)
            grids[points] = (RelativeModePropagator(x, mass, angular), fields)
        owner, fields = grids[points]
        result = owner.evolve(E, fields, np.linspace(*TIME, steps+1),
                             SuppliedKSHarmonicMetric(tuple(AMPLITUDES), OMEGA, TIME))
        key = f'{points}_{steps}'
        arrays[key+'/projection'] = result['mode_projection']
        diagnostics[key] = {k: v for k, v in result.items() if not isinstance(v, np.ndarray)}
        if (points, steps) == CONFIGS[-1]:
            arrays['x'] = x
            arrays['reference_field'] = result['reference_characteristic_field']
            arrays['difference_field'] = result['difference_characteristic_field']
            arrays['difference_L2_norms'] = result['difference_L2_norms']
    arrays['linear_prediction'] = linear_prediction(data, vertex)
    # One mixed massive family compares the new PDE derivative independently
    # against the already verified analytic source-mode Duhamel insertion.
    if name == '14_1':
        owner, fields = grids[CONFIGS[-1][0]]
        for sign in (-1, 1):
            result = owner.evolve(E, fields, np.linspace(*TIME, CONFIGS[-1][1]+1),
                SuppliedKSHarmonicMetric(tuple(sign*.05*AMPLITUDES), OMEGA, TIME))
            arrays[f'linear_limit_{sign}'] = result['mode_projection']
    arrays['diagnostics_json'] = np.frombuffer(json.dumps(diagnostics, sort_keys=True).encode(), np.uint8)
    return name, arrays


def prepare(workers):
    modes, linear = payload(read(MODE)), payload(read(LINEAR))
    names = sorted({k.split('/')[0] for k in modes if '/' in k})
    jobs = [(name, {k[len(name)+1:]: v for k, v in modes.items() if k.startswith(name+'/')},
             linear[name+'/vertex']) for name in names]
    arrays = {}
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for n, (name, result) in enumerate(pool.map(one_family, jobs), 1):
            arrays.update({name+'/'+k: v for k, v in result.items()})
            if n % 8 == 0 or n == len(names):
                print(f'{n}/{len(names)} new spatial history controls complete', flush=True)
    arrays['metadata_json'] = np.frombuffer(json.dumps({
        'signature': {p: sha(p) for p in (*SOURCES, *INPUTS)},
        'mode_payload': read(MODE)['payload'], 'linear_payload': read(LINEAR)['payload'],
    }, sort_keys=True).encode(), np.uint8)
    raw = deterministic_npz_bytes(arrays)
    digest = hashlib.sha256(raw).hexdigest()
    out = ROOT/f'results/development/artifacts/nsc-transmitting-history-modes.{digest}.npz'
    out.write_bytes(raw)
    return out


def make_record(path):
    a = load(path)
    metadata = json.loads(a['metadata_json'].tobytes())
    for p, digest in metadata['signature'].items():
        if sha(p) != digest:
            raise ValueError('history source dependency changed: '+p)
    rows, maxima = {}, {}
    for name in sorted({k.split('/')[0] for k in a if '/' in k}):
        get = lambda key: a[name+'/'+key]
        diagnostics = json.loads(get('diagnostics_json').tobytes())
        fine = diagnostics['801_128']
        n = len(get('x'))
        weights = np.full(n, (get('x')[-1]-get('x')[0])/(n-1))
        weights[[0, -1]] *= .5
        weights = np.tile(weights, 2)
        rebuilt = get('reference_field').conj().T@(weights[:, None]*get('difference_field'))
        group, sign = map(int, name.split('_'))
        row = {
            'group': group, 'angular_sign': sign,
            'spatial_response_refinement': float(np.max(abs(get('801_64/projection')-get('401_64/projection')))),
            'temporal_response_refinement': float(np.max(abs(get('801_128/projection')-get('801_64/projection')))),
            'projection_reconstruction': float(np.max(abs(rebuilt-get('801_128/projection')))),
            'norm_flux_algebra': max(d['norm_flux_algebra'] for d in diagnostics.values()),
            'forcing_phase_residual': max(d['forcing_phase_residual'] for d in diagnostics.values()),
            'initial_metric_match': fine['initial_metric_match'],
            'initial_time_jet_match': fine['initial_time_jet_match'],
            'fine_left_trace': fine['left_difference_trace'],
            'fine_right_trace': fine['right_difference_trace'],
            'sampled_causal_buffer': fine['sampled_speed_buffer'],
            'maximum_right_characteristic_speed': fine['maximum_right_characteristic_speed'],
            'maximum_difference_L2_norm': float(get('difference_L2_norms').max()),
            'finite_history_minus_linear_prediction': float(np.max(abs(get('801_128/projection')-get('linear_prediction')))),
        }
        if name == '14_1':
            row['linear_limit_vs_locked_Duhamel'] = float(np.max(abs(
                (get('linear_limit_1')-get('linear_limit_-1'))/.1-get('linear_prediction'))))
        rows[name] = row
        for k in ('spatial_response_refinement', 'temporal_response_refinement', 'projection_reconstruction',
                  'norm_flux_algebra', 'forcing_phase_residual', 'initial_metric_match',
                  'fine_left_trace', 'fine_right_trace', 'linear_limit_vs_locked_Duhamel'):
            if k in row:
                maxima[k] = max(maxima.get(k, 0.), row[k])
    if len(rows) != 63 or {r['group'] for r in rows.values()} != set(range(33)):
        raise ValueError('retained group coverage changed')
    tols = {k: (3e-8 if k in ('spatial_response_refinement', 'temporal_response_refinement',
                              'linear_limit_vs_locked_Duhamel') else 3e-11) for k in maxima}
    failures = {k: v for k, v in maxima.items() if v > tols[k]}
    if any(r['sampled_causal_buffer'] <= 0 or r['maximum_right_characteristic_speed'] >= 0 for r in rows.values()):
        failures['control_propagation_domain'] = 'buffer or trapped characteristic direction violated'
    return {
        'schema': 'NSC-TRANSMITTING-HISTORY-MODES-v1', 'accountable_author': 'Douglas Ek',
        'status': 'PASS: conditional nonlinear spatial mode propagation controls' if not failures else 'OPEN: numerical history-control gate failed',
        'source_hashes': {p: sha(p) for p in SOURCES}, 'input_hashes': {p: sha(p) for p in INPUTS},
        'payload': {'path': str(path.relative_to(ROOT)), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
                    'bytes': path.stat().st_size},
        'locked_inputs': read(MODE)['locked_inputs'],
        'domain': {
            'operator': 'same full transmitting canonical Dirac field on PG slices',
            'new_unknown': 'w_E=(U[g]-U_ref)Phi_E as a spatial field',
            'equation': 'i d_tau w_E=H[g]w_E+(H[g]-H_ref)exp(-i E (tau-ti))Phi_E',
            'initial_data': 'zero difference; authenticated reference mode columns, not seed covariance',
            'computational_edges': 'zero right incoming DIFFERENCE, left outflow; no reflecting carrier',
            'auxiliary_exponential_amplitudes': 'deterministic forcing integration only; not new physical degrees of freedom',
            'energy_samples': 'inherited derivative-control real frequencies; not a complete covariance integral',
            'state_preparation': 'same C_H plus inherited occupation; preparation unchanged at pulse start',
            'metric_selection': 'none; a supplied response history is an argument of the propagator',
        },
        'numerical_control': {
            'time_coordinates': list(TIME), 'raw_KS_amplitudes': AMPLITUDES.tolist(), 'omega': OMEGA,
            'time_envelope': 'inherited C-infinity compact profile; all endpoint derivatives vanish',
            'spatial_interval': [-2., 1.2], 'grid_and_steps': [list(c) for c in CONFIGS],
            'spatial_method': 'owned second-order diagonal-norm split SBP/SAT convention',
            'time_method': 'midpoint coefficients; sparse exponential with exact exponential forcing within each step',
            'refinement_scope': 'measured differences of sampled response kernels, not rigorous whole-spectrum error bounds',
        },
        'per_family': rows, 'maxima': maxima, 'tolerances': tols, 'failures': failures,
        'gate': {
            'conditional_mode_propagation': 'PASS' if not failures else 'OPEN',
            'full_EndpointBranchJets': 'OPEN', 'complete_evolved_C_PG': 'not integrated',
            'nonlinear_history_metric_derivatives': 'not computed',
            'physical_endpoint_family': 'not selected by these response histories',
            'remaining_state_integral': 'retain all real panels, analytic subgap representation and tail; no complex contour node promoted to canonical frequency',
            'Gamma_rest_extra_term': 'none introduced; full Gaussian action counted once',
            'B2': 'OPEN', 'V_c': None, 'extended_stationarity': 'OPEN',
            'Weyl_time_node_diagnostic': 93.54264532195464,
            'absolute_stress': None, 'new_nulls': None, 'updated_constraints': None,
            'metric_timestep': False, 'history_selected': False, 'Z3': 'OUT OF SCOPE',
            'homogeneous_nonexistence': 'imported in its original scope, not rerun',
        },
        'comparison': {'fields': 'all', 'float_atol': 3e-12, 'float_rtol': 3e-10,
                       'exact': 'source/input hashes, structure and scope labels'},
        'reproducer': 'python3 scripts/derive_nsc_transmitting_history_modes.py --check',
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--prepare', action='store_true')
    p.add_argument('--check', action='store_true')
    p.add_argument('--workers', type=int, default=2)
    args = p.parse_args()
    if args.prepare == args.check:
        p.error('select --prepare or --check')
    path = prepare(args.workers) if args.prepare else ROOT/read(OUTPUT)['payload']['path']
    record = make_record(path)
    if args.check:
        from derive_nsc_transmitting_boundary_binding import compare
        compare(read(OUTPUT), record)
    else:
        (ROOT/OUTPUT).write_text(json.dumps(record, indent=2, sort_keys=True)+'\n')
    print(json.dumps({'status': record['status'], 'maxima': record['maxima'], 'failures': record['failures']}, indent=2), flush=True)
    if record['failures']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
