#!/usr/bin/env python3
"""Prepare/replay the paired incoming tail; never run old field generators."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys

import mpmath as mp
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from recursive_horizons.nsc_incoming_source_tail import (
    incoming_source_tail, transformed_tail_quadrature, evaluate_tail_series,
    leading_reference_matching_identity,
)
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes

OUTPUT = 'results/development/nsc-incoming-source-tail.json'
BASELINE = 'e8a77c69df1c6b272a58d1596a77b8f302cc8203'
FINITE = 'results/development/nsc-incoming-spectral-source.json'
SOURCES = ('src/recursive_horizons/nsc_incoming_source_tail.py',
    'tests/test_nsc_incoming_source_tail.py', 'scripts/derive_nsc_incoming_source_tail.py',
    'docs/nsc-incoming-source-tail.md')
INPUTS = (FINITE, 'results/development/nsc-mode-resolved-cauchy-state.json',
    'results/development/nsc-compact-matched-restart.json',
    'src/recursive_horizons/nsc_compact_ctp_neck.py',
    'src/recursive_horizons/nsc_pg_high_energy.py',
    'src/recursive_horizons/nsc_mode_resolved_cauchy_state.py')


def read(path): return json.loads((ROOT/path).read_text())
def sha(path): return hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
def signature(): return {p: sha(p) for p in (*SOURCES, *INPUTS)}
def load(path):
    with np.load(path, allow_pickle=False) as payload:
        return {k: payload[k].copy() for k in payload.files}


def inputs():
    finite = read(FINITE)
    published_bytes = subprocess.run(['git', 'cat-file', 'blob', BASELINE+':'+FINITE],
        cwd=ROOT, capture_output=True, check=True).stdout
    if hashlib.sha256(published_bytes).hexdigest() != sha(FINITE):
        raise ValueError('finite-source input no longer matches its declared local checkpoint')
    for key in ('source_hashes', 'input_hashes', 'numerical_owner_and_input_hashes'):
        for path, expected in finite[key].items():
            if sha(path) != expected: raise ValueError('finite-source dependency changed: '+path)
    for payload in [finite['payload'], *finite['input_payloads']]:
        if sha(payload['path']) != payload['sha256']:
            raise ValueError('authenticated source artifact changed: '+payload['path'])
    channels = read(INPUTS[1])['channels']
    config = read(INPUTS[2])['scattering_provenance']['config']
    return finite, channels, config, load(ROOT/finite['payload']['path'])


def prepare():
    finite, channels, _, _ = inputs()
    before = signature(); arrays = {}; controls = {}
    for group in range(1, 33):
        channel = channels[group]
        lower = 320. if group in (10, 11, 12, 31, 32) else 160.
        tail = incoming_source_tail(channel, lower, maximum_order=13, precision=50)
        if not tail['nonintegrable_coefficients_exactly_zero']:
            raise ArithmeticError('nonintegrable incoming-source coefficient was not structurally zero')
        prefix = f'group{group}/'
        arrays[prefix+'coefficients'] = tail['coefficients_rho_parallel_sphere']
        arrays[prefix+'quadrature16_precision50'] = transformed_tail_quadrature(channel, lower, points=16, precision=50)
        arrays[prefix+'quadrature24_precision60'] = transformed_tail_quadrature(channel, lower, points=24, precision=60)
        for order, value in tail['corrections_by_order'].items():
            arrays[prefix+f'correction_order{order}'] = value
        controls[str(group)] = {k: v for k, v in tail.items()
            if k not in ('coefficients_rho_parallel_sphere', 'corrections_by_order', 'series_order_difference')}
        if group % 8 == 0:
            print(f'paired incoming tail {group}/32; local algebra only', flush=True)
    if signature() != before: raise ValueError('tail dependency changed during preparation')
    meta = {'signature': before, 'controls': controls, 'baseline_commit': BASELINE,
            'finite_source_payload': finite['payload'],
            'mpmath_version': mp.__version__, 'numpy_version': np.__version__}
    arrays['metadata_json'] = np.frombuffer(json.dumps(meta, sort_keys=True).encode(), np.uint8)
    raw = deterministic_npz_bytes(arrays); digest = hashlib.sha256(raw).hexdigest()
    path = ROOT/f'results/development/artifacts/nsc-incoming-source-tail.{digest}.npz'
    if path.exists() and path.read_bytes() != raw: raise ValueError('content-addressed collision')
    if not path.exists(): path.write_bytes(raw)
    return path


def thermal_bound(channel, lower, config):
    """Same source-law trace-norm inequality plus elementary exponential integrals."""
    with mp.workdps(60):
        r, a = mp.sqrt(2), mp.sqrt(3*mp.pi/2-4)
        m, ell = mp.mpf(float(channel['compact_mass'])), mp.mpf(float(channel['angular_eigenvalue']))
        massnorm = mp.sqrt(m*m+(ell/r)**2)
        kappa, omega = mp.mpf(float(config['surface_gravity'])), mp.mpf(float(config['omega']))
        moments0 = moments1 = mp.mpf(0)
        for coefficient, alpha in ((2, mp.pi/kappa), (1, 2*mp.pi/(omega*kappa))):
            exponential = coefficient*mp.exp(-alpha*lower)
            moments0 += exponential/alpha
            moments1 += exponential*(lower/alpha+1/alpha**2)
        # Actual angular-sign sum cancels n_signs in the group factor.
        factor = mp.mpf(int(channel['copy_count']))*int(channel['degeneracy'])/(4*mp.pi**2*r*r*a)
        bounds = factor*mp.matrix([moments1/a+massnorm*moments0, moments1/a,
                                  moments1/a, abs(ell)/(2*r)*moments0])
        return bounds


def bound_record(values):
    with mp.workdps(60):
        return {'decimal_upper_bounds': [mp.nstr(v, 24) for v in values],
                'natural_log_upper_bounds': [float(mp.log(v)) if v else None for v in values],
                'exact_zero_components': [i for i, v in enumerate(values) if not v]}


def make_record(path):
    finite, channels, config, panels = inputs()
    arrays = load(path); meta = json.loads(arrays['metadata_json'].tobytes())
    if meta['signature'] != signature(): raise ValueError('tail producer or scientific inputs changed')
    if meta['baseline_commit'] != BASELINE or meta['finite_source_payload'] != finite['payload']:
        raise ValueError('tail has a different finite-source parent')
    if set(meta['controls']) != {str(g) for g in range(1, 33)}:
        raise ValueError('all retained non-LLL groups are required')
    if leading_reference_matching_identity()['vacuum_minus_ad1'] != ['0', '0', '0']:
        raise ArithmeticError('regular endpoint identity failed')
    groups = {}; sums = {n: np.zeros(4) for n in (7, 9, 11, 13)}
    absolute_order = np.zeros(4); thermal_total = mp.matrix([0, 0, 0, 0]); failures = {}
    maxima = {'coefficient_integral_replay': 0., 'order11_13': 0.,
              'quadrature16_24': 0., 'series13_quadrature24': 0., 'stored_band_absolute': 0.,
              'stored_band_excess_over_declared_indicator': 0.}
    tolerances = {'coefficient_integral_replay': 3e-24, 'order11_13': 3e-15,
                  'quadrature16_24': 3e-20, 'series13_quadrature24': 3e-18,
                  'stored_band_excess_over_declared_indicator': 3e-20}
    for group in range(1, 33):
        control = meta['controls'][str(group)]; prefix = f'group{group}/'
        coefficients = arrays[prefix+'coefficients']; lower = control['lower']
        if np.count_nonzero(coefficients[:2]) or not control['nonintegrable_coefficients_exactly_zero'] or control['numerical_differentiation_chop']:
            raise ArithmeticError('nonintegrable coefficients may not be clipped or ignored')
        if control['regular_endpoint_identity'] != leading_reference_matching_identity():
            raise ValueError('different regularized endpoint definition')
        stored = {order: arrays[prefix+f'correction_order{order}'] for order in sums}
        replay_error = 0.
        for order in sums:
            powers = np.arange(2, order+1)
            values = control['multiplicity_factor']*np.sum(
                coefficients[2:order+1]/((powers-1)*lower**(powers-1))[:, None], axis=0)
            replay = np.array([values[0], values[1], 0., values[2]])
            replay_error = max(replay_error, float(np.max(abs(replay-stored[order]))))
            sums[order] += stored[order]
        q16, q24 = arrays[prefix+'quadrature16_precision50'], arrays[prefix+'quadrature24_precision60']
        difference = abs(stored[13]-stored[11]); absolute_order += difference
        # Check already evaluated source bands, including actual angular pairs.
        names = ['group13/mid'] if group == 13 else ([f'group14/mid24_{s}' for s in (1, -1)] if group == 14
                else [f'mid/{group}_{s}' for s in ((1,) if channels[group]['angular_eigenvalue'] == 0 else (1, -1))])
        E = panels[names[0]+'/energies']; mask = (E > .75*lower) & (E < lower)
        if not mask.any() or any(not np.array_equal(E, panels[name+'/energies']) for name in names):
            raise ValueError('paired source-band frequency grids differ')
        expected = control['multiplicity_factor']*sum(panels[name+'/kernels'][mask] for name in names)
        indicator = control['multiplicity_factor']*sum(panels[name+'/arithmetic_indicator'][mask] for name in names)
        series_control = dict(control, coefficients_rho_parallel_sphere=coefficients)
        band13 = evaluate_tail_series(series_control, E[mask])
        band11 = evaluate_tail_series(series_control, E[mask], order=11)
        error = abs(band13-expected)
        allowance = 16*indicator+2*abs(band13-band11)
        residuals = {'coefficient_integral_replay': replay_error,
            'order11_13': float(difference.max()), 'quadrature16_24': float(np.max(abs(q24-q16))),
            'series13_quadrature24': float(np.max(abs(stored[13]-q24))),
            'stored_band_absolute': float(error.max()),
            'stored_band_excess_over_declared_indicator': float(np.maximum(error-allowance, 0).max())}
        bad = {k: v for k, v in residuals.items() if k in tolerances and (not np.isfinite(v) or v > tolerances[k])}
        if bad: failures[str(group)] = bad
        for key, value in residuals.items(): maxima[key] = max(maxima[key], value)
        bound = thermal_bound(channels[group], lower, config)
        with mp.workdps(60): thermal_total += bound
        groups[str(group)] = {'lower': lower, 'actual_angular_signs': [1] if channels[group]['angular_eigenvalue'] == 0 else [1, -1],
            'multiplicity_factor': control['multiplicity_factor'],
            'tail_by_order': {str(n): stored[n].tolist() for n in stored},
            'quadrature16_precision50': q16.tolist(), 'quadrature24_precision60': q24.tolist(),
            'order11_13_absolute_change': difference.tolist(), 'residuals': residuals,
            'stored_source_band': {'panels': names, 'rows': int(mask.sum()), 'energies': [float(E[mask].min()), float(E[mask].max())],
                'allowance': '16*stored arithmetic indicator + 2*local order11/13 change; comparison control, not rigorous error bound'},
            'nonintegrable_coefficients_exactly_zero': True,
            'unfactored_first_derivative_numerical_diagnostic': control['unfactored_first_derivative_numerical_diagnostic'],
            'thermal_source_bound': bound_record(bound), 'status': 'PASS' if not bad else 'OPEN'}
    return {'schema': 'NSC-INCOMING-SOURCE-TAIL-v1', 'accountable_author': 'Douglas Ek',
        'status': 'PASS: all retained source tails in the original Riccati/ad4 approximation; physical remainder OPEN' if not failures else 'OPEN',
        'source_hashes': {p: sha(p) for p in SOURCES}, 'input_hashes': {p: sha(p) for p in INPUTS},
        'parent_finite_source_commit': BASELINE, 'input_payload': finite['payload'],
        'payload': {'path': str(path.relative_to(ROOT)), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest(), 'bytes': path.stat().st_size},
        'locked_inputs': finite['locked_inputs'], 'kernel_order': finite['kernel_order'],
        'groups': groups, 'group_count': len(groups), 'total_tail_by_order': {str(n): sums[n].tolist() for n in sums},
        'total_order11_13_absolute_change': absolute_order.tolist(),
        'finite_source_plus_approximate_vacuum_tail': (np.array(finite['finite_quantum_source_approximant'])+sums[13]).tolist(),
        'thermal_source_bound': dict(bound_record(thermal_total),
            source_trace_norm='2*exp(-pi*E/kappa)+exp(-2*pi*E/(Omega*kappa))',
            scope='same horizon/incoming state correction with current-normalized compression; includes horizon coherence, not vacuum mode error'),
        'regular_endpoint_identity': leading_reference_matching_identity(),
        'maxima': maxima, 'tolerances': tolerances, 'failures': failures,
        'domain': {'surface': 'same rho=1 incoming KS half-density', 'riccati_order': 16, 'reference_order': 4,
            'Taylor_orders': [7, 9, 11, 13], 'Taylor_precision': 50, 'differentiation_chop': False,
            'exact_zero_source': 'source=x^2*g(x), g(0)=0 by the explicit c1/c2/ad1 identity',
            'base_and_refinement_panels': 'no reassembly; only missing tails beyond the existing group endpoints added',
            'state_changed': False, 'old_generators_rerun': False, 'mode_or_scattering_solve': False,
            'history_selected': False, 'error_scope': 'measured truncation/quadrature controls for the fixed approximation, not a rigorous physical remainder'},
        'gate': {'specified_approximation_tail': 'PASS' if not failures else 'OPEN',
            'physical_Riccati_mode_remainder': 'OPEN', 'subgap_source_refinement': 'OPEN',
            'retained_angular_compact_remainder': 'OPEN', 'full_renormalized_source': 'OPEN',
            'incoming_constraints': 'OPEN', 'extended_stationarity': 'OPEN',
            'updated_constraints': None, 'selected_nulls': None, 'metric_timestep': False,
            'Z3': 'OUT OF SCOPE', 'PDF_bumped': False, 'homogeneous_nonexistence': 'preserved; not rerun',
            'Weyl_time_node_diagnostic': 93.54264532195464},
        'comparison': {'fields': 'all', 'float_atol': 3e-24, 'float_rtol': 3e-12,
            'exact': 'hashes, labels, structural zeros, strings, coverage and scope'},
        'reproducer': 'python3 scripts/derive_nsc_incoming_source_tail.py --check'}


def compare(expected, actual, path='$'):
    if type(expected) is not type(actual): raise AssertionError('type mismatch: '+path)
    if isinstance(expected, dict):
        if expected.keys() != actual.keys(): raise AssertionError('keys differ: '+path)
        for key in expected: compare(expected[key], actual[key], path+'/'+key)
    elif isinstance(expected, list):
        if len(expected) != len(actual): raise AssertionError('length mismatch: '+path)
        for index, (left, right) in enumerate(zip(expected, actual)): compare(left, right, path+'/'+str(index))
    elif isinstance(expected, float):
        if not math.isfinite(expected) or not math.isfinite(actual) or not math.isclose(expected, actual, abs_tol=3e-24, rel_tol=3e-12):
            raise AssertionError('number mismatch: '+path)
    elif expected != actual: raise AssertionError('value mismatch: '+path)


def main():
    parser = argparse.ArgumentParser(); parser.add_argument('--prepare', action='store_true'); parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    if args.prepare and args.check: parser.error('choose prepare or check')
    if args.prepare: path = prepare()
    elif args.check:
        receipt = read(OUTPUT); path = ROOT/receipt['payload']['path']
        if sha(receipt['payload']['path']) != receipt['payload']['sha256']: raise ValueError('tail artifact changed')
    else: parser.error('use --prepare or --check')
    result = make_record(path)
    if args.check: compare(read(OUTPUT), result)
    else: (ROOT/OUTPUT).write_text(json.dumps(result, sort_keys=True, indent=2)+'\n')
    print(json.dumps({key: result[key] for key in ('status', 'group_count', 'total_tail_by_order', 'maxima', 'failures')}, indent=2), flush=True)
    if result['failures']: raise SystemExit(1)


if __name__ == '__main__': main()
