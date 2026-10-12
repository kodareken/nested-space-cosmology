#!/usr/bin/env python3
"""Two source-window updates from saved radial certificates; no radial runs."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

import mpmath as mp
import numpy as np

ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT/'src'))
from recursive_horizons.nsc_incoming_window_refinement import authenticated_window, prepare_window, bind_window_certificate
from recursive_horizons.nsc_incoming_middle_bound import thermal_middle_difference_bound, constraint_action_bound
from recursive_horizons.nsc_incoming_vacuum_tail_bound import _precision, _up_float
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes
from recursive_horizons.nsc_incoming_joint_constraints import source_action_gradient
from recursive_horizons.nsc_incoming_cauchy_jets import incoming_cauchy_jets

OUTPUT = 'results/development/nsc-incoming-window-reuse.json'
CASES = (
    ('group12_low', 12, 'low_vacuum', 32., 40., 'results/development/nsc-incoming-retained-order24-bound.json'),
    ('group22_middle', 22, 'middle_correction', 40., 160., 'results/development/nsc-incoming-low-high-source.json'),
)
SOURCES = ('src/recursive_horizons/nsc_incoming_window_refinement.py',
           'tests/test_nsc_incoming_window_refinement.py',
           'scripts/derive_nsc_incoming_window_reuse.py', 'docs/nsc-incoming-window-reuse.md')
INPUTS = tuple(c[-1] for c in CASES)+(
    'results/development/nsc-incoming-spectral-source.json',
    'results/development/nsc-pg-retained-covariance.json',
    'results/development/nsc-mode-resolved-cauchy-state.json',
    'results/development/nsc-compact-matched-restart.json',
    'src/recursive_horizons/nsc_incoming_source_assembly.py',
    'src/recursive_horizons/nsc_incoming_source_tail.py',
    'src/recursive_horizons/nsc_incoming_middle_order_correction.py',
    'src/recursive_horizons/nsc_incoming_state_moments.py',
    'src/recursive_horizons/nsc_incoming_projector_energy_bound.py',
    'src/recursive_horizons/nsc_incoming_middle_bound.py',
    'src/recursive_horizons/nsc_incoming_vacuum_tail_bound.py',
    'src/recursive_horizons/nsc_incoming_joint_constraints.py',
    'src/recursive_horizons/nsc_incoming_cauchy_jets.py',
    'src/recursive_horizons/nsc_mode_resolved_cauchy_state.py',
)


def digest(path): return sha256((ROOT/path).read_bytes()).hexdigest()
def signature(): return {p: digest(p) for p in (*SOURCES, *INPUTS)}


def certificate(path):
    record = json.loads((ROOT/path).read_text())
    for field in ('source_hashes', 'input_hashes'):
        for p, expected in record[field].items():
            if digest(p) != expected: raise ValueError('reused certificate owner/input changed: '+p)
    spec = record['payload']
    if digest(spec['path']) != spec['sha256']: raise ValueError('reused coefficient artifact changed')
    if spec['path'].endswith('.json'): payload = json.loads((ROOT/spec['path']).read_text())
    else:
        with np.load(ROOT/spec['path'], allow_pickle=False) as a: payload = json.loads(a['metadata_json'].tobytes())
    radial, cross = payload['radial'], payload['local_cross']
    if radial != record.get('radial', record.get('radial_certificate')) or cross != record.get('local_cross', record.get('local_cross_certificate')):
        raise ValueError('certificate record differs from authenticated coefficient payload')
    return radial, cross, spec


def prepare():
    before = signature(); arrays = {}; windows = {}
    for name, group, kind, left, right, prior in CASES:
        radial, cross, spec = certificate(prior)
        values, report, owned = prepare_window(ROOT, group, kind, left, right)
        arrays.update({name+'/'+k: v for k, v in values.items()})
        windows[name] = {'report': report, 'certificate_record': prior, 'certificate_payload': spec,
                         'binding': bind_window_certificate(owned['channel'], radial, cross, left, right)}
    if signature() != before: raise ValueError('window source or reused dependency changed during preparation')
    arrays['metadata_json'] = np.frombuffer(json.dumps({'schema': 'NSC-INCOMING-WINDOW-REUSE-ARTIFACT-v1',
        'signature': before, 'windows': windows}, sort_keys=True, allow_nan=False).encode(), np.uint8)
    raw = deterministic_npz_bytes(arrays); h = sha256(raw).hexdigest()
    path = ROOT/f'results/development/artifacts/nsc-incoming-window-reuse.{h}.npz'
    if path.exists() and path.read_bytes() != raw: raise ValueError('content-addressed collision')
    if not path.exists(): path.write_bytes(raw)
    return path


def make_record(path):
    with np.load(path, allow_pickle=False) as a: arrays = {k: a[k].copy() for k in a.files}
    metadata = json.loads(arrays['metadata_json'].tobytes())
    if metadata['signature'] != signature(): raise ValueError('window producer or input changed')
    if set(metadata['windows']) != {c[0] for c in CASES}: raise ValueError('exactly the two declared windows required')
    windows = {}; all_failures = {}
    for name, group, kind, left, right, prior in CASES:
        item = metadata['windows'][name]; report = item['report']; owned = authenticated_window(ROOT, group, kind, left, right)
        radial, cross, spec = certificate(prior); binding = bind_window_certificate(owned['channel'], radial, cross, left, right)
        if item['binding'] != binding or item['certificate_payload'] != spec or item['certificate_record'] != prior:
            raise ValueError('fresh interval view differs from immutable coefficient replay')
        if (report['group'], report['kind'], report['interval'], report['source_order']) != (group, kind, [left, right], 24):
            raise ValueError('source window or numerical order changed')
        if report['selected_panels'] != owned['descriptors'] or report['input_payloads'] != owned['input_payloads']:
            raise ValueError('selected archived source panels changed')
        old = np.zeros(4); old_thermal = np.zeros(4); numerical = np.zeros(4); original = np.zeros(4)
        maximum = prefix = measure = convention = 0.
        for sign in owned['signs']:
            selected = owned['selected'][str(sign)]
            for key, value in selected.items():
                if not np.array_equal(arrays[f'{name}/{sign}/archived/{key}'], value): raise ValueError('old source values changed')
            old += owned['factor']*np.einsum('n,nv->v', selected['weights'], selected['kernels'])
            if kind == 'middle_correction': old_thermal += owned['factor']*np.einsum('n,nv->v', selected['weights'], selected['thermal_kernels'])
            for label, dps in (('original', 70), ('refined48', 70), ('refined64', 50), ('refined64', 70)):
                key = f'{sign}/{label}/dps{dps}'; row = report['reports'][key]
                E, w, K = [arrays[f'{name}/{key}/'+k] for k in ('energies', 'weights', 'kernels')]
                if np.any(w <= 0) or np.any(E <= left) or np.any(E >= right) or K.shape != (len(E), 4) or np.count_nonzero(K[:, 2]):
                    raise ValueError('positive source quadrature and exact vacuum-current identity required')
                if label == 'original' and (not np.array_equal(E, selected['energies']) or not np.array_equal(w, selected['weights'])):
                    raise ValueError('old source quadrature was replaced')
                replay = owned['factor']*np.einsum('n,nv->v', w, K)
                maximum = max(maximum, float(np.max(abs(replay-row['physical_value']))))
                prefix = max(prefix, row['prefix16_residual']); measure = max(measure, row['measure_residual'])
                if label == 'refined64' and dps == 70: numerical += row['physical_value']
                if label == 'original': original += row['physical_value']
            value = report['reports'][str(sign)]['old16_source_convention_residual']
            if value is not None: convention = max(convention, value)
        delta = numerical-old if kind == 'low_vacuum' else numerical
        updated = numerical if kind == 'low_vacuum' else old+numerical
        thermal = thermal_middle_difference_bound(owned['channel'], owned['config'], left, right)
        thermal_action = constraint_action_bound(thermal); energy = binding['energy']
        with _precision(40):
            total_density = _up_float(mp.iv.mpf(energy['density_error_upper'])+mp.iv.mpf(thermal[0]))
            total_lapse = _up_float(mp.iv.mpf(energy['lapse_action_error_upper'])+mp.iv.mpf(thermal_action[0]))
        residuals = {'quadrature_replay': maximum, 'prefix16': prefix, 'quadrature_measure': measure,
                     'original_measure': report['original_measure_residual'], 'middle_source16_convention': convention,
                     **{k: max(v) for k, v in report['indicators_sum_absolute_signs'].items()}}
        tolerances = {'quadrature_replay': 3e-22, 'prefix16': 3e-40, 'quadrature_measure': 3e-11,
                      'original_measure': 3e-11, 'middle_source16_convention': 3e-11,
                      'original_to48': 3e-20, 'quadrature48_to64': 3e-20, 'precision50_to70': 3e-35}
        failures = {k: v for k, v in residuals.items() if not np.isfinite(v) or v > tolerances[k]}
        if failures: all_failures[name] = failures
        fits = max(total_lapse, thermal_action[1]) <= 3e-11
        windows[name] = {'group': group, 'kind': kind, 'interval': [left, right], 'source_order': 24,
            'selected_panels': owned['descriptors'], 'factor_per_sign': owned['factor'],
            'old_source': old.tolist(), 'new_numerical_value': numerical.tolist(),
            'original_grid_numerical_value': original.tolist(), 'explicit_source_delta': delta.tolist(),
            'updated_source_approximant': updated.tolist(),
            'retained_archived_thermal': old_thermal.tolist() if kind == 'middle_correction' else None,
            'thermal_policy': ('old insertion retained, difference to exact thermal bounded' if kind == 'middle_correction' else
                               'physical thermal retained as nonzero bounded remainder'),
            'thermal_scattering_stress_upper': thermal,
            'additive_action_gradient': source_action_gradient(incoming_cauchy_jets(), delta).tolist(),
            'certificate_record': prior, 'certificate_payload': spec, 'interval_binding': binding,
            'physical_error_budget': {'density_error_upper': total_density, 'lapse_action_error_upper': total_lapse,
                'shift_action_error_upper': thermal_action[1], 'stationarity_tolerance': 3e-11,
                'fits_component_tolerance': fits, 'rigorous_source_quadrature_error_included': False},
            'numerical_indicators_sum_absolute_signs': report['indicators_sum_absolute_signs'],
            'residuals': residuals, 'numerical_tolerances': tolerances, 'failures': failures,
            'status': 'PASS: component physical budget and numerical controls; full source OPEN' if fits and not failures else 'OPEN'}
    return {'schema': 'NSC-INCOMING-WINDOW-REUSE-v1', 'accountable_author': 'Douglas Ek',
        'status': 'PASS: two explicit source-window updates; full source OPEN' if not all_failures and all(w['physical_error_budget']['fits_component_tolerance'] for w in windows.values()) else 'OPEN',
        'kernel_order': ['rho', 'p_parallel', 'T01', 'p_perp'], 'raw_metric_order': ['N', 'beta', 'a', 'r'],
        'windows': windows, 'failures': all_failures,
        'scope': {'radial_recurrences_rerun': 0, 'modes_or_scattering_solved': 0, 'cases_calculated': 2,
            'original_certificates_or_source_arrays_changed': False, 'physical_state_action_scales_changed': False,
            'thermal_source_physically_zero': False, 'rigorous_source_quadrature_error_certified': False,
            'pressure_mode_error_certified': False, 'other_groups_calculated': False,
            'full_source_converged': False, 'metric_step': False, 'Gamma_rest_assigned': False},
        'source_hashes': {p: digest(p) for p in SOURCES}, 'input_hashes': {p: digest(p) for p in INPUTS},
        'payload': {'path': str(path.relative_to(ROOT)), 'sha256': sha256(path.read_bytes()).hexdigest(), 'bytes': path.stat().st_size},
        'reproducer': 'python3 scripts/derive_nsc_incoming_window_reuse.py --check'}


def main():
    parser = argparse.ArgumentParser(); mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--prepare', action='store_true'); mode.add_argument('--check', action='store_true'); args = parser.parse_args()
    if args.prepare:
        record = make_record(prepare()); (ROOT/OUTPUT).write_text(json.dumps(record, indent=2, sort_keys=True, allow_nan=False)+'\n')
    else:
        prior = json.loads((ROOT/OUTPUT).read_text()); spec = prior['payload']
        if digest(spec['path']) != spec['sha256']: raise ValueError('window source artifact changed')
        record = make_record(ROOT/spec['path'])
        if record != prior: raise ValueError('source-window record differs from authenticated replay')
    print(json.dumps({'status': record['status'], 'windows': {k: {name: v[name] for name in
        ('explicit_source_delta', 'physical_error_budget', 'residuals')} for k, v in record['windows'].items()}, 'failures': record['failures']}, indent=2))
    return 1 if record['failures'] else 0


if __name__ == '__main__': raise SystemExit(main())
