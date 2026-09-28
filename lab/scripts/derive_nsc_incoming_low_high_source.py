#!/usr/bin/env python3
"""Prepare/replay one group22 LOW[32,40] source and matching energy certificate."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

import mpmath as mp
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from recursive_horizons.nsc_incoming_low_high_source import (
    prepare_source, prepare_group22_defect, authenticated_low_inputs, SIGNS,
)
from recursive_horizons.nsc_incoming_middle_bound import (
    INPUT_RECORDS, thermal_middle_difference_bound, constraint_action_bound,
)
from recursive_horizons.nsc_incoming_projector_energy_bound import (
    local_cross_product_coefficients, projected_energy_bound,
)
from recursive_horizons.nsc_incoming_vacuum_tail_bound import _precision, _up_float
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes
from recursive_horizons.nsc_incoming_joint_constraints import source_action_gradient
from recursive_horizons.nsc_incoming_cauchy_jets import incoming_cauchy_jets

OUTPUT = 'results/development/nsc-incoming-low-high-source.json'
SOURCES = ('src/recursive_horizons/nsc_incoming_low_high_source.py',
           'tests/test_nsc_incoming_low_high_source.py',
           'scripts/derive_nsc_incoming_low_high_source.py',
           'docs/nsc-incoming-low-high-source.md')
OWNERS = ('src/recursive_horizons/nsc_incoming_centered_order24.py',
          'src/recursive_horizons/nsc_incoming_defect_taylor_bound.py',
          'src/recursive_horizons/nsc_incoming_projector_energy_bound.py',
          'src/recursive_horizons/nsc_incoming_vacuum_tail_bound.py',
          'src/recursive_horizons/nsc_incoming_middle_bound.py',
          'src/recursive_horizons/nsc_incoming_source_tail.py',
          'src/recursive_horizons/nsc_incoming_state_moments.py',
          'src/recursive_horizons/nsc_incoming_joint_constraints.py',
          'src/recursive_horizons/nsc_incoming_cauchy_jets.py',
          'src/recursive_horizons/nsc_mode_resolved_cauchy_state.py')


def digest(path): return sha256((ROOT/path).read_bytes()).hexdigest()
def signature(): return {p: digest(p) for p in (*SOURCES, *INPUT_RECORDS, *OWNERS)}


def prepare():
    before = signature(); arrays, report, owned = prepare_source(ROOT)
    radial = prepare_group22_defect(owned['config'], owned['channel'],
        progress=lambda n: print(f'group22 order24 centered radial enclosure {n}/128', flush=True))
    cross = local_cross_product_coefficients(owned['channel'], order=24)
    if signature() != before: raise ValueError('frozen source or helper changed during preparation')
    metadata = {'schema': 'NSC-INCOMING-LOW-HIGH-SOURCE-ARTIFACT-v1', 'signature': before,
                'report': report, 'radial': radial, 'local_cross': cross}
    arrays['metadata_json'] = np.frombuffer(json.dumps(metadata, sort_keys=True, allow_nan=False).encode(), np.uint8)
    raw = deterministic_npz_bytes(arrays); h = sha256(raw).hexdigest()
    path = ROOT/f'results/development/artifacts/nsc-incoming-low-high-source.{h}.npz'
    if path.exists() and path.read_bytes() != raw: raise ValueError('content-addressed collision')
    if not path.exists(): path.write_bytes(raw)
    return path


def make_record(path):
    owned = authenticated_low_inputs(ROOT)
    with np.load(path, allow_pickle=False) as data: arrays = {k: data[k].copy() for k in data.files}
    metadata = json.loads(arrays['metadata_json'].tobytes()); report = metadata['report']
    if metadata['signature'] != signature(): raise ValueError('source/certificate owner or input changed')
    if report['input_payloads'] != owned['input_payloads'] or report['source_factor_per_sign'] != owned['factor']:
        raise ValueError('same archived source and factor required')
    radial, cross = metadata['radial'], metadata['local_cross']
    if (radial['group'], radial['physical_Riccati_order'], radial['intervals'], radial['centered_remainder_depth']) != (22, 24, 128, 4):
        raise ValueError('one matching group22 order24 centered128/depth4 certificate required')
    physical = projected_energy_bound(owned['channel'], radial, cross, 32., 40.)
    thermal = thermal_middle_difference_bound(owned['channel'], owned['config'], 32., 40.)
    thermal_action = constraint_action_bound(thermal)
    with _precision(40):
        density_error = _up_float(mp.iv.mpf(physical['density_error_upper'])+mp.iv.mpf(thermal[0]))
        lapse_error = _up_float(mp.iv.mpf(physical['lapse_action_error_upper'])+mp.iv.mpf(thermal_action[0]))
    factor = owned['factor']; replay_error = prefix_error = measure_error = 0.
    old = {n: np.zeros(4) for n in (24, 32)}; source_by_grid = {n: np.zeros(4) for n in (24, 32, 48, 64)}
    per_sign = {}
    for sign in SIGNS:
        for label, precision in (('old24', 70), ('old32', 70), ('new48', 70), ('new64', 50), ('new64', 70)):
            key = f'{sign}/{label}/dps{precision}'; row = report['reports'][key]
            E, w, K = [arrays[key+'/'+name] for name in ('energies', 'weights', 'vacuum_kernels')]
            if np.any(w <= 0) or np.any(E <= 32) or np.any(E >= 40) or K.shape != (len(E), 4) or np.count_nonzero(K[:, 2]):
                raise ValueError('positive source24 quadrature with exact vacuum current zero required')
            replay = factor*np.einsum('n,nv->v', w, K)
            replay_error = max(replay_error, float(np.max(abs(replay-row['physical_vacuum_source']))))
            prefix_error = max(prefix_error, row['coefficient16_prefix_residual'])
            measure_error = max(measure_error, row['quadrature_measure_residual'])
            if precision == 70: source_by_grid[int(label[-2:])] += row['physical_vacuum_source']
        for n in (24, 32):
            original = owned['selected'][f'{sign}/old{n}']
            for key, value in original.items():
                if not np.array_equal(arrays[f'{sign}/archived{n}/'+key], value):
                    raise ValueError('archived LOW contribution changed')
                if key in ('energies', 'weights') and not np.array_equal(arrays[f'{sign}/old{n}/dps70/'+key], value):
                    raise ValueError('original source-grid comparison differs from archive')
            old[n] += factor*np.einsum('n,nv->v', original['weights'], original['kernels'])
        per_sign[str(sign)] = report['reports'][f'{sign}/new64/dps70']['physical_vacuum_source']
    new = source_by_grid[64]; delta = new-old[32]
    indicators = report['indicators_sum_absolute_signs']
    residuals = {'quadrature_contraction_replay': replay_error, 'coefficient16_prefix': prefix_error,
                 'original_cell_measure': report['owned_interval_measure_residual'],
                 'source_quadrature_measure': measure_error,
                 **{name: max(values) for name, values in indicators.items()}}
    tolerances = {'quadrature_contraction_replay': 3e-22, 'coefficient16_prefix': 3e-40,
                  'original_cell_measure': 3e-13, 'source_quadrature_measure': 3e-13,
                  'quadrature32_to48': 3e-20, 'quadrature48_to64': 3e-20, 'precision50_to70': 3e-35}
    failures = {k: v for k, v in residuals.items() if not np.isfinite(v) or v > tolerances[k]}
    fits = max(lapse_error, thermal_action[1]) <= 3e-11
    return {'schema': 'NSC-INCOMING-LOW-HIGH-SOURCE-v1', 'accountable_author': 'Douglas Ek',
        'status': ('PASS: group22 LOW[32,40] vacuum/thermal energy budget and numerical source controls; full source OPEN'
                   if fits and not failures else 'OPEN: group22 LOW[32,40] error budget or numerical source controls unresolved'),
        'group': 22, 'numerical_Riccati_order': 24, 'interval': [32., 40.], 'actual_angular_signs': [1, -1],
        'kernel_order': ['rho', 'p_parallel', 'T01', 'p_perp'], 'group_factor_per_sign': factor,
        'archived_LOW24_source': old[24].tolist(), 'archived_LOW32_source': old[32].tolist(),
        'archived_LOW32_minus_LOW24': (old[32]-old[24]).tolist(),
        'new_order24_vacuum_source': new.tolist(), 'new_vacuum_minus_archived_LOW32': delta.tolist(),
        'new_vacuum_minus_archived_LOW24': (new-old[24]).tolist(),
        'new_source_by_positive_quadrature': {str(n): v.tolist() for n, v in source_by_grid.items()},
        'new_source_per_sign': per_sign,
        'raw_metric_order': ['N', 'beta', 'a', 'r'],
        'additive_action_gradient_correction': source_action_gradient(incoming_cauchy_jets(), delta).tolist(),
        'physical_mode_energy_bound': physical,
        'omitted_thermal_stress_upper': thermal,
        'thermal_interpretation': 'physical contribution retained as bounded remainder, not assigned zero; conservatively uses existing2*b(E) bound',
        'physical_error_budget': {'density_error_upper': density_error,
            'lapse_action_error_upper': lapse_error, 'shift_action_error_upper': thermal_action[1],
            'stationarity_tolerance': 3e-11, 'fits_stationarity_tolerance': fits,
            'rigorous_numerical_source_quadrature_error_included': False},
        'numerical_indicators_sum_absolute_signs': indicators,
        'residuals': residuals, 'numerical_tolerances': tolerances, 'failures': failures,
        'radial_certificate': radial, 'local_cross_certificate': cross,
        'scope': {'same_physical_operator_affine_horizon_source_and_ad4': True,
            'physical_state_selected_or_refitted': False, 'old_LOW_source_artifacts_overwritten': False,
            'old_finite_offset_accuracy_certified': False, 'old_difference_identified_as_inverse_noise': False,
            'thermal_source_physically_zero': False, 'source_and_certificate_order_match': True,
            'pressure_mode_error_certified': False, 'full_source_convergence_claimed': False,
            'other_groups_or_intervals_changed': False, 'mode_or_scattering_solve': False,
            'metric_step': False, 'Gamma_rest_assigned': False,
            'seed_used_as_state': False, 'numerical_indicators_are_rigorous_bounds': False},
        'input_payloads': owned['input_payloads'], 'input_hashes': {p: digest(p) for p in (*INPUT_RECORDS, *OWNERS)},
        'source_hashes': {p: digest(p) for p in SOURCES},
        'payload': {'path': str(path.relative_to(ROOT)), 'sha256': sha256(path.read_bytes()).hexdigest(), 'bytes': path.stat().st_size},
        'reproducer': 'python3 scripts/derive_nsc_incoming_low_high_source.py --check'}


def main():
    parser = argparse.ArgumentParser(description=__doc__); mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--prepare', action='store_true'); mode.add_argument('--check', action='store_true'); args = parser.parse_args()
    if args.prepare:
        record = make_record(prepare()); (ROOT/OUTPUT).write_text(json.dumps(record, indent=2, sort_keys=True, allow_nan=False)+'\n')
    else:
        prior = json.loads((ROOT/OUTPUT).read_text()); spec = prior['payload']
        if digest(spec['path']) != spec['sha256']: raise ValueError('source/certificate artifact changed')
        record = make_record(ROOT/spec['path'])
        if record != prior: raise ValueError('record differs from authenticated source/certificate replay')
    print(json.dumps({k: record[k] for k in ('status', 'new_order24_vacuum_source', 'new_vacuum_minus_archived_LOW32',
                                          'physical_error_budget', 'residuals', 'failures')}, indent=2))
    return 1 if record['failures'] else 0


if __name__ == '__main__': raise SystemExit(main())
