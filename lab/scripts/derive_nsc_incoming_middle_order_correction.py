#!/usr/bin/env python3
"""Prepare/replay explicit order24-minus16 group14 finite-middle correction."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from recursive_horizons.nsc_incoming_middle_order_correction import prepare_correction, authenticated_group14, SIGNS
from recursive_horizons.nsc_incoming_middle_bound import INPUT_RECORDS
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes
from recursive_horizons.nsc_incoming_cauchy_jets import incoming_cauchy_jets
from recursive_horizons.nsc_incoming_joint_constraints import source_action_gradient

OUTPUT = 'results/development/nsc-incoming-middle-order-correction.json'
SOURCES = ('src/recursive_horizons/nsc_incoming_middle_order_correction.py',
           'tests/test_nsc_incoming_middle_order_correction.py',
           'scripts/derive_nsc_incoming_middle_order_correction.py',
           'docs/nsc-incoming-middle-order-correction.md')
INPUTS = (*INPUT_RECORDS, 'results/development/nsc-incoming-middle-bound.json',
          'src/recursive_horizons/nsc_incoming_middle_bound.py',
          'src/recursive_horizons/nsc_incoming_source_tail.py',
          'src/recursive_horizons/nsc_incoming_state_moments.py',
          'src/recursive_horizons/nsc_incoming_joint_constraints.py',
          'src/recursive_horizons/nsc_incoming_cauchy_jets.py',
          'scripts/derive_nsc_pg_mid_inputs.py',
          'src/recursive_horizons/nsc_mode_resolved_cauchy_state.py')


def digest(path): return sha256((ROOT/path).read_bytes()).hexdigest()
def signature(): return {p: digest(p) for p in (*SOURCES, *INPUTS)}


def prepare():
    before = signature(); arrays, report = prepare_correction(ROOT)
    if signature() != before: raise ValueError('correction source or input changed during preparation')
    metadata = {'schema': 'NSC-INCOMING-MIDDLE-ORDER-CORRECTION-ARTIFACT-v1',
                'signature': before, 'report': report}
    arrays['metadata_json'] = np.frombuffer(json.dumps(metadata, sort_keys=True, allow_nan=False).encode(), np.uint8)
    raw = deterministic_npz_bytes(arrays); h = sha256(raw).hexdigest()
    path = ROOT/f'results/development/artifacts/nsc-incoming-middle-order-correction.{h}.npz'
    if path.exists() and path.read_bytes() != raw: raise ValueError('content-addressed collision')
    if not path.exists(): path.write_bytes(raw)
    return path


def make_record(path):
    owned = authenticated_group14(ROOT)
    with np.load(path, allow_pickle=False) as data: arrays = {k: data[k].copy() for k in data.files}
    metadata = json.loads(arrays['metadata_json'].tobytes()); report = metadata['report']
    if metadata['signature'] != signature(): raise ValueError('correction producer or dependency changed')
    if report['group'] != 14 or report['old_order'] != 16 or report['new_order'] != 24:
        raise ValueError('explicit group14 order24-minus16 correction required')
    if report['input_payloads'] != owned['input_payloads'] or report['interval'] != owned['interval']:
        raise ValueError('correction baseline or selected interval changed')
    factor = owned['factor']; replay_max = 0.; prefix_max = 0.; measure_max = 0.
    original = np.zeros(4); refined = np.zeros(4); old = np.zeros(4); oldvac = np.zeros(4); thermal = np.zeros(4)
    per_sign = {}
    for sign in SIGNS:
        for label, dps in (('original24', 70), ('refined48', 50), ('refined48', 70)):
            key = f'{sign}/{label}/dps{dps}'; row = report['reports'][key]
            E, weights, kernels = [arrays[key+'/'+k] for k in ('energies', 'weights', 'correction_kernels')]
            if (len(E) != row['rows'] or kernels.shape != (len(E), 4) or np.any(weights <= 0)
                    or np.any(E <= 16) or np.any(E >= 160) or np.count_nonzero(kernels[:, 2])):
                raise ValueError('owned positive finite-band correction with exact zero current required')
            value = factor*np.einsum('n,nv->v', weights, kernels)
            replay_max = max(replay_max, float(np.max(abs(value-row['physical_correction']))))
            prefix_max = max(prefix_max, row['coefficient_prefix_residual'])
            measure_max = max(measure_max, row['integrated_dE_measure_residual'])
            if label == 'original24':
                archive = owned['arrays']; source_prefix = f'group14/mid24_{sign}/'
                if not np.array_equal(E, archive[source_prefix+'energies']) or not np.array_equal(weights, archive[source_prefix+'weights']):
                    raise ValueError('original correction nodes or weights differ from authenticated source')
                original += row['physical_correction']
            elif dps == 70: refined += row['physical_correction']
        weights = arrays[f'{sign}/archived16/weights']
        for key in ('energies', 'weights', 'kernels', 'vacuum_kernels', 'thermal_kernels'):
            if not np.array_equal(arrays[f'{sign}/archived16/'+key], owned['arrays'][f'group14/mid24_{sign}/'+key]):
                raise ValueError('order16 archived values cannot be replaced by this correction')
        contribution = lambda key: factor*np.einsum('n,nv->v', weights, arrays[f'{sign}/archived16/'+key])
        old += contribution('kernels'); oldvac += contribution('vacuum_kernels'); thermal += contribution('thermal_kernels')
        per_sign[str(sign)] = {'original24_correction': report['reports'][f'{sign}/original24/dps70']['physical_correction'],
                              'refined48_correction': report['reports'][f'{sign}/refined48/dps70']['physical_correction']}
    residuals = {'quadrature_contraction_replay': replay_max, 'coefficient16_prefix': prefix_max,
        'old16_vacuum_kernel_vs_archive': max(report['reports'][str(s)]['order16_vacuum_kernel_vs_archive'] for s in SIGNS),
        'source16_decomposition': max(report['reports'][str(s)]['source16_decomposition_residual'] for s in SIGNS),
        'source24_grid_recipe': report['grid_recipe_residual'], 'quadrature_measure': measure_max,
        'quadrature24_to48_physical': max(report['sum_absolute_sign_quadrature_indicator']),
        'precision50_to70_physical': max(report['sum_absolute_sign_precision_indicator'])}
    tolerances = {'quadrature_contraction_replay': 3e-24, 'coefficient16_prefix': 3e-40,
        'old16_vacuum_kernel_vs_archive': 3e-11, 'source16_decomposition': 3e-11,
        'source24_grid_recipe': 3e-13, 'quadrature_measure': 3e-11,
        'quadrature24_to48_physical': 3e-20, 'precision50_to70_physical': 3e-35}
    failures = {k: v for k, v in residuals.items() if not np.isfinite(v) or v > tolerances[k]}
    gradient = source_action_gradient(incoming_cauchy_jets(), refined)
    return {'schema': 'NSC-INCOMING-MIDDLE-ORDER-CORRECTION-v1', 'accountable_author': 'Douglas Ek',
        'status': ('PASS: explicit group14 order24-minus16 numerical correction; physical order24 error OPEN' if not failures else
                   'OPEN: correction numerical controls failed'),
        'group': 14, 'old_order': 16, 'new_order': 24, 'kernel_order': ['rho', 'p_parallel', 'T01', 'p_perp'],
        'interval': report['interval'], 'group_factor_per_angular_sign': factor, 'per_sign': per_sign,
        'old16_middle_source': old.tolist(), 'old16_middle_vacuum': oldvac.tolist(),
        'unchanged_archived_thermal_insertion': thermal.tolist(),
        'original24_node_vacuum_correction': original.tolist(),
        'refined48_node_vacuum_correction': refined.tolist(),
        'additively_corrected_middle_approximation': (old+refined).tolist(),
        'raw_metric_order': ['N', 'beta', 'a', 'r'], 'source_action_gradient_correction': gradient.tolist(),
        'quadrature_indicator_sum_absolute_signs': report['sum_absolute_sign_quadrature_indicator'],
        'precision_indicator_sum_absolute_signs': report['sum_absolute_sign_precision_indicator'],
        'thermal_scattering_stress_error_upper': report['thermal_scattering_stress_error_upper'],
        'physical_order24_defect_bound': None, 'physical_order24_defect_status': 'OPEN',
        'residuals': residuals, 'numerical_tolerances': tolerances, 'failures': failures,
        'scope': {'physical_operator_horizon_state_and_ad4_unchanged': True,
            'old_source_artifact_replaced': False, 'vacuum_only_additive_correction': True,
            'archived_thermal_insertion_changed': False, 'current_correction_exactly_zero': True,
            'quadrature_or_precision_indicator_is_rigorous_bound': False,
            'numerical_order_difference_is_physical_error_bound': False,
            'source_law_refit': False, 'full_source_convergence_claimed': False,
            'field_mode_scattering_or_horizon_solve': False, 'metric_step': False,
            'seed_used_as_state': False, 'Gamma_rest_assigned': False,
            'other_groups_low_subgap_and_existing_tail_changed': False},
        'input_payloads': report['input_payloads'], 'input_hashes': {p: digest(p) for p in INPUTS},
        'source_hashes': {p: digest(p) for p in SOURCES},
        'payload': {'path': str(path.relative_to(ROOT)), 'sha256': sha256(path.read_bytes()).hexdigest(), 'bytes': path.stat().st_size},
        'reproducer': 'python3 scripts/derive_nsc_incoming_middle_order_correction.py --check'}


def main():
    parser = argparse.ArgumentParser(description=__doc__); mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--prepare', action='store_true'); mode.add_argument('--check', action='store_true'); args = parser.parse_args()
    if args.prepare:
        result = make_record(prepare()); (ROOT/OUTPUT).write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False)+'\n')
    else:
        prior = json.loads((ROOT/OUTPUT).read_text()); spec = prior['payload']
        if digest(spec['path']) != spec['sha256']: raise ValueError('correction artifact changed')
        result = make_record(ROOT/spec['path'])
        if result != prior: raise ValueError('record differs from authenticated correction replay')
    print(json.dumps({k: result[k] for k in ('status', 'refined48_node_vacuum_correction', 'source_action_gradient_correction', 'residuals', 'failures')}, indent=2))
    return 1 if result['failures'] else 0


if __name__ == '__main__': raise SystemExit(main())
