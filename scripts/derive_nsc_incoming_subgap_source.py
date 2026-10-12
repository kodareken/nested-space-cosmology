#!/usr/bin/env python3
"""Read-only numerical replay of the prepared group13 incoming source panel.

Neither mode preparation nor any horizon, scattering or history solve is
called. --write creates the new result record; --check compares every field.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT/'src'))
from recursive_horizons.nsc_incoming_subgap_source import read_prepared_incoming_subgap, recompose_incoming_subgap


OUTPUT = 'results/development/nsc-incoming-subgap-source.json'
DIGEST = '97054f5d41d4fc9711b111d0482760f3fe59d0eee8b092a28603da6c3cad8817'
PAYLOAD = 'results/development/artifacts/nsc-incoming-subgap-source.'+DIGEST+'.npz'
SOURCES = ('src/recursive_horizons/nsc_incoming_subgap_source.py',
           'tests/test_nsc_incoming_subgap_source.py',
           'scripts/derive_nsc_incoming_subgap_source.py', 'docs/nsc-incoming-subgap-source.md')
CURRENT_IDENTITY = 'results/development/nsc-incoming-constraint-gate.json'


def sha(path): return hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
def read(path): return json.loads((ROOT/path).read_text())
def encoded(value):
    value = np.asarray(value)
    return {'real': value.real.tolist(), 'imag': value.imag.tolist()} if np.iscomplexobj(value) else value.tolist()


def make_record():
    if sha(PAYLOAD) != DIGEST: raise ValueError('pinned incoming subgap artifact changed')
    arrays = read_prepared_incoming_subgap(ROOT, ROOT/PAYLOAD)
    meta = json.loads(arrays['metadata_json'].tobytes())
    report = recompose_incoming_subgap(arrays)
    identities = read(CURRENT_IDENTITY)['exact_identities']
    if any(value != '0' for value in identities.values()): raise ValueError('inherited sewing-current identity changed')
    maxima = dict(report['diagnostics'])
    maxima.update({name: float(np.max(value)) for name, value in report['error_indicators'].items()})
    maxima['integrated_T01_zero_identity'] = float(abs(report['physical_group13_contribution'][2]))
    tolerances = {name: 3e-11 for name in maxima}
    # Raw analytic bilinears and physical integrated stress have distinct
    # scales and gates. These tolerances were declared in the focused tests.
    tolerances.update(complex_interpolation=3e-10, real8_bilinear_interpolation=3e-10)
    failures = {key: value for key, value in maxima.items() if not np.isfinite(value) or value > tolerances[key]}
    if report['reflection_samples_used'] != 360: failures['reflection_inventory'] = report['reflection_samples_used']
    if report['short_stationary_solves'] != 19: failures['preparation_solve_budget'] = report['short_stationary_solves']
    tables = {name: {key: encoded(value) if isinstance(value, np.ndarray) else value for key, value in row.items()}
              for name, row in report['tables'].items()}
    return {
        'schema': 'NSC-INCOMING-GROUP13-SUBGAP-SOURCE-v1', 'accountable_author': 'Douglas Ek',
        'status': 'PASS: refined group13 incoming local subgap source panel' if not failures else 'OPEN',
        'source_hashes': {p: sha(p) for p in SOURCES},
        'prepared_dependency_hashes': meta['source_hashes'],
        'input_payloads': meta['input_payloads'],
        'current_identity_pointer': {'path': CURRENT_IDENTITY, 'sha256': sha(CURRENT_IDENTITY),
                                     'exact_identities': identities, 'old_generator_rerun': False},
        'payload': {'path': PAYLOAD, 'sha256': DIGEST, 'bytes': (ROOT/PAYLOAD).stat().st_size},
        'source': {'kernel_order': list(report['kernel_order']), 'group_factor': meta['group_factor'],
                   'raw_positive_energy_integral': report['raw_integral'].tolist(),
                   'refined_physical_contribution': report['physical_group13_contribution'].tolist(),
                   'coarse8_positive_energy_integral': report['coarse8_raw_integral'].tolist(),
                   'coarse8_physical_contribution': report['coarse8_physical_contribution'].tolist(),
                   'refined_minus_coarse8': report['refined_minus_coarse8'].tolist(),
                   'exact_physical_subgap_T01': 0.,
                   'raw_T01_interpretation': 'nonzero floating-point zero-identity residual, not a physical flux'},
        'error_indicators_physical': {k: v.tolist() for k, v in report['error_indicators'].items()},
        'diagnostics': report['diagnostics'], 'maxima': maxima, 'tolerances': tolerances, 'failures': failures,
        'quadrature_tables': tables,
        'domain': {'group': 13, 'mass': meta['mass'], 'angular': 0., 'positive_energy_panel': [1., meta['mass']],
                   'incoming_surface_rho': 1., 'canonical_frame': 'exact owned trace inverse and exp(i E S(1))',
                   'clock_shift': meta['clock_shift'], 'state': 'unchanged C_H; closed incoming infinity channel',
                   'basis_order': ['v', 'w'], 'physical_real_columns': '(R v,w,0)',
                   'source_law': 'd(Bvv-Bww)+2 Re[-i s R Bwv]+Tr(V)/2-Tr(Pad4 V)',
                   'normalization': 'raw dE integral times incoming_group_factor; no history-action -1/pi',
                   'interpolation': '17 archived real unsewn bilinear nodes, nested 9 control; never interpolate reflection',
                   'complex_control': 'separately archived conjugate-frequency bases; analytic bilinear, never a Gaussian state',
                   'reference_subtraction': 'existing incoming ad4 evaluated only on the real axis',
                   'current_identity': 'Tr(S Csrc Sdagger)-1=T(n_in-f_H); T=0 for the whole open subgap interval',
                   'error_scope': 'measured interpolation/quadrature/deformation controls, not a rigorous uniform remainder bound'},
        'computation': {'preparation_short_stationary_solves_0_to_1': 19,
                        'preparation_function_evaluations': meta['function_evaluations'],
                        'replay_stationary_solves': 0, 'history_horizon_scattering_solves': 0,
                        'archived_contour_reflections_used': report['reflection_samples_used'],
                        'old_generators_rerun': False},
        'gate': {'group13_local_subgap_refinement': 'PASS' if not failures else 'OPEN',
                 'other_37_compact_signed_family_refinements': 'OPEN: not evaluated here',
                 'full_33_group_source_convergence': 'OPEN', 'full_stress_tensor': None,
                 'local_restoration_added': False, 'compact_induced_added': False,
                 'state_changed': False, 'constraints_solved': False, 'metric_evolution': False,
                 'PDF_or_public_release_changed': False},
        'comparison': {'fields': 'all', 'float_atol': 3e-13, 'float_rtol': 3e-12,
                       'exact': 'hashes, labels, structure, counting and scope'},
        'reproducer': 'python3 scripts/derive_nsc_incoming_subgap_source.py --check',
    }


def compare(expected, actual, path='$'):
    if type(expected) is not type(actual): raise AssertionError('type mismatch '+path)
    if isinstance(expected, dict):
        if expected.keys() != actual.keys(): raise AssertionError('keys differ '+path)
        for key in expected: compare(expected[key], actual[key], path+'.'+key)
    elif isinstance(expected, list):
        if len(expected) != len(actual): raise AssertionError('length differs '+path)
        for i, (left, right) in enumerate(zip(expected, actual)): compare(left, right, path+f'[{i}]')
    elif isinstance(expected, float):
        if not np.isfinite(actual) or not np.isclose(expected, actual, atol=3e-13, rtol=3e-12):
            raise AssertionError('numeric mismatch '+path)
    elif expected != actual: raise AssertionError('value mismatch '+path)


def main():
    parser = argparse.ArgumentParser(); mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write', action='store_true'); mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    if args.write and (ROOT/OUTPUT).exists(): raise ValueError('record already exists; use --check')
    result = make_record()
    if result['failures']: raise ArithmeticError(json.dumps(result['failures']))
    if args.check: compare(read(OUTPUT), result)
    else: (ROOT/OUTPUT).write_text(json.dumps(result, sort_keys=True, indent=2, allow_nan=False)+'\n')
    print(json.dumps({k: result[k] for k in ('status', 'source', 'error_indicators_physical', 'maxima', 'failures')}, indent=2), flush=True)


if __name__ == '__main__': main()
