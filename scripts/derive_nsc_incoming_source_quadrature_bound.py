#!/usr/bin/env python3
"""Group6 quadrature pilot; interval source algebra only, no radial solves."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from recursive_horizons.nsc_incoming_source_quadrature_bound import (
    certify_source_window, gauss_rule, replay_window,
)
from recursive_horizons.nsc_incoming_window_refinement import authenticated_window

OUTPUT = ROOT/'results/development/nsc-incoming-source-quadrature-bound.json'
SOURCES = ('src/recursive_horizons/nsc_incoming_source_quadrature_bound.py',
           'scripts/derive_nsc_incoming_source_quadrature_bound.py',
           'tests/test_nsc_incoming_source_quadrature_bound.py',
           'docs/nsc-incoming-source-quadrature-bound.md')
INPUTS = ('results/development/nsc-incoming-spectral-source.json',
          'results/development/nsc-pg-retained-covariance.json',
          'results/development/nsc-mode-resolved-cauchy-state.json',
          'results/development/nsc-compact-matched-restart.json',
          'src/recursive_horizons/nsc_incoming_window_refinement.py',
          'src/recursive_horizons/nsc_incoming_source_tail.py',
          'src/recursive_horizons/nsc_incoming_source_assembly.py',
          'src/recursive_horizons/nsc_incoming_middle_order_correction.py',
          'src/recursive_horizons/nsc_incoming_projector_energy_bound.py',
          'src/recursive_horizons/nsc_incoming_middle_bound.py',
          'src/recursive_horizons/nsc_incoming_vacuum_tail_bound.py',
          'src/recursive_horizons/nsc_incoming_state_moments.py',
          'src/recursive_horizons/nsc_compact_ctp_neck.py')
CASES = (('group6_low', 'low_vacuum', 32, 40),
         ('group6_middle_correction', 'middle_correction', 40, 160))


def digest(path):
    return sha256((ROOT/path).read_bytes()).hexdigest()


def signature():
    return {p: digest(p) for p in (*SOURCES, *INPUTS)}


def prepare():
    before = signature()
    rule = gauss_rule(48, 80)
    windows = {}
    for name, kind, left, right in CASES:
        owned = authenticated_window(ROOT, 6, kind, left, right)
        data = certify_source_window(owned['channel'], kind, left, right, rule=rule,
                                     progress=lambda row: print(json.dumps(row), flush=True))
        data['selected_panels'] = owned['descriptors']
        data['input_payloads'] = owned['input_payloads']
        windows[name] = data
    if before != signature():
        raise ValueError('quadrature producer or input changed during preparation')
    raw = json.dumps({'schema': 'NSC-INCOMING-SOURCE-QUADRATURE-ARTIFACT-v1',
                      'signature': before, 'windows': windows}, sort_keys=True, allow_nan=False).encode()
    h = sha256(raw).hexdigest()
    target = ROOT/f'results/development/artifacts/nsc-incoming-source-quadrature-bound.{h}.json'
    if target.exists() and target.read_bytes() != raw:
        raise ValueError('content-addressed collision')
    if not target.exists():
        target.write_bytes(raw)
    return target


def make_record(path):
    payload = json.loads(path.read_text())
    if payload['signature'] != signature():
        raise ValueError('quadrature owner or dependency changed')
    if set(payload['windows']) != {row[0] for row in CASES}:
        raise ValueError('exact two pilot windows required')
    reports = {}
    for name, kind, left, right in CASES:
        item = payload['windows'][name]
        owned = authenticated_window(ROOT, 6, kind, left, right)
        if (item['channel'] != owned['channel'] or item['selected_panels'] != owned['descriptors']
                or item['input_payloads'] != owned['input_payloads'] or item['group'] != 6
                or item['kind'] != kind or item['interval'] != [left, right] or item['source_order'] != 24):
            raise ValueError('selected archived source window changed')
        expected = {(sign, start, start+8) for sign in owned['signs'] for start in range(left, right, 8)}
        actual = [(cell['sign'], cell['left'], cell['right']) for cell in item['cells']]
        if len(actual) != len(expected) or set(actual) != expected:
            raise ValueError('source cells missing, repeated or outside the owned window')
        reports[name] = {'group': 6, 'kind': kind, 'interval': [left, right],
                         'source_order': 24, 'cells_per_sign': (right-left)//8,
                         'selected_panels': owned['descriptors'], **replay_window(item)}
    passed = all(item['numerical_accuracy_pass'] for item in reports.values())
    return {'schema': 'NSC-INCOMING-SOURCE-QUADRATURE-BOUND-v1', 'accountable_author': 'Douglas Ek',
            'status': 'PASS: pilot source quadrature and arithmetic enclosed; physical source OPEN' if passed else 'OPEN: pilot numerical enclosure exceeds tolerance',
            'source_hashes': {p: digest(p) for p in SOURCES}, 'input_hashes': {p: digest(p) for p in INPUTS},
            'payload': {'path': str(path.relative_to(ROOT)), 'sha256': digest(path)},
            'windows': reports, 'kernel_order': ['rho', 'p_parallel', 'T01', 'p_perp'],
            'method': {'Gauss_points_per_cell': 48, 'directed_precision': 80,
                       'cell_halfwidth': 4, 'analytic_disk_radius': 8,
                       'analytic_bound': '4*h*M*(h/R)^(2*n)/(1-h/R)',
                       'circle_arcs': 32, 'exact_rule': 'certified Legendre root brackets and interval derivative weights',
                       'parameter_policy': 'exact definitions and stored binary labels/factor enclosed together',
                       'roots_source': 'https://dlmf.nist.gov/3.5#v',
                       'analytic_source': 'https://dlmf.nist.gov/1.9#iv'},
            'scope': {'radial_recurrences_rerun': 0, 'modes_or_scattering_solved': 0,
                      'point_Riccati_coefficients_enclosed': True,
                      'rigorous_quadrature_and_arithmetic': True,
                      'physical_projector_thermal_errors_included': False,
                      'middle_certificate_covers': 'order24-minus16 correction only; unchanged order16 baseline quadrature remains separate',
                      'source_arrays_or_old_records_changed': False,
                      'physical_IV_selected': False, 'constraint_root': False,
                      'metric_timestep': False, 'scales_or_source_law_changed': False,
                      'push_or_PDF': False},
            'reproducer': 'python3 scripts/derive_nsc_incoming_source_quadrature_bound.py --check'}


def main():
    parser = argparse.ArgumentParser()
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument('--prepare', action='store_true')
    modes.add_argument('--check', action='store_true')
    args = parser.parse_args()
    if args.prepare:
        if OUTPUT.exists():
            raise FileExistsError('existing quadrature record is not overwritten')
        result = make_record(prepare())
        OUTPUT.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False)+'\n')
    else:
        old = json.loads(OUTPUT.read_text())
        target = ROOT/old['payload']['path']
        if digest(target) != old['payload']['sha256']:
            raise ValueError('quadrature payload bytes changed')
        result = make_record(target)
        if result != old:
            raise ValueError('quadrature record replay differs')
    print(json.dumps({'status': result['status'], 'lapse_error_upper':
                      {name: item['lapse_error_upper'] for name, item in result['windows'].items()}}, indent=2))


if __name__ == '__main__':
    main()
