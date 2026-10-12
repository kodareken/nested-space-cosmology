#!/usr/bin/env python3
"""Bound one finite-collar endpoint, without propagation or spectral generation."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from recursive_horizons.nsc_incoming_horizon_endpoint_bound import (
    endpoint_bound, frobenius_residual_identity, archived_formula_controls,
)
from recursive_horizons.nsc_incoming_source_quadrature_bound import unpack, float_enclosure
from recursive_horizons.nsc_incoming_vacuum_tail_bound import _precision, _hi

OUTPUT = ROOT/'results/development/nsc-incoming-horizon-endpoint-bound.json'
SOURCES = ('src/recursive_horizons/nsc_incoming_horizon_endpoint_bound.py',
           'scripts/derive_nsc_incoming_horizon_endpoint_bound.py',
           'tests/test_nsc_incoming_horizon_endpoint_bound.py',
           'docs/nsc-incoming-horizon-endpoint-bound.md')
INPUTS = ('results/development/nsc-mode-resolved-cauchy-state.json',
          'results/development/nsc-compact-matched-restart.json',
          'src/recursive_horizons/nsc_paired_horizon_preparation.py',
          'src/recursive_horizons/nsc_unruh_state.py',
          'src/recursive_horizons/nsc_lorentzian.py',
          'src/recursive_horizons/nsc_angular_stress.py',
          'src/recursive_horizons/nsc_incoming_vacuum_tail_bound.py',
          'src/recursive_horizons/nsc_incoming_middle_bound.py',
          'src/recursive_horizons/nsc_incoming_source_quadrature_bound.py',
          'src/recursive_horizons/nsc_pg_massive_modes.py')


def digest(path):
    return sha256((ROOT/path).read_bytes()).hexdigest()


def signature():
    return {path: digest(path) for path in (*SOURCES, *INPUTS)}


def prepare():
    before = signature()
    state = json.loads((ROOT/INPUTS[0]).read_text())
    restart = json.loads((ROOT/INPUTS[1]).read_text())
    # The endpoint inputs are explicit immutable values from these recorded
    # source owners. Nothing in the actual initializer is changed or rerun.
    channel = state['channels'][32]; config = restart['scattering_provenance']['config']
    data = {'signature': before, 'channel': channel, 'config': config,
            'identity': frobenius_residual_identity(), 'certificate': endpoint_bound(channel, config),
            'numerical_controls': archived_formula_controls(channel, config)}
    if signature() != before:
        raise ValueError('endpoint owner or input changed during bounded evaluation')
    raw = json.dumps(data, sort_keys=True, allow_nan=False).encode()
    h = sha256(raw).hexdigest()
    relative = f'results/development/artifacts/nsc-incoming-horizon-endpoint-bound.{h}.json'
    target = ROOT/relative
    if target.exists() and target.read_bytes() != raw:
        raise ValueError('content-addressed collision')
    if not target.exists():
        target.write_bytes(raw)
    return target


def make_record(path):
    data = json.loads(path.read_text())
    if data['signature'] != signature():
        raise ValueError('endpoint producer or input changed')
    if data['identity'] != frobenius_residual_identity():
        raise ValueError('finite-series algebraic identity changed')
    certificate = data['certificate']; bounds = certificate['bounds']
    with _precision(certificate['precision']):
        for value in bounds.values():
            if float_enclosure(unpack(value['binary_interval'])) != value['float_enclosure']:
                raise ValueError('directed bound endpoint serialization changed')
        total = unpack(bounds['total_mathematical_initializer_lapse_error_upper']['binary_interval'])
        passed = _hi(total) <= certificate['contextual_remaining_lapse_budget']
        if passed != certificate['mathematical_endpoint_budget_pass']:
            raise ValueError('endpoint gate changed')
    control = data['numerical_controls']
    if max(control['sampled_working_frame_covariance_difference'], control['sampled_raw_frame_norm_residual']) > 3e-14:
        raise ArithmeticError('phase-free formula no longer matches original initializer')
    return {'schema': 'NSC-INCOMING-HORIZON-ENDPOINT-BOUND-v1', 'accountable_author': 'Douglas Ek',
            'status': ('PASS: mathematical endpoint bias bound; propagation and arithmetic OPEN' if passed else
                       'OPEN: finite-collar endpoint enclosure exceeds remaining lapse budget'),
            'source_hashes': {p: digest(p) for p in SOURCES}, 'input_hashes': {p: digest(p) for p in INPUTS},
            'payload': {'path': str(path.relative_to(ROOT)), 'sha256': digest(path)},
            'identity': data['identity'], 'certificate': certificate, 'numerical_controls': control,
            'numerical_control_tolerance': 3e-14,
            'interpretation': 'upper bounds are not measured stress errors or a non-existence certificate; ten-term truncation is negligible, geometric/chart mismatch requires a better representation',
            'next_connection': 'same frozen source input, exact-profile coordinate and geometric Frobenius exponent; include first omitted delta^(3/2) profile correction before any validated propagation',
            'scope': {'archived_initial_modes_changed': False, 'source_kappa_changed': False,
                      'source_occupation_or_coherence_changed': False,
                      'radial_mode_or_source_solves': 0, 'metric_timestep': False,
                      'physical_initial_data_selected': False, 'constraint_solution': False,
                      'push_or_PDF': False},
            'reproducer': 'python3 scripts/derive_nsc_incoming_horizon_endpoint_bound.py --check'}


def main():
    parser = argparse.ArgumentParser()
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument('--write', action='store_true')
    modes.add_argument('--check', action='store_true')
    args = parser.parse_args()
    if args.write:
        if OUTPUT.exists():
            raise FileExistsError('endpoint receipt is not overwritten')
        result = make_record(prepare())
        OUTPUT.write_text(json.dumps(result, sort_keys=True, indent=2, allow_nan=False)+'\n')
    else:
        old = json.loads(OUTPUT.read_text()); target = ROOT/old['payload']['path']
        if digest(target) != old['payload']['sha256']:
            raise ValueError('endpoint artifact bytes changed')
        result = make_record(target)
        if result != old:
            raise ValueError('endpoint record replay differs')
    print(json.dumps({'status': result['status'], 'lapse_bound': result['certificate']['bounds']['total_mathematical_initializer_lapse_error_upper']['float_enclosure'],
                      'controls': result['numerical_controls']}, indent=2))


if __name__ == '__main__':
    main()
