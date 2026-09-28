#!/usr/bin/env python3
"""v3 control for finite-history UV coefficients; frozen v2 bytes are historical.

--record writes the v3 successor only. The v2 JSON is never overwritten.
"""
import argparse
import json
import subprocess
import time
from hashlib import sha256
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))

from recursive_horizons.nsc_ks_finite_history_uv_coefficients import (
    ARCHIVE_RHO_UP,
    CURRENT_HISTORY_IDENTITY_PREFIX,
    V2_RECORD,
    coefficient_identities,
    e_minus2_current_bookkeeping,
    e_minus2_sigma_contraction,
    finite_history_uv_coefficient_report,
    integrate_imag_major_A2_z,
    named_gaps,
    original_angular_pair_ledger,
    paired_leading_e_minus2_weight,
)
from recursive_horizons.nsc_ks_profile_identity import profile_identity
from recursive_horizons.nsc_local_incoming_family import LocalIncomingFamily

OUTPUT_V2 = ROOT / V2_RECORD
OUTPUT = ROOT / 'results/development/nsc-ks-current-uv-coefficients-v3.json'
HISTORY = ROOT / 'results/development/nsc-ks-gate-history-lm-broyden.json'
OWNERS = (
    'src/recursive_horizons/nsc_ks_finite_history_uv_coefficients.py',
    'docs/nsc-ks-current-uv-coefficients-v2.md',
    'scripts/derive_nsc_ks_current_uv_coefficient_control_v2.py',
    'tests/test_nsc_ks_finite_history_uv_coefficients.py',
)
INPUTS = (
    'docs/nsc-source-cutoff-bridge.md',
    'src/recursive_horizons/nsc_incoming_fixed_transfer.py',
    'src/recursive_horizons/nsc_ks_linear_uv_contractions.py',
    'src/recursive_horizons/nsc_ks_chebyshev_jet_bound.py',
    'src/recursive_horizons/nsc_common_subtracted_ks_source.py',
    'results/development/nsc-ks-gate-history-lm-broyden.json',
    'results/development/nsc-ks-source-inventory.json',
    'results/development/nsc-ks-cutoff-bridge-control.json',
    V2_RECORD,
)
PILOT_MASS = 1.0
PILOT_ANGULAR = float(np.sqrt(5.0))
PILOT_GAUSS = 48
FROZEN_V2_SHA256 = 'de8e67e4cd89eb0331aff4e18d577220b83ddca94c0c6bac6d0c5a0fe32585df'


def _digest(path):
    return sha256((ROOT / path).read_bytes()).hexdigest()


def current_history_pilot():
    payload = json.loads(HISTORY.read_text())
    family = LocalIncomingFamily(np.array(payload['history']['coefficients']))
    identity = profile_identity(family)
    if not identity.startswith(CURRENT_HISTORY_IDENTITY_PREFIX):
        raise ValueError('expected live history 0b0e4ced')
    z = family.collocation_nodes(3)
    started = time.process_time()
    result = integrate_imag_major_A2_z(
        family, z, mass=PILOT_MASS, angular=PILOT_ANGULAR,
        rho_up=ARCHIVE_RHO_UP, gauss_nodes=PILOT_GAUSS)
    cpu = time.process_time() - started
    finer = integrate_imag_major_A2_z(
        family, z, mass=PILOT_MASS, angular=PILOT_ANGULAR,
        rho_up=ARCHIVE_RHO_UP, gauss_nodes=64)
    reference = integrate_imag_major_A2_z(
        LocalIncomingFamily(np.zeros((2, 32))), z, mass=PILOT_MASS,
        angular=PILOT_ANGULAR, rho_up=ARCHIVE_RHO_UP, gauss_nodes=8)
    if cpu >= 120.0:
        raise TimeoutError('geometry pilot exceeded the 120 CPU-second target')
    if not np.array_equal(reference['h_z'], np.zeros((2, 3))):
        raise ArithmeticError('reference history must give vanishing h_z')
    return {
        'profile_identity': identity,
        'z': [float(v) for v in z],
        'mass': PILOT_MASS,
        'angular': PILOT_ANGULAR,
        'rho_up': ARCHIVE_RHO_UP,
        'gauss_nodes': PILOT_GAUSS,
        'h_z': result['h_z'].tolist(),
        'source_signs': list(result['source_signs']),
        'quadrature_n_nodes': result['quadrature']['n_nodes'],
        'certified_quadrature_error_bound': None,
        'gauss_48_minus_64_max_abs': float(np.max(np.abs(result['h_z'] - finer['h_z']))),
        'reference_h_z_max_abs': float(np.max(np.abs(reference['h_z']))),
        'cpu_seconds': cpu,
        'cpu_target_seconds': 120.0,
        'is_full_e_minus2_coefficient': False,
        'numerical_C_M': None,
        'field_or_source_evolutions': 0,
    }


def _historical_commit():
    return subprocess.check_output(
        ['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()


def record():
    if not OUTPUT_V2.exists():
        raise FileNotFoundError('frozen v2 record required as historical context')
    v2_sha = _digest(V2_RECORD)
    if v2_sha != FROZEN_V2_SHA256:
        raise ValueError('frozen v2 JSON bytes changed')
    started = time.process_time()
    identities = coefficient_identities()
    current = e_minus2_current_bookkeeping()
    sigma = e_minus2_sigma_contraction()
    report = finite_history_uv_coefficient_report()
    pairs = original_angular_pair_ledger()
    mu = pairs[0]['multiplicity_per_signed_family']
    pilot = current_history_pilot()
    h = np.asarray(pilot['h_z'][0])
    equal = paired_leading_e_minus2_weight(h, -h, mu, mu)
    broken = paired_leading_e_minus2_weight(h, -h, mu, 0.5 * mu)
    cpu = time.process_time() - started
    return {
        'schema': 'NSC-KS-CURRENT-UV-COEFFICIENTS-v3',
        'accountable_author': 'Douglas Ek',
        'status': (
            'OPEN: leading paired e^{-2} vacuum coefficient cancels for original '
            'equal-mu same-cutoff angular pairs; C_M, C4, thermal remainder and '
            'local gate unbounded'
        ),
        'historical_context': {
            'commit': _historical_commit(),
            'predecessor_record': V2_RECORD,
            'predecessor_sha256': v2_sha,
            'predecessor_schema': 'NSC-KS-CURRENT-UV-COEFFICIENTS-v2',
        },
        'identities': {
            'V': identities['V'],
            'T_s': identities['T_s'],
            'A0': identities['A0'],
            'minor_A1': identities['minor_A1'],
            'T_s_major_A1': identities['T_s_major_A1'],
            'Im_T_s_major_A2_pure_imag_major_A1': identities[
                'Im_T_s_major_A2_pure_imag_major_A1'],
            'dummy_symbol_L0_A_j': identities['dummy_symbol_L0_A_j'],
            'jet_residuals': dict(identities['residuals']),
            'sigma_residuals': dict(sigma['residuals']),
        },
        'e_minus2_current': {
            'n3': sigma['n3'],
            'Sigma_I_momentum': sigma['Sigma_I_momentum_e_minus2'],
            'Sigma_S3_momentum': sigma['Sigma_S3_momentum_e_minus2'],
            'action_N2': sigma['action_N2'],
            'action_beta2': sigma['action_beta2'],
            'local_delta_r_rho_cancels_in_N': True,
            'A3_solved_componentwise': False,
            'raw_expansion_still_contains_n3': True,
            'linear_fixed_transfer_e_minus2_zero': current[
                'linear_fixed_transfer_e_minus2_zero'],
        },
        'angular_pairs': {
            'count': len(pairs),
            'equal_mu_and_same_cutoff': True,
            'equal_mu_weighted_h_z_max_abs': float(np.max(np.abs(equal))),
            'broken_mu_weighted_h_z_max_abs': float(np.max(np.abs(broken))),
        },
        'upstream': dict(report['upstream']),
        'vertices': dict(report['vertices']),
        'gaps': dict(named_gaps()),
        'pilot': {**pilot, 'indicator_only': True},
        'cpu_seconds_total': cpu,
        'physical_local_gate': report['physical_local_gate'],
        'certificate_from_asymptotic_orders_alone': False,
        'leading_paired_e_minus2_cancels': True,
        'C_M_or_physical_gate_from_leading_cancellation': False,
        'new_action_term': False,
        'existing_source_ad759_evaluation_changed': False,
        'source_hashes': {path: _digest(path) for path in OWNERS},
        'input_hashes': {path: _digest(path) for path in INPUTS},
        'reproducer': 'python3 scripts/derive_nsc_ks_current_uv_coefficient_control_v2.py --check',
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--check', action='store_true')
    mode.add_argument('--record', action='store_true')
    args = parser.parse_args()
    result = record()

    def frozen(payload):
        copy = json.loads(json.dumps(payload))
        copy.pop('cpu_seconds_total', None)
        copy.get('pilot', {}).pop('cpu_seconds', None)
        return copy

    if args.record:
        if OUTPUT.exists():
            raise FileExistsError('v3 UV coefficient control already recorded; use --check')
        if _digest(V2_RECORD) != FROZEN_V2_SHA256:
            raise ValueError('refusing to record v3 after v2 bytes changed')
        OUTPUT.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    elif OUTPUT.exists() and frozen(json.loads(OUTPUT.read_text())) != frozen(result):
        raise ValueError('v3 UV coefficient control replay differs')
    summary = {
        'schema': result['schema'],
        'status': result['status'],
        'predecessor_sha256': result['historical_context']['predecessor_sha256'],
        'leading_paired_e_minus2_cancels': True,
        'numerical_C_M': None,
        'physical_local_gate': 'OPEN',
        'cpu_seconds_total': result['cpu_seconds_total'],
        'angular_pairs': result['angular_pairs']['count'],
        'broken_mu_control': result['angular_pairs']['broken_mu_weighted_h_z_max_abs'],
    }
    print(json.dumps(summary, indent=2))
