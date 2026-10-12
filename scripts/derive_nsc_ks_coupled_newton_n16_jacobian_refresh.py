#!/usr/bin/env python3
"""Name the n=16 leftover stall from the already refreshed iterate4 Jacobian."""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np
from numpy.polynomial.chebyshev import chebfit
from scipy.interpolate import BarycentricInterpolator

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
sys.path.insert(0, str(ROOT/'scripts'))
import derive_nsc_ks_coupled_newton_n16_damped_iterate4 as I4
import derive_nsc_ks_coupled_newton_n16_damped_next5 as N5
from recursive_horizons.nsc_ks_profile_identity import profile_identity
from recursive_horizons.nsc_local_incoming_family import LocalIncomingFamily

OUTPUT = ROOT/'results/development/nsc-ks-coupled-newton-n16-jacobian-refresh.json'
DENSE_COUNT = 129
OWNERS = (
    'scripts/derive_nsc_ks_coupled_newton_n16_jacobian_refresh.py',
    'docs/nsc-ks-coupled-newton-n16-jacobian-refresh.md',
    'src/recursive_horizons/nsc_local_incoming_family.py',
)
STALL = 'n16_damped_frozen_jacobian_leftover_cannot_reach_existence'


def digest(path):
    path = Path(path)
    return sha256((path if path.is_absolute() else ROOT/path).read_bytes()).hexdigest()


def chebyshev_tail(z, values, center, halfwidth):
    x = (np.asarray(z, float) - center) / halfwidth
    coeff = chebfit(x, np.asarray(values, float), deg=min(len(z) - 1, 46))
    return float(np.max(abs(coeff[-4:])))


def compute():
    trial = json.loads(I4.OUTPUT.read_text())
    proposal = json.loads(N5.OUTPUT.read_text())
    for path, expected in {**trial['source_hashes'], **proposal['source_hashes'],
                           **proposal['input_hashes']}.items():
        if digest(path) != expected:
            raise ValueError('jacobian refresh input changed: '+path)
    if trial['profile_identity'] != I4.EXPECTED_IDENTITY:
        raise ValueError('refresh must use the fifth clipped history')
    if trial['new_operator_solves'] != 60 or trial['reused_operator_solves'] != 0:
        raise ValueError('refresh requires the already rebuilt sixty-family Jacobian')
    if proposal['named_search_stall'] != STALL:
        raise ValueError('next5 must already name the leftover stall')
    leftover = np.asarray(
        proposal['unclipped_linear_predicted_all_node_residual_maxima'], float)
    measured = np.asarray(trial['constraint_maxima'], float)
    if proposal['unclipped_leftover_reaches_existence_tolerance']:
        raise ValueError('leftover that reaches 3e-11 would continue Newton')
    if np.any(leftover <= N5.EXISTENCE_TOLERANCE):
        raise ValueError('both leftover components must stay above the gate')
    family = LocalIncomingFamily(np.asarray(trial['history']['coefficients'], float))
    if profile_identity(family, include_normal_window=True) != trial['profile_identity']:
        raise ValueError('refresh history identity changed')
    z = np.asarray(trial['z'], float)
    residual = np.asarray(trial['action_gradient'], float)
    dense = family.collocation_nodes(DENSE_COUNT)
    interpolants = [BarycentricInterpolator(z, residual[:, k]) for k in range(2)]
    residual_dense = np.column_stack([fn(dense) for fn in interpolants])
    tails = [round(chebyshev_tail(z, residual[:, k], family.center, family.axial_inner), 12)
             for k in range(2)]
    interpolant_maxima = np.round(np.max(abs(residual_dense), axis=0), 12)
    rect = [row for row in proposal['proposals'] if row['layout'] == 'rectangular'][0]
    n32_motivated = bool(np.any(np.asarray(tails) > 0.1 * leftover))
    return {
        'schema': 'NSC-KS-COUPLED-NEWTON-N16-JACOBIAN-REFRESH-v1',
        'accountable_author': 'Douglas Ek',
        'status': (
            'OPEN: Jacobian at the fifth clipped g is already rebuilt; '
            'unclipped leftover cannot reach EXISTENCE'),
        'profile_identity': trial['profile_identity'],
        'jacobian_reused_from': str(I4.OUTPUT.relative_to(ROOT)),
        'new_operator_solves_in_this_owner': 0,
        'already_rebuilt_operator_solves': trial['new_operator_solves'],
        'measured_constraint_maxima': measured.tolist(),
        'unclipped_linear_predicted_all_node_residual_maxima': leftover.tolist(),
        'existence_tolerance': N5.EXISTENCE_TOLERANCE,
        'unclipped_leftover_reaches_existence_tolerance': False,
        'rectangular_condition_number': rect['condition_number'],
        'unclipped_step_norm': rect['unclipped_step_norm'],
        'residual_chebyshev_tail_max_last4': tails,
        'interpolated_residual_maxima_on_dense_nodes': interpolant_maxima.tolist(),
        'n32_basis_motivated_by_residual_tail': n32_motivated,
        'named_search_stall': STALL,
        'scope': {
            'physical_local_gate': 'OPEN',
            'jacobian_already_rebuilt_at_accepted_g': True,
            'no_repeated_sixty_family_evolution': True,
            'prediction_is_not_a_nonlinear_residual': True,
            'failed_optimization_is_not_NONEXISTENCE': True,
            'physical_EXISTENCE_certificate': False,
            'physical_NONEXISTENCE_claimed': False,
            'n32_not_opened': not n32_motivated,
            'new_field_or_source_evolutions': 0,
            'metric_timestep': False,
            'forbidden_identities_not_evolved': True,
        },
        'source_hashes': {p: digest(p) for p in OWNERS},
        'input_hashes': {
            str(I4.OUTPUT.relative_to(ROOT)): digest(I4.OUTPUT),
            str(N5.OUTPUT.relative_to(ROOT)): digest(N5.OUTPUT),
        },
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--record', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = compute()
    if args.record:
        if OUTPUT.exists():
            raise FileExistsError('n=16 jacobian refresh register exists; use --check')
        OUTPUT.write_text(json.dumps(result, sort_keys=True, indent=2, allow_nan=False)+'\n')
    elif result != json.loads(OUTPUT.read_text()):
        raise ValueError('n=16 jacobian refresh replay differs')
    print(json.dumps({
        'status': result['status'],
        'leftover': result['unclipped_linear_predicted_all_node_residual_maxima'],
        'condition': result['rectangular_condition_number'],
        'chebyshev_tail': result['residual_chebyshev_tail_max_last4'],
        'n32_motivated': result['n32_basis_motivated_by_residual_tail'],
        'named_search_stall': result['named_search_stall'],
        'physical_NONEXISTENCE_claimed': result['scope']['physical_NONEXISTENCE_claimed'],
    }, indent=2))
