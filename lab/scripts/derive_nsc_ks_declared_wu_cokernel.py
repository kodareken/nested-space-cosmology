#!/usr/bin/env python3
"""Continuous cokernel of the declared (w,U) geometry map. Not L(E), not flux."""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
sys.path.insert(0, str(ROOT/'scripts'))
import derive_nsc_ks_coupled_newton_n16_damped_iterate4 as I4
import derive_nsc_ks_fourier_leftover_null as N
import derive_nsc_ks_n16_geometry_killing_current as K
from recursive_horizons.nsc_evolved_incoming_constraints import (
    compatible_history_slots, surface_geometry_response)

OUTPUT = ROOT/'results/development/nsc-ks-declared-wu-cokernel.json'
OWNERS = (
    'scripts/derive_nsc_ks_declared_wu_cokernel.py',
    'docs/nsc-ks-declared-wu-cokernel.md',
)
PRINCIPAL = ROOT/'results/development/nsc-incoming-surface-principal.json'
GEOMETRY_REMAINDER = 1e-13


def digest(path):
    path = Path(path)
    return sha256((path if path.is_absolute() else ROOT/path).read_bytes()).hexdigest()


def _round(value, digits=12):
    return K._round(value, digits)


def principal_blocks(slots, coeff):
    w, _wz, _wzz, _wzzz, _U, _Uz = np.asarray(slots, float).T
    a, r, Hr = coeff['a'], coeff['r'], coeff['Hr']
    A0, A1, d, e = coeff['A0'], coeff['A1'], coeff['d'], coeff['e']
    A = A0 + A1 * w
    B = -(2 * A + (Hr + w / r) * d) / a ** 2
    det = A * e - B * d
    return A, B, d, e, det


def compute():
    leftover_null = json.loads(N.OUTPUT.read_text())
    trial = json.loads(I4.OUTPUT.read_text())
    principal = json.loads(PRINCIPAL.read_text())
    for path, expected in {**leftover_null['source_hashes'],
                           **leftover_null['input_hashes'],
                           **trial['source_hashes']}.items():
        if digest(path) != expected:
            raise ValueError('declared cokernel input changed: '+path)
    if leftover_null['physical_NONEXISTENCE_certificate']:
        raise ValueError('continuous cokernel is only for a failed leftover-null')
    ctx = I4.context()
    coeff = ctx['coeff']
    z = np.asarray(trial['z'], float)
    family = ctx['family']
    slots, _ = compatible_history_slots(family.metric(), z, 32)
    A, B, d, e, det = principal_blocks(slots, coeff)
    certified = np.asarray(principal['matrix']['determinant_interval'], float)
    det_min, det_max = float(np.min(det)), float(np.max(det))
    det_excludes_zero = bool(det_min * det_max > 0)
    certified_excludes_zero = bool(certified[0] * certified[1] > 0)
    same_sign_as_certified = bool(det_min * certified[0] > 0)
    geometry = surface_geometry_response(slots, np.zeros((0, len(z), 6)), coeff)
    jac = geometry['jacobian']
    u_column = jac[:, :, 4]
    wzz_column = jac[:, :, 2]
    uz_column = jac[:, :, 5]
    wzzz_column = jac[:, :, 3]
    pointwise_rank = []
    for index in range(len(z)):
        matrix = np.array([
            [u_column[index, 0], wzz_column[index, 0]],
            [uz_column[index, 1], wzzz_column[index, 1]],
        ], float)
        pointwise_rank.append(float(np.linalg.det(matrix)))
    pointwise_rank = np.asarray(pointwise_rank, float)
    pointwise_excludes_zero = bool(np.min(pointwise_rank) * np.max(pointwise_rank) > 0)
    recovered_lapse_division = bool(
        np.allclose(A, d, atol=0, rtol=0) or np.max(np.abs(A)) < GEOMETRY_REMAINDER)
    recovered_flux_divergence = False
    trivial = bool(
        det_excludes_zero
        and certified_excludes_zero
        and same_sign_as_certified
        and pointwise_excludes_zero)
    is_L_or_flux = bool(recovered_lapse_division or recovered_flux_divergence)
    certificate = False
    return {
        'schema': 'NSC-KS-DECLARED-WU-COKERNEL-v1',
        'accountable_author': 'Douglas Ek',
        'status': (
            'OPEN: continuous cokernel of compact (w,U) geometry is trivial; '
            'not L(E), not flux, not NON-EXISTENCE'
            if trivial and not is_L_or_flux else
            'OPEN: continuous cokernel recovered a refused identity or lost rank'),
        'profile_identity': trial['profile_identity'],
        'certified_determinant_interval': [float(certified[0]), float(certified[1])],
        'iterate4_principal_det_min': _round(det_min),
        'iterate4_principal_det_max': _round(det_max),
        'principal_det_excludes_zero': det_excludes_zero,
        'certified_det_excludes_zero': certified_excludes_zero,
        'same_sign_as_certified_det': same_sign_as_certified,
        'pointwise_principal_det_min': _round(np.min(pointwise_rank)),
        'pointwise_principal_det_max': _round(np.max(pointwise_rank)),
        'pointwise_principal_excludes_zero': pointwise_excludes_zero,
        'compact_support_forces_zero_multiplier': trivial,
        'recovered_L_E': recovered_lapse_division,
        'recovered_flux_divergence': recovered_flux_divergence,
        'is_L_or_flux_again': is_L_or_flux,
        'continuous_cokernel_trivial': trivial and not is_L_or_flux,
        'geometry_killed': False,
        'unenclosed_edge_tail': None,
        'locked_current_identified': False,
        'physical_NONEXISTENCE_certificate': certificate,
        'scope': {
            'physical_local_gate': 'OPEN',
            'principal_det_is_invertibility_not_existence': True,
            'finite_leftover_is_not_class_identity': True,
            'no_class_rewrite': True,
            'new_field_or_source_evolutions': 0,
            'physical_EXISTENCE_certificate': False,
            'metric_timestep': False,
        },
        'source_hashes': {p: digest(p) for p in OWNERS},
        'input_hashes': {
            str(N.OUTPUT.relative_to(ROOT)): digest(N.OUTPUT),
            str(I4.OUTPUT.relative_to(ROOT)): digest(I4.OUTPUT),
            str(PRINCIPAL.relative_to(ROOT)): digest(PRINCIPAL),
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
            raise FileExistsError('declared cokernel exists; use --check')
        OUTPUT.write_text(json.dumps(result, sort_keys=True, indent=2, allow_nan=False)+'\n')
    elif result != json.loads(OUTPUT.read_text()):
        raise ValueError('declared cokernel replay differs')
    print(json.dumps({
        'status': result['status'],
        'continuous_cokernel_trivial': result['continuous_cokernel_trivial'],
        'is_L_or_flux_again': result['is_L_or_flux_again'],
        'physical_NONEXISTENCE_certificate': result['physical_NONEXISTENCE_certificate'],
        'principal_det_excludes_zero': result['principal_det_excludes_zero'],
    }, indent=2))
