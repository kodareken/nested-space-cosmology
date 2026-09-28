#!/usr/bin/env python3
"""Fourier leftover-null current on the walk. Not leftover v, not L(E)."""
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
import derive_nsc_ks_fourier_wu_leftover as F
import derive_nsc_ks_n16_geometry_killing_current as K
import derive_nsc_ks_n16_leftover_anatomy as A
import derive_nsc_ks_n32_geometry_leftover as N32
from recursive_horizons.nsc_evolved_incoming_constraints import (
    compatible_history_slots, surface_geometry_response)
from recursive_horizons.nsc_local_incoming_family import LocalIncomingFamily
from recursive_horizons.nsc_local_incoming_fourier_family import LocalFourierIncomingFamily

OUTPUT = ROOT/'results/development/nsc-ks-fourier-leftover-null.json'
OWNERS = (
    'scripts/derive_nsc_ks_fourier_leftover_null.py',
    'docs/nsc-ks-fourier-leftover-null.md',
)
GEOMETRY_REMAINDER = 1e-13
NEIGHBORHOOD = (-1e-3, 0.0, 1e-3)


def digest(path):
    path = Path(path)
    return sha256((path if path.is_absolute() else ROOT/path).read_bytes()).hexdigest()


def _round(value, digits=12):
    return N32._round(value, digits)


def flatten_field(field):
    field = np.asarray(field, float)
    return np.concatenate([field[:, 0], field[:, 1]])


def leftover_null(jacobian, residual):
    residual = np.asarray(residual, float)
    jacobian = np.asarray(jacobian, float)
    delta, *_ = np.linalg.lstsq(jacobian, residual, rcond=None)
    leftover = residual - jacobian @ delta
    norm = float(np.linalg.norm(leftover))
    if not np.isfinite(norm) or norm <= 0:
        raise ValueError('Fourier leftover-null requires a nonzero leftover')
    return leftover / norm, leftover, delta


def current_of(null, field):
    return float(np.dot(null, flatten_field(field)))


def compute():
    leftover = json.loads(F.OUTPUT.read_text())
    trial = json.loads(I4.OUTPUT.read_text())
    killing = json.loads(K.OUTPUT.read_text())
    for register in (leftover, trial, killing):
        for path, expected in {**register.get('source_hashes', {}),
                               **register.get('input_hashes', {})}.items():
            if digest(path) != expected:
                raise ValueError('Fourier leftover-null input changed: '+path)
    if leftover['evolve_authorized']:
        raise ValueError('leftover-null is only for a refused Fourier evolve')
    ctx = I4.context()
    family = ctx['family']
    coeff = ctx['coeff']
    z = np.asarray(trial['z'], float)
    gradient = np.asarray(trial['action_gradient'], float)
    matter = np.asarray(trial['summed_matter_change'], float)
    edge = np.asarray(trial['source_cutoff_edge_gradient'], float)
    raw = np.asarray(trial['raw_action_gradient'], float)
    slots, _ = compatible_history_slots(family.metric(), z, 32)
    count = int(leftover['selected_mode_count'])
    fourier = LocalFourierIncomingFamily(np.zeros((2, count)))
    _, fourier_ds = compatible_history_slots(fourier.metric(), z, 2 * count)
    if leftover['orthonormalized_columns']:
        fourier_ds = F.orthonormalize_tangents(fourier_ds)
    geometry = surface_geometry_response(slots, fourier_ds, coeff)
    residual, jacobian, _n_nodes = N32.flatten_system(
        gradient, geometry['action_gradient_tangent'])
    null, leftover_vec, _delta = leftover_null(jacobian, residual)
    geometry_change = surface_geometry_response(
        slots, np.zeros((0, len(z), 6)), coeff)['action_gradient_change']
    baseline = raw - geometry_change - matter
    pieces = {
        'baseline': current_of(null, baseline),
        'geometry': current_of(null, geometry_change),
        'matter': current_of(null, matter),
        'edge': current_of(null, edge),
        'assembled': current_of(null, gradient),
    }
    walk = []
    for path in A.WALK:
        record = json.loads(path.read_text())
        hist = LocalIncomingFamily(np.asarray(record['history']['coefficients'], float))
        hist_z = np.asarray(record['z'], float)
        hist_residual = np.asarray(record['action_gradient'], float)
        if hist_z.shape != z.shape or not np.allclose(hist_z, z):
            raise ValueError('walk residual must share the iterate4 node set')
        walk.append({
            'path': str(path.relative_to(ROOT)),
            'profile_identity': record['profile_identity'],
            'current': _round(current_of(null, hist_residual)),
        })
    walk_values = np.array([row['current'] for row in walk], float)
    walk_sign_stable = bool(np.all(walk_values > 0) or np.all(walk_values < 0))
    neighborhood = []
    for direction, tangent in enumerate(fourier_ds):
        for scale in NEIGHBORHOOD:
            perturbed = slots + scale * tangent
            geometry_shift = surface_geometry_response(
                perturbed, fourier_ds, coeff)['action_gradient_change']
            field = gradient - geometry_change + geometry_shift
            neighborhood.append({
                'direction': int(direction),
                'scale': scale,
                'current': _round(current_of(null, field)),
            })
    neighborhood_values = np.array([row['current'] for row in neighborhood], float)
    neighborhood_sign_stable = bool(
        np.all(neighborhood_values > 0) or np.all(neighborhood_values < 0))
    gap = abs(pieces['assembled'])
    edge_can_eat = bool(abs(pieces['edge']) >= gap)
    missing_tail = True
    geometry_killed = bool(abs(pieces['geometry']) <= 1000 * GEOMETRY_REMAINDER)
    certificate = bool(
        walk_sign_stable
        and neighborhood_sign_stable
        and geometry_killed
        and gap > GEOMETRY_REMAINDER
        and not edge_can_eat
        and not missing_tail)
    return {
        'schema': 'NSC-KS-FOURIER-LEFTOVER-NULL-v1',
        'accountable_author': 'Douglas Ek',
        'status': (
            'OPEN: Fourier leftover-null is not a class identity'
            if not certificate else
            'scoped NON-EXISTENCE: Fourier leftover-null has an error-controlled gap'),
        'profile_identity': trial['profile_identity'],
        'not_n16_walk_v': True,
        'L_E_not_reevaluated': True,
        'selected_mode_count': count,
        'leftover_null_norm': 1.0,
        'leftover_euclidean_norm': _round(float(np.linalg.norm(leftover_vec))),
        'pieces': {name: _round(value) for name, value in pieces.items()},
        'walk': walk,
        'walk_sign_stable': walk_sign_stable,
        'walk_min': _round(np.min(walk_values)),
        'walk_max': _round(np.max(walk_values)),
        'neighborhood_scales': list(NEIGHBORHOOD),
        'neighborhood_sign_stable': neighborhood_sign_stable,
        'neighborhood_min': _round(np.min(neighborhood_values)),
        'neighborhood_max': _round(np.max(neighborhood_values)),
        'geometry_killed': geometry_killed,
        'edge_can_eat_gap': edge_can_eat,
        'unenclosed_edge_tail': None,
        'locked_current_identified': certificate,
        'physical_NONEXISTENCE_certificate': certificate,
        'existence_tolerance': A.EXISTENCE_TOLERANCE,
        'scope': {
            'physical_local_gate': 'OPEN' if not certificate else 'NON-EXISTENCE',
            'finite_span_is_not_the_class': True,
            'failed_optimization_is_not_NONEXISTENCE': True,
            'missing_tail_blocks_gap': missing_tail,
            'new_field_or_source_evolutions': 0,
            'physical_EXISTENCE_certificate': False,
            'metric_timestep': False,
        },
        'source_hashes': {p: digest(p) for p in OWNERS},
        'input_hashes': {
            str(F.OUTPUT.relative_to(ROOT)): digest(F.OUTPUT),
            str(I4.OUTPUT.relative_to(ROOT)): digest(I4.OUTPUT),
            str(K.OUTPUT.relative_to(ROOT)): digest(K.OUTPUT),
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
            raise FileExistsError('Fourier leftover-null exists; use --check')
        OUTPUT.write_text(json.dumps(result, sort_keys=True, indent=2, allow_nan=False)+'\n')
    elif result != json.loads(OUTPUT.read_text()):
        raise ValueError('Fourier leftover-null replay differs')
    print(json.dumps({
        'status': result['status'],
        'physical_NONEXISTENCE_certificate': result['physical_NONEXISTENCE_certificate'],
        'walk_sign_stable': result['walk_sign_stable'],
        'neighborhood_sign_stable': result['neighborhood_sign_stable'],
        'walk_min': result['walk_min'],
        'walk_max': result['walk_max'],
        'edge_can_eat_gap': result['edge_can_eat_gap'],
    }, indent=2))
