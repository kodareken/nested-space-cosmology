#!/usr/bin/env python3
"""Score the unused flux first integral on iterate4/walk. No L(E), no leftover v."""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np
from numpy.polynomial.polynomial import polyder, polyval

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
sys.path.insert(0, str(ROOT/'scripts'))
import derive_nsc_ks_coupled_newton_n16_damped_iterate4 as I4
import derive_nsc_ks_fourier_wu_evolve_refusal as R
import derive_nsc_ks_n16_geometry_killing_current as K
import derive_nsc_ks_n16_leftover_anatomy as A
from recursive_horizons.nsc_evolved_incoming_constraints import (
    compatible_history_slots, surface_geometry_response)
from recursive_horizons.nsc_local_incoming_family import LocalIncomingFamily
from recursive_horizons.nsc_local_incoming_fourier_family import LocalFourierIncomingFamily

OUTPUT = ROOT/'results/development/nsc-ks-flux-first-integral.json'
OWNERS = (
    'scripts/derive_nsc_ks_flux_first_integral.py',
    'docs/nsc-ks-flux-first-integral.md',
)
GEOMETRY_REMAINDER = 1e-13
NEIGHBORHOOD = (-1e-3, 0.0, 1e-3)


def digest(path):
    path = Path(path)
    return sha256((path if path.is_absolute() else ROOT/path).read_bytes()).hexdigest()


def _round(value, digits=12):
    return K._round(value, digits)


def summary(values):
    return K.summary(values)


def first_integral(slots, coeff):
    w, _wz, wzz, _wzzz, U, _Uz = np.asarray(slots, float).T
    return coeff['d'] * U + coeff['e'] * wzz + polyval(w, coeff['F'])


def first_integral_derivative(slots, coeff):
    w, wz, _wzz, wzzz, _U, Uz = np.asarray(slots, float).T
    return coeff['d'] * Uz + coeff['e'] * wzzz + polyval(w, polyder(coeff['F'])) * wz


def flux_divergence(slots, coeff):
    w, wz, wzz, _wzzz, U, _Uz = np.asarray(slots, float).T
    a, r, Hr = coeff['a'], coeff['r'], coeff['Hr']
    A0, A1, d, e, c = coeff['A0'], coeff['A1'], coeff['d'], coeff['e'], coeff['C']
    A = A0 + A1 * w
    B = -(2 * A + (Hr + w / r) * d) / a ** 2
    delta = A * e - B * d
    F = polyval(w, coeff['F'])
    Dval = polyval(w, coeff['D'])
    Pi = d * U + e * wzz + F
    return -d * c * wz * wz + delta * wzz - A * (Pi - F) - d * Dval


def balance(slots, residual, z, coeff):
    dPi = first_integral_derivative(slots, coeff)
    beta = np.asarray(residual, float)[:, 1]
    Pi = first_integral(slots, coeff)
    integrated = float(Pi[-1] - Pi[0] + np.trapezoid(beta, np.asarray(z, float)))
    return dPi + beta, integrated, flux_divergence(slots, coeff)


def compute():
    killing = json.loads(K.OUTPUT.read_text())
    refusal = json.loads(R.OUTPUT.read_text())
    trial = json.loads(I4.OUTPUT.read_text())
    for register in (killing, refusal, trial):
        for path, expected in {**register.get('source_hashes', {}),
                               **register.get('input_hashes', {})}.items():
            if digest(path) != expected:
                raise ValueError('flux first-integral input changed: '+path)
    if refusal['evolve_authorized'] or killing['physical_NONEXISTENCE_certificate']:
        raise ValueError('flux owner is only for a failed leftover / failed L(E)')
    ctx = I4.context()
    family = ctx['family']
    coeff = ctx['coeff']
    z = np.asarray(trial['z'], float)
    gradient = np.asarray(trial['action_gradient'], float)
    matter = np.asarray(trial['summed_matter_change'], float)
    edge = np.asarray(trial['source_cutoff_edge_gradient'], float)
    raw = np.asarray(trial['raw_action_gradient'], float)
    slots, _ = compatible_history_slots(family.metric(), z, 32)
    geometry_change = surface_geometry_response(
        slots, np.zeros((0, len(z), 6)), coeff)['action_gradient_change']
    baseline = raw - geometry_change - matter
    pieces = {
        'baseline': baseline,
        'geometry': geometry_change,
        'matter': matter,
        'edge': edge,
        'assembled': gradient,
    }
    integrals = {}
    divergences = {}
    for name, field in pieces.items():
        pointwise, integrated, divergence = balance(slots, field, z, coeff)
        integrals[name] = {**summary(pointwise), 'integrated': _round(integrated)}
        divergences[name] = summary(divergence)
    walk = []
    for path in A.WALK:
        record = json.loads(path.read_text())
        hist = LocalIncomingFamily(np.asarray(record['history']['coefficients'], float))
        hist_z = np.asarray(record['z'], float)
        residual = np.asarray(record['action_gradient'], float)
        hist_slots, _ = compatible_history_slots(hist.metric(), hist_z, 32)
        pointwise, integrated, divergence = balance(hist_slots, residual, hist_z, coeff)
        walk.append({
            'path': str(path.relative_to(ROOT)),
            'profile_identity': record['profile_identity'],
            'integral': {**summary(pointwise), 'integrated': _round(integrated)},
            'divergence': summary(divergence),
        })
    walk_int_mins = np.array([row['integral']['min'] for row in walk], float)
    walk_int_maxs = np.array([row['integral']['max'] for row in walk], float)
    walk_div_mins = np.array([row['divergence']['min'] for row in walk], float)
    walk_div_maxs = np.array([row['divergence']['max'] for row in walk], float)
    walk_integral_sign_stable = bool(
        (np.all(walk_int_mins > 0) and np.all(walk_int_maxs > 0))
        or (np.all(walk_int_mins < 0) and np.all(walk_int_maxs < 0)))
    walk_divergence_sign_stable = bool(
        (np.all(walk_div_mins > 0) and np.all(walk_div_maxs > 0))
        or (np.all(walk_div_mins < 0) and np.all(walk_div_maxs < 0)))
    fourier = LocalFourierIncomingFamily(np.zeros((2, 8)))
    _, fourier_ds = compatible_history_slots(fourier.metric(), z, 16)
    neighborhood = []
    for direction, tangent in enumerate(fourier_ds):
        for scale in NEIGHBORHOOD:
            perturbed = slots + scale * tangent
            geometry_shift = surface_geometry_response(
                perturbed, fourier_ds, coeff)['action_gradient_change']
            residual = gradient - geometry_change + geometry_shift
            pointwise, integrated, divergence = balance(perturbed, residual, z, coeff)
            neighborhood.append({
                'direction': direction,
                'scale': scale,
                'integral': summary(pointwise) | {'integrated': _round(integrated)},
                'divergence': summary(divergence),
            })
    neighborhood_integral_stable = bool(
        all(row['integral']['sign_stable'] for row in neighborhood)
        and (
            all(row['integral']['min'] > 0 for row in neighborhood)
            or all(row['integral']['max'] < 0 for row in neighborhood)))
    neighborhood_divergence_stable = bool(
        all(row['divergence']['sign_stable'] for row in neighborhood)
        and (
            all(row['divergence']['min'] > 0 for row in neighborhood)
            or all(row['divergence']['max'] < 0 for row in neighborhood)))
    assembled = integrals['assembled']
    gap = assembled['abs_min']
    edge_can_eat = bool(integrals['edge']['abs_max'] >= gap)
    missing_tail = True
    geometry_killed = bool(integrals['geometry']['abs_max'] <= 1000 * GEOMETRY_REMAINDER)
    divergence_geometry_killed = bool(
        divergences['geometry']['abs_max'] <= 1000 * GEOMETRY_REMAINDER)
    integral_certificate = bool(
        walk_integral_sign_stable
        and assembled['sign_stable']
        and neighborhood_integral_stable
        and geometry_killed
        and gap is not None
        and gap > GEOMETRY_REMAINDER
        and not edge_can_eat
        and not missing_tail)
    divergence_certificate = bool(
        walk_divergence_sign_stable
        and divergences['assembled']['sign_stable']
        and neighborhood_divergence_stable
        and divergence_geometry_killed
        and divergences['assembled']['abs_min'] > GEOMETRY_REMAINDER
        and divergences['edge']['abs_max'] < divergences['assembled']['abs_min']
        and not missing_tail)
    certificate = bool(integral_certificate or divergence_certificate)
    return {
        'schema': 'NSC-KS-FLUX-FIRST-INTEGRAL-v1',
        'accountable_author': 'Douglas Ek',
        'status': (
            'OPEN: flux first integral and A-free divergence are not locked currents'
            if not certificate else
            'scoped NON-EXISTENCE: flux current has an error-controlled gap'),
        'profile_identity': trial['profile_identity'],
        'leftover_v_not_reused': True,
        'L_E_not_reevaluated': True,
        'pieces': integrals,
        'divergences': divergences,
        'walk': walk,
        'walk_integral_sign_stable': walk_integral_sign_stable,
        'walk_divergence_sign_stable': walk_divergence_sign_stable,
        'walk_integral_min': _round(np.min(walk_int_mins)),
        'walk_integral_max': _round(np.max(walk_int_maxs)),
        'walk_divergence_min': _round(np.min(walk_div_mins)),
        'walk_divergence_max': _round(np.max(walk_div_maxs)),
        'neighborhood_scales': list(NEIGHBORHOOD),
        'neighborhood_integral_sign_stable': neighborhood_integral_stable,
        'neighborhood_divergence_sign_stable': neighborhood_divergence_stable,
        'geometry_killed_by_integral': geometry_killed,
        'geometry_killed_by_divergence': divergence_geometry_killed,
        'edge_can_eat_gap': edge_can_eat,
        'unenclosed_edge_tail': None,
        'locked_current_identified': certificate,
        'physical_NONEXISTENCE_certificate': certificate,
        'existence_tolerance': A.EXISTENCE_TOLERANCE,
        'scope': {
            'physical_local_gate': 'OPEN' if not certificate else 'NON-EXISTENCE',
            'failed_optimization_is_not_NONEXISTENCE': True,
            'missing_tail_blocks_gap': missing_tail,
            'new_field_or_source_evolutions': 0,
            'physical_EXISTENCE_certificate': False,
            'metric_timestep': False,
            'principal_det_not_promoted': True,
        },
        'source_hashes': {p: digest(p) for p in OWNERS},
        'input_hashes': {
            str(K.OUTPUT.relative_to(ROOT)): digest(K.OUTPUT),
            str(R.OUTPUT.relative_to(ROOT)): digest(R.OUTPUT),
            str(I4.OUTPUT.relative_to(ROOT)): digest(I4.OUTPUT),
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
            raise FileExistsError('flux first integral exists; use --check')
        OUTPUT.write_text(json.dumps(result, sort_keys=True, indent=2, allow_nan=False)+'\n')
    elif result != json.loads(OUTPUT.read_text()):
        raise ValueError('flux first integral replay differs')
    print(json.dumps({
        'status': result['status'],
        'physical_NONEXISTENCE_certificate': result['physical_NONEXISTENCE_certificate'],
        'walk_integral_sign_stable': result['walk_integral_sign_stable'],
        'walk_divergence_sign_stable': result['walk_divergence_sign_stable'],
        'neighborhood_integral_sign_stable': result['neighborhood_integral_sign_stable'],
        'edge_can_eat_gap': result['edge_can_eat_gap'],
        'walk_integral_min': result['walk_integral_min'],
        'walk_integral_max': result['walk_integral_max'],
    }, indent=2))
