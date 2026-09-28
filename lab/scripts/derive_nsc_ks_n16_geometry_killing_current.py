#!/usr/bin/env python3
"""Evaluate L(E)=E_beta-(d/A(w)) d_z E_N on saved n=16 pieces. No Dirac."""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np
from numpy.polynomial.chebyshev import chebder, chebfit, chebval

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
sys.path.insert(0, str(ROOT/'scripts'))
import derive_nsc_ks_coupled_newton_n16_damped_iterate4 as I4
import derive_nsc_ks_n16_leftover_anatomy as A
import derive_nsc_ks_n16_leftover_identification as ID
from recursive_horizons.nsc_evolved_incoming_constraints import (
    compatible_history_slots, surface_geometry_response)
from recursive_horizons.nsc_local_incoming_family import LocalIncomingFamily

OUTPUT = ROOT/'results/development/nsc-ks-n16-geometry-killing-current.json'
OWNERS = (
    'scripts/derive_nsc_ks_n16_geometry_killing_current.py',
    'docs/nsc-ks-n16-geometry-killing-current.md',
)
BRANCH = ROOT/'results/development/nsc-incoming-surface-regular-branch.json'
EXISTENCE = A.EXISTENCE_TOLERANCE
GEOMETRY_REMAINDER = 1e-13


def digest(path):
    path = Path(path)
    return sha256((path if path.is_absolute() else ROOT/path).read_bytes()).hexdigest()


def _round(value, digits=12):
    array = np.asarray(value, float)
    if array.shape == ():
        return None if not np.isfinite(array) else round(float(array), digits)
    return [None if not np.isfinite(item) else round(float(item), digits)
            for item in array.ravel()]


def axial_derivative(z, values, center, halfwidth):
    x = (np.asarray(z, float) - center) / halfwidth
    values = np.asarray(values, float)
    coeff = chebfit(x, values, deg=min(len(z) - 1, 46))
    return chebval(x, chebder(coeff)) / halfwidth


def A_of_w(w, A0, A1):
    return A0 + A1 * np.asarray(w, float)


def killing(E, A, d, dN):
    if np.any(A == 0) or not np.isfinite(A).all():
        raise ValueError('A(w) must stay finite and nonzero on I')
    return np.asarray(E[:, 1], float) - (d / A) * dN


def interval_A(w, A0, A1):
    terms = np.vstack((A1[0] * w, A1[1] * w))
    return A0[0] + np.min(terms, axis=0), A0[1] + np.max(terms, axis=0)


def interval_ratio(d_int, A_lo, A_hi):
    if np.any(A_lo <= 0) and np.any(A_hi >= 0):
        if np.any((A_lo <= 0) & (A_hi >= 0)):
            raise ValueError('A(w) interval contains zero; elimination pivot fails')
    corners = np.vstack((
        d_int[0] / A_lo, d_int[0] / A_hi, d_int[1] / A_lo, d_int[1] / A_hi))
    return np.min(corners, axis=0), np.max(corners, axis=0)


def summary(values):
    values = np.asarray(values, float)
    return {
        'min': _round(np.min(values)),
        'max': _round(np.max(values)),
        'abs_min': _round(np.min(np.abs(values))),
        'abs_max': _round(np.max(np.abs(values))),
        'sign_stable': bool(np.all(values > 0) or np.all(values < 0)),
    }


def compute():
    identification = json.loads(ID.OUTPUT.read_text())
    trial = json.loads(I4.OUTPUT.read_text())
    branch = json.loads(BRANCH.read_text())
    for path, expected in {**identification['source_hashes'],
                           **identification['input_hashes'],
                           **trial['source_hashes']}.items():
        if digest(path) != expected:
            raise ValueError('geometry-killing input changed: '+path)
    ctx = I4.context()
    family = ctx['family']
    z = np.asarray(trial['z'], float)
    gradient = np.asarray(trial['action_gradient'], float)
    matter = np.asarray(trial['summed_matter_change'], float)
    edge = np.asarray(trial['source_cutoff_edge_gradient'], float)
    raw = np.asarray(trial['raw_action_gradient'], float)
    slots, slot_tangents = compatible_history_slots(family.metric(), z, 32)
    geometry = surface_geometry_response(
        slots, slot_tangents, ctx['coeff'])['action_gradient_change']
    baseline = raw - geometry - matter
    imported = branch['branch']['imported_unchanged']
    A0_int = np.asarray(imported['A0'], float)
    d_int = np.asarray(imported['d'], float)
    A1_int = np.asarray(branch['branch']['intervals']['A1_total'], float)
    A0 = float(ctx['coeff']['A0'])
    A1 = float(ctx['coeff']['A1'])
    d = float(ctx['coeff']['d'])
    w = slots[:, 0]
    A_values = A_of_w(w, A0, A1)
    A_lo, A_hi = interval_A(w, A0_int, A1_int)
    ratio_lo, ratio_hi = interval_ratio(d_int, A_lo, A_hi)
    pieces = {
        'baseline': baseline,
        'geometry': geometry,
        'matter': matter,
        'edge': edge,
        'assembled': gradient,
    }
    L = {}
    for name, field in pieces.items():
        dN = axial_derivative(z, field[:, 0], family.center, family.axial_inner)
        values = killing(field, A_values, d, dN)
        L[name] = {
            **summary(values),
            'axial_range': _round([np.min(field[:, 0]), np.max(field[:, 0])]),
            'beta_range': _round([np.min(field[:, 1]), np.max(field[:, 1])]),
        }
        if name == 'assembled':
            L[name]['interval_min'] = _round(np.min(
                field[:, 1] - np.where(dN >= 0, ratio_hi, ratio_lo) * dN))
            L[name]['interval_max'] = _round(np.max(
                field[:, 1] - np.where(dN >= 0, ratio_lo, ratio_hi) * dN))
    walk = []
    for path in A.WALK:
        record = json.loads(path.read_text())
        hist = LocalIncomingFamily(np.asarray(record['history']['coefficients'], float))
        residual = np.asarray(record['action_gradient'], float)
        hist_z = np.asarray(record['z'], float)
        hist_slots, hist_ds = compatible_history_slots(hist.metric(), hist_z, 32)
        hist_A = A_of_w(hist_slots[:, 0], A0, A1)
        dN = axial_derivative(
            hist_z, residual[:, 0], hist.center, hist.axial_inner)
        values = killing(residual, hist_A, d, dN)
        walk.append({
            'path': str(path.relative_to(ROOT)),
            'profile_identity': record['profile_identity'],
            **summary(values),
        })
    walk_mins = np.array([row['min'] for row in walk], float)
    walk_maxs = np.array([row['max'] for row in walk], float)
    walk_sign_stable = bool(
        (np.all(walk_mins > 0) and np.all(walk_maxs > 0))
        or (np.all(walk_mins < 0) and np.all(walk_maxs < 0)))
    baseline_axial_span = float(np.ptp(baseline[:, 0]))
    baseline_beta_span = float(np.ptp(baseline[:, 1]))
    baseline_constant = bool(
        baseline_axial_span <= 100 * GEOMETRY_REMAINDER
        and baseline_beta_span <= 100 * GEOMETRY_REMAINDER)
    geometry_killed = bool(L['geometry']['abs_max'] <= 1000 * GEOMETRY_REMAINDER)
    assembled = L['assembled']
    gap = assembled['abs_min']
    edge_can_eat = bool(L['edge']['abs_max'] >= gap)
    missing_tail = True
    certificate = bool(
        walk_sign_stable
        and assembled['sign_stable']
        and geometry_killed
        and gap is not None
        and gap > GEOMETRY_REMAINDER
        and not edge_can_eat
        and not missing_tail
        and baseline_constant)
    return {
        'schema': 'NSC-KS-N16-GEOMETRY-KILLING-CURRENT-v1',
        'accountable_author': 'Douglas Ek',
        'status': (
            'OPEN: L(E) is not a locked P_K-style current; '
            'no scoped NON-EXISTENCE certificate'
            if not certificate else
            'scoped NON-EXISTENCE: L(E) has an error-controlled gap'),
        'profile_identity': trial['profile_identity'],
        'A0': _round(A0),
        'A1': _round(A1),
        'd': _round(d),
        'A_on_I': {'min': _round(np.min(A_values)), 'max': _round(np.max(A_values))},
        'A_interval_on_I': {'min': _round(np.min(A_lo)), 'max': _round(np.max(A_hi))},
        'L_pieces': L,
        'baseline_axially_constant': baseline_constant,
        'baseline_N_span': _round(baseline_axial_span),
        'baseline_beta_span': _round(baseline_beta_span),
        'L_baseline_equals_baseline_beta': bool(
            baseline_constant and abs(L['baseline']['min'] - np.min(baseline[:, 1])) <= 1e-12),
        'geometry_killed': geometry_killed,
        'ell0_matter_change_exact': trial['scope']['ell0_matter_change_exact'],
        'ell0_baseline_vector': None,
        'ell0_lives_inside_L_baseline': True,
        'walk': walk,
        'walk_sign_stable': walk_sign_stable,
        'walk_L_min': _round(np.min(walk_mins)),
        'walk_L_max': _round(np.max(walk_maxs)),
        'assembled_interval_contains_zero': bool(
            assembled.get('interval_min', 0) <= 0 <= assembled.get('interval_max', 0)),
        'edge_can_eat_gap': edge_can_eat,
        'unenclosed_edge_tail': None,
        'locked_current_identified': certificate,
        'physical_NONEXISTENCE_certificate': certificate,
        'existence_tolerance': EXISTENCE,
        'scope': {
            'physical_local_gate': 'OPEN' if not certificate else 'NON-EXISTENCE',
            'leftover_v_dot_E_not_reused': True,
            'failed_optimization_is_not_NONEXISTENCE': True,
            'missing_tail_blocks_gap': missing_tail,
            'new_field_or_source_evolutions': 0,
            'physical_EXISTENCE_certificate': False,
            'metric_timestep': False,
        },
        'source_hashes': {p: digest(p) for p in OWNERS},
        'input_hashes': {
            str(ID.OUTPUT.relative_to(ROOT)): digest(ID.OUTPUT),
            str(I4.OUTPUT.relative_to(ROOT)): digest(I4.OUTPUT),
            str(BRANCH.relative_to(ROOT)): digest(BRANCH),
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
            raise FileExistsError('geometry-killing register exists; use --check')
        OUTPUT.write_text(json.dumps(result, sort_keys=True, indent=2, allow_nan=False)+'\n')
    elif result != json.loads(OUTPUT.read_text()):
        raise ValueError('geometry-killing replay differs')
    print(json.dumps({
        'status': result['status'],
        'physical_NONEXISTENCE_certificate': result['physical_NONEXISTENCE_certificate'],
        'walk_sign_stable': result['walk_sign_stable'],
        'walk_L_min': result['walk_L_min'],
        'walk_L_max': result['walk_L_max'],
        'geometry_killed': result['geometry_killed'],
        'L_pieces': {name: {k: row[k] for k in ('min', 'max', 'abs_min', 'sign_stable')}
                     for name, row in result['L_pieces'].items()},
        'edge_can_eat_gap': result['edge_can_eat_gap'],
    }, indent=2))
