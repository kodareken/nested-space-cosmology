#!/usr/bin/env python3
"""Dirac-free sign search for L(E) and the flux defect on the declared class.

The center is ad759424. The neighborhood is the n=128 coefficient tangents
at the leftover-null scales. The five saved n=16 histories are scored with
their stored residuals. No family is evolved, and the beta-primary ray is
not. A sign that holds only on this probe is not scoped NON-EXISTENCE.
"""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['VECLIB_MAXIMUM_THREADS'] = '1'
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np
from numpy.polynomial.polynomial import polyval

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
sys.path.insert(0, str(ROOT/'scripts'))
import derive_nsc_ks_flux_first_integral as FLUX
import derive_nsc_ks_n128_high_modes as N128
import derive_nsc_ks_n16_geometry_killing_current as KILL
import derive_nsc_ks_n16_leftover_anatomy as ANATOMY
from derive_nsc_evolved_incoming_constraints import coefficients_from_records
from recursive_horizons.nsc_evolved_incoming_constraints import (
    compatible_history_slots, surface_geometry_response)
from recursive_horizons.nsc_ks_profile_identity import profile_identity
from recursive_horizons.nsc_local_incoming_family import LocalIncomingFamily

OUTPUT = ROOT/'results/development/nsc-ks-class-current-neighborhood.json'
IMAGE = ROOT/'results/development/nsc-ks-n128-image-channels.json'
BRANCH = ROOT/'results/development/nsc-incoming-surface-regular-branch.json'
COEFFICIENT = ROOT/'results/development/nsc-incoming-surface-coefficients.json'
PARENT = N128.PARENT
SCALES = (-1e-3, 1e-3)
ROOM_FACTOR = 10.0
EXISTENCE = 3e-11
OWNERS = (
    'scripts/derive_nsc_ks_class_current_neighborhood.py',
    'docs/nsc-ks-class-current-neighborhood.md',
)
GEOMETRY_RAW = (2.0172752357439094e-13, 6.527672152811803e-11)
GEOMETRY_ENCLOSED = (2.1172752357439098e-13, 6.528672152811804e-11)


def digest(path):
    path = Path(path)
    return sha256((path if path.is_absolute() else ROOT/path).read_bytes()).hexdigest()


def empty_tangents(count):
    return np.zeros((0, count, 6))


def killing_on_residual(z, residual, slots, family, coeff, branch):
    imported = branch['branch']['imported_unchanged']
    A0_int = np.asarray(imported['A0'], float)
    d_int = np.asarray(imported['d'], float)
    A1_int = np.asarray(branch['branch']['intervals']['A1_total'], float)
    w = slots[:, 0]
    A_values = KILL.A_of_w(w, float(coeff['A0']), float(coeff['A1']))
    A_lo, A_hi = KILL.interval_A(w, A0_int, A1_int)
    ratio_lo, ratio_hi = KILL.interval_ratio(d_int, A_lo, A_hi)
    dN = KILL.axial_derivative(z, residual[:, 0], family.center, family.axial_inner)
    values = KILL.killing(residual, A_values, float(coeff['d']), dN)
    interval_min = float(np.min(
        residual[:, 1] - np.where(dN >= 0, ratio_hi, ratio_lo) * dN))
    interval_max = float(np.max(
        residual[:, 1] - np.where(dN >= 0, ratio_lo, ratio_hi) * dN))
    return values, interval_min, interval_max


def score(z, residual, slots, family, coeff, branch):
    values, interval_min, interval_max = killing_on_residual(
        z, residual, slots, family, coeff, branch)
    pointwise, integrated, divergence = FLUX.balance(slots, residual, z, coeff)
    return {
        'L_min': float(np.min(values)),
        'L_max': float(np.max(values)),
        'L_interval_min': float(interval_min),
        'L_interval_max': float(interval_max),
        'flux_min': float(np.min(pointwise)),
        'flux_max': float(np.max(pointwise)),
        'integrated': float(integrated),
        'div_min': float(np.min(divergence)),
        'div_max': float(np.max(divergence)),
    }


def span_of(rows, key_min, key_max=None):
    if key_max is None:
        values = np.array([row[key_min] for row in rows], float)
        lo = float(np.min(values))
        hi = float(np.max(values))
    else:
        lo = float(min(row[key_min] for row in rows))
        hi = float(max(row[key_max] for row in rows))
    return {
        'min': lo,
        'max': hi,
        'sign_stable': bool((lo > 0.0 and hi > 0.0) or (lo < 0.0 and hi < 0.0)),
        'abs_min': float(min(abs(lo), abs(hi))) if lo * hi > 0.0 else 0.0,
    }


def cancelling_u(slots, coeff):
    """Algebraic U that zeros the coded flux-divergence expression."""
    w, wz, wzz, _wzzz, _U, _Uz = np.asarray(slots, float).T
    a = float(coeff['a'])
    r = float(coeff['r'])
    Hr = float(coeff['Hr'])
    d = float(coeff['d'])
    A = float(coeff['A0']) + float(coeff['A1']) * w
    if np.any(A == 0.0) or not np.isfinite(A).all():
        raise ValueError('cancelling U needs A(w) nonzero on I')
    B = -(2.0 * A + (Hr + w / r) * d) / a ** 2
    Dval = polyval(w, coeff['D'])
    return -(B * wzz + float(coeff['C']) * wz * wz + Dval) / A


def relation_rows():
    return (
        ('L(E)', 'L_min', 'L_max', 'pointwise', True),
        ('Pi_prime_plus_E_beta', 'flux_min', 'flux_max', 'pointwise', True),
        ('integrated_Pi_plus_beta', 'integrated', None, 'scalar', True),
        ('flux_divergence', 'div_min', 'div_max', 'slot', False),
    )


def compute():
    parent = json.loads(PARENT.read_text())
    image = json.loads(IMAGE.read_text())
    branch = json.loads(BRANCH.read_text())
    coefficient = json.loads(COEFFICIENT.read_text())
    if parent['profile_identity'] != N128.PARENT_IDENTITY or not parent['newton_step_accepted']:
        raise ValueError('neighborhood search requires the accepted ad759424 history')
    if parent['solver'] != N128.SOLVER:
        raise ValueError('neighborhood search requires primal-block control')
    if image['stall'] != 'n128_sup_norm_image_outside_existence_ball':
        raise ValueError('image stall changed before the current search')
    if image['physical_NONEXISTENCE_certificate'] or image['probe_authorized']:
        raise ValueError('image register is not an open Dirac-free input')
    raw = tuple(image['channels']['geometry_between_node_raw'])
    enclosed = tuple(image['channels']['geometry_between_node_remainder'])
    if raw != GEOMETRY_RAW or enclosed != GEOMETRY_ENCLOSED:
        raise ValueError('geometry remainder binding changed')
    coeff = coefficients_from_records(coefficient, branch)
    family = N128.LocalIncomingFamily128(N128.lifted_coefficients(parent))
    if profile_identity(family, include_normal_window=True) != N128.LIFTED_IDENTITY:
        raise ValueError('n=128 family is not the zero-pad lift')
    z = np.asarray(parent['z'], float)
    gradient = np.asarray(parent['action_gradient'], float)
    metric = family.metric()
    slots, tangents = compatible_history_slots(metric, z, len(metric.directions))
    if tangents.shape[0] != 256:
        raise ValueError('declared n=128 class no longer has 256 tangents')
    base_geometry = surface_geometry_response(
        slots, empty_tangents(len(z)), coeff)['action_gradient_change']
    center = score(z, gradient, slots, family, coeff, branch)
    center['kind'] = 'center'
    raw = np.asarray(parent['raw_action_gradient'], float)
    matter = np.asarray(parent['summed_matter_change'], float)
    edge = np.asarray(parent['source_cutoff_edge_gradient'], float)
    baseline = raw - base_geometry - matter
    if np.max(np.abs(baseline + base_geometry + matter - raw)) > 1e-12:
        raise ValueError('baseline plus geometry plus matter left the raw residual')
    piece_fields = {
        'baseline': baseline,
        'geometry': base_geometry,
        'matter': matter,
        'edge': edge,
        'assembled': gradient,
    }
    pieces = {
        name: score(z, field, slots, family, coeff, branch)
        for name, field in piece_fields.items()}
    rows = [center]
    for path in ANATOMY.WALK:
        record = json.loads(path.read_text())
        hist = LocalIncomingFamily(np.asarray(record['history']['coefficients'], float))
        hist_z = np.asarray(record['z'], float)
        hist_slots, _hist_tangent = compatible_history_slots(hist.metric(), hist_z, 0)
        row = score(
            hist_z,
            np.asarray(record['action_gradient'], float),
            hist_slots,
            hist,
            coeff,
            branch,
        )
        row['kind'] = 'walk'
        row['profile_identity'] = record['profile_identity']
        rows.append(row)
    for direction, tangent in enumerate(tangents):
        for scale in SCALES:
            perturbed = slots + scale * tangent
            shifted = surface_geometry_response(
                perturbed, empty_tangents(len(z)), coeff)['action_gradient_change']
            residual = gradient - base_geometry + shifted
            row = score(z, residual, perturbed, family, coeff, branch)
            row['kind'] = 'neighborhood'
            row['direction'] = int(direction)
            row['scale'] = float(scale)
            rows.append(row)
    substituted = slots.copy()
    substituted[:, 4] = cancelling_u(slots, coeff)
    substituted_divergence = FLUX.flux_divergence(substituted, coeff)
    divergence_after_U = {
        'max_abs_U_change': float(np.max(np.abs(substituted[:, 4] - slots[:, 4]))),
        'div_min': float(np.min(substituted_divergence)),
        'div_max': float(np.max(substituted_divergence)),
        'abs_max': float(np.max(np.abs(substituted_divergence))),
    }
    if divergence_after_U['abs_max'] > 1e-9:
        raise ValueError('documented flux divergence did not zero under its algebraic U')
    length = float(z[-1] - z[0])
    if length <= 0.0:
        raise ValueError('declared interval needs positive length')
    beta_room = float(np.nextafter(GEOMETRY_ENCLOSED[1] * ROOM_FACTOR, np.inf))
    integral_room = float(np.nextafter(GEOMETRY_ENCLOSED[1] * length * ROOM_FACTOR, np.inf))
    relations = {}
    for name, key_min, key_max, kind, necessary in relation_rows():
        span = span_of(rows, key_min, key_max)
        threshold = integral_room if kind == 'scalar' else beta_room
        if kind == 'slot':
            threshold = beta_room
            necessary = False
        gap_survives = bool(span['sign_stable'] and span['abs_min'] > threshold)
        relations[name] = {
            **span,
            'kind': kind,
            'necessary_when_constraints_vanish': necessary,
            'survival_threshold': threshold,
            'gap_survives_geometry_remainder': gap_survives,
            'class_wide': False,
        }
    relations['L(E)']['interval_min'] = float(min(row['L_interval_min'] for row in rows))
    relations['L(E)']['interval_max'] = float(max(row['L_interval_max'] for row in rows))
    relations['L(E)']['interval_excludes_zero'] = bool(
        relations['L(E)']['interval_min'] * relations['L(E)']['interval_max'] > 0.0)
    relations['L(E)']['axial_derivative_error_enclosed'] = False
    center_sign = np.sign(center['integrated'])
    crossings = []
    for direction in range(tangents.shape[0]):
        minus = next(
            row['integrated'] for row in rows
            if row.get('direction') == direction and row.get('scale') == SCALES[0])
        plus = next(
            row['integrated'] for row in rows
            if row.get('direction') == direction and row.get('scale') == SCALES[1])
        slope = (plus - minus) / (SCALES[1] - SCALES[0])
        if slope == 0.0 or not np.isfinite(slope):
            continue
        crossings.append(float(-center['integrated'] / slope))
    edge = np.asarray(parent['source_cutoff_edge_gradient'], float)
    edge_beta_integral = float(np.trapezoid(edge[:, 1], z))
    certificate = False
    sign_flips = [
        name for name, row in relations.items()
        if row['necessary_when_constraints_vanish'] and not row['sign_stable']]
    subclass = [
        name for name, row in relations.items()
        if row['sign_stable'] and not row['class_wide']]
    if 'L(E)' not in sign_flips or 'Pi_prime_plus_E_beta' not in sign_flips:
        raise ValueError('pointwise current changed stability; refuse a silent certificate')
    if 'integrated_Pi_plus_beta' not in sign_flips:
        raise ValueError('integrated defect stayed sign-stable; refuse a silent certificate')
    saved_rows = [row for row in rows if row['kind'] in ('center', 'walk')]
    neighborhood_rows = [row for row in rows if row['kind'] == 'neighborhood']
    saved_integrated = span_of(saved_rows, 'integrated')
    if not saved_integrated['sign_stable'] or saved_integrated['min'] <= 0.0:
        raise ValueError('saved histories no longer share a positive integrated defect')
    witness = min(
        (row for row in neighborhood_rows if row['integrated'] < 0.0),
        key=lambda row: (row['direction'], row['scale']))
    if not pieces['baseline']['L_min'] < 0.0 or not pieces['baseline']['L_max'] < 0.0:
        raise ValueError('L(baseline) is no longer negative on ad759424')
    if pieces['geometry']['L_min'] * pieces['geometry']['L_max'] > 0.0:
        raise ValueError('geometry piece of L(E) no longer changes sign')
    return {
        'schema': 'NSC-KS-CLASS-CURRENT-NEIGHBORHOOD-v1',
        'accountable_author': 'Douglas Ek',
        'status': (
            'OPEN: L(E) and Pi-prime plus E_beta change sign on the declared '
            'probe; no scoped NON-EXISTENCE'),
        'stall': 'n128_sup_norm_image_outside_existence_ball',
        'search_result': 'class_current_neighborhood_not_sign_stable',
        'profile_identity': N128.PARENT_IDENTITY,
        'lifted_profile_identity': N128.LIFTED_IDENTITY,
        'families_evolved': 0,
        'probe_authorized': False,
        'beta_primary_evolved': False,
        'higher_n_authorized': False,
        'newton_step_accepted': False,
        'new_best_residual': False,
        'physical_EXISTENCE_certificate': False,
        'physical_NONEXISTENCE_certificate': certificate,
        'existence_tolerance': EXISTENCE,
        'neighborhood_scales': list(SCALES),
        'direction_count': int(tangents.shape[0]),
        'room_factor': ROOM_FACTOR,
        'geometry_between_node_raw': list(GEOMETRY_RAW),
        'geometry_between_node_remainder': list(GEOMETRY_ENCLOSED),
        'interval_length': length,
        'pointwise_survival_threshold': beta_room,
        'integral_survival_threshold': integral_room,
        'center': {key: center[key] for key in (
            'L_min', 'L_max', 'L_interval_min', 'L_interval_max',
            'flux_min', 'flux_max', 'integrated', 'div_min', 'div_max')},
        'center_pieces': pieces,
        'saved_history_integrated': saved_integrated,
        'integrated_flip_witness': {
            'direction': witness['direction'],
            'scale': witness['scale'],
            'integrated': witness['integrated'],
            'L_min': witness['L_min'],
            'L_max': witness['L_max'],
            'flux_min': witness['flux_min'],
            'flux_max': witness['flux_max'],
        },
        'walk': [
            {
                'profile_identity': row['profile_identity'],
                'L_min': row['L_min'],
                'L_max': row['L_max'],
                'flux_min': row['flux_min'],
                'flux_max': row['flux_max'],
                'integrated': row['integrated'],
                'div_min': row['div_min'],
                'div_max': row['div_max'],
            }
            for row in rows if row['kind'] == 'walk'],
        'relations': relations,
        'sign_flip_relations': sign_flips,
        'subclass_only_relations': subclass,
        'integrated_linear_crossing_scales': {
            'count': len(crossings),
            'min_abs': None if not crossings else float(np.min(np.abs(crossings))),
            'center_sign': float(center_sign),
        },
        'edge_beta_integral': edge_beta_integral,
        'edge_can_eat_integrated_gap': bool(
            abs(edge_beta_integral) >= abs(center['integrated'])),
        'flux_divergence_zeroed_by_algebraic_U': divergence_after_U,
        'class_argument': (
            'The probe is the zero-pad lift, its 256 coefficient tangents at '
            'scales -1e-3 and 1e-3 with matter held fixed, and the five saved '
            'n=16 residuals. L(E) and Pi-prime plus E_beta change sign on '
            'ad759424, which is inside the declared class. L(baseline) stays '
            'negative there, while the geometry, matter and edge pieces of '
            'L(E) change sign and are larger, so the baseline piece does not '
            'control the assembled residual. The six saved histories keep a '
            'positive integrated defect, and the w-mode tangent named in '
            'integrated_flip_witness changes that sign inside the same '
            'scales. The slot divergence does not use the residual; '
            'substituting the algebraic U of the coded expression drops it '
            'to roundoff.'),
        'next_owner_action': (
            'Do not evolve the beta-primary ray and do not open n=256. '
            'L(E), Pi-prime plus E_beta, and the slot divergence are spent '
            'on this probe. A later NON-EXISTENCE certificate still needs a '
            'different necessary relation with an error-controlled gap on '
            'every declared (w,U) history.'),
        'samples': rows,
        'source_hashes': {path: digest(path) for path in OWNERS},
        'input_hashes': {
            str(PARENT.relative_to(ROOT)): digest(PARENT),
            str(IMAGE.relative_to(ROOT)): digest(IMAGE),
            str(BRANCH.relative_to(ROOT)): digest(BRANCH),
            str(COEFFICIENT.relative_to(ROOT)): digest(COEFFICIENT),
            **{str(path.relative_to(ROOT)): digest(path) for path in ANATOMY.WALK},
            'scripts/derive_nsc_ks_flux_first_integral.py': digest(
                'scripts/derive_nsc_ks_flux_first_integral.py'),
            'scripts/derive_nsc_ks_n16_geometry_killing_current.py': digest(
                'scripts/derive_nsc_ks_n16_geometry_killing_current.py'),
            'scripts/derive_nsc_ks_n128_high_modes.py': digest(
                'scripts/derive_nsc_ks_n128_high_modes.py'),
        },
    }


def check():
    saved = json.loads(OUTPUT.read_text())
    fresh = json.loads(json.dumps(compute(), sort_keys=True))
    saved_norm = json.loads(json.dumps(saved, sort_keys=True))
    if saved_norm != fresh:
        raise ValueError('class-current neighborhood register changed')
    if saved['physical_NONEXISTENCE_certificate'] or saved['families_evolved']:
        raise ValueError('this probe must stay OPEN and unevolved')
    if saved['beta_primary_evolved'] or saved['higher_n_authorized']:
        raise ValueError('beta-primary evolution and n=256 stay unauthorized')
    return saved


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    if args.write:
        if OUTPUT.exists():
            raise SystemExit('class-current neighborhood register already exists')
        payload = compute()
        OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+'\n')
        brief = {
            name: {
                'min': row['min'],
                'max': row['max'],
                'sign_stable': row['sign_stable'],
                'gap_survives_geometry_remainder': row['gap_survives_geometry_remainder'],
            }
            for name, row in payload['relations'].items()
        }
        print(json.dumps({
            'search_result': payload['search_result'],
            'physical_NONEXISTENCE_certificate': payload['physical_NONEXISTENCE_certificate'],
            'relations': brief,
            'integrated_linear_crossing_scales': payload['integrated_linear_crossing_scales'],
            'flux_divergence_zeroed_by_algebraic_U': payload[
                'flux_divergence_zeroed_by_algebraic_U'],
            'edge_can_eat_integrated_gap': payload['edge_can_eat_integrated_gap'],
        }, indent=2))
    else:
        saved = check()
        print(json.dumps({
            'search_result': saved['search_result'],
            'physical_NONEXISTENCE_certificate': saved['physical_NONEXISTENCE_certificate'],
            'sign_flip_relations': saved['sign_flip_relations'],
        }, indent=2))
