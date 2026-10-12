#!/usr/bin/env python3
"""Channels and sup-norm image of the n=128 truncation at ad759424.

No Dirac evolution. Geometry between-node is evaluated on I. UV, field,
low/subgap and the full residual between-node bound stay None when the
owning tool does not accept this history. The stable-rank sup-norm image
is screened for a beta-primary step. A trivial cokernel and a positive
joint linear gain are not scoped NON-EXISTENCE.
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

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
sys.path.insert(0, str(ROOT/'scripts'))
import derive_nsc_ks_declared_wu_cokernel as COK
import derive_nsc_ks_flux_first_integral as FLUX
import derive_nsc_ks_n128_high_modes as N128
import derive_nsc_ks_n16_geometry_killing_current as KILL
from derive_nsc_ks_coupled_newton_n64_truncated_iterate import LocalIncomingFamily64
from derive_nsc_evolved_incoming_constraints import coefficients_from_records
from recursive_horizons.nsc_evolved_incoming_constraints import (
    compatible_history_slots, surface_geometry_response)
from recursive_horizons.nsc_ks_profile_identity import profile_identity
from recursive_horizons.nsc_local_history_newton import (
    LocalHistoryNewtonSettings, _linear_step)

PARENT = N128.PARENT
SOURCE = N128.OUTPUT
OUTPUT = ROOT/'results/development/nsc-ks-n128-image-channels.json'
BRANCH = ROOT/'results/development/nsc-incoming-surface-regular-branch.json'
COEFFICIENT = ROOT/'results/development/nsc-incoming-surface-coefficients.json'
PRINCIPAL = COK.PRINCIPAL
FINE = ROOT/'results/development/nsc-ks-fine-matter-accuracy.json'
FINE_TRAJECTORY = ROOT/'results/development/nsc-ks-fine-trajectory.json'
DENSE_COUNT = 129
EXISTENCE = 3e-11
RATIO = 1e6
BETA_HS = (1e-6, 1e-4, 1e-3)
JOINT_HS = (1e-8, 1e-6, 1e-4)
BALL_HS = (1e-6, 1e-4, 1e-2)
REQUIRED_BETA_DROP = 1e-8
OWNERS = (
    'scripts/derive_nsc_ks_n128_image_channels.py',
    'docs/nsc-ks-n128-image-channels.md',
)


def digest(path):
    path = Path(path)
    return sha256((path if path.is_absolute() else ROOT/path).read_bytes()).hexdigest()


def pair(values):
    return [float(values[0]), float(values[1])]


def barycentric_weights(nodes):
    """Scaled barycentric weights in node order. Pure products, no LAPACK."""
    nodes = np.asarray(nodes, np.float64)
    weights = np.empty(len(nodes), np.float64)
    for index, node in enumerate(nodes):
        product = 1.0
        for other_index, other in enumerate(nodes):
            if other_index != index:
                product *= (node - other)
        weights[index] = 1.0 / product
    return weights / np.max(np.abs(weights))


def barycentric_values(nodes, values, weights, query):
    values = np.asarray(values, np.float64)
    query = np.asarray(query, np.float64)
    out = np.empty(len(query), np.float64)
    for index, point in enumerate(query):
        numerator = 0.0
        denominator = 0.0
        exact = None
        for node_index, node in enumerate(nodes):
            difference = point - node
            if difference == 0.0:
                exact = values[node_index]
                break
            term = weights[node_index] / difference
            numerator += term * values[node_index]
            denominator += term
        out[index] = exact if exact is not None else numerator / denominator
    return out


def slots_and_geometry(family, z, coeff):
    slots, _tangent = compatible_history_slots(family.metric(), z, 0)
    empty = np.zeros((0, len(z), 6))
    value = surface_geometry_response(slots, empty, coeff)['action_gradient_change']
    return slots, np.asarray(value, float)


def geometry_remainder(family, z, coeff):
    _slots, geo_nodes = slots_and_geometry(family, z, coeff)
    dense = family.collocation_nodes(DENSE_COUNT)
    if dense[0] < family.interval[0] - 1e-15 or dense[-1] > family.interval[1] + 1e-15:
        raise ValueError('dense nodes must stay inside the declared interval I')
    _dense_slots, geo_dense = slots_and_geometry(family, dense, coeff)
    weights = barycentric_weights(z)
    geo_interp = np.column_stack([
        barycentric_values(z, geo_nodes[:, k], weights, dense) for k in range(2)])
    raw = np.max(np.abs(geo_dense - geo_interp), axis=0)
    if np.any(~np.isfinite(raw)) or np.any(raw < 0):
        raise ValueError('finite nonnegative geometry interpolation remainder required')
    if np.all(raw < 1e-13):
        enclosed = np.full(2, 1e-13)
    else:
        enclosed = np.nextafter(raw + 1e-14, np.inf)
    return raw, enclosed


def stable_basis(matrix):
    row = np.linalg.norm(matrix, axis=1)
    col = np.linalg.norm(matrix, axis=0)
    row = np.where(row > 0, row, 1.0)
    col = np.where(col > 0, col, 1.0)
    _left, spectrum, right = np.linalg.svd(
        matrix / row[:, None] / col[None, :], full_matrices=False)
    kept = np.flatnonzero(spectrum[0] / spectrum < RATIO)
    if kept.size == 0 or int(kept[0]) != 0 or int(kept[-1]) != kept.size - 1:
        raise ValueError('stable singular values are not a leading block')
    last = int(kept[-1]) + 1
    basis = N128.whiten(right, col, 0, last)
    return basis, last, pair_spectrum(spectrum)


def pair_spectrum(spectrum):
    spectrum = np.asarray(spectrum, float)
    return {
        's0': float(spectrum[0]),
        'smin': float(spectrum[-1]),
        'ratio': float(spectrum[0] / spectrum[-1]),
        'count': int(spectrum.size),
    }


def max_beta_gain(response, residual, measured, h, n_cap):
    lo = 0.0
    hi = float(measured[1])
    best = None
    for _step in range(48):
        mid = 0.5 * (lo + hi)
        ok, direction = N128.feasible(
            response, residual, n_cap, float(measured[1]) - mid, h)
        if ok:
            lo = mid
            best = direction
        else:
            hi = mid
    return lo, best


def screen(response, residual, n_cap, beta_cap, h_values):
    return [{
        'h': float(h),
        'feasible': bool(N128.feasible(response, residual, n_cap, beta_cap, h)[0]),
    } for h in h_values]


def pointwise_principal(slots, coeff):
    geometry = surface_geometry_response(
        slots, np.zeros((0, len(slots), 6)), coeff)
    jac = geometry['jacobian']
    dets = []
    for index in range(len(slots)):
        matrix = np.array([
            [jac[index, 0, 4], jac[index, 0, 2]],
            [jac[index, 1, 5], jac[index, 1, 3]],
        ], float)
        dets.append(float(np.linalg.det(matrix)))
    dets = np.asarray(dets, float)
    return float(np.min(dets)), float(np.max(dets))


def killing_on_residual(z, residual, slots, family, coeff, branch):
    imported = branch['branch']['imported_unchanged']
    A0_int = np.asarray(imported['A0'], float)
    d_int = np.asarray(imported['d'], float)
    A1_int = np.asarray(branch['branch']['intervals']['A1_total'], float)
    A0 = float(coeff['A0'])
    A1 = float(coeff['A1'])
    d = float(coeff['d'])
    w = slots[:, 0]
    A_values = KILL.A_of_w(w, A0, A1)
    A_lo, A_hi = KILL.interval_A(w, A0_int, A1_int)
    ratio_lo, ratio_hi = KILL.interval_ratio(d_int, A_lo, A_hi)
    dN = KILL.axial_derivative(z, residual[:, 0], family.center, family.axial_inner)
    values = KILL.killing(residual, A_values, d, dN)
    interval_min = float(np.min(
        residual[:, 1] - np.where(dN >= 0, ratio_hi, ratio_lo) * dN))
    interval_max = float(np.max(
        residual[:, 1] - np.where(dN >= 0, ratio_lo, ratio_hi) * dN))
    return {
        'min': float(np.min(values)),
        'max': float(np.max(values)),
        'sign_stable': bool(np.all(values > 0) or np.all(values < 0)),
        'interval_min': interval_min,
        'interval_max': interval_max,
        'interval_excludes_zero': bool(interval_min * interval_max > 0),
    }


def compute():
    parent = json.loads(PARENT.read_text())
    saved = json.loads(SOURCE.read_text())
    branch = json.loads(BRANCH.read_text())
    coefficient = json.loads(COEFFICIENT.read_text())
    principal = json.loads(PRINCIPAL.read_text())
    fine = json.loads(FINE.read_text())
    fine_trajectory = json.loads(FINE_TRAJECTORY.read_text())
    if parent['profile_identity'] != N128.PARENT_IDENTITY or not parent['newton_step_accepted']:
        raise ValueError('image register requires the accepted ad759424 history')
    if parent['solver'] != N128.SOLVER or saved['solver'] != N128.SOLVER:
        raise ValueError('image register requires primal-block control')
    if saved['profile_identity'] != N128.LIFTED_IDENTITY:
        raise ValueError('lifted n=128 identity changed')
    if saved['probe_authorized'] or saved['newton_step_accepted']:
        raise ValueError('n=128 column owner must not already authorize a step')
    coeff = coefficients_from_records(coefficient, branch)
    family64 = LocalIncomingFamily64(np.asarray(parent['history']['coefficients'], float))
    family128 = N128.LocalIncomingFamily128(N128.lifted_coefficients(parent))
    if profile_identity(family64, include_normal_window=True) != N128.PARENT_IDENTITY:
        raise ValueError('n=64 family is not ad759424')
    if profile_identity(family128, include_normal_window=True) != N128.LIFTED_IDENTITY:
        raise ValueError('n=128 family is not the zero-pad lift')
    z = np.asarray(parent['z'], float)
    gradient = np.asarray(parent['action_gradient'], float)
    measured = np.asarray(saved['best_constraint_maxima'], float)
    if not np.allclose(np.max(np.abs(gradient), axis=0), measured, rtol=0, atol=0):
        raise ValueError('saved maxima left the parent residual')
    jacobian = np.asarray(saved['history_jacobian'], float)
    if jacobian.shape != (256, len(z), 2) or gradient.shape != (len(z), 2):
        raise ValueError('n=128 Jacobian shape changed')
    slots64, geo64 = slots_and_geometry(family64, z, coeff)
    slots128, geo128 = slots_and_geometry(family128, z, coeff)
    if np.max(np.abs(slots64 - slots128)) != 0.0 or np.max(np.abs(geo64 - geo128)) != 0.0:
        raise ValueError('zero-pad lift changed the owned slots or geometry')
    raw64, enclosed64 = geometry_remainder(family64, z, coeff)
    raw128, enclosed128 = geometry_remainder(family128, z, coeff)
    if np.max(np.abs(raw64 - raw128)) != 0.0 or np.max(np.abs(enclosed64 - enclosed128)) != 0.0:
        raise ValueError('zero-pad lift changed the geometry between-node remainder')
    _A, _B, _d, _e, det = COK.principal_blocks(slots64, coeff)
    det_min, det_max = float(np.min(det)), float(np.max(det))
    certified = np.asarray(principal['matrix']['determinant_interval'], float)
    point_min, point_max = pointwise_principal(slots64, coeff)
    det_excludes_zero = bool(det_min * det_max > 0 and point_min * point_max > 0)
    same_sign = bool(det_min * certified[0] > 0 and point_min * certified[0] > 0)
    if not det_excludes_zero or not same_sign:
        raise ValueError('principal cokernel at ad759424 is no longer the trivial certified sign')
    residual = np.concatenate([gradient[:, 0], gradient[:, 1]])
    matrix = np.vstack([jacobian[:, :, 0].T, jacobian[:, :, 1].T])
    raw_spectrum = pair_spectrum(np.linalg.svd(matrix, compute_uv=False))
    basis, stable_rank, equilibrated = stable_basis(matrix)
    response = matrix @ basis
    delta, diagnostics = _linear_step(
        matrix, residual, matrix, LocalHistoryNewtonSettings())
    beta_rows = []
    for h in BETA_HS:
        _gain, direction = max_beta_gain(response, residual, measured, h, float(measured[0]))
        if direction is None:
            raise ValueError('beta-primary screen lost the zero step')
        step = h * (basis @ direction)
        predicted = N128.predicted_maxima(residual, matrix, step)
        beta_rows.append({
            'h': float(h),
            'predicted_maxima': pair(predicted),
            'n_gain': float(measured[0] - predicted[0]),
            'beta_gain': float(measured[1] - predicted[1]),
            'step_norm': float(np.linalg.norm(step)),
            'n_strictly_improved': bool(predicted[0] < measured[0]),
        })
    floor = min(beta_rows, key=lambda row: row['predicted_maxima'][1])
    joint_rows = []
    for h in JOINT_HS:
        gain, _direction = N128.max_n_gain(
            response, residual, measured, h, REQUIRED_BETA_DROP)
        joint_rows.append({'h': float(h), 'n_gain': float(gain), 'required_beta_drop': REQUIRED_BETA_DROP})
    both_3e11 = screen(response, residual, EXISTENCE, EXISTENCE, BALL_HS)
    both_1e9 = screen(response, residual, 1e-9, 1e-9, BALL_HS)
    beta_3e11_n_held = screen(response, residual, float(measured[0]), EXISTENCE, BALL_HS)
    beta_1e9_n_held = screen(response, residual, float(measured[0]), 1e-9, BALL_HS)
    useful = screen(response, residual, float(measured[0]) - 1e-12, 2e-9, (1e-4,))
    ball_closed = (
        not any(row['feasible'] for row in both_3e11)
        and not any(row['feasible'] for row in both_1e9)
        and not any(row['feasible'] for row in beta_3e11_n_held)
        and not any(row['feasible'] for row in beta_1e9_n_held))
    if not ball_closed or floor['predicted_maxima'][1] <= EXISTENCE:
        raise ValueError('stable image reached a ball this register treats as closed')
    if floor['n_strictly_improved'] or floor['beta_gain'] <= 0.0:
        raise ValueError('beta-primary floor does not move beta while holding N')
    n_abs = np.abs(gradient[:, 0])
    beta_abs = np.abs(gradient[:, 1])
    killing = killing_on_residual(z, gradient, slots64, family64, coeff, branch)
    dpi_plus_beta, integrated, divergence = FLUX.balance(slots64, gradient, z, coeff)
    flux_values = KILL.summary(dpi_plus_beta)
    divergence_values = KILL.summary(divergence)
    fine_identity = fine_trajectory['axial_profile_identity']
    if fine_identity == N128.PARENT_IDENTITY:
        raise ValueError('field-accuracy owner unexpectedly matches ad759424')
    if fine['input_hashes']['results/development/nsc-ks-fine-trajectory.json'] != digest(FINE_TRAJECTORY):
        raise ValueError('fine-matter trajectory binding changed')
    n_short = float(floor['predicted_maxima'][0] - EXISTENCE)
    beta_short = float(floor['predicted_maxima'][1] - EXISTENCE)
    joint_gain = float(joint_rows[-1]['n_gain'])
    certificate = False
    return {
        'schema': 'NSC-KS-N128-IMAGE-CHANNELS-v1',
        'accountable_author': 'Douglas Ek',
        'status': (
            'OPEN: n=128 sup-norm image stays outside the 3e-11 ball; '
            'cokernel trivial; not NON-EXISTENCE'),
        'stall': 'n128_sup_norm_image_outside_existence_ball',
        'profile_identity': N128.PARENT_IDENTITY,
        'lifted_profile_identity': N128.LIFTED_IDENTITY,
        'best_constraint_maxima': pair(measured),
        'existence_tolerance': EXISTENCE,
        'n_gap': float(measured[0] - EXISTENCE),
        'beta_gap': float(measured[1] - EXISTENCE),
        'channels': {
            'geometry_between_node_raw': pair(raw64),
            'geometry_between_node_remainder': pair(enclosed64),
            'geometry_remainder_shared_by_zero_pad': True,
            'full_between_node_remainder': None,
            'changed_history_UV_tail': None,
            'field_error_bound': None,
            'baseline_low_subgap': None,
            'missing_paths': {
                'full_between_node_remainder': (
                    'saved operators are the 47-node target; dense matter is not stored'),
                'changed_history_UV_tail': (
                    'linear UV contraction is a coefficient identity with no numerical tail on I'),
                'field_error_bound': (
                    'fine-matter accuracy is bound to axial profile '
                    + fine_identity + ', not ad759424'),
                'baseline_low_subgap': (
                    'group13 subgap panel and the Pauli matrix identity are not a '
                    'constraint error at ad759424'),
            },
        },
        'cokernel': {
            'principal_det_min': det_min,
            'principal_det_max': det_max,
            'pointwise_principal_det_min': point_min,
            'pointwise_principal_det_max': point_max,
            'certified_determinant_interval': [float(certified[0]), float(certified[1])],
            'principal_det_excludes_zero': det_excludes_zero,
            'same_sign_as_certified_det': same_sign,
            'continuous_cokernel_trivial': True,
            'principal_det_is_invertibility_not_existence': True,
        },
        'jacobian_spectrum': {
            'raw': raw_spectrum,
            'equilibrated': equilibrated,
            'stable_rank': stable_rank,
            'singular_ratio_cut': RATIO,
            'whitened_width': int(basis.shape[1]),
            'linear_step_delta_is_none': delta is None,
            'linear_step_rank': int(diagnostics['rank']),
            'linear_step_n_unknowns': int(diagnostics['n_unknowns']),
            'linear_step_n_equations': int(diagnostics['n_equations']),
            'linear_step_rank_or_conditioning': bool(diagnostics['rank_or_conditioning']),
            'left_cokernel_trivial': bool(
                equilibrated['smin'] > 0.0 and equilibrated['count'] == matrix.shape[0]),
            'fat_system_not_a_left_null': bool(
                delta is None and diagnostics['rank'] < diagnostics['n_unknowns']),
        },
        'nodal_flatness': {
            'n_second_over_first': float(np.partition(n_abs, -2)[-2] / np.max(n_abs)),
            'beta_second_over_first': float(np.partition(beta_abs, -2)[-2] / np.max(beta_abs)),
        },
        'beta_primary_rows': beta_rows,
        'beta_primary_floor': floor,
        'image_short_of_existence_ball': {'N': n_short, 'beta': beta_short},
        'joint_n_gain_at_beta_drop_1e-8': joint_rows,
        'ball_screens': {
            'both_under_3e-11': both_3e11,
            'both_under_1e-9': both_1e9,
            'beta_under_3e-11_n_held': beta_3e11_n_held,
            'beta_under_1e-9_n_held': beta_1e9_n_held,
            'n_down_1e-12_and_beta_at_most_2e-9': useful,
        },
        'killing_current': killing,
        'flux_first_integral': {
            'dpi_plus_beta': flux_values,
            'integrated_pi_plus_beta': float(integrated),
            'flux_divergence': divergence_values,
        },
        'joint_linear_gain_positive': bool(joint_gain > 0.0),
        'probe_authorized': False,
        'newton_step_accepted': False,
        'new_best_residual': False,
        'families_evolved': 0,
        'physical_EXISTENCE_certificate': False,
        'physical_NONEXISTENCE_certificate': certificate,
        'higher_n_authorized': False,
        'higher_n_note': 'n=256 is the same subclass pattern and is not authorized by this image',
        'newton_rank_note': (
            'Rank 90 of 94 because the solver floor is 1e-14 when the leading '
            'singular value is below 1. All 94 equilibrated singular values stay '
            'positive. The absent step is the fat shape: 256 unknowns, 94 equations.'),
        'extra_certificate': (
            'Scoped NON-EXISTENCE needs a necessary relation on the evolved state '
            'with an error-controlled gap above zero for every history in the '
            'declared (w,U) class on I. The truncation image gap is representation '
            'evidence. The geometry principal cokernel is trivial. The Fourier '
            'leftover-null is a finite-span vector of an older geometry Jacobian. '
            'This coupled Jacobian has 94 positive equilibrated singular values; '
            'Newton rank 90 is the absolute floor, not a certified left null. '
            'L(E) changes sign on this residual. Pi-prime plus E_beta changes sign. '
            'The slot flux divergence is positive here and does not use the residual. '
            'None of these is the class relation.'),
        'source_hashes': {path: digest(path) for path in OWNERS},
        'input_hashes': {
            str(PARENT.relative_to(ROOT)): digest(PARENT),
            str(SOURCE.relative_to(ROOT)): digest(SOURCE),
            str(BRANCH.relative_to(ROOT)): digest(BRANCH),
            str(COEFFICIENT.relative_to(ROOT)): digest(COEFFICIENT),
            str(PRINCIPAL.relative_to(ROOT)): digest(PRINCIPAL),
            str(FINE.relative_to(ROOT)): digest(FINE),
            str(FINE_TRAJECTORY.relative_to(ROOT)): digest(FINE_TRAJECTORY),
            'scripts/derive_nsc_ks_n128_high_modes.py': digest(
                'scripts/derive_nsc_ks_n128_high_modes.py'),
            'scripts/derive_nsc_ks_n16_geometry_killing_current.py': digest(
                'scripts/derive_nsc_ks_n16_geometry_killing_current.py'),
            'scripts/derive_nsc_ks_flux_first_integral.py': digest(
                'scripts/derive_nsc_ks_flux_first_integral.py'),
            'scripts/derive_nsc_ks_declared_wu_cokernel.py': digest(
                'scripts/derive_nsc_ks_declared_wu_cokernel.py'),
        },
    }


def check():
    saved = json.loads(OUTPUT.read_text())
    fresh = json.loads(json.dumps(compute(), sort_keys=True))
    saved_norm = json.loads(json.dumps(saved, sort_keys=True))
    if saved_norm != fresh:
        raise ValueError('n=128 image-channel register changed')
    if saved['physical_NONEXISTENCE_certificate'] or saved['probe_authorized']:
        raise ValueError('image register must stay OPEN without a probe')
    return saved


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    if args.write:
        if OUTPUT.exists():
            raise SystemExit('image-channel register already exists')
        payload = compute()
        OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+'\n')
        print(json.dumps({
            'stall': payload['stall'],
            'geometry_between_node_raw': payload['channels']['geometry_between_node_raw'],
            'geometry_between_node_remainder': payload['channels']['geometry_between_node_remainder'],
            'beta_primary_floor': payload['beta_primary_floor'],
            'joint_n_gain_at_beta_drop_1e-8': payload['joint_n_gain_at_beta_drop_1e-8'],
            'killing_current': payload['killing_current'],
            'flux_sign_stable': payload['flux_first_integral']['dpi_plus_beta']['sign_stable'],
            'left_cokernel_trivial': payload['jacobian_spectrum']['left_cokernel_trivial'],
            'stable_rank': payload['jacobian_spectrum']['stable_rank'],
        }, indent=2))
    else:
        saved = check()
        print(json.dumps({
            'stall': saved['stall'],
            'beta_primary_floor': saved['beta_primary_floor']['predicted_maxima'],
        }, indent=2))
