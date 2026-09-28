#!/usr/bin/env python3
"""Assembled Fourier-on-I EXISTENCE search at iterate4. No new matter columns."""
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
import derive_nsc_ks_coupled_newton_n16_damped as D
import derive_nsc_ks_coupled_newton_n16_damped_iterate4 as I4
import derive_nsc_ks_coupled_newton_n16_damped_next5 as N5
import derive_nsc_ks_n16_leftover_anatomy as A
import derive_nsc_ks_n32_geometry_leftover as N32
from recursive_horizons.nsc_evolved_incoming_constraints import compatible_history_slots
from recursive_horizons.nsc_local_history_newton import LocalHistoryNewtonSettings
from recursive_horizons.nsc_local_incoming_fourier_family import LocalFourierIncomingFamily

OUTPUT = ROOT/'results/development/nsc-ks-fourier-assembled-search.json'
OWNERS = (
    'scripts/derive_nsc_ks_fourier_assembled_search.py',
    'docs/nsc-ks-fourier-assembled-search.md',
    'src/recursive_horizons/nsc_local_incoming_fourier_family.py',
)
TRANSFER_RELATIVE = 1e-8
FORBIDDEN = D.FORBIDDEN_IDENTITIES + (I4.EXPECTED_IDENTITY,)


def digest(path):
    path = Path(path)
    return sha256((path if path.is_absolute() else ROOT/path).read_bytes()).hexdigest()


def _round(value, digits=12):
    return N32._round(value, digits)


def flatten_slots(tangents):
    tangents = np.asarray(tangents, float)
    return tangents.reshape(len(tangents), -1).T


def transfer_assembled(cheb_slots, fourier_slots, cheb_jacobian, gradient):
    design = flatten_slots(cheb_slots)
    target = flatten_slots(fourier_slots)
    change, *_ = np.linalg.lstsq(design, target, rcond=None)
    slot_residual = target - design @ change
    relative = float(
        np.linalg.norm(slot_residual) / max(float(np.linalg.norm(target)), 1e-30))
    _residual, cheb_flat, n_nodes = N32.flatten_system(gradient, cheb_jacobian)
    transferred = cheb_flat @ change
    tangent = np.zeros((transferred.shape[1], n_nodes, 2), float)
    tangent[:, :, 0] = transferred[:n_nodes].T
    tangent[:, :, 1] = transferred[n_nodes:].T
    return tangent, relative, change.shape


def leftover_row(count, cheb_tangents, z, gradient, cheb_jacobian, settings):
    family = LocalFourierIncomingFamily(np.zeros((2, count)))
    _slots, fourier_slots = compatible_history_slots(family.metric(), z, 2 * count)
    tangent, relative, shape = transfer_assembled(
        cheb_tangents, fourier_slots, cheb_jacobian, gradient)
    leftover = N32.leftover_of(gradient, tangent, settings)
    unclipped = np.array([
        leftover['unclipped']['N_max'], leftover['unclipped']['beta_max']], float)
    clipped = np.array([
        leftover['clipped']['N_max'], leftover['clipped']['beta_max']], float)
    return {
        'count': count,
        'direction_count': int(tangent.shape[0]),
        'change_of_basis_shape': [int(shape[0]), int(shape[1])],
        'slot_transfer_relative_residual': _round(relative),
        'assembled_transfer_faithful': bool(relative <= TRANSFER_RELATIVE),
        'all_node_unclipped_leftover': _round(unclipped),
        'all_node_clipped_leftover': _round(clipped),
        'all_node_unclipped_step_norm': _round(leftover['unclipped']['step_norm']),
        'all_node_clipped_step_norm': _round(leftover['clipped']['step_norm']),
        'all_node_condition_number': _round(leftover['condition_number']),
        'unclipped_reaches_existence': bool(np.all(unclipped <= A.EXISTENCE_TOLERANCE)),
        'delta_is_none': leftover['delta'] is None,
    }


def compute():
    trial = json.loads(I4.OUTPUT.read_text())
    next5 = json.loads(N5.OUTPUT.read_text())
    for register in (trial, next5):
        for path, expected in {**register.get('source_hashes', {}),
                               **register.get('input_hashes', {})}.items():
            if digest(path) != expected:
                raise ValueError('assembled Fourier search input changed: '+path)
    if next5['evolve_authorized'] or next5['unclipped_leftover_reaches_existence_tolerance']:
        raise ValueError('n=16 assembled leftover already decides; this owner is the Fourier path')
    if trial['profile_identity'] in FORBIDDEN[:2]:
        raise ValueError('search must not start from a forbidden identity')
    ctx = I4.context()
    family16 = ctx['family']
    z = np.asarray(trial['z'], float)
    gradient = np.asarray(trial['action_gradient'], float)
    cheb_jacobian = np.asarray(trial['history_jacobian'], float)
    measured = np.asarray(trial['constraint_maxima'], float)
    _slots, cheb_tangents = compatible_history_slots(family16.metric(), z, 32)
    settings = LocalHistoryNewtonSettings()
    n8 = leftover_row(8, cheb_tangents, z, gradient, cheb_jacobian, settings)
    selected = n8
    n16 = None
    n8_leftover = np.asarray(n8['all_node_unclipped_leftover'], float)
    if (n8['assembled_transfer_faithful']
            and not n8['unclipped_reaches_existence']
            and n8['all_node_condition_number'] is not None
            and n8['all_node_condition_number'] < 1e12):
        n16 = leftover_row(16, cheb_tangents, z, gradient, cheb_jacobian, settings)
        n16_leftover = np.asarray(n16['all_node_unclipped_leftover'], float)
        if (n16['assembled_transfer_faithful']
                and n16['all_node_condition_number'] < 1e12
                and np.all(n16_leftover <= 0.5 * n8_leftover)):
            selected = n16
    leftover = np.asarray(selected['all_node_unclipped_leftover'], float)
    clipped = np.asarray(selected['all_node_clipped_leftover'], float)
    reaches = bool(np.all(leftover <= A.EXISTENCE_TOLERANCE))
    clipped_beats = bool(np.all(clipped < measured))
    clip_norm = selected['all_node_clipped_step_norm']
    evolve_authorized = bool(
        selected['assembled_transfer_faithful']
        and not selected['delta_is_none']
        and reaches
        and clipped_beats
        and clip_norm is not None and clip_norm <= 1.0 + 1e-12)
    cheb_leftover = np.asarray(next5['unclipped_linear_predicted_all_node_residual_maxima'], float)
    return {
        'schema': 'NSC-KS-FOURIER-ASSEMBLED-SEARCH-v1',
        'accountable_author': 'Douglas Ek',
        'status': (
            'OPEN: assembled Fourier leftover cannot reach EXISTENCE; no family evolved'
            if not evolve_authorized else
            'assembled Fourier leftover authorizes one retarded evolve'),
        'profile_identity': trial['profile_identity'],
        'representation': 'Fourier-on-I times owned axial plateau',
        'acts_on_assembled_residual': True,
        'not_geometry_leftover_alone': True,
        'no_invented_matter_columns': True,
        'slot_transfer_of_iterate4_assembled_jacobian': True,
        'selected_mode_count': selected['count'],
        'n8': {k: n8[k] for k in n8},
        'n16': None if n16 is None else {k: n16[k] for k in n16},
        'selected': {k: selected[k] for k in selected},
        'measured_constraint_maxima': measured.tolist(),
        'chebyshev_assembled_unclipped_leftover': _round(cheb_leftover),
        'all_node_unclipped_leftover': selected['all_node_unclipped_leftover'],
        'all_node_clipped_leftover': selected['all_node_clipped_leftover'],
        'all_node_condition_number': selected['all_node_condition_number'],
        'assembled_transfer_faithful': selected['assembled_transfer_faithful'],
        'unclipped_reaches_existence': reaches,
        'clipped_prediction_beats_measured': clipped_beats,
        'existence_tolerance': A.EXISTENCE_TOLERANCE,
        'forbidden_identities': list(FORBIDDEN),
        'evolve_authorized': evolve_authorized,
        'families_evolved': 0,
        'named_search_stall': (
            None if evolve_authorized else
            'finite_wu_geometry_image_cannot_reach_existence'),
        'fas_c_cause': (
            None if evolve_authorized else
            'insufficient representation of the declared class: assembled leftover stays residual-scale'),
        'scope': {
            'physical_local_gate': 'OPEN',
            'finite_image_is_not_class_identity': True,
            'no_class_rewrite': True,
            's5_not_next': True,
            'prediction_is_not_a_nonlinear_residual': True,
            'new_field_or_source_evolutions': 0,
            'physical_EXISTENCE_certificate': False,
            'physical_NONEXISTENCE_claimed': False,
            'metric_timestep': False,
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
            raise FileExistsError('assembled Fourier search exists; use --check')
        OUTPUT.write_text(json.dumps(result, sort_keys=True, indent=2, allow_nan=False)+'\n')
    elif result != json.loads(OUTPUT.read_text()):
        raise ValueError('assembled Fourier search replay differs')
    print(json.dumps({
        'status': result['status'],
        'selected_mode_count': result['selected_mode_count'],
        'all_node_unclipped_leftover': result['all_node_unclipped_leftover'],
        'all_node_condition_number': result['all_node_condition_number'],
        'assembled_transfer_faithful': result['assembled_transfer_faithful'],
        'evolve_authorized': result['evolve_authorized'],
        'families_evolved': result['families_evolved'],
        'named_search_stall': result['named_search_stall'],
    }, indent=2))
