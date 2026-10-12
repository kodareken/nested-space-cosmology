#!/usr/bin/env python3
"""Dirac-free Fourier-on-I geometry leftover at iterate4 g. No matter columns."""
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
import derive_nsc_ks_n16_geometry_killing_current as K
import derive_nsc_ks_n16_leftover_anatomy as A
import derive_nsc_ks_n32_geometry_leftover as N32
from recursive_horizons.nsc_evolved_incoming_constraints import (
    compatible_history_slots, surface_geometry_response)
from recursive_horizons.nsc_ks_profile_identity import profile_identity
from recursive_horizons.nsc_local_history_newton import (
    HistoryCollocation, LocalHistoryNewtonSettings)
from recursive_horizons.nsc_local_incoming_fourier_family import (
    LocalFourierIncomingFamily, project_callable_on_fourier)

OUTPUT = ROOT/'results/development/nsc-ks-fourier-wu-leftover.json'
OWNERS = (
    'scripts/derive_nsc_ks_fourier_wu_leftover.py',
    'docs/nsc-ks-fourier-wu-leftover.md',
    'src/recursive_horizons/nsc_local_incoming_fourier_family.py',
)


def digest(path):
    path = Path(path)
    return sha256((path if path.is_absolute() else ROOT/path).read_bytes()).hexdigest()


def _round(value, digits=12):
    return N32._round(value, digits)


def orthonormalize_tangents(tangents):
    """One L2 orthonormalization of the slot columns. Same span."""
    flat = np.asarray(tangents, float).reshape(len(tangents), -1).T
    q, r = np.linalg.qr(flat, mode='reduced')
    rank = int(np.sum(np.abs(np.diag(r)) > 1e-14 * (abs(r[0, 0]) if r.size else 1.)))
    if rank == 0:
        raise ValueError('Fourier slot tangents have no numerical rank')
    return q[:, :rank].T.reshape(rank, *np.asarray(tangents).shape[1:])


def leftover_for(count, slots, z, gradient, coeff, settings, *, orthonormalize=False):
    family = LocalFourierIncomingFamily(np.zeros((2, count)))
    _, tangent = compatible_history_slots(family.metric(), z, 2 * count)
    if orthonormalize:
        tangent = orthonormalize_tangents(tangent)
    geometry = surface_geometry_response(slots, tangent, coeff)
    all_node = N32.leftover_of(gradient, geometry['action_gradient_tangent'], settings)
    residual, jacobian, n_nodes = N32.flatten_system(
        gradient, geometry['action_gradient_tangent'])
    return {
        'count': count,
        'direction_count': int(tangent.shape[0]),
        'orthonormalized': orthonormalize,
        'all_node': all_node,
        'residual': residual,
        'jacobian': jacobian,
        'n_nodes': n_nodes,
        'truncated': N32.truncated(
            residual, jacobian, all_node['singular'], n_nodes, settings),
    }


def summarize(row):
    leftover = np.array([
        row['all_node']['unclipped']['N_max'],
        row['all_node']['unclipped']['beta_max']], float)
    return {
        'mode_count': row['count'],
        'direction_count': row['direction_count'],
        'orthonormalized': row['orthonormalized'],
        'all_node_unclipped_leftover': _round(leftover),
        'all_node_clipped_leftover': _round([
            row['all_node']['clipped']['N_max'], row['all_node']['clipped']['beta_max']]),
        'all_node_unclipped_step_norm': _round(row['all_node']['unclipped']['step_norm']),
        'all_node_clipped_step_norm': _round(row['all_node']['clipped']['step_norm']),
        'all_node_condition_number': _round(row['all_node']['condition_number']),
        'all_node_rank_or_conditioning': row['all_node']['rank_or_conditioning'],
        'truncated_svd_all_node': row['truncated'],
        'unclipped_reaches_existence': bool(np.all(leftover <= A.EXISTENCE_TOLERANCE)),
    }


def compute():
    killing = json.loads(K.OUTPUT.read_text())
    trial = json.loads(I4.OUTPUT.read_text())
    n32 = json.loads(N32.OUTPUT.read_text())
    for path, expected in {**killing['source_hashes'], **killing['input_hashes'],
                           **trial['source_hashes'], **n32['source_hashes'],
                           **n32['input_hashes']}.items():
        if digest(path) != expected:
            raise ValueError('Fourier leftover input changed: '+path)
    ctx = I4.context()
    family16 = ctx['family']
    z = np.asarray(trial['z'], float)
    gradient = np.asarray(trial['action_gradient'], float)
    cheb_slots, _ = compatible_history_slots(family16.metric(), z, 32)
    settings = LocalHistoryNewtonSettings()
    n8 = leftover_for(8, cheb_slots, z, gradient, ctx['coeff'], settings)
    n8_ortho = leftover_for(
        8, cheb_slots, z, gradient, ctx['coeff'], settings, orthonormalize=True)
    selected = n8
    if ((n8['all_node']['condition_number'] > 1e16 or n8['all_node']['delta'] is None)
            and n8_ortho['all_node']['delta'] is not None
            and n8_ortho['all_node']['condition_number'] < n8['all_node']['condition_number']):
        selected = n8_ortho
    n8_leftover = np.array([
        selected['all_node']['unclipped']['N_max'],
        selected['all_node']['unclipped']['beta_max']], float)
    n8_condition = float(selected['all_node']['condition_number'])
    resolution_like = bool(
        np.all(n8_leftover > 10 * A.EXISTENCE_TOLERANCE) and n8_condition < 1e12)
    n16 = leftover_for(16, cheb_slots, z, gradient, ctx['coeff'], settings) if resolution_like else None
    if n16 is not None:
        n16_leftover = np.array([
            n16['all_node']['unclipped']['N_max'],
            n16['all_node']['unclipped']['beta_max']], float)
        if (np.all(n16_leftover <= 0.5 * n8_leftover)
                and float(n16['all_node']['condition_number']) < 1e12):
            selected = n16
    leftover = np.array([
        selected['all_node']['unclipped']['N_max'],
        selected['all_node']['unclipped']['beta_max']], float)
    measured = np.asarray(trial['constraint_maxima'], float)
    reaches = bool(np.all(leftover <= A.EXISTENCE_TOLERANCE))
    clipped_beats = bool(np.all(np.array([
        selected['all_node']['clipped']['N_max'],
        selected['all_node']['clipped']['beta_max']]) < measured))
    clip_norm = selected['all_node']['clipped']['step_norm']
    evolve_authorized = bool(
        selected['all_node']['delta'] is not None
        and reaches
        and clipped_beats
        and clip_norm is not None and clip_norm <= 1.0 + 1e-12)
    w, U = family16.functions
    projected = np.array([
        project_callable_on_fourier(lambda point, fn=fn: fn(point, 0),
                                    family16.center, family16.axial_inner,
                                    selected['count'])
        for fn in (w, U)
    ], float)
    fourier = LocalFourierIncomingFamily(projected)
    collocation = HistoryCollocation(16, layout='rectangular')
    solve_nodes = collocation.solve_nodes(family16)
    indices = np.searchsorted(z, solve_nodes)
    if not np.array_equal(z[indices], solve_nodes):
        raise ValueError('exact owned n=16 solve-node subset required')
    solve_basis = compatible_history_slots(
        LocalFourierIncomingFamily(np.zeros((2, selected['count']))).metric(),
        z, 2 * selected['count'])[1]
    if selected['orthonormalized']:
        solve_basis = orthonormalize_tangents(solve_basis)
    solve = N32.leftover_of(
        gradient[indices],
        surface_geometry_response(
            cheb_slots[indices], solve_basis[:, indices], ctx['coeff'])['action_gradient_tangent'],
        settings)
    return {
        'schema': 'NSC-KS-FOURIER-WU-LEFTOVER-v1',
        'accountable_author': 'Douglas Ek',
        'status': (
            'OPEN: Fourier-on-I geometry leftover cannot reach EXISTENCE; '
            'no family evolved'
            if not evolve_authorized else
            'Fourier-on-I geometry leftover authorizes one retarded evolve'),
        'current_profile_identity': profile_identity(family16),
        'representation': 'Fourier-on-I times owned axial plateau',
        'selected_mode_count': selected['count'],
        'selected_direction_count': selected['direction_count'],
        'orthonormalized_columns': selected['orthonormalized'],
        'raised_to_16': selected['count'] == 16,
        'resolution_like_at_8': resolution_like,
        'n8': summarize(n8),
        'n8_orthonormalized': summarize(n8_ortho),
        'n16': None if n16 is None else summarize(n16),
        'all_node_unclipped_leftover': _round(leftover),
        'all_node_clipped_leftover': _round([
            selected['all_node']['clipped']['N_max'],
            selected['all_node']['clipped']['beta_max']]),
        'all_node_unclipped_step_norm': _round(selected['all_node']['unclipped']['step_norm']),
        'all_node_clipped_step_norm': _round(selected['all_node']['clipped']['step_norm']),
        'all_node_condition_number': _round(selected['all_node']['condition_number']),
        'all_node_rank_or_conditioning': selected['all_node']['rank_or_conditioning'],
        'solve_unclipped_leftover': _round([
            solve['unclipped']['N_max'], solve['unclipped']['beta_max']]),
        'solve_unclipped_step_norm': _round(solve['unclipped']['step_norm']),
        'solve_condition_number': _round(solve['condition_number']),
        'truncated_svd_all_node': selected['truncated'],
        'projected_fourier_coefficients': [list(map(float, row)) for row in projected],
        'projected_radius_lower_bound': fourier.radius_lower_bound(),
        'unclipped_reaches_existence': reaches,
        'clipped_beats_measured': clipped_beats,
        'evolve_authorized': evolve_authorized,
        'existence_tolerance': A.EXISTENCE_TOLERANCE,
        'forbidden_identities': list(D.FORBIDDEN_IDENTITIES),
        'named_search_stall': (
            None if evolve_authorized else 'fourier_wu_leftover_cannot_reach_existence'),
        'scope': {
            'physical_local_gate': 'OPEN',
            'prediction_is_not_a_nonlinear_residual': True,
            'high_mode_matter_columns_invented': False,
            'new_field_or_source_evolutions': 0,
            'physical_EXISTENCE_certificate': False,
            'physical_NONEXISTENCE_claimed': False,
            'metric_timestep': False,
            'full_scale_n16_not_evolved': True,
            'n8_next_not_evolved': True,
            'chebyshev_family_not_widened': True,
        },
        'source_hashes': {p: digest(p) for p in OWNERS},
        'input_hashes': {
            str(K.OUTPUT.relative_to(ROOT)): digest(K.OUTPUT),
            str(I4.OUTPUT.relative_to(ROOT)): digest(I4.OUTPUT),
            str(A.OUTPUT.relative_to(ROOT)): digest(A.OUTPUT),
            str(N32.OUTPUT.relative_to(ROOT)): digest(N32.OUTPUT),
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
            raise FileExistsError('Fourier leftover exists; use --check')
        OUTPUT.write_text(json.dumps(result, sort_keys=True, indent=2, allow_nan=False)+'\n')
    elif result != json.loads(OUTPUT.read_text()):
        raise ValueError('Fourier leftover replay differs')
    print(json.dumps({
        'status': result['status'],
        'evolve_authorized': result['evolve_authorized'],
        'selected_mode_count': result['selected_mode_count'],
        'all_node_unclipped_leftover': result['all_node_unclipped_leftover'],
        'all_node_condition_number': result['all_node_condition_number'],
        'named_search_stall': result['named_search_stall'],
    }, indent=2))
