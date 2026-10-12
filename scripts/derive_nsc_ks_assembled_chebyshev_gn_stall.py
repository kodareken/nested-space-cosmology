#!/usr/bin/env python3
"""Confirm the assembled Chebyshev n=16 leftover cannot decide EXISTENCE."""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
# Do not import iterate4/next5 modules: they pull the Dirac evolution stack.
I4_OUTPUT = ROOT/'results/development/nsc-ks-coupled-newton-n16-damped-iterate4.json'
N5_OUTPUT = ROOT/'results/development/nsc-ks-coupled-newton-n16-damped-next5.json'
OUTPUT = ROOT/'results/development/nsc-ks-assembled-chebyshev-gn-stall.json'
OWNERS = (
    'scripts/derive_nsc_ks_assembled_chebyshev_gn_stall.py',
    'docs/nsc-ks-assembled-chebyshev-gn-stall.md',
)
EXPECTED_IDENTITY = 'e15982e241bdc2244481491f38878b501b34ad3bbb208a54f40415d862482898'
FORBIDDEN = (
    '2d2588c3c9b1c83fe1a3f291314b8030300c4f14d3b12c24eadd2664eb20231b',
    '949a188181bd9ba0df43bc8ea00db2c89ce35a458017052b7af928059f819680',
    EXPECTED_IDENTITY,
)
EXISTENCE_TOLERANCE = 3e-11


def digest(path):
    path = Path(path)
    return sha256((path if path.is_absolute() else ROOT/path).read_bytes()).hexdigest()


def predicted_maxima(gradient, jacobian, step):
    prediction = gradient + np.tensordot(step, jacobian, axes=(0, 0))
    return np.max(abs(prediction), axis=0)


def unclipped_leftover(gradient, jacobian, direction, unclipped_norm):
    step = np.asarray(direction, float).reshape(-1) * float(unclipped_norm)
    return predicted_maxima(gradient, jacobian, step)


def bind_register(register):
    for path, expected in {**register.get('source_hashes', {}),
                           **register.get('input_hashes', {})}.items():
        if digest(path) != expected:
            raise ValueError('assembled Chebyshev GN stall input changed: '+path)


def compute():
    trial = json.loads(I4_OUTPUT.read_text())
    next5 = json.loads(N5_OUTPUT.read_text())
    # Iterate4 source_hashes name the retained family npz tree; bind those
    # bytes only through the iterate4 JSON digest already stored in next5.
    bind_register(next5)
    if trial['profile_identity'] != EXPECTED_IDENTITY:
        raise ValueError('stall must start from the accepted iterate4 identity')
    if trial['profile_identity'] != next5['current_profile_identity']:
        raise ValueError('next5 must sit on the same iterate4 history')
    if trial['profile_identity'] in FORBIDDEN[:2]:
        raise ValueError('search must not start from a forbidden identity')
    if next5['evolve_authorized'] or next5['unclipped_leftover_reaches_existence_tolerance']:
        raise ValueError('next5 leftover already decides; this owner records the leftover floor')
    gradient = np.asarray(trial['action_gradient'], float)
    jacobian = np.asarray(trial['history_jacobian'], float)
    measured = np.asarray(trial['constraint_maxima'], float)
    prediction = np.asarray(trial['linear_predicted_all_node_residual_maxima'], float)
    if jacobian.shape != (32,)+gradient.shape or gradient.ndim != 2 or gradient.shape[1] != 2:
        raise ValueError('iterate4 assembled Chebyshev n=16 Jacobian required')
    if not np.isfinite(gradient).all() or not np.isfinite(jacobian).all():
        raise ValueError('iterate4 residual/Jacobian must be finite')
    rect = [row for row in next5['proposals'] if row['layout'] == 'rectangular'][0]
    if rect['candidate_profile_identity'] in FORBIDDEN:
        raise ValueError('forbidden identity selected')
    if rect.get('unclipped_candidate_profile_identity') in FORBIDDEN[:2]:
        raise ValueError('forbidden unclipped identity selected')
    leftover = unclipped_leftover(
        gradient, jacobian, rect['direction'], rect['unclipped_step_norm'])
    recorded = np.asarray(next5['unclipped_linear_predicted_all_node_residual_maxima'], float)
    if not np.array_equal(leftover, recorded):
        raise ValueError('recomputed assembled leftover differs from next5')
    clipped = np.asarray(rect['predicted_all_node_residual_maxima'], float)
    reaches = bool(np.all(leftover <= EXISTENCE_TOLERANCE))
    clipped_beats = bool(np.all(clipped < measured))
    evolve_authorized = bool(reaches and clipped_beats and rect.get('selected_scale') is not None)
    if evolve_authorized:
        raise ValueError('leftover that reaches 3e-11 would require a fresh assembled evolve')
    return {
        'schema': 'NSC-KS-ASSEMBLED-CHEBYSHEV-GN-STALL-v1',
        'accountable_author': 'Douglas Ek',
        'status': (
            'OPEN: assembled Chebyshev n=16 leftover cannot reach EXISTENCE; '
            'frozen Jacobian leftover floor; no family evolved'),
        'profile_identity': trial['profile_identity'],
        'representation': 'Chebyshev n=16 assembled (baseline+geometry+matter+edge)',
        'acts_on_assembled_residual': True,
        'not_geometry_leftover_alone': True,
        'no_invented_matter_columns': True,
        'no_slot_transfer': True,
        'full_retarded_history_jacobian_consumed': True,
        'full_minus_geometric_jacobian_max': rect['full_minus_geometric_jacobian_max'],
        'measured_constraint_maxima': measured.tolist(),
        'iterate4_linear_predicted_all_node_residual_maxima': prediction.tolist(),
        'iterate4_prediction_minus_measured': (prediction - measured).tolist(),
        'clipped_linear_predicted_all_node_residual_maxima': clipped.tolist(),
        'clipped_prediction_beats_measured': clipped_beats,
        'clipped_candidate_profile_identity': rect['candidate_profile_identity'],
        'clipped_step_norm': rect['clipped_step_norm'],
        'unclipped_step_norm': rect['unclipped_step_norm'],
        'assembled_unclipped_leftover': leftover.tolist(),
        'unclipped_leftover_reaches_existence_tolerance': reaches,
        'existence_tolerance': EXISTENCE_TOLERANCE,
        'forbidden_identities': list(FORBIDDEN),
        'forbidden_identity_selected': False,
        'prediction_is_not_a_nonlinear_residual': True,
        'evolve_authorized': evolve_authorized,
        'families_evolved': 0,
        'named_search_stall': 'finite_wu_geometry_image_cannot_reach_existence',
        'fas_c_cause': (
            'insufficient representation / leftover floor of the declared class'),
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
            str(I4_OUTPUT.relative_to(ROOT)): digest(I4_OUTPUT),
            str(N5_OUTPUT.relative_to(ROOT)): digest(N5_OUTPUT),
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
            raise FileExistsError('assembled Chebyshev GN stall exists; use --check')
        OUTPUT.write_text(json.dumps(result, sort_keys=True, indent=2, allow_nan=False)+'\n')
    elif result != json.loads(OUTPUT.read_text()):
        raise ValueError('assembled Chebyshev GN stall replay differs')
    print(json.dumps({
        'status': result['status'],
        'assembled_unclipped_leftover': result['assembled_unclipped_leftover'],
        'unclipped_leftover_reaches_existence_tolerance': result[
            'unclipped_leftover_reaches_existence_tolerance'],
        'evolve_authorized': result['evolve_authorized'],
        'families_evolved': result['families_evolved'],
        'named_search_stall': result['named_search_stall'],
    }, indent=2))
