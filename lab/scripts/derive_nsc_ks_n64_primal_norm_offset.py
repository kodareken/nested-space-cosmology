#!/usr/bin/env python3
"""Explain the n=64 N offset from saved family solves, then test its removal.

The two rank-64 evolutions share one fewer DOP853 step than iterate6 on the
same twenty families. Their step-independent N piece is that common-mode
matter shift. Beta stays with the Jacobian. This is the joint RMS norm
letting tangent slots dilute the primal tolerance, not a missing constant
mode and not a class obstruction. A zero coefficient step is not evolved here.
"""
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
import derive_nsc_ks_coupled_newton_n16 as N16
import derive_nsc_ks_coupled_retained_control as R
import recursive_horizons.nsc_ks_energy_propagator as P
import recursive_horizons.nsc_ks_history_evaluator as H
from recursive_horizons.nsc_ks_primal_step_control import evolve_primal_controlled
from recursive_horizons.nsc_ks_profile_identity import profile_identity
from recursive_horizons.nsc_local_incoming_family import LocalIncomingFamily

import derive_nsc_ks_coupled_newton_n64_truncated_iterate as I64

ARTIFACTS = ROOT/'results/development/artifacts'
CAMPAIGNS = {
    'iterate6': 'nsc-ks-coupled-newton-n32-truncated-iterate6',
    'iterate8': 'nsc-ks-coupled-newton-n64-truncated-iterate8',
    'iterate9': 'nsc-ks-coupled-newton-n64-truncated-iterate9',
}
PARENT = ROOT/'results/development/nsc-ks-coupled-newton-n32-truncated-iterate6.json'
OUTPUT = ROOT/'results/development/nsc-ks-n64-primal-norm-offset.json'
OWNERS = (
    'scripts/derive_nsc_ks_n64_primal_norm_offset.py',
    'docs/nsc-ks-primal-step-control.md',
    'src/recursive_horizons/nsc_ks_primal_step_control.py',
)
PROBE_KEYS = ((32, -1), (1, 1))
AGREEMENT = 1e-12


def digest(path):
    return sha256((ROOT/path).read_bytes()).hexdigest()


def _family_rows():
    folder = {name: ARTIFACTS/stem for name, stem in CAMPAIGNS.items()}
    stems = sorted(path.name[:-4] for path in folder['iterate6'].glob('family-*.npz'))
    if len(stems) != 60:
        raise ValueError('iterate6 family set is not the retained sixty')
    rows = []
    for stem in stems:
        matter, steps, ell = {}, {}, None
        for name, directory in folder.items():
            with np.load(directory/f'{stem}.npz') as payload:
                matter[name] = np.asarray(payload['paired_matter_change'], float)
            record = json.loads((directory/f'{stem}.json').read_text())
            steps[name] = int(record['operator']['diagnostics']['accepted_steps'])
            if ell is None:
                ell = float(next(iter(record['family_records'].values()))['angular_eigenvalue'])
        if matter['iterate8'].shape != matter['iterate6'].shape:
            raise ValueError('family matter nodes changed')
        delta8 = matter['iterate8'] - matter['iterate6']
        delta9 = matter['iterate9'] - matter['iterate6']
        common = (32.0 * delta9 - delta8) / 31.0
        linear = delta8 - common
        same = steps['iterate6'] == steps['iterate8'] == steps['iterate9']
        rows.append({
            'family': stem,
            'angular_eigenvalue': ell,
            'accepted_steps': [steps['iterate6'], steps['iterate8'], steps['iterate9']],
            'same_step_count': same,
            'common_mode_N_mean': float(common[:, 0].mean()),
            'linear_N_mean': float(linear[:, 0].mean()),
            'common_mode_beta_mean': float(common[:, 1].mean()),
            'common_mode_N_flatness': float(np.max(np.abs(common[:, 0] - common[:, 0].mean()))),
        })
    return rows


def decompose():
    rows = _family_rows()
    same = [row for row in rows if row['same_step_count']]
    changed = [row for row in rows if not row['same_step_count']]
    if any(row['accepted_steps'][1] != row['accepted_steps'][2] for row in rows):
        raise ValueError('iterate8 and iterate9 must share a step count on every family')
    if any(row['accepted_steps'][0] - row['accepted_steps'][1] not in (0, 1) for row in rows):
        raise ValueError('the n=64 runs drop at most one accepted step')
    common_n = float(sum(row['common_mode_N_mean'] for row in rows))
    changed_n = float(sum(row['common_mode_N_mean'] for row in changed))
    same_n = float(sum(row['common_mode_N_mean'] for row in same))
    return {
        'families': rows,
        'common_mode_N_sum': common_n,
        'common_mode_beta_sum': float(sum(row['common_mode_beta_mean'] for row in rows)),
        'linear_N_sum': float(sum(row['linear_N_mean'] for row in rows)),
        'same_step_family_count': len(same),
        'same_step_common_mode_N_sum': same_n,
        'step_drop_family_count': len(changed),
        'step_drop_common_mode_N_sum': changed_n,
        'step_drop_fraction_of_common_mode_N': changed_n / common_n,
    }


def _context(base, coefficients):
    parent = json.loads(PARENT.read_text())
    if coefficients.shape == (2, 32):
        family = LocalIncomingFamily(coefficients)
        directions = 64
    elif coefficients.shape == (2, 64):
        family = I64.LocalIncomingFamily64(coefficients)
        directions = 128
    else:
        raise ValueError('probe history must be 32 or 64 coefficients per function')
    identity = profile_identity(family, include_normal_window=True)
    solve = family.collocation_nodes(N16.SOLVE_COUNT)
    verify = family.collocation_nodes(N16.VERIFY_COUNT)
    target = H.fixed_target_union(solve, verify)
    if target.tolist() != parent['z']:
        raise ValueError('probe must stay on the iterate6 nodes')
    ctx = dict(base)
    ctx.update({
        'family': family, 'identity': identity, 'solve': solve,
        'verify': verify, 'target': target,
    })
    return ctx, directions


def _probe_one(base, coefficients, key):
    ctx, directions = _context(base, coefficients)
    N16.RETARDED = directions
    P.evolve_ks_difference_envelope = evolve_primal_controlled
    arrays, meta, _signed, _operator = N16.evaluate_family(ctx, key, saved=None, allow_new=True)
    matter = np.asarray(arrays['paired_matter_change'], float)
    return {
        'positive_family': [int(key[0]), int(key[1])],
        'directions': directions,
        'accepted_steps': int(meta['operator']['diagnostics']['accepted_steps']),
        'function_evaluations': int(meta['operator']['diagnostics']['function_evaluations']),
        'matter_N_mean': float(matter[:, 0].mean()),
        'matter_beta_mean': float(matter[:, 1].mean()),
        'matter_change_maxima': np.max(np.abs(matter), axis=0).tolist(),
        'profile_identity': ctx['identity'],
    }


def probe(keys):
    parent = json.loads(PARENT.read_text())
    coefficients = np.asarray(parent['history']['coefficients'], float)
    if coefficients.shape != (2, 32):
        raise ValueError('iterate6 coefficient layout changed')
    padded = np.zeros((2, 64), float)
    padded[:, :32] = coefficients
    base = R.context()
    runs = []
    for key in keys:
        narrow = _probe_one(base, coefficients, key)
        wide = _probe_one(base, padded, key)
        runs.append({'narrow': narrow, 'wide': wide})
        print(json.dumps({
            'family': [int(key[0]), int(key[1])],
            'steps': [narrow['accepted_steps'], wide['accepted_steps']],
            'N_mean': [narrow['matter_N_mean'], wide['matter_N_mean']],
        }), flush=True)
    return runs


def build(runs=None):
    record = {
        'schema': 'NSC-KS-N64-PRIMAL-NORM-OFFSET-v1',
        'accountable_author': 'Douglas Ek',
        'status': 'OPEN: N offset is a joint-norm step drop; primal-block control is measured on two families',
        'named_search_stall': 'n64_joint_rms_norm_drops_one_primal_step',
        'previous_stall': 'n64_truncated_N_offset_independent_of_step',
        'classification': 'bookkeeping: joint DOP853 RMS norm dilutes the primal tolerance by the tangent block',
        'not_a_missing_constant_mode': True,
        'not_an_edge_or_baseline_piece': True,
        'not_class_nonexistence': True,
        'physical_EXISTENCE_certificate': False,
        'physical_NONEXISTENCE_certificate': False,
        'gate': 'OPEN',
        'existence_tolerance': 3e-11,
        'best_measured_profile_identity': json.loads(PARENT.read_text())['profile_identity'],
        'best_measured_constraint_maxima': json.loads(PARENT.read_text())['constraint_maxima'],
        'decomposition': decompose(),
        'probe': runs,
        'probe_agreement_tolerance': AGREEMENT,
        'coefficient_step_evolved': False,
        'source_hashes': {path: digest(path) for path in OWNERS},
    }
    if runs is not None:
        for pair in runs:
            gap = abs(pair['narrow']['matter_N_mean'] - pair['wide']['matter_N_mean'])
            if pair['narrow']['accepted_steps'] != pair['wide']['accepted_steps'] or gap > AGREEMENT:
                record['status'] = 'OPEN: primal-block control did not remove the family N split'
    return record


def check(record):
    fresh = decompose()
    stored = record['decomposition']
    for key in (
            'common_mode_N_sum', 'common_mode_beta_sum', 'linear_N_sum',
            'same_step_common_mode_N_sum', 'step_drop_common_mode_N_sum'):
        if abs(fresh[key] - stored[key]) > 1e-15 * max(1.0, abs(fresh[key])):
            raise ValueError('stored step decomposition does not replay')
    if fresh['step_drop_family_count'] != 20 or fresh['same_step_family_count'] != 40:
        raise ValueError('step-drop family count changed')
    if abs(fresh['step_drop_fraction_of_common_mode_N']) < 0.999:
        raise ValueError('the common-mode N offset is not carried by the step-drop families')
    if abs(fresh['same_step_common_mode_N_sum']) > 1e-6:
        raise ValueError('same-step families still carry a visible common-mode N')
    if abs(fresh['common_mode_beta_sum']) > 1e-10:
        raise ValueError('the common mode is not confined to N')
    if record['source_hashes'] != {path: digest(path) for path in OWNERS}:
        raise ValueError('primal-norm owner hashes changed')
    if not record['probe']:
        raise ValueError('probe measurements are part of this register')
    for pair in record['probe']:
        gap = abs(pair['narrow']['matter_N_mean'] - pair['wide']['matter_N_mean'])
        if pair['narrow']['accepted_steps'] != pair['wide']['accepted_steps'] or gap > AGREEMENT:
            raise ValueError('primal-controlled 64- and 128-direction matters disagree')
        if pair['narrow']['directions'] != 64 or pair['wide']['directions'] != 128:
            raise ValueError('probe must compare 64 and 128 directions')
    if record['physical_EXISTENCE_certificate'] or record['physical_NONEXISTENCE_certificate']:
        raise ValueError('this diagnosis is not a gate certificate')
    return record


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--record', action='store_true')
    mode.add_argument('--probe', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    if args.check:
        check(json.loads(OUTPUT.read_text()))
        print('PASS')
    elif args.record:
        if OUTPUT.exists():
            raise SystemExit('offset register already exists')
        OUTPUT.write_text(json.dumps(build(), indent=2, sort_keys=True)+'\n')
        print(json.dumps({
            'common_mode_N_sum': build()['decomposition']['common_mode_N_sum'],
            'step_drop_common_mode_N_sum': build()['decomposition']['step_drop_common_mode_N_sum'],
        }, indent=2))
    else:
        runs = probe(PROBE_KEYS)
        OUTPUT.write_text(json.dumps(build(runs), indent=2, sort_keys=True)+'\n')
        check(json.loads(OUTPUT.read_text()))
        print('PASS')
