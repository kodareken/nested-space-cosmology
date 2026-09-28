#!/usr/bin/env python3
"""Re-measure the iterate6 history with primal-block DOP853 control.

High Chebyshev coefficients are zero, so the radius polynomial is the
iterate6 history. 128 retarded directions are carried, and the step
controller ignores them. This is not a coefficient step and not a Newton
acceptance against the joint-norm iterate6 residual.
"""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
sys.path.insert(0, str(ROOT/'scripts'))
import derive_nsc_ks_coupled_newton_n16 as N16
import derive_nsc_ks_coupled_newton_n16_damped as D
import derive_nsc_ks_coupled_retained_control as R
import recursive_horizons.nsc_ks_energy_propagator as P
import recursive_horizons.nsc_ks_history_evaluator as H
from recursive_horizons.nsc_ks_primal_step_control import evolve_primal_controlled
from recursive_horizons.nsc_ks_profile_identity import profile_identity

import derive_nsc_ks_coupled_newton_n64_truncated_iterate as I64

PARENT = ROOT/'results/development/nsc-ks-coupled-newton-n32-truncated-iterate6.json'
DIRECTORY = ROOT/'results/development/artifacts/nsc-ks-coupled-newton-n64-primal-norm-zeropad'
OUTPUT = ROOT/'results/development/nsc-ks-coupled-newton-n64-primal-norm-zeropad.json'
OWNERS = (
    'scripts/derive_nsc_ks_n64_primal_norm_zeropad.py',
    'docs/nsc-ks-n64-primal-norm-zeropad.md',
    'src/recursive_horizons/nsc_ks_primal_step_control.py',
)


def digest(path):
    path = Path(path)
    return sha256((path if path.is_absolute() else ROOT/path).read_bytes()).hexdigest()


def context():
    parent = json.loads(PARENT.read_text())
    coefficients = np.asarray(parent['history']['coefficients'], float)
    if coefficients.shape != (2, 32):
        raise ValueError('iterate6 coefficient layout changed')
    padded = np.zeros((2, 64), float)
    padded[:, :32] = coefficients
    if np.max(np.abs(padded[:, 32:])) != 0.0:
        raise ValueError('zero pad must leave modes 32..63 unused')
    family = I64.LocalIncomingFamily64(padded)
    identity = profile_identity(family, include_normal_window=True)
    if identity in D.FORBIDDEN_IDENTITIES:
        raise ValueError('forbidden history must not be evolved')
    if identity == parent['profile_identity']:
        raise ValueError('64-coefficient identity unexpectedly matches the 32-coefficient profile')
    base = R.context()
    solve = family.collocation_nodes(N16.SOLVE_COUNT)
    verify = family.collocation_nodes(N16.VERIFY_COUNT)
    target = H.fixed_target_union(solve, verify)
    if target.tolist() != parent['z']:
        raise ValueError('zero pad must stay on the iterate6 nodes')
    return {
        **base,
        'family': family,
        'identity': identity,
        'solve': solve,
        'verify': verify,
        'target': target,
        'hashes': {
            **base['hashes'],
            **{path: digest(path) for path in OWNERS},
            str(PARENT.relative_to(ROOT)): digest(PARENT),
        },
        'proposal': {
            'constraint_maxima_with_verified_phase_values': parent['constraint_maxima'],
        },
        'selected': {'step_scale': 0.0, 'step_norm': 0.0},
        'parent_profile_identity': parent['profile_identity'],
        'parent_constraint_maxima': parent['constraint_maxima'],
    }


def annotate(result, ctx):
    measured = np.asarray(result['constraint_maxima'], float)
    parent = np.asarray(ctx['parent_constraint_maxima'], float)
    result['schema'] = 'NSC-KS-N64-PRIMAL-NORM-ZEROPAD-v1'
    result['status'] = (
        'OPEN: iterate6 radius remeasured with primal-block DOP853 control; uncertified physical gate'
        if result['all_retained_nonzero_angular_families_included']
        else 'OPEN: partial primal-block remeasure of the iterate6 radius')
    result['solver'] = 'primal-block DOP853 RMS on (A, D); tangent directions excluded'
    result['history_class'] = 'zero pad of iterate6 inside declared (w, U); no s^5'
    result['coefficient_step_evolved'] = False
    result['newton_step_accepted'] = False
    result['comparison_to_iterate6_is_a_solver_change'] = True
    result['parent_profile_identity'] = ctx['parent_profile_identity']
    result['iterate6_constraint_maxima'] = parent.tolist()
    result['measured_minus_iterate6'] = (measured - parent).tolist()
    result['residual_improved_on_all_components'] = False
    result['physical_EXISTENCE_certificate'] = False
    result['physical_NONEXISTENCE_certificate'] = False
    result['reproducer'] = 'python3 scripts/derive_nsc_ks_n64_primal_norm_zeropad.py --check'
    return result


def run(cpu_budget, max_new):
    P.evolve_ks_difference_envelope = evolve_primal_controlled
    N16.RETARDED = I64.DIRECTIONS
    N16.DIRECTORY = DIRECTORY
    N16.OUTPUT = OUTPUT
    R.FAMILY_CAP = 240.0
    R.CPU_CAP = 14400.0
    ctx = context()
    DIRECTORY.mkdir(parents=True, exist_ok=True)
    start = time.process_time()
    new = 0
    for key in sorted(ctx['archived']):
        path, _payload = N16.family_paths(key)
        if path.exists():
            continue
        if new >= max_new:
            continue
        remaining = cpu_budget - (time.process_time() - start)
        if remaining <= 1.0:
            break
        record = N16.new_family(ctx, key, cpu_limit=remaining - 1.0)
        new += record['new_operator_solves']
        print(json.dumps({
            'positive_family': record['positive_family'],
            'CPU_seconds': record['CPU_seconds'],
            'accepted_steps': record['operator']['diagnostics']['accepted_steps'],
            'matter_change_maxima': record['matter_change_maxima'],
        }), flush=True)
    result = annotate(N16.assemble(ctx), ctx)
    if result['all_retained_nonzero_angular_families_included']:
        OUTPUT.write_text(json.dumps(R.jsonable(result), indent=2, sort_keys=True)+'\n')
        progress = DIRECTORY/'progress.json'
        if progress.exists():
            progress.unlink()
    else:
        (DIRECTORY/'progress.json').write_text(json.dumps(R.jsonable(result), indent=2, sort_keys=True)+'\n')
    return result


def check():
    ctx = context()
    P.evolve_ks_difference_envelope = evolve_primal_controlled
    N16.RETARDED = I64.DIRECTIONS
    N16.DIRECTORY = DIRECTORY
    N16.OUTPUT = OUTPUT
    result = annotate(N16.assemble(ctx, replay=True), ctx)
    saved = json.loads(OUTPUT.read_text())
    if saved['profile_identity'] != result['profile_identity']:
        raise ValueError('zero-pad identity changed')
    if saved['constraint_maxima'] != result['constraint_maxima']:
        raise ValueError('zero-pad residual changed')
    if saved['newton_step_accepted'] or saved['coefficient_step_evolved']:
        raise ValueError('zero pad is not a coefficient step')
    if saved['physical_EXISTENCE_certificate'] or saved['physical_NONEXISTENCE_certificate']:
        raise ValueError('zero pad is not a gate certificate')
    return saved


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--run', action='store_true')
    mode.add_argument('--check', action='store_true')
    parser.add_argument('--cpu-budget', type=float, default=14400.0)
    parser.add_argument('--max-new', type=int, default=60)
    args = parser.parse_args()
    if args.check:
        saved = check()
        print(json.dumps({
            'constraint_maxima': saved['constraint_maxima'],
            'iterate6_constraint_maxima': saved['iterate6_constraint_maxima'],
            'measured_minus_iterate6': saved['measured_minus_iterate6'],
        }, indent=2))
    else:
        result = run(args.cpu_budget, args.max_new)
        print(json.dumps({
            'status': result['status'],
            'completed_positive_families': result['completed_positive_families'],
            'constraint_maxima': result['constraint_maxima'],
            'iterate6_constraint_maxima': result['iterate6_constraint_maxima'],
            'measured_minus_iterate6': result['measured_minus_iterate6'],
            'profile_identity': result['profile_identity'],
        }, indent=2))
