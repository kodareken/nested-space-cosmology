#!/usr/bin/env python3
"""One ledger and one step record for the declared (w, U) search.

The merit is fixed before any step uses it:

    m(g) = max over nodes of max(|E_N|, |E_beta|).

Rows are append-only. A prediction is stored as a prediction.
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
from recursive_horizons.nsc_ks_profile_identity import profile_identity
from recursive_horizons.nsc_ks_value_evaluator import MERIT_DEFINITION, merit
from recursive_horizons.nsc_local_incoming_family import LocalIncomingFamily

LEDGER = ROOT/'results/development/nsc-ks-gate-ledger.json'
FORBIDDEN = (
    '2d2588c3c9b1c83fe1a3f291314b8030300c4f14d3b12c24eadd2664eb20231b',
    '949a188181bd9ba0df43bc8ea00db2c89ce35a458017052b7af928059f819680',
)
OWNERS = (
    'scripts/derive_nsc_ks_gate_step.py',
    'docs/nsc-ks-gate-step.md',
    'src/recursive_horizons/nsc_ks_value_evaluator.py',
)


def digest(path):
    path = Path(path)
    return sha256((path if path.is_absolute() else ROOT/path).read_bytes()).hexdigest()


def empty_ledger():
    return {
        'schema': 'NSC-KS-GATE-LEDGER-v1',
        'accountable_author': 'Douglas Ek',
        'merit_definition': MERIT_DEFINITION,
        'merit_declared_before_use': True,
        'rows': [],
        'source_hashes': {p: digest(p) for p in OWNERS},
    }


def load_ledger(path=LEDGER):
    if not Path(path).exists():
        raise FileNotFoundError('gate ledger is not declared')
    ledger = json.loads(Path(path).read_text())
    if ledger.get('merit_definition') != MERIT_DEFINITION or not ledger.get('merit_declared_before_use'):
        raise ValueError('gate ledger merit was not declared before use')
    return ledger


def declare_ledger(path=LEDGER):
    path = Path(path)
    if path.exists():
        return load_ledger(path)
    ledger = empty_ledger()
    path.write_text(json.dumps(ledger, sort_keys=True, indent=2)+'\n')
    return ledger


def append_row(row, path=LEDGER):
    ledger = load_ledger(path)
    required = (
        'parent', 'rule', 'parameters', 'step_size',
        'predicted_merit', 'measured_merit', 'accepted')
    missing = [name for name in required if name not in row]
    if missing:
        raise ValueError('ledger row is missing ' + missing[0])
    if row['accepted'] not in (True, False):
        raise ValueError('ledger acceptance must be true or false')
    identity = row.get('candidate_profile_identity')
    if identity in FORBIDDEN:
        raise ValueError('forbidden history cannot be appended')
    previous = json.dumps(ledger['rows'], sort_keys=True)
    ledger['rows'].append(row)
    if not json.dumps(json.loads(previous), sort_keys=True) == previous:
        raise ValueError('ledger history was rewritten')
    Path(path).write_text(json.dumps(ledger, sort_keys=True, indent=2)+'\n')
    reloaded = json.loads(Path(path).read_text())
    if reloaded['rows'][:-1] != ledger['rows'][:-1]:
        raise ValueError('append changed an existing ledger row')
    return reloaded


def proposal(coefficients, step, parent_identity, parent_merit):
    """A step inside the declared class. It does not evolve and is not a residual."""
    coefficients = np.asarray(coefficients, float)
    step = np.asarray(step, float)
    if coefficients.shape[0] != 2 or coefficients.shape[1] not in (32, 64):
        raise ValueError('step owner accepts 32 or 64 coefficients per function')
    if step.shape != coefficients.shape:
        raise ValueError('step and coefficient shapes differ')
    family = LocalIncomingFamily(coefficients + step)
    identity = profile_identity(family, include_normal_window=True)
    if identity in FORBIDDEN or identity == parent_identity:
        raise ValueError('step candidate is forbidden or unchanged')
    if family.radius_lower_bound() <= 0:
        raise ValueError('step candidate must keep a positive radius')
    return {
        'schema': 'NSC-KS-GATE-STEP-PROPOSAL-v1',
        'merit_definition': MERIT_DEFINITION,
        'parent_profile_identity': parent_identity,
        'parent_merit': None if parent_merit is None else float(parent_merit),
        'candidate_profile_identity': identity,
        'candidate_coefficients': family.description()['coefficients'],
        'step_size': float(np.linalg.norm(step)),
        'radius_lower_bound': family.radius_lower_bound(),
        'predicted_merit': None,
        'prediction_is_not_a_residual': True,
        'evolve_authorized': False,
    }


def measured_merit(gradient):
    return merit(gradient)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--declare', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    if args.declare:
        ledger = declare_ledger()
    else:
        ledger = load_ledger()
        fresh = empty_ledger()
        fresh['rows'] = ledger['rows']
        fresh['source_hashes'] = ledger['source_hashes']
        if {k: ledger[k] for k in ('merit_definition', 'merit_declared_before_use')} != {
                'merit_definition': MERIT_DEFINITION, 'merit_declared_before_use': True}:
            raise ValueError('gate ledger replay differs')
    print(json.dumps({
        'merit_definition': ledger['merit_definition'],
        'rows': len(ledger['rows']),
    }, indent=2))
