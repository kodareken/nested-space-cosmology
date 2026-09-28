#!/usr/bin/env python3
"""Feasibility of the history-dependent bounds on the value-only best.

No new field enclosure is run. The value-only ruler still moves by more
than the allocation, so a pilot at certificate precision cannot decide
EXISTENCE. ESCALATE is not INFEASIBLE and not NON-EXISTENCE.
"""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
sys.path.insert(0, str(ROOT/'scripts'))

OUTPUT = ROOT/'results/development/nsc-ks-gate-budget-feasibility.json'
REANCHOR = ROOT/'results/development/nsc-ks-gate-reanchor.json'
BUDGET = ROOT/'results/development/nsc-ks-gate-budget-independent.json'
OWNERS = (
    'scripts/derive_nsc_ks_gate_budget_feasibility.py',
    'docs/nsc-ks-gate-budget-feasibility.md',
)
COMPONENTS = (
    'finite_energy_field',
    'energy_interpolation',
    'changed_history_UV_tail',
    'phase_value',
    'between_node_on_I',
)


def digest(path):
    path = Path(path)
    return sha256((path if path.is_absolute() else ROOT/path).read_bytes()).hexdigest()


def compute():
    reanchor = json.loads(REANCHOR.read_text())
    budget = json.loads(BUDGET.read_text())
    merit = float(reanchor['campaign_best_merit'])
    if merit <= 1e-11:
        raise ValueError('a merit at the existence line would need a fresh feasibility pilot')
    rows = {}
    for name in COMPONENTS:
        rows[name] = {
            'method': 'not piloted on the value-only campaign best',
            'pilot_bound': None,
            'projected_cost': 'withheld until the value-only floor is below the allocation',
            'verdict': 'ESCALATE',
        }
    return {
        'schema': 'NSC-KS-GATE-BUDGET-FEASIBILITY-v1',
        'accountable_author': 'Douglas Ek',
        'status': (
            'OPEN: history-dependent bounds are ESCALATE; '
            'the value-only campaign best is above 1e-4'),
        'campaign_best_name': reanchor['campaign_best_name'],
        'campaign_best_merit': merit,
        'campaign_best_profile_identity': reanchor['campaign_best_profile_identity'],
        'components': rows,
        'infeasible_components': [],
        'physical_EXISTENCE_certificate': False,
        'physical_NONEXISTENCE_certificate': False,
        'named_gap': budget['named_gap'],
        'source_hashes': {path: digest(path) for path in OWNERS},
        'input_hashes': {
            str(REANCHOR.relative_to(ROOT)): digest(REANCHOR),
            str(BUDGET.relative_to(ROOT)): digest(BUDGET),
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
            raise FileExistsError('feasibility register exists; use --check')
        OUTPUT.write_text(json.dumps(result, sort_keys=True, indent=2)+'\n')
    elif result != json.loads(OUTPUT.read_text()):
        raise ValueError('feasibility replay differs')
    print(json.dumps({
        'status': result['status'],
        'campaign_best_merit': result['campaign_best_merit'],
        'named_gap': result['named_gap'],
    }, indent=2))
