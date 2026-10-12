#!/usr/bin/env python3
"""History-independent budget against the declared allocation.

A component is REUSE only when an existing certificate is cited by hash.
An indicator is not a bound. A missing bound is an escalation, not a zero.
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

OUTPUT = ROOT/'results/development/nsc-ks-gate-budget-independent.json'
FLOOR = ROOT/'results/development/nsc-ks-gate-numerical-floor.json'
V5 = ROOT/'results/development/nsc-incoming-source-update-v5.json'
BUFFER = ROOT/'results/development/nsc-incoming-boundary-buffer.json'
CF4 = ROOT/'results/development/nsc-local-cf4-response.json'
OWNERS = (
    'scripts/derive_nsc_ks_gate_budget_independent.py',
    'docs/nsc-ks-gate-budget-independent.md',
)
ALLOCATION = {
    'field_space_time': 5e-12,
    'changed_history_UV_tail': 5e-12,
    'baseline_low_subgap': 3e-12,
    'upstream': 2e-12,
    'energy_interpolation': 1e-12,
    'covered_regions': 1e-12,
    'phase_value': 1e-12,
    'between_node': 1e-12,
    'arithmetic': 1e-12,
}


def digest(path):
    path = Path(path)
    return sha256((path if path.is_absolute() else ROOT/path).read_bytes()).hexdigest()


def compute():
    floor = json.loads(FLOOR.read_text())
    movement = float(floor['largest_one_knob_movement'])
    if movement <= 1e-10:
        raise ValueError('the recorded floor no longer forces an escalation')
    cited = {
        str(path.relative_to(ROOT)): digest(path)
        for path in (V5, BUFFER, CF4, FLOOR)
    }
    components = {}
    for name, allocation in ALLOCATION.items():
        components[name] = {
            'allocation': allocation,
            'bound': None,
            'verdict': 'ESCALATE',
            'reason': (
                'value-only one-knob movement %.3e exceeds this allocation; '
                'no changed-history enclosure is claimed' % movement),
        }
    components['covered_regions']['cited_register'] = str(V5.relative_to(ROOT))
    components['covered_regions']['cited_status'] = json.loads(V5.read_text())['status']
    components['field_space_time']['cited_registers'] = [
        str(BUFFER.relative_to(ROOT)), str(CF4.relative_to(ROOT))]
    return {
        'schema': 'NSC-KS-GATE-BUDGET-INDEPENDENT-v1',
        'accountable_author': 'Douglas Ek',
        'status': (
            'OPEN: every history-independent component is ESCALATE; '
            'the value-only floor exceeds the allocations'),
        'allocation': ALLOCATION,
        'components': components,
        'largest_one_knob_movement': movement,
        'bounds_spent': False,
        'physical_EXISTENCE_certificate': False,
        'physical_NONEXISTENCE_certificate': False,
        'named_gap': 'value_only_indicator_floor_above_1e-10',
        'source_hashes': {path: digest(path) for path in OWNERS},
        'input_hashes': cited,
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
            raise FileExistsError('budget register exists; use --check')
        OUTPUT.write_text(json.dumps(result, sort_keys=True, indent=2)+'\n')
    elif result != json.loads(OUTPUT.read_text()):
        raise ValueError('budget replay differs')
    print(json.dumps({
        'status': result['status'],
        'verdicts': {name: row['verdict'] for name, row in result['components'].items()},
        'named_gap': result['named_gap'],
    }, indent=2))
