#!/usr/bin/env python3
"""One-knob value-only ladder for the iterate6 history.

Each row is an indicator: the change in max |E_N| and max |E_beta| when one
solver setting moves. It is not an error bound. The base row is the measured
47-node value-only residual. Later knobs are evolved only when no other
derive_nsc_ --run is live.
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
import derive_nsc_ks_gate_value as V

OUTPUT = ROOT/'results/development/nsc-ks-gate-numerical-floor.json'
BASE = ROOT/'results/development/nsc-ks-gate-value-iterate6-47.json'
OWNERS = (
    'scripts/derive_nsc_ks_gate_numerical_floor.py',
    'docs/nsc-ks-gate-numerical-floor.md',
    'scripts/derive_nsc_ks_gate_value.py',
)
# One change from the base solver. Order is the ladder order.
LADDER = (
    ('rtol-2e-12-n47', 47, {'rtol': 2e-12}),
    ('rtol-2e-13-n47', 47, {'rtol': 2e-13}),
    ('max-step-1e-3-n47', 47, {'max_step': 0.001}),
    ('max-step-5e-4-n47', 47, {'max_step': 0.0005}),
    ('grid-128-n47', 47, {'grid_nodes': 128}),
    ('grid-256-n47', 47, {'grid_nodes': 256}),
    ('degree-48-n47', 47, {'degree': 48}),
    ('phase-384-n47', 47, {'phase_nodes': 384}),
    ('base-n129', 129, {}),
)


def digest(path):
    path = Path(path)
    return sha256((path if path.is_absolute() else ROOT/path).read_bytes()).hexdigest()


def base_row():
    measured = json.loads(BASE.read_text())
    maxima = measured['constraint_maxima']
    return {
        'label': 'base-n47',
        'knob': 'base',
        'nodes': 47,
        'solver': V.solver_settings(None),
        'constraint_maxima': maxima,
        'merit': measured['merit'],
        'delta_max_abs_from_base': [0.0, 0.0],
        'wall_seconds': measured.get('wall_seconds'),
        'indicator_not_a_bound': True,
        'families_evolved_for_this_row': 0,
    }


def completed_row(label, nodes, overrides):
    path = V.output_paths(label)[1]
    if not path.exists():
        return None
    measured = json.loads(path.read_text())
    base = base_row()['constraint_maxima']
    maxima = measured['constraint_maxima']
    return {
        'label': label,
        'knob': label,
        'nodes': nodes,
        'solver': measured.get('solver', V.solver_settings(overrides)),
        'constraint_maxima': maxima,
        'merit': measured['merit'],
        'delta_max_abs_from_base': [abs(maxima[0] - base[0]), abs(maxima[1] - base[1])],
        'wall_seconds': measured.get('wall_seconds'),
        'indicator_not_a_bound': True,
        'families_evolved_for_this_row': 0,
    }


def compute():
    rows = [base_row()]
    pending = []
    for label, nodes, overrides in LADDER:
        row = completed_row(label, nodes, overrides)
        if row is None:
            pending.append(label)
        else:
            rows.append(row)
    within = [
        row['label'] for row in rows
        if row['label'] != 'base-n47'
        and row['delta_max_abs_from_base'][0] <= 1e-10
        and row['delta_max_abs_from_base'][1] <= 1e-10]
    movements = [max(row['delta_max_abs_from_base']) for row in rows if row['label'] != 'base-n47']
    return {
        'schema': 'NSC-KS-GATE-NUMERICAL-FLOOR-v1',
        'accountable_author': 'Douglas Ek',
        'status': (
            'OPEN: no one-knob setting keeps both components inside 1e-10; '
            'loop and certificate numerics are not declared'),
        'indicator_not_a_bound': True,
        'drift_subtracted_as_physics': False,
        'loop_numerics': None,
        'certificate_numerics': None,
        'knobs_within_1e-10_on_both_components': within,
        'largest_one_knob_movement': None if not movements else max(movements),
        'fas_a_required_before_search': True,
        'rows': rows,
        'pending_labels': pending,
        'physical_EXISTENCE_certificate': False,
        'physical_NONEXISTENCE_certificate': False,
        'named_gap': pending[0] if pending else 'value_only_indicator_floor_above_1e-10',
        'source_hashes': {path: digest(path) for path in OWNERS},
        'input_hashes': {str(BASE.relative_to(ROOT)): digest(BASE)},
    }


def run_next(workers):
    result = compute()
    if not result['pending_labels']:
        return result
    label = result['pending_labels'][0]
    nodes, overrides = next((n, o) for name, n, o in LADDER if name == label)
    V.run(nodes, workers, 14400.0, label, overrides or None)
    return compute()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--record', action='store_true')
    mode.add_argument('--check', action='store_true')
    mode.add_argument('--next', action='store_true')
    parser.add_argument('--workers', type=int, default=9)
    args = parser.parse_args()
    if args.next:
        result = run_next(args.workers)
        OUTPUT.write_text(json.dumps(result, sort_keys=True, indent=2, allow_nan=False)+'\n')
    else:
        result = compute()
        if args.record:
            if OUTPUT.exists():
                raise FileExistsError('numerical-floor register exists; use --check or --next')
            OUTPUT.write_text(json.dumps(result, sort_keys=True, indent=2, allow_nan=False)+'\n')
        elif result != json.loads(OUTPUT.read_text()):
            raise ValueError('numerical-floor replay differs')
    print(json.dumps({
        'status': result['status'],
        'named_gap': result['named_gap'],
        'pending_labels': result['pending_labels'],
        'base_merit': result['rows'][0]['merit'],
    }, indent=2))
