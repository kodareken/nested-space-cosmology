#!/usr/bin/env python3
"""Demonstrate the cubic leading addition on saved finite-cutoff value nodes.

Coefficient transport is the existing shared channel basis at production
resolution. Missing nodes stay empty. The center value is not broadcast, a
tiny cell count is not substituted, and C_M blocks certification.
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

from flint import arb

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parent
sys.path.insert(0, str(ROOT / 'src'))
sys.path.insert(0, str(ROOT / 'scripts'))
from recursive_horizons.nsc_ks_cubic_channel_basis import (
    CubicChannelBudgetExceeded, enclose_cubic_channel_basis, unit_cubic_geometry)
from recursive_horizons.nsc_ks_cubic_residual import (
    VALUE_LABEL, apply_cubic_leading, certification_blockers, certify,
    load_saved_finite_cutoff_value)
from recursive_horizons.nsc_ks_finite_history_uv_coefficients import ARCHIVE_RHO_UP
from recursive_horizons.nsc_ks_massive_cubic_uv import require_original_history

CPU_BUDGET = 180.0
CELLS = 1024
ORDER = 8
BORROWED = (
    'lab/src/recursive_horizons/nsc_ks_cubic_inventory.py',
    'lab/tests/test_nsc_ks_cubic_inventory.py',
    'lab/docs/nsc-ks-cubic-inventory.md',
)
CHANNEL_NAMES = ('N2', 'N4', 'NM', 'J2', 'J4', 'JM')


def _borrowed_hashes():
    hashes = {}
    for relative in BORROWED:
        path = REPO / relative
        if not path.is_file():
            raise ValueError('borrowed cubic inventory file missing')
        hashes[relative] = sha256(path.read_bytes()).hexdigest()
    return hashes


def _stack(pairs, identity, rho_up):
    bases = {}
    for sign in (1, -1):
        channels = {name: [] for name in CHANNEL_NAMES}
        domains = []
        for pair in pairs:
            basis = pair[sign]
            if (basis.get('sign') != sign or basis.get('massless_A') is not True
                    or basis.get('profile_identity') != identity
                    or basis.get('C_M') is not None
                    or basis.get('cells') != CELLS or basis.get('order') != ORDER
                    or basis.get('accepted_cells') != CELLS
                    or basis.get('entire_incoming_interval')
                    or basis.get('z_box_is_full_incoming_interval')):
                raise ValueError('tiny-box pilot is not a target-node enclosure')
            for name in CHANNEL_NAMES:
                channels[name].append(basis[name])
            domains.append(basis['z_domain'])
        bases[sign] = {
            'sign': sign, 'massless_A': True, 'profile_identity': identity,
            'C_M': None, 'higher_uv_remainder_C_M': None,
            'rho_up': arb(float(rho_up)), 'z_domain': domains, 'channels': channels,
        }
    return bases


def _node_record(before, report, position):
    after = report['residual'][position]
    delta = report['approximate_leading'][position]
    enclosure = report['leading_enclosures'][position]
    return {
        'node_index': enclosure['node_index'],
        'z_hex': enclosure['z_hex'],
        'residual_before_hex': [float(before[0]).hex(), float(before[1]).hex()],
        'residual_after_hex': [float(after[0]).hex(), float(after[1]).hex()],
        'approximate_leading_hex': [float(delta[0]).hex(), float(delta[1]).hex()],
        'approximate_leading_decimal': [float(delta[0]), float(delta[1])],
        'approximate_is_enclosure': False,
        'leading_enclosure': {'N': enclosure['N'], 'beta': enclosure['beta']},
    }


def demonstrate(root=None, cpu_budget=CPU_BUDGET):
    if cpu_budget != CPU_BUDGET:
        raise ValueError('demonstration CPU budget is 180 seconds')
    root = ROOT if root is None else Path(root)
    started = time.process_time()
    # Stop the march before the last cell and the nodal assembly can cross 180s.
    deadline = started + cpu_budget - 2.0
    saved = load_saved_finite_cutoff_value(root)
    family, identity = require_original_history(root)
    if identity != saved['profile_identity']:
        raise ValueError('original history 0b0e4ced required')
    if float(saved['rho_up']) != float(ARCHIVE_RHO_UP):
        raise ValueError('original upstream rho changed')
    model = unit_cubic_geometry(family, bits=160)
    nodes = saved['nodes']
    order = [saved['center_index'], 0, nodes.size - 1]
    order.extend(index for index in range(nodes.size) if index not in order)
    completed = []
    transport_cpu = []
    for index in order:
        if time.process_time() >= deadline:
            break
        row = {}
        spent = []
        aborted = False
        for sign in (1, -1):
            if time.process_time() >= deadline:
                aborted = True
                break
            tick = time.process_time()
            try:
                basis = enclose_cubic_channel_basis(
                    model, arb(float(nodes[index])), sign, saved['rho_up'],
                    cells=CELLS, max_depth=16, order=ORDER, cpu_deadline=deadline)
            except CubicChannelBudgetExceeded:
                aborted = True
                break
            spent.append(time.process_time() - tick)
            row[sign] = basis
        if aborted or set(row) != {1, -1}:
            break
        completed.append((index, row))
        transport_cpu.append({'node_index': index, 'cpu_seconds': spent})
    if not completed:
        raise RuntimeError('no exact target node was enclosed within the CPU budget')
    indices = [index for index, _row in completed]
    pairs = [row for _index, row in completed]
    report = apply_cubic_leading(
        saved['residual'][indices], _stack(pairs, identity, saved['rho_up']), root,
        saved['nodes'], node_indices=indices,
        coefficient_source='enclose_cubic_channel_basis')
    if report['cubic_leading_applications'] != 1 or report['center_was_broadcast']:
        raise ValueError('cubic leading was not added once on its own nodes')
    refusal = 'certification blocked: higher remainder C_M is None'
    certified = False
    try:
        certify(report)
    except ValueError as error:
        refusal = str(error)
    else:
        certified = True
    nodes_out = [
        _node_record(saved['residual'][index], report, position)
        for position, index in enumerate(indices)]
    deltas = report['approximate_leading']
    cpu = time.process_time() - started
    return {
        'status': 'OPEN: cubic leading added on completed exact target nodes; not a UV certificate',
        'value_label': VALUE_LABEL,
        'profile_identity': identity,
        'history_path': saved['history_path'],
        'node_count': saved['node_count'],
        'completed_node_count': len(indices),
        'all_exact_target_nodes': report['all_exact_target_nodes'],
        'missing_node_count': saved['node_count'] - len(indices),
        'center_was_broadcast': False,
        'tiny_box_substituted': False,
        'cells': CELLS,
        'order': ORDER,
        'energy_folding_factor': 1,
        'large_cutoff_groups': list(report['large_cutoff_groups']),
        'saved_merit_hex': float(saved['merit']).hex(),
        'saved_merit_decimal': saved['merit'],
        'saved_constraint_maxima': saved['constraint_maxima'],
        'full_successor_merit': None,
        'max_abs_approximate_change_on_completed_nodes': float(np_max(deltas)),
        'completed_nodes': nodes_out,
        'transport_cpu_seconds': transport_cpu,
        'cpu_seconds': cpu,
        'cpu_budget': cpu_budget,
        'cpu_budget_kept': bool(cpu <= cpu_budget),
        'higher_uv_remainder_C_M': None,
        'uniform_C4_on_I': None,
        'complete_UV_tail': None,
        'certified': certified,
        'certification_refusal': refusal,
        'certification_blockers': certification_blockers(report),
        'physical_local_gate': 'OPEN',
        'physical_EXISTENCE_certificate': False,
        'physical_NONEXISTENCE_certificate': False,
        'field_or_source_evolutions': 0,
        'approximate_is_enclosure': False,
        'borrowed_inventory_sha256': _borrowed_hashes(),
        'source_identity_digest': saved['source_identity_digest'],
        'source_identity_count': saved['source_identity_count'],
        'evaluation_identity': saved['evaluation_identity'],
        'replay_merit_matches_register': saved['merit'] == saved['register_merit'],
    }


def np_max(values):
    import numpy as np
    return np.max(np.abs(values))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true', required=True)
    args = parser.parse_args()
    if not args.check:
        raise ValueError('cubic residual demonstration is check-only')
    report = demonstrate()
    if not report['cpu_budget_kept'] or report['certified'] or report['tiny_box_substituted']:
        print(json.dumps(report, indent=2, sort_keys=True, allow_nan=False))
        raise SystemExit(2)
    if report['max_abs_approximate_change_on_completed_nodes'] <= 0:
        raise RuntimeError('completed nodes did not change the value residual')
    print(json.dumps(report, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
