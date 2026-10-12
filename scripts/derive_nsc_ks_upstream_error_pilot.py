#!/usr/bin/env python3
"""Directed family14 numerical-continuation pilot; physical rho=1 error stays OPEN."""
import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import sys
import time

import numpy as np
from flint import arb, ctx

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
sys.path.insert(0, str(ROOT / 'scripts'))
import derive_nsc_ks_source_control_v2 as C
from recursive_horizons.evidence_io import publish_exclusive_file
from recursive_horizons.nsc_ks_ball_trajectory import exact_upper
from recursive_horizons.nsc_ks_retained_upstream_archive import RetainedUpstreamArchive
from recursive_horizons.nsc_ks_source_inventory import RetainedSourceInventory
from recursive_horizons.nsc_ks_upstream_error import (
    ASSUMPTIONS, DIRECTED_FORMULA, SCHEMA_CONTINUATION,
    combine_upstream_row_error, continue_columns, independent_constant_matrix_control,
    nine_component_upstream_feed, pack_upper, packed_upper_fraction,
    restore_packed_upper, weighted_column_errors)

OUTPUT = ROOT / 'results/development/nsc-ks-upstream-error-pilot.json'
OWNERS = (
    'scripts/derive_nsc_ks_upstream_error_pilot.py',
    'src/recursive_horizons/nsc_ks_upstream_error.py',
    'tests/test_nsc_ks_upstream_error.py',
    'docs/nsc-ks-upstream-error.md',
)


def digest(path):
    return sha256((ROOT / path).read_bytes()).hexdigest()


def json_bounds(value):
    if isinstance(value, Fraction):
        return exact_upper(arb(value.numerator) / value.denominator)
    if isinstance(value, dict):
        return {key: json_bounds(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_bounds(item) for item in value]
    return value


def inputs():
    archive = RetainedUpstreamArchive(ROOT)
    source, _negative, archived, first, channel, selections = C.select_source(
        archive, (14, 1))
    unique = source.energies.reshape(-1, 3)[:, 0]
    at_one = []
    panel_names = []
    with RetainedSourceInventory(ROOT) as inventory:
        panels = {panel.name: panel for panel in inventory.positive_panels(groups=(14,))
                  if panel.angular_sign == 1}
        for energy in unique:
            panel = next(value for value in panels.values()
                         if np.any(value.energies == energy))
            index = np.flatnonzero(panel.energies == energy)
            if len(index) != 1:
                raise ValueError('unique family14 rho=1 source row required')
            at_one.append(panel.amplitudes_at_one[index[0]])
            panel_names.append(panel.name)
        inventory_hashes = {item['path']: item['sha256'] for item in inventory.provenance}
    return (archive, source, archived, first, channel, selections, unique,
            np.asarray(at_one), panel_names, inventory_hashes)


def endpoint_difference_columns(calculated, archived, *, bits=100):
    calculated = np.asarray(calculated, complex).reshape(2, -1)
    archived = np.asarray(archived, complex)
    if calculated.shape != archived.shape:
        raise ValueError('matching calculated and archived endpoint columns required')
    with ctx.workprec(bits):
        result = []
        for column in range(calculated.shape[1]):
            total = arb(0)
            for spin in range(2):
                value = calculated[spin, column] - archived[spin, column]
                from recursive_horizons.nsc_ks_ball_trajectory import complex_ball
                total += complex_ball(value).abs_upper() ** 2
            result.append(total.sqrt().upper())
        return result


def calculate():
    started_cpu, started_wall = time.process_time(), time.monotonic()
    archive, source, archived, first, channel, selections, energies, at_one, panel_names, inventory_hashes = inputs()
    validation = continue_columns(
        at_one, energies, first.mass, first.angular,
        max_step=1/4096, panels=1, bits=100)
    archive_delta = endpoint_difference_columns(
        validation['amplitudes_at_rho_up'], archived, bits=100)
    per_column = []
    for proof, delta in zip(validation['transported_column_errors'], archive_delta):
        per_column.append(pack_upper(restore_packed_upper(proof) + delta))
    with ctx.workprec(100):
        total = sum((restore_packed_upper(value) ** 2 for value in per_column), arb(0)).sqrt().upper()
    weights = source.column_weights
    labels = source.energies
    weighted_field, weighted_axial = weighted_column_errors(
        per_column, weights, labels, bits=100)
    endpoint_indicator = float(np.max(abs(
        validation['amplitudes_at_rho_up'].reshape(2, -1) - archived)))
    numerical_row = packed_upper_fraction(pack_upper(total))
    return {
        'schema': SCHEMA_CONTINUATION,
        'status': (
            'OPEN: directed numerical continuation bound complete for family14 control; '
            'physical rho=1 source error and history-correlated constraint transport remain missing'),
        'family': [14, 1],
        'source_rows': selections,
        'rho1_panels': panel_names,
        'energies': energies.tolist(),
        'rho_up': first.rho_up,
        'mass': first.mass,
        'angular': first.angular,
        'signed_channel_multiplicity': channel['copy_count'] * channel['degeneracy'],
        'directed_formula': DIRECTED_FORMULA,
        'assumptions': list(ASSUMPTIONS),
        'solver': {'integrator': 'DOP853', 'rtol': 5e-14, 'atol': 5e-16,
                   'max_step': 1/4096, 'defect_method': 'Bernstein after cancellation'},
        'accepted_steps': validation['accepted_steps'],
        'function_evaluations': validation['nfev'],
        'validation_endpoint_error_upper': validation['unweighted_frobenius_error_upper'],
        'archive_endpoint_difference_indicator': endpoint_indicator,
        'archive_column_error_uppers': per_column,
        'archive_unweighted_frobenius_error_upper': pack_upper(total),
        'weighted_field_error_upper': exact_upper(arb(weighted_field.numerator) / weighted_field.denominator),
        'weighted_axial_error_upper': exact_upper(arb(weighted_axial.numerator) / weighted_axial.denominator),
        'physical_rho1_source_error_upper': None,
        'combined_upstream_row_error': json_bounds(
            combine_upstream_row_error(numerical_row, None)),
        'nine_component_feed': nine_component_upstream_feed(None, None),
        'constraint_transport': (
            'must use the correlated reference/difference equations on the final history; '
            'a full-reference matter bound is intentionally not substituted'),
        'constant_matrix_control': independent_constant_matrix_control(),
        'refinement_used_as_bound': False,
        'physical_EXISTENCE_certificate': False,
        'source_hashes': {path: digest(path) for path in OWNERS},
        'input_hashes': {
            **archive.input_hashes,
            **inventory_hashes,
            C.HISTORY: digest(C.HISTORY),
        },
        'runtime': {'CPU_seconds': time.process_time() - started_cpu,
                    'wall_seconds': time.monotonic() - started_wall,
                    'CPU_cap': 180},
    }


def check(record):
    if record.get('schema') != SCHEMA_CONTINUATION:
        raise ValueError('unexpected upstream pilot schema')
    for path, expected in {**record['source_hashes'], **record['input_hashes']}.items():
        if digest(path) != expected:
            raise ValueError('upstream pilot dependency changed: ' + path)
    if record['physical_rho1_source_error_upper'] is not None:
        raise ValueError('pilot may not invent the physical rho=1 source error')
    if record['nine_component_feed']['combined_upstream_bound_N_beta'] is not None:
        raise ValueError('uncorrelated pilot cannot fill the gate upstream slot')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument('--record', action='store_true')
    modes.add_argument('--check', action='store_true')
    modes.add_argument('--replay', action='store_true')
    args = parser.parse_args()
    if args.check:
        record = json.loads(OUTPUT.read_text())
        check(record)
    else:
        record = calculate()
        if args.record:
            publish_exclusive_file(ROOT, str(OUTPUT.relative_to(ROOT)),
                (json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + '\n').encode())
        elif args.replay:
            saved = json.loads(OUTPUT.read_text())
            stable = lambda value: {key: item for key, item in value.items() if key != 'runtime'}
            if stable(record) != stable(saved):
                raise ValueError('upstream pilot replay differs')
    print(json.dumps({key: record[key] for key in (
        'status', 'archive_unweighted_frobenius_error_upper',
        'weighted_field_error_upper', 'weighted_axial_error_upper',
        'physical_rho1_source_error_upper', 'runtime')}, indent=2))
