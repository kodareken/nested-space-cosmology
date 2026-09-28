#!/usr/bin/env python3
"""Phase 2A rho=1 source-bound core: coverage, validator and a small sample.

This recorder authenticates the V5 inventory and the existing
upstream/low-subgap owners. It does not reconstruct all 102 remaining
panels / 29356 rows and does not write a finite physical N,beta budget.
"""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'

import argparse
import ctypes
from hashlib import sha256
import json
from pathlib import Path
import shutil
import sys
import time

try:
    import resource
except ModuleNotFoundError:  # Windows has no resource module.
    resource = None

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from recursive_horizons.nsc_ks_evaluation_binding import write_bytes_atomic
from recursive_horizons.nsc_ks_rho1_source_bound import (
    SCHEMA,
    audit_row_threshold_partition,
    build_rho1_coverage,
    inspect_authenticated_sample,
    load_authenticated_inventory_payload,
    physical_source_error_for_upstream,
    validate_rho1_source_bound_record,
)

OUTPUT = ROOT / 'results/development/nsc-ks-rho1-source-bound-v1.json'
OWNERS = (
    'src/recursive_horizons/nsc_ks_rho1_source_bound.py',
    'scripts/derive_nsc_ks_rho1_source_bound_v1.py',
    'tests/test_nsc_ks_rho1_source_bound.py',
    'docs/nsc-ks-rho1-source-bound-v1.md',
)
INPUTS = (
    'results/development/nsc-ks-source-inventory.json',
    'results/development/nsc-incoming-source-update-v5.json',
    'results/development/nsc-ks-low-subgap-decision-v2.json',
    'results/development/nsc-ks-upstream-error-pilot.json',
    'results/development/nsc-incoming-retained-order24-bound.json',
    'results/development/nsc-incoming-group32-order24-bound.json',
    'results/development/nsc-mode-resolved-cauchy-state.json',
)
EXTERNAL_PAYLOAD = os.environ.get('NSC_SOURCE_INVENTORY_PAYLOAD')


def peak_rss():
    """Return a process-lifetime resident-memory peak and its unit/source."""
    if resource is not None:
        value = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        unit = 'bytes' if sys.platform == 'darwin' else 'kilobytes'
        return value, unit, 'resource.ru_maxrss'
    if os.name == 'nt':
        class ProcessMemoryCounters(ctypes.Structure):
            _fields_ = (
                ('cb', ctypes.c_ulong),
                ('PageFaultCount', ctypes.c_ulong),
                ('PeakWorkingSetSize', ctypes.c_size_t),
                ('WorkingSetSize', ctypes.c_size_t),
                ('QuotaPeakPagedPoolUsage', ctypes.c_size_t),
                ('QuotaPagedPoolUsage', ctypes.c_size_t),
                ('QuotaPeakNonPagedPoolUsage', ctypes.c_size_t),
                ('QuotaNonPagedPoolUsage', ctypes.c_size_t),
                ('PagefileUsage', ctypes.c_size_t),
                ('PeakPagefileUsage', ctypes.c_size_t),
            )
        counters = ProcessMemoryCounters()
        counters.cb = ctypes.sizeof(counters)
        get_process = ctypes.windll.kernel32.GetCurrentProcess
        get_process.argtypes = []
        get_process.restype = ctypes.c_void_p
        get_memory = ctypes.windll.psapi.GetProcessMemoryInfo
        get_memory.argtypes = (
            ctypes.c_void_p, ctypes.POINTER(ProcessMemoryCounters),
            ctypes.c_ulong)
        get_memory.restype = ctypes.c_int
        process = get_process()
        ok = get_memory(process, ctypes.byref(counters), counters.cb)
        if not ok:
            raise OSError('GetProcessMemoryInfo failed')
        return int(counters.PeakWorkingSetSize), 'bytes', (
            'GetProcessMemoryInfo.PeakWorkingSetSize')
    return None, 'unavailable', 'unavailable'


def digest(path):
    return sha256((ROOT / path).read_bytes()).hexdigest()


def ensure_inventory_payload(inventory_record):
    spec = inventory_record['payload']
    dest = ROOT / spec['path']
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.is_file() and sha256(dest.read_bytes()).hexdigest() == spec['sha256']:
        return dest
    external = Path(EXTERNAL_PAYLOAD) if EXTERNAL_PAYLOAD else None
    if external is not None and external.is_file():
        shutil.copyfile(external, dest)
    if not dest.is_file() or sha256(dest.read_bytes()).hexdigest() != spec['sha256']:
        raise FileNotFoundError('authenticated source inventory payload required')
    return dest


def stable(record):
    return {key: value for key, value in record.items() if key != 'runtime'}


def calculate():
    started_cpu, started_wall = time.process_time(), time.monotonic()
    started_rss, rss_units, rss_source = peak_rss()
    inventory_record = json.loads((ROOT / INPUTS[0]).read_text())
    v5 = json.loads((ROOT / INPUTS[1]).read_text())
    decision = json.loads((ROOT / INPUTS[2]).read_text())
    pilot = json.loads((ROOT / INPUTS[3]).read_text())
    group12 = json.loads((ROOT / INPUTS[4]).read_text())
    group32 = json.loads((ROOT / INPUTS[5]).read_text())
    channels = json.loads((ROOT / INPUTS[6]).read_text())['channels']
    ensure_inventory_payload(inventory_record)
    arrays, meta, payload_digest = load_authenticated_inventory_payload(
        ROOT, inventory_record)
    coverage = build_rho1_coverage(
        decision['coverage']['windows'],
        inventory_record=inventory_record,
        v5_coverage=v5['coverage'],
        channels=channels)
    row_partition = audit_row_threshold_partition(arrays, coverage)
    sample = inspect_authenticated_sample(
        arrays, meta, coverage['remaining_windows'], meta['config'],
        panel_names=('group14/low32_1', 'group13/low'), rows_per_panel=2)
    sample['payload_sha256'] = payload_digest
    upstream = physical_source_error_for_upstream(None)
    record = {
        'schema': SCHEMA,
        'accountable_author': 'Douglas Ek',
        'status': (
            'OPEN: rho=1 source-bound core and validator; complete 164-panel '
            'physical column reconstruction and directed N,beta budget unevaluated; '
            '102-panel low/subgap subset identified'),
        'owners_used': {
            'source_inventory': INPUTS[0],
            'v5_source_update': INPUTS[1],
            'low_subgap_decision': INPUTS[2],
            'upstream_continuation_pilot': INPUTS[3],
            'group12_order24_action': INPUTS[4],
            'group32_order24_action': INPUTS[5],
            'channel_ledger': INPUTS[6],
            'source_inventory_payload': inventory_record['payload']['path'],
        },
        'coverage': {
            key: coverage[key] for key in (
                'inventory_positive_panels', 'inventory_positive_rows',
                'v5_action_covered_panels', 'v5_action_covered_rows',
                'remaining_reconstruction_panels', 'remaining_reconstruction_rows',
                'remaining_by_kind', 'angular_families_with_columns',
                'ell0_zero_response_groups',
                'ell0_remaining_reconstruction_panels',
                'ell0_remaining_reconstruction_rows',
                'remaining_panel_identities', 'v5_action_panel_identities',
                'remaining_identity_digest', 'v5_action_identity_digest',
                'double_counted', 'v5_action_is_not_column_reconstruction',
                'baseline_ell0_allocations_retained',
                'identities_that_are_not_error_bounds',
                'remaining_physical_rho1_source_error_bound',
                'remaining_low_subgap_source_error_bound',
                'directed_N_beta_contribution_bound')
        },
        'row_threshold_partition': row_partition,
        'sample': sample,
        'reused_successors': {
            'group12_order24_lapse_error_upper':
                group12['energy']['lapse_action_error_upper'],
            'group32_combined_action_error_upper':
                group32['windows']['combined']['conditional_constraint_action_error_upper'],
            'group32_validated_lower_window': v5['coverage']['32']['validated_lower_window'],
            'v5_covered_action_bound_N_beta':
                v5['partial_error_budget']['covered_spectral_regions_action_error_upper'],
            'action_bounds_are_not_column_reconstruction': True,
        },
        'upstream_physical_slot': {
            'physical_rho1_source_error_upper': pilot['physical_rho1_source_error_upper'],
            'combined_upstream_bound_N_beta':
                upstream['nine_component_feed']['combined_upstream_bound_N_beta'],
            'combined_status': upstream['nine_component_feed']['combined_status'],
        },
        'remaining_physical_rho1_source_error_bound': None,
        'remaining_low_subgap_source_error_bound': None,
        'directed_N_beta_contribution_bound': coverage['directed_N_beta_contribution_bound'],
        'historical_coarse_4p998e_minus9_is_current_remaining_bound': False,
        'redo_group12_or_group32_order24_middle_bounds': False,
        'physical_EXISTENCE_certificate': False,
        'physical_NONEXISTENCE_certificate': False,
        'next_method': (
            'split the 59 straddling panel-level windows at their V5 joined '
            'thresholds, then produce authenticated directed rho=1 physical-column '
            'reconstruction errors for the complete 164-panel / 49372-row source '
            'universe; V5 action bounds are not column-accuracy bounds'),
        'source_hashes': {path: digest(path) for path in OWNERS},
        'input_hashes': {path: digest(path) for path in INPUTS} | {
            inventory_record['payload']['path']: payload_digest,
        },
        'reproducer': 'python3 scripts/derive_nsc_ks_rho1_source_bound_v1.py --check',
    }
    validate_rho1_source_bound_record(record)
    ended_rss, ended_units, ended_source = peak_rss()
    if (ended_units, ended_source) != (rss_units, rss_source):
        raise RuntimeError('resident-memory measurement source changed')
    record['runtime'] = {
        'CPU_seconds': time.process_time() - started_cpu,
        'wall_seconds': time.monotonic() - started_wall,
        'ru_maxrss': ended_rss,
        'ru_maxrss_at_start': started_rss,
        'rss_units': rss_units,
        'rss_source': rss_source,
    }
    return record


def check(record):
    if record.get('schema') != SCHEMA:
        raise ValueError('unexpected rho=1 source-bound schema')
    inventory_record = json.loads((ROOT / INPUTS[0]).read_text())
    ensure_inventory_payload(inventory_record)
    for path, expected in {**record['source_hashes'], **record['input_hashes']}.items():
        if digest(path) != expected:
            raise ValueError('rho=1 source-bound dependency changed: ' + path)
    validate_rho1_source_bound_record(record)
    if record['remaining_physical_rho1_source_error_bound'] is not None:
        raise ValueError('recorder may not invent the missing physical rho=1 bound')
    if record['directed_N_beta_contribution_bound']['N'] is not None:
        raise ValueError('recorder may not claim a finite directed N budget')
    if record['historical_coarse_4p998e_minus9_is_current_remaining_bound']:
        raise ValueError('historical 4.998e-9 coarse aggregate cannot be current')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--record', action='store_true')
    mode.add_argument('--check', action='store_true')
    mode.add_argument('--replay', action='store_true')
    args = parser.parse_args()
    if args.check:
        record = json.loads(OUTPUT.read_text())
        check(record)
        fresh = calculate()
        if stable(record) != stable(fresh):
            raise ValueError('rho=1 source-bound replay differs')
    else:
        record = calculate()
        if args.record:
            output = (json.dumps(
                record, indent=2, sort_keys=True, allow_nan=False) + '\n').encode()
            if OUTPUT.is_file():
                if OUTPUT.read_bytes() != output:
                    raise ValueError('rho=1 source-bound record already differs')
            else:
                write_bytes_atomic(OUTPUT, output, exclusive=True)
        elif args.replay:
            saved = json.loads(OUTPUT.read_text())
            if stable(record) != stable(saved):
                raise ValueError('rho=1 source-bound replay differs')
        check(record)
    print(json.dumps({
        'status': record['status'],
        'remaining_reconstruction_panels':
            record['coverage']['remaining_reconstruction_panels'],
        'remaining_reconstruction_rows':
            record['coverage']['remaining_reconstruction_rows'],
        'v5_action_covered_panels': record['coverage']['v5_action_covered_panels'],
        'v5_action_covered_rows': record['coverage']['v5_action_covered_rows'],
        'inspected_sample_panels': record['sample']['inspected_sample_panels'],
        'inspected_sample_rows': record['sample']['inspected_sample_rows'],
        'directed_N_beta_contribution_bound':
            record['directed_N_beta_contribution_bound'],
        'remaining_physical_rho1_source_error_bound':
            record['remaining_physical_rho1_source_error_bound'],
        'runtime': record.get('runtime'),
    }, indent=2, sort_keys=True, allow_nan=False))
