#!/usr/bin/env python3
"""Physical rho=1 source-bound v2: split windows, traced columns, directed maps.

This recorder authenticates the v1 OPEN foundation and the stored source
inventory. It splits the 59 straddling panels, writes row/panel manifests,
and binds directed generator/axial Lipschitz coefficients. It does not
enclose stored columns against exact physical incoming modes and does not
write a finite N,beta reconstruction budget. Historical v1 bytes stay
immutable. --check rebuilds the stable record.
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
    INVENTORY_POSITIVE_PANELS,
    INVENTORY_POSITIVE_ROWS,
    build_rho1_coverage,
    inspect_authenticated_sample,
    load_authenticated_inventory_payload,
)
from recursive_horizons.nsc_ks_rho1_source_bound_v2 import (
    CONSTRUCTION_OWNERS,
    MISSING_PRIMITIVE,
    SCHEMA,
    SUCCESSOR_STATUS,
    V1_RECORD,
    assert_v1_immutable,
    audit_row_threshold_partition,
    correlated_budget_feeds,
    coverage_digest,
    manifest_payload_bytes,
    reject_straddling_misclassification,
    resource_forecast,
    split_straddling_windows,
    trace_physical_rho1_universe,
    validate_rho1_source_bound_v2_record,
)
from recursive_horizons.nsc_ks_upstream_error import REMAINING_N_AFTER_KNOWN

OUTPUT = ROOT / 'results/development/nsc-ks-rho1-source-bound-v2.json'
PAYLOAD = ROOT / 'results/development/artifacts/nsc-ks-rho1-source-bound-v2.npz'
OWNERS = (
    'src/recursive_horizons/nsc_ks_rho1_source_bound_v2.py',
    'scripts/derive_nsc_ks_rho1_source_bound_v2.py',
    'tests/test_nsc_ks_rho1_source_bound_v2.py',
    'docs/nsc-ks-rho1-source-bound-v2.md',
)
INPUTS = (
    'results/development/nsc-ks-source-inventory.json',
    'results/development/nsc-incoming-source-update-v5.json',
    'results/development/nsc-ks-low-subgap-decision-v2.json',
    'results/development/nsc-ks-upstream-error-pilot.json',
    'results/development/nsc-incoming-retained-order24-bound.json',
    'results/development/nsc-incoming-group32-order24-bound.json',
    'results/development/nsc-mode-resolved-cauchy-state.json',
    V1_RECORD,
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
    v1_sha256 = assert_v1_immutable(ROOT)
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
    split = split_straddling_windows(arrays, coverage)
    reject_straddling_misclassification(arrays, coverage, split)
    trace = trace_physical_rho1_universe(arrays, meta, split, channels)
    sample = inspect_authenticated_sample(
        arrays, meta, coverage['remaining_windows'], meta['config'],
        panel_names=('group14/low32_1', 'group13/low'), rows_per_panel=2)
    sample['payload_sha256'] = payload_digest
    sample['full_reconstruction_enclosure_run'] = False
    sample['reconstruction_bounds_evaluated_panels'] = 0
    sample['reconstruction_bounds_evaluated_rows'] = 0
    sample['propagation_coefficients_evaluated_rows'] = trace[
        'propagation_coefficients_evaluated_rows']
    feeds = correlated_budget_feeds()
    payload_bytes = manifest_payload_bytes(trace, split, extra={
        'inventory_payload_sha256': payload_digest,
        'v1_record_sha256': v1_sha256,
    })
    payload_sha256 = sha256(payload_bytes).hexdigest()
    jsonable_windows = [
        {key: window[key] for key in (
            'group', 'sign', 'panel', 'kind', 'side', 'n_rows',
            'energy_interval', 'v5_joined', 'radius_history_response',
            'straddling_parent', 'source_reconstruction', 'error_bound',
            'v5_action_is_not_column_reconstruction', 'baseline_retained')}
        for window in split['windows']]
    if coverage_digest(split['window_identities']) != split['window_identity_digest']:
        raise ValueError('split-window identity digest changed')
    directed = {
        'N': None,
        'beta': None,
        'status': (
            'OPEN: directed N,beta requires an enclosure payload on every '
            'authenticated inventory identity'),
        'declared_universe_size': INVENTORY_POSITIVE_ROWS,
        'identity_digest': trace['identity_digest'],
        'independent_column_max_used': False,
        'reconstruction_enclosures_bound': 0,
        'propagation_coefficients_bound': INVENTORY_POSITIVE_ROWS,
        'gate_2_closed': False,
    }
    record = {
        'schema': SCHEMA,
        'accountable_author': 'Douglas Ek',
        'status': SUCCESSOR_STATUS,
        'supersedes_record_only': V1_RECORD,
        'supersedes_schema_only': 'NSC-KS-RHO1-SOURCE-BOUND-v1',
        'v1_record_sha256': v1_sha256,
        'v1_rewritten': False,
        'missing_primitive': MISSING_PRIMITIVE,
        'owners_used': {
            'source_inventory': INPUTS[0],
            'v5_source_update': INPUTS[1],
            'low_subgap_decision': INPUTS[2],
            'upstream_continuation_pilot': INPUTS[3],
            'group12_order24_action': INPUTS[4],
            'group32_order24_action': INPUTS[5],
            'channel_ledger': INPUTS[6],
            'v1_source_bound': V1_RECORD,
            'source_inventory_payload': inventory_record['payload']['path'],
        },
        'coverage': {
            'inventory_positive_panels': INVENTORY_POSITIVE_PANELS,
            'inventory_positive_rows': INVENTORY_POSITIVE_ROWS,
            'v5_action_covered_panels': coverage['v5_action_covered_panels'],
            'v5_action_covered_rows': coverage['v5_action_covered_rows'],
            'remaining_reconstruction_panels': coverage['remaining_reconstruction_panels'],
            'remaining_reconstruction_rows': coverage['remaining_reconstruction_rows'],
            'straddling_remaining_panels': split['straddling_remaining_panels'],
            'split_window_count': split['split_window_count'],
            'straddling_panels_were_split': True,
            'below_joined_rows': split['below_joined_rows'],
            'at_or_above_joined_uncovered_rows':
                split['at_or_above_joined_uncovered_rows'],
            'v5_action_rows': split['v5_action_rows'],
            'construction_census': trace['construction_census'],
            'construction_owners': CONSTRUCTION_OWNERS,
            'angular_families_with_columns': coverage['angular_families_with_columns'],
            'ell0_zero_response_groups': coverage['ell0_zero_response_groups'],
            'ell0_remaining_reconstruction_panels':
                coverage['ell0_remaining_reconstruction_panels'],
            'ell0_remaining_reconstruction_rows':
                coverage['ell0_remaining_reconstruction_rows'],
            'double_counted': False,
            'v5_action_is_not_column_reconstruction': True,
            'baseline_ell0_allocations_retained': True,
            'identities_that_are_not_error_bounds': coverage[
                'identities_that_are_not_error_bounds'],
            'row_identity_digest': trace['identity_digest'],
            'window_identity_digest': split['window_identity_digest'],
            'remaining_physical_rho1_source_error_bound': None,
            'remaining_low_subgap_source_error_bound': None,
            'directed_N_beta_contribution_bound': directed,
        },
        'row_threshold_partition': row_partition,
        'split': {
            'split_window_count': split['split_window_count'],
            'straddling_remaining_panels': split['straddling_remaining_panels'],
            'straddling_below_windows': split['straddling_below_windows'],
            'straddling_above_uncovered_windows':
                split['straddling_above_uncovered_windows'],
            'nonstraddling_remaining_windows': split['nonstraddling_remaining_windows'],
            'v5_action_windows': split['v5_action_windows'],
            'below_joined_rows': split['below_joined_rows'],
            'at_or_above_joined_uncovered_rows':
                split['at_or_above_joined_uncovered_rows'],
            'v5_action_rows': split['v5_action_rows'],
            'window_identities': split['window_identities'],
            'window_identity_digest': split['window_identity_digest'],
            'windows': jsonable_windows,
            'v5_action_is_not_column_reconstruction': True,
            'straddling_panels_were_split': True,
        },
        'propagation': {
            'axial_restriction': (
                'F_z = G(rho=1) F with homogeneous ks_generator; source energy is a label'),
            'minus_energy_factor': 2,
            'minus_energy_retained': True,
            'weights_retained': True,
            'coherence_retained': True,
            'independent_column_max_used': False,
            'v5_action_is_not_column_reconstruction': True,
            'sums': trace['propagation_sums'],
            'linear_indicator_uniform_weighted_fiber_threshold':
                trace['linear_indicator_uniform_weighted_fiber_threshold'],
            'linear_indicator_is_not_a_certificate': True,
            'documented_remaining_N_after_known_components': REMAINING_N_AFTER_KNOWN,
            'reconstruction_enclosures_evaluated_rows': 0,
            'propagation_coefficients_evaluated_rows':
                trace['propagation_coefficients_evaluated_rows'],
        },
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
        'upstream_physical_slot': feeds,
        'remaining_physical_rho1_source_error_bound': None,
        'remaining_low_subgap_source_error_bound': None,
        'directed_N_beta_contribution_bound': directed,
        'historical_coarse_4p998e_minus9_is_current_remaining_bound': False,
        'redo_group12_or_group32_order24_middle_bounds': False,
        'physical_EXISTENCE_certificate': False,
        'physical_NONEXISTENCE_certificate': False,
        'gate_2_closed': False,
        'payload': {
            'path': 'results/development/artifacts/nsc-ks-rho1-source-bound-v2.npz',
            'sha256': payload_sha256,
            'bytes': len(payload_bytes),
        },
        'resource_forecast': resource_forecast(trace),
        'next_method': MISSING_PRIMITIVE,
        'source_hashes': {path: digest(path) for path in OWNERS},
        'input_hashes': {path: digest(path) for path in INPUTS} | {
            inventory_record['payload']['path']: payload_digest,
        },
        'reproducer': 'python3 scripts/derive_nsc_ks_rho1_source_bound_v2.py --check',
        '_payload_bytes': payload_bytes,
    }
    validate_rho1_source_bound_v2_record(
        {key: value for key, value in record.items() if key != '_payload_bytes'})
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
        raise ValueError('unexpected rho=1 source-bound v2 schema')
    inventory_record = json.loads((ROOT / INPUTS[0]).read_text())
    ensure_inventory_payload(inventory_record)
    assert_v1_immutable(ROOT)
    for path, expected in {**record['source_hashes'], **record['input_hashes']}.items():
        if digest(path) != expected:
            raise ValueError('rho=1 source-bound v2 dependency changed: ' + path)
    validate_rho1_source_bound_v2_record(record)
    if record['remaining_physical_rho1_source_error_bound'] is not None:
        raise ValueError('recorder may not invent the missing physical rho=1 bound')
    if record['directed_N_beta_contribution_bound']['N'] is not None:
        raise ValueError('recorder may not claim a finite directed N budget')
    if record['historical_coarse_4p998e_minus9_is_current_remaining_bound']:
        raise ValueError('historical 4.998e-9 coarse aggregate cannot be current')
    if record['v1_rewritten'] or record['v1_record_sha256'] != digest(V1_RECORD):
        raise ValueError('historical v1 source-bound bytes changed')
    spec = record['payload']
    path = ROOT / spec['path']
    if sha256(path.read_bytes()).hexdigest() != spec['sha256']:
        raise ValueError('rho=1 source-bound v2 payload changed')


def _public(record):
    return {key: value for key, value in record.items() if key != '_payload_bytes'}


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
        if stable(_public(record)) != stable(_public(fresh)):
            raise ValueError('rho=1 source-bound v2 replay differs')
        if sha256(fresh['_payload_bytes']).hexdigest() != record['payload']['sha256']:
            raise ValueError('rho=1 source-bound v2 payload replay differs')
    else:
        record = calculate()
        public = _public(record)
        output = (json.dumps(
            public, indent=2, sort_keys=True, allow_nan=False) + '\n').encode()
        if args.record:
            PAYLOAD.parent.mkdir(parents=True, exist_ok=True)
            if OUTPUT.is_file() or PAYLOAD.is_file():
                if OUTPUT.is_file() and OUTPUT.read_bytes() != output:
                    raise ValueError('rho=1 source-bound v2 record already differs')
                if PAYLOAD.is_file() and PAYLOAD.read_bytes() != record['_payload_bytes']:
                    raise ValueError('rho=1 source-bound v2 payload already differs')
            else:
                write_bytes_atomic(PAYLOAD, record['_payload_bytes'], exclusive=True)
                write_bytes_atomic(OUTPUT, output, exclusive=True)
        elif args.replay:
            saved = json.loads(OUTPUT.read_text())
            if stable(public) != stable(saved):
                raise ValueError('rho=1 source-bound v2 replay differs')
        check(public)
    shown = record if args.check else public
    print(json.dumps({
        'status': shown['status'],
        'schema': SCHEMA,
        'inventory_positive_panels': INVENTORY_POSITIVE_PANELS,
        'inventory_positive_rows': INVENTORY_POSITIVE_ROWS,
        'straddling_remaining_panels': shown['coverage']['straddling_remaining_panels'],
        'split_window_count': shown['coverage']['split_window_count'],
        'directed_N_beta_contribution_bound': shown['directed_N_beta_contribution_bound'],
        'remaining_physical_rho1_source_error_bound': None,
        'gate_2_closed': False,
        'runtime': shown.get('runtime'),
    }, indent=2, sort_keys=True, allow_nan=False))
