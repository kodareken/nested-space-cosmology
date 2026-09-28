#!/usr/bin/env python3
"""Materialize existing finite source inputs; no field ODE or history solve."""
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
sys.path.insert(0, str(ROOT / 'src'))
from recursive_horizons.nsc_ks_source_inventory import RetainedSourceInventory
from recursive_horizons.nsc_incoming_spectral_panels import NUMERICAL_OWNERS, RECORDS
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes

OUTPUT = 'results/development/nsc-ks-source-inventory.json'
PAYLOAD = 'results/development/artifacts/nsc-ks-source-inventory.npz'
OWNERS = tuple(dict.fromkeys((
    'src/recursive_horizons/nsc_ks_source_inventory.py', 'docs/nsc-ks-source-inventory.md',
    'scripts/derive_nsc_ks_source_inventory.py', 'src/recursive_horizons/nsc_incoming_source_assembly.py',
    *NUMERICAL_OWNERS, *RECORDS)))


def digest(path):
    return sha256((ROOT / path).read_bytes()).hexdigest()


def summarize(arrays, metadata, payload):
    counts, measures = {}, {}
    for name, desc in metadata['panels'].items():
        key = f"{desc['group']}_{desc['angular_sign']}"
        counts[key] = counts.get(key, 0) + len(arrays[name + '/energies'])
        measures[key] = measures.get(key, 0.) + float(arrays[name + '/weights'].sum())
    measure_residual = max(abs(value - (320. if int(k.split('_')[0]) in (10, 11, 12, 31, 32) else 160.))
                           for k, value in measures.items())
    if len(counts) != 62 or measure_residual > 1e-11:
        raise ValueError('finite disjoint source inventory is incomplete or double counted')
    return {
        'schema': 'NSC-KS-SOURCE-INVENTORY-v1', 'accountable_author': 'Douglas Ek',
        'status': 'PASS: finite runtime input inventory; physical local gate OPEN',
        'positive_energy_rows': sum(counts.values()), 'positive_panels': len(metadata['panels']),
        'angular_families_with_columns': len(counts), 'rows_by_family': counts,
        'quadrature_measure_by_family': measures, 'quadrature_measure_residual': measure_residual,
        'quadrature_measure_tolerance': 1e-11,
        'negative_energy_rule': 'same weights; S3 conjugation; opposite angular sign; owned negative source law',
        'negative_source_identity_residual': metadata['negative_source_identity_residual'],
        'negative_source_identity_tolerance': 3e-13,
        'radius_response_exactly_zero_groups': [0, 13, 23],
        'scope': {'full_source_error_bound': None, 'physical_local_gate': 'OPEN',
                  'low_subgap_error_bound': None, 'changed_history_UV_bound': None,
                  'incoming_state_frozen': False, 'history_solves': 0, 'source_ODE_solves': 0,
                  'reference_continuation_performed': False, 'new_source_scales': False,
                  'baseline_ell0_allocations_retained': True, 'metric_timestep': False},
        'source_hashes': metadata['signature'], 'input_payloads': metadata['input_payloads'],
        'payload': payload, 'reproducer': 'python3 scripts/derive_nsc_ks_source_inventory.py --check',
    }


def prepare():
    if (ROOT / OUTPUT).exists() or (ROOT / PAYLOAD).exists():
        raise FileExistsError('source inventory exists; use --check')
    signature = {p: digest(p) for p in OWNERS}
    arrays, panels, residual = {}, {}, 0.
    with RetainedSourceInventory(ROOT) as inventory:
        for panel in inventory.positive_panels():
            for key in ('energies', 'weights', 'covariance', 'amplitudes_at_one'):
                arrays[panel.name + '/' + key] = getattr(panel, key)
            partner = panel.negative_partner(inventory.config)
            residual = max(residual, float(np.max(abs(partner.covariance - (np.eye(3) - panel.covariance.conj())))))
            panels[panel.name] = {k: getattr(panel, k) for k in ('group', 'angular_sign', 'mass', 'angular_magnitude', 'provenance')}
        metadata = {'signature': signature, 'panels': panels, 'input_payloads': inventory.provenance,
                    'config': inventory.config, 'channels': inventory.channels,
                    'negative_source_identity_residual': residual}
    if residual > 3e-13 or signature != {p: digest(p) for p in OWNERS}:
        raise ValueError('source identity/owner changed')
    arrays['metadata_json'] = np.frombuffer(json.dumps(metadata, sort_keys=True).encode(), np.uint8)
    raw = deterministic_npz_bytes(arrays)
    item = {'path': PAYLOAD, 'sha256': sha256(raw).hexdigest(), 'bytes': len(raw)}
    result = summarize(arrays, metadata, item)
    (ROOT / PAYLOAD).write_bytes(raw)
    (ROOT / OUTPUT).write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + '\n')
    return result


def check():
    saved = json.loads((ROOT / OUTPUT).read_text())
    if digest(PAYLOAD) != saved['payload']['sha256']:
        raise ValueError('source inventory payload changed')
    with np.load(ROOT / PAYLOAD, allow_pickle=False) as loaded:
        arrays = {k: loaded[k] for k in loaded.files}
    meta = json.loads(arrays['metadata_json'].tobytes())
    if meta['signature'] != {p: digest(p) for p in OWNERS}:
        raise ValueError('source inventory implementation/input changed')
    for item in meta['input_payloads']:
        if digest(item['path']) != item['sha256']:
            raise ValueError('source input payload changed')
    if summarize(arrays, meta, saved['payload']) != saved:
        raise ValueError('source inventory array replay differs')
    return saved


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--prepare', action='store_true')
    group.add_argument('--check', action='store_true')
    result = prepare() if parser.parse_args().prepare else check()
    print(json.dumps({k: result[k] for k in ('status', 'positive_energy_rows', 'positive_panels',
        'angular_families_with_columns', 'quadrature_measure_residual', 'negative_source_identity_residual')}, indent=2))
