#!/usr/bin/env python3
"""Authenticate and contract the new incoming source panels, never old modes."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from recursive_horizons.nsc_incoming_source_assembly import assemble, PANEL_PATH, PANEL_DIGEST
from recursive_horizons.nsc_incoming_spectral_panels import sha, record, compute_panels

OUTPUT = 'results/development/nsc-incoming-spectral-source.json'
SOURCES = ('src/recursive_horizons/nsc_incoming_source_assembly.py',
    'src/recursive_horizons/nsc_incoming_spectral_panels.py',
    'src/recursive_horizons/nsc_incoming_state_moments.py',
    'src/recursive_horizons/nsc_incoming_high_energy_source.py',
    'tests/test_nsc_incoming_source_assembly.py', 'tests/test_nsc_incoming_state_moments.py',
    'tests/test_nsc_incoming_high_energy_source.py', 'docs/nsc-incoming-spectral-source.md',
    'scripts/derive_nsc_incoming_spectral_source.py')
INPUTS = ('results/development/nsc-incoming-constraint-gate.json',
    'results/development/nsc-magnetic-light-restoration.json',
    'src/recursive_horizons/nsc_magnetic_light_reference.py',
    'src/recursive_horizons/nsc_horizon_source.py',
    'scripts/derive_nsc_transmitting_boundary_binding.py')


def make_record():
    result = assemble(ROOT)
    meta = result.pop('payload_metadata')
    checks = dict(result['residuals'])
    checks.update({k: v for k, v in result['low_modal_diagnostics'].items() if k != 'raw_flux_discrepancy'})
    tolerance = 3e-11
    failures = {k: v for k, v in checks.items() if not np.isfinite(v) or v > tolerance}
    return {'schema': 'NSC-INCOMING-SPECTRAL-SOURCE-v1', 'accountable_author': 'Douglas Ek',
        'status': 'PASS: finite incoming source assembly; converged full source OPEN' if not failures else 'OPEN',
        'source_hashes': {p: sha(ROOT, p) for p in SOURCES},
        'input_hashes': {p: sha(ROOT, p) for p in INPUTS},
        'numerical_owner_and_input_hashes': meta['signature'], 'input_payloads': meta['input_payloads'],
        'payload': {'path': PANEL_PATH, 'sha256': PANEL_DIGEST, 'bytes': (ROOT/PANEL_PATH).stat().st_size},
        'locked_inputs': record(ROOT, INPUTS[0])['locked_inputs'],
        'panel_inventory': meta['panels'], **result, 'acceptance_tolerance': tolerance,
        'gate_residuals': checks, 'failures': failures,
        'domain': {'geometry': meta['geometry'], 'canonical_frame': 'same incoming KS normal and half-density',
            'state': 'unchanged global horizon/incoming source; seed covariance control only',
            'selected_history': False, 'changed_normal_jets_evaluated': False,
            'old_generators_rerun': False, 'compact_induced_terms': 'excluded, remain on local-action side',
            'refined_panel_selection': 'explicit new observable; frozen packet result not changed',
            'error_indicators': 'arithmetic and measured quadrature differences, not total or rigorous source bounds'},
        'gate': {'finite_panel_source': 'PASS' if not failures else 'OPEN',
            'subgap_source_refinement': 'OPEN', 'infinite_energy_tail': 'OPEN',
            'retained_angular_compact_remainder': 'OPEN', 'full_renormalized_source': 'OPEN',
            'incoming_constraints': 'OPEN', 'extended_stationarity': 'OPEN',
            'physical_EndpointBranchJets': 'OPEN', 'updated_constraints': None, 'selected_nulls': None,
            'metric_timestep': False, 'Z3': 'OUT OF SCOPE', 'PDF_bumped': False,
            'homogeneous_nonexistence': 'preserved in original scope; not rerun',
            'Weyl_time_node_diagnostic': 93.54264532195464},
        'comparison': {'fields': 'all', 'float_atol': 3e-13, 'float_rtol': 3e-13,
                       'exact': 'hashes, masks, labels, inventory and scope'},
        'reproducer': 'python3 scripts/derive_nsc_incoming_spectral_source.py --check'}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--prepare-panels', action='store_true')
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    if args.prepare_panels:
        from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes
        arrays = compute_panels(ROOT, progress=lambda i, n: print(f'new source arithmetic {i}/{n}', flush=True))
        raw = deterministic_npz_bytes(arrays); digest = hashlib.sha256(raw).hexdigest()
        path = ROOT/f'results/development/artifacts/nsc-incoming-spectral-panels.{digest}.npz'
        if path.exists() and path.read_bytes() != raw:
            raise ValueError('existing content-addressed artifact differs')
        if not path.exists(): path.write_bytes(raw)
        print(json.dumps({'path': str(path), 'sha256': digest, 'bytes': len(raw)}))
        return
    result = make_record()
    if args.check:
        from derive_nsc_transmitting_boundary_binding import compare
        compare(record(ROOT, OUTPUT), result)
    else:
        (ROOT/OUTPUT).write_text(json.dumps(result, sort_keys=True, indent=2)+'\n')
    print(json.dumps({k: result[k] for k in ('status', 'finite_quantum_source_approximant',
        'net_quadrature_change', 'residuals', 'failures')}, indent=2))
    if result['failures']: raise SystemExit(1)


if __name__ == '__main__': main()
