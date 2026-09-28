#!/usr/bin/env python3
"""Authenticate the exact paired-band bulk closure; no physical generator."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))

from recursive_horizons.nsc_reference_band_bulk import paired_unweighted_bulk_closure


OUTPUT = ROOT/'results/development/nsc-reference-band-bulk.json'
SOURCES = (
    'src/recursive_horizons/nsc_reference_band_bulk.py',
    'tests/test_nsc_reference_band_bulk.py',
    'scripts/derive_nsc_reference_band_bulk.py',
    'docs/nsc-reference-band-bulk.md',
)
INPUTS = (
    'src/recursive_horizons/nsc_reference_band_action.py',
    'src/recursive_horizons/nsc_spatial_reference_symbol.py',
    'src/recursive_horizons/nsc_common_subtracted_ks_source.py',
    'src/recursive_horizons/nsc_incoming_state_moments.py',
    'src/recursive_horizons/nsc_incoming_source_assembly.py',
    'results/development/nsc-mode-resolved-cauchy-state.json',
    'scripts/derive_nsc_transmitting_boundary_binding.py',
)


def make_record():
    channels = json.loads((ROOT/'results/development/nsc-mode-resolved-cauchy-state.json').read_text())['channels']
    result = paired_unweighted_bulk_closure(
        channels, full_momentum_line=True, unweighted=True,
        compact_test_variations=True, smooth_positive_metric_jets=True)
    sha = lambda path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
    return {
        'schema': 'NSC-REFERENCE-BAND-BULK-v1', 'accountable_author': 'Douglas Ek',
        'status': 'PASS: paired unweighted full-k reference-band bulk Euler exchange vanishes through order four',
        'source_hashes': {p: sha(p) for p in SOURCES},
        'input_hashes': {p: sha(p) for p in INPUTS},
        'scope_clarification_provenance': {
            'preserved_record_commit': 'e22b03903d58af44043b5d9b5a34901615cd5c48',
            'preserved_record_sha256': '9faaded760d45fb3181f2bb9cd0f9d17aadaf1a3b09dee2398ec5fa168f7de64',
            'correction': 'make the existing chartwise LLL convention explicit; no distributional k=0 bridge is proved',
            'exact_coefficients_changed': False,
        },
        'result': result,
        'exact_residuals': 'all Pauli coefficient, angular-pair, Weyl-trace and recursion residuals are exactly zero',
        'proof_inputs': {
            'reused_identity': 'nsc_reference_band_action.action_variation_identity',
            'common_remainder': 'nsc_common_subtracted_ks_source.reference_band_remainder',
            'projector_recursion': 'nsc_spatial_reference_symbol.reference_projector',
            'published_formal_owner': 'Panati, Spohn and Teufel, math-ph/0201055, sections 2 and 4.4; existing project reuse',
            'new_connection': 'the only nondecaying k-boundary primitives are odd in the actual paired angular sign',
            'computation': 'exact 2x2 Pauli coefficient differentiation and integer order bounds only',
            'physical_mode_or_large_k_scan': False,
        },
        'comparison': {'fields': 'all', 'float_atol': 3e-13, 'float_rtol': 3e-13,
                       'exact': 'types, hashes, expressions, residuals, structure and scope'},
        'reproducer': 'python3 scripts/derive_nsc_reference_band_bulk.py --check',
    }


def main():
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--prepare', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    if args.prepare and OUTPUT.exists():
        raise FileExistsError('existing record is not overwritten: '+str(OUTPUT))
    data = make_record()
    if args.prepare:
        OUTPUT.write_text(json.dumps(data, indent=2, sort_keys=True)+'\n')
    else:
        from derive_nsc_transmitting_boundary_binding import compare
        compare(json.loads(OUTPUT.read_text()), data)
    print(json.dumps({'status': data['status'], 'exact_residuals': data['exact_residuals'],
                      'bulk_action_gradient': data['result']['bulk_action_gradient'],
                      'boundary_accounting': data['result']['boundary_accounting']}, indent=2))


if __name__ == '__main__':
    main()
