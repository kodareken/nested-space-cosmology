#!/usr/bin/env python3
"""Record the canonical/measure allocation using authenticated existing data."""
import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from recursive_horizons.nsc_canonical_spectral_bridge import allocate
from check_nsc_compact_boundary_action import authenticated_record
from check_nsc_compact_casimir import compare
from check_nsc_vacuum_charge_matching import hashes

OUTPUT = ROOT/'results/development/canonical-spectral-bridge.json'
INPUTS = ('results/development/compact-matching.json',)
SOURCES = ('src/recursive_horizons/nsc_canonical_spectral_bridge.py',
           'scripts/check_nsc_canonical_spectral_bridge.py',
           'docs/nsc-canonical-spectral-bridge.md',
           'docs/nsc-compact-mass-map.md', 'docs/nsc-warped-source.md',
           'src/recursive_horizons/nsc_compact_matching.py',
           'scripts/check_nsc_compact_boundary_action.py',
           'scripts/check_nsc_compact_casimir.py',
           'scripts/check_nsc_vacuum_charge_matching.py')


def calculate():
    result = allocate(authenticated_record(INPUTS[0]))
    for row in result['rows']:
        # Binding to the stored coefficients, not a rerun of their generator.
        assert max(abs(x) for x in row['moment_reconstruction_residuals'].values()) < 5e-11
        assert max(abs(x) for x in row['coefficient_reconstruction_residuals'].values()) < 5e-13
    return {'schema': 'NSC-CANONICAL-SPECTRAL-BRIDGE-v1',
            'status': 'canonical tower and covariant-conversion origins of the stored Euclidean source are evaluated',
            'source_hashes': hashes(SOURCES), 'input_hashes': hashes(INPUTS),
            'allocation': result,
            'scope': {'new_spectral_or_source_integrals': False,
                      'new_independent_coefficients': False,
                      'physical_canonical_field_content_identified': True,
                      'full_real_time_covariant_conversion_derived': False,
                      'absolute_self_sourced_solution_derived': False,
                      'local_a4_truncation_valid_at_neck_established': False},
            'comparison': {'fields': 'all', 'float_atol': 3e-9, 'float_rtol': 3e-8,
                           'exact': 'structure, non-float values, source and input hashes',
                           'exceptions': []}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--output', type=Path, default=OUTPUT)
    args = parser.parse_args()
    result = calculate()
    if args.check:
        compare(json.loads(args.output.read_text()), result)
        print('Canonical/spectral allocation reproduced: all fields agree; no old generator ran.')
    else:
        # Refuse to overwrite a scientific record.
        with args.output.open('x') as handle:
            json.dump(result, handle, indent=2, sort_keys=True, allow_nan=False)
            handle.write('\n')
        print(f'Wrote {args.output}')


if __name__ == '__main__':
    main()
