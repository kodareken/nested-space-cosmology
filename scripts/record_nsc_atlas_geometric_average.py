"""Save the fixed-domain geometric scalar average and quadrature differences."""
import argparse
import hashlib
import json
from pathlib import Path
import time

from recursive_horizons.nsc_atlas_geometric_average import geometric_slice_average


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    started = time.process_time()
    rows = {scheme: geometric_slice_average(quadrature=scheme)
            for scheme in ('gauss32', 'gauss64', 'adaptive')}
    reference = rows['adaptive']
    keys = tuple(reference['normalized_sources'])
    moments = reference['integrals_in_ell_units']
    I = moments['theta']
    J = 2*moments['theta_squared']/3-2*moments['shear_squared']
    IR = moments['curvature_total']
    root = Path(__file__).resolve().parents[1]
    paths = ('src/recursive_horizons/nsc_atlas_geometric_average.py',
             'src/recursive_horizons/nsc_atlas_connection_slice.py',
             'src/recursive_horizons/nsc_atlas_moving_boundaries.py',
             'scripts/record_nsc_atlas_geometric_average.py',
             'tests/test_nsc_atlas_geometric_average.py')
    record = {
        'schema': 'nsc-atlas-geometric-average-summary-v1',
        'reference': reference,
        'maximum_absolute_quadrature_difference': {
            key: max(abs(row['normalized_sources'][key]-reference['normalized_sources'][key])
                     for row in rows.values()) for key in keys},
        'maximum_volume_quadrature_difference_over_ell3': max(
            abs(row['volume_over_ell3']-reference['volume_over_ell3']) for row in rows.values()),
        'time_symmetric_vacuum_extension_condition': {
            'Itheta_over_ell2': I, 'J_over_ell': J, 'IR_over_ell': IR,
            'required_volume_per_body_over_ell3': 2*I**2/(3*J+IR),
            'minimum_retained_volume_per_body_over_ell3': reference['volume_over_ell3'],
            'scope': 'Conditional algebra for disjoint copies and a K=0 vacuum connector; no such connector metric is evaluated here.',
        },
        'uncertainty_scope': 'Quadrature differences only, not a certified bound or a cosmological/model uncertainty.',
        'identification': 'The scalar geometric correction fails p=-rho on this slice. Its trace is not a measured Lambda.',
        'unresolved': ['Physical common-parent averaging domain and clock',
                       'Common-parent metric, evolution and population data',
                       'Tensor vacuum form, conservation and observational comparison'],
        'source_sha256': {name: hashlib.sha256((root/name).read_bytes()).hexdigest() for name in paths},
        'reproduce': 'python scripts/run_python.py scripts/record_nsc_atlas_geometric_average.py --output results/summaries/nsc-atlas-geometric-average-v1.json',
    }
    if time.process_time()-started > 30:
        raise RuntimeError('summary exceeded its combined 30 CPU-second budget')
    args.output.write_text(json.dumps(record, indent=2, sort_keys=True, allow_nan=False)+'\n',
                          encoding='utf-8', newline='\n')
    print(f'Wrote {args.output}; {time.process_time()-started:.3f} CPU seconds; no evolution.')


if __name__ == '__main__':
    main()
