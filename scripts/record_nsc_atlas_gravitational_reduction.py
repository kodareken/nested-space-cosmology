"""Record the quiet body's parent-end metric response; no field evolution."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import time

from recursive_horizons.nsc_atlas_connection_slice import wall_slice_data
from recursive_horizons.nsc_atlas_gravitational_reduction import (
    boundary_stress_from_extrinsic, exact_tube_identities,
    homogeneous_action_variations, schwarzschild_boundary_response,
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    started = time.process_time()
    wall = wall_slice_data()
    ell = wall['ell']
    mass = wall['parameters'].mass
    radius = 4*ell  # The outer cut already used by geometric_slice_average.
    response = schwarzschild_boundary_response(mass, radius)
    tensor = boundary_stress_from_extrinsic(mass, radius)
    identities = exact_tube_identities()
    action = homogeneous_action_variations()
    action_checks = {name: action[name] for name in (
        'pressure_plus_dE', 'fixed_pressure', 'extensive_vacuum_defect')}
    if any(value != 0 for value in (*identities.values(), *action_checks.values())):
        raise ArithmeticError('parent metric-response identity failed')
    root = Path(__file__).resolve().parents[1]
    paths = (
        'scripts/record_nsc_atlas_gravitational_reduction.py',
        'src/recursive_horizons/nsc_atlas_gravitational_reduction.py',
        'src/recursive_horizons/nsc_atlas_connection_slice.py',
        'src/recursive_horizons/nsc_atlas_moving_boundaries.py',
        'tests/test_nsc_atlas_gravitational_reduction.py',
    )
    record = {
        'schema': 'nsc-atlas-gravitational-reduction-v1',
        'domain': 'Quiet right Schwarzschild tube; leading dilute conserved-body response to a retained parent metric.',
        'reference': {'mass_over_ell': mass/ell, 'tube_radius_over_ell': radius/ell},
        'tube': {
            'proper_lapse_factor': response['lapse_factor'],
            'energy_over_ell': response['quasilocal_energy']/ell,
            'surface_energy_density_times_ell': response['surface_energy_density']*ell,
            'surface_pressure_times_ell': response['surface_pressure']*ell,
            'adm_mass_over_ell': response['adm_mass']/ell,
            'extrinsic_tensor_residuals_times_ell': {
                name: tensor[name]*ell for name in (
                    'time_time_plus_surface_density', 'angular_minus_surface_pressure')},
        },
        'leading_parent_source': {
            'action': '-sum_i M_i integral d_tau_i, in G=c=1 units',
            'conserved_homogeneous_energy': 'sum_i M_i, independent of physical volume',
            'rho_geom': 'n*M for identical bodies',
            'p_geom': '0',
            'physical_number_density_n_ell3': None,
            'lambda_trace_ell2_coefficient_of_n_ell3': 2*math.pi*mass/ell,
            'vacuum_form_defect_ell2_coefficient_of_n_ell3': 8*math.pi*mass/ell,
            'generated_volume_lambda_ell2_at_this_order': 0,
            'population_or_common_parent_constructed': False,
        },
        'vacuum_energy_condition': 'dE/dV=E/V, hence E=rho_Lambda*V with constant rho_Lambda for a conserved homogeneous sector',
        'exact_zero_checks': {
            'tube': {name: str(value) for name, value in identities.items()},
            'homogeneous_action': {name: str(value) for name, value in action_checks.items()},
        },
        'uncertainty_scope': 'Exact leading identities and floating reference evaluation; no total physical or observational error bound.',
        'unresolved': [
            'Compatible parent geometry, physical population and parent observer domain',
            'Finite-size, tidal, mutual-interaction and velocity corrections',
            'Further physical parent-metric functional, state and conserved exchange needed for a nonzero vacuum term',
        ],
        'source': 'https://arxiv.org/abs/gr-qc/9209012',
        'reading': 'Section III equations 3.14-3.19; Section VI equations 6.2-6.20; independently checked tube curvatures, derivatives and homogeneous action',
        'source_sha256': {name: hashlib.sha256((root/name).read_bytes()).hexdigest() for name in paths},
        'reproduce': 'python scripts/run_python.py scripts/record_nsc_atlas_gravitational_reduction.py --output .build/nsc-atlas-gravitational-reduction-reproduced.json',
        'physical_time_steps': 0,
        'cpu_seconds_limit': 30,
    }
    if time.process_time()-started > 30:
        raise RuntimeError('parent metric-response record exceeded 30 CPU seconds')
    args.output.write_text(json.dumps(record, indent=2, sort_keys=True, allow_nan=False)+'\n',
                           encoding='utf-8', newline='\n')
    print(f'Wrote {args.output}; {time.process_time()-started:.3f} CPU seconds; no evolution.')


if __name__ == '__main__':
    main()
