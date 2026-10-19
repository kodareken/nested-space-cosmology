"""Connect the existing nested slice to associative averaging and clock rates.

Only scalar summaries are printed. The existing quadratures supply all
regional moments; there is no new field grid or physical evolution.
"""
from __future__ import annotations

import json
import argparse
import hashlib
import math
from pathlib import Path
import time

import numpy as np

from recursive_horizons.nsc_atlas_geometric_average import (
    KAPPA, expansion_source_inference, partitioned_scalar_average, partition_volume_rates,
)
from recursive_horizons.nsc_atlas_nested_slice import nested_slice_record
from recursive_horizons.nsc_infinite_inheritance import (
    geometric_inherited_tree, root_spectral_moments, scale_free_response_inverse,
)


def physical_groups(regions):
    grouped = {name: np.zeros(7) for name in ('cap', 'host', 'exterior')}
    for name, moments in regions.items():
        group = 'cap' if name == 'cap' else 'host' if name.startswith('host_') else 'exterior'
        grouped[group] += np.asarray(moments)
    return grouped


def calculation():
    started = time.process_time()
    record = nested_slice_record(include_partition=True)
    part = record['partition']
    regions = part['bulk_regions_in_ell_units']
    surfaces = part['surfaces_in_ell_units']
    wall_integrals = tuple(math.fsum(wall[key] for wall in surfaces.values())
                           for key in ('curvature', 'energy', 'stress_trace'))
    fine = partitioned_scalar_average(regions, surface_integrals=wall_integrals)
    groups = physical_groups(regions)
    composed = partitioned_scalar_average(groups, surface_integrals=wall_integrals)
    coarse = partitioned_scalar_average(physical_groups(part['coarse_bulk_regions_in_ell_units']),
                                        surface_integrals=wall_integrals)
    expected = record['average_diagnostic']['scalar_outputs_in_ell_units']
    error = max(abs(fine['sources'][key]-expected[key]) for key in expected)
    regrouping = max(abs(composed['sources'][key]-fine['sources'][key]) for key in expected)
    if max(error, regrouping) > 2e-12:
        raise RuntimeError('composition disagrees with the original geometry')
    Finner, Fouter = (surfaces[name]['volume_flux'] for name in ('inner', 'outer'))
    fluxes = dict(cap=Finner, host=Fouter-Finner, exterior=-Fouter)
    rates = partition_volume_rates({name: row[0] for name, row in groups.items()},
                                    composed['regional_H'], fluxes)
    inner = record['walls']['inner']
    direct_cap_rate = 4*math.pi*inner['R_over_ell']**2*inner['Rdot']/inner['gamma']
    if abs(direct_cap_rate-rates['regional_volume_rates']['cap']) > 2e-12:
        raise RuntimeError('wall proper clock and Reynolds volume rate disagree')
    source = composed['sources']
    H = source['H_normal']
    H_rate = source['three_volume_acceleration']/3-H*H
    inferred = expansion_source_inference(H=H, H_rate=H_rate,
        mean_energy=source['kappa_mean_energy']/KAPPA,
        mean_stress_trace=3*source['kappa_mean_pressure']/KAPPA)
    if max(abs(inferred['Lambda_trace']-source['geometric_trace_projection']),
           abs(inferred['vacuum_condition']-source['geometric_vacuum_condition'])) > 2e-12:
        raise RuntimeError('independent expansion/source reconstruction disagrees')
    elapsed = time.process_time()-started
    if elapsed > 30:
        raise RuntimeError('compositional gravity application exceeded30 CPU seconds')
    return dict(
        inputs=dict(geometry='existing two-level cap/wall/Kottler/wall/Schwarzschild slice',
                    H_cap_ell=1., k_ell=2., outer_M_over_ell=.4, Rcut_over_ell=4.,
                    clock='declared initial Gaussian normal; not observed redshift time'),
        group_volumes_over_ell3={name: float(row[0]) for name, row in groups.items()},
        volume_fractions=composed['volume_fractions'],
        regional_H_times_ell=composed['regional_H'], regional_Q_times_ell2=composed['regional_Q'],
        within_region_Q_times_ell2=composed['within_region_Q'],
        between_region_Q_times_ell2=composed['between_region_Q'],
        total_Q_times_ell2=source['Q'],
        vacuum_condition_times_ell2=source['geometric_vacuum_condition'],
        geometric_Lambda_trace_times_ell2=source['geometric_trace_projection'],
        bulk_sources_and_wall_sources_separate=True,
        wall_integrals=dict(zip(('curvature_over_ell', 'energy_over_ell', 'stress_trace_over_ell'), wall_integrals)),
        boundary_fluxes_over_ell2=fluxes, initial_partition_rates=rates,
        expansion_inference=dict(H_times_ell=H, Raychaudhuri_H_rate_times_ell2=H_rate,
                                **inferred),
        controls=dict(original_moment_map_max_error=error,
                      associative_regrouping_max_error=regrouping,
                      direct_cap_volume_rate_error=direct_cap_rate-rates['regional_volume_rates']['cap'],
                      Q_composition_residual=composed['Q_composition_residual']),
        quadrature_indicators=dict(
            between_region_Q_64_minus_32=composed['between_region_Q']-coarse['between_region_Q'],
            within_region_Q_64_minus_32=composed['within_region_Q']-coarse['within_region_Q'],
            not_error_bounds=True),
        counts=dict(existing_spatial_values=record['geometric_evaluations'],
                    additional_quadrature_values=0, new_evolution_grids=0, evolution_steps=0),
        process_cpu_seconds=elapsed,
        scope='Composition and instantaneous common-normal kinematics. No parent history, Lambda drift, inherited H, or light-cone observable is inferred.')


def proof_record():
    root = Path(__file__).resolve().parents[1]
    examples = {}
    for q in (.5, .4):
        tree = geometric_inherited_tree(3, q, spatial_dimension=3)
        moments = root_spectral_moments(tree)
        coefficients = [(-1)**k*mu for k, mu in enumerate(moments['moments'][:4])]
        examples[str(q)] = dict(root_moments=moments,
            inference=scale_free_response_inverse(coefficients))
    paths = ('scripts/check_nsc_atlas_compositional_gravity.py',
             'src/recursive_horizons/nsc_atlas_geometric_average.py',
             'src/recursive_horizons/nsc_atlas_nested_slice.py',
             'src/recursive_horizons/nsc_infinite_inheritance.py',
             'tests/test_nsc_atlas_geometric_average.py',
             'tests/test_nsc_atlas_nested_slice.py',
             'tests/test_nsc_infinite_inheritance.py')
    return dict(schema='nsc-compositional-proofs-v1',
        atlas_scalar_composition=calculation(), inherited_operator_response=examples,
        source_sha256={path: hashlib.sha256((root/path).read_bytes()).hexdigest() for path in paths},
        reproduce='python scripts/run_python.py scripts/check_nsc_atlas_compositional_gravity.py --output <fresh-path>',
        identification='Two scoped mathematical bridges; the inherited tree is not identified with the cap-wall Einstein operator.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, help='optional compact summary; never overwrite evidence')
    args = parser.parse_args()
    if args.output is not None and (args.output.exists() or not args.output.parent.is_dir()):
        raise ValueError('use a fresh output in an existing directory')
    text = json.dumps(proof_record(), indent=2, allow_nan=False)+'\n'
    if len(text.encode('utf-8')) > 24000:
        raise RuntimeError('proof summary exceeds24KB')
    if args.output is None:
        print(text, end='')
    else:
        args.output.write_text(text, encoding='utf-8', newline='\n')
        print('Wrote compact proof summary:', args.output)
