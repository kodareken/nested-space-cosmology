"""Invert the evolved lens for a fixed receiver and check its delay clock."""
import argparse
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict, replace
import hashlib
import json
from pathlib import Path
import time

from recursive_horizons.nsc_atlas_core_optics import OpticalRun
from recursive_horizons.nsc_atlas_core_images import run_images, reduce_run, EPOCH, OFFSET
from recursive_horizons.nsc_atlas_core_reservoir import CoreReservoir
from check_nsc_atlas_core_optics import SOURCE_INPUTS

ROOT = Path(__file__).resolve().parents[1]
INPUTS = SOURCE_INPUTS+('src/recursive_horizons/nsc_atlas_core_images.py',)


def worker(item):
    name, config = item
    key = hashlib.sha256(json.dumps(dict(config=asdict(config),
        source_sha256={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in INPUTS}),
        sort_keys=True).encode()).hexdigest()
    folder = ROOT/'.build/core-images'/key; path = folder/'completed.json'
    record = json.loads(path.read_text(encoding='utf-8')) if path.is_file() else run_images(config, cache_dir=folder)
    return name, record


def calculation():
    wall = time.perf_counter(); cpu = time.process_time(); base = OpticalRun()
    configs = dict(coarse=replace(base, spacing=.001, step=.0001), medium=base,
        fine=replace(base, spacing=.00025, step=.000025), optical_time_control=replace(base, block_steps=2))
    with ProcessPoolExecutor(max_workers=4) as pool:
        raw = dict(pool.map(worker, configs.items()))
    reduced = {name:reduce_run(r) for name,r in raw.items()}
    old = json.loads((ROOT/'results/summaries/nsc-atlas-core-evolution-v1.json').read_text(encoding='utf-8'))
    expected = dict(coarse='diffuse_coarse', medium='diffuse_medium', fine='diffuse_fine', optical_time_control='diffuse_medium')
    matches = {name:r['final_field_sha256']==old['runs'][expected[name]]['final_state_sha256'] for name,r in raw.items()}
    centre = {name:r['pairs'][2] for name,r in reduced.items()}
    values = centre['fine']; envelopes = {}
    for key,value in values.items():
        if key=='emission_proper_age': continue
        envelopes[key] = 2*max(abs(value-centre['medium'][key]),
            abs(centre['medium'][key]-centre['optical_time_control'][key]))
    # At the fixed target, curvature and the differentiated landing map
    # give independently calculated beam areas. Include their discrepancy
    # in the ratio uncertainty rather than relying only on grid differences.
    images = [x for x in reduced['fine']['images'] if x['emission_proper_age']==EPOCH]
    alt_areas = [(im['jacobi_perpendicular']-im['rotation_jacobi_gap'])*
        (im['jacobi_parallel']-im['parallel_jacobi_gap']) for im in images]
    for key,power in (('fluence_ratio_negative_over_positive',1), ('luminosity_flux_ratio_negative_over_positive',2)):
        alt=(images[1]['received_frequency']/images[0]['received_frequency'])**power*abs(alt_areas[0]/alt_areas[1])
        envelopes[key]=max(envelopes[key],2*abs(alt-values[key]))
    checks = {}
    for name,r in reduced.items():
        checks[name] = {key:max(abs(im[key]) for im in r['images']) for key in
            ('independent_root_gap', 'landing_residual', 'endpoint_variation_gap',
             'rotation_jacobi_gap', 'parallel_jacobi_gap', 'frequency_transport_gap')}
        checks[name]['maximum_nested_root_gap'] = max(abs(im['nested_angular_gaps']['launch_angle']) for im in r['images'])
        checks[name]['delay_drift'] = r['delay_drift_check']
    drift_fine=checks['fine']['delay_drift']; drift_medium=checks['medium']['delay_drift']
    drift_control=checks['optical_time_control']['delay_drift']
    drift_fine['derivative_comparison_envelope']=2*max(
        abs(drift_fine['delay_derivative']-drift_medium['delay_derivative']),
        abs(drift_control['delay_derivative']-drift_medium['delay_derivative']),
        abs(drift_fine['identity_gap']), abs(drift_fine['derivative_refinement_gap']))
    errors = []
    if not all(matches.values()): errors.append('source replay hash differs from frozen evolution')
    rotation_convergence=[checks[name]['rotation_jacobi_gap'] for name in ('coarse','medium','fine')]
    if not rotation_convergence[0]>rotation_convergence[1]>rotation_convergence[2]:
        errors.append('rotation/Jacobi error does not decrease under refinement')
    # Acceptance is on independent identities, not agreement with a desired lens signal.
    if checks['fine']['independent_root_gap']>2e-7: errors.append('curvature/root mismatch')
    if checks['fine']['maximum_nested_root_gap']>2e-8: errors.append('angular interpolation unresolved')
    if checks['fine']['endpoint_variation_gap']>2e-7: errors.append('endpoint/Fermat check failed')
    if checks['fine']['rotation_jacobi_gap']>4e-6: errors.append('rotation/Jacobi check failed')
    if checks['fine']['parallel_jacobi_gap']>8e-6: errors.append('parallel/Jacobi check failed')
    if checks['fine']['frequency_transport_gap']>1e-7: errors.append('frequency transport failed')
    if abs(checks['fine']['delay_drift']['identity_gap'])>1e-6: errors.append('fixed-observer delay/redshift identity failed')
    if values['delay_negative_minus_positive']<=envelopes['delay_negative_minus_positive']:
        errors.append('image arrival order unresolved')
    if images[0]['parity']!=1 or images[1]['parity']!=-1:
        errors.append('selected opposite image parities unresolved')
    total_cpu = time.process_time()-cpu+sum(r['cpu_seconds'] for r in raw.values())
    parent=CoreReservoir().reference
    vacuum=parent.scalar.kappa*parent.scalar.vacuum(parent.parent_flux)
    Lref=parent.cap_geometry['parent_maximum_radius']
    delay=values['delay_negative_minus_positive']; delay_error=envelopes['delay_negative_minus_positive']
    scale_connection=dict(vacuum_Lambda_in_atlas_units=vacuum, reference_length_in_atlas_units=Lref,
        vacuum_Lambda_reference_length_squared=vacuum*Lref**2,
        Lambda_c_delay_squared=vacuum*delay**2,
        numerical_envelope=vacuum*(2*abs(delay)*delay_error+delay_error**2),
        delay_over_reference_length=delay/Lref,
        receiver_clock_delay_drift=1/values['frequency_ratio_negative_over_positive']-1,
        interpretation='Lambda is the specified inherited parent vacuum, not a newly selected observed cosmological constant. The same dimensional embedding implies Lambda_phys*(c*delay_seconds)^2 equals this number. A delay-based scale calibration leaves the parent vacuum, angles and ratios as linked tests; they cannot be independently rescaled.')
    files = INPUTS+('scripts/check_nsc_atlas_core_images.py', 'scripts/plot_nsc_atlas_core_images.py',
        'tests/test_nsc_atlas_core_images.py', 'results/summaries/nsc-atlas-core-evolution-v1.json')
    record = dict(schema='nsc-atlas-core-images-v1',
        result='Two resolved images of the same emitted pulse on one held receiver in the unchanged evolved scalar parent.',
        protocol=dict(screen_radius=.6, source_proper_age=EPOCH, receiver_azimuth='pi - offset',
            offset=OFFSET, opposite_orientation_target='pi + offset',
            claim='Two resolved branches in declared patches; not an exhaustive image census.'),
        central_pair=values, numerical_comparison_envelope=envelopes, parent_vacuum_connection=scale_connection,
        central_images=[x for x in reduced['fine']['images'] if x['emission_proper_age']==EPOCH],
        refinement={name:r['pairs'][2] for name,r in reduced.items()}, independent_checks=checks,
        source_replay_exact_matches=matches,
        runs={name:dict(config=r['config'], final_field_sha256=r['final_field_sha256'], cpu_seconds=r['cpu_seconds']) for name,r in raw.items()},
        identities=['d tau_o = d tau_e / nu_o + (L/nu_o) d phi_o',
            'd(tau_minus-tau_plus)/d tau_e = 1/nu_minus - 1/nu_plus',
            'sin xi = sign(L) sin psi / nu_o',
            'fluence per emitted isotropic energy = nu_o/(4*pi*abs(D_perp*D_parallel))'],
        limits='Conditional supplied cloud, source pulse, held source/receiver and geometric optics. Two angular patches do not exclude other images. No astronomical object identification, physical scale, population selection, finite-source smoothing or nonspherical stability is inferred. Numerical envelopes are comparison estimates, not interval proofs or preparation uncertainties.',
        numerical_validation_errors=errors, process_cpu_seconds_including_workers=total_cpu,
        wall_seconds=time.perf_counter()-wall,
        longest_worker_original_wall_seconds=max(r['wall_seconds'] for r in raw.values()),
        compute_budget_seconds=4500, budget_exhausted=total_cpu>4500,
        timing_scope='Worker CPU includes original cached computations; wall time is this invocation.',
        source_sha256={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in files},
        reproduce='python scripts/run_python.py scripts/check_nsc_atlas_core_images.py --output results/summaries/nsc-atlas-core-images-v1.json')
    return record, raw, reduced


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output', type=Path)
    p.add_argument('--figure',action='store_true');args=p.parse_args()
    record,raw,reduced=calculation(); encoded=json.dumps(record, indent=2, allow_nan=False)+'\n'
    if args.output:args.output.write_text(encoded, encoding='utf-8', newline='\n')
    else:print(encoded)
    print(json.dumps(dict(pair=record['central_pair'], errors=record['numerical_validation_errors'],
        cpu_seconds=record['process_cpu_seconds_including_workers'],wall_seconds=record['wall_seconds'])))
    if record['numerical_validation_errors']:raise RuntimeError(record['numerical_validation_errors'])
    if record['budget_exhausted']:raise RuntimeError('4500 CPU second budget exhausted; diagnose before another campaign')
    if args.figure:
        from plot_nsc_atlas_core_images import render
        render(record,raw,reduced)
