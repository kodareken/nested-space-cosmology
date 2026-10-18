"""Frozen finite Dirac mean-metric/radar-clock reference; no P001 evolution."""
from __future__ import annotations

import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['NUMBA_NUM_THREADS'] = '1'

import argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
import hashlib
import json
import multiprocessing
from pathlib import Path
import subprocess
import time

import numpy as np
import scipy
from scipy import constants
from threadpoolctl import threadpool_limits

from recursive_horizons.nsc_atlas_dirac_radar import WallRadarChart
from recursive_horizons.nsc_atlas_dirac_mean import PreparedDiracSource
from recursive_horizons.nsc_atlas_dirac_backreaction import MeanBackreaction

ROOT = Path(__file__).resolve().parents[1]
SOURCES = (
    'src/recursive_horizons/nsc_atlas_dirac_mean.py',
    'src/recursive_horizons/nsc_atlas_dirac_backreaction.py',
    'src/recursive_horizons/nsc_atlas_dirac_gravity.py',
    'src/recursive_horizons/nsc_atlas_dirac_stress.py',
    'src/recursive_horizons/nsc_atlas_dirac_scattering.py',
    'src/recursive_horizons/nsc_atlas_dirac_radar.py',
    'src/recursive_horizons/nsc_atlas_dirac_evolution.py',
    'src/recursive_horizons/nsc_atlas_cap_wall_response.py',
    'src/recursive_horizons/nsc_atlas_causal_reference.py',
    'tests/test_nsc_atlas_dirac_mean.py',
    'tests/test_nsc_atlas_dirac_backreaction.py',
    'scripts/measure_nsc_atlas_dirac_mean_clock.py',
    'pyproject.toml',
)


def prepare_pair(cells):
    with threadpool_limits(limits=1):
        started = time.process_time()
        def remaining():
            value = 20-(time.process_time()-started)
            if value<=0:
                raise RuntimeError('mean-clock setup exhausted its 20 CPU seconds')
            return value
        chart = WallRadarChart()
        U,V = [float(chart.maps(r)[0]) for r in (.475,.525)]
        nested_source = PreparedDiracSource(chart,cells)
        nested = MeanBackreaction(nested_source,V,wall_cpu_limit=remaining())
        cap_source = PreparedDiracSource(chart,cells,region='cap')
        cap = MeanBackreaction(cap_source,V,shared_worldtube=nested,wall_cpu_limit=remaining())
        metadata = dict(cells=cells,wall_clock_duration=chart.duration,
            receiver_U=U,receiver_V=V,receiver_radar_range=(V-U)/2,
            nested_input_flux=nested_source.grid.incoming_flux,
            cap_input_flux=cap_source.grid.incoming_flux,
            nested_flux_defect=nested_source.grid.flux_error,
            cap_flux_defect=cap_source.grid.flux_error,
            packet_norm_quadrature_error=nested_source.packet.norm_quadrature_error,
            nested_wall_rhs=nested.wall_nfev,cap_wall_rhs=cap.wall_nfev)
        state = nested.wall_state(V)
        metadata['nested_wall_endpoint'] = {k:float(state[k]) for k in
            ('radius','displacement','displacement_rate','exterior_mass','cap_mass')}
        return dict(cells=cells,nested=nested,cap=cap,metadata=metadata,
                    pid=os.getpid(),process_CPU_total=time.process_time())


def evaluate_clock(cells,geometry,model,order,U,V):
    with threadpool_limits(limits=1):
        result = model.clock_response(U,V,order=order,cpu_limit=20)
        result.update(cells=cells,geometry=geometry,quadrature_order=order)
        return dict(result=result,pid=os.getpid(),process_CPU_total=time.process_time())


def compare(results,metadata):
    issues = []
    rows = []
    for data in metadata:
        for geometry in ('nested','cap'):
            norm = data[geometry+'_input_flux']
            if abs(norm-1)>1e-5:
                issues.append(f'{geometry} N={data["cells"]}: unresolved input-flux quadrature')
            if data[geometry+'_flux_defect']>32*data['cells']*np.finfo(float).eps*max(1.,norm):
                issues.append(f'{geometry} N={data["cells"]}: propagation flux defect')
    for geometry in ('nested','cap'):
        by_level = {}
        for cells in (64,128):
            lo,hi = [next(r for r in results if r['cells']==cells and
                r['geometry']==geometry and r['quadrature_order']==q) for q in (8,12)]
            quadrature_spread = abs(hi['frequency_fractional_change']-lo['frequency_fractional_change'])
            normalized_mass_gap = abs(hi['mass_path_gap'])/max(1.,abs(hi['mass']),abs(hi['mass_other_path']))
            if quadrature_spread > 2e-6*max(1.,abs(hi['frequency_fractional_change'])):
                issues.append(f'{geometry} N={cells}: unresolved quadrature')
            if normalized_mass_gap>1e-3:
                issues.append(f'{geometry} N={cells}: unresolved mass-path constraint')
            if abs(hi['wall_clock_log_C_at_corner'])>2e-9:
                issues.append(f'{geometry} N={cells}: inconsistent proper-wall clock')
            by_level[cells] = dict(clock_coefficient=hi['frequency_fractional_change'],
                normalized_mass_path_gap=normalized_mass_gap,
                quadrature_spread=quadrature_spread)
        if (by_level[128]['normalized_mass_path_gap']>1e-7 and
                by_level[128]['normalized_mass_path_gap']>by_level[64]['normalized_mass_path_gap']/1.3):
            issues.append(geometry+': mass-path discrepancy did not contract')
        rows.append(dict(geometry=geometry,levels=by_level,
            absolute_64_128_clock_spread=abs(by_level[128]['clock_coefficient']-
                                            by_level[64]['clock_coefficient'])))
    fine = {r['geometry']:r for r in results if r['cells']==128 and r['quadrature_order']==12}
    coarse = {r['geometry']:r for r in results if r['cells']==64 and r['quadrature_order']==12}
    delta_fine = fine['nested']['frequency_fractional_change']-fine['cap']['frequency_fractional_change']
    delta_coarse = coarse['nested']['frequency_fractional_change']-coarse['cap']['frequency_fractional_change']
    return rows,dict(fine_clock_coefficient_difference=delta_fine,
        absolute_64_128_difference_spread=abs(delta_fine-delta_coarse)),issues


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    if args.output.exists() or not args.output.parent.is_dir():
        raise ValueError('use an unused compact-summary path with an existing parent')
    commit = subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    bindings = {}
    for name in SOURCES:
        raw = (ROOT/name).read_bytes()
        if subprocess.check_output(['git','show',f'{commit}:{name}'],cwd=ROOT) != raw:
            raise ValueError('freeze the mean-clock method before its reference: '+name)
        bindings[name] = hashlib.sha256(raw).hexdigest()
    started = time.perf_counter()
    workers = min(8,os.cpu_count() or 1)
    worker_totals = {}
    results,metadata = [],[]
    budget_failure = None
    pool = ProcessPoolExecutor(max_workers=workers,mp_context=multiprocessing.get_context('spawn'))
    try:
        preparations = {pool.submit(prepare_pair,n):n for n in (64,128)}
        observations = []
        for future in as_completed(preparations):
            prepared = future.result()
            worker_totals[prepared['pid']] = prepared['process_CPU_total']
            data = prepared['metadata']
            metadata.append(data)
            for geometry in ('nested','cap'):
                for order in (8,12):
                    observations.append(pool.submit(evaluate_clock,prepared['cells'],geometry,
                        prepared[geometry],order,data['receiver_U'],data['receiver_V']))
            print(json.dumps(dict(prepared_cells=prepared['cells'],
                wall_rhs=data['nested_wall_rhs']+data['cap_wall_rhs'])),flush=True)
        for future in as_completed(observations):
            value = future.result()
            worker_totals[value['pid']] = max(worker_totals.get(value['pid'],0.),
                                             value['process_CPU_total'])
            row = value['result'];results.append(row)
            charged_cpu = time.process_time()+sum(worker_totals.values())+8
            if charged_cpu>240:
                budget_failure='240 process-CPU-second shared allowance exhausted'
            print(json.dumps(dict(cells=row['cells'],geometry=row['geometry'],
                quadrature=row['quadrature_order'],clock=row['frequency_fractional_change'],
                mass_path_gap=row['mass_path_gap'],charged_CPU_seconds=charged_cpu)),flush=True)
    except Exception as error:
        if hasattr(pool,'terminate_workers'):
            pool.terminate_workers()
        failed = dict(schema='nsc-atlas-dirac-mean-clock-v1',status='execution_exit_recorded',
            source_commit=commit,source_sha256=bindings,error=str(error),
            completed_results=results,metadata=metadata,
            scope='Stopped before a complete mean-clock result; no new P001 grids.',
            runtime=dict(wall_seconds=time.perf_counter()-started,
                completed_worker_CPU_seconds=sum(worker_totals.values()),
                interrupted_worker_CPU_unavailable=True))
        args.output.write_text(json.dumps(failed,indent=2,sort_keys=True)+'\n',
                               encoding='utf-8',newline='\n')
        raise SystemExit(str(error))
    finally:
        pool.shutdown(wait=True,cancel_futures=True)
    results.sort(key=lambda r:(r['cells'],r['geometry'],r['quadrature_order']))
    metadata.sort(key=lambda r:r['cells'])
    comparison,difference,issues = compare(results,metadata)
    if budget_failure:
        issues.append(budget_failure)
    G,Gunit,Gerror = constants.physical_constants['Newtonian constant of gravitation']
    ell = 10000.
    coupling = 8*constants.pi*constants.hbar*G/(constants.c**3*ell**2)
    physical = {r['geometry']:dict(
        frequency_fractional_change=coupling*r['frequency_fractional_change'],
        absolute_frequency_ratio_change=coupling*r['received_over_wall_frequency_change'],
        observer_areal_radius_change_m=ell*coupling*r['radius_change'],
        G_standard_uncertainty_fractional=Gerror/G)
        for r in results if r['cells']==128 and r['quadrature_order']==12}
    record = dict(schema='nsc-atlas-dirac-mean-clock-v1',
        status='unresolved_mean_clock_reference' if issues else 'finite_prepared_mean_clock_evaluated',
        source_commit=commit,source_sha256=bindings,validation_issues=issues,
        model=dict(H_ell=1.,k_wall_ell=2.,M_over_ell=.4,angular_k=1,degeneracy=2,
            field='minimal massless Dirac probe, continuum-normalized finite occupation difference',
            gravity='original Einstein-Hilbert plus fixed regional vacuum terms and material tension',
            order='first order in 8*pi*(l_P/ell)^2 times occupation',
            root_vacuum_selected=False,vacuum_angular_sum=False),
        preparation=dict(incoming='eta_minus on U=0, eta_plus=0',
            envelope_center_over_duration=.6,half_width_over_duration=.3,
            carrier_phase='exp(-12i*V/duration)',continuum_flux=1.,
            covariance_difference='nu*|f><f| after emptying that mode in the supplied parent',
            energy_generator='per multiplet <i*dV> = 2*12/duration in hbar=c=ell=1'),
        observer=dict(emission_wall_radius_label=.475,reception_wall_radius_label=.525,
            phase='-U',worldline='fixed proper-wall radar position',
            cap_control='smooth cap with the nested wall areal-radius history explicitly supplied; no membrane'),
        metadata=metadata,results=results,comparison=comparison,nesting_difference=difference,
        physical_scale=dict(ell_m=ell,origin='independently supplied nominal scale, not a fitted sky value',
            hbar_J_s=constants.hbar,c_m_s=constants.c,G=G,G_unit=Gunit,G_standard_uncertainty=Gerror,
            constants_source='https://physics.nist.gov/cuu/Constants/Table/allascii.txt',
            constant_set='CODATA 2022 in SciPy '+scipy.__version__,quantum_gravitational_factor=coupling,
            wall_clock_duration_s=metadata[0]['wall_clock_duration']*ell/constants.c,
            radar_range_m=metadata[0]['receiver_radar_range']*ell,
            multiplet_incoming_clock_generator_energy_J=2*constants.hbar*constants.c*12/(
                metadata[0]['wall_clock_duration']*ell),predictions=physical),
        uncertainty='64/128 and 8/12 quadrature spreads are numerical diagnostics; '
            'G uncertainty is separate. No proven continuum bound or observational confidence interval.',
        runtime=dict(wall_seconds=time.perf_counter()-started,
            observed_process_CPU_seconds=time.process_time()+sum(worker_totals.values()),
            unobserved_serialization_overhead_CPU_charge=8,
            charged_CPU_seconds=time.process_time()+sum(worker_totals.values())+8,
            process_workers=workers,BLAS_threads_per_worker=1,Numba_threads_per_worker=1,
            new_P001_grids=0,new_P001_time_steps=0,baseline_vacuum_polarizations=0),
        limits=dict(charged_CPU_seconds=240,summary_bytes=24000),
        nested_Psi_over_Phi=None,
        scope='Mean metric and finite radar-clock response to one specified prepared multiplet; '
              'the mixed retarded source/tadpole response and full baseline-vacuum Psi/Phi are not computed.')
    raw = json.dumps(record,indent=2,sort_keys=True,allow_nan=False)+'\n'
    if len(raw.encode())>24000:
        raise RuntimeError('compact mean-clock summary exceeds its admitted size')
    args.output.write_text(raw,encoding='utf-8',newline='\n')
    print(json.dumps(dict(status=record['status'],validation_issues=issues,comparison=comparison,
        nesting_difference=difference,physical=physical,charged_CPU_seconds=record['runtime']['charged_CPU_seconds'])),flush=True)
    if issues:
        raise SystemExit(1)


if __name__=='__main__':
    main()
