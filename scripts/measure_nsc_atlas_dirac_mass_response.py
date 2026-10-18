"""Prepared Dirac mean response to a conserved compact mass on the nested cap."""
from __future__ import annotations

import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
os.environ['NUMBA_NUM_THREADS']='1'

import argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
import hashlib
import json
import multiprocessing
from pathlib import Path
import subprocess
import time

from scipy import constants
from threadpoolctl import threadpool_limits

from recursive_horizons.nsc_atlas_dirac_radar import WallRadarChart
from recursive_horizons.nsc_atlas_dirac_compact_source import CompactMassChart
from recursive_horizons.nsc_atlas_dirac_mean import PreparedDiracSource
from recursive_horizons.nsc_atlas_dirac_backreaction import MeanBackreaction
from recursive_horizons.nsc_atlas_dirac_tidal import mean_tidal_response
from measure_nsc_atlas_dirac_mean_clock import SOURCES as MEAN_SOURCES

ROOT=Path(__file__).resolve().parents[1]
SOURCES=tuple(dict.fromkeys((*MEAN_SOURCES,
    'src/recursive_horizons/nsc_atlas_dirac_compact_source.py',
    'src/recursive_horizons/nsc_atlas_dirac_tidal.py',
    'tests/test_nsc_atlas_dirac_compact_source.py',
    'tests/test_nsc_atlas_dirac_tidal.py',
    'scripts/measure_nsc_atlas_dirac_mass_response.py')))
MASSES=(0.,5e-5,1e-4)
CPU_LIMIT=480
DIAGNOSIS='results/summaries/nsc-atlas-dirac-source-precision-v1.json'


def prepare(cells,mass):
    with threadpool_limits(limits=1):
        started=time.process_time()
        chart=CompactMassChart(WallRadarChart(),mass)
        U,V=.8*chart.duration,.7*chart.duration
        def remaining():
            value=30-(time.process_time()-started)
            if value<=0:
                raise RuntimeError('compact-source setup exhausted 30 CPU seconds')
            return value
        nested_source=PreparedDiracSource(chart,cells)
        nested=MeanBackreaction(nested_source,U,wall_rtol=2e-11,wall_atol=2e-13,
            wall_cpu_limit=remaining())
        cap_source=PreparedDiracSource(chart,cells,region='cap')
        cap=MeanBackreaction(cap_source,U,shared_worldtube=nested,
            wall_rtol=2e-11,wall_atol=2e-13,wall_cpu_limit=remaining())
        metadata=dict(cells=cells,central_mass=mass,parent_mass=chart.parent_mass,
            duration=chart.duration,wall_end_radius=chart.wall_end_radius,
            receiver_U=U,receiver_V=V,
            incoming_flux=nested_source.grid.incoming_flux,
            flux_defect=max(nested_source.grid.flux_error,cap_source.grid.flux_error),
            actual_wall_rhs=nested.wall_nfev+cap.wall_nfev)
        return dict(models=dict(nested=nested,cap=cap),metadata=metadata,
            pid=os.getpid(),process_CPU_total=time.process_time())


def observe(model,cells,mass,geometry,order,U,V):
    with threadpool_limits(limits=1):
        result=model.clock_response(U,V,order=order,cpu_limit=30)
        source=model.source.evaluate(U,V)
        fields,stress=source['fields'],source['orthonormal']
        radius=float(fields['radius'])
        tidal=mean_tidal_response(radius=radius,background_mass=mass,
            mass_change=result['mass'],radius_change=result['radius_change'],stress=stress)
        tidal={k:None if v is None else float(v) for k,v in tidal.items()}
        mass_scale=max(1.,abs(result['mass']),abs(result['mass_other_path']))
        row=dict(cells=cells,central_mass=mass,geometry=geometry,order=order,
            clock=result['frequency_fractional_change'],radius_change=result['radius_change'],
            mass_change=result['mass'],normalized_mass_path_gap=abs(result['mass_path_gap'])/mass_scale,
            wall_clock_defect=abs(result['wall_clock_log_C_at_corner']),
            receiver_radius=radius,anisotropic_pressure=float(stress['radial_pressure']-stress['angular_pressure']),
            tidal=tidal)
        return dict(result=row,pid=os.getpid(),process_CPU_total=time.process_time())


def comparisons(rows,metadata):
    issues=[]
    table={(r['cells'],r['central_mass'],r['geometry'],r['order']):r for r in rows}
    for data in metadata:
        if abs(data['incoming_flux']-1)>1e-5 or data['flux_defect']>2e-12:
            issues.append('Unresolved normalization/charge flux: '+str(data))
    for r in rows:
        if r['normalized_mass_path_gap']>1e-3 or r['wall_clock_defect']>2e-9:
            issues.append('Unresolved mass or proper-wall clock: '+str((r['cells'],r['central_mass'],r['geometry'])))
    summary={}
    for geometry in ('nested','cap'):
        by_level={}
        for cells in (64,128):
            fine=[table[cells,b,geometry,12] for b in MASSES]
            derivatives=[(fine[i]['clock']-fine[0]['clock'])/MASSES[i] for i in (1,2)]
            gaps=[r['normalized_mass_path_gap'] for r in fine]
            qdiff=max(abs(table[cells,b,geometry,8]['clock']-table[cells,b,geometry,12]['clock'])
                for b in MASSES)
            if qdiff>2e-6*max(1.,*(abs(r['clock']) for r in fine)):
                issues.append(f'{geometry} N={cells}: unresolved clock quadrature')
            by_level[cells]=dict(clock_zero_mass=fine[0]['clock'],
                clock_mass_derivative_h_half=derivatives[0],clock_mass_derivative_h=derivatives[1],
                extrapolated_clock_mass_derivative=2*derivatives[0]-derivatives[1],
                derivative_step_spread=abs(derivatives[0]-derivatives[1]),
                max_normalized_mass_path_gap=max(gaps),max_clock_quadrature_spread=qdiff,
                tidal_ratio_coefficient_at_b_1e_4=fine[-1]['tidal']['ratio_change'])
            if by_level[cells]['derivative_step_spread']>.02*max(
                    1.,abs(by_level[cells]['extrapolated_clock_mass_derivative'])):
                issues.append(f'{geometry} N={cells}: unresolved source derivative step')
        if (by_level[128]['max_normalized_mass_path_gap']>1e-7 and
            by_level[128]['max_normalized_mass_path_gap']>
                by_level[64]['max_normalized_mass_path_gap']/1.3):
            issues.append(geometry+': mass-path gap did not contract')
        summary[geometry]=dict(levels=by_level,
            derivative_grid_spread=abs(by_level[128]['extrapolated_clock_mass_derivative']-
                                      by_level[64]['extrapolated_clock_mass_derivative']),
            tidal_ratio_grid_spread=abs(by_level[128]['tidal_ratio_coefficient_at_b_1e_4']-
                                       by_level[64]['tidal_ratio_coefficient_at_b_1e_4']))
        if summary[geometry]['derivative_grid_spread']>.02*max(
                1.,abs(by_level[128]['extrapolated_clock_mass_derivative'])):
            issues.append(geometry+': unresolved source derivative field resolution')
    delta={}
    for cells in (64,128):
        n,c=summary['nested']['levels'][cells],summary['cap']['levels'][cells]
        delta[str(cells)]=dict(clock_derivative=n['extrapolated_clock_mass_derivative']-c['extrapolated_clock_mass_derivative'],
            tidal_ratio_coefficient=n['tidal_ratio_coefficient_at_b_1e_4']-c['tidal_ratio_coefficient_at_b_1e_4'])
    delta['derivative_grid_spread']=abs(delta['128']['clock_derivative']-delta['64']['clock_derivative'])
    delta['tidal_ratio_grid_spread']=abs(delta['128']['tidal_ratio_coefficient']-delta['64']['tidal_ratio_coefficient'])
    return summary,delta,issues


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists() or not args.output.parent.is_dir():
        raise ValueError('use an unused compact output path with an existing parent')
    commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    hashes={}
    for name in SOURCES:
        raw=(ROOT/name).read_bytes()
        if subprocess.check_output(['git','show',f'{commit}:{name}'],cwd=ROOT)!=raw:
            raise ValueError('freeze the mass-response method before execution: '+name)
        hashes[name]=hashlib.sha256(raw).hexdigest()
    diagnosis_raw=(ROOT/DIAGNOSIS).read_bytes()
    if subprocess.check_output(['git','show',f'{commit}:{DIAGNOSIS}'],cwd=ROOT)!=diagnosis_raw:
        raise ValueError('freeze the prior precision diagnosis before execution')
    diagnosis_hash=hashlib.sha256(diagnosis_raw).hexdigest()
    started=time.perf_counter()
    workers=min(8,os.cpu_count() or 1)
    totals,metadata,rows={},[],[]
    def account(value):
        totals[value['pid']]=max(totals.get(value['pid'],0.),value['process_CPU_total'])
        if time.process_time()+sum(totals.values())+8>CPU_LIMIT:
            raise RuntimeError('Shared mass-response CPU allowance exhausted')
    pool=ProcessPoolExecutor(max_workers=workers,mp_context=multiprocessing.get_context('spawn'))
    try:
        preparation=[pool.submit(prepare,n,b) for n in (64,128) for b in MASSES]
        pending=[]
        for future in as_completed(preparation):
            value=future.result();account(value)
            data=value['metadata'];metadata.append(data)
            for geometry,model in value['models'].items():
                for order in (8,12):
                    pending.append(pool.submit(observe,model,data['cells'],data['central_mass'],
                        geometry,order,data['receiver_U'],data['receiver_V']))
            print(json.dumps(dict(prepared_cells=data['cells'],central_mass=data['central_mass'])),flush=True)
        for future in as_completed(pending):
            value=future.result();account(value)
            row=value['result'];rows.append(row)
            print(json.dumps(dict(cells=row['cells'],central_mass=row['central_mass'],geometry=row['geometry'],
                order=row['order'],clock=row['clock'],mass_gap=row['normalized_mass_path_gap'])),flush=True)
    except Exception as error:
        if hasattr(pool,'terminate_workers'):
            pool.terminate_workers()
        failure=dict(schema='nsc-atlas-dirac-mass-response-v2',status='execution_exit_recorded',
            source_commit=commit,source_sha256=hashes,error=str(error),results=rows,metadata=metadata,
            completed_worker_CPU_seconds=sum(totals.values()),new_P001_grids=0)
        args.output.write_text(json.dumps(failure,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
        raise SystemExit(str(error))
    finally:
        pool.shutdown(wait=True,cancel_futures=True)
    rows.sort(key=lambda r:(r['cells'],r['central_mass'],r['geometry'],r['order']))
    metadata.sort(key=lambda r:(r['cells'],r['central_mass']))
    summary,delta,issues=comparisons(rows,metadata)
    ell=10000.
    quantum=8*constants.pi*constants.hbar*constants.G/(constants.c**3*ell**2)
    physical={g:dict(tidal_ratio_state_difference=quantum*s['levels'][128]['tidal_ratio_coefficient_at_b_1e_4'],
        tidal_ratio_grid_spread=quantum*s['tidal_ratio_grid_spread'],
        clock_derivative_per_m=quantum*s['levels'][128]['extrapolated_clock_mass_derivative']/ell)
        for g,s in summary.items()}
    record=dict(schema='nsc-atlas-dirac-mass-response-v2',
        status='unresolved_prepared_mass_response' if issues else 'prepared_parametric_mass_response_evaluated',
        source_commit=commit,source_sha256=hashes,validation_issues=issues,
        model=dict(H_ell=1.,k_wall_ell=2.,initial_wall_radius=.45,
            source='central mass present before the window, initial proper wall speed fixed; parent mass follows Israel',
            central_mass_values=list(MASSES),angular_k=1,multiplet_degeneracy=2,
            preparation='same continuum unit-flux exterior mode and occupation difference as the mean-clock owner',
            order='leading 8*pi*(l_P/ell)^2 times occupation; parametric derivative in b/ell at zero',
            root_vacuum_selected=False,vacuum_angular_sum=False),
        observer=dict(U_over_duration=.8,V_over_duration=.7,phase='-V, emitted toward cap at the earlier wall time V',
            frame='fixed proper-wall radar position',
            cap_control='smooth SdS continuation with same classical and mean nested areal-worldtube motion; no membrane'),
        precision=dict(wall_rtol=2e-11,wall_atol=2e-13,ray_rtol=2e-8,ray_atol=2e-10,
            diagnosis=DIAGNOSIS,diagnosis_sha256=diagnosis_hash,
            reason='same-field wall mass/clock error floor; source, field grids, ray method and acceptance unchanged'),
        metadata=metadata,results=rows,comparison=summary,nesting_difference=delta,
        physical_scale=dict(ell_m=ell,central_mass_length_m=1.,
            scale_origin='supplied nominal scale, not calibrated to data',
            quantum_gravitational_factor=quantum,predictions=physical,
            constants_source='https://physics.nist.gov/cuu/Constants/Table/allascii.txt'),
        uncertainty='64/128 field, 8/12 quadrature and one-sided h/h2 derivative spreads are diagnostics, '
            'not proven continuum bounds or observational confidence intervals.',
        runtime=dict(wall_seconds=time.perf_counter()-started,
            observed_process_CPU_seconds=time.process_time()+sum(totals.values()),
            overhead_CPU_charge=8,charged_CPU_seconds=time.process_time()+sum(totals.values())+8,
            process_workers=workers,BLAS_threads_per_worker=1,Numba_threads_per_worker=1,
            new_P001_grids=0,new_P001_time_steps=0,baseline_vacuum_polarizations=0),
        limits=dict(charged_CPU_seconds=CPU_LIMIT,summary_bytes=32768),
        full_Psi_over_Phi=None,
        scope='Finite prepared-state parametric mean gravitational response to a physical mass family; '
              'local tidal ratio B/A is an operational curvature observable, not full vacuum Psi/Phi or retarded switching.')
    raw=json.dumps(record,indent=2,sort_keys=True,allow_nan=False)+'\n'
    if len(raw.encode())>32768:
        raise RuntimeError('Compact source-response summary exceeded 32 KiB')
    args.output.write_text(raw,encoding='utf-8',newline='\n')
    print(json.dumps(dict(status=record['status'],validation_issues=issues,comparison=summary,
        nesting_difference=delta,physical=physical,charged_CPU_seconds=record['runtime']['charged_CPU_seconds'])),flush=True)
    if issues:
        raise SystemExit(1)


if __name__=='__main__':
    main()
