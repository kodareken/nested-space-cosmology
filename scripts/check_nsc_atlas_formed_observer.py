"""Fixed-clock observations of the retained forming concentration."""
import argparse
from concurrent.futures import ProcessPoolExecutor
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import time

import numpy as np

from recursive_horizons.nsc_atlas_formed_observer import SourceSlab,Observer,comparison,clock_age

ROOT=Path(__file__).resolve().parents[1]
INPUTS=('src/recursive_horizons/nsc_atlas_formed_observer.py','scripts/check_nsc_atlas_formed_observer.py',
        'scripts/plot_nsc_atlas_formed_observer.py','tests/test_nsc_atlas_formed_observer.py')


def retained_source(name):
    from check_nsc_atlas_seed_collapse import SIMULATION_INPUTS
    path=ROOT/'results/summaries/nsc-atlas-seed-collapse-v1.json';record=json.loads(path.read_text(encoding='utf-8'))
    row=record['runs'][name];source={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in SIMULATION_INPUTS}
    for f,h in source.items():
        if h!=record['source_sha256'][f]:raise RuntimeError('Source changed: '+f)
    key=hashlib.sha256(json.dumps(dict(args=row['config'],sources=source),sort_keys=True).encode()).hexdigest()
    raw=ROOT/'.build/seed-collapse'/key/'fields.npz'
    if not raw.is_file():raise FileNotFoundError('Regenerate the retained source with scripts/check_nsc_atlas_seed_collapse.py: '+str(raw))
    with np.load(raw) as fields:r=fields['radius'];u=fields['final_u']
    digest=hashlib.sha256(np.asarray(np.column_stack((r,u.real,u.imag)),dtype='<f8').tobytes()).hexdigest()
    if digest!=row['final_field_sha256']:raise RuntimeError('Retained final field mismatch: '+name)
    return r,u,row,digest


def worker(job):
    name,source_name,spacing,field_step,ray_tight=job;cpu=time.process_time();wall=time.perf_counter()
    r,u,row,digest=retained_source(source_name)
    slab=SourceSlab(r,u,length=row['config']['length'],spacing=spacing,field_step=field_step)
    obs=Observer(ray_rtol=2e-13,ray_atol=2e-15) if ray_tight else Observer()
    result=comparison(slab,observer=obs)
    result.update(source_name=source_name,source_final_field_sha256=digest,
                  slab_echo_relative_error=slab.echo_relative_error,slab_slices=len(slab.times),
                  field_step=field_step,slab_spacing=spacing,ray_tight=ray_tight)
    if name=='fine':
        age=result['signal']['source_proper_age'];derivatives={}
        derivative_obs=replace(obs,ray_rtol=2e-13,ray_atol=2e-15)
        for h in (1e-4,5e-5):
            plus=comparison(slab,age+h,derivative_obs);minus=comparison(slab,age-h,derivative_obs)
            value=(plus['proper_delay_excess']-minus['proper_delay_excess'])/(2*h)
            derivatives[str(h)]=dict(derivative=value,gap=value-result['endpoint_delay_derivative_prediction'],
                                     ray_rtol=derivative_obs.ray_rtol,ray_atol=derivative_obs.ray_atol)
        result['proper_age_derivative_checks']=derivatives
        angular={}
        for h in (.002,.001):
            plus=comparison(slab,age,replace(derivative_obs,separation=obs.separation+h))
            minus=comparison(slab,age,replace(derivative_obs,separation=obs.separation-h))
            value=(plus['proper_delay_excess']-minus['proper_delay_excess'])/(2*h)
            angular[str(h)]=dict(derivative=value,gap=value-result['endpoint_angular_derivative_prediction'])
        result['receiver_angle_derivative_checks']=angular
        result['frozen_comoving_potential']=comparison(slab,age,obs,frozen=True)
        fan=[]
        for offset in (.04,.08,.16):
            fan.append(dict(receiver_offset=offset,**comparison(slab,age,replace(obs,separation=np.pi-offset))))
        result['fixed_receiver_fan']=fan
    result['cpu_seconds']=time.process_time()-cpu;result['wall_seconds']=time.perf_counter()-wall
    return name,result


def calculation():
    start=time.process_time();wall=time.perf_counter()
    jobs=[('coarse','coarse',1e-5,2e-6,False),('medium','medium',1e-5,2e-6,False),
          ('fine','fine',1e-5,2e-6,False),('time_control','time_control',1e-5,2e-6,False),
          ('outer_control','outer_control',1e-5,2e-6,False),
          ('slab_control','fine',5e-6,1e-6,False),('ray_control','fine',1e-5,2e-6,True)]
    with ProcessPoolExecutor(max_workers=4) as pool:rows=dict(pool.map(worker,jobs))
    f=rows['fine'];m=rows['medium'];ranges={}
    for key in ('proper_delay_excess','image_direction_change','log_frequency_change','expansion_removed_log_frequency'):
        ranges[key]=2*max(abs(f[key]-m[key]),abs(rows['time_control'][key]-m[key]),
            abs(rows['outer_control'][key]-m[key]),abs(rows['slab_control'][key]-f[key]),abs(rows['ray_control'][key]-f[key]))
    units=json.loads((ROOT/'results/summaries/nsc-atlas-scale-hierarchy-v1.json').read_text(encoding='utf-8'))['illustrative_physical_conversion']
    lightunit=units['length_unit_metres']/units['constants']['c_m_per_s']
    mas=180/np.pi*3600*1000
    physical=dict(assumed_Lambda_per_m2=units['assumed_Lambda_per_m2'],
        image_shift_milliarcseconds=f['image_direction_change']*mas,
        image_shift_comparison_milliarcseconds=ranges['image_direction_change']*mas,
        proper_delay_days=f['proper_delay_excess']*lightunit/86400,
        proper_delay_comparison_days=ranges['proper_delay_excess']*lightunit/86400,
        equivalent_frequency_velocity_m_per_s=f['log_frequency_change']*units['constants']['c_m_per_s'],
        scope='Conditional finite-comoving-source/receiver clock experiment with the previous supplied Lambda ruler. Not an observed astronomical lens or a physical velocity measurement.')
    errors=[]
    for n,row in rows.items():
        if abs(row['signal']['fixed_receiver_residual'])>2e-10:errors.append(n+': receiver mismatch')
        if row['slab_echo_relative_error']>1e-8:errors.append(n+': source echo failed')
        if abs(row['signal']['endpoint_frequency_identity_gap']-row['homogeneous']['endpoint_frequency_identity_gap'])>5e-14:
            errors.append(n+': paired photon energy identity failed')
    for key in ranges:
        if ranges[key]>.1*abs(f[key]):errors.append(key+': signal unresolved')
    if min(abs(v['gap']) for v in f['proper_age_derivative_checks'].values())>5e-13:
        errors.append('fixed-clock arrival derivative disagrees with frequency')
    if min(abs(v['gap']) for v in f['receiver_angle_derivative_checks'].values())>2e-13:
        errors.append('angular arrival derivative disagrees with observed direction')
    total=time.process_time()-start+sum(r['cpu_seconds'] for r in rows.values())
    return dict(schema='nsc-atlas-formed-observer-v1',
        result='Fixed-observer direction, proper flight-time excess and frequency from the actual forming scalar source in a declared time-dependent weak metric.',
        observer=dict(comoving_radius_over_L=.25,source_azimuth=0.,receiver_azimuth=float(np.pi-.08),
            initial_clock_slice_N_offset=-.0002,emission_age='Homogeneous clock age from the initial slab slice to A=1.2; matched proper emitter age in every metric.'),
        runs=rows,numerical_comparison_ranges=ranges,conditional_physical_conversion=physical,
        dimensionless_vacuum_delay_relation=dict(Lambda_times_delay_squared=(9/32)*f['proper_delay_excess']**2,
            numerical_comparison=2*(9/32)*abs(f['proper_delay_excess'])*ranges['proper_delay_excess'],
            interpretation='Lambda_phys*(c*proper_delay)^2 for the fixed model and observer protocol. This relation does not select the absolute observed Lambda.'),
        numerical_validation_errors=errors,process_cpu_seconds_including_workers=total,wall_seconds=time.perf_counter()-wall,
        compute_budget_seconds=120,budget_exhausted=total>120,
        limits='Leading averaged scalar metric with Phi=Psi and retained parent spatial curvature. Null propagation is solved in that specified metric; omitted full-GR and carrier terms have no certified total error enclosure. Seed, ruler and observer locations are supplied; no observational fit or multiple-image census.',
        source_sha256={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in INPUTS},
        evidence_sha256={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in ('results/summaries/nsc-atlas-seed-collapse-v1.json','results/summaries/nsc-atlas-scale-hierarchy-v1.json')},
        reproduce='python scripts/run_python.py scripts/check_nsc_atlas_formed_observer.py --output results/summaries/nsc-atlas-formed-observer-v1.json --figure')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path)
    parser.add_argument('--pilot',action='store_true');parser.add_argument('--figure',action='store_true');args=parser.parse_args()
    if args.pilot:
        print(json.dumps(worker(('pilot','medium',1e-5,2e-6,False))[1],indent=2));raise SystemExit(0)
    result=calculation()
    if args.output:args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({k:result[k] for k in ('conditional_physical_conversion','numerical_comparison_ranges','numerical_validation_errors','process_cpu_seconds_including_workers')}))
    if result['numerical_validation_errors'] or result['budget_exhausted']:raise RuntimeError('Inspect retained observer result before further calculations')
    if args.figure:
        from plot_nsc_atlas_formed_observer import render
        render(result)
