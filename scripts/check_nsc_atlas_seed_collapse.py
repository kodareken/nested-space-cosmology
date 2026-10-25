"""Nonlinear compact-seed formation with source and conservation checks."""
import argparse
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict,replace
import hashlib
import json
from pathlib import Path
import time

import numpy as np

from recursive_horizons.nsc_atlas_seed_collapse import SeedRun,evolve_seed,initial_constraint
from recursive_horizons.nsc_atlas_formation_window import compact_packet

ROOT=Path(__file__).resolve().parents[1]
SIMULATION_INPUTS=('src/recursive_horizons/nsc_atlas_seed_collapse.py',
    'src/recursive_horizons/nsc_atlas_hierarchy_evolution.py','src/recursive_horizons/nsc_atlas_formation_window.py',
    'src/recursive_horizons/nsc_atlas_core_reservoir.py','src/recursive_horizons/nsc_atlas_scale_hierarchy.py',
    'src/recursive_horizons/nsc_atlas_shared_scalar_slice.py','src/recursive_horizons/nsc_atlas_flux_conversion.py',
    'src/recursive_horizons/nsc_atlas_flux_transition.py','src/recursive_horizons/nsc_atlas_scalar_core.py')
INPUTS=SIMULATION_INPUTS+('scripts/check_nsc_atlas_seed_collapse.py','scripts/plot_nsc_atlas_seed_collapse.py',
    'tests/test_nsc_atlas_seed_collapse.py')


def worker(job):
    name,args=job
    source={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in SIMULATION_INPUTS}
    key=hashlib.sha256(json.dumps(dict(args=args,sources=source),sort_keys=True).encode()).hexdigest()
    folder=ROOT/'.build/seed-collapse'/key;path=folder/'completed.json'
    if path.is_file():row=json.loads(path.read_text(encoding='utf-8'))
    else:
        folder.mkdir(parents=True,exist_ok=True)
        row,raw=evolve_seed(SeedRun(**args),retain=True)
        np.savez_compressed(folder/'fields.npz',radius=raw['radius'],initial_u=raw['initial_u'],final_u=raw['final_u'])
        (folder/'history.json').write_text(json.dumps(raw['history'],allow_nan=False),encoding='utf-8')
        path.write_text(json.dumps(row,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    return name,row,str(folder)


def calculation():
    start=time.process_time();wall=time.perf_counter();base=SeedRun(points=8191,log_step=.001,phase_cap=.04)
    configs=dict(coarse=replace(base,points=4095,log_step=.002,phase_cap=.08),medium=base,
        fine=replace(base,points=16383,log_step=.0005,phase_cap=.02),time_control=replace(base,log_step=.0005,phase_cap=.02),
        outer_control=replace(base,outer=12.,points=12287),no_gravity=replace(base,gravity=False),
        small_amplitude=replace(base,amplitude=.0005),half_small_amplitude=replace(base,amplitude=.00025),
        weak_coarse=replace(base,amplitude=.05,final_scale=16.,points=4095,log_step=.002,phase_cap=.08),
        weak_medium=replace(base,amplitude=.05,final_scale=16.),
        weak_fine=replace(base,amplitude=.05,final_scale=16.,points=16383,log_step=.0005,phase_cap=.02))
    with ProcessPoolExecutor(max_workers=6) as pool:outputs=list(pool.map(worker,[(n,asdict(c)) for n,c in configs.items()]))
    rows={n:r for n,r,_ in outputs};folders={n:Path(f) for n,_,f in outputs}
    from recursive_horizons.nsc_atlas_hierarchy_evolution import enclosed
    for name,row in rows.items():
        if row['final_central_density_contrast']<=0:
            row['final_physical_half_height_radius']=None
            row['scalar_mass_inside_half_height']=None
            continue
        with np.load(folders[name]/'fields.npz') as data:
            r=data['radius'];u=data['final_u'];dr=r[1]-r[0]
            digest=hashlib.sha256(np.asarray(np.column_stack((r,u.real,u.imag)),dtype='<f8').tobytes()).hexdigest()
            if digest!=row['final_field_sha256']:raise RuntimeError('retained field hash failed: '+name)
            cut=row['final_physical_half_height_radius']/(row['config']['final_scale']*row['config']['length'])
            row['scalar_mass_inside_half_height']=row['preparation']['mass_unit']*enclosed(r,r*r*abs(1+u/r)**2,dr,cut)
    fine=rows['fine'];medium=rows['medium'];errors=[]
    ranges={}
    for key in ('final_central_density_contrast','final_physical_mass_shell_radii','final_over_maximum_shell_radius','final_cut_excess_norms','final_physical_half_height_radius','scalar_mass_inside_half_height'):
        value=2*np.maximum.reduce([abs(np.array(fine[key])-np.array(medium[key])),
            abs(np.array(rows['time_control'][key])-np.array(medium[key])),
            abs(np.array(rows['outer_control'][key])-np.array(medium[key]))])
        ranges[key]=value.tolist()
    linear,raw=compact_packet(base.length,initial_scale=base.initial_scale,order=256,qmax=60.,rtol=2e-12)
    at_end=float(raw['weights']@raw['solution'].sol(1/base.final_scale)[:256])
    linear_checks={}
    for name in ('small_amplitude','half_small_amplitude'):
        row=rows[name];gain=row['final_central_density_contrast']/row['config']['amplitude']
        linear_checks[name]=dict(central_gain=gain,linear_gain=at_end,difference=gain-at_end)
    initial=fine['preparation'];tight=initial_constraint(base.amplitude,rtol=2e-13)[0]
    initial_control=dict(correction_difference=tight['correction']-initial['correction'],
        proper_mass_endpoint=tight['proper_distance_endpoint_norm'],
        corrected_mass_endpoint=tight['corrected_GR_endpoint_norm'])
    for name,row in rows.items():
        if row['maximum_total_flux_gap_over_initial_positive']>1e-5:errors.append(name+': total mass flux failed')
        if row['maximum_energy_work_gap_over_initial_energy']>.01:errors.append(name+': expansion work failed')
    if fine['maximum_mass_flux_gap_over_initial_positive']>.005:errors.append('fine local flux unresolved')
    if ranges['final_central_density_contrast']>.03*max(abs(fine['final_central_density_contrast']),1):errors.append('central density unresolved')
    if abs(linear_checks['half_small_amplitude']['difference'])>=abs(linear_checks['small_amplitude']['difference']):errors.append('linear limit does not improve')
    history=json.loads((folders['fine']/'history.json').read_text(encoding='utf-8'))
    # Keep only a small public plotting subset; complete fields/history regenerate locally.
    curves={}
    for name in ('fine','no_gravity','small_amplitude','weak_fine'):
        h=json.loads((folders[name]/'history.json').read_text(encoding='utf-8'))
        if name=='fine':h=[r for i,r in enumerate(h) if r['scale_factor']>=.3 or i%5==0]+[h[-1]]
        else:h=h[::max(1,len(h)//120)]+[h[-1]]
        curves[name]=dict(A=[r['scale_factor'] for r in h],contrast=[r['central_density_contrast'] for r in h],
            half_height_radius=[r['physical_half_height_radius'] if r['central_density_contrast']>0 else None for r in h],
            shell_radii=[r['physical_mass_shell_radii'] for r in h])
    units=json.loads((ROOT/'results/summaries/nsc-atlas-scale-hierarchy-v1.json').read_text(encoding='utf-8'))['illustrative_physical_conversion']
    Lunit=units['length_unit_metres'];pc=units['constants']['parsec_metres'];cc=units['constants']['c_m_per_s'];GM=units['constants']['nominal_solar_GM_m3_per_s2']
    physical=dict(assumed_Lambda_per_m2=units['assumed_Lambda_per_m2'],
        shell_label_initial_comoving_radius_over_L=.25,
        conserved_shell_mass_nominal_solar=cc*cc*Lunit*fine['conserved_shell_geometric_masses'][2]/GM,
        final_shell_radius_kpc=Lunit*fine['final_physical_mass_shell_radii'][2]/(1000*pc),
        maximum_shell_radius_kpc=Lunit*fine['maximum_physical_mass_shell_radii'][2]/(1000*pc),
        final_half_height_radius_kpc=Lunit*fine['final_physical_half_height_radius']/(1000*pc),
        final_shell_radius_comparison_kpc=Lunit*ranges['final_physical_mass_shell_radii'][2]/(1000*pc),
        final_half_height_radius_comparison_kpc=Lunit*ranges['final_physical_half_height_radius']/(1000*pc),
        scalar_mass_inside_half_height_nominal_solar=cc*cc*Lunit*fine['scalar_mass_inside_half_height']/GM,
        scalar_mass_inside_half_height_comparison_nominal_solar=cc*cc*Lunit*ranges['scalar_mass_inside_half_height']/GM,
        scope='Conditional conversion using the existing supplied Lambda ruler. The selected shell mass includes ambient material; no galaxy mass or size was fitted.')
    total=time.process_time()-start+sum(r['cpu_seconds'] for r in rows.values())
    return dict(schema='nsc-atlas-seed-collapse-v1',
        result='Finite nonlinear envelope evolution of an initially weak compact compensated density seed on the declared earlier parent.',
        runs=rows,numerical_comparison_ranges=ranges,initial_constraint_control=initial_control,
        conditional_physical_conversion=physical,
        finite_formation_assessment=dict(
            tested_strong_amplitude=.1,tested_weak_amplitude=.05,
            strong_finite_contraction=bool(fine['final_over_maximum_shell_radius'][2]+ranges['final_over_maximum_shell_radius'][2]<.2),
            weak_sampled_shells_still_expanding=bool(min(rows['weak_fine']['final_over_maximum_shell_radius'])==1.),
            no_gravity_sampled_shells_still_expanding=bool(min(rows['no_gravity']['final_over_maximum_shell_radius'])==1.),
            mass_definition='Fixed integral of the oscillatory scalar density in the envelope, in geometric units; excludes vacuum volume mass and higher-order GR gradient terms.',
            limits='Finite regular concentration and rebound for a separately supplied 10% seed; not an amplitude threshold, asymptotic bound state or exact full-GR evolution.'),
        linear_limit_checks=linear_checks,curves=curves,
        numerical_validation_errors=errors,process_cpu_seconds_including_workers=total,
        wall_seconds=time.perf_counter()-wall,compute_budget_seconds=1200,budget_exhausted=total>1200,
        limits='Earlier homogeneous parent and compact seed remain supplied. Exact local initial turning-point constraints and approximate nonlinear envelope evolution are distinguished. No shared-wall past, nonspherical or asymptotic stability, abundance or observational confirmation.',
        source_sha256={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in INPUTS},
        evidence_sha256={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in ('results/summaries/nsc-atlas-formation-window-v1.json',)},
        reproduce='python scripts/run_python.py scripts/check_nsc_atlas_seed_collapse.py --output results/summaries/nsc-atlas-seed-collapse-v1.json --figure')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path)
    p.add_argument('--pilot',action='store_true');p.add_argument('--figure',action='store_true');args=p.parse_args()
    if args.pilot:
        name,row,folder=worker(('coarse',asdict(SeedRun())))
        print(json.dumps(dict(result=row,folder=folder),indent=2));raise SystemExit(0)
    r=calculation()
    if args.output:args.output.write_text(json.dumps(r,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(errors=r['numerical_validation_errors'],cpu=r['process_cpu_seconds_including_workers'],
        central=r['runs']['fine']['final_central_density_contrast'],ranges=r['numerical_comparison_ranges'])))
    if r['numerical_validation_errors'] or r['budget_exhausted']:raise RuntimeError('Inspect retained diagnosis before more calculations')
    if args.figure:
        from plot_nsc_atlas_seed_collapse import render
        render(r)
