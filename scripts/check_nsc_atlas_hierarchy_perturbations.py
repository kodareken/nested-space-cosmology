"""Finite perturbation response on the unchanged evolving source."""
import argparse
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict,replace
import hashlib
import json
from pathlib import Path
import time

import numpy as np

from recursive_horizons.nsc_atlas_hierarchy_perturbations import (
    PerturbationRun,run_perturbations,nonlinear_radial_check,multipole_potential)
from check_nsc_atlas_hierarchy_evolution import INPUTS as SOURCE_INPUTS

ROOT=Path(__file__).resolve().parents[1]
INPUTS=SOURCE_INPUTS+('src/recursive_horizons/nsc_atlas_hierarchy_perturbations.py',)


def worker(job):
    label,config=job
    key=hashlib.sha256(json.dumps(dict(config=asdict(config),sources={
        f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in INPUTS}),sort_keys=True).encode()).hexdigest()
    folder=ROOT/'.build/hierarchy-perturbations'/key;path=folder/'completed.json'
    if path.is_file():record=json.loads(path.read_text(encoding='utf-8'))
    else:
        folder.mkdir(parents=True,exist_ok=True);record,raw=run_perturbations(config,retain=True)
        fields=dict(radius=raw['radius'],initial_u=raw['initial_u'],final_u=raw['final_u'])
        for name,value in raw['initial_bases'].items():fields[name+'_seed']=value
        for name,value in raw['final_states'].items():fields[name+'_final']=value
        np.savez_compressed(folder/'fields.npz',**fields)
        (folder/'history.json').write_text(json.dumps(raw['histories'],allow_nan=False),encoding='utf-8')
        path.write_text(json.dumps(record,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    return label,record,str(folder)


def nonlinear_worker(job):
    config,folder,coefficients,epsilon=job;folder=Path(folder)
    key=hashlib.sha256(json.dumps(dict(config=asdict(config),coefficients=coefficients,epsilon=epsilon),sort_keys=True).encode()).hexdigest()
    path=folder/('nonlinear-'+key+'.json')
    if path.is_file():return json.loads(path.read_text(encoding='utf-8'))
    result=nonlinear_radial_check(config,folder/'fields.npz',coefficients,epsilon)
    path.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    return result


def calculation():
    wall=time.perf_counter();cpu=time.process_time();base=PerturbationRun()
    configs=dict(fine=replace(base,points=16383,step=.01),medium=base,
        coarse=replace(base,points=4095,step=.04),time_control=replace(base,step=.01),
        outer_control=replace(base,outer=1536.,points=12287))
    with ProcessPoolExecutor(max_workers=5) as pool:out=list(pool.map(worker,configs.items()))
    raw={name:r for name,r,_ in out};folders={name:f for name,_,f in out}
    previous=json.loads((ROOT/'results/summaries/nsc-atlas-hierarchy-evolution-v1.json').read_text(encoding='utf-8'))
    source_match={name:r['source_final_field_sha256']==previous['runs'][name]['final_field_sha256'] for name,r in raw.items()}
    # An outward-facing response: potential at twice the initial reservoir
    # radius, on a fixed physical sphere. It is an absolute quadrupolar
    # metric-potential coefficient, not an already integrated clock signal.
    prep=previous['preparation'];physical_cut=2*prep['reservoir_outer_radius']
    probe_scaled=physical_cut*1e6;probe={}
    for name,folder in folders.items():
        field=np.load(Path(folder)/'fields.npz');r=field['radius'];u=field['final_u']
        A=previous['runs'][name]['final_scale_factor'];b=previous['runs'][name]['background_envelope']
        dr=configs[name].outer/(configs[name].points+1);psi=b+u/r
        q=field['quadrupole_final'];v=multipole_potential(psi,q,r,dr,A,2)
        coefficient=np.array([np.interp(probe_scaled/A,r,v[:,j]) for j in range(4)])
        core_norm=raw[name]['responses']['quadrupole']['initial_core_field_norm']
        probe[name]=float(1e-6*core_norm*np.linalg.norm(coefficient))
    probe_error=2*max(abs(probe['fine']-probe['medium']),abs(probe['time_control']-probe['medium']),
        abs(probe['outer_control']-probe['medium']))
    centre=raw['fine']['responses'];ranges={}
    for mode in ('radial','quadrupole','late_quadrupole'):
        ranges[mode]={}
        for key in ('maximum_field_norm_gain','maximum_core_density_gain','maximum_quadrupole_moment_gain'):
            if centre[mode][key] is None:continue
            f=centre[mode][key];m=raw['medium']['responses'][mode][key]
            ranges[mode][key]=2*max(abs(f-m),abs(raw['time_control']['responses'][mode][key]-m),
                abs(raw['outer_control']['responses'][mode][key]-m))
    coefficients=raw['medium']['responses']['radial']['final']['worst_density_coefficients']
    jobs=[(base,folders['medium'],coefficients,e) for e in (1e-4,5e-5)]
    with ProcessPoolExecutor(max_workers=2) as pool:nonlinear=list(pool.map(nonlinear_worker,jobs))
    errors=[]
    if not all(source_match.values()):errors.append('spherical source replay changed')
    if max(x['maximum_symplectic_gap'] for x in centre.values())>1e-7:errors.append('symplectic transport unresolved')
    passive=raw['fine']['responses']['passive_quadrupole']
    if abs(passive['final']['field_norm_gain']-1)>1e-8:errors.append('passive angular propagation is not unitary')
    for key in ('relative_field_derivative_error','relative_core_density_derivative_error'):
        if nonlinear[1][key]>max(nonlinear[0][key]/3,2e-8):errors.append('nonlinear/tangent refinement failed: '+key)
        if nonlinear[1][key]>1e-3:errors.append('nonlinear derivative unresolved: '+key)
    total=time.process_time()-cpu+sum(r['cpu_seconds'] for r in raw.values())+sum(r['cpu_seconds'] for r in nonlinear)
    files=INPUTS+('scripts/check_nsc_atlas_hierarchy_perturbations.py','tests/test_nsc_atlas_hierarchy_perturbations.py',
        'results/summaries/nsc-atlas-hierarchy-evolution-v1.json')
    result=dict(schema='nsc-atlas-hierarchy-perturbations-v1',
        result='Finite real-linear radial and quadrupole response on the authenticated evolving envelope source.',
        source_replay_exact_matches=source_match,responses=centre,numerical_comparison_envelope=ranges,
        external_quadrupole_potential=dict(physical_probe_radius_atlas_units=physical_cut,
            maximum_Phi_over_c_squared_per_initial_fractional_field_norm=probe['fine'],
            numerical_comparison_envelope=probe_error,by_resolution=probe,
            meaning='Largest first-order angular-RMS quadrupolar Newtonian metric-potential coefficient at the final fixed physical probe sphere, per unit initial core-field-norm perturbation in the declared radial span. Multiply by epsilon for a chosen infinitesimal perturbation. This is not a redshift or lens image.'),
        runs=raw,nonlinear_radial_derivative_checks=nonlinear,
        scope='Two smooth compact radial profiles with widths2 and4, plus phase partners, real-orthonormalized. Radial seeds preserve total mass to first order. All five quadrupole orientations share the same radial propagator, so the four-column Gram bound extends to their20-dimensional real span. Late quadrupole seeds start at480; primary seeds start at0. Gains are maxima on saved unit-slow-time samples, not rigorous continuous-time enclosures or all-profile/asymptotic/nonlinear/full-GR stability.',
        field_versus_density='Field-norm amplification and density deformation are distinct. A large field norm can be mainly phase sensitivity; it is not alone evidence of breakup.',
        numerical_validation_errors=errors,process_cpu_seconds_including_workers=total,
        wall_seconds=time.perf_counter()-wall,compute_budget_seconds=4500,budget_exhausted=total>4500,
        timing_scope='CPU figures include underlying cached runs; wall time measures this invocation. Raw histories and tangent fields remain local.',
        source_sha256={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in files},
        reproduce='python scripts/run_python.py scripts/check_nsc_atlas_hierarchy_perturbations.py --output results/summaries/nsc-atlas-hierarchy-perturbations-v1.json')
    return result,folders


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path);args=p.parse_args()
    record,folders=calculation();encoded=json.dumps(record,indent=2,allow_nan=False)+'\n'
    if args.output:args.output.write_text(encoded,encoding='utf-8',newline='\n')
    else:print(encoded)
    print(json.dumps(dict(responses=record['responses'],errors=record['numerical_validation_errors'],cpu=record['process_cpu_seconds_including_workers'])))
    if record['numerical_validation_errors']:raise RuntimeError(record['numerical_validation_errors'])
    if record['budget_exhausted']:raise RuntimeError('4500 CPU-second budget exhausted; preserve completed cases and diagnose')
