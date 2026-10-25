"""Evolve the conserved core/reservoir and check the reduced equations."""
import argparse
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict,replace
import hashlib
import json
from pathlib import Path
import time

import numpy as np
import sympy as sp

from recursive_horizons.nsc_atlas_hierarchy_evolution import HierarchyRun,run_hierarchy,full_equation_comparison

ROOT=Path(__file__).resolve().parents[1]
INPUTS=('src/recursive_horizons/nsc_atlas_hierarchy_evolution.py',
    'src/recursive_horizons/nsc_atlas_scale_hierarchy.py','src/recursive_horizons/nsc_atlas_core_evolution.py',
    'src/recursive_horizons/nsc_atlas_core_reservoir.py','src/recursive_horizons/nsc_atlas_scalar_core.py',
    'src/recursive_horizons/nsc_atlas_shared_scalar_slice.py','src/recursive_horizons/nsc_atlas_flux_conversion.py',
    'src/recursive_horizons/nsc_atlas_flux_transition.py')


def symbolic_normalizations():
    A,w,z,k,d,C,L,E=sp.symbols('A w z k d C L E',positive=True)
    H2=L/3+k*E/(3*A**3)-C/A**2
    Hdot=C/A**2-k*E/(2*A**3)
    residual=-sp.Rational(3,2)*Hdot-sp.Rational(9,4)*H2
    amp=sp.sqrt(2/k)*z*z*A**(-sp.Rational(3,2))
    density=w*w*amp*amp/2
    mass_prefactor=sp.simplify(k*density/2*(A/(w*z))**3)
    background_amp=sp.sqrt(k/2)*d/(w*z*z)
    return dict(background_dust_cancellation=str(sp.simplify(residual-sp.Rational(3,4)*(C/A**2-L))),
        radial_geometric_mass_prefactor=str(sp.simplify(mass_prefactor-z/(2*w))),
        real_background_amplitude=str(sp.simplify(amp*background_amp-d/(w*A**sp.Rational(3,2)))))


def worker(job):
    name,kind,args=job
    key=hashlib.sha256(json.dumps(dict(kind=kind,args=args,
        sources={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in INPUTS}),sort_keys=True).encode()).hexdigest()
    folder=ROOT/'.build/hierarchy-evolution'/key;path=folder/'completed.json'
    if path.is_file():record=json.loads(path.read_text(encoding='utf-8'))
    else:
        folder.mkdir(parents=True,exist_ok=True)
        if kind=='envelope':
            record,raw=run_hierarchy(HierarchyRun(**args),retain=True)
            np.savez_compressed(folder/'fields.npz',radius=raw['radius'],initial_u=raw['initial_u'],final_u=raw['final_u'])
            (folder/'history.json').write_text(json.dumps(raw['history'],allow_nan=False),encoding='utf-8')
        else:record=full_equation_comparison(**args)
        path.write_text(json.dumps(record,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    return name,record,str(folder)


def calculation():
    wall=time.perf_counter();cpu=time.process_time();base=HierarchyRun();identities=symbolic_normalizations()
    configs=dict(fine=replace(base,points=16383,step=.01),medium=base,
        coarse=replace(base,points=4095,step=.04),time_control=replace(base,step=.01),
        outer_control=replace(base,outer=1536.,points=12287),no_gravity=replace(base,gravity=False),
        comparison_envelope_large=replace(base,frequency=2000.,weak_amplitude=.08,outer=160.,points=4095,step=.002,end=5.12),
        comparison_envelope_small=replace(base,frequency=8000.,weak_amplitude=.04,outer=320.,points=8191,step=.002,end=5.12))
    jobs=[(name,'envelope',asdict(config)) for name,config in configs.items()]
    jobs += [('comparison_GR_large','full',dict(zeta=.08,refinement=1)),
        ('comparison_GR_small','full',dict(zeta=.04,refinement=1)),
        ('comparison_GR_small_refined','full',dict(zeta=.04,refinement=2))]
    with ProcessPoolExecutor(max_workers=6) as pool:outputs=list(pool.map(worker,jobs))
    raw={name:r for name,r,_ in outputs};folders={name:f for name,_,f in outputs}
    fine=raw['fine'];medium=raw['medium'];envelopes={}
    history=json.loads((Path(folders['fine'])/'history.json').read_text(encoding='utf-8'))
    breathing={}
    for label,lo,hi in (('early',0.,120.),('late',480.,600.)):
        values=[r['physical_half_mass_radius'] for r in history if lo<=r['tau']<=hi]
        breathing[label]=dict(slow_time_interval=[lo,hi],minimum=min(values),maximum=max(values),
            peak_to_peak=max(values)-min(values))
    for key in ('late_mean_core_mass_ratio','late_minimum_core_mass_ratio','late_mean_physical_half_mass_radius','final_core_mass_ratio'):
        envelopes[key]=2*max(abs(fine[key]-medium[key]),
            abs(raw['time_control'][key]-medium[key]),abs(raw['outer_control'][key]-medium[key]))
    comparisons={}
    for label in ('large','small'):
        exact=raw['comparison_GR_'+label];approx=raw['comparison_envelope_'+label]
        comparisons[label]=dict(zeta=exact['zeta'],omega=exact['omega'],
            full_core_mass_ratio=exact['core_mass_ratio'],envelope_core_mass_ratio=approx['final_core_mass_ratio'],
            difference=exact['core_mass_ratio']-approx['final_core_mass_ratio'],
            scale_factor_difference=exact['final_scale_factor']-approx['final_scale_factor'],
            full_mass_flux_error=exact['maximum_mass_flux_gap_over_initial_core'])
    comparisons['small_full_refinement_gap']=raw['comparison_GR_small_refined']['core_mass_ratio']-raw['comparison_GR_small']['core_mass_ratio']
    errors=[]
    if any(v!='0' for v in identities.values()):errors.append('symbolic envelope normalization failed')
    if fine['maximum_energy_work_gap_over_initial_energy']>1e-4:errors.append('energy/expansion work unresolved')
    if max(fine['maximum_core_mass_flux_gaps_over_initial_core'])>1e-4:errors.append('fixed-physical-radius mass flux unresolved')
    if fine['maximum_total_mass_flux_gap_over_initial_core']>1e-6:errors.append('outer mass transport unresolved')
    if abs(comparisons['small']['difference'])>=abs(comparisons['large']['difference']):
        errors.append('short full-equation comparison does not improve in the weaker-field member')
    if abs(comparisons['small_full_refinement_gap'])>abs(comparisons['small']['difference'])/4:
        errors.append('full-equation comparison dominated by its numerical resolution')
    if min(raw[k]['finite_boundary_margin'] for k in ('comparison_GR_large','comparison_GR_small','comparison_GR_small_refined'))<=0:
        errors.append('full-equation comparison boundary margin unresolved')
    cpu_total=time.process_time()-cpu+sum(r['cpu_seconds'] for r in raw.values())
    short={name:{k:v for k,v in r.items() if k!='preparation'} for name,r in raw.items()}
    physical=json.loads((ROOT/'results/summaries/nsc-atlas-scale-hierarchy-v1.json').read_text(encoding='utf-8'))['illustrative_physical_conversion']
    ruler=physical['length_unit_metres'];c=physical['constants']['c_m_per_s'];pc=physical['constants']['parsec_metres']
    dimensional=dict(assumed_Lambda_per_m2=physical['assumed_Lambda_per_m2'],
        duration_gigayears=fine['end_atlas_proper_time']*ruler/c/(365.25*86400*1e9),
        fixed_core_mass_cut_kpc=20/(base.frequency*base.weak_amplitude)*ruler/(1000*pc),
        late_mean_half_mass_radius_kpc=fine['late_mean_physical_half_mass_radius']/(base.frequency*base.weak_amplitude)*ruler/(1000*pc),
        numerical_half_mass_radius_range_kpc=envelopes['late_mean_physical_half_mass_radius']/(base.frequency*base.weak_amplitude)*ruler/(1000*pc),
        interpretation='Same supplied dimensional reference as the initial hierarchy example; not an independent observed age or galaxy fit. The long run is an envelope calculation, with no rigorous full-GR model-error range.')
    source_files=INPUTS+('scripts/check_nsc_atlas_hierarchy_evolution.py','tests/test_nsc_atlas_hierarchy_evolution.py',
        'scripts/plot_nsc_atlas_hierarchy_evolution.py','results/summaries/nsc-atlas-scale-hierarchy-v1.json')
    record=dict(schema='nsc-atlas-hierarchy-evolution-v1',
        result='Finite classical-envelope localization of the retained high-frequency core and solved reservoir in its expanding inherited parent.',
        runs=short,numerical_comparison_envelope=envelopes,short_full_equation_comparison=comparisons,
        symbolic_zero_identities=identities,conditional_physical_interpretation=dimensional,
        resolved_radial_breathing=breathing,
        preparation=fine['preparation'],numerical_validation_errors=errors,
        process_cpu_seconds_including_workers=cpu_total,wall_seconds=time.perf_counter()-wall,
        longest_worker_original_wall_seconds=max(r['wall_seconds'] for r in raw.values()),
        compute_budget_seconds=1200,budget_exhausted=cpu_total>1200,
        timing_scope='Worker CPU includes original computations reused from authenticated local checkpoints. Wall time measures this invocation.',
        equations='i psi_tau=-Delta psi/(2 A^2)+V psi; Delta V=(abs(psi)^2-b^2)/(2 A); dE/dtau=-(A_prime/A)*(2T+W). The real scalar is reconstructed with its fast frequency and A^-3/2 amplitude.',
        limits='Leading weak-field, nonrelativistic and locally flat spatial envelope, with an averaged excited closed parent. This is not the full long-time Einstein-scalar evolution. Small expansion/curvature/weak-field parameters and a short full-equation comparison do not establish a rigorous model-error bound over the long run. Supplied initial core and reservoir; finite spherical localization is not formation from primordial data, asymptotic or nonspherical stability, a galaxy fit, or a dark-matter abundance.',
        source_sha256={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in source_files},
        reproduce='python scripts/run_python.py scripts/check_nsc_atlas_hierarchy_evolution.py --output results/summaries/nsc-atlas-hierarchy-evolution-v1.json')
    return record,folders


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path)
    p.add_argument('--figure',action='store_true');args=p.parse_args()
    record,folders=calculation();encoded=json.dumps(record,indent=2,allow_nan=False)+'\n'
    if args.output:args.output.write_text(encoded,encoding='utf-8',newline='\n')
    else:print(encoded)
    print(json.dumps(dict(result=record['runs']['fine'],errors=record['numerical_validation_errors'],cpu=record['process_cpu_seconds_including_workers'])))
    if record['numerical_validation_errors']:raise RuntimeError(record['numerical_validation_errors'])
    if record['budget_exhausted']:raise RuntimeError('1200 CPU-second budget exhausted; completed cases retained for diagnosis')
    if args.figure:
        from plot_nsc_atlas_hierarchy_evolution import render
        render(record,folders)
