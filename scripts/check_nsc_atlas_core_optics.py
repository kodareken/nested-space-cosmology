"""Lensing and curvature focusing from the frozen scalar-parent evolution."""
import argparse
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict,replace
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np
import sympy as sp

from recursive_horizons.nsc_atlas_core_optics import OpticalRun,run_optics,homogeneous_screen

ROOT=Path(__file__).resolve().parents[1]


def identities():
    R,m,k,U,P,F,b,E,q,L=sp.symbols('R m k U P F b E q L',real=True)
    rho=(P*P+b*b*F*F)/2+U;pr=rho-2*U
    hnn=-m/R**2-k*R*pr/2;hss=m/R**2-k*R*rho/2;hns=-k*R*P*b*F/2
    warped=-(E*E*hnn+2*E*q*hns+q*q*hss)/R+2*m*L*L/R**5
    W=m/R**3-k*(U+(P*P-b*b*F*F)/2)/6
    tidal=k*(E*P+q*b*F)**2/2+3*W*L*L/R**2
    # Terms proportional to E^2-q^2-L^2/R^2 vanish on the null mass shell.
    difference=sp.expand(warped-tidal).subs(E**2,q*q+L*L/R**2)
    v,gam,mu=sp.symbols('v gam mu',real=True)
    w=-v/b;gamma2=1/(1-w*w)
    En=gam*(1+w*mu);qn=gam*(mu+w)
    raw=dict(warped_product_optical_eigenvalue=difference,
        vacuum_null_contraction=U*(-E*E+q*q+L*L/R**2).subs(E**2,q*q+L*L/R**2),
        source_boost_mass_shell=sp.expand(En*En-qn*qn-(1-mu*mu)).subs(gam**2,gamma2),
        source_energy_normalization=sp.expand(gam*(En-w*qn)-1).subs(gam**2,gamma2))
    dU=sp.symbols('dU')
    raw['uniform_vacuum_Weyl_cancellation']=W.subs({U:U+dU,m:m+k*dU*R**3/6},simultaneous=True)-W
    answer={name:str(sp.simplify(value)) for name,value in raw.items()}
    if any(value!='0' for value in answer.values()):raise RuntimeError(answer)
    return answer


def index(record):return {(row['emission_proper_age'],round(row['sky_angle'],12)):row for row in record}


def map_checks(record):
    lookup=index(record['rays']);rows=[];delta=.0002
    for age in (2.8,3.4):
        for angle in (.02,.04,.08):
            r=dict(lookup[(age,angle)])
            phi=lambda x:lookup[(age,round(x,12))]['landing_angle']
            d1=(phi(angle+delta)-phi(angle-delta))/(2*delta)
            d2=(phi(angle+delta/2)-phi(angle-delta/2))/delta
            derivative=(4*d2-d1)/3
            angular=-.6*r['received_radial_cosine']*derivative
            r.update(angular_jacobi_parallel=angular,
                angular_jacobi_gap=r['jacobi_parallel']-angular,
                angular_difference_estimate=abs(.6*r['received_radial_cosine']*(d2-d1)/3))
            rows.append(r)
    return rows


def worker(item):
    label,config=item
    names=('src/recursive_horizons/nsc_atlas_core_optics.py','src/recursive_horizons/nsc_atlas_core_evolution.py',
           'src/recursive_horizons/nsc_atlas_core_reservoir.py','src/recursive_horizons/nsc_atlas_scalar_core.py',
           'src/recursive_horizons/nsc_atlas_shared_scalar_slice.py')
    key=hashlib.sha256(json.dumps(dict(config=asdict(config),
        sources={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in names}),sort_keys=True).encode()).hexdigest()
    folder=ROOT/'.build/core-optics'/key
    path=folder/'completed.json'
    if path.is_file():result=json.loads(path.read_text(encoding='utf-8'))
    else:result=run_optics(config,cache_dir=folder)
    return label,result


def calculation():
    cpu0=time.process_time();wall=time.perf_counter();symbolic=identities()
    base=OpticalRun()
    configs={
        'fine':replace(base,spacing=.00025,step=.000025),
        'medium':base,
        'coarse':replace(base,spacing=.001,step=.0001),
        'optical_time_control':replace(base,block_steps=2),
        'outer_control':replace(base,outer=4.5),
        'homogeneous':replace(base,kind='homogeneous',spacing=.002,step=.0001),
    }
    with ProcessPoolExecutor(max_workers=6) as pool:outputs=list(pool.map(worker,configs.items()))
    raw=dict(outputs);full_hom=homogeneous_screen();hom=index(full_hom)
    previous=json.loads((ROOT/'results/summaries/nsc-atlas-core-evolution-v1.json').read_text(encoding='utf-8'))
    expected={'fine':'diffuse_fine','medium':'diffuse_medium','coarse':'diffuse_coarse',
              'optical_time_control':'diffuse_medium','outer_control':'diffuse_outer_control','homogeneous':'homogeneous_late'}
    exact={name:raw[name]['final_field_sha256']==previous['runs'][old]['final_state_sha256'] for name,old in expected.items()}
    probes={name:map_checks(r) for name,r in raw.items()}
    predictions=[]
    for j,p in enumerate(probes['fine']):
        ref=hom[(p['emission_proper_age'],p['sky_angle'])]
        def observables(row):
            area=row['signed_beam_area']/ref['signed_beam_area']
            frequency=row['received_frequency']/ref['received_frequency']
            return dict(deflection_radians=row['landing_angle']-ref['landing_angle'],
                relative_signed_area=area,relative_beam_flux=frequency**2/abs(area),
                frequency_departure=frequency-1,arrival_proper_delay=row['arrival_proper_age']-ref['arrival_proper_age'],
                axis_ratio=abs(row['jacobi_parallel']/row['jacobi_perpendicular']))
        values=observables(p);medium=observables(probes['medium'][j]);control=[observables(probes[name][j]) for name in ('optical_time_control','outer_control')]
        errors={key:2*max(abs(values[key]-medium[key]),*[abs(v[key]-medium[key]) for v in control]) for key in values}
        # Optical identities provide a second, independently computed area.
        independent_area=p['rotation_jacobi']*p['angular_jacobi_parallel']/ref['signed_beam_area']
        independent_flux=(p['received_frequency']/ref['received_frequency'])**2/abs(independent_area)
        errors['relative_signed_area']=max(errors['relative_signed_area'],2*abs(independent_area-values['relative_signed_area']))
        errors['relative_beam_flux']=max(errors['relative_beam_flux'],2*abs(independent_flux-values['relative_beam_flux']))
        predictions.append(dict(emission_proper_age=p['emission_proper_age'],sky_angle=p['sky_angle'],
            **values,numerical_comparison_envelope=errors,
            parity_relative_to_homogeneous=1 if values['relative_signed_area']>0 else -1))
    homogeneous_checks={key:max(abs(row[key]-hom[(row['emission_proper_age'],round(row['sky_angle'],12))][key]) for row in raw['homogeneous']['rays'])
        for key in ('arrival_proper_age','landing_angle','received_frequency','jacobi_perpendicular')}
    fan={}
    for age in (2.8,3.4):
        rows=[r for r in raw['fine']['rays'] if r['emission_proper_age']==age]
        brackets=[]
        for a,b in zip(rows[:-1],rows[1:]):
            if a['jacobi_perpendicular']*b['jacobi_perpendicular']<0:brackets.append([a['sky_angle'],b['sky_angle']])
        fan[str(age)]=dict(tangential_caustic_angle_brackets=brackets,
            minimum_signed_beam_area=min(r['signed_beam_area'] for r in rows),
            angle_range=[rows[0]['sky_angle'],rows[-1]['sky_angle']])
    compact={name:dict(config=r['config'],final_field_sha256=r['final_field_sha256'],cpu_seconds=r['cpu_seconds'],
        wall_seconds=r['wall_seconds'],probes=probes[name],
        maximum_rotation_jacobi_gap=max(abs(p['rotation_jacobi_gap']) for p in r['rays']),
        maximum_frequency_transport_gap=max(abs(p['frequency_transport_gap']) for p in r['rays'])) for name,r in raw.items()}
    cpu=time.process_time()-cpu0+sum(r['cpu_seconds'] for r in raw.values())
    files=('src/recursive_horizons/nsc_atlas_core_optics.py','src/recursive_horizons/nsc_atlas_core_evolution.py',
        'src/recursive_horizons/nsc_atlas_core_reservoir.py','src/recursive_horizons/nsc_atlas_shared_scalar_slice.py',
        'scripts/check_nsc_atlas_core_optics.py','scripts/plot_nsc_atlas_core_optics.py','tests/test_nsc_atlas_core_optics.py',
        'results/summaries/nsc-atlas-core-evolution-v1.json')
    errors=[]
    if not all(exact.values()):errors.append('streamed source does not reproduce the frozen field bit for bit')
    convergence=[max(abs(p['rotation_jacobi_gap']) for p in probes[name]) for name in ('coarse','medium','fine')]
    if not convergence[0]>convergence[1]>convergence[2]:errors.append('curvature/ray-map convergence failed')
    if convergence[2]>3e-6:errors.append('rotation Jacobi discrepancy unresolved')
    if max(abs(p['angular_jacobi_gap']) for p in probes['fine'])>5e-6:errors.append('parallel beam scale unresolved')
    if max(abs(p['frequency_transport_gap']) for p in probes['fine'])>1e-7:errors.append('photon frequency mismatch')
    if homogeneous_checks['jacobi_perpendicular']>2e-6:errors.append('homogeneous analytic beam check failed')
    record=dict(schema='nsc-atlas-core-optics-v1',result='Nonradial light bending, anisotropic focusing and received beam flux from the same frozen evolved scalar geometry.',
        symbolic_zero_identities=symbolic,source_replay_exact_matches=exact,runs=compact,
        predictions=predictions,angular_fan=fan,homogeneous_analytic_checks=homogeneous_checks,
        protocol='Held areal-radius .6 source and receiving sphere, clocks zero on initial slice. Emit unit rest-frame energy at proper ages2.8 and3.4, inward sky angle psi. Follow each ray to its first outward sphere recrossing. Comparisons identify rays by the same source event protocol and launch angle, not by an assumed common receiving worldline.',
        flux_definition='Received bolometric beam flux for fixed emitted luminosity is proportional to nu_received^2/abs(D_perpendicular*D_parallel). Signed area records parity separately. This finite-screen beam comparison is not a fit to an astronomical lens catalogue.',
        direct_focusing='R_ab k^a k^b=kappa(k phi)^2; potential cancels from this contraction. It still affects field evolution and geometry. Weyl focusing supplies the opposite-sign optical shear eigenvalues.',
        limits='Same supplied nonlinear cloud, action and parent; no new physical scale or abundance selection. Geometric optics and finite spherical evolution. Caustic brackets do not predict infinite physical brightness or resolve finite-source/wavelength smoothing. Error envelopes are finite numerical comparisons, not rigorous bounds or preparation uncertainty.',
        numerical_validation_errors=errors,worker_count=6,process_cpu_seconds_including_workers=cpu,
        wall_seconds=time.perf_counter()-wall,compute_budget_seconds=4500,budget_exhausted=cpu>4500,
        source_sha256={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in files},
        reproduce='python scripts/run_python.py scripts/check_nsc_atlas_core_optics.py --figure',
        next='Use the same ray map for fixed source/receiver imaging and compare its frequency/focusing relation against alternative admissible preparations; population selection and nonspherical longevity remain separate.')
    return record,raw,full_hom


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path);p.add_argument('--figure',action='store_true');args=p.parse_args()
    result,raw,hom=calculation();encoded=json.dumps(result,indent=2,allow_nan=False)+'\n'
    if len(encoded.encode())>100000:raise RuntimeError('100 KB compact optical record limit exceeded')
    if args.output:args.output.write_text(encoded,encoding='utf-8',newline='\n')
    else:print(encoded,end='')
    # Completed checkpoints and the diagnosis survive any failed check or budget.
    if result['numerical_validation_errors']:raise RuntimeError('; '.join(result['numerical_validation_errors']))
    if result['budget_exhausted']:raise RuntimeError('optical campaign exceeded4500 CPU seconds; completed checkpoints retained for diagnosis')
    if args.figure:
        from plot_nsc_atlas_core_optics import render
        render(result,raw,hom)
