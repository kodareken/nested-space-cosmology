"""Identical-covariance inherited-parent/LambdaCDM retrospective comparison.

No new field evolution, raw numerical history or observational holdout claim.
The pinned publisher input loader belongs to the original contrast producer.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
import time

import numpy as np
from scipy.stats import chi2
from threadpoolctl import threadpool_limits

from measure_nsc_atlas_distance_contrast import cached_file, FILES, RELEASE
from recursive_horizons.nsc_atlas_cosmological_distance import distance_contrast, project_distance_contrast
from recursive_horizons.nsc_atlas_inherited_matter import exact_initial_geometry, uniform_bound
from recursive_horizons.nsc_atlas_parent_supernova import (
    SupernovaShapeLikelihood, fit_shape, flrw_distance, magnitude_error_bound,
    parent_distance, parent_fractions, parent_ray_angle_upper,
)

ROOT=Path(__file__).resolve().parents[1]
SOURCES=(
    'src/recursive_horizons/nsc_atlas_parent_supernova.py',
    'src/recursive_horizons/nsc_atlas_inherited_matter.py',
    'src/recursive_horizons/nsc_atlas_cosmological_distance.py',
    'scripts/measure_nsc_atlas_distance_contrast.py',
    'scripts/measure_nsc_atlas_parent_supernova.py',
    'tests/test_nsc_atlas_parent_supernova.py',
)


def measure(cache, fetch=False, threads=4):
    started=time.process_time()
    bindings={name:cached_file(cache,name,fetch) for name in FILES}
    table=np.genfromtxt(cache/'distances.dat',names=True,dtype=None,encoding='utf-8')
    with (cache/'stat-plus-systematic.cov').open(encoding='ascii') as stream:
        dimension=int(stream.readline())
    C=np.loadtxt(cache/'stat-plus-systematic.cov',skiprows=1).reshape(dimension,dimension)
    if table.size!=1701 or dimension!=1701:
        raise ValueError('unexpected publisher dimensions')
    z,hel,mag,flags=(table[k] for k in ('zHD','zHEL','m_b_corr','IS_CALIBRATOR'))
    with threadpool_limits(limits=threads):
        primary,w,stat=distance_contrast(z,hel,mag,C,flags)
        extended=np.flatnonzero((flags==0)&(z>=.025))
        collections={}
        for name,indices in (('original_two_bands',primary),('all_noncalibrators_z_ge_0_025',extended)):
            zz,hh=z[indices],hel[indices]
            likelihood=SupernovaShapeLikelihood(mag[indices],C[np.ix_(indices,indices)])
            zmax=float(np.max(zz)); upper=1/(1+zmax)
            rows={}
            def put(label,distance,parameters,fitted=0,parent_s=None):
                predicted=distance(z[primary],hel[primary],64)
                contrast=project_distance_contrast(w,z[primary],hel[primary],predicted)
                dist=distance(zz,hh,64)
                numerical=float(np.max(abs(5*np.log10(distance(zz,hh,32)/dist))))
                row=dict(parameters=parameters, fitted_shape_parameters=fitted,
                    fitted_intercept_parameters=1, **likelihood.evaluate(dist),
                    original_contrast_prediction_magnitudes=contrast,
                    original_contrast_standardized_residual=(stat['contrast_magnitudes']-contrast)/stat['standard_uncertainty_magnitudes'],
                    quadrature_32_vs_64_max_magnitude_difference=numerical)
                dof=len(indices)-1-fitted
                row.update(degrees_of_freedom=dof,chi_squared_survival_probability=float(chi2.sf(row['chi_squared'],dof)))
                if parent_s is not None:
                    angle=parent_ray_angle_upper(zmax,parent_s)
                    if angle>=1:
                        raise ValueError('parent ray-angle comparison not established')
                    error=magnitude_error_bound()
                    row.update(parent_comparison_fractions=parent_fractions(parent_s),
                        minimum_emission_scale=(None if parent_s==0 else 1/(parent_s*(1+zmax))),
                        scalar_uniform_magnitude_error_up_to_common_offset=error,
                        original_contrast_scalar_model_error_bound=float(np.sum(abs(w)))*error,
                        ray_angle_upper=angle)
                rows[label]=row
            put('quiet_de_sitter',lambda z,h,o:flrw_distance(z,h,matter=0,order=o),{})
            for matter in (.3,.315,.334):
                put('flat_LCDM_fixed_'+str(matter),lambda z,h,o,m=matter:flrw_distance(z,h,matter=m,order=o),
                    dict(Omega_m=matter,Omega_k=0))
            lcdm=fit_shape(likelihood,lambda m:flrw_distance(zz,hh,matter=m),1.)
            m=lcdm['parameter']
            put('flat_LCDM_calibrated',lambda z,h,o:flrw_distance(z,h,matter=m,order=o),
                dict(Omega_m=m,Omega_k=0,at_boundary=lcdm['at_boundary']),1)
            rows['flat_LCDM_calibrated']['objective_evaluations']=lcdm['objective_evaluations']
            if .5<=upper:
                put('parent_fixed_Ao_2',lambda z,h,o:parent_distance(z,h,.5,order=o),dict(Ao=2.),parent_s=.5)
            parent=fit_shape(likelihood,lambda s:parent_distance(zz,hh,s),upper)
            s=parent['parameter']
            put('parent_calibrated_epoch',lambda z,h,o:parent_distance(z,h,s,order=o),
                dict(inverse_Ao=s,Ao=None if s==0 else 1/s,at_boundary=parent['at_boundary']),1,parent_s=s)
            rows['parent_calibrated_epoch']['objective_evaluations']=parent['objective_evaluations']
            d=exact_initial_geometry(); U=float(uniform_bound(10**9)['upper'])
            h2,matter,C0=map(float,(d['vacuum']/3,d['excitation']/3,d['curvature']))
            max_fraction=matter*U*upper**3/(h2+matter*U*upper**3-C0*upper**2)
            collections[name]=dict(light_curves=len(indices),distinct_CIDs=int(np.unique(table['CID'][indices]).size),
                selected_index_sha256=hashlib.sha256(indices.astype('<i8').tobytes()).hexdigest(),
                zHD_range=[float(np.min(zz)),zmax],parent_minimum_observer_scale=1+zmax,
                parent_scalar_fraction_future_upper_bound=max_fraction,models=rows)
            if time.process_time()-started>120:
                raise RuntimeError('comparison exceeded 120 process CPU seconds')
    return dict(schema='nsc-atlas-parent-supernova-v1',status='retrospective_common_likelihood_comparison',
        release_commit=RELEASE,sources=bindings,original_contrast=stat,comparisons=collections,
        physical_inputs=dict(omega=10**9,**{k:float(v) for k,v in exact_initial_geometry().items()}),
        scope=dict(data_usage='Retrospective calibration and model comparison, no withheld validation.',
            matter='Parent scalar excitation only; not a separately derived baryon/CDM mixture.',
            scale='Magnitude intercept removes common SI distance/Hubble scale; Lambda SI is not inferred.',
            history='Every parent emitter A>=1 in protected homogeneous future; no pre-slice extrapolation.',
            likelihood='Published standardized m_b_corr, joint STAT+SYS covariance, common intercept profiled.',
            errors='Quadrature difference is numerical indicator; scalar density/light-cone envelope is analytic. Neither covers missing matter/radiation/population physics.',
            sigma='Standardized residual of one fixed linear contrast; not a global-model sigma.'),
        computation=dict(process_CPU_seconds=time.process_time()-started,CPU_limit=120,BLAS_threads=threads,
                         parameter_grids=0,new_field_evolutions=0))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--fetch',action='store_true')
    parser.add_argument('--threads',type=int,default=4)
    args=parser.parse_args()
    out=(ROOT/args.output).resolve()
    if not out.is_relative_to(ROOT) or out.exists() or not out.parent.is_dir() or not 1<=args.threads<=8:
        raise ValueError('new repository output and 1..8 BLAS threads required')
    if subprocess.check_output(['git','status','--porcelain','--',*SOURCES],cwd=ROOT,text=True).strip():
        raise ValueError('commit method sources before comparing observations')
    cache=ROOT/'.build/data/pantheonplus-c447f0fea703'
    cache.mkdir(parents=True,exist_ok=True)
    result=measure(cache,args.fetch,args.threads)
    result['method_commit']=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    result['source_sha256']={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in SOURCES}
    result['reproduce']='python scripts/run_python.py scripts/measure_nsc_atlas_parent_supernova.py --fetch --output .build/parent-supernova-reproduced.json'
    with out.open('x',encoding='utf-8',newline='\n') as stream:
        stream.write(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n')
    print(json.dumps(dict(output=out.relative_to(ROOT).as_posix(),**result['computation'],comparisons={
        key:{model:dict(chi_squared=row['chi_squared'],contrast_sigma=row['original_contrast_standardized_residual'],
                        parameters=row['parameters']) for model,row in block['models'].items()}
        for key,block in result['comparisons'].items()})))


if __name__=='__main__':
    main()
