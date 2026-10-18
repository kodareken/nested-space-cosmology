"""Bounded full-radial, finite-split vacuum comparison on the actual atlas.

The state is the declared correlated cylinder preparation. The first signed
angular sector is complete radially but is not the angular sum, a coincident
renormalized stress, or the source-induced gravitational potential ratio.
"""
from __future__ import annotations

import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'):
    os.environ[name]='1'

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import time

import numpy as np
from scipy import constants
from scipy.integrate import quad_vec
from threadpoolctl import threadpool_limits

from recursive_horizons.nsc_atlas_dirac_radar import WallRadarChart
from recursive_horizons.nsc_atlas_dirac_preparation import CylinderPreparation
from recursive_horizons.nsc_atlas_dirac_retarded import retarded_null_responses,retarded_response
from recursive_horizons.nsc_atlas_dirac_root_pairing import NullResponse,pair_cylinder_responses

ROOT=Path(__file__).resolve().parents[1]
CPU_LIMIT=240.
SOURCES=tuple('src/recursive_horizons/'+name+'.py' for name in (
    'nsc_atlas_dirac_retarded','nsc_atlas_dirac_root_pairing','nsc_atlas_dirac_preparation',
    'nsc_atlas_dirac_mean','nsc_atlas_dirac_scattering','nsc_atlas_dirac_radar',
    'nsc_atlas_cap_wall_response','nsc_atlas_causal_reference','nsc_atlas_dirac_evolution'))+(
    'tests/test_nsc_atlas_dirac_retarded.py','tests/test_nsc_atlas_dirac_root_pairing.py',
    'scripts/measure_nsc_atlas_dirac_root_response.py')


def difference(a,b):
    if a.lengths!=b.lengths:
        raise ValueError('relative root pairing requires identical null endpoints')
    return NullResponse(a.lengths,tuple(x-y for x,y in zip(a.green,b.green)),
        a.direct-b.direct,a.direct_derivative-b.direct_derivative)


def relative_pair(a,b,c,d,mu):
    # A C_root B^dagger-S C_root T^dagger
    #       = (A-S) C_root B^dagger + S C_root (B-T)^dagger.
    # Cancel common direct distributions before doing numerical integration.
    first=pair_cylinder_responses(difference(a,c),b,mu,epsabs=2e-12)
    second=pair_cylinder_responses(c,difference(b,d),mu,epsabs=2e-12)
    return (first['kernel']+second['kernel'],
        first['pv_quadrature_estimate']+second['pv_quadrature_estimate'],
        first['bessel_series_tail_bound']+second['bessel_series_tail_bound'])


def endpoint(preparation,U,V,n,integral_cache):
    du,dv=U/(n-.5),V/(n-.5)
    u=(np.arange(n)+.5)*du;v=(np.arange(n)+.5)*dv
    fields=preparation.fields(*np.meshgrid(u,v,indexing='ij'))
    if float(preparation.fields(U,V)['preparation_weight'])!=1:
        raise ValueError('physical endpoints must lie in the matched operator region')
    cache_key=(id(preparation),U,V)
    integrals=[];errors=[]
    for leg,length in enumerate((V,U)) if cache_key not in integral_cache else ():
        def integrand(s):
            f=preparation.fields(U,s) if leg==0 else preparation.fields(s,V)
            W=float(f['potential'])
            derivative=float((f['potential_t']-f['potential_x'])/2 if leg==0 else
                             (f['potential_t']+f['potential_x'])/2)
            return np.array([W*W,W*derivative])
        points=[t for t in (.1*preparation.duration,.2*preparation.duration,
            .9*preparation.duration,U if leg==0 else V) if 0<t<length]
        value,error=quad_vec(integrand,0,length,points=sorted(set(points)),
            epsabs=2e-12,epsrel=2e-11,limit=160)
        integrals.append(value);errors.append(float(error))
    if cache_key not in integral_cache:
        integral_cache[cache_key]=(integrals,errors)
    integrals,errors=integral_cache[cache_key]
    result=retarded_null_responses(fields['potential'],du,dv,
        potential_u=(fields['potential_t']-fields['potential_x'])/2,
        potential_v=(fields['potential_t']+fields['potential_x'])/2,
        boundary_square_integrals=[v[0] for v in integrals],
        boundary_mixed_integrals=[v[1] for v in integrals],order=24,workers=8)
    result['endpoint_quadrature_estimate']=max(errors)
    return result


def packed(matrix):
    return dict(real=matrix.real.tolist(),imag=matrix.imag.tolist())


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists() or not args.output.parent.is_dir():
        raise ValueError('unused output path and existing parent required')
    commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    hashes={}
    for name in SOURCES:
        raw=(ROOT/name).read_bytes()
        if raw!=subprocess.check_output(['git','show',f'{commit}:{name}'],cwd=ROOT):
            raise ValueError('freeze the response method before execution: '+name)
        hashes[name]=hashlib.sha256(raw).hexdigest()
    cpu_start=time.process_time();wall_start=time.perf_counter()
    def budget():
        if time.process_time()-cpu_start>CPU_LIMIT:
            raise RuntimeError('Full-root response exhausted its 240 process-CPU seconds; stop before another calculation')
    chart=WallRadarChart();S=chart.duration;U,V=.8*S,.7*S
    T,x=(U+V)/2,(V-U)/2;ell=10000.;mu=1/.45
    metric=chart.fields(U,V,region='cap');a=float(metric['lapse']);r=float(metric['radius'])
    preparations={g:CylinderPreparation(chart,duration=S,root_radius=.45,region=g,angular_k=1)
                  for g in ('nested','cap')}
    rows=[];endpoint_evidence=[];stopped=None;integral_cache={}
    try:
        with threadpool_limits(limits=1):
            for n in (64,128):
                for fraction in (.02,.01):
                    eta=fraction*S
                    points=((U-eta/2,V+eta/2),(U+eta/2,V-eta/2))
                    maps={}
                    for g,preparation in preparations.items():
                        for side,point in enumerate(points):
                            budget();sample=endpoint(preparation,*point,n,integral_cache);budget()
                            maps[g,side]=sample
                            endpoint_evidence.append(dict(cells=n,separation_over_S=fraction,
                                geometry=g,side=side,quadrature_estimate=sample['endpoint_quadrature_estimate'],
                                reconstruction_projection_residual=sample['cell_projection_residual'].tolist(),
                                response_sha256=hashlib.sha256(b''.join(
                                    c.tobytes() for key in ('cell_field','cell_spatial')
                                    for c in sample[key].green)).hexdigest()))
                            print(json.dumps(dict(cells=n,separation_over_S=fraction,geometry=g,side=side)),flush=True)
                    separation,separation_error=quad_vec(lambda y:float(chart.fields(T-y,T+y,region='cap')['lapse']),
                        x-eta/2,x+eta/2,epsabs=1e-13,epsrel=1e-12)
                    for order in (12,24):
                        budget()
                        lf,rf,ls,rs=[retarded_response(maps[g,side],'field',order)
                            for g,side in (('nested',0),('nested',1),('cap',0),('cap',1))]
                        ld,rd,sd,td=[retarded_response(maps[g,side],'spatial',order)
                            for g,side in (('nested',0),('nested',1),('cap',0),('cap',1))]
                        kernel,ke,kb=relative_pair(lf,rf,ls,rs,mu)
                        dl,de,db=relative_pair(ld,rf,sd,rs,mu)
                        dr,re,rb=relative_pair(lf,rd,ls,td,mu)
                        derivative=(dl-dr)/2
                        q=float((-1j*np.trace(np.diag([1.,-1.])@derivative)).real)
                        s=float(np.trace(np.array([[0.,1.],[1.,0.]])@kernel).real)
                        radial=2/(4*np.pi*a*a*r*r);angular=2/(8*np.pi*a*r**3)
                        pressure=radial*q-angular*s
                        rows.append(dict(cells=n,polynomial_order=order,separation_over_S=fraction,
                            proper_separation_m=float(separation)*ell,
                            proper_separation_quadrature_estimate_m=float(separation_error)*ell,
                            nested_minus_cap_kernel=packed(kernel),
                            nested_minus_cap_split_derivative=packed(derivative),
                            finite_split_anisotropy_difference=pressure,
                            nominal_finite_split_pressure_Pa=pressure*constants.hbar*constants.c/ell**4,
                            pv_pressure_quadrature_estimate=2*radial*(de+re)/2+2*angular*ke,
                            bessel_pressure_tail_bound=2*radial*(db+rb)/2+2*angular*kb))
                        budget()
    except Exception as error:
        stopped=type(error).__name__+': '+str(error)
    comparisons={}
    if stopped is None:
        lookup={(v['cells'],v['polynomial_order'],v['separation_over_S']):v['finite_split_anisotropy_difference'] for v in rows}
        comparisons=dict(grid_spread={f'{p}:{f}':abs(lookup[128,p,f]-lookup[64,p,f])
            for p in (12,24) for f in (.02,.01)},
            polynomial_spread={f'{n}:{f}':abs(lookup[n,24,f]-lookup[n,12,f]) for n in (64,128) for f in (.02,.01)},
            separation_spread={f'{n}:{p}':abs(lookup[n,p,.01]-lookup[n,p,.02]) for n in (64,128) for p in (12,24)})
    result=dict(schema='nsc-atlas-dirac-full-root-finite-split-v1',method_commit=commit,
        sources_sha256=hashes,status='finite_split_reference_only' if stopped is None else 'stopped_before_complete_comparison',
        stopped_reason=stopped,domain=dict(wall_clock_duration=S,root_radius=.45,
            detector_U_over_S=.8,detector_V_over_S=.7,physical_reference_length_m=ell,
            detector_radius_over_ell=r,physical_radar_lapse=a,central_mass=0.,angular_sectors=[1],
            magnetic_degeneracy=2,root_radial_spectrum='complete analytic covariance, no momentum cutoff',
            state='same correlated product-cylinder parent and fixed auxiliary switch in both geometries'),
        cases=rows,endpoint_evidence=endpoint_evidence,numerical_spreads=comparisons,
        limits=['One signed angular sector only; the angular tail is not computed.',
            'Finite nonzero split in the declared canonical radar frame; no coincident renormalized stress.',
            'The pressure-shaped bilocal statistic is not full Psi/Phi or its compact-source derivative.',
            'Grid, polynomial and separation spreads are diagnostics, not certified continuum error bounds.',
            'Bessel tail and PV estimates exclude field representation and floating arithmetic.',
            'The physical formation of the supplied preparation remains open; no P001 grid ran.'],
        compute=dict(mode_workers=8,BLAS_threads=1,process_CPU_seconds=time.process_time()-cpu_start,
            CPU_limit_seconds=CPU_LIMIT,wall_seconds=time.perf_counter()-wall_start,
            live_field_array_ceiling_bytes=256*1024**2,retained_full_arrays=False))
    raw=json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n'
    if len(raw.encode())>24000:
        raise RuntimeError('Compact full-root summary exceeds its 24000-byte ceiling')
    with args.output.open('x',encoding='utf-8',newline='\n') as stream:
        stream.write(raw)
    print(json.dumps(dict(output=str(args.output),bytes=len(raw.encode()),status=result['status'],
        stopped_reason=stopped,compute=result['compute'],spreads=comparisons)),flush=True)
    if stopped is not None:
        raise SystemExit(1)


if __name__=='__main__':
    main()
