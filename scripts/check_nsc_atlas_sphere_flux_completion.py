"""Check a smooth source completion and its actual frozen-slice shape operator.

The old directional background is unchanged. The new higher-dimensional field
content matches its reduced scalar action; no inherited state, SI ruler or
child is inferred from this calculation.
"""
import argparse
from dataclasses import asdict
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import time

import numpy as np
from scipy.linalg import eigh
from scipy.special import roots_chebyu
from threadpoolctl import threadpool_limits

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from continue_nsc_atlas_valve_directional import authentic_run
from recursive_horizons.nsc_atlas_valve_spherical import SphericalValve,gegenbauer_jets
from recursive_horizons.nsc_atlas_sphere_flux_completion import (
    action_reduction_record,flux_stress_check,shape_mass_squared,shape_threshold,integer_flux_embedding,
)


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def source_pins():
    package=(ROOT/'src/recursive_horizons').resolve();paths={Path(__file__).resolve(),
        ROOT/'scripts/continue_nsc_atlas_valve_directional.py'}
    for module in tuple(sys.modules.values()):
        name=getattr(module,'__file__',None)
        if name:
            p=Path(name).resolve()
            if p.suffix=='.py' and p.is_relative_to(package):paths.add(p)
    return {p.relative_to(ROOT).as_posix():sha(p) for p in sorted(paths)}


def matrices(system,state,degree,order):
    u,w=roots_chebyu(order);z=system.pieces(state,u)
    Y,Yu,_=gegenbauer_jets(u,degree+1);Y/=math.sqrt(2*math.pi**2);Yu/=math.sqrt(2*math.pi**2)
    measure=4*math.pi*w*z['q']*z['b']**2
    local_s=system.sstar+z['y'][4]/math.sqrt(84)
    mass=shape_mass_squared(local_s)
    M=Y.T@(measure[:,None]*Y)
    stiffness=Yu.T@((measure*(1-u*u)/z['q']**2)[:,None]*Yu)+Y.T@((measure*mass)[:,None]*Y)
    return .5*(stiffness+stiffness.T),.5*(M+M.T),dict(
        sampled_mass_range=[float(min(mass)),float(max(mass))],
        proper_volume_fraction_negative_mass=float(measure@(mass<0)/sum(measure)),
        constant_trial_Rayleigh=float(measure@mass/sum(measure)),
        Gram_condition=float(np.linalg.cond(M)))


def pointwise_trial(system,state,coefficients):
    u=np.linspace(-1.,1.,2001);z=system.pieces(state,u)
    Y,Yu,Yuu=gegenbauer_jets(u,len(coefficients));normal=math.sqrt(2*math.pi**2)
    phi=Y@coefficients/normal;du=Yu@coefficients/normal;duu=Yuu@coefficients/normal
    if np.sum(phi)<0:phi=-phi;du=-du;duu=-duu
    v=1-u*u;bu=z['yu'][0];cu=z['yu'][1];qu=bu-2*u*z['y'][1]+v*cu
    drift=(v*(2*bu/z['b']-qu/z['q'])-3*u)/z['q']**2
    lap=v*duu/z['q']**2+drift*du
    mass=shape_mass_squared(system.sstar+z['y'][4]/math.sqrt(84))
    positive=bool(min(phi)>1e-10*max(abs(phi)))
    ratio=(-lap+mass*phi)/phi if positive else None
    return dict(sampled_trial_positive=positive,sampled_trial_minimum=float(min(phi)),
        sampled_Aphi_over_phi_range=[float(min(ratio)),float(max(ratio))] if positive else None,
        scope='Point samples including the smooth poles, not a whole-interval Barta certificate.')


def calculation(args):
    started=time.perf_counter();before=source_pins();histories={};origins=[]
    for path in args.run:
        record,_,_,origin=authentic_run(path);origins.append(origin)
        with np.load(ROOT/record['state_cache'],allow_pickle=False) as raw:
            states=raw['states'].copy();times=raw['times'].copy()
        for row,state,t in zip(record['snapshots'],states,times):
            if not row['sampled_admitted']:continue
            key=round(float(t),10)
            if key in histories and not np.array_equal(histories[key][0],state):
                raise ValueError('duplicate time has different inherited state')
            histories[key]=(state,record['modes'])
    if not histories:raise ValueError('admitted source history required')
    modes={v[1] for v in histories.values()}
    if len(modes)!=1:raise ValueError('one fixed background resolution required')
    system=SphericalValve(modes.pop(),reflection_symmetry=False)
    rows=[];vectors=[]
    with threadpool_limits(limits=4):
        for t,(state,_) in sorted(histories.items()):
            if time.perf_counter()-started>60:raise TimeoutError('one60s frozen-shape operator allowance exhausted')
            K,M,info=matrices(system,state,48,601)
            low,vec=eigh(K,M,subset_by_index=(0,3),check_finite=True)
            low32=eigh(K[:33,:33],M[:33,:33],subset_by_index=(0,3),eigvals_only=True)
            Ka,Ma,_=matrices(system,state,48,401)
            check=eigh(Ka,Ma,subset_by_index=(0,3),eigvals_only=True)
            trial=vec[:,0];rayleigh=float(trial@Ka@trial/(trial@Ma@trial))
            if low[0]>low32[0]+1e-10:raise RuntimeError('nested Ritz variational ordering violated')
            rows.append(dict(proper_time=t,background_modes=system.modes,
                degrees=[32,48],quadrature_orders=[401,601],
                first_four_eigenvalues=low.tolist(),first_four_eigenvalues_degree32=low32.tolist(),
                first_four_eigenvalues_quadrature401=check.tolist(),
                ground_eigenvalue_over_H2=float(low[0]/system.H**2),
                degree_ground_difference=float(abs(low[0]-low32[0])),
                quadrature_ground_difference=float(abs(low[0]-check[0])),
                independent_trial_Rayleigh=rayleigh,**info,
                positive_trial_samples=pointwise_trial(system,state,trial)))
            vectors.append(trial)
    if before!=source_pins():raise RuntimeError('probe sources changed during calculation')
    folder=(ROOT/args.raw_directory).resolve()
    if not folder.is_relative_to(ROOT/'.build') or folder.exists():raise ValueError('fresh raw directory under .build required')
    folder.mkdir(parents=True);raw=folder/'ground-trials.npz'
    np.savez_compressed(raw,times=np.array([r['proper_time'] for r in rows]),coefficients=np.array(vectors))
    stress=flux_stress_check();action=asdict(action_reduction_record())
    action=json.loads(json.dumps(action,default=lambda v:str(v) if isinstance(v,Fraction) else v))
    minimum=min(rows,key=lambda r:r['first_four_eigenvalues'][0])
    embeddings=[]
    for name,unit,vector,step in (
        ('diagonal_six_cycle',Fraction(1,96),(5,)*6,(1,)*6),
        ('single_pair_cycle',Fraction(1,16),(5,0,0,0,0,0),(1,0,0,0,0,0))):
        embedding=integer_flux_embedding(unit,vector)
        payload=json.loads(json.dumps(asdict(embedding),default=str))
        payload.update(name=name,inherited_step_norm_squared=str(unit*sum(v*v for v in step)),
            initial_shape_accelerations=[float(v)*math.exp(-24*system.sstar) for v in embedding.equal_radius_shape_acceleration_coefficients])
        embeddings.append(payload)
    return dict(schema='nsc-atlas-sphere-flux-completion-v1',
        scientific_owner='docs/nsc-atlas-valve-formation.md',
        action=action,exact_unit_flux=dict(norm_squared=str(stress.norm_squared),
            internal_contraction_diagonal=str(stress.internal_contraction[0][0]),
            internal_stress_zero=all(v==0 for row in stress.internal_stress for v in row),
            external_metric_coefficient=str(stress.external_metric_coefficient)),
        shape=dict(canonical_kinetic_prefactor=3,local_mass_squared='-5exp(-14s)/2+25exp(-24s)/8',
            threshold_s=shape_threshold(),linear_channels=3,
            spatial_operator='A=-w^-1 d_u[p d_u]+m²(s); w=4pi sqrt(1-u²)q b²,p=w(1-u²)/q²',
            scope='Full radial parity; smooth poles. Higher external angular sectors add a nonnegative term. Frozen-slice operator, not physical seed/state selection or nonlinear growth.'),
        inheritance_embedding_comparison=dict(examples=embeddings,
            equal_reduced_T=all(x['reduced_flux_energy']=='25/32' for x in embeddings),
            equal_inherited_step_norm=all(x['inherited_step_norm_squared']=='1/16' for x in embeddings),
            missing_geometric_datum='Which internal six-cycle or combination is wrapped by the inherited unit charge.',
            scope='Two conditional microscopic embeddings of the same effective flux unit and energy, with different initial shape forces. The actual trajectory probed here belongs only to the diagonal equal-flux completion; it is not reused as a single-pair history.'),
        origins=origins,slices=rows,
        verdict=dict(minimum_retained_ground_eigenvalue=minimum['first_four_eigenvalues'][0],
            at_proper_time=minimum['proper_time'],negative_trial_times=[r['proper_time'] for r in rows if r['independent_trial_Rayleigh']<0],
            full_continuum_stability_proved=False,physical_shape_seed_selected=False,
            child_proved=False,SI_scale_selected=False,
            scope='Positive Ritz minima alone cannot exclude omitted negative directions; sampled positive trials are not certified Barta bounds. Exact zero shapes remain zero.'),
        provenance=dict(source_sha256=before,source_unchanged=True,
            producer_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
            raw_cache=raw.relative_to(ROOT).as_posix(),raw_cache_sha256=sha(raw)),
        compute=dict(threads=4,budget_seconds=60,elapsed_seconds=time.perf_counter()-started))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run',type=Path,action='append',required=True)
    parser.add_argument('--raw-directory',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();path=(ROOT/args.output).resolve()
    if not path.is_relative_to(ROOT) or path.exists():raise ValueError('fresh in-repository output required')
    result=calculation(args);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(output=str(path),verdict=result['verdict'],compute=result['compute']),indent=2))


if __name__=='__main__':main()
