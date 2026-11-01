"""Coupled same-action counterrotation test; actual GH fields, no spherical RHS.

Dense solver polynomials remain under .build. Admission is a sampled numerical
criterion, not a continuum error certificate or an automatic child diagnosis.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import time

import numpy as np
from numba import set_num_threads
from scipy.fft import set_workers
from scipy.integrate import DOP853
from scipy.optimize import brentq
from threadpoolctl import threadpool_limits

from recursive_horizons.nsc_atlas_parent_gh import ParentGH,PAIRS
from recursive_horizons.nsc_atlas_s3_hopf import build_hopf_frame_basis
from recursive_horizons.nsc_atlas_valve_counterrotation import CounterrotatingData
from recursive_horizons.nsc_atlas_gh_constraints import adm_constraints
from recursive_horizons.nsc_atlas_gh_observer import observer_readout
from recursive_horizons.nsc_atlas_gh_history import DenseGHHistory
from recursive_horizons.nsc_atlas_passive_dop853 import PassiveDOP853Quadrature
from recursive_horizons.nsc_atlas_trapped_parent import GHSurfaceGeometry
from recursive_horizons.nsc_atlas_trapped_surfaces import evaluate_surface,round_surface

ROOT=Path(__file__).resolve().parents[1]
V3=2*math.pi**2


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def project_initial_state(system,data):
    """Explicit preparation adapter; never invoke the mixed-source initial()."""
    if system.backend!='hopf' or system.kick_amplitude!=data.kick_amplitude or system.kick_source!=2:
        raise ValueError('matching Hopf counterrotation system required')
    for name,left,right in (('H',system.r,data.parent.H),('mass',system.parent.mass_squared,data.parent.mass_squared),
                            ('vacuum',system.parent.vacuum,data.parent.vacuum),('kappa',system.parent.kappa,data.parent.kappa)):
        if not math.isclose(left,right,rel_tol=2e-14):raise ValueError('parent mismatch: '+name)
    D=max(4*data.degree+system.basis.degree,2*system.basis.degree)
    startup=build_hopf_frame_basis(system.basis.degree,(math.ceil(D/4)+2,2*math.ceil((D+1)/2)))
    psi=np.empty(len(startup.points))
    for lo in range(0,len(psi),8192):psi[lo:lo+8192]=data.values(startup.points[lo:lo+8192])
    Fcoef=startup.project(psi**4)
    Fcoef[~np.all(system.basis.characters==[1,1],axis=1)]=0.
    F=startup.evaluate(Fcoef)
    if np.min(F)<=0:raise ValueError('represented initial metric is not positive')
    C=np.zeros((24,len(system.basis.labels)));C[0,0]=-math.sqrt(V3)
    for k,(a,b) in enumerate(PAIRS):
        if a==b and a:C[k]=Fcoef
    j=system.basis.labels.index((1,0,1,'sin'))
    C[21,j]=-data.A*math.pi/math.sqrt(2)
    C[23]=startup.project(data.kick_amplitude*data.f(startup.points)*F**-1.5)
    pairs=[(k,a-1,b-1) for k,(a,b) in enumerate(PAIRS) if a and b]
    velocity=np.empty((len(F),6))
    for lo in range(0,len(F),8192):
        tensor=data.tensor(startup.points[lo:lo+8192])
        for col,(_,a,b) in enumerate(pairs):
            velocity[lo:lo+8192,col]=-2*data.kick_amplitude*F[lo:lo+8192]**-.5*tensor[:,a,b]
    for col,(k,_,_) in enumerate(pairs):C[10+k]=startup.project(velocity[:,col])
    for i,Di in enumerate(system.basis.generators):
        k=10+PAIRS.index((0,i+1))
        C[k]=startup.project(-.5*startup.evaluate(system.r*(Di@Fcoef))/F)
    leakage=0.
    for k,char in enumerate(system.chars):
        wrong=~np.all(system.basis.characters==char,axis=1)
        leakage=max(leakage,float(np.linalg.norm(C[k,wrong])));C[k,wrong]=0.
    # This declared inflow-lobe geodesic defines the proper endpoint.
    system.X0=np.array([0.,2**-.5,0.,2**-.5]);system.update_point(system.X0)
    metadata=dict(preparation_id='one-mode-counterrotation-v1',chi='A X1',Pi='lambda f psi^-6',
        K_reference='lambda psi^-2 Abar',A=data.A,kick_amplitude=data.kick_amplitude,
        source_harmonic=2,geometry_degree=data.degree,GH_degree=system.basis.degree,
        startup_quadrature_shape=list(startup.quadrature_shape),min_original_psi=float(min(psi)),
        min_represented_F=float(min(F)),disallowed_character_leakage=leakage,
        observer_X=system.X0.tolist(),observer_role='Initially normal inflow-lobe geodesic; not the old mixed-source anchor.',
        freely_specified_TT=0.,input_provenance=data.provenance)
    system.initial_kick_metadata=metadata
    return np.r_[system.pack(C),system.X0,[1.,0.,0.,0.],0.,0.],metadata


def gauge_covector(system,t,g,v,dg,lapse):
    """First-jet GH constraint, independently of metric acceleration/source."""
    inverse=np.linalg.inv(g);D=np.empty((len(g),4,4,4));D[:,0]=v
    for i in range(3):
        G=system.Gamma[:,i+1,:]
        D[:,i+1]=dg[:,i]-G.T@g-g@G
    C=np.einsum('nbc,nbac->na',inverse,D)-.5*np.einsum('nbc,nabc->na',inverse,D)
    base,_=system.gauge_source(t)
    C[:,0]+=base+(system.r*(lapse*lapse-lapse**-2) if system.gauge=='lapse_controlled' else 0.)
    return C


def constraint_ratios(system,adm,h,C):
    hi=np.linalg.inv(h);Ku=np.einsum('nij,njk->nik',hi,adm['Kij'])
    trace=np.trace(Ku,axis1=1,axis2=2);K2=np.einsum('nij,nji->n',Ku,Ku)
    scale=system.r**2+abs(adm['R3'])+trace*trace+K2+2*system.parent.kappa*abs(adm['rho'])
    normM=np.sqrt(np.maximum(0.,np.einsum('ni,nij,nj->n',adm['M'],hi,adm['M'])))
    normalC=(C[:,0]-np.sum(adm['shift']*C[:,1:],axis=1))/adm['lapse']
    normC=np.sqrt(np.maximum(0.,normalC**2+np.einsum('ni,nij,nj->n',C[:,1:],hi,C[:,1:])))
    return dict(H=abs(adm['H'])/scale,M=normM/scale,GH=normC/np.sqrt(scale)),trace,K2


def diagnostic(system,t,state,*,limit=.02):
    J=system.jets(state);W=system.potential(J.q)
    adm=adm_constraints(J.g,J.v,J.dg,J.d2g,J.dv,J.q,J.p,J.dq,r=system.r,
                        kappa=system.parent.kappa,potential=W)
    h=J.g[:,1:,1:];hi=np.linalg.inv(h)
    C=gauge_covector(system,t,J.g,J.v,J.dg,adm['lapse'])
    ratios,trace,K2=constraint_ratios(system,adm,h,C)
    weights=system.basis.weights*np.sqrt(np.linalg.det(h[:system.N]))/system.r**3
    volume=float(np.sum(weights));energy=float(weights@adm['rho'][:system.N])
    pi=(J.p-np.einsum('ni,nia->na',adm['shift'],J.dq))/adm['lapse'][:,None]
    grad2=np.einsum('nij,nia,nja->n',hi,J.dq,J.dq)
    pressure=.5*np.sum(pi*pi,axis=1)-grad2/6-W
    rms=lambda a:float(np.sqrt(weights@(a[:system.N]**2)/volume))
    provider=GHSurfaceGeometry(system,state,coordinate_t=t,provenance={'preparation_id':'one-mode-counterrotation-v1'})
    surface=round_surface([1.,0.,0.,0.],math.pi/2)
    cut=evaluate_surface(surface,provider,frame_scale=system.r,quadrature_order=(14,29))
    geometry=provider(cut.points);sgrad=np.einsum('ni,nia->na',cut.outward_normal,geometry['dq'])
    js=-np.sum(geometry['Pi']*sgrad,axis=1);material=-geometry['lapse']*js
    # Neighborhood controls use actual second spatial/first temporal metric jets.
    patch=np.concatenate([evaluate_surface(round_surface([1.,0.,0.,0.],math.pi/2+d),provider,
        frame_scale=system.r,quadrature_order=(8,17)).points for d in (-.04,0.,.04)])
    pj=provider._jets(patch,second=True)
    pa=adm_constraints(**pj,r=system.r,kappa=system.parent.kappa,potential=system.potential(pj['q']))
    pc=gauge_covector(system,t,pj['g'],pj['v'],pj['dg'],pa['lapse'])
    pr,_,_=constraint_ratios(system,pa,pj['g'][:,1:,1:],pc)
    obs=observer_readout(J.g[-1],J.q[-1],J.p[-1],J.dq[-1],state[system.dof+4:system.dof+8],
                         potential=W[-1],omega_ref=system.parent.reference_frequency)
    maxima={key:float(np.max(val[:system.N])) for key,val in ratios.items()}
    patchmax={key:float(np.max(val)) for key,val in pr.items()}
    admitted=max(*maxima.values(),*patchmax.values())<=limit and float(obs['input_norm_error'])<=1e-5
    coef=system.unpack(state);ell=np.array([v[0] for v in system.basis.labels]);high=ell>=system.basis.degree-1
    power=coef[21]**2+coef[23]**2/system.parent.mass_squared
    return dict(coordinate_t=float(t),observer_proper_duration=float(state[system.dof+8]),
        observer_X=state[system.dof:system.dof+4].tolist(),observer_unit_norm_error=float(obs['input_norm_error']),
        observer_total_density=float(obs['density_total']),observer_normal_expansion=-float(trace[-1]),
        observer_normal_shear_norm=math.sqrt(max(0.,float(K2[-1]-trace[-1]**2/3))),
        proper_slice_volume=volume,total_slice_normal_matter_energy=energy,
        slice_mean_density=energy/volume,slice_integrated_isotropic_pressure=float(weights@pressure[:system.N]),
        slice_w=float(weights@pressure[:system.N])/energy,
        normal_energy_minus_fixed_vacuum=energy-system.parent.vacuum*volume,
        Hamiltonian_RMS_over_H2=rms(adm['H'])/system.r**2,
        momentum_RMS_over_H2=rms(np.sqrt(np.maximum(0.,np.einsum('ni,nij,nj->n',adm['M'],hi,adm['M']))))/system.r**2,
        global_max_relative_constraints=maxima,neck_neighborhood_max_relative_constraints=patchmax,
        relative_scale='H_ref²+|R3|+(trK)²+KijKij+2kappa|rho|, an ADM-term scale, not a spacetime-curvature certificate.',
        sampled_admitted=admitted,relative_constraint_limit=limit,admission_is_continuum_error_bound=False,
        min_lapse=float(min(adm['lapse'])),max_lapse=float(max(adm['lapse'])),
        last_two_scalar_shell_power_fraction=float(np.sum(power[high])/np.sum(power)) if np.sum(power)>0 else 0.,
        neck=dict(cut.summary(),area_radius=math.sqrt(cut.area/(4*math.pi)),
            theta_plus_over_H_range=[float(min(cut.theta_plus)/system.r),float(max(cut.theta_plus)/system.r)],
            theta_minus_over_H_range=[float(min(cut.theta_minus)/system.r),float(max(cut.theta_minus)/system.r)],
            net_material_influx=float(cut.area_weights@material),
            positive_material_influx=float(cut.area_weights@np.maximum(material,0.)),
            material_outflux=float(cut.area_weights@np.maximum(-material,0.)),
            orientation='Outward from center+X1; positive inflow enters that hemisphere. The fixed round interface is not yet a tracked receiving-child boundary.'),
        regular_child_proved=False)


def source_pins():
    package=(ROOT/'src/recursive_horizons').resolve();paths={Path(__file__).resolve()}
    for module in tuple(sys.modules.values()):
        name=getattr(module,'__file__',None)
        if name:
            path=Path(name).resolve()
            if path.suffix=='.py' and path.is_relative_to(package):paths.add(path)
    return {p.relative_to(ROOT).as_posix():sha(p) for p in sorted(paths)}


def run(args):
    folder=(ROOT/args.output).resolve()
    if not folder.is_relative_to(ROOT/'.build') or folder.exists():raise ValueError('fresh output directory under .build required')
    folder.mkdir(parents=True)
    started=time.perf_counter();set_num_threads(args.threads)
    input_path=(ROOT/args.initial_record).resolve();input_sha=sha(input_path)
    result=dict(schema='nsc-atlas-parent-gh-v3',degree=args.degree,epsilon=1.,
        prescribed_GH_gauge='lapse_controlled',spatial_transform='hopf',
        requested_proper_endpoint=args.proper_endpoint,threads=args.threads,budget_seconds=args.budget,
        initial_record=args.initial_record.as_posix(),initial_record_sha256=input_sha,
        snapshots=[],dense_history_chunks=[],regular_child_proved=False,status='initializing',
        temporal_tolerances=dict(rtol=args.rtol,atol=args.atol,max_step=args.max_step))
    saved=[];history=DenseGHHistory();chunk=0
    def save_record():
        result['elapsed_seconds']=time.perf_counter()-started
        tmp=folder/'live.tmp';tmp.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
        tmp.replace(folder/'run.json')
    with threadpool_limits(limits=args.threads),set_workers(args.threads):
        data=CounterrotatingData.from_record(input_path)
        system=ParentGH(args.degree,1.,gauge='lapse_controlled',backend='hopf',kick_amplitude=data.kick_amplitude,kick_source=2)
        state,prep=project_initial_state(system,data)
        pins=source_pins();commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
        bound=True
        for name,digest in pins.items():
            old=subprocess.run(['git','show',commit+':'+name],cwd=ROOT,capture_output=True)
            bound=bound and old.returncode==0 and hashlib.sha256(old.stdout).hexdigest()==digest
        result.update(initial_preparation=prep,source_sha256=pins,producer_commit=commit,source_bound_to_commit=bool(bound),
            model_inputs=dict(H_reference=system.r,vacuum=system.parent.vacuum,scalar_mass_squared=system.parent.mass_squared,
                volume_kinetic=system.K,sstar=system.sstar,kick_amplitude=data.kick_amplitude,kick_source=2))
        meta=dict(degree=args.degree,epsilon=1.,gauge=system.gauge,backend='hopf',
            quadrature_order=list(system.basis.quadrature_shape[:2]),kick_amplitude=data.kick_amplitude,kick_source=2,
            source_sha256=pins,producer_commit=commit,preparation_id='one-mode-counterrotation-v1')
        def flush():
            nonlocal history,chunk
            if not history.bounds:return
            path=folder/f'dense-{chunk:04d}.npz';history.save(path,meta)
            result['dense_history_chunks'].append(dict(path=path.relative_to(ROOT).as_posix(),sha256=sha(path),
                coordinate_bounds=[history.bounds[0][0],history.bounds[-1][1]],intervals=len(history.bounds)))
            history=DenseGHHistory();chunk+=1
        initial=diagnostic(system,0.,state);initial['integrated_excitation_work']=0.;initial['excitation_work_residual']=0.
        result['snapshots'].append(initial);saved.append((0.,state.copy()));save_record()
        print(json.dumps(dict(event='initial',degree=args.degree,H_res=initial['Hamiltonian_RMS_over_H2'],
            M_res=initial['momentum_RMS_over_H2'],admitted=initial['sampled_admitted'],elapsed=time.perf_counter()-started)),flush=True)
        if args.initial_only or not initial['sampled_admitted']:
            result['status']='initial_only' if args.initial_only else 'initial_not_admitted'
        else:
            system.started=started;system.budget=args.budget
            ledger=PassiveDOP853Quadrature();system.passive_work=ledger
            solver=DOP853(system.rhs,0.,state,t_bound=8/system.r,rtol=args.rtol,atol=args.atol,
                          first_step=.002,max_step=args.max_step)
            ledger.initialize(solver);normalization=system.parent.vacuum*V3/system.r**3
            targets=sorted(set([.02]+[i/10 for i in range(1,math.ceil(10*args.proper_endpoint)+1)]+[args.proper_endpoint]))
            targets=[v for v in targets if v<=args.proper_endpoint];steps=0;last_admitted=0.;status='incomplete';reason=None
            heartbeat=time.perf_counter()
            try:
                while solver.status=='running' and targets:
                    if time.perf_counter()-started>args.budget:raise TimeoutError('declared formation run budget exhausted')
                    message=solver.step();steps+=1
                    if solver.status=='failed':raise RuntimeError(message or 'DOP853 failed')
                    ledger.accepted_step(solver);dense=solver.dense_output()
                    history.record(solver.t_old,solver.t,dense)
                    if len(history.bounds)>=32:flush()
                    if time.perf_counter()-heartbeat>=30:
                        result['current_progress']=dict(coordinate_t=float(solver.t),proper=float(solver.y[system.dof+8]),
                            steps=steps,rhs_calls=system.calls)
                        save_record();print(json.dumps(dict(event='progress',degree=args.degree,
                            **result['current_progress'],elapsed=result['elapsed_seconds'])),flush=True)
                        heartbeat=time.perf_counter()
                    while targets and solver.y[system.dof+8]>=targets[0]:
                        target=targets.pop(0)
                        t=brentq(lambda at:dense(at)[system.dof+8]-target,solver.t_old,solver.t,xtol=2e-12)
                        point=dense(t);diag=diagnostic(system,t,point)
                        work=float(ledger.dense_value(solver,t))*normalization
                        diag['integrated_excitation_work']=work
                        diag['excitation_work_residual']=diag['normal_energy_minus_fixed_vacuum']-initial['normal_energy_minus_fixed_vacuum']-work
                        result['snapshots'].append(diag);saved.append((t,point.copy()))
                        if diag['sampled_admitted']:last_admitted=t
                        result.update(steps=steps,rhs_calls=system.calls,last_admitted_coordinate_t=last_admitted)
                        save_record()
                        print(json.dumps(dict(event='checkpoint',degree=args.degree,proper=diag['observer_proper_duration'],
                            area_radius=diag['neck']['area_radius'],theta_plus=diag['neck']['theta_plus_over_H_range'],
                            theta_minus=diag['neck']['theta_minus_over_H_range'],admitted=diag['sampled_admitted'],
                            work_residual=diag['excitation_work_residual'],elapsed=result['elapsed_seconds'])),flush=True)
                        if not diag['sampled_admitted']:raise RuntimeError('sampled constraint/gauge/observer admission ended; no physical endpoint is inferred')
                if not targets:status='endpoint_reached'
                elif solver.status=='finished':reason='coordinate allowance ended before proper target'
            except (TimeoutError,ValueError,RuntimeError,FloatingPointError) as exc:reason=str(exc)
            flush();result.update(status=status,stop_reason=reason,steps=steps,rhs_calls=system.calls,
                actual_last_accepted_coordinate_t=float(solver.t),actual_last_accepted_proper=float(solver.y[system.dof+8]),
                last_admitted_coordinate_t=last_admitted,invalid_solver_stages=system.invalid_stages,passive_work_valid=system.work_valid)
            # Keep the sparse cache aligned with its diagnostic snapshots.
            # The actual accepted endpoint may lie beyond the last controlled
            # event; save that tail separately without promoting its admission.
            terminal=folder/'accepted-endpoint.npz'
            np.savez_compressed(terminal,state=solver.y,coordinate_t=float(solver.t),metadata=json.dumps(meta,sort_keys=True))
            result['accepted_endpoint_cache']=dict(path=terminal.relative_to(ROOT).as_posix(),sha256=sha(terminal),
                admitted_at_this_exact_time=bool(result['snapshots'][-1]['coordinate_t']==solver.t and result['snapshots'][-1]['sampled_admitted']))
        cache=folder/'sparse.npz';np.savez_compressed(cache,times=np.array([v[0] for v in saved]),states=np.array([v[1] for v in saved]),metadata=json.dumps(meta,sort_keys=True))
        result.update(state_cache=cache.relative_to(ROOT).as_posix(),state_cache_sha256=sha(cache),
            producer_source_unchanged_during_run=pins=={name:sha(ROOT/name) for name in pins},
            initial_record_unchanged=input_sha==sha(input_path))
        save_record();print(json.dumps(dict(event='terminal',degree=args.degree,status=result['status'],stop_reason=result.get('stop_reason'),
            elapsed=result['elapsed_seconds'],output=str(folder/'run.json'))),flush=True)
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--degree',type=int,required=True);p.add_argument('--threads',type=int,default=4)
    p.add_argument('--proper-endpoint',type=float,default=2.);p.add_argument('--budget',type=float,default=3600.)
    p.add_argument('--rtol',type=float,default=3e-9);p.add_argument('--atol',type=float,default=3e-11)
    p.add_argument('--max-step',type=float,default=.025);p.add_argument('--initial-only',action='store_true')
    p.add_argument('--initial-record',type=Path,default=Path('results/summaries/nsc-atlas-valve-counterrotation-v1.json'))
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    if args.threads<1 or args.threads>8 or args.proper_endpoint<=0 or args.budget<=0:raise ValueError('positive endpoint/budget and threads1..8 required')
    run(args)


if __name__=='__main__':main()
