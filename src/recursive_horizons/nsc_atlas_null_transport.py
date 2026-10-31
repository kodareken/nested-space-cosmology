"""Null geodesics, parallel screens and regular Sachs Jacobi transport.

Components use (partial_t,E_i), E_i=r*A_i*X on unit S3. Gamma^a_bc has
the derivative index b: nabla_{E_b}E_c=Gamma^a_bc E_a. Jacobi D and its
affine derivative P are evolved directly; no Riccati expansion enters
the adaptive state or stops integration at a focal point.
"""
from dataclasses import dataclass
import math
import time

import numpy as np
from scipy.integrate import DOP853


STATE_SIZE=27


def _real(value,name,shape):
    raw=np.asarray(value)
    if raw.dtype.kind not in 'iuf' or raw.shape!=shape or not np.all(np.isfinite(raw)):
        raise ValueError(name+' must be finite real shape '+str(shape))
    return np.array(raw,dtype=float,copy=True)


def pack_null_state(initial):
    """Copy a bundle while retaining the caller's exact ray ordering."""
    X=np.asarray(initial['X'])
    if X.ndim!=2 or X.shape[1]!=4 or not len(X):raise ValueError('X must have shape(rays,4)')
    n=len(X);out=np.empty((n,STATE_SIZE));out[:,:4]=_real(X,'X',(n,4))
    if np.any(abs(np.sum(out[:,:4]**2,axis=1)-1)>2e-10):raise ValueError('initial X must lie on unit S3')
    out[:,4:8]=_real(initial['k'],'k',(n,4))
    if np.any(out[:,4]<=0):raise ValueError('future null rays require positive initial k0')
    screen=np.asarray(initial['screen'])
    if screen.shape==(n,2,3):
        spatial=_real(screen,'spatial screen',(n,2,3));screen=np.zeros((n,2,4));screen[:,:,1:]=spatial
    out[:,8:16]=_real(screen,'screen',(n,2,4)).reshape(n,8)
    out[:,16:20]=_real(initial.get('D',np.broadcast_to(np.eye(2),(n,2,2))),'D',(n,2,2)).reshape(n,4)
    if 'P' not in initial:raise ValueError('initial affine Jacobi derivative P is required')
    out[:,20:24]=_real(initial['P'],'P',(n,2,2)).reshape(n,4)
    for name,column in (('affine_lambda',24),('integrated_Rkk',25),('integrated_source_Rkk',26)):
        value=np.asarray(initial.get(name,np.zeros(n)))
        if value.ndim==0:value=np.full(n,float(value))
        out[:,column]=_real(value,name,(n,))
    return out


def unpack_null_state(packed):
    y=np.asarray(packed)
    if y.ndim!=2 or y.shape[1]!=STATE_SIZE:raise ValueError('packed bundle must have shape(rays,27)')
    return dict(X=y[:,:4],k=y[:,4:8],screen=y[:,8:16].reshape(-1,2,4),
        D=y[:,16:20].reshape(-1,2,2),P=y[:,20:24].reshape(-1,2,2),
        affine_lambda=y[:,24],integrated_Rkk=y[:,25],integrated_source_Rkk=y[:,26])


def _frame(X):
    x,y,z,w=X.T
    return np.stack((np.c_[w,z,-y,-x],np.c_[-z,w,x,-y],np.c_[y,-x,w,-z]),axis=1)


def _geometry(sampler,t,X):
    norm=np.linalg.norm(X,axis=1)
    if np.any(~np.isfinite(norm)) or np.any(norm==0):raise ValueError('invalid null position direction')
    # Smooth directional extension for intermediate RK stages only. Stored
    # X, k, screens and optical data are never projected or renormalized.
    sample=sampler(float(t),X/norm[:,None]);n=len(X);out={}
    for key,shape in (('g',(n,4,4)),('Gamma',(n,4,4,4)),('R_cov',(n,4,4,4,4)),('Ricci',(n,4,4))):
        alternate=sample.get('fullGamma') if key=='Gamma' else (sample.get('Riemann') if key=='R_cov' else None)
        value=sample.get(key,alternate)
        if value is None:raise ValueError('sampler lacks '+key)
        array=np.asarray(value)
        if array.shape==shape[1:]:array=np.broadcast_to(array,shape)
        out[key]=_real(array,'sampled '+key,shape)
    if 'matter_Ricci' in sample:
        source=np.asarray(sample['matter_Ricci'])
        if source.shape==(4,4):source=np.broadcast_to(source,(n,4,4))
        out['matter_Ricci']=_real(source,'sampled matter_Ricci',(n,4,4))
        out['matter_Ricci_role']='provided_matter_Ricci'
    else:
        # Analytic geometry-only controls may omit a matter model. This
        # explicit fallback is not an independent Einstein-source check.
        out['matter_Ricci']=out['Ricci']
        out['matter_Ricci_role']='geometric_Ricci_fallback_control'
    if np.any(np.linalg.eigvalsh(out['g'][:,1:,1:])<=0):raise ValueError('sampled spatial metric is not positive')
    inverse=np.linalg.inv(out['g'])
    if not np.all(np.isfinite(inverse)) or np.any(inverse[:,0,0]>=0):raise ValueError('sampled lapse is not real positive')
    out['lapse']=1/np.sqrt(-inverse[:,0,0])
    return out


def optical_diagnostics(state,geometry,*,initial_wronskian=None,D_reference_scale=None,caustic_tolerance=1e-9):
    """Expansion is undefined near singular D; D and P themselves remain valid."""
    y=unpack_null_state(state);g=geometry['g'];k=y['k'];e=y['screen'];D=y['D'];P=y['P'];n=len(k)
    energy=geometry['lapse']*k[:,0]
    if np.any(energy<=0):raise ValueError('accepted ray is not future directed')
    null=np.einsum('nab,na,nb->n',g,k,k)
    gram=np.einsum('nab,nAa,nBb->nAB',g,e,e)
    ek=np.einsum('nab,nAa,nb->nA',g,e,k)
    tidal=np.einsum('nabcd,nAa,nb,nBc,nd->nAB',geometry['R_cov'],e,k,e,k,optimize=True)
    Rkk=np.einsum('nab,na,nb->n',geometry['Ricci'],k,k)
    source_Rkk=np.einsum('nab,na,nb->n',geometry['matter_Ricci'],k,k)
    singular=np.linalg.svd(D,compute_uv=False)
    reference=np.ones(n) if D_reference_scale is None else np.asarray(D_reference_scale,dtype=float)
    defined=singular[:,-1]>caustic_tolerance*reference
    theta=np.full(n,np.nan)
    if np.any(defined):
        # P*D^-1, computed through an actual solve; no pseudoinverse or
        # regulated expansion is substituted at a caustic.
        B=np.linalg.solve(D[defined].swapaxes(1,2),P[defined].swapaxes(1,2)).swapaxes(1,2)
        theta[defined]=np.trace(B,axis1=1,axis2=2)
    wronskian=D.swapaxes(1,2)@P-P.swapaxes(1,2)@D
    initial=np.zeros_like(wronskian) if initial_wronskian is None else initial_wronskian
    return dict(theta=theta,theta_defined=defined,normal_energy=energy,
        theta_over_normal_energy=theta/energy,detD=np.linalg.det(D),
        min_singular_value=singular[:,-1],null_norm_over_normal_energy_squared=null/energy**2,
        screen_gram_error=gram-np.eye(2),screen_k_inner_over_normal_energy=ek/energy[:,None],
        X_norm_error=np.linalg.norm(y['X'],axis=1)-1,
        wronskian=wronskian,wronskian_drift=wronskian-initial,
        tidal=tidal,Ricci_kk=Rkk,source_Rkk=source_Rkk,Ricci_minus_source_Rkk=Rkk-source_Rkk,
        matter_Ricci_role=geometry['matter_Ricci_role'],
        tidal_trace_minus_Ricci_kk=np.trace(tidal,axis1=1,axis2=2)-Rkk,
        tidal_antisymmetry=tidal-tidal.swapaxes(1,2))


@dataclass(frozen=True)
class NullBundleResult:
    times: np.ndarray
    packed_states: np.ndarray
    diagnostics: tuple
    status: str
    reason: str|None
    requested_interval: tuple
    effective_interval: tuple
    steps: int
    rhs_calls: int
    invalid_trials: int
    elapsed_seconds: float
    sampled_caustic_rays: np.ndarray

    @property
    def states(self):return [unpack_null_state(y) for y in self.packed_states]

    @property
    def final(self):return unpack_null_state(self.packed_states[-1])

    @property
    def caustic_encountered(self):return bool(np.any(self.sampled_caustic_rays))

    def summary(self):
        points=[]
        for t,y,d in zip(self.times,self.states,self.diagnostics):
            defined=d['theta_defined'];theta=d['theta'][defined]
            unit_theta=d['theta_over_normal_energy'][defined]
            points.append(dict(coordinate_t=float(t),theta_defined_rays=int(np.sum(defined)),
                theta_undefined_rays=int(np.sum(~defined)),
                caustic_rays_so_far=int(np.sum(d['caustic_rays_so_far'])),
                caustic_prefix_checked_through_t=float(d['caustic_prefix_checked_through_t']),
                theta_min=float(np.min(theta)) if len(theta) else None,
                theta_max=float(np.max(theta)) if len(theta) else None,
                unit_normalized_theta_min=float(np.min(unit_theta)) if len(unit_theta) else None,
                unit_normalized_theta_max=float(np.max(unit_theta)) if len(unit_theta) else None,
                determinant_range=[float(np.min(d['detD'])),float(np.max(d['detD']))],
                smallest_Jacobi_singular_value=float(np.min(d['min_singular_value'])),
                max_null_norm_over_normal_energy_squared=float(np.max(abs(d['null_norm_over_normal_energy_squared']))),
                max_screen_gram_error=float(np.max(abs(d['screen_gram_error']))),
                max_screen_k_inner_over_normal_energy=float(np.max(abs(d['screen_k_inner_over_normal_energy']))),
                max_X_norm_error=float(np.max(abs(d['X_norm_error']))),
                max_wronskian_drift=float(np.max(abs(d['wronskian_drift']))),
                max_tidal_trace_Ricci_gap=float(np.max(abs(d['tidal_trace_minus_Ricci_kk']))),
                matter_Ricci_role=d['matter_Ricci_role'],
                source_Rkk_range=[float(np.min(d['source_Rkk'])),float(np.max(d['source_Rkk']))],
                max_Ricci_minus_source_Rkk=float(np.max(abs(d['Ricci_minus_source_Rkk']))),
                affine_range=[float(np.min(y['affine_lambda'])),float(np.max(y['affine_lambda']))],
                integrated_Rkk_range=[float(np.min(y['integrated_Rkk'])),float(np.max(y['integrated_Rkk']))],
                integrated_source_Rkk_range=[float(np.min(y['integrated_source_Rkk'])),float(np.max(y['integrated_source_Rkk']))],
                integrated_Ricci_minus_source_range=[float(np.min(y['integrated_Rkk']-y['integrated_source_Rkk'])),
                                                     float(np.max(y['integrated_Rkk']-y['integrated_source_Rkk']))]))
        return dict(status=self.status,reason=self.reason,rays=self.packed_states.shape[1],
            requested_interval=list(self.requested_interval),effective_interval=list(self.effective_interval),
            steps=self.steps,rhs_calls=self.rhs_calls,invalid_trials=self.invalid_trials,
            elapsed_seconds=self.elapsed_seconds,checkpoints=points,
            sampled_caustic_rays=int(np.sum(self.sampled_caustic_rays)),
            caustic_detection_scope='Near-singular accepted/checkpoint Jacobi matrices and determinant sign changes only; no certificate excluding caustics between samples.',
            scope='Coordinate-time null and Jacobi transport on supplied data. A caustic or data endpoint is not a spacetime singularity or child continuation.')


def transport_null_bundle(sampler,initial,*,coordinate_interval,frame_scale,checkpoints=(),
        rtol=2e-9,atol=2e-11,max_step=.03,budget_seconds=300.,data_end=None,
        caustic_tolerance=1e-9,initial_constraint_tolerance=2e-7,progress=False,dense_sink=None):
    """Evolve a future bundle through caustics within a declared data domain.

    Optional dense_sink(t_old,t,dense) receives the actual DOP853 dense
    polynomial once per accepted step. Its extra dense-stage RHS calls do
    not enter adaptive state/error control. Checkpoints reuse that same
    dense output; no sparse-state interpolation is introduced.
    """
    r=float(frame_scale);t0,requested_end=map(float,coordinate_interval)
    for value in (r,rtol,atol,max_step,budget_seconds,caustic_tolerance,initial_constraint_tolerance):
        if not math.isfinite(float(value)) or value<=0:raise ValueError('positive finite transport controls required')
    if not math.isfinite(t0+requested_end) or requested_end<=t0:
        raise ValueError('future increasing coordinate interval required')
    end=requested_end
    if data_end is not None:
        data_end=float(data_end)
        if not math.isfinite(data_end) or data_end<=t0:raise ValueError('data endpoint must lie to the future')
        end=min(end,data_end)
    requested=np.asarray(checkpoints,dtype=float)
    if requested.ndim!=1 or not np.all(np.isfinite(requested)) or np.any(requested<t0) or np.any(requested>requested_end):
        raise ValueError('checkpoints must lie inside the requested interval')
    targets=sorted(set(map(float,requested[requested>t0]))|{end})
    targets=[t for t in targets if t<=end]
    y0=pack_null_state(initial);rays=len(y0);started=time.perf_counter();calls=0;invalid=0;failure=None
    D0=y0[:,16:20].reshape(-1,2,2);P0=y0[:,20:24].reshape(-1,2,2)
    W0=D0.swapaxes(1,2)@P0-P0.swapaxes(1,2)@D0
    Dscale=np.linalg.svd(D0,compute_uv=False)[:,0];Dscale=np.where(Dscale>0,Dscale,1.)
    initial_geometry=_geometry(sampler,t0,y0[:,:4])
    source_role=initial_geometry['matter_Ricci_role']
    first=optical_diagnostics(y0,initial_geometry,initial_wronskian=W0,D_reference_scale=Dscale,caustic_tolerance=caustic_tolerance)
    if max(float(np.max(abs(first[name]))) for name in ('null_norm_over_normal_energy_squared','screen_gram_error','screen_k_inner_over_normal_energy'))>initial_constraint_tolerance:
        raise ValueError('initial null/screen normalization is inconsistent with the sampled metric')
    caustic=~first['theta_defined'];previous_det=first['detD'];times=[t0];states=[y0.copy()];diagnostics=[first]
    first['caustic_rays_so_far']=caustic.copy()
    first['caustic_prefix_checked_through_t']=t0

    def rhs(t,flat):
        nonlocal calls,invalid
        if time.perf_counter()-started>budget_seconds:raise TimeoutError('null transport compute budget exceeded')
        calls+=1;y=flat.reshape(rays,STATE_SIZE)
        if not np.all(np.isfinite(y)) or np.any(y[:,4]<=0) or np.any(np.linalg.norm(y[:,:4],axis=1)==0):
            invalid+=1;return np.full_like(flat,np.nan)
        data=unpack_null_state(y);k=data['k'];e=data['screen'];geometry=_geometry(sampler,t,data['X'])
        if geometry['matter_Ricci_role']!=source_role:
            raise ValueError('matter_Ricci source role changed during null transport')
        inverse_k0=1/k[:,0];Gamma=geometry['Gamma']
        result=np.empty_like(y)
        result[:,:4]=r*np.einsum('ni,nia->na',k[:,1:],_frame(data['X']))*inverse_k0[:,None]
        result[:,4:8]=-np.einsum('nabc,nb,nc->na',Gamma,k,k)*inverse_k0[:,None]
        result[:,8:16]=(-np.einsum('nabc,nb,nAc->nAa',Gamma,k,e)*inverse_k0[:,None,None]).reshape(rays,8)
        tidal=np.einsum('nabcd,nAa,nb,nBc,nd->nAB',geometry['R_cov'],e,k,e,k,optimize=True)
        result[:,16:20]=(data['P']*inverse_k0[:,None,None]).reshape(rays,4)
        result[:,20:24]=(-np.einsum('nAB,nBC->nAC',tidal,data['D'])*inverse_k0[:,None,None]).reshape(rays,4)
        result[:,24]=inverse_k0
        result[:,25]=np.einsum('nab,na,nb->n',geometry['Ricci'],k,k)*inverse_k0
        result[:,26]=np.einsum('nab,na,nb->n',geometry['matter_Ricci'],k,k)*inverse_k0
        return result.ravel()

    solver=None;steps=0
    try:
        solver=DOP853(rhs,t0,y0.ravel(),t_bound=end,rtol=rtol,atol=atol,max_step=max_step)
        while solver.status=='running':
            message=solver.step();steps+=1
            if solver.status=='failed':failure=message or 'null integrator failed';break
            accepted=solver.y.reshape(rays,STATE_SIZE);D=accepted[:,16:20].reshape(-1,2,2)
            singular=np.linalg.svd(D,compute_uv=False)[:,-1];determinant=np.linalg.det(D)
            caustic|=(singular<=caustic_tolerance*Dscale)|(previous_det*determinant<0)
            previous_det=determinant
            dense=None
            if dense_sink is not None:
                dense=solver.dense_output()
                dense_sink(solver.t_old,solver.t,dense)
            if targets and solver.t>=targets[0]:
                if dense is None:dense=solver.dense_output()
                while targets and solver.t>=targets[0]:
                    t=targets.pop(0);snapshot=dense(t).reshape(rays,STATE_SIZE)
                    geometry=_geometry(sampler,t,snapshot[:,:4])
                    check=optical_diagnostics(snapshot,geometry,initial_wronskian=W0,D_reference_scale=Dscale,caustic_tolerance=caustic_tolerance)
                    caustic|=~check['theta_defined']
                    check['caustic_rays_so_far']=caustic.copy()
                    # Accepted-step monitoring can conservatively include
                    # an event later than this interpolated checkpoint.
                    check['caustic_prefix_checked_through_t']=float(solver.t)
                    times.append(t);states.append(snapshot.copy());diagnostics.append(check)
                    if progress:print(f'Null bundle: t={t:.7g}/{end:.7g}, rays={rays}',flush=True)
    except (TimeoutError,ValueError,FloatingPointError,np.linalg.LinAlgError) as error:
        failure=str(error)
    complete=solver is not None and solver.status=='finished' and failure is None and not targets
    if not complete and solver is not None and solver.t>times[-1] and np.all(np.isfinite(solver.y)):
        snapshot=solver.y.reshape(rays,STATE_SIZE)
        try:
            geometry=_geometry(sampler,solver.t,snapshot[:,:4])
            check=optical_diagnostics(snapshot,geometry,initial_wronskian=W0,D_reference_scale=Dscale,caustic_tolerance=caustic_tolerance)
            caustic|=~check['theta_defined']
            check['caustic_rays_so_far']=caustic.copy()
            check['caustic_prefix_checked_through_t']=float(solver.t)
            times.append(float(solver.t));states.append(snapshot.copy());diagnostics.append(check)
        except (ValueError,FloatingPointError,np.linalg.LinAlgError):pass
    status=('data_end_reached' if end<requested_end else 'endpoint_reached') if complete else 'incomplete'
    reason=('Declared data ended before the requested coordinate endpoint; no physical singularity is inferred.' if complete and end<requested_end else failure)
    return NullBundleResult(np.array(times),np.array(states),tuple(diagnostics),status,reason,
        (t0,requested_end),(t0,float(times[-1])),steps,calls,invalid,time.perf_counter()-started,caustic)
