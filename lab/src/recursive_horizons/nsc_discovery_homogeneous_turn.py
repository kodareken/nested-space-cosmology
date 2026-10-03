"""Exact homogeneous leading-Einstein/three-AP-sector turn reduction.

This is the existing four-geometry action restricted to its symmetry-invariant
uniform covariance, not a new dynamics. The independent held-out measurement
backend is invoked only after its source-radius prediction is locked.
"""
from __future__ import annotations
import hashlib
import importlib
import io
import json
import os
from pathlib import Path
import subprocess
import time

import numpy as np
from scipy.integrate import solve_ivp

from . import nsc_discovery_dynamic_preparation as preparation
from . import nsc_discovery_leading_einstein as leading
from . import provenance

LAB=Path(__file__).resolve().parents[2]
ROOT=LAB.parent
INPUT=LAB/'results/development/nsc-discovery-dynamic-preparation-v2'
INITIAL=LAB/'results/development/nsc-discovery-leading-einstein-v1/nf128_uniform_dt0.0005-000000'
REFERENCE=LAB/'results/development/nsc-discovery-leading-einstein-v2/nf128_uniform_dt0.0005-000001'
OUTPUT=LAB/'results/development/nsc-discovery-homogeneous-turn-v1'
SCHEMA='NSC-DISCOVERY-HOMOGENEOUS-TURN-v1'
ALPHA=.002
DIRECTION=np.array([1.,0.,-1.])
BASE_C=np.array([.75,.5,.25])
CPU_LIMIT=30.
CHUNK_LIMIT=64*1024*1024
S1=np.array([[0.,1.],[1.,0.]],complex)
S2=np.array([[0.,-1j],[1j,0.]],complex)
S3=np.diag([1.,-1.])
BACKEND_MODULE='recursive_horizons.nsc_discovery_turn_fullfield'


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def producers():
    paths=[Path(__file__),LAB/'scripts/derive_nsc_discovery_homogeneous_turn.py',LAB/'tests/test_nsc_discovery_homogeneous_turn.py',
           LAB/'docs/nsc-discovery-homogeneous-turn.md',
           LAB/'src/recursive_horizons/nsc_discovery_turn_fullfield.py',LAB/'tests/test_nsc_discovery_turn_fullfield.py',
           LAB/'src/recursive_horizons/nsc_discovery_leading_step_control.py',Path(provenance.__file__),
           Path(leading.__file__),Path(preparation.__file__),Path(leading.coupling.__file__),Path(leading.galerkin.__file__),
           LAB/'src/recursive_horizons/nsc_spherical_feedback_action.py',LAB/'src/recursive_horizons/nsc_discovery_backend.py']
    return dict(leading.source_hashes(),**{str(path.relative_to(ROOT)):sha(path) for path in paths})


def _commit(pins):
    try:
        commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT).decode().strip()
        if all(hashlib.sha256(subprocess.check_output(['git','cat-file','blob',commit+':'+path],cwd=ROOT,stderr=subprocess.DEVNULL)).hexdigest()==digest for path,digest in pins.items()):return commit
    except (OSError,subprocess.SubprocessError):pass
    return None


def _check_inputs(hashes):
    if any(sha(path)!=digest for path,digest in hashes.items()):raise ValueError('frozen turn input changed')


def unpack(y):return np.asarray(y[:4],float),(np.asarray(y[4:10])+1j*np.asarray(y[10:16])).reshape(3,2)


def pack(geometry,u):return np.r_[geometry,np.asarray(u).real.ravel(),np.asarray(u).imag.ravel()]


def spin(u,matrix):return np.einsum('ji,ik,jk->j',u.conj(),matrix,u).real


def quantities(y,c,k):
    (q,r,p,v),u=unpack(y);sx=spin(u,S1);sy=spin(u,S2)
    E=2.*k['M']*np.dot(c,k['momenta']*sy+k['kappa']*q*sx)
    S=2.*k['M']*k['kappa']*np.dot(c,sx)
    a=k['a'];V=-a*r*r-k['mag'];Z=k['Z']
    Hg=k['period']*(q*p*v/(2*a*r)-Z*q*q*p*p/(4*a*a*r*r)-q*q*V)
    return {'field_energy':float(E),'geometry_energy':float(Hg),'total_energy':float(E+Hg),'source_Q_force':float(S)}


def rhs(t,y,c,k):
    (q,r,p,v),u=unpack(y);a=k['a'];ell=k['period'];kap=k['kappa'];mom=k['momenta']
    V=-a*r*r-k['mag'];S=2*k['M']*kap*np.dot(c,spin(u,S1))
    geo=[q*v/(2*a*r)-3*q*q*p/(2*a*r*r),q*p/(2*a*r),
         -p*v/(2*a*r)+3*q*p*p/(2*a*r*r)+2*q*V-S/ell,
         q*p*v/(2*a*r*r)-3*q*q*p*p/(2*a*r**3)-2*a*r*q*q]
    du=np.stack(((-mom-1j*kap*q)*u[:,1],(mom-1j*kap*q)*u[:,0]),axis=1)
    return pack(geo,du)


def jvp(y,delta,c,dc,k):
    """Analytic full state/source directional image; no finite differences."""
    (q,r,p,v),u=unpack(y);(dq,dr,dp,dv),du=unpack(delta)
    a=k['a'];ell=k['period'];kap=k['kappa'];mom=k['momenta'];V=-a*r*r-k['mag']
    dsx=2*np.einsum('ji,ik,jk->j',du.conj(),S1,u).real
    dS=2*k['M']*kap*(np.dot(dc,spin(u,S1))+np.dot(c,dsx))
    geo=[(dq*v+q*dv)/(2*a*r)-q*v*dr/(2*a*r*r)
          -3/(2*a)*((2*q*dq*p+q*q*dp)/r**2-2*q*q*p*dr/r**3),
         (dq*p+q*dp)/(2*a*r)-q*p*dr/(2*a*r*r),
         -(dp*v+p*dv)/(2*a*r)+p*v*dr/(2*a*r*r)
          +3/(2*a)*((dq*p*p+2*q*p*dp)/r**2-2*q*p*p*dr/r**3)+2*dq*V-4*a*q*r*dr-dS/ell,
         (dq*p*v+q*dp*v+q*p*dv)/(2*a*r*r)-q*p*v*dr/(a*r**3)
          -3/(2*a)*((2*q*dq*p*p+2*q*q*p*dp)/r**3-3*q*q*p*p*dr/r**4)-2*a*dr*q*q-4*a*r*q*dq]
    field=np.stack(((-mom-1j*kap*q)*du[:,1]-1j*kap*dq*u[:,1],
                    (mom-1j*kap*q)*du[:,0]-1j*kap*dq*u[:,0]),axis=1)
    return pack(geo,field)


def initial_state(c,k,Q0):
    epsilon=np.sqrt(k['momenta']**2+(k['kappa']*Q0)**2)
    u=np.stack((np.ones(3),(k['kappa']*Q0+1j*k['momenta'])/epsilon),axis=1)/np.sqrt(2.)
    E0=2*k['M']*np.dot(c,epsilon)
    square=k['mag']/(-k['a'])+E0/(-k['a']*k['period']*Q0*Q0)
    if square<=0. or np.min(c)<0. or np.max(c)>1.:raise ValueError('nonadmissible homogeneous source/constraint preparation')
    return pack([Q0,np.sqrt(square),0.,0.],u)


def initial_tangent(c,dc,k,Q0):
    y=initial_state(c,k,Q0);r=y[1];eps=np.sqrt(k['momenta']**2+(k['kappa']*Q0)**2)
    dE=2*k['M']*np.dot(dc,eps)
    dr=dE/(2*(-k['a'])*k['period']*Q0*Q0*r)
    return np.r_[0.,dr,np.zeros(14)]


def reconstruct_columns(y,frame,k,nf):
    """Recover the original standing columns from fixed AP frame coefficients."""
    _geo,u=unpack(y);x=np.arange(nf)*k['period']/nf
    result=np.zeros((nf,2,6),complex)
    for j,p in enumerate(k['momenta']):
        result+=np.exp(1j*p*x)[:,None,None]*u[j][None,:,None]*frame[2*j][None,None,:]/np.sqrt(nf)
        result+=np.exp(-1j*p*x)[:,None,None]*(S1@u[j])[None,:,None]*frame[2*j+1][None,None,:]/np.sqrt(nf)
    return result[:,0,:],result[:,1,:]


def load_preparation():
    prep_binding=preparation.authenticate_preparation(INPUT)
    record=json.loads(INPUT.with_suffix('.json').read_text())
    paths=[INPUT.with_suffix('.json'),INPUT.with_suffix('.npz'),Path(str(INITIAL)+'.json'),Path(str(INITIAL)+'.npz'),
           Path(str(REFERENCE)+'.json'),Path(str(REFERENCE)+'.npz'),
           INITIAL.parent/'manifest.json',REFERENCE.parent/'manifest.json',
           REFERENCE.parent/'nf128_uniform_dt0.0005-observations.jsonl']
    hashes={str(path):sha(path) for path in paths}
    historical_sources={}
    for directory in (INITIAL.parent,REFERENCE.parent):
        manifest=json.loads((directory/'manifest.json').read_text())
        for name,digest in manifest['source_hashes'].items():
            provenance.resolve_pinned_source_bytes(ROOT,name,digest,commit=manifest['producing_commit'])
        historical_sources[str(directory)]={'producing_commit':manifest['producing_commit'],'source_hashes':manifest['source_hashes']}
    for stem in (INITIAL,REFERENCE):
        meta=json.loads(Path(str(stem)+'.json').read_text())
        if sha(Path(str(stem)+'.npz'))!=meta['arrays_sha256']:raise ValueError('leading checkpoint payload binding failed')
    with np.load(Path(str(INITIAL)+'.npz'),allow_pickle=False) as saved:
        nf=len(saved['phi0']);W=saved['W'];dx=record['period']/len(saved['Q'])
        geo=[float(np.mean(W@saved[name]/dx if name.startswith('pi') else W@saved[name])) for name in ('Q','r','pi_Q','pi_r')]
        original0=saved['phi0'].copy();original1=saved['phi1'].copy();weights=saved['source_weights'].copy()
        geom_spread=[float(np.ptp(W@saved[name]/dx if name.startswith('pi') else W@saved[name])) for name in ('Q','r','pi_Q','pi_r')]
    coeff=record['locked_coefficients'];A=coeff['A'];a=-8*np.pi*A
    k={'A':A,'C_F':coeff['C_F'],'flux':coeff['flux'],'a':a,'Z':3*a,
       'mag':2*np.pi*coeff['C_F']*coeff['flux']**2,'period':record['period'],'kappa':record['multiplicity']/4.,
       'M':record['multiplicity'],'momenta':2*np.pi/record['period']*np.array([.5,1.5,2.5])}
    y0=initial_state(BASE_C,k,geo[0]);_geom,u=unpack(y0)
    demod=np.exp(-1j*np.pi*np.arange(nf)/nf)[:,None]
    modal=np.stack((np.fft.fft(original0*demod,axis=0),np.fft.fft(original1*demod,axis=0)),axis=1)/np.sqrt(nf)
    ids=[0,nf-1,1,nf-2,2,nf-3];frame=np.zeros((6,6),complex)
    for j in range(3):frame[2*j]=u[j].conj()@modal[ids[2*j]];frame[2*j+1]=(S1@u[j]).conj()@modal[ids[2*j+1]]
    rebuilt0,rebuilt1=reconstruct_columns(y0,frame,k,nf)
    tail=float(np.linalg.norm(modal[[index for index in range(nf) if index not in ids]]))
    cross=max(float(np.max(np.abs((modal[i]*weights)@modal[j].conj().T))) for i in ids for j in ids if i!=j)
    reconstruction=max(float(np.max(np.abs(rebuilt0-original0))),float(np.max(np.abs(rebuilt1-original1))))
    if max(geom_spread)>1e-10 or reconstruction>1e-10 or cross>1e-10 or tail>1e-10:
        raise ValueError('stored uniform source/geometry does not close on the homogeneous symmetry reduction')
    if abs(geo[1]-y0[1])>1e-10 or any(abs(value)>1e-12 for value in geo[2:]):raise ValueError('leading initial constraint/momenta do not match the homogeneous preparation')
    with np.load(Path(str(REFERENCE)+'.npz'),allow_pickle=False) as final:
        W=final['W'];dx=k['period']/len(final['Q'])
        final_geo=[float(np.mean(W@final[name]/dx if name.startswith('pi') else W@final[name])) for name in ('Q','r','pi_Q','pi_r')]
    rows=[json.loads(line) for line in (REFERENCE.parent/'nf128_uniform_dt0.0005-observations.jsonl').read_text().splitlines()]
    minimum=min(rows,key=lambda row:row['r']['min'])
    actual_final=next(row for row in rows if abs(row['time']-3.)<1e-10)
    return {'constants':k,'Q0':geo[0],'initial_state':y0,'frame_coefficients':frame,'input_hashes':hashes,'preparation_binding':prep_binding,
            'historical_sources':historical_sources,
            'closure':{'geometry_spreads':geom_spread,'inactive_modal_amplitude_tail':tail,'cross_momentum_covariance_max':cross,
                       'standing_column_reconstruction_max':reconstruction,'frame_Gram_gap':float(np.max(np.abs(frame@frame.conj().T-np.eye(6)))),
                       'exact_invariance_is_symmetry_statement':True,'actual_stored_input_agrees_up_to_recorded_roundoff':True},
            'fullband_reference':{'time':3.,'geometry':final_geo,'energy':actual_final['energy'],
                                 'sampled_minimum_time':minimum['time'],'sampled_minimum_radius':minimum['r']['min'],
                                 'sampled_minimum_field_energy':minimum['energy']['field'],'sampled_minimum_is_not_rooted_event':True}}


def integrate(c,k,Q0,*,direction=None,cpu_remaining=CPU_LIMIT):
    start=time.process_time();y0=initial_state(c,k,Q0)
    if direction is not None:y0=np.r_[y0,initial_tangent(c,direction,k,Q0)]
    y0=np.r_[y0,0.]  # work ledger is an auxiliary integral, not a dynamical variable.
    def fun(t,y):
        if time.process_time()-start>cpu_remaining:raise RuntimeError('homogeneous turn CPU budget exhausted')
        value=rhs(t,y[:16],c,k)
        work=quantities(y[:16],c,k)['source_Q_force']*value[0]
        return np.r_[value,work] if direction is None else np.r_[value,jvp(y[:16],y[16:32],c,direction,k),work]
    def turn(t,y):return y[2]
    turn.direction=-1
    def zero(t,y):return quantities(y[:16],c,k)['field_energy']
    zero.direction=-1
    solved=solve_ivp(fun,(0.,3.),y0,method='DOP853',rtol=2e-12,atol=2e-14,events=(turn,zero),dense_output=False)
    if not solved.success or len(solved.t_events[0])!=1 or not len(solved.t_events[1]):raise RuntimeError('required homogeneous event did not occur in T<=3')
    event=solved.y_events[0][0];state=event[:16];t=solved.t_events[0][0];rate=rhs(t,state,c,k)
    q,r,p,v=state[:4];values=quantities(state,c,k);S=values['source_Q_force'];E=values['field_energy']
    acceleration=q*rate[2]/(2*k['a']*r)
    predicted_square=k['mag']/(-k['a'])+E/(-k['a']*k['period']*q*q)
    _g,u=unpack(state);eps=np.sqrt(k['momenta']**2+(k['kappa']*q)**2)
    minus=(1-(k['momenta']*spin(u,S2)+k['kappa']*q*spin(u,S1))/eps)/2
    report={'time':float(t),'radius':float(r),'Q':float(q),'p_Q':float(p),'p_Q_dot':float(rate[2]),
            'event_time_scope':'coordinate T; equal-proper-time comparison is not asserted',
            'r_ddot':float(acceleration),'field_energy_zero_time':float(solved.t_events[1][0]),
            'event_energy':values,'turn_radius_squared_from_constraint':float(predicted_square),
            'turn_constraint_radius_squared_gap':float(r*r-predicted_square),'negative_sector_probabilities':minus.tolist(),
            'spinor_norm_gap':float(np.max(np.abs(np.sum(abs(u)**2,axis=1)-1.))),
            'final_geometry':solved.y[:4,-1].tolist(),'final_energy':quantities(solved.y[:16,-1],c,k),
            'CPU_seconds':time.process_time()-start,'nfev':solved.nfev}
    initial_energy=quantities(y0[:16],c,k)
    final_energy=report['final_energy']
    report.update(initial_energy=initial_energy,integrated_fieldwork=float(solved.y[-1,-1]),
         fieldwork_closure_gap=float(final_energy['field_energy']-initial_energy['field_energy']-solved.y[-1,-1]),
         geometry_transfer_closure_gap=float(final_energy['geometry_energy']-initial_energy['geometry_energy']+solved.y[-1,-1]),
         maximum_total_energy_residual=max(abs(quantities(row,c,k)['total_energy']) for row in solved.y[:16].T),
         source_CAR_eigenvalues=np.repeat(c,2).tolist(),source_trace=float(2*np.sum(c)),
         Q_positive=bool(np.min(solved.y[0])>0),r_positive=bool(np.min(solved.y[1])>0),
         turn_acceleration_from_energy=float((2*E-q*S)/(2*k['a']*r*k['period'])),
         submagnetic_turn_requires_negative_field_energy=bool(r*r<k['mag']/(-k['a']) and E<0))
    if direction is not None:
        delta=event[16:32];report.update(radius_derivative=float(delta[1]),event_time_derivative=float(-delta[2]/rate[2]),
            event_time_correction_in_radius=0.,radius_formula='rstar(alpha)=rstar(0)+alpha delta_r(tstar), since rdot(tstar)=0')
    return report,{'event_state':state,'event_tangent':event[16:32] if direction is not None else np.zeros(16),'final_state':solved.y[:16,-1]}


def _json(value):
    if isinstance(value,np.ndarray):return value.tolist()
    if isinstance(value,np.generic):return value.item()
    if isinstance(value,dict):return {key:_json(item) for key,item in value.items()}
    if isinstance(value,(list,tuple)):return [_json(item) for item in value]
    return value


def _write(directory,prefix,record,arrays):
    directory=Path(directory);jp=directory/(prefix+'.json');npz=directory/(prefix+'.npz')
    if jp.exists() or npz.exists():raise FileExistsError('immutable homogeneous-turn prefix already exists')
    stream=io.BytesIO();np.savez_compressed(stream,**arrays);payload=stream.getvalue()
    record=dict(record,payload_sha256=hashlib.sha256(payload).hexdigest(),array_sha256={key:hashlib.sha256(np.ascontiguousarray(value).tobytes()).hexdigest() for key,value in arrays.items()})
    text=(json.dumps(_json(record),sort_keys=True,indent=2,allow_nan=False)+'\n').encode()
    if len(payload)+len(text)>CHUNK_LIMIT:raise RuntimeError('turn record exceeds 64 MiB')
    directory.mkdir(parents=True,exist_ok=True)
    for path,blob in ((npz,payload),(jp,text)):
        descriptor=os.open(path,os.O_CREAT|os.O_EXCL|os.O_WRONLY,0o444)
        with os.fdopen(descriptor,'wb') as output:output.write(blob)
    return record


def _read(directory,prefix):
    directory=Path(directory);record=json.loads((directory/(prefix+'.json')).read_text());path=directory/(prefix+'.npz')
    if sha(path)!=record['payload_sha256']:raise ValueError('turn payload binding failed')
    with np.load(path,allow_pickle=False) as saved:arrays={key:saved[key].copy() for key in saved.files}
    for key,digest in record['array_sha256'].items():
        if hashlib.sha256(np.ascontiguousarray(arrays[key]).tobytes()).hexdigest()!=digest:raise ValueError('turn array binding failed')
    _check_inputs(record['input_hashes'])
    commit=record.get('producing_commit')
    if record['producers']!=producers() and not commit:raise ValueError('historical uncommitted turn producers do not authenticate')
    if commit:
        for path,digest in record['producers'].items():
            provenance.resolve_pinned_source_bytes(ROOT,path,digest,commit=commit)
    return record,arrays


def preview():
    data=load_preparation();k=data['constants']
    return {'schema':SCHEMA,'mode':'preview','evolved':False,'constants':_json(k),'Q0':data['Q0'],'closure':data['closure'],
            'source_pair_occupations':BASE_C.tolist(),'held_pair_occupations':np.repeat(BASE_C+ALPHA*DIRECTION,2).tolist(),
            'alpha':ALPHA,'direction_per_pair':DIRECTION.tolist(),'reference':data['fullband_reference'],
            'selected_scale_requires_actual_state_history':True,'physical_vacuum_Gamma_matching':'unresolved','EFT_strong_curvature_domain':'unresolved'}


def prepare(directory):
    if (Path(directory)/'prepare.json').exists():raise FileExistsError('turn preparation already exists')
    start=time.process_time();pins=producers();data=load_preparation();k=data['constants']
    initial_spinors=unpack(data['initial_state'])[1]
    description={'input_record':str(INPUT.with_suffix('.json')),'input_payload':str(INPUT.with_suffix('.npz')),
       'leading_initial_record':str(Path(str(INITIAL)+'.json')),'leading_initial_payload':str(Path(str(INITIAL)+'.npz')),
       'source_case':'uniform','Q0':data['Q0'],'period':k['period'],'kappa':k['kappa'],'multiplicity':k['M'],
       'baseline_pair_occupations':np.repeat(BASE_C,2).tolist(),'held_pair_occupations':np.repeat(BASE_C+ALPHA*DIRECTION,2).tolist(),
       'direction_per_pair':DIRECTION.tolist(),'initial_momenta_zero':True,'constants':_json(k),'input_hashes':data['input_hashes'],
       'frame_mode_order':[.5,-.5,1.5,-1.5,2.5,-2.5],'prepared_payload':str(Path(directory).resolve()/'prepare.npz'),
       'frame_coefficients_key':'frame_coefficients','initial_spinors_key':'initial_spinors','momenta_key':'momenta'}
    return _write(directory,'prepare',{'schema':SCHEMA,'mode':'prepared_invariant_source','producers':pins,'producing_commit':_commit(pins),
        'input_hashes':data['input_hashes'],'preparation':description,'closure':data['closure'],'fullband_reference':data['fullband_reference'],
        'preparation_binding':data['preparation_binding'],'historical_sources':data['historical_sources'],'CPU_seconds':time.process_time()-start,'held_data_generated':False},
        {'initial_state':data['initial_state'],'initial_spinors':initial_spinors,'frame_coefficients':data['frame_coefficients'],'momenta':k['momenta']})


def predict(directory):
    if (Path(directory)/'prediction.json').exists():raise FileExistsError('turn prediction already exists')
    start=time.process_time();prepared,arrays=_read(directory,'prepare');pins=producers()
    if pins!=prepared['producers']:raise ValueError('prediction needs current frozen producers')
    k=dict(prepared['preparation']['constants'],momenta=arrays['momenta']);Q0=prepared['preparation']['Q0']
    baseline,payload=integrate(BASE_C,k,Q0,direction=DIRECTION,cpu_remaining=CPU_LIMIT-prepared['CPU_seconds'])
    centred=[]
    for h in (.0001,.00005):
        plus,_=integrate(BASE_C+h*DIRECTION,k,Q0,cpu_remaining=CPU_LIMIT-prepared['CPU_seconds']-(time.process_time()-start))
        minus,_=integrate(BASE_C-h*DIRECTION,k,Q0,cpu_remaining=CPU_LIMIT-prepared['CPU_seconds']-(time.process_time()-start))
        centred.append({'h':h,'event_radius_derivative':(plus['radius']-minus['radius'])/(2*h),
                        'analytic_difference':(plus['radius']-minus['radius'])/(2*h)-baseline['radius_derivative'],
                        'event_time_derivative':(plus['time']-minus['time'])/(2*h),'held_alpha_not_used':True})
    reference=prepared['fullband_reference'];reference_gaps=np.asarray(baseline['final_geometry'])-np.asarray(reference['geometry'])
    field_gap=baseline['final_energy']['field_energy']-reference['energy']['field']
    record={key:value for key,value in prepared.items() if key not in ('payload_sha256','array_sha256')}
    record.update(mode='locked_prediction',alpha=ALPHA,locked_before_held_data=True,held_data_generated=False,
         baseline=baseline,centred_event_checks=centred,
         predicted_held_radius=baseline['radius']+ALPHA*baseline['radius_derivative'],
         predicted_radius_change=ALPHA*baseline['radius_derivative'],
         predicted_event_time=baseline['time']+ALPHA*baseline['event_time_derivative'],
         fullband_T3_comparison={'geometry_gaps':reference_gaps.tolist(),'field_energy_gap':float(field_gap),'indicator_not_error_bound':True},
         physical_domain='same e6 leading action; full finite CAR source; magnetic term unchanged',
         physical_vacuum_Gamma_matching='unresolved',strong_curvature_EFT_validity='unresolved',
         event_is_not_bounce_novelty_or_holding=True,scale_is_selected_by_actual_source_history=True,
         prepare_json_sha256=sha(Path(directory)/'prepare.json'),CPU_seconds=prepared['CPU_seconds']+time.process_time()-start)
    _check_inputs(prepared['input_hashes'])
    if pins!=producers():raise ValueError('producer source changed during turn prediction')
    return _write(directory,'prediction',record,payload)


def measure(directory):
    if (Path(directory)/'measurement.json').exists():raise FileExistsError('turn measurement already exists')
    start=time.process_time();locked,_arrays=_read(directory,'prediction');pins=producers()
    if not locked['locked_before_held_data'] or locked['mode']!='locked_prediction':raise ValueError('held measurement requires a locked prediction')
    if locked['producers']!=pins:raise ValueError('held measurement requires current frozen producers')
    backend=importlib.import_module(BACKEND_MODULE)
    records,arrays=backend.measure_locked_prediction(Path(directory)/'prediction.json',cpu_limit_seconds=min(20.,CPU_LIMIT-locked['CPU_seconds']))
    comparisons=[]
    for event in records['events']:
        if event['source']!='held' or event.get('status') not in (None,'FIRST_POSITIVE_TO_NEGATIVE_P_Q_EVENT'):continue
        value=event['radius'];change=value-locked['baseline']['radius'];error=value-locked['predicted_held_radius']
        comparisons.append({'step_cap':event.get('step_cap'),'radius':value,'time':event['time'],'measured_change':change,'predicted_change':locked['predicted_radius_change'],
                            'prediction_error':error,'error_over_effect':None if change==0 else abs(error/change),
                            'event_time_prediction_error':event['time']-locked['predicted_event_time']})
    _check_inputs(locked['input_hashes'])
    if pins!=producers():raise ValueError('source changed during independent turn measurement')
    record={'schema':SCHEMA,'mode':'independent_fullfield_measurement','producers':pins,'producing_commit':locked['producing_commit'],
            'input_hashes':locked['input_hashes'],'prediction_json_sha256':sha(Path(directory)/'prediction.json'),
            'backend':records,'comparisons':comparisons,'complete_matched_events':len(comparisons)==2,
            'finite_checkpoint_retained_if_deferred':len(comparisons)!=2,'aggregate_CPU_seconds':locked['CPU_seconds']+time.process_time()-start,
            'CPU_limit_seconds':CPU_LIMIT,'held_measurement_not_fit':True,'physical_vacuum_Gamma_matching':'unresolved',
            'strong_curvature_EFT_validity':'unresolved','holding_or_bounce_novelty_claimed':False}
    if record['aggregate_CPU_seconds']>CPU_LIMIT:raise RuntimeError('independent measurement exceeded aggregate 30 CPU budget')
    return _write(directory,'measurement',record,arrays)


def check(directory):
    prepared,_source=_read(directory,'prepare')
    if not (Path(directory)/'prediction.json').exists():return {'schema':SCHEMA,'mode':'read_only_check','ok':True,'evolved':False}
    locked,arrays=_read(directory,'prediction')
    if locked['prepare_json_sha256']!=sha(Path(directory)/'prepare.json'):raise ValueError('turn prediction preparation binding failed')
    k=dict(locked['preparation']['constants'],momenta=np.asarray(locked['preparation']['constants']['momenta']))
    y=arrays['event_state'];delta=arrays['event_tangent']
    if abs(y[1]+ALPHA*delta[1]-locked['predicted_held_radius'])>1e-12:raise ValueError('turn radius prediction replay failed')
    rate=rhs(locked['baseline']['time'],y,BASE_C,k)
    if abs(-delta[2]/rate[2]-locked['baseline']['event_time_derivative'])>1e-10:raise ValueError('turn event-time tangent replay failed')
    report={'schema':SCHEMA,'mode':'read_only_check','ok':True,'evolved':False,'prediction_replayed':True,'held_data_used_in_prediction':False}
    if (Path(directory)/'measurement.json').exists():
        measured,_endpoints=_read(directory,'measurement')
        if measured['prediction_json_sha256']!=sha(Path(directory)/'prediction.json'):raise ValueError('independent measurement prediction binding failed')
        report.update(measurement_authenticated=True,comparisons=measured['comparisons'])
    return report
