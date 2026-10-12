"""Finite static leading-action collar: source response and required traction.

A local real standing Dirac solution is prescribed at the symmetric center.
Its amplitude weight is not assumed to be a globally normalized occupation.
"""
from __future__ import annotations
import hashlib
import io
import json
import os
from pathlib import Path
import subprocess
import time

import numpy as np
from scipy.integrate import solve_ivp
from scipy.integrate import simpson

from . import nsc_discovery_dynamic_preparation as preparation
from . import nsc_discovery_leading_einstein as leading
from . import provenance

LAB=Path(__file__).resolve().parents[2]
ROOT=LAB.parent
INPUT=LAB/'results/development/nsc-discovery-dynamic-preparation-v2'
OUTPUT=LAB/'results/development/nsc-discovery-parent-traction-v1'
SCHEMA='NSC-DISCOVERY-PARENT-TRACTION-v1'
HELD_WEIGHT=.0013
CPU_LIMIT=30.
CHUNK_LIMIT=64*1024**2
CUTS={'child':.5,'parent':1.}


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def source_hashes():
    paths=[Path(__file__),LAB/'scripts/derive_nsc_discovery_parent_traction.py',
           LAB/'tests/test_nsc_discovery_parent_traction.py',LAB/'docs/nsc-discovery-parent-traction.md',
           Path(leading.__file__),Path(preparation.__file__),Path(provenance.__file__)]
    return {str(path.relative_to(ROOT)):sha(path) for path in paths}


def producing_commit(pins):
    try:
        commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT).decode().strip()
        for path,digest in pins.items():provenance.resolve_pinned_source_bytes(ROOT,path,digest,commit=commit)
        return commit
    except (OSError,subprocess.SubprocessError,RuntimeError):return None


def load_inputs():
    binding=preparation.authenticate_preparation(INPUT)
    record=json.loads(INPUT.with_suffix('.json').read_text());locked=record['locked_coefficients']
    g=8*np.pi*locked['A'];mag=2*np.pi*locked['C_F']*locked['flux']**2
    if not np.isclose(mag/g,1.,rtol=0.,atol=5e-14):raise ValueError('reviewed finite BR collar requires locked mag/g=1')
    k={'A':locked['A'],'C_F':locked['C_F'],'flux':locked['flux'],'g':g,'mag':mag,'M':4.,'kappa':1.,'epsilon':1.5}
    if record['multiplicity']!=4:raise ValueError('unchanged multiplicity four required')
    return k,{str(INPUT.with_suffix(extension)):sha(INPUT.with_suffix(extension)) for extension in ('.json','.npz')},binding


def source_values(y,k):
    r,s,h,hx,u,v=y[:6];Q=np.exp(h);n=u*u+v*v;S=2*u*v
    return Q,n,S,k['epsilon']*n-k['kappa']*Q*S


def rhs(x,y,w,k):
    r,s,h,hx,u,v=y[:6];Q,n,S,K=source_values(y,k)
    rxx=Q*Q*(r*r-k['mag']/k['g'])/r-s*s/r-k['M']*w*k['kappa']*Q*S/(2*k['g']*r)
    hxx=Q*Q-3*rxx/r
    return np.array([s,rxx,hx,hxx,k['epsilon']*v-k['kappa']*Q*u,
                     k['kappa']*Q*v-k['epsilon']*u,n,Q*Q,r*Q])


def first_integral(y,w,k):
    r,s,h,hx=y[:4];Q,n,S,K=source_values(y,k)
    return Q*Q*(r*r-k['mag']/k['g'])-3*s*s-2*r*s*hx+k['M']*w*K/k['g']


def first_integral_gradient(y,w,k):
    r,s,h,hx,u,v=y[:6];Q,n,S,K=source_values(y,k);mw=k['M']*w/k['g'];kap=k['kappa'];eps=k['epsilon']
    return np.array([2*Q*Q*r-2*s*hx,-6*s-2*r*hx,2*Q*Q*(r*r-k['mag']/k['g'])-mw*kap*Q*S,
                     -2*r*s,mw*(2*eps*u-2*kap*Q*v),mw*(2*eps*v-2*kap*Q*u),0.,0.,0.])


def center_state(w,k):
    if not np.isfinite(w) or w<0:raise ValueError('source amplitude weight must be finite and nonnegative')
    square=k['mag']/k['g']-k['M']*w/(k['g']*np.pi)
    if square<=0:raise ValueError('source center has no positive areal radius')
    return np.array([np.sqrt(square),0.,0.,0.,1/np.sqrt(np.pi),1/np.sqrt(np.pi),0.,0.,0.])


def radius_response(x,k):
    x=np.asarray(x);c=k['M']/(4*np.pi*k['g']);sec2=1/np.cos(x)**2;tan=np.tan(x)
    a=-c*(2+np.sin(x)**2+3*x*tan)
    ax=-c*(np.sin(2*x)+3*tan+3*x*sec2)
    axx=-c*(2*np.cos(2*x)+6*sec2+6*x*sec2*tan)
    return a,ax,axx


def baseline_spinor(x):
    x=np.asarray(x);cos=np.cos(x);sin=np.sin(x)
    return cos*np.sqrt(1+sin)/np.sqrt(np.pi),cos*np.sqrt(1-sin)/np.sqrt(np.pi)


def _grid(points):
    if int(points)!=points or points<3:raise ValueError('collar grid needs at least three integer points')
    x=np.linspace(0.,1.,int(points))
    if any(np.min(abs(x-cut))>1e-12 for cut in CUTS.values()):raise ValueError('collar grid must contain every declared cut')
    return x


def linear_solution(k,*,points=257,cpu_remaining=CPU_LIMIT):
    """Predictor solves only the reviewed first-order q equation, never held IVP."""
    started=time.process_time();x=_grid(points)
    def equation(x,y):
        if time.process_time()-started>cpu_remaining:raise RuntimeError('linear collar CPU budget exhausted')
        a,ax,axx=radius_response(x,k);sec=1/np.cos(x)
        return [y[1],2*sec*sec*y[0]-3*axx,sec,sec*(a+y[0]),sec*sec*y[0]]
    result=solve_ivp(equation,(0.,1.),np.zeros(5),method='DOP853',t_eval=x,rtol=2e-12,atol=2e-14)
    if not result.success:raise RuntimeError('linear collar integration failed')
    q,qx,length0,dlength,magnetic_response=result.y;a,ax,axx=radius_response(x,k);u,v=baseline_spinor(x)
    return {'x':x,'a':a,'ax':ax,'axx':axx,'q':q,'qx':qx,'baseline_u':u,'baseline_v':v,
            'length0_half':length0,'dlength_half':dlength,'magnetic_response_half':magnetic_response}


def linear_readouts(arrays,k,w):
    x=arrays['x'];out={};a0=float(arrays['a'][0])
    for name,cut in CUTS.items():
        index=int(np.argmin(abs(x-cut)))
        if abs(x[index]-cut)>1e-12:raise ValueError('linear grid must contain both declared cuts')
        a=arrays['a'][index];ax=arrays['ax'][index];q=arrays['q'][index];qx=arrays['qx'][index]
        sec=1/np.cos(cut);tan=np.tan(cut);dB=k['g']*(ax+2*a*tan+qx)
        norm0=(2*cut+np.sin(2*cut))/np.pi;energy=k['M']*k['epsilon']*w*norm0
        magnetic=2*k['mag']*tan+4*k['mag']*w*arrays['magnetic_response_half'][index]
        out[name]={'cut':cut,'radius':float(1+w*a),'radius_change':float(w*a),
             'Q':float(sec*(1+w*q)),'Q_change':float(w*sec*q),
             'proper_lapse_N':float(sec*(1+w*(a+q))),'proper_lapse_N_change':float(w*sec*(a+q)),
             'center_clock_lapse_ratio':float(sec*(1+w*(a+q-a0))),
             'boundary_gradient_r_x':float(w*ax),'boundary_gradient_h_x':float(tan+w*qx),
             'boundary_gradient_Q_x':float(sec*(tan+w*(qx+q*tan))),
             'right_endpoint_charge':float(k['g']*tan+w*dB),'right_endpoint_charge_change':float(w*dB),
             'full_collar_boundary_charge':float(2*(k['g']*tan+w*dB)),
             'full_collar_boundary_charge_change':float(2*w*dB),
             'finite_collar_covariance_eigenvalue':float(w*norm0),'partial_source_energy':float(energy),
             'magnetic_integral':float(magnetic),'linear_static_boundary_identity_residual':float(2*(k['g']*tan+w*dB)-energy-magnetic),
             'proper_radial_length':float(2*(arrays['length0_half'][index]+w*arrays['dlength_half'][index])),
             'proper_radial_length_change':float(2*w*arrays['dlength_half'][index])}
    return out


def exact_local(w,k,*,rtol=2e-11,atol=2e-13,points=257,cpu_remaining=CPU_LIMIT):
    """Independent nonlinear local IVP; no linear predictor or predicted endpoint."""
    started=time.process_time();x=_grid(points);y0=center_state(w,k)
    def equation(x,y):
        if time.process_time()-started>cpu_remaining:raise RuntimeError('nonlinear collar CPU budget exhausted')
        if not np.isfinite(y).all() or y[0]<=0:raise ValueError('local positive-radius chart exited')
        return rhs(x,y,w,k)
    solved=solve_ivp(equation,(0.,1.),y0,method='DOP853',t_eval=x,rtol=rtol,atol=atol)
    if not solved.success:raise RuntimeError('local collar IVP failed')
    state=solved.y.T;Q=np.exp(state[:,2]);n=state[:,4]**2+state[:,5]**2
    integrals=np.array([first_integral(row,w,k) for row in state]);readouts={}
    for name,cut in CUTS.items():
        index=int(np.argmin(abs(x-cut)));r,s,h,hx,u,v,jn,jQ,jlength=state[index];N=r*np.exp(h)
        norm=2*jn;energy=k['M']*w*k['epsilon']*norm;magnetic=2*k['mag']*jQ
        endpoint=k['g']*(r*s+r*r*hx);charge=2*endpoint
        quadrature_norm=2*simpson(n[:index+1],x=x[:index+1]);quadrature_mag=2*k['mag']*simpson(Q[:index+1]**2,x=x[:index+1])
        readouts[name]={'cut':cut,'radius':float(r),'Q':float(np.exp(h)),'proper_lapse_N':float(N),
            'center_clock_lapse_ratio':float(N/y0[0]),'boundary_gradient_r_x':float(s),'boundary_gradient_h_x':float(hx),
            'boundary_gradient_Q_x':float(np.exp(h)*hx),
            'right_endpoint_charge':float(endpoint),'left_endpoint_charge':float(-endpoint),'full_collar_boundary_charge':float(charge),
            'proper_radial_length':float(2*jlength),'finite_collar_mode_norm':float(norm),
            'finite_collar_covariance_eigenvalue':float(w*norm),'partial_source_energy':float(energy),'magnetic_integral':float(magnetic),
            'static_boundary_identity_residual':float(charge-energy-magnetic),
            'mode_norm_quadrature_gap':float(quadrature_norm-norm),'magnetic_quadrature_gap':float(quadrature_mag-magnetic),
            'covariance_admissible_on_this_collar':bool(0<=w*norm<=1)}
    report={'w':w,'readouts':readouts,'first_integral_center':float(first_integral(y0,w,k)),
            'first_integral_max':float(np.max(abs(integrals))),
            'r_min':float(np.min(state[:,0])),'Q_min':float(np.min(Q)),
            'center_radius':float(y0[0]),'center_proper_lapse_N':float(y0[0]),
            'real_standing_current_zero':True,'source_amplitude_globally_normalized':False,
            'epsilon_scope':'fixed center-coordinate frequency; center-proper frequency epsilon/N(0)',
            'center_proper_frequency':float(k['epsilon']/y0[0]),'CPU_seconds':time.process_time()-started,
            'solver':{'method':'DOP853','rtol':rtol,'atol':atol,'nfev':solved.nfev,'points':points}}
    return report,{'x':x,'state':state,'first_integral':integrals}


def _check_inputs(hashes):
    for path,digest in hashes.items():
        if sha(path)!=digest:raise ValueError('frozen collar input changed')


def _write(directory,name,record,arrays):
    directory=Path(directory);jp=directory/(name+'.json');ap=directory/(name+'.npz')
    if jp.exists() or ap.exists():raise FileExistsError('immutable collar stage already exists')
    commit=record.get('producing_commit')
    if not commit:raise ValueError('record creation requires a frozen non-None producing commit')
    if record['producers']!=source_hashes():raise ValueError('record creation requires current frozen producer bytes')
    for path,digest in record['producers'].items():provenance.resolve_pinned_source_bytes(ROOT,path,digest,commit=commit)
    stream=io.BytesIO();np.savez_compressed(stream,**arrays);payload=stream.getvalue()
    record=dict(record,payload_sha256=hashlib.sha256(payload).hexdigest(),
       array_sha256={key:hashlib.sha256(np.ascontiguousarray(value).tobytes()).hexdigest() for key,value in arrays.items()})
    text=(json.dumps(record,sort_keys=True,indent=2,allow_nan=False)+'\n').encode()
    if len(payload)+len(text)>CHUNK_LIMIT:raise RuntimeError('collar stage exceeds 64 MiB')
    directory.mkdir(parents=True,exist_ok=True)
    for path,blob in ((ap,payload),(jp,text)):
        descriptor=os.open(path,os.O_CREAT|os.O_EXCL|os.O_WRONLY,0o444)
        with os.fdopen(descriptor,'wb') as handle:handle.write(blob)
    return record


def _read(directory,name):
    directory=Path(directory);record=json.loads((directory/(name+'.json')).read_text());ap=directory/(name+'.npz')
    if record['schema']!=SCHEMA or sha(ap)!=record['payload_sha256']:raise ValueError('collar schema/payload binding failed')
    with np.load(ap,allow_pickle=False) as saved:arrays={key:saved[key].copy() for key in saved.files}
    for key,digest in record['array_sha256'].items():
        if hashlib.sha256(np.ascontiguousarray(arrays[key]).tobytes()).hexdigest()!=digest:raise ValueError('collar array binding failed')
    _check_inputs(record['input_hashes']);pins=record['producers'];commit=record.get('producing_commit')
    if pins!=source_hashes() and not commit:raise ValueError('historical collar source requires original producing commit')
    if commit:
        for path,digest in pins.items():provenance.resolve_pinned_source_bytes(ROOT,path,digest,commit=commit)
    return record,arrays


def _scope():
    return {'physical_domain':'unchanged leading Einstein/Dirac local static finite BR collar',
            'rank_one_source_normalization':'w is local amplitude weight; eigenvalue on collar = w integral n dx',
            'required_exterior_data':['endpoint r,Q and their outward gradients','Dirac continuation with compatible phase/frequency',
                                      'global mode normalization and CAR state','matching exterior boundary charge without an imposed wall'],
            'global_parent_completed':False,'bound_spectrum_completed':False,'autonomous_renewal_claimed':False,
            'prior_w_0p002_scratch_scope':'development calculation observed before lock; not held-out evidence',
            'physical_vacuum_Gamma_matching':'unresolved','strong_curvature_EFT_validity':'unresolved'}


def preview():
    k,hashes,binding=load_inputs()
    return {'schema':SCHEMA,'mode':'preview','evolved':False,'constants':k,'cuts':CUTS,'held_w':HELD_WEIGHT,
            'center_data':center_state(0.,k)[:6].tolist(),'center_constraint':'r0²=mag/g−M w/(g pi)',**_scope()}


def prepare(directory):
    start=time.process_time();k,hashes,binding=load_inputs();pins=source_hashes()
    return _write(directory,'prepare',{'schema':SCHEMA,'mode':'prepared_local_collar','producers':pins,'producing_commit':producing_commit(pins),
        'input_hashes':hashes,'input_binding':binding,'constants':k,'cuts':CUTS,'held_w':HELD_WEIGHT,
        'preparation':{'center_Q':1.,'center_r_x':0.,'center_h_x':0.,'center_u':1/np.sqrt(np.pi),'center_v':1/np.sqrt(np.pi),
                       'epsilon':1.5,'kappa':1.,'M':4.,'center_constraint_applied_for_each_w':True},
        'CPU_seconds':time.process_time()-start,'nonlinear_held_generated':False,**_scope()},
        {'baseline_center_state':center_state(0.,k),'held_center_state':center_state(HELD_WEIGHT,k)})


def predict(directory):
    start=time.process_time();prepared,_initial=_read(directory,'prepare');pins=source_hashes()
    if pins!=prepared['producers']:raise ValueError('prediction requires current frozen producers')
    arrays=linear_solution(prepared['constants'],cpu_remaining=CPU_LIMIT-prepared['CPU_seconds'])
    record={key:value for key,value in prepared.items() if key not in ('array_sha256','payload_sha256')}
    record.update(mode='locked_linear_prediction',locked_before_held_nonlinear=True,
       prediction=linear_readouts(arrays,prepared['constants'],prepared['held_w']),
       baseline=linear_readouts(arrays,prepared['constants'],0.),
       prepare_json_sha256=sha(Path(directory)/'prepare.json'),CPU_seconds=prepared['CPU_seconds']+time.process_time()-start,
       response_method='analytic a and linear q IVP on exact BR; no held nonlinear data used')
    _check_inputs(prepared['input_hashes'])
    if pins!=source_hashes():raise ValueError('producer changed during collar prediction')
    return _write(directory,'prediction',record,arrays)


def comparisons(locked,measured):
    out={}
    for name in CUTS:
        out[name]={}
        for key in ('radius','Q','proper_lapse_N','center_clock_lapse_ratio','boundary_gradient_r_x','boundary_gradient_h_x','boundary_gradient_Q_x',
                    'right_endpoint_charge','full_collar_boundary_charge','proper_radial_length',
                    'finite_collar_covariance_eigenvalue','partial_source_energy','magnetic_integral'):
            baseline=locked['baseline'][name][key];forecast=locked['prediction'][name][key];actual=measured['readouts'][name][key]
            out[name][key]={'measured':actual,'predicted':forecast,'measured_change':actual-baseline,
                          'predicted_change':forecast-baseline,'prediction_residual':actual-forecast,
                          'residual_over_effect':None if actual==baseline else abs((actual-forecast)/(actual-baseline))}
    return out


def measure(directory):
    start=time.process_time();prediction_path=Path(directory)/'prediction.json';prediction_hash=sha(prediction_path)
    locked,_linear=_read(directory,'prediction');pins=source_hashes()
    if sha(prediction_path)!=prediction_hash:raise ValueError('prediction changed during authentication')
    _read(directory,'prepare')
    if locked['prepare_json_sha256']!=sha(Path(directory)/'prepare.json'):raise ValueError('collar preparation binding failed before measurement')
    if locked.get('locked_before_held_nonlinear') is not True or locked['mode']!='locked_linear_prediction':raise ValueError('measurement needs a locked linear prediction')
    if pins!=locked['producers']:raise ValueError('measurement requires current frozen producers')
    remaining=lambda:CPU_LIMIT-locked['CPU_seconds']-(time.process_time()-start)
    primary,arrays=exact_local(locked['held_w'],locked['constants'],cpu_remaining=remaining())
    tighter,fine=exact_local(locked['held_w'],locked['constants'],rtol=2e-13,atol=2e-15,points=513,cpu_remaining=remaining())
    indicators={name:{key:tighter['readouts'][name][key]-primary['readouts'][name][key]
        for key in ('radius','Q','full_collar_boundary_charge','partial_source_energy','proper_radial_length')} for name in CUTS}
    aggregate=locked['CPU_seconds']+time.process_time()-start
    if aggregate>CPU_LIMIT:raise RuntimeError('aggregate collar budget exceeded')
    _check_inputs(locked['input_hashes'])
    if pins!=source_hashes():raise ValueError('producer changed during collar measurement')
    if sha(prediction_path)!=prediction_hash:raise ValueError('prediction changed during collar measurement')
    arrays.update({'fine_'+key:value for key,value in fine.items()})
    record={'schema':SCHEMA,'mode':'held_nonlinear_measurement','producers':pins,'producing_commit':locked['producing_commit'],
        'input_hashes':locked['input_hashes'],'prediction_json_sha256':prediction_hash,
        'measurement':primary,'tighter_measurement':tighter,'comparisons':comparisons(locked,primary),
        'tolerance_and_grid_indicators':indicators,'indicator_not_continuous_error_bound':True,
        'aggregate_CPU_seconds':aggregate,'CPU_limit_seconds':CPU_LIMIT,'prediction_not_fit':True,**_scope()}
    return _write(directory,'measurement',record,arrays)


def _verify_nonlinear(result,arrays,k,w,*,prefix=''):
    x=arrays[prefix+'x'];states=arrays[prefix+'state']
    if any(np.min(abs(x-cut))>1e-12 for cut in CUTS.values()):raise ValueError('saved collar grid misses a declared cut')
    residual=np.array([first_integral(row,w,k) for row in states])
    if not np.array_equal(residual,arrays[prefix+'first_integral']):raise ValueError('collar first-integral replay differs')
    top={'center_radius':states[0,0],'center_proper_lapse_N':states[0,0],
         'center_proper_frequency':k['epsilon']/states[0,0],'first_integral_center':residual[0],
         'first_integral_max':np.max(abs(residual)),'r_min':np.min(states[:,0]),'Q_min':np.min(np.exp(states[:,2]))}
    if any(abs(value-result[key])>1e-13 for key,value in top.items()):raise ValueError('collar center/clock/conservation replay differs')
    for name,cut in CUTS.items():
        row=states[np.argmin(abs(x-cut))];r,s,h,hx=row[:4];Q=np.exp(h);B=k['g']*(r*s+r*r*hx)
        norm=2*row[6];eigenvalue=w*norm;energy=k['M']*k['epsilon']*eigenvalue;magnetic=2*k['mag']*row[7]
        replay={'radius':r,'Q':Q,'proper_lapse_N':r*Q,'center_clock_lapse_ratio':r*Q/states[0,0],
                'boundary_gradient_r_x':s,'boundary_gradient_h_x':hx,'boundary_gradient_Q_x':Q*hx,
                'right_endpoint_charge':B,'left_endpoint_charge':-B,'full_collar_boundary_charge':2*B,
                'finite_collar_mode_norm':norm,'finite_collar_covariance_eigenvalue':eigenvalue,'partial_source_energy':energy,
                'magnetic_integral':magnetic,'proper_radial_length':2*row[8],'static_boundary_identity_residual':2*B-energy-magnetic}
        stored=result['readouts'][name]
        if any(abs(value-stored[key])>1e-13 for key,value in replay.items()):raise ValueError('collar endpoint/source/clock replay differs')
        if stored['covariance_admissible_on_this_collar']!=bool(0<=eigenvalue<=1):raise ValueError('collar covariance admissibility replay differs')


def check(directory):
    prepared,_initial=_read(directory,'prepare');report={'schema':SCHEMA,'mode':'read_only_check','ok':True,'evolved':False,
                                                       'numerical_replay':False}
    if not (Path(directory)/'prediction.json').exists():return report
    locked,linear=_read(directory,'prediction')
    if locked['prepare_json_sha256']!=sha(Path(directory)/'prepare.json'):raise ValueError('collar preparation binding failed')
    if not locked['locked_before_held_nonlinear']:raise ValueError('missing collar prediction lock')
    measured=arrays=None
    if (Path(directory)/'measurement.json').exists():
        measured,arrays=_read(directory,'measurement')
        if measured['prediction_json_sha256']!=sha(Path(directory)/'prediction.json'):raise ValueError('collar measurement lock binding failed')
    records=(prepared,locked) if measured is None else (prepared,locked,measured)
    if any(record['producers']!=source_hashes() for record in records):
        report.update(mode='historical_authentication_only',historical_sources_authenticated=True,
                      prediction_replayed=False,measurement_replayed=False,
                      numerical_replay_reason='current producers differ from authenticated historical producers')
        return report
    if linear_readouts(linear,locked['constants'],locked['held_w'])!=locked['prediction']:raise ValueError('collar prediction replay differs')
    report.update(prediction_replayed=True,numerical_replay=True)
    if measured is not None:
        if comparisons(locked,measured['measurement'])!=measured['comparisons']:raise ValueError('collar comparison replay differs')
        for prefix,result in (('',measured['measurement']),('fine_',measured['tighter_measurement'])):
            _verify_nonlinear(result,arrays,locked['constants'],locked['held_w'],prefix=prefix)
        report['measurement_replayed']=True;report['comparisons']=measured['comparisons']
    return report
