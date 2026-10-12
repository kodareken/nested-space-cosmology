"""One fixed-k exterior-population response on the actual coupled minus parent."""
from dataclasses import replace
from pathlib import Path
import json
import numpy as np
from threadpoolctl import threadpool_limits
from . import nsc_discovery_parent as parent
from . import nsc_discovery_parent_response as response
from . import nsc_discovery_leading_einstein as leading
from . import nsc_discovery_extent as extent
from . import provenance

LAB=parent.LAB;ROOT=parent.ROOT
INPUT=parent.OUTPUT/'strong-balanced-v1/minus'
OUTPUT=LAB/'results/development/nsc-discovery-parent-cut-response-v1'
SCHEMA='NSC-DISCOVERY-PARENT-CUT-RESPONSE-v1'
STATIONS=(0.,.25,.5,1.,1.5,2.25)
ALPHA=.05;CPU_LIMIT=300.;STEP_CAP=.001
CUTS=(3.,5.);NORMALS=(-1.,1.)
S2=np.array([[0.,-1j],[1j,0.]])


def pins():
    result=dict(parent.source_hashes(),**response.source_closure())
    for path in (Path(__file__),LAB/'scripts/derive_nsc_discovery_parent_cut_response.py',LAB/'tests/test_nsc_discovery_parent_cut_response.py',LAB/'docs/nsc-discovery-parent-cut-response.md'):
        result[str(path.relative_to(ROOT))]=parent._sha_file(path)
    return result


def reconstruct(record,a):
    grid=response.install_native_fft(parent.galerkin.build_grid(int(record['nf']),gauge='conformal'))
    grid=replace(grid,fine=replace(grid.fine,occupations=a['weights'].copy()))
    pair=parent.ParentPair(grid,a['W'],a['weights'],a['reference_columns'],a['source_columns'],
        record.get('weight_metadata',{}),dict(response.regional_map()),clock_locations=(grid.length/2,))
    state=leading.State(*(a[key].copy() for key in ('Q','r','pi_Q','pi_r','phi0','phi1')))
    return pair,state


def load_baseline():
    record,a=parent._load(INPUT)
    if record['nf']!=128 or record['sign']!=-1 or record['population']!=1 or record['fallback']['used'] or not record['converged']:
        raise ValueError('requires authenticated converged strong balanced NF128 minus without fallback')
    if abs(record['k']-.45388971484879426)>1e-14 or record['profile']['initial_radius']!='magnetic':raise ValueError('baseline physical preparation differs')
    if abs(record['source_strength']-173.16013550038755)>1e-10:raise ValueError('baseline source strength differs')
    return record,a,*reconstruct(record,a)


def direction(pair):return np.array([0.,float(pair.weights[1])])


def car(pair,weights,state=None):
    columns=pair.source_columns if state is None else np.vstack((state.phi0,state.phi1))
    eigen=parent._weighted_eigenvalues(columns.conj().T@columns,weights)
    if np.min(weights)<0 or np.min(eigen)<-1e-12 or np.max(eigen)>1+1e-10:raise ValueError('perturbed source left CAR interval')
    return eigen.tolist()


def problem(pair,state,k,sign):
    nodal=leading.decode(pair,state)
    seed=parent._offset_seed(pair,nodal.r,nodal.Q,state.phi0,state.phi1,.05);seed['suggested_k']=float(k)
    p=parent._MomentumProblem(pair,nodal.r,nodal.Q,state.phi0,state.phi1,sign,0.,seed['g'])
    seeded=parent._seed_theta(p,seed,sign)
    theta=np.r_[p.even.T@nodal.p_Q,p.even.T@nodal.p_r]
    return p,seed,seeded,theta


def reprepare(pair,state,k,sign,weights,*,cpu_limit=30.):
    """Fixed Q,r,Phi and k: solve only actual finite momenta and seed anchor."""
    car(pair,weights);held=response.with_weights(pair,weights)
    p,seed,guess,theta=problem(held,state,k,sign)
    result,report=parent.correct_momenta(p,theta,parent._cpu_time()+cpu_limit)
    if not report['converged'] or not report['hard_anchor']['satisfied']:raise ValueError('perturbed finite momentum preparation unresolved')
    return held,parent._state_from_theta(held,p,result),report


def initial_tangent(pair,state,k,sign=-1,*,width=1e-3,checks=True,cpu_limit=30.):
    started=parent._cpu_time();p,seed,guess,theta=problem(pair,state,k,sign);dc=direction(pair)
    zero=response.ParentTangent(*(np.zeros_like(getattr(state,n)) for n in leading.FIELDS),dc)
    dL,dQ,dB,_op0,_op1=response._force_and_operator(pair,state,zero)
    dCs=dL/pair.grid.dx_q;dj=dB/pair.grid.dx_q
    integrand=4*seed['g']*(seed['system'].derivative@seed['fine'].r)*dCs/(seed['fine'].r**2*seed['fine'].Q)
    primitive,removed=parent._primitive(integrand,pair.grid.length,parent.CENTER);dJ=primitive(pair.grid.xi_q)
    P=sign*np.sqrt(seed['fine'].r**3*(k+seed['J']));dP=P*dJ/(2*(k+seed['J']))
    dnode=p.even@(p.even.T@parent.galerkin.pull_geometry(pair.grid,dP))
    seededP,_=p.fine_momenta(guess);danchor=float(np.mean(2*seededP*(pair.grid.A_g@dnode)/p.fine_r**3))
    source_rows=np.r_[p.even.T@parent.galerkin.pull_geometry(pair.grid,p.fine_Q*dCs),
                      p.odd.T@parent.galerkin.pull_geometry(pair.grid,dj),-danchor]
    jac=p.jacobian(theta);step,rank,singular=parent._hard_anchor_step(jac,source_rows,p.whitener,p.row_scales(theta))
    dtheta=p.whitener@step;dp,dv=p.nodal(dtheta)
    encoded=leading.encode(pair,leading.State(np.zeros(pair.grid.ng),np.zeros(pair.grid.ng),dp,dv,
                                           np.zeros_like(state.phi0),np.zeros_like(state.phi1)))
    tangent=response.ParentTangent(*(getattr(encoded,n) for n in leading.FIELDS),dc)
    report={'method':'analytic implicit actual finite momentum C/D plus fixed-k seed anchor',
        'rank':rank,'required_rank':theta.size,'analytic_rank_usable':rank==theta.size,
        'linear_residual_max':parent._maximum(jac@dtheta+source_rows),'baseline_root_residual_max':parent._maximum(p.residual(theta)),
        'delta_seed_anchor':danchor,'fixed_k_is_not_fixed_anchor':True,'delta_J_max':parent._maximum(dJ),
        'primitive_removed':removed,'direction':dc.tolist(),'Q_r_Phi_fixed':True,'finite_difference_checks':[]}
    if rank!=theta.size:
        report['scope']='unresolved truncated initial derivative; not an analytic full-rank preparation'
    if checks and rank==theta.size:
        for h in (width,width/2):
            try:
                plus,up,_=reprepare(pair,state,k,sign,pair.weights+h*dc,cpu_limit=cpu_limit-(parent._cpu_time()-started))
                minus,down,_=reprepare(pair,state,k,sign,pair.weights-h*dc,cpu_limit=cpu_limit-(parent._cpu_time()-started))
                report['finite_difference_checks'].append({'width':h,'status':'RESOLVED','momentum_max_gap':max(parent._maximum((getattr(up,n)-getattr(down,n))/(2*h)-getattr(tangent,n)) for n in ('p_Q','p_r'))})
            except ValueError as error:
                report['finite_difference_checks'].append({'width':h,'status':'UNRESOLVED','reason':str(error)})
    return tangent,report


def ap_values(phi,length,x):
    """AP Fourier cut values in continuum density units, no periodic wrap."""
    nf=len(phi);m=np.fft.fftfreq(nf)*nf
    coeff=np.fft.fft(phi*np.exp(-1j*np.pi*np.arange(nf)/nf)[:,None],axis=0)/nf
    return np.exp(2j*np.pi*np.asarray(x)[:,None]*(m[None,:]+.5)/length)@coeff/np.sqrt(length/nf)


def observables(pair,state,tangent=None):
    rate,bundle=leading.rates(pair,state,return_bundle=True);fine=bundle['fine_state']
    velocity=leading.prolong(pair.grid,leading.decode(pair,rate));D=bundle['fine_system'].derivative
    ev=lambda a:extent.real_periodic_values(pair.grid,a,CUTS)
    N=ev(fine.r*fine.Q);r=ev(fine.r);Q=ev(fine.Q)
    ut={n:ev(getattr(velocity,n)) for n in ('r','Q')};ux={n:ev(D@getattr(fine,n)) for n in ('r','Q')}
    psi=np.stack((ap_values(state.phi0,pair.grid.length,CUTS),ap_values(state.phi1,pair.grid.length,CUTS)),axis=1)
    values={'child_content':response.child_regional_content(pair,state),
            'child_proper_length':extent.real_interval_integral(pair.grid,fine.r*fine.Q,pair.child_interval),
            'centre_clock_rate':response.centre_clock(pair,state)};delta={}
    if tangent is not None:
        lifted=leading.prolong(pair.grid,leading.decode(pair,response._as_state(tangent)))
        image=response.response_rates(pair,state,tangent);dvel=leading.prolong(pair.grid,leading.decode(pair,response._as_state(image)))
        dN=ev(lifted.r*fine.Q+fine.r*lifted.Q)
        dpsi=np.stack((ap_values(tangent.phi0,pair.grid.length,CUTS),ap_values(tangent.phi1,pair.grid.length,CUTS)),axis=1)
        delta['child_content']=response.child_regional_content(pair,state,tangent)[1]
        delta['child_proper_length']=extent.real_interval_integral(pair.grid,lifted.r*fine.Q+fine.r*lifted.Q,pair.child_interval)
        delta['centre_clock_rate']=response.centre_clock(pair,state,tangent)[1]
    for i,(x,n) in enumerate(zip(CUTS,NORMALS)):
        label='left' if n<0 else 'right';incoming_r=ut['r'][i]+n*ux['r'][i];incoming_Q=ut['Q'][i]+n*ux['Q'][i]
        projected=(np.eye(2)-n*S2)/2@psi[i];power=np.sum(abs(projected)**2,axis=0);incoming=float(np.dot(pair.weights,power))
        values.update({label+'_r_incoming_normal':float(incoming_r/N[i]),label+'_Q_incoming_coordinate':float(incoming_Q),
                       label+'_Dirac_incoming_normal_current':incoming/N[i]})
        if tangent is not None:
            dr=ev(dvel.r)[i]+n*ev(D@lifted.r)[i];dq=ev(dvel.Q)[i]+n*ev(D@lifted.Q)[i]
            dprojected=(np.eye(2)-n*S2)/2@dpsi[i]
            dpower=2*np.real(np.sum(projected.conj()*dprojected,axis=0))
            di=float(np.dot(tangent.occupations,power)+np.dot(pair.weights,dpower))
            delta.update({label+'_r_incoming_normal':float(dr/N[i]-incoming_r*dN[i]/N[i]**2),
                label+'_Q_incoming_coordinate':float(dq),label+'_Dirac_incoming_normal_current':di/N[i]-incoming*dN[i]/N[i]**2})
    return values,delta


def row(pair,state,tangent,time,tau,dtau):
    values,delta=observables(pair,state,tangent)
    _,slope=observables(pair,state,response.velocity_tangent(pair,state))
    fixed={key:response.centre_clock_variation(delta[key],slope[key],dtau,values['centre_clock_rate']) for key in values}
    return {'time':time,'tau':tau,'delta_tau':dtau,'values':values,'delta_t':delta,'delta_equal_tau':fixed,
            'constraints':leading.constraints(pair,state),'CAR':car(pair,pair.weights,state)}


def evolve(pair,state,*,tangent=None,stations=STATIONS,proper_target=None,cpu_limit=300.):
    """One serial bounded loop using existing admission/shared stage steppers."""
    start=parent._cpu_time();time=tau=dtau=0.;rows=[];steps=0;stop='COMPLETE';index=0
    while True:
        if proper_target is None and index<len(stations) and abs(time-stations[index])<1e-11:
            rows.append(row(pair,state,tangent,time,tau,dtau));index+=1
        if (proper_target is None and index==len(stations)) or (proper_target is not None and tau>=proper_target-1e-11):break
        if parent._cpu_time()-start>cpu_limit-1.:stop='CPU_LIMIT';break
        dt,_=response._admitted_step(pair,state,STEP_CAP)
        if proper_target is None:dt=min(dt,stations[index]-time)
        else:
            full=response._nonlinear_centre_increment(pair,state,dt)
            if tau+full>proper_target:
                lo=0.;hi=dt
                for _ in range(32):
                    mid=(lo+hi)/2
                    if tau+response._nonlinear_centre_increment(pair,state,mid)<proper_target:lo=mid
                    else:hi=mid
                dt=(lo+hi)/2
        try:
            if tangent is not None:
                state,tangent,clock=response.advance_tangent(pair,state,tangent,dt);tau+=clock['tau'];dtau+=clock['delta_tau']
            else:
                increment=response._nonlinear_centre_increment(pair,state,dt);state=leading.rk4_step(pair,state,dt);tau+=increment
        except leading.coupling.PositiveChartExit:stop='POSITIVE_CHART_EXIT';break
        time+=dt;steps+=1
    return state,tangent,{'rows':rows,'time':time,'tau':tau,'delta_tau':dtau,'steps':steps,'stop':stop,'CPU_seconds':parent._cpu_time()-start}


def _arrays(pair,state,tangent=None):
    a={n:getattr(state,n) for n in leading.FIELDS};a.update(weights=pair.weights,W=pair.geometry_map,
          reference_columns=pair.reference_columns,source_columns=pair.source_columns)
    if tangent is not None:a.update({'delta_'+n:getattr(tangent,n) for n in leading.FIELDS},delta_c=tangent.occupations)
    return a


def _state(a):return leading.State(*(a[n].copy() for n in leading.FIELDS))

def _stage(directory,name,record,a):
    if record['producers']!=pins() or not record.get('producing_commit'):raise ValueError('explicit stage needs frozen current producers')
    for path,digest in record['producers'].items():provenance.resolve_pinned_source_bytes(ROOT,path,digest,commit=record['producing_commit'])
    for path,digest in record['inputs'].items():
        if parent._sha_file(path)!=digest:raise ValueError('response input changed before write')
    return response._write_exclusive(directory,name,dict(record,source_closure=response.source_closure(),frozen_producer=True,physical_binding=True),a,input_paths=(INPUT,))


def _read(directory,name):
    record,a=response._load_stage(directory,name)
    for path,digest in record['producers'].items():provenance.resolve_pinned_source_bytes(ROOT,path,digest,commit=record['producing_commit'])
    if record['producers']!=pins():raise ValueError('current source differs; historical authentication only, no numerical replay')
    for path,digest in record['inputs'].items():
        if parent._sha_file(path)!=digest:raise ValueError('response input binding changed')
    return record,a


def preview():
    record,a,pair,state=load_baseline()
    return {'schema':SCHEMA,'mode':'preview','evolved':False,'input':str(INPUT),'direction':direction(pair).tolist(),'stations':STATIONS,
            'cuts':CUTS,'outward_normals':NORMALS,'child_window':pair.child_interval,'alpha':ALPHA,'k_fixed':record['k']}


def prepare(directory):
    start=parent._cpu_time();closure=pins();commit=parent.producing_commit(closure)
    if not commit:raise ValueError('prepare requires frozen producer')
    record,a,pair,state=load_baseline();tangent,derivative=initial_tangent(pair,state,record['k'],record['sign'])
    for sign in (-1,1):car(pair,pair.weights+sign*ALPHA*direction(pair))
    inputs={str(INPUT/name):parent._sha_file(INPUT/name) for name in ('parent.json','parent.npz')}
    values= row(pair,state,tangent,0.,0.,0.)
    return _stage(directory,'prepare',{'schema':SCHEMA,'mode':'prepared','producers':closure,'producing_commit':commit,'inputs':inputs,
       'baseline_record':record,'initial_derivative':derivative,'initial_row':values,'alpha':ALPHA,'CPU_seconds':parent._cpu_time()-start,
       'declared_held_amplitudes':[-ALPHA,ALPHA],
       'state_hashes':{n:parent._sha_array(getattr(state,n)) for n in leading.FIELDS},'direction_hash':parent._sha_array(tangent.occupations),
       'finite_band_tails_and_initial_constraint_response_retained':True,'exact_compact_delay_claimed':False},_arrays(pair,state,tangent))


def predict(directory):
    start=parent._cpu_time()
    prepared,a=_read(directory,'prepare');pair,_=reconstruct(prepared['baseline_record'],dict(a,pi_Q=a['p_Q'],pi_r=a['p_r']))
    if not prepared['initial_derivative']['analytic_rank_usable'] or any(item['status']!='RESOLVED' for item in prepared['initial_derivative']['finite_difference_checks']):
        raise ValueError('initial derivative/preparation scope unresolved; forecast not authorized')
    tangent=response.ParentTangent(*(a['delta_'+n] for n in leading.FIELDS),a['delta_c'])
    state,tangent,forecast=evolve(pair,_state(a),tangent=tangent,cpu_limit=CPU_LIMIT-prepared['CPU_seconds']-(parent._cpu_time()-start))
    for item in forecast['rows']:
        item['predicted_held_equal_tau']={label:{n:item['values'][n]+sign*ALPHA*item['delta_equal_tau'][n] for n in item['values']}
                                        for label,sign in (('plus',1),('minus',-1))}
        item['evolved_directional_change']={n:item['delta_t'][n]-prepared['initial_row']['delta_t'][n] for n in item['values']}
    return _stage(directory,'prediction',dict(prepared,mode='locked_prediction',forecast_locked=True,forecast=forecast,
        prepare_json_sha256=parent._sha_file(Path(directory)/'prepare.json'),CPU_seconds=prepared['CPU_seconds']+parent._cpu_time()-start),_arrays(pair,state,tangent))


def run(directory,*,station=2.25,amplitude=ALPHA):
    start=parent._cpu_time()
    if not np.isfinite(amplitude):raise ValueError('measurement amplitude must be finite')
    name=f'measurement-a{amplitude:g}-T{station:g}'
    if any((Path(directory)/(name+suffix)).exists() for suffix in ('.json','.npz')):raise FileExistsError('immutable measurement already exists')
    locked,_a=_read(directory,'prediction');prepared,a=_read(directory,'prepare')
    if not locked.get('forecast_locked') or parent._sha_file(Path(directory)/'prepare.json')!=locked['prepare_json_sha256']:raise ValueError('held run requires unchanged locked forecast')
    prediction_hash=parent._sha_file(Path(directory)/'prediction.json')
    target=next((item for item in locked['forecast']['rows'] if abs(item['time']-station)<1e-11),None)
    if target is None:raise ValueError('requested station not reached by locked forecast')
    used=sum(json.loads(path.read_text())['run_CPU_seconds'] for path in Path(directory).glob('measurement-*.json'))
    pair,state=reconstruct(prepared['baseline_record'],dict(a,pi_Q=a['p_Q'],pi_r=a['p_r']))
    remaining=CPU_LIMIT-locked['CPU_seconds']-used-(parent._cpu_time()-start)
    if remaining<=1:raise ValueError('aggregate batch CPU budget exhausted before held preparation')
    pair,state,prep=reprepare(pair,state,prepared['baseline_record']['k'],-1,pair.weights+amplitude*a['delta_c'],cpu_limit=min(30.,remaining-1))
    state,_t,measured=evolve(pair,state,proper_target=target['tau'],cpu_limit=CPU_LIMIT-locked['CPU_seconds']-used-(parent._cpu_time()-start))
    elapsed=parent._cpu_time()-start;aggregate=locked['CPU_seconds']+used+elapsed
    if aggregate>CPU_LIMIT:measured['stop']='AGGREGATE_CPU_LIMIT'
    values,_=observables(pair,state);matched=measured['stop']=='COMPLETE' and abs(measured['tau']-target['tau'])<1e-9
    comparison={n:{'measured_change':values[n]-target['values'][n],
          'prediction':amplitude*target['delta_equal_tau'][n],
          'residual':values[n]-target['values'][n]-amplitude*target['delta_equal_tau'][n]} for n in values} if matched else None
    if parent._sha_file(Path(directory)/'prediction.json')!=prediction_hash:raise ValueError('forecast changed during held run')
    return _stage(directory,name,dict(locked,mode='held_nonlinear_measurement',prediction_json_sha256=prediction_hash,
        amplitude=amplitude,station=station,measurement=measured,comparison=comparison,held_preparation=prep,
        values=values,matched_proper_clock=matched,comparison_valid=matched,aggregate_CPU_seconds=aggregate,
        run_CPU_seconds=parent._cpu_time()-start,finite_difference_control=abs(amplitude)!=ALPHA),_arrays(pair,state))
