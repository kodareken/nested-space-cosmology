"""Conditional prescribed-pulse work on the actual prepared uniform CAR source.

Local AP momentum blocks replace neither the gravity equations nor the state
interpretation. M is outside the column Hamiltonian. The filled-negative
reference is a mathematical frozen-field control, not a physical vacuum.
"""
from __future__ import annotations

import hashlib
import io
import json
import os
from pathlib import Path
import re
import subprocess
import time

import numpy as np
from scipy.integrate import solve_ivp

from . import nsc_discovery_dynamic_preparation as preparation
from . import nsc_spherical_coupling as coupling
from . import nsc_spherical_galerkin_coupling as galerkin

LAB=Path(__file__).resolve().parents[2]
INPUT=LAB/'results/development/nsc-discovery-dynamic-preparation-v2'
OUTPUT=LAB/'results/development/nsc-discovery-spectral-work-v1'
SCHEMA='NSC-DISCOVERY-SPECTRAL-WORK-v1'
CPU_LIMIT=30.
CHUNK_LIMIT=64*1024*1024
SIGMA1=np.array([[0.,1.],[1.,0.]],complex)
SIGMA2=np.array([[0.,-1j],[1j,0.]],complex)


def sha256(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def producer_hashes():
    paths=[Path(__file__),LAB/'scripts/derive_nsc_discovery_spectral_work.py',Path(preparation.__file__),
           Path(coupling.__file__),Path(galerkin.__file__),LAB/'src/recursive_horizons/nsc_spherical_feedback_action.py']
    return {str(path.relative_to(LAB)):sha256(path) for path in paths}


def _git_hash(commit,path):
    if not re.fullmatch(r'[0-9a-f]{40}',str(commit)) or Path(path).is_absolute() or '..' in Path(path).parts:
        raise ValueError('invalid frozen producer commit/path')
    blob=subprocess.run(['git','cat-file','blob',f'{commit}:lab/{path}'],cwd=LAB.parent,capture_output=True,check=True,timeout=5).stdout
    return hashlib.sha256(blob).hexdigest()


def producing_commit(pins):
    try:
        commit=subprocess.run(['git','rev-parse','HEAD'],cwd=LAB.parent,capture_output=True,check=True,timeout=5).stdout.decode().strip()
        if all(_git_hash(commit,path)==digest for path,digest in pins.items()):return commit
    except (ValueError,OSError,subprocess.SubprocessError):pass
    return None


def _input_check(hashes):
    for path,digest in hashes.items():
        if sha256(LAB/path)!=digest:raise ValueError('sealed pulse input changed: '+path)


def pulse_parameters(gap,*,amplitude=.001,frequency_ratio=.8):
    gap=float(gap);amplitude=float(amplitude);frequency_ratio=float(frequency_ratio)
    if gap<=0. or not 0.<frequency_ratio<1. or not np.isfinite(amplitude) or amplitude<=0.:
        raise ValueError('positive gap/amplitude and sub-gap frequency ratio required')
    omega=frequency_ratio*gap
    return {'amplitude':amplitude,'omega':omega,'T':4.*np.pi/omega,'frequency_ratio':frequency_ratio,
            'shape':'A sin^2(pi t/T) cos(omega t)','endpoint_J':0.,'endpoint_Jdot':0.,
            'analytic_integral_J':0.,'coordinate_pulse':True}


def pulse(time_value,settings):
    t=float(time_value);T=settings['T'];A=settings['amplitude'];omega=settings['omega']
    if t<=0. or t>=T:return 0.,0.
    s=np.sin(np.pi*t/T);c=np.cos(np.pi*t/T)
    value=A*s*s*np.cos(omega*t)
    slope=A*(2.*s*c*np.pi/T*np.cos(omega*t)-s*s*omega*np.sin(omega*t))
    return float(value),float(slope)


def pulse_transform(gap,settings):
    """Integral J(t) exp(i gap t), stable also at coincident frequencies."""
    T=settings['T'];omega=settings['omega'];A=settings['amplitude'];envelope=2.*np.pi/T
    integral=lambda frequency:T*np.exp(.5j*frequency*T)*np.sinc(frequency*T/(2.*np.pi))
    value=A/4.*(integral(gap+omega)+integral(gap-omega))
    value-=A/8.*sum(integral(gap+sign*omega+side*envelope) for sign in (-1.,1.) for side in (-1.,1.))
    return complex(value)


def block_operators(momentum,Q,kappa=1.):
    H=float(momentum)*SIGMA2+float(kappa)*float(Q)*SIGMA1
    eigenvalues,vectors=np.linalg.eigh(H)
    G=float(kappa)*SIGMA1
    return eigenvalues,vectors,vectors.conj().T@G@vectors


def _modal_columns(phi0,phi1):
    nf=phi0.shape[0]
    demod=np.exp(-1j*np.pi*np.arange(nf)/nf)[:,None]
    return np.stack((np.fft.fft(phi0*demod,axis=0)/np.sqrt(nf),
                     np.fft.fft(phi1*demod,axis=0)/np.sqrt(nf)),axis=1)


def _noise_and_characteristic(columns_by_mode,weights,momenta,Q,kappa,*,angle=.2):
    """Rank-six determinant update; no dense ambient C or propagator."""
    dim=columns_by_mode.shape[-1];I=np.eye(2,dtype=complex)
    U=np.cos(kappa*angle)*I+1j*np.sin(kappa*angle)*SIGMA1
    G=kappa*SIGMA1
    gramG=sum(B.conj().T@G@B for B in columns_by_mode)
    full_noise=float(kappa*kappa*np.sum(weights)-np.sum(weights[:,None]*weights[None,:]*np.abs(gramG)**2))
    full_update=sum(B.conj().T@(U-I)@B for B in columns_by_mode)
    reference_update=np.zeros((dim,dim),complex)
    vacuum_noise=0.;cross_noise=0.
    for p,B in zip(momenta,columns_by_mode):
        levels,V,Ge=block_operators(p,Q,kappa)
        lower=V[:,0:1];P=lower@lower.conj().T
        vacuum_noise+=abs(Ge[0,1])**2
        C=(B*weights[None,:])@B.conj().T
        cross_noise+=float(np.trace(P@G@C@G).real)
        Avac=I-P+P@U
        reference_update+=B.conj().T@(U-I)@np.linalg.solve(Avac,B)
    Zfull=np.linalg.det(np.eye(dim)+weights[:,None]*full_update)
    Zratio=np.linalg.det(np.eye(dim)+weights[:,None]*reference_update)
    difference_noise=full_noise-2.*cross_noise
    return {'full_C_noise':full_noise,'vacuum_reference_noise':float(vacuum_noise),
            'reference_minus_vacuum_noise':float(difference_noise),
            'noise_subtraction_minus_full':float(difference_noise-full_noise),
            'characteristic_angle':angle,'Zfull_real':float(Zfull.real),'Zfull_imag':float(Zfull.imag),
            'Zreference_over_Zvacuum_real':float(Zratio.real),'Zreference_over_Zvacuum_imag':float(Zratio.imag),
            'characteristic_ratio_minus_full_abs':float(abs(Zratio-Zfull)),
            'normalized_Gaussian_ratio_is_excitation_Gaussian':False,
            'noise_and_characteristic_scope':'one representative canonical CAR sector; M multiplies the mean/action work outside H',
            'covariance_interpretation':'full finite CAR versus mathematical filled-negative reference; physical filled-sea Gamma matching unresolved'}


def load_uniform():
    binding=preparation.authenticate_preparation(INPUT)
    jp=INPUT.with_suffix('.json');npz=INPUT.with_suffix('.npz');record=json.loads(jp.read_text())
    with np.load(npz,allow_pickle=False) as saved:
        Q=np.array(saved['uniform_nodal_Q'],copy=True);r=np.array(saved['uniform_nodal_r'],copy=True)
        phi0=np.array(saved['uniform_nodal_phi0'],copy=True);phi1=np.array(saved['uniform_nodal_phi1'],copy=True)
        stored_levels=np.array(saved['uniform_eigenvalues'],copy=True)
    Q0=float(np.mean(Q));r0=float(np.mean(r));nf=int(record['nf']);period=float(record['period'])
    weights=np.asarray(record['occupations'],float);M=float(record['multiplicity']);kappa=M/4.
    if record['gauge']!='conformal' or np.ptp(Q)>1e-11 or np.ptp(r)>1e-10 or Q0<=0. or r0<=0.:
        raise ValueError('pulse control requires the actual positive uniform conformal preparation')
    if kappa!=1. or M!=4. or weights.shape!=(6,):raise ValueError('prepared source is not the declared kappa=1, M=4 rank-six sector')
    B=_modal_columns(phi0,phi1)
    momenta=2.*np.pi/period*(np.fft.fftfreq(nf)*nf+.5)
    populations=[];offdiag=0.;negative=0.;eigen_defect=0.
    for p,b in zip(momenta,B):
        levels,V,G=block_operators(p,Q0,kappa)
        amplitudes=V.conj().T@b
        C=(amplitudes*weights[None,:])@amplitudes.conj().T
        populations.append(float(C[1,1].real));negative+=float(max(0.,C[0,0].real))
        offdiag=max(offdiag,float(abs(C[0,1])))
        eigen_defect=max(eigen_defect,float(np.max(np.abs((p*SIGMA2+kappa*Q0*SIGMA1)@b-b*stored_levels[None,:]))))
    populations=np.asarray(populations)
    active=np.flatnonzero(populations>1e-10)
    cross=max(float(np.max(np.abs((B[a]*weights)@B[b].conj().T))) for a in active for b in active if a!=b)
    positive=np.flatnonzero(momenta>0.);order=positive[np.argsort(momenta[positive])]
    p=momenta[order];cp=populations[order]
    for index in order:
        opposite=int(np.argmin(np.abs(momenta+momenta[index])))
        if abs(populations[index]-populations[opposite])>1e-10:raise ValueError('source breaks the paired AP momentum occupations')
    if len(active)!=6 or np.max(np.abs(cp[:3]-weights[::2]))>1e-10 or np.max(cp[3:])>1e-10:
        raise ValueError('actual uniform source does not have the expected six positive-only occupied modes')
    energies=[];Gblocks=[]
    for momentum in p:
        levels,V,G=block_operators(momentum,Q0,kappa);energies.append(levels);Gblocks.append(G)
    energies=np.asarray(energies);Gblocks=np.asarray(Gblocks)
    source_energy=M*float(np.dot(weights,stored_levels))
    block_energy=M*float(2.*np.dot(cp,energies[:,1]))
    if eigen_defect>1e-9 or abs(source_energy-block_energy)>1e-9:raise ValueError('actual source/operator block reconstruction failed')
    hashes={str(jp.relative_to(LAB)):sha256(jp),str(npz.relative_to(LAB)):sha256(npz)}
    hashes.update({str(Path(path).relative_to('lab')):digest for path,digest in record['input_hashes'].items() if Path(path).parts[0]=='lab'})
    return {'nf':nf,'period':period,'Q0':Q0,'r0':r0,'N0':r0*Q0,'kappa':kappa,'M':M,
            'momenta':p,'energies':energies,'Gblocks':Gblocks,'positive_occupations':cp,'degeneracy':np.full(len(p),2.),
            'input_hashes':hashes,'preparation_binding':binding,
            'source_checks':{'Gram_gap':float(np.max(np.abs(np.vstack((phi0,phi1)).conj().T@np.vstack((phi0,phi1))-np.eye(6)))),
                'negative_occupation_roundoff':negative,'spectral_offdiagonal_max':offdiag,'cross_momentum_covariance_max':cross,
                'eigen_equation_max':eigen_defect,'source_energy':source_energy,'block_energy':block_energy,
                'uniform_Q_spread':float(np.ptp(Q)),'uniform_r_spread':float(np.ptp(r))},
            'Gaussian_reference_gap':_noise_and_characteristic(B,weights,momenta,Q0,kappa)}


def forecast(data,settings=None):
    gaps=data['energies'][:,1]-data['energies'][:,0]
    settings=pulse_parameters(gaps[0]) if settings is None else dict(settings)
    if data['Q0']-settings['amplitude']<=0.:raise ValueError('prescribed Q pulse can leave the positive chart')
    matrix2=np.abs(data['Gblocks'][:,0,1])**2
    transform=np.array([pulse_transform(gap,settings) for gap in gaps])
    transitions=matrix2*np.abs(transform)**2
    c=data['positive_occupations'];d=data['degeneracy'];M=data['M']
    full=M*float(np.sum(d*(-c)*gaps*transitions))
    vacuum=M*float(np.sum(d*gaps*transitions))
    reference=M*float(np.sum(d*(1.-c)*gaps*transitions))
    channels=[{'momentum':float(data['momenta'][i]),'gap':float(gaps[i]),'positive_occupation':float(c[i]),
               'degeneracy':int(d[i]),'G_lower_upper_squared':float(matrix2[i]),
               'positive_frequency_residue_per_channel':float(-c[i]*matrix2[i]),
               'reference_positive_frequency_residue_per_channel':float((1.-c[i])*matrix2[i]),
               'Jtilde_real':float(transform[i].real),'Jtilde_imag':float(transform[i].imag),
               'second_order_transition_probability':float(transitions[i])} for i in range(len(gaps))]
    return {'pulse':settings,'channels':channels,'second_order_work_full_C':full,'second_order_work_passive_reference':reference,
            'second_order_work_vacuum_reference':vacuum,'mean_reference_subtraction_gap':reference-vacuum-full,
            'reference_proper_endpoint':data['N0']*settings['T'],'pulse_proper_endpoint':data['N0']*settings['T'],
            'endpoint_clocks_match_by_integral_J_zero':True,'clocks_are_not_equal_at_all_intermediate_times':True,
            'formula':'M sum_{a<b}(c_a-c_b) Delta_epsilon |G_ab|^2 |Jtilde(Delta_epsilon)|^2',
            'M_outside_H_and_column_law':True,'no_pole_evaluated':True,'nonlinear_measurement_performed':False,
            'scope':'prescribed uniform frozen-field pulse; gain is conditional on the declared full finite CAR covariance',
            'physical_filled_sea_Gamma_matching':'unresolved','autonomous_regeneration_claimed':False}


def preview():
    data=load_uniform();predicted=forecast(data)
    return {'schema':SCHEMA,'mode':'preview','evolved':False,'input_hashes':data['input_hashes'],
            'Q0':data['Q0'],'r0':data['r0'],'kappa':data['kappa'],'M':data['M'],'source_checks':data['source_checks'],
            'prediction':predicted,'Gaussian_reference_gap':data['Gaussian_reference_gap'],'CPU_limit_seconds':CPU_LIMIT}


def evolve_blocks(energies,Gblocks,settings,*,rtol=2e-11,atol=2e-13,cpu_limit=CPU_LIMIT):
    """Interaction-picture DOP853 for local 2x2 blocks and centered work integrals."""
    started=time.process_time();energies=np.asarray(energies,float);Gblocks=np.asarray(Gblocks,complex)
    count=len(energies);initial=np.zeros(6*count+1,complex)
    initial[:4*count]=np.tile(np.eye(2,dtype=complex),(count,1,1)).ravel()
    initial_means=np.diagonal(Gblocks,axis1=1,axis2=2).real
    differences=energies[:,:,None]-energies[:,None,:]
    def rhs(t,y):
        if time.process_time()-started>cpu_limit:raise RuntimeError('spectral pulse exceeded its 30 CPU-second cap')
        J,Jdot=pulse(t,settings)
        GI=Gblocks*np.exp(1j*differences*t)
        U=y[:4*count].reshape(count,2,2)
        derivative=-1j*J*np.einsum('nab,nbc->nac',GI,U)
        means=np.einsum('nai,nab,nbi->ni',U.conj(),GI,U).real
        work=(means-initial_means)*Jdot
        return np.concatenate((derivative.ravel(),work.ravel(),np.array([J],complex)))
    max_gap=float(np.max(energies[:,1]-energies[:,0]))
    max_step=min(settings['T']/32.,np.pi/(4.*max_gap))
    solved=solve_ivp(rhs,(0.,settings['T']),initial,method='DOP853',rtol=rtol,atol=atol,max_step=max_step,t_eval=[settings['T']])
    if not solved.success:raise RuntimeError('spectral pulse integration failed: '+solved.message)
    final=solved.y[:,-1];U=final[:4*count].reshape(count,2,2);work=final[4*count:6*count].reshape(count,2).real
    transition_up=np.abs(U[:,1,0])**2;transition_down=np.abs(U[:,0,1])**2
    gram=np.einsum('nai,naj->nij',U.conj(),U)
    return {'U_interaction_endpoint':U,'integrated_work_per_initial_branch':work,
            'transition_up':transition_up,'transition_down':transition_down,
            'unitarity_max':float(np.max(np.abs(gram-np.eye(2)))),
            'up_down_probability_gap':float(np.max(np.abs(transition_up-transition_down))),
            'pulse_area_integral':float(final[-1].real),
            'nfev':solved.nfev,'CPU_seconds':time.process_time()-started,'max_step':max_step,'rtol':rtol,'atol':atol}


def _write(directory,prefix,record,arrays):
    directory=Path(directory);jp=directory/(prefix+'.json');npz=directory/(prefix+'.npz')
    if jp.exists() or npz.exists():raise FileExistsError('immutable pulse prefix already exists: '+prefix)
    stream=io.BytesIO();np.savez_compressed(stream,**arrays);blob=stream.getvalue()
    bound=dict(record,payload_sha256=hashlib.sha256(blob).hexdigest(),array_sha256={key:hashlib.sha256(np.ascontiguousarray(value).tobytes()).hexdigest() for key,value in arrays.items()})
    text=(json.dumps(bound,indent=2,sort_keys=True,allow_nan=False)+'\n').encode()
    if len(blob)+len(text)>CHUNK_LIMIT:raise RuntimeError('spectral pulse record exceeds 64 MiB')
    directory.mkdir(parents=True,exist_ok=True)
    for path,content in ((npz,blob),(jp,text)):
        descriptor=os.open(path,os.O_CREAT|os.O_EXCL|os.O_WRONLY,0o444)
        with os.fdopen(descriptor,'wb') as output:output.write(content)
    return bound


def _read(directory,prefix):
    directory=Path(directory);record=json.loads((directory/(prefix+'.json')).read_text());npz=directory/(prefix+'.npz')
    if sha256(npz)!=record['payload_sha256']:raise ValueError('pulse payload binding failed')
    with np.load(npz,allow_pickle=False) as saved:arrays={key:saved[key].copy() for key in saved.files}
    for key,digest in record['array_sha256'].items():
        if hashlib.sha256(np.ascontiguousarray(arrays[key]).tobytes()).hexdigest()!=digest:raise ValueError('pulse array binding failed')
    _input_check(record['input_hashes'])
    if record['producers']!=producer_hashes():
        commit=record.get('producing_commit')
        if not commit or any(_git_hash(commit,path)!=digest for path,digest in record['producers'].items()):
            raise ValueError('frozen spectral work producers do not authenticate')
    return record,arrays


def prepare(directory):
    if (Path(directory)/'prepare.json').exists() or (Path(directory)/'prepare.npz').exists():raise FileExistsError('pulse preparation already exists')
    started=time.process_time();pins=producer_hashes();data=load_uniform();predicted=forecast(data)
    arrays={key:data[key] for key in ('momenta','energies','Gblocks','positive_occupations','degeneracy')}
    _input_check(data['input_hashes'])
    if pins!=producer_hashes():raise ValueError('spectral work source changed during prediction')
    return _write(directory,'prepare',{'schema':SCHEMA,'mode':'frozen_source_operator','producers':pins,'producing_commit':producing_commit(pins),
           'input_hashes':data['input_hashes'],'preparation_binding':data['preparation_binding'],'Q0':data['Q0'],'r0':data['r0'],
           'N0':data['N0'],'kappa':data['kappa'],'M':data['M'],'source_checks':data['source_checks'],'pulse':predicted['pulse'],
           'Gaussian_reference_gap':data['Gaussian_reference_gap'],'CPU_seconds':time.process_time()-started,
           'source_operator_frozen_before_prediction':True,'nonlinear_measurement_performed':False},arrays)


def predict(directory):
    if (Path(directory)/'prediction.json').exists() or (Path(directory)/'prediction.npz').exists():raise FileExistsError('pulse prediction already exists')
    started=time.process_time();prepared,arrays=_read(directory,'prepare')
    pins=producer_hashes()
    if prepared['producers']!=pins:raise ValueError('prediction requires current frozen source/operator producers')
    predicted=forecast(dict(arrays,Q0=prepared['Q0'],M=prepared['M'],N0=prepared['N0']),prepared['pulse'])
    record={key:value for key,value in prepared.items() if key not in ('payload_sha256','array_sha256','source_operator_frozen_before_prediction')}
    record.update(mode='locked_prediction',prediction=predicted,prepare_json_sha256=sha256(Path(directory)/'prepare.json'),
                  locked_before_nonlinear_pulse=True,CPU_seconds=prepared['CPU_seconds']+time.process_time()-started)
    return _write(directory,'prediction',record,arrays)


def measurement(data,endpoint):
    gaps=data['energies'][:,1]-data['energies'][:,0];c=data['positive_occupations'];d=data['degeneracy'];M=data['M']
    down=endpoint['transition_down'];up=endpoint['transition_up'];work=endpoint['integrated_work_per_initial_branch']
    full=-M*float(np.sum(d*c*gaps*down));vacuum=M*float(np.sum(d*gaps*up));reference=vacuum+full
    integrated_full=M*float(np.sum(d*c*work[:,1]));integrated_vacuum=M*float(np.sum(d*work[:,0]))
    integrated_reference=integrated_vacuum+integrated_full
    U=endpoint['U_interaction_endpoint'];CAR=[]
    for matrix,occupation in zip(U,c):
        CAR.extend(np.linalg.eigvalsh((matrix*np.array([0.,occupation])[None,:])@matrix.conj().T).tolist())
        CAR.extend(np.linalg.eigvalsh((matrix*np.array([1.,occupation])[None,:])@matrix.conj().T).tolist())
    return {'work_full_C':full,'work_passive_reference':reference,'work_vacuum_reference':vacuum,
            'integrated_work_full_C':integrated_full,'integrated_work_passive_reference':integrated_reference,
            'integrated_work_vacuum_reference':integrated_vacuum,
            'energy_work_closure_full':integrated_full-full,'energy_work_closure_reference':integrated_reference-reference,
            'mean_reference_subtraction_gap':reference-vacuum-full,
            'CAR_min':float(min(CAR)),'CAR_max':float(max(CAR)),'CAR_preserved_with_unitarity_error':endpoint['unitarity_max'],
            'energy_from_transition_probabilities_not_large_subtraction':True,
            'integrated_work_definition':'M integral Tr C(t) G Jdot; initial constant force integrates exactly to its zero endpoint term',
            'passive_reference_not_fed_to_gravity':True,'physical_vacuum_claimed':False,
            'active_source_is_not_a_universal_model_or_antimatter_verdict':True,'autonomous_regeneration_claimed':False}


def run(directory):
    if (Path(directory)/'measurement.json').exists() or (Path(directory)/'measurement.npz').exists():raise FileExistsError('pulse measurement already exists')
    started=time.process_time();locked,arrays=_read(directory,'prediction');pins=producer_hashes()
    if pins!=locked['producers']:raise ValueError('nonlinear pulse needs current frozen producer bytes')
    remaining=CPU_LIMIT-locked['CPU_seconds']-(time.process_time()-started)
    endpoint=evolve_blocks(arrays['energies'],arrays['Gblocks'],locked['prediction']['pulse'],cpu_limit=remaining)
    measured=measurement(dict(arrays,M=locked['M']),endpoint)
    predicted=locked['prediction']
    for state in ('full_C','passive_reference','vacuum_reference'):
        measured['prediction_error_'+state]=measured['work_'+state]-predicted['second_order_work_'+state]
    settings=predicted['pulse'];area=endpoint['pulse_area_integral'];clock_gap=locked['r0']*area
    measured.update(Q_lower_bound=locked['Q0']-settings['amplitude'],
        Q_endpoint=locked['Q0'],J_endpoint=0.,Jdot_endpoint=0.,integral_J_numeric=area,integral_J_analytic=0.,
        proper_clock_endpoint=locked['N0']*settings['T']+clock_gap,reference_proper_clock_endpoint=locked['N0']*settings['T'],
        proper_clock_endpoint_gap=clock_gap,proper_clock_gap_computed_from_integrated_pulse_area=True)
    _input_check(locked['input_hashes'])
    if pins!=producer_hashes():raise ValueError('pulse producer changed during evolution')
    record={'schema':SCHEMA,'mode':'nonlinear_pulse_measurement','producers':pins,'producing_commit':locked['producing_commit'],
            'input_hashes':locked['input_hashes'],'prediction_json_sha256':sha256(Path(directory)/'prediction.json'),
            'measurement':measured,'solver':{key:value for key,value in endpoint.items() if not isinstance(value,np.ndarray)},
            'aggregate_CPU_seconds':locked['CPU_seconds']+time.process_time()-started,'CPU_limit_seconds':CPU_LIMIT,
            'source_full_finite_CAR_interpretation':True,'physical_filled_sea_Gamma_matching':'unresolved',
            'Gaussian_reference_gap':locked['Gaussian_reference_gap'],'no_dense_propagator_or_covariance_history':True}
    payload={key:value for key,value in endpoint.items() if isinstance(value,np.ndarray)}
    payload['pulse_area_integral_saved']=np.asarray(area)
    return _write(directory,'measurement',record,payload)


def check(directory):
    prepared,source=_read(directory,'prepare')
    if not (Path(directory)/'prediction.json').exists():
        return {'schema':SCHEMA,'mode':'read_only_check','ok':True,'evolved':False,'source_operator_preparation_authenticated':True}
    locked,arrays=_read(directory,'prediction')
    if locked['prepare_json_sha256']!=sha256(Path(directory)/'prepare.json'):raise ValueError('prediction preparation binding failed')
    predicted=forecast(dict(arrays,Q0=locked['Q0'],M=locked['M'],N0=locked['N0']),locked['prediction']['pulse'])
    for key in ('second_order_work_full_C','second_order_work_passive_reference','second_order_work_vacuum_reference'):
        if abs(predicted[key]-locked['prediction'][key])>1e-14:raise ValueError('locked work forecast replay failed')
    report={'schema':SCHEMA,'mode':'read_only_check','ok':True,'evolved':False,'original_prediction_preserved':True}
    if (Path(directory)/'measurement.json').exists():
        measured,final=_read(directory,'measurement')
        if measured['prediction_json_sha256']!=sha256(Path(directory)/'prediction.json'):raise ValueError('measurement prediction binding failed')
        U=final['U_interaction_endpoint']
        final['unitarity_max']=float(np.max(np.abs(np.einsum('nai,naj->nij',U.conj(),U)-np.eye(2))))
        if max(np.max(np.abs(final['transition_down']-np.abs(U[:,0,1])**2)),np.max(np.abs(final['transition_up']-np.abs(U[:,1,0])**2)))>1e-14:
            raise ValueError('stored transition probabilities disagree with local block endpoints')
        replay=measurement(dict(arrays,M=locked['M']),final)
        for key in ('work_full_C','work_passive_reference','energy_work_closure_full','energy_work_closure_reference'):
            if abs(replay[key]-measured['measurement'][key])>1e-14:raise ValueError('pulse endpoint work replay failed')
        area=float(final['pulse_area_integral_saved'])
        if abs(area-measured['measurement']['integral_J_numeric'])>1e-14 or abs(locked['r0']*area-measured['measurement']['proper_clock_endpoint_gap'])>1e-14:
            raise ValueError('pulse proper-clock integral replay failed')
        report['measurement_replayed']=True;report['measurement']=replay
    return report
