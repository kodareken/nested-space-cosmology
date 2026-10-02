"""Finite inherited grandchild forecast and driven nested elimination.

The only dynamics are the owned common-action rates, full analytic Jv and RK4.
Observer bases are fixed at the saved opening. Spatial content is distinct from
modal covariance. Forecast creation precedes all new held-out-arm evolution.
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

from . import nsc_discovery_prediction as prediction
from . import nsc_discovery_response as response
from . import nsc_discovery_backend as backend
from . import nsc_discovery_episode as episode
from . import nsc_nested_parent_child as model
from . import nsc_spherical_galerkin_coupling as galerkin
from . import nsc_spherical_coupling as coupling

LAB = Path(__file__).resolve().parents[2]
INPUT = LAB / 'results/development/nsc-discovery-prediction-v1'
OUTPUT = LAB / 'results/development/nsc-discovery-grandchild-v1'
SCHEMA = 'NSC-DISCOVERY-GRANDCHILD-v1'
DELTA_TAU = .05
MAX_PROPER_INCREMENT = 1.
ALPHA = .01
CPU_BUDGET = 300.
CHUNK_CAP = 64 * 1024 * 1024
IG = (1.5, 2.5)
IC = (1., 3.)
IP = (0., 4.)
SPACE_PIECES = (('left',(0.,1.)), ('child_left',(1.,1.5)), ('grandchild',IG),
                ('child_right',(2.5,3.)), ('right',(3.,4.)), ('ambient',(4.,8.)))
GEOMETRY = model.GEOMETRY_NAMES + model.MOMENTUM_NAMES


def producers():
    modules = (prediction, response, backend, episode, model, galerkin, coupling)
    paths = [Path(module.__file__) for module in modules] + [Path(__file__), LAB/'scripts/derive_nsc_discovery_grandchild.py',
             LAB/'src/recursive_horizons/nsc_spherical_feedback_action.py', LAB/'src/recursive_horizons/nsc_conformal_adm_source.py']
    return {str(path.relative_to(LAB)): episode.file_sha256(path) for path in paths}


def _git_blob_sha(commit,path):
    """Read a frozen local Git blob; no checkout or repository mutation."""
    if not re.fullmatch(r'[0-9a-f]{40}',str(commit)):
        raise ValueError('historical producing commit must be a full hexadecimal SHA')
    relative=Path(path)
    if relative.is_absolute() or '..' in relative.parts:
        raise ValueError('producer path must stay relative to lab')
    result=subprocess.run(['git','cat-file','blob',f'{commit}:lab/{relative.as_posix()}'],
                          cwd=LAB.parent,capture_output=True,check=True,timeout=5)
    return hashlib.sha256(result.stdout).hexdigest()


def _current_producing_commit(pins):
    """Only call HEAD the producer when every recorded file matches its blob."""
    try:
        result=subprocess.run(['git','rev-parse','HEAD'],cwd=LAB.parent,capture_output=True,check=True,timeout=5)
        commit=result.stdout.decode().strip()
        if all(_git_blob_sha(commit,path)==digest for path,digest in pins.items()):return commit
    except (OSError,subprocess.SubprocessError,ValueError):
        pass
    return None


def _authenticate_historical_producers(directory,prefix,record):
    binding_path=Path(directory)/'observed-run-binding.json'
    if not binding_path.is_file():raise ValueError('historical producer mismatch has no observed run binding')
    binding=json.loads(binding_path.read_text())
    if binding.get('schema')!='NSC-DISCOVERY-OBSERVED-RUN-BINDING-v1':raise ValueError('unsupported historical run binding')
    for extension in ('.json','.npz'):
        path=Path(directory)/(prefix+extension)
        key=str(path.resolve().relative_to(LAB.parent))
        entry=binding.get('artifact_paths',{}).get(key)
        if not entry or episode.file_sha256(path)!=entry['sha256'] or path.stat().st_size!=entry['bytes']:
            raise ValueError('historical grandchild artifact differs from its observed run binding')
    commit=binding['producing_commit']
    try:
        for path,digest in record['producers'].items():
            if _git_blob_sha(commit,path)!=digest:raise ValueError('historical producer differs from its frozen Git blob: '+path)
    except (OSError,subprocess.SubprocessError) as error:
        raise ValueError('historical producing Git blobs are unavailable') from error
    return commit


def validate_proper_increment(value):
    value=float(value)
    if not np.isfinite(value) or not 0.<value<=MAX_PROPER_INCREMENT:
        raise ValueError('finite future proper increment must lie in (0,1]')
    return value


def _guard(started, used=0.):
    if used + time.process_time() - started >= CPU_BUDGET:
        raise RuntimeError('grandchild aggregate 300 CPU-second budget exhausted')


def columns(state):
    return np.vstack((state.phi0,state.phi1))


def with_columns(state, field):
    result = state.copy()
    nf = len(state.phi0)
    result.phi0, result.phi1 = field[:nf].copy(),field[nf:].copy()
    return result


def content(pair,state,tangent=None,interval=IG):
    nodal = model.reconstruct_state(pair,state)
    fine = galerkin.prolong_state(pair.grid,nodal)
    weights = np.asarray(pair.weights)
    mass = np.sum(weights[None,:]*(np.abs(fine.phi0)**2+np.abs(fine.phi1)**2),axis=1)/pair.grid.dx_q
    value = model.interval_integral(pair.grid,mass,interval)
    if tangent is None:
        return value
    lifted = response._prolong_tangent(pair.grid,prediction._nodal_tangent(pair,tangent))
    delta = response._density_tangent(fine.phi0,fine.phi1,lifted.phi0,lifted.phi1,weights,lifted.occupations,pair.grid.dx_q)
    return value,model.interval_integral(pair.grid,delta,interval)


def content_velocity(pair,state):
    rate = model.rates(pair,state)
    velocity = response.StateTangent(*(getattr(rate,name) for name in model.STATE_NAMES),np.zeros(6))
    return content(pair,state,velocity)[1]


def recover_opening_clock(pair,state,tangent,case):
    raw = response.nested_readout(pair,state,tangent,coordinate_duration=0.)
    rows = []
    for key,delta,dot in [('primary','delta_child_regional_content_t','child_regional_content_dot'),
                          ('secondary','delta_child_proper_mean_r_t','child_proper_mean_r_dot')]:
        slope = float(raw[dot])
        if abs(slope) < 1e-10:
            raise RuntimeError('saved analytic clock recovery is ill-conditioned; unchanged-prefix replay required')
        value = (raw[delta]-case['analytic_'+key])*raw['tau_dot']/slope
        rows.append(float(value))
    gap = abs(rows[0]-rows[1])
    if gap > 1e-8 * max(1.,abs(rows[0]),abs(rows[1])):
        raise RuntimeError('primary/secondary opening clock tangents disagree; unchanged-prefix replay required')
    rate = model.rates(pair,state)
    aligned = prediction._shift_tangent(tangent,rate,-rows[0]/raw['tau_dot'])
    return aligned, {'delta_tau_primary':rows[0], 'delta_tau_secondary':rows[1], 'cross_check_gap':gap,
                     'tau_dot':raw['tau_dot'], 'formula':'delta_y_tau = delta_y_t - F delta_tau / tau_dot',
                     'new_segment_delta_tau':0., 'occupation_direction_preserved':np.array_equal(aligned.occupations,tangent.occupations)}


def load_opening(nf=256):
    if nf not in (128,256):
        raise ValueError('saved grandchild openings are nf128 or nf256')
    paths = (INPUT/'prediction-run.json',INPUT/'prediction-run.npz',episode.V1_JSON,episode.V1_NPZ,episode.BASIS_JSON,episode.BASIS_NPZ)
    hashes = {str(path.relative_to(LAB)):episode.file_sha256(path) for path in paths}
    record = json.loads(paths[0].read_text())
    if hashes[str(paths[1].relative_to(LAB))] != record['payload_sha256']:
        raise ValueError('prediction payload binding failed')
    case = next(case for case in record['cases'] if case['nf']==nf)
    if abs(case['duration']-.3)>1e-12:
        raise ValueError('inherited opening must be the saved T=.3 prediction endpoint')
    arm = next(arm for arm in case['arms'] if arm['alpha']==ALPHA)
    if not arm['equal_tau']['reached'] or arm['equal_tau']['uses_linear_clock_correction']:
        raise ValueError('held-out opening must be independently rooted at equal proper time')
    with np.load(paths[1],allow_pickle=False) as payload:
        state = prediction._state_from_prefix(payload,case['analytic_state_prefix']).copy()
        held = prediction._state_from_prefix(payload,arm['event_prefix']).copy()
        prefix = case['analytic_tangent_prefix']
        tangent = response.StateTangent(*(np.array(payload[prefix+'_'+name],copy=True) for name in (*model.STATE_NAMES,'occupations')))
        saved_W = np.array(payload[case['W_key']],copy=True)
    seed = prediction.load_saved_baseline(case['case_name'],time=0.)
    if not np.array_equal(saved_W,seed['pair'].geometry_map):
        raise ValueError('prediction W differs from the frozen inherited basis')
    pair = backend.make_fft_pair(seed['pair'])
    aligned,clock = recover_opening_clock(pair,state,tangent,case)
    held_pair = prediction.pair_with_occupations(pair,arm['weights'])
    prediction.assert_admissible_source(held.phi0,held.phi1,held_pair.weights)
    return {'pair':pair,'held_pair':held_pair,'state':state,'held':held,'tangent':aligned,
            'clock_alignment':clock,'tau_opening':float(case['tau_target']),
            'held_tau_opening':float(arm['equal_tau']['tau']),'held_coordinate_time':float(arm['equal_tau']['time']),
            'input_hashes':hashes,'nf':nf}


def _ap_values(values,coordinates,period=8.):
    count = values.shape[0]
    demod = np.exp(-1j*np.pi*np.arange(count)/count)[:,None]
    coefficients = np.fft.fft(values*demod,axis=0)/count
    modes = np.fft.fftfreq(count)*count+.5
    return np.exp(2j*np.pi*np.asarray(coordinates)[:,None]*modes[None,:]/period) @ coefficients


def observer_basis(pair):
    """One QR: dilated grandchild seeds, middle pair, original four outer modes."""
    grid = pair.grid
    coordinate = grid.xi_q
    mapped = 2.+2.*(coordinate-2.)
    mask = ((coordinate>=IG[0]) & (coordinate<IG[1]))[:,None]
    middle = pair.original_columns[:,[2,3]]
    high = []
    low = []
    for block in (middle[:grid.nf],middle[grid.nf:]):
        samples = np.sqrt(2.*grid.nf/grid.nq)*_ap_values(block,mapped)*mask
        projected = grid.U_f.conj().T @ samples
        high.append(samples); low.append(projected)
    seeds = np.vstack(low)
    fine_seeds = np.vstack(high)
    prolonged = np.vstack((grid.U_f@low[0],grid.U_f@low[1]))
    projection_tail = float(np.linalg.norm(fine_seeds-prolonged)/np.linalg.norm(fine_seeds))
    packed = np.concatenate((seeds,middle,pair.original_columns[:,[0,1,4,5]]),axis=1)
    basis,R = np.linalg.qr(packed,mode='reduced')
    if np.min(np.abs(np.diag(R))) < 1e-10:
        raise ValueError('grandchild/child/parent seeds do not have ranks 2/4/8')
    basis.setflags(write=False)
    tails = {}
    for name,width,interval in [('grandchild',2,IG),('child',4,IC),('parent',8,IP)]:
        fine0 = grid.U_f @ basis[:grid.nf,:width]
        fine1 = grid.U_f @ basis[grid.nf:,:width]
        density = np.abs(fine0)**2+np.abs(fine1)**2
        inside = [model.interval_integral(grid,density[:,column]/grid.dx_q,interval) for column in range(width)]
        tails[name] = {'inside_mass':inside,'outside_mass':[1.-value for value in inside], 'exact_compact_support_claimed':False}
    return basis, {'ranks':[2,4,8], 'single_opening_QR':True,'dilation_factor':2.,
                   'dilated_seed_projection_tail':projection_tail,'spatial_tails':tails,
                   'gram_gap':float(np.max(np.abs(basis.conj().T@basis-np.eye(8)))),
                   'primary_is_spatial_content_not_projector_trace':True}


def initial_routes(field,basis):
    G,D,F = basis[:,:2],basis[:,2:4],basis[:,4:8]
    direct = {'a':G.conj().T@field,'drive':field-G@(G.conj().T@field),'memory':np.zeros_like(field)}
    outer = field-basis@(basis.conj().T@field)
    sequential = {'a':G.conj().T@field,'d_drive':D.conj().T@field,'d_memory':np.zeros((2,field.shape[1]),complex),
                  'f_drive':F.conj().T@field,'f_memory':np.zeros((4,field.shape[1]),complex),
                  'drive':outer,'memory':np.zeros_like(field)}
    return direct,sequential


def reconstruct_route(route,basis,sequential=False):
    result = basis[:,:2]@route['a']+route['drive']+route['memory']
    if sequential:
        result = result + basis[:,2:4]@(route['d_drive']+route['d_memory'])+basis[:,4:8]@(route['f_drive']+route['f_memory'])
    return result


def route_rates(image,basis,direct,sequential,*,omit_outer_detail=False):
    """Streaming causal splits. Detail memory includes the outer reconstructed drive."""
    G,D,F = basis[:,:2],basis[:,2:4],basis[:,4:8]
    parts = [G@direct['a'],direct['drive'],direct['memory'],G@sequential['a'],
             D@sequential['d_drive'],D@sequential['d_memory'],F@sequential['f_drive'],F@sequential['f_memory'],
             sequential['drive'],sequential['memory']]
    images = np.split(image(np.concatenate(parts,axis=1)),10,axis=1)
    dr,dd,dm,g,ddrive,dmem,fdrive,fmem,odrive,omem = images
    qG = lambda value:value-G@(G.conj().T@value)
    qP = lambda value:value-basis@(basis.conj().T@value)
    first = {'a':-1j*G.conj().T@(dr+dd+dm),'drive':-1j*qG(dd),'memory':-1j*qG(dm+dr)}
    outside = fdrive+fmem+odrive+omem
    second = {'a':-1j*G.conj().T@(g+ddrive+dmem+outside),
              'd_drive':-1j*D.conj().T@ddrive,
              'd_memory':-1j*D.conj().T@(dmem+g+(0. if omit_outer_detail else outside)),
              'f_drive':-1j*F.conj().T@fdrive,
              'f_memory':-1j*F.conj().T@(fmem+g+ddrive+dmem+odrive+omem),
              'drive':-1j*qP(odrive),'memory':-1j*qP(omem+g+ddrive+dmem+fdrive+fmem)}
    return first,second


def _combine_route(route,rate,factor):
    return {key:value+factor*rate[key] for key,value in route.items()}


def routed_step(pair,state,direct,sequential,basis,dt):
    """Same full RK4 geometry; reconstruct both sources at every stage before forces."""
    states=[state]; directs=[direct]; sequences=[sequential]
    rates=[]; drates=[]; srates=[]; clocks=[]; source_gaps=[]
    zero=response.zero_tangent(pair.grid)
    for index,factor in enumerate((.5,.5,1.,None)):
        current,dr,sr=states[index],directs[index],sequences[index]
        dfield=reconstruct_route(dr,basis); sfield=reconstruct_route(sr,basis,True)
        dfull,dbundle=model.rates(pair,with_columns(current,dfield),return_bundle=True)
        sfull,sbundle=model.rates(pair,with_columns(current,sfield),return_bundle=True)
        full,bundle=model.rates(pair,current,return_bundle=True)
        source_gaps.append(max(float(np.max(np.abs(dbundle['source'][key]-bundle['source'][key]))) for key in ('force_Q','force_L','force_beta')))
        source_gaps.append(max(float(np.max(np.abs(sbundle['source'][key]-bundle['source'][key]))) for key in ('force_Q','force_L','force_beta')))
        dnext,snext=route_rates(lambda field:model.apply_hamiltonian(pair,current,field),basis,dr,sr)
        rates.append(full); drates.append(dnext); srates.append(snext)
        clocks.append(prediction.clock_rates(pair,current,zero)[0])
        if factor is not None:
            states.append(model._combine(state,full,dt*factor))
            directs.append(_combine_route(direct,dnext,dt*factor)); sequences.append(_combine_route(sequential,snext,dt*factor))
    out=model.NestedState(*(getattr(state,name)+dt/6.*(getattr(rates[0],name)+2*getattr(rates[1],name)+2*getattr(rates[2],name)+getattr(rates[3],name)) for name in model.STATE_NAMES))
    average=lambda items:{key:(items[0][key]+2*items[1][key]+2*items[2][key]+items[3][key])/6. for key in items[0]}
    dr=_combine_route(direct,average(drates),dt); sr=_combine_route(sequential,average(srates),dt)
    prediction._chart(pair,out)
    return {'state':out,'direct':dr,'sequential':sr,'tau':dt*prediction._rk_quadrature(clocks),
            'stage_source_gap':max(source_gaps), 'reconstruction_gap':max(float(np.max(np.abs(reconstruct_route(dr,basis)-columns(out)))),float(np.max(np.abs(reconstruct_route(sr,basis,True)-columns(out)))))}


def march_analytic(pair,state,tangent,target,cap,started,used=0.):
    current=state.copy(); direction=tangent; tau=0.; delta_tau=0.; mark=0.; steps=0
    while tau<target:
        _guard(started,used)
        dt=model.stable_timestep(pair,current,cap)[0]
        old=current.copy(); old_direction=direction; old_tau=tau; old_delta=delta_tau; old_mark=mark
        advanced=prediction.coupled_rk4_step(pair,current,direction,dt)
        current,direction=advanced['state'],advanced['tangent']; tau+=advanced['tau'];delta_tau+=advanced['delta_tau'];mark+=dt;steps+=1
        if tau>=target:
            _guard(started,used)
            event=prediction._root_step_to_tau(pair,{'state':old,'tau':old_tau,'time':old_mark,'dt_hi':dt},target)
            rooted=prediction.coupled_rk4_step(pair,old,old_direction,event['dt'])
            current,direction=rooted['state'],rooted['tangent'];tau=old_tau+rooted['tau'];delta_tau=old_delta+rooted['delta_tau'];mark=event['time']
            break
    value,raw=content(pair,current,direction)
    dot=content_velocity(pair,current); tau_dot=prediction.clock_rates(pair,current,direction)[0]
    derivative=raw-dot*delta_tau/tau_dot
    return current,direction,{'baseline_NG':value,'derivative_NG_tau':derivative,'raw_derivative_NG_t':raw,'NG_dot':dot,
        'delta_tau':delta_tau,'tau_dot':tau_dot,'proper_increment':tau,'tau_residual':tau-target,
        'coordinate_duration':mark,'steps':steps,'root_bisections':event['bisections'],
        'predicted_heldout_NG':value+ALPHA*derivative,'predicted_change':ALPHA*derivative,
        'forecast_frozen_before_new_heldout_evolution':True}


def march_held(pair,state,basis,target,cap,started,used):
    current=state.copy(); direct,sequential=initial_routes(columns(current),basis)
    tau=0.;mark=0.;steps=0;max_force=0.;max_rebuild=0.
    while tau<target:
        _guard(started,used)
        dt=model.stable_timestep(pair,current,cap)[0]
        old=current.copy();old_d=direct;old_s=sequential;old_tau=tau;old_mark=mark
        advanced=routed_step(pair,current,direct,sequential,basis,dt)
        if tau+advanced['tau']>=target:
            _guard(started,used)
            event=prediction._root_step_to_tau(pair,{'state':old,'tau':old_tau,'time':old_mark,'dt_hi':dt},target)
            dt=event['dt'];advanced=routed_step(pair,old,old_d,old_s,basis,dt)
        current,direct,sequential=advanced['state'],advanced['direct'],advanced['sequential']
        tau+=advanced['tau'];mark+=dt;steps+=1
        max_force=max(max_force,advanced['stage_source_gap']);max_rebuild=max(max_rebuild,advanced['reconstruction_gap'])
        if 'event' in locals():break
    return current,direct,sequential,{'proper_increment':tau,'tau_residual':tau-target,'coordinate_duration':mark,
           'steps':steps,'max_stage_force_gap':max_force,'max_reconstruction_gap':max_rebuild,
           'root_bisections':event['bisections'],'uses_linear_clock_correction':False}


def endpoint_diagnostics(pair,state,basis,direct=None,sequential=None):
    field=columns(state);G=basis[:,:2]
    amplitudes=G.conj().T@field
    covariance=(amplitudes*pair.weights[None,:])@amplitudes.conj().T
    exterior=field-G@amplitudes
    cross_norm=float(np.sqrt(max(0.,np.trace((amplitudes*pair.weights)@(exterior.conj().T@exterior)@(amplitudes*pair.weights).conj().T).real)))
    fine,system=model._active(pair,state);source=coupling.source_from_columns(system,fine)
    energy={name:model.interval_integral(pair.grid,source['force_L']/fine.r/pair.grid.dx_q,interval) for name,interval in SPACE_PIECES}
    gram=field.conj().T@field;w=np.sqrt(pair.weights)
    eig=np.linalg.eigvalsh(w[:,None]*gram*w[None,:])
    report={'NG':content(pair,state),'grand_covariance_real':covariance.real.tolist(),'grand_covariance_imag':covariance.imag.tolist(),
            'modal_trace_is_secondary':float(np.trace(covariance).real),'grand_exterior_cross_frobenius':cross_norm,
            'source_gram_gap':float(np.max(np.abs(gram-np.eye(6)))),'source_CAR_eigenvalues':eig.tolist(),
            'source_CAR_admissible':bool(np.min(eig)>=-1e-8 and np.max(eig)<=1.+1e-8),'disjoint_normal_energies':energy,
            'normal_energy_closure_claimed':False,'clock_rate':prediction.clock_rates(pair,state,response.zero_tangent(pair.grid))[0]}
    # A diagonal source partition alone drops physical cross terms. Compare
    # its forces against those from the reconstructed full six-column source.
    pieces=[basis[:,:2]@(basis[:,:2].conj().T@field),basis[:,2:4]@(basis[:,2:4].conj().T@field),
            basis[:,4:8]@(basis[:,4:8].conj().T@field),field-basis@(basis.conj().T@field)]
    diagonal={key:np.zeros_like(source[key]) for key in ('force_Q','force_L','force_beta')}
    for piece in pieces:
        piece_fine,piece_system=model._active(pair,with_columns(state,piece))
        piece_source=coupling.source_from_columns(piece_system,piece_fine)
        for key in diagonal:diagonal[key]+=piece_source[key]
    report['cross_force_max']={key:float(np.max(np.abs(source[key]-diagonal[key]))) for key in diagonal}
    report['forces_use_full_reconstructed_source_with_cross']=True
    if direct is not None:
        dfield=reconstruct_route(direct,basis);sfield=reconstruct_route(sequential,basis,True)
        da=G.conj().T@dfield;sa=G.conj().T@sfield
        dcov=(da*pair.weights)@da.conj().T;scov=(sa*pair.weights)@sa.conj().T
        report.update(direct_sequential_Phi_gap=float(np.max(np.abs(dfield-sfield))),grand_covariance_composition_gap=float(np.max(np.abs(dcov-scov))),
                      direct_NG=content(pair,with_columns(state,dfield)),sequential_NG=content(pair,with_columns(state,sfield)),
                      direct_clock_rate=prediction.clock_rates(pair,with_columns(state,dfield),response.zero_tangent(pair.grid))[0],
                      sequential_clock_rate=prediction.clock_rates(pair,with_columns(state,sfield),response.zero_tangent(pair.grid))[0],
                      geometry_shared_between_routes=True,outer_detail_drive_included=True)
    return report


def _pack_state(arrays,prefix,state):
    prediction._pack_state(arrays,prefix,state)


def _pack_direction(arrays,prefix,tangent):
    prediction._pack_tangent(arrays,prefix,tangent)


def _check_bindings(hashes):
    if any(episode.file_sha256(LAB/path)!=digest for path,digest in hashes.items()):
        raise ValueError('saved input hashes changed')


def _write(directory,prefix,record,arrays):
    directory=Path(directory)
    jpath=directory/(prefix+'.json');npath=directory/(prefix+'.npz')
    if jpath.exists() or npath.exists():raise FileExistsError('immutable grandchild prefix already exists: '+prefix)
    stream=io.BytesIO();np.savez_compressed(stream,**arrays);blob=stream.getvalue()
    if len(blob)>CHUNK_CAP:raise RuntimeError('grandchild checkpoint exceeds 64 MiB')
    directory.mkdir(parents=True,exist_ok=True)
    record=dict(record,payload_sha256=hashlib.sha256(blob).hexdigest(),array_sha256={key:episode.array_sha256(value) for key,value in arrays.items()},payload_bytes=len(blob))
    json_blob=(json.dumps(record,indent=2,sort_keys=True,allow_nan=False)+'\n').encode()
    if len(blob)+len(json_blob)>CHUNK_CAP:raise RuntimeError('grandchild record exceeds 64 MiB')
    for path,data in ((npath,blob),(jpath,json_blob)):
        descriptor=os.open(path,os.O_CREAT|os.O_EXCL|os.O_WRONLY,0o444)
        with os.fdopen(descriptor,'wb') as handle:handle.write(data)
    return record


def _read(directory,prefix,*,allow_historical=False):
    directory=Path(directory);record=json.loads((directory/(prefix+'.json')).read_text())
    npath=directory/(prefix+'.npz')
    if episode.file_sha256(npath)!=record['payload_sha256']:raise ValueError('grandchild payload binding failed')
    with np.load(npath,allow_pickle=False) as stored:arrays={key:stored[key].copy() for key in stored.files}
    if any(episode.array_sha256(arrays[key])!=digest for key,digest in record['array_sha256'].items()):raise ValueError('grandchild array binding failed')
    _check_bindings(record['input_hashes'])
    if record['producers']!=producers():
        if not allow_historical:raise ValueError('grandchild producer sources differ from the frozen record')
        record=dict(record,historical_producing_commit=_authenticate_historical_producers(directory,prefix,record),
                    producer_binding_method='observed immutable artifacts and frozen local Git blobs')
    return record,arrays


def preview(nf=256,*,proper_increment=DELTA_TAU):
    proper_increment=validate_proper_increment(proper_increment)
    opening=load_opening(nf);basis,basis_report=observer_basis(opening['pair'])
    return {'schema':SCHEMA,'mode':'preview','nf':nf,'opening_clock_alignment':opening['clock_alignment'],'basis':basis_report,
            'opening':endpoint_diagnostics(opening['pair'],opening['state'],basis),'heldout_opening':endpoint_diagnostics(opening['held_pair'],opening['held'],basis),
            'future_proper_increment':proper_increment,'absolute_future_tau_target':opening['tau_opening']+proper_increment,'evolved':False,
            'physical_claim_from_composition_identity':False,'CPU_budget_seconds':CPU_BUDGET}


def prepare(directory,nf=256,*,proper_increment=DELTA_TAU,matched=True):
    directory=Path(directory)
    if (directory/'forecast.json').exists() or (directory/'forecast.npz').exists():raise FileExistsError('forecast prefix already exists')
    proper_increment=validate_proper_increment(proper_increment)
    started=time.process_time();pins=producers();opening=load_opening(nf)
    producing_commit=_current_producing_commit(pins)
    pair,state,tangent=opening['pair'],opening['state'],opening['tangent']
    basis,basis_report=observer_basis(pair)
    cap=.0005
    t=time.process_time();prediction.coupled_rk4_step(pair,state,tangent,cap);analytic_cpu=time.process_time()-t
    # Admission probes use the baseline only. No future held-out state is
    # generated before the forecast has been written and locked.
    dr,sr=initial_routes(columns(state),basis)
    t=time.process_time();routed_step(pair,state,dr,sr,basis,cap);routed_cpu=time.process_time()-t
    tau_dot=prediction.clock_rates(pair,state,tangent)[0]
    estimate_steps=int(np.ceil(proper_increment/max(tau_dot,1e-8)/cap))
    forecast_one=1.75*((estimate_steps+22)*analytic_cpu+(estimate_steps+22)*routed_cpu)
    remaining=CPU_BUDGET-(time.process_time()-started)
    if forecast_one>remaining:raise RuntimeError('grandchild continuation not admitted under 300 CPU seconds')
    arrays={'basis':basis}
    _pack_state(arrays,'opening_analytic',state);_pack_direction(arrays,'opening_aligned_tangent',tangent);_pack_state(arrays,'opening_heldout',opening['held'])
    forecasts=[]
    heldout_costs={}
    caps=[cap]
    matched_status='not_requested' if not matched else 'deferred_budget_requires_later_confirmation'
    for index,step_cap in enumerate(caps):
        future,direction,forecast=march_analytic(pair,state,tangent,proper_increment,step_cap,started)
        _pack_state(arrays,f'future_{index}',future);_pack_direction(arrays,f'future_tangent_{index}',direction)
        forecasts.append(dict(forecast,step_cap=step_cap,prefix=f'future_{index}'))
        heldout_costs[f'{step_cap:g}']=1.5*(forecast['steps']+22)*routed_cpu
        elapsed=time.process_time()-started
        if elapsed+sum(heldout_costs.values())>CPU_BUDGET:
            raise RuntimeError('actual baseline forecast leaves insufficient CPU for the primary heldout measurement')
        if index==0 and matched:
            fine_steps=2*forecast['steps']+22
            fine_reserve=1.5*fine_steps*(analytic_cpu+routed_cpu)
            if elapsed+sum(heldout_costs.values())+fine_reserve<CPU_BUDGET:
                caps.append(cap/2.);matched_status='admitted_after_actual_primary_forecast'
    _guard(started);_check_bindings(opening['input_hashes'])
    if pins!=producers():raise ValueError('producer sources changed while preparing the forecast')
    record={'schema':SCHEMA,'mode':'locked_forecast','nf':nf,'alpha':ALPHA,'IG':list(IG),'IC':list(IC),'IP':list(IP),
        'basis':basis_report,'opening_clock_alignment':opening['clock_alignment'],'tau_opening':opening['tau_opening'],
        'held_tau_opening':opening['held_tau_opening'],'absolute_future_tau_target':opening['tau_opening']+proper_increment,
        'proper_increment':proper_increment,'forecasts':forecasts,'input_hashes':opening['input_hashes'],'producers':pins,
        'producing_commit':producing_commit,'producing_commit_matches_declared_files':producing_commit is not None,
        'forecast_admission':{'probe_analytic_step_CPU':analytic_cpu,'probe_routed_step_CPU':routed_cpu,'primary_full_run_forecast_CPU':forecast_one,
                             'probes_use_baseline_only':True,
                             'heldout_forecast_CPU_by_cap':heldout_costs,
                             'matched_step_status':matched_status,
                             'later_step_confirmation_required':matched and len(caps)==1,
                             'actual_baseline_steps':[item['steps'] for item in forecasts],
                             'matched_step_admitted':len(caps)==2,'CPU_budget_seconds':CPU_BUDGET,'checkpoint_cap_bytes':CHUNK_CAP},
        'CPU_used_before_measurement':time.process_time()-started,'new_heldout_future_measured':False,
        'forecast_locked_before_measurement':True,'initializer_called':False,'source_reset':False,'geometry_reset':False,
        'physical_claim_from_composition_identity':False,'continuous_error_bound':None,
        'opening':endpoint_diagnostics(pair,state,basis)}
    return _write(directory,'forecast',record,arrays)


def run(directory):
    directory=Path(directory)
    if (directory/'measurement.json').exists() or (directory/'measurement.npz').exists():raise FileExistsError('measurement prefix already exists')
    started=time.process_time();forecast,arrays=_read(directory,'forecast');used=forecast['CPU_used_before_measurement']
    check(directory)  # Replay the locked numerical forecast before new held-out evolution.
    opening=load_opening(forecast['nf']);basis=arrays['basis'];pair=opening['held_pair'];state=prediction._state_from_prefix(arrays,'opening_heldout')
    target=forecast['absolute_future_tau_target']-forecast['held_tau_opening']
    results=[];endpoints={}
    for index,locked in enumerate(forecast['forecasts']):
        estimated=forecast['forecast_admission'].get('heldout_forecast_CPU_by_cap',{}).get(f"{locked['step_cap']:g}",forecast['forecast_admission']['primary_full_run_forecast_CPU']*(2 if index else 1))
        if index and used+time.process_time()-started+estimated>CPU_BUDGET:
            results.append({'step_cap':locked['step_cap'],'status':'matched_measurement_deferred_by_budget'});break
        final,direct,sequential,event=march_held(pair,state,basis,target,locked['step_cap'],started,used)
        measured=endpoint_diagnostics(pair,final,basis,direct,sequential)
        measured_change=measured['NG']-locked['baseline_NG']
        error=measured['NG']-locked['predicted_heldout_NG']
        results.append({'step_cap':locked['step_cap'],'status':'heldout_future_measured','event':event,'measurement':measured,
                        'baseline_NG':locked['baseline_NG'],'predicted_NG':locked['predicted_heldout_NG'],
                        'predicted_change':locked['predicted_change'],'measured_change':measured_change,'prediction_error':error,
                        'error_over_measured_effect':None if measured_change==0 else abs(error/measured_change),
                        'absolute_future_tau':forecast['held_tau_opening']+event['proper_increment']})
        _pack_state(endpoints,f'heldout_{index}',final)
        for label,route in [('direct',direct),('sequential',sequential)]:
            for key,value in route.items():endpoints[f'{label}_{index}_{key}']=value
    _guard(started,used);_check_bindings(forecast['input_hashes'])
    if producers()!=forecast['producers']:raise ValueError('producer source changed during heldout measurement')
    comparison=None
    if len(results)==2 and results[1]['status']=='heldout_future_measured':
        comparison={'NG_gap':results[1]['measurement']['NG']-results[0]['measurement']['NG'],
                    'prediction_error_gap':results[1]['prediction_error']-results[0]['prediction_error'],
                    'fine_primary_effect':results[1]['measured_change']}
    record={'schema':SCHEMA,'mode':'heldout_measurement','nf':forecast['nf'],'input_hashes':forecast['input_hashes'],'producers':forecast['producers'],
        'forecast_json_sha256':episode.file_sha256(directory/'forecast.json'),'forecast_payload_sha256':forecast['payload_sha256'],
        'results':results,'matched_step':comparison,'aggregate_CPU_seconds':used+time.process_time()-started,'CPU_budget_seconds':CPU_BUDGET,
        'producing_commit':forecast.get('producing_commit'),
        'later_step_confirmation_required':forecast['forecast_admission'].get('later_step_confirmation_required',False) or any(item['status']=='matched_measurement_deferred_by_budget' for item in results),
        'forecast_preceded_future_measurement':True,'new_future_NG_measurement':True,
        'composition_identity_is_separate_from_prediction':True,'conditional_finite_handoff_only':True,
        'universal_Omega_scaling_claimed':False,'normal_energy_closure_claimed':False,'continuous_error_bound':None,
        'dense_exterior_propagators_stored':False,'dense_ambient_covariance_stored':False,'initializer_called':False}
    return _write(directory,'measurement',record,endpoints)


def check(directory):
    forecast,arrays=_read(directory,'forecast',allow_historical=True)
    opening=load_opening(forecast['nf']);basis=arrays['basis'];G=basis[:,:2]
    if np.max(np.abs(basis.conj().T@basis-np.eye(8)))>1e-10:raise ValueError('frozen nested observer is not orthonormal')
    for index,locked in enumerate(forecast['forecasts']):
        state=prediction._state_from_prefix(arrays,f'future_{index}')
        direction=response.StateTangent(*(arrays[f'future_tangent_{index}_{name}'] for name in (*model.STATE_NAMES,'occupations')))
        value,raw=content(opening['pair'],state,direction)
        tau_dot=prediction.clock_rates(opening['pair'],state,direction)[0]
        derivative=raw-content_velocity(opening['pair'],state)*locked['delta_tau']/tau_dot
        if abs(value-locked['baseline_NG'])>1e-10 or abs(value+ALPHA*derivative-locked['predicted_heldout_NG'])>1e-10:
            raise ValueError('locked forecast readout replay disagrees')
    reconstructed=[]
    path=Path(directory)/'measurement.json'
    if path.exists():
        measured,states=_read(directory,'measurement',allow_historical=True)
        if measured['forecast_json_sha256']!=episode.file_sha256(Path(directory)/'forecast.json'):raise ValueError('measurement forecast binding failed')
        for index,result in enumerate(measured['results']):
            if result['status']!='heldout_future_measured':continue
            state=prediction._state_from_prefix(states,f'heldout_{index}')
            value=content(opening['held_pair'],state)
            if abs(value-result['measurement']['NG'])>1e-10:raise ValueError('saved endpoint NG replay disagrees')
            routes=[]
            for label,seq in [('direct',False),('sequential',True)]:
                keys=('a','drive','memory') if not seq else ('a','d_drive','d_memory','f_drive','f_memory','drive','memory')
                route={key:states[f'{label}_{index}_{key}'] for key in keys}
                routes.append(reconstruct_route(route,basis,seq))
            if max(np.max(np.abs(route-columns(state))) for route in routes)>1e-8:raise ValueError('nested endpoint reconstruction replay failed')
            reconstructed.append(value)
    return {'schema':SCHEMA,'ok':True,'mode':'read_only_replay','future_measurement_replayed':bool(reconstructed),'NG':reconstructed,
            'evolved':False,'producer_and_input_hashes_match':True,'composition_only_is_not_a_physical_claim':True,
            'historical_producing_commit':forecast.get('historical_producing_commit'),
            'producer_binding_method':forecast.get('producer_binding_method','current frozen producer bytes'),
            'replay_uses_current_readout_code':True,'historical_trajectory_recomputed':False}
