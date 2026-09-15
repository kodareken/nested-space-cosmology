#!/usr/bin/env python3
"""Full-energy PG covariance for the first retained massive compact group.

Cold preparation computes a new, hash-pinned input artifact from physical
mode owners. Verification recombines that artifact without old generators.
Other massive groups and the complete retained C1b gate remain distinct.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import inspect
import json
from pathlib import Path
import sys
import tempfile

import numpy as np
from numpy.polynomial.legendre import leggauss

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'scripts'))
from recursive_horizons.nsc_paired_horizon_preparation import PairedHorizonSeedMap
from recursive_horizons.nsc_massive_jost_modes import solve_jost
from recursive_horizons.nsc_pg_massive_modes import MassivePGModeResolution
from recursive_horizons.nsc_pg_packet_modes import project_computed_modes,assemble_signed_window,SIGNED_PACKET_MAP
from recursive_horizons.nsc_pg_fast_packets import FastVacuumPacketProjector
from recursive_horizons.nsc_pg_tail_packets import MassiveTailPackets
from recursive_horizons.nsc_pg_tail_series import EndpointTailSeries
from recursive_horizons.nsc_compact_ctp_neck import _frequency_grid
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes
from recursive_horizons.nsc_pg_state_covariance import equal_time_ctp_blocks
import derive_nsc_pg_threshold_projection as threshold

OUTPUT='results/development/nsc-pg-group13-covariance.json'
PHYSICAL_SOURCES=('src/recursive_horizons/nsc_paired_horizon_preparation.py',
 'src/recursive_horizons/nsc_massive_jost_modes.py','src/recursive_horizons/nsc_pg_massive_modes.py',
 'src/recursive_horizons/nsc_pg_packet_modes.py','src/recursive_horizons/nsc_pg_fast_packets.py',
 'src/recursive_horizons/nsc_pg_high_energy.py','src/recursive_horizons/nsc_pg_tail_packets.py',
 'src/recursive_horizons/nsc_pg_tail_series.py','src/recursive_horizons/nsc_pg_threshold_projection.py',
 'src/recursive_horizons/nsc_complex_horizon_modes.py','src/recursive_horizons/nsc_lorentzian.py',
 'src/recursive_horizons/nsc_pg_lll_preparation.py','src/recursive_horizons/nsc_unruh_state.py',
 'src/recursive_horizons/nsc_gauge_source.py','src/recursive_horizons/nsc_common_time_bulk_split.py',
 'src/recursive_horizons/nsc_transmitting_dirac_domain.py','src/recursive_horizons/nsc_compact_ctp_neck.py',
 'src/recursive_horizons/nsc_mode_resolved_cauchy_state.py')

def sha(path):return hashlib.sha256((ROOT/path).read_bytes()).hexdigest()

def context():
    state=json.loads((ROOT/'results/development/nsc-mode-resolved-cauchy-state.json').read_text())
    restart=json.loads((ROOT/'results/development/nsc-compact-matched-restart.json').read_text())
    channel=next(c for c in state['channels'] if c['index']==13)
    if channel['compact_level']!=1 or channel['angular_eigenvalue']!=0:raise ValueError('retained group13 changed')
    return channel,restart['scattering_provenance']['config'],restart['locked_inputs']

def low_node(job):
    index,E,w,config,cache=job;path=Path(cache)/f'low-{index:04d}.npz'
    if path.exists():return
    p=PairedHorizonSeedMap(config['horizon_rho'],config['surface_gravity'],config['omega'],config['horizon_offset'],config['scattering_tolerance'])
    j=solve_jost(p.background,E,np.pi/2,0.)
    modes=MassivePGModeResolution(p).section(E,np.pi/2,0.,[0.],[3.],jost=j)
    F,_=project_computed_modes(E,np.pi/2,0.,modes.interior[0],modes.exterior[0])
    path.write_bytes(deterministic_npz_bytes({'energy':np.array(E),'weight':np.array(w),'projection':F,
                    'source':modes.source_covariance,'projector':modes.source_projector}))

def _threshold_cache(channel,config):
    signature={'source_hashes':{p:threshold.sha(p) for p in threshold.DEPENDENCIES},
       'node_function_sha256':hashlib.sha256(inspect.getsource(threshold.evaluate_node).encode()).hexdigest(),
       'config':config,'channel':channel,'angular_sign':1}
    tag=hashlib.sha256(json.dumps(signature,sort_keys=True).encode()).hexdigest()[:16]
    return Path(tempfile.gettempdir())/('nsc-threshold-contour-'+tag)

def prepare(workers=2):
    channel,config,locked=context();source_hashes={p:sha(p) for p in PHYSICAL_SOURCES}
    signature={'sources':source_hashes,'config':config,'low_node':hashlib.sha256(inspect.getsource(low_node).encode()).hexdigest()}
    tag=hashlib.sha256(json.dumps(signature,sort_keys=True).encode()).hexdigest()[:16]
    cache=Path(tempfile.gettempdir())/('nsc-group13-input-'+tag);cache.mkdir(exist_ok=True)
    E,w,_=_frequency_grid(np.pi/2,40.,32);mask=(E<1.)|((E>np.pi/2)&(E<16.));E=E[mask];w=w[mask]
    print(f'Canonical group13 low input: {len(E)} nodes; {cache}',flush=True)
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for i,_ in enumerate(pool.map(low_node,((i,e,v,config,str(cache)) for i,(e,v) in enumerate(zip(E,w))),chunksize=1),1):
            if i%32==0:print(f'Low input {i}/{len(E)}',flush=True)
    low=[]
    for i in range(len(E)):
        with np.load(cache/f'low-{i:04d}.npz',allow_pickle=False) as a:low.append({k:a[k].copy() for k in a.files})
    arrays={'low_energy':E,'low_weight':w,'low_projection':np.array([a['projection'] for a in low]),
            'low_source':np.array([a['source'] for a in low]),'low_projector':np.array([a['projector'] for a in low])}
    p=PairedHorizonSeedMap(config['horizon_rho'],config['surface_gravity'],config['omega'],config['horizon_offset'],config['scattering_tolerance'])
    fast=FastVacuumPacketProjector(np.pi/2,0.,p,minimum_energy=16.)
    for points in (16,24):
        x,q=leggauss(points);energy=np.concatenate([a+4*(x+1) for a in range(16,160,8)])
        weight=np.tile(4*q,18);file=cache/f'fast-{points}.npz'
        if not file.exists():
            F=fast.project_many(energy);file.write_bytes(deterministic_npz_bytes({'projection':F}))
        with np.load(file) as a:F=a['projection'].copy()
        arrays[f'fast_energy_{points}']=energy;arrays[f'fast_weight_{points}']=weight;arrays[f'fast_projection_{points}']=F
        print(f'Fast input {points}: {len(energy)} nodes',flush=True)
    subcache=_threshold_cache(channel,config)
    for points,height in ((48,3),(48,4),(24,3)):
        path=subcache/f'panel-p{points}-h{height}.npz'
        if not path.exists():raise ValueError('verified subgap panel must be prepared first: '+str(path))
        with np.load(path,allow_pickle=False) as a:
            arrays[f'subgap_gram_{points}_{height}']=a['gram_signed'].copy()
            arrays[f'subgap_centered_{points}_{height}']=a['centered_signed'].copy()
    for order,grid in ((16,64),(12,64),(16,128)):
        tail=EndpointTailSeries(MassiveTailPackets(p,np.pi/2,0.,order=order,points_per_unit=grid),160.)
        K,bK=tail.centered_integral();G,bG=tail.integral((1.,1.,1.));S=SIGNED_PACKET_MAP
        arrays[f'tail_centered_{order}_{grid}']=K-S@K.conj()@S
        arrays[f'tail_gram_{order}_{grid}']=G+S@G.conj()@S
        arrays[f'tail_bound_{order}_{grid}']=np.array([2*bK,2*bG])
        print(f'Tail input order {order}, grid {grid}',flush=True)
    metadata={'schema':'NSC-PG-GROUP13-INPUT-v1','source_hashes':source_hashes,'config':config,
              'locked_inputs':locked,'channel':channel,'seed_covariance_used_as_input':False,
              'subgap_scope':'existing verified contour panels, before any full-state declaration',
              'tail_scope':'mode approximation plus explicit Fourier/normalization remainder; no stress subtraction'}
    arrays['metadata_json']=np.frombuffer(json.dumps(metadata,sort_keys=True).encode(),np.uint8)
    raw=deterministic_npz_bytes(arrays);digest=hashlib.sha256(raw).hexdigest()
    path=ROOT/f'results/development/artifacts/nsc-pg-group13-input.{digest}.npz';path.write_bytes(raw)
    print(json.dumps({'input':str(path.relative_to(ROOT)),'sha256':digest,'bytes':len(raw)},indent=2))

def assemble(arrays):
    F=arrays['low_projection'];low=assemble_signed_window(arrays['low_energy'],arrays['low_weight'],F,F,arrays['low_source'],arrays['low_projector'])
    S=SIGNED_PACKET_MAP;D=np.diag([-.5,.5,-.5]);windows={}
    for points in (16,24):
        F=arrays[f'fast_projection_{points}'];w=arrays[f'fast_weight_{points}']/(2*np.pi)
        K=np.einsum('n,nai,ij,nbj->ab',w,F,D,F.conj());G=np.einsum('n,nai,nbi->ab',w,F,F.conj())
        windows[points]=(K-S@K.conj()@S,G+S@G.conj()@S)
    K=low.centered+arrays['subgap_centered_48_3']+windows[24][0]+arrays['tail_centered_16_64']
    G=low.gram+arrays['subgap_gram_48_3']+windows[24][1]+arrays['tail_gram_16_64']
    C=.5*G+K
    result={'covariance':C,'car_gram':G,**equal_time_ctp_blocks(C,G)}
    residuals={'CAR':float(np.linalg.norm(G-np.eye(8),2)),'Hermiticity':float(np.linalg.norm(C-C.conj().T,2)),
       'finite_energy_refinement_C':float(np.linalg.norm(windows[24][0]-windows[16][0],2)),
       'finite_energy_refinement_Gram':float(np.linalg.norm(windows[24][1]-windows[16][1],2)),
       'subgap_quadrature_C':float(np.linalg.norm(arrays['subgap_centered_48_3']-arrays['subgap_centered_24_3'],2)),
       'subgap_height_C':float(np.linalg.norm(arrays['subgap_centered_48_3']-arrays['subgap_centered_48_4'],2)),
       'tail_mode_order_C':float(np.linalg.norm(arrays['tail_centered_16_64']-arrays['tail_centered_12_64'],2)),
       'tail_grid_C':float(np.linalg.norm(arrays['tail_centered_16_64']-arrays['tail_centered_16_128'],2)),
       'probe_bulk_correlation_norm':float(np.linalg.norm(C[:4,4:],2)),
       'covariance_eigenvalue_min':float(np.linalg.eigvalsh(C).min()),'covariance_eigenvalue_max':float(np.linalg.eigvalsh(C).max())}
    return result,residuals

def encoded_matrix(value):return np.stack((value.real,value.imag),axis=-1).tolist()

def make_record(path,arrays,metadata):
    result,residuals=assemble(arrays)
    tolerance=3e-9
    checked=('CAR','Hermiticity','finite_energy_refinement_C','finite_energy_refinement_Gram',
             'subgap_quadrature_C','subgap_height_C','tail_mode_order_C','tail_grid_C')
    if any(residuals[k]>tolerance for k in checked):raise ArithmeticError('group13 numerical covariance gate failed')
    if residuals['covariance_eigenvalue_min'] < -tolerance or residuals['covariance_eigenvalue_max'] > 1+tolerance:
        raise ArithmeticError('computed covariance violates CAR bounds')
    if float(arrays['tail_bound_16_64'][0])>tolerance:raise ArithmeticError('Fourier tail representation error is too large')
    source_names=(*PHYSICAL_SOURCES,'src/recursive_horizons/nsc_pg_state_covariance.py',
                  'src/recursive_horizons/nsc_pg_retarded_packets.py',
                  'scripts/derive_nsc_pg_group13_covariance.py','docs/nsc-pg-group13-covariance.md',
                  'tests/test_nsc_pg_group13_covariance.py','scripts/derive_nsc_pg_threshold_projection.py',
                  'scripts/derive_nsc_transmitting_boundary_binding.py','docs/nsc-pg-threshold-projection.md',
                  'tests/test_nsc_pg_threshold_projection.py','tests/test_nsc_complex_horizon_modes.py')
    record={'schema':'NSC-PG-GROUP13-COVARIANCE-v1',
       'status':'PASS: full-energy common-PG covariance of retained group13 on the declared probes; full retained C1b OPEN',
       'source_hashes':{p:sha(p) for p in source_names},
       'input_hashes':{p:sha(p) for p in ('results/development/nsc-pg-massive-mode-resolution.json',
                     'results/development/nsc-compact-matched-restart.json','results/development/nsc-mode-resolved-cauchy-state.json')},
       'payload':{'path':str(path.relative_to(ROOT)),'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'bytes':path.stat().st_size},
       'locked_inputs':metadata['locked_inputs'],'channel':metadata['channel'],
       'domain':{'surface':'same common PG Cauchy surface on the fixed transmitted Dirac background',
                 'basis':['parent spin2','child spin2','orthogonal child-bulk spin2','exterior-bulk spin2'],
                 'field_state':'C_PG = F (C_H direct_sum n_in) F dagger on the full mode space',
                 'probe_complement':'Q=I-J J dagger; C_JQ(E;rho)=F_J C_src (Phi(rho)-J(rho)F_J) dagger',
                 'closed_eight_mode_system':False,'angular_zero_mode_doubled':False,
                 'negative_frequency':'same locked charge-conjugate radial map, not an assigned partner covariance'},
       'energy_account':{'low':[0.,16.],'subgap':[1.,float(np.pi/2)],'fast_mode_window':[16.,160.],
          'tail':[160.,'infinity'],'negative_energy':'included by the physical signed mode/source identity',
          'covariance_assembly':'C_raw=Gram_raw/2+integral centered density; Gram was integrated independently, not filled with I',
          'seed_covariance_used_as_input':False,'LLL_covariance_used_as_substitute':False,
          'reference_or_local_allocation_changed':False},
       'computed':{name:encoded_matrix(value) for name,value in result.items()},
       'verification':residuals,
       'error_account':{'numerical_tolerance':tolerance,
          'tail_Fourier_and_normalization_bound_C':float(arrays['tail_bound_16_64'][0]),
          'tail_Fourier_and_normalization_bound_Gram':float(arrays['tail_bound_16_64'][1]),
          'warm_source_log_norm_bound_above_16':float(np.log(2)-min(np.pi*16/metadata['config']['surface_gravity'],2*np.pi*16/(metadata['config']['omega']*metadata['config']['surface_gravity']))),
          'mode_approximation':'checked by exact-mode comparison at the existing boundary controls and order/grid convergence; not claimed as a uniform theorem on every nonperturbative remainder',
          'eight_probe_compression':'does not replace the full source-to-field mode operator or its eliminated bulk'},
       'gate':{'group13_covariance':'PASS','massive_groups_remaining':31,'full_C1b':'OPEN',
               'transmitting_EndpointBranchJets':'OPEN','Gamma_rest_Weyl':'OPEN','extended_stationarity':'OPEN',
               'physical_stress_computed':False,'metric_timestep':False,'physical_history_selected':False,'PDF_bump':False},
       'comparison':{'fields':'all','float_atol':3e-13,'float_rtol':3e-13,
                     'exact':'structure, strings, integers, booleans, source/input/payload hashes',
                     'routine':'recombine the authenticated full-energy input artifact; no old generator rerun'}}
    return record

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--prepare-inputs',action='store_true');parser.add_argument('--workers',type=int,default=2)
    parser.add_argument('--input');parser.add_argument('--record',action='store_true');parser.add_argument('--check',action='store_true');args=parser.parse_args()
    if args.prepare_inputs:prepare(args.workers);return
    stored=json.loads((ROOT/OUTPUT).read_text()) if args.check else None
    input_name=stored['payload']['path'] if stored else args.input
    if not input_name:raise ValueError('explicit authenticated input artifact required')
    path=ROOT/input_name;raw=path.read_bytes();digest=hashlib.sha256(raw).hexdigest()
    if path.name!=f'nsc-pg-group13-input.{digest}.npz':raise ValueError('input artifact digest differs')
    with np.load(path,allow_pickle=False) as a:arrays={k:a[k].copy() for k in a.files}
    metadata=json.loads(arrays['metadata_json'].tobytes())
    for name,expected in metadata['source_hashes'].items():
        if sha(name)!=expected:raise ValueError('physical input dependency changed: '+name)
    if args.record or args.check:
        record=make_record(path,arrays,metadata)
        if args.check:
            from derive_nsc_transmitting_boundary_binding import compare
            compare(stored,record)
        else:(ROOT/OUTPUT).write_text(json.dumps(record,indent=2,sort_keys=True)+'\n')
        print(json.dumps({'status':record['status'],'verification':record['verification'],'gate':record['gate']},indent=2))
    else:
        _,residuals=assemble(arrays);print(json.dumps(residuals,indent=2))

if __name__=='__main__':main()
