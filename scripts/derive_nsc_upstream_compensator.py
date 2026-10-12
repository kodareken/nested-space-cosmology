#!/usr/bin/env python3
"""Six bounded slab solves, one real geometry triple, saved-array replay."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import argparse
import json
from pathlib import Path
import signal
import sys
import time

import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import derive_nsc_retarded_radial_response as D
from recursive_horizons.nsc_retarded_radial_response import archived_amplitudes,pulse_fourier
from recursive_horizons.nsc_upstream_compensator import (
    operator_response,fit_selected,combine,covariance_response,small_amplitude_range,vertex_identities,DIRECTIONS,MAX_STEP)

P=D.R.P;F=D.F
OUTPUT='results/development/nsc-upstream-compensator.json'
ARTIFACT='results/development/artifacts/nsc-upstream-compensator.npz'
OWNED=('src/recursive_horizons/nsc_upstream_compensator.py','scripts/derive_nsc_upstream_compensator.py',
       'tests/test_nsc_upstream_compensator.py','docs/nsc-upstream-compensator.md')
EXTRA=(D.OUTPUT,'results/development/nsc-upstream-compensation-rank.json','docs/nsc-upstream-metric-compensation.md',
       'docs/nsc-upstream-compensation-rank.md','src/recursive_horizons/nsc_incoming_cauchy_jets.py')
FAMILIES=('14_1','14_-1','13_1','0_1','24_1')
CASES=(('14_1','coarse',2e-10,2e-12),('14_1','fine',2e-13,2e-15),
       *((name,'fine',2e-13,2e-15) for name in FAMILIES[1:]))
TOL=3e-11;CPU_CAP=60.


def signature():return {**D.signature(),**{p:P.digest(p) for p in (*OWNED,*EXTRA)}}


def inputs():
    radial=D.check()
    history=json.loads((ROOT/P.RECORDS[0]).read_text());modes=json.loads((ROOT/P.RECORDS[1]).read_text())
    for record in (history,modes):
        if P.digest(record['payload']['path'])!=record['payload']['sha256']:raise ValueError('original source/history payload changed')
    result={}
    with np.load(ROOT/history['payload']['path'],allow_pickle=False) as h,np.load(ROOT/modes['payload']['path'],allow_pickle=False) as m:
        for name in FAMILIES:
            data={k:m[name+'/'+k] for k in ('energies','mass','angular','source','projector','weights')}
            x=h[name+'/x'];phi=h[name+'/reference_field'];indices=np.flatnonzero(x==1.)
            if len(x)!=801 or len(indices)!=1 or phi.shape!=(1602,12):raise ValueError('owned reference layout changed')
            f=archived_amplitudes(x,phi,data['energies'],int(indices[0]))
            C=f@data['source']@f.swapaxes(-1,-2).conj()
            result[name]={**data,'stationary':f,'covariance':C}
    with np.load(ROOT/radial['payload']['path'],allow_pickle=False) as a:
        B=a['fine/response'];f=a['reference/incoming']
        baseline=np.array([[np.linalg.solve(f[i,:,:2].T,B[o,i,:,:2].T).T for i in range(4)] for o in range(4)])
    return radial,result,baseline,[history['payload'],modes['payload'],radial['payload']]


def support_guards(original):
    left,right=203/200,41/40;lo,hi=original['radial_support']
    if not lo<1<left<right<hi<1.2:raise ValueError('compensator union does not retain prior support')
    initial=original['axial_center']-.06-P.chart_coordinates(hi)[1]
    final=original['axial_center']+.06-P.chart_coordinates(lo)[1]
    if not 0<initial<final<.3:raise ValueError('union support changes the declared reference past or time coverage')
    return {'compensator_support':[left,right],'support_exact_rationals':[[203,200],[41,40]],
            'union_radial_support':[lo,hi],'union_PG_time_support':[initial,final],
            'initial_time':0.,'final_time':.3,'upstream_right_edge':1.2,'incoming_jet_neighborhood_rho':[1.,left],
            'all_compensator_incoming_jets':'exact zero by vanishing on an open neighborhood of rho1',
            'unchanged_initial_inflow_source_tangents':True}


def row_arrays(a,case):return {k[len(case)+1:]:v for k,v in a.items() if k.startswith(case+'/')}
def family_arrays(a,name):return {k[len(name)+7:]:v for k,v in a.items() if k.startswith('input/'+name+'/')}


def evaluate_case(a,case,c):
    name=case.split('/')[0];data=family_arrays(a,name);U=a[case+'/U'];K=a[case+'/K']
    total=combine(K,c);original=K[0];C=data['covariance']
    before=covariance_response(original,C);after=covariance_response(total,C)
    fullF=data['stationary'];Cs=data['source']
    source_response=np.array([[total[o,i]@fullF[i]@Cs[i]@fullF[i].conj().T+
                               fullF[o]@Cs[o]@fullF[o].conj().T@total[i,o].conj().T for i in range(4)] for o in range(4)])
    algebra={'unitarity':F.maxabs(U@U.swapaxes(-1,-2).conj()-np.eye(2)),
             'all_direction_CAR':F.maxabs(K+K.transpose(0,2,1,4,3).conj()),
             'compensated_CAR':F.maxabs(total+total.transpose(1,0,3,2).conj()),
             'covariance_pair_Hermiticity':F.maxabs(after-after.transpose(1,0,3,2).conj()),
             'full_source_vs_slab_reduction':F.maxabs(source_response-after)}
    pairs=[]
    for o in range(4):
        for i in range(4):
            pairs.append({'output_index':o,'input_index':i,'output_energy':float(data['energies'][o]),'input_energy':float(data['energies'][i]),
                          'original_covariance_norm':F.maxabs(before[o,i]),'compensated_covariance_norm':F.maxabs(after[o,i]),
                          'residual_above_numerical_tolerance':F.maxabs(after[o,i])>=TOL})
    row={'original_K_max':F.maxabs(original),'compensated_K_max':F.maxabs(total),
         'original_covariance_max':F.maxabs(before),'compensated_covariance_max':F.maxabs(after),
         'original_radius_identically_zero':bool(data['angular'].item()==0.),
         'created_response_in_original_zero_channel':bool(data['angular'].item()==0. and F.maxabs(after)>=TOL),
         'algebra':algebra,'stationary_source_coisometry_residual':F.maxabs(fullF@fullF.swapaxes(-1,-2).conj()-np.eye(2)),
         'pairs':pairs,'original_covariance_pairs':F.pack(before),'compensated_covariance_pairs':F.pack(after),
         'compensated_operator_pairs':F.pack(total)}
    return row,total,after


def analyze(arrays,meta):
    rows={};checks={};selected={};c=arrays.get('coefficients')
    if c is None:return rows,checks,selected,False
    outputs={}
    for case in meta['completed']:
        row,K,C=evaluate_case(arrays,case,c);rows[case]=row;outputs[case]=(K,C)
        checks.update({case+'/'+k:v for k,v in row['algebra'].items()})
    fine=outputs['14_1/fine'];coarse=outputs['14_1/coarse']
    checks.update(selected_operator_residual=F.maxabs(fine[0][3,3]),selected_covariance_residual=F.maxabs(fine[1][3,3]),
        fixed_fine_coefficient_coarse_operator=F.maxabs(coarse[0]-fine[0]),
        fixed_fine_coefficient_coarse_covariance=F.maxabs(coarse[1]-fine[1]),
        original_radius_vs_b622632=F.maxabs(arrays['14_1/fine/K'][0]-arrays['baseline/K']))
    fit=meta['fit'];checks.update({key:fit[key] for key in ('selected_antiHermitian','selected_trace','Pauli_coefficient_imaginary','linear_system_residual')})
    selected={'family':'14_1','energy_index':3,'energy':float(arrays['input/14_1/energies'][3]),'coefficients_N_a_r':c.tolist(),
              'fit':fit,'same_coefficients_all_cases':True,'coefficient_frozen_before_validation':meta['coefficient_frozen_before_validation'],
              'small_amplitude_chart_range':small_amplitude_range(c)}
    local=all(v=='0' for v in meta['vertex_proof']['raw_vertex_residuals']) and meta['vertex_proof']['Weyl_midpoint_residual']=='0'
    passed=len(meta['completed'])==6 and max(checks.values())<TOL and local and not meta['cpu_budget_exceeded']
    return rows,checks,selected,passed


def record_from_arrays(arrays,meta,digest,size):
    rows,checks,selected,passed=analyze(arrays,meta)
    validation={case:row['compensated_covariance_max'] for case,row in rows.items() if case.endswith('/fine')}
    unresolved=any(value>=TOL for value in validation.values())
    return {'schema':'NSC-UPSTREAM-COMPENSATOR-v1','accountable_author':'Douglas Ek',
        'numerical_status':'PASS: common geometry calibration and numerical controls' if passed else 'OPEN: bounded numerical controls unresolved',
        'matching_status':'OPEN: frozen common triple leaves uncalibrated covariance responses' if unresolved else 'OPEN: finite sampled checks do not establish global matching',
        'control':{'directions':list(DIRECTIONS),'families':list(FAMILIES),'cases_completed':meta['completed'],'maximum_radial_step':MAX_STEP,
                   'tolerance':TOL,'upstream':meta.get('upstream'),'support':meta.get('support'),'energy_weights_used':False},
        'selected_compensation':selected,'per_case':rows,'residuals':checks,'vertex_proof':meta.get('vertex_proof'),
        'fine_channel_residual_maxima':validation,
        'runtime':{k:meta[k] for k in ('CPU_seconds','wall_seconds','CPU_cap_seconds','cpu_budget_exceeded','solver','stop_reason')},
        'scope':{'one_common_geometric_triple':True,'channel_or_energy_dependent_fit':False,'new_matrix_force':False,
                 'source_or_coupling_changed':False,'archived_columns_normalized':False,'incoming_compensator_jets':'zero by support',
                 'global_C0_matching':'OPEN','physical_preparation_accuracy':'inherited qualifications; no new preparation bound',
                 'action_constraints_and_stationarity':'OPEN','metric_timestep':False,'old_producers_rerun':False},
        'source_hashes':{p:meta['signature'][p] for p in OWNED},'input_hashes':{p:h for p,h in meta['signature'].items() if p not in OWNED},
        'input_payloads':meta.get('input_payloads',[]),'payload':{'path':ARTIFACT,'sha256':digest,'bytes':size},
        'reproducer':'python3 scripts/derive_nsc_upstream_compensator.py --check'}


def publish(arrays,meta):
    content={**arrays,'metadata_json':np.frombuffer(json.dumps(meta,sort_keys=True).encode(),np.uint8)}
    raw=P.deterministic_npz_bytes(content);(ROOT/ARTIFACT).write_bytes(raw)
    record=record_from_arrays(content,meta,P.sha256(raw).hexdigest(),len(raw))
    (ROOT/OUTPUT).write_text(json.dumps(record,indent=2,sort_keys=True,allow_nan=False)+'\n')
    return record


def run():
    begin=time.process_time();wall=time.monotonic();before=signature();arrays={}
    meta={'signature':before,'completed':[],'solver':{},'CPU_cap_seconds':CPU_CAP,'cpu_budget_exceeded':False,
          'coefficient_frozen_before_validation':False,'stop_reason':'checkpoint'}
    def stop(signum,frame):raise P.CPUExceeded('common compensator CPU budget exhausted')
    previous=signal.signal(signal.SIGPROF,stop);signal.setitimer(signal.ITIMER_PROF,CPU_CAP)
    def checkpoint(reason):
        meta.update(CPU_seconds=time.process_time()-begin,wall_seconds=time.monotonic()-wall,stop_reason=reason)
        if signature()!=before:raise ValueError('compensator owner/input changed during run')
        return publish(arrays,meta)
    try:
        radial,families,baseline,payloads=inputs();meta['input_payloads']=payloads;meta['upstream']=radial['control']['upstream']
        guard=radial['control']['support'];meta['support']=support_guards(guard);meta['vertex_proof']=vertex_identities()
        arrays['baseline/K']=baseline
        for name,data in families.items():arrays.update({'input/'+name+'/'+k:v for k,v in data.items()})
        for name,label,rtol,atol in CASES:
            if name!='14_1' and not meta['coefficient_frozen_before_validation']:raise ValueError('fine coefficients must freeze before validation')
            data=families[name];E=data['energies'];what=pulse_fourier(E[:,None]-E[None,:],guard['axial_center'],guard['axial_halfwidth'],192)
            case=name+'/'+label;t0=time.process_time()
            U,B,K,info=operator_response(meta['upstream'],E,data['mass'].item(),data['angular'].item(),what,rtol=rtol,atol=atol)
            info['CPU_seconds']=time.process_time()-t0;meta['solver'][case]=info
            arrays.update({case+'/U':U,case+'/B':B,case+'/K':K,case+'/what':what});meta['completed'].append(case)
            if case=='14_1/fine':
                c,fit=fit_selected(K);arrays['coefficients']=c;meta['fit']=fit;meta['coefficient_frozen_before_validation']=True
            checkpoint('completed '+case);print(case+' CPU '+str(info['CPU_seconds']),flush=True)
        meta['stop_reason']='completed exactly six prescribed short coupled solves'
    except P.CPUExceeded:
        meta['cpu_budget_exceeded']=True;meta['stop_reason']='OPEN: CPU cap reached; no automatic extension'
    finally:
        signal.setitimer(signal.ITIMER_PROF,0);signal.signal(signal.SIGPROF,previous)
    return checkpoint(meta['stop_reason'])


def check():
    old=json.loads((ROOT/OUTPUT).read_text());radial,families,baseline,payloads=inputs()
    if P.digest(ARTIFACT)!=old['payload']['sha256']:raise ValueError('compensator artifact changed')
    with np.load(ROOT/ARTIFACT,allow_pickle=False) as a:arrays={k:a[k] for k in a.files}
    meta=json.loads(arrays['metadata_json'].tobytes())
    if meta['signature']!=signature() or meta['input_payloads']!=payloads:raise ValueError('compensator input binding changed')
    for name,data in families.items():
        for k,v in data.items():
            if not np.array_equal(arrays['input/'+name+'/'+k],v):raise ValueError('saved original source changed')
    if not np.array_equal(arrays['baseline/K'],baseline):raise ValueError('original radial comparison changed')
    if 'coefficients' in arrays:
        c,fit=fit_selected(arrays['14_1/fine/K'])
        if not np.array_equal(c,arrays['coefficients']) or fit!=meta['fit']:raise ValueError('single real calibration replay differs')
    for case in meta['completed']:
        K=np.einsum('doiab,icb->doiac',arrays[case+'/B'],arrays[case+'/U'].conj())
        if not np.array_equal(K,arrays[case+'/K']):raise ValueError('operator response contraction differs')
    record=record_from_arrays(arrays,meta,P.digest(ARTIFACT),(ROOT/ARTIFACT).stat().st_size)
    if record!=old:raise ValueError('compensator saved receipt differs')
    return record


def main():
    parser=argparse.ArgumentParser();mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--run',action='store_true');mode.add_argument('--check',action='store_true');args=parser.parse_args()
    record=run() if args.run else check()
    print(json.dumps({k:record[k] for k in ('numerical_status','matching_status','selected_compensation','fine_channel_residual_maxima','runtime')},indent=2))


if __name__=='__main__':main()
