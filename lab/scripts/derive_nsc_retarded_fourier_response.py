#!/usr/bin/env python3
"""Full first-order Fourier-pair reduction from authenticated saved traces.

No field propagation, energy quadrature or finite-window background transform.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import simpson, trapezoid
import derive_nsc_retarded_response_resolution as R

ROOT=Path(__file__).resolve().parents[1]
OUTPUT=ROOT/'results/development/nsc-retarded-fourier-response.json'
OWNED=('scripts/derive_nsc_retarded_fourier_response.py',
       'tests/test_nsc_retarded_fourier_response.py',
       'docs/nsc-retarded-fourier-response.md')
INPUTS=(R.OUTPUT,'scripts/derive_nsc_retarded_response_resolution.py',
        'docs/nsc-incoming-fourier-matching.md')


def transform_modes(energies,times,delta_trace,clock_shift,*,rule='simpson'):
    """dA[o,i,spin,source]=integral dz exp(i Eo z) dF(z,Ei).

    All three input-source columns are unweighted. Only the perturbation is
    integrated numerically; the analytic background plane wave is not windowed.
    """
    E=np.asarray(energies);t=np.asarray(times);dF=np.asarray(delta_trace,complex)
    if (np.iscomplexobj(E) or np.iscomplexobj(t) or E.ndim!=1 or t.ndim!=1
        or not np.isfinite(E).all() or not np.isfinite(t).all()
        or len(t)<3 or np.any(np.diff(t)<=0) or not np.isfinite(clock_shift)):
        raise ValueError('finite real frequencies, clock and ordered time nodes required')
    if dF.ndim!=3 or dF.shape[0]!=len(t) or dF.shape[1]!=2 or dF.shape[2]%len(E) or not np.isfinite(dF).all():
        raise ValueError('unweighted time/spin/source columns with matching energies required')
    if rule not in ('simpson','trapezoid'):raise ValueError('declared time integration rule required')
    integrate=simpson if rule=='simpson' else trapezoid
    columns=dF.shape[2]//len(E)
    result=[]
    for energy in E:
        value=integrate(np.exp(1j*energy*(t+clock_shift))[:,None,None]*dF,x=t,axis=0)
        result.append(value.reshape(2,len(E),columns).transpose(1,0,2))
    return np.asarray(result)


def covariance_pairs(A,stationary,source):
    """Exact pair algebra after analytic source-energy delta collapse."""
    A=np.asarray(A,complex);f=np.asarray(stationary,complex);C=np.asarray(source,complex)
    if f.ndim!=3 or f.shape[1]!=2 or A.shape!=(len(f),len(f),2,f.shape[2]) or C.shape!=(len(f),f.shape[2],f.shape[2]):
        raise ValueError('matching pair, stationary and full source-fiber shapes required')
    if not all(np.isfinite(x).all() for x in (A,f,C)):
        raise ValueError('finite full coherent matrices required')
    return np.asarray([[A[o,i]@C[i]@f[i].conj().T+f[o]@C[o]@A[i,o].conj().T
                        for i in range(len(f))] for o in range(len(f))])


def maxabs(a):return float(np.max(abs(a)))
def pack(a):return {'real':a.real.tolist(),'imag':a.imag.tolist()}
def signature():return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in (*OWNED,*INPUTS)}


def make_record():
    parent=R.check()  # authenticated array replay only, never calls a producer
    if not parent['status'].startswith('PASS:'):
        raise ValueError('completed declared resolution indicators required')
    with np.load(ROOT/parent['payload']['path'],allow_pickle=False) as data:
        a={k:data[k] for k in data.files if k!='metadata_json'}
    E=a['source/energies'];C=a['source/source'];P=a['source/projector']
    S1=R.P.chart_coordinates(1.)[1]
    rows={};pairs={};CAR={};phase_error=0.
    for case in ('1601_128','3201_128','3201_256'):
        row=R.case_arrays(a,case);t=row['times'];dF=row['trace_tangent']
        x,phi=R.subgrid(a['reference/x'],a['reference/phi'],len(row['x']))
        initial=row['restriction']@phi[row['sample_rows']]
        f=initial.reshape(2,len(E),3).transpose(1,0,2)*np.exp(1j*E*S1)[:,None,None]
        A=transform_modes(E,t,dF,S1)
        B=covariance_pairs(A,f,C);G=covariance_pairs(A,f,P)
        Bt=covariance_pairs(transform_modes(E,t,dF,S1,rule='trapezoid'),f,C)
        # Independent origin choice: tau convention and exact pair phase law.
        f_tau=initial.reshape(2,len(E),3).transpose(1,0,2)
        Btau=covariance_pairs(transform_modes(E,t,dF,0.),f_tau,C)
        expected=np.exp(1j*(E[:,None]-E[None,:])*S1)[:,:,None,None]*Btau
        phase=maxabs(B-expected);phase_error=max(phase_error,phase)
        peak=np.unravel_index(abs(B).argmax(),B.shape);norm=maxabs(B)
        rows[case]={'covariance_pair_max':norm,'CAR_pair_residual':maxabs(G),
            'CAR_over_response':maxabs(G)/norm,'time_quadrature_difference':maxabs(B-Bt),
            'quadrature_difference_over_response':maxabs(B-Bt)/norm,
            'pair_hermiticity_residual':maxabs(B-B.transpose(1,0,3,2).conj()),
            'phase_origin_residual':phase,'stationary_initial_column_residual':maxabs(initial-row['trace'][0]),
            'closed_input_tangent_column':maxabs(A[:,:,:,2]),
            'maximum_pair':{'output_energy':float(E[peak[0]]),'input_energy':float(E[peak[1]]),
                'output_spin':int(peak[2]),'input_spin':int(peak[3]),
                'value_real':float(B[peak].real),'value_imag':float(B[peak].imag)},
            'covariance_pairs':pack(B),'CAR_pairs':pack(G),'mode_pair_response':pack(A)}
        pairs[case]=B;CAR[case]=G
    comparison={}
    for name,coarse,fine,limit in (('spatial','1601_128','3201_128',.25),
                                  ('temporal','3201_128','3201_256',.1)):
        error=maxabs(pairs[coarse]-pairs[fine]);ratio=error/maxabs(pairs[fine])
        comparison[name]={'difference':error,'relative_difference':ratio,'indicator_limit':limit,'passes':ratio<limit}
    fine=rows['3201_256']
    algebra=max(r[k] for r in rows.values() for k in ('pair_hermiticity_residual','phase_origin_residual',
                 'stationary_initial_column_residual','closed_input_tangent_column'))
    passed=(all(v['passes'] for v in comparison.values()) and algebra<3e-11
            and fine['CAR_over_response']<.01 and fine['quadrature_difference_over_response']<.01)
    return {'schema':'NSC-RETARDED-FOURIER-RESPONSE-v1','accountable_author':'Douglas Ek',
        'status':'PASS: numerical first-order pair indicators; rigorous matching certificate OPEN' if passed else 'OPEN: first-order pair indicators unresolved',
        'source_sha256':signature(),'input_payload':parent['payload'],'energies':E.tolist(),
        'convention':{'z':'tau+S1','S1':float(S1),
            'double_transform':'integral dz dz_prime exp(i Eo z) deltaC(z,z_prime) exp(-i Ei z_prime)',
            'deltaA':'integral dz exp(i Eo z) deltaF(z,Ei)',
            'source_energy_integral':'analytically collapsed at first order, not discretely approximated',
            'energy_quadrature_weights_used':False,'angular_multiplicity_used':False,
            'background':'analytic stationary amplitude from authenticated reference columns',
            'response_window':[0.,.3],'exact_response_support_contained_in':[.0525,.2475],
            'numerical_trace_clipped':False},
        'per_case':rows,'comparison':comparison,
        'indicators':{'algebra_max':algebra,'algebra_tolerance':3e-11,
            'CAR_over_response_limit':.01,'quadrature_over_response_limit':.01,
            'differences_are_error_bounds':False},
        'scope':{'single_angular_block':'14_1','fixed_original_source_law':True,
            'linearized_C0_numerical_response':'nonzero' if fine['covariance_pair_max'] else 'zero at selected pairs',
            'rigorous_continuum_response_error_bound':None,'exact_fixed_C0_exclusion_certificate':'OPEN',
            'matching_for_other_variations':'OPEN','global_constraints':'OPEN','extended_stationarity':'OPEN',
            'field_or_radial_solves':0,'metric_timestep':False,'stress_computed':False,
            'added_action_or_fit':False,'CAR_residual_subtracted':False}}


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--check',action='store_true');args=parser.parse_args()
    result=make_record()
    if args.check:
        if json.loads(OUTPUT.read_text())!=result:raise ValueError('Fourier pair receipt differs from authenticated replay')
    else:OUTPUT.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'status':result['status'],'comparison':result['comparison'],
        'fine':{k:v for k,v in result['per_case']['3201_256'].items() if k not in ('covariance_pairs','CAR_pairs','mode_pair_response')}},indent=2))


if __name__=='__main__':main()
