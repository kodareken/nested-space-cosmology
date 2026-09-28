#!/usr/bin/env python3
"""Saved-array localization and local frozen-generator diagnostics; no evolution."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys
import numpy as np
from scipy.integrate import quad

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
sys.path.insert(0,str(ROOT/'scripts'))
import derive_nsc_retarded_compatible_response as pilot
from recursive_horizons.nsc_transmitting_history_jets import generator_direction

OUTPUT='results/development/nsc-retarded-spatial-diagnosis.json'
SOURCES=('scripts/derive_nsc_retarded_spatial_diagnosis.py','tests/test_nsc_retarded_spatial_diagnosis.py','docs/nsc-retarded-spatial-diagnosis.md')
INPUTS=('results/development/nsc-retarded-compatible-response.json','scripts/derive_nsc_retarded_compatible_response.py',
    'results/development/nsc-transmitting-history-modes.json','results/development/nsc-pg-ctp-mode-jets.json',
    'src/recursive_horizons/nsc_compatible_history_geometry.py','src/recursive_horizons/nsc_transmitting_history_jets.py',
    'src/recursive_horizons/nsc_transmitting_history_modes.py','src/recursive_horizons/nsc_prepared_history_jets.py',
    'src/recursive_horizons/nsc_ks_spacetime_variation.py','src/recursive_horizons/nsc_transmitting_dirac_domain.py',
    'src/recursive_horizons/nsc_pg_ks_metric_pullback.py','src/recursive_horizons/nsc_transmitting_resolvent.py',
    'src/recursive_horizons/nsc_lorentzian.py')


def digest(path):return sha256((ROOT/path).read_bytes()).hexdigest()
def peak(value):return tuple(int(v) for v in np.unravel_index(np.argmax(abs(value)),value.shape))
def norm(value):return float(np.linalg.norm(value))
def maxabs(value):return float(np.max(abs(value)))
def complex_pair(value):return [float(value.real),float(value.imag)]


def sbp_l2(value,x):
    """Owned diagonal SBP(4,2) norm, summed over all finite source columns."""
    w=np.full(len(x),float(x[1]-x[0]));edge=np.array([17,59,43,49])/48
    w[:4]*=edge;w[-4:]*=edge[::-1]
    return float(np.sqrt(np.sum(np.tile(w,2)[:,None]*abs(value)**2)))


def localize(arrays):
    E=arrays['source/energies'];result={};cases=('401_64','801_64','801_128')
    def label(row,column,x):return {'rho':float(x[row%len(x)]),'characteristic_spin':row//len(x),
                                  'energy':float(E[column//3]),'source_column':column%3}
    for case in cases:
        x=arrays[case+'/x'];n=len(x);t=arrays[case+'/times'];Y=arrays[case+'/field_tangent']
        drift=arrays[case+'/field']-arrays[case+'/stationary_field']
        T=arrays[case+'/trace_tangent'];tdrift=arrays[case+'/trace']-arrays[case+'/stationary_trace']
        row={}
        for name,value in [('field_tangent',Y),('stationary_field_drift',drift)]:
            i,j=peak(value);row[name]={'maximum':float(abs(value[i,j])),'location':label(i,j,x),'SBP_L2':sbp_l2(value,x)}
        for name,value in [('trace_tangent',T),('stationary_trace_drift',tdrift)]:
            i,s,j=peak(value);row[name]={'maximum':float(abs(value[i,s,j])),
                'location':{'PG_time':float(t[i]),'canonical_spin':s,'energy':float(E[j//3]),'source_column':j%3}}
        row['final_boundary_tangents']={'left':maxabs(Y[[0,n]]),'right':maxabs(Y[[n-1,2*n-1]])}
        row['spin_shape']=[]
        for s in range(2):
            square=np.sum(abs(Y[s*n:(s+1)*n])**2,axis=1)
            row['spin_shape'].append({'characteristic_spin':s,'squared_column_norm_peak_rho':float(x[np.argmax(square)]),
                                      'unweighted_squared_norm_centroid':float(np.sum(x*square)/np.sum(square))})
        result[case]=row
    x=arrays['401_64/x'];t=arrays['401_64/times'];n=len(x)
    coarse=arrays['401_64/field_tangent'].reshape(2,n,12)
    fine=arrays['801_64/field_tangent'].reshape(2,801,12)[:,::2]
    difference=coarse-fine;i,j=peak(difference.reshape(2*n,12));total=float(np.sum(abs(difference)**2))
    spin_squares=np.sum(abs(difference)**2,axis=(1,2))
    spin0_density=np.sum(abs(difference[0])**2,axis=1);mask=(x>=.95)&(x<=1.05)
    spatial={'maximum':float(abs(difference.reshape(2*n,12)[i,j])),'location':label(i,j,x),
        'coarse_at_maximum':complex_pair(coarse.reshape(2*n,12)[i,j]),'fine_at_maximum':complex_pair(fine.reshape(2*n,12)[i,j]),
        'SBP_L2_difference':sbp_l2(difference.reshape(2*n,12),x),'SBP_L2_fine':sbp_l2(fine.reshape(2*n,12),x),
        'SBP_L2_relative':sbp_l2(difference.reshape(2*n,12),x)/sbp_l2(fine.reshape(2*n,12),x),
        'unweighted_error_squared_norm_fraction_by_spin':(spin_squares/total).tolist(),
        'spin0_error_fraction_in_rho_0_95_to_1_05':float(np.sum(spin0_density[mask])/np.sum(spin0_density)),
        'unweighted_error_fraction_by_energy':(np.sum(abs(difference.reshape(2,n,4,3))**2,axis=(0,1,3))/total).tolist(),
        'unweighted_error_fraction_by_source':(np.sum(abs(difference.reshape(2,n,4,3))**2,axis=(0,1,2))/total).tolist()}
    tc=arrays['401_64/trace_tangent'];tf=arrays['801_64/trace_tangent'];dt=tc-tf
    ti,s,j=peak(dt);weights=np.full(len(t),float(t[1]-t[0]));weights[[0,-1]]/=2
    tl2=lambda value:float(np.sqrt(np.sum(weights[:,None,None]*abs(value)**2)))
    trace={'maximum':float(abs(dt[ti,s,j])),'location':{'PG_time':float(t[ti]),'canonical_spin':s,'energy':float(E[j//3]),'source_column':j%3},
        'coarse_norm_at_maximum_time':maxabs(tc[ti]),'fine_norm_at_maximum_time':maxabs(tf[ti]),
        'time_L2_difference':tl2(dt),'time_L2_fine':tl2(tf),'time_L2_relative':tl2(dt)/tl2(tf),
        'squared_time_norm_fraction_after_0_21':float(np.sum(weights[t>.21,None,None]*abs(dt[t>.21])**2)/tl2(dt)**2),
        'unweighted_error_fraction_by_spin':(np.sum(abs(dt)**2,axis=(0,2))/np.sum(abs(dt)**2)).tolist()}
    Kc=arrays['401_64/covariance_tangent'];Kf=arrays['801_64/covariance_tangent'];dk=Kc-Kf;ki,kj=peak(dk)
    kernel={'maximum':float(abs(dk[ki,kj])),'time1':float(t[ki//2]),'time2':float(t[kj//2]),
            'canonical_spin1':ki%2,'canonical_spin2':kj%2,'Frobenius_relative':norm(dk)/norm(Kf)}
    bd=(arrays['401_64/field']-arrays['401_64/stationary_field']).reshape(2,n,12)-(arrays['801_64/field']-arrays['801_64/stationary_field']).reshape(2,801,12)[:,::2]
    btd=(arrays['401_64/trace']-arrays['401_64/stationary_trace'])-(arrays['801_64/trace']-arrays['801_64/stationary_trace'])
    correlation=lambda a,b:float(abs(np.vdot(a,b))/(norm(a)*norm(b)))
    baseline={'field_discrepancy_normalized_correlation':correlation(difference,bd),'trace_discrepancy_normalized_correlation':correlation(dt,btd),
        'baseline_drift_difference_at_tangent_maximum':complex_pair(bd.reshape(2*n,12)[i,j]),
        'baseline_drift_difference_norm_at_tangent_maximum':float(abs(bd.reshape(2*n,12)[i,j])),
        'baseline_field_discrepancy_squared_norm_fraction_by_spin':(np.sum(abs(bd)**2,axis=(1,2))/np.sum(abs(bd)**2)).tolist(),
        'baseline_trace_discrepancy_squared_norm_fraction_by_spin':(np.sum(abs(btd)**2,axis=(0,2))/np.sum(abs(btd)**2)).tolist(),
        'temporal_field_tangent_max_difference':maxabs(arrays['801_64/field_tangent']-arrays['801_128/field_tangent']),
        'temporal_baseline_field_max_difference':maxabs(arrays['801_64/field']-arrays['801_128/field'])}
    # Saved-row comparison only: this does not shift, interpolate or correct a field.
    roi=np.flatnonzero((x>.8)&(x<1.12));shifts=[]
    for shift in range(-10,11):
        c=coarse[0,roi+shift];f=fine[0,roi]
        shifts.append({'coarse_index_shift':shift,'relative_shape_difference':norm(c-f)/norm(f)})
    displacement={'unshifted':shifts[10],'best_integer_shift':min(shifts,key=lambda r:r['relative_shape_difference']),
                  'interpretation':'diagnostic index comparison only; no output, coordinate or physical profile adjusted'}
    return {'per_case':result,'spatial_field':spatial,'spatial_trace':trace,'spatial_covariance':kernel,
            'baseline_discriminator':baseline,'translation_discriminator':displacement,
            'norm_scope':'finite unweighted source-column norms; SBP spatial and trapezoidal time norms are descriptive diagnostics, not physical state norms or error bounds'}


def local_forcing_diagnostic(x801,phi801,data,support):
    direction=pilot.CompatibleRadiusDirection(pilot.axial_bump(support['axial_center']),lambda z,n:0.,.007,.03)
    grids={}
    for points in (401,801):
        stride=800//(points-1);x=x801[::stride];phi=phi801.reshape(2,801,12)[:,::stride].reshape(2*points,12)
        grids[points]=(x,phi,pilot.FourthOrderModePropagator(x,data['mass'].item(),data['angular'].item()),pilot.CachedZeroRadiusProvider(x,direction))
    result=[]
    for t in (.11,.15,.19):
        row={'PG_time':t,'grids':{}};samples={}
        for n,(x,phi,owner,provider) in grids.items():
            metric=provider.values(t,x);directions=provider.log_directions(t,x,metric)[0]
            L,checks=owner.generator(metric);dL,tangent_residual=generator_direction(owner,metric,directions)
            phase=np.exp(-1j*np.repeat(data['energies'],3)*t)
            F=dL@(phi*phase[None,:]);LF=L@F;samples[n]=(F,LF)
            geometry_basis=directions[3]/metric[3]  # delta r/r^2, not delta log r
            theta=2*np.pi*np.fft.fftfreq(n);group_factor=(4*np.cos(theta)-np.cos(2*theta))/3
            power=abs(np.fft.fft(geometry_basis))**2;total=float(np.sum(power))
            row['grids'][str(n)]={'F_maximum':maxabs(F),'LF_maximum':maxabs(LF),
                'geometry_power_fraction_group_error_above_0_25':float(np.sum(power[abs(group_factor-1)>.25])/total),
                'geometry_power_fraction_reversed_group_factor':float(np.sum(power[group_factor<0])/total),
                'generator_flux_algebra_residual':float(checks['norm_flux_algebra']),'generator_tangent_flux_algebra_residual':float(tangent_residual)}
        common=lambda value:value.reshape(2,801,12)[:,::2].reshape(802,12)
        f,lf=map(common,samples[801]);delta=samples[401][1]-lf
        row.update(F_common_node_difference=maxabs(samples[401][0]-f),
                   LF_common_node_difference=maxabs(delta),LF_fine_common_norm=maxabs(lf),
                   LF_relative_difference=maxabs(delta)/maxabs(lf))
        result.append(row)
    return {'times':result,'interior_group_factor':'(4*cos(theta)-cos(2*theta))/3',
        'forcing':'F=dL*phi_archive*exp(-iE*t); LF applies one owned spatial generator, no time step',
        'interpretation':'local frozen-coefficient dispersion discriminator; not a continuum/global causal proof',
        'geometry_spectrum':'DFT of the same full-grid delta_r/r^2; finite-domain diagnostic, not new periodic boundary data'}


def causal_support_diagnostic(arrays,support):
    """Radius-only fixed cones; numerical edge estimates, no directed bound."""
    upper=support['radial_support'][1]
    X,error=quad(lambda rho:1/(float(pilot.geometry(rho)[0])**2-1),1.,upper,epsabs=2e-14,epsrel=2e-14)
    start,stop=.09-X,.21+X
    cases={};E=arrays['source/energies']
    for case in ('401_64','801_64','801_128'):
        t=arrays[case+'/times'];v=arrays[case+'/trace_tangent'];mask=(t<start)|(t>stop)
        i,spin,col=peak(v[mask]);square=abs(v)**2
        weights=np.full(len(t),float(t[1]-t[0]));weights[[0,-1]]/=2
        cases[case]={'outside_window_maximum':maxabs(v[mask]),
            'location':{'PG_time':float(t[mask][i]),'canonical_spin':spin,'energy':float(E[col//3]),'source_column':col%3},
            'outside_unweighted_squared_norm_fraction':float(np.sum(square[mask])/np.sum(square)),
            'outside_trapezoidal_time_norm_fraction':float(np.sum(weights[mask,None,None]*square[mask])/np.sum(weights[:,None,None]*square)),
            'nearest_sample_distance_from_estimated_window_edge':float(min(np.min(abs(t-start)),np.min(abs(t-stop))))}
    return {'X':float(X),'quadrature_error_estimate_not_bound':float(error),'PG_time_window':[float(start),float(stop)],
        'formula':'X=integral_1^rho_plus d_rho/(beta²-1); tau in [0.09-X,0.21+X]',
        'fast_arrival_identity':'integral[1/(beta+1)-beta/(beta²-1)] = -X',
        'slow_arrival_identity':'integral[1/(beta-1)-beta/(beta²-1)] = +X',
        'conditions':['initial tangent zero before the source','source compact with z-S1 in[0.09,0.21]',
            'both characteristic velocities left-going; only upstream rho>1 can influence Sigma',
            'radius-only variation leaves principal cones fixed; lower-order spin coupling stays inside those cones'],
        'per_case':cases,'interpretation':'outside-window saved response is numerical leakage under these continuum assumptions',
        'directed_edge_certificate':False,'discrete_finite_propagation_claimed':False}


def make_record():
    receipt=pilot.check()  # authenticated saved arrays only, never pilot.run
    with np.load(ROOT/receipt['payload']['path'],allow_pickle=False) as a:arrays={k:a[k] for k in a.files}
    x,phi,data,records=pilot.inputs()  # authenticate existing global fields/source labels
    localization=localize(arrays);local=local_forcing_diagnostic(x,phi,data,receipt['control']['support'])
    F_residual=max(row['F_common_node_difference'] for row in local['times'])
    algebra=max(row['grids'][str(n)][key] for row in local['times'] for n in (401,801)
                for key in ('generator_flux_algebra_residual','generator_tangent_flux_algebra_residual'))
    if F_residual!=0 or algebra>3e-11:raise ArithmeticError('local forcing identity or flux algebra failed')
    return {'schema':'NSC-RETARDED-SPATIAL-DIAGNOSIS-v1','accountable_author':'Douglas Ek',
        'status':'PASS: bounded diagnosis; inherited retarded response remains spatially OPEN',
        'localization':localization,'local_discriminator':local,
        'conditional_causal_support':causal_support_diagnostic(arrays,receipt['control']['support']),
        'residuals':{'common_node_forcing_identity':F_residual,'local_generator_flux_algebra':algebra},
        'verification_tolerances':{'common_node_forcing_identity':0.,'local_generator_flux_algebra':3e-11},
        'decision':{'next_uncertainty':'spatial representation of the same narrow forcing and slow characteristic transport',
            'selected_comparison':'finer spatial evaluation of SAME pulse, with sufficiently accurate reference columns and temporal control',
            'pulse_widening_or_source_alteration_selected':False,'new_resolution_or_field_run_selected':False,
            'existing_three_solve_response_status':receipt['status']},
        'scope':{'field_propagations':0,'radial_solves':0,'old_producer_reruns':0,'saved_data_only_plus_local_matrix_algebra':True,
            'continuum_error_bound':None,'global_causality_proof':False,'physical_C0_change_proved':False,
            'sampled_CAR_response_subtracted':False,'stress_or_counterforce_added':False,
            'physical_profile_or_metric_evolution':False},
        'source_hashes':{p:digest(p) for p in SOURCES},'input_hashes':{p:digest(p) for p in INPUTS},
        'reused_artifacts':[receipt['payload'],*[r['payload'] for r in records]],
        'reproducer':'python3 scripts/derive_nsc_retarded_spatial_diagnosis.py --check'}


def main():
    p=argparse.ArgumentParser();mode=p.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write',action='store_true');mode.add_argument('--check',action='store_true');args=p.parse_args()
    record=make_record();path=ROOT/OUTPUT
    if args.write:path.write_text(json.dumps(record,indent=2,sort_keys=True,allow_nan=False)+'\n')
    elif json.loads(path.read_text())!=record:raise ValueError('spatial diagnosis differs from authenticated replay')
    print(json.dumps({'status':record['status'],'residuals':record['residuals'],
        'local_LF_relative':[r['LF_relative_difference'] for r in record['local_discriminator']['times']],
        'field_L2_relative':record['localization']['spatial_field']['SBP_L2_relative']},indent=2))

if __name__=='__main__':main()
