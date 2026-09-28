#!/usr/bin/env python3
"""Authenticate exact compatible-surface polynomials and replay their intervals."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from recursive_horizons.nsc_incoming_surface_coefficients import polynomial_proof,directed_coefficients,new_algebra_controls
OUTPUT='results/development/nsc-incoming-surface-coefficients.json'
SOURCES=('src/recursive_horizons/nsc_incoming_surface_coefficients.py','scripts/derive_nsc_incoming_surface_coefficients.py',
         'tests/test_nsc_incoming_surface_coefficients.py','docs/nsc-incoming-surface-coefficients.md')
RECORDS=tuple('results/development/nsc-incoming-'+x+'.json' for x in ('shift-coefficient','lapse-coefficient','surface-principal','constraint-germ','local-constraints','joint-constraints'))+('results/development/nsc-mode-resolved-cauchy-state.json',)
OWNERS=tuple('src/recursive_horizons/'+x+'.py' for x in ('nsc_incoming_surface_quadratic_lapse','nsc_incoming_shift_coefficient',
    'nsc_incoming_lapse_coefficient','nsc_incoming_cauchy_jets','nsc_spatial_reference_symbol','nsc_incoming_local_constraints',
    'nsc_light_restoration_action','nsc_spherical_local_history','nsc_magnetic_light_reference','nsc_incoming_vacuum_tail_bound'))+(
    'docs/nsc-incoming-surface-quadratic-lapse.md','tests/test_nsc_incoming_surface_quadratic_lapse.py','docs/nsc-incoming-compatible-surface.md')
C_SHA='98a1a52275bfb198d671d493541548774589f09db57f9fbeabcc3a7ba635b529'


def digest(path):return sha256((ROOT/path).read_bytes()).hexdigest()
def signature():
    if digest(OWNERS[0])!=C_SHA:raise ValueError('independent C reference owner changed')
    return {p:digest(p) for p in (*SOURCES,*RECORDS,*OWNERS)}
def inputs():
    records=[json.loads((ROOT/p).read_text()) for p in RECORDS]
    for record in records:
        for key in ('source_hashes','input_hashes'):
            for p,h in record.get(key,{}).items():
                if digest(p)!=h:raise ValueError('surface polynomial input changed: '+p)
        if 'payload' in record and digest(record['payload']['path'])!=record['payload']['sha256']:
            raise ValueError('surface polynomial archived payload changed')
    return records

def stored_controls(directed,records):
    shift,cu,principal,germ,local,joint,inventory=records
    midpoint=lambda bounds:sum(bounds)/2
    C=midpoint(directed['C_interval']);A=midpoint(cu['coefficient']['total_coefficient_interval'])
    cv=midpoint(shift['coefficient']['total_coefficient_interval']);Cref=midpoint(directed['C_reference_interval'])
    errors={'v_lapse_action':abs(C*.01**2-germ['v_probe_response'][0]),
            'mixed_lapse_action':abs(A*.02+C*(-.005)**2-germ['mixed_control_response'][0]),
            'v_shift_action':abs(cv*.01-germ['v_probe_response'][1]),
            'mixed_shift_action':abs(cv*(-.005)-germ['mixed_control_response'][1]),
            'reference_v_lapse_action':abs(Cref*.01**2-joint['reference_controls']['n24']['reference_action_gradient_change'][0])}
    localerrors={n:abs(midpoint(bounds)*.01**2-(joint['local_probe']['channel_action_gradients'][n][0]-local['baseline']['channel_action_gradients'][n][0])) for n,bounds in directed['local_C_channels'].items()}
    return {'action_errors':errors,'local_C_channel_action_errors':localerrors,
            'maximum':max(*errors.values(),*localerrors.values()),'tolerance':3e-11,
            'old_coefficient_difference_is_not_an_action_error':abs(C-germ['numerical_coefficients']['c_v2']),
            'new_probe_runs':0,'control_parameters':{'v':.01,'mixed_u':.02,'mixed_v':-.005}}

def prepare():
    records=inputs();before=signature();proof=polynomial_proof(records[0])
    ledger=records[4]['locked_inputs'];channels=records[6]['channels'];cv=records[0]['coefficient']['total_coefficient_interval']
    directed=directed_coefficients(proof,channels,ledger,cv)
    fine=directed_coefficients(proof,channels,ledger,cv,precision=70)
    algebra=new_algebra_controls(ledger);controls=stored_controls(directed,records)
    if signature()!=before:raise ValueError('surface coefficient producer changed during preparation')
    payload={'schema':'NSC-INCOMING-SURFACE-COEFFICIENTS-ARTIFACT-v1','signature':before,'proof':proof,
        'directed':directed,'directed70':fine,'new_algebra_controls':algebra,'stored_controls':controls}
    raw=(json.dumps(payload,indent=2,sort_keys=True,allow_nan=False)+'\n').encode();h=sha256(raw).hexdigest()
    path=ROOT/f'results/development/artifacts/nsc-incoming-surface-coefficients.{h}.json'
    if path.exists() and path.read_bytes()!=raw:raise ValueError('content addressed collision')
    if not path.exists():path.write_bytes(raw)
    return path

def make_record(path):
    records=inputs();payload=json.loads(path.read_text())
    if payload['signature']!=signature():raise ValueError('surface coefficient producer/input changed')
    proof=payload['proof'];ledger=records[4]['locked_inputs'];channels=records[6]['channels'];cv=records[0]['coefficient']['total_coefficient_interval']
    directed=directed_coefficients(proof,channels,ledger,cv)
    if directed!=payload['directed']:raise ValueError('surface coefficient directed replay mismatch')
    fine=directed_coefficients(proof,channels,ledger,cv,precision=70)
    if fine!=payload['directed70']:raise ValueError('surface coefficient refined replay mismatch')
    if set(proof['identities'].values())!={'0'}:raise ArithmeticError('surface coefficient symbolic identity failed')
    controls=stored_controls(directed,records)
    if controls!=payload['stored_controls']:raise ValueError('stored control replay mismatch')
    # Both precision runs enclose the same input hull; measure intersection gaps.
    pairs=[(directed['C_interval'],fine['C_interval'])]
    pairs += [(directed['polynomial_intervals'][f][str(j)],fine['polynomial_intervals'][f][str(j)]) for f,n in [('D',4),('F',2)] for j in range(1,n+1)]
    gap=max(max(0.,x[0]-y[1],y[0]-x[1]) for x,y in pairs)
    algebra=payload['new_algebra_controls']
    residuals={'exact_identities':0.,'new_w_point_energy':algebra['new_w_point_energy_maximum'],
        'local_density_partial_derivatives':algebra['local_density_maximum'],'stored_action_controls':controls['maximum'],
        'directed_precision_overlap_gap':gap,'directed_replay':0.}
    tolerances={'exact_identities':0.,'new_w_point_energy':3e-13,'local_density_partial_derivatives':3e-13,
                'stored_action_controls':3e-11,'directed_precision_overlap_gap':0.,'directed_replay':0.}
    if any(value>tolerances[k] for k,value in residuals.items()):raise ArithmeticError('surface coefficient verification failed: '+str(residuals))
    return {'schema':'NSC-INCOMING-SURFACE-COEFFICIENTS-v1','accountable_author':'Douglas Ek',
        'status':'PASS: fixed included compatible-surface polynomial coefficients',
        'operator':{'E_N':'A(w)*U+B(w)*w_zz+C*w_z²+S_N+sum(Dj*w^j,j=1..4)',
                    'E_beta':'d*U_z+e*w_zzz+(cv+F1*w+F2*w²)*w_z+S_beta',
                    'source_constants':'S_N and S_beta remain exact formal finite baseline constraints; no numerical source substitution',
                    'primitive':'G(w)=cv*w+F1*w²/2+F2*w³/3; E_beta=d_z[d*U+e*w_zz+G(w)]+S_beta',
                    'fixed_raw_derivatives':['r_TT','r_TTT','all a derivatives'],'family':'delta r=T*w(z)+T³*U(z)/6, delta a=0'},
        'coefficients':directed,'exact_polynomials':proof,'stored_controls':controls,'new_algebra_controls':algebra,
        'residuals':residuals,'verification_tolerances':tolerances,
        'reused':{'cv_interval':cv,'cv_producer_rerun':False,'cu_producer_rerun':False,
            'C_reference_owner_commit':'b74f5fc','C_reference_source_sha256':C_SHA,
            'principal_coefficients_record':RECORDS[2]},
        'scope':{'included_N_beta_operators':True,'pressure_equations_claimed':False,'source_state_scale_action_parameters_changed':False,
            'old_reference_quadrature_or_probe_campaign_rerun':False,'source_numerical_accuracy_inferred':False,
            'physical_normal_data_or_boundary_selected':False,'w_profile_integrated':False,'metric_evolution':False,'Gamma_rest_assigned':False},
        'source_hashes':{p:digest(p) for p in SOURCES},'input_hashes':{p:digest(p) for p in (*RECORDS,*OWNERS)},
        'payload':{'path':str(path.relative_to(ROOT)),'sha256':sha256(path.read_bytes()).hexdigest(),'bytes':path.stat().st_size},
        'reproducer':'python3 scripts/derive_nsc_incoming_surface_coefficients.py --check'}

def main():
    parser=argparse.ArgumentParser();mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--prepare',action='store_true');mode.add_argument('--check',action='store_true');args=parser.parse_args()
    if args.prepare:
        record=make_record(prepare());(ROOT/OUTPUT).write_text(json.dumps(record,indent=2,sort_keys=True,allow_nan=False)+'\n')
    else:
        old=json.loads((ROOT/OUTPUT).read_text());p=old['payload']
        if digest(p['path'])!=p['sha256']:raise ValueError('surface coefficient artifact changed')
        record=make_record(ROOT/p['path'])
        if record!=old:raise ValueError('surface coefficient record differs')
    print(json.dumps({'status':record['status'],'polynomials':record['coefficients']['polynomial_intervals'],'C':record['coefficients']['C_interval'],'residuals':record['residuals']},indent=2))

if __name__=='__main__':main()
