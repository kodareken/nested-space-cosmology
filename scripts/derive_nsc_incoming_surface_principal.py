#!/usr/bin/env python3
"""Prepare/replay the compatible-family principal matrix, reusing saved cu."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from recursive_horizons.nsc_incoming_surface_principal import (
    local_principal_expressions,local_euler_principal_identity,local_principal_controls,
    reference_shift_principal,linear_weyl_response_identity,combined_principal_identities,
    directed_principal_matrix,
)
from recursive_horizons.nsc_incoming_surface_reference_lapse import reference_lapse_principal

OUTPUT='results/development/nsc-incoming-surface-principal.json'
SOURCES=('src/recursive_horizons/nsc_incoming_surface_principal.py',
         'scripts/derive_nsc_incoming_surface_principal.py','tests/test_nsc_incoming_surface_principal.py',
         'docs/nsc-incoming-surface-principal.md')
RECORDS=('results/development/nsc-incoming-lapse-coefficient.json',
         'results/development/nsc-incoming-local-constraints.json',
         'results/development/nsc-mode-resolved-cauchy-state.json',
         'results/development/nsc-incoming-surface-integrability.json')
OWNERS=('src/recursive_horizons/nsc_incoming_surface_reference_lapse.py',
        'tests/test_nsc_incoming_surface_reference_lapse.py','docs/nsc-incoming-surface-reference-lapse.md',
        'src/recursive_horizons/nsc_incoming_lapse_coefficient.py',
        'src/recursive_horizons/nsc_incoming_cauchy_jets.py',
        'src/recursive_horizons/nsc_incoming_local_constraints.py',
        'src/recursive_horizons/nsc_spatial_reference_symbol.py',
        'src/recursive_horizons/nsc_spherical_local_history.py',
        'src/recursive_horizons/nsc_light_restoration_action.py',
        'src/recursive_horizons/nsc_magnetic_light_reference.py',
        'src/recursive_horizons/nsc_incoming_vacuum_tail_bound.py')
B_OWNER_SHA='afa03c896559f9772409ec4a8fcad81c69414fbae7ae3b520ed7746c62626f3f'
CU_EQUATION='-d*a*L^2/(120*pi*r)*(4*m^2*H_r/M^6+(H_r-H_a)/M^4)'


def digest(path):return sha256((ROOT/path).read_bytes()).hexdigest()
def signature():
    if digest(OWNERS[0])!=B_OWNER_SHA:raise ValueError('frozen independently proved reference B helper changed')
    return {p:digest(p) for p in (*SOURCES,*RECORDS,*OWNERS)}


def inputs():
    records=[json.loads((ROOT/p).read_text()) for p in RECORDS]
    for record in records:
        for key in ('source_hashes','input_hashes'):
            for p,h in record.get(key,{}).items():
                if digest(p)!=h:raise ValueError('principal input changed: '+p)
        if 'payload' in record and digest(record['payload']['path'])!=record['payload']['sha256']:
            raise ValueError('principal input artifact changed')
    if records[0]['identity']['reference_group_coefficient']!=CU_EQUATION:
        raise ValueError('reused exact cu equation changed')
    if not records[0]['coefficient']['strictly_positive']:
        raise ValueError('positive existing cu certificate required')
    return records


def prepare():
    records=inputs();before=signature();cu=records[0]['coefficient']['total_coefficient_interval']
    ledger=records[1]['locked_inputs'];channels=records[2]['channels']
    local=local_principal_expressions();beta=reference_shift_principal();lapse=reference_lapse_principal()
    proof={'local_entries':{name:{k:str(v) for k,v in row.items()} for name,row in local['entries'].items()},
        'local_beta_wave_ratio_residuals':local['beta_wave_ratio_residuals'],
        'full_local_Euler_principal_residuals':local_euler_principal_identity(),
        'reference_beta_transport':linear_weyl_response_identity(),
        'reference_beta_principal':{'kernels':{k:str(v) for k,v in beta['kernels'].items()},
            'closed':{k:str(v) for k,v in beta['closed'].items()},
            'kernel_identity_residuals':beta['kernel_identity_residuals'],
            'moment_identity_residuals':beta['moment_identity_residuals'],
            'beta_wave_ratio_residual':beta['beta_wave_ratio_residual']},
        'reference_lapse_principal':{'closed':str(lapse['closed']),'kernel':str(lapse['kernel']),
            'third_order_odd_kernel':str(lapse['third_order_odd_kernel']),'residuals':lapse['residuals'],
            'trace_G4_pieces':{k:str(v) for k,v in lapse['trace_G4_pieces'].items()},
            'owner_source_sha256':B_OWNER_SHA,'owner_commit':'4680388'},
        'combined':combined_principal_identities(),
        'reused_cu_equation':records[0]['identity']['reference_group_coefficient']}
    controls={str(step):local_principal_controls(ledger,step) for step in (.5,1.)}
    matrix=directed_principal_matrix(channels,ledger,cu)
    if signature()!=before:raise ValueError('principal producer or input changed during preparation')
    payload={'schema':'NSC-INCOMING-SURFACE-PRINCIPAL-ARTIFACT-v1','signature':before,
             'proof':proof,'local_density_controls':controls,'directed_matrix':matrix}
    raw=(json.dumps(payload,indent=2,sort_keys=True,allow_nan=False)+'\n').encode();h=sha256(raw).hexdigest()
    path=ROOT/f'results/development/artifacts/nsc-incoming-surface-principal.{h}.json'
    if path.exists() and path.read_bytes()!=raw:raise ValueError('content-addressed collision')
    if not path.exists():path.write_bytes(raw)
    return path


def make_record(path):
    records=inputs();payload=json.loads(path.read_text())
    if payload['signature']!=signature():raise ValueError('principal producer/input changed')
    matrix=directed_principal_matrix(records[2]['channels'],records[1]['locked_inputs'],records[0]['coefficient']['total_coefficient_interval'])
    if matrix!=payload['directed_matrix']:raise ValueError('directed principal replay differs')
    proof=payload['proof'];zero=[]
    zero+=list(proof['local_beta_wave_ratio_residuals'].values())
    zero+=[x for row in proof['full_local_Euler_principal_residuals'].values() for x in row.values()]
    zero+=[v for k,v in proof['reference_beta_transport'].items() if k.endswith('_residual')]
    zero+=proof['reference_beta_principal']['kernel_identity_residuals']+proof['reference_beta_principal']['moment_identity_residuals']
    zero.append(proof['reference_beta_principal']['beta_wave_ratio_residual'])
    zero+=[v for row in proof['reference_lapse_principal']['residuals'].values() for v in row]
    zero+=[proof['combined']['reference_reused_cu_relation'],proof['combined']['determinant_identity_residual'],proof['combined']['elimination_identity_residual']]
    zero+=list(proof['combined']['local_reused_cu_relations'].values())
    if set(zero)!={'0'}:raise ArithmeticError('exact principal identity failed')
    maximum=0.
    for rows in payload['local_density_controls'].values():
        for name,row in rows.items():
            for key,value in row.items():
                bounds=matrix['components']['local_entries'].get(name,{}).get(key+'_interval',[0.,0.])
                maximum=max(maximum,max(0.,bounds[0]-value,value-bounds[1]))
    first=matrix['determinant_interval'];second=matrix['direct_matrix_determinant_interval']
    overlap=max(0.,first[0]-second[1],second[0]-first[1])
    if maximum>3e-13 or overlap or not matrix['strictly_negative']:
        raise ArithmeticError('principal control or nonzero directed determinant failed')
    return {'schema':'NSC-INCOMING-SURFACE-PRINCIPAL-v1','accountable_author':'Douglas Ek',
        'status':'PASS: nonzero compatible-family principal matrix; no physical initial data or global surface selected',
        'family':{'delta_r':'T*w(z)+T^3*U(z)/6','delta_a':0,
                  'slots':['w','w_z','w_zz','w_zzz','U','U_z'],'evaluation':'all six changes zero at the original baseline'},
        'matrix':matrix,'proof':proof,
        'reused_cu':{'record':RECORDS[0],'interval':records[0]['coefficient']['total_coefficient_interval'],'producer_rerun':False},
        'local_density_controls':payload['local_density_controls'],
        'residuals':{'new_local_density_stencils':maximum,'exact_symbolic_identities':0.,'directed_replay':0.,'independent_determinant_overlap_gap':overlap},
        'verification_tolerances':{'new_local_density_stencils':3e-13,'exact_symbolic_identities':0.,'directed_replay':0.,'independent_determinant_overlap_gap':0.},
        'decision':{'third_spatial_principal_term_cancels':False,'eliminated_w3_coefficient_nonzero':True,
                    'next_connection':'apply the local compatible-function ODE theorem with included finite analytic coefficients; global domain/parent/endpoints remain separate'},
        'scope':{'included_local_and_reference_constraints_only':True,'physical_source_cancels_in_new_coefficients':True,
            'cu_cv_producers_rerun':False,'old_probe_or_reference_quadrature_rerun':False,
            'physical_state_or_parameters_changed':False,'spatial_domain_boundary_length_selected':False,
            'lower_coefficient_numerical_expansion_performed':False,'physical_initial_data_or_root_selected':False,
            'global_Cauchy_surface_proved':False,'parent_preparation_or_endpoints_complete':False,
            'extended_stationarity':'OPEN','metric_evolution':False,'Gamma_rest_assigned':False,'push_or_PDF':False},
        'source_hashes':{p:digest(p) for p in SOURCES},'input_hashes':{p:digest(p) for p in (*RECORDS,*OWNERS)},
        'payload':{'path':str(path.relative_to(ROOT)),'sha256':sha256(path.read_bytes()).hexdigest(),'bytes':path.stat().st_size},
        'reproducer':'python3 scripts/derive_nsc_incoming_surface_principal.py --check'}


def main():
    parser=argparse.ArgumentParser();mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--prepare',action='store_true');mode.add_argument('--check',action='store_true');args=parser.parse_args()
    if args.prepare:
        record=make_record(prepare());(ROOT/OUTPUT).write_text(json.dumps(record,indent=2,sort_keys=True,allow_nan=False)+'\n')
    else:
        old=json.loads((ROOT/OUTPUT).read_text());spec=old['payload']
        if digest(spec['path'])!=spec['sha256']:raise ValueError('principal artifact changed')
        record=make_record(ROOT/spec['path'])
        if record!=old:raise ValueError('principal record differs from authenticated replay')
    print(json.dumps({'status':record['status'],'matrix':record['matrix']['matrix_intervals'],
        'determinant':record['matrix']['determinant_interval'],
        'eliminated_w3':record['matrix']['eliminated_w3_coefficient_interval'],'residuals':record['residuals']},indent=2))


if __name__=='__main__':main()
