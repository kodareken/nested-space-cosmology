#!/usr/bin/env python3
"""Extract intrinsic Weyl boundary jets and reconstruct the frozen end nodes."""
import hashlib
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from recursive_horizons.nsc_weyl_boundary_jets import nodal_decomposition,FIELDS

OUTPUT='results/development/nsc-weyl-boundary-jets.json'
INPUTS=('results/development/nsc-general-ks-local-history.json',
        'results/development/nsc-landau-cauchy-isometry.json',
        'results/development/nsc-transmitting-boundary-remainder.json',
        'results/development/compact-boundary-action.json',
        'src/recursive_horizons/nsc_general_ks_local_history.py',
        'src/recursive_horizons/nsc_boundary_state.py',
        'src/recursive_horizons/nsc_causal_common.py',
        'src/recursive_horizons/nsc_transmitting_boundary_history.py',
        'scripts/derive_nsc_transmitting_boundary_binding.py')
SOURCES=('src/recursive_horizons/nsc_weyl_boundary_jets.py',
         'scripts/derive_nsc_weyl_boundary_jets.py','tests/test_nsc_weyl_boundary_jets.py',
         'docs/nsc-weyl-boundary-jets.md')
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def rec(p):return json.loads((ROOT/p).read_text())

def calculate():
    local=rec(INPUTS[0]);isometry=rec(INPUTS[1]);binding=rec(INPUTS[2])
    payload=isometry['payload']
    if sha(payload['path'])!=payload['sha256']:raise ValueError('authenticated history changed')
    with np.load(ROOT/payload['path'],allow_pickle=False) as a:
        index=np.linspace(0,len(a['short_time'])-1,33,dtype=int)
        h=SimpleNamespace(time=a['short_time'][index],lapse=a['short_lapse'][index],
            shift=a['short_shift'][index],a_parallel=a['short_a'][index],radius=a['short_r'][index])
    result=nodal_decomposition(h,local['locked_inputs']['C_Weyl'])
    original=np.array([local['diagnostic_evaluation']['action_gradients']['weyl_bulk'][f] for f in FIELDS])
    aliases=('N','beta','q_ADM','r')
    old_endpoints={f:original[i,[0,-1]].tolist() for i,f in enumerate(aliases)}
    if old_endpoints!=binding['locked_local_components']:raise ValueError('the shared eight-coefficient object changed')
    residuals=dict(result['residuals'])
    residuals['full_gradient_vs_locked']=float(np.max(abs(result['gradient']-original)))
    residuals['legacy_eight_coefficient_binding']=0.
    tests={k:v for k,v in residuals.items() if k!='stencil_SBP_defect_norm'}
    if any(not np.isfinite(v) or v>3e-11 for v in tests.values()):raise ArithmeticError('boundary-jet reconstruction failed')
    proper={k:result[k][[0,-1]].tolist() for k in ('proper_boundary_N','proper_boundary_a','proper_boundary_r','proper_boundary_shear')}
    raw={k:{f:result[k][i,[0,-1]].tolist() for i,f in enumerate(aliases)} for k in ('P','Q')}
    parts={k:{f:result[k][i,[0,-1]].tolist() for i,f in enumerate(aliases)} for k in
           ('bulk_component','canonical_boundary_component','stencil_component','reconstructed_gradient')}
    record={'schema':'NSC-WEYL-BOUNDARY-JETS-v1','accountable_author':'Douglas Ek',
      'status':'PASS: intrinsic boundary-jet extraction and discrete reconstruction; physical interface stationarity OPEN',
      'source_hashes':{p:sha(p) for p in SOURCES},'input_hashes':{p:sha(p) for p in INPUTS},
      'history_payload':payload,'locked_inputs':local['locked_inputs'],
      'legacy_endpoint_basis':binding['endpoint_basis'],'legacy_eight_coefficients':old_endpoints,
      'preserved_Weyl_value':binding['preserved_Weyl_value'],
      'boundary_jet_basis':{'raw_values':['N','beta','a=q_ADM','r'],'independent_coordinate_derivatives':['adot','rdot'],
        'orientation':'[boundary form] at final minus initial; P and Q below precede this orientation',
        'proper_normal_rates':'H_a=adot/(N*a), H_r=rdot/(N*r)',
        'proper_shear':'sigma=H_a-H_r','normal_lapse_coefficient':'zero by the same Weyl density identity',
        'not_a_Cauchy_state_isometry':True},
      'raw_boundary_coefficients':raw,'proper_boundary_coefficients':proper,
      'endpoint_nodal_decomposition':parts,
      'full_nodal_decomposition':{k:result[k].tolist() for k in
        ('gradient','bulk_component','canonical_boundary_component','stencil_component','reconstructed_gradient')},
      'identities':{'first_variation':'delta S=bulk+[P_N delta N+P_a delta a+P_r delta r+Q_a delta adot+Q_r delta rdot]',
        'proper_boundary':'[Pbar_a delta a+Pbar_r delta r+Pi_sigma delta(H_a-H_r)]',
        'highest_jet':'a*Q_a+r*Q_r=0; N*P_N+adot*Q_a+rdot*Q_r=0',
        'discrete_gradient':'G=W E+S P+D^T S Q+(B-S)P+D^T(B-S)Q; B=WD+D^TW',
        'reason_for_stencil_term':'the existing edge_order=2 gradient with trapezoid weights is not an exact endpoint-only summation-by-parts operator'},
      'residuals':residuals,'tolerance':3e-11,
      'variational_domain':{
        'fixed_intrinsic_metric_and_normal_rates':'canonical Weyl boundary form vanishes on these variations; not a self-sourcing claim',
        'fixed_intrinsic_metric_only':'free normal shear still has coefficient Pi_sigma and cannot be discarded',
        'free_or_glued_interface':'requires the actual matching of intrinsic metric, normal shear and their oriented conjugates',
        'physical_matching_data_supplied':False,
        'time_end_nodes_are_not_the_two_sides_of_rho_zero':True},
      'remaining_action_ownership':{
        'Gaussian_complement':'GaussianBoundaryState.overlap_schur factorizes the same Gamma_G; it is not a second added action',
        'geometric_terms':'Euler/boxR endpoints already belong to GeneralKSLocalInducedHistory and are not introduced again',
        'Gamma_rest':'the declared remainder excludes geometric pi but its full transmitting functional and permitted boundary variations are not specified by these records',
        'Gamma_rest_assigned_zero':False,'arbitrary_B_introduced':False},
      'gate':{'B1_full':'OPEN: physical history/endpoint Cauchy embedding remains distinct from nodal bookkeeping',
        'B2':'OPEN: intrinsic boundary-jet matching and remaining same-action derivative not evaluated',
        'physical_two_sided_mismatch':None,'stationarity':'OPEN','history_selected':False,
        'metric_timestep':False,'new_physical_term':False,'stress_computed':False,'Z3':'OUT OF SCOPE'},
      'comparison':{'fields':'all','float_atol':3e-13,'float_rtol':3e-13,'exact':'hashes, original eight coefficients, schemas and scope'},
      'reproducer':'python3 scripts/derive_nsc_weyl_boundary_jets.py --check'}
    return record

def main():
    actual=calculate()
    if '--check' in sys.argv:
        from derive_nsc_transmitting_boundary_binding import compare
        compare(rec(OUTPUT),actual)
    else:(ROOT/OUTPUT).write_text(json.dumps(actual,sort_keys=True,indent=2)+'\n')
    print(json.dumps({'status':actual['status'],'residuals':actual['residuals'],
         'proper_boundary_coefficients':actual['proper_boundary_coefficients']},indent=2))

if __name__=='__main__':main()
