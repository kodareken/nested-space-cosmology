#!/usr/bin/env python3
"""Bind the existing local action to the actual shared rho=0 surface."""
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import sympy as sp

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from recursive_horizons.nsc_smooth_seam_variation import smooth_seam_match
from recursive_horizons.nsc_weyl_boundary_jets import density_partials
from recursive_horizons.nsc_lorentzian import geometry
from recursive_horizons.nsc_transmitting_dirac_domain import TransmittingDiracSeamDomain

OUTPUT='results/development/nsc-smooth-seam-variation.json'
SOURCES=('src/recursive_horizons/nsc_smooth_seam_variation.py',
 'scripts/derive_nsc_smooth_seam_variation.py','tests/test_nsc_smooth_seam_variation.py',
 'docs/nsc-smooth-seam-variation.md')
INPUTS=('results/development/nsc-weyl-boundary-jets.json',
 'results/development/nsc-common-time-bulk-split.json',
 'results/development/nsc-transmitting-dirac-domain.json',
 'src/recursive_horizons/nsc_weyl_boundary_jets.py','src/recursive_horizons/nsc_lorentzian.py',
 'src/recursive_horizons/nsc_general_ks_local_history.py','src/recursive_horizons/nsc_boundary_state.py',
 'src/recursive_horizons/nsc_transmitting_dirac_domain.py',
 'docs/nsc-common-time-bulk-split.md','docs/nsc-adm-neck-source-map.md',
 'scripts/derive_nsc_transmitting_boundary_binding.py')
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def rec(p):return json.loads((ROOT/p).read_text())

def reference_jet():
    # Existing geometry() expression; differentiate along its owned future
    # proper child normal. No Einstein equation or geometry solution is run.
    x=sp.Symbol('rho',real=True)
    beta2=3*((1+x*x)*(sp.pi/2-sp.atan(x))-x)
    a=sp.sqrt(beta2-1);r=sp.sqrt(1+x*x)
    DT=lambda f:-a*sp.diff(f,x)
    def at0(f):return sp.simplify(f.subs(x,0))
    symbols={'N':sp.Integer(1),'Nd':sp.Integer(0),'Ndd':sp.Integer(0)}
    for names,value in ((('a','ad','add','addd'),a),(('r','rd','rdd','rddd'),r)):
        for name in names:
            symbols[name]=at0(value);value=DT(value)
    return {k:float(v) for k,v in symbols.items()},{k:str(v) for k,v in symbols.items()}

def plain(value):
    if isinstance(value,np.ndarray):return value.tolist()
    if isinstance(value,dict):return {k:plain(v) for k,v in value.items()}
    if isinstance(value,np.floating):return float(value)
    return value

def calculate():
    old=rec(INPUTS[0]);c=old['locked_inputs'];jet,symbols=reference_jet()
    result=smooth_seam_match(jet,c);p=result['parent']
    independent=density_partials(np.array([[jet[k] for k in ('N','Nd','a','ad','add','r','rd','rdd')]]),c['C_Weyl'])[0]
    result['residuals']['highest_jet_against_owned_density']=float(max(abs(independent[4]-p['Q_a']),abs(independent[7]-p['Q_r'])))
    b,bp,A=map(float,geometry(0.))
    normal=TransmittingDiracSeamDomain(1.,1.,b,1.).normal_geometry()
    result['normal_vectors_PG']={'parent':normal['parent_outward_normal'].tolist(),
                                 'child':normal['child_outward_normal'].tolist()}
    result['residuals']['normal_unit']=normal['normal_unit_residual']
    result['residuals']['normal_tangent']=normal['normal_tangent_residual']
    result['residuals']['reference_geometry_a']=float(abs(jet['a']**2+A))
    result['residuals']['reference_geometry_normal_adot']=float(abs(jet['ad']+b*bp))
    tol=3e-11
    if any(not np.isfinite(v) or v>tol for v in result['residuals'].values()):raise ArithmeticError('physical geometric seam binding failed')
    return {'schema':'NSC-SMOOTH-SEAM-VARIATION-v1','accountable_author':'Douglas Ek',
      'status':'PASS: geometric variation cancels on the shared smooth seam; quantum/full stationarity OPEN',
      'source_hashes':{p:sha(p) for p in SOURCES},'input_hashes':{p:sha(p) for p in INPUTS},
      'locked_inputs':c,'reference_metric_normal_jet':jet,'exact_reference_jet':symbols,
      'domain':{'seam':'one rho=0 surface of the existing smooth global metric',
        'parent':'rho>0, outward normal is the future child normal',
        'child':'rho<0, outward normal is the opposite normal',
        'admissible_local_variations':'common intrinsic metric and matched normal metric jets through order three, expressed with opposite outward normals',
        'justification':'Z1 uses one smooth global metric and smooth compactly supported coefficient variations; a cut does not create two independently variable metrics',
        'scope':'the local reduced action on the reference KS geometry and its smooth shared-jet gluing; not an arbitrary nonsmooth interface',
        'independent_interfacial_field_added':False,'old_control_time_endpoints_identified_with_seam_sides':False},
      'matching':plain(result),'tolerance':tol,
      'legacy_diagnostic':{'Weyl_value':old['preserved_Weyl_value'],'original_eight_coefficients':old['legacy_eight_coefficients'],
          'role':'unchanged discrete initial/final history gradient; not the covector on this common rho=0 surface'},
      'quantum_account':{'full_Gaussian':'the already owned full CTP determinant, with its complementary bulk and initial correlations retained',
        'Schur_factor_added_as_extra_action':False,'Gamma_rest_assigned_zero':False,
        'quantum_or_other_boundary_variation_evaluated':False,
        'warning_scope':'geometric seam cancellation does not remove the bulk expectation-value stress or certify self-sourcing'},
      'gate':{'geometric_seam':'PASS','B1_full':'OPEN: complete physical history/source functional and endpoint embedding',
        'B2_full':'OPEN: remaining same-action quantum/field variation must be identified and evaluated on this same domain',
        'physical_Vc_selected':False,'stationarity':'OPEN','new_stress_assigned':False,
        'metric_timestep':False,'Z3':'OUT OF SCOPE','physical_history_selected':False},
      'comparison':{'fields':'all','float_atol':3e-13,'float_rtol':3e-13,'exact':'hashes, schemas, symbolic identities and scope'},
      'reproducer':'python3 scripts/derive_nsc_smooth_seam_variation.py --check'}

def main():
    result=calculate()
    if '--check' in sys.argv:
        from derive_nsc_transmitting_boundary_binding import compare
        compare(rec(OUTPUT),result)
    else:(ROOT/OUTPUT).write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'parent':result['matching']['parent'],
                      'residuals':result['matching']['residuals']},indent=2))

if __name__=='__main__':main()
