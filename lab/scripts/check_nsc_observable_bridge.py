#!/usr/bin/env python3
"""Verify spinor identification and source-to-observable requirements."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from recursive_horizons.nsc_spinor_bridge import exact_checks, chiral_evolution, weyl_matrices
from check_nsc_plateau_conditions import compare
OUTPUT=ROOT/'results/nsc-7-observable-bridge.json'


def conservation_checks():
    b,d,h=sp.symbols('b d H',positive=True)
    wb,wd,q,jb,jd=sp.symbols('w_b w_d Q_b J_b J_d',real=True)
    bp=-3*(1+wb)*b+(q+jb)/h
    dp=-3*(1+wd)*d+(-q+jd)/h
    rho=b+d;f=d/rho
    fp=sp.diff(f,b)*bp+sp.diff(f,d)*dp
    total=sp.simplify(bp+dp+3*((1+wb)*b+(1+wd)*d))
    expected=3*f*(1-f)*(wb-wd)-q/(h*rho)+jd/(h*rho)-f*(jb+jd)/(h*rho)
    residuals={'internal_transfer_cancels_total':sp.simplify(total-(jb+jd)/h),
               'fraction_with_external_source_partition':sp.simplify(fp-expected)}
    # Applying both usual Einstein equations to a nonconserved complete source
    # is inconsistent. Differentiating the stipulated flat Friedmann constraint
    # alone changes the inferred deceleration by -J/(2Hrho).
    wt=(wb*b+wd*d)/rho
    q_from_constraint=-1-(bp+dp)/(2*rho)
    residuals['sourced_Friedmann_derivative']=sp.simplify(q_from_constraint-(1+3*wt-(jb+jd)/(h*rho))/2)
    dust_residual=fp.subs(wb,0)-(1-f)*(1-2*q_from_constraint.subs(wb,0))+q/(h*rho)+jb/(h*rho)
    residuals['dust_fraction_if_only_sourced_Friedmann_retained']=sp.simplify(dust_residual)
    assert all(v==0 for v in residuals.values())
    controls=[]
    for wd_value in (-1.,-.5):
        rows=[]
        for a in (1.,2.,4.,16.,256.):
            rb=.05*a**-3;rd=.95*a**(-3*(1+wd_value))
            fval=rd/(rb+rd)
            rows.append({'a':a,'H_over_H0':float(np.sqrt(rb+rd)),'dark_fraction':fval,
                         'q_dec':(1+3*wd_value*fval)/2})
        controls.append({'internal_Q':0.,'external_J':0.,'stipulated_w_dark':wd_value,'rows':rows})
    return {'exact_residuals':{k:str(v) for k,v in residuals.items()},
        'closed_room_total':'dot rho+3H(rho+p)=0; internal Q_b cancels',
        'open_room_total':'dot rho+3H(rho+p)=J_b+J_d; external source requires reservoir/boundary gravity accounting',
        'general_fraction':'f_prime=3f(1-f)(w_b-w_d)-epsilon_b+j_d-f*(j_b+j_d)',
        'normalizations':'epsilon_b=Q_b/(H*rho); j_i=J_i/(H*rho)',
        'zero_Q_counterexamples':controls,
        'consequence':'Q=0 alone neither forces a finite plateau nor selects w,H; fixed finite plateau is an additional hypothesis',
        'boundary_dimensions':{'local_stress':'energy/volume','surface_flux_integral':'energy/time',
            'worldtube_flux_integral':'energy','cosmological_Q':'energy/(proper volume*proper time)'},
        'missing_growth_data':['rest-frame pressure response','anisotropic stress','perturbation of energy transfer',
                               'momentum transfer','metric constraints','initial fluctuation spectrum']}


def calculate():
    algebra=exact_checks();cosmology=conservation_checks()
    samples=[chiral_evolution(p,m,h) for p,m,h in ((1.2,.7,1),(1.2,.7,-1),(0.,.7,1),(1.2,0.,1),(1.2,0.,-1))]
    for sample in samples:
        assert abs(sample['positive_energy_chirality']-sample['expected_positive_energy_chirality'])<1e-12
        values=[]
        for row in sample['controls']:
            assert abs(row['pure_left_to_right']-row['analytic_transition'])<1e-12
            values.append(row['positive_energy_right_probability'])
        assert np.ptp(values)<1e-12
    w=weyl_matrices();alpha=np.array(w['alpha'][2],complex);beta=np.array(w['gamma'][0],complex)
    tau1=np.array([[0,1],[1,0]])
    # At |p|=Phi, a spinor-identity sheet link is gapless; the scalar beta link
    # has the Lorentzian mass shell +/-sqrt(p²+Phi²).
    scalar=np.kron(np.eye(2),alpha)+np.kron(tau1,beta)
    identity=np.kron(np.eye(2),alpha)+np.kron(tau1,np.eye(4))
    return {'schema':'nsc-observable-bridge-v1','artifact_id':'NSC-7-OBSERVABLE-BRIDGE',
        'classification':'exact_identification_requirements_and_kinematic_Lorentz_scalar_sheet_bridge_not_particle_or_cosmology_closure',
        'source_hashes':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in (
            'src/recursive_horizons/nsc_spinor_bridge.py','scripts/check_nsc_observable_bridge.py',
            'scripts/check_nsc_plateau_conditions.py')},
        'input_hashes':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in (
            'results/nsc-6-energy-transfer.json','results/nsc-6-vacuum-work.json','results/nsc-5-plateau-conditions.json')},
        'spinor_algebra':algebra,'chiral_evolution':samples,'cosmological_source_accounting':cosmology,
        'scalar_vs_identity_link_control':{'momentum':1.,'Phi':1.,'scalar_link_spectrum':np.linalg.eigvalsh(scalar).tolist(),
                                           'identity_link_spectrum':np.linalg.eigvalsh(identity).tolist()},
        'conventions':{'signature':'+---','basis':'Weyl for4D gamma matrices; separate Dirac basis for spherical expansion',
            'sheet':'tau matrices, independent of spinor chirality',
            'angular_radial':'eta acts on(+kappa,-kappa), rho on(F,G)',
            'spherical_expansion':'psi_kappa=(F Omega_kappa, i G Omega_-kappa)/r; sigma.rhat Omega_kappa=-Omega_-kappa',
            'paired_gamma5':'-eta1 tensor rho2, from full angular expansion',
            'radial_domain_scope':'continuum paired angular operator and compatible domain; no exact naive gamma5 intertwiner claimed for the staggered node/edge regulator',
            'C_scope':'psi^c=i gamma2 psi* is an antilinear c-number map; the quantized charge-conjugation symmetry and physical antiparticle labels are distinct',
            'scalar_sheet_bridge':'L=barPsi(i slash_partial - M)Psi, M=Phi*tau1; a candidate local Lorentz scalar interaction, not derived from current B(E)',
            'quantum_Ward_caution':'classical mass axial identity can acquire gauge/gravitational axial anomaly; no anomaly cancellation claimed'},
        'nonclaims':{'sheet_equals_chirality':False,'sheet_equals_antiparticle':False,
            'radial_angular_gap_is_a_4D_chiral_mass':False,'scalar_sheet_vertex_derived_from_throat':False,
            'Q_alone_predicts_background_or_growth':False,'cosmological_tensions_resolved':False,
            'CPT_of_full_recursive_domain_proved':False,'new_to_world_theorem':False},'terminal':True}


def main():
    parser=argparse.ArgumentParser(description=__doc__);group=parser.add_mutually_exclusive_group()
    group.add_argument('--check',action='store_true');group.add_argument('--output',type=Path)
    args=parser.parse_args()
    if args.output and args.output.exists():raise FileExistsError('refusing to overwrite result')
    record=calculate()
    if args.check:compare(json.loads(OUTPUT.read_text()),record);print('Spinor identification and cosmological source requirements reproduced; every field checked.')
    elif args.output:
        with args.output.open('x') as stream:json.dump(record,stream,indent=2,sort_keys=True,allow_nan=False);stream.write('\n')
        print(args.output)
    else:print(json.dumps(record,indent=2,allow_nan=False))


if __name__=='__main__':main()
