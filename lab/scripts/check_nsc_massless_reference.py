#!/usr/bin/env python3
"""Bind an absolute massless reference to the unchanged finite determinant."""
import argparse
import json
from math import pi
from pathlib import Path
import sys

import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
from recursive_horizons.nsc_massless_reference import (
    reference_extension,canonical_reference_density,light_conversion,matched_reference_budget,
)
from check_nsc_compact_boundary_action import authenticated_record
from check_nsc_compact_casimir import compare
from check_nsc_vacuum_charge_matching import hashes

OUTPUT=ROOT/"results/development/massless-reference.json"
SOURCES=("scripts/check_nsc_massless_reference.py","src/recursive_horizons/nsc_massless_reference.py",
         "docs/nsc-massless-reference.md","scripts/check_nsc_compact_boundary_action.py",
         "scripts/check_nsc_compact_casimir.py","scripts/check_nsc_vacuum_charge_matching.py")
INPUTS=("results/development/child-state.json","results/development/spectral-endpoint.json",
        "results/development/compact-matching.json","results/nsc-9-covariant-source.json")


def calculate():
    child,endpoint,matching,_=[authenticated_record(p) for p in INPUTS]
    exact={}
    def zero(name,value):
        assert sp.simplify(value)==0,name
        exact[name]="0"
    h,a=sp.symbols("H a",positive=True)
    R=-12*h*h;ric=-3*h*h
    H3=ric*ric-sp.Rational(2,3)*R*ric-(4*ric*ric)/2+R*R/4
    H1=-R*R/2+2*R*ric
    zero("imported_H3_on_limiting_geometry",H3-3*h**4)
    zero("finite_R2_variation_on_limiting_geometry",H1)
    zero("Dirac_Euler_to_reference_density",sp.Rational(11,360)*H3/(8*sp.pi**2)-11*h**4/(960*sp.pi**2))
    volume=8*sp.pi**2*a**4/3
    zero("reference_trace_to_sphere_log_derivative",4*volume*11/(960*sp.pi**2*a**4)-sp.Rational(11,90))

    # General conversion is an action identity, not an adjustable source.
    F,L,Ren=sp.symbols("Gamma5 Gamma_light_nu Gamma_light_ren")
    zero("same_functional_partition",Ren+(F-L)+(L-Ren)-F)
    nu,mu,t,A0,A2,A4=sp.symbols("nu mu t A0 A2 A4",positive=True)
    divergent=(A0*nu**4/2+A2*nu**2+A4*sp.log(nu**2/mu**2))/(32*sp.pi**2)
    heat_local=(A0*nu**4+A2*nu**2+A4)/(16*sp.pi**2)
    zero("proper_time_divergent_primitive",nu*sp.diff(divergent,nu)-heat_local)
    sphere_heat=(sp.Rational(2,3)*(a**4/t**2-a*a/t)+sp.Rational(11,90))
    zero("sphere_scale_and_heat_derivative",a*sp.diff(sphere_heat,a)+2*t*sp.diff(sphere_heat,t))

    # Apply the imported conformal-stress tensors to the recorded geometry.
    # No full Christoffel/Riemann derivation or old generator is repeated.
    u=sp.Symbol("u",positive=True)
    W=3*(sp.pi-sp.atan(u))+(3*u-u*u)/(1+u*u)
    radius=sp.sqrt(1+u*u)/u;lapse=-(1+u*u)*W/u**2
    prime=lambda v:u*u*sp.diff(v,u)
    lp=prime(lapse);lpp=prime(lp);rp=-1/sp.sqrt(1+u*u)
    sec=(lpp/2,lp*rp/(2*radius),lapse/radius**4+lp*rp/(2*radius),(1-lapse*rp*rp)/radius**2)
    aa,bb,cc,dd=sec
    ricci=[sp.series(v,u,0,4).removeO().expand() for v in (aa+2*bb,aa+2*cc,bb+cc-dd,bb+cc-dd)]
    scalar=sp.expand(sum(ricci))
    zero("actual_R_series",scalar-(-36*sp.pi+18*sp.pi*u*u))
    zero("actual_Weyl_amplitude_order",sp.series(aa-bb-cc-dd,u,0,4).removeO()-6*u**3)
    scalar_p=prime(scalar);scalar_pp=prime(scalar_p)
    hessian=[-lp*scalar_p/2,-lapse*scalar_pp-lp*scalar_p/2,
             -lapse*rp/radius*scalar_p,-lapse*rp/radius*scalar_p]
    hessian=[sp.series(v,u,0,3).removeO().expand() for v in hessian]
    box=sum(hessian);ricci2=sum(v*v for v in ricci)
    third=[sp.series(v*v-sp.Rational(2,3)*scalar*v-ricci2/2+scalar*scalar/4,u,0,3).removeO() for v in ricci]
    first=[sp.series(-scalar*scalar/2+2*scalar*v+2*box-2*hh,u,0,3).removeO() for v,hh in zip(ricci,hessian)]
    # In the inherited heat scheme alpha=-1/20, beta=11/360, gamma=-1/30.
    tensor=[sp.simplify((11*t3-t1)/(2880*sp.pi**2)) for t3,t1 in zip(third,first)]
    wanted=[sp.Rational(33,320)-u*u/32,sp.Rational(33,320)-3*u*u/32,
            sp.Rational(33,320)-u*u/32,sp.Rational(33,320)-u*u/32]
    for name,v,target in zip(("z","T","theta","phi"),tensor,wanted):zero("reference_mixed_"+name,v-target)
    zero("reference_radial_null_coefficient",tensor[1]-tensor[0]+u*u/16)
    euler=sp.series(8*(-aa*dd+2*bb*cc),u,0,3).removeO()
    anomaly=(sp.Rational(11,360)*euler-sp.Rational(1,30)*box)/(16*sp.pi**2)
    zero("actual_tail_trace_matches_inherited_heat_scheme",sum(tensor)-anomaly)
    zero("actual_tail_work_through_u_squared",-sp.sqrt(3*sp.pi)*u*sp.diff(tensor[1],u)
         +3*sp.sqrt(3*sp.pi)*(tensor[1]-tensor[0]))

    # Euler normalization is consumed from the authenticated coefficient map.
    coeff=next(c for c in matching["matched_coefficients"] if c["matching_cutoff"]==1.)
    imported_heat_constant=128*pi*pi*coeff["C_Euler_Dirac"]/coeff["matched_weight_at_zero"]
    assert abs(imported_heat_constant-11/90)<2e-15
    full=endpoint["endpoint_comparisons"][0]["full"]
    assert child["unchanged_endpoint_input"]["Dirac_log_radius_derivative"]==full["d_log_radius"]
    budgets=[matched_reference_budget(full["isotropic_action_response"],full["radius"],nu) for nu in (.5,1.,1.4)]
    for row in budgets:assert abs(row["matching_residual"])<3e-16
    heat_controls=[light_conversion(full["radius"],x/full["radius"]) for x in (2.,4.,8.,16.)]
    remainder=[abs(r["finite_remainder_log_derivative"]) for r in heat_controls]
    assert all(a>b for a,b in zip(remainder,remainder[1:]))
    assert remainder[-1]<7e-5
    extension=[{"u":u,"W_extended":reference_extension(u),"width":.1} for u in (-.2,-.1,-.075,-.05,0.,.1)]
    assert all(r["W_extended"]>0 for r in extension)
    assert extension[0]["W_extended"]==extension[1]["W_extended"]==3*pi
    return {"schema":"NSC-MASSLESS-REFERENCE-v1","status":"massless reference asymptotics and exact regulator-allocation map",
       "source_hashes":hashes(SOURCES),"input_hashes":hashes(INPUTS),"exact_residuals":exact,
       "reference_specification":{"field":"one free massless four-component Dirac field; untwisted full angular spectrum, no link mass",
         "construction":"smooth extension of g/r² across u=0 to an auxiliary ultrastatic R_eta x R_z x S² region; its ground covariance propagated back by the conformal Dirac operator, then psi=r^-3/2 psi_tilde",
         "width":.1,"width_role":"auxiliary reference choice, not a derived physical parameter or a stress-fitting knob",
         "ground_covariance":"C_out=1_(H_ultrastatic<0); omega_k,kappa²=k²/(3*pi)+kappa², |kappa|>=1; no zero eigenvalue",
         "time_direction":"physical future decreases u; in the constant region use conformal proper eta and extend that ultrastatic time freely",
         "transport":"C_ref=U_barD C_out U_barD^dagger; all channels retained in the reference definition",
         "parent_covariance_matched":False,"parent_flux_compatibility_verified":False,
         "finite_radius_absolute_stress_evaluated":False,"extension_probes":extension},
       "renormalization":{"signature":"+---; physical density is T^T_T",
         "reference_scheme":"covariant proper-time subtraction matching the inherited heat coefficients; no extra dimensionful vacuum/Einstein term allocated to renormalized light",
         "anomaly_coefficients":{"C2":"-1/20","Euler":"11/360","boxR":"-1/30"},
         "finite_R2_scheme":"must be matched when converting a Hadamard implementation; the u² coefficient below is in this heat scheme",
         "normalization_mu":"held fixed and distinct from matching nu and physical Lambda",
         "complete_physical_finite_coefficients_set_to_zero":False},
       "actual_tail_reference":{"mixed_order":"z,T,theta,phi","Ricci_series":[str(v) for v in ricci],
         "R_L":"-36*pi+18*pi*u²+O(u4)","physical_Weyl_amplitude":"6*u³+O(u4)",
         "barred_Weyl_amplitude":"O(u), because C_orth[g]=r^-2 C_orth[g_bar]",
         "conformal_Weyl_term":"r^-4 * bar(nabla nabla)(log(r)*bar(C)) = O(u3); it must not be discarded at finite u",
         "barred_state_and_curvature_terms":"bounded at u=0, hence O(u4) after physical rescaling",
         "reference_mixed_series":[str(v)+" + O(u³)" for v in tensor],
         "rho":"33/320-3*u²/32+O(u³)","p_parallel_and_sphere":"-33/320+u²/32+O(u³)",
         "radial_null":"-u²/16+O(u³), for this massless reference and heat scheme",
         "limiting_density":canonical_reference_density(full["radius"]),
         "finite_neck_or_full_NSC_source":False},
       "conversion_functional":{"definition":"C_nu,mu=Gamma_light,nu-Gamma_light,ren(mu)",
         "heat_representation":"I_div(nu,mu)-(1/2)int_0^(nu^-2) ds/s [K_light(s)-(4*pi*s)^-2(A0+A2*s+A4*s²)]",
         "I_div":"[A0*nu^4/2+A2*nu²+A4*log(nu²/mu²)]/(32*pi²)",
         "partition":"Gamma5=Gamma_light,ren + [Gamma_H,nu+C_nu,mu]",
         "general_metric_variation":"delta C=delta I_div-(1/2)int ds/s delta(K-K_asym); finite nonlocal remainder retained",
         "spacetime_zero_modes_and_boundaries":"the representation needs the same prescribed IR and boundary treatment; S4 control has no Dirac zero mode",
         "causal_conversion_on_actual_child_computed":False},
       "invariant_endpoint_budgets":budgets,"imported_light_heat_constant":imported_heat_constant,
       "cutoff_removal_controls":heat_controls,
       "required_next_source":{"full_angular_covariance_on_actual_child":True,"parent_reference_correlations_and_flux":True,
         "spacetime_CTP_conversion_kernel":True,"massive_compact_and_interaction_contributions":True,
         "actual_neck_independent_metric_residuals":True},
       "scope":{"sphere_or_compact_generator_rerun":False,"instantaneous_spatial_cutoff_substituted":False,
         "retardation_rule_applied_without_state_hypotheses":False,"normalization_conversion_changes_Gamma5":False,
         "negative_reference_null_stress_proves_full_self_sourcing":False,"full_angular_finite_child_stress_computed":False,
         "source_budget_at_neck_closed":False,"cosmological_Q_or_dark_fraction_derived":False},
       "comparison":{"fields":"all","float_atol":3e-9,"float_rtol":3e-8,"exact":"symbolic formulas, identities, keys, source/input hashes and scope","exceptions":[]}}


def main():
    parser=argparse.ArgumentParser(description=__doc__);group=parser.add_mutually_exclusive_group()
    group.add_argument("--check",action="store_true");group.add_argument("--output",type=Path)
    args=parser.parse_args()
    if args.output and args.output.exists():raise FileExistsError("refusing to overwrite recorded evidence")
    expected=json.loads(OUTPUT.read_text()) if args.check else None
    if args.check:compare(expected["source_hashes"],hashes(SOURCES));compare(expected["input_hashes"],hashes(INPUTS))
    result=calculate()
    if args.check:compare(expected,result);print("massless reference: all fields reproduced; prior endpoint and compact generators were not run")
    elif args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        with args.output.open("x") as stream:json.dump(result,stream,indent=2,sort_keys=True,allow_nan=False);stream.write("\n")
        print(args.output)
    else:print(json.dumps(result,indent=2,allow_nan=False))


if __name__=="__main__":main()
