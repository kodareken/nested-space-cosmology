#!/usr/bin/env python3
"""Bind the known curvature coefficients to a first-order matter/metric EFT."""
import argparse
import json
from pathlib import Path
import sys

import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
from recursive_horizons.nsc_curvature_eft import (
    tensor_trace, tensor_square, inverse_metric_shift, stress_contact,
    spherical_contact, corrected_constraints,
)
from recursive_horizons.nsc_spherical_action import canonical_constraints, potentials
from check_nsc_compact_boundary_action import authenticated_record
from check_nsc_vacuum_charge_matching import compare, hashes

OUTPUT=ROOT/"results/development/curvature-eft.json"
SOURCES=("scripts/check_nsc_curvature_eft.py","src/recursive_horizons/nsc_curvature_eft.py",
         "docs/nsc-curvature-eft.md","scripts/check_nsc_compact_boundary_action.py",
         "scripts/check_nsc_vacuum_charge_matching.py")
INPUTS=("results/development/spherical-action.json","results/development/compact-matching.json",
        "results/development/gauge-source.json")


def calculate():
    for path in INPUTS:authenticated_record(path)
    exact={}
    def zero(name,value):
        if isinstance(value,sp.MatrixBase):
            assert value.applyfunc(sp.simplify)==sp.zeros(*value.shape),name
        else:assert sp.factor(value)==0,name
        exact[name]="0"

    eta=sp.diag(1,-1,-1,-1)
    A,CF,X,mag=sp.symbols("A C_F X q_mag",positive=True)
    cW,cR,V,eps=sp.symbols("cW cR V eps",real=True)
    def symmetric(prefix):
        fields=sp.symbols(prefix+"0:10",real=True)
        matrix=sp.zeros(4);k=0
        for i in range(4):
            for j in range(i,4):matrix[i,j]=matrix[j,i]=fields[k];k+=1
        return matrix
    E,T=symmetric("Einstein_residual_"),symmetric("stress_")
    G=E+T/(2*A);R=-tensor_trace(G,eta);Ric=G+R*eta/2
    shift=inverse_metric_shift(Ric,T,eta,A,cW,cR)
    S1=-2*cW*(tensor_square(Ric,eta)-R*R/3)-cR*R*R
    deltaS0=-A*sp.trace(E*shift)
    zero("off_shell_order_reduction_identity",S1+deltaS0-stress_contact(T,eta,A,cW,cR))
    Riem,Ric2,scalar,Euler=sp.symbols("Riemann2 Ricci2 R Euler4")
    zero("fourD_Weyl_Euler_identity",Riem-2*Ric2+scalar**2/3-(Riem-4*Ric2+scalar**2)-2*(Ric2-scalar**2/3))

    # Matter and vacuum are included together before squaring the stress.
    area,u,v,rho=sp.symbols("area u v rho_EM",real=True,nonzero=True)
    K=u*u-v*v;kin=(u*u+v*v)/(2*area)
    Ts=sp.Matrix([[kin+rho+V,u*v/area,0,0],[u*v/area,kin-rho-V,0,0],
                  [0,0,rho-V,0],[0,0,0,rho-V]])
    zero("spherical_total_trace_includes_vacuum",tensor_trace(Ts,eta)-4*V)
    targetQ=K*K/(2*area**2)+4*rho*rho-4*V*V/3
    zero("spherical_stress_contraction",tensor_square(Ts,eta)-tensor_trace(Ts,eta)**2/3-targetQ)
    chi=sp.Symbol("chi",real=True)
    pars=potentials(X,A,V,CF,mag)
    H=X/(2*A);r2=X/(8*sp.pi*A)
    rhoEM=CF*mag**2/(2*r2**2)+pars["mu_squared"]*chi**2/(2*H)
    spherical=spherical_contact(X,chi,K,A,V,CF,mag,cW,cR)
    projected=area*stress_contact(Ts,eta,A,cW,cR)
    zero("contact_to_existing_spherical_source",projected.subs({area:H,rho:rhoEM})-spherical)

    # Leading electric elimination commutes with the first-order insertion.
    electric,e2,field_shift=sp.symbols("electric e2_squared field_shift",real=True,nonzero=True)
    E0=e2*sp.sqrt(mag/sp.pi)*chi
    L0=electric**2/(2*e2)-sp.sqrt(mag/sp.pi)*electric*chi
    zero("leading_electric_elimination_at_first_order",sp.expand(L0.subs(electric,E0+eps*field_shift)).coeff(eps,1))

    # A vacuum shift in redefined variables is not a physical Lambda shift.
    lam=V/(2*A);T_vac=V*eta;Ric_vac=-lam*eta
    map_vac=inverse_metric_shift(Ric_vac,T_vac,eta,A,cW,cR)
    shift_vac=sp.factor(map_vac[0,0])
    contact_vac=sp.factor(stress_contact(T_vac,eta,A,cW,cR))
    lam_new=(V-eps*contact_vac)/(2*A)
    zero("constant_curvature_observable_pullback",sp.expand((1+eps*shift_vac)*lam_new-lam).coeff(eps,1))

    # New canonical correction is the perturbative Legendre transform.
    p=(X,*sp.symbols("Xplus Xminus",real=True))
    q=sp.symbols("omega emin eplus",nonzero=True);dp=sp.symbols("X_x Xplus_x Xminus_x")
    dq=sp.symbols("omega_x emin_x eplus_x")
    fx,fxx,pi,pix=sp.symbols("chi_x chi_xx pi_chi pi_chi_x",real=True)
    b2,b3,velocity=sp.symbols("bar2 bar3 velocity",real=True)
    volume=b2*q[2]-b3*q[1]
    g00=-2*q[1]*q[2]/volume**2;g01=(b2*q[2]+b3*q[1])/volume**2;g11=-2*b2*b3/volume**2
    velocity0=(pi-volume*g01*fx)/(volume*g00)
    kinetic=(g00*velocity**2+2*g01*velocity*fx+g11*fx*fx)
    K0=(fx*fx-pi*pi)/(2*q[1]*q[2])
    zero("leading_kinetic_in_canonical_variables",kinetic.subs(velocity,velocity0)-K0)
    G0=canonical_constraints(p,q,dp,pars["U"],pars["V"],chi,fx,pi,pars["mu_squared"])
    corrected=corrected_constraints(p,q,dp,chi,fx,pi,A,V,CF,mag,cW,cR,eps)
    C1=sp.simplify((corrected-G0)/eps)
    L1=volume*spherical_contact(X,chi,kinetic,A,V,CF,mag,cW,cR)
    zero("contact_Legendre_to_corrected_constraints",-L1.subs(velocity,velocity0)+b2*C1[1]+b3*C1[2])

    # Verify the constraint algebra for any d(X) K²+F(X,chi), then bind
    # the actual contact to that class. This avoids assuming its closure.
    d=sp.Function("d")(X);f=sp.Function("F")(X,chi)
    F=d*K0*K0+f
    actualF=spherical_contact(X,chi,K0,A,V,CF,mag,cW,cR)
    actual_d=-cW/(2*A*X)
    actual_f=spherical_contact(X,chi,0,A,V,CF,mag,cW,cR)
    zero("actual_contact_in_verified_constraint_class",actualF-actual_d*K0*K0-actual_f)
    Gc=G0+eps*sp.Matrix([0,q[2]*F,-q[1]*F])
    a,ax,b,bx=sp.symbols("a a_x b b_x")
    variables=(*p,*q,chi,fx,pi);jets=(*dp,*dq,fx,fxx,pix)
    def derivative(value):return sp.expand(sum(sp.diff(value,z)*dz for z,dz in zip(variables,jets)))
    def functional(i,z,zx):
        value=Gc[i]
        return ([z*sp.diff(value,t) for t in q],
                [z*sp.diff(value,t)-zx*int(k==i) for k,t in enumerate(p)],
                z*sp.diff(value,chi)-zx*sp.diff(value,fx)-z*derivative(sp.diff(value,fx)),
                z*sp.diff(value,pi))
    structures={}
    for i,j in ((0,1),(0,2),(1,2)):
        aq,ap,af,am=functional(i,a,ax);bq,bp,bf,bm=functional(j,b,bx)
        B=sp.expand(sum(x*y-z*w for x,y,z,w in zip(aq,bp,ap,bq))+af*bm-am*bf)
        density=sp.expand(B.coeff(a).coeff(b)-derivative(B.coeff(ax).coeff(b)))
        delta=sp.expand(B.coeff(a).coeff(bx)-B.coeff(ax).coeff(b))
        structure=[]
        for order in (0,1):
            zero(f"constraint_{i+1}{j+1}_delta_order_{order}",delta.coeff(eps,order))
            structure.append(sp.Matrix([sp.factor(sp.diff(density.coeff(eps,order),z)) for z in dp]))
        residual=sp.expand(density-(structure[0]+eps*structure[1]).dot(Gc))
        for order in (0,1):zero(f"constraint_{i+1}{j+1}_closure_order_{order}",residual.coeff(eps,order))
        structures[f"{i+1}{j+1}"]={"leading":[str(z) for z in structure[0]],"first_order":[str(z) for z in structure[1]]}

    return {
        "schema":"NSC-CURVATURE-EFT-v1",
        "status":"local curvature-squared corrections mapped to source contacts and constraints at first EFT order",
        "source_hashes":hashes(SOURCES),"input_hashes":hashes(INPUTS),
        "conventions":{"leading_action":"sqrt|g|[-A R_L+L_m]", "leading_variation":"-A*(G_mu_nu-T_mu_nu/(2A))*delta g^{mu nu}",
                       "signature":"+---; curvature and stress conventions inherited from spherical-action",
                       "original_curvature_terms":"-cW C²-cR R²-cEuler E4-cBox boxR",
                       "coefficients":"full A,cW,cR,C_F,V symbolic; Dirac cR=0 does not set the full cR to zero",
                       "order":"first order in a common formal parameter multiplying cW,cR; O(eps²) not claimed",
                       "frame":"canonical X=8*pi*A*r² belongs to the redefined EFT metric; original observables require pullback"},
        "field_map":"g_old^{mu nu}=g_new^{mu nu}-(2*cW/A)*(R^{mu nu}-R*g^{mu nu}/6+(T^{mu nu}-T*g^{mu nu}/3)/(2*A))+(cR/A)*(R-T/(2*A))*g^{mu nu}",
        "stress_contact":"-cW/(2*A²)*(T_mu_nu*T^{mu nu}-T²/3)-cR*T²/(4*A²)",
        "Euler_coefficient_after_map":"cEuler+cW, in the original static-energy sign convention",
        "exact_residuals":exact,
        "spherical_source":{"area":"X/(2*A)","kinetic":"K=(grad chi)²","stress_trace":"4*V before quantum anomaly/renormalization",
                            "Q_T":"K²/(2*area²)+4*rho_EM²-4*V²/3",
                            "rho_EM":"rho_B+mu²*chi²/(2*area)",
                            "contact_density":str(sp.factor(spherical_contact(X,chi,sp.Symbol('K'),A,V,CF,mag,cW,cR))),
                            "extra_light_fields":"include all sectors in total T BEFORE squaring; cross terms must not be lost"},
        "canonical":{"first_order_correction":[str(sp.factor(x)) for x in C1],
                     "kinetic0":str(K0),"constraint_structure_check":structures,
                     "new_local_matter_density":"d(X)*K²+F(X,chi), with d=-cW/(2*A*X)",
                     "first_order_initial_data_enlarged":False,
                     "old_exact_linear_geometry_integration_reused_without_correction":False},
        "constant_curvature_control":{"inverse_metric_scale_shift":str(shift_vac),"apparent_vacuum_Lagrangian_shift":str(contact_vac),
                                       "physical_curvature_shift_after_pullback":"0 through first EFT order"},
        "quantum_and_boundary_requirements":{"Jacobian":"J_F=det[delta g_old/delta g_new], with other field, ghost, cutoff and state transformations included",
                                             "Jacobian_set_to_one":False,
                                             "Euler_box_and_induced_boundary_terms_discarded":False,
                                             "stress_contact_expectation_factorized_into_mean_stresses":False,
                                             "composite_products_and_state_matched":False},
        "scope":{"fourD_topological_identity_applied_to_fiveD":False,
                 "correction_added_on_top_of_the_same_unremoved_curvature_term":False,
                 "full_nonlocal_spectral_action_replaced_by_truncation":False,
                 "derivative_expansion_validity_at_original_neck_established":False,
                 "full_theory_ghost_or_causality_verdict_derived":False,
                 "complete_couplings_state_or_self_sourced_solution_derived":False,
                 "physical_vacuum_energy_removed_by_change_of_variables":False},
        "regime":"requires |c_i R|/A and |c_i T|/A² small, external scales below the matched expansion scales; boundaries/state/regulator must be transformed consistently",
        "comparison":{"fields":"all","float_atol":3e-13,"float_rtol":3e-13,
                      "exact":"identities, symbolic formulas, structures, strings and source/input hashes","exceptions":[]},
    }


def main():
    parser=argparse.ArgumentParser(description=__doc__);group=parser.add_mutually_exclusive_group()
    group.add_argument("--check",action="store_true");group.add_argument("--output",type=Path)
    args=parser.parse_args()
    if args.output and args.output.exists():raise FileExistsError("refusing to overwrite recorded evidence")
    expected=json.loads(OUTPUT.read_text()) if args.check else None
    if args.check:compare(expected["source_hashes"],hashes(SOURCES),"$/source_hashes");compare(expected["input_hashes"],hashes(INPUTS),"$/input_hashes")
    result=calculate()
    if args.check:compare(expected,result);print("curvature EFT: every first-order identity and record field reproduced; prior generators were not run")
    elif args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        with args.output.open("x") as stream:json.dump(result,stream,indent=2,sort_keys=True,allow_nan=False);stream.write("\n")
        print(args.output)
    else:print(json.dumps(result,indent=2,allow_nan=False))


if __name__=="__main__":main()
