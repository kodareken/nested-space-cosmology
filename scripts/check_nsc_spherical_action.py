#!/usr/bin/env python3
"""Match the local NSC source action to the established canonical geometry owner."""
import argparse
import json
from pathlib import Path
import sys

import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
from recursive_horizons.nsc_spherical_action import (
    potentials, poisson_tensor, canonical_constraints, primitive, casimir, pg_null_coframe,
)
from check_nsc_compact_boundary_action import authenticated_record
from check_nsc_vacuum_charge_matching import compare, hashes

OUTPUT=ROOT/"results/development/spherical-action.json"
SOURCES=("scripts/check_nsc_spherical_action.py","src/recursive_horizons/nsc_spherical_action.py",
         "docs/nsc-spherical-action.md","docs/nsc-compact-matching.md",
         "scripts/check_nsc_compact_boundary_action.py","scripts/check_nsc_vacuum_charge_matching.py")
INPUTS=("results/development/compact-matching.json","results/development/gauge-source.json")


def calculate():
    for path in INPUTS:authenticated_record(path)
    exact={}

    def zero(name,value):
        if isinstance(value,sp.MatrixBase):
            assert value.applyfunc(sp.simplify)==sp.zeros(*value.shape),name
        else:assert sp.simplify(value)==0,name
        exact[name]="0"

    X,A4,C_F,mag,r=sp.symbols("X A4 C_F q_mag r",positive=True)
    Xp,Xm,V4=sp.symbols("Xplus Xminus V4",real=True)
    pars=potentials(X,A4,V4,C_F,mag)
    R2,grad_r2,grad_X2=sp.symbols("R2 grad_r_squared grad_X_squared",real=True)
    reduced_r=-4*sp.pi*A4*(r*r*R2+2*grad_r2-2)-4*sp.pi*V4*r*r-2*sp.pi*C_F*mag**2/r**2
    reduced_X=-X*R2/2+pars["U"]*grad_X2/2-pars["V"]
    zero("spherical_bulk_normalization",reduced_X.subs({X:8*sp.pi*A4*r*r,grad_X2:(16*sp.pi*A4*r)**2*grad_r2})-reduced_r)
    g4_squared=1/(4*C_F)
    zero("charge_mass_same_coefficient",pars["mu_squared"].subs(X,8*sp.pi*A4*r*r)-mag*g4_squared/(4*sp.pi**2*r*r))
    chi=sp.Symbol("chi",real=True)
    Vm=pars["mu_squared"]*chi**2/2
    zero("area_source_matches_Maxwell_pressure",-X*sp.diff(Vm,X)-Vm)

    # The known spherical IBP term cancels the radius part of the same 4D GHY.
    K2,nr,epsilon=sp.symbols("K2 normal_derivative_r epsilon",real=True)
    spherical_IBP=16*sp.pi*A4*epsilon*r*nr
    GHY4=-8*sp.pi*A4*epsilon*r*r*(K2+2*nr/r)
    GHY2=-8*sp.pi*A4*epsilon*r*r*K2
    zero("same_GHY_radius_cancellation",spherical_IBP+GHY4-GHY2)

    # Check the canonical component sign, not just an overall vacuum sign.
    n,h,xt,xx,nt,nx,ht,hx,xtx=sp.symbols("N h X_t X_x N_t N_x h_t h_x X_tx",real=True)
    nxx,htt=sp.symbols("N_xx h_tt",real=True)
    wc_t,wc_x=nx/h,ht/n
    xplus=(xt/n-xx/h)/sp.sqrt(2)
    xminus=(-xt/n-xx/h)/sp.sqrt(2)
    e2,e3,b2,b3=-h/sp.sqrt(2),h/sp.sqrt(2),n/sp.sqrt(2),n/sp.sqrt(2)
    base=(X,n,h,xt,xx,nx,ht)
    dx_values=(xx,nx,hx,xtx,sp.Symbol('X_xx'),nxx,sp.Symbol('h_tx'))
    dx=lambda value:sum(sp.diff(value,a)*b for a,b in zip(base,dx_values))
    curvature=2/(n*h)*(nxx/h-nx*hx/h**2-htt/n+ht*nt/n**2)
    total=pars["U"]*xplus*xminus+pars["V"]
    canonical=(X*(htt/n-ht*nt/n**2)+xplus*(-ht/sp.sqrt(2))+xminus*(ht/sp.sqrt(2))
               +wc_t*(xx+xminus*e3-xplus*e2)
               +b2*(dx(xplus)+wc_x*xplus-e3*total)
               +b3*(dx(xminus)-wc_x*xminus+e2*total))
    boundary=X*wc_t+xplus*b2+xminus*b3
    bulk=n*h*(-X*curvature/2+pars["U"]*(xt*xt/n**2-xx*xx/h**2)/2-pars["V"])
    zero("canonical_component_action_and_boundary",canonical-dx(boundary)-bulk)

    p=(X,Xp,Xm);q=sp.symbols("omega emin eplus",nonzero=True)
    dp=sp.symbols("X_x Xplus_x Xminus_x");dq=sp.symbols("omega_x emin_x eplus_x")
    cp,cpp,momentum,mp=sp.symbols("chi_x chi_xx pi_chi pi_chi_x",real=True)
    G=canonical_constraints(p,q,dp,pars["U"],pars["V"],chi,cp,momentum,pars["mu_squared"])
    P=poisson_tensor(*p,pars["U"],pars["V"])
    Pmatter=poisson_tensor(*p,pars["U"],pars["V"]+Vm)
    # Independent Legendre transform fixes the positive-energy matter sign.
    bar2,bar3,velocity=sp.symbols("bar2 bar3 chi_velocity",real=True)
    volume=bar2*q[2]-bar3*q[1]
    inverse00=-2*q[1]*q[2]/volume**2
    inverse01=(bar2*q[2]+bar3*q[1])/volume**2
    inverse11=-2*bar2*bar3/volume**2
    scalar_L=volume*((inverse00*velocity**2+2*inverse01*velocity*cp+inverse11*cp**2)/2-Vm)
    solved_velocity=(momentum-volume*inverse01*cp)/(volume*inverse00)
    scalar_H=(momentum*velocity-scalar_L).subs(velocity,solved_velocity)
    matter_G=G-sp.Matrix(dp)-P*sp.Matrix(q)
    zero("matter_Legendre_to_same_constraints",scalar_H+bar2*matter_G[1]+bar3*matter_G[2])
    for i,j,k in ((0,1,2),):
        jac=sum(P[i,l]*sp.diff(P[j,k],p[l])+P[j,l]*sp.diff(P[k,i],p[l])+P[k,l]*sp.diff(P[i,j],p[l]) for l in range(3))
        zero("geometric_Poisson_Jacobi",jac)
    cas=casimir(*p,A4,V4,C_F,mag)
    zero("Casimir_Poisson_kernel",P*sp.Matrix([sp.diff(cas,x) for x in p]))
    zero("Casimir_primitive",sp.diff(primitive(X,A4,V4,C_F,mag),X)-pars["V"]/sp.sqrt(X))

    # Smeared functional Poisson brackets: derivative-delta terms are checked,
    # then integrated by parts. This includes the X-dependent mass term.
    a,ax,b,bx=sp.symbols("a a_x b b_x")
    variables=(*p,*q,chi,cp,momentum);jets=(*dp,*dq,cp,cpp,mp)
    derivative=lambda f:sp.expand(sum(sp.diff(f,z)*dz for z,dz in zip(variables,jets)))
    def functional(i,z,zx):
        f=G[i]
        return ([z*sp.diff(f,x) for x in q],
                [z*sp.diff(f,x)-zx*int(k==i) for k,x in enumerate(p)],
                z*sp.diff(f,chi)-zx*sp.diff(f,cp)-z*derivative(sp.diff(f,cp)),
                z*sp.diff(f,momentum))
    for i,j in ((0,1),(0,2),(1,2)):
        Fq,Fp,Fc,Fm=functional(i,a,ax);Hq,Hp,Hc,Hm=functional(j,b,bx)
        B=sp.expand(sum(u*v-w*z for u,v,w,z in zip(Fq,Hp,Fp,Hq))+Fc*Hm-Fm*Hc)
        B0=B.coeff(a).coeff(b);B10=B.coeff(ax).coeff(b);B01=B.coeff(a).coeff(bx)
        zero(f"coupled_constraint_{i+1}{j+1}_derivative_delta",B01-B10)
        zero(f"coupled_constraint_{i+1}{j+1}_closure",B0-derivative(B10)-sum(sp.diff(Pmatter[i,j],x)*g for x,g in zip(p,G)))

    # Pure geometry and fixed matter histories have the same full Jacobian
    # operator as the temporal-gauge ghost operator. Keep the full operator.
    d0=sp.Symbol("partial_0",commutative=True)
    T=pars["U"]*Xp*Xm+pars["V"]+Vm
    M=sp.diag(d0,d0,d0)+sp.Matrix([[sp.diff(Pmatter[i,1],x) for x in p] for i in range(3)])
    J=sp.Matrix([[d0,-1,0],[0,d0,0],[sp.diff(T,X),sp.diff(T,Xp),d0+sp.diff(T,Xm)]])
    zero("full_ghost_and_constraint_Jacobian",M-J)
    # Noncommuting finite time blocks: determinant factorization is a finite
    # identity, not an independent continuum regulator prescription.
    D=sp.Matrix([[1,0,0],[-1,1,0],[0,-1,1]])
    B=sp.diag(sp.Rational(1,3),sp.Rational(1,5),sp.Rational(1,7))
    K=sp.Matrix([[1,2,0],[3,0,1],[0,4,2]]);L=sp.diag(2,3,5);Z=sp.zeros(3)
    full=sp.BlockMatrix([[D,-sp.eye(3),Z],[Z,D,Z],[K,L,D+B]]).as_explicit()
    assert D*B-B*D!=sp.zeros(3)
    zero("finite_noncommuting_FP_factorization",full.det()-D.det()**2*(D+B).det())

    # Geometry responds to the same matter history, rather than a new weight.
    charge_dot=sp.Symbol("chi_dot",real=True)
    p_dot=sp.Matrix([Xp,-charge_dot**2,-T])
    cas_dot=sum(sp.diff(cas,x)*v for x,v in zip(p,p_dot))
    zero("matter_Casimir_balance",cas_dot+(Xp*Vm+Xm*charge_dot**2)/sp.sqrt(X))
    q3=sp.Symbol("q3",positive=True)
    zero("coframe_measure_evolution",q3*Xp/(2*sp.sqrt(X))+pars["U"]*Xp*q3*sp.sqrt(X))
    C0,Mgeo=sp.symbols("C0 M_geom",real=True)
    metric=2/sp.sqrt(X)*(C0-primitive(X,A4,V4,C_F,mag))/(32*sp.pi*A4)
    conversion=-32*sp.pi*A4*sp.sqrt(8*sp.pi*A4)*Mgeo
    expected=1-2*Mgeo/r+C_F*mag**2/(4*A4*r*r)-V4*r*r/(6*A4)
    zero("known_vacuum_metric_normalization",metric.subs({X:8*sp.pi*A4*r*r,C0:conversion})-expected)
    grad_r=sp.Symbol("grad_r_squared",real=True)
    mass=r*(1+grad_r+C_F*mag**2/(4*A4*r*r)-V4*r*r/(6*A4))/2
    mass_casimir=cas.subs(Xp*Xm,-(16*sp.pi*A4*r)**2*grad_r/2).subs(X,8*sp.pi*A4*r*r)
    zero("Casimir_to_spherical_mass_function",mass_casimir+32*sp.pi*A4*sp.sqrt(8*sp.pi*A4)*mass)

    beta=sp.Symbol("beta",real=True);coframe=pg_null_coframe(beta)
    zero("PG_null_coframe_metric",coframe.T*sp.Matrix([[0,1],[1,0]])*coframe-sp.Matrix([[1-beta**2,-beta],[-beta,-1]]))
    zero("PG_coframe_nondegenerate",coframe.det()-1)
    neck_P=P.subs({Xp:0,Xm:0})
    assert neck_P[1,2]==-pars["V"]

    return {
        "schema":"NSC-SPHERICAL-ACTION-v1",
        "status":"local spherical action and coupled canonical constraints matched to existing first-order gravity",
        "source_hashes":hashes(SOURCES),"input_hashes":hashes(INPUTS),
        "conventions":{"fourD_action":"sqrt|g4|[-A4 R_L - V4 - C_F F²]","signature":"+---; R_E=-R_L",
                       "X":"8*pi*A4*r²; area field, not a resolution compensator",
                       "coefficients":"full effective A4,V4,C_F kept symbolic; Dirac-only values not substituted",
                       "q":"absolute magnetic flux for the displayed charged mode; extra neutral CFT sectors required when q>1",
                       "GKV_sign":"physical bulk action is minus GKV (2.9); chapter-7 canonical component convention retained",
                       "range":"X>0, A4>0, C_F>0; oriented canonical patch q3>0,q2<0 has spacelike slices; local two-derivative spherical action only"},
        "potentials":{k:str(v) for k,v in pars.items()},
        "bulk_action_density_without_measure":str(reduced_X),
        "charge_potential":str(Vm),"exact_residuals":exact,
        "canonical":{"q":[str(x) for x in q],"p":[str(x) for x in p],
                     "constraints":[str(sp.factor(x)) for x in G],
                     "smeared_bracket":"{G_i[a],G_j[b]}=integral a*b*(partial Pmatter_ij/partial p_k)*G_k",
                     "Hamiltonian":"-bar_q_i G_i",
                     "Liouville_measure":"Dq Dp Dchi Dpi_chi, with constraint multipliers and BRST ghosts",
                     "temporal_gauge":"bar_q=(omega_0,e^-_0,e^+_0)=(0,1,0)",
                     "ghost_matrix":[[str(x) for x in row] for row in M.tolist()],
                     "formal_ghost_determinant":"(det partial_0)^2 det(partial_0+Xplus*U(X)); domains/zero modes must agree",
                     "geometry_momentum_equations":[str(x) for x in p_dot],
                     "q3_equation":"partial_0 q3=U(X)*Xplus*q3; q3/q3_initial=sqrt(X_initial/X) with auxiliary generating sources set to zero, retaining the charge history",
                     "matter_measure":"standard covariant scalar measure contains sqrt(q3); its bosonization/UV matching and zero-mode prescription remain required",
                     "local_gravity_configuration_degrees":0},
        "boundary":{"spherical_IBP":"integral_boundary epsilon*n.dot(grad X)",
                    "remaining_GHY":"-integral_boundary sqrt|h1| epsilon X K2",
                    "canonical_to_bulk_spatial_boundary":"[X omega_0+Xplus e^-_0+Xminus e^+_0]_ends",
                    "corners_null_boundaries_and_transmission_matched":False},
        "mass_and_domain":{"Casimir":str(cas),"Casimir_to_geometric_mass":str(conversion),
                           "parent_energy":"M_geom/G, G=1/(16*pi*A4), in the asymptotically flat parent limit with unit-normalized Killing time and the stated reference subtraction",
                           "outgoing_power_balance":"dM_geom/dtau=-G*P_infinity under that parent normalization with the complete source and fixed magnetic/vacuum reference terms",
                           "canonical_time_is_automatically_parent_Killing_time":False,
                           "radius_as_coordinate_valid_at_static_neck":False,
                           "area_as_field_valid_at_positive_radius_neck":True,
                           "Poisson_rank_at_Xplus_Xminus_zero":"2 if V(X)!=0; 0 only when V(X)=0; not every neck is a rank-drop point"},
        "scope":{"new_spherical_GR_or_quantization_theorem_claimed":False,
                 "GKV_auxiliary_torsion_is_the_5D_physical_three_form":False,
                 "pure_geometry_ghost_cancellation_removes_matter_measure":False,
                 "full_nonlocal_or_higher_curvature_action_quantized":False,
                 "full_5D_or_resolution_field_measure_derived":False,
                 "quantum_regulated_constraint_algebra_verified":False,
                 "global_gauge_state_or_source_boundary_selected":False,
                 "self_sourced_geometry_or_physical_couplings_solved":False},
        "comparison":{"fields":"all","float_atol":3e-13,"float_rtol":3e-13,
                      "exact":"all symbolic identities, structures, strings and hashes","exceptions":[]},
    }


def main():
    parser=argparse.ArgumentParser(description=__doc__);group=parser.add_mutually_exclusive_group()
    group.add_argument("--check",action="store_true");group.add_argument("--output",type=Path)
    args=parser.parse_args()
    if args.output and args.output.exists():raise FileExistsError("refusing to overwrite recorded evidence")
    expected=json.loads(OUTPUT.read_text()) if args.check else None
    if args.check:
        compare(expected["source_hashes"],hashes(SOURCES),"$/source_hashes");compare(expected["input_hashes"],hashes(INPUTS),"$/input_hashes")
    result=calculate()
    if args.check:compare(expected,result);print("spherical action: every field and exact identity reproduced; prior generators were not run")
    elif args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        with args.output.open("x") as stream:json.dump(result,stream,indent=2,sort_keys=True,allow_nan=False);stream.write("\n")
        print(args.output)
    else:print(json.dumps(result,indent=2,allow_nan=False))


if __name__=="__main__":main()
