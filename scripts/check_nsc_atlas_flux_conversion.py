"""Reachability, isolated bounce controls and one conditional conversion history."""
import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np
import sympy as sp

from recursive_horizons.nsc_atlas_flux_transition import FluxLaw
from recursive_horizons.nsc_atlas_flux_conversion import (
    ScalarFlux, sector_minimum, isolated_bounce, evolve_prepared,
)

ROOT=Path(__file__).resolve().parents[1]


def symbolic_checks():
    m2,mu,q,e,b,phi=sp.symbols("m2 mu q e b phi",real=True)
    w2=m2+mu**2
    eps=m2/w2
    minimum=lambda Q:-mu*Q/w2
    vacuum=lambda Q:b+eps*Q*Q/2
    U=lambda P,Q:b+m2*P*P/2+(Q+mu*P)**2/2
    excitation=mu**2*e**2/(2*w2)
    release=eps*e*q-e*e/2
    delta,velocity0=sp.symbols("delta velocity0",real=True)
    excited_parent=vacuum(q)+(velocity0**2+w2*delta**2)/2
    excited_child=vacuum(q-e)+(velocity0**2+w2*(delta-mu*e/w2)**2)/2
    x,v,H,a,k,W,rv=sp.symbols("x v H a k W rv",real=True)
    rho=rv+(v*v+W*W*x*x)/2
    pressure=(v*v-W*W*x*x)/2
    derivative=lambda expr: (sp.diff(expr,x)*v+sp.diff(expr,v)*(-3*H*v-W*W*x)
                             +sp.diff(expr,H)*(-k*v*v/2)+sp.diff(expr,a)*H*a)
    checks={
        "completed_square":U(phi,q)-vacuum(q)-w2*(phi-minimum(q))**2/2,
        "prepared_excitation":U(minimum(q),q-e)-vacuum(q-e)-excitation,
        "transition_energy":vacuum(q)-vacuum(q-e)-excitation-release,
        "wall_force_uses_actual_field":U(phi,q)-U(phi,q-e)-e*(q+mu*phi-e/2),
        "excited_parent_transfer":excited_parent-excited_child-e*(eps*q+mu*delta-e/2),
        "scalar_conservation":derivative(rho)+3*H*v*v,
        "Friedmann_propagation":derivative(3*H*H-k*rho),
        "comoving_virial":derivative(a**3*x*v)-2*a**3*pressure,
    }
    result={name:str(sp.factor(expr)) for name,expr in checks.items()}
    if any(x!="0" for x in result.values()):
        raise RuntimeError(result)
    return result


def calculation():
    start=time.process_time()
    exact=symbolic_checks()
    model=ScalarFlux()
    coarse,mc=evolve_prepared()
    fine,mf=evolve_prepared(rtol=2e-12,atol=2e-14)
    alternative,ma=evolve_prepared(clock="efolds",rtol=2e-12,atol=2e-14)
    end=fine.t[-1]
    nodes=np.linspace(0,min(coarse.t[-1],end),257)
    state=fine.sol(nodes)
    x,v,N,H,work=state
    energy=(v*v+model.omega**2*x*x)/2
    rv=mf["preparation"]["child_vacuum"]
    friedmann=float(np.max(np.abs(3*H*H-model.kappa*(rv+energy))))
    conservation=float(np.max(np.abs(energy+work-mf["preparation"]["oscillation_energy"])))
    refinement=float(np.max(np.abs(coarse.sol(nodes)-state)))
    endpoint=fine.y[:,-1]
    alt=alternative.y[:,-1]
    clock_gap=max(abs(endpoint[0]-alt[0]),abs(endpoint[1]-alt[1]),abs(end-alt[2]))
    final_energy=(endpoint[1]**2+model.omega**2*endpoint[0]**2)/2
    comoving_ratio=final_energy*np.exp(3*endpoint[2])/mf["preparation"]["oscillation_energy"]
    if max(friedmann,conservation,clock_gap,refinement)>3e-8:
        raise RuntimeError("independent evolution/conservation check failed")
    paths=("src/recursive_horizons/nsc_atlas_flux_conversion.py",
           "src/recursive_horizons/nsc_atlas_flux_transition.py",
           "scripts/check_nsc_atlas_flux_conversion.py","tests/test_nsc_atlas_flux_conversion.py")
    record=dict(
        schema="nsc-atlas-flux-conversion-v1",
        decision="Test selection capability, then connect a specified flux transfer to residual vacuum and an evolving scalar. Full P002 nucleation is deferred.",
        pure_flux=dict(aligned=sector_minimum("1.5",".5"),offset=sector_minimum("1.55",".5"),
                       exact_aligned_smallest_positive_Lambda="1/8",
                       invariant="f modulo e is unchanged by any integer-charge sequence; rates cannot create absent levels.",
                       local_degeneracy="A constant-flux bulk supplies only a cosmological-constant source; away from transitions it has no additional propagating mode.",
                       isolated_controls=[isolated_bounce(FluxLaw(),f) for f in (1.5,1.,.5)],
                       bounce_scale="Actions listed are Bhat; physical B=Bhat*Lref^2/kappa_phys in normalized variables. No determinant, absolute rate or semiclassical validity is asserted.",
                       boundary_limit="P002 has nonzero Kottler masses and a finite dust parent; these isolated vacuum instantons are not its preparation or rates."),
        scalar_candidate=dict(
            potential="rho_b + m^2*phi^2/2 + (q+mu*phi)^2/2",
            epsilon=model.epsilon,omega=model.omega,bare_density=model.bare_density,
            normalization="phi_hat=sqrt(kappa_phys)*phi, q_hat=sqrt(kappa_phys)*Lref*q, omega_hat=Lref*omega, rho_hat=kappa_phys*Lref^2*rho. Arithmetic kappa=1 selects no SI scale.",
            imported="Kaloper-Sorbo0811.1989 Eqs1-3 for scalar/four-form mixing; additional positive m^2 is declared, without an axion symmetry or ultraviolet-completion claim.",
            preparation=mf["preparation"],
            no_work_condition="e <= 2*epsilon*q_parent for e>0 in this rapid ground-state transfer, not a condition on every possible tunnelling history",
            initial_excitation_to_vacuum_ratio=mf["preparation"]["oscillation_energy"]/rv,
            strong_suppression_control=ScalarFlux(epsilon=.01).prepare(),
            excited_parent_control=model.prepare(parent_displacement=.03,parent_velocity=.4),
            ground_state_no_work_bound="For rho_b=0, q_parent>0 and 0<epsilon<1/2, Eosc/rho_v <= 4 epsilon(1-epsilon)/(1-2epsilon)^2 if the rapid ground-state transfer releases nonnegative energy. This bound does not apply to an excited parent.",
            scope="A prescribed rapid-transfer preparation and a homogeneous post-transfer patch. No moving scalar/membrane junction or causal creation of a whole homogeneous slice has been solved."),
        evolution=dict(efolds=2.,scale_factor=float(np.exp(endpoint[2])),proper_time=float(end),
                       initial_H=mf["initial_H"],initial_omega_over_H=mf["initial_omega_over_H"],
                       final_H=float(endpoint[3]),vacuum_density=rv,
                       final_oscillation_density=float(final_energy),final_oscillation_to_vacuum=float(final_energy/rv),
                       rho_osc_a_cubed_over_initial=float(comoving_ratio),
                       statement="The oscillating component dilutes approximately as dust; no perturbation/clustering or thermalization result is implied.",
                       pressure_identity="d(a^3*chi*chi_dot)/dt=2*a^3*p_osc; a zero-endpoint full-cycle comoving pressure integral vanishes exactly."),
        checks=dict(symbolic_zero_identities=exact,Friedmann_max_residual=friedmann,
                    integrated_energy_max_residual=conservation,tolerance_max_state_difference=refinement,
                    independent_efold_clock_max_endpoint_gap=float(clock_gap),
                    rhs_evaluations=[mc["rhs_evaluations"],mf["rhs_evaluations"],ma["rhs_evaluations"]]),
        remaining="Boundary-driven preparation and coupled scalar-wall evolution using the derived actual-field normal force; perturbations, measured couplings/scale, thermalization and observed density fractions.",
        numerical_uncertainty="Tolerance and clock comparisons only, not a rigorous interval enclosure, model error or observational uncertainty.",
        field_grids=0,process_cpu_seconds=time.process_time()-start,
        source_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        reproduce="python scripts/run_python.py scripts/check_nsc_atlas_flux_conversion.py")
    if record["process_cpu_seconds"]>60:
        raise RuntimeError("finite calculation exceeded its CPU budget")
    return record


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path)
    args=parser.parse_args()
    encoded=json.dumps(calculation(),indent=2,allow_nan=False)+"\n"
    if len(encoded.encode())>16000:
        raise RuntimeError("compact summary budget exceeded")
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(encoded,encoding="utf-8",newline="\n")
    else:
        print(encoded,end="")
