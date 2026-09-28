#!/usr/bin/env python3
"""Apply known free-state evolution to the actual expanding-child volume."""
import argparse
import json
from math import pi,sqrt
from pathlib import Path
import sys

import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
from recursive_horizons.nsc_child_state import (
    interior_factors,relative_energy_bound,channel_energy_norm,PARALLEL_RATIO_MINIMUM,
    relative_null_state_source,
)
from recursive_horizons.nsc_warped_source import angular_modes
from recursive_horizons.nsc_spinor_bridge import weyl_matrices
from check_nsc_compact_boundary_action import authenticated_record
from check_nsc_vacuum_charge_matching import compare,hashes

OUTPUT=ROOT/"results/development/child-state.json"
SOURCES=("scripts/check_nsc_child_state.py","src/recursive_horizons/nsc_child_state.py",
         "docs/nsc-child-state.md","src/recursive_horizons/nsc_spinor_bridge.py",
         "docs/nsc-compact-mass-map.md","scripts/check_nsc_compact_mass_map.py",
         "scripts/check_nsc_compact_boundary_action.py","scripts/check_nsc_vacuum_charge_matching.py")
INPUTS=("results/nsc-5-clock-horizon.json","results/nsc-10-influence.json",
        "results/development/spectral-endpoint.json","results/development/horizon-source.json")


def calculate():
    clock,_,endpoint,horizon=[authenticated_record(p) for p in INPUTS]
    assert clock["child_geometry"]["asymptotic_directional_H"]=="sqrt(3pi)"
    exact={}
    def zero(name,value):
        if isinstance(value,sp.MatrixBase):
            assert value.applyfunc(sp.simplify)==sp.zeros(*value.shape),name
        else:assert sp.simplify(value)==0,name
        exact[name]="0"

    x=sp.Symbol("rho",real=True);u=sp.Symbol("u",positive=True)
    r=sp.sqrt(1+x*x)
    A=1+3*x+3*(1+x*x)*(sp.atan(x)-sp.pi/2)
    W=3*(sp.pi/2-sp.atan(x))-(1+3*x)/(1+x*x)
    zero("parallel_ratio_from_recorded_metric",W+A/r**2)
    zero("parallel_ratio_monotonicity",sp.diff(W,x)-2*(x-3)/(1+x*x)**2)
    Wu=3*(sp.pi-sp.atan(u))+(3*u-u*u)/(1+u*u)
    zero("future_ratio_series",sp.series(Wu,u,0,4).removeO()-(3*sp.pi-u*u-2*u**3))
    zero("future_conformal_time_metric",(1/(u*u)/sp.sqrt((1+u*u)*Wu/u**2)/(sp.sqrt(1+u*u)/u))**2-1/((1+u*u)**2*Wu))

    # Canonical half-density removes the existing expansion connection.
    hp,hs,a,R=sp.symbols("H_parallel H_sphere a_parallel r",positive=True)
    v=a*R**2
    dv=sp.diff(v,a)*hp*a+sp.diff(v,R)*hs*R
    zero("canonical_volume_connection",dv/(2*v)-hp/2-hs)
    frame=weyl_matrices();beta=frame["gamma"][0];az,angular=frame["alpha"][:2]
    m,k,kap=sp.symbols("m k kappa",real=True)
    H=m*beta+k/a*az+kap/R*angular
    zero("channel_norm_for_interior_scaling",H*H-(m*m+k*k/a**2+kap*kap/R**2)*sp.eye(4))
    dotH=-hp*k/a*az-hs*kap/R*angular
    # Operator identity implies work/conservation for any common-evolution
    # covariance difference, with two equal angular pressure components.
    pparallel=k/a*az
    psphere=kap/(2*R)*angular
    zero("same_action_anisotropic_work",dotH+hp*pparallel+2*hs*psphere)
    dm=sp.Symbol("m_dot",real=True)
    zero("variable_mass_is_extra_work",dotH+dm*beta+hp*pparallel+2*hs*psphere-dm*beta)

    # New proper-interior projection of the ALREADY computed null constants.
    du,dv=sp.symbols("delta_tu delta_tv",real=True)
    unull,vnull=sp.Matrix([-1/a,1]),sp.Matrix([1/a,1])
    difference=du*unull*unull.T+dv*vnull*vnull.T
    zero("null_constants_to_interior_density",difference[0,0]-(du+dv)/a**2)
    zero("null_constants_to_interior_flux",difference[0,1]/a-(dv-du)/a**2)
    rho_delta=(du+dv)/(4*sp.pi*a*a*R*R)
    drho=sp.diff(rho_delta,a)*hp*a+sp.diff(rho_delta,R)*hs*R
    zero("null_state_difference_work",drho+(hp+2*hs)*rho_delta+hp*rho_delta)

    factors=[interior_factors(x) for x in (0.,-1.,-4.,-16.)]
    for g in factors:
        assert PARALLEL_RATIO_MINIMUM<=g["parallel_over_radius"]<sqrt(3*pi)
        mass,momentum,kap=pi/2,1.,1.
        norm=channel_energy_norm(g["radius"],g["a_parallel"],mass,momentum,kap)
        envelope=mass+(abs(momentum)/PARALLEL_RATIO_MINIMUM+abs(kap))/g["radius"]
        assert norm<=envelope
    # These rows display only normalized moment bounds, not evolved states
    # or a prediction of their initial occupations.
    bound_rows=[{"radius":r,"massless_M0_K1":relative_energy_bound(r,0.,1.),
                 "massive_M1_K1":relative_energy_bound(r,1.,1.)} for r in (1.,2.,4.,16.,64.)]
    coefficients=[row["massless_M0_K1"]["total"]*row["radius"]**4 for row in bound_rows]
    assert max(coefficients)-min(coefficients)<1e-15
    state_constant=horizon["benchmark_Unruh_per_abs_q"]["parent_Killing_power"]
    state_rows=[relative_null_state_source(x,0.,-state_constant) for x in (0.,-1.,-4.,-16.,-64.)]
    for row in state_rows:assert abs(row["delta_density"])<=row["absolute_density_bound"]
    gap_rows=[]
    for radius in (1.,2.,4.,16.):
        gap=angular_modes(1,radius,1)[1][0]
        gap_rows.append({"radius":radius,"abs_flux":1,"first_nonzero_angular_squared":gap,
                         "gap_squared_over_child_H_squared":gap/(3*pi)})
    return {"schema":"NSC-CHILD-STATE-v1","status":"conditional dilution of homogeneous free Dirac state differences on the actual child",
      "source_hashes":hashes(SOURCES),"input_hashes":hashes(INPUTS),"exact_residuals":exact,
      "conventions":{"metric":"dT²-a_parallel(T)² dz²-r(T)² dOmega²; z is the parent's static time used as an interior spatial coordinate",
                     "units":"hbar=c=L_throat=1; restoring L gives s_min=sqrt(3*pi/2-1)/L, k and m inverse length, kappa dimensionless; mode measure includes longitudinal momentum density",
                     "free_tower":"fixed ell, time-independent compact domain/warp; m_n=n*pi/ell, before additional link/gauge interactions",
                     "canonical_field":"chi=(a_parallel*r²)^(1/2)*psi; unit S² measure retained",
                     "channel_H":"m*beta+(k/a_parallel)*alpha_z+(kappa/r)*alpha_sphere; physical angular sectors and degeneracies retained in channel measure",
                     "volume":"4*pi*a_parallel*r² per unit coordinate z; no finite global spatial volume assumed",
                     "stress":"difference between two states on the same operator, geometry and renormalization prescription; local state-independent terms cancel"},
      "future_geometry":{"squared_ratio_series":"a_parallel²/r²=3*pi-u²-2*u³+O(u4), u=-1/rho",
                         "conformal_metric":"g/r²=du²/[(1+u²)² W(u)]-W(u) dz²-dOmega²",
                         "W_at_future":3*pi,"s_min":PARALLEL_RATIO_MINIMUM,
                         "remaining_conformal_time_bound":"0<eta_future-eta(u)<=atan(u)/s_min",
                         "ratios":factors},
      "state_bound":{"norm":"N_j=||Delta C_j||_1; conserved for Delta C_j=U_j Delta C_j(initial) U_j^dagger",
                     "mass_moment":"M=sum_j/integral dnu_j |m_j| N_j < infinity",
                     "momentum_moment":"K=sum_j/integral dnu_j (|k_j|/s_min+|kappa_j|) N_j < infinity",
                     "density_bound":"|Delta rho| <= (M+K/r)/(4*pi*s_min*r³)",
                     "massless_rate":"O(r^-4) if M=0", "massive_envelope":"O(r^-3) for finite M,K",
                     "pressure_bounds":"|Delta p_parallel|+2|Delta p_sphere| <= K/(4*pi*s_min*r^4)",
                     "energy_ledger":"Delta rho_dot+(H_parallel+2H_sphere)Delta rho+H_parallel Delta p_parallel+2H_sphere Delta p_sphere=0",
                     "mass_change_work":"m_dot*tr(beta Delta C)/(4*pi*a_parallel*r²) must be retained when m varies",
                     "normalized_bound_controls":bound_rows},
      "required_hypotheses":["common closed free evolution on the actual expanding interior",
                              "translation-invariant/mode-diagonal and angularly symmetric states, or the specified spatially averaged density",
                              "finite weighted trace-norm moments for the complete retained tower and continuous mode measure",
                              "fixed compact size, transmission domain and masses; no unaccounted reservoir injection",
                              "one common renormalization prescription for both states"],
      "existing_parent_state_control":{"states":"free LLL Unruh minus Hartle-Hawking, per |q|; only their recorded null constants differ",
                                       "delta_tu":0.,"delta_tv":-state_constant,"interior_projection":state_rows,
                                       "asymptotic_r4_density":-state_constant/(12*pi*pi),
                                       "absolute_reference_stress_identified":False},
      "angular_applicability":{"imported_gap_formula":"kappa_1²=(1+|q|)/r²; at fixed finite q its ratio to H_child² tends to zero",
                               "normalized_flux_one_controls":gap_rows,
                               "LLL_only_approximation_uniform_at_late_child":False,
                               "vacuum_limit_may_be_interchanged_with_unsummed_angular_modes":False,
                               "state_bound_limit_interchange":"permitted only by the finite M,K dominating moments; not a bound on vacuum zero-point sums"},
      "reference_state_gap":{"difference_decay_selects_absolute_vacuum":False,
                             "Hadamard_condition_alone_proves_uniform_future_bounds":False,
                             "round_S4_is_global_continuation_of_actual_interior":False,
                             "actual_parent_state_summability_verified":False,
                             "canonical_future_instantaneous_vacuum_assumed":False,
                             "full_finite_regulator_causal_continuation_established":False,
                             "required":"construct or import an applicable Hadamard reference on the actual domain whose absolute stress approaches the same-regulator invariant endpoint; prove parent/reference moment conditions"},
      "unchanged_endpoint_input":{"radius":endpoint["endpoint_comparisons"][0]["full"]["radius"],
                                  "Dirac_log_radius_derivative":endpoint["endpoint_comparisons"][0]["full"]["d_log_radius"],
                                  "spectral_sum_rerun":False},
      "scope":{"interacting_gauge_or_metric_state_attractor_proved":False,"persistent_vacuum_energy_ruled_out":False,
               "zero_flux_implies_transient_acceleration":False,"cosmological_Q_or_dark_fraction_predicted":False,
               "self_sourcing_completed":False,"new_general_cosmic_no_hair_theorem_claimed":False},
      "comparison":{"fields":"all","float_atol":3e-13,"float_rtol":3e-13,"exact":"symbolic identities, formulas, structures, strings and source/input hashes","exceptions":[]}}


def main():
    parser=argparse.ArgumentParser(description=__doc__);group=parser.add_mutually_exclusive_group()
    group.add_argument("--check",action="store_true");group.add_argument("--output",type=Path)
    args=parser.parse_args()
    if args.output and args.output.exists():raise FileExistsError("refusing to overwrite recorded evidence")
    expected=json.loads(OUTPUT.read_text()) if args.check else None
    if args.check:compare(expected["source_hashes"],hashes(SOURCES));compare(expected["input_hashes"],hashes(INPUTS))
    result=calculate()
    if args.check:compare(expected,result);print("child state: all bounds, identities and record fields reproduced; previous generators were not run")
    elif args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        with args.output.open("x") as stream:json.dump(result,stream,indent=2,sort_keys=True,allow_nan=False);stream.write("\n")
        print(args.output)
    else:print(json.dumps(result,indent=2,allow_nan=False))


if __name__=="__main__":main()
