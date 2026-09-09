#!/usr/bin/env python3
"""Apply the known projected QED2 interaction to the recorded NSC exterior."""
import argparse
import json
from pathlib import Path
import sys

import numpy as np
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
from recursive_horizons.nsc_gauge_source import Exterior, outgoing_power, strength, transmission
from recursive_horizons.nsc_lorentzian import geometry
from check_nsc_compact_boundary_action import authenticated_record
from check_nsc_compact_casimir import compare, native
from check_nsc_vacuum_charge_matching import hashes

OUTPUT=ROOT/"results/development/gauge-source.json"
SOURCES=("scripts/check_nsc_gauge_source.py","src/recursive_horizons/nsc_gauge_source.py",
         "docs/nsc-gauge-source.md","src/recursive_horizons/nsc_lorentzian.py",
         "scripts/check_nsc_compact_boundary_action.py","scripts/check_nsc_compact_casimir.py",
         "scripts/check_nsc_vacuum_charge_matching.py")
INPUTS=("results/development/horizon-source.json","results/development/compact-matching.json")


def calculate():
    previous=authenticated_record(INPUTS[0])
    authenticated_record(INPUTS[1])
    source=previous["benchmark_Unruh_per_abs_q"]
    bg=Exterior(source["reused_horizon_rho"],source["reused_surface_gravity"])
    radius_h=source["reused_horizon_areal_radius"]
    exact={}

    def zero(name,value):
        assert sp.simplify(value)==0,name
        exact[name]="0"

    r,g,q,X,E=sp.symbols("r g4 q X E",positive=True)
    e2=g*g/(4*sp.pi*r*r)
    mu2=q*e2/sp.pi
    electric=e2*sp.sqrt(q/sp.pi)*X
    gauge_lagrangian=E*E/(2*e2)-sp.sqrt(q/sp.pi)*E*X
    zero("gauge_elimination_gives_charge_potential",gauge_lagrangian.subs(E,electric)+mu2*X*X/2)
    rho_e=electric**2/(2*g*g)
    p_sphere=mu2*X*X/(8*sp.pi*r*r)
    zero("charge_potential_is_radial_Maxwell_energy",rho_e-p_sphere)
    zero("sphere_pressure_from_radius_variation",sp.diff(-mu2*X*X/2,r)/(8*sp.pi*r)-p_sphere)
    # On-shell div2 t^a_b = X² partial_b(mu²)/2. The angular pressure
    # restores the full 4D conservation law after division by sphere area.
    rp=sp.Symbol("r_prime",real=True)
    zero("spherical_Ward_with_Maxwell_pressure",X*X*sp.diff(mu2,r)*rp/(8*sp.pi*r*r)+2*rp/r*p_sphere)
    kappa=sp.Symbol("kappa",positive=True)
    temperature=kappa/(2*sp.pi)
    zero("boson_CFT_flux_matches_one_complex_Dirac",temperature**2/(2*sp.pi)*(sp.pi**2/6)-kappa*kappa/(48*sp.pi))
    rho=sp.Symbol("rho",positive=True)
    c=sp.Symbol("c",positive=True)
    zero("potential_tail_primitive",sp.diff(c*sp.atan(rho),rho)-c/(1+rho*rho))
    J,Ttt,Ttr,A,beta=sp.symbols("J Ttt Ttr A beta",real=True)
    zero("future_regular_PG_horizon_null_source",(-beta*Ttt-A*Ttr).subs({A:0,beta:1})+Ttt)

    # Stable evaluation is used only where the original formula loses digits.
    overlaps=[]
    for rho in (bg.horizon_rho+5e-5,19.9,20.1,40.):
        difference=bg.A(rho)-float(geometry(rho)[2])
        assert abs(difference)<2e-11
        overlaps.append({"rho":rho,"difference_from_recorded_formula":difference})

    methods=[]
    for frequency in (.004,.04,.12):
        a=transmission(bg,frequency,.5)
        b=transmission(bg,frequency,.5,method="field")
        assert abs(a["transmission"]-b["transmission"])<1e-7
        assert a["relative_current_residual"]<1e-8 and b["relative_current_residual"]<1e-8
        methods.append({"frequency":frequency,"amplitudes":a,"direct_field":b,
                        "transmission_difference":b["transmission"]-a["transmission"]})

    cases=[]
    for coupling in (.25,.5,1.):
        result=outgoing_power(bg,coupling,nodes=64)
        assert 0<result["charged_ratio_lower_with_tail_bounds"]<=result["charged_power_over_one_free_channel"]<1
        assert result["charged_power_over_one_free_channel"]<=result["charged_ratio_upper_with_tail_bounds"]
        assert result["maximum_relative_current_residual"]<1e-8
        cases.append({"g4":coupling,"flux":1,"strength":strength(coupling),"power":result,
                      "conditional_horizon_null_stress":-result["parent_Killing_power"]/(4*np.pi*radius_h**2),
                      "mu_over_first_Landau_gap":coupling/(2*np.pi*np.sqrt(2)),
                      "fourD_alpha":coupling**2/(4*np.pi)})
    base=cases[1]["power"]
    refinements={}
    controls={"frequency_nodes_48":{"nodes":48},"frequency_nodes_96":{"nodes":96},
              "outer_phase_120":{"phase_extent":120.},"horizon_offset_1e-9":{"horizon_offset":1e-9},
              "ODE_rtol_2e-10":{"rtol":2e-10}}
    for name,options in controls.items():
        settings={"nodes":64,**options}
        result=outgoing_power(bg,.5,**settings)
        difference=result["charged_power_over_one_free_channel"]-base["charged_power_over_one_free_channel"]
        assert abs(difference)<2e-6,name
        refinements[name]={"power":result,"charged_ratio_difference":difference}
    free=outgoing_power(bg,0.,flux=1)
    assert free["parent_Killing_power"]==source["parent_Killing_power"]

    return native({
        "schema":"NSC-GAUGE-SOURCE-v1",
        "status":"known projected Dirac-gauge interaction mapped to the unwrapped source and transmission",
        "source_hashes":hashes(SOURCES),"input_hashes":hashes(INPUTS),
        "conventions":{
            "units":"hbar=c_light=1, original benchmark throat radius unit",
            "domain":"recorded spherical exterior, same reflecting transverse compact candidate; photon and fermion angular projection stated separately",
            "retained_sector":"one canonical massless 4D Dirac -> abs(q_mag) LLL flavours; spherical radial U(1) gauge field only",
            "state":"specified projected Unruh state, no parent incoming bath, zero electric-flux sector",
            "scattering_orientation":"unit ingoing horizon amplitude; real-potential reciprocity supplies the outgoing transmission",
            "statistics":"q-1 neutral CFT units and one Bose charged collective mode",
            "e2_squared":"g4^2/(4*pi*r^2)","collective_mu_squared":"abs(q_mag)*g4^2/(4*pi^2*r^2)",
            "potential":"A(rho)*abs(q_mag)*g4^2/[4*pi^2*(1+rho^2)]",
            "collective_field":"electric-displacement/current variable from bosonization; not a new fundamental scalar or the resolution compensator",
            "g4":"unselected full matched coupling; the Dirac-only C_F is not its complete value",
            "base_parameters":{"flux":1,"frequency_nodes":64,"minimum_omega_over_T":1e-4,
                               "maximum_omega_over_T":24.,"outer_phase_extent":60.,
                               "horizon_offset":1e-8,"ODE_rtol":2e-9},
        },
        "exact_mapping_residuals":exact,
        "reused_horizon":{"rho":bg.horizon_rho,"surface_gravity":bg.surface_gravity,"areal_radius":radius_h,
                          "temperature":source["Hawking_temperature"],"one_free_channel_power":source["parent_Killing_power"]},
        "geometry_evaluation_overlap":overlaps,
        "independent_scattering_methods":methods,
        "coupling_controls":cases,"independent_refinements":refinements,"free_limit_reused":free,
        "source_accounting":{
            "interaction_correction":"Gamma_scalar[mu(r)]-Gamma_scalar[0], plus required local/zero-mode matching terms",
            "sphere_pressure_operator":"mu(r)^2*X^2/(8*pi*r^2), with composite expectation values consistently renormalized",
            "horizon_null_relation":"T(K,K)_h=-P_infinity/(4*pi*r_h^2), conditional on stationary conserved current and future-regular PG components",
            "bare_Dirac_UV_Maxwell_term_added_again":False,
            "full_absolute_renormalized_stress_computed":False,
        },
        "error_scope":{
            "domain_tail_formula":"S_tail=c*[atan(rho_left)-atan(rho_h)+atan(1/rho_right)]/(2*omega)",
            "transmission_bounds":"sech^2(theta+S_tail) <= T_full <= sech^2(max(theta-S_tail,0))",
            "assumptions":"real stationary integrable potential and disjoint-tail transfer-matrix composition",
            "numerical_caveat":"comparison bounds evaluated in floating point using the computed finite transfer; ODE and quadrature errors assessed separately, not rigorous interval arithmetic",
            "transmission_extraction":"T=1/(1+abs(b)^2) uses the exact Wronskian; raw 1/abs(a)^2 and current defects retained",
        },
        "scope":{
            "new_Schwinger_or_Hawking_law_claimed":False,
            "full_4D_reduction_or_magnetic_catalysis_controlled":False,
            "higher_photon_Landau_and_compact_modes_included":False,
            "physical_g4_selected_or_fitted":False,
            "collective_scale_is_fourD_fermion_or_photon_rest_mass":False,
            "asymptotic_mass_gap_present":False,
            "global_unified_state_or_self_sourced_geometry_solved":False,
            "cosmological_density_rate_or_dark_fraction_predicted":False,
        },
        "comparison":{"fields":"all","float_atol":3e-9,"float_rtol":3e-8,
                      "exact":"keys, structures, strings, non-float types and source/input hashes","exceptions":[]},
    })


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group();group.add_argument("--check",action="store_true");group.add_argument("--output",type=Path)
    args=parser.parse_args()
    if args.output and args.output.exists():raise FileExistsError("refusing to overwrite recorded evidence")
    expected=json.loads(OUTPUT.read_text()) if args.check else None
    if args.check:
        compare(expected["source_hashes"],hashes(SOURCES),"$/source_hashes")
        compare(expected["input_hashes"],hashes(INPUTS),"$/input_hashes")
    result=calculate()
    if args.check:
        compare(expected,result);print("gauge source: every field reproduced; prior generators were not run")
    elif args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        with args.output.open("x") as stream:json.dump(result,stream,indent=2,sort_keys=True,allow_nan=False);stream.write("\n")
        print(args.output)
    else:print(json.dumps(result,indent=2,allow_nan=False))


if __name__=="__main__":main()
