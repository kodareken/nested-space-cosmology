#!/usr/bin/env python3
"""Execute the homogeneous GeneralKSSameActionHistoryFunctional gate."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys


ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
from recursive_horizons.nsc_general_ks_same_action_history import GeneralKSSameActionHistoryFunctional


OUTPUT=ROOT/"results/development/nsc-general-ks-same-action-history.json"
INPUTS=(
    "results/development/nsc-joint-history-state-bvp.json",
    "results/development/nsc-mode-resolved-cauchy-state.json",
    "results/development/nsc-landau-cauchy-isometry.json",
    "results/development/nsc-background-projection.json",
    "results/development/charged-ctp-neck-source.json",
    "results/development/scale-binding.json",
    "results/development/gauge-source.json",
    "results/development/causal-common-functional.json",
    "results/development/spherical-action.json",
)
SOURCES=(
    "src/recursive_horizons/nsc_general_ks_same_action_history.py",
    "scripts/derive_nsc_general_ks_same_action_history.py",
    "docs/nsc-general-ks-same-action-history.md",
)


def hashes(paths):return {path:hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in paths}


def compare(expected,actual,path="$",atol=3e-13,rtol=3e-13):
    if isinstance(expected,dict):
        if not isinstance(actual,dict) or expected.keys()!=actual.keys():raise AssertionError(f"keys differ at {path}")
        for key in expected:compare(expected[key],actual[key],f"{path}/{key}",atol,rtol)
    elif isinstance(expected,list):
        if not isinstance(actual,list) or len(expected)!=len(actual):raise AssertionError(f"list differs at {path}")
        for i,(a,b) in enumerate(zip(expected,actual)):compare(a,b,f"{path}/{i}",atol,rtol)
    elif isinstance(expected,float):
        if isinstance(actual,bool) or not isinstance(actual,(int,float)) or abs(expected-actual)>atol+rtol*abs(expected):raise AssertionError(f"number differs at {path}")
    elif type(expected) is not type(actual) or expected!=actual:raise AssertionError(f"value differs at {path}")


def calculate():
    bvp,state,isometry,background,charged,scale,gauge,causal,spherical=(json.loads((ROOT/p).read_text()) for p in INPUTS)
    branch=scale["development_branch"];locked=state["locked_inputs"]
    if not (locked["magnetic_flux"]==branch["magnetic_flux"]==4 and locked["Omega"]==branch["omega"] and locked["zeta"]==branch["zeta"] and locked["V_full"]==0):raise ValueError("locked branch changed")
    if not state["gate"]["seed_serialization_pass"] or state["gate"]["full_contract_pass"]:raise ValueError("unexpected state gate")
    if isometry["gate"]["physical_U_L0_pass"] or bvp["gate"]["full_joint_BVP_pass"]:raise ValueError("upstream history gate changed")
    if causal["gate"]["actual_physical_CausalSource_complete"] or spherical["scope"]["self_sourced_geometry_or_physical_couplings_solved"]:raise ValueError("upstream action owner changed")
    tensor=background["completed_tensor_input"]
    compact_local=charged["runs"]["base"]["compact_local_source"]["child_frame"]
    if compact_local["T01"]!=0 or gauge["source_accounting"]["full_absolute_renormalized_stress_computed"]:raise ValueError("local ledger assumptions changed")
    owner=GeneralKSSameActionHistoryFunctional(
        branch["coefficients"]["A_Wilsonian_at_magnetic_scale"],
        branch["magnetic_flux"],branch["omega"],branch["zeta"],
        branch["coefficients"]["V_full_relational"],
    )
    tolerance=3e-11
    certificate=owner.homogeneous_nonexistence_certificate(state_T01=tensor["T01"],tolerance=tolerance)
    if certificate["stationary_history_exists_in_declared_class"]:raise ArithmeticError("homogeneous momentum obstruction disappeared")
    variation=certificate["variation"]
    return {
        "schema":"NSC-GENERAL-KS-SAME-ACTION-HISTORY-v1",
        "status":"NON-EXISTENCE PASS: no smooth stationary history exists in the declared homogeneous KS class with the locked nonzero-momentum seed state",
        "source_hashes":hashes(SOURCES),"input_hashes":hashes(INPUTS),
        "owner":{
            "name":"GeneralKSSameActionHistoryFunctional",
            "history_class":"smooth frequency-diagonal homogeneous Kantowski-Sachs ADM fields N(tau), beta(tau), q_ADM(tau)=a_parallel(tau), r(tau); no spatial gradients and no transmitting tilted interface",
            "declared_variation":"delta Gamma_one_CTP[g,C[g]] / delta g_Delta^A(tau)=0 for A=N,beta,q_ADM,r",
            "execution_order":"evaluate lapse/shift constraints on fixed boundary data before interior collocation",
        },
        "locked_inputs":{
            "A":locked["A"],"magnetic_flux":locked["magnetic_flux"],"Omega":locked["Omega"],"zeta":locked["zeta"],"V_full":locked["V_full"],"occupation_or_parameter_refit":False,
        },
        "same_action_momentum_ledger":{
            **variation,
            "physical_mode_state":"serialized ModeResolvedCauchyState",
            "magnetic_gauge_reason":"homogeneous static magnetic F_theta_phi has no normal-axial Poynting component",
            "local_induced_reason":"local scalar metric invariants and their fourth-order reference are axial-reflection even and have no homogeneous T01",
            "compact_local_verified_zero":compact_local["T01"],
            "relational_vacuum_verified_zero":locked["V_full"],
            "boundary_condition":"fixed smooth Cauchy endpoints; no declared momentum shell",
        },
        "nonexistence_certificate":certificate,
        "functional_component_ledger":{
            "ModeResolved_state_evolution":"implemented blockwise for a supplied history",
            "canonical_CTP_determinant":"existing CausalCommonFunctional; supplied-history evaluator",
            "Einstein_A":"locked and owns the geometric constraints",
            "magnetic_gauge_C":"already owned; diagonal in homogeneous KS",
            "V_full":"zero by the relational law",
            "fourth_order_reference_history":"not needed to decide nonexistence because its homogeneous T01 is zero by axial reflection",
            "local_induced_history":"not needed to decide nonexistence because its homogeneous T01 is zero by axial reflection",
            "boundary_interface":"no momentum shell is declared",
            "locked_local_coefficients":{
                "C_gauge":branch["coefficients"]["C_Wilsonian_at_magnetic_scale"],
                "C_Weyl":branch["coefficients"]["C_Weyl"],
                "C_Euler":branch["coefficients"]["C_Euler"],
                "C_boxR":branch["coefficients"]["C_boxR"],
            },
            "reference_momentum_policy":"use the parity-completed renormalized reference; the one-sided stored positive-frequency rows are never interpreted as a vacuum T01",
        },
        "variations":{
            "beta_seed_residual":variation["normalized_shift_variation"],
            "N_seed_residual":bvp["route_A_dynamic_seed"]["Hamiltonian_constraint_residual"],
            "q_ADM_and_r_interior_variations_evaluated":False,
            "reason":"one failed independent constraint is sufficient to exclude every smooth history in the scoped class",
        },
        "finite_solution":{
            "g_star":None,"C_star":None,"finite_stress":None,"null_signs":None,"updated_constraints":None,
            "reason":"no admissible homogeneous history reaches an interior stationarity solve from the locked seed",
        },
        "escape_conditions_outside_scope":[
            "derive a same-action counterflow Delta T01=-0.0012238701558050947",
            "allow spatially inhomogeneous spherical data whose momentum constraint has gradient terms",
            "declare a physical boundary momentum shell",
            "change the initial quantum state",
        ],
        "gate":{
            "variation_residual_below_tolerance":False,"nonexistence_in_declared_class_proved":True,"full_history_search_needed":False,"finite_stress_fabricated":False,"coupled_evolution_reopened":False,
        },
        "nonexistence_scope":"frequency-diagonal homogeneous no-interface truncation only; a tilted frequency-mixing interface or inhomogeneous spherical geometry is outside this certificate",
        "full_class_unimplemented_history_terms":[
            "GeneralKSFourthOrderReferenceHistory with all four ADM variations and parity-completed momentum subtraction",
            "GeneralKSLocalInducedHistory with nodal action forces rather than constant endpoint forces",
            "transmitting tilted-Landau interface variation with its frequency-mixing Cauchy map",
        ],
        "scope":{
            "diagnostic_duration_or_profile_selected":False,"metric_timestep_started":False,"A_q_Omega_zeta_or_V_full_changed":False,"new_source_compensator_or_dark_fluid":False,"w_of_a_H_of_z_or_Hessian":False,"standard_GR_rederived":False,
        },
        "comparison":{"fields":"all","exact":"schema, strings, booleans and source/input hashes","float_atol":3e-13,"float_rtol":3e-13,"exceptions":[]},
    }


def main():
    parser=argparse.ArgumentParser(description=__doc__);group=parser.add_mutually_exclusive_group();group.add_argument("--check",action="store_true");group.add_argument("--output",type=Path);args=parser.parse_args();result=calculate()
    if args.check:compare(json.loads(OUTPUT.read_text()),result);print("general-KS same-action homogeneous nonexistence certificate reproduced")
    elif args.output:
        if args.output.exists():raise FileExistsError("refusing to overwrite recorded evidence")
        args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n");print(args.output)
    else:print(json.dumps(result,indent=2,sort_keys=True))


if __name__=="__main__":main()
