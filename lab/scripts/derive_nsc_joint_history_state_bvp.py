#!/usr/bin/env python3
"""Audit the joint history/state BVP and record its exact missing selector."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src"))

from recursive_horizons.nsc_joint_history_state_bvp import nonuniqueness_witness


OUTPUT = ROOT/"results/development/nsc-joint-history-state-bvp.json"
INPUTS = (
    "results/development/nsc-landau-cauchy-isometry.json",
    "results/development/nsc-mode-resolved-cauchy-state.json",
    "results/development/nsc-child-metric-backreaction.json",
    "results/development/nsc-constraint-complete-neck.json",
    "results/development/causal-common-functional.json",
    "results/development/spherical-action.json",
    "results/development/relational-vacuum-normalization.json",
)
SOURCES = (
    "src/recursive_horizons/nsc_joint_history_state_bvp.py",
    "scripts/derive_nsc_joint_history_state_bvp.py",
    "docs/nsc-joint-history-state-bvp.md",
)


def hashes(paths):
    return {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in paths}


def compare(expected, actual, path="$", atol=3e-13, rtol=3e-13):
    if isinstance(expected, dict):
        if not isinstance(actual, dict) or expected.keys() != actual.keys():
            raise AssertionError(f"keys differ at {path}")
        for key in expected: compare(expected[key], actual[key], f"{path}/{key}", atol, rtol)
    elif isinstance(expected, list):
        if not isinstance(actual, list) or len(expected) != len(actual):
            raise AssertionError(f"list differs at {path}")
        for index,(left,right) in enumerate(zip(expected,actual)):
            compare(left,right,f"{path}/{index}",atol,rtol)
    elif isinstance(expected,float):
        if isinstance(actual,bool) or not isinstance(actual,(int,float)):
            raise AssertionError(f"numeric type differs at {path}")
        if abs(expected-actual)>atol+rtol*abs(expected):
            raise AssertionError(f"number differs at {path}")
    elif type(expected) is not type(actual) or expected!=actual:
        raise AssertionError(f"value differs at {path}")


def calculate():
    isometry,state,seed_fail,selected,causal,spherical,vacuum = (
        json.loads((ROOT/path).read_text()) for path in INPUTS
    )
    payload_path=ROOT/isometry["payload"]["path"]
    if hashlib.sha256(payload_path.read_bytes()).hexdigest()!=isometry["payload"]["sha256"]:
        raise ValueError("conditional-isometry payload hash mismatch")
    with np.load(payload_path,allow_pickle=False) as payload:
        witness=nonuniqueness_witness({name:payload[name] for name in payload.files})
    if witness["maximum_absolute_U_difference"]<=0 or witness["maximum_absolute_C_difference"]<=0:
        raise ArithmeticError("nonuniqueness witness collapsed")
    if causal["gate"]["actual_physical_CausalSource_complete"]:
        raise ValueError("causal owner now contains the missing source")
    if spherical["scope"]["self_sourced_geometry_or_physical_couplings_solved"]:
        raise ValueError("spherical owner now contains a physical solution")
    if vacuum["gate"]["V_full"]!=0:
        raise ValueError("locked relational vacuum changed")
    seed_constraints=seed_fail["initial_constraint_gate"]["normalized_constraints"]
    locked=selected["locked_inputs"]
    return {
        "schema":"NSC-JOINT-HISTORY-STATE-BVP-v1",
        "status":(
            "FAIL / ROUTE C: dynamic seed evolution fails its initial constraints, "
            "while the two-point BVP has multiple unitary state histories and no "
            "executable same-action history selector"
        ),
        "source_hashes":hashes(SOURCES),"input_hashes":hashes(INPUTS),
        "selected_route":{
            "choice":"C_nonuniqueness_bound_and_stop",
            "A_dynamic_seed_attempted":False,
            "A_stop_is_algebraic_before_timestep":True,
            "B_physical_history_selected":False,
        },
        "locked_inputs":{
            "magnetic_flux":locked["q"],"Omega":locked["Omega"],
            "zeta":locked["zeta"],"A":locked["A"],"V_full":locked["V_full"],
            "occupation_or_parameter_refit":False,
        },
        "route_A_dynamic_seed":{
            "seed_surface":"serialized r=1 ModeResolvedCauchyState",
            "Hamiltonian_constraint_residual":seed_constraints["hamiltonian_geometry_minus_source"],
            "momentum_constraint_residual":seed_constraints["momentum_geometry_minus_source"],
            "can_start_ADM_evolution":False,
            "reason":"constraints are initial-data equations and fail before a state or metric step",
        },
        "route_B_two_point_BVP":{
            "seed_state_serialized":state["gate"]["seed_serialization_pass"],
            "conditional_U_API_pass":isometry["gate"]["conditional_U_unitarity_pass"],
            "same_endpoints":isometry["control_results"]["same_endpoint_data"],
            "nonuniqueness_witness":witness,
            "unique_history_selected":False,
            "finite_rstar_stress_available":False,
        },
        "existing_selector_audit":{
            "H_sphere_zero": "endpoint condition satisfied by both controls; does not select the interior history",
            "CTP_equal_history_normalization": "normalizes Gamma[g,g]=0; does not extremize the physical average history",
            "CTP_energy_work_Ward": "holds for each declared unitary history; does not distinguish them",
            "relational_zero_tadpole": "selects V_full=0, not lapse/shift/radius duration or profile",
            "recursive_scale_stationarity": "locks Omega/zeta but supplies no KS history functional",
            "spherical_constraints": "restrict each Cauchy slice but do not supply the missing quantum stress history",
            "declared_gauge_invariant_mismatch_with_stationary_equation": None,
        },
        "exact_missing_selector":{
            "name":"GeneralKSSameActionHistoryFunctional",
            "equation":"delta Gamma_one_CTP[g,C[g]] / delta g_Delta^A(tau)=0 for A=N,beta,q_ADM,r",
            "must_include":[
                "blockwise ModeResolvedCauchyState Dirac evolution",
                "general-KS fourth-order reference history",
                "general-KS local induced coefficient history counted once",
                "boundary/interface variation and the declared endpoint conditions",
            ],
            "current_causal_owner_role":"evaluates a supplied finite history; it does not solve or select one",
            "new_physical_law_required":False,
            "missing_executable_part_of_declared_action":True,
        },
        "finite_stress_and_constraints":{
            "finite_rstar_stress":None,
            "null_signs":None,
            "updated_constraint_residuals":None,
            "reason":"different allowed histories give different C_L before the absent history-dependent reference/local terms are evaluated",
            "fabricated_value_used":False,
        },
        "gate":{
            "route_A_pass":False,"route_B_unique_selector_pass":False,
            "nonuniqueness_lower_bound_pass":True,
            "full_joint_BVP_pass":False,"coupled_evolution_reopened":False,
            "missing_selector_named_exactly":True,
        },
        "scope":{
            "diagnostic_duration_selected_physical":False,"metric_timestep_started":False,
            "A_q_Omega_zeta_or_V_full_changed":False,"pointwise_spin_boost_used":False,
            "unit_radius_tensor_or_linearized_density_reused":False,
            "w_of_a_H_of_z_Hessian_or_dark_fit":False,"new_source_channel":False,
        },
        "comparison":{
            "fields":"all","exact":"schema, strings, booleans and source/input/artifact hashes",
            "float_atol":3e-13,"float_rtol":3e-13,"exceptions":[],
        },
    }


def main():
    parser=argparse.ArgumentParser(description=__doc__);group=parser.add_mutually_exclusive_group()
    group.add_argument("--check",action="store_true");group.add_argument("--output",type=Path)
    args=parser.parse_args();result=calculate()
    if args.check:
        compare(json.loads(OUTPUT.read_text()),result);print("joint history/state BVP uniqueness failure reproduced from authenticated controls")
    elif args.output:
        if args.output.exists():raise FileExistsError("refusing to overwrite recorded evidence")
        args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n");print(args.output)
    else:print(json.dumps(result,indent=2,sort_keys=True))


if __name__=="__main__":main()
