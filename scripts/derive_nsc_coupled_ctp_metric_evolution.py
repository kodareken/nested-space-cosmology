#!/usr/bin/env python3
"""Gate the first same-covariance stress evaluation on the selected neck."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src"))

from recursive_horizons.nsc_coupled_ctp_metric_evolution import (
    source_revaluation_preflight,
)


OUTPUT = ROOT/"results/development/nsc-coupled-ctp-metric-evolution.json"
INPUTS = (
    "results/development/nsc-constraint-complete-neck.json",
    "results/development/nsc-child-metric-backreaction.json",
    "results/development/nsc-background-projection.json",
    "results/development/charged-compact-ctp-completion.json",
    "results/development/charged-ctp-neck-source.json",
    "results/development/causal-common-functional.json",
    "results/development/adm-source-constraints.json",
    "results/development/child-state.json",
)
SOURCES = (
    "src/recursive_horizons/nsc_coupled_ctp_metric_evolution.py",
    "scripts/derive_nsc_coupled_ctp_metric_evolution.py",
    "docs/nsc-coupled-ctp-metric-evolution.md",
)


def hashes(paths):
    return {
        path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
        for path in paths
    }


def compare(expected, actual, path="$", atol=2e-13, rtol=2e-13):
    if isinstance(expected, dict):
        if not isinstance(actual, dict) or expected.keys() != actual.keys():
            raise AssertionError(f"keys differ at {path}")
        for key in expected:
            compare(expected[key], actual[key], f"{path}/{key}", atol, rtol)
    elif isinstance(expected, list):
        if not isinstance(actual, list) or len(expected) != len(actual):
            raise AssertionError(f"list differs at {path}")
        for index, (left, right) in enumerate(zip(expected, actual)):
            compare(left, right, f"{path}/{index}", atol, rtol)
    elif isinstance(expected, float):
        if isinstance(actual, bool) or not isinstance(actual, (int, float)):
            raise AssertionError(f"numeric type differs at {path}")
        if abs(expected-actual) > atol+rtol*abs(expected):
            raise AssertionError(f"number differs at {path}: {expected} != {actual}")
    elif type(expected) is not type(actual) or expected != actual:
        raise AssertionError(f"value differs at {path}: {expected!r} != {actual!r}")


def calculate():
    selected, failed, projection, compact, angular, causal, adm, child = (
        json.loads((ROOT/path).read_text()) for path in INPUTS
    )
    if not selected["gate"]["passes"]:
        raise ValueError("constraint-complete initial geometry is required")
    if failed["initial_constraint_gate"]["passes"]:
        raise ValueError("historical unit-radius gate changed")
    if causal["gate"]["actual_physical_CausalSource_complete"]:
        raise ValueError("physical causal source record changed")
    if adm["scope"]["joint_geometry_solution_found"]:
        raise ValueError("ADM owner already contains a joint solution")

    # The accepted source records contain integrated rows/totals and covariance
    # eigenvalue diagnostics.  They intentionally do not persist C_j(k,k').
    forbidden_state_keys = {
        "mode_resolved_cauchy_covariance",
        "mode_resolved_cauchy_state",
        "cauchy_covariance_artifact",
    }
    if forbidden_state_keys & compact.keys() or forbidden_state_keys & angular.keys():
        raise ValueError("source records now expose a mode-resolved state; rerun this gate")

    frame = selected["source_selected_frame"]
    stress = frame["stress"]
    geometry = selected["source_selected_initial_geometry"]
    preflight = source_revaluation_preflight(
        seed_radius=geometry["seed_r"],
        selected_radius=geometry["r"],
        rho_landau=stress["rho"],
        p_sphere_landau=stress["p_sphere"],
    )
    if abs(preflight["fixed_covariance_density_log_radius_vertex"]) < 1e-12:
        raise ArithmeticError("metric source vertex unexpectedly vanished")

    locked = selected["locked_inputs"]
    initial_row = {
        "T": 0.0,
        "a_parallel": geometry["a_parallel"],
        "r": geometry["r"],
        "H_parallel": geometry["H_parallel"],
        "H_sphere": geometry["H_sphere"],
        "rho_fixed_tensor_gate": stress["rho"],
        "T01_fixed_tensor_gate": stress["T01"],
        "p_parallel_fixed_tensor_gate": stress["p_parallel"],
        "p_sphere_fixed_tensor_gate": stress["p_sphere"],
        "fixed_tensor_Hamiltonian_constraint": selected["constraints"]["residuals"]["Hamiltonian"],
        "fixed_tensor_momentum_constraint": selected["constraints"]["residuals"]["momentum"],
        "null_plus_fixed_tensor_gate": frame["radial_nulls"]["Landau_plus"],
        "null_minus_fixed_tensor_gate": frame["radial_nulls"]["Landau_minus"],
        "same_covariance_stress_on_selected_geometry": None,
    }
    return {
        "schema": "NSC-COUPLED-CTP-METRIC-EVOLUTION-v1",
        "status": (
            "STOP AT T=0 SOURCE RE-EVALUATION: the integrated neck tensor does "
            "not contain the mode-resolved Gaussian covariance required by "
            "unitary CTP evolution on the source-selected geometry"
        ),
        "source_hashes": hashes(SOURCES),
        "input_hashes": hashes(INPUTS),
        "locked_inputs": {
            "q": locked["q"],
            "Omega": locked["Omega"],
            "zeta": locked["zeta"],
            "A": locked["A"],
            "V_full": locked["V_full"],
            "occupations_changed": False,
            "seed_retained_a_parallel": geometry["a_parallel"],
            "seed_retained_H_parallel": geometry["H_parallel"],
        },
        "prestep_metric_vertex": preflight,
        "state_information_gate": {
            "completed_records_store": [
                "integrated rho, T01, p_parallel and p_sphere",
                "mode labels, quadrature metadata and aggregate row stresses",
                "minimum/maximum covariance eigenvalue diagnostics",
            ],
            "completed_records_do_not_store": [
                "complex C_jk for every compact/angular/frequency mode",
                "canonical isometry from the old neck slice to the Landau slice",
                "mode-resolved renormalized reference on the selected geometry",
            ],
            "four_integrated_stress_moments_determine_future_commutators": False,
            "reason": (
                "future pressure uses C_dot=-i[H[g],C]; traces of C with four "
                "current vertices do not determine the mode commutators"
            ),
            "actual_physical_CausalSource_complete_in_reused_record": causal["gate"][
                "actual_physical_CausalSource_complete"
            ],
        },
        "trajectory": {
            "rows": [initial_row],
            "time_steps_started": 0,
            "break_time": 0.0,
            "break_stage": "same-covariance stress evaluation before first metric step",
            "constraint_break_time": None,
            "null_sign_break_time": None,
            "reason": (
                "the fixed-tensor constraint row is not a CTP stress "
                "re-evaluation on r_star"
            ),
        },
        "decision": {
            "fixed_tensor_initial_constraints_pass": True,
            "fixed_tensor_initial_null_signs_negative": True,
            "coupled_CTP_metric_trajectory_defined": False,
            "fake_fluid_or_constant_w_closure_used": False,
            "chosen_next_owner": "ModeResolvedCauchyState",
            "required_contract": {
                "basis": "canonical half-density chi=sqrt(a_parallel*r^2)*psi",
                "payload": (
                    "frequency weights and complex Gaussian covariance for "
                    "every retained compact/angular channel"
                ),
                "slice_map": "old child frame to source-selected Landau normal",
                "renormalization": (
                    "the existing fourth-order reference and local induced "
                    "allocation evaluated on the same metric history"
                ),
                "consumer": "CausalCommonFunctional plus ADM evolution",
            },
        },
        "rethinking_decision": {
            "discarded_moment_fluid": "would invent a pressure law from four moments",
            "discarded_fixed_tensor_reuse": "nonzero radius vertex proves it is not the same source evaluation",
            "discarded_old_generator_replay_only": "would reconstruct the old slice but not define the Landau slice/history",
            "selected": "persist and transport the actual mode-resolved covariance before integration",
        },
        "scope": {
            "new_metric_ansatz_or_source_channel": False,
            "old_source_generator_rerun": False,
            "A_q_Omega_zeta_or_V_full_changed": False,
            "unit_radius_Bronnikov_restored": False,
            "w_of_a_or_fluid_closure_introduced": False,
            "H_of_z_likelihood_or_Hessian_computed": False,
            "trajectory_claimed": False,
        },
        "comparison": {
            "fields": "all",
            "exact": "keys, types, strings and source/input hashes",
            "float_atol": 2e-13,
            "float_rtol": 2e-13,
            "exceptions": [],
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--check", action="store_true")
    group.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = calculate()
    if args.check:
        expected = json.loads(OUTPUT.read_text())
        compare(expected["source_hashes"], hashes(SOURCES), "$/source_hashes")
        compare(expected["input_hashes"], hashes(INPUTS), "$/input_hashes")
        compare(expected, result)
        print("coupled evolution T=0 source-state gate reproduced; no source was rerun")
    elif args.output:
        if args.output.exists():
            raise FileExistsError("refusing to overwrite recorded evidence")
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True)+"\n")
        print(args.output)
    else:
        print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
