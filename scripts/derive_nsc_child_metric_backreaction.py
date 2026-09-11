#!/usr/bin/env python3
"""Gate locked child backreaction before evolving the completed CTP state."""
from __future__ import annotations

import argparse
import hashlib
import json
from math import pi
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src"))

from recursive_horizons.nsc_background_projection import ChildFrameTensor
from recursive_horizons.nsc_child_metric_backreaction import (
    locked_neck_constraint_gate,
)


OUTPUT = ROOT/"results/development/nsc-child-metric-backreaction.json"
INPUTS = (
    "results/development/nsc-background-projection.json",
    "results/development/charged-compact-ctp-completion.json",
    "results/development/scale-binding.json",
    "results/development/adm-source-constraints.json",
    "results/development/spherical-action.json",
    "results/development/unruh-state.json",
)
SOURCES = (
    "src/recursive_horizons/nsc_child_metric_backreaction.py",
    "scripts/derive_nsc_child_metric_backreaction.py",
    "docs/nsc-child-metric-backreaction.md",
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
    projection, completion, scale, adm, spherical, parent = (
        json.loads((ROOT/path).read_text()) for path in INPUTS
    )
    locked = projection["locked_inputs"]
    branch = scale["development_branch"]
    if not (
        locked["q"] == branch["magnetic_flux"] == 4
        and locked["Omega"] == branch["omega"]
        and locked["zeta"] == branch["zeta"]
        and locked["V_full"] == branch["coefficients"]["V_full_relational"] == 0.0
    ):
        raise ValueError("locked scale branch changed")
    if adm["scope"]["joint_geometry_solution_found"]:
        raise ValueError("ADM input no longer describes the open evolution gate")
    if spherical["scope"]["self_sourced_geometry_or_physical_couplings_solved"]:
        raise ValueError("spherical input no longer describes the imported owner")
    required = parent["neck_source_budget"]["required_two_derivative_tensor"]
    if required["rho"] != "2*a_EH" or required["T_hatT_hatz"] != "0":
        raise ValueError("stored Bronnikov source convention changed")

    data = projection["completed_tensor_input"]
    tensor = ChildFrameTensor(
        rho=data["rho"], T01=data["T01"],
        p_parallel=data["p_parallel"], p_sphere=data["p_sphere"],
        parent_Killing_power=data["parent_Killing_power"],
    )
    geometry = projection["geometry_at_neck"]
    coefficient = branch["coefficients"]["A_Wilsonian_at_magnetic_scale"]
    newton = branch["mmp_map"]["G_N"]
    newton_residual = newton-1.0/(16.0*pi*coefficient)
    if abs(newton_residual) > 2e-14:
        raise ArithmeticError("scale-record Einstein normalization changed")
    result = locked_neck_constraint_gate(
        tensor,
        einstein_coefficient=coefficient,
        a_parallel=geometry["a_parallel"],
        sphere_radius=geometry["sphere_radius"],
        H_parallel=geometry["H_parallel"],
        H_sphere=geometry["H_sphere"],
    )
    constraints = result["normalized_constraints"]
    tolerance = 2e-10
    maximum = max(abs(
        constraints["hamiltonian_geometry_minus_source"]
    ), abs(constraints["momentum_geometry_minus_source"]))
    nulls = result["radial_nulls"]
    initial_constraints_pass = maximum < tolerance
    if initial_constraints_pass:
        raise AssertionError(
            "this record is the stored failing branch; a changed dependency "
            "requires a new evolution record"
        )
    if not (nulls["plus"] < 0 and nulls["minus"] < 0):
        raise AssertionError("completed source null signs changed")

    initial_row = {
        "T": 0.0,
        "a_parallel": geometry["a_parallel"],
        "r": geometry["sphere_radius"],
        "H_parallel": geometry["H_parallel"],
        "H_sphere": geometry["H_sphere"],
        "rho": tensor.rho,
        "T01": tensor.T01,
        "p_parallel": tensor.p_parallel,
        "p_sphere": tensor.p_sphere,
        "null_plus": nulls["plus"],
        "null_minus": nulls["minus"],
        "hamiltonian_constraint_residual": constraints[
            "hamiltonian_geometry_minus_source"
        ],
        "momentum_constraint_residual": constraints[
            "momentum_geometry_minus_source"
        ],
    }
    return {
        "schema": "NSC-CHILD-METRIC-BACKREACTION-v1",
        "status": (
            "FAIL AT T=0: the locked Bronnikov neck and completed charged CTP "
            "tensor do not lie on the same Hamiltonian/momentum constraint surface"
        ),
        "source_hashes": hashes(SOURCES),
        "input_hashes": hashes(INPUTS),
        "locked_inputs": {
            "q": locked["q"],
            "Omega": locked["Omega"],
            "zeta": locked["zeta"],
            "V_full": locked["V_full"],
            "A_at_magnetic_scale": coefficient,
            "metric_Cauchy_jet_changed": False,
            "state_or_scale_refitted": False,
        },
        "imported_owners": {
            "gravity": "existing ADM/spherical two-derivative owner",
            "geometry": "stored Bronnikov child neck and child clock",
            "state": "completed charged angular+compact CTP tensor",
            "standard_GR_rederived": False,
        },
        "initial_constraint_gate": {
            **result,
            "normalization_crosscheck": {
                "G_N": newton,
                "G_N_minus_one_over_16piA": newton_residual,
                "eight_pi_G_times_CTP_rho": 8.0*pi*newton*tensor.rho,
            },
            "declared_tolerance": tolerance,
            "maximum_absolute_normalized_constraint_residual": maximum,
            "passes": initial_constraints_pass,
        },
        "trajectory": {
            "rows": [initial_row],
            "requested_time_steps_started": 0,
            "stopped_at_T": 0.0,
            "stop_reason": (
                "ADM constraints precede evolution; stepping this slice would "
                "propagate invalid initial data rather than produce backreaction"
            ),
            "same_CTP_covariance_replayed": False,
            "reason_covariance_not_replayed": (
                "the algebraic initial constraint fails before a pressure or "
                "covariance update can affect the decision"
            ),
        },
        "decision": {
            "negative_null_signs_present_at_T0": True,
            "null_sign_break_time": None,
            "short_background_trajectory_defined": False,
            "failure_owner": (
                "constraint-complete same-action neck assembly: the source or "
                "source-determined initial geometry must supply the recorded "
                "density and momentum residuals before evolution"
            ),
            "required_constraint_completion": {
                "delta_rho": result["component_residual_required_minus_CTP"]["rho"],
                "delta_T01": result["component_residual_required_minus_CTP"]["T01"],
            },
            "longer_background_or_CTP_Hessian_authorized": False,
        },
        "scope": {
            "new_metric_ansatz": False,
            "Einstein_or_ADM_equations_rederived": False,
            "old_CTP_source_generator_rerun": False,
            "new_vacuum_subtraction_or_compensator": False,
            "q_Omega_zeta_or_V_full_changed": False,
            "dark_fraction_or_fourteen_over_nineteen_fit": False,
            "H_of_z_or_likelihood_computed": False,
            "CTP_Hessian_or_perturbations_computed": False,
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
        print("child backreaction T=0 obstruction reproduced; no CTP source was rerun")
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
