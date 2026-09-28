#!/usr/bin/env python3
"""Select constraint-complete neck data from the completed NSC CTP tensor."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src"))

from recursive_horizons.nsc_background_projection import ChildFrameTensor
from recursive_horizons.nsc_constraint_complete_neck import (
    source_selected_neck,
)


OUTPUT = ROOT/"results/development/nsc-constraint-complete-neck.json"
INPUTS = (
    "results/development/nsc-child-metric-backreaction.json",
    "results/development/nsc-background-projection.json",
    "results/development/scale-binding.json",
    "results/development/adm-source-constraints.json",
    "results/development/spherical-action.json",
    "results/development/compact-interaction.json",
    "results/development/gauge-source.json",
    "results/development/relational-vacuum-normalization.json",
)
SOURCES = (
    "src/recursive_horizons/nsc_constraint_complete_neck.py",
    "scripts/derive_nsc_constraint_complete_neck.py",
    "docs/nsc-constraint-complete-neck.md",
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
    failed, projection, scale, adm, spherical, interaction, gauge, vacuum = (
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
    if failed["initial_constraint_gate"]["passes"]:
        raise ValueError("stored Bronnikov branch no longer needs completion")
    if interaction["normalization"]["s_T"] is not None:
        raise ValueError("interacting compact channel is no longer undetermined")
    if vacuum["gate"]["V_full"] != 0.0:
        raise ValueError("relational vacuum law changed")
    if adm["scope"]["joint_geometry_solution_found"]:
        raise ValueError("ADM owner changed")
    if spherical["scope"]["self_sourced_geometry_or_physical_couplings_solved"]:
        raise ValueError("spherical owner changed")

    data = projection["completed_tensor_input"]
    tensor = ChildFrameTensor(
        rho=data["rho"], T01=data["T01"],
        p_parallel=data["p_parallel"], p_sphere=data["p_sphere"],
        parent_Killing_power=data["parent_Killing_power"],
    )
    seed = projection["geometry_at_neck"]
    coefficient = branch["coefficients"]["A_Wilsonian_at_magnetic_scale"]
    selected = source_selected_neck(
        tensor,
        einstein_coefficient=coefficient,
        retained_a_parallel=seed["a_parallel"],
        retained_H_parallel=seed["H_parallel"],
    )
    frame = selected["frame"]
    geometry = selected["geometry"]
    constraints = selected["constraints"]
    residuals = {
        **constraints,
        **frame["residuals"],
    }
    maximum = max(abs(value) for value in residuals.values())
    tolerance = 2e-12
    nulls = frame["radial_nulls"]
    gate = (
        maximum < tolerance
        and nulls["Landau_plus"] < 0
        and nulls["Landau_minus"] < 0
    )
    if not gate:
        raise ArithmeticError("source-selected initial geometry failed its gate")

    return {
        "schema": "NSC-CONSTRAINT-COMPLETE-NECK-v1",
        "status": (
            "PASS: geometry-from-source selects a Landau-frame Kantowski-Sachs "
            "neck whose Hamiltonian and momentum constraints close"
        ),
        "source_hashes": hashes(SOURCES),
        "input_hashes": hashes(INPUTS),
        "selected_route": {
            "choice": "B_geometry_from_source",
            "reason": (
                "no already evaluated same-action source channel has a fixed "
                "counterflow; the completed type-I tensor instead determines "
                "its constraint-compatible normal and radius"
            ),
        },
        "source_completion_A_audit": {
            "magnetic_gauge_Maxwell": (
                "already owned but diagonal on the homogeneous magnetic "
                "background; it does not supply the required T01 counterflow"
            ),
            "relational_vacuum": "V_full=0 is locked and supplies no flux",
            "interacting_compact_vertex": (
                "not selected because s_T and its interacting state remain undetermined"
            ),
            "unevaluated_nonlocal_remainder": (
                "no recorded finite tensor may be assigned to the residual"
            ),
            "new_term_added": False,
            "gauge_record_absolute_stress_complete": gauge["source_accounting"][
                "full_absolute_renormalized_stress_computed"
            ],
        },
        "locked_inputs": {
            "q": locked["q"],
            "Omega": locked["Omega"],
            "zeta": locked["zeta"],
            "A": coefficient,
            "V_full": locked["V_full"],
            "state_occupations_changed": False,
            "continuous_parameter_fitted": False,
        },
        "completed_tensor_old_child_frame": {
            "rho": tensor.rho,
            "T01": tensor.T01,
            "p_parallel": tensor.p_parallel,
            "p_sphere": tensor.p_sphere,
        },
        "source_selected_frame": frame,
        "source_selected_initial_geometry": {
            **geometry,
            "seed_r": seed["sphere_radius"],
            "r_minus_seed_r": geometry["r"]-seed["sphere_radius"],
            "a_parallel_and_H_parallel_owner": (
                "retained seed data; the two constraints do not fix the axial "
                "coordinate normalization or H_parallel at H_sphere=0"
            ),
            "Bronnikov_profile_retained": False,
            "MMP_unit_radius_imposed_on_new_geometry": False,
        },
        "constraints": {
            "equations_imported": True,
            "Hamiltonian": (
                "H_sphere^2+2*H_parallel*H_sphere+1/r^2-rho_L/(2*A)=0"
            ),
            "momentum": "T01_L=0 for homogeneous Kantowski-Sachs data",
            "residuals": residuals,
            "maximum_absolute_residual": maximum,
            "declared_tolerance": tolerance,
        },
        "gate": {
            "constraints_below_tolerance": maximum < tolerance,
            "both_radial_null_signs_negative": (
                nulls["Landau_plus"] < 0 and nulls["Landau_minus"] < 0
            ),
            "passes": gate,
            "covariance_metric_evolution_authorized_next": True,
            "next_owner": (
                "evolve the same CTP covariance and the source-selected metric "
                "together, re-evaluating stress at every step"
            ),
        },
        "scope": {
            "new_source_channel_or_metric_ansatz": False,
            "Einstein_or_ADM_equations_rederived": False,
            "A_q_Omega_zeta_or_V_full_changed": False,
            "fake_H_sphere_adjustment": False,
            "old_source_generator_rerun": False,
            "time_evolution_started": False,
            "H_of_z_Hessian_or_dark_fraction_fit": False,
            "source_selected_geometry_is_the_stored_Bronnikov_profile": False,
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
        print("constraint-complete source-selected neck reproduced; no source was rerun")
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
