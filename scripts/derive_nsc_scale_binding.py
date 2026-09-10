#!/usr/bin/env python3
"""Bind the common spectral coefficients to one recursive throat scale."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from recursive_horizons.nsc_scale_binding import (
    ScaleBindingConventions,
    ScaleBindingSolver,
)


OUTPUT = ROOT / "results/development/scale-binding.json"
INPUTS = (
    "results/development/relational-vacuum-normalization.json",
    "results/development/charged-sector.json",
    "results/development/compact-matching.json",
    "results/nsc-5-tail-limit.json",
    "results/nsc-5-clock-horizon.json",
)
SOURCES = (
    "src/recursive_horizons/nsc_scale_binding.py",
    "scripts/derive_nsc_scale_binding.py",
    "docs/nsc-scale-binding.md",
)


def hashes(paths):
    return {path: hashlib.sha256((ROOT / path).read_bytes()).hexdigest() for path in paths}


def compare(expected, actual, path="$", atol=4e-9, rtol=4e-8):
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
        if abs(expected - actual) > atol + rtol * abs(expected):
            raise AssertionError(f"numeric value differs at {path}")
    elif type(expected) is not type(actual) or expected != actual:
        raise AssertionError(f"exact value differs at {path}")


def branch(flux, bracket, modes, points, tolerance, transfer):
    solver = ScaleBindingSolver(ScaleBindingConventions(
        magnetic_flux=flux,
        compact_modes=modes,
        quadrature_points=points,
        quadrature_tolerance=tolerance,
    ))
    return solver.solve(bracket, recursive_transfer=transfer).to_dict()


def calculate():
    vacuum, charged, matching, tail, clock = (
        json.loads((ROOT / path).read_text()) for path in INPUTS
    )
    if not vacuum["gate"]["V_full_determined"] or vacuum["gate"]["V_full"] != 0:
        raise ValueError("the relational vacuum gate must be applied first")
    if not any(row["anomaly_free_zero_sector"] for row in charged["parity_and_anomaly_rows"]):
        raise ValueError("an anomaly-free charged zero sector is required")
    transfer = tail["families"][0]["q"]
    coarse = branch(4, (3.0, 5.0), 26, 96, 3e-8, transfer)
    refined = branch(4, (3.0, 5.0), 34, 112, 1e-9, transfer)
    lower = [
        branch(2, (1.5, 2.2), 34, 112, 1e-9, transfer),
        branch(3, (2.5, 3.5), 34, 112, 1e-9, transfer),
    ]
    if not all(abs(row["radius_residual"]) < 1e-9 for row in [coarse, refined, *lower]):
        raise AssertionError("radius solve did not close")
    if lower[0]["child_H_over_cutoff"] <= 1 or lower[1]["child_H_over_cutoff"] <= 1:
        raise AssertionError("lower-flux cutoff controls changed")
    if refined["child_H_over_cutoff"] >= 1:
        raise AssertionError("selected development branch does not resolve the child curvature")
    if refined["recursive_endpoint_product"] >= 1:
        raise AssertionError("selected branch left the stored limit-point regime")
    return {
        "schema": "NSC-SCALE-BINDING-v1",
        "status": (
            "first checked fixed-flux spectral/MMP radius branch that also places the "
            "stored child curvature below the cutoff"
        ),
        "source_hashes": hashes(SOURCES),
        "input_hashes": hashes(INPUTS),
        "equations": {
            "magnetic_matching_scale": "mu_B=sqrt(abs(q))/r_e",
            "Wilsonian_coefficients": (
                "A=(Q1-mu_B^2)/(6*(4*pi)^2), "
                "C=H_(Lambda,mu_B)(0)/(3*(4*pi)^2)"
            ),
            "relational_vacuum": "V_full=0",
            "MMP_radius": "r_e^2=q^2*C/(4*A)",
            "room_scale": "Omega=Lambda*r_e; zeta=Omega^2 in r_e=L_star units",
            "root": "q^2*C(Omega)/(4*A(Omega))-1=0",
        },
        "field_allocation": {
            "explicit_light_sector": "one charged four-dimensional Dirac zero field",
            "integrated_complement": "compact spectral modes above mu_B",
            "double_counted_light_loop": False,
            "matching_cutoff_fitted_to_root": False,
            "magnetic_scale_source": "MMP equation (2.3) running-coupling scale",
        },
        "development_branch": refined,
        "resolution_control": {
            "coarse": coarse,
            "absolute_Omega_change": abs(refined["omega"] - coarse["omega"]),
            "absolute_zeta_change": abs(refined["zeta"] - coarse["zeta"]),
        },
        "adjacent_lower_flux_branches": lower,
        "selection": {
            "q_checked": [2, 3, 4],
            "first_checked_integer_with_child_H_below_cutoff": 4,
            "q_selected_by_complete_action": False,
            "low_q_failure": "the radius roots for q=2 and q=3 lie below the stored child curvature scale",
            "MMP_large_q_control_satisfied_by_q4": False,
            "physical_branch_claimed": False,
        },
        "recursive_compatibility": {
            "stored_carrier": tail["conventions"]["geometry"],
            "endpoint_condition": tail["endpoint_theorem"]["limit_point_condition"],
            "Omega_times_transfer": refined["recursive_endpoint_product"],
            "no_remote_endpoint_datum_in_stored_carrier": refined["recursive_endpoint_product"] <= 1,
            "actual_charged_throat_recursive_map_evaluated": False,
        },
        "child_control": {
            "asymptotic_H_times_throat": clock["child_geometry"]["asymptotic_directional_H"],
            "curvature_resolved_by_cutoff": refined["child_H_over_cutoff"] < 1,
            "conditional_scalar_gap_real": refined["conditional_gap_squared_over_cutoff"] > 0,
            "scalar_gap_applied_to_massless_MMP_LLL": False,
        },
        "remaining_gate": (
            "evaluate the actual charged recursive boundary/state map and the full "
            "four-component CTP metric source before promoting this fixed-q branch"
        ),
        "old_scientific_generators_rerun": False,
        "comparison": {
            "fields": "all",
            "exact": "keys, types, strings and source/input hashes",
            "float_atol": 4e-9,
            "float_rtol": 4e-8,
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
        compare(json.loads(OUTPUT.read_text()), result)
        print("scale binding reproduced from authenticated inputs")
    elif args.output:
        if args.output.exists():
            raise FileExistsError("refusing to overwrite recorded evidence")
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
        print(args.output)
    else:
        print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
