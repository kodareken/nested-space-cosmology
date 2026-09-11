#!/usr/bin/env python3
"""Deposit the completed charged CTP neck tensor into the child ledger."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from recursive_horizons.nsc_background_projection import (
    ChildFrameTensor,
    bronnikov_child_slice,
    conserved_local_jet,
)


OUTPUT = ROOT / "results/development/nsc-background-projection.json"
INPUTS = (
    "results/development/charged-compact-ctp-completion.json",
    "results/development/charged-ctp-neck-source.json",
    "results/development/scale-binding.json",
    "results/development/adm-neck-source-map.json",
    "results/development/child-state.json",
    "results/nsc-5-clock-horizon.json",
)
SOURCES = (
    "src/recursive_horizons/nsc_background_projection.py",
    "scripts/derive_nsc_background_projection.py",
    "docs/nsc-background-projection.md",
)


def hashes(paths):
    return {
        path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
        for path in paths
    }


def compare(expected, actual, path="$", atol=2e-12, rtol=2e-12):
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


def tensor_dict(tensor):
    return {
        "rho": tensor.rho,
        "T01": tensor.T01,
        "p_parallel": tensor.p_parallel,
        "p_sphere": tensor.p_sphere,
        "parent_Killing_power": tensor.parent_Killing_power,
    }


def calculate():
    completion, angular, scale, adm, child_state, clock = (
        json.loads((ROOT/path).read_text()) for path in INPUTS
    )
    locked = completion["locked_inputs"]
    branch = scale["development_branch"]
    if not (
        locked["q"] == branch["magnetic_flux"] == 4
        and locked["Omega"] == branch["omega"]
        and locked["zeta"] == branch["zeta"]
        and locked["V_full"] == branch["coefficients"]["V_full_relational"] == 0.0
    ):
        raise ValueError("the locked q=4 scale branch changed")
    if adm["coframe_map"]["child_time"] != "dT=-d rho/sqrt(f)":
        raise ValueError("unexpected child clock convention")
    if "proper_volume_per_dz" not in json.dumps(child_state):
        raise ValueError("authenticated child proper-volume owner is missing")

    base = angular["runs"]["base"]
    massive = base["massive_angular_source"]
    lll = base["lll_source"]
    zero_level = ChildFrameTensor(
        rho=massive["source"]["rho"]+lll["child_frame"]["rho"],
        T01=(
            base["diagnostics"]["massive_child_T01"]
            + lll["child_frame"]["T01"]
        ),
        p_parallel=(
            massive["source"]["p_parallel"]
            + lll["child_frame"]["p_parallel"]
        ),
        p_sphere=(
            massive["source"]["p_sphere"]
            + lll["child_frame"]["p_sphere_state"]
        ),
        parent_Killing_power=(
            massive["parent_Killing_power"]
            + lll["PG_coordinate"]["parent_Killing_power"]
        ),
    )

    compact_nonlocal = completion["runs"]["base"]["totals"]
    compact_local = base["compact_local_source"]
    positive_compact = ChildFrameTensor(
        rho=compact_nonlocal["rho"]+compact_local["child_frame"]["rho"],
        T01=compact_nonlocal["T01"]+compact_local["child_frame"]["T01"],
        p_parallel=(
            compact_nonlocal["p_parallel"]
            + compact_local["child_frame"]["p_parallel"]
        ),
        p_sphere=(
            compact_nonlocal["p_sphere"]
            + compact_local["child_frame"]["p_sphere"]
        ),
        parent_Killing_power=(
            compact_nonlocal["parent_Killing_power"]
            + compact_local["Killing_power"]
        ),
    )
    reconstructed = zero_level.plus(positive_compact)
    completed = completion["completed_tensor"]
    total = ChildFrameTensor(
        rho=completed["rho"],
        T01=completed["T01"],
        p_parallel=completed["p_parallel"],
        p_sphere=completed["p_sphere"],
        parent_Killing_power=completed["parent_Killing_power"],
    )
    reconstruction = {
        name: getattr(reconstructed, name)-getattr(total, name)
        for name in tensor_dict(total)
    }
    if max(abs(value) for value in reconstruction.values()) > 2e-14:
        raise ArithmeticError("spectral allocation does not reconstruct tensor")

    geometry = bronnikov_child_slice(0.0)
    if abs(geometry["a_parallel"]**2-adm["geometry"]["f_minus_A"]) > 2e-14:
        raise ArithmeticError("stored child geometry conventions differ")
    if abs(clock["observer_map"]["geometry_probes"][2]["rho"]) > 0:
        raise ArithmeticError("stored neck probe moved")

    # The accepted realization is free and block diagonal in compact level.
    # There is no interaction vertex from which a regular inter-block Q could
    # arise.  Apply Q=0 to the total; component rows expose the allocation but
    # are not relabelled as visible or dark fluids.
    component_rows = []
    for name, tensor, allocation in (
        (
            "compact_zero_charged_angular_sector",
            zero_level,
            "LLL plus nonzero magnetic angular modes at compact level zero",
        ),
        (
            "positive_compact_and_Wilsonian_complement",
            positive_compact,
            "j=1,2 nonlocal CTP covariance plus the locked local compact response",
        ),
    ):
        component_rows.append({
            "name": name,
            "allocation": allocation,
            "tensor_at_neck": tensor_dict(tensor),
            "no_exchange_local_jet": conserved_local_jet(tensor, geometry),
        })
    jet = conserved_local_jet(total, geometry)
    tolerance = 2e-12
    residuals = {
        "tensor_reconstruction": reconstruction,
        "energy_conservation": jet["energy_conservation_residual"],
        "integrated_energy_ledger": jet["integrated_energy_ledger_residual"],
        "momentum_conservation": jet["momentum_conservation_residual"],
        "Killing_power_to_child_spatial_charge": jet["stored_Killing_power_map_residual"],
    }
    maximum = max(
        abs(value)
        for group in residuals.values()
        for value in (group.values() if isinstance(group, dict) else (group,))
    )
    if maximum >= tolerance:
        raise ArithmeticError("local child ledger does not close")

    rho_rate = jet["density_rate_child_time"]
    scale_rate = geometry["volume_scale_H"]
    return {
        "schema": "NSC-BACKGROUND-PROJECTION-v1",
        "status": (
            "PASS: the completed neck tensor has a covariantly conserved "
            "proper-volume/child-clock first jet; global background evolution "
            "requires metric backreaction with the same CTP source"
        ),
        "source_hashes": hashes(SOURCES),
        "input_hashes": hashes(INPUTS),
        "locked_inputs": {
            "q": locked["q"],
            "Omega": locked["Omega"],
            "zeta": locked["zeta"],
            "V_full": locked["V_full"],
            "dark_fraction_fit_used": False,
            "fourteen_over_nineteen_fit_used": False,
        },
        "imported_ledger": {
            "geometry": (
                "ds^2=dT^2-a_parallel(T)^2 dz^2-r(T)^2 dOmega_2^2"
            ),
            "child_clock": "dT=-d rho/sqrt(-A)",
            "proper_volume_per_coordinate_dz": "V=4*pi*a_parallel*r^2",
            "energy_conservation": (
                "rho_dot+H_parallel*(rho+p_parallel)"
                "+2*H_sphere*(rho+p_sphere)=Q"
            ),
            "momentum_conservation": (
                "j_dot+2*(H_parallel+H_sphere)*j=S_j"
            ),
            "derivation_repeated": False,
        },
        "geometry_at_neck": geometry,
        "completed_tensor_input": tensor_dict(total),
        "proper_volume_clock_projection": jet,
        "spectrum_component_ledger": {
            "rows": component_rows,
            "explicit_free_compact_block_transfer_Q": 0.0,
            "reason": (
                "the accepted free Gaussian Hamiltonian is block diagonal in "
                "compact level and the interacting compact vertex is disabled"
            ),
            "visible_dark_component_Q": None,
            "visible_dark_assignment_made": False,
            "allocation_is_a_renormalized_functional_split_not_a_dark_fluid_split": True,
            "reconstruction_residuals": reconstruction,
        },
        "homogeneous_output": {
            "type": "local first-order Cauchy jet at child time T=0",
            "rho_T": f"{total.rho:.17g}+{rho_rate:.17g}*T+O(T^2)",
            "normalized_volume_scale_aV_T": (
                f"1+{scale_rate:.17g}*T+O(T^2), "
                "aV=(V/V_neck)^(1/3)"
            ),
            "rho_aV": (
                f"{total.rho:.17g}+"
                f"{jet['density_derivative_per_log_volume_scale']:.17g}*log(aV)"
                "+O(log(aV)^2)"
            ),
            "p_parallel_at_neck": total.p_parallel,
            "p_sphere_at_neck": total.p_sphere,
            "pressure_time_derivatives": None,
            "global_rho_or_pressure_history_determined": False,
            "scope_boundary": (
                "one neck tensor fixes the conserved first density derivative "
                "but not the later pressure functions; evolve the same CTP "
                "state and metric together before constructing a background history"
            ),
        },
        "deposition_decision": {
            "regular_child_bulk_Q": 0.0,
            "proper_cell_energy_change_is_pressure_work": True,
            "parent_Killing_power_is_child_energy_injection": False,
            "parent_Killing_power_maps_to": (
                "minus the conserved child spatial-Killing momentum charge: "
                "P_parent=-P_z=-a_parallel*V*T01, with the sign fixed by "
                "dT=-d rho/a_parallel"
            ),
            "formation_surface_impulse_determined_by_neck_snapshot": False,
        },
        "residuals": {**residuals, "maximum_absolute": maximum},
        "gate": {
            "declared_tolerance": tolerance,
            "local_conservation_residual_below_tolerance": maximum < tolerance,
            "proper_volume_and_clock_projection_complete_at_neck": True,
            "global_H_of_z_authorized": False,
            "next_owner": (
                "metric backreaction and same-state CTP tensor evolution on "
                "the locked Bronnikov child branch"
            ),
        },
        "scope": {
            "new_source_computed": False,
            "new_metric_ansatz_or_equation": False,
            "old_scientific_generator_rerun": False,
            "vacuum_subtraction_changed": False,
            "q_Omega_or_zeta_changed": False,
            "dark_fraction_or_H_of_z_fitted": False,
            "CTP_Hessian_or_perturbations_computed": False,
            "global_background_claimed": False,
        },
        "comparison": {
            "fields": "all",
            "exact": "keys, types, strings and source/input hashes",
            "float_atol": 2e-12,
            "float_rtol": 2e-12,
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
        print("background projection reproduced from the completed tensor; no source generator was run")
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
